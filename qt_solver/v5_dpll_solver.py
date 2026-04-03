"""
V5 DPLL Solver: branch on products, propagate via Gaussian elimination.

Works in ORIGINAL variable space (not projected) where each quad eq
has exactly 1 product. This gives maximum sparsity for propagation.

Algorithm:
  1. Start with linear system + list of quad equations
  2. Pick unresolved quad equation with smallest variable count
  3. Try product = 0: add linear eq L_k(x) = 0, Gaussian eliminate
     - Check: does this determine any variables?
     - If yes: check if determined vars resolve other products → CASCADE
  4. If contradiction: try product = 1
  5. If both contradict: backtrack
  6. Repeat until all products resolved or contradiction

Key insight: in the original space, each quad eq has 1 product involving
variables from specific rounds (T_DEP). The triangular structure means
early-round products can be resolved first without affecting later rounds.
"""

import random
import time

from qt_solver.sha256_traced import (
    MASK32, IV256, sha256_compress, sha256_compress_traced,
    get_all_carry_chains, get_bit,
)
from qt_solver.gf2 import gf2_gaussian_eliminate, gf2_solve, gf2_rank
from qt_solver.qt_system import build_qt_system
from qt_solver.v3_solver import _iv_constant_value


def build_sparse_system(msg_words, num_rounds):
    """Build system keeping quad equations in original sparse form."""
    system, trace = build_qt_system(msg_words, num_rounds)
    vm = system.vm
    n = vm.num_vars

    linear_rows = list(system.to_linear_rows())
    quad_eqs = []  # (lin_row_int, (v1,v2), const)

    for eq in system.q_quadratic:
        new_lin = set(eq.linear)
        new_quad = set()
        new_const = eq.constant

        for v1, v2 in eq.quadratic:
            c1 = _iv_constant_value(vm, v1)
            c2 = _iv_constant_value(vm, v2)
            if c1 is not None and c2 is not None:
                new_const ^= (c1 & c2)
            elif c1 is not None:
                if c1 == 1:
                    new_lin ^= {v2}
            elif c2 is not None:
                if c2 == 1:
                    new_lin ^= {v1}
            else:
                new_quad.add((min(v1, v2), max(v1, v2)))

        if not new_quad:
            row = 0
            for v in new_lin:
                row ^= (1 << v)
            if new_const:
                row |= (1 << n)
            linear_rows.append(row)
        else:
            row = 0
            for v in new_lin:
                row ^= (1 << v)
            if new_const:
                row |= (1 << n)
            # Should have exactly 1 product after IV substitution
            quad_eqs.append((row, list(new_quad), new_const))

    return linear_rows, quad_eqs, system, trace, n


def dpll_solve(num_rounds, seed=42, max_nodes=500000, verbose=True):
    """
    DPLL solver with Gaussian propagation.
    """
    rng = random.Random(seed)
    msg = [rng.randint(0, MASK32) for _ in range(16)]

    if verbose:
        print(f"\n{'='*60}")
        print(f"V5 DPLL Solver: R={num_rounds}")
        print(f"{'='*60}")

    linear_rows, quad_eqs, system, trace, n = build_sparse_system(msg, num_rounds)
    vm = system.vm
    target = trace.hash_words

    if verbose:
        print(f"Linear rows: {len(linear_rows)}, Quad eqs: {len(quad_eqs)}")
        products_per_eq = [len(qp) for _, qp, _ in quad_eqs]
        if products_per_eq:
            print(f"Products per eq: min={min(products_per_eq)}, max={max(products_per_eq)}")

    # Preprocess: identify which quad eqs have exactly 1 product (most of them)
    single_product = []
    multi_product = []
    for lr, qp, const in quad_eqs:
        if len(qp) == 1:
            single_product.append((lr, qp[0], const))
        else:
            multi_product.append((lr, qp, const))

    if verbose:
        print(f"Single-product eqs: {len(single_product)}")
        print(f"Multi-product eqs: {len(multi_product)}")

    t0 = time.time()
    nodes = 0
    cascades = []

    # State: current linear rows + resolved products
    def solve_with_products(base_rows, product_assignments, remaining_singles):
        """
        Given product assignments, add them as linear equations,
        solve, check consistency and product validity.
        """
        rows = list(base_rows)
        for lr, (v1, v2), const, pval in product_assignments:
            var_mask = (1 << n) - 1
            row = lr & var_mask
            new_const = const ^ pval
            if new_const:
                row |= (1 << n)
            rows.append(row)

        result = gf2_solve(rows, n)
        if result is None:
            return None, None

        particular, kernel = result

        # Check if any remaining product is determined
        newly_determined = []
        still_remaining = []

        for lr, (v1, v2), const in remaining_singles:
            # Check if v1 or v2 is determined (not in any kernel vector)
            v1_free = any((kvec >> v1) & 1 for kvec in kernel)
            v2_free = any((kvec >> v2) & 1 for kvec in kernel)

            if not v1_free and not v2_free:
                # Both determined → product determined
                pval = ((particular >> v1) & 1) & ((particular >> v2) & 1)
                newly_determined.append((lr, (v1, v2), const, pval))
            elif not v1_free:
                v1_val = (particular >> v1) & 1
                if v1_val == 0:
                    # Product = 0 regardless
                    newly_determined.append((lr, (v1, v2), const, 0))
                else:
                    still_remaining.append((lr, (v1, v2), const))
            elif not v2_free:
                v2_val = (particular >> v2) & 1
                if v2_val == 0:
                    newly_determined.append((lr, (v1, v2), const, 0))
                else:
                    still_remaining.append((lr, (v1, v2), const))
            else:
                still_remaining.append((lr, (v1, v2), const))

        return (particular, kernel, newly_determined, still_remaining), rows

    def dpll(rows, assigned, remaining, depth=0):
        nonlocal nodes
        nodes += 1

        if nodes > max_nodes:
            return None

        if not remaining:
            # All products assigned → check full solution
            result = gf2_solve(rows, n)
            if result is None:
                return None
            particular, kernel = result

            # Verify all product constraints
            for lr, (v1, v2), const, pval in assigned:
                actual = ((particular >> v1) & 1) & ((particular >> v2) & 1)
                if actual != pval:
                    # Try sampling from kernel
                    for _ in range(min(100, 2 ** min(len(kernel), 7))):
                        sample = particular
                        for kvec in kernel:
                            if random.getrandbits(1):
                                sample ^= kvec
                        ok = True
                        for lr2, (va, vb), c2, pv2 in assigned:
                            if ((sample >> va) & 1) & ((sample >> vb) & 1) != pv2:
                                ok = False
                                break
                        if ok:
                            return sample
                    return None
            return particular

        # Pick equation to branch on (smallest variable count heuristic)
        best_idx = 0
        # Just take first for speed
        lr, (v1, v2), const = remaining[0]
        rest = remaining[1:]

        for pval in [0, 1]:
            # Add this product as linear equation
            new_rows = list(rows)
            var_mask = (1 << n) - 1
            row = lr & var_mask
            new_const = const ^ pval
            if new_const:
                row |= (1 << n)
            new_rows.append(row)

            # Quick consistency check
            result = gf2_solve(new_rows, n)
            if result is None:
                continue  # Inconsistent

            particular, kernel = result

            # Propagation: check if any remaining products are now determined
            new_assigned = list(assigned) + [(lr, (v1, v2), const, pval)]
            new_remaining = []
            cascade = 0

            for lr2, (va, vb), c2 in rest:
                va_free = any((kvec >> va) & 1 for kvec in kernel)
                vb_free = any((kvec >> vb) & 1 for kvec in kernel)

                if not va_free and not vb_free:
                    pv = ((particular >> va) & 1) & ((particular >> vb) & 1)
                    new_assigned.append((lr2, (va, vb), c2, pv))
                    # Add as linear constraint
                    row2 = lr2 & var_mask
                    nc2 = c2 ^ pv
                    if nc2:
                        row2 |= (1 << n)
                    new_rows.append(row2)
                    cascade += 1
                elif not va_free and ((particular >> va) & 1) == 0:
                    new_assigned.append((lr2, (va, vb), c2, 0))
                    row2 = lr2 & var_mask
                    if c2:
                        row2 |= (1 << n)
                    new_rows.append(row2)
                    cascade += 1
                elif not vb_free and ((particular >> vb) & 1) == 0:
                    new_assigned.append((lr2, (va, vb), c2, 0))
                    row2 = lr2 & var_mask
                    if c2:
                        row2 |= (1 << n)
                    new_rows.append(row2)
                    cascade += 1
                else:
                    new_remaining.append((lr2, (va, vb), c2))

            if cascade > 0:
                cascades.append(cascade)
                # Re-check consistency after cascade
                result2 = gf2_solve(new_rows, n)
                if result2 is None:
                    continue

            sol = dpll(new_rows, new_assigned, new_remaining, depth + 1)
            if sol is not None:
                return sol

        return None

    # DPLL on single-product eqs (sparse, cascade-friendly)
    # Then verify multi-product eqs on the solution
    solution = dpll(linear_rows, [], single_product)
    elapsed = time.time() - t0

    if verbose:
        print(f"\nDPLL: {nodes} nodes, {elapsed:.2f}s")
        if cascades:
            print(f"Cascades: {len(cascades)}, mean={sum(cascades)/len(cascades):.1f}, "
                  f"max={max(cascades)}")

    if solution is not None:
        found_msg = []
        for w in range(16):
            word = 0
            for b in range(32):
                if (solution >> vm.w(w, b)) & 1:
                    word |= (1 << b)
            found_msg.append(word)

        found_hash = sha256_compress(found_msg, num_rounds)
        success = found_hash == target and found_msg != msg

        if verbose:
            print(f"Solution found! Hash match: {found_hash == target}")
            if success:
                print(f"*** SECOND PREIMAGE CONFIRMED ***")

        return {
            'success': success,
            'nodes': nodes,
            'cascades': cascades,
            'elapsed': elapsed,
            'msg': found_msg,
        }
    else:
        if verbose:
            print(f"No solution (nodes={nodes}, timeout={nodes >= max_nodes})")
        return {
            'success': False,
            'nodes': nodes,
            'elapsed': elapsed,
        }


def benchmark_dpll(max_R=8, seeds=3):
    """Benchmark DPLL solver."""
    print(f"\n{'═'*60}")
    print(f"V5 DPLL BENCHMARK")
    print(f"{'═'*60}")

    for R in range(2, max_R + 1):
        successes = 0
        total_nodes = 0
        for s in range(seeds):
            r = dpll_solve(R, seed=s * 100 + 42, max_nodes=200000, verbose=False)
            if r.get('success'):
                successes += 1
            total_nodes += r.get('nodes', 0)
        avg_nodes = total_nodes // seeds
        print(f"R={R}: {successes}/{seeds} success, avg_nodes={avg_nodes}")
