"""
Six parallel mathematical directions for Q∩T solving.

Dir 1: Product Enumeration — enumerate p_k = x_i·x_j, solve linear
Dir 2: Decomposition — solve 29 components independently
Dir 3: Kernel-Product Interaction — measure rank change per product
Dir 4: Carry Window Partitioning — group by carry locality
Dir 5: Incremental Telescope — extend R=2 solution to R=3,4,...
Dir 6: Algebraic DPLL — fix product → propagate → cascade
"""

import random
import time
import math
from collections import defaultdict

from qt_solver.sha256_traced import (
    MASK32, sha256_compress, sha256_compress_traced, get_all_carry_chains,
)
from qt_solver.gf2 import gf2_gaussian_eliminate, gf2_solve, gf2_rank
from qt_solver.qt_system import build_qt_system, VarMap, QTSystem
from qt_solver.relinearize import iterative_relinearize, iter_bits


def _extract_msg(bv, vm):
    msg = []
    for w in range(16):
        word = 0
        for b in range(32):
            if (bv >> vm.w(w, b)) & 1:
                word |= (1 << b)
        msg.append(word)
    return msg


def _get_remaining_quad(msg, R):
    """Get the remaining quad system after relinearization."""
    system, trace = build_qt_system(msg, R)
    n = system.vm.num_vars
    linear_rows = system.to_linear_rows()
    # Convert GF2Equation objects to (lin_row, quad_pairs, const) tuples
    quad_tuples = []
    for eq in system.q_quadratic:
        lr = 0
        for v in eq.linear:
            lr ^= (1 << v)
        if eq.constant:
            lr |= (1 << n)
        quad_tuples.append((lr, list(eq.quadratic), eq.constant))
    relin = iterative_relinearize(
        linear_rows, quad_tuples, n, verbose=False
    )
    return system, trace, relin, n


# ════════════════════════════════════════════════════════════
# Dir 1: Product Enumeration
# ════════════════════════════════════════════════════════════

def product_enumeration_analysis(num_rounds=3, seed=42, verbose=True):
    """
    For each product assignment (p_1,...,p_K):
      - Add K linear equations to existing system
      - Check consistency
      - If consistent: check product constraints on solution

    For small K: exhaustive. For large K: sample.
    """
    rng = random.Random(seed)
    msg = [rng.randint(0, MASK32) for _ in range(16)]
    system, trace, relin, n = _get_remaining_quad(msg, num_rounds)
    remaining = relin['quad_equations']
    vm = system.vm

    if verbose:
        print(f"\n{'='*60}")
        print(f"Dir 1: Product Enumeration (R={num_rounds})")
        print(f"{'='*60}")
        print(f"Remaining quad: {len(remaining)}, Free vars: {relin['final_kernel_dim']}")

    if not remaining:
        print("No quadratic equations!")
        return {}

    # Extract products and their linear equations
    products = []  # (linear_row_in_free_space, quad_pair, const)
    for lr, qpairs, const in remaining:
        assert len(qpairs) == 1, f"Expected 1 product, got {len(qpairs)}"
        products.append((lr, qpairs[0], const))

    K = len(products)
    if verbose:
        print(f"Products: {K}")

    # Get base linear system
    base_rows = list(system.to_linear_rows())
    # Add linearized quad eqs from relinearization
    from qt_solver.relinearize import relinearize_once
    new_lin, rem_quad, _ = relinearize_once(
        system.linear_rows, system.quad_equations, n, verbose=False
    )
    base_rows = new_lin  # Already includes relinearized eqs

    base_rank = gf2_rank(base_rows, n)

    if verbose:
        print(f"Base rank: {base_rank}, Base kernel: {n - base_rank}")

    # Test: how many product assignments are consistent?
    max_enum = min(2**K, 2**18)  # Cap at 2^18 for speed
    sample_mode = (2**K > max_enum)

    consistent = 0
    product_valid = 0
    t0 = time.time()

    for trial in range(max_enum):
        if sample_mode:
            pvec = rng.getrandbits(K)
        else:
            pvec = trial

        # Add product-as-linear equations
        trial_rows = list(base_rows)
        for k in range(K):
            pk = (pvec >> k) & 1
            lr, (v1, v2), const = products[k]
            # Original eq: lin ⊕ v1*v2 = 0
            # With p_k = v1*v2: lin = p_k, i.e., lin ⊕ p_k = 0
            var_mask = (1 << n) - 1
            row = lr & var_mask
            new_const = const ^ pk
            if new_const:
                row |= (1 << n)
            trial_rows.append(row)

        # Check consistency
        result = gf2_solve(trial_rows, n)
        if result is not None:
            consistent += 1
            particular, kernel = result

            # Check product constraints
            assignment = {}
            for v in range(n):
                assignment[v] = (particular >> v) & 1

            all_valid = True
            for k in range(K):
                _, (v1, v2), _ = products[k]
                expected = (pvec >> k) & 1
                actual = assignment.get(v1, 0) & assignment.get(v2, 0)
                if actual != expected:
                    all_valid = False
                    break

            if all_valid:
                product_valid += 1
                # VERIFY: does this give correct hash?
                found_msg = _extract_msg(particular, vm)
                found_hash = sha256_compress(found_msg, num_rounds)
                if found_hash == trace.hash_words:
                    if verbose:
                        print(f"\n  *** PREIMAGE at trial {trial}! ***")
                        elapsed = time.time() - t0
                        print(f"  Time: {elapsed:.2f}s")
                    return {
                        'success': True,
                        'trial': trial,
                        'msg': found_msg,
                    }

    elapsed = time.time() - t0
    if verbose:
        mode = f"sampled {max_enum}" if sample_mode else f"exhaustive {max_enum}"
        print(f"\n  {mode} product vectors in {elapsed:.1f}s:")
        print(f"  Consistent: {consistent}/{max_enum} ({100*consistent/max_enum:.1f}%)")
        print(f"  Product-valid: {product_valid}/{max_enum}")
        if consistent > 0:
            print(f"  P(consistent) ≈ 2^{math.log2(consistent/max_enum):.1f}")

    return {
        'K': K,
        'tested': max_enum,
        'consistent': consistent,
        'product_valid': product_valid,
    }


# ════════════════════════════════════════════════════════════
# Dir 3: Kernel-Product Interaction
# ════════════════════════════════════════════════════════════

def kernel_product_interaction(num_rounds=3, seed=42, verbose=True):
    """
    Measure how adding each product equation changes the rank.
    Some products may be redundant (already implied by linear system).
    """
    rng = random.Random(seed)
    msg = [rng.randint(0, MASK32) for _ in range(16)]
    system, trace, relin, n = _get_remaining_quad(msg, num_rounds)
    remaining = relin['quad_equations']

    if verbose:
        print(f"\n{'='*60}")
        print(f"Dir 3: Kernel-Product Interaction (R={num_rounds})")
        print(f"{'='*60}")

    from qt_solver.relinearize import relinearize_once
    base_rows, _, _ = relinearize_once(
        system.linear_rows, system.quad_equations, n, verbose=False
    )
    base_rank = gf2_rank(base_rows, n)

    # Add products one by one (with p=0) and measure rank change
    cumulative_rows = list(base_rows)
    rank_progression = [base_rank]

    for k, (lr, qpairs, const) in enumerate(remaining):
        var_mask = (1 << n) - 1
        row = lr & var_mask  # linear part (with p=0: lin ⊕ 0 = 0 → lin = 0)
        if const:
            row |= (1 << n)
        cumulative_rows.append(row)
        new_rank = gf2_rank(cumulative_rows, n)
        rank_progression.append(new_rank)

    if verbose:
        total_increase = rank_progression[-1] - rank_progression[0]
        print(f"Base rank: {base_rank}")
        print(f"After all {len(remaining)} products (p=0): rank={rank_progression[-1]}")
        print(f"Total rank increase: {total_increase}")
        print(f"Kernel: {n - base_rank} → {n - rank_progression[-1]}")
        print(f"Compression: {total_increase} bits")

        # How many products are redundant?
        increments = [rank_progression[i+1] - rank_progression[i]
                      for i in range(len(remaining))]
        redundant = sum(1 for inc in increments if inc == 0)
        print(f"Independent products: {sum(1 for i in increments if i > 0)}/{len(remaining)}")
        print(f"Redundant products: {redundant}/{len(remaining)}")

    return {
        'base_rank': base_rank,
        'final_rank': rank_progression[-1],
        'progression': rank_progression,
    }


# ════════════════════════════════════════════════════════════
# Dir 5: Incremental Telescope
# ════════════════════════════════════════════════════════════

def incremental_telescope(max_rounds=8, seed=42, verbose=True):
    """
    Solve R=2 → extend to R=3 → extend to R=4 → ...

    At each step: use the previous solution as starting point
    for the greedy combined walk on the extended system.

    Key idea: the R-round solution is a GOOD starting point for
    R+1 rounds because most of the state is already correct.
    """
    rng = random.Random(seed)

    if verbose:
        print(f"\n{'='*60}")
        print(f"Dir 5: Incremental Telescope (R=2→{max_rounds})")
        print(f"{'='*60}")

    # Start: find R=2 preimage
    base_msg = [rng.randint(0, MASK32) for _ in range(16)]

    current_msg = list(base_msg)
    target_hashes = {}

    for R in range(2, max_rounds + 1):
        # Build system for R rounds with ORIGINAL message's hash as target
        trace_orig = sha256_compress_traced(base_msg, R)
        target = trace_orig.hash_words
        target_hashes[R] = target

        system, trace = build_qt_system(current_msg, R, target)
        vm = system.vm
        ref_carries = get_all_carry_chains(trace)

        rows = system.to_linear_rows()
        result = gf2_solve(rows, vm.num_vars)
        if result is None:
            if verbose:
                print(f"  R={R}: linear system inconsistent — rebuilding")
            # Rebuild from base_msg
            system, trace = build_qt_system(base_msg, R)
            ref_carries = get_all_carry_chains(trace)
            rows = system.to_linear_rows()
            result = gf2_solve(rows, vm.num_vars)
            if result is None:
                print(f"  R={R}: FATAL inconsistency")
                break

        particular, kernel = result
        quad_eqs = system.q_quadratic

        # Greedy combined walk
        current_bv = particular
        cost, cmm, qviol = _eval_cost(current_bv, vm, ref_carries, R, quad_eqs)
        best_cost = cost
        best_bv = current_bv

        for step in range(min(500, 100 + R * 50)):
            improved = False
            for ki, kvec in enumerate(kernel):
                if kvec & ((1 << 512) - 1) == 0:
                    continue
                cand = current_bv ^ kvec
                c, cm, qv = _eval_cost(cand, vm, ref_carries, R, quad_eqs)
                if c < cost:
                    cost, cmm, qviol = c, cm, qv
                    current_bv = cand
                    improved = True
                    if c < best_cost:
                        best_cost = c
                        best_bv = cand
                    break

            if not improved:
                for _ in range(rng.randint(1, 5)):
                    current_bv ^= kernel[rng.randint(0, len(kernel) - 1)]
                cost, cmm, qviol = _eval_cost(current_bv, vm, ref_carries, R, quad_eqs)

            if best_cost == 0:
                break

        found_msg = _extract_msg(best_bv, vm)
        found_hash = sha256_compress(found_msg, R)
        success = (found_hash == target)

        if verbose:
            print(f"  R={R}: cost={best_cost} (carry={cmm},quad={qviol}), "
                  f"hash_match={success}, steps={step+1}")

        if success:
            current_msg = found_msg  # Use as starting point for R+1
        else:
            if verbose:
                print(f"    → failed, keeping previous message")
            break

    return target_hashes


def _eval_cost(bv, vm, ref_carries, R, quad_eqs):
    msg = _extract_msg(bv, vm)
    trace = sha256_compress_traced(msg, R)
    carries = get_all_carry_chains(trace)
    cmm = sum(bin(ref_carries[i] ^ carries[i]).count('1')
              for i in range(len(ref_carries)))
    assignment = {v: (bv >> v) & 1 for v in range(vm.num_vars)}
    qv = sum(1 for eq in quad_eqs if eq.evaluate(assignment) != 0)
    return cmm + qv * 5, cmm, qv


# ════════════════════════════════════════════════════════════
# Dir 6: Algebraic DPLL with Cascade Propagation
# ════════════════════════════════════════════════════════════

def algebraic_dpll(num_rounds=3, seed=42, verbose=True):
    """
    Fix products one by one with propagation.

    1. Pick a product p_k
    2. Try p_k = 0: add linear eq, Gaussian eliminate
       - If new rank determines variables that appear in other products
         → those products become determined → cascade!
    3. If contradiction: try p_k = 1
    4. If both contradict: backtrack

    Measure cascade lengths and search tree depth.
    """
    rng = random.Random(seed)
    msg = [rng.randint(0, MASK32) for _ in range(16)]
    system, trace, relin, n = _get_remaining_quad(msg, num_rounds)
    remaining = relin['quad_equations']
    vm = system.vm
    target = trace.hash_words

    if verbose:
        print(f"\n{'='*60}")
        print(f"Dir 6: Algebraic DPLL (R={num_rounds})")
        print(f"{'='*60}")
        print(f"Products: {len(remaining)}")

    from qt_solver.relinearize import relinearize_once
    base_rows, _, _ = relinearize_once(
        system.linear_rows, system.quad_equations, n, verbose=False
    )

    # Extract product info
    products = []
    for lr, qpairs, const in remaining:
        var_mask = (1 << n) - 1
        lin_bits = lr & var_mask
        products.append({
            'lin_row': lr,
            'pair': qpairs[0],
            'const': const,
            'v1': qpairs[0][0],
            'v2': qpairs[0][1],
        })

    # DPLL search
    assigned_products = {}  # k -> p_value
    current_rows = list(base_rows)
    cascade_lengths = []
    nodes_explored = 0
    stack = []  # for backtracking: (k, old_rows_len, old_assigned)

    def propagate(rows, assigned):
        """After adding equations, check if any products become determined."""
        result = gf2_solve(rows, n)
        if result is None:
            return None, 0  # Inconsistent

        particular, kernel = result
        cascade = 0

        # Check unassigned products
        for k, prod in enumerate(products):
            if k in assigned:
                continue
            v1, v2 = prod['v1'], prod['v2']

            # Check if v1 or v2 is determined (appears in no kernel vector)
            v1_free = any((kvec >> v1) & 1 for kvec in kernel)
            v2_free = any((kvec >> v2) & 1 for kvec in kernel)

            if not v1_free:
                # v1 is determined
                v1_val = (particular >> v1) & 1
                if v1_val == 0:
                    # x_i = 0 → x_i·x_j = 0 regardless of x_j
                    assigned[k] = 0
                    var_mask = (1 << n) - 1
                    row = prod['lin_row'] & var_mask
                    if prod['const']:
                        row |= (1 << n)
                    rows.append(row)
                    cascade += 1
                elif not v2_free:
                    # Both determined
                    v2_val = (particular >> v2) & 1
                    pval = v1_val & v2_val
                    assigned[k] = pval
                    var_mask = (1 << n) - 1
                    row = prod['lin_row'] & var_mask
                    new_const = prod['const'] ^ pval
                    if new_const:
                        row |= (1 << n)
                    rows.append(row)
                    cascade += 1

            elif not v2_free:
                v2_val = (particular >> v2) & 1
                if v2_val == 0:
                    assigned[k] = 0
                    var_mask = (1 << n) - 1
                    row = prod['lin_row'] & var_mask
                    if prod['const']:
                        row |= (1 << n)
                    rows.append(row)
                    cascade += 1

        return result, cascade

    # Initial propagation
    result0, cascade0 = propagate(current_rows, assigned_products)
    if verbose:
        print(f"Initial cascade: {cascade0} products determined")
        print(f"Assigned: {len(assigned_products)}/{len(products)}")

    # DPLL loop
    unassigned = [k for k in range(len(products)) if k not in assigned_products]

    if verbose and unassigned:
        print(f"Remaining to assign: {len(unassigned)}")

    t0 = time.time()

    def dpll_solve(rows, assigned, depth=0):
        nonlocal nodes_explored
        nodes_explored += 1

        if nodes_explored > 100000:
            return None  # Timeout

        unassigned = [k for k in range(len(products)) if k not in assigned]
        if not unassigned:
            # All products assigned — check full solution
            result = gf2_solve(rows, n)
            if result is None:
                return None
            particular, kernel = result
            # Verify product constraints
            for k, prod in enumerate(products):
                v1_val = (particular >> prod['v1']) & 1
                v2_val = (particular >> prod['v2']) & 1
                expected = assigned[k]
                if (v1_val & v2_val) != expected:
                    return None
            return particular

        # Pick next product to assign (first unassigned)
        k = unassigned[0]
        prod = products[k]

        for pval in [0, 1]:
            new_assigned = dict(assigned)
            new_assigned[k] = pval
            new_rows = list(rows)

            # Add linear equation
            var_mask = (1 << n) - 1
            row = prod['lin_row'] & var_mask
            new_const = prod['const'] ^ pval
            if new_const:
                row |= (1 << n)
            new_rows.append(row)

            # Propagate
            result, cascade = propagate(new_rows, new_assigned)
            if result is None:
                continue  # Inconsistent, try other value

            if cascade > 0:
                cascade_lengths.append(cascade)

            # Recurse
            solution = dpll_solve(new_rows, new_assigned, depth + 1)
            if solution is not None:
                return solution

        return None  # Both values failed

    solution = dpll_solve(current_rows, assigned_products)
    elapsed = time.time() - t0

    if verbose:
        print(f"\nDPLL result ({elapsed:.2f}s):")
        print(f"  Nodes explored: {nodes_explored}")
        if cascade_lengths:
            print(f"  Cascades: {len(cascade_lengths)}, "
                  f"mean length={sum(cascade_lengths)/len(cascade_lengths):.1f}, "
                  f"max={max(cascade_lengths)}")

        if solution is not None:
            found_msg = _extract_msg(solution, vm)
            found_hash = sha256_compress(found_msg, num_rounds)
            match = found_hash == target
            print(f"  Solution found! Hash match: {match}")
            if match:
                print(f"  *** PREIMAGE CONFIRMED ***")
        else:
            print(f"  No solution found (timeout or inconsistent)")

    return {
        'nodes': nodes_explored,
        'cascades': cascade_lengths,
        'solution': solution is not None,
        'elapsed': elapsed,
    }


# ════════════════════════════════════════════════════════════
# Run all directions
# ════════════════════════════════════════════════════════════

def run_all(R=3, seed=42):
    """Run all 6 directions for round R."""
    print(f"\n{'═'*60}")
    print(f"ALL 6 DIRECTIONS: R={R}")
    print(f"{'═'*60}")

    product_enumeration_analysis(R, seed)
    kernel_product_interaction(R, seed)
    algebraic_dpll(R, seed)
    incremental_telescope(max_rounds=R+2, seed=seed)
