"""
V6 Hybrid Solver: DPLL cascade + greedy walk.

Phase 1 (DPLL): resolve single-product equations via cascade.
  - 32 single-product eqs from round 1 (IV-adjacent)
  - Cascade determines ~31 from just 1 branch
  - Adds ~32 linear equations to the system
  - Kernel shrinks by ~32, but 32 quad eqs eliminated

Phase 2 (Greedy): walk the REDUCED Q-variety.
  - Kernel is smaller but cleaner (round-1 products resolved)
  - Only multi-product eqs remain
  - Combined cost: carry_mismatch + quad_violations
"""

import random
import time
import math

from qt_solver.sha256_traced import (
    MASK32, IV256, sha256_compress, sha256_compress_traced,
    get_all_carry_chains, get_bit,
)
from qt_solver.gf2 import gf2_gaussian_eliminate, gf2_solve, gf2_rank
from qt_solver.qt_system import build_qt_system
from qt_solver.v3_solver import _iv_constant_value


def _extract_msg(bv, vm):
    msg = []
    for w in range(16):
        word = 0
        for b in range(32):
            if (bv >> vm.w(w, b)) & 1:
                word |= (1 << b)
        msg.append(word)
    return msg


def build_hybrid_system(msg_words, num_rounds):
    """
    Build system and separate into:
    - linear_rows: all linear equations + IV-linearized quads
    - single_quads: equations with exactly 1 product (DPLL targets)
    - multi_quads: equations with 2+ products (greedy targets)
    """
    system, trace = build_qt_system(msg_words, num_rounds)
    vm = system.vm
    n = vm.num_vars

    linear_rows = list(system.to_linear_rows())
    single_quads = []
    multi_quads = []

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
        elif len(new_quad) == 1:
            row = 0
            for v in new_lin:
                row ^= (1 << v)
            pair = list(new_quad)[0]
            single_quads.append((row, pair, new_const))
        else:
            multi_quads.append((new_lin, new_quad, new_const))

    return linear_rows, single_quads, multi_quads, system, trace, n


def dpll_cascade(linear_rows, single_quads, n, max_nodes=10000):
    """
    Phase 1: DPLL on single-product equations.
    Returns augmented linear_rows with products resolved.
    """
    result = gf2_solve(linear_rows, n)
    if result is None:
        return None, []

    particular, kernel = result

    # Try to resolve single_quads via cascade
    remaining = list(single_quads)
    resolved = []
    rows = list(linear_rows)

    # Try each single-product eq as seed
    for seed_idx in range(min(len(remaining), 4)):
        if not remaining:
            break

        lr, (v1, v2), const = remaining[0]

        for pval in [0, 1]:
            trial_rows = list(rows)
            trial_resolved = list(resolved)
            trial_remaining = list(remaining[1:])

            # Add seed product
            var_mask = (1 << n) - 1
            row = lr & var_mask
            nc = const ^ pval
            if nc:
                row |= (1 << n)
            trial_rows.append(row)
            trial_resolved.append((lr, (v1, v2), const, pval))

            # Cascade: check if other products are now determined
            changed = True
            while changed:
                changed = False
                res = gf2_solve(trial_rows, n)
                if res is None:
                    break  # Inconsistent
                part, kern = res

                new_remaining = []
                for lr2, (va, vb), c2 in trial_remaining:
                    va_free = any((k >> va) & 1 for k in kern)
                    vb_free = any((k >> vb) & 1 for k in kern)

                    determined = False
                    det_pval = 0

                    if not va_free and not vb_free:
                        det_pval = ((part >> va) & 1) & ((part >> vb) & 1)
                        determined = True
                    elif not va_free and ((part >> va) & 1) == 0:
                        det_pval = 0
                        determined = True
                    elif not vb_free and ((part >> vb) & 1) == 0:
                        det_pval = 0
                        determined = True

                    if determined:
                        row2 = lr2 & var_mask
                        nc2 = c2 ^ det_pval
                        if nc2:
                            row2 |= (1 << n)
                        trial_rows.append(row2)
                        trial_resolved.append((lr2, (va, vb), c2, det_pval))
                        changed = True
                    else:
                        new_remaining.append((lr2, (va, vb), c2))

                trial_remaining = new_remaining

            # Check consistency
            if res is not None:
                rows = trial_rows
                resolved = trial_resolved
                remaining = trial_remaining
                break
        else:
            # Both values failed for seed, skip it
            remaining = remaining[1:]

    return rows, resolved


def v6_hybrid(num_rounds, max_steps=1500, seed=42, verbose=True):
    """
    V6 Hybrid: DPLL cascade → greedy walk.
    """
    rng = random.Random(seed)
    msg = [rng.randint(0, MASK32) for _ in range(16)]

    if verbose:
        print(f"\n{'='*60}")
        print(f"V6 Hybrid Solver: R={num_rounds}")
        print(f"{'='*60}")

    linear_rows, single_quads, multi_quads, system, trace, n = \
        build_hybrid_system(msg, num_rounds)
    vm = system.vm
    target = trace.hash_words
    ref_carries = get_all_carry_chains(trace)

    if verbose:
        print(f"Phase 0: {len(linear_rows)} linear, "
              f"{len(single_quads)} single-prod, {len(multi_quads)} multi-prod")

    # Phase 1: DPLL cascade
    t0 = time.time()
    augmented_rows, resolved = dpll_cascade(linear_rows, single_quads, n)

    if augmented_rows is None:
        if verbose:
            print("Phase 1: DPLL inconsistent!")
        return {'success': False}

    if verbose:
        print(f"Phase 1 (DPLL): resolved {len(resolved)}/{len(single_quads)} "
              f"single-product eqs via cascade ({time.time()-t0:.2f}s)")

    # Phase 2: Greedy walk on reduced system
    result = gf2_solve(augmented_rows, n)
    if result is None:
        if verbose:
            print("Phase 2: augmented system inconsistent!")
        return {'success': False}

    particular, kernel = result

    if verbose:
        print(f"Phase 2: kernel={len(kernel)}, multi-quads={len(multi_quads)}")
        print(f"  Effective DOF: {len(kernel) - len(multi_quads)}")

    if not multi_quads:
        # All quads resolved by DPLL! Just verify hash.
        found_msg = _extract_msg(particular, vm)
        found_hash = sha256_compress(found_msg, num_rounds)
        success = found_hash == target and found_msg != msg
        if verbose:
            print(f"  No multi-quads remaining! Hash match: {found_hash == target}")
            if success:
                print(f"  *** SECOND PREIMAGE FOUND (DPLL only)! ***")
        return {'success': success, 'steps': 0, 'kernel_dim': len(kernel),
                'remaining_quad': 0, 'dpll_resolved': len(resolved)}

    # Greedy + SA walk on multi-quads
    def evaluate(bv):
        m = _extract_msg(bv, vm)
        t = sha256_compress_traced(m, num_rounds)
        c = get_all_carry_chains(t)
        cmm = sum(bin(ref_carries[i] ^ c[i]).count('1') for i in range(len(ref_carries)))
        assignment = {v: (bv >> v) & 1 for v in range(n)}
        qv = 0
        for lin_terms, quad_terms, const in multi_quads:
            val = const
            for v in lin_terms:
                val ^= assignment.get(v, 0)
            for v1, v2 in quad_terms:
                val ^= (assignment.get(v1, 0) & assignment.get(v2, 0))
            if val != 0:
                qv += 1
        return cmm + qv * 5, cmm, qv

    current = particular
    cost, cmm, qviol = evaluate(current)
    best_cost, best_bv = cost, current

    if verbose:
        print(f"  Start: cost={cost} (carry={cmm}, quad={qviol})")

    use_sa = num_rounds >= 5
    temperature = 80.0 if use_sa else 0
    cooling = 0.997

    for step in range(max_steps):
        improved = False

        if use_sa:
            trials = min(60, len(kernel))
            indices = rng.sample(range(len(kernel)), trials) if len(kernel) > trials else range(len(kernel))
        else:
            indices = range(len(kernel))

        for ki in indices:
            kvec = kernel[ki]
            if kvec & ((1 << 512) - 1) == 0:
                continue
            cand = current ^ kvec
            c, cm, qv = evaluate(cand)
            if c < cost:
                cost, cmm, qviol = c, cm, qv
                current = cand
                improved = True
                if c < best_cost:
                    best_cost = c
                    best_bv = cand
                if not use_sa:
                    break
            elif use_sa and temperature > 0.1:
                delta = c - cost
                if delta < temperature and rng.random() < math.exp(-delta / temperature):
                    cost, cmm, qviol = c, cm, qv
                    current = cand
                    improved = True

        if use_sa:
            temperature *= cooling

        if not improved:
            nf = rng.randint(2, min(10, len(kernel)))
            for _ in range(nf):
                current ^= kernel[rng.randint(0, len(kernel) - 1)]
            cost, cmm, qviol = evaluate(current)

        if best_cost == 0:
            found_msg = _extract_msg(best_bv, vm)
            found_hash = sha256_compress(found_msg, num_rounds)
            success = found_hash == target and found_msg != msg
            if verbose:
                print(f"\n  cost=0 at step {step}! Hash: {found_hash == target}")
                if success:
                    print(f"  *** SECOND PREIMAGE FOUND! ***")
            return {'success': success, 'steps': step,
                    'kernel_dim': len(kernel),
                    'remaining_quad': len(multi_quads),
                    'dpll_resolved': len(resolved)}

        if verbose and step % 200 == 0:
            sa_info = f", T={temperature:.1f}" if use_sa else ""
            print(f"  step {step}: cost={cost} (c={cmm},q={qviol}), "
                  f"best={best_cost}{sa_info}, {time.time()-t0:.1f}s")

    elapsed = time.time() - t0
    if verbose:
        print(f"\n  Best: cost={best_cost} ({elapsed:.1f}s)")
    return {'success': False, 'best_cost': best_cost,
            'kernel_dim': len(kernel),
            'remaining_quad': len(multi_quads),
            'dpll_resolved': len(resolved)}


def benchmark_v6(max_R=8, seeds=5):
    """Benchmark V6 hybrid solver."""
    print(f"\n{'═'*60}")
    print(f"V6 HYBRID BENCHMARK (DPLL + Greedy)")
    print(f"{'═'*60}")
    print(f"{'R':>3} | {'DPLL':>5} | {'Kern':>5} | {'MQuad':>5} | {'DOF':>4} | {'Success':>7}")
    print(f"{'-'*3}-+-{'-'*5}-+-{'-'*5}-+-{'-'*5}-+-{'-'*4}-+-{'-'*7}")

    for R in range(2, max_R + 1):
        successes = 0
        dpll_res = kern = mq = 0
        for s in range(seeds):
            r = v6_hybrid(R, max_steps=1500, seed=s * 100 + 42, verbose=False)
            if r.get('success'):
                successes += 1
            dpll_res = r.get('dpll_resolved', 0)
            kern = r.get('kernel_dim', 0)
            mq = r.get('remaining_quad', 0)
        dof = kern - mq
        print(f"{R:>3} | {dpll_res:>5} | {kern:>5} | {mq:>5} | {dof:>4} | {successes}/{seeds}")
