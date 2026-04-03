"""
V9 Adaptive Score Solver.

Key insight: 55/192 kernel vectors are IMPROVING at any given point.
But which 55 changes as we move. Pre-compute score per vector per step,
take the BEST improving move, repeat.

This is "smart greedy" — instead of trying random vectors and evaluating
SHA-256 (expensive), we compute an algebraic SCORE that predicts
improvement WITHOUT full evaluation.

Score(k) = #{violated eqs k would fix} - #{satisfied eqs k would break}
Positive score → guaranteed improvement of quad violations.
"""

import random
import time

from qt_solver.sha256_traced import (
    MASK32, sha256_compress, sha256_compress_traced, get_all_carry_chains,
)
from qt_solver.gf2 import gf2_solve
from qt_solver.v3_solver import build_v3_system


def _extract_msg(bv, vm):
    msg = []
    for w in range(16):
        word = 0
        for b in range(32):
            if (bv >> vm.w(w, b)) & 1:
                word |= (1 << b)
        msg.append(word)
    return msg


def v9_adaptive(num_rounds, max_steps=3000, seed=42, verbose=True):
    rng = random.Random(seed)
    msg = [rng.randint(0, MASK32) for _ in range(16)]

    if verbose:
        print(f"\n{'='*60}")
        print(f"V9 Adaptive Score: R={num_rounds}")
        print(f"{'='*60}")

    lr, rq, system, trace, info = build_v3_system(msg, num_rounds)
    vm = system.vm
    n = info['n']
    target = trace.hash_words
    ref_carries = get_all_carry_chains(trace)

    result = gf2_solve(lr, n)
    if not result:
        return {'success': False}
    particular, kernel = result

    # Filter to relevant kernel vectors
    rel_mask = 0
    for w in range(num_rounds):
        for b in range(32):
            rel_mask |= (1 << vm.w(w, b))
    filtered = [k for k in kernel if (k & rel_mask) != 0]

    if verbose:
        print(f"Filtered kernel: {len(filtered)}, Quad eqs: {len(rq)}")

    if not filtered:
        found = _extract_msg(particular, vm)
        h = sha256_compress(found, num_rounds)
        return {'success': h == target and found != msg}

    # Precompute: for each eq, which vars does it use?
    eq_var_masks = []
    for lin_terms, quad_terms, const in rq:
        mask = 0
        for v in lin_terms:
            mask |= (1 << v)
        for v1, v2 in quad_terms:
            mask |= (1 << v1)
            mask |= (1 << v2)
        eq_var_masks.append(mask)

    # For each kernel vector: which equations it affects
    kvec_eq_sets = []
    for kvec in filtered:
        affected = []
        for i, mask in enumerate(eq_var_masks):
            if kvec & mask:
                affected.append(i)
        kvec_eq_sets.append(affected)

    def eval_equations(bv):
        """Return list of 0/1 per equation (0=satisfied, 1=violated)."""
        assignment = {}
        for v in range(n):
            assignment[v] = (bv >> v) & 1
        results = []
        for lin_terms, quad_terms, const in rq:
            val = const
            for v in lin_terms:
                val ^= assignment.get(v, 0)
            for v1, v2 in quad_terms:
                val ^= (assignment.get(v1, 0) & assignment.get(v2, 0))
            results.append(val)
        return results

    def compute_scores(bv, eq_status):
        """
        For each kernel vector, compute score = fixes - breaks.
        Only checks AFFECTED equations (sparse!).
        """
        scores = []
        for ki, affected_eqs in enumerate(kvec_eq_sets):
            if not affected_eqs:
                scores.append(0)
                continue

            # Evaluate only affected equations for bv ^ kvec
            new_bv = bv ^ filtered[ki]
            new_assignment = {}
            # Only need vars that appear in affected eqs
            # For speed: compute full assignment lazily
            fix = 0
            brk = 0
            for ei in affected_eqs:
                old_val = eq_status[ei]
                # Compute new value
                lin_terms, quad_terms, const = rq[ei]
                new_val = const
                for v in lin_terms:
                    new_val ^= (new_bv >> v) & 1
                for v1, v2 in quad_terms:
                    new_val ^= ((new_bv >> v1) & 1) & ((new_bv >> v2) & 1)

                if old_val == 1 and new_val == 0:
                    fix += 1
                elif old_val == 0 and new_val == 1:
                    brk += 1

            scores.append(fix - brk)
        return scores

    def full_cost(bv):
        m = _extract_msg(bv, vm)
        t = sha256_compress_traced(m, num_rounds)
        c = get_all_carry_chains(t)
        cmm = sum(bin(ref_carries[i] ^ c[i]).count('1') for i in range(len(ref_carries)))
        eq_st = eval_equations(bv)
        qv = sum(eq_st)
        return cmm + qv * 5, cmm, qv, eq_st

    # Initialize
    current = particular
    cost, cmm, qviol, eq_status = full_cost(current)
    best_cost, best_bv = cost, current

    if verbose:
        print(f"Start: cost={cost} (carry={cmm}, quad={qviol})")

    t0 = time.time()
    score_evals = 0

    for step in range(max_steps):
        # Compute algebraic scores (CHEAP: only touches affected eqs)
        scores = compute_scores(current, eq_status)
        score_evals += 1

        # Find best positive-score vector
        best_idx = -1
        best_score = 0
        for i, s in enumerate(scores):
            if s > best_score:
                best_score = s
                best_idx = i

        if best_idx >= 0:
            # Apply improving move
            current ^= filtered[best_idx]
            cost, cmm, qviol, eq_status = full_cost(current)
            if cost < best_cost:
                best_cost = cost
                best_bv = current
        else:
            # No improving move found — random perturbation
            nf = rng.randint(2, min(8, len(filtered)))
            for _ in range(nf):
                current ^= filtered[rng.randint(0, len(filtered) - 1)]
            cost, cmm, qviol, eq_status = full_cost(current)

        if best_cost == 0:
            found_msg = _extract_msg(best_bv, vm)
            found_hash = sha256_compress(found_msg, num_rounds)
            success = found_hash == target and found_msg != msg
            if verbose:
                print(f"\n  cost=0 at step {step}! Hash: {found_hash == target}")
                if success:
                    print(f"  *** SECOND PREIMAGE FOUND! ***")
            return {'success': success, 'steps': step,
                    'score_evals': score_evals}

        if verbose and step % 500 == 0:
            improving = sum(1 for s in scores if s > 0)
            print(f"  step {step}: cost={cost} (c={cmm},q={qviol}), "
                  f"best={best_cost}, improving={improving}, "
                  f"{time.time()-t0:.1f}s")

    elapsed = time.time() - t0
    if verbose:
        print(f"\nBest: cost={best_cost} ({elapsed:.1f}s)")
    return {'success': False, 'best_cost': best_cost}


def benchmark_v9(max_R=8, seeds=5):
    print(f"\n{'═'*60}")
    print(f"V9 ADAPTIVE SCORE BENCHMARK")
    print(f"{'═'*60}")
    for R in range(2, max_R + 1):
        successes = 0
        for s in range(seeds):
            r = v9_adaptive(R, max_steps=3000, seed=s * 100 + 42, verbose=False)
            if r.get('success'):
                successes += 1
        print(f"R={R}: {successes}/{seeds}")
