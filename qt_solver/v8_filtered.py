"""
V8 Filtered Kernel Solver.

Key discovery: for R=6, only 192/512 kernel vectors (38%) change
relevant message bits W[0..R-1]. The other 320 are junk (change
W[R..15] which don't enter any round).

By filtering to relevant-only kernel vectors:
  - Every move matters (no wasted evaluations)
  - Search space is 192-dimensional (not 512)
  - Greedy/SA explores the RIGHT subspace
"""

import random
import time
import math

from qt_solver.sha256_traced import (
    MASK32, sha256_compress, sha256_compress_traced, get_all_carry_chains,
)
from qt_solver.gf2 import gf2_solve
from qt_solver.v3_solver import build_v3_system, _iv_constant_value


def _extract_msg(bv, vm):
    msg = []
    for w in range(16):
        word = 0
        for b in range(32):
            if (bv >> vm.w(w, b)) & 1:
                word |= (1 << b)
        msg.append(word)
    return msg


def v8_filtered(num_rounds, max_steps=3000, seed=42, verbose=True):
    rng = random.Random(seed)
    msg = [rng.randint(0, MASK32) for _ in range(16)]

    if verbose:
        print(f"\n{'='*60}")
        print(f"V8 Filtered Kernel: R={num_rounds}")
        print(f"{'='*60}")

    lr, rq, system, trace, info = build_v3_system(msg, num_rounds)
    vm = system.vm
    n = info['n']
    target = trace.hash_words
    ref_carries = get_all_carry_chains(trace)

    result = gf2_solve(lr, n)
    if result is None:
        return {'success': False}
    particular, kernel = result

    # Filter kernel: keep only vectors that change W[0..R-1]
    relevant_mask = 0
    for w in range(num_rounds):
        for b in range(32):
            relevant_mask |= (1 << vm.w(w, b))

    filtered = [k for k in kernel if (k & relevant_mask) != 0]

    if verbose:
        pct = 100*len(filtered)//max(1,len(kernel))
        print(f"Total kernel: {len(kernel)}, Filtered: {len(filtered)} ({pct}%)")
        print(f"Remaining quad: {len(rq)}, Eff DOF: {len(filtered) - len(rq)}")

    if not filtered or not rq:
        # Fully linear or no relevant moves — verify directly
        found_msg = _extract_msg(particular, vm)
        found_hash = sha256_compress(found_msg, num_rounds)
        success = found_hash == target and found_msg != msg
        if verbose:
            print(f"Direct solve: hash_match={found_hash == target}")
        return {'success': success, 'filtered_kern': len(filtered), 'remaining_quad': len(rq)}

    def evaluate(bv):
        m = _extract_msg(bv, vm)
        t = sha256_compress_traced(m, num_rounds)
        c = get_all_carry_chains(t)
        cmm = sum(bin(ref_carries[i] ^ c[i]).count('1') for i in range(len(ref_carries)))
        assignment = {v: (bv >> v) & 1 for v in range(n)}
        qv = 0
        for lin_terms, quad_terms, const in rq:
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
        print(f"Start: cost={cost} (carry={cmm}, quad={qviol})")

    # SA with filtered kernel
    temperature = 120.0
    cooling = 0.998
    t0 = time.time()

    for step in range(max_steps):
        improved = False
        trials = min(len(filtered), 100)
        indices = rng.sample(range(len(filtered)), trials)

        for ki in indices:
            kvec = filtered[ki]
            cand = current ^ kvec
            c, cm, qv = evaluate(cand)
            if c < cost:
                cost, cmm, qviol = c, cm, qv
                current = cand
                improved = True
                if c < best_cost:
                    best_cost = c
                    best_bv = cand
            elif temperature > 0.1:
                delta = c - cost
                if delta < temperature and rng.random() < math.exp(-delta / temperature):
                    cost, cmm, qviol = c, cm, qv
                    current = cand
                    improved = True

        temperature *= cooling

        if not improved:
            nf = rng.randint(1, max(2, min(8, len(filtered))))
            for _ in range(nf):
                current ^= filtered[rng.randint(0, len(filtered) - 1)]
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
                    'filtered_kern': len(filtered), 'remaining_quad': len(rq)}

        if verbose and step % 500 == 0:
            print(f"  step {step}: cost={cost} (c={cmm},q={qviol}), "
                  f"best={best_cost}, T={temperature:.1f}, {time.time()-t0:.1f}s")

    elapsed = time.time() - t0
    if verbose:
        print(f"\nBest: cost={best_cost} ({elapsed:.1f}s)")
    return {'success': False, 'best_cost': best_cost,
            'filtered_kern': len(filtered), 'remaining_quad': len(rq)}


def benchmark_v8(max_R=8, seeds=5):
    print(f"\n{'═'*60}")
    print(f"V8 FILTERED KERNEL BENCHMARK")
    print(f"{'═'*60}")
    print(f"{'R':>3} | {'Filt':>5} | {'RemQ':>4} | {'DOF':>4} | {'Success':>7}")
    print(f"{'-'*3}-+-{'-'*5}-+-{'-'*4}-+-{'-'*4}-+-{'-'*7}")

    for R in range(2, max_R + 1):
        successes = 0
        filt = remq = 0
        for s in range(seeds):
            r = v8_filtered(R, max_steps=3000, seed=s * 100 + 42, verbose=False)
            if r.get('success'):
                successes += 1
            filt = r.get('filtered_kern', 0)
            remq = r.get('remaining_quad', 0)
        dof = filt - remq
        print(f"{R:>3} | {filt:>5} | {remq:>4} | {dof:>4} | {successes}/{seeds}")
