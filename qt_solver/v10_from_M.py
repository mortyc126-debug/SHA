"""
V10: Start from M (cost=0), search for second preimage.

Key insight: M itself is in Q-variety at cost=0. Starting from M
instead of random particular (cost~1200) dramatically reduces distance.

Strategy: greedy walk from M, minimizing cost. If we find M' ≠ M
with cost=0 → second preimage. If we only find M → try perturbation.
"""

import random, time, math
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


def _encode_M(msg, trace, vm, R):
    """Encode actual message M as bitvec in Q-variety."""
    bv = 0
    for w in range(16):
        for b in range(32):
            if (msg[w] >> b) & 1:
                bv |= (1 << vm.w(w, b))
    for r in range(1, R + 1):
        for reg in range(8):
            for b in range(32):
                if (trace.states[r][reg] >> b) & 1:
                    bv |= (1 << vm.s(r, reg, b))
    for w in range(16, R):
        for b in range(32):
            if (trace.schedule[w] >> b) & 1:
                bv |= (1 << vm.w(w, b))
    return bv


def v10_from_M(num_rounds, max_steps=5000, seed=42, verbose=True):
    rng = random.Random(seed)
    msg = [rng.randint(0, MASK32) for _ in range(16)]

    if verbose:
        print(f"\n{'='*60}")
        print(f"V10 From-M Solver: R={num_rounds}")
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

    # Encode M
    M_bv = _encode_M(msg, trace, vm, num_rounds)

    # Filter kernel
    rel_mask = 0
    for w in range(num_rounds):
        for b in range(32):
            rel_mask |= (1 << vm.w(w, b))
    filtered = [k for k in kernel if (k & rel_mask) != 0]

    if verbose:
        print(f"Filtered kernel: {len(filtered)}, Quad eqs: {len(rq)}")

    def full_cost(bv):
        m = _extract_msg(bv, vm)
        t = sha256_compress_traced(m, num_rounds)
        c = get_all_carry_chains(t)
        cmm = sum(bin(ref_carries[i] ^ c[i]).count('1') for i in range(len(ref_carries)))
        assignment = {v: (bv >> v) & 1 for v in range(n)}
        qv = 0
        for lin_terms, quad_terms, const in rq:
            val = const
            for v in lin_terms: val ^= assignment.get(v, 0)
            for v1, v2 in quad_terms: val ^= (assignment.get(v1, 0) & assignment.get(v2, 0))
            if val != 0: qv += 1
        return cmm + qv * 5, cmm, qv

    if not filtered:
        return {'success': False}

    # Phase 1: FROM M, explore outward
    # M has cost=0. Search for another cost=0 point.
    # Strategy: 2-hop and 3-hop from M (combine kernel vectors)
    if verbose:
        print(f"Phase 1: Multi-hop search from M (cost=0)")

    t0 = time.time()
    best_cost = 9999
    best_bv = None

    # 1-hop
    for kvec in filtered:
        candidate = M_bv ^ kvec
        if candidate == M_bv:
            continue
        c, cm, qv = full_cost(candidate)
        if c < best_cost:
            best_cost = c
            best_bv = candidate
        if c == 0:
            found = _extract_msg(candidate, vm)
            if sha256_compress(found, num_rounds) == target and found != msg:
                if verbose:
                    print(f"  1-hop: SECOND PREIMAGE FOUND!")
                return {'success': True, 'hops': 1}

    if verbose:
        print(f"  1-hop best: {best_cost} ({time.time()-t0:.1f}s)")

    # 2-hop: combine pairs of kernel vectors
    # Smart: only combine vectors that individually gave low cost
    scored = []
    for i, kvec in enumerate(filtered):
        c, _, _ = full_cost(M_bv ^ kvec)
        scored.append((c, i))
    scored.sort()
    top_k = min(40, len(scored))
    top_indices = [idx for _, idx in scored[:top_k]]

    for i in range(top_k):
        for j in range(i + 1, top_k):
            combo = filtered[top_indices[i]] ^ filtered[top_indices[j]]
            candidate = M_bv ^ combo
            if candidate == M_bv:
                continue
            c, cm, qv = full_cost(candidate)
            if c < best_cost:
                best_cost = c
                best_bv = candidate
            if c == 0:
                found = _extract_msg(candidate, vm)
                if sha256_compress(found, num_rounds) == target and found != msg:
                    if verbose:
                        print(f"  2-hop: SECOND PREIMAGE FOUND!")
                    return {'success': True, 'hops': 2}

    if verbose:
        print(f"  2-hop best: {best_cost} ({time.time()-t0:.1f}s)")

    # 3-hop from top candidates
    top3 = min(15, len(scored))
    top3_indices = [idx for _, idx in scored[:top3]]
    for i in range(top3):
        for j in range(i + 1, top3):
            for k in range(j + 1, top3):
                combo = filtered[top3_indices[i]] ^ filtered[top3_indices[j]] ^ filtered[top3_indices[k]]
                candidate = M_bv ^ combo
                if candidate == M_bv:
                    continue
                c, cm, qv = full_cost(candidate)
                if c < best_cost:
                    best_cost = c
                    best_bv = candidate
                if c == 0:
                    found = _extract_msg(candidate, vm)
                    if sha256_compress(found, num_rounds) == target and found != msg:
                        if verbose:
                            print(f"  3-hop: SECOND PREIMAGE FOUND!")
                        return {'success': True, 'hops': 3}

    if verbose:
        print(f"  3-hop best: {best_cost} ({time.time()-t0:.1f}s)")

    # Phase 2: Greedy descent from best point found near M
    if verbose:
        print(f"\nPhase 2: Greedy descent from best={best_cost}")

    current = best_bv if best_bv else M_bv
    cost, cmm, qviol = full_cost(current)

    temperature = 60.0
    cooling = 0.998

    for step in range(max_steps):
        improved = False
        trials = min(len(filtered), 80)
        indices = rng.sample(range(len(filtered)), trials)

        for ki in indices:
            kvec = filtered[ki]
            cand = current ^ kvec
            c, cm, qv = full_cost(cand)
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
            nf = rng.randint(1, max(2, min(6, len(filtered))))
            for _ in range(nf):
                current ^= filtered[rng.randint(0, len(filtered) - 1)]
            cost, cmm, qviol = full_cost(current)

        if best_cost == 0:
            found = _extract_msg(best_bv, vm)
            fh = sha256_compress(found, num_rounds)
            success = fh == target and found != msg
            if verbose:
                print(f"  cost=0 at step {step}! Hash: {fh == target}")
                if success:
                    print(f"  *** SECOND PREIMAGE FOUND! ***")
            return {'success': success, 'steps': step, 'best_cost': 0}

        if verbose and step % 500 == 0:
            print(f"  step {step}: cost={cost} (c={cmm},q={qviol}), "
                  f"best={best_cost}, T={temperature:.1f}, {time.time()-t0:.1f}s")

    if verbose:
        print(f"\nBest: {best_cost}")
    return {'success': False, 'best_cost': best_cost}


def benchmark_v10(max_R=8, seeds=5):
    print(f"\n{'═'*60}")
    print(f"V10 FROM-M BENCHMARK")
    print(f"{'═'*60}")
    for R in range(2, max_R + 1):
        successes = 0
        best = 9999
        for s in range(seeds):
            r = v10_from_M(R, max_steps=3000, seed=s * 100 + 42, verbose=False)
            if r.get('success'): successes += 1
            bc = r.get('best_cost', 9999)
            if bc < best: best = bc
        print(f"R={R}: {successes}/{seeds}, best_cost={best}")
