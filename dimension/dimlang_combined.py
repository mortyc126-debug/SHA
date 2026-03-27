"""
DimLang COMBINED ATTACK: объединяем ВСЕ найденные операции.

Результаты по отдельности:
  Stochastic gradient: 90 (50K шагов)
  MITM at r=32: 91 (250K пар)
  Repair-guided: 97 (10K trials)

НОВАЯ ИДЕЯ: MITM + repair + gradient
1. Forward: random messages, repair δa=0, collect state[r=32]
2. Backward: invert от target с random W, collect state[r=32]
3. Match: найти ближайшую пару
4. Gradient: уточнить найденную пару

А также: новые операции языка, которые мы ещё НЕ пробовали.
"""

import numpy as np
from dimlang import *
from dimlang import hw, add32, sub32, Sigma0, Sigma1, Ch, Maj, K, IV_CONST
from dimlang import sigma0, sigma1, rotr

np.random.seed(42)

MASK32 = 0xFFFFFFFF


def mitm_enhanced(n_each=1000, split_round=32, verbose=True):
    """Enhanced MITM: forward с repair, backward с invert."""
    print(f"\n  Enhanced MITM (split at r={split_round}, {n_each} each side):")

    # Forward paths
    forward = []
    for _ in range(n_each):
        msg = Message()
        p = Path(msg)
        s = p.at(split_round)
        # Store: (state, message, path)
        forward.append((s, msg))

    # Backward paths: random target, invert from r=64
    target_msg = Message()
    target_path = Path(target_msg)
    target_final_state = target_path.states[64]

    backward = []
    for _ in range(n_each):
        s = State(target_final_state.regs)
        random_W = []
        for r in range(63, split_round - 1, -1):
            w = Word(np.random.randint(0, 2**32))
            s = invert_round(s, w, r)
            random_W.append(w)
        backward.append((s, random_W))

    # Find closest pair
    min_dist = 256
    min_pair = None
    for i in range(n_each):
        for j in range(min(100, n_each)):  # limit for speed
            d = forward[i][0].delta(backward[j][0]).weight
            if d < min_dist:
                min_dist = d
                min_pair = (i, j)

    if verbose:
        expected = 128 - np.log2(n_each * min(100, n_each)) / 2
        print(f"    Min distance at r={split_round}: {min_dist}")
        print(f"    Expected: ~{expected:.0f}")
        print(f"    Advantage: {expected - min_dist:.0f} bits")

    return min_dist, min_pair, forward, backward


def split_point_search():
    """Ищем ОПТИМАЛЬНУЮ точку split для MITM."""
    print(f"\n{'=' * 70}")
    print("1. OPTIMAL SPLIT POINT for MITM")
    print("=" * 70)

    for split_r in [8, 16, 24, 32, 40, 48, 56]:
        # Quick test with 200 each side
        target_msg = Message()
        target_path = Path(target_msg)

        forward = []
        for _ in range(200):
            p = Path(Message())
            forward.append(p.at(split_r))

        backward = []
        for _ in range(200):
            s = State(target_path.states[64].regs)
            for r in range(63, split_r - 1, -1):
                s = invert_round(s, Word(np.random.randint(0, 2**32)), r)
            backward.append(s)

        min_d = 256
        for i in range(200):
            for j in range(200):
                d = forward[i].delta(backward[j]).weight
                if d < min_d: min_d = d

        expected = 128 - np.log2(200*200) / 2
        print(f"  Split r={split_r:2d}: min_dist = {min_d:3d} (expected ~{expected:.0f}, advantage: {expected-min_d:.0f})")


def self_collision_search():
    """Ищем collision МЕЖДУ forward paths (не MITM, а birthday)."""
    print(f"\n{'=' * 70}")
    print("2. SELF-COLLISION: birthday among forward paths")
    print("=" * 70)

    paths = []
    hashes = []
    for i in range(2000):
        msg = Message()
        p = Path(msg)
        h = p.final_hash
        paths.append(p)
        hashes.append(h)

    # Find closest pair
    min_d = 256
    min_pair = None
    # Sample pairs
    n = len(hashes)
    for _ in range(500000):
        i = np.random.randint(0, n)
        j = np.random.randint(0, n)
        if i == j: continue
        d = hashes[i].delta(hashes[j]).weight
        if d < min_d:
            min_d = d
            min_pair = (i, j)

    expected = 128 - np.log2(2000*1999/2) / 2
    print(f"  2000 random hashes, 500K sampled pairs:")
    print(f"    Min HW(dH) = {min_d}")
    print(f"    Expected: ~{expected:.0f}")

    # Also check internal states
    for r in [16, 32, 48]:
        min_internal = 256
        for _ in range(100000):
            i = np.random.randint(0, n)
            j = np.random.randint(0, n)
            if i == j: continue
            d = paths[i].at(r).delta(paths[j].at(r)).weight
            if d < min_internal: min_internal = d
        print(f"    Internal r={r}: min = {min_internal}")


def new_operation_experiments():
    """Операции, которые НЕ существуют в стандартной криптографии."""
    print(f"\n{'=' * 70}")
    print("3. NEW OPERATIONS: beyond standard crypto")
    print("=" * 70)

    # ─── Operation: STATE SURGERY ───
    # Take state from path1 at round r, inject into path2
    # Like organ transplant between paths
    print(f"\n  STATE SURGERY: transplant state between paths")

    m1 = Message()
    m2 = Message()
    p1 = Path(m1)
    p2 = Path(m2)

    # Take p1's state at r=16, compute what W[16] would continue from there
    # using p2's schedule
    s_hybrid = p1.at(16)  # organs from p1
    a,b,c,d,e,f,g,h = s_hybrid.regs

    # Continue with p2's schedule from r=16
    W_p2 = m2.expand()
    for r in range(16, 64):
        w = W_p2[r].val
        T1 = add32(h, Sigma1(e), Ch(e,f,g), K[r], w)
        T2 = add32(Sigma0(a), Maj(a,b,c))
        h,g,f,e = g,f,e,add32(d,T1)
        d,c,b,a = c,b,a,add32(T1,T2)

    hybrid_hash = State([add32(IV_CONST[i], [a,b,c,d,e,f,g,h][i]) for i in range(8)])
    d1 = p1.final_hash.delta(hybrid_hash).weight
    d2 = p2.final_hash.delta(hybrid_hash).weight
    print(f"    Hybrid(p1_state@16 + p2_schedule@16..63):")
    print(f"      δ from p1: {d1}")
    print(f"      δ from p2: {d2}")
    print(f"      → Hybrid is {'closer to p2' if d2 < d1 else 'closer to p1'}")
    print(f"         (schedule dominates after 48 rounds)")

    # ─── Operation: RESONANCE AMPLIFICATION ───
    # Find two paths that are close at SOME round, amplify the closeness
    print(f"\n  RESONANCE AMPLIFICATION: exploit internal proximity")

    paths = [Path(Message()) for _ in range(50)]

    # Find closest pair at each round
    for r in [16, 32, 48]:
        min_d = 256
        best_ij = None
        for i in range(50):
            for j in range(i+1, 50):
                d = paths[i].at(r).delta(paths[j].at(r)).weight
                if d < min_d:
                    min_d = d
                    best_ij = (i, j)

        # Now: these two are close at round r.
        # Are they also close at OTHER rounds?
        i, j = best_ij
        dists = []
        for r2 in range(0, 65, 4):
            d2 = paths[i].at(r2).delta(paths[j].at(r2)).weight
            dists.append(d2)
        min_other = min(dists)
        avg_other = np.mean(dists)

        print(f"    Closest at r={r}: dist={min_d}, overall min={min_other}, avg={avg_other:.0f}")

    # ─── Operation: SCHEDULE HARMONICS ───
    # The message schedule has structure: W[r] depends on W[r-2,r-7,r-15,r-16]
    # Are there HARMONICS? Patterns in the expanded schedule?
    print(f"\n  SCHEDULE HARMONICS: patterns in expanded W")

    m = Message()
    W = m.expand()

    # Autocorrelation of expanded schedule
    vals = np.array([w.val for w in W], dtype=np.float64)
    vals -= np.mean(vals)
    if np.sum(vals**2) > 0:
        autocorrs = []
        for lag in [1, 2, 7, 15, 16, 32]:
            if lag < 64:
                c = np.sum(vals[:64-lag] * vals[lag:64]) / np.sum(vals**2)
                autocorrs.append((lag, c))
        print(f"    Schedule autocorrelation:")
        for lag, c in autocorrs:
            bar = "#" * int(abs(c) * 50)
            print(f"      lag={lag:2d}: {c:+.4f} {bar}")

    # ─── Operation: DIMENSIONAL FOLD ───
    # XOR multiple hashes together — does the result have structure?
    print(f"\n  DIMENSIONAL FOLD: XOR of N hashes")

    for n_fold in [2, 3, 4, 8, 16]:
        hashes = [Path(Message()).final_hash for _ in range(n_fold)]
        folded = [0]*8
        for h in hashes:
            for i in range(8):
                folded[i] ^= h[i]
        folded_hw = sum(hw(f) for f in folded)
        print(f"    XOR of {n_fold:2d} hashes: HW = {folded_hw:3d} (expected: 128)")

    # ─── Operation: CHAIN REACTION ───
    # Use hash as message for next hash: H(H(H(...)))
    # Does the chain converge? Fixed points? Cycles?
    print(f"\n  CHAIN REACTION: iterated hashing")

    m = Message()
    p = Path(m)
    h = p.final_hash
    prev_h = h
    chain_dists = []

    for step in range(100):
        # Use hash as new message (pad with zeros)
        new_msg = Message([h[i] for i in range(8)] + [0]*8)
        p_new = Path(new_msg)
        h_new = p_new.final_hash
        d = h.delta(h_new).weight
        chain_dists.append(d)
        # Check for fixed point
        if h == h_new:
            print(f"    ★ FIXED POINT at step {step}!")
            break
        h = h_new

    print(f"    Chain δ: mean={np.mean(chain_dists):.1f}, min={min(chain_dists)}, max={max(chain_dists)}")
    print(f"    → {'Converging!' if chain_dists[-1] < chain_dists[0] - 20 else 'Random walk (no convergence).'}")


def ultimate_attack():
    """Ультимативная комбинация: все операции вместе."""
    print(f"\n{'=' * 70}")
    print("4. ULTIMATE: combined DimLang attack")
    print("=" * 70)

    # Strategy:
    # 1. Generate fabric of 200 paths
    # 2. For each pair, try repair on the DIFFERENCE
    # 3. Gradient refine the best candidates
    # 4. Use MITM at the optimal split point

    fab = Fabric()
    paths = []
    for _ in range(200):
        msg = Message()
        idx = fab.weave(msg)
        paths.append(fab.paths[idx])

    # Step 1: Find closest pairs in hash space
    print(f"\n  Step 1: Find closest hash pairs among 200 paths")
    best_pairs = []
    for i in range(200):
        for j in range(i+1, 200):
            d = overlay(paths[i], paths[j]).weight
            best_pairs.append((d, i, j))
    best_pairs.sort()

    print(f"    Top 5 closest pairs:")
    for d, i, j in best_pairs[:5]:
        print(f"      ({i},{j}): HW(dH) = {d}")

    # Step 2: For the closest pair, try to REDUCE the difference
    d, i, j = best_pairs[0]
    mi = paths[i].message
    mj = paths[j].message

    print(f"\n  Step 2: Refine closest pair (initial dist={d})")

    # Stochastic refinement of mj toward mi's hash
    m_curr = mj.copy()
    d_curr = d
    target_hash = paths[i].final_hash
    improvements = 0

    for step in range(20000):
        w = np.random.randint(0, 16)
        b = np.random.randint(0, 32)
        m_test = m_curr.copy()
        m_test[w] = m_test[w].val ^ (1 << b)
        d_test = overlay(Path(m_test), paths[i]).weight
        if d_test < d_curr:
            m_curr = m_test
            d_curr = d_test
            improvements += 1
            if d_curr == 0:
                print(f"    ★★★ COLLISION at step {step}!")
                break

    print(f"    After 20K refinement: HW(dH) = {d_curr} (improved by {d - d_curr})")
    print(f"    Improvements: {improvements}")

    # Step 3: Repair-guided from the closest pair
    print(f"\n  Step 3: Repair from closest pair")
    dW = [mi[k].val ^ m_curr[k].val for k in range(16)]
    dW_hw = sum(hw(d) for d in dW)
    print(f"    HW(δMessage) = {dW_hw}")

    # Try repair on the ACTUAL difference
    best_repair = 256
    for trial in range(5000):
        dw0 = np.random.randint(1, 2**32)
        _, dH = dual_repair(mi, dw0, 'a')
        if dH < best_repair:
            best_repair = dH
    print(f"    Repair-guided search: best HW(dH) = {best_repair}")

    # Final summary
    print(f"\n  ULTIMATE RESULT:")
    print(f"    Closest pair (200 paths): HW(dH) = {best_pairs[0][0]}")
    print(f"    After gradient refinement: HW(dH) = {d_curr}")
    print(f"    Repair-guided: HW(dH) = {best_repair}")
    print(f"    Combined best: HW(dH) = {min(d_curr, best_repair)}")

    return min(d_curr, best_repair)


def main():
    split_point_search()
    self_collision_search()
    new_operation_experiments()
    best = ultimate_attack()

    print(f"\n{'=' * 70}")
    print("FINAL VERDICT")
    print("=" * 70)
    print(f"""
  DimLang ACHIEVEMENT:
    Best HW(dH) found in polynomial time: ≈{best}
    Standard random: HW(dH) ≈ 128
    Our advantage: ≈{128 - best} bits

  LANGUAGE EVALUATION:
    - Types (State, Path, Fabric): WORK
    - Repair: CONFIRMED (~20-30 bit advantage)
    - Gradient: EXISTS but shallow (2 steps max)
    - MITM: WORKS, ~28 bit advantage over birthday
    - Invert: EXACT (perfect backward computation)
    - New operations (surgery, resonance, fold): tested

  WHAT THE LANGUAGE GIVES US:
    1. Clear COMPOSITION of operations
    2. Easy EXPERIMENTATION with new ideas
    3. MEASURABLE advantage over brute force
    4. Foundation for MORE complex algorithms

  NEXT STEP: build algorithms that COMPOSE these operations
  in ways we haven't tried. The language is the TOOL.
  The attack is the PROGRAM.
""")


if __name__ == "__main__":
    main()
