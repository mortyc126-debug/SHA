"""
DimLang GRADIENT ATTACK: используем тот факт, что gradient = 33 бита за 1 бит.
Поверхность hash НЕ плоская. Есть направления.

Стандартный gradient descent застрял на 1 шаге.
НО: мы не ограничены стандартным descent.

НОВЫЕ АЛГОРИТМЫ:
1. Multi-bit gradient: менять 2-3 бита одновременно
2. Stochastic gradient: случайные шаги в направлении снижения
3. REPAIR-GUIDED gradient: repair определяет 15 из 16 слов, gradient ищет δW[0]
4. FABRIC gradient: параллельный поиск из многих стартовых точек
5. RESONANCE gradient: использовать БЛИЗОСТЬ между путями как трамплин
"""

import numpy as np
from dimlang import *
from dimlang import hw

np.random.seed(42)


def gradient_full(m_start, m_target, max_steps=100, verbose=True):
    """Full gradient descent: пробуем ВСЕ 512 однобитных flip, выбираем лучший."""
    m_current = m_start.copy()
    p_target = Path(m_target)
    d_current = overlay(Path(m_current), p_target).weight
    history = [d_current]

    for step in range(max_steps):
        best_d = d_current
        best_move = None

        for word in range(16):
            for bit in range(32):
                m_test = m_current.copy()
                m_test[word] = m_test[word].val ^ (1 << bit)
                d_test = overlay(Path(m_test), p_target).weight
                if d_test < best_d:
                    best_d = d_test
                    best_move = (word, bit)

        if best_move is None:
            break

        m_current[best_move[0]] = m_current[best_move[0]].val ^ (1 << best_move[1])
        d_current = best_d
        history.append(d_current)

    if verbose:
        print(f"    Steps: {len(history)-1}, final dist: {d_current}")
        print(f"    History: {history}")
    return m_current, d_current, history


def stochastic_gradient(m_start, m_target, max_steps=10000, verbose=True):
    """Stochastic: случайный бит, принимаем если лучше (hill climbing)."""
    m_current = m_start.copy()
    p_target = Path(m_target)
    d_current = overlay(Path(m_current), p_target).weight
    best_d = d_current
    best_m = m_current.copy()
    history = [d_current]
    improvements = 0

    for step in range(max_steps):
        word = np.random.randint(0, 16)
        bit = np.random.randint(0, 32)
        m_test = m_current.copy()
        m_test[word] = m_test[word].val ^ (1 << bit)
        d_test = overlay(Path(m_test), p_target).weight

        if d_test < d_current:
            m_current = m_test
            d_current = d_test
            improvements += 1
            if d_current < best_d:
                best_d = d_current
                best_m = m_current.copy()
                history.append(d_current)
                if d_current == 0:
                    break

        # Simulated annealing: sometimes accept WORSE moves
        elif np.random.random() < 0.01:
            m_current = m_test
            d_current = d_test

    if verbose:
        print(f"    Steps: {max_steps}, improvements: {improvements}")
        print(f"    Best dist: {best_d}")
        if len(history) > 1:
            print(f"    Progress: {history[:15]}...")
    return best_m, best_d, history


def multi_bit_gradient(m_start, m_target, n_bits=2, trials=5000, verbose=True):
    """Multi-bit: flip N bits randomly, keep if improves."""
    m_current = m_start.copy()
    p_target = Path(m_target)
    d_current = overlay(Path(m_current), p_target).weight
    best_d = d_current
    improvements = 0

    for _ in range(trials):
        m_test = m_current.copy()
        for __ in range(n_bits):
            w = np.random.randint(0, 16)
            b = np.random.randint(0, 32)
            m_test[w] = m_test[w].val ^ (1 << b)

        d_test = overlay(Path(m_test), p_target).weight
        if d_test < d_current:
            m_current = m_test
            d_current = d_test
            improvements += 1
            if d_current < best_d:
                best_d = d_current

    if verbose:
        print(f"    Trials: {trials}, improvements: {improvements}, best: {best_d}")
    return m_current, best_d


def repair_guided_gradient(m_base, delta_w0_range, strategy='a', trials=10000, verbose=True):
    """Repair-guided: repair fixes 15 words, gradient searches δW[0]."""
    best_dH = 256
    best_dw0 = 0

    for _ in range(trials):
        dw0 = np.random.randint(1, 2**32)
        _, dH = dual_repair(m_base, dw0, strategy)
        if dH < best_dH:
            best_dH = dH
            best_dw0 = dw0

    if verbose:
        print(f"    Trials: {trials}, best HW(dH): {best_dH}, best dW[0]: {hex(best_dw0)}")
    return best_dw0, best_dH


def fabric_parallel_search(n_starts, m_target, steps_per_start=1000, verbose=True):
    """Parallel: много стартовых точек, каждая делает stochastic gradient."""
    p_target = Path(m_target)
    best_overall = 256
    best_m = None

    for i in range(n_starts):
        m_start = Message()
        d = overlay(Path(m_start), p_target).weight

        # Quick stochastic descent
        m_curr = m_start.copy()
        d_curr = d
        for _ in range(steps_per_start):
            w = np.random.randint(0, 16)
            b = np.random.randint(0, 32)
            m_test = m_curr.copy()
            m_test[w] = m_test[w].val ^ (1 << b)
            d_test = overlay(Path(m_test), p_target).weight
            if d_test < d_curr:
                m_curr = m_test
                d_curr = d_test

        if d_curr < best_overall:
            best_overall = d_curr
            best_m = m_curr

    if verbose:
        print(f"    Starts: {n_starts}, steps each: {steps_per_start}")
        print(f"    Best overall dist: {best_overall}")
    return best_m, best_overall


def meet_in_middle(m_target, n_forward=200, n_backward=200, verbose=True):
    """Meet-in-the-middle using DimLang operations.
    Forward: random messages → compute Path, collect state at r=32
    Backward: from target hash, invert 32 rounds with random W[32..63]
    Find: closest pair between forward and backward states at r=32
    """
    p_target = Path(m_target)

    # Forward paths
    forward_states = []  # (state_at_32, message)
    for _ in range(n_forward):
        msg = Message()
        p = Path(msg)
        forward_states.append((p.at(32), msg))

    # Backward paths: invert from target's final state
    # We need the state BEFORE feedforward: state[64] = final_state - IV
    final_state = p_target.states[64]  # before feedforward

    backward_states = []
    for _ in range(n_backward):
        s = State(final_state.regs)
        # Generate random W[32..63]
        random_W = [Word(np.random.randint(0, 2**32)) for _ in range(32)]
        for r in range(63, 31, -1):
            w = random_W[r - 32]
            s = invert_round(s, w, r)
        backward_states.append((s, random_W))

    # Find closest pair
    min_dist = 256
    min_pair = None
    for i, (sf, _) in enumerate(forward_states):
        for j, (sb, _) in enumerate(backward_states):
            d = sf.delta(sb).weight
            if d < min_dist:
                min_dist = d
                min_pair = (i, j)

    if verbose:
        print(f"    Forward: {n_forward}, Backward: {n_backward}")
        print(f"    Min distance at r=32: {min_dist}")
        expected = 128 - np.log2(n_forward * n_backward) / 2
        print(f"    Expected (birthday on 256-bit state): ~{expected:.0f}")
        print(f"    → {'MITM ADVANTAGE!' if min_dist < expected - 10 else 'No MITM advantage.'}")
    return min_dist, min_pair


def main():
    print("=" * 70)
    print("DimLang GRADIENT ALGORITHMS")
    print("=" * 70)

    m_target = Message()
    m_start = Message()

    # =========================================================
    print(f"\n{'=' * 70}")
    print("1. FULL GRADIENT DESCENT")
    print("=" * 70)
    _, d, hist = gradient_full(m_start, m_target, max_steps=20)
    print(f"  Gradient is {'exploitable' if len(hist) > 3 else 'flat (1-2 steps max)'}")

    # =========================================================
    print(f"\n{'=' * 70}")
    print("2. STOCHASTIC GRADIENT (hill climbing + annealing)")
    print("=" * 70)
    _, d_stoch, hist_stoch = stochastic_gradient(m_start, m_target, max_steps=50000)

    # =========================================================
    print(f"\n{'=' * 70}")
    print("3. MULTI-BIT GRADIENT (2-bit, 3-bit, 4-bit flips)")
    print("=" * 70)
    for n_bits in [1, 2, 3, 4, 8]:
        print(f"  {n_bits}-bit flips:")
        _, d_multi = multi_bit_gradient(m_start, m_target, n_bits=n_bits, trials=5000)

    # =========================================================
    print(f"\n{'=' * 70}")
    print("4. REPAIR-GUIDED SEARCH")
    print("=" * 70)
    m_base = Message()
    print(f"  Strategy 'a':")
    _, dH_a = repair_guided_gradient(m_base, None, 'a', trials=10000)
    print(f"  Strategy 'e':")
    _, dH_e = repair_guided_gradient(m_base, None, 'e', trials=10000)
    print(f"  Strategy 'smart':")
    _, dH_s = repair_guided_gradient(m_base, None, 'smart', trials=10000)

    # =========================================================
    print(f"\n{'=' * 70}")
    print("5. FABRIC PARALLEL SEARCH")
    print("=" * 70)
    _, d_fab = fabric_parallel_search(20, m_target, steps_per_start=2000)

    # =========================================================
    print(f"\n{'=' * 70}")
    print("6. MEET-IN-THE-MIDDLE at r=32")
    print("=" * 70)
    min_dist, _ = meet_in_middle(m_target, n_forward=500, n_backward=500)

    # =========================================================
    print(f"\n{'=' * 70}")
    print("7. COMBINED: repair + gradient + fabric")
    print("=" * 70)

    # Best strategy: use repair to get close, then gradient to optimize
    best_overall = 256
    for trial in range(500):
        m_trial = Message()
        dw0 = np.random.randint(1, 2**32)
        m_repaired, dH = dual_repair(m_trial, dw0, 'e')
        if dH < best_overall:
            best_overall = dH

        # Now: gradient on the repaired message pair
        # Can we improve the repaired pair?
        if dH < 115:  # only refine promising candidates
            # Try flipping non-repair words? We used W[0]..W[15] for repair.
            # No free words left! Repair consumed all.
            pass

    print(f"  500 repair trials, best HW(dH): {best_overall}")

    # =========================================================
    print(f"\n{'=' * 70}")
    print("8. ANALYSIS")
    print("=" * 70)

    print(f"""
  ALGORITHM COMPARISON:
    Full gradient:        ~2 steps, then flat
    Stochastic (50K):     best ≈ {d_stoch}
    Multi-bit (4-bit):    fast but no better
    Repair-guided (10K):  best ≈ {min(dH_a, dH_e, dH_s)}
    Fabric parallel:      best ≈ {d_fab}
    MITM at r=32:         min_dist = {min_dist}

  KEY FINDINGS:

  1. GRADIENT EXISTS but is VERY SHALLOW.
     One bit flip can improve by up to 33 bits.
     But second step: only ~5 bits improvement.
     Third step: 0. Surface is almost flat after 2 moves.

  2. STOCHASTIC BEATS DETERMINISTIC.
     Random search + hill climbing does better than exhaustive
     single-bit gradient. Because: landscape has many local minima.

  3. REPAIR GIVES CONSISTENT EDGE.
     ~20-25 bits below random 128 consistently.
     But: repair consumes ALL 16 free words.
     No degrees of freedom left for further optimization.

  4. THE 2^128 BARRIER:
     All algorithms converge to HW(dH) ≈ 100-110.
     Never below ~95.
     This is ~20-30 bits below random.
     But for COLLISION need dH = 0 (256 bits more).
     Gap: ~95-100 bits. Still astronomical.

  WHAT THIS MEANS FOR ATTACK:
     We can find messages with HW(dH) ≈ 100 in polynomial time.
     Standard birthday would give HW(dH) = 0 in 2^128 time.
     Our gradient gives 25-bit "head start".
     Effective cost: ~2^(100/2) = 2^50 with gradient???
     NO! HW(dH)=100 is NEAR-collision, not collision.
     To go from 100 to 0 still needs birthday on 100-bit space.
     Cost: 2^50. Total: polynomial + 2^50 = 2^50???

  WAIT. Let's verify this.
""")


def verify_near_collision_to_collision():
    """Can we use near-collisions as stepping stones?"""
    print(f"\n{'=' * 70}")
    print("9. NEAR-COLLISION TO COLLISION: stepping stone?")
    print("=" * 70)

    # If we can find pairs with HW(dH) = 100 easily,
    # and HW(dH) = 100 means 100 specific bits differ,
    # can we find a THIRD message that matches those 100 bits?

    # Multi-collision approach:
    # Find m1, m2 with small dH.
    # Then find m3 close to m1, m4 close to m2.
    # If m3 is closer to m4 than m1 is to m2...

    m1 = Message()
    p1 = Path(m1)

    # Find near-collision partner
    best_d = 256
    best_m2 = None
    for _ in range(1000):
        m2 = Message()
        d = overlay(p1, Path(m2)).weight
        if d < best_d:
            best_d = d
            best_m2 = m2

    print(f"  Near-collision: HW(dH) = {best_d}")

    # Now: can we REFINE this near-collision?
    # Try: perturb m2 to reduce dH
    m_curr = best_m2.copy()
    d_curr = best_d
    for _ in range(10000):
        w = np.random.randint(0, 16)
        b = np.random.randint(0, 32)
        m_test = m_curr.copy()
        m_test[w] = m_test[w].val ^ (1 << b)
        d_test = overlay(p1, Path(m_test)).weight
        if d_test < d_curr:
            m_curr = m_test
            d_curr = d_test
            if d_curr == 0:
                print(f"  ★★★ COLLISION FOUND!")
                break

    print(f"  After 10K refinement: HW(dH) = {d_curr}")
    print(f"  Improvement: {best_d - d_curr} bits")

    # Multi-target: find MANY near-collision pairs, then birthday among them
    print(f"\n  Multi-near-collision birthday:")
    pairs = []
    for _ in range(200):
        ma = Message()
        mb = Message()
        pa = Path(ma)
        pb = Path(mb)
        ha = pa.final_hash
        hb = pb.final_hash
        d = ha.delta(hb).weight
        pairs.append((d, ha, hb))

    pairs.sort(key=lambda x: x[0])
    print(f"  200 random pairs, 5 closest:")
    for d, ha, hb in pairs[:5]:
        print(f"    HW(dH) = {d}")

    # Birthday on hashes themselves
    hashes = {}
    for i in range(5000):
        m = Message()
        p = Path(m)
        h = p.final_hash
        # Truncate to first word for birthday
        key = h[0]
        if key in hashes:
            d = h.delta(hashes[key][1]).weight
            print(f"  32-bit birthday collision at i={i}: HW(full dH) = {d}")
            break
        hashes[key] = (m, h)
    else:
        print(f"  No 32-bit birthday in 5K (expected: ~65K needed)")


if __name__ == "__main__":
    main()
    verify_near_collision_to_collision()
