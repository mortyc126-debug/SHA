"""
DimLang TEST SUITE: тестируем каждую операцию.
Ищем то, чего НЕ ДОЛЖНО БЫТЬ.
Границ нет — экспериментируем.
"""

import numpy as np
import hashlib, struct
from dimlang import *
from dimlang import rotr, hw, add32, sub32, Sigma0, Sigma1, Ch, Maj, K, IV_CONST

np.random.seed(42)


def test_basics():
    print("=" * 70)
    print("TEST 1: BASIC TYPES")
    print("=" * 70)

    s = State()
    print(f"  IV State: {s}")
    print(f"  a = {hex(s.a)}, e = {hex(s.e)}")

    m = Message()
    print(f"  Random message: {[hex(w.val) for w in m.words[:4]]}...")

    p = Path(m)
    print(f"  Path: {len(p.states)} states, {len(p.schedule)} schedule words")
    print(f"  State[0]: {p.at(0)}")
    print(f"  State[64]: {p.at(64)}")
    print(f"  Hash: {p.final_hash}")

    # Verify against hashlib
    raw = b''
    for w in m.words:
        raw += w.val.to_bytes(4, 'big')
    # pad manually
    raw += b'\x80' + b'\x00' * (55) + (512).to_bytes(8, 'big')
    h = hashlib.sha256(raw[:64]).digest()
    import struct
    h_words = struct.unpack('>8I', h)
    h_state = State(h_words)
    print(f"  hashlib: {h_state}")
    print(f"  Match: {p.final_hash == h_state}")


def test_overlay():
    print(f"\n{'=' * 70}")
    print("TEST 2: OVERLAY")
    print("=" * 70)

    m1 = Message()
    m2 = Message()
    p1 = Path(m1)
    p2 = Path(m2)

    d = overlay(p1, p2)
    print(f"  Random overlay: δ = {d}, weight = {d.weight}")

    # Same message
    d_same = overlay(p1, Path(Message(m1.raw())))
    print(f"  Same message: δ = {d_same}, weight = {d_same.weight}, zero = {d_same.is_zero}")

    # Per-round overlay
    print(f"\n  Per-round overlay (random vs random):")
    for r in [0, 1, 4, 8, 16, 32, 64]:
        if r <= 64:
            d_r = overlay(p1, p2, r)
            print(f"    r={r:2d}: weight = {d_r.weight}")


def test_repair():
    print(f"\n{'=' * 70}")
    print("TEST 3: REPAIR OPERATIONS")
    print("=" * 70)

    m = Message()

    for strategy in ['a', 'e', 'alternate', 'smart']:
        _, dH = dual_repair(m, 1, strategy)
        print(f"  Strategy '{strategy:>9}': HW(δH) = {dH}")

    # Try multiple initial differences
    print(f"\n  Varying δW[0] with strategy 'a':")
    for dw0 in [1, 0x80000000, 0xFFFFFFFF, 0x12345678]:
        _, dH = dual_repair(m, dw0, 'a')
        print(f"    δW[0]={hex(dw0):>12}: HW(δH) = {dH}")

    # Try many random differences, find BEST
    print(f"\n  Random δW[0] search (1000 trials):")
    best = 256
    best_dw = 0
    for _ in range(1000):
        dw = np.random.randint(1, 2**32)
        _, dH = dual_repair(m, dw, 'a')
        if dH < best:
            best = dH
            best_dw = dw
    print(f"    Best: HW(δH) = {best}, δW[0] = {hex(best_dw)}")


def test_invert():
    print(f"\n{'=' * 70}")
    print("TEST 4: INVERT (backward round)")
    print("=" * 70)

    m = Message()
    p = Path(m)

    # Invert round 63 → should get state[63] from state[64]
    s_recovered = invert_round(p.at(64), p.schedule[63], 63)
    s_actual = p.at(63)
    d = s_recovered.delta(s_actual)
    print(f"  Invert round 63: δ = {d.weight} (should be 0)")

    # Invert multiple rounds
    print(f"  Multi-round inversion:")
    s = p.at(64)
    for r in range(63, -1, -1):
        s = invert_round(s, p.schedule[r], r)
    d = s.delta(p.at(0))
    print(f"    Invert all 64 rounds: δ = {d.weight} (should be 0)")
    print(f"    Recovered IV: {s}")
    print(f"    Actual IV:    {p.at(0)}")
    print(f"    Match: {s == p.at(0)}")


def test_tunnel():
    print(f"\n{'=' * 70}")
    print("TEST 5: TUNNEL (skip rounds)")
    print("=" * 70)

    m = Message()
    p = Path(m)

    # Tunnel: δ between different rounds of SAME path
    print(f"  Self-tunnel (δ between rounds of same path):")
    for gap in [1, 2, 4, 8, 16, 32]:
        weights = []
        for start_r in range(0, 64 - gap, 4):
            d = tunnel(p, start_r, start_r + gap)
            weights.append(d.weight)
        print(f"    gap={gap:2d}: mean δ = {np.mean(weights):.1f}, std = {np.std(weights):.1f}")


def test_warp():
    print(f"\n{'=' * 70}")
    print("TEST 6: WARP (state deformation)")
    print("=" * 70)

    m = Message()
    p = Path(m)

    # Warp 1: bit rotation of all registers
    def warp_rotate(state, r):
        return State([rotr(state[i], 1) for i in range(8)])

    warped = warp(p, warp_rotate, range(64))
    # Compare warped states to original
    for r in [0, 16, 32, 64]:
        d = warped[r].delta(p.at(r))
        print(f"  Warp(rotate) r={r:2d}: δ = {d.weight}")

    # Warp 2: XOR all registers with constant
    def warp_xor_const(state, r):
        return State([state[i] ^ 0xDEADBEEF for i in range(8)])

    warped2 = warp(p, warp_xor_const, range(64))
    for r in [0, 16, 32, 64]:
        d = warped2[r].delta(p.at(r))
        print(f"  Warp(XOR const) r={r:2d}: δ = {d.weight}")

    # Warp 3: swap a↔e registers
    def warp_swap_ae(state, r):
        regs = list(state.regs)
        regs[0], regs[4] = regs[4], regs[0]
        return State(regs)

    warped3 = warp(p, warp_swap_ae, range(64))
    dists = []
    for r in range(1, 65):
        d = warped3[r].delta(p.at(r))
        dists.append(d.weight)
    print(f"\n  Warp(swap a<->e): mean δ = {np.mean(dists):.1f}")
    print(f"    This tells us: a and e are {'similar' if np.mean(dists) < 64 else 'different'}")


def test_resonate():
    print(f"\n{'=' * 70}")
    print("TEST 7: RESONATE (find closest paths)")
    print("=" * 70)

    # Create fabric with many paths
    paths = []
    for _ in range(100):
        paths.append(Path(Message()))

    # Find resonance at different rounds
    for r in [0, 16, 32, 64]:
        res = resonate(paths[:20], r)  # smaller set for speed
        print(f"  Resonate at r={r:2d}: min_dist = {res['min_dist']}, pair = {res['min_pair']}")

    # KEY QUESTION: is min_dist at r=64 (output) EVER less than expected?
    # Birthday: for 100 random 256-bit values, min dist ≈ 128 - log2(100*99/2) ≈ 128 - 13 ≈ 115
    print(f"\n  Full 100-path resonance at final hash:")
    res_full = resonate(paths, 64)
    print(f"    min_dist = {res_full['min_dist']}")
    print(f"    Expected (birthday): ≈ 115")

    # Now: resonance at INTERNAL rounds — is it different?
    print(f"\n  Internal resonance (20 paths):")
    internal_min = []
    for r in range(0, 65, 4):
        res = resonate(paths[:20], r)
        internal_min.append((r, res['min_dist']))
    for r, d in internal_min:
        bar = "#" * (d // 4)
        print(f"    r={r:2d}: min_dist = {d:3d} {bar}")


def test_new_operations():
    print(f"\n{'=' * 70}")
    print("TEST 8: NEW UNEXPLORED OPERATIONS")
    print("=" * 70)

    m = Message()
    p = Path(m)

    # Operation: FOLD two paths through XOR of their states
    print(f"  FOLD: XOR states of two paths at each round")
    m2 = Message()
    p2 = Path(m2)

    fold_states = []
    for r in range(65):
        s1 = p.at(r)
        s2 = p2.at(r)
        folded = State([s1[i] ^ s2[i] for i in range(8)])
        fold_states.append(folded)

    # Is folded state trajectory structured or random?
    fold_deltas = []
    for r in range(64):
        d = fold_states[r].delta(fold_states[r+1])
        fold_deltas.append(d.weight)
    print(f"    Fold trajectory δ(r,r+1): mean={np.mean(fold_deltas):.1f}, std={np.std(fold_deltas):.1f}")
    print(f"    Compare: single path δ(r,r+1): ", end="")
    single_deltas = [p.at(r).delta(p.at(r+1)).weight for r in range(64)]
    print(f"mean={np.mean(single_deltas):.1f}, std={np.std(single_deltas):.1f}")

    # Operation: COMPOSE two paths (use output of p1 as IV for p2)
    print(f"\n  COMPOSE: use hash(m1) as 'IV' for second computation")
    h1 = p.final_hash
    # Can we start from h1 as IV?
    a,b,c,d,e,f,g,h = h1.regs
    states_composed = [h1]
    for r in range(64):
        w = p2.schedule[r].val
        T1 = add32(h, Sigma1(e), Ch(e,f,g), K[r], w)
        T2 = add32(Sigma0(a), Maj(a,b,c))
        h,g,f,e = g,f,e,add32(d,T1)
        d,c,b,a = c,b,a,add32(T1,T2)
        states_composed.append(State((a,b,c,d,e,f,g,h)))
    composed_hash = State([add32(h1[i], states_composed[64][i]) for i in range(8)])
    print(f"    Composed hash: {composed_hash}")

    # Operation: MIRROR — run path BACKWARDS
    print(f"\n  MIRROR: run path backwards (invert all rounds)")
    s = p.at(64)
    mirror_states = [s]
    for r in range(63, -1, -1):
        s = invert_round(s, p.schedule[r], r)
        mirror_states.append(s)
    mirror_states.reverse()
    # Check
    print(f"    Mirror[0] == Path[0]: {mirror_states[0] == p.at(0)}")
    # What's the trajectory δ of mirror?
    mirror_deltas = [mirror_states[r].delta(mirror_states[r+1]).weight for r in range(64)]
    print(f"    Mirror δ(r,r+1): mean={np.mean(mirror_deltas):.1f}")
    print(f"    Forward δ(r,r+1): mean={np.mean(single_deltas):.1f}")
    print(f"    → Mirror = Forward? {np.allclose(mirror_deltas, single_deltas)}")

    # Operation: INTERLEAVE — alternate words from two messages
    print(f"\n  INTERLEAVE: m3[r] = m1[r] if r even, m2[r] if r odd")
    m3 = Message([m.words[r].val if r % 2 == 0 else m2.words[r].val for r in range(16)])
    p3 = Path(m3)
    d13 = overlay(p, p3)
    d23 = overlay(p2, p3)
    print(f"    δ(m1, interleaved) = {d13.weight}")
    print(f"    δ(m2, interleaved) = {d23.weight}")
    print(f"    Expected (random): ~128")

    # Operation: GRADIENT — which input bit reduces δH most?
    print(f"\n  GRADIENT: which W[i] bit MOST reduces δH toward target?")
    target = p2.final_hash  # arbitrary target
    d_base = p.final_hash.delta(target).weight
    print(f"    Base distance to target: {d_base}")

    best_improvement = 0
    best_bit = None
    for word in range(16):
        for bit in range(32):
            m_test = m.copy()
            m_test[word] = m_test[word].val ^ (1 << bit)
            p_test = Path(m_test)
            d_test = p_test.final_hash.delta(target).weight
            improvement = d_base - d_test
            if improvement > best_improvement:
                best_improvement = improvement
                best_bit = (word, bit, d_test)
    print(f"    Best single-bit improvement: {best_improvement} bits")
    if best_bit:
        print(f"    Bit: W[{best_bit[0]}] bit {best_bit[1]} → dist = {best_bit[2]}")
    print(f"    → {'GRADIENT WORKS!' if best_improvement > 10 else 'Gradient ≈ flat (random)'}")


def test_algorithms():
    print(f"\n{'=' * 70}")
    print("TEST 9: DIMLANG ALGORITHMS (composite operations)")
    print("=" * 70)

    m = Message()

    # Algorithm 1: GREEDY GRADIENT DESCENT toward collision
    print(f"\n  ALGORITHM: Greedy gradient descent")
    p_target = Path(m)
    m_current = Message()  # random start
    d_current = overlay(Path(m_current), p_target).weight
    print(f"    Initial distance: {d_current}")

    steps = 0
    history = [d_current]
    for iteration in range(50):
        best_d = d_current
        best_move = None

        # Try flipping each bit
        for word in range(16):
            for bit in range(0, 32, 4):  # sample every 4th bit for speed
                m_test = m_current.copy()
                m_test[word] = m_test[word].val ^ (1 << bit)
                d_test = overlay(Path(m_test), p_target).weight
                if d_test < best_d:
                    best_d = d_test
                    best_move = (word, bit)

        if best_move is None:
            print(f"    Stuck at iteration {iteration}, dist = {d_current}")
            break

        m_current[best_move[0]] = m_current[best_move[0]].val ^ (1 << best_move[1])
        d_current = best_d
        history.append(d_current)
        steps += 1

    print(f"    After {steps} steps: distance = {d_current}")
    print(f"    History: {history[:10]}...")
    print(f"    → {'CONVERGING!' if d_current < history[0] - 20 else 'Gradient too flat.'}")

    # Algorithm 2: REPAIR + PROPAGATE
    print(f"\n  ALGORITHM: Multi-strategy repair comparison")
    results = {}
    for strat in ['a', 'e', 'alternate', 'smart']:
        dHs = []
        for trial in range(100):
            m_trial = Message()
            dw = np.random.randint(1, 2**32)
            _, dH = dual_repair(m_trial, dw, strat)
            dHs.append(dH)
        results[strat] = (np.mean(dHs), min(dHs), max(dHs))
        print(f"    {strat:>9}: mean={results[strat][0]:.1f}, min={results[strat][1]}, max={results[strat][2]}")

    # Algorithm 3: FABRIC RESONANCE SEARCH
    print(f"\n  ALGORITHM: Fabric resonance at internal rounds")
    fab = Fabric()
    paths_list = []
    for i in range(50):
        msg = Message()
        fab.weave(msg)
        paths_list.append(fab.paths[-1])

    # Find: at which INTERNAL round are paths CLOSEST?
    min_per_round = []
    for r in range(0, 65, 8):
        min_d = 256
        for i in range(50):
            for j in range(i+1, min(i+10, 50)):  # limited pairs for speed
                d = overlay(paths_list[i], paths_list[j], r).weight
                if d < min_d: min_d = d
        min_per_round.append((r, min_d))

    print(f"    Min distance per round (50 paths, limited pairs):")
    for r, d in min_per_round:
        bar = "#" * (d // 4)
        print(f"      r={r:2d}: {d:3d} {bar}")


if __name__ == "__main__":
    test_basics()
    test_overlay()
    test_repair()
    test_invert()
    test_tunnel()
    test_warp()
    test_resonate()
    test_new_operations()
    test_algorithms()

    print(f"\n{'=' * 70}")
    print("ALL TESTS COMPLETE")
    print("=" * 70)
