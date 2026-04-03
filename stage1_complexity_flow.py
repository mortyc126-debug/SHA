"""
Stage 1, Part 3: Complexity Flow — tracking algebraic depth through rounds.

We proved: within ONE round, bit k of a[r+1] has algebraic depth proportional to k
(bit 0 = degree 2 exact, bit 31 = degree ~33).

Rotations in Sig0/Sig1 move bits between positions:
  Sig0: positions {2, 13, 22} → bit 0
  Sig1: positions {6, 11, 25} → bit 0

So after one round, bit 0 of a[r+1] depends on bits at positions
{0, 2, 6, 11, 13, 22, 25} of the PREVIOUS state.

Each of THOSE bits has its own complexity from the round before that.

Question: Can we define a "complexity tag" for each bit that tracks
how many carry-levels deep it is? And does this tag have a pattern
that survives across rounds?

New object: for each bit position, assign a "carry depth" =
number of carry chain levels that contributed to its current value.
- After an addition: bit k has carry depth = k
- After XOR: carry depth = max of inputs (XOR doesn't add carries)
- After rotation by j: carry depth of position i = carry depth of position (i+j)%32
"""

import random
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, schedule
)


def track_carry_depth_one_round(depths_a, depths_e, depths_prev_a):
    """
    Given carry-depth arrays for registers a, e (and previous a's for b,c,d),
    compute carry-depth of new a and new e after one SHA-256 round.

    depths_a[i] = carry depth of bit i of register a at round r
    Similarly for e.

    Returns (new_depths_a, new_depths_e) after one round.
    """
    # T2 = Sig0(a) + Maj(a,b,c)
    # Sig0(a) = ROTR(a,2) XOR ROTR(a,13) XOR ROTR(a,22)
    # Bit i of Sig0(a) = a[(i+2)%32] XOR a[(i+13)%32] XOR a[(i+22)%32]
    # Carry depth = max of the three rotated depths (XOR doesn't add carries)

    sig0_depth = [0] * 32
    for i in range(32):
        sig0_depth[i] = max(depths_a[(i+2)%32], depths_a[(i+13)%32], depths_a[(i+22)%32])

    # Maj(a,b,c) = (a AND b) XOR (a AND c) XOR (b AND c)
    # Bit i: degree 2, depth = max of inputs at bit i
    # b = a_{r-1}, c = a_{r-2}
    depths_b = depths_prev_a[0]  # a_{r-1}
    depths_c = depths_prev_a[1]  # a_{r-2}

    maj_depth = [0] * 32
    for i in range(32):
        maj_depth[i] = max(depths_a[i], depths_b[i], depths_c[i])

    # T2 = Sig0(a) + Maj(a,b,c): ADDITION adds carry chain
    # Bit i of T2: depth = max(sig0_depth[i], maj_depth[i]) + i (carry from bits 0..i-1)
    # But carry depth from addition is more subtle:
    # carry into bit i comes from bits 0..i-1, each with their own depth
    # The carry depth at position i = max over j<i of (depth[j]) + 1
    # Simplification: just use i as the additional depth from carry chain

    t2_input_depth = [max(sig0_depth[i], maj_depth[i]) for i in range(32)]
    t2_depth = [0] * 32
    max_below = 0
    for i in range(32):
        carry_depth = max_below + 1 if i > 0 else 0
        t2_depth[i] = max(t2_input_depth[i], carry_depth)
        max_below = max(max_below, t2_input_depth[i])

    # T1 = h + Sig1(e) + Ch(e,f,g) + K + W
    # Sig1(e) = ROTR(e,6) XOR ROTR(e,11) XOR ROTR(e,25)
    sig1_depth = [0] * 32
    for i in range(32):
        sig1_depth[i] = max(depths_e[(i+6)%32], depths_e[(i+11)%32], depths_e[(i+25)%32])

    # Ch(e,f,g) = (e AND f) XOR (NOT e AND g)
    # f = e_{r-1}, g = e_{r-2}
    depths_f = depths_prev_a[2]  # We'll approximate: e depths from 1 round ago
    depths_g = depths_prev_a[3]  # e depths from 2 rounds ago

    ch_depth = [0] * 32
    for i in range(32):
        ch_depth[i] = max(depths_e[i], depths_f[i], depths_g[i])

    # h = e_{r-3}
    depths_h = depths_prev_a[4]  # e depths from 3 rounds ago

    # T1 = h + Sig1(e) + Ch(e,f,g) + K + W  (4 additions, K and W have depth 0)
    t1_input_depth = [0] * 32
    for i in range(32):
        t1_input_depth[i] = max(depths_h[i], sig1_depth[i], ch_depth[i])

    t1_depth = [0] * 32
    max_below = 0
    for i in range(32):
        carry_depth = max_below + 1 if i > 0 else 0
        t1_depth[i] = max(t1_input_depth[i], carry_depth)
        max_below = max(max_below, t1_input_depth[i])

    # a_new = T1 + T2 (one more addition)
    a_input_depth = [max(t1_depth[i], t2_depth[i]) for i in range(32)]
    new_depths_a = [0] * 32
    max_below = 0
    for i in range(32):
        carry_depth = max_below + 1 if i > 0 else 0
        new_depths_a[i] = max(a_input_depth[i], carry_depth)
        max_below = max(max_below, a_input_depth[i])

    # e_new = d + T1 (d = a_{r-3})
    depths_d = depths_prev_a[5]  # a_{r-3}
    e_input_depth = [max(depths_d[i], t1_depth[i]) for i in range(32)]
    new_depths_e = [0] * 32
    max_below = 0
    for i in range(32):
        carry_depth = max_below + 1 if i > 0 else 0
        new_depths_e[i] = max(e_input_depth[i], carry_depth)
        max_below = max(max_below, e_input_depth[i])

    return new_depths_a, new_depths_e


def experiment_25_complexity_flow():
    """
    Track carry-depth through 64 rounds.
    Start: all bits at depth 0 (message/IV bits are "fresh").
    Each round: update depths according to the rules above.
    """
    print("=" * 80)
    print("EXPERIMENT 25: Carry-depth flow through 64 rounds")
    print("=" * 80)

    # Initialize: all state bits at depth 0
    depth_a = [0] * 32
    depth_e = [0] * 32

    # History for b,c,d,f,g,h (shifted copies)
    # We need a_depths from rounds r-1, r-2, r-3
    # and e_depths from rounds r-1, r-2, r-3
    history_a = [[0]*32 for _ in range(6)]  # [a_{r-1}, a_{r-2}, a_{r-3} depths, e_{r-1}, e_{r-2}, e_{r-3}]
    history_e = [[0]*32 for _ in range(3)]

    print(f"\n{'r':>3} | {'depth_a[0]':>10} {'depth_a[15]':>11} {'depth_a[31]':>11} | {'depth_e[0]':>10} {'depth_e[31]':>11} | {'max(a)':>6} {'max(e)':>6} | {'mean(a)':>7} {'mean(e)':>7}")
    print("-" * 110)

    for r in range(64):
        # Pack history: [depths_b, depths_c, depths_d_from_a, depths_f, depths_g, depths_h_from_e]
        # b = a_{r-1}, c = a_{r-2}, d = a_{r-3}
        # f = e_{r-1}, g = e_{r-2}, h = e_{r-3}
        prev = [
            history_a[0] if len(history_a) > 0 else [0]*32,  # b = a_{r-1}
            history_a[1] if len(history_a) > 1 else [0]*32,  # c = a_{r-2}
            history_e[0] if len(history_e) > 0 else [0]*32,  # f = e_{r-1}
            history_e[1] if len(history_e) > 1 else [0]*32,  # g = e_{r-2}
            history_e[2] if len(history_e) > 2 else [0]*32,  # h = e_{r-3}
            history_a[2] if len(history_a) > 2 else [0]*32,  # d = a_{r-3}
        ]

        new_a, new_e = track_carry_depth_one_round(depth_a, depth_e, prev)

        max_a = max(new_a)
        max_e = max(new_e)
        mean_a = sum(new_a) / 32
        mean_e = sum(new_e) / 32

        print(f"{r:>3} | {new_a[0]:>10} {new_a[15]:>11} {new_a[31]:>11} | {new_e[0]:>10} {new_e[31]:>11} | {max_a:>6} {max_e:>6} | {mean_a:>7.1f} {mean_e:>7.1f}")

        # Update history
        history_a = [list(depth_a)] + history_a[:5]
        history_e = [list(depth_e)] + history_e[:2]

        depth_a = new_a
        depth_e = new_e

    print(f"\nFinal carry-depth spectrum of a[64]:")
    print(f"  Bit: {' '.join(f'{i:>4}' for i in range(32))}")
    print(f"  Dep: {' '.join(f'{depth_a[i]:>4}' for i in range(32))}")

    print(f"\nDoes depth SATURATE or keep growing?")
    print(f"  If saturates → there's a finite 'complexity ceiling'")
    print(f"  If grows linearly → complexity is unbounded (= random)")


def experiment_26_actual_bit_degree():
    """
    Instead of TRACKING carry depth theoretically,
    MEASURE actual algebraic degree of bit 0 at each round.

    Method: for bit 0 of a[r], test degree-1 (linearity) and degree-2:
    - Linear: f(x⊕y⊕z) = f(x)⊕f(y)⊕f(z)⊕f(0) for all x,y,z
    - Degree 2: sum over all x,y of f(x⊕y)⊕f(x)⊕f(y)⊕f(0) = 0

    Actually simpler: just count how often flipping 2 bits is predicted
    by XOR of single-bit effects (this tests degree > 1).
    """
    print("\n" + "=" * 80)
    print("EXPERIMENT 26: Actual nonlinearity of bit 0 across rounds")
    print("=" * 80)

    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]

    for target_r in [1, 2, 4, 8, 16, 32, 64]:
        states_base, _ = sha256_round_trace(M, rounds=target_r)
        base_bit0 = states_base[target_r][0] & 1

        # Single-bit effects
        single_effects = {}
        for word in range(min(4, 16)):  # Just first 4 words for speed
            for bit in range(32):
                M_flip = list(M)
                M_flip[word] ^= (1 << bit)
                states_f, _ = sha256_round_trace(M_flip, rounds=target_r)
                effect = (states_f[target_r][0] & 1) ^ base_bit0
                single_effects[(word, bit)] = effect

        # Two-bit effects: flip bits i and j, check if effect = effect_i XOR effect_j
        # If yes → degree 1 (linear). If no → degree ≥ 2.
        N_pairs = 200
        nonlinear_count = 0
        total_pairs = 0

        pairs_tested = set()
        rng2 = random.Random(target_r)
        for _ in range(N_pairs):
            w1, b1 = rng2.randint(0, min(3, 15)), rng2.randint(0, 31)
            w2, b2 = rng2.randint(0, min(3, 15)), rng2.randint(0, 31)
            if (w1, b1) == (w2, b2):
                continue
            key = (min((w1,b1),(w2,b2)), max((w1,b1),(w2,b2)))
            if key in pairs_tested:
                continue
            pairs_tested.add(key)

            M_flip2 = list(M)
            M_flip2[w1] ^= (1 << b1)
            M_flip2[w2] ^= (1 << b2)
            states_f2, _ = sha256_round_trace(M_flip2, rounds=target_r)
            actual_effect = (states_f2[target_r][0] & 1) ^ base_bit0

            predicted_linear = single_effects[(w1, b1)] ^ single_effects[(w2, b2)]

            if actual_effect != predicted_linear:
                nonlinear_count += 1
            total_pairs += 1

        nl_rate = nonlinear_count / total_pairs if total_pairs > 0 else 0
        print(f"  Round {target_r:>2}: nonlinearity of bit 0 = {nonlinear_count}/{total_pairs} = {nl_rate:.3f} "
              f"({'LINEAR' if nl_rate == 0 else 'DEGREE≥2' if nl_rate < 0.4 else 'HIGH DEGREE'})")


if __name__ == "__main__":
    experiment_25_complexity_flow()
    experiment_26_actual_bit_degree()
