"""
Stage 1, Part 4: Searching for a foothold — multiple angles.

What we know WON'T work:
- Differentials (random after 6 rounds)
- Tracking complexity across rounds (destroyed by rotation)
- Simple invariants (none found)
- Any single algebra (GF(2) or Z/2^32 each see only 50%)

Angles to try:
A. The INVERSE direction: what does F^{-1} look like?
B. Fixed points and cycles: does F have structure as a DYNAMICAL SYSTEM?
C. Bit-slice: treat each of 32 bit-planes independently
D. The SCHEDULE as constraint: W[16..63] = f(W[0..15]) — what does this do?
E. Two-message JOINT trajectory: not diff, but the PAIR together
"""

import random
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, schedule
)


def experiment_A_inverse():
    """
    SHA-256 round function is invertible (given the full state).
    Given state[r+1] and W[r], K[r], we can recover state[r].

    The inverse uses: a_prev = b, b_prev = c, c_prev = d, e_prev = f, etc.
    Then: T1 = e_new - d_prev, T2 = a_new - T1

    Question: Is F^{-1} "simpler" than F in any sense?
    Specifically: does it have fewer carries? Simpler bit dependencies?
    """
    print("=" * 80)
    print("ANGLE A: The inverse — is F^{-1} simpler than F?")
    print("=" * 80)

    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]
    states, W = sha256_round_trace(M)

    # Verify inverse
    violations = 0
    for r in range(63, 0, -1):
        a, b, c, d, e, f, g, h = states[r]

        # From state[r], recover state[r-1]:
        # b = a_{r-1}, c = b_{r-1} = a_{r-2}, d = c_{r-1} = a_{r-3}
        # f = e_{r-1}, g = f_{r-1} = e_{r-2}, h = g_{r-1} = e_{r-3}
        a_prev = b
        # e_prev = f
        # T1 = e - d_prev = e - c_{r-1} = e - d
        # Wait: e_r = d_{r-1} + T1_{r-1}, so T1_{r-1} = e_r - d_{r-1}
        # But d_{r-1} = c_{r-2} = ... we need state[r-1]'s d

        # Actually: from state[r]:
        # a_r was computed as T1_{r-1} + T2_{r-1}
        # e_r was computed as d_{r-1} + T1_{r-1}
        # So: T1_{r-1} = e_r - d_{r-1}
        #     T2_{r-1} = a_r - T1_{r-1}
        # d_{r-1} = state[r-1][3]

        # We have state[r-1] to verify:
        d_prev = states[r-1][3]
        T1_recovered = (e - d_prev) & MASK
        T2_recovered = (a - T1_recovered) & MASK

        # Verify T2 = Sig0(a_{r-1}) + Maj(a_{r-1}, b_{r-1}, c_{r-1})
        a_p, b_p, c_p = states[r-1][0], states[r-1][1], states[r-1][2]
        T2_expected = (Sig0(a_p) + Maj(a_p, b_p, c_p)) & MASK

        if T2_recovered != T2_expected:
            violations += 1

    print(f"  Inverse verification: {violations}/63 violations")
    if violations == 0:
        print(f"  F is exactly invertible ✓")

    # The inverse requires subtraction (a - T1) which has the SAME carry structure
    # as addition. So F^{-1} is NOT simpler than F.
    print(f"  F^{{-1}} uses subtraction = same carry structure as addition")
    print(f"  → F^{{-1}} is NOT simpler than F")


def experiment_B_dynamics():
    """
    View SHA-256 as a dynamical system: state → state.
    Question: Are there periodic orbits? Near-periodic behavior?

    If we run SHA-256 with FIXED W (same message word every round),
    does the state eventually cycle?

    This removes the schedule complication and isolates the round function.
    """
    print("\n" + "=" * 80)
    print("ANGLE B: Dynamical system — periodicity with fixed W")
    print("=" * 80)

    # Use fixed W = 0 and K = 0 (strip away constants)
    # Run the round function many times and look for cycles
    state = list(H0)
    seen = {}
    cycle_found = False

    for step in range(10000):
        key = tuple(state)
        if key in seen:
            cycle_len = step - seen[key]
            print(f"  CYCLE found! Length = {cycle_len}, starts at step {seen[key]}")
            cycle_found = True
            break
        seen[key] = step

        a, b, c, d, e, f, g, h = state
        T1 = (h + Sig1(e) + Ch(e, f, g)) & MASK  # K=0, W=0
        T2 = (Sig0(a) + Maj(a, b, c)) & MASK
        a_new = (T1 + T2) & MASK
        e_new = (d + T1) & MASK
        state = [a_new, a, b, c, e_new, e, f, g]

    if not cycle_found:
        print(f"  No cycle in 10000 steps")

        # Check: is state approaching a fixed point?
        # Measure distance between consecutive states
        dists = []
        state = list(H0)
        for step in range(100):
            a, b, c, d, e, f, g, h = state
            T1 = (h + Sig1(e) + Ch(e, f, g)) & MASK
            T2 = (Sig0(a) + Maj(a, b, c)) & MASK
            a_new = (T1 + T2) & MASK
            e_new = (d + T1) & MASK
            new_state = [a_new, a, b, c, e_new, e, f, g]

            dist = sum(bin(s1 ^ s2).count('1') for s1, s2 in zip(state, new_state))
            dists.append(dist)
            state = new_state

        print(f"  Distance between consecutive states:")
        print(f"    First 5:  {dists[:5]}")
        print(f"    Last 5:   {dists[-5:]}")
        print(f"    Average:  {sum(dists)/len(dists):.1f} (random=128)")


def experiment_C_bitslice():
    """
    Treat SHA-256 as 32 PARALLEL 1-bit systems.

    At bit position k, the round function is:
      a_new[k] = T1[k] XOR T2[k] XOR carry_into_k(T1+T2)
      T1[k] = h[k] XOR Sig1_k(e) XOR Ch_k(e,f,g) XOR K[r][k] XOR W[r][k]
              XOR carries from T1 additions
      T2[k] = Sig0_k(a) XOR Maj_k(a,b,c) XOR carry from Sig0+Maj addition

    For k=0: NO carries at all → the 1-bit system at position 0 is pure GF(2).
    For k>0: carries couple position k to positions 0..k-1.

    Question: If we extract just the k=0 bit-slice across all 64 rounds,
    is it a COMPLETE GF(2) system (closed, no coupling to other bits)?
    """
    print("\n" + "=" * 80)
    print("ANGLE C: Bit-slice at position 0 — is it a closed GF(2) system?")
    print("=" * 80)

    # Bit 0 of a[r+1] depends on:
    # h[0], e[6], e[11], e[25], e[0], f[0], g[0], a[2], a[13], a[22], a[0], b[0], c[0]
    # K[r][0], W[r][0]
    #
    # Through register shifts: b[0] = a_{r-1}[0], c[0] = a_{r-2}[0], etc.
    # But e[6] = bit 6 of e, which is NOT bit 0 of anything!
    # e[6] has carry contamination from the round that produced it.
    #
    # So the bit-0 slice is NOT closed — it depends on bits 2, 6, 11, 13, 22, 25
    # of OTHER positions.

    print(f"  Bit 0 of a[r+1] depends on bits at positions: 0, 2, 6, 11, 13, 22, 25")
    print(f"  of the current state registers.")
    print(f"  Positions 2, 6, 11, 13, 22, 25 are NOT bit-0 → slice is NOT closed.")
    print(f"")
    print(f"  These positions come from Sig0 rotations {{2,13,22}} and Sig1 rotations {{6,11,25}}.")
    print(f"  They couple bit 0 to bits that HAVE carry contamination.")
    print(f"")

    # But: what if we define a LARGER slice that IS closed?
    # Bit 0 needs bits {0, 2, 6, 11, 13, 22, 25}.
    # Each of those needs their own dependencies...
    # This grows exponentially → eventually covers all 32 bits.

    # Let's measure: starting from bit 0, how fast does the dependency set grow?
    needed = {0}
    round_deps = []
    for step in range(8):
        new_needed = set()
        for bit in needed:
            # Each bit position k depends on:
            # From Sig0: (k+2)%32, (k+13)%32, (k+22)%32
            # From Sig1: (k+6)%32, (k+11)%32, (k+25)%32
            # From Ch/Maj: k itself
            # Plus carry: bits 0..k-1 (but only within one addition)
            for offset in [0, 2, 6, 11, 13, 22, 25]:
                new_needed.add((bit + offset) % 32)
        needed = needed | new_needed
        round_deps.append(len(needed))
        print(f"  After {step+1} expansion(s): {len(needed)}/32 bit positions needed")
        if len(needed) == 32:
            print(f"  → ALL 32 positions needed after {step+1} expansions")
            break

    print(f"\n  Conclusion: No closed bit-slice exists smaller than all 32 bits.")
    print(f"  The rotation constants {{2,6,11,13,22,25}} generate all of Z/32Z.")


def experiment_D_schedule_constraint():
    """
    The message schedule W[16..63] = f(W[0..15]) is a HARD CONSTRAINT.
    It's GF(2)-linear (sig0, sig1 are linear over GF(2) for XOR).
    But it's nonlinear over Z/2^32 (because sig0/sig1 use addition).

    Wait — actually sig0 and sig1 in the schedule are:
      sig0(x) = ROTR(x,7) XOR ROTR(x,18) XOR (x >> 3)
      sig1(x) = ROTR(x,17) XOR ROTR(x,19) XOR (x >> 10)

    These are PURE XOR/shift — NO addition! So the schedule IS GF(2)-linear.

    W[i] = sig1(W[i-2]) XOR W[i-7] XOR sig0(W[i-15]) XOR W[i-16]
    (over GF(2))

    But wait: the actual schedule uses + not XOR:
    W[i] = (sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16]) mod 2^32

    So the schedule has ADDITIONS too. Let's check how much carry matters.
    """
    print("\n" + "=" * 80)
    print("ANGLE D: Schedule — XOR vs ADD deviation")
    print("=" * 80)

    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]
    W_real = schedule(M)

    # XOR-only schedule
    W_xor = list(M)
    for i in range(16, 64):
        s1 = rotr(W_xor[i-2], 17) ^ rotr(W_xor[i-2], 19) ^ (W_xor[i-2] >> 10)
        s0 = rotr(W_xor[i-15], 7) ^ rotr(W_xor[i-15], 18) ^ (W_xor[i-15] >> 3)
        W_xor.append(s1 ^ W_xor[i-7] ^ s0 ^ W_xor[i-16])

    print(f"  Schedule deviation (real ADD vs XOR-only):")
    for i in [16, 17, 20, 32, 48, 63]:
        diff = W_real[i] ^ W_xor[i]
        hw = bin(diff).count('1')
        print(f"    W[{i:>2}]: HW(diff) = {hw}")

    # Average
    avg_hw = sum(bin(W_real[i] ^ W_xor[i]).count('1') for i in range(16, 64)) / 48
    print(f"    Average HW(diff) for W[16..63]: {avg_hw:.1f}/32")
    print(f"    Schedule carry contribution ≈ {avg_hw:.0f}/32 = {avg_hw/32:.0%}")


def experiment_E_joint_trajectory():
    """
    Instead of comparing M and M' (differential), look at them TOGETHER.

    For two messages M, M', define the "joint state" at round r as:
      J[r] = (state(M)[r], state(M')[r])  — a 512-bit object

    The joint state evolves under a JOINT round function:
      J[r+1] = (F(state_1, K, W_1), F(state_2, K, W_2))

    The two halves are INDEPENDENT (they share K but have different W).
    So the joint state carries no extra information... unless we look at
    specific FUNCTIONS of the joint state.

    The differential δ[r] = state_1 XOR state_2 is one such function.
    We proved it goes random. Are there OTHERS?

    For example: state_1[r] AND state_2[r]? Or state_1[r] + state_2[r]?
    """
    print("\n" + "=" * 80)
    print("ANGLE E: Joint trajectory — functions of (state_1, state_2)")
    print("=" * 80)

    N = 300
    rng = random.Random(42)
    M1 = [rng.randint(0, MASK) for _ in range(16)]
    M2 = [rng.randint(0, MASK) for _ in range(16)]

    s1, _ = sha256_round_trace(M1)
    s2, _ = sha256_round_trace(M2)

    print(f"\nJoint functions of two INDEPENDENT messages:")
    print(f"{'r':>3} | {'HW(XOR)':>8} {'HW(AND)':>8} {'HW(OR)':>8} {'(a1+a2)%2^32':>14} {'HW(sum)':>8}")
    print("-" * 60)

    for r in [0, 1, 4, 8, 16, 32, 64]:
        xor_hw = sum(bin(s1[r][i] ^ s2[r][i]).count('1') for i in range(8))
        and_hw = sum(bin(s1[r][i] & s2[r][i]).count('1') for i in range(8))
        or_hw = sum(bin(s1[r][i] | s2[r][i]).count('1') for i in range(8))
        sum_a = (s1[r][0] + s2[r][0]) & MASK
        sum_hw = bin(sum_a).count('1')

        print(f"{r:>3} | {xor_hw:>8} {and_hw:>8} {or_hw:>8} {sum_a:>14x} {sum_hw:>8}")

    # For independent random 256-bit vectors:
    # E[HW(XOR)] = 128, E[HW(AND)] = 64, E[HW(OR)] = 192
    print(f"\n  Random expectation: XOR=128, AND=64, OR=192")
    print(f"  All joint functions → random by round 4")
    print(f"  Independent messages carry no joint structure.")

    # What about RELATED messages (M2 = M1 with one bit flipped)?
    print(f"\n  Related messages (M2 = M1 XOR bit 0 of W[0]):")
    M2_related = list(M1)
    M2_related[0] ^= 1
    s2r, _ = sha256_round_trace(M2_related)

    for r in [0, 1, 4, 8, 16, 32, 64]:
        xor_hw = sum(bin(s1[r][i] ^ s2r[r][i]).count('1') for i in range(8))
        and_hw = sum(bin(s1[r][i] & s2r[r][i]).count('1') for i in range(8))
        print(f"  r={r:>2}: HW(XOR)={xor_hw:>3}, HW(AND)={and_hw:>3}")


def experiment_F_round_function_as_permutation():
    """
    For FIXED W[r] and K[r], the round function F_r is a map
    from (Z/2^32)^8 → (Z/2^32)^8.

    Is it a permutation? (Invertible = yes, as we showed.)
    What is its ORDER? (Smallest k with F_r^k = identity)

    If the order is small → there's exploitable periodicity.
    If the order is astronomical → no structure.
    """
    print("\n" + "=" * 80)
    print("ANGLE F: Round function as permutation — what is its order?")
    print("=" * 80)

    # We can't compute the exact order (state space = 2^256).
    # But we can check: starting from a random state, how long until
    # F_r^k(state) returns close to state?

    # Use W=0, K=0 for simplicity
    state = list(H0)
    initial = tuple(state)

    # Track HW of (state XOR initial) — if it drops to 0, we cycled
    print(f"  Iterating F with W=0, K=0 from IV:")
    min_dist = 256
    min_dist_step = 0

    for step in range(1, 100001):
        a, b, c, d, e, f, g, h = state
        T1 = (h + Sig1(e) + Ch(e, f, g)) & MASK
        T2 = (Sig0(a) + Maj(a, b, c)) & MASK
        a_new = (T1 + T2) & MASK
        e_new = (d + T1) & MASK
        state = [a_new, a, b, c, e_new, e, f, g]

        dist = sum(bin(s ^ i).count('1') for s, i in zip(state, initial))
        if dist < min_dist:
            min_dist = dist
            min_dist_step = step

        if dist == 0:
            print(f"  CYCLE at step {step}!")
            break

    else:
        print(f"  No cycle in 100000 steps")
        print(f"  Closest approach: HW={min_dist} at step {min_dist_step}")
        print(f"  (random distance = 128)")


if __name__ == "__main__":
    experiment_A_inverse()
    experiment_B_dynamics()
    experiment_C_bitslice()
    experiment_D_schedule_constraint()
    experiment_E_joint_trajectory()
    experiment_F_round_function_as_permutation()
