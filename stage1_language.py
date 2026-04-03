"""
Stage 1: WHY do existing languages fail for F?

F = one round of SHA-256 a-recurrence.
F mixes: rotations (bit permutation), XOR (GF(2) addition),
         AND (GF(2) multiplication), mod 2^32 addition (carry chain).

Each of these operations alone is SIMPLE:
  - Rotations: trivial permutation
  - XOR: linear over GF(2)
  - AND: quadratic over GF(2)
  - mod 2^32 add: linear over Z/2^32, but nonlinear over GF(2)

The problem: F uses ALL of them TOGETHER.
- GF(2) sees XOR perfectly but mod-add creates degree explosion
- Z/2^32 sees mod-add perfectly but XOR creates "bit mixing"
- Neither sees the WHOLE picture

Question: What does ONE round of F do to the information
in a way that NEITHER GF(2) NOR Z/2^32 can describe?
"""

import random
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, schedule
)


def experiment_17_operation_interference():
    """
    The core issue: XOR and ADD interfere through carries.

    (a XOR b) != (a + b) mod 2^32 because of carries.
    The DIFFERENCE is exactly the carry contribution:
      a + b = (a XOR b) + 2*(a AND b)  ... but this is recursive (carry propagation)

    Formally: a + b = a XOR b XOR carry_chain(a, b)
    where carry_chain itself depends on a AND b at each bit.

    This means: the carry is the COUPLING between GF(2) and Z/2^32.
    The carry is what makes neither language sufficient alone.

    Let's measure: in F, how much of the output comes from
    "pure XOR" vs "carry corrections"?
    """
    print("=" * 80)
    print("EXPERIMENT 17: XOR vs carry decomposition of F")
    print("=" * 80)

    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]
    states, W = sha256_round_trace(M)

    print("\nFor each round: compute a[r+1] via mod-add vs XOR-only")
    print("The DIFFERENCE = total carry contribution\n")

    print(f"{'r':>3} | {'a[r+1] real':>10} | {'a[r+1] XOR-only':>15} | {'carry contribution':>18} | {'HW(carry)':>9}")
    print("-" * 75)

    for r in range(20):
        a, b, c, d, e, f, g, h = states[r]

        # Real computation (mod 2^32 add)
        T1_real = (h + Sig1(e) + Ch(e, f, g) + K[r] + W[r]) & MASK
        T2_real = (Sig0(a) + Maj(a, b, c)) & MASK
        a_real = (T1_real + T2_real) & MASK

        # XOR-only computation (replace all + with XOR)
        T1_xor = h ^ Sig1(e) ^ Ch(e, f, g) ^ K[r] ^ W[r]
        T2_xor = Sig0(a) ^ Maj(a, b, c)
        a_xor = T1_xor ^ T2_xor

        carry_contrib = a_real ^ a_xor
        hw_carry = bin(carry_contrib).count('1')

        print(f"{r:>3} | {a_real:>10x} | {a_xor:>15x} | {carry_contrib:>18x} | {hw_carry:>9}")

    # Statistics over many messages
    print(f"\n--- Statistics over 500 messages, rounds 8-63 ---")
    hw_carries = []
    for seed in range(500):
        rng_s = random.Random(seed)
        M = [rng_s.randint(0, MASK) for _ in range(16)]
        states, W = sha256_round_trace(M)
        for r in range(8, 64):
            a, b, c, d, e, f, g, h = states[r]
            T1_real = (h + Sig1(e) + Ch(e, f, g) + K[r] + W[r]) & MASK
            T2_real = (Sig0(a) + Maj(a, b, c)) & MASK
            a_real = (T1_real + T2_real) & MASK
            T1_xor = h ^ Sig1(e) ^ Ch(e, f, g) ^ K[r] ^ W[r]
            T2_xor = Sig0(a) ^ Maj(a, b, c)
            a_xor = T1_xor ^ T2_xor
            hw_carries.append(bin(a_real ^ a_xor).count('1'))

    avg_hw = sum(hw_carries) / len(hw_carries)
    print(f"E[HW(carry contribution)] = {avg_hw:.2f} out of 32 bits")
    print(f"Fraction of output determined by carries: {avg_hw/32:.1%}")
    print(f"Carries are NOT small corrections — they affect ~{avg_hw:.0f}/32 bits!")


def experiment_18_what_F_preserves():
    """
    If F doesn't preserve any simple quantity (exp 11 showed this),
    maybe it preserves something MORE COMPLEX.

    Idea: instead of looking for g(state) = const,
    look for g(state[r], state[r+1]) = const.
    I.e., a relation between CONSECUTIVE states that always holds.

    We already know: a[r+1] = F(state[r], K[r], W[r]).
    That IS a relation — but it's the trivial one (the recurrence itself).

    Are there OTHERS? I.e., non-trivial functions of (state[r], state[r+1])
    that are independent of M?

    If yes → these are the "laws of motion" beyond the obvious recurrence.
    """
    print("\n" + "=" * 80)
    print("EXPERIMENT 18: Relations between consecutive states")
    print("=" * 80)

    N = 300

    all_data = []
    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        states, W = sha256_round_trace(M)
        all_data.append((states, W))

    # We know a[r+1] = F(state[r], K[r], W[r]). This uses K[r] and W[r].
    # Are there relations NOT involving W[r]?

    # Key insight: a[r+1] and e[r+1] together encode ALL of F.
    # But: a[r+1] - e[r+1] = T2[r] - d[r] (from experiment 10).
    # T2[r] = Sig0(a[r]) + Maj(a[r], b[r], c[r])
    # d[r] = a[r-3]
    # So: a[r+1] - e[r+1] = Sig0(a[r]) + Maj(a[r], a[r-1], a[r-2]) - a[r-3]
    #
    # This relation does NOT involve K[r] or W[r]!
    # It's a PURE structural relation of the a/e sequences.
    #
    # But we verified this in experiment 10. Let's go DEEPER.

    # Question: Does (a[r+1] - e[r+1]) have any pattern across rounds?
    # If it's random → no additional structure.
    # If it's NOT random → there's structure in the a-trajectory itself.

    print(f"\nPattern of T2[r] - d[r] = a[r+1] - e[r+1] across rounds:")
    print(f"(This quantity depends only on a[r], a[r-1], a[r-2], a[r-3])\n")

    # For one message: show the sequence
    states, W = all_data[0]
    vals = []
    for r in range(4, 64):
        v = (states[r+1][0] - states[r+1][4]) & MASK
        vals.append(v)

    # Is this sequence more structured than random?
    # Test: autocorrelation
    mean_hw = sum(bin(v).count('1') for v in vals) / len(vals)
    print(f"E[HW(a[r+1]-e[r+1])] = {mean_hw:.2f} (random=16)")

    # Unique values
    print(f"Unique values: {len(set(vals))}/{len(vals)}")

    # XOR consecutive
    xor_consecutive = [vals[i] ^ vals[i+1] for i in range(len(vals)-1)]
    mean_xor = sum(bin(v).count('1') for v in xor_consecutive) / len(xor_consecutive)
    print(f"E[HW(v[r] XOR v[r+1])] = {mean_xor:.2f} (random=16)")


def experiment_19_the_gap():
    """
    The GAP between GF(2) and Z/2^32 is the carry.

    For a + b (mod 2^32):
      result = a XOR b XOR carry_vector
      carry_vector[0] = 0
      carry_vector[i] = MAJ(a[i-1], b[i-1], carry_vector[i-1])

    The carry propagation is a SEQUENTIAL process — bit i depends on bit i-1.
    This is fundamentally different from both XOR (parallel) and the
    nonlinear functions Ch/Maj (parallel per bit).

    The carry is the ONLY sequential element in F.
    Everything else (Sig0, Sig1, Ch, Maj, XOR) is parallel per bit.

    What if the right language describes carries as a FIRST-CLASS object,
    not as a side effect of addition?
    """
    print("\n" + "=" * 80)
    print("EXPERIMENT 19: Carry as first-class object")
    print("=" * 80)

    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]
    states, W = sha256_round_trace(M)

    # For round r, decompose each addition into (XOR part, carry part)
    r = 20
    a, b, c, d, e, f, g, h = states[r]

    def add_decompose(x, y):
        """Decompose x + y into XOR part and carry vector."""
        result = (x + y) & MASK
        xor_part = x ^ y
        carry_part = result ^ xor_part  # carry contribution

        # Also compute carry VECTOR (which positions had carries)
        carry_vec = 0
        c_bit = 0
        for i in range(32):
            xi = (x >> i) & 1
            yi = (y >> i) & 1
            # carry into position i is c_bit
            if c_bit:
                carry_vec |= (1 << i)
            # carry out of position i
            c_bit = (xi & yi) | (xi & c_bit) | (yi & c_bit)

        return result, xor_part, carry_part, carry_vec

    print(f"\nRound {r} addition decomposition:")
    print(f"{'operation':>20} | {'result':>10} | {'XOR part':>10} | {'carry part':>10} | {'HW(carry)':>9}")
    print("-" * 75)

    # T1 = h + Sig1(e) + Ch(e,f,g) + K[r] + W[r]
    s1e = Sig1(e)
    che = Ch(e, f, g)
    kw = (K[r] + W[r]) & MASK

    additions = [
        ("h + Sig1(e)", h, s1e),
        ("prev + Ch", (h + s1e) & MASK, che),
        ("prev + K+W", (h + s1e + che) & MASK, kw),
    ]

    total_carry_hw = 0
    for name, x, y in additions:
        result, xor_p, carry_p, carry_v = add_decompose(x, y)
        hw = bin(carry_v).count('1')
        total_carry_hw += hw
        print(f"{name:>20} | {result:>10x} | {xor_p:>10x} | {carry_p:>10x} | {hw:>9}")

    # T2
    s0a = Sig0(a)
    maja = Maj(a, b, c)
    result_t2, xor_t2, carry_t2, carry_v_t2 = add_decompose(s0a, maja)
    hw_t2 = bin(carry_v_t2).count('1')
    total_carry_hw += hw_t2
    print(f"{'Sig0(a) + Maj':>20} | {result_t2:>10x} | {xor_t2:>10x} | {carry_t2:>10x} | {hw_t2:>9}")

    # T1 + T2
    T1 = (h + s1e + che + kw) & MASK
    result_final, xor_final, carry_final, carry_v_final = add_decompose(T1, result_t2)
    hw_final = bin(carry_v_final).count('1')
    total_carry_hw += hw_final
    print(f"{'T1 + T2':>20} | {result_final:>10x} | {xor_final:>10x} | {carry_final:>10x} | {hw_final:>9}")

    print(f"\nTotal carry positions in one round: {total_carry_hw}")
    print(f"These carries are the COUPLING between XOR-world and ADD-world")

    # KEY OBSERVATION: Each carry bit is determined by a CHAIN of previous bits.
    # Carry at position i depends on bits 0..i-1 of BOTH operands.
    # This creates a DIRECTED DEPENDENCY: low bits → high bits.
    # There is NO reverse dependency: high bits never affect low bits.

    print(f"\n--- CRITICAL STRUCTURAL OBSERVATION ---")
    print(f"In mod 2^32 addition:")
    print(f"  - Bit 0 of result: NEVER has carry (always XOR)")
    print(f"  - Bit i of result: carry from bits 0..i-1")
    print(f"  - Bit 31 of result: carry from ALL 31 lower bits")
    print(f"")
    print(f"This means: there is a NATURAL ORDERING of bits.")
    print(f"Low bits are 'causes', high bits are 'effects'.")
    print(f"XOR ignores this ordering. Addition CREATES it.")
    print(f"")
    print(f"Sig0/Sig1 (rotations+XOR) SCRAMBLE this ordering.")
    print(f"Then the next addition RECREATES a new ordering.")
    print(f"")
    print(f"SHA-256 = repeated cycle of:")
    print(f"  1. Create ordering (addition, carry propagation)")
    print(f"  2. Destroy ordering (rotation + XOR)")
    print(f"  3. Create new ordering (next addition)")
    print(f"")
    print(f"No fixed language can describe this because the 'coordinate system'")
    print(f"CHANGES at every operation.")


def experiment_20_dual_representation():
    """
    Every 32-bit word x can be viewed simultaneously as:
      - A vector in GF(2)^32 (for XOR, AND, rotations)
      - An element of Z/2^32 (for addition)

    These are the SAME bits, but with DIFFERENT algebraic structure.

    The carry is what TRANSLATES between them:
      x +_{Z} y = x XOR y XOR carry(x,y)

    What if our new "language" keeps BOTH representations simultaneously?
    Not choosing GF(2) or Z/2^32, but tracking the PAIR.

    For any word x, define its "dual" as (x_gf2, x_int) where both
    are the same bits but we track which OPERATIONS were applied.

    Then: XOR operates on the gf2 component, ADD on the int component,
    and carry is the COUPLING between them.
    """
    print("\n" + "=" * 80)
    print("EXPERIMENT 20: Dual tracking — GF(2) shadow vs Z-shadow")
    print("=" * 80)

    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]
    states, W = sha256_round_trace(M)

    # For each round, compute a[r+1] TWO ways:
    # 1. "GF(2) shadow": replace all mod-add with XOR
    # 2. "Z shadow": the real computation
    # Track how far apart they drift

    a_gf2 = H0[:]  # GF(2) shadow state
    a_real = H0[:]  # Real state

    # For GF(2) shadow: use XOR for everything
    W_full = schedule(M)

    gf2_state = list(H0)  # [a,b,c,d,e,f,g,h]
    real_state = list(H0)

    print(f"{'r':>3} | {'HW(a_real^a_gf2)':>16} | {'HW(e_real^e_gf2)':>16} | {'comment':>30}")
    print("-" * 75)

    for r in range(64):
        a_r, b_r, c_r, d_r, e_r, f_r, g_r, h_r = real_state
        ag, bg, cg, dg, eg, fg, gg, hg = gf2_state

        # Real round
        T1_real = (h_r + Sig1(e_r) + Ch(e_r, f_r, g_r) + K[r] + W_full[r]) & MASK
        T2_real = (Sig0(a_r) + Maj(a_r, b_r, c_r)) & MASK
        a_new_real = (T1_real + T2_real) & MASK
        e_new_real = (d_r + T1_real) & MASK

        # GF(2) round (XOR instead of ADD)
        T1_gf2 = hg ^ Sig1(eg) ^ Ch(eg, fg, gg) ^ K[r] ^ W_full[r]
        T2_gf2 = Sig0(ag) ^ Maj(ag, bg, cg)
        a_new_gf2 = T1_gf2 ^ T2_gf2
        e_new_gf2 = dg ^ T1_gf2

        # Update states
        real_state = [a_new_real, a_r, b_r, c_r, e_new_real, e_r, f_r, g_r]
        gf2_state = [a_new_gf2, ag, bg, cg, e_new_gf2, eg, fg, gg]

        drift_a = bin(a_new_real ^ a_new_gf2).count('1')
        drift_e = bin(e_new_real ^ e_new_gf2).count('1')

        comment = ""
        if r == 0:
            comment = "← first divergence"
        elif r == 5:
            comment = "← approaching saturation"
        elif r == 63:
            comment = "← final round"

        if r < 10 or r >= 60 or r % 10 == 0:
            print(f"{r:>3} | {drift_a:>16} | {drift_e:>16} | {comment:>30}")

    print(f"\nThe GF(2) shadow diverges from reality within ~5 rounds.")
    print(f"After that: drift ≈ 16 bits = random.")
    print(f"This is EXACTLY the carry contribution — it accumulates.")
    print(f"")
    print(f"KEY INSIGHT: The 'language gap' between GF(2) and Z/2^32")
    print(f"is not a small perturbation. It's a FULL 50% of the bits.")
    print(f"Neither language captures even half of what F does.")
    print(f"")
    print(f"A new language must treat XOR and ADD as EQUALLY FUNDAMENTAL,")
    print(f"not one as primary and the other as correction.")


if __name__ == "__main__":
    experiment_17_operation_interference()
    experiment_18_what_F_preserves()
    experiment_19_the_gap()
    experiment_20_dual_representation()
