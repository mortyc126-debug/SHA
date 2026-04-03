"""
Stage 1, Part 2: The Bit-Position Asymmetry.

Key hypothesis: Different bit positions in SHA-256 have fundamentally
different algebraic complexity. Bit 0 = simple (pure XOR in additions).
Bit 31 = complex (depends on 31-bit carry chain).

Rotations MIX these complexities. This create-destroy-recreate cycle
is the CORE of what SHA-256 does.

If true, this gives us a NATURAL "coordinate" for a new language:
not the bit position itself, but the COMPLEXITY LEVEL of each bit.
"""

import random
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, schedule
)


def experiment_21_bit0_linearity():
    """
    Bit 0 of (x + y) mod 2^32 = x[0] XOR y[0] (no carry into bit 0).
    This means: bit 0 of a[r+1] is a GF(2)-linear function of bit 0 of inputs.

    Let's verify: does bit 0 of a[r+1] satisfy superposition?
    I.e., for messages M1, M2, M3 = M1 XOR M2:
      a[r+1](M3)[bit0] = a[r+1](M1)[bit0] XOR a[r+1](M2)[bit0] ?

    If yes for all rounds → bit 0 carries a LINEAR signal through ALL 64 rounds.
    """
    print("=" * 80)
    print("EXPERIMENT 21: Is bit 0 of a[r+1] linear across messages?")
    print("=" * 80)

    N_tests = 200
    max_violations_by_round = [0] * 65

    for trial in range(N_tests):
        rng = random.Random(trial * 3)
        M1 = [rng.randint(0, MASK) for _ in range(16)]
        M2 = [rng.randint(0, MASK) for _ in range(16)]
        M3 = [m1 ^ m2 for m1, m2 in zip(M1, M2)]  # M3 = M1 XOR M2

        s1, _ = sha256_round_trace(M1)
        s2, _ = sha256_round_trace(M2)
        s3, _ = sha256_round_trace(M3)

        for r in range(65):
            a1_bit0 = s1[r][0] & 1
            a2_bit0 = s2[r][0] & 1
            a3_bit0 = s3[r][0] & 1

            # Linear over GF(2) would mean: a3[bit0] = a1[bit0] XOR a2[bit0]
            # But we need to account for IV: a(M)[bit0] = L(M)[bit0] XOR IV_contribution
            # For XOR linearity: a(M1 XOR M2) = a(M1) XOR a(M2) XOR a(0)
            # Let's just check if a3 = a1 XOR a2 (simplest test)
            if a3_bit0 != (a1_bit0 ^ a2_bit0):
                max_violations_by_round[r] += 1

    print(f"\nTest: a(M1 XOR M2)[bit0] = a(M1)[bit0] XOR a(M2)[bit0]")
    print(f"(Violations out of {N_tests} trials per round)\n")

    for r in range(min(20, 65)):
        v = max_violations_by_round[r]
        pct = v / N_tests * 100
        status = "LINEAR ✓" if v == 0 else f"NONLINEAR ({pct:.0f}%)"
        print(f"  r={r:>2}: {v:>3}/{N_tests} violations → {status}")

    # Even if bit 0 is not XOR-linear (because of IV),
    # it might still be simpler than other bits.
    # Let's compare: how many input bits does bit k of a[r] depend on?


def experiment_22_bit_sensitivity():
    """
    For each bit position k of a[r+1], how sensitive is it to input bits?

    Measure: flip each of the 512 message bits, count how often
    bit k of a[round] changes. This gives sensitivity per bit position.

    If bit 0 is truly simpler: it should be sensitive to FEWER input bits.
    """
    print("\n" + "=" * 80)
    print("EXPERIMENT 22: Sensitivity of each bit position of a[r] to input bits")
    print("=" * 80)

    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]
    states_base, _ = sha256_round_trace(M)

    # For round 32 (middle of computation) and round 64 (output):
    for target_round in [1, 4, 8, 16, 32, 64]:
        a_base = states_base[target_round][0]

        sensitivity = [0] * 32  # how many input bits affect each output bit

        for word in range(16):
            for bit in range(32):
                M_flip = list(M)
                M_flip[word] ^= (1 << bit)
                states_flip, _ = sha256_round_trace(M_flip, rounds=target_round)
                a_flip = states_flip[target_round][0]

                diff = a_base ^ a_flip
                for k in range(32):
                    if (diff >> k) & 1:
                        sensitivity[k] += 1

        # Show sensitivity by bit position
        print(f"\n  Round {target_round}: sensitivity of each bit of a[r] to 512 input bits")
        print(f"  (max possible = 512)")

        # Group by 8
        for group_start in [0, 8, 16, 24]:
            bits_str = " ".join(f"{sensitivity[k]:>3}" for k in range(group_start, group_start + 8))
            print(f"    bits {group_start:>2}-{group_start+7:>2}: {bits_str}")

        # Is there a gradient? (low bits less sensitive than high bits?)
        avg_low = sum(sensitivity[0:8]) / 8
        avg_mid = sum(sensitivity[8:24]) / 16
        avg_high = sum(sensitivity[24:32]) / 8
        print(f"    Average: low[0-7]={avg_low:.1f}  mid[8-23]={avg_mid:.1f}  high[24-31]={avg_high:.1f}")


def experiment_23_carry_free_bit():
    """
    Bit 0 of ANY mod-2^32 addition has no carry-in.
    But after rotation by k, what was bit 0 is now at position k.

    In Sig1 = ROTR(e, 6) XOR ROTR(e, 11) XOR ROTR(e, 25):
      Bit 0 of Sig1 = e[6] XOR e[11] XOR e[25]
      (pure XOR of three specific bits of e)

    In Ch(e,f,g): bit 0 = (e[0] AND f[0]) XOR (NOT e[0] AND g[0])
      (depends on bit 0 of e, f, g only)

    In the addition h + Sig1(e):
      Bit 0 = h[0] XOR Sig1(e)[0] = h[0] XOR e[6] XOR e[11] XOR e[25]
      (pure XOR of 4 bits!)

    Let's trace bit 0 through the ENTIRE T1 computation.
    """
    print("\n" + "=" * 80)
    print("EXPERIMENT 23: Tracing bit 0 through one round (no carries)")
    print("=" * 80)

    # T1 = h + Sig1(e) + Ch(e,f,g) + K[r] + W[r]
    # Bit 0 of T1 (no carry into bit 0 at ANY addition):
    # Step 1: h[0] XOR Sig1(e)[0]
    #        = h[0] XOR e[6] XOR e[11] XOR e[25]
    # Step 2: XOR Ch(e,f,g)[0]
    #        = h[0] XOR e[6] XOR e[11] XOR e[25] XOR (e[0]&f[0] XOR ~e[0]&g[0])
    # Step 3: XOR K[r][0]
    # Step 4: XOR W[r][0]
    #
    # T2 = Sig0(a) + Maj(a,b,c)
    # Bit 0 of T2 (no carry):
    #        = Sig0(a)[0] XOR Maj(a,b,c)[0]
    #        = (a[2] XOR a[13] XOR a[22]) XOR (a[0]&b[0] XOR a[0]&c[0] XOR b[0]&c[0])
    #
    # a_new[0] = T1[0] XOR T2[0] (no carry!)
    #
    # So a_new[0] is an EXACT GF(2) function of:
    #   h[0], e[6], e[11], e[25], e[0], f[0], g[0], K[r][0], W[r][0],
    #   a[2], a[13], a[22], a[0], b[0], c[0]
    #
    # That's 15 input bits, degree 2 (from Ch and Maj), NO CARRIES.

    print(f"\nBit 0 of a[r+1]:")
    print(f"  = h[0] XOR e[6] XOR e[11] XOR e[25]")
    print(f"    XOR (e[0]&f[0] XOR ~e[0]&g[0])")
    print(f"    XOR K[r][0] XOR W[r][0]")
    print(f"    XOR (a[2] XOR a[13] XOR a[22])")
    print(f"    XOR (a[0]&b[0] XOR a[0]&c[0] XOR b[0]&c[0])")
    print(f"")
    print(f"  = EXACT degree-2 GF(2) polynomial in 15 state bits + K[r][0] + W[r][0]")
    print(f"  NO carries involved. This is TRUE for ALL messages, ALL rounds.")

    # Verify numerically
    violations = 0
    total = 0
    for seed in range(500):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        states, W = sha256_round_trace(M)

        for r in range(64):
            a, b, c, d, e, f, g, h = states[r]

            # Compute bit 0 of a[r+1] using our formula
            sig1_bit0 = ((e >> 6) ^ (e >> 11) ^ (e >> 25)) & 1
            ch_bit0 = ((e & f) ^ ((~e) & g)) & 1  # bit 0
            sig0_bit0 = ((a >> 2) ^ (a >> 13) ^ (a >> 22)) & 1
            maj_bit0 = ((a & b) ^ (a & c) ^ (b & c)) & 1  # bit 0

            t1_bit0 = (h & 1) ^ sig1_bit0 ^ ch_bit0 ^ (K[r] & 1) ^ (W[r] & 1)
            t2_bit0 = sig0_bit0 ^ maj_bit0
            a_new_bit0_formula = t1_bit0 ^ t2_bit0

            a_new_bit0_actual = states[r+1][0] & 1

            if a_new_bit0_formula != a_new_bit0_actual:
                violations += 1
            total += 1

    print(f"\n  Verification: {violations}/{total} violations")
    if violations == 0:
        print(f"  CONFIRMED: Bit 0 of a[r+1] is an exact degree-2 GF(2) polynomial ✓")
        print(f"  This holds for ALL 64 rounds, ALL messages.")
        print(f"  Bit 0 NEVER sees a carry. It lives entirely in GF(2).")

    # Now: what about bit 1?
    print(f"\n--- Bit 1 of a[r+1]: ---")
    print(f"  Bit 1 has ONE carry-in: from bit 0 of each addition.")
    print(f"  Carry from bit 0 of (x+y) = x[0] AND y[0]")
    print(f"  So bit 1 = XOR of inputs at bit 1 + carry = (degree 2 from bit 1) XOR (degree 2 from bit 0)")
    print(f"  Total degree: 2 (but more terms than bit 0)")

    # And bit k in general:
    print(f"\n--- General: bit k of a[r+1] ---")
    print(f"  Carry into bit k depends on bits 0..k-1 of both operands")
    print(f"  For k additions in T1 and T2:")
    print(f"    Bit 0:  degree 2, depends on 15 state bits (NO carry)")
    print(f"    Bit 1:  degree 2+, depends on ~30 state bits (1 carry level)")
    print(f"    Bit k:  degree ≤ k+2, depends on O(15·k) state bits")
    print(f"    Bit 31: degree ≤ 33, depends on nearly ALL state bits")
    print(f"")
    print(f"  This is the BIT-POSITION COMPLEXITY GRADIENT.")
    print(f"  Low bits = algebraically simple.")
    print(f"  High bits = algebraically complex.")
    print(f"  Rotations SHUFFLE this gradient.")
    print(f"  Each addition RECREATES it.")


def experiment_24_bit0_across_rounds():
    """
    Since bit 0 of a[r+1] is an exact GF(2) polynomial of 15 state bits,
    and the state bits themselves are functions of previous a-values,
    we can trace bit 0 through MULTIPLE rounds.

    After 2 rounds: bit 0 of a[r+2] depends on bit 0 of a[r+1]
    (which we know is degree 2) plus other bits. But those other bits
    include bits from rotated positions (e.g., a[2], a[13], a[22])
    which MAY have carries in their history.

    Key question: Does bit 0's "carry-free" property compose across rounds?
    I.e., is bit 0 of a[r+k] still a "manageable" function for large k?
    """
    print("\n" + "=" * 80)
    print("EXPERIMENT 24: Bit 0 across multiple rounds — does simplicity survive?")
    print("=" * 80)

    # Method: For bit 0 of a[r], compute its dependency on M bits
    # at different round depths. If it stays narrow → simplicity survives.

    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]
    states_base, _ = sha256_round_trace(M)

    print(f"\nFlipping each of 512 message bits, counting which flip bit 0 of a[r]:\n")

    for r in [1, 2, 3, 4, 5, 6, 8, 10, 16, 32, 64]:
        a_base_bit0 = states_base[r][0] & 1

        flips = 0
        for word in range(16):
            for bit in range(32):
                M_flip = list(M)
                M_flip[word] ^= (1 << bit)
                states_flip, _ = sha256_round_trace(M_flip, rounds=r)
                if (states_flip[r][0] & 1) != a_base_bit0:
                    flips += 1

        print(f"  r={r:>2}: bit 0 of a[r] sensitive to {flips:>3}/512 message bits ({flips/512*100:.1f}%)")


if __name__ == "__main__":
    experiment_21_bit0_linearity()
    experiment_22_bit_sensitivity()
    experiment_23_carry_free_bit()
    experiment_24_bit0_across_rounds()
