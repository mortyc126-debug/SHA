"""
Stage 1, Part 8: Derive from formulas. What doesn't fit?

Write the COMPLETE algebra of SHA-256.
Look at what SHOULD work but DOESN'T.
Look at what SHOULDN'T exist but DOES.
The gaps are the fingerprint of what we need to invent.
"""

import random
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, schedule
)


def formulas():
    """
    THE COMPLETE FORMULAS OF SHA-256 (one round)

    Given state (a, b, c, d, e, f, g, h) and constants K[r], W[r]:

      a_new = (h + Sig1(e) + Ch(e,f,g) + K + W) + (Sig0(a) + Maj(a,b,c))    [mod 2^32]
      e_new = d + (h + Sig1(e) + Ch(e,f,g) + K + W)                           [mod 2^32]
      b_new = a, c_new = b, d_new = c
      f_new = e, g_new = f, h_new = g

    Where:
      Sig0(x) = ROTR(x,2) ⊕ ROTR(x,13) ⊕ ROTR(x,22)     [linear over GF(2)]
      Sig1(x) = ROTR(x,6) ⊕ ROTR(x,11) ⊕ ROTR(x,25)     [linear over GF(2)]
      Ch(e,f,g) = (e ∧ f) ⊕ (¬e ∧ g)                      [degree 2 over GF(2)]
      Maj(a,b,c) = (a ∧ b) ⊕ (a ∧ c) ⊕ (b ∧ c)           [degree 2 over GF(2)]

    Using створочное число: e[r] = a[r] + a[r-4] - T2[r-1]
    where T2[r-1] = Sig0(a[r-1]) + Maj(a[r-1], a[r-2], a[r-3])

    FULL SUBSTITUTION (a-only recurrence):
      Let E(r) = a[r] + a[r-4] - Sig0(a[r-1]) - Maj(a[r-1], a[r-2], a[r-3])

      a[r+1] = E(r-3)                                    [= h in original]
             + Sig1(E(r))                                  [Sig1 of e]
             + Ch(E(r), E(r-1), E(r-2))                   [Ch of e, f, g]
             + K[r] + W[r]                                 [constants]
             + Sig0(a[r])                                  [direct]
             + Maj(a[r], a[r-1], a[r-2])                   [direct]

    This is the MASTER EQUATION. Everything follows from this.
    """
    pass


def anomaly_analysis():
    """
    ANOMALY 1: The operation mismatch.

    Sig0, Sig1 are GF(2)-linear (XOR of rotations).
    Ch, Maj are GF(2)-quadratic (AND + XOR).
    The additions (+) are Z/2^32-linear but GF(2)-nonlinear.

    Key: a_new = [GF(2)-LINEAR] + [GF(2)-QUADRATIC] + [CONSTANTS]
                  ↑ Sig0, Sig1      ↑ Ch, Maj           ↑ K, W

    The + here is mod 2^32. If we decompose + = XOR + carry:

    a_new = Sig0(a) ⊕ Maj(a,b,c)                              [T2 XOR part]
          ⊕ h ⊕ Sig1(e) ⊕ Ch(e,f,g) ⊕ K ⊕ W                 [T1 XOR part]
          ⊕ carry_T2(Sig0(a), Maj(a,b,c))                      [T2 carry]
          ⊕ carry_T1(h, Sig1(e), Ch(e,f,g), K, W)             [T1 carries]
          ⊕ carry_final(T1, T2)                                 [final carry]

    The XOR parts are degree ≤ 2 over GF(2). KNOWN.
    The carry parts are degree up to 32 over GF(2). UNKNOWN.

    ANOMALY: The carry depends on LOWER bits of both operands.
    At bit position k:
      carry_k = f(operand1[0..k-1], operand2[0..k-1])

    This function f is the carry CHAIN. It's a specific Boolean function.
    It's not random — it has EXACT structure.

    What IS this structure?
    """
    print("=" * 80)
    print("ANOMALY ANALYSIS: What the formulas predict but we haven't exploited")
    print("=" * 80)

    print("""
    FORMULA DECOMPOSITION:

    a[r+1] = XOR_PART(state[r], K[r], W[r])  ⊕  CARRY_PART(state[r], K[r], W[r])

    Where:
      XOR_PART = Sig0(a) ⊕ Maj(a,b,c) ⊕ h ⊕ Sig1(e) ⊕ Ch(e,f,g) ⊕ K ⊕ W
               = degree 2 GF(2) polynomial in 256 state bits
               = KNOWN, EXACT, computable

      CARRY_PART = XOR of all carry chains from ~6 additions
                 = degree up to 32, but with SPECIFIC structure
                 = carry_chain(x, y)[k] = MAJ(x[k-1], y[k-1], carry[k-1])

    The carry chain is a SEQUENTIAL RECURSION:
      c[0] = 0
      c[k] = MAJ(x[k-1], y[k-1], c[k-1])
      c[k] = x[k-1]·y[k-1] + x[k-1]·c[k-1] + y[k-1]·c[k-1]  (over GF(2))

    THIS IS THE KEY FORMULA. The carry at position k is:
      c[k] = ∨_{j=0}^{k-1} (G[j] · ∏_{i=j+1}^{k-1} P[i])

    Where G[i] = x[i]·y[i] (generate) and P[i] = x[i]⊕y[i] (propagate).

    Written out: c[k] = G[k-1] + P[k-1]·G[k-2] + P[k-1]·P[k-2]·G[k-3] + ...

    This is a SUM OF PRODUCTS where each product is a CHAIN of propagates
    terminated by a generate. The LENGTH of the chain determines the degree.
    """)

    # Verify the carry formula
    print("  Verifying carry chain formula...")
    N = 1000
    violations = 0
    for _ in range(N):
        x = random.randint(0, MASK)
        y = random.randint(0, MASK)
        result = (x + y) & MASK
        xor_part = x ^ y

        # Compute carry chain using formula
        carry = 0
        carry_word = 0
        for k in range(32):
            if carry:
                carry_word |= (1 << k)
            xk = (x >> k) & 1
            yk = (y >> k) & 1
            carry = (xk & yk) | (xk & carry) | (yk & carry)

        # result should equal xor_part XOR carry_word
        if result != (xor_part ^ carry_word):
            violations += 1

    print(f"    x + y = (x ⊕ y) ⊕ carry_chain(x,y): {violations}/{N} violations")

    print("""
    ═══════════════════════════════════════════════════════════════
    WHAT DOESN'T FIT — THE THREE ANOMALIES
    ═══════════════════════════════════════════════════════════════

    ANOMALY 1: CARRY IS BOTH SEQUENTIAL AND ALGEBRAIC
      - Sequentially: c[k] depends on c[k-1] (chain)
      - Algebraically: c[k] = sum of products of G's and P's (polynomial)
      - These are TWO EQUIVALENT descriptions of the SAME object
      - GF(2) uses the algebraic view → degree explosion (degree k at bit k)
      - Sequential view never explodes → always depth 1 recursion
      - NO EXISTING FRAMEWORK uses BOTH views simultaneously
      - The carry is the SIMPLEST POSSIBLE nonlinear recursion (MAJ of 3 bits)
        that GF(2) algebra makes look maximally complex (degree 32)

    ANOMALY 2: SHA-256 HAS TWO "TIMES"
      - Time 1: ROUND index r = 0, 1, ..., 63 (computation progress)
      - Time 2: BIT index k = 0, 1, ..., 31 (carry propagation)
      - Round function advances Time 1 by 1
      - Carry chain advances Time 2 by 1 (within each addition)
      - Rotation COUPLES Time 1 and Time 2 (bit k at round r
        becomes bit (k+rot) at the same round)
      - No framework models two coupled "times" in one system

    ANOMALY 3: THE FUNCTION IS SIMPLER THAN IT LOOKS
      - Over GF(2): degree 32 after 6 rounds → looks maximally complex
      - Over Z/2^32: linear (just additions) → looks trivially simple
      - Neither view is "right" — the function has intermediate complexity
      - ACTUAL information: ~107 carry bits per round (experiment 3)
        out of 224 possible = 48%
      - The function is HALF as complex as GF(2) thinks, and INFINITELY
        more complex than Z/2^32 thinks
      - There should be a description where the complexity is EXACTLY 48%
        — not 0% (Z/2^32) and not 100% (GF(2))
    """)


def anomaly_2_two_times():
    """
    ANOMALY 2 concrete test: SHA-256 has two "times".

    Time 1 (round): a[r+1] = F(a[r..r-7], ...)
    Time 2 (bit/carry): carry[k] = MAJ(x[k-1], y[k-1], carry[k-1])

    If we VIEW the round function at bit level:
      a_new[k] = xor_of_inputs[k] ⊕ carry[k]

    The carry at position k is:
      carry[k] = recursion over positions 0..k-1

    This means: a_new[0] is determined FIRST (no carry).
    Then a_new[1] (needs carry from 0).
    Then a_new[2] (needs carry from 0,1).
    ...

    There is a NATURAL ORDERING within one round:
    bit 0 is "past", bit 31 is "future" (in carry-time).

    Rotation BREAKS this ordering. ROTR(x, 6) sends bit 6 to bit 0.
    Bit 6 (which was "future") becomes bit 0 (which is "past").

    This means: rotation is a TIME MACHINE in carry-time.
    It sends future-bits to the past.

    What if we track this "carry-time" explicitly?
    """
    print("\n" + "=" * 80)
    print("ANOMALY 2: Two times — verify carry propagation ordering")
    print("=" * 80)

    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]
    states, W = sha256_round_trace(M)

    r = 20
    a, b, c, d, e, f, g, h = states[r]

    # Compute T1 step by step, tracking carry at each bit
    operands_T1 = [h, Sig1(e), Ch(e, f, g), K[r], W[r]]
    T2 = (Sig0(a) + Maj(a, b, c)) & MASK

    # Sequential addition of T1 components
    running = 0
    all_carries = []
    for op in operands_T1:
        carry_chain = []
        c_bit = 0
        for k in range(32):
            carry_chain.append(c_bit)
            xk = (running >> k) & 1
            yk = (op >> k) & 1
            c_bit = (xk & yk) | (xk & c_bit) | (yk & c_bit)
        all_carries.append(carry_chain)
        running = (running + op) & MASK

    T1 = running
    # Final addition: T1 + T2
    final_carry = []
    c_bit = 0
    for k in range(32):
        final_carry.append(c_bit)
        xk = (T1 >> k) & 1
        yk = (T2 >> k) & 1
        c_bit = (xk & yk) | (xk & c_bit) | (yk & c_bit)
    all_carries.append(final_carry)

    # The "carry time" at each bit position:
    # = how many carry-levels deep is this bit's value?
    # Bit 0: always 0 carry levels (no carry-in)
    # Bit k: up to k carry levels

    # But MORE importantly: where does the carry chain BREAK?
    # A "break" = position where carry = 0 regardless of lower bits.
    # This happens when G[k] = P[k] = 0 (both operand bits = 0).
    # At a break, bit k+1 becomes independent of bits 0..k.

    # Count break points in each addition
    print(f"\n  Carry chain breaks in round {r} additions:")
    add_names = ['h+Sig1', '+Ch', '+K', '+W', 'T1+T2']

    running = 0
    operands = [h, Sig1(e), Ch(e, f, g), (K[r] + W[r]) & MASK]
    for name, op in zip(['h+Sig1(e)', '+(Ch)', '+(K+W)'], [Sig1(e), Ch(e, f, g), (K[r] + W[r]) & MASK]):
        running_before = running
        running = (running + op) & MASK

        breaks = 0
        for k in range(1, 32):
            xk = (running_before >> k) & 1
            yk = (op >> k) & 1
            # Break if both bits are 0 AND carry is 0
            # Actually, a carry chain breaks when carry_out = 0
            # carry_out[k] = G[k] | (P[k] & carry_in[k])
            # If G[k] = 0 and P[k] = 0 → carry_out = 0 → break
            if xk == 0 and yk == 0:
                breaks += 1

        print(f"    {name:>12}: {breaks}/31 break points ({breaks/31*100:.0f}%)")

    # The average: about 25% of positions are break points.
    # At a break, the carry chain RESETS → independence.
    # This means: bits are NOT all coupled through carries.
    # There are ~8 INDEPENDENT SEGMENTS of the carry chain.

    print(f"""
    ═══════════════════════════════════════════════════════════════
    THE "ANTI-PARTICLE": What the formulas predict should NOT exist
    ═══════════════════════════════════════════════════════════════

    The formulas say:
      carry[k] depends on ALL bits 0..k-1

    But in PRACTICE:
      carry chains BREAK at ~25% of positions
      → Average carry chain length ≈ 4 bits
      → Bits are grouped into ~8 INDEPENDENT segments

    This means: GF(2) degree = 32 is the WORST CASE.
    AVERAGE degree ≈ 4-5 (length of typical carry chain).

    The "effective algebra" of SHA-256 is NOT degree-32 polynomials.
    It's degree-4 polynomials CONNECTED by carry propagation.

    This is the GAP:
      - GF(2) says: degree 32 → 2^32 monomials → impossible
      - Reality says: degree ~4 connected by chains → much simpler
      - Neither GF(2) nor Z/2^32 captures this intermediate structure

    THE ANTI-PARTICLE:
      A framework where "degree" is not a global number (32)
      but a LOCAL property that varies by position (0 at breaks, 4 on average).
      The total complexity = PRODUCT of segment complexities, not the
      global degree.

      If segments have degree d and there are s segments:
        GF(2) says: complexity = 2^(d*s) = 2^(4*8) = 2^32
        Segmented view: complexity = 2^d * s = 2^4 * 8 = 128
                                              or (2^d)^s with structure

      The truth is between these extremes.
    """)


def measure_carry_segments():
    """
    Measure carry chain segment statistics across many rounds and messages.
    """
    print("\n" + "=" * 80)
    print("MEASUREMENT: Carry chain segment statistics")
    print("=" * 80)

    N = 500
    all_segment_lengths = []

    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        states, W = sha256_round_trace(M)

        for r in range(64):
            a, b, c, d, e, f, g, h = states[r]
            T1 = (h + Sig1(e) + Ch(e, f, g) + K[r] + W[r]) & MASK
            T2 = (Sig0(a) + Maj(a, b, c)) & MASK

            # Carry chain for T1 + T2
            c_bit = 0
            current_len = 0
            for k in range(32):
                xk = (T1 >> k) & 1
                yk = (T2 >> k) & 1
                new_carry = (xk & yk) | (xk & c_bit) | (yk & c_bit)

                if new_carry:
                    current_len += 1
                else:
                    if current_len > 0:
                        all_segment_lengths.append(current_len)
                    current_len = 0
                c_bit = new_carry

            if current_len > 0:
                all_segment_lengths.append(current_len)

    # Statistics
    avg = sum(all_segment_lengths) / len(all_segment_lengths)
    max_len = max(all_segment_lengths)

    # Distribution
    dist = {}
    for l in all_segment_lengths:
        dist[l] = dist.get(l, 0) + 1

    print(f"  Total segments measured: {len(all_segment_lengths)}")
    print(f"  Average segment length: {avg:.2f}")
    print(f"  Max segment length: {max_len}")
    print(f"\n  Distribution of carry chain segment lengths:")
    for length in sorted(dist.keys())[:15]:
        count = dist[length]
        pct = count / len(all_segment_lengths) * 100
        bar = '█' * int(pct)
        print(f"    len={length:>2}: {count:>6} ({pct:>5.1f}%) {bar}")

    # P(segment ≥ k)
    print(f"\n  Survival function P(segment ≥ k):")
    total = len(all_segment_lengths)
    for k in [1, 2, 3, 4, 5, 8, 10, 16, 20, 31]:
        count = sum(1 for l in all_segment_lengths if l >= k)
        print(f"    P(len ≥ {k:>2}) = {count/total:.4f}")


if __name__ == "__main__":
    anomaly_analysis()
    anomaly_2_two_times()
    measure_carry_segments()
