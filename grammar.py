"""
SHA-256 GRAMMAR EXTRACTION.

SHA-256 speaks a language. We need its grammar.

What we know:
- Alphabet: {0, 1} (bits)
- Words: 32-bit values
- Sentences: 64-round trajectories
- Grammar rule 1 (93.7%): fixed skeleton (rotations, shifts, Ch, Maj XOR parts)
- Grammar rule 2 (6.3%): P-mask determines carry coupling

The grammar is HOW P-masks at round r constrain P-masks at round r+1.

A "grammatical sentence" = a valid SHA-256 trajectory.
Not all P-mask sequences are grammatical — only those arising from
actual SHA-256 computations on some message M.

The GRAMMAR = constraints on P-mask sequences that make them valid.

These constraints ARE the SHA-256 algorithm, written in P-mask language.
Extracting them = translating SHA-256 from value-space to P-mask-space.
"""

import random
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, schedule
)


def extract_p_masks(M):
    """
    Extract ALL P-masks from a complete SHA-256 computation.

    For each of the ~6 additions per round, record the P-mask.
    Returns list of 64 tuples, each containing P-masks for that round.
    """
    states, W = sha256_round_trace(M)
    all_pmasks = []

    for r in range(64):
        a, b, c, d, e, f, g, h = states[r]

        # T1 = h + Sig1(e) + Ch(e,f,g) + K[r] + W[r]
        sig1e = Sig1(e)
        che = Ch(e, f, g)

        # Addition 1: h + Sig1(e)
        sum1 = (h + sig1e) & MASK
        p1 = h ^ sig1e

        # Addition 2: sum1 + Ch
        sum2 = (sum1 + che) & MASK
        p2 = sum1 ^ che

        # Addition 3: sum2 + (K[r] + W[r])
        kw = (K[r] + W[r]) & MASK
        T1 = (sum2 + kw) & MASK
        p3 = sum2 ^ kw

        # T2 = Sig0(a) + Maj(a,b,c)
        sig0a = Sig0(a)
        maja = Maj(a, b, c)
        T2 = (sig0a + maja) & MASK
        p4 = sig0a ^ maja

        # a_new = T1 + T2
        p5 = T1 ^ T2

        # e_new = d + T1
        p6 = d ^ T1

        all_pmasks.append((p1, p2, p3, p4, p5, p6))

    return all_pmasks, states, W


def experiment_pmask_constraints():
    """
    THE GRAMMAR: How does P-mask at round r+1 depend on P-masks at round r?

    P5[r] = T1[r] XOR T2[r]
    T1[r] = f(state[r], K[r], W[r])
    T2[r] = g(state[r])
    state[r+1] = derived from T1[r], T2[r], and shifts of state[r]

    So P5[r+1] = T1[r+1] XOR T2[r+1]
               = f(state[r+1], K[r+1], W[r+1]) XOR g(state[r+1])

    state[r+1] depends on state[r] through T1[r] and T2[r].
    And T1[r], T2[r] are related to P5[r] by: T1+T2 = result, T1 XOR T2 = P5.

    The GRAMMAR RULE: P5[r+1] is a function of P5[r] (and other P-masks)
    AND the state values (which carry the context).

    Let's measure: given P5[r], how much information does it give about P5[r+1]?
    """
    print("=" * 80)
    print("GRAMMAR EXTRACTION: P-mask → P-mask transition rules")
    print("=" * 80)

    N = 500
    # For each (P5[r], P5[r+1]) pair, compute mutual information
    # Simplified: compute XOR distance between P5[r] and P5[r+1]
    # and compare with random baseline

    xor_dists_consecutive = []
    xor_dists_random = []

    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        pmasks, _, _ = extract_p_masks(M)

        for r in range(8, 60):
            p5_r = pmasks[r][4]     # P5 at round r
            p5_r1 = pmasks[r+1][4]  # P5 at round r+1
            p5_rand = pmasks[(r+17) % 64][4]  # Some distant round

            xor_dists_consecutive.append(bin(p5_r ^ p5_r1).count('1'))
            xor_dists_random.append(bin(p5_r ^ p5_rand).count('1'))

    avg_cons = sum(xor_dists_consecutive) / len(xor_dists_consecutive)
    avg_rand = sum(xor_dists_random) / len(xor_dists_random)

    print(f"  E[HW(P5[r] XOR P5[r+1])] = {avg_cons:.3f} (consecutive rounds)")
    print(f"  E[HW(P5[r] XOR P5[r+k])] = {avg_rand:.3f} (distant rounds)")
    print(f"  Random expectation: 16.000")
    print(f"  Difference: {abs(avg_cons - avg_rand):.3f}")


def experiment_pmask_joint():
    """
    Instead of P5 alone, look at ALL 6 P-masks together.

    The 6 P-masks per round = 6 × 32 = 192 bits.
    These 192 bits encode the COMPLETE nonlinear state of one round.

    How much do the 192 bits at round r determine the 192 bits at round r+1?
    """
    print("\n" + "=" * 80)
    print("GRAMMAR: Joint P-mask (192 bits) round-to-round")
    print("=" * 80)

    N = 300
    # Compute: how many of the 192 output P-mask bits at r+1
    # can be predicted from the 192 input P-mask bits at r?
    # (Plus the known skeleton)

    # Method: collect (pmask[r], pmask[r+1]) pairs, compute GF(2) rank
    # of the mapping. If rank < 192 → P-masks at r+1 are partially
    # determined by P-masks at r.

    def flatten_pmask(pm_tuple):
        """Flatten 6 × 32-bit P-masks to 192-bit vector."""
        bits = []
        for pm in pm_tuple:
            for k in range(32):
                bits.append((pm >> k) & 1)
        return bits

    # Compute Jacobian: ∂(pmask[r+1])/∂(pmask[r])
    # But P-masks aren't free variables — they're derived from state.
    # Instead: just measure correlation.

    # Bit-level: for each of 192 output bits, what fraction of messages
    # have that bit = 1?
    bit_freqs_r = [0] * 192
    bit_freqs_r1 = [0] * 192
    joint_11 = [[0]*192 for _ in range(192)]  # P(bit_i[r]=1 AND bit_j[r+1]=1)

    target_r = 20
    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        pmasks, _, _ = extract_p_masks(M)

        bits_r = flatten_pmask(pmasks[target_r])
        bits_r1 = flatten_pmask(pmasks[target_r + 1])

        for i in range(192):
            if bits_r[i]:
                bit_freqs_r[i] += 1
            if bits_r1[i]:
                bit_freqs_r1[i] += 1
            if bits_r[i] and bits_r1[i]:
                joint_11[i][i] += 1  # Just diagonal for speed

    # Check diagonal correlations
    print(f"  Diagonal P-mask correlations (bit i at round r vs bit i at round r+1):")
    max_corr = 0
    max_corr_bit = -1
    significant = 0
    for i in range(192):
        p_r = bit_freqs_r[i] / N
        p_r1 = bit_freqs_r1[i] / N
        p_11 = joint_11[i][i] / N
        expected = p_r * p_r1
        if abs(p_r - 0.5) < 0.01 and abs(p_r1 - 0.5) < 0.01:
            corr = (p_11 - expected)
            if abs(corr) > max_corr:
                max_corr = abs(corr)
                max_corr_bit = i
            if abs(corr) > 0.02:
                significant += 1

    print(f"  Max |correlation excess| = {max_corr:.4f} at bit {max_corr_bit}")
    print(f"  Significant correlations (|excess|>0.02): {significant}/192")
    print(f"  (Random: max excess ≈ 0.03 from {N} samples)")


def experiment_grammar_rules():
    """
    Extract EXPLICIT grammar rules.

    For SHA-256: state[r+1] = F(state[r], K[r], W[r]).
    In P-mask language:
      T1[r] and T2[r] are determined by state[r] and W[r].
      P5[r] = T1[r] XOR T2[r].
      state[r+1][a] = T1[r] + T2[r] = P5[r] XOR carry_chain(T1[r] AND T2[r], P5[r]).

    The grammar rule connecting P-masks across rounds:
      P5[r+1] depends on state[r+1] which depends on (T1[r], T2[r])
      which is related to P5[r] by: knowing P5[r] and the result = T1+T2
      determines the carry vector, which determines T1 AND T2 = generate mask.

    Specifically: given a_new = T1 + T2 and P5 = T1 XOR T2:
      carry_vector = a_new XOR P5
      generate_mask = (a_new - P5) >> 1... no, more complex.

    Actually: a_new[k] = P5[k] XOR carry[k]
    So: carry[k] = a_new[k] XOR P5[k]

    And: T1[k] = (P5[k] AND mask[k]) OR generate[k]... complicated.

    Let me just directly measure: given (a_new, P5) can we recover (T1, T2)?
    """
    print("\n" + "=" * 80)
    print("GRAMMAR RULES: Can (a_new, P5) recover (T1, T2)?")
    print("=" * 80)

    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]
    states, W = sha256_round_trace(M)

    # For each round: given a_new and P5, recover T1 and T2
    recovery_ok = 0
    total = 0

    for r in range(64):
        a, b, c, d, e, f, g, h = states[r]
        T1 = (h + Sig1(e) + Ch(e, f, g) + K[r] + W[r]) & MASK
        T2 = (Sig0(a) + Maj(a, b, c)) & MASK
        a_new = states[r+1][0]
        P5 = T1 ^ T2

        # carry[k] = a_new[k] XOR P5[k]
        carry = a_new ^ P5

        # From P5 and carry, can we get T1 and T2?
        # P5[k] = T1[k] XOR T2[k]
        # carry[k] was generated by bits at position k-1
        # generate[k-1] = T1[k-1] AND T2[k-1]
        # propagate[k-1] = P5[k-1]

        # At bit 0: carry = 0 always. T1[0] XOR T2[0] = P5[0]. T1[0] AND T2[0] = ?
        # We need ONE more bit of info per position to separate T1 from T2.

        # But we HAVE that info: T2 = Sig0(a[r]) + Maj(a[r], b[r], c[r])
        # And a[r] = state[r+1][1] (b of next state), which we know.
        # b[r] = state[r+1][2], c[r] = state[r+1][3].

        # So T2 IS KNOWN from the next state's b,c,d registers!
        T2_recovered = (Sig0(states[r+1][1]) + Maj(states[r+1][1], states[r+1][2], states[r+1][3])) & MASK
        T1_recovered = (a_new - T2_recovered) & MASK

        if T1_recovered == T1 and T2_recovered == T2:
            recovery_ok += 1
        total += 1

    print(f"  Recovery of (T1, T2) from next state: {recovery_ok}/{total}")

    # So: at each round, T2 is KNOWN (from b,c,d of next state).
    # This means: P5 + knowledge of T2 → T1 = a_new - T2.
    # And T1 = h + Sig1(e) + Ch(e,f,g) + K + W.
    # Since h,e,f,g are known (from shifts of next state):
    #   W[r] = T1 - h - Sig1(e) - Ch(e,f,g) - K[r]
    # For r < 16: W[r] = M[r] → message recovered!
    # For r >= 16: W[r] = schedule(M) → constrains M.

    print(f"\n  CONSEQUENCE: Given the FULL trajectory state[0..64],")
    print(f"  every P-mask, T1, T2, and W[r] is EXACTLY recoverable.")
    print(f"  The trajectory IS the message in a different encoding.")
    print(f"")
    print(f"  The GRAMMAR is: the constraints on W[16..63] = schedule(W[0..15]).")
    print(f"  In P-mask language: P-masks at rounds 16-63 must be CONSISTENT")
    print(f"  with P-masks at rounds 0-15 through the message schedule.")
    print(f"")
    print(f"  THE SCHEDULE IS THE GRAMMAR.")
    print(f"  The round function is the alphabet.")
    print(f"  The schedule constrains which 'sentences' are valid.")


def experiment_schedule_in_pmask():
    """
    The message schedule W[16..63] = f(W[0..15]) constrains the P-masks.

    In value space: W[i] = sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16]
    This is a hard constraint on 48 of the 64 W-values.

    In P-mask space: what does this constraint look like?
    Does it constrain P-masks in a way that's simpler than in value space?

    Key test: For two different W[0..15] giving the same W[16] (or same P5[16]),
    how constrained are the later P-masks?
    """
    print("\n" + "=" * 80)
    print("THE GRAMMAR = SCHEDULE: How schedule constrains P-masks")
    print("=" * 80)

    # The schedule links W[r] values across rounds.
    # In the P-mask framework:
    # P5[r] = T1[r] XOR T2[r]
    # T1[r] depends on W[r], and W[r] for r >= 16 depends on W[0..15].
    # So: P5[r] for r >= 16 is constrained by W[0..15] through the schedule.

    # Measure: how many "free" P-masks are there?
    # W[0..15] = 16 words = 512 bits of freedom.
    # P5[0..63] = 64 words = 2048 bits.
    # So 2048 - 512 = 1536 bits of P5 are CONSTRAINED by the schedule.
    # (In addition to the constraints from state evolution.)

    # But P5 isn't determined by W alone — it also depends on state.
    # And state depends on ALL previous W values.
    # So the constraint structure is complex.

    # Let's measure: for varying ONLY W[0], how many P5 values change?
    N = 500
    rng = random.Random(42)
    base_msg = [rng.randint(0, MASK) for _ in range(16)]

    base_pmasks, _, _ = extract_p_masks(base_msg)

    changes_per_round = [0] * 64
    for trial in range(N):
        msg = list(base_msg)
        msg[0] = random.Random(trial).randint(0, MASK)
        pmasks, _, _ = extract_p_masks(msg)

        for r in range(64):
            if pmasks[r][4] != base_pmasks[r][4]:
                changes_per_round[r] += 1

    print(f"  Varying W[0] only: fraction of messages that change P5[r]:")
    for r in [0, 1, 4, 8, 12, 15, 16, 20, 32, 48, 63]:
        frac = changes_per_round[r] / N
        print(f"    r={r:>2}: {frac:.3f}")

    # T_DEP tells us: e_r depends on W[0..r-2].
    # So state at round r depends on W[0..r-2].
    # P5[r] depends on state[r] and W[r].
    # For r >= 2: P5[r] changes when W[0] changes (through state dependency).
    # For r = 0: P5[0] doesn't depend on W[0]... wait, T1[0] includes W[0]!
    # Actually: T1[0] = h_iv + Sig1(e_iv) + Ch(e_iv,f_iv,g_iv) + K[0] + W[0]
    # So P5[0] = (T1[0]) XOR T2[0], and T1[0] depends on W[0].
    # → P5[0] DOES depend on W[0]. Every round's P5 depends on W[0].


if __name__ == "__main__":
    experiment_pmask_constraints()
    experiment_pmask_joint()
    experiment_grammar_rules()
    experiment_schedule_in_pmask()
