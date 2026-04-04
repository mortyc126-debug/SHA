"""
Step 2: The Coin Mechanism — WHY does each monomial survive or die?

The "coin flip" for monomial m at round r is not random.
It's determined by the algebraic interaction of ANFs during the round.

When we compute a_new = T1 + T2 (mod 2^n):
  a_new = T1 XOR T2 XOR carry(T1, T2)

In ANF terms:
  ANF(a_new)[m] = ANF(T1)[m] ⊕ ANF(T2)[m] ⊕ ANF(carry)[m]

Monomial m SURVIVES in a_new iff this XOR = 1.
Monomial m DIES iff this XOR = 0.

The "coin" = ANF(T1)[m] ⊕ ANF(T2)[m] ⊕ ANF(carry)[m].

Let's compute this exactly and see WHAT DETERMINES the coin.
"""

import numpy as np
from collections import Counter, defaultdict
from step1_monomial_genealogy import (
    N, MASK, N_MSG, N_INPUT, N_TOTAL, IV, K,
    make_constant, make_input_bit, make_const_word,
    tt_xor, tt_and, tt_not,
    mobius_transform, monomial_set,
    symbolic_rotr, symbolic_shr, symbolic_xor_word,
    symbolic_ch, symbolic_maj, symbolic_sigma0, symbolic_sigma1,
    symbolic_add, symbolic_add_no_carry
)


def analyze_coin_mechanism(R_max=6):
    """
    For each round, decompose the ANF of a_new[bit] into contributions
    from T1, T2, and carry. Determine what makes each monomial live or die.
    """

    # Initialize state
    state_words = []
    for reg in range(8):
        word = []
        for bit in range(N):
            val = (IV[reg] >> bit) & 1
            word.append(make_constant(val))
        state_words.append(word)

    W_msg = []
    for w in range(N_MSG):
        word = []
        for bit in range(N):
            var_idx = w * N + bit
            word.append(make_input_bit(var_idx))
        W_msg.append(word)

    a, b, c, d, e, f, g, h = state_words

    for r in range(R_max):
        w_r = W_msg[r] if r < len(W_msg) else [make_constant(0)] * N
        k_r = make_const_word(K[r % len(K)])

        # Intermediates
        sig1_e = symbolic_sigma1(e)
        ch_efg = symbolic_ch(e, f, g)
        sig0_a = symbolic_sigma0(a)
        maj_abc = symbolic_maj(a, b, c)

        # T1 = h + sig1(e) + ch(e,f,g) + K[r] + W[r]
        t1_1, _ = symbolic_add(h, sig1_e)
        t1_2, _ = symbolic_add(t1_1, ch_efg)
        t1_3, _ = symbolic_add(t1_2, k_r)
        T1, _ = symbolic_add(t1_3, w_r)

        # T2 = sig0(a) + maj(a,b,c)
        T2, _ = symbolic_add(sig0_a, maj_abc)

        # a_new = T1 + T2 (with carry)
        a_new, carry_a = symbolic_add(T1, T2)

        # e_new = d + T1
        e_new, carry_e = symbolic_add(d, T1)

        # ============================================================
        # THE COIN: for each bit, decompose a_new ANF
        # ============================================================
        print(f"\n{'='*70}")
        print(f"ROUND {r}: The Coin Mechanism")
        print(f"{'='*70}")

        for bit in [0, 1, N-1]:
            # ANFs of the three components
            anf_T1 = mobius_transform(T1[bit])
            anf_T2 = mobius_transform(T2[bit])
            anf_a_new = mobius_transform(a_new[bit])

            # For bit > 0, carry from previous bit matters
            if bit > 0:
                anf_carry = mobius_transform(carry_a[bit - 1])
            else:
                anf_carry = np.zeros(N_TOTAL, dtype=np.uint8)

            # The algebraic identity:
            # a_new[bit] = T1[bit] XOR T2[bit] XOR carry[bit]
            # Verify:
            reconstructed = anf_T1 ^ anf_T2 ^ anf_carry
            # This should equal anf_a_new... but only for truth tables!
            # For ANFs, the relationship is different because ANF(f XOR g) = ANF(f) XOR ANF(g)
            # but ANF(f AND g) ≠ ANF(f) AND ANF(g).
            # However, since a_new[bit] = T1[bit] ⊕ T2[bit] ⊕ carry[bit] at truth table level,
            # and Möbius transform is linear over GF(2):
            # ANF(a_new) = ANF(T1) ⊕ ANF(T2) ⊕ ANF(carry)

            match = np.array_equal(reconstructed, anf_a_new)

            monos_T1 = monomial_set(anf_T1)
            monos_T2 = monomial_set(anf_T2)
            monos_carry = monomial_set(anf_carry)
            monos_result = monomial_set(anf_a_new)

            print(f"\n  a_new[{bit}]: T1⊕T2⊕carry = a_new? {match}")
            print(f"    |T1| = {len(monos_T1)}, |T2| = {len(monos_T2)}, "
                  f"|carry| = {len(monos_carry)}, |a_new| = {len(monos_result)}")

            # Categorize each monomial in a_new by its source
            # A monomial m is in a_new iff it's in ODD number of {T1, T2, carry}
            only_T1 = monos_T1 - monos_T2 - monos_carry  # in T1 only
            only_T2 = monos_T2 - monos_T1 - monos_carry  # in T2 only
            only_carry = monos_carry - monos_T1 - monos_T2  # in carry only
            all_three = monos_T1 & monos_T2 & monos_carry  # in all three

            # Monomials in exactly two of three → cancelled (XOR = 0)
            t1_and_t2 = (monos_T1 & monos_T2) - monos_carry
            t1_and_carry = (monos_T1 & monos_carry) - monos_T2
            t2_and_carry = (monos_T2 & monos_carry) - monos_T1

            print(f"    SOURCES of monomials in a_new[{bit}]:")
            print(f"      Only T1:    {len(only_T1):>6}  ({len(only_T1)/max(len(monos_result),1)*100:>5.1f}%)")
            print(f"      Only T2:    {len(only_T2):>6}  ({len(only_T2)/max(len(monos_result),1)*100:>5.1f}%)")
            print(f"      Only carry: {len(only_carry):>6}  ({len(only_carry)/max(len(monos_result),1)*100:>5.1f}%)")
            print(f"      All three:  {len(all_three):>6}  ({len(all_three)/max(len(monos_result),1)*100:>5.1f}%)")
            print(f"    CANCELLED (in exactly 2 of 3):")
            print(f"      T1∩T2:      {len(t1_and_t2):>6}")
            print(f"      T1∩carry:   {len(t1_and_carry):>6}")
            print(f"      T2∩carry:   {len(t2_and_carry):>6}")

            # Degree analysis of each source
            if len(monos_result) > 0 and len(monos_result) < 40000:
                print(f"    DEGREE by source:")

                for name, mono_set in [("Only_T1", only_T1), ("Only_T2", only_T2),
                                        ("Only_carry", only_carry), ("All_three", all_three)]:
                    if mono_set:
                        degrees = [bin(m).count('1') for m in mono_set]
                        avg_deg = sum(degrees) / len(degrees)
                        max_deg = max(degrees)
                        print(f"      {name:>12}: count={len(mono_set):>5}, "
                              f"avg_deg={avg_deg:.1f}, max_deg={max_deg}")

        # ============================================================
        # KEY ANALYSIS: Carry contribution vs T1/T2 as rounds increase
        # ============================================================
        print(f"\n  CARRY CONTRIBUTION SUMMARY (round {r}):")
        for bit in range(N):
            anf_T1 = mobius_transform(T1[bit])
            anf_T2 = mobius_transform(T2[bit])
            anf_a_new = mobius_transform(a_new[bit])

            if bit > 0:
                anf_carry = mobius_transform(carry_a[bit - 1])
            else:
                anf_carry = np.zeros(N_TOTAL, dtype=np.uint8)

            monos_carry = monomial_set(anf_carry)
            monos_result = monomial_set(anf_a_new)

            # What fraction of a_new comes from carry?
            carry_in_result = monos_carry & monos_result
            carry_exclusive = monos_carry - monomial_set(anf_T1) - monomial_set(anf_T2)
            carry_contribution = len(carry_exclusive & monos_result)

            total = max(len(monos_result), 1)
            print(f"    bit {bit}: |carry|={len(monos_carry):>6}, "
                  f"carry_exclusive_in_result={carry_contribution:>6} "
                  f"({carry_contribution/total*100:>5.1f}% of output)")

        # Update state for next round
        h = [g[k].copy() for k in range(N)]
        g = [f[k].copy() for k in range(N)]
        f = [e[k].copy() for k in range(N)]
        e_new_word, _ = symbolic_add(d, T1)
        e = e_new_word
        d = [c[k].copy() for k in range(N)]
        c = [b[k].copy() for k in range(N)]
        b = [a[k].copy() for k in range(N)]
        a = a_new

    # ============================================================
    # FINAL: The anatomy of "randomness"
    # ============================================================
    print(f"\n{'='*70}")
    print(f"ANATOMY OF RANDOMNESS")
    print(f"{'='*70}")
    print(f"""
The 'coin flip' for monomial m at each round is:

    survive(m) = ANF(T1)[m] ⊕ ANF(T2)[m] ⊕ ANF(carry)[m]

This is NOT random. It is a DETERMINISTIC function of:
  1. The current state polynomials (which determine T1, T2)
  2. The carry propagation (which is a product of AND operations)

The reason it LOOKS like 50/50:
  - T1 and T2 are HIGH-DEGREE polynomials (after a few rounds)
  - Whether a specific monomial appears in a high-degree polynomial
    depends on GLOBAL interactions across all variables
  - These interactions are so complex that the outcome is
    PSEUDORANDOM — not random, but unpredictable without
    computing the full polynomial

The 'coin' IS the round function. Understanding it = understanding SHA-256.
""")


if __name__ == "__main__":
    analyze_coin_mechanism(R_max=6)
