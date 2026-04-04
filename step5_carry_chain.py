"""
Step 5: The Carry Chain — building 100% prediction bit by bit

We know:
  a_new[0] = T1[0] ⊕ T2[0]                              → 100% (carry-free)
  a_new[1] = T1[1] ⊕ T2[1] ⊕ carry[1]
  a_new[2] = T1[2] ⊕ T2[2] ⊕ carry[2]
  a_new[3] = T1[3] ⊕ T2[3] ⊕ carry[3]

Where the carry chain is:
  carry[0] = 0
  carry[k] = MAJ(T1[k-1], T2[k-1], carry[k-1])
           = T1[k-1]·T2[k-1] ⊕ T1[k-1]·carry[k-1] ⊕ T2[k-1]·carry[k-1]

Each carry = known function of previous carries + T1, T2.
So carry is a CHAIN of degree multiplications.

Question: What is the COST of understanding each bit?
  - Bit 0: FREE (T1⊕T2)
  - Bit 1: cost = compute T1[0]·T2[0] (one AND = degree doubling)
  - Bit 2: cost = compute MAJ(T1[1], T2[1], carry[1]) (chain deepens)
  - Bit k: cost = k ANDs deep (degree ≈ 2^k)

This shows exactly WHERE the complexity lives and HOW it grows.
"""

import numpy as np
from collections import Counter
from step1_monomial_genealogy import (
    N, MASK, N_MSG, N_INPUT, N_TOTAL, IV, K,
    make_constant, make_input_bit, make_const_word,
    tt_xor, tt_and, tt_not,
    mobius_transform, monomial_set,
    symbolic_rotr, symbolic_shr, symbolic_xor_word,
    symbolic_ch, symbolic_maj, symbolic_sigma0, symbolic_sigma1,
    symbolic_add, symbolic_add_no_carry
)
from step3_breaking_equilibrium import compute_round_state, degree_of


def analyze_carry_chain(R):
    """Analyze the carry chain at round R, bit by bit."""
    T1, T2, carry_a, a_new, e_new, carry_e = compute_round_state(R)

    print(f"\n{'='*80}")
    print(f"CARRY CHAIN ANATOMY — Round {R}")
    print(f"{'='*80}")

    # Rebuild carry chain from T1, T2
    carries = [make_constant(0)]  # carry[0] = 0

    for k in range(1, N):
        # carry[k] = MAJ(T1[k-1], T2[k-1], carry[k-1])
        t1_prev = T1[k-1]
        t2_prev = T2[k-1]
        c_prev = carries[k-1]

        new_carry = tt_xor(tt_xor(
            tt_and(t1_prev, t2_prev),
            tt_and(t1_prev, c_prev)),
            tt_and(t2_prev, c_prev))
        carries.append(new_carry)

    # Verify carries match
    for k in range(N):
        original = carry_a[k-1] if k > 0 else make_constant(0)
        if k == 0:
            assert np.array_equal(carries[0], make_constant(0))
        else:
            rebuilt_anf = mobius_transform(carries[k])
            original_anf = mobius_transform(carry_a[k-1])
            match = np.array_equal(rebuilt_anf, original_anf)
            if not match:
                print(f"  WARNING: carry[{k}] mismatch!")

    # For each bit, show the cumulative prediction model
    print(f"\n  CUMULATIVE PREDICTION MODEL:")
    print(f"  {'Bit':>4} {'Model':>30} {'Accuracy':>9} {'Carry ANF deg':>14} "
          f"{'Carry #mono':>12} {'Total #mono':>12}")

    cumulative_accuracy = []

    for bit in range(N):
        # Model: a_new[bit] = T1[bit] ⊕ T2[bit] ⊕ carry[bit]
        t1_bit = T1[bit]
        t2_bit = T2[bit]
        carry_bit = carries[bit]

        # Prediction using T1⊕T2⊕carry (should be 100%)
        pred = tt_xor(tt_xor(t1_bit, t2_bit), carry_bit)
        actual = a_new[bit]
        acc = np.mean(mobius_transform(pred) == mobius_transform(actual))

        # ANF analysis of carry
        anf_carry = mobius_transform(carry_bit)
        carry_monos = monomial_set(anf_carry)
        carry_degrees = [degree_of(m) for m in carry_monos] if carry_monos else [0]

        # ANF of result
        anf_result = mobius_transform(actual)
        result_monos = monomial_set(anf_result)

        if bit == 0:
            model_str = "T1⊕T2 (no carry)"
        else:
            model_str = f"T1⊕T2⊕carry[{bit}]"

        max_carry_deg = max(carry_degrees) if carry_monos else 0
        avg_carry_deg = sum(carry_degrees)/len(carry_degrees) if carry_degrees else 0

        print(f"  {bit:>4} {model_str:>30} {acc*100:>8.1f}% "
              f"  deg≤{max_carry_deg:<3} avg={avg_carry_deg:.1f}"
              f"  {len(carry_monos):>10}"
              f"  {len(result_monos):>10}")

    # ================================================================
    # The KEY: How does carry ANF complexity grow with bit position?
    # ================================================================
    print(f"\n  CARRY COMPLEXITY GROWTH:")
    print(f"  {'Bit':>4} {'#Monomials':>11} {'Max Degree':>11} {'Avg Degree':>11} "
          f"{'Density':>8}")

    for k in range(N):
        anf = mobius_transform(carries[k])
        monos = monomial_set(anf)
        if monos:
            degrees = [degree_of(m) for m in monos]
            max_deg = max(degrees)
            avg_deg = sum(degrees) / len(degrees)
            density = len(monos) / N_TOTAL * 100
        else:
            max_deg = 0
            avg_deg = 0
            density = 0

        print(f"  {k:>4} {len(monos):>11} {max_deg:>11} {avg_deg:>11.1f} {density:>7.1f}%")

    # ================================================================
    # Degree distribution of each carry level
    # ================================================================
    print(f"\n  CARRY DEGREE DISTRIBUTION:")
    for k in range(N):
        anf = mobius_transform(carries[k])
        monos = monomial_set(anf)
        if not monos:
            print(f"  carry[{k}]: (empty)")
            continue

        deg_dist = Counter(degree_of(m) for m in monos)
        total = len(monos)
        degs_str = ", ".join(f"d{d}={c}" for d, c in sorted(deg_dist.items()) if c > 0)
        print(f"  carry[{k}]: {total} monomials — {degs_str}")

    # ================================================================
    # The critical question: carry as function of T1, T2
    # ================================================================
    print(f"\n  CARRY DECOMPOSITION into T1, T2 components:")
    print(f"  carry[1] = T1[0] AND T2[0]")

    if N > 1:
        # Verify carry[1] = T1[0] AND T2[0]
        carry1_expected = tt_and(T1[0], T2[0])
        match = np.array_equal(carries[1], carry1_expected)
        print(f"  Verified: {match}")

        # ANF of T1[0] AND T2[0]
        anf_t1_0 = mobius_transform(T1[0])
        anf_t2_0 = mobius_transform(T2[0])
        t1_monos = monomial_set(anf_t1_0)
        t2_monos = monomial_set(anf_t2_0)
        carry1_monos = monomial_set(mobius_transform(carries[1]))

        print(f"  |T1[0]| = {len(t1_monos)}, |T2[0]| = {len(t2_monos)}, "
              f"|carry[1]| = {len(carry1_monos)}")
        print(f"  Degree: T1[0] max={max(degree_of(m) for m in t1_monos) if t1_monos else 0}, "
              f"T2[0] max={max(degree_of(m) for m in t2_monos) if t2_monos else 0}, "
              f"carry[1] max={max(degree_of(m) for m in carry1_monos) if carry1_monos else 0}")

    # ================================================================
    # The equation of understanding
    # ================================================================
    print(f"\n{'='*80}")
    print(f"THE EQUATION OF UNDERSTANDING (Round {R})")
    print(f"{'='*80}")

    # For each bit: what % of output monomials come from T1⊕T2 vs carry?
    print(f"\n  For each output bit: origin of monomials")
    print(f"  {'Bit':>4} {'From T1⊕T2':>12} {'From carry':>12} {'Cancelled':>12} {'Total':>8}")

    for bit in range(N):
        anf_t1 = mobius_transform(T1[bit])
        anf_t2 = mobius_transform(T2[bit])
        anf_carry = mobius_transform(carries[bit])
        anf_result = mobius_transform(a_new[bit])

        t1_xor_t2 = anf_t1 ^ anf_t2  # T1⊕T2 in ANF

        # Monomials in result that are in T1⊕T2 but not carry
        from_xor = monomial_set(t1_xor_t2) & monomial_set(anf_result)
        from_carry = monomial_set(anf_carry) & monomial_set(anf_result)
        # Some are in both (cancelled by XOR of all three)
        from_both = from_xor & from_carry

        total = len(monomial_set(anf_result))
        only_xor = len(from_xor - from_carry)
        only_carry = len(from_carry - from_xor)
        from_all = len(from_both)

        print(f"  {bit:>4} {only_xor:>12} {only_carry:>12} {from_all:>12} {total:>8}")
        print(f"       {only_xor/max(total,1)*100:>11.1f}% {only_carry/max(total,1)*100:>11.1f}% "
              f"{from_all/max(total,1)*100:>11.1f}%")


def main():
    for R in [2, 4, 5, 6]:
        analyze_carry_chain(R)

    # Final summary
    print(f"\n{'='*80}")
    print(f"SUMMARY: The Price of Each Bit of Understanding")
    print(f"{'='*80}")
    print(f"""
  Bit 0: FREE — T1[0] ⊕ T2[0], degree ≤ max(deg(T1), deg(T2))
         No carry. 100% prediction. Degree = degree of state polynomials.

  Bit 1: CHEAP — add carry[1] = T1[0]·T2[0]
         One AND. Degree roughly doubles (product of two polynomials).

  Bit 2: MODERATE — add carry[2] = MAJ(T1[1], T2[1], T1[0]·T2[0])
         Carry depends on previous carry. Degree grows further.

  Bit k: cost = chain of k AND operations
         Degree ≈ 2^k (each AND can double degree)

  For SHA-256 (n=32):
    Bit 0:  degree ≈ d (degree of state)
    Bit 1:  degree ≈ 2d
    Bit 16: degree ≈ 2^16 · d (ASTRONOMICAL)
    Bit 31: degree ≈ 2^31 · d (IMPOSSIBLE to represent)

  THIS is why SHA-256 is hard: the carry chain creates
  EXPONENTIAL degree growth across bit positions.

  But it's STRUCTURED exponential growth — each step is one MAJ.
  The question: can we exploit the structure of the chain
  without computing the full degree-2^31 polynomial?
""")


if __name__ == "__main__":
    main()
