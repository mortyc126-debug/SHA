"""
Step 12: Cocycle-Based Algebra — formulate collision through T5

The carry cocycle (BTE T5):
  E(a,b) = (a+b) ⊕ (a⊕b)     [carry correction of one addition]
  E(a,b,c) = E(a,b) ⊕ E(a+b,c)  [composition law]

SHA round: a_new = T1 + T2 = T1 ⊕ T2 ⊕ E(T1, T2)

Collision at bit k: a_new(M)[k] = a_new(M')[k]
  where M' = M ⊕ δ

This means:
  (T1⊕T2⊕E(T1,T2))(M)[k] = (T1⊕T2⊕E(T1,T2))(M')[k]

Let Δf = f(M) ⊕ f(M'). Then:
  ΔT1[k] ⊕ ΔT2[k] ⊕ ΔE(T1,T2)[k] = 0

The XOR parts (ΔT1, ΔT2) are "easy" — they propagate linearly
through XOR and rotation, with degree-2 from Ch/Maj.

The carry part ΔE is the HARD part. Let's study it specifically.

KEY IDEA: Express ΔE as a function of ΔT1, ΔT2, and the BASE values
T1(M), T2(M). If ΔE has special structure in this formulation,
we might solve the collision without degree explosion.
"""

import numpy as np
from step0_exact_algebra import mini_sha, N, MASK

N_MSG = 4
N_INPUT = N * N_MSG
N_TOTAL = 1 << N_INPUT


def carry_correction(a, b):
    """E(a,b) = (a+b) ⊕ (a⊕b) — the carry correction."""
    return ((a + b) & MASK) ^ (a ^ b)


def carry_correction_bits(a, b, n=N):
    """Return carry correction as list of bits."""
    e = carry_correction(a, b)
    return [(e >> k) & 1 for k in range(n)]


def delta_carry(T1_base, T2_base, T1_mod, T2_mod):
    """
    ΔE = E(T1(M), T2(M)) ⊕ E(T1(M'), T2(M'))

    This is the carry DIFFERENCE between base and modified computations.
    """
    E_base = carry_correction(T1_base, T2_base)
    E_mod = carry_correction(T1_mod, T2_mod)
    return E_base ^ E_mod


def analyze_delta_carry():
    """
    Study ΔE(T1,T2) as a function of (ΔT1, ΔT2) given fixed (T1_base, T2_base).

    If ΔE depends ONLY on (ΔT1, ΔT2) and NOT on the base values
    → we can solve the carry equation independently.

    If ΔE depends on base values → state-dependent (harder).
    """
    import random
    rng = random.Random(42)

    print(f"{'='*80}")
    print(f"DELTA-CARRY ANALYSIS")
    print(f"{'='*80}")

    # Test: for FIXED ΔT1, ΔT2, does ΔE depend on the base T1, T2?
    print(f"\n  Test: is ΔE(T1,T2) determined by (ΔT1, ΔT2) alone?")
    print(f"  Fix ΔT1=1, ΔT2=0. Vary T1_base over all values.\n")

    for dT1 in [1, 2, 5]:
        for dT2 in [0, 1, 3]:
            values = set()
            for T1_base in range(1 << N):
                for T2_base in range(1 << N):
                    T1_mod = T1_base ^ dT1
                    T2_mod = T2_base ^ dT2
                    de = delta_carry(T1_base, T2_base, T1_mod, T2_mod)
                    values.add(de)

            print(f"  ΔT1={dT1:>2}, ΔT2={dT2:>2}: |ΔE values| = {len(values):>3} "
                  f"/ {1 << N} possible  "
                  f"{'→ STATE-DEPENDENT' if len(values) > 1 else '→ DETERMINED!'}")

    # ΔE is state-dependent. But HOW state-dependent?
    print(f"\n  Characterizing state-dependence:")
    print(f"  For fixed (ΔT1, ΔT2), how many distinct ΔE values?")
    print(f"  And what's the DISTRIBUTION?")

    for dT1 in range(1 << N):
        for dT2 in range(1 << N):
            if dT1 == 0 and dT2 == 0:
                continue
            de_count = {}
            for T1_base in range(1 << N):
                for T2_base in range(1 << N):
                    T1_mod = T1_base ^ dT1
                    T2_mod = T2_base ^ dT2
                    de = delta_carry(T1_base, T2_base, T1_mod, T2_mod)
                    de_count[de] = de_count.get(de, 0) + 1
            # Only print interesting cases
            n_values = len(de_count)
            if n_values <= 4 or (dT1 < 4 and dT2 < 4):
                dist = sorted(de_count.items())
                dist_str = ", ".join(f"{v}:{c}" for v, c in dist[:8])
                print(f"  ΔT1={dT1:>2x}, ΔT2={dT2:>2x}: {n_values:>2} values — {dist_str}")

    # ================================================================
    # KEY: Carry difference as function of BASE state, given δ
    # ================================================================
    print(f"\n{'='*80}")
    print(f"COCYCLE DECOMPOSITION OF COLLISION")
    print(f"{'='*80}")

    # For a specific round: collision requires
    #   ΔT1[k] ⊕ ΔT2[k] ⊕ ΔE[k] = 0 for all k
    #
    # ΔE[k] depends on T1, T2 (base state) AND ΔT1, ΔT2.
    #
    # Rewrite: ΔE[k] = ΔT1[k] ⊕ ΔT2[k]
    # This means: the carry difference must EXACTLY CANCEL the XOR difference.
    #
    # E(T1,T2)[k] ⊕ E(T1⊕ΔT1, T2⊕ΔT2)[k] = ΔT1[k] ⊕ ΔT2[k]

    print(f"\n  Collision condition per bit:")
    print(f"  ΔE(T1,T2)[k] = ΔT1[k] ⊕ ΔT2[k]")
    print(f"  = carry_diff must cancel XOR_diff exactly")

    # For bit 0: ΔE[0] = 0 always (no carry into bit 0)
    # So collision at bit 0: ΔT1[0] ⊕ ΔT2[0] = 0 ← pure XOR condition
    print(f"\n  Bit 0: ΔE[0] = 0 → condition: ΔT1[0] = ΔT2[0] (XOR only)")

    # For bit 1: ΔE[1] = Δcarry[1] = Δ(T1[0]·T2[0])
    #   = (T1[0]⊕ΔT1[0])·(T2[0]⊕ΔT2[0]) ⊕ T1[0]·T2[0]
    #   = T1[0]·ΔT2[0] ⊕ ΔT1[0]·T2[0] ⊕ ΔT1[0]·ΔT2[0]
    print(f"  Bit 1: ΔE[1] = T1[0]·ΔT2[0] ⊕ ΔT1[0]·T2[0] ⊕ ΔT1[0]·ΔT2[0]")
    print(f"         condition: ΔT1[1] ⊕ ΔT2[1] = T1[0]·ΔT2[0] ⊕ ΔT1[0]·T2[0] ⊕ ΔT1[0]·ΔT2[0]")
    print(f"         This is DEGREE 1 in (ΔT1, ΔT2) when base (T1, T2) is fixed!")
    print(f"         (the ΔT1·ΔT2 term is degree 2, but T1·ΔT2 is degree 1)")

    # For bit k: ΔE[k] involves carry chain, but carry chain is RECURSIVE
    # ΔE[k] = function of (T1[0..k-1], T2[0..k-1], ΔT1[0..k-1], ΔT2[0..k-1])

    # VERIFY the formula for bit 1
    print(f"\n  VERIFICATION of bit-1 formula:")
    errors = 0
    total = 0
    for T1 in range(1 << N):
        for T2 in range(1 << N):
            for dT1 in range(1 << N):
                for dT2 in range(1 << N):
                    T1m = T1 ^ dT1
                    T2m = T2 ^ dT2

                    # Actual ΔE[1]
                    E1 = carry_correction(T1, T2)
                    E2 = carry_correction(T1m, T2m)
                    actual_DE1 = ((E1 >> 1) & 1) ^ ((E2 >> 1) & 1)

                    # Formula: T1[0]·ΔT2[0] ⊕ ΔT1[0]·T2[0] ⊕ ΔT1[0]·ΔT2[0]
                    t1_0 = T1 & 1
                    t2_0 = T2 & 1
                    dt1_0 = dT1 & 1
                    dt2_0 = dT2 & 1
                    formula_DE1 = (t1_0 & dt2_0) ^ (dt1_0 & t2_0) ^ (dt1_0 & dt2_0)

                    if actual_DE1 != formula_DE1:
                        errors += 1
                    total += 1

    print(f"  Errors: {errors}/{total} {'✓ VERIFIED' if errors == 0 else '✗ FAILED'}")

    # ================================================================
    # The STRUCTURE of ΔE per bit
    # ================================================================
    print(f"\n{'='*80}")
    print(f"STRUCTURE OF ΔE[k] — degree in (ΔT1, ΔT2) for fixed base")
    print(f"{'='*80}")

    # For each bit k, ΔE[k] is a function of (T1[0..k-1], T2[0..k-1], ΔT1[0..k-1], ΔT2[0..k-1])
    # When base (T1, T2) is FIXED, ΔE[k] is a polynomial in (ΔT1[0..k-1], ΔT2[0..k-1])
    # What degree?

    # Measure empirically: for fixed T1_base, T2_base, what's the degree of ΔE[k]
    # as a function of (ΔT1, ΔT2)?

    T1_base = 0b0110  # fixed
    T2_base = 0b1011

    print(f"\n  Base: T1={T1_base:04b}, T2={T2_base:04b}")

    for k in range(N):
        # Truth table of ΔE[k] as function of (ΔT1, ΔT2) — 2^(2N) = 256 entries
        from step0_exact_algebra import mobius_transform
        tt = np.zeros(1 << (2*N), dtype=np.uint8)

        for idx in range(1 << (2*N)):
            dT1 = idx & MASK
            dT2 = (idx >> N) & MASK
            T1m = T1_base ^ dT1
            T2m = T2_base ^ dT2
            E1 = carry_correction(T1_base, T2_base)
            E2 = carry_correction(T1m, T2m)
            de_k = ((E1 >> k) & 1) ^ ((E2 >> k) & 1)
            tt[idx] = de_k

        anf = mobius_transform(tt, 2*N)
        # Max degree
        max_deg = 0
        n_monos = 0
        for m in range(1 << (2*N)):
            if anf[m]:
                deg = bin(m).count('1')
                max_deg = max(max_deg, deg)
                n_monos += 1

        print(f"  ΔE[{k}]: degree={max_deg}, monomials={n_monos}/{1<<(2*N)} "
              f"({n_monos/(1<<(2*N))*100:.1f}%)")

    # Compare: degree of full a_new collision vs degree of ΔE alone
    print(f"\n  COMPARISON:")
    print(f"  Full collision polynomial (Step 7): degree 9-16 in 16 message bits")
    print(f"  ΔE[k] for fixed base: degree in (ΔT1, ΔT2) — see above")
    print(f"  If ΔE[k] has LOW degree → carry equation solvable!")
    print(f"  If ΔE[k] has HIGH degree → same barrier as before")


def main():
    analyze_delta_carry()


if __name__ == "__main__":
    main()
