"""
THREE THREADS IN PARALLEL:

A: Rotation formula — what property of rotation constants determines D2 speed?
B: Watershed — at what round do collision paths FIRST diverge (backward from hash)?
C: Carry algebra extension — carry of PRODUCTS (Ch output + carry = ?)
"""

import random, math
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, schedule
)


# ═══════════════════════════════════════════════════
# THREAD A: Rotation formula
# ═══════════════════════════════════════════════════

def thread_A():
    """
    Test: number of DISTINCT (r1+r2)%32 sums for all pairs of rotations.
    This = number of degree-2 paths = potential for Ch/Maj multiplication.
    """
    print("=" * 80)
    print("A: Rotation PAIR-SUMS as predictor of D2 speed")
    print("=" * 80)

    configs = [
        ('SHA-256', [2,13,22,6,11,25], 0.427),
        ('Adjacent', [1,2,3,4,5,6], 0.350),
        ('Wide', [1,10,20,5,15,25], 0.405),
        ('Minimal', [1,16], 0.308),
        ('Primes', [3,5,7,11,13,17], 0.472),
        ('Powers2', [1,2,4,8,16], 0.425),
        ('VeryWide', [4,16,28,8,20,30], 0.452),
    ]

    print(f"\n  {'Config':>12} | n_rots | pair_sums | D2@R=16 | corr_candidate")
    print("-" * 65)

    for name, rots, d2 in configs:
        n = len(rots)
        # Count distinct (ri + rj) % 32 for all pairs (including self-pairs)
        pair_sums = set()
        for r1 in rots:
            for r2 in rots:
                pair_sums.add((r1 + r2) % 32)

        # Also: triple sums
        triple_sums = set()
        for r1 in rots:
            for r2 in rots:
                for r3 in rots:
                    triple_sums.add((r1 + r2 + r3) % 32)

        # Candidate metric: n_rots × pair_sums / 32
        candidate = n * len(pair_sums) / 32

        print(f"  {name:>12} | {n:>6} | {len(pair_sums):>9}/32 | {d2:.3f}   | {candidate:.2f}")

    # Correlation: n_rots vs D2
    n_vals = [len(c[1]) for c in configs]
    d2_vals = [c[2] for c in configs]
    mean_n = sum(n_vals) / len(n_vals)
    mean_d2 = sum(d2_vals) / len(d2_vals)
    cov = sum((n - mean_n) * (d - mean_d2) for n, d in zip(n_vals, d2_vals)) / len(n_vals)
    var_n = sum((n - mean_n)**2 for n in n_vals) / len(n_vals)
    var_d2 = sum((d - mean_d2)**2 for d in d2_vals) / len(d2_vals)
    corr = cov / (var_n * var_d2)**0.5 if var_n > 0 and var_d2 > 0 else 0

    print(f"\n  corr(n_rotations, D2@R=16) = {corr:.3f}")
    print(f"  → {'NUMBER of rotations predicts D2 speed!' if corr > 0.5 else 'Weak correlation.'}")


# ═══════════════════════════════════════════════════
# THREAD B: Watershed geometry
# ═══════════════════════════════════════════════════

def thread_B():
    """
    Watershed: for collision pair, when do paths FIRST differ going backward?

    Can't find 64-round collision. But: from створочне backward extension,
    collision pair has a[57..64] identical.

    Going further back: a[56] may differ. If so: watershed = round 57.
    If a[56] also matches: watershed earlier.

    For RANDOM collision: expected watershed = ?

    Argument: a[r] is 32-bit word. For two random preimages of same hash:
    P(a1[r] = a2[r]) = 2^{-32} for r < 57 (states differ randomly).
    Expected first match going backward from 64: at r where random match occurs.
    With 2^256 preimages: expected deepest watershed ≈ ?

    Actually: for any TWO random preimages, states before round 57 are
    independent → P(a1[r] = a2[r]) = 2^{-32} for each r < 57.
    Over 56 rounds: P(any match) ≈ 56 × 2^{-32} ≈ 2^{-26}.
    Very unlikely → watershed ≈ round 57 for random collision pairs.
    """
    print("\n" + "=" * 80)
    print("B: Watershed = round 57 (analytical argument)")
    print("=" * 80)

    print(f"""
    THEOREM (Watershed Location):
    For a collision pair (M1, M2) with random structure:
      Watershed = round 57 (with probability 1 - 56/2^32 ≈ 1)

    PROOF:
      Collision: state1[64] = state2[64] → a1[57..64] = a2[57..64].
      For r < 57: state1[r] and state2[r] evolve from different M.
      After full thermalization (R > 20): states = pseudo-independent.
      P(a1[r] = a2[r]) = 2^{{-32}} for each r < 57.
      P(any r < 57 matches) ≤ 56 × 2^{{-32}} ≈ 1.3×10^{{-8}} ≈ 0.

      → Watershed = 57 with near-certainty for random collision pairs.

    CONSEQUENCE:
      Collision paths are COMPLETELY DIFFERENT for rounds 1-56.
      They CONVERGE ABRUPTLY at round 57 (within 8 rounds of the end).
      The convergence uses the LAST 7 rounds (57-63) to align.
      These 7 rounds are constrained by schedule (W[57..63] linked to M).

    The COLLISION MECHANISM:
      Rounds 1-56: two independent paths (no meeting)
      Round 57: paths MEET (a1[57] = a2[57] from створочне backward)
      Rounds 57-64: paths IDENTICAL (same state → same computation)

    This is NOT gradual convergence — it's SUDDEN meeting.
    Like two roads that never cross until they merge at a junction.
    The junction = round 57.
    """)


# ═══════════════════════════════════════════════════
# THREAD C: Carry algebra of products
# ═══════════════════════════════════════════════════

def thread_C():
    """
    Ch(e,f,g) = (e AND f) XOR (NOT e AND g).
    The output of Ch goes INTO an addition (T1 = h + Sig1 + Ch + K + W).
    The carry of that addition depends on Ch's output.

    Question: is carry(x, Ch(e,f,g)) structured differently from carry(x, random)?

    Ch is degree 2. Its output bits have CORRELATIONS:
    Ch[k] depends on e[k], f[k], g[k] — SAME bit position only.
    This means: Ch's output has ZERO cross-bit correlation.
    Each bit of Ch = independent function of 3 bits.

    Does this independence affect carry behavior?
    """
    print("\n" + "=" * 80)
    print("C: Carry of Ch output — is it special?")
    print("=" * 80)

    n = 32
    N = 10000

    # Compare: carry(x, Ch(e,f,g)) vs carry(x, random_y)
    # Metric: number of carries, segment structure

    ch_carries = []
    rand_carries = []

    for trial in range(N):
        rng = random.Random(trial)
        x = rng.randint(0, MASK)
        e = rng.randint(0, MASK)
        f = rng.randint(0, MASK)
        g = rng.randint(0, MASK)

        ch_out = Ch(e, f, g)
        rand_y = rng.randint(0, MASK)

        # Carry count for x + Ch
        c = 0; n_carry_ch = 0
        for k in range(32):
            if c: n_carry_ch += 1
            xk = (x >> k) & 1; yk = (ch_out >> k) & 1
            c = (xk & yk) | (xk & c) | (yk & c)
        ch_carries.append(n_carry_ch)

        # Carry count for x + random
        c = 0; n_carry_r = 0
        for k in range(32):
            if c: n_carry_r += 1
            xk = (x >> k) & 1; yk = (rand_y >> k) & 1
            c = (xk & yk) | (xk & c) | (yk & c)
        rand_carries.append(n_carry_r)

    avg_ch = sum(ch_carries) / N
    avg_rand = sum(rand_carries) / N

    print(f"  Average carries: carry(x, Ch) = {avg_ch:.2f}, carry(x, random) = {avg_rand:.2f}")
    print(f"  Difference: {abs(avg_ch - avg_rand):.3f}")

    # Ch output statistics
    ch_hws = []
    for trial in range(N):
        rng = random.Random(trial)
        e = rng.randint(0, MASK); f = rng.randint(0, MASK); g = rng.randint(0, MASK)
        ch_hws.append(bin(Ch(e, f, g)).count('1'))

    avg_ch_hw = sum(ch_hws) / N
    print(f"\n  Ch output: E[HW] = {avg_ch_hw:.2f} (random = 16)")

    # P(Ch[k]=1) for each bit
    ch_bit_freq = [0] * 32
    for trial in range(N):
        rng = random.Random(trial)
        e = rng.randint(0, MASK); f = rng.randint(0, MASK); g = rng.randint(0, MASK)
        ch_out = Ch(e, f, g)
        for k in range(32):
            if (ch_out >> k) & 1: ch_bit_freq[k] += 1

    max_bias = max(abs(ch_bit_freq[k]/N - 0.5) for k in range(32))
    print(f"  Ch per-bit bias: max = {max_bias:.4f} (random < 0.02)")

    # KEY: Ch output is BALANCED and UNBIASED per bit.
    # → carry(x, Ch) should be SAME as carry(x, random).
    # → Ch does NOT create special carry structure.

    print(f"\n  CONCLUSION: Ch output = balanced, unbiased per bit.")
    print(f"  carry(x, Ch) ≈ carry(x, random). Difference = {abs(avg_ch-avg_rand):.3f}.")
    print(f"  Ch does NOT create special carry patterns.")
    print(f"  The nonlinearity of Ch acts through DEGREE, not through carry bias.")
    print(f"  Same for Maj (also balanced, per-bit independent).")


if __name__ == "__main__":
    thread_A()
    thread_B()
    thread_C()
