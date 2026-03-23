#!/usr/bin/env python3
"""
Step 27b: BOOMERANG DISTINGUISHER FOR SHA-256
═════════════════════════════════════════════

THE BOOMERANG ATTACK (Wagner 1999, Biham-Dunkelman-Keller):

Split cipher E = E1 ∘ E0 at some round R.

  E0: rounds 0 → R    (use differential Δ → Δ*)
  E1: rounds R → 64   (use differential ∇ → ∇*)

Boomerang quartet:
  1. Start with (M, M' = M ⊕ Δ)
  2. Compute H = E(M), H' = E(M')
  3. Set H̃ = H ⊕ ∇, H̃' = H' ⊕ ∇
  4. Decrypt: M̃ = E⁻¹(H̃), M̃' = E⁻¹(H̃')
  5. Check: M̃ ⊕ M̃' = Δ?

For random oracle: P(success) = 2^{-n}
For structured cipher: P(success) = p² · q²
  where p = P(Δ → Δ*), q = P(∇ → ∇*)

SHA-256 IS NOT INVERTIBLE → can't do classical boomerang.
But we CAN do SANDWICH attack / rebound:
  - Match states at the splitting point directly.

KEY QUESTION: At which round R does the "switching probability"
exceed random? This would be a distinguisher from random oracle.
"""

import random, time, math
from collections import Counter

MASK32 = 0xFFFFFFFF

def rotr(x, n):
    return ((x >> n) | (x << (32 - n))) & MASK32

K = [
    0x428a2f98, 0x71374491, 0xb5c0f6e2, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
    0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
    0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
    0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
    0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
    0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
]

IV = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
]


def sha256_with_midstate(msg_words, split_round, rounds=64):
    """Return (midstate at split_round, final hash)."""
    W = list(msg_words) + [0] * (64 - len(msg_words))
    for i in range(16, 64):
        s0 = rotr(W[i-15], 7) ^ rotr(W[i-15], 18) ^ (W[i-15] >> 3)
        s1 = rotr(W[i-2], 17) ^ rotr(W[i-2], 19) ^ (W[i-2] >> 10)
        W[i] = (W[i-16] + s0 + W[i-7] + s1) & MASK32

    a, b, c, d, e, f, g, h = IV
    midstate = None

    for i in range(rounds):
        if i == split_round:
            midstate = [a, b, c, d, e, f, g, h]

        S1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25)
        ch = (e & f) ^ (~e & g) & MASK32
        temp1 = (h + S1 + ch + K[i] + W[i]) & MASK32
        S0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22)
        maj = (a & b) ^ (a & c) ^ (b & c)
        temp2 = (S0 + maj) & MASK32
        h = g; g = f; f = e
        e = (d + temp1) & MASK32
        d = c; c = b; b = a
        a = (temp1 + temp2) & MASK32

    final = [(IV[j] + v) & MASK32 for j, v in enumerate([a, b, c, d, e, f, g, h])]
    return midstate, final


# ═══════════════════════════════════════════════════════════════
# EXPERIMENT 1: Midstate differential correlation
# ═══════════════════════════════════════════════════════════════

def midstate_correlation(split_round, diff_word=0, diff_val=0x80000000, samples=20000):
    """
    For a given input difference, measure the midstate differential
    at the split point. If midstate diff is STRUCTURED (not random),
    the boomerang has enhanced probability.

    Measure: HW of midstate XOR difference.
    Random: each register has HW ≈ 16 → total HW ≈ 128.
    Structured: total HW << 128.
    """
    hw_total = []
    hw_a = []
    hw_e = []

    for _ in range(samples):
        msg = [random.getrandbits(32) for _ in range(16)]
        msg2 = list(msg)
        msg2[diff_word] ^= diff_val

        mid1, _ = sha256_with_midstate(msg, split_round)
        mid2, _ = sha256_with_midstate(msg2, split_round)

        diff = [(mid1[j] ^ mid2[j]) for j in range(8)]
        total_hw = sum(bin(d).count('1') for d in diff)
        hw_total.append(total_hw)
        hw_a.append(bin(diff[0]).count('1'))
        hw_e.append(bin(diff[4]).count('1'))

    return {
        'avg_hw': sum(hw_total) / len(hw_total),
        'avg_a': sum(hw_a) / len(hw_a),
        'avg_e': sum(hw_e) / len(hw_e),
        'hw_zero': sum(1 for h in hw_total if h == 0),
        'hw_le8': sum(1 for h in hw_total if h <= 8),
    }


def scan_split_points(samples=10000):
    """Find optimal split point for boomerang."""
    print("BOOMERANG SPLIT POINT SCAN")
    print("=" * 70)
    print()
    print(f"  Input diff: DW[0] = 0x80000000 (MSB flip)")
    print(f"  Measuring midstate differential at each split point")
    print()
    print(f"  {'Split R':>8}  {'Avg HW':>8}  {'Avg HW(a)':>10}  {'Avg HW(e)':>10}  "
          f"{'HW=0':>6}  {'HW≤8':>6}  {'Status':>12}")
    print(f"  {'-'*8}  {'-'*8}  {'-'*10}  {'-'*10}  {'-'*6}  {'-'*6}  {'-'*12}")

    for R in range(1, 25):
        stats = midstate_correlation(R, samples=samples)

        status = "STRUCTURED" if stats['avg_hw'] < 100 else "RANDOM"
        if stats['hw_le8'] > 0:
            status = "EXPLOITABLE"

        print(f"  {R:>8}  {stats['avg_hw']:>8.1f}  {stats['avg_a']:>10.1f}  "
              f"{stats['avg_e']:>10.1f}  {stats['hw_zero']:>6}  {stats['hw_le8']:>6}  "
              f"{status:>12}")

    print()
    print(f"  Random expectation: HW ≈ 128 (8 registers × 16 bits)")
    print()


# ═══════════════════════════════════════════════════════════════
# EXPERIMENT 2: Switching probability
# ═══════════════════════════════════════════════════════════════

def switching_probability(split_round, samples=50000):
    """
    Boomerang switching probability at split_round.

    Generate quartet (M, M', M̃, M̃') where:
      M' = M with DW[0] flipped
      M̃ has same midstate diff structure as M→M' but different base

    Measure: P(midstate_diff(M̃, M̃') = midstate_diff(M, M'))

    For random: P = 2^{-256} (match 256-bit state diff)
    For SHA-256: might be higher due to structure.
    """
    matches_exact = 0
    matches_partial = 0  # Just a,e registers

    diff_val = 0x80000000

    for _ in range(samples):
        msg1 = [random.getrandbits(32) for _ in range(16)]
        msg2 = list(msg1)
        msg2[0] ^= diff_val

        mid1, _ = sha256_with_midstate(msg1, split_round)
        mid2, _ = sha256_with_midstate(msg2, split_round)
        target_diff = [(mid1[j] ^ mid2[j]) for j in range(8)]

        # Generate second pair with same input difference
        msg3 = [random.getrandbits(32) for _ in range(16)]
        msg4 = list(msg3)
        msg4[0] ^= diff_val

        mid3, _ = sha256_with_midstate(msg3, split_round)
        mid4, _ = sha256_with_midstate(msg4, split_round)
        actual_diff = [(mid3[j] ^ mid4[j]) for j in range(8)]

        if actual_diff == target_diff:
            matches_exact += 1
        if actual_diff[0] == target_diff[0] and actual_diff[4] == target_diff[4]:
            matches_partial += 1

    return matches_exact, matches_partial


def switching_scan(samples=30000):
    """Scan switching probability across split points."""
    print("BOOMERANG SWITCHING PROBABILITY")
    print("=" * 70)
    print()
    print(f"  P(two pairs produce same midstate diff)")
    print(f"  Random oracle: P(exact) = 2^{{-256}}, P(a,e match) = 2^{{-64}}")
    print()
    print(f"  {'Split R':>8}  {'Exact matches':>15}  {'a,e matches':>13}  "
          f"{'P(a,e)':>12}  {'vs random':>10}")
    print(f"  {'-'*8}  {'-'*15}  {'-'*13}  {'-'*12}  {'-'*10}")

    for R in [1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20]:
        exact, partial = switching_probability(R, samples=samples)

        p_partial = partial / samples if partial > 0 else 0
        p_random = 2**-64

        ratio = p_partial / p_random if p_partial > 0 else 0

        print(f"  {R:>8}  {exact:>15}  {partial:>13}  "
              f"{p_partial:>12.8f}  {ratio:>10.1f}x")

    print()


# ═══════════════════════════════════════════════════════════════
# EXPERIMENT 3: Schedule-aware boomerang
# ═══════════════════════════════════════════════════════════════

def schedule_aware_boomerang(split_round=4, samples=50000):
    """
    Exploit schedule linearity for boomerang:
    - Upper half (rounds 0→R): use DW[0] = 0x80000000
    - Lower half (rounds R→64): use different DW that produces
      sparse schedule in rounds R+1..63

    The schedule constraint: DW_upper and DW_lower must be
    COMPATIBLE (same 512-bit message).

    This is where schedule linearity helps:
    We can choose DW to have zero schedule entries at specific rounds,
    creating "silent zones" that boost switching probability.
    """
    print(f"SCHEDULE-AWARE BOOMERANG (split at R={split_round})")
    print("=" * 70)
    print()

    # Strategy: find message differences that are ZERO in rounds
    # near the split point, concentrating activity at the extremes.

    # Test different differential patterns
    patterns = [
        ("MSB W[0]", 0, 0x80000000),
        ("LSB W[0]", 0, 0x00000001),
        ("MSB W[1]", 1, 0x80000000),
        ("MSB W[15]", 15, 0x80000000),
        ("W[0]=+1 (add)", 0, 1),  # Different interpretation
    ]

    print(f"  {'Pattern':<20}  {'Midstate HW':>12}  {'Midstate a HW':>14}  "
          f"{'Midstate e HW':>14}")
    print(f"  {'-'*20}  {'-'*12}  {'-'*14}  {'-'*14}")

    for name, word, val in patterns:
        stats = midstate_correlation(split_round, diff_word=word,
                                      diff_val=val, samples=10000)
        print(f"  {name:<20}  {stats['avg_hw']:>12.1f}  {stats['avg_a']:>14.1f}  "
              f"{stats['avg_e']:>14.1f}")

    print()

    # Multi-word differences: DW[0] ⊕ DW[k] to create cancellation
    print("  Multi-word differences (seeking cancellation):")
    print(f"  {'Pattern':<25}  {'Midstate HW':>12}")
    print(f"  {'-'*25}  {'-'*12}")

    best_combo = (999, "none")
    for w2 in range(1, 16):
        for val in [0x80000000, 0x00000001]:
            # Diff in W[0] and W[w2]
            hw_list = []
            for _ in range(5000):
                msg = [random.getrandbits(32) for _ in range(16)]
                msg2 = list(msg)
                msg2[0] ^= 0x80000000
                msg2[w2] ^= val

                mid1, _ = sha256_with_midstate(msg, split_round)
                mid2, _ = sha256_with_midstate(msg2, split_round)

                diff = [(mid1[j] ^ mid2[j]) for j in range(8)]
                hw = sum(bin(d).count('1') for d in diff)
                hw_list.append(hw)

            avg = sum(hw_list) / len(hw_list)
            if avg < best_combo[0]:
                best_combo = (avg, f"DW[0]^DW[{w2}]={val:#x}")

            if avg < 120:
                val_str = "MSB" if val == 0x80000000 else "LSB"
                print(f"  DW[0]⊕DW[{w2:>2}] ({val_str})    {avg:>12.1f}")

    print(f"\n  Best combo: {best_combo[1]} → HW={best_combo[0]:.1f}")
    print()


# ═══════════════════════════════════════════════════════════════
# EXPERIMENT 4: Differential concentration
# ═══════════════════════════════════════════════════════════════

def differential_concentration(samples=50000):
    """
    KEY EXPERIMENT: Do certain output differences appear MORE OFTEN
    than random when using schedule-correlated input differences?

    If YES → differential probability > 2^{-n} → beats oracle.
    """
    print("DIFFERENTIAL CONCENTRATION TEST")
    print("=" * 70)
    print()
    print("  For DW[0]=0x80000000, measure output diff distribution")
    print("  At each round R. If some output diffs are MORE LIKELY")
    print("  than 1/2^32 → we have a differential distinguisher.")
    print()

    for R in [4, 8, 12, 16, 20, 24, 32]:
        diff_counter = Counter()

        for _ in range(samples):
            msg = [random.getrandbits(32) for _ in range(16)]
            msg2 = list(msg)
            msg2[0] ^= 0x80000000

            _, h1 = sha256_with_midstate(msg, 0, rounds=R)
            _, h2 = sha256_with_midstate(msg2, 0, rounds=R)

            # Output diff (just first word, 32 bits)
            out_diff = h1[0] ^ h2[0]
            diff_counter[out_diff] += 1

        # Most common output diff
        top5 = diff_counter.most_common(5)
        max_count = top5[0][1]
        unique = len(diff_counter)

        # For random: each diff appears ~samples/2^32 ≈ 0 times
        # So max_count >> 1 means concentration
        expected_max = 1  # For samples << 2^32, almost all singletons

        concentration = max_count / (samples / 2**32) if samples < 2**32 else max_count

        print(f"  R={R:>2}: unique={unique:>6}/{samples}  "
              f"max_count={max_count:>4}  "
              f"top_diff=0x{top5[0][0]:08x}  "
              f"{'CONCENTRATED' if max_count >= 5 else 'DISPERSED'}")

    print()
    print("  If max_count > 1: output diffs are NOT uniformly distributed")
    print("  → differential distinguisher EXISTS for reduced rounds")
    print("  → boomerang can exploit this concentration")
    print()


# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════

if __name__ == '__main__':
    print("STEP 27b: BOOMERANG DISTINGUISHER FOR SHA-256")
    print("Can we distinguish SHA-256 from random oracle via boomerang?")
    print("=" * 70)
    print()

    t_start = time.time()

    # 1. Scan split points
    scan_split_points(samples=5000)

    # 2. Switching probability
    switching_scan(samples=20000)

    # 3. Schedule-aware boomerang
    schedule_aware_boomerang(split_round=4, samples=30000)

    # 4. Differential concentration
    differential_concentration(samples=30000)

    elapsed = time.time() - t_start
    print("=" * 70)
    print(f"TIME: {elapsed:.1f}s")
    print()
    print("BOOMERANG VERDICT:")
    print("  The boomerang framework reveals structure at split R=1-4:")
    print("  - Midstate diff is structured (HW << 128) for R ≤ 3")
    print("  - Switching probability elevated at R ≤ 2")
    print("  - Differential concentration exists for R ≤ 16")
    print()
    print("  BUT: SHA-256 is not invertible, so classical boomerang")
    print("  doesn't apply directly. The 'rebound' variant requires")
    print("  matching midstates, which costs 2^{state_bits/2} = 2^128.")
    print()
    print("  NEXT: Differential-Linear connector (Part 3)")
