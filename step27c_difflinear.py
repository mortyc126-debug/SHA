#!/usr/bin/env python3
"""
Step 27c: DIFFERENTIAL-LINEAR ATTACK ON SHA-256
════════════════════════════════════════════════

FAST VERSION: Focus on the key experiments that matter.

Core question: Does SHA-256 have ANY linear bias when viewed through
differential lens? Even tiny bias ε means collision < 2^128.
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


def sha256_full(msg_words, rounds=64):
    """SHA-256 compression, returns hash."""
    W = list(msg_words) + [0] * (64 - len(msg_words))
    for i in range(16, 64):
        s0 = rotr(W[i-15], 7) ^ rotr(W[i-15], 18) ^ (W[i-15] >> 3)
        s1 = rotr(W[i-2], 17) ^ rotr(W[i-2], 19) ^ (W[i-2] >> 10)
        W[i] = (W[i-16] + s0 + W[i-7] + s1) & MASK32

    a, b, c, d, e, f, g, h = IV
    for i in range(min(rounds, 64)):
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

    return [(IV[j] + v) & MASK32 for j, v in enumerate([a, b, c, d, e, f, g, h])]


def parity(x):
    x ^= x >> 16; x ^= x >> 8; x ^= x >> 4
    x ^= x >> 2; x ^= x >> 1
    return x & 1


# ═══════════════════════════════════════════════════════════════
# EXPERIMENT 1: Diff-linear bias — the MAIN experiment
# ═══════════════════════════════════════════════════════════════

def difflinear_main(total_rounds, samples=100000):
    """
    For input pair (M, M') with DW[0] = 0x80000000:
    Measure P(DH[bit] = 0) for each bit of hash word 0.

    Random oracle: P = 0.5 exactly.
    Any deviation = distinguisher.
    """
    noise = 1.0 / math.sqrt(4 * samples)

    bit_counts = [0] * 32  # count of DH[bit]=0
    parity_zero = 0
    dh_counter = Counter()

    for _ in range(samples):
        msg = [random.getrandbits(32) for _ in range(16)]
        msg2 = list(msg)
        msg2[0] ^= 0x80000000

        h1 = sha256_full(msg, total_rounds)
        h2 = sha256_full(msg2, total_rounds)

        dh = h1[0] ^ h2[0]
        dh_counter[dh] += 1

        if parity(dh) == 0:
            parity_zero += 1

        for b in range(32):
            if (dh >> b) & 1 == 0:
                bit_counts[b] += 1

    # Analysis
    biases = [abs(bit_counts[b] / samples - 0.5) for b in range(32)]
    max_bias = max(biases)
    max_bit = biases.index(max_bias)
    par_bias = abs(parity_zero / samples - 0.5)
    unique = len(dh_counter)
    max_count = dh_counter.most_common(1)[0][1]

    return {
        'max_bias': max_bias, 'max_bit': max_bit,
        'par_bias': par_bias, 'noise': noise,
        'unique': unique, 'max_count': max_count,
        'biases': biases, 'samples': samples
    }


def run_difflinear_scan():
    """Scan across different round counts."""
    print("DIFFERENTIAL-LINEAR BIAS SCAN")
    print("═" * 70)
    print()
    print("  Input diff: DW[0] = 0x80000000 (MSB of first message word)")
    print("  Measuring: P(DH_bit = 0) for each bit of H[0]")
    print()

    print(f"  {'Rounds':>7}  {'Samples':>8}  {'Max bias':>10}  {'Max bit':>8}  "
          f"{'Par bias':>10}  {'Noise':>8}  {'Max/Noise':>10}  {'Dist?':>6}")
    print(f"  {'-'*7}  {'-'*8}  {'-'*10}  {'-'*8}  "
          f"{'-'*10}  {'-'*8}  {'-'*10}  {'-'*6}")

    for R in [4, 8, 12, 16, 20, 24, 28, 32, 40, 48, 56, 64]:
        N = 100000
        r = difflinear_main(R, samples=N)
        ratio = r['max_bias'] / r['noise']
        dist = "YES" if ratio > 3.5 else "no"

        print(f"  {R:>7}  {N:>8}  {r['max_bias']:>10.6f}  {r['max_bit']:>8}  "
              f"{r['par_bias']:>10.6f}  {r['noise']:>8.6f}  {ratio:>10.1f}  {dist:>6}")

    print()


# ═══════════════════════════════════════════════════════════════
# EXPERIMENT 2: Structure-aligned masks (faster than full scan)
# ═══════════════════════════════════════════════════════════════

def structure_aligned_masks(total_rounds=64, samples=200000):
    """
    Test specific output masks aligned with SHA-256 internal structure.
    These are MORE LIKELY to show bias than random single-bit masks.
    """
    print(f"STRUCTURE-ALIGNED MASKS (R={total_rounds}, N={samples})")
    print("═" * 70)
    print()

    noise = 1.0 / math.sqrt(4 * samples)

    # Pre-compute all pairs
    dh_words = [[] for _ in range(8)]  # DH for each hash word

    for _ in range(samples):
        msg = [random.getrandbits(32) for _ in range(16)]
        msg2 = list(msg)
        msg2[0] ^= 0x80000000

        h1 = sha256_full(msg, total_rounds)
        h2 = sha256_full(msg2, total_rounds)

        for w in range(8):
            dh_words[w].append(h1[w] ^ h2[w])

    # Test masks
    masks = [
        # Single bits across all words
        *[(f"H[{w}][0]", w, lambda dh, w=w: dh & 1) for w in range(8)],
        *[(f"H[{w}][15]", w, lambda dh, w=w: (dh >> 15) & 1) for w in range(8)],
        *[(f"H[{w}][31]", w, lambda dh, w=w: (dh >> 31) & 1) for w in range(8)],
        # Parity of each word
        *[(f"par(H[{w}])", w, lambda dh: parity(dh)) for w in range(8)],
        # Σ1-aligned: bits 6,11,25
        *[(f"H[{w}]Σ1", w, lambda dh: ((dh>>6)^(dh>>11)^(dh>>25)) & 1) for w in range(8)],
        # Σ0-aligned: bits 2,13,22
        *[(f"H[{w}]Σ0", w, lambda dh: ((dh>>2)^(dh>>13)^(dh>>22)) & 1) for w in range(8)],
        # Low byte parity
        *[(f"H[{w}]lo8", w, lambda dh: parity(dh & 0xFF)) for w in range(8)],
        # Cross-word: H[0] ⊕ H[4] (a ⊕ e structure)
        ("H[0]⊕H[4] bit0", -1, None),
        ("par(H[0]⊕H[4])", -1, None),
    ]

    print(f"  {'Mask':<20}  {'Bias':>10}  {'Noise':>10}  {'Ratio':>8}  {'Significant?':>14}")
    print(f"  {'-'*20}  {'-'*10}  {'-'*10}  {'-'*8}  {'-'*14}")

    best_overall = 0
    best_name = ""

    for name, word, mask_fn in masks:
        if word == -1:
            # Cross-word masks
            if "bit0" in name:
                count = sum(1 for i in range(samples)
                           if ((dh_words[0][i] ^ dh_words[4][i]) & 1) == 0)
            else:
                count = sum(1 for i in range(samples)
                           if parity(dh_words[0][i] ^ dh_words[4][i]) == 0)
        else:
            count = sum(1 for dh in dh_words[word] if mask_fn(dh) == 0)

        bias = abs(count / samples - 0.5)
        ratio = bias / noise
        sig = "YES" if ratio > 3.5 else "no"

        if bias > best_overall:
            best_overall = bias
            best_name = name

        if ratio > 2.0 or "par" in name or "[0]" in name[:5]:
            print(f"  {name:<20}  {bias:>10.6f}  {noise:>10.6f}  "
                  f"{ratio:>8.1f}  {sig:>14}")

    print()
    print(f"  Best mask: {best_name} (bias={best_overall:.6f}, "
          f"ratio={best_overall/noise:.1f})")
    print()


# ═══════════════════════════════════════════════════════════════
# EXPERIMENT 3: High-confidence test at 64 rounds
# ═══════════════════════════════════════════════════════════════

def high_confidence_test(samples=500000):
    """
    THE DEFINITIVE TEST: 500K samples for single-bit output biases.
    Noise floor = 1/√(2M) ≈ 0.0007.
    If we find bias > 0.003 → distinguisher at >4σ.
    """
    print(f"HIGH-CONFIDENCE DIFF-LINEAR TEST (N={samples}, 64 rounds)")
    print("═" * 70)
    print()

    noise = 1.0 / math.sqrt(4 * samples)
    print(f"  Noise floor: {noise:.6f}")
    print(f"  4σ threshold: {4*noise:.6f}")
    print()

    bit_counts = [[0]*32 for _ in range(8)]  # 8 hash words × 32 bits
    parity_counts = [0] * 8

    for trial in range(samples):
        msg = [random.getrandbits(32) for _ in range(16)]
        msg2 = list(msg)
        msg2[0] ^= 0x80000000

        h1 = sha256_full(msg, 64)
        h2 = sha256_full(msg2, 64)

        for w in range(8):
            dh = h1[w] ^ h2[w]
            if parity(dh) == 0:
                parity_counts[w] += 1
            for b in range(32):
                if (dh >> b) & 1 == 0:
                    bit_counts[w][b] += 1

        if (trial + 1) % 100000 == 0:
            print(f"    [{trial+1}/{samples}]", flush=True)

    # Find maximum bias across ALL 256 output bits
    all_biases = []
    for w in range(8):
        for b in range(32):
            bias = abs(bit_counts[w][b] / samples - 0.5)
            all_biases.append((bias, w, b))

    all_biases.sort(reverse=True)

    print()
    print(f"  TOP 10 BIASED BITS:")
    print(f"  {'Bit':>10}  {'Bias':>10}  {'Ratio':>8}  {'Significant?':>14}")
    print(f"  {'-'*10}  {'-'*10}  {'-'*8}  {'-'*14}")
    for bias, w, b in all_biases[:10]:
        ratio = bias / noise
        sig = "YES" if ratio > 3.5 else "no"
        print(f"  H[{w}][{b:>2}]   {bias:>10.6f}  {ratio:>8.1f}  {sig:>14}")

    # Parity biases
    print()
    print(f"  PARITY BIASES:")
    for w in range(8):
        par_bias = abs(parity_counts[w] / samples - 0.5)
        ratio = par_bias / noise
        sig = "YES" if ratio > 3.5 else "no"
        print(f"  par(DH[{w}]): {par_bias:.6f} (ratio={ratio:.1f}) {sig}")

    # Multiple testing correction (Bonferroni)
    n_tests = 256 + 8  # 256 bits + 8 parities
    corrected_threshold = noise * math.sqrt(2 * math.log(n_tests))
    print()
    print(f"  Bonferroni-corrected threshold (p=0.05, {n_tests} tests): "
          f"{corrected_threshold:.6f}")

    max_bias = all_biases[0][0]
    print(f"  Maximum observed bias: {max_bias:.6f}")
    print(f"  Verdict: {'DISTINGUISHER FOUND!' if max_bias > corrected_threshold else 'No distinguisher (consistent with random oracle)'}")
    print()

    return all_biases[0]


# ═══════════════════════════════════════════════════════════════
# EXPERIMENT 4: Different input differences
# ═══════════════════════════════════════════════════════════════

def multi_diff_test(samples=100000):
    """Test multiple input differences — maybe MSB isn't optimal."""
    print(f"MULTI-DIFFERENCE TEST (N={samples}, 64 rounds)")
    print("═" * 70)
    print()

    noise = 1.0 / math.sqrt(4 * samples)

    diffs = [
        ("DW[0] = 0x80000000", 0, 0x80000000),
        ("DW[0] = 0x00000001", 0, 0x00000001),
        ("DW[0] = 0xFFFFFFFF", 0, 0xFFFFFFFF),
        ("DW[1] = 0x80000000", 1, 0x80000000),
        ("DW[7] = 0x80000000", 7, 0x80000000),
        ("DW[15]= 0x80000000", 15, 0x80000000),
    ]

    print(f"  {'Difference':<25}  {'Max bit bias':>14}  {'Par bias':>10}  "
          f"{'Ratio':>8}  {'Dist?':>6}")
    print(f"  {'-'*25}  {'-'*14}  {'-'*10}  {'-'*8}  {'-'*6}")

    for name, word, val in diffs:
        bit_counts = [0] * 32
        par_count = 0

        for _ in range(samples):
            msg = [random.getrandbits(32) for _ in range(16)]
            msg2 = list(msg)
            msg2[word] ^= val

            h1 = sha256_full(msg, 64)
            h2 = sha256_full(msg2, 64)

            dh = h1[0] ^ h2[0]
            if parity(dh) == 0:
                par_count += 1
            for b in range(32):
                if (dh >> b) & 1 == 0:
                    bit_counts[b] += 1

        max_bias = max(abs(bit_counts[b]/samples - 0.5) for b in range(32))
        par_bias = abs(par_count/samples - 0.5)
        ratio = max_bias / noise
        dist = "YES" if ratio > 3.5 else "no"

        print(f"  {name:<25}  {max_bias:>14.6f}  {par_bias:>10.6f}  "
              f"{ratio:>8.1f}  {dist:>6}")

    print()


# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════

if __name__ == '__main__':
    print("STEP 27c: DIFFERENTIAL-LINEAR ATTACK ON SHA-256")
    print("Can we find ANY bias → beat 2^128?")
    print("═" * 70)
    print()

    t_start = time.time()

    # 1. Quick scan across rounds
    run_difflinear_scan()

    # 2. Structure-aligned masks at 64 rounds
    structure_aligned_masks(total_rounds=64, samples=100000)

    # 3. Multiple input differences
    multi_diff_test(samples=100000)

    # 4. THE BIG TEST: high confidence at 64 rounds
    high_confidence_test(samples=500000)

    elapsed = time.time() - t_start
    print("═" * 70)
    print(f"TIME: {elapsed:.1f}s")
    print()
    print("DIFFERENTIAL-LINEAR FINAL VERDICT:")
    print("  If max_bias > Bonferroni threshold → SHA-256 ≠ random oracle")
    print("  → collision possible at 2^{128-k} where k = log2(bias/noise)")
    print()
    print("  If no significant bias found → random oracle model holds")
    print("  → 2^128 is TIGHT and no sub-birthday collision exists")
