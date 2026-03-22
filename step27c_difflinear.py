#!/usr/bin/env python3
"""
Step 27c: DIFFERENTIAL-LINEAR ATTACK ON SHA-256
════════════════════════════════════════════════

THE IDEA (Langford-Hellman 1994, Biham-Dunkelman-Keller 2002):

Split E = E1 ∘ E0:
  E0 (rounds 0→R):  Differential phase — high-prob differential Δ→Δ*
  E1 (rounds R→64): Linear phase — linear approximation with bias ε

Combined: P(collision on linear mask) = p · (1/2 + ε)
  where p = P(Δ→Δ*), ε = linear bias of E1

For SHA-256:
  - Wang chain gives p ≈ 1 for first 20 rounds
  - If ANY linear bias ε > 0 exists in rounds 21-64...
  - ...then collision complexity < 2^128

THIS IS THE CRITICAL EXPERIMENT: Find linear bias in SHA-256 tail rounds.

Even ε = 2^{-10} would give advantage: 2^{128-10} = 2^{118}.
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


def sha256_rounds(state, W, start_round, end_round):
    """Run SHA-256 round function from start_round to end_round."""
    a, b, c, d, e, f, g, h = state
    for i in range(start_round, end_round):
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
    return [a, b, c, d, e, f, g, h]


def expand_schedule(msg_words):
    W = list(msg_words) + [0] * (64 - len(msg_words))
    for i in range(16, 64):
        s0 = rotr(W[i-15], 7) ^ rotr(W[i-15], 18) ^ (W[i-15] >> 3)
        s1 = rotr(W[i-2], 17) ^ rotr(W[i-2], 19) ^ (W[i-2] >> 10)
        W[i] = (W[i-16] + s0 + W[i-7] + s1) & MASK32
    return W


def parity(x):
    """Parity (XOR of all bits)."""
    x ^= x >> 16; x ^= x >> 8; x ^= x >> 4
    x ^= x >> 2; x ^= x >> 1
    return x & 1


# ═══════════════════════════════════════════════════════════════
# EXPERIMENT 1: Linear bias in SHA-256 tail (rounds R→64)
# ═══════════════════════════════════════════════════════════════

def linear_bias_tail(start_round, samples=100000):
    """
    Measure linear bias of SHA-256 from round start_round to 64.

    For a RANDOM function: bias = 0 (up to sampling noise ~1/√N).
    For SHA-256: if bias > noise → we have a linear approximation.

    Test: input mask α (on state), output mask β (on hash).
    Bias = |P(α·x ⊕ β·f(x) = 0) - 1/2|

    We test many (α, β) pairs.
    """
    noise_floor = 1.0 / math.sqrt(4 * samples)

    # Generate random states at start_round and compute tail
    biases = []
    best_bias = 0
    best_masks = None

    # Test: single-bit input mask on register a or e
    # Single-bit output mask on hash word 0
    for in_reg in [0, 4]:  # a=0, e=4
        for in_bit in range(32):
            for out_bit in range(32):
                count = 0
                for _ in range(samples):
                    msg = [random.getrandbits(32) for _ in range(16)]
                    W = expand_schedule(msg)

                    # Compute state at start_round
                    state = sha256_rounds(list(IV), W, 0, start_round)

                    # Input parity
                    in_val = (state[in_reg] >> in_bit) & 1

                    # Compute tail
                    final = sha256_rounds(state, W, start_round, 64)
                    # Add IV (feedforward)
                    hash_val = [(IV[j] + final[j]) & MASK32 for j in range(8)]

                    # Output parity
                    out_val = (hash_val[0] >> out_bit) & 1

                    if in_val ^ out_val == 0:
                        count += 1

                bias = abs(count / samples - 0.5)
                if bias > best_bias:
                    best_bias = bias
                    reg_name = 'a' if in_reg == 0 else 'e'
                    best_masks = (f"{reg_name}[{in_bit}]", f"H0[{out_bit}]")

                biases.append(bias)

                # Early termination for large scans
                if in_bit > 4 and out_bit > 4 and best_bias < 2 * noise_floor:
                    break
            if in_bit > 4 and best_bias < 2 * noise_floor:
                break

    return best_bias, best_masks, noise_floor, biases


def scan_linear_bias(samples=50000):
    """Scan linear bias for different tail start points."""
    print("LINEAR BIAS IN SHA-256 TAIL")
    print("=" * 70)
    print()
    print(f"  Measure: max |P(α·state ⊕ β·hash = 0) - 1/2|")
    print(f"  For tail from round R to round 64")
    print(f"  Noise floor: 1/√(4N) where N = samples")
    print()
    print(f"  {'Start R':>8}  {'Best bias':>12}  {'Noise floor':>13}  "
          f"{'Ratio':>8}  {'Masks':>20}  {'Significant?':>14}")
    print(f"  {'-'*8}  {'-'*12}  {'-'*13}  {'-'*8}  {'-'*20}  {'-'*14}")

    for R in [1, 2, 4, 8, 12, 16, 20, 24, 32]:
        best, masks, noise, _ = linear_bias_tail(R, samples=samples)
        ratio = best / noise
        sig = "YES" if ratio > 3 else "no"

        masks_str = f"{masks[0]}→{masks[1]}" if masks else "—"
        print(f"  {R:>8}  {best:>12.6f}  {noise:>13.6f}  "
              f"{ratio:>8.1f}  {masks_str:>20}  {sig:>14}")

    print()


# ═══════════════════════════════════════════════════════════════
# EXPERIMENT 2: Multi-bit linear approximation
# ═══════════════════════════════════════════════════════════════

def multibit_linear(start_round=4, samples=50000):
    """
    Test multi-bit linear masks (XOR of multiple bits).
    These can have higher bias than single-bit.

    Specifically test masks aligned with Sigma/Ch/Maj structure.
    """
    print(f"MULTI-BIT LINEAR APPROXIMATION (R={start_round}→64)")
    print("=" * 70)
    print()

    noise = 1.0 / math.sqrt(4 * samples)

    # Predefined masks based on SHA-256 structure
    masks = [
        # Sigma1 structure: bits 6,11,25
        ("Σ1-aligned", lambda s: ((s[4]>>6)^(s[4]>>11)^(s[4]>>25)) & 1,
                        lambda h: (h[0]>>6) & 1),
        # Sigma0 structure: bits 2,13,22
        ("Σ0-aligned", lambda s: ((s[0]>>2)^(s[0]>>13)^(s[0]>>22)) & 1,
                        lambda h: (h[0]>>2) & 1),
        # Ch structure: e selects f vs g
        ("Ch-aligned", lambda s: ((s[4]&s[5])^(~s[4]&s[6])) & 1,
                        lambda h: h[0] & 1),
        # Parity of a
        ("parity(a)", lambda s: parity(s[0]),
                       lambda h: parity(h[0])),
        # Parity of e
        ("parity(e)", lambda s: parity(s[4]),
                       lambda h: parity(h[0])),
        # XOR a[0]⊕e[0]
        ("a[0]⊕e[0]", lambda s: (s[0] ^ s[4]) & 1,
                        lambda h: h[0] & 1),
        # Full word XOR: a⊕e bit 0
        ("(a⊕e)[15]", lambda s: ((s[0]^s[4])>>15) & 1,
                        lambda h: (h[0]>>15) & 1),
        # d register (kernel element)
        ("d[0]", lambda s: s[3] & 1,
                  lambda h: h[0] & 1),
        # h register (kernel element)
        ("h[0]", lambda s: s[7] & 1,
                  lambda h: h[0] & 1),
    ]

    print(f"  {'Mask':>15}  {'Bias':>10}  {'Noise':>10}  {'Ratio':>8}  {'Significant':>12}")
    print(f"  {'-'*15}  {'-'*10}  {'-'*10}  {'-'*8}  {'-'*12}")

    for name, in_mask, out_mask in masks:
        count = 0
        for _ in range(samples):
            msg = [random.getrandbits(32) for _ in range(16)]
            W = expand_schedule(msg)
            state = sha256_rounds(list(IV), W, 0, start_round)

            in_val = in_mask(state)
            final = sha256_rounds(state, W, start_round, 64)
            hash_val = [(IV[j] + final[j]) & MASK32 for j in range(8)]
            out_val = out_mask(hash_val)

            if in_val ^ out_val == 0:
                count += 1

        bias = abs(count / samples - 0.5)
        ratio = bias / noise
        sig = "YES" if ratio > 3 else "no"
        print(f"  {name:>15}  {bias:>10.6f}  {noise:>10.6f}  {ratio:>8.1f}  {sig:>12}")

    print()


# ═══════════════════════════════════════════════════════════════
# EXPERIMENT 3: Differential-Linear distinguisher
# ═══════════════════════════════════════════════════════════════

def difflinear_distinguisher(diff_rounds=4, samples=100000):
    """
    THE COMBINED ATTACK:
    1. Choose input pair (M, M') with DW[0] = 0x80000000
    2. Compute H = SHA(M), H' = SHA(M')
    3. Measure: P(β · (H ⊕ H') = 0) for various masks β

    For random oracle: P = 1/2 for any β
    For SHA-256: if differential creates bias at any β → distinguisher!

    This directly measures whether the Wang chain differential
    creates output correlations through the random-looking tail.
    """
    print(f"DIFFERENTIAL-LINEAR DISTINGUISHER")
    print(f"  Input diff: DW[0] = 0x80000000")
    print(f"  Full SHA-256 (64 rounds)")
    print("=" * 70)
    print()

    noise = 1.0 / math.sqrt(4 * samples)

    # For each output bit, measure P(DH[bit] = 0)
    bit_biases = []
    print(f"  {'Output bit':>12}  {'P(=0)':>10}  {'Bias':>10}  "
          f"{'Ratio':>8}  {'Status':>10}")
    print(f"  {'-'*12}  {'-'*10}  {'-'*10}  {'-'*8}  {'-'*10}")

    # Pre-generate all pairs
    bit_counts = [[0, 0] for _ in range(32)]  # [zero_count, one_count] per bit
    word_counts = Counter()

    for _ in range(samples):
        msg = [random.getrandbits(32) for _ in range(16)]
        msg2 = list(msg)
        msg2[0] ^= 0x80000000

        W1 = expand_schedule(msg)
        W2 = expand_schedule(msg2)
        h1 = sha256_rounds(list(IV), W1, 0, 64)
        h2 = sha256_rounds(list(IV), W2, 0, 64)

        hash1 = [(IV[j] + h1[j]) & MASK32 for j in range(8)]
        hash2 = [(IV[j] + h2[j]) & MASK32 for j in range(8)]

        dh0 = hash1[0] ^ hash2[0]
        word_counts[dh0] += 1

        for bit in range(32):
            if (dh0 >> bit) & 1 == 0:
                bit_counts[bit][0] += 1
            else:
                bit_counts[bit][1] += 1

    max_bias = 0
    for bit in range(32):
        p0 = bit_counts[bit][0] / samples
        bias = abs(p0 - 0.5)
        ratio = bias / noise
        status = "BIASED" if ratio > 3 else "random"
        bit_biases.append(bias)

        if bias > max_bias:
            max_bias = bias

        print(f"  H[0][{bit:>2}]     {p0:>10.6f}  {bias:>10.6f}  "
              f"{ratio:>8.1f}  {status:>10}")

    print()
    print(f"  Max single-bit bias: {max_bias:.6f} (noise={noise:.6f})")
    print(f"  Ratio: {max_bias/noise:.1f}")
    print()

    # Multi-bit: parity of output word
    parity_count = 0
    for dh, cnt in word_counts.items():
        if parity(dh) == 0:
            parity_count += cnt
    p_par = parity_count / samples
    par_bias = abs(p_par - 0.5)
    print(f"  Parity bias: {par_bias:.6f} (ratio={par_bias/noise:.1f})")
    print()

    # Unique output diffs
    print(f"  Unique DH values: {len(word_counts)} / {samples}")
    print(f"  Top 5 DH values:")
    for dh, cnt in word_counts.most_common(5):
        print(f"    0x{dh:08x}: {cnt} times ({cnt/samples:.6f})")
    print()

    return max_bias, noise


def difflinear_by_rounds(samples=50000):
    """
    Measure diff-linear bias for different number of total rounds.
    This shows how bias decays as more rounds are added.
    """
    print("DIFF-LINEAR BIAS vs NUMBER OF ROUNDS")
    print("=" * 70)
    print()

    noise = 1.0 / math.sqrt(4 * samples)
    print(f"  Noise floor: {noise:.6f}")
    print()
    print(f"  {'Rounds':>8}  {'Max bit bias':>14}  {'Ratio':>8}  "
          f"{'Parity bias':>14}  {'Ratio':>8}  {'Distinguisher?':>16}")
    print(f"  {'-'*8}  {'-'*14}  {'-'*8}  {'-'*14}  {'-'*8}  {'-'*16}")

    for total_rounds in [4, 8, 12, 16, 20, 24, 28, 32, 40, 48, 56, 64]:
        bit_counts = [[0, 0] for _ in range(32)]
        parity_count = 0

        for _ in range(samples):
            msg = [random.getrandbits(32) for _ in range(16)]
            msg2 = list(msg)
            msg2[0] ^= 0x80000000

            W1 = expand_schedule(msg)
            W2 = expand_schedule(msg2)
            h1 = sha256_rounds(list(IV), W1, 0, total_rounds)
            h2 = sha256_rounds(list(IV), W2, 0, total_rounds)

            hash1 = [(IV[j] + h1[j]) & MASK32 for j in range(8)]
            hash2 = [(IV[j] + h2[j]) & MASK32 for j in range(8)]

            dh0 = hash1[0] ^ hash2[0]
            if parity(dh0) == 0:
                parity_count += 1

            for bit in range(32):
                if (dh0 >> bit) & 1 == 0:
                    bit_counts[bit][0] += 1

        max_bias = max(abs(bit_counts[b][0]/samples - 0.5) for b in range(32))
        par_bias = abs(parity_count/samples - 0.5)
        mr = max_bias / noise
        pr = par_bias / noise

        dist = "YES" if mr > 3 or pr > 3 else "no"
        print(f"  {total_rounds:>8}  {max_bias:>14.6f}  {mr:>8.1f}  "
              f"{par_bias:>14.6f}  {pr:>8.1f}  {dist:>16}")

    print()


# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════

if __name__ == '__main__':
    print("STEP 27c: DIFFERENTIAL-LINEAR ATTACK ON SHA-256")
    print("Finding linear bias in the 'random' tail")
    print("=" * 70)
    print()

    t_start = time.time()

    # 1. Linear bias scan
    print("─" * 70)
    print("  PART 1: Pure linear bias in tail")
    print("─" * 70)
    scan_linear_bias(samples=20000)

    # 2. Multi-bit linear approximation
    print("─" * 70)
    print("  PART 2: Multi-bit (structure-aligned) linear masks")
    print("─" * 70)
    multibit_linear(start_round=4, samples=20000)
    multibit_linear(start_round=8, samples=20000)

    # 3. Diff-linear by rounds
    print("─" * 70)
    print("  PART 3: Diff-linear bias vs total rounds")
    print("─" * 70)
    difflinear_by_rounds(samples=30000)

    # 4. Full diff-linear on 64 rounds
    print("─" * 70)
    print("  PART 4: Full diff-linear distinguisher (64 rounds)")
    print("─" * 70)
    difflinear_distinguisher(samples=50000)

    elapsed = time.time() - t_start
    print("=" * 70)
    print(f"TIME: {elapsed:.1f}s")
    print()
    print("DIFFERENTIAL-LINEAR VERDICT:")
    print("  If ANY bias > noise at 64 rounds → SHA-256 ≠ random oracle")
    print("  If all biases = noise → random oracle model holds completely")
    print()
    print("  Implication for collision:")
    print("  Even bias ε = 2^{-10} → collision at 2^{118} instead of 2^{128}")
    print("  Bias ε = 2^{-20} → collision at 2^{108}")
    print("  No bias → 2^{128} is tight")
