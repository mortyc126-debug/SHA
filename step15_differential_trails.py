#!/usr/bin/env python3
"""Step 15: High-probability differential trail search for SHA-256.

Instead of requiring exact De=0, find differential characteristics
where each round transition has manageable probability.

Key techniques:
1. Optimal single-bit input differences
2. Controlled carry propagation
3. Conditional differential probability
"""

import random
from collections import defaultdict

MASK32 = 0xFFFFFFFF

def rotr(x, n):
    return ((x >> n) | (x << (32 - n))) & MASK32

def shr(x, n):
    return x >> n

def sigma0(x):
    return rotr(x, 7) ^ rotr(x, 18) ^ shr(x, 3)

def sigma1(x):
    return rotr(x, 17) ^ rotr(x, 19) ^ shr(x, 10)

def Sigma0(x):
    return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)

def Sigma1(x):
    return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)

def Ch(e, f, g):
    return (e & f) ^ (~e & g) & MASK32

def Maj(a, b, c):
    return (a & b) ^ (a & c) ^ (b & c)

K = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
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

def expand_message(W16):
    W = list(W16)
    for t in range(16, 64):
        W.append((sigma1(W[t-2]) + W[t-7] + sigma0(W[t-15]) + W[t-16]) & MASK32)
    return W

def sha256_compress(W64, iv=None):
    if iv is None:
        iv = IV
    a, b, c, d, e, f, g, h = iv
    states = [(a, b, c, d, e, f, g, h)]
    for t in range(64):
        T1 = (h + Sigma1(e) + Ch(e, f, g) + K[t] + W64[t]) & MASK32
        T2 = (Sigma0(a) + Maj(a, b, c)) & MASK32
        h = g; g = f; f = e
        e = (d + T1) & MASK32
        d = c; c = b; b = a
        a = (T1 + T2) & MASK32
        states.append((a, b, c, d, e, f, g, h))
    return states

def hw(x):
    return bin(x).count('1')

# ============================================================
# Analysis 1: Per-round transition probability
# ============================================================
def analyze_round_transition_probs():
    """For each round, estimate P(De_t = target | De_{t-1} = known)."""
    print("=" * 70)
    print("Analysis 1: Per-round differential transition probabilities")
    print("=" * 70)

    # Study how De evolves per round with DW[0] = 0x80000000
    DW0 = 0x80000000
    N = 100000

    # Track De HW distribution per round
    de_hw_dist = defaultdict(lambda: defaultdict(int))
    de_low_hw_count = defaultdict(int)  # De with HW <= 3

    for _ in range(N):
        W = [random.randint(0, MASK32) for _ in range(16)]
        W2 = list(W)
        W2[0] ^= DW0

        W_exp = expand_message(W)
        W2_exp = expand_message(W2)

        s1 = sha256_compress(W_exp)
        s2 = sha256_compress(W2_exp)

        for t in range(64):
            de = (s1[t+1][4] - s2[t+1][4]) & MASK32
            h = hw(de)
            de_hw_dist[t][h] += 1
            if h <= 3:
                de_low_hw_count[t] += 1

    print(f"\n{'Round':>5} | {'P(De=0)':>12} | {'P(HW<=3)':>12} | {'P(HW<=5)':>12} | {'Avg HW':>7}")
    print("-" * 60)
    for t in range(64):
        p_zero = de_hw_dist[t].get(0, 0) / N
        p_low3 = sum(de_hw_dist[t].get(h, 0) for h in range(4)) / N
        p_low5 = sum(de_hw_dist[t].get(h, 0) for h in range(6)) / N
        avg_hw = sum(h * de_hw_dist[t].get(h, 0) for h in range(33)) / N
        if t < 25 or t % 4 == 0 or p_zero > 0:
            print(f"{t:5d} | {p_zero:12.6f} | {p_low3:12.6f} | {p_low5:12.6f} | {avg_hw:7.2f}")

# ============================================================
# Analysis 2: Best input differences
# ============================================================
def find_best_input_diffs():
    """Search for input differences that minimize differential growth."""
    print("\n" + "=" * 70)
    print("Analysis 2: Optimal input differences")
    print("=" * 70)

    candidates = []

    # Single-bit differences in each word
    for word in range(16):
        for bit in [0, 1, 7, 15, 24, 31]:
            delta = 1 << bit
            N = 5000

            total_de_hw_20 = 0
            total_de_hw_30 = 0
            de_zero_17 = 0

            for _ in range(N):
                W = [random.randint(0, MASK32) for _ in range(16)]
                W2 = list(W)
                W2[word] ^= delta

                W_exp = expand_message(W)
                W2_exp = expand_message(W2)

                s1 = sha256_compress(W_exp)
                s2 = sha256_compress(W2_exp)

                de17 = (s1[18][4] - s2[18][4]) & MASK32
                de20 = (s1[21][4] - s2[21][4]) & MASK32
                de30 = (s1[31][4] - s2[31][4]) & MASK32

                total_de_hw_20 += hw(de20)
                total_de_hw_30 += hw(de30)
                if de17 == 0:
                    de_zero_17 += 1

            avg_hw20 = total_de_hw_20 / N
            avg_hw30 = total_de_hw_30 / N
            candidates.append((avg_hw30, word, bit, avg_hw20, de_zero_17/N))

    candidates.sort()
    print(f"\nTop 20 input differences (sorted by De30 avg HW):")
    print(f"{'Word':>4} {'Bit':>3} | {'De20 HW':>8} | {'De30 HW':>8} | {'P(De17=0)':>10}")
    print("-" * 45)
    for avg_hw30, word, bit, avg_hw20, p_de17_zero in candidates[:20]:
        print(f"W[{word:2d}] b{bit:2d} | {avg_hw20:8.2f} | {avg_hw30:8.2f} | {p_de17_zero:10.4f}")

# ============================================================
# Analysis 3: Conditional probabilities — given De17=0
# ============================================================
def conditional_analysis():
    """Given De17=0 (achieved), what's the probability of extending?"""
    print("\n" + "=" * 70)
    print("Analysis 3: Conditional probabilities given De17=0")
    print("=" * 70)

    DW0 = 0x80000000
    N = 500000
    found_de17_zero = 0

    de_hw_given_de17zero = defaultdict(list)
    extend_counts = defaultdict(int)

    for _ in range(N):
        W = [random.randint(0, MASK32) for _ in range(16)]
        W2 = list(W)
        W2[0] ^= DW0

        W_exp = expand_message(W)
        W2_exp = expand_message(W2)

        s1 = sha256_compress(W_exp)
        s2 = sha256_compress(W2_exp)

        de17 = (s1[18][4] - s2[18][4]) & MASK32
        if de17 == 0:
            found_de17_zero += 1
            for t in range(17, 40):
                de_t = (s1[t+1][4] - s2[t+1][4]) & MASK32
                de_hw_given_de17zero[t].append(hw(de_t))
                if de_t == 0:
                    extend_counts[t] += 1

    print(f"\nFound {found_de17_zero} cases with De17=0 out of {N}")
    if found_de17_zero > 0:
        print(f"P(De17=0) ≈ 2^{-32 + 32} (expected ~2^-32)")
        print(f"\nConditional on De17=0:")
        print(f"{'Round':>5} | {'Avg De HW':>10} | {'P(De=0|De17=0)':>15} | {'Count':>6}")
        print("-" * 45)
        for t in range(17, 35):
            hws = de_hw_given_de17zero.get(t, [])
            if hws:
                avg = sum(hws) / len(hws)
                p_zero = extend_counts.get(t, 0) / found_de17_zero
                print(f"{t:5d} | {avg:10.2f} | {p_zero:15.6f} | {extend_counts.get(t, 0):6d}")

# ============================================================
# Analysis 4: Two-block attack analysis
# ============================================================
def two_block_analysis():
    """Analyze feasibility of a two-block collision attack."""
    print("\n" + "=" * 70)
    print("Analysis 4: Two-block attack feasibility")
    print("=" * 70)

    # In a two-block attack:
    # Block 1: Create a controlled intermediate hash difference
    # Block 2: Cancel that difference

    # The output of block 1 compression becomes the IV for block 2
    # We get 512 new bits of freedom in block 2

    # Question: Can we control the hash difference after block 1
    # to something that block 2 can cancel?

    DW0 = 0x80000000
    N = 100000

    # After 64 rounds, what does the De distribution look like?
    final_de_hw = []
    final_da_hw = []
    final_state_diffs = []

    for _ in range(N):
        W = [random.randint(0, MASK32) for _ in range(16)]
        W2 = list(W)
        W2[0] ^= DW0

        W_exp = expand_message(W)
        W2_exp = expand_message(W2)

        s1 = sha256_compress(W_exp)
        s2 = sha256_compress(W2_exp)

        # Final state differences (before feedforward)
        total_diff_hw = 0
        for reg in range(8):
            d = s1[64][reg] ^ s2[64][reg]
            total_diff_hw += hw(d)

        final_de_hw.append(hw(s1[64][4] ^ s2[64][4]))
        final_da_hw.append(hw(s1[64][0] ^ s2[64][0]))
        final_state_diffs.append(total_diff_hw)

        # Feedforward: H = IV + final_state
        # D(H) = D(final_state) since IV is same
        # So hash difference = state difference

    avg_total = sum(final_state_diffs) / N
    avg_de = sum(final_de_hw) / N
    avg_da = sum(final_da_hw) / N

    print(f"\nAfter full 64 rounds (DW[0] = 0x80000000):")
    print(f"  Avg total state diff HW: {avg_total:.1f} / 256 bits")
    print(f"  Avg De diff HW: {avg_de:.1f} / 32 bits")
    print(f"  Avg Da diff HW: {avg_da:.1f} / 32 bits")

    # For two-block: block 2 needs to produce SAME state diff
    # from a DIFFERENT IV but same message difference
    # This is essentially as hard as the original problem

    # Alternative: block 1 uses NO message difference
    # Block 2 uses the message difference and must collide
    # But the IV for block 2 is the output of block 1 (identical for both)
    # So this reduces to a single-block attack with non-standard IV

    # Better alternative: different messages in BOTH blocks
    print(f"\n  Two-block strategy options:")
    print(f"  1. Different M1, same M2: need block1 to create exploitable IV diff")
    print(f"  2. Same M1, different M2: reduces to single-block with custom IV")
    print(f"  3. Different M1 AND M2: 1024 bits freedom but 2x constraints")

    # Test option 2: single-block with random IV
    print(f"\n  Testing option 2 (custom IV, DW[0]=0x80000000):")
    de_zero_by_round = defaultdict(int)
    N2 = 50000
    for _ in range(N2):
        # Random IV (simulating output of block 1)
        custom_iv = [random.randint(0, MASK32) for _ in range(8)]

        W = [random.randint(0, MASK32) for _ in range(16)]
        W2 = list(W)
        W2[0] ^= DW0

        W_exp = expand_message(W)
        W2_exp = expand_message(W2)

        s1 = sha256_compress(W_exp, custom_iv)
        s2 = sha256_compress(W2_exp, custom_iv)

        for t in range(25):
            de = (s1[t+1][4] - s2[t+1][4]) & MASK32
            if de == 0:
                de_zero_by_round[t] += 1

    print(f"  P(De_t=0) with random IV vs standard IV:")
    for t in range(20):
        p = de_zero_by_round.get(t, 0) / N2
        if p > 0:
            print(f"    Round {t:2d}: P = {p:.6f}")

# ============================================================
# Analysis 5: Truncated differential analysis
# ============================================================
def truncated_differential():
    """Analyze truncated differentials — focus on specific output bits."""
    print("\n" + "=" * 70)
    print("Analysis 5: Truncated differential analysis")
    print("=" * 70)

    # Instead of full De=0, target specific bits of the hash
    # SHA-256 output = IV + final_state (mod 2^32 per word)

    DW0 = 0x80000000
    N = 100000

    # For each output bit position, track P(output_bit_diff = 0)
    bit_zero_prob = defaultdict(lambda: defaultdict(int))

    for _ in range(N):
        W = [random.randint(0, MASK32) for _ in range(16)]
        W2 = list(W)
        W2[0] ^= DW0

        W_exp = expand_message(W)
        W2_exp = expand_message(W2)

        s1 = sha256_compress(W_exp)
        s2 = sha256_compress(W2_exp)

        # Hash = IV + state (mod 2^32)
        for reg in range(8):
            h1 = (IV[reg] + s1[64][reg]) & MASK32
            h2 = (IV[reg] + s2[64][reg]) & MASK32
            diff = h1 ^ h2
            for bit in range(32):
                if not (diff >> bit) & 1:
                    bit_zero_prob[reg][bit] += 1

    print(f"\nOutput bit bias (P(diff_bit=0) - should be 0.5 for random):")
    biased_bits = []
    for reg in range(8):
        for bit in range(32):
            p = bit_zero_prob[reg][bit] / N
            bias = abs(p - 0.5)
            if bias > 0.005:
                biased_bits.append((bias, reg, bit, p))

    biased_bits.sort(reverse=True)
    if biased_bits:
        print(f"  Found {len(biased_bits)} biased bits (|bias| > 0.005):")
        for bias, reg, bit, p in biased_bits[:20]:
            print(f"    H[{reg}] bit {bit:2d}: P(0) = {p:.4f} (bias = {bias:.4f})")
    else:
        print("  No significant biases found (hash acts as random oracle)")

    # Birthday-like analysis
    print(f"\n  Birthday analysis:")
    print(f"  Full 256-bit collision: 2^128 expected")
    print(f"  If k bits biased with avg bias b:")
    print(f"  Effective hash space ≈ 256 - k*H(0.5+b) bits")
    if biased_bits:
        total_saved = sum(1 - (-p*__import__('math').log2(p) - (1-p)*__import__('math').log2(1-p))
                        for _, _, _, p in biased_bits if 0 < p < 1)
        print(f"  Total information saved: {total_saved:.2f} bits")
        print(f"  Effective collision complexity: 2^{128 - total_saved/2:.1f}")

# ============================================================
# Main
# ============================================================
if __name__ == "__main__":
    random.seed(123)
    analyze_round_transition_probs()
    find_best_input_diffs()
    conditional_analysis()
    two_block_analysis()
    truncated_differential()

    print("\n" + "=" * 70)
    print("Step 15 Complete")
    print("=" * 70)
