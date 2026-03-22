#!/usr/bin/env python3
"""Step 14: Message schedule analysis for SHA-256.

Analyze sigma0/sigma1 dependencies in the message schedule
to find exploitable structure for extending differential control
beyond round 20.

Key questions:
1. What is the differential propagation through sigma0/sigma1?
2. Can we find message schedule differentials that cancel?
3. What constraints does the schedule impose on extended attacks?
"""

import random
import numpy as np
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

# SHA-256 round constants
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
    """Expand 16-word message to 64 words."""
    W = list(W16)
    for t in range(16, 64):
        W.append((sigma1(W[t-2]) + W[t-7] + sigma0(W[t-15]) + W[t-16]) & MASK32)
    return W

def sha256_compress(W64, iv=None):
    """Full SHA-256 compression with 64-word expanded message."""
    if iv is None:
        iv = IV
    a, b, c, d, e, f, g, h = iv
    states = [(a, b, c, d, e, f, g, h)]
    for t in range(64):
        T1 = (h + Sigma1(e) + Ch(e, f, g) + K[t] + W64[t]) & MASK32
        T2 = (Sigma0(a) + Maj(a, b, c)) & MASK32
        h = g
        g = f
        f = e
        e = (d + T1) & MASK32
        d = c
        c = b
        b = a
        a = (T1 + T2) & MASK32
        states.append((a, b, c, d, e, f, g, h))
    return states

def hw(x):
    return bin(x).count('1')

# ============================================================
# Analysis 1: Differential propagation through sigma0/sigma1
# ============================================================
def analyze_sigma_differential():
    """Study how single-bit differences propagate through sigma0/sigma1."""
    print("=" * 70)
    print("Analysis 1: Sigma differential propagation")
    print("=" * 70)

    for name, func in [("sigma0", sigma0), ("sigma1", sigma1)]:
        print(f"\n{name} single-bit differential HW:")
        hw_map = []
        for bit in range(32):
            delta = 1 << bit
            # Average over random x
            total_hw = 0
            N = 10000
            for _ in range(N):
                x = random.randint(0, MASK32)
                diff = (func((x + delta) & MASK32) - func(x)) & MASK32
                total_hw += hw(diff)
            avg_hw = total_hw / N
            hw_map.append(avg_hw)

        print(f"  Min HW: {min(hw_map):.1f} (bit {hw_map.index(min(hw_map))})")
        print(f"  Max HW: {max(hw_map):.1f} (bit {hw_map.index(max(hw_map))})")
        print(f"  Avg HW: {sum(hw_map)/32:.1f}")

        # XOR differential (GF2-linear part)
        print(f"\n{name} XOR differential HW (GF2-linear):")
        for bit in [0, 1, 15, 31]:
            delta = 1 << bit
            xor_diff = func(0 ^ delta) ^ func(0)  # Linear part
            print(f"  bit {bit:2d}: XOR diff HW = {hw(xor_diff)}, value = 0x{xor_diff:08x}")

# ============================================================
# Analysis 2: Message schedule dependency graph
# ============================================================
def analyze_schedule_dependencies():
    """Map which initial words W[0..15] affect each expanded word W[t]."""
    print("\n" + "=" * 70)
    print("Analysis 2: Message schedule dependency graph")
    print("=" * 70)

    # For each W[t], track which W[i] (i<16) it depends on
    deps = {}
    for i in range(16):
        deps[i] = {i}

    for t in range(16, 64):
        deps[t] = deps[t-2] | deps[t-7] | deps[t-15] | deps[t-16]

    print("\nW[t] depends on these initial words:")
    for t in [16, 17, 18, 19, 20, 24, 28, 32, 48, 63]:
        dep_list = sorted(deps[t])
        print(f"  W[{t:2d}]: {dep_list} ({len(dep_list)} words)")

    # Find when all 16 words are mixed
    for t in range(16, 64):
        if len(deps[t]) == 16:
            print(f"\n  Full mixing at W[{t}] — all 16 initial words contribute")
            break

    return deps

# ============================================================
# Analysis 3: XOR-linear differential through message schedule
# ============================================================
def analyze_linear_schedule_diff():
    """Analyze how a single-word difference in W[i] propagates through schedule."""
    print("\n" + "=" * 70)
    print("Analysis 3: XOR-linear schedule differential propagation")
    print("=" * 70)

    # sigma0 and sigma1 are GF2-linear, so XOR differential is exact
    # But addition is not GF2-linear — carries make it nonlinear

    for src_word in [0, 1, 14, 15]:
        print(f"\n  DW[{src_word}] = 0x80000000 (MSB flip):")

        # Compute schedule with and without difference
        N_trials = 5000
        hw_stats = defaultdict(list)

        for _ in range(N_trials):
            W = [random.randint(0, MASK32) for _ in range(16)]
            W2 = list(W)
            W2[src_word] ^= 0x80000000

            W_exp = expand_message(W)
            W2_exp = expand_message(W2)

            for t in range(16, 40):
                diff = (W_exp[t] - W2_exp[t]) & MASK32
                xdiff = W_exp[t] ^ W2_exp[t]
                hw_stats[t].append(hw(xdiff))

        for t in range(16, 35):
            hws = hw_stats[t]
            avg = sum(hws) / len(hws)
            if avg > 0.01:
                print(f"    W[{t:2d}]: avg XOR HW = {avg:.1f}")

# ============================================================
# Analysis 4: Can we find schedule differentials that stay sparse?
# ============================================================
def find_sparse_schedule_diffs():
    """Search for message differences that keep schedule differentials sparse."""
    print("\n" + "=" * 70)
    print("Analysis 4: Sparse schedule differentials")
    print("=" * 70)

    best_results = []

    for trial in range(100000):
        # Try random single-word XOR difference
        word_idx = random.randint(0, 15)
        # Try low HW differences
        n_bits = random.choice([1, 1, 1, 2, 2, 3])
        delta = 0
        for _ in range(n_bits):
            delta ^= 1 << random.randint(0, 31)

        # Evaluate over several random messages
        total_expansion_hw = 0
        N = 50
        for _ in range(N):
            W = [random.randint(0, MASK32) for _ in range(16)]
            W2 = list(W)
            W2[word_idx] = (W2[word_idx] + delta) & MASK32  # Additive difference

            W_exp = expand_message(W)
            W2_exp = expand_message(W2)

            # Total HW of XOR differences from round 16 to 40
            for t in range(16, 40):
                total_expansion_hw += hw(W_exp[t] ^ W2_exp[t])

        avg_hw = total_expansion_hw / N
        best_results.append((avg_hw, word_idx, delta, n_bits))

    best_results.sort()

    print("\nTop 10 sparsest additive schedule differentials (rounds 16-40):")
    seen = set()
    count = 0
    for avg_hw, word_idx, delta, n_bits in best_results:
        key = (word_idx, delta)
        if key in seen:
            continue
        seen.add(key)
        print(f"  DW[{word_idx:2d}] = 0x{delta:08x} (HW={hw(delta)}): "
              f"avg expansion XOR HW = {avg_hw:.1f}")
        count += 1
        if count >= 10:
            break

# ============================================================
# Analysis 5: Multi-word differences for schedule cancellation
# ============================================================
def find_schedule_cancellation():
    """Find multi-word differences where schedule terms cancel."""
    print("\n" + "=" * 70)
    print("Analysis 5: Multi-word schedule cancellation")
    print("=" * 70)

    # Key insight: W[t] = sigma1(W[t-2]) + W[t-7] + sigma0(W[t-15]) + W[t-16]
    # If we set differences in multiple words, terms can cancel

    # Strategy: set DW[i] and DW[j] such that their contributions
    # to later W[t] cancel out

    # First: purely linear analysis (ignore carries)
    # sigma0 and sigma1 are GF2-linear
    # We want: sigma0(DW[i]) XOR DW[j] = 0 for some relationship

    print("\nLinear cancellation pairs:")
    print("If DW[t-16] = sigma0(DW[t-15]), schedule diff at t cancels one term")

    # For DW[1] = delta, DW[0] = sigma0(delta):
    for bit in [0, 1, 7, 15, 31]:
        delta = 1 << bit
        cancel_val = sigma0(delta)
        print(f"  DW[1] = 1<<{bit:2d}, DW[0] = sigma0(DW[1]) = 0x{cancel_val:08x} (HW={hw(cancel_val)})")

    # Test actual cancellation with carries
    print("\nActual cancellation test (with carries):")
    for bit in [0, 15, 31]:
        delta = 1 << bit
        s0_delta = sigma0(delta)

        N = 10000
        cancel_count = 0
        total_hw = 0

        for _ in range(N):
            W = [random.randint(0, MASK32) for _ in range(16)]
            W2 = list(W)
            W2[0] = (W2[0] + s0_delta) & MASK32  # DW[0] = sigma0(delta)
            W2[1] = (W2[1] + delta) & MASK32      # DW[1] = delta

            W_exp = expand_message(W)
            W2_exp = expand_message(W2)

            # Check W[16] = sigma1(W[14]) + W[9] + sigma0(W[1]) + W[0]
            # DW[16] should have sigma0(delta) from W[1] and s0_delta from W[0]
            diff16 = (W_exp[16] - W2_exp[16]) & MASK32
            xdiff16 = W_exp[16] ^ W2_exp[16]
            total_hw += hw(xdiff16)
            if diff16 == 0:
                cancel_count += 1

        avg_hw = total_hw / N
        print(f"  bit {bit:2d}: DW[16] cancel rate = {cancel_count/N:.4f}, avg XOR HW = {avg_hw:.1f}")

# ============================================================
# Analysis 6: Full differential path probability estimation
# ============================================================
def estimate_path_probability():
    """Estimate probability of extending De=0 beyond round 20."""
    print("\n" + "=" * 70)
    print("Analysis 6: Differential path probability (rounds 20-64)")
    print("=" * 70)

    # Start from our known De17=0 solution setup
    # DW[0] = 0x80000000, rest adjustable

    N = 50000
    de_hw_by_round = defaultdict(list)
    da_hw_by_round = defaultdict(list)

    for _ in range(N):
        W = [random.randint(0, MASK32) for _ in range(16)]
        W[0] = random.randint(0, MASK32)

        W2 = list(W)
        W2[0] ^= 0x80000000  # MSB difference

        W_exp = expand_message(W)
        W2_exp = expand_message(W2)

        states1 = sha256_compress(W_exp)
        states2 = sha256_compress(W2_exp)

        for t in range(64):
            de = (states1[t+1][4] - states2[t+1][4]) & MASK32
            da = (states1[t+1][0] - states2[t+1][0]) & MASK32
            de_hw_by_round[t].append(hw(de))
            da_hw_by_round[t].append(hw(da))

    print("\nAverage differential HW by round (DW[0] = 0x80000000):")
    print(f"{'Round':>5} | {'De avg HW':>9} | {'Da avg HW':>9} | {'De=0 count':>10}")
    print("-" * 45)

    for t in range(64):
        de_avg = sum(de_hw_by_round[t]) / N
        da_avg = sum(da_hw_by_round[t]) / N
        de_zero = de_hw_by_round[t].count(0)
        if t < 25 or t % 4 == 0 or de_zero > 0:
            print(f"{t:5d} | {de_avg:9.2f} | {da_avg:9.2f} | {de_zero:10d}")

# ============================================================
# Analysis 7: Neutral bit candidates
# ============================================================
def find_neutral_bits():
    """Find message bits that don't affect the differential path."""
    print("\n" + "=" * 70)
    print("Analysis 7: Neutral bit candidates")
    print("=" * 70)

    # For each bit of each message word, check if flipping it
    # changes the De differential at rounds 1-20

    DW0 = 0x80000000  # Our known difference
    N = 2000

    neutral_scores = {}  # (word, bit) -> rounds where neutral

    for word_idx in range(16):
        for bit in range(32):
            flip = 1 << bit
            neutral_rounds = 0

            for _ in range(N):
                W = [random.randint(0, MASK32) for _ in range(16)]
                W2 = list(W)
                W2[0] ^= DW0

                # Flip bit in both W and W2 (same message modification)
                W_flip = list(W)
                W_flip[word_idx] ^= flip
                W2_flip = list(W2)
                W2_flip[word_idx] ^= flip

                W_exp = expand_message(W)
                W2_exp = expand_message(W2)
                Wf_exp = expand_message(W_flip)
                W2f_exp = expand_message(W2_flip)

                states1 = sha256_compress(W_exp)
                states2 = sha256_compress(W2_exp)
                states1f = sha256_compress(Wf_exp)
                states2f = sha256_compress(W2f_exp)

                # Check if differential is preserved
                max_neutral_round = 0
                for t in range(20):
                    de_orig = (states1[t+1][4] - states2[t+1][4]) & MASK32
                    de_flip = (states1f[t+1][4] - states2f[t+1][4]) & MASK32
                    if de_orig == de_flip:
                        max_neutral_round = t + 1
                    else:
                        break
                neutral_rounds += max_neutral_round

            avg_neutral = neutral_rounds / N
            neutral_scores[(word_idx, bit)] = avg_neutral

    # Sort by neutrality
    ranked = sorted(neutral_scores.items(), key=lambda x: -x[1])

    print("\nTop 30 neutral bits (preserve De differential longest):")
    print(f"{'Word':>4} {'Bit':>3} | {'Avg neutral rounds':>18}")
    print("-" * 30)
    for (word, bit), score in ranked[:30]:
        if score > 1.0:
            print(f"W[{word:2d}] b{bit:2d} | {score:18.1f}")

    # Count bits neutral for at least N rounds
    for threshold in [5, 10, 15, 17, 20]:
        count = sum(1 for s in neutral_scores.values() if s >= threshold)
        print(f"\n  Bits neutral for >= {threshold} rounds: {count}")

# ============================================================
# Main
# ============================================================
if __name__ == "__main__":
    random.seed(42)

    # Run analyses
    analyze_sigma_differential()
    deps = analyze_schedule_dependencies()
    analyze_linear_schedule_diff()
    find_sparse_schedule_diffs()
    find_schedule_cancellation()
    estimate_path_probability()
    find_neutral_bits()

    print("\n" + "=" * 70)
    print("Step 14 Complete")
    print("=" * 70)
