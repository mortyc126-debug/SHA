#!/usr/bin/env python3
"""
ЗАДАНИЕ 13, Часть 2: H[7]-CASCADE — масштабная проверка
Generate 2M random W[0] (W[1..15]=0), find H[7] collisions,
measure HW(ΔH[i]) compression. Null control with random pairs.
"""
import random, time, math
from collections import defaultdict

MASK = 0xFFFFFFFF

K = [
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2,
]
H0 = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]

def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK
def _sig0(x):   return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def _sig1(x):   return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)
def Sig0(x):    return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sig1(x):    return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return ((e & f) ^ ((~e & MASK) & g)) & MASK
def Maj(a, b, c): return ((a & b) ^ (a & c) ^ (b & c)) & MASK
def add(*args):
    s = 0
    for x in args: s = (s + x) & MASK
    return s
def hw(x): return bin(x & MASK).count('1')

def schedule(M):
    W = list(M[:16])
    for i in range(16, 64):
        W.append(add(_sig1(W[i-2]), W[i-7], _sig0(W[i-15]), W[i-16]))
    return W

def sha256_hash(W):
    """Full SHA-256 of 16-word message, returns 8 words."""
    s = list(H0)
    Ws = schedule(W)
    for r in range(64):
        T1 = add(s[7], Sig1(s[4]), Ch(s[4], s[5], s[6]), K[r], Ws[r])
        T2 = add(Sig0(s[0]), Maj(s[0], s[1], s[2]))
        s = [add(T1, T2), s[0], s[1], s[2], add(s[3], T1), s[4], s[5], s[6]]
    return [add(H0[i], s[i]) for i in range(8)]

def main():
    random.seed(42)
    print("=" * 70)
    print("ЗАДАНИЕ 13, Часть 2: H[7]-CASCADE")
    print("=" * 70)

    # ── Generate samples and find H[7] collisions ──
    N = 2000000
    print(f"\nGenerating {N} random W[0] values (W[1..15]=0)...")
    t0 = time.time()

    # Store H[7] → list of (w0, full_hash)
    h7_table = defaultdict(list)

    batch_size = 100000
    total_collisions = 0
    collision_pairs = []  # list of (Ha, Hb) full hashes

    for batch_start in range(0, N, batch_size):
        batch_end = min(batch_start + batch_size, N)
        for i in range(batch_start, batch_end):
            w0 = random.getrandbits(32)
            W = [w0] + [0] * 15
            H = sha256_hash(W)
            h7 = H[7]

            # Check for collision with existing entries (exclude same W[0])
            if h7 in h7_table:
                for (prev_w0, prev_H) in h7_table[h7]:
                    if prev_w0 != w0:  # Skip trivial W[0] duplicates
                        collision_pairs.append((prev_H, H))
                        total_collisions += 1

            h7_table[h7].append((w0, H))

        elapsed = time.time() - t0
        rate = (batch_end) / elapsed if elapsed > 0 else 0
        print(f"  Processed {batch_end}/{N} ({rate:.0f}/s), collisions so far: {total_collisions}")

    elapsed = time.time() - t0
    print(f"\nTotal H[7] collisions found: {total_collisions}")
    print(f"Time: {elapsed:.1f}s")
    expected = N * (N - 1) / (2 * (2**32))
    print(f"Expected (birthday): ~{expected:.1f}")

    if total_collisions == 0:
        print("ERROR: No H[7] collisions found. Cannot proceed.")
        return

    # ── Analyze HW(ΔH[i]) for collision pairs ──
    print(f"\n{'─'*70}")
    print(f"Analyzing HW(ΔH[i]) for {len(collision_pairs)} H[7]-collision pairs")
    print(f"{'─'*70}")

    hw_sums = [0.0] * 8
    hw_sq_sums = [0.0] * 8
    n_pairs = len(collision_pairs)

    # Track best (min HW(ΔH[4])) pair
    best_h4_hw = 33
    best_pair = None

    for (Ha, Hb) in collision_pairs:
        for i in range(8):
            d = hw(Ha[i] ^ Hb[i])
            hw_sums[i] += d
            hw_sq_sums[i] += d * d
            if i == 4 and d < best_h4_hw:
                best_h4_hw = d
                best_pair = (Ha, Hb)

    print(f"\n{'i':>3} | {'mean HW(ΔH[i])':>15} | {'std':>8} | {'compression':>12} | N_pairs")
    print(f"{'─'*3}-+-{'─'*15}-+-{'─'*8}-+-{'─'*12}-+-{'─'*8}")
    for i in range(8):
        mean_hw = hw_sums[i] / n_pairs
        var_hw = hw_sq_sums[i] / n_pairs - mean_hw ** 2
        std_hw = math.sqrt(max(0, var_hw))
        compression = 16.0 - mean_hw
        marker = " ← H[7]=0 (collision)" if i == 7 else ""
        print(f"  {i} | {mean_hw:15.4f} | {std_hw:8.4f} | {compression:12.4f} | {n_pairs}{marker}")

    # ── NULL CONTROL ──
    print(f"\n{'─'*70}")
    print("NULL CONTROL: 500 random (non-H[7]-collisional) pairs")
    print(f"{'─'*70}")

    null_hw_sums = [0.0] * 8
    null_hw_sq_sums = [0.0] * 8
    n_null = 500

    for _ in range(n_null):
        w0a = random.getrandbits(32)
        w0b = random.getrandbits(32)
        while w0a == w0b:
            w0b = random.getrandbits(32)
        Wa = [w0a] + [0] * 15
        Wb = [w0b] + [0] * 15
        Ha = sha256_hash(Wa)
        Hb = sha256_hash(Wb)
        for i in range(8):
            d = hw(Ha[i] ^ Hb[i])
            null_hw_sums[i] += d
            null_hw_sq_sums[i] += d * d

    print(f"\n{'i':>3} | {'mean HW(ΔH[i])':>15} | {'std':>8} | {'deviation from 16':>18}")
    print(f"{'─'*3}-+-{'─'*15}-+-{'─'*8}-+-{'─'*18}")
    for i in range(8):
        mean_hw = null_hw_sums[i] / n_null
        var_hw = null_hw_sq_sums[i] / n_null - mean_hw ** 2
        std_hw = math.sqrt(max(0, var_hw))
        dev = mean_hw - 16.0
        print(f"  {i} | {mean_hw:15.4f} | {std_hw:8.4f} | {dev:18.4f}")

    # ── Best H[4] near-collision ──
    if best_pair is not None:
        Ha, Hb = best_pair
        print(f"\n{'─'*70}")
        print(f"Best two-word near-collision: H[7] exact + H[4] near")
        print(f"  HW(ΔH[4]) = {best_h4_hw}")
        print(f"  Ha = {' '.join(f'{x:08x}' for x in Ha)}")
        print(f"  Hb = {' '.join(f'{x:08x}' for x in Hb)}")
        print(f"  ΔH  = ", end="")
        for i in range(8):
            d = Ha[i] ^ Hb[i]
            print(f"{d:08x}(hw={hw(d):2d}) ", end="")
        print()

    # ── Triple-check if compression > 2 bits (Rule 14) ──
    compression_h4 = 16.0 - hw_sums[4] / n_pairs
    if abs(compression_h4) > 2.0:
        print(f"\n{'='*70}")
        print(f"RULE 14 ALERT: H[4] compression = {compression_h4:.4f} > 2 bits")
        print(f"Running TRIPLE CHECK with 3 independent seeds...")
        print(f"{'='*70}")

        for check_idx, seed in enumerate([1337, 2024, 9999]):
            random.seed(seed)
            check_n = 500000
            ch7 = defaultdict(list)
            cpairs = []
            for _ in range(check_n):
                w0 = random.getrandbits(32)
                W = [w0] + [0] * 15
                H = sha256_hash(W)
                h7 = H[7]
                if h7 in ch7:
                    for (prev_w0, pH) in ch7[h7]:
                        if prev_w0 != w0:
                            cpairs.append((pH, H))
                ch7[h7].append((w0, H))

            if cpairs:
                mean_all = [sum(hw(a[i] ^ b[i]) for a, b in cpairs) / len(cpairs) for i in range(8)]
                print(f"  Check {check_idx+1} (seed={seed}): {len(cpairs)} genuine pairs")
                print(f"    mean HW(ΔH[i]): {' '.join(f'{m:.2f}' for m in mean_all)}")
            else:
                print(f"  Check {check_idx+1} (seed={seed}): 0 genuine pairs found in {check_n} samples")

    print(f"\n{'='*70}")
    print("Part 2 COMPLETE")
    print(f"{'='*70}")

if __name__ == "__main__":
    main()
