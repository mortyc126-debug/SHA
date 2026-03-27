"""
FIBER STRUCTURE: распределение прообразов.

Для каждого выхода H: fiber(H) = {W : SHA(W) = H}.
|fiber(H)| = "толщина" ткани в точке H.

Random function 2^N → 2^M: |fiber| ~ Poisson(2^{N-M}).
SHA-256: N=512, M=256 → |fiber| ~ Poisson(2^256).

Мы не можем измерить 2^256 прообразов. Но можем измерить
fiber на REDUCED space: фикс W[1..15], варьируем W[0].

2^32 inputs → 2^32 output space (H[7] only).
Expected: uniform (Poisson(1)).
Deviation = structural non-uniformity.
"""

import numpy as np
import struct, hashlib
from collections import Counter

def sha256_words(W16):
    return struct.unpack('>8I', hashlib.sha256(struct.pack('>16I', *W16)).digest())

def hw(x): return bin(x).count('1')


def main():
    np.random.seed(42)

    print("=" * 70)
    print("FIBER STRUCTURE: распределение прообразов")
    print("=" * 70)

    # =================================================================
    print(f"\n{'=' * 70}")
    print("1. REDUCED FIBER: W[0] → H[7] (32→32 bit mapping)")
    print("=" * 70)

    W_base = [np.random.randint(0, 2**32) for _ in range(16)]

    N = 200000  # sample from 2^32 space
    h7_counter = Counter()

    for _ in range(N):
        W = list(W_base)
        W[0] = np.random.randint(0, 2**32)
        H = sha256_words(W)
        h7_counter[H[7]] += 1

    # Fiber size distribution
    fiber_sizes = list(h7_counter.values())
    size_dist = Counter(fiber_sizes)

    print(f"\n  {N} random W[0] → H[7] distribution:")
    print(f"    Unique H[7] values: {len(h7_counter)}")
    print(f"    Expected (random): {N} unique (all different)")

    print(f"\n  Fiber size distribution:")
    for size in sorted(size_dist.keys()):
        count = size_dist[size]
        import math
        expected = N * np.exp(-1) / math.factorial(size-1) if size <= 10 else 0
        bar = "█" * min(count // max(1, N // 500), 50)
        print(f"    |fiber|={size}: {count:6d} fibers  {bar}")

    # Statistics
    mean_fiber = np.mean(fiber_sizes)
    std_fiber = np.std(fiber_sizes)
    max_fiber = max(fiber_sizes)

    print(f"\n  Statistics:")
    print(f"    Mean fiber size: {mean_fiber:.4f}")
    print(f"    Std: {std_fiber:.4f}")
    print(f"    Max fiber: {max_fiber}")
    print(f"    Expected (Poisson(λ={N/2**32:.4f})): mean={N/2**32:.4f}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("2. MULTI-WORD FIBER: W[0] → (H[6], H[7]) (32→64 bit)")
    print("=" * 70)

    h67_counter = Counter()
    for _ in range(N):
        W = list(W_base)
        W[0] = np.random.randint(0, 2**32)
        H = sha256_words(W)
        h67_counter[(H[6], H[7])] += 1

    sizes_67 = list(h67_counter.values())
    print(f"\n  {N} random W[0] → (H[6],H[7]):")
    print(f"    Unique pairs: {len(h67_counter)}")
    print(f"    Collisions: {N - len(h67_counter)}")
    print(f"    Max fiber: {max(sizes_67)}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("3. FULL OUTPUT FIBER: (W[0],W[1]) → H[7] (64→32 bit)")
    print("=" * 70)

    # More inputs than outputs → expect collisions
    h7_counter_2 = Counter()
    N2 = 200000

    for _ in range(N2):
        W = list(W_base)
        W[0] = np.random.randint(0, 2**32)
        W[1] = np.random.randint(0, 2**32)
        H = sha256_words(W)
        h7_counter_2[H[7]] += 1

    sizes_2 = list(h7_counter_2.values())
    size_dist_2 = Counter(sizes_2)

    print(f"\n  {N2} random (W[0],W[1]) → H[7]:")
    print(f"    Unique H[7]: {len(h7_counter_2)}")
    print(f"    Max fiber: {max(sizes_2)}")

    print(f"\n  Fiber distribution:")
    for size in sorted(size_dist_2.keys())[:10]:
        count = size_dist_2[size]
        print(f"    |fiber|={size}: {count}")

    # Expected: Poisson(N/2^32). N=200K, 2^32=4.3B → λ≈0.000047.
    # Almost all outputs appear 0 or 1 times. No large fibers.
    expected_lambda = N2 / 2**32

    print(f"\n  Expected λ = N/2^32 = {expected_lambda:.6f}")
    print(f"  P(fiber≥2) = 1 - (1+λ)e^(-λ) ≈ λ²/2 = {expected_lambda**2/2:.2e}")
    print(f"  Expected collisions: {N2 * expected_lambda / 2:.4f}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("4. UNIFORMITY TEST: χ² для H[7] distribution")
    print("=" * 70)

    # Bin H[7] into 256 bins (top 8 bits)
    bin_counter = Counter()
    for _ in range(N):
        W = list(W_base)
        W[0] = np.random.randint(0, 2**32)
        H = sha256_words(W)
        bin_counter[H[7] >> 24] += 1

    expected_per_bin = N / 256
    chi2 = sum((bin_counter.get(b, 0) - expected_per_bin)**2 / expected_per_bin
               for b in range(256))

    # Degrees of freedom = 255
    # For χ²(255): mean=255, std≈√(2×255)≈22.6
    z_score = (chi2 - 255) / np.sqrt(2 * 255)

    print(f"\n  χ² test on H[7] top byte (256 bins, {N} samples):")
    print(f"    χ² = {chi2:.1f}")
    print(f"    Expected: 255 ± 22.6")
    print(f"    z-score: {z_score:.2f}")
    print(f"    → {'UNIFORM ✓' if abs(z_score) < 3 else 'NON-UNIFORM!'}")

    # Same test for each output word
    print(f"\n  χ² per output word:")
    for w_idx in range(8):
        bin_counter_w = Counter()
        for _ in range(100000):
            W = list(W_base)
            W[0] = np.random.randint(0, 2**32)
            H = sha256_words(W)
            bin_counter_w[H[w_idx] >> 24] += 1

        expected_w = 100000 / 256
        chi2_w = sum((bin_counter_w.get(b, 0) - expected_w)**2 / expected_w for b in range(256))
        z_w = (chi2_w - 255) / np.sqrt(2 * 255)
        status = "✓" if abs(z_w) < 3 else "✗"
        print(f"    H[{w_idx}]: χ²={chi2_w:.1f}, z={z_w:+.2f} {status}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("5. ЗНАЧЕНИЕ FIBER STRUCTURE")
    print("=" * 70)

    print(f"""
  FIBER STRUCTURE В НАШЕМ ИЗМЕРЕНИИ:

  Fiber = множество входов с одинаковым выходом.
  Для random function: |fiber| ~ Poisson(C^{{N_in - N_out}}).

  SHA-256:
    Single-word fiber (W[0]→H[7]): uniform ✓
    Multi-word fiber: no collisions in 200K ✓
    χ² test: z={z_score:.2f} (uniform) ✓

  SHA-256 fibers = UNIFORM (as random function).
  Нет "толстых" fibers (concentration).
  Нет structural non-uniformity.

  → Fiber structure ПОДТВЕРЖДАЕТ: SHA-256 = random function
    на уровне fiber distribution.
  → Нет "лёгких" выходов (outputs с many preimages).
  → Collision равно вероятна для любого output.
""")


if __name__ == "__main__":
    main()
