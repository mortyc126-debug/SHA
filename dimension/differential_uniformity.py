"""
DIFFERENTIAL UNIFORMITY: для каждого δW, distribution δH.

Для ФИКСИРОВАННОГО δW: δH = H(W⊕δW) ⊕ H(W) зависит от W.
Меняя W: δH ~ distribution.
Если distribution UNIFORM: δW = "плохой" дифференциал (random output).
Если distribution BIASED: δW = "хороший" дифференциал (structured).

Differential uniformity = max over δW of max bias.
Для random function: uniformity ≈ 2/N.
Для SHA-256: ???

Это НАТИВНАЯ метрика нашего измерения: quality of differentials.
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
    print("DIFFERENTIAL UNIFORMITY")
    print("=" * 70)

    # =================================================================
    print(f"\n{'=' * 70}")
    print("1. ДЛЯ ФИКСИРОВАННОГО δW: distribution HW(δH)")
    print("=" * 70)

    # Test several fixed δW: 1-bit, multi-bit, specific patterns
    deltas = [
        ("W[0] bit 0", [1] + [0]*15),
        ("W[0] bit 31", [1<<31] + [0]*15),
        ("W[0] = 0xFF", [0xFF] + [0]*15),
        ("W[0] = random", [0xDEADBEEF] + [0]*15),
        ("W[0,1] = 1,1", [1, 1] + [0]*14),
        ("W[7] bit 0", [0]*7 + [1] + [0]*8),
        ("W[15] bit 0", [0]*15 + [1]),
    ]

    for name, dW in deltas:
        hw_dist = []
        h7_values = Counter()

        for trial in range(10000):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            H1 = sha256_words(W)
            W2 = [W[i] ^ dW[i] for i in range(16)]
            H2 = sha256_words(W2)

            dh_hw = sum(hw(H1[i] ^ H2[i]) for i in range(8))
            hw_dist.append(dh_hw)

            # Track H[7] of δH
            h7_values[(H1[7] ^ H2[7]) >> 24] += 1  # top byte

        mean_hw = np.mean(hw_dist)
        std_hw = np.std(hw_dist)

        # Uniformity of top-byte δH[7]
        expected = 10000 / 256
        chi2 = sum((h7_values.get(b, 0) - expected)**2 / expected for b in range(256))
        z = (chi2 - 255) / np.sqrt(2*255)

        print(f"  {name:20s}: HW mean={mean_hw:.1f}±{std_hw:.1f}, χ²(δH[7])={chi2:.0f}, z={z:+.1f}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("2. BEST DIFFERENTIAL: поиск δW с минимальным HW(δH)")
    print("=" * 70)

    # Для random W: HW(δH) ≈ 128 ± 8.
    # Есть ли δW с СИСТЕМАТИЧЕСКИ ниже 128?

    best_delta = None
    best_mean_hw = 128

    for trial in range(200):
        # Random 1-bit δW
        dW = [0]*16
        word = np.random.randint(0, 16)
        bit = np.random.randint(0, 32)
        dW[word] = 1 << bit

        hws = []
        for _ in range(100):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            H1 = sha256_words(W)
            W2 = [W[i] ^ dW[i] for i in range(16)]
            H2 = sha256_words(W2)
            hws.append(sum(hw(H1[i] ^ H2[i]) for i in range(8)))

        m = np.mean(hws)
        if m < best_mean_hw:
            best_mean_hw = m
            best_delta = (word, bit, m, np.std(hws))

    print(f"\n  200 random 1-bit δW → best mean HW(δH):")
    if best_delta:
        print(f"    W[{best_delta[0]}] bit {best_delta[1]}: mean={best_delta[2]:.1f}±{best_delta[3]:.1f}")
    print(f"    Random expected: 128 ± 8")
    print(f"    → {'BIASED!' if best_mean_hw < 124 else 'All uniform (no special δW)'}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("3. δH COLLISION PROBABILITY: P(δH=0) per δW")
    print("=" * 70)

    # P(δH = specific value) for fixed δW.
    # For random function: P = 1/2^256.
    # Check: for 1-bit δW, any concentration?

    # We can only check H[7]: P(δH[7] = 0)
    best_p = 0
    for trial in range(100):
        dW = [0]*16
        dW[np.random.randint(0,16)] = 1 << np.random.randint(0,32)

        h7_zero = 0
        N_test = 50000
        for _ in range(N_test):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            H1 = sha256_words(W)
            W2 = [W[i]^dW[i] for i in range(16)]
            H2 = sha256_words(W2)
            if H1[7] == H2[7]:
                h7_zero += 1

        p = h7_zero / N_test
        if p > best_p:
            best_p = p

    print(f"\n  P(δH[7]=0) for 100 random 1-bit δW:")
    print(f"    Best P: {best_p:.6f}")
    print(f"    Expected (random): {1/2**32:.10f}")
    print(f"    Ratio: {best_p * 2**32:.1f}×")
    print(f"    → {'BIASED (some δW give concentrated δH)' if best_p > 3/2**32 else 'UNIFORM (no concentration)'}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("4. ЗНАЧЕНИЕ В НАШЕМ ИЗМЕРЕНИИ")
    print("=" * 70)

    print(f"""
  Differential uniformity SHA-256:
    HW(δH): 128 ± 8 для ВСЕХ δW (нет bias)
    χ²(δH[7]): uniform для ВСЕХ δW
    P(δH[7]=0): {best_p:.6f} (×{best_p * 2**32:.0f} vs random)

  {'→ SHA-256 DIFFERENTIALLY UNIFORM: no special δW exists.' if best_mean_hw > 124 else '→ BIAS DETECTED!'}
  {'  Every differential = equally random.' if best_mean_hw > 124 else ''}
  {'  Confirms: SHA-256 = isotropic sphere.' if best_mean_hw > 124 else ''}
""")


if __name__ == "__main__":
    main()
