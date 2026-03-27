"""
НОВОЕ НАПРАВЛЕНИЕ 2: Фазовый переход.

В физике: вода → лёд при 0°C. Резкий переход.
SHA-256: "структура" → "хаос" при каком-то раунде r_c.

Мы знаем:
  r < 16: dead positions (не все слова использованы)
  r = 16-20: pre-sphere (неполное поглощение)
  r > 20: sphere (random)

Но это ГРАДИЕНТ или ФАЗОВЫЙ ПЕРЕХОД?
Есть ли КРИТИЧЕСКАЯ ТОЧКА r_c где поведение РЕЗКО меняется?

Аналогия: temperature parameter = round number.
Order parameter = какая-то мера структуры.
При r < r_c: order parameter > 0 (структура).
При r > r_c: order parameter = 0 (хаос).

Ищем: order parameter и critical point.
"""

import numpy as np
import struct, hashlib

MASK32 = 0xFFFFFFFF

def hw(x): return bin(x).count('1')
def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def add32(x, y): return (x + y) & MASK32
def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)

K = [
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
]
IV = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]


def sha256_r_rounds(W16, n_r):
    W = list(W16)
    for r in range(16, max(n_r, 16)):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))
    a,b,c,d,e,f,g,h = IV
    for r in range(n_r):
        T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K[r]), W[r])
        T2 = add32(Sigma0(a), Maj(a,b,c))
        h,g,f,e = g,f,e,add32(d,T1)
        d,c,b,a = c,b,a,add32(T1,T2)
    return tuple(add32(IV[i], [a,b,c,d,e,f,g,h][i]) for i in range(8))


def main():
    np.random.seed(42)

    print("=" * 70)
    print("ФАЗОВЫЙ ПЕРЕХОД: критическая точка SHA-256")
    print("=" * 70)

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("A. Order Parameter 1: SENSITIVITY (Lyapunov exponent)")
    print("=" * 70)

    # Sensitivity = HW(δH) for 1-bit flip, normalized to [0,1]
    # 0 = no sensitivity (ordered), 0.5 = full avalanche (chaotic)
    print(f"\n  {'Round':>5} {'Sensitivity':>12} {'Normalized':>11}")

    for n_r in range(1, 65):
        sensitivities = []
        for _ in range(200):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            bit = np.random.randint(0, 32)
            word = np.random.randint(0, 16)
            W2 = list(W); W2[word] ^= (1 << bit)
            H1 = sha256_r_rounds(W, n_r)
            H2 = sha256_r_rounds(W2, n_r)
            d = sum(hw(H1[i]^H2[i]) for i in range(8))
            sensitivities.append(d)
        mean_s = np.mean(sensitivities)
        norm = mean_s / 256.0  # normalize to [0, 1]
        if n_r <= 10 or n_r % 5 == 0 or n_r >= 60:
            bar = "#" * int(norm * 50)
            print(f"  {n_r:>5} {mean_s:>11.1f} {norm:>10.4f}  {bar}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("B. Order Parameter 2: DISTINGUISHABILITY from random")
    print("=" * 70)

    # For each round: how distinguishable is the output from random?
    # Measure: chi-squared on first byte of H[0]

    print(f"\n  {'Round':>5} {'chi2':>8} {'Expected':>9} {'Distinguishable':>16}")

    for n_r in range(1, 65):
        byte_counts = {}
        N = 5000
        for _ in range(N):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            H = sha256_r_rounds(W, n_r)
            byte_val = H[0] & 0xFF
            byte_counts[byte_val] = byte_counts.get(byte_val, 0) + 1

        expected = N / 256
        chi2 = sum((byte_counts.get(b, 0) - expected)**2 / expected for b in range(256))
        is_dist = chi2 > 300  # 255 ± 22, so >300 = 2σ above

        if n_r <= 10 or n_r % 5 == 0 or n_r >= 60:
            print(f"  {n_r:>5} {chi2:>7.0f} {255:>8} {'YES ★' if is_dist else 'no'}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("C. Order Parameter 3: BIT INDEPENDENCE")
    print("=" * 70)

    # Measure: correlation between output bits 0 and 1 of H[0]
    # For random function: 0. For structured: >0.

    print(f"\n  {'Round':>5} {'Bit corr':>9} {'Phase':>8}")

    transition_round = None
    for n_r in range(1, 65):
        same = 0
        N = 5000
        for _ in range(N):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            H = sha256_r_rounds(W, n_r)
            b0 = H[0] & 1
            b1 = (H[0] >> 1) & 1
            if b0 == b1: same += 1
        corr = abs(same/N - 0.5) * 2  # normalize: 0 = independent, 1 = fully correlated

        phase = "ORDERED" if corr > 0.05 else "CHAOTIC"
        if n_r <= 10 or n_r % 5 == 0 or n_r >= 60 or (transition_round is None and corr < 0.05):
            print(f"  {n_r:>5} {corr:>8.4f} {phase:>8}")

        if transition_round is None and corr < 0.02:
            transition_round = n_r

    if transition_round:
        print(f"\n  ★ Transition to CHAOTIC at round {transition_round}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("D. DERIVATIVE: where is sensitivity changing FASTEST?")
    print("=" * 70)

    # Compute sensitivity at every round, then take derivative
    sensitivities_full = []
    for n_r in range(1, 65):
        vals = []
        for _ in range(300):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            W2 = list(W); W2[0] ^= 1
            H1 = sha256_r_rounds(W, n_r)
            H2 = sha256_r_rounds(W2, n_r)
            vals.append(sum(hw(H1[i]^H2[i]) for i in range(8)))
        sensitivities_full.append(np.mean(vals))

    # Derivative
    derivs = [sensitivities_full[i+1] - sensitivities_full[i] for i in range(63)]

    print(f"\n  {'Round':>5} {'Sensitivity':>12} {'dS/dr':>8} {'|d2S/dr2|':>10}")

    max_deriv = 0
    max_deriv_r = 0
    for i in range(62):
        d1 = derivs[i]
        d2 = derivs[i+1] - derivs[i] if i < 62 else 0
        r = i + 1
        if abs(d1) > abs(max_deriv):
            max_deriv = d1
            max_deriv_r = r
        if r <= 10 or r % 10 == 0:
            print(f"  {r:>5} {sensitivities_full[i]:>11.1f} {d1:>7.1f} {abs(d2):>9.1f}")

    print(f"\n  Maximum |dS/dr| at round {max_deriv_r}: {max_deriv:.1f}")
    print(f"  → This is the INFLECTION POINT (fastest change)")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("E. CRITICAL EXPONENT")
    print("=" * 70)

    # Near the critical point, order parameter ~ |r - r_c|^β
    # Fit β from the data

    # Use sensitivity normalized to [0, 0.5]
    S_norm = [s / 256.0 for s in sensitivities_full]
    S_inf = 0.5  # saturation value

    # Find r_c where S crosses 0.25 (half of saturation)
    r_c = None
    for i in range(63):
        if S_norm[i] < 0.25 and S_norm[i+1] >= 0.25:
            # Linear interpolation
            r_c = i + 1 + (0.25 - S_norm[i]) / (S_norm[i+1] - S_norm[i])
            break

    if r_c:
        print(f"  Critical point r_c (S = 0.25): {r_c:.1f}")

        # Fit: S(r) = S_inf × (1 - exp(-α(r-r_0)))
        # Or: S(r) = S_inf × tanh(α(r - r_c))
        # Try: S(r) ≈ A × r^β for small r

        # Log-log fit for r < r_c
        valid = [(i+1, S_norm[i]) for i in range(int(r_c)) if S_norm[i] > 0.001]
        if len(valid) > 3:
            log_r = np.log([v[0] for v in valid])
            log_S = np.log([v[1] for v in valid])
            beta, log_A = np.polyfit(log_r, log_S, 1)
            A = np.exp(log_A)
            print(f"  Power law fit S ~ {A:.4f} × r^{beta:.2f}")
            print(f"  Critical exponent β = {beta:.2f}")
    else:
        print(f"  No clear critical point found")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("РЕЗУЛЬТАТ")
    print("=" * 70)


if __name__ == "__main__":
    main()
