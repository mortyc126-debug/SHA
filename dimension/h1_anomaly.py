"""
H[1] ANOMALY: λ/√N = 0.05 — циклы 20× короче. Масштабируется?

H[1] = b[64] + IV[1] = a[63] + IV[1] (pipe a→b).
Почему H[1] даёт аномально короткие циклы?

Тестируем на 2^8, 2^12, 2^16, 2^20, 2^24.
"""

import numpy as np
import struct, hashlib

MASK32 = 0xFFFFFFFF

def sha256_h_pos(x, n_bits, position):
    mask = (1 << n_bits) - 1
    W = [x & mask] + [0] * 15
    h = hashlib.sha256(struct.pack('>16I', *W)).digest()
    words = struct.unpack('>8I', h)
    return words[position] & mask


def floyd(x0, n_bits, position):
    f = lambda x: sha256_h_pos(x, n_bits, position)
    t = f(x0)
    h = f(f(x0))
    steps = 1
    while t != h:
        t = f(t)
        h = f(f(h))
        steps += 1
        if steps > 5 * (2**n_bits):
            return None, None
    t = x0
    mu = 0
    while t != h:
        t = f(t)
        h = f(h)
        mu += 1
    lam = 1
    h = f(t)
    while t != h:
        h = f(h)
        lam += 1
    return mu, lam


def main():
    np.random.seed(42)

    print("=" * 70)
    print("H[1] ANOMALY VERIFICATION")
    print("=" * 70)

    positions = [0, 1, 3, 4, 7]  # NODE, PIPE-anomaly, PIPE, NODE, PIPE
    types = {0: 'NODE', 1: 'PIPE★', 3: 'PIPE', 4: 'NODE', 7: 'PIPE'}

    for n_bits in [8, 12, 16, 20, 24]:
        N = 2**n_bits
        expected = np.sqrt(np.pi * N / 8)
        n_trials = min(20, max(3, 2**(min(n_bits, 16) // 3)))

        print(f"\n  2^{n_bits} (N={N:,}, √N≈{int(expected)}, {n_trials} trials):")

        for pos in positions:
            mus = []
            lams = []

            for trial in range(n_trials):
                x0 = np.random.randint(0, 2**n_bits)
                mu, lam = floyd(x0, n_bits, pos)
                if mu is not None:
                    mus.append(mu)
                    lams.append(lam)

            if lams:
                r_lam = np.mean(lams) / expected
                marker = " ★★★" if r_lam < 0.1 else " ★★" if r_lam < 0.3 else " ★" if r_lam < 0.7 else ""
                print(f"    H[{pos}] ({types[pos]:>5}): λ={np.mean(lams):>8.0f}  λ/√N={r_lam:.3f}{marker}")
            else:
                print(f"    H[{pos}] ({types[pos]:>5}): TIMEOUT")

    print(f"\n{'=' * 70}")
    print("ANALYSIS")
    print("=" * 70)

    print(f"""
  Если H[1] anomaly МАСШТАБИРУЕТСЯ:
    λ/√N ≈ 0.05 для все размеры
    ρ-cost = μ + λ ≈ √N × (ratio_μ + ratio_λ)
    При ratio_λ = 0.05: ρ-cost ≈ √N × 1.05 ≈ √N

  Hmm — ρ-cost = μ + λ. Если λ мал но μ нормальный:
    Total cost ≈ μ ≈ √N. Короткий цикл не помогает с ХВОСТОМ.

  Для collision: нужен вход в ЦИКЛ (μ шагов) + цикл (λ шагов).
  Collision точка = на стыке хвоста и цикла.
  Стоимость = μ ≈ √N независимо от λ.

  Короткий λ полезен для:
    - Определения collision (нужен λ чтобы найти cycle entry)
    - ПАМЯТИ: λ определяет размер хранимого
    - НО НЕ ДЛЯ СКОРОСТИ: μ доминирует

  ИСКЛЮЧЕНИЕ: если μ ТОЖЕ короче.
""")


if __name__ == "__main__":
    main()
