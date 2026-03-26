"""
ГЛОБАЛЬНАЯ СТРАТЕГИЯ: все инструменты вместе.

Инструменты:
  1. R обратим (backward computation)
  2. T1 управляем через W
  3. Forward fix: δstate[1..16]=0
  4. Schedule линеен над GF(2)
  5. δe — максимальный усилитель (14.4 бит)
  6. Carry предсказуем при известном state

Стратегия А: Hill Climbing в δW-пространстве
  — 512 бит, flip по одному, минимизируем HW(δH)

Стратегия B: Partial Fix + Compatible Schedule
  — НЕ обнуляем δstate полностью, а направляем δ СОВМЕСТИМО с schedule

Стратегия C: Backward-Forward Pincer
  — Forward от IV, backward от target H, ищем meeting point

Стратегия D: GF(2)-Schedule Nullspace
  — Ищем δW где GF(2)-часть schedule diff МИНИМАЛЬНА
"""

import numpy as np
import time

MASK32 = 0xFFFFFFFF
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
IV = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]

def hw(x): return bin(x).count('1')
def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
def add32(x, y): return (x + y) & MASK32
def sub32(x, y): return (x - y) & MASK32
def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)

def sha256_compress(W16):
    W = list(W16)
    for r in range(16, 64):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))
    a, b, c, d, e, f, g, h = IV
    for r in range(64):
        T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e, f, g)), K[r]), W[r])
        T2 = add32(Sigma0(a), Maj(a, b, c))
        h, g, f, e = g, f, e, add32(d, T1)
        d, c, b, a = c, b, a, add32(T1, T2)
    return tuple(add32(IV[i], [a,b,c,d,e,f,g,h][i]) for i in range(8))

def hash_diff(W1, W2):
    H1 = sha256_compress(W1)
    H2 = sha256_compress(W2)
    return sum(hw(H1[i] ^ H2[i]) for i in range(8))


def strategy_A_hill_climbing(time_limit=30):
    """Hill climbing: flip бит в δW, принимаем если HW(δH) уменьшился."""
    print("\n" + "=" * 70)
    print("СТРАТЕГИЯ A: HILL CLIMBING (512-бит ландшафт)")
    print("=" * 70)

    np.random.seed(42)
    best_global = 256
    best_W1 = None
    best_W2 = None
    total_evals = 0
    improvements = []

    start = time.time()

    for restart in range(50):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        # Начальный δ: 1 случайный бит
        delta = [0] * 16
        delta[0] = 1 << np.random.randint(0, 32)
        W2 = [W1[i] ^ delta[i] for i in range(16)]

        current_hd = hash_diff(W1, W2)
        total_evals += 1

        # Hill climb: flip биты в delta
        no_improve = 0
        while no_improve < 200 and (time.time() - start) < time_limit:
            # Выбираем случайный бит для flip
            word_idx = np.random.randint(0, 16)
            bit_idx = np.random.randint(0, 32)

            delta[word_idx] ^= (1 << bit_idx)
            # Если delta стал 0 — пропускаем
            if all(d == 0 for d in delta):
                delta[word_idx] ^= (1 << bit_idx)
                continue

            W2_new = [W1[i] ^ delta[i] for i in range(16)]
            new_hd = hash_diff(W1, W2_new)
            total_evals += 1

            if new_hd < current_hd:
                W2 = W2_new
                current_hd = new_hd
                no_improve = 0
                if current_hd < best_global:
                    best_global = current_hd
                    best_W1 = list(W1)
                    best_W2 = list(W2)
                    improvements.append((total_evals, best_global))
            else:
                # Revert
                delta[word_idx] ^= (1 << bit_idx)
                no_improve += 1

        if (time.time() - start) >= time_limit:
            break

    elapsed = time.time() - start
    print(f"\n  Результат ({elapsed:.1f}с, {total_evals:,} evals, {len(improvements)} improvements):")
    print(f"    Best HW(δH) = {best_global}")
    print(f"    Прогресс: {' → '.join(str(x[1]) for x in improvements[:15])}")

    return best_global, total_evals


def strategy_B_partial_fix(time_limit=30):
    """Partial fix: не обнуляем δstate полностью, оставляем совместимый с schedule."""
    print("\n" + "=" * 70)
    print("СТРАТЕГИЯ B: PARTIAL FIX + COMPATIBLE SCHEDULE")
    print("=" * 70)

    np.random.seed(123)
    best_global = 256
    total_evals = 0
    start = time.time()

    # Идея: вместо δstate=0, ищем δstate который МИНИМИЗИРУЕТ
    # эффект на раундах 17-64.
    # Подход: пробуем разные "глубины" fix'а (0, 4, 8, 12, 16 раундов)
    # и для каждой — несколько δW

    for fix_depth in [0, 4, 8, 12, 16]:
        best_for_depth = 256
        trials = 0

        while trials < 5000 and (time.time() - start) < time_limit:
            W1 = [np.random.randint(0, 2**32) for _ in range(16)]
            W2 = list(W1)
            W2[0] ^= (1 << np.random.randint(0, 32))

            if fix_depth > 0:
                # Forward fix для первых fix_depth раундов
                def R_func(state, W_r, r_idx):
                    a, b, c, d, e, f, g, h = state
                    T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e, f, g)), K[r_idx]), W_r)
                    T2 = add32(Sigma0(a), Maj(a, b, c))
                    return (add32(T1, T2), a, b, c, add32(d, T1), e, f, g)

                # Compute target states from W1
                s1 = tuple(IV)
                states1 = [s1]
                for r in range(fix_depth):
                    s1 = R_func(s1, W1[r], r)
                    states1.append(s1)

                # Fix W2 to match states1
                s2 = tuple(IV)
                s2 = R_func(s2, W2[0], 0)
                for r in range(1, min(fix_depth, 16)):
                    a, b, c, d, e, f, g, h = s2
                    ap = states1[r + 1][0]
                    T2 = add32(Sigma0(a), Maj(a, b, c))
                    T1 = sub32(ap, T2)
                    W2[r] = sub32(sub32(sub32(sub32(T1, h), Sigma1(e)), Ch(e, f, g)), K[r])
                    s2 = R_func(s2, W2[r], r)

            hd = hash_diff(W1, W2)
            total_evals += 1
            trials += 1
            if hd < best_for_depth:
                best_for_depth = hd
            if hd < best_global:
                best_global = hd

        print(f"    fix_depth={fix_depth:2d}: best HW(δH) = {best_for_depth:3d} ({trials} trials)")

    print(f"\n  Best overall: {best_global}")
    return best_global, total_evals


def strategy_C_simulated_annealing(time_limit=30):
    """Simulated annealing: hill climbing с вероятностью принять ухудшение."""
    print("\n" + "=" * 70)
    print("СТРАТЕГИЯ C: SIMULATED ANNEALING")
    print("=" * 70)

    np.random.seed(456)
    best_global = 256
    best_W1 = None
    best_delta = None
    total_evals = 0
    start = time.time()

    for restart in range(20):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        delta = [0] * 16
        # Начальный δ: несколько бит
        for _ in range(np.random.randint(1, 5)):
            delta[np.random.randint(0, 16)] ^= (1 << np.random.randint(0, 32))
        if all(d == 0 for d in delta):
            delta[0] = 1

        W2 = [W1[i] ^ delta[i] for i in range(16)]
        current_hd = hash_diff(W1, W2)
        total_evals += 1
        temperature = 20.0

        for step in range(5000):
            if (time.time() - start) >= time_limit:
                break

            temperature *= 0.998
            word_idx = np.random.randint(0, 16)
            bit_idx = np.random.randint(0, 32)

            delta[word_idx] ^= (1 << bit_idx)
            if all(d == 0 for d in delta):
                delta[word_idx] ^= (1 << bit_idx)
                continue

            W2_new = [W1[i] ^ delta[i] for i in range(16)]
            new_hd = hash_diff(W1, W2_new)
            total_evals += 1

            diff = new_hd - current_hd
            if diff < 0 or (temperature > 0.1 and np.random.random() < np.exp(-diff / temperature)):
                W2 = W2_new
                current_hd = new_hd
                if current_hd < best_global:
                    best_global = current_hd
                    best_W1 = list(W1)
                    best_delta = list(delta)
            else:
                delta[word_idx] ^= (1 << bit_idx)

        if (time.time() - start) >= time_limit:
            break

    elapsed = time.time() - start
    print(f"\n  Результат ({elapsed:.1f}с, {total_evals:,} evals):")
    print(f"    Best HW(δH) = {best_global}")

    if best_delta:
        total_delta_hw = sum(hw(d) for d in best_delta)
        active_words = sum(1 for d in best_delta if d != 0)
        print(f"    HW(δW total) = {total_delta_hw}, active words = {active_words}")

    return best_global, total_evals


def strategy_D_multi_delta(time_limit=30):
    """Multi-word δ: вместо δ в одном слове, используем δ в нескольких."""
    print("\n" + "=" * 70)
    print("СТРАТЕГИЯ D: MULTI-WORD δ (оптимальная структура δW)")
    print("=" * 70)

    np.random.seed(789)
    best_global = 256
    total_evals = 0
    start = time.time()

    configs = [
        ("1 word, 1 bit", 1, 1),
        ("1 word, 4 bits", 1, 4),
        ("2 words, 1 bit each", 2, 1),
        ("4 words, 1 bit each", 4, 1),
        ("8 words, 1 bit each", 8, 1),
        ("16 words, 1 bit each", 16, 1),
        ("2 words, 8 bits each", 2, 8),
    ]

    for desc, n_words, n_bits in configs:
        best_for_config = 256
        trials = 0

        while trials < 10000 and (time.time() - start) < time_limit:
            W1 = [np.random.randint(0, 2**32) for _ in range(16)]
            W2 = list(W1)

            # Выбираем n_words случайных слов
            words = np.random.choice(16, n_words, replace=False)
            for w in words:
                bits = np.random.choice(32, n_bits, replace=False)
                for b in bits:
                    W2[w] ^= (1 << b)

            hd = hash_diff(W1, W2)
            total_evals += 1
            trials += 1
            if hd < best_for_config:
                best_for_config = hd
            if hd < best_global:
                best_global = hd

        print(f"    {desc:30s}: best = {best_for_config:3d} ({trials} trials)")

    print(f"\n  Best overall: {best_global}")
    return best_global, total_evals


def main():
    print("=" * 70)
    print("ГЛОБАЛЬНАЯ СТРАТЕГИЯ: ВСЕ ИНСТРУМЕНТЫ ВМЕСТЕ")
    print("Цель: минимизировать HW(δH)")
    print("=" * 70)

    time_per_strategy = 15  # секунд на стратегию

    results = {}

    r_a, e_a = strategy_A_hill_climbing(time_per_strategy)
    results['A: Hill Climbing'] = (r_a, e_a)

    r_b, e_b = strategy_B_partial_fix(time_per_strategy)
    results['B: Partial Fix'] = (r_b, e_b)

    r_c, e_c = strategy_C_simulated_annealing(time_per_strategy)
    results['C: Sim. Annealing'] = (r_c, e_c)

    r_d, e_d = strategy_D_multi_delta(time_per_strategy)
    results['D: Multi-δ'] = (r_d, e_d)

    # Random baseline
    print("\n" + "=" * 70)
    print("BASELINE: RANDOM SEARCH")
    print("=" * 70)
    np.random.seed(999)
    best_random = 256
    for trial in range(50000):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W1)
        W2[0] ^= (1 << np.random.randint(0, 32))
        hd = hash_diff(W1, W2)
        if hd < best_random:
            best_random = hd
    results['E: Random'] = (best_random, 50000)
    print(f"    Best HW(δH) = {best_random} (50K evals)")

    # === СВОДКА ===
    print("\n" + "=" * 70)
    print("СВОДКА ГЛОБАЛЬНОЙ СТРАТЕГИИ")
    print("=" * 70)

    print(f"\n  {'Стратегия':<25} {'Best HW(δH)':>12} {'Evals':>10} {'Efficiency':>12}")
    print(f"  {'-'*25} {'-'*12} {'-'*10} {'-'*12}")

    for name, (best, evals) in sorted(results.items(), key=lambda x: x[1][0]):
        eff = (128 - best) / max(np.log2(evals), 1)
        print(f"  {name:<25} {best:>11} {evals:>10,} {eff:>11.2f}")

    winner = min(results.items(), key=lambda x: x[1][0])
    print(f"\n  ЛУЧШАЯ СТРАТЕГИЯ: {winner[0]} → HW(δH) = {winner[1][0]}")
    print(f"  Birthday random (теор.): HW(δH) ≈ 0 при 2^128 evals")
    print(f"  Наш лучший результат:    HW(δH) = {winner[1][0]} при {winner[1][1]:,} evals")

    # Экстраполяция
    best_hw = winner[1][0]
    best_evals = winner[1][1]
    print(f"""
  ЭКСТРАПОЛЯЦИЯ:
    При {best_evals:,} evals: HW = {best_hw}
    Random при {best_evals:,}: HW ≈ {128 - int(np.sqrt(np.log2(best_evals)) * 8)}
    Gain over random: {128 - int(np.sqrt(np.log2(best_evals)) * 8) - best_hw:+d} бит

  ДЛЯ HW=0 (коллизия):
    Нужно ≈ 2^128 evals (birthday bound)
    Наша лучшая стратегия ≈ 2^{128 - 0:.0f} evals
    → Преимущество: ~0 бит
""")


if __name__ == "__main__":
    main()
