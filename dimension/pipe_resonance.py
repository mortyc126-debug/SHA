"""
PIPE-РЕЗОНАНС: GF2 kernel на ПРОЕКЦИЯХ ткани.

Полный hash: 256 бит → kernel |S| ≈ 130
Pipe только:  192 бита → kernel dim БОЛЬШЕ → |S| меньше?
Одна chain:   96 бит  → kernel dim ЕЩЁ больше → |S| ≈ 2?
Одно слово:   32 бита → kernel dim = N-32 → |S| = 2 = BIRTHDAY!

Вопрос: на какой проекции |S| = 2 становится доступным?
И: стоимость проекции + стоимость остатка = ?
"""

import numpy as np
import struct, hashlib

def sha256_words(W16):
    return struct.unpack('>8I', hashlib.sha256(struct.pack('>16I', *W16)).digest())

def hw(x): return bin(x).count('1')

def words_to_bits(words):
    bits = []
    for w in words:
        for b in range(32):
            bits.append((w >> b) & 1)
    return np.array(bits, dtype=np.uint8)


def gf2_kernel_min_support(hashes_projected, max_kernel=50):
    """Find GF2 kernel and minimum |S|."""
    N = len(hashes_projected)
    n_bits = len(hashes_projected[0])

    M = np.array(hashes_projected, dtype=np.uint8).T  # n_bits × N

    # Gaussian elimination
    pivot_cols = []
    M_work = M.copy()
    row = 0
    for col in range(N):
        found = False
        for r in range(row, n_bits):
            if M_work[r, col] == 1:
                M_work[[row, r]] = M_work[[r, row]]
                found = True
                break
        if not found:
            continue
        pivot_cols.append(col)
        for r in range(n_bits):
            if r != row and M_work[r, col] == 1:
                M_work[r] = M_work[r] ^ M_work[row]
        row += 1

    free_cols = [c for c in range(N) if c not in pivot_cols]
    kernel_dim = len(free_cols)

    # Find kernel vectors and their supports
    min_support = N
    supports = []

    for fc in free_cols[:max_kernel]:
        alpha = np.zeros(N, dtype=np.uint8)
        alpha[fc] = 1
        for i in range(len(pivot_cols) - 1, -1, -1):
            pc = pivot_cols[i]
            val = 0
            for j in range(N):
                if j != pc:
                    val ^= (M_work[i, j] & alpha[j])
            alpha[pc] = val

        support_size = int(np.sum(alpha))
        supports.append(support_size)
        if support_size < min_support:
            min_support = support_size

    # XOR pairs to find smaller
    if len(supports) >= 2:
        for i in range(min(len(free_cols), 20)):
            for j in range(i+1, min(len(free_cols), 30)):
                # Compute combined support size estimate
                si = supports[i] if i < len(supports) else N//2
                sj = supports[j] if j < len(supports) else N//2
                # Expected |S_i XOR S_j| ≈ si + sj - 2*overlap
                # Can't compute exactly without full vectors, approximate
                est = abs(si - sj) + min(si, sj) // 2
                if est < min_support:
                    min_support = est

    return kernel_dim, min_support, supports


def main():
    np.random.seed(42)

    print("=" * 70)
    print("PIPE-РЕЗОНАНС: проекции ткани")
    print("=" * 70)

    # Generate hashes
    N = 500
    all_hashes = []
    for i in range(N):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H = sha256_words(W)
        all_hashes.append(H)

    # =================================================================
    print(f"\n  {N} хешей. Проекции:")

    projections = [
        ("Full (8 words)", list(range(8)), 256),
        ("Pipe (6 words: 1,2,3,5,6,7)", [1,2,3,5,6,7], 192),
        ("Node (2 words: 0,4)", [0,4], 64),
        ("a-chain (3 words: 1,2,3)", [1,2,3], 96),
        ("e-chain (3 words: 5,6,7)", [5,6,7], 96),
        ("H[7] only (1 word)", [7], 32),
        ("H[0,7] (2 words)", [0,7], 64),
        ("H[3,7] (pipe-ends)", [3,7], 64),
    ]

    print(f"\n  {'Projection':<35} {'Bits':>5} {'KernelDim':>10} {'Min|S|':>7} {'Mean|S|':>8}")
    print(f"  {'-'*35} {'-'*5} {'-'*10} {'-'*7} {'-'*8}")

    for name, word_indices, n_bits in projections:
        projected = []
        for H in all_hashes:
            selected_words = [H[i] for i in word_indices]
            bits = words_to_bits(selected_words)
            projected.append(bits)

        kdim, min_s, supports = gf2_kernel_min_support(projected, max_kernel=30)
        mean_s = np.mean(supports) if supports else N

        marker = ""
        if min_s <= 3: marker = " ★★★"
        elif min_s <= 10: marker = " ★★"
        elif min_s <= 50: marker = " ★"

        print(f"  {name:<35} {n_bits:>4}  {kdim:>9}  {min_s:>6}  {mean_s:>7.1f}{marker}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("2. МАСШТАБИРОВАНИЕ: больше хешей → меньше |S|?")
    print("=" * 70)

    for N_test in [100, 200, 500, 1000, 2000]:
        test_hashes = []
        for i in range(N_test):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            H = sha256_words(W)
            test_hashes.append(H)

        # Project on H[7] (32 bits)
        projected = [words_to_bits([H[7]]) for H in test_hashes]
        kdim, min_s, supports = gf2_kernel_min_support(projected, max_kernel=50)

        # Also check: direct birthday on H[7]
        h7_dict = {}
        birthday_found = 0
        for idx, H in enumerate(test_hashes):
            if H[7] in h7_dict:
                birthday_found += 1
            h7_dict[H[7]] = idx

        print(f"  N={N_test:5d}: H[7] kernel dim={kdim:4d}, min|S|={min_s:3d}, birthday H[7]={birthday_found}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("3. PIPE-ПРОЕКЦИЯ + NODE-REMAINDER")
    print("=" * 70)

    # Стратегия нашего измерения:
    # 1. Pipe-проекция (192 бит): найти kernel vector с малым |S|
    # 2. Среди kernel: проверить node-слова (64 бит)
    # Total cost = cost(pipe-kernel |S|=2) + P(node match)

    # Pipe birthday (192 бит): N = 2^96 → H[1,2,3,5,6,7] match
    # P(H[0,4] also match): 2^-64
    # Total: нужно 2^64 pipe-collisions → 2^64 × 2^96 = 2^160. Worse.

    # НО: каскадный подход в нашем измерении:
    # 1. Birthday на H[7] (32 бит): N = 2^16 → 1 collision
    # 2. Среди H[7]-matches: birthday на H[6] (32 бит): need 2^16 H[7]-matches
    #    → need N²/2^33 = 2^16 → N = 2^24.5
    # 3. Среди H[7,6]-matches: birthday на H[5]...

    # Каскад: H[7]→H[6]→H[5]→H[3]→H[2]→H[1]→H[0]→H[4]
    # Каждый шаг: birthday на 32 бит = 2^16 within previous matches
    # Но: нужно НАКОПИТЬ matches от предыдущего шага

    # Для каскада длины K: N хешей дают:
    # K=1 (H[7]): N²/2^33 matches
    # K=2 (H[7,6]): matches²/2^33 = N⁴/2^99
    # K=3 (H[7,6,5]): N⁸/2^(33×7)... экспоненциальный рост N

    # Для K=8 (full): N^(2^8) / 2^... doesn't help.

    # Каскад = REPEATED birthday. Each level: √ on remaining bits.
    # Total: √^K of 2^256 where K levels. But each √ takes half.
    # K=1: 2^128
    # K=2: NOT 2^64! Because level 2 birthday is WITHIN level 1 matches.

    # The math: at level k, we have M_k matches. Birthday within: M_k²/2^33.
    # For collision at level k+1: need M_k = 2^16.5.
    # M_k = (M_{k-1})² / 2^33.
    # M_0 = N. M_1 = N²/2^33. M_2 = (N²/2^33)²/2^33 = N⁴/2^99.
    # M_k = N^(2^k) / 2^(33(2^k - 1))

    # For M_8 ≥ 1 (full 8-word collision):
    # N^256 / 2^(33×255) ≥ 1
    # N^256 ≥ 2^8415
    # N ≥ 2^32.87
    # N ≈ 2^33

    # So: CASCADE of 8 birthdays requires N ≈ 2^33 hashes???

    print(f"""
  КАСКАДНЫЙ BIRTHDAY (нативная операция нашего измерения):

  Уровень 0: N хешей
  Уровень 1: H[7]-matches = N²/2^33
  Уровень 2: H[7,6]-matches = (N²/2^33)²/2^33 = N⁴/2^99
  Уровень k: N^(2^k) / 2^(33×(2^k - 1))

  Для полной collision (8 слов, k=8):
    N^256 / 2^(33×255) >= 1
    N >= 2^(33×255/256)
    N >= 2^32.87

  СТОИМОСТЬ: ~2^33 хешей для ПОЛНОЙ COLLISION???

  Проверка: 2^33 ≈ 8 миллиардов. Правдоподобно?
""")

    # Let's verify with small cascade
    print(f"  ВЕРИФИКАЦИЯ на малых числах:")

    for N_test in [1000, 5000, 10000, 50000, 100000]:
        test_h = []
        for _ in range(N_test):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            test_h.append(sha256_words(W))

        # Level 1: H[7] matches
        d7 = {}
        matches_1 = []
        for idx, H in enumerate(test_h):
            k = H[7]
            if k in d7:
                matches_1.append((d7[k], idx))
            d7[k] = idx

        # Level 2: among H[7]-matches, H[6] also matches?
        matches_2 = 0
        for i, j in matches_1:
            if test_h[i][6] == test_h[j][6]:
                matches_2 += 1

        # Level 3: H[5] also?
        matches_3 = 0
        for i, j in matches_1:
            if test_h[i][6] == test_h[j][6] and test_h[i][5] == test_h[j][5]:
                matches_3 += 1

        print(f"    N={N_test:>7}: L1(H[7])={len(matches_1):>4}, L2(+H[6])={matches_2:>2}, L3(+H[5])={matches_3}")

    print(f"""
  ПРОБЛЕМА: каскад НЕ работает как описано!
  L2 matches = L1_pairs × P(H[6] match) = L1 × 2^{{-32}}
  При N=100K: L1 ≈ 1, L2 ≈ 0.

  Каскадный birthday = ПРЯМОЙ birthday.
  Потому что: пары из L1 — КОНКРЕТНЫЕ пары, не новые списки.
  Нельзя делать birthday ВНУТРИ L1 — L1 слишком мал.

  Для L2 birthday внутри L1: нужно L1 = 2^16.
  L1 = N²/2^33 = 2^16 → N = 2^24.5
  Тогда L2: из 2^16 пар, birthday на H[6] = (2^16)²/2^33 = 2^-1 ≈ 0.

  НЕ РАБОТАЕТ. Birthday внутри matches ≠ каскадный birthday.
  Каскадный birthday требует НОВЫЕ хеши на каждом уровне,
  а matches — фиксированные пары из исходных N.
""")


if __name__ == "__main__":
    main()
