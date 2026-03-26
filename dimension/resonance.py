"""
РЕЗОНАНС (⊛): N-арная операция нашего измерения.

Идея: 257 хешей → Гауссово исключение над GF(2) → ядро.
Ядро содержит подмножество S с XOR = 0.
Если |S| = 2: collision!
Если |S| > 2: K-XOR-collision. Может конвертировать?

Также: в нашем измерении ткань имеет 16 ПОЗИЦИЙ входа.
Каждая позиция — независимый "голос".
Что если мы варьируем РАЗНЫЕ позиции НЕЗАВИСИМО?

Wagner's idea: K списков, K-XOR = 0.
  K = 16 позиций: cost = 2^(256/(1+4)) = 2^51 ???

Тестируем всё.
"""

import numpy as np
import struct, hashlib

def sha256_words(W16):
    return struct.unpack('>8I', hashlib.sha256(struct.pack('>16I', *W16)).digest())

def hw(x): return bin(x).count('1')

def hash_to_bits(H):
    """8 words → 256 bits as numpy array."""
    bits = np.zeros(256, dtype=np.uint8)
    for w in range(8):
        for b in range(32):
            bits[w * 32 + b] = (H[w] >> b) & 1
    return bits


def main():
    np.random.seed(42)

    print("=" * 70)
    print("РЕЗОНАНС: N-арная операция нашего измерения")
    print("=" * 70)

    # =================================================================
    print(f"\n{'=' * 70}")
    print("1. GF(2) ЯДРО: 257 хешей → линейная зависимость")
    print("=" * 70)

    N = 300  # > 256
    hashes = []
    W_list = []

    for i in range(N):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H = sha256_words(W)
        hashes.append(H)
        W_list.append(W)

    # Build matrix: 256 rows (bits) × N columns (hashes)
    M = np.zeros((256, N), dtype=np.uint8)
    for j in range(N):
        M[:, j] = hash_to_bits(hashes[j])

    # Gaussian elimination over GF(2) to find kernel
    # Kernel = vectors α where M·α = 0 (mod 2)

    # Augmented row reduction
    M_aug = M.copy()
    pivot_cols = []
    row = 0

    for col in range(N):
        # Find pivot
        found = False
        for r in range(row, 256):
            if M_aug[r, col] == 1:
                M_aug[[row, r]] = M_aug[[r, row]]
                found = True
                break
        if not found:
            continue  # free variable

        pivot_cols.append(col)

        # Eliminate
        for r in range(256):
            if r != row and M_aug[r, col] == 1:
                M_aug[r] = M_aug[r] ^ M_aug[row]
        row += 1

    free_cols = [c for c in range(N) if c not in pivot_cols]
    print(f"\n  {N} хешей, 256 бит каждый:")
    print(f"  Rank: {len(pivot_cols)}")
    print(f"  Free variables: {len(free_cols)}")
    print(f"  Kernel dimension: {len(free_cols)}")

    # Find kernel vectors
    kernel_vectors = []
    for fc in free_cols[:10]:  # max 10
        alpha = np.zeros(N, dtype=np.uint8)
        alpha[fc] = 1
        # Back-substitute
        for i in range(len(pivot_cols) - 1, -1, -1):
            pc = pivot_cols[i]
            val = 0
            for j in range(N):
                if j != pc:
                    val ^= (M_aug[i, j] & alpha[j])
            alpha[pc] = val

        # Verify
        check = np.zeros(256, dtype=np.uint8)
        for j in range(N):
            if alpha[j]:
                check = check ^ hash_to_bits(hashes[j])
        is_zero = np.sum(check) == 0

        support = [j for j in range(N) if alpha[j]]
        kernel_vectors.append((support, is_zero))

    print(f"\n  Kernel vectors (первые {len(kernel_vectors)}):")
    sizes = []
    for support, valid in kernel_vectors:
        sizes.append(len(support))
        if len(support) <= 10:
            print(f"    |S| = {len(support):3d}, valid = {valid}, indices = {support}")
        else:
            print(f"    |S| = {len(support):3d}, valid = {valid}")

    print(f"\n  Размеры подмножеств: min={min(sizes)}, max={max(sizes)}, mean={np.mean(sizes):.1f}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("2. ИЩЕМ МАЛЫЕ ПОДМНОЖЕСТВА (|S| = 2 = collision!)")
    print("=" * 70)

    # Все kernel vectors имеют |S| >> 2. Можем ли найти |S|=2?
    # |S|=2 → H[i] = H[j] → collision. P очень мала для 300 samples.

    # Но: можем комбинировать kernel vectors (XOR of kernel vectors = kernel)
    # Если v₁ и v₂ — kernel vectors, v₁ ⊕ v₂ тоже kernel vector.
    # |support(v₁ ⊕ v₂)| может быть МЕНЬШЕ!

    min_found = min(sizes)
    best_support = None

    # Try XOR of pairs of kernel vectors
    improvements = 0
    for i in range(len(kernel_vectors)):
        for j in range(i + 1, len(kernel_vectors)):
            si = set(kernel_vectors[i][0])
            sj = set(kernel_vectors[j][0])
            combined = si.symmetric_difference(sj)
            if len(combined) < min_found:
                min_found = len(combined)
                best_support = combined
                improvements += 1

    print(f"\n  XOR пар kernel vectors:")
    print(f"  Improvements: {improvements}")
    print(f"  Minimum |S| found: {min_found}")

    if min_found <= 5 and best_support:
        print(f"  Best support: {sorted(best_support)}")
        # Verify
        check = np.zeros(256, dtype=np.uint8)
        for j in best_support:
            check = check ^ hash_to_bits(hashes[j])
        print(f"  Valid (XOR=0): {np.sum(check) == 0}")

    # Triple XOR
    if min_found > 3:
        print(f"\n  Trying triple XOR...")
        for i in range(min(len(kernel_vectors), 5)):
            for j in range(i+1, min(len(kernel_vectors), 10)):
                for k in range(j+1, min(len(kernel_vectors), 15)):
                    si = set(kernel_vectors[i][0])
                    sj = set(kernel_vectors[j][0])
                    sk = set(kernel_vectors[k][0])
                    combined = si.symmetric_difference(sj).symmetric_difference(sk)
                    if len(combined) < min_found:
                        min_found = len(combined)
                        improvements += 1

        print(f"  After triple XOR: min |S| = {min_found}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("3. СТОИМОСТЬ РЕЗОНАНСА")
    print("=" * 70)

    print(f"""
  РЕЗОНАНС (GF2 kernel):
    N = 257 хешей → rank = 256 → kernel dim = {N - 256}
    Стоимость: 257 свёрток + O(256³) = O(2^24) операций

    Kernel vectors: |S| ≈ {np.mean(sizes):.0f} (не 2!)

    Для |S| = 2 (collision):
      Нужно МНОГО kernel vectors и XOR-комбинирование
      Или: N >> 257 хешей

    Оценка N для P(|S|=2):
      Kernel dim = N - 256.
      Число kernel vectors ≈ 2^{{N-256}}.
      Среди 2^{{N-256}} vectors: P(∃ |S|=2) ≈ ???

      Каждый vector: |S| ≈ N/2 (random).
      XOR двух: |S| ≈ N/2 (не уменьшается для random).
      Нужно |S|=2 из O(N) elements → NP-hard (subset sum-like).

  ВЕРДИКТ:
    Резонанс находит K-XOR-collision (K >> 2) за O(2^24).
    Но K-XOR → 2-collision конвертация = NP-hard.
    Резонанс ≠ collision shortcut.
""")

    # =================================================================
    print(f"{'=' * 70}")
    print("4. НАТИВНЫЙ ПОДХОД: ПОЗИЦИОННЫЙ РЕЗОНАНС")
    print("=" * 70)

    # Вместо random W: варьируем ОДНУ позицию за раз
    # 16 позиций × 2^32 значений = 16 списков по 2^32

    # Для collision: H(W_A) = H(W_B)
    # H зависит от 16 позиций НЕЛИНЕЙНО.

    # НО: в нашем измерении первые 3 позиции фиксированы,
    # break на W[3], a-repair на W[4..15].
    # Реальная свобода = W₁ (все 16 позиций).
    # W₂ ОПРЕДЕЛЕНА из W₁ через a-repair.

    # Collision = H(W₁) = H(a-repair(W₁)).
    # Это функция ОДНОГО аргумента W₁.
    # f(W₁) = H(W₁) ⊕ H(a-repair(W₁)) = δH(W₁)
    # Collision: f(W₁) = 0.

    # Это PREIMAGE задача: найти W₁ с f(W₁)=0.
    # f: 512 бит → 256 бит. Preimage: 2^256 (worse than birthday!)

    # НО: f(W₁) = δH не random! A-repair создаёт СТРУКТУРУ.
    # Мы показали: при a-repair, 58/64 раундов = 0.
    # f зависит от 5 amplification раундов.

    # Может f(W₁) имеет НИЗКИЙ АЛГЕБРАИЧЕСКИЙ DEGREE?

    print(f"""
  f(W1) = H(W1) XOR H(a-repair(W1)) = deltaH(W1)

  Collision: f(W1) = 0 (preimage, 256 бит)

  f зависит от:
    - W1[0..3]: fixed/break (не варьируем)
    - W1[4..15]: ОПРЕДЕЛЯЮТ a-repair → ОПРЕДЕЛЯЮТ W2
    - Amplification zone: 5 раундов

  В нашем измерении: f = свёртка 5 amplification слоёв.
  Если 5 слоёв имеют НИЗКИЙ РАНГ...
  Rank(5 amp layers) = ? Измеряем!
""")

    # Measure: что такое f(W1)?
    # Варьируем только W1[4] (1 позиция, 32 бит)
    # и смотрим как меняется δH.

    W_base = [np.random.randint(0, 2**32) for _ in range(16)]
    dh_by_w4 = []

    for _ in range(1000):
        W1 = list(W_base)
        W1[4] = np.random.randint(0, 2**32)
        H1 = sha256_words(W1)

        # A-repair для W1
        W1_saved = list(W1)
        H2 = sha256_words(W1_saved)  # same W → same H (need a-repair)
        # Simplified: just compare H(W1) for different W1[4]
        dh_by_w4.append(hash_to_bits(H1))

    # How many independent bits does varying W1[4] produce?
    mat = np.array(dh_by_w4[:256]).T  # 256 × 256
    rank_w4 = np.linalg.matrix_rank(mat.astype(float))

    print(f"  Varying W1[4] only (1000 samples):")
    print(f"  Rank of H output matrix: {rank_w4}")
    print(f"  → W1[4] controls {rank_w4} independent output bits")

    # All positions
    print(f"\n  Rank by position (varying one at a time):")
    for pos in range(16):
        samples = []
        for _ in range(260):
            W1 = list(W_base)
            W1[pos] = np.random.randint(0, 2**32)
            H1 = sha256_words(W1)
            samples.append(hash_to_bits(H1))

        mat = np.array(samples[:256]).T
        rk = np.linalg.matrix_rank(mat.astype(float))
        bar = "█" * (rk // 8)
        print(f"    W[{pos:2d}]: rank = {rk:3d}/256  {bar}")


if __name__ == "__main__":
    main()
