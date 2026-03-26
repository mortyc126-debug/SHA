"""
АЛГЕБРА МЕТОК: какая математическая структура у δ?

Метка δ = разность двух следов. Но какие ОПЕРАЦИИ определены на метках?
Замкнуты ли они? Есть ли нейтральный элемент? Обратный?

Определяем и проверяем:
  1. δ₁ ⊕ δ₂ = ? (overlay двух меток)
  2. δ · α = ? (масштабирование метки)
  3. Ассоциативность, коммутативность
  4. Нейтральный элемент (δ=0)
  5. Обратный элемент (-δ)
  6. ЗАМКНУТОСТЬ: δ₁ ⊕ δ₂ — тоже метка? (т.е. существует пара следов с такой разностью?)

А также:
  7. Ткань как ОПЕРАТОР на метках: T(δ) = δ_out
  8. Ядро оператора: ker(T) = {δ : T(δ) = 0} = collision set
  9. Образ оператора: im(T) = {T(δ)} = достижимые δ на выходе
"""

import numpy as np
import struct, hashlib

MASK32 = 0xFFFFFFFF

def sha256_words(W16):
    return struct.unpack('>8I', hashlib.sha256(struct.pack('>16I', *W16)).digest())

def hw(x): return bin(x).count('1')


def main():
    np.random.seed(42)

    print("=" * 70)
    print("АЛГЕБРА МЕТОК: структура пространства δ")
    print("=" * 70)

    # =================================================================
    print(f"\n{'=' * 70}")
    print("1. ОПЕРАЦИЯ НА МЕТКАХ: δ₁ ⊕ δ₂")
    print("=" * 70)

    # δ₁ = H(W_A) ⊕ H(W_B) — метка пары (A,B)
    # δ₂ = H(W_C) ⊕ H(W_D) — метка пары (C,D)
    # δ₁ ⊕ δ₂ = H(W_A) ⊕ H(W_B) ⊕ H(W_C) ⊕ H(W_D)

    # Замкнутость: δ₁ ⊕ δ₂ — это метка какой пары?
    # Если H(W_A) ⊕ H(W_B) ⊕ H(W_C) ⊕ H(W_D) = H(W_E) ⊕ H(W_F)?
    # Не обязательно. XOR четырёх хешей ≠ XOR двух.

    # НО: δ₁ ⊕ δ₂ = (H(A) ⊕ H(C)) ⊕ (H(B) ⊕ H(D))
    # = δ(A,C) ⊕ δ(B,D)
    # Это СУММА двух других меток. Замкнутость: XOR меток = XOR меток ✓

    # Формально: множество всех меток = {H(W₁) ⊕ H(W₂) : W₁,W₂ ∈ input}
    # Это ПОДПРОСТРАНСТВО GF(2)^256? Или всё GF(2)^256?

    N = 10000
    marks = set()
    for _ in range(N):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = [np.random.randint(0, 2**32) for _ in range(16)]
        H1 = sha256_words(W1)
        H2 = sha256_words(W2)
        delta = tuple(H1[i] ^ H2[i] for i in range(8))
        marks.add(delta)

    print(f"\n  {N} random пар → {len(marks)} уникальных меток")
    print(f"  Если бы всё GF(2)^256: все уникальны ✓" if len(marks) == N else
          f"  Коллизии меток: {N - len(marks)}")

    # Нулевая метка?
    zero_mark = tuple(0 for _ in range(8))
    has_zero = zero_mark in marks
    print(f"  Нулевая метка (collision): {'НАЙДЕНА!' if has_zero else 'не найдена'}")

    # Является ли множество меток ВСЕМ GF(2)^256?
    # Проверяем: можем ли сгенерировать ЛЮБУЮ заданную метку?
    target = tuple(np.random.randint(0, 2**32) for _ in range(8))
    found_target = target in marks
    print(f"  Случайная target метка в множестве: {'да' if found_target else 'нет (N слишком мал)'}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("2. СТРУКТУРА: (Marks, ⊕) — группа?")
    print("=" * 70)

    print(f"""
  Проверяем аксиомы группы для (Marks, ⊕):

  1. Замкнутость: δ₁ ⊕ δ₂ ∈ Marks?
     δ₁ ⊕ δ₂ = (H(A)⊕H(B)) ⊕ (H(C)⊕H(D))
     = H(A) ⊕ H(B) ⊕ H(C) ⊕ H(D)
     Это XOR 4 хешей. Не обязательно XOR 2 хешей.
     → НЕ ЗАМКНУТО в общем случае.

  2. Ассоциативность: (δ₁⊕δ₂)⊕δ₃ = δ₁⊕(δ₂⊕δ₃)?
     XOR ассоциативен → ДА.

  3. Нейтральный элемент: δ⊕0 = δ?
     0 = H(W)⊕H(W) (пара из одинаковых следов) → ДА.

  4. Обратный: δ⊕δ = 0?
     (H(A)⊕H(B)) ⊕ (H(A)⊕H(B)) = 0 → ДА. Каждая метка самообратна.

  Итого: (GF(2)^256, ⊕) — ГРУППА (Абелева).
  Marks ⊂ GF(2)^256 — ПОДМНОЖЕСТВО (не обязательно подгруппа).

  Ткань T: Marks → Marks? Нет!
  T: Inputs → Outputs. T(δ_input) = δ_output.
  T — ОПЕРАТОР на GF(2)^256 → GF(2)^256.
""")

    # =================================================================
    print(f"{'=' * 70}")
    print("3. ОПЕРАТОР ТКАНИ T: свойства")
    print("=" * 70)

    # T(δW) = δH = H(W₁) ⊕ H(W₁ ⊕ δW)
    # T: GF(2)^512 → GF(2)^256

    # Линейность: T(δW₁ ⊕ δW₂) = T(δW₁) ⊕ T(δW₂)?
    # Мы проверили (закон 3): НЕЛИНЕЙНО. diff ≈ 128.

    # НО: может T ПОЧТИ линеен? Для МАЛЫХ δW?

    print(f"\n  T(δW) = H(W) ⊕ H(W ⊕ δW)")
    print(f"  Линейность: T(a⊕b) = T(a) ⊕ T(b)?")

    # Малые δW: 1 бит
    linearity_diffs = []
    for _ in range(2000):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H_base = sha256_words(W)

        # δW₁ = 1 бит в W[0]
        bit1 = np.random.randint(0, 32)
        W_a = list(W); W_a[0] ^= (1 << bit1)
        H_a = sha256_words(W_a)
        T_a = tuple(H_base[i] ^ H_a[i] for i in range(8))

        # δW₂ = 1 бит в W[1]
        bit2 = np.random.randint(0, 32)
        W_b = list(W); W_b[1] ^= (1 << bit2)
        H_b = sha256_words(W_b)
        T_b = tuple(H_base[i] ^ H_b[i] for i in range(8))

        # δW₁ ⊕ δW₂
        W_ab = list(W); W_ab[0] ^= (1 << bit1); W_ab[1] ^= (1 << bit2)
        H_ab = sha256_words(W_ab)
        T_ab = tuple(H_base[i] ^ H_ab[i] for i in range(8))

        # T(a) ⊕ T(b)
        T_sum = tuple(T_a[i] ^ T_b[i] for i in range(8))

        # Linearity error
        error = sum(hw(T_ab[i] ^ T_sum[i]) for i in range(8))
        linearity_diffs.append(error)

    print(f"\n  |T(a⊕b) - T(a)⊕T(b)| для 1-бит δW:")
    print(f"    Mean: {np.mean(linearity_diffs):.1f}")
    print(f"    Min:  {min(linearity_diffs)}")
    print(f"    Max:  {max(linearity_diffs)}")
    print(f"    P(=0): {sum(1 for d in linearity_diffs if d == 0)/len(linearity_diffs)*100:.2f}%")

    if np.mean(linearity_diffs) < 10:
        print(f"    → ПОЧТИ ЛИНЕЙНО!")
    elif np.mean(linearity_diffs) < 64:
        print(f"    → ЧАСТИЧНО ЛИНЕЙНО")
    else:
        print(f"    → НЕЛИНЕЙНО")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("4. РАНГ ОПЕРАТОРА T (над GF(2))")
    print("=" * 70)

    # T: GF(2)^512 → GF(2)^256.
    # rank(T) = число независимых output bits.
    # Если rank < 256: image(T) < GF(2)^256 → не все метки достижимы.
    # Ker(T) > 0 → collisions exist.

    # Из T_GF2_BIJECTION (методичка): rank(J) = 256 для free-start.
    # Проверяем сами:

    # Генерируем 512 "базисных" δW (1 бит каждый)
    # и смотрим rank матрицы T.

    W_base = [np.random.randint(0, 2**32) for _ in range(16)]
    H_base = sha256_words(W_base)

    T_matrix = []  # 512 rows (input bits) × 256 columns (output bits)

    for word in range(16):
        for bit in range(32):
            W_mod = list(W_base)
            W_mod[word] ^= (1 << bit)
            H_mod = sha256_words(W_mod)

            row = []
            for w in range(8):
                delta = H_base[w] ^ H_mod[w]
                for b in range(32):
                    row.append((delta >> b) & 1)
            T_matrix.append(row)

    T_arr = np.array(T_matrix, dtype=np.float64)  # 512 × 256
    rank = np.linalg.matrix_rank(T_arr)

    print(f"\n  Матрица T: 512 входных бит → 256 выходных бит")
    print(f"  Rank(T) = {rank}")
    print(f"  Ker(T) dimension = 512 - {rank} = {512 - rank}")
    print(f"  → {'ПОЛНОРАНГОВЫЙ (rank=256)' if rank == 256 else f'НЕПОЛНОРАНГОВЫЙ (rank={rank})'}")

    if rank == 256:
        print(f"\n  Image(T) = GF(2)^256 (все метки достижимы)")
        print(f"  Ker(T) dimension = {512 - 256} = 256 (collision space)")
        print(f"  → Collision space 256-мерный над GF(2)")
        print(f"  → ~2^256 различных δW дают δH = 0 (в GF(2))")
        print(f"  → НО: GF(2) collision ≠ Z/2^32 collision (carry!)")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("5. GF(2)-COLLISION: решаем T·x = 0")
    print("=" * 70)

    # T_arr: 512×256. Ищем x ∈ {0,1}^512 с T·x = 0 (mod 2).
    # Ker(T) = 256-мерный. Базис ядра.

    # Transpose: T^T is 256×512. Solve T^T · x = 0.
    # Actually: T_matrix has rows = input bits, cols = output bits.
    # We need: for which input bit combinations, output = 0?
    # This is: null space of T_arr^T (as 256×512).

    T_gf2 = np.array(T_matrix, dtype=np.uint8)  # 512 × 256

    # Gaussian elimination over GF(2) on T_gf2^T (256 × 512)
    M = T_gf2.T.copy()  # 256 × 512
    pivots = []
    row = 0
    for col in range(512):
        found = False
        for r in range(row, 256):
            if M[r, col] == 1:
                M[[row, r]] = M[[r, row]]
                found = True
                break
        if not found:
            continue
        pivots.append(col)
        for r in range(256):
            if r != row and M[r, col] == 1:
                M[r] = M[r] ^ M[row]
        row += 1

    free_vars = [c for c in range(512) if c not in pivots]
    print(f"\n  GF(2) null space of T:")
    print(f"    Pivots: {len(pivots)}")
    print(f"    Free variables: {len(free_vars)}")
    print(f"    Kernel dimension: {len(free_vars)}")

    # Construct one kernel vector
    if free_vars:
        fc = free_vars[0]
        x = np.zeros(512, dtype=np.uint8)
        x[fc] = 1
        for i in range(len(pivots)-1, -1, -1):
            pc = pivots[i]
            val = 0
            for j in range(512):
                if j != pc:
                    val ^= (M[i, j] & x[j])
            x[pc] = val

        # Verify: T_gf2^T · x should be 0
        check = T_gf2.T @ x % 2
        is_valid = np.sum(check) == 0

        # Convert x to δW (16 words)
        delta_W = [0] * 16
        for word in range(16):
            for bit in range(32):
                if x[word * 32 + bit]:
                    delta_W[word] ^= (1 << bit)

        support = sum(1 for w in delta_W if w != 0)
        hw_total = sum(hw(w) for w in delta_W)

        print(f"\n  GF(2) kernel vector:")
        print(f"    Valid: {is_valid}")
        print(f"    Support: {support}/16 words nonzero")
        print(f"    Total HW: {hw_total}")

        # Apply to real SHA-256
        W_test = list(W_base)
        W_test2 = [W_test[i] ^ delta_W[i] for i in range(16)]
        H1 = sha256_words(W_test)
        H2 = sha256_words(W_test2)
        real_dh = sum(hw(H1[i] ^ H2[i]) for i in range(8))

        print(f"\n  Реальный SHA-256 с GF(2)-kernel δW:")
        print(f"    HW(δH) = {real_dh}")
        print(f"    {'GF(2) COLLISION WORKS!' if real_dh == 0 else f'GF(2)≠Z/2^32: carry adds {real_dh} бит'}")

        # Try multiple kernel vectors
        print(f"\n  10 kernel vectors → HW(δH):")
        real_dhs = []
        for fc in free_vars[:10]:
            x = np.zeros(512, dtype=np.uint8)
            x[fc] = 1
            for i in range(len(pivots)-1, -1, -1):
                pc = pivots[i]
                val = 0
                for j in range(512):
                    if j != pc:
                        val ^= (M[i, j] & x[j])
                x[pc] = val

            dW = [0] * 16
            for word in range(16):
                for bit in range(32):
                    if x[word * 32 + bit]:
                        dW[word] ^= (1 << bit)

            W2 = [W_base[i] ^ dW[i] for i in range(16)]
            H1 = sha256_words(W_base)
            H2 = sha256_words(W2)
            dh = sum(hw(H1[i] ^ H2[i]) for i in range(8))
            real_dhs.append(dh)
            print(f"    kernel[{fc}]: HW(δH) = {dh}")

        print(f"\n  Mean HW(δH) from GF(2) kernel: {np.mean(real_dhs):.1f}")
        print(f"  Min: {min(real_dhs)}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("6. ЗНАЧЕНИЕ: GF(2)-kernel в нашем измерении")
    print("=" * 70)

    print(f"""
  GF(2)-kernel dimension = {len(free_vars)}
  Каждый kernel vector = δW с T(δW) = 0 над GF(2).
  В реальной SHA-256 (Z/2^32): carry добавляет ≈{np.mean(real_dhs):.0f} бит.

  В нашем измерении:
    GF(2)-collision = ЛИНЕЙНОЕ ядро ткани.
    Z/2^32-collision = нелинейное ядро (carry).
    Расстояние: {np.mean(real_dhs):.0f} бит (carry noise).

  Carry noise = ЕДИНСТВЕННОЕ что отделяет GF(2)-collision от real collision.
  Если бы carry = 0: collision за O(256³) = 2^24 (Gaussian elimination).
  С carry: GF(2)-kernel + carry compensation = ???

  Стоимость carry compensation:
    Mean carry error = {np.mean(real_dhs):.0f} бит.
    Нужно найти kernel vector с carry error = 0.
    Среди 2^{len(free_vars)} kernel vectors: P(carry=0) ≈ ???
""")


if __name__ == "__main__":
    main()
