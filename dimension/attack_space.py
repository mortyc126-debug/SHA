"""
Построение пространства атаки на минимальное ядро SHA-256.

Ядро: CARRY + РОТАЦИИ. Всё остальное — избыточно.

Стратегия: разделяем SHA-256 на:
  - Линейную часть L: ротации, XOR, сдвиги, Ch/Maj (над GF(2))
  - Нелинейную часть N: carry от каждого ADD

В GF(2)-пространстве L — полностью вычислима.
N — это "ошибка", которую вносит carry.

Collision = найти δW такой, что L(δW) + N(δW) = 0.
Т.е. N(δW) = -L(δW) = L(δW) (в GF(2)).

Эксперимент:
  1. Для reduced SHA-256 (8 раундов) вычисляем L и N раздельно
  2. Измеряем: L предсказывает δH? N случайна?
  3. Ищем δW где N компенсирует L
  4. Масштабируем: 8 → 16 → 32 → 64 раунда
"""

import numpy as np

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
def shr(x, n): return x >> n
def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ shr(x, 3)
def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ shr(x, 10)
def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
def add32(x, y): return (x + y) & MASK32

def get_carry_word(x, y):
    """Carry-слово: бит i = carry out of position i."""
    c = 0
    carry = 0
    for i in range(32):
        s = ((x >> i) & 1) + ((y >> i) & 1) + c
        c = s >> 1
        carry |= (c << i)
    return carry


def sha256_decomposed(W_input, num_rounds=64):
    """
    SHA-256 с РАЗДЕЛЬНЫМ вычислением линейной и нелинейной частей.

    Для каждого ADD: result = x XOR y XOR carry_word(x,y)
    Линейная часть: x XOR y
    Нелинейная часть: carry_word(x,y)

    Возвращает:
      H_real: реальный хеш
      H_linear: хеш если бы carry=0 (XOR вместо ADD)
      carry_total: суммарное carry-слово (XOR всех промежуточных)
    """
    W = list(W_input[:16])
    for r in range(16, 64):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))

    # Реальное вычисление
    a, b, c, d, e, f, g, h = IV
    for r in range(num_rounds):
        T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e, f, g)), K[r]), W[r])
        T2 = add32(Sigma0(a), Maj(a, b, c))
        h, g, f, e = g, f, e, add32(d, T1)
        d, c, b, a = c, b, a, add32(T1, T2)
    H_real = [add32(IV[i], [a,b,c,d,e,f,g,h][i]) for i in range(8)]

    # Линейное вычисление (XOR вместо ADD)
    a, b, c, d, e, f, g, h = IV
    for r in range(num_rounds):
        T1 = h ^ Sigma1(e) ^ Ch(e, f, g) ^ K[r] ^ W[r]
        T2 = Sigma0(a) ^ Maj(a, b, c)
        h, g, f, e = g, f, e, (d ^ T1) & MASK32
        d, c, b, a = c, b, a, (T1 ^ T2) & MASK32
    H_linear = [IV[i] ^ [a,b,c,d,e,f,g,h][i] for i in range(8)]

    # Carry = разница между реальным и линейным
    carry_effect = [H_real[i] ^ H_linear[i] for i in range(8)]

    return H_real, H_linear, carry_effect


def experiment_attack_space():
    np.random.seed(42)

    print("=" * 70)
    print("ПРОСТРАНСТВО АТАКИ: L + N ДЕКОМПОЗИЦИЯ")
    print("=" * 70)

    # ===================================================================
    # ЭТАП 1: Декомпозиция δH = δL ⊕ δN
    # ===================================================================
    print("\n" + "=" * 70)
    print("1. ДЕКОМПОЗИЦИЯ: δH = δL ⊕ δN")
    print("   δL = линейная часть (XOR), δN = carry-эффект")
    print("=" * 70)

    for num_rounds in [8, 16, 32, 64]:
        dl_hws = []
        dn_hws = []
        dh_hws = []
        dn_predicts_dh = 0

        for trial in range(2000):
            W1 = [np.random.randint(0, 2**32) for _ in range(16)]
            W2 = list(W1)
            W2[0] ^= (1 << np.random.randint(0, 32))

            H1_r, H1_l, C1 = sha256_decomposed(W1, num_rounds)
            H2_r, H2_l, C2 = sha256_decomposed(W2, num_rounds)

            # δ в каждом слое
            dH = [H1_r[i] ^ H2_r[i] for i in range(8)]
            dL = [H1_l[i] ^ H2_l[i] for i in range(8)]
            dN = [C1[i] ^ C2[i] for i in range(8)]

            hw_dH = sum(hw(d) for d in dH)
            hw_dL = sum(hw(d) for d in dL)
            hw_dN = sum(hw(d) for d in dN)

            dl_hws.append(hw_dL)
            dn_hws.append(hw_dN)
            dh_hws.append(hw_dH)

        print(f"\n  {num_rounds} раундов:")
        print(f"    HW(δH real):    {np.mean(dh_hws):6.1f} ± {np.std(dh_hws):.1f}")
        print(f"    HW(δL linear):  {np.mean(dl_hws):6.1f} ± {np.std(dl_hws):.1f}")
        print(f"    HW(δN carry):   {np.mean(dn_hws):6.1f} ± {np.std(dn_hws):.1f}")
        print(f"    Доля carry в δH: {np.mean(dn_hws)/max(np.mean(dh_hws),1)*100:.1f}%")

    # ===================================================================
    # ЭТАП 2: Корреляция δL и δN
    # ===================================================================
    print("\n" + "=" * 70)
    print("2. КОРРЕЛЯЦИЯ δL и δN: независимы или связаны?")
    print("=" * 70)

    for num_rounds in [8, 16, 32, 64]:
        dl_list = []
        dn_list = []

        for trial in range(2000):
            W1 = [np.random.randint(0, 2**32) for _ in range(16)]
            W2 = list(W1)
            W2[0] ^= (1 << np.random.randint(0, 32))

            H1_r, H1_l, C1 = sha256_decomposed(W1, num_rounds)
            H2_r, H2_l, C2 = sha256_decomposed(W2, num_rounds)

            dL = sum(hw(H1_l[i] ^ H2_l[i]) for i in range(8))
            dN = sum(hw(C1[i] ^ C2[i]) for i in range(8))

            dl_list.append(dL)
            dn_list.append(dN)

        corr = np.corrcoef(dl_list, dn_list)[0, 1]
        print(f"    {num_rounds:2d} раундов: corr(δL, δN) = {corr:+.4f} {'← СВЯЗАНЫ!' if abs(corr) > 0.1 else ''}")

    # ===================================================================
    # ЭТАП 3: Можно ли предсказать δN из δW?
    # ===================================================================
    print("\n" + "=" * 70)
    print("3. ПРЕДСКАЗУЕМОСТЬ δN: знаем δW → знаем δN?")
    print("=" * 70)

    # Для фиксированного W1, меняем 1 бит → δN зависит от ПОЗИЦИИ бита?
    for num_rounds in [8, 64]:
        W_base = [np.random.randint(0, 2**32) for _ in range(16)]
        H_base_r, H_base_l, C_base = sha256_decomposed(W_base, num_rounds)

        bit_to_dn = {}
        for bit in range(32):
            dns = []
            for trial in range(100):
                W_base = [np.random.randint(0, 2**32) for _ in range(16)]
                H_base_r, H_base_l, C_base = sha256_decomposed(W_base, num_rounds)

                W2 = list(W_base)
                W2[0] ^= (1 << bit)
                H2_r, H2_l, C2 = sha256_decomposed(W2, num_rounds)
                dn = sum(hw(C_base[i] ^ C2[i]) for i in range(8))
                dns.append(dn)
            bit_to_dn[bit] = (np.mean(dns), np.std(dns))

        print(f"\n  {num_rounds} раундов — HW(δN) по позиции бита δW₀:")
        for bit in [0, 1, 4, 8, 15, 16, 24, 31]:
            m, s = bit_to_dn[bit]
            bar = "█" * int(m / 4)
            print(f"    бит {bit:2d}: {m:5.1f} ± {s:.1f}  {bar}")

        # Разброс между битами
        means = [bit_to_dn[b][0] for b in range(32)]
        print(f"    Разброс mean(δN): {np.std(means):.2f}")
        print(f"    → {'ПОЗИЦИЯ БИТА ВЛИЯЕТ!' if np.std(means) > 2.0 else 'Все биты эквивалентны'}")

    # ===================================================================
    # ЭТАП 4: Поиск в carry-пространстве
    # ===================================================================
    print("\n" + "=" * 70)
    print("4. ПОИСК В ПРОСТРАНСТВЕ: δN = δL (условие коллизии)")
    print("=" * 70)

    # Collision: δH = δL ⊕ δN = 0, т.е. δL = δN
    # Ищем δW такой, что линейная и нелинейная части СОВПАДАЮТ
    # Мера близости: HW(δL ⊕ δN) = HW(δH) → минимизируем

    for num_rounds in [8, 16, 24, 32, 64]:
        best_hw = 256
        trials_done = 0
        hw_distribution = []

        for trial in range(5000):
            W1 = [np.random.randint(0, 2**32) for _ in range(16)]
            W2 = list(W1)
            W2[0] ^= np.random.randint(1, 2**32)

            H1_r, _, _ = sha256_decomposed(W1, num_rounds)
            H2_r, _, _ = sha256_decomposed(W2, num_rounds)

            dH_hw = sum(hw(H1_r[i] ^ H2_r[i]) for i in range(8))
            hw_distribution.append(dH_hw)

            if dH_hw < best_hw:
                best_hw = dH_hw
            trials_done += 1

        print(f"    {num_rounds:2d} раундов: best HW(δH) = {best_hw:3d}/256 "
              f"[mean={np.mean(hw_distribution):.1f}, "
              f"min теор.={max(0, int(np.mean(hw_distribution) - 4*np.std(hw_distribution)))}, "
              f"5000 проб]")

    # ===================================================================
    # ЭТАП 5: Направленный поиск — градиент в L-пространстве
    # ===================================================================
    print("\n" + "=" * 70)
    print("5. ГРАДИЕНТНЫЙ ПОИСК В L-ПРОСТРАНСТВЕ")
    print("=" * 70)

    # Идея: линейная часть L предсказуема. Если δL(δW) "указывает"
    # куда двигаться, можем ли мы выбрать δW так, чтобы δL был мал?
    # А потом надеяться что δN тоже мал?

    for num_rounds in [8, 16, 32, 64]:
        # Шаг 1: найти δW с минимальным HW(δL)
        best_dl = 256
        best_dw = 0
        best_dh_at_best_dl = 256

        for trial in range(5000):
            W1 = [np.random.randint(0, 2**32) for _ in range(16)]
            delta = np.random.randint(1, 2**32)
            W2 = list(W1)
            W2[0] ^= delta

            H1_r, H1_l, _ = sha256_decomposed(W1, num_rounds)
            H2_r, H2_l, _ = sha256_decomposed(W2, num_rounds)

            dl = sum(hw(H1_l[i] ^ H2_l[i]) for i in range(8))
            dh = sum(hw(H1_r[i] ^ H2_r[i]) for i in range(8))

            if dl < best_dl:
                best_dl = dl
                best_dw = delta
                best_dh_at_best_dl = dh

        print(f"    {num_rounds:2d} раундов: min HW(δL) = {best_dl:3d}, "
              f"при этом HW(δH real) = {best_dh_at_best_dl:3d}")

    # ===================================================================
    # ЭТАП 6: Корреляция размера δW с HW(δH)
    # ===================================================================
    print("\n" + "=" * 70)
    print("6. ОПТИМАЛЬНЫЙ ТИП δW")
    print("=" * 70)

    for num_rounds in [8, 64]:
        results_by_hw = {}
        for target_hw in [1, 2, 4, 8, 16]:
            dh_list = []
            for trial in range(1000):
                W1 = [np.random.randint(0, 2**32) for _ in range(16)]
                W2 = list(W1)
                # Генерим δW с заданным HW
                delta = 0
                bits = np.random.choice(32, target_hw, replace=False)
                for b in bits:
                    delta ^= (1 << b)
                W2[0] ^= delta

                H1, _, _ = sha256_decomposed(W1, num_rounds)
                H2, _, _ = sha256_decomposed(W2, num_rounds)
                dh = sum(hw(H1[i] ^ H2[i]) for i in range(8))
                dh_list.append(dh)

            results_by_hw[target_hw] = (np.mean(dh_list), np.min(dh_list))

        print(f"\n  {num_rounds} раундов:")
        for thw in [1, 2, 4, 8, 16]:
            m, mn = results_by_hw[thw]
            print(f"    HW(δW)={thw:2d}: mean HW(δH)={m:6.1f}, min={mn:3d}")

    # ===================================================================
    # ФИНАЛ
    # ===================================================================
    print("\n" + "=" * 70)
    print("ФИНАЛ: МОЖНО ЛИ ПОСТРОИТЬ ПРОСТРАНСТВО АТАКИ?")
    print("=" * 70)

    print(f"""
  ДЕКОМПОЗИЦИЯ δH = δL ⊕ δN:

    δL (линейная, XOR):  предсказуема, но HW ≈ 128 (полная диффузия)
    δN (carry-эффект):   непредсказуема, HW ≈ 128 (тоже полная)
    δH = δL ⊕ δN:        HW ≈ 128 (два случайных 128-бит XOR ≈ 128)

  ПРОБЛЕМА:
    Для collision нужно δL = δN (чтобы δH = 0).
    Но δL и δN — НЕКОРРЕЛИРОВАНЫ.
    Это два НЕЗАВИСИМЫХ случайных 256-бит вектора.
    P(δL = δN) = 2^{{-256}}.

  ПОЧЕМУ L-N ПРОСТРАНСТВО НЕ РАБОТАЕТ:
    Линейная часть L использует значения (a,b,c,d,e,f,g,h),
    которые САМИ зависят от carry.
    L(W) вычисляется как "XOR вместо ADD",
    но ПРОМЕЖУТОЧНЫЕ СОСТОЯНИЯ другие.
    Поэтому δL ≠ "линейная часть δH".
    δL — это δH ДРУГОЙ ФУНКЦИИ (XOR-SHA-256).

  КОРЕНЬ ПРОБЛЕМЫ:
    Carry меняет не только РЕЗУЛЬТАТ сложения,
    но и ВХОДЫ следующего сложения.
    Каждое carry изменяет state → изменяет входы ротаций →
    → изменяет входы следующего carry.

    CARRY И РОТАЦИИ НЕ РАЗДЕЛИМЫ.
    Они связаны через state. Декомпозиция L+N невозможна
    без полного знания промежуточного state.

  ВЫВОД:
    Пространство, разделяющее carry и ротации, НЕ СУЩЕСТВУЕТ.
    Это два компонента ОДНОГО неразделимого механизма.
    Их взаимодействие через state — и есть защита SHA-256.
""")


if __name__ == "__main__":
    experiment_attack_space()
