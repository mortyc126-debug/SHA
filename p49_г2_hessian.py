"""
П-49А: ГИПОТЕЗА Г2 — РАНГ ГЕССИАНА КВАДРАТИЧНОГО БАРЬЕРА

Вопрос: Почему mod-4 барьер абсолютен (0/393K решений)?

Гипотеза Г2: Матрица Гессиана (билинейная форма квадратичной поправки Q)
  над GF(2) имеет полный ранг 15.
  Если да → Q(x₀,δ) ≢ 0 для любого δ → лифтинг невозможен структурно.

Метод:
  1. Для каждого GF(2)-семени x₀ (F(x₀)≡0 mod 2):
     a. Вычислить матрицу Гессиана второй разности:
        B^(i)_{jk} = [F_i(x₀+2eⱼ+2eₖ) - F_i(x₀+2eⱼ) - F_i(x₀+2eₖ) + F_i(x₀)] mod 2
        (15 компонент × 105 пар j<k) → матрица 15×105 над GF(2)
     b. Ранг этой матрицы → если 15, то Q полноранговая (Г2 подтверждена).

  2. Полный перебор δ ∈ {0,1}^15 (32768):
     - Вычислить F(x₀+2δ) mod 4 для каждого δ
     - Найти распределение значений (image coverage)
     - Проверить: есть ли δ с F(x₀+2δ) ≡ 0 mod 4?

  3. Аналитика лифтинга:
     - Для δ Hensel (из решения J̃·δ = q̃ mod 2) проверить F(x₀+2δ) mod 4
     - Посмотреть на Q(x₀, δ_Hensel) — ненулевой ли он?
"""

import random

MASK = 0xFFFFFFFF

K_SHA = [
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
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
]
H0 = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]


def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK
def Sig0(x):    return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sig1(x):    return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def sig0(x):    return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sig1(x):    return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)


def _expand(W16):
    W = list(W16)
    for i in range(16, 18):
        W.append((sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16]) & MASK)
    return W


def sha256_states_17(W16):
    """Один проход 17 раундов, возвращает состояния после раундов 3..17.
    Ключи: r=3..17 → (a,b,c,d,e,f,g,h).  7× быстрее чем 15 отдельных вызовов."""
    W = _expand(W16)
    a, b, c, d, e, f, g, h = H0
    states = {}
    for r in range(17):
        T1 = (h + Sig1(e) + Ch(e, f, g) + K_SHA[r] + W[r]) & MASK
        T2 = (Sig0(a) + Maj(a, b, c)) & MASK
        h = g; g = f; f = e; e = (d + T1) & MASK
        d = c; c = b; b = a; a = (T1 + T2) & MASK
        if r >= 2:   # после раундов 0,1,2 → состояние r+1 = 3,4,...,17
            states[r + 1] = (a, b, c, d, e, f, g, h)
    return states   # states[3]..states[17]


# кэш базового блока — избегаем повторного вычисления sha256_states_17(W_base)
_cache_base = {}


def compute_f_fast(W_base, W2):
    """F = (Da₃..Da₁₆, De₁₇).  Принимает уже сложенный W2 (список 16 u32)."""
    key = tuple(W_base)
    if key not in _cache_base:
        _cache_base.clear()
        _cache_base[key] = sha256_states_17(W_base)
    s1 = _cache_base[key]
    s2 = sha256_states_17(W2)
    result = [(s2[r][0] - s1[r][0]) & MASK for r in range(3, 17)]   # Da
    result.append((s2[17][4] - s1[17][4]) & MASK)                    # De17
    return result


def make_W2(W_base, x16):
    return [(W_base[i] + x16[i]) & MASK for i in range(16)]


def rank_gf2(M):
    """Ранг матрицы M (список строк, каждая строка — список 0/1) над GF(2)."""
    if not M or not M[0]:
        return 0
    nrows, ncols = len(M), len(M[0])
    A = [[M[i][j] & 1 for j in range(ncols)] for i in range(nrows)]
    rank = 0
    pivot_col = 0
    for row in range(nrows):
        found = None
        for c in range(pivot_col, ncols):
            for r in range(row, nrows):
                if A[r][c] == 1:
                    found = (r, c)
                    break
            if found:
                break
        if not found:
            break
        r, c = found
        A[row], A[r] = A[r], A[row]
        pivot_col = c + 1
        for r2 in range(nrows):
            if r2 != row and A[r2][c] == 1:
                for k in range(ncols):
                    A[r2][k] ^= A[row][k]
        rank += 1
    return rank


def gauss_mod2_solve(J15, b):
    """Решить J·δ = b mod 2. Возвращает δ или None."""
    n = 15
    aug = [[J15[i][j] & 1 for j in range(n)] + [b[i] & 1] for i in range(n)]
    for col in range(n):
        pivot = None
        for row in range(col, n):
            if aug[row][col] == 1:
                pivot = row
                break
        if pivot is None:
            return None
        aug[col], aug[pivot] = aug[pivot], aug[col]
        for row in range(n):
            if row != col and aug[row][col] == 1:
                for k in range(n + 1):
                    aug[row][k] ^= aug[col][k]
    return [aug[i][n] for i in range(n)]


def find_seeds_mod2(W_base, dw0, n_trials=20000):
    """Найти ВСЕ GF(2)-семена x₀ ∈ {0,1}^15 с F(x₀)≡0 mod 2."""
    seeds = []
    seen = set()
    for _ in range(n_trials):
        bits = random.randrange(1 << 15)
        if bits in seen:
            continue
        seen.add(bits)
        x16 = [dw0] + [(bits >> j) & 1 for j in range(15)]
        W2 = make_W2(W_base, x16)
        f = compute_f_fast(W_base, W2)
        if all(v % 2 == 0 for v in f):
            seeds.append(x16)
    return seeds


# ═══════════════════════════════════════════════════════════════════════════════
# ЧАСТЬ A: МАТРИЦА ГЕССИАНА (ВТОРАЯ РАЗНОСТЬ)
# ═══════════════════════════════════════════════════════════════════════════════

def compute_hessian_matrix(W_base, x0_16):
    """
    Вычислить матрицу Гессиана второй разности для GF(2)-семени x₀.

    Возвращает: матрицу H размером 15×105 над GF(2), где
      строка i = компонента F_i,
      столбец (j,k) = B^(i)_{jk} = [F_i(x₀+2eⱼ+2eₖ) - F_i(x₀+2eⱼ) - F_i(x₀+2eₖ) + F_i(x₀)] mod 2
    """
    # Базовое значение F(x₀)
    W2_base = make_W2(W_base, x0_16)
    f0 = compute_f_fast(W_base, W2_base)

    # F(x₀ + 2eⱼ) для каждого j = 0..14  (соответствует DW_{j+1})
    f_ej = []
    for j in range(15):
        x_j = list(x0_16)
        x_j[j + 1] = (x_j[j + 1] + 2) & MASK
        W2_j = make_W2(W_base, x_j)
        f_ej.append(compute_f_fast(W_base, W2_j))

    # Матрица Гессиана: строки = компоненты F, столбцы = пары (j,k)
    H = []  # 15 строк
    for i in range(15):
        row = []
        for j in range(15):
            for k in range(j + 1, 15):
                # F(x₀ + 2eⱼ + 2eₖ)
                x_jk = list(x0_16)
                x_jk[j + 1] = (x_jk[j + 1] + 2) & MASK
                x_jk[k + 1] = (x_jk[k + 1] + 2) & MASK
                W2_jk = make_W2(W_base, x_jk)
                f_jk = compute_f_fast(W_base, W2_jk)

                # B = F(x₀+2eⱼ+2eₖ) - F(x₀+2eⱼ) - F(x₀+2eₖ) + F(x₀) mod 2
                B = (f_jk[i] - f_ej[j][i] - f_ej[k][i] + f0[i]) & MASK
                row.append(B & 1)
        H.append(row)

    return H, f0, f_ej


# ═══════════════════════════════════════════════════════════════════════════════
# ЧАСТЬ B: ПОЛНЫЙ ПЕРЕБОР δ ∈ {0,1}^15 — ОБРАЗ F(x₀ + 2δ) mod 4
# ═══════════════════════════════════════════════════════════════════════════════

def full_scan_mod4(W_base, x0_16, max_delta=32768):
    """
    Перебрать все δ ∈ {0,1}^15 (2^15 = 32768) или max_delta из них.
    Для каждого: вычислить F(x₀+2δ) mod 4.
    Вернуть: (список результатов, количество нулей mod 4, размер образа).
    """
    zeros_mod4 = 0
    image_set = set()
    zero_vectors = []

    limit = min(max_delta, 1 << 15)
    for bits in range(limit):
        delta = [(bits >> j) & 1 for j in range(15)]
        x_delta = list(x0_16)
        for j in range(15):
            x_delta[j + 1] = (x_delta[j + 1] + 2 * delta[j]) & MASK
        W2 = make_W2(W_base, x_delta)
        f = compute_f_fast(W_base, W2)

        # mod 4
        f_mod4 = tuple(v % 4 for v in f)
        image_set.add(f_mod4)

        if all(v == 0 for v in f_mod4):
            zeros_mod4 += 1
            zero_vectors.append(bits)

    return zeros_mod4, len(image_set), zero_vectors


# ═══════════════════════════════════════════════════════════════════════════════
# ЧАСТЬ C: АНАЛИЗ ЛИФТИНГА — HENSEL δ и квадратичный остаток Q
# ═══════════════════════════════════════════════════════════════════════════════

def compute_jacobian_mod2(W_base, x0_16):
    """Якобиан J (15×15) mod 2 в точке x₀ (шаг +1 по DW₁..DW₁₅)."""
    W2_base = make_W2(W_base, x0_16)
    f0 = compute_f_fast(W_base, W2_base)
    J = []
    for i in range(15):
        row = []
        for j in range(15):
            x_j = list(x0_16)
            x_j[j + 1] = (x_j[j + 1] + 1) & MASK
            W2_j = make_W2(W_base, x_j)
            fj = compute_f_fast(W_base, W2_j)
            row.append((fj[i] - f0[i]) & MASK)
        J.append(row)
    return J, f0


def analyze_hensel_lift(W_base, x0_16):
    """
    Проверить шаг Хенселя: найти δ s.t. J̃·δ ≡ q̃ (mod 2)
    где q̃ = (F(x₀)/2) mod 2, затем проверить F(x₀+2δ) mod 4.
    """
    J, f0 = compute_jacobian_mod2(W_base, x0_16)

    # q̃ = (F(x₀)/2) mod 2 = bit 1 of F(x₀)
    q_tilde = [(v >> 1) & 1 for v in f0]

    # Решить J̃·δ ≡ q̃ mod 2
    J_mod2 = [[J[i][j] & 1 for j in range(15)] for i in range(15)]
    delta_h = gauss_mod2_solve(J_mod2, q_tilde)

    result = {
        'rank_J': rank_gf2(J_mod2),
        'q_tilde': q_tilde,
        'delta_hensel': delta_h,
        'f0_mod4': [v % 4 for v in f0],
    }

    if delta_h is not None:
        # Проверить F(x₀ + 2·δ_h) mod 4
        x_lift = list(x0_16)
        for j in range(15):
            x_lift[j + 1] = (x_lift[j + 1] + 2 * delta_h[j]) & MASK
        W2_lift = make_W2(W_base, x_lift)
        f_lift = compute_f_fast(W_base, W2_lift)
        result['f_hensel_lift_mod4'] = [v % 4 for v in f_lift]
        result['lift_success'] = all(v == 0 for v in result['f_hensel_lift_mod4'])

        # Квадратичный остаток Q: разница между реальным F и Hensel-предсказанием
        # Hensel предсказывает F(x₀+2δ) ≡ F(x₀) + 2J·δ mod 4 = 0 mod 4
        # Реально: F(x₀+2δ) = F(x₀) + 2J·δ + 2Q + O(4) — вот Q:
        Jd = [sum(J[i][j] * delta_h[j] for j in range(15)) & MASK for i in range(15)]
        Q = [((f_lift[i] - f0[i] - 2 * Jd[i]) & MASK) % 4 for i in range(15)]
        result['Q_residual'] = Q
        result['Q_nonzero_count'] = sum(1 for q in Q if q != 0)

    return result


# ═══════════════════════════════════════════════════════════════════════════════
# ГЛАВНЫЙ ЗАПУСК
# ═══════════════════════════════════════════════════════════════════════════════

def run_г2_analysis(n_seeds=10):
    print("=" * 72)
    print("П-49А: ГИПОТЕЗА Г2 — РАНГ ГЕССИАНА КВАДРАТИЧНОГО БАРЬЕРА")
    print("=" * 72)
    print(f"Анализируем {n_seeds} GF(2)-семян.\n")

    hessian_ranks = []
    image_sizes = []
    zeros_found_total = 0
    seeds_with_zero = 0
    hensel_successes = 0
    Q_nonzero_counts = []
    rank_J_list = []

    for seed_idx in range(n_seeds):
        W_base = [random.randint(0, MASK) for _ in range(16)]

        # Найти GF(2) семя (пробуем разные DW₀)
        x0_16 = None
        for dw0_try in range(1, 200):
            dw0 = random.randint(1, MASK)
            seeds = find_seeds_mod2(W_base, dw0, n_trials=5000)
            if seeds:
                x0_16 = seeds[0]
                break

        if x0_16 is None:
            print(f"  Семя {seed_idx+1}: не найдено, пропуск.")
            continue

        print(f"\n─── Семя {seed_idx+1} / {n_seeds} ─────────────────────────────")
        print(f"  DW₀ = {x0_16[0]:08x}")
        dw_bits = sum((x0_16[j+1] & 1) << j for j in range(15))
        print(f"  x₁..x₁₅ (bits) = {dw_bits:#017b}")

        # A. Гессиан
        print("  [A] Вычисляю матрицу Гессиана (15×105 над GF(2))...")
        H, f0, _ = compute_hessian_matrix(W_base, x0_16)
        r_H = rank_gf2(H)
        hessian_ranks.append(r_H)
        print(f"      rank(H) = {r_H} / 15")
        print(f"      F(x₀) mod 4 = {[v%4 for v in f0[:5]]}... (первые 5)")
        nz = sum(1 for row in H for v in row if v)
        print(f"      Ненулевых бит в H: {nz}/{15*105} ({nz*100//(15*105)}%)")

        # B. Полный перебор δ
        print("  [B] Полный перебор δ ∈ {0,1}^15 (32768 вариантов)...")
        n_zeros, img_size, zero_vecs = full_scan_mod4(W_base, x0_16, max_delta=32768)
        image_sizes.append(img_size)
        zeros_found_total += n_zeros
        if n_zeros > 0:
            seeds_with_zero += 1
        print(f"      Образ |Im(g)| = {img_size} из 4^15 = {4**15} возможных")
        print(f"      Образ/пространство: {img_size}/{1<<15} (из {1<<15} δ)")
        print(f"      Нулей mod 4: {n_zeros}   {'← !!!РЕШЕНИЕ!!!' if n_zeros > 0 else '← нет решений'}")

        # C. Hensel лифт
        print("  [C] Анализ шага Хенселя...")
        h_res = analyze_hensel_lift(W_base, x0_16)
        rank_J_list.append(h_res['rank_J'])
        print(f"      rank(J̃) = {h_res['rank_J']}")
        if h_res['delta_hensel'] is not None:
            print(f"      Hensel δ найден, F(x₀+2δ) mod 4 = {h_res['f_hensel_lift_mod4'][:5]}...")
            print(f"      Успех лифта: {h_res['lift_success']}")
            print(f"      Q ненулевых компонент: {h_res['Q_nonzero_count']}/15  ← квадратичный остаток")
            Q_nonzero_counts.append(h_res['Q_nonzero_count'])
            if h_res['lift_success']:
                hensel_successes += 1
        else:
            print("      Hensel δ: нет (якобиан сингулярен, rank<15)")

    # ── Итоговая статистика ──
    print("\n" + "=" * 72)
    print("ИТОГ")
    print("=" * 72)
    print(f"Семян проанализировано: {len(hessian_ranks)}")
    print()
    print(f"[Г2 — ГЕССИАН]")
    if hessian_ranks:
        print(f"  rank(H) по семенам: {hessian_ranks}")
        print(f"  E[rank(H)] = {sum(hessian_ranks)/len(hessian_ranks):.2f}")
        print(f"  P(rank=15) = {sum(1 for r in hessian_ranks if r==15)/len(hessian_ranks):.2f}")
        if all(r == 15 for r in hessian_ranks):
            print("  → Г2 ПОДТВЕРЖДЕНА: Гессиан полноранговый во ВСЕХ семенах!")
        elif any(r == 15 for r in hessian_ranks):
            print("  → Г2 ЧАСТИЧНО: Гессиан полноранговый в некоторых семенах.")
        else:
            print("  → Г2 НЕ ПОДТВЕРЖДЕНА: Гессиан не имеет полного ранга.")
    print()
    print(f"[Г1 — ПОИСК MOD-4 РЕШЕНИЙ]")
    if image_sizes:
        print(f"  Размер образа Im(g): {image_sizes}")
        print(f"  E[|Im(g)|] = {sum(image_sizes)/len(image_sizes):.1f} / 32768 δ проверено")
        coverage = [s/(1<<15) for s in image_sizes]
        print(f"  Покрытие δ-пространства: {[f'{c:.3f}' for c in coverage]}")
    print(f"  Нулей mod 4 найдено: {zeros_found_total} / ({len(image_sizes)}×32768)")
    print(f"  Семян с решением: {seeds_with_zero}/{len(image_sizes)}")
    if zeros_found_total == 0:
        print("  → Г1 НЕ ОПРОВЕРГНУТА: решений mod 4 не найдено (в данной выборке).")
    else:
        print(f"  → Г1 ОПРОВЕРГНУТА: найдено {zeros_found_total} решений mod 4!")
    print()
    print(f"[КВАДРАТИЧНЫЙ ОСТАТОК Q]")
    if Q_nonzero_counts:
        print(f"  Q ненулевых компонент: {Q_nonzero_counts}")
        print(f"  E[|Q ≠ 0|] = {sum(Q_nonzero_counts)/len(Q_nonzero_counts):.1f} / 15")
        print(f"  Hensel успехов: {hensel_successes}/{len(Q_nonzero_counts)}")
    print()
    print(f"[РАНГ ЯКОБИАНА J]")
    if rank_J_list:
        print(f"  rank(J̃) по семенам: {rank_J_list}")
        print(f"  E[rank(J̃)] = {sum(rank_J_list)/len(rank_J_list):.2f}")
    print()
    print("ВЫВОД:")
    if hessian_ranks and all(r == 15 for r in hessian_ranks) and zeros_found_total == 0:
        print("  ✓ Г2 ПОДТВЕРЖДЕНА: Гессиан полноранговый.")
        print("  ✓ Г1 НЕ ОПРОВЕРГНУТА: решений mod 4 нет.")
        print("  Теорема P1 (mod-4 барьер) получает структурное объяснение:")
        print("  квадратичный барьер Q полноранговый → лифтинг невозможен.")
    elif zeros_found_total > 0:
        print("  ! Г1 ОПРОВЕРГНУТА: решения mod 4 существуют!")
    else:
        print("  Данные неоднозначны, нужно больше семян.")


if __name__ == '__main__':
    random.seed(42)
    run_г2_analysis(n_seeds=10)
