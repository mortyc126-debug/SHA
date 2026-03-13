"""
П-49Б: ГИПОТЕЗА Г1 — СУЩЕСТВОВАНИЕ РЕШЕНИЙ F(x) ≡ 0 mod 4

Вопрос: Существует ли x ∈ (Z/4Z)^15 (с DW₀=const) с F(x) ≡ 0 mod 4?

Три стратегии поиска:

  Стратегия A (случайный поиск):
    Случайные x ∈ {0,1,2,3}^15 → вычислить F(x) mod 4.
    Если F равномерна → P(F=0 mod 4) = 4^{-15} ≈ 10^{-9}.
    Задача: проверить, есть ли структурное смещение.

  Стратегия B (жадный mod-4 каскад):
    Последовательно выбирать DW₂,DW₃,...,DW₁₅ ∈ {0,1,2,3}
    чтобы каждый Da_{k} ≡ 0 mod 4 (перебор 4 вариантов на шаге).
    Цель: проверить, существует ли mod-4 каскад и какой De₁₇ он даёт.

  Стратегия C (GF(2) → mod-4 подъём, полный):
    Для каждого (W_base, DW₀): найти ВСЕ GF(2)-семена (не только первое).
    Для каждого семени: проверить ВСЕ 2^15 лифтов δ.
    Статистика: total_attempts / solutions.

Ключевой вопрос Б.A vs Б.B:
  - Стратегия A: mod-4 "маленькие" x (близкие к нулю по mod-2³²)
  - Стратегия B: mod-4 "каскадные" x (близкие к структуре adap cascade)
  Разница важна: в Серии II мы изучаем МАЛЕНЬКИЕ x (GF(2) lifting).
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


def _expand17(W16):
    W = list(W16)
    for i in range(16, 18):
        W.append((sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16]) & MASK)
    return W


def _sha256_states_17(W16):
    W = _expand17(W16)
    a, b, c, d, e, f, g, h = H0
    states = {}
    for r in range(17):
        T1 = (h + Sig1(e) + Ch(e, f, g) + K_SHA[r] + W[r]) & MASK
        T2 = (Sig0(a) + Maj(a, b, c)) & MASK
        h = g; g = f; f = e; e = (d + T1) & MASK
        d = c; c = b; b = a; a = (T1 + T2) & MASK
        if r >= 2:
            states[r + 1] = (a, b, c, d, e, f, g, h)
    return states


_base_cache = {}


def compute_f_vec(W_base, x16):
    """F(x) = (Da₃..Da₁₆, De₁₇). x16 — список из 16 значений DW."""
    W2 = [(W_base[i] + x16[i]) & MASK for i in range(16)]
    key = tuple(W_base)
    if key not in _base_cache:
        _base_cache.clear()
        _base_cache[key] = _sha256_states_17(W_base)
    s1 = _base_cache[key]
    s2 = _sha256_states_17(W2)
    result = [(s2[r][0] - s1[r][0]) & MASK for r in range(3, 17)]
    result.append((s2[17][4] - s1[17][4]) & MASK)
    return result


def hamming_weight_mod4(f_mod4):
    """Число ненулевых компонент вектора mod 4."""
    return sum(1 for v in f_mod4 if v != 0)


# ═══════════════════════════════════════════════════════════════════════════════
# СТРАТЕГИЯ A: СЛУЧАЙНЫЙ ПОИСК В {0,1,2,3}^15
# ═══════════════════════════════════════════════════════════════════════════════

def strategy_a_random(W_base, dw0, n_trials=50000):
    """
    Случайный поиск x ∈ {0,1,2,3}^15 с F(x) ≡ 0 mod 4.
    Дополнительно: минимальный L1-вес F(x) mod 4, распределение числа нулей.
    """
    best_zeros = 0
    best_x = None
    zeros_hist = {}  # histogram: zeros_hist[k] = count of x with k components=0 mod 4
    solutions = []

    for _ in range(n_trials):
        x16 = [dw0] + [random.randint(0, 3) for _ in range(15)]
        f = compute_f_vec(W_base, x16)
        f4 = [v % 4 for v in f]
        n_zero = sum(1 for v in f4 if v == 0)
        zeros_hist[n_zero] = zeros_hist.get(n_zero, 0) + 1
        if n_zero > best_zeros:
            best_zeros = n_zero
            best_x = x16[:]
        if n_zero == 15:
            solutions.append(x16[:])

    return best_zeros, best_x, zeros_hist, solutions


# ═══════════════════════════════════════════════════════════════════════════════
# СТРАТЕГИЯ B: ЖАДНЫЙ MOD-4 КАСКАД
# ═══════════════════════════════════════════════════════════════════════════════

def strategy_b_greedy_cascade(W_base, dw0, dw1=0, n_tries=50):
    """
    Жадный mod-4 каскад: последовательно выбирать DW₂..DW₁₅ ∈ {0,1,2,3}
    чтобы максимизировать число Da_k ≡ 0 mod 4.

    Ключевой вопрос: насколько глубоко можно зайти?
    (В Серии I каскад работал ТОЧНО для любого DW₀ — за счёт больших значений DW)
    Здесь DW ограничены 4-мя значениями — сможет ли это работать?
    """
    results = []

    for attempt in range(n_tries):
        # Варьируем DW₁ по попыткам чтобы исследовать разные ветки каскада
        dw1_try = attempt % 4
        x16 = [dw0, dw1_try] + [0] * 14
        zeros_achieved = 0

        for pos in range(2, 16):  # DW₂..DW₁₅
            best_v = 0
            best_zeros = -1
            best_da_abs = MASK

            for v in range(4):  # Перебрать {0,1,2,3}
                x_try = list(x16)
                x_try[pos] = v
                # Da_{pos+1}: индекс в F = pos - 2 (Da3=idx0, Da4=idx1, ...)
                da = compute_f_vec(W_base, x_try)[pos - 2]
                da_mod4 = da % 4
                z = 1 if da_mod4 == 0 else 0
                da_abs = min(da_mod4, 4 - da_mod4)  # расстояние до нуля mod 4

                if z > best_zeros or (z == best_zeros and da_abs < best_da_abs):
                    best_zeros = z
                    best_v = v
                    best_da_abs = da_abs

            x16[pos] = best_v
            if best_zeros == 1:
                zeros_achieved += 1

        # Проверить полный результат
        f = compute_f_vec(W_base, x16)
        f4 = [v % 4 for v in f]
        n_total_zero = sum(1 for v in f4 if v == 0)

        results.append({
            'zeros_in_cascade': zeros_achieved,
            'total_zeros_F_mod4': n_total_zero,
            'De17_mod4': f4[-1],
            'x16': x16[:],
        })

    # Лучший результат
    best = max(results, key=lambda r: r['total_zeros_F_mod4'])
    avg_zeros = sum(r['total_zeros_F_mod4'] for r in results) / n_tries

    return best, avg_zeros, results


# ═══════════════════════════════════════════════════════════════════════════════
# СТРАТЕГИЯ C: GF(2) СЕМЕНА + ПОЛНЫЙ ПЕРЕБОР ЛИФТОВ
# ═══════════════════════════════════════════════════════════════════════════════

def find_all_seeds_mod2_exhaustive(W_base, dw0):
    """
    ИСЧЕРПЫВАЮЩИЙ поиск всех x₀ ∈ {0,1}^15 с F(x₀) ≡ 0 mod 2.
    Проверяет все 2^15 = 32768 кандидатов.
    """
    seeds = []
    for bits in range(1 << 15):
        x16 = [dw0] + [(bits >> j) & 1 for j in range(15)]
        W2 = [(W_base[i] + x16[i]) & MASK for i in range(16)]
        f = compute_f_vec(W_base, x16)
        if all(v % 2 == 0 for v in f):
            seeds.append(x16[:])
    return seeds


def strategy_c_full_lift(W_base, dw0, exhaustive_seeds=True, n_seed_trials=5000, max_seeds_for_lift=5):
    """
    Для каждого GF(2)-семени: проверить все 2^15 лифтов δ.
    Вернуть: total_attempts, solutions_found, seeds_count.
    """
    if exhaustive_seeds:
        seeds = find_all_seeds_mod2_exhaustive(W_base, dw0)
    else:
        # Случайный поиск семян
        seeds = []
        for _ in range(n_seed_trials):
            bits = random.randrange(1 << 15)
            x16 = [dw0] + [(bits >> j) & 1 for j in range(15)]
            f = compute_f_vec(W_base, x16)
            if all(v % 2 == 0 for v in f):
                seeds.append(x16[:])

    total_attempts = 0
    solutions = []

    for seed in seeds[:max_seeds_for_lift]:
        for dbits in range(1 << 15):
            delta = [(dbits >> j) & 1 for j in range(15)]
            x_lift = [seed[0]] + [(seed[j+1] + 2*delta[j]) & MASK for j in range(15)]
            f = compute_f_vec(W_base, x_lift)
            total_attempts += 1
            if all(v % 4 == 0 for v in f):
                solutions.append({'seed': seed[:], 'delta': delta[:], 'x': x_lift[:]})

    return len(seeds), total_attempts, solutions


# ═══════════════════════════════════════════════════════════════════════════════
# ДОПОЛНИТЕЛЬНО: РАСПРЕДЕЛЕНИЕ |F(x) mod 4| ДЛЯ СЛУЧАЙНЫХ x
# ═══════════════════════════════════════════════════════════════════════════════

def analyze_f_mod4_distribution(W_base, dw0, n_samples=10000):
    """
    Распределение ненулевых компонент F(x) mod 4 для случайных x ∈ {0,1,2,3}^15.
    Позволяет оценить, насколько далеко F от нуля mod 4.
    """
    dist = {}
    min_nonzero = 15
    best_x = None

    for _ in range(n_samples):
        x16 = [dw0] + [random.randint(0, 3) for _ in range(15)]
        f = compute_f_vec(W_base, x16)
        f4 = tuple(v % 4 for v in f)
        n_nz = sum(1 for v in f4 if v != 0)
        dist[n_nz] = dist.get(n_nz, 0) + 1
        if n_nz < min_nonzero:
            min_nonzero = n_nz
            best_x = x16[:]

    return dist, min_nonzero, best_x


# ═══════════════════════════════════════════════════════════════════════════════
# ГЛАВНЫЙ ЗАПУСК
# ═══════════════════════════════════════════════════════════════════════════════

def run_г1_analysis(n_bases=5):
    print("=" * 72)
    print("П-49Б: ГИПОТЕЗА Г1 — СУЩЕСТВОВАНИЕ РЕШЕНИЙ F(x) ≡ 0 mod 4")
    print("=" * 72)
    print(f"Анализируем {n_bases} базовых блоков W_base.\n")

    all_solutions = []
    total_attempts_c = 0
    total_seeds_c = 0
    global_best_zeros = 0
    cascade_de17_mod4_hist = {}

    for base_idx in range(n_bases):
        W_base = [random.randint(0, MASK) for _ in range(16)]
        dw0 = random.randint(1, MASK)

        print(f"\n{'─'*72}")
        print(f"БАЗА {base_idx+1}/{n_bases}  DW₀={dw0:#010x}")

        # ── Стратегия A ──
        print("  [A] Случайный поиск x ∈ {0,1,2,3}^15 (50K попыток)...")
        best_z, best_x, zhist, sols_a = strategy_a_random(W_base, dw0, n_trials=50000)
        global_best_zeros = max(global_best_zeros, best_z)
        print(f"      Лучший результат: {best_z}/15 компонент = 0 mod 4")
        top5 = sorted(zhist.items(), reverse=True)[:5]
        print(f"      Распределение (zeros: count): {top5}")
        if sols_a:
            print(f"      !!! НАЙДЕНО {len(sols_a)} РЕШЕНИЙ (A)! !!!")
            all_solutions.extend(sols_a)

        # ── Стратегия B ──
        print("  [B] Жадный mod-4 каскад (50 попыток)...")
        best_res, avg_z, cascade_results = strategy_b_greedy_cascade(W_base, dw0, n_tries=50)
        print(f"      Лучший: {best_res['total_zeros_F_mod4']}/15 нулей F mod 4")
        print(f"      Среднее: {avg_z:.1f}/15")
        print(f"      Нулей в шагах каскада: {best_res['zeros_in_cascade']}/14")
        de17_v = best_res['De17_mod4']
        cascade_de17_mod4_hist[de17_v] = cascade_de17_mod4_hist.get(de17_v, 0) + 1
        print(f"      De₁₇ mod 4 (лучший): {de17_v}")
        if best_res['total_zeros_F_mod4'] == 15:
            print(f"      !!! НАЙДЕНО РЕШЕНИЕ (B)! !!!")
            all_solutions.append(best_res['x16'])

        # ── Стратегия C ──
        print("  [C] GF(2) семена → полный перебор лифтов (исчерпывающий по δ, макс 5 семян)...")
        n_seeds, n_att, sols_c = strategy_c_full_lift(W_base, dw0, exhaustive_seeds=True, max_seeds_for_lift=5)
        total_seeds_c += n_seeds
        total_attempts_c += n_att
        print(f"      GF(2) семян найдено: {n_seeds} из 32768 проверенных")
        print(f"      Всего лифтов проверено: {n_att}")
        print(f"      Решений найдено: {len(sols_c)}")
        if sols_c:
            print(f"      !!! НАЙДЕНО {len(sols_c)} РЕШЕНИЙ (C)! !!!")
            all_solutions.extend([s['x'] for s in sols_c])

        # ── Распределение F mod 4 ──
        print("  [D] Распределение |F(x) mod 4| (10K случ.)...")
        dist, min_nz, _ = analyze_f_mod4_distribution(W_base, dw0, n_samples=10000)
        sorted_dist = sorted(dist.items())
        print(f"      #нулей : #случаев (из 10K)")
        for k, cnt in sorted_dist:
            bar = "█" * (cnt // 100)
            print(f"         {15-k:2d} ненул. ({k:2d} нул.): {cnt:5d}  {bar}")
        print(f"      Минимум ненулевых: {min_nz}")

    # ── Итоговая статистика ──
    print("\n" + "=" * 72)
    print("ИТОГ — ГИПОТЕЗА Г1")
    print("=" * 72)
    print(f"Общая статистика:")
    print(f"  Стратегия C: семян = {total_seeds_c}, лифтов = {total_attempts_c}")
    print(f"  Лучший результат стратегии A: {global_best_zeros}/15 нулей mod 4")
    print(f"  De₁₇ mod 4 в жадном каскаде: {cascade_de17_mod4_hist}")
    print()

    if all_solutions:
        print(f"!!! НАЙДЕНО {len(all_solutions)} РЕШЕНИЙ F(x) ≡ 0 mod 4 !!!")
        print("Г1 ОПРОВЕРГНУТА — решения СУЩЕСТВУЮТ.")
        for i, sol in enumerate(all_solutions[:3]):
            f = compute_f_vec([0]*16, [0]*16)  # dummy
            print(f"  Решение {i+1}: DW₀={sol[0]:#010x}, DW₁..₁₅={sol[1:]}")
    else:
        print(f"Решений mod 4 не найдено.")
        print(f"  Суммарно проверено: {total_attempts_c} + {n_bases*50000} = {total_attempts_c + n_bases*50000}")
        print()
        print("ОЦЕНКА:")
        if global_best_zeros >= 12:
            print(f"  Максимум {global_best_zeros}/15 нулей → близко к решению,")
            print(f"  но не достигнуто. Вероятно, структурный барьер.")
        else:
            print(f"  Максимум {global_best_zeros}/15 нулей → далеко от решения.")
            print(f"  Сильный структурный барьер подтверждается.")
        print()
        print("ВЫВОД: Г1 пока не опровергнута.")
        print("  Все три стратегии дают 0 решений mod 4.")
        print("  Согласуется с T_MOD4_BARRIER_ABSOLUTE и Г2 (если подтверждена).")

    print()
    print("Г1 и Г2 взаимосвязь:")
    print("  Если Г2 верна (Гессиан ранг 15):")
    print("    → квадратичный барьер фундаментален")
    print("    → решений mod 4 из GF(2)-семян нет НИКОГДА")
    print("    → Г1 (в контексте GF(2)-лифтинга) решена отрицательно")
    print("  Если Г1 позитивна (решения найдены):")
    print("    → Г2 ложна, или лифт идёт не из GF(2)-семян")


if __name__ == '__main__':
    random.seed(1337)
    run_г1_analysis(n_bases=5)
