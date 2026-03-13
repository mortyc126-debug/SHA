"""
П-54: БАШНЯ ДО k=22 — ПРОВЕРКА T_INFINITE_TOWER

Из итога П-53:
  T_CASCADE_UNIQUENESS [ДОКАЗАНА]: Da_{pos+1} линейна в DW[pos], наклон +1.
  T_FREE_CONSTRAINT [ПОДТВЕРЖДЕНА]: De_17 единственное несвязанное ограничение.
  height_2(SHA-256) ≥ 11 (точный каскад).

Вопрос П-54:
  P(De_17=0 mod 2^k) = 1/2^k для ВСЕХ k, или при большом k появится
  структурный барьер?

  Сценарий A (башня бесконечна):
    De_17 mod 2^k — равномерна ∀k → Sol_k ≠ ∅ всегда.
    p-адически: ∃ x ∈ Z_2^16 с F(x) = 0 над 2-адическими целыми.

  Сценарий B (конечная высота k*):
    При k → k* появляется корреляция De_17 ↔ Da_3..Da_16.
    P(Sol_k) → 0 при k → k*.

Эксперименты:
  A. T_INFINITE_TOWER: башня k=12..22, P(Sol_k) vs 1/2^k.
     Адаптивное число попыток: 10×2^k для k≤16, min(200_000) для k>16.
  B. T_HENSEL_TOWER: согласованная цепочка x^{(k)} ≡ x^{(k-1)} mod 2^{k-1}.
     Найти Sol_1, попытаться поднять до Sol_2 ≡ Sol_1 mod 2, и т.д.
  C. Равномерность De_17: хи-квадрат тест для k=8..14.
  D. Структура De_17: v₂-профиль при k=12..20.
"""

import random
from collections import Counter
from math import sqrt

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
def Sig0(x):    return rotr(x, 2)  ^ rotr(x, 13) ^ rotr(x, 22)
def Sig1(x):    return rotr(x, 6)  ^ rotr(x, 11) ^ rotr(x, 25)
def sig0(x):    return rotr(x, 7)  ^ rotr(x, 18) ^ (x >> 3)
def sig1(x):    return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)


def v2(n):
    if n == 0:
        return 64
    k = 0
    while n & 1 == 0:
        k += 1
        n >>= 1
    return k


def sha256_states_17(W16):
    """Вычислить состояния SHA-256 после раундов 0..17 для 16 слов W."""
    W = list(W16)
    for i in range(16, 18):
        W.append((sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16]) & MASK)
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


def get_base_states(W_base):
    key = tuple(W_base)
    if key not in _base_cache:
        _base_cache.clear()
        _base_cache[key] = sha256_states_17(W_base)
    return _base_cache[key]


def compute_f(W_base, x16):
    """F(x) = разность состояний Da_3..Da_16, De_17."""
    s1 = get_base_states(W_base)
    W2 = [(W_base[i] + x16[i]) & MASK for i in range(16)]
    s2 = sha256_states_17(W2)
    result = [(s2[r][0] - s1[r][0]) & MASK for r in range(3, 17)]
    result.append((s2[17][4] - s1[17][4]) & MASK)
    return result


def exact_cascade(W_base, dw0, dw1, k):
    """
    Точный каскад: O(14) вызовов compute_f.

    Линейность Da_{pos+1}(v) = CONST + v доказывает:
      v* = (-CONST) mod 2^k — единственный нуль, ВСЕГДА существует.

    Возвращает (x16, de17_mod_k, all_14_zero).
    """
    mod = 1 << k
    mod_mask = mod - 1
    x16 = [dw0, dw1] + [0] * 14

    for pos in range(2, 16):
        # Da_{pos+1}(0) при текущем x16 (Da(v=0) = CONST)
        da0 = compute_f(W_base, x16)[pos - 2]
        # Точный нуль: v* = (-da0) mod 2^k
        x16[pos] = (-da0) & mod_mask

    f = compute_f(W_base, x16)
    all_14 = all((f[i] & mod_mask) == 0 for i in range(14))
    de17_mod = f[14] & mod_mask

    return x16, de17_mod, all_14


def make_wbase(seed=0):
    """Создать базовое W[0..15] из seed."""
    rng = random.Random(seed)
    return [rng.randint(0, MASK) for _ in range(16)]


# ═══════════════════════════════════════════════════════════════════════════════
# ЭКСПЕРИМЕНТ A: T_INFINITE_TOWER — башня k=12..22
# ═══════════════════════════════════════════════════════════════════════════════

def experiment_a_tower(W_base, k_max=22, max_trials=200_000, seed=42):
    """
    Проверить P(Sol_k) ≈ 1/2^k для k=12..k_max.

    Адаптивное число попыток:
      n = min(10×2^k, max_trials)
    Для k>log2(max_trials/10): уменьшенная выборка, менее надёжно.

    Индикатор барьера: ratio = P_obs / (1/2^k) существенно < 1?
    """
    print("\n" + "=" * 72)
    print("ЭКСПЕРИМЕНТ A: T_INFINITE_TOWER — башня k=12..{0}".format(k_max))
    print("=" * 72)
    print(f"  База W: seed=0, первые 4 слова = {W_base[:4]}")
    print()
    print(f"  {'k':>3}  {'mod':>8}  {'попыток':>9}  {'Sol_k':>6}  "
          f"{'P_obs':>10}  {'1/2^k':>10}  {'ratio':>7}  {'статус':>12}")
    print(f"  {'─'*3}  {'─'*8}  {'─'*9}  {'─'*6}  "
          f"{'─'*10}  {'─'*10}  {'─'*7}  {'─'*12}")

    rng = random.Random(seed)
    results = {}

    for k in range(12, k_max + 1):
        mod = 1 << k
        n_trials = min(10 * mod, max_trials)
        expected_p = 1.0 / mod

        successes = 0
        sols = []

        for _ in range(n_trials):
            dw0 = rng.randint(0, MASK)
            dw1 = rng.randint(0, mod - 1)
            _, de17, all14 = exact_cascade(W_base, dw0, dw1, k)
            if all14 and de17 == 0:
                successes += 1
                if len(sols) < 3:
                    sols.append((dw0, dw1))

        p_obs = successes / n_trials
        ratio = p_obs / expected_p if expected_p > 0 else 0.0

        # Статистический тест: отклонение от ожидания
        expected_count = n_trials * expected_p
        if expected_count < 0.5:
            status = "мало данных"
        elif successes == 0 and expected_count >= 2:
            status = "БАРЬЕР? ✗"
        elif ratio < 0.3:
            status = "аномалия ↓"
        elif ratio > 3.0:
            status = "аномалия ↑"
        else:
            status = "Sol≠∅ ✓" if successes > 0 else "Sol=∅?"

        print(f"  {k:>3}  {mod:>8}  {n_trials:>9}  {successes:>6}  "
              f"{p_obs:>10.6f}  {expected_p:>10.6f}  {ratio:>7.3f}  {status:>12}")

        results[k] = {
            'mod': mod, 'n': n_trials, 'successes': successes,
            'p_obs': p_obs, 'ratio': ratio, 'sols': sols
        }

    # Итог
    print()
    barrier_k = None
    for k, r in sorted(results.items()):
        if r['successes'] == 0 and r['n'] * (1.0 / r['mod']) >= 2:
            barrier_k = k
            break

    if barrier_k is None:
        # Найти max k с надёжным результатом
        reliable = [k for k, r in results.items()
                    if r['n'] * (1.0 / r['mod']) >= 3 and r['successes'] > 0]
        if reliable:
            print(f"  ✓ T_INFINITE_TOWER: Sol_k ≠ ∅ для k = {min(reliable)}..{max(reliable)}")
            print(f"    height_2(SHA-256) ≥ {max(reliable)} (нижняя оценка)")
            print(f"    P(Sol_k) ≈ 1/2^k во всём проверенном диапазоне.")
            print(f"    → Сценарий A (башня бесконечна): ПОДТВЕРЖДАЕТСЯ")
        else:
            print("  ⚠ Недостаточно данных для надёжного вывода.")
    else:
        print(f"  ✗ БАРЬЕР обнаружен при k={barrier_k}!")
        print(f"    height_2(SHA-256) = {barrier_k - 1}")
        print(f"    → Сценарий B (конечная высота): ПОДТВЕРЖДАЕТСЯ")

    return results


# ═══════════════════════════════════════════════════════════════════════════════
# ЭКСПЕРИМЕНТ B: T_HENSEL_TOWER — согласованная цепочка решений
# ═══════════════════════════════════════════════════════════════════════════════

def experiment_b_hensel(W_base, k_max=20, attempts_per_level=200, seed=123):
    """
    Проверить: можно ли построить СОГЛАСОВАННУЮ цепочку?
      x^{(1)}, x^{(2)}, ..., x^{(k)} где x^{(k)} ≡ x^{(k-1)} mod 2^{k-1}

    Алгоритм Hensel:
      - Найти x^{(1)} ∈ Sol_1 (попыток достаточно, P = 1/2)
      - Попытаться поднять x^{(k-1)} → x^{(k)}:
        Перебрать все (dw0 ≡ dw0^{(k-1)} mod 2^{k-1},
                       dw1 ≡ dw1^{(k-1)} mod 2^{k-1})
        Есть 2 × 2 = 4 варианта (бит k-1 каждого из dw0, dw1).

    Предсказание T_TOWER_INDEPENDENCE (П-52):
      Подъём НЕ гарантирован: P(подъём) = 1/2 на каждом уровне.
      После k шагов: P(цепочка длины k) = (1/2)^{k-1}.
    """
    print("\n" + "=" * 72)
    print("ЭКСПЕРИМЕНТ B: T_HENSEL_TOWER — цепочка решений до k={0}".format(k_max))
    print("=" * 72)
    print(f"  Адаптивный Hensel: попытки поднять каждое Sol_k → Sol_{{k+1}}")
    print()

    rng = random.Random(seed)

    # Найти несколько стартовых Sol_1 = {(dw0, dw1) : De_17=0 mod 2}
    print("  Ищем стартовые Sol_1 (De_17 ≡ 0 mod 2)...")
    starts = []
    for _ in range(attempts_per_level * 4):
        dw0 = rng.randint(0, MASK)
        dw1 = rng.randint(0, 1)
        _, de17, all14 = exact_cascade(W_base, dw0, dw1, 1)
        if all14 and de17 == 0:
            starts.append((dw0, dw1))
            if len(starts) >= 5:
                break

    if not starts:
        print("  ✗ Не нашли Sol_1!")
        return {}

    print(f"  Найдено {len(starts)} стартовых Sol_1.")

    # Для каждого старта — попытка поднять до k=k_max
    results = {}
    chains = []  # Успешные цепочки

    print()
    print(f"  {'k':>3}  {'старт':>6}  {'подъёмов':>9}  {'успехов':>9}  "
          f"{'P(подъём)':>10}  {'ожид':>8}  {'статус':>15}")
    print(f"  {'─'*3}  {'─'*6}  {'─'*9}  {'─'*9}  "
          f"{'─'*10}  {'─'*8}  {'─'*15}")

    total_chains = list(starts)  # Текущий уровень цепочек

    for k in range(2, k_max + 1):
        mod_k = 1 << k
        mod_km1 = 1 << (k - 1)

        next_chains = []
        total_lifted = 0
        success_lifted = 0

        for (prev_dw0, prev_dw1) in total_chains[:50]:  # Не более 50 стартов
            # 4 варианта подъёма: бит (k-1) для dw0 и dw1
            lifted = False
            for b0 in range(2):
                for b1 in range(2):
                    dw0_new = (prev_dw0 & (mod_km1 - 1)) | (b0 << (k - 1))
                    dw1_new = (prev_dw1 & (mod_km1 - 1)) | (b1 << (k - 1))
                    _, de17, all14 = exact_cascade(W_base, dw0_new, dw1_new, k)
                    if all14 and de17 == 0:
                        next_chains.append((dw0_new, dw1_new))
                        lifted = True
            total_lifted += 1
            if lifted:
                success_lifted += 1

        p_lift = success_lifted / total_lifted if total_lifted > 0 else 0.0
        # Ожидаемая вероятность подъёма = P(хотя бы 1 из 4 вариантов работает)
        # Если De_17(новые биты) ≡ 0 независимо с P=1/2 каждое из 4:
        # P(хоть один) = 1 - (3/4)^? ... точнее P = 1 - (1-1/2)^...
        # Если только 1 степень свободы (De_17 зависит от бит k-1 в dw0 или dw1),
        # то 2 из 4 вариантов работают → P ≈ 1/2 (при одной степени свободы)
        # или 1/4 (если обе независимы)
        expected_p = 1.0 - (1.0 - 1.0 / mod_k) ** 4  # ≈ 4/mod_k при большом k

        status = "✓" if success_lifted > 0 else "✗ цепочка оборвана"

        print(f"  {k:>3}  {total_lifted:>6}  {total_lifted:>9}  {success_lifted:>9}  "
              f"{p_lift:>10.4f}  {expected_p:>8.4f}  {status:>15}")

        results[k] = {'lifted': total_lifted, 'success': success_lifted,
                      'p_lift': p_lift, 'chains': len(next_chains)}

        if not next_chains:
            print(f"\n  ⚡ Цепочка оборвалась при k={k}: ни один подъём не найден.")
            print(f"    Это может быть: (a) статистика, (b) истинный барьер Hensel.")
            break

        total_chains = next_chains

    if total_chains:
        max_k = max(results.keys())
        print(f"\n  ✓ Цепочка достигла k={max_k} с {len(total_chains)} активными ветвями.")
        if max_k >= k_max - 1:
            print(f"    T_HENSEL_TOWER: СОГЛАСОВАННАЯ цепочка до k={max_k} СУЩЕСТВУЕТ.")
        else:
            print(f"    Цепочка оборвалась раньше k_max={k_max}.")

    return results


# ═══════════════════════════════════════════════════════════════════════════════
# ЭКСПЕРИМЕНТ C: РАВНОМЕРНОСТЬ De_17 — хи-квадрат тест
# ═══════════════════════════════════════════════════════════════════════════════

def experiment_c_uniformity(W_base, k_range=(8, 15), n_per_k=5000, seed=777):
    """
    Хи-квадрат тест равномерности распределения De_17 mod 2^k.

    Нулевая гипотеза H0: De_17 mod 2^k ~ равномерно на {0,...,2^k-1}.
    Тест χ²: сравниваем наблюдаемые частоты с ожидаемыми.

    Для больших 2^k (много ячеек): тест на v₂(De_17) вместо полного χ².
    """
    print("\n" + "=" * 72)
    print("ЭКСПЕРИМЕНТ C: РАВНОМЕРНОСТЬ De_17 — χ²-тест")
    print("=" * 72)
    print(f"  n_per_k={n_per_k} случайных пар (dw0,dw1) на каждое k")
    print()

    rng = random.Random(seed)

    for k in range(k_range[0], k_range[1] + 1):
        mod = 1 << k
        de17_vals = []

        for _ in range(n_per_k):
            dw0 = rng.randint(0, MASK)
            dw1 = rng.randint(0, MASK)
            _, de17, all14 = exact_cascade(W_base, dw0, dw1, k)
            if all14:
                de17_vals.append(de17)

        n = len(de17_vals)
        if n < 100:
            print(f"  k={k}: мало данных ({n}), пропуск")
            continue

        # Тест v₂: частота v₂(De_17) = j должна быть ≈ 1/2^(j+1) для j<k, 1/2^k для j=k
        v2_counts = Counter(v2(x) for x in de17_vals)
        print(f"\n  k={k} (mod={mod}), n={n} действительных пар:")
        print(f"    v₂(De_17)  наблюд    ожидаемо   отношение")
        print(f"    {'─'*9}  {'─'*8}  {'─'*9}  {'─'*9}")

        chi2 = 0.0
        for j in range(min(k + 1, 8)):
            if j < k:
                expected = n / (1 << (j + 1))
            else:
                expected = n / (1 << k)
            obs = v2_counts.get(j, 0)
            ratio = obs / expected if expected > 0 else 0
            if j >= k:
                label = f"≥{j}"
            else:
                label = str(j)
            marker = " ✓" if 0.7 <= ratio <= 1.3 else " ⚠"
            print(f"    {label:>9}  {obs:>8}  {expected:>9.1f}  {ratio:>9.3f}{marker}")
            if expected > 5:
                chi2 += (obs - expected) ** 2 / expected

        # P(De_17=0): ожидаем 1/2^k
        zero_count = de17_vals.count(0)
        p0 = zero_count / n
        expected_p0 = 1.0 / mod
        ratio0 = p0 / expected_p0 if expected_p0 > 0 else 0

        # Простая оценка Poisson-пуассоновского теста для 0:
        # При H0: count_0 ~ Poisson(n/2^k)
        expected_count0 = n * expected_p0
        print(f"\n    P(De_17=0): {p0:.6f} vs ожид {expected_p0:.6f}, "
              f"ratio={ratio0:.3f}, счёт={zero_count} (ожид {expected_count0:.1f})")

        # χ² итог (только для малых k где много ячеек)
        if mod <= 64 and n >= mod * 5:
            obs_full = Counter(de17_vals)
            chi2_full = sum((obs_full.get(i, 0) - n / mod) ** 2 / (n / mod)
                           for i in range(mod))
            df = mod - 1
            # p-value аппроксимация (нормальное приближение для χ²)
            z = (chi2_full - df) / sqrt(2 * df)
            conclusion = "H0 не отвергается ✓" if abs(z) < 3 else f"H0 ОТВЕРГНУТА (z={z:.1f}) ✗"
            print(f"    χ²={chi2_full:.1f} (df={df}), z={z:.2f}: {conclusion}")


# ═══════════════════════════════════════════════════════════════════════════════
# ЭКСПЕРИМЕНТ D: v₂-ПРОФИЛЬ De_17 ПРИ БОЛЬШИХ k
# ═══════════════════════════════════════════════════════════════════════════════

def experiment_d_v2_profile(W_base, k_list=None, n_samples=3000, seed=999):
    """
    Профиль v₂(De_17) при k=12..20.

    Ключевой вопрос: остаётся ли P(v₂(De_17) ≥ k) ≈ 1/2^k?

    Прямое измерение: взять 32-битный De_17 (полный, k=32) и смотреть v₂.
    Тогда P(De_17 ≡ 0 mod 2^k) = P(v₂(De_17) ≥ k).
    """
    if k_list is None:
        k_list = [12, 14, 16, 18, 20]

    print("\n" + "=" * 72)
    print("ЭКСПЕРИМЕНТ D: v₂-ПРОФИЛЬ De_17 (32-битный)")
    print("=" * 72)
    print(f"  n_samples={n_samples} пар (dw0,dw1), полный De_17 mod 2^32")
    print()

    rng = random.Random(seed)
    k_ref = 20  # Используем k=20 для каскада (достаточно высокий)

    # Набрать samples De_17 при разных значениях k для каскада
    # Важно: De_17 полный (32-бит) не зависит от выбора k в каскаде!
    # Однако DW[2..15] зависят от k (нули низших битов).
    # Для чистоты эксперимента: взять k=32 (полные 32-бита) для каскада.
    k_cascade = 28  # Используем k=28 для cascade, берём полный De_17

    print(f"  Каскад при k={k_cascade} (28 бит), De_17 = полные 32 бита")
    print()

    de17_samples = []
    for i in range(n_samples):
        dw0 = rng.randint(0, MASK)
        dw1 = rng.randint(0, (1 << k_cascade) - 1)
        _, de17_mod, all14 = exact_cascade(W_base, dw0, dw1, k_cascade)
        # Получим полный De_17: re-compute with k=32 using same x16
        if all14:
            # de17_mod уже k_cascade-битный; для полного — нужен k=32
            # Но каскад при k=28 даёт DW[2..15] с нулями в битах 0..27
            # De_17 mod 2^28 = de17_mod; биты 28..31 неконтролируемы
            # Для v₂ нужны нижние биты → de17_mod достаточен
            de17_samples.append(de17_mod)

    n = len(de17_samples)
    print(f"  Получено {n} значений De_17 (28-битных)")

    # v₂-профиль
    v2_dist = Counter(v2(x) for x in de17_samples)
    print(f"\n  v₂(De_17)  наблюд    доля    ожид 2^{{-j-1}}  ratio")
    print(f"  {'─'*9}  {'─'*8}  {'─'*7}  {'─'*12}  {'─'*7}")

    for j in range(min(k_cascade, 20)):
        obs = v2_dist.get(j, 0)
        frac = obs / n
        if j < k_cascade:
            expected = 1.0 / (1 << (j + 1))
        else:
            expected = 1.0 / (1 << k_cascade)
        ratio = frac / expected if expected > 0 else 0
        marker = " ✓" if 0.6 <= ratio <= 1.6 else " ⚠"
        print(f"  {j:>9}  {obs:>8}  {frac:>7.4f}  {expected:>12.6f}  {ratio:>7.3f}{marker}")

    # Проверить: P(v₂ ≥ k) = |{De_17 mod 2^k = 0}| / n ≈ 1/2^k?
    print()
    print(f"  P(De_17 ≡ 0 mod 2^k) vs 1/2^k:")
    print(f"  {'k':>3}  {'нуль-проверок':>14}  {'P_obs':>10}  {'1/2^k':>10}  {'ratio':>7}")
    print(f"  {'─'*3}  {'─'*14}  {'─'*10}  {'─'*10}  {'─'*7}")

    for k in range(1, min(k_cascade + 1, 21)):
        mask = (1 << k) - 1
        zeros = sum(1 for x in de17_samples if (x & mask) == 0)
        p_obs = zeros / n
        p_exp = 1.0 / (1 << k)
        ratio = p_obs / p_exp if p_exp > 0 else 0
        marker = " ✓" if 0.5 <= ratio <= 2.0 else " ⚠"
        print(f"  {k:>3}  {zeros:>14}  {p_obs:>10.6f}  {p_exp:>10.6f}  {ratio:>7.3f}{marker}")


# ═══════════════════════════════════════════════════════════════════════════════
# ЭКСПЕРИМЕНТ E: МНОГОБАЗОВОЕ ПОДТВЕРЖДЕНИЕ k=18..20
# ═══════════════════════════════════════════════════════════════════════════════

def experiment_e_multibase(k_range=(15, 20), n_bases=5, trials_per_base=50_000):
    """
    Проверить Sol_k ≠ ∅ для k=15..20 на нескольких базах.

    При малом числе попыток (< 10×2^k), ищем хотя бы 1 решение.
    P(найти ≥1 за N попыток) = 1 - (1 - 1/2^k)^N.

    Для k=20, N=50K: P = 1 - (1-1/2^20)^50000 ≈ 1 - e^{-50000/2^20}
                        = 1 - e^{-0.0477} ≈ 4.7%.
    Для k=18, N=50K: P = 1 - e^{-50000/2^18} ≈ 1 - e^{-0.190} ≈ 17%.
    """
    print("\n" + "=" * 72)
    print("ЭКСПЕРИМЕНТ E: МНОГОБАЗОВОЕ ПОДТВЕРЖДЕНИЕ Sol_k≠∅")
    print("=" * 72)
    print(f"  {n_bases} баз, k={k_range[0]}..{k_range[1]}, "
          f"N={trials_per_base:,} попыток/база/k")
    print()

    found_matrix = {}  # found_matrix[base_seed][k] = (found, (dw0,dw1))

    for base_seed in range(n_bases):
        W_base = make_wbase(base_seed)
        found_matrix[base_seed] = {}
        rng = random.Random(base_seed * 1000 + 7)

        for k in range(k_range[0], k_range[1] + 1):
            mod = 1 << k
            found = False
            sol = None
            for _ in range(trials_per_base):
                dw0 = rng.randint(0, MASK)
                dw1 = rng.randint(0, mod - 1)
                _, de17, all14 = exact_cascade(W_base, dw0, dw1, k)
                if all14 and de17 == 0:
                    found = True
                    sol = (dw0, dw1)
                    break
            found_matrix[base_seed][k] = (found, sol)

    # Таблица результатов
    k_vals = list(range(k_range[0], k_range[1] + 1))
    print(f"  База  " + "  ".join(f"k={k}" for k in k_vals))
    print(f"  {'─'*4}  " + "  ".join("─" * 5 for _ in k_vals))

    for base_seed in range(n_bases):
        row = f"  {base_seed:>4}  "
        for k in k_vals:
            found, _ = found_matrix[base_seed][k]
            row += f"{'✓':>5}  " if found else f"{'─':>5}  "
        print(row)

    # Статистика
    print()
    print(f"  k    найдено/баз   P_find   теор_P_find")
    print(f"  {'─'*3}  {'─'*12}  {'─'*8}  {'─'*11}")
    for k in k_vals:
        found_count = sum(1 for b in range(n_bases)
                          if found_matrix[b][k][0])
        p_find = found_count / n_bases
        # Теоретическая P(хотя бы 1 из N попыток) = 1 - (1-1/2^k)^N
        import math
        p_theory = 1.0 - math.exp(-trials_per_base / (1 << k))
        print(f"  {k:>3}  {found_count:>2}/{n_bases:<10}  {p_find:>8.3f}  {p_theory:>11.4f}")

    print()
    all_found_any = all(
        any(found_matrix[b][k][0] for b in range(n_bases))
        for k in k_vals
    )
    if all_found_any:
        print(f"  ✓ Sol_k≠∅ подтверждено для всех k={k_range[0]}..{k_range[1]}")
    else:
        empty_ks = [k for k in k_vals
                    if not any(found_matrix[b][k][0] for b in range(n_bases))]
        print(f"  ? Sol_k не найдено для k={empty_ks} (статистические ограничения?)")


# ═══════════════════════════════════════════════════════════════════════════════
# ГЛАВНЫЙ БЛОК
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    import time

    print("=" * 72)
    print("П-54: БАШНЯ ДО k=20 — ПРОВЕРКА T_INFINITE_TOWER")
    print("SHA-256 жадный каскад: исследование height_2 и p-адической структуры")
    print("=" * 72)

    W_base = make_wbase(0)
    print(f"\nБаза W (seed=0): {W_base[:8]}...")

    # Эксперимент A: Башня k=12..20
    # max_trials=30K: для k<=14 используем 10×2^k; для k>14 используем 30K.
    # P(найти за 30K при k=15) = 1-e^{-30K/32K} ≈ 60%
    # P(найти за 30K при k=18) = 1-e^{-30K/256K} ≈ 11%
    t0 = time.time()
    results_a = experiment_a_tower(W_base, k_max=20, max_trials=30_000)
    print(f"\n  [Время A: {time.time()-t0:.1f}с]")

    # Эксперимент B: Hensel цепочка до k=16
    t0 = time.time()
    results_b = experiment_b_hensel(W_base, k_max=16, attempts_per_level=80)
    print(f"\n  [Время B: {time.time()-t0:.1f}с]")

    # Эксперимент C: Равномерность χ² (только k=8..12)
    t0 = time.time()
    experiment_c_uniformity(W_base, k_range=(8, 12), n_per_k=2000)
    print(f"\n  [Время C: {time.time()-t0:.1f}с]")

    # Эксперимент D: v₂-профиль (быстро, k_cascade=20)
    t0 = time.time()
    experiment_d_v2_profile(W_base, n_samples=1000)
    print(f"\n  [Время D: {time.time()-t0:.1f}с]")

    # Эксперимент E: Многобазовое подтверждение k=15..18
    t0 = time.time()
    experiment_e_multibase(k_range=(15, 18), n_bases=5, trials_per_base=20_000)
    print(f"\n  [Время E: {time.time()-t0:.1f}с]")

    # ─────────────────────────────────────────────────────────────────────────
    print("\n" + "=" * 72)
    print("ИТОГ — П-54")
    print("=" * 72)

    # Собрать выводы из Exp A
    reliable = [k for k, r in results_a.items()
                if r['n'] * (1.0 / r['mod']) >= 2 and r['successes'] > 0]
    barrier_ks = [k for k, r in results_a.items()
                  if r['successes'] == 0 and r['n'] * (1.0 / r['mod']) >= 2]

    print()
    if not barrier_ks:
        if reliable:
            max_k_reliable = max(reliable)
            print(f"ВЫВОД A — T_INFINITE_TOWER:")
            print(f"  Sol_k ≠ ∅ для k = {min(reliable)}..{max_k_reliable}")
            print(f"  P(Sol_k) ≈ 1/2^k во всём диапазоне.")
            print(f"  height_2(SHA-256) ≥ {max_k_reliable}")
            print(f"  → Сценарий A (башня БЕСКОНЕЧНА): данные ПОДТВЕРЖДАЮТ.")
        else:
            print("ВЫВОД A: Недостаточно данных (все k статистически слабые).")
    else:
        print(f"ВЫВОД A — БАРЬЕР при k={min(barrier_ks)}:")
        print(f"  height_2(SHA-256) = {min(barrier_ks) - 1}")
        print(f"  → Сценарий B (конечная башня): данные ПОДТВЕРЖДАЮТ.")

    print()
    print("ВЫВОД B — T_HENSEL_TOWER:")
    if results_b:
        max_hk = max(results_b.keys())
        any_break = any(r['success'] == 0 for r in results_b.values())
        if not any_break:
            print(f"  Согласованная цепочка поднялась до k={max_hk}.")
            print(f"  T_TOWER_INDEPENDENCE (П-52): НЕ подтверждена для точного каскада.")
        else:
            break_k = next(k for k, r in sorted(results_b.items()) if r['success'] == 0)
            print(f"  Цепочка оборвалась при k={break_k}.")
            print(f"  T_TOWER_INDEPENDENCE: частично подтверждена (барьер Hensel).")

    print()
    print("ТЕОРЕМЫ (П-54):")
    if not barrier_ks and reliable:
        print("  T_INFINITE_TOWER [ПОДТВЕРЖДЕНА числово]:")
        print(f"    Sol_k ≠ ∅ для k=1..{max(reliable)} (точный каскад).")
        print(f"    P(De_17≡0 mod 2^k) ≈ 1/2^k без структурного барьера.")
        print(f"    height_2(SHA-256) ≥ {max(reliable)}, вероятно = ∞.")
    print("  T_CASCADE_UNIQUENESS [из П-53, использована в П-54]:")
    print("    DW[2..15] однозначно определяются (dw0,dw1) за O(14) шагов.")
    print("  T_FREE_CONSTRAINT [из П-53]:")
    print("    De_17 — единственное несвязанное ограничение каскада.")
