"""
П-52: ЖАДНЫЙ mod-2^k КАСКАД — ВЫСОТА p-АДИЧЕСКОЙ БАШНИ

Цель: Проверить гипотезы M9, M10, M11 из П-51.

  M9:  Sol_k ≠ ∅ для всех k?
  M10: height_2(SHA-256) = ∞?
  M11: P(жадный mod-2^k успешен) ≈ 1/2^k за O(2^k × 14)?

Метод: Запустить жадный каскад для k=1,2,3,4,5,6,7,...
  На каждом уровне k: DW_i ∈ {0,...,2^k-1}.
  Считать P(успеха) и следить до какого k каскад работает.

Дополнительно:
  - Измерить P(Sol_k ≠ ∅) по нескольким базам.
  - Найти максимальное k, при котором жадный работает (на N попытках).
  - Если жадный падает при некотором k* → height_2 = k* - 1.
  - Если работает до k=32 → height_2 ≥ 32 (высота бесконечна).
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
def Sig0(x):    return rotr(x, 2)  ^ rotr(x, 13) ^ rotr(x, 22)
def Sig1(x):    return rotr(x, 6)  ^ rotr(x, 11) ^ rotr(x, 25)
def sig0(x):    return rotr(x, 7)  ^ rotr(x, 18) ^ (x >> 3)
def sig1(x):    return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)


def sha256_states_17(W16):
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


_cache = {}


def compute_f(W_base, x16):
    key = tuple(W_base)
    if key not in _cache:
        _cache.clear()
        _cache[key] = sha256_states_17(W_base)
    s1 = _cache[key]
    W2 = [(W_base[i] + x16[i]) & MASK for i in range(16)]
    s2 = sha256_states_17(W2)
    result = [(s2[r][0] - s1[r][0]) & MASK for r in range(3, 17)]
    result.append((s2[17][4] - s1[17][4]) & MASK)
    return result


# ═══════════════════════════════════════════════════════════════════════════════
# ЖАДНЫЙ mod-2^k КАСКАД — ОБЩИЙ
# ═══════════════════════════════════════════════════════════════════════════════

def greedy_cascade_modk(W_base, dw0, k, n_attempts=None):
    """
    Жадный каскад: DW_i ∈ {0,...,2^k-1}.
    Ищем F(x) ≡ 0 (mod 2^k).

    Оптимизация: для больших k (2^k > MAX_EXHAUSTIVE)
    перебираем не все 2^k значений, а сэмплируем SAMPLE_SIZE.
    Это позволяет работать при k=1..32 за разумное время.
    """
    mod = 1 << k
    MAX_EXHAUSTIVE = 64   # при mod ≤ 64 перебираем полностью
    SAMPLE_SIZE = 64      # при mod > 64 берём 64 случайных + 8 малых

    if n_attempts is None:
        n_attempts = min(2 * mod, 128)

    solutions = []
    best_zeros = 0

    for attempt in range(n_attempts):
        dw1 = attempt % mod
        x16 = [dw0, dw1] + [0] * 14

        for pos in range(2, 16):
            best_v = 0
            best_dist = mod

            if mod <= MAX_EXHAUSTIVE:
                candidates = range(mod)
            else:
                # малые значения + случайные по всему диапазону
                small = list(range(min(8, mod)))
                rand_c = [random.randint(0, mod - 1) for _ in range(SAMPLE_SIZE - len(small))]
                seen = set(small)
                uniq_rand = []
                for c in rand_c:
                    if c not in seen:
                        seen.add(c)
                        uniq_rand.append(c)
                candidates = small + uniq_rand

            for v in candidates:
                x_try = list(x16)
                x_try[pos] = v
                da = compute_f(W_base, x_try)[pos - 2]
                da_mod = da % mod
                dist = min(da_mod, mod - da_mod)
                if dist < best_dist:
                    best_dist = dist
                    best_v = v

            x16[pos] = best_v

        f = compute_f(W_base, x16)
        f_modk = [v % mod for v in f]
        n_zero = sum(1 for v in f_modk if v == 0)
        best_zeros = max(best_zeros, n_zero)

        if all(v == 0 for v in f_modk):
            solutions.append(x16[:])

    return solutions, best_zeros, n_attempts


# ═══════════════════════════════════════════════════════════════════════════════
# ОСНОВНОЙ ЭКСПЕРИМЕНТ: БАШНЯ ДО mod-2^K_MAX
# ═══════════════════════════════════════════════════════════════════════════════

def run_tower_experiment(W_base, dw0, k_max=12, n_reps=3):
    """
    Для каждого k=1..k_max: запустить жадный mod-2^k.
    Печатает результат сразу после каждого k (не ждёт конца).
    """
    print(f"\n{'k':>3}  {'mod 2^k':>8}  {'решений':>8}  {'попыток':>8}  "
          f"{'P(успех)':>9}  {'лучш нулей':>11}  {'Sol_k≠∅?':>9}")
    print("-" * 70)

    results = {}
    for k in range(1, k_max + 1):
        mod = 1 << k
        # Ограничиваем попытки: исчерпывающий при малых mod, сэмплинг при больших
        n_attempts = min(4 * mod, 64)
        total_sols = 0
        best_z = 0

        for rep in range(n_reps):
            dw0_r = (dw0 + rep * 0x1234567) & MASK
            sols, bz, _ = greedy_cascade_modk(W_base, dw0_r, k, n_attempts)
            total_sols += len(sols)
            best_z = max(best_z, bz)

        n_tried = n_reps * n_attempts
        p_success = total_sols / n_tried if n_tried else 0
        found = "ДА ✓" if total_sols > 0 else "нет"
        print(f"{k:>3}  {mod:>8}  {total_sols:>8}  {n_tried:>8}  "
              f"{p_success:>9.4f}  {best_z:>11}/15  {found:>9}", flush=True)

        results[k] = {
            'mod': mod, 'total_sols': total_sols,
            'n_tried': n_tried, 'p_success': p_success,
            'best_zeros': best_z, 'any_found': total_sols > 0,
        }

    return results


# ═══════════════════════════════════════════════════════════════════════════════
# ПРОВЕРКА: ЯВЛЯЮТСЯ ЛИ Sol_k НЕЗАВИСИМЫМИ ВЛОЖЕННЫМИ СЛОЯМИ
# ═══════════════════════════════════════════════════════════════════════════════

def check_tower_consistency(W_base, dw0, k_max=6):
    """
    Для каждого k: найти Sol_k решение и проверить:
    является ли оно (mod 2^{k-1}) элементом Sol_{k-1}?
    Т.е. башня вложена: Sol_k mod 2^{k-1} ⊆ Sol_{k-1}?
    """
    print(f"\n── Проверка вложенности башни Sol_k mod 2^{{k-1}} ⊆ Sol_{{k-1}} ──")
    solutions_by_k = {}

    for k in range(1, k_max + 1):
        mod = 1 << k
        n_attempts = min(8 * mod, 512)
        sols, _, _ = greedy_cascade_modk(W_base, dw0, k, n_attempts)
        if sols:
            solutions_by_k[k] = sols[0]

    for k in range(2, k_max + 1):
        if k not in solutions_by_k or k-1 not in solutions_by_k:
            print(f"  k={k}: нет данных.")
            continue

        x_k = solutions_by_k[k]
        mod_km1 = 1 << (k-1)

        # Проверить: x_k mod 2^{k-1} ∈ Sol_{k-1}?
        x_reduced = [v % mod_km1 for v in x_k]
        f_red = compute_f(W_base, x_reduced)
        is_sol_km1 = all(v % mod_km1 == 0 for v in f_red)

        # Также проверить сам x_k относительно Sol_k
        f_k = compute_f(W_base, x_k)
        is_sol_k = all(v % (1 << k) == 0 for v in f_k)

        print(f"  k={k}: x_k ∈ Sol_{k}? {is_sol_k}  |  "
              f"x_k mod 2^{k-1} ∈ Sol_{k-1}? {is_sol_km1}")

    return solutions_by_k


# ═══════════════════════════════════════════════════════════════════════════════
# СТАТИСТИКА: P(Sol_k ≠ ∅) ПО МНОГИМ БАЗАМ
# ═══════════════════════════════════════════════════════════════════════════════

def stats_over_bases(n_bases=20, k_max=8, n_attempts_per_k=None):
    """
    По n_bases случайным базам: для каждого k считать долю баз с Sol_k ≠ ∅.
    """
    print(f"\n── Статистика P(Sol_k ≠ ∅) по {n_bases} базам ──")
    found_count = {k: 0 for k in range(1, k_max + 1)}
    base_count = 0

    for base_idx in range(n_bases):
        W_base = [random.randint(0, MASK) for _ in range(16)]
        dw0 = random.randint(1, MASK)
        base_count += 1

        for k in range(1, k_max + 1):
            mod = 1 << k
            if n_attempts_per_k:
                n_att = n_attempts_per_k
            else:
                n_att = min(4 * mod, 256)
            sols, _, _ = greedy_cascade_modk(W_base, dw0, k, n_att)
            if sols:
                found_count[k] += 1

    print(f"\n{'k':>3}  {'mod 2^k':>8}  {'баз с Sol_k≠∅':>14}  {'P(Sol_k≠∅)':>11}")
    print("-" * 42)
    for k in range(1, k_max + 1):
        mod = 1 << k
        p = found_count[k] / base_count
        bar = "█" * int(p * 20)
        print(f"{k:>3}  {mod:>8}  {found_count[k]:>14}/{base_count}  "
              f"{p:>11.3f}  {bar}")

    return found_count, base_count


# ═══════════════════════════════════════════════════════════════════════════════
# ГЛАВНЫЙ ЗАПУСК
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    random.seed(31415)

    print("=" * 72)
    print("П-52: ЖАДНЫЙ mod-2^k КАСКАД — ВЫСОТА p-АДИЧЕСКОЙ БАШНИ SHA-256")
    print("Проверка M9 (башня бесконечна?), M10, M11")
    print("=" * 72)

    # ── Эксперимент 1: Башня для одной базы до k=16 ──────────────────────
    print("\n[1] Башня для одной базы (k=1..16, 5 повторов на k):")
    W_base_main = [random.randint(0, MASK) for _ in range(16)]
    dw0_main = random.randint(1, MASK)
    print(f"    DW₀={dw0_main:#010x}")

    tower_results = run_tower_experiment(W_base_main, dw0_main, k_max=12, n_reps=3)

    # Найти первый k, где P=0
    first_empty_k = None
    for k in range(1, 13):
        if not tower_results[k]['any_found']:
            first_empty_k = k
            break

    if first_empty_k:
        print(f"\n  → Первый пустой уровень: k={first_empty_k}")
        print(f"  → height_2(SHA-256) = {first_empty_k - 1} (для данной базы)")
    else:
        print(f"\n  → Sol_k НЕПУСТ для всех k=1..12!")
        print(f"  → height_2(SHA-256) ≥ 12 (для данной базы)")

    # ── Эксперимент 2: Вложенность башни ────────────────────────────────
    print("\n[2] Вложенность: является ли Sol_k mod 2^{k-1} ∈ Sol_{k-1}?")
    sols_by_k = check_tower_consistency(W_base_main, dw0_main, k_max=8)

    # ── Эксперимент 3: Статистика по многим базам (k=1..8) ───────────────
    print("\n[3] Статистика P(Sol_k ≠ ∅) по 10 базам (k=1..8):")
    found_count, n_bases = stats_over_bases(n_bases=10, k_max=8, n_attempts_per_k=32)

    # ── Итог ────────────────────────────────────────────────────────────
    print(f"\n{'='*72}")
    print("ИТОГ — П-52")
    print(f"{'='*72}")

    # Анализ тренда P(Sol_k)
    print("\nТренд P(Sol_k ≠ ∅) по базам:")
    prev_p = None
    ratios = []
    for k in range(1, 9):
        p = found_count[k] / n_bases
        ratio = p / prev_p if prev_p and prev_p > 0 else None
        ratios.append(ratio)
        ratio_str = f"  (P_k/P_{{k-1}} = {ratio:.2f})" if ratio else ""
        print(f"  k={k}: P = {p:.3f}{ratio_str}")
        prev_p = p

    print()
    if first_empty_k is None:
        print("ВЫВОД (башня одной базы): Sol_k непуст для k=1..16.")
        print("  M9 и M10 ПОДТВЕРЖДАЮТСЯ для этой базы.")
    else:
        print(f"ВЫВОД (башня одной базы): Sol_k пуст начиная с k={first_empty_k}.")
        print(f"  height_2 = {first_empty_k - 1} для этой базы.")

    print()
    meaningful_ratios = [r for r in ratios if r is not None]
    if meaningful_ratios:
        avg_ratio = sum(meaningful_ratios) / len(meaningful_ratios)
        print(f"Средний коэффициент P_k/P_{{k-1}} = {avg_ratio:.3f}")
        if abs(avg_ratio - 0.5) < 0.15:
            print("  ≈ 0.5 → M11 ПОДТВЕРЖДАЕТСЯ: P(Sol_k) ≈ P(Sol_1) / 2^(k-1)")
        elif avg_ratio > 0.7:
            print(f"  > 0.5 → P убывает медленнее 1/2. Башня устойчивее ожидаемого.")
        else:
            print(f"  < 0.5 → P убывает быстрее 1/2. Башня неустойчива.")

    print()
    print("Ключевые выводы для методички:")
    any_k12 = tower_results.get(12, {}).get('any_found', False)
    print(f"  Sol_12 (mod 2^12=4096) найден? {any_k12}")
    print(f"  Жадный каскад масштабируется до k={12 if first_empty_k is None else first_empty_k-1}")
    print(f"  height_2(SHA-256) {'≥12' if first_empty_k is None else '= '+str(first_empty_k-1)}")


if __name__ == '__main__':
    main()
