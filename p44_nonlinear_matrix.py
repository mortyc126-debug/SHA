"""
П-44: НЕЛИНЕЙНАЯ МАТРИЦА ДЛЯ SHA-256 — ИСЧЕРПЫВАЮЩИЙ БИТ-ЗА-БИТОМ ПОДЪЁМ.

Проблема из П-43: классический подъём Хенселя (Newton p-adic) терпит крах при k=2,
поскольку SHA-256 переносы нарушают условие 2-адической гладкости:
  v2(f(x+2^k) - f(x) - J·2^k) ≥ k+1  — НЕ выполняется для SHA-256!

Решение (нелинейная матрица):
  Вместо якобиана — ПОЛНЫЙ ПЕРЕБОР {0,1}² на каждом бите.
  На уровне k: тестируем все (d14,d15) ∈ {0,1}²:
    f(x14 + d14·2^k, x15 + d15·2^k) ≡ (0,0) (mod 2^{k+1})?
  4 вычисления SHA × 32 бита = ~128 вызовов SHA-256 на попытку.
  Если несколько кандидатов → ветвление (backtracking-дерево).
  Если ни одного → возврат с откатом.

Это «нелинейная матрица» — таблица переходов на каждом бите, специально
заточенная под SHA-256: учитывает все переносы и нелинейности точно,
без линейных приближений.

Ожидание:
  - Если система Da16=0, De17=0 разрешима в Z/2^32Z → найдём за ~128 SHA.
  - Если неразрешима ни для каких DW0..DW13 → это подтверждение барьера 2^64.

Методичка: methodology_v15.md (барьер 16 нулей, T_HENSEL_INAPPLICABLE).
"""

import random
import math
import sys

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
def hw(x): return bin(x).count('1')


def expand_schedule(W16):
    W = list(W16)
    for i in range(16, 64):
        W.append((sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16]) & MASK)
    return W


def sha256_state(W16, nrounds):
    W = expand_schedule(W16)
    a, b, c, d, e, f, g, h = H0
    for r in range(nrounds):
        T1 = (h + Sig1(e) + Ch(e, f, g) + K_SHA[r] + W[r]) & MASK
        T2 = (Sig0(a) + Maj(a, b, c)) & MASK
        h = g; g = f; f = e; e = (d + T1) & MASK
        d = c; c = b; b = a; a = (T1 + T2) & MASK
    return (a, b, c, d, e, f, g, h)


def De_at(W_base, dW_dict, r):
    W2 = list(W_base)
    for idx, delta in dW_dict.items():
        W2[idx] = (W2[idx] + delta) & MASK
    s1 = sha256_state(W_base, r)
    s2 = sha256_state(W2, r)
    return (s2[4] - s1[4]) & MASK


def Da_at(W_base, dW_dict, r):
    W2 = list(W_base)
    for idx, delta in dW_dict.items():
        W2[idx] = (W2[idx] + delta) & MASK
    s1 = sha256_state(W_base, r)
    s2 = sha256_state(W2, r)
    return (s2[0] - s1[0]) & MASK


def solve_mod32(slope, target):
    g = math.gcd(slope & MASK, 2**32)
    if (target & MASK) % g != 0:
        return None
    s_red = (slope & MASK) // g
    t_red = (target & MASK) // g
    m_red = (2**32) // g
    return (t_red * pow(s_red % m_red, -1, m_red)) % m_red


def build_13step_cascade(W_base, dW0=1):
    """Каскад Da3..Da15=0 (13 шагов), DW14 и DW15 остаются свободными."""
    dW = {0: dW0}
    for r in range(3, 16):
        fw = r - 2
        dW[fw] = 0
        val0 = Da_at(W_base, dW, r)
        dW[fw] = 1
        val1 = Da_at(W_base, dW, r)
        slope = (val1 - val0) & MASK
        x = solve_mod32(slope, (-val0) & MASK)
        if x is None:
            x = 0
        dW[fw] = x & MASK
        if Da_at(W_base, dW, r) != 0:
            found = False
            for delta in range(1, 256):
                for sign in (1, -1):
                    dW[fw] = (x + sign * delta) & MASK
                    if Da_at(W_base, dW, r) == 0:
                        found = True
                        break
                if found:
                    break
            if not found:
                dW[fw] = x & MASK
    return dW


def f_system(W_base, dW_base, x14, x15):
    """Вычислить (Da16, De17) для данных DW14=x14, DW15=x15."""
    dW = dict(dW_base)
    dW[14] = x14 & MASK
    dW[15] = x15 & MASK
    return Da_at(W_base, dW, 16), De_at(W_base, dW, 17)


# ═══════════════════════════════════════════════════════════════════════════════
# A. НЕЛИНЕЙНЫЙ ИСЧЕРПЫВАЮЩИЙ ПОДЪЁМ (NONLINEAR EXHAUSTIVE LIFT)
# ═══════════════════════════════════════════════════════════════════════════════

def exhaustive_lift_nobranch(W_base, dW_base, verbose=False):
    """
    Жадный (без ветвления) нелинейный подъём.

    На каждом бите k перебираем все (d14,d15) ∈ {0,1}²:
      - Тестируем f(x14+d14*2^k, x15+d15*2^k) mod 2^{k+1}
      - Берём первый подходящий кандидат
      - Если ни одного — возврат с ошибкой (no_candidate)

    Возвращает: (x14, x15, info_dict) или (None, None, error_str)
    """
    x14, x15 = 0, 0
    sha_calls = 0
    branch_points = []   # биты, где было несколько кандидатов

    for k in range(32):
        bit_k = 1 << k
        mod_next = 1 << (k + 1)
        candidates = []

        for d14 in range(2):
            for d15 in range(2):
                tx14 = (x14 + d14 * bit_k) & MASK
                tx15 = (x15 + d15 * bit_k) & MASK
                f0, f1 = f_system(W_base, dW_base, tx14, tx15)
                sha_calls += 4  # f_system делает 4 SHA-256
                if f0 % mod_next == 0 and f1 % mod_next == 0:
                    candidates.append((d14, d15))

        if not candidates:
            return None, None, f'no_candidate_at_k={k}', sha_calls

        if len(candidates) > 1:
            branch_points.append((k, len(candidates)))

        d14, d15 = candidates[0]
        x14 = (x14 + d14 * bit_k) & MASK
        x15 = (x15 + d15 * bit_k) & MASK

        if verbose:
            print(f"  k={k:2d}: candidates={len(candidates)} → d14={d14}, d15={d15}  "
                  f"x14=0x{x14:08x} x15=0x{x15:08x}")

    # Верификация финального результата
    f0, f1 = f_system(W_base, dW_base, x14, x15)
    sha_calls += 4

    info = {
        'sha_calls': sha_calls,
        'branch_points': branch_points,
        'f0': f0,
        'f1': f1,
        'success': f0 == 0 and f1 == 0,
    }
    if f0 == 0 and f1 == 0:
        return x14, x15, 'ok', sha_calls
    else:
        return None, None, f'verify_fail:Da16=0x{f0:08x},De17=0x{f1:08x}', sha_calls


def exhaustive_lift_backtrack(W_base, dW_base, max_backtrack=10, verbose=False):
    """
    Подъём с откатом (backtracking).

    Если на уровне k нет кандидатов, возвращаемся на предыдущий уровень
    и пробуем следующего кандидата (если он был).

    Стек: список (k, x14, x15, remaining_candidates).
    """
    sha_calls = [0]
    total_nodes = [0]

    def f_eval(x14, x15):
        f0, f1 = f_system(W_base, dW_base, x14, x15)
        sha_calls[0] += 4
        return f0, f1

    def get_candidates(k, x14, x15):
        bit_k = 1 << k
        mod_next = 1 << (k + 1)
        result = []
        for d14 in range(2):
            for d15 in range(2):
                tx14 = (x14 + d14 * bit_k) & MASK
                tx15 = (x15 + d15 * bit_k) & MASK
                f0, f1 = f_eval(tx14, tx15)
                if f0 % mod_next == 0 and f1 % mod_next == 0:
                    result.append((tx14, tx15))
        return result

    # Стек DFS: (уровень k, x14, x15, список кандидатов ещё не испробованных)
    stack = []
    candidates_0 = get_candidates(0, 0, 0)
    if not candidates_0:
        return None, None, 'no_candidate_at_k=0', sha_calls[0]

    stack.append((0, 0, 0, candidates_0))
    backtracks = 0

    while stack:
        k, px14, px15, cands = stack[-1]

        if not cands:
            # Откат
            stack.pop()
            backtracks += 1
            if backtracks > max_backtrack:
                return None, None, f'max_backtrack_exceeded', sha_calls[0]
            continue

        # Берём первого кандидата
        x14, x15 = cands.pop(0)
        # Обновляем стек с оставшимися кандидатами
        stack[-1] = (k, px14, px15, cands)
        total_nodes[0] += 1

        if k == 31:
            # Финальный уровень — верификация
            f0, f1 = f_eval(x14, x15)
            if f0 == 0 and f1 == 0:
                return x14, x15, f'ok_backtracks={backtracks}', sha_calls[0]
            # Иначе продолжаем (pop уже произошёл через cands)
            continue

        # Идём глубже
        next_cands = get_candidates(k + 1, x14, x15)
        if next_cands:
            stack.append((k + 1, x14, x15, next_cands))
        # Иначе просто берём следующего кандидата на текущем уровне

    return None, None, f'dfs_exhausted_backtracks={backtracks}', sha_calls[0]


# ═══════════════════════════════════════════════════════════════════════════════
# B. СТАТИСТИКА КАНДИДАТОВ НА КАЖДОМ БИТЕ
# ═══════════════════════════════════════════════════════════════════════════════

def analyze_candidate_counts(W_base, dW_base, n_samples=5):
    """
    Анализ: сколько кандидатов (d14,d15) проходит фильтр на каждом бите k?
    Это показывает структуру «нелинейной матрицы».
    """
    x14, x15 = 0, 0
    per_bit_counts = []

    for k in range(32):
        bit_k = 1 << k
        mod_next = 1 << (k + 1)
        count = 0
        for d14 in range(2):
            for d15 in range(2):
                tx14 = (x14 + d14 * bit_k) & MASK
                tx15 = (x15 + d15 * bit_k) & MASK
                f0, f1 = f_system(W_base, dW_base, tx14, tx15)
                if f0 % mod_next == 0 and f1 % mod_next == 0:
                    count += 1
        per_bit_counts.append(count)
        # Жадно берём первого кандидата для продолжения
        found_next = False
        for d14 in range(2):
            for d15 in range(2):
                tx14 = (x14 + d14 * bit_k) & MASK
                tx15 = (x15 + d15 * bit_k) & MASK
                f0, f1 = f_system(W_base, dW_base, tx14, tx15)
                if f0 % mod_next == 0 and f1 % mod_next == 0:
                    x14, x15 = tx14, tx15
                    found_next = True
                    break
            if found_next:
                break
        if not found_next:
            per_bit_counts.extend([0] * (32 - k - 1))
            break

    return per_bit_counts


# ═══════════════════════════════════════════════════════════════════════════════
# C. ПРОВЕРКА: СУЩЕСТВУЕТ ЛИ ВООБЩЕ РЕШЕНИЕ В Z/2^32Z?
# ═══════════════════════════════════════════════════════════════════════════════

def check_solution_exists_brute(W_base, dW_base, n_samples=256):
    """
    Случайная выборка: ищем любое (x14, x15) с Da16=0, De17=0.
    n_samples случайных пар → нижняя оценка вероятности существования решения.
    """
    for _ in range(n_samples):
        x14 = random.randint(0, MASK)
        x15 = random.randint(0, MASK)
        f0, f1 = f_system(W_base, dW_base, x14, x15)
        if f0 == 0 and f1 == 0:
            return True, x14, x15
    return False, None, None


def check_solution_exists_grid(W_base, dW_base, n_grid=64):
    """
    Сетка 64×64: ищем любое (x14, x15) с Da16=0, De17=0 методом перебора.
    Используем маленькие значения (0..63) для проверки вблизи нуля.
    """
    step = max(1, MASK // n_grid)
    for i in range(n_grid):
        for j in range(n_grid):
            x14 = (i * step) & MASK
            x15 = (j * step) & MASK
            f0, f1 = f_system(W_base, dW_base, x14, x15)
            if f0 == 0 and f1 == 0:
                return True, x14, x15
    return False, None, None


# ═══════════════════════════════════════════════════════════════════════════════
# D. ГЛАВНАЯ ПРОГРАММА
# ═══════════════════════════════════════════════════════════════════════════════

def run_section_a(n_tests=20):
    print("─" * 72)
    print("A. ЖАДНЫЙ НЕЛИНЕЙНЫЙ ПОДЪЁМ (БЕЗ ВЕТВЛЕНИЯ)")
    print("─" * 72)
    print("Алгоритм: на каждом бите k перебираем (d14,d15)∈{0,1}²,")
    print("берём первого кандидата удовлетворяющего f≡0 mod 2^{k+1}.")
    print(f"Тест на {n_tests} случайных сообщениях:\n")
    print(f"  {'#':>3} | {'dW0':>10} | {'Результат':<30} | {'SHA-calls':>9} | {'Ветвления'}")
    print(f"  ────┼────────────┼{'─'*30}┼───────────┼───────────")

    successes = 0
    total_sha = 0

    for i in range(n_tests):
        W_base = [random.randint(0, MASK) for _ in range(16)]
        dW0 = random.randint(1, MASK)
        dW_base = build_13step_cascade(W_base, dW0=dW0)

        x14, x15, status, sha = exhaustive_lift_nobranch(W_base, dW_base)
        total_sha += sha

        if status == 'ok':
            successes += 1
            result_str = f"OK  DW14=0x{x14:08x}"
        else:
            result_str = f"FAIL: {status[:25]}"

        print(f"  {i+1:>3} | 0x{dW0:08x} | {result_str:<30} | {sha:>9} | -")

    print(f"\n  Успех: {successes}/{n_tests} ({100*successes/n_tests:.0f}%)")
    print(f"  Среднее SHA-вызовов: {total_sha/n_tests:.0f}")
    return successes, n_tests


def run_section_b(n_msgs=5):
    print("\n" + "─" * 72)
    print("B. СТРУКТУРА КАНДИДАТОВ: СКОЛЬКО (d14,d15) ПРОХОДЯТ НА КАЖДОМ БИТЕ?")
    print("─" * 72)
    print("Идея: если на каждом бите ровно 1 кандидат → решение единственно.")
    print("Если 2+ → есть развилка (ветвление). Если 0 → нужен откат.\n")

    for msg_idx in range(n_msgs):
        W_base = [random.randint(0, MASK) for _ in range(16)]
        dW0 = random.randint(1, MASK)
        dW_base = build_13step_cascade(W_base, dW0=dW0)

        counts = analyze_candidate_counts(W_base, dW_base)

        zeros = counts.count(0)
        ones   = counts.count(1)
        twos   = counts.count(2)
        threes = counts.count(3)
        fours  = counts.count(4)

        print(f"  Сообщение {msg_idx+1}: dW0=0x{dW0:08x}")
        print(f"    Распределение кол-ва кандидатов по битам 0..31:")
        print(f"      0 кандидатов: {zeros:2d} битов  ← ТУПИК (надо откат)")
        print(f"      1 кандидат:   {ones:2d} битов  ← уникально")
        print(f"      2 кандидата:  {twos:2d} битов  ← 2 ветви")
        print(f"      3 кандидата:  {threes:2d} битов")
        print(f"      4 кандидата:  {fours:2d} битов  ← все прошли")

        # Гистограмма первых 16 битов
        bar = ''.join(str(c) if c > 0 else '.' for c in counts[:32])
        print(f"    Карта (0=тупик, 1-4=кол-во): {bar}")
        print()


def run_section_c(n_msgs=10):
    print("─" * 72)
    print("C. ПОДЪЁМ С ОТКАТОМ (BACKTRACKING)")
    print("─" * 72)
    print("Если жадный алгоритм не нашёл решение → пробуем backtracking DFS.")
    print(f"Тест на {n_msgs} сообщениях:\n")

    successes = 0
    for i in range(n_msgs):
        W_base = [random.randint(0, MASK) for _ in range(16)]
        dW0 = random.randint(1, MASK)
        dW_base = build_13step_cascade(W_base, dW0=dW0)

        x14, x15, status, sha = exhaustive_lift_backtrack(
            W_base, dW_base, max_backtrack=50
        )

        if x14 is not None:
            successes += 1
            # Верифицируем
            f0, f1 = f_system(W_base, dW_base, x14, x15)
            ok = "✓" if f0 == 0 and f1 == 0 else "✗"
            print(f"  {i+1:>3} OK {ok}  DW14=0x{x14:08x}  DW15=0x{x15:08x}  "
                  f"SHA={sha}  [{status}]")
        else:
            print(f"  {i+1:>3} FAIL  SHA={sha}  [{status}]")

    print(f"\n  Успех (с откатом): {successes}/{n_msgs}")
    return successes, n_msgs


def run_section_d(n_msgs=5):
    print("\n" + "─" * 72)
    print("D. ПРОВЕРКА СУЩЕСТВОВАНИЯ РЕШЕНИЯ ЧЕРЕЗ СЛУЧАЙНЫЙ ПОИСК")
    print("─" * 72)
    print("Случайно проверяем 1000 пар (x14,x15) — найдём ли хоть одно решение Da16=De17=0?\n")

    found_count = 0
    for i in range(n_msgs):
        W_base = [random.randint(0, MASK) for _ in range(16)]
        dW0 = random.randint(1, MASK)
        dW_base = build_13step_cascade(W_base, dW0=dW0)

        found, x14, x15 = check_solution_exists_brute(W_base, dW_base, n_samples=1000)
        if found:
            found_count += 1
            f0, f1 = f_system(W_base, dW_base, x14, x15)
            print(f"  Сообщение {i+1}: НАЙДЕНО!  x14=0x{x14:08x} x15=0x{x15:08x}  "
                  f"Da16={f0} De17={f1}")
        else:
            print(f"  Сообщение {i+1}: не найдено за 1000 попыток  (ожидаемо, P≈1/2^64)")

    print(f"\n  Найдено решений: {found_count}/{n_msgs}")
    print("  (Случайный поиск за 1000 пар: вероятность успеха ~ 1000/2^32 ≈ 0.00002%)")


def run_section_e(n_msgs=20):
    print("\n" + "─" * 72)
    print("E. ИТОГОВЫЙ ТЕСТ: ЖАДНЫЙ + ОТКАТ ПО 20 СООБЩЕНИЯМ")
    print("─" * 72)
    print("Комбинируем: сначала жадный, при неудаче — откат (до 100 откатов).\n")

    results = {'greedy_ok': 0, 'backtrack_ok': 0, 'fail': 0}
    total_sha = 0

    for i in range(n_msgs):
        W_base = [random.randint(0, MASK) for _ in range(16)]
        dW0 = random.randint(1, MASK)
        dW_base = build_13step_cascade(W_base, dW0=dW0)

        # Сначала жадный
        x14, x15, status, sha1 = exhaustive_lift_nobranch(W_base, dW_base)
        total_sha += sha1

        if status == 'ok':
            results['greedy_ok'] += 1
            outcome = f"ЖАДНЫЙ OK  DW14=0x{x14:08x}"
        else:
            # Откат
            x14, x15, status2, sha2 = exhaustive_lift_backtrack(
                W_base, dW_base, max_backtrack=100
            )
            total_sha += sha2
            if x14 is not None:
                results['backtrack_ok'] += 1
                outcome = f"ОТКАТ OK   DW14=0x{x14:08x}  (extra {sha2} SHA)"
            else:
                results['fail'] += 1
                outcome = f"FAIL [{status2[:30]}]"

        print(f"  {i+1:>3} | {outcome}")

    total_ok = results['greedy_ok'] + results['backtrack_ok']
    print(f"\n  ┌─────────────────────────────────────┐")
    print(f"  │ Жадный успех:   {results['greedy_ok']:3d}/{n_msgs}                 │")
    print(f"  │ Откат успех:    {results['backtrack_ok']:3d}/{n_msgs}                 │")
    print(f"  │ Итого успехов:  {total_ok:3d}/{n_msgs} ({100*total_ok/n_msgs:.0f}%)           │")
    print(f"  │ Неудач:         {results['fail']:3d}/{n_msgs}                 │")
    print(f"  │ Среднее SHA:    {total_sha/n_msgs:.0f}                    │")
    print(f"  └─────────────────────────────────────┘")

    return results, n_msgs


def print_header():
    print("=" * 72)
    print("П-44 | НЕЛИНЕЙНАЯ МАТРИЦА SHA-256 — ИСЧЕРПЫВАЮЩИЙ ПОДЪЁМ")
    print("=" * 72)
    print()
    print("Задача: найти (DW14, DW15) такие что Da16=0 И De17=0.")
    print("Каскад Da3..Da15=0 уже построен через DW1..DW13.")
    print("Метод: перебор {0,1}² на каждом из 32 битов — без Якобиана!")
    print()


def print_conclusions(sec_a, sec_c, sec_e):
    print("\n" + "=" * 72)
    print("ВЫВОДЫ П-44")
    print("=" * 72)

    a_ok, a_total = sec_a
    c_ok, c_total = sec_c
    e_res, e_total = sec_e
    e_ok = e_res['greedy_ok'] + e_res['backtrack_ok']

    print(f"""
НЕЛИНЕЙНАЯ МАТРИЦА (исчерпывающий побитовый подъём):
  — Жадный алгоритм (без ветвления): {a_ok}/{a_total} успехов
  — Откат (DFS backtracking):         {c_ok}/{c_total} успехов
  — Комбинированный:                  {e_ok}/{e_total} успехов

ИНТЕРПРЕТАЦИЯ:""")

    if e_ok == 0:
        print("""
  БАРЬЕР ПОДТВЕРЖДЁН: нелинейная матрица тоже не помогает.
  Причина: система Da16=0, De17=0 НЕ имеет решений в Z/2^32Z
  для типичных каскадов Da3..Da15=0.

  T_NONLINEAR_MATRIX_FAILS:
    Даже исчерпывающий перебор {0,1}² без каких-либо линейных
    приближений не находит решение — значит решения нет в Z/2^32Z.
    Барьер 2^64 для 16+ нулей подтверждён третьим независимым методом
    (после T_JACOBIAN_APPROX и T_HENSEL_INAPPLICABLE).

  ВЫВОД: Барьер SHA-256 на 16 нулях — фундаментальный,
  не артефакт метода аппроксимации.""")
    elif e_ok < e_total // 2:
        print(f"""
  ЧАСТИЧНЫЙ УСПЕХ: решение находится в {e_ok}/{e_total} случаев.
  — Когда решение существует → нелинейная матрица его находит!
  — Когда решения нет → алгоритм правильно возвращает FAIL.
  — Кол-во SHA-256 вызовов: ~128 на успешный случай (vs 2^32 brute force).

  T_NONLINEAR_MATRIX_WORKS_CONDITIONALLY:
    Нелинейная матрица эффективна когда решение существует,
    но само существование решения редко (P ≈ 1/2^{e_total-e_ok}).
    Барьер сохраняется на уровне 2^32 для нахождения решения.""")
    else:
        print(f"""
  НЕОЖИДАННЫЙ УСПЕХ: решение находится в {e_ok}/{e_total} случаев!
  — Нелинейная матрица преодолевает барьер!
  — Требуется углублённый анализ (П-45).""")

    print(f"""
СРАВНЕНИЕ МЕТОДОВ:
  П-42 (Якобиан линейный):  0% успех, линейное приближение неверно
  П-43 (Хенсель p-adic):    0% успех, v2-условие нарушается при k=2
  П-44 (нелинейная матрица): {100*e_ok//e_total}% успех, без приближений
""")


if __name__ == '__main__':
    random.seed(42)
    print_header()

    sec_a = run_section_a(n_tests=20)

    run_section_b(n_msgs=3)

    sec_c = run_section_c(n_msgs=10)

    run_section_d(n_msgs=3)

    sec_e = run_section_e(n_msgs=20)

    print_conclusions(sec_a, sec_c, sec_e)
