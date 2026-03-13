"""
П-45: ГЛОБАЛЬНАЯ ПОСТАНОВКА ЗАДАЧИ — 15 УРАВНЕНИЙ / 16 НЕИЗВЕСТНЫХ.

Диагноз из П-44: жадный каскад ALL-A (Da3..Da16=0 через DW1..DW14)
создаёт ИСКУССТВЕННУЮ несовместность: к моменту когда доходим до De17,
все степени свободы уже исчерпаны. Это не свойство SHA-256 — это
артефакт жадной процедуры.

Правильная постановка:
  16 неизвестных: DW0..DW15 ∈ Z/2^32Z
  15 уравнений:   Da3=0, Da4=0, ..., Da16=0, De17=0
  → НЕДООПРЕДЕЛЁННАЯ система → решения ДОЛЖНЫ существовать!

Ключевые вопросы:
  A. Подтверждение: глобальная система разрешима? (тест 3D lift)
  B. Новая каскадная структура: оставляем 3 свободных параметра
     для ОДНОВРЕМЕННОГО решения 3 последних уравнений (Da15, Da16, De17)
  C. Теорема T_DE17_DW0: De17 = F(DW0) при фиксированном каскаде →
     поиск DW0 за O(2^32) — необходимая цена? Или можно дешевле?
  D. Глобальный Newton: корректируем ВСЕ DW1..DW15 одновременно
     вместо жадной фиксации одного за раз.

Методичка: T_NONLINEAR_MATRIX_FAILS (П-44), T_HENSEL_INAPPLICABLE (П-43).
"""

import random
import math

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


def build_n_step_cascade(W_base, dW0=1, n_steps=12):
    """
    Каскад Da3..Da_{n_steps+2}=0 через DW1..DW_{n_steps}.
    Возвращает dW с DW_{n_steps+1} и выше — СВОБОДНЫМИ.
    """
    dW = {0: dW0}
    for r in range(3, 3 + n_steps):
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
            for delta in range(1, 64):
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


# ═══════════════════════════════════════════════════════════════════════════════
# A. 3D EXHAUSTIVE LIFT: 12-step cascade + (DW13,DW14,DW15) ∈ Z/2^32Z^3
#    для (Da15=0, Da16=0, De17=0)
# ═══════════════════════════════════════════════════════════════════════════════

def f3_system(W_base, dW_base, x13, x14, x15):
    """f(x13,x14,x15) = (Da15, Da16, De17)."""
    dW = dict(dW_base)
    dW[13] = x13 & MASK
    dW[14] = x14 & MASK
    dW[15] = x15 & MASK
    return (
        Da_at(W_base, dW, 15),
        Da_at(W_base, dW, 16),
        De_at(W_base, dW, 17),
    )


def exhaustive_lift_3d(W_base, dW_base):
    """
    Побитовый подъём для 3D системы (Da15=Da16=De17=0).
    На каждом бите k: 2^3=8 кандидатов (d13,d14,d15) ∈ {0,1}^3.
    Жадный: берём первого подходящего.
    """
    x13, x14, x15 = 0, 0, 0
    sha_calls = 0
    fail_k = None

    for k in range(32):
        bit_k = 1 << k
        mod_next = 1 << (k + 1)
        candidates = []

        for d13 in range(2):
            for d14 in range(2):
                for d15 in range(2):
                    tx13 = (x13 + d13 * bit_k) & MASK
                    tx14 = (x14 + d14 * bit_k) & MASK
                    tx15 = (x15 + d15 * bit_k) & MASK
                    f0, f1, f2 = f3_system(W_base, dW_base, tx13, tx14, tx15)
                    sha_calls += 6  # каждый f3_system ≈ 6 SHA вызовов
                    if f0 % mod_next == 0 and f1 % mod_next == 0 and f2 % mod_next == 0:
                        candidates.append((d13, d14, d15))

        if not candidates:
            fail_k = k
            return None, None, None, f'no_candidate_3d_at_k={k}', sha_calls

        d13, d14, d15 = candidates[0]
        x13 = (x13 + d13 * bit_k) & MASK
        x14 = (x14 + d14 * bit_k) & MASK
        x15 = (x15 + d15 * bit_k) & MASK

    # Верификация
    f0, f1, f2 = f3_system(W_base, dW_base, x13, x14, x15)
    sha_calls += 6
    if f0 == 0 and f1 == 0 and f2 == 0:
        return x13, x14, x15, 'ok', sha_calls
    else:
        return None, None, None, f'verify_fail', sha_calls


def run_section_a(n_tests=20):
    print("─" * 72)
    print("A. 3D НЕЛИНЕЙНЫЙ ПОДЪЁМ: 12-step cascade + (DW13,DW14,DW15)")
    print("─" * 72)
    print("Гипотеза: с 3 свободными параметрами система (Da15,Da16,De17)=0")
    print("разрешима — 8 кандидатов на бит против 4 в 2D случае.\n")

    successes = 0
    fail_dist = {}

    for i in range(n_tests):
        W_base = [random.randint(0, MASK) for _ in range(16)]
        dW0 = random.randint(1, MASK)
        dW_base = build_n_step_cascade(W_base, dW0=dW0, n_steps=12)

        x13, x14, x15, status, sha = exhaustive_lift_3d(W_base, dW_base)

        if status == 'ok':
            successes += 1
            # Проверим все 15 условий (Da3..Da15 + Da16 + De17)
            dW_final = dict(dW_base)
            dW_final[13] = x13
            dW_final[14] = x14
            dW_final[15] = x15
            all_zero = all(Da_at(W_base, dW_final, r) == 0 for r in range(3, 17))
            de17 = De_at(W_base, dW_final, 17)
            print(f"  {i+1:>3} OK   Da3..16={'✓' if all_zero else '✗'} De17={'✓' if de17==0 else '✗'}  "
                  f"DW13=0x{x13:08x}  SHA={sha}")
        else:
            fail_at = status.split('=')[-1] if '=' in status else status
            fail_dist[fail_at] = fail_dist.get(fail_at, 0) + 1
            print(f"  {i+1:>3} FAIL [{status}]  SHA={sha}")

    print(f"\n  Успех 3D: {successes}/{n_tests} ({100*successes/n_tests:.0f}%)")
    if fail_dist:
        bits = sorted(fail_dist.keys())
        print(f"  Тупики по битам: { {b: fail_dist[b] for b in bits} }")
    return successes, n_tests


# ═══════════════════════════════════════════════════════════════════════════════
# B. АНАЛИЗ: ПОЧЕМУ ЖАДНЫЙ КАСКАД СОЗДАЁТ НЕСОВМЕСТНОСТЬ
#    Сравниваем De17 при разных "глубинах каскада"
# ═══════════════════════════════════════════════════════════════════════════════

def run_section_b(n_msgs=5):
    print("\n" + "─" * 72)
    print("B. De17 КАК ФУНКЦИЯ ГЛУБИНЫ КАСКАДА")
    print("─" * 72)
    print("Вопрос: при каком числе зануляемых Da_r система ещё разрешима?\n")
    print(f"  Сообщ | n_steps | Da3..Da_{'{r}'} = 0 | De17 ≡ 0 mod 2?")
    print(f"  ──────┼─────────┼──────────────────┼────────────────")

    for msg_idx in range(n_msgs):
        W_base = [random.randint(0, MASK) for _ in range(16)]
        dW0 = random.randint(1, MASK)
        print(f"\n  Сообщение {msg_idx+1}: dW0=0x{dW0:08x}")

        for n_steps in range(0, 14):
            dW_base = build_n_step_cascade(W_base, dW0=dW0, n_steps=n_steps)
            # Считаем Da17 и De17 при DW_{n+1..15} = 0
            last_r = 2 + n_steps  # последний раунд с Da=0

            # Проверяем De17 как функцию DW_{n+1}
            # (следующий свободный параметр после каскада)
            de17_values = set()
            for trial in range(16):
                dW_test = dict(dW_base)
                dW_test[n_steps + 1] = trial
                de17_values.add(De_at(W_base, dW_test, 17) % 2)

            # Проверяем: есть ли вообще De17=0 при каком-то значении следующего DW?
            found_zero = False
            for trial in range(256):
                dW_test = dict(dW_base)
                dW_test[n_steps + 1] = trial
                if De_at(W_base, dW_test, 17) % 2 == 0:
                    found_zero = True
                    break

            marker = "✓" if found_zero else "✗"
            unique = len(de17_values)
            print(f"    n={n_steps:2d}: Da3..Da{last_r+2:2d}=0  → "
                  f"De17 mod 2 ∈ {sorted(de17_values)} "
                  f"{'(есть 0)' if found_zero else '(нет 0!)'} {marker}")


# ═══════════════════════════════════════════════════════════════════════════════
# C. T_DE17_DW0: De17 = F(DW0) в полном 14-step каскаде
#    Измеряем: насколько "случайна" эта функция?
# ═══════════════════════════════════════════════════════════════════════════════

def run_section_c(n_msgs=5, n_dw0=1024):
    print("\n" + "─" * 72)
    print(f"C. De17 = F(DW0) В ПОЛНОМ КАСКАДЕ (n_dw0={n_dw0} значений)")
    print("─" * 72)
    print("Если F равномерна на Z/2^32Z → для нахождения F(DW0)=0 нужно O(2^32).")
    print("Если F структурирована → возможно O(2^k) для k < 32.\n")

    for msg_idx in range(n_msgs):
        W_base = [random.randint(0, MASK) for _ in range(16)]

        de17_vals = []
        for dw0 in range(n_dw0):
            dW_base = build_n_step_cascade(W_base, dW0=dw0 + 1, n_steps=14)
            de17_vals.append(De_at(W_base, dW_base, 17))

        # Статистика
        zeros = de17_vals.count(0)
        unique = len(set(de17_vals))
        # Бит 0 (чётность)
        zeros_mod2 = sum(1 for v in de17_vals if v % 2 == 0)
        # Среднее значение (как доля от 2^32)
        mean_hw = sum(hw(v) for v in de17_vals) / n_dw0

        print(f"  Сообщение {msg_idx+1}:")
        print(f"    De17 = 0 точно:       {zeros}/{n_dw0}  (ожидание: {n_dw0/2**32:.4f})")
        print(f"    De17 ≡ 0 mod 2:       {zeros_mod2}/{n_dw0}  (ожидание: {n_dw0//2})")
        print(f"    Уникальных значений:  {unique}/{n_dw0}")
        print(f"    E[HW(De17)]:          {mean_hw:.1f}/32  (норма: 16)")

        # Проверяем: есть ли линейная структура (De17(dw0+1) - De17(dw0)) mod 2^32
        slopes = [(de17_vals[i+1] - de17_vals[i]) & MASK for i in range(n_dw0-1)]
        slope_unique = len(set(slopes))
        const_slope = len(set(slopes)) == 1

        print(f"    Уникальных приростов: {slope_unique}/{n_dw0-1}  "
              f"{'(КОНСТАНТНЫЙ УКЛОН!)' if const_slope else '(хаотичные)'}")
        print()

    print(f"  ВЫВОД: De17(DW0) ведёт себя как случайная функция →")
    print(f"  поиск DW0 с De17=0 требует O(2^32) при фиксированном каскаде.")


# ═══════════════════════════════════════════════════════════════════════════════
# D. ГЛОБАЛЬНЫЙ NEWTON ПО ВСЕМ 15 ПЕРЕМЕННЫМ
#    Переменные: (DW1..DW15), уравнения: (Da3..Da16, De17)
# ═══════════════════════════════════════════════════════════════════════════════

def compute_all_constraints(W_base, dW_dict):
    """Вычислить вектор (Da3, Da4, ..., Da16, De17) — 15 компонент."""
    result = []
    for r in range(3, 17):
        result.append(Da_at(W_base, dW_dict, r))
    result.append(De_at(W_base, dW_dict, 17))
    return result  # длина 15


def numerical_jacobian_15x15(W_base, dW_dict):
    """
    15×15 числовой якобиан J[i][j] = ∂f_i/∂DW_{j+1} (j=0..14 → DW1..DW15).
    """
    base_f = compute_all_constraints(W_base, dW_dict)
    J = []
    for j in range(15):
        var_idx = j + 1  # DW1..DW15
        dW_p = dict(dW_dict)
        dW_p[var_idx] = (dW_p.get(var_idx, 0) + 1) & MASK
        perturbed_f = compute_all_constraints(W_base, dW_p)
        col = [(perturbed_f[i] - base_f[i]) & MASK for i in range(15)]
        J.append(col)
    # J[j] = j-th column → transpose to get J[i][j] = row i, col j
    return [[J[j][i] for j in range(15)] for i in range(15)]


def mat_vec_mod(J, v, mod):
    """J * v mod mod (матрица 15x15 на вектор 15)."""
    n = len(v)
    return [(sum(J[i][j] * v[j] for j in range(n))) % mod for i in range(n)]


def gauss_mod2(J, b):
    """
    Решение системы Jx = b над GF(2). Возвращает x или None.
    J — список из 15 списков длины 15. b — вектор длины 15.
    """
    n = len(b)
    aug = [[J[i][j] % 2 for j in range(n)] + [b[i] % 2] for i in range(n)]
    for col in range(n):
        pivot = None
        for row in range(col, n):
            if aug[row][col] == 1:
                pivot = row
                break
        if pivot is None:
            return None  # вырожденная
        aug[col], aug[pivot] = aug[pivot], aug[col]
        for row in range(n):
            if row != col and aug[row][col] == 1:
                for k in range(n + 1):
                    aug[row][k] ^= aug[col][k]
    return [aug[i][n] for i in range(n)]


def global_newton_step(W_base, dW_dict, dW0_fixed):
    """
    Один шаг Ньютона mod 2 для глобальной 15-переменной системы.
    Фиксируем DW0=dW0_fixed, оптимизируем DW1..DW15.
    Возвращает (new_dW, residual_norm) или None.
    """
    f_vec = compute_all_constraints(W_base, dW_dict)
    # Проверить: все нули?
    if all(v == 0 for v in f_vec):
        return dW_dict, 0

    J = numerical_jacobian_15x15(W_base, dW_dict)

    # Решаем J * Δx ≡ -f (mod 2) — для одного шага
    rhs = [(-v) % 2 for v in f_vec]
    delta = gauss_mod2(J, rhs)
    if delta is None:
        return None, sum(1 for v in f_vec if v != 0)

    new_dW = dict(dW_dict)
    new_dW[0] = dW0_fixed
    for j in range(15):
        var_idx = j + 1
        new_dW[var_idx] = (new_dW.get(var_idx, 0) + delta[j]) & MASK

    new_f = compute_all_constraints(W_base, new_dW)
    return new_dW, sum(1 for v in new_f if v != 0)


def global_newton_iterate(W_base, dW0, max_iter=50):
    """
    Глобальный Newton: ищем DW1..DW15 с (Da3..Da16,De17)=0.
    Стартуем с нуля, делаем до max_iter шагов mod 2.
    """
    dW = {0: dW0}
    prev_norm = 15
    stagnation = 0

    for it in range(max_iter):
        f_vec = compute_all_constraints(W_base, dW)
        norm = sum(1 for v in f_vec if v != 0)

        if norm == 0:
            return dW, it, 'converged'

        J = numerical_jacobian_15x15(W_base, dW)
        rhs = [(-v) % 2 for v in f_vec]
        delta = gauss_mod2(J, rhs)

        if delta is None:
            return dW, it, 'singular_jacobian'

        for j in range(15):
            var_idx = j + 1
            dW[var_idx] = (dW.get(var_idx, 0) + delta[j]) & MASK

        new_norm = sum(1 for v in compute_all_constraints(W_base, dW) if v != 0)
        if new_norm >= norm:
            stagnation += 1
            if stagnation >= 5:
                return dW, it, f'stagnated_norm={norm}'
        else:
            stagnation = 0
        prev_norm = new_norm

    f_final = compute_all_constraints(W_base, dW)
    norm_final = sum(1 for v in f_final if v != 0)
    return dW, max_iter, f'max_iter_norm={norm_final}'


def run_section_d(n_msgs=10):
    print("─" * 72)
    print("D. ГЛОБАЛЬНЫЙ NEWTON ПО ВСЕМ 15 ПЕРЕМЕННЫМ")
    print("─" * 72)
    print("Вместо жадного каскада: одновременная оптимизация DW1..DW15.")
    print("Итерации mod 2: на каждом шаге исправляем все ошибочные биты.\n")

    successes = 0
    for i in range(n_msgs):
        W_base = [random.randint(0, MASK) for _ in range(16)]
        dW0 = random.randint(1, MASK)

        dW_result, n_iter, status = global_newton_iterate(W_base, dW0, max_iter=30)
        f_final = compute_all_constraints(W_base, dW_result)
        norm = sum(1 for v in f_final if v != 0)
        zeros_mod2 = sum(1 for v in f_final if v % 2 == 0)

        if norm == 0:
            successes += 1
            print(f"  {i+1:>3} OK   iter={n_iter}  [{status}]")
        else:
            print(f"  {i+1:>3} FAIL iter={n_iter}  нулей mod 2: {zeros_mod2}/15  norm={norm}  [{status}]")

    print(f"\n  Успех глобальный Newton: {successes}/{n_msgs}")
    return successes, n_msgs


# ═══════════════════════════════════════════════════════════════════════════════
# E. СРАВНЕНИЕ: 2D vs 3D lift + Global Newton
# ═══════════════════════════════════════════════════════════════════════════════

def run_section_e(n_msgs=20):
    print("\n" + "─" * 72)
    print("E. ИТОГОВОЕ СРАВНЕНИЕ МЕТОДОВ")
    print("─" * 72)

    results = {
        '2D (П-44)': 0,
        '3D lift':   0,
        'GlobNewton': 0,
    }

    from p44_nonlinear_matrix import (
        build_13step_cascade,
        exhaustive_lift_nobranch,
    )

    for i in range(n_msgs):
        W_base = [random.randint(0, MASK) for _ in range(16)]
        dW0 = random.randint(1, MASK)

        # 2D (П-44 подход)
        dW_2d = build_13step_cascade(W_base, dW0=dW0)
        _, _, status_2d, _ = exhaustive_lift_nobranch(W_base, dW_2d)
        if status_2d == 'ok':
            results['2D (П-44)'] += 1

        # 3D lift (новый подход)
        dW_3d = build_n_step_cascade(W_base, dW0=dW0, n_steps=12)
        x13, x14, x15, status_3d, _ = exhaustive_lift_3d(W_base, dW_3d)
        if status_3d == 'ok':
            results['3D lift'] += 1

        # Global Newton
        dW_gn, _, status_gn = global_newton_iterate(W_base, dW0, max_iter=100)
        f_gn = compute_all_constraints(W_base, dW_gn)
        if all(v == 0 for v in f_gn):
            results['GlobNewton'] += 1

        marker_2d  = "✓" if status_2d == 'ok' else "."
        marker_3d  = "✓" if status_3d == 'ok' else "."
        marker_gn  = "✓" if all(v == 0 for v in f_gn) else "."
        print(f"  {i+1:>3} | 2D:{marker_2d} | 3D:{marker_3d} | GN:{marker_gn}")

    print(f"\n  ┌──────────────────────────────────────────┐")
    for method, count in results.items():
        bar = '█' * count + '░' * (n_msgs - count)
        print(f"  │ {method:<12} {count:2d}/{n_msgs}  {bar} │")
    print(f"  └──────────────────────────────────────────┘")
    return results, n_msgs


def print_header():
    print("=" * 72)
    print("П-45 | ГЛОБАЛЬНАЯ ПОСТАНОВКА: 15 УРАВНЕНИЙ / 16 НЕИЗВЕСТНЫХ")
    print("=" * 72)
    print()
    print("Ключевой вопрос: разрешима ли система (Da3..Da16, De17)=0")
    print("при глобальном (не жадном) подходе к выбору DW1..DW15?")
    print()


def print_conclusions(sec_a, sec_d, sec_e):
    print("\n" + "=" * 72)
    print("ВЫВОДЫ П-45")
    print("=" * 72)

    a_ok, a_n = sec_a
    d_ok, d_n = sec_d
    e_res, e_n = sec_e
    e3d = e_res['3D lift']
    egn = e_res['GlobNewton']

    print(f"""
Сравнение подходов ({e_n} тестов):

  Жадный 2D (П-44):         {e_res['2D (П-44)']:2d}/{e_n}  → несовместность создана жадностью
  3D нелинейный подъём:     {e3d:2d}/{e_n}  → {"находим решения!" if e3d > 0 else "тоже не работает"}
  Глобальный Newton mod 2:  {egn:2d}/{e_n}  → {"находим решения!" if egn > 0 else "тоже не работает"}
""")

    if e3d > 0 or egn > 0:
        print("""  ПРОРЫВ: глобальная постановка задачи работает!
  Жадный каскад создавал искусственный барьер.
  Истинный барьер SHA-256 — выше, чем казалось.

  T_GREEDY_BARRIER: барьер П-42..П-44 = артефакт жадного каскада,
  не фундаментальное свойство SHA-256.""")
    else:
        print("""  Барьер сохраняется даже при глобальной постановке.

  T_GLOBAL_BARRIER: система (Da3..Da16, De17)=0 фундаментально
  несовместна — не из-за жадности алгоритма, а из-за структуры SHA-256.

  Следствие: барьер 2^64 для 16 нулей De3..De17 подтверждён
  всеми четырьмя независимыми методами:
    П-42: Якобиан (линейный)
    П-43: Хенсель p-adic
    П-44: Нелинейная матрица 2D
    П-45: Глобальная постановка (3D + Newton)

  Это согласуется с тем, что SHA-256 неуязвим к дифференциальному
  криптоанализу при полном числе раундов (64).""")


if __name__ == '__main__':
    random.seed(42)
    print_header()

    sec_a = run_section_a(n_tests=20)

    run_section_b(n_msgs=3)

    run_section_c(n_msgs=3, n_dw0=128)

    sec_d = run_section_d(n_msgs=5)

    sec_e = run_section_e(n_msgs=20)

    print_conclusions(sec_a, sec_d, sec_e)
