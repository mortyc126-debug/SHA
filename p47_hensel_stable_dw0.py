"""
П-47: HENSEL-STABLE DW0 — ПРАВИЛЬНЫЙ ПОБИТОВЫЙ ПОДЪЁМ.

Диагноз П-46:
  rank=15 при старте (DW1..DW15=0) — необходимое, но НЕ достаточное условие.
  Причина: Hensel lifting требует СТАРТОВОГО РЕШЕНИЯ MOD 2, а не mod 2^32.
  Мы начинали с x=0, которое НЕ является решением mod 2 (f(0) ≢ 0 mod 2).
  Поэтому "lifting" не работал — не было базы для подъёма.

Правильный алгоритм (Hensel lifting):
  Шаг 0: найти x₀ ∈ {0,1}^15 такое, что f(x₀) ≡ 0 mod 2
  Шаг k: из x_k ∈ (Z/2^k)^15 получить x_{k+1} ∈ (Z/2^{k+1})^15:
          J(x_k) * delta ≡ -(f(x_k) / 2^k) mod 2
          x_{k+1} = x_k + delta * 2^k

Ключевые вопросы П-47:
  A. Для фиксированного DW0: существует ли x₀ ∈ {0,1}^15 с f(x₀) ≡ 0 mod 2?
     Стратегия поиска: линейная алгебра над GF(2) + проверка нелинейных поправок.
  B. Если x₀ найден: поднимается ли он до решения mod 2^32?
     Трассировка Hensel по всем 32 битам.
  C. Зависимость от DW0: для каких DW0 алгоритм работает?
  D. Итог: существует ли вообще решение системы (Da3..Da16, De17)=0?
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


def compute_all_constraints(W_base, dW_dict):
    """Вектор (Da3..Da16, De17) — 15 компонент."""
    result = [Da_at(W_base, dW_dict, r) for r in range(3, 17)]
    result.append(De_at(W_base, dW_dict, 17))
    return result


def numerical_jacobian_15x15(W_base, dW_dict):
    base_f = compute_all_constraints(W_base, dW_dict)
    cols = []
    for j in range(15):
        var_idx = j + 1
        dW_p = dict(dW_dict)
        dW_p[var_idx] = (dW_p.get(var_idx, 0) + 1) & MASK
        pf = compute_all_constraints(W_base, dW_p)
        cols.append([(pf[i] - base_f[i]) & MASK for i in range(15)])
    return [[cols[j][i] for j in range(15)] for i in range(15)]


def rank_mod2(J):
    n = len(J)
    A = [[J[i][j] % 2 for j in range(n)] for i in range(n)]
    rank = 0
    pivot_col = 0
    for row in range(n):
        found = (-1, -1)
        for c in range(pivot_col, n):
            for r in range(row, n):
                if A[r][c] == 1:
                    found = (r, c)
                    break
            if found[0] != -1:
                break
        if found[0] == -1:
            break
        r, c = found
        A[row], A[r] = A[r], A[row]
        pivot_col = c + 1
        for r2 in range(n):
            if r2 != row and A[r2][c] == 1:
                for k in range(n):
                    A[r2][k] ^= A[row][k]
        rank += 1
    return rank


def gauss_mod2(J, b):
    n = len(b)
    aug = [[J[i][j] % 2 for j in range(n)] + [b[i] % 2] for i in range(n)]
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


# ═══════════════════════════════════════════════════════════════════════════════
# A. ПОИСК СТАРТОВОГО СЕМЕНИ mod 2
#    Для фиксированного DW0 ищем x₀ ∈ {0,1}^15 с f(x₀) ≡ 0 mod 2.
#
#    Стратегия: GF(2) линейная аппроксимация как стартовая точка,
#    затем локальный перебор (flip) до нахождения настоящего нулевого семени.
# ═══════════════════════════════════════════════════════════════════════════════

def find_seed_mod2_linear(W_base, dw0):
    """
    Стратегия 1: линейная аппроксимация J*x ≡ -f(0) mod 2 → кандидат x₀.
    Затем проверяем и корректируем нелинейными поправками.
    """
    dW = {0: dw0}
    for j in range(1, 16):
        dW[j] = 0

    f0 = compute_all_constraints(W_base, dW)
    J = numerical_jacobian_15x15(W_base, dW)
    rhs = [(-v) % 2 for v in f0]
    x0 = gauss_mod2(J, rhs)
    if x0 is None:
        return None  # система несовместна линейно

    # Проверяем нелинейное решение
    dW_cand = {0: dw0}
    for j in range(15):
        dW_cand[j + 1] = x0[j]
    f_cand = compute_all_constraints(W_base, dW_cand)
    if all(v % 2 == 0 for v in f_cand):
        return dW_cand

    return None  # линейный кандидат не удовлетворяет нелинейной системе


def find_seed_mod2_random(W_base, dw0, n_trials=2000, rng=None):
    """
    Стратегия 2: случайный перебор x₀ ∈ {0,1}^15.
    Эффективнее с ранним выходом по первым несовпавшим компонентам.
    """
    if rng is None:
        rng = random
    for _ in range(n_trials):
        bits = rng.randrange(1 << 15)
        dW_cand = {0: dw0}
        for j in range(15):
            dW_cand[j + 1] = (bits >> j) & 1

        # Проверяем компоненты с ранним выходом
        ok = True
        for r in range(3, 17):
            if Da_at(W_base, dW_cand, r) % 2 != 0:
                ok = False
                break
        if ok and De_at(W_base, dW_cand, 17) % 2 != 0:
            ok = False
        if ok:
            return dW_cand
    return None


def count_seeds_mod2_exhaustive(W_base, dw0):
    """
    Полный перебор {0,1}^15: считаем, сколько x₀ дают f(x₀) ≡ 0 mod 2.
    Возвращает (count, first_seed_or_None).
    """
    count = 0
    first = None
    for bits in range(1 << 15):
        dW_cand = {0: dw0}
        for j in range(15):
            dW_cand[j + 1] = (bits >> j) & 1

        ok = True
        for r in range(3, 17):
            if Da_at(W_base, dW_cand, r) % 2 != 0:
                ok = False
                break
        if ok and De_at(W_base, dW_cand, 17) % 2 != 0:
            ok = False
        if ok:
            count += 1
            if first is None:
                first = dict(dW_cand)
    return count, first


# ═══════════════════════════════════════════════════════════════════════════════
# B. HENSEL LIFTING ОТ СЕМЕНИ mod 2 → mod 2^32
# ═══════════════════════════════════════════════════════════════════════════════

def hensel_lift_from_seed(W_base, seed_dW, verbose=False):
    """
    Правильный Hensel lifting:
      x_0 = seed_dW (решение mod 2)
      На каждом бите k:
        J(x_k) * delta ≡ -(f(x_k) / 2^k) mod 2
        x_{k+1} = x_k + delta * 2^k
    Возвращает (dW_result, status, rank_trace).
    """
    dW = dict(seed_dW)
    rank_trace = []

    # Проверяем стартовое условие: f ≡ 0 mod 2
    f_start = compute_all_constraints(W_base, dW)
    if not all(v % 2 == 0 for v in f_start):
        return dW, 'bad_seed', []

    for bit in range(32):
        # Проверяем текущее решение mod 2^(bit+1)
        f_vec = compute_all_constraints(W_base, dW)
        mod_k = 1 << bit
        mod_kp1 = 1 << (bit + 1)

        # Все компоненты должны быть ≡ 0 mod 2^bit
        if not all(v % mod_k == 0 for v in f_vec):
            return dW, f'lift_fail_mod_check_bit={bit}', rank_trace

        # Готово?
        if all(v == 0 for v in f_vec):
            return dW, 'converged_early', rank_trace

        # Вычисляем сокращённый RHS: r_i = f_i / 2^bit mod 2
        rhs = [(v // mod_k) % 2 for v in f_vec]

        # Якобиан и ранг
        J = numerical_jacobian_15x15(W_base, dW)
        r = rank_mod2(J)
        rank_trace.append(r)

        if verbose:
            nonzero_rhs = sum(1 for v in rhs if v != 0)
            print(f"    bit={bit}: rank(J)={r}  nonzero_rhs={nonzero_rhs}")

        if r < 15:
            # Проверяем: совместна ли система?
            delta = gauss_mod2(J, rhs)
            if delta is None:
                # Может быть несовместна из-за сингулярности
                # Попробуем backtrack: delta=0 (не поднимаем этот бит)
                delta = [0] * 15
                dW_test = dict(dW)
                for j in range(15):
                    dW_test[j + 1] = (dW_test.get(j + 1, 0) + 0 * (1 << bit)) & MASK
                f_test = compute_all_constraints(W_base, dW_test)
                if all(v % mod_kp1 == 0 for v in f_test):
                    if verbose:
                        print(f"    bit={bit}: rank<15, delta=0 works!")
                    continue  # delta=0 работает
                return dW, f'singular_inconsistent_bit={bit}', rank_trace
            # delta найдена — используем её
        else:
            delta = gauss_mod2(J, rhs)
            if delta is None:
                return dW, f'gauss_fail_bit={bit}', rank_trace

        # Обновляем x_{k+1} = x_k + delta * 2^k
        for j in range(15):
            dW[j + 1] = (dW.get(j + 1, 0) + delta[j] * (1 << bit)) & MASK

        # Проверяем: f(x_{k+1}) ≡ 0 mod 2^{k+1}?
        f_new = compute_all_constraints(W_base, dW)
        if not all(v % mod_kp1 == 0 for v in f_new):
            # Это неожиданно — нелинейная поправка испортила результат
            return dW, f'nonlinear_correction_fail_bit={bit}', rank_trace

    f_final = compute_all_constraints(W_base, dW)
    if all(v == 0 for v in f_final):
        return dW, 'converged', rank_trace
    norm = sum(1 for v in f_final if v != 0)
    return dW, f'done_norm={norm}', rank_trace


# ═══════════════════════════════════════════════════════════════════════════════
# A. СЕКЦИЯ A: статистика семян mod 2
# ═══════════════════════════════════════════════════════════════════════════════

def run_section_a(n_msgs=5, n_dw0_per_msg=50):
    print("─" * 72)
    print("A. ПОИСК СЕМЯН mod 2: x₀ ∈ {0,1}^15 с f(x₀) ≡ 0 mod 2")
    print("─" * 72)
    print("Метод 1: линейная аппроксимация J*x ≡ -f(0) mod 2")
    print("Метод 2: случайный перебор (2000 попыток)\n")

    linear_ok = 0
    random_ok = 0
    total = 0

    for msg_idx in range(n_msgs):
        W_base = [random.randint(0, MASK) for _ in range(16)]
        lin_this = 0
        rnd_this = 0

        for _ in range(n_dw0_per_msg):
            dw0 = random.randint(1, MASK)

            # Метод 1: линейный
            seed = find_seed_mod2_linear(W_base, dw0)
            if seed is not None:
                lin_this += 1

            # Метод 2: случайный
            seed2 = find_seed_mod2_random(W_base, dw0, n_trials=2000)
            if seed2 is not None:
                rnd_this += 1

            total += 1

        linear_ok += lin_this
        random_ok += rnd_this
        print(f"  Сообщение {msg_idx+1}: "
              f"линейный={lin_this}/{n_dw0_per_msg}  "
              f"случайный={rnd_this}/{n_dw0_per_msg}")

    print(f"\n  Итого: линейный={linear_ok}/{total}={100*linear_ok/total:.1f}%  "
          f"случайный={random_ok}/{total}={100*random_ok/total:.1f}%")
    print(f"\n  Ключевой вопрос: если линейный < случайный → нелинейная система")
    print(f"  имеет решения, но их не находит линейная аппроксимация.")
    return linear_ok, random_ok, total


# ═══════════════════════════════════════════════════════════════════════════════
# B. EXHAUSTIVE count: сколько семян mod 2 существует для фикс. DW0?
# ═══════════════════════════════════════════════════════════════════════════════

def run_section_b(n_msgs=3, n_dw0_per_msg=5):
    print("\n" + "─" * 72)
    print("B. ПОЛНЫЙ ПЕРЕБОР {0,1}^15: КОЛИЧЕСТВО СЕМЯН mod 2")
    print("─" * 72)
    print("Сколько x₀ из 2^15=32768 дают f(x₀) ≡ 0 mod 2?\n")

    for msg_idx in range(n_msgs):
        W_base = [random.randint(0, MASK) for _ in range(16)]
        print(f"  Сообщение {msg_idx+1}:")

        for _ in range(n_dw0_per_msg):
            dw0 = random.randint(1, MASK)
            count, first = count_seeds_mod2_exhaustive(W_base, dw0)
            log2_count = math.log2(count) if count > 0 else -float('inf')
            print(f"    DW0=0x{dw0:08x}:  {count} семян  "
                  f"(≈2^{log2_count:.1f})  "
                  f"первое={list(first.values())[1:6] if first else 'None'}")
        print()


# ═══════════════════════════════════════════════════════════════════════════════
# C. HENSEL LIFTING от найденных семян
# ═══════════════════════════════════════════════════════════════════════════════

def run_section_c(n_msgs=5, n_dw0_per_msg=10):
    print("─" * 72)
    print("C. HENSEL LIFTING: seed mod 2 → решение mod 2^32")
    print("─" * 72)
    print("Для каждого DW0 находим семя (случайный перебор),")
    print("затем поднимаем по Hensel до мод 2^32.\n")

    full_lifts = 0
    seeds_found = 0
    total = 0
    fail_dist = {}

    for msg_idx in range(n_msgs):
        W_base = [random.randint(0, MASK) for _ in range(16)]
        print(f"  Сообщение {msg_idx+1}:")

        for _ in range(n_dw0_per_msg):
            dw0 = random.randint(1, MASK)
            total += 1

            # Ищем семя
            seed = find_seed_mod2_random(W_base, dw0, n_trials=5000)
            if seed is None:
                seed = find_seed_mod2_linear(W_base, dw0)

            if seed is None:
                print(f"    DW0=0x{dw0:08x}: семя НЕ найдено")
                fail_dist['no_seed'] = fail_dist.get('no_seed', 0) + 1
                continue

            seeds_found += 1

            # Hensel lifting
            dW_result, status, rank_trace = hensel_lift_from_seed(W_base, seed)
            f_final = compute_all_constraints(W_base, dW_result)
            success = all(v == 0 for v in f_final)

            if success:
                full_lifts += 1
                all_da = all(Da_at(W_base, dW_result, r) == 0 for r in range(3, 17))
                de17 = De_at(W_base, dW_result, 17)
                print(f"    DW0=0x{dw0:08x}: РЕШЕНИЕ НАЙДЕНО! "
                      f"Da3..16={'✓' if all_da else '✗'}  De17={'✓' if de17==0 else '✗'}")
            else:
                key = status.split('_')[0] + '_' + status.split('_')[1] if '_' in status else status
                fail_dist[key] = fail_dist.get(key, 0) + 1
                min_r = min(rank_trace) if rank_trace else '?'
                print(f"    DW0=0x{dw0:08x}: подъём={status}  min_rank={min_r}")

        print()

    print(f"  Семян найдено: {seeds_found}/{total}")
    print(f"  Полных подъёмов (решений): {full_lifts}/{seeds_found}")
    if fail_dist:
        print(f"  Статус отказов: {fail_dist}")
    return full_lifts, seeds_found, total


# ═══════════════════════════════════════════════════════════════════════════════
# D. ДЕТАЛЬНАЯ ТРАССИРОВКА для 3 лучших кандидатов
# ═══════════════════════════════════════════════════════════════════════════════

def run_section_d(n_msgs=3):
    print("─" * 72)
    print("D. ДЕТАЛЬНАЯ ТРАССИРОВКА HENSEL LIFTING")
    print("─" * 72)
    print("Для каждого сообщения: находим лучший DW0 и трассируем по битам.\n")

    for msg_idx in range(n_msgs):
        W_base = [random.randint(0, MASK) for _ in range(16)]
        print(f"  Сообщение {msg_idx+1}:")

        best_bits = -1
        best_dw0 = None
        best_seed = None

        # Ищем DW0 с наилучшим подъёмом
        for _ in range(30):
            dw0 = random.randint(1, MASK)
            seed = find_seed_mod2_random(W_base, dw0, n_trials=3000)
            if seed is None:
                continue

            # Трассируем
            dW_result, status, rank_trace = hensel_lift_from_seed(W_base, seed)

            # Считаем сколько бит успешно поднялось
            if 'converged' in status:
                bits_ok = 32
            elif 'bit=' in status:
                bits_ok = int(status.split('bit=')[1])
            else:
                bits_ok = len(rank_trace)

            if bits_ok > best_bits:
                best_bits = bits_ok
                best_dw0 = dw0
                best_seed = seed
                best_status = status
                best_rank_trace = rank_trace

        if best_seed is None:
            print(f"    Семян не найдено\n")
            continue

        print(f"    Лучший DW0=0x{best_dw0:08x}, биты поднято: {best_bits}/32")
        print(f"    Статус: {best_status}")
        if best_rank_trace:
            ranks_str = ''.join(str(r) for r in best_rank_trace[:16])
            print(f"    Ранги J (биты 0..{min(15,len(best_rank_trace)-1)}): {ranks_str}")

        # Полная трассировка
        print(f"    Детальная трассировка:")
        dW_result, status, _ = hensel_lift_from_seed(W_base, best_seed, verbose=True)
        print()


# ═══════════════════════════════════════════════════════════════════════════════
# E. ИТОГОВЫЙ ДИАГНОЗ
# ═══════════════════════════════════════════════════════════════════════════════

def run_section_e(full_lifts, seeds_found, total):
    print("─" * 72)
    print("E. ИТОГОВЫЙ ДИАГНОЗ")
    print("─" * 72)
    print()

    if full_lifts > 0:
        print(f"  ПРОРЫВ! Найдено {full_lifts} полных решений системы!")
        print(f"  → T_GLOBAL_BARRIER ОПРОВЕРГНУТ.")
        print(f"  → Система (Da3..Da16, De17)=0 РАЗРЕШИМА.")
        print(f"  → Правильный Hensel lifting работает.")
        print(f"\n  Сложность: найти семя (random ~{total//max(seeds_found,1)} попыток) + 32 бита подъёма.")
    elif seeds_found > 0:
        print(f"  Семян mod 2 найдено: {seeds_found}/{total}.")
        print(f"  Но подъём до mod 2^32 не удался ни разу.")
        print(f"\n  Диагноз: система разрешима mod 2 (семена существуют),")
        print(f"  но разрешима ли mod 2^32 — неизвестно.")
        print(f"\n  → П-48: детальный анализ барьера при подъёме (на каком бите?),")
        print(f"    возможно с backtracking по ветвям rank=14.")
    else:
        print(f"  Семян mod 2 не найдено ни для одного DW0.")
        print(f"\n  Сильный диагноз: система f(x) ≡ 0 mod 2 НЕСОВМЕСТНА")
        print(f"  для всех протестированных DW0.")
        print(f"\n  → T_GLOBAL_BARRIER ПОДТВЕРЖДЁН на уровне mod 2.")
        print(f"  → Нет решений даже в GF(2)^15.")


def print_header():
    print("=" * 72)
    print("П-47 | HENSEL-STABLE DW0: ПРАВИЛЬНЫЙ ПОБИТОВЫЙ ПОДЪЁМ")
    print("=" * 72)
    print()
    print("Ключевая коррекция: Hensel lifting требует стартового решения mod 2.")
    print("П-46 не имел его → якобиан был бессмысленен.")
    print()


if __name__ == '__main__':
    random.seed(42)
    print_header()

    lin_ok, rnd_ok, total_a = run_section_a(n_msgs=5, n_dw0_per_msg=50)

    run_section_b(n_msgs=3, n_dw0_per_msg=5)

    full_lifts, seeds_found, total_c = run_section_c(n_msgs=5, n_dw0_per_msg=10)

    run_section_d(n_msgs=3)

    run_section_e(full_lifts, seeds_found, total_c)
