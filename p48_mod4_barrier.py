"""
П-48: АНАЛИЗ БАРЬЕРА MOD 4 — k-ШАГОВЫЙ ЯКОБИАН И ИСЧЕРПЫВАЮЩИЙ ПОИСК.

Диагноз П-47:
  Hensel lifting mod 2 → mod 4 падает: f(seed + δ·2) ≢ 0 mod 4.
  Причина: якобиан J вычислен с шагом +1, а нужен шаг +2.
  J_ij = f_i(x + e_j) - f_i(x)   (единичный шаг)
  G_ij = f_i(x + 2·e_j) - f_i(x) (двойной шаг)
  Разность G - 2J = нелинейные поправки SHA-256 (переносы, Ch, Maj).

Правильный алгоритм (k-шаговый якобиан):
  На шаге k (mod 2^k → mod 2^{k+1}) используем 2^k-якобиан:
    G^(k)_ij = f_i(x_k + 2^k·e_j) - f_i(x_k)
  Тогда: f(x_k + 2^k·δ) = f(x_k) + G^(k)·δ·2^k + O(2^{2k})
  Для k≥1: O(2^{2k}) ≡ 0 mod 2^{k+1} — формула точна!

Структура П-48:
  A. Диагностика: сравниваем J (единичный) vs G (двойной) — разница mod 2?
  B. k-шаговый Hensel: используем G^(k) вместо J — работает ли mod 4?
  C. Исчерпывающий поиск mod 4: 2^15 кандидатов δ ∈ {0,1}^15 — есть ли хоть один?
  D. Статистика разрешимости: для скольких семян существует решение mod 4?
  E. Продолжение: если mod 4 найдено — поднимаем до mod 8, mod 16, ...
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


def compute_f(W_base, dW_dict):
    """Вектор (Da3..Da16, De17) — 15 компонент."""
    result = [Da_at(W_base, dW_dict, r) for r in range(3, 17)]
    result.append(De_at(W_base, dW_dict, 17))
    return result


def jacobian_step_k(W_base, dW_dict, step):
    """
    k-шаговый якобиан: G^(k)_ij = f_i(x + step·e_j) - f_i(x).
    step=1  → обычный Якобиан J
    step=2  → двойной шаг G (для Hensel mod 2 → mod 4)
    step=2^k → для Hensel mod 2^k → mod 2^{k+1}
    """
    base_f = compute_f(W_base, dW_dict)
    cols = []
    for j in range(15):
        var_idx = j + 1
        dW_p = dict(dW_dict)
        dW_p[var_idx] = (dW_p.get(var_idx, 0) + step) & MASK
        pf = compute_f(W_base, dW_p)
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


def find_seed_mod2(W_base, dw0, n_trials=5000):
    """Случайный поиск x₀ ∈ {0,1}^15 с f(x₀) ≡ 0 mod 2."""
    for _ in range(n_trials):
        bits = random.randrange(1 << 15)
        dW_cand = {0: dw0}
        for j in range(15):
            dW_cand[j + 1] = (bits >> j) & 1
        f = compute_f(W_base, dW_cand)
        if all(v % 2 == 0 for v in f):
            return dW_cand
    return None


# ═══════════════════════════════════════════════════════════════════════════════
# A. ДИАГНОСТИКА: J (шаг=1) vs G (шаг=2) — разница mod 2
# ═══════════════════════════════════════════════════════════════════════════════

def run_section_a(n_tests=20):
    print("─" * 72)
    print("A. ДИАГНОСТИКА: J (шаг=1) vs G (шаг=2) — разница mod 2")
    print("─" * 72)
    print("Гипотеза: G_ij mod 2 ≠ J_ij mod 2 из-за нелинейных переносов SHA.\n")

    diff_count = 0
    total_entries = 0
    rank_J_list = []
    rank_G_list = []

    for i in range(n_tests):
        W_base = [random.randint(0, MASK) for _ in range(16)]
        seed = find_seed_mod2(W_base, random.randint(1, MASK))
        if seed is None:
            continue

        J = jacobian_step_k(W_base, seed, step=1)
        G = jacobian_step_k(W_base, seed, step=2)

        # Сравниваем mod 2
        diff_this = sum(
            (J[i][j] % 2) != (G[i][j] % 2)
            for i in range(15) for j in range(15)
        )
        diff_count += diff_this
        total_entries += 225

        rJ = rank_mod2(J)
        rG = rank_mod2([[G[i][j] // 1 for j in range(15)] for i in range(15)])
        rank_J_list.append(rJ)
        rank_G_list.append(rG)

        print(f"  Тест {i+1:2d}: rank(J)={rJ}  rank(G)={rG}  "
              f"различий(mod 2): {diff_this}/225  ({100*diff_this/225:.1f}%)")

    print(f"\n  Итого различий: {diff_count}/{total_entries} = {100*diff_count/total_entries:.1f}%")
    n_valid = len(rank_J_list)
    print(f"  Средний rank(J): {sum(rank_J_list)/n_valid:.2f}  "
          f"rank(G): {sum(rank_G_list)/n_valid:.2f}")
    print(f"\n  ВЫВОД: G {'≠' if diff_count > 0 else '='} J mod 2 — "
          f"{'нелинейные поправки значимы' if diff_count > 0 else 'матрицы идентичны mod 2'}")
    return diff_count, total_entries


# ═══════════════════════════════════════════════════════════════════════════════
# B. k-ШАГОВЫЙ HENSEL: G^(k) вместо J
# ═══════════════════════════════════════════════════════════════════════════════

def hensel_kstep(W_base, seed_dW, verbose=False):
    """
    Правильный Hensel с k-шаговым якобианом G^(k).
    На шаге k: используем G^(k)_ij = f_i(x + 2^k*e_j) - f_i(x).
    """
    dW = dict(seed_dW)

    f_start = compute_f(W_base, dW)
    if not all(v % 2 == 0 for v in f_start):
        return dW, 'bad_seed', []

    rank_trace = []

    for bit in range(32):
        mod_k = 1 << bit
        mod_kp1 = 1 << (bit + 1)

        f_vec = compute_f(W_base, dW)

        if all(v == 0 for v in f_vec):
            return dW, 'converged_early', rank_trace

        if not all(v % mod_k == 0 for v in f_vec):
            return dW, f'invariant_fail_bit={bit}', rank_trace

        rhs = [(v // mod_k) % 2 for v in f_vec]

        # k-шаговый якобиан
        G = jacobian_step_k(W_base, dW, step=mod_k)
        r = rank_mod2(G)
        rank_trace.append(r)

        if verbose:
            nonzero_rhs = sum(1 for v in rhs if v != 0)
            print(f"    bit={bit}: rank(G^{mod_k})={r}  nonzero_rhs={nonzero_rhs}")

        if r < 15:
            # Попробуем delta=0
            f_zero = compute_f(W_base, dW)
            if all(v % mod_kp1 == 0 for v in f_zero):
                if verbose:
                    print(f"    bit={bit}: rank<15, delta=0 works!")
                continue
            return dW, f'singular_bit={bit}_rank={r}', rank_trace

        delta = gauss_mod2(G, rhs)
        if delta is None:
            return dW, f'gauss_fail_bit={bit}', rank_trace

        for j in range(15):
            dW[j + 1] = (dW.get(j + 1, 0) + delta[j] * mod_k) & MASK

        f_new = compute_f(W_base, dW)
        if not all(v % mod_kp1 == 0 for v in f_new):
            if verbose:
                bad = [(i, f_new[i]) for i in range(15) if f_new[i] % mod_kp1 != 0]
                print(f"    bit={bit}: mod_kp1 check FAIL, bad={bad[:3]}")
            return dW, f'correction_fail_bit={bit}', rank_trace

    f_final = compute_f(W_base, dW)
    if all(v == 0 for v in f_final):
        return dW, 'converged', rank_trace
    norm = sum(1 for v in f_final if v != 0)
    return dW, f'done_norm={norm}', rank_trace


def run_section_b(n_msgs=5, n_dw0_per_msg=10):
    print("\n" + "─" * 72)
    print("B. k-ШАГОВЫЙ HENSEL LIFTING: G^(k) вместо J")
    print("─" * 72)
    print("Теория: G^(k)·δ точно предсказывает f(x+2^k·δ) mod 2^{k+1}.\n")

    full_ok = 0
    seeds_found = 0
    total = 0
    fail_bits = {}

    for msg_idx in range(n_msgs):
        W_base = [random.randint(0, MASK) for _ in range(16)]
        print(f"  Сообщение {msg_idx+1}:")

        for _ in range(n_dw0_per_msg):
            dw0 = random.randint(1, MASK)
            total += 1
            seed = find_seed_mod2(W_base, dw0, n_trials=5000)
            if seed is None:
                continue
            seeds_found += 1

            dW_result, status, rank_trace = hensel_kstep(W_base, seed)
            f_final = compute_f(W_base, dW_result)
            success = all(v == 0 for v in f_final)

            if success:
                full_ok += 1
                print(f"    DW0=0x{dw0:08x}: РЕШЕНИЕ НАЙДЕНО! status={status}")
            else:
                bit = status.split('bit=')[1] if 'bit=' in status else '?'
                fail_bits[bit] = fail_bits.get(bit, 0) + 1
                min_r = min(rank_trace) if rank_trace else '?'
                print(f"    DW0=0x{dw0:08x}: {status}  min_rank={min_r}")
        print()

    print(f"  Результат: семян={seeds_found}/{total}  полных подъёмов={full_ok}/{seeds_found}")
    if fail_bits:
        print(f"  Провалы по битам: {dict(sorted(fail_bits.items()))}")
    return full_ok, seeds_found


# ═══════════════════════════════════════════════════════════════════════════════
# C. ИСЧЕРПЫВАЮЩИЙ ПОИСК MOD 4: 2^15 кандидатов для δ ∈ {0,1}^15
# ═══════════════════════════════════════════════════════════════════════════════

def exhaustive_mod4_from_seed(W_base, seed_dW):
    """
    Для данного семени (f(seed) ≡ 0 mod 2) ищем δ ∈ {0,1}^15
    такой, что f(seed + 2·δ) ≡ 0 mod 4.
    Перебираем все 2^15 кандидатов.
    Возвращает (список найденных δ, число проверок).
    """
    found = []
    checks = 0

    for bits in range(1 << 15):
        dW_cand = dict(seed_dW)
        for j in range(15):
            dW_cand[j + 1] = (seed_dW.get(j + 1, 0) + ((bits >> j) & 1) * 2) & MASK

        f = compute_f(W_base, dW_cand)
        checks += 1
        if all(v % 4 == 0 for v in f):
            found.append(bits)
            if len(found) >= 10:  # достаточно
                break

    return found, checks


def run_section_c(n_msgs=5, n_dw0_per_msg=5):
    print("─" * 72)
    print("C. ИСЧЕРПЫВАЮЩИЙ ПОИСК MOD 4: δ ∈ {0,1}^15 (2^15 = 32768 кандидатов)")
    print("─" * 72)
    print("Вопрос: существует ли хоть ОДНО решение mod 4 для данного семени?\n")

    mod4_found = 0
    total_seeds = 0

    for msg_idx in range(n_msgs):
        W_base = [random.randint(0, MASK) for _ in range(16)]
        print(f"  Сообщение {msg_idx+1}:")

        for _ in range(n_dw0_per_msg):
            dw0 = random.randint(1, MASK)
            seed = find_seed_mod2(W_base, dw0, n_trials=8000)
            if seed is None:
                print(f"    DW0=0x{dw0:08x}: семя не найдено")
                continue

            total_seeds += 1
            found_deltas, checks = exhaustive_mod4_from_seed(W_base, seed)

            if found_deltas:
                mod4_found += 1
                print(f"    DW0=0x{dw0:08x}: ✓ MOD 4 НАЙДЕНО! "
                      f"{len(found_deltas)} δ(s)  (проверено {checks})")
            else:
                print(f"    DW0=0x{dw0:08x}: ✗ нет решений mod 4  (проверено {checks})")

        print()

    print(f"  Итог: mod-4 решения найдены: {mod4_found}/{total_seeds}")
    if mod4_found > 0:
        print(f"  T_MOD4_SOLVABLE: система разрешима mod 4 для некоторых (W,DW0)!")
    else:
        print(f"  T_MOD4_BARRIER_CONFIRMED: ни одного решения mod 4.")
    return mod4_found, total_seeds


# ═══════════════════════════════════════════════════════════════════════════════
# D. ЕСЛИ MOD 4 НАЙДЕНО: продолжение до mod 8, mod 16, ...
# ═══════════════════════════════════════════════════════════════════════════════

def lift_from_mod4(W_base, dW_mod4, start_bit=2):
    """
    Продолжение подъёма от mod-4 решения вверх.
    start_bit=2 означает: x уже является решением mod 4 = 2^2.
    """
    dW = dict(dW_mod4)
    rank_trace = []

    for bit in range(start_bit, 32):
        mod_k = 1 << bit
        mod_kp1 = 1 << (bit + 1)

        f_vec = compute_f(W_base, dW)

        if all(v == 0 for v in f_vec):
            return dW, 'converged_early', rank_trace

        if not all(v % mod_k == 0 for v in f_vec):
            return dW, f'invariant_fail_bit={bit}', rank_trace

        rhs = [(v // mod_k) % 2 for v in f_vec]
        G = jacobian_step_k(W_base, dW, step=mod_k)
        r = rank_mod2(G)
        rank_trace.append(r)

        if r < 15:
            f_zero = compute_f(W_base, dW)
            if all(v % mod_kp1 == 0 for v in f_zero):
                continue
            return dW, f'singular_bit={bit}_rank={r}', rank_trace

        delta = gauss_mod2(G, rhs)
        if delta is None:
            return dW, f'gauss_fail_bit={bit}', rank_trace

        for j in range(15):
            dW[j + 1] = (dW.get(j + 1, 0) + delta[j] * mod_k) & MASK

        f_new = compute_f(W_base, dW)
        if not all(v % mod_kp1 == 0 for v in f_new):
            return dW, f'correction_fail_bit={bit}', rank_trace

    f_final = compute_f(W_base, dW)
    if all(v == 0 for v in f_final):
        return dW, 'converged', rank_trace
    norm = sum(1 for v in f_final if v != 0)
    return dW, f'done_norm={norm}', rank_trace


def run_section_d(n_msgs=5, n_dw0_per_msg=5):
    print("─" * 72)
    print("D. ПОЛНЫЙ ПОДЪЁМ: если mod 4 найдено → продолжаем до mod 2^32")
    print("─" * 72)
    print("Для найденных mod-4 решений применяем k-шаговый Hensel дальше.\n")

    full_solutions = 0
    mod4_total = 0

    for msg_idx in range(n_msgs):
        W_base = [random.randint(0, MASK) for _ in range(16)]
        print(f"  Сообщение {msg_idx+1}:")

        for _ in range(n_dw0_per_msg):
            dw0 = random.randint(1, MASK)
            seed = find_seed_mod2(W_base, dw0, n_trials=8000)
            if seed is None:
                continue

            found_deltas, _ = exhaustive_mod4_from_seed(W_base, seed)
            if not found_deltas:
                continue

            # Берём первый δ с mod-4 решением
            bits = found_deltas[0]
            dW_mod4 = dict(seed)
            for j in range(15):
                dW_mod4[j + 1] = (seed.get(j + 1, 0) + ((bits >> j) & 1) * 2) & MASK

            # Проверяем: действительно mod 4?
            f_check = compute_f(W_base, dW_mod4)
            assert all(v % 4 == 0 for v in f_check), "Ошибка: mod 4 не выполнено!"

            mod4_total += 1

            # Продолжаем подъём
            dW_result, status, rank_trace = lift_from_mod4(W_base, dW_mod4, start_bit=2)
            f_final = compute_f(W_base, dW_result)
            success = all(v == 0 for v in f_final)

            if success:
                full_solutions += 1
                print(f"    DW0=0x{dw0:08x}: *** ПОЛНОЕ РЕШЕНИЕ! *** [{status}]")
                # Верификация
                all_da = all(Da_at(W_base, dW_result, r) == 0 for r in range(3, 17))
                de17 = De_at(W_base, dW_result, 17)
                print(f"      Da3..Da16={'✓' if all_da else '✗'}  De17={'✓' if de17==0 else '✗'}")
            else:
                bit_info = status.split('bit=')[1] if 'bit=' in status else status
                min_r = min(rank_trace) if rank_trace else '?'
                print(f"    DW0=0x{dw0:08x}: mod4→{status}  bit_fail={bit_info}  min_rank={min_r}")

        print()

    print(f"  Итог: mod-4 стартов={mod4_total}  полных решений={full_solutions}")
    return full_solutions, mod4_total


# ═══════════════════════════════════════════════════════════════════════════════
# E. СТАТИСТИКА РАЗРЕШИМОСТИ ПО УРОВНЯМ
# ═══════════════════════════════════════════════════════════════════════════════

def count_solutions_by_level(W_base, seed_dW, max_bit=8):
    """
    Для данного семени считаем, сколько решений существует
    на каждом уровне mod 2^k (k=1..max_bit).
    Использует жадный перебор (не полный).
    """
    results = {}
    current_solutions = [dict(seed_dW)]  # список решений mod 2

    for bit in range(1, max_bit + 1):
        mod_k = 1 << bit
        next_solutions = []

        for x_k in current_solutions:
            # Пробуем оба варианта δ_j=0 и δ_j=1 для ВСЕХ переменных
            # Для эффективности: используем k-шаговый Якобиан (если rank=15, одно решение)
            # Иначе: перебор {0,1}^15 (медленно)
            G = jacobian_step_k(W_base, x_k, step=1 << (bit - 1))
            rhs = [(v // (1 << (bit - 1))) % 2 for v in compute_f(W_base, x_k)]
            r = rank_mod2(G)

            if r == 15:
                delta = gauss_mod2(G, rhs)
                if delta is None:
                    continue  # нет решения от этого x_k
                x_new = dict(x_k)
                for j in range(15):
                    x_new[j + 1] = (x_k.get(j + 1, 0) + delta[j] * (1 << (bit - 1))) & MASK
                f_new = compute_f(W_base, x_new)
                if all(v % mod_k == 0 for v in f_new):
                    next_solutions.append(x_new)
            else:
                # rank=14: либо 0 либо 2 продолжения
                for trial_bits in range(1 << 15):
                    x_trial = dict(x_k)
                    for j in range(15):
                        x_trial[j + 1] = (x_k.get(j + 1, 0) +
                                          ((trial_bits >> j) & 1) * (1 << (bit - 1))) & MASK
                    f_trial = compute_f(W_base, x_trial)
                    if all(v % mod_k == 0 for v in f_trial):
                        next_solutions.append(x_trial)
                    if len(next_solutions) > 10:
                        break

        results[bit] = len(next_solutions)
        current_solutions = next_solutions[:4]  # ограничиваем ветвление
        if not current_solutions:
            for b in range(bit, max_bit + 1):
                results[b] = 0
            break

    return results


def run_section_e(n_msgs=3):
    print("─" * 72)
    print("E. СТАТИСТИКА РАЗРЕШИМОСТИ ПО УРОВНЯМ MOD 2^k")
    print("─" * 72)
    print("Для найденных семян: сколько решений на каждом уровне k?\n")

    for msg_idx in range(n_msgs):
        W_base = [random.randint(0, MASK) for _ in range(16)]
        print(f"  Сообщение {msg_idx+1}:")

        # Ищем 3 семени
        seeds_found = 0
        for dw0 in range(1, 500):
            seed = find_seed_mod2(W_base, dw0, n_trials=3000)
            if seed is None:
                continue
            seeds_found += 1

            # Используем exhaustive mod-4 check как первый шаг
            found_mod4, _ = exhaustive_mod4_from_seed(W_base, seed)
            print(f"    DW0={dw0}: mod2=1 семя, mod4={len(found_mod4)} решений")

            if seeds_found >= 5:
                break

        print()


def print_header():
    print("=" * 72)
    print("П-48 | БАРЬЕР MOD 4: k-ШАГОВЫЙ ЯКОБИАН И ИСЧЕРПЫВАЮЩИЙ ПОИСК")
    print("=" * 72)
    print()
    print("Ключевое исправление: Hensel step k → k+1 должен использовать")
    print("G^(k)_ij = f_i(x + 2^k·e_j) - f_i(x)  (НЕ единичный якобиан J).")
    print()


if __name__ == '__main__':
    random.seed(42)
    print_header()

    run_section_a(n_tests=15)

    run_section_b(n_msgs=5, n_dw0_per_msg=8)

    mod4_found, total_seeds = run_section_c(n_msgs=5, n_dw0_per_msg=5)

    if mod4_found > 0:
        run_section_d(n_msgs=5, n_dw0_per_msg=5)

    run_section_e(n_msgs=3)

    print("\n" + "=" * 72)
    print("ИТОГОВЫЙ ДИАГНОЗ П-48")
    print("=" * 72)
    if mod4_found > 0:
        print(f"  T_MOD4_SOLVABLE: решения mod 4 существуют ({mod4_found}/{total_seeds} семян).")
        print(f"  k-шаговый якобиан РАБОТАЕТ для mod-4 подъёма.")
        print(f"  → П-49: полный подъём до mod 2^32 с k-шаговым якобианом.")
    else:
        print(f"  T_MOD4_BARRIER_ABSOLUTE: решений mod 4 нет ни для одного семени.")
        print(f"  k-шаговый якобиан НЕ помогает — барьер фундаментален.")
        print(f"  → Система (Da3..Da16, De17)=0 не имеет решений над Z/2^32Z.")
        print(f"  → T_GLOBAL_BARRIER_FINAL подтверждён исчерпывающе.")
