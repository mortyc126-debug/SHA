"""
П-46: ПРОВЕРКА ГИПОТЕЗЫ DW0.

Диагноз П-45 (секция D):
  Якобиан 15×15 сингулярен mod 2 при ФИКСИРОВАННОМ DW0.
  → Система Da3..Da16=0, De17=0 не имеет решений для большинства DW0.
  → DW0 — единственная степень свободы, которая управляет разрешимостью.

Гипотеза:
  Существует ~1/2^15 значений DW0 (из 2^32), при которых якобиан
  15×15 НЕВЫРОЖДЕН mod 2 → система разрешима.
  Поиск: O(2^15) перебор DW0, а не O(2^64).

Структура П-46:
  A. Empirical: для N случайных DW0 считаем rank(J mod 2).
     Проверяем: rank < 15 почти всегда? Какой типичный дефект?
  B. Search experiment: перебираем DW0 = 1..2^16, ищем ранг = 15.
     Считаем плотность "удачных" DW0.
  C. Verification: для найденных DW0 запускаем глобальный Newton
     (из П-45) и проверяем сходимость.
  D. Bit-profile: для удачных DW0 — есть ли битовые паттерны?
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
    """
    15×15 якобиан: J[i][j] = ∂f_i / ∂DW_{j+1}.
    Переменные: DW1..DW15 (индексы j=0..14 → DW1..DW15).
    """
    base_f = compute_all_constraints(W_base, dW_dict)
    cols = []
    for j in range(15):
        var_idx = j + 1
        dW_p = dict(dW_dict)
        dW_p[var_idx] = (dW_p.get(var_idx, 0) + 1) & MASK
        pf = compute_all_constraints(W_base, dW_p)
        cols.append([(pf[i] - base_f[i]) & MASK for i in range(15)])
    # cols[j] = j-th column
    return [[cols[j][i] for j in range(15)] for i in range(15)]


def rank_mod2(J):
    """Ранг матрицы 15×15 над GF(2)."""
    n = len(J)
    m = len(J[0])
    A = [[J[i][j] % 2 for j in range(m)] for i in range(n)]
    rank = 0
    pivot_col = 0
    for row in range(n):
        # Найти ненулевой элемент в текущем столбце
        found = -1
        for c in range(pivot_col, m):
            for r in range(row, n):
                if A[r][c] == 1:
                    found = (r, c)
                    break
            if found != -1:
                break
        if found == -1:
            break
        r, c = found
        A[row], A[r] = A[r], A[row]
        pivot_col = c + 1
        for r2 in range(n):
            if r2 != row and A[r2][c] == 1:
                for k in range(m):
                    A[r2][k] ^= A[row][k]
        rank += 1
    return rank


def gauss_mod2(J, b):
    """Решение Jx = b над GF(2). Возвращает x или None."""
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


def global_newton_iterate(W_base, dW0, max_iter=64):
    """
    Глобальный Newton для 15-переменной системы (DW1..DW15).
    DW0 = dW0 фиксировано.
    Побитовый подъём: на каждом бите k решаем систему mod 2.
    """
    dW = {0: dW0}
    for j in range(1, 16):
        dW[j] = 0

    for bit in range(32):
        f_vec = compute_all_constraints(W_base, dW)
        norm = sum(1 for v in f_vec if v != 0)
        if norm == 0:
            return dW, 'converged_early'

        J = numerical_jacobian_15x15(W_base, dW)

        # Проверяем ранг
        r = rank_mod2(J)
        if r < 15:
            return dW, f'singular_bit={bit}_rank={r}'

        rhs = [(-v) % 2 for v in f_vec]
        delta = gauss_mod2(J, rhs)
        if delta is None:
            return dW, f'gauss_fail_bit={bit}'

        for j in range(15):
            var_idx = j + 1
            dW[var_idx] = (dW.get(var_idx, 0) + delta[j] * (1 << bit)) & MASK

    f_final = compute_all_constraints(W_base, dW)
    if all(v == 0 for v in f_final):
        return dW, 'converged'
    norm = sum(1 for v in f_final if v != 0)
    return dW, f'failed_norm={norm}'


# ═══════════════════════════════════════════════════════════════════════════════
# A. РАНГ ЯКОБИАНА КАК ФУНКЦИЯ DW0
# ═══════════════════════════════════════════════════════════════════════════════

def run_section_a(n_msgs=5, n_dw0_per_msg=200):
    print("─" * 72)
    print("A. РАНГ ЯКОБИАНА 15×15 mod 2 КАК ФУНКЦИЯ DW0")
    print("─" * 72)
    print("Гипотеза: rank < 15 для большинства DW0. Дефект = 15 - rank.\n")

    rank_hist = {}  # rank → count
    total = 0

    for msg_idx in range(n_msgs):
        W_base = [random.randint(0, MASK) for _ in range(16)]
        ranks_this = []

        for _ in range(n_dw0_per_msg):
            dw0 = random.randint(1, MASK)
            dW = {0: dw0}
            for j in range(1, 16):
                dW[j] = 0
            J = numerical_jacobian_15x15(W_base, dW)
            r = rank_mod2(J)
            ranks_this.append(r)
            rank_hist[r] = rank_hist.get(r, 0) + 1
            total += 1

        mean_r = sum(ranks_this) / len(ranks_this)
        full = sum(1 for r in ranks_this if r == 15)
        print(f"  Сообщение {msg_idx+1}: mean_rank={mean_r:.2f}  rank=15: {full}/{n_dw0_per_msg}")

    print(f"\n  Распределение рангов по {total} тестам:")
    for r in sorted(rank_hist):
        pct = 100.0 * rank_hist[r] / total
        bar = '█' * int(pct / 2)
        print(f"    rank={r:2d}: {rank_hist[r]:5d}  ({pct:5.1f}%)  {bar}")

    full_rank = rank_hist.get(15, 0)
    print(f"\n  Доля rank=15: {full_rank}/{total} = {100.0*full_rank/total:.2f}%")
    expected_exp = math.log2(total / full_rank) if full_rank > 0 else float('inf')
    print(f"  Ожидаемый поиск DW0: ~2^{expected_exp:.1f} попыток")
    return rank_hist, total


# ═══════════════════════════════════════════════════════════════════════════════
# B. ПОИСК: ПЕРЕБОР DW0 = 1..2^N, ИЩЕМ RANK=15
# ═══════════════════════════════════════════════════════════════════════════════

def run_section_b(n_msgs=3, search_limit=65536):
    print("\n" + "─" * 72)
    print(f"B. ПОИСК DW0 С RANK=15  (перебор до {search_limit} = 2^{math.log2(search_limit):.0f})")
    print("─" * 72)
    print("Для каждого сообщения: ищем первый DW0 с невырожденным якобианом.\n")

    for msg_idx in range(n_msgs):
        W_base = [random.randint(0, MASK) for _ in range(16)]
        found_list = []
        gaps = []
        prev = 0

        for dw0 in range(1, search_limit + 1):
            dW = {0: dw0}
            for j in range(1, 16):
                dW[j] = 0
            J = numerical_jacobian_15x15(W_base, dW)
            r = rank_mod2(J)
            if r == 15:
                found_list.append(dw0)
                gaps.append(dw0 - prev)
                prev = dw0
                if len(found_list) >= 20:
                    break

        if found_list:
            mean_gap = sum(gaps) / len(gaps) if gaps else 0
            print(f"  Сообщение {msg_idx+1}: найдено {len(found_list)} DW0 до {found_list[-1]}")
            print(f"    Первые DW0: {found_list[:8]}")
            print(f"    Средний интервал: {mean_gap:.1f}  (~2^{math.log2(mean_gap):.1f})")
            density = len(found_list) / found_list[-1]
            print(f"    Плотность: {density:.4f}  (~1/2^{-math.log2(density):.1f})")
        else:
            print(f"  Сообщение {msg_idx+1}: rank=15 НЕ найден в первых {search_limit}!")
        print()


# ═══════════════════════════════════════════════════════════════════════════════
# C. ВЕРИФИКАЦИЯ: NEWTON ДЛЯ НАЙДЕННЫХ DW0
# ═══════════════════════════════════════════════════════════════════════════════

def run_section_c(n_msgs=5, search_limit=32768):
    print("─" * 72)
    print("C. ВЕРИФИКАЦИЯ: NEWTON ДЛЯ УДАЧНЫХ DW0 (rank=15)")
    print("─" * 72)
    print("Если rank=15 → система разрешима. Проверяем сходимость Newton.\n")

    newton_ok = 0
    newton_total = 0

    for msg_idx in range(n_msgs):
        W_base = [random.randint(0, MASK) for _ in range(16)]
        found_dw0 = None

        # Ищем DW0 с rank=15
        for dw0 in range(1, search_limit + 1):
            dW = {0: dw0}
            for j in range(1, 16):
                dW[j] = 0
            J = numerical_jacobian_15x15(W_base, dW)
            r = rank_mod2(J)
            if r == 15:
                found_dw0 = dw0
                break

        if found_dw0 is None:
            print(f"  Сообщение {msg_idx+1}: DW0 с rank=15 не найден в {search_limit}")
            continue

        # Запускаем Newton
        dW_result, status = global_newton_iterate(W_base, found_dw0, max_iter=64)
        f_final = compute_all_constraints(W_base, dW_result)
        success = all(v == 0 for v in f_final)
        norm = sum(1 for v in f_final if v != 0)

        newton_total += 1
        if success:
            newton_ok += 1
            # Верификация: проверим Da3..Da16 и De17 все = 0
            all_da = all(Da_at(W_base, dW_result, r) == 0 for r in range(3, 17))
            de17 = De_at(W_base, dW_result, 17)
            print(f"  Сообщение {msg_idx+1}: DW0=0x{found_dw0:08x}  Newton=OK  "
                  f"Da3..16={'✓' if all_da else '✗'}  De17={'✓' if de17==0 else '✗'}  [{status}]")
        else:
            print(f"  Сообщение {msg_idx+1}: DW0=0x{found_dw0:08x}  Newton=FAIL  "
                  f"norm={norm}  [{status}]")

    print(f"\n  Newton успех: {newton_ok}/{newton_total}")
    return newton_ok, newton_total


# ═══════════════════════════════════════════════════════════════════════════════
# D. БИТОВЫЙ ПРОФИЛЬ УДАЧНЫХ DW0
# ═══════════════════════════════════════════════════════════════════════════════

def run_section_d(n_msgs=3, search_limit=65536, max_found=100):
    print("\n" + "─" * 72)
    print("D. БИТОВЫЙ ПРОФИЛЬ УДАЧНЫХ DW0")
    print("─" * 72)
    print("Есть ли структура в DW0 с rank=15? Битовые частоты vs равномерное.\n")

    all_found = []

    for msg_idx in range(n_msgs):
        W_base = [random.randint(0, MASK) for _ in range(16)]
        found = []

        for dw0 in range(1, search_limit + 1):
            dW = {0: dw0}
            for j in range(1, 16):
                dW[j] = 0
            J = numerical_jacobian_15x15(W_base, dW)
            r = rank_mod2(J)
            if r == 15:
                found.append(dw0)
                if len(found) >= max_found:
                    break

        all_found.extend(found)
        print(f"  Сообщение {msg_idx+1}: найдено {len(found)} удачных DW0 до {search_limit}")

    if not all_found:
        print("  Нет данных для анализа.")
        return

    # Битовые частоты для малых бит (0..15)
    print(f"\n  Частоты бит (0..15) для {len(all_found)} удачных DW0:")
    print(f"  {'Бит':<5} {'Частота':<10} {'Ожидание':<10} {'Отклонение'}")
    total = len(all_found)
    expected = total / 2
    for bit in range(16):
        cnt = sum(1 for v in all_found if (v >> bit) & 1)
        dev = (cnt - expected) / (total ** 0.5)
        print(f"  {bit:<5} {cnt:<10} {expected:<10.0f} {dev:+.2f}σ")

    # Распределение по модулям
    print(f"\n  DW0 mod 8 распределение:")
    mod8 = {}
    for v in all_found:
        k = v % 8
        mod8[k] = mod8.get(k, 0) + 1
    for k in range(8):
        cnt = mod8.get(k, 0)
        bar = '█' * int(20 * cnt / total)
        print(f"    DW0 mod 8 = {k}: {cnt:4d}  {bar}")


# ═══════════════════════════════════════════════════════════════════════════════
# E. ИТОГ: ОЦЕНКА СЛОЖНОСТИ ПОИСКА
# ═══════════════════════════════════════════════════════════════════════════════

def run_section_e(rank_hist, total):
    print("\n" + "─" * 72)
    print("E. ИТОГОВАЯ ОЦЕНКА И УТОЧНЁННЫЙ ДИАГНОЗ")
    print("─" * 72)

    full_rank = rank_hist.get(15, 0)
    if full_rank == 0:
        print("  rank=15 не наблюдался! Гипотеза о 2^15 — НЕВЕРНА.")
        return

    fraction = full_rank / total
    log2_inv = -math.log2(fraction)

    print(f"  Доля DW0 с rank=15 при старте: {100*fraction:.1f}%  (~1/2^{log2_inv:.2f})")
    print()
    print(f"  Секция C показала: даже при rank=15 при старте,")
    print(f"  Newton ПАДАЕТ на бите 1-2 (rank J падает до 14).")
    print()
    print(f"  Уточнённый диагноз:")
    print(f"    rank=15 при DW=0 — НЕОБХОДИМОЕ, но НЕ ДОСТАТОЧНОЕ условие.")
    print(f"    Нужен DW0 такой, что rank(J) = 15 на ВСЕХ 32 битах подъёма.")
    print(f"    Это более сильное условие → реальная плотность ещё ниже.")
    print()
    print(f"  → Для П-47: искать DW0 с rank=15 через ВСЕ биты (Hensel-stable DW0).")
    print(f"    Такой DW0 и есть истинный «ключ» к решению системы.")


def print_header():
    print("=" * 72)
    print("П-46 | ГИПОТЕЗА DW0: ПОИСК В ПРОСТРАНСТВЕ ~2^15")
    print("=" * 72)
    print()
    print("Ключевой вопрос: доля DW0 с невырожденным якобианом 15×15 mod 2.")
    print("Если ~1/2^15 → поиск коллизии за O(2^15), а не O(2^64).")
    print()


if __name__ == '__main__':
    random.seed(42)
    print_header()

    rank_hist, total = run_section_a(n_msgs=5, n_dw0_per_msg=200)

    run_section_b(n_msgs=3, search_limit=65536)

    run_section_c(n_msgs=5, search_limit=32768)

    run_section_d(n_msgs=3, search_limit=32768, max_found=50)

    run_section_e(rank_hist, total)
