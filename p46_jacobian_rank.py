"""
П-46: РАНГ ЯКОБИАНА mod 2 — ФИНАЛЬНЫЙ ДИАГНОЗ.

Вопрос: является ли система (Da3..Da16, De17) = 0 над DW0..DW15
фундаментально несовместной, или есть шанс найти решение при правильном DW0?

Метод: вычислить ранг 15×16 якобиана (по всем 16 переменным DW0..DW15)
над GF(2). Если rank < 15 → система несовместна для "почти всех" правых частей.
Если rank = 15 → глобальный Newton должен сходиться.

Дополнительно: вычислить ранг 15×15 якобиана (только по DW1..DW15)
при фиксированном DW0 — показывает, насколько "управляема" система.
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


def compute_constraints(W_base, dW_dict):
    """15 компонент: Da3..Da16 + De17."""
    result = [Da_at(W_base, dW_dict, r) for r in range(3, 17)]
    result.append(De_at(W_base, dW_dict, 17))
    return result


def jacobian_15x16_mod2(W_base, dW_dict):
    """
    15×16 якобиан по DW0..DW15 над GF(2).
    J[i][j] = (∂f_i/∂DW_j) mod 2, j=0..15.
    """
    base_f = compute_constraints(W_base, dW_dict)
    J = []
    for j in range(16):
        dW_p = dict(dW_dict)
        dW_p[j] = (dW_p.get(j, 0) + 1) & MASK
        perturbed_f = compute_constraints(W_base, dW_p)
        col = [(perturbed_f[i] - base_f[i]) % 2 for i in range(15)]
        J.append(col)
    # J[j] = column j → row-major:
    return [[J[j][i] for j in range(16)] for i in range(15)]


def jacobian_15x15_mod2(W_base, dW_dict):
    """15×15 якобиан по DW1..DW15 над GF(2)."""
    base_f = compute_constraints(W_base, dW_dict)
    J = []
    for j in range(1, 16):
        dW_p = dict(dW_dict)
        dW_p[j] = (dW_p.get(j, 0) + 1) & MASK
        perturbed_f = compute_constraints(W_base, dW_p)
        col = [(perturbed_f[i] - base_f[i]) % 2 for i in range(15)]
        J.append(col)
    return [[J[j][i] for j in range(15)] for i in range(15)]


def rank_mod2(M):
    """Ранг матрицы над GF(2) методом Гаусса."""
    if not M:
        return 0
    nrows = len(M)
    ncols = len(M[0])
    A = [row[:] for row in M]
    rank = 0
    pivot_col = 0
    for row in range(nrows):
        # Найти ненулевой элемент в столбце >= pivot_col
        found = False
        for col in range(pivot_col, ncols):
            for r in range(row, nrows):
                if A[r][col] == 1:
                    A[row], A[r] = A[r], A[row]
                    pivot_col = col
                    found = True
                    break
            if found:
                break
        if not found:
            break
        rank += 1
        for r in range(nrows):
            if r != row and A[r][pivot_col] == 1:
                A[r] = [A[r][k] ^ A[row][k] for k in range(ncols)]
        pivot_col += 1
    return rank


def solvability_mod2(J, b):
    """
    Проверка совместности Jx=b над GF(2).
    Возвращает (solvable, rank_J, rank_aug).
    """
    n = len(b)
    ncols = len(J[0])
    # Расширенная матрица
    aug = [J[i][:] + [b[i] % 2] for i in range(n)]
    rank_J = rank_mod2([row[:ncols] for row in aug])
    rank_aug = rank_mod2([row[:] for row in aug])
    return rank_J == rank_aug, rank_J, rank_aug


# ════════════════════════════════════════════════════════════════════════════
# A. РАНГ 15×16 ЯКОБИАНА (все 16 переменных)
# ════════════════════════════════════════════════════════════════════════════

def run_section_a(n_msgs=20):
    print("─" * 68)
    print("A. РАНГ 15×16 ЯКОБИАНА (DW0..DW15) mod 2")
    print("─" * 68)
    print("Если ранг = 15 → система потенциально разрешима (недоопределённая).")
    print("Если ранг < 15 → система несовместна для почти всех сообщений.\n")
    print(f"  {'#':>3}  {'ранг_15x16':>10}  {'ранг_15x15(DW1..15)':>20}")
    print(f"  {'─'*3}  {'─'*10}  {'─'*20}")

    ranks_16 = []
    ranks_15 = []
    for i in range(n_msgs):
        W_base = [random.randint(0, MASK) for _ in range(16)]
        dW = {}  # нулевые возмущения

        J16 = jacobian_15x16_mod2(W_base, dW)
        J15 = jacobian_15x15_mod2(W_base, dW)

        r16 = rank_mod2(J16)
        r15 = rank_mod2(J15)
        ranks_16.append(r16)
        ranks_15.append(r15)
        print(f"  {i+1:>3}  {r16:>10}  {r15:>20}")

    print(f"\n  Среднее ранг 15×16: {sum(ranks_16)/n_msgs:.1f}/15")
    print(f"  Среднее ранг 15×15: {sum(ranks_15)/n_msgs:.1f}/15")
    print(f"  Ранг = 15 (полный): {sum(1 for r in ranks_16 if r==15)}/{n_msgs}")
    return ranks_16, ranks_15


# ════════════════════════════════════════════════════════════════════════════
# B. СОВМЕСТНОСТЬ СИСТЕМЫ mod 2 для разных DW0
#    При фиксированном W_base, варьируем DW0 → ищем DW0, при котором
#    система Jx = -f над GF(2) СОВМЕСТНА
# ════════════════════════════════════════════════════════════════════════════

def run_section_b(n_msgs=5, n_dw0=512):
    print("\n" + "─" * 68)
    print(f"B. СОВМЕСТНОСТЬ J·Δ = -f mod 2 при вариации DW0 (n={n_dw0})")
    print("─" * 68)
    print("Ключевой вопрос: существует ли DW0, при котором система над GF(2)")
    print("имеет решение?\n")

    for msg_idx in range(n_msgs):
        W_base = [random.randint(0, MASK) for _ in range(16)]
        n_solvable = 0
        n_full_rank = 0

        for dw0_idx in range(n_dw0):
            dW = {0: dw0_idx + 1}
            f = compute_constraints(W_base, dW)
            # Якобиан 15×15 по DW1..DW15
            J15 = jacobian_15x15_mod2(W_base, dW)
            rhs = [(-v) % 2 for v in f]

            solvable, rj, ra = solvability_mod2(J15, rhs)
            if solvable:
                n_solvable += 1
            if rj == 15:
                n_full_rank += 1

        print(f"  Сообщение {msg_idx+1}:")
        print(f"    Полный ранг J (=15): {n_full_rank}/{n_dw0} ({100*n_full_rank/n_dw0:.0f}%)")
        print(f"    Система совместна:   {n_solvable}/{n_dw0} ({100*n_solvable/n_dw0:.0f}%)")
        print()


# ════════════════════════════════════════════════════════════════════════════
# C. ГЛОБАЛЬНЫЙ NEWTON с перебором DW0
#    Для каждого DW0 из {1..N}: один шаг Newton mod 2, проверяем норму
# ════════════════════════════════════════════════════════════════════════════

def gauss_mod2_solve(J, b):
    """Решить Jx=b над GF(2). Возвращает x или None."""
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


def newton_mod2_from_zero(W_base, dw0, max_iter=32):
    """
    Итерации Newton mod 2, начиная с DW1..DW15 = 0.
    Возвращает финальную норму (число ненулевых f_i mod 2).
    """
    dW = {0: dw0}
    for _ in range(max_iter):
        f = compute_constraints(W_base, dW)
        if all(v % 2 == 0 for v in f):
            return 0
        J15 = jacobian_15x15_mod2(W_base, dW)
        rhs = [(-v) % 2 for v in f]
        delta = gauss_mod2_solve(J15, rhs)
        if delta is None:
            break
        for j in range(15):
            dW[j + 1] = (dW.get(j + 1, 0) + delta[j]) & MASK
    f_final = compute_constraints(W_base, dW)
    return sum(1 for v in f_final if v % 2 != 0)


def run_section_c(n_msgs=3, n_dw0=256):
    print("─" * 68)
    print(f"C. NEWTON mod 2 с перебором DW0 (n_dw0={n_dw0})")
    print("─" * 68)
    print("Ищем DW0 при котором Newton mod 2 сходится к норме=0.\n")

    for msg_idx in range(n_msgs):
        W_base = [random.randint(0, MASK) for _ in range(16)]
        best_norm = 15
        best_dw0 = None
        n_zero_norm = 0
        norm_hist = {}

        for dw0_idx in range(n_dw0):
            dw0 = dw0_idx + 1
            norm = newton_mod2_from_zero(W_base, dw0, max_iter=16)
            norm_hist[norm] = norm_hist.get(norm, 0) + 1
            if norm < best_norm:
                best_norm = norm
                best_dw0 = dw0
            if norm == 0:
                n_zero_norm += 1

        print(f"  Сообщение {msg_idx+1}:")
        print(f"    Лучшая норма: {best_norm}  (DW0=0x{best_dw0:08x})")
        print(f"    Норма=0:      {n_zero_norm}/{n_dw0}")
        hist_str = "  ".join(f"norm={k}:{v}" for k, v in sorted(norm_hist.items()))
        print(f"    Распределение: {hist_str}")
        print()


# ════════════════════════════════════════════════════════════════════════════
# D. СТРУКТУРА ЯКОБИАНА: какие переменные влияют на какие уравнения?
# ════════════════════════════════════════════════════════════════════════════

def run_section_d(n_msgs=5):
    print("─" * 68)
    print("D. ПАТТЕРН ЯКОБИАНА mod 2 (визуализация)")
    print("─" * 68)
    print("J[i][j] = ∂(Da_{i+3} или De17)/∂DW_{j+1}")
    print("Строки: уравнения (Da3..Da16, De17), Столбцы: DW1..DW15\n")

    W_base = [random.randint(0, MASK) for _ in range(16)]
    dW = {}
    J15 = jacobian_15x15_mod2(W_base, dW)

    labels_row = [f"Da{r}" for r in range(3, 17)] + ["De17"]
    labels_col = [f"W{j}" for j in range(1, 16)]

    print(f"       " + " ".join(f"{c:>4}" for c in labels_col))
    for i, row in enumerate(J15):
        print(f"  {labels_row[i]:>5}: " + " ".join("  ██" if v else "  ░░" for v in row))

    r = rank_mod2(J15)
    print(f"\n  Ранг = {r}/15")

    # Нулевые столбцы (DW_j который ни на что не влияет)
    zero_cols = [j for j in range(15) if all(J15[i][j] == 0 for i in range(15))]
    zero_rows = [i for i in range(15) if all(J15[i][j] == 0 for j in range(15))]
    if zero_cols:
        print(f"  Нулевые столбцы (DW без влияния): {[labels_col[j] for j in zero_cols]}")
    if zero_rows:
        print(f"  Нулевые строки (уравнения без зависимостей): {[labels_row[i] for i in zero_rows]}")


def main():
    random.seed(42)
    print("=" * 68)
    print("П-46 | РАНГ ЯКОБИАНА mod 2 — ФИНАЛЬНЫЙ ДИАГНОЗ")
    print("=" * 68)
    print()

    ranks_16, ranks_15 = run_section_a(n_msgs=20)

    run_section_b(n_msgs=3, n_dw0=256)

    run_section_c(n_msgs=3, n_dw0=128)

    run_section_d(n_msgs=1)

    print("\n" + "=" * 68)
    print("ИТОГ П-46")
    print("=" * 68)

    avg16 = sum(ranks_16) / len(ranks_16)
    avg15 = sum(ranks_15) / len(ranks_15)
    full16 = sum(1 for r in ranks_16 if r == 15)

    print(f"""
Ранг 15×16 якобиана (DW0..DW15) mod 2: среднее = {avg16:.1f}/15
Ранг 15×15 якобиана (DW1..DW15)  mod 2: среднее = {avg15:.1f}/15
Случаев полного ранга (15/15): {full16}/{len(ranks_16)}
""")

    if avg15 < 14:
        print("""  T_JACOBIAN_DEFICIENT: якобиан 15×15 систематически вырожден mod 2.
  Система Da3..Da16,De17=0 фундаментально несовместна.
  Барьер подтверждён: 2^(32*(16-rank)) сложность поиска коллизии.
""")
    elif full16 == len(ranks_16):
        print("""  T_FULL_RANK: якобиан 15×16 имеет полный ранг!
  Система потенциально разрешима при правильном DW0.
  Нужен поиск DW0 за O(2^(32*(16-15))) = O(2^32).
""")
    else:
        print(f"""  T_PARTIAL_RANK: ранг = {avg16:.1f} < 15 в среднем.
  Для {full16}/{len(ranks_16)} сообщений полный ранг — есть шанс найти решение.
""")


if __name__ == '__main__':
    main()
