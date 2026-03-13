"""
П-47: HENSEL LIFTING — подъём решения mod 2 → mod 2^32.

Оптимизация: sha256_trace вычисляет все раунды за один проход,
вместо отдельного вызова для каждого раунда r=3..17.
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


def expand_schedule(W16):
    W = list(W16)
    for i in range(16, 64):
        W.append((sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16]) & MASK)
    return W


def sha256_trace(W16):
    """
    Вычисляет трассу SHA-256 за один проход (17 раундов).
    Возвращает list[(a_r, e_r)] для r=1..17 (индекс 0 = после раунда 1).
    """
    W = expand_schedule(W16)
    a, b, c, d, e, f, g, h = H0
    trace_a = []
    trace_e = []
    for r in range(17):
        T1 = (h + Sig1(e) + Ch(e, f, g) + K_SHA[r] + W[r]) & MASK
        T2 = (Sig0(a) + Maj(a, b, c)) & MASK
        h = g; g = f; f = e; e = (d + T1) & MASK
        d = c; c = b; b = a; a = (T1 + T2) & MASK
        trace_a.append(a)
        trace_e.append(e)
    return trace_a, trace_e


def apply_dW(W_base, dW):
    W2 = list(W_base)
    for idx, delta in dW.items():
        W2[idx] = (W2[idx] + delta) & MASK
    return W2


def compute_constraints(W_base, dW, base_trace=None):
    """
    15 значений: Da_r = a_r(W+dW) - a_r(W) для r=3..16, De17 = e17 diff.
    base_trace: (trace_a, trace_e) для W_base — кэш для ускорения.
    """
    if base_trace is None:
        base_a, base_e = sha256_trace(W_base)
    else:
        base_a, base_e = base_trace

    W2 = apply_dW(W_base, dW)
    pert_a, pert_e = sha256_trace(W2)

    # Da_r = a после r раундов. r=3 → index 2; r=16 → index 15
    result = [(pert_a[r - 1] - base_a[r - 1]) & MASK for r in range(3, 17)]
    result.append((pert_e[16] - base_e[16]) & MASK)  # De17 = e после 17 раундов
    return result


def jacobian_and_f(W_base, dW):
    """
    Вычисляет f и 15×15 якобиан по DW1..DW15 mod 2.
    Один базовый trace + 15 возмущённых → итого 16 SHA trace.
    """
    base_trace = sha256_trace(W_base)
    f = compute_constraints(W_base, dW, base_trace)

    # Вычислить 15 возмущённых
    J_cols = []
    for j in range(15):
        dW_p = dict(dW)
        dW_p[j + 1] = (dW_p.get(j + 1, 0) + 1) & MASK
        W2 = apply_dW(W_base, dW_p)
        pert_a, pert_e = sha256_trace(W2)
        base_a, base_e = base_trace
        col = []
        for r in range(3, 17):
            col.append((pert_a[r-1] - base_a[r-1]) % 2)
        col.append((pert_e[16] - base_e[16]) % 2)
        J_cols.append(col)

    J = [[J_cols[j][i] for j in range(15)] for i in range(15)]
    return J, f


def gauss_mod2(J, b):
    """Решение Jx=b над GF(2). Возвращает (x, rank) или (None, rank)."""
    n = len(b)
    aug = [[J[i][j] for j in range(n)] + [b[i] % 2] for i in range(n)]
    rank = 0
    pivot_info = []
    pivot_row = 0
    for col in range(n):
        found = -1
        for row in range(pivot_row, n):
            if aug[row][col] == 1:
                found = row
                break
        if found == -1:
            continue
        aug[pivot_row], aug[found] = aug[found], aug[pivot_row]
        for row in range(n):
            if row != pivot_row and aug[row][col] == 1:
                for k in range(n + 1):
                    aug[row][k] ^= aug[pivot_row][k]
        pivot_info.append((pivot_row, col))
        rank += 1
        pivot_row += 1
    if rank < n:
        return None, rank
    x = [0] * n
    for row, col in pivot_info:
        x[col] = aug[row][n]
    return x, rank


def find_dw0_mod2(W_base, max_search=512, max_newton=16):
    """
    Поиск DW0 ∈ [1..max_search] где система (Da3..Da16,De17)=0 над GF(2)
    достижима итеративным Newton (до max_newton шагов).
    Возвращает (dw0, dW) или (None, None).
    """
    for dw0 in range(1, max_search + 1):
        dW = {0: dw0}
        for _ in range(max_newton):
            J, f = jacobian_and_f(W_base, dW)
            if all(v % 2 == 0 for v in f):
                return dw0, dict(dW)
            rhs = [(-v) % 2 for v in f]
            x, rank = gauss_mod2(J, rhs)
            if x is None:
                break  # вырожденный якобиан — пробуем следующий DW0
            # XOR-накопление поправки (Newton mod 2)
            for j in range(15):
                dW[j + 1] = (dW.get(j + 1, 0) + x[j]) & MASK
    return None, None


def hensel_step(W_base, dW, k):
    """
    Шаг Хенселя: поднять x (≡0 mod 2^k) до mod 2^{k+1}.
    Ищет δ ∈ {0,1}^15: J·δ ≡ -f(x)/2^k (mod 2).
    """
    mod_k = 1 << k
    mod_k1 = 1 << (k + 1)

    J, f = jacobian_and_f(W_base, dW)

    if any(v % mod_k != 0 for v in f):
        return dW, False, -1, 'precondition_fail'

    rhs = [(-(v >> k)) % 2 for v in f]
    delta, rank = gauss_mod2(J, rhs)
    if delta is None:
        return dW, False, rank, f'singular_rank={rank}'

    new_dW = dict(dW)
    for j in range(15):
        new_dW[j + 1] = (new_dW.get(j + 1, 0) + delta[j] * mod_k) & MASK

    base_trace = sha256_trace(W_base)
    f_new = compute_constraints(W_base, new_dW, base_trace)
    bad = sum(1 for v in f_new if v % mod_k1 != 0)
    if bad > 0:
        return dW, False, rank, f'verify_fail_{bad}_eqs'

    return new_dW, True, rank, 'ok'


def full_hensel_lift(W_base, dW_start, verbose=True):
    """Полный подъём mod 2 → mod 2^32."""
    dW = dict(dW_start)
    for k in range(1, 32):
        new_dW, ok, rank, msg = hensel_step(W_base, dW, k)
        if not ok:
            if verbose:
                print(f"    k={k:2d}: FAIL [{msg}]")
            return dW, False, k
        dW = new_dW
        if verbose and (k <= 3 or k % 8 == 0 or k == 31):
            f = compute_constraints(W_base, dW)
            norm_k1 = sum(1 for v in f if v % (1 << (k + 1)) != 0)
            print(f"    k={k:2d}: ok  rank={rank}  норма_mod_2^{k+1}={norm_k1}")
    f_final = compute_constraints(W_base, dW)
    return dW, all(v == 0 for v in f_final), 32


# ─────────────────────────────────────────────────────────────────────────────

def run_section_a(n_msgs=5):
    print("─" * 64)
    print("A. ПОЛНЫЙ HENSEL LIFT: mod 2 → mod 2^32")
    print("─" * 64)

    total_mod2 = 0
    total_ok = 0

    for i in range(n_msgs):
        W_base = [random.randint(0, MASK) for _ in range(16)]
        print(f"\n  Сообщение {i+1}:")

        dw0, dW = find_dw0_mod2(W_base, max_search=512)
        if dW is None:
            print(f"    DW0 не найден за 512 итераций")
            continue

        total_mod2 += 1
        print(f"    DW0=0x{dw0:08x} — решение mod 2 найдено")

        dW_final, success, k_stop = full_hensel_lift(W_base, dW, verbose=True)

        if success:
            total_ok += 1
            f = compute_constraints(W_base, dW_final)
            print(f"    ★ УСПЕХ! Все 15 ур-й = 0: {'✓' if all(v==0 for v in f) else '✗'}")
            dw_vals = [dW_final.get(j, 0) for j in range(16)]
            print(f"    DW[0..7]:  {' '.join(f'{v:08x}' for v in dw_vals[:8])}")
            print(f"    DW[8..15]: {' '.join(f'{v:08x}' for v in dw_vals[8:])}")
        else:
            print(f"    FAIL на шаге k={k_stop}")

    print(f"\n  Итог A: mod2={total_mod2}/{n_msgs}, полный подъём={total_ok}/{n_msgs}")
    return total_ok, n_msgs


def run_section_b(n_msgs=20):
    print("\n" + "─" * 64)
    print("B. СТАТИСТИКА (20 сообщений, тихий режим)")
    print("─" * 64)

    fail_at = {}
    success = 0
    no_mod2 = 0

    for i in range(n_msgs):
        W_base = [random.randint(0, MASK) for _ in range(16)]
        dw0, dW = find_dw0_mod2(W_base, max_search=256)
        if dW is None:
            no_mod2 += 1
            continue
        dW_final, ok, k_stop = full_hensel_lift(W_base, dW, verbose=False)
        if ok:
            success += 1
        else:
            fail_at[k_stop] = fail_at.get(k_stop, 0) + 1

    print(f"  DW0 не найден (mod 2): {no_mod2}/{n_msgs}")
    print(f"  Успехов (полный подъём): {success}/{n_msgs}")
    if fail_at:
        print(f"  Тупики по шагам k: {dict(sorted(fail_at.items()))}")
    return success, n_msgs


def main():
    random.seed(42)
    print("=" * 64)
    print("П-47 | HENSEL LIFTING: mod 2 → mod 2^32")
    print("=" * 64)
    print()

    ok_a, n_a = run_section_a(n_msgs=5)
    ok_b, n_b = run_section_b(n_msgs=20)

    print("\n" + "=" * 64)
    print("ИТОГ П-47")
    print("=" * 64)

    if ok_a + ok_b > 0:
        print(f"""
  ПРОРЫВ: Hensel lifting работает!
  Секция A: {ok_a}/{n_a}, Секция B: {ok_b}/{n_b}

  Алгоритм решения:
    1. Перебрать DW0 ∈ [1..512] — O(512 × 16 SHA-traces)
    2. Hensel: 31 шаг × 16 SHA-traces = O(496 SHA-traces)
  Итого: ~8000 SHA вычислений на сообщение.
""")
    else:
        print(f"""
  Hensel lifting не сошёлся: {ok_a+ok_b}/{n_a}.
  Якобиан теряет ранг на промежуточных шагах.
  T_GLOBAL_BARRIER подтверждён: SHA-256 устойчив на каждом уровне.
""")


if __name__ == '__main__':
    main()
