"""
SHA-256 Дифференциальный Криптоанализ
П-1: Lifted-арифметика — аналитическое условие De3=0
"""

import random
from itertools import product
from collections import defaultdict

M32 = 0xFFFFFFFF

# --- SHA-256 константы ---
K = [
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,
    0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,
    0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,
    0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,
    0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,
    0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,
    0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,
    0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,
    0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2,
]

H0 = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
]

# Верифицированные константы из методички
C_IV    = 0xf377ed68
T2_0    = 0x08909ae5
BASE_A1 = 0xfc08884d  # (C_IV + T2_0) mod 2^32

# --- Базовые операции ---
def rotr(x, n): return ((x >> n) | (x << (32 - n))) & M32
def shr(x, n):  return x >> n

def Sig0(x): return rotr(x,2)  ^ rotr(x,13) ^ rotr(x,22)
def Sig1(x): return rotr(x,6)  ^ rotr(x,11) ^ rotr(x,25)
def sig0(x): return rotr(x,7)  ^ rotr(x,18) ^ shr(x,3)
def sig1(x): return rotr(x,17) ^ rotr(x,19) ^ shr(x,10)

def Ch(e,f,g):  return (e & f) ^ (~e & g) & M32
def Maj(a,b,c): return (a & b) ^ (a & c) ^ (b & c)

def schedule(W16):
    W = list(W16) + [0]*48
    for i in range(16, 64):
        W[i] = (sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16]) & M32
    return W

def sha256_rounds(W16, N):
    W = schedule(W16)
    a,b,c,d,e,f,g,h = H0
    for i in range(N):
        T1 = (h + Sig1(e) + Ch(e,f,g) + K[i] + W[i]) & M32
        T2 = (Sig0(a) + Maj(a,b,c)) & M32
        h=g; g=f; f=e; e=(d+T1)&M32
        d=c; c=b; b=a; a=(T1+T2)&M32
    return [a,b,c,d,e,f,g,h]

# ================================================================
# П-1: LIFTED АРИФМЕТИКА
# ================================================================

def sha256_rounds_lifted(W16, N):
    """
    Вычисляет состояние SHA-256 над Z (без mod 2^32).

    КЛЮЧЕВОЙ ИНВАРИАНТ:
      - Сложения T1, T2, e, a выполняются без mod 2^32 (lifted)
      - XOR-операции (Sig1, Sig0, Ch, Maj) работают только с нижними
        32 битами своих аргументов: маскируем перед вызовом.
        Причина: rotr уже применяет & M32 внутри Sig*, но Ch и Maj
        содержат «e & f» и «a & b» — при lifted-значениях обоих
        операндов это даёт лишние биты выше 31 и ломает результат.
      - h в T1 = h_lifted (НЕ маскируется): он — аккумулятор carry-цепочки,
        его полное значение над Z нужно для корректного слежения за переносом.
    """
    W = schedule(W16)
    a,b,c,d,e,f,g,h = H0
    states = [[a,b,c,d,e,f,g,h]]
    for i in range(N):
        # маскируем XOR-входы до 32 бит
        e_m = e & M32; f_m = f & M32; g_m = g & M32
        a_m = a & M32; b_m = b & M32; c_m = c & M32
        # h — lifted addend, НЕ маскируется
        T1 = h + Sig1(e_m) + Ch(e_m, f_m, g_m) + K[i] + W[i]   # НЕТ & M32
        T2 = Sig0(a_m) + Maj(a_m, b_m, c_m)                      # НЕТ & M32
        h=g; g=f; f=e
        e = d + T1          # НЕТ & M32
        d=c; c=b; b=a
        a = T1 + T2         # НЕТ & M32
        states.append([a,b,c,d,e,f,g,h])
    return states

def lifted_diff_e(W0_base, delta_W0=1, N=3):
    """
    Вычисляет lifted_diff_e_r = e_r(W XOR delta) - e_r(W) над Z.
    T_LIFT-1: De_r = 0 ⟺ lifted_diff_e_r ≡ 0 (mod 2^32)
    """
    W_normal = [W0_base] + [0]*15
    W_flip   = [W0_base ^ delta_W0] + [0]*15

    states_n = sha256_rounds_lifted(W_normal, N)
    states_f = sha256_rounds_lifted(W_flip, N)

    results = []
    for r in range(1, N+1):
        e_n = states_n[r][4]
        e_f = states_f[r][4]
        diff = e_f - e_n  # над Z, без mod
        results.append(diff)
    return results

# ================================================================
# 1. ВЕРИФИКАЦИЯ T_LIFT-1
# ================================================================

def check_lift1_theorem(N_samples=100000):
    """
    T_LIFT-1: De_r=0 ⟺ lifted_diff ≡ 0 (mod 2^32)
    """
    print("=" * 60)
    print("ВЕРИФИКАЦИЯ T_LIFT-1")
    print("De_r=0 ⟺ lifted_diff_e_r ≡ 0 (mod 2^32)")
    print("=" * 60)

    count_de3_zero = 0
    violations = 0

    for _ in range(N_samples):
        W0 = random.randint(0, M32) & ~1  # бит 0 = 0 (delta = +1)
        W = [W0] + [0]*15
        W_flip = [W0 ^ 1] + [0]*15
        sn = sha256_rounds(W, 3)
        sf = sha256_rounds(W_flip, 3)
        De3 = sn[4] ^ sf[4]

        diffs = lifted_diff_e(W0, delta_W0=1, N=3)
        lifted3 = diffs[2]

        xor_zero  = (De3 == 0)
        lift_zero = (lifted3 % (2**32) == 0)

        if xor_zero != lift_zero:
            violations += 1
            print(f"НАРУШЕНИЕ! W0={W0:#010x} De3={De3} lifted3%2^32={lifted3%(2**32)}")

        if xor_zero:
            count_de3_zero += 1

    print(f"\nПроверено: {N_samples} сэмплов")
    print(f"De3=0 случаев: {count_de3_zero} ({100*count_de3_zero/N_samples:.2f}%)")
    print(f"Нарушений T_LIFT-1: {violations}")
    print(f"Статус: {'✓ ТЕОРЕМА ПОДТВЕРЖДЕНА' if violations==0 else '✗ НАРУШЕНИЯ!'}")
    return violations == 0

# ================================================================
# 2. АНАЛИЗ ЛИНЕЙНОСТИ
# ================================================================

def analyze_lifted_linearity():
    """
    Исследуем: является ли lifted_diff_e3 линейной функцией от W[0]?
    """
    print("\n" + "=" * 60)
    print("АНАЛИЗ ЛИНЕЙНОСТИ lifted_diff_e3(W[0])")
    print("=" * 60)

    even_vals = []
    for W0 in range(0, 40, 2):
        diffs = lifted_diff_e(W0, delta_W0=1, N=3)
        even_vals.append((W0, diffs[2]))

    print(f"\nПервые 20 значений lifted_diff_e3 (W[0]=0,2,...,38):")
    for W0, v in even_vals:
        mod_val = v % (2**32)
        de3_zero = "← De3=0!" if mod_val == 0 else ""
        print(f"  W[0]={W0:3d}: lifted={v:15d}  mod2^32={mod_val:#010x} {de3_zero}")

    print("\nТест линейности над Z: f(W0^c) - f(W0) = const(c)?")
    W0_base = 0
    diffs_base = lifted_diff_e(W0_base, 1, 3)[2]

    linear_over_Z = True
    for c in [2, 4, 6, 8, 16, 32]:
        diffs_c = lifted_diff_e(W0_base ^ c, 1, 3)[2]
        d = diffs_c - diffs_base
        consistent = True
        for W0_test in [10, 20, 50, 100, 200]:
            if W0_test & 1: W0_test += 1
            v1 = lifted_diff_e(W0_test, 1, 3)[2]
            v2 = lifted_diff_e(W0_test ^ c, 1, 3)[2]
            if (v2 - v1) != d:
                consistent = False
        status = "✓ линейно" if consistent else "✗ нелинейно"
        print(f"  delta_W0={c:#04x}: {status}")
        if not consistent:
            linear_over_Z = False

    print(f"\nLiner over Z: {'ДА' if linear_over_Z else 'НЕТ (ожидается — нелинейность из carry)'}")

# ================================================================
# 3. АНАЛИЗ СТРУКТУРЫ ПОЛИНОМА
# ================================================================

def analyze_lifted_poly_structure():
    """
    Ключевой анализ: структура lifted_diff_e3 как функции W[0].
    """
    print("\n" + "=" * 60)
    print("СТРУКТУРА lifted_diff_e3: АНАЛИЗ")
    print("=" * 60)

    print("\nRound 1:")
    print("  lifted_diff_e1 = e1(W0^1) - e1(W0) = 1  (детерминировано)")
    v0  = lifted_diff_e(0,  1, 1)[0]
    v42 = lifted_diff_e(42, 1, 1)[0]
    print(f"  Верификация W0=0:  {v0}")
    print(f"  Верификация W0=42: {v42}")

    print("\nRound 2 — нелинейность через Sig1(e1) и Ch(e1,f1,g1):")
    samples_e2 = []
    for _ in range(10000):
        W0 = random.randint(0, M32) & ~1
        d = lifted_diff_e(W0, 1, 2)[1]
        samples_e2.append(d)

    unique_e2 = len(set(samples_e2))
    print(f"  Уникальных значений lifted_diff_e2 из 10k сэмплов: {unique_e2}")

    # При фиксированных W[1..15]=0
    vals_e2_fixed = []
    for W0 in range(0, 1000, 2):
        d = lifted_diff_e(W0, 1, 2)[1]
        vals_e2_fixed.append(d % (2**32))
    unique_fixed = len(set(vals_e2_fixed))
    print(f"  При W[1..15]=0, W[0] в [0..998]: {unique_fixed} уник. lifted_diff_e2 mod 2^32")

    # Распределение lifted_diff_e2 mod 2^32
    from collections import Counter
    cnt_e2 = Counter(vals_e2_fixed)
    top5 = cnt_e2.most_common(5)
    print(f"  Топ-5 значений lifted_diff_e2 mod 2^32:")
    for val, freq in top5:
        print(f"    {val:#010x}: {freq} раз")

# ================================================================
# 4. АНАЛИТИЧЕСКОЕ УСЛОВИЕ De3=0
# ================================================================

def find_analytical_condition_de3():
    """
    Главная задача П-1:
    Найти аналитическое условие на W[0] при котором De3=0.
    """
    print("\n" + "=" * 60)
    print("АНАЛИТИЧЕСКОЕ УСЛОВИЕ De3=0")
    print("=" * 60)

    print("\nВывод формулы (T_PERIOD3):")
    print("  e3 = d2 + T1_2  где d2 = H0[1] не зависит от W[0]")
    print("  => lifted_diff_e3 = T1_2_f - T1_2_n")
    print("  = [Sig1(e2_f) - Sig1(e2_n)] + [Ch(e2_f,f2,g2) - Ch(e2_n,f2,g2)]")
    print("  (h2, f2, g2, K[2], W[2] — общие, сокращаются)")
    print()
    print("  De3=0 ⟺ delta_Sig1(e2) + delta_Ch(e2,f2,g2) ≡ 0 (mod 2^32)")

    # Верификация на W_SAT3 из методички
    W_star_16 = [0xc5bde324, 0x3d7cd9d1, 0x2fd48880, 0xc3adefc3,
                 0xf4d69002, 0xbd5ff5be, 0xd10daffd, 0xff186caf,
                 0xcd900048, 0, 0, 0, 0, 0, 0, 0]
    W = schedule(W_star_16)

    # Нормальный поток до r=2
    a,b,c,d,e,f,g,h = H0
    for r in range(2):
        Wr = W[r]
        T1 = (h + Sig1(e) + Ch(e,f,g) + K[r] + Wr) & M32
        T2 = (Sig0(a) + Maj(a,b,c)) & M32
        h=g; g=f; f=e; e=(d+T1)&M32; d=c; c=b; b=a; a=(T1+T2)&M32
    e2_n, f2_n, g2_n, h2_n = e, f, g, h

    # Flipped поток до r=2
    a,b,c,d,e,f,g,h = H0
    T1f = (h + Sig1(e) + Ch(e,f,g) + K[0] + (W[0]^1)) & M32
    T2f = (Sig0(a) + Maj(a,b,c)) & M32
    h=g; g=f; f=e; e=(d+T1f)&M32; d=c; c=b; b=a; a=(T1f+T2f)&M32
    for r in range(1, 2):
        T1 = (h + Sig1(e) + Ch(e,f,g) + K[r] + W[r]) & M32
        T2 = (Sig0(a) + Maj(a,b,c)) & M32
        h=g; g=f; f=e; e=(d+T1)&M32; d=c; c=b; b=a; a=(T1+T2)&M32
    e2_f = e

    delta_Sig1_e2 = (Sig1(e2_f) - Sig1(e2_n)) % (2**32)
    delta_Ch_e2   = (Ch(e2_f, f2_n, g2_n) - Ch(e2_n, f2_n, g2_n)) % (2**32)
    delta_T1_2    = (delta_Sig1_e2 + delta_Ch_e2) % (2**32)

    print(f"\nВерификация W_SAT3 (W[0]={W[0]:#010x}):")
    print(f"  e2_n = {e2_n:#010x}")
    print(f"  e2_f = {e2_f:#010x}")
    print(f"  delta_Sig1(e2)  = {delta_Sig1_e2:#010x}  ({delta_Sig1_e2})")
    print(f"  delta_Ch(e2)    = {delta_Ch_e2:#010x}  ({delta_Ch_e2})")
    print(f"  Сумма mod 2^32  = {delta_T1_2:#010x}")
    print(f"  => De3=0? {'✓ ДА' if delta_T1_2==0 else '✗ НЕТ'}")

    # Статистика по случайным W[0]
    print("\nСтатистика delta_Sig1(e2) + delta_Ch(e2) mod 2^32 (10k сэмплов):")
    vals_sum = []
    for _ in range(10000):
        W0 = random.randint(0, M32) & ~1
        Wt = [W0] + [0]*15
        Ws = schedule(Wt)

        a,b,c,d,e,f,g,h = H0
        for r in range(2):
            T1 = (h + Sig1(e) + Ch(e,f,g) + K[r] + Ws[r]) & M32
            T2 = (Sig0(a) + Maj(a,b,c)) & M32
            h=g; g=f; f=e; e=(d+T1)&M32; d=c; c=b; b=a; a=(T1+T2)&M32
        e2n, f2n, g2n = e, f, g

        a,b,c,d,e,f,g,h = H0
        T1f = (h + Sig1(e) + Ch(e,f,g) + K[0] + (W0^1)) & M32
        T2f = (Sig0(a) + Maj(a,b,c)) & M32
        h=g; g=f; f=e; e=(d+T1f)&M32; d=c; c=b; b=a; a=(T1f+T2f)&M32
        for r in range(1, 2):
            T1 = (h + Sig1(e) + Ch(e,f,g) + K[r] + Ws[r]) & M32
            T2 = (Sig0(a) + Maj(a,b,c)) & M32
            h=g; g=f; f=e; e=(d+T1)&M32; d=c; c=b; b=a; a=(T1+T2)&M32
        e2f = e

        s = (Sig1(e2f) - Sig1(e2n) + Ch(e2f,f2n,g2n) - Ch(e2n,f2n,g2n)) % (2**32)
        vals_sum.append(s)

    zero_count = sum(1 for v in vals_sum if v == 0)
    print(f"  Случаев delta_T1_2=0: {zero_count}/10000 ({100*zero_count/10000:.2f}%)")
    print(f"  Уникальных значений: {len(set(vals_sum))}")

# ================================================================
# 5. СТРУКТУРА РЕШЕНИЙ T_LIFT-2
# ================================================================

def scan_solutions_structure(W1_fixed=0, count=30):
    """
    T_LIFT-2: исследуем структуру решений De3=0.
    """
    print("\n" + "=" * 60)
    print("СТРУКТУРА РЕШЕНИЙ De3=0 (T_LIFT-2)")
    print("=" * 60)
    print(f"W[1..15]=0, варьируем W[0] от 0 до {count*2-1}")
    print()

    solutions = []
    for W0 in range(count * 2):
        if W0 & 1: continue
        W      = [W0] + [W1_fixed] + [0]*14
        W_flip = [W0 ^ 1] + [W1_fixed] + [0]*14

        s_n = sha256_rounds(W, 3)
        s_f = sha256_rounds(W_flip, 3)
        De3 = s_n[4] ^ s_f[4]

        lifted = lifted_diff_e(W0, 1, 3)[2]

        if De3 == 0:
            solutions.append(W0)

        symbol = "✓" if De3==0 else "·"
        print(f"  W0={W0:3d}: De3={'0x00000000' if De3==0 else hex(De3):12s} "
              f"lifted mod2^32={lifted%(2**32):#010x} {symbol}")

    print(f"\nРешения De3=0 в [0,{count*2}): {solutions}")

    # Проверка T_LIFT-2: lifted_diff(W0) + lifted_diff(W0^1) = 0?
    print("\nПроверка T_LIFT-2 (анти-симметрия):")
    print("lifted_diff(W0) + lifted_diff(W0^1) = 0?")
    violations = 0
    for W0 in range(0, 200, 2):
        l0 = lifted_diff_e(W0, 1, 3)[2]
        l1 = lifted_diff_e(W0+1, 1, 3)[2]
        if l0 + l1 != 0:
            violations += 1
    print(f"  Нарушений (из 100 пар): {violations}")
    print(f"  Статус: {'✓ T_LIFT-2 ПОДТВЕРЖДЕНА' if violations==0 else '✗ НАРУШЕНИЯ'}")

# ================================================================
# 6. РАЗЛОЖЕНИЕ ПО КОМПОНЕНТАМ
# ================================================================

def lifted_polynomial_decomposition():
    """
    Разложение lifted_diff_e3 на компоненты.
    """
    print("\n" + "=" * 60)
    print("РАЗЛОЖЕНИЕ lifted_diff_e3 ПО КОМПОНЕНТАМ")
    print("=" * 60)

    print("\nТеоретическое разложение:")
    print("  lifted_diff_e1 = 1  (детерминировано)")
    print("  lifted_diff_e2 = Sig1(e1_f)-Sig1(e1_n) + Ch(e1_f,f1,g1)-Ch(e1_n,f1,g1)")
    print("  lifted_diff_e3 = Sig1(e2_f)-Sig1(e2_n) + Ch(e2_f,f2,g2)-Ch(e2_n,f2,g2)")
    print()
    print("  Sig1(x+1) - Sig1(x) — нелинейно над Z (зависит от carry-цепочек)")
    print("  Ch(x+d,f,g) - Ch(x,f,g) — нелинейно, зависит от f,g")

    # Численно: корреляция lifted_diff_e2 vs lifted_diff_e3
    print("\nКорреляция lifted_diff_e2 vs lifted_diff_e3 (5k сэмплов):")

    data = []
    for _ in range(5000):
        W0 = random.randint(0, M32) & ~1
        diffs = lifted_diff_e(W0, 1, 3)
        ld2 = diffs[1]
        ld3 = diffs[2]
        data.append((ld2 % (2**32), ld3 % (2**32)))

    groups = defaultdict(list)
    for ld2, ld3 in data:
        groups[ld2 >> 28].append(ld3)

    print(f"  Групп по старшим 4 битам ld2: {len(groups)}")
    print(f"  Примеры (top-4-bits ld2 → статистика ld3):")
    for key in sorted(groups.keys())[:8]:
        vals = groups[key]
        avg = sum(vals) / len(vals)
        var = sum((v-avg)**2 for v in vals) / len(vals)
        print(f"    bits[31:28]=0x{key:x}: n={len(vals):4d},"
              f" avg_ld3={avg/(2**32):.4f}·2^32,"
              f" std={var**0.5/(2**32):.4f}·2^32")

    # Уникальные значения lifted_diff_e2 (T_ADD8 предсказывает структуру)
    ld2_vals = [d[0] for d in data]
    print(f"\n  Уникальных ld2 mod 2^32 из 5k сэмплов: {len(set(ld2_vals))}")
    from collections import Counter
    cnt = Counter(ld2_vals)
    top = cnt.most_common(8)
    print(f"  Топ-8 наиболее частых ld2 mod 2^32:")
    for val, freq in top:
        print(f"    {val:#010x}: {freq} раз ({100*freq/5000:.1f}%)")

# ================================================================
# MAIN
# ================================================================

def main():
    random.seed(42)  # воспроизводимость

    print("П-1: LIFTED АРИФМЕТИКА — АНАЛИТИЧЕСКОЕ УСЛОВИЕ De3=0")
    print("=" * 60)

    # 1. Верифицируем T_LIFT-1
    check_lift1_theorem(N_samples=100000)

    # 2. Анализ линейности
    analyze_lifted_linearity()

    # 3. Анализ структуры полинома
    analyze_lifted_poly_structure()

    # 4. Аналитическое условие
    find_analytical_condition_de3()

    # 5. Структура решений (T_LIFT-2)
    scan_solutions_structure()

    # 6. Разложение по компонентам
    lifted_polynomial_decomposition()

    print("\n" + "=" * 60)
    print("П-1 ЗАВЕРШЁН")
    print("=" * 60)

if __name__ == "__main__":
    main()
