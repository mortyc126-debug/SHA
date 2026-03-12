"""
p27_verify_theorems.py — Верификация теорем v27 (Раздел 13)

V2: T_Da_GENERAL — Da_k = De_k - Da_{k-4} + ΔT2_{k-1}  (точная идентичность)
V1: T_P23_PHASE1 — измерение реальной P(3-раундовый XOR-след)
    с фильтрацией δe1 = 2^j и без (обе версии)
V3: T_JOINT_ZERO — аналитическое подтверждение
"""

import random

# ─── базовые операции SHA-256 ────────────────────────────────────────────────
MASK = 0xFFFFFFFF

def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK
def Sig0(x):  return rotr(x,2)  ^ rotr(x,13) ^ rotr(x,22)
def Sig1(x):  return rotr(x,6)  ^ rotr(x,11) ^ rotr(x,25)
def Ch(e,f,g):  return ((e & f) ^ (~e & g)) & MASK
def Maj(a,b,c): return ((a & b) ^ (a & c) ^ (b & c)) & MASK

K = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
    0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
    0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
]

IV = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
]

def sig0(x):  return rotr(x,7)  ^ rotr(x,18) ^ (x >> 3)
def sig1(x):  return rotr(x,17) ^ rotr(x,19) ^ (x >> 10)

def schedule(W16):
    W = list(W16)
    for i in range(16, 64):
        W.append((sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16]) & MASK)
    return W

def sha_rounds(W16, R, iv=None):
    """R раундов. Возвращает состояния [round0..roundR]."""
    W = schedule(W16)
    a,b,c,d,e,f,g,h = iv if iv else IV
    states = [(a,b,c,d,e,f,g,h)]
    for r in range(R):
        T1 = (h + Sig1(e) + Ch(e,f,g) + K[r] + W[r]) & MASK
        T2 = (Sig0(a) + Maj(a,b,c)) & MASK
        h=g; g=f; f=e
        e=(d+T1) & MASK
        d=c; c=b; b=a
        a=(T1+T2) & MASK
        states.append((a,b,c,d,e,f,g,h))
    return states

# ─── V2: T_Da_GENERAL (точная идентичность) ──────────────────────────────────
def verify_T_Da_GENERAL(N=5000):
    """Da_k = De_k - Da_{k-4} + ΔT2_{k-1}  (mod 2^32) — должна быть ТОЧНОЙ."""
    print("=" * 62)
    print("V2: T_Da_GENERAL   Da_k = De_k − Da_{k−4} + ΔT2_{k−1}")
    print("=" * 62)
    failures = 0
    for _ in range(N):
        W1 = [random.randint(0, MASK) for _ in range(16)]
        dw  = random.randint(1, MASK)
        W2  = list(W1); W2[0] = (W1[0] + dw) & MASK
        st1 = sha_rounds(W1, 17)
        st2 = sha_rounds(W2, 17)
        for k in range(4, 17):
            Da_k   = (st2[k][0] - st1[k][0]) & MASK
            De_k   = (st2[k][4] - st1[k][4]) & MASK
            Da_km4 = (st2[k-4][0] - st1[k-4][0]) & MASK
            T2_1 = (Sig0(st1[k-1][0]) + Maj(st1[k-1][0], st1[k-1][1], st1[k-1][2])) & MASK
            T2_2 = (Sig0(st2[k-1][0]) + Maj(st2[k-1][0], st2[k-1][1], st2[k-1][2])) & MASK
            DT2   = (T2_2 - T2_1) & MASK
            rhs   = (De_k - Da_km4 + DT2) & MASK
            if Da_k != rhs:
                failures += 1
    total = N * 13
    print(f"  Проверено: {total:,} (k=4..16, N={N})")
    print(f"  Нарушений: {failures}")
    if failures == 0:
        print(f"  ✓ ТОЧНАЯ ИДЕНТИЧНОСТЬ — 0 нарушений из {total:,}")
    else:
        print(f"  ✗ ОШИБКА: {failures}/{total} нарушений")
    return failures == 0

# ─── V1: T_P23_PHASE1 — XOR-дифференциал через 3 раунда ─────────────────────
def verify_T_P23_PHASE1_exact(N=20000, j=15):
    """
    Версия с ФИЛЬТРАЦИЕЙ: только пары где δe1 = 2^j (XOR, точно).
    Выбираем δW1 = Sig1(2^j) для компенсации.
    Меряем P(δe2 = 0) и P(HW(δe3) = 0..1).
    """
    bit = 1 << j
    Sig1_bit = Sig1(bit)

    # Сначала посмотрим: какова P(δe1 = 2^j) при δW0 = 2^j
    count_de1_exact = 0
    for _ in range(N):
        W1 = [random.randint(0, MASK) for _ in range(16)]
        W2 = list(W1); W2[0] ^= bit
        st1 = sha_rounds(W1, 1); st2 = sha_rounds(W2, 1)
        if (st1[1][4] ^ st2[1][4]) == bit:
            count_de1_exact += 1
    P_de1_exact = count_de1_exact / N
    print(f"  P(δe1=2^{j} | δW0=2^{j}) = {P_de1_exact:.4f}  (теория: ~0.5 из-за переносов)")

    # Основной эксперимент с фильтрацией δe1 = 2^j
    hits_de2_zero = 0
    hits_de3_hw0  = 0
    hits_de3_hw1  = 0
    total_filtered = 0
    de3_vals = {}

    attempts = 0
    while total_filtered < 5000:
        W1 = [random.randint(0, MASK) for _ in range(16)]
        W2 = list(W1); W2[0] ^= bit; W2[1] ^= Sig1_bit
        st1 = sha_rounds(W1, 1); st2 = sha_rounds(W2, 1)
        de1 = st1[1][4] ^ st2[1][4]
        attempts += 1
        if de1 != bit:
            continue  # фильтрация: только чистые δe1 = 2^j
        total_filtered += 1

        st1f = sha_rounds(W1, 3)
        st2f = sha_rounds(W2, 3)
        de2 = st1f[2][4] ^ st2f[2][4]
        de3 = st1f[3][4] ^ st2f[3][4]
        hw3 = bin(de3).count('1')

        if de2 == 0:   hits_de2_zero += 1
        if hw3 == 0:   hits_de3_hw0  += 1
        if hw3 <= 1:   hits_de3_hw1  += 1
        de3_vals[de3] = de3_vals.get(de3, 0) + 1

    p2 = hits_de2_zero / total_filtered
    p3_hw0 = hits_de3_hw0 / total_filtered
    p3_hw1 = hits_de3_hw1 / total_filtered
    efficiency = total_filtered / attempts

    print(f"  Фильтрация: {total_filtered} пар из {attempts} попыток  (P_filter={efficiency:.3f})")
    print(f"  P(δe2=0)       = {p2:.4f}  (теория GF(2): 0.5)")
    print(f"  P(δe3=0)       = {p3_hw0:.4f}  (теория GF(2): 0.25)")
    print(f"  P(HW(δe3)≤1)   = {p3_hw1:.4f}")
    top5 = sorted(de3_vals.items(), key=lambda x: -x[1])[:5]
    print(f"  Топ-5 δe3:")
    for v, c in top5:
        print(f"    0x{v:08x}  HW={bin(v).count('1')}  P={c/total_filtered:.4f}")

    return p2, p3_hw0, total_filtered

def verify_T_P23_PHASE1_additive(N=5000, j=15):
    """
    Версия с АДДИТИВНЫМ дифференциалом δW0 = +2^j (не XOR).
    В аддитивной модели: δe1 = (C+W0+2^j)-(C+W0) = 2^j (если нет переполнения).
    """
    bit = 1 << j

    # P(δe1_add = 2^j | δW0_add = 2^j)
    count_exact = 0; count_nonzero = 0
    for _ in range(N):
        W1 = [random.randint(0, MASK) for _ in range(16)]
        W2 = list(W1); W2[0] = (W1[0] + bit) & MASK
        st1 = sha_rounds(W1, 1); st2 = sha_rounds(W2, 1)
        de1_add = (st2[1][4] - st1[1][4]) & MASK
        if de1_add == bit: count_exact += 1
        if de1_add != 0:   count_nonzero += 1
    P_exact_add = count_exact / N
    print(f"  [Аддитивный] P(δe1_add=2^{j}) = {P_exact_add:.4f}  (ожидаем: ~1.0 из-за линейности)")

    # Аддитивные эксперименты
    hits2 = 0; hits3_hw0 = 0; hits3_hw1 = 0
    total = 0; attempts = 0
    while total < min(N, 3000):
        W1 = [random.randint(0, MASK) for _ in range(16)]
        W2 = list(W1); W2[0] = (W1[0] + bit) & MASK
        # Оптимальный дельта W1 для аддитивной модели — закрытая формула сложна
        # Просто берём δW1=0 чтобы измерить "естественный" переход
        st1 = sha_rounds(W1, 3); st2 = sha_rounds(W2, 3)
        de1 = (st2[1][4] - st1[1][4]) & MASK
        attempts += 1
        if de1 != bit:
            continue
        total += 1
        de2 = (st2[2][4] - st1[2][4]) & MASK
        de3 = (st2[3][4] - st1[3][4]) & MASK
        hw3 = bin(de3 if de3 <= MASK//2 else MASK+1-de3).count('1')
        if de2 == 0:   hits2    += 1
        if de3 == 0:   hits3_hw0 += 1
        if hw3 <= 1:   hits3_hw1 += 1
    if total > 0:
        print(f"  [Аддитивный δW1=0] P(δe2=0)={hits2/total:.4f}  P(δe3=0)={hits3_hw0/total:.4f}  P(HW≤1)={hits3_hw1/total:.4f}  N={total}")
    return P_exact_add

# ─── АНАЛИЗ ПЕРЕНОСОВ ────────────────────────────────────────────────────────
def analyze_carry_structure(j=15, N=50000):
    """
    Анализируем: когда XOR δW0 = 2^j приводит к δe1 = 2^j (vs multi-bit)?
    Результат: как часто перенос портит δe1.
    """
    bit = 1 << j
    good = 0; bad_hw = {}
    for _ in range(N):
        W1 = [random.randint(0, MASK) for _ in range(16)]
        W2 = list(W1); W2[0] ^= bit
        st1 = sha_rounds(W1, 1); st2 = sha_rounds(W2, 1)
        de1 = st1[1][4] ^ st2[1][4]
        if de1 == bit:
            good += 1
        else:
            hw = bin(de1).count('1')
            bad_hw[hw] = bad_hw.get(hw, 0) + 1

    print(f"  Carry-анализ j={j}: P(δe1=2^j)={good/N:.4f}")
    print(f"  Распределение HW(δe1) при δe1≠2^j:")
    for hw in sorted(bad_hw.keys())[:6]:
        print(f"    HW={hw}: {bad_hw[hw]/N:.4f}")

# ─── MAIN ─────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    random.seed(42)
    print()
    print("╔══════════════════════════════════════════════════════════════╗")
    print("║  p27_verify_theorems.py — SHA-256 Верификация v27           ║")
    print("╚══════════════════════════════════════════════════════════════╝")
    print()

    # V2: T_Da_GENERAL — должна быть ТОЧНОЙ
    ok_v2 = verify_T_Da_GENERAL(N=5000)
    print()

    # Carry-структура для понимания XOR-дифференциала
    print("=" * 62)
    print("Анализ carry-структуры (почему XOR-след сложнее GF(2)-теории)")
    print("=" * 62)
    analyze_carry_structure(j=15, N=50000)
    analyze_carry_structure(j=0,  N=50000)
    print()

    # V1: T_P23_PHASE1 с фильтрацией
    print("=" * 62)
    print("V1: T_P23_PHASE1 — XOR-дифференциал (с фильтрацией δe1=2^j)")
    print("=" * 62)
    for j in [0, 7, 15, 24]:
        print(f"\n  --- j={j} ---")
        p2, p3, n = verify_T_P23_PHASE1_exact(N=40000, j=j)
        print(f"  Итог j={j}: P(δe2=0)={p2:.4f}  P(δe3=0)={p3:.4f}  n={n}")
    print()

    # V1b: аддитивная версия
    print("=" * 62)
    print("V1b: T_P23_PHASE1 — аддитивный дифференциал (δW0=+2^j)")
    print("=" * 62)
    for j in [15, 0]:
        print(f"\n  --- j={j} ---")
        verify_T_P23_PHASE1_additive(N=5000, j=j)
    print()

    # Итоговый отчёт
    print("=" * 62)
    print("ИТОГОВЫЙ ОТЧЁТ v27:")
    print(f"  V2 T_Da_GENERAL : {'✓ ТОЧНАЯ ИДЕНТИЧНОСТЬ' if ok_v2 else '✗ ОШИБКА'}")
    print(f"  V1 T_P23_PHASE1 : см. таблицу выше")
    print(f"  V3 T_JOINT_ZERO : P(De=Da=0) ~2^{{-64}} — прямая верификация невозможна")
    print(f"                    аналитически верна (следствие T_Da_GENERAL)")
    print("=" * 62)
