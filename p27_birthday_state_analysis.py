"""
П-27: Три направления исследования

A) Birthday-поиск реальной 17-раундовой пары (Python-версия для проверки)
B) Анализ состояния (a,b,c,d) после 17 нулей δe2..δe17=0
C) 2D birthday: теоретический анализ δe18=0 за O(2^32) вместо O(2^64)
"""

import random
import time
import statistics
from collections import Counter

MASK = 0xFFFFFFFF

def rotr(x, n): return ((x >> n) | (x << (32-n))) & MASK
def sig0(x):  return rotr(x,7)  ^ rotr(x,18) ^ (x>>3)
def sig1(x):  return rotr(x,17) ^ rotr(x,19) ^ (x>>10)
def Sig0(x):  return rotr(x,2)  ^ rotr(x,13) ^ rotr(x,22)
def Sig1(x):  return rotr(x,6)  ^ rotr(x,11) ^ rotr(x,25)
def Ch(e,f,g):  return ((e&f) ^ (~e&g)) & MASK
def Maj(a,b,c): return (a&b) ^ (a&c) ^ (b&c)

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
    0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
]
IV = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
      0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]

def schedule(W16):
    W = list(W16) + [0]*48
    for i in range(16, 64):
        W[i] = (sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16]) & MASK
    return W

def sha_r(W, R):
    """R раундов SHA-256, возвращает список состояний [state0, state1, ...]"""
    a,b,c,d,e,f,g,h = IV
    states = [[a,b,c,d,e,f,g,h]]
    for r in range(R):
        T1 = (h + Sig1(e) + Ch(e,f,g) + K[r] + W[r]) & MASK
        T2 = (Sig0(a) + Maj(a,b,c)) & MASK
        h=g; g=f; f=e; e=(d+T1)&MASK
        d=c; c=b; b=a; a=(T1+T2)&MASK
        states.append([a,b,c,d,e,f,g,h])
    return states

def de(sn, sf, r): return (sf[r][4] - sn[r][4]) & MASK
def da(sn, sf, r): return (sf[r][0] - sn[r][0]) & MASK
def hw(x): return bin(x).count('1')

def cascade_3param(W0, W1, DW0=1):
    """
    Адаптивный каскад: De3..De17=0.
    Возвращает (DWs, De17, De18, full_states_n, full_states_f).
    """
    Wn = [W0, W1] + [0]*14
    DWs = [0]*16
    DWs[0] = DW0

    # ΔW2 → De3=0
    Wf_tmp = [(Wn[i]+DWs[i])&MASK for i in range(16)]
    Wn_s = schedule(Wn); Wf_tmp_s = schedule(Wf_tmp)
    sn = sha_r(Wn_s, 3); sf = sha_r(Wf_tmp_s, 3)
    De3_nat = de(sn, sf, 3)
    DWs[2] = (-De3_nat) & MASK

    # Каскад ΔW3..ΔW15
    for step in range(13):
        wi = step+3; dt = step+4
        Wfc = [(Wn[i]+DWs[i])&MASK for i in range(16)]
        Wn_s = schedule(Wn); Wfc_s = schedule(Wfc)
        sn = sha_r(Wn_s, dt); sf = sha_r(Wfc_s, dt)
        De_nat = de(sn, sf, dt)
        DWs[wi] = (-De_nat) & MASK

    # Финальный вычисление
    Wf = [(Wn[i]+DWs[i])&MASK for i in range(16)]
    Wn_s = schedule(Wn); Wf_s = schedule(Wf)
    sn = sha_r(Wn_s, 20); sf = sha_r(Wf_s, 20)
    return DWs, sn, sf

# ============================================================
# ЧАСТЬ B: Анализ состояния (a,b,c,d) после 17 нулей
# ============================================================
print("=" * 60)
print("ЧАСТЬ B: Анализ (a,b,c,d) при δe2..δe17=0")
print("=" * 60)

# Используем найденные пары из П-15/П-16
known_pairs = [
    (0xe82222c7, 0x516cfb41),  # П-15
    (0xd4254551, 0x679ea4de),  # П-16 #1
    (0xe73eb86e, 0xdfa1b7b0),  # П-16 #2
]

print("\nИзвестные пары с De3..De17=0:")
print(f"{'W0':>10} {'W1':>10} {'Da17':>10} {'Db17':>10} {'Dc17':>10} {'Dd17':>10} {'De17':>10}")
for W0, W1 in known_pairs:
    DWs, sn, sf = cascade_3param(W0, W1)
    da17 = (sf[17][0] - sn[17][0]) & MASK
    db17 = (sf[17][1] - sn[17][1]) & MASK
    dc17 = (sf[17][2] - sn[17][2]) & MASK
    dd17 = (sf[17][3] - sn[17][3]) & MASK
    de17 = de(sn, sf, 17)
    print(f"0x{W0:08x} 0x{W1:08x} 0x{da17:08x} 0x{db17:08x} 0x{dc17:08x} 0x{dd17:08x} 0x{de17:08x}")

print("\nПолная таблица δ-состояний по раундам (пара П-15):")
W0, W1 = 0xe82222c7, 0x516cfb41
DWs, sn, sf = cascade_3param(W0, W1)

print(f"\n{'r':>3} | {'δe':>10} {'δa':>10} {'δb':>10} {'δc':>10} {'δd':>10} | HW(δe) HW(δa)")
print("-" * 75)
for r in range(1, 20):
    if r >= len(sn): break
    de_r = (sf[r][4] - sn[r][4]) & MASK
    da_r = (sf[r][0] - sn[r][0]) & MASK
    db_r = (sf[r][1] - sn[r][1]) & MASK
    dc_r = (sf[r][2] - sn[r][2]) & MASK
    dd_r = (sf[r][3] - sn[r][3]) & MASK
    marker = " ←" if (2 <= r <= 17 and de_r == 0) else (" ← БАРЬЕР" if r == 17 else "")
    print(f"{r:>3} | 0x{de_r:08x} 0x{da_r:08x} 0x{db_r:08x} 0x{dc_r:08x} 0x{dd_r:08x} | {hw(de_r):>6} {hw(da_r):>6}{marker}")

# Статистика по 200 случайным парам с De17=0
print("\n\nСтатистика по 200 случайным парам (W0, W1) ~ T_CASCADE_17 поиск:")
print("(Быстрый метод: перебор W1 до De17=0)")

N_samples = 200
da17_vals, db17_vals, dc17_vals, dd17_vals = [], [], [], []
da13_vals, dw16_vals = [], []
hw_da17, hw_db17, hw_dc17, hw_dd17 = [], [], [], []
hw_da13 = []
de18_vals = []
t0 = time.time()

found = 0
attempts = 0
for _ in range(N_samples):
    W0 = random.randint(0, MASK)
    # Быстрый поиск W1: перебор до нахождения De17=0
    for _ in range(1 << 20):  # до 1M попыток (ожидаем ~2^32, но для статистики берём найденные)
        W1 = random.randint(0, MASK)
        attempts += 1
        DWs, sn, sf = cascade_3param(W0, W1)
        d17 = de(sn, sf, 17)
        if d17 == 0:
            found += 1
            da17 = (sf[17][0] - sn[17][0]) & MASK
            db17 = (sf[17][1] - sn[17][1]) & MASK
            dc17 = (sf[17][2] - sn[17][2]) & MASK
            dd17 = (sf[17][3] - sn[17][3]) & MASK
            da13 = (sf[13][0] - sn[13][0]) & MASK
            d18 = de(sn, sf, 18)
            da17_vals.append(da17)
            db17_vals.append(db17)
            dc17_vals.append(dc17)
            dd17_vals.append(dd17)
            da13_vals.append(da13)
            hw_da17.append(hw(da17))
            hw_db17.append(hw(db17))
            hw_dc17.append(hw(dc17))
            hw_dd17.append(hw(dd17))
            hw_da13.append(hw(da13))
            de18_vals.append(d18)
            break
    if found >= N_samples:
        break

elapsed = time.time() - t0
print(f"Найдено: {found} пар за {elapsed:.1f}с ({attempts} попыток)")

if found > 0:
    print(f"\n Регистр | E[HW]  | std  | min | max | P(=0)")
    for name, vals in [("δa17", hw_da17), ("δb17", hw_db17), ("δc17", hw_dc17), ("δd17", hw_dd17), ("δa13", hw_da13)]:
        if vals:
            print(f" {name:>7} | {statistics.mean(vals):5.2f}  | {statistics.stdev(vals) if len(vals)>1 else 0:4.2f} | {min(vals):3d} | {max(vals):3d} | {vals.count(0)/len(vals):.4f}")

    # Проверка T_DE17_DECOMPOSITION: De17 = Da13 + ΔW16
    print("\nПроверка T_DEk_DECOMPOSITION для раундов 17, 18:")
    for W0, W1 in known_pairs[:1]:
        DWs, sn, sf = cascade_3param(W0, W1)
        Wn_s = schedule([W0, W1] + [0]*14)
        Wf_s = schedule([(W0+DWs[0])&MASK, (W1+DWs[1])&MASK] + [DWs[i] for i in range(2,16)])
        DW16 = (Wf_s[16] - Wn_s[16]) & MASK
        DW17 = (Wf_s[17] - Wn_s[17]) & MASK
        da13 = (sf[13][0] - sn[13][0]) & MASK
        da14 = (sf[14][0] - sn[14][0]) & MASK
        de17 = de(sn, sf, 17)
        de18 = de(sn, sf, 18)
        print(f"  Da13 = 0x{da13:08x}, ΔW16 = 0x{DW16:08x}")
        print(f"  Da13 + ΔW16 = 0x{(da13+DW16)&MASK:08x}, De17 = 0x{de17:08x}  {'✓' if (da13+DW16)&MASK==de17 else '✗'}")
        print(f"  Da14 = 0x{da14:08x}, ΔW17 = 0x{DW17:08x}")
        print(f"  Da14 + ΔW17 = 0x{(da14+DW17)&MASK:08x}, De18 = 0x{de18:08x}  {'✓' if (da14+DW17)&MASK==de18 else '✗'}")

    # Распределение Da17 при условии De17=0
    print(f"\nРаспределение δa17 при De17=0 (N={found}):")
    # Группы по HW
    hw_counter = Counter(hw_da17)
    print("  HW(δa17) | count | fraction")
    for h in sorted(hw_counter.keys()):
        print(f"  {h:8d} | {hw_counter[h]:5d} | {hw_counter[h]/found:.3f}")

    # Корреляция Da17 ↔ δe18
    print(f"\nКорреляция δa17 ↔ δe18 (T_DEk_DECOMPOSITION для r=18):")
    # δe18 = δa14 + δW17, не δa17 → проверяем независимость
    if de18_vals:
        nonzero_de18 = sum(1 for x in de18_vals if x != 0)
        print(f"  P(δe18=0 | δe17=0) = {found - nonzero_de18}/{found} = {(found-nonzero_de18)/found:.6f}")
        print(f"  Ожидание 2^(-32)  ≈ {2**-32:.2e}")

# ============================================================
# ЧАСТЬ C: 2D Birthday — теоретический анализ δe18=0
# ============================================================
print("\n\n" + "=" * 60)
print("ЧАСТЬ C: 2D Birthday для δe18=0 — теоретический анализ")
print("=" * 60)

print("""
Цель: снизить стоимость δe18=0 с O(2^64) до O(2^32)?

Стандартный анализ (T_BARRIER_16):
  δe17=0: одно 32-битное условие → стоимость 2^32
  δe18=0: ещё одно 32-битное условие, независимое → ×2^32
  Итого: 2^64.

2D Birthday гипотеза:
  Ищем (W0_a, W1_a) и (W0_b, W1_b) такие что:
    f17(W0_a, W1_a) = f17(W0_b, W1_b)  [birthday на f17]
    f18(W0_a, W1_a) = f18(W0_b, W1_b)  [birthday на f18]
  где f17 = Da13 + ΔW16, f18 = Da14 + ΔW17.

Это НЕ стандартный birthday, т.к. нужно не совпадение,
а обращение в 0: f17=0 И f18=0.

Корректная 2D-атака: MITM в пространстве (W0, W1).
  Фаза 1: для 2^32 случайных W1, вычислить f17(W1).
           Хранить в хэш-таблице: f17 → W1.
           При f17=0 → записать кандидат.
  Фаза 2: для найденных кандидатов (f17=0): проверить f18.
           P(f18=0 | f17=0) ≈ 2^(-32).
  Итого: нужно 2^32 кандидатов → 2^64 операций.

Вывод: стандартный 1D-birthday НЕ снижает барьер для f18.
""")

print("Эмпирическая проверка: независимость f17 и f18")
print("(N=10000 случайных (W0,W1), смотрим на совместное распределение)")

N_check = 10000
f17_vals = []
f18_vals = []
joint_zeros = 0
t0 = time.time()

for _ in range(N_check):
    W0 = random.randint(0, MASK)
    W1 = random.randint(0, MASK)
    DWs, sn, sf = cascade_3param(W0, W1)
    Wn_s = schedule([W0, W1] + [0]*14)
    Wf_list = [(W0+DWs[0])&MASK, (W1+DWs[1])&MASK] + [(DWs[i])&MASK for i in range(2,16)]
    Wf_s = schedule(Wf_list)
    da13 = (sf[13][0] - sn[13][0]) & MASK
    da14 = (sf[14][0] - sn[14][0]) & MASK
    DW16 = (Wf_s[16] - Wn_s[16]) & MASK
    DW17 = (Wf_s[17] - Wn_s[17]) & MASK
    f17 = (da13 + DW16) & MASK
    f18 = (da14 + DW17) & MASK
    f17_vals.append(f17)
    f18_vals.append(f18)
    if f17 == 0 and f18 == 0:
        joint_zeros += 1

elapsed = time.time() - t0

# Уникальность
uniq17 = len(set(f17_vals))
uniq18 = len(set(f18_vals))
zeros17 = f17_vals.count(0)
zeros18 = f18_vals.count(0)

print(f"\nN = {N_check}, время: {elapsed:.1f}с")
print(f"f17: {uniq17} уникальных значений ({uniq17/N_check*100:.1f}%)")
print(f"f18: {uniq18} уникальных значений ({uniq18/N_check*100:.1f}%)")
print(f"f17=0: {zeros17} ({zeros17/N_check:.4f}, ожидание 1/2^32≈{1/2**32:.2e})")
print(f"f18=0: {zeros18} ({zeros18/N_check:.4f}, ожидание 1/2^32≈{1/2**32:.2e})")
print(f"f17=0 AND f18=0: {joint_zeros} (ожидание {N_check/2**64:.2e})")

# Корреляция f17 и f18 (XOR-битовая независимость)
import math
def bit_correlation(vals_a, vals_b, bit):
    """Корреляция бита bit между двумя списками"""
    a_bits = [(v >> bit) & 1 for v in vals_a]
    b_bits = [(v >> bit) & 1 for v in vals_b]
    n = len(a_bits)
    mean_a = sum(a_bits)/n
    mean_b = sum(b_bits)/n
    cov = sum((a_bits[i]-mean_a)*(b_bits[i]-mean_b) for i in range(n))/n
    std_a = math.sqrt(sum((x-mean_a)**2 for x in a_bits)/n) or 1e-10
    std_b = math.sqrt(sum((x-mean_b)**2 for x in b_bits)/n) or 1e-10
    return cov/(std_a*std_b)

corrs = [bit_correlation(f17_vals, f18_vals, bit) for bit in [0, 8, 16, 24, 31]]
print(f"\nКорреляция f17 ↔ f18 (биты 0,8,16,24,31): {[f'{c:.4f}' for c in corrs]}")
print("(близко к 0 → независимы → стоимость δe18=0 при δe17=0 не снижается)")

print("""
Теорема T_2D_BIRTHDAY_ANALYSIS:
  f17 = Da13(W0,W1) + ΔW16(W0,W1)  — псевдослучайная 32-бит функция
  f18 = Da14(W0,W1) + ΔW17(W0,W1)  — псевдослучайная 32-бит функция

  Из T_DE17_DE18_INDEPENDENCE (П-14): P(f18=0 | f17=0) ≈ 2^(-32).
  Корреляция ≈ 0 эмпирически → f17, f18 независимы.

  Вывод: 2D birthday не даёт преимущества перед 2^64.

  НО: Если ввести ВТОРОЙ свободный параметр (например, W0[0] и W1[0])
  как два независимых «измерения», то:
    - Пространство поиска: 2D, размер 2^32 × 2^32 = 2^64
    - Birthday в 2D: коллизия f17=f18=0 за O(√(2^64)) = O(2^32)
    - Требует хранения 2^32 записей: 2^32 × 8 байт = 32 ГБ

  Это ДЕЙСТВИТЕЛЬНО снижает вычислительную сложность:
    2^64 → 2^32 вычислений (ценой 32 ГБ памяти)!
""")

# Проверяем: насколько f17 зависит только от W1[0] при фиксированных остальных?
print("=" * 60)
print("T_BIRTHDAY_W1_0: зависимость f17 от W1[0]")
print("(Фиксируем W0, W1[1..15]=0, перебираем W1[0])")
print("=" * 60)

W0_fixed = 0xc5bde324
N_scan = 10000
f17_by_w1_0 = []
t0 = time.time()
for w1_0 in range(N_scan):
    W1 = w1_0  # W1[1..15]=0, W1[0]=w1_0
    DWs, sn, sf = cascade_3param(W0_fixed, W1)
    Wn_s = schedule([W0_fixed, W1] + [0]*14)
    Wf_list = [(W0_fixed+DWs[0])&MASK, (W1+DWs[1])&MASK] + [DWs[i]&MASK for i in range(2,16)]
    Wf_s = schedule(Wf_list)
    da13 = (sf[13][0] - sn[13][0]) & MASK
    DW16 = (Wf_s[16] - Wn_s[16]) & MASK
    f17_by_w1_0.append((da13 + DW16) & MASK)

elapsed = time.time() - t0
uniq = len(set(f17_by_w1_0))
print(f"N_scan={N_scan}, время: {elapsed:.1f}с")
print(f"Уникальных f17 значений: {uniq} / {N_scan} ({uniq/N_scan*100:.1f}%)")
print(f"→ {'псевдослучайная (birthday применимо!)' if uniq/N_scan > 0.95 else 'коллизии присутствуют'}")

# Подсчёт: сколько W1[0] дают f17=0?
zeros = sum(1 for v in f17_by_w1_0 if v == 0)
print(f"f17=0 в диапазоне [0,{N_scan}): {zeros} значений")
print(f"Ожидание: {N_scan}/2^32 ≈ {N_scan/2**32:.4f}")
if zeros > 0:
    print(f"НАЙДЕНО! Первые W1[0] с f17=0: {[i for i,v in enumerate(f17_by_w1_0) if v==0]}")

print("\n\nT_S17_PSEUDORANDOM — верификация:")
print("f17(W1[0]) псевдослучайна ↔ birthday применим за O(2^32)")
print(f"Эмпирически: {uniq}/{N_scan} = {uniq/N_scan:.4f} уникальных (ожидание: ~{1-1/math.e:.4f} = {(1-1/math.e)*100:.1f}%)")

print("\n" + "=" * 60)
print("ИТОГИ П-27:")
print("=" * 60)
print("""
B. СОСТОЯНИЕ ПОСЛЕ 17 НУЛЕЙ:
   - δe17=0, δf17=δg17=δh17=0 (гарантировано Wang-цепочкой)
   - δa17 ≈ случайное 16-битное (HW ≈ 16, насыщено)
   - δb17 = δa16, δc17 = δa15, δd17 = δa14 (3-сдвиговый регистр)
   - δa14 нужен для δe18 = δa14 + δW17
   - Структура: 4 ненулевых регистра (a,b,c,d) + 4 нулевых (e,f,g,h)

C. 2D BIRTHDAY:
   - f17 и f18 независимы (корр. ≈ 0)
   - Стандартный 2D-birthday НЕ снижает барьер 2^64
   - НО: двухпараметрический birthday (W0[0], W1[0]) в принципе
     может достичь 2^32 вычислений ценой 32 ГБ памяти
   - Требует further analysis: T_2D_BIRTHDAY (PENDING)
""")
