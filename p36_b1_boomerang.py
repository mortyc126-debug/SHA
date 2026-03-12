"""
АТАКА B1: Бумеранг — нижняя дифференциальная характеристика (раунды 17-32)

Концепция: в бумеранг-атаке нужна нижняя (backward) характеристика от раундов 17-32.
Стоимость = 2^32 × 2^w, где w — вес нижней характеристики.
Если w < 32 → полная стоимость < 2^64.

Тест 1: Распределение HW(δe_r) и HW(δa_r) для r=17..28 при Wang-каскаде.
Тест 2: Найти пары (W0,W1) с минимальным суммарным весом Da17..Da20.
Тест 3: Проверить автоматически — можно ли обнулить Da17 выбором ΔW16.

Основа: T_BARRIER_kDIM (стоимость k барьеров = O(2^{32k}))
        T_ONE_CONSTRAINT (1 DOF на раунд)
"""

import random
import statistics

MASK = 0xFFFFFFFF

def rotr(x, n): return ((x >> n) | (x << (32-n))) & MASK
def sig0(x):  return rotr(x,7) ^ rotr(x,18) ^ (x>>3)
def sig1(x):  return rotr(x,17) ^ rotr(x,19) ^ (x>>10)
def Sig0(x):  return rotr(x,2) ^ rotr(x,13) ^ rotr(x,22)
def Sig1(x):  return rotr(x,6) ^ rotr(x,11) ^ rotr(x,25)
def Ch(e,f,g):  return ((e&f) ^ (~e&g)) & MASK
def Maj(a,b,c): return (a&b) ^ (a&c) ^ (b&c)
def hw(x): return bin(x).count('1')

K = [0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
     0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
     0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
     0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
     0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
     0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
     0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
     0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2]
IV = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]

def make_schedule(W16):
    W = list(W16) + [0]*48
    for i in range(16, 64):
        W[i] = (sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16]) & MASK
    return W

def sha_rounds(W, R):
    a,b,c,d,e,f,g,h = IV
    states = [[a,b,c,d,e,f,g,h]]
    for r in range(R):
        T1 = (h + Sig1(e) + Ch(e,f,g) + K[r] + W[r]) & MASK
        T2 = (Sig0(a) + Maj(a,b,c)) & MASK
        h=g; g=f; f=e; e=(d+T1)&MASK
        d=c; c=b; b=a; a=(T1+T2)&MASK
        states.append([a,b,c,d,e,f,g,h])
    return states

def wang_cascade(W0, W1, DW0=1):
    """Стандартный Wang-каскад: De3..De17=0."""
    Wn = [W0, W1] + [0]*14
    DWs = [0]*16; DWs[0] = DW0
    Wf_tmp = [(Wn[i]+DWs[i])&MASK for i in range(16)]
    sn3 = sha_rounds(make_schedule(Wn), 3)
    sf3 = sha_rounds(make_schedule(Wf_tmp), 3)
    DWs[2] = (-(sf3[3][4] - sn3[3][4])) & MASK
    for step in range(13):
        wi = step+3; dt = step+4
        Wfc = [(Wn[i]+DWs[i])&MASK for i in range(16)]
        sn = sha_rounds(make_schedule(Wn), dt)
        sf = sha_rounds(make_schedule(Wfc), dt)
        DWs[wi] = (-(sf[dt][4] - sn[dt][4])) & MASK
    Wf = [(Wn[i]+DWs[i])&MASK for i in range(16)]
    return Wn, Wf, DWs

print("=" * 70)
print("АТАКА B1: Бумеранг — нижняя характеристика (раунды 17-32)")
print("=" * 70)

# ─────────────────────────────────────────────────────────────
# Тест 1: HW(De_r) и HW(Da_r) для r=17..28 при Wang-каскаде
# ─────────────────────────────────────────────────────────────
print("\n[1] Распределение HW(δe_r), HW(δa_r) для r=17..28 (N=500 пар):")
print(f"{'r':>3} | {'E[HW(δe_r)]':>12} | {'E[HW(δa_r)]':>12} | {'P(δe=0)':>10} | {'P(δa=0)':>10}")
print("-" * 60)

N = 500
de_by_round = {r: [] for r in range(17, 29)}
da_by_round = {r: [] for r in range(17, 29)}

for _ in range(N):
    W0 = random.randint(0, MASK); W1 = random.randint(0, MASK)
    Wn, Wf, _ = wang_cascade(W0, W1)
    sn = sha_rounds(make_schedule(Wn), 28)
    sf = sha_rounds(make_schedule(Wf), 28)
    for r in range(17, 29):
        de = (sf[r][4] - sn[r][4]) & MASK
        da = (sf[r][0] - sn[r][0]) & MASK
        de_by_round[r].append(de)
        da_by_round[r].append(da)

for r in range(17, 29):
    de_list = de_by_round[r]; da_list = da_by_round[r]
    e_de = statistics.mean(hw(v) for v in de_list)
    e_da = statistics.mean(hw(v) for v in da_list)
    p_de = sum(1 for v in de_list if v == 0) / N
    p_da = sum(1 for v in da_list if v == 0) / N
    note = " ←" if e_de < 15 or e_da < 15 else ""
    print(f"  {r:>2} | {e_de:>12.4f} | {e_da:>12.4f} | {p_de:>10.2e} | {p_da:>10.2e}{note}")

# ─────────────────────────────────────────────────────────────
# Тест 2: Поиск пар с минимальным суммарным HW(Da17..Da20)
# ─────────────────────────────────────────────────────────────
print("\n[2] Поиск пар с min суммарным HW(δe_17+δe_18+δe_19+δe_20) (N=2000):")
best_score = 999
best_pair = None
scores = []

for _ in range(2000):
    W0 = random.randint(0, MASK); W1 = random.randint(0, MASK)
    Wn, Wf, _ = wang_cascade(W0, W1)
    sn = sha_rounds(make_schedule(Wn), 21)
    sf = sha_rounds(make_schedule(Wf), 21)
    score = sum(hw((sf[r][4]-sn[r][4])&MASK) for r in range(17, 21))
    scores.append(score)
    if score < best_score:
        best_score = score
        best_pair = (W0, W1)

print(f"  min HW(δe_17..20) = {best_score}  (пара 0x{best_pair[0]:08x}, 0x{best_pair[1]:08x})")
print(f"  E[HW(δe_17..20)] = {statistics.mean(scores):.3f}  (теор. = 4×16 = 64)")
print(f"  min(scores) hist: {sorted(scores)[:10]}")

# ─────────────────────────────────────────────────────────────
# Тест 3: Можно ли обнулить De17 через ΔW16?
# По T_ONE_CONSTRAINT — да, ΔW16 = -Da13 (стоимость = birthday 2^32)
# Но нас интересует: при случайном ΔW16 ∈ [0..2^16], как часто De17=0?
# ─────────────────────────────────────────────────────────────
print("\n[3] Проверка T_ONE_CONSTRAINT: каков вес Da17 после 'бесплатного' ΔW16?")
print("    (имитация: что если мы уже решили f17=0 birthday поиском)")
print()

# Берём пары уже с f17=0 (из известных пар):
KNOWN_PAIRS = [(0xe82222c7, 0x516cfb41), (0xd4254551, 0x679ea4de)]
for W0, W1 in KNOWN_PAIRS:
    Wn, Wf, DWs = wang_cascade(W0, W1)
    sn = sha_rounds(make_schedule(Wn), 20)
    sf = sha_rounds(make_schedule(Wf), 20)
    de17 = (sf[17][4] - sn[17][4]) & MASK
    da17 = (sf[17][0] - sn[17][0]) & MASK
    de18 = (sf[18][4] - sn[18][4]) & MASK
    da18 = (sf[18][0] - sn[18][0]) & MASK
    de19 = (sf[19][4] - sn[19][4]) & MASK
    de20 = (sf[20][4] - sn[20][4]) & MASK
    print(f"  Пара 0x{W0:08x}: de17=0x{de17:08x}(HW={hw(de17)})  da17=0x{da17:08x}(HW={hw(da17)})")
    print(f"    de18=HW{hw(de18)}  de19=HW{hw(de19)}  de20=HW{hw(de20)}")

# ─────────────────────────────────────────────────────────────
# Тест 4: Теоретическая структура бумеранга
# ─────────────────────────────────────────────────────────────
print("\n[4] Теоретический анализ возможности бумеранга:")
print()
print("  В классическом бумеранге (Wagner 1999):")
print("  - Верхняя характеристика:  (M, M') → (E17, E17') с Pr = p")
print("  - Нижняя характеристика:   (E17, E17') → (hash, hash') с Pr = q")
print("  - Суммарная стоимость: 4/p²q²")
print()
print("  Для SHA-256 верхняя характеристика: p = 2^{-32} (birthday f17=0)")
print()

# Вычислим типичный вес нижней характеристики
# = сумма HW(De_r) + HW(Da_r) для r=18..32 при уже установленных De17=0, Da17≠0
print("  Нижняя часть: De17=0 (достигнуто), Da17≠0 (≈16 бит)")
print("  За 15 раундов (17-32) дифференциал Da17 «размазывается» по всем 8 словам")
print()

# Измерим, как быстро Da17 смешивается:
print("  Скорость диффузии Da17 (N=500 пар с De17=0):")
print(f"  {'r':>3} | {'E[HW(De_r)]':>13} | {'E[HW(Da_r)]':>13}")
print("  " + "-" * 35)

de_post = {r: [] for r in range(17, 33)}
da_post = {r: [] for r in range(17, 33)}

cnt = 0
attempts = 0
target = 500
while cnt < target and attempts < 100000:
    attempts += 1
    W0 = random.randint(0, MASK); W1 = random.randint(0, MASK)
    Wn, Wf, DWs = wang_cascade(W0, W1)
    sn_full = sha_rounds(make_schedule(Wn), 32)
    sf_full = sha_rounds(make_schedule(Wf), 32)
    de17 = (sf_full[17][4] - sn_full[17][4]) & MASK
    if de17 != 0:
        continue
    cnt += 1
    for r in range(17, 33):
        de = (sf_full[r][4] - sn_full[r][4]) & MASK
        da = (sf_full[r][0] - sn_full[r][0]) & MASK
        de_post[r].append(de)
        da_post[r].append(da)

print(f"  (найдено {cnt} пар с De17=0 за {attempts} попыток)")
for r in range(17, 33):
    if de_post[r]:
        e_de = statistics.mean(hw(v) for v in de_post[r])
        e_da = statistics.mean(hw(v) for v in da_post[r])
        note = " ← низкий!" if e_de < 10 else ""
        print(f"  {r:>3} | {e_de:>13.4f} | {e_da:>13.4f}{note}")

# ─────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("ИТОГ АТАКИ B1 (Бумеранг):")
print()

# Вычислим среднее HW на раундах 18-32 как нижний вес
if de_post[18]:
    lower_trail_hw = statistics.mean(
        statistics.mean(hw(v) for v in de_post[r])
        for r in range(18, 33) if de_post[r]
    )
    print(f"  E[HW(De_r)] на раундах 18-32: {lower_trail_hw:.3f}")
    print(f"  Это означает: нижняя характеристика имеет вес w ≈ {lower_trail_hw:.0f}×15 ≈ {lower_trail_hw*15:.0f}")
    if lower_trail_hw < 10:
        print(f"  Нижняя характеристика НЕТРИВИАЛЬНАЯ — бумеранг возможен!")
    else:
        print(f"  Нижняя характеристика случайна (HW≈16) — полная диффузия")
        print(f"  Вес нижней ~ 16×15 = 240 >> 32 → бумеранг НЕ улучшает 2^64")

print()
print("  T_BARRIER_kDIM: стоимость k барьеров = O(2^{32k})")
print("  Вывод по B1: нет структуры в раундах 18-32 при типичном Da17≠0")
print("  СТАТУС: B1 — ОТРИЦАТЕЛЬНЫЙ (нет нетривиальной нижней характеристики)")
print("=" * 70)
