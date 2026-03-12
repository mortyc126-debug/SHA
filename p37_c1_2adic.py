"""
АТАКА C1: 2-адический анализ Da13 и f17

Идея: Если v_2(Da13) ≥ k (т.е. Da13 ≡ 0 mod 2^k) с вероятностью > 2^{-k},
то birthday для f17 = Da13 + δW16 = 0 обходится дешевле:
  - Нужно δW16 ≡ -Da13 mod 2^32
  - Если Da13 всегда чётна (v_2≥1): δW16 тоже чётна → поиск в 2^31 значениях

Тест 1: Распределение v_2(Da13) — сравнение с теоретическим Геом(1/2).
         Геом(1/2): P(v_2=k) = 2^{-(k+1)}, E[v_2] = 1.

Тест 2: Распределение v_2(f17) для случайных пар.
         Если E[v_2(f17)] > 1 → Da13 и δW16 взаимодействуют нетривиально.

Тест 3: Условное распределение v_2(Da13) при фиксированном W0 mod 2.
         Влияет ли чётность W0 на чётность Da13?

Тест 4: v_2(Da13) при Da1=1 (аддитивный) — есть ли связь?
         Da1 = a1(W0+1) - a1(W0). Проверить v_2(Da13) | Da1=1.

Основа: T_CARRY_DIST_GEOMETRIC, T_BARRIER_UNIFORM, T_MOD_ADD_DIFF_PROB
"""

import random
import statistics
from collections import Counter
import math

MASK = 0xFFFFFFFF

def rotr(x, n): return ((x >> n) | (x << (32-n))) & MASK
def sig0(x):  return rotr(x,7) ^ rotr(x,18) ^ (x>>3)
def sig1(x):  return rotr(x,17) ^ rotr(x,19) ^ (x>>10)
def Sig0(x):  return rotr(x,2) ^ rotr(x,13) ^ rotr(x,22)
def Sig1(x):  return rotr(x,6) ^ rotr(x,11) ^ rotr(x,25)
def Ch(e,f,g):  return ((e&f) ^ (~e&g)) & MASK
def Maj(a,b,c): return (a&b) ^ (a&c) ^ (b&c)
def hw(x): return bin(x).count('1')
def v2(x):
    """2-адическая валюация: число trailing нулей."""
    if x == 0: return 32
    k = 0
    while (x >> k) & 1 == 0: k += 1
    return k

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
    sts = [[a,b,c,d,e,f,g,h]]
    for r in range(R):
        T1=(h+Sig1(e)+Ch(e,f,g)+K[r]+W[r])&MASK; T2=(Sig0(a)+Maj(a,b,c))&MASK
        h=g;g=f;f=e;e=(d+T1)&MASK;d=c;c=b;b=a;a=(T1+T2)&MASK
        sts.append([a,b,c,d,e,f,g,h])
    return sts

def compute_da13_f17(W0, W1, j=0):
    DW0 = 1 << j
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
    Wn_s = make_schedule(Wn); Wf_s = make_schedule(Wf)
    sn = sha_rounds(Wn_s, 14); sf = sha_rounds(Wf_s, 14)
    da13 = (sf[13][0] - sn[13][0]) & MASK
    dw16 = (Wf_s[16] - Wn_s[16]) & MASK
    f17 = (da13 + dw16) & MASK
    return da13, dw16, f17

print("=" * 70)
print("АТАКА C1: 2-адический анализ Da13 и f17")
print("=" * 70)

N = 100_000

# ─────────────────────────────────────────────────────────────
# Тест 1: Распределение v_2(Da13) vs геометрическое
# ─────────────────────────────────────────────────────────────
print(f"\n[1] Распределение v_2(Da13) и v_2(f17) (N={N:,}):")
print(f"    Теоретическое (случайное): P(v_2=k) = 2^{{-(k+1)}}")
print(f"{'k':>3} | {'P(v2(Da13)=k)':>15} | {'Теор.':>8} | {'P(v2(f17)=k)':>15} | {'Теор.':>8}")
print("-" * 60)

v2_da13 = Counter()
v2_f17 = Counter()
v2_dw16 = Counter()
da13_odd = 0
da13_even = 0

for _ in range(N):
    W0 = random.randint(0, MASK); W1 = random.randint(0, MASK)
    da13, dw16, f17 = compute_da13_f17(W0, W1)
    v2_da13[v2(da13)] += 1
    v2_f17[v2(f17)] += 1
    v2_dw16[v2(dw16)] += 1
    if da13 & 1: da13_odd += 1
    else: da13_even += 1

for k in range(8):
    p_da13 = v2_da13[k] / N
    p_f17 = v2_f17[k] / N
    p_theor = 0.5**(k+1)
    d1 = p_da13 - p_theor
    d2 = p_f17 - p_theor
    note = " ←!" if abs(d1) > 0.005 or abs(d2) > 0.005 else ""
    print(f"  {k:>2} | {p_da13:>15.6f} | {p_theor:>8.6f} | {p_f17:>15.6f} | {p_theor:>8.6f}{note}")

# v2=32 (zero)
k = 32
p_da13 = v2_da13[k] / N
p_f17 = v2_f17[k] / N
print(f"  32 | {p_da13:>15.6f} | {'2^-32':>8} | {p_f17:>15.6f} | {'2^-32':>8}")

print(f"\n  P(Da13 нечётна) = {da13_odd/N:.6f}  (теор. = 0.5)")
print(f"  P(Da13 чётна)  = {da13_even/N:.6f}  (теор. = 0.5)")

# ─────────────────────────────────────────────────────────────
# Тест 2: Зависимость v_2(Da13) от v_2(W0) и v_2(W1)
# ─────────────────────────────────────────────────────────────
print(f"\n[2] E[v_2(Da13)] при фиксированных v_2(W0)=m и v_2(W1)=n:")
print(f"    Если E[v_2(Da13)] > 1 при некоторых m,n → можно выбрать W0,W1 для ускорения")
print(f"{'v2(W0)':>7} | {'E[v2(Da13)]':>13} | {'N пар':>7}")
print("-" * 35)

cond_v2 = {m: [] for m in range(6)}
for _ in range(200_000):
    W0 = random.randint(0, MASK); W1 = random.randint(0, MASK)
    m = v2(W0)
    if m < 6:
        da13, _, _ = compute_da13_f17(W0, W1)
        cond_v2[m].append(v2(da13))

for m in range(6):
    if cond_v2[m]:
        e_v2 = statistics.mean(cond_v2[m])
        note = " ← СТРУКТУРА!" if e_v2 > 1.2 else ""
        print(f"  v2(W0)={m:>1} | {e_v2:>13.4f} | {len(cond_v2[m]):>7}{note}")

# ─────────────────────────────────────────────────────────────
# Тест 3: v_2(Da13) при W0 ≡ 0 mod 2^k (контролируем младшие биты W0)
# ─────────────────────────────────────────────────────────────
print(f"\n[3] Условное E[v_2(Da13)] при W0 ≡ 0 mod 2^k:")
print(f"    (Устанавливаем k младших битов W0 в 0)")
print(f"{'k':>3} | {'E[v2(Da13)]':>13} | {'E[v2(f17)]':>13} | {'P(f17=0)':>10}")
print("-" * 50)

N2 = 50_000
for k in range(0, 9):
    mask_lo = MASK ^ ((1 << k) - 1)  # зануляем k младших бит W0
    v2_da_list = []; v2_f17_list = []; zeros = 0
    for _ in range(N2):
        W0_base = random.randint(0, MASK) & mask_lo
        W1 = random.randint(0, MASK)
        da13, dw16, f17 = compute_da13_f17(W0_base, W1)
        v2_da_list.append(v2(da13))
        v2_f17_list.append(v2(f17))
        if f17 == 0: zeros += 1
    e_da = statistics.mean(v2_da_list)
    e_f17 = statistics.mean(v2_f17_list)
    p0 = zeros / N2
    note = " ← СТРУКТУРА!" if e_da > 1.5 or e_f17 > 1.5 else ""
    print(f"  {k:>2} | {e_da:>13.4f} | {e_f17:>13.4f} | {p0:>10.2e}{note}")

# ─────────────────────────────────────────────────────────────
# Тест 4: v_2(Da_r) по раундам — когда теряется 2-адическая структура?
# ─────────────────────────────────────────────────────────────
print(f"\n[4] E[v_2(Da_r)] по раундам r=1..13 (N=20000):")
print(f"    (Отслеживаем, как 2-адическая структура Da1=1 распространяется)")
print(f"{'r':>3} | {'E[v2(Da_r)]':>13} | {'P(v2≥2)':>10} | Теор.геом")
print("-" * 45)

N3 = 20_000
da_vals_by_round = {r: [] for r in range(1, 14)}

for _ in range(N3):
    W0 = random.randint(0, MASK); W1 = random.randint(0, MASK)
    Wn = [W0, W1] + [0]*14
    Wf = [(W0+1)&MASK, W1] + [0]*14  # DW0=1
    sn = sha_rounds(make_schedule(Wn), 13)
    sf = sha_rounds(make_schedule(Wf), 13)
    for r in range(1, 14):
        da_r = (sf[r][0] - sn[r][0]) & MASK
        da_vals_by_round[r].append(da_r)

for r in range(1, 14):
    vals = da_vals_by_round[r]
    e_v2 = statistics.mean(v2(x) for x in vals)
    p_ge2 = sum(1 for x in vals if v2(x) >= 2) / N3
    note = " ← структура!" if e_v2 > 1.2 else ""
    print(f"  {r:>2} | {e_v2:>13.4f} | {p_ge2:>10.4f} | 1.000 (геом){note}")

# ─────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("ИТОГ АТАКИ C1 (2-адика):")
# Проверим общий E[v_2(Da13)]
e_v2_da13 = statistics.mean(v2_da13[k] * k for k in v2_da13) / N  # = E[v_2]
# Нет, правильно:
e_v2_total = sum(k * v2_da13[k] for k in v2_da13) / N
print(f"  E[v_2(Da13)] = {e_v2_total:.4f}  (теор. геом(1/2) = 1.000)")
if abs(e_v2_total - 1.0) > 0.05:
    print(f"  ОТКЛОНЕНИЕ {e_v2_total-1.0:+.4f} — 2-адическая структура Da13!")
    print(f"  СТАТУС: C1 — ПЕРСПЕКТИВНЫЙ")
else:
    print(f"  Da13 имеет стандартное 2-адическое распределение (геом(1/2))")
    print(f"  f17 = Da13 + δW16: сумма двух равномерных → равномерная")
    print(f"  СТАТУС: C1 — ОТРИЦАТЕЛЬНЫЙ")
print("=" * 70)
