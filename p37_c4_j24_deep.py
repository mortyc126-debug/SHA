"""
АТАКА C4: Глубокое тестирование j=24 каскада

Гипотеза C4: j=24 даёт δe2 с P(δe2=0)=2^{-2} (лучший carry для XOR-следа).
Вопрос: влияет ли это на аддитивный Da13?

SA-sweep показал j=2 оптимален для 8-раундовой XOR-стоимости.
j-sweep показал разброс 0.42 бита (незначимо) для E[HW(Da13)].

Но: T_CARRY_GEOMETRY специфична для XOR-следа.
Для аддитивного следа нужно проверить непосредственно P(Da13+δW16=0).

Тест 1: Прямой birthday-поиск f17=0 с DW0=2^24 vs DW0=1 (N=10^6 пар).
         Сколько раз f17=0 встречается? Должно быть ≈ 10^6 / 2^32 ≈ 0.23.

Тест 2: Гистограмма Da13 для j=24 vs j=0 — есть ли пики?

Тест 3: Корреляция бит f17 для j=24 vs j=0 — распределение бит.

Тест 4: Проверить P(δe2=0) для j=24 (XOR-след):
         δe2 = e2(W0+2^24,W1) XOR e2(W0,W1) — должно быть P≈2^{-2}?

Основа: T_CARRY_GEOMETRY (j=24 лучший для XOR), T_J_INDEPENDENCE (аддитивные)
"""

import random
import statistics
from collections import Counter
import time

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
    sts = [[a,b,c,d,e,f,g,h]]
    for r in range(R):
        T1=(h+Sig1(e)+Ch(e,f,g)+K[r]+W[r])&MASK; T2=(Sig0(a)+Maj(a,b,c))&MASK
        h=g;g=f;f=e;e=(d+T1)&MASK;d=c;c=b;b=a;a=(T1+T2)&MASK
        sts.append([a,b,c,d,e,f,g,h])
    return sts

def compute_f17_j(W0, W1, j):
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
    return (da13 + dw16) & MASK, da13, dw16

print("=" * 70)
print("АТАКА C4: Глубокий анализ j=24 каскада vs j=0")
print("=" * 70)

# ─────────────────────────────────────────────────────────────
# Тест 1: Прямой birthday count для j=0, j=24, j=31
# ─────────────────────────────────────────────────────────────
N_BIRTHDAY = 1_000_000
print(f"\n[1] Birthday поиск f17=0: N={N_BIRTHDAY:,} пар для j∈{{0, 24, 31}}")
print(f"    Ожидается ≈{N_BIRTHDAY/2**32:.3f} нулей на j")

for j in [0, 24, 31]:
    t0 = time.time()
    zeros = 0
    hw_f17 = []
    for _ in range(N_BIRTHDAY):
        W0 = random.randint(0, MASK); W1 = random.randint(0, MASK)
        f17, _, _ = compute_f17_j(W0, W1, j)
        if f17 == 0:
            zeros += 1
        hw_f17.append(hw(f17))
    elapsed = time.time() - t0
    rate = N_BIRTHDAY / elapsed
    mean_hw = statistics.mean(hw_f17)
    print(f"  j={j:2d}: найдено {zeros} нулей  E[HW(f17)]={mean_hw:.4f}  "
          f"({elapsed:.1f}с, {rate:.0f} пар/с)")

# ─────────────────────────────────────────────────────────────
# Тест 2: XOR-след δe2 для j=24 (T_CARRY_GEOMETRY)
# ─────────────────────────────────────────────────────────────
print("\n[2] XOR-след δe2 для DW0=2^j (T_CARRY_GEOMETRY):")
print(f"    δe2_XOR = e2(W0+2^j, W1) XOR e2(W0, W1)")
print(f"    T_CARRY_GEOMETRY: j=24 → P(δe2=0)≈2^{{-2}}, т.е. ≈1/4 пар")
print(f"{'j':>3} | {'P(δe2_XOR=0)':>14} | {'E[HW(δe2_XOR)]':>16}")
print("-" * 40)

N_XOR = 50000
for j in [0, 2, 7, 24, 28, 31]:
    DW0 = 1 << j
    zeros = 0
    hw_list = []
    for _ in range(N_XOR):
        W0 = random.randint(0, MASK); W1 = random.randint(0, MASK)
        Wn = [W0, W1] + [0]*14
        Wf = [(W0 + DW0) & MASK, W1] + [0]*14
        sn2 = sha_rounds(make_schedule(Wn), 2)
        sf2 = sha_rounds(make_schedule(Wf), 2)
        de2_xor = sn2[2][4] ^ sf2[2][4]   # XOR-разность
        hw_list.append(hw(de2_xor))
        if de2_xor == 0: zeros += 1
    p0 = zeros / N_XOR
    e_hw = statistics.mean(hw_list)
    note = " ← T_CARRY_GEOMETRY!" if abs(p0 - 0.25) < 0.05 else ""
    print(f"  {j:>2} | {p0:>14.4f} | {e_hw:>16.4f}{note}")

# ─────────────────────────────────────────────────────────────
# Тест 3: Гистограмма Da13 для j=0 vs j=24 (детальная)
# ─────────────────────────────────────────────────────────────
print("\n[3] Сравнение распределения HW(Da13) для j=0 и j=24 (N=10000 каждый):")
N_HIST = 10000
for j in [0, 24]:
    hw_da13 = []
    for _ in range(N_HIST):
        W0 = random.randint(0, MASK); W1 = random.randint(0, MASK)
        _, da13, _ = compute_f17_j(W0, W1, j)
        hw_da13.append(hw(da13))
    m = statistics.mean(hw_da13)
    s = statistics.stdev(hw_da13)
    cnt = Counter(hw_da13)
    # Биномиальное B(32,0.5): E=16, peak at 16
    p_low = sum(cnt[k] for k in range(0, 13)) / N_HIST
    p_mid = sum(cnt[k] for k in range(13, 20)) / N_HIST
    print(f"  j={j}: E[HW(Da13)]={m:.4f}  std={s:.4f}  P(HW<13)={p_low:.4f}  P(13≤HW<20)={p_mid:.4f}")

# ─────────────────────────────────────────────────────────────
# Тест 4: Битовое распределение f17 для j=0 и j=24
# ─────────────────────────────────────────────────────────────
print("\n[4] Битовое распределение f17 — P(бит b = 1) для j=0 и j=24:")
print(f"    (ожидается P(bit=1)=0.5 для всех b)")
print(f"{'бит':>4} | {'P(b=1) j=0':>12} | {'P(b=1) j=24':>12} | Разница")
print("-" * 48)

N_BIT = 5000
bit_count_j0 = [0]*32
bit_count_j24 = [0]*32
for _ in range(N_BIT):
    W0 = random.randint(0, MASK); W1 = random.randint(0, MASK)
    f17_j0, _, _ = compute_f17_j(W0, W1, 0)
    f17_j24, _, _ = compute_f17_j(W0, W1, 24)
    for b in range(32):
        if (f17_j0 >> b) & 1: bit_count_j0[b] += 1
        if (f17_j24 >> b) & 1: bit_count_j24[b] += 1

max_diff = 0
for b in range(0, 32, 4):
    p0 = bit_count_j0[b] / N_BIT
    p24 = bit_count_j24[b] / N_BIT
    diff = abs(p0 - p24)
    max_diff = max(max_diff, diff)
    note = " ←!" if diff > 0.05 else ""
    print(f"  {b:>3} | {p0:>12.4f} | {p24:>12.4f} | {diff:.4f}{note}")

print(f"\n  max |P_j0 - P_j24| = {max_diff:.4f}  (>0.05 → структурное отличие)")

# ─────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("ИТОГ АТАКИ C4 (j=24 каскад):")
# birthday count сравнение уже выведено выше
print("  XOR-след δe2: проверка T_CARRY_GEOMETRY")
print("  Da13: сравнение j=0 vs j=24")
if max_diff > 0.05:
    print(f"  СТРУКТУРНОЕ ОТЛИЧИЕ в битовом распределении! max_diff={max_diff:.4f}")
    print(f"  СТАТУС: C4 — ПЕРСПЕКТИВНЫЙ")
else:
    print(f"  Нет структурных отличий j=24 vs j=0 для аддитивного следа")
    print(f"  T_J_INDEPENDENCE подтверждена для f17")
    print(f"  СТАТУС: C4 — ОТРИЦАТЕЛЬНЫЙ")
print("=" * 70)
