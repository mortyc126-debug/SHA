"""
АТАКА A3: Малое δa1 через T_HW_DA1 → проверка влияния на Da13

Гипотеза: Выбор W0 с HW(Da1)=1 (trailing_ones=0) → Da1 маленький →
  рекуррентность Da13 = ΔT2_12 - ΔT2_8 + ΔT2_4 - Da1 (T_Da13_ALTERNATING)
  → HW(Da13) < 16 → P(f17=0) > 2^{-32}?

T_HW_DA1: HW(Da1) = trailing_ones_from_bit_j(BASE_A1 + W0) + 1
  BASE_A1 = 0xfc08884d
  trailing_ones = 0 → HW(Da1) = 1  (однобитовый Da1)
  P(trailing_ones=0) ≈ 1/2

Основа: T_HW_DA1, T_Da13_ALTERNATING, T_BARRIER_UNIFORM, T_CARRY_DIST_GEOMETRIC
"""

import random
import statistics
import time

MASK = 0xFFFFFFFF

def rotr(x, n): return ((x >> n) | (x << (32-n))) & MASK
def sig0(x):  return rotr(x,7)  ^ rotr(x,18) ^ (x>>3)
def sig1(x):  return rotr(x,17) ^ rotr(x,19) ^ (x>>10)
def Sig0(x):  return rotr(x,2)  ^ rotr(x,13) ^ rotr(x,22)
def Sig1(x):  return rotr(x,6)  ^ rotr(x,11) ^ rotr(x,25)
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

BASE_A1 = 0xfc08884d  # S_base из T_HW_DA1

def trailing_ones_from_bit(x, j):
    """Число последовательных единиц начиная с бита j в x."""
    count = 0
    while (x >> (j + count)) & 1:
        count += 1
    return count

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

def compute_da1(W0, DW0=1):
    """Вычислить Da1 = a1(W0+DW0) - a1(W0) (mod 2^32)."""
    Wn = [W0] + [0]*15
    Wf = [(W0+DW0)&MASK] + [0]*15
    sn = sha_rounds(make_schedule(Wn), 1)
    sf = sha_rounds(make_schedule(Wf), 1)
    return (sf[1][0] - sn[1][0]) & MASK

def compute_f17_da13(W0, W1, DW0=1):
    """Адаптивный каскад → f17 = Da13 + δW16."""
    Wn = [W0, W1] + [0]*14
    DWs = [0]*16
    DWs[0] = DW0

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
    DW16 = (Wf_s[16] - Wn_s[16]) & MASK
    f17 = (da13 + DW16) & MASK
    return f17, da13, DW16

# ─────────────────────────────────────────────────────────────
print("=" * 70)
print("АТАКА A3: Малое HW(Da1) → влияние на Da13 и f17")
print("=" * 70)

# Шаг 1: Верификация T_HW_DA1 и доля W0 с HW(Da1)=1
print("\n[1] Верификация T_HW_DA1: распределение HW(Da1) для j=0:")
DW0 = 1  # ΔW0 = 2^0 = 1
hw_dist = {}
N_verify = 10000
t0 = time.time()
for _ in range(N_verify):
    W0 = random.randint(0, MASK)
    da1 = compute_da1(W0, DW0)
    h = hw(da1)
    hw_dist[h] = hw_dist.get(h, 0) + 1

print(f"  N={N_verify}, время: {time.time()-t0:.1f}с")
print(f"  HW(Da1) | count | фракция | теория 2^(-HW)")
for h in sorted(hw_dist)[:8]:
    frac = hw_dist[h] / N_verify
    theory = 0.5**h if h > 0 else 0
    print(f"  {h:>7} | {hw_dist[h]:5d} | {frac:.4f}  | {theory:.4f}")
frac_hw1 = hw_dist.get(1, 0) / N_verify
print(f"\n  P(HW(Da1)=1) = {frac_hw1:.4f}  (теория: ~0.5)")

# Шаг 2: Метод выбора W0 с HW(Da1)=1
print("\n[2] Метод фильтрации W0 с trailing_ones(BASE_A1+W0, j=0)=0:")
count_hw1 = 0
N_filter = 5000
for _ in range(N_filter):
    W0 = random.randint(0, MASK)
    to = trailing_ones_from_bit(BASE_A1 + W0, 0)  # от бита 0
    if to == 0:
        count_hw1 += 1
print(f"  P(trailing_ones=0) = {count_hw1/N_filter:.4f}  (ожидание: ~0.5)")

# Шаг 3: ГЛАВНЫЙ ТЕСТ — сравнение E[HW(Da13)] для случайных vs HW(Da1)=1
print("\n[3] ГЛАВНЫЙ ТЕСТ: E[HW(Da13)] при случайных W0 vs HW(Da1)=1")
print("    Гипотеза: HW(Da1)=1 → HW(Da13) < 16?")

N_test = 2000
hw_da13_random = []  # случайные W0
hw_da13_hw1    = []  # W0 с HW(Da1)=1
f17_random_zero = 0
f17_hw1_zero    = 0

t0 = time.time()
collected_random = 0
collected_hw1 = 0

W1_fixed = random.randint(0, MASK)  # зафиксируем W1 для чистоты

# Собираем два набора
W0_list_rand = [random.randint(0, MASK) for _ in range(N_test)]
W0_list_hw1  = []
tries = 0
while len(W0_list_hw1) < N_test:
    W0 = random.randint(0, MASK)
    tries += 1
    if trailing_ones_from_bit(BASE_A1 + W0, 0) == 0:
        W0_list_hw1.append(W0)

print(f"  Найдено {N_test} W0 с trailing_ones=0 за {tries} попыток (эффективность: {N_test/tries:.2%})")

for W0 in W0_list_rand:
    W1 = random.randint(0, MASK)
    f17, da13, dw16 = compute_f17_da13(W0, W1)
    hw_da13_random.append(hw(da13))
    if f17 == 0: f17_random_zero += 1

for W0 in W0_list_hw1:
    W1 = random.randint(0, MASK)
    f17, da13, dw16 = compute_f17_da13(W0, W1)
    hw_da13_hw1.append(hw(da13))
    if f17 == 0: f17_hw1_zero += 1

elapsed = time.time() - t0
print(f"\n  Время: {elapsed:.1f}с")
print(f"\n  Метрика         | Случайные W0  | W0 с HW(Da1)=1")
print(f"  ────────────────────────────────────────────────")
m_rand = statistics.mean(hw_da13_random)
m_hw1  = statistics.mean(hw_da13_hw1)
s_rand = statistics.stdev(hw_da13_random)
s_hw1  = statistics.stdev(hw_da13_hw1)
print(f"  E[HW(Da13)]     | {m_rand:7.4f} ± {s_rand:.2f} | {m_hw1:7.4f} ± {s_hw1:.2f}")
print(f"  P(Da13=0)       | {hw_da13_random.count(0)/N_test:.2e}  | {hw_da13_hw1.count(0)/N_test:.2e}")
print(f"  P(f17=0)        | {f17_random_zero/N_test:.2e}  | {f17_hw1_zero/N_test:.2e}")
print(f"  Ожидание P(f17=0): {1/2**32:.2e}")

# Шаг 4: Тест при разных HW(Da1) — HW=1,2,3,4,5+
print("\n[4] Зависимость E[HW(Da13)] от HW(Da1):")
print(f"  HW(Da1) | N пар | E[HW(Da13)] | std  | P(f17=0)")
groups = {1:[], 2:[], 3:[], 4:[], 5:[]}

N_group = 3000
t0 = time.time()
for _ in range(N_group * 10):
    W0 = random.randint(0, MASK)
    W1 = random.randint(0, MASK)
    da1 = compute_da1(W0, DW0)
    h_da1 = hw(da1)
    key = min(h_da1, 5)
    if key in groups and len(groups[key]) < N_group:
        f17, da13, _ = compute_f17_da13(W0, W1)
        groups[key].append((hw(da13), f17==0))
    if all(len(v) >= N_group for v in groups.values()):
        break

elapsed = time.time() - t0
print(f"  (время сбора: {elapsed:.1f}с)")
for k in sorted(groups):
    vals = groups[k]
    if not vals: continue
    hw13_list = [v[0] for v in vals]
    zero_f17  = sum(v[1] for v in vals)
    theory_hw1 = 0.5**k
    label = f"HW={k}" if k < 5 else "HW≥5"
    m = statistics.mean(hw13_list)
    s = statistics.stdev(hw13_list) if len(hw13_list)>1 else 0
    print(f"  {label:>7} | {len(vals):5d} | {m:11.4f} | {s:4.2f} | {zero_f17/len(vals):.2e}")

# Итог
print("\n" + "=" * 70)
print("ИТОГ АТАКИ A3:")
delta = m_rand - m_hw1
if abs(delta) > 0.5:
    direction = "СНИЖЕНИЕ" if delta > 0 else "ПОВЫШЕНИЕ"
    print(f"  ОБНАРУЖЕН ЭФФЕКТ: {direction} E[HW(Da13)] на {delta:.3f} бит при HW(Da1)=1")
    print(f"  СТАТУС: Структурное влияние существует — гипотеза A3 ЧАСТИЧНО ПОДТВЕРЖДЕНА")
else:
    print(f"  E[HW(Da13)] не изменяется значимо: Δ={delta:.4f} бит")
    print(f"  T_BARRIER_UNIFORM подтверждена: Da13 не зависит от HW(Da1)")
    print(f"  СТАТУС: Гипотеза A3 ОПРОВЕРГНУТА")
print("=" * 70)
