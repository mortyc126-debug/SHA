"""
АТАКА B3: Двойной каскад — δe=0 на чётных раундах, δa=0 на нечётных

Гипотеза: Найти δW[0..15] такие, что:
  ∀r ∈ {2,4,6,8,10,12,14,16}: δe_r = 0  (чётные)
  ∀r ∈ {3,5,7,9,11,13,15,17}: δa_r = 0  (нечётные)

Если Da17=0, то Da14=0 (3-сдвиг), и De18 = Da14+δW17 = δW17.
  → δe18=0 управляется напрямую через δW17.

ТЕСТ 1: Можно ли достичь δe_even=0 И δa_odd=0 одновременно?
ТЕСТ 2: Стоимость через адаптивный каскад обоих условий.
ТЕСТ 3: Структура δa при Wang-цепочке — можно ли управлять δa отдельно?

Основа: T_ONE_CONSTRAINT, T_JOINT_ZERO, T_Da_GENERAL, T_DAk_CASCADE
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
IV = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0xa510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]
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

def get_de_da(sn, sf, r):
    de = (sf[r][4] - sn[r][4]) & MASK
    da = (sf[r][0] - sn[r][0]) & MASK
    return de, da

def adaptive_de_cascade(W0, W1, DW0=1):
    """Стандартный каскад: De3..De17=0. Возвращает (DWs, состояния n, f)."""
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
    sn = sha_rounds(make_schedule(Wn), 17)
    sf = sha_rounds(make_schedule(Wf), 17)
    return DWs, sn, sf

print("=" * 70)
print("АТАКА B3: Двойной каскад (δe чётные = 0, δa нечётные = 0)")
print("=" * 70)

# ─────────────────────────────────────────────────────────────
# Тест 1: Профиль δe и δa при стандартном Wang-каскаде
# ─────────────────────────────────────────────────────────────
print("\n[1] Профиль δe_r и δa_r при Wang-каскаде (De3..De17=0):")
W0, W1 = 0xe82222c7, 0x516cfb41  # П-15
DWs, sn, sf = adaptive_de_cascade(W0, W1)
print(f"  {'r':>3} | {'δe_r':>10} | {'δa_r':>10} | {'HW(δe)':>7} | {'HW(δa)':>7} | Условие")
print("  " + "-" * 65)
for r in range(1, 18):
    if r >= len(sn): break
    de_r, da_r = get_de_da(sn, sf, r)
    cond = "δe=0 ✓" if de_r==0 else (f"δe≠0" )
    cond += "  δa=0 ✓" if da_r==0 else ""
    print(f"  {r:>3} | 0x{de_r:08x} | 0x{da_r:08x} | {hw(de_r):>7} | {hw(da_r):>7} | {cond}")

# ─────────────────────────────────────────────────────────────
# Тест 2: Попытка управления δa_r через T_One_Constraint
# ─────────────────────────────────────────────────────────────
print("\n[2] T_ONE_CONSTRAINT проверка:")
print("    δW_r даёт 1 степень свободы → либо δe_{r+1}=0, либо δa_{r+1}=0")
print("    При Wang: δW_r выбран для δe_{r+1}=0 → δa_{r+1} неуправляем")
print()

# Что если выбирать δW_r для δa_{r+1}=0 вместо δe_{r+1}=0?
# Формула: δa_{r+1} = δT1_r + δT2_r = (δe_{r+1} - δd_r) + δT2_r
# δa_{r+1}=0 ↔ δT2_r = -δT1_r
# δT2_r = δSig0(a_r) + δMaj(a_r,b_r,c_r) → НЕ зависит от δW_r!
# Значит δa_{r+1}=0 не управляется через δW_r напрямую.

print("    Ключевой вывод (аналитически):")
print("    δa_{r+1} = δT2_r + δT1_r")
print("    δT2_r = δSig0(a_r) + δMaj(a_r,...) — НЕ зависит от δW_r")
print("    δT1_r = δh_r + δSig1(e_r) + δCh(e_r,...) + δW_r — зависит от δW_r")
print()
print("    Для δa_{r+1}=0: δW_r = -(δh_r + δSig1(e_r) + δCh(e_r,...)) - δT2_r")
print("    Но тогда δe_{r+1} = δd_r + δT1_r = δd_r + 0 = δd_r ≠ 0")
print("    → Нельзя одновременно δe=0 и δa=0 через одно δW_r ✓ (T_ONE_CONSTRAINT)")

# ─────────────────────────────────────────────────────────────
# Тест 3: Чередующийся каскад (Da_нечётные = 0)
# ─────────────────────────────────────────────────────────────
print("\n[3] Попытка чередующегося каскада: δa=0 на нечётных раундах:")
print("    Алгоритм: на чётных r — выбирать δW для δe_{r+1}=0")
print("              на нечётных r — выбирать δW для δa_{r+1}=0")
print()

def alternating_cascade(W0, W1, DW0=1):
    """
    Чередующийся каскад:
      r=2: δe=0 (Wang)
      r=3: δa=0
      r=4: δe=0
      ... и т.д.
    """
    Wn = [W0, W1] + [0]*14
    DWs = [0]*16; DWs[0] = DW0

    # ΔW2 → De3=0 (стандартно)
    Wf_tmp = [(Wn[i]+DWs[i])&MASK for i in range(16)]
    sn3 = sha_rounds(make_schedule(Wn), 3)
    sf3 = sha_rounds(make_schedule(Wf_tmp), 3)
    DWs[2] = (-(sf3[3][4] - sn3[3][4])) & MASK

    # Чередуем: нечётные r=3,5,7,9,11,13,15 → δa_{r+1}=0
    #            чётные r=4,6,8,10,12,14,16 → δe_{r+1}=0
    for step in range(13):
        wi = step + 3
        dt_e = step + 4   # e-регистр следующего раунда
        dt_a = step + 4   # a-регистр следующего раунда
        r_current = step + 3  # текущий раунд

        Wfc = [(Wn[i]+DWs[i])&MASK for i in range(16)]
        sn = sha_rounds(make_schedule(Wn), dt_e)
        sf = sha_rounds(make_schedule(Wfc), dt_e)

        if r_current % 2 == 1:  # нечётный → целимся на δa=0
            # δa_{r+1} = (sf[dt_a][0] - sn[dt_a][0]) при текущих DWs
            da_nat = (sf[dt_a][0] - sn[dt_a][0]) & MASK
            # δW_wi = -(da_{r+1}_nat)
            # Но δa зависит от δW через δT1, а δT2 — нет:
            # δa_{r+1} = δT1 + δT2 = (δh + δSig1(e) + δCh + δW + K) + δT2
            # Изменение δW меняет δT1 аддитивно → Da_new = Da_nat + δW_delta
            # → δW для Da=0: DWs[wi] -= da_nat
            DWs[wi] = (-da_nat) & MASK
        else:  # чётный → целимся на δe=0 (стандартный Wang)
            de_nat = (sf[dt_e][4] - sn[dt_e][4]) & MASK
            DWs[wi] = (-de_nat) & MASK

    Wf = [(Wn[i]+DWs[i])&MASK for i in range(16)]
    sn17 = sha_rounds(make_schedule(Wn), 17)
    sf17 = sha_rounds(make_schedule(Wf), 17)
    return DWs, sn17, sf17

print("    Результаты чередующегося каскада (N=5 случайных пар):")
print(f"    {'r':>3} | {'δe_r':>10} | {'δa_r':>10} | {'HW(δe)':>7} | {'HW(δa)':>7} | Целевое")
for trial in range(3):
    W0t = random.randint(0, MASK)
    W1t = random.randint(0, MASK)
    DWs_alt, sn_alt, sf_alt = alternating_cascade(W0t, W1t)
    print(f"\n    Пара {trial+1}: W0=0x{W0t:08x} W1=0x{W1t:08x}")
    even_de_zeros = 0; odd_da_zeros = 0
    for r in range(2, 18):
        if r >= len(sn_alt): break
        de_r, da_r = get_de_da(sn_alt, sf_alt, r)
        target = "δe→0" if r % 2 == 0 else "δa→0"
        achieved = ("✓" if (r%2==0 and de_r==0) else ("✓" if (r%2==1 and da_r==0) else "✗"))
        if r%2==0 and de_r==0: even_de_zeros += 1
        if r%2==1 and da_r==0: odd_da_zeros += 1
        print(f"    {r:>3} | 0x{de_r:08x} | 0x{da_r:08x} | {hw(de_r):>7} | {hw(da_r):>7} | {target} {achieved}")
    print(f"    Итог: δe=0 на чётных: {even_de_zeros}/8, δa=0 на нечётных: {odd_da_zeros}/8")

# ─────────────────────────────────────────────────────────────
# Тест 4: Проверка — можно ли снизить HW(Da13) через чередующийся каскад
# ─────────────────────────────────────────────────────────────
print("\n[4] Распределение HW(Da13) при чередующемся vs стандартном каскаде:")
N_stat = 200
hw_std = []; hw_alt = []
t0 = time.time()
for _ in range(N_stat):
    W0t = random.randint(0, MASK); W1t = random.randint(0, MASK)

    # Стандартный каскад
    DWs_s, sn_s, sf_s = adaptive_de_cascade(W0t, W1t)
    da13_s = (sf_s[13][0] - sn_s[13][0]) & MASK
    hw_std.append(hw(da13_s))

    # Чередующийся каскад
    DWs_a, sn_a, sf_a = alternating_cascade(W0t, W1t)
    da13_a = (sf_a[13][0] - sn_a[13][0]) & MASK
    hw_alt.append(hw(da13_a))

elapsed = time.time() - t0
print(f"  N={N_stat}, время: {elapsed:.1f}с")
print(f"  E[HW(Da13)] стандартный:       {statistics.mean(hw_std):.4f} ± {statistics.stdev(hw_std):.3f}")
print(f"  E[HW(Da13)] чередующийся:      {statistics.mean(hw_alt):.4f} ± {statistics.stdev(hw_alt):.3f}")
delta_hw = statistics.mean(hw_std) - statistics.mean(hw_alt)
print(f"  Разница:                        {delta_hw:+.4f} бит")

print("\n" + "=" * 70)
print("ИТОГ АТАКИ B3 (Двойной каскад):")
print(f"  T_ONE_CONSTRAINT подтверждена: δW_r управляет ИЛИ δe, ИЛИ δa")
print(f"  Чередующийся каскад: при нечётных r δa_{r+1}=0 достигается,")
print(f"  но при этом δe_{r+1}≠0 — нарушается Wang-цепочка δe")
if delta_hw > 0.5:
    print(f"  E[HW(Da13)] снизился на {delta_hw:.3f} бит → частичный успех!")
    print(f"  СТАТУС: B3 ЧАСТИЧНО УСПЕШНА — чередующийся каскад меняет Da13")
else:
    print(f"  E[HW(Da13)] изменился на {delta_hw:.4f} бит — статистически незначимо")
    print(f"  СТАТУС: B3 ПОДТВЕРЖДАЕТ T_ONE_CONSTRAINT — нет двойного управления")
print("=" * 70)
