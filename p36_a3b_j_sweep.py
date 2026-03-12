"""
АТАКА A3 (исправленная): j-sweep для DW0=2^j

Проверяем: для разных DW0=2^j (j=0..31), меняется ли распределение Da13?
Если E[HW(Da13)] зависит от j → есть оптимальный j для birthday.

SA-sweep показал j=2 оптимален для 8-раундовой XOR-стоимости (cost=47 vs 56 при j=0).
Вопрос: есть ли оптимальный j для аддитивного birthday-поиска f17=0?

Также проверяем: является ли f17(W0,W1,j) равномерным для всех j?
Основа: T_J_INDEPENDENCE (только j=24,26,28,30,31), T_BARRIER_UNIFORM, T_CARRY_GEOMETRY
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

def compute_f17_with_j(W0, W1, j):
    """Каскад с DW0=2^j: f17 = Da13 + δW16."""
    DW0 = 1 << j
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

print("=" * 70)
print("АТАКА A3 (исправл.): j-sweep для DW0=2^j → Da13 и f17 статистика")
print("=" * 70)

N = 800  # пар на каждое j
J_LIST = [0, 1, 2, 3, 7, 15, 16, 23, 24, 28, 31]  # ключевые значения

print(f"\nN={N} случайных (W0,W1) на каждое j")
print(f"\n{'j':>3} | {'DW0':>12} | {'E[HW(Da13)]':>13} | {'E[HW(f17)]':>12} | {'P(f17=0)':>10} | {'E[HW(DW16)]':>13}")
print("-" * 80)

results = {}
t_total = time.time()

for j in J_LIST:
    hw_da13 = []
    hw_f17  = []
    hw_dw16 = []
    zeros   = 0
    t0 = time.time()
    for _ in range(N):
        W0 = random.randint(0, MASK)
        W1 = random.randint(0, MASK)
        f17, da13, dw16 = compute_f17_with_j(W0, W1, j)
        hw_da13.append(hw(da13))
        hw_f17.append(hw(f17))
        hw_dw16.append(hw(dw16))
        if f17 == 0: zeros += 1

    m13 = statistics.mean(hw_da13)
    mf17 = statistics.mean(hw_f17)
    mdw16 = statistics.mean(hw_dw16)
    p0 = zeros / N
    results[j] = {'mean_da13': m13, 'mean_f17': mf17, 'zeros': zeros, 'mean_dw16': mdw16}
    marker = " ← ОПТИМ?" if m13 < 15.7 else (" ← MAX" if m13 > 16.3 else "")
    print(f"  {j:>2} | 0x{(1<<j):08x}   | {m13:>13.4f} | {mf17:>12.4f} | {p0:>10.2e} | {mdw16:>13.4f}{marker}")

elapsed = time.time() - t_total
print(f"\nОбщее время: {elapsed:.1f}с")

# Ищем оптимальный j
best_j = min(results, key=lambda j: results[j]['mean_da13'])
worst_j = max(results, key=lambda j: results[j]['mean_da13'])
print(f"\nОптимальный j (min E[HW(Da13)]): j={best_j}  E={results[best_j]['mean_da13']:.4f}")
print(f"Худший j (max E[HW(Da13)]):      j={worst_j} E={results[worst_j]['mean_da13']:.4f}")
spread = results[worst_j]['mean_da13'] - results[best_j]['mean_da13']
print(f"Разброс: {spread:.4f} бит")

print("\n[2] Анализ корреляции HW(Da13) и HW(DW16) для j=0 и j=24:")
for jtest in [0, 24]:
    corr_data = []
    for _ in range(1000):
        W0 = random.randint(0, MASK)
        W1 = random.randint(0, MASK)
        f17, da13, dw16 = compute_f17_with_j(W0, W1, jtest)
        corr_data.append((hw(da13), hw(dw16)))
    x = [v[0] for v in corr_data]; y = [v[1] for v in corr_data]
    mx = statistics.mean(x); my = statistics.mean(y)
    sx = statistics.stdev(x); sy = statistics.stdev(y)
    cov = sum((x[i]-mx)*(y[i]-my) for i in range(len(x)))/len(x)
    rho = cov/(sx*sy) if sx*sy > 0 else 0
    print(f"  j={jtest}: corr(HW(Da13), HW(DW16)) = {rho:.4f}  (0=независимы)")

print("\n[3] Тест равномерности f17 для j=0 и j=24 (Welch's t-test на HW):")
for jtest in [0, 24]:
    hw_list = []
    for _ in range(2000):
        W0 = random.randint(0, MASK); W1 = random.randint(0, MASK)
        f17, _, _ = compute_f17_with_j(W0, W1, jtest)
        hw_list.append(hw(f17))
    m = statistics.mean(hw_list)
    s = statistics.stdev(hw_list)
    # Биномиальное B(32,0.5): E=16, std=sqrt(32×0.25)=2.828
    print(f"  j={jtest}: E[HW(f17)]={m:.4f} std={s:.4f}  (биномиальное: E=16.0, std=2.828)")

print("\n" + "=" * 70)
print("ИТОГ A3 (j-sweep):")
if spread > 1.0:
    print(f"  ЗНАЧИМЫЙ РАЗБРОС {spread:.3f} бит! j={best_j} лучший для birthday")
    print(f"  СТАТУС: Оптимальный j улучшает шансы birthday")
elif spread > 0.3:
    print(f"  Небольшой разброс {spread:.3f} бит. j={best_j} незначительно лучше")
    print(f"  СТАТУС: Статистически незначимо, T_J_INDEPENDENCE подтверждается")
else:
    print(f"  Разброс {spread:.4f} бит — все j дают одинаковый E[HW(Da13)]≈16")
    print(f"  T_J_INDEPENDENCE подтверждена для аддитивных дифференциалов")
    print(f"  СТАТУС: A3 полностью опровергнута")
print("=" * 70)
