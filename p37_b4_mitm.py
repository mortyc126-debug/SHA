"""
АТАКА B4: 2-блочный MITM (Meet-in-the-Middle)

Концепция: Если после блока 1 с Wang-каскадом (f17=0) финальное состояние
H = (a64+IV_a, b64+IV_b, ...) имеет структуру (смещение, корреляцию) —
то 2-блочный birthday может быть дешевле 2^128.

Тест 1: Распределение HW(H[i]) для i=0..7 при f17=0 vs случайные пары.
         Если E[HW] ≠ 16 → есть структура → MITM имеет смысл.

Тест 2: Корреляция между словами H[i] и H[j] при f17=0.
         corr(H[0], H[1]) ≠ 0 → зависимость → можно эксплуатировать.

Тест 3: Разность δH = H(M')-H(M) при f17=0. Если δH имеет низкий вес →
         различить коллизионные пары легче → помогает в MITM.

Тест 4: Проверить T_FULL_DIFFUSION — распределение финального хэша
         при Wang-каскаде идентично равномерному?

Основа: T_BARRIER_kDIM, T_2D_BIRTHDAY_NEGATIVE, T_BIRTHDAY_COST17
"""

import random
import statistics
from collections import Counter

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

def sha256_block(W16, state=None):
    """Один блок SHA-256: 64 раунда + добавление к состоянию (или IV)."""
    s = state if state else IV
    a,b,c,d,e,f,g,h = s
    W = make_schedule(W16)
    for r in range(64):
        T1 = (h + Sig1(e) + Ch(e,f,g) + K[r] + W[r]) & MASK
        T2 = (Sig0(a) + Maj(a,b,c)) & MASK
        h=g; g=f; f=e; e=(d+T1)&MASK
        d=c; c=b; b=a; a=(T1+T2)&MASK
    return [(a+s[0])&MASK,(b+s[1])&MASK,(c+s[2])&MASK,(d+s[3])&MASK,
            (e+s[4])&MASK,(f+s[5])&MASK,(g+s[6])&MASK,(h+s[7])&MASK]

def wang_cascade(W0, W1, DW0=1):
    Wn = [W0, W1] + [0]*14
    DWs = [0]*16; DWs[0] = DW0
    # round 3: zero de3
    from_sha = lambda W, R: __sha_rounds(make_schedule(W), R)
    def __sha_rounds(W, R):
        a,b,c,d,e,f,g,h = IV
        sts = [[a,b,c,d,e,f,g,h]]
        for r in range(R):
            T1=(h+Sig1(e)+Ch(e,f,g)+K[r]+W[r])&MASK; T2=(Sig0(a)+Maj(a,b,c))&MASK
            h=g;g=f;f=e;e=(d+T1)&MASK;d=c;c=b;b=a;a=(T1+T2)&MASK
            sts.append([a,b,c,d,e,f,g,h])
        return sts
    Wf_tmp = [(Wn[i]+DWs[i])&MASK for i in range(16)]
    sn3 = __sha_rounds(make_schedule(Wn), 3)
    sf3 = __sha_rounds(make_schedule(Wf_tmp), 3)
    DWs[2] = (-(sf3[3][4] - sn3[3][4])) & MASK
    for step in range(13):
        wi = step+3; dt = step+4
        Wfc = [(Wn[i]+DWs[i])&MASK for i in range(16)]
        sn = __sha_rounds(make_schedule(Wn), dt)
        sf = __sha_rounds(make_schedule(Wfc), dt)
        DWs[wi] = (-(sf[dt][4] - sn[dt][4])) & MASK
    Wf = [(Wn[i]+DWs[i])&MASK for i in range(16)]
    # check f17
    sn = __sha_rounds(make_schedule(Wn), 17)
    sf = __sha_rounds(make_schedule(Wf), 17)
    de17 = (sf[17][4] - sn[17][4]) & MASK
    return Wn, Wf, de17

print("=" * 70)
print("АТАКА B4: 2-блочный MITM — структура финального состояния при f17=0")
print("=" * 70)

# ─────────────────────────────────────────────────────────────
# Тест 1: HW(δH[i]) для пар Wang-каскада с f17≠0 (обычные пары)
# ─────────────────────────────────────────────────────────────
print("\n[1] Распределение HW(δH[i]) = HW(H'[i] - H[i]) при Wang-каскаде (N=500):")
print(f"    (i=0..7 — 8 слов итогового хэша блока 1)")
print(f"{'слово':>6} | {'E[HW(δH)]':>11} | {'P(δH=0)':>10} | Теор. E=16")
print("-" * 45)

N = 500
dH_by_word = [[] for _ in range(8)]
for _ in range(N):
    W0 = random.randint(0, MASK); W1 = random.randint(0, MASK)
    Wn, Wf, de17 = wang_cascade(W0, W1)
    Hn = sha256_block(Wn)
    Hf = sha256_block(Wf)
    for i in range(8):
        dH_by_word[i].append((Hf[i] - Hn[i]) & MASK)

for i in range(8):
    lst = dH_by_word[i]
    e_hw = statistics.mean(hw(v) for v in lst)
    p0 = sum(1 for v in lst if v == 0) / N
    note = " ← структура!" if abs(e_hw - 16) > 0.5 else ""
    print(f"  H[{i}] | {e_hw:>11.4f} | {p0:>10.2e} |{note}")

# ─────────────────────────────────────────────────────────────
# Тест 2: Корреляция δH[i] и δH[j] (должна быть = 0)
# ─────────────────────────────────────────────────────────────
print("\n[2] Корреляция HW(δH[i]) ↔ HW(δH[j]) (независимость слов хэша):")
pairs_to_check = [(0,1),(0,4),(1,5),(2,6),(3,7),(4,5)]
for i, j in pairs_to_check:
    xi = [hw(v) for v in dH_by_word[i]]
    xj = [hw(v) for v in dH_by_word[j]]
    mx, mj = statistics.mean(xi), statistics.mean(xj)
    sx, sj = statistics.stdev(xi), statistics.stdev(xj)
    cov = sum((xi[k]-mx)*(xj[k]-mj) for k in range(N))/N
    rho = cov/(sx*sj) if sx*sj > 0 else 0
    note = " ← КОРРЕЛЯЦИЯ!" if abs(rho) > 0.05 else ""
    print(f"  corr(H[{i}], H[{j}]) = {rho:+.4f}{note}")

# ─────────────────────────────────────────────────────────────
# Тест 3: Разность финального хэша при de17=0 vs de17≠0
# ─────────────────────────────────────────────────────────────
print("\n[3] Сравнение распределения δH[0] при de17=0 (известные пары) vs de17≠0:")
KNOWN_PAIRS = [(0xe82222c7, 0x516cfb41), (0xd4254551, 0x679ea4de), (0xe73eb86e, 0xdfa1b7b0)]
print("  Пары с de17=0:")
for W0, W1 in KNOWN_PAIRS:
    Wn, Wf, de17 = wang_cascade(W0, W1)
    Hn = sha256_block(Wn)
    Hf = sha256_block(Wf)
    dH_hw = [hw((Hf[i]-Hn[i])&MASK) for i in range(8)]
    print(f"    W0=0x{W0:08x}: de17=0x{de17:08x}  HW(δH)=[{','.join(map(str,dH_hw))}]  sum={sum(dH_hw)}")

print("\n  Типичные пары (de17≠0), первые 5:")
cnt = 0
for _ in range(10000):
    W0 = random.randint(0, MASK); W1 = random.randint(0, MASK)
    Wn, Wf, de17 = wang_cascade(W0, W1)
    if de17 != 0 and cnt < 5:
        Hn = sha256_block(Wn)
        Hf = sha256_block(Wf)
        dH_hw = [hw((Hf[i]-Hn[i])&MASK) for i in range(8)]
        print(f"    W0=0x{W0:08x}: de17=0x{de17:08x}  HW(δH)=[{','.join(map(str,dH_hw))}]  sum={sum(dH_hw)}")
        cnt += 1
    if cnt >= 5: break

# ─────────────────────────────────────────────────────────────
# Тест 4: Минимальный суммарный вес δH (N=5000) — есть ли малые δH?
# ─────────────────────────────────────────────────────────────
print("\n[4] Распределение суммарного HW(δH) = sum_i HW(δH[i]) (N=5000 пар):")
total_hw = []
for i in range(N):
    total_hw.append(sum(hw(v) for v in [dH_by_word[j][i] for j in range(8)]))

# дополнительные 4500 пар
for _ in range(4500):
    W0 = random.randint(0, MASK); W1 = random.randint(0, MASK)
    Wn, Wf, _ = wang_cascade(W0, W1)
    Hn = sha256_block(Wn); Hf = sha256_block(Wf)
    total_hw.append(sum(hw((Hf[i]-Hn[i])&MASK) for i in range(8)))

print(f"  E[sum HW(δH)] = {statistics.mean(total_hw):.3f}  (теор. = 8×16 = 128)")
print(f"  std = {statistics.stdev(total_hw):.3f}")
print(f"  min = {min(total_hw)}, max = {max(total_hw)}")
print(f"  Гистограмма (bin=10):")
from collections import Counter
bins = Counter(v // 10 * 10 for v in total_hw)
for b in sorted(bins):
    bar = "#" * (bins[b] // 20)
    print(f"    [{b:3d}-{b+9:3d}]: {bins[b]:5d}  {bar}")

# ─────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("ИТОГ АТАКИ B4 (2-блочный MITM):")
mean_dH = statistics.mean(statistics.mean(hw(v) for v in dH_by_word[i]) for i in range(8))
print(f"  E[HW(δH)] = {mean_dH:.4f} для каждого слова (теор.=16)")
if abs(mean_dH - 16) > 0.5:
    print(f"  СТРУКТУРА НАЙДЕНА — E[HW(δH)] ≠ 16, отклонение {mean_dH-16:.2f} бит!")
    print(f"  СТАТУС: B4 — ПЕРСПЕКТИВНЫЙ")
else:
    print(f"  Нет структуры в δH — финальный хэш случайный")
    print(f"  MITM стоимость = 2^128 (birthday на 256 бит) — не улучшает 2^64")
    print(f"  СТАТУС: B4 — ОТРИЦАТЕЛЬНЫЙ")
print("=" * 70)
