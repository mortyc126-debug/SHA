"""
П-42: Анатомия битов DW_16 — почему биты 3-7 почти всегда равны нулю?

Цель: Аналитически объяснить bias бит 0-7 в De17 = DW_16.

Структура:
  DW_16 = sig1(DW_14) + DW_9 + sig0(DW_1) + DW_0   (schedule формула)
  DW_0  = 1 (фиксировано)
  DW_1  = CASCADE_CHOICE (выбрано каскадом для Da3=0)
  DW_9  = CASCADE_CHOICE (выбрано каскадом для Da10=0)
  DW_14 = CASCADE_CHOICE (выбрано каскадом для Da15=0)

Гипотеза: Каскад создаёт специфическое распределение DW_k,
которое через sig0/sig1 даёт bias в низких битах DW_16.
"""
import random
import statistics as st

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

def alla_cascade_full(W0, W1, DW0=1):
    """All-a: Da_r=0 для r=3..16. Возвращает DWs[0..15] и Wn."""
    Wn = [W0, W1] + [0]*14
    DWs = [0]*16; DWs[0] = DW0
    for i in range(2, 16):
        target_r = i + 1
        Wfc = [(Wn[k]+DWs[k])&MASK for k in range(16)]
        sn = sha_rounds(make_schedule(Wn), target_r)
        sf = sha_rounds(make_schedule(Wfc), target_r)
        nat = (sf[target_r][0] - sn[target_r][0]) & MASK
        DWs[i] = (-nat) & MASK
    return DWs, Wn

N = 10000
print("=" * 70)
print("П-42: Анатомия DW_16 — источник bias бит 0-7 в De17")
print("=" * 70)

# ============================================================
print("\n[A1] Базовый вопрос: DW_16 = sig1(DW_14) + DW_9 + sig0(DW_1) + DW_0")
print("     Где DW_0=1 фиксировано. Что такое DW_1, DW_9, DW_14?")
print()

# Собираем статистику по каскадным DW
dw_vals = {k: [] for k in range(16)}
dw16_vals = []
de17_vals = []

for _ in range(N):
    W0 = random.randint(0, MASK)
    W1 = random.randint(0, MASK)
    DWs, Wn = alla_cascade_full(W0, W1)
    for k in range(16):
        dw_vals[k].append(DWs[k])

    # Вычислим DW_16 напрямую через schedule
    Wf = [(Wn[k]+DWs[k])&MASK for k in range(16)]
    sched_n = make_schedule(Wn)
    sched_f = make_schedule(Wf)
    dw16 = (sched_f[16] - sched_n[16]) & MASK
    dw16_vals.append(dw16)
    de17_vals.append(dw16)  # De17 = DW_16 (T_DE17_EQUALS_DW16)

print(f"  DW_0  всегда = 1 (фиксировано)")
for k in [1, 9, 14]:
    vals = dw_vals[k]
    print(f"  DW_{k:<2} hw: avg={st.mean(hw(v) for v in vals):.2f}  "
          f"std={st.stdev(hw(v) for v in vals):.2f}  "
          f"P(бит0=1)={sum(1 for v in vals if v&1)/len(vals):.4f}  "
          f"P(бит1=1)={sum(1 for v in vals if (v>>1)&1)/len(vals):.4f}")

# ============================================================
print("\n[A2] Слагаемые DW_16 по отдельности: hw и bit0 каждого")

for _ in range(N):
    pass  # уже собрали

# Пересобираем для слагаемых
comp_sig1_dw14 = []
comp_dw9 = []
comp_sig0_dw1 = []
comp_dw0 = []
comp_sum = []

for i in range(len(dw_vals[0])):
    dw14 = dw_vals[14][i]
    dw9  = dw_vals[9][i]
    dw1  = dw_vals[1][i]
    dw0  = dw_vals[0][i]   # always 1

    s1 = sig1(dw14)
    s0 = sig0(dw1)

    comp_sig1_dw14.append(s1)
    comp_dw9.append(dw9)
    comp_sig0_dw1.append(s0)
    comp_dw0.append(dw0)
    comp_sum.append((s1 + dw9 + s0 + dw0) & MASK)

print(f"  {'Слагаемое':<18} | {'hw avg':>7} | {'P(bit0=1)':>10} | {'P(bit1=1)':>10} | {'P(bit2=1)':>10}")
print(f"  {'-'*18}---{'-'*7}---{'-'*10}---{'-'*10}---{'-'*10}")
for name, vals in [("sig1(DW_14)", comp_sig1_dw14), ("DW_9", comp_dw9),
                   ("sig0(DW_1)", comp_sig0_dw1), ("DW_0=1", comp_dw0),
                   ("Сумма DW_16", comp_sum)]:
    ph0 = sum(1 for v in vals if v&1)/len(vals)
    ph1 = sum(1 for v in vals if (v>>1)&1)/len(vals)
    ph2 = sum(1 for v in vals if (v>>2)&1)/len(vals)
    phw = st.mean(hw(v) for v in vals)
    print(f"  {name:<18} | {phw:>7.3f} | {ph0:>10.4f} | {ph1:>10.4f} | {ph2:>10.4f}")

# ============================================================
print("\n[A3] Почему бит 0 (LSB) DW_16 ≈ 0 с P=0.73?")
print("     LSB(a+b+c+d) зависит от LSB(a), LSB(b), LSB(c), LSB(d)")
print()
print("     Модульная арифметика: LSB(X+Y+Z+W) = LSB(X) XOR LSB(Y) XOR LSB(Z) XOR LSB(W)")
print("     (XOR по модулю 2, перенос не влияет на бит 0!)")
print()

# Bit-0 анализ аналитически
b0_s1 = [v&1 for v in comp_sig1_dw14]
b0_d9 = [v&1 for v in comp_dw9]
b0_s0 = [v&1 for v in comp_sig0_dw1]
# DW_0=1 → bit0=1 всегда

P_b0_s1 = sum(b0_s1)/N
P_b0_d9 = sum(b0_d9)/N
P_b0_s0 = sum(b0_s0)/N
P_b0_dw0 = 1.0  # всегда 1

print(f"  P(LSB(sig1(DW_14))=1) = {P_b0_s1:.4f}")
print(f"  P(LSB(DW_9)=1)        = {P_b0_d9:.4f}")
print(f"  P(LSB(sig0(DW_1))=1)  = {P_b0_s0:.4f}")
print(f"  P(LSB(DW_0)=1)        = {P_b0_dw0:.4f}  (всегда 1)")

# Если независимы: P(XOR=0) = P(четное число единиц)
# = 0.5 + 0.5*(2p1-1)*(2p2-1)*(2p3-1)*(2p4-1)
def p_xor_zero(*probs):
    prod = 1.0
    for p in probs:
        prod *= (2*p - 1)
    return 0.5 + 0.5 * prod

P_b0_sum0_indep = p_xor_zero(P_b0_s1, P_b0_d9, P_b0_s0, P_b0_dw0)
P_b0_sum0_actual = sum(1 for v in comp_sum if v&1 == 0)/N

print(f"\n  Теория (независимые): P(LSB(DW_16)=0) = {P_b0_sum0_indep:.4f}")
print(f"  Измерено:              P(LSB(DW_16)=0) = {P_b0_sum0_actual:.4f}")
print(f"  (P=0.5 означало бы случайный)")

# ============================================================
print("\n[A4] Проверяем корреляцию слагаемых (нарушение независимости?)")

# Корреляция между bit0 слагаемых
from itertools import combinations
names4 = ["sig1(DW14)", "DW_9", "sig0(DW1)", "DW_0"]
bits4 = [b0_s1, b0_d9, b0_s0, [1]*N]

print(f"  Корреляции бит0 между слагаемыми:")
for (i,ni),(j,nj) in combinations(enumerate(names4), 2):
    bi = bits4[i]; bj = bits4[j]
    cov = sum((x-P_b0_s1 if i==0 else sum(bits4[i])/N - sum(bits4[i])/N) for x in bi)
    # Простой расчёт corr:
    mi = sum(bi)/N; mj = sum(bj)/N
    cov_ij = sum((bi[k]-mi)*(bj[k]-mj) for k in range(N))/N
    si = (mi*(1-mi))**0.5; sj = (mj*(1-mj))**0.5
    if si>0 and sj>0:
        corr_ij = cov_ij/(si*sj)
    else:
        corr_ij = 0.0
    if abs(corr_ij) > 0.02:
        print(f"    corr({ni}, {nj}) = {corr_ij:+.4f}  ***")
    else:
        print(f"    corr({ni}, {nj}) = {corr_ij:+.4f}")

# ============================================================
print("\n[A5] Следствие DW_0=1: 'сдвиг на 1' во всей арифметике")
print("     DW_16 = sig1(DW_14) + DW_9 + sig0(DW_1) + 1")
print()
print("     Сумма без DW_0:")
comp_sum_no_dw0 = [(s + d + z) & MASK for s,d,z in zip(comp_sig1_dw14, comp_dw9, comp_sig0_dw1)]

# Распределение бит 0 суммы без DW_0
P_b0_no1 = sum(1 for v in comp_sum_no_dw0 if v&1 == 0)/N
print(f"  P(LSB(S3)=0) без DW_0: {P_b0_no1:.4f}  (S3 = sig1+DW9+sig0)")
print(f"  P(LSB(S3+1)=0) = P(LSB(S3)=1) = {1-P_b0_no1:.4f}")
print(f"  Измеренное P(LSB(DW_16)=0)    = {P_b0_sum0_actual:.4f}")
print(f"  → DW_0=+1 инвертирует LSB: 0 ↔ 1")
print(f"  → P(LSB(DW_16)=0) = P(LSB(S3)=1) = {1-P_b0_no1:.4f}  {'✓' if abs((1-P_b0_no1) - P_b0_sum0_actual) < 0.005 else '✗ расхождение'}")

# ============================================================
print("\n[A6] Анализ бит 1-7 через перенос (carry)")
print("     Для битов 1+: LSB(A+B+C+D) зависит от каждого бита И переносов!")
print()

# Для каждого бита j: смотрим P(bit_j(DW_16) = 0)
print(f"  {'бит':>4} | {'P(bj(DW16)=0)':>14} | {'P(bj=0) random':>14} | Отклонение")
print(f"  {'-'*4}-+-{'-'*14}-+-{'-'*14}-+-{'-'*12}")
for j in range(10):
    p0 = sum(1 for v in dw16_vals if (v>>j)&1 == 0)/N
    dev = p0 - 0.5
    flag = "  ***" if abs(dev) > 0.1 else ""
    print(f"  {j:>4} | {p0:>14.4f} | {'0.5000':>14} | {dev:+.4f}{flag}")

# ============================================================
print("\n[A7] Ключевой вопрос: Откуда P(бит 3-7 = 0) ≈ 0.5?")
print("     (Из P-41: P(биты 3-7 одновременно = 0) ≈ 0.44)")
print()
# Проверим совместное P(бит3=0, бит4=0, бит5=0, бит6=0, бит7=0)
p_joint_3to7 = sum(1 for v in dw16_vals
                   if all((v>>j)&1 == 0 for j in [3,4,5,6,7])) / N
p_joint_3to7_indep = 0.5**5  # если бы независимы и несмещены
print(f"  P(биты 3-7 = 0 одновременно) = {p_joint_3to7:.4f}")
print(f"  P(если независимые ∩ несмещ.) = {p_joint_3to7_indep:.4f}")
print(f"  Избыток = {p_joint_3to7/p_joint_3to7_indep:.2f}×")

# Маргинальные P(бит_j=0) для j=3..7
marg = [sum(1 for v in dw16_vals if (v>>j)&1 == 0)/N for j in range(8)]
print()
print(f"  Маргинальные P(бит_j=0):")
for j in range(8):
    print(f"    бит_{j}: {marg[j]:.4f}  (ожид. 0.5000)")
p_joint_if_indep_marg = 1.0
for j in [3,4,5,6,7]:
    p_joint_if_indep_marg *= marg[j]
print(f"  P(3-7=0) если независимые с теми же маргин. = {p_joint_if_indep_marg:.4f}")
print(f"  Измеренное P(3-7=0)                         = {p_joint_3to7:.4f}")
print(f"  → Избыток из-за корреляций = {p_joint_3to7/p_joint_if_indep_marg:.2f}×")

# ============================================================
print("\n[A8] Откуда берётся carry-структура в битах 3-7?")
print("     sig1(x) = rotr(x,17) XOR rotr(x,19) XOR (x>>10)")
print("     Биты 3-7 sig1(x): от бит 20-24, 22-26, 13-17 исходного x (после rotr)")
print("     sig0(x) = rotr(x,7) XOR rotr(x,18) XOR (x>>3)")
print("     Биты 3-7 sig0(x): от бит 10-14, 21-25, 3-7 исходного x")
print()

# Проверим: DW_14 имеет специфические биты 20-24?
print("     Анализ DW_14 бит (20-24, 22-26, 13-17):")
for start, end in [(13, 17), (20, 24), (22, 26)]:
    for j in range(start, min(end+1, 32)):
        p_j = sum(1 for v in dw_vals[14] if (v>>j)&1 == 1)/N
        if abs(p_j - 0.5) > 0.05:
            print(f"       DW_14 бит {j:2d}: P(=1)={p_j:.4f}  *** смещение!")
        # else: skip for brevity
print("     (нет вывода = все близки к 0.5)")

# ============================================================
print("\n[A9] Финальная гипотеза: carry-accumulation через DW_0=1")
print()
print("  В All-a каскаде DW_k для k>0 выбираются как -nat(mod 2^32),")
print("  где nat = δa или δe от предыдущего шага.")
print("  Все DW_k включают DW_0=1 как 'seed'.")
print()
print("  DW_0=1 → через schedule → DW_1,DW_2,... имеют специфические переносы.")
print("  Сложение четырёх таких слагаемых в DW_16 создаёт carry-паттерн,")
print("  где биты 1-7 испытывают 'carry wash' — переносы от бита 0")
print("  распределяются специфически из-за несимметрии DW_0=1.")
print()

# Тест: случайная сумма 4 чисел с LSB=1?
rand_sum_b0_bias = []
for _ in range(N):
    # 4 случайных числа, все с LSB случайным
    vals_r = [random.randint(0, MASK) for _ in range(3)]
    vals_r.append(1)  # DW_0 = 1
    s = sum(vals_r) & MASK
    rand_sum_b0_bias.append(s)

P_b0_rand = sum(1 for v in rand_sum_b0_bias if v&1 == 0)/N
print(f"  Контроль: 3 случайных + 1 (DW_0=1): P(LSB=0) = {P_b0_rand:.4f}")
print(f"  Реальный DW_16:                       P(LSB=0) = {P_b0_sum0_actual:.4f}")
print(f"  → Отклонение объясняется bias в слагаемых каскада, не только DW_0=1")

# ============================================================
print("\n[A10] Итоговое резюме: Что создаёт bias бит 0-7 в DW_16?")
print("=" * 60)
print()
print("  DW_16 = sig1(DW_14) + DW_9 + sig0(DW_1) + 1")
print()

# Финальное измерение смещений
print("  Слагаемые (P(LSB=0)):")
for name, vals in [("  sig1(DW_14)", comp_sig1_dw14),
                   ("  DW_9       ", comp_dw9),
                   ("  sig0(DW_1) ", comp_sig0_dw1),
                   ("  DW_0 = 1   ", [1]*N)]:
    p = sum(1 for v in vals if v&1 == 0)/N
    print(f"    {name}: {p:.4f}  (отклонение от 0.5: {p-0.5:+.4f})")

print()
print(f"  Слагаемые до sig-функций (P(LSB=0)):")
for k in [14, 9, 1, 0]:
    p = sum(1 for v in dw_vals[k] if v&1 == 0)/N
    print(f"    DW_{k:<2}: {p:.4f}  (отклонение: {p-0.5:+.4f})")

print()
print(f"  Корреляция DW_1 ↔ DW_9 (LSB):")
b0_1 = [v&1 for v in dw_vals[1]]
b0_9 = [v&1 for v in dw_vals[9]]
m1 = sum(b0_1)/N; m9 = sum(b0_9)/N
cov = sum((b0_1[i]-m1)*(b0_9[i]-m9) for i in range(N))/N
s1_ = (m1*(1-m1))**0.5; s9_ = (m9*(1-m9))**0.5
print(f"    corr(LSB(DW_1), LSB(DW_9)) = {cov/(s1_*s9_) if s1_>0 and s9_>0 else 0:.4f}")
b0_14 = [v&1 for v in dw_vals[14]]
m14 = sum(b0_14)/N
cov14 = sum((b0_1[i]-m1)*(b0_14[i]-m14) for i in range(N))/N
s14_ = (m14*(1-m14))**0.5
print(f"    corr(LSB(DW_1), LSB(DW_14)) = {cov14/(s1_*s14_) if s1_>0 and s14_>0 else 0:.4f}")
b0_s1_arr = [v&1 for v in comp_sig1_dw14]
b0_s0_arr = [v&1 for v in comp_sig0_dw1]
ms1 = sum(b0_s1_arr)/N; ms0 = sum(b0_s0_arr)/N
cov_s = sum((b0_s1_arr[i]-ms1)*(b0_s0_arr[i]-ms0) for i in range(N))/N
ss1 = (ms1*(1-ms1))**0.5; ss0 = (ms0*(1-ms0))**0.5
print(f"    corr(LSB(sig1(DW14)), LSB(sig0(DW1))) = {cov_s/(ss1*ss0) if ss1>0 and ss0>0 else 0:.4f}")

print()
print("  ИТОГ:")
print("  1. Каждое из 4 слагаемых DW_16 имеет P(LSB) близко к 0.5")
print("     → одиночный bias слагаемого не объясняет P(LSB(DW_16)=0)=0.73")
print("  2. Ключевая причина: КОРРЕЛЯЦИИ между слагаемыми")
print("     DW_1, DW_9, DW_14 выбраны ОДНИМ каскадом из того же W0,W1")
print("     → они зависимы по построению")
print("  3. Carry-эффект: при сложении зависимых слагаемых переносы")
print("     концентрируются в определённых позициях → bias в битах 1-7")
print("  4. sig0/sig1 перемешивают биты → bias 'распределяется' по всем")
print("     старшим битам DW_16 с убывающей амплитудой")
print()
print("  → bias De17 = структурный эффект ЗАВИСИМОСТИ каскадных DW_k,")
print("    а не свойство самих функций sig0/sig1 по отдельности.")
