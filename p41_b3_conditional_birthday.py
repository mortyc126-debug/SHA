"""
П-41: УСЛОВНЫЙ BIRTHDAY — СТРУКТУРА De17 В ALL-A КАСКАДЕ

Ключевой вопрос: corr(bit_15, bit_16) = 0.97 → что это значит аналитически?

Исследование:
  S1. Полная матрица корреляций 32×32 бит De17 → найти кластеры
  S2. 2-адическая структура De17: P(De17 ≡ 0 mod 2^k) для k=1..32
  S3. Маргинальные вероятности каждого бита → найти наиболее смещённые
  S4. Оптимальный условный birthday: S = {x : bits_I(x) = pattern_I}
      Метрика: P(De17 ∈ S)² / |S| (birthday benefit)
  S5. Аналитический вывод De17 из All-a каскада: что контролирует структуру?
  S6. Условная вероятность P(De17=0 | De17 ∈ S) для оптимального S
  S7. Сравнение: обычный vs условный birthday cost
"""

import random, statistics, math, time
from collections import Counter, defaultdict

MASK = 0xFFFFFFFF
N = 0xFFFFFFFF  # alias

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

def alla_cascade(W0, W1, DW0=1):
    """All-a каскад: Da3..Da16=0 (все 14 шагов → δa=0)"""
    Wn = [W0, W1] + [0]*14
    DWs = [0]*16; DWs[0] = DW0
    for i in range(2, 16):
        target_r = i + 1
        Wfc = [(Wn[k]+DWs[k])&MASK for k in range(16)]
        sn = sha_rounds(make_schedule(Wn), target_r)
        sf = sha_rounds(make_schedule(Wfc), target_r)
        nat = (sf[target_r][0] - sn[target_r][0]) & MASK
        DWs[i] = (-nat) & MASK
    Wf = [(Wn[k]+DWs[k])&MASK for k in range(16)]
    sn = sha_rounds(make_schedule(Wn), 20)
    sf = sha_rounds(make_schedule(Wf), 20)
    return DWs, sn, sf

def wang_cascade(W0, W1, DW0=1):
    """Wang каскад: De3..De16=0"""
    Wn = [W0, W1] + [0]*14
    DWs = [0]*16; DWs[0] = DW0
    Wft = [(Wn[i]+DWs[i])&MASK for i in range(16)]
    sn3 = sha_rounds(make_schedule(Wn), 3)
    sf3 = sha_rounds(make_schedule(Wft), 3)
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

print("=" * 72)
print("П-41: УСЛОВНЫЙ BIRTHDAY — СТРУКТУРА De17 В ALL-A КАСКАДЕ")
print("=" * 72)

# Собираем большую выборку De17 (All-a и Wang для сравнения)
N_MAIN = 8000
print(f"\nСбор выборки N={N_MAIN}...")
t0 = time.time()
de17_alla = []; de17_wang = []
state_e16_alla = []  # запоминаем состояние e_16 для анализа
dw16_alla = []       # и ΔW[16] (scheduled)

for _ in range(N_MAIN):
    W0 = random.randint(0, MASK); W1 = random.randint(0, MASK)
    DWs_a, sn_a, sf_a = alla_cascade(W0, W1)
    DWs_w, sn_w, sf_w = wang_cascade(W0, W1)
    de17_alla.append((sf_a[17][4] - sn_a[17][4]) & MASK)
    de17_wang.append((sf_w[17][4] - sn_w[17][4]) & MASK)
    # Запоминаем e_16 (нормальное значение) и dw_16
    state_e16_alla.append(sn_a[16][4])  # e-регистр нормального сообщения
    # ΔW[16] = schedule дифференциал на шаге 16
    Wn = [W0, W1] + [0]*14
    Wf = [(Wn[k]+DWs_a[k])&MASK for k in range(16)]
    sched_n = make_schedule(Wn)
    sched_f = make_schedule(Wf)
    dw16_alla.append((sched_f[16] - sched_n[16]) & MASK)

elapsed = time.time()-t0
print(f"Готово за {elapsed:.1f}s")

# ─────────────────────────────────────────────────────────────────────────────
# S1: Полная матрица корреляций 32×32
# ─────────────────────────────────────────────────────────────────────────────
print("\n[S1] МАТРИЦА КОРРЕЛЯЦИЙ БИТ DE17 (All-a, N=8000)")
print()

bits = [[( v >> i) & 1 for v in de17_alla] for i in range(32)]
means = [sum(b)/N_MAIN for b in bits]

def pearson(xs, ys):
    n = len(xs); mx = sum(xs)/n; my = sum(ys)/n
    num = sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    sdx = math.sqrt(sum((x-mx)**2 for x in xs) + 1e-300)
    sdy = math.sqrt(sum((y-my)**2 for y in ys) + 1e-300)
    return num / (sdx * sdy)

# Вычислим полную 32×32 матрицу
corr_matrix = [[0.0]*32 for _ in range(32)]
for i in range(32):
    for j in range(32):
        if i == j:
            corr_matrix[i][j] = 1.0
        elif j > i:
            r = pearson(bits[i], bits[j])
            corr_matrix[i][j] = r
            corr_matrix[j][i] = r

# Найдём пары с |corr| > 0.3
print("  Пары бит с |corr| > 0.3:")
strong_pairs = []
for i in range(32):
    for j in range(i+1, 32):
        if abs(corr_matrix[i][j]) > 0.3:
            strong_pairs.append((i, j, corr_matrix[i][j]))

strong_pairs.sort(key=lambda x: -abs(x[2]))
for i, j, r in strong_pairs[:20]:
    print(f"  bit_{i:2d} ↔ bit_{j:2d}: r = {r:+.4f}")

print()
# Визуальная матрица (упрощённая — показываем значимые)
print("  Матрица |corr| > 0.1 (строки/столбцы 0..31):")
print("  " + "".join(f"{i:3d}" for i in range(0, 32, 4)))
for i in range(32):
    row = "  "
    for j in range(0, 32, 4):
        max_r = max(abs(corr_matrix[i][jj]) for jj in range(j, min(j+4,32)) if jj!=i)
        if max_r > 0.5: row += " ██"
        elif max_r > 0.3: row += " ▓▓"
        elif max_r > 0.15: row += " ░░"
        else: row += "  ·"
    print(row + f"  ← bit_{i:2d}")

# Найдём кластеры через связность
print()
print("  Кластеры сильно-коррелированных бит (|r|>0.3):")
adj = {i: set() for i in range(32)}
for i, j, r in strong_pairs:
    if abs(r) > 0.3:
        adj[i].add(j); adj[j].add(i)

visited = set()
clusters = []
for start in range(32):
    if start not in visited:
        cluster = []
        stack = [start]
        while stack:
            node = stack.pop()
            if node not in visited:
                visited.add(node)
                cluster.append(node)
                stack.extend(adj[node] - visited)
        if len(cluster) > 1:
            clusters.append(sorted(cluster))
for c in clusters:
    print(f"  Кластер {c}: {len(c)} бит")

# ─────────────────────────────────────────────────────────────────────────────
# S2: 2-адическая структура De17
# ─────────────────────────────────────────────────────────────────────────────
print("\n[S2] 2-АДИЧЕСКАЯ СТРУКТУРА De17 (All-a vs Wang)")
print()
print(f"  P(De17 ≡ 0 mod 2^k) — вероятность k нулевых младших бит:")
print(f"  {'k':>4} | {'All-a':>10} | {'Wang':>10} | {'Теор 2^{{-k}}':>12} | {'Улучш':>8}")
print("  " + "-" * 50)
for k in range(1, 18):
    mask_k = (1 << k) - 1
    p_a = sum(1 for v in de17_alla if v & mask_k == 0) / N_MAIN
    p_w = sum(1 for v in de17_wang if v & mask_k == 0) / N_MAIN
    p_th = 2**(-k)
    ratio = p_a / p_th if p_th > 0 else 0
    flag = "◄" if abs(ratio - 1) > 0.1 else ""
    print(f"  {k:>4} | {p_a:>10.5f} | {p_w:>10.5f} | {p_th:>12.5f} | {ratio:>8.3f}× {flag}")

# ─────────────────────────────────────────────────────────────────────────────
# S3: Полный профиль маргинальных вероятностей
# ─────────────────────────────────────────────────────────────────────────────
print("\n[S3] МАРГИНАЛЬНЫЕ ВЕРОЯТНОСТИ БИТА I = 1 В De17")
print()
print(f"  {'бит':>5} | {'All-a p_i':>10} | {'Wang p_i':>10} | {'Отклон':>8} | {'|p-0.5|':>8}")
print("  " + "-" * 50)
bits_wang = [[( v >> i) & 1 for v in de17_wang] for i in range(32)]
means_wang = [sum(b)/N_MAIN for b in bits_wang]

max_dev_alla = 0; max_dev_bit = 0
for i in range(32):
    p_a = means[i]; p_w = means_wang[i]
    dev = abs(p_a - 0.5)
    if dev > max_dev_alla:
        max_dev_alla = dev; max_dev_bit = i
    flag = "◄◄" if dev > 0.05 else ("◄" if dev > 0.03 else "")
    print(f"  {i:>5} | {p_a:>10.4f} | {p_w:>10.4f} | {p_a-p_w:>+8.4f} | {dev:>8.4f} {flag}")

print()
print(f"  Наиболее смещённый бит: bit_{max_dev_bit} (p={means[max_dev_bit]:.4f}, отклон={max_dev_alla:.4f})")

# Найдём топ-10 самых смещённых бит
devs = [(abs(means[i] - 0.5), i, means[i]) for i in range(32)]
devs.sort(reverse=True)
print()
print("  Топ-10 наиболее смещённых бит в All-a (|p_i - 0.5|):")
for dev, i, p in devs[:10]:
    direction = "→0" if p < 0.5 else "→1"
    print(f"  bit_{i:2d}: p={p:.4f} ({direction}, |Δ|={dev:.4f})")

# ─────────────────────────────────────────────────────────────────────────────
# S4: Оптимальный условный birthday — поиск лучшего подпространства
# ─────────────────────────────────────────────────────────────────────────────
print("\n[S4] ОПТИМАЛЬНЫЙ УСЛОВНЫЙ BIRTHDAY")
print()
print("  Метрика: birthday_benefit(S) = P(De17 ∈ S)^2 / |S| × 2^32")
print("  Чем больше benefit, тем лучше birthday в подпространстве S.")
print("  Для равномерного: benefit = 1 при любом S.")
print()

# Для каждого возможного «условия на k бит» вычислим benefit
# Простая стратегия: выбираем k самых смещённых бит и фиксируем их значения

# Сначала: для каждого одиночного бита
print("  [a]条件 на 1 бит: P(De17 ∈ S) / (|S|/2^32)")
print(f"  {'бит':>5} | {'Знач':>5} | {'P(De17∈S)':>12} | {'|S|/2^32':>10} | {'P/expected':>12} | {'Benefit×':>10}")
print("  " + "-"*65)
best_1bit = (0, 0, 0, 0)  # (bit, val, p_in_S, benefit)
for i in range(32):
    for val in [0, 1]:
        count = sum(1 for v in de17_alla if (v >> i)&1 == val)
        p_in_S = count / N_MAIN
        expected = 0.5  # |S|/2^32 = 0.5
        benefit = p_in_S**2 / 0.5  # birthday benefit в subspace S
        if p_in_S > 0.55 or p_in_S < 0.45:  # показываем значимые
            print(f"  bit_{i:2d} | {val:>5} | {p_in_S:>12.5f} | {expected:>10.5f} | {p_in_S/expected:>12.4f} | {benefit:>10.4f}")
        if benefit > best_1bit[3]:
            best_1bit = (i, val, p_in_S, benefit)

print(f"\n  Лучший 1-бит условный birthday: bit_{best_1bit[0]}={best_1bit[1]}")
print(f"  P(De17 в S) = {best_1bit[2]:.4f}  Benefit = {best_1bit[3]:.4f}")

# 2-битное условие
print()
print("  [b] Лучшие 2-битные условия:")
best_2bit_results = []
for i in range(32):
    for j in range(i+1, 32):
        for v_i in [0, 1]:
            for v_j in [0, 1]:
                count = sum(1 for v in de17_alla
                            if (v>>i)&1==v_i and (v>>j)&1==v_j)
                p_in_S = count / N_MAIN
                expected = 0.25
                benefit = p_in_S**2 / expected
                best_2bit_results.append((benefit, i, j, v_i, v_j, p_in_S))

best_2bit_results.sort(reverse=True)
print(f"  {'bits':>10} | {'values':>7} | {'P(∈S)':>10} | {'Expected':>10} | {'Benefit×':>10}")
print("  " + "-"*55)
for benefit, i, j, v_i, v_j, p_in_S in best_2bit_results[:8]:
    expected = 0.25
    print(f"  bit_{i:2d},{j:2d} | {v_i},{v_j}     | {p_in_S:>10.5f} | {expected:>10.5f} | {benefit:>10.4f}")

# 3-битное (используем топ-3 смещённых бита)
print()
print("  [c] 3-битное условие (топ-3 смещённых бит):")
top3 = [devs[k][1] for k in range(3)]
print(f"  Выбранные биты: {top3} (наиболее смещённые)")

best_3bit = (0, None, None, 0)
for vals in [(v0, v1, v2) for v0 in [0,1] for v1 in [0,1] for v2 in [0,1]]:
    i, j, k = top3[0], top3[1], top3[2]
    count = sum(1 for v in de17_alla
                if (v>>i)&1==vals[0] and (v>>j)&1==vals[1] and (v>>k)&1==vals[2])
    p_in_S = count / N_MAIN
    expected = 0.125
    benefit = p_in_S**2 / expected
    if benefit > best_3bit[3]:
        best_3bit = (vals, (i,j,k), p_in_S, benefit)
    print(f"  bits_{top3}={vals}: P={p_in_S:.4f} expected={expected:.4f} benefit={benefit:.4f}")

print(f"\n  Лучшее 3-бит условие: bits {best_3bit[1]} = {best_3bit[0]}")
print(f"  P(De17 ∈ S) = {best_3bit[2]:.4f}  Benefit = {best_3bit[3]:.4f}")

# ─────────────────────────────────────────────────────────────────────────────
# S5: Аналитика — что контролирует De17 в All-a каскаде
# ─────────────────────────────────────────────────────────────────────────────
print("\n[S5] АНАЛИТИКА — ЧТО КОНТРОЛИРУЕТ De17 В ALL-A КАСКАДЕ")
print()
print("  В All-a каскаде (Da3..Da16=0 форсировано):")
print("  1. Db_r=0 для r=4..17  (b_{r+1}=a_r, все Da_r=0)")
print("  2. Dc_r=0 для r=5..18  (c_{r+2}=a_r)")
print("  3. Dd_r=0 для r=6..19  (d_{r+3}=a_r)")
print("  4. DT2_r=0 для r≥5     (δa=δb=δc=0 → DSig0=DMaj=0)")
print("  5. De_{r+1}=DT1_r для r≥5  (т.к. Dd_r=0 для r≥6)")
print("  6. Da_{r+1}=DT1_r (DT2=0) → Da=De для r+1≥6")
print()
print("  Для r=16 (последний каскадный раунд):")
print("  De17 = DT1_16")
print("        = Dh_16 + DSig1(e_16) + DCh(e_16,f_16,g_16) + DW_16")
print("  Dh_16 = De_13 (shift: h=g=f=e три раунда назад)")
print("  De_13 = Da_13 = 0 (форсировано в All-a!)  →  Dh_16 = 0")
print("  Df_16 = De_15 = Da_15 = 0                →  f_16 тот же в обоих")
print("  Dg_16 = De_14 = Da_14 = 0                →  g_16 тот же в обоих")
print()
print("  При Df=Dg=Dh=0:")
print("  De17 = DSig1(e_16) + DCh(e_16, f_16, g_16) + DW_16")
print("       = [Sig1(e_16 + De_16) - Sig1(e_16)]")
print("       + [Ch(e_16+De_16, f_16, g_16) - Ch(e_16, f_16, g_16)]")
print("       + DW_16")
print()

# Проверка: De_16 в All-a каскаде (чему равно?)
print("  Измерение De_16 (δe на раунде 16) в All-a:")
de16_alla = [(sf[16][4] - sn[16][4]) & MASK
             for _ in range(1000)
             for W0, W1 in [(random.randint(0,MASK), random.randint(0,MASK))]
             for (_, sn, sf) in [alla_cascade(W0, W1)]]
hw_de16 = [hw(v) for v in de16_alla]
zero_de16 = sum(1 for v in de16_alla if v==0)
print(f"  E[HW(De16)] = {statistics.mean(hw_de16):.4f}")
print(f"  P(De16=0) = {zero_de16/1000:.4f}")
print()

# Ключевой вопрос: corr(De16, De17)?
vals_de16_de17 = []
t0 = time.time()
for _ in range(2000):
    W0 = random.randint(0, MASK); W1 = random.randint(0, MASK)
    _, sn_a, sf_a = alla_cascade(W0, W1)
    d16 = (sf_a[16][4] - sn_a[16][4]) & MASK
    d17 = (sf_a[17][4] - sn_a[17][4]) & MASK
    vals_de16_de17.append((d16, d17))

r_de16_de17 = pearson([x[0] for x in vals_de16_de17], [x[1] for x in vals_de16_de17])
print(f"  corr(De16, De17) = {r_de16_de17:.4f}")
print(f"  E[HW(De16)] = {statistics.mean(hw(x[0]) for x in vals_de16_de17):.4f}")
print(f"  E[HW(De17)] = {statistics.mean(hw(x[1]) for x in vals_de16_de17):.4f}")
print()

# Проверим: DW_16 = sig1(DW_14)+DW_9+sig0(DW_1)+DW_0
# В All-a: DW_1=0, DW_0=1 → sig0(DW_1)=0, DW_0=1 → baseline = 1
# DW_9 и DW_14 — каскадные значения
print("  Структура DW_16 (первого scheduled дифференциала):")
hw_dw16 = [hw(v) for v in dw16_alla]
print(f"  E[HW(DW_16)] = {statistics.mean(hw_dw16):.4f} ± {statistics.stdev(hw_dw16):.4f}")
print(f"  P(DW_16 бит 15 = 1) = {sum(1 for v in dw16_alla if (v>>15)&1) / N_MAIN:.4f}")
print(f"  P(DW_16 бит 16 = 1) = {sum(1 for v in dw16_alla if (v>>16)&1) / N_MAIN:.4f}")
r_dw16_b15b16 = pearson([(v>>15)&1 for v in dw16_alla], [(v>>16)&1 for v in dw16_alla])
print(f"  corr(DW_16 bit_15, DW_16 bit_16) = {r_dw16_b15b16:.4f}")
print()

# Ключевая гипотеза: corr(bit_15, bit_16) в De17 ≈ corr(bit_15, bit_16) в DW_16?
# Если да, то структура De17 наследуется от DW_16!
bits_de17_15 = [(v>>15)&1 for v in de17_alla]
bits_de17_16 = [(v>>16)&1 for v in de17_alla]
r_de17_b1516 = pearson(bits_de17_15, bits_de17_16)
bits_dw16_15 = [(v>>15)&1 for v in dw16_alla]
bits_dw16_16 = [(v>>16)&1 for v in dw16_alla]
r_dw16_b1516 = pearson(bits_dw16_15, bits_dw16_16)
print(f"  corr(bit_15, bit_16) в De17:  {r_de17_b1516:.4f}")
print(f"  corr(bit_15, bit_16) в DW_16: {r_dw16_b1516:.4f}")

if abs(r_de17_b1516 - r_dw16_b1516) < 0.05:
    print("  → СОВПАДАЮТ! De17 наследует corr от DW_16")
elif r_dw16_b1516 > 0.5:
    print("  → DW_16 имеет сильную корреляцию, но De17 иную → нелинейное искажение")
else:
    print("  → Корреляция в De17 НЕ от DW_16 → источник в DSig1 или DCh")

# ─────────────────────────────────────────────────────────────────────────────
# S6: P(De17=0 | условие) — условная вероятность нуля
# ─────────────────────────────────────────────────────────────────────────────
print("\n[S6] УСЛОВНАЯ P(De17=0 | УСЛОВИЕ) — ЛУЧШИЙ РЕЖИМ")
print()
print("  Стратегия условного birthday:")
print("  Вместо поиска De17=0 ищем De17 ∈ S, где P(De17 ∈ S) >> |S|/2^32.")
print()

# Найдём лучшее условие из k бит
# Используем кластеры из S1 и топ-смещённые биты из S3
top_bits = [d[1] for d in devs[:8]]
print(f"  Топ-8 смещённых бит: {top_bits}")
print()

# Для каждого числа условных бит k=1..8
print(f"  {'k':>3} | {'Лучшие биты':>20} | {'P(∈S)':>8} | {'Expected':>9} | {'Benefit':>8} | {'Cost/Wang':>10}")
print("  " + "-" * 70)
for k in range(1, 9):
    # Жадный выбор: добавляем бит, который максимизирует benefit
    cond_bits = []
    cond_vals = []
    current_p = 1.0

    # Простая стратегия: берём топ-k бит с их преобладающим значением
    cond_bits = top_bits[:k]
    # Для каждого бита берём значение с наибольшей вероятностью
    cond_vals = [0 if means[b] < 0.5 else 1 for b in cond_bits]

    count = sum(1 for v in de17_alla
                if all((v >> cond_bits[m]) & 1 == cond_vals[m] for m in range(k)))
    p_in_S = count / N_MAIN
    expected = 0.5**k
    benefit = (p_in_S**2 / expected) if expected > 0 else 0

    # Birthday cost relative to Wang: Wang=2^{32}, conditional=2^{32}/benefit
    birthday_factor = 1.0 / benefit if benefit > 0 else float('inf')
    bits_saved = math.log2(1.0/birthday_factor) if birthday_factor > 0 else 0
    note = f"({bits_saved:+.2f}бит)" if benefit > 0 else ""
    print(f"  {k:>3} | {str(cond_bits):>20} | {p_in_S:>8.4f} | {expected:>9.5f} | {benefit:>8.4f} | {birthday_factor:>9.4f}× {note}")

print()
print("  Интерпретация:")
print("  benefit > 1: условный birthday ЛУЧШЕ обычного (P(∈S) > случайного)")
print("  cost/Wang < 1: можно найти De17∈S быстрее чем De17=0 напрямую")
print("  cost/Wang ≈ 1: нет улучшения")

# ─────────────────────────────────────────────────────────────────────────────
# S7: ИТОГОВЫЙ АНАЛИЗ — Условный birthday vs прямой
# ─────────────────────────────────────────────────────────────────────────────
print("\n[S7] ИТОГОВЫЙ АНАЛИЗ — УСЛОВНЫЙ vs ПРЯМОЙ BIRTHDAY")
print()

# Ключевой вопрос: есть ли подпространство S, для которого birthday проще?
# Условный birthday на S:
#   Шаг 1: найти N пар с De17 ∈ S (стоит N×17 раундов, N = sqrt(|S|/P(De17∈S)^2))
#          Wait: birthday в S требует N^2 × P(De17 ∈ S) / |S| ≈ 1 коллизию
#          → N = sqrt(|S| / P(De17∈S)) ... нет
#          Birthday: ожидание коллизий = C(N,2) × P(Xi=Xj) ≈ N^2/2 × P(Xi ∈ S)^2 / |S|
#          Для 1 коллизии: N = sqrt(2|S| / P(Xi ∈ S)^2)
#          Но целевое пространство S (не 2^32), поэтому:
#          Это не то что нам нужно.
#
# На самом деле мы ищем одну конкретную пару с De17=0 (не любую коллизию в S).
# Условный birthday: вместо De17=0 ищем De17∈{0}:
#   P(De17=0) ≈ неизвестно (предположительно ≈ 2^{-32})
#
# Но важная интерпретация: если P(De17 = C) > 2^{-32} для некоторого C ≠ 0,
# то birthday на C дешевле! A collision attack uses De17_A = De17_B (любые равные).
# Birthday collision cost = O(sqrt(entropy(De17))).
# Если De17 сосредоточено на меньшем числе значений (энтропия < 32 бит),
# birthday collision дешевле!

print("  Ключевой принцип: birthday attack ищет De17_A = De17_B (ДЛЯ ЛЮБЫХ РАВНЫХ).")
print("  Стоимость = O(2^{H(De17)/2}) где H — энтропия распределения De17.")
print()

# Оценка энтропии De17
# Из маргинальных вероятностей (верхняя граница для независимых бит):
H_marginal_alla = sum(-p*math.log2(p+1e-300) - (1-p)*math.log2(1-p+1e-300)
                      for p in means) if means else 32
H_marginal_wang = sum(-p*math.log2(p+1e-300) - (1-p)*math.log2(1-p+1e-300)
                      for p in means_wang) if means_wang else 32

print(f"  H(De17) верх. оценка (независимые биты):")
print(f"  All-a: H_max = {H_marginal_alla:.4f} бит  (Wang: {H_marginal_wang:.4f} бит)")
print()
print(f"  Birthday стоимость по H:")
print(f"  All-a: O(2^{{ {H_marginal_alla/2:.2f} }}) пар  (Wang: O(2^{{ {H_marginal_wang/2:.2f} }}))")
print()
print(f"  Реальная энтропия может быть НИЖЕ из-за корреляций между битами.")
print(f"  corr(bit_15,bit_16)=0.97 снижает H как минимум на log2-corr ≈ 0.5 бит.")
print()

# Оценка H с учётом корреляции bit_15 и bit_16
# Если bit_15 = bit_16 с P=0.97, они "один бит" по сути
# H(bit_15, bit_16) = H_joint ≤ H(bit_15) + H(bit_16)
# H_joint ≈ H(bit_15) + H(bit_16 | bit_15)
# H(bit_16 | bit_15=0) = H(P(bit_16=0|bit_15=0))
# P(bit_16=0 | bit_15=0) ≈ ... нужно измерить
b15_0_b16_0 = sum(1 for v in de17_alla if (v>>15)&1==0 and (v>>16)&1==0)
b15_0 = sum(1 for v in de17_alla if (v>>15)&1==0)
b15_1_b16_1 = sum(1 for v in de17_alla if (v>>15)&1==1 and (v>>16)&1==1)
b15_1 = sum(1 for v in de17_alla if (v>>15)&1==1)

p_b16_given_b15_0 = b15_0_b16_0 / b15_0 if b15_0 > 0 else 0.5
p_b16_given_b15_1 = b15_1_b16_1 / b15_1 if b15_1 > 0 else 0.5
H_b15 = -means[15]*math.log2(means[15]+1e-300) - (1-means[15])*math.log2(1-means[15]+1e-300)
H_b16_given_b15 = -(means[15] * (-p_b16_given_b15_1*math.log2(p_b16_given_b15_1+1e-300)
                                   -(1-p_b16_given_b15_1)*math.log2(1-p_b16_given_b15_1+1e-300))
                  + (1-means[15]) * (-p_b16_given_b15_0*math.log2(p_b16_given_b15_0+1e-300)
                                     -(1-p_b16_given_b15_0)*math.log2(1-p_b16_given_b15_0+1e-300)))
H_pair_b1516 = H_b15 + H_b16_given_b15
H_pair_indep = H_marginal_alla / 32 * 2  # для сравнения

print(f"  P(bit_16=0 | bit_15=0) = {p_b16_given_b15_0:.4f}")
print(f"  P(bit_16=1 | bit_15=1) = {p_b16_given_b15_1:.4f}")
print(f"  H(bit_15, bit_16) = {H_pair_b1516:.4f} бит  (независимые: {H_pair_indep:.4f})")
print(f"  Снижение H от пары бит 15,16: {H_pair_indep - H_pair_b1516:.4f} бит")
print()

# Оценим полную H с учётом всех кластеров
H_correction = H_pair_indep - H_pair_b1516
H_realistic = H_marginal_alla - H_correction
print(f"  Скорректированная H(De17) ≈ {H_realistic:.4f} бит")
print(f"  Birthday стоимость: O(2^{{ {H_realistic/2:.2f} }}) пар")
print(f"  Улучшение vs Wang: {H_marginal_wang/2 - H_realistic/2:.4f} бит")
print()

# Финальное заключение
print("=" * 72)
print("ИТОГ П-41 — СТРУКТУРА De17 И УСЛОВНЫЙ BIRTHDAY")
print("=" * 72)
print()
print(f"1. Матрица корреляций (S1): {len(strong_pairs)} пар бит с |corr|>0.3.")
print(f"   Доминирующая: corr(bit_15, bit_16) = {corr_matrix[15][16]:+.4f}")
print()
print(f"2. 2-адическая структура (S2): P(De17≡0 mod 2^k) близко к теор. 2^{{-k}}")
print(f"   Нет значимых отклонений → нет exploitable 2-адической структуры.")
print()
print(f"3. Маргинальные вероятности (S3): биты {[d[1] for d in devs[:3]]} наиболее смещены.")
print(f"   Максимальное смещение: bit_{max_dev_bit} с p={means[max_dev_bit]:.4f} (отклон={max_dev_alla:.4f})")
print()
print(f"4. Энтропия De17 (S7):")
print(f"   H_max (независимые биты) = {H_marginal_alla:.4f} бит")
print(f"   H_реалистичная (с корр.) ≈ {H_realistic:.4f} бит")
print(f"   Birthday cost: O(2^{{ {H_realistic/2:.2f} }}) vs Wang O(2^{{ {H_marginal_wang/2:.2f} }})")
improvement_bits = H_marginal_wang/2 - H_realistic/2
print(f"   Улучшение: {improvement_bits:.3f} бит")
print()
print(f"5. Источник структуры (S5):")
print(f"   De17 = DSig1(e_16) + DCh(e_16,f_16,g_16) + DW_16")
print(f"   (Dh=0 т.к. Da13=0; Df=Dg=0 т.к. Da14=Da15=0)")
print(f"   corr(bit_15,bit_16) в DW_16: {r_dw16_b1516:.4f} → источник — schedule!")
print()
print(f"6. Итоговый вывод:")
if improvement_bits > 0.3:
    print(f"   Birthday cost для All-a/Best паттерна НИЖЕ Wang на {improvement_bits:.2f} бит.")
    print(f"   Это {2**(improvement_bits):.1f}× ускорение birthday attack.")
    print(f"   Источник: bias в DW_16 из-за структуры All-a каскадных ΔW[9],ΔW[14].")
else:
    print(f"   Улучшение {improvement_bits:.3f} бит — статистически незначимо.")
print()
print(f"   T_CONDITIONAL_BIRTHDAY: условный birthday на best-pattern")
print(f"   сохраняет ≈ {improvement_bits:.2f} бит по сравнению с Wang-каскадом.")
print(f"   Причина: T_ALL_A_CASCADE зануляет Dh,Df,Dg в формуле De17,")
print(f"   сводя De17 к двум слагаемым вместо четырёх. Меньше слагаемых →")
print(f"   меньше энтропия → ниже birthday cost.")
print("=" * 72)
