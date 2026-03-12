"""
П-43: DW_9 — источник всего bias в De17.

Из П-42 установлено:
  DW_1 = 0 (не устанавливается каскадом)
  DW_14 = 0 (nat=0 после 12 предыдущих коррекций)
  sig1(0) = sig0(0) = 0

  → DW_16 = DW_9 + 1  (точно!)
  → De17 = DW_9 + 1   (по T_DE17_EQUALS_DW16)

Цель: Объяснить откуда P(LSB(DW_9)=1) = 0.7378 и бит-структуру DW_9.

DW_9 = -nat9 (mod 2^32), где nat9 = δa на раунде 10 без коррекции DW_9.
Коррекция: каскад при i=9 вычисляет nat9 и ставит DWs[9] = -nat9.
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

def alla_cascade_partial(W0, W1, stop_at, DW0=1):
    """All-a каскад до шага stop_at (включительно). Возвращает DWs."""
    Wn = [W0, W1] + [0]*14
    DWs = [0]*16; DWs[0] = DW0
    for i in range(2, min(stop_at+1, 16)):
        target_r = i + 1
        Wfc = [(Wn[k]+DWs[k])&MASK for k in range(16)]
        sn = sha_rounds(make_schedule(Wn), target_r)
        sf = sha_rounds(make_schedule(Wfc), target_r)
        nat = (sf[target_r][0] - sn[target_r][0]) & MASK
        DWs[i] = (-nat) & MASK
    return DWs, Wn

N = 10000

print("=" * 70)
print("П-43: DW_9 как единственный источник bias De17")
print("=" * 70)

print("\n[B1] Подтверждение: De17 = DW_9 + 1")

dw9_vals = []
de17_vals = []

for _ in range(N):
    W0 = random.randint(0, MASK)
    W1 = random.randint(0, MASK)
    DWs, Wn = alla_cascade_partial(W0, W1, 15)

    # De17 через schedule
    Wf = [(Wn[k]+DWs[k])&MASK for k in range(16)]
    sched_n = make_schedule(Wn)
    sched_f = make_schedule(Wf)
    dw16 = (sched_f[16] - sched_n[16]) & MASK
    de17 = dw16  # T_DE17_EQUALS_DW16

    dw9 = DWs[9]
    dw9_vals.append(dw9)
    de17_vals.append(de17)

    # Проверка De17 = DW_9 + 1
    expected = (dw9 + 1) & MASK
    if de17 != expected:
        print(f"  ОШИБКА: De17={de17:08x} != DW_9+1={expected:08x}")

match = sum(1 for i in range(N) if de17_vals[i] == (dw9_vals[i]+1)&MASK)
print(f"  De17 == DW_9 + 1: {match}/{N} ({match/N*100:.1f}%)")
print(f"  (Ожидается 100%)")

print("\n[B2] Структура DW_9: откуда P(LSB=1)=0.737?")
print()

# DW_9 = -nat9, где nat9 = δa_10 без коррекции DW_9
# Сбор nat9 напрямую
nat9_vals = []
for _ in range(N):
    W0 = random.randint(0, MASK)
    W1 = random.randint(0, MASK)
    # Коррекции DW_2..DW_8 (без DW_9)
    DWs_no9, Wn = alla_cascade_partial(W0, W1, 8)  # до i=8 включительно

    # nat9 = δa_10 с DWs[0..8], без DWs[9]
    Wfc = [(Wn[k]+DWs_no9[k])&MASK for k in range(16)]
    sn = sha_rounds(make_schedule(Wn), 10)
    sf = sha_rounds(make_schedule(Wfc), 10)
    nat9 = (sf[10][0] - sn[10][0]) & MASK
    nat9_vals.append(nat9)

# DW_9 = -nat9
dw9_from_nat = [(-v) & MASK for v in nat9_vals]

print(f"  nat9 (δa_10 перед коррекцией DW_9):")
print(f"    hw:  avg={st.mean(hw(v) for v in nat9_vals):.3f}  std={st.stdev(hw(v) for v in nat9_vals):.3f}")
print(f"    P(LSB(nat9)=0) = {sum(1 for v in nat9_vals if v&1==0)/N:.4f}")
print(f"    P(LSB(nat9)=1) = {sum(1 for v in nat9_vals if v&1==1)/N:.4f}")
print()
print(f"  DW_9 = -nat9 (mod 2^32):")
print(f"    P(LSB(DW_9)=0) = {sum(1 for v in dw9_from_nat if v&1==0)/N:.4f}")
print(f"    P(LSB(DW_9)=1) = {sum(1 for v in dw9_from_nat if v&1==1)/N:.4f}")
print()
print(f"  Ключевое наблюдение:")
print(f"  LSB(-x) = LSB(x) для нечётных x, LSB(-x) = 0 для чётных x.")
print(f"  Поэтому P(LSB(DW_9)=1) = P(LSB(nat9)=1) = P(nat9 нечётное)")
# Verify:
p_lsb_nat9_1 = sum(1 for v in nat9_vals if v&1==1)/N
p_lsb_dw9_1  = sum(1 for v in dw9_from_nat if v&1==1)/N
print(f"  P(LSB(nat9)=1) = {p_lsb_nat9_1:.4f}")
print(f"  P(LSB(DW_9)=1) = {p_lsb_dw9_1:.4f}  {'✓ совпадает' if abs(p_lsb_nat9_1-p_lsb_dw9_1)<0.01 else '✗ не совпадает'}")

print("\n[B3] Откуда P(LSB(nat9)=1) = 0.737?")
print("     nat9 = δa_10 = f(DW_0=1, DW_2..DW_8, W_0, W_1)")
print()
# nat9 аналитически раскладывается через раунды 2..9
# Ключевые вклады: DW_0=1 проходит через 9 раундов
# В линейном приближении: δa_r ≈ DW_{r-1} + cross-terms
# LSB(δa_r) = XOR-сумма LSB(DW_{r-9..r-1}) + нелинейность

# Проверим: зависит ли nat9 только от DW_0, DW_2..DW_8 и W_0, W_1?
# Да — это именно те коррекции которые применены

# Разложим: вклад от одного DW_0=1
print("  Тест изолированного вклада DW_0=1:")
nat9_seed_only = []
for _ in range(N):
    W0 = random.randint(0, MASK)
    W1 = random.randint(0, MASK)
    Wn = [W0, W1] + [0]*14
    # Только DW_0=1, остальные=0
    DWs_only0 = [0]*16; DWs_only0[0] = 1
    Wfc = [(Wn[k]+DWs_only0[k])&MASK for k in range(16)]
    sn = sha_rounds(make_schedule(Wn), 10)
    sf = sha_rounds(make_schedule(Wfc), 10)
    nat9_only = (sf[10][0] - sn[10][0]) & MASK
    nat9_seed_only.append(nat9_only)

print(f"    nat9 (только DW_0=1, остальные DW=0):")
print(f"    hw:  avg={st.mean(hw(v) for v in nat9_seed_only):.3f}")
print(f"    P(LSB=0) = {sum(1 for v in nat9_seed_only if v&1==0)/N:.4f}")
print(f"    P(LSB=1) = {sum(1 for v in nat9_seed_only if v&1==1)/N:.4f}")
print()

# Разница: реальный nat9 vs seed-only nat9
nat9_diff = [(nat9_vals[i] - nat9_seed_only[i]) & MASK for i in range(N)]
print(f"  nat9_real - nat9_seed_only (вклад DW_2..DW_8 коррекций):")
print(f"    hw avg = {st.mean(hw(v) for v in nat9_diff):.3f}")
print(f"    P(LSB_diff=0) = {sum(1 for v in nat9_diff if v&1==0)/N:.4f}")

print("\n[B4] Прямое тестирование LSB(nat9) через единственный вклад")
print("     LSB не зависит от переносов в предыдущих битах!")
print("     LSB(A+B+C...) = XOR(LSB(A), LSB(B), LSB(C)...)")
print()

# Разложим каждую коррекцию на вклад в LSB(δa_10)
# Сначала измерим вклад каждого DW_i (i=0,2..8) изолированно
contributions = {}
for target_k in [0, 2, 3, 4, 5, 6, 7, 8]:
    lsb_contrib = []
    for _ in range(2000):  # меньше для скорости
        W0 = random.randint(0, MASK)
        W1 = random.randint(0, MASK)
        Wn = [W0, W1] + [0]*14

        # Вычислим DW_k для нужного k (через partial cascade до k)
        DWs_partial, _ = alla_cascade_partial(W0, W1, target_k)
        # Изолированный вклад только DW_k
        DWs_single = [0]*16
        DWs_single[target_k] = DWs_partial[target_k]

        Wfc = [(Wn[k]+DWs_single[k])&MASK for k in range(16)]
        sn = sha_rounds(make_schedule(Wn), 10)
        sf = sha_rounds(make_schedule(Wfc), 10)
        nat = (sf[10][0] - sn[10][0]) & MASK
        lsb_contrib.append(nat & 1)

    p_lsb1 = sum(lsb_contrib)/len(lsb_contrib)
    contributions[target_k] = p_lsb1
    print(f"  P(LSB(вклад от DW_{target_k:2d} в δa_10) = 1) = {p_lsb1:.4f}")

print()
print("  Ключевые наблюдения:")
print("  - P=0.5 → вклад случаен (не объясняет bias)")
print("  - P≠0.5 → стабильный вклад в LSB(nat9)")
for k, p in contributions.items():
    if abs(p - 0.5) > 0.05:
        dev = p - 0.5
        print(f"    DW_{k}: P={p:.4f} (bias {dev:+.4f}) ← значимый вклад")

print("\n[B5] ФИНАЛЬНЫЙ РЕЗУЛЬТАТ")
print("=" * 60)
print()

# Итоговые P(бит_j DW_9 = 0)
print(f"  Распределение DW_9 (N={N}):")
print(f"  {'бит':>4} | {'P(=0)':>8} | Откл. | Примечание")
print(f"  {'-'*4}-+-{'-'*8}-+-{'-'*5}-+-{'-'*20}")
for j in range(12):
    p0 = sum(1 for v in dw9_vals if (v>>j)&1==0)/N
    dev = p0-0.5
    note = "***" if abs(dev)>0.1 else ""
    print(f"  {j:>4} | {p0:>8.4f} | {dev:+.3f} | {note}")

print()
print("  Распределение De17 = DW_9 + 1:")
print(f"  {'бит':>4} | {'P(=0)':>8} | Откл. | Примечание")
print(f"  {'-'*4}-+-{'-'*8}-+-{'-'*5}-+-{'-'*20}")
for j in range(12):
    p0 = sum(1 for v in de17_vals if (v>>j)&1==0)/N
    dev = p0-0.5
    note = "***" if abs(dev)>0.1 else ""
    print(f"  {j:>4} | {p0:>8.4f} | {dev:+.3f} | {note}")

print()
print("  Связь DW_9 ↔ De17 через +1:")
print("  bit_0(De17) = NOT(bit_0(DW_9)) когда DW_9 нечётное")
print(f"  P(bit_0(DW_9)=1) = {sum(1 for v in dw9_vals if v&1==1)/N:.4f}")
print(f"  P(bit_0(De17)=0) = {sum(1 for v in de17_vals if v&1==0)/N:.4f}")
print(f"  Проверка: равны? {'✓' if abs(sum(1 for v in dw9_vals if v&1==1)/N - sum(1 for v in de17_vals if v&1==0)/N) < 0.01 else '✗'}")

print("\n[B6] ИТОГОВАЯ ЦЕПОЧКА ПРИЧИННОСТИ")
print("=" * 60)
print()
print("  DW_0 = 1  (фиксировано в атаке)")
print("  DW_1 = 0  (loop starts at i=2, никогда не устанавливается)")
print("  DW_2..DW_8: каскадные коррекции для Da3..Da9=0")
print("  DW_9 = -nat9 (коррекция для Da10=0)")
print("         nat9 = δa_10 с только DW_0..DW_8 как inputs")
print("         → зависит от W_0, W_1 (random) И DW_0=1 (fixed)")
print()
print("  DW_10..DW_13: коррекции для Da11..Da14=0")
print("  DW_14 = 0  (nat14=0 из-за накопленных предыдущих коррекций)")
print("  DW_15: коррекция для Da16=0")
print()
print("  schedule формула: DW_16 = sig1(DW_14) + DW_9 + sig0(DW_1) + DW_0")
print("                           = sig1(0)    + DW_9 + sig0(0)    + 1")
print("                           = 0          + DW_9 + 0          + 1")
print("                           = DW_9 + 1")
print()
print("  T_DE17_EQUALS_DW16: De17 = DW_16 = DW_9 + 1")
print()
print(f"  Следствие: bias De17 = bias(DW_9 + 1)")
print(f"  P(LSB(DW_9)=1) = {sum(1 for v in dw9_vals if v&1)/N:.4f}")
print(f"  → P(LSB(De17)=0) = P(DW_9 нечётное) = {sum(1 for v in dw9_vals if v&1)/N:.4f}")
print()

# Ключевой вопрос: почему LSB(DW_9) = 1 с P=0.74?
# nat9 = δa_10. В модуле 2: LSB(δa_10) = XOR(LSB(DW_0), ...) через раунды
# DW_0=1 → его вклад в LSB(δa_10) определяет это
print("  Корень: P(LSB(nat9) нечётное) = P(LSB(DW_9) нечётное)")
p_nat9_odd = sum(1 for v in nat9_vals if v&1)/N
print(f"  P(nat9 нечётное) = {p_nat9_odd:.4f}")
print()
print("  Почему nat9 нечётное с P=0.74?")
print("  LSB(δa_10) = XOR(LSB(вклады от каждого DW_k)), k=0..8")
print(f"  В модуле 2: вклад DW_0=1 фиксирован (его LSB=1 всегда).")
print(f"  Вклады DW_2..DW_8 добавляют XOR со своими LSB.")
print(f"  Если большинство вкладов от DW_2..DW_8 имеют LSB=0 (чётные),")
print(f"  то LSB(nat9) = 1 (вклад от DW_0=1) → P ≈ 0.74.")

p_seed_only_odd = sum(1 for v in nat9_seed_only if v&1)/N
print(f"\n  nat9 при только DW_0=1 (остальные DW=0): P(нечётное)={p_seed_only_odd:.4f}")
print(f"  nat9 реальный (все DW_0..DW_8):           P(нечётное)={p_nat9_odd:.4f}")
print()
if abs(p_seed_only_odd - p_nat9_odd) < 0.05:
    print(f"  → DW_2..DW_8 корр. НЕ МЕНЯЮТ LSB(nat9)!")
    print(f"    Их суммарный XOR-вклад в LSB ≈ 0 (чётный)")
    print(f"    → ВЕСЬ bias P(LSB(nat9)=1)={p_nat9_odd:.4f} объясняется")
    print(f"       ОДНИМ слагаемым DW_0=1!")
else:
    print(f"  → DW_2..DW_8 меняют LSB(nat9): разница {abs(p_seed_only_odd-p_nat9_odd):.4f}")
