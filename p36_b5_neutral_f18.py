"""
АТАКА B5: Нейтральные биты для f18=0 при уже выполненном f17=0

После нахождения пары (W0*, W1*) с f17=0 (De3..De17=0):
  Ищем биты b в (W0,W1), флип которых сохраняет f17=0 И также даёт f18=0.

Если такой бит существует → стоимость 2 барьеров = 2^32 вместо 2^64!

Также тестируем:
  - Нейтральные биты для f18=0 при f17≠0 (независимо)
  - Парная структура (f17,f18) при флипе бита

Основа: T_2D_BIRTHDAY_NEGATIVE, T_DE18_DECOMPOSITION, T_STATE17
"""

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

def compute_f17_f18_full(W0, W1, DW0=1):
    """Полный каскад: возвращает (f17, f18, da13, da14, dw16, dw17)."""
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
    sn = sha_rounds(Wn_s, 15); sf = sha_rounds(Wf_s, 15)
    da13 = (sf[13][0] - sn[13][0]) & MASK
    da14 = (sf[14][0] - sn[14][0]) & MASK
    dw16 = (Wf_s[16] - Wn_s[16]) & MASK
    dw17 = (Wf_s[17] - Wn_s[17]) & MASK
    f17 = (da13 + dw16) & MASK
    f18 = (da14 + dw17) & MASK
    return f17, f18, da13, da14, dw16, dw17

KNOWN_PAIRS = [
    (0xe82222c7, 0x516cfb41, "П-15"),
    (0xd4254551, 0x679ea4de, "П-16 #1"),
    (0xe73eb86e, 0xdfa1b7b0, "П-16 #2"),
]

print("=" * 70)
print("АТАКА B5: Нейтральные биты для f18=0 (при f17=0 и независимо)")
print("=" * 70)

# ─────────────────────────────────────────────────────────────
# Шаг 1: Базовые f18 для известных пар
# ─────────────────────────────────────────────────────────────
print("\n[1] Базовые значения f17, f18 для известных пар:")
t0 = time.time()
for W0, W1, name in KNOWN_PAIRS:
    f17, f18, da13, da14, dw16, dw17 = compute_f17_f18_full(W0, W1)
    print(f"  {name}: f17=0x{f17:08x} ({'✓' if f17==0 else '✗'})  "
          f"f18=0x{f18:08x}  Da13=0x{da13:08x}  Da14=0x{da14:08x}")
print(f"  Время: {time.time()-t0:.2f}с")

# ─────────────────────────────────────────────────────────────
# Шаг 2: Нейтральные биты для f18=0 при сохранении f17=0
# ─────────────────────────────────────────────────────────────
print("\n[2] Поиск битов b: флип сохраняет f17=0 И даёт f18=0:")
print(f"{'Бит':>4} | {'Пара':>8} | {'f17_new':>12} | {'f18_new':>12} | Статус")
print("-" * 60)

double_neutral = []  # биты нейтральные для f17 И f18=0
t0 = time.time()
for W0, W1, name in KNOWN_PAIRS:
    for w_idx in range(2):   # w_idx=0 → W0, w_idx=1 → W1
        label = "W0" if w_idx == 0 else "W1"
        for b in range(32):
            W0_m = W0 ^ (1 << b) if w_idx == 0 else W0
            W1_m = W1 ^ (1 << b) if w_idx == 1 else W1
            f17n, f18n, _, _, _, _ = compute_f17_f18_full(W0_m, W1_m)
            if f17n == 0 and f18n == 0:
                double_neutral.append((b, w_idx, name))
                print(f"  {b:>2} | {name:>8}.{label} | 0x{f17n:08x}★ | 0x{f18n:08x}★ | ДВОЙНОЙ НЕЙТРАЛЬНЫЙ!")
            elif f17n == 0:
                print(f"  {b:>2} | {name:>8}.{label} | 0x{f17n:08x}★ | 0x{f18n:08x}  | f17 сохранён, f18≠0")

elapsed = time.time() - t0
print(f"Время поиска: {elapsed:.2f}с  (итого 192 вызовов каскада)")

# ─────────────────────────────────────────────────────────────
# Шаг 3: Нейтральные биты независимо для f18=0 (без f17=0)
# ─────────────────────────────────────────────────────────────
print("\n[3] Нейтральные биты для f18=0 независимо (независимо от f17):")
print("    (флип b: f18(W0^2^b, W1) = f18(W0, W1)?  — нет, ищем f18=0)")
print()
print(f"{'Бит':>4} | {'П-15 f18_new':>14} | {'П-16#1 f18_new':>14} | {'П-16#2 f18_new':>14} | {'f18=0?':>6}")
print("-" * 65)

neutral_f18_w0 = []
neutral_f18_w1 = []
for b in range(32):
    row_w0 = []
    for W0, W1, name in KNOWN_PAIRS:
        _, f18n, _, _, _, _ = compute_f17_f18_full(W0 ^ (1<<b), W1)
        row_w0.append(f18n)
    if any(v == 0 for v in row_w0):
        neutral_f18_w0.append(b)
        print(f"  {b:>2}W0| {f'0x{row_w0[0]:08x}':>14} | {f'0x{row_w0[1]:08x}':>14} | {f'0x{row_w0[2]:08x}':>14} | {'★ НОЛЬ' if any(v==0 for v in row_w0) else ''}")

for b in range(32):
    row_w1 = []
    for W0, W1, name in KNOWN_PAIRS:
        _, f18n, _, _, _, _ = compute_f17_f18_full(W0, W1 ^ (1<<b))
        row_w1.append(f18n)
    if any(v == 0 for v in row_w1):
        neutral_f18_w1.append(b)
        print(f"  {b:>2}W1| {f'0x{row_w1[0]:08x}':>14} | {f'0x{row_w1[1]:08x}':>14} | {f'0x{row_w1[2]:08x}':>14} | {'★ НОЛЬ' if any(v==0 for v in row_w1) else ''}")

# ─────────────────────────────────────────────────────────────
# Шаг 4: Анализ чувствительности f18 к битам — HW(Δf18) при флипе
# ─────────────────────────────────────────────────────────────
print("\n[4] Чувствительность f17 и f18 к флипу каждого бита (пара П-15):")
W0r, W1r, _ = KNOWN_PAIRS[0]
_, f18_base, _, _, _, _ = compute_f17_f18_full(W0r, W1r)
f17_base = 0  # знаем f17=0

print(f"  Базовые: f17=0x{f17_base:08x} f18=0x{f18_base:08x} HW(f18)={hw(f18_base)}")
print(f"\n  {'Бит':>3} | {'HW(Δf17_W0)':>12} {'HW(Δf18_W0)':>12} | {'HW(Δf17_W1)':>12} {'HW(Δf18_W1)':>12}")
print("  " + "-" * 58)
for b in range(0, 32, 2):
    f17_w0, f18_w0, _, _, _, _ = compute_f17_f18_full(W0r ^ (1<<b), W1r)
    f17_w1, f18_w1, _, _, _, _ = compute_f17_f18_full(W0r, W1r ^ (1<<b))
    df17_w0 = hw(f17_w0 ^ f17_base)
    df18_w0 = hw((f18_w0 - f18_base) & MASK)
    df17_w1 = hw(f17_w1 ^ f17_base)
    df18_w1 = hw((f18_w1 - f18_base) & MASK)
    low_w0 = " ←" if df18_w0 < 8 else ""
    low_w1 = " ←" if df18_w1 < 8 else ""
    print(f"  {b:>3} | {df17_w0:>12} {df18_w0:>12}{low_w0} | {df17_w1:>12} {df18_w1:>12}{low_w1}")

# ─────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("ИТОГ АТАКИ B5:")
if double_neutral:
    print(f"  НАЙДЕНО {len(double_neutral)} ДВОЙНЫХ нейтральных бит!")
    for b, widx, name in double_neutral:
        print(f"  Бит {b} в {'W0' if widx==0 else 'W1'} (пара {name})")
    print(f"  СТАТУС: ПРОРЫВ — 2 барьера за стоимость одного!")
elif neutral_f18_w0 or neutral_f18_w1:
    print(f"  Нейтральных для f18 (независимо): W0={neutral_f18_w0}, W1={neutral_f18_w1}")
    print(f"  СТАТУС: f18 чувствительна, но есть биты с нулём f18 при флипе")
else:
    print(f"  Двойных нейтральных бит: 0")
    print(f"  Независимых нейтральных для f18: 0")
    print(f"  СТАТУС: B5 ОТРИЦАТЕЛЬНЫЙ — f17 и f18 полностью независимы")
print("=" * 70)
