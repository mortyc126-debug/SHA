"""
АТАКА A2: Поиск нейтральных битов для барьерной функции f17

Нейтральный бит b для f17=0:
  Если f17(W0, W1) = 0, то f17(W0 ^ 2^b, W1) = 0  (или f17(W0, W1 ^ 2^(b-32)) = 0)
  при полном пересчёте адаптивного каскада.

Аналогия: Biham-Chen (SHA-1, 2004) — 14 нейтральных бит → ускорение 2^14x.
Если m нейтральных бит → birthday стоимость 2^{32-m}.

Основа: T_NEUTRAL_STEP, T_CASCADE_17, T_BIRTHDAY_COST17
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

def compute_f17_f18(W0, W1, DW0=1):
    """Вычислить f17 = Da13+δW16 и f18 = Da14+δW17 через адаптивный каскад."""
    Wn = [W0, W1] + [0]*14
    DWs = [0]*16
    DWs[0] = DW0

    # ΔW2 → De3=0
    Wf_tmp = [(Wn[i]+DWs[i])&MASK for i in range(16)]
    sn3 = sha_rounds(make_schedule(Wn), 3)
    sf3 = sha_rounds(make_schedule(Wf_tmp), 3)
    DWs[2] = (-(sf3[3][4] - sn3[3][4])) & MASK

    # Каскад ΔW3..ΔW15
    for step in range(13):
        wi = step+3; dt = step+4
        Wfc = [(Wn[i]+DWs[i])&MASK for i in range(16)]
        sn = sha_rounds(make_schedule(Wn), dt)
        sf = sha_rounds(make_schedule(Wfc), dt)
        DWs[wi] = (-(sf[dt][4] - sn[dt][4])) & MASK

    Wf = [(Wn[i]+DWs[i])&MASK for i in range(16)]
    Wn_s = make_schedule(Wn)
    Wf_s = make_schedule(Wf)
    sn = sha_rounds(Wn_s, 15)
    sf = sha_rounds(Wf_s, 15)

    da13 = (sf[13][0] - sn[13][0]) & MASK
    da14 = (sf[14][0] - sn[14][0]) & MASK
    DW16 = (Wf_s[16] - Wn_s[16]) & MASK
    DW17 = (Wf_s[17] - Wn_s[17]) & MASK
    f17 = (da13 + DW16) & MASK
    f18 = (da14 + DW17) & MASK
    return f17, f18, da13, da14, DW16, DW17

# ─────────────────────────────────────────────────────────────
# Известные пары с De3..De17=0 (f17=0 верифицировано)
# ─────────────────────────────────────────────────────────────
KNOWN_PAIRS = [
    (0xe82222c7, 0x516cfb41, "П-15"),
    (0xd4254551, 0x679ea4de, "П-16 #1"),
    (0xe73eb86e, 0xdfa1b7b0, "П-16 #2"),
]

print("=" * 70)
print("АТАКА A2: Поиск нейтральных битов для f17=0")
print("=" * 70)

# ─────────────────────────────────────────────────────────────
# Шаг 0: Верификация исходных пар
# ─────────────────────────────────────────────────────────────
print("\n[0] Верификация исходных пар:")
t0 = time.time()
for W0, W1, name in KNOWN_PAIRS:
    f17, f18, da13, da14, dw16, dw17 = compute_f17_f18(W0, W1)
    ok = "✓" if f17 == 0 else "✗"
    print(f"  {name}: W0=0x{W0:08x} W1=0x{W1:08x}  f17=0x{f17:08x} {ok}  f18=0x{f18:08x}")
print(f"  Время верификации: {time.time()-t0:.2f}с")

# ─────────────────────────────────────────────────────────────
# Шаг 1: Тест нейтральных битов — W0 пространство (биты 0..31)
# ─────────────────────────────────────────────────────────────
print("\n[1] Тест нейтральных битов W0 (биты 0..31) для всех 3 пар:")
print(f"{'Бит':>4} | {'П-15':>12} {'П-16#1':>12} {'П-16#2':>12} | {'Нейт.':>6}")
print("-" * 55)

neutral_w0 = []
t0 = time.time()
for b in range(32):
    row = []
    neutral_count = 0
    for W0, W1, name in KNOWN_PAIRS:
        W0_mod = W0 ^ (1 << b)
        f17_new, _, _, _, _, _ = compute_f17_f18(W0_mod, W1)
        is_neutral = (f17_new == 0)
        if is_neutral:
            neutral_count += 1
        row.append((f17_new, is_neutral))
    neutral_in_all = all(r[1] for r in row)
    neutral_in_any = any(r[1] for r in row)
    vals = [f"0x{r[0]:08x}{'★' if r[1] else ' '}" for r in row]
    marker = " ← НЕЙТРАЛЬНЫЙ" if neutral_in_all else (" ← частично" if neutral_in_any else "")
    print(f"  {b:>2} | {vals[0]} {vals[1]} {vals[2]} |{marker}")
    if neutral_in_all:
        neutral_w0.append(b)

print(f"\nВремя W0-тестирования: {time.time()-t0:.2f}с")
print(f"Нейтральных W0-битов (во всех 3 парах): {len(neutral_w0)} → {neutral_w0}")

# ─────────────────────────────────────────────────────────────
# Шаг 2: Тест нейтральных битов — W1 пространство (биты 0..31)
# ─────────────────────────────────────────────────────────────
print("\n[2] Тест нейтральных битов W1 (биты 0..31) для всех 3 пар:")
print(f"{'Бит':>4} | {'П-15':>12} {'П-16#1':>12} {'П-16#2':>12} | {'Нейт.':>6}")
print("-" * 55)

neutral_w1 = []
t0 = time.time()
for b in range(32):
    row = []
    for W0, W1, name in KNOWN_PAIRS:
        W1_mod = W1 ^ (1 << b)
        f17_new, _, _, _, _, _ = compute_f17_f18(W0, W1_mod)
        is_neutral = (f17_new == 0)
        row.append((f17_new, is_neutral))
    neutral_in_all = all(r[1] for r in row)
    neutral_in_any = any(r[1] for r in row)
    vals = [f"0x{r[0]:08x}{'★' if r[1] else ' '}" for r in row]
    marker = " ← НЕЙТРАЛЬНЫЙ" if neutral_in_all else (" ← частично" if neutral_in_any else "")
    print(f"  {b:>2} | {vals[0]} {vals[1]} {vals[2]} |{marker}")
    if neutral_in_all:
        neutral_w1.append(b)

print(f"\nВремя W1-тестирования: {time.time()-t0:.2f}с")
print(f"Нейтральных W1-битов (во всех 3 парах): {len(neutral_w1)} → {neutral_w1}")

# ─────────────────────────────────────────────────────────────
# Шаг 3: «Мягкий» нейтральный тест — хотя бы в 1 паре
# ─────────────────────────────────────────────────────────────
print("\n[3] Нейтральные биты хотя бы в 1 из 3 пар (W0 и W1):")
soft_w0 = []
soft_w1 = []
for b in range(32):
    for W0, W1, name in KNOWN_PAIRS:
        W0_mod = W0 ^ (1 << b)
        f17_new, _, _, _, _, _ = compute_f17_f18(W0_mod, W1)
        if f17_new == 0:
            soft_w0.append((b, name))
            break
for b in range(32):
    for W0, W1, name in KNOWN_PAIRS:
        W1_mod = W1 ^ (1 << b)
        f17_new, _, _, _, _, _ = compute_f17_f18(W0, W1_mod)
        if f17_new == 0:
            soft_w1.append((b, name))
            break

print(f"  W0 биты (хотя бы 1 пара): {soft_w0}")
print(f"  W1 биты (хотя бы 1 пара): {soft_w1}")

# ─────────────────────────────────────────────────────────────
# Шаг 4: Анализ чувствительности — HW(f17_new) при флипе бита
# ─────────────────────────────────────────────────────────────
print("\n[4] Средний HW(f17) после флипа бита (для пары П-15):")
W0_ref, W1_ref, _ = KNOWN_PAIRS[0]
print(f"{'Бит':>4} | {'HW(f17_W0)':>12} {'HW(f17_W1)':>12} | Интерпретация")
print("-" * 55)
for b in range(0, 32, 4):  # каждые 4 бита для краткости
    f17_w0, _, _, _, _, _ = compute_f17_f18(W0_ref ^ (1<<b), W1_ref)
    f17_w1, _, _, _, _, _ = compute_f17_f18(W0_ref, W1_ref ^ (1<<b))
    interp = "← малое влияние!" if (hw(f17_w0) < 8 or hw(f17_w1) < 8) else ""
    print(f"  {b:>2} | {hw(f17_w0):>12} {hw(f17_w1):>12} | {interp}")

# ─────────────────────────────────────────────────────────────
# Итог
# ─────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("ИТОГ АТАКИ A2:")
total_neutral = len(neutral_w0) + len(neutral_w1)
if total_neutral > 0:
    print(f"  НАЙДЕНО {total_neutral} нейтральных бит!")
    print(f"  W0: {neutral_w0}")
    print(f"  W1: {neutral_w1}")
    print(f"  Ускорение birthday: 2^32 → 2^{32-total_neutral}")
    print(f"  СТАТУС: ПРОРЫВ — барьер снижается!")
else:
    print(f"  Строгих нейтральных бит (во всех 3 парах): 0")
    print(f"  Мягких нейтральных W0 (хотя бы 1 пара): {len(soft_w0)}")
    print(f"  Мягких нейтральных W1 (хотя бы 1 пара): {len(soft_w1)}")
    if soft_w0 or soft_w1:
        print(f"  СТАТУС: Частичные нейтральные биты — возможна пара-специфичная оптимизация")
    else:
        print(f"  СТАТУС: Нейтральных бит нет — f17 полностью чувствительна ко всем битам входа")
print("=" * 70)
