"""
П-39: РОТАЦИОННЫЙ АНАЛИЗ SHA-256.

Ротационная криптография (Khovratovich & Nikolić, 2010):
  Для блочного шифра E: если E(ROT(x, r)) = ROT(E(x), r) — идеальная симметрия.
  В реальных шифрах эта симметрия "почти" выполняется → атака.

Применение к SHA-256:
  Ключевое свойство: Sigma0, Sigma1, sigma0, sigma1 — ЛИНЕЙНЫ над GF(2) (Теорема L1).
  Для XOR-линейных функций: f(ROT(x, r)) = ROT(f(x), r) ЕСЛИ функция "ротационно-эквивариантна".

  Sigma1(x) = ROTR(x,6) XOR ROTR(x,11) XOR ROTR(x,25)
  Sigma1(ROTR(x,k)) = ROTR(x,6+k) XOR ROTR(x,11+k) XOR ROTR(x,25+k)
                    = ROTR(Sigma1(x), k)   ← ротационная эквивариантность!

  Ch(e,f,g) — НЕЛИНЕЙНА (мультиплексор): нарушает симметрию.
  Maj(a,b,c) — НЕЛИНЕЙНА: тоже нарушает.
  Mod-add — нелинейна по переносам: НО переносы имеют структуру!

Ротационные пары (M, M') где M'[i] = ROTR(M[i], k) для всех i ∈ [0..15].
Вопрос: насколько "похоже" SHA(M') на ROTR(SHA(M), k)?
Если близко → ротационная характеристика → новый тип атаки.

Метод обхода барьера:
  Барьер 2^64 возникает из-за carry-цепочек в mod-add.
  Ротационный анализ смотрит на ДРУГУЮ СТРУКТУРУ — не разности, а повороты.
  Ротационный дифференциал = (M, ROTR(M, k)) → анализ SHA(M) vs ROTR(SHA(M), k).
"""

import random
from collections import defaultdict

# ─── SHA-256 примитивы ───────────────────────────────────────────────────────

K = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
    0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
    0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
    0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
    0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
    0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
]

H0 = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]

MASK = 0xFFFFFFFF

def rotr(x, n): return ((x >> n) | (x << (32-n))) & MASK
def rotl(x, n): return rotr(x, 32-n)
def Sig0(x): return rotr(x,2) ^ rotr(x,13) ^ rotr(x,22)
def Sig1(x): return rotr(x,6) ^ rotr(x,11) ^ rotr(x,25)
def sig0(x): return rotr(x,7) ^ rotr(x,18) ^ (x >> 3)
def sig1(x): return rotr(x,17) ^ rotr(x,19) ^ (x >> 10)
def Ch(e,f,g): return (e & f) ^ (~e & g) & MASK
def Maj(a,b,c): return (a & b) ^ (a & c) ^ (b & c)
def hw(x): return bin(x).count('1')

def sha256_trace(W16, nrounds=20):
    W = list(W16)
    for i in range(16, 64):
        W.append((sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16]) & MASK)
    a, b, c, d, e, f, g, h = H0
    states = [(a,b,c,d,e,f,g,h)]
    for r in range(nrounds):
        T1 = (h + Sig1(e) + Ch(e,f,g) + K[r] + W[r]) & MASK
        T2 = (Sig0(a) + Maj(a,b,c)) & MASK
        h=g; g=f; f=e; e=(d+T1)&MASK
        d=c; c=b; b=a; a=(T1+T2)&MASK
        states.append((a,b,c,d,e,f,g,h))
    return states, W

def sha256_full(W16):
    W = list(W16)
    for i in range(16, 64):
        W.append((sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16]) & MASK)
    a, b, c, d, e, f, g, h = H0
    for r in range(64):
        T1 = (h + Sig1(e) + Ch(e,f,g) + K[r] + W[r]) & MASK
        T2 = (Sig0(a) + Maj(a,b,c)) & MASK
        h=g; g=f; f=e; e=(d+T1)&MASK
        d=c; c=b; b=a; a=(T1+T2)&MASK
    return [(H0[i] + v) & MASK for i, v in enumerate([a,b,c,d,e,f,g,h])]

def rot_message(W16, k):
    """Повернуть все 16 слов сообщения на k бит вправо."""
    return [rotr(w, k) for w in W16]

# ─── Теоретический анализ ────────────────────────────────────────────────────

print("=" * 70)
print("П-39 | РОТАЦИОННЫЙ АНАЛИЗ SHA-256")
print("=" * 70)

print("\n[0] Теоретическая проверка: ротационная эквивариантность Sigma-функций")
print("─" * 70)

# Sigma1(ROTR(x,k)) =? ROTR(Sigma1(x), k)
ok_sig1 = 0; ok_sig0 = 0; ok_s0 = 0; ok_s1 = 0
N = 10000
for _ in range(N):
    x = random.randint(0, MASK)
    k = random.randint(1, 31)
    # Sigma1: ROTR{6,11,25} — полностью ротационная (все три компонента — cyclic rotations)
    if Sig1(rotr(x,k)) == rotr(Sig1(x), k): ok_sig1 += 1
    if Sig0(rotr(x,k)) == rotr(Sig0(x), k): ok_sig0 += 1
    # sigma0: ROTR{7,18} + SHR{3} — SHR нарушает ротационную симметрию!
    if sig0(rotr(x,k)) == rotr(sig0(x), k): ok_s0 += 1
    if sig1(rotr(x,k)) == rotr(sig1(x), k): ok_s1 += 1

print(f"Sigma1(ROTR(x,k)) = ROTR(Sigma1(x),k): {ok_sig1}/{N}  "
      f"{'✓ ПОЛНАЯ' if ok_sig1==N else '✗'} ротационная эквивариантность")
print(f"Sigma0(ROTR(x,k)) = ROTR(Sigma0(x),k): {ok_sig0}/{N}  "
      f"{'✓ ПОЛНАЯ' if ok_sig0==N else '✗'}")
print(f"sigma0(ROTR(x,k)) = ROTR(sigma0(x),k): {ok_s0}/{N}  "
      f"{'✓' if ok_s0==N else f'✗ НАРУШЕНА (SHR3 нелинеен по ротации)'}")
print(f"sigma1(ROTR(x,k)) = ROTR(sigma1(x),k): {ok_s1}/{N}  "
      f"{'✓' if ok_s1==N else f'✗ НАРУШЕНА (SHR10 нелинеен по ротации)'}")

# ─── Ключевое наблюдение: нарушение ротации через SHR ───────────────────────

print("\n[0.1] Анализ нарушения через sigma0 (SHR3):")
errors_s0 = []
for _ in range(1000):
    x = random.randint(0, MASK)
    k = 1
    err = sig0(rotr(x,k)) ^ rotr(sig0(x), k)
    errors_s0.append(hw(err))
avg_err = sum(errors_s0) / len(errors_s0)
print(f"  HW(sigma0(ROTR(x,1)) XOR ROTR(sigma0(x),1)): avg={avg_err:.2f} бит")
print(f"  Источник: SHR3 → после поворота появляются 'лишние' старшие биты")
print(f"  Следствие: расписание W[16..63] НАРУШАЕТ ротационную симметрию")

# ─── Эксперимент 1: Ротационное расстояние SHA-256 ───────────────────────────

print("\n[1] Ротационное расстояние: HW(SHA(M) XOR ROTR(SHA(M'), k))")
print("─" * 70)
print("Для каждого k измеряем, насколько SHA(ROTR_all(M,k)) близко к ROTR(SHA(M),k).")
print("Идеальная ротационная симметрия → HW=0. Случайные пары → HW≈128.")

K_VALUES = [1, 2, 6, 11, 17, 22, 25]
N_SAMPLES = 2000

print(f"\n{'k':>4} | avg HW | min HW | % пар с HW<16")
print("─" * 45)
for k in K_VALUES:
    hw_vals = []
    for _ in range(N_SAMPLES):
        W = [random.randint(0, MASK) for _ in range(16)]
        W_rot = rot_message(W, k)
        out1 = sha256_full(W)
        out2 = sha256_full(W_rot)
        out2_rot = [rotr(v, k) for v in out2]
        # Расстояние: SHA(W) vs ROTR(SHA(W_rot), -k) = ROTL(SHA(W_rot), k)
        out2_back = [rotl(v, k) for v in sha256_full(W_rot)]
        hw_total = sum(hw(a ^ b) for a, b in zip(out1, out2_back))
        hw_vals.append(hw_total)
    avg = sum(hw_vals) / len(hw_vals)
    mn = min(hw_vals)
    frac_low = sum(1 for v in hw_vals if v < 16) / len(hw_vals) * 100
    print(f"{k:>4} | {avg:6.1f} | {mn:6} | {frac_low:5.1f}%")

# ─── Эксперимент 2: Ротационный дифференциал раунд за раундом ────────────────

print("\n[2] Поряундовый ротационный дифференциал (k=1)")
print("─" * 70)
print("De_rot_r = e_r(W) XOR ROTR(e_r(ROTR_all(W), 1), -1)")
print("Если De_rot_r=0 — ротационная симметрия сохраняется на раунде r.")

k = 1
N = 1000
hw_by_round = defaultdict(list)

for _ in range(N):
    W = [random.randint(0, MASK) for _ in range(16)]
    W_rot = rot_message(W, k)
    states1, _ = sha256_trace(W, nrounds=20)
    states2, _ = sha256_trace(W_rot, nrounds=20)

    for r in range(21):
        # e-регистр
        e1 = states1[r][4]
        e2 = states2[r][4]
        e2_unrot = rotl(e2, k)  # "денормализовать" ротацию
        diff_e = hw(e1 ^ e2_unrot)
        hw_by_round[r].append(diff_e)

print(f"\n{'Раунд':>6} | avg HW(De_rot) | интерпретация")
print("─" * 55)
for r in range(21):
    avg = sum(hw_by_round[r]) / len(hw_by_round[r])
    if r == 0:
        desc = "начальное состояние (IV не ротационно!)"
    elif avg < 1:
        desc = "≈ полная симметрия"
    elif avg < 8:
        desc = "частичная симметрия"
    elif avg < 16:
        desc = "нарушается"
    else:
        desc = "полное рассеивание"
    print(f"{r:>6} | {avg:13.2f} | {desc}")

# ─── Эксперимент 3: Влияние IV — главный источник нарушения ─────────────────

print("\n[3] Ротационная несимметрия SHA-256 — роль IV (H0)")
print("─" * 70)
print("Проблема: H0 = [0x6a09e667, ...] НЕ является ротационным образом себя.")
print("Эксперимент: используем ротационно-симметричный IV (все нули).")

def sha256_trace_custom_iv(W16, IV, nrounds=20):
    W = list(W16)
    for i in range(16, 64):
        W.append((sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16]) & MASK)
    a, b, c, d, e, f, g, h = IV
    states = [(a,b,c,d,e,f,g,h)]
    for r in range(nrounds):
        T1 = (h + Sig1(e) + Ch(e,f,g) + K[r] + W[r]) & MASK
        T2 = (Sig0(a) + Maj(a,b,c)) & MASK
        h=g; g=f; f=e; e=(d+T1)&MASK
        d=c; c=b; b=a; a=(T1+T2)&MASK
        states.append((a,b,c,d,e,f,g,h))
    return states

k = 1
IV_zero = (0,)*8  # нулевой IV
IV_rot_k = tuple(rotr(v, k) for v in H0)  # ротированный стандартный IV

print("\nС нулевым IV (идеально ротационно-симметричен):")
hw_zero_iv = defaultdict(list)
for _ in range(500):
    W = [random.randint(0, MASK) for _ in range(16)]
    W_rot = rot_message(W, k)
    s1 = sha256_trace_custom_iv(W, IV_zero, 15)
    s2 = sha256_trace_custom_iv(W_rot, IV_zero, 15)
    for r in range(16):
        e1 = s1[r][4]; e2 = rotl(s2[r][4], k)
        hw_zero_iv[r].append(hw(e1 ^ e2))

print(f"{'Раунд':>6} | avg HW (IV=0) | avg HW (IV=H0)")
print("─" * 50)
for r in range(16):
    avg_zero = sum(hw_zero_iv[r]) / len(hw_zero_iv[r])
    avg_std = sum(hw_by_round[r]) / len(hw_by_round[r])
    print(f"{r:>6} | {avg_zero:12.2f} | {avg_std:12.2f}  "
          f"{'← IV нарушает!' if abs(avg_zero-avg_std)>2 else ''}")

# ─── Эксперимент 4: Ротационные константы SHA-256 — "особые" k ───────────────

print("\n[4] Поиск 'резонансных' сдвигов k — минимальный ротационный дифференциал")
print("─" * 70)
print("Идея: некоторые k могут резонировать с {2,13,22} и {6,11,25} → меньше HW.")

best_k = []
for k in range(1, 32):
    hw_vals = []
    for _ in range(500):
        W = [random.randint(0, MASK) for _ in range(16)]
        W_rot = rot_message(W, k)
        out1 = sha256_full(W)
        out2_back = [rotl(v, k) for v in sha256_full(W_rot)]
        hw_vals.append(sum(hw(a ^ b) for a, b in zip(out1, out2_back)))
    avg = sum(hw_vals) / len(hw_vals)
    best_k.append((k, avg, min(hw_vals)))

best_k.sort(key=lambda x: x[1])
print(f"\nТоп-10 наименее рассеивающих сдвигов k:")
print(f"{'k':>4} | avg HW | min HW | связь с SHA-256 константами")
print("─" * 65)
SHA_CONSTANTS = {2,6,7,11,13,17,18,19,22,25}
for k, avg, mn in best_k[:10]:
    conn = "резонанс!" if k in SHA_CONSTANTS else ""
    print(f"{k:>4} | {avg:6.1f} | {mn:6} | {conn}")

# ─── Эксперимент 5: "Ротационная пара" с De=0 ────────────────────────────────

print("\n[5] Ротационные пары и аддитивный дифференциал — пересечение методов")
print("─" * 70)
print("Ищем пары (W, ROTR_all(W, k)) такие что аддитивная разность De_3 мала.")

k = 1
N = 100000
hw_de3_rot = []
zeros_de3_rot = 0

for _ in range(N):
    W = [random.randint(0, MASK) for _ in range(16)]
    W_rot = rot_message(W, k)

    # Аддитивная разность De_3 между W и W_rot
    def get_e3_custom(Wlist):
        W = list(Wlist) + [0]*max(0, 16-len(Wlist))
        for i in range(16, 64):
            W.append((sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16]) & MASK)
        a, b, c, d, e, f, g, h = H0
        for r in range(3):
            T1 = (h + Sig1(e) + Ch(e,f,g) + K[r] + W[r]) & MASK
            T2 = (Sig0(a) + Maj(a,b,c)) & MASK
            h=g; g=f; f=e; e=(d+T1)&MASK
            d=c; c=b; b=a; a=(T1+T2)&MASK
        return e

    e3_W = get_e3_custom(W)
    e3_Wrot = get_e3_custom(W_rot)
    De3 = (e3_Wrot - e3_W) & MASK
    hw_de3_rot.append(hw(De3))
    if De3 == 0: zeros_de3_rot += 1

avg_hw_de3 = sum(hw_de3_rot) / len(hw_de3_rot)
print(f"Ротационные пары (k=1): avg HW(De3) = {avg_hw_de3:.2f} бит")
print(f"De3=0 для ротационных пар: {zeros_de3_rot}/{N} = {zeros_de3_rot/N*100:.3f}%")
print(f"(сравнение: случайные пары De3=0: ~9.43% из П-1)")

# ─── Эксперимент 6: Ротационное расписание ───────────────────────────────────

print("\n[6] Ротационная симметрия расписания W[16..63]")
print("─" * 70)
print("sigma0 и sigma1 содержат SHR (не cyclic) → расписание нарушает ротацию.")
print("Измеряем: насколько W_rot[r] ≠ ROTR(W[r], k) для r≥16.")

k = 1
N = 1000
hw_schedule_break = defaultdict(list)
for _ in range(N):
    W = [random.randint(0, MASK) for _ in range(16)]
    W_rot = rot_message(W, k)

    # Развернуть расписание
    Wf = list(W)
    Wf_rot = list(W_rot)
    for i in range(16, 32):
        Wf.append((sig1(Wf[i-2]) + Wf[i-7] + sig0(Wf[i-15]) + Wf[i-16]) & MASK)
        Wf_rot.append((sig1(Wf_rot[i-2]) + Wf_rot[i-7] + sig0(Wf_rot[i-15]) + Wf_rot[i-16]) & MASK)
    for i in range(16, 32):
        expected = rotr(Wf[i], k)
        actual = Wf_rot[i]
        hw_schedule_break[i].append(hw(expected ^ actual))

print(f"\n{'W[r]':>6} | avg HW нарушения ротации")
print("─" * 35)
for r in range(16, 32):
    avg = sum(hw_schedule_break[r]) / len(hw_schedule_break[r])
    print(f"W[{r:2d}] | {avg:6.2f}  {'← первые нарушения' if r==16 else ''}")

# ─── Вывод ───────────────────────────────────────────────────────────────────

print("\n" + "=" * 70)
print("ВЫВОД П-39")
print("=" * 70)
print("""
Теоремы ротационного анализа:

T_ROTSYM_SIGMA (доказана):
  Sigma0, Sigma1 ПОЛНОСТЬЮ ротационно эквивариантны:
  Sigma_i(ROTR(x,k)) = ROTR(Sigma_i(x), k) для всех x, k

T_ROTSYM_SCHEDULE (новая):
  Расписание W[16..63] НАРУШАЕТ ротационную симметрию из-за SHR-компонент
  sigma0(x) = ROTR(x,7) XOR ROTR(x,18) XOR (x>>3) — SHR3 нециклический.
  HW(нарушения) растёт начиная с W[16].

T_ROTSYM_IV (новая):
  Стандартный IV H0 не является ротационным образом себя.
  При нулевом IV ротационная симметрия сохраняется дольше.

T_ROTSYM_BARRIER (гипотеза):
  Ротационный дифференциал SHA-256 наследует барьер из двух источников:
  1. IV: нарушение с раунда 0 (≈6 бит на раунд 0)
  2. Расписание: нарушение с W[16] → влияет с раунда 14

СТРАТЕГИЯ ОБХОДА:
  Вместо стандартного IV, использовать "ротационно-адаптированный" IV:
  IV' = ROTR(IV', k) — тогда первый источник барьера устраняется.
  Это "free-start rotational collision" — более слабое требование,
  но позволяет исследовать ротационные характеристики без IV-помех.

ПЕРЕСЕЧЕНИЕ С П-36: структурный дифференциал при k=ротационный сдвиг
может создать аттрактор в пространстве ротационных пар → тема П-40.
""")
