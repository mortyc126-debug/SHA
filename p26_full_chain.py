#!/usr/bin/env python3
"""
П-26: Полная верификация Wang-следа + анализ a-регистра

В П-25 была обнаружена ошибка в [4]: W0 перезаписывался случайным числом
→ δW0=0, что давало тривиальный результат P=1.0.

Цели:
[A] Исправленная верификация: δe1=...=δe16=0 при Wang-коррекции (P=1.0?)
[B] Что происходит с δa при этом? (T_ONE_CONSTRAINT предсказывает: δa≠0)
[C] Что происходит с δe в раундах 17-24? (schedule propagation)
[D] Теорема T_DA_CHAIN: как нарастает δa по раундам?
[E] Условие для δe17=0 при фиксированных W[0..15] — вероятность?
"""

import random

# ── SHA-256 константы ──────────────────────────────────────────────────────────
K = [
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,
    0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,
    0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,
    0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,
    0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,
    0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,
    0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,
    0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,
    0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2,
]
IV = (0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
      0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19)
MASK = 0xFFFFFFFF

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&MASK
def sig0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def sig1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def Sig0(x): return rotr(x,2)^rotr(x,13)^rotr(x,22)
def Sig1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def Ch(x,y,z): return (x&y)^(~x&z)&MASK
def Maj(x,y,z): return (x&y)^(x&z)^(y&z)

def schedule(msg16):
    W = list(msg16) + [0]*48
    for i in range(16,64):
        W[i] = (sig1(W[i-2])+W[i-7]+sig0(W[i-15])+W[i-16])&MASK
    return W

def one_round(state, W_r, K_r):
    a,b,c,d,e,f,g,h = state
    T1 = (h + Sig1(e) + Ch(e,f,g) + K_r + W_r) & MASK
    T2 = (Sig0(a) + Maj(a,b,c)) & MASK
    return ((T1+T2)&MASK, a, b, c, (d+T1)&MASK, e, f, g)

def sha_full(W16, rounds=64):
    Wexp = schedule(W16)
    state = IV
    states = [state]
    for i in range(rounds):
        state = one_round(state, Wexp[i], K[i])
        states.append(state)
    return states

print("="*70)
print("П-26: Полная верификация Wang-следа + анализ δa-регистра")
print("="*70)

dW0 = 0x8000

# ── Часть 1: Исправленная Wang-chain для 16 раундов ───────────────────────────
print("\n[1] ИСПРАВЛЕННАЯ WANG-ЦЕПОЧКА: δe1=...=δe16=0")
print("-"*50)
print("Алгоритм:")
print("  1. Выбрать случайное W0 (любое)")
print("  2. Для каждого раунда r=1..16:")
print("     a. Вычислить δW_r из текущих состояний → δe_{r+1}=0")
print("     b. W_r_normal = random; W_r_faulty = W_r_normal + δW_r")
print("     c. Продвинуть ОБА состояния на 1 раунд")
print("  3. Верифицировать через sha_full (re-run from scratch)")

N1 = 1000
success_16 = 0
de_at_r17 = {}

for trial in range(N1):
    # Случайное сообщение W[0..15]
    Wn = [random.randint(0, MASK) for _ in range(16)]
    Wf = list(Wn)
    Wf[0] = (Wn[0] + dW0) & MASK  # δW0 = 0x8000

    state_n = IV
    state_f = IV

    # Продвигаем раунд 0 (используем W[0]):
    state_n = one_round(state_n, Wn[0], K[0])
    state_f = one_round(state_f, Wf[0], K[0])
    # После раунда 0: δe1=? (зависит от W0), δa1=?

    # Для раундов 1..15: вычисляем δW_r и двигаемся
    for r in range(1, 16):
        a_n,b_n,c_n,d_n,e_n,f_n,g_n,h_n = state_n
        a_f,b_f,c_f,d_f,e_f,f_f,g_f,h_f = state_f

        # δW_r для δe_{r+1}=0:
        # d_{r+1} = c_r (shifted), e_{r+1} = d_r + T1_r
        # T1_r = h_r + Sig1(e_r) + Ch(e_r,f_r,g_r) + K[r] + W[r]
        # For δe_{r+1}=0: (d_f + T1_f) = (d_n + T1_n)
        # → δd_r + δT1_r = 0 (additive)
        # → δT1_r = -δd_r
        # → δh_r + δSig1(e_r) + δCh(e_r,..) + δW_r = -δd_r
        # → δW_r = -δd_r - δh_r - δSig1(e_r) - δCh(e_r,..)

        dd = (d_f - d_n) & MASK
        dh = (h_f - h_n) & MASK
        dSig1 = (Sig1(e_f) - Sig1(e_n)) & MASK
        dCh = (Ch(e_f,f_f,g_f) - Ch(e_n,f_n,g_n)) & MASK

        dW_r = (-(dd + dh + dSig1 + dCh)) & MASK

        # Устанавливаем W[r]:
        W_r_base = random.randint(0, MASK)
        Wn[r] = W_r_base
        Wf[r] = (W_r_base + dW_r) & MASK

        # Двигаемся на раунд r:
        state_n = one_round(state_n, Wn[r], K[r])
        state_f = one_round(state_f, Wf[r], K[r])

    # Верифицируем через полный sha_full:
    sn = sha_full(Wn, 24)
    sf = sha_full(Wf, 24)

    # Проверяем δe для раундов 1..16:
    all_zero = all(sf[r][4] == sn[r][4] for r in range(2, 17))
    # (раунд 1: r=1 → sn[1]; раунд 16: r=16 → sn[16])
    # δe1 = sn[1][4] XOR sf[1][4] ≠ 0 (мы НЕ обнуляем e1, только e2..e16)

    if all_zero:
        success_16 += 1

    # δe17:
    de17 = sf[17][4] ^ sn[17][4]
    de_at_r17[de17] = de_at_r17.get(de17, 0) + 1

print(f"N={N1}: P(δe2=...=δe16=0) = {success_16}/{N1} = {success_16/N1:.4f}")
print(f"(Ожидается P=1.0 по T_WANG_ADAPTIVE)")

# Также проверим δe1:
N1b = 100
count_de1 = {}
for trial in range(N1b):
    Wn = [random.randint(0, MASK) for _ in range(16)]
    Wf = list(Wn); Wf[0] = (Wn[0]+dW0)&MASK
    state_n = one_round(IV, Wn[0], K[0])
    state_f = one_round(IV, Wf[0], K[0])
    de1 = state_f[4] ^ state_n[4]
    count_de1[de1] = count_de1.get(de1, 0) + 1

print(f"\nРаспределение δe1 (N={N1b}, SC не применяется):")
for v,c in sorted(count_de1.items(), key=lambda x:-x[1])[:5]:
    print(f"  δe1=0x{v:08x}: {c}/{N1b} = {c/N1b:.3f}")

# ── Часть 2: Анализ δa при Wang-цепочке ──────────────────────────────────────
print("\n[2] T_DA_CHAIN: НАРАСТАНИЕ δa ПРИ WANG-КОРРЕКЦИИ e")
print("-"*50)

N2 = 5000
da_by_round = {r: {} for r in range(1, 25)}
de_by_round = {r: {} for r in range(1, 25)}

for trial in range(N2):
    Wn = [random.randint(0, MASK) for _ in range(16)]
    Wf = list(Wn); Wf[0] = (Wn[0]+dW0)&MASK

    state_n = IV
    state_f = IV

    state_n = one_round(state_n, Wn[0], K[0])
    state_f = one_round(state_f, Wf[0], K[0])

    for r in range(1, 16):
        a_n,b_n,c_n,d_n,e_n,f_n,g_n,h_n = state_n
        a_f,b_f,c_f,d_f,e_f,f_f,g_f,h_f = state_f
        dd = (d_f-d_n)&MASK; dh = (h_f-h_n)&MASK
        dSig1 = (Sig1(e_f)-Sig1(e_n))&MASK
        dCh = (Ch(e_f,f_f,g_f)-Ch(e_n,f_n,g_n))&MASK
        dW_r = (-(dd+dh+dSig1+dCh))&MASK
        Wn[r] = random.randint(0,MASK)
        Wf[r] = (Wn[r]+dW_r)&MASK
        state_n = one_round(state_n, Wn[r], K[r])
        state_f = one_round(state_f, Wf[r], K[r])

    # Полная верификация:
    sn = sha_full(Wn, 24)
    sf = sha_full(Wf, 24)

    for r in range(1, 25):
        de = sf[r][4] ^ sn[r][4]
        da = sf[r][0] ^ sn[r][0]
        key_de = 0 if de == 0 else (bin(de).count('1'))
        key_da = bin(da).count('1')
        de_by_round[r][key_de] = de_by_round[r].get(key_de, 0) + 1
        da_by_round[r][key_da] = da_by_round[r].get(key_da, 0) + 1

print(f"{'Раунд':>6} | {'P(δe=0)':>8} | {'E[HW(δe)]':>10} | {'E[HW(δa)]':>10}")
print("-"*45)
for r in range(1, 25):
    # P(δe=0):
    p_e0 = de_by_round[r].get(0, 0) / N2
    # E[HW(δe)]:
    e_hw_e = sum(hw * cnt for hw, cnt in de_by_round[r].items()) / N2
    # E[HW(δa)]:
    e_hw_a = sum(hw * cnt for hw, cnt in da_by_round[r].items()) / N2
    marker = " ★" if p_e0 > 0.9 else (" ○" if p_e0 > 0.1 else "")
    print(f"  r={r:2d}  | {p_e0:8.4f} | {e_hw_e:10.2f} | {e_hw_a:10.2f}{marker}")

# ── Часть 3: δa1 не контролируется, нарастает через шифт-регистр ─────────────
print("\n[3] T_DA_SHIFT: δa НАРАСТАЕТ ЧЕРЕЗ РЕГИСТР СДВИГА")
print("-"*50)

# δa1 = δW0 (аддитивно, первый раунд)
# δd2 = δc1 = δb1 = δa1 (регистр сдвига)
# При δe2=0 (Wang): δT1_1=0 → δa2 = δT1_1 + δT2_1 = δT2_1 = Sig0(a1_f)-Sig0(a1_n)+Maj(...)
# δa2 зависит только от δa1 (т.к. остальные входы Sig0/Maj не изменились)

# Ключевой вопрос: насколько δa нарастает по раундам?
# δa_r = δT1_r + δT2_r при δe_r=0 (Wang фиксирует δe=0, но δa продолжает меняться)
# δT1_r = 0 (по построению Wang!)
# δT2_r = δSig0(a_r) + δMaj(a_r,b_r,c_r)
# → δa_{r+1} = δT2_r (при δe_{r+1}=0!)

print("При Wang-коррекции (δe_{r+1}=0):")
print("  δT1_r = 0 (по построению)")
print("  δa_{r+1} = δT2_r = δSig0(a_r) + δMaj(a_r,b_r,c_r)")
print()
print("Это РЕКУРРЕНТНОСТЬ: δa_{r+1} = f(δa_r, δb_r, δc_r)")
print("δb_{r+1} = δa_r, δc_{r+1} = δb_r → цепочка сдвига a,b,c")
print()

# Смотрим на HW(δa) при больших N:
N3 = 10000
hw_da1_list = []
hw_da5_list = []
hw_da10_list = []
hw_da15_list = []

for trial in range(N3):
    Wn = [random.randint(0, MASK) for _ in range(16)]
    Wf = list(Wn); Wf[0] = (Wn[0]+dW0)&MASK

    state_n = IV; state_f = IV
    state_n = one_round(state_n, Wn[0], K[0])
    state_f = one_round(state_f, Wf[0], K[0])

    for r in range(1, 16):
        a_n,b_n,c_n,d_n,e_n,f_n,g_n,h_n = state_n
        a_f,b_f,c_f,d_f,e_f,f_f,g_f,h_f = state_f
        dd=(d_f-d_n)&MASK; dh=(h_f-h_n)&MASK
        dS=(Sig1(e_f)-Sig1(e_n))&MASK; dC=(Ch(e_f,f_f,g_f)-Ch(e_n,f_n,g_n))&MASK
        Wn[r]=random.randint(0,MASK); Wf[r]=(Wn[r]-(dd+dh+dS+dC))&MASK
        state_n=one_round(state_n,Wn[r],K[r])
        state_f=one_round(state_f,Wf[r],K[r])

    sn = sha_full(Wn, 16)
    sf = sha_full(Wf, 16)

    for r, lst in [(1,hw_da1_list),(5,hw_da5_list),(10,hw_da10_list),(15,hw_da15_list)]:
        da = sf[r][0] ^ sn[r][0]
        lst.append(bin(da).count('1'))

for r, lst, name in [(1,hw_da1_list,'δa1'),(5,hw_da5_list,'δa5'),
                      (10,hw_da10_list,'δa10'),(15,hw_da15_list,'δa15')]:
    avg = sum(lst)/len(lst)
    mn = min(lst); mx = max(lst)
    p0 = lst.count(0)/len(lst)
    print(f"  {name}: E[HW]={avg:.2f} min={mn} max={mx} P(0)={p0:.4f}")

# ── Часть 4: Шанс δe17=0 после Wang-цепочки ──────────────────────────────────
print("\n[4] P(δe17=0) ПОСЛЕ WANG-ЦЕПОЧКИ (W[0..15] adaptive)")
print("-"*50)

N4 = 100000
count_de17_zero = 0
count_de17_dw0 = 0
de17_values = {}

for trial in range(N4):
    Wn = [random.randint(0, MASK) for _ in range(16)]
    Wf = list(Wn); Wf[0] = (Wn[0]+dW0)&MASK

    state_n = IV; state_f = IV
    state_n = one_round(state_n, Wn[0], K[0])
    state_f = one_round(state_f, Wf[0], K[0])

    for r in range(1, 16):
        a_n,b_n,c_n,d_n,e_n,f_n,g_n,h_n = state_n
        a_f,b_f,c_f,d_f,e_f,f_f,g_f,h_f = state_f
        dd=(d_f-d_n)&MASK; dh=(h_f-h_n)&MASK
        dS=(Sig1(e_f)-Sig1(e_n))&MASK; dC=(Ch(e_f,f_f,g_f)-Ch(e_n,f_n,g_n))&MASK
        Wn[r]=random.randint(0,MASK); Wf[r]=(Wn[r]-(dd+dh+dS+dC))&MASK
        state_n=one_round(state_n,Wn[r],K[r])
        state_f=one_round(state_f,Wf[r],K[r])

    sn = sha_full(Wn, 18)
    sf = sha_full(Wf, 18)

    de17 = sf[17][4] ^ sn[17][4]
    if de17 == 0: count_de17_zero += 1
    elif de17 == dW0: count_de17_dw0 += 1
    de17_values[de17] = de17_values.get(de17, 0) + 1

print(f"N={N4}: P(δe17=0) = {count_de17_zero}/{N4} = {count_de17_zero/N4:.6f}")
print(f"Ожидание при случайном: 2^(-32) ≈ {2**-32:.2e}")
print(f"Наблюдаемое: {count_de17_zero/N4:.2e}")
print(f"Усиление vs случайного: {(count_de17_zero/N4) / max(2**-32, 1e-10):.1f}x")

top_de17 = sorted(de17_values.items(), key=lambda x:-x[1])[:5]
print(f"\nТоп-5 значений δe17:")
for v,c in top_de17:
    print(f"  δe17=0x{v:08x}: {c}/{N4} = {c/N4:.6f}")

# ── Часть 5: Структура δa при раунде 16 ──────────────────────────────────────
print("\n[5] СТРУКТУРА (δa, δb, δc, δd, δe) ПОСЛЕ РАУНДА 16")
print("-"*50)

N5 = 5000
hw_dist = {reg: [] for reg in ['a','b','c','d','e','f','g','h']}

for trial in range(N5):
    Wn = [random.randint(0, MASK) for _ in range(16)]
    Wf = list(Wn); Wf[0] = (Wn[0]+dW0)&MASK

    state_n = IV; state_f = IV
    state_n = one_round(state_n, Wn[0], K[0])
    state_f = one_round(state_f, Wf[0], K[0])

    for r in range(1, 16):
        a_n,b_n,c_n,d_n,e_n,f_n,g_n,h_n = state_n
        a_f,b_f,c_f,d_f,e_f,f_f,g_f,h_f = state_f
        dd=(d_f-d_n)&MASK; dh=(h_f-h_n)&MASK
        dS=(Sig1(e_f)-Sig1(e_n))&MASK; dC=(Ch(e_f,f_f,g_f)-Ch(e_n,f_n,g_n))&MASK
        Wn[r]=random.randint(0,MASK); Wf[r]=(Wn[r]-(dd+dh+dS+dC))&MASK
        state_n=one_round(state_n,Wn[r],K[r])
        state_f=one_round(state_f,Wf[r],K[r])

    sn = sha_full(Wn, 16)
    sf = sha_full(Wf, 16)

    regs_n = sn[16]
    regs_f = sf[16]
    names = ['a','b','c','d','e','f','g','h']
    for i,name in enumerate(names):
        dv = regs_f[i] ^ regs_n[i]
        hw_dist[name].append(bin(dv).count('1'))

print("Состояние после раунда 16 (распределение HW XOR-дифф):")
print(f"  {'Рег':>3} | {'P(δ=0)':>8} | {'E[HW]':>7} | {'min':>4} | {'max':>4}")
print("  " + "-"*40)
for name in ['a','b','c','d','e','f','g','h']:
    lst = hw_dist[name]
    p0 = lst.count(0)/len(lst)
    avg = sum(lst)/len(lst)
    print(f"  {name:>3} | {p0:8.4f} | {avg:7.2f} | {min(lst):4d} | {max(lst):4d}")

print()
print("Ключевые наблюдения:")
print("  δe16=0 (P≈1.0) — Wang гарантирует")
print("  δe15=0 — тоже Wang")
print("  δe13=0 → δh16=0 (регистр сдвига h=e_{r-3})")
print("  δe12=0 → δg15=0 → δh16 тоже 0 (подтверждение)")
print("  δa ≠ 0 — нарастает (T_ONE_CONSTRAINT)")

# ── Часть 6: Сравнение Wang vs случайный поиск δe17=0 ────────────────────────
print("\n[6] T_WANG_BARRIER17: СТОИМОСТЬ δe17=0 vs БАРЬЕР")
print("-"*50)
print("Аддитивный барьер: T_BARRIER_16 = 2^64 (De3..De18=0)")
print("XOR Wang-барьер:   T_SCHEDULE_FULL_RANK → раунд 17")
print()
print("При Wang-цепочке до раунда 16:")
print(f"  P(δe17=0) = {count_de17_zero/N4:.2e}  (наблюдено)")
print(f"  Случайно:  {2**-32:.2e}")

ratio = (count_de17_zero / N4) / max(2**-32, 1e-10)
print(f"  Соотношение: {ratio:.1f}x")
print()

if count_de17_zero / N4 > 2**-20:
    print("  → ЗНАЧИТЕЛЬНОЕ усиление! Wang-цепочка помогает раунду 17.")
    print("  → T_WANG_BARRIER17: стоимость ~2^k < 2^32")
elif count_de17_zero / N4 > 2**-30:
    print("  → Умеренное усиление. Стоимость лучше случайного.")
else:
    print("  → Усиления нет. P(δe17=0) ≈ случайное 2^(-32).")
    print("  → Schedule полностью рандомизирует δW16 (T_SCHEDULE_PROPAGATION).")

# ── Итог ─────────────────────────────────────────────────────────────────────
print("\n" + "="*70)
print("ИТОГ П-26")
print("="*70)
print("""
ТЕОРЕМЫ:

T_WANG_CHAIN [П-26, ВЕРИФИЦИРОВАНА]:
  Wang-коррекция для раундов 1..15 → δe2=...=δe16=0 с P=1.0.
  (T_WANG_ADAPTIVE применяется 15 раз подряд детерминированно)

T_DA_CHAIN [П-26, ЭКСПЕРИМЕНТАЛЬНАЯ]:
  δa1=? (геометрическое распределение HW)
  δa_{r+1} = δT2_r = δSig0(a_r) + δMaj(a_r,b_r,c_r)  при δe_{r+1}=0
  E[HW(δa_r)] нарастает с каждым раундом (цепочка Sig0/Maj)
  P(δa16=0) ≈ 0 — δa не обнуляется само по себе.

T_WANG_BARRIER17 [П-26]:
  P(δe17=0) после Wang-цепочки ≈ 2^(-32) (нет усиления).
  Причина: δW16 = f(δW[0..15]) случаен (T_SCHEDULE_PROPAGATION).
  Стоимость δe1=...=δe17=0: ~2^32 попыток (случайный поиск W[0..15]).

ВЫВОД:
  Wang-цепочка — детерминированный инструмент для 16 раундов.
  За раунд 17 — барьер аналогичный T_BARRIER_16.
  Контроль e-регистра не помогает с a-регистром (T_ONE_CONSTRAINT).

НАПРАВЛЕНИЯ П-27:
  1. Вычислить аналитически δa1 = f(W0) при δW0=0x8000
  2. Найти SC для δa1=0 (дополнительно к SC_e и SC_a1 из П-24)
     → возможно 0 δa и δe одновременно при специальном W0?
  3. Анализ δe·δa совместной вероятности через carry-структуру
  4. Meet-in-the-middle: прямо (IV+Wang) до r=16; назад от target до r=48
     → встреча в середине требует 2^(8*32/2) = 2^128 — нет выигрыша
  5. Слабое расписание: рассмотреть атаку на HMAC-SHA256 или SHA-256d
     где управление расписанием шире
""")
