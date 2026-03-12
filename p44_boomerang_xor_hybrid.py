"""
П-44A: БУМЕРАНГ — XOR-ГИБРИДНЫЙ ПОДХОД

Проблема с T_BOOMERANG_INFEASIBLE (П-29):
  Прежний анализ рассматривал бумеранг в ADD-дифференциалах.
  HW нижней ADD-характеристики ≈ 64 → C ≈ 2^{96} (провал).

Новый угол: XOR-дифференциалы + сплит на раунде 8

Теория:
  Расписание SHA-256 ЛИНЕЙНО над GF(2) (L1):
    ΔW[i] = sig1(ΔW[i-2]) XOR ΔW[i-7] XOR sig0(ΔW[i-15]) XOR ΔW[i-16]
  Значит XOR-дифференциалы в W распространяются ТОЧНО (без вероятностных потерь).

  ВЕРХНЯЯ ЧАСТЬ (раунды 1-8):
    ΔW → Δstate_8 с вероятностью p₈ (через нелинейность Ch/Maj)

  НИЖНЯЯ ЧАСТЬ (раунды 9-17):
    Δstate_8 → Δe17=0 с вероятностью q₁₇ (снизу вверх)

  КЛАССИЧЕСКАЯ СТОИМОСТЬ: O(1/(p₈·q₁₇))
  БУМЕРАНГ (quartet): O(4/√(p₈·q₁₇)) — квадратное ускорение
    При p₈ = q₁₇ = 2⁻¹⁶: классика=2³², бумеранг=2¹⁷

  КЛЮЧЕВОЙ ВОПРОС: каковы реальные p₈ и q₁₇ для XOR-характеристики?

Тест 1: P(Δe_r=0) для каждого раунда r=1..17 при фиксированном XOR ΔW₀=1
Тест 2: Сплит-анализ — независимость P(Δe8=0) и P(Δe17=0|Δe8=0)
Тест 3: Поиск оптимального ΔW ≠ e₀ с минимальным HW-накоплением к раунду 8
Тест 4: Quartet-стоимость vs birthday-стоимость сравнение
Тест 5: Условный бумеранг — P(Δe17=0 | Δe8=0) vs P(Δe17=0)
"""

import random
import statistics
import math
from itertools import combinations

MASK = 0xFFFFFFFF

def rotr(x, n): return ((x >> n) | (x << (32-n))) & MASK
def sig0(x):  return rotr(x,7) ^ rotr(x,18) ^ (x>>3)
def sig1(x):  return rotr(x,17) ^ rotr(x,19) ^ (x>>10)
def Sig0(x):  return rotr(x,2) ^ rotr(x,13) ^ rotr(x,22)
def Sig1(x):  return rotr(x,6) ^ rotr(x,11) ^ rotr(x,25)
def Ch(e,f,g): return ((e&f) ^ (~e&g)) & MASK
def Maj(a,b,c): return (a&b) ^ (a&c) ^ (b&c)
def hw(x): return bin(x).count('1')

K = [0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,
     0x923f82a4,0xab1c5ed5,0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,
     0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,0xe49b69c1,0xefbe4786,
     0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
     0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,
     0x06ca6351,0x14292967,0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,
     0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,0xa2bfe8a1,0xa81a664b,
     0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
     0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,
     0x5b9cca4f,0x682e6ff3,0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,
     0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2]
IV = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
      0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]

def make_schedule(W16):
    W = list(W16) + [0]*48
    for i in range(16, 64):
        W[i] = (sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16]) & MASK
    return W

def make_schedule_xor(W16):
    """Расписание над XOR (линейное, точное для XOR-дифференциалов)."""
    W = list(W16) + [0]*48
    for i in range(16, 64):
        W[i] = sig1(W[i-2]) ^ W[i-7] ^ sig0(W[i-15]) ^ W[i-16]
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

def xor_diff_state(W_base, DW_xor, R):
    """XOR-дифференциал: вычислить Δstate_r = state_r(M XOR ΔM) XOR state_r(M)."""
    W2 = [W_base[i] ^ DW_xor[i] for i in range(16)]
    sn = sha_rounds(make_schedule(W_base), R)
    sf = sha_rounds(make_schedule(W2), R)
    # XOR-дифференциалы для всех регистров
    diff = [(sf[r][i] ^ sn[r][i]) for r in range(R+1) for i in range(8)]
    # Вернуть только e-регистр
    de = [(sf[r][4] ^ sn[r][4]) for r in range(R+1)]
    da = [(sf[r][0] ^ sn[r][0]) for r in range(R+1)]
    return de, da, sn, sf

print("=" * 72)
print("П-44A: БУМЕРАНГ — XOR-ГИБРИДНЫЙ АНАЛИЗ")
print("=" * 72)

# ─────────────────────────────────────────────────────────────
# Тест 1: P(Δe_r = 0) для каждого раунда при XOR ΔW[0]=1
# ─────────────────────────────────────────────────────────────
print("\n[1] XOR-дифференциал ΔW[0]=1 → P(Δe_r=0), P(Δa_r=0) для r=1..17")
print(f"{'r':>3} | {'P(Δe=0)':>10} | {'P(Δa=0)':>10} | {'E[HW(Δe)]':>10} | {'E[HW(Δa)]':>10}")
print("-" * 55)

N = 2000
DW_BASE = [0]*16
DW_BASE[0] = 1  # XOR-дифференциал: только бит 0 в W[0]

de_zero_cnt = {r: 0 for r in range(1, 18)}
da_zero_cnt = {r: 0 for r in range(1, 18)}
de_hw_sum   = {r: 0 for r in range(1, 18)}
da_hw_sum   = {r: 0 for r in range(1, 18)}

for _ in range(N):
    W0_val = random.randint(0, MASK)
    W1_val = random.randint(0, MASK)
    W_base = [W0_val, W1_val] + [0]*14
    de, da, sn, sf = xor_diff_state(W_base, DW_BASE, 17)
    for r in range(1, 18):
        if de[r] == 0: de_zero_cnt[r] += 1
        if da[r] == 0: da_zero_cnt[r] += 1
        de_hw_sum[r] += hw(de[r])
        da_hw_sum[r] += hw(da[r])

# Сохраняем для сплит-анализа
split_data_r8 = {}  # (de8, da8) → de17 статистика

for r in range(1, 18):
    p_de = de_zero_cnt[r] / N
    p_da = da_zero_cnt[r] / N
    e_hw_de = de_hw_sum[r] / N
    e_hw_da = da_hw_sum[r] / N
    marker = " ←SPLIT" if r == 8 else (" ←TARGET" if r == 17 else "")
    print(f"{r:>3} | {p_de:>10.6f} | {p_da:>10.6f} | {e_hw_de:>10.4f} | {e_hw_da:>10.4f}{marker}")

# ─────────────────────────────────────────────────────────────
# Тест 2: Сплит-анализ: независимость p₈ и q₁₇
# ─────────────────────────────────────────────────────────────
print("\n[2] СПЛИТ-АНАЛИЗ при r=8: независимость верхней/нижней части")
print("    Гипотеза: P(Δe17=0 | Δe8=0) != P(Δe17=0) → бумеранг-структура")

N2 = 5000
DW_SPLIT = [0]*16
DW_SPLIT[0] = 1

cnt_e8_zero = 0
cnt_e17_zero = 0
cnt_both_zero = 0
cnt_e17_given_e8 = 0
hw_de8_list = []
hw_de17_list = []
corr_data = []  # (de8, de17) для анализа корреляции

for _ in range(N2):
    W0_val = random.randint(0, MASK)
    W1_val = random.randint(0, MASK)
    W_base = [W0_val, W1_val] + [0]*14
    de, da, _, _ = xor_diff_state(W_base, DW_SPLIT, 17)
    de8  = de[8]
    de17 = de[17]
    hw_de8_list.append(hw(de8))
    hw_de17_list.append(hw(de17))
    corr_data.append((de8, de17))
    if de8  == 0: cnt_e8_zero += 1
    if de17 == 0: cnt_e17_zero += 1
    if de8 == 0 and de17 == 0: cnt_both_zero += 1
    if de8  == 0: cnt_e17_given_e8 += (1 if de17 == 0 else 0)

p_e8      = cnt_e8_zero / N2
p_e17     = cnt_e17_zero / N2
p_both    = cnt_both_zero / N2
p_e17_given_e8 = (cnt_e17_given_e8 / cnt_e8_zero) if cnt_e8_zero > 0 else 0

print(f"\n    P(Δe8=0)  = {p_e8:.6f}  (≈2^{{{math.log2(p_e8):.2f}}})" if p_e8 > 0 else "    P(Δe8=0) = 0")
print(f"    P(Δe17=0) = {p_e17:.6f}  (≈2^{{{math.log2(p_e17):.2f}}})" if p_e17 > 0 else "    P(Δe17=0) = 0")
if p_both > 0:
    print(f"    P(Δe8=0 AND Δe17=0) = {p_both:.6f}  (≈2^{{{math.log2(p_both):.2f}}})")
else:
    print(f"    P(Δe8=0 AND Δe17=0) = 0/{N2}")

# Если P(Δe8=0)·P(Δe17=0) = P(оба=0) → независимы → нет бумеранг-преимущества
if p_e8 > 0 and p_e17 > 0:
    independence_ratio = p_both / (p_e8 * p_e17) if (p_e8 * p_e17) > 0 else float('nan')
    print(f"    P(оба)/[P(e8)·P(e17)] = {independence_ratio:.4f}  "
          f"(1.0=независимы, >>1=бумеранг)")
    print(f"    P(Δe17=0 | Δe8=0) = {p_e17_given_e8:.6f}  "
          f"(vs безусловное P(Δe17=0) = {p_e17:.6f})")
    if p_e17 > 0:
        cond_ratio = p_e17_given_e8 / p_e17
        print(f"    Условное/безусловное = {cond_ratio:.4f}  "
              f"({'ЗАВИСИМОСТЬ → бумеранг!' if abs(cond_ratio-1)>0.1 else 'независимые'})")

# ─────────────────────────────────────────────────────────────
# Тест 3: Оптимальный однобитовый XOR ΔW — поиск минимального p₈
# ─────────────────────────────────────────────────────────────
print("\n[3] ОПТИМИЗАЦИЯ: поиск ΔW (1 бит в позиции k, раунд r) с мин. HW(Δe8)")
print("    (Минимальный HW(Δe8) → максимальное p₈ для верхней части бумеранга)")

N3 = 500
best_results = []
# Перебираем: 1 бит в каждом из W[0..7] (8 слов × 32 бита = 256 вариантов, берём часть)
test_positions = [(w_idx, bit) for w_idx in range(8) for bit in [0, 1, 7, 15, 16, 17, 31]]

for w_idx, bit in test_positions:
    DW_test = [0]*16
    DW_test[w_idx] = 1 << bit

    hw_e8_sum = 0
    e8_zero = 0
    for _ in range(N3):
        W_base = [random.randint(0, MASK), random.randint(0, MASK)] + [0]*14
        de, _, _, _ = xor_diff_state(W_base, DW_test, 8)
        hw_e8_sum += hw(de[8])
        if de[8] == 0: e8_zero += 1

    e_hw = hw_e8_sum / N3
    p_zero = e8_zero / N3
    best_results.append((e_hw, p_zero, w_idx, bit))

best_results.sort(key=lambda x: x[0])  # сортировка по E[HW] возрастанию
print(f"\n    {'W[i],bit':>12} | {'E[HW(Δe8)]':>12} | {'P(Δe8=0)':>12} | Оценка log₂(1/p₈)")
print("    " + "-" * 58)
for e_hw, p_zero, w_idx, bit in best_results[:10]:
    log2_inv = -math.log2(p_zero) if p_zero > 0 else float('inf')
    marker = " ← ЛУЧШИЙ" if (e_hw, p_zero, w_idx, bit) == best_results[0] else ""
    print(f"    W[{w_idx}], b{bit:02d}   | {e_hw:>12.4f} | {p_zero:>12.6f} | {log2_inv:>8.2f}{marker}")

best_ew, best_pz, best_wi, best_bit = best_results[0]
print(f"\n    → Лучший: W[{best_wi}], бит {best_bit}: E[HW(Δe8)]={best_ew:.4f}, P(Δe8=0)≈2^{{{-math.log2(best_pz):.2f} if best_pz > 0 else 'inf'}}")

# ─────────────────────────────────────────────────────────────
# Тест 4: Стоимость бумеранга vs birthday
# ─────────────────────────────────────────────────────────────
print("\n[4] СРАВНЕНИЕ СТОИМОСТИ: классика vs бумеранг")

# Используем результаты из теста 1
p8_xor  = de_zero_cnt[8]  / N   # P(Δe8=0 | ΔW[0]=1)
p17_xor = de_zero_cnt[17] / N   # P(Δe17=0 | ΔW[0]=1)

print(f"\n    XOR-дифференциал ΔW[0]=1:")
print(f"    p₈  = P(Δe8=0)  = {p8_xor:.6f}")
print(f"    p₁₇ = P(Δe17=0) = {p17_xor:.6f}")

if p8_xor > 0 and p17_xor > 0:
    # Классическая стоимость (прямая характеристика 1→17):
    C_classic = 1.0 / p17_xor
    log2_classic = math.log2(C_classic)

    # Бумеранг-стоимость (quartet): 4 / (p8 * q17)
    # где q17 = P(Δe17=0 | Δe8=0) = p_e17_given_e8 из теста 2
    q17 = p_e17_given_e8 if p_e17_given_e8 > 0 else p17_xor  # fallback

    C_boomerang = 4.0 / (p8_xor * q17)
    log2_boom = math.log2(C_boomerang)

    # Теоретический идеальный бумеранг (независимые p, q):
    C_boom_ideal = 4.0 / math.sqrt(p8_xor * p17_xor) if p8_xor * p17_xor > 0 else float('inf')
    log2_boom_ideal = math.log2(C_boom_ideal) if C_boom_ideal < float('inf') else float('inf')

    print(f"\n    Классическая стоимость:     O(2^{{{log2_classic:.2f}}})")
    print(f"    Бумеранг (реальный):        O(2^{{{log2_boom:.2f}}})")
    print(f"    Бумеранг (идеальный, √):    O(2^{{{log2_boom_ideal:.2f}}})")
    print(f"    Birthday (Wang ADD):         O(2^32.00)")

    if log2_boom < log2_classic:
        gain = log2_classic - log2_boom
        print(f"\n    → Бумеранг БЫСТРЕЕ классики на {gain:.2f} бит!")
    else:
        print(f"\n    → Бумеранг НЕ ЛУЧШЕ классики (XOR характеристика слабее ADD)")
else:
    print("    → Недостаточно данных (p8 или p17 = 0)")

# ─────────────────────────────────────────────────────────────
# Тест 5: Условный бумеранг — анализ зависимости Δe8 ↔ Δe17
# ─────────────────────────────────────────────────────────────
print("\n[5] КЛЮЧЕВОЙ ВОПРОС: КОРРЕЛИРУЮТ ли Δe8 и Δe17 (XOR)?")
print("    (Корреляция → есть смысл в бумеранге, нет → бесполезен)")

N5 = 3000
de8_vals  = []
de17_vals = []

for _ in range(N5):
    W_base = [random.randint(0, MASK), random.randint(0, MASK)] + [0]*14
    de, _, _, _ = xor_diff_state(W_base, [1]+[0]*15, 17)
    de8_vals.append(de[8])
    de17_vals.append(de[17])

# Битовая корреляция: для каждого бита b8, b17 измерить корреляцию
max_corr = 0.0
max_pair = (0, 0)
corr_matrix = []
for b8 in range(32):
    row = []
    for b17 in range(32):
        x8  = [(v >> b8)  & 1 for v in de8_vals]
        x17 = [(v >> b17) & 1 for v in de17_vals]
        m8  = sum(x8)  / N5
        m17 = sum(x17) / N5
        num = sum((x8[i]-m8)*(x17[i]-m17) for i in range(N5))
        std8  = math.sqrt(sum((x-m8)**2  for x in x8))
        std17 = math.sqrt(sum((x-m17)**2 for x in x17))
        corr  = num / (std8 * std17) if std8 * std17 > 0 else 0.0
        row.append(abs(corr))
        if abs(corr) > max_corr:
            max_corr = abs(corr)
            max_pair = (b8, b17)
    corr_matrix.append(row)

print(f"\n    Максимальная |corr(Δe8[b8], Δe17[b17])| = {max_corr:.6f}")
print(f"    Пара: бит {max_pair[0]} (Δe8) vs бит {max_pair[1]} (Δe17)")
print(f"    (Случайное ожидаемое: ≈0)")

# Топ-5 коррелированных пар
all_corr = [(corr_matrix[b8][b17], b8, b17)
            for b8 in range(32) for b17 in range(32)]
all_corr.sort(reverse=True)
print("\n    Топ-5 пар (|corr(Δe8[i], Δe17[j])|):")
for corr, b8, b17 in all_corr[:5]:
    print(f"      бит Δe8[{b8:2d}] ↔ Δe17[{b17:2d}]: |corr|={corr:.6f}")

# ─────────────────────────────────────────────────────────────
# ИТОГОВАЯ ТЕОРЕМА
# ─────────────────────────────────────────────────────────────
print("\n" + "=" * 72)
print("ИТОГ П-44A: ТЕОРЕМЫ")
print("=" * 72)

print("""
T_BOOMERANG_XOR_SPLIT:
  При XOR-дифференциале ΔW[0]=1:
    p₈  = P(Δe8=0)  ≈ 2^{log2(p8_xor):.2f}   (верхняя часть раунды 1-8)
    p₁₇ = P(Δe17=0) ≈ 2^{log2(p17_xor):.2f}  (полная характеристика)

  Если Δe8 и Δe17 НЕЗАВИСИМЫ (corr≈0):
    → Нет бумеранг-преимущества: P(оба)=p₈·p₁₇
    → Классика = бумеранг

  Если Δe8 и Δe17 ЗАВИСИМЫ (corr≠0):
    → Бумеранг экономит квадратный корень по стоимости
    → При max_corr>{0:.2f}: потенциальный выигрыш
""".format(max_corr))

if max_corr > 0.1:
    print("  СТАТУС: max_corr > 0.1 → ЕСТЬ ЗАВИСИМОСТЬ → исследовать дальше!")
elif max_corr > 0.05:
    print("  СТАТУС: слабая зависимость → возможен малый бумеранг-выигрыш")
else:
    print("  СТАТУС: Δe8 и Δe17 практически независимы при XOR-дифференциале")
    print("  T_BOOMERANG_XOR_INDEPENDENT: бумеранг не даёт преимущества в XOR-схеме")

print(f"""
T_XOR_CASCADE_DECAY:
  XOR-дифференциал ΔW[0]=1 через 8 раундов:
    E[HW(Δe8)]  = {de_hw_sum[8]/N:.4f} бит (лавина, ожидаемые 16)
    E[HW(Δe17)] = {de_hw_sum[17]/N:.4f} бит
  Нелинейность Ch/Maj полностью рассеивает XOR-дифференциал за r≈5 раундов.
  Это согласуется с T_DEGREE_BARRIER_8 и T_BARRIER_H5.
""")
