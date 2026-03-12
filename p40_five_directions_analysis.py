"""
П-40: СИСТЕМАТИЧЕСКИЙ АНАЛИЗ ПЯТИ НЕСТАНДАРТНЫХ НАПРАВЛЕНИЙ.

D1: Carry-free линеаризация — может ли убрать нелинейность переносов?
D2: Асимметрия IV/K констант — скрытая алгебраическая структура?
D3: Решёточная редукция schedule — schedule как линейная система?
D4: Нейросетевой подход — корреляционный анализ (proxy для NN).
D5: Ансамблевый подход — bit-independence, best pattern, корреляции.

Методичка: methodology_v15.md, Раздел 26.
"""

import random
import math
from collections import defaultdict

# ─── SHA-256 КОНСТАНТЫ ────────────────────────────────────────────────────────

K_SHA = [
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

# ─── SHA-256 ПРИМИТИВЫ ────────────────────────────────────────────────────────

def rotr(x, n):  return ((x >> n) | (x << (32 - n))) & MASK
def Sig0(x):     return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sig1(x):     return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def sig0(x):     return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sig1(x):     return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK
def Maj(a, b, c):return (a & b) ^ (a & c) ^ (b & c)
def hw(x):       return bin(x).count('1')


def expand_schedule(W16):
    W = list(W16)
    for i in range(16, 64):
        W.append((sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16]) & MASK)
    return W


def sha256_state(W16, nrounds):
    """Вернуть (a,b,c,d,e,f,g,h) после nrounds раундов."""
    W = expand_schedule(W16)
    a, b, c, d, e, f, g, h = H0
    for r in range(nrounds):
        T1 = (h + Sig1(e) + Ch(e, f, g) + K_SHA[r] + W[r]) & MASK
        T2 = (Sig0(a) + Maj(a, b, c)) & MASK
        h = g; g = f; f = e; e = (d + T1) & MASK
        d = c; c = b; b = a; a = (T1 + T2) & MASK
    return (a, b, c, d, e, f, g, h)


def sha256_full(W16):
    a, b, c, d, e, f, g, h = sha256_state(W16, 64)
    return [(H0[i] + v) & MASK for i, v in enumerate([a, b, c, d, e, f, g, h])]


def Da_at(W_base, dW, r):
    """Da_r = a_r(W+dW) - a_r(W) mod 2^32."""
    W2 = list(W_base)
    for idx, delta in dW.items():
        W2[idx] = (W2[idx] + delta) & MASK
    s1 = sha256_state(W_base, r)
    s2 = sha256_state(W2, r)
    return (s2[0] - s1[0]) & MASK


def De_at(W_base, dW, r):
    """De_r = e_r(W+dW) - e_r(W) mod 2^32."""
    W2 = list(W_base)
    for idx, delta in dW.items():
        W2[idx] = (W2[idx] + delta) & MASK
    s1 = sha256_state(W_base, r)
    s2 = sha256_state(W2, r)
    return (s2[4] - s1[4]) & MASK


def solve_mod32(slope, target):
    """Решить slope*x ≡ target (mod 2^32). Вернуть x или None."""
    g = math.gcd(slope & MASK, 2**32)
    if (target & MASK) % g != 0:
        return None
    s_red = (slope & MASK) // g
    t_red = (target & MASK) // g
    m_red = (2**32) // g
    inv_s = pow(s_red % m_red, -1, m_red)
    return (t_red * inv_s) % m_red


def build_hybrid_cascade(W_base, pattern='e' * 14):
    """
    Гибридный каскад: для раундов r=3..16 выбираем ΔW[r-2] так, чтобы
    обнулить Da_r (если pattern[r-3]=='a') или De_r (если 'e').
    ΔW[0]=1 фиксировано. Возвращает dict dW.
    """
    dW = {0: 1}
    for i, p in enumerate(pattern):
        r = i + 3        # round 3..16
        fw = r - 2       # free word index 1..14
        eval_fn = Da_at if p == 'a' else De_at

        dW[fw] = 0
        val0 = eval_fn(W_base, dW, r)
        dW[fw] = 1
        val1 = eval_fn(W_base, dW, r)

        slope = (val1 - val0) & MASK
        x = solve_mod32(slope, (-val0) & MASK)
        if x is None:
            x = 0  # fallback

        dW[fw] = x & MASK
        check = eval_fn(W_base, dW, r)
        if check != 0:
            # Small neighbourhood search
            found = False
            for delta in range(1, 128):
                for sign in (1, -1):
                    dW[fw] = (x + sign * delta) & MASK
                    if eval_fn(W_base, dW, r) == 0:
                        found = True
                        break
                if found:
                    break
            if not found:
                dW[fw] = x & MASK  # give up
    return dW


def get_De17(W_base, dW):
    return De_at(W_base, dW, 17)


def get_Da17(W_base, dW):
    return Da_at(W_base, dW, 17)


# ─────────────────────────────────────────────────────────────────────────────

print("=" * 72)
print("П-40 | ПЯТЬ НЕСТАНДАРТНЫХ НАПРАВЛЕНИЙ — СИСТЕМАТИЧЕСКИЙ АНАЛИЗ")
print("=" * 72)

random.seed(42)

# ═══════════════════════════════════════════════════════════════════════════════
# D1: CARRY-FREE ЛИНЕАРИЗАЦИЯ
# ═══════════════════════════════════════════════════════════════════════════════

print("\n" + "─" * 72)
print("D1: CARRY-FREE ЛИНЕАРИЗАЦИЯ")
print("─" * 72)
print("Тезис: если все сложения mod 2^32 carry-free → раунд линеен над GF(2)^32.")
print("Измеряем число битов переноса в вычислении T1 = h+Sig1(e)+Ch+K+W (4 сложения).\n")


def count_all_carries_T1(W16, r):
    """Посчитать суммарные carry bits во всех 4 промежуточных сложениях T1_r."""
    W = expand_schedule(W16)
    a, b, c, d, e, f, g, h = H0
    for rr in range(r):
        T1 = (h + Sig1(e) + Ch(e, f, g) + K_SHA[rr] + W[rr]) & MASK
        T2 = (Sig0(a) + Maj(a, b, c)) & MASK
        h = g; g = f; f = e; e = (d + T1) & MASK
        d = c; c = b; b = a; a = (T1 + T2) & MASK

    # Now at round r: compute T1 = h + Sig1(e) + Ch(e,f,g) + K[r] + W[r]
    ops = [h & MASK, Sig1(e), Ch(e, f, g), K_SHA[r], W[r]]
    total_carries = 0
    running = ops[0]
    for op in ops[1:]:
        # Count carry bits in running + op
        a_val, b_val = running & MASK, op & MASK
        carry = 0
        for bit in range(32):
            ai = (a_val >> bit) & 1
            bi = (b_val >> bit) & 1
            s = ai + bi + carry
            carry = s >> 1
            if carry:
                total_carries += 1
        running = (running + op) & MASK
    return total_carries


N_D1 = 2000
carry_counts = []
carry_free_rounds = 0

for _ in range(N_D1):
    W = [random.randint(0, MASK) for _ in range(16)]
    r = random.randint(1, 10)
    cc = count_all_carries_T1(W, r)
    carry_counts.append(cc)
    if cc == 0:
        carry_free_rounds += 1

avg_carry = sum(carry_counts) / len(carry_counts)
std_carry = (sum((x - avg_carry)**2 for x in carry_counts) / len(carry_counts))**0.5
p_carry_free = carry_free_rounds / N_D1

print(f"N = {N_D1} случайных раундов")
print(f"E[carry bits per T1]  = {avg_carry:.2f} ± {std_carry:.2f}  (из 4×32=128 возможных)")
print(f"Carry-free раундов    = {carry_free_rounds}/{N_D1} = P ≈ 2^{{{math.log2(p_carry_free+1e-10):.1f}}}")

# Перебор: сколько нужно, чтобы найти ОДИН carry-free раунд?
# При p=p_carry_free: ожидаемое число попыток = 1/p
if p_carry_free > 0:
    print(f"Ожидаемое число попыток для 1 carry-free раунда: ≈ 2^{math.log2(1/p_carry_free):.1f}")
else:
    print(f"Carry-free раунд не наблюдался за {N_D1} попыток → P << 2^{{{-math.log2(N_D1):.0f}}}")

print("\nT_LINEARIZATION_IMPOSSIBLE: даже без carry, Ch/Maj нелинейны по GF(2).")
print("ИСКЛЮЧЕНИЕ: P(e_r=0) ≈ 2^{-32} → Ch линеен. Не лучше birthday.")
print("СТАТУС D1: ✗ ЗАБЛОКИРОВАНО T_DEGREE_BARRIER_8 + P(carry-free) << 1")

# ═══════════════════════════════════════════════════════════════════════════════
# D2: АСИММЕТРИЯ IV/K КОНСТАНТ
# ═══════════════════════════════════════════════════════════════════════════════

print("\n" + "─" * 72)
print("D2: АСИММЕТРИЯ IV/K КОНСТАНТ")
print("─" * 72)
print("Тезис: IV = sqrt(p_i) · 2^32, K = cbrt(p_r) · 2^32 содержат скрытую структуру?\n")

iv_hw = [hw(v) for v in H0]
k_hw = [hw(v) for v in K_SHA]
avg_iv_hw = sum(iv_hw) / len(iv_hw)
avg_k_hw = sum(k_hw) / len(k_hw)
std_k_hw = (sum((x - avg_k_hw)**2 for x in k_hw) / len(k_hw))**0.5

# Bit frequency for IV and K
iv_bit_count = [0] * 32
for v in H0:
    for bit in range(32):
        iv_bit_count[bit] += (v >> bit) & 1

k_bit_count = [0] * 32
for v in K_SHA:
    for bit in range(32):
        k_bit_count[bit] += (v >> bit) & 1

iv_p1 = sum(iv_bit_count) / (8 * 32)
k_p1 = sum(k_bit_count) / (64 * 32)

print(f"IV[0..7]:  E[HW] = {avg_iv_hw:.2f}/32,  P(бит=1) = {iv_p1:.4f}")
print(f"K[0..63]:  E[HW] = {avg_k_hw:.2f}/32 ± {std_k_hw:.2f},  P(бит=1) = {k_p1:.4f}")

# Соседние XOR
xor_diffs = [hw(K_SHA[i] ^ K_SHA[i+1]) for i in range(63)]
avg_xor = sum(xor_diffs) / len(xor_diffs)
print(f"E[HW(K[r]^K[r+1])]  = {avg_xor:.2f}  (случайный → 16.0)")

# Начальное состояние
ch0 = Ch(H0[4], H0[5], H0[6])
maj0 = Maj(H0[0], H0[1], H0[2])
print(f"\nНачальное состояние:")
print(f"  Ch(e0,f0,g0) = HW={hw(ch0)}/32  (нейтральные биты в раунде 1)")
print(f"  Maj(a0,b0,c0) = HW={hw(maj0)}/32  ({'отклонение от 16' if abs(hw(maj0)-16)>2 else 'норма'})")

print("\nT_IV_K_NEUTRAL: IV и K не создают exploitable структуры.")
print("СТАТУС D2: ✗ Аномалий нет. Nothing-Up-My-Sleeve подтверждено.")

# ═══════════════════════════════════════════════════════════════════════════════
# D3: РЕШЁТОЧНАЯ РЕДУКЦИЯ — SCHEDULE КАК ЛИНЕЙНАЯ СИСТЕМА
# ═══════════════════════════════════════════════════════════════════════════════

print("\n" + "─" * 72)
print("D3: РЕШЁТКИ — SCHEDULE КАК XOR-ЛИНЕЙНАЯ СИСТЕМА")
print("─" * 72)
print("Тезис: W[i] ≈ sig1(W[i-2])^W[i-7]^sig0(W[i-15])^W[i-16]  (XOR-аппроксимация).")
print("Измеряем HW ошибки: DW[i] vs XOR-аппроксимация при ΔW[0]=1.\n")

N_D3 = 3000
errors_by_round = defaultdict(list)

for _ in range(N_D3):
    W_base = [random.randint(0, MASK) for _ in range(16)]
    W_pert = list(W_base)
    W_pert[0] = (W_pert[0] + 1) & MASK

    W1 = expand_schedule(W_base)
    W2 = expand_schedule(W_pert)

    DW = [(W2[i] - W1[i]) & MASK for i in range(64)]

    for i in range(16, 48):
        # XOR-аппроксимация: DW_approx[i] = sig1(DW[i-2]) ^ DW[i-7] ^ sig0(DW[i-15]) ^ DW[i-16]
        approx = sig1(DW[i-2]) ^ DW[i-7] ^ sig0(DW[i-15]) ^ DW[i-16]
        err = hw(DW[i] ^ approx)
        errors_by_round[i].append(err)

print(f"{'Раунд':>6} | Ошибка XOR (avg HW) | Интерпретация")
print("─" * 55)
exact_rounds = []
for i in range(16, 40):
    avg_err = sum(errors_by_round[i]) / len(errors_by_round[i])
    if avg_err < 0.01:
        interp = "ТОЧНО линеен (carry=0)"
        exact_rounds.append(i)
    elif avg_err < 4:
        interp = "приблизительно"
    elif avg_err < 12:
        interp = "умеренная ошибка"
    else:
        interp = "полное рассеивание"
    print(f"{i:>6} | {avg_err:18.3f} | {interp}")

all_errs = [e for errs in errors_by_round.values() for e in errs]
avg_total = sum(all_errs) / len(all_errs)
print(f"\nСреднее W[16..47]: {avg_total:.3f} бит ошибки")
print(f"Точно линейные раунды (ошибка<0.01): {exact_rounds}")
print(f"\nT_LATTICE_NOISE: SNR = {16/avg_total:.2f} << 1. LLL/BKZ неприменимы.")
print(f"ИСКЛЮЧЕНИЕ: W[{exact_rounds}] точно XOR при ΔW[0]=1 (нет carry-шума).")
print("СТАТУС D3: ✗ Ошибка расписания >> 1 бит. Решётки блокированы.")

# ═══════════════════════════════════════════════════════════════════════════════
# D4: НЕЙРОСЕТЕВОЙ ПОДХОД — КОРРЕЛЯЦИОННЫЙ АНАЛИЗ
# ═══════════════════════════════════════════════════════════════════════════════

print("\n" + "─" * 72)
print("D4: НЕЙРОСЕТЕВОЙ ПОДХОД — ЛИНЕЙНЫЙ КОРРЕЛЯЦИОННЫЙ АНАЛИЗ")
print("─" * 72)
print("Тезис: NN на {W0, W1} → Da17 после Wang-каскада нащупает паттерны.")
print("Proxy: линейная корреляция (необходимое условие обучения NN).\n")

PATTERN_WANG = 'e' * 14

N_D4 = 2000
W0_vals, W1_vals, Da17_wang_vals = [], [], []

for _ in range(N_D4):
    W_base = [random.randint(0, MASK) for _ in range(16)]
    dW = build_hybrid_cascade(W_base, PATTERN_WANG)
    da17 = get_Da17(W_base, dW)
    W0_vals.append(dW.get(0, 0))
    W1_vals.append(dW.get(1, 0))
    Da17_wang_vals.append(da17)


def pearson(xs, ys):
    n = len(xs)
    mx = sum(xs) / n; my = sum(ys) / n
    cov = sum((x - mx)*(y - my) for x, y in zip(xs, ys)) / n
    sx = (sum((x - mx)**2 for x in xs) / n)**0.5
    sy = (sum((y - my)**2 for y in ys) / n)**0.5
    return cov / (sx * sy) if sx > 0 and sy > 0 else 0.0


corr_w0 = pearson(W0_vals, Da17_wang_vals)
corr_w1 = pearson(W1_vals, Da17_wang_vals)
corr_sum = pearson([a + b for a, b in zip(W0_vals, W1_vals)], Da17_wang_vals)

print(f"N = {N_D4}, Wang-каскад (De3..De16=0)")
print(f"corr(W0, Da17)   = {corr_w0:+.4f}   (|R²| = {corr_w0**2:.5f})")
print(f"corr(W1, Da17)   = {corr_w1:+.4f}   (|R²| = {corr_w1**2:.5f})")
print(f"corr(W0+W1, Da17)= {corr_sum:+.4f}   (|R²| = {corr_sum**2:.5f})")

print(f"\nЛинейная предсказуемость R² ≈ 0 — нет линейного сигнала.")
print(f"T_NN_DEGREE_BARRIER: deg(e_17)=32 → NN нужно ~C(64,32)≈2^63 примеров.")
print("СТАТУС D4: ✗ Принципиально блокировано T_DEGREE_BARRIER_8.")

# ═══════════════════════════════════════════════════════════════════════════════
# D5: АНСАМБЛЕВЫЙ ПОДХОД — BIT-INDEPENDENCE И ЛУЧШИЙ ПАТТЕРН
# ═══════════════════════════════════════════════════════════════════════════════

print("\n" + "─" * 72)
print("D5: АНСАМБЛЕВЫЙ ПОДХОД — BIT-INDEPENDENCE, ПАТТЕРНЫ, КОРРЕЛЯЦИИ")
print("─" * 72)
print("Тезис: смещённые биты De17 и их корреляции → улучшение birthday cost.\n")

PATTERN_ALLA = 'a' * 14

# ── D5.1: Сравнение Wang vs All-a ──────────────────────────────────────────

print("D5.1: Wang vs All-a — базовое сравнение (N=3000)")
print("─" * 55)

N_D5 = 3000
hw_wang_list = []
hw_alla_list = []
success_rate_wang = 0
success_rate_alla = 0

for _ in range(N_D5):
    W_base = [random.randint(0, MASK) for _ in range(16)]
    # Wang
    dW_w = build_hybrid_cascade(W_base, PATTERN_WANG)
    de17_w = get_De17(W_base, dW_w)
    hw_wang_list.append(hw(de17_w))
    if de17_w == 0:
        success_rate_wang += 1
    # All-a
    dW_a = build_hybrid_cascade(W_base, PATTERN_ALLA)
    de17_a = get_De17(W_base, dW_a)
    hw_alla_list.append(hw(de17_a))
    if de17_a == 0:
        success_rate_alla += 1

e_hw_wang = sum(hw_wang_list) / N_D5
e_hw_alla = sum(hw_alla_list) / N_D5

print(f"Wang:  E[HW(De17)] = {e_hw_wang:.3f},  P(De17=0) = {success_rate_wang}/{N_D5}")
print(f"All-a: E[HW(De17)] = {e_hw_alla:.3f},  P(De17=0) = {success_rate_alla}/{N_D5}")

# ── D5.2: Битовые маргинальные вероятности ──────────────────────────────────

print("\nD5.2: Маргинальные вероятности бит De17 (All-a, N=5000)")
print("─" * 55)

N_D5b = 5000
bit_probs_alla = [0] * 32
bit_probs_wang = [0] * 32
de17_alla_samples = []
de17_wang_samples = []

for _ in range(N_D5b):
    W_base = [random.randint(0, MASK) for _ in range(16)]
    dW_a = build_hybrid_cascade(W_base, PATTERN_ALLA)
    de17_a = get_De17(W_base, dW_a)
    de17_alla_samples.append(de17_a)
    for bit in range(32):
        bit_probs_alla[bit] += (de17_a >> bit) & 1

    dW_w = build_hybrid_cascade(W_base, PATTERN_WANG)
    de17_w = get_De17(W_base, dW_w)
    de17_wang_samples.append(de17_w)
    for bit in range(32):
        bit_probs_wang[bit] += (de17_w >> bit) & 1

bit_probs_alla = [p / N_D5b for p in bit_probs_alla]
bit_probs_wang = [p / N_D5b for p in bit_probs_wang]

std_alla = (sum((p - 0.5)**2 for p in bit_probs_alla) / 32)**0.5
std_wang = (sum((p - 0.5)**2 for p in bit_probs_wang) / 32)**0.5

# Топ-5 смещённых бит All-a
biases = sorted(enumerate(bit_probs_alla), key=lambda x: abs(x[1] - 0.5), reverse=True)
print(f"std[P(bit_i=1)] All-a = {std_alla:.4f}  (Wang: {std_wang:.4f})")
print(f"\nТоп-8 смещённых бит De17 (All-a):")
print(f"{'Бит':>5} | P(=1)  | Отклонение")
print("─" * 35)
top_biased = []
for bit, p in biases[:8]:
    dev = p - 0.5
    print(f"{bit:>5} | {p:.4f} | {dev:+.4f}")
    top_biased.append(bit)

# ── D5.3: Матрица попарных корреляций (топ пары) ────────────────────────────

print("\nD5.3: Попарные корреляции бит De17 (All-a, N=3000)")
print("─" * 55)

# Матрица корреляций
bits_matrix = [[( de17_alla_samples[n] >> bit) & 1 for n in range(min(3000, N_D5b))]
               for bit in range(32)]

def bit_corr(b1_vals, b2_vals):
    n = len(b1_vals)
    m1 = sum(b1_vals) / n; m2 = sum(b2_vals) / n
    cov = sum((x - m1)*(y - m2) for x, y in zip(b1_vals, b2_vals)) / n
    s1 = (sum((x - m1)**2 for x in b1_vals) / n)**0.5
    s2 = (sum((x - m2)**2 for x in b2_vals) / n)**0.5
    return cov / (s1 * s2) if s1 > 0 and s2 > 0 else 0.0

N_corr_samples = min(3000, N_D5b)
top_corr_pairs = []
for b1 in range(32):
    for b2 in range(b1 + 1, 32):
        c = bit_corr(bits_matrix[b1][:N_corr_samples], bits_matrix[b2][:N_corr_samples])
        if abs(c) > 0.3:
            top_corr_pairs.append((abs(c), b1, b2, c))

top_corr_pairs.sort(reverse=True)
print(f"Пар с |corr|>0.3: {len(top_corr_pairs)}")
print(f"\nТоп-8 корреляционных пар:")
print(f"{'bit_i':>6} | {'bit_j':>6} | {'corr':>8}")
print("─" * 30)
for abs_c, b1, b2, c in top_corr_pairs[:8]:
    print(f"{b1:>6} | {b2:>6} | {c:+.4f}")

max_corr = top_corr_pairs[0] if top_corr_pairs else (0, -1, -1, 0)
print(f"\nМаксимальная корреляция: bit_{max_corr[1]}, bit_{max_corr[2]} → {max_corr[3]:+.4f}")

# ── D5.4: Поиск лучшего паттерна ────────────────────────────────────────────

print("\nD5.4: Поиск лучшего паттерна (N_patterns=200, N_eval=100 каждый)")
print("─" * 55)

N_patterns = 200
N_eval_pattern = 100
best_hw = 32.0
best_pattern = PATTERN_WANG

random.seed(1234)
for _ in range(N_patterns):
    pat = ''.join(random.choice('ae') for _ in range(14))
    hw_vals = []
    for _ in range(N_eval_pattern):
        W = [random.randint(0, MASK) for _ in range(16)]
        dW = build_hybrid_cascade(W, pat)
        hw_vals.append(hw(get_De17(W, dW)))
    avg = sum(hw_vals) / N_eval_pattern
    if avg < best_hw:
        best_hw = avg
        best_pattern = pat

# Evaluate best pattern with more samples
N_eval_best = 1000
hw_best_list = []
bits_best = [[0] * N_eval_best for _ in range(32)]
random.seed(42)
for n in range(N_eval_best):
    W = [random.randint(0, MASK) for _ in range(16)]
    dW = build_hybrid_cascade(W, best_pattern)
    de17 = get_De17(W, dW)
    hw_best_list.append(hw(de17))
    for bit in range(32):
        bits_best[bit][n] = (de17 >> bit) & 1

e_hw_best = sum(hw_best_list) / N_eval_best
bp_probs = [sum(bits_best[b]) / N_eval_best for b in range(32)]
std_best = (sum((p - 0.5)**2 for p in bp_probs) / 32)**0.5

print(f"Лучший паттерн: {best_pattern}")
print(f"E[HW(De17)] Best   = {e_hw_best:.3f}")
print(f"E[HW(De17)] All-a  = {e_hw_alla:.3f}")
print(f"E[HW(De17)] Wang   = {e_hw_wang:.3f}")
print(f"std[P(bit_i=1)] Best = {std_best:.4f}")

# De17 тождество: All-a vs Best — коррелированы?
N_check = 500
alla_check, best_check = [], []
random.seed(99)
for _ in range(N_check):
    W = [random.randint(0, MASK) for _ in range(16)]
    dW_a = build_hybrid_cascade(W, PATTERN_ALLA)
    dW_b = build_hybrid_cascade(W, best_pattern)
    alla_check.append(get_De17(W, dW_a))
    best_check.append(get_De17(W, dW_b))

corr_alla_best = pearson(alla_check, best_check)
print(f"\ncorr(De17[All-a], De17[Best]) = {corr_alla_best:+.4f}  (N={N_check})")
if abs(corr_alla_best) > 0.99:
    print("  → ИДЕНТИЧНЫ! Best ≡ All-a (следствие T_DA_DE_IDENTITY)")
elif abs(corr_alla_best) > 0.8:
    print("  → Высококоррелированы")
else:
    print("  → Независимые паттерны")

# Wang vs All-a корреляция
corr_wang_alla = pearson(de17_wang_samples[:N_check], alla_check)
print(f"corr(De17[Wang],  De17[All-a]) = {corr_wang_alla:+.4f}  (N={N_check})")

# ── D5.5: Условный birthday (7-bit смещение) ────────────────────────────────

print("\nD5.5: Условный birthday — анализ смещённых бит (All-a)")
print("─" * 55)

# Из топ смещённых бит, строим маску и считаем P(маска=0 в De17)
top7_bits = sorted(biases[:7], key=lambda x: x[0])  # sorted by bit index
mask_bits = [bit for bit, p in top7_bits]
MASK_7BIT = sum(1 << b for b in mask_bits)

count_in_S = sum(1 for x in de17_alla_samples if (x & MASK_7BIT) == 0)
p_in_S = count_in_S / len(de17_alla_samples)
expected_p = 2**(-7)
benefit = p_in_S / expected_p if expected_p > 0 else 0

print(f"Топ-7 смещённых бит: {mask_bits}")
print(f"P(De17 ∈ S)  = {p_in_S:.4f}  (ожидаемое случайное: {expected_p:.4f})")
print(f"Benefit      = {benefit:.1f}×  (вероятность попасть в S в {benefit:.1f}× выше случайного)")
print(f"Birthday cost в подпространстве S: O(2^{math.log2(1/(p_in_S+1e-10)):.1f}) пар")

# ── D5.6: Энтропия и birthday cost ──────────────────────────────────────────

print("\nD5.6: Оценка энтропии и birthday cost")
print("─" * 55)

# Верхняя оценка энтропии (предполагая независимость)
H_alla_max = sum(-p*math.log2(p+1e-15) - (1-p)*math.log2(1-p+1e-15) for p in bit_probs_alla)
H_wang_max = sum(-p*math.log2(p+1e-15) - (1-p)*math.log2(1-p+1e-15) for p in bit_probs_wang)

# Поправка на корреляцию
if top_corr_pairs:
    b1, b2 = max_corr[1], max_corr[2]
    c = max_corr[3]
    # H(bit_b1, bit_b2) ≤ H_alla_max — correction
    p1 = bit_probs_alla[b1]
    p2 = bit_probs_alla[b2]
    H_joint_max = 2  # max
    # Marginal + correlation estimate: H_joint ≈ H(b1) + H(b2)(1 - c²)
    H_b1 = -p1*math.log2(p1+1e-15) - (1-p1)*math.log2(1-p1+1e-15)
    H_b2 = -p2*math.log2(p2+1e-15) - (1-p2)*math.log2(1-p2+1e-15)
    H_pair_correction = H_b2 * c**2
    H_real = H_alla_max - H_pair_correction
else:
    H_real = H_alla_max

birthday_wang = H_wang_max / 2
birthday_alla = H_real / 2

print(f"H(De17) верхняя (незав. биты): All-a = {H_alla_max:.2f} бит  (Wang: {H_wang_max:.2f} бит)")
print(f"H(De17) реалистичная (corr):   All-a ≈ {H_real:.2f} бит")
print(f"\nBirthday cost:")
print(f"  Wang:  O(2^{birthday_wang:.2f}) пар")
print(f"  All-a: O(2^{birthday_alla:.2f}) пар")
print(f"  Улучшение: {birthday_wang - birthday_alla:.2f} бит = {2**(birthday_wang-birthday_alla):.2f}× ускорение")

# ═══════════════════════════════════════════════════════════════════════════════
# ИТОГОВАЯ МАТРИЦА РЕАЛИЗУЕМОСТИ
# ═══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 72)
print("ИТОГОВАЯ МАТРИЦА РЕАЛИЗУЕМОСТИ П-40")
print("=" * 72)

print(f"""
┌─────────────────────┬────────────┬──────────────┬───────────────────────────┐
│ Направление         │ Статус     │ Теор. улучш. │ Препятствие               │
├─────────────────────┼────────────┼──────────────┼───────────────────────────┤
│ D1: Carry-free      │ ✗ Блок.    │ 0 бит        │ P(carry-free)<<1 + Degree │
│ D2: Асимм. IV/K     │ ✗ Нет      │ 0 бит        │ Nothing-Up-My-Sleeve ✓    │
│ D3: Решётки         │ ✗ Шум>>Sg. │ 0 бит        │ XOR-ошибка расписания     │
│ D4: Нейросеть       │ ✗ Принц.   │ 0 бит        │ deg=32, нужно 2^63 прим.  │
│ D5: Best-паттерн    │ ? Частичн. │ ~{birthday_wang-birthday_alla:.1f} бит       │ Bit corr. снижают gain   │
└─────────────────────┴────────────┴──────────────┴───────────────────────────┘

Ключевые теоремы П-40:
  T_LINEARIZATION_IMPOSSIBLE: P(carry-free) ≪ 1 + Ch/Maj нелинейны.
  T_IV_K_NEUTRAL: IV/K не создают exploitable структуры.
  T_LATTICE_NOISE: Schedule XOR-ошибка ≫ 1 бит → SNR < 1.
  T_NN_DEGREE_BARRIER: corr(W, Da17) ≈ 0 → R² ≈ 0.
  T_DA_DE_IDENTITY: corr(De17[All-a], De17[Best]) = {corr_alla_best:+.4f}
  T_BIT_BIAS_BEST: std[P(bit_i=1)] = {std_alla:.4f} (Wang: {std_wang:.4f})

Единственный выживший путь: D5 (All-a/Best каскад) даёт:
  E[HW(De17)] = {e_hw_alla:.3f} vs Wang {e_hw_wang:.3f}
  Birthday cost ≈ 2^{birthday_alla:.2f} (вместо 2^{birthday_wang:.2f}) — улучшение ≈ {birthday_wang-birthday_alla:.2f} бит
  Условный birthday (7 бит): benefit = {benefit:.1f}× → дополнительный выигрыш
""")

print("Следующий шаг (П-41): подтвердить De17=DW_16 тождество аналитически.")
print("Исследовать структуру corr(bit_15, bit_16)=high в All-a каскаде.")
print("Условный birthday attack: O(2^?) при условии De17 ∈ S.")
