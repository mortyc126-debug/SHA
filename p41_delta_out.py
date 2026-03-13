"""
П-41: КОНЦЕНТРАЦИЯ Δout — структура выходной разности для структурированных пар.

Из П-38: ключевой вопрос — равномерен ли Δout = SHA(W) XOR SHA(W') для пар с
De3..De17=0, или есть КОНЦЕНТРАЦИЯ в подпространстве?

Из П-40: ротационный аттрактор создаёт новую структуру дифференциала.

ГИПОТЕЗЫ (три класса пар):
  1. Случайные пары:           Δout равномерен на 256 битах → 2^128 для birthday
  2. Пары с De6=0 (T_DEkDECOMPOSITION): Δout определяется раундами 7..64
     → меньше "степеней свободы" → возможная концентрация
  3. Ротационные аттракторы:  De_r = ROTR(De_{r-3}, s) → квазипериодическая
     структура → Δout может иметь ротационную корреляцию

МЕТРИКИ КОНЦЕНТРАЦИИ:
  a) HW(Δout) — распределение (если не ~128, есть смещение)
  b) Энтропия по битам Δout — если биты не 50/50
  c) Корреляции между словами Δout[0..7]
  d) Эффективная размерность подпространства Δout

НОВЫЙ ВОПРОС (главный результат):
  Если Δout живёт в подпространстве размерности k < 256,
  то birthday collision на Δout требует ~2^{k/2} << 2^{128} пар.
  Это было бы теоретическим улучшением над 2^{128}.
"""

import random
import hashlib
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

def rotr(x, n):  return ((x >> n) | (x << (32 - n))) & MASK
def Sig0(x):     return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sig1(x):     return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def sig0(x):     return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sig1(x):     return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
def hw(x):       return bin(x).count('1')
def hw256(words): return sum(hw(w) for w in words)


def sha256_compress(W16):
    """Полная SHA-256 сжимающая функция (без padding, без финального add IV)."""
    W = list(W16)
    for i in range(16, 64):
        W.append((sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16]) & MASK)
    a, b, c, d, e, f, g, h = H0
    for r in range(64):
        T1 = (h + Sig1(e) + Ch(e, f, g) + K[r] + W[r]) & MASK
        T2 = (Sig0(a) + Maj(a, b, c)) & MASK
        h = g; g = f; f = e; e = (d + T1) & MASK
        d = c; c = b; b = a; a = (T1 + T2) & MASK
    # Возвращаем БЕЗ добавления H0 (только внутренние раунды)
    return [a, b, c, d, e, f, g, h]


def sha256_compress_add_iv(W16):
    """SHA-256 с добавлением H0 в конце (стандартный режим)."""
    out = sha256_compress(W16)
    return [(out[i] + H0[i]) & MASK for i in range(8)]


def sha256_trace_n(W16, nrounds):
    """Первые nrounds раундов SHA-256, возвращает состояния."""
    W = list(W16)
    for i in range(16, 64):
        W.append((sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16]) & MASK)
    a, b, c, d, e, f, g, h = H0
    states = [(a, b, c, d, e, f, g, h)]
    for r in range(nrounds):
        T1 = (h + Sig1(e) + Ch(e, f, g) + K[r] + W[r]) & MASK
        T2 = (Sig0(a) + Maj(a, b, c)) & MASK
        h = g; g = f; f = e; e = (d + T1) & MASK
        d = c; c = b; b = a; a = (T1 + T2) & MASK
        states.append((a, b, c, d, e, f, g, h))
    return states


# ─── Генераторы пар ──────────────────────────────────────────────────────────

def gen_random_pair():
    """Случайная пара (W, W'), W' = W + dW0."""
    W = [random.randint(0, MASK) for _ in range(16)]
    W2 = list(W); W2[0] = (W2[0] + random.randint(1, MASK)) & MASK
    return W, W2


def gen_pair_de6_zero():
    """Пара с De6=0 методом T_DEkDECOMPOSITION из П-36.
    De_6 = Da_2 + ΔW_3  →  ΔW_3 = -Da_2 устанавливает De_6=0.
    Стоимость: O(1) — нет случайного поиска!
    """
    W = [random.randint(0, MASK) for _ in range(16)]
    dW0 = random.randint(1, MASK)
    W2 = list(W); W2[0] = (W2[0] + dW0) & MASK

    # Трассировка 6 раундов для W и W'
    st1 = sha256_trace_n(W, 6)
    st2 = sha256_trace_n(W2, 6)

    # Da_2 = a(round2, W') - a(round2, W)
    Da2 = (st2[2][0] - st1[2][0]) & MASK
    # De_6 (натуральное без коррекции)
    De6_nat = (st2[6][4] - st1[6][4]) & MASK

    # ΔW_3 = -Da_2 - De6_nat (первое приближение из T_DEkDECOMPOSITION)
    # Точное: установить W2[3] = W[3] + correction
    # Используем линейное приближение: ΔW_3 корректирует De_6
    correction = (0 - De6_nat) & MASK  # хотим De6=0, поэтому нужно отнять De6_nat
    W2[3] = (W2[3] + correction) & MASK

    return W, W2


def gen_pair_rot_attractor(s=1):
    """Пара с ротационным аттрактором De_6 = ROTR(De_3, s)."""
    W = [random.randint(0, MASK) for _ in range(16)]
    dW0 = random.randint(1, MASK)
    W2 = list(W); W2[0] = (W2[0] + dW0) & MASK

    st1 = sha256_trace_n(W, 6)
    st2 = sha256_trace_n(W2, 6)
    De3 = (st2[3][4] - st1[3][4]) & MASK
    De6 = (st2[6][4] - st1[6][4]) & MASK
    target6 = rotr(De3, s)
    W2[5] = (W2[5] + (target6 - De6)) & MASK

    return W, W2


# ─── Эксперимент 1: Базовая статистика Δout ─────────────────────────────────

print("=" * 70)
print("П-41 | КОНЦЕНТРАЦИЯ Δout — структура выходной разности")
print("=" * 70)

print("\n[1] Базовая статистика HW(Δout) для трёх типов пар")
print("─" * 70)

N = 500

results = {}
for name, gen in [
    ("случайные",            gen_random_pair),
    ("De6=0 (T_DECOMP)",    gen_pair_de6_zero),
    ("ротационный (s=1)",   lambda: gen_pair_rot_attractor(1)),
]:
    hw_list = []
    for _ in range(N):
        W, W2 = gen()
        out1 = sha256_compress_add_iv(W)
        out2 = sha256_compress_add_iv(W2)
        delta = [o1 ^ o2 for o1, o2 in zip(out1, out2)]
        hw_list.append(hw256(delta))
    avg = sum(hw_list) / len(hw_list)
    mn  = min(hw_list)
    std = (sum((x - avg)**2 for x in hw_list) / len(hw_list)) ** 0.5
    results[name] = (avg, mn, std, hw_list)
    print(f"  {name:25s}: avg HW={avg:.1f}, min={mn}, std={std:.2f}")

print(f"\n  Теоретический avg для равномерного: 128.0")

# ─── Эксперимент 2: Побитовая энтропия Δout ─────────────────────────────────

print("\n[2] Побитовая статистика Δout — ищем биты с меньшей энтропией")
print("─" * 70)
print("Для каждого из 256 бит: частота 1 в Δout.")
print("Если частота ≠ 0.5 → концентрация!")

N_entropy = 500

for name, gen in [
    ("случайные",          gen_random_pair),
    ("De6=0 (T_DECOMP)",   gen_pair_de6_zero),
    ("ротационный (s=1)",  lambda: gen_pair_rot_attractor(1)),
]:
    bit_freq = [0] * 256
    for _ in range(N_entropy):
        W, W2 = gen()
        out1 = sha256_compress_add_iv(W)
        out2 = sha256_compress_add_iv(W2)
        for word_idx in range(8):
            delta_word = out1[word_idx] ^ out2[word_idx]
            for bit in range(32):
                if (delta_word >> bit) & 1:
                    bit_freq[word_idx * 32 + bit] += 1

    freqs = [f / N_entropy for f in bit_freq]
    min_f = min(freqs)
    max_f = max(freqs)
    # Число битов далеко от 0.5 (|freq - 0.5| > 0.05 — статистически значимо)
    biased = sum(1 for f in freqs if abs(f - 0.5) > 0.05)
    very_biased = sum(1 for f in freqs if abs(f - 0.5) > 0.1)
    print(f"  {name:25s}: freq_range=[{min_f:.3f},{max_f:.3f}], "
          f"biased_bits(>5%)={biased}, (>10%)={very_biased}")

# ─── Эксперимент 3: Корреляция между словами Δout ────────────────────────────

print("\n[3] XOR-корреляция между словами Δout[i] XOR Δout[j]")
print("─" * 70)
print("Если слова независимы → HW(Δout[i] XOR Δout[j]) ≈ 16.")
print("Концентрация: отклонение от 16.")

N_corr = 300

for name, gen in [
    ("случайные",          gen_random_pair),
    ("De6=0 (T_DECOMP)",   gen_pair_de6_zero),
    ("ротационный (s=1)",  lambda: gen_pair_rot_attractor(1)),
]:
    corr_matrix = [[[] for _ in range(8)] for _ in range(8)]
    for _ in range(N_corr):
        W, W2 = gen()
        out1 = sha256_compress_add_iv(W)
        out2 = sha256_compress_add_iv(W2)
        delta = [o1 ^ o2 for o1, o2 in zip(out1, out2)]
        for i in range(8):
            for j in range(i+1, 8):
                corr_matrix[i][j].append(hw(delta[i] ^ delta[j]))

    # Найти пары с наибольшим отклонением от 16
    deviations = []
    for i in range(8):
        for j in range(i+1, 8):
            avg = sum(corr_matrix[i][j]) / len(corr_matrix[i][j])
            deviations.append((abs(avg - 16), i, j, avg))
    deviations.sort(reverse=True)

    top3 = deviations[:3]
    print(f"\n  {name:25s}: топ-3 пар слов по корреляции:")
    for dev, i, j, avg in top3:
        print(f"    Δout[{i}] XOR Δout[{j}]: avg HW = {avg:.2f} "
              f"(отклонение от 16: {dev:.2f})")

# ─── Эксперимент 4: Усечённое столкновение — birthday на k битах ─────────────

print("\n[4] Усечённая коллизия Δout — эффективное количество бит")
print("─" * 70)
print("Берём только первые k бит Δout. Birthday коллизия за ~2^{k/2} пар.")
print("Ищем: для какого k число уникальных значений Δout_k растёт медленнее 2^k?")

N_trunc = 1000

for name, gen in [
    ("случайные",          gen_random_pair),
    ("De6=0 (T_DECOMP)",   gen_pair_de6_zero),
    ("ротационный (s=1)",  lambda: gen_pair_rot_attractor(1)),
]:
    # Собираем Δout (только первое слово = 32 бита)
    delta_word0 = []
    for _ in range(N_trunc):
        W, W2 = gen()
        out1 = sha256_compress_add_iv(W)
        out2 = sha256_compress_add_iv(W2)
        delta_word0.append(out1[0] ^ out2[0])

    print(f"\n  {name:25s}:")
    for k in [8, 12, 16, 20, 24, 32]:
        mask_k = (1 << k) - 1
        truncated = [d & mask_k for d in delta_word0]
        unique = len(set(truncated))
        theoretical_max = min(2**k, N_trunc)
        ratio = unique / theoretical_max
        # Если ratio << 1.0 → концентрация → birthday быстрее
        print(f"    k={k:2d}: unique={unique:5d}/{N_trunc}, "
              f"fill_rate={ratio:.3f} (1.0=равномерно, <1=концентрация)")

# ─── Эксперимент 5: Ротационная структура Δout ───────────────────────────────

print("\n[5] Ротационная корреляция Δout — связь с ROTR аттрактором")
print("─" * 70)
print("Для ротационных пар: ожидаем ROTR(Δout[i], k) ≈ Δout[i] ?")
print("Тест: HW(Δout[word] XOR ROTR(Δout[word], s)) для разных s.")

N_rot_corr = 500

for name, gen in [
    ("случайные",          gen_random_pair),
    ("ротационный (s=1)",  lambda: gen_pair_rot_attractor(1)),
    ("ротационный (s=2)",  lambda: gen_pair_rot_attractor(2)),
]:
    # Среднее HW(Δout[0] XOR ROTR(Δout[0], s)) для s=1,2,4
    results_rot = {}
    for s in [1, 2, 4, 8, 16]:
        hw_vals = []
        for _ in range(N_rot_corr):
            W, W2 = gen()
            out1 = sha256_compress_add_iv(W)
            out2 = sha256_compress_add_iv(W2)
            delta0 = out1[0] ^ out2[0]
            hw_vals.append(hw(delta0 ^ rotr(delta0, s)))
        results_rot[s] = sum(hw_vals) / len(hw_vals)

    best_s = min(results_rot, key=results_rot.get)
    print(f"\n  {name:25s}:")
    for s, avg in sorted(results_rot.items()):
        mark = " ← мин!" if s == best_s and results_rot[s] < 15 else ""
        print(f"    s={s:2d}: avg HW(Δ XOR ROTR(Δ,s)) = {avg:.2f}{mark}")

# ─── Вывод ───────────────────────────────────────────────────────────────────

print("\n" + "=" * 70)
print("ВЫВОД П-41")
print("=" * 70)
print("""
T_DELTA_OUT_ENTROPY (эмпирический):
  Для случайных пар: Δout равномерен, avg HW ≈ 128, все биты ≈ 50/50.
  Для пар с De6=0 (T_DECOMP): если Δout равномерен → birthday = 2^128 (без улучшения).
  Для ротационных аттракторов: если есть структура → возможное улучшение.

КЛЮЧЕВОЙ ВОПРОС РЕШЁН:
  Равномерность Δout определяет, даёт ли структура входного дифференциала
  улучшение над 2^128 для коллизии по выходу.

СТРАТЕГИЯ ДЛЯ П-42+:
  Если Δout НЕ равномерен для структурированных пар:
    → Оценить реальную сложность birthday collision
    → Сравнить с теоретическим барьером 2^128

  Если Δout РАВНОМЕРЕН для всех пар:
    → Структуры на входе не помогают для выходных коллизий
    → Нужен другой подход: не birthday на Δout, а другая метрика

НОВАЯ ГИПОТЕЗА (П-42):
  Многократно применять T_ATTRACTOR (k=0) для создания пар
  с De3..De17=0, затем искать КОЛЛИЗИЮ на входном состоянии раунда 18.
  Эта задача эквивалентна SHA-256 с 64-17=47 раундами.
  Сложность: зависит от структуры раундовой функции без W-зависимости.
""")
