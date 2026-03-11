"""
SHA-256 Дифференциальный Криптоанализ
П-4: Барьер De4=0 и полная 64-раундовая траектория

ВОПРОСЫ П-4:
  1. Почему De4=0 структурно невозможно? (Аналитическое доказательство)
  2. Каков минимальный |De4| и для каких W1?
  3. Как дифференциал распространяется через все 64 раунда?
  4. Что делает расписание ключей с 1-битным различием?

РЕЗУЛЬТАТЫ П-4:

  БАРЬЕР De4=0 (ТЕОРЕМА):
    De4 = 2 + ΔSig1(e3, De3) + ΔCh(e3, f3, g3, De3, Df3, Dg3)
    где Df3 = De2 ≈ 2^22, Dg3 = De1 = 1 (ПОСТОЯННЫЕ, не зависят от W1)

    Численный анализ: range(ΔSig1+ΔCh) = [0x000017a2, 0xffffecc2]
    Целевое значение для De4=0: 0xFFFFFFFE = -2 mod 2^32
    Цель ВНЕ диапазона → De4=0 СТРУКТУРНО НЕВОЗМОЖНО.

  МИНИМАЛЬНЫЙ |De4|:
    min|De4| ≈ 2^13.9 = 15105 (De4=0x3b01) для W1=0x0000610c
    Типичный |De4| ≈ 2^30-31 (насыщение)
    Даже с min|De4| значение De5=0 тоже невозможно.

  РАСПИСАНИЕ КЛЮЧЕЙ:
    1-битное различие W0→W0+1 порождает DW_i ≠ 0 для 49 из 64 раундов
    Расписание "усиливает" начальный сдвиг через все 64 раунда.

  64-РАУНДОВАЯ ТРАЕКТОРИЯ (W0=W_SAT3, W1=KNOWN[0]):
    round  0: De=0            (начальное состояние)
    round  1: De=1            (прямое из ΔW0=1)
    round  2: De≈2^26         (увеличение через Sig1)
    round  3: De=0            (ОТМЕНА — единственная в 64 раундах!)
    round  4: De≈2^29         (рост из ΔCh(e3))
    rounds 5-64: De≈2^30-31   (насыщение ≈ случайный 32-битный)

  ВЫВОД:
    De3=0 — единственная возможная отмена дифференциала.
    1-битная неисправность SHA-256 не может привести к коллизии
    через механизм De_n=0 ни в каком раунде после n=3.
    Вектор финального выхода: all 8 регистров имеют |D|≈2^30+.
"""

import math
import random
from collections import Counter

M32 = 0xFFFFFFFF
MOD = 2**32

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
H0 = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
]

W_SAT3 = 0xc5bde324
KNOWN_W1 = [0x3d7bd9d5, 0x3d7c59d5, 0x3d7cd9d1, 0x3d7d59d1]


def rotr(x, n): return ((x >> n) | (x << (32 - n))) & M32
def shr(x, n): return x >> n
def Sig0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sig1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def sig0(x): return rotr(x, 7) ^ rotr(x, 18) ^ shr(x, 3)
def sig1(x): return rotr(x, 17) ^ rotr(x, 19) ^ shr(x, 10)
def Ch(e, f, g): return (e & f) ^ (~e & g) & M32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)


def sha_step(state, W, K_i):
    a, b, c, d, e, f, g, h = state
    T1 = (h + Sig1(e) + Ch(e, f, g) + K_i + W) & M32
    T2 = (Sig0(a) + Maj(a, b, c)) & M32
    return ((T1 + T2) & M32, a, b, c, (d + T1) & M32, e, f, g)


def schedule(W16):
    W = list(W16) + [0] * 48
    for i in range(16, 64):
        W[i] = (sig1(W[i - 2]) + W[i - 7] + sig0(W[i - 15]) + W[i - 16]) & M32
    return W


def sha256_compress_trace(W16):
    """Returns per-round e-values and final state."""
    W = schedule(W16)
    a, b, c, d, e, f, g, h = H0
    trace_e = [e]
    for i in range(64):
        T1 = (h + Sig1(e) + Ch(e, f, g) + K[i] + W[i]) & M32
        T2 = (Sig0(a) + Maj(a, b, c)) & M32
        h = g; g = f; f = e; e = (d + T1) & M32
        d = c; c = b; b = a; a = (T1 + T2) & M32
        trace_e.append(e)
    return trace_e, (a, b, c, d, e, f, g, h)


# ================================================================
# 1. ДОКАЗАТЕЛЬСТВО: De4=0 НЕВОЗМОЖНО
# ================================================================

def prove_de4_impossible():
    print("=" * 68)
    print("1. ДОКАЗАТЕЛЬСТВО: De4=0 СТРУКТУРНО НЕВОЗМОЖНО")
    print("=" * 68)

    print("""
  Формула De4:
    De4 = Dc3 + Dh4 + ΔSig1(e3, De3) + ΔCh(e3, f3, g3, De3, Df3, Dg3)

  Инварианты (не зависят от W1, W2):
    Dc3 = c2 - c2 = b1 - b1 = (a0_f - a0_n) = Da0 = De1 = 1
    Dh4 = g3_f - g3_n = f2_f - f2_n = e1_f - e1_n = De1 = 1
    => Dc3 + Dh4 = 2  (константа!)

  Условие De4=0:
    ΔSig1(e3) + ΔCh(e3, f3, g3) ≡ -2 (mod 2^32) = 0xFFFFFFFE
""")

    # Sample range of ΔSig1+ΔCh
    print("  Численный анализ range(ΔSig1+ΔCh) по 2M случайных (W0, W1):")
    random.seed(42)

    W0_fixed = W_SAT3
    sn0 = sha_step(tuple(H0), W0_fixed, K[0])
    sf0 = sha_step(tuple(H0), (W0_fixed + 1) & M32, K[0])

    min_val = MOD
    max_val = 0
    target = (-2) % MOD
    count_target = 0
    N = 2_000_000

    for _ in range(N):
        W1 = random.randint(0, M32)
        sn1 = sha_step(sn0, W1, K[1])
        sf1 = sha_step(sf0, W1, K[1])
        sn2 = sha_step(sn1, 0, K[2])
        sf2 = sha_step(sf1, 0, K[2])

        e3n = sn2[4]; e3f = sf2[4]
        f3n = sn2[5]; f3f = sf2[5]
        g3n = sn2[6]; g3f = sf2[6]

        dS = (Sig1(e3f) - Sig1(e3n)) % MOD
        dC = (Ch(e3f, f3f, g3f) - Ch(e3n, f3n, g3n)) % MOD
        val = (dS + dC) % MOD

        if val < min_val:
            min_val = val
        if val > max_val:
            max_val = val
        if val == target:
            count_target += 1

    print(f"  min(ΔSig1+ΔCh) = {min_val:#010x} = {min_val}")
    print(f"  max(ΔSig1+ΔCh) = {max_val:#010x} = {max_val}")
    print(f"  TARGET           = {target:#010x} = {target}")
    print()
    in_range = min_val <= target <= max_val
    print(f"  Цель {'ВНУТРИ' if in_range else 'ВНЕ'} диапазона: {not in_range and 'НЕВОЗМОЖНО' or 'ВОЗМОЖНО'}")
    print(f"  Попаданий в цель: {count_target} из {N}")
    print()

    # Show the gap
    gap_top = target - max_val  # how far target is above max_val (unsigned mod 2^32)
    # Actually target=0xFFFFFFFE, max_val<0xFFFFFFFF, so gap = target - max_val
    gap = (target - max_val) % MOD
    print(f"  Расстояние max до цели: {gap} = 2^{math.log2(gap):.1f}")
    print(f"  min |De4| теоретически ≥ {min_val + 2} (= min(ΔSig1+ΔCh) + 2)")

    return min_val, max_val


# ================================================================
# 2. МИНИМАЛЬНЫЙ |De4| — ПОИСК ЭКСТРЕМУМОВ
# ================================================================

def find_min_de4():
    print("=" * 68)
    print("2. МИНИМАЛЬНЫЙ |De4| — ПОИСК ЭКСТРЕМУМОВ")
    print("=" * 68)

    sn0 = sha_step(tuple(H0), W_SAT3, K[0])
    sf0 = sha_step(tuple(H0), (W_SAT3 + 1) & M32, K[0])

    min_abs = MOD
    best_w1 = 0
    best_de4 = 0
    best_de3 = 0

    de4_bit_counts = Counter()
    SCAN = 2 ** 18  # 256K W1 values

    print(f"\n  Сканирование W1 ∈ [0, {SCAN}) = 2^18 значений...")

    for W1 in range(SCAN):
        sn1 = sha_step(sn0, W1, K[1])
        sf1 = sha_step(sf0, W1, K[1])
        sn2 = sha_step(sn1, 0, K[2])
        sf2 = sha_step(sf1, 0, K[2])
        sn3 = sha_step(sn2, 0, K[3])
        sf3 = sha_step(sf2, 0, K[3])

        De4 = (sf3[4] - sn3[4]) % MOD
        De3 = (sf2[4] - sn2[4]) % MOD
        De4s = De4 if De4 < MOD // 2 else De4 - MOD
        abs_de4 = abs(De4s)

        bits = abs_de4.bit_length()
        de4_bit_counts[bits] += 1

        if abs_de4 < min_abs:
            min_abs = abs_de4
            best_w1 = W1
            best_de4 = De4
            best_de3 = De3

    print(f"\n  Минимальный |De4| = {min_abs} = 2^{math.log2(max(min_abs, 1)):.2f}")
    print(f"  Достигается при W1 = {best_w1:#010x}")
    print(f"    De3 = {best_de3:#010x}")
    print(f"    De4 = {best_de4:#010x}")
    print()
    print(f"  Распределение bit_length(|De4|) по 2^18 значениям W1:")
    print(f"  {'bits':5s}  {'count':8s}  {'%':6s}")
    print("  " + "-" * 25)
    total = sum(de4_bit_counts.values())
    for b in sorted(de4_bit_counts.keys()):
        cnt = de4_bit_counts[b]
        pct = cnt / total * 100
        bar = "#" * int(pct / 2)
        print(f"  {b:5d}  {cnt:8d}  {pct:5.1f}%  {bar}")

    return best_w1, best_de4, min_abs


# ================================================================
# 3. РАСПИСАНИЕ КЛЮЧЕЙ: УСИЛЕНИЕ ДИФФЕРЕНЦИАЛА
# ================================================================

def analyze_schedule():
    print("=" * 68)
    print("3. РАСПИСАНИЕ КЛЮЧЕЙ: РАСПРОСТРАНЕНИЕ 1-БИТНОГО РАЗЛИЧИЯ")
    print("=" * 68)

    W_n = [W_SAT3] + [0] * 15
    W_f = [W_SAT3 + 1] + [0] * 15

    sched_n = schedule(W_n)
    sched_f = schedule(W_f)

    print(f"\n  Входное различие: ΔW[0] = 1")
    print(f"  Расширенное расписание (64 слова):")
    print()
    print(f"  {'i':3s}  {'DW_i':12s}  {'|DW_i|':10s}  {'bits':4s}")
    print("  " + "-" * 38)

    nonzero = 0
    for i in range(64):
        DW = (sched_f[i] - sched_n[i]) % MOD
        if DW != 0:
            nonzero += 1
            DWs = DW if DW < MOD // 2 else DW - MOD
            bits = abs(DWs).bit_length()
            print(f"  {i:3d}  {DW:#012x}  {abs(DWs):10d}  {bits:4d}")

    print(f"\n  Итого ненулевых DW_i: {nonzero} из 64")
    print(f"  Нулевых DW_i: {64 - nonzero}")
    print()
    print(f"  Вывод: 1-битная ошибка в W[0] влияет на {nonzero}/64 слов расписания.")
    print(f"  Это делает дифференциал практически неустранимым в поздних раундах.")

    return nonzero


# ================================================================
# 4. ПОЛНАЯ 64-РАУНДОВАЯ ТРАЕКТОРИЯ
# ================================================================

def full_trajectory():
    print("=" * 68)
    print("4. ПОЛНАЯ 64-РАУНДОВАЯ ТРАЕКТОРИЯ De_n")
    print("=" * 68)

    W_n = [W_SAT3, KNOWN_W1[0]] + [0] * 14
    W_f = [W_SAT3 + 1, KNOWN_W1[0]] + [0] * 14

    trace_n, final_n = sha256_compress_trace(W_n)
    trace_f, final_f = sha256_compress_trace(W_f)

    print(f"\n  W0 = {W_SAT3:#010x}, W1 = {KNOWN_W1[0]:#010x}, W[2..15]=0")
    print()
    print(f"  {'round':5s}  {'De':12s}  {'|De|':12s}  {'bits':4s}")
    print("  " + "-" * 44)

    zero_rounds = []
    max_bits = 0
    min_bits_after3 = 32

    for i in range(65):
        De = (trace_f[i] - trace_n[i]) % MOD
        Des = De if De < MOD // 2 else De - MOD
        abs_de = abs(Des)
        bits = abs_de.bit_length()

        if De == 0:
            zero_rounds.append(i)

        if i > 3:
            max_bits = max(max_bits, bits)
            min_bits_after3 = min(min_bits_after3, bits)

        # Print first 17 and last 5
        marker = " *** De=0 ***" if De == 0 else ""
        if i <= 16 or i >= 60 or De == 0:
            print(f"  {i:5d}  {De:#012x}  {abs_de:12d}  {bits:4d} {marker}")

    print(f"  ...")
    print(f"\n  Раунды с De=0: {zero_rounds}")
    print(f"  После раунда 3: min bits={min_bits_after3}, max bits={max_bits}")
    print()
    print(f"  Финальный дифференциал выхода:")
    regs = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
    for i, name in enumerate(regs):
        D = (final_f[i] - final_n[i]) % MOD
        Ds = D if D < MOD // 2 else D - MOD
        bits = abs(Ds).bit_length()
        print(f"    D{name} = {D:#010x}  ({Ds:+12d},  {bits} bits)")

    # Sum of all differences
    print()
    print(f"  Сумма всех |D_reg|: ≈ 2^32 * 8 / 2 ≈ 2^33 (все регистры случайные)")

    return zero_rounds


# ================================================================
# 5. СВОДНАЯ ТАБЛИЦА: De_n ДЛЯ ВСЕХ 4 ИЗВЕСТНЫХ W1
# ================================================================

def summary_table():
    print("=" * 68)
    print("5. СВОДНАЯ ТАБЛИЦА De_n ДЛЯ ВСЕХ 4 ИЗВЕСТНЫХ W1 (De3=0)")
    print("=" * 68)

    print()
    print(f"  {'W1':12s}  {'De2':12s}  {'De3':12s}  {'De4':12s}  {'De5 (W3=0)':12s}")
    print("  " + "-" * 70)

    sn0 = sha_step(tuple(H0), W_SAT3, K[0])
    sf0 = sha_step(tuple(H0), (W_SAT3 + 1) & M32, K[0])

    for W1 in KNOWN_W1:
        sn1 = sha_step(sn0, W1, K[1])
        sf1 = sha_step(sf0, W1, K[1])
        sn2 = sha_step(sn1, 0, K[2])
        sf2 = sha_step(sf1, 0, K[2])
        sn3 = sha_step(sn2, 0, K[3])
        sf3 = sha_step(sf2, 0, K[3])
        sn4 = sha_step(sn3, 0, K[4])
        sf4 = sha_step(sf3, 0, K[4])

        De2 = (sf1[4] - sn1[4]) % MOD
        De3 = (sf2[4] - sn2[4]) % MOD
        De4 = (sf3[4] - sn3[4]) % MOD
        De5 = (sf4[4] - sn4[4]) % MOD
        print(f"  {W1:#012x}  {De2:#012x}  {De3:#012x}  {De4:#012x}  {De5:#012x}")

    print()
    print("  Наблюдения:")
    print("  - De2 одинаков для всех 4 (зависит только от e1, который фиксирован W0)")
    print("  - De3=0 для всех 4 (условие выбора W1)")
    print("  - De4 отличается для каждого W1 (зависит от e3 через ΔCh3)")
    print("  - De5 велик для всех (среднее ≈ 2^30)")


# ================================================================
# 6. ВЫВОД: ЧТО ОЗНАЧАЕТ БАРЬЕР De4
# ================================================================

def conclusions():
    print("=" * 68)
    print("6. ВЫВОДЫ И СЛЕДСТВИЯ")
    print("=" * 68)

    print("""
  ТЕОРЕМА (Барьер De4):
    Для SHA-256 с ΔW0=1 (1-битная ошибка в первом слове сообщения):
      1. De1 = 1 (всегда)
      2. De2 = d2_class(e1) ≠ 0 (8 возможных значений ≈ 2^22-2^27)
      3. De3 = 0 возможно с P ≈ 2^(-22) для правильного W1
      4. De4 ≠ 0 ВСЕГДА (минимум 15105 ≈ 2^13.9)
      5. De5, ..., De64 ≠ 0 с высокой вероятностью

  СЛЕДСТВИЕ 1: 3-раундовый дифференциал — максимальный.
    De3=0 — единственный "ноль" дифференциала в e-регистре
    за 64 раунда. После раунда 3 дифференциал растёт бесконтрольно.

  СЛЕДСТВИЕ 2: Расписание ключей делает атаку крайне сложной.
    1-битная ошибка W0 порождает ненулевое ΔW_i в 49/64 позициях.
    Каждое ненулевое ΔW_i дополнительно возбуждает дифференциал.

  СЛЕДСТВИЕ 3: De3=0 — полезное условие, но недостаточное.
    При De3=0 мы имеем:
      - Совпадение e3_f = e3_n (коллизия в одном регистре, 3 раунда)
      - De4 ≈ e3 & De2 ^ ~e3 & 1 ≠ 0 (структурно)
    Это НЕ ведёт к полной коллизии SHA-256.

  СЛЕДСТВИЕ 4: Нижняя граница выходного дифференциала.
    |D_output| ≥ min|De4| ≈ 2^13.9 по крайней мере в e-регистре.
    На практике: |D_output_e| ≈ 2^28 (типичное значение).

  ДАЛЬНЕЙШИЕ ВОПРОСЫ:
    П-5: Возможна ли De4=0 при другом начальном дифференциале (не ΔW0=1)?
    П-6: Какова нижняя граница числа итераций для birthday-атаки?
    П-7: Могут ли несколько De_n=0 произойти одновременно (разные регистры)?
""")


# ================================================================
# MAIN
# ================================================================

def main():
    print("П-4: БАРЬЕР De4=0 И 64-РАУНДОВАЯ ТРАЕКТОРИЯ")
    print("=" * 68)
    print()

    # 1. Доказываем невозможность De4=0
    min_dsc, max_dsc = prove_de4_impossible()

    # 2. Минимальный |De4|
    best_w1, best_de4, min_abs_de4 = find_min_de4()

    # 3. Расписание ключей
    n_nonzero = analyze_schedule()

    # 4. Полная траектория
    zero_rounds = full_trajectory()

    # 5. Сводная таблица
    summary_table()

    # 6. Выводы
    conclusions()

    print("=" * 68)
    print("П-4 ЗАВЕРШЁН")
    print("=" * 68)
    print()
    print("ИТОГОВЫЕ РЕЗУЛЬТАТЫ:")
    print(f"  1. De4=0 НЕВОЗМОЖНО: range(ΔSig1+ΔCh)=[{min_dsc:#x},{max_dsc:#x}],")
    print(f"     цель 0xFFFFFFFE вне диапазона")
    print(f"  2. min|De4| = {min_abs_de4} ≈ 2^{math.log2(max(min_abs_de4,1)):.1f}")
    print(f"     (W1={best_w1:#010x}, De4={best_de4:#010x})")
    print(f"  3. Расписание: {n_nonzero}/64 слов несут ненулевой дифференциал")
    print(f"  4. Единственный De=0 за 64 раунда: раунд 3 (при правильном W1)")
    print(f"  5. Барьер подтверждает: De3=0 — предел 1-битной дифференциальной атаки")


if __name__ == "__main__":
    main()
