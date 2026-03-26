"""
ФИНАЛЬНАЯ СТОИМОСТЬ COLLISION через трубный каскад.

Из pipe_cascade.py:
  6 слов (H[1,2,3,5,6,7]) → бесплатные через трубы
  2 слова (H[0,4]) → уравнения:
    δH[0] = δe[60] + δW[63] + carry
    δH[4] = δa[60] + δe[60] + δW[63] + carry
    → δH[4] - δH[0] = δa[60] + carry

Вопросы:
  1. H[0] и H[4] ЗАВИСИМЫ? (δH[4]-δH[0] = δa[60])
  2. Можно ли сделать birthday на 7 слов + использовать зависимость?
  3. Точная стоимость с учётом трубной структуры
"""

import numpy as np
import struct
import hashlib

MASK32 = 0xFFFFFFFF

def sha256_words(W16):
    h = hashlib.sha256(struct.pack('>16I', *W16)).digest()
    return struct.unpack('>8I', h)

def hw(x): return bin(x).count('1')


def experiment():
    np.random.seed(42)

    print("=" * 70)
    print("ФИНАЛЬНАЯ СТОИМОСТЬ: трубный каскад + зависимость H[0]↔H[4]")
    print("=" * 70)

    # =================================================================
    print("\n" + "=" * 70)
    print("1. ЗАВИСИМОСТЬ H[0] и H[4]: δH[4] ⊕ δH[0] = ?")
    print("=" * 70)

    # Теория: δH[4] - δH[0] = δa[60] (mod 2^32, + carry terms)
    # Проверяем: при случайных парах, corr(δH[0], δH[4])?

    dh0_list = []
    dh4_list = []
    dh04_xor_list = []  # HW(δH[0] ⊕ δH[4])

    for _ in range(100000):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W1)
        W2[0] ^= (1 << np.random.randint(0, 32))

        H1 = sha256_words(W1)
        H2 = sha256_words(W2)

        dh0 = hw(H1[0] ^ H2[0])
        dh4 = hw(H1[4] ^ H2[4])
        dh04 = hw((H1[0] ^ H2[0]) ^ (H1[4] ^ H2[4]))  # XOR of two deltas

        dh0_list.append(dh0)
        dh4_list.append(dh4)
        dh04_xor_list.append(dh04)

    corr_04 = np.corrcoef(dh0_list, dh4_list)[0, 1]

    print(f"\n  100K случайных пар:")
    print(f"    Mean HW(δH[0]):         {np.mean(dh0_list):.2f}")
    print(f"    Mean HW(δH[4]):         {np.mean(dh4_list):.2f}")
    print(f"    Mean HW(δH[0]⊕δH[4]):  {np.mean(dh04_xor_list):.2f}")
    print(f"    corr(δH[0], δH[4]):     {corr_04:+.4f}")
    print(f"    Если независимы: corr ≈ 0, HW(XOR) ≈ 16")

    # =================================================================
    print("\n" + "=" * 70)
    print("2. УСЛОВНАЯ ВЕРОЯТНОСТЬ: P(δH[4]=0 | δH[0]=0)")
    print("=" * 70)

    # Ищем пары где δH[0]=0 (birthday на 1 слово)
    h0_dict = {}
    h0_collisions = []

    N = 200000
    for i in range(N):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H = sha256_words(W)
        h0 = H[0]
        if h0 in h0_dict:
            h0_collisions.append((H, h0_dict[h0]))
        h0_dict[h0] = H

    print(f"\n  {N:,} хешей → {len(h0_collisions)} H[0]-коллизий")

    if h0_collisions:
        # При H[0]-collision: каково δH[4]?
        dh4_at_h0coll = [hw(h1[4] ^ h2[4]) for h1, h2 in h0_collisions]
        # И другие слова для сравнения
        dh1_at_h0coll = [hw(h1[1] ^ h2[1]) for h1, h2 in h0_collisions]
        dh7_at_h0coll = [hw(h1[7] ^ h2[7]) for h1, h2 in h0_collisions]
        dh4_zero = sum(1 for d in dh4_at_h0coll if d == 0)

        print(f"\n  При δH[0]=0:")
        print(f"    Mean HW(δH[4]): {np.mean(dh4_at_h0coll):.2f}  (random: 16)")
        print(f"    Mean HW(δH[1]): {np.mean(dh1_at_h0coll):.2f}  (random: 16)")
        print(f"    Mean HW(δH[7]): {np.mean(dh7_at_h0coll):.2f}  (random: 16)")
        print(f"    P(δH[4]=0 | δH[0]=0): {dh4_zero}/{len(h0_collisions)}")
        print(f"    Если независимы: P ≈ 2^-32 ≈ {1/2**32:.2e}")

        if np.mean(dh4_at_h0coll) < 15.0:
            compression = 16 - np.mean(dh4_at_h0coll)
            print(f"    >>> СЖАТИЕ: {compression:.1f} бит! H[4] ЗАВИСИТ от H[0]!")
        else:
            print(f"    >>> Нет сжатия: H[4] НЕЗАВИСИМ от H[0]")

    # =================================================================
    print("\n" + "=" * 70)
    print("3. МНОГОСЛОВНЫЙ BIRTHDAY: стоимость по комбинациям")
    print("=" * 70)

    # Для каждой комбинации k слов: birthday cost = 2^(k*16)
    combos = [
        ("H[7] (1 word = 32 bit)",              1, 32),
        ("H[7,6] (2 words = 64 bit)",           2, 64),
        ("H[7,6,5] (3 words = 96 bit)",         3, 96),
        ("H[7,6,5,3] (4 words = 128 bit)",      4, 128),
        ("H[7,6,5,3,2] (5 words = 160 bit)",    5, 160),
        ("H[7,6,5,3,2,1] (6 pipe = 192 bit)",   6, 192),
        ("H[7,6,5,3,2,1,0] (7 = 224 bit)",      7, 224),
        ("H[0..7] (full = 256 bit)",             8, 256),
    ]

    print(f"\n  {'Комбинация':<42} {'Бит':>5} {'Birthday cost':>15} {'Остаток':>10}")
    print(f"  {'-'*42} {'-'*5} {'-'*15} {'-'*10}")

    for desc, n_words, bits in combos:
        cost = bits // 2
        remaining = 256 - bits
        marker = ""
        if n_words == 6:
            marker = " ← ТРУБНЫЙ ОПТИМУМ"
        if n_words == 8:
            marker = " ← ПОЛНАЯ COLLISION"

        print(f"  {desc:<42} {bits:>4}  2^{cost:<12d} {remaining:>5} бит{marker}")

    # =================================================================
    print("\n" + "=" * 70)
    print("4. СТОИМОСТЬ ТРУБНОГО КАСКАДА — ТОЧНАЯ")
    print("=" * 70)

    print(f"""
  ПОДХОД 1: Прямой birthday (стандарт)
    256 бит выхода → 2^128 пар → 2^128 SHA-256

  ПОДХОД 2: Трубный каскад (наш)
    Шаг 1: Birthday на H[1,2,3,5,6,7] = 192 бит
            Стоимость: 2^96 SHA-256 + O(2^96) память
            Результат: пара (W1,W2) с δH[1,2,3,5,6,7]=0

    Шаг 2: Проверка δH[0]=0 и δH[4]=0
            P(δH[0]=0 | 6-word match) = 2^-32
            P(δH[4]=0 | 6-word match, δH[0]=0):""")

    # Вычисляем: сколько 6-word collisions нужно чтобы одна стала 8-word
    # Из наших уравнений: δH[4]-δH[0] ≈ δa[60]
    # Если δH[0]=0: δH[4] ≈ δa[60] ≈ random 32-bit
    # Значит P(δH[4]=0 | δH[0]=0, 6-word) ≈ 2^-32

    # Но если H[0] и H[4] коррелированы (corr ≠ 0):
    conditional_p = 2**-32  # baseline
    if len(h0_collisions) > 0:
        if np.mean(dh4_at_h0coll) < 15.5:
            # Some compression — adjust probability
            effective_bits = np.mean(dh4_at_h0coll) * 2  # rough estimate
            conditional_p = 2**(-effective_bits)

    print(f"              P ≈ 2^-32 (из эксперимента: δH[4] independent)")
    print(f"""
    Шаг 3: Повторять Шаг 1 пока Шаг 2 не сработает
            Нужно: 2^32 × 2^32 = 2^64 проверок (для обоих H[0] и H[4])

    ИТОГО:
      Вариант A: 2^64 birthday поисков × 2^96 каждый = 2^160
      → ХУЖЕ стандарта!

      Вариант B: Один birthday на ВСЕ 8 слов = 2^128
      → Стандарт

      Вариант C: Birthday на 7 слов (включая H[0]) = 224 бит
                 Cost: 2^112
                 P(δH[4]=0 | 7-word): 2^-32
                 Нужно: 2^32 seven-word collisions
                 Каждая стоит 2^112
                 → 2^32 × 2^112 = 2^144 ХУЖЕ

      Вариант D: Один ОГРОМНЫЙ birthday
                 N хешей в таблице
                 6-word коллизий: N²/2^193
                 Полных: N²/2^193 × 2^-64 = N²/2^257
                 Для 1 полной: N = 2^128.5
                 → Стандарт""")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("5. ПОЧЕМУ ТРУБЫ НЕ СНИЖАЮТ BIRTHDAY BOUND")
    print("=" * 70)

    print(f"""
  ПРИЧИНА: 192 + 64 = 256

  Трубный каскад разбивает 256 бит на:
    192 бит (трубных, κ=1) + 64 бит (узловых, κ≈0.05)

  Но birthday bound зависит от СУММАРНОГО числа бит, не от разбивки.
  Неважно что 192 бита "бесплатны" — их всё равно нужно обнулить.

  АНАЛОГИЯ: замок с 8 цилиндрами.
  6 цилиндров — простые (трубы). 2 — сложные (узлы).
  Но открыть замок = совместить ВСЕ 8.
  Простота 6-ти не помогает если нужны все одновременно.

  ЧТО ТРУБЫ РЕАЛЬНО ДАЮТ:
    ✓ ПОНИМАНИЕ: почему H[7]↔e[61], H[3]↔a[61]
    ✓ УРАВНЕНИЯ: δH[0] и δH[4] через δa[60], δe[60], δW[63]
    ✓ INSIGHT: Ch и Maj ВЫКЛЮЧАЮТСЯ при 6-word collision
    ✗ УСКОРЕНИЕ: нет (192+64=256, birthday = 2^128)

  ЕДИНСТВЕННЫЙ ПУТЬ К УСКОРЕНИЮ:
    Найти ЗАВИСИМОСТЬ между трубными и узловыми битами.
    Наш эксперимент: corr(δH[0], δH[4]) = {corr_04:+.4f}
    → {'ЕСТЬ ЗАВИСИМОСТЬ!' if abs(corr_04) > 0.05 else 'НЕТ зависимости.'}
    → Стоимость = 2^128 (стандарт)
""")


if __name__ == "__main__":
    experiment()
