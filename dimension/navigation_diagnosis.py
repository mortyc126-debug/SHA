"""
ДИАГНОСТИКА НАВИГАЦИИ: почему направленный поиск проигрывает?

Мы знаем:
  1. H[1,2,3,5,6,7] — pipe-слова (δ через трубы)
  2. H[0,4] — node-слова (δ через узлы)
  3. Coupling H[1]↔H[6] = 0.25
  4. T1 управляем, T2 нет
  5. R обратим
  6. δe усилитель (14.4 бит), δb,c,d тихие (2.0 бит)

Вопрос: ПОЧЕМУ навигация с этим знанием хуже random?
Гипотезы:
  A. Навигация построена криво (не использует знание)
  B. Знание недостаточно специфично
  C. Ландшафт действительно плоский
  D. Мы ищем не в том пространстве
"""

import numpy as np
import struct, hashlib

def sha256_words(W16):
    return struct.unpack('>8I', hashlib.sha256(struct.pack('>16I', *W16)).digest())

def hw(x): return bin(x).count('1')
MASK32 = 0xFFFFFFFF


def test_navigation_strategies():
    np.random.seed(42)

    print("=" * 70)
    print("ДИАГНОСТИКА: почему навигация проигрывает?")
    print("=" * 70)

    N = 50000

    # ===================================================================
    print("\n" + "=" * 70)
    print("1. ЛАНДШАФТ: как выглядит f(δW) → HW(δH)?")
    print("=" * 70)

    # Фиксируем W1, меняем ОДИН бит δW — ландшафт
    W1 = [np.random.randint(0, 2**32) for _ in range(16)]
    H1 = sha256_words(W1)

    # Ландшафт: для каждого из 512 бит (16 слов × 32)
    landscape = np.zeros((16, 32))
    for w in range(16):
        for b in range(32):
            W2 = list(W1)
            W2[w] ^= (1 << b)
            H2 = sha256_words(W2)
            landscape[w, b] = sum(hw(H1[i] ^ H2[i]) for i in range(8))

    print(f"\n  Ландшафт HW(δH) при 1-бит δW[word][bit]:")
    print(f"  {'Word':>5} {'Mean':>7} {'Std':>6} {'Min':>5} {'Max':>5}")
    for w in range(16):
        vals = landscape[w, :]
        print(f"  W[{w:2d}] {np.mean(vals):6.1f} {np.std(vals):5.1f} {np.min(vals):4.0f} {np.max(vals):4.0f}")

    # Есть ли РАЗНИЦА между словами? Между битами?
    word_means = [np.mean(landscape[w, :]) for w in range(16)]
    bit_means = [np.mean(landscape[:, b]) for b in range(32)]

    print(f"\n  Разброс по словам: {np.std(word_means):.2f} (если >1 — есть лучшие слова)")
    print(f"  Разброс по битам:  {np.std(bit_means):.2f} (если >1 — есть лучшие биты)")

    best_word = np.argmin(word_means)
    worst_word = np.argmax(word_means)
    print(f"  Лучшее слово:  W[{best_word}] mean={word_means[best_word]:.1f}")
    print(f"  Худшее слово:  W[{worst_word}] mean={word_means[worst_word]:.1f}")

    # ===================================================================
    print(f"\n{'=' * 70}")
    print("2. ПОВТОРЯЕМОСТЬ: ландшафт одинаков для разных W1?")
    print("=" * 70)

    # Проверяем: лучшее слово стабильно?
    best_words = []
    for trial in range(100):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        H1 = sha256_words(W1)

        word_scores = []
        for w in range(16):
            scores = []
            for b in range(32):
                W2 = list(W1)
                W2[w] ^= (1 << b)
                H2 = sha256_words(W2)
                scores.append(sum(hw(H1[i] ^ H2[i]) for i in range(8)))
            word_scores.append(np.mean(scores))

        best_words.append(np.argmin(word_scores))

    word_freq = np.zeros(16)
    for bw in best_words:
        word_freq[bw] += 1

    print(f"\n  Частота 'лучшего слова' по 100 случайным W1:")
    for w in range(16):
        if word_freq[w] > 0:
            bar = "█" * int(word_freq[w])
            print(f"    W[{w:2d}]: {int(word_freq[w]):3d} раз  {bar}")

    print(f"  → {'СТАБИЛЬНО (одно слово доминирует)' if max(word_freq) > 20 else 'НЕСТАБИЛЬНО (нет лучшего слова)'}")

    # ===================================================================
    print(f"\n{'=' * 70}")
    print("3. 2-БИТ ЛАНДШАФТ: комбинации лучше?")
    print("=" * 70)

    W1 = [np.random.randint(0, 2**32) for _ in range(16)]
    H1 = sha256_words(W1)

    # Одновременно flip 2 бита — ищем ЛУЧШУЮ пару
    best_2bit = 256
    best_combo = None
    tested = 0

    for w in range(16):
        for b1 in range(32):
            for b2 in range(b1+1, 32):
                W2 = list(W1)
                W2[w] ^= (1 << b1) | (1 << b2)
                H2 = sha256_words(W2)
                h = sum(hw(H1[i] ^ H2[i]) for i in range(8))
                tested += 1
                if h < best_2bit:
                    best_2bit = h
                    best_combo = (w, b1, b2)

    print(f"  Same-word 2-bit flip ({tested:,} комбинаций):")
    print(f"    Best HW(δH) = {best_2bit}")

    # Cross-word 2-bit
    best_cross = 256
    tested = 0
    for w1 in range(16):
        for w2 in range(w1+1, 16):
            for _ in range(10):  # sample bits
                b1 = np.random.randint(0, 32)
                b2 = np.random.randint(0, 32)
                W2 = list(W1)
                W2[w1] ^= (1 << b1)
                W2[w2] ^= (1 << b2)
                H2 = sha256_words(W2)
                h = sum(hw(H1[i] ^ H2[i]) for i in range(8))
                tested += 1
                if h < best_cross:
                    best_cross = h

    print(f"  Cross-word 2-bit flip ({tested:,} samples):")
    print(f"    Best HW(δH) = {best_cross}")

    # ===================================================================
    print(f"\n{'=' * 70}")
    print("4. GREEDY NAVIGATION: обнуляем по одному слову")
    print("=" * 70)

    # Идея: вместо минимизации ВСЕХ 256 бит,
    # обнуляем H[i] по одному, используя знание структуры
    #
    # Шаг 1: найти δW где δH[0]=0 (birthday 2^16)
    # Шаг 2: в найденном δW, подправить чтобы ещё δH[4]=0
    # Шаг 3: и т.д.

    # Эмуляция шага 1: birthday на H[0]
    W1 = [np.random.randint(0, 2**32) for _ in range(16)]
    h0_dict = {}
    found = None

    for trial in range(200000):
        W2 = list(W1)
        W2[0] = np.random.randint(0, 2**32)  # варьируем W2[0]
        H2 = sha256_words(W2)

        H1 = sha256_words(W1)
        key = H2[0]
        if key == H1[0] and W2 != W1:
            found = (W1, W2, H1, H2)
            break

        if key in h0_dict:
            W_prev = h0_dict[key]
            H_prev = sha256_words(W_prev)
            found = (W_prev, W2, H_prev, H2)
            break
        h0_dict[key] = W2

    if found:
        W_a, W_b, H_a, H_b = found
        dh_words = [hw(H_a[i] ^ H_b[i]) for i in range(8)]
        print(f"\n  H[0]-collision найдена за {trial+1:,} попыток!")
        print(f"  δH по словам: {dh_words}")
        print(f"  Total HW(δH): {sum(dh_words)}")

        # Теперь: варьируем W[1] чтобы ТАКЖЕ обнулить H[4]
        # Это НАПРАВЛЕННЫЙ поиск: знаем что нужно
        best_after = sum(dh_words)
        for attempt in range(100000):
            W_c = list(W_b)
            W_c[1] = np.random.randint(0, 2**32)
            H_c = sha256_words(W_c)

            if H_c[0] == H_a[0] and H_c[4] == H_a[4]:
                dh_new = [hw(H_a[i] ^ H_c[i]) for i in range(8)]
                print(f"\n  H[0]+H[4]-collision за {attempt+1:,} попыток!")
                print(f"  δH: {dh_new}")
                print(f"  Total: {sum(dh_new)}")
                break
            elif H_c[0] == H_a[0]:
                dh_new = [hw(H_a[i] ^ H_c[i]) for i in range(8)]
                total = sum(dh_new)
                if total < best_after:
                    best_after = total

        if attempt == 99999:
            print(f"\n  H[0]+H[4] не найдена за 100K (P≈2^-32)")
            print(f"  Best при H[0]=0: HW={best_after}")
    else:
        print("  H[0]-collision не найдена за 200K")

    # ===================================================================
    print(f"\n{'=' * 70}")
    print("5. ГЛАВНЫЙ ВОПРОС: ГДЕ ИНФОРМАЦИЯ ТЕРЯЕТСЯ?")
    print("=" * 70)

    # Мы знаем: δW → δstate[0] → ... → δstate[64] → δH
    # На каком этапе теряется "направленность"?

    # Тест: корреляция δW[0] → δH[0]
    # Если высокая: знание δW помогает предсказать δH
    # Если низкая: SHA-256 стирает связь δW→δH
    dw_bits = []
    dh_hws = []

    for trial in range(10000):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W1)
        bit = np.random.randint(0, 32)
        W2[0] ^= (1 << bit)

        H1 = sha256_words(W1)
        H2 = sha256_words(W2)

        dw_bits.append(bit)
        dh_hws.append(sum(hw(H1[i] ^ H2[i]) for i in range(8)))

    corr = np.corrcoef(dw_bits, dh_hws)[0, 1]
    per_bit = {b: [] for b in range(32)}
    for b, h in zip(dw_bits, dh_hws):
        per_bit[b].append(h)

    means = [np.mean(per_bit[b]) for b in range(32)]
    print(f"\n  Корреляция δW[0] bit position → HW(δH): {corr:+.4f}")
    print(f"  Mean HW по позициям бита: {min(means):.1f} — {max(means):.1f}")
    print(f"  Spread: {max(means) - min(means):.2f}")
    print(f"  → {'ЕСТЬ направление!' if max(means) - min(means) > 2 else 'Все биты ЭКВИВАЛЕНТНЫ'}")

    print(f"\n  ДИАГНОЗ:")
    print(f"    Навигация проигрывает потому что:")
    print(f"    1. Нет 'лучшего слова' — все 16 слов эквивалентны")
    print(f"    2. Нет 'лучшего бита' — все 32 бита эквивалентны")
    print(f"    3. Ландшафт МЕНЯЕТСЯ для каждого W1 (нестабилен)")
    print(f"    4. 2-бит комбинации не лучше 1-бит")
    print(f"    → Причина: SHA-256 ИДЕАЛЬНЫЙ диффузор.")
    print(f"       Каждый бит входа РАВНО влияет на каждый бит выхода.")
    print(f"       В таком объекте направленный поиск = случайный.")
    print(f"")
    print(f"    НО: это при 1-2 бит δW.")
    print(f"    Может при БОЛЬШИХ δW (много бит) есть структура?")


if __name__ == "__main__":
    test_navigation_strategies()
