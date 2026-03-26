"""
Каскад H[7]→H[4]: измеряем +17 бит сжатия из методички (П-1253).

Теория:
  1. Найти H[7]-collision (birthday: ~2^16.6 ≈ 100K сэмплов)
  2. При H[7]-collision: H[4] сжимается (E[HW(ΔH[4])]=12 вместо 16)
  3. Каскад: +17 бит суммарного сжатия → 2^119.5 вместо 2^128

План:
  Этап 1: Набираем H[7]-collisions через birthday
  Этап 2: Измеряем ΔH по словам внутри H[7]-collision пар
  Этап 3: Оцениваем реальное сжатие
  Этап 4: Экстраполяция на полную collision
"""

import hashlib
import struct
import numpy as np
import time

def sha256_compress(W16):
    """SHA-256 compression function (single block)."""
    return hashlib.sha256(struct.pack('>16I', *W16)).digest()

def h_to_words(h):
    return struct.unpack('>8I', h)

def hw(x):
    return bin(x).count('1')


def find_h7_collisions(n_samples=500000):
    """Birthday search на H[7] (32-бит): ожидаем ~N²/2^33 коллизий."""
    np.random.seed(42)

    print(f"  Генерация {n_samples:,} хешей...")
    t0 = time.time()

    h7_dict = {}  # H[7] → (index, W, full_hash)
    collisions = []
    all_data = []

    for i in range(n_samples):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        h = sha256_compress(W)
        words = h_to_words(h)
        h7 = words[7]

        entry = (i, W, h, words)
        all_data.append(entry)

        if h7 in h7_dict:
            j, W_prev, h_prev, words_prev = h7_dict[h7]
            # Проверяем что это реальная H[7]-collision, не дубликат
            if W != W_prev:
                collisions.append((entry, h7_dict[h7]))
        else:
            h7_dict[h7] = entry

    elapsed = time.time() - t0
    expected = n_samples**2 / 2**33

    print(f"  Время: {elapsed:.1f}с")
    print(f"  H[7]-коллизий: {len(collisions)} (ожидалось ≈{expected:.1f})")

    return collisions, all_data


def analyze_collisions(collisions):
    """Анализируем ΔH по словам для каждой H[7]-collision."""
    if not collisions:
        return None

    per_word_hw = [[] for _ in range(8)]
    total_hw = []

    for (i1, W1, h1, w1), (i2, W2, h2, w2) in collisions:
        for k in range(8):
            d = hw(w1[k] ^ w2[k])
            per_word_hw[k].append(d)

        dh = sum(hw(w1[k] ^ w2[k]) for k in range(8))
        total_hw.append(dh)

    return per_word_hw, total_hw


def experiment():
    print("=" * 70)
    print("КАСКАД H[7]→H[4]: Birthday + Сжатие")
    print("=" * 70)

    # ===================================================================
    print("\n" + "=" * 70)
    print("1. BIRTHDAY SEARCH: H[7]-collisions")
    print("=" * 70)

    collisions, all_data = find_h7_collisions(n_samples=500000)

    if not collisions:
        print("  Нет коллизий. Увеличиваем выборку...")
        collisions, all_data = find_h7_collisions(n_samples=2000000)

    if not collisions:
        print("  Всё ещё нет. Нужно ≈2^16.6 ≈ 100K для ожидания 1 collision.")
        print("  При N=2M ожидалось ≈465. Что-то не так.")
        return

    # ===================================================================
    print("\n" + "=" * 70)
    print("2. АНАЛИЗ: ΔH по словам при H[7]-collision")
    print("=" * 70)

    per_word, total = analyze_collisions(collisions)

    reg = ['H[0](a)', 'H[1](b)', 'H[2](c)', 'H[3](d)',
           'H[4](e)', 'H[5](f)', 'H[6](g)', 'H[7](h)']

    print(f"\n  {len(collisions)} H[7]-коллизий:")
    print(f"\n  {'Слово':<10} {'Mean HW(Δ)':>12} {'Std':>8} {'Expected':>10} {'Compression':>13}")
    print(f"  {'-'*10} {'-'*12} {'-'*8} {'-'*10} {'-'*13}")

    total_compression = 0
    for k in range(8):
        m = np.mean(per_word[k])
        s = np.std(per_word[k])
        expected = 0 if k == 7 else 16
        comp = expected - m
        total_compression += max(comp, 0)

        marker = ""
        if k == 7:
            marker = " ← COLLISION"
        elif comp > 1:
            marker = f" ← +{comp:.1f} бит"

        print(f"  {reg[k]:<10} {m:11.2f} {s:7.2f} {expected:9d} {comp:+12.1f}{marker}")

    print(f"\n  Суммарное сжатие: +{total_compression:.1f} бит")
    print(f"  Методичка (П-1253): +17 бит")

    # Полный HW(ΔH) без H[7]
    remaining_hw = [sum(per_word[k][i] for k in range(7)) for i in range(len(collisions))]
    print(f"\n  HW(ΔH[0..6]) — без H[7] (7 слов × 32 = 224 бит):")
    print(f"    Mean: {np.mean(remaining_hw):.1f}")
    print(f"    Std:  {np.std(remaining_hw):.1f}")
    print(f"    Min:  {min(remaining_hw)}")
    print(f"    Max:  {max(remaining_hw)}")
    print(f"    Random expected: 112 (7 × 16)")

    # ===================================================================
    print("\n" + "=" * 70)
    print("3. H[4] ДЕТАЛЬНО: механизм сжатия")
    print("=" * 70)

    # H[4] = e[64] + IV[4]. H[7] = h[64] + IV[7].
    # h[63] = g[62] = f[61] = e[60] (трубы!)
    # H[7]-collision: h1[64] + IV[7] = h2[64] + IV[7] → h1[64] = h2[64]
    # h[64] = g[63] (труба) = f[62] (труба) = e[61] (труба)
    # Значит: e[61] одинаков!
    # e[64] = e[60+4] — но e[64] определяется раундом 63.
    # Связь H[4]↔H[7]: через e-chain трубы.

    h4_vals = per_word[4]
    h7_vals = per_word[7]  # все = 0

    print(f"\n  H[7]-collision → h1[64] = h2[64]")
    print(f"  h[64] = g[63] = f[62] = e[61] (3 трубы)")
    print(f"  → δe[61] = 0!")
    print(f"  → H[4] = e[64] + IV[4], и e[64] зависит от e[61]")
    print(f"  → δe[64] частично подавлен через δe[61]=0")

    print(f"\n  Распределение HW(ΔH[4]):")
    h4_hist = np.zeros(33)
    for v in h4_vals:
        h4_hist[v] += 1
    for hw_val in range(0, 33, 2):
        count = int(h4_hist[hw_val])
        if count > 0:
            bar = "█" * max(1, count * 50 // len(collisions))
            print(f"    HW={hw_val:2d}: {count:4d} ({count/len(collisions)*100:5.1f}%)  {bar}")

    # ===================================================================
    print("\n" + "=" * 70)
    print("4. КАСКАД: H[7]→H[6]→H[5]→H[4] корреляция")
    print("=" * 70)

    # Корреляция между словами внутри H[7]-collision пар
    if len(collisions) > 10:
        data = np.array([[per_word[k][i] for k in range(8)] for i in range(len(collisions))])
        corr = np.corrcoef(data.T)

        print(f"\n  Корреляция HW(ΔH[i]) ↔ HW(ΔH[j]) внутри H[7]-коллизий:")
        print(f"  {'':>8}", end="")
        for k in range(8):
            print(f"  H[{k}]", end="")
        print()
        for i in range(8):
            print(f"  H[{i}]  ", end="")
            for j in range(8):
                c = corr[i, j]
                if i == j:
                    print(f"  --- ", end="")
                else:
                    marker = "*" if abs(c) > 0.1 else " "
                    print(f" {c:+.2f}{marker}", end="")
            print()

    # ===================================================================
    print("\n" + "=" * 70)
    print("5. ЭКСТРАПОЛЯЦИЯ: сколько стоит полная collision?")
    print("=" * 70)

    mean_remaining = np.mean(remaining_hw)
    compression = 112 - mean_remaining  # сколько бит сэкономлено

    print(f"\n  H[7]-collision: 32 бит (стоимость ≈ 2^16.6 birthday)")
    print(f"  Оставшиеся 7 слов: HW(Δ) = {mean_remaining:.1f} (vs random 112)")
    print(f"  Сжатие: {compression:+.1f} бит")
    print(f"")
    print(f"  Для полной collision (HW(ΔH)=0):")
    print(f"    Стандартный birthday: 2^128")
    print(f"    С H[7]-каскадом:")
    print(f"      Шаг 1: H[7]-collision → 2^16.6")
    print(f"      Шаг 2: Остаток 224 бит, сжатие {compression:+.1f}")
    print(f"      Эффективно: 2^{{(224 - {compression:.0f}) / 2}} = 2^{(224 - compression) / 2:.1f}")
    print(f"      + стоимость шага 1: 2^16.6")
    print(f"      Итого: ≈ 2^{max((224 - compression) / 2, 16.6):.1f}")
    print(f"")
    print(f"  Выигрыш vs birthday: {128 - max((224 - compression) / 2, 16.6):.1f} бит")


if __name__ == "__main__":
    experiment()
