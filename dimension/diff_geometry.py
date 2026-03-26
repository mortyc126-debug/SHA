"""
ДИФФЕРЕНЦИАЛЬНАЯ ГЕОМЕТРИЯ ТКАНИ.

Ткань = отображение F: input → output.
В нашем измерении: F "изгибает" пространство.

Кривизна = мера того, насколько F отличается от линейного.
Метрика = расстояние между следами ПОСЛЕ прохождения ткани.
Геодезические = кратчайшие пути в пространстве следов.

Определения:
  1. МЕТРИЧЕСКИЙ ТЕНЗОР: g(δ₁, δ₂) = <F(δ₁), F(δ₂)>
     Как ткань деформирует скалярное произведение меток.

  2. КРИВИЗНА: K = |F(δ₁⊕δ₂) - F(δ₁) - F(δ₂)| / |δ₁||δ₂|
     Мера нелинейности (= carry effect, normalized).

  3. ЛОКАЛЬНАЯ КРИВИЗНА: K(W, δ) — кривизна в точке W в направлении δ.
     Разная для разных точек и направлений?

  4. ГЕОДЕЗИЧЕСКАЯ: путь δ(t) в пространстве меток с минимальной
     "длиной" (= минимальной carry error).
"""

import numpy as np
import struct, hashlib

MASK32 = 0xFFFFFFFF

def sha256_words(W16):
    return struct.unpack('>8I', hashlib.sha256(struct.pack('>16I', *W16)).digest())

def hw(x): return bin(x).count('1')

def hash_xor_hw(H1, H2):
    return sum(hw(H1[i] ^ H2[i]) for i in range(8))


def main():
    np.random.seed(42)

    print("=" * 70)
    print("ДИФФЕРЕНЦИАЛЬНАЯ ГЕОМЕТРИЯ ТКАНИ")
    print("=" * 70)

    # =================================================================
    print(f"\n{'=' * 70}")
    print("1. МЕТРИЧЕСКИЙ ТЕНЗОР: как ткань деформирует расстояния")
    print("=" * 70)

    # d_input(W₁, W₂) = HW(W₁ ⊕ W₂) (расстояние на входе)
    # d_output(W₁, W₂) = HW(H₁ ⊕ H₂) (расстояние на выходе)
    # Метрика: d_output / d_input = "растяжение"

    stretches = []
    for _ in range(5000):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        # Малый δ: 1 бит
        bit_pos = np.random.randint(0, 512)
        W2 = list(W1)
        W2[bit_pos // 32] ^= (1 << (bit_pos % 32))

        H1 = sha256_words(W1)
        H2 = sha256_words(W2)

        d_in = 1  # 1 бит
        d_out = hash_xor_hw(H1, H2)
        stretches.append(d_out / d_in)

    print(f"\n  1-бит δW → HW(δH):")
    print(f"    Mean stretch: {np.mean(stretches):.1f}×")
    print(f"    Std: {np.std(stretches):.1f}")
    print(f"    Min: {min(stretches):.0f}×, Max: {max(stretches):.0f}×")

    # Для БОЛЬШИХ δ:
    for target_din in [1, 10, 50, 128, 256]:
        stretches_d = []
        for _ in range(1000):
            W1 = [np.random.randint(0, 2**32) for _ in range(16)]
            W2 = list(W1)
            bits = np.random.choice(512, target_din, replace=False)
            for b in bits:
                W2[b // 32] ^= (1 << (b % 32))

            d_in = target_din
            H1 = sha256_words(W1)
            H2 = sha256_words(W2)
            d_out = hash_xor_hw(H1, H2)
            stretches_d.append(d_out)

        print(f"  d_in={target_din:3d}: d_out mean={np.mean(stretches_d):.1f}, stretch={np.mean(stretches_d)/target_din:.2f}×")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("2. КРИВИЗНА: K(W, δ₁, δ₂)")
    print("=" * 70)

    # K = HW(F(δ₁⊕δ₂) ⊕ F(δ₁) ⊕ F(δ₂)) / (HW(δ₁) × HW(δ₂))
    # F(δ) = H(W⊕δ) ⊕ H(W)

    curvatures = []
    for _ in range(2000):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H_base = sha256_words(W)

        b1 = np.random.randint(0, 512)
        b2 = np.random.randint(0, 512)
        if b1 == b2:
            continue

        W1 = list(W); W1[b1//32] ^= (1 << (b1%32))
        W2 = list(W); W2[b2//32] ^= (1 << (b2%32))
        W12 = list(W); W12[b1//32] ^= (1 << (b1%32)); W12[b2//32] ^= (1 << (b2%32))

        H1 = sha256_words(W1)
        H2 = sha256_words(W2)
        H12 = sha256_words(W12)

        F1 = tuple(H_base[i] ^ H1[i] for i in range(8))
        F2 = tuple(H_base[i] ^ H2[i] for i in range(8))
        F12 = tuple(H_base[i] ^ H12[i] for i in range(8))
        F_sum = tuple(F1[i] ^ F2[i] for i in range(8))

        nonlin = sum(hw(F12[i] ^ F_sum[i]) for i in range(8))
        curvatures.append(nonlin)

    print(f"\n  Кривизна K = |F(δ₁⊕δ₂) ⊕ F(δ₁) ⊕ F(δ₂)| (1-бит δ):")
    print(f"    Mean: {np.mean(curvatures):.1f}")
    print(f"    Std: {np.std(curvatures):.1f}")
    print(f"    Min: {min(curvatures)}, Max: {max(curvatures)}")
    print(f"    Flat (K=0): {sum(1 for k in curvatures if k == 0)} / {len(curvatures)}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("3. ЛОКАЛЬНАЯ КРИВИЗНА: зависит ли K от W?")
    print("=" * 70)

    # Для РАЗНЫХ W: кривизна K одинакова?
    K_by_W = []
    for _ in range(100):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H_base = sha256_words(W)

        Ks = []
        for _ in range(50):
            b1, b2 = np.random.choice(512, 2, replace=False)
            W1 = list(W); W1[b1//32] ^= (1 << (b1%32))
            W2 = list(W); W2[b2//32] ^= (1 << (b2%32))
            W12 = list(W); W12[b1//32] ^= (1 << (b1%32)); W12[b2//32] ^= (1 << (b2%32))

            H1 = sha256_words(W1)
            H2 = sha256_words(W2)
            H12 = sha256_words(W12)

            F1 = tuple(H_base[i] ^ H1[i] for i in range(8))
            F2 = tuple(H_base[i] ^ H2[i] for i in range(8))
            F12 = tuple(H_base[i] ^ H12[i] for i in range(8))
            nonlin = sum(hw(F12[i] ^ (F1[i]^F2[i])) for i in range(8))
            Ks.append(nonlin)

        K_by_W.append(np.mean(Ks))

    print(f"\n  K по 100 random W:")
    print(f"    Mean: {np.mean(K_by_W):.1f}")
    print(f"    Std: {np.std(K_by_W):.1f}")
    print(f"    Min: {min(K_by_W):.1f}, Max: {max(K_by_W):.1f}")
    print(f"    → {'K ПОСТОЯННА (ткань однородно изогнута)' if np.std(K_by_W) < 5 else 'K ВАРЬИРУЕТСЯ (ткань неоднородна)'}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("4. КРИВИЗНА ПО НАПРАВЛЕНИЯМ: same word vs cross word")
    print("=" * 70)

    K_same = []  # δ₁, δ₂ в одном слове
    K_cross = []  # δ₁, δ₂ в разных словах
    K_close = []  # δ₁, δ₂ в соседних словах
    K_far = []    # δ₁, δ₂ в далёких словах

    for _ in range(3000):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H_base = sha256_words(W)

        w1 = np.random.randint(0, 16)
        w2 = np.random.randint(0, 16)
        b1 = np.random.randint(0, 32)
        b2 = np.random.randint(0, 32)

        if w1 == w2 and b1 == b2:
            continue

        pos1 = w1 * 32 + b1
        pos2 = w2 * 32 + b2

        W1 = list(W); W1[w1] ^= (1 << b1)
        W2 = list(W); W2[w2] ^= (1 << b2)
        W12 = list(W); W12[w1] ^= (1 << b1); W12[w2] ^= (1 << b2)

        H1 = sha256_words(W1)
        H2 = sha256_words(W2)
        H12 = sha256_words(W12)

        F1 = tuple(H_base[i] ^ H1[i] for i in range(8))
        F2 = tuple(H_base[i] ^ H2[i] for i in range(8))
        F12 = tuple(H_base[i] ^ H12[i] for i in range(8))
        nonlin = sum(hw(F12[i] ^ (F1[i]^F2[i])) for i in range(8))

        if w1 == w2:
            K_same.append(nonlin)
        else:
            K_cross.append(nonlin)
            if abs(w1 - w2) <= 1:
                K_close.append(nonlin)
            elif abs(w1 - w2) >= 8:
                K_far.append(nonlin)

    print(f"\n  Кривизна по типу пары:")
    print(f"    Same word:    K = {np.mean(K_same):.1f}")
    print(f"    Cross word:   K = {np.mean(K_cross):.1f}")
    if K_close:
        print(f"    Close words:  K = {np.mean(K_close):.1f}")
    if K_far:
        print(f"    Far words:    K = {np.mean(K_far):.1f}")

    diff = np.mean(K_same) - np.mean(K_cross)
    print(f"\n    Same - Cross = {diff:+.1f}")
    print(f"    → {'АНИЗОТРОПНАЯ кривизна!' if abs(diff) > 3 else 'ИЗОТРОПНАЯ (нет предпочтительного направления)'}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("5. МЕТРИЧЕСКИЙ ТЕНЗОР: РАСТЯЖЕНИЕ по позициям")
    print("=" * 70)

    # Для каждой входной позиции: среднее растяжение
    # g_ii = E[d_out | δ in position i]

    print(f"\n  Диагональ метрического тензора g_ii:")
    g_diag = []
    for word in range(16):
        stretches_w = []
        for _ in range(200):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            bit = np.random.randint(0, 32)
            W2 = list(W); W2[word] ^= (1 << bit)
            H1 = sha256_words(W)
            H2 = sha256_words(W2)
            stretches_w.append(hash_xor_hw(H1, H2))
        g_ii = np.mean(stretches_w)
        g_diag.append(g_ii)
        print(f"    g[{word:2d},{word:2d}] = {g_ii:.1f}")

    print(f"\n  Spread: {max(g_diag) - min(g_diag):.2f}")
    print(f"  → {'ПЛОСКАЯ метрика (все позиции равны)' if max(g_diag) - min(g_diag) < 3 else 'ИСКРИВЛЁННАЯ метрика'}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("6. ГЕОМЕТРИЯ НАШЕЙ ТКАНИ: РЕЗЮМЕ")
    print("=" * 70)

    print(f"""
  МЕТРИКА:
    Stretch = {np.mean(stretches):.0f}× (1-бит δ → 128 бит δH)
    Для всех d_in: d_out ≈ 128 (constant output distance)
    Метрический тензор: ПЛОСКИЙ (g_ii ≈ 128 для всех i)
    → Ткань = "сфера" (все точки эквивалентны, все направления равны)

  КРИВИЗНА:
    K = {np.mean(curvatures):.0f} (carry nonlinearity)
    K постоянна по W: std = {np.std(K_by_W):.1f} (однородная)
    K изотропна по направлениям: same ≈ cross ({abs(diff):.1f} разница)
    → Ткань = "сфера с постоянной кривизной"

  СЛЕДСТВИЯ:
    Нет "плоских" направлений (K > 0 везде).
    Нет "мягких" точек (K одинакова для всех W).
    Нет привилегированных путей (метрика плоская).
    → SHA-256 = ОДНОРОДНОЕ ИЗОТРОПНОЕ пространство.
    → Как сфера: нет кратчайшего пути, все пути эквивалентны.
    → Collision = случайная точка на сфере. Нет геодезической.
""")


if __name__ == "__main__":
    main()
