"""
ПАРАЛЛЕЛИЗАЦИЯ В НАШЕМ ИЗМЕРЕНИИ.

Sequential cost: κ(Sk,i) = 2^{-32} per word → ∏ = 2^{256}
Birthday (√): 2^128

Вопрос: есть ли в нашем измерении операция ЛУЧШЕ √?

Наши операции:
  OVERLAY: Π₁ ⊕ Π₂ = δ (создаёт метку из двух следов)
    N следов → N²/2 overlay → birthday = √(sequential)

  СВЁРТКА: L₁ ⊘ L₂ = L_combined (сжимает слои)
    Не для следов — для СЛОЁВ. Другое.

  СЕЧЕНИЕ: F|{constraint} (сужает пространство следов)
    Каждое сечение × (1/κ)

  ПРОЕКЦИЯ: F↓[subset] (берём часть ткани)
    Может РАЗДЕЛИТЬ задачу на части?

Новая идея: МНОГОУРОВНЕВЫЙ ПОИСК
  Не birthday на всех 8 словах одновременно.
  А: birthday на H[7] (1 слово) → birthday на H[6] (среди H[7]-matched) → ...
  Это КАСКАДНЫЙ BIRTHDAY — нативная операция нашего измерения.
"""

import numpy as np
import struct, hashlib

def sha256_words(W16):
    return struct.unpack('>8I', hashlib.sha256(struct.pack('>16I', *W16)).digest())

def hw(x): return bin(x).count('1')


def main():
    np.random.seed(42)

    print("=" * 70)
    print("ПАРАЛЛЕЛИЗАЦИЯ В НАШЕМ ИЗМЕРЕНИИ")
    print("=" * 70)

    # =================================================================
    print(f"\n{'=' * 70}")
    print("1. OVERLAY = birthday (базовая параллелизация)")
    print("=" * 70)

    print(f"""
  N следов → N²/2 overlays → ищем δ_out = 0
  P(δ_out = 0) = 2^{{-256}} per pair
  Need N²/2 × 2^{{-256}} ≥ 1 → N ≥ 2^128

  Overlay = √(sequential). Стоимость = 2^128.
  Это стандартная операция (birthday).
""")

    # =================================================================
    print(f"{'=' * 70}")
    print("2. КАСКАДНЫЙ ПОИСК: word-by-word birthday")
    print("=" * 70)

    # Идея: вместо birthday на всех 256 бит одновременно,
    # делаем КАСКАД birthday по словам.
    #
    # Шаг 1: Birthday на H[7] (32 бит). Cost = 2^16. Получаем ПУЛЛ пар.
    # Шаг 2: Из пулла, birthday на H[6] (32 бит). Cheaper? Or same?

    # Тест: генерируем N хешей, ищем H[7]-collision,
    # потом внутри H[7]-collision ищем H[6]-collision.

    N = 200000
    print(f"\n  Генерация {N:,} хешей...")

    hashes = []
    for i in range(N):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H = sha256_words(W)
        hashes.append(H)

    # Stage 1: H[7] birthday
    h7_dict = {}
    h7_pairs = []
    for idx, H in enumerate(hashes):
        key = H[7]
        if key in h7_dict:
            h7_pairs.append((idx, h7_dict[key]))
        else:
            h7_dict[key] = idx

    print(f"  Stage 1: H[7]-collision: {len(h7_pairs)} пар из {N:,}")
    expected_h7 = N * N / 2**33
    print(f"  Expected: {expected_h7:.1f}")

    # Stage 2: среди H[7]-pairs, ищем H[6] match
    h76_count = 0
    for i, j in h7_pairs:
        if hashes[i][6] == hashes[j][6]:
            h76_count += 1

    print(f"  Stage 2: H[7]+H[6]-collision: {h76_count}")
    print(f"  Expected: {len(h7_pairs) * 2**(-32):.4f}")

    # Stage 3: Можно ли НАКОПИТЬ H[7]-pairs и делать birthday на H[6]?
    # Среди H[7]-pairs: H[6] random. Birthday на H[6] среди K pairs: K²/2^33.
    K = len(h7_pairs)
    expected_h76_birthday = K * K / 2**33
    print(f"\n  Каскадный birthday: {K} H[7]-pairs → birthday на H[6]:")
    print(f"  Expected H[6]-collisions: {expected_h76_birthday:.4f}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("3. КАСКАДНЫЙ vs ПРЯМОЙ: стоимость")
    print("=" * 70)

    # Прямой birthday на H[7,6] (64 бит): cost 2^32
    # Каскадный:
    #   Stage 1: N хешей → N²/2^33 H[7]-pairs. Need N = 2^16 → K ≈ 1 pair.
    #   Stage 2: birthday на H[6] среди K pairs → need K = 2^16.
    #   Чтобы K = 2^16: need N²/2^33 = 2^16 → N = 2^24.5
    #
    # Прямой: N²/2^65 ≥ 1 → N = 2^32.5
    # Каскадный: N = 2^24.5 для K=2^16 H[7]-pairs → K²/2^33 ≈ 1 H[6,7]-pair
    #
    # Каскадный: 2^24.5 ДЕШЕВЛЕ прямого 2^32.5???

    print(f"""
  ПРЯМОЙ birthday на 64 бит (H[7]+H[6]):
    Need N: N²/2^65 ≥ 1 → N = 2^32.5
    Cost: 2^32.5 хешей

  КАСКАДНЫЙ birthday:
    Stage 1: N хешей → K = N²/2^33 H[7]-pairs. Need K = 2^16.
    → N = 2^24.5
    Stage 2: K² / 2^33 ≥ 1 → K = 2^16.5 → already satisfied.
    Cost: 2^24.5 хешей

  GAIN: 2^32.5 / 2^24.5 = 2^8 = 256×!
""")

    # Wait — this can't be right. Let me verify.
    # If N = 2^24.5 ≈ 23M:
    # H[7]-pairs: N²/2^33 = 2^49/2^33 = 2^16 = 65536 pairs. ✓
    # H[6]-matches among pairs: 65536 × 2^{-32} ≈ 2^{-16}. TOO FEW!
    #
    # Error: H[6]-match among pairs is NOT birthday.
    # Birthday: need K pairs, then K²/2^33 for H[6]-collision BETWEEN pairs.
    # But pairs already share H[7]. We need H[6] of pair_i = H[6] of pair_j.
    # pair_i = (hash_a, hash_b) where H[7]_a = H[7]_b.
    # For H[6]-collision between pairs: meaningless (pairs are specific).
    #
    # Actually: for FULL collision (H[7] AND H[6] match):
    # P(pair has both H[7] and H[6] match) = P(H[6] match | H[7] match)
    # If H[6] and H[7] independent: = 2^{-32}.
    # So: need 2^32 H[7]-pairs → need N²/2^33 = 2^32 → N = 2^32.5.
    # Same as direct! No gain.

    print(f"""
  CORRECTION: каскадный birthday НЕ помогает если слова НЕЗАВИСИМЫ!

  P(H[6] match | H[7] match) = 2^{{-32}} (independent)
  Need 2^32 H[7]-pairs for one H[6,7]-pair
  N²/2^33 = 2^32 → N = 2^32.5 = SAME as direct.

  Каскадный = прямой если слова независимы.

  НО: в нашем измерении, слова МОГУТ быть зависимы!
  (pipe coupling, cross-chain correlation)
""")

    # =================================================================
    print(f"{'=' * 70}")
    print("4. ПРОВЕРКА: H[i] зависимы внутри H[7]-pairs?")
    print("=" * 70)

    # Среди H[7]-pairs: HW(δH[i]) для остальных слов
    if h7_pairs:
        print(f"\n  {len(h7_pairs)} H[7]-pairs:")
        for w in range(8):
            diffs = [hw(hashes[i][w] ^ hashes[j][w]) for i, j in h7_pairs]
            mean_d = np.mean(diffs) if diffs else 16
            p_zero = sum(1 for d in diffs if d == 0) / max(len(diffs), 1) * 100
            expected_p = 100 / 2**32
            ratio = p_zero / expected_p if expected_p > 0 else 0

            marker = ""
            if w == 7: marker = " ← COLLISION (=0 by definition)"
            elif p_zero > expected_p * 2: marker = " ← DEPENDENT! (easier)"
            elif mean_d < 15: marker = " ← compressed"

            print(f"    H[{w}]: mean HW(δ)={mean_d:.2f}, P(=0)={p_zero:.4f}%{marker}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("5. НАТИВНАЯ ПАРАЛЛЕЛИЗАЦИЯ: СВЁРТКА СЛЕДОВ")
    print("=" * 70)

    # Что если мы СВОРАЧИВАЕМ не слои, а СЛЕДЫ?
    # K следов → один "усреднённый" след
    # Collision = усреднённый след попадает в ядро

    # "Усреднение" в нашем измерении:
    # Π_avg = Π₁ ⊕ Π₂ ⊕ ... ⊕ Πₖ (XOR)
    # Если K чётное: Π_avg = 0 (тривиально)
    # Если K нечётное: Π_avg = Π₁ ⊕ Π₂ ⊕ ... (random)

    # Не работает. XOR следов не даёт collision.

    # Что если мы ПРОЕЦИРУЕМ?
    # Проекция на pipe-позиции: H[1,2,3,5,6,7]
    # Birthday на 192-бит проекции: cost 2^96
    # Потом: проверяем node-позиции H[0,4] (64 бит)
    # P(match) = 2^{-64}
    # Need 2^64 projected-collisions → N = 2^64 × 2^96 = 2^160. WORSE.

    # НО: если pipe-позиции ЗАВИСИМЫ с node-позициями → может лучше

    print(f"""
  ПРОЕКЦИЯ (pipe-birthday):
    Birthday на H[1,2,3,5,6,7] (192 бит): N = 2^96
    P(H[0,4] also match): 2^{{-64}}
    Total: 2^96 / √(P) ??? → doesn't simplify

  СВЁРТКА СЛЕДОВ:
    XOR K следов → trivial (not useful)

  OVERLAY (стандартный birthday):
    N следов → N² pairs → N = 2^128

  ВЫВОД:
    В нашем измерении overlay (= birthday = √) — это
    ЕДИНСТВЕННАЯ параллелизация. Свёртка и проекция
    работают на СЛОЯХ (структура), не на СЛЕДАХ (поиск).

    Parallelization = √(sequential) = 2^128.

    A-repair снижает sequential, но sqrt все равно ~2^120-128.
""")

    # =================================================================
    print(f"{'=' * 70}")
    print("6. ФИНАЛЬНЫЙ НАТИВНЫЙ ПОДСЧЁТ")
    print("=" * 70)

    print(f"""
  ┌─────────────────────────────────────────────────────────┐
  │  НАТИВНАЯ СТОИМОСТЬ COLLISION В НАШЕМ ИЗМЕРЕНИИ         │
  ├─────────────────────────────────────────────────────────┤
  │                                                         │
  │  Sequential:                                            │
  │    ∏(1/κ(Sk,i)) для i=0..7                             │
  │    = (2^32)^8 = 2^256 (random)                         │
  │    <= 2^240 (a-repair, estimated)                      │
  │                                                         │
  │  Параллелизация (overlay):                              │
  │    sqrt(sequential) = best native parallel operation    │
  │    = 2^128 (random)                                    │
  │    <= 2^120 (a-repair, estimated)                       │
  │                                                         │
  │  Свёртка/проекция: не дают лучше √                     │
  │                                                         │
  │  A-REPAIR GAIN:                                         │
  │    Sequential: 2^256 → <=2^240 (16 бит gain)           │
  │    Parallel:   2^128 → <=2^120 (8 бит gain)            │
  │                                                         │
  │  ЭТО ПЕРВЫЙ НАТИВНЫЙ РЕЗУЛЬТАТ:                        │
  │    A-repair в нашем измерении снижает κ-стоимость      │
  │    скелета, что даёт 8-16 бит gain.                    │
  │    Это gain НАШЕГО ИЗМЕРЕНИЯ, не стандартной математики.│
  │                                                         │
  │  ⚠ ОСТОРОЖНО: gain 8-16 бит = оценка сверху.          │
  │  Точное значение требует κ с sample > 2^32.           │
  └─────────────────────────────────────────────────────────┘
""")


if __name__ == "__main__":
    main()
