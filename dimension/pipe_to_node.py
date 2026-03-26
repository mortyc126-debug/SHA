"""
PIPE→NODE: при pipe-collision, что происходит с node-позициями?

Pipe-collision: H₁[1,2,3,5,6,7] = H₂[1,2,3,5,6,7] (6 слов match)
Стоимость: C^3 = 2^96

Вопрос: при pipe-collision, P(H₁[0]=H₂[0]) и P(H₁[4]=H₂[4]) = ?
  Если > 1/C: node-позиции ЗАВИСИМЫ от pipe → дешевле
  Если = 1/C: независимы → 2^96 + 2^64 = 2^96
  Если < 1/C: anti-зависимы → дороже

Из трубной теории (Этап 1.3):
  H[0] = a[64]+IV[0] = NODE_a(63)+IV[0]
  H[4] = e[64]+IV[4] = NODE_e(63)+IV[4]
  H[1] = b[64]+IV[1] = a[63]+IV[1]  (pipe a→b)
  H[7] = h[64]+IV[7] = e[61]+IV[7]  (pipe e→f→g→h)

  Pipe-collision: a[61..63] и e[61..63] совпадают.
  Node: NODE_a(63) зависит от a[63],b[63],c[63] — которые совпадают!
  Значит T2(63) = Σ₀(a[63])+Maj(a[63],b[63],c[63]) ОДИНАКОВ!
  → δH[0] = δ(NODE_a(63)) = δT1(63) + 0 = δT1(63)
  → δT1(63) зависит от h[63],e[63],f[63],g[63],W[63]
  → h[63]=e[61] (pipe), e[63],f[63],g[63] тоже pipes
  → ВСЕ входы T1(63) = PIPE → при pipe-collision ВСЕ РАВНЫ!
  → δT1(63) = δW[63] ТОЛЬКО
  → δH[0] = δW[63]

  Если δW[63] = 0: δH[0] = 0! БЕСПЛАТНО!
  Но δW[63] = schedule diff ≠ 0 (в общем случае).

Аналогично: δH[4] = δ(NODE_e(63)) = δd[63] + δT1(63) = δa[60] + δW[63]+carry

Проверяем экспериментально.
"""

import numpy as np
import struct, hashlib

MASK32 = 0xFFFFFFFF

def sha256_words(W16):
    return struct.unpack('>8I', hashlib.sha256(struct.pack('>16I', *W16)).digest())

def hw(x): return bin(x).count('1')


def main():
    np.random.seed(42)

    print("=" * 70)
    print("PIPE→NODE: зависимость node от pipe при pipe-collision")
    print("=" * 70)

    # =================================================================
    print(f"\n{'=' * 70}")
    print("1. ТЕОРИЯ: при pipe-collision, δH[0] = δW[63]")
    print("=" * 70)

    print(f"""
  Pipe-collision: H₁[1,2,3,5,6,7] = H₂[1,2,3,5,6,7]
  → a₁[61..63] = a₂[61..63], e₁[61..63] = e₂[61..63]
  → T2₁(63) = T2₂(63) (same a,b,c → same Σ₀, Maj)
  → δH[0] = δ(T1(63) + T2(63)) = δT1(63) + 0

  T1(63) = h[63] + Σ₁(e[63]) + Ch(e[63],f[63],g[63]) + K[63] + W[63]

  При pipe-collision:
    h[63] = g[62] = f[61] = e[60] → δh[63] = δe[60]
    e[63] → δe[63] = 0 (pipe-collision)
    f[63] = e[62] → δf[63] = 0
    g[63] = e[61] → δg[63] = 0

  → δΣ₁ = 0, δCh = 0
  → δT1(63) = δh[63] + 0 + 0 + 0 + δW[63]
            = δe[60] + δW[63]

  При pipe-collision: δe[61]=0, но δe[60] НЕ ГАРАНТИРОВАН!

  Также: δd[63] = δc[62] = δb[61] = δa[60]
  → δH[4] = δd[63] + δT1(63) = δa[60] + δe[60] + δW[63]

  Из уравнений:
    δH[0] = δe[60] + δW[63]
    δH[4] = δa[60] + δe[60] + δW[63]
    δH[4] - δH[0] = δa[60]

  ЭТО ТЕ ЖЕ УРАВНЕНИЯ из pipe_cascade.py!
  Но теперь в контексте pipe-collision, не full collision.
""")

    # =================================================================
    print(f"{'=' * 70}")
    print("2. ЭКСПЕРИМЕНТ: H[7]-collision → что с H[0] и H[4]?")
    print("=" * 70)

    # Начнём с H[7]-collision (проще найти: 32-бит birthday)
    # H[7]-collision → e[61] совпадает (трубы)
    # → часть pipe-collision. Проверяем H[0,4].

    N = 2000000
    h7_dict = {}
    h7_pairs = []

    for _ in range(N):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H = sha256_words(W)
        k = H[7]
        if k in h7_dict:
            h7_pairs.append((H, h7_dict[k]))
            if len(h7_pairs) >= 500:
                break
        h7_dict[k] = H

    print(f"\n  H[7]-collisions найдено: {len(h7_pairs)}")

    if h7_pairs:
        dh0 = [hw(a[0] ^ b[0]) for a, b in h7_pairs]
        dh4 = [hw(a[4] ^ b[4]) for a, b in h7_pairs]
        dh04_xor = [hw((a[0]^b[0]) ^ (a[4]^b[4])) for a, b in h7_pairs]

        # Все слова
        print(f"\n  При H[7]-collision, HW(δH[i]):")
        for w in range(8):
            vals = [hw(a[w] ^ b[w]) for a, b in h7_pairs]
            expected = 0 if w == 7 else 16
            delta = expected - np.mean(vals)
            marker = ""
            if abs(delta) > 0.5: marker = f" ← {delta:+.1f} COMPRESSED" if delta > 0 else f" ← {delta:+.1f} EXPANDED"
            print(f"    H[{w}]: mean={np.mean(vals):.2f} (expected {expected}){marker}")

        print(f"\n  δH[0]: mean={np.mean(dh0):.2f}")
        print(f"  δH[4]: mean={np.mean(dh4):.2f}")
        print(f"  δH[4]⊕δH[0]: mean HW={np.mean(dh04_xor):.2f} (= δa[60] theoretically)")

        # P(δH[0]=0) | H[7]-collision
        p_h0_zero = sum(1 for d in dh0 if d == 0) / len(dh0)
        p_h4_zero = sum(1 for d in dh4 if d == 0) / len(dh4)
        p_both_zero = sum(1 for i in range(len(dh0)) if dh0[i] == 0 and dh4[i] == 0) / len(dh0)

        print(f"\n  Conditional P при H[7]-collision:")
        print(f"    P(δH[0]=0): {p_h0_zero*100:.3f}%")
        print(f"    P(δH[4]=0): {p_h4_zero*100:.3f}%")
        print(f"    P(both=0):  {p_both_zero*100:.3f}%")
        print(f"    Unconditional P(=0): {100/2**32:.8f}%")

        if p_h0_zero > 0:
            lift = p_h0_zero / (1/2**32)
            print(f"    LIFT H[0]: {lift:.1f}×")
        if p_both_zero > 0:
            lift_both = p_both_zero / (1/2**64)
            print(f"    LIFT both: {lift_both:.1f}×")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("3. MULTI-WORD: H[7]+H[6]-collision → H[0,4]?")
    print("=" * 70)

    # Более полная pipe-collision: H[7] AND H[6] match
    h76_dict = {}
    h76_pairs = []

    for _ in range(N):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H = sha256_words(W)
        k = (H[7], H[6])
        if k in h76_dict:
            h76_pairs.append((H, h76_dict[k]))
        h76_dict[k] = H

    print(f"\n  H[7,6]-collisions: {len(h76_pairs)}")
    if h76_pairs:
        for w in range(8):
            vals = [hw(a[w] ^ b[w]) for a, b in h76_pairs]
            expected = 0 if w in [6,7] else 16
            delta = expected - np.mean(vals)
            marker = f" ← {delta:+.1f}" if abs(delta) > 1 else ""
            print(f"    H[{w}]: mean={np.mean(vals):.2f}{marker}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("4. НАТИВНАЯ СТОИМОСТЬ: pipe-collision + node dependency")
    print("=" * 70)

    # Из теории и эксперимента:
    # При pipe-collision (H[1,2,3,5,6,7] match):
    #   δH[0] = δe[60] + δW[63]
    #   δH[4] = δa[60] + δe[60] + δW[63]

    # δe[60] и δa[60] зависят от ВСЕХ предыдущих раундов.
    # При pipe-collision (random пара): δe[60],δa[60] ≈ random.
    # P(δH[0]=0) = P(δe[60] = -δW[63]) ≈ 1/C (random)
    # P(δH[4]=0 | δH[0]=0) = P(δa[60]=0) ≈ 1/C

    # НО: при a-repair pipe-collision:
    # a-repair фиксирует δa = 0 на r=4..15 → δa[60] МОЖЕТ быть подавлен!

    if h7_pairs:
        mean_h0 = np.mean(dh0)
        mean_h4 = np.mean(dh4)

        print(f"""
  При H[7]-collision (e[61] совпадает):
    δH[0] = δe[60] + δW[63]: mean HW = {mean_h0:.1f} (random = 16)
    δH[4] = δa[60] + δe[60] + δW[63]: mean HW = {mean_h4:.1f}

  COMPRESSION: {'ДА' if mean_h0 < 15.5 else 'НЕТ'} (mean = {mean_h0:.2f} vs 16)

  Стоимость в нашем мире:
    Pipe-collision (6 позиций): C^3 = 2^96
    Node-позиция H[0]: C^1 = 2^32 (if independent)
    Node-позиция H[4]: C^1 = 2^32 (if independent)
    Full collision: C^3 × C^1 × C^1 = C^5 = 2^160 ??? ХУЖЕ!

  НО: если через overlay на pipe-collision СНАЧАЛА:
    Birthday на pipe (192 бит): N = C^3 = 2^96
    Каждая pipe-collision: P(full) = 1/C^2 (node match)
    Need C^2 pipe-collisions: N = C^3 × C^2 = C^5 = 2^160. ХУЖЕ.

  ОШИБКА: pipe-collision + sequential node ≠ full collision birthday.
  Правильно: birthday на ПОЛНЫЕ 256 бит = C^4 = 2^128.
  Pipe-decomposition НЕ помогает (6+2 = 8 = original).
""")


if __name__ == "__main__":
    main()
