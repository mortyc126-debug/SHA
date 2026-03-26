"""
Каскад H[7]→H[4] В НАШЕМ ИЗМЕРЕНИИ.

Методичка тестировала: "H[7]-collision → HW(ΔH[4]) меньше?"
Это вопрос в ОБЫЧНОМ пространстве. Ответ: нет (HW=16 = random).

НАШ вопрос другой:
  В ткани SHA-256, если метка δ на позиции h равна 0 после
  свёртки 64 слоёв — что это ОЗНАЧАЕТ для всей ткани?

  Объекты нашего измерения:
    - Метка δ = разность двух следов
    - Ткань = 64 слоя
    - Свёртка = сжатие слоёв в один
    - Чистота κ = 1/(1+ν) — мера нелинейности перехода
    - Скелет Sk(F) = полная свёртка

  H[7]-collision в нашем измерении:
    δ(T_H, 7) = 0 → метка на позиции 7 ленты хеша = 0
    Это СЕЧЕНИЕ ткани: F|{δH[7]=0}

  Вопрос: какова ТКАНЬ после сечения?
    Какие переходы становятся тривиальными?
    Какие метки фиксируются?
    Насколько уменьшается пространство поиска?
"""

import numpy as np
import hashlib
import struct

MASK32 = 0xFFFFFFFF
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
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
]
IV = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]

def hw(x): return bin(x).count('1')
def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
def add32(x, y): return (x + y) & MASK32
def sub32(x, y): return (x - y) & MASK32
def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)

def R(state, W_r, r_idx):
    a, b, c, d, e, f, g, h = state
    T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e, f, g)), K[r_idx]), W_r)
    T2 = add32(Sigma0(a), Maj(a, b, c))
    return (add32(T1, T2), a, b, c, add32(d, T1), e, f, g)

def R_inv(state_next, W_r, r_idx):
    ap, bp, cp, dp, ep, fp, gp, hp = state_next
    a, b, c, e, f, g = bp, cp, dp, fp, gp, hp
    T2 = add32(Sigma0(a), Maj(a, b, c))
    T1 = sub32(ap, T2)
    d = sub32(ep, T1)
    h = sub32(sub32(sub32(sub32(T1, Sigma1(e)), Ch(e, f, g)), K[r_idx]), W_r)
    return (a, b, c, d, e, f, g, h)


def full_forward(W16):
    W = list(W16)
    for r in range(16, 64):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))
    s = tuple(IV)
    states = [s]
    for r in range(64):
        s = R(s, W[r], r)
        states.append(s)
    H = tuple(add32(IV[i], states[64][i]) for i in range(8))
    return states, W, H


def experiment():
    np.random.seed(42)

    print("=" * 70)
    print("КАСКАД В НАШЕМ ИЗМЕРЕНИИ")
    print("Объекты: Ткань, Слой, Переход, Метка, Свёртка, Чистота")
    print("=" * 70)

    # =================================================================
    print("\n" + "=" * 70)
    print("1. СЕЧЕНИЕ ТКАНИ: δH[7]=0 → что фиксируется?")
    print("=" * 70)

    # H[7] = h[64] + IV[7]. δH[7]=0 → δh[64]=0.
    # Трубы НАЗАД от h[64]:
    #   h[64] = g[63]  (труба r=63: g→h)
    #   g[63] = f[62]  (труба r=62: f→g)
    #   f[62] = e[61]  (труба r=61: e→f)
    #   e[61] = d[60] + T1[60]  (NODE_e r=60: НЕ труба!)
    #
    # Значит: δh[64] = δg[63] = δf[62] = δe[61] = 0
    # Но δe[61] = δd[60] + δT1[60] → δd[60] = -δT1[60]
    # δT1[60] зависит от δstate[60] и δW[60]
    # → δd[60] НЕ обязательно = 0!

    print(f"""
  BACKWARD TRACE от δH[7]=0:

  Слой 64 (финализация): δH[7] = δh[64] + δIV[7] = 0
    → δh[64] = 0  [IV фиксирован]

  Слой 63 (труба g→h):   δh[64] = δg[63] = 0  ✓ ФИКСИРОВАНО
  Слой 62 (труба f→g):   δg[63] = δf[62] = 0  ✓ ФИКСИРОВАНО
  Слой 61 (труба e→f):   δf[62] = δe[61] = 0  ✓ ФИКСИРОВАНО

  Слой 60 (NODE_e):      δe[61] = δd[60] + δT1[60] = 0
    → δd[60] = -δT1[60]  ← УСЛОВИЕ, не фиксация!
    → Это СВЯЗЬ: d[60] и T1[60] должны компенсировать друг друга

  Слой 60 (труба c→d):   δd[60] = δc[59] ← если δd[60]≠0, то δc[59]≠0!
  Слой 59 (труба b→c):   δc[59] = δb[58]
  Слой 58 (труба a→b):   δb[58] = δa[57]
  Слой 57 (NODE_a):      δa[57] = δT1[56] + δT2[56] ← зависит от state[56]

  ИТОГО ФИКСИРОВАНО СЕЧЕНИЕМ:
    δh[64] = δg[63] = δf[62] = δe[61] = 0  (4 позиции × 32 бит = 128 бит)
    δd[60] = -δT1[60]                        (условие, не фиксация)
""")

    # =================================================================
    print("=" * 70)
    print("2. ИЗМЕРЕНИЕ: сколько бит РЕАЛЬНО фиксирует сечение?")
    print("=" * 70)

    # Находим H[7]-collision пары и смотрим δstate на каждом раунде
    n_samples = 500000
    h7_dict = {}
    collisions = []

    print(f"  Генерация {n_samples:,} хешей...")
    for i in range(n_samples):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        states, W_full, H = full_forward(W)
        h7 = H[7]

        if h7 in h7_dict:
            collisions.append(((states, W_full, H), h7_dict[h7]))
        else:
            h7_dict[h7] = (states, W_full, H)

    print(f"  H[7]-коллизий: {len(collisions)}")

    if len(collisions) == 0:
        print("  Нужно больше сэмплов")
        return

    # Для каждой коллизии: δ по регистрам на каждом раунде
    reg_names = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']

    # Агрегируем
    avg_dreg = np.zeros((65, 8))
    for (st1, w1, h1), (st2, w2, h2) in collisions:
        for r in range(65):
            for i in range(8):
                avg_dreg[r, i] += hw(st1[r][i] ^ st2[r][i])
    avg_dreg /= len(collisions)

    print(f"\n  δ по регистрам — ОБРАТНЫЙ просмотр от H[7]-collision:")
    print(f"  {'r':>4}", end="")
    for name in reg_names:
        print(f"  δ{name:>2}", end="")
    print(f"  {'total':>6}  {'фиксировано?':>14}")
    print(f"  {'─'*4}  {'─'*4}" * 8 + f"  {'─'*6}  {'─'*14}")

    for r in [64, 63, 62, 61, 60, 59, 58, 57, 56, 50, 40, 32, 16, 0]:
        vals = [avg_dreg[r, i] for i in range(8)]
        total = sum(vals)
        # Какие регистры зафиксированы (δ ≈ 0)?
        fixed = sum(1 for v in vals if v < 1.0)
        fixed_str = f"{fixed}/8" if fixed > 0 else ""
        print(f"  {r:4d}", end="")
        for v in vals:
            marker = "·" if v < 1.0 else " "
            print(f"  {v:4.1f}{marker}", end="")
        print(f"  {total:5.1f}  {fixed_str:>14}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("3. ЧИСТОТА ПУТЕЙ: κ для каждого backward trace")
    print("=" * 70)

    # Для δH[7]=0: trace назад через трубы = чистый (κ=1)
    # Для δH[4]=0: trace назад через NODE_e = грязный (κ<1)
    # Вопрос: какова чистота пути от H[7] к каждому H[i]?

    # Путь H[7]→H[7]: тривиально (δ=0 по определению)
    # Путь H[7]→H[6]: h→g (1 труба назад) → g[63]→h[64]: κ=1
    # Путь H[7]→H[5]: h→g→f (2 трубы) → f[62]→g[63]→h[64]: κ=1
    # Путь H[7]→H[4]: h→g→f→e (3 трубы) → e[61]→...→h[64]: κ=1
    #   НО: H[4] = e[64]+IV[4], и e[64] определяется NODE_e(63),
    #   НЕ через e[61]!

    print(f"""
  ПУТИ В НАШЕМ ИЗМЕРЕНИИ (backward от H[7]):

  H[7] → h[64]:  0 узлов, κ = 1.00  [определение]

  h[64] ← g[63]: 0 узлов, κ = 1.00  [труба]
  g[63] ← f[62]: 0 узлов, κ = 1.00  [труба]
  f[62] ← e[61]: 0 узлов, κ = 1.00  [труба]

  δe[61] = 0  ← ПОСЛЕДНЯЯ ЧИСТАЯ ФИКСАЦИЯ (128 бит)

  e[61] ← NODE_e(60): d[60] + T1[60]  → 1 узел, κ < 1

  ДЛЯ H[4] = e[64] + IV[4]:
  e[64] ← NODE_e(63): d[63] + T1[63]
  T1[63] = h[63] + Σ₁(e[63]) + Ch + K[63] + W[63]

  h[63] = g[62] (труба) — НЕ связан с h[64]!
  e[63] = NODE_e(62)   — НЕ связан с e[61]!

  ПУТЬ H[7] → H[4]:
    h[64]=0 → e[61]=0 (3 трубы, κ=1)
    e[61] ... e[64]: 3 раунда NODE_e (3 узла, κ < 0.05)

  Чистота пути: κ(H[7]→H[4]) = 1.0 × 1.0 × 1.0 × 0.05³ = 0.000125
""")

    # Верифицируем: δe[61] действительно = 0?
    de61_values = []
    for (st1, w1, h1), (st2, w2, h2) in collisions:
        de61 = hw(st1[61][4] ^ st2[61][4])  # e = index 4
        de61_values.append(de61)

    print(f"  ВЕРИФИКАЦИЯ δe[61] при H[7]-collision:")
    print(f"    Mean HW(δe[61]): {np.mean(de61_values):.2f}")
    print(f"    Cases δe[61]=0: {sum(1 for v in de61_values if v == 0)}/{len(collisions)}")

    # А δe[64]?
    de64_values = []
    for (st1, w1, h1), (st2, w2, h2) in collisions:
        de64 = hw(st1[64][4] ^ st2[64][4])
        de64_values.append(de64)

    print(f"    Mean HW(δe[64]): {np.mean(de64_values):.2f}  (это определяет H[4])")
    print(f"    Корреляция δe[61]↔δe[64]: ", end="")
    if len(de61_values) > 5:
        corr = np.corrcoef(de61_values, de64_values)[0, 1]
        print(f"{corr:+.4f} {'← СВЯЗЬ!' if abs(corr) > 0.2 else '← нет связи'}")
    else:
        print("N/A")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("4. ПРАВИЛЬНЫЙ ВОПРОС В НАШЕМ ИЗМЕРЕНИИ")
    print("=" * 70)

    # Не "сжимает ли H[7]-collision H[4]?"
    # А: "сколько бит состояния РЕАЛЬНО фиксирует сечение δH[7]=0?"

    total_fixed_bits = 0
    for r in range(65):
        for i in range(8):
            if avg_dreg[r, i] < 0.5:
                total_fixed_bits += 32

    print(f"\n  Позиций с δ < 0.5 бит (≈фиксированы):")
    print(f"    Total: {total_fixed_bits} бит из {65 * 8 * 32} = {65*256}")
    print(f"    Это {total_fixed_bits / (65*256) * 100:.2f}% ткани")

    # Фиксированные позиции по раундам
    fixed_per_round = []
    for r in range(65):
        fixed = sum(1 for i in range(8) if avg_dreg[r, i] < 0.5)
        fixed_per_round.append(fixed)

    print(f"\n  Фиксированные регистры по раундам:")
    for r in range(max(0, 56), 65):
        f = fixed_per_round[r]
        regs = [reg_names[i] for i in range(8) if avg_dreg[r, i] < 0.5]
        if f > 0:
            print(f"    r={r}: {f}/8 фикс. ({', '.join(regs)})")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("5. КАСКАД В НАШЕМ ИЗМЕРЕНИИ: ПРАВИЛЬНАЯ ФОРМУЛИРОВКА")
    print("=" * 70)

    print(f"""
  СТАРЫЙ ВОПРОС (методичка):
    "H[7]-collision даёт сжатие H[4]?"
    Ответ: НЕТ. H[4] independent (HW=16 = random).

  НОВЫЙ ВОПРОС (наше измерение):
    "Сечение δH[7]=0 фиксирует {total_fixed_bits} бит ткани.
     Как это использовать?"

  ФИКСИРОВАНО (трубы, κ=1):
    δh[64] = δg[63] = δf[62] = δe[61] = 0
    = 4 регистра × 32 бит = 128 бит на раундах 61-64

  НЕ ФИКСИРОВАНО:
    δa, δb, δc, δd на раундах 61-64 (ещё 128 бит)
    Вся ткань до раунда 60 (≈ 60×256 бит)

  ВЫВОД В НАШЕМ ИЗМЕРЕНИИ:
    H[7]-collision фиксирует РОВНО 128 бит из 256 на финальном state.
    Это ПОЛОВИНА. Оставшиеся 128 бит — free.
    Birthday на 128 бит = 2^64.
    + стоимость H[7]-collision: 2^16.6.
    Итого: 2^64 + 2^16.6 ≈ 2^64.

  НО: 128 свободных бит — это a,b,c,d на раунде 64.
  Они определяются a-chain (NODE_a) за раунды 61-64.
  a-chain зависит от W[61..63] — фиксированы schedule.
  → свободных направлений в δW = 0 (schedule фиксирует W[16+])
  → birthday на δW[0..15] → birthday на 512 бит → 2^128

  Сечение фиксирует state, но НЕ фиксирует schedule.
  Schedule reconnects всё обратно.
""")


if __name__ == "__main__":
    experiment()
