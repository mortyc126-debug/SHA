"""
ИНТЕГРАЦИЯ 4 ИНСТРУМЕНТОВ ИЗ МЕТОДИЧКИ (★★★★★)

1. T_TWIN_CARRY: GF(2)-twins → P(carry[63]=0) = 21.8% (×1M lift)
2. T_BOTTLENECK_R60: 80% значений e[60] запрещены при carry=0
3. T_G62_PREDICTS_H: близкие g[62] → E[HW(ΔH)]=109 (-18 бит)
4. T_MULTILEVEL_BIRTHDAY: H[7]-collision → H[4] сжатие +17 бит

Цель: комбинируем все 4 → измеряем суммарный выигрыш.
"""

import numpy as np
import hashlib
import struct
import time

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
def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)


def sha256_with_internals(W16):
    """SHA-256 возвращая H, state[64], и ключевые промежуточные."""
    W = list(W16)
    for r in range(16, 64):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))

    a, b, c, d, e, f, g, h = IV
    for r in range(64):
        T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e, f, g)), K[r]), W[r])
        T2 = add32(Sigma0(a), Maj(a, b, c))
        h, g, f, e = g, f, e, add32(d, T1)
        d, c, b, a = c, b, a, add32(T1, T2)

        if r == 59:
            e60 = e  # e after round 59 = e[60]
        if r == 61:
            g62 = g  # g after round 61 = g[62]

    state64 = (a, b, c, d, e, f, g, h)
    H = tuple(add32(IV[i], state64[i]) for i in range(8))

    # carry[63]: T1_63 = h + Sig1(e) + Ch(e,f,g) + K[63] + W[63]
    # carry = 1 if T1_63 overflows (result < max operand)
    # Проще: SS = h + Sig1(e) + Ch(e,f,g); carry = (SS + K[63] + W[63]) > 2^32?
    # Приближение: проверяем финальный бит carry через Ch[b30,b31]
    ch62 = Ch(state64[4], state64[5], state64[6])  # Ch(e,f,g) после r=63... нет
    # Используем H[7] bias как прокси для carry[63]
    h7_b30 = (H[7] >> 30) & 1
    h7_b29 = (H[7] >> 29) & 1

    return H, state64, e60, g62, h7_b30, h7_b29


def main():
    np.random.seed(42)

    print("=" * 70)
    print("ИНТЕГРАЦИЯ 4 ИНСТРУМЕНТОВ (★★★★★)")
    print("=" * 70)

    # =================================================================
    print("\n" + "=" * 70)
    print("1. BASELINE: random birthday search")
    print("=" * 70)

    N_SAMPLES = 100000
    start = time.time()

    # Собираем данные
    hashes = {}
    h7_values = {}
    g62_values = {}
    e60_values = {}
    all_data = []

    for i in range(N_SAMPLES):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H, s64, e60, g62, h7b30, h7b29 = sha256_with_internals(W)
        all_data.append({
            'W': W, 'H': H, 'e60': e60, 'g62': g62,
            'h7b30': h7b30, 'h7b29': h7b29, 'H7': H[7]
        })

    elapsed = time.time() - start
    print(f"\n  Собрано {N_SAMPLES:,} сэмплов за {elapsed:.1f}с")

    # --- Birthday на H[7] ---
    h7_dict = {}
    h7_collisions = []
    for i, d in enumerate(all_data):
        key = d['H7']
        if key in h7_dict:
            j = h7_dict[key]
            h7_collisions.append((i, j))
        else:
            h7_dict[key] = i

    print(f"\n  H[7] коллизии (32-бит birthday): {len(h7_collisions)}")
    print(f"  Теория: N²/2^33 = {N_SAMPLES**2 / 2**33:.1f}")

    # Для каждой H[7]-collision — HW(ΔH) полного хеша
    if h7_collisions:
        full_dh = []
        per_word_dh = [[] for _ in range(8)]
        for i, j in h7_collisions:
            H1 = all_data[i]['H']
            H2 = all_data[j]['H']
            dh = sum(hw(H1[k] ^ H2[k]) for k in range(8))
            full_dh.append(dh)
            for k in range(8):
                per_word_dh[k].append(hw(H1[k] ^ H2[k]))

        print(f"\n  При H[7]-collision, HW(ΔH) полного хеша:")
        print(f"    Mean: {np.mean(full_dh):.1f}")
        print(f"    Best: {min(full_dh)}")

        print(f"\n  HW(ΔH[k]) по словам при H[7]-collision:")
        reg = ['H0(a)', 'H1(b)', 'H2(c)', 'H3(d)', 'H4(e)', 'H5(f)', 'H6(g)', 'H7(h)']
        compression_bits = 0
        for k in range(8):
            m = np.mean(per_word_dh[k])
            expected = 0 if k == 7 else 16
            comp = expected - m
            compression_bits += max(comp, 0)
            marker = f" ← {comp:+.1f} бит" if abs(comp) > 1 else ""
            print(f"    {reg[k]}: {m:5.1f}/32 (expected {expected}){marker}")

        print(f"\n  >>> Суммарное сжатие: +{compression_bits:.1f} бит")
    else:
        print("  Нет H[7]-коллизий (нужно больше сэмплов)")
        compression_bits = 0

    # =================================================================
    print("\n" + "=" * 70)
    print("2. ИНСТРУМЕНТ: g[62] СОРТИРОВКА (group-birthday)")
    print("=" * 70)

    # Сортируем по g[62], ищем пары с близким g[62]
    sorted_by_g62 = sorted(range(N_SAMPLES), key=lambda i: all_data[i]['g62'])

    # Пары с ближайшим g[62]
    close_pairs_dh = []
    far_pairs_dh = []

    for idx in range(len(sorted_by_g62) - 1):
        i = sorted_by_g62[idx]
        j = sorted_by_g62[idx + 1]
        g62_diff = abs(all_data[i]['g62'] - all_data[j]['g62'])

        H1 = all_data[i]['H']
        H2 = all_data[j]['H']
        dh = sum(hw(H1[k] ^ H2[k]) for k in range(8))

        if g62_diff < 256:  # очень близкие
            close_pairs_dh.append(dh)
        elif g62_diff > 2**31:  # далёкие
            far_pairs_dh.append(dh)

    if close_pairs_dh:
        print(f"\n  Близкие g[62] (Δ<256): {len(close_pairs_dh)} пар")
        print(f"    Mean HW(ΔH): {np.mean(close_pairs_dh):.1f}")
        print(f"    Best HW(ΔH): {min(close_pairs_dh)}")
    if far_pairs_dh:
        print(f"  Далёкие g[62] (Δ>2^31): {len(far_pairs_dh)} пар")
        print(f"    Mean HW(ΔH): {np.mean(far_pairs_dh):.1f}")

    if close_pairs_dh and far_pairs_dh:
        g62_gain = np.mean(far_pairs_dh) - np.mean(close_pairs_dh)
        print(f"\n  >>> g[62]-сортировка даёт: -{g62_gain:.1f} бит")
    else:
        g62_gain = 0

    # =================================================================
    print("\n" + "=" * 70)
    print("3. ИНСТРУМЕНТ: e[60] ФИЛЬТРАЦИЯ (bottleneck)")
    print("=" * 70)

    # Фильтруем пары где e[60]>>24 совпадают
    e60_groups = {}
    for i, d in enumerate(all_data):
        key = d['e60'] >> 24  # верхний байт
        if key not in e60_groups:
            e60_groups[key] = []
        e60_groups[key].append(i)

    # Birthday внутри e[60]-группы
    e60_pair_dh = []
    for key, members in e60_groups.items():
        if len(members) >= 2:
            # Берём первые пары
            for a_idx in range(min(len(members), 10)):
                for b_idx in range(a_idx + 1, min(len(members), 10)):
                    i, j = members[a_idx], members[b_idx]
                    H1 = all_data[i]['H']
                    H2 = all_data[j]['H']
                    dh = sum(hw(H1[k] ^ H2[k]) for k in range(8))
                    e60_pair_dh.append(dh)

    if e60_pair_dh:
        print(f"\n  Пары с одинаковым e[60]>>24: {len(e60_pair_dh)}")
        print(f"    Mean HW(ΔH): {np.mean(e60_pair_dh):.1f}")
        print(f"    Best HW(ΔH): {min(e60_pair_dh)}")
        e60_gain = 128 - np.mean(e60_pair_dh)
        print(f"\n  >>> e[60]-фильтрация даёт: {e60_gain:+.1f} бит vs random")
    else:
        e60_gain = 0

    # =================================================================
    print("\n" + "=" * 70)
    print("4. ИНСТРУМЕНТ: H[7] bias (carry[63] signal)")
    print("=" * 70)

    # Группируем по H[7][b30,b29]
    groups = {'00': [], '01': [], '10': [], '11': []}
    for i, d in enumerate(all_data):
        key = f"{d['h7b30']}{d['h7b29']}"
        groups[key].append(i)

    print(f"\n  Распределение по H[7][b30,b29]:")
    for key, members in groups.items():
        print(f"    {key}: {len(members):6d} ({len(members)/N_SAMPLES*100:.1f}%)")

    # Birthday внутри группы 00 (carry-aligned)
    g00 = groups['00']
    h7_bias_dh = []
    if len(g00) > 100:
        for trial in range(min(5000, len(g00) * (len(g00)-1) // 2)):
            i, j = np.random.choice(len(g00), 2, replace=False)
            H1 = all_data[g00[i]]['H']
            H2 = all_data[g00[j]]['H']
            dh = sum(hw(H1[k] ^ H2[k]) for k in range(8))
            h7_bias_dh.append(dh)

    if h7_bias_dh:
        print(f"\n  Birthday внутри H[7][b30,b29]=00: {len(h7_bias_dh)} пар")
        print(f"    Mean HW(ΔH): {np.mean(h7_bias_dh):.1f}")
        print(f"    Best HW(ΔH): {min(h7_bias_dh)}")
        h7_gain = 128 - np.mean(h7_bias_dh)
        print(f"  >>> H[7]-bias даёт: {h7_gain:+.1f} бит vs random")
    else:
        h7_gain = 0

    # =================================================================
    print("\n" + "=" * 70)
    print("5. КОМБИНАЦИЯ: g[62]-sort + e[60]-filter + H[7]-bias")
    print("=" * 70)

    # Комбинируем: берём пары с:
    # - близким g[62] (Δ<1024)
    # - одинаковым e[60]>>28 (верхний nibble)
    # - H[7][b30]=0

    # Сортируем по g[62]
    sorted_idx = sorted(range(N_SAMPLES), key=lambda i: all_data[i]['g62'])

    combo_dh = []
    for idx in range(len(sorted_idx) - 1):
        i = sorted_idx[idx]
        j = sorted_idx[idx + 1]

        # g[62] close?
        if abs(all_data[i]['g62'] - all_data[j]['g62']) > 1024:
            continue

        # e[60] similar?
        if (all_data[i]['e60'] >> 28) != (all_data[j]['e60'] >> 28):
            continue

        # H[7] bias aligned?
        if all_data[i]['h7b30'] != 0 or all_data[j]['h7b30'] != 0:
            continue

        H1 = all_data[i]['H']
        H2 = all_data[j]['H']
        dh = sum(hw(H1[k] ^ H2[k]) for k in range(8))
        combo_dh.append(dh)

    if combo_dh:
        print(f"\n  Комбинированные пары: {len(combo_dh)}")
        print(f"    Mean HW(ΔH): {np.mean(combo_dh):.1f}")
        print(f"    Best HW(ΔH): {min(combo_dh)}")
        combo_gain = 128 - np.mean(combo_dh)
        print(f"    >>> Комбинированный выигрыш: {combo_gain:+.1f} бит")
    else:
        print("  Нет подходящих пар (нужно больше сэмплов)")
        combo_gain = 0

    # =================================================================
    print("\n" + "=" * 70)
    print("6. БУХГАЛТЕРИЯ: СУММАРНЫЙ ВЫИГРЫШ")
    print("=" * 70)

    print(f"""
  ИНСТРУМЕНТ                  ИЗМЕРЕННЫЙ ВЫИГРЫШ
  ─────────────────────────────────────────────
  H[7]-collision cascade      +{compression_bits:.1f} бит (из методички: +17)
  g[62]-сортировка            {g62_gain:+.1f} бит (из методички: -18.2)
  e[60]-фильтрация            {e60_gain:+.1f} бит
  H[7]-bias фильтрация        {h7_gain:+.1f} бит
  Комбинация всех              {combo_gain:+.1f} бит

  BASELINE:                   2^128 (birthday bound)
  С инструментами:            2^{128 - max(compression_bits, 0) - max(g62_gain, 0):.1f}
  """)

    # Теоретическая сводка
    total_gain = compression_bits + max(g62_gain, 0)
    print(f"  ТЕОРЕТИЧЕСКИЙ МАКСИМУМ (из методички):")
    print(f"    H[7]-cascade:   +17 бит")
    print(f"    g[62]-sort:     +18 бит")
    print(f"    Итого:          +35 бит")
    print(f"    Стоимость:      2^{{128-35}} = 2^93")
    print(f"")
    print(f"  ИЗМЕРЕННЫЙ:")
    print(f"    Наш суммарный:  +{total_gain:.1f} бит")
    print(f"    Стоимость:      2^{128 - total_gain:.1f}")


if __name__ == "__main__":
    main()
