"""
МОСТИК В ХАОС: не обнуляем δ, а ТОРМОЗИМ.

Текущая ситуация:
  r=0-17:  δstate=0 (a-repair + schedule transparency, cost ~2^2)
  r=18:    δstate=10 (δW[18]≈4-8 бит, неустранимый)
  r=19-21: δstate взрывается 0→32→64→128 за 4 раунда
  r=22-60: δstate≈128 (плато, хаос)
  r=61-64: backward δ < 50

Идея: вместо обнуления r=18, ЗАМЕДЛИТЬ рост.
Если δstate[18]=4 бит вместо 10 → насыщение на r=23 вместо r=22
Каждый спасённый раунд = меньше хаоса к r=60.

Эксперимент:
  1. При каких W1 рост МЕДЛЕННЕЕ? (δstate[22] < 128)
  2. Связан ли малый δW[18] с медленным ростом?
  3. Сколько стоит "торможение" на 1 раунд?
  4. Связь δstate на meeting point (r=60) с δH
"""

import numpy as np
import struct, hashlib

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

def expand_schedule(W16):
    W = list(W16)
    for r in range(16, 64):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))
    return W

def a_repair_W(state2, target_a, r_idx):
    a, b, c, d, e, f, g, h = state2
    T2 = add32(Sigma0(a), Maj(a, b, c))
    T1_needed = sub32(target_a, T2)
    return sub32(sub32(sub32(sub32(T1_needed, h), Sigma1(e)), Ch(e, f, g)), K[r_idx])

def state_diff(s1, s2):
    return sum(hw(s1[i] ^ s2[i]) for i in range(8))


def run_full(W1, break_bit=31):
    """A-repair + full forward, return profile and schedule."""
    W1f = expand_schedule(W1)
    s1 = tuple(IV)
    states1 = [s1]
    for r in range(64):
        s1 = R(s1, W1f[r], r)
        states1.append(s1)

    W2 = list(W1)
    W2[3] ^= (1 << break_bit)
    s2 = tuple(IV)
    for r in range(3):
        s2 = R(s2, W2[r], r)
    s2 = R(s2, W2[3], 3)
    for r in range(4, 16):
        W2[r] = a_repair_W(s2, states1[r + 1][0], r)
        s2 = R(s2, W2[r], r)

    W2f = expand_schedule(W2)
    s2_full = tuple(IV)
    profile = []
    for r in range(64):
        w = W2[r] if r < 16 else W2f[r]
        s2_full = R(s2_full, w, r)
        profile.append(state_diff(states1[r + 1], s2_full))

    dw = [hw(W1f[r] ^ W2f[r]) for r in range(64)]
    H1 = tuple(add32(IV[i], states1[64][i]) for i in range(8))
    H2 = tuple(add32(IV[i], s2_full[i]) for i in range(8))
    dh = sum(hw(H1[i] ^ H2[i]) for i in range(8))

    return profile, dw, dh


def main():
    np.random.seed(42)

    print("=" * 70)
    print("МОСТИК В ХАОС: торможение вместо обнуления")
    print("=" * 70)

    # =================================================================
    print(f"\n{'=' * 70}")
    print("1. РАСПРЕДЕЛЕНИЕ δstate[18] и скорость роста")
    print("=" * 70)

    N = 10000
    ds18_list = []
    ds20_list = []
    ds25_list = []
    ds30_list = []
    dw18_list = []
    dh_list = []
    profiles_all = []

    for trial in range(N):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        profile, dw, dh = run_full(W1)
        ds18_list.append(profile[17])  # index 17 = after round 18
        ds20_list.append(profile[19])
        ds25_list.append(profile[24])
        ds30_list.append(profile[29])
        dw18_list.append(dw[18])
        dh_list.append(dh)
        profiles_all.append(profile)

    print(f"\n  {N} trials, a-repair break=3 bit=31:")
    print(f"  {'Metric':>12} {'Mean':>8} {'Std':>7} {'Min':>6} {'5th%':>7} {'Median':>7}")
    for name, vals in [('δstate[18]', ds18_list), ('δstate[20]', ds20_list),
                       ('δstate[25]', ds25_list), ('δstate[30]', ds30_list),
                       ('δW[18]', dw18_list), ('HW(δH)', dh_list)]:
        arr = np.array(vals)
        print(f"  {name:>12} {np.mean(arr):7.1f} {np.std(arr):6.1f} {np.min(arr):5d} "
              f"{np.percentile(arr, 5):6.1f} {np.median(arr):6.1f}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("2. КОРРЕЛЯЦИЯ: малый δW[18] → медленный рост?")
    print("=" * 70)

    corr_18_20 = np.corrcoef(dw18_list, ds20_list)[0, 1]
    corr_18_25 = np.corrcoef(dw18_list, ds25_list)[0, 1]
    corr_18_30 = np.corrcoef(dw18_list, ds30_list)[0, 1]
    corr_18_dh = np.corrcoef(dw18_list, dh_list)[0, 1]
    corr_ds18_dh = np.corrcoef(ds18_list, dh_list)[0, 1]

    print(f"\n  corr(δW[18], δstate[20]):  {corr_18_20:+.4f}")
    print(f"  corr(δW[18], δstate[25]):  {corr_18_25:+.4f}")
    print(f"  corr(δW[18], δstate[30]):  {corr_18_30:+.4f}")
    print(f"  corr(δW[18], HW(δH)):     {corr_18_dh:+.4f}")
    print(f"  corr(δstate[18], HW(δH)): {corr_ds18_dh:+.4f}")

    # Условно: при δW[18] ≤ 4 (best 7%) vs δW[18] > 10
    small_18 = [i for i in range(N) if dw18_list[i] <= 4]
    large_18 = [i for i in range(N) if dw18_list[i] > 10]

    if small_18 and large_18:
        print(f"\n  При δW[18] ≤ 4 ({len(small_18)} пар):")
        print(f"    δstate[20]: {np.mean([ds20_list[i] for i in small_18]):.1f}")
        print(f"    δstate[25]: {np.mean([ds25_list[i] for i in small_18]):.1f}")
        print(f"    HW(δH):     {np.mean([dh_list[i] for i in small_18]):.1f}, min={min(dh_list[i] for i in small_18)}")

        print(f"\n  При δW[18] > 10 ({len(large_18)} пар):")
        print(f"    δstate[20]: {np.mean([ds20_list[i] for i in large_18]):.1f}")
        print(f"    δstate[25]: {np.mean([ds25_list[i] for i in large_18]):.1f}")
        print(f"    HW(δH):     {np.mean([dh_list[i] for i in large_18]):.1f}, min={min(dh_list[i] for i in large_18)}")

        gain = np.mean([dh_list[i] for i in large_18]) - np.mean([dh_list[i] for i in small_18])
        print(f"\n  Gain from small δW[18]: {gain:+.1f} бит HW(δH)")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("3. СТОИМОСТЬ ТОРМОЖЕНИЯ: δW[18] → δstate рост")
    print("=" * 70)

    # Сколько стоит (в попытках) получить δW[18] ≤ threshold?
    for threshold in [4, 6, 8, 10, 12]:
        count = sum(1 for x in dw18_list if x <= threshold)
        rate = count / N
        cost = 1.0 / max(rate, 1e-10)
        print(f"  P(δW[18] ≤ {threshold:2d}): {rate*100:6.2f}%  cost=2^{np.log2(max(cost,1)):.1f}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("4. КАСКАД ТОРМОЖЕНИЯ: δW[18] мал → δW[19] тоже?")
    print("=" * 70)

    # δW[19] = σ₁(δW[17]) + δW[12] + σ₀(δW[4]) + δW[3]
    # При a-repair: δW[3]=1bit, δW[4]=a-repair, δW[12]=0(?), δW[17]=0.5bit
    # Есть ли каскадное торможение?

    dw19_list = []
    dw20_list = []
    for trial in range(N):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        _, dw, _ = run_full(W1)
        dw19_list.append(dw[19])
        dw20_list.append(dw[20])

    print(f"\n  Schedule cascade:")
    print(f"    δW[18]: {np.mean(dw18_list):.1f} бит")
    print(f"    δW[19]: {np.mean(dw19_list):.1f} бит")
    print(f"    δW[20]: {np.mean(dw20_list):.1f} бит")

    # При малом δW[18] → δW[19] тоже мал?
    corr_18_19 = np.corrcoef(dw18_list[:N], dw19_list)[0, 1]
    print(f"    corr(δW[18], δW[19]): {corr_18_19:+.4f}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("5. БУХГАЛТЕРИЯ: цена каждого раунда r=18..22")
    print("=" * 70)

    # Для каждого раунда: P(δstate[r]=0) и стоимость
    for r in range(17, 25):
        zero_count = sum(1 for p in profiles_all if p[r] == 0)
        small_count = sum(1 for p in profiles_all if p[r] <= 4)
        rate = zero_count / N
        cost = 1.0 / max(rate, 1e-10)

        print(f"  r={r+1}: P(δstate=0)={rate*100:6.2f}% (cost=2^{np.log2(max(cost,1)):.1f}), "
              f"P(≤4)={small_count/N*100:.1f}%, mean={np.mean([p[r] for p in profiles_all]):.1f}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("6. СВОДКА: стоимость в нашем измерении")
    print("=" * 70)

    # Подсчитаем: стоимость прохода через каждый раунд
    costs = {}
    for r in range(17, 25):
        rate = sum(1 for p in profiles_all if p[r] == 0) / N
        costs[r+1] = -np.log2(max(rate, 1/N))

    print(f"\n  Стоимость ОБНУЛЕНИЯ каждого раунда (в битах):")
    total = 0
    for r in sorted(costs):
        total += costs[r]
        print(f"    r={r}: 2^{costs[r]:.1f}  (cumulative: 2^{total:.1f})")

    print(f"\n  Стандартный birthday: 2^128")
    print(f"  A-repair до r=17 + birthday per round: 2^{total:.1f}")

    # Альтернатива: birthday на δstate[meeting_point]
    # Meeting point = где forward δ ≈ backward δ
    # Forward: δstate[60] ≈ 128
    # Backward: δstate[60] ≈ 64
    # Birthday на min(128, 64) × 8 бит... нет, δstate = HW(XOR по 8 словам)

    # Реальная стоимость: birthday на δH (256 бит) = 2^128
    # С a-repair: first 18 rounds free → δ injected at r=18 → 46 rounds gap
    # Growth: 4 бит/раунд на первые 4 раунда, потом saturation
    # Не даёт преимущества для ПОЛНОЙ collision

    # НО: для NEAR-collision (HW(δH) < target):
    near_targets = [120, 110, 100, 90, 80]
    print(f"\n  Near-collision стоимость:")
    for target in near_targets:
        count = sum(1 for d in dh_list if d <= target)
        rate = count / N
        cost = -np.log2(max(rate, 1/N))
        print(f"    HW(δH) ≤ {target}: P={rate*100:.2f}%  cost=2^{cost:.1f}")


if __name__ == "__main__":
    main()
