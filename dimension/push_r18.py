"""
ПРОБИВАЕМ R=18: минимизация δW[18].

δW[18] = σ₁(δW[16]) + δW[11] + σ₀(δW[3]) + δW[2]
δW[2] = 0 (by design)

Три слагаемых. Каждое зависит от a-repair.
Идея: a-repair задаёт δW[4..15]. Но δW[9,11,14] входят в schedule.
Может для НЕКОТОРЫХ W1 a-repair создаёт δW[9]=δW[11]=δW[14]=0?

Тогда δW[16] = 0, δW[18] = σ₀(δW[3]) ≈ 2-3 бита.
Или даже δW[18] = 0 если σ₀(δW[3]) тоже компенсируется?

Стратегия: перебираем W1, для каждого запускаем a-repair,
смотрим HW(δW[9]), HW(δW[11]), HW(δW[14]).
Ищем W1 где все три малы.
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


def full_a_repair_analysis(W1, break_round=3, break_bit=31):
    """A-repair + полный анализ schedule."""
    W1f = expand_schedule(W1)

    s1 = tuple(IV)
    states1 = [s1]
    for r in range(64):
        s1 = R(s1, W1f[r], r)
        states1.append(s1)

    W2 = list(W1)
    W2[break_round] ^= (1 << break_bit)

    s2 = tuple(IV)
    for r in range(break_round):
        s2 = R(s2, W2[r], r)
    s2 = R(s2, W2[break_round], break_round)

    for r in range(break_round + 1, 16):
        W2[r] = a_repair_W(s2, states1[r + 1][0], r)
        s2 = R(s2, W2[r], r)

    W2f = expand_schedule(W2)

    # δW per word
    dw = [hw(W1f[r] ^ W2f[r]) for r in range(64)]
    dw_raw = [W1f[r] ^ W2f[r] for r in range(64)]

    # Forward profile
    s2_full = tuple(IV)
    profile = []
    for r in range(64):
        w = W2[r] if r < 16 else W2f[r]
        s2_full = R(s2_full, w, r)
        ds = state_diff(states1[r + 1], s2_full)
        profile.append(ds)

    H1 = tuple(add32(IV[i], states1[64][i]) for i in range(8))
    H2 = tuple(add32(IV[i], s2_full[i]) for i in range(8))
    dh = sum(hw(H1[i] ^ H2[i]) for i in range(8))

    return dw, dw_raw, profile, dh


def main():
    np.random.seed(42)

    print("=" * 70)
    print("ПРОБИВАЕМ R=18: поиск W1 с минимальным δW[18]")
    print("=" * 70)

    # =================================================================
    print(f"\n{'=' * 70}")
    print("1. РАСПРЕДЕЛЕНИЕ δW[9], δW[11], δW[14] от a-repair")
    print("=" * 70)

    N = 5000
    dw9s, dw11s, dw14s, dw16s, dw17s, dw18s = [], [], [], [], [], []
    dw18_zeros = 0
    best_dw18 = 32
    best_W1_for_dw18 = None

    for trial in range(N):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        dw, dw_raw, profile, dh = full_a_repair_analysis(W1)

        dw9s.append(dw[9])
        dw11s.append(dw[11])
        dw14s.append(dw[14])
        dw16s.append(dw[16])
        dw17s.append(dw[17])
        dw18s.append(dw[18])

        if dw[18] == 0:
            dw18_zeros += 1
        if dw[18] < best_dw18:
            best_dw18 = dw[18]
            best_W1_for_dw18 = list(W1)

    print(f"\n  {N} random W1, break r=3 bit=31:")
    print(f"  {'Word':>7} {'Mean HW':>9} {'P(=0)':>8} {'Min':>5}")
    for name, vals in [('δW[9]', dw9s), ('δW[11]', dw11s), ('δW[14]', dw14s),
                       ('δW[16]', dw16s), ('δW[17]', dw17s), ('δW[18]', dw18s)]:
        p0 = sum(1 for v in vals if v == 0) / N * 100
        print(f"  {name:>7} {np.mean(vals):8.2f} {p0:7.1f}% {min(vals):4d}")

    print(f"\n  Best δW[18] = {best_dw18}")
    print(f"  P(δW[18]=0) = {dw18_zeros/N*100:.2f}%")

    # Joint probabilities
    both_9_14_zero = sum(1 for i in range(N) if dw9s[i] == 0 and dw14s[i] == 0) / N * 100
    all_three_zero = sum(1 for i in range(N) if dw9s[i] == 0 and dw11s[i] == 0 and dw14s[i] == 0) / N * 100
    print(f"\n  P(δW[9]=0 AND δW[14]=0) = {both_9_14_zero:.2f}%")
    print(f"  P(δW[9]=0 AND δW[11]=0 AND δW[14]=0) = {all_three_zero:.2f}%")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("2. ФИЛЬТРАЦИЯ: W1 с малым δW[18]")
    print("=" * 70)

    # Среди 5K: сколько имеют δW[18] ≤ 2?
    small_dw18 = [(i, dw18s[i]) for i in range(N) if dw18s[i] <= 4]
    print(f"\n  δW[18] ≤ 4: {len(small_dw18)} / {N} = {len(small_dw18)/N*100:.2f}%")

    if small_dw18:
        # Для этих: какой r=18 profile?
        for idx, dw18_val in small_dw18[:10]:
            W1 = [np.random.randint(0, 2**32) for _ in range(16)]  # need to regenerate
            # Пересобираем с тем же seed
            pass

    # =================================================================
    print(f"\n{'=' * 70}")
    print("3. МАССОВЫЙ ПОИСК: δW[18]=0")
    print("=" * 70)

    # Перебираем W1, ищем δW[18]=0 после a-repair
    found_18 = 0
    found_1718 = 0
    found_161718 = 0
    best_chain = 0
    best_profile = None

    for trial in range(50000):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        dw, dw_raw, profile, dh = full_a_repair_analysis(W1)

        # Longest BOTH=0 chain from start
        chain = 0
        for r in range(64):
            if profile[r] == 0:
                chain = r + 1
            else:
                break

        # Check schedule transparency
        if dw[18] == 0:
            found_18 += 1
        if dw[17] == 0 and dw[18] == 0:
            found_1718 += 1
        if dw[16] == 0 and dw[17] == 0 and dw[18] == 0:
            found_161718 += 1

        # Best chain after reboot
        all_zeros = [r for r in range(64) if profile[r] == 0]
        max_round_zero = max(all_zeros) if all_zeros else 0

        if max_round_zero > best_chain:
            best_chain = max_round_zero
            best_profile = (profile[:25], dw[:25], dh, trial)

        if trial % 10000 == 9999:
            print(f"    {trial+1}K: δW[18]=0: {found_18}, "
                  f"δW[17,18]=0: {found_1718}, "
                  f"δW[16,17,18]=0: {found_161718}, "
                  f"best chain to r={best_chain}")

    print(f"\n  ИТОГО за 50K:")
    print(f"    P(δW[18]=0):          {found_18} = {found_18/50000*100:.3f}%")
    print(f"    P(δW[17,18]=0):       {found_1718} = {found_1718/50000*100:.3f}%")
    print(f"    P(δW[16,17,18]=0):    {found_161718} = {found_161718/50000*100:.3f}%")
    print(f"    Best BOTH=0 to round: {best_chain}")

    if best_profile:
        prof, dw_prof, dh_best, trial_best = best_profile
        print(f"\n  Лучший профиль (trial {trial_best}):")
        print(f"    {'r':>4} {'δstate':>8} {'δW':>6}")
        for r in range(min(25, len(prof))):
            marker = " ★" if prof[r] == 0 else ""
            print(f"    {r:4d} {prof[r]:7d} {dw_prof[r]:5d}{marker}")
        print(f"    HW(δH) = {dh_best}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("4. СТОИМОСТЬ В НАШЕМ ИЗМЕРЕНИИ")
    print("=" * 70)

    r17_rate = sum(1 for i in range(N) if dw16s[i] == 0 and dw17s[i] == 0) / N
    r18_rate = found_161718 / 50000 if found_161718 > 0 else 0

    print(f"""
  A-REPAIR + BREAK BIT=31 + W1 SEARCH:

    P(δW[16]=0): ~50%     → стоимость 2^1
    P(δW[17]=0): ~50%     → стоимость 2^1
    P(δW[18]=0): {found_18/50000*100:.3f}%   → стоимость 2^{-np.log2(max(found_18/50000, 1e-10)):.1f}
    P(все три=0): {found_161718/50000*100:.3f}% → стоимость 2^{-np.log2(max(found_161718/50000, 1e-10)):.1f}

    Schedule прозрачен до r=17 за ~4 попытки (2^2).
    Расширение до r=18: {'возможно' if found_18 > 0 else 'не найдено'} (стоимость 2^{-np.log2(max(found_18/50000, 1e-10)):.0f}).

    Если r=18 пройден → gap начинается с r=19.
    Forward fix до r=18 + backward 4 = 22 раунда controlled.
    Gap: 64 - 22 = 42 раунда (vs Wang 48, vs basic a-repair 44).
""")


if __name__ == "__main__":
    main()
