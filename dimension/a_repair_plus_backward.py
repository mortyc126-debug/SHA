"""
A-REPAIR + BACKWARD ASYMMETRY: комбинированная стратегия.

A-repair даёт:
  - Reboot на r=12 (100%)
  - BOTH=0 на r=1-3, 12-16 (8 раундов)
  - δW[16] = 1.9 бит (малый!)
  - r=17 BOTH=0 = 20.8%

Backward даёт:
  - r=61-64: δstate < 50 (4 раунда, стабильно)
  - Backward СТАБИЛЕН (не зависит от forward)

Комбинация:
  Front:  a-repair r=1-16 (+ 20.8% проход r=17)
  Back:   backward asymmetry r=61-64
  Gap:    r=17(18)..60 = 43(42) раунда
  vs Wang: gap = 48 раундов
  vs a-repair only: gap = 44 раунда

Ключевой вопрос: a-repair + backward = сколько TOTAL controlled?
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

def R_inv(sn, W_r, r_idx):
    ap, bp, cp, dp, ep, fp, gp, hp = sn
    a, b, c, e, f, g = bp, cp, dp, fp, gp, hp
    T2 = add32(Sigma0(a), Maj(a, b, c))
    T1 = sub32(ap, T2)
    d = sub32(ep, T1)
    h = sub32(sub32(sub32(sub32(T1, Sigma1(e)), Ch(e, f, g)), K[r_idx]), W_r)
    return (a, b, c, d, e, f, g, h)

def expand_schedule(W16):
    W = list(W16)
    for r in range(16, 64):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))
    return W

def state_diff(s1, s2):
    return sum(hw(s1[i] ^ s2[i]) for i in range(8))

def a_repair_W(state2, target_a, r_idx):
    a, b, c, d, e, f, g, h = state2
    T2 = add32(Sigma0(a), Maj(a, b, c))
    T1_needed = sub32(target_a, T2)
    return sub32(sub32(sub32(sub32(T1_needed, h), Sigma1(e)), Ch(e, f, g)), K[r_idx])


def main():
    np.random.seed(42)

    print("=" * 70)
    print("A-REPAIR + BACKWARD ASYMMETRY")
    print("=" * 70)

    N = 1000

    # Полный профиль: forward a-repair + backward trace
    fwd_profiles = np.zeros((N, 65))
    bwd_profiles = np.zeros((N, 65))
    combined = np.zeros((N, 65))
    both_zero = np.zeros((N, 65))
    r17_pass = 0

    for trial in range(N):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W1f = expand_schedule(W1)

        # Trace 1: full forward
        s1 = tuple(IV)
        states1 = [s1]
        for r in range(64):
            s1 = R(s1, W1f[r], r)
            states1.append(s1)

        # Trace 2: a-repair (break r=3)
        W2 = list(W1)
        W2[3] ^= (1 << (trial % 32))

        s2 = tuple(IV)
        for r in range(3):
            s2 = R(s2, W2[r], r)
        s2 = R(s2, W2[3], 3)  # break
        for r in range(4, 16):
            W2[r] = a_repair_W(s2, states1[r + 1][0], r)
            s2 = R(s2, W2[r], r)

        W2f = expand_schedule(W2)

        # Forward full
        states2_fwd = list(states1[:4])  # r=0..3 same (except we need actual)
        s2_fwd = tuple(IV)
        for r in range(64):
            w = W2[r] if r < 16 else W2f[r]
            s2_fwd = R(s2_fwd, w, r)
            if r >= 3:
                states2_fwd.append(s2_fwd)

        # Rebuild states2 properly
        s2_full = tuple(IV)
        s2_states = [s2_full]
        for r in range(64):
            w = W2[r] if r < 16 else W2f[r]
            s2_full = R(s2_full, w, r)
            s2_states.append(s2_full)

        # Forward profile
        for r in range(65):
            fwd_profiles[trial, r] = state_diff(states1[r], s2_states[r])
            if fwd_profiles[trial, r] == 0:
                both_zero[trial, r] = 1

        # r=17 pass?
        if fwd_profiles[trial, 17] == 0:
            r17_pass += 1

        # Backward: from states1[64] with W2f
        sb = states1[64]
        bwd_states = [None] * 65
        bwd_states[64] = sb
        for r in range(63, -1, -1):
            sb = R_inv(sb, W2f[r], r)
            bwd_states[r] = sb

        for r in range(65):
            bwd_profiles[trial, r] = state_diff(states1[r], bwd_states[r])

        # Combined: min of forward and backward δ at each round
        for r in range(65):
            combined[trial, r] = min(fwd_profiles[trial, r], bwd_profiles[trial, r])

    # =================================================================
    print(f"\n{'=' * 70}")
    print("1. ПОЛНЫЙ ПРОФИЛЬ: forward + backward + combined")
    print("=" * 70)

    mean_fwd = np.mean(fwd_profiles, axis=0)
    mean_bwd = np.mean(bwd_profiles, axis=0)
    mean_comb = np.mean(combined, axis=0)
    mean_both = np.mean(both_zero, axis=0) * 100

    print(f"\n  {'r':>4} {'Fwd':>6} {'Bwd':>6} {'Min':>6} {'BOTH%':>7}  Visual")
    for r in range(65):
        f, b, c, bz = mean_fwd[r], mean_bwd[r], mean_comb[r], mean_both[r]
        # Visual: F=forward zone, B=backward zone, G=gap
        if f < 2:
            zone = "F"
        elif b < 50:
            zone = "B"
        elif c < 64:
            zone = "~"
        else:
            zone = "G"

        bar_f = "█" * min(int(f / 8), 16)
        bar_b = "░" * min(int(b / 8), 16)

        if r <= 20 or r >= 58 or r in [25, 30, 40, 50]:
            print(f"  {r:4d} {f:5.1f} {b:5.1f} {c:5.1f} {bz:6.1f}%  [{zone}] {bar_f}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("2. ЗОНЫ КОНТРОЛЯ")
    print("=" * 70)

    fwd_fixed = sum(1 for r in range(65) if mean_fwd[r] < 2)
    bwd_ctrl = sum(1 for r in range(65) if mean_bwd[r] < 50)
    gap_rounds = sum(1 for r in range(65) if mean_fwd[r] >= 2 and mean_bwd[r] >= 50)

    print(f"\n  Forward fixed (δ<2):     {fwd_fixed} раундов")
    print(f"  Backward controlled (δ<50): {bwd_ctrl} раундов")
    print(f"  Gap (обе > threshold):    {gap_rounds} раундов")
    print(f"  r=17 BOTH=0:              {r17_pass}/{N} = {r17_pass/N*100:.1f}%")

    # r=17 pass: what happens AFTER?
    print(f"\n  При r=17 BOTH=0: что дальше?")
    r17_pass_profiles = [fwd_profiles[t] for t in range(N) if fwd_profiles[t, 17] == 0]
    if r17_pass_profiles:
        r17_arr = np.array(r17_pass_profiles)
        for r in [17, 18, 19, 20, 24, 32]:
            m = np.mean(r17_arr[:, r])
            p = np.mean(r17_arr[:, r] == 0) * 100
            print(f"    r={r}: δstate={m:.1f}, BOTH=0={p:.1f}%")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("3. СВОДКА: что даёт комбинация")
    print("=" * 70)

    print(f"""
  СТРАТЕГИЯ               FwdFixed  BwdCtrl  Gap    r17pass
  ─────────────────────── ──────── ──────── ───── ─────────
  Wang (forward fix)         17       3.6    43.4    0%
  a-repair (break=3)          8       3.6    52.4   {r17_pass/N*100:.1f}%
  a-repair + backward         8       3.6    52.4   {r17_pass/N*100:.1f}%

  a-repair r=17 PASS:
    {r17_pass/N*100:.1f}% пар проходят schedule barrier!
    Это {r17_pass} из {N} пар.
    Стоимость одного прохода: ~{int(N/max(r17_pass,1))} попыток = 2^{np.log2(max(N/max(r17_pass,1),1)):.1f}

  Если r=17 пройден, следующий barrier = r=18:
    δW[17] зависит от σ₁(δW[15]) + δW[10] + σ₀(δW[2]) + δW[1]
    δW[2]=0 (a-repair), δW[1]=0 (a-repair)
    → δW[17] = σ₁(δW[15]) + δW[10]
    Два ненулевых слагаемых → δW[17] ≈ random
    → P(δe[18]=0 | δe[17]=0) ≈ 2^{{-32}}
""")

    # =================================================================
    print(f"{'=' * 70}")
    print("4. КАСКАД SCHEDULE: какие δW[16..] малы при a-repair?")
    print("=" * 70)

    dw_profiles = np.zeros((N, 64))
    for trial in range(N):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W1f = expand_schedule(W1)

        W2 = list(W1)
        W2[3] ^= (1 << (trial % 32))

        s1 = tuple(IV)
        states1_local = [s1]
        for r in range(64):
            s1 = R(s1, W1f[r], r)
            states1_local.append(s1)

        s2 = tuple(IV)
        for r in range(3):
            s2 = R(s2, W2[r], r)
        s2 = R(s2, W2[3], 3)
        for r in range(4, 16):
            W2[r] = a_repair_W(s2, states1_local[r + 1][0], r)
            s2 = R(s2, W2[r], r)

        W2f = expand_schedule(W2)
        for r in range(64):
            dw_profiles[trial, r] = hw(W1f[r] ^ W2f[r])

    mean_dw = np.mean(dw_profiles, axis=0)
    p_zero_dw = np.mean(dw_profiles == 0, axis=0) * 100

    print(f"\n  {'r':>4} {'HW(δW)':>8} {'P(δW=0)':>9}")
    for r in range(16, 25):
        print(f"  {r:4d} {mean_dw[r]:7.1f} {p_zero_dw[r]:8.1f}%")

    print(f"\n  δW[0..2]=0 по дизайну (a-repair)")
    print(f"  δW[3]=1 бит (break)")
    print(f"  δW[4..15]: a-repair подбирает → schedule зависит")

    # W[16] formula: sig1(W[14]) + W[9] + sig0(W[1]) + W[0]
    # δW[16] = sig1(δW[14]) + δW[9] + sig0(δW[1]) + δW[0]
    # δW[0]=0, δW[1]=0, δW[9] и δW[14] = a-repair values
    print(f"\n  δW[16] = σ₁(δW[14]) + δW[9] + σ₀(δW[1]) + δW[0]")
    print(f"         = σ₁(δW[14]) + δW[9] + 0 + 0")
    print(f"  δW[9] и δW[14] — из a-repair. Обычно малые.")
    print(f"  Mean HW(δW[16]) = {mean_dw[16]:.1f}")
    print(f"  P(δW[16]=0) = {p_zero_dw[16]:.1f}%")


if __name__ == "__main__":
    main()
