"""
АНАТОМИЯ GAP-ЗОНЫ (r=17..60): что происходит внутри?

Мы знаем: forward фиксирует r=0..16, backward — r=61..64.
Gap = r=17..60 = 44 раунда. Там δstate растёт 0→128.

Вопросы:
  1. КАК растёт δ внутри gap? Линейно? Экспоненциально? С плато?
  2. ВСЕ ЛИ раунды gap одинаковы? Есть тихие/громкие?
  3. КАКИЕ РЕГИСТРЫ несут δ в gap? a-chain vs e-chain?
  4. δW[16..63] — все ли одинаково вредны? Есть ли тихие W?
  5. Можно ли ЧАСТИЧНО контролировать gap через выбор W[0..15]?
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

def compute_needed_W(state, target_next, r_idx):
    a, b, c, d, e, f, g, h = state
    T2 = add32(Sigma0(a), Maj(a, b, c))
    T1 = sub32(target_next[0], T2)
    return sub32(sub32(sub32(sub32(T1, h), Sigma1(e)), Ch(e, f, g)), K[r_idx])

reg_names = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']


def main():
    np.random.seed(42)
    N = 500

    print("=" * 70)
    print("АНАТОМИЯ GAP-ЗОНЫ (r=17..60)")
    print("=" * 70)

    # Собираем profiles: forward fix (depth=15), потом gap
    all_ds = np.zeros((N, 65))    # δstate total
    all_dreg = np.zeros((N, 65, 8))  # δ per register
    all_dw = np.zeros((N, 64))     # δW schedule
    all_gap_growth = np.zeros((N, 44))  # Δ(δstate) per gap round

    for trial in range(N):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W1f = expand_schedule(W1)

        s1 = tuple(IV)
        states1 = [s1]
        for r in range(64):
            s1 = R(s1, W1f[r], r)
            states1.append(s1)

        W2 = list(W1)
        W2[0] ^= (1 << (trial % 32))
        s2 = tuple(IV)
        s2 = R(s2, W2[0], 0)
        for r in range(1, 16):
            W2[r] = compute_needed_W(s2, states1[r + 1], r)
            s2 = R(s2, W2[r], r)

        W2f = expand_schedule(W2)

        # Forward through gap
        states2 = list(states1[:17])  # forward fix
        for r in range(16, 64):
            s2 = R(s2, W2f[r], r)
            states2.append(s2)

        for r in range(65):
            ds = sum(hw(states1[r][i] ^ states2[r][i]) for i in range(8))
            all_ds[trial, r] = ds
            for i in range(8):
                all_dreg[trial, r, i] = hw(states1[r][i] ^ states2[r][i])

        for r in range(64):
            all_dw[trial, r] = hw(W1f[r] ^ W2f[r])

        for idx, r in enumerate(range(17, 61)):
            all_gap_growth[trial, idx] = all_ds[trial, r] - all_ds[trial, r - 1]

    # =================================================================
    print(f"\n{'=' * 70}")
    print("1. ПРОФИЛЬ GAP: δstate по раундам")
    print("=" * 70)

    mean_ds = np.mean(all_ds, axis=0)
    std_ds = np.std(all_ds, axis=0)

    print(f"\n  {'r':>4} {'δstate':>8} {'±':>5} {'δW[r]':>7} {'Growth':>7}")
    for r in range(16, 61):
        growth = mean_ds[r] - mean_ds[r - 1] if r > 0 else 0
        dw = np.mean(all_dw[:, r])
        bar = "+" * max(0, int(growth))
        neg_bar = "-" * max(0, int(-growth))
        print(f"  {r:4d} {mean_ds[r]:7.1f} {std_ds[r]:4.1f} {dw:6.1f} {growth:+6.1f}  {bar}{neg_bar}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("2. РОСТ δ ПО РЕГИСТРАМ внутри gap")
    print("=" * 70)

    mean_dreg = np.mean(all_dreg, axis=0)

    print(f"\n  {'r':>4}", end="")
    for n in reg_names:
        print(f" δ{n:>3}", end="")
    print(f" {'Σ':>5}  {'Fastest':>8}")

    for r in range(16, 30):
        vals = mean_dreg[r]
        total = sum(vals)
        fastest = reg_names[np.argmax(vals)]
        print(f"  {r:3d}", end="")
        for v in vals:
            marker = "█" if v > 10 else "░" if v > 3 else " "
            print(f" {v:4.1f}{marker}", end="")
        print(f" {total:5.1f}  {fastest:>8}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("3. δW SCHEDULE: тихие и громкие раунды")
    print("=" * 70)

    mean_dw = np.mean(all_dw, axis=0)

    quiet_rounds = []
    loud_rounds = []

    print(f"\n  {'r':>4} {'HW(δW)':>8} {'Тип':>6}")
    for r in range(16, 64):
        dw = mean_dw[r]
        if dw < 4:
            quiet_rounds.append(r)
            typ = "QUIET"
        elif dw > 20:
            loud_rounds.append(r)
            typ = "LOUD"
        else:
            typ = ""
        bar = "█" * int(dw)
        print(f"  {r:4d} {dw:7.1f}  {typ:>5}  {bar}")

    print(f"\n  Тихие (δW<4): {quiet_rounds}")
    print(f"  Громкие (δW>20): {loud_rounds}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("4. КОРРЕЛЯЦИЯ: δW[r] → рост δstate[r]")
    print("=" * 70)

    # Для каждого раунда в gap: corr(HW(δW[r]), growth δstate[r])
    print(f"\n  {'r':>4} {'corr(δW, growth)':>18} {'Значимо?':>10}")
    for idx, r in enumerate(range(17, 55)):
        dw_vals = all_dw[:, r]
        growth_vals = all_gap_growth[:, idx]
        if np.std(dw_vals) > 0 and np.std(growth_vals) > 0:
            corr = np.corrcoef(dw_vals, growth_vals)[0, 1]
        else:
            corr = 0
        sig = "★" if abs(corr) > 0.1 else ""
        print(f"  {r:4d} {corr:+17.4f}  {sig:>8}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("5. ВЫБОР W[0..15]: можно ли сделать gap ТИШЕ?")
    print("=" * 70)

    # Варьируем W1[0] (32 бит свободы) и смотрим δstate[25] (середина gap)
    # При РАЗНЫХ W[0]: schedule РАЗНЫЙ → δW[16..63] РАЗНЫЙ → gap РАЗНЫЙ?

    W_base = [np.random.randint(0, 2**32) for _ in range(16)]

    ds25_by_w0 = []
    dh_by_w0 = []

    for w0_trial in range(2000):
        W1 = list(W_base)
        W1[0] = np.random.randint(0, 2**32)
        W1f = expand_schedule(W1)

        W2 = list(W1)
        W2[0] ^= 0x80000000

        s1 = tuple(IV)
        states1 = [s1]
        for r in range(64):
            s1 = R(s1, W1f[r], r)
            states1.append(s1)

        s2 = tuple(IV)
        s2 = R(s2, W2[0], 0)
        for r in range(1, 16):
            W2[r] = compute_needed_W(s2, states1[r + 1], r)
            s2 = R(s2, W2[r], r)

        W2f = expand_schedule(W2)
        for r in range(16, 64):
            s2 = R(s2, W2f[r], r)

        # δstate at r=25 (mid-gap)
        s1_25 = states1[25]
        s2_fwd = tuple(IV)
        s2_fwd = R(s2_fwd, W2[0], 0)
        for r in range(1, 16):
            s2_fwd = R(s2_fwd, W2[r], r)
        for r in range(16, 25):
            s2_fwd = R(s2_fwd, W2f[r], r)
        ds25 = sum(hw(s1_25[i] ^ s2_fwd[i]) for i in range(8))

        H1 = tuple(add32(IV[i], states1[64][i]) for i in range(8))
        H2 = tuple(add32(IV[i], s2[i]) for i in range(8))
        dh = sum(hw(H1[i] ^ H2[i]) for i in range(8))

        ds25_by_w0.append(ds25)
        dh_by_w0.append(dh)

    print(f"\n  2000 вариантов W[0] (δbit=MSB, forward fix depth=15):")
    print(f"    δstate[25] (mid-gap): mean={np.mean(ds25_by_w0):.1f}, min={min(ds25_by_w0)}, max={max(ds25_by_w0)}")
    print(f"    HW(δH):               mean={np.mean(dh_by_w0):.1f}, min={min(dh_by_w0)}")
    print(f"    corr(δstate[25], δH): {np.corrcoef(ds25_by_w0, dh_by_w0)[0, 1]:+.4f}")

    # Можно ли найти W[0] где gap тише?
    best_idx = np.argmin(ds25_by_w0)
    print(f"\n    Лучший W[0]: δstate[25]={ds25_by_w0[best_idx]}, δH={dh_by_w0[best_idx]}")
    print(f"    Худший:      δstate[25]={max(ds25_by_w0)}")
    print(f"    Spread: {max(ds25_by_w0) - min(ds25_by_w0)}")
    if max(ds25_by_w0) - min(ds25_by_w0) > 20:
        print(f"    → W[0] ВЛИЯЕТ на gap! Можно оптимизировать.")
    else:
        print(f"    → W[0] мало влияет на gap.")


if __name__ == "__main__":
    main()
