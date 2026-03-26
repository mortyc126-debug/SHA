"""
A-REPAIR ОПТИМИЗИРОВАННЫЙ: минимизация δW[18].

Три виновника δW[18] = σ₁(δW[16]) + δW[11] + σ₀(δW[3]):

Оптимизации:
  1. Break bit = бит 0 → σ₀(1) = HW=2 (минимум из всех позиций)
  2. Варьируем W1[0..2] чтобы минимизировать δW[9,11,14] от a-repair
  3. Тестируем разные break rounds (r=3,4,5)
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


def run_a_repair(W1, break_round, break_bit):
    """Полный a-repair: возвращает W2, states, schedule stats."""
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
    for r in range(16, 64):
        s2 = R(s2, W2f[r], r)

    # Full profile
    s2_check = tuple(IV)
    profile = []
    for r in range(64):
        w = W2[r] if r < 16 else W2f[r]
        s2_check = R(s2_check, w, r)
        ds = state_diff(states1[r + 1], s2_check)
        profile.append(ds)

    # Schedule stats
    dw = [hw(W1f[r] ^ W2f[r]) for r in range(64)]

    # Hash diff
    H1 = tuple(add32(IV[i], states1[64][i]) for i in range(8))
    H2 = tuple(add32(IV[i], s2[i]) for i in range(8))
    dh = sum(hw(H1[i] ^ H2[i]) for i in range(8))

    both_zero = sum(1 for d in profile if d == 0)

    return profile, dw, dh, both_zero


def main():
    np.random.seed(42)

    print("=" * 70)
    print("A-REPAIR ОПТИМИЗИРОВАННЫЙ")
    print("=" * 70)

    # =================================================================
    print(f"\n{'=' * 70}")
    print("1. BREAK BIT ОПТИМИЗАЦИЯ: какой бит break минимизирует δW[18]?")
    print("=" * 70)

    break_round = 3
    N = 500

    print(f"\n  Break round={break_round}, варьируем break bit 0..31:")
    print(f"  {'Bit':>4} {'σ₀ HW':>6} {'δW16':>6} {'δW17':>6} {'δW18':>6} {'r17%':>6} {'r18%':>6} {'BOTH0':>6} {'minδH':>6}")

    best_bit = 0
    best_r17 = 0

    for break_bit in range(32):
        r17_pass = 0
        r18_pass = 0
        dw16s = []
        dw17s = []
        dw18s = []
        both_zeros = []
        dhs = []

        for trial in range(N):
            W1 = [np.random.randint(0, 2**32) for _ in range(16)]
            profile, dw, dh, bz = run_a_repair(W1, break_round, break_bit)

            dw16s.append(dw[16])
            dw17s.append(dw[17])
            dw18s.append(dw[18])
            both_zeros.append(bz)
            dhs.append(dh)

            if profile[16] == 0 and profile[17] == 0:
                r17_pass += 1
                if profile[18] == 0:
                    r18_pass += 1

        r17p = r17_pass / N * 100
        r18p = r18_pass / N * 100
        sig0_hw = hw(sigma0(1 << break_bit))

        marker = ""
        if r17p > best_r17:
            best_r17 = r17p
            best_bit = break_bit
            marker = " ★"

        if break_bit < 8 or break_bit > 28 or r17p > 15:
            print(f"  {break_bit:4d} {sig0_hw:5d} {np.mean(dw16s):5.1f} {np.mean(dw17s):5.1f} "
                  f"{np.mean(dw18s):5.1f} {r17p:5.1f} {r18p:5.1f} {np.mean(both_zeros):5.1f} "
                  f"{min(dhs):5d}{marker}")

    print(f"\n  ЛУЧШИЙ break bit: {best_bit} (r17 pass = {best_r17:.1f}%)")

    # =================================================================
    print(f"\n{'=' * 70}")
    print(f"2. ДЕТАЛЬНО: break bit={best_bit}")
    print("=" * 70)

    profiles_best = []
    dw_best = []
    dh_best = []
    for trial in range(2000):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        profile, dw, dh, bz = run_a_repair(W1, break_round, best_bit)
        profiles_best.append(profile)
        dw_best.append(dw)
        dh_best.append(dh)

    mean_profile = np.mean(profiles_best, axis=0)
    both_pct = [np.mean([p[r] == 0 for p in profiles_best]) * 100 for r in range(64)]
    mean_dw = np.mean(dw_best, axis=0)

    print(f"\n  {'r':>4} {'δstate':>8} {'BOTH%':>7} {'δW':>6}")
    for r in range(25):
        status = ""
        if both_pct[r] > 99: status = " ★ZERO"
        elif both_pct[r] > 10: status = " ~pass"
        print(f"  {r:4d} {mean_profile[r]:7.1f} {both_pct[r]:6.1f}% {mean_dw[r]:5.1f}{status}")

    # Schedule chain
    print(f"\n  Schedule chain δW[16..24]:")
    for r in range(16, 25):
        p0 = sum(1 for dw in dw_best if dw[r] == 0) / len(dw_best) * 100
        print(f"    δW[{r}]: {mean_dw[r]:.1f} бит, P(=0)={p0:.1f}%")

    print(f"\n  HW(δH): mean={np.mean(dh_best):.1f}, min={min(dh_best)}, std={np.std(dh_best):.1f}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("3. BREAK ROUND ОПТИМИЗАЦИЯ: r=3 vs r=4 vs r=5")
    print("=" * 70)

    for br in [2, 3, 4, 5]:
        r17_pass = 0
        r18_pass = 0
        r19_pass = 0
        dhs = []
        both_list = []

        for trial in range(2000):
            W1 = [np.random.randint(0, 2**32) for _ in range(16)]
            profile, dw, dh, bz = run_a_repair(W1, br, 0)  # bit 0

            dhs.append(dh)
            both_list.append(bz)

            # Count passes
            all_zero_to = 0
            for r in range(64):
                if profile[r] == 0:
                    all_zero_to = r + 1
                else:
                    break

            passed = [r for r in range(64) if profile[r] == 0]
            max_consec_after_break = 0
            consec = 0
            for r in range(64):
                if profile[r] == 0:
                    consec += 1
                    max_consec_after_break = max(max_consec_after_break, consec)
                else:
                    consec = 0

            if len(passed) > 0 and max(passed) >= 17 and profile[17] == 0:
                r17_pass += 1
            if len(passed) > 0 and max(passed) >= 18 and profile[17] == 0 and profile[18] == 0:
                r18_pass += 1

        print(f"  break r={br}, bit=0:")
        print(f"    BOTH=0: mean={np.mean(both_list):.1f}")
        print(f"    r=17 pass: {r17_pass/20:.1f}%")
        print(f"    r=18 pass: {r18_pass/20:.1f}%")
        print(f"    HW(δH): mean={np.mean(dhs):.1f}, min={min(dhs)}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("4. КОМБИНАЦИЯ: лучший break + W1 search")
    print("=" * 70)

    # Ищем W1 где a-repair + best break даёт максимальный BOTH=0
    best_both = 0
    best_dh = 256
    best_W1 = None

    for trial in range(10000):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        profile, dw, dh, bz = run_a_repair(W1, break_round, best_bit)

        if bz > best_both or (bz == best_both and dh < best_dh):
            best_both = bz
            best_dh = dh
            best_W1 = list(W1)

            # Показываем прогресс
            if bz > 10:
                passed = [r + 1 for r in range(64) if profile[r] == 0]
                print(f"    trial {trial}: BOTH={bz}, δH={dh}, zeros at: {passed[:20]}...")

    print(f"\n  Лучший результат из 10K:")
    print(f"    BOTH=0 раундов: {best_both}")
    print(f"    HW(δH): {best_dh}")


if __name__ == "__main__":
    main()
