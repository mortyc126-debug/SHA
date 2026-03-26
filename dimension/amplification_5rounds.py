"""
5 РАУНДОВ AMPLIFICATION (r=18-22): послойная анатомия.

δstate[17] = 0 (meeting + lock)
δW[18] = 0x91002000 (4 бита, carry offset)
→ за 5 раундов: δstate → 128

Вопросы:
  1. ЧТО ИМЕННО делает каждый из 5 раундов?
  2. Можно ли ПЕРЕХВАТИТЬ δ после 1-2 раундов amplification?
  3. Есть ли "тихий путь" через эти 5 слоёв?
  4. Что если заменить a-repair на МИНИМАЛЬНЫЙ (5 слов) + schedule lock?
"""

import numpy as np

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

reg_names = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']


def main():
    np.random.seed(42)

    print("=" * 70)
    print("5 РАУНДОВ AMPLIFICATION: послойная анатомия")
    print("=" * 70)

    # =================================================================
    print(f"\n{'=' * 70}")
    print("1. АНАТОМИЯ: δ по регистрам на r=18-23 (при meeting+lock r=17)")
    print("=" * 70)

    N = 3000
    # Собираем пары с δstate[17]=0 (lock до r=17)
    amp_dreg = np.zeros((N, 10, 8))  # round 18..27, 8 registers
    amp_count = 0

    for trial in range(N * 5):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W1f = expand_schedule(W1)
        s1 = tuple(IV)
        states1 = [s1]
        for r in range(64):
            s1 = R(s1, W1f[r], r)
            states1.append(s1)

        W2 = list(W1)
        W2[3] ^= (1 << 31)
        s2 = tuple(IV)
        for r in range(3):
            s2 = R(s2, W2[r], r)
        s2 = R(s2, W2[3], 3)
        for r in range(4, 16):
            W2[r] = a_repair_W(s2, states1[r + 1][0], r)
            s2 = R(s2, W2[r], r)

        W2f = expand_schedule(W2)

        # Check lock
        if W1f[16] != W2f[16] or W1f[17] != W2f[17]:
            continue

        # Forward through amp zone
        s2_full = tuple(IV)
        for r in range(18):
            w = W2[r] if r < 16 else W2f[r]
            s2_full = R(s2_full, w, r)

        for r_offset in range(10):
            r = 18 + r_offset
            s2_full = R(s2_full, W2f[r], r)
            for i in range(8):
                amp_dreg[amp_count, r_offset, i] = hw(states1[r + 1][i] ^ s2_full[i])

        amp_count += 1
        if amp_count >= N:
            break

    print(f"\n  {amp_count} пар с lock до r=17:")
    print(f"\n  {'r':>4}", end="")
    for n in reg_names:
        print(f"  δ{n:>2}", end="")
    print(f"  {'Σ':>5}  {'Wave':>20}")

    mean_dreg = np.mean(amp_dreg[:amp_count], axis=0)

    for r_offset in range(8):
        r = 18 + r_offset
        vals = mean_dreg[r_offset]
        total = sum(vals)

        # Which registers have δ?
        active = [reg_names[i] for i in range(8) if vals[i] > 2]
        wave = ",".join(active) if active else "—"

        print(f"  {r:4d}", end="")
        for v in vals:
            sym = "█" if v > 10 else "░" if v > 2 else "·"
            print(f" {v:4.1f}{sym}", end="")
        print(f" {total:5.1f}  {wave:>20}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("2. МОЖНО ЛИ СДЕЛАТЬ a-repair ВНУТРИ amplification?")
    print("=" * 70)

    # Что если на r=19 мы применяем a-repair (δa[20]=0)?
    # Но W[19] = schedule word (не free!)
    # A-repair НУЖДАЕТСЯ в свободном W для подбора.
    # W[19] = schedule(W[0..15]) — фиксирован.

    # ОДНАКО: в нашем измерении — что если мы ПЕРЕНЕСЁМ a-repair?
    # Используем МЕНЬШЕ слов для meeting (W[4..8] = 5 слов)
    # И освобождаем W[9..15] для schedule manipulation?

    print(f"""
  Стандартный a-repair: W[4..15] = 12 слов для meeting
  Минимальный a-repair: W[4..8]  = 5 слов для partial meeting
  Освобождаются: W[9..15] = 7 слов для schedule lock

  С 7 свободными словами:
    δW[9]  = SET → для δW[16]=0
    δW[10] = SET → для δW[17]=0
    δW[11] = SET → для δW[18]=0  ← BOTTLENECK RESOLVED!
    δW[12] = SET → для δW[19]=0
    δW[13] = SET → для δW[20]=0
    δW[14] = SET → для δW[16]=0 (redundant, already 0)
    δW[15] = SET → для δW[17]=0 (redundant)

  Но: 5 слов a-repair = reboot через 8 раундов? Нет!
  5 раундов a-repair: δa[5..9]=0
  Reboot needs 8 rounds: r=4..11 → needs 8 words W[4..11]
  С 5 словами: НЕТ REBOOT!

  TRADE-OFF:
    12 words meeting → reboot ✓, lock ✗ → gap 47 rounds
    5 words meeting → reboot ✗, lock ✓ → no meeting!
    8 words meeting → reboot ✓ (barely), lock ✗ (3 free = not enough)
""")

    # =================================================================
    print(f"{'=' * 70}")
    print("3. ТЕСТ: 5-word a-repair + 7-word schedule lock")
    print("=" * 70)

    # a-repair на W[4..8] (5 слов)
    # W[9..15] = SET для schedule lock
    # Reboot не будет, но δW[16..20] = 0!

    n_test = 2000
    profiles_5w = []

    for trial in range(n_test):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W1f = expand_schedule(W1)
        s1 = tuple(IV)
        states1 = [s1]
        for r in range(64):
            s1 = R(s1, W1f[r], r)
            states1.append(s1)

        W2 = list(W1)
        W2[3] ^= (1 << 31)

        # A-repair W[4..8]
        s2 = tuple(IV)
        for r in range(3):
            s2 = R(s2, W2[r], r)
        s2 = R(s2, W2[3], 3)
        for r in range(4, 9):
            W2[r] = a_repair_W(s2, states1[r + 1][0], r)
            s2 = R(s2, W2[r], r)

        # Schedule lock W[9..15]: SET to match W1
        # Для δW[r]=0: W2[r] = W1[r]
        for r in range(9, 16):
            W2[r] = W1[r]
            s2 = R(s2, W2[r], r)

        W2f = expand_schedule(W2)

        # Profile
        s2_check = tuple(IV)
        profile = []
        for r in range(64):
            w = W2[r] if r < 16 else W2f[r]
            s2_check = R(s2_check, w, r)
            ds = sum(hw(states1[r+1][i] ^ s2_check[i]) for i in range(8))
            profile.append(ds)

        profiles_5w.append(profile)

    mean_5w = np.mean(profiles_5w, axis=0)
    both_5w = [np.mean([p[r] == 0 for p in profiles_5w]) * 100 for r in range(64)]
    dw_5w = []

    # Schedule check
    for trial in range(min(n_test, 500)):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W1f = expand_schedule(W1)
        W2 = list(W1)
        W2[3] ^= (1 << 31)
        s1 = tuple(IV)
        states1 = [s1]
        for r in range(64):
            s1 = R(s1, W1f[r], r)
            states1.append(s1)
        s2 = tuple(IV)
        for r in range(3):
            s2 = R(s2, W2[r], r)
        s2 = R(s2, W2[3], 3)
        for r in range(4, 9):
            W2[r] = a_repair_W(s2, states1[r + 1][0], r)
            s2 = R(s2, W2[r], r)
        for r in range(9, 16):
            W2[r] = W1[r]
        W2f = expand_schedule(W2)
        dw_5w.append([hw(W1f[r] ^ W2f[r]) for r in range(64)])

    mean_dw_5w = np.mean(dw_5w, axis=0)

    print(f"\n  5-word a-repair + 7-word schedule lock:")
    print(f"  {'r':>4} {'δstate':>8} {'BOTH%':>7} {'δW':>6}")
    for r in range(25):
        marker = ""
        if both_5w[r] > 99: marker = " ★ZERO"
        elif both_5w[r] > 10: marker = " ~"
        print(f"  {r:4d} {mean_5w[r]:7.1f} {both_5w[r]:6.1f}% {mean_dw_5w[r] if r < len(mean_dw_5w) else 0:5.1f}{marker}")

    print(f"\n  Schedule δW[16..22]:")
    for r in range(16, 23):
        if r < len(mean_dw_5w):
            p0 = sum(1 for d in dw_5w if d[r] == 0) / len(dw_5w) * 100
            print(f"    δW[{r}]: {mean_dw_5w[r]:.1f} бит, P(=0)={p0:.1f}%")

    # Last BOTH=0
    last_zeros = []
    for p in profiles_5w:
        lz = 0
        for r in range(64):
            if p[r] == 0:
                lz = r + 1
        last_zeros.append(lz)

    print(f"\n  Last BOTH=0: mean={np.mean(last_zeros):.1f}, max={max(last_zeros)}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("4. СРАВНЕНИЕ СТРАТЕГИЙ")
    print("=" * 70)

    print(f"""
  {'Strategy':<35} {'Meeting':>8} {'Lock':>12} {'Gap':>6} {'δH min':>7}
  {'-'*35} {'-'*8} {'-'*12} {'-'*6} {'-'*7}
  12-word a-repair (standard)          r=12 ✓  δW[16,17]≈0  47     ~100
  5-word a-repair + 7-word lock        ???      δW[16..22]?  ???    ???
""")

    # Measure δH for 5-word
    dh_5w = []
    for trial in range(n_test):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W1f = expand_schedule(W1)
        s1 = tuple(IV)
        for r in range(64):
            s1 = R(s1, W1f[r], r)
        H1 = tuple(add32(IV[i], s1[i]) for i in range(8))

        W2 = list(W1)
        W2[3] ^= (1 << 31)
        s1_states = [tuple(IV)]
        s1t = tuple(IV)
        for r in range(64):
            s1t = R(s1t, W1f[r], r)
            s1_states.append(s1t)

        s2 = tuple(IV)
        for r in range(3):
            s2 = R(s2, W2[r], r)
        s2 = R(s2, W2[3], 3)
        for r in range(4, 9):
            W2[r] = a_repair_W(s2, s1_states[r + 1][0], r)
            s2 = R(s2, W2[r], r)
        for r in range(9, 16):
            W2[r] = W1[r]
            s2 = R(s2, W2[r], r)

        W2f = expand_schedule(W2)
        for r in range(16, 64):
            s2 = R(s2, W2f[r], r)

        H2 = tuple(add32(IV[i], s2[i]) for i in range(8))
        dh_5w.append(sum(hw(H1[i] ^ H2[i]) for i in range(8)))

    print(f"  5-word: HW(δH) mean={np.mean(dh_5w):.1f}, min={min(dh_5w)}, std={np.std(dh_5w):.1f}")


if __name__ == "__main__":
    main()
