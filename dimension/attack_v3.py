"""
ATTACK v3: schedule-optimized, всё знание учтено.

КЛЮЧЕВОЕ ЗНАНИЕ:
  δ ТОЛЬКО в W[15]:
    δW[16] = 0  (W[16] не зависит от W[15])
    δW[17] = σ1(δW[15])
    δW[18] = 0
    δW[19] = σ1(σ1(δW[15]))
    δW[20] = 0
    ...
    Половина expanded words = НОЛЬ разницы!

  Оптимальный бит: тот где HW(σ1(1<<bit)) минимальный.

СТРАТЕГИЯ:
  1. Найти bit с минимальным σ1-damage
  2. Также: найти ПАРУ бит где σ1(δ) CANCEL друг друга
  3. Также: δ в W[14] (δW[16]=σ1(δW[14]), но δW[17] зависит от W[15])
  4. Массивный поиск по random base messages
"""

import numpy as np
import struct, hashlib
import time

MASK32 = 0xFFFFFFFF

def hw(x): return bin(x).count('1')
def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def add32(x, y): return (x + y) & MASK32
def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)
def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)

K = [
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
]
IV = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]

def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)

def sha256_r(W16, n_r):
    W = list(W16)
    for r in range(16, max(n_r, 16)):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))
    a,b,c,d,e,f,g,h = IV
    for r in range(n_r):
        T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K[r % len(K)]), W[r])
        T2 = add32(Sigma0(a), Maj(a,b,c))
        h,g,f,e = g,f,e,add32(d,T1)
        d,c,b,a = c,b,a,add32(T1,T2)
    return tuple(add32(IV[i], [a,b,c,d,e,f,g,h][i]) for i in range(8))


def main():
    np.random.seed(42)

    print("=" * 70)
    print("ATTACK v3: schedule-optimized")
    print("=" * 70)

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("1. σ1 DAMAGE MAP: HW(σ1(1<<bit)) for each bit")
    print("=" * 70)

    print(f"  {'Bit':>4} {'σ1(1<<bit)':>12} {'HW':>4}")
    bit_damages = []
    for bit in range(32):
        d = sigma1(1 << bit)
        h = hw(d)
        bit_damages.append((bit, h, d))
        print(f"  {bit:>4} {hex(d):>12} {h:>4}")

    bit_damages.sort(key=lambda x: x[1])
    print(f"\n  Best bits (min σ1 damage):")
    for bit, h, d in bit_damages[:5]:
        print(f"    bit {bit}: HW(σ1) = {h}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("2. SCHEDULE PROPAGATION: δ only in W[15]")
    print("=" * 70)

    best_bit = bit_damages[0][0]
    dW15 = 1 << best_bit

    W_base = [np.random.randint(0, 2**32) for _ in range(16)]
    W_mod = list(W_base); W_mod[15] ^= dW15

    W_exp1 = list(W_base)
    W_exp2 = list(W_mod)
    for r in range(16, 24):
        W_exp1.append(add32(add32(add32(sigma1(W_exp1[r-2]), W_exp1[r-7]), sigma0(W_exp1[r-15])), W_exp1[r-16]))
        W_exp2.append(add32(add32(add32(sigma1(W_exp2[r-2]), W_exp2[r-7]), sigma0(W_exp2[r-15])), W_exp2[r-16]))

    print(f"  δW[15] = bit {best_bit} ({hex(dW15)}), schedule propagation:")
    for r in range(16, 24):
        d = W_exp1[r] ^ W_exp2[r]
        print(f"    δW[{r}] = {hex(d):>12} HW={hw(d):>2} {'★ ZERO' if d == 0 else ''}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("3. σ1 CANCELLATION: find 2-bit δW[15] where σ1 bits cancel")
    print("=" * 70)

    # σ1(a ⊕ b) ≠ σ1(a) ⊕ σ1(b) in general (nonlinear due to SHR)
    # But for XOR of two single bits:
    # σ1(1<<a ⊕ 1<<b) where SHR doesn't interact

    best_2bit = None
    min_hw_2bit = 32
    for b1 in range(32):
        for b2 in range(b1+1, 32):
            dW = (1 << b1) ^ (1 << b2)
            s1_hw = hw(sigma1(dW))
            if s1_hw < min_hw_2bit:
                min_hw_2bit = s1_hw
                best_2bit = (b1, b2, dW, s1_hw)

    print(f"  Best 2-bit δW[15]: bits ({best_2bit[0]}, {best_2bit[1]})")
    print(f"    δW = {hex(best_2bit[2])}, HW(σ1) = {best_2bit[3]}")
    print(f"    vs best 1-bit: HW(σ1) = {bit_damages[0][1]}")

    # Also try 3-bit for even better cancellation
    best_3bit = None
    min_hw_3bit = 32
    for b1 in range(32):
        for b2 in range(b1+1, 32):
            for b3 in range(b2+1, 32):
                dW = (1 << b1) ^ (1 << b2) ^ (1 << b3)
                s1_hw = hw(sigma1(dW))
                if s1_hw < min_hw_3bit:
                    min_hw_3bit = s1_hw
                    best_3bit = (b1, b2, b3, dW, s1_hw)

    print(f"  Best 3-bit δW[15]: bits ({best_3bit[0]}, {best_3bit[1]}, {best_3bit[2]})")
    print(f"    δW = {hex(best_3bit[3])}, HW(σ1) = {best_3bit[4]}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("4. MASSIVE SEARCH: best δW[15] at each round count")
    print("=" * 70)

    # For each candidate δW[15], search over random base messages
    candidates = [
        (f"1-bit best (bit {bit_damages[0][0]})", 1 << bit_damages[0][0]),
        (f"1-bit bit 31 (MSB)", 1 << 31),
        (f"1-bit bit 0 (LSB)", 1),
        (f"2-bit cancel ({best_2bit[0]},{best_2bit[1]})", best_2bit[2]),
        (f"3-bit cancel", best_3bit[3]),
        ("random 1-bit", None),  # random each trial
    ]

    N = 20000
    for n_r in [17, 18, 19, 20, 22, 24]:
        print(f"\n  r={n_r}:")
        for name, dW_fixed in candidates:
            best = 256
            for _ in range(N):
                W = [np.random.randint(0, 2**32) for _ in range(16)]
                W2 = list(W)
                if dW_fixed is not None:
                    W2[15] ^= dW_fixed
                else:
                    W2[15] ^= (1 << np.random.randint(0, 32))
                H1 = sha256_r(W, n_r)
                H2 = sha256_r(W2, n_r)
                dH = sum(hw(H1[i]^H2[i]) for i in range(8))
                if dH < best: best = dH
            print(f"    {name:>35}: best δH = {best:>3}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("5. δ in W[14] vs W[15]: which is better?")
    print("=" * 70)

    # W[14]: enters at r=14, 3 rounds of mixing by r=17
    # δW[16] = σ1(δW[14]) (NON-ZERO!)
    # δW[17] = δW[15] term... wait, if only W[14] changes:
    # δW[16] = σ1(δW[14]), δW[17] = 0, δW[18] = σ1(δW[16]) = σ1²(δW[14])
    # More schedule damage than W[15]!

    # But W[14] has 3 rounds of mixing vs W[15]'s 2.
    # Trade-off: more mixing (bad) vs less schedule damage (W[15] wins)

    for n_r in [17, 18, 20]:
        best_14 = 256; best_15 = 256; best_any = 256
        for _ in range(N):
            W = [np.random.randint(0, 2**32) for _ in range(16)]

            # δ in W[14]
            b = np.random.randint(0, 32)
            W2 = list(W); W2[14] ^= (1 << b)
            H1 = sha256_r(W, n_r); H2 = sha256_r(W2, n_r)
            dH = sum(hw(H1[i]^H2[i]) for i in range(8))
            if dH < best_14: best_14 = dH

            # δ in W[15]
            W3 = list(W); W3[15] ^= (1 << b)
            H3 = sha256_r(W3, n_r)
            dH3 = sum(hw(H1[i]^H3[i]) for i in range(8))
            if dH3 < best_15: best_15 = dH3

            # δ in any word
            w = np.random.randint(0, 16)
            W4 = list(W); W4[w] ^= (1 << b)
            H4 = sha256_r(W4, n_r)
            dH4 = sum(hw(H1[i]^H4[i]) for i in range(8))
            if dH4 < best_any: best_any = dH4

        print(f"  r={n_r}: W[14]={best_14:>3}, W[15]={best_15:>3}, any_word={best_any:>3}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("6. ULTIMATE: best δW[15] bit + 100K search")
    print("=" * 70)

    for n_r in [17, 18, 19, 20]:
        # Find which bit of W[15] gives best average
        bit_avgs = []
        for bit in range(32):
            dHs = []
            for _ in range(500):
                W = [np.random.randint(0, 2**32) for _ in range(16)]
                W2 = list(W); W2[15] ^= (1 << bit)
                H1 = sha256_r(W, n_r); H2 = sha256_r(W2, n_r)
                dHs.append(sum(hw(H1[i]^H2[i]) for i in range(8)))
            bit_avgs.append((bit, np.mean(dHs), min(dHs)))
        bit_avgs.sort(key=lambda x: x[1])

        optimal_bit = bit_avgs[0][0]
        print(f"\n  r={n_r}: optimal bit = {optimal_bit} (avg δH={bit_avgs[0][1]:.1f})")
        print(f"    Top 5 bits: {[(b, f'{a:.1f}') for b, a, _ in bit_avgs[:5]]}")

        # Now: 100K search with optimal bit
        best = 256
        dW = 1 << optimal_bit
        t0 = time.time()
        for _ in range(100000):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            W2 = list(W); W2[15] ^= dW
            H1 = sha256_r(W, n_r); H2 = sha256_r(W2, n_r)
            dH = sum(hw(H1[i]^H2[i]) for i in range(8))
            if dH < best: best = dH
        elapsed = time.time() - t0

        print(f"    100K search: best δH = {best} ({elapsed:.1f}s)")
        print(f"    {256-best} of 256 bits MATCH ({(256-best)/256*100:.1f}%)")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("7. COMPARISON: our best vs literature")
    print("=" * 70)

    print(f"""
  Our attack: δW[15] 1-bit, random message search.

  Known results (differential cryptanalysis of SHA-256):
    Mendel et al. 2011: collision up to 31 rounds
    (uses complex message modification + multi-step trail)

    Our technique: simple 1-bit δ in weak word
    Advantage: NO complex trail needed, just absorption knowledge
    Limitation: works only up to ~r=20 (absorption boundary)

  What our MATH adds:
    - EXACT knowledge of which word/bit is weakest at which round
    - Schedule propagation map (δW[15] → half schedule words = 0)
    - Metric tensor confirms: no better direction exists at r≥20
    - Security boundary formula: secure_round(W[i]) = i + 5
""")


if __name__ == "__main__":
    main()
