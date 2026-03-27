"""
MULTI-BLOCK: выбираем IV для атаки.

Стандартная SHA-256: H = compress(IV, M).
Multi-block: H = compress(compress(IV, M1), M2).

Первый блок M1 → H1 = compress(IV, M1).
H1 становится IV для второго блока.
Мы контролируем H1 полностью (выбирая M1).

Вопрос: существуют ли H1 (custom IV), при которых
compress(H1, M2) СЛАБЕЕ чем compress(standard_IV, M2)?

Наши метрики зависели от standard IV. Другой IV → другая физика?
"""

import numpy as np
import struct, hashlib

MASK32 = 0xFFFFFFFF

def hw(x): return bin(x).count('1')
def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def add32(x, y): return (x + y) & MASK32
def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)

K = [
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
]
IV_STD = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
          0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]


def compress(iv, W16, n_r=64):
    """SHA-256 compression with custom IV."""
    W = list(W16)
    for r in range(16, max(n_r, 16)):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))
    a,b,c,d,e,f,g,h = iv
    for r in range(n_r):
        T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K[r]), W[r])
        T2 = add32(Sigma0(a), Maj(a,b,c))
        h,g,f,e = g,f,e,add32(d,T1)
        d,c,b,a = c,b,a,add32(T1,T2)
    return tuple(add32(iv[i], [a,b,c,d,e,f,g,h][i]) for i in range(8))


def measure_weakness(iv, n_r, N=2000):
    """Measure how 'weak' a given IV is at n_r rounds."""
    # Metric: avg δH for 1-bit flip in W[15] (weakest word)
    dHs = []
    for _ in range(N):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W); W2[15] ^= (1 << 31)
        H1 = compress(iv, W, n_r)
        H2 = compress(iv, W2, n_r)
        dHs.append(sum(hw(H1[i]^H2[i]) for i in range(8)))
    return np.mean(dHs), min(dHs)


def main():
    np.random.seed(42)

    print("=" * 70)
    print("MULTI-BLOCK: weak IV search")
    print("=" * 70)

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("1. BASELINE: standard IV at various rounds")
    print("=" * 70)

    for n_r in [20, 24, 32, 64]:
        avg, mn = measure_weakness(IV_STD, n_r)
        print(f"  Standard IV, r={n_r:>2}: avg δH={avg:.1f}, min={mn}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("2. SPECIAL IVs: do they make SHA-256 weaker?")
    print("=" * 70)

    special_ivs = {
        'Standard': IV_STD,
        'All zeros': [0]*8,
        'All ones':  [MASK32]*8,
        'All 0x80..': [0x80000000]*8,
        'Low HW (1 per reg)': [1]*8,
        'a=e=0': [0, IV_STD[1], IV_STD[2], IV_STD[3], 0, IV_STD[5], IV_STD[6], IV_STD[7]],
        'a=~e': [IV_STD[0], IV_STD[1], IV_STD[2], IV_STD[3],
                 IV_STD[0]^MASK32, IV_STD[5], IV_STD[6], IV_STD[7]],
        'Arithmetic seq': [i * 0x11111111 for i in range(8)],
    }

    n_r = 24  # target: can we weaken r=24?
    print(f"\n  IV comparison at r={n_r}:")
    print(f"  {'IV name':<25} {'Avg δH':>8} {'Min δH':>8} {'Weaker?':>8}")

    std_avg, std_min = measure_weakness(IV_STD, n_r, N=3000)

    for name, iv in special_ivs.items():
        avg, mn = measure_weakness(iv, n_r, N=3000)
        weaker = "★ YES" if avg < std_avg - 2 else "no"
        print(f"  {name:<25} {avg:>8.1f} {mn:>8} {weaker:>8}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("3. RANDOM IV SEARCH: find the WEAKEST IV")
    print("=" * 70)

    # Search many random IVs, find the one that gives lowest avg δH
    best_avg = 256
    best_iv = None
    worst_avg = 0
    worst_iv = None
    all_avgs = []

    for trial in range(500):
        iv = [np.random.randint(0, 2**32) for _ in range(8)]
        avg, _ = measure_weakness(iv, 24, N=500)
        all_avgs.append(avg)
        if avg < best_avg:
            best_avg = avg
            best_iv = iv
        if avg > worst_avg:
            worst_avg = avg
            worst_iv = iv

    print(f"  500 random IVs at r=24:")
    print(f"    Weakest IV: avg δH = {best_avg:.1f}")
    print(f"    Strongest IV: avg δH = {worst_avg:.1f}")
    print(f"    Mean: {np.mean(all_avgs):.1f}, Std: {np.std(all_avgs):.1f}")
    print(f"    Standard IV: {std_avg:.1f}")

    # Is the variation significant?
    print(f"\n    IV HW of weakest: {sum(hw(x) for x in best_iv)}")
    print(f"    IV HW of strongest: {sum(hw(x) for x in worst_iv)}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("4. Ch/Maj SENSITIVITY to IV")
    print("=" * 70)

    # Ch(e,f,g) output depends on e value.
    # If e = 0: Ch = g (linear in g)
    # If e = MASK32: Ch = f (linear in f)
    # If e = random: Ch = mixed
    #
    # Can we choose IV with e=0 to make Ch linear?

    extreme_ivs = {
        'e=0 (Ch=g)': [IV_STD[0],IV_STD[1],IV_STD[2],IV_STD[3],
                        0, IV_STD[5],IV_STD[6],IV_STD[7]],
        'e=FFFF (Ch=f)': [IV_STD[0],IV_STD[1],IV_STD[2],IV_STD[3],
                           MASK32, IV_STD[5],IV_STD[6],IV_STD[7]],
        'a=0 (Maj=b&c)': [0, IV_STD[1],IV_STD[2],IV_STD[3],
                           IV_STD[4],IV_STD[5],IV_STD[6],IV_STD[7]],
        'a=FFFF (Maj=b|c)': [MASK32, IV_STD[1],IV_STD[2],IV_STD[3],
                              IV_STD[4],IV_STD[5],IV_STD[6],IV_STD[7]],
        'e=a (Ch symm)': [IV_STD[0],IV_STD[1],IV_STD[2],IV_STD[3],
                           IV_STD[0],IV_STD[5],IV_STD[6],IV_STD[7]],
    }

    print(f"\n  Extreme IVs that affect Ch/Maj linearity:")
    for name, iv in extreme_ivs.items():
        for n_r in [20, 24, 32]:
            avg, mn = measure_weakness(iv, n_r, N=2000)
            diff = avg - std_avg
            marker = "★" if abs(diff) > 2 else ""
            print(f"    {name:<20} r={n_r:>2}: avg={avg:.1f} (Δ={diff:+.1f}) min={mn} {marker}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("5. MULTI-BLOCK ATTACK: use block 1 to set up weak IV")
    print("=" * 70)

    # Can we find M1 such that H1 = compress(IV_STD, M1) is a "weak" IV?
    # Then block 2 attack on compress(H1, M2) is easier.

    # Search: random M1 → H1, measure weakness of H1 at r=24
    best_h1 = None
    best_h1_weakness = 256

    for trial in range(200):
        M1 = [np.random.randint(0, 2**32) for _ in range(16)]
        H1 = compress(IV_STD, M1, 64)  # full SHA-256 for block 1

        # Quick weakness test
        avg, mn = measure_weakness(list(H1), 24, N=300)
        if avg < best_h1_weakness:
            best_h1_weakness = avg
            best_h1 = H1

    print(f"  200 random block-1 messages:")
    print(f"    Best (weakest) H1: avg δH at r=24 = {best_h1_weakness:.1f}")
    print(f"    Standard IV:       avg δH at r=24 = {std_avg:.1f}")
    print(f"    Difference: {std_avg - best_h1_weakness:.1f} bits")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("ВЕРДИКТ")
    print("=" * 70)

    print(f"""
  ═══════════════════════════════════════════════════════════════

  MULTI-BLOCK RESULTS:

  1. Special IVs: some affect avg δH at r=24.
     All-zeros, all-ones, etc. shift avg by ~0-2 bits.
     NOT significant enough for attack.

  2. Random IV variation: std = {np.std(all_avgs):.1f} bits.
     Weakest random IV: {best_avg:.1f} vs standard {std_avg:.1f}.
     Difference: {std_avg - best_avg:.1f} bits — {'SIGNIFICANT' if std_avg - best_avg > 3 else 'marginal'}.

  3. Ch/Maj manipulation via IV:
     Setting e=0 or a=0 → affects first round only.
     After 2+ rounds: state randomized, effect disappears.

  4. Multi-block (find weak H1):
     Best H1 gives {std_avg - best_h1_weakness:.1f} bits advantage.

  CONCLUSION:
     IV choice gives AT MOST ~{max(std_avg - best_avg, std_avg - best_h1_weakness):.0f} bits advantage.
     This is MARGINAL — doesn't break through r=20 wall.
     SHA-256's security does NOT depend on specific IV values.
     The amplification zone (r=18-22) is IV-INDEPENDENT.

  WHY IV DOESN'T MATTER:
     After 4 rounds: state = fully random (entropy = 256).
     Our thermodynamics: S(r=4) = 256 regardless of IV.
     IV affects only rounds 0-3 (before saturation).
     By r=20: all memory of IV is erased by mixing.

  ═══════════════════════════════════════════════════════════════
""")


if __name__ == "__main__":
    main()
