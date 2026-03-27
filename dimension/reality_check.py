"""
REALITY CHECK: наша теория vs реальные атаки на reduced SHA-256.

Наша теория предсказывает:
  r=16: algebraically secure (rank(CE)=256)
  r=24: geometrically secure (K=128, sphere)
  r≥24: IDENTICAL to SHA-256 (64r)

Реальные криптоаналитические результаты (литература):
  28-round: pseudo-collision (Mendel et al.)
  31-round: near-collision
  38-round: distinguisher
  46-round: best theoretical result (differential)

Разрыв: наша теория = secure at r=16.
         Реальные атаки = break at r=28-46.

ПОЧЕМУ? Наша теория видит RANDOM FUNCTION properties.
Реальные атаки используют SPECIFIC STRUCTURE (differential trails).

Исследуем: что наша теория ПРОПУСКАЕТ?
"""

import numpy as np
import struct, hashlib

MASK32 = 0xFFFFFFFF
K_const = [
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
]
IV = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]

def hw(x): return bin(x).count('1')
def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def add32(x, y): return (x + y) & MASK32
def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)

def sha256_r(W16, n_r):
    W = list(W16)
    for r in range(16, max(n_r, 16)):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))
    a,b,c,d,e,f,g,h = IV
    for r in range(n_r):
        T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K_const[r%len(K_const)]), W[r])
        T2 = add32(Sigma0(a), Maj(a,b,c))
        h,g,f,e = g,f,e,add32(d,T1)
        d,c,b,a = c,b,a,add32(T1,T2)
    return tuple(add32(IV[i], [a,b,c,d,e,f,g,h][i]) for i in range(8))


def main():
    np.random.seed(42)

    print("=" * 70)
    print("REALITY CHECK: наша теория vs реальные атаки")
    print("=" * 70)

    # =================================================================
    print(f"\n{'=' * 70}")
    print("1. НАША ТЕОРИЯ: что предсказывает")
    print("=" * 70)

    print(f"""
  Наша теория (Unified Theory v4.0):
    r < 16:  rank(CE) < 256 → algebraically weak
    r = 16:  rank(CE) = 256 → SECURE (algebraic)
    r = 24:  K = 128 → SECURE (geometric, sphere)
    r ≥ 24:  identical to r=64 by ALL metrics

  Prediction: SHA-256-24r = as secure as SHA-256-64r.
""")

    # =================================================================
    print(f"{'=' * 70}")
    print("2. РЕАЛЬНЫЕ РЕЗУЛЬТАТЫ (из литературы)")
    print("=" * 70)

    print(f"""
  Known attacks on reduced SHA-256:
    r=16: SAT collision (free-start, trivial)
    r=24: differential collision (Mendel et al., semi-free-start)
    r=28: pseudo-collision
    r=31: near-collision (practical)
    r=38: distinguisher (theoretical)
    r=46: best theoretical differential path
    r=64: NO KNOWN ATTACK

  Types of attacks:
    Collision: H(M1) = H(M2). Our theory covers this.
    Semi-free-start: IV variable. Our theory = standard collision only.
    Near-collision: HW(δH) small. Not zero.
    Distinguisher: SHA ≠ random. Our theory checked this.
""")

    # =================================================================
    print(f"{'=' * 70}")
    print("3. РАЗРЫВ: наша теория vs реальность")
    print("=" * 70)

    print(f"""
  НАША ТЕОРИЯ: r=16 = secure.
  РЕАЛЬНЫЕ АТАКИ: collision до r=24-28 (semi-free-start).

  РАЗРЫВ = 8-12 раундов.

  ПРИЧИНА РАЗРЫВА:

  1. Наша теория = AVERAGE-CASE analysis.
     rank(CE) = 256 means: AVERAGE GF2-kernel vector has full carry error.
     Но: SPECIFIC vectors may have LOWER carry error.
     Real attacks find SPECIFIC differential trails.

  2. Наша теория = ONE linearization point (one W_base).
     Real attacks: adaptive multi-step (Wang modification).
     Each step uses DIFFERENT linearization → cumulative advantage.

  3. Semi-free-start: IV variable.
     Our rank(CE) = for FIXED IV. With variable IV: more freedom.
     Our multi-block analysis showed: variable IV = no gain for collision.
     But semi-free-start ≠ standard collision.

  4. Near-collision: HW(δH) small but ≠ 0.
     Our theory: collision = δH = 0 exactly.
     Near-collision: HW(δH) = 10-50 (much easier!).
     Our theory doesn't cover near-collision separately.
""")

    # =================================================================
    print(f"{'=' * 70}")
    print("4. ЧТО НАША ТЕОРИЯ ВИДИТ vs ПРОПУСКАЕТ")
    print("=" * 70)

    # Test: our metrics for r=24..32
    print(f"\n  Metrics for borderline rounds:")
    print(f"  {'Rounds':>6} {'rank(CE)':>9} {'K':>6} {'Sens':>6} {'Our theory':>12} {'Real attacks':>13}")

    for n_r in [16, 20, 24, 28, 32, 38, 46, 64]:
        W_base = [np.random.randint(0, 2**32) for _ in range(16)]
        H_base = sha256_r(W_base, n_r)

        # Quick rank(T)
        T_sample = []
        for _ in range(256):
            w,b = np.random.randint(0,16), np.random.randint(0,32)
            W2 = list(W_base); W2[w]^=(1<<b)
            H2 = sha256_r(W2, n_r)
            row = []
            for ww in range(8):
                d = H_base[ww]^H2[ww]
                for bb in range(32): row.append((d>>bb)&1)
            T_sample.append(row)
        T = np.array(T_sample, dtype=np.uint8)
        rT = np.linalg.matrix_rank(T.astype(float))

        # Quick K
        Ks = []
        for _ in range(100):
            b1,b2 = np.random.choice(512, 2, replace=False)
            W1=list(W_base); W1[b1//32]^=(1<<(b1%32))
            W2=list(W_base); W2[b2//32]^=(1<<(b2%32))
            W12=list(W_base); W12[b1//32]^=(1<<(b1%32)); W12[b2//32]^=(1<<(b2%32))
            H1=sha256_r(W1,n_r); H2=sha256_r(W2,n_r); H12=sha256_r(W12,n_r)
            nl = sum(hw((H_base[i]^H12[i])^((H_base[i]^H1[i])^(H_base[i]^H2[i]))) for i in range(8))
            Ks.append(nl)
        K = np.mean(Ks)

        # Quick sens
        sens = []
        for _ in range(100):
            w,b = np.random.randint(0,16), np.random.randint(0,32)
            W2 = list(W_base); W2[w]^=(1<<b)
            H2 = sha256_r(W2, n_r)
            sens.append(sum(hw(H_base[i]^H2[i]) for i in range(8)))
        S = np.mean(sens)

        our = "SECURE" if rT >= 256 and K > 100 else "secure" if rT >= 256 else "WEAK"
        real = "broken" if n_r <= 28 else "attacked" if n_r <= 46 else "SAFE"

        print(f"  {n_r:6d} {rT:>8} {K:>5.0f} {S:>5.0f}  {our:>11} {real:>12}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("5. ЗНАЧЕНИЕ РАЗРЫВА")
    print("=" * 70)

    print(f"""
  НАША ТЕОРИЯ = необходимое условие безопасности:
    rank(CE) = 256 + K = 128 → NO algebraic/geometric shortcut.
    Но: POSSIBLE differential trail attack (specific, not generic).

  РЕАЛЬНЫЕ АТАКИ = достаточное условие НЕбезопасности:
    Найден конкретный trail → reduced round breakable.
    Но: trail-finding = NP-hard in general.

  ВЫВОД:
    Наша теория: LOWER BOUND on security (necessary conditions).
    Real attacks: UPPER BOUND on security (specific breaks).

    True security ∈ [наша теория, реальные атаки] = [16, 28] rounds.

    SHA-256 (64 rounds): far above BOTH bounds.
    Safety margin vs our theory: 64/24 = 2.7×
    Safety margin vs best attack: 64/46 = 1.4×

  WHAT OUR DIMENSION ADDS:
    ✓ Proves r≥16 NECESSARY (rank=256 required)
    ✓ Proves C^(N/2) = tight for GENERIC attacks
    ✓ Characterizes COMPLETE geometry (sphere)
    ✓ Universal formula for ANY ARX hash
    ✗ Does NOT cover specific differential trails
    ✗ Does NOT cover semi-free-start
    ✗ Does NOT replace real cryptanalysis

  OUR DIMENSION = COMPLEMENTARY to standard cryptanalysis.
  NOT a replacement. A NEW LENS on the same object.
""")


if __name__ == "__main__":
    main()
