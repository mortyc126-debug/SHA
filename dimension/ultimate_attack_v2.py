"""
ULTIMATE ATTACK v2: fix the bug.

Bug in v1: a-repair overwrites W[1..15], so δW[13-15] gets erased → W2=W1.
Fix: δ ONLY in W[0] (not overwritten by repair).
     a-repair controls W[1..15] to minimize δ in rounds 1-15.
     δW[0] survives → actual difference propagates.

Also try: NO repair, just smart δW selection using metric knowledge.
"""

import numpy as np
import struct, hashlib
import time

MASK32 = 0xFFFFFFFF

def hw(x): return bin(x).count('1')
def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def add32(x, y): return (x + y) & MASK32
def sub32(x, y): return (x - y) & MASK32
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
IV = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]


def sha256_r(W16, n_r):
    W = list(W16)
    for r in range(16, max(n_r, 16)):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))
    a,b,c,d,e,f,g,h = IV
    for r in range(n_r):
        T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K[r]), W[r])
        T2 = add32(Sigma0(a), Maj(a,b,c))
        h,g,f,e = g,f,e,add32(d,T1)
        d,c,b,a = c,b,a,add32(T1,T2)
    return tuple(add32(IV[i], [a,b,c,d,e,f,g,h][i]) for i in range(8))


def a_repair_w0(W_base, dW0, n_r):
    """a-repair with δ ONLY in W[0]. W[1..15] modified by repair."""
    W2 = list(W_base)
    W2[0] ^= dW0

    a1,b1,c1,d1,e1,f1,g1,h1 = IV
    a2,b2,c2,d2,e2,f2,g2,h2 = IV

    # Round 0: both use their OWN W[0]
    T1=add32(add32(add32(add32(h1,Sigma1(e1)),Ch(e1,f1,g1)),K[0]),W_base[0])
    T2=add32(Sigma0(a1),Maj(a1,b1,c1))
    h1,g1,f1,e1=g1,f1,e1,add32(d1,T1); d1,c1,b1,a1=c1,b1,a1,add32(T1,T2)

    T1=add32(add32(add32(add32(h2,Sigma1(e2)),Ch(e2,f2,g2)),K[0]),W2[0])
    T2=add32(Sigma0(a2),Maj(a2,b2,c2))
    h2,g2,f2,e2=g2,f2,e2,add32(d2,T1); d2,c2,b2,a2=c2,b2,a2,add32(T1,T2)

    # Rounds 1-15: repair δa=0
    for r in range(1, 16):
        T1_1=add32(add32(add32(add32(h1,Sigma1(e1)),Ch(e1,f1,g1)),K[r]),W_base[r])
        T2_1=add32(Sigma0(a1),Maj(a1,b1,c1))
        a1_new=add32(T1_1,T2_1); e1_new=add32(d1,T1_1)

        T2_2=add32(Sigma0(a2),Maj(a2,b2,c2))
        T1_n=sub32(a1_new,T2_2)
        W2[r]=sub32(sub32(sub32(sub32(T1_n,h2),Sigma1(e2)),Ch(e2,f2,g2)),K[r])

        h1,g1,f1,e1=g1,f1,e1,e1_new; d1,c1,b1,a1=c1,b1,a1,a1_new
        h2,g2,f2,e2=g2,f2,e2,add32(d2,T1_n); d2,c2,b2,a2=c2,b2,a2,a1_new

    # VERIFY: W2 ≠ W_base (non-trivial)
    diff_words = sum(1 for i in range(16) if W_base[i] != W2[i])
    if diff_words == 0:
        return 0, W2, True  # trivial!

    H1 = sha256_r(W_base, n_r)
    H2 = sha256_r(W2, n_r)
    dH = sum(hw(H1[i]^H2[i]) for i in range(8))
    return dH, W2, False


def main():
    np.random.seed(42)

    print("=" * 70)
    print("ULTIMATE ATTACK v2: fixed (δ only in W[0])")
    print("=" * 70)

    N = 10000

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("1. RANDOM BASELINE (10K trials)")
    print("=" * 70)

    for n_r in [17, 18, 19, 20, 22, 24, 28, 32]:
        best = 256
        for _ in range(N):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            b = np.random.randint(0, 32)
            W2 = list(W); W2[0] ^= (1 << b)
            H1 = sha256_r(W, n_r); H2 = sha256_r(W2, n_r)
            dH = sum(hw(H1[i]^H2[i]) for i in range(8))
            if dH < best: best = dH
        print(f"  r={n_r:>2}: random best = {best:>3}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("2. A-REPAIR (δW[0] = 1 bit, repair W[1..15])")
    print("=" * 70)

    for n_r in [17, 18, 19, 20, 22, 24, 28, 32]:
        best = 256
        trivial_count = 0
        for _ in range(N):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            dW0 = 1 << np.random.randint(0, 32)
            dH, W2, is_trivial = a_repair_w0(W, dW0, n_r)
            if is_trivial:
                trivial_count += 1
                continue
            if dH < best: best = dH
        print(f"  r={n_r:>2}: repair best = {best:>3} (trivial: {trivial_count}/{N})")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("3. METRIC-GUIDED: δW in weakest bit of weakest word")
    print("=" * 70)

    # From metric analysis: W[15] bit 31 is weakest at r=17
    # W[0] bit 0 at r=20+ is same as random
    # Use absorption law: target W[r_target - 5] for maximum effect

    for n_r in [17, 18, 19, 20, 22, 24, 28, 32]:
        best = 256
        # Target word: the one that's WEAKEST at this round count
        # Weakest = last one to be absorbed = W[n_r - 6] roughly
        target_word = min(15, max(0, n_r - 6))

        for _ in range(N):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            b = np.random.randint(0, 32)
            W2 = list(W); W2[target_word] ^= (1 << b)
            H1 = sha256_r(W, n_r); H2 = sha256_r(W2, n_r)
            dH = sum(hw(H1[i]^H2[i]) for i in range(8))
            if dH < best: best = dH
        print(f"  r={n_r:>2}: metric-guided W[{target_word}] best = {best:>3}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("4. COMBINED: repair(W[0]) + metric-guided search")
    print("=" * 70)

    # a-repair with δW[0] = specific high bits (MSB more effective?)
    for n_r in [17, 18, 19, 20, 22, 24, 28, 32]:
        best_any = 256
        best_lsb = 256
        best_msb = 256

        for _ in range(N):
            W = [np.random.randint(0, 2**32) for _ in range(16)]

            # Try bit 31 (MSB)
            dH, _, triv = a_repair_w0(W, 1 << 31, n_r)
            if not triv and dH < best_msb: best_msb = dH

            # Try bit 0 (LSB)
            dH, _, triv = a_repair_w0(W, 1, n_r)
            if not triv and dH < best_lsb: best_lsb = dH

            # Try random bit
            dH, _, triv = a_repair_w0(W, 1 << np.random.randint(0, 32), n_r)
            if not triv and dH < best_any: best_any = dH

        print(f"  r={n_r:>2}: repair+MSB={best_msb:>3}, repair+LSB={best_lsb:>3}, "
              f"repair+random={best_any:>3}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("5. SUMMARY TABLE")
    print("=" * 70)

    print(f"\n  {'Round':>5} {'Random':>8} {'Repair':>8} {'Metric':>8} {'Advantage':>10}")

    for n_r in [17, 18, 19, 20, 22, 24, 28, 32]:
        # Random
        best_r = 256
        for _ in range(N):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            W2 = list(W); W2[np.random.randint(0,16)] ^= (1 << np.random.randint(0,32))
            H1 = sha256_r(W, n_r); H2 = sha256_r(W2, n_r)
            dH = sum(hw(H1[i]^H2[i]) for i in range(8))
            if dH < best_r: best_r = dH

        # Repair
        best_rep = 256
        for _ in range(N):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            dH, _, triv = a_repair_w0(W, 1 << np.random.randint(0, 32), n_r)
            if not triv and dH < best_rep: best_rep = dH

        # Metric-guided
        target_word = min(15, max(0, n_r - 6))
        best_m = 256
        for _ in range(N):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            W2 = list(W); W2[target_word] ^= (1 << np.random.randint(0, 32))
            H1 = sha256_r(W, n_r); H2 = sha256_r(W2, n_r)
            dH = sum(hw(H1[i]^H2[i]) for i in range(8))
            if dH < best_m: best_m = dH

        best_ours = min(best_rep, best_m)
        adv = best_r - best_ours
        bar = "█" * max(0, adv)
        print(f"  {n_r:>5} {best_r:>8} {best_rep:>8} {best_m:>8} {adv:>+9} {bar}")

    print(f"\n  Advantage = random_best - min(repair, metric)")
    print(f"  Positive = our techniques find BETTER near-collisions")


if __name__ == "__main__":
    main()
