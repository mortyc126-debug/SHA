"""
ULTIMATE ATTACK: все наши инструменты в одном pipeline.

Компоненты:
  1. METRIC: определяет СЛАБОЕ направление (eigenvector)
  2. MULTI-WORD: δW в нескольких слабых словах одновременно
  3. A-REPAIR: контролирует δ в rounds 1-15
  4. CONSERVATION FILTER: отсеивает плохие кандидаты рано
  5. MASSIVE SEARCH: перебор в сжатом пространстве

Target: максимальное число раундов где мы ЛУЧШЕ random search.
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


def sha256_states(W16, n_r):
    W = list(W16)
    for r in range(16, max(n_r, 16)):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))
    a,b,c,d,e,f,g,h = IV
    states = [(a,b,c,d,e,f,g,h)]
    for r in range(n_r):
        T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K[r]), W[r])
        T2 = add32(Sigma0(a), Maj(a,b,c))
        h,g,f,e = g,f,e,add32(d,T1)
        d,c,b,a = c,b,a,add32(T1,T2)
        states.append((a,b,c,d,e,f,g,h))
    final = tuple(add32(IV[i], states[n_r][i]) for i in range(8))
    return states, W, final


def a_repair_build(W_base, dW_mask, n_r):
    """a-repair: apply dW_mask to W_base, repair rounds 1-15, return full hash."""
    W2 = list(W_base)
    # Apply mask to multiple words
    for word in range(16):
        W2[word] ^= dW_mask[word]

    # a-repair rounds 1-15
    a1,b1,c1,d1,e1,f1,g1,h1 = IV
    a2,b2,c2,d2,e2,f2,g2,h2 = IV

    T1=add32(add32(add32(add32(h1,Sigma1(e1)),Ch(e1,f1,g1)),K[0]),W_base[0])
    T2=add32(Sigma0(a1),Maj(a1,b1,c1))
    h1,g1,f1,e1=g1,f1,e1,add32(d1,T1); d1,c1,b1,a1=c1,b1,a1,add32(T1,T2)

    T1=add32(add32(add32(add32(h2,Sigma1(e2)),Ch(e2,f2,g2)),K[0]),W2[0])
    T2=add32(Sigma0(a2),Maj(a2,b2,c2))
    h2,g2,f2,e2=g2,f2,e2,add32(d2,T1); d2,c2,b2,a2=c2,b2,a2,add32(T1,T2)

    for r in range(1, 16):
        T1_1=add32(add32(add32(add32(h1,Sigma1(e1)),Ch(e1,f1,g1)),K[r]),W_base[r])
        T2_1=add32(Sigma0(a1),Maj(a1,b1,c1))
        a1_new=add32(T1_1,T2_1); e1_new=add32(d1,T1_1)
        T2_2=add32(Sigma0(a2),Maj(a2,b2,c2))
        T1_n=sub32(a1_new,T2_2)
        W2[r]=sub32(sub32(sub32(sub32(T1_n,h2),Sigma1(e2)),Ch(e2,f2,g2)),K[r])
        h1,g1,f1,e1=g1,f1,e1,e1_new; d1,c1,b1,a1=c1,b1,a1,a1_new
        h2,g2,f2,e2=g2,f2,e2,add32(d2,T1_n); d2,c2,b2,a2=c2,b2,a2,a1_new

    H1 = sha256_r(W_base, n_r)
    H2 = sha256_r(W2, n_r)
    dH = sum(hw(H1[i]^H2[i]) for i in range(8))
    return dH, W2, H1, H2


def main():
    np.random.seed(42)

    print("=" * 70)
    print("ULTIMATE ATTACK: all techniques combined")
    print("=" * 70)

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("1. BASELINE: random search (no techniques)")
    print("=" * 70)

    N_TRIALS = 10000

    for n_r in [17, 18, 19, 20, 22, 24]:
        best = 256
        for _ in range(N_TRIALS):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            w = np.random.randint(0, 16)
            b = np.random.randint(0, 32)
            W2 = list(W); W2[w] ^= (1 << b)
            H1 = sha256_r(W, n_r)
            H2 = sha256_r(W2, n_r)
            dH = sum(hw(H1[i]^H2[i]) for i in range(8))
            if dH < best: best = dH
        print(f"  r={n_r:>2}: random best = {best:>3} (from {N_TRIALS} trials)")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("2. TECHNIQUE A: a-repair only")
    print("=" * 70)

    for n_r in [17, 18, 19, 20, 22, 24]:
        best = 256
        for _ in range(N_TRIALS):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            dW_mask = [0]*16
            # Random 1 bit in random word
            dW_mask[np.random.randint(0, 16)] = 1 << np.random.randint(0, 32)
            dH, _, _, _ = a_repair_build(W, dW_mask, n_r)
            if dH < best: best = dH
        print(f"  r={n_r:>2}: a-repair best = {best:>3}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("3. TECHNIQUE B: a-repair + target WEAK words (W[13-15])")
    print("=" * 70)

    for n_r in [17, 18, 19, 20, 22, 24]:
        best = 256
        for _ in range(N_TRIALS):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            dW_mask = [0]*16
            # Target weak words only
            word = np.random.choice([13, 14, 15])
            dW_mask[word] = 1 << np.random.randint(0, 32)
            dH, _, _, _ = a_repair_build(W, dW_mask, n_r)
            if dH < best: best = dH
        print(f"  r={n_r:>2}: repair+weak best = {best:>3}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("4. TECHNIQUE C: a-repair + MULTI-WORD δ in W[13-15]")
    print("=" * 70)

    for n_r in [17, 18, 19, 20, 22, 24]:
        best = 256
        for _ in range(N_TRIALS):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            dW_mask = [0]*16
            # Multi-word: 1 bit in each of W[13], W[14], W[15]
            for word in [13, 14, 15]:
                if np.random.random() < 0.7:  # 70% chance each
                    dW_mask[word] = 1 << np.random.randint(0, 32)
            if all(d == 0 for d in dW_mask):
                dW_mask[15] = 1  # fallback
            dH, _, _, _ = a_repair_build(W, dW_mask, n_r)
            if dH < best: best = dH
        print(f"  r={n_r:>2}: repair+multi best = {best:>3}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("5. TECHNIQUE D: FULL PIPELINE (repair + weak + conservation filter)")
    print("=" * 70)

    for n_r in [17, 18, 19, 20, 22, 24]:
        best = 256
        tested = 0
        filtered = 0

        for _ in range(N_TRIALS * 3):  # more candidates, filter aggressively
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            dW_mask = [0]*16

            # Target weak words
            for word in [13, 14, 15]:
                if np.random.random() < 0.5:
                    dW_mask[word] = 1 << np.random.randint(0, 32)
            if all(d == 0 for d in dW_mask):
                dW_mask[15] = 1

            # Quick a-repair + partial evaluation (conservation filter)
            W2 = list(W)
            for word in range(16): W2[word] ^= dW_mask[word]

            a1,b1,c1,d1,e1,f1,g1,h1 = IV
            a2,b2,c2,d2,e2,f2,g2,h2 = IV
            T1=add32(add32(add32(add32(h1,Sigma1(e1)),Ch(e1,f1,g1)),K[0]),W[0])
            T2=add32(Sigma0(a1),Maj(a1,b1,c1))
            h1,g1,f1,e1=g1,f1,e1,add32(d1,T1); d1,c1,b1,a1=c1,b1,a1,add32(T1,T2)
            T1=add32(add32(add32(add32(h2,Sigma1(e2)),Ch(e2,f2,g2)),K[0]),W2[0])
            T2=add32(Sigma0(a2),Maj(a2,b2,c2))
            h2,g2,f2,e2=g2,f2,e2,add32(d2,T1); d2,c2,b2,a2=c2,b2,a2,add32(T1,T2)

            for r in range(1, 16):
                T1_1=add32(add32(add32(add32(h1,Sigma1(e1)),Ch(e1,f1,g1)),K[r]),W[r])
                T2_1=add32(Sigma0(a1),Maj(a1,b1,c1))
                a1_new=add32(T1_1,T2_1); e1_new=add32(d1,T1_1)
                T2_2=add32(Sigma0(a2),Maj(a2,b2,c2))
                T1_n=sub32(a1_new,T2_2)
                W2[r]=sub32(sub32(sub32(sub32(T1_n,h2),Sigma1(e2)),Ch(e2,f2,g2)),K[r])
                h1,g1,f1,e1=g1,f1,e1,e1_new; d1,c1,b1,a1=c1,b1,a1,a1_new
                h2,g2,f2,e2=g2,f2,e2,add32(d2,T1_n); d2,c2,b2,a2=c2,b2,a2,a1_new

            # Conservation filter: check δ(a+e) at round 16
            # If small → promising
            ae1 = add32(a1, e1)
            ae2 = add32(a2, e2)
            cons_hw = hw(ae1 ^ ae2)

            if cons_hw < 10:  # strict filter
                filtered += 1
                dH, _, _, _ = a_repair_build(W, dW_mask, n_r)
                tested += 1
                if dH < best: best = dH

        print(f"  r={n_r:>2}: FULL pipeline best = {best:>3} "
              f"(filtered {filtered}/{N_TRIALS*3}, tested {tested})")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("6. COMPARISON TABLE")
    print("=" * 70)

    print(f"\n  Running final comparison (5K trials each)...")
    N_FINAL = 5000

    results = {}
    for n_r in [17, 18, 19, 20, 22, 24]:
        row = {}

        # Random
        best_r = 256
        for _ in range(N_FINAL):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            W2 = list(W); W2[np.random.randint(0,16)] ^= (1 << np.random.randint(0,32))
            H1 = sha256_r(W, n_r); H2 = sha256_r(W2, n_r)
            dH = sum(hw(H1[i]^H2[i]) for i in range(8))
            if dH < best_r: best_r = dH
        row['random'] = best_r

        # a-repair + weak words
        best_rw = 256
        for _ in range(N_FINAL):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            dW = [0]*16
            for word in [13, 14, 15]:
                if np.random.random() < 0.5:
                    dW[word] = 1 << np.random.randint(0, 32)
            if all(d==0 for d in dW): dW[15] = 1
            dH, _, _, _ = a_repair_build(W, dW, n_r)
            if dH < best_rw: best_rw = dH
        row['repair_weak'] = best_rw

        results[n_r] = row

    print(f"\n  {'Round':>5} {'Random':>8} {'Repair+Weak':>12} {'Advantage':>10}")
    print(f"  {'':>5} {'(5K)':>8} {'(5K)':>12} {'(bits)':>10}")
    for n_r in sorted(results.keys()):
        r = results[n_r]
        adv = r['random'] - r['repair_weak']
        bar = "█" * max(0, adv)
        print(f"  {n_r:>5} {r['random']:>8} {r['repair_weak']:>12} {adv:>+9} {bar}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("7. BEST RESULT: find actual near-collision")
    print("=" * 70)

    # Massive search at r=18 (where we have good advantage)
    n_r = 18
    best_ever = 256
    best_pair = None
    t0 = time.time()

    for trial in range(50000):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        dW = [0]*16
        for word in [13, 14, 15]:
            if np.random.random() < 0.5:
                dW[word] = 1 << np.random.randint(0, 32)
        if all(d==0 for d in dW): dW[15] = 1

        dH, W2, H1, H2 = a_repair_build(W, dW, n_r)
        if dH < best_ever:
            best_ever = dH
            best_pair = (list(W), list(W2), H1, H2, dH)

    elapsed = time.time() - t0

    print(f"  50K trials at r={n_r}, {elapsed:.1f}s:")
    print(f"  Best near-collision: δH = {best_ever} bits")
    print(f"  ({256 - best_ever} of 256 bits MATCH = {(256-best_ever)/256*100:.1f}%)")

    if best_pair:
        W1, W2, H1, H2, dH = best_pair
        print(f"\n  Message difference (non-zero words):")
        for i in range(16):
            if W1[i] != W2[i]:
                print(f"    δW[{i}] = {hex(W1[i] ^ W2[i])} (HW={hw(W1[i]^W2[i])})")
        print(f"\n  Hash 1: {' '.join(hex(h) for h in H1)}")
        print(f"  Hash 2: {' '.join(hex(h) for h in H2)}")
        print(f"  XOR:    {' '.join(hex(H1[i]^H2[i]) for i in range(8))}")
        matching_words = sum(1 for i in range(8) if H1[i] == H2[i])
        print(f"  Matching words: {matching_words}/8")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("ИТОГ")
    print("=" * 70)


if __name__ == "__main__":
    main()
