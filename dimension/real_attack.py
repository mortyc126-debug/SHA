"""
РЕАЛЬНАЯ АТАКА: применяем знания измерения в стандартной математике.

Из измерения мы знаем:
  1. W[15] при r=16 влияет на 4 бита (H[0] и H[4] only)
  2. a-repair: δa=0 forcing → δe converges to 0 in 4 rounds
  3. Pipe cascade: a→b→c→d, e→f→g→h (1 round each)
  4. Schedule: W[16]=σ1(W[14])+W[9]+σ0(W[1])+W[0]

СТРАТЕГИЯ (Wang-style message modification + наши знания):
  1. Выбираем δW с МИНИМАЛЬНЫМ hamming weight
  2. Используем a-repair чтобы δ не росло в первых 16 раундах
  3. Контролируем schedule чтобы δW[16+] были малыми
  4. Target: reduced-round collision (r=20-24)
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


def sha256_n_rounds(W16, n_rounds):
    """SHA-256 truncated to n_rounds."""
    W = list(W16)
    for r in range(16, max(n_rounds, 16)):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))
    a,b,c,d,e,f,g,h = IV
    for r in range(n_rounds):
        T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K[r]), W[r])
        T2 = add32(Sigma0(a), Maj(a,b,c))
        h,g,f,e = g,f,e,add32(d,T1)
        d,c,b,a = c,b,a,add32(T1,T2)
    return tuple(add32(IV[i], [a,b,c,d,e,f,g,h][i]) for i in range(8))


def expand_schedule(W16):
    W = list(W16)
    for r in range(16, 64):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))
    return W


def main():
    np.random.seed(42)

    print("=" * 70)
    print("РЕАЛЬНАЯ АТАКА: знания измерения → differential trails")
    print("=" * 70)

    # ═══════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("1. СТРАТЕГИЯ 1: δW concentrated in W[15], target r=17")
    print("   (absorption law: W[15] at r=16 affects only 4 bits)")
    print("=" * 70)

    # At r=17: W[15] has ~24 bits influence (from absorption law)
    # δW[15] = 1 bit → δH at r=17 should be ~12 bits
    # Search: which bit of W[15] gives LOWEST δH?

    best_results = {}
    for n_r in [17, 18, 19, 20]:
        best_hw_r = 256
        best_bit = -1
        for bit in range(32):
            # Try many random base messages
            total_hw = 0
            n_tries = 200
            for _ in range(n_tries):
                W = [np.random.randint(0, 2**32) for _ in range(16)]
                W2 = list(W); W2[15] ^= (1 << bit)
                H1 = sha256_n_rounds(W, n_r)
                H2 = sha256_n_rounds(W2, n_r)
                total_hw += sum(hw(H1[i]^H2[i]) for i in range(8))

            avg = total_hw / n_tries
            if avg < best_hw_r:
                best_hw_r = avg
                best_bit = bit

        best_results[n_r] = (best_bit, best_hw_r)
        print(f"  r={n_r}: best bit={best_bit}, avg HW(δH)={best_hw_r:.1f}")

    # ═══════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("2. СТРАТЕГИЯ 2: a-repair + specific δW[0]")
    print("   (build differential path using repair)")
    print("=" * 70)

    # Choose δW[0] = 1 bit
    # Use a-repair to control first rounds
    # Let schedule propagate naturally
    # Measure: how far can we keep δ small?

    for start_bit in [0, 15, 31]:
        dW0 = 1 << start_bit
        controlled_rounds = []

        for trial in range(100):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            W2 = list(W)
            W2[0] ^= dW0

            # Forward both, track δ per round
            a1,b1,c1,d1,e1,f1,g1,h1 = IV
            a2,b2,c2,d2,e2,f2,g2,h2 = IV

            # Round 0: free, both use their own W[0]
            for state_pair, w_val in [
                ([a1,b1,c1,d1,e1,f1,g1,h1], W[0]),
            ]:
                pass

            T1 = add32(add32(add32(add32(h1,Sigma1(e1)),Ch(e1,f1,g1)),K[0]),W[0])
            T2 = add32(Sigma0(a1),Maj(a1,b1,c1))
            h1,g1,f1,e1 = g1,f1,e1,add32(d1,T1)
            d1,c1,b1,a1 = c1,b1,a1,add32(T1,T2)

            T1 = add32(add32(add32(add32(h2,Sigma1(e2)),Ch(e2,f2,g2)),K[0]),W2[0])
            T2 = add32(Sigma0(a2),Maj(a2,b2,c2))
            h2,g2,f2,e2 = g2,f2,e2,add32(d2,T1)
            d2,c2,b2,a2 = c2,b2,a2,add32(T1,T2)

            # Rounds 1-15: a-repair (force δa=0)
            for r in range(1, 16):
                T1_1 = add32(add32(add32(add32(h1,Sigma1(e1)),Ch(e1,f1,g1)),K[r]),W[r])
                T2_1 = add32(Sigma0(a1),Maj(a1,b1,c1))
                a1_new = add32(T1_1, T2_1)
                e1_new = add32(d1, T1_1)

                # Force a2_new = a1_new
                T2_2 = add32(Sigma0(a2),Maj(a2,b2,c2))
                T1_needed = sub32(a1_new, T2_2)
                W2[r] = sub32(sub32(sub32(sub32(T1_needed,h2),Sigma1(e2)),Ch(e2,f2,g2)),K[r])

                h1,g1,f1,e1 = g1,f1,e1,e1_new
                d1,c1,b1,a1 = c1,b1,a1,a1_new

                e2_new = add32(d2, T1_needed)
                h2,g2,f2,e2 = g2,f2,e2,e2_new
                d2,c2,b2,a2 = c2,b2,a2,a1_new

            # Now: W2[0..15] are set. Compute full hash for various rounds.
            for n_r in [17, 18, 20, 24]:
                H1 = sha256_n_rounds(W, n_r)
                H2 = sha256_n_rounds(W2, n_r)
                dH = sum(hw(H1[i]^H2[i]) for i in range(8))
                if trial == 0:
                    controlled_rounds.append((n_r, dH))

        # Average over 100 trials for each round count
        results_by_round = {}
        for trial in range(100):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            W2 = list(W); W2[0] ^= dW0

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
                e2_n=add32(d2,T1_n)
                h2,g2,f2,e2=g2,f2,e2,e2_n; d2,c2,b2,a2=c2,b2,a2,a1_new

            for n_r in [17, 18, 20, 24, 32, 64]:
                H1 = sha256_n_rounds(W, n_r)
                H2 = sha256_n_rounds(W2, n_r)
                dH = sum(hw(H1[i]^H2[i]) for i in range(8))
                if n_r not in results_by_round:
                    results_by_round[n_r] = []
                results_by_round[n_r].append(dH)

        print(f"\n  δW[0] = bit {start_bit} + a-repair(r=1..15):")
        for n_r in sorted(results_by_round.keys()):
            vals = results_by_round[n_r]
            print(f"    r={n_r:>2}: mean HW(δH)={np.mean(vals):>6.1f}, "
                  f"min={min(vals):>3}, max={max(vals):>3}")

    # ═══════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("3. СТРАТЕГИЯ 3: schedule-aware δW")
    print("   (choose δW[0..15] to minimize δW[16..23])")
    print("=" * 70)

    # δW[16] = σ1(δW[14]) + δW[9] + σ0(δW[1]) + δW[0]  (XOR approx)
    # To minimize δW[16]: set δW only in words that DON'T affect W[16]
    # W[16] depends on W[14], W[9], W[1], W[0]
    # W[17] depends on W[15], W[10], W[2], W[1]
    # ...
    # "Safe" words for W[16]: W[2..8, 10..13, 15]
    # But W[17] uses W[15], W[10], W[2] → removes them
    # Intersection of "safe" for W[16] AND W[17]: W[3..8, 11..13]

    print(f"  Schedule dependencies:")
    for r in range(16, 24):
        deps = [r-2, r-7, r-15, r-16]
        print(f"    W[{r}] depends on W[{deps}]")

    # Find: which input words affect FEWEST schedule words 16-23?
    impact = {}
    for w in range(16):
        affected = set()
        for r in range(16, 24):
            if w in [r-2, r-7, r-15, r-16]:
                affected.add(r)
        impact[w] = affected

    print(f"\n  Input word impact on W[16..23]:")
    for w in range(16):
        print(f"    W[{w:>2}]: affects {sorted(impact[w])} ({len(impact[w])}/8)")

    # Best word: affects fewest schedule words in 16-23
    best_word = min(impact, key=lambda w: len(impact[w]))
    print(f"\n  Best word for minimal schedule impact: W[{best_word}] (affects {len(impact[best_word])})")

    # Test: δ in best_word vs worst_word
    worst_word = max(impact, key=lambda w: len(impact[w]))

    for word, label in [(best_word, "BEST"), (worst_word, "WORST")]:
        dH_list = []
        for trial in range(500):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            W2 = list(W); W2[word] ^= 1
            H1 = sha256_n_rounds(W, 20)
            H2 = sha256_n_rounds(W2, 20)
            dH_list.append(sum(hw(H1[i]^H2[i]) for i in range(8)))

        print(f"  {label} word W[{word}] at r=20: mean δH={np.mean(dH_list):.1f}, min={min(dH_list)}")

    # ═══════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("4. КОМБИНИРОВАННАЯ АТАКА: repair + schedule-aware + search")
    print("=" * 70)

    # Combined:
    # 1. δ in W[best_word] (minimal schedule impact)
    # 2. a-repair for rounds 1-15 (converge δ to 0)
    # 3. Random search over base message for best δH

    target_rounds = 20
    best_global = 256
    best_W = None
    best_W2 = None
    n_search = 5000

    t0 = time.time()
    for trial in range(n_search):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W)
        W2[best_word] ^= (1 << np.random.randint(0, 32))

        # a-repair
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

        H1 = sha256_n_rounds(W, target_rounds)
        H2 = sha256_n_rounds(W2, target_rounds)
        dH = sum(hw(H1[i]^H2[i]) for i in range(8))

        if dH < best_global:
            best_global = dH
            best_W = list(W)
            best_W2 = list(W2)

    elapsed = time.time() - t0

    print(f"  Combined attack at r={target_rounds}:")
    print(f"    Search: {n_search} trials, {elapsed:.1f}s")
    print(f"    Best HW(δH): {best_global}")
    print(f"    Random baseline at r={target_rounds}: ~128")

    if best_global < 128:
        print(f"    → {128 - best_global} bits below random! Near-collision.")
        # Verify
        H1 = sha256_n_rounds(best_W, target_rounds)
        H2 = sha256_n_rounds(best_W2, target_rounds)
        print(f"    H1: {' '.join(hex(h) for h in H1)}")
        print(f"    H2: {' '.join(hex(h) for h in H2)}")
        print(f"    δH: {' '.join(hex(H1[i]^H2[i]) for i in range(8))}")

    # Also test without repair (pure random)
    best_random = 256
    for trial in range(n_search):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W); W2[best_word] ^= (1 << np.random.randint(0, 32))
        H1 = sha256_n_rounds(W, target_rounds)
        H2 = sha256_n_rounds(W2, target_rounds)
        dH = sum(hw(H1[i]^H2[i]) for i in range(8))
        if dH < best_random: best_random = dH

    print(f"\n    Without repair (pure random search): best={best_random}")
    print(f"    With repair+schedule-aware:           best={best_global}")
    print(f"    Advantage: {best_random - best_global} bits")

    # ═══════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("5. SCALING: attack quality vs rounds")
    print("=" * 70)

    for n_r in [17, 18, 19, 20, 22, 24, 28, 32]:
        best_rep = 256
        best_rnd = 256
        for trial in range(2000):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            W2 = list(W); W2[best_word] ^= (1 << np.random.randint(0, 32))

            # With repair
            a1,b1,c1,d1,e1,f1,g1,h1 = IV
            a2,b2,c2,d2,e2,f2,g2,h2 = IV
            T1=add32(add32(add32(add32(h1,Sigma1(e1)),Ch(e1,f1,g1)),K[0]),W[0])
            T2=add32(Sigma0(a1),Maj(a1,b1,c1))
            h1,g1,f1,e1=g1,f1,e1,add32(d1,T1); d1,c1,b1,a1=c1,b1,a1,add32(T1,T2)
            T1=add32(add32(add32(add32(h2,Sigma1(e2)),Ch(e2,f2,g2)),K[0]),W2[0])
            T2=add32(Sigma0(a2),Maj(a2,b2,c2))
            h2,g2,f2,e2=g2,f2,e2,add32(d2,T1); d2,c2,b2,a2=c2,b2,a2,add32(T1,T2)
            W2r = list(W2)
            for r in range(1, 16):
                T1_1=add32(add32(add32(add32(h1,Sigma1(e1)),Ch(e1,f1,g1)),K[r]),W[r])
                T2_1=add32(Sigma0(a1),Maj(a1,b1,c1))
                a1_new=add32(T1_1,T2_1); e1_new=add32(d1,T1_1)
                T2_2=add32(Sigma0(a2),Maj(a2,b2,c2))
                T1_n=sub32(a1_new,T2_2)
                W2r[r]=sub32(sub32(sub32(sub32(T1_n,h2),Sigma1(e2)),Ch(e2,f2,g2)),K[r])
                h1,g1,f1,e1=g1,f1,e1,e1_new; d1,c1,b1,a1=c1,b1,a1,a1_new
                h2,g2,f2,e2=g2,f2,e2,add32(d2,T1_n); d2,c2,b2,a2=c2,b2,a2,a1_new

            H1=sha256_n_rounds(W,n_r); H2=sha256_n_rounds(W2r,n_r)
            dH_rep = sum(hw(H1[i]^H2[i]) for i in range(8))
            if dH_rep < best_rep: best_rep = dH_rep

            # Without repair
            H1=sha256_n_rounds(W,n_r); H2=sha256_n_rounds(W2,n_r)
            dH_rnd = sum(hw(H1[i]^H2[i]) for i in range(8))
            if dH_rnd < best_rnd: best_rnd = dH_rnd

        advantage = best_rnd - best_rep
        print(f"  r={n_r:>2}: repair={best_rep:>3}, random={best_rnd:>3}, advantage={advantage:>+3} bits")

    # ═══════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("ИТОГ")
    print("=" * 70)


if __name__ == "__main__":
    main()
