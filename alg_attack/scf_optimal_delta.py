#!/usr/bin/env python3
"""
SCF: OPTIMAL δ SEARCH — какой начальный δW даёт лучший near-collision?

До сих пор: δW[0] = 0x80000000 (MSB flip). Почему именно этот?
Wang использовал его, но мы можем ВЫЧИСЛИТЬ оптимальный.

Пространство поиска:
A. Single-bit: δW[word][bit] для всех 512 бит
B. Multi-bit: δW[0] = random 32-bit value
C. Multi-word: δW distributed across several words

Для каждого δ: запускаем v4 solver и измеряем best HW.
"""
import os, sys, math

MASK32 = 0xFFFFFFFF
K = [
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2,
]
IV = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
      0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&MASK32
def Sig0(x): return rotr(x,2)^rotr(x,13)^rotr(x,22)
def Sig1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def sig0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def sig1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def Ch(e,f,g): return (e&f)^(~e&g)&MASK32
def Maj(a,b,c): return (a&b)^(a&c)^(b&c)
def add32(*args):
    s=0
    for x in args: s=(s+x)&MASK32
    return s
def sub32(a,b): return (a-b)&MASK32
def hw(x): return bin(x).count('1')

def expand_real(W16):
    W=list(W16)
    for i in range(16,64):
        W.append(add32(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    return W

def sha_compress(W16):
    W=expand_real(W16); s=list(IV)
    for r in range(64):
        a,b,c,d,e,f,g,h=s
        T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add32(Sig0(a),Maj(a,b,c))
        s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
    return [add32(IV[i],s[i]) for i in range(8)]

def hash_hw(H1, W2):
    H2 = sha_compress(W2)
    return sum(hw(H1[i]^H2[i]) for i in range(8))

def quick_sa(W1, W2_init, H1, budget=500):
    """Quick SA polish."""
    W2 = list(W2_init)
    best = hash_hw(H1, W2)
    cur = list(W2); cur_s = best
    for it in range(budget):
        t = list(cur)
        w = int.from_bytes(os.urandom(1),'big') % 16
        b = int.from_bytes(os.urandom(1),'big') % 32
        t[w] ^= (1<<b)
        if t == W1: continue
        s = hash_hw(H1, t)
        T = max(0.01, 1-it/budget)
        if s < cur_s or math.exp(-(s-cur_s)/(T*2)) > int.from_bytes(os.urandom(4),'big')/(1<<32):
            cur = t; cur_s = s
            if s < best: best = s; W2 = list(t)
    return W2, best


# ============================================================
# EXP 1: Single-bit δ — which bit position gives best HW?
# ============================================================
def exp1_single_bit(W1, H1):
    print("="*70)
    print("EXP 1: SINGLE-BIT δ — best bit position")
    print("="*70)

    results = []
    for word in range(16):
        for bit in range(32):
            W2 = list(W1)
            W2[word] ^= (1 << bit)
            raw_hw = hash_hw(H1, W2)

            # Quick SA polish
            _, polished = quick_sa(W1, W2, H1, budget=200)
            results.append((polished, raw_hw, word, bit))

    results.sort()
    print(f"\n  Top 10 single-bit δ (with SA polish):")
    for polished, raw, word, bit in results[:10]:
        print(f"    W[{word:2d}][{bit:2d}]: polished={polished}, raw={raw}")

    print(f"\n  Bottom 3:")
    for polished, raw, word, bit in results[-3:]:
        print(f"    W[{word:2d}][{bit:2d}]: polished={polished}, raw={raw}")

    print(f"\n  MSB of W[0] (our default): ", end="")
    msb = next((p,r,w,b) for p,r,w,b in results if w==0 and b==31)
    print(f"polished={msb[0]}, raw={msb[1]}, rank={results.index(msb)+1}/512")

    best = results[0]
    return best[2], best[3], best[0]  # word, bit, hw


# ============================================================
# EXP 2: Multi-bit δW[0] — random 32-bit deltas
# ============================================================
def exp2_multibit(W1, H1, N):
    print("\n" + "="*70)
    print("EXP 2: MULTI-BIT δW[0] — random 32-bit values")
    print("="*70)

    results = []
    # Test specific patterns
    test_deltas = [
        0x80000000, 0x00000001, 0x00008000, 0x80008000,
        0x00000003, 0x00000007, 0x0000000F, 0x000000FF,
        0xFFFFFFFF, 0x55555555, 0xAAAAAAAA,
    ]
    # Add random deltas
    for _ in range(N):
        test_deltas.append(int.from_bytes(os.urandom(4), 'big') | 1)

    for delta in test_deltas:
        if delta == 0: continue
        W2 = list(W1); W2[0] ^= delta
        _, polished = quick_sa(W1, W2, H1, budget=200)
        results.append((polished, delta, hw(delta)))

    results.sort()
    print(f"\n  Top 10 δW[0] values:")
    for pol, delta, delta_hw in results[:10]:
        print(f"    δ=0x{delta:08x} (HW={delta_hw:2d}): polished={pol}")

    print(f"\n  0x80000000 (default): polished={next(p for p,d,_ in results if d==0x80000000)}")

    best_delta = results[0][1]
    best_hw_result = results[0][0]
    return best_delta, best_hw_result


# ============================================================
# EXP 3: Multi-word δ — spread across 2-4 words
# ============================================================
def exp3_multiword(W1, H1, N):
    print("\n" + "="*70)
    print("EXP 3: MULTI-WORD δ — spread across 2-4 words")
    print("="*70)

    results = []

    for trial in range(N):
        # Random 2-word delta
        W2 = list(W1)
        w1 = int.from_bytes(os.urandom(1),'big') % 16
        w2 = (w1 + 1 + int.from_bytes(os.urandom(1),'big') % 15) % 16
        W2[w1] ^= int.from_bytes(os.urandom(4),'big')
        W2[w2] ^= int.from_bytes(os.urandom(4),'big')
        if W2 == W1: continue

        _, polished = quick_sa(W1, W2, H1, budget=200)
        n_words = sum(1 for i in range(16) if W1[i] != W2[i])
        results.append((polished, n_words, w1, w2))

    results.sort()
    print(f"\n  Top 10 multi-word δ:")
    for pol, nw, w1, w2 in results[:10]:
        print(f"    words=[{w1},{w2}] ({nw} changed): polished={pol}")

    # Compare with single-word
    single_results = []
    for trial in range(N):
        W2 = list(W1)
        w = int.from_bytes(os.urandom(1),'big') % 16
        W2[w] ^= int.from_bytes(os.urandom(4),'big')
        if W2 == W1: continue
        _, pol = quick_sa(W1, W2, H1, budget=200)
        single_results.append(pol)

    avg_multi = sum(r[0] for r in results[:len(single_results)])/max(len(single_results),1)
    avg_single = sum(single_results)/max(len(single_results),1)

    print(f"\n  Average: multi-word={avg_multi:.1f}, single-word={avg_single:.1f}")
    print(f"  Multi-word advantage: {avg_single-avg_multi:+.1f}")

    return results[0]


# ============================================================
if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 50

    W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    H1 = sha_compress(W1)

    print("="*70)
    print("SCF: OPTIMAL δ SEARCH")
    print("="*70)

    best_word, best_bit, best_hw1 = exp1_single_bit(W1, H1)
    best_delta, best_hw2 = exp2_multibit(W1, H1, N)
    best_multi = exp3_multiword(W1, H1, N)

    print(f"\n{'='*70}")
    print(f"SUMMARY")
    print(f"  Best single-bit: W[{best_word}][{best_bit}] → HW={best_hw1}")
    print(f"  Best multi-bit δW[0]: 0x{best_delta:08x} → HW={best_hw2}")
    print(f"  Best multi-word: HW={best_multi[0]}")
    print(f"  Default (MSB W[0]): HW≈100-104")
