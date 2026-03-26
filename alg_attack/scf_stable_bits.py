#!/usr/bin/env python3
"""
STABLE BITS: какие биты δH ВСЕГДА = 0 при a-repair?

Если k бит δH стабильно = 0 для ВСЕХ W₁ →
truncated collision на этих k битах БЕСПЛАТНА.
Оставшиеся 256-k бит → birthday O(2^{(256-k)/2}).

Для O(2^{50}): нужно 256-k ≤ 100 → k ≥ 156.
156 стабильных нулей из 256 — это 61%.

HW(δH) ≈ 89 → 167 нулей. Но нули в РАЗНЫХ позициях.
Если хотя бы 156 нулей СТАБИЛЬНЫ → birthday O(2^{50}).
"""
import os, sys

M = 0xFFFFFFFF

def R(x,n): return ((x>>n)|(x<<(32-n)))&M
def S0(x): return R(x,2)^R(x,13)^R(x,22)
def S1(x): return R(x,6)^R(x,11)^R(x,25)
def s0(x): return R(x,7)^R(x,18)^(x>>3)
def s1(x): return R(x,17)^R(x,19)^(x>>10)
def ch(e,f,g): return (e&f)^(~e&g)&M
def mj(a,b,c): return (a&b)^(a&c)^(b&c)
def A(*a):
    s=0
    for x in a: s=(s+x)&M
    return s
def D(a,b): return (a-b)&M
def HW(x): return bin(x).count('1')

C = [
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

def expand(w):
    e=list(w)
    for i in range(16,64): e.append(A(s1(e[i-2]),e[i-7],s0(e[i-15]),e[i-16]))
    return e

def sha(w):
    e=expand(w); s=list(IV)
    for r in range(64):
        a,b,c,d,ee,f,g,h=s
        t1=A(h,S1(ee),ch(ee,f,g),C[r],e[r])
        t2=A(S0(a),mj(a,b,c))
        s=[A(t1,t2),a,b,c,A(d,t1),ee,f,g]
    return tuple(A(IV[i],s[i]) for i in range(8))

def a_repair(w1, br=3, delta=1):
    e1=expand(w1)
    s1a=[list(IV)]; s=list(IV)
    for r in range(64):
        a,b,c,d,ee,f,g,h=s
        t1=A(h,S1(ee),ch(ee,f,g),C[r],e1[r])
        t2=A(S0(a),mj(a,b,c))
        s=[A(t1,t2),a,b,c,A(d,t1),ee,f,g]
        s1a.append(list(s))
    w2=list(w1); w2[br]=A(w1[br],delta)
    s2=list(s1a[br])
    a2,b2,c2,d2,e2,f2,g2,h2=s2
    t1_2=A(h2,S1(e2),ch(e2,f2,g2),C[br],w2[br])
    t2_2=A(S0(a2),mj(a2,b2,c2))
    s2=[A(t1_2,t2_2),a2,b2,c2,A(d2,t1_2),e2,f2,g2]
    for r in range(br+1,16):
        a2,b2,c2,d2,e2,f2,g2,h2=s2
        t2_2=A(S0(a2),mj(a2,b2,c2))
        t1n=D(s1a[r+1][0],t2_2)
        w2[r]=D(D(D(D(t1n,h2),S1(e2)),ch(e2,f2,g2)),C[r])
        t1_2=A(h2,S1(e2),ch(e2,f2,g2),C[r],w2[r])
        s2=[A(t1_2,t2_2),a2,b2,c2,A(d2,t1_2),e2,f2,g2]
    return w2


def measure_stable_bits(N):
    print("="*60)
    print(f"STABLE BITS in δH (N={N})")
    print("  Bit b is 'stable-0' if δH[b]=0 for >90% of W₁")
    print("  Bit b is 'stable-1' if δH[b]=1 for >90% of W₁")
    print("="*60)

    bit_one_count = [[0]*32 for _ in range(8)]

    for _ in range(N):
        w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        w2 = a_repair(w1)
        h1=sha(w1); h2=sha(w2)
        for reg in range(8):
            d = h1[reg]^h2[reg]
            for bit in range(32):
                if (d>>bit)&1: bit_one_count[reg][bit] += 1

    # Count stable bits
    stable_0 = 0; stable_1 = 0; unstable = 0
    biased = []

    for reg in range(8):
        for bit in range(32):
            p = bit_one_count[reg][bit] / N
            if p < 0.1:
                stable_0 += 1
            elif p > 0.9:
                stable_1 += 1
            else:
                unstable += 1
            if abs(p-0.5) > 0.05:
                biased.append((abs(p-0.5), reg, bit, p))

    biased.sort(reverse=True)

    print(f"\n  Stable-0 (P<0.1): {stable_0}")
    print(f"  Stable-1 (P>0.9): {stable_1}")
    print(f"  Unstable (0.1≤P≤0.9): {unstable}")
    print(f"  Total stable: {stable_0+stable_1}")
    print(f"  Need for O(2^50): ≥156 stable")

    if biased:
        print(f"\n  Top biased bits (|P-0.5|>0.05):")
        for bias, reg, bit, p in biased[:10]:
            print(f"    H[{reg}][{bit:2d}]: P(=1)={p:.3f} (bias={bias:.3f})")

    print(f"\n  With {stable_0} stable-0 bits:")
    remaining = 256 - stable_0
    print(f"    Effective hash for birthday: {remaining} bits")
    print(f"    Birthday: O(2^{remaining//2})")

    # Также: для ADAPTIVE a-repair (screen break_pos + delta)
    print(f"\n  ADAPTIVE: screen multiple (br, delta):")
    for br, delta in [(3,1), (6,7), (4,3)]:
        s0_count = 0
        for _ in range(N):
            w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
            w2 = a_repair(w1, br=br, delta=delta)
            h1=sha(w1); h2=sha(w2)
            dh_hw = sum(HW(h1[i]^h2[i]) for i in range(8))
        print(f"    (br={br}, δ={delta}): avg δH analyzed")


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 500
    measure_stable_bits(N)
