#!/usr/bin/env python3
"""
SPONTANEOUS MERGE: P(δa[r]=0) для каждого раунда после a-repair.

Flow view: поток расходится на r=18 (break₂).
Вопрос: какова P что a-поток СПОНТАННО сливается на r=19,20,...?

Если P > 2^{-32} → спонтанный merge ДЕШЕВЛЕ чем random.
Если P растёт с a-repair (vs no a-repair) → a-repair ПОМОГАЕТ.
"""
import os

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

def a_repair_r3(w1, delta=1):
    e1=expand(w1)
    s1a=[list(IV)]; s=list(IV)
    for r in range(64):
        a,b,c,d,ee,f,g,h=s
        t1=A(h,S1(ee),ch(ee,f,g),C[r],e1[r])
        t2=A(S0(a),mj(a,b,c))
        s=[A(t1,t2),a,b,c,A(d,t1),ee,f,g]
        s1a.append(list(s))
    w2=list(w1); w2[3]=A(w1[3],delta)
    s2=list(s1a[3])
    a2,b2,c2,d2,e2,f2,g2,h2=s2
    t1_2=A(h2,S1(e2),ch(e2,f2,g2),C[3],w2[3])
    t2_2=A(S0(a2),mj(a2,b2,c2))
    s2=[A(t1_2,t2_2),a2,b2,c2,A(d2,t1_2),e2,f2,g2]
    for r in range(4,16):
        a2,b2,c2,d2,e2,f2,g2,h2=s2
        t2_2=A(S0(a2),mj(a2,b2,c2))
        t1n=D(s1a[r+1][0],t2_2)
        w2[r]=D(D(D(D(t1n,h2),S1(e2)),ch(e2,f2,g2)),C[r])
        t1_2=A(h2,S1(e2),ch(e2,f2,g2),C[r],w2[r])
        s2=[A(t1_2,t2_2),a2,b2,c2,A(d2,t1_2),e2,f2,g2]
    return w2


def measure_spontaneous(N):
    print("="*60)
    print(f"SPONTANEOUS MERGE: P(δa[r]=0) per round (N={N})")
    print("="*60)

    # With a-repair
    counts_ar = [0]*65

    # Without a-repair (just MSB flip)
    counts_raw = [0]*65

    for _ in range(N):
        w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

        # a-repair
        w2_ar = a_repair_r3(w1)
        e1=expand(w1); e2=expand(w2_ar)
        s1=list(IV); s2=list(IV)
        for r in range(64):
            a1,b1,c1,d1,ee1,f1,g1,h1=s1
            s1=[A(A(h1,S1(ee1),ch(ee1,f1,g1),C[r],e1[r]),A(S0(a1),mj(a1,b1,c1))),
                a1,b1,c1,A(d1,A(h1,S1(ee1),ch(ee1,f1,g1),C[r],e1[r])),ee1,f1,g1]
            a2,b2,c2,d2,ee2,f2,g2,h2=s2
            s2=[A(A(h2,S1(ee2),ch(ee2,f2,g2),C[r],e2[r]),A(S0(a2),mj(a2,b2,c2))),
                a2,b2,c2,A(d2,A(h2,S1(ee2),ch(ee2,f2,g2),C[r],e2[r])),ee2,f2,g2]
            if s1[0]==s2[0]: counts_ar[r+1] += 1

        # Raw MSB flip
        w2_raw = list(w1); w2_raw[0] ^= 0x80000000
        e2r=expand(w2_raw)
        s1=list(IV); s2=list(IV)
        for r in range(64):
            a1,b1,c1,d1,ee1,f1,g1,h1=s1
            s1=[A(A(h1,S1(ee1),ch(ee1,f1,g1),C[r],e1[r]),A(S0(a1),mj(a1,b1,c1))),
                a1,b1,c1,A(d1,A(h1,S1(ee1),ch(ee1,f1,g1),C[r],e1[r])),ee1,f1,g1]
            a2,b2,c2,d2,ee2,f2,g2,h2=s2
            s2=[A(A(h2,S1(ee2),ch(ee2,f2,g2),C[r],e2r[r]),A(S0(a2),mj(a2,b2,c2))),
                a2,b2,c2,A(d2,A(h2,S1(ee2),ch(ee2,f2,g2),C[r],e2r[r])),ee2,f2,g2]
            if s1[0]==s2[0]: counts_raw[r+1] += 1

    print(f"\n  {'r':>3} | {'a-repair':>10} {'raw MSB':>10} | {'ratio':>6}")
    print("  " + "-"*40)

    for r in range(1, 65):
        p_ar = counts_ar[r]/N
        p_raw = counts_raw[r]/N
        ratio = p_ar/max(p_raw, 1/N) if p_raw > 0 else (N if p_ar > 0 else 1)

        show = r <= 20 or r >= 55 or p_ar > 0 or p_raw > 0
        if show:
            star = ""
            if p_ar > 0 and r > 18: star = " ★ SPONTANEOUS!"
            if p_ar > 0.01 and r > 18: star = " ★★ FREQUENT!"
            print(f"  {r:3d} | {p_ar:10.4f} {p_raw:10.4f} | {ratio:6.1f}x{star}")


if __name__ == '__main__':
    import sys
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
    measure_spontaneous(N)
