#!/usr/bin/env python3
"""
LATE BREAK: break на r=8 вместо r=0.

Rounds 0..7: W₂ = W₁ (идентичны, state одинаков)
Round 8: break (δW[8] = 1)
Rounds 9..15: a-repair (7 раундов)

Преимущество: δW[0..7] = 0! Значит schedule deps от W[0..7] = 0!
δW[16] = σ1(δW[14]) + δW[9] + σ0(δW[1]) + δW[0]
       = σ1(δW[14]) + δW[9] + 0 + 0
       = σ1(δW[14]) + δW[9] + carry!

Только 2 слагаемых (vs 4 при break на r=0)!
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

def sha(w):
    e=expand(w); s=list(IV)
    for r in range(64):
        a,b,c,d,ee,f,g,h=s
        t1=A(h,S1(ee),ch(ee,f,g),C[r],e[r])
        t2=A(S0(a),mj(a,b,c))
        s=[A(t1,t2),a,b,c,A(d,t1),ee,f,g]
    return tuple(A(IV[i],s[i]) for i in range(8))


def late_break_repair(w1, break_round, break_delta=1):
    """Break на раунде break_round. a-repair на break_round+1..15."""

    e1 = expand(w1)
    s1_all = [list(IV)]
    s = list(IV)
    for r in range(64):
        a,b,c,d,ee,f,g,h = s
        t1=A(h,S1(ee),ch(ee,f,g),C[r],e1[r])
        t2=A(S0(a),mj(a,b,c))
        s=[A(t1,t2),a,b,c,A(d,t1),ee,f,g]
        s1_all.append(list(s))

    w2 = list(w1)
    # Break: change W[break_round]
    w2[break_round] = A(w1[break_round], break_delta)

    # Run stream 2: identical to stream 1 until break_round
    s2 = list(IV)
    for r in range(break_round):
        a2,b2,c2,d2,e2,f2,g2,h2 = s2
        t1_2=A(h2,S1(e2),ch(e2,f2,g2),C[r],e1[r])
        t2_2=A(S0(a2),mj(a2,b2,c2))
        s2=[A(t1_2,t2_2),a2,b2,c2,A(d2,t1_2),e2,f2,g2]

    # Break round
    a2,b2,c2,d2,e2,f2,g2,h2 = s2
    t1_2=A(h2,S1(e2),ch(e2,f2,g2),C[break_round],w2[break_round])
    t2_2=A(S0(a2),mj(a2,b2,c2))
    s2=[A(t1_2,t2_2),a2,b2,c2,A(d2,t1_2),e2,f2,g2]

    # a-repair: break_round+1 .. 15
    for r in range(break_round+1, 16):
        a2,b2,c2,d2,e2,f2,g2,h2 = s2
        a_target = s1_all[r+1][0]
        t2_2 = A(S0(a2),mj(a2,b2,c2))
        t1_need = D(a_target, t2_2)
        w2[r] = D(D(D(D(t1_need,h2),S1(e2)),ch(e2,f2,g2)),C[r])
        t1_2=A(h2,S1(e2),ch(e2,f2,g2),C[r],w2[r])
        s2=[A(t1_2,t2_2),a2,b2,c2,A(d2,t1_2),e2,f2,g2]

    return w2


def compare_break_positions(N):
    print("="*60)
    print("LATE BREAK: сравнение позиций break")
    print("="*60)

    for trial in range(N):
        w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        e1 = expand(w1)

        print(f"\n  Trial {trial}:")
        print(f"  {'break_r':>7} | δW[0..15] zeros | ��W[16] | ��W[17] | dH")
        print("  " + "-"*60)

        for break_r in range(16):
            w2 = late_break_repair(w1, break_r, break_delta=1)
            e2 = expand(w2)

            dw = [HW(e1[r]^e2[r]) for r in range(64)]
            n_zeros_015 = sum(1 for r in range(16) if e1[r]==e2[r])

            h1=sha(w1); h2=sha(w2)
            dh=sum(HW(h1[i]^h2[i]) for i in range(8))
            ndiff = sum(1 for i in range(16) if w1[i]!=w2[i])

            # BOTH=0 раунды?
            s1=list(IV); s2=list(IV)
            both0 = 0
            for r in range(64):
                a1,b1,c1,d1,ee1,f1,g1,h1r = s1
                t1_1=A(h1r,S1(ee1),ch(ee1,f1,g1),C[r],e1[r])
                t2_1=A(S0(a1),mj(a1,b1,c1))
                s1=[A(t1_1,t2_1),a1,b1,c1,A(d1,t1_1),ee1,f1,g1]
                a2,b2,c2,d2,ee2,f2,g2,h2r = s2
                t1_2=A(h2r,S1(ee2),ch(ee2,f2,g2),C[r],e2[r])
                t2_2=A(S0(a2),mj(a2,b2,c2))
                s2=[A(t1_2,t2_2),a2,b2,c2,A(d2,t1_2),ee2,f2,g2]
                if s1==s2: both0+=1

            star = ""
            if dw[16] <= 5: star += " ★16"
            if dw[17] <= 5: star += " ★17"
            if both0 > 20: star += f" ★★BOTH={both0}"

            print(f"  r={break_r:2d}    | {n_zeros_015:2d}/16 zeros   "
                  f"| {dw[16]:6d} | {dw[17]:6d} | {dh:3d} (diff={ndiff}){star}")


if __name__ == '__main__':
    import sys
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    compare_break_positions(N)
