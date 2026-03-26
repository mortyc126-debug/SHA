#!/usr/bin/env python3
"""
Break r=3 — ГЛУБОКОЕ исследование. Как далеко цепочка BOTH=0?

δW[0..2]=0, δW[3]=break, δW[4..15]=a-repair → cascade.
δW[16]=0, δW[17]=0 (все 4 deps нулевые!).
δW[18]=σ0(δW[3]) (единственное ненулевое слагаемое).

Вопрос: ВСЯ карта δW[16..63] для break r=3.
И: СКОЛЬКО раундов BOTH=0?
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


def break_r3(w1, break_delta=1):
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
    w2[3] = A(w1[3], break_delta)  # Break at r=3

    # Rounds 0..2: identical
    s2 = list(IV)
    for r in range(3):
        a2,b2,c2,d2,e2,f2,g2,h2 = s2
        t1_2=A(h2,S1(e2),ch(e2,f2,g2),C[r],e1[r])
        t2_2=A(S0(a2),mj(a2,b2,c2))
        s2=[A(t1_2,t2_2),a2,b2,c2,A(d2,t1_2),e2,f2,g2]

    # Round 3: break
    a2,b2,c2,d2,e2,f2,g2,h2 = s2
    t1_2=A(h2,S1(e2),ch(e2,f2,g2),C[3],w2[3])
    t2_2=A(S0(a2),mj(a2,b2,c2))
    s2=[A(t1_2,t2_2),a2,b2,c2,A(d2,t1_2),e2,f2,g2]

    # Rounds 4..15: a-repair
    for r in range(4, 16):
        a2,b2,c2,d2,e2,f2,g2,h2 = s2
        a_target = s1_all[r+1][0]
        t2_2 = A(S0(a2),mj(a2,b2,c2))
        t1_need = D(a_target, t2_2)
        w2[r] = D(D(D(D(t1_need,h2),S1(e2)),ch(e2,f2,g2)),C[r])
        t1_2=A(h2,S1(e2),ch(e2,f2,g2),C[r],w2[r])
        s2=[A(t1_2,t2_2),a2,b2,c2,A(d2,t1_2),e2,f2,g2]

    return w2


def deep_analysis(N):
    print("="*60)
    print("BREAK r=3 — DEEP ANALYSIS")
    print("="*60)

    for trial in range(N):
        w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

        for break_delta in [1, 3, 7]:
            w2 = break_r3(w1, break_delta)
            e1 = expand(w1); e2 = expand(w2)

            # Full δW map
            dw = [HW(e1[r]^e2[r]) for r in range(64)]
            is_zero = [e1[r]==e2[r] for r in range(64)]

            # Full BOTH=0 check
            s1=list(IV); s2=list(IV)
            both0_rounds = []
            for r in range(64):
                a1,b1,c1,d1,ee1,f1,g1,h1r = s1
                t1_1=A(h1r,S1(ee1),ch(ee1,f1,g1),C[r],e1[r])
                t2_1=A(S0(a1),mj(a1,b1,c1))
                s1=[A(t1_1,t2_1),a1,b1,c1,A(d1,t1_1),ee1,f1,g1]
                a2,b2,c2,d2,ee2,f2,g2,h2r = s2
                t1_2=A(h2r,S1(ee2),ch(ee2,f2,g2),C[r],e2[r])
                t2_2=A(S0(a2),mj(a2,b2,c2))
                s2=[A(t1_2,t2_2),a2,b2,c2,A(d2,t1_2),ee2,f2,g2]
                if s1 == s2:
                    both0_rounds.append(r+1)

            h1=sha(w1); h2=sha(w2)
            dh=sum(HW(h1[i]^h2[i]) for i in range(8))

            # Consecutive BOTH=0 from start
            consec = 0
            for r in range(1, 65):
                if r in both0_rounds:
                    consec += 1
                else:
                    break

            print(f"\n  Trial {trial}, δ={break_delta}:")
            print(f"    BOTH=0 rounds: {len(both0_rounds)}/64 (consecutive from start: {consec})")
            if len(both0_rounds) <= 30:
                print(f"    Which: {both0_rounds}")
            print(f"    dH = {dh}")

            # δW map summary
            zeros = [r for r in range(64) if is_zero[r]]
            print(f"    δW=0: {zeros}")
            print(f"    δW[16..20]: {dw[16:21]}")
            print(f"    δW[21..25]: {dw[21:26]}")

            if consec >= 18:
                print(f"    ★★★ {consec} consecutive BOTH=0!")


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    deep_analysis(N)
