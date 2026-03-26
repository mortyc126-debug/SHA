#!/usr/bin/env python3
"""
Birthday на W₁: найти W₁ где a-repair с break r=3 даёт δW[9]=0.

δW[9]=0 → δW[16] = σ0(0) + 0 = carry only.
При δW[9]=0 И δW[14]=0 И δW[1]=0 И δW[0]=0:
  carry от четырёх нулей = 0!
  → δW[16] = 0 ТОЧНО!

Потом δW[17] deps: δW[15]=0, δW[10]=?, δW[2]=0, δW[1]=0.
Если δW[10]=0 тоже → δW[17]=0 → chain extends дальше!

Перебираю W₁ до нахождения δW[9]=0.
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


def do_a_repair_r3(w1, break_delta=1):
    """Break r=3, a-repair r=4..15. Return w2."""
    e1 = expand(w1)
    s1_all = [list(IV)]
    s = list(IV)
    for r in range(64):
        a,b,c,d,ee,f,g,h = s
        t1=A(h,S1(ee),ch(ee,f,g),C[r],e1[r])
        t2=A(S0(a),mj(a,b,c))
        s=[A(t1,t2),a,b,c,A(d,t1),ee,f,g]
        s1_all.append(list(s))

    w2 = list(w1); w2[3] = A(w1[3], break_delta)
    s2 = list(IV)
    for r in range(3):
        a2,b2,c2,d2,e2,f2,g2,h2 = s2
        t1_2=A(h2,S1(e2),ch(e2,f2,g2),C[r],e1[r])
        t2_2=A(S0(a2),mj(a2,b2,c2))
        s2=[A(t1_2,t2_2),a2,b2,c2,A(d2,t1_2),e2,f2,g2]

    a2,b2,c2,d2,e2,f2,g2,h2 = s2
    t1_2=A(h2,S1(e2),ch(e2,f2,g2),C[3],w2[3])
    t2_2=A(S0(a2),mj(a2,b2,c2))
    s2=[A(t1_2,t2_2),a2,b2,c2,A(d2,t1_2),e2,f2,g2]

    for r in range(4, 16):
        a2,b2,c2,d2,e2,f2,g2,h2 = s2
        a_target = s1_all[r+1][0]
        t2_2 = A(S0(a2),mj(a2,b2,c2))
        t1_need = D(a_target, t2_2)
        w2[r] = D(D(D(D(t1_need,h2),S1(e2)),ch(e2,f2,g2)),C[r])
        t1_2=A(h2,S1(e2),ch(e2,f2,g2),C[r],w2[r])
        s2=[A(t1_2,t2_2),a2,b2,c2,A(d2,t1_2),e2,f2,g2]

    return w2


def search_w1_for_zero_dw9(N):
    """Перебираю W₁: ищу такой что δW[9]=0 при break r=3."""
    print("="*60)
    print(f"BIRTHDAY: ищу W₁ с δW[9]=0 (N={N})")
    print("="*60)

    best_dw9 = 32
    best_w1 = None
    best_full = None

    for trial in range(N):
        w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        w2 = do_a_repair_r3(w1, break_delta=1)
        e1 = expand(w1); e2 = expand(w2)

        dw9 = HW(e1[9] ^ e2[9])

        if dw9 < best_dw9:
            best_dw9 = dw9
            best_w1 = list(w1)

            # Full analysis
            dw = [HW(e1[r]^e2[r]) for r in range(64)]
            zeros = [r for r in range(20) if e1[r]==e2[r]]

            if dw9 <= 1:
                # Check chain extension
                h1=sha(w1); h2=sha(w2)
                dh=sum(HW(h1[i]^h2[i]) for i in range(8))
                n_diff=sum(1 for i in range(16) if w1[i]!=w2[i])

                print(f"  [{trial:6d}] δW[9]={dw9} ★")
                print(f"    zeros: {zeros}")
                print(f"    δW[16]={dw[16]}, δW[17]={dw[17]}, δW[18]={dw[18]}")
                print(f"    dH={dh}, diff={n_diff}")

                if dw9 == 0:
                    # Full BOTH=0 check
                    s1c=list(IV); s2c=list(IV)
                    both0=[]
                    for r in range(64):
                        a1,b1,c1,d1,ee1,f1,g1,h1r=s1c
                        t1_1=A(h1r,S1(ee1),ch(ee1,f1,g1),C[r],e1[r])
                        t2_1=A(S0(a1),mj(a1,b1,c1))
                        s1c=[A(t1_1,t2_1),a1,b1,c1,A(d1,t1_1),ee1,f1,g1]
                        a2r,b2r,c2r,d2r,ee2r,f2r,g2r,h2r=s2c
                        t1_2=A(h2r,S1(ee2r),ch(ee2r,f2r,g2r),C[r],e2[r])
                        t2_2=A(S0(a2r),mj(a2r,b2r,c2r))
                        s2c=[A(t1_2,t2_2),a2r,b2r,c2r,A(d2r,t1_2),ee2r,f2r,g2r]
                        if s1c==s2c: both0.append(r+1)
                    print(f"    BOTH=0: {both0[:30]}")
                    print(f"    Total BOTH=0: {len(both0)}")

    print(f"\n  Best δW[9] = {best_dw9} in {N} trials")

    if best_dw9 == 0:
        print(f"  ★★★ δW[9]=0 FOUND!")
    elif best_dw9 <= 2:
        print(f"  ★★ Close: δW[9]={best_dw9}")
    else:
        print(f"  δW[9]={best_dw9} — cascade not deep enough")
        # Estimate: P(δW[9]=0) from distribution
        print(f"  Need more trials or different strategy")


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10000
    search_w1_for_zero_dw9(N)
