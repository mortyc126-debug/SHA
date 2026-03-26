#!/usr/bin/env python3
"""
КАСКАДНАЯ КАРТА: какие δW[r] обнуляются через a-repair?

Schedule[r] = σ1(W[r-2]) + W[r-7] + σ0(W[r-15]) + W[r-16]

Для r=16..63: δW[r] з��висит от δW[r-2], δW[r-7], δW[r-15], δW[r-16].
Если какие-то из них = 0 → меньше слагаемых → проще обнулить.

Я построю ПОЛНУЮ карту: δW[0..15] → какие нулевые → как это
каскадирует в δW[16..63].
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

def a_repair(w1, break_delta=1):
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
    w2[0] = A(w1[0], break_delta)
    s2 = list(IV)
    a2,b2,c2,d2,e2,f2,g2,h2 = s2
    t1_2=A(h2,S1(e2),ch(e2,f2,g2),C[0],w2[0])
    t2_2=A(S0(a2),mj(a2,b2,c2))
    s2=[A(t1_2,t2_2),a2,b2,c2,A(d2,t1_2),e2,f2,g2]

    for r in range(1, 16):
        a2,b2,c2,d2,e2,f2,g2,h2 = s2
        a_target = s1_all[r+1][0]
        t2_2 = A(S0(a2),mj(a2,b2,c2))
        t1_need = D(a_target, t2_2)
        w2[r] = D(D(D(D(t1_need,h2),S1(e2)),ch(e2,f2,g2)),C[r])
        t1_2=A(h2,S1(e2),ch(e2,f2,g2),C[r],w2[r])
        s2=[A(t1_2,t2_2),a2,b2,c2,A(d2,t1_2),e2,f2,g2]

    return w2


def full_map(N):
    print("="*60)
    print("КАСКАДНАЯ КАРТА δW[0..63] после a-repair")
    print("="*60)

    for trial in range(N):
        w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        w2 = a_repair(w1, break_delta=1)

        e1 = expand(w1); e2 = expand(w2)
        dw = [HW(e1[r]^e2[r]) for r in range(64)]
        is_zero = [e1[r]==e2[r] for r in range(64)]

        print(f"\n  Trial {trial}:")
        print(f"  r  : HW(δW) zero? | schedule deps")
        print("  " + "-"*55)

        for r in range(64):
            z = "=0" if is_zero[r] else f"={dw[r]:2d}"

            deps = ""
            if r >= 16:
                d1, d2, d3, d4 = r-2, r-7, r-15, r-16
                z1 = "0" if is_zero[d1] else f"{dw[d1]}"
                z2 = "0" if is_zero[d2] else f"{dw[d2]}"
                z3 = "0" if is_zero[d3] else f"{dw[d3]}"
                z4 = "0" if is_zero[d4] else f"{dw[d4]}"
                n_zero_deps = sum([is_zero[d1],is_zero[d2],is_zero[d3],is_zero[d4]])
                deps = f"σ1(W[{d1}]:{z1})+W[{d2}]:{z2}+σ0(W[{d3}]:{z3})+W[{d4}]:{z4} [{n_zero_deps}/4 zero]"

            star = ""
            if is_zero[r]: star = " ★"
            elif r >= 16:
                n_zero = sum([is_zero[r-2],is_zero[r-7],is_zero[r-15] if r>=15 else False,is_zero[r-16] if r>=16 else False])
                if n_zero >= 3: star = " ◆ (3/4 deps zero!)"
                elif n_zero >= 2: star = " ○ (2/4 deps zero)"

            if r < 16 or dw[r] < 10 or is_zero[r] or star:
                print(f"  {r:2d} : δW{z:4s}    | {deps}{star}")

        # Считаю: сколько δW[16..63] = 0?
        n_zero_sched = sum(1 for r in range(16,64) if is_zero[r])
        n_low = sum(1 for r in range(16,64) if dw[r] <= 5)
        print(f"\n  δW=0 in schedule: {n_zero_sched}/48")
        print(f"  δW≤5 in schedule: {n_low}/48")


if __name__ == '__main__':
    full_map(3)
