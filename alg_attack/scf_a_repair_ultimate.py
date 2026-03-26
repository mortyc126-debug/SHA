#!/usr/bin/env python3
"""
A-REPAIR как стартовая точка для оптимизации.

a-repair даёт: 10 BOTH=0, 18 NEAR≤10, dH≈120.
SA от random start даёт: dH≈100.

Вопрос: SA от a-repair start → dH < 100?

a-repair создаёт СТРУКТУРИРОВАННЫЙ W₂ (не random):
  - δW[0,1,2]=0 (до break)
  - δW[3]=1 (маленький break)
  - δW[4..15]=a-repair corrections (структурированные!)

Может быть a-repair start ЛУЧШЕ для SA чем random start?
"""
import os, sys, math

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

def dh(w1, w2):
    h1=sha(w1); h2=sha(w2)
    return sum(HW(h1[i]^h2[i]) for i in range(8))


def do_a_repair_r3(w1, break_delta=1):
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
    for r in range(3): s2 = list(s1_all[r+1])  # Same as stream 1
    # Break round 3
    a2,b2,c2,d2,e2,f2,g2,h2 = s1_all[3]
    t1_2=A(h2,S1(e2),ch(e2,f2,g2),C[3],w2[3])
    t2_2=A(S0(a2),mj(a2,b2,c2))
    s2=[A(t1_2,t2_2),a2,b2,c2,A(d2,t1_2),e2,f2,g2]
    for r in range(4, 16):
        a2,b2,c2,d2,e2,f2,g2,h2 = s2
        t2_2 = A(S0(a2),mj(a2,b2,c2))
        t1_need = D(s1_all[r+1][0], t2_2)
        w2[r] = D(D(D(D(t1_need,h2),S1(e2)),ch(e2,f2,g2)),C[r])
        t1_2=A(h2,S1(e2),ch(e2,f2,g2),C[r],w2[r])
        s2=[A(t1_2,t2_2),a2,b2,c2,A(d2,t1_2),e2,f2,g2]
    return w2


def sa_polish(w1, w2_start, budget):
    """SA от заданной стартовой точки."""
    best = dh(w1, w2_start); cur = list(w2_start); best_w = list(w2_start)
    for it in range(budget):
        t = list(cur)
        w = int.from_bytes(os.urandom(1),'big') % 16
        b = int.from_bytes(os.urandom(1),'big') % 32
        t[w] ^= (1<<b)
        if t == w1: continue
        s = dh(w1, t)
        T = max(0.01, 1-it/budget)
        if s < best or math.exp(-(s-best)/(T*2)) > int.from_bytes(os.urandom(4),'big')/(1<<32):
            cur = t
            if s < best: best = s; best_w = list(t)
    return best_w, best


def experiment(N):
    print("="*60)
    print("A-REPAIR + SA vs PURE SA")
    print("="*60)

    a_repair_results = []
    pure_sa_results = []
    a_repair_sa_results = []

    budget = 2000

    for trial in range(N):
        w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

        # 1. Pure a-repair (no SA)
        w2_ar = do_a_repair_r3(w1)
        dh_ar = dh(w1, w2_ar)
        a_repair_results.append(dh_ar)

        # 2. a-repair + SA polish
        w2_ar_sa, dh_ar_sa = sa_polish(w1, w2_ar, budget)
        a_repair_sa_results.append(dh_ar_sa)

        # 3. Pure SA (MSB delta start)
        w2_sa_start = list(w1); w2_sa_start[0] ^= 0x80000000
        w2_sa, dh_sa = sa_polish(w1, w2_sa_start, budget)
        pure_sa_results.append(dh_sa)

        winner = "A+SA ★" if dh_ar_sa < dh_sa else ("SA" if dh_sa < dh_ar_sa else "tie")
        if trial < 10 or dh_ar_sa < 90:
            print(f"  [{trial:3d}] a-repair={dh_ar:3d} → a-repair+SA={dh_ar_sa:3d}  pure_SA={dh_sa:3d}  {winner}")

    avg_ar = sum(a_repair_results)/N
    avg_ar_sa = sum(a_repair_sa_results)/N
    avg_sa = sum(pure_sa_results)/N

    print(f"\n{'='*60}")
    print(f"  a-repair only:    avg={avg_ar:.1f}  min={min(a_repair_results)}")
    print(f"  a-repair + SA:    avg={avg_ar_sa:.1f}  min={min(a_repair_sa_results)}")
    print(f"  pure SA:          avg={avg_sa:.1f}  min={min(pure_sa_results)}")
    print(f"")
    wins = sum(1 for i in range(N) if a_repair_sa_results[i] < pure_sa_results[i])
    print(f"  a-repair+SA wins: {wins}/{N}")
    print(f"  Advantage: {avg_sa - avg_ar_sa:+.1f}")


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    experiment(N)
