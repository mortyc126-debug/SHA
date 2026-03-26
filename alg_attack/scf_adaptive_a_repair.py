#!/usr/bin/env python3
"""
ADAPTIVE A-REPAIR: screen break_position × break_delta → SA polish.
Объединяю ВСЕ наши инструменты в одном.
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

def hdiff(w1,w2):
    h1=sha(w1);h2=sha(w2)
    return sum(HW(h1[i]^h2[i]) for i in range(8))


def a_repair_general(w1, break_round, break_delta):
    """a-repair с произвольной позицией break и delta."""
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
    w2[break_round] = A(w1[break_round], break_delta)

    s2 = list(IV)
    for r in range(break_round):
        s2 = list(s1_all[r+1])  # Same as stream 1

    # Break round
    a2,b2,c2,d2,e2,f2,g2,h2 = s1_all[break_round]
    t1_2=A(h2,S1(e2),ch(e2,f2,g2),C[break_round],w2[break_round])
    t2_2=A(S0(a2),mj(a2,b2,c2))
    s2=[A(t1_2,t2_2),a2,b2,c2,A(d2,t1_2),e2,f2,g2]

    # a-repair rounds
    for r in range(break_round+1, 16):
        a2,b2,c2,d2,e2,f2,g2,h2 = s2
        t2_2 = A(S0(a2),mj(a2,b2,c2))
        t1_need = D(s1_all[r+1][0], t2_2)
        w2[r] = D(D(D(D(t1_need,h2),S1(e2)),ch(e2,f2,g2)),C[r])
        t1_2=A(h2,S1(e2),ch(e2,f2,g2),C[r],w2[r])
        s2=[A(t1_2,t2_2),a2,b2,c2,A(d2,t1_2),e2,f2,g2]

    return w2


def adaptive_a_repair(w1, budget=3000):
    """Screen break positions and deltas, SA polish best."""

    # Phase 1: Screen (200 evals)
    candidates = []
    for br in range(1, 15):
        for delta in [1, 3, 7, 15, 0xFF]:
            w2 = a_repair_general(w1, br, delta)
            d = hdiff(w1, w2)
            candidates.append((d, br, delta, list(w2)))

    # Also screen plain single-bit deltas (no a-repair)
    for word in range(16):
        for bit in [0, 7, 15, 23, 31]:
            w2 = list(w1); w2[word] ^= (1<<bit)
            d = hdiff(w1, w2)
            candidates.append((d, -1, (word, bit), list(w2)))

    candidates.sort()
    screen_evals = len(candidates)

    # Phase 2: SA polish top-5 (remaining budget)
    sa_budget = (budget - screen_evals) // 5
    global_best = 256; global_w2 = None; global_type = ""

    for d, br, delta, w2 in candidates[:5]:
        cur = list(w2); best = hdiff(w1, cur)
        for it in range(sa_budget):
            t = list(cur)
            w = int.from_bytes(os.urandom(1),'big') % 16
            b = int.from_bytes(os.urandom(1),'big') % 32
            t[w] ^= (1<<b)
            if t == w1: continue
            s = hdiff(w1, t)
            T = max(0.01, 1-it/sa_budget)
            if s < best or math.exp(-(s-best)/(T*2)) > int.from_bytes(os.urandom(4),'big')/(1<<32):
                cur = t
                if s < best: best = s

        if best < global_best:
            global_best = best
            global_w2 = list(cur)
            global_type = f"a-repair(r={br},δ={delta})" if br >= 0 else f"bit({delta})"

    return global_w2, global_best, global_type


def compare(N, budget=3000):
    print("="*60)
    print(f"ADAPTIVE A-REPAIR vs ALL METHODS (budget={budget})")
    print("="*60)

    adaptive_ar = []; pure_sa = []; adaptive_bit = []

    for trial in range(N):
        w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

        # Adaptive a-repair
        _, d_ar, tp = adaptive_a_repair(w1, budget)
        adaptive_ar.append(d_ar)

        # Pure SA (MSB)
        w2_sa = list(w1); w2_sa[0] ^= 0x80000000
        cur = list(w2_sa); best_sa = hdiff(w1, cur)
        for it in range(budget):
            t = list(cur)
            w = int.from_bytes(os.urandom(1),'big') % 16
            b = int.from_bytes(os.urandom(1),'big') % 32
            t[w] ^= (1<<b)
            if t == w1: continue
            s = hdiff(w1, t)
            T = max(0.01,1-it/budget)
            if s<best_sa or math.exp(-(s-best_sa)/(T*2))>int.from_bytes(os.urandom(4),'big')/(1<<32):
                cur=t
                if s<best_sa: best_sa=s
        pure_sa.append(best_sa)

        # Adaptive bit (no a-repair, just screen bits + SA)
        screen = []
        for word in range(8):
            for bit in range(32):
                w2=list(w1); w2[word]^=(1<<bit)
                screen.append((hdiff(w1,w2), word, bit))
        screen.sort()
        sa_b = (budget-256)//3
        best_ab = 256
        for _,wd,bt in screen[:3]:
            w2=list(w1); w2[wd]^=(1<<bt)
            cur=list(w2); b=hdiff(w1,cur)
            for it in range(sa_b):
                t=list(cur); w=int.from_bytes(os.urandom(1),'big')%16
                bb=int.from_bytes(os.urandom(1),'big')%32
                if w==wd and bb==bt: continue
                t[w]^=(1<<bb)
                if t==w1: continue
                s=hdiff(w1,t)
                T=max(0.01,1-it/sa_b)
                if s<b or math.exp(-(s-b)/(T*2))>int.from_bytes(os.urandom(4),'big')/(1<<32):
                    cur=t
                    if s<b: b=s
            if b < best_ab: best_ab = b
        adaptive_bit.append(best_ab)

        winner = "AR" if d_ar <= best_sa and d_ar <= best_ab else ("SA" if best_sa <= best_ab else "AB")
        print(f"  [{trial:2d}] adaptive_a-repair={d_ar:3d}({tp[:15]:15s}) SA={best_sa:3d} adaptive_bit={best_ab:3d} → {winner}")

    print(f"\n{'='*60}")
    print(f"  Adaptive a-repair: avg={sum(adaptive_ar)/N:.1f}  min={min(adaptive_ar)}")
    print(f"  Pure SA:           avg={sum(pure_sa)/N:.1f}  min={min(pure_sa)}")
    print(f"  Adaptive bit:      avg={sum(adaptive_bit)/N:.1f}  min={min(adaptive_bit)}")


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    compare(N)
