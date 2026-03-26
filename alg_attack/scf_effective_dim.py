#!/usr/bin/env python3
"""
EFFECTIVE DIMENSION of δH from a-repair.

Per-bit: 256. Closest pair: 97 (vs 239 random).
Вопрос: РЕАЛЬНАЯ размерность подпространства δH.

Методы:
1. GF(2) rank матрицы δH векторов
2. HW distribution (если d<256 → narrower distribution)
3. Closest pair scaling (d_closest vs N)
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

def dH_vec(w1):
    """δH as 256-bit integer."""
    w2 = a_repair(w1)
    h1=sha(w1); h2=sha(w2)
    v = 0
    for i in range(8):
        v |= (h1[i]^h2[i]) << (i*32)
    return v


# ============================================================
# Method 1: GF(2) rank of δH matrix
# ============================================================
def method1_gf2_rank(N):
    print("="*60)
    print(f"METHOD 1: GF(2) RANK of δH matrix (N={N})")
    print("="*60)

    rows = []
    for _ in range(N):
        w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        v = dH_vec(w1)
        rows.append(v)

    # GF(2) rank
    m = list(rows); rank = 0
    for bp in range(255, -1, -1):
        mask = 1 << bp
        pv = -1
        for i in range(rank, len(m)):
            if m[i] & mask: pv=i; break
        if pv == -1: continue
        m[rank],m[pv] = m[pv],m[rank]
        for i in range(len(m)):
            if i!=rank and m[i]&mask: m[i]^=m[rank]
        rank += 1

    print(f"  GF(2) rank of {N} δH vectors: {rank} / 256")

    if rank < 256:
        print(f"  ★ RANK DEFECT = {256-rank}!")
        print(f"  δH lives in a {rank}-dimensional subspace!")
        print(f"  Birthday: O(2^{rank//2}) instead of O(2^128)!")
    else:
        print(f"  Full rank — no GF(2) rank defect")

    return rank


# ============================================================
# Method 2: HW distribution of δH
# ============================================================
def method2_hw_dist(N):
    print(f"\n{'='*60}")
    print(f"METHOD 2: HW distribution of δH (N={N})")
    print("="*60)

    hws = []
    for _ in range(N):
        w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        v = dH_vec(w1)
        hws.append(HW(v))

    avg = sum(hws)/len(hws)
    std = (sum((h-avg)**2 for h in hws)/len(hws))**0.5

    # For random 256-bit: avg=128, std=8
    # For d-dimensional: avg=d/2, std≈√(d/4)
    est_d_from_avg = avg * 2
    est_d_from_std = (std**2) * 4

    print(f"  avg(HW) = {avg:.2f}  (random=128)")
    print(f"  std(HW) = {std:.2f}  (random=8.0)")
    print(f"  Effective d from avg: {est_d_from_avg:.0f}")
    print(f"  Effective d from std: {est_d_from_std:.0f}")


# ============================================================
# Method 3: Closest pair scaling
# ============================================================
def method3_closest_pair(N):
    print(f"\n{'='*60}")
    print(f"METHOD 3: Closest pair of δH (N={N})")
    print("="*60)

    samples = []
    for _ in range(min(N, 2000)):
        w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        samples.append(dH_vec(w1))

    # Closest pair (subsample for speed)
    best = 256
    n = len(samples)
    for i in range(min(n, 500)):
        for j in range(i+1, min(i+200, n)):
            d = HW(samples[i] ^ samples[j])
            if d < best: best = d

    print(f"  Closest pair (among {n} samples): dist = {best}")

    # Random 256-bit: expected closest ≈ 256 - 2*log₂(N*(N-1)/2)
    expected = 256 - int(2*math.log2(max(n*(n-1)//2, 1)))
    print(f"  Expected for random 256-bit: ~{expected}")
    print(f"  Ratio: {best/max(expected,1):.2f}")

    if best < expected * 0.7:
        print(f"  ★ SIGNIFICANTLY closer than random!")

    # Also: closest pair of random SHA hashes (control)
    ctrl_samples = []
    for _ in range(min(N, 2000)):
        w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        h = sha(w1)
        v = 0
        for i in range(8): v |= h[i]<<(i*32)
        ctrl_samples.append(v)

    ctrl_best = 256
    for i in range(min(len(ctrl_samples), 500)):
        for j in range(i+1, min(i+200, len(ctrl_samples))):
            d = HW(ctrl_samples[i]^ctrl_samples[j])
            if d < ctrl_best: ctrl_best = d

    print(f"\n  Control (random SHA hashes): closest = {ctrl_best}")
    print(f"  a-repair δH closest: {best}")
    print(f"  Difference: {ctrl_best - best:+d}")


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 500
    rank = method1_gf2_rank(min(N, 300))
    method2_hw_dist(N)
    method3_closest_pair(N)

    print(f"\n{'='*60}")
    print(f"CONCLUSION")
    print(f"  GF(2) rank: {rank}")
    if rank < 256:
        print(f"  ★ Effective dimension = {rank}")
        print(f"  ★ Birthday on a-repair residuals: O(2^{rank//2})")
    else:
        print(f"  Full rank 256 — but closest pair may reveal more")
