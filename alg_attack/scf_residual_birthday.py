#!/usr/bin/env python3
"""
BIRTHDAY НА RESIDUALS: δH от a-repair — сколько бит энтропии?

a-repair с фиксированным break_pos и break_delta создаёт δH.
Разные W₁ → разные δH. Но δH СТРУКТУРИРОВАН (18 NEAR раундов).

Если H(δH) < 256 → birthday на δH дешевле чем 2^{128}.

И: если два a-repair residuals СОВПАДАЮТ (δH₁ = δH₂) →
multiblock collision: M₁||M₃ и M₂||M₃ дают одинаковый хеш
(если block 2 одинаков для обоих).

Нет! Для collision: нужно sha(M₁)=sha(M₂). Не δH₁=δH₂.

Но: если мы найдём M₁,M₂ с sha(M₁)=sha(M₂):
  → a-repair: W₂ = a_repair(W₁)
  → sha(W₁) и sha(W₂) отличаются на δH (≈89 бит)
  → НЕ collision!

Birthday нужен на ПОЛНОМ хеше, не на δH.
НО: a-repair определяет W₂ из W₁. Значит sha(W₂) = f(W₁).
Collision = sha(W₂) = sha(W₂') для ДВУХ РАЗНЫХ W₁.
  = f(W₁) = f(W₁'). Birthday на f!

f: 512 бит → 256 бит. Birthday на f = 2^{128}. Стандартно.

НО: a-repair W₂ зависит от W₁ ДЕТЕРМИНИРОВАНО.
sha(a_repair(W₁)) — это ДРУГАЯ функция g: 512 → 256.
Birthday на g тоже 2^{128}. Не помогает.

ОДНАКО: что если g(W₁) НЕ имеет 256 бит энтропии?
Что если a-repair ограничивает выход? Тогда birthday дешевле.

Измерю: H(sha(a_repair(W₁))) vs H(sha(W₁)).
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


def measure_residual_entropy(N):
    """Измерю: per-bit entropy δH при фиксированном break."""
    print("="*60)
    print("ENTROPY OF δH FROM a-repair")
    print("="*60)

    def binary_entropy(p):
        if p<=0 or p>=1: return 0
        return -p*math.log2(p)-(1-p)*math.log2(1-p)

    bit_counts = [[0]*32 for _ in range(8)]

    for _ in range(N):
        w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        w2 = a_repair(w1, br=3, delta=1)
        h1=sha(w1); h2=sha(w2)
        dH = [h1[i]^h2[i] for i in range(8)]
        for reg in range(8):
            for bit in range(32):
                if (dH[reg]>>bit)&1: bit_counts[reg][bit] += 1

    total_h = 0
    for reg in range(8):
        for bit in range(32):
            p = bit_counts[reg][bit]/N
            total_h += binary_entropy(p)

    print(f"  N = {N}")
    print(f"  Per-bit entropy of δH: {total_h:.2f} / 256")
    print(f"  Deficit: {256-total_h:.2f} bits")
    print(f"  Birthday on δH: O(2^{total_h/2:.1f})")

    # Also: what's the closest pair of δH values?
    print(f"\n  Closest pair search among {min(N,500)} δH values:")
    residuals = []
    for trial in range(min(N, 500)):
        w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        w2 = a_repair(w1, br=3, delta=1)
        h1=sha(w1); h2=sha(w2)
        residuals.append(tuple(h1[i]^h2[i] for i in range(8)))

    best_dist = 256
    for i in range(len(residuals)):
        for j in range(i+1, min(i+100, len(residuals))):
            d = sum(HW(residuals[i][k]^residuals[j][k]) for k in range(8))
            if d < best_dist:
                best_dist = d

    print(f"  Closest δH pair: dist = {best_dist}")
    print(f"  Expected for random: ~{256 - int(2*math.log2(max(len(residuals),2)))}")

    # Birthday on sha(a_repair(W₁)) — collision within a-repair family
    print(f"\n  Birthday on sha(W₂) where W₂=a_repair(W₁):")
    seen = {}
    for trial in range(min(N, 5000)):
        w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        w2 = a_repair(w1, br=3, delta=1)
        h2 = sha(w2)
        trunc = h2[0]  # First 32 bits
        if trunc in seen:
            old_h = seen[trunc][1]
            full_d = sum(HW(h2[i]^old_h[i]) for i in range(8))
            if full_d < 200:
                print(f"    Partial match at trial {trial}: full dist={full_d}")
        else:
            seen[trunc] = (w1, h2)
    print(f"    {len(seen)} unique 32-bit prefixes from {min(N,5000)} trials")


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
    measure_residual_entropy(N)
