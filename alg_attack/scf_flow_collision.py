#!/usr/bin/env python3
"""
FLOW COLLISION — collision в измерении потоков.

Поток P(W) проходит через SHA: P: W → H.
Collision: P(W₁) = P(W₂).

В потоковом измерении: два потока СЛИВАЮТСЯ на выходе.
Слияние = потоки стали НЕРАЗЛИЧИМЫ.

Подпотоки: H = (H₀, H₁, ..., H₇). Каждый — отдельный подпоток.
Полная collision = 8 подпотоков совпали.

Мой подход: заставить подпотоки СЛИВАТЬСЯ ПО ОДНОМУ.

Подпоток H₀ зависит от state[64][0] = a[64].
a[64] — конец a-sequence.
a-sequence — одномерный поток.

Если два a-потока СЛИЛИСЬ к r=64 → H₀ совпадает.
a-repair показал: a-потоки сливаются на r=12 (reboot).
Если они остаются слитыми до r=64 → H₀ collision!

Остальные H₁..H₇: зависят от a[63..61] и e[64..61].
e = a - D (D-coupling). Если a-поток слит → e тоже
(если D одинаковые, что зависит от a[-4..0]).

Мой инструмент: FLOW MERGER.
Заставляю a-поток₂ = a-поток₁ начиная с r=R_merge.
Если R_merge ≤ 52 (= 64 - 12 для D-coupling) → collision.
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


def flow_observe(w1, w2):
    """Наблюдаю два потока: ГДЕ они сливаются и расходятся.
    Не считаю HW — смотрю на a-поток и e-поток ОТДЕЛЬНО."""

    e1 = expand(w1); e2 = expand(w2)
    s1 = list(IV); s2 = list(IV)

    a_merged = []  # раунды где a₁ = a₂
    e_merged = []  # раунды где e₁ = e₂
    full_merged = []  # раунды где state₁ = state₂

    for r in range(64):
        a1,b1,c1,d1,ee1,f1,g1,h1 = s1
        a2,b2,c2,d2,ee2,f2,g2,h2 = s2

        t1_1=A(h1,S1(ee1),ch(ee1,f1,g1),C[r],e1[r])
        t2_1=A(S0(a1),mj(a1,b1,c1))
        s1=[A(t1_1,t2_1),a1,b1,c1,A(d1,t1_1),ee1,f1,g1]

        t1_2=A(h2,S1(ee2),ch(ee2,f2,g2),C[r],e2[r])
        t2_2=A(S0(a2),mj(a2,b2,c2))
        s2=[A(t1_2,t2_2),a2,b2,c2,A(d2,t1_2),ee2,f2,g2]

        if s1[0] == s2[0]: a_merged.append(r+1)
        if s1[4] == s2[4]: e_merged.append(r+1)
        if s1 == s2: full_merged.append(r+1)

    # Hash comparison per sub-stream
    h1 = sha(w1); h2 = sha(w2)
    sub_match = [h1[i] == h2[i] for i in range(8)]

    return {
        'a_merged': a_merged,
        'e_merged': e_merged,
        'full_merged': full_merged,
        'sub_match': sub_match,
        'n_sub_match': sum(sub_match),
    }


def flow_merger():
    """Инструмент: заставляю потоки сливаться."""

    print("="*60)
    print("FLOW MERGER: наблюдаю слияние подпотоков")
    print("="*60)

    # a-repair break r=3
    for trial in range(5):
        w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

        # a-repair
        e1 = expand(w1)
        s1a = [list(IV)]; s = list(IV)
        for r in range(64):
            a,b,c,d,ee,f,g,h = s
            t1=A(h,S1(ee),ch(ee,f,g),C[r],e1[r])
            t2=A(S0(a),mj(a,b,c))
            s=[A(t1,t2),a,b,c,A(d,t1),ee,f,g]
            s1a.append(list(s))

        w2 = list(w1); w2[3] = A(w1[3], 1)
        s2 = list(s1a[3])
        a2,b2,c2,d2,e2,f2,g2,h2 = s2
        t1_2=A(h2,S1(e2),ch(e2,f2,g2),C[3],w2[3])
        t2_2=A(S0(a2),mj(a2,b2,c2))
        s2=[A(t1_2,t2_2),a2,b2,c2,A(d2,t1_2),e2,f2,g2]
        for r in range(4, 16):
            a2,b2,c2,d2,e2,f2,g2,h2 = s2
            t2_2 = A(S0(a2),mj(a2,b2,c2))
            t1n = D(s1a[r+1][0], t2_2)
            w2[r] = D(D(D(D(t1n,h2),S1(e2)),ch(e2,f2,g2)),C[r])
            t1_2=A(h2,S1(e2),ch(e2,f2,g2),C[r],w2[r])
            s2=[A(t1_2,t2_2),a2,b2,c2,A(d2,t1_2),e2,f2,g2]

        obs = flow_observe(w1, w2)

        print(f"\n  Trial {trial}:")
        print(f"    a-merged rounds: {obs['a_merged'][:20]}")
        print(f"    e-merged rounds: {obs['e_merged'][:20]}")
        print(f"    full-merged:     {obs['full_merged'][:20]}")
        print(f"    sub-stream hash match: {obs['sub_match']}")
        print(f"    matched sub-streams: {obs['n_sub_match']}/8")

        if obs['n_sub_match'] > 0:
            print(f"    ★ {obs['n_sub_match']} sub-streams MERGED at output!")
            # Which ones?
            h1=sha(w1); h2=sha(w2)
            for i in range(8):
                if h1[i]==h2[i]:
                    print(f"      H[{i}] = 0x{h1[i]:08x} (MATCH!)")


if __name__ == '__main__':
    flow_merger()
