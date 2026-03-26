#!/usr/bin/env python3
"""
14 раундов BOTH=0. Ломается на r=17.
Причина: W₂[16] ≠ W₁[16] (schedule из разных W[0..15]).

Вопрос: НАСКОЛЬКО W₂[16] отличается от W₁[16]?
Если отличие маленькое → может продлить.

И главное: δW[16] = schedule(W₂)[16] - schedule(W₁)[16].
schedule[16] = σ1(W[14]) + W[9] + σ0(W[1]) + W[0].

Все компоненты (W[0], W[1], W[9], W[14]) — из a-repair.
Мы ЗНАЕМ δW[0] (break), δW[1], δW[9], δW[14] (a-repair).
Значит δW[16] = σ1(δW[14]) + δW[9] + σ0(δW[1]) + δW[0]
(с поправкой на carry).

Это ВЫЧИСЛИМО. И если δW[16] = 0 → раунд 17 тоже BOTH=0!
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


def a_repair_and_measure(w1, break_delta=1):
    """a-repair с break_delta, вернуть W₂ и δschedule."""

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
    # Round 0: break
    a2,b2,c2,d2,e2,f2,g2,h2 = s2
    t1_2=A(h2,S1(e2),ch(e2,f2,g2),C[0],w2[0])
    t2_2=A(S0(a2),mj(a2,b2,c2))
    s2=[A(t1_2,t2_2),a2,b2,c2,A(d2,t1_2),e2,f2,g2]

    # Rounds 1..15: a-repair
    for r in range(1, 16):
        a2,b2,c2,d2,e2,f2,g2,h2 = s2
        a_target = s1_all[r+1][0]
        t2_2 = A(S0(a2),mj(a2,b2,c2))
        t1_need = D(a_target, t2_2)
        w2[r] = D(D(D(D(t1_need,h2),S1(e2)),ch(e2,f2,g2)),C[r])
        t1_2=A(h2,S1(e2),ch(e2,f2,g2),C[r],w2[r])
        s2=[A(t1_2,t2_2),a2,b2,c2,A(d2,t1_2),e2,f2,g2]

    # Schedule δ для r=16..63
    e1_full = expand(w1)
    e2_full = expand(w2)

    return w1, w2, e1_full, e2_full


def measure_extension(N):
    print("="*60)
    print("EXTENSION: δschedule[16..63] после a-repair")
    print("="*60)

    for trial in range(N):
        w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

        # Разные break_delta
        for delta in [1, 2, 7, 0x100, 0x10000, 0x80000000]:
            _, w2, e1, e2 = a_repair_and_measure(w1, break_delta=delta)

            # δW для первых нескольких schedule rounds
            dw = [HW(e1[r] ^ e2[r]) for r in range(64)]

            # Сколько δW[r]=0 подряд с r=16?
            consec_zero = 0
            for r in range(16, 64):
                if e1[r] == e2[r]:
                    consec_zero += 1
                else:
                    break

            # Полный hash
            h1=sha(w1); h2=sha(w2)
            dh=sum(HW(h1[i]^h2[i]) for i in range(8))
            ndiff = sum(1 for i in range(16) if w1[i]!=w2[i])

            if trial == 0:
                print(f"\n  δ=0x{delta:08x}: schedule zeros from r=16: {consec_zero}, "
                      f"dH={dh}, words_diff={ndiff}")
                # Show δW[16..25]
                print(f"    δW: ", end="")
                for r in range(16, min(28, 64)):
                    print(f"r{r}={dw[r]:2d} ", end="")
                print()


def search_for_zero_dw16(N):
    """Ищем break_delta где δW[16]=0.
    δW[16] = schedule(W₂)[16] ⊕ schedule(W₁)[16].
    Зависит от δW[0], δW[1], δW[9], δW[14]."""

    print(f"\n{'='*60}")
    print("SEARCH: break_delta where δW[16]=0")
    print("  If found → r=17 also BOTH=0 → chain extends!")
    print("="*60)

    w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

    best_dw16 = 32
    best_delta = 0
    best_consec = 0

    for delta in range(1, N+1):
        _, w2, e1, e2 = a_repair_and_measure(w1, break_delta=delta)
        dw16 = HW(e1[16] ^ e2[16])

        consec = 0
        for r in range(16, 64):
            if e1[r] == e2[r]: consec += 1
            else: break

        if dw16 < best_dw16:
            best_dw16 = dw16
            best_delta = delta
            best_consec = consec
            if dw16 <= 3:
                h1=sha(w1); h2=sha(w2)
                dh=sum(HW(h1[i]^h2[i]) for i in range(8))
                print(f"  δ={delta:8d}: δW[16]={dw16} bits, consec_zeros={consec}, dH={dh}")

        if dw16 == 0:
            print(f"  ★★★ δW[16]=0! Chain extends!")
            # Check further
            for r in range(16, 30):
                dwr = HW(e1[r]^e2[r])
                print(f"    δW[{r}]={dwr}", end="")
                if dwr == 0: print(" ✓", end="")
                print()
            h1=sha(w1); h2=sha(w2)
            dh=sum(HW(h1[i]^h2[i]) for i in range(8))
            print(f"    dH={dh}")
            break

    print(f"\n  Best: δ={best_delta}, δW[16]={best_dw16} bits, consec={best_consec}")
    print(f"  Expected: ~16 bits (random)")
    print(f"  If best < 10 → a-repair partially controls schedule!")

    # Birthday estimate for δW[16]=0
    p = 2**(-32)  # P(δW[16]=0) ≈ 2^{-32} if random
    print(f"\n  P(δW[16]=0) ≈ 2^{{-32}} → need ~{1<<32:.0e} trials")


if __name__ == '__main__':
    import sys
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 100000

    measure_extension(1)
    search_for_zero_dw16(N)
