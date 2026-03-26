#!/usr/bin/env python3
"""
FLOW MATHEMATICS — математика потоков.

Поток = δ-данные перемещающиеся через SHA.
Не значения — МАРШРУТЫ.

Каждый поток имеет:
  - Birth: (round, bit) — откуда начался
  - Width: сколько бит-позиций затронут на каждом раунде
  - Direction: carry идёт LSB→MSB, rotate перенаправляет

Collision = ВСЕ потоки скомпенсированы к раунду 64.

Мои уравнения потоков:
  W(r, b) = width потока начавшегося в (r, b) через K раундов
  W(r, b, 0) = 1 (рождение = 1 бит)
  W(r, b, 1) = ? (через 1 раунд)
  W(r, b, K) = ? (через K раундов)

Измерю W(r, b, K) для всех (r, b, K).
Это ХАРАКТЕРИСТИКА SHA как потоковой машины.
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


def measure_flow_width(w1, flip_word, flip_bit):
    """Ширина потока: сколько state-бит отличаются на каждом раунде."""
    e1 = expand(w1)
    w2 = list(w1); w2[flip_word] ^= (1 << flip_bit)
    e2 = expand(w2)

    s1 = list(IV); s2 = list(IV)
    widths = []  # (round, width_a, width_e, total_state_width)

    for r in range(64):
        a1,b1,c1,d1,ee1,f1,g1,h1 = s1
        a2,b2,c2,d2,ee2,f2,g2,h2 = s2

        t1_1=A(h1,S1(ee1),ch(ee1,f1,g1),C[r],e1[r])
        t2_1=A(S0(a1),mj(a1,b1,c1))
        s1=[A(t1_1,t2_1),a1,b1,c1,A(d1,t1_1),ee1,f1,g1]

        t1_2=A(h2,S1(ee2),ch(ee2,f2,g2),C[r],e2[r])
        t2_2=A(S0(a2),mj(a2,b2,c2))
        s2=[A(t1_2,t2_2),a2,b2,c2,A(d2,t1_2),ee2,f2,g2]

        wa = HW(s1[0]^s2[0])  # δa width
        we = HW(s1[4]^s2[4])  # δe width
        total = sum(HW(s1[i]^s2[i]) for i in range(8))

        widths.append((r+1, wa, we, total))

    return widths


def flow_equations():
    """Мои уравнения потоков."""
    print("="*60)
    print("FLOW EQUATIONS")
    print("="*60)

    # Измерю W(word, bit, K) для разных starting points
    w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

    # Flow from different bit positions
    print(f"\n  Flow width from different birth points:")
    print(f"  {'birth':>12} | K=1  K=2  K=3  K=4  K=5  K=8  K=16 K=32 K=64")
    print("  " + "-"*75)

    for word in [0, 3, 7, 15]:
        for bit in [0, 15, 31]:
            widths = measure_flow_width(w1, word, bit)
            # Extract total width at specific K values
            vals = []
            for K in [1, 2, 3, 4, 5, 8, 16, 32, 64]:
                if K <= len(widths):
                    vals.append(widths[K-1][3])
                else:
                    vals.append(0)

            # Flow RATE = width[K+1] / width[K]
            rates = []
            for i in range(min(5, len(widths)-1)):
                if widths[i][3] > 0:
                    rates.append(widths[i+1][3] / widths[i][3])
                else:
                    rates.append(0)

            print(f"  W[{word:2d}][{bit:2d}]  | " +
                  " ".join(f"{v:4d}" for v in vals))

    # Flow RATE universal?
    print(f"\n  Flow growth rates (width[K+1]/width[K]):")
    rates_all = {k: [] for k in range(10)}

    for _ in range(20):
        w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        word = int.from_bytes(os.urandom(1),'big') % 16
        bit = int.from_bytes(os.urandom(1),'big') % 32
        widths = measure_flow_width(w1, word, bit)
        for k in range(min(10, len(widths)-1)):
            if widths[k][3] > 0:
                rates_all[k].append(widths[k+1][3] / widths[k][3])

    for k in range(8):
        if rates_all[k]:
            avg_rate = sum(rates_all[k]) / len(rates_all[k])
            mn = min(rates_all[k]); mx = max(rates_all[k])
            print(f"    K={k}→{k+1}: avg rate = {avg_rate:.2f} (range {mn:.1f}-{mx:.1f})")

    # FLOW EQUATION:
    # width(K) = min(256, width(0) × R^K)
    # where R = growth rate ≈ ? (from data)

    print(f"\n  Flow equation: width(K) = min(256, 1 × R^K)")
    print(f"  where R = average growth rate per round")
    print(f"  Saturation at width = 256 (all state bits affected)")
    print(f"  Saturation round = log_R(256) = ?")

    if rates_all[0]:
        R_avg = sum(rates_all[0])/len(rates_all[0])
        import math
        if R_avg > 1:
            sat_round = math.log(256) / math.log(R_avg)
            print(f"  R = {R_avg:.2f}")
            print(f"  Saturation round = {sat_round:.1f}")
            print(f"  After saturation: width = 256 (constant)")

    # COMPENSATION equation:
    # Two flows F₁(birth₁) and F₂(birth₂) CANCEL if:
    # At every point where they overlap: F₁ + F₂ = 0 (mod 2)
    # Overlap happens when both flows reach the same (round, bit)

    print(f"\n  COMPENSATION:")
    print(f"  Two flows overlap where both affect the same state bits.")
    print(f"  If overlap is TOTAL (both affect same bits same way) → cancel!")
    print(f"  a-repair = anti-flow: cancels δ on a-register.")
    print(f"  Carry = flow leakage: part of δ escapes through carry chain.")
    print(f"  Carry leakage rate = how much flow escapes per round.")

    # Measure: carry leakage rate
    # = (actual δ) - (GF(2) predicted δ) = carry contribution
    # = flow that escapes the XOR-cancellation

    print(f"\n  Carry leakage measurement:")
    for trial in range(3):
        w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        word = 3; bit = 0
        widths = measure_flow_width(w1, word, bit)

        # GF(2) width (XOR-SHA)
        e1x = list(w1)
        for i in range(16,64):
            e1x.append(s1(e1x[i-2])^e1x[i-7]^s0(e1x[i-15])^e1x[i-16])
        w2x = list(w1); w2x[word] ^= (1<<bit)
        e2x = list(w2x)
        for i in range(16,64):
            e2x.append(s1(e2x[i-2])^e2x[i-7]^s0(e2x[i-15])^e2x[i-16])

        s1x=list(IV); s2x=list(IV)
        xor_widths = []
        for r in range(64):
            a1,b1,c1,d1,ee1,f1,g1,h1=s1x
            t1_1=h1^S1(ee1)^ch(ee1,f1,g1)^C[r]^e1x[r]
            t2_1=S0(a1)^mj(a1,b1,c1)
            s1x=[t1_1^t2_1,a1,b1,c1,d1^t1_1,ee1,f1,g1]
            a2,b2,c2,d2,ee2,f2,g2,h2=s2x
            t1_2=h2^S1(ee2)^ch(ee2,f2,g2)^C[r]^e2x[r]
            t2_2=S0(a2)^mj(a2,b2,c2)
            s2x=[t1_2^t2_2,a2,b2,c2,d2^t1_2,ee2,f2,g2]
            xw = sum(HW(s1x[i]^s2x[i]) for i in range(8))
            xor_widths.append(xw)

        if trial == 0:
            print(f"    {'K':>3} | {'real':>5} {'xor':>5} {'leak':>5}")
            print("    " + "-"*25)
            for K in [1,2,3,4,5,8,12,16,32,64]:
                real = widths[K-1][3]
                xor = xor_widths[K-1]
                leak = real - xor
                print(f"    {K:3d} | {real:5d} {xor:5d} {leak:5d}")


if __name__ == '__main__':
    flow_equations()
