#!/usr/bin/env python3
"""
Откуда 4 бита residual на δW[16]?

W[16] = σ1(W[14]) + W[9] + σ0(W[1]) + W[0]

δW[16] = δ(σ1(W[14])) + δW[9] + δ(σ0(W[1])) + δW[0] + carry_corrections

Каждое слагаемое: σ1 и σ0 — линейны над XOR.
Значит: δW[16]_xor = σ1(δW[14]) ⊕ δW[9] ⊕ σ0(δW[1]) ⊕ δW[0]
И: δW[16]_real = δW[16]_xor ⊕ carry_в_сумме

XOR-часть ВЫЧИСЛИМА. Carry-часть зависит от ЗНАЧЕНИЙ.

Если XOR-часть = 0 → δW[16] = carry only (≈ 16 бит).
Если XOR-часть ≠ 0 но мала → δW[16] мал.

Я посмотрю: чему равна XOR-часть при a-repair?
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
    """a-repair: вернуть w2 и все промежуточные данные."""
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

    return w1, w2


def decompose_dw16(N):
    print("="*60)
    print("DECOMPOSITION: δW[16] = σ1(δW[14]) ⊕ δW[9] ⊕ σ0(δW[1]) ⊕ δW[0]")
    print("                        + carry_correction")
    print("="*60)

    for trial in range(N):
        w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        _, w2 = a_repair(w1, break_delta=1)

        dw = [w1[i] ^ w2[i] for i in range(16)]

        # XOR-часть δW[16]
        dw16_xor = s1(dw[14]) ^ dw[9] ^ s0(dw[1]) ^ dw[0]

        # Реальная δW[16]
        e1 = expand(w1); e2 = expand(w2)
        dw16_real = e1[16] ^ e2[16]

        # Carry = разница
        carry = dw16_xor ^ dw16_real

        print(f"\n  Trial {trial}:")
        print(f"    δW[0]:  HW={HW(dw[0]):2d}  (break)")
        print(f"    δW[1]:  HW={HW(dw[1]):2d}  (a-repair round 1)")
        print(f"    δW[9]:  HW={HW(dw[9]):2d}  (a-repair round 9)")
        print(f"    δW[14]: HW={HW(dw[14]):2d}  (a-repair round 14)")
        print(f"    ---")
        print(f"    σ1(δW[14]): HW={HW(s1(dw[14])):2d}")
        print(f"    σ0(δW[1]):  HW={HW(s0(dw[1])):2d}")
        print(f"    ---")
        print(f"    δW[16]_xor:  HW={HW(dw16_xor):2d}  (XOR part)")
        print(f"    δW[16]_real: HW={HW(dw16_real):2d}  (real)")
        print(f"    carry:       HW={HW(carry):2d}")

        # Какие компоненты ДОМИНИРУЮТ?
        components = [
            ("σ1(δW[14])", s1(dw[14])),
            ("δW[9]", dw[9]),
            ("σ0(δW[1])", s0(dw[1])),
            ("δW[0]", dw[0]),
        ]

        # Попарная отмена?
        print(f"    ---")
        print(f"    Pairwise cancellation:")
        for i in range(4):
            for j in range(i+1, 4):
                cancel = HW(components[i][1] ^ components[j][1])
                name = f"{components[i][0]} ⊕ {components[j][0]}"
                if cancel < 16:
                    print(f"      {name}: HW={cancel} ★ (partial cancel!)")

        # Все 4 вместе
        all_xor = 0
        for _, v in components:
            all_xor ^= v
        print(f"    All 4 ⊕: HW={HW(all_xor)}")


def search_cancellation(N):
    """Ищем break_delta где XOR-часть δW[16] мала."""
    print(f"\n{'='*60}")
    print("SEARCH: break_delta where XOR-part of δW[16] is small")
    print("="*60)

    w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

    best_xor_hw = 32
    best_real_hw = 32
    best_delta = 0

    for delta in range(1, N+1):
        _, w2 = a_repair(w1, break_delta=delta)
        dw = [w1[i]^w2[i] for i in range(16)]

        dw16_xor = s1(dw[14]) ^ dw[9] ^ s0(dw[1]) ^ dw[0]
        e1=expand(w1); e2=expand(w2)
        dw16_real = e1[16]^e2[16]

        hw_xor = HW(dw16_xor)
        hw_real = HW(dw16_real)

        if hw_xor < best_xor_hw:
            best_xor_hw = hw_xor
            if hw_xor <= 2:
                print(f"  δ={delta}: XOR HW={hw_xor}, real HW={hw_real} ★")

        if hw_real < best_real_hw:
            best_real_hw = hw_real
            best_delta = delta

        if hw_real == 0:
            print(f"  ★★★ δW[16]=0 at δ={delta}!")
            # Check full chain
            for r in range(16, 30):
                print(f"    δW[{r}]={HW(e1[r]^e2[r])}")
            break

    print(f"\n  Best XOR: HW={best_xor_hw}")
    print(f"  Best real: HW={best_real_hw} at δ={best_delta}")


if __name__ == '__main__':
    decompose_dw16(3)
    search_cancellation(20000)
