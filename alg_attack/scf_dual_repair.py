#!/usr/bin/env python3
"""
DUAL REPAIR: чиним И δe И δa одновременно.

Wang чинит δe (через T1, управляемый W).
Wang НЕ чинит δa (через T2, автономный).

Новое: a-matching. Вместо T1₂=T1₁ → a₂[r+1]=a₁[r+1].
a_new = T1 + T2. Значит T1₂ = a₁[r+1] - T2₂.
W₂[r] = T1₂_needed - h₂ - S1(e₂) - Ch₂ - K[r]

Если a₁[r+1] = a₂[r+1]:
  → b₂[r+2] = a₂[r+1] = a₁[r+1] = b₁[r+2] ✓
  → через 3 шага: d₂ = a₂ = a₁ = d₁
  → T2 начинает совпадать!

a-matching создаёт КАСКАД: a→b→c→d совпадают через 3 раунда.
Потом T2 совпадает (зависит от a,b,c).
Потом δa_new = δT1 (управляемый) + δT2=0 = δT1.
И δe_new = δd + δT1. δd уже = 0 (через 3 раунда a-matching).
Значит δe_new = δT1. Управляемый!

ПОСЛЕ 3 раундов a-matching: ВСЁ управляемо через W!
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


def dual_repair(w1):
    """Строю W2 с a-matching: a₂[r+1] = a₁[r+1] для r=0..14.

    Это ДРУГАЯ стратегия чем Wang (который делает e-matching).
    a-matching → b,c,d совпадают через 3 раунда → T2 совпадает.
    """

    e1 = expand(w1)
    w2 = list(w1)
    w2[0] = A(w1[0], 1)  # Начальный break

    # Прогоняю поток 1 полностью, запоминаю все состояния
    s1_all = [list(IV)]
    s = list(IV)
    for r in range(64):
        a,b,c,d,ee,f,g,h = s
        t1=A(h,S1(ee),ch(ee,f,g),C[r],e1[r])
        t2=A(S0(a),mj(a,b,c))
        s=[A(t1,t2),a,b,c,A(d,t1),ee,f,g]
        s1_all.append(list(s))

    # Строю W2[0..15] для a-matching
    s2 = list(IV)
    for r in range(16):
        a2,b2,c2,d2,e2,f2,g2,h2 = s2

        # Цель: a₂[r+1] = a₁[r+1] = s1_all[r+1][0]
        a_target = s1_all[r+1][0]

        # a_new = T1 + T2
        t2_2 = A(S0(a2), mj(a2,b2,c2))
        t1_needed = D(a_target, t2_2)

        # W₂[r] = T1_needed - h₂ - S1(e₂) - Ch(e₂,f₂,g₂) - K[r]
        w2[r] = D(D(D(D(t1_needed, h2), S1(e2)), ch(e2,f2,g2)), C[r])

        # Шагаю поток 2
        t1_2 = A(h2,S1(e2),ch(e2,f2,g2),C[r],w2[r])
        s2 = [A(t1_2,t2_2),a2,b2,c2,A(d2,t1_2),e2,f2,g2]

    # Проверяю: a₂ совпадает с a₁ на каких раундах?
    e2_full = expand(w2)
    s2_check = list(IV)

    print(f"  {'r':>3} | {'δa':>4} {'δe':>4} {'δb':>4} {'δc':>4} {'δd':>4} {'δf':>4} {'δg':>4} {'δh':>4} | notes")
    print("  " + "-"*70)

    for r in range(22):
        a2,b2,c2,d2,e2,f2,g2,h2 = s2_check
        t1=A(h2,S1(e2),ch(e2,f2,g2),C[r],e2_full[r])
        t2=A(S0(a2),mj(a2,b2,c2))
        s2_check=[A(t1,t2),a2,b2,c2,A(d2,t1),e2,f2,g2]

        # Дельты
        da = HW(s2_check[0]^s1_all[r+1][0])
        de = HW(s2_check[4]^s1_all[r+1][4])
        db = HW(s2_check[1]^s1_all[r+1][1])
        dc = HW(s2_check[2]^s1_all[r+1][2])
        dd = HW(s2_check[3]^s1_all[r+1][3])
        df = HW(s2_check[5]^s1_all[r+1][5])
        dg = HW(s2_check[6]^s1_all[r+1][6])
        dh = HW(s2_check[7]^s1_all[r+1][7])

        notes = ""
        if da == 0: notes += " a=0!"
        if de == 0: notes += " e=0!"
        if da == 0 and de == 0: notes = " FULL MATCH!"

        free = r < 16
        print(f"  {r:3d} | {da:4d} {de:4d} {db:4d} {dc:4d} {dd:4d} {df:4d} {dg:4d} {dh:4d} | "
              f"{'[F]' if free else '[S]'}{notes}")

    # Hash
    h1 = sha(w1); h2 = sha(w2)
    dh = sum(HW(h1[i]^h2[i]) for i in range(8))
    n_diff = sum(1 for i in range(16) if w1[i]!=w2[i])
    print(f"\n  Hash diff: {dh}, words different: {n_diff}")

    return w2, dh


def wang_repair(w1):
    """Стандартный Wang (e-matching) для сравнения."""
    e1 = expand(w1)
    w2 = list(w1)
    w2[0] = A(w1[0], 1)

    s1_all = [list(IV)]
    s = list(IV)
    for r in range(64):
        a,b,c,d,ee,f,g,h = s
        t1=A(h,S1(ee),ch(ee,f,g),C[r],e1[r])
        t2=A(S0(a),mj(a,b,c))
        s=[A(t1,t2),a,b,c,A(d,t1),ee,f,g]
        s1_all.append(list(s))

    # Wang: e-matching (T1₂ = T1₁)
    s2 = list(IV)
    for r in range(16):
        a2,b2,c2,d2,e2,f2,g2,h2 = s2
        a1,b1,c1,d1,ee1,f1,g1,h1 = s1_all[r]

        t1_1 = A(h1,S1(ee1),ch(ee1,f1,g1),C[r],e1[r])
        w2[r] = D(D(D(D(t1_1,h2),S1(e2)),ch(e2,f2,g2)),C[r])

        t1_2 = t1_1  # forced
        t2_2 = A(S0(a2),mj(a2,b2,c2))
        s2 = [A(t1_2,t2_2),a2,b2,c2,A(d2,t1_2),e2,f2,g2]

    h1 = sha(w1); h2 = sha(w2)
    dh = sum(HW(h1[i]^h2[i]) for i in range(8))
    return w2, dh


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 3

    for trial in range(N):
        w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

        print(f"\n{'='*70}")
        print(f"Trial {trial}: DUAL REPAIR (a-matching)")
        print(f"{'='*70}")
        w2_dual, dh_dual = dual_repair(w1)

        print(f"\n--- WANG (e-matching) для сравнения ---")
        w2_wang, dh_wang = wang_repair(w1)
        print(f"  Wang hash diff: {dh_wang}")

        print(f"\n  DUAL: {dh_dual} vs WANG: {dh_wang}")
