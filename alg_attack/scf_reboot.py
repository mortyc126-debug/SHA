#!/usr/bin/env python3
"""
REBOOT: после BOTH=0 на r=18, δW[18]≈4 бит = МАЛЕНЬКИЙ автоbreak.

Вопрос: маленький δW[18] → маленький δstate[19] →
может ли schedule СЛУЧАЙНО продолжить BOTH=0?

Если δstate[19] маленький (4 бита) И δW[19] СЛУЧАЙНО ≈ нужный
→ BOTH=0 на r=20 → НОВЫЙ cascade?

P(BOTH=0 на r=19 | BOTH=0 на r=18, δW[18]=4 bits):
  δstate[19] = step(same_state, W+4bits) - step(same_state, W)
  ≈ 4 бита разницы (маленький break)
  Для BOTH=0 на r=19: δW[19] должен КОМПЕНСИРОВАТЬ δstate[19].
  P ≈ 2^{-4}? (4 бита нужно совпасть)

Если P ≈ 2^{-4} на каждый "reboot" → за O(2^4)=16 проб W₁
можно продлить цепочку на 1 раунд.

За K раундов: O(2^{4K}). Для 46 раундов gap: O(2^{184}).
Хуже birthday. НО: δW может быть < 4 бит для правильных W₁!
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

def step(state, w, k):
    a,b,c,d,e,f,g,h = state
    t1 = A(h,S1(e),ch(e,f,g),k,w)
    t2 = A(S0(a),mj(a,b,c))
    return [A(t1,t2),a,b,c,A(d,t1),e,f,g]

def sha(w):
    e=expand(w); s=list(IV)
    for r in range(64):
        s = step(s, e[r], C[r])
    return tuple(A(IV[i],s[i]) for i in range(8))


def do_a_repair_r3(w1, break_delta=1):
    e1 = expand(w1)
    s1_all = [list(IV)]
    s = list(IV)
    for r in range(64):
        s = step(s, e1[r], C[r])
        s1_all.append(list(s))

    w2 = list(w1); w2[3] = A(w1[3], break_delta)
    s2 = list(IV)
    for r in range(3):
        s2 = step(s2, e1[r], C[r])
    s2 = step(s2, w2[3], C[3])
    for r in range(4, 16):
        a2,b2,c2,d2,e2,f2,g2,h2 = s2
        a_target = s1_all[r+1][0]
        t2_2 = A(S0(a2),mj(a2,b2,c2))
        t1_need = D(a_target, t2_2)
        w2[r] = D(D(D(D(t1_need,h2),S1(e2)),ch(e2,f2,g2)),C[r])
        s2 = step(s2, w2[r], C[r])
    return w2


def measure_reboot(N):
    """После a-repair, смотрю что происходит ПОСЛЕ r=18."""
    print("="*60)
    print("REBOOT: что происходит после BOTH=0 обрывается?")
    print("="*60)

    for trial in range(N):
        w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        w2 = do_a_repair_r3(w1)
        e1 = expand(w1); e2 = expand(w2)

        # Track BOTH=0 AND δstate size AND δW size
        s1 = list(IV); s2 = list(IV)

        print(f"\n  Trial {trial}:")
        print(f"  {'r':>3} | {'same':>4} | {'δa':>3} {'δe':>3} {'δ_tot':>5} | {'δW':>3} | notes")
        print("  " + "-"*50)

        prev_same = True
        reboots = 0

        for r in range(64):
            s1 = step(s1, e1[r], C[r])
            s2 = step(s2, e2[r], C[r])

            same = (s1 == s2)
            da = HW(s1[0]^s2[0])
            de = HW(s1[4]^s2[4])
            dtot = sum(HW(s1[i]^s2[i]) for i in range(8))
            dw = HW(e1[r]^e2[r])

            notes = ""
            if same: notes = "BOTH=0"
            elif dtot <= 10: notes = f"NEAR (δ={dtot})"

            if prev_same and not same:
                notes += " ← BREAK"
                reboots += 1
            if not prev_same and same:
                notes += " ← REBOOT!"

            if r < 25 or same or dtot <= 15 or (not prev_same and same):
                print(f"  {r+1:3d} | {'Y' if same else 'N':>4} | {da:3d} {de:3d} {dtot:5d} | {dw:3d} | {notes}")

            prev_same = same

        h1=sha(w1); h2=sha(w2)
        dh=sum(HW(h1[i]^h2[i]) for i in range(8))
        print(f"\n  dH={dh}, reboots={reboots}")


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    measure_reboot(N)
