#!/usr/bin/env python3
"""
Carry как ПОТОК. Не число. Поток.

Carry течёт от bit 0 к bit 31 в КАЖДОМ сложении.
Между сложениями: rotate перемешивает позиции.
Результат: carry-путь ИЗВИВАЕТСЯ через 64 раунда.

Я прослежу: для конкретного δW (flip одного бита),
ГДЕ carry возникает, КУДА течёт, и ЧЕГО КАСАЕТСЯ.

Не "сколько бит carry" а "КАКОЙ ПУТЬ carry проходит".
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


def trace_carry(w1, flip_word, flip_bit):
    """Прослеживаю carry-поток для одного bit-flip.

    Запускаю SHA с W1 и W1⊕flip. На КАЖДОМ сложении:
    смотрю δ-carry (разница carry между двумя потоками).

    carry(a,b) = ((a+b) ^ (a^b)) >> 1

    Для каждого раунда и каждого сложения:
    записываю δ-carry (32-бит маска: какие позиции carry отличаются).
    """

    e1 = expand(w1)
    w2 = list(w1); w2[flip_word] ^= (1 << flip_bit)
    e2 = expand(w2)

    s1 = list(IV); s2 = list(IV)

    carry_trace = []

    for r in range(64):
        a1,b1,c1,d1,ee1,f1,g1,h1 = s1
        a2,b2,c2,d2,ee2,f2,g2,h2 = s2

        # T1 = h + S1(e) + Ch(e,f,g) + K + W
        # Trace carry at each addition step
        round_carries = []

        # Step 1: h + S1(e)
        p1_1 = A(h1, S1(ee1)); p1_2 = A(h2, S1(ee2))
        c1_1 = (p1_1 ^ (h1 ^ S1(ee1))); c1_2 = (p1_2 ^ (h2 ^ S1(ee2)))
        dc1 = c1_1 ^ c1_2
        round_carries.append(('h+S1', HW(dc1), dc1))

        # Step 2: + Ch
        q1_1 = A(p1_1, ch(ee1,f1,g1)); q1_2 = A(p1_2, ch(ee2,f2,g2))
        dc2 = (q1_1 ^ (p1_1 ^ ch(ee1,f1,g1))) ^ (q1_2 ^ (p1_2 ^ ch(ee2,f2,g2)))
        round_carries.append(('...+Ch', HW(dc2), dc2))

        # Step 3: + K
        r1_1 = A(q1_1, C[r]); r1_2 = A(q1_2, C[r])
        dc3 = (r1_1 ^ (q1_1 ^ C[r])) ^ (r1_2 ^ (q1_2 ^ C[r]))
        round_carries.append(('...+K', HW(dc3), dc3))

        # Step 4: + W
        t1_1 = A(r1_1, e1[r]); t1_2 = A(r1_2, e2[r])
        dc4 = (t1_1 ^ (r1_1 ^ e1[r])) ^ (t1_2 ^ (r1_2 ^ e2[r]))
        round_carries.append(('...+W=T1', HW(dc4), dc4))

        # T2 = S0(a) + Maj
        t2_1 = A(S0(a1), mj(a1,b1,c1)); t2_2 = A(S0(a2), mj(a2,b2,c2))
        dc5 = (t2_1 ^ (S0(a1)^mj(a1,b1,c1))) ^ (t2_2 ^ (S0(a2)^mj(a2,b2,c2)))
        round_carries.append(('S0+Maj=T2', HW(dc5), dc5))

        # a_new = T1 + T2
        an1 = A(t1_1,t2_1); an2 = A(t1_2,t2_2)
        dc6 = (an1 ^ (t1_1^t2_1)) ^ (an2 ^ (t1_2^t2_2))
        round_carries.append(('T1+T2=a', HW(dc6), dc6))

        # e_new = d + T1
        en1 = A(d1,t1_1); en2 = A(d2,t1_2)
        dc7 = (en1 ^ (d1^t1_1)) ^ (en2 ^ (d2^t1_2))
        round_carries.append(('d+T1=e', HW(dc7), dc7))

        # Total δ-carry this round
        total_dc = sum(c[1] for c in round_carries)
        carry_trace.append((r, total_dc, round_carries))

        s1 = [an1,a1,b1,c1,en1,ee1,f1,g1]
        s2 = [an2,a2,b2,c2,en2,ee2,f2,g2]

    return carry_trace


def look_at_flow():
    print("="*60)
    print("CARRY FLOW: tracing δ-carry through SHA")
    print("="*60)

    w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

    # Flip bit 0 of W[3] (our a-repair break position)
    trace = trace_carry(w1, 3, 0)

    print(f"\n  Flip: W[3] bit 0")
    print(f"  {'r':>3} | {'total_δcarry':>12} | per-addition breakdown")
    print("  " + "-"*70)

    for r, total, adds in trace:
        if r < 20 or total < 5 or r > 58:
            breakdown = " ".join(f"{name}:{hw}" for name, hw, _ in adds if hw > 0)
            mark = ""
            if total == 0: mark = " ← NO CARRY DIFF"
            elif total <= 3: mark = " ← tiny"
            print(f"  {r:3d} | {total:12d} | {breakdown}{mark}")

    # Summarize: where is carry flow concentrated?
    print(f"\n  Carry flow by phase:")
    phases = [(0,4,"r=0..3"), (4,12,"r=4..11"), (12,18,"r=12..17"),
              (18,32,"r=18..31"), (32,48,"r=32..47"), (48,64,"r=48..63")]
    for start, end, name in phases:
        total = sum(trace[r][1] for r in range(start, end))
        avg = total / (end - start)
        print(f"    {name}: total δ-carry = {total}, avg = {avg:.1f}/round")

    # KEY: does carry flow STOP or CONTINUE?
    # If it stops after round K → no carry influence after K → can exploit!
    zero_carry_rounds = [r for r, total, _ in trace if total == 0]
    print(f"\n  Zero δ-carry rounds: {zero_carry_rounds[:20]}")
    print(f"  Total: {len(zero_carry_rounds)}/64")

    # The DIRECTION of carry: which bit positions are affected?
    print(f"\n  δ-carry bit position heatmap (first 10 rounds):")
    for r in range(min(10, len(trace))):
        _, _, adds = trace[r]
        all_dc = 0
        for _, _, dc_val in adds:
            all_dc |= dc_val
        # Show which bits
        bits_affected = [b for b in range(32) if (all_dc>>b)&1]
        if bits_affected:
            print(f"    r={r}: bits {bits_affected[:10]}{'...' if len(bits_affected)>10 else ''}")
        else:
            print(f"    r={r}: NO carry diff")


if __name__ == '__main__':
    look_at_flow()
