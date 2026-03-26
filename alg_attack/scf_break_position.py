#!/usr/bin/env python3
"""
BREAK POSITION: carry-free break (MSB) vs carry-heavy break (LSB).

Flow prediction: MSB break → no carry → faster convergence → more merged rounds.

Measure: для каждой bit position break:
  - convergence speed (rounds to first BOTH=0)
  - merged duration (how many consecutive BOTH=0)
  - total BOTH=0 rounds
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

def a_repair_with_bit(w1, break_word, break_bit):
    """a-repair с break = flip одного бита W[break_word][break_bit]."""
    e1 = expand(w1)
    s1a = [list(IV)]; s = list(IV)
    for r in range(64):
        a,b,c,d,ee,f,g,h = s
        t1=A(h,S1(ee),ch(ee,f,g),C[r],e1[r])
        t2=A(S0(a),mj(a,b,c))
        s=[A(t1,t2),a,b,c,A(d,t1),ee,f,g]
        s1a.append(list(s))

    w2 = list(w1)
    w2[break_word] ^= (1 << break_bit)

    s2 = list(s1a[break_word])
    a2,b2,c2,d2,e2,f2,g2,h2 = s2
    t1_2=A(h2,S1(e2),ch(e2,f2,g2),C[break_word],w2[break_word])
    t2_2=A(S0(a2),mj(a2,b2,c2))
    s2=[A(t1_2,t2_2),a2,b2,c2,A(d2,t1_2),e2,f2,g2]

    for r in range(break_word+1, 16):
        a2,b2,c2,d2,e2,f2,g2,h2 = s2
        t2_2 = A(S0(a2),mj(a2,b2,c2))
        t1n = D(s1a[r+1][0], t2_2)
        w2[r] = D(D(D(D(t1n,h2),S1(e2)),ch(e2,f2,g2)),C[r])
        t1_2=A(h2,S1(e2),ch(e2,f2,g2),C[r],w2[r])
        s2=[A(t1_2,t2_2),a2,b2,c2,A(d2,t1_2),e2,f2,g2]

    return w2


def measure_convergence(w1, w2):
    """Измеряю: convergence speed и merge duration."""
    e1=expand(w1); e2=expand(w2)
    s1=list(IV); s2=list(IV)

    both0_rounds = []
    for r in range(64):
        a1,b1,c1,d1,ee1,f1,g1,h1=s1
        t1_1=A(h1,S1(ee1),ch(ee1,f1,g1),C[r],e1[r])
        t2_1=A(S0(a1),mj(a1,b1,c1))
        s1=[A(t1_1,t2_1),a1,b1,c1,A(d1,t1_1),ee1,f1,g1]
        a2,b2,c2,d2,ee2,f2,g2,h2=s2
        t1_2=A(h2,S1(ee2),ch(ee2,f2,g2),C[r],e2[r])
        t2_2=A(S0(a2),mj(a2,b2,c2))
        s2=[A(t1_2,t2_2),a2,b2,c2,A(d2,t1_2),ee2,f2,g2]
        if s1==s2: both0_rounds.append(r+1)

    return both0_rounds


def compare_bit_positions(N):
    print("="*60)
    print("BREAK BIT POSITION: carry-free (MSB) vs carry-heavy (LSB)")
    print("="*60)

    # Break word = 3 (our optimal position)
    break_word = 3

    results = {}

    for bit_pos in range(32):
        total_both0 = 0
        total_last = 0
        total_first_merge = 0

        for _ in range(N):
            w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
            w2 = a_repair_with_bit(w1, break_word, bit_pos)
            b0 = measure_convergence(w1, w2)

            total_both0 += len(b0)

            # First merge AFTER break (not counting before-break rounds)
            after_break = [r for r in b0 if r > break_word + 1]
            if after_break:
                total_first_merge += after_break[0]
            else:
                total_first_merge += 64

            total_last += (max(b0) if b0 else 0)

        avg_both0 = total_both0 / N
        avg_first = total_first_merge / N
        avg_last = total_last / N
        convergence = avg_first - break_word - 1  # rounds to converge

        results[bit_pos] = (avg_both0, convergence, avg_last)

    # Sort by both0 (more = better)
    print(f"\n  break_word={break_word}, N={N}")
    print(f"  {'bit':>4} | {'both0':>6} {'convergence':>12} {'last_b0':>8} | carrier")
    print("  " + "-"*55)

    for bit in sorted(results.keys(), key=lambda b: -results[b][0]):
        b0, conv, last = results[bit]
        carrier = "LSB" if bit == 0 else ("MSB" if bit == 31 else f"mid({bit})")
        star = ""
        if b0 > 9: star = " ★"
        if b0 > 10: star = " ★★"
        print(f"  {bit:4d} | {b0:6.1f} {conv:12.1f} {last:8.1f} | {carrier}{star}")

    # Summary: LSB vs MSB
    lsb = results[0]; msb = results[31]
    print(f"\n  LSB (bit 0): both0={lsb[0]:.1f}, convergence={lsb[1]:.1f}")
    print(f"  MSB (bit 31): both0={msb[0]:.1f}, convergence={msb[1]:.1f}")
    print(f"  MSB advantage: {msb[0]-lsb[0]:+.1f} both0, {lsb[1]-msb[1]:+.1f} faster convergence")


if __name__ == '__main__':
    import sys
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 50
    compare_bit_positions(N)
