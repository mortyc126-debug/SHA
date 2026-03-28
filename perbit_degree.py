#!/usr/bin/env python3
"""
Per-Bit Algebraic Degree — Axis 4.2 Deep Dive
=================================================

N3 found: H[7][b0] degree = 9/10 at 64 rounds (10-bit input).
Question: is degree < n REAL, or artifact of small input?

If degree < n persists at n=16,20 → algebraic weakness in specific bits.
If degree = n at larger n → small-n artifact.

Also: which OUTPUT BITS have LOWEST degree? These are attack targets.

Method: ANF via Möbius transform on truth tables.
Limitation: 2^n evaluations needed. n=16 → 65536 evals (fast).
n=20 → 1M evals (~3 min). n=24 → 16M (too slow in Python).
"""

import numpy as np
import time, sys

MASK = 0xFFFFFFFF
H0 = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
      0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]
K = [0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
     0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
     0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
     0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
     0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
     0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
     0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
     0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2]

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&MASK
def sig0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def sig1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def Sig0(x): return rotr(x,2)^rotr(x,13)^rotr(x,22)
def Sig1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def Ch(e,f,g): return (e&f)^(~e&g)&MASK
def Maj(a,b,c): return (a&b)^(a&c)^(b&c)
def hw(x): return bin(x).count('1')

def msg_sched(W16):
    W=list(W16)+[0]*48
    for i in range(16,64): W[i]=(sig1(W[i-2])+W[i-7]+sig0(W[i-15])+W[i-16])&MASK
    return W

def sha256_reduced(W16, num_rounds):
    W=msg_sched(W16); a,b,c,d,e,f,g,h=H0
    for r in range(num_rounds):
        T1=(h+Sig1(e)+Ch(e,f,g)+K[r]+W[r])&MASK; T2=(Sig0(a)+Maj(a,b,c))&MASK
        h,g,f=g,f,e; e=(d+T1)&MASK; d,c,b=c,b,a; a=(T1+T2)&MASK
    return [(v+iv)&MASK for v,iv in zip([a,b,c,d,e,f,g,h],H0)]


def compute_degree(n_bits, num_rounds, target_word, target_bit, base_W=None):
    """
    Compute algebraic degree of H[target_word][target_bit]
    as function of W[0][0..n_bits-1], with W[0][n_bits..31]=0, W[1..15]=0.

    If base_W provided, use it as base (non-zero higher bits).
    """
    N = 1 << n_bits
    tt = np.zeros(N, dtype=np.int8)

    for x in range(N):
        if base_W is not None:
            W16 = list(base_W)
            W16[0] = (W16[0] & ~((1 << n_bits) - 1)) | x
        else:
            W16 = [x] + [0]*15
        H = sha256_reduced(W16, num_rounds)
        tt[x] = (H[target_word] >> target_bit) & 1

    # Möbius transform → ANF
    anf = tt.copy()
    for i in range(n_bits):
        step = 1 << i
        for j in range(0, N, step << 1):
            for k in range(step):
                anf[j + k + step] ^= anf[j + k]

    # Degree = max HW of index with anf[idx]=1
    max_deg = 0
    for idx in range(N):
        if anf[idx] == 1:
            d = hw(idx)
            if d > max_deg:
                max_deg = d

    # Count monomials per degree
    mono_counts = {}
    for idx in range(N):
        if anf[idx] == 1:
            d = hw(idx)
            mono_counts[d] = mono_counts.get(d, 0) + 1

    return max_deg, mono_counts


def experiment_degree_scan(n_bits=14, seed=13000):
    print("="*70)
    print(f"PER-BIT DEGREE SCAN — n={n_bits} input bits")
    print(f"Truth table size: 2^{n_bits} = {1<<n_bits}")
    print("="*70)

    # Scan all 8 output words × selected bits × selected rounds
    target_bits = [0, 7, 15, 23, 31]  # 5 bits per word
    rounds_to_test = [5, 8, 16, 64]

    print(f"\n  {'rounds':>6} | {'H[w][b]':>10} | {'degree':>8} | {'monomials':>10} | {'note':>20}")
    print(f"  {'-'*6}-+-{'-'*10}-+-{'-'*8}-+-{'-'*10}-+-{'-'*20}")

    t0 = time.time()
    results = []

    for num_r in rounds_to_test:
        for w in [0, 4, 5, 6, 7]:
            for b in target_bits:
                deg, monos = compute_degree(n_bits, num_r, w, b)
                total_monos = sum(monos.values())
                note = ""
                if deg < n_bits:
                    note = f"★ BELOW MAX ({deg}<{n_bits})"
                elif deg == n_bits and total_monos < 10:
                    note = "sparse"

                results.append((num_r, w, b, deg, total_monos, note))
                if note or num_r in [5, 64]:
                    print(f"  {num_r:6d} | H[{w}][b{b:2d}] | {deg:4d}/{n_bits} | {total_monos:10d} | {note}")

        elapsed = time.time() - t0
        print(f"  ... {num_r} rounds done ({elapsed:.1f}s)")

    # Summary: which bits have degree < n at 64 rounds?
    print(f"\n  === SUMMARY: Bits with degree < {n_bits} at 64 rounds ===")
    low_deg = [(w, b, d, m) for r, w, b, d, m, _ in results if r == 64 and d < n_bits]
    if low_deg:
        for w, b, d, m in low_deg:
            print(f"    H[{w}][b{b:2d}]: degree = {d}/{n_bits} ★ (monomials: {m})")
        print(f"\n    ★★★ {len(low_deg)} bits have degree BELOW maximum!")
        print(f"    These bits are algebraically simpler than expected.")
        print(f"    Potential linearization attack on these specific bits.")
    else:
        print(f"    All tested bits have degree = {n_bits} (full)")
        print(f"    No algebraic weakness detected at n={n_bits}")

    # Summary at 5 rounds
    print(f"\n  === Bits with degree < {n_bits} at 5 rounds ===")
    low_5 = [(w, b, d, m) for r, w, b, d, m, _ in results if r == 5 and d < n_bits]
    for w, b, d, m in low_5:
        print(f"    H[{w}][b{b:2d}]: degree = {d}/{n_bits} (monomials: {m})")

    return results


def experiment_degree_vs_n(target_word=7, target_bit=0, seed=13001):
    """How does degree scale with n for a fixed output bit?"""
    print(f"\n{'='*70}")
    print(f"DEGREE vs n — H[{target_word}][b{target_bit}] at 64 rounds")
    print("="*70)

    print(f"\n  {'n':>3} | {'degree':>7} | {'deg/n':>6} | {'monomials':>10} | {'time':>6}")
    print(f"  {'-'*3}-+-{'-'*7}-+-{'-'*6}-+-{'-'*10}-+-{'-'*6}")

    for n in [8, 10, 12, 14, 16]:
        t0 = time.time()
        deg, monos = compute_degree(n, 64, target_word, target_bit)
        elapsed = time.time() - t0
        total_monos = sum(monos.values())
        ratio = deg / n

        marker = " ★" if deg < n else ""
        print(f"  {n:3d} | {deg:4d}/{n:<2d} | {ratio:6.2f} | {total_monos:10d} | {elapsed:6.1f}s{marker}")

    # Same for H[0][b0] as comparison
    print(f"\n  Comparison: H[0][b0] at 64 rounds:")
    for n in [8, 10, 12, 14]:
        t0 = time.time()
        deg, monos = compute_degree(n, 64, 0, 0)
        elapsed = time.time() - t0
        print(f"  n={n:2d}: degree={deg}/{n} ({time.time()-t0:.1f}s)")


def experiment_degree_randomized(n_bits=12, num_rounds=64, n_bases=5, seed=13002):
    """Test degree with RANDOM base message (not just W[0][0..n-1])."""
    print(f"\n{'='*70}")
    print(f"RANDOMIZED DEGREE — n={n_bits}, {num_rounds} rounds, {n_bases} random bases")
    print("="*70)

    rng = np.random.RandomState(seed)

    print(f"\n  Does degree change with different base messages?")
    print(f"  (n=10 test used W[1..15]=0. Real messages have W[1..15]≠0)")

    for target_w, target_b in [(7, 0), (7, 31), (0, 0), (6, 31)]:
        degrees = []
        for base_idx in range(n_bases):
            base_W = [int(rng.randint(0, 1<<32)) for _ in range(16)]
            deg, _ = compute_degree(n_bits, num_rounds, target_w, target_b, base_W=base_W)
            degrees.append(deg)

        min_d = min(degrees)
        max_d = max(degrees)
        marker = " ★★" if min_d < n_bits else ""
        print(f"  H[{target_w}][b{target_b:2d}]: degrees = {degrees}, range [{min_d},{max_d}]{marker}")

    return


# ============================================================
# MAIN
# ============================================================

if __name__ == '__main__':
    print("Per-Bit Algebraic Degree — Axis 4.2 Deep Dive")
    print(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")

    n_bits = 14
    if len(sys.argv) > 1 and sys.argv[1] == '--fast':
        n_bits = 12

    t_start = time.time()
    results = experiment_degree_scan(n_bits=n_bits)
    experiment_degree_vs_n(target_word=7, target_bit=0)
    experiment_degree_randomized(n_bits=min(n_bits, 12))

    total = time.time() - t_start
    print(f"\n{'='*70}")
    print(f"TOTAL TIME: {total:.1f}s")
    print(f"{'='*70}")
