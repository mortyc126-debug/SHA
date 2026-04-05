#!/usr/bin/env python3
"""
Word-level differential structure.

All previous tests: bit-level (XOR, HW).
This test: WORD-level (mod 2^32 arithmetic).

δa[r] = a₁[r] - a₂[r] mod 2^32 — the ADDITIVE difference.
Not XOR, not HW — the actual NUMBER.

Questions:
Q1: Is δa[r] uniformly distributed mod 2^32?
    Or does it cluster around specific values?
Q2: Is δa[r] related to δa[r-1] arithmetically?
    (not bit-correlation, but mod-arithmetic relation)
Q3: Does δa[r] mod small numbers have structure?
    (mod 8, mod 256, mod 2^k for various k)
Q4: Is the CARRY of δa[r] + δe[r] structured?
"""

import random, math
from collections import defaultdict

MASK32 = 0xFFFFFFFF
K=[0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2]
IV=[0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]
def rotr(x,n):return((x>>n)|(x<<(32-n)))&MASK32
def ssig0(x):return rotr(x,7)^rotr(x,18)^(x>>3)
def ssig1(x):return rotr(x,17)^rotr(x,19)^(x>>10)
def sig0(x):return rotr(x,2)^rotr(x,13)^rotr(x,22)
def sig1(x):return rotr(x,6)^rotr(x,11)^rotr(x,25)
def ch(e,f,g):return(e&f)^(~e&g)&MASK32
def maj(a,b,c):return(a&b)^(a&c)^(b&c)
def hw(x):return bin(x&MASK32).count('1')

def sha_states(M):
    W=list(M)+[0]*(64-len(M))
    for i in range(16,64):W[i]=(ssig1(W[i-2])+W[i-7]+ssig0(W[i-15])+W[i-16])&MASK32
    a,b,c,d,e,f,g,h=IV
    ss=[(a,b,c,d,e,f,g,h)]
    for r in range(64):
        T1=(h+sig1(e)+ch(e,f,g)+K[r]+W[r])&MASK32;T2=(sig0(a)+maj(a,b,c))&MASK32
        h,g,f=g,f,e;e=(d+T1)&MASK32;d,c,b=c,b,a;a=(T1+T2)&MASK32
        ss.append((a,b,c,d,e,f,g,h))
    return ss

# ================================================================
# Q1: Distribution of additive δa[r] mod 2^32
# Is it uniform or clustered?
# ================================================================

print("=" * 70)
print("Q1: Distribution of additive word differences")
print("=" * 70)

N = 5000
random.seed(42)

# Collect δa and δe (additive) at each round
da_vals = [[] for _ in range(65)]
de_vals = [[] for _ in range(65)]

for _ in range(N):
    M1 = [random.randint(0,MASK32) for _ in range(16)]
    M2 = list(M1); M2[0] ^= 1

    ss1 = sha_states(M1)
    ss2 = sha_states(M2)

    for r in range(65):
        da = (ss2[r][0] - ss1[r][0]) & MASK32  # additive diff of a
        de = (ss2[r][4] - ss1[r][4]) & MASK32  # additive diff of e
        da_vals[r].append(da)
        de_vals[r].append(de)

# Check: is δa[r] uniform mod 2^32?
# If uniform: δa mod m should be uniform mod m for any m.
# Test with m = 8, 256, 65536

print(f"\n  N={N}, δM = W[0] bit 0 flip")
print(f"\n  Uniformity test: chi² of δa[r] mod m")
print(f"\n  {'r':>3} | {'mod 2':>6} | {'mod 4':>6} | {'mod 8':>6} | {'mod 256':>7} | note")

for r in [1, 2, 3, 4, 5, 8, 16, 32, 48, 63, 64]:
    results = []
    for m in [2, 4, 8, 256]:
        counts = defaultdict(int)
        for v in da_vals[r]:
            counts[v % m] += 1
        expected = N / m
        chi2 = sum((counts[k] - expected)**2 / expected for k in range(m))
        # p-value approx: chi2 with m-1 DOF
        # Significant if chi2 > 2*(m-1) roughly
        pval_approx = chi2 / (m - 1) if m > 1 else 0
        results.append(pval_approx)

    note = ""
    if any(p > 3 for p in results): note = "★ NON-UNIFORM"

    print(f"  {r:>3} | {results[0]:>6.2f} | {results[1]:>6.2f} | {results[2]:>6.2f} | {results[3]:>7.2f} | {note}")

# ================================================================
# Q2: Additive recurrence — is δa[r+1] predicted by δa[r]?
# From SHA-256: a[r+1] = T1[r] + T2[r]
# δa[r+1] = δT1[r] + δT2[r] (mod 2^32)
# δT2[r] = δ(Σ₀(a[r]) + Maj(a[r], a[r-1], a[r-2]))
#
# If δa is small: δΣ₀(a[r]) ≈ Σ₀(δa[r]) (linearization)
# Does δa[r+1] ≈ f(δa[r]) for some f?
# ================================================================

print()
print("=" * 70)
print("Q2: Additive recurrence δa[r+1] = f(δa[r])?")
print("=" * 70)

# Correlation: δa[r+1] vs δa[r] mod 2^k
print(f"\n  Correlation of δa[r+1] with δa[r] (mod 2^k)")
print(f"\n  {'r':>3} | {'corr_raw':>8} | {'corr_mod256':>11} | {'corr_mod8':>9} | note")

def mod_corr(x_list, y_list, m):
    """Pearson correlation of x mod m vs y mod m."""
    xm = [x % m for x in x_list]
    ym = [y % m for y in y_list]
    n = len(xm)
    mx = sum(xm)/n; my = sum(ym)/n
    sx = math.sqrt(sum((a-mx)**2 for a in xm)/n)
    sy = math.sqrt(sum((a-my)**2 for a in ym)/n)
    if sx*sy == 0: return 0
    return sum((a-mx)*(b-my) for a,b in zip(xm,ym))/(n*sx*sy)

for r in [1, 2, 3, 4, 5, 8, 16, 32, 48, 63]:
    if r >= 64: continue
    c_raw = mod_corr(da_vals[r], da_vals[r+1], 2**32)
    c_256 = mod_corr(da_vals[r], da_vals[r+1], 256)
    c_8 = mod_corr(da_vals[r], da_vals[r+1], 8)

    note = ""
    if abs(c_raw) > 0.05 or abs(c_256) > 0.05: note = "★ CORRELATED"
    print(f"  {r:>3} | {c_raw:>+8.4f} | {c_256:>+11.4f} | {c_8:>+9.4f} | {note}")

# ================================================================
# Q3: Створочное число (a-e difference) — additive view
# створочное = a[r] - e[r] mod 2^32
# From methodology: a[r]-e[r] = T2[r-1] - a[r-4]
# This is an EXACT identity. How does it look additively?
# ================================================================

print()
print("=" * 70)
print("Q3: Створочне (a-e) additive structure")
print("=" * 70)

# For PAIRS: δ(a-e)[r] = δa[r] - δe[r] mod 2^32
# This should = δT2[r-1] - δa[r-4] (exact)
# Does δ(a-e) have structure that δa and δe separately don't?

print(f"\n  δ(a-e) = δa - δe (mod 2^32) at each round:")
print(f"\n  {'r':>3} | {'E[δ(a-e)]':>12} | {'std':>12} | {'mod8_chi2':>9} | note")

for r in [1, 2, 3, 4, 5, 8, 16, 32, 48, 64]:
    dae = [(da_vals[r][i] - de_vals[r][i]) & MASK32 for i in range(N)]

    avg = sum(dae) / N
    std = math.sqrt(sum((x - avg)**2 for x in dae) / N)

    # mod 8 chi2
    counts = defaultdict(int)
    for v in dae:
        counts[v % 8] += 1
    chi2 = sum((counts[k] - N/8)**2 / (N/8) for k in range(8))
    chi2_norm = chi2 / 7

    note = ""
    if chi2_norm > 3: note = "★ STRUCTURED mod 8"
    if std < 2**31 * 0.5: note += " LOW STD"

    print(f"  {r:>3} | {avg:>12.0f} | {std:>12.0f} | {chi2_norm:>9.2f} | {note}")


# ================================================================
# Q4: v₂ (2-adic valuation) of δa[r]
# v₂(x) = number of trailing zeros in binary representation
# If v₂(δa[r]) is structured → p-adic structure at WORD level
# ================================================================

print()
print("=" * 70)
print("Q4: 2-adic valuation v₂(δa[r])")
print("=" * 70)

def v2(x):
    """2-adic valuation: number of trailing zeros."""
    if x == 0: return 32  # convention
    v = 0
    while (x >> v) & 1 == 0:
        v += 1
    return v

print(f"\n  {'r':>3} | {'E[v₂(δa)]':>10} | {'P(v₂≥5)':>8} | {'P(v₂=0)':>8} | {'E[v₂] expected':>14} | note")

# For uniform random 32-bit: E[v₂] = Σ k/2^(k+1) = 1 - 2^{-32} ≈ 1.0
# P(v₂ = k) = 1/2^(k+1) for k=0..30, P(v₂ ≥ 31) = 1/2^31
# P(v₂ ≥ 5) = 1/2^5 = 0.03125
# P(v₂ = 0) = 1/2

expected_v2 = sum(k / 2**(k+1) for k in range(32))

for r in [1, 2, 3, 4, 5, 8, 16, 32, 48, 63, 64]:
    v2_vals = [v2(da_vals[r][i]) for i in range(N)]
    avg_v2 = sum(v2_vals) / N
    p_ge5 = sum(1 for v in v2_vals if v >= 5) / N
    p_eq0 = sum(1 for v in v2_vals if v == 0) / N

    note = ""
    if abs(avg_v2 - expected_v2) > 0.3: note = "★ ANOMALOUS v₂"
    if p_ge5 > 0.05: note += " HIGH v₂"
    if r <= 3 and avg_v2 > 2: note += " (early round, large δa)"

    print(f"  {r:>3} | {avg_v2:>10.3f} | {p_ge5:>8.4f} | {p_eq0:>8.4f} | {expected_v2:>14.3f} | {note}")


# ================================================================
# Q5: Additive relation between δa and δe across ROUNDS
# δe[r+1] = δd[r] + δT1[r] = δa[r-3] + δT1[r]
# Does knowing δa[r] help predict δe[r+4]?
# (This is the створочне relation in DIFFERENCE form)
# ================================================================

print()
print("=" * 70)
print("Q5: Cross-chain additive relation δa[r] → δe[r+4]")
print("=" * 70)

print(f"\n  Correlation: δa[r] vs δe[r+4] (additive, mod 2^32)")
print(f"  From SHA-256 structure: δe[r+4] = δa[r] + δT1[r+3]")
print(f"  If δT1 independent of δa → no correlation")
print(f"\n  {'r':>3} | {'corr(δa[r], δe[r+4])':>21} | note")

for r in range(1, 56):
    c = mod_corr(da_vals[r], de_vals[r+4], 2**32)
    if r <= 5 or abs(c) > 0.03 or r % 16 == 0:
        note = "★" if abs(c) > 0.05 else ""
        print(f"  {r:>3} | {c:>+21.4f} | {note}")


# ================================================================
# Q6: KEY TEST — is δa[r] mod 2^k predictable for small k?
# If δa[r] mod 8 is structured → carry structure at word level
# ================================================================

print()
print("=" * 70)
print("Q6: δa[r] mod 2^k — per-round distribution")
print("=" * 70)

# Focus on mod 8 (3 low bits of additive difference)
print(f"\n  δa[r] mod 8 distribution (should be uniform = 12.5% each):")
print(f"  N={N}")

for r in [1, 2, 3, 5, 8, 16, 32, 64]:
    counts = [0]*8
    for v in da_vals[r]:
        counts[v % 8] += 1
    probs = [c/N for c in counts]
    max_dev = max(abs(p - 0.125) for p in probs)

    bar = ""
    for i in range(8):
        bar += f" {i}:{probs[i]:.3f}"
    note = ""
    if max_dev > 0.03: note = " ★ BIASED"
    if r <= 2 and da_vals[r][0] < 256: note += " (small δa)"

    print(f"  r={r:>2}: {bar} {note}")


# ================================================================
# SYNTHESIS
# ================================================================

print()
print("=" * 70)
print("SYNTHESIS: Word-level differential structure")
print("=" * 70)

print("""
  Looking at WORDS as NUMBERS (mod 2^32), not as bit vectors.

  Results above show whether:
  - δa[r] is uniformly distributed mod small numbers
  - consecutive δa values are correlated arithmetically
  - створочне (a-e) has word-level structure
  - 2-adic valuation reveals p-adic patterns
  - cross-chain δa→δe relation is exploitable

  Key: if ANY mod-k test shows structure at round > 8,
  that's information surviving the thermostat —
  invisible to bit-level analysis but visible at word level.
""")
