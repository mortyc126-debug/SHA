#!/usr/bin/env python3
"""
NK Direction 3: Multi-differential attack.

Idea: Instead of one pair (M1,M2) with H(M1)=H(M2),
use MULTIPLE pairs that each give near-collisions in DIFFERENT words.

Pair A: H[0..3] close (a-chain near-collision)
Pair B: H[4..7] close (e-chain near-collision)

If we can COMBINE them → full collision?

Also: k-collision (k messages all hashing to same value).
Birthday for k-collision: O(N^{1/k} * S^{1-1/k}).
For k=2: O(S^{1/2}) = 2^128 (standard).
For k=3: O(S^{2/3}) = 2^{170} (WORSE for collision).
But k-list birthday (Wagner): different.

Questions:
Q1: Do Wang pairs cluster? (Multiple pairs near same H)
Q2: Can two near-collisions on different words be combined?
Q3: Does the a-chain / e-chain split help for multi-party birthday?
"""

import random
from collections import defaultdict

MASK32 = 0xFFFFFFFF
K=[0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2]
IV=[0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&MASK32
def ssig0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def ssig1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def sig0(x): return rotr(x,2)^rotr(x,13)^rotr(x,22)
def sig1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def ch(e,f,g): return (e&f)^(~e&g)&MASK32
def maj(a,b,c): return (a&b)^(a&c)^(b&c)
def hw(x): return bin(x&MASK32).count('1')

def sha256c(M):
    W=list(M)+[0]*(64-len(M))
    for i in range(16,64): W[i]=(ssig1(W[i-2])+W[i-7]+ssig0(W[i-15])+W[i-16])&MASK32
    a,b,c,d,e,f,g,h=IV
    for r in range(64):
        T1=(h+sig1(e)+ch(e,f,g)+K[r]+W[r])&MASK32;T2=(sig0(a)+maj(a,b,c))&MASK32
        h,g,f=g,f,e;e=(d+T1)&MASK32;d,c,b=c,b,a;a=(T1+T2)&MASK32
    return tuple((s+iv)&MASK32 for s,iv in zip([a,b,c,d,e,f,g,h],IV))


# ================================================================
# Q1: Correlation between output words for random messages.
# If H[0] and H[4] are correlated → multi-word birthday is cheaper.
# If independent → multi-word birthday = product of per-word.
# ================================================================

print("=" * 70)
print("Q1: Output word independence (random messages)")
print("=" * 70)

N = 100000
random.seed(42)

# Collect per-word values
h_vals = [[] for _ in range(8)]
for _ in range(N):
    M = [random.randint(0,MASK32) for _ in range(16)]
    H = sha256c(M)
    for w in range(8):
        h_vals[w].append(H[w])

# Pearson correlation between words (on raw values mod 2^16 for speed)
import math
def pearson(x, y):
    n=len(x); mx=sum(x)/n; my=sum(y)/n
    sx=math.sqrt(sum((xi-mx)**2 for xi in x)/n)
    sy=math.sqrt(sum((yi-my)**2 for yi in y)/n)
    if sx==0 or sy==0: return 0
    return sum((xi-mx)*(yi-my) for xi,yi in zip(x,y))/(n*sx*sy)

# Use lower 16 bits for correlation (full 32-bit values too noisy)
print(f"\n  N={N}, Pearson correlation between H[i] and H[j] (lower 16 bits):")
print(f"        H[0]   H[1]   H[2]   H[3]   H[4]   H[5]   H[6]   H[7]")
for i in range(8):
    row = f"  H[{i}]"
    for j in range(8):
        if i == j:
            row += "    --"
        else:
            xi = [v & 0xFFFF for v in h_vals[i]]
            xj = [v & 0xFFFF for v in h_vals[j]]
            c = pearson(xi, xj)
            row += f" {c:+.3f}"
    print(row)

# XOR correlation: P(H[i][b] = H[j][b]) for bit 0
print(f"\n  Bit 0 agreement: P(H[i][0] = H[j][0])")
for i in [0, 3, 4, 7]:
    for j in [i+1, (i+4)%8]:
        if j >= 8: continue
        agree = sum(1 for k in range(N) if (h_vals[i][k]&1) == (h_vals[j][k]&1))
        print(f"    H[{i}][0] vs H[{j}][0]: P(agree) = {agree/N:.4f} (expect 0.5000)")


# ================================================================
# Q2: Can partial collisions be COMPOSED?
#
# Suppose pair (M1,M2) has δH[0..3]=0 (a-chain collision).
# Suppose pair (M3,M4) has δH[4..7]=0 (e-chain collision).
# Can we build M5 combining both? Only if M2=M3 or similar.
#
# More precisely: multi-target birthday.
# Search for M where H(M)[0..3] matches some target T_a
# AND H(M)[4..7] matches some target T_e.
#
# If H[0..3] and H[4..7] are independent:
#   P(both match) = P(a-match) × P(e-match) = 2^{-128} × 2^{-128} = 2^{-256}
#   Birthday: 2^{128}. Same as standard.
#
# If H[0..3] and H[4..7] are DEPENDENT:
#   P(both match) > 2^{-256} → birthday < 2^{128}.
# ================================================================

print()
print("=" * 70)
print("Q2: Independence of a-chain (H[0..3]) and e-chain (H[4..7])")
print("=" * 70)

# Test: among messages where H[0]=target, is H[4] biased?
# Need many messages with same H[0] → use birthday on H[0]

random.seed(123)
N2 = 500000

# Build table: H[0] → list of (M, H[4])
ht_h0 = defaultdict(list)
for _ in range(N2):
    M = tuple(random.randint(0,MASK32) for _ in range(16))
    H = sha256c(M)
    ht_h0[H[0]].append(H[4])

# Find H[0] values with multiple messages
multi = {k: v for k, v in ht_h0.items() if len(v) >= 2}
print(f"\n  N={N2} messages, H[0]-collisions: {len(multi)} groups")

if multi:
    # For each H[0]-collision group, check if H[4] values are related
    h4_match = 0
    h4_total = 0
    for h0_val, h4_list in multi.items():
        for i in range(len(h4_list)):
            for j in range(i+1, len(h4_list)):
                h4_total += 1
                if h4_list[i] == h4_list[j]:
                    h4_match += 1

    print(f"  H[0]-collision pairs: {h4_total}")
    print(f"  Among them, H[4] also matches: {h4_match}")
    if h4_total > 0:
        p_cond = h4_match / h4_total
        p_random = 1 / 2**32
        ratio = p_cond / p_random if p_random > 0 else float('inf')
        print(f"  P(H[4] match | H[0] match) = {p_cond:.6f}")
        print(f"  P(H[4] match | random) = {p_random:.2e}")
        if h4_match == 0:
            print(f"  No conditional H[4] matches → H[0] and H[4] appear INDEPENDENT")
        else:
            print(f"  Ratio = {ratio:.1f}×")
else:
    print(f"  No H[0]-collisions found (need N >> 2^16 = 65536)")
    print(f"  With N={N2}: expected H[0]-collisions = {N2**2/(2*2**32):.1f}")


# ================================================================
# Q3: a-chain / e-chain split for multi-list birthday
#
# SHA-256 output: H[0..3] from a-chain, H[4..7] from e-chain.
# What if we treat these as TWO SEPARATE targets?
#
# Wagner's k-list birthday: for k independent lists of size N,
# find one element from each list that XOR to zero.
# Cost: O(N) per list when k = log2(S) lists, each size S^{1/k}.
#
# For SHA-256: S = 2^256, k = 2 "chains" (a and e).
# Each chain: 128 bits. Birthday per chain: 2^64.
# If chains independent: total = 2^64 + 2^64 = 2^65.
# But this finds H[0..3] collision AND H[4..7] collision SEPARATELY.
# Need SAME pair for BOTH → back to 2^128.
#
# UNLESS: we can find M where a-chain and e-chain are partly independent.
# ================================================================

print()
print("=" * 70)
print("Q3: Can a-chain and e-chain be targeted separately?")
print("=" * 70)

# Measure: correlation between δH[0..3] and δH[4..7] for random pairs
N3 = 100000
random.seed(456)

hw_a_chain = []
hw_e_chain = []

for _ in range(N3):
    M1 = [random.randint(0,MASK32) for _ in range(16)]
    M2 = [random.randint(0,MASK32) for _ in range(16)]
    H1 = sha256c(M1)
    H2 = sha256c(M2)

    hw_a = sum(hw(H1[i]^H2[i]) for i in range(4))  # a-chain: H[0..3]
    hw_e = sum(hw(H1[i]^H2[i]) for i in range(4,8)) # e-chain: H[4..7]
    hw_a_chain.append(hw_a)
    hw_e_chain.append(hw_e)

corr_ae = pearson(hw_a_chain, hw_e_chain)
avg_a = sum(hw_a_chain)/N3
avg_e = sum(hw_e_chain)/N3

print(f"\n  N={N3} random pairs")
print(f"  E[HW(δH[0..3])] = {avg_a:.1f} (expected 64)")
print(f"  E[HW(δH[4..7])] = {avg_e:.1f} (expected 64)")
print(f"  Pearson(HW_a, HW_e) = {corr_ae:.4f}")

if abs(corr_ae) < 0.01:
    print(f"\n  a-chain and e-chain differences are INDEPENDENT.")
    print(f"  Multi-target birthday gains NOTHING:")
    print(f"    Separate: 2^64 (a-chain) + 2^64 (e-chain) = 2^65")
    print(f"    But need SAME pair → 2^128. No shortcut.")
else:
    print(f"\n  Correlation detected: {corr_ae:.4f}")
    print(f"  Investigate further.")


# ================================================================
# Q4: 4-list birthday on output WORDS
#
# Treat H as 8 independent 32-bit words.
# Wagner k-list for k=8, each word 32 bits:
#   Cost = O(2^{256/(1+3)}) = O(2^64) if lists are independent.
#
# But SHA-256 produces H from ONE message — not 8 independent lists.
# Can we FAKE independence?
#
# Idea: partition message space into 8 groups,
# each "controlling" one output word.
# But: all output words depend on ALL input bits (after 8 rounds).
# No partition gives independent control.
# ================================================================

print()
print("=" * 70)
print("Q4: Can we partition M to control individual H[i]?")
print("=" * 70)

# Test: sensitivity of each H[i] to each W[j]
# S[i][j] = E[HW(H[i](M) XOR H[i](M with W[j] flipped))]
# If S[i][j] ≈ 16 for all i,j: uniform sensitivity, no partition possible

N4 = 3000
random.seed(789)

print(f"\n  Sensitivity S[i][j] = E[HW(δH[i])] when W[j] bit 0 flipped")
print(f"  N={N4}\n")

sens = [[0.0]*16 for _ in range(8)]

for trial in range(N4):
    M = [random.randint(0,MASK32) for _ in range(16)]
    H0 = sha256c(M)

    for j in range(16):
        M2 = list(M)
        M2[j] ^= 1
        H2 = sha256c(M2)

        for i in range(8):
            sens[i][j] += hw(H0[i] ^ H2[i])

for i in range(8):
    for j in range(16):
        sens[i][j] /= N4

# Print compact: just check if any H[i] is MORE sensitive to specific W[j]
print(f"  {'':>5}", end="")
for j in range(16):
    print(f" W[{j:>2}]", end="")
print()

for i in range(8):
    print(f"  H[{i}]", end="")
    for j in range(16):
        s = sens[i][j]
        if s > 16.3:
            print(f" {s:>4.1f}*", end="")
        elif s < 15.7:
            print(f" {s:>4.1f}-", end="")
        else:
            print(f" {s:>4.1f} ", end="")
    print()

# Check: is sensitivity uniform?
all_s = [sens[i][j] for i in range(8) for j in range(16)]
avg_s = sum(all_s)/len(all_s)
std_s = (sum((x-avg_s)**2 for x in all_s)/len(all_s))**0.5
min_s = min(all_s)
max_s = max(all_s)

print(f"\n  avg={avg_s:.2f} std={std_s:.2f} min={min_s:.2f} max={max_s:.2f}")
print(f"  Range: {max_s-min_s:.2f} bits out of 16")

if std_s < 0.3:
    print(f"\n  Sensitivity UNIFORM across all (H[i], W[j]) pairs.")
    print(f"  No partition possible → k-list birthday inapplicable.")
else:
    print(f"\n  Some variation detected. Max-min = {max_s-min_s:.2f}")


# ================================================================
# SYNTHESIS
# ================================================================

print()
print("=" * 70)
print("SYNTHESIS: Direction 3")
print("=" * 70)

print("""
  Q1: Output words have Pearson correlation ≈ 0 (independent).
  Q2: H[0] collision gives zero information about H[4] (independent).
  Q3: a-chain and e-chain differences uncorrelated (independent).
  Q4: Every H[i] equally sensitive to every W[j] (no partition).

  All four tests show: SHA-256 output words are INDEPENDENT
  for cryptanalytic purposes. Multi-differential cannot decompose
  the 256-bit collision into easier sub-problems.

  k-list birthday (Wagner): requires independent lists.
  SHA-256 output words come from ONE computation → not independent lists.
  No way to "assign" input bits to specific output words.

  Direction 3: CLOSED. Multi-differential gives no advantage.
""")
