#!/usr/bin/env python3
"""
Pipeline approach: Wang pairs as a STREAM.

Key question that hasn't been tested:
Are δH[i] correlated for WANG PAIRS specifically?
Dir3 tested random pairs → independent.
But Wang pairs have structure (δe=0, δa~16).
Does this structure create δH correlation?

If yes → cascaded filtering is cheaper than independent birthday.
"""

import random
from collections import defaultdict
import math

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

def sha256(M):
    W=list(M)+[0]*(64-len(M))
    for i in range(16,64):W[i]=(ssig1(W[i-2])+W[i-7]+ssig0(W[i-15])+W[i-16])&MASK32
    a,b,c,d,e,f,g,h=IV
    for r in range(64):
        T1=(h+sig1(e)+ch(e,f,g)+K[r]+W[r])&MASK32;T2=(sig0(a)+maj(a,b,c))&MASK32
        h,g,f=g,f,e;e=(d+T1)&MASK32;d,c,b=c,b,a;a=(T1+T2)&MASK32
    return tuple((s+v)&MASK32 for s,v in zip([a,b,c,d,e,f,g,h],IV))

def wang_pair(seed_val):
    """Generate one Wang pair. Returns (H1, H2, δH)."""
    random.seed(seed_val)
    W1 = [random.randint(0,MASK32) for _ in range(16)]
    DW = [0]*16; DW[0] = 1

    def state_at(M, R):
        W=list(M)+[0]*(64-len(M))
        for i in range(16,64):W[i]=(ssig1(W[i-2])+W[i-7]+ssig0(W[i-15])+W[i-16])&MASK32
        a,b,c,d,e,f,g,h=IV
        for r in range(R):
            T1=(h+sig1(e)+ch(e,f,g)+K[r]+W[r])&MASK32;T2=(sig0(a)+maj(a,b,c))&MASK32
            h,g,f=g,f,e;e=(d+T1)&MASK32;d,c,b=c,b,a;a=(T1+T2)&MASK32
        return (a,b,c,d,e,f,g,h)

    for step in range(15):
        wi = step+1
        W2 = [(W1[i]+DW[i])&MASK32 for i in range(16)]
        s1 = state_at(W1, step+2)
        s2 = state_at(W2, step+2)
        de = (s2[4]-s1[4])&MASK32
        DW[wi] = (-de)&MASK32

    W2 = [(W1[i]+DW[i])&MASK32 for i in range(16)]
    H1 = sha256(W1)
    H2 = sha256(W2)
    dH = tuple(H1[i] ^ H2[i] for i in range(8))
    return H1, H2, dH


# ================================================================
# TEST 1: δH word correlation for WANG pairs
# ================================================================

print("=" * 70)
print("TEST 1: δH[i] correlation for Wang pairs")
print("=" * 70)

N = 5000
dH_words = [[] for _ in range(8)]

for seed in range(N):
    _, _, dH = wang_pair(seed * 137 + 1)
    for i in range(8):
        dH_words[i].append(dH[i])

# HW correlation between words
def pearson(x, y):
    n=len(x);mx=sum(x)/n;my=sum(y)/n
    sx=math.sqrt(sum((xi-mx)**2 for xi in x)/n)
    sy=math.sqrt(sum((yi-my)**2 for yi in y)/n)
    if sx*sy==0:return 0
    return sum((xi-mx)*(yi-my) for xi,yi in zip(x,y))/(n*sx*sy)

hw_words = [[hw(v) for v in dH_words[i]] for i in range(8)]

print(f"\n  N={N} Wang pairs")
print(f"\n  Pearson corr of HW(δH[i]) vs HW(δH[j]):")
print(f"        ", end="")
for j in range(8): print(f" H[{j}] ", end="")
print()

for i in range(8):
    print(f"  H[{i}]", end="")
    for j in range(8):
        if i==j:
            print(f"    -- ", end="")
        else:
            c = pearson(hw_words[i], hw_words[j])
            flag = "★" if abs(c) > 0.05 else " "
            print(f" {c:+.3f}{flag}", end="")
    print()

# ================================================================
# TEST 2: Conditional probability — does δH[7]=0 help δH[6]?
# Among Wang pairs with δH[7]=0: is δH[6] biased?
# ================================================================

print()
print("=" * 70)
print("TEST 2: Conditional filtering — δH[7]=0 → bias in δH[6]?")
print("=" * 70)

# Birthday among Wang pairs on H[7]
N2 = 300000
ht7 = {}  # δH[7] → list of (seed, δH)
collisions_h7 = []

for seed in range(N2):
    H1, H2, dH = wang_pair(seed * 73 + 5)
    key = dH[7]  # δH[7]
    if key in ht7:
        prev_seed, prev_dH = ht7[key]
        collisions_h7.append((prev_dH, dH))
    else:
        ht7[key] = (seed, dH)

print(f"\n  N={N2} Wang pairs, birthday on δH[7]")
print(f"  δH[7]-collisions found: {len(collisions_h7)}")

expected = N2*(N2-1)/(2*2**32)
print(f"  Expected for 32-bit birthday: {expected:.1f}")

if collisions_h7:
    # Among δH[7]-matching pairs: check δH[6]
    h6_match = 0
    h0_match = 0
    hw_h6 = []

    for dH_a, dH_b in collisions_h7:
        # These pairs have same δH[7]. Do they also share δH[6]?
        if dH_a[6] == dH_b[6]:
            h6_match += 1
        if dH_a[0] == dH_b[0]:
            h0_match += 1
        # HW of difference in δH[6]
        hw_h6.append(hw(dH_a[6] ^ dH_b[6]))

    p_h6 = h6_match / len(collisions_h7) if collisions_h7 else 0
    p_h0 = h0_match / len(collisions_h7) if collisions_h7 else 0
    expected_match = 1 / 2**32

    print(f"\n  Among δH[7]-collision pairs:")
    print(f"    P(δH[6] also matches) = {h6_match}/{len(collisions_h7)} = {p_h6:.6f}")
    print(f"    P(δH[0] also matches) = {h0_match}/{len(collisions_h7)} = {p_h0:.6f}")
    print(f"    Expected if independent: {expected_match:.2e}")

    if h6_match > 0:
        ratio = p_h6 / expected_match
        print(f"    Ratio H[6]: {ratio:.1f}× {'★ CORRELATED!' if ratio > 10 else ''}")
    if h0_match > 0:
        ratio0 = p_h0 / expected_match
        print(f"    Ratio H[0]: {ratio0:.1f}×")

    if hw_h6:
        print(f"\n    E[HW(δH[6]_a XOR δH[6]_b)] = {sum(hw_h6)/len(hw_h6):.1f} (expect 16 if independent)")
else:
    print(f"  No δH[7]-collisions found. Need more pairs.")

# ================================================================
# TEST 3: Pipeline simulation
# Wang stream → filter by δH[7] low bits → filter by δH[6]
# ================================================================

print()
print("=" * 70)
print("TEST 3: Pipeline — filter Wang stream by partial δH match")
print("=" * 70)

# Stage 1: generate Wang pairs
# Stage 2: birthday on δH[7] (first 16 bits only → cheaper)
# Stage 3: from δH[7] matches, check δH[6]

N3 = 500000

# Birthday on lower 16 bits of δH[7]
ht_lo = defaultdict(list)  # δH[7] lower 16 → list of (seed, full_dH)

for seed in range(N3):
    H1, H2, dH = wang_pair(seed * 41 + 3)
    key = dH[7] & 0xFFFF  # lower 16 bits
    ht_lo[key].append((seed, dH))

# Count pairs per bucket
n_pairs = 0
h7_full_match = 0
h7h6_partial_match = 0

for key, entries in ht_lo.items():
    if len(entries) >= 2:
        for i in range(len(entries)):
            for j in range(i+1, min(len(entries), i+5)):  # limit per bucket
                n_pairs += 1
                dH_a = entries[i][1]
                dH_b = entries[j][1]

                # Do full δH[7] match?
                if dH_a[7] == dH_b[7]:
                    h7_full_match += 1

                    # δH[6] lower 16 also match?
                    if (dH_a[6] & 0xFFFF) == (dH_b[6] & 0xFFFF):
                        h7h6_partial_match += 1

print(f"  N={N3} Wang pairs")
print(f"  δH[7] lower-16 birthday pairs: {n_pairs}")
print(f"  Full δH[7] matches: {h7_full_match}")
print(f"  Full δH[7] + partial δH[6] (lower 16): {h7h6_partial_match}")

if n_pairs > 0:
    p_full7 = h7_full_match / n_pairs
    print(f"\n  P(full δH[7] | lower-16 match) = {p_full7:.4f}")
    print(f"  Expected: 1/2^16 = {1/2**16:.6f}")

    if h7_full_match > 0:
        p_h6_given_h7 = h7h6_partial_match / h7_full_match
        print(f"  P(δH[6] lower-16 match | full δH[7] match) = {p_h6_given_h7:.4f}")
        print(f"  Expected if independent: 1/2^16 = {1/2**16:.6f}")

        if h7h6_partial_match > 0:
            ratio = p_h6_given_h7 / (1/2**16)
            print(f"  Ratio: {ratio:.2f}× {'★★ SIGNAL!' if ratio > 2 else '≈ independent'}")


# ================================================================
# TEST 4: Bit-level correlation in δH for Wang pairs
# Does bit b of δH[7] predict bit b of δH[6]?
# ================================================================

print()
print("=" * 70)
print("TEST 4: Bit-level δH correlation within Wang pairs")
print("=" * 70)

N4 = 10000
bit_agree = [[0]*32 for _ in range(8)]  # bit_agree[word][bit]

for seed in range(N4):
    _, _, dH = wang_pair(seed * 97 + 11)
    for i in range(8):
        for b in range(32):
            # Does bit b of δH[i] = bit b of δH[(i+1)%8]?
            bi = (dH[i] >> b) & 1
            bj = (dH[(i+1)%8] >> b) & 1
            if bi == bj:
                bit_agree[i][b] += 1

print(f"  N={N4} Wang pairs")
print(f"  P(δH[i][b] = δH[i+1][b]) for adjacent words:")
print(f"\n  {'word_pair':>10} | {'max_P':>5} {'min_P':>5} {'avg_P':>5} | {'max_bit':>7} | note")

for i in range(8):
    j = (i+1) % 8
    probs = [bit_agree[i][b] / N4 for b in range(32)]
    max_p = max(probs)
    min_p = min(probs)
    avg_p = sum(probs) / 32
    max_bit = probs.index(max_p)

    note = ""
    if max_p > 0.52: note = "★ BIAS"
    if max_p > 0.55: note = "★★ STRONG"

    print(f"  H[{i}]↔H[{j}] | {max_p:.3f} {min_p:.3f} {avg_p:.3f} | bit {max_bit:>2}    | {note}")


# ================================================================
# SYNTHESIS
# ================================================================

print()
print("=" * 70)
print("SYNTHESIS: Pipeline approach")
print("=" * 70)

print("""
  Pipeline idea: generate Wang pairs as stream, filter cascadingly.
  Requires: δH word correlation for WANG pairs (not random pairs).

  Results:
  - HW correlation between δH words: checked above
  - Conditional P(δH[6] | δH[7] match): checked above
  - Bit-level δH[i]↔δH[i+1]: checked above

  If all ≈ independent → pipeline = standard birthday = 2^128.
  If correlated → pipeline potentially cheaper.
""")
