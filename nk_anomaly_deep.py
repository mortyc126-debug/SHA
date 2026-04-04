#!/usr/bin/env python3
"""
NK: Deep verification of W[1..15]=0 birthday anomaly.

Questions:
1. Is it only H[7] or ALL output words?
2. How much total image compression?
3. Other structured subsets — do they also compress?
4. Connection to GBP: can multi-word collisions be found faster?
"""

import random
from collections import defaultdict
import time, math

MASK32 = 0xFFFFFFFF
K = [0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2]
IV = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&MASK32
def sig0(x): return rotr(x,2)^rotr(x,13)^rotr(x,22)
def sig1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def ssig0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def ssig1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def ch(e,f,g): return (e&f)^(~e&g)&MASK32
def maj(a,b,c): return (a&b)^(a&c)^(b&c)

def sha256c(M):
    W=list(M)+[0]*(64-len(M))
    for i in range(16,64): W[i]=(ssig1(W[i-2])+W[i-7]+ssig0(W[i-15])+W[i-16])&MASK32
    a,b,c,d,e,f,g,h=IV
    for r in range(64):
        T1=(h+sig1(e)+ch(e,f,g)+K[r]+W[r])&MASK32
        T2=(sig0(a)+maj(a,b,c))&MASK32
        h,g,f=g,f,e;e=(d+T1)&MASK32;d,c,b=c,b,a;a=(T1+T2)&MASK32
    return tuple((s+iv)&MASK32 for s,iv in zip([a,b,c,d,e,f,g,h],IV))

# ================================================================
# TEST 1: Per-word collision rate for W[1..15]=0 vs random
# Check ALL 8 output words independently
# ================================================================

print("=" * 70)
print("TEST 1: Per-word H[i] collision rate (5 seeds × 300K)")
print("=" * 70)

N = 300000
SEEDS = 5

for label, gen_fn in [
    ("random",     lambda: [random.randint(0,MASK32) for _ in range(16)]),
    ("W[1..15]=0", lambda: [random.randint(0,MASK32)]+[0]*15),
]:
    word_colls = [0]*8
    for seed in range(SEEDS):
        random.seed(seed*777+1)
        hts = [dict() for _ in range(8)]
        for _ in range(N):
            H = sha256c(gen_fn())
            for w in range(8):
                if H[w] in hts[w]:
                    word_colls[w] += 1
                else:
                    hts[w][H[w]] = 1

    expected = N*(N-1)/(2*2**32) * SEEDS
    print(f"\n  {label}:  (expected per word = {expected:.1f})")
    for w in range(8):
        ratio = word_colls[w] / max(expected, 0.1)
        bar = "█" * int(ratio * 20)
        flag = " ★" if ratio > 1.5 else ""
        print(f"    H[{w}]: {word_colls[w]:>4} coll, ratio={ratio:.2f}× {bar}{flag}")


# ================================================================
# TEST 2: Pair-word collisions (H[i]=H[j] for i≠j simultaneously)
# Do words collide TOGETHER more often?
# ================================================================

print()
print("=" * 70)
print("TEST 2: Joint 2-word collision rate")
print("=" * 70)

N2 = 300000

for label, gen_fn in [
    ("random",     lambda: [random.randint(0,MASK32) for _ in range(16)]),
    ("W[1..15]=0", lambda: [random.randint(0,MASK32)]+[0]*15),
]:
    random.seed(42)
    # Track pairs of words
    pair_ht = {}  # (H[i], H[j]) -> count for selected pairs
    pair_colls = defaultdict(int)

    ht_pairs = {}  # (H[0], H[7]) as key for joint collision test
    joint_coll = 0

    ht_h07 = {}
    for trial in range(N2):
        H = sha256c(gen_fn())
        key = (H[0], H[7])
        if key in ht_h07:
            joint_coll += 1
        else:
            ht_h07[key] = 1

    expected_joint = N2*(N2-1)/(2*2**64)  # 64-bit key
    ratio = joint_coll / max(expected_joint, 1e-10)
    print(f"  {label}: H[0]∧H[7] joint: {joint_coll} coll, expect={expected_joint:.4f}, ratio={ratio:.1f}×")


# ================================================================
# TEST 3: Other structured message subsets
# Which patterns give the most compression?
# ================================================================

print()
print("=" * 70)
print("TEST 3: Collision rate for different message structures")
print("=" * 70)

N3 = 300000

structures = [
    ("random",            lambda: [random.randint(0,MASK32) for _ in range(16)]),
    ("W[1..15]=0",        lambda: [random.randint(0,MASK32)]+[0]*15),
    ("W[0,2..15]=0",      lambda: [0,random.randint(0,MASK32)]+[0]*14),
    ("W[0..14]=0,W[15]",  lambda: [0]*15+[random.randint(0,MASK32)]),
    ("W[0]=W[1]=...=W[15]", lambda: (lambda w:[w]*16)(random.randint(0,MASK32))),
    ("W[1..15]=1",        lambda: [random.randint(0,MASK32)]+[1]*15),
    ("W[1..15]=0xFFFFFFFF", lambda: [random.randint(0,MASK32)]+[MASK32]*15),
    ("W[8..15]=0",        lambda: [random.randint(0,MASK32) for _ in range(8)]+[0]*8),
    ("only W[0],W[1]",    lambda: [random.randint(0,MASK32),random.randint(0,MASK32)]+[0]*14),
]

print(f"  N={N3}, 3 seeds per structure, H[7] collisions\n")
print(f"  {'structure':>25} | {'avg_coll':>8} | {'ratio':>6} | note")

results = []

for label, gen_fn in structures:
    total_coll = 0
    for seed in range(3):
        random.seed(seed*333+5)
        ht = {}
        coll = 0
        for _ in range(N3):
            h7 = sha256c(gen_fn())[7]
            if h7 in ht: coll += 1
            else: ht[h7] = 1
        total_coll += coll

    avg = total_coll / 3
    expected = N3*(N3-1)/(2*2**32)
    ratio = avg / max(expected, 0.1)
    note = ""
    if ratio > 1.8: note = "★★ STRONG"
    elif ratio > 1.3: note = "★ elevated"
    results.append((ratio, label, avg, note))

for ratio, label, avg, note in sorted(results, reverse=True):
    print(f"  {label:>25} | {avg:>8.1f} | {ratio:>5.1f}× | {note}")


# ================================================================
# TEST 4: WHY does W[1..15]=0 compress?
# Hypothesis: degenerate schedule reduces effective dimension.
#
# W[1..15]=0 means:
#   W[16] = σ₁(W[14]) + W[9] + σ₀(W[1]) + W[0]
#         = σ₁(0) + 0 + σ₀(0) + W[0] = W[0]
#   W[17] = σ₁(W[15]) + W[10] + σ₀(W[2]) + W[1] = 0
#   W[18] = σ₁(W[16]) + W[11] + σ₀(W[3]) + W[2] = σ₁(W[0])
#   W[19] = 0, W[20] = σ₁(σ₁(W[0])), W[21] = 0, ...
#
# Only EVEN-indexed W[16+2k] are nonzero, and all are functions of W[0]!
# This means: 32-bit input (W[0]) → 256-bit output (H).
# The function is 32→256, so H[7] has at most 2^32 distinct values.
# Birthday on 2^32 values: N*2/2^32 expected collisions.
# At N=300K: 300000^2/2^33 ≈ 10.5. This matches baseline!
#
# BUT: the observed ratio is 2× → only 2^31 distinct H[7] values.
# This means: the mapping W[0] → H[7] is NOT injective.
# On average, 2 different W[0] values map to the same H[7].
# ================================================================

print()
print("=" * 70)
print("TEST 4: Image size estimation for W[1..15]=0")
print("=" * 70)

# Direct measurement: how many unique H[7] values for N distinct W[0]?
random.seed(42)
N4 = 500000
seen_h7 = set()
seen_w0 = set()
dupes_h7 = 0

for _ in range(N4):
    w0 = random.randint(0, MASK32)
    if w0 in seen_w0:
        continue
    seen_w0.add(w0)

    M = [w0] + [0]*15
    h7 = sha256c(M)[7]

    if h7 in seen_h7:
        dupes_h7 += 1
    seen_h7.add(h7)

unique_w0 = len(seen_w0)
unique_h7 = len(seen_h7)
collision_rate = dupes_h7 / unique_w0

# Birthday estimate of image size:
# If |Im| = S, then expected duplicates ≈ N²/(2S)
# dupes = N²/(2S) → S = N²/(2*dupes)
if dupes_h7 > 0:
    estimated_S = unique_w0**2 / (2 * dupes_h7)
else:
    estimated_S = float('inf')

print(f"  Unique W[0] tested: {unique_w0}")
print(f"  Unique H[7] seen:  {unique_h7}")
print(f"  H[7] duplicates:   {dupes_h7}")
print(f"  Collision rate:    {collision_rate:.6f}")
print(f"  Estimated |Im(H[7])| ≈ {estimated_S:.0f} ≈ 2^{math.log2(estimated_S):.1f}")
print(f"  Full space: 2^32 = {2**32}")
print(f"  Compression ratio: {2**32 / estimated_S:.2f}×")

# Same for H[0]
seen_h0 = set()
dupes_h0 = 0
random.seed(42)
seen_w0_2 = set()
for _ in range(N4):
    w0 = random.randint(0, MASK32)
    if w0 in seen_w0_2: continue
    seen_w0_2.add(w0)
    h0 = sha256c([w0]+[0]*15)[0]
    if h0 in seen_h0: dupes_h0 += 1
    seen_h0.add(h0)

if dupes_h0 > 0:
    est_S0 = len(seen_w0_2)**2 / (2*dupes_h0)
    print(f"\n  H[0]: dupes={dupes_h0}, estimated |Im| ≈ 2^{math.log2(est_S0):.1f}, compression={2**32/est_S0:.2f}×")
else:
    print(f"\n  H[0]: dupes=0 (image ≈ full)")


# ================================================================
# TEST 5: Is this just birthday statistics of a 32→32 function?
#
# For ANY function f: {0,1}^32 → {0,1}^32, birthday expects
# N²/(2·2^32) collisions for N distinct inputs.
# A RANDOM function has |Im| ≈ (1-1/e)·2^32 ≈ 0.632·2^32.
# So birthday collision rate for random function:
# N²/(2·0.632·2^32) = 1.58 × N²/(2·2^32).
#
# Factor 1.58 comes from image compression of random function!
# Is our 2.09× just this 1.58× plus noise?
# ================================================================

print()
print("=" * 70)
print("TEST 5: Is 2× = random function image compression?")
print("=" * 70)

# For a random function f: {0,1}^32 → {0,1}^32:
# Expected unique outputs for N inputs: S = 2^32 × (1 - (1-1/2^32)^N)
# For N = 300000: S ≈ 300000 (N << 2^32, so almost all unique)
# Birthday collisions ≈ N²/(2·2^32) regardless of image size.
#
# But WAIT: our hash table only stores FIRST occurrence.
# So collision = second W[0] with same H[7].
# This is exactly birthday for a 32→32 function.
#
# For N=300K << 2^32: expected collisions = N(N-1)/(2·2^32) ≈ 10.5.
# We see 21.9. So ratio = 21.9/10.5 = 2.09.
#
# For a TRULY RANDOM f: {0,1}^32 → {0,1}^32, birthday also gives 10.5.
# (Because N << 2^32, image compression doesn't matter.)
#
# So 2.09× means: H[7] as function of W[0] (with W[1..15]=0) is
# NOT a random 32→32 function. It has MORE collisions than random.

# Verify with a control: random 32→32 function
random.seed(42)
coll_control = 0
ht_ctrl = {}
for _ in range(300000):
    x = random.randint(0, MASK32)
    y = random.randint(0, MASK32)  # truly random output
    if y in ht_ctrl:
        coll_control += 1
    else:
        ht_ctrl[y] = 1

print(f"  Control (random 32→32): {coll_control} collisions (expected ~10.5)")

# And: actual SHA-256 with W[1..15]=0
random.seed(42)
coll_sha = 0
ht_sha = {}
for _ in range(300000):
    w0 = random.randint(0, MASK32)
    h7 = sha256c([w0]+[0]*15)[7]
    if h7 in ht_sha:
        coll_sha += 1
    else:
        ht_sha[h7] = 1

print(f"  SHA-256 W[1..15]=0:     {coll_sha} collisions")
print(f"  Ratio SHA/control:      {coll_sha/max(coll_control,1):.2f}×")

if coll_sha > coll_control * 1.5:
    print(f"\n  ★ SHA-256 with sparse schedule has MORE collisions")
    print(f"    than a random function. This is a STRUCTURAL effect.")
    print(f"    The schedule degeneracy creates non-random clustering in H[7].")
elif coll_sha > coll_control * 1.1:
    print(f"\n  ~ Marginal effect. Needs more data to confirm.")
else:
    print(f"\n  ✗ No effect. SHA-256 ≈ random function even with sparse schedule.")

# ================================================================
# FINAL: What does this mean for birthday attack cost?
# ================================================================

print()
print("=" * 70)
print("IMPLICATIONS")
print("=" * 70)

ratio_measured = coll_sha / max(coll_control, 1)
bits_saved = math.log2(ratio_measured) if ratio_measured > 1 else 0

print(f"""
  Measured effect: {ratio_measured:.2f}× more H[7] collisions for W[1..15]=0
  Bits saved: {bits_saved:.2f} bits on H[7] birthday

  For FULL collision (all 8 words):
    Standard birthday: 2^128
    If ALL 8 words have {ratio_measured:.1f}× compression:
      Birthday = 2^(128 - 8×{bits_saved:.2f}) = 2^{128 - 8*bits_saved:.1f}
    If ONLY H[7] compressed:
      Birthday = 2^(128 - {bits_saved:.2f}) = 2^{128 - bits_saved:.1f}

  Connection to GBP (Generalized Birthday):
    If we can treat each word as a separate list,
    k-list birthday could further reduce cost.
    Wagner k-tree for k=8 lists: O(2^(256/(1+log2(8)))) = O(2^{256/4:.0f})
    But: SHA-256 words are NOT independent lists.
    Correlation between words limits GBP applicability.
""")
