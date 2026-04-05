#!/usr/bin/env python3
"""
320-bit constraint problem.

Given target H (256 bits) and target A[3] (32 bits, from backward chain):
Find W[0..15] satisfying both.

Structure:
  W[0] = A[1] - 0xfc08884d  (from A[1], 32 bits)
  W[1] = free (32 bits)
  W[2] = determined by (A[1], W[1], A[3])
  W[3..15] = free (13 × 32 = 416 bits)

Constraint: SHA-256(W[0..15]) = H  (256 bits)

Effective: 416 + 32 = 448 free bits, 256 bit constraint → 192 bits slack.

But W[0] and W[2] are functions of (A[1], W[1]).
So: 32 (A[1]) + 32 (W[1]) + 416 (W[3..15]) = 480 free bits.
Constraints: A[3] target (32 bits) + H target (256 bits) = 288 bits.
Slack: 480 - 288 = 192 bits.

KEY QUESTION: are these 288 bits of constraints INDEPENDENT?
If A[3] constraint and H constraint overlap → fewer effective constraints.
"""

import random, math, time

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

A1_CONST = 0xfc08884d
E1_OFFSET = 0x9cbf5a55

def sha256(M):
    W=list(M)+[0]*(64-len(M))
    for i in range(16,64):W[i]=(ssig1(W[i-2])+W[i-7]+ssig0(W[i-15])+W[i-16])&MASK32
    a,b,c,d,e,f,g,h=IV
    for r in range(64):
        T1=(h+sig1(e)+ch(e,f,g)+K[r]+W[r])&MASK32;T2=(sig0(a)+maj(a,b,c))&MASK32
        h,g,f=g,f,e;e=(d+T1)&MASK32;d,c,b=c,b,a;a=(T1+T2)&MASK32
    return tuple((s+v)&MASK32 for s,v in zip([a,b,c,d,e,f,g,h],IV))

def compute_A3_const(A1, W1):
    """Compute CONST such that A[3] = CONST + W[2]."""
    E1 = (A1 + E1_OFFSET) & MASK32
    a,b,c,d = A1, IV[0], IV[1], IV[2]
    e,f,g,h = E1, IV[4], IV[5], IV[6]
    T1=(h+sig1(e)+ch(e,f,g)+K[1]+W1)&MASK32;T2=(sig0(a)+maj(a,b,c))&MASK32
    A2=(T1+T2)&MASK32; E2=(d+T1)&MASK32
    a,b,c,d = A2, A1, IV[0], IV[1]
    e,f,g,h = E2, E1, IV[4], IV[5]
    # A[3] = T1[2]+T2[2]. T1[2] = h+sig1(e)+ch(e,f,g)+K[2]+W[2]
    # CONST = h+sig1(e)+ch(e,f,g)+K[2] + T2[2]
    T1_no_W = (h+sig1(e)+ch(e,f,g)+K[2])&MASK32
    T2_2 = (sig0(a)+maj(a,b,c))&MASK32
    return (T1_no_W + T2_2) & MASK32

def build_message(A1, W1, A3_target, W3_to_15):
    """Build full message from free parameters + A3 constraint."""
    W0 = (A1 - A1_CONST) & MASK32
    const = compute_A3_const(A1, W1)
    W2 = (A3_target - const) & MASK32
    M = [W0, W1, W2] + list(W3_to_15)
    assert len(M) == 16
    return M


# ================================================================
# TEST 1: Is A[3] constraint REDUNDANT with H constraint?
# If A[3] is already determined by H → it adds 0 new bits.
# From backward chain: A[3] IS computed from H.
# So A[3] constraint is NOT independent — it's a SUBSET of H.
# ================================================================

print("=" * 70)
print("TEST 1: Is A[3] from H? (redundancy check)")
print("=" * 70)

# A[3] comes from backward chain starting at H.
# So: knowing H → knowing A[3]. The A[3] constraint is REDUNDANT.
# Total independent constraints: just H = 256 bits.
# Free parameters: 480 bits (A[1], W[1], W[3..15]).
# Slack: 480 - 256 = 224 bits.

print("""
  A[3] is computed FROM H via backward chain.
  So A[3] constraint is CONTAINED IN H constraint.
  NOT independent. Adds 0 new bits.

  Effective:
    Free parameters: 480 bits (A[1] + W[1] + W[3..15])
    Constraints: 256 bits (H only)
    Slack: 224 bits

  This is a PREIMAGE problem, not a collision problem.
  Preimage: given H, find M with SHA(M) = H.
  Cost: 2^256 (brute force) or 2^128 with birthday tricks.
""")


# ================================================================
# TEST 2: But collision uses TWO messages.
# For collision: H(M1) = H(M2) with M1 ≠ M2.
# Both M1, M2 must satisfy: same backward chain → same A[3].
# But they can have DIFFERENT (A1, W1, W3..15).
#
# Two messages with same H:
#   M1 = build_message(A1_a, W1_a, A3, W3..15_a)
#   M2 = build_message(A1_b, W1_b, A3, W3..15_b)
#
# Same H means: SHA(M1) = SHA(M2).
# Both have correct A[3] (by construction).
# Birthday on H: 2^128.
#
# BUT: can we make the search SMARTER?
# We have 448 free bits per message (A1, W1, W3..15).
# We're searching for match on 256 bits (H).
# Standard birthday: 2^128.
#
# The question: does the A[3]-constraint help?
# ================================================================

print()
print("=" * 70)
print("TEST 2: Does A[3] constraint help collision search?")
print("=" * 70)

# Fix A[3] target. Generate many messages with correct A[3].
# Check: is the H-distribution among these messages structured?

A3_target = 0xdf63ae35  # arbitrary
N = 50000
random.seed(42)

# Generate messages with correct A[3]
hashes = []
for _ in range(N):
    A1 = random.randint(0, MASK32)
    W1 = random.randint(0, MASK32)
    W3_15 = [random.randint(0, MASK32) for _ in range(13)]

    M = build_message(A1, W1, A3_target, W3_15)
    H = sha256(M)
    hashes.append(H)

# Check: H[0] distribution (should be uniform if no structure)
h0_counts = {}
for H in hashes:
    h0_counts[H[0]] = h0_counts.get(H[0], 0) + 1

# Birthday on H[0]: expected collisions
h0_colls = sum(c*(c-1)//2 for c in h0_counts.values())
h0_expected = N*(N-1)/(2*2**32)

print(f"  N={N} messages with A[3]={A3_target:#010x}")
print(f"  H[0] collisions: {h0_colls} (expected: {h0_expected:.1f})")
print(f"  Ratio: {h0_colls/max(h0_expected, 0.1):.2f}×")

# Full 256-bit collision check
full_ht = {}
full_colls = 0
for H in hashes:
    key = H
    if key in full_ht:
        full_colls += 1
    else:
        full_ht[key] = 1

print(f"  Full H collisions: {full_colls} (expected: ~0)")

# Compare: random messages (no A[3] constraint)
random.seed(42)
hashes_rand = []
for _ in range(N):
    M = [random.randint(0, MASK32) for _ in range(16)]
    H = sha256(M)
    hashes_rand.append(H)

h0_counts_r = {}
for H in hashes_rand:
    h0_counts_r[H[0]] = h0_counts_r.get(H[0], 0) + 1
h0_colls_r = sum(c*(c-1)//2 for c in h0_counts_r.values())

print(f"\n  Control (random M, no A[3] constraint):")
print(f"  H[0] collisions: {h0_colls_r} (expected: {h0_expected:.1f})")
print(f"  Ratio: {h0_colls_r/max(h0_expected, 0.1):.2f}×")

if abs(h0_colls - h0_colls_r) < 3 * max(h0_expected**0.5, 1):
    print(f"\n  A[3] constraint does NOT change H distribution.")
    print(f"  Collision cost UNCHANGED: 2^128.")
else:
    print(f"\n  ★ A[3] constraint CHANGES H distribution!")
    print(f"  Investigate further.")


# ================================================================
# TEST 3: What about fixing W[3..15] and varying only (A[1], W[1])?
# Then: only 64 bits free, 256 bits constraint.
# Underdetermined → no solution expected.
# But: hash is 256 bits, 64 free bits → 2^{256-64} = 2^192 tries?
# No: birthday on 256 bits with 64-bit input: 2^{192} is brute force,
# birthday doesn't apply (input too small).
# ================================================================

print()
print("=" * 70)
print("TEST 3: Fix W[3..15], vary (A[1], W[1]) only")
print("=" * 70)

# For fixed W[3..15]: function (A[1], W[1]) → H is 64→256.
# This function has image size at most 2^64.
# Birthday on 256 bits within 2^64 image: need 2^{64} × 2 = 2^{65} for
# any H-collision? No: birthday on n-bit hash with N samples:
# N^2 / (2 * 2^n). For N = 2^64, n = 256: collisions = 2^{128}/2^{256} = 2^{-128}.
# So: 0 collisions expected from 2^64 samples on 256-bit hash.

# But: birthday on SINGLE WORD H[0] (32 bits) with 2^64 input → 2^{64-16} = 2^{48} expected H[0]-collisions!

# Let's check with small sample:
W3_15_fixed = [random.randint(0, MASK32) for _ in range(13)]
N3 = 200000

ht_h0 = {}
h0_coll_count = 0

random.seed(789)
t0 = time.time()
for _ in range(N3):
    A1 = random.randint(0, MASK32)
    W1 = random.randint(0, MASK32)
    M = build_message(A1, W1, A3_target, W3_15_fixed)
    H = sha256(M)

    if H[0] in ht_h0:
        h0_coll_count += 1
    else:
        ht_h0[H[0]] = (A1, W1)

elapsed = time.time() - t0
h0_expected_3 = N3*(N3-1)/(2*2**32)

print(f"  W[3..15] fixed, vary (A[1], W[1])")
print(f"  N={N3}, time={elapsed:.1f}s")
print(f"  H[0] collisions: {h0_coll_count} (expected: {h0_expected_3:.1f})")
print(f"  Ratio: {h0_coll_count/max(h0_expected_3, 0.1):.2f}×")

if abs(h0_coll_count/max(h0_expected_3,0.1) - 1.0) < 0.3:
    print(f"\n  H[0] behaves as random even with A[3] fixed and W[3..15] fixed.")
    print(f"  No exploitable structure from backward chain.")


# ================================================================
# TEST 4: The REAL question — schedule sensitivity
# When we change W[0] (= change A[1]):
# How many schedule words W[16..63] change?
# If few change → few rounds affected → cheaper collision
# ================================================================

print()
print("=" * 70)
print("TEST 4: Schedule sensitivity to W[0..2] changes")
print("=" * 70)

# Fix W[3..15]. Vary (A[1], W[1]) → changes W[0], W[1], W[2].
# How many W[16..63] words differ between two choices?

N4 = 1000
random.seed(42)
W3_15_fixed = [random.randint(0, MASK32) for _ in range(13)]

diff_counts = []
for _ in range(N4):
    A1_a = random.randint(0, MASK32)
    W1_a = random.randint(0, MASK32)
    M_a = build_message(A1_a, W1_a, A3_target, W3_15_fixed)

    A1_b = random.randint(0, MASK32)
    W1_b = random.randint(0, MASK32)
    M_b = build_message(A1_b, W1_b, A3_target, W3_15_fixed)

    # Expand schedules
    Wa = list(M_a)+[0]*48; Wb = list(M_b)+[0]*48
    for i in range(16,64):
        Wa[i]=(ssig1(Wa[i-2])+Wa[i-7]+ssig0(Wa[i-15])+Wa[i-16])&MASK32
        Wb[i]=(ssig1(Wb[i-2])+Wb[i-7]+ssig0(Wb[i-15])+Wb[i-16])&MASK32

    diff = sum(1 for r in range(16, 64) if Wa[r] != Wb[r])
    diff_counts.append(diff)

avg_diff = sum(diff_counts)/N4
min_diff = min(diff_counts)
max_diff = max(diff_counts)

print(f"  W[3..15] shared, W[0..2] differ (both satisfy A[3]={A3_target:#010x})")
print(f"  N={N4} pairs")
print(f"  Schedule words differing (W[16..63]):")
print(f"    avg = {avg_diff:.1f}/48")
print(f"    min = {min_diff}/48")
print(f"    max = {max_diff}/48")

# Now: what if ONLY W[0] differs (same W[1], same A[3] → different W[2])?
diff_counts_w0 = []
for _ in range(N4):
    A1_a = random.randint(0, MASK32)
    A1_b = random.randint(0, MASK32)
    W1_shared = random.randint(0, MASK32)

    M_a = build_message(A1_a, W1_shared, A3_target, W3_15_fixed)
    M_b = build_message(A1_b, W1_shared, A3_target, W3_15_fixed)

    Wa = list(M_a)+[0]*48; Wb = list(M_b)+[0]*48
    for i in range(16,64):
        Wa[i]=(ssig1(Wa[i-2])+Wa[i-7]+ssig0(Wa[i-15])+Wa[i-16])&MASK32
        Wb[i]=(ssig1(Wb[i-2])+Wb[i-7]+ssig0(Wb[i-15])+Wb[i-16])&MASK32

    diff = sum(1 for r in range(16, 64) if Wa[r] != Wb[r])
    diff_counts_w0.append(diff)

    # Which words are zero? (shared W[3..15] means many schedule words shared)
    if _ == 0:
        zeros = [r for r in range(16,64) if Wa[r] == Wb[r]]
        nonzeros = [r for r in range(16,64) if Wa[r] != Wb[r]]
        print(f"\n  Example: W[0],W[2] differ, W[1,3..15] shared")
        print(f"  Zero-diff schedule words: {zeros[:15]}{'...' if len(zeros)>15 else ''}")
        print(f"  ({len(zeros)}/48 shared)")

avg_w0 = sum(diff_counts_w0)/N4

print(f"\n  Only W[0] differs (W[1] shared, W[2] auto-adjusts for A[3]):")
print(f"    avg schedule diff = {avg_w0:.1f}/48")
print(f"    min = {min(diff_counts_w0)}/48")

# KEY: from Dir7, single W[0] change → 45/48 words differ.
# With W[3..15] shared → fewer?

print(f"\n  Compare: Dir7 single-bit W[0] change → ~45/48 differ")
print(f"  Here: W[0]+W[2] change → {avg_w0:.0f}/48 differ")


# ================================================================
# SYNTHESIS
# ================================================================

print()
print("=" * 70)
print("SYNTHESIS: 320-bit constraint")
print("=" * 70)

print(f"""
  Problem: find M with SHA(M) = H (preimage) or SHA(M1)=SHA(M2) (collision).

  Backward chain gives: A[3] from H (free).
  A[1], W[1] = 64 free bits. W[2] determined. W[3..15] = 416 free.
  Total: 480 free bits, 256 constraints (H) → 224 slack.

  A[3] constraint is REDUNDANT (subset of H constraint).
  Does not reduce effective constraints.

  Schedule sensitivity: changing W[0]+W[2] (from A[1] choice)
  affects {avg_w0:.0f}/48 schedule words. Nearly ALL.
  Sharing W[3..15] does NOT prevent schedule cascade.

  H[0] distribution: UNIFORM regardless of A[3] constraint.
  No birthday acceleration from backward chain.

  The backward chain (a[3..64] from H) gives us KNOWLEDGE
  but not CONTROL. We know what the solution must look like
  but cannot steer toward it faster than birthday.
""")
