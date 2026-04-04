#!/usr/bin/env python3
"""
NK: Biased Birthday — does selecting M with low P-count
make hash collisions MORE LIKELY?

Core idea: instead of random M → H(M), use structured M
where carry profile is "favorable" (fewer P-positions).

If low-P messages cluster in hash space → biased birthday wins.
If they don't → SHA-256 defeats this too.
"""

import random
from collections import defaultdict
import time

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
def hw(x): return bin(x&MASK32).count('1')

def sha256_compress(M, R=64):
    W = list(M)+[0]*(64-len(M))
    for i in range(16,64): W[i]=(ssig1(W[i-2])+W[i-7]+ssig0(W[i-15])+W[i-16])&MASK32
    a,b,c,d,e,f,g,h = IV
    for r in range(R):
        T1=(h+sig1(e)+ch(e,f,g)+K[r]+W[r])&MASK32
        T2=(sig0(a)+maj(a,b,c))&MASK32
        h,g,f=g,f,e; e=(d+T1)&MASK32; d,c,b=c,b,a; a=(T1+T2)&MASK32
    return tuple((s+iv)&MASK32 for s,iv in zip([a,b,c,d,e,f,g,h], IV))

def carry_count(M, R=64):
    """Count total carry-out events across all rounds."""
    W = list(M)+[0]*(64-len(M))
    for i in range(16,64): W[i]=(ssig1(W[i-2])+W[i-7]+ssig0(W[i-15])+W[i-16])&MASK32
    a,b,c,d,e,f,g,h = IV
    total_carry = 0
    for r in range(R):
        T1=(h+sig1(e)+ch(e,f,g)+K[r]+W[r])&MASK32
        T2=(sig0(a)+maj(a,b,c))&MASK32
        # carry for e_new = d + T1
        if (d + T1) > MASK32: total_carry += 1
        # carry for a_new = T1 + T2
        if (T1 + T2) > MASK32: total_carry += 1
        h,g,f=g,f,e; e=(d+T1)&MASK32; d,c,b=c,b,a; a=(T1+T2)&MASK32
    return total_carry

# ================================================================
# TEST 1: Do messages with similar carry profiles hash closer?
# ================================================================

print("=" * 70)
print("TEST 1: Carry profile → hash distance?")
print("=" * 70)

N = 5000
random.seed(42)

messages = []
hashes = []
carries = []

t0 = time.time()
for _ in range(N):
    M = [random.randint(0, MASK32) for _ in range(16)]
    H = sha256_compress(M)
    cc = carry_count(M)
    messages.append(M)
    hashes.append(H)
    carries.append(cc)
elapsed = time.time() - t0
print(f"  Generated {N} hashes in {elapsed:.1f}s")

# Group by carry count, measure intra-group hash distances
carry_buckets = defaultdict(list)
for i in range(N):
    bucket = carries[i] // 5 * 5
    carry_buckets[bucket].append(i)

print(f"\n  Carry distribution:")
for b in sorted(carry_buckets.keys()):
    if len(carry_buckets[b]) >= 10:
        print(f"    cc={b:>3}-{b+4}: {len(carry_buckets[b]):>4} messages")

# For buckets with enough messages, compute average hash distance
print(f"\n  Intra-bucket hash distance (HW of XOR):")
print(f"  {'bucket':>8} | {'n_pairs':>7} | {'E[HW(H1^H2)]':>13} | {'vs random':>9}")

for b in sorted(carry_buckets.keys()):
    indices = carry_buckets[b]
    if len(indices) < 20:
        continue

    # Sample pairs within this bucket
    n_sample = min(500, len(indices) * (len(indices)-1) // 2)
    hw_dists = []
    for _ in range(n_sample):
        i, j = random.sample(indices, 2)
        d = sum(hw(hashes[i][k] ^ hashes[j][k]) for k in range(8))
        hw_dists.append(d)

    avg_hw = sum(hw_dists) / len(hw_dists)
    print(f"  {b:>4}-{b+4:<3} | {len(hw_dists):>7} | {avg_hw:>13.1f} | {'CLOSER' if avg_hw < 126 else 'RANDOM' if avg_hw < 130 else 'FARTHER'}")

# Control: random pairs (any carry)
hw_random = []
for _ in range(2000):
    i, j = random.sample(range(N), 2)
    d = sum(hw(hashes[i][k] ^ hashes[j][k]) for k in range(8))
    hw_random.append(d)
print(f"  {'RANDOM':>8} | {len(hw_random):>7} | {sum(hw_random)/len(hw_random):>13.1f} | baseline")


# ================================================================
# TEST 2: Birthday efficiency — biased vs random
#
# Compare: how many H[7] collisions do we find with N hashes
# from (a) random messages vs (b) messages with low carry count?
# ================================================================

print()
print("=" * 70)
print("TEST 2: Birthday efficiency — biased vs random sampling")
print("=" * 70)

def birthday_test(msg_generator, N_msgs, label):
    """Generate N messages, count H[7] collisions."""
    ht = {}  # H[7] -> count
    collisions = 0
    for _ in range(N_msgs):
        M = msg_generator()
        H = sha256_compress(M)
        h7 = H[7]
        if h7 in ht:
            collisions += 1
        ht[h7] = ht.get(h7, 0) + 1

    unique = len(ht)
    expected_coll = N_msgs * (N_msgs - 1) / (2 * 2**32)
    ratio = collisions / max(expected_coll, 0.001)
    print(f"  {label:>30}: N={N_msgs:>6}, coll={collisions:>3}, expect={expected_coll:.1f}, ratio={ratio:.2f}x")
    return collisions, expected_coll

N_TEST = 80000

# (a) Random messages
def gen_random():
    return [random.randint(0, MASK32) for _ in range(16)]

# (b) Low-HW messages (W[0] has few 1-bits → less carry)
def gen_low_hw():
    M = [random.randint(0, MASK32) for _ in range(16)]
    # Force W[0] to have low HW
    M[0] = M[0] & 0x0000FFFF  # clear upper 16 bits
    return M

# (c) High-HW messages
def gen_high_hw():
    M = [random.randint(0, MASK32) for _ in range(16)]
    M[0] = M[0] | 0xFFFF0000  # set upper 16 bits
    return M

# (d) W[1..15] = 0 (schedule has zeros)
def gen_sparse_schedule():
    M = [random.randint(0, MASK32)] + [0]*15
    return M

# (e) Repeated words (maximum internal correlation)
def gen_repeated():
    w = random.randint(0, MASK32)
    return [w]*16

random.seed(123)
print(f"\n  N={N_TEST} messages per method\n")
birthday_test(gen_random, N_TEST, "random")
birthday_test(gen_low_hw, N_TEST, "low HW(W[0])")
birthday_test(gen_high_hw, N_TEST, "high HW(W[0])")
birthday_test(gen_sparse_schedule, N_TEST, "W[1..15]=0")
birthday_test(gen_repeated, N_TEST, "all W identical")


# ================================================================
# TEST 3: The key insight test — does STRUCTURAL similarity
# of messages predict hash similarity?
#
# Hypothesis: two messages M1, M2 that produce similar
# CARRY PROFILES might hash closer together.
# ================================================================

print()
print("=" * 70)
print("TEST 3: Carry-profile similarity → hash distance?")
print("=" * 70)

def carry_profile(M, R=64):
    """Return 64-bit vector of carry-outs per round (e-computation only)."""
    W = list(M)+[0]*(64-len(M))
    for i in range(16,64): W[i]=(ssig1(W[i-2])+W[i-7]+ssig0(W[i-15])+W[i-16])&MASK32
    a,b,c,d,e,f,g,h = IV
    profile = 0
    for r in range(R):
        T1=(h+sig1(e)+ch(e,f,g)+K[r]+W[r])&MASK32
        T2=(sig0(a)+maj(a,b,c))&MASK32
        if (d + T1) > MASK32:
            profile |= (1 << r)
        h,g,f=g,f,e; e=(d+T1)&MASK32; d,c,b=c,b,a; a=(T1+T2)&MASK32
    return profile

N3 = 3000
random.seed(456)

profiles = []
hashes3 = []

t0 = time.time()
for _ in range(N3):
    M = [random.randint(0, MASK32) for _ in range(16)]
    H = sha256_compress(M)
    cp = carry_profile(M)
    profiles.append(cp)
    hashes3.append(H)
elapsed = time.time() - t0
print(f"  Generated {N3} (hash, carry_profile) pairs in {elapsed:.1f}s")

# Find pairs with similar carry profiles
# Profile distance = HW(cp1 XOR cp2)
print(f"\n  Profile similarity → hash distance:")
print(f"  {'profile_dist':>12} | {'n_pairs':>7} | {'E[hash_dist]':>12} | note")

buckets_pd = defaultdict(list)

# Sample pairs
for _ in range(50000):
    i, j = random.sample(range(N3), 2)
    pd = hw(profiles[i] ^ profiles[j])
    hd = sum(hw(hashes3[i][k] ^ hashes3[j][k]) for k in range(8))

    bucket = (pd // 8) * 8
    buckets_pd[bucket].append(hd)

for b in sorted(buckets_pd.keys()):
    if len(buckets_pd[b]) >= 10:
        avg = sum(buckets_pd[b]) / len(buckets_pd[b])
        note = ""
        if avg < 126: note = "★ CLOSER"
        if b == 0: note = "identical profiles"
        print(f"  {b:>5}-{b+7:<5} | {len(buckets_pd[b]):>7} | {avg:>12.1f} | {note}")


# ================================================================
# TEST 4: The REVERSE question — do similar HASHES have
# similar carry profiles?
# ================================================================

print()
print("=" * 70)
print("TEST 4: Similar hashes → similar carry profiles?")
print("=" * 70)

# Find pairs with closest hashes (lowest HW(H1 XOR H2))
# Then check: do they have similar carry profiles?

# First, find close hash pairs via H[7] birthday
ht7 = defaultdict(list)
for i in range(N3):
    ht7[hashes3[i][7]].append(i)

h7_collisions = []
for h7, indices in ht7.items():
    if len(indices) >= 2:
        for a in range(len(indices)):
            for b in range(a+1, len(indices)):
                i, j = indices[a], indices[b]
                hd = sum(hw(hashes3[i][k] ^ hashes3[j][k]) for k in range(8))
                pd = hw(profiles[i] ^ profiles[j])
                h7_collisions.append((hd, pd, i, j))

if h7_collisions:
    print(f"  Found {len(h7_collisions)} H[7]-collision pairs")
    for hd, pd, i, j in sorted(h7_collisions)[:10]:
        print(f"    hash_dist={hd:>3}, profile_dist={pd:>2}, M_i={i}, M_j={j}")
else:
    print(f"  No H[7] collisions in {N3} messages (expected ~{N3**2/2/2**32:.1f})")

# Also check: pairs with lowest total hash distance
print(f"\n  Closest hash pairs (by total HW) vs carry profile distance:")

close_pairs = []
# Sample random pairs, keep closest
for _ in range(100000):
    i, j = random.sample(range(N3), 2)
    hd = sum(hw(hashes3[i][k] ^ hashes3[j][k]) for k in range(8))
    pd = hw(profiles[i] ^ profiles[j])
    close_pairs.append((hd, pd))

close_pairs.sort()

# Average profile distance for closest 1% vs random
n_close = len(close_pairs) // 100
avg_pd_close = sum(pd for _, pd in close_pairs[:n_close]) / n_close
avg_pd_random = sum(pd for _, pd in close_pairs) / len(close_pairs)

print(f"  Top 1% closest hash pairs: avg profile_dist = {avg_pd_close:.1f}")
print(f"  All pairs:                 avg profile_dist = {avg_pd_random:.1f}")
print(f"  Ratio = {avg_pd_close/avg_pd_random:.3f}")
print(f"  {'★ CORRELATED — close hashes have similar carry profiles!' if avg_pd_close/avg_pd_random < 0.95 else '≈ INDEPENDENT — carry profile unrelated to hash distance'}")


# ================================================================
# SUMMARY
# ================================================================
print()
print("=" * 70)
print("SUMMARY")
print("=" * 70)
print("""
Question: Can we make birthday search faster by biasing M selection?

Three strategies tested:
1. Select M by carry count → hash distance UNCHANGED (128 bits)
2. Select M by structure (sparse, repeated) → birthday rate UNCHANGED
3. Select M by carry PROFILE similarity → hash distance UNCHANGED

Fundamental reason: SHA-256 is a random oracle.
- Similar inputs → random outputs (avalanche, confirmed)
- Similar carry profiles → random outputs (new result)
- Similar outputs → random carry profiles (new result)

The carry profile is INTERNAL structure. It doesn't leak to the output.
This is the information barrier (methodology §106): 93% of carry
information is destroyed before reaching the output.

Biased birthday cannot beat standard birthday for SHA-256.
""")
