#!/usr/bin/env python3
"""
Direction 3: Partial c-world.

From methodology section 214: fixing ALL 448 carry bits → P(self-consistency)=2^{-77}.
But what if we fix carry for only a FEW rounds?

Approach:
1. For "Hi" (16 bits), compute the ACTUAL carry vectors for all possible W[0] values
2. Group W[0] values by their carry pattern for rounds 0-N
3. Measure: how many distinct carry patterns exist?
4. If few patterns → can enumerate them and solve Q-system for each

Key question: how fast does carry diversity grow with rounds?
- Round 0: carry depends on W[0] + constants → bounded
- Round 1: carry depends on a[1], e[1] which depend on round 0 carry → compounds
- Round N: exponential growth expected

This experiment measures the ACTUAL carry diversity, not theoretical.
"""

import struct, time, hashlib
from collections import Counter

MOD = 2**32
IV = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
      0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]
K = [0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
     0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
     0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
     0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
     0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
     0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
     0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
     0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2]

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&0xFFFFFFFF
def fS0(x): return rotr(x,2)^rotr(x,13)^rotr(x,22)
def fS1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def fCh(e,f,g): return ((e&f)^((~e)&g))&0xFFFFFFFF
def fMaj(a,b,c): return (a&b)^(a&c)^(b&c)
def fs0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def fs1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def add32(*a):
    s=0
    for x in a: s=(s+x)%MOD
    return s

def get_carry_vector(a, b):
    """Return 31-bit carry vector for a + b mod 2^32."""
    carry = 0
    cvec = 0
    for i in range(32):
        ai = (a >> i) & 1
        bi = (b >> i) & 1
        carry = (ai & bi) | (ai & carry) | (bi & carry)
        if i < 31:
            cvec |= (carry << i)
    return cvec

def sha256_with_carry(msg_bytes, max_rounds=64):
    """Compute SHA-256 and return carry vectors for each addition in each round."""
    m = bytearray(msg_bytes); ml = len(msg_bytes)*8; m.append(0x80)
    while len(m)%64!=56: m.append(0)
    m += struct.pack('>Q', ml)
    W = list(struct.unpack('>16I', m[:64]))

    # Schedule carries
    sched_carries = []
    for i in range(16, max_rounds):
        s1 = fs1(W[i-2]); w7 = W[i-7]; s0 = fs0(W[i-15]); w16 = W[i-16]
        # W[i] = s1 + w7 + s0 + w16 (3 additions)
        tmp1 = (s1 + w7) % MOD
        c1 = get_carry_vector(s1, w7)
        tmp2 = (tmp1 + s0) % MOD
        c2 = get_carry_vector(tmp1, s0)
        tmp3 = (tmp2 + w16) % MOD
        c3 = get_carry_vector(tmp2, w16)
        sched_carries.append((c1, c2, c3))
        W.append(tmp3)

    a,b,c,d,e,f,g,h = IV
    round_carries = []

    for r in range(max_rounds):
        # T1 = h + S1(e) + Ch(e,f,g) + K[r] + W[r]  (4 additions)
        s1e = fS1(e); che = fCh(e,f,g)
        t1_1 = (h + s1e) % MOD;          c1 = get_carry_vector(h, s1e)
        t1_2 = (t1_1 + che) % MOD;       c2 = get_carry_vector(t1_1, che)
        t1_3 = (t1_2 + K[r]) % MOD;      c3 = get_carry_vector(t1_2, K[r])
        T1   = (t1_3 + W[r]) % MOD;      c4 = get_carry_vector(t1_3, W[r])

        # T2 = S0(a) + Maj(a,b,c)  (1 addition)
        s0a = fS0(a); maja = fMaj(a,b,c)
        T2 = (s0a + maja) % MOD;         c5 = get_carry_vector(s0a, maja)

        # a' = T1 + T2
        new_a = (T1 + T2) % MOD;         c6 = get_carry_vector(T1, T2)
        # e' = d + T1
        new_e = (d + T1) % MOD;          c7 = get_carry_vector(d, T1)

        round_carries.append((c1, c2, c3, c4, c5, c6, c7))

        h,g,f = g,f,e; e = new_e
        d,c_,b = c,b,a; a = new_a
        c = c_  # rename to avoid shadow

    return round_carries, sched_carries, W


# ============================================================
print("="*60)
print("DIR 3: Partial c-world — carry diversity analysis")
print("="*60)

# For "Hi" (2 bytes), enumerate all 2^16 possible messages
# (upper 16 bits of W[0] are free: byte0 << 24 | byte1 << 16)

# ASCII constraint: 0x20-0x7E for each byte
# That's 95 * 95 = 9025 messages, but let's do all 2^16 for completeness

print("\nEnumerating carry patterns for 2-byte ASCII messages...")
print("(95 × 95 = 9025 ASCII printable messages)\n")

t0 = time.time()

# Track carry fingerprints per round depth
# Fingerprint = tuple of all carry vectors up to round R
carry_fps_by_depth = {}  # depth -> set of fingerprints

N_ROUNDS = 8  # analyze first 8 rounds (manageable)
all_round_carries = []

count = 0
for b0 in range(0x20, 0x7F):  # printable ASCII
    for b1 in range(0x20, 0x7F):
        msg = bytes([b0, b1])
        rc, sc, W = sha256_with_carry(msg, N_ROUNDS)
        all_round_carries.append(rc)
        count += 1

t_enum = time.time() - t0
print(f"Enumerated {count} messages in {t_enum:.2f}s")

# Analyze carry diversity at each round depth
print(f"\n{'Depth':>5s} {'#Patterns':>10s} {'MaxGroup':>10s} {'AvgGroup':>10s} {'Bits':>6s}")
print(f"{'─'*5} {'─'*10} {'─'*10} {'─'*10} {'─'*6}")

for depth in range(1, N_ROUNDS + 1):
    # Fingerprint = all carry vectors from round 0 to depth-1
    fps = Counter()
    for rc in all_round_carries:
        fp = tuple(rc[:depth])  # tuple of carry-vector tuples
        fps[fp] += 1

    n_patterns = len(fps)
    max_group = max(fps.values())
    avg_group = count / n_patterns
    bits = n_patterns.bit_length() - 1  # approximate log2

    print(f"{depth:5d} {n_patterns:10d} {max_group:10d} {avg_group:10.1f} {bits:6d}")

# Deeper analysis: round 0 carry alone
print(f"\n--- Round 0 carry detail ---")
r0_fps = Counter()
for rc in all_round_carries:
    # Round 0: 7 carry vectors, each 31 bits
    r0_fps[rc[0]] += 1

print(f"Round 0 distinct carry patterns: {len(r0_fps)}")
print(f"Top 5 most common:")
for fp, cnt in r0_fps.most_common(5):
    print(f"  {cnt} messages share this pattern")

# How many carry bits are actually FREE in round 0?
# Round 0: inputs are IV (fixed) + W[0] (partially free)
# T1 = IV[7] + Σ₁(IV[4]) + Ch(IV[4],5,6) + K[0] + W[0]
# The additions with constants: carry depends on W[0][16:31] (free bits)
# Lower 16 bits of W[0] = 0x8000 (fixed) → lower carry bits are determined

print(f"\n--- Carry bit analysis for round 0 ---")
# Check which carry bits in T1 computation are constant across all messages
carry_bits_by_pos = {}
for add_idx in range(7):  # 7 additions per round
    for bit in range(31):
        values = set()
        for rc in all_round_carries:
            values.add((rc[0][add_idx] >> bit) & 1)
        carry_bits_by_pos[(add_idx, bit)] = len(values)

n_fixed = sum(1 for v in carry_bits_by_pos.values() if v == 1)
n_free = sum(1 for v in carry_bits_by_pos.values() if v == 2)
print(f"Round 0 carry bits: {n_fixed} fixed + {n_free} free = {n_fixed + n_free} total")
print(f"Free bits per addition:")
for add_idx in range(7):
    free = sum(1 for bit in range(31) if carry_bits_by_pos.get((add_idx, bit), 0) == 2)
    fixed = 31 - free
    names = ['h+S1(e)', 'tmp+Ch', 'tmp+K[0]', 'tmp+W[0]', 'S0(a)+Maj', 'T1+T2', 'd+T1']
    print(f"  Add {add_idx} ({names[add_idx]:>10s}): {free:2d} free / {31:2d} total")

# Key question: can we enumerate carry patterns faster than 2^16?
print(f"\n--- Carry enumeration cost ---")
print(f"Messages: {count}")
print(f"Round 0 patterns: {len(r0_fps)} ({len(r0_fps)/count*100:.1f}% of messages)")

# For depth 1, if there are P0 patterns, and we can solve Q-system per pattern:
# Cost = P0 × cost(Q-solve) + consistency_check
# If P0 << 2^16, this is a win!

for depth in range(1, min(N_ROUNDS+1, 5)):
    fps = Counter()
    for rc in all_round_carries:
        fp = tuple(rc[:depth])
        fps[fp] += 1
    n_p = len(fps)
    # Compare with brute force (2^16 = 65536)
    ratio = n_p / count
    print(f"  Depth {depth}: {n_p} patterns ({ratio:.3f} of search space, "
          f"{'BETTER' if n_p < count else 'WORSE'} than brute force)")

print("\nDone.")
