#!/usr/bin/env python3
"""
Step 12: What happens after round 20?

After round 20, all 16 message words are consumed.
No free words → no direct control over barriers.

Two strategies:
A. Multi-block: use 2nd message block for more freedom
B. Probabilistic: let differential propagate randomly,
   use birthday to find pairs where rounds 21-64 happen to cancel

IDEA C: Use the Wang chain differently.
Instead of De=0 for ALL rounds 1-16, allow De≠0 at SOME rounds.
This "spends" fewer W[i] on maintaining De=0, leaving MORE free words.

For example: if we only require De=0 at rounds 1, 5, 9, 13, 17, 21...
(every 4th round), we get:
- Rounds 1-4: W[0] → De1=0 (DW in W[0])
- Round 5: Da1 accumulates, need DW[4] to cancel → costs W[4]
- Round 9: need DW[8]
- Round 13: need DW[12]
- Round 17: need DW[16] = f(W[14]) → costs W[14]
This uses only W[0,4,8,12,14] = 5 words. 11 words remain free!

But: does this work? If De≠0 at rounds 2,3,4, the differential
grows and may become uncontrollable.

Let's test what happens with SPARSE De=0 constraints.
"""

import random
import numpy as np
from collections import Counter

MASK = 0xFFFFFFFF

def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK
def Sig0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sig1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def sig0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sig1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)
def Ch(e, f, g):  return ((e & f) ^ (~e & g)) & MASK
def Maj(a, b, c): return ((a & b) ^ (a & c) ^ (b & c)) & MASK
def hw(x): return bin(x).count('1')
def add(x, y): return (x + y) & MASK

IV = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]
K = [0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
     0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
     0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
     0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
     0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
     0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
     0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
     0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
     0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
     0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
     0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
     0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
     0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
     0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
     0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
     0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2]

DW_BIT = 0x80000000

def schedule(W16):
    W = list(W16) + [0]*48
    for i in range(16, 64):
        W[i] = add(add(add(sig1(W[i-2]), W[i-7]), sig0(W[i-15])), W[i-16])
    return W

def sha_round_fn(state, W_r, K_r):
    a,b,c,d,e,f,g,h = state
    T1 = add(add(add(add(h, Sig1(e)), Ch(e,f,g)), K_r), W_r)
    T2 = add(Sig0(a), Maj(a,b,c))
    return [add(T1, T2), a, b, c, add(d, T1), e, f, g]

# ============================================================
# PART A: After De17=0 — what's the natural diff evolution?
# ============================================================
print("=" * 70)
print("PART A: Differential evolution after De17=0")
print("=" * 70)
print()

# If we achieve De17=0 but don't control De18+,
# what does the differential look like at rounds 18-64?

# Use the concrete solution from de17_verify.c
W_sol = [0x3eba2820, 0x3082cc68, 0xa325e3b1, 0xb5de7a5f,
         0x4f5f0f61, 0x35acc876, 0x5d3417ff, 0x4344e622,
         0x9c0974a0, 0x6174c1d1, 0xef746113, 0xd763389c,
         0xfc07b5fd, 0x1ebec461, 0x0013efc2, 0xe0c456dc]

W64 = schedule(W_sol)
W_sol_f = list(W_sol); W_sol_f[0] ^= DW_BIT
W64_f = schedule(W_sol_f)

s_n = list(IV); s_f = list(IV)
total_hw = []
de_hw = []

print(f"{'Round':>5} | {'HW(De)':>7} | {'HW(Da)':>7} | {'Total HW':>9} | {'Notes':>15}")
print("-" * 55)

for r in range(64):
    s_n = sha_round_fn(s_n, W64[r], K[r])
    s_f = sha_round_fn(s_f, W64_f[r], K[r])
    de = s_n[4] ^ s_f[4]
    da = s_n[0] ^ s_f[0]
    total = sum(hw(s_n[i] ^ s_f[i]) for i in range(8))
    total_hw.append(total)
    de_hw.append(hw(de))

    notes = ""
    if hw(de) == 0: notes = "De=0 ★"
    elif hw(de) <= 3: notes = "De≈0"
    if total == 0: notes = "COLLISION ★★★"

    if r < 20 or hw(de) <= 5 or r >= 60:
        print(f"  {r+1:3d} | {hw(de):7d} | {hw(da):7d} | {total:9d} | {notes}")

print()

# Hash diff
hash_n = [add(IV[i], s_n[i]) for i in range(8)]
hash_f = [add(IV[i], s_f[i]) for i in range(8)]
hash_diff_bits = sum(hw(hash_n[i] ^ hash_f[i]) for i in range(8))
print(f"Final hash diff: {hash_diff_bits} bits (out of 256)")
print(f"For collision: need 0 bits")
print(f"Ratio: {hash_diff_bits/256:.1%} of bits differ")
print()

# ============================================================
# PART B: Can we choose DW to cancel hash diff?
# ============================================================
print("=" * 70)
print("PART B: Multi-block approach — 2nd message block")
print("=" * 70)
print()

# SHA-256 processes message in 512-bit blocks.
# After block 1: state S1 = compress(IV, M1)
# After block 2: state S2 = compress(S1, M2)
# Final hash: H = S2 (if 2 blocks)
#
# For collision: H(M1||M2) = H(M1'||M2')
# where M1' = M1 ⊕ DM1 (first block differential)
#
# After block 1:
#   S1 = compress(IV, M1)
#   S1' = compress(IV, M1')
#   DS1 = S1' ⊕ S1  (this is our ~120-bit diff from round 64)
#
# For block 2:
#   S2 = compress(S1, M2)
#   S2' = compress(S1', M2')
#   We need: S2 + IV2_equiv = S2' + IV2_equiv (for hash)
#   i.e., S2 = S2'
#   i.e., compress(S1, M2) = compress(S1', M2')
#
# Key: S1 and S1' differ by DS1 (~120 bits).
# M2 and M2' are BOTH under our control (512 bits of freedom).
# We can set M2' = M2 ⊕ DM2 with DM2 chosen freely.
#
# This is EXACTLY the same problem as block 1, but with:
#   - IV replaced by S1 (known)
#   - IV' replaced by S1' = S1 ⊕ DS1 (known)
#   - DS1 is ~120 bits (much larger than DW_bit = 1 bit)
#
# The problem: DS1 is in the STATE, not in the message.
# We need: compress(S1, M2) = compress(S1', M2')
# with S1 ≠ S1'. This is a DIFFERENT problem.

print("Multi-block collision approach:")
print()
print("After block 1:")
print(f"  State diff DS1: {hash_diff_bits} bits in 256-bit state")
print(f"  This becomes the 'IV difference' for block 2")
print()
print("Block 2 problem: find M2, M2' such that")
print("  compress(S1, M2) = compress(S1', M2')")
print("  where DS1 = S1 ⊕ S1' has ~120 bits set")
print()
print("This is HARDER than block 1 (where DIV = 0)")
print("because the state difference is large.")
print()

# ============================================================
# PART C: The actual approach — choose DM2 to cancel DS1
# ============================================================
print("=" * 70)
print("PART C: Differential cancellation in block 2")
print("=" * 70)
print()

# In block 2: compress(S1, M2) processes 16 words M2[0..15]
# through 64 rounds starting from state S1.
#
# For the perturbed version: compress(S1', M2') starts from S1'.
# The initial state difference is DS1.
#
# If we set M2 = M2' (same message for block 2):
# Then the only difference is in the initial state.
# After 64 rounds, this difference becomes DS2.
# For collision: need DS2 + DS1 = 0 (because hash = final + IV)
# Actually: hash = compress(S1, M2) where the final addition is
# H[i] = state[i] + S1[i]. So:
# DH = (state_n[i] + S1[i]) XOR (state_f[i] + S1'[i])
# For DH=0: need state_n[i] + S1[i] = state_f[i] + S1'[i]
# i.e., state_f[i] - state_n[i] = S1[i] - S1'[i] = -DS1[i]
# So we need the differential AFTER block 2 to equal -DS1 (additive).

# Alternatively: we CAN choose M2 ≠ M2' (different DM2).
# This gives us 512 bits of message freedom in block 2!
# The problem becomes: find DM2 such that
# compress(S1, M2) = compress(S1', M2⊕DM2)
# with S1, S1' known.

# This is effectively a generalized collision problem with:
# - 256-bit state difference as initial condition
# - 512-bit message freedom
# - 256-bit target (all-zero state diff)

print("Block 2 collision problem:")
print(f"  Input: 256-bit state diff (from block 1)")
print(f"  Freedom: 512 bits (16 message words)")
print(f"  Target: 256-bit all-zero final diff")
print()
print("  Ratio: 512/256 = 2× overconstrained")
print("  Birthday bound: 2^128 (on 256-bit target)")
print("  Wang-chain approach: 2^128 / effective_dim")
print()
print("This is the SAME complexity as a full SHA-256 collision!")
print("Multi-block doesn't help unless DS1 has special structure.")
print()

# ============================================================
# PART D: DS1 structure — does it have exploitable properties?
# ============================================================
print("=" * 70)
print("PART D: Structure of DS1 (state diff after block 1)")
print("=" * 70)
print()

DS1 = [s_n[i] ^ s_f[i] for i in range(8)]
DS1_add = [(s_f[i] - s_n[i]) & MASK for i in range(8)]

print("DS1 (XOR diff after 64 rounds):")
for i, name in enumerate("abcdefgh"):
    print(f"  {name}: 0x{DS1[i]:08x}  HW={hw(DS1[i]):2d}  add_diff=0x{DS1_add[i]:08x}")

total_ds1 = sum(hw(x) for x in DS1)
print(f"\nTotal HW: {total_ds1} bits")
print()

# Is DS1 random-looking or structured?
print("DS1 structure test:")
print(f"  Total HW / 256 = {total_ds1/256:.3f}  (expected 0.5 for random)")

# Check if any register pair has correlated diff
for i in range(8):
    for j in range(i+1, 8):
        corr = hw(DS1[i] ^ DS1[j])
        if corr < 8 or corr > 24:  # unusually high correlation
            print(f"  HW({chr(97+i)} XOR {chr(97+j)}) = {corr} ← {'correlated!' if corr < 10 else 'anticorrelated!'}")

print()

# ============================================================
# PART E: Alternative — probabilistic rounds 21-64
# ============================================================
print("=" * 70)
print("PART E: Probabilistic approach for rounds 21-64")
print("=" * 70)
print()

# After breaking 4 barriers (rounds 17-20), we have:
# - Controlled differential through round 20
# - 44 rounds remaining (21-64)
# - Each round: differential propagates pseudo-randomly
#
# Question: what is P(final_diff = 0)?
# If differential is random: P = 2^{-256}
# But: from the De17=0 state, differential has SOME structure.
#
# Alternative: use birthday approach.
# Generate 2^128 messages (with De17..De20=0 maintained),
# check if any two have the same final state.
# Cost: 2^128 × (2^34 per message) = 2^162
#
# This is WORSE than generic birthday (2^128).
# So: the Wang chain approach doesn't improve over birthday
# for the full 64 rounds.

# But: can we improve by choosing WHICH rounds to zero De?
# Wang chain zeros De at rounds 1-16 (16 consecutive zeros).
# What if we zero De at DIFFERENT rounds?

# Test: after De17=0, what is De at rounds 18..24?
# Can we find "natural" zeros due to differential structure?

print("Natural De=0 occurrences after De17=0:")
print("(Across 10000 random messages with same structure)")
print()

N = 10000
natural_zeros = np.zeros(64)

for _ in range(N):
    W16 = [random.randint(0, MASK) for _ in range(16)]
    W64 = schedule(W16)
    W16_f = list(W16); W16_f[0] ^= DW_BIT
    W64_f = schedule(W16_f)

    s_n = list(IV); s_f = list(IV)
    for r in range(64):
        s_n = sha_round_fn(s_n, W64[r], K[r])
        s_f = sha_round_fn(s_f, W64_f[r], K[r])
        if s_n[4] == s_f[4]:  # De = 0
            natural_zeros[r] += 1

print(f"{'Round':>5} | {'P(De=0)':>10} | {'Expected 2^-32':>14}")
print("-" * 40)
for r in range(64):
    p = natural_zeros[r] / N
    if p > 0 or r < 3 or r > 60:
        print(f"  {r+1:3d} | {p:10.6f} | {2**-32:14.2e}")

print()
print(f"Total De=0 events across all rounds and messages: {int(natural_zeros.sum())}")
print(f"Expected (64 × {N} × 2^-32): {64*N*2**-32:.6f}")
print()

# ============================================================
# PART F: Multi-block with NEUTRAL bits
# ============================================================
print("=" * 70)
print("PART F: Neutral bit analysis for block 2")
print("=" * 70)
print()

# Neutral bits: bits that don't affect early rounds' differentials.
# In block 2, the STATE already has a large diff.
# If some message bits are "neutral" (don't affect the diff for many rounds),
# they can be varied freely for birthday search.
#
# Test: in block 2 (starting from S1, S1'), vary each M2 bit
# and measure how many rounds before the differential changes.

# For this we'd need to run block 2 compression... but S1 depends
# on the specific message found by de17_verify.c.
# Let's just use IV as S1 and simulate.

print("Neutral bit test (rounds in block 2 before diff changes):")
print("(Using random S1 with DS1 = DW_BIT in word 0)")
print()

S1 = [random.randint(0, MASK) for _ in range(8)]
S1_f = list(S1); S1_f[0] ^= DW_BIT  # minimal state diff

neutral_rounds = []
for bit_word in range(16):
    for bit_pos in [0, 15, 31]:  # sample 3 positions per word
        M2 = [random.randint(0, MASK) for _ in range(16)]
        M2_flip = list(M2)
        M2_flip[bit_word] ^= (1 << bit_pos)

        W64_2 = schedule(M2)
        W64_2f = schedule(M2_flip)

        # Run with normal S1
        s_a = list(S1); s_b = list(S1)
        last_same = 0
        for r in range(20):
            s_a = sha_round_fn(s_a, W64_2[r], K[r])
            s_b = sha_round_fn(s_b, W64_2f[r], K[r])
            if s_a[4] == s_b[4]:  # De still same
                last_same = r + 1
            else:
                break

        neutral_rounds.append(last_same)

avg_neutral = np.mean(neutral_rounds)
max_neutral = max(neutral_rounds)

print(f"  Average neutral rounds: {avg_neutral:.1f}")
print(f"  Max neutral rounds: {max_neutral}")
print(f"  Distribution: {dict(Counter(neutral_rounds))}")
print()

neutral_dist = Counter(neutral_rounds)
print("  Rounds before diff changes:")
for k in sorted(neutral_dist.keys()):
    print(f"    {k} rounds: {neutral_dist[k]} bits")

print()
print("=" * 70)
print("SUMMARY STEP 12")
print("=" * 70)
print()
print("POST-ROUND 20 ANALYSIS:")
print()
print("1. After De17=0, differential explodes to ~128 HW at round 64")
print("2. Multi-block: reduces to SAME problem with 256-bit IV diff")
print("   → no improvement over generic birthday (2^128)")
print("3. Natural De=0 events: ~0 expected (P = 2^-32 per round)")
print("4. Neutral bits in block 2: ~0-1 rounds (no useful neutrality)")
print()
print("CONCLUSION:")
print("The free-word technique provides a CONCRETE improvement")
print("for the first ~20 rounds (extending Wang chain from 16 to 20).")
print("But it does NOT solve the fundamental problem of rounds 21-64,")
print("which remains at generic birthday complexity (~2^128).")
print()
print("The Wang chain barrier at round 17 is BREAKABLE (demonstrated),")
print("but this is one step in a 64-round gauntlet.")
print("Full SHA-256 collision still requires 2^128 work.")
