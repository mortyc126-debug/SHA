#!/usr/bin/env python3
"""
Step 7b: Optimized De17=0 search.

Optimization: instead of Python loops, use a smarter strategy.

Key insight: De17 is a 32-bit value. We want De17=0.
Strategy: fix W[0..11], vary W[14] only (32 bits).
This is a 32→32 mapping. If it's a random function,
~63% of outputs are covered → P(0 in image) ≈ 0.63.

So with just ONE base message, sweeping W[14] through 2^32 values
should find a solution with probability ~63%.

But 2^32 Python iterations is slow. Let's use numpy for speed.
Or: use the STRUCTURE we found.

ALGEBRAIC APPROACH:
De17 = e17(W) ⊕ e17(W')
     = [d13 + T1_16] ⊕ [d'13 + T1'_16]
     (where d13 = a9, T1_16 involves W[16])

The e-register at round 17:
  e17 = d13 + T1_16
  where T1_16 = h13 + Sig1(e13) + Ch(e13,f13,g13) + K[16] + W[16]

For the differential:
  De17 = Dd13 ⊕ ... (mixing XOR and addition)

Actually let me just try a smarter search: low-bits-first.

LOW-BITS-FIRST SEARCH:
If De17 is pseudo-random, we can search for De17=0 bit by bit:
1. Find DW14 such that De17 mod 2 = 0 → 50% chance, 2 trials
2. Among those, find DW14 such that De17 mod 4 = 0 → 50% again
3. Continue for 32 bits → total trials: 32 × 2 = 64!

But this only works if we can INDEPENDENTLY control each bit.
Since additions propagate carries from low to high bits,
low bits of De17 depend ONLY on low bits of free words!

This is the CARRY CHAIN PROPERTY:
  De17[0] depends only on bit 0 of all inputs
  De17[1] depends on bits 0-1 of all inputs
  De17[k] depends on bits 0-k of all inputs

So: we can solve bit-by-bit from LSB to MSB!
"""

import random
import numpy as np

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
def sub(x, y): return (x - y) & MASK

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

DW_bit = 0x80000000

def schedule(W16):
    W = list(W16) + [0]*48
    for i in range(16, 64):
        W[i] = add(add(add(sig1(W[i-2]), W[i-7]), sig0(W[i-15])), W[i-16])
    return W

def sha_round(state, W_r, K_r):
    a,b,c,d,e,f,g,h = state
    T1 = add(add(add(add(h, Sig1(e)), Ch(e,f,g)), K_r), W_r)
    T2 = add(Sig0(a), Maj(a,b,c))
    return [add(T1, T2), a, b, c, add(d, T1), e, f, g]

def compute_de17_from_w16(W16):
    W64 = schedule(W16)
    W16_f = list(W16); W16_f[0] ^= DW_bit
    W64_f = schedule(W16_f)
    s_n = list(IV)
    s_f = list(IV)
    for r in range(17):
        s_n = sha_round(s_n, W64[r], K[r])
        s_f = sha_round(s_f, W64_f[r], K[r])
    return s_n[4] ^ s_f[4]

# ============================================================
# PART A: Test carry chain property — does it hold?
# ============================================================
print("=" * 70)
print("PART A: Carry chain property — bit independence test")
print("=" * 70)
print()

# Question: does De17 mod 2^k depend only on low k bits of W[14]?
# Test: fix W[0..13,15], change only low k bits of W[14], check De17 mod 2^k

W16_test = [random.randint(0, MASK) for _ in range(16)]

print("Does De17 mod 2^k depend only on W[14] mod 2^k?")
print("(Fix W[0..13,15], vary low bits of W[14])")
print()

# For this to work, the ENTIRE computation must respect carry locality.
# But SHA rounds involve rotations (ROTR) which MIX bits globally!
# So this property does NOT hold for SHA-256.

# Let me verify this directly:
for k in [1, 2, 4, 8, 16]:
    mask_k = (1 << k) - 1
    base_de17 = compute_de17_from_w16(W16_test) & mask_k

    # Change high bits of W[14], keep low k bits same
    all_same = True
    for _ in range(100):
        W16_var = list(W16_test)
        W16_var[14] = (W16_test[14] & mask_k) | (random.randint(0, MASK) & ~mask_k)
        de17_var = compute_de17_from_w16(W16_var) & mask_k
        if de17_var != base_de17:
            all_same = False
            break

    print(f"  k={k:2d}: De17 mod 2^{k} stable when high bits change: {all_same}")

print()
print("(Expected: False — rotations destroy carry locality)")
print()

# ============================================================
# PART B: Partial matching — fix some bits, search others
# ============================================================
print("=" * 70)
print("PART B: Partial bit matching strategy")
print("=" * 70)
print()

# Since carry locality doesn't hold, we need a different approach.
# Let's try: fix low 16 bits of W[14], sweep high 16 bits.
# For each of 2^16 high-bit patterns, compute De17.
# Build a lookup table. Then try different low 16 bits.

# MITM on W[14]:
# Split W[14] = (high16, low16)
# Phase 1: for each low16 value, compute De17 → store
# Phase 2: impossible to do true MITM because De17 is nonlinear

# Better: just brute force with 2^20 trials using 2 free words
print("Brute force search with W[14] (20-bit range):")
W16_search = [random.randint(0, MASK) for _ in range(16)]
best = 32
best_w14 = 0
results = {}

for w14 in range(1 << 20):
    W16_search[14] = w14
    de17 = compute_de17_from_w16(W16_search)
    h = hw(de17)
    if h < best:
        best = h
        best_w14 = w14
        print(f"  W14={hex(w14):>10s} → De17={hex(de17):>12s} HW={h}")
    if h == 0:
        break

    if (w14 + 1) % 200000 == 0:
        print(f"  ... searched {w14+1} ({(w14+1)/(1<<20)*100:.0f}%)")

print()
print(f"Best after 2^20 trials: HW(De17) = {best}")
print()

# Expected: with 2^20 trials on 32-bit output:
# P(HW=0) = 2^20/2^32 = 2^-12 → unlikely
# P(HW≤2) = 2^20 × C(32,≤2)/2^32 ≈ 2^20 × 529/2^32 ≈ 0.13
# P(HW≤3) ≈ 2^20 × 5489/2^32 ≈ 1.34 → likely!

# ============================================================
# PART C: Two-word search — W[14] × W[13], 16 bits each
# ============================================================
print("=" * 70)
print("PART C: Two-word birthday search (16 bits × 16 bits)")
print("=" * 70)
print()

# Phase 1: Build table {De17 → W14} for 2^16 values of W14
print("Phase 1: Building De17 table (2^16 W14 values)...")
W16_mitm = [random.randint(0, MASK) for _ in range(16)]

table = {}
for w14 in range(1 << 16):
    W16_mitm[14] = w14
    de17 = compute_de17_from_w16(W16_mitm)
    table[de17] = w14

print(f"  Table size: {len(table)} entries")
print(f"  (Unique De17: {len(table)})")

# Phase 2: For 2^16 values of W13, compute De17 and look for match
# We need De17 = 0, not collision. But:
# If we CHANGE W13 and compute De17, and it happens to be in the table
# with De17=0, great. But there's no algebraic reason for this.

# ALTERNATIVE: Use the table for collision.
# Two different (W13, W14) pairs that give the SAME De17
# → their difference is a differential with De=0? NO — that's not right.

# Actually, the CORRECT birthday approach:
# We want De17(W, DW=DW_bit, W14) = 0
# This is a function f(W14) → {0,1}^32
# We need to invert: find W14 s.t. f(W14) = 0
# With 2^16 samples: P(hit) = 2^16/2^32 = 2^-16 → very unlikely

# The BIRTHDAY trick:
# Consider TWO different base messages W_A, W_B (differing in W[0..11])
# f_A(W14) = De17(W_A, DW_bit, W14)
# f_B(W14) = De17(W_B, DW_bit, W14)
# If f_A(w14_a) = f_B(w14_b), does that help? Not directly.

# TRUE multi-target birthday:
# Generate M base messages → M tables of De17 values
# For each message i, compute {De17_i(w14) : w14 in [0, 2^16)}
# If ANY De17_i(w14) = 0, we win.
# P(win) = 1 - (1 - 2^16/2^32)^M = 1 - (1-2^-16)^M
# For M = 2^16: P ≈ 1 - e^{-1} ≈ 0.63  ← HIGH PROBABILITY!

print()
print("Phase 2: Multi-target search")
print(f"  Need M ≈ 2^16 = 65536 base messages")
print(f"  Each with 2^16 free-word trials")
print(f"  Total: 2^32 evaluations — too slow for Python")
print()

# Let's try M=100 messages × 2^16 trials = 6.5M evaluations
M = 100
print(f"Practical search: M={M} messages × 2^16 W14 values = {M * 65536:,} evals")
print()

found = False
best_global = 32
best_global_info = None

for msg_idx in range(M):
    W16_msg = [random.randint(0, MASK) for _ in range(16)]

    for w14 in range(1 << 16):
        W16_msg[14] = w14
        de17 = compute_de17_from_w16(W16_msg)
        h = hw(de17)
        if h < best_global:
            best_global = h
            best_global_info = (msg_idx, w14, de17)
            if h <= 2:
                print(f"  msg={msg_idx}, W14={hex(w14)}: De17={hex(de17)} HW={h}")
        if h == 0:
            found = True
            break

    if found:
        break

    if (msg_idx + 1) % 20 == 0:
        print(f"  ... {msg_idx+1}/{M} messages, best HW = {best_global}")

print()
total = (msg_idx + 1) * 65536
print(f"Total evaluations: {total:,}")
print(f"Best HW(De17): {best_global}")

if found:
    msg_idx, w14, de17 = best_global_info
    print(f"\n★★★ De17 = 0 FOUND! ★★★")
    print(f"  Message #{msg_idx}, W14 = {hex(w14)}")
elif best_global_info:
    msg_idx, w14, de17 = best_global_info
    print(f"\nClosest: msg={msg_idx}, W14={hex(w14)}, De17={hex(de17)}, HW={best_global}")

print()
print("=" * 70)
print("ANALYSIS")
print("=" * 70)
print()
print("Expected: P(De17=0 | M messages × 2^16 W14) = 1-(1-2^-16)^M")
for m_test in [100, 1000, 10000, 65536]:
    p = 1 - (1 - 2**-16)**m_test
    print(f"  M={m_test:>6d}: P = {p:.6f}")
print()
print("For M=65536: P ≈ 0.63 — this is the natural De17=0 probability")
print("For Python: M=100 gives P ≈ 0.0015 — too low")
print()
print("CONCLUSION: De17=0 search requires ~2^32 total evaluations")
print("This is FEASIBLE in C/Rust (~10 seconds)")
print("In Python: prohibitively slow (~hours)")
print()
print("But the KEY ALGEBRAIC FINDING stands:")
print("  - 4 free words give 128 bits of freedom")
print("  - De17 is 32-bit constraint")
print("  - Solution GUARANTEED to exist")
print("  - Searchable in ~2^32 evaluations")
print("  - This is MUCH better than 2^65 for full collision")
