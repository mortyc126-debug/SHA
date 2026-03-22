#!/usr/bin/env python3
"""
Step 7: Find ACTUAL De17=0 using MITM (Meet-in-the-Middle).

We have 4 free words: W[12], W[13], W[14], W[15] (128 bits freedom).
De17 is 32-bit target. This is massively overconstrained.

MITM approach:
  Group 1: vary W[13] → compute De17 → store in dict
  Group 2: vary W[14] → compute De17 → look for collision where
           De17_group1 = De17_group2 (both give same De17 → their
           difference gives us a way to zero out De17)

Wait, actually we need De17(W, DW) = 0 directly.
With DW = {0: DW_bit, 13: DW13, 14: DW14}:

MITM on De17:
  Phase 1: For 2^16 random DW13 values, compute De17_partial
           (with DW14=0). Store: {De17_partial → DW13}
  Phase 2: For 2^16 random DW14 values, compute De17_partial
           Look for: De17(DW13, DW14) = 0

But De17 is not separable into DW13 and DW14 parts (nonlinear mixing).
So we can't do standard MITM.

Better approach: since we have 128 bits and need 32 bits = 0,
just do random search with 2^16 trials on 64 bits.
With 4 free words, expected ~2^(128-32) = 2^96 solutions.
Random search with 2^32 trials should find one.

Actually simplest: we have freedom in W[12..15] (128 bits).
Fix W[0..11]. For each random (W12, W13, W14, W15), compute De17.
Need De17 = 0 → 32-bit condition → birthday with 2^16 trials.

Wait — it's NOT birthday. It's direct search: P(De17=0) = 1/2^32.
With 128 bits of freedom: 2^32 trials expected.

BUT: if we use two groups (birthday on 2 halves):
  Group A: vary (W12, W13) → compute De17 → store
  Group B: vary (W14, W15) → compute De17 → look for match
  If De17(A) = De17(B) → the XOR differential of A and B gives De17=0?
  NO — that's not how it works because De17 is nonlinear.

The CORRECT birthday approach:
  De17 depends on W[0..15] and DW[0..15].
  Fix everything except DW[14].
  De17(DW14) maps 32 bits → 32 bits.
  Find DW14 such that De17 = 0: ~2^32 trials.

  With TWO free parameters (DW13, DW14):
  De17(DW13, DW14) maps 64 bits → 32 bits.
  ~2^32 images per 2^32 target values.
  Random search: P(hit 0) ≈ 1/2^32 per trial → 2^32 trials.

  Birthday DOES apply if we can split:
  Fix DW14, vary DW13: get list L1 of De17 values
  Fix DW13, vary DW14: get list L2 of De17 values
  Find De17_1 in L1 that equals De17_2 in L2?
  NO — we need De17 = 0, not a collision.

The REAL trick: for each random (W12, W13, W14, W15),
De17 is approximately random 32-bit → need 2^32 trials.
500K wasn't enough. Let's try smarter.

IDEA: Multi-target. For MANY different base messages W[0..11],
compute the required Da13 values. Then for each free-word choice,
check if ANY of the base messages is satisfied.

If we have N base messages and M free-word trials:
P(success) ≈ 1 - (1 - 1/2^32)^(N*M) ≈ N*M/2^32

With N=M=2^16: P ≈ 2^32/2^32 = 1  ← EXPECTED SUCCESS!
"""

import random
import numpy as np
from collections import defaultdict

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

def sha_rounds_fast(state, W64, r_start, r_end):
    a,b,c,d,e,f,g,h = state
    for r in range(r_start, r_end):
        T1 = add(add(add(add(h, Sig1(e)), Ch(e,f,g)), K[r]), W64[r])
        T2 = add(Sig0(a), Maj(a,b,c))
        a,b,c,d,e,f,g,h = add(T1,T2),a,b,c,add(d,T1),e,f,g
    return [a,b,c,d,e,f,g,h]

DW_bit = 0x80000000

# ============================================================
# PART A: Precompute states at round 12 (before free words)
# ============================================================
print("=" * 70)
print("PART A: Precompute states through round 12")
print("=" * 70)
print()

# Since Da13 depends on W[0..11], and W[12..15] are free,
# we can precompute the state after round 11 (uses W[0..11])
# Then rounds 12..16 use W[12..16] where W[12..15] are free

# Actually: round 12 uses W[12], so state after round 12 depends on W[12]
# But Da13 = state[0] after round 13, which depends on W[12] via round 12
# Wait — Step 5 showed W[12] is FREE (Da13 doesn't depend on it)
# This means round 12 processes W[12] but Da doesn't change?!
# That seems wrong. Let me recheck.

# Actually, Da13 is defined relative to two messages: W and W' = W ⊕ DW
# where DW = {0: 0x80000000}. Round 12 processes W[12] identically
# for both messages (same W[12]), so the DIFFERENTIAL doesn't change.
# But the STATE does depend on W[12].

print("Approach: fix W[0..11], precompute state at round 12.")
print("Then vary W[12..15] to search for De17=0.")
print()

# Precompute for one base message
W_base = [random.randint(0, MASK) for _ in range(16)]
W64_base = schedule(W_base)
W_base_f = list(W_base); W_base_f[0] ^= DW_bit
W64_base_f = schedule(W_base_f)

# State after round 11 (uses W[0..11])
state_n_12 = sha_rounds_fast(list(IV), W64_base, 0, 12)
state_f_12 = sha_rounds_fast(list(IV), W64_base_f, 0, 12)

print(f"State at r=12 (normal):    a={hex(state_n_12[0])}")
print(f"State at r=12 (perturbed): a={hex(state_f_12[0])}")
print(f"Da at r=12: {hex(sub(state_f_12[0], state_n_12[0]))}")
print()

# Now: for each choice of W[12..15], we need to:
# 1. Recompute schedule W[16] = f(W[14], W[9], W[1], W[0])
#    W[16] depends on W[14] — which is free!
# 2. Run rounds 12..16 with the new W[12..16]
# 3. Check De[17] = e_n[17] XOR e_f[17]

def compute_de17_fast(state_n_12, state_f_12, W_full, W_full_f, W12, W13, W14, W15):
    """Compute De17 given precomputed states at r=12 and free words."""
    # Update schedule: only W[12..15] change, affecting W[16+]
    # W[16] = sig1(W[14]) + W[9] + sig0(W[1]) + W[0]
    # W[17] = sig1(W[15]) + W[10] + sig0(W[2]) + W[1]
    # We need W[12..16] for rounds 12..16

    W_n = list(W_full)
    W_f = list(W_full_f)
    W_n[12] = W12; W_n[13] = W13; W_n[14] = W14; W_n[15] = W15
    W_f[12] = W12; W_f[13] = W13; W_f[14] = W14; W_f[15] = W15
    # Recompute W[16]
    W_n[16] = add(add(add(sig1(W_n[14]), W_n[9]), sig0(W_n[1])), W_n[0])
    W_f[16] = add(add(add(sig1(W_f[14]), W_f[9]), sig0(W_f[1])), W_f[0])

    # Run rounds 12..16
    s_n = sha_rounds_fast(list(state_n_12), W_n, 12, 17)
    s_f = sha_rounds_fast(list(state_f_12), W_f, 12, 17)

    return s_n[4] ^ s_f[4]  # De17

# Verify the fast computation matches slow:
de17_slow_n = sha_rounds_fast(list(IV), W64_base, 0, 17)
de17_slow_f = sha_rounds_fast(list(IV), W64_base_f, 0, 17)
de17_slow = de17_slow_n[4] ^ de17_slow_f[4]

de17_fast = compute_de17_fast(state_n_12, state_f_12, W64_base, W64_base_f,
                               W_base[12], W_base[13], W_base[14], W_base[15])
print(f"Slow De17: {hex(de17_slow)}")
print(f"Fast De17: {hex(de17_fast)}")
print(f"Match: {de17_slow == de17_fast} ✓" if de17_slow == de17_fast else "MISMATCH!")
print()

# ============================================================
# PART B: MITM-style search using W[12] and W[14]
# ============================================================
print("=" * 70)
print("PART B: Birthday/multi-target search")
print("=" * 70)
print()

# Strategy: fix W[13]=W[15]=0. Vary W[12] and W[14].
# Phase 1: for N1 random W[14] values, compute De17 → store in dict
# Phase 2: for N2 random W[12] values, compute De17 → check if in dict
# If a match exists → the W[14] from phase 1 and W[12] from phase 2
# give the same De17 value → but we need De17=0, not a collision!

# CORRECT approach: just random search with 2 free words.
# With 64 bits of freedom, 32-bit target: need 2^32 / 2^32 = 1 on average.
# But in PRACTICE each trial costs O(1), so we want P(hit) ≈ 1/2^32 per trial.

# Multi-target: fix M base messages, for each find De17.
# Then for each free-word trial, check against ALL M targets.
# P(hit any) ≈ M/2^32 per trial.
# With M=2^16 and 2^16 trials: P ≈ 2^32/2^32 = 1

print("Multi-target approach:")
print("  Phase 1: Generate M base messages (vary W[0..11])")
print("           Each gives a unique De17 for default free words")
print("  Phase 2: Vary free words, check if De17 hits ANY target")
print()

# Phase 1: generate base messages and their "default" states
M = 2000  # number of base messages
base_messages = []
for _ in range(M):
    W0_11 = [random.randint(0, MASK) for _ in range(12)]
    W16 = W0_11 + [0, 0, 0, 0]  # default free words = 0
    W64 = schedule(W16)
    W16_f = list(W16); W16_f[0] ^= DW_bit
    W64_f = schedule(W16_f)

    s_n_12 = sha_rounds_fast(list(IV), W64, 0, 12)
    s_f_12 = sha_rounds_fast(list(IV), W64_f, 0, 12)

    base_messages.append({
        'W0_11': W0_11,
        'W64': W64,
        'W64_f': W64_f,
        'state_n_12': s_n_12,
        'state_f_12': s_f_12,
    })

print(f"  Generated {M} base messages")

# Phase 2: for each free-word choice, compute De17 for ALL bases
# and check if any = 0
print(f"  Searching with random free words...")
print()

N_trials = 50000  # number of free-word trials
found = False
best_hw = 32
best_info = None

for trial in range(N_trials):
    W12 = random.randint(0, MASK)
    W13 = random.randint(0, MASK)
    W14 = random.randint(0, MASK)
    W15 = random.randint(0, MASK)

    for idx, bm in enumerate(base_messages):
        de17 = compute_de17_fast(bm['state_n_12'], bm['state_f_12'],
                                  bm['W64'], bm['W64_f'],
                                  W12, W13, W14, W15)
        h = hw(de17)
        if h < best_hw:
            best_hw = h
            best_info = {
                'msg_idx': idx,
                'W12': W12, 'W13': W13, 'W14': W14, 'W15': W15,
                'De17': de17,
            }
        if h == 0:
            found = True
            break

    if found:
        break

    if (trial + 1) % 10000 == 0:
        print(f"    Trial {trial+1}/{N_trials}: best HW = {best_hw}")

print()
total_checks = (trial + 1) * M
print(f"Total De17 evaluations: {total_checks:,}")
print(f"Best De17 HW: {best_hw}")

if found:
    print(f"\n★★★ De17 = 0 FOUND! ★★★")
    print(f"  Message index: {best_info['msg_idx']}")
    print(f"  W[12] = {hex(best_info['W12'])}")
    print(f"  W[13] = {hex(best_info['W13'])}")
    print(f"  W[14] = {hex(best_info['W14'])}")
    print(f"  W[15] = {hex(best_info['W15'])}")
    print(f"  De17 = {hex(best_info['De17'])}")

    # VERIFY
    bm = base_messages[best_info['msg_idx']]
    W16_sol = bm['W0_11'] + [best_info['W12'], best_info['W13'],
                               best_info['W14'], best_info['W15']]
    W64_sol = schedule(W16_sol)
    W16_sol_f = list(W16_sol); W16_sol_f[0] ^= DW_bit
    W64_sol_f = schedule(W16_sol_f)

    s_n_full = sha_rounds_fast(list(IV), W64_sol, 0, 17)
    s_f_full = sha_rounds_fast(list(IV), W64_sol_f, 0, 17)

    de17_verify = s_n_full[4] ^ s_f_full[4]
    print(f"\n  FULL VERIFICATION: De17 = {hex(de17_verify)}")
    if de17_verify == 0:
        print(f"  ★★★ VERIFIED — De17 = 0 CONFIRMED! ★★★")
        print(f"\n  Full message W[0..15]:")
        for i, w in enumerate(W16_sol):
            print(f"    W[{i:2d}] = {hex(w)}")
        print(f"\n  Perturbation: W'[0] = W[0] ^ 0x80000000 = {hex(W16_sol[0] ^ DW_bit)}")

        # Check all De values
        print(f"\n  De profile through 17 rounds:")
        s_n = list(IV)
        s_f = list(IV)
        for r in range(17):
            s_n = sha_round(s_n, W64_sol[r], K[r])
            s_f = sha_round(s_f, W64_sol_f[r], K[r])
            de = s_n[4] ^ s_f[4]
            print(f"    r={r+1:2d}: De = {hex(de):>12s}  HW={hw(de):2d}")
    else:
        print(f"  VERIFICATION FAILED — De17 = {hex(de17_verify)}")
elif best_info:
    print(f"\n  Closest solution:")
    print(f"  Message index: {best_info['msg_idx']}")
    print(f"  W[14] = {hex(best_info['W14'])}")
    print(f"  De17 = {hex(best_info['De17'])}  HW = {best_hw}")

print()
print("=" * 70)
print("SUMMARY STEP 7")
print("=" * 70)
print()
print(f"Total evaluations: {total_checks:,}")
print(f"Best HW(De17): {best_hw}")
expected_evals = 2**32 / M
print(f"Expected evaluations for hit: ~{expected_evals:,.0f} "
      f"(with {M} base messages)")
if found:
    print(f"\n★ SUCCESS: Found a message pair with De17=0!")
    print(f"  This breaks the Wang chain barrier at round 17")
    print(f"  using multi-word differential + free word exploitation")
else:
    print(f"\nNot found yet — need more trials or larger M")
    print(f"Estimated total needed: ~{expected_evals:,.0f}")
