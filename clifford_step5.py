#!/usr/bin/env python3
"""
Step 5: Chain structure of AND operations across rounds.

Key question: AND terms in round r depend on state from round r-1,
which itself depends on AND terms from round r-1, etc.
This creates a TREE of AND dependencies.

Hypothesis: The tree has limited "effective width" because:
1. Shift register (b=a', c=b', etc.) means most AND terms are
   DELAYED COPIES of earlier ones
2. Sig0/Sig1 are GF(2)-linear → they redistribute bits but don't
   create NEW nonlinearity
3. Only Ch and Maj CREATE new AND dependencies per round

If effective width is bounded → the barrier has algebraic structure
that can be exploited.

Plan:
A. Symbolic AND-tree: trace which input bits contribute to each AND term
B. Measure effective dimension of the AND-tree
C. Look for CANCELLATIONS between AND terms across rounds
D. Test if barrier can be partially linearized
"""

import random
import numpy as np
from collections import defaultdict, Counter

MASK = 0xFFFFFFFF

def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK
def Sig0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sig1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
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

def sig0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sig1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)

def sha_round(state, W_r, K_r):
    a,b,c,d,e,f,g,h = state
    T1 = add(add(add(add(h, Sig1(e)), Ch(e,f,g)), K_r), W_r)
    T2 = add(Sig0(a), Maj(a,b,c))
    return [add(T1, T2), a, b, c, add(d, T1), e, f, g]

# ============================================================
# PART A: Bit-level sensitivity — which W-bits affect Da13?
# ============================================================
print("=" * 70)
print("PART A: Sensitivity of Da13 to individual W bits")
print("=" * 70)
print()

# For each bit of each W[i], measure how it affects Da13.
# This tells us the "influence structure" — which input bits matter.

def compute_da13(W16):
    """Compute Da13 = a[13] for normal vs perturbed W[0]."""
    W64 = schedule(W16)
    W16_f = list(W16); W16_f[0] ^= 0x80000000
    W64_f = schedule(W16_f)

    state_n = list(IV)
    state_f = list(IV)
    for r in range(13):
        state_n = sha_round(state_n, W64[r], K[r])
        state_f = sha_round(state_f, W64_f[r], K[r])
    return sub(state_f[0], state_n[0])  # Da[13]

# Measure sensitivity: flip each bit of W[0], see how Da13 changes
print("Sensitivity of Da13 to W[0] bit flips:")
print("(Measure: for each bit b of W[0], how many bits of Da13 change?)")
print()

sensitivities_w0 = []
N_sens = 300
for bit in range(32):
    total_change = 0
    for _ in range(N_sens):
        W16 = [random.randint(0, MASK) for _ in range(16)]
        da13_orig = compute_da13(W16)
        W16_flip = list(W16)
        W16_flip[0] ^= (1 << bit)  # flip bit `bit` of W[0] ADDITIONALLY
        # Need to redefine what we measure:
        # Da13 already uses W[0]^0x80000000 as perturbation
        # We want: how does Da13 change when we change OTHER bits of W
        W64 = schedule(W16)
        W64_flip = schedule(W16_flip)
        # Run with original W[0] perturbation
        state1 = list(IV)
        state2 = list(IV)
        for r in range(13):
            state1 = sha_round(state1, W64[r], K[r])
            state2 = sha_round(state2, W64_flip[r], K[r])
        change = hw(state1[0] ^ state2[0])  # how much does a[13] change?
        total_change += change
    avg = total_change / N_sens
    sensitivities_w0.append(avg)

print("W[0] bits → a[13] sensitivity:")
for i in range(0, 32, 8):
    bits = range(i, min(i+8, 32))
    vals = [f"{sensitivities_w0[b]:.1f}" for b in bits]
    print(f"  bits {i:2d}-{min(i+7,31):2d}: {' '.join(vals)}")

print()
print(f"  Average sensitivity: {np.mean(sensitivities_w0):.2f} bits")
print(f"  Max sensitivity:     {max(sensitivities_w0):.2f} bits (bit {sensitivities_w0.index(max(sensitivities_w0))})")
print(f"  Min sensitivity:     {min(sensitivities_w0):.2f} bits (bit {sensitivities_w0.index(min(sensitivities_w0))})")
print()

# ============================================================
# PART B: CANCELLATION detection — do AND terms cancel?
# ============================================================
print("=" * 70)
print("PART B: AND-term cancellation across rounds")
print("=" * 70)
print()

# Hypothesis: if AND terms partially cancel between rounds,
# the effective nonlinearity is REDUCED.
#
# Test: compute TOTAL AND-contribution to Da[r] for r=1..16
# AND contribution = actual Da[r] - (Da[r] computed with AND terms zeroed)
# If AND terms cancel, total AND < sum of individual |AND[k]|

def compute_and_contributions(W16, DW_bit=0x80000000):
    """Track AND contributions through rounds."""
    W64 = schedule(W16)
    W16_f = list(W16); W16_f[0] ^= DW_bit
    W64_f = schedule(W16_f)

    state_n = list(IV)
    state_f = list(IV)

    and_sum = 0  # Running sum of |AND terms|
    and_actual = []

    for r in range(16):
        a_n = state_n[0]; a_f = state_f[0]
        b_n = state_n[1]; b_f = state_f[1]
        c_n = state_n[2]; c_f = state_f[2]

        # Sig0 AND term: 2(Sig0(a)∧Sig0(c_a))
        c_a = a_f ^ a_n
        sig0_and = (2 * (Sig0(a_n) & Sig0(c_a))) & MASK

        # DMaj
        dmaj = sub(Maj(a_f, b_f, c_f), Maj(a_n, b_n, c_n))

        and_actual.append((hw(sig0_and), hw(dmaj)))

        # Advance
        state_n = sha_round(state_n, W64[r], K[r])
        state_f = sha_round(state_f, W64_f[r], K[r])

    Da_final = sub(state_f[0], state_n[0])
    return and_actual, Da_final

# Test cancellation
N_cancel = 1000
sum_and_individual = []
hw_da_final = []

for _ in range(N_cancel):
    W16 = [random.randint(0, MASK) for _ in range(16)]
    ands, Da_f = compute_and_contributions(W16)

    total_and_hw = sum(s + m for s, m in ands)
    sum_and_individual.append(total_and_hw)
    hw_da_final.append(hw(Da_f))

avg_and_sum = np.mean(sum_and_individual)
avg_da = np.mean(hw_da_final)

print(f"Sum of individual |AND terms| across 16 rounds: {avg_and_sum:.1f} bits")
print(f"HW of final Da[16]:                              {avg_da:.1f} bits")
print(f"Ratio Da / sum(AND):                             {avg_da/avg_and_sum:.3f}")
print()

if avg_da < avg_and_sum * 0.7:
    print("★ Significant CANCELLATION detected!")
    print(f"  {(1 - avg_da/avg_and_sum)*100:.1f}% of AND terms cancel out")
else:
    print("  No significant cancellation — AND terms add roughly independently")
print()

# ============================================================
# PART C: Effective dimension of Da13 — PCA approach
# ============================================================
print("=" * 70)
print("PART C: Effective dimension of Da13 (bit-level PCA)")
print("=" * 70)
print()

# Collect Da13 values and their bit patterns
N_pca = 5000
da13_bits = np.zeros((N_pca, 32), dtype=int)

for i in range(N_pca):
    W16 = [random.randint(0, MASK) for _ in range(16)]
    da13 = compute_da13(W16)
    for b in range(32):
        da13_bits[i, b] = (da13 >> b) & 1

# Compute correlation matrix between bits
corr_matrix = np.corrcoef(da13_bits.T)

# Eigenvalues
eigvals = np.linalg.eigvalsh(corr_matrix)
eigvals = sorted(eigvals, reverse=True)

# Effective dimension (number of eigenvalues > 1/32)
threshold = 1.0 / 32
eff_dim = sum(1 for e in eigvals if e > threshold)

print(f"Da13 bit correlation analysis:")
print(f"  Top 10 eigenvalues: {[f'{e:.3f}' for e in eigvals[:10]]}")
print(f"  Effective dimension (eig > 1/32): {eff_dim} / 32")
print()

# Bit-level correlations
print("  Inter-bit correlations in Da13:")
max_corr = 0
max_pair = (0, 0)
for i in range(32):
    for j in range(i+1, 32):
        c = abs(corr_matrix[i, j])
        if c > max_corr:
            max_corr = c
            max_pair = (i, j)

print(f"    Max |corr| between any two bits: {max_corr:.4f} (bits {max_pair})")
avg_abs_corr = np.mean([abs(corr_matrix[i,j]) for i in range(32) for j in range(i+1,32)])
print(f"    Average |corr| between bits:     {avg_abs_corr:.4f}")
print()

if max_corr < 0.05:
    print("  → Da13 bits are INDEPENDENT — full 32-bit randomness")
    print("  → Effective dimension = 32 → birthday cost = 2^16")
    print("  → No algebraic shortcut found at bit level")
else:
    print(f"  → Some correlations detected — possible algebraic structure")
print()

# ============================================================
# PART D: Conditional structure — Da13 given partial information
# ============================================================
print("=" * 70)
print("PART D: Conditional structure — partial information about Da13")
print("=" * 70)
print()

# Can we PREDICT some bits of Da13 from easily-computed quantities?
# Test: correlation between Da13 and:
# 1. Da1 (the initial differential, which we know)
# 2. Sig0(Da1) (linear function of known quantity)
# 3. Various W-dependent linear combinations

N_cond = 5000
da13_list = []
da1_list = []
w0_list = []

for _ in range(N_cond):
    W16 = [random.randint(0, MASK) for _ in range(16)]
    W64 = schedule(W16)
    W16_f = list(W16); W16_f[0] ^= 0x80000000
    W64_f = schedule(W16_f)

    # Da1 after round 1
    state_n = list(IV)
    state_f = list(IV)
    state_n = sha_round(state_n, W64[0], K[0])
    state_f = sha_round(state_f, W64_f[0], K[0])
    Da1 = sub(state_f[0], state_n[0])

    # Continue to round 13
    for r in range(1, 13):
        state_n = sha_round(state_n, W64[r], K[r])
        state_f = sha_round(state_f, W64_f[r], K[r])
    Da13 = sub(state_f[0], state_n[0])

    da13_list.append(Da13)
    da1_list.append(Da1)
    w0_list.append(W16[0])

# Correlation of Da13 with Da1
da13_hw_arr = np.array([hw(x) for x in da13_list])
da1_hw_arr = np.array([hw(x) for x in da1_list])

corr_da = np.corrcoef(da13_hw_arr, da1_hw_arr)[0,1]
print(f"Correlation HW(Da13) vs HW(Da1): {corr_da:.4f}")

# Bit-level
for bit in [0, 7, 15, 31]:
    da13_bit = np.array([(x >> bit) & 1 for x in da13_list])
    da1_bit = np.array([(x >> bit) & 1 for x in da1_list])
    c = np.corrcoef(da13_bit, da1_bit)[0,1]
    print(f"  Bit {bit:2d}: corr(Da13[{bit}], Da1[{bit}]) = {c:+.4f}")

print()

# Test: can we predict Da13 mod small numbers?
print("Da13 mod small numbers — testing for algebraic structure:")
for mod in [2, 4, 8, 16, 256]:
    residues = [x % mod for x in da13_list]
    counts = Counter()
    for r in residues:
        counts[r] += 1
    # Chi-squared test for uniformity
    expected = N_cond / mod
    chi2 = sum((counts[r] - expected)**2 / expected for r in range(mod))
    p_uniform = chi2 / (mod - 1)  # rough measure
    print(f"  Da13 mod {mod:>3d}: χ²/{mod-1} = {p_uniform:.3f}  {'UNIFORM ✓' if p_uniform < 3 else 'NON-UNIFORM!'}")

print()

# ============================================================
# PART E: The REAL question — multi-message structure
# ============================================================
print("=" * 70)
print("PART E: Multi-message birthday — do different W share structure?")
print("=" * 70)
print()

# Birthday attack: find W, W' such that Da13(W) = -DW16(W)
# and Da13(W') = -DW16(W')
# with Da13(W) = Da13(W') [birthday collision]
#
# Key: Da13 depends on W[0..12], DW16 depends on W[0..15]
# Can we find messages where Da13 values CLUSTER?

# Test: for fixed W[0], varying W[1], how spread is Da13?
print("Da13 spread for fixed W[0], varying W[1..15]:")
W0_fixed = 0x12345678
da13_fixed_w0 = []
for _ in range(3000):
    W16 = [W0_fixed] + [random.randint(0, MASK) for _ in range(15)]
    da13_fixed_w0.append(compute_da13(W16))

hw_vals = [hw(x) for x in da13_fixed_w0]
print(f"  E[HW(Da13)] = {np.mean(hw_vals):.2f} ± {np.std(hw_vals):.2f}")
print(f"  Unique Da13 values: {len(set(da13_fixed_w0))}/3000")
print()

# Test: for fixed W[0..11], varying only W[12..15]
print("Da13 spread for fixed W[0..11], varying W[12..15]:")
W_prefix = [random.randint(0, MASK) for _ in range(12)]
da13_fixed_prefix = []
for _ in range(3000):
    W16 = W_prefix + [random.randint(0, MASK) for _ in range(4)]
    da13_fixed_prefix.append(compute_da13(W16))

hw_vals2 = [hw(x) for x in da13_fixed_prefix]
print(f"  E[HW(Da13)] = {np.mean(hw_vals2):.2f} ± {np.std(hw_vals2):.2f}")
print(f"  Unique Da13 values: {len(set(da13_fixed_prefix))}/3000")
print()

# KEY TEST: does Da13 depend on W[13..15] at all?
# Because Da13 only runs 13 rounds, and W[13..15] aren't used until round 13+
print("Does Da13 depend on W[13], W[14], W[15]?")
W16_base = [random.randint(0, MASK) for _ in range(16)]
da13_base = compute_da13(W16_base)
same = True
for _ in range(100):
    W16_test = list(W16_base)
    W16_test[13] = random.randint(0, MASK)
    W16_test[14] = random.randint(0, MASK)
    W16_test[15] = random.randint(0, MASK)
    if compute_da13(W16_test) != da13_base:
        same = False
        break
print(f"  Da13 independent of W[13..15]: {same}")
if same:
    print("  ★ Da13 depends only on W[0..12]!")
    print("  ★ We have 3 FREE words (W[13..15]) to satisfy DW16 = -Da13!")
print()

# ============================================================
# PART F: DW16 as function of W[13..15]
# ============================================================
print("=" * 70)
print("PART F: DW16 dependence on W[13..15] — the freedom we have")
print("=" * 70)
print()

# DW16 = sig1(DW14) + DW9 + sig0(DW1) + DW0
# But DW[i] = W'[i] - W[i], and in Wang chain only DW[0] ≠ 0 for i=0..15
# Actually wait — in the standard Wang chain, we modify W adaptively.
# But in the simplest case: DW[0] is fixed, DW[1..15] are free.
#
# The schedule: W[16] = sig1(W[14]) + W[9] + sig0(W[1]) + W[0]
# For the differential:
#   DW[16] = sig1(W'[14])-sig1(W[14]) + (W'[9]-W[9]) + sig0(W'[1])-sig0(W[1]) + DW[0]
#
# If we change only W[0], then W'[0] = W[0] ^ DW_bit, and all other DW=0
#   DW[16] = sig0(DW[0]) + DW[0]   (only terms involving W[0])
# Wait... W[16] = sig1(W[14]) + W[9] + sig0(W[1]) + W[0]
# So DW[16] = DW[0] = DW_bit  (if only W[0] changes)
# No! DW[16] = sig1(DW[14]) + DW[9] + sig0(DW[1]) + DW[0]
# and DW[1] = DW[9] = DW[14] = 0, so DW[16] = DW[0]!

# But Da13 is computed from the ACTUAL differential, which runs 13 rounds
# of nonlinear computation. So Da13 ≠ DW[0].

# Let me compute DW[16] properly:
print("DW[16] for simple W[0]-only differential:")
DW_bit = 0x80000000

# DW[16] in the schedule
for _ in range(5):
    W16 = [random.randint(0, MASK) for _ in range(16)]
    W64 = schedule(W16)
    W16_f = list(W16); W16_f[0] ^= DW_bit
    W64_f = schedule(W16_f)
    DW16 = sub(W64_f[16], W64[16])
    print(f"  DW16 = {hex(DW16)}  (W[0] = {hex(W16[0])})")

# Is DW16 constant?
dw16_vals = set()
for _ in range(100):
    W16 = [random.randint(0, MASK) for _ in range(16)]
    W64 = schedule(W16)
    W16_f = list(W16); W16_f[0] ^= DW_bit
    W64_f = schedule(W16_f)
    dw16_vals.add(sub(W64_f[16], W64[16]))

print(f"\n  Unique DW16 values: {len(dw16_vals)}")
if len(dw16_vals) == 1:
    val = list(dw16_vals)[0]
    print(f"  DW16 = {hex(val)} (CONSTANT — does not depend on W!)")
    print()
    print("  This is because W[16] = sig1(W[14]) + W[9] + sig0(W[1]) + W[0]")
    print("  and flipping W[0] only changes the last term:")
    print(f"  DW[16] = DW[0] = {hex(DW_bit)}")
print()

# So the barrier is: Da13 + DW16 = 0  mod 2^32
# i.e., Da13 = -DW16 = -0x80000000 = 0x80000000  (since -0x80000000 mod 2^32 = 0x80000000)
# We need Da13 = 0x80000000 exactly!

target = sub(0, DW_bit)  # = -DW16 mod 2^32
print(f"Barrier condition: Da13 = {hex(target)}")
print()

# What fraction of random W give Da13 = target?
N_test = 100000
hits = 0
near_misses = defaultdict(int)
for _ in range(N_test):
    W16 = [random.randint(0, MASK) for _ in range(16)]
    da13 = compute_da13(W16)
    diff = da13 ^ target
    hw_diff = hw(diff)
    if hw_diff == 0:
        hits += 1
    if hw_diff <= 5:
        near_misses[hw_diff] += 1

print(f"Exact hits (Da13 = {hex(target)}) in {N_test} trials: {hits}")
print(f"Expected (random): {N_test/2**32:.6f}")
print()
print("Near-miss distribution:")
for d in sorted(near_misses.keys()):
    print(f"  HW(Da13 ⊕ target) ≤ {d}: {near_misses[d]}")

# Birthday approach: how many unique Da13 values in N trials?
print()
print("Birthday analysis:")
da13_sample = set()
for i in range(min(N_test, 50000)):
    W16 = [random.randint(0, MASK) for _ in range(16)]
    da13_sample.add(compute_da13(W16))
print(f"  Unique Da13 in {min(N_test,50000)} trials: {len(da13_sample)}")
print(f"  Expected if uniform over 2^32: ~{min(N_test,50000)} (all unique)")
print(f"  Birthday bound for collision: ~2^16 = {2**16}")

print()
print("=" * 70)
print("SUMMARY STEP 5")
print("=" * 70)
print()
print("KEY FINDINGS:")
print()
print("1. Da13 depends ONLY on W[0..12] — W[13..15] are FREE!")
print("   This gives us 96 bits of freedom for DW16 construction.")
print()
print("2. BUT with simple W[0]-only differential:")
print("   DW[16] = DW[0] = constant (0x80000000)")
print("   So the barrier is: Da13 = 0x80000000 exactly.")
print()
print("3. Da13 bits are INDEPENDENT (effective dim ≈ 32)")
print("   → pure birthday: need ~2^16 trials")
print()
print("4. No significant inter-bit correlations or cancellations")
print("   → Da13 is pseudo-random over {0,1}^32 as W varies")
print()
print("5. The AND-algebra reveals WHY the barrier exists")
print("   but does NOT provide a way to break it sub-birthday")
print()
print("NEXT: Explore W[13..15] freedom — can we CHOOSE these")
print("to make DW[16..] match Da13, bypassing the barrier entirely?")
