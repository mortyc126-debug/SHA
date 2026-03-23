#!/usr/bin/env python3
"""
Step 3: Carry as structural operator — the W_32 algebra.

Key insight from Step 2: GF(2) model fails (gap ~30 bits).
But the gap SATURATES — carry is bounded, not chaotic.

New approach: work in W_32 = (Z/2^32Z, +) where:
  - Addition is LINEAR (no carry noise!)
  - Ch, Maj are the ONLY nonlinearities
  - Sig0, Sig1, sig0, sig1 are NONLINEAR (rotations ≠ Z-linear)

Question 1: In W_32, how does the differential look?
  δa' = δT1 + δT2  (exact in Z, no carry!)
  where δT1 = δh + δSig1 + δCh + δW  (but Sig1 is NOT Z-linear!)

Question 2: What if we work in ADDITIVE differentials (Z/2^32Z)?
  Da' = Da(T1) + Da(T2)
  Da(T1) = Dh + Da(Sig1(e)) + Da(Ch(e,f,g)) + DW

  Addition is linear: Da(x+y) = Da(x) + Da(y) exactly!
  Sig1 is NOT additive: Da(Sig1(e)) ≠ Sig1(Da(e))

  BUT: Da(Sig1(e)) = Sig1(e⊕De) - Sig1(e) + CARRY_CORRECTION

  The key: can we characterize CARRY_CORRECTION?

Question 3: What is the ALGEBRA of additive differentials?
  Da operates on Z/2^32Z. Ch, Maj are quadratic over GF(2).
  But over Z/2^32Z, Ch and Maj have a different algebraic structure.
  Ch(e,f,g) = (e & f) ^ (~e & g) — this mixes XOR and AND.

  Over W_32: Ch(e+de, f, g) - Ch(e, f, g) = ???
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
def sub(x, y): return (x - y) & MASK

IV = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]
K = [0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
     0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
     0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
     0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174]

# ============================================================
# PART A: Additive differential of Ch
# ============================================================
print("=" * 70)
print("PART A: Additive differential of Ch")
print("=" * 70)
print()

# Ch(e,f,g) = (e & f) ^ (~e & g)
# DCh = Ch(e+de, f+df, g+dg) - Ch(e, f, g)  mod 2^32
#
# For de only (df=dg=0):
# DCh(e, de) = Ch(e+de, f, g) - Ch(e, f, g)
#
# Key: Ch is NOT additive over Z/2^32Z.
# But what is the STRUCTURE of DCh?

N = 10000
dch_hw = []
de_hw = []
# Fix f, g from IV
f_val, g_val = IV[5], IV[6]

for _ in range(N):
    e = random.randint(0, MASK)
    de = 1  # minimal additive perturbation

    ch1 = Ch(e, f_val, g_val)
    ch2 = Ch(add(e, de), f_val, g_val)
    dch = sub(ch2, ch1)

    dch_hw.append(hw(dch))
    de_hw.append(hw(de))

print(f"DCh(e, de=1, f=IV[5], g=IV[6]):")
print(f"  E[HW(DCh)]:  {np.mean(dch_hw):.2f} ± {np.std(dch_hw):.2f}")
print(f"  P(DCh=0):    {dch_hw.count(0)/N:.4f}")
print(f"  P(DCh=1):    {sum(1 for x in dch_hw if x==1)/N:.4f}")
print()

# More interesting: DCh as function of e for fixed de=1
# Is there structure?
dch_values = Counter()
for e in range(0, 2**16, 1):  # sample 2^16 values of e
    ch1 = Ch(e, f_val, g_val)
    ch2 = Ch(add(e, 1), f_val, g_val)
    dch = sub(ch2, ch1)
    dch_values[dch] += 1

print(f"Unique DCh values (e ∈ [0, 2^16), de=1): {len(dch_values)}")
print(f"Top DCh values:")
for val, cnt in dch_values.most_common(5):
    print(f"  DCh = {hex(val):>12s}  count={cnt:5d}  P={cnt/65536:.4f}  HW={hw(val)}")
print()

# ============================================================
# PART B: Additive differential of Sig1
# ============================================================
print("=" * 70)
print("PART B: Additive differential of Sig1 — the CRITICAL function")
print("=" * 70)
print()

# Sig1 is XOR of three rotations: rotr(6) ^ rotr(11) ^ rotr(25)
# Over GF(2): Sig1 is LINEAR → DSig1_xor = Sig1(de)  (constant!)
# Over Z/2^32Z: DSig1_add = Sig1(e+de) - Sig1(e) ≠ Sig1(de) in general
#
# The GAP: Sig1(e+de) - Sig1(e) - Sig1(de) = ???
# This gap IS the carry correction.

N = 10000
sig1_gap_hw = []
sig1_linear_hw = []

for _ in range(N):
    e = random.randint(0, MASK)
    de = 1

    actual = sub(Sig1(add(e, de)), Sig1(e))  # Sig1(e+1) - Sig1(e)
    linear_pred = Sig1(de)  # = Sig1(1)  (GF2-linear prediction)
    gap = actual ^ linear_pred  # XOR difference between actual and linear

    sig1_gap_hw.append(hw(gap))
    sig1_linear_hw.append(hw(actual))

print(f"DSig1(e, de=1):")
print(f"  E[HW(actual DSig1)]:     {np.mean(sig1_linear_hw):.2f} ± {np.std(sig1_linear_hw):.2f}")
print(f"  HW(Sig1(1)) [linear]:    {hw(Sig1(1))}")
print(f"  E[HW(gap = actual⊕lin)]: {np.mean(sig1_gap_hw):.2f} ± {np.std(sig1_gap_hw):.2f}")
print(f"  P(gap=0):                {sig1_gap_hw.count(0)/N:.4f}")
print()

# KEY: The gap comes from carry propagation through the addition e+de.
# carry(e, 1) = (e+1) XOR e XOR 1 = carry_xor(e, 0)
# The carry chain length is geometric with P=1/2.
# So: DSig1 = Sig1(carry_xor(e,0))  (the carry pattern through Sig1)

# Verify this interpretation:
print("Verification: DSig1_add = Sig1(e⊕de) ⊕ carry_correction")
print("  carry(e, de=1) = (e+1) ⊕ e ⊕ 1")
print()

N2 = 5000
exact = 0
for _ in range(N2):
    e = random.randint(0, MASK)
    # (e+1) = e XOR carry_xor(e,0)
    carry_pattern = (add(e, 1)) ^ e  # = 1 if no carry, longer if carry
    actual_dsig = sub(Sig1(add(e, 1)), Sig1(e))
    # Sig1 is GF2-linear: Sig1(a XOR b) = Sig1(a) XOR Sig1(b)
    # So Sig1(e+1) = Sig1(e XOR carry_pattern) = Sig1(e) XOR Sig1(carry_pattern)
    # → DSig1_add = Sig1(e+1) - Sig1(e)
    #            = [Sig1(e) XOR Sig1(carry_pattern)] - Sig1(e)
    #
    # The question is: a XOR b - a = b + carry_correction_2
    # Not the same as b in general!

    # Let's check: is DSig1 = Sig1(carry_pattern) + correction?
    xor_part = Sig1(carry_pattern)
    add_dsig = sub(Sig1(add(e, 1)), Sig1(e))
    # add_dsig = xor_part + second_order_carry???
    gap2 = sub(add_dsig, xor_part)  # should this be structured?
    if gap2 == 0:
        exact += 1

print(f"  P(DSig1_add = Sig1(carry_pattern)):  {exact/N2:.4f}")
print(f"  (If 1.0 → Sig1 preserves additive structure through carry)")
print()

# ============================================================
# PART C: THE FUNDAMENTAL IDENTITY
# ============================================================
print("=" * 70)
print("PART C: Searching for the fundamental identity")
print("=" * 70)
print()

# IDEA: Since Sig1 is GF(2)-linear:
#   Sig1(e+de) = Sig1(e XOR carry_xor(e, de))
#              = Sig1(e) XOR Sig1(carry_xor(e, de))
#
# Therefore:
#   DSig1_ADD = Sig1(e+de) - Sig1(e)
#            = [Sig1(e) XOR Sig1(carry_xor(e,de))] - Sig1(e)
#
# Let A = Sig1(e), B = Sig1(carry_xor(e,de))
# Then: DSig1_ADD = (A XOR B) - A
#
# For any X, Y: (X XOR Y) - X = Y - 2*(X AND Y)  (mod 2^32)
# Proof: X XOR Y = X + Y - 2*(X AND Y)
#        → (X XOR Y) - X = Y - 2*(X AND Y)
#
# So: DSig1_ADD = Sig1(carry_xor(e,de)) - 2*(Sig1(e) AND Sig1(carry_xor(e,de)))
#
# THIS IS THE FUNDAMENTAL IDENTITY!

print("FUNDAMENTAL IDENTITY:")
print()
print("  DSig1_add(e, de) = Sig1(carry(e,de)) - 2·(Sig1(e) ∧ Sig1(carry(e,de)))")
print()
print("where carry(e,de) = (e+de) ⊕ e ⊕ de")
print()

# Verify
N3 = 10000
verified = 0
for _ in range(N3):
    e = random.randint(0, MASK)
    de = random.randint(0, MASK)

    carry_pat = add(e, de) ^ e ^ de
    A = Sig1(e)
    B = Sig1(carry_pat)

    actual = sub(Sig1(add(e, de)), Sig1(e))
    predicted = sub(B, (2 * (A & B)) & MASK)

    if actual == predicted:
        verified += 1

print(f"Verification: {verified}/{N3}")
if verified == N3:
    print("★ IDENTITY CONFIRMED — 100% exact! ★")
else:
    print(f"  Failures: {N3-verified}")
print()

# ============================================================
# PART D: Generalize — the (A XOR B) - A = B - 2(A∧B) identity
# ============================================================
print("=" * 70)
print("PART D: The universal carry identity")
print("=" * 70)
print()

# For ANY GF(2)-linear function L:
#   L(x+y) = L(x XOR carry(x,y))
#           = L(x) XOR L(carry(x,y))
#
#   D_add L(x,y) = L(x+y) - L(x)
#                = L(carry(x,y)) - 2·(L(x) ∧ L(carry(x,y)))
#
# This works for L = Sig0, Sig1, sig0, sig1.
# And carry(x,y) has geometric distribution: E[HW] = 2.

print("Universal identity for GF(2)-linear functions:")
print()
print("  D_add[L](x, y) = L(carry(x,y)) - 2·(L(x) ∧ L(carry(x,y)))")
print()
print("Verified for:")

for name, L in [("Sig0", Sig0), ("Sig1", Sig1), ("sig0", sig0), ("sig1", sig1)]:
    ok = 0
    for _ in range(5000):
        x = random.randint(0, MASK)
        y = random.randint(0, MASK)
        c = add(x, y) ^ x ^ y
        actual = sub(L(add(x, y)), L(x))
        pred = sub(L(c), (2 * (L(x) & L(c))) & MASK)
        if actual == pred:
            ok += 1
    print(f"  {name}: {ok}/5000 {'✓' if ok==5000 else '✗'}")

print()

# ============================================================
# PART E: One round in the NEW algebra
# ============================================================
print("=" * 70)
print("PART E: SHA-256 round in the carry-aware algebra")
print("=" * 70)
print()

# In the new algebra, one round differential (additive) decomposes as:
#
# DT1 = Dh + DSig1(e,De) + DCh(e,f,g,De,Df,Dg) + DW
#
# Using the fundamental identity:
#   DSig1(e,De) = Sig1(c_e) - 2·(Sig1(e) ∧ Sig1(c_e))
#   where c_e = carry(e, De) = (e+De) XOR e XOR De
#
# The carry pattern c_e depends on (e, De) — this is the state-dependence.
# But c_e has KNOWN distribution: geometric, E[HW]=2.
#
# Key question: can we separate the state-independent and state-dependent parts?

# State-independent part of DSig1:
#   When carry is minimal (c_e = De, i.e., no actual carry):
#     DSig1 = Sig1(De) - 2·(Sig1(e) ∧ Sig1(De))
#   = Sig1(De) [GF2-linear] - 2·[correction term depending on e]
#
# The correction 2·(Sig1(e) ∧ Sig1(c_e)) is the ONLY state-dependent part!

# How big is this correction?
N4 = 10000
correction_hw = []
for _ in range(N4):
    e = random.randint(0, MASK)
    De = 1  # minimal perturbation
    c_e = add(e, De) ^ e ^ De
    correction = (2 * (Sig1(e) & Sig1(c_e))) & MASK
    correction_hw.append(hw(correction))

print(f"Correction term 2·(Sig1(e) ∧ Sig1(carry(e,1))):")
print(f"  E[HW]: {np.mean(correction_hw):.2f} ± {np.std(correction_hw):.2f}")
print(f"  P(correction=0): {correction_hw.count(0)/N4:.4f}")
print(f"  This is SMALL compared to the 32-bit word!")
print()

# For De = random:
correction_hw2 = []
for _ in range(N4):
    e = random.randint(0, MASK)
    De = random.randint(0, MASK)
    c_e = add(e, De) ^ e ^ De
    correction = (2 * (Sig1(e) & Sig1(c_e))) & MASK
    correction_hw2.append(hw(correction))

print(f"Correction for random De:")
print(f"  E[HW]: {np.mean(correction_hw2):.2f} ± {np.std(correction_hw2):.2f}")
print()

# ============================================================
# PART F: The algebra structure
# ============================================================
print("=" * 70)
print("PART F: Structure of the carry-aware round algebra")
print("=" * 70)
print()
print("SHA-256 round in carry-aware decomposition:")
print()
print("  DT1 = Dh                                    [linear, W_32]")
print("      + Sig1(carry(e,De))                      [GF2-linear of carry]")
print("      - 2·(Sig1(e) ∧ Sig1(carry(e,De)))        [correction: quadratic in e]")
print("      + DCh_add(e,f,g,De,Df,Dg)                [quadratic in efg]")
print("      + DW                                      [linear, W_32]")
print()
print("  DT2 = Sig0(carry(a,Da))                      [GF2-linear of carry]")
print("      - 2·(Sig0(a) ∧ Sig0(carry(a,Da)))        [correction: quadratic in a]")
print("      + DMaj_add(a,b,c,Da,Db,Dc)               [quadratic in abc]")
print()
print("Structure: ROUND = LINEAR + QUADRATIC(state) + CARRY_CORRECTION(state)")
print()
print("The carry correction depends on state through:")
print("  carry(x, Dx) = (x+Dx) ⊕ x ⊕ Dx")
print("  This is a DETERMINISTIC function of (x, Dx)")
print("  With geometric distribution: P(HW(carry)=k) = (1/2)^k")
print()
print("NEXT STEPS:")
print("  1. The 2·(L(x)∧L(c)) term is the ONLY barrier-creating term")
print("  2. It involves AND — this is the same nonlinearity as Ch/Maj!")
print("  3. ALL nonlinearity in SHA-256 reduces to AND operations")
print("  4. Can we build an algebra where AND is tractable?")
print()

# Verify: how many AND operations per round?
print("AND operations per round:")
print("  Ch(e,f,g):  2 ANDs (e∧f, ¬e∧g)")
print("  Maj(a,b,c): 3 ANDs (a∧b, a∧c, b∧c)")
print("  Sig1 carry correction: 1 AND (Sig1(e) ∧ Sig1(carry))")
print("  Sig0 carry correction: 1 AND (Sig0(a) ∧ Sig0(carry))")
print("  Addition carries (T1 = h+Sig1+Ch+K+W): ~4 carry chains")
print("  Addition carries (T2 = Sig0+Maj): ~1 carry chain")
print()
print("  Total AND-type operations per round: ~12")
print("  Over 64 rounds: ~768 AND operations")
print("  This is the TRUE complexity of SHA-256")
