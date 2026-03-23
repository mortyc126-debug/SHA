#!/usr/bin/env python3
"""
Step 3b: Fix the identity. Find the CORRECT algebraic structure.

Error in Step 3: (A XOR B) - A ≠ B - 2(A∧B) in general over Z/2^32Z.
Actually: A XOR B = A + B - 2(A∧B), so (A XOR B) - A = B - 2(A∧B) IS correct!

But the ERROR was: Sig1(e+de) ≠ Sig1(e) XOR Sig1(carry(e,de))
Because Sig1 involves ROTATIONS, and e+de = e XOR carry(e,de)
only when the XOR and carry don't overlap.

Let me re-derive carefully:
  e + de = e XOR de XOR 2·carry_bits(e, de)
  where carry_bits(e,de) = the actual carry chain

So: Sig1(e+de) = Sig1(e XOR de XOR 2·carry_bits)

Since Sig1 is GF(2)-linear on individual bits (rotation+XOR):
  Sig1(a XOR b) = Sig1(a) XOR Sig1(b)  ALWAYS (GF2 linearity)

But: Sig1(a + b) ≠ Sig1(a) + Sig1(b) because + and rotation don't commute.

NEW APPROACH: Don't try to make additive identity.
Instead, DIRECTLY measure the algebraic structure of the round map.
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
     0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174]

# ============================================================
# PART A: The CORRECT identity for XOR-linear functions
# ============================================================
print("=" * 70)
print("PART A: Correct identity verification")
print("=" * 70)
print()

# Identity: For any A, B: A XOR B = A + B - 2·(A AND B)  (mod 2^32)
# Verify this first:
ok = 0
for _ in range(10000):
    A = random.randint(0, MASK)
    B = random.randint(0, MASK)
    lhs = A ^ B
    rhs = sub(add(A, B), (2 * (A & B)) & MASK)
    if lhs == rhs:
        ok += 1
print(f"A⊕B = A+B-2(A∧B) mod 2^32: {ok}/10000 {'✓' if ok==10000 else '✗'}")

# So: Sig1(e+de) = Sig1(e ⊕ de ⊕ 2·c(e,de))  where c = carry bits
# But: Sig1 involves ROTATIONS which don't commute with multiplication by 2.
# rotr(2x, n) ≠ 2·rotr(x, n) because bit-rotation ≠ multiplication.

# Therefore: we need to work with the EXACT XOR identity:
# e + de = e ⊕ (carry_chain(e, de))
# where carry_chain(e, de) = (e+de) ⊕ e ⊕ de... wait, that's wrong too.
# Actually: (e + de) = e XOR de XOR (2·carries)
# But 2·carries shifts carry bits left by 1.

# Let me just verify the CORRECT decomposition:
# carry_xor(e, de) = (e + de) XOR e
# This is the total XOR effect of adding de to e.
# Then: Sig1(e+de) = Sig1(e XOR carry_xor(e,de))
#                   = Sig1(e) XOR Sig1(carry_xor(e,de))  [GF2 linearity]

ok = 0
for _ in range(10000):
    e = random.randint(0, MASK)
    de = random.randint(0, MASK)
    carry_xor = add(e, de) ^ e  # = de XOR 2·carries
    lhs = Sig1(add(e, de))
    rhs = Sig1(e) ^ Sig1(carry_xor)
    if lhs == rhs:
        ok += 1
print(f"Sig1(e+de) = Sig1(e) ⊕ Sig1((e+de)⊕e): {ok}/10000 {'✓' if ok==10000 else '✗'}")
print()

# GREAT! This DOES work. Because:
# Let c = (e+de) XOR e. Then e+de = e XOR c.
# Sig1(e XOR c) = Sig1(e) XOR Sig1(c)  [GF2 linearity]  ✓
#
# Now for the ADDITIVE differential:
# D_add[Sig1](e, de) = Sig1(e+de) - Sig1(e)
#                     = [Sig1(e) XOR Sig1(c)] - Sig1(e)
#                     = Sig1(c) - 2·(Sig1(e) AND Sig1(c))  [using A⊕B-A = B-2(A∧B)]

ok = 0
for _ in range(10000):
    e = random.randint(0, MASK)
    de = random.randint(0, MASK)
    c = add(e, de) ^ e  # carry_xor
    actual = sub(Sig1(add(e, de)), Sig1(e))
    A = Sig1(e)
    B = Sig1(c)
    predicted = sub(B, (2 * (A & B)) & MASK)
    if actual == predicted:
        ok += 1
print(f"D_add[Sig1] = Sig1(c) - 2(Sig1(e)∧Sig1(c)): {ok}/10000 {'✓' if ok==10000 else '✗'}")
print()

# ============================================================
# PART B: The CORRECTED fundamental identity — now for all functions
# ============================================================
print("=" * 70)
print("PART B: Universal carry-aware identity (CORRECTED)")
print("=" * 70)
print()
print("For any GF(2)-linear function L and Z-addition x+y:")
print("  D_add[L](x, y) = L(x+y) - L(x)")
print("                  = L(c) - 2·(L(x) ∧ L(c))")
print("  where c = (x+y) ⊕ x  (the XOR-effect of adding y to x)")
print()

# c = (x+y) XOR x. Note: c ≠ y in general! c = y when there are no carries.
# c depends on BOTH x and y — this is the state-dependence.

for name, L in [("Sig0", Sig0), ("Sig1", Sig1), ("sig0", sig0), ("sig1", sig1)]:
    ok = 0
    for _ in range(5000):
        x = random.randint(0, MASK)
        y = random.randint(0, MASK)
        c = add(x, y) ^ x  # NOTE: c = (x+y)⊕x, NOT (x+y)⊕x⊕y !
        actual = sub(L(add(x, y)), L(x))
        A = L(x)
        B = L(c)
        pred = sub(B, (2 * (A & B)) & MASK)
        if actual == pred:
            ok += 1
    print(f"  {name}: {ok}/5000 {'✓ EXACT' if ok==5000 else '✗'}")

print()
print("★ FUNDAMENTAL IDENTITY CONFIRMED FOR ALL 4 SIGMA FUNCTIONS ★")
print()

# ============================================================
# PART C: What IS carry_xor c = (x+y)⊕x?
# ============================================================
print("=" * 70)
print("PART C: Structure of c = (x+y)⊕x")
print("=" * 70)
print()

# c = (x+y) XOR x
# Note: x + y = x XOR y XOR 2·carry(x,y)
# So: c = (x+y) XOR x = y XOR 2·carry(x,y)
#
# carry(x,y) at bit i = MAJ(x[i], y[i], carry[i-1])
# 2·carry = carry shifted left by 1
#
# So: c[0] = y[0]  (no carry into bit 0)
# c[i] = y[i] XOR carry[i-1]  for i >= 1
# c[32] = carry[31]  (but we're mod 2^32, so this is lost)

# Verify c = y XOR 2·carry(x,y):
print("c = y ⊕ 2·carry(x,y):")
ok = 0
for _ in range(10000):
    x = random.randint(0, MASK)
    y = random.randint(0, MASK)
    c = add(x, y) ^ x
    # carry = ((x+y) XOR x XOR y) >> 1... no.
    # (x + y) = x XOR y XOR (carry_chain << 1)... approximately
    # Actually: x + y = x XOR y XOR 2·carry_bits
    # where carry_bits[i] = MAJ(x[i], y[i], carry_bits[i-1])
    # So: (x+y) XOR x = y XOR 2·carry_bits
    # Therefore: c = y XOR 2·carry_bits
    # carry_bits = ((x+y) XOR x XOR y) >> 1  ... let's check:
    carry_2 = (add(x, y) ^ x ^ y)  # = 2·carry_bits
    c_pred = y ^ carry_2
    if c == c_pred:
        ok += 1
print(f"  Verification c = y ⊕ ((x+y)⊕x⊕y): {ok}/10000 {'✓' if ok==10000 else '✗'}")
print()

# Distribution of HW(c - y) = HW(2·carry)
print("Distribution of carry effect |c - y| = |2·carry(x,y)|:")
carry_hw = []
for _ in range(50000):
    x = random.randint(0, MASK)
    y = 1  # minimal y
    c = add(x, y) ^ x
    carry_eff = c ^ y  # = 2·carry
    carry_hw.append(hw(carry_eff))

print(f"  y=1: E[HW(carry_effect)] = {np.mean(carry_hw):.2f}")
print(f"  Distribution:")
for k in range(8):
    cnt = carry_hw.count(k)
    print(f"    HW={k}: P={cnt/50000:.4f}  (≈ {cnt/50000*100:.1f}%)")

print()

# ============================================================
# PART D: The FULL round identity
# ============================================================
print("=" * 70)
print("PART D: Full SHA-256 round in carry-aware algebra")
print("=" * 70)
print()

# One round: given state (a,b,c,d,e,f,g,h), W, K:
#   T1 = h + Sig1(e) + Ch(e,f,g) + K + W
#   T2 = Sig0(a) + Maj(a,b,c)
#   a' = T1 + T2
#   e' = d + T1
#
# For the ADDITIVE differential (D = state' - state for perturbed input):
#   DT1 = Dh + D_add[Sig1](e, De) + D_add[Ch](e,f,g, De,Df,Dg) + DW
#   DT2 = D_add[Sig0](a, Da) + D_add[Maj](a,b,c, Da,Db,Dc)
#   Da' = DT1 + DT2  (EXACT in Z/2^32Z — addition IS linear!)
#   De' = Dd + DT1   (EXACT!)
#
# Using our identity:
#   D_add[Sig1](e,De) = Sig1(c_e) - 2(Sig1(e)∧Sig1(c_e))
#     where c_e = (e+De)⊕e = De ⊕ 2·carry(e,De)
#
# The FULL differential is:
#   Da' = Dh + [Sig1(c_e) - 2(Sig1(e)∧Sig1(c_e))]
#        + [D_add_Ch] + DW
#        + [Sig0(c_a) - 2(Sig0(a)∧Sig0(c_a))]
#        + [D_add_Maj]

# Now verify the FULL round identity:
def sha_one_round(state, W, K_r):
    a,b,c,d,e,f,g,h = state
    T1 = add(add(add(add(h, Sig1(e)), Ch(e,f,g)), K_r), W)
    T2 = add(Sig0(a), Maj(a,b,c))
    return [add(T1, T2), a, b, c, add(d, T1), e, f, g]

N5 = 5000
exact_round = 0
for _ in range(N5):
    state = [random.randint(0, MASK) for _ in range(8)]
    Dstate = [random.randint(0, MASK) for _ in range(8)]
    W = random.randint(0, MASK)

    state2 = [add(state[i], Dstate[i]) for i in range(8)]

    out1 = sha_one_round(state, W, K[0])
    out2 = sha_one_round(state2, W, K[0])

    # Actual additive differential
    Dout_actual = [sub(out2[i], out1[i]) for i in range(8)]

    # Predicted: Da' = DT1 + DT2, De' = Dd + DT1
    a,b,c,d,e,f,g,h = state
    Da,Db,Dc,Dd,De,Df,Dg,Dh = Dstate

    # D_add[Sig1] using identity
    c_e = add(e, De) ^ e
    DSig1 = sub(Sig1(c_e), (2 * (Sig1(e) & Sig1(c_e))) & MASK)

    # D_add[Ch] — direct computation
    DCh = sub(Ch(add(e,De), add(f,Df), add(g,Dg)), Ch(e,f,g))

    # D_add[Sig0] using identity
    c_a = add(a, Da) ^ a
    DSig0 = sub(Sig0(c_a), (2 * (Sig0(a) & Sig0(c_a))) & MASK)

    # D_add[Maj]
    DMaj = sub(Maj(add(a,Da), add(b,Db), add(c,Dc)), Maj(a,b,c))

    DT1 = add(add(add(Dh, DSig1), DCh), 0)  # DW = 0 (same W)
    DT2 = add(DSig0, DMaj)

    Dout_pred = [add(DT1, DT2), Da, Db, Dc, add(Dd, DT1), De, Df, Dg]

    if all(Dout_actual[i] == Dout_pred[i] for i in range(8)):
        exact_round += 1

print(f"Full round differential identity: {exact_round}/{N5} {'✓ EXACT' if exact_round==N5 else '✗'}")
print()

if exact_round == N5:
    print("★★★ SHA-256 ROUND DIFFERENTIAL IDENTITY CONFIRMED ★★★")
    print()
    print("Da' = Dh + [Sig1(c_e) - 2(Sig1(e)∧Sig1(c_e))] + DCh + DW")
    print("    + [Sig0(c_a) - 2(Sig0(a)∧Sig0(c_a))] + DMaj")
    print()
    print("De' = Dd + Dh + [Sig1(c_e) - 2(Sig1(e)∧Sig1(c_e))] + DCh + DW")
    print()
    print("where:")
    print("  c_e = (e+De)⊕e       [carry-XOR of e-register]")
    print("  c_a = (a+Da)⊕a       [carry-XOR of a-register]")
    print("  DCh, DMaj computed directly (quadratic in state)")
    print()
    print("NONLINEAR TERMS (all involving AND):")
    print("  1. Sig1(e) ∧ Sig1(c_e)    [carry correction for Sig1]")
    print("  2. Sig0(a) ∧ Sig0(c_a)    [carry correction for Sig0]")
    print("  3. Ch(e+De,f+Df,g+Dg)-Ch(e,f,g)  [Ch nonlinearity]")
    print("  4. Maj(a+Da,b+Db,c+Dc)-Maj(a,b,c) [Maj nonlinearity]")
    print()
    print("ALL 4 terms are products of GF(2)-linear functions!")
    print("This is the algebraic basis for Clifford-like structure.")

print()
print("=" * 70)
print("SUMMARY STEP 3b")
print("=" * 70)
print()
print("CONFIRMED IDENTITIES:")
print("  1. A⊕B = A + B - 2(A∧B)  mod 2^32")
print("  2. Sig1(e+de) = Sig1(e) ⊕ Sig1((e+de)⊕e)  [GF2 linearity]")
print("  3. D_add[L](x,y) = L(c) - 2(L(x)∧L(c))  where c=(x+y)⊕x")
print("  4. Full round identity: Da', De' decompose EXACTLY")
print()
print("ALGEBRA STRUCTURE:")
print("  - Linear part: additions, shift registers (EXACT)")
print("  - GF2-linear part: Sig0, Sig1 applied to carry patterns")
print("  - Quadratic part: 4 AND terms per round")
print("  - ALL nonlinearity = AND operations on GF2-linear images")
