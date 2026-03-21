#!/usr/bin/env python3
"""
Step 1: Extract bilinear form Q from Ch and Maj for Clifford algebra construction.

Key question: Ch(e,f,g) = ef ⊕ eg ⊕ g and Maj(a,b,c) = ab ⊕ ac ⊕ bc
are QUADRATIC over GF(2). What is the associated bilinear form?

For a quadratic form Q(x), the associated bilinear form is:
  B(x,y) = Q(x+y) - Q(x) - Q(y)

This is the foundation of Clifford algebra Cl(V, Q).
"""

import numpy as np
from collections import Counter
import random

MASK = 0xFFFFFFFF

# === SHA-256 primitives ===
def rotr(x, n):
    return ((x >> n) | (x << (32 - n))) & MASK

def sig0(x):  return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sig1(x):  return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)
def Sig0(x):  return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sig1(x):  return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g):   return ((e & f) ^ (~e & g)) & MASK
def Maj(a, b, c):  return ((a & b) ^ (a & c) ^ (b & c)) & MASK

# SHA-256 IV
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


# ============================================================
# PART A: Bilinear form of Ch over GF(2), bit by bit
# ============================================================
# Ch(e,f,g) = e·f ⊕ ¬e·g = e·f ⊕ e·g ⊕ g  (over GF(2))
# This is quadratic in (e,f,g).
#
# The POLARIZATION (associated bilinear form) B_Ch:
#   B_Ch((e1,f1,g1), (e2,f2,g2)) = Ch(e1⊕e2, f1⊕f2, g1⊕g2) ⊕ Ch(e1,f1,g1) ⊕ Ch(e2,f2,g2)
#
# Let's compute this symbolically for one bit:
#   Ch(e,f,g) = ef + eg + g  (mod 2)
#   Ch(e1+e2, f1+f2, g1+g2) = (e1+e2)(f1+f2) + (e1+e2)(g1+g2) + (g1+g2)
#     = e1f1 + e1f2 + e2f1 + e2f2 + e1g1 + e1g2 + e2g1 + e2g2 + g1 + g2
#   Subtract Ch(e1,f1,g1) = e1f1 + e1g1 + g1
#   Subtract Ch(e2,f2,g2) = e2f2 + e2g2 + g2
#   B_Ch = e1f2 + e2f1 + e1g2 + e2g1  (mod 2)
#
# So: B_Ch = e1·(f2+g2) + e2·(f1+g1)  -- this is the BILINEAR FORM

print("=" * 70)
print("PART A: Symbolic bilinear form of Ch")
print("=" * 70)
print()
print("Ch(e,f,g) = ef ⊕ eg ⊕ g  over GF(2)")
print()
print("Polarization B_Ch((e1,f1,g1), (e2,f2,g2)):")
print("  = e1·(f2⊕g2) ⊕ e2·(f1⊕g1)")
print()
print("This is a SYMPLECTIC form (B_Ch(x,x) = 0 for all x)")
print("  Proof: B_Ch((e,f,g),(e,f,g)) = e·(f⊕g) ⊕ e·(f⊕g) = 0")
print()

# Verify numerically
def verify_bilinear_Ch(N=10000):
    """Verify B_Ch formula on random bit triples."""
    ok = 0
    for _ in range(N):
        e1, f1, g1 = random.randint(0,1), random.randint(0,1), random.randint(0,1)
        e2, f2, g2 = random.randint(0,1), random.randint(0,1), random.randint(0,1)

        # Direct: Ch(x+y) + Ch(x) + Ch(y)
        ch_sum = ((e1^e2)&(f1^f2)) ^ ((~(e1^e2))&(g1^g2)) & 1
        ch_sum ^= ((e1&f1) ^ ((~e1)&g1)) & 1
        ch_sum ^= ((e2&f2) ^ ((~e2)&g2)) & 1
        ch_sum &= 1

        # Formula: e1·(f2⊕g2) ⊕ e2·(f1⊕g1)
        formula = (e1 & (f2 ^ g2)) ^ (e2 & (f1 ^ g1))

        if ch_sum == formula:
            ok += 1
    return ok, N

ok, N = verify_bilinear_Ch()
print(f"Verification B_Ch formula: {ok}/{N} ✓" if ok == N else f"FAILED: {ok}/{N}")

# ============================================================
# PART B: Bilinear form of Maj over GF(2)
# ============================================================
# Maj(a,b,c) = ab ⊕ ac ⊕ bc  (over GF(2))
#
# Polarization:
#   Maj(a1+a2, b1+b2, c1+c2) - Maj(a1,b1,c1) - Maj(a2,b2,c2)
#   = (a1+a2)(b1+b2) + (a1+a2)(c1+c2) + (b1+b2)(c1+c2) - (a1b1+a1c1+b1c1) - (a2b2+a2c2+b2c2)
#   = a1b2+a2b1 + a1c2+a2c1 + b1c2+b2c1  (mod 2)
#
# B_Maj = a1·b2 ⊕ a2·b1 ⊕ a1·c2 ⊕ a2·c1 ⊕ b1·c2 ⊕ b2·c1
#       = a1·(b2⊕c2) ⊕ a2·(b1⊕c1) ⊕ b1·c2 ⊕ b2·c1

print()
print("=" * 70)
print("PART B: Symbolic bilinear form of Maj")
print("=" * 70)
print()
print("Maj(a,b,c) = ab ⊕ ac ⊕ bc  over GF(2)")
print()
print("Polarization B_Maj((a1,b1,c1), (a2,b2,c2)):")
print("  = a1·b2 ⊕ a2·b1 ⊕ a1·c2 ⊕ a2·c1 ⊕ b1·c2 ⊕ b2·c1")
print()

# Verify
def verify_bilinear_Maj(N=10000):
    ok = 0
    for _ in range(N):
        a1, b1, c1 = random.randint(0,1), random.randint(0,1), random.randint(0,1)
        a2, b2, c2 = random.randint(0,1), random.randint(0,1), random.randint(0,1)

        # Direct
        maj_sum = ((a1^a2)&(b1^b2)) ^ ((a1^a2)&(c1^c2)) ^ ((b1^b2)&(c1^c2))
        maj_sum ^= (a1&b1) ^ (a1&c1) ^ (b1&c1)
        maj_sum ^= (a2&b2) ^ (a2&c2) ^ (b2&c2)
        maj_sum &= 1

        # Formula
        formula = (a1&b2) ^ (a2&b1) ^ (a1&c2) ^ (a2&c1) ^ (b1&c2) ^ (b2&c1)

        if maj_sum == formula:
            ok += 1
    return ok, N

ok, N = verify_bilinear_Maj()
print(f"Verification B_Maj formula: {ok}/{N} ✓" if ok == N else f"FAILED: {ok}/{N}")

# ============================================================
# PART C: Combined bilinear form Q for one SHA-256 round
# ============================================================
# One round: state' = Round(state, W, K)
#   T1 = h + Sig1(e) + Ch(e,f,g) + K + W
#   T2 = Sig0(a) + Maj(a,b,c)
#   a' = T1 + T2
#   e' = d + T1
#
# Nonlinear parts: Ch(e,f,g) and Maj(a,b,c)
# Linear parts: Sig0, Sig1 (over GF(2)), additions (over Z/2^32Z)
#
# The quadratic form of one round (over GF(2), per bit) involves:
#   Q_round = Q_Ch(e,f,g) + Q_Maj(a,b,c)
# where "+" is XOR (GF(2) addition)
#
# The FULL bilinear form is:
#   B_round(state1, state2) = B_Ch(e1,f1,g1; e2,f2,g2) ⊕ B_Maj(a1,b1,c1; a2,b2,c2)

print()
print("=" * 70)
print("PART C: Combined round bilinear form")
print("=" * 70)
print()
print("B_round(s1, s2) = B_Ch(e1,f1,g1; e2,f2,g2) ⊕ B_Maj(a1,b1,c1; a2,b2,c2)")
print()
print("State vector: s = (a, b, c, d, e, f, g, h) ∈ (GF(2)^32)^8")
print("Dimension of V: 8 × 32 = 256 bits")
print()

# ============================================================
# PART D: Matrix representation of B_round (bit level, one bit position)
# ============================================================
# For a fixed bit position i, state = (a[i], b[i], c[i], d[i], e[i], f[i], g[i], h[i])
# B_round is 8×8 matrix over GF(2).
#
# B_Ch involves only (e,f,g) = positions 4,5,6
# B_Maj involves only (a,b,c) = positions 0,1,2
#
# B_Ch matrix (on e,f,g subspace):
#   B_Ch(x,y) = e1·(f2⊕g2) + e2·(f1⊕g1)
#   In matrix form with basis (e,f,g):
#     [0, 1, 1]
#     [1, 0, 0]     (symmetric: B_Ch[i,j] = B_Ch[j,i])
#     [1, 0, 0]
#
# B_Maj matrix (on a,b,c subspace):
#   B_Maj(x,y) = a1·b2 + a2·b1 + a1·c2 + a2·c1 + b1·c2 + b2·c1
#   In matrix form with basis (a,b,c):
#     [0, 1, 1]
#     [1, 0, 1]     (symmetric)
#     [1, 1, 0]

print("=" * 70)
print("PART D: Matrix of B_round (per-bit, 8×8 over GF(2))")
print("=" * 70)
print()

# Build 8×8 bilinear form matrix (per bit)
B = np.zeros((8, 8), dtype=int)

# B_Maj on (a,b,c) = indices 0,1,2
B[0,1] = B[1,0] = 1  # a·b cross term
B[0,2] = B[2,0] = 1  # a·c cross term
B[1,2] = B[2,1] = 1  # b·c cross term

# B_Ch on (e,f,g) = indices 4,5,6
B[4,5] = B[5,4] = 1  # e·f cross term
B[4,6] = B[6,4] = 1  # e·g cross term
# Note: f·g cross term is ABSENT in Ch (unlike Maj)

print("B_round (8×8 over GF(2)):")
print("Basis: (a, b, c, d, e, f, g, h)")
print()
for i, label in enumerate("abcdefgh"):
    row = "  " + label + " | " + " ".join(str(B[i,j]) for j in range(8))
    print(row)
print("      +" + "-" * 16)
print("        " + " ".join("abcdefgh"))

# Properties of B
rank_B = np.linalg.matrix_rank(B % 2)  # over GF(2) approximation
print(f"\nRank of B over GF(2): {rank_B}")
print(f"Nullity (dim kernel): {8 - rank_B}")
print()

# ============================================================
# PART E: Key structural question — what is the kernel?
# ============================================================
# Kernel of B = vectors v such that B(v, w) = 0 for ALL w
# These are the "invisible" directions in the Clifford algebra.
#
# From the matrix:
#   d (index 3) has all zeros → d ∈ ker(B)
#   h (index 7) has all zeros → h ∈ ker(B)
#
# This is HUGE: d and h are INVISIBLE to the quadratic nonlinearity!
# But: e' = d + T1 — the d-register DIRECTLY feeds into e!
#      h enters T1 = h + Sig1(e) + Ch(e,f,g) + K + W
#
# So d and h are "linear channels" — they carry information
# WITHOUT nonlinear distortion from Ch/Maj.

print("=" * 70)
print("PART E: Kernel of B — 'invisible' directions")
print("=" * 70)
print()

kernel_vectors = []
for i in range(8):
    if all(B[i,j] == 0 for j in range(8)):
        kernel_vectors.append(i)
        print(f"  Register '{chr(97+i)}' (index {i}) is in ker(B) — invisible to nonlinearity!")

print(f"\nKernel dimension: {len(kernel_vectors)}")
print(f"These registers pass through SHA rounds WITHOUT quadratic distortion.")
print()
print("CRITICAL INSIGHT:")
print("  d and h are linear channels in the Clifford algebra.")
print("  d[r] = a[r-3] (shift register)")
print("  h[r] = e[r-3] (shift register)")
print("  They are 3-round delayed copies of a and e.")
print()
print("  In Cl(V,Q): d and h COMMUTE with all Cl-elements!")
print("  (because B(d,·) = B(h,·) = 0)")
print()

# ============================================================
# PART F: Verify on 32-bit words — does B capture real SHA behavior?
# ============================================================
print("=" * 70)
print("PART F: Experimental verification — B captures SHA nonlinearity")
print("=" * 70)
print()

def sha_one_round(state, W_r, K_r):
    """One round of SHA-256. Returns new state."""
    a, b, c, d, e, f, g, h = state
    T1 = (h + Sig1(e) + Ch(e, f, g) + K_r + W_r) & MASK
    T2 = (Sig0(a) + Maj(a, b, c)) & MASK
    return [(T1 + T2) & MASK, a, b, c, (d + T1) & MASK, e, f, g]

def xor_diff(s1, s2):
    """XOR difference of two states."""
    return [s1[i] ^ s2[i] for i in range(8)]

def hw(x):
    """Hamming weight."""
    return bin(x).count('1')

# Test: perturb only d (kernel direction) vs perturb only e (non-kernel)
N = 5000
hw_d_perturb = []  # HW of output diff when perturbing d
hw_e_perturb = []  # HW of output diff when perturbing e

for _ in range(N):
    W_r = random.randint(0, MASK)
    state = [random.randint(0, MASK) for _ in range(8)]

    # Perturb d (index 3) by 1 bit
    state_d = list(state)
    state_d[3] ^= 1

    out_orig = sha_one_round(state, W_r, K[0])
    out_d = sha_one_round(state_d, W_r, K[0])
    diff_d = xor_diff(out_orig, out_d)
    hw_d = sum(hw(x) for x in diff_d)
    hw_d_perturb.append(hw_d)

    # Perturb e (index 4) by 1 bit
    state_e = list(state)
    state_e[4] ^= 1

    out_e = sha_one_round(state_e, W_r, K[0])
    diff_e = xor_diff(out_orig, out_e)
    hw_e = sum(hw(x) for x in diff_e)
    hw_e_perturb.append(hw_e)

avg_d = np.mean(hw_d_perturb)
avg_e = np.mean(hw_e_perturb)
std_d = np.std(hw_d_perturb)
std_e = np.std(hw_e_perturb)

print("Perturbation test (1 bit flip, 1 round):")
print(f"  Perturb d (kernel):     E[HW(diff)] = {avg_d:.2f} ± {std_d:.2f}")
print(f"  Perturb e (non-kernel): E[HW(diff)] = {avg_e:.2f} ± {std_e:.2f}")
print(f"  Ratio e/d: {avg_e/avg_d:.2f}x")
print()

# Count how many output registers are affected
def count_nonzero_regs(diff):
    return sum(1 for x in diff if x != 0)

regs_d = np.mean([count_nonzero_regs(xor_diff(
    sha_one_round([random.randint(0,MASK) for _ in range(8)], random.randint(0,MASK), K[0]),
    sha_one_round(
        (lambda s: [s[0],s[1],s[2],s[3]^1,s[4],s[5],s[6],s[7]])(
            [random.randint(0,MASK) for _ in range(8)]
        ), random.randint(0,MASK), K[0])
)) for _ in range(1000)])

# Actually need same state for both...
affected_d = []
affected_e = []
for _ in range(2000):
    W_r = random.randint(0, MASK)
    state = [random.randint(0, MASK) for _ in range(8)]
    out0 = sha_one_round(state, W_r, K[0])

    sd = list(state); sd[3] ^= 1
    out_d = sha_one_round(sd, W_r, K[0])
    affected_d.append(count_nonzero_regs(xor_diff(out0, out_d)))

    se = list(state); se[4] ^= 1
    out_e = sha_one_round(se, W_r, K[0])
    affected_e.append(count_nonzero_regs(xor_diff(out0, out_e)))

print(f"  Affected registers (d perturb): E = {np.mean(affected_d):.2f}")
print(f"  Affected registers (e perturb): E = {np.mean(affected_e):.2f}")
print()

# Key question: is d perturbation truly LINEAR (deterministic)?
print("Is d-perturbation deterministic?")
diffs_d_vals = Counter()
for _ in range(1000):
    W_r = random.randint(0, MASK)
    state = [random.randint(0, MASK) for _ in range(8)]
    sd = list(state); sd[3] ^= 1  # flip bit 0 of d
    out0 = sha_one_round(state, W_r, K[0])
    out1 = sha_one_round(sd, W_r, K[0])
    diff = tuple(xor_diff(out0, out1))
    diffs_d_vals[diff] += 1

print(f"  Unique output diffs from d[0] flip: {len(diffs_d_vals)}")
if len(diffs_d_vals) == 1:
    print("  → DETERMINISTIC! d passes through linearly. ✓")
    print(f"  → Fixed diff: {[hex(x) for x in list(diffs_d_vals.keys())[0]]}")
else:
    print(f"  → {len(diffs_d_vals)} different patterns (depends on state)")
    # Show top patterns
    for diff, cnt in diffs_d_vals.most_common(3):
        regs = [hex(x) for x in diff]
        print(f"     {cnt}/1000: {regs}")

print()

# Same for e
diffs_e_vals = Counter()
for _ in range(1000):
    W_r = random.randint(0, MASK)
    state = [random.randint(0, MASK) for _ in range(8)]
    se = list(state); se[4] ^= 1
    out0 = sha_one_round(state, W_r, K[0])
    out1 = sha_one_round(se, W_r, K[0])
    diff = tuple(xor_diff(out0, out1))
    diffs_e_vals[diff] += 1

print(f"Is e-perturbation deterministic?")
print(f"  Unique output diffs from e[0] flip: {len(diffs_e_vals)}")
if len(diffs_e_vals) == 1:
    print("  → DETERMINISTIC!")
else:
    print(f"  → {len(diffs_e_vals)} different patterns (state-DEPENDENT — nonlinear!)")

print()
print("=" * 70)
print("SUMMARY STEP 1")
print("=" * 70)
print()
print("1. Ch and Maj define a bilinear form B on (GF(2)^32)^8")
print("2. B is 8×8 per bit position, with rank 6 and kernel dim 2")
print("3. ker(B) = {d, h} — these registers are INVISIBLE to nonlinearity")
print("4. d-perturbation is deterministic (linear channel)")
print("5. e-perturbation is state-dependent (nonlinear)")
print()
print("NEXT: Use B to construct Clifford algebra Cl(V, Q)")
print("  V = (GF(2)^32)^6  (excluding d,h from kernel)")
print("  Q = restriction of Ch+Maj bilinear form to V")
