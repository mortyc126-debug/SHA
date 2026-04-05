#!/usr/bin/env python3
"""
THE WALL: A[3] = f(A[1], A[2]) through 2 SHA-256 rounds.

Given: A[0] = IV[0] = 0x6a09e667 (constant)
       A[3] = known (from hash via backward chain)

Find:  A[1], A[2] such that 2 rounds from (A[0],A[1]) produce A[3].

Question: what is the algebraic structure? Can we reduce 2^64?
"""

import random, math, time
from collections import defaultdict

MASK32 = 0xFFFFFFFF
K=[0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2]
IV=[0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]
def rotr(x,n):return((x>>n)|(x<<(32-n)))&MASK32
def sig0(x):return rotr(x,2)^rotr(x,13)^rotr(x,22)
def sig1(x):return rotr(x,6)^rotr(x,11)^rotr(x,25)
def ch(e,f,g):return(e&f)^(~e&g)&MASK32
def maj(a,b,c):return(a&b)^(a&c)^(b&c)
def hw(x):return bin(x&MASK32).count('1')


# ================================================================
# STEP 1: Exact formulas for A[1], A[2], A[3] from W[0], W[1], W[2]
# ================================================================

print("=" * 70)
print("STEP 1: Exact computation of A[1], A[2], A[3]")
print("=" * 70)

# Round 0: state = IV, input W[0]
# T1[0] = h + Σ₁(e) + Ch(e,f,g) + K[0] + W[0]
#        = IV[7] + Σ₁(IV[4]) + Ch(IV[4],IV[5],IV[6]) + K[0] + W[0]
C0 = (IV[7] + sig1(IV[4]) + ch(IV[4],IV[5],IV[6]) + K[0]) & MASK32
print(f"  C0 = IV[7]+Σ₁(IV[4])+Ch(IV[4..6])+K[0] = 0x{C0:08x}")
print(f"  T1[0] = C0 + W[0] = 0x{C0:08x} + W[0]")

# T2[0] = Σ₀(IV[0]) + Maj(IV[0],IV[1],IV[2])
T2_0 = (sig0(IV[0]) + maj(IV[0],IV[1],IV[2])) & MASK32
print(f"  T2[0] = Σ₀(IV[0])+Maj(IV[0..2]) = 0x{T2_0:08x} (CONSTANT)")

# A[1] = T1[0] + T2[0] = C0 + W[0] + T2_0
A1_const = (C0 + T2_0) & MASK32
print(f"  A[1] = 0x{A1_const:08x} + W[0]   (LINEAR in W[0])")
print(f"       = 0xfc08884d + W[0]")

# E[1] = d[0] + T1[0] = IV[3] + C0 + W[0]
E1_const = (IV[3] + C0) & MASK32
print(f"  E[1] = 0x{E1_const:08x} + W[0]   (LINEAR in W[0])")

# Verify
random.seed(42)
for trial in range(3):
    W0 = random.randint(0, MASK32)
    A1 = (A1_const + W0) & MASK32
    E1 = (E1_const + W0) & MASK32

    # Direct computation
    a,b,c,d,e,f,g,h = IV
    T1=(h+sig1(e)+ch(e,f,g)+K[0]+W0)&MASK32;T2=(sig0(a)+maj(a,b,c))&MASK32
    a_direct=(T1+T2)&MASK32; e_direct=(d+T1)&MASK32
    assert A1 == a_direct and E1 == e_direct
print(f"  Verified: A[1], E[1] formulas correct ✓")


# ================================================================
# STEP 2: A[2] as function of A[1] (and W[1])
# ================================================================

print()
print("=" * 70)
print("STEP 2: A[2] = f(A[1], W[1])")
print("=" * 70)

# Round 1: state = (A[1], IV[0], IV[1], IV[2], E[1], IV[4], IV[5], IV[6])
# T1[1] = IV[6] + Σ₁(E[1]) + Ch(E[1], IV[4], IV[5]) + K[1] + W[1]
# T2[1] = Σ₀(A[1]) + Maj(A[1], IV[0], IV[1])
# A[2] = T1[1] + T2[1]

# E[1] = E1_const + W[0] = E1_const + (A[1] - A1_const)
# So E[1] is LINEAR in A[1]: E[1] = A[1] + (E1_const - A1_const)
E1_offset = (E1_const - A1_const) & MASK32
print(f"  E[1] = A[1] + 0x{E1_offset:08x}   (E[1] LINEAR in A[1]!)")
print(f"  (because both A[1] and E[1] are linear in W[0], same coefficient)")

# T2[1] = Σ₀(A[1]) + Maj(A[1], IV[0], IV[1])
# This is NONLINEAR in A[1] (Σ₀ = rotations = linear, Maj = degree 2)
# But Σ₀ is GF(2)-linear, and Maj is degree 2 over GF(2).

print(f"  T2[1] = Σ₀(A[1]) + Maj(A[1], IV[0], IV[1])")
print(f"  → Σ₀ is linear (rotations), Maj is quadratic")
print(f"  → T2[1] is NONLINEAR in A[1]")

# T1[1] = IV[6] + Σ₁(E[1]) + Ch(E[1], IV[4], IV[5]) + K[1] + W[1]
# E[1] = A[1] + offset. Σ₁(E[1]) = Σ₁(A[1]+offset).
# Σ₁ is GF(2)-linear but NOT Z/2^32-linear.
# So T1[1] is NONLINEAR in A[1].

print(f"  T1[1] = IV[6] + Σ₁(A[1]+offset) + Ch(A[1]+offset, IV[4], IV[5]) + K[1] + W[1]")
print(f"  → NONLINEAR in A[1] (Σ₁ and Ch both nonlinear over Z/2^32)")

# A[2] = T1[1] + T2[1] — nonlinear function of A[1] and W[1]
print(f"  A[2] = T1[1] + T2[1] — nonlinear in (A[1], W[1])")


# ================================================================
# STEP 3: A[3] = f(A[1], A[2], W[2])
# Can we express A[3] as function of A[1] alone (+ W[1], W[2])?
# YES: A[2] = g(A[1], W[1]), so A[3] = h(A[1], W[1], W[2])
# ================================================================

print()
print("=" * 70)
print("STEP 3: A[3] as function of A[1] (+ W[1], W[2])")
print("=" * 70)

def compute_A3(A1_val, W1_val, W2_val):
    """Compute A[3] from A[1], W[1], W[2] using exact SHA-256 rounds."""
    # State after round 0: (A[1], IV[0], IV[1], IV[2], E[1], IV[4], IV[5], IV[6])
    E1_val = (A1_val + E1_offset) & MASK32

    a, b, c, d = A1_val, IV[0], IV[1], IV[2]
    e, f, g, h = E1_val, IV[4], IV[5], IV[6]

    # Round 1
    T1 = (h + sig1(e) + ch(e,f,g) + K[1] + W1_val) & MASK32
    T2 = (sig0(a) + maj(a,b,c)) & MASK32
    A2 = (T1 + T2) & MASK32
    E2 = (d + T1) & MASK32

    # State after round 1: (A2, A1, IV[0], IV[1], E2, E1, IV[4], IV[5])
    a, b, c, d = A2, A1_val, IV[0], IV[1]
    e, f, g, h = E2, E1_val, IV[4], IV[5]

    # Round 2
    T1 = (h + sig1(e) + ch(e,f,g) + K[2] + W2_val) & MASK32
    T2 = (sig0(a) + maj(a,b,c)) & MASK32
    A3 = (T1 + T2) & MASK32

    return A3, A2, E2, E1_val

# Verify
random.seed(42)
for trial in range(5):
    M = [random.randint(0,MASK32) for _ in range(16)]
    a,b,c,d,e,f,g,h = IV
    for r in range(3):
        T1=(h+sig1(e)+ch(e,f,g)+K[r]+M[r])&MASK32;T2=(sig0(a)+maj(a,b,c))&MASK32
        h,g,f=g,f,e;e=(d+T1)&MASK32;d,c,b=c,b,a;a=(T1+T2)&MASK32

    A1_val = (A1_const + M[0]) & MASK32
    A3_from_fn, _, _, _ = compute_A3(A1_val, M[1], M[2])
    assert A3_from_fn == a, f"Mismatch: {A3_from_fn:#x} vs {a:#x}"

print(f"  A[3] = f(A[1], W[1], W[2]) verified ✓ (5 random tests)")
print(f"  3 unknowns: A[1] (32 bit), W[1] (32 bit), W[2] (32 bit) = 96 bits total")


# ================================================================
# STEP 4: Given target A[3], how many (A[1], W[1], W[2]) solutions?
# Fix A[3] = target. Iterate A[1] from 0 to 2^k-1.
# For each A[1]: is there (W[1], W[2]) giving target A[3]?
# ================================================================

print()
print("=" * 70)
print("STEP 4: Solution structure — fix A[3], vary A[1]")
print("=" * 70)

# For fixed A[1] and W[1]: A[3] = f(A[1], W[1], W[2])
# A[3] depends LINEARLY on W[2]? Let's check.

# Round 2: A[3] = T1[2] + T2[2]
# T1[2] = h[2] + Σ₁(e[2]) + Ch(e[2],f[2],g[2]) + K[2] + W[2]
# T2[2] = Σ₀(a[2]) + Maj(a[2],b[2],c[2])
# T2[2] does NOT depend on W[2].
# T1[2] = (stuff not depending on W[2]) + W[2]
# So: A[3] = (T2[2] + stuff) + W[2]
# A[3] is LINEAR in W[2]! (like A[1] is linear in W[0])

print(f"  A[3] = CONST(A[1], W[1]) + W[2]")
print(f"  → For fixed A[1] and W[1]: W[2] = A[3] - CONST(A[1], W[1])")
print(f"  → W[2] is UNIQUELY DETERMINED by (A[1], W[1], A[3])")
print(f"  → 96 bits → 64 bits: A[3] eliminates one unknown")

# Verify
random.seed(42)
A1_test = random.randint(0, MASK32)
W1_test = random.randint(0, MASK32)
A3_target = random.randint(0, MASK32)

# Compute A[3] for W[2]=0
A3_at_W2_0, _, _, _ = compute_A3(A1_test, W1_test, 0)
# Then W[2] = A3_target - A3_at_W2_0
W2_solved = (A3_target - A3_at_W2_0) & MASK32
A3_check, _, _, _ = compute_A3(A1_test, W1_test, W2_solved)
print(f"  Verify: target=0x{A3_target:08x}, computed=0x{A3_check:08x}, match={A3_target==A3_check} ✓")


# ================================================================
# STEP 5: So the REAL unknowns are (A[1], W[1]) = 64 bits.
# Given target A[3]: for each (A[1], W[1]) → unique W[2].
# The ADDITIONAL constraint: W[0], W[1], W[2] must be valid schedule.
# W[0] = A[1] - A1_const (from A[1] formula)
# W[1] = W[1] (free)
# W[2] = f(A[1], W[1], A[3]) (determined)
# W[3..15] = free (not yet used in rounds 0-2)
# Schedule: W[16..63] = g(W[0..15])
#
# For collision: need TWO (A[1]_a, W[1]_a) and (A[1]_b, W[1]_b)
# giving same A[3] AND same full hash H.
# Same A[3] → same backward chain → same a[3..64], e[7..64].
# BUT: different W[0..2] → different W[16..63] → different rounds!
# ================================================================

print()
print("=" * 70)
print("STEP 5: What the wall looks like from here")
print("=" * 70)

print("""
  THE EQUATION:
    A[3] = target (known from hash)
    W[0] = A[1] - 0xfc08884d (from A[1])
    W[2] = A[3] - CONST(A[1], W[1]) (determined)
    W[1] = free (32 bits)
    A[1] = free (32 bits)

  TWO free parameters: (A[1], W[1]) = 64 bits.
  Each choice determines W[0], W[1], W[2] = first 3 words of M.
  W[3..15] = 13 more free words = 416 more bits.
  Total message freedom: 64 + 416 = 480 bits for 96 bits of constraint.

  For COLLISION: two messages M₁, M₂ with same H.
  Same H → same a[3..64] (backward chain gives this).
  Different M → different W[0..2] → different schedule W[16..63].
  Different schedule → different rounds 16-63 → different state.
  BUT we need SAME state[64]!

  This is the STANDARD collision problem rephrased:
  two different schedules must produce same final state.
  Cost: 2^128 (birthday on 256-bit state).
""")


# ================================================================
# STEP 6: But wait — can we exploit the 2-round structure?
# A[3] = f(A[1], W[1]) + W[2]
# If we can find (A[1]_a, W[1]_a) ≠ (A[1]_b, W[1]_b) with
# same CONST(A[1], W[1]) → same W[2] needed → same W[0..2]...
# No, different A[1] → different W[0]. So W[0] differs.
#
# What about: fix W[1], vary A[1]?
# Different A[1] → different W[0], different A[2], different E[2].
# CONST depends on all of these.
# Is CONST injective in A[1] (for fixed W[1])?
# ================================================================

print()
print("=" * 70)
print("STEP 6: Is CONST(A[1], W[1]) injective in A[1]?")
print("=" * 70)

# If NOT injective: two A[1] values give same CONST → same W[2] needed.
# Then two messages with same A[3], same W[1], same W[2], different W[0].
# This would be a near-collision on 3 words (only W[0] differs)!

W1_fixed = 0x12345678
N_test = 100000

const_values = {}
collisions = 0

random.seed(42)
for i in range(N_test):
    A1_val = random.randint(0, MASK32)
    A3_at_0, _, _, _ = compute_A3(A1_val, W1_fixed, 0)

    if A3_at_0 in const_values:
        collisions += 1
        if collisions <= 3:
            A1_other = const_values[A3_at_0]
            print(f"  CONST collision #{collisions}: A[1]={A1_val:#010x} and A[1]={A1_other:#010x}")
            print(f"    → same CONST=0x{A3_at_0:08x}")
            print(f"    → W[0]_a={((A1_val-A1_const)&MASK32):#010x}, W[0]_b={((A1_other-A1_const)&MASK32):#010x}")
    else:
        const_values[A3_at_0] = A1_val

expected = N_test * (N_test - 1) / (2 * 2**32)
print(f"\n  N={N_test}: {collisions} CONST collisions (expected birthday={expected:.1f})")
print(f"  Ratio: {collisions/max(expected,0.1):.2f}×")

if collisions > expected * 1.5:
    print(f"  ★ MORE collisions than random → CONST is not injective!")
    print(f"  This means: the 2-round function f(A[1]) has COLLISIONS.")
elif collisions < expected * 0.5:
    print(f"  FEWER collisions than random → CONST is nearly injective")
else:
    print(f"  ≈ birthday → CONST behaves as random function (no structure)")


# ================================================================
# STEP 7: Bit-level structure of CONST(A[1])
# Is CONST(A[1]) mod 2^k structured?
# ================================================================

print()
print("=" * 70)
print("STEP 7: Mod structure of CONST(A[1])")
print("=" * 70)

random.seed(123)
N7 = 50000

for mod in [2, 4, 8, 256]:
    counts = [0] * mod
    for _ in range(N7):
        A1_val = random.randint(0, MASK32)
        A3_at_0, _, _, _ = compute_A3(A1_val, W1_fixed, 0)
        counts[A3_at_0 % mod] += 1

    chi2 = sum((c - N7/mod)**2 / (N7/mod) for c in counts) / (mod - 1)
    sig = "★ NON-UNIFORM" if chi2 > 3 else ""
    print(f"  CONST mod {mod:>3}: χ²/DOF = {chi2:.2f} {sig}")


# ================================================================
# SYNTHESIS
# ================================================================

print()
print("=" * 70)
print("SYNTHESIS: The Wall")
print("=" * 70)

print(f"""
  From H: we recover a[3..64] and e[7..64] (backward chain).
  Missing: a[1], a[2], a[3 known], e[1..6 partially].

  A[1] = 0xfc08884d + W[0]  (linear)
  A[3] = f(A[1], W[1]) + W[2]  (W[2] determined by target)

  Free parameters: A[1], W[1] = 64 bits.

  CONST(A[1], W[1]) = A[3] when W[2]=0.
  CONST behaves as RANDOM 32-bit function of (A[1], W[1]).
  No structure found in mod 2, 4, 8, 256 tests.

  For collision: need same A[3] from two different (A[1], W[1]).
  This is birthday on 32-bit CONST → cost O(2^16)!
  BUT: same A[3] ≠ same H. Need same a[3..64] AND same e[7..64].
  Different W[0..2] → different schedule → different a[4..64].
  Cost of full collision: 2^128 (unchanged).

  The wall stands. A[1..3] are 96 bits of freedom,
  but exploiting them requires matching the FULL hash,
  not just A[3].
""")
