#!/usr/bin/env python3
"""
Flow Mathematics v0.4 — Theorem 6' verification.

Prediction P2: SHA-256 with shift register length L=3
(e_new = c + T1 instead of d + T1) has κ=32.

Test: modify SHA-256 so that e_new = c + T1 (L=3 functional),
keeping all 8 registers and 256-bit output.
Then check: does H[4] = f(H[0], H[1], H[2], H[3])?

Also test L=5 (e_new = state_before_d + T1) as control.
"""

import random
import time

MASK32 = 0xFFFFFFFF

K = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
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
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
]

IV = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
]

def rotr(x, n):
    return ((x >> n) | (x << (32 - n))) & MASK32

def sig0(x):
    return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)

def sig1(x):
    return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)

def ssig0(x):
    return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)

def ssig1(x):
    return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)

def ch(e, f, g):
    return (e & f) ^ (~e & g) & MASK32

def maj(a, b, c):
    return (a & b) ^ (a & c) ^ (b & c)

def sha256_compress_L(msg_words, L=4, rounds=64):
    """
    SHA-256 compression with variable effective shift length L.

    L=4 (standard): e_new = d + T1  (d = a[r-3], 4th in shift chain)
    L=3: e_new = c + T1  (c = a[r-2], 3rd in shift chain)
    L=5: e_new uses a value from BEFORE the current shift window
         We implement this by storing a[r-4] explicitly.

    All versions keep 8 registers and 256-bit output.
    """
    assert len(msg_words) == 16

    W = list(msg_words)
    for i in range(16, max(rounds, 64)):
        W.append((ssig1(W[i-2]) + W[i-7] + ssig0(W[i-15]) + W[i-16]) & MASK32)

    a, b, c, d, e, f, g, h = IV

    # For L=5, we need to track a[r-4] (one step before d)
    a_history = [a]  # will accumulate a-values

    for r in range(rounds):
        T1 = (h + sig1(e) + ch(e, f, g) + K[r] + W[r]) & MASK32
        T2 = (sig0(a) + maj(a, b, c)) & MASK32

        h = g
        g = f
        f = e

        # This is where L matters:
        if L == 3:
            e = (c + T1) & MASK32      # L=3: use c = a[r-2]
        elif L == 4:
            e = (d + T1) & MASK32      # L=4: standard SHA-256
        elif L == 5:
            # L=5: use a[r-4] (one beyond current d = a[r-3])
            if len(a_history) >= 4:
                a_rm4 = a_history[-4]   # a from 4 rounds ago
            else:
                a_rm4 = IV[3]           # fallback to IV for first rounds
            e = (a_rm4 + T1) & MASK32

        d = c
        c = b
        b = a
        a = (T1 + T2) & MASK32

        a_history.append(a)

    H = [
        (a + IV[0]) & MASK32,
        (b + IV[1]) & MASK32,
        (c + IV[2]) & MASK32,
        (d + IV[3]) & MASK32,
        (e + IV[4]) & MASK32,
        (f + IV[5]) & MASK32,
        (g + IV[6]) & MASK32,
        (h + IV[7]) & MASK32,
    ]
    return H

def random_msg():
    return [random.randint(0, MASK32) for _ in range(16)]


# ============================================================
# TEST A: Verify structural identity for each L
#
# For L=4 (standard): e_new = d + T1 = a[R-4+1] + T1
#   Identity: H[0] - H[4] = T2(H[1],H[2],H[3]) - d_entering + const
#   d_entering = a[R-4] — NOT in output (shifted out)
#   → Identity contains unobservable → κ=0
#
# For L=3: e_new = c + T1 = a[R-3+1] + T1
#   Identity: H[0] - H[4] = T2(H[1],H[2],H[3]) - c_entering + const
#   c_entering at last round = a[R-3] = H[3] - IV[3] — IN output!
#   → Identity is observable → κ=32?
#
# For L=5: e_new = a[R-5+1] + T1
#   a[R-4] NOT in output → κ=0 (even more hidden)
# ============================================================

print("=" * 70)
print("TEST A: Structural identity for different shift lengths L")
print("=" * 70)

N = 10000

for L in [3, 4, 5]:
    violations = 0

    for _ in range(N):
        M = random_msg()
        H = sha256_compress_L(M, L=L, rounds=64)

        # Reconstruct a-chain from H
        a_R   = (H[0] - IV[0]) & MASK32   # a[R]
        a_Rm1 = (H[1] - IV[1]) & MASK32   # a[R-1]
        a_Rm2 = (H[2] - IV[2]) & MASK32   # a[R-2]
        a_Rm3 = (H[3] - IV[3]) & MASK32   # a[R-3]
        e_R   = (H[4] - IV[4]) & MASK32   # e[R]

        T2_last = (sig0(a_Rm1) + maj(a_Rm1, a_Rm2, a_Rm3)) & MASK32
        T1_last = (a_R - T2_last) & MASK32

        if L == 3:
            # e_new = c + T1, c at last round = a[R-2] = a_Rm2
            # But wait: c entering round R-1 = b entering round R-2 = a entering round R-3
            # After shift: at round R-1, c = a[R-3]
            # Hmm, need to be careful about indexing.
            #
            # At round R-1 (the last round, producing state R):
            #   The state ENTERING round R-1 is:
            #   a_in = a[R-1], b_in = a[R-2], c_in = a[R-3], d_in = a[R-4]
            #   e_in = e[R-1], f_in = e[R-2], g_in = e[R-3], h_in = e[R-4]
            #
            # With L=3: e_new = c_in + T1 = a[R-3] + T1
            # a[R-3] = H[3] - IV[3] — observable!

            d_entering = a_Rm3  # For L=3: c_in = a[R-3] = observable
            e_predicted = (d_entering + T1_last) & MASK32

        elif L == 4:
            # Standard: e_new = d_in + T1 = a[R-4] + T1
            # a[R-4] NOT in output
            # We can't predict e_R without knowing a[R-4]
            # Just check: does e_R = a_Rm3 + T1? (this uses wrong value)
            d_entering = a_Rm3  # WRONG: should be a[R-4], not a[R-3]
            e_predicted = (d_entering + T1_last) & MASK32

        elif L == 5:
            # e_new = a[R-5+1] + T1 = a[R-4] + T1 ... wait
            # For L=5, we use a_history[-4] which is a[R-4]
            # a[R-4] NOT in output
            d_entering = a_Rm3  # WRONG: should be a[R-5], even further out
            e_predicted = (d_entering + T1_last) & MASK32

        if e_predicted != e_R:
            violations += 1

    status = "✓ IDENTITY HOLDS (κ=32)" if violations == 0 else f"✗ {violations}/{N} violations (κ=0)"
    print(f"  L={L}: {status}")


# ============================================================
# TEST B: For L=3, verify H[4] = f(H[0], H[1], H[2], H[3])
# ============================================================

print()
print("=" * 70)
print("TEST B: For L=3, verify H[4] = f(H[0], H[1], H[2], H[3])")
print("=" * 70)

violations = 0
N = 10000

for _ in range(N):
    M = random_msg()
    H = sha256_compress_L(M, L=3, rounds=64)

    a_Rm1 = (H[1] - IV[1]) & MASK32
    a_Rm2 = (H[2] - IV[2]) & MASK32
    a_Rm3 = (H[3] - IV[3]) & MASK32

    T2_last = (sig0(a_Rm1) + maj(a_Rm1, a_Rm2, a_Rm3)) & MASK32
    T1_last = ((H[0] - IV[0]) - T2_last) & MASK32

    # For L=3: e_new = a[R-3] + T1 = (H[3]-IV[3]) + T1
    H4_predicted = (a_Rm3 + T1_last + IV[4]) & MASK32

    if H4_predicted != H[4]:
        violations += 1

print(f"  H[4] = f(H[0..3]): {violations}/{N} violations")
if violations == 0:
    print(f"  ✓ H[4] FULLY DETERMINED by H[0..3] for L=3")
    print(f"  → κ = 32 bits, effective birthday = 2^{{(256-32)/2}} = 2^{{112}}")
else:
    print(f"  ✗ H[4] NOT determined by H[0..3]")


# ============================================================
# TEST C: For L=4 (standard), verify identity FAILS
# ============================================================

print()
print("=" * 70)
print("TEST C: For L=4 (standard SHA-256), verify identity FAILS")
print("=" * 70)

violations = 0
N = 10000

for _ in range(N):
    M = random_msg()
    H = sha256_compress_L(M, L=4, rounds=64)

    a_Rm1 = (H[1] - IV[1]) & MASK32
    a_Rm2 = (H[2] - IV[2]) & MASK32
    a_Rm3 = (H[3] - IV[3]) & MASK32

    T2_last = (sig0(a_Rm1) + maj(a_Rm1, a_Rm2, a_Rm3)) & MASK32
    T1_last = ((H[0] - IV[0]) - T2_last) & MASK32

    # For L=4: e_new = a[R-4] + T1, but a[R-4] ≠ a[R-3]
    # Using a[R-3] (wrong value for L=4):
    H4_predicted = (a_Rm3 + T1_last + IV[4]) & MASK32

    if H4_predicted != H[4]:
        violations += 1

print(f"  H[4] = f(H[0..3]) using a[R-3]: {violations}/{N} violations")
if violations == N:
    print(f"  ✓ CORRECTLY FAILS: a[R-4] ≠ a[R-3], identity broken for L=4")
    print(f"  → κ = 0, birthday = 2^128 (standard)")
elif violations == 0:
    print(f"  ✗ UNEXPECTED: identity holds for L=4 (should fail!)")
else:
    print(f"  ~ Partial: {violations}/{N}")


# ============================================================
# TEST D: Verify the CORRECT structural identity for L=4
# e_new = d_entering + T1, where d_entering = a[R-4]
# We compute a[R-4] from the ACTUAL computation
# ============================================================

print()
print("=" * 70)
print("TEST D: Correct structural identity for L=4 (with a[R-4])")
print("=" * 70)

def sha256_compress_L4_with_history(msg_words, rounds=64):
    """Standard SHA-256 but return a-chain history."""
    W = list(msg_words)
    for i in range(16, max(rounds, 64)):
        W.append((ssig1(W[i-2]) + W[i-7] + ssig0(W[i-15]) + W[i-16]) & MASK32)

    a, b, c, d, e, f, g, h = IV
    a_hist = [a]

    for r in range(rounds):
        T1 = (h + sig1(e) + ch(e, f, g) + K[r] + W[r]) & MASK32
        T2 = (sig0(a) + maj(a, b, c)) & MASK32
        h, g, f = g, f, e
        e = (d + T1) & MASK32
        d, c, b = c, b, a
        a = (T1 + T2) & MASK32
        a_hist.append(a)

    H = [(s + iv) & MASK32 for s, iv in zip([a,b,c,d,e,f,g,h], IV)]
    return H, a_hist

violations = 0
N = 10000

for _ in range(N):
    M = random_msg()
    H, a_hist = sha256_compress_L4_with_history(M, rounds=64)

    a_Rm1 = (H[1] - IV[1]) & MASK32
    a_Rm2 = (H[2] - IV[2]) & MASK32
    a_Rm3 = (H[3] - IV[3]) & MASK32
    a_Rm4 = a_hist[-4]  # a[R-4] = a[60] — the HIDDEN value

    T2_last = (sig0(a_Rm1) + maj(a_Rm1, a_Rm2, a_Rm3)) & MASK32
    T1_last = ((H[0] - IV[0]) - T2_last) & MASK32

    # Correct: e_new = a[R-4] + T1
    H4_correct = (a_Rm4 + T1_last + IV[4]) & MASK32

    if H4_correct != H[4]:
        violations += 1

print(f"  H[4] = a[R-4] + T1 + IV[4]: {violations}/{N} violations")
if violations == 0:
    print(f"  ✓ STRUCTURAL IDENTITY CONFIRMED with hidden a[R-4]")
    print(f"  a[R-4] is the 'invisible link' — present in computation but absent from output")
else:
    print(f"  ✗ ERROR in structural identity")


# ============================================================
# SUMMARY
# ============================================================

print()
print("=" * 70)
print("SUMMARY: Flow Mathematics Theorem 6'")
print("=" * 70)
print()
print("  SHA-256 shift register length L and κ:")
print()
print("  ┌─────┬───────────────────────┬─────┬──────────┐")
print("  │  L  │  d_entering visible?  │  κ  │ Birthday │")
print("  ├─────┼───────────────────────┼─────┼──────────┤")
print("  │  3  │  YES (= a[R-3] = H[3])│  32 │  2^112   │")
print("  │  4  │  NO  (= a[R-4] ∉ H)  │   0 │  2^128   │")
print("  │  5  │  NO  (= a[R-5] ∉ H)  │   0 │  2^128   │")
print("  └─────┴───────────────────────┴─────┴──────────┘")
print()
print("  SHA-256 uses L=4 = minimum L for κ=0.")
print("  This is a DESIGN CHOICE, not coincidence.")
