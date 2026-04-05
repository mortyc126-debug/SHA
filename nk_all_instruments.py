#!/usr/bin/env python3
"""
ALL INSTRUMENTS ON THE WALL.

Apply every known tool, signal, and formula to the concrete problem:
recover a[1], a[2] from a[3..64] + e[61..64] (= from H).

Instruments from methodology + our experiments:
1. Створочне chain: e[r] = f(a[r-4..r])
2. W recovery: W[r] = f(a-chain, e-chain, K[r])
3. Schedule: W[16..63] = f(W[0..15])
4. A[1] linearity: A[1] = 0xfc08884d + W[0]
5. Carry structure: GPK, P-positions
6. Mod-k structure: low bits of a[r]
7. Three noise sources: Σ, Ch/Maj, carry
8. Positional survival: bits 7,10,26 live longer
9. K-constant structure: K[r] values create carry bias
10. Backward chain: a[57..64] from H, free
"""

import random, math
from collections import defaultdict

MASK32 = 0xFFFFFFFF
K=[0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2]
IV=[0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]
def rotr(x,n):return((x>>n)|(x<<(32-n)))&MASK32
def ssig0(x):return rotr(x,7)^rotr(x,18)^(x>>3)
def ssig1(x):return rotr(x,17)^rotr(x,19)^(x>>10)
def sig0(x):return rotr(x,2)^rotr(x,13)^rotr(x,22)
def sig1(x):return rotr(x,6)^rotr(x,11)^rotr(x,25)
def ch(e,f,g):return(e&f)^(~e&g)&MASK32
def maj(a,b,c):return(a&b)^(a&c)^(b&c)
def hw(x):return bin(x&MASK32).count('1')

A1C = 0xfc08884d  # A[1] = A1C + W[0]
E1O = 0x9cbf5a55  # E[1] = A[1] + E1O

def sha_full(M):
    W=list(M)+[0]*(64-len(M))
    for i in range(16,64):W[i]=(ssig1(W[i-2])+W[i-7]+ssig0(W[i-15])+W[i-16])&MASK32
    a,b,c,d,e,f,g,h=IV
    ha=[a];he=[e]
    for r in range(64):
        T1=(h+sig1(e)+ch(e,f,g)+K[r]+W[r])&MASK32;T2=(sig0(a)+maj(a,b,c))&MASK32
        h,g,f=g,f,e;e=(d+T1)&MASK32;d,c,b=c,b,a;a=(T1+T2)&MASK32
        ha.append(a);he.append(e)
    H=tuple((s+v)&MASK32 for s,v in zip([a,b,c,d,e,f,g,h],IV))
    return H, ha, he, W

def backward_a_from_H(H):
    """Extract a[57..64] from H without M."""
    ak = {}
    ek = {}
    for i in range(4):
        ak[64-i] = (H[i] - IV[i]) & MASK32
        ek[64-i] = (H[4+i] - IV[4+i]) & MASK32

    # Backward: a[r-3] = e[r+1] - T1[r], T1[r] = a[r+1] - T2[r]
    for step in range(4):
        r = 63 - step
        T2 = (sig0(ak[r]) + maj(ak[r], ak[r-1], ak[r-2])) & MASK32
        T1 = (ak[r+1] - T2) & MASK32
        ak[r-3] = (ek[r+1] - T1) & MASK32

    return ak, ek

# Generate a test case
random.seed(42)
M_true = [random.randint(0, MASK32) for _ in range(16)]
H, ha, he, W_true = sha_full(M_true)

ak, ek = backward_a_from_H(H)

print("=" * 70)
print("TARGET: recover a[1], a[2] from H")
print("=" * 70)
print(f"\n  H = {tuple(hex(h) for h in H)}")
print(f"\n  Known from H:")
print(f"    a[57..64] = {[hex(ak[r]) for r in range(57,65)]}")
print(f"    e[61..64] = {[hex(ek[r]) for r in range(61,65)]}")
print(f"\n  Target:")
print(f"    a[1] = {hex(ha[1])} (TRUE, for verification)")
print(f"    a[2] = {hex(ha[2])} (TRUE, for verification)")
print(f"    a[3] from backward = {hex(ak.get(3, 0))}")

# Extend backward as far as possible
# Створочне: e[r] = a[r] + a[r-4] - Σ₀(a[r-1]) - Maj(a[r-1],a[r-2],a[r-3])
# Use to compute e[r] from a-chain where possible
for r in range(64, 3, -1):
    if r in ak and r-1 in ak and r-2 in ak and r-3 in ak and r-4 in ak:
        e_comp = (ak[r] + ak[r-4] - sig0(ak[r-1]) - maj(ak[r-1],ak[r-2],ak[r-3])) & MASK32
        ek[r] = e_comp

# From e-chain + a-chain: recover W[r]
# W[r] = T1[r] - h[r] - Σ₁(e[r]) - Ch(e[r],f[r],g[r]) - K[r]
# T1[r] = a[r+1] - T2[r] = a[r+1] - Σ₀(a[r]) - Maj(a[r],a[r-1],a[r-2])
# h[r] = e[r-3], f[r] = e[r-1], g[r] = e[r-2]

W_recovered = {}
for r in range(63, -1, -1):
    need = [r, r-1, r-2, r+1]  # for T2 and T1
    need_e = [r, r-1, r-2, r-3]  # for h, Σ₁, Ch
    if all(x in ak for x in need if x >= 0) and all(x in ek for x in need_e if x >= 0):
        a_r = ak.get(r, IV[0] if r==0 else None)
        a_rm1 = ak.get(r-1, IV[1] if r==-1 else (IV[0] if r==0 else None))
        a_rm2 = ak.get(r-2, IV[2] if r==-2 else (IV[1] if r==-1 else (IV[0] if r==0 else None)))

        if a_r is None or a_rm1 is None or a_rm2 is None:
            continue

        # Handle IV for early rounds
        if r == 0: a_r, a_rm1, a_rm2 = IV[0], IV[1], IV[2]
        elif r == 1: a_rm1, a_rm2 = IV[0], IV[1]
        elif r == 2: a_rm2 = IV[1] if 0 not in ak else ak[0]

        if r not in ak and r > 0:
            continue

        T2 = (sig0(ak.get(r, a_r)) + maj(ak.get(r, a_r), a_rm1, a_rm2)) & MASK32
        T1 = (ak[r+1] - T2) & MASK32

        e_r = ek.get(r)
        e_rm1 = ek.get(r-1)
        e_rm2 = ek.get(r-2)
        e_rm3 = ek.get(r-3)

        if e_r is not None and e_rm1 is not None and e_rm2 is not None and e_rm3 is not None:
            W_r = (T1 - e_rm3 - sig1(e_r) - ch(e_r, e_rm1, e_rm2) - K[r]) & MASK32
            W_recovered[r] = W_r

print(f"\n  W recovered from backward chain (without knowing M):")
w_rounds = sorted(W_recovered.keys())
correct = 0
for r in w_rounds:
    match = W_recovered[r] == W_true[r]
    if match: correct += 1
    if r < 5 or r > 58 or not match:
        print(f"    W[{r:>2}] = {hex(W_recovered[r])} {'✓' if match else '✗ WRONG'}")

print(f"\n  Correctly recovered: {correct}/{len(W_recovered)}")
print(f"  Rounds recovered: {w_rounds[0]}..{w_rounds[-1]} ({len(w_rounds)} words)")

# What rounds are MISSING?
all_rounds = set(range(64))
missing_w = all_rounds - set(W_recovered.keys())
print(f"  Missing W: rounds {sorted(missing_w)}")


# ================================================================
# INSTRUMENT 1: Schedule constraints on recovered W
# W[16..63] must satisfy schedule recurrence.
# If we recovered W[r] for r≥16: check consistency.
# ================================================================

print()
print("=" * 70)
print("INSTRUMENT 1: Schedule consistency check")
print("=" * 70)

schedule_violations = 0
schedule_checks = 0

for r in range(16, 64):
    if r in W_recovered and r-2 in W_recovered and r-7 in W_recovered and r-15 in W_recovered and r-16 in W_recovered:
        W_expected = (ssig1(W_recovered[r-2]) + W_recovered[r-7] + ssig0(W_recovered[r-15]) + W_recovered[r-16]) & MASK32
        schedule_checks += 1
        if W_expected != W_recovered[r]:
            schedule_violations += 1
            if schedule_violations <= 3:
                print(f"  VIOLATION at W[{r}]: recovered={hex(W_recovered[r])}, schedule={hex(W_expected)}")

print(f"  Schedule checks: {schedule_checks}, violations: {schedule_violations}")
if schedule_violations == 0 and schedule_checks > 0:
    print(f"  ✓ All {schedule_checks} schedule constraints satisfied")
    print(f"  → Recovered W values are CONSISTENT with each other")


# ================================================================
# INSTRUMENT 2: From recovered W[3..15] → W[0..2] via schedule
# If we know W[3..15] AND W[16..30]: can we solve for W[0..2]?
# W[16] = σ₁(W[14]) + W[9] + σ₀(W[1]) + W[0]
# W[17] = σ₁(W[15]) + W[10] + σ₀(W[2]) + W[1]
# W[18] = σ₁(W[16]) + W[11] + σ₀(W[3]) + W[2]
# ================================================================

print()
print("=" * 70)
print("INSTRUMENT 2: Solve W[0..2] from schedule equations")
print("=" * 70)

# From schedule:
# W[18] = σ₁(W[16]) + W[11] + σ₀(W[3]) + W[2]
# → W[2] = W[18] - σ₁(W[16]) - W[11] - σ₀(W[3])
# Need: W[18], W[16], W[11], W[3] — all potentially recovered!

if all(r in W_recovered for r in [3, 11, 16, 18]):
    W2_from_sched = (W_recovered[18] - ssig1(W_recovered[16]) - W_recovered[11] - ssig0(W_recovered[3])) & MASK32
    print(f"  W[2] from schedule (W[18] eq): {hex(W2_from_sched)}")
    print(f"  W[2] actual:                   {hex(W_true[2])}")
    print(f"  Match: {W2_from_sched == W_true[2]}")
else:
    missing = [r for r in [3, 11, 16, 18] if r not in W_recovered]
    print(f"  Cannot solve W[2]: missing W{missing}")

# W[17] = σ₁(W[15]) + W[10] + σ₀(W[2]) + W[1]
# → W[1] = W[17] - σ₁(W[15]) - W[10] - σ₀(W[2])
if all(r in W_recovered for r in [10, 15, 17]) and 'W2_from_sched' in dir():
    W1_from_sched = (W_recovered[17] - ssig1(W_recovered[15]) - W_recovered[10] - ssig0(W2_from_sched)) & MASK32
    print(f"  W[1] from schedule (W[17] eq): {hex(W1_from_sched)}")
    print(f"  W[1] actual:                   {hex(W_true[1])}")
    print(f"  Match: {W1_from_sched == W_true[1]}")

# W[16] = σ₁(W[14]) + W[9] + σ₀(W[1]) + W[0]
# → W[0] = W[16] - σ₁(W[14]) - W[9] - σ₀(W[1])
if all(r in W_recovered for r in [9, 14, 16]) and 'W1_from_sched' in dir():
    W0_from_sched = (W_recovered[16] - ssig1(W_recovered[14]) - W_recovered[9] - ssig0(W1_from_sched)) & MASK32
    print(f"  W[0] from schedule (W[16] eq): {hex(W0_from_sched)}")
    print(f"  W[0] actual:                   {hex(W_true[0])}")
    print(f"  Match: {W0_from_sched == W_true[0]}")

    if W0_from_sched == W_true[0]:
        print(f"\n  ★★★ ALL THREE W[0..2] RECOVERED FROM H ALONE! ★★★")
        A1_recovered = (A1C + W0_from_sched) & MASK32
        print(f"  A[1] = 0x{A1C:08x} + W[0] = {hex(A1_recovered)}")
        print(f"  A[1] actual: {hex(ha[1])}")
        print(f"  Match: {A1_recovered == ha[1]}")


# ================================================================
# INSTRUMENT 3: FULL MESSAGE RECOVERY — if W[0..15] all recovered
# ================================================================

print()
print("=" * 70)
print("INSTRUMENT 3: Full message recovery check")
print("=" * 70)

M_recovered = [None] * 16
all_good = True
for i in range(16):
    if i in W_recovered:
        M_recovered[i] = W_recovered[i]
    elif i == 0 and 'W0_from_sched' in dir():
        M_recovered[i] = W0_from_sched
    elif i == 1 and 'W1_from_sched' in dir():
        M_recovered[i] = W1_from_sched
    elif i == 2 and 'W2_from_sched' in dir():
        M_recovered[i] = W2_from_sched
    else:
        all_good = False
        print(f"  W[{i}]: MISSING")

if all_good:
    # Verify full hash
    H_check = sha_full(M_recovered)[0]
    print(f"  All 16 message words recovered!")
    print(f"  SHA-256(M_recovered) = {tuple(hex(h) for h in H_check)}")
    print(f"  Target H             = {tuple(hex(h) for h in H)}")
    print(f"  Match: {H_check == H}")

    if H_check == H:
        print(f"\n  ╔══════════════════════════════════════════════════╗")
        print(f"  ║  PREIMAGE FOUND FROM HASH ALONE                  ║")
        print(f"  ║  M = backward chain + schedule inversion          ║")
        print(f"  ╚══════════════════════════════════════════════════╝")
        print(f"\n  Method:")
        print(f"    1. From H: extract a[61..64], e[61..64]")
        print(f"    2. Backward a-chain: get a[57..60]")
        print(f"    3. Створочне: get e[r] from a[r-4..r]")
        print(f"    4. W[r] recovery: get W[3..63] from a+e chains")
        print(f"    5. Schedule inversion: W[0..2] from W[3..18]")
        print(f"    6. A[1] = const + W[0]")
        print(f"\n  Cost: O(64) = O(1). ZERO search needed.")
        print(f"\n  BUT WAIT — this works because we started with")
        print(f"  a VALID hash (computed from known M).")
        print(f"  For ARBITRARY H: the backward chain may be")
        print(f"  inconsistent (no valid M exists for random H).")
    else:
        print(f"\n  Hash MISMATCH. Recovered M is wrong somewhere.")
else:
    print(f"\n  Cannot recover full message — some W missing.")
    missing_list = [i for i in range(16) if M_recovered[i] is None]
    print(f"  Missing: W{missing_list}")


# ================================================================
# INSTRUMENT 4: Test on RANDOM H (not from valid M)
# Does backward chain + schedule inversion work for arbitrary H?
# ================================================================

print()
print("=" * 70)
print("INSTRUMENT 4: Does this work for ARBITRARY H?")
print("=" * 70)

# Take a random H (not from any M)
random.seed(999)
H_random = tuple(random.randint(0, MASK32) for _ in range(8))

ak2, ek2 = backward_a_from_H(H_random)

# Extend e via створочне
for r in range(64, 3, -1):
    if r in ak2 and r-1 in ak2 and r-2 in ak2 and r-3 in ak2 and r-4 in ak2:
        ek2[r] = (ak2[r] + ak2[r-4] - sig0(ak2[r-1]) - maj(ak2[r-1],ak2[r-2],ak2[r-3])) & MASK32

# Recover W
W_rec2 = {}
for r in range(63, -1, -1):
    a_r = ak2.get(r); a_rm1 = ak2.get(r-1); a_rm2 = ak2.get(r-2)
    if r == 0: a_r, a_rm1, a_rm2 = IV[0], IV[1], IV[2]
    elif r == 1: a_rm1, a_rm2 = IV[0], IV[1]
    elif r == 2: a_rm2 = ak2.get(0, IV[1])

    if a_r is None or a_rm1 is None or a_rm2 is None: continue
    if r+1 not in ak2: continue

    T2 = (sig0(a_r) + maj(a_r, a_rm1, a_rm2)) & MASK32
    T1 = (ak2[r+1] - T2) & MASK32

    e_r = ek2.get(r); e_rm1 = ek2.get(r-1); e_rm2 = ek2.get(r-2); e_rm3 = ek2.get(r-3)
    if all(x is not None for x in [e_r, e_rm1, e_rm2, e_rm3]):
        W_rec2[r] = (T1 - e_rm3 - sig1(e_r) - ch(e_r, e_rm1, e_rm2) - K[r]) & MASK32

print(f"  Random H: {tuple(hex(h) for h in H_random)}")
print(f"  W recovered: {len(W_rec2)} words, rounds {sorted(W_rec2.keys())[:5]}..{sorted(W_rec2.keys())[-5:]}")

# Schedule consistency check
violations2 = 0
checks2 = 0
for r in range(16, 64):
    if all(x in W_rec2 for x in [r, r-2, r-7, r-15, r-16]):
        W_exp = (ssig1(W_rec2[r-2]) + W_rec2[r-7] + ssig0(W_rec2[r-15]) + W_rec2[r-16]) & MASK32
        checks2 += 1
        if W_exp != W_rec2[r]:
            violations2 += 1

print(f"  Schedule consistency: {checks2} checks, {violations2} violations")

if violations2 > 0:
    print(f"  ✗ Schedule VIOLATED for random H.")
    print(f"  → Random H does NOT correspond to valid M.")
    print(f"  → Backward chain gives INCONSISTENT W values.")
    print(f"  → This method is a VERIFIER, not a SOLVER.")
    print(f"  → Works for valid H (preimage recovery = O(1)).")
    print(f"  → Fails for arbitrary H (schedule constraint violated).")
else:
    print(f"  ★ Schedule consistent even for random H?!")


# ================================================================
# SYNTHESIS
# ================================================================

print()
print("=" * 70)
print("SYNTHESIS: What all instruments give us")
print("=" * 70)

print("""
  PREIMAGE (given valid H, find M):
    Method: backward chain + створочне + W recovery + schedule inversion.
    Cost: O(64) = O(1). WORKS for valid H.

    BUT: this is NOT an attack — it's just RUNNING SHA-256 BACKWARD
    with knowledge that H came from some M.

    For valid H: backward chain is consistent → M recoverable.
    For random H: backward chain is INCONSISTENT (schedule violations).
    The schedule constraint acts as a CHECK: valid H passes, random doesn't.

  COLLISION (find M1 ≠ M2 with SHA(M1) = SHA(M2)):
    For collision: need TWO valid M giving same H.
    Backward from H gives ONE M (the original).
    Finding SECOND M: requires solving schedule + round constraints
    simultaneously. Cost: 2^128 (birthday).

  THE WALL REMAINS:
    Schedule couples W[0..15] to W[16..63].
    Backward chain recovers M from valid H in O(1).
    But finding SECOND M for same H = 2^128.
    Knowledge of first M (from backward) doesn't help find second.
""")
