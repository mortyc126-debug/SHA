#!/usr/bin/env python3
"""
Carry-Conditional Solver (CCS) v0.1

The algorithm no one tested:
  Step A: Fix carry profile c. Solve SHA-256|c (quadratic, no carry nonlinearity).
  Step B: Compute real carry c'(M). Update c → c'.
  Repeat until convergence or budget exhausted.

Key insight from 1300+ experiments:
  - Fix carry → SHA-256 becomes degree-2 (Ch + Maj only)
  - Carry match for Wang pairs = 60.5%
  - Self-healing τ = 1.8 rounds
  - P(carry=1) = 99% per round (T_CARRY_DOMINANT)

This is NOT Hensel (single lift, fails at mod 4).
This is NOT Newton GF(2) (doesn't converge).
This is NOT SA (no metric optimization).
This is ITERATIVE ALTERNATION between quadratic solve and carry update.
"""

import random, time, math

MASK32 = 0xFFFFFFFF
K=[0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2]
IV=[0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&MASK32
def ssig0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def ssig1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def sig0(x): return rotr(x,2)^rotr(x,13)^rotr(x,22)
def sig1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def ch(e,f,g): return (e&f)^(~e&g)&MASK32
def maj(a,b,c): return (a&b)^(a&c)^(b&c)
def hw(x): return bin(x&MASK32).count('1')

def add_with_carry(a, b):
    """Return (sum mod 2^32, carry vector)."""
    s = (a + b) & MASK32
    c = 0; cv = 0
    for k in range(32):
        ak=(a>>k)&1; bk=(b>>k)&1
        c=(ak&bk)|(ak&c)|(bk&c)
        cv|=(c<<k)
    return s, cv

def add_fixed_carry(a, b, carry_vec):
    """Addition where carry is FIXED (not computed from a,b).
    Result: a XOR b XOR (carry_vec << 1), truncated to 32 bits.
    This is what addition looks like when carry is predetermined."""
    # Real addition: sum[k] = a[k] XOR b[k] XOR carry[k-1]
    # With fixed carry: same formula but carry doesn't propagate from inputs
    result = 0
    for k in range(32):
        ak = (a >> k) & 1
        bk = (b >> k) & 1
        ck_prev = (carry_vec >> (k-1)) & 1 if k > 0 else 0
        result |= ((ak ^ bk ^ ck_prev) << k)
    return result & MASK32

def sha256_real(M):
    """Standard SHA-256, returns H and all carry vectors."""
    W=list(M)+[0]*(64-len(M))
    for i in range(16,64): W[i]=(ssig1(W[i-2])+W[i-7]+ssig0(W[i-15])+W[i-16])&MASK32
    a,b,c,d,e,f,g,h=IV
    all_carries = []
    for r in range(64):
        # T1 = h + sig1(e) + ch(e,f,g) + K[r] + W[r]  (4 additions)
        s1, c1 = add_with_carry(h, sig1(e))
        s2, c2 = add_with_carry(s1, ch(e,f,g))
        s3, c3 = add_with_carry(s2, K[r])
        T1, c4 = add_with_carry(s3, W[r])
        # T2 = sig0(a) + maj(a,b,c)  (1 addition)
        T2, c5 = add_with_carry(sig0(a), maj(a,b,c))
        # e_new = d + T1  (1 addition)
        e_new, c6 = add_with_carry(d, T1)
        # a_new = T1 + T2  (1 addition)
        a_new, c7 = add_with_carry(T1, T2)

        all_carries.append((c1,c2,c3,c4,c5,c6,c7))
        h,g,f=g,f,e; e=e_new; d,c,b=c,b,a; a=a_new

    H = tuple((s+v)&MASK32 for s,v in zip([a,b,c,d,e,f,g,h],IV))
    return H, all_carries, W

def sha256_fixed_carry(M, carry_profile):
    """SHA-256 with FIXED carry vectors (not computed from values).
    carry_profile[r] = (c1..c7) for round r."""
    W=list(M)+[0]*(64-len(M))
    # Schedule also has carries — for simplicity use real schedule
    for i in range(16,64): W[i]=(ssig1(W[i-2])+W[i-7]+ssig0(W[i-15])+W[i-16])&MASK32

    a,b,c,d,e,f,g,h=IV
    for r in range(64):
        c1,c2,c3,c4,c5,c6,c7 = carry_profile[r]

        s1 = add_fixed_carry(h, sig1(e), c1)
        s2 = add_fixed_carry(s1, ch(e,f,g), c2)
        s3 = add_fixed_carry(s2, K[r], c3)
        T1 = add_fixed_carry(s3, W[r], c4)
        T2 = add_fixed_carry(sig0(a), maj(a,b,c), c5)
        e_new = add_fixed_carry(d, T1, c6)
        a_new = add_fixed_carry(T1, T2, c7)

        h,g,f=g,f,e; e=e_new; d,c,b=c,b,a; a=a_new

    H = tuple((s+v)&MASK32 for s,v in zip([a,b,c,d,e,f,g,h],IV))
    return H, W

def carry_distance(c1, c2):
    """Number of carry bits that differ between two profiles."""
    total = 0
    for r in range(len(c1)):
        for op in range(7):
            total += hw(c1[r][op] ^ c2[r][op])
    return total

def carry_total_bits(c):
    return len(c) * 7 * 32


# ================================================================
# EXPERIMENT 1: Verify that fixed-carry SHA gives same result
# when carry is correct
# ================================================================

print("=" * 70)
print("EXP 1: Verify fixed-carry = real SHA when carry is correct")
print("=" * 70)

random.seed(42)
N_verify = 1000
ok = 0
for _ in range(N_verify):
    M = [random.randint(0,MASK32) for _ in range(16)]
    H_real, carries, _ = sha256_real(M)
    H_fixed, _ = sha256_fixed_carry(M, carries)
    if H_real == H_fixed:
        ok += 1

print(f"  {ok}/{N_verify} match → {'✓ VERIFIED' if ok==N_verify else '✗ BUG'}")


# ================================================================
# EXPERIMENT 2: CCS Algorithm — iterative carry-conditional solve
#
# Goal: find M2 such that H(M2) = H(M1) (collision with known M1).
#
# Algorithm:
#   1. Start with M1, compute carry_profile(M1)
#   2. Perturb M1 → M2 (small change)
#   3. Compute H2_approx = SHA_fixed_carry(M2, carry_profile(M1))
#      (use M1's carry, not M2's)
#   4. Compute H2_real = SHA_real(M2)
#   5. Measure carry_distance(carry(M1), carry(M2))
#   6. Update: use carry(M2) as new profile
#   7. Adjust M2 to reduce |H2_real - H_target|
# ================================================================

print()
print("=" * 70)
print("EXP 2: CCS v0.1 — Iterative carry-conditional solver")
print("=" * 70)

def ccs_solve(M_target, max_iter=100, verbose=True):
    """
    Try to find M2 with H(M2) = H(M_target).

    Strategy:
    1. Start from M_current = random M
    2. Compute carry(M_current)
    3. Use fixed-carry to estimate which W[i] changes would help
    4. Apply best change
    5. Update carry profile
    6. Repeat
    """
    H_target, carries_target, _ = sha256_real(M_target)

    # Start from random M
    M_current = [random.randint(0, MASK32) for _ in range(16)]
    H_current, carries_current, _ = sha256_real(M_current)

    hw_best = sum(hw(H_current[i] ^ H_target[i]) for i in range(8))
    history = [hw_best]

    for iteration in range(max_iter):
        # Try flipping each bit of each W[i] and pick the best
        # using FIXED-CARRY approximation (cheap: no carry recomputation)
        best_flip = None
        best_hw_approx = hw_best

        # Sample random flips (full search = 512 bits, too slow per iter)
        for _ in range(64):
            w = random.randint(0, 15)
            b = random.randint(0, 31)
            M_try = list(M_current)
            M_try[w] ^= (1 << b)

            # Approximate: use current carry profile
            H_approx, _ = sha256_fixed_carry(M_try, carries_current)
            hw_approx = sum(hw(H_approx[i] ^ H_target[i]) for i in range(8))

            if hw_approx < best_hw_approx:
                best_hw_approx = hw_approx
                best_flip = (w, b)

        if best_flip is not None:
            w, b = best_flip
            M_current[w] ^= (1 << b)

            # NOW: recompute real hash and real carry
            H_current, carries_current, _ = sha256_real(M_current)
            hw_real = sum(hw(H_current[i] ^ H_target[i]) for i in range(8))

            # Measure: how good was the approximation?
            carry_dist = carry_distance(carries_current, carries_target)

            history.append(hw_real)

            if verbose and (iteration < 5 or iteration % 20 == 0 or hw_real < hw_best):
                print(f"    iter {iteration:>3}: HW(δH)={hw_real:>3} (approx={best_hw_approx:>3}) carry_dist={carry_dist:>5}/{carry_total_bits(carries_current)}")

            if hw_real < hw_best:
                hw_best = hw_real

            if hw_real == 0:
                print(f"    ★★★ COLLISION FOUND at iteration {iteration}!")
                return M_current, H_current, history
        else:
            # No improving flip found
            history.append(hw_best)

    return None, None, history


# Run CCS multiple times
print(f"\n  Running CCS v0.1 (10 trials, 200 iterations each)")

random.seed(42)
all_best = []

for trial in range(10):
    M_target = [random.randint(0, MASK32) for _ in range(16)]

    print(f"\n  Trial {trial}:")
    result, H, history = ccs_solve(M_target, max_iter=200, verbose=True)

    best = min(history)
    all_best.append(best)
    print(f"    Best HW(δH) = {best}")

avg_best = sum(all_best) / len(all_best)
overall_best = min(all_best)

print(f"\n  Summary:")
print(f"    Avg best HW(δH) = {avg_best:.1f}")
print(f"    Overall best = {overall_best}")
print(f"    Random baseline = 128")


# ================================================================
# EXPERIMENT 3: Compare CCS with pure random search
# Same number of SHA evaluations → which gets lower HW?
# ================================================================

print()
print("=" * 70)
print("EXP 3: CCS vs random search (same budget)")
print("=" * 70)

N_budget = 200 * 64  # CCS: 200 iters × 64 samples = 12800 SHA evals

random.seed(123)
M_target = [random.randint(0, MASK32) for _ in range(16)]
H_target = sha256_real(M_target)[0]

# CCS
random.seed(123)
_, _, ccs_history = ccs_solve(M_target, max_iter=200, verbose=False)
ccs_best = min(ccs_history)

# Random search: same budget
random.seed(456)
rand_best = 256
for _ in range(N_budget):
    M = [random.randint(0, MASK32) for _ in range(16)]
    H = sha256_real(M)[0]
    hw_d = sum(hw(H[i] ^ H_target[i]) for i in range(8))
    if hw_d < rand_best:
        rand_best = hw_d

# HC search: same budget
random.seed(789)
M_hc = [random.randint(0, MASK32) for _ in range(16)]
H_hc = sha256_real(M_hc)[0]
hc_best = sum(hw(H_hc[i] ^ H_target[i]) for i in range(8))

for _ in range(N_budget):
    w = random.randint(0,15); b = random.randint(0,31)
    M_try = list(M_hc); M_try[w] ^= (1<<b)
    H_try = sha256_real(M_try)[0]
    hw_try = sum(hw(H_try[i] ^ H_target[i]) for i in range(8))
    if hw_try < hc_best:
        M_hc = M_try; hc_best = hw_try

print(f"\n  Budget: {N_budget} SHA-256 evaluations each")
print(f"  Target: H(M_target)")
print(f"\n  {'Method':>15} | {'best HW(δH)':>11}")
print(f"  {'CCS v0.1':>15} | {ccs_best:>11}")
print(f"  {'Random search':>15} | {rand_best:>11}")
print(f"  {'Hill climbing':>15} | {hc_best:>11}")
print(f"  {'Birthday est.':>15} | {'~128':>11}")

if ccs_best < rand_best and ccs_best < hc_best:
    print(f"\n  ★ CCS BEATS both random and HC!")
    print(f"    Improvement over random: {rand_best - ccs_best} bits")
    print(f"    Improvement over HC: {hc_best - ccs_best} bits")
elif ccs_best < rand_best:
    print(f"\n  CCS beats random (+{rand_best-ccs_best}) but not HC ({hc_best-ccs_best})")
else:
    print(f"\n  CCS does NOT beat alternatives.")


# ================================================================
# EXPERIMENT 4: CCS convergence analysis
# Does carry profile CONVERGE? I.e., does carry_distance decrease?
# ================================================================

print()
print("=" * 70)
print("EXP 4: Does carry profile converge?")
print("=" * 70)

random.seed(42)
M1 = [random.randint(0, MASK32) for _ in range(16)]
H1, c1, _ = sha256_real(M1)

# Start with M2 = M1 + small perturbation
M2 = list(M1)
M2[0] ^= 0x80000000  # flip MSB of W[0]

carry_dists = []
hw_dists = []

for iteration in range(50):
    H2, c2, _ = sha256_real(M2)
    cdist = carry_distance(c1, c2)
    hdist = sum(hw(H1[i]^H2[i]) for i in range(8))
    carry_dists.append(cdist)
    hw_dists.append(hdist)

    # CCS step: use c1 (target carry) to guide M2 adjustment
    # Try 32 random flips, pick best under fixed-carry approximation
    best_w, best_b, best_score = -1, -1, hdist
    for _ in range(32):
        w = random.randint(0,15); b = random.randint(0,31)
        M_try = list(M2); M_try[w] ^= (1<<b)
        H_approx, _ = sha256_fixed_carry(M_try, c1)  # USE TARGET CARRY
        score = sum(hw(H_approx[i]^H1[i]) for i in range(8))
        if score < best_score:
            best_w, best_b, best_score = w, b, score

    if best_w >= 0:
        M2[best_w] ^= (1 << best_b)

print(f"  Iter | carry_dist | HW(δH)")
for i in range(min(50, len(carry_dists))):
    if i < 5 or i % 10 == 0 or i == len(carry_dists)-1:
        print(f"  {i:>4} | {carry_dists[i]:>10} | {hw_dists[i]:>6}")

if len(carry_dists) >= 2:
    cd_start = carry_dists[0]
    cd_end = carry_dists[-1]
    hd_start = hw_dists[0]
    hd_end = hw_dists[-1]

    print(f"\n  Carry distance: {cd_start} → {cd_end} ({'↓ converging' if cd_end < cd_start else '→ NOT converging'})")
    print(f"  Hash distance:  {hd_start} → {hd_end} ({'↓ converging' if hd_end < hd_start else '→ NOT converging'})")


# ================================================================
# SYNTHESIS
# ================================================================

print()
print("=" * 70)
print("SYNTHESIS: CCS v0.1")
print("=" * 70)

print(f"""
  CCS (Carry-Conditional Solver) tested:
  - Step A: fix carry, evaluate SHA under fixed carry (approximation)
  - Step B: apply best flip, recompute real carry
  - Iterate

  Results:
  - CCS best HW(δH) = {overall_best} (10 trials × 200 iterations)
  - Random search: {rand_best}
  - Hill climbing: {hc_best}
  - Birthday baseline: 128

  Key question: does fixed-carry approximation HELP guide the search?
  If CCS ≈ HC → carry info doesn't help (approximation = noise)
  If CCS > HC → carry info HURTS (wrong carry = wrong direction)
  If CCS < HC → carry info HELPS (approximation = useful signal)
""")
