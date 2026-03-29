"""
PARALLEL EXPERIMENTS — All remaining directions simultaneously.

1. Invariant subspace of Σ₀·Σ₁ at reduced rounds (2D subspace exploitation)
2. Joint algebra dimension (did it finish?)
3. Da₁₃ as function of MULTIPLE words (not just W[0])
4. Carry-profile clustering for collision search
5. Schedule linear algebra: σ₀·σ₁ eigenstructure for message freedom
6. Quantitative carry-security complementarity verification
"""
import struct
import random
import time
import numpy as np

M32 = 0xFFFFFFFF

def ror(x, n, bits=32):
    return ((x >> n) | (x << (bits - n))) & ((1 << bits) - 1)

K = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0xfc19dc68, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
]
IV = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]

def sha256_R(W, R):
    """SHA-256 reduced to R rounds with feedforward"""
    w = list(W)
    for i in range(16, max(R, 16)):
        s0 = (ror(w[i-15], 7) ^ ror(w[i-15], 18) ^ (w[i-15] >> 3)) & M32
        s1 = (ror(w[i-2], 17) ^ ror(w[i-2], 19) ^ (w[i-2] >> 10)) & M32
        w.append((w[i-16] + s0 + w[i-7] + s1) & M32)
    a,b,c,d,e,f,g,h = IV
    for i in range(R):
        S1 = (ror(e,6)^ror(e,11)^(e>>25))&M32
        ch = ((e&f)^(~e&g))&M32
        t1 = (h+S1+ch+K[i]+w[i])&M32
        S0 = (ror(a,2)^ror(a,13)^ror(a,22))&M32
        maj = ((a&b)^(a&c)^(b&c))&M32
        t2 = (S0+maj)&M32
        h=g;g=f;f=e;e=(d+t1)&M32;d=c;c=b;b=a;a=(t1+t2)&M32
    return tuple((s+iv)&M32 for s,iv in zip((a,b,c,d,e,f,g,h), IV))

def hw(x):
    return bin(x).count('1')

def hw_hash(h1, h2):
    return sum(hw(a^b) for a,b in zip(h1,h2))

random.seed(42)

# ============================================================
# EXP 1: 2D invariant subspace of Σ₀·Σ₁ at reduced rounds
# ============================================================
print("=" * 70)
print("EXP 1: Delta along Σ₀·Σ₁ invariant subspace at reduced rounds")
print("=" * 70)

# The 2D invariant vectors
INV1 = 0x2a811e81
INV2 = 0x80955cba

for R in [17, 18, 19, 20, 22, 24, 64]:
    hw_inv1 = 0
    hw_inv2 = 0
    hw_inv_combo = 0
    hw_msb = 0
    N = 5000
    for _ in range(N):
        W = [random.randint(0, M32) for _ in range(16)]
        H0 = sha256_R(W, R)

        W1 = list(W); W1[0] ^= INV1
        hw_inv1 += hw_hash(H0, sha256_R(W1, R))

        W2 = list(W); W2[0] ^= INV2
        hw_inv2 += hw_hash(H0, sha256_R(W2, R))

        W3 = list(W); W3[0] ^= (INV1 ^ INV2)
        hw_inv_combo += hw_hash(H0, sha256_R(W3, R))

        W4 = list(W); W4[0] ^= 0x80000000
        hw_msb += hw_hash(H0, sha256_R(W4, R))

    print(f"  R={R:2d}: inv1={hw_inv1/N:.1f} inv2={hw_inv2/N:.1f} combo={hw_inv_combo/N:.1f} MSB={hw_msb/N:.1f}")

# ============================================================
# EXP 2: Carry-security complementarity QUANTITATIVE verification
# ============================================================
print("\n" + "=" * 70)
print("EXP 2: Carry-security complementarity S(r) + C(r) = 128")
print("=" * 70)

# For each round r: measure how much information about the collision
# is determined by knowing state[r].
# S(r) ≈ HW(delta_H) averaged over random delta with known state[r]
# C(r) ≈ cost to determine state[r] from H ≈ 128 - S(r)

# Approximation: for delta on W[15] MSB, measure HW(delta_H) at round R
# This gives "security AFTER round R": how many bits still differ
# S(R) = HW(delta_H_fullrounds) given state[R] is known

print(f"{'R':>4} | {'mean HW(ΔH[R:64])':>20} | {'S(R)=128-HW':>15} | {'C(R)=HW':>10}")
print("-" * 60)

N_TEST = 5000
for R_known in [0, 4, 8, 12, 16, 17, 18, 19, 20, 24, 32, 48, 64]:
    # Measure: if attacker knows state at round R_known,
    # how different is ΔH from what can be predicted?
    # Proxy: HW(ΔH) for delta_W[15]=MSB at R_known..64 rounds
    hw_sum = 0
    for _ in range(N_TEST):
        W = [random.randint(0, M32) for _ in range(16)]
        H1 = sha256_R(W, 64)
        W2 = list(W); W2[15] ^= 0x80000000
        H2 = sha256_R(W2, 64)
        hw_sum += hw_hash(H1, H2)

    mean_hw = hw_sum / N_TEST
    S = 128  # Security of full hash (constant, by definition)
    # What changes: per-round contribution
    # Actually compute HW(ΔH) for REDUCED rounds R_known
    hw_reduced = 0
    for _ in range(N_TEST):
        W = [random.randint(0, M32) for _ in range(16)]
        H1 = sha256_R(W, R_known) if R_known > 0 else tuple(IV)
        W2 = list(W); W2[15] ^= 0x80000000
        H2 = sha256_R(W2, R_known) if R_known > 0 else tuple(IV)
        hw_reduced += hw_hash(H1, H2)

    hw_r = hw_reduced / N_TEST
    # S(R) = bits that would need to be "fixed" after round R
    # C(R) = bits already determined by round R
    # hw_r = bits different at round R output
    print(f"{R_known:>4} | {hw_r:>20.2f} | {128 - hw_r/2:>15.2f} | {hw_r/2:>10.2f}")

# ============================================================
# EXP 3: Da₁₃ dependence on W[1] (neutral bit test)
# ============================================================
print("\n" + "=" * 70)
print("EXP 3: Does Da₁₃ depend on W[1]? (Wang independence test)")
print("=" * 70)

def wang_simple(W, dW0=0x80000000):
    """Simplified Wang: just compute delta_a at round 16"""
    W1 = list(W); W2 = list(W)
    W2[0] ^= dW0

    a1,b1,c1,d1,e1,f1,g1,h1 = list(IV)
    a2,b2,c2,d2,e2,f2,g2,h2 = list(IV)

    for r in range(16):
        S1_1=(ror(e1,6)^ror(e1,11)^(e1>>25))&M32
        ch1=((e1&f1)^(~e1&g1))&M32
        t1_1=(h1+S1_1+ch1+K[r]+W1[r])&M32
        S0_1=(ror(a1,2)^ror(a1,13)^ror(a1,22))&M32
        maj1=((a1&b1)^(a1&c1)^(b1&c1))&M32
        t2_1=(S0_1+maj1)&M32

        e1_new=(d1+t1_1)&M32
        a1_new=(t1_1+t2_1)&M32

        if r >= 1:
            S1_2=(ror(e2,6)^ror(e2,11)^(e2>>25))&M32
            ch2=((e2&f2)^(~e2&g2))&M32
            needed=(e1_new-d2-h2-S1_2-ch2-K[r])&M32
            W2[r]=needed

        S1_2=(ror(e2,6)^ror(e2,11)^(e2>>25))&M32
        ch2=((e2&f2)^(~e2&g2))&M32
        t1_2=(h2+S1_2+ch2+K[r]+W2[r])&M32
        S0_2=(ror(a2,2)^ror(a2,13)^ror(a2,22))&M32
        maj2=((a2&b2)^(a2&c2)^(b2&c2))&M32
        t2_2=(S0_2+maj2)&M32

        e2_new=(d2+t1_2)&M32
        a2_new=(t1_2+t2_2)&M32

        h1=g1;g1=f1;f1=e1;e1=e1_new;d1=c1;c1=b1;b1=a1;a1=a1_new
        h2=g2;g2=f2;f2=e2;e2=e2_new;d2=c2;c2=b2;b2=a2;a2=a2_new

    return (a1-a2)&M32  # Da at round 16

# Test: fix W[0], vary W[1], measure how Da changes
W_base = [random.randint(0, M32) for _ in range(16)]
da_base = wang_simple(W_base)

n_same = 0
n_test = 10000
da_by_w1 = set()
for _ in range(n_test):
    W_test = list(W_base)
    W_test[1] = random.randint(0, M32)  # Vary W[1] only
    da = wang_simple(W_test)
    da_by_w1.add(da)
    if da == da_base:
        n_same += 1

print(f"  Fix W[0], vary W[1]: {len(da_by_w1)} unique Da out of {n_test}")
print(f"  Da unchanged: {n_same}/{n_test} ({100*n_same/n_test:.1f}%)")

# Test W[12] (neutral bit)
da_by_w12 = set()
n_same_12 = 0
for _ in range(n_test):
    W_test = list(W_base)
    W_test[12] = random.randint(0, M32)
    da = wang_simple(W_test)
    da_by_w12.add(da)
    if da == da_base:
        n_same_12 += 1

print(f"  Fix W[0], vary W[12]: {len(da_by_w12)} unique Da out of {n_test}")
print(f"  Da unchanged: {n_same_12}/{n_test} ({100*n_same_12/n_test:.1f}%)")

# Test W[15] (should be very influential via schedule)
da_by_w15 = set()
n_same_15 = 0
for _ in range(n_test):
    W_test = list(W_base)
    W_test[15] = random.randint(0, M32)
    da = wang_simple(W_test)
    da_by_w15.add(da)
    if da == da_base:
        n_same_15 += 1

print(f"  Fix W[0], vary W[15]: {len(da_by_w15)} unique Da out of {n_test}")
print(f"  Da unchanged: {n_same_15}/{n_test} ({100*n_same_15/n_test:.1f}%)")

# ============================================================
# EXP 4: Carry profile correlation between messages
# ============================================================
print("\n" + "=" * 70)
print("EXP 4: Carry profile similarity for structured vs random pairs")
print("=" * 70)

def get_carry_profile(W, R=64):
    """Get carry overflow bits for each round"""
    w = list(W)
    for i in range(16, max(R, 16)):
        s0 = (ror(w[i-15], 7) ^ ror(w[i-15], 18) ^ (w[i-15] >> 3)) & M32
        s1 = (ror(w[i-2], 17) ^ ror(w[i-2], 19) ^ (w[i-2] >> 10)) & M32
        w.append((w[i-16] + s0 + w[i-7] + s1) & M32)

    a,b,c,d,e,f,g,h = IV
    carries = []
    for i in range(min(R, 64)):
        S1 = (ror(e,6)^ror(e,11)^(e>>25))&M32
        ch = ((e&f)^(~e&g))&M32
        raw = h + S1 + ch + K[i] + w[i]
        carries.append(1 if raw >= (1 << 32) else 0)
        t1 = raw & M32
        S0 = (ror(a,2)^ror(a,13)^ror(a,22))&M32
        maj = ((a&b)^(a&c)^(b&c))&M32
        t2 = (S0+maj)&M32
        h=g;g=f;f=e;e=(d+t1)&M32;d=c;c=b;b=a;a=(t1+t2)&M32
    return carries

# Compare carry profiles
N_PAIRS = 5000

# Same W[0], different W[1..15]
same_w0_match = 0
for _ in range(N_PAIRS):
    w0 = random.randint(0, M32)
    W1 = [w0] + [random.randint(0, M32) for _ in range(15)]
    W2 = [w0] + [random.randint(0, M32) for _ in range(15)]
    c1 = get_carry_profile(W1)
    c2 = get_carry_profile(W2)
    same_w0_match += sum(a == b for a, b in zip(c1, c2))

# Random pairs
random_match = 0
for _ in range(N_PAIRS):
    W1 = [random.randint(0, M32) for _ in range(16)]
    W2 = [random.randint(0, M32) for _ in range(16)]
    c1 = get_carry_profile(W1)
    c2 = get_carry_profile(W2)
    random_match += sum(a == b for a, b in zip(c1, c2))

# Adjacent W[0]
adj_match = 0
for _ in range(N_PAIRS):
    W1 = [random.randint(0, M32) for _ in range(16)]
    W2 = list(W1); W2[0] = (W1[0] + 1) & M32
    c1 = get_carry_profile(W1)
    c2 = get_carry_profile(W2)
    adj_match += sum(a == b for a, b in zip(c1, c2))

total_bits = N_PAIRS * 64
print(f"  Same W[0], diff rest: {same_w0_match}/{total_bits} match ({100*same_w0_match/total_bits:.1f}%)")
print(f"  Random pairs:         {random_match}/{total_bits} match ({100*random_match/total_bits:.1f}%)")
print(f"  Adjacent W[0]:        {adj_match}/{total_bits} match ({100*adj_match/total_bits:.1f}%)")
print(f"  Expected random:      {50.0:.1f}%")

# ============================================================
# EXP 5: Schedule freedom — how many message bits are truly free?
# ============================================================
print("\n" + "=" * 70)
print("EXP 5: Schedule independence structure")
print("=" * 70)

# Which W[j] does delta_e[17] actually depend on?
# Vary each W[j] individually, check if delta_e17 changes

W_base = [random.randint(0, M32) for _ in range(16)]
da_base = wang_simple(W_base)

print("  Sensitivity of Da₁₆ to each W[j]:")
for j in range(16):
    n_diff = 0
    n_test_j = 1000
    for _ in range(n_test_j):
        W_test = list(W_base)
        W_test[j] = random.randint(0, M32)
        da = wang_simple(W_test)
        if da != da_base:
            n_diff += 1
    print(f"    W[{j:2d}]: Da changes in {n_diff}/{n_test_j} ({100*n_diff/n_test_j:5.1f}%) cases")

# ============================================================
# EXP 6: Multi-delta: XOR delta on MULTIPLE words simultaneously
# ============================================================
print("\n" + "=" * 70)
print("EXP 6: Multi-word delta — HW(ΔH) for delta on multiple words")
print("=" * 70)

N_TEST = 5000
for n_words in [1, 2, 4, 8, 16]:
    hw_sum = 0
    best = 256
    for _ in range(N_TEST):
        W = [random.randint(0, M32) for _ in range(16)]
        H1 = sha256_R(W, 64)

        W2 = list(W)
        # XOR MSB on first n_words
        for j in range(n_words):
            W2[j] ^= 0x80000000
        H2 = sha256_R(W2, 64)
        d = hw_hash(H1, H2)
        hw_sum += d
        best = min(best, d)

    print(f"  Delta on W[0..{n_words-1}] MSB: mean={hw_sum/N_TEST:.1f}, best={best}")
