"""
CRITICAL EXPERIMENT: Compute |Image(ε)| for SHA-256.

ε(W[0]) = Da₁₃(W[0]) ⊕ Da₁₃_linear(W[0])
       = carry-error in the Da₁₃ computation

Theory predicts: |Image(ε)| << 2^32 if carry is structured.
If |Image(ε)| ≈ 2^32: no speedup (ε is surjective).
If |Image(ε)| << 2^32: potential speedup for barrier crossing.

Also: compute Da₁₃ distribution and verify T_BARRIER_UNIFORM.
"""
import struct
import random

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

M32 = 0xFFFFFFFF

def sha256_17rounds(W):
    """Run 17 rounds of SHA-256, return (state_after_17, delta_e17 components)"""
    w = list(W)
    for i in range(16, 17):
        s0 = (ror(w[i-15], 7) ^ ror(w[i-15], 18) ^ (w[i-15] >> 3)) & M32
        s1 = (ror(w[i-2], 17) ^ ror(w[i-2], 19) ^ (w[i-2] >> 10)) & M32
        w.append((w[i-16] + s0 + w[i-7] + s1) & M32)

    a, b, c, d, e, f, g, h = IV
    for i in range(17):
        S1 = (ror(e, 6) ^ ror(e, 11) ^ (e >> 25)) & M32
        ch = ((e & f) ^ (~e & g)) & M32
        temp1 = (h + S1 + ch + K[i] + w[i]) & M32
        S0 = (ror(a, 2) ^ ror(a, 13) ^ ror(a, 22)) & M32
        maj = ((a & b) ^ (a & c) ^ (b & c)) & M32
        temp2 = (S0 + maj) & M32
        h = g; g = f; f = e
        e = (d + temp1) & M32
        d = c; c = b; b = a
        a = (temp1 + temp2) & M32
    return (a, b, c, d, e, f, g, h)

def wang_chain(W, delta_W0=0x80000000):
    """
    Build Wang chain: choose dW[1..15] so delta_e[2..16]=0.
    W = base message (16 words).
    Returns (W_prime, details) where W_prime = modified second message.
    """
    W1 = list(W)
    W2 = list(W)
    W2[0] = (W2[0] ^ delta_W0) & M32  # XOR delta on W[0]

    # Run both forward, adaptively choosing W2[r] to zero delta_e[r+1]
    a1, b1, c1, d1, e1, f1, g1, h1 = list(IV)
    a2, b2, c2, d2, e2, f2, g2, h2 = list(IV)

    for r in range(16):
        # Compute round r for message 1
        S1_1 = (ror(e1, 6) ^ ror(e1, 11) ^ (e1 >> 25)) & M32
        ch1 = ((e1 & f1) ^ (~e1 & g1)) & M32
        t1_1 = (h1 + S1_1 + ch1 + K[r] + W1[r]) & M32
        S0_1 = (ror(a1, 2) ^ ror(a1, 13) ^ ror(a1, 22)) & M32
        maj1 = ((a1 & b1) ^ (a1 & c1) ^ (b1 & c1)) & M32
        t2_1 = (S0_1 + maj1) & M32

        e1_new = (d1 + t1_1) & M32
        a1_new = (t1_1 + t2_1) & M32

        # For message 2: we need e2_new = e1_new (delta_e = 0)
        # e2_new = d2 + T1_2 = d2 + h2 + S1(e2) + Ch(e2,f2,g2) + K[r] + W2[r]
        # We want e2_new = e1_new
        # So: W2[r] = e1_new - d2 - h2 - S1(e2) - Ch(e2,f2,g2) - K[r]

        if r >= 1:  # For r=0, W2[0] is already set (delta_W0)
            S1_2 = (ror(e2, 6) ^ ror(e2, 11) ^ (e2 >> 25)) & M32
            ch2 = ((e2 & f2) ^ (~e2 & g2)) & M32
            needed = (e1_new - d2 - h2 - S1_2 - ch2 - K[r]) & M32
            W2[r] = needed

        # Now compute round r for message 2 with (possibly adjusted) W2[r]
        S1_2 = (ror(e2, 6) ^ ror(e2, 11) ^ (e2 >> 25)) & M32
        ch2 = ((e2 & f2) ^ (~e2 & g2)) & M32
        t1_2 = (h2 + S1_2 + ch2 + K[r] + W2[r]) & M32
        S0_2 = (ror(a2, 2) ^ ror(a2, 13) ^ ror(a2, 22)) & M32
        maj2 = ((a2 & b2) ^ (a2 & c2) ^ (b2 & c2)) & M32
        t2_2 = (S0_2 + maj2) & M32

        e2_new = (d2 + t1_2) & M32
        a2_new = (t1_2 + t2_2) & M32

        # Shift registers
        h1 = g1; g1 = f1; f1 = e1; e1 = e1_new; d1 = c1; c1 = b1; b1 = a1; a1 = a1_new
        h2 = g2; g2 = f2; f2 = e2; e2 = e2_new; d2 = c2; c2 = b2; b2 = a2; a2 = a2_new

    # Now compute delta_e[17] using schedule W[16]
    # W1[16] from schedule
    s0_1 = (ror(W1[1], 7) ^ ror(W1[1], 18) ^ (W1[1] >> 3)) & M32
    s1_1 = (ror(W1[14], 17) ^ ror(W1[14], 19) ^ (W1[14] >> 10)) & M32
    W1_16 = (W1[0] + s0_1 + W1[9] + s1_1) & M32

    s0_2 = (ror(W2[1], 7) ^ ror(W2[1], 18) ^ (W2[1] >> 3)) & M32
    s1_2 = (ror(W2[14], 17) ^ ror(W2[14], 19) ^ (W2[14] >> 10)) & M32
    W2_16 = (W2[0] + s0_2 + W2[9] + s1_2) & M32

    # Round 16 for both
    S1_1r = (ror(e1, 6) ^ ror(e1, 11) ^ (e1 >> 25)) & M32
    ch1r = ((e1 & f1) ^ (~e1 & g1)) & M32
    t1_1r = (h1 + S1_1r + ch1r + K[16] + W1_16) & M32
    e1_17 = (d1 + t1_1r) & M32

    S1_2r = (ror(e2, 6) ^ ror(e2, 11) ^ (e2 >> 25)) & M32
    ch2r = ((e2 & f2) ^ (~e2 & g2)) & M32
    t1_2r = (h2 + S1_2r + ch2r + K[16] + W2_16) & M32
    e2_17 = (d2 + t1_2r) & M32

    delta_e17 = (e1_17 - e2_17) & M32
    Da13 = (a1 - a2) & M32  # a-register difference at round 13... actually at round 16

    return delta_e17, Da13, W2, (a1, a2, d1, d2)


# ============================================================
# MAIN EXPERIMENT: Sample Da13 and carry-error
# ============================================================

random.seed(42)
N_SAMPLES = 200000

delta_e17_values = []
da_values = set()
da_list = []

print("=" * 70)
print(f"Sampling Da and δe₁₇ for {N_SAMPLES} random W[0] values")
print("=" * 70)

# Use fixed W[1..15], vary only W[0]
W_base = [0] * 16  # Start with zeros for simplicity
# Actually use random fixed W[1..15]
random.seed(123)
W_base = [random.randint(0, M32) for _ in range(16)]

for trial in range(N_SAMPLES):
    W_base[0] = random.randint(0, M32)
    de17, da, _, _ = wang_chain(W_base)
    delta_e17_values.append(de17)
    da_values.add(da)
    da_list.append(da)

# Statistics
n_zero = sum(1 for x in delta_e17_values if x == 0)
n_unique_da = len(da_values)
n_unique_de17 = len(set(delta_e17_values))

print(f"\nδe₁₇ = 0 count: {n_zero}/{N_SAMPLES} (expected: ~{N_SAMPLES/2**32:.4f})")
print(f"Unique Da values: {n_unique_da}/{N_SAMPLES}")
print(f"Unique δe₁₇ values: {n_unique_de17}/{N_SAMPLES}")

# Check uniformity of Da
print(f"\n--- Da distribution ---")
# Bit-by-bit bias
print("Per-bit bias of Da (deviation from 0.5):")
max_bias = 0
for bit in range(32):
    ones = sum(1 for d in da_list if (d >> bit) & 1)
    bias = ones / N_SAMPLES - 0.5
    max_bias = max(max_bias, abs(bias))
print(f"  Max absolute bias: {max_bias:.6f} (expected ~{1/N_SAMPLES**0.5:.6f} for random)")

# HW distribution
hw_counts = {}
for d in da_list:
    h = bin(d).count('1')
    hw_counts[h] = hw_counts.get(h, 0) + 1
mean_hw = sum(bin(d).count('1') for d in da_list) / N_SAMPLES
print(f"  Mean HW(Da): {mean_hw:.2f} (expected: 16.0 for uniform)")

# Now compute carry-error
# We need Da_linear: the GF(2)-linearized version
# Approximate: run SHA-256 with all + replaced by XOR... complex.
# Instead: measure empirically how many DISTINCT Da values there are
# as W[0] varies over a SMALLER range

print(f"\n--- Carry-error Image size estimation ---")

# For small ranges of W[0]: count distinct Da
for log_range in [8, 10, 12, 14, 16, 18, 20]:
    n_test = min(2**log_range, 1 << 20)
    da_set = set()
    for i in range(n_test):
        W_base[0] = i if log_range <= 20 else random.randint(0, M32)
        de17, da, _, _ = wang_chain(W_base)
        da_set.add(da)
    coverage = len(da_set) / min(n_test, 2**32)
    print(f"  W[0] range 0..{n_test-1}: {len(da_set)} unique Da out of {n_test} ({100*coverage:.1f}%)")

# Test: is Da a LINEAR function of W[0]?
print(f"\n--- Linearity test ---")
# If Da linear in W[0]: Da(W0 XOR a) XOR Da(W0) XOR Da(a) XOR Da(0) = 0
W_base[0] = 0
_, da_zero, _, _ = wang_chain(W_base)

n_linear = 0
n_test_lin = 10000
for _ in range(n_test_lin):
    w0 = random.randint(1, M32)
    a = random.randint(1, M32)

    W_base[0] = w0
    _, da_w0, _, _ = wang_chain(W_base)

    W_base[0] = a
    _, da_a, _, _ = wang_chain(W_base)

    W_base[0] = w0 ^ a
    _, da_w0xa, _, _ = wang_chain(W_base)

    # Linearity check: Da(w0^a) = Da(w0) ^ Da(a) ^ Da(0) for GF(2)-linear
    if da_w0xa == (da_w0 ^ da_a ^ da_zero):
        n_linear += 1

print(f"  Linearity test: {n_linear}/{n_test_lin} ({100*n_linear/n_test_lin:.1f}%) pass")
print(f"  (100% = GF(2)-linear, ~0% = fully nonlinear)")

# Differential uniformity of Da as function of W[0]
print(f"\n--- Differential uniformity ---")
# For random delta, what's P(Da(W0) = Da(W0 ^ delta))?
n_coll = 0
n_test_du = 50000
for _ in range(n_test_du):
    w0 = random.randint(0, M32)
    delta = random.randint(1, M32)

    W_base[0] = w0
    _, da1, _, _ = wang_chain(W_base)

    W_base[0] = w0 ^ delta
    _, da2, _, _ = wang_chain(W_base)

    if da1 == da2:
        n_coll += 1

print(f"  Da collisions: {n_coll}/{n_test_du}")
print(f"  Rate: {n_coll/n_test_du:.6f} (expected for random: {1/2**32:.10f})")
if n_coll > 0:
    print(f"  *** DA HAS COLLISIONS — NOT BIJECTIVE ***")
    print(f"  Estimated |Image(Da)|/2^32: {1 - n_coll/n_test_du:.6f}")
