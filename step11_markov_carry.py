#!/usr/bin/env python3
"""
Step 11: Markov structure of carry propagation — M(0.313) connection.

From the methodology: carry probability was measured as 0.313.
Our fundamental identity shows: carry(x, y) = (x+y) ⊕ x ⊕ y
And carry has geometric distribution: P(carry_chain_length ≥ k) = (1/2)^k.

But 0.313 ≠ 1/2. Where does 0.313 come from?

Hypothesis: 0.313 is the probability that a SPECIFIC bit position
has carry=1 AFTER the full SHA-256 round (not just one addition).
The round involves ~5 additions, and their carries interact.

Key questions:
1. What is the per-bit carry probability in SHA-256?
2. Does it form a Markov chain?
3. What is the transition matrix?
4. How does this relate to barrier structure?
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

DW_BIT = 0x80000000

def sha_round_fn(state, W_r, K_r):
    a,b,c,d,e,f,g,h = state
    T1 = add(add(add(add(h, Sig1(e)), Ch(e,f,g)), K_r), W_r)
    T2 = add(Sig0(a), Maj(a,b,c))
    return [add(T1, T2), a, b, c, add(d, T1), e, f, g]

# ============================================================
# PART A: Per-bit carry probability in T1 computation
# ============================================================
print("=" * 70)
print("PART A: Carry probability in T1 = h + Sig1(e) + Ch(e,f,g) + K + W")
print("=" * 70)
print()

# T1 involves 4 additions. The total carry from all 4 additions:
# total_carry = (T1_actual) XOR (h XOR Sig1(e) XOR Ch(e,f,g) XOR K XOR W)
# where the XOR sum is the "no-carry" answer.

N = 100000
carry_per_bit = np.zeros(32)

for _ in range(N):
    state = [random.randint(0, MASK) for _ in range(8)]
    W_r = random.randint(0, MASK)
    a,b,c,d,e,f,g,h = state

    # XOR sum (no carry):
    xor_sum = h ^ Sig1(e) ^ Ch(e,f,g) ^ K[0] ^ W_r
    # Actual sum (with carry):
    T1_actual = add(add(add(add(h, Sig1(e)), Ch(e,f,g)), K[0]), W_r)
    # Carry pattern:
    carry = T1_actual ^ xor_sum

    for bit in range(32):
        carry_per_bit[bit] += (carry >> bit) & 1

carry_prob = carry_per_bit / N

print(f"Per-bit carry probability in T1 (N={N}):")
for i in range(0, 32, 8):
    bits = range(i, min(i+8, 32))
    vals = [f"{carry_prob[b]:.3f}" for b in bits]
    print(f"  bits {i:2d}-{min(i+7,31):2d}: {' '.join(vals)}")

avg_carry = np.mean(carry_prob)
print(f"\nAverage per-bit carry probability: {avg_carry:.4f}")
print(f"Methodology value: 0.313")
print(f"Match: {'YES' if abs(avg_carry - 0.313) < 0.01 else 'NO'} (diff = {abs(avg_carry-0.313):.4f})")
print()

# ============================================================
# PART B: Carry probability by number of additions
# ============================================================
print("=" * 70)
print("PART B: Carry probability vs number of additions")
print("=" * 70)
print()

# Theory: for k additions of random 32-bit numbers:
# E[carry at bit i] ≈ (k-1) × P(carry per addition)
# where P(carry per addition) starts at 0 for bit 0 and
# approaches 1/2 for higher bits.
#
# More precisely: for random x, y:
#   P(carry at bit i of x+y) = P(x[0..i]+y[0..i] ≥ 2^(i+1))
#   = sum_{j=0}^{i} P(bit j generates carry and it propagates)

# Test carry probability for different numbers of additions:
print(f"{'# additions':>12} | {'E[carry prob]':>13} | {'Theoretical':>12}")
print("-" * 45)

for k in [1, 2, 3, 4, 5]:
    carry_bits = np.zeros(32)
    for _ in range(50000):
        terms = [random.randint(0, MASK) for _ in range(k+1)]
        xor_sum = 0
        for t in terms:
            xor_sum ^= t
        add_sum = 0
        for t in terms:
            add_sum = add(add_sum, t)
        carry = add_sum ^ xor_sum
        for bit in range(32):
            carry_bits[bit] += (carry >> bit) & 1
    avg = np.mean(carry_bits / 50000)
    # Theoretical: for k additions of uniform random,
    # the carry probability converges to k/(2k+2) for high bits
    # Actually for 2 terms: P(carry at bit i) → 1/2 as i→∞
    # For k+1 terms: P ≈ 1 - 2^(-k) / (k+1)... not simple.
    # Let's compute: exact P for bit 31 (MSB)
    print(f"  {k+1:>10d} | {avg:>13.4f} |")

print()

# ============================================================
# PART C: T1 has 5 terms, T2 has 2 terms
# ============================================================
print("=" * 70)
print("PART C: Carry in full round — T1 (5 terms) + T2 (2 terms)")
print("=" * 70)
print()

# T1 = h + Sig1(e) + Ch(e,f,g) + K + W  → 5 terms, 4 additions
# T2 = Sig0(a) + Maj(a,b,c)              → 2 terms, 1 addition
# a' = T1 + T2                            → 1 more addition
# e' = d + T1                             → 1 more addition
#
# Total: 7 additions per round
# But they're CHAINED, not independent.

# Measure: carry in a' = T1 + T2
N2 = 100000
carry_a_prime = np.zeros(32)
carry_e_prime = np.zeros(32)

for _ in range(N2):
    state = [random.randint(0, MASK) for _ in range(8)]
    W_r = random.randint(0, MASK)
    a,b,c,d,e,f,g,h = state

    T1 = add(add(add(add(h, Sig1(e)), Ch(e,f,g)), K[0]), W_r)
    T2 = add(Sig0(a), Maj(a,b,c))

    a_prime = add(T1, T2)
    e_prime = add(d, T1)

    # Carry in a' = T1 + T2
    carry_a = (T1 + T2) ^ (T1 ^ T2)  # would give wrong result mod 2^32
    # Better: carry = (T1+T2 mod 2^32) XOR T1 XOR T2
    # Wait: (x+y) XOR x XOR y = carry pattern (shifted)
    carry_a_val = a_prime ^ T1 ^ T2
    carry_e_val = e_prime ^ d ^ T1

    for bit in range(32):
        carry_a_prime[bit] += (carry_a_val >> bit) & 1
        carry_e_prime[bit] += (carry_e_val >> bit) & 1

avg_carry_a = np.mean(carry_a_prime / N2)
avg_carry_e = np.mean(carry_e_prime / N2)

print(f"Carry probability in a' = T1 + T2: {avg_carry_a:.4f}")
print(f"Carry probability in e' = d + T1:  {avg_carry_e:.4f}")
print()

# ============================================================
# PART D: The TOTAL carry probability per round
# ============================================================
print("=" * 70)
print("PART D: Total per-bit carry — GF(2) vs Z gap per round")
print("=" * 70)
print()

# For the DIFFERENTIAL: what matters is carry in the DIFFERENCE.
# Da' = Da'(actual) = a'_f - a'_n (additive difference)
# vs Da'(GF2) = a'_f XOR a'_n (XOR difference)
#
# The GAP between these is the carry effect.
# Gap = (a'_f - a'_n) XOR (a'_f XOR a'_n)

N3 = 50000
gap_per_bit = np.zeros(32)

for _ in range(N3):
    W_r = random.randint(0, MASK)
    state_n = [random.randint(0, MASK) for _ in range(8)]
    # Small perturbation
    state_f = list(state_n)
    state_f[4] ^= (1 << random.randint(0, 31))  # flip one bit of e

    out_n = sha_round_fn(state_n, W_r, K[0])
    out_f = sha_round_fn(state_f, W_r, K[0])

    # Additive diff
    Da_add = (out_f[0] - out_n[0]) & MASK
    # XOR diff
    Da_xor = out_f[0] ^ out_n[0]
    # Gap
    gap = Da_add ^ Da_xor

    for bit in range(32):
        gap_per_bit[bit] += (gap >> bit) & 1

gap_prob = gap_per_bit / N3

print(f"Per-bit gap probability (additive vs XOR differential):")
for i in range(0, 32, 8):
    bits = range(i, min(i+8, 32))
    vals = [f"{gap_prob[b]:.3f}" for b in bits]
    print(f"  bits {i:2d}-{min(i+7,31):2d}: {' '.join(vals)}")

avg_gap = np.mean(gap_prob)
print(f"\nAverage per-bit gap probability: {avg_gap:.4f}")
print(f"Methodology value M(0.313): {0.313:.3f}")
print()

# ============================================================
# PART E: Markov chain analysis — bit-to-bit carry transition
# ============================================================
print("=" * 70)
print("PART E: Markov chain — carry transition probabilities")
print("=" * 70)
print()

# For T1 = h + Sig1(e) + Ch + K + W:
# carry[i] depends on carry[i-1] and the input bits at position i.
# If inputs are random: P(carry[i]=1 | carry[i-1]=0) = p_01
#                        P(carry[i]=1 | carry[i-1]=1) = p_11
#
# For sum of 5 random terms:
# At bit i, the sum of 5 bits + carry_in can be 0..6.
# carry_out = 1 iff sum ≥ 2 (since we're doing binary addition)
# Wait, this isn't right. We're doing sequential additions, not parallel.

# For sequential addition: a+b+c+d+e = ((((a+b)+c)+d)+e)
# The carries are CORRELATED between additions.
# But the FINAL carry (after all 4 additions) can be modeled.

# Let's measure the Markov transition directly:
# carry[i] = bit i of ((x+y) XOR x XOR y) >> 1
# For the sum of 5 random terms, track carry at each bit position.

N_mk = 200000
# Transition counts for the ROUND carry
# State: (carry_in at bit i) → (carry_out at bit i)
transitions = np.zeros((2, 2))

for _ in range(N_mk):
    terms = [random.randint(0, MASK) for _ in range(5)]
    xor_sum = 0
    for t in terms:
        xor_sum ^= t
    add_sum = 0
    for t in terms:
        add_sum = add(add_sum, t)

    carry_pattern = add_sum ^ xor_sum  # bits where carry had net effect

    # Look at consecutive bits
    for bit in range(1, 31):  # avoid edges
        c_prev = (carry_pattern >> (bit-1)) & 1
        c_curr = (carry_pattern >> bit) & 1
        transitions[c_prev, c_curr] += 1

# Normalize
for i in range(2):
    row_sum = transitions[i].sum()
    if row_sum > 0:
        transitions[i] /= row_sum

print("Markov transition matrix for carry propagation (5-term sum):")
print(f"  P(carry[i]=0 | carry[i-1]=0) = {transitions[0,0]:.4f}")
print(f"  P(carry[i]=1 | carry[i-1]=0) = {transitions[0,1]:.4f}")
print(f"  P(carry[i]=0 | carry[i-1]=1) = {transitions[1,0]:.4f}")
print(f"  P(carry[i]=1 | carry[i-1]=1) = {transitions[1,1]:.4f}")
print()

# Stationary distribution
# π₁ = P(carry=1) = P(0→1) / (P(0→1) + P(1→0))
p01 = transitions[0, 1]
p10 = transitions[1, 0]
stationary_carry = p01 / (p01 + p10)
print(f"Stationary carry probability: {stationary_carry:.4f}")
print(f"Methodology value: 0.313")
print(f"Match: {abs(stationary_carry - 0.313) < 0.02}")
print()

# ============================================================
# PART F: Does 0.313 come from 5-term addition?
# ============================================================
print("=" * 70)
print("PART F: Carry probability by number of terms in addition")
print("=" * 70)
print()

for n_terms in range(2, 8):
    carry_count = 0
    total_bits = 0
    for _ in range(50000):
        terms = [random.randint(0, MASK) for _ in range(n_terms)]
        xor_sum = 0
        for t in terms:
            xor_sum ^= t
        add_sum = 0
        for t in terms:
            add_sum = add(add_sum, t)
        carry_pattern = add_sum ^ xor_sum
        carry_count += hw(carry_pattern)
        total_bits += 32
    p = carry_count / total_bits
    print(f"  {n_terms} terms: P(carry) = {p:.4f}  {'← M(0.313)?' if abs(p-0.313) < 0.01 else ''}")

print()

# ============================================================
# PART G: Carry in the DIFFERENTIAL round function
# ============================================================
print("=" * 70)
print("PART G: Carry in differential propagation")
print("=" * 70)
print()

# The carry that matters for the barrier is in the DIFFERENTIAL:
# Da[r+1] = DT1 + DT2  (additive)
# The carry in DT1 + DT2 is the carry of ONE addition (2 terms).
# But DT1 itself involves multiple additions → carry accumulated.
#
# What is the EFFECTIVE number of carry-producing additions
# that the DIFFERENTIAL sees?

# Measure: per-bit carry in the additive differential of a full round
N4 = 100000
diff_carry = np.zeros(32)

for _ in range(N4):
    state = [random.randint(0, MASK) for _ in range(8)]
    W = random.randint(0, MASK)

    # Perturbed state: flip bit 31 of e
    state_f = list(state)
    state_f[4] ^= DW_BIT

    out_n = sha_round_fn(state, W, K[0])
    out_f = sha_round_fn(state_f, W, K[0])

    # Da' additive
    Da_add = (out_f[0] - out_n[0]) & MASK
    # Da' XOR
    Da_xor = out_f[0] ^ out_n[0]
    # Carry = difference between additive and XOR
    carry = Da_add ^ Da_xor

    for bit in range(32):
        diff_carry[bit] += (carry >> bit) & 1

diff_carry_prob = diff_carry / N4

print(f"Per-bit carry in differential Da[r+1]:")
print(f"(Perturbation: δe[bit 31] = 1)")
print()
for i in range(0, 32, 8):
    bits = range(i, min(i+8, 32))
    vals = [f"{diff_carry_prob[b]:.3f}" for b in bits]
    print(f"  bits {i:2d}-{min(i+7,31):2d}: {' '.join(vals)}")

avg_diff_carry = np.mean(diff_carry_prob)
print(f"\nAverage: {avg_diff_carry:.4f}")
print(f"This is the carry probability in the DIFFERENTIAL,")
print(f"not in the absolute computation.")
print()

# Now measure across multiple rounds — is there a Markov structure?
print("Carry probability evolution across rounds:")
for R in [1, 2, 4, 8, 16]:
    carry_total = 0
    bit_total = 0
    for _ in range(10000):
        state_n = list(IV)
        state_f = list(IV)
        state_f[4] ^= DW_BIT
        W16 = [random.randint(0, MASK) for _ in range(16)]
        from functools import reduce

        for r in range(R):
            state_n = sha_round_fn(state_n, W16[r] if r < 16 else 0, K[r])
            state_f = sha_round_fn(state_f, W16[r] if r < 16 else 0, K[r])

        Da_add = (state_f[0] - state_n[0]) & MASK
        Da_xor = state_f[0] ^ state_n[0]
        carry = Da_add ^ Da_xor
        carry_total += hw(carry)
        bit_total += 32

    p = carry_total / bit_total
    print(f"  After {R:2d} rounds: P(carry in Da) = {p:.4f}")

print()
print("=" * 70)
print("SUMMARY STEP 11")
print("=" * 70)
print()
print("CARRY PROBABILITY FINDINGS:")
print(f"  T1 (5-term sum): P(carry) ≈ {avg_carry:.4f}")
print(f"  Differential Da: P(carry) ≈ {avg_diff_carry:.4f}")
print(f"  Markov stationary: {stationary_carry:.4f}")
print(f"  Methodology M(0.313): 0.313")
print()
print("The 0.313 value appears to correspond to the carry probability")
print("in a specific SHA-256 sub-computation, not the full round.")
print("The Markov model captures carry propagation between bit positions.")
