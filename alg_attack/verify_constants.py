"""
ALG Theory Verification - Step 1: Fundamental Constants
Proves: h_carry = 0.811, P(carry) = 1/2, red spectrum, N(c) = 3^{n-transitions}
"""
import random
import math
from collections import Counter

def carry_chain(a, b, n=32):
    """Compute carry chain for a + b (n-bit)."""
    c = [0] * n
    carry = 0
    for i in range(n):
        ai = (a >> i) & 1
        bi = (b >> i) & 1
        # MAJ(ai, bi, carry)
        c[i] = (ai & bi) | (ai & carry) | (bi & carry)
        carry = c[i]
    return c

def transitions(c):
    """Count number of 0->1 and 1->0 transitions in carry chain."""
    t = 0
    prev = 0  # c[-1] = 0
    for bit in c:
        if bit != prev:
            t += 1
        prev = bit
    return t

def entropy(p):
    """Binary entropy H(p)."""
    if p <= 0 or p >= 1:
        return 0.0
    return -p * math.log2(p) - (1 - p) * math.log2(1 - p)

# ============================================================
# TEST 1: h_carry = 0.811
# ============================================================
print("=" * 60)
print("TEST 1: h_carry = 0.811 (carry entropy per bit)")
print("=" * 60)

N = 1_000_000
n = 32
carry_ones = [0] * n

for _ in range(N):
    a = random.getrandbits(n)
    b = random.getrandbits(n)
    c = carry_chain(a, b, n)
    for i in range(n):
        carry_ones[i] += c[i]

print(f"\nSamples: {N:,}")
print(f"\nBit | P(carry=1) | H(carry)")
print("-" * 40)
for i in [0, 1, 2, 3, 4, 15, 31]:
    p = carry_ones[i] / N
    h = entropy(p)
    marker = " <-- transient" if i < 2 else ""
    print(f"  {i:2d} | {p:.6f}   | {h:.6f}{marker}")

# Stationary bits (i >= 2)
p_stationary = sum(carry_ones[2:]) / (N * (n - 2))

# h_carry = CONDITIONAL entropy H(c[i] | c[i-1])
# NOT marginal H(c[i]) which → 1.0
# Measure: P(c[i]=1 | c[i-1]=0) and P(c[i]=1 | c[i-1]=1)
cond_count = [[0,0],[0,0]]  # cond_count[prev][curr]
cond_total = [0, 0]
for _ in range(N):
    a = random.getrandbits(n)
    b = random.getrandbits(n)
    c = carry_chain(a, b, n)
    for i in range(3, n):  # stationary
        prev = c[i-1]
        curr = c[i]
        cond_count[prev][curr] += 1
        cond_total[prev] += 1

p_1_given_0 = cond_count[0][1] / cond_total[0]  # P(c=1|prev=0) = 1/4
p_1_given_1 = cond_count[1][1] / cond_total[1]  # P(c=1|prev=1) = 3/4
h_given_0 = entropy(p_1_given_0)
h_given_1 = entropy(p_1_given_1)
h_conditional = 0.5 * h_given_0 + 0.5 * h_given_1  # weighted by P(prev)

print(f"\nStationary P(carry=1)        = {p_stationary:.6f}  (expected: 0.5000)")
print(f"P(c[i]=1 | c[i-1]=0)        = {p_1_given_0:.6f}  (expected: 0.2500 = 1/4)")
print(f"P(c[i]=1 | c[i-1]=1)        = {p_1_given_1:.6f}  (expected: 0.7500 = 3/4)")
print(f"H(c|prev=0)                 = {h_given_0:.6f}  (expected: 0.8113)")
print(f"H(c|prev=1)                 = {h_given_1:.6f}  (expected: 0.8113)")
print(f"h_carry = H(c[i]|c[i-1])    = {h_conditional:.6f}  (expected: 0.8113)")
print(f"Theoretical h_carry          = {-0.75*math.log2(0.75) - 0.25*math.log2(0.25):.6f}")
h_stationary = h_conditional

# ============================================================
# TEST 2: Red spectrum S(f) = 1/(5/4 - cos(2*pi*f))
# ============================================================
print("\n" + "=" * 60)
print("TEST 2: Carry spectrum (red noise, not white)")
print("=" * 60)

N2 = 100_000
spectrum = [0.0] * (n // 2 + 1)

for _ in range(N2):
    a = random.getrandbits(n)
    b = random.getrandbits(n)
    c = carry_chain(a, b, n)
    # DFT
    for f in range(n // 2 + 1):
        real = sum(c[i] * math.cos(2 * math.pi * f * i / n) for i in range(n))
        imag = sum(c[i] * math.sin(2 * math.pi * f * i / n) for i in range(n))
        spectrum[f] += real * real + imag * imag

for f in range(len(spectrum)):
    spectrum[f] /= N2

print(f"\nSamples: {N2:,}")
print(f"\nFreq | Measured    | Theory 1/(5/4-cos) | Ratio")
print("-" * 60)
for f in [0, 1, 2, 4, 8, 16]:
    if f >= len(spectrum):
        break
    freq = f / n
    theory = 1.0 / (1.25 - math.cos(2 * math.pi * freq)) if f > 0 else 4.0
    ratio = spectrum[f] / theory if theory > 0 else 0
    print(f"  {f:2d}  | {spectrum[f]:10.2f} | {theory:10.4f}          | {ratio:.2f}")

ratio_low_high = spectrum[0] / spectrum[n // 2] if spectrum[n // 2] > 0 else float('inf')
print(f"\nS(0)/S(n/2) = {ratio_low_high:.2f}  (expected: ~9.0 = RED spectrum)")
print(f"White noise would give ratio = 1.0")

# ============================================================
# TEST 3: N(c) = 3^{n - transitions(c)}
# ============================================================
print("\n" + "=" * 60)
print("TEST 3: N(c) = 3^{n - transitions(c)}  [n=8, exhaustive]")
print("=" * 60)

n3 = 8
pair_count = Counter()

for a in range(2 ** n3):
    for b in range(2 ** n3):
        c = carry_chain(a, b, n3)
        c_tuple = tuple(c)
        pair_count[c_tuple] += 1

print(f"\nExhaustive search: {2**n3}×{2**n3} = {2**(2*n3):,} pairs")
print(f"Unique carry patterns: {len(pair_count)}")
print(f"\nCarry pattern       | Trans | N(c) measured | 3^(n-trans) | Match?")
print("-" * 75)

correct = 0
total = 0
for c_tuple in sorted(pair_count.keys(), key=lambda x: pair_count[x], reverse=True)[:10]:
    t = transitions(list(c_tuple))
    measured = pair_count[c_tuple]
    predicted = 3 ** (n3 - t)
    match = "YES" if measured == predicted else "NO"
    if measured == predicted:
        correct += 1
    total += 1
    c_str = ''.join(map(str, c_tuple))
    print(f"  {c_str} | {t:5d} | {measured:13d} | {predicted:11d} | {match}")

# Check ALL patterns
all_correct = all(
    pair_count[c] == 3 ** (n3 - transitions(list(c)))
    for c in pair_count
)
print(f"\nALL {len(pair_count)} patterns match N(c) = 3^(n-trans): {all_correct}")

# ============================================================
# TEST 4: Markov correlation rho_1 = 1/2
# ============================================================
print("\n" + "=" * 60)
print("TEST 4: Markov correlation rho_1 = P(c[i]=c[i+1]) = 3/4")
print("=" * 60)

N4 = 1_000_000
same_count = 0
total_pairs = 0

for _ in range(N4):
    a = random.getrandbits(n)
    b = random.getrandbits(n)
    c = carry_chain(a, b, n)
    for i in range(2, n - 1):  # stationary region
        if c[i] == c[i + 1]:
            same_count += 1
        total_pairs += 1

p_same = same_count / total_pairs
print(f"\nP(c[i] = c[i+1]) = {p_same:.6f}  (expected: 0.7500 = 3/4)")
print(f"Correlation rho_1  = {p_same - 0.5:.6f}  (expected: 0.2500 = 1/4)")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 60)
print("SUMMARY: ALG Constants Verification")
print("=" * 60)
print(f"""
  h_carry = {h_stationary:.4f}     (theory: 0.8113)  {'PASS' if abs(h_stationary - 0.8113) < 0.01 else 'FAIL'}
  P(c=1)  = {p_stationary:.4f}     (theory: 0.5000)  {'PASS' if abs(p_stationary - 0.5) < 0.01 else 'FAIL'}
  S(0)/S(n/2) = {ratio_low_high:.2f}  (theory: ~9.0)    {'PASS' if 7 < ratio_low_high < 12 else 'FAIL'}
  N(c) = 3^(n-t)           (exact match)  {'PASS' if all_correct else 'FAIL'}
  P(same) = {p_same:.4f}     (theory: 0.7500)  {'PASS' if abs(p_same - 0.75) < 0.01 else 'FAIL'}
""")
