#!/usr/bin/env python3
"""
NK Direction 7: Schedule as error-correcting code.

SHA-256 schedule: W[0..15] → W[16..63] (linear over GF(2) for XOR diffs).
This is a [1536, 512] linear code over GF(2).

Key question: what is the MINIMUM DISTANCE of this code?
min_dist = min HW(Σ(δM)) for nonzero δM.

If min_dist is small → there exist δM producing very sparse δW[16..63].
Sparse δW → few active rounds → potentially cheap differential path.

From NK §3: Λ(path) depends on number of ACTIVE rounds.
Active round = round where δW[r] ≠ 0.
min_dist / 32 ≈ minimum number of active WORDS.

If min_dist < 128: there exist δM with < 4 active words in W[16..63].
Combined with Wang (16 free rounds): potentially < 128 total Λ.
"""

import random, math, time

MASK32 = 0xFFFFFFFF

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&MASK32
def ssig0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def ssig1(x): return rotr(x,17)^rotr(x,19)^(x>>10)

def expand_schedule(M):
    W=list(M)+[0]*(64-len(M))
    for i in range(16,64): W[i]=(ssig1(W[i-2])+W[i-7]+ssig0(W[i-15])+W[i-16])&MASK32
    return W

def expand_schedule_xor(M):
    """XOR-linearized schedule (no carry in additions)."""
    W=list(M)+[0]*(64-len(M))
    for i in range(16,64):
        W[i] = ssig1(W[i-2]) ^ W[i-7] ^ ssig0(W[i-15]) ^ W[i-16]
    return W

def hw(x): return bin(x & MASK32).count('1')

def hw_schedule(M, xor_model=True):
    """Total HW of W[16..63] for given M (as δM from zero)."""
    if xor_model:
        W = expand_schedule_xor(M)
    else:
        W = expand_schedule(M)
    return sum(hw(W[r]) for r in range(16, 64))

def active_words(M, xor_model=True):
    """Count nonzero words in W[16..63]."""
    if xor_model:
        W = expand_schedule_xor(M)
    else:
        W = expand_schedule(M)
    return sum(1 for r in range(16, 64) if W[r] != 0)


# ================================================================
# TEST 1: Single-bit inputs — HW of schedule expansion
# For each input bit (word w, bit b): what's HW(δW[16..63])?
# ================================================================

print("=" * 70)
print("TEST 1: HW(δW[16..63]) for single-bit δM")
print("=" * 70)

results_1bit = []
for w in range(16):
    for b in range(32):
        M = [0]*16
        M[w] = 1 << b
        h = hw_schedule(M, xor_model=True)
        aw = active_words(M, xor_model=True)
        results_1bit.append((h, aw, w, b))

results_1bit.sort()

print(f"\n  512 single-bit inputs, XOR model")
print(f"\n  Top 10 LOWEST HW (sparsest schedule diffs):")
print(f"  {'HW':>4} {'active_words':>12} | input")
for h, aw, w, b in results_1bit[:10]:
    print(f"  {h:>4} {aw:>12} | W[{w}] bit {b}")

print(f"\n  Top 5 HIGHEST HW:")
for h, aw, w, b in results_1bit[-5:]:
    print(f"  {h:>4} {aw:>12} | W[{w}] bit {b}")

print(f"\n  Stats: min={results_1bit[0][0]}, max={results_1bit[-1][0]}, avg={sum(r[0] for r in results_1bit)/512:.1f}")
print(f"  Active words: min={results_1bit[0][1]}, max={max(r[1] for r in results_1bit)}")


# ================================================================
# TEST 2: Two-bit inputs — can XOR of two sparse expand. be sparser?
# If schedule is linear (XOR model): δW(a XOR b) = δW(a) XOR δW(b)
# Possible cancellation: two high-HW expansions may XOR to low-HW
# ================================================================

print()
print("=" * 70)
print("TEST 2: Two-bit inputs (XOR cancellation)")
print("=" * 70)

# Try all pairs of the 10 sparsest single-bit inputs
sparse_10 = results_1bit[:20]  # take 20 sparsest

best_2bit = []
for i in range(len(sparse_10)):
    for j in range(i+1, len(sparse_10)):
        _, _, w1, b1 = sparse_10[i]
        _, _, w2, b2 = sparse_10[j]
        M = [0]*16
        M[w1] ^= (1 << b1)
        M[w2] ^= (1 << b2)
        h = hw_schedule(M, xor_model=True)
        aw = active_words(M, xor_model=True)
        best_2bit.append((h, aw, w1, b1, w2, b2))

best_2bit.sort()

print(f"\n  Pairs from 20 sparsest single-bit inputs:")
print(f"\n  Top 10 LOWEST HW:")
print(f"  {'HW':>4} {'active':>6} | inputs")
for h, aw, w1, b1, w2, b2 in best_2bit[:10]:
    print(f"  {h:>4} {aw:>6} | W[{w1}]b{b1} + W[{w2}]b{b2}")

# Also: random 2-bit pairs
random.seed(42)
rand_2bit = []
for _ in range(10000):
    w1, b1 = random.randint(0,15), random.randint(0,31)
    w2, b2 = random.randint(0,15), random.randint(0,31)
    if (w1,b1) == (w2,b2): continue
    M = [0]*16
    M[w1] ^= (1 << b1)
    M[w2] ^= (1 << b2)
    h = hw_schedule(M, xor_model=True)
    aw = active_words(M, xor_model=True)
    rand_2bit.append((h, aw))

rand_2bit.sort()
print(f"\n  Random 2-bit pairs (N=10000):")
print(f"  min HW = {rand_2bit[0][0]}, min active = {rand_2bit[0][1]}")
print(f"  Top 5: {[r[0] for r in rand_2bit[:5]]}")


# ================================================================
# TEST 3: Random low-HW δM — search for minimum distance
# For HW(δM) = k, search for min HW(Σ(δM))
# ================================================================

print()
print("=" * 70)
print("TEST 3: Min distance search (random low-HW δM)")
print("=" * 70)

random.seed(123)

for input_hw in [1, 2, 3, 4, 8, 16, 32]:
    N_trials = 50000 if input_hw <= 4 else 20000

    min_hw = 9999
    min_aw = 99
    best_M = None

    for _ in range(N_trials):
        M = [0]*16
        # Set input_hw random bits
        positions = set()
        while len(positions) < input_hw:
            w = random.randint(0, 15)
            b = random.randint(0, 31)
            positions.add((w, b))
        for w, b in positions:
            M[w] ^= (1 << b)

        h = hw_schedule(M, xor_model=True)
        aw = active_words(M, xor_model=True)

        if h < min_hw:
            min_hw = h
            min_aw = aw
            best_M = list(M)

    print(f"  HW(δM)={input_hw:>2}: min HW(δW[16..63])={min_hw:>4}, min active_words={min_aw:>2}  (N={N_trials})")


# ================================================================
# TEST 4: REAL schedule (with carries) vs XOR model
# Does carry make it sparser or denser?
# ================================================================

print()
print("=" * 70)
print("TEST 4: Real schedule (ADD) vs XOR model")
print("=" * 70)

random.seed(456)
N4 = 5000

xor_hws = []
real_hws = []

for _ in range(N4):
    # Random single-bit δM
    w = random.randint(0, 15)
    b = random.randint(0, 31)
    M = [0]*16
    M[w] = 1 << b

    xor_hws.append(hw_schedule(M, xor_model=True))

    # For real: need base message + perturbed
    M_base = [random.randint(0, MASK32) for _ in range(16)]
    M_pert = list(M_base)
    M_pert[w] ^= (1 << b)

    W_base = expand_schedule(M_base)
    W_pert = expand_schedule(M_pert)
    real_hw = sum(hw(W_base[r] ^ W_pert[r]) for r in range(16, 64))
    real_hws.append(real_hw)

avg_xor = sum(xor_hws)/N4
avg_real = sum(real_hws)/N4
min_xor = min(xor_hws)
min_real = min(real_hws)

print(f"\n  N={N4} single-bit inputs")
print(f"  XOR model: avg HW={avg_xor:.1f}, min={min_xor}")
print(f"  Real (ADD): avg HW={avg_real:.1f}, min={min_real}")
print(f"  Difference: {avg_real - avg_xor:+.1f} (carry adds noise)")


# ================================================================
# TEST 5: Zero-word pattern — which δM gives most zero W[16..63]?
# Each zero word = one "free" round (no schedule difference).
# ================================================================

print()
print("=" * 70)
print("TEST 5: Maximum zero words in W[16..63]")
print("=" * 70)

# From TEST 1, we know single-bit min active = some number.
# Search more broadly.

random.seed(789)
best_zeros = 0
best_zeros_M = None

# Strategy 1: single-bit inputs
for w in range(16):
    for b in range(32):
        M = [0]*16; M[w] = 1 << b
        W = expand_schedule_xor(M)
        zeros = sum(1 for r in range(16, 64) if W[r] == 0)
        if zeros > best_zeros:
            best_zeros = zeros
            best_zeros_M = (f"W[{w}] bit {b}", list(M))

print(f"  Best single-bit: {best_zeros} zero words ({best_zeros_M[0]})")

# Strategy 2: random multi-bit with hill climbing
for trial in range(20):
    M = [0]*16
    # Start with random HW=3 input
    for _ in range(3):
        w = random.randint(0,15); b = random.randint(0,31)
        M[w] ^= (1<<b)

    W = expand_schedule_xor(M)
    zeros = sum(1 for r in range(16,64) if W[r]==0)

    # Hill climb: flip bits to increase zeros
    for step in range(500):
        w = random.randint(0,15); b = random.randint(0,31)
        M_try = list(M); M_try[w] ^= (1<<b)
        if M_try == [0]*16: continue
        W_try = expand_schedule_xor(M_try)
        zeros_try = sum(1 for r in range(16,64) if W_try[r]==0)
        if zeros_try > zeros:
            M = M_try; zeros = zeros_try

    if zeros > best_zeros:
        best_zeros = zeros
        input_hw = sum(hw(m) for m in M)
        best_zeros_M = (f"HW={input_hw}, HC trial {trial}", list(M))

print(f"  Best after HC: {best_zeros} zero words ({best_zeros_M[0]})")

# List the zero rounds for best
if best_zeros_M:
    M = best_zeros_M[1]
    W = expand_schedule_xor(M)
    zero_rounds = [r for r in range(16, 64) if W[r] == 0]
    nonzero_rounds = [r for r in range(16, 64) if W[r] != 0]
    print(f"  Zero rounds: {zero_rounds}")
    print(f"  Nonzero rounds: {len(nonzero_rounds)} out of 48")
    print(f"  Nonzero HW per word: {[hw(W[r]) for r in nonzero_rounds[:10]]}...")


# ================================================================
# TEST 6: Implications for Λ
# ================================================================

print()
print("=" * 70)
print("TEST 6: Implications for path cost Λ")
print("=" * 70)

# Best case: X zero words out of 48 in W[16..63]
# Active words: 48 - X
# Each active word contributes ~16 P-bits (from P2 results)
# But: word with low HW contributes fewer P-bits

# Lower bound on active words:
print(f"\n  Min active words in W[16..63] (from tests above):")
print(f"    Single-bit δM: min active = {results_1bit[0][1]}")
print(f"    Best found: {48 - best_zeros} active words")

min_active = 48 - best_zeros
lambda_schedule = min_active * 16  # 16 P-bits per active word (from P2)

print(f"\n  Λ lower bound from schedule alone:")
print(f"    Active words: {min_active}")
print(f"    P-bits per word: ~16 (from P2 experiment)")
print(f"    Λ_schedule ≥ {min_active} × 16 = {lambda_schedule}")
print(f"    Add Wang (16 free rounds): Λ_total ≥ {lambda_schedule}")
print(f"    Birthday: 2^128")
print(f"    Ratio: 2^{lambda_schedule}/2^128 = 2^{lambda_schedule - 128}")

if lambda_schedule > 128:
    print(f"\n  Schedule minimum distance PREVENTS path cheaper than birthday.")
    print(f"  Even the sparsest δM gives {min_active} active words = {lambda_schedule} bits.")
    print(f"  Need active ≤ 8 words for Λ ≤ 128. Have ≥ {min_active}.")
else:
    print(f"\n  ★ Schedule allows path with Λ ≤ 128!")
    print(f"  This is a potential attack direction!")
