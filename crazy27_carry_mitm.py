#!/usr/bin/env python3
"""
crazy27_carry_mitm.py — Meet-in-the-middle on carry chain for De17(W[14])

Key insight: In modular addition a+b, carry at bit k depends ONLY on bits 0..k-1.
Split W[14] into LOW (bits 0-15) and HIGH (bits 16-31).
Carries propagate LOW→HIGH but not HIGH→LOW.
This creates a natural MITM split reducing 2^32 to ~2^17.
"""

import random
import time
from collections import defaultdict

# --- Standard SHA-256 primitives ---
MASK = 0xFFFFFFFF
K = [0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
     0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
     0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc]
H0 = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&MASK
def Sig0(x): return rotr(x,2)^rotr(x,13)^rotr(x,22)
def Sig1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def sig0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def sig1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def Ch(e,f,g): return ((e&f)^(~e&g))&MASK
def Maj(a,b,c): return ((a&b)^(a&c)^(b&c))&MASK
def add32(*a):
    s=0
    for x in a: s=(s+x)&MASK
    return s
def hw(x): return bin(x&MASK).count('1')

def sha_round(st,w,k):
    a,b,c,d,e,f,g,h=st
    T1=add32(h,Sig1(e),Ch(e,f,g),k,w)
    T2=add32(Sig0(a),Maj(a,b,c))
    return [add32(T1,T2),a,b,c,add32(d,T1),e,f,g]

def expand_W(W16):
    W=list(W16)
    for i in range(16,20):
        W.append(add32(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    return W

def get_de17(msg, iv=None):
    if iv is None: iv = list(H0)
    W=list(msg); Wp=list(W); Wp[0]^=0x80000000
    s=list(iv); sp=list(iv)
    s=sha_round(s,W[0],K[0]); sp=sha_round(sp,Wp[0],K[0])
    for t in range(1,16):
        a,b,c,d,e,f,g,h=s
        a2,b2,c2,d2,e2,f2,g2,h2=sp
        tp=add32(h,Sig1(e),Ch(e,f,g),K[t])
        tp2=add32(h2,Sig1(e2),Ch(e2,f2,g2),K[t])
        target=add32(d,tp,W[t])
        Wp[t]=(target-d2-tp2)&MASK
        s=sha_round(s,W[t],K[t]); sp=sha_round(sp,Wp[t],K[t])
    We=expand_W(msg); Wpe=expand_W(Wp)
    s1=list(iv); s2=list(iv)
    for t in range(17):
        s1=sha_round(s1,We[t],K[t])
        s2=sha_round(s2,Wpe[t],K[t])
    return s1[4]^s2[4]

# ============================================================
# SETUP: Base message from seed 0xC027
# ============================================================
print("="*70)
print("CARRY-CHAIN MEET-IN-THE-MIDDLE FOR De17(W[14])")
print("="*70)

rng = random.Random(0xC027)
base_msg = [rng.randint(0, MASK) for _ in range(16)]
print(f"\nBase message (seed 0xC027):")
for i in range(0, 16, 4):
    print(f"  W[{i:2d}..{i+3:2d}] = {' '.join(f'{base_msg[j]:08x}' for j in range(i,i+4))}")

base_de17 = get_de17(base_msg)
print(f"\nBase De17 = 0x{base_de17:08x}  (hw={hw(base_de17)})")

# ============================================================
# PHASE 1: Understand carry coupling structure
# Sample 10000 random W[14] values, measure how LOW affects HIGH of De17
# ============================================================
print("\n" + "="*70)
print("PHASE 1: CARRY COUPLING ANALYSIS (10000 random W[14])")
print("="*70)

t0 = time.time()
N_SAMPLE = 10000
samples = []
rng2 = random.Random(42)

for _ in range(N_SAMPLE):
    w14 = rng2.randint(0, MASK)
    msg = list(base_msg)
    msg[14] = w14
    de17 = get_de17(msg)
    samples.append((w14, de17))

dt = time.time() - t0
print(f"Computed {N_SAMPLE} De17 values in {dt:.2f}s")

# Analyze: for same HIGH half of W[14], how much does LOW affect De17_high?
# Group by HIGH half
high_groups = defaultdict(list)
for w14, de17 in samples:
    high_half = (w14 >> 16) & 0xFFFF
    low_half = w14 & 0xFFFF
    de17_high = (de17 >> 16) & 0xFFFF
    de17_low = de17 & 0xFFFF
    high_groups[high_half].append((low_half, de17_low, de17_high))

# Find groups with multiple entries
multi_groups = {k: v for k, v in high_groups.items() if len(v) >= 2}
print(f"\nHIGH groups with 2+ entries: {len(multi_groups)}")

if multi_groups:
    # For groups with same HIGH(W[14]), measure variation in De17_high
    carry_leak_bits = []
    for hval, entries in list(multi_groups.items())[:100]:
        de17_highs = set(dh for _, _, dh in entries)
        # How many distinct De17_high values?
        variation = len(de17_highs)
        if variation > 1:
            # XOR all pairs to see which bits differ
            xor_acc = 0
            vals = list(de17_highs)
            for i in range(len(vals)):
                for j in range(i+1, len(vals)):
                    xor_acc |= vals[i] ^ vals[j]
            carry_leak_bits.append(hw(xor_acc))

    if carry_leak_bits:
        avg_leak = sum(carry_leak_bits) / len(carry_leak_bits)
        print(f"Avg carry leakage into De17_high: {avg_leak:.1f} bits (from {len(carry_leak_bits)} groups)")
    else:
        print("Not enough variation to measure carry leakage in this sample")

# Also measure: correlation between LOW(W14) change and De17 bits
print("\n--- Bit-by-bit sensitivity analysis ---")
# For each bit of W[14], flip it and see which De17 bits change
bit_affects = [0] * 32  # bit_affects[k] = OR of all De17 changes when bit k of W14 flips
msg_test = list(base_msg)
de17_ref = get_de17(msg_test)
for bit in range(32):
    msg_test2 = list(base_msg)
    msg_test2[14] ^= (1 << bit)
    de17_flip = get_de17(msg_test2)
    diff = de17_ref ^ de17_flip
    bit_affects[bit] = diff

print("W[14] bit → De17 bits affected (hw of diff):")
print("  LOW bits (0-15) → De17 change:")
for b in range(16):
    d = bit_affects[b]
    low_bits = hw(d & 0xFFFF)
    high_bits = hw(d >> 16)
    print(f"    bit {b:2d}: De17_low={low_bits:2d} bits, De17_high={high_bits:2d} bits  (total hw={hw(d):2d})")
print("  HIGH bits (16-31) → De17 change:")
for b in range(16, 32):
    d = bit_affects[b]
    low_bits = hw(d & 0xFFFF)
    high_bits = hw(d >> 16)
    print(f"    bit {b:2d}: De17_low={low_bits:2d} bits, De17_high={high_bits:2d} bits  (total hw={hw(d):2d})")

# Key question: do LOW bits of W[14] affect HIGH bits of De17?
low_to_high_leakage = 0
for b in range(16):
    low_to_high_leakage |= (bit_affects[b] >> 16)
print(f"\nCombined LOW→HIGH leakage mask: 0x{low_to_high_leakage:04x} (hw={hw(low_to_high_leakage)})")

high_to_low_leakage = 0
for b in range(16, 32):
    high_to_low_leakage |= (bit_affects[b] & 0xFFFF)
print(f"Combined HIGH→LOW leakage mask: 0x{high_to_low_leakage:04x} (hw={hw(high_to_low_leakage)})")

# ============================================================
# PHASE 2: MITM IMPLEMENTATION
# ============================================================
print("\n" + "="*70)
print("PHASE 2: MEET-IN-THE-MIDDLE (2^16 + 2^16 = 2^17 evaluations)")
print("="*70)

base_w14 = base_msg[14]
base_low = base_w14 & 0xFFFF
base_high = (base_w14 >> 16) & 0xFFFF

print(f"\nBase W[14] = 0x{base_w14:08x}  (LOW=0x{base_low:04x}, HIGH=0x{base_high:04x})")
print(f"Base De17  = 0x{base_de17:08x}")

# --- Forward phase: enumerate all 2^16 LOW values (HIGH fixed to base) ---
print("\n[Forward] Enumerating 2^16 LOW values...")
t0 = time.time()
low_table = {}  # de17 → list of LOW values
low_results = []  # (low_val, de17)

for low_val in range(0x10000):
    msg = list(base_msg)
    msg[14] = (base_high << 16) | low_val
    de17 = get_de17(msg)
    low_results.append((low_val, de17))
    if de17 not in low_table:
        low_table[de17] = []
    low_table[de17].append(low_val)

dt_low = time.time() - t0
print(f"  Done in {dt_low:.1f}s. Unique De17 values: {len(low_table)}")

# Check for De17=0 in forward phase alone
if 0 in low_table:
    print(f"  *** De17=0 found in LOW sweep! {len(low_table[0])} solutions ***")
    for lv in low_table[0][:5]:
        w14 = (base_high << 16) | lv
        print(f"      W[14]=0x{w14:08x}")

# --- Backward phase: enumerate all 2^16 HIGH values (LOW fixed to base) ---
print("\n[Backward] Enumerating 2^16 HIGH values...")
t0 = time.time()
high_table = {}  # de17 → list of HIGH values
high_results = []

for high_val in range(0x10000):
    msg = list(base_msg)
    msg[14] = (high_val << 16) | base_low
    de17 = get_de17(msg)
    high_results.append((high_val, de17))
    if de17 not in high_table:
        high_table[de17] = []
    high_table[de17].append(high_val)

dt_high = time.time() - t0
print(f"  Done in {dt_high:.1f}s. Unique De17 values: {len(high_table)}")

if 0 in high_table:
    print(f"  *** De17=0 found in HIGH sweep! {len(high_table[0])} solutions ***")
    for hv in high_table[0][:5]:
        w14 = (hv << 16) | base_low
        print(f"      W[14]=0x{w14:08x}")

# --- MITM Matching ---
# Naive MITM: De17 is NOT linear, so we can't just XOR.
# But we can try: for combined W[14] = (HIGH << 16) | LOW,
# De17(combined) ≈ De17(LOW_only) ⊕ De17(HIGH_only) ⊕ De17(base)
# if carries didn't couple. Let's test this approximation.

print("\n[Matching] Attempting MITM match...")
print("  Strategy 1: Linear approximation De17(L,H) ≈ De17(L,H0) XOR De17(L0,H) XOR De17(L0,H0)")

t0 = time.time()
mitm_candidates = []

# For each LOW value, we need De17(LOW, base_HIGH) XOR base_de17
# to match De17(base_LOW, HIGH) for some HIGH
# Under linear approx: De17(LOW, HIGH) = De17(LOW, H0) XOR De17(L0, HIGH) XOR De17(L0, H0)
# Want De17(LOW, HIGH) = 0
# So need: De17(LOW, H0) XOR De17(L0, HIGH) = De17(L0, H0) = base_de17

# Build lookup: for each LOW, compute target = De17(LOW, H0) XOR base_de17
# Then look up target in high_table

match_count = 0
for low_val, de17_low in low_results:
    target = de17_low ^ base_de17
    if target in high_table:
        for high_val in high_table[target]:
            mitm_candidates.append((low_val, high_val))
            match_count += 1

dt_match = time.time() - t0
print(f"  Found {match_count} MITM candidates in {dt_match:.3f}s")

# --- Verify MITM candidates ---
print(f"\n[Verification] Checking {min(match_count, 50000)} MITM candidates...")
t0 = time.time()
mitm_exact_solutions = []
mitm_near_solutions = []  # hw(De17) <= 4

checked = 0
for low_val, high_val in mitm_candidates[:50000]:
    w14 = (high_val << 16) | low_val
    msg = list(base_msg)
    msg[14] = w14
    de17 = get_de17(msg)
    checked += 1
    if de17 == 0:
        mitm_exact_solutions.append(w14)
    elif hw(de17) <= 4:
        mitm_near_solutions.append((w14, de17))

dt_verify = time.time() - t0
print(f"  Checked {checked} candidates in {dt_verify:.1f}s")
print(f"  Exact De17=0: {len(mitm_exact_solutions)}")
print(f"  Near (hw<=4):  {len(mitm_near_solutions)}")

if mitm_exact_solutions:
    print("\n  MITM exact solutions:")
    for w14 in mitm_exact_solutions[:10]:
        print(f"    W[14] = 0x{w14:08x}")

# Distribution of De17 hamming weight in candidates
hw_dist = defaultdict(int)
for low_val, high_val in mitm_candidates[:50000]:
    w14 = (high_val << 16) | low_val
    msg = list(base_msg)
    msg[14] = w14
    de17 = get_de17(msg)
    hw_dist[hw(de17)] += 1

print("\n  Hamming weight distribution of De17 for MITM candidates:")
for h in sorted(hw_dist.keys())[:20]:
    bar = "#" * min(hw_dist[h], 60)
    print(f"    hw={h:2d}: {hw_dist[h]:6d}  {bar}")

# ============================================================
# PHASE 3: Strategy 2 — Carry-aware MITM
# ============================================================
print("\n" + "="*70)
print("PHASE 3: CARRY-AWARE MITM")
print("="*70)
print("Since De17 is nonlinear (carry coupling), we try a refined approach:")
print("For each LOW, record De17_low_bits (bits 0-15) AND carry signature.")
print("The 'carry signature' = De17_high XOR De17_high_from_base_low")

# For carry-aware matching, we index by De17_low (16 bits) plus carry info
# We know: changing LOW affects all 32 bits of De17 (due to nonlinearity)
# But the coupling might be structured.

# Let's measure: how well does De17_low predict De17 for the combined value?
# For the combined W[14], De17_low_bits are determined mainly by LOW half.

# Actually, let's try a different MITM: match on De17 directly with tolerance.
# Use the LOW sweep to build a table indexed by De17 value.
# Then for each HIGH value, compute what De17 "correction" is needed.

# Better approach: partial-bit matching
# Match on bits 0-15 of De17 first (controlled mainly by LOW),
# then check bits 16-31.

print("\n  Strategy 2: Match on De17_low (bits 0-15), verify De17_high")

t0 = time.time()
# Index LOW results by De17 lower 16 bits
low_by_de17low = defaultdict(list)
for low_val, de17 in low_results:
    de17_low_bits = de17 & 0xFFFF
    low_by_de17low[de17_low_bits].append((low_val, de17))

# For each HIGH, compute what De17_low we need: want combined De17_low = 0
# Under approximation: De17_low(combined) ≈ De17_low(L,H0) XOR De17_low(L0,H) XOR De17_low(base)
strategy2_candidates = []
for high_val, de17_high_sweep in high_results:
    target_low = (de17_high_sweep ^ base_de17) & 0xFFFF
    if target_low in low_by_de17low:
        for low_val, de17_low_sweep in low_by_de17low[target_low]:
            strategy2_candidates.append((low_val, high_val))

print(f"  Strategy 2 candidates: {len(strategy2_candidates)}")

# Verify
s2_exact = []
s2_near = []
for low_val, high_val in strategy2_candidates[:50000]:
    w14 = (high_val << 16) | low_val
    msg = list(base_msg)
    msg[14] = w14
    de17 = get_de17(msg)
    if de17 == 0:
        s2_exact.append(w14)
    elif hw(de17) <= 4:
        s2_near.append((w14, de17))

dt_s2 = time.time() - t0
print(f"  Checked in {dt_s2:.1f}s")
print(f"  Exact De17=0: {len(s2_exact)}")
print(f"  Near (hw<=4):  {len(s2_near)}")

if s2_exact:
    print("  Solutions:")
    for w14 in s2_exact[:10]:
        print(f"    W[14] = 0x{w14:08x}")

# ============================================================
# PHASE 4: BRUTE-FORCE VERIFICATION (2^20 random)
# ============================================================
print("\n" + "="*70)
print("PHASE 4: BRUTE-FORCE VERIFICATION (2^20 = 1048576 random W[14])")
print("="*70)

t0 = time.time()
rng3 = random.Random(0xBF)
bf_solutions = []
bf_near = []
N_BF = 1 << 20
hw_bf_dist = defaultdict(int)

for i in range(N_BF):
    w14 = rng3.randint(0, MASK)
    msg = list(base_msg)
    msg[14] = w14
    de17 = get_de17(msg)
    h = hw(de17)
    hw_bf_dist[h] += 1
    if de17 == 0:
        bf_solutions.append(w14)
    elif h <= 4:
        bf_near.append((w14, de17, h))

dt_bf = time.time() - t0
print(f"Brute-force done in {dt_bf:.1f}s")
print(f"De17=0 solutions: {len(bf_solutions)}")
print(f"Near (hw<=4):     {len(bf_near)}")
print(f"Empirical P(De17=0) ≈ {len(bf_solutions)/N_BF:.2e}  (expected ~2^-32 = {2**-32:.2e})")
if bf_solutions:
    print(f"  → That's ~2^{-32 + len(bf_solutions).bit_length():.1f} ... way more than expected!")
    for w14 in bf_solutions[:5]:
        print(f"    W[14] = 0x{w14:08x}")

print("\nHamming weight distribution (brute-force):")
total = 0
cumulative = 0
for h in range(33):
    c = hw_bf_dist.get(h, 0)
    cumulative += c
    if c > 0:
        expected_frac = 1.0  # just show raw
        bar_len = int(c / max(hw_bf_dist.values()) * 50)
        bar = "#" * bar_len
        print(f"  hw={h:2d}: {c:7d} ({100*c/N_BF:6.3f}%)  {bar}")

peak_hw = max(hw_bf_dist, key=hw_bf_dist.get)
mean_hw = sum(h * c for h, c in hw_bf_dist.items()) / N_BF
print(f"\n  Peak at hw={peak_hw}, mean hw={mean_hw:.2f}")
print(f"  (Random 32-bit: expected mean=16.0, peak=16)")

# ============================================================
# PHASE 5: COMBINED ANALYSIS
# ============================================================
print("\n" + "="*70)
print("PHASE 5: COMBINED ANALYSIS & CARRY COUPLING FACTOR")
print("="*70)

# All unique De17 values from LOW sweep
low_de17_set = set(de17 for _, de17 in low_results)
high_de17_set = set(de17 for _, de17 in high_results)
print(f"\nUnique De17 from LOW sweep:  {len(low_de17_set)} / 65536")
print(f"Unique De17 from HIGH sweep: {len(high_de17_set)} / 65536")

# How many De17_low (bits 0-15) values are covered?
low_de17_low_set = set(de17 & 0xFFFF for _, de17 in low_results)
low_de17_high_set = set((de17 >> 16) & 0xFFFF for _, de17 in low_results)
high_de17_low_set = set(de17 & 0xFFFF for _, de17 in high_results)
high_de17_high_set = set((de17 >> 16) & 0xFFFF for _, de17 in high_results)

print(f"\nDe17 half-coverage:")
print(f"  LOW sweep  → De17_low:  {len(low_de17_low_set):5d}/65536  De17_high: {len(low_de17_high_set):5d}/65536")
print(f"  HIGH sweep → De17_low:  {len(high_de17_low_set):5d}/65536  De17_high: {len(high_de17_high_set):5d}/65536")

# Carry coupling factor estimation
# If De17 were linear in W[14], then MITM would find all solutions.
# The "coupling factor" C means MITM misses solutions by a factor of ~2^C.
# We can estimate C from the MITM match rate vs expected.

total_mitm_solutions = len(mitm_exact_solutions) + len(s2_exact)
all_mitm_solutions = set(mitm_exact_solutions) | set(s2_exact)
print(f"\nTotal unique MITM solutions (both strategies): {len(all_mitm_solutions)}")

# What fraction of the 2^32 space did each phase search?
print(f"\nSearch effort:")
print(f"  LOW sweep:    2^16 = {1<<16}")
print(f"  HIGH sweep:   2^16 = {1<<16}")
print(f"  MITM matching:  {len(mitm_candidates)} + {len(strategy2_candidates)} candidates checked")
print(f"  Brute-force:  2^20 = {N_BF}")
print(f"  Total De17 evals: ~{(1<<16)*2 + min(len(mitm_candidates),50000) + min(len(strategy2_candidates),50000) + N_BF}")

# Linearity test: is De17(L1,H1) XOR De17(L1,H0) XOR De17(L0,H1) XOR De17(L0,H0) = 0?
print("\n--- Linearity Test ---")
print("Testing: De17(L,H) XOR De17(L,H0) XOR De17(L0,H) XOR De17(L0,H0) = 0 ?")
print("(If zero, MITM is exact; nonzero = carry coupling residual)")

rng4 = random.Random(123)
residuals = []
for _ in range(1000):
    l = rng4.randint(0, 0xFFFF)
    h = rng4.randint(0, 0xFFFF)
    msg_lh = list(base_msg); msg_lh[14] = (h << 16) | l
    msg_lh0 = list(base_msg); msg_lh0[14] = (base_high << 16) | l
    msg_l0h = list(base_msg); msg_l0h[14] = (h << 16) | base_low
    de17_lh = get_de17(msg_lh)
    de17_lh0 = get_de17(msg_lh0)
    de17_l0h = get_de17(msg_l0h)
    residual = de17_lh ^ de17_lh0 ^ de17_l0h ^ base_de17
    residuals.append(residual)

nonzero = sum(1 for r in residuals if r != 0)
mean_hw_residual = sum(hw(r) for r in residuals) / len(residuals)
print(f"  Out of 1000 random (L,H) pairs:")
print(f"  Nonzero residuals: {nonzero}/1000 ({100*nonzero/1000:.1f}%)")
print(f"  Mean hw of residual: {mean_hw_residual:.2f}")
print(f"  → Carry coupling factor ≈ {mean_hw_residual:.1f} bits of uncertainty")

if mean_hw_residual < 1:
    print("  *** Nearly linear! MITM should work very well. ***")
elif mean_hw_residual < 8:
    print(f"  Moderate coupling. MITM reduces search by ~2^{{32-17-{mean_hw_residual:.0f}}} = 2^{32-17-mean_hw_residual:.0f}")
else:
    print(f"  Heavy coupling ({mean_hw_residual:.1f} bits). MITM provides limited speedup.")
    print(f"  Effective reduction: 2^32 → 2^{17+mean_hw_residual:.0f}")

# Residual bit distribution
residual_or = 0
for r in residuals:
    residual_or |= r
print(f"\n  OR of all residuals: 0x{residual_or:08x} (hw={hw(residual_or)})")
print(f"  → {32 - hw(residual_or)} bits are perfectly linear (carry-free)")
print(f"  → {hw(residual_or)} bits have carry coupling")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "="*70)
print("SUMMARY")
print("="*70)
print(f"""
Base message seed: 0xC027
Target: De17(W[14]) = 0  (32-bit XOR of SHA-256 register e after round 17)

De17 statistics (from 2^20 brute-force):
  Mean Hamming weight: {mean_hw:.2f} (random: 16.0)
  P(De17=0): {len(bf_solutions)}/{N_BF} = {len(bf_solutions)/N_BF:.2e}
  Brute-force solutions: {len(bf_solutions)}

MITM results (2^17 evaluations):
  Strategy 1 (full-word XOR match): {len(mitm_exact_solutions)} exact solutions
  Strategy 2 (half-word match):     {len(s2_exact)} exact solutions
  Total unique MITM solutions:      {len(all_mitm_solutions)}
  MITM candidates generated:        {len(mitm_candidates)} (S1) + {len(strategy2_candidates)} (S2)

Carry coupling analysis:
  Linearity residual: {100*nonzero/1000:.1f}% nonzero
  Mean coupling: {mean_hw_residual:.2f} bits
  Linear bits: {32 - hw(residual_or)}/32
  Coupled bits: {hw(residual_or)}/32

Conclusion: {"MITM is effective! Carry coupling is low." if mean_hw_residual < 4 else "Carry coupling is significant. MITM provides partial speedup but not full 2^15 reduction."}
  Effective search reduction: 2^32 → ~2^{17 + mean_hw_residual:.0f}
  Speedup factor: ~2^{max(0, 32 - 17 - mean_hw_residual):.0f}x over brute-force
""")
