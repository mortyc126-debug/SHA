#!/usr/bin/env python3
"""
Positional survival analysis.

Thermostat controls TOTAL HW → 128. But are all POSITIONS equal?
If some bit positions "survive" (stay different) longer than others,
the thermostat is homogeneous in quantity but INHOMOGENEOUS in position.

That would be a crack: specific bits carry information longer.

Experiment: 1-bit input diff, track WHICH bits of δstate are 1
at each round. Measure per-position survival time.
"""

import random
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

def sha_states(M):
    W=list(M)+[0]*(64-len(M))
    for i in range(16,64):W[i]=(ssig1(W[i-2])+W[i-7]+ssig0(W[i-15])+W[i-16])&MASK32
    a,b,c,d,e,f,g,h=IV
    ss=[(a,b,c,d,e,f,g,h)]
    for r in range(64):
        T1=(h+sig1(e)+ch(e,f,g)+K[r]+W[r])&MASK32;T2=(sig0(a)+maj(a,b,c))&MASK32
        h,g,f=g,f,e;e=(d+T1)&MASK32;d,c,b=c,b,a;a=(T1+T2)&MASK32
        ss.append((a,b,c,d,e,f,g,h))
    return ss

REG = ['a','b','c','d','e','f','g','h']

# ================================================================
# TEST 1: Per-position survival — P(bit (reg,k) is 1 in δstate[r])
# For 1-bit input diff δW[0] bit 0.
# ================================================================

print("=" * 70)
print("TEST 1: Per-position δ-bit frequency at each round")
print("=" * 70)

N = 3000
random.seed(42)

# bit_freq[r][reg][bit] = count of times δstate[r] has bit 1 at (reg, bit)
bit_freq = [[[ 0 for _ in range(32)] for _ in range(8)] for _ in range(65)]

for trial in range(N):
    M1 = [random.randint(0,MASK32) for _ in range(16)]
    M2 = list(M1); M2[0] ^= 1

    ss1 = sha_states(M1)
    ss2 = sha_states(M2)

    for r in range(65):
        for reg in range(8):
            diff = ss1[r][reg] ^ ss2[r][reg]
            for bit in range(32):
                if (diff >> bit) & 1:
                    bit_freq[r][reg][bit] += 1

# Convert to probability
bit_prob = [[[bit_freq[r][reg][bit]/N for bit in range(32)] for reg in range(8)] for r in range(65)]

# At equilibrium (r=32): all should be ≈ 0.50
# Deviation from 0.50 = positional bias

print(f"\n  N={N}, δM = W[0] bit 0")

# Round 1: where does the diff appear?
print(f"\n  Round 1 — initial diff pattern:")
for reg in range(8):
    nonzero = [(bit, bit_prob[1][reg][bit]) for bit in range(32) if bit_prob[1][reg][bit] > 0.01]
    if nonzero:
        top = sorted(nonzero, key=lambda x: -x[1])[:5]
        print(f"    {REG[reg]}: {', '.join(f'b{b}={p:.3f}' for b,p in top)}")

# Track specific positions through rounds
print(f"\n  Tracking initial diff bits through rounds:")
# At round 1, which positions have P > 0.5? Track them.
initial_hot = []
for reg in range(8):
    for bit in range(32):
        if bit_prob[1][reg][bit] > 0.1:
            initial_hot.append((reg, bit, bit_prob[1][reg][bit]))

print(f"  {len(initial_hot)} initially hot positions (P > 0.1 at r=1)")

# Track each hot position
print(f"\n  {'position':>10} | r=1    r=2    r=3    r=4    r=5    r=8    r=16   r=32   r=64")
for reg, bit, p0 in sorted(initial_hot, key=lambda x: -x[2])[:15]:
    row = f"  {REG[reg]}[{bit:>2}]    |"
    for r in [1, 2, 3, 4, 5, 8, 16, 32, 64]:
        p = bit_prob[r][reg][bit]
        if p > 0.55:
            row += f" {p:.3f}★"
        elif p < 0.45:
            row += f" {p:.3f}↓"
        else:
            row += f" {p:.3f} "
    print(row)


# ================================================================
# TEST 2: Survival time per position
# Define: survival(reg, bit) = last round r where P(diff=1) > 0.55
# (significantly above 0.50)
# ================================================================

print()
print("=" * 70)
print("TEST 2: Survival time per position")
print("=" * 70)

# For EVERY position: find last round with P > threshold
threshold = 0.55
survival = {}

for reg in range(8):
    for bit in range(32):
        last_r = 0
        for r in range(1, 65):
            if bit_prob[r][reg][bit] > threshold:
                last_r = r
        survival[(reg, bit)] = last_r

# Distribution of survival times
surv_dist = defaultdict(int)
for (reg, bit), s in survival.items():
    surv_dist[s] += 1

print(f"\n  Threshold: P > {threshold}")
print(f"\n  Survival time distribution:")
for s in sorted(surv_dist.keys()):
    if surv_dist[s] > 0:
        bar = "█" * surv_dist[s]
        print(f"    r={s:>2}: {surv_dist[s]:>3} positions {bar}")

# Longest surviving positions
top_survivors = sorted(survival.items(), key=lambda x: -x[1])
print(f"\n  Top 20 longest surviving positions:")
for (reg, bit), s in top_survivors[:20]:
    print(f"    {REG[reg]}[{bit:>2}]: survives until round {s}, P(r={s})={bit_prob[s][reg][bit]:.3f}")

# Average survival
avg_surv = sum(survival.values()) / len(survival)
max_surv = max(survival.values())
print(f"\n  Average survival = {avg_surv:.1f} rounds")
print(f"  Maximum survival = {max_surv} rounds")
print(f"  Expected if uniform: ~3-4 rounds (shift register length)")


# ================================================================
# TEST 3: Is survival related to register or bit position?
# ================================================================

print()
print("=" * 70)
print("TEST 3: Survival by register and bit position")
print("=" * 70)

# By register
print(f"\n  Average survival by register:")
for reg in range(8):
    avg = sum(survival[(reg, bit)] for bit in range(32)) / 32
    mx = max(survival[(reg, bit)] for bit in range(32))
    print(f"    {REG[reg]}: avg={avg:.1f}, max={mx}")

# By bit position
print(f"\n  Average survival by bit position:")
bit_avg = []
for bit in range(32):
    avg = sum(survival[(reg, bit)] for reg in range(8)) / 8
    bit_avg.append((bit, avg))

bit_avg.sort(key=lambda x: -x[1])
print(f"    Top 5 (longest surviving bits):")
for bit, avg in bit_avg[:5]:
    print(f"      bit {bit:>2}: avg survival = {avg:.1f}")
print(f"    Bottom 5:")
for bit, avg in bit_avg[-5:]:
    print(f"      bit {bit:>2}: avg survival = {avg:.1f}")


# ================================================================
# TEST 4: CRITICAL — do surviving positions CORRELATE with output?
# If a bit survives to round r and correlates with H → useful!
# ================================================================

print()
print("=" * 70)
print("TEST 4: Do surviving bits correlate with output H?")
print("=" * 70)

# For each position that survives > 4 rounds:
# Check: P(H[j][k] = 1 | δstate[r][reg][bit] = 1) vs P(H[j][k] = 1)

# Collect: for pairs with 1-bit diff, record δstate at round 5
# and δH (output diff)
N4 = 5000
random.seed(123)

# We want: δstate[r] bits that predict δH bits.
# If such predictors exist at r > 4 → info survives longer than expected.

d_state_5 = []  # δstate at round 5 for each trial
d_hash = []     # δH for each trial

for trial in range(N4):
    M1 = [random.randint(0,MASK32) for _ in range(16)]
    M2 = list(M1); M2[0] ^= 1

    ss1 = sha_states(M1)
    ss2 = sha_states(M2)

    ds5 = tuple(ss1[5][i] ^ ss2[5][i] for i in range(8))
    d_state_5.append(ds5)

    # δH via feedforward
    dH = tuple((ss1[64][i]+IV[i])&MASK32 ^ (ss2[64][i]+IV[i])&MASK32 for i in range(8))
    # Fix: compute H properly
    H1 = tuple((ss1[64][i]+IV[i])&MASK32 for i in range(8))
    H2 = tuple((ss2[64][i]+IV[i])&MASK32 for i in range(8))
    dH = tuple(H1[i] ^ H2[i] for i in range(8))
    d_hash.append(dH)

# For top surviving positions at round 5: check correlation with δH
print(f"\n  N={N4}, checking δstate[5] → δH correlation")

# Pick positions that survive past round 4
long_survivors = [(reg, bit) for (reg, bit), s in survival.items() if s >= 4]
print(f"  {len(long_survivors)} positions survive to round 4+")

# For each surviving position: compute max correlation with any δH bit
best_corrs = []

for reg, bit in long_survivors[:50]:  # check top 50
    # δstate[5][reg][bit] for each trial
    ds_bits = [(d_state_5[t][reg] >> bit) & 1 for t in range(N4)]

    # Check correlation with each δH bit
    max_corr = 0
    max_target = (0, 0)

    for h_reg in range(8):
        for h_bit in range(32):
            dh_bits = [(d_hash[t][h_reg] >> h_bit) & 1 for t in range(N4)]

            # Correlation = P(agree) - 0.5
            agree = sum(1 for a, b in zip(ds_bits, dh_bits) if a == b) / N4
            corr = abs(agree - 0.5)

            if corr > max_corr:
                max_corr = corr
                max_target = (h_reg, h_bit)

    best_corrs.append((max_corr, reg, bit, max_target))

best_corrs.sort(reverse=True)

print(f"\n  Top 10 correlations (δstate[5] bit → δH bit):")
print(f"  {'source':>10} | {'target':>10} | {'corr':>6} | note")

for corr, reg, bit, (h_reg, h_bit) in best_corrs[:10]:
    note = ""
    # Expected max for N=5000, 256 targets: ~0.02 (noise)
    expected_max = 2.5 / (N4**0.5)  # ~3σ for 256 tests
    if corr > expected_max * 2:
        note = "★ ABOVE NOISE"
    else:
        note = "noise"
    print(f"  {REG[reg]}[{bit:>2}]   | H[{h_reg}][{h_bit:>2}] | {corr:.4f} | {note}")

print(f"\n  Expected noise level: ~{expected_max:.4f} (3σ for {N4} samples)")


# ================================================================
# TEST 5: Same analysis for DIFFERENT input positions
# Does the choice of which W[0] bit to flip affect survival?
# ================================================================

print()
print("=" * 70)
print("TEST 5: Survival depends on INPUT bit position?")
print("=" * 70)

N5 = 1000

for input_bit in [0, 7, 15, 16, 25, 31]:
    random.seed(input_bit * 100 + 42)
    local_surv = defaultdict(int)

    for trial in range(N5):
        M1 = [random.randint(0,MASK32) for _ in range(16)]
        M2 = list(M1); M2[0] ^= (1 << input_bit)

        ss1 = sha_states(M1)
        ss2 = sha_states(M2)

        # For each position: is it diff at round 5?
        for reg in range(8):
            diff = ss1[5][reg] ^ ss2[5][reg]
            for bit in range(32):
                if (diff >> bit) & 1:
                    local_surv[(reg, bit)] += 1

    # Which positions have highest freq at round 5?
    top5 = sorted(local_surv.items(), key=lambda x: -x[1])[:3]
    avg_count = sum(local_surv.values()) / max(len(local_surv), 1)

    print(f"  δW[0] bit {input_bit:>2}: {len(local_surv)} active positions at r=5, "
          f"top={', '.join(f'{REG[r]}[{b}]={c/N5:.2f}' for (r,b),c in top5)}")


# ================================================================
# SYNTHESIS
# ================================================================

print()
print("=" * 70)
print("SYNTHESIS")
print("=" * 70)

print(f"""
  Thermostat controls TOTAL HW → 128.
  Question: is it HOMOGENEOUS across positions?

  Results:
  - Max survival at threshold P>0.55: {max_surv} rounds
  - Average survival: {avg_surv:.1f} rounds
  - Expected from shift register: 3-4 rounds

  If max_survival >> 4: POSITIONAL INHOMOGENEITY (crack!)
  If max_survival ≈ 4: HOMOGENEOUS thermostat (no crack)

  Correlation of surviving bits with output:
  - Check results above for ★ ABOVE NOISE signals
""")
