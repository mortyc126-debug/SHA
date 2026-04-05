#!/usr/bin/env python3
"""
Amplification law: is 1.12× per round EXACT?

If amplification α is constant: HW(δstate[r]) = HW(δstate[0]) × α^r
After 64 rounds: HW = α^64

If α = 1.12: α^64 = 1.12^64 = 1556. But HW caps at 256.
So: amplification MUST saturate. Where? At what round?

Key questions:
1. Is α constant across rounds?
2. Does α depend on current HW(δstate)?
3. What is the exact amplification curve HW(δstate[r])?
4. Where does saturation happen?
5. Is there a FORMULA connecting input δ and output δ?
"""

import random, math

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

def state_hw_diff(s1, s2):
    return sum(hw(a^b) for a,b in zip(s1, s2))


# ================================================================
# TEST 1: Exact amplification curve — HW(δstate[r]) vs r
# For 1-bit input difference (δW[0] = 1, rest same)
# ================================================================

print("=" * 70)
print("TEST 1: Amplification curve HW(δstate[r]) for 1-bit δM")
print("=" * 70)

N = 3000
random.seed(42)

hw_by_round = [[] for _ in range(65)]

for _ in range(N):
    M1 = [random.randint(0,MASK32) for _ in range(16)]
    M2 = list(M1); M2[0] ^= 1  # 1-bit difference

    ss1 = sha_states(M1)
    ss2 = sha_states(M2)

    for r in range(65):
        hw_by_round[r].append(state_hw_diff(ss1[r], ss2[r]))

print(f"\n  N={N}, δM = flip W[0] bit 0")
print(f"\n  {'r':>3} | {'E[HW]':>7} | {'std':>5} | {'α(r)':>6} | {'α_cum':>7} | model")

prev_hw = 0
cum_alpha = 1.0

for r in range(65):
    avg = sum(hw_by_round[r]) / N
    std = (sum((x-avg)**2 for x in hw_by_round[r])/N)**0.5

    if r > 0 and prev_hw > 0.5:
        alpha = avg / prev_hw
        cum_alpha *= alpha
    else:
        alpha = 0
        cum_alpha = avg if avg > 0 else 1

    # Theoretical model: logistic growth
    # HW(r) = 256 / (1 + (256/HW0 - 1) * exp(-k*r))
    # At saturation HW → 128 (half of 256)

    if r <= 4 or r % 4 == 0 or r >= 60:
        prev_hw_safe = prev_hw if prev_hw > 0 else 1
        print(f"  {r:>3} | {avg:>7.1f} | {std:>5.1f} | {alpha:>6.3f} | {cum_alpha:>7.1f} | ")

    prev_hw = avg

# ================================================================
# TEST 2: Per-round amplification α(r) — is it CONSTANT?
# α(r) = E[HW(δstate[r+1])] / E[HW(δstate[r])]
# ================================================================

print()
print("=" * 70)
print("TEST 2: Per-round amplification α(r)")
print("=" * 70)

alphas = []
for r in range(1, 64):
    avg_r = sum(hw_by_round[r]) / N
    avg_r1 = sum(hw_by_round[r+1]) / N
    if avg_r > 0.5:
        alpha = avg_r1 / avg_r
        alphas.append((r, alpha, avg_r, avg_r1))

print(f"\n  {'r':>3} | {'α = HW[r+1]/HW[r]':>18} | {'HW[r]':>6} → {'HW[r+1]':>7}")

for r, alpha, hw_r, hw_r1 in alphas:
    flag = ""
    if alpha > 1.5: flag = " ★ GROWTH"
    elif alpha < 0.95 and alpha > 0: flag = " ↓ DECAY"
    elif abs(alpha - 1.0) < 0.02: flag = " = STABLE"

    if r <= 8 or r % 8 == 0 or r >= 56 or abs(alpha-1) > 0.1:
        print(f"  {r:>3} | {alpha:>18.4f} | {hw_r:>6.1f} → {hw_r1:>7.1f}{flag}")

# Split into phases
growth = [a for r,a,_,_ in alphas if a > 1.05]
stable = [a for r,a,_,_ in alphas if 0.95 <= a <= 1.05]
decay = [a for r,a,_,_ in alphas if a < 0.95]

print(f"\n  Growth phases (α > 1.05): {len(growth)} rounds, avg α = {sum(growth)/max(len(growth),1):.3f}")
print(f"  Stable phases (0.95-1.05): {len(stable)} rounds, avg α = {sum(stable)/max(len(stable),1):.3f}")
print(f"  Decay phases (α < 0.95): {len(decay)} rounds")


# ================================================================
# TEST 3: α as function of CURRENT HW (not round number)
# α(HW) = E[HW_next | HW_current = HW] / HW
# This reveals: is amplification state-dependent or round-dependent?
# ================================================================

print()
print("=" * 70)
print("TEST 3: α as function of current HW(δstate)")
print("=" * 70)

# Collect (HW_current, HW_next) pairs from ALL rounds
hw_pairs = []
for trial in range(N):
    for r in range(64):
        hw_cur = hw_by_round[r][trial]
        hw_nxt = hw_by_round[r+1][trial]
        if hw_cur > 0:
            hw_pairs.append((hw_cur, hw_nxt))

# Bin by current HW
bins = {}
for hw_cur, hw_nxt in hw_pairs:
    bucket = (hw_cur // 8) * 8
    if bucket not in bins:
        bins[bucket] = []
    bins[bucket].append(hw_nxt)

print(f"\n  {'HW_cur':>7} | {'E[HW_next]':>10} | {'α':>6} | {'n':>6} | note")

for bucket in sorted(bins.keys()):
    if len(bins[bucket]) < 10: continue
    avg_next = sum(bins[bucket]) / len(bins[bucket])
    mid = bucket + 4
    alpha = avg_next / mid if mid > 0 else 0

    note = ""
    if alpha > 1.5: note = "STRONG GROWTH"
    elif alpha > 1.05: note = "growth"
    elif alpha < 0.95: note = "DECAY"
    elif abs(alpha - 1.0) < 0.03: note = "equilibrium"

    print(f"  {bucket:>3}-{bucket+7:<3} | {avg_next:>10.1f} | {alpha:>6.3f} | {len(bins[bucket]):>6} | {note}")


# ================================================================
# TEST 4: Fit a MODEL to the amplification curve
# Candidates:
# 1. Logistic: HW(r) = L / (1 + exp(-k*(r-r0)))
# 2. Exponential + cap: HW(r) = min(L, HW0 * α^r)
# 3. Thermostat: HW(r+1) = β*HW(r) + γ
# ================================================================

print()
print("=" * 70)
print("TEST 4: Model fitting")
print("=" * 70)

avg_hw = [sum(hw_by_round[r])/N for r in range(65)]

# Model 3: Thermostat HW[r+1] = β*HW[r] + γ
# Linear regression on (HW[r], HW[r+1]) for r where HW > 10
x_data = []
y_data = []
for r in range(64):
    if avg_hw[r] > 10:
        x_data.append(avg_hw[r])
        y_data.append(avg_hw[r+1])

if len(x_data) > 2:
    n = len(x_data)
    mx = sum(x_data)/n; my = sum(y_data)/n
    sxx = sum((x-mx)**2 for x in x_data)
    sxy = sum((x-mx)*(y-my) for x,y in zip(x_data, y_data))

    beta = sxy / sxx if sxx > 0 else 0
    gamma = my - beta * mx

    # Residual
    residuals = [(y - beta*x - gamma)**2 for x,y in zip(x_data, y_data)]
    rmse = math.sqrt(sum(residuals)/n)

    # Equilibrium: HW* = γ/(1-β)
    if abs(1 - beta) > 0.001:
        hw_eq = gamma / (1 - beta)
    else:
        hw_eq = float('inf')

    print(f"  Thermostat model: HW[r+1] = {beta:.4f} × HW[r] + {gamma:.2f}")
    print(f"  β = {beta:.4f}, γ = {gamma:.2f}")
    print(f"  RMSE = {rmse:.2f}")
    print(f"  Equilibrium HW* = γ/(1-β) = {hw_eq:.1f}")
    print(f"  Actual equilibrium: {avg_hw[20]:.1f} (round 20+)")

    # Convergence rate
    print(f"\n  Starting from HW=1 (1-bit diff):")
    hw_model = 1.0
    for r in range(20):
        hw_model = beta * hw_model + gamma
        actual = avg_hw[r+1] if r+1 < 65 else 128
        err = abs(hw_model - actual)
        if r < 10 or r == 19:
            print(f"    r={r+1:>2}: model={hw_model:>7.1f}, actual={actual:>7.1f}, error={err:>5.1f}")


# ================================================================
# TEST 5: The EXACT equation
# From thermostat: HW[r+1] = β*HW[r] + γ
# Solution: HW[r] = HW* + (HW[0] - HW*) × β^r
# where HW* = γ/(1-β)
#
# This PREDICTS HW at any round from HW at round 0.
# If exact → we have a closed-form equation through all 64 rounds!
# ================================================================

print()
print("=" * 70)
print("TEST 5: Closed-form prediction")
print("=" * 70)

hw_star = gamma / (1 - beta) if abs(1-beta) > 0.001 else 128
hw0 = avg_hw[1]  # HW after round 0 (first nonzero)

print(f"  HW* = {hw_star:.1f}")
print(f"  β = {beta:.4f}")
print(f"  HW[0] = {hw0:.1f}")
print(f"\n  Prediction: HW[r] = {hw_star:.1f} + ({hw0:.1f} - {hw_star:.1f}) × {beta:.4f}^r")
print(f"            = {hw_star:.1f} + {hw0 - hw_star:.1f} × {beta:.4f}^r")
print(f"\n  {'r':>3} | {'predicted':>9} | {'actual':>7} | {'error':>6} | {'error%':>7}")

max_err = 0
for r in range(1, 65):
    predicted = hw_star + (hw0 - hw_star) * (beta ** (r-1))
    actual = avg_hw[r]
    err = predicted - actual
    pct = abs(err) / max(actual, 1) * 100
    max_err = max(max_err, abs(err))

    if r <= 10 or r % 8 == 0 or r >= 56:
        print(f"  {r:>3} | {predicted:>9.1f} | {actual:>7.1f} | {err:>+6.1f} | {pct:>6.1f}%")

print(f"\n  Max error: {max_err:.1f} bits")
print(f"  Model {'ACCURATE' if max_err < 5 else 'APPROXIMATE' if max_err < 15 else 'POOR'}")


# ================================================================
# TEST 6: Does the thermostat equation help for collision?
#
# Collision: HW[64] = 0. From model: 0 = HW* + (HW[0]-HW*) × β^63
# → HW[0] = HW* × (1 - 1/β^63)
# → requires HW[0] = HW* (1 - β^{-63})
#
# But HW* ≈ 128 and β < 1 → β^{-63} huge → HW[0] huge negative.
# Nonsensical: HW is always ≥ 0.
#
# The thermostat PREVENTS HW from reaching 0.
# Starting from any HW[0] > 0: HW[r] → HW* = 128.
# Never reaches 0 (collision).
# ================================================================

print()
print("=" * 70)
print("TEST 6: What the thermostat says about collision")
print("=" * 70)

print(f"  Model: HW[r] = {hw_star:.1f} + (HW[0] - {hw_star:.1f}) × {beta:.4f}^r")
print(f"  For collision at r=64: need HW[64] = 0")
print(f"  0 = {hw_star:.1f} + (HW[0] - {hw_star:.1f}) × {beta:.4f}^63")
print(f"  HW[0] = {hw_star:.1f} × (1 - {beta:.4f}^{{-63}})")
print(f"  β^{{-63}} = {(1/beta)**63:.2e}")
print(f"  HW[0] = {hw_star:.1f} × (1 - {(1/beta)**63:.2e}) = {hw_star * (1 - (1/beta)**63):.2e}")
print(f"\n  Required HW[0] is NEGATIVE ({hw_star * (1 - (1/beta)**63):.0f}).")
print(f"  The thermostat says: starting from ANY positive HW,")
print(f"  the system converges to HW*={hw_star:.0f}, NEVER to 0.")
print(f"\n  To reach HW=0 (collision): need to ESCAPE the thermostat.")
print(f"  Probability of escape by random fluctuation:")
print(f"  P(HW[r+1] < HW[r] - k) depends on std of HW.")

# What's the std of HW at equilibrium?
eq_stds = [math.sqrt(sum((x - avg_hw[r])**2 for x in hw_by_round[r])/N)
           for r in range(20, 65)]
avg_std = sum(eq_stds) / len(eq_stds)

print(f"  std(HW) at equilibrium = {avg_std:.1f}")
print(f"  For HW=0 from HW*={hw_star:.0f}: need {hw_star/avg_std:.1f}σ deviation")
print(f"  P({hw_star/avg_std:.1f}σ) ≈ 2^{{-{hw_star**2/(2*avg_std**2)/math.log(2):.0f}}}")
print(f"  Compare with birthday: 2^128")

bits_needed = hw_star**2 / (2 * avg_std**2) / math.log(2)
print(f"\n  Thermostat escape cost: 2^{bits_needed:.0f}")
print(f"  Birthday cost: 2^128")
print(f"  Ratio: 2^{bits_needed - 128:.0f}")

if bits_needed > 128:
    print(f"\n  Thermostat escape HARDER than birthday.")
    print(f"  The thermostat is THE reason birthday is optimal.")
else:
    print(f"\n  ★ Thermostat escape EASIER than birthday!")
