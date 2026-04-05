#!/usr/bin/env python3
"""
The Doctor: tracking δe through arithmetic lens.

Bit-level doctor: sees HW(δe XOR) — loses signal at round ~5.
Arithmetic doctor: sees δe mod k — loses signal at round ???

If mod k sees further → arithmetic axis survives rotation mixing.
This would mean: Ch/Maj is the LAST barrier, not rotation.
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

N = 5000
random.seed(42)

# ================================================================
# Collect δe (additive) and δe (XOR) at every round
# ================================================================

de_add = [[] for _ in range(65)]  # additive: e2 - e1 mod 2^32
de_xor = [[] for _ in range(65)]  # XOR: e2 ^ e1
da_add = [[] for _ in range(65)]
da_xor = [[] for _ in range(65)]

for _ in range(N):
    M1 = [random.randint(0,MASK32) for _ in range(16)]
    M2 = list(M1); M2[0] ^= 1

    ss1 = sha_states(M1)
    ss2 = sha_states(M2)

    for r in range(65):
        de_add[r].append((ss2[r][4] - ss1[r][4]) & MASK32)
        de_xor[r].append(ss2[r][4] ^ ss1[r][4])
        da_add[r].append((ss2[r][0] - ss1[r][0]) & MASK32)
        da_xor[r].append(ss2[r][0] ^ ss1[r][0])

# ================================================================
# The Doctor's instruments: mod k uniformity + bit-level HW
# ================================================================

print("=" * 70)
print("THE DOCTOR: δe tracking — bit-level vs arithmetic")
print("=" * 70)

def chi2_mod(values, m):
    """Chi-squared statistic for values mod m (normalized by DOF)."""
    counts = [0]*m
    for v in values:
        counts[v % m] += 1
    expected = len(values) / m
    chi2 = sum((c - expected)**2 / expected for c in counts)
    return chi2 / (m - 1)

def hw_deviation(xor_values):
    """How far is E[HW(XOR)] from 16 (random)?"""
    avg_hw = sum(hw(v) for v in xor_values) / len(xor_values)
    return abs(avg_hw - 16)

# Track δe through ALL metrics simultaneously
print(f"\n  N={N}, δM = W[0] bit 0")
print(f"\n  {'r':>3} | {'HW_dev':>6} | {'mod2':>5} | {'mod4':>5} | {'mod8':>5} | {'mod16':>5} | {'mod256':>6} | best_metric")

last_round = {
    'hw_dev': 0, 'mod2': 0, 'mod4': 0, 'mod8': 0, 'mod16': 0, 'mod256': 0
}

THRESH = 3.0  # chi2/DOF threshold for significance

for r in range(65):
    hw_dev = hw_deviation(de_xor[r])
    m2 = chi2_mod(de_add[r], 2)
    m4 = chi2_mod(de_add[r], 4)
    m8 = chi2_mod(de_add[r], 8)
    m16 = chi2_mod(de_add[r], 16)
    m256 = chi2_mod(de_add[r], 256)

    # Update last structured round
    if hw_dev > 1.0: last_round['hw_dev'] = r
    if m2 > THRESH: last_round['mod2'] = r
    if m4 > THRESH: last_round['mod4'] = r
    if m8 > THRESH: last_round['mod8'] = r
    if m16 > THRESH: last_round['mod16'] = r
    if m256 > THRESH: last_round['mod256'] = r

    # Best metric = the one still showing structure
    metrics = {'HW': hw_dev > 1.0, 'mod2': m2 > THRESH, 'mod4': m4 > THRESH,
               'mod8': m8 > THRESH, 'mod16': m16 > THRESH, 'mod256': m256 > THRESH}
    active = [k for k, v in metrics.items() if v]
    best = ','.join(active) if active else '-'

    if r <= 10 or r % 8 == 0 or r >= 56 or active:
        print(f"  {r:>3} | {hw_dev:>6.2f} | {m2:>5.1f} | {m4:>5.1f} | {m8:>5.1f} | {m16:>5.1f} | {m256:>6.1f} | {best}")

# ================================================================
# SAME for δa
# ================================================================

print()
print("=" * 70)
print("THE DOCTOR: δa tracking — same metrics")
print("=" * 70)

last_round_a = {
    'hw_dev': 0, 'mod2': 0, 'mod4': 0, 'mod8': 0, 'mod16': 0, 'mod256': 0
}

print(f"\n  {'r':>3} | {'HW_dev':>6} | {'mod2':>5} | {'mod4':>5} | {'mod8':>5} | {'mod16':>5} | {'mod256':>6} | best_metric")

for r in range(65):
    hw_dev = hw_deviation(da_xor[r])
    m2 = chi2_mod(da_add[r], 2)
    m4 = chi2_mod(da_add[r], 4)
    m8 = chi2_mod(da_add[r], 8)
    m16 = chi2_mod(da_add[r], 16)
    m256 = chi2_mod(da_add[r], 256)

    if hw_dev > 1.0: last_round_a['hw_dev'] = r
    if m2 > THRESH: last_round_a['mod2'] = r
    if m4 > THRESH: last_round_a['mod4'] = r
    if m8 > THRESH: last_round_a['mod8'] = r
    if m16 > THRESH: last_round_a['mod16'] = r
    if m256 > THRESH: last_round_a['mod256'] = r

    metrics = {'HW': hw_dev > 1.0, 'mod2': m2 > THRESH, 'mod4': m4 > THRESH,
               'mod8': m8 > THRESH, 'mod16': m16 > THRESH, 'mod256': m256 > THRESH}
    active = [k for k, v in metrics.items() if v]
    best = ','.join(active) if active else '-'

    if r <= 10 or r % 8 == 0 or r >= 56 or active:
        print(f"  {r:>3} | {hw_dev:>6.2f} | {m2:>5.1f} | {m4:>5.1f} | {m8:>5.1f} | {m16:>5.1f} | {m256:>6.1f} | {best}")

# ================================================================
# SUMMARY: which metric sees furthest?
# ================================================================

print()
print("=" * 70)
print("SUMMARY: Last structured round per metric")
print("=" * 70)

print(f"\n  δe (e-register):")
for metric, r in sorted(last_round.items(), key=lambda x: -x[1]):
    bar = "█" * r
    print(f"    {metric:>7}: round {r:>2} {bar}")

print(f"\n  δa (a-register):")
for metric, r in sorted(last_round_a.items(), key=lambda x: -x[1]):
    bar = "█" * r
    print(f"    {metric:>7}: round {r:>2} {bar}")

# ================================================================
# CROSS-CHECK: is mod 8 structure at late rounds REAL or noise?
# Run 5 independent seeds for the latest structured round
# ================================================================

print()
print("=" * 70)
print("CROSS-CHECK: verify latest structure with 5 seeds")
print("=" * 70)

# Find the metric+register with latest structure
best_metric_r = max(max(last_round.values()), max(last_round_a.values()))
if best_metric_r > 4:
    # Which metric and register?
    for reg_name, lr in [("δe", last_round), ("δa", last_round_a)]:
        for metric, r in lr.items():
            if r == best_metric_r:
                print(f"\n  Checking: {reg_name} {metric} at round {r}")

                for seed in range(5):
                    random.seed(seed * 9999 + 7)
                    vals = []
                    for _ in range(N):
                        M1 = [random.randint(0,MASK32) for _ in range(16)]
                        M2 = list(M1); M2[0] ^= 1
                        ss1 = sha_states(M1); ss2 = sha_states(M2)
                        if reg_name == "δe":
                            d = (ss2[r][4] - ss1[r][4]) & MASK32
                        else:
                            d = (ss2[r][0] - ss1[r][0]) & MASK32

                        if 'mod' in metric:
                            m = int(metric.replace('mod',''))
                            vals.append(d % m)
                        else:
                            vals.append(hw(ss1[r][0 if 'a' in reg_name else 4] ^ ss2[r][0 if 'a' in reg_name else 4]))

                    if 'mod' in metric:
                        m = int(metric.replace('mod',''))
                        counts = [0]*m
                        for v in vals: counts[v%m]+=1
                        chi2 = sum((c-N/m)**2/(N/m) for c in counts)/(m-1)
                        sig = "★ CONFIRMED" if chi2 > THRESH else "noise"
                        print(f"    seed {seed}: χ²/DOF = {chi2:.2f} {sig}")
                    else:
                        avg = sum(vals)/N
                        sig = "★" if abs(avg-16) > 1 else "noise"
                        print(f"    seed {seed}: E[HW] = {avg:.2f} {sig}")
else:
    print(f"\n  No structure beyond round 4. Nothing to cross-check.")
