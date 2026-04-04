#!/usr/bin/env python3
"""
NK P2 extension: Can we find M with LOW P-count on ALL barrier rounds simultaneously?

Key question: are P-counts at rounds 17,18,19,20 INDEPENDENT or CORRELATED?
If independent: joint min = product of individual mins → still expensive
If correlated: joint min << product → possible shortcut

Also: can we SEARCH for low-P configurations via hill climbing on W[0]?
"""

import random
from collections import defaultdict

MASK32 = 0xFFFFFFFF
K = [0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2]
IV = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&MASK32
def sig0(x): return rotr(x,2)^rotr(x,13)^rotr(x,22)
def sig1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def ssig0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def ssig1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def ch(e,f,g): return (e&f)^(~e&g)&MASK32
def maj(a,b,c): return (a&b)^(a&c)^(b&c)
def hw(x): return bin(x&MASK32).count('1')

def run_sha(M, R):
    W = list(M)+[0]*(64-len(M))
    for i in range(16,64): W[i]=(ssig1(W[i-2])+W[i-7]+ssig0(W[i-15])+W[i-16])&MASK32
    a,b,c,d,e,f,g,h = IV
    states = [{'a':a,'b':b,'c':c,'d':d,'e':e,'f':f,'g':g,'h':h}]
    T1s = []
    for r in range(R):
        T1=(h+sig1(e)+ch(e,f,g)+K[r]+W[r])&MASK32
        T2=(sig0(a)+maj(a,b,c))&MASK32
        T1s.append(T1)
        h,g,f=g,f,e; e=(d+T1)&MASK32; d,c,b=c,b,a; a=(T1+T2)&MASK32
        states.append({'a':a,'b':b,'c':c,'d':d,'e':e,'f':f,'g':g,'h':h})
    return states, T1s, W

def wang_cascade(W_base, dW0, R):
    W1 = list(W_base)
    DW = [0]*16; DW[0] = dW0
    for step in range(min(R-1,15)):
        wi = step+1
        W2 = [(W1[i]+DW[i])&MASK32 for i in range(16)]
        s1,_,_ = run_sha(W1, step+2)
        s2,_,_ = run_sha(W2, step+2)
        de = (s2[-1]['e']-s1[-1]['e'])&MASK32
        DW[wi] = (-de)&MASK32
    W2 = [(W1[i]+DW[i])&MASK32 for i in range(16)]
    return W1, W2, DW

def count_p(a_val, b_val):
    """Count P-positions in GPK string of a+b."""
    p = 0
    for bit in range(32):
        ab = (a_val >> bit) & 1
        bb = (b_val >> bit) & 1
        if ab != bb:
            p += 1
    return p

def get_p_counts_r17_20(W_base, dW0=1):
    """Get P-counts at rounds 17-20 for a Wang pair."""
    W1, W2, DW = wang_cascade(W_base, dW0, 21)
    s1, t1_1, _ = run_sha(W1, 21)
    s2, t1_2, _ = run_sha(W2, 21)

    p_counts = []
    for r in range(17, 21):
        if r < len(t1_1):
            d1 = s1[r]['d']
            T1_v1 = t1_1[r]
            p_counts.append(count_p(d1, T1_v1))
        else:
            p_counts.append(16)
    return p_counts, sum(p_counts)

# ================================================================
# TEST 1: Correlation between P-counts at different barrier rounds
# ================================================================

print("=" * 70)
print("TEST 1: Correlation of P-counts across barrier rounds")
print("=" * 70)

N = 2000
random.seed(42)

all_p = {17: [], 18: [], 19: [], 20: []}
all_sums = []

for trial in range(N):
    W_base = [random.randint(0, MASK32) for _ in range(16)]
    p_counts, p_sum = get_p_counts_r17_20(W_base)
    for i, r in enumerate(range(17, 21)):
        all_p[r].append(p_counts[i])
    all_sums.append(p_sum)

# Pearson correlation between rounds
import math

def pearson(x, y):
    n = len(x)
    mx, my = sum(x)/n, sum(y)/n
    sx = math.sqrt(sum((xi-mx)**2 for xi in x)/n)
    sy = math.sqrt(sum((yi-my)**2 for yi in y)/n)
    if sx == 0 or sy == 0: return 0
    return sum((xi-mx)*(yi-my) for xi, yi in zip(x, y)) / (n * sx * sy)

print(f"\nN={N} Wang pairs")
print(f"\nPearson correlation matrix:")
print(f"       r=17   r=18   r=19   r=20")
for r1 in range(17, 21):
    row = f"r={r1}"
    for r2 in range(17, 21):
        c = pearson(all_p[r1], all_p[r2])
        row += f"  {c:+.3f}"
    print(row)

print(f"\nE[P_sum(17-20)] = {sum(all_sums)/N:.1f}")
print(f"std[P_sum]      = {(sum((x-sum(all_sums)/N)**2 for x in all_sums)/N)**0.5:.1f}")
print(f"min[P_sum]      = {min(all_sums)}")
print(f"max[P_sum]      = {max(all_sums)}")

# If independent: E[sum] = 4*16 = 64, std = sqrt(4*var_single)
var_single = sum((x - 16)**2 for x in all_p[17]) / N
print(f"\nvar(P_single) = {var_single:.1f}")
print(f"std(P_single) = {var_single**0.5:.1f}")
print(f"If independent: std(sum) = {(4*var_single)**0.5:.1f}")

# Distribution of P_sum
print(f"\nP_sum distribution:")
buckets = defaultdict(int)
for s in all_sums:
    buckets[(s//4)*4] += 1
for b in sorted(buckets.keys()):
    pct = buckets[b]/N*100
    bar = '#' * int(pct)
    print(f"  {b:>3}-{b+3:<3}: {buckets[b]:>4} ({pct:>5.1f}%) {bar}")

# ================================================================
# TEST 2: Joint minimum — how often are ALL rounds low?
# ================================================================

print()
print("=" * 70)
print("TEST 2: Joint low-P events")
print("=" * 70)

thresholds = [14, 12, 10, 8]

for thresh in thresholds:
    joint_count = 0
    for trial in range(N):
        if all(all_p[r][trial] <= thresh for r in range(17, 21)):
            joint_count += 1

    # If independent: P(all ≤ thresh) = P(single ≤ thresh)^4
    single_probs = {}
    for r in range(17, 21):
        single_probs[r] = sum(1 for x in all_p[r] if x <= thresh) / N

    p_indep = 1.0
    for r in range(17, 21):
        p_indep *= single_probs[r]

    p_joint = joint_count / N
    ratio = p_joint / max(p_indep, 1e-10)

    print(f"\n  Threshold P ≤ {thresh}:")
    print(f"    Per-round: r17={single_probs[17]:.3f} r18={single_probs[18]:.3f} r19={single_probs[19]:.3f} r20={single_probs[20]:.3f}")
    print(f"    P(joint)  = {p_joint:.4f} ({joint_count}/{N})")
    print(f"    P(indep)  = {p_indep:.4f}")
    print(f"    Ratio     = {ratio:.2f}x {'★ CORRELATED' if ratio > 1.5 else '≈ INDEPENDENT' if ratio > 0.67 else '★ ANTI-CORRELATED'}")


# ================================================================
# TEST 3: Hill climbing on W[0] to minimize P_sum
# ================================================================

print()
print("=" * 70)
print("TEST 3: Hill climbing to minimize joint P-count")
print("=" * 70)

N_HC = 20  # number of HC runs
HC_STEPS = 300

best_overall = 64
best_W = None

for run in range(N_HC):
    W_base = [random.randint(0, MASK32) for _ in range(16)]
    _, current_sum = get_p_counts_r17_20(W_base)
    best_local = current_sum

    for step in range(HC_STEPS):
        # Flip random bit of W[0]
        bit = random.randint(0, 31)
        W_try = list(W_base)
        W_try[0] ^= (1 << bit)

        _, new_sum = get_p_counts_r17_20(W_try)

        if new_sum < current_sum:
            W_base = W_try
            current_sum = new_sum
            if new_sum < best_local:
                best_local = new_sum

    if best_local < best_overall:
        best_overall = best_local
        best_W = list(W_base)
        p_detail, _ = get_p_counts_r17_20(W_base)
        print(f"  Run {run:>2}: P_sum={best_local:>2} (P17={p_detail[0]:>2} P18={p_detail[1]:>2} P19={p_detail[2]:>2} P20={p_detail[3]:>2}) ★ NEW BEST")
    elif run < 5 or best_local <= best_overall + 4:
        p_detail, _ = get_p_counts_r17_20(W_base)
        print(f"  Run {run:>2}: P_sum={best_local:>2} (P17={p_detail[0]:>2} P18={p_detail[1]:>2} P19={p_detail[2]:>2} P20={p_detail[3]:>2})")

print(f"\n  Best P_sum = {best_overall}")
print(f"  Random E[P_sum] = {sum(all_sums)/N:.1f}")
print(f"  HC improvement = {sum(all_sums)/N - best_overall:.1f} bits")

# ================================================================
# TEST 4: Expand search — also vary W[1..15]
# ================================================================

print()
print("=" * 70)
print("TEST 4: HC on ALL W[0..15] to minimize P_sum")
print("=" * 70)

N_HC2 = 10
HC_STEPS2 = 500

best_overall2 = 64

for run in range(N_HC2):
    W_base = [random.randint(0, MASK32) for _ in range(16)]
    _, current_sum = get_p_counts_r17_20(W_base)

    for step in range(HC_STEPS2):
        # Flip random bit of random W[k]
        word = random.randint(0, 15)
        bit = random.randint(0, 31)
        W_try = list(W_base)
        W_try[word] ^= (1 << bit)

        _, new_sum = get_p_counts_r17_20(W_try)

        if new_sum < current_sum:
            W_base = W_try
            current_sum = new_sum

    if current_sum < best_overall2:
        best_overall2 = current_sum
        p_detail, _ = get_p_counts_r17_20(W_base)
        print(f"  Run {run:>2}: P_sum={current_sum:>2} (P17={p_detail[0]:>2} P18={p_detail[1]:>2} P19={p_detail[2]:>2} P20={p_detail[3]:>2}) ★ BEST")
    elif run < 3:
        p_detail, _ = get_p_counts_r17_20(W_base)
        print(f"  Run {run:>2}: P_sum={current_sum:>2} (P17={p_detail[0]:>2} P18={p_detail[1]:>2} P19={p_detail[2]:>2} P20={p_detail[3]:>2})")

print(f"\n  Best P_sum (all W) = {best_overall2}")
print(f"  Best P_sum (W[0] only) = {best_overall}")
print(f"  Random E[P_sum] = {sum(all_sums)/N:.1f}")

# ================================================================
# TEST 5: What's the Λ (total NK cost) for the best path found?
# ================================================================

print()
print("=" * 70)
print("TEST 5: Λ computation for best found path")
print("=" * 70)

# Best path: Wang cascade r=1-16 (free) + barrier r=17-20 (P_sum bits)
# Plus rounds 21-63 (43 rounds at ~16 P each)

lambda_cascade = 0  # rounds 1-16: free (Wang)
lambda_barrier = best_overall2  # rounds 17-20: measured P_sum
lambda_post = 43 * 16  # rounds 21-63: assume random (16 P/round)

lambda_total = lambda_cascade + lambda_barrier + lambda_post

print(f"  Λ_cascade (r=1-16):    {lambda_cascade:>5} bits (Wang, free)")
print(f"  Λ_barrier (r=17-20):   {lambda_barrier:>5} bits (HC-optimized)")
print(f"  Λ_post    (r=21-63):   {lambda_post:>5} bits (random, 43×16)")
print(f"  ─────────────────────────────")
print(f"  Λ_total:               {lambda_total:>5} bits")
print(f"  C(path) = 2^{lambda_total}")
print(f"  Birthday = 2^128")
print(f"  Ratio: path/birthday = 2^{lambda_total - 128}")

if lambda_total > 128:
    print(f"\n  Path WORSE than birthday by {lambda_total - 128} bits")
    print(f"  Bottleneck: Λ_post = {lambda_post} ({lambda_post*100/lambda_total:.0f}% of total)")
    print(f"  Even with perfect barrier (Λ_barrier=0): Λ_total = {lambda_post} >> 128")
    print(f"\n  CONCLUSION: The post-barrier zone (r=21-63) dominates.")
    print(f"  43 rounds × 16 P-bits = 688 >> 128.")
    print(f"  No barrier optimization can overcome the post-barrier cost.")
    print(f"  This confirms: the REAL barrier is Zone 2 (r=21-59), not round 17.")
else:
    print(f"\n  Path BETTER than birthday by {128 - lambda_total} bits!")

# ================================================================
# TEST 6: Reality check — are post-barrier P-counts really 16?
# Or does the path structure help?
# ================================================================

print()
print("=" * 70)
print("TEST 6: Post-barrier P-counts for Wang pairs (are they really 16?)")
print("=" * 70)

N6 = 300
random.seed(999)

post_p = defaultdict(list)

for trial in range(N6):
    W_base = [random.randint(0, MASK32) for _ in range(16)]
    W1, W2, DW = wang_cascade(W_base, 1, 32)
    s1, t1_1, _ = run_sha(W1, 32)

    for r in range(17, 32):
        if r < len(t1_1):
            d1 = s1[r]['d']
            T1_v1 = t1_1[r]
            p = count_p(d1, T1_v1)
            post_p[r].append(p)

print(f"\nN={N6}")
print(f"\n{'r':>3} | {'E[P]':>5} | {'min':>4} {'max':>4} | {'std':>5}")
for r in range(17, 32):
    if post_p[r]:
        ps = post_p[r]
        ep = sum(ps)/len(ps)
        mn = min(ps)
        mx = max(ps)
        std = (sum((x-ep)**2 for x in ps)/len(ps))**0.5
        print(f"{r:>3} | {ep:>5.1f} | {mn:>4} {mx:>4} | {std:>5.1f}")
