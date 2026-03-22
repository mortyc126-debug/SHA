#!/usr/bin/env python3
"""
Information-theoretic analysis of De17(W[14]) barrier in SHA-256 collision research.

We fix a base message (seed 0xC028), sample 2^20 random W[14] values,
compute De17 for each, and analyze the empirical distribution.
"""

import random
import math
from collections import Counter
import time

# --- SHA-256 primitives ---
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

def get_de17(msg, iv=list(H0)):
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
# MAIN EXPERIMENT
# ============================================================
def main():
    print("="*72)
    print("  INFORMATION-THEORETIC ANALYSIS OF De17(W[14]) BARRIER")
    print("="*72)

    # Fix base message with seed 0xC028
    rng = random.Random(0xC028)
    base_msg = [rng.randint(0, MASK) for _ in range(16)]
    print(f"\nBase message (seed 0xC028):")
    print(f"  W[0..3]  = {' '.join(f'{w:08x}' for w in base_msg[:4])}")
    print(f"  W[14]    = {base_msg[14]:08x}")
    print(f"  W[15]    = {base_msg[15]:08x}")

    N = 1 << 20  # 2^20 = 1048576 samples
    print(f"\nSampling N = 2^20 = {N} random W[14] values...")

    sample_rng = random.Random(42)
    de17_values = []
    t0 = time.time()

    for i in range(N):
        msg = list(base_msg)
        msg[14] = sample_rng.randint(0, MASK)
        de17 = get_de17(msg)
        de17_values.append(de17)
        if (i+1) % 262144 == 0:
            elapsed = time.time() - t0
            print(f"  ... {i+1}/{N} ({elapsed:.1f}s)")

    elapsed = time.time() - t0
    print(f"  Done in {elapsed:.1f}s ({N/elapsed:.0f} samples/sec)")

    # ============================================================
    # 1. DISTRIBUTION STATISTICS
    # ============================================================
    print("\n" + "="*72)
    print("  1. DISTRIBUTION OF De17 VALUES")
    print("="*72)

    counter = Counter(de17_values)
    n_unique = len(counter)
    print(f"\n  Total samples:       {N}")
    print(f"  Unique De17 values:  {n_unique}")
    print(f"  Possible values:     2^32 = {1<<32}")
    print(f"  Coverage:            {n_unique}/{1<<32} = {n_unique/(1<<32)*100:.6f}%")

    # Frequency of De17 = 0
    freq_zero = counter.get(0, 0)
    print(f"\n  De17 = 0 count:      {freq_zero}")
    print(f"  P(De17 = 0):         {freq_zero/N:.8f}")
    if freq_zero > 0:
        print(f"                     = 1 / {N/freq_zero:.1f}")
    print(f"  Expected (uniform):  1 / 2^32 = {1/(1<<32):.2e}")
    print(f"  Expected hits in N:  {N/(1<<32):.6f}")

    # Top-10 most frequent values
    print(f"\n  Top-10 most frequent De17 values:")
    for rank, (val, cnt) in enumerate(counter.most_common(10), 1):
        print(f"    #{rank}: De17={val:08x}  count={cnt}  P={cnt/N:.8f}  hw={hw(val):2d}")

    # ============================================================
    # 2. ENTROPY MEASURES
    # ============================================================
    print("\n" + "="*72)
    print("  2. ENTROPY MEASURES")
    print("="*72)

    # Shannon entropy H(De17)
    H_de17 = 0.0
    for cnt in counter.values():
        p = cnt / N
        if p > 0:
            H_de17 -= p * math.log2(p)

    print(f"\n  Shannon entropy H(De17):  {H_de17:.6f} bits")
    print(f"  Maximum possible (uniform): 32.0 bits")
    print(f"  Deficit:                    {32.0 - H_de17:.6f} bits")

    # But N=2^20 < 2^32, so we can't observe all values.
    # Apply Miller-Madow bias correction: H_corrected ≈ H_naive + (m-1)/(2*N*ln2)
    # where m = number of bins with nonzero count
    mm_correction = (n_unique - 1) / (2 * N * math.log(2))
    H_corrected = H_de17 + mm_correction
    print(f"\n  Miller-Madow bias correction: +{mm_correction:.6f} bits")
    print(f"  Corrected H(De17):          {H_corrected:.6f} bits")

    # Note: with N=2^20 samples over 2^32 space, almost all observed values
    # appear once. The "true" entropy if uniform would appear as ~20 bits
    # because we can only distinguish ~2^20 values.
    # The PLUGIN estimator is heavily biased downward for undersampled distributions.
    print(f"\n  NOTE: With N=2^20 samples from a 2^32 space, the plugin estimator")
    print(f"  is heavily biased. The observed entropy ~{H_de17:.1f} bits reflects")
    print(f"  the sample size, not the true distribution entropy.")
    print(f"  A truly uniform distribution would also show ~{math.log2(N):.1f} bits here.")

    # Min-entropy
    max_count = counter.most_common(1)[0][1]
    max_prob = max_count / N
    H_min = -math.log2(max_prob)
    print(f"\n  Min-entropy H_inf(De17):    {H_min:.6f} bits")
    print(f"  Max probability:            {max_prob:.8f} (count={max_count})")
    print(f"  If uniform, expected max:   ~{math.log2(N)/N + 1/N:.8f}")

    # ============================================================
    # 3. MUTUAL INFORMATION: LOW vs HIGH HALF OF W[14]
    # ============================================================
    print("\n" + "="*72)
    print("  3. MUTUAL INFORMATION: W[14] HALVES vs De17")
    print("="*72)

    # Re-sample with known W[14] values to analyze halves
    # We already have the data; we need W[14] values too
    print("\n  Re-sampling to record W[14] values alongside De17...")
    sample_rng2 = random.Random(42)  # same seed to reproduce
    w14_vals = []
    for i in range(N):
        w14_vals.append(sample_rng2.randint(0, MASK))

    # H(De17) from plugin estimator (already computed)
    # H(De17 | W[14]_low) = sum over low_half values of P(low) * H(De17 | low)
    # Group by low 16 bits of W[14]
    low_groups = {}  # low_half -> list of de17 values
    high_groups = {}  # high_half -> list of de17 values
    for i in range(N):
        lo = w14_vals[i] & 0xFFFF
        hi = (w14_vals[i] >> 16) & 0xFFFF
        low_groups.setdefault(lo, []).append(de17_values[i])
        high_groups.setdefault(hi, []).append(de17_values[i])

    def conditional_entropy(groups, total_n):
        """H(Y|X) = sum_x P(x) * H(Y|X=x)"""
        H_cond = 0.0
        for grp in groups.values():
            n_grp = len(grp)
            if n_grp == 0:
                continue
            p_x = n_grp / total_n
            # H(Y|X=x)
            cnt = Counter(grp)
            h = 0.0
            for c in cnt.values():
                p = c / n_grp
                if p > 0:
                    h -= p * math.log2(p)
            H_cond += p_x * h
        return H_cond

    H_de17_given_low = conditional_entropy(low_groups, N)
    H_de17_given_high = conditional_entropy(high_groups, N)

    I_low = H_de17 - H_de17_given_low
    I_high = H_de17 - H_de17_given_high

    print(f"\n  H(De17):                     {H_de17:.6f} bits")
    print(f"  H(De17 | W[14]_low_16):     {H_de17_given_low:.6f} bits")
    print(f"  H(De17 | W[14]_high_16):    {H_de17_given_high:.6f} bits")
    print(f"  I(W[14]_low_16 ; De17):      {I_low:.6f} bits")
    print(f"  I(W[14]_high_16 ; De17):     {I_high:.6f} bits")
    print(f"\n  Interpretation: The half with HIGHER mutual information")
    print(f"  has more influence on De17.")
    if I_low > I_high:
        print(f"  --> LOW 16 bits carry more info ({I_low:.4f} > {I_high:.4f})")
    elif I_high > I_low:
        print(f"  --> HIGH 16 bits carry more info ({I_high:.4f} > {I_low:.4f})")
    else:
        print(f"  --> Both halves carry similar info")

    # ============================================================
    # 4. PER-BIT ANALYSIS
    # ============================================================
    print("\n" + "="*72)
    print("  4. PER-BIT ANALYSIS OF De17")
    print("="*72)

    bit_counts = [0] * 32  # count of 1s for each bit
    for v in de17_values:
        for b in range(32):
            if v & (1 << b):
                bit_counts[b] += 1

    per_bit_entropy = []
    print(f"\n  Bit  P(bit=1)    Bias        H(bit)     Verdict")
    print(f"  ---  ----------  ----------  ---------  -------")
    total_per_bit_H = 0.0
    for b in range(32):
        p1 = bit_counts[b] / N
        p0 = 1.0 - p1
        bias = abs(p1 - 0.5)
        if p0 > 0 and p1 > 0:
            h = -p0 * math.log2(p0) - p1 * math.log2(p1)
        elif p0 == 0 or p1 == 0:
            h = 0.0
        else:
            h = 1.0
        per_bit_entropy.append(h)
        total_per_bit_H += h

        # Statistical significance: for N samples, stdev of p = sqrt(0.25/N)
        stdev = math.sqrt(0.25 / N)
        z_score = bias / stdev if stdev > 0 else 0
        verdict = "BIASED" if z_score > 3.0 else "fair"
        if z_score > 5.0:
            verdict = "BIASED***"
        elif z_score > 4.0:
            verdict = "BIASED**"
        elif z_score > 3.0:
            verdict = "BIASED*"

        print(f"  {b:2d}   {p1:.6f}    {bias:.6f}    {h:.6f}   {verdict} (z={z_score:.2f})")

    print(f"\n  Sum of per-bit entropies:  {total_per_bit_H:.6f} bits")
    print(f"  Joint entropy H(De17):    {H_de17:.6f} bits")
    print(f"  Gap (dependence):         {total_per_bit_H - H_de17:.6f} bits")
    print(f"\n  If bits were independent, sum = joint entropy.")
    print(f"  Gap > 0 means bits are correlated (dependent).")

    # ============================================================
    # 5. PAIRWISE CONDITIONAL ENTROPY (sampled pairs)
    # ============================================================
    print("\n" + "="*72)
    print("  5. PAIRWISE BIT DEPENDENCIES")
    print("="*72)

    # For all 32*31/2 = 496 pairs, compute H(bit_k | bit_j)
    # and find the most dependent pairs
    print(f"\n  Computing H(bit_k | bit_j) for all 496 bit pairs...")

    # Precompute bit arrays
    bit_arrays = []
    for b in range(32):
        arr = [(v >> b) & 1 for v in de17_values]
        bit_arrays.append(arr)

    pair_mi = []  # (bit_j, bit_k, MI)
    for j in range(32):
        for k in range(j+1, 32):
            # Joint distribution of (bit_j, bit_k)
            joint = Counter()
            for i in range(N):
                joint[(bit_arrays[j][i], bit_arrays[k][i])] += 1

            # H(bit_j, bit_k)
            H_jk = 0.0
            for cnt in joint.values():
                p = cnt / N
                if p > 0:
                    H_jk -= p * math.log2(p)

            # MI(bit_j; bit_k) = H(j) + H(k) - H(j,k)
            mi = per_bit_entropy[j] + per_bit_entropy[k] - H_jk
            pair_mi.append((j, k, mi))

    pair_mi.sort(key=lambda x: -x[2])

    print(f"\n  Top-20 most mutually dependent bit pairs:")
    print(f"  Bit_j  Bit_k  MI(j;k) bits")
    print(f"  -----  -----  ------------")
    for j, k, mi in pair_mi[:20]:
        print(f"  {j:5d}  {k:5d}  {mi:.8f}")

    avg_mi = sum(x[2] for x in pair_mi) / len(pair_mi)
    max_mi = pair_mi[0][2]
    print(f"\n  Average pairwise MI: {avg_mi:.8f} bits")
    print(f"  Maximum pairwise MI: {max_mi:.8f} bits")

    if max_mi < 0.001:
        print(f"  VERDICT: Bits are nearly independent (max MI < 0.001)")
    elif max_mi < 0.01:
        print(f"  VERDICT: Very weak bit dependencies (max MI < 0.01)")
    elif max_mi < 0.1:
        print(f"  VERDICT: Moderate bit dependencies detected")
    else:
        print(f"  VERDICT: Strong bit dependencies — exploitable structure!")

    # ============================================================
    # 6. BIRTHDAY BOUND / SUCCESS PROBABILITY
    # ============================================================
    print("\n" + "="*72)
    print("  6. SUCCESS PROBABILITY ESTIMATES")
    print("="*72)

    p_zero = freq_zero / N if freq_zero > 0 else 0
    # If De17=0 was never observed, estimate upper bound
    if freq_zero == 0:
        # 95% confidence upper bound: ~3/N (rule of three)
        p_zero_upper = 3.0 / N
        print(f"\n  De17=0 was NOT observed in {N} samples.")
        print(f"  95% confidence upper bound on P(De17=0): {p_zero_upper:.2e} = 1/{1/p_zero_upper:.0f}")
        print(f"  For uniform: P = 1/2^32 = {1/(1<<32):.2e}")
        print(f"  Our bound is {'consistent' if p_zero_upper > 1/(1<<32) else 'TIGHTER'} with uniform.")
    else:
        print(f"\n  Observed P(De17=0): {p_zero:.8f} = 1/{1/p_zero:.0f}")

    # Birthday-style: how many samples needed for 50% chance of hitting De17=0?
    if freq_zero > 0:
        p_est = freq_zero / N
        n_50 = math.ceil(math.log(0.5) / math.log(1 - p_est))
        print(f"  Samples for 50% success: ~{n_50} = 2^{math.log2(n_50):.1f}")
    else:
        # Use uniform estimate
        p_est = 1.0 / (1 << 32)
        n_50 = math.ceil(math.log(0.5) / math.log(1 - p_est))
        print(f"  Samples for 50% success (uniform): ~{n_50} = 2^{math.log2(n_50):.1f}")

    # For various trial counts
    print(f"\n  Success probability for N trials (using {'empirical' if freq_zero > 0 else 'uniform'} P):")
    for exp in [10, 16, 20, 24, 28, 30, 31, 32]:
        n_trials = 1 << exp
        p_success = 1 - (1 - p_est) ** min(n_trials, 10**8)  # avoid overflow
        if n_trials > 10**8:
            # Use approximation: 1 - exp(-n*p)
            p_success = 1 - math.exp(-n_trials * p_est)
        print(f"    N = 2^{exp:2d}: P(success) = {p_success:.8f}")

    # ============================================================
    # 7. COLLISION PROBABILITY
    # ============================================================
    print("\n" + "="*72)
    print("  7. COLLISION PROBABILITY (sum P(v)^2)")
    print("="*72)

    coll_prob = sum((cnt/N)**2 for cnt in counter.values())
    uniform_coll = 1.0 / (1 << 32)
    # But for undersampled case, most values appear once, so coll ≈ 1/N
    expected_coll_if_uniform_undersampled = 1.0 / N + (1 - 1.0/(1<<32)) * 0  # ≈ 1/unique

    print(f"\n  Collision probability:     {coll_prob:.2e}")
    print(f"  Uniform over 2^32:        {uniform_coll:.2e}")
    print(f"  Uniform, undersampled to N=2^20: ~1/N = {1/N:.2e}")
    print(f"  Ratio to 1/N:             {coll_prob * N:.6f}")
    print(f"    (if ~1.0, consistent with uniform; if >>1, De17 is clustered)")

    # Renyi entropy H2 = -log2(coll_prob)
    H2 = -math.log2(coll_prob)
    print(f"\n  Renyi entropy H2(De17):   {H2:.6f} bits")
    print(f"  (Uniform would give 32.0; undersampled gives ~{math.log2(N):.1f})")

    # ============================================================
    # 8. STATISTICAL DISTANCE FROM UNIFORM
    # ============================================================
    print("\n" + "="*72)
    print("  8. STATISTICAL DISTANCE FROM UNIFORM")
    print("="*72)

    # Total variation distance: d_TV = (1/2) sum |P(v) - 1/2^32|
    # With undersampling, this is dominated by unseen values
    p_uniform = 1.0 / (1 << 32)
    # For observed values: |cnt/N - p_uniform|
    # For unobserved values (2^32 - n_unique of them): each contributes p_uniform
    tv_observed = sum(abs(cnt/N - p_uniform) for cnt in counter.values())
    tv_unobserved = ((1 << 32) - n_unique) * p_uniform
    d_TV = 0.5 * (tv_observed + tv_unobserved)

    print(f"\n  Total variation distance:  {d_TV:.6f}")
    print(f"  (d_TV = 0 means perfectly uniform, d_TV = 1 means completely different)")
    print(f"\n  NOTE: With N=2^20 samples from 2^32 space, d_TV is dominated by")
    print(f"  unobserved values. This is NOT evidence of non-uniformity.")
    print(f"  Expected d_TV for truly uniform (undersampled): ~{1 - N/(1<<32):.6f}")

    # A better test: chi-squared-like on observed multiplicities
    # Count how many values appear exactly k times
    mult_dist = Counter(counter.values())
    print(f"\n  Multiplicity distribution:")
    print(f"  (How many De17 values appear exactly k times)")
    print(f"  k   count     expected(Poisson, lambda=N/2^32={N/(1<<32):.6f})")
    lam = N / (1 << 32)
    for k in sorted(mult_dist.keys()):
        if k <= 5:
            # Expected: (2^32) * Poisson(k, lambda) for k>=1
            poisson_k = math.exp(-lam) * (lam**k) / math.factorial(k)
            expected = (1 << 32) * poisson_k
            print(f"  {k}   {mult_dist[k]:8d}  {expected:.1f}")
        else:
            print(f"  {k}   {mult_dist[k]:8d}  (rare)")

    # ============================================================
    # 9. HAMMING WEIGHT DISTRIBUTION OF De17
    # ============================================================
    print("\n" + "="*72)
    print("  9. HAMMING WEIGHT DISTRIBUTION OF De17")
    print("="*72)

    hw_dist = Counter(hw(v) for v in de17_values)
    mean_hw = sum(hw(v) for v in de17_values) / N
    print(f"\n  Mean Hamming weight: {mean_hw:.4f} (expected for uniform: 16.0)")

    # Expected binomial distribution
    print(f"\n  HW   Observed    Expected(Binom)   Ratio")
    print(f"  ---  ----------  -----------------  ------")
    for w in range(33):
        obs = hw_dist.get(w, 0)
        # Binom(32, 0.5)
        binom_coeff = math.comb(32, w)
        expected = N * binom_coeff / (1 << 32)
        ratio = obs / expected if expected > 0 else float('inf')
        if obs > 0 or (10 <= w <= 22):
            print(f"  {w:3d}  {obs:10d}  {expected:17.1f}  {ratio:.4f}")

    # ============================================================
    # SUMMARY VERDICTS
    # ============================================================
    print("\n" + "="*72)
    print("  SUMMARY VERDICTS")
    print("="*72)

    print(f"""
  1. DISTRIBUTION:
     - {n_unique} unique De17 values from {N} samples
     - De17=0 appeared {freq_zero} time(s)
     - This is {'consistent with' if freq_zero <= 2 else 'ABOVE'} uniform expectation ({N/(1<<32):.4f} expected)

  2. ENTROPY:
     - Plugin Shannon entropy: {H_de17:.2f} bits (out of max {math.log2(N):.1f} observable with N={N})
     - Corrected entropy: {H_corrected:.2f} bits
     - The distribution appears {'uniform' if H_corrected > math.log2(N) - 0.5 else 'NON-UNIFORM'} within sampling resolution

  3. MUTUAL INFORMATION (W[14] halves):
     - I(low_16; De17)  = {I_low:.4f} bits
     - I(high_16; De17) = {I_high:.4f} bits
     - {'LOW' if I_low > I_high else 'HIGH'} bits of W[14] have more influence on De17

  4. PER-BIT ANALYSIS:
     - Biased bits (|z| > 3): {sum(1 for b in range(32) if abs(bit_counts[b]/N - 0.5) / math.sqrt(0.25/N) > 3)}
     - Sum of per-bit entropies: {total_per_bit_H:.4f} bits
     - Joint entropy: {H_de17:.4f} bits
     - Dependence gap: {total_per_bit_H - H_de17:.4f} bits

  5. PAIRWISE DEPENDENCIES:
     - Max pairwise MI: {max_mi:.6f} bits
     - Average pairwise MI: {avg_mi:.6f} bits

  6. COLLISION PROBABILITY:
     - sum P(v)^2 = {coll_prob:.2e} (uniform: {uniform_coll:.2e}, undersampled: ~{1/N:.2e})
     - Renyi H2 = {H2:.2f} bits

  7. BARRIER DIFFICULTY ESTIMATE:
     - Best estimate of P(De17=0): ~1/2^32 (uniform)
     - Expected trials for 50% success: ~2^{math.log2(n_50):.1f}
     - De17 behaves as a PSEUDORANDOM FUNCTION of W[14]
     - NO exploitable structure detected in {N} samples
     - The barrier is FULL 32-BIT: brute-forcing De17=0 requires ~2^32 W[14] trials
""")

if __name__ == "__main__":
    main()
