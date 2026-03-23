#!/usr/bin/env python3
"""
CRAZY-30: Bit-position algebraic gradient of De17.

Hypothesis: De17 bit k has algebraic degree that INCREASES with k,
because carry chains accumulate from LSB to MSB.

Plan:
1. Fix base message (seed 0xC030), vary 16-bit subset of W[14].
2. For each output bit k (0..31) of De17, compute full truth table on 2^16 inputs.
3. Apply Möbius transform to get ANF for each output bit.
4. Record: degree, number of terms, density at each degree level.
5. Cross-validate with bits 16-31 of W[14].
6. Compute entropy and bias for each output bit.
"""

import math
import time
from collections import defaultdict

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
def sha_round(st,w,k):
    a,b,c,d,e,f,g,h=st
    T1=add32(h,Sig1(e),Ch(e,f,g),k,w); T2=add32(Sig0(a),Maj(a,b,c))
    return [add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
def expand_W(W16):
    W=list(W16)
    for i in range(16,20): W.append(add32(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    return W

def get_de17(msg, iv=list(H0)):
    W=list(msg); Wp=list(W); Wp[0]^=0x80000000
    s=list(iv); sp=list(iv)
    s=sha_round(s,W[0],K[0]); sp=sha_round(sp,Wp[0],K[0])
    for t in range(1,16):
        a,b,c,d,e,f,g,h=s; a2,b2,c2,d2,e2,f2,g2,h2=sp
        tp=add32(h,Sig1(e),Ch(e,f,g),K[t]); tp2=add32(h2,Sig1(e2),Ch(e2,f2,g2),K[t])
        target=add32(d,tp,W[t]); Wp[t]=(target-d2-tp2)&MASK
        s=sha_round(s,W[t],K[t]); sp=sha_round(sp,Wp[t],K[t])
    We=expand_W(msg); Wpe=expand_W(Wp)
    s1=list(iv); s2=list(iv)
    for t in range(17): s1=sha_round(s1,We[t],K[t]); s2=sha_round(s2,Wpe[t],K[t])
    return s1[4]^s2[4]

def mobius_transform(tt, n):
    """In-place Möbius transform. tt is list of 2^n GF(2) values."""
    for i in range(n):
        for j in range(1 << n):
            if j & (1 << i):
                tt[j] ^= tt[j ^ (1 << i)]
    return tt

def popcount(x):
    return bin(x).count('1')

def entropy(p):
    """Binary entropy of a bit with P(1) = p."""
    if p <= 0 or p >= 1:
        return 0.0
    return -p * math.log2(p) - (1-p) * math.log2(1-p)

def run_experiment(base_msg, input_word_idx, input_bit_lo, n_input_bits, label):
    """
    Vary n_input_bits of W[input_word_idx] starting at bit input_bit_lo.
    For each of 32 output bits of De17, compute ANF via Möbius transform.
    """
    n = n_input_bits
    N = 1 << n
    print(f"\n{'='*80}")
    print(f"EXPERIMENT: {label}")
    print(f"Varying W[{input_word_idx}] bits {input_bit_lo}..{input_bit_lo+n-1} ({n} bits, {N} inputs)")
    print(f"{'='*80}")

    # Step 1: Build truth tables for all 32 output bits
    print(f"Building truth tables ({N} evaluations)...")
    t0 = time.time()

    # truth_tables[k] = list of 2^n GF(2) values for output bit k
    truth_tables = [[0]*N for _ in range(32)]
    ones_count = [0]*32  # for bias computation

    for idx in range(N):
        msg = list(base_msg)
        # Set the input bits
        # Clear the target bits first, then set from idx
        clear_mask = ((1 << n) - 1) << input_bit_lo
        msg[input_word_idx] = (msg[input_word_idx] & ~clear_mask) & MASK
        msg[input_word_idx] = (msg[input_word_idx] | ((idx << input_bit_lo) & clear_mask)) & MASK

        de17 = get_de17(msg)
        for k in range(32):
            bit_val = (de17 >> k) & 1
            truth_tables[k][idx] = bit_val
            ones_count[k] += bit_val

    elapsed = time.time() - t0
    print(f"Truth tables built in {elapsed:.1f}s")

    # Step 2: Möbius transform each truth table
    print("Computing Möbius transforms...")
    t0 = time.time()

    results = []
    for k in range(32):
        tt = list(truth_tables[k])  # copy
        bias = ones_count[k] / N
        ent = entropy(bias)

        mobius_transform(tt, n)

        # Analyze ANF
        max_degree = 0
        num_terms = 0
        degree_counts = defaultdict(int)

        for s in range(N):
            if tt[s]:
                deg = popcount(s)
                num_terms += 1
                degree_counts[deg] += 1
                if deg > max_degree:
                    max_degree = deg

        # Constant term (s=0) has degree 0
        # Total possible terms at each degree: C(n, d)
        density = num_terms / N if N > 0 else 0

        results.append({
            'bit': k,
            'degree': max_degree,
            'num_terms': num_terms,
            'density': density,
            'bias': bias,
            'entropy': ent,
            'degree_counts': dict(degree_counts),
        })

    elapsed = time.time() - t0
    print(f"Möbius transforms done in {elapsed:.1f}s")

    # Step 3: Print results table
    print(f"\n{'─'*100}")
    print(f"{'Bit':>4} {'Degree':>7} {'#Terms':>7} {'Density':>8} {'Bias':>7} {'Entropy':>8} │ "
          f"{'Deg0':>5} {'Deg1':>5} {'Deg2':>5} {'Deg3':>5} {'Deg4':>5} {'Deg5+':>6}")
    print(f"{'─'*100}")

    for r in results:
        dc = r['degree_counts']
        d0 = dc.get(0, 0)
        d1 = dc.get(1, 0)
        d2 = dc.get(2, 0)
        d3 = dc.get(3, 0)
        d4 = dc.get(4, 0)
        d5plus = sum(v for k2,v in dc.items() if k2 >= 5)
        print(f"{r['bit']:>4} {r['degree']:>7} {r['num_terms']:>7} {r['density']:>8.4f} "
              f"{r['bias']:>7.4f} {r['entropy']:>8.4f} │ "
              f"{d0:>5} {d1:>5} {d2:>5} {d3:>5} {d4:>5} {d5plus:>6}")

    # Step 4: Summary statistics
    print(f"\n{'─'*60}")
    print("SUMMARY:")
    degrees = [r['degree'] for r in results]
    print(f"  Min degree: {min(degrees)} (at bits {[r['bit'] for r in results if r['degree']==min(degrees)]})")
    print(f"  Max degree: {max(degrees)} (at bits {[r['bit'] for r in results if r['degree']==max(degrees)]})")
    print(f"  Mean degree: {sum(degrees)/len(degrees):.2f}")

    # Check monotonicity
    mono_increases = sum(1 for i in range(31) if degrees[i+1] >= degrees[i])
    strict_increases = sum(1 for i in range(31) if degrees[i+1] > degrees[i])
    print(f"  Monotonic non-decreasing steps: {mono_increases}/31")
    print(f"  Strict increases: {strict_increases}/31")

    # Correlation between bit position and degree
    mean_k = sum(range(32)) / 32
    mean_d = sum(degrees) / 32
    cov = sum((k - mean_k) * (degrees[k] - mean_d) for k in range(32)) / 32
    var_k = sum((k - mean_k)**2 for k in range(32)) / 32
    var_d = sum((d - mean_d)**2 for d in degrees) / 32
    if var_k > 0 and var_d > 0:
        corr = cov / math.sqrt(var_k * var_d)
    else:
        corr = 0.0
    print(f"  Correlation(bit_pos, degree): {corr:.4f}")

    # Degree distribution summary
    print(f"\n  Degree distribution across all bits:")
    all_degree_counts = defaultdict(int)
    for r in results:
        for d, c in r['degree_counts'].items():
            all_degree_counts[d] += c
    for d in sorted(all_degree_counts.keys()):
        print(f"    Degree {d:>2}: {all_degree_counts[d]:>7} terms total")

    # Bias statistics
    biases = [r['bias'] for r in results]
    print(f"\n  Bias range: [{min(biases):.4f}, {max(biases):.4f}]")
    print(f"  Mean bias: {sum(biases)/len(biases):.4f}")
    print(f"  Mean entropy: {sum(r['entropy'] for r in results)/32:.4f}")

    # Degree profile: fraction of terms at each degree for each bit
    print(f"\n  Fraction of ANF terms by degree (selected bits):")
    print(f"  {'Bit':>4} {'Frac_d1':>8} {'Frac_d2':>8} {'Frac_d3':>8} {'Frac_d4+':>8}")
    for k in [0, 1, 2, 3, 4, 7, 15, 23, 29, 31]:
        r = results[k]
        nt = max(r['num_terms'], 1)
        dc = r['degree_counts']
        f1 = dc.get(1, 0) / nt
        f2 = dc.get(2, 0) / nt
        f3 = dc.get(3, 0) / nt
        f4p = sum(v for k2,v in dc.items() if k2 >= 4) / nt
        print(f"  {k:>4} {f1:>8.4f} {f2:>8.4f} {f3:>8.4f} {f4p:>8.4f}")

    return results


def main():
    print("CRAZY-30: Bit-position algebraic gradient of De17")
    print("=" * 60)

    # Generate base message from seed 0xC030
    import random
    rng = random.Random(0xC030)
    base_msg = [rng.randint(0, MASK) for _ in range(16)]
    print(f"Base message W[0..3]: {[hex(w) for w in base_msg[:4]]}")

    # Experiment A: Vary low 16 bits of W[14]
    results_lo = run_experiment(base_msg, 14, 0, 16, "W[14] bits 0-15 (LOW)")

    # Experiment B: Vary high 16 bits of W[14]
    results_hi = run_experiment(base_msg, 14, 16, 16, "W[14] bits 16-31 (HIGH)")

    # Cross-comparison
    print(f"\n{'='*80}")
    print("CROSS-VALIDATION: Comparing LOW vs HIGH input bits")
    print(f"{'='*80}")
    print(f"{'Bit':>4} {'Deg_LO':>7} {'Deg_HI':>7} {'Match':>6} │ {'Bias_LO':>8} {'Bias_HI':>8}")
    print(f"{'─'*60}")
    matches = 0
    for k in range(32):
        dl = results_lo[k]['degree']
        dh = results_hi[k]['degree']
        m = 'YES' if dl == dh else 'no'
        if dl == dh: matches += 1
        bl = results_lo[k]['bias']
        bh = results_hi[k]['bias']
        print(f"{k:>4} {dl:>7} {dh:>7} {m:>6} │ {bl:>8.4f} {bh:>8.4f}")
    print(f"\nDegree matches: {matches}/32")

    # Final gradient analysis
    print(f"\n{'='*80}")
    print("GRADIENT ANALYSIS")
    print(f"{'='*80}")
    deg_lo = [r['degree'] for r in results_lo]
    deg_hi = [r['degree'] for r in results_hi]

    # Split into lower half (bits 0-15) and upper half (bits 16-31) of OUTPUT
    lo_out_deg_lo = deg_lo[:16]
    hi_out_deg_lo = deg_lo[16:]
    print(f"Using W[14] bits 0-15 as input:")
    print(f"  Output bits 0-15  mean degree: {sum(lo_out_deg_lo)/16:.2f}")
    print(f"  Output bits 16-31 mean degree: {sum(hi_out_deg_lo)/16:.2f}")
    print(f"  Gradient (high - low): {sum(hi_out_deg_lo)/16 - sum(lo_out_deg_lo)/16:.2f}")

    lo_out_deg_hi = deg_hi[:16]
    hi_out_deg_hi = deg_hi[16:]
    print(f"Using W[14] bits 16-31 as input:")
    print(f"  Output bits 0-15  mean degree: {sum(lo_out_deg_hi)/16:.2f}")
    print(f"  Output bits 16-31 mean degree: {sum(hi_out_deg_hi)/16:.2f}")
    print(f"  Gradient (high - low): {sum(hi_out_deg_hi)/16 - sum(lo_out_deg_hi)/16:.2f}")

    # Hypothesis test: is degree correlated with bit position?
    print(f"\nHYPOTHESIS: De17 bit k has degree increasing with k")
    for label, degs in [("LOW input", deg_lo), ("HIGH input", deg_hi)]:
        mean_k = sum(range(32)) / 32
        mean_d = sum(degs) / 32
        cov = sum((k - mean_k) * (degs[k] - mean_d) for k in range(32)) / 32
        var_k = sum((k - mean_k)**2 for k in range(32)) / 32
        var_d = sum((d - mean_d)**2 for d in degs) / 32
        if var_k > 0 and var_d > 0:
            corr = cov / math.sqrt(var_k * var_d)
        else:
            corr = 0.0
        verdict = "CONFIRMED" if corr > 0.3 else ("WEAK" if corr > 0.1 else "REJECTED")
        print(f"  {label}: correlation = {corr:.4f} → {verdict}")

    print(f"\nDone.")


if __name__ == '__main__':
    main()
