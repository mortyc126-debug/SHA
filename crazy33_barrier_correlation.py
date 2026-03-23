#!/usr/bin/env python3
"""
crazy33_barrier_correlation.py — Barrier Correlation Analysis

Measure statistical dependencies between De17, De18, De19, De20 barriers.
If barriers are positively correlated, satisfying one helps satisfy others.
"""

import random
import math
import numpy as np
from collections import defaultdict

MASK = 0xFFFFFFFF
K = [0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
     0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
     0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da]
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

def expand_W(W16, num_rounds=24):
    W=list(W16)
    for i in range(16, num_rounds):
        W.append(add32(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    return W

def get_all_de(msg, max_round=21, iv=None):
    if iv is None:
        iv = list(H0)
    W=list(msg); Wp=list(W); Wp[0]^=0x80000000
    # Force De1..De16 = 0 by adjusting Wp[t]
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
    # Now expand and run full chain
    We=expand_W(msg, max_round)
    Wpe=expand_W(Wp, max_round)
    s1=list(iv); s2=list(iv)
    des = {}
    for t in range(max_round):
        s1=sha_round(s1,We[t],K[t])
        s2=sha_round(s2,Wpe[t],K[t])
        des[t+1] = {'e': s1[4]^s2[4], 'a': s1[0]^s2[0]}
    return des

def bits_to_vec(x, nbits=32):
    """Convert 32-bit int to numpy array of bits."""
    return np.array([(x >> i) & 1 for i in range(nbits)], dtype=np.float64)

def mutual_info_bits(x_arr, y_arr):
    """Mutual information between two arrays of 32-bit values, treating each as binary."""
    n = len(x_arr)
    # Bin by (hw_x, hw_y) is approximate; instead compute per-bit MI and sum
    # For efficiency, compute MI between each bit pair and average
    # Actually compute MI(X;Y) where X,Y are 32-bit. Too expensive directly.
    # Instead: estimate from HW correlation as proxy, and compute per-bit MI sum.
    total_mi = 0.0
    for bi in range(32):
        bx = np.array([(v >> bi) & 1 for v in x_arr], dtype=np.float64)
        for bj in range(32):
            by = np.array([(v >> bj) & 1 for v in y_arr], dtype=np.float64)
            # 2x2 contingency
            p11 = np.mean(bx * by)
            p10 = np.mean(bx * (1-by))
            p01 = np.mean((1-bx) * by)
            p00 = np.mean((1-bx) * (1-by))
            px1 = p11 + p10
            px0 = p01 + p00
            py1 = p11 + p01
            py0 = p10 + p00
            mi = 0.0
            for pxy, px, py in [(p11,px1,py1),(p10,px1,py0),(p01,px0,py1),(p00,px0,py0)]:
                if pxy > 1e-12 and px > 1e-12 and py > 1e-12:
                    mi += pxy * math.log2(pxy / (px * py))
            total_mi += mi
    return total_mi

def main():
    print("="*78)
    print("CRAZY-33: Barrier Correlation Analysis")
    print("="*78)
    print()

    # ── Phase 1: Collect samples ──
    random.seed(0xC033)
    base_msg = [random.getrandbits(32) for _ in range(16)]

    N_SMALL = 10000
    N_LARGE = 2**20  # ~1M for joint probability estimation

    rounds_of_interest = [17, 18, 19, 20]
    labels_e = [f"De{r}" for r in rounds_of_interest]
    labels_a = [f"Da{r}" for r in rounds_of_interest]

    print(f"Base message seed: 0xC033")
    print(f"Varying: W[14] randomly")
    print(f"Small sample: {N_SMALL}, Large sample: {N_LARGE}")
    print()

    # ── Collect small sample for detailed analysis ──
    print("Phase 1: Collecting small sample (N={})...".format(N_SMALL))
    de_samples = {r: [] for r in rounds_of_interest}
    da_samples = {r: [] for r in rounds_of_interest}
    hw_e = {r: [] for r in rounds_of_interest}
    hw_a = {r: [] for r in rounds_of_interest}

    for i in range(N_SMALL):
        msg = list(base_msg)
        msg[14] = random.getrandbits(32)
        des = get_all_de(msg, max_round=21)
        for r in rounds_of_interest:
            de_samples[r].append(des[r]['e'])
            da_samples[r].append(des[r]['a'])
            hw_e[r].append(hw(des[r]['e']))
            hw_a[r].append(hw(des[r]['a']))

    # ── Section 2: Hamming Weight Statistics ──
    print("\n" + "="*78)
    print("SECTION 2: Hamming Weight Statistics")
    print("="*78)
    for r in rounds_of_interest:
        arr = np.array(hw_e[r])
        print(f"  HW(De{r}): mean={arr.mean():.2f}, std={arr.std():.2f}, "
              f"min={arr.min()}, max={arr.max()}")
    print()
    for r in rounds_of_interest:
        arr = np.array(hw_a[r])
        print(f"  HW(Da{r}): mean={arr.mean():.2f}, std={arr.std():.2f}, "
              f"min={arr.min()}, max={arr.max()}")

    # ── Section 3: HW Correlation Matrix (Pearson) ──
    print("\n" + "="*78)
    print("SECTION 3: Pearson Correlation of Hamming Weights")
    print("="*78)

    all_labels = labels_e + labels_a
    all_hw = [np.array(hw_e[r]) for r in rounds_of_interest] + \
             [np.array(hw_a[r]) for r in rounds_of_interest]

    n = len(all_labels)
    corr_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            corr_matrix[i][j] = np.corrcoef(all_hw[i], all_hw[j])[0, 1]

    # Print matrix
    header = "         " + "  ".join(f"{l:>7s}" for l in all_labels)
    print(header)
    for i in range(n):
        row = f"{all_labels[i]:>7s}  " + "  ".join(f"{corr_matrix[i][j]:+7.4f}" for j in range(n))
        print(row)

    # ── Section 4: Conditional HW Analysis ──
    print("\n" + "="*78)
    print("SECTION 4: Conditional Hamming Weight Analysis")
    print("="*78)
    print("  When HW(DeR) <= threshold, what is avg HW(DeS)?")
    print()

    for threshold in [5, 8, 10, 12]:
        print(f"  --- Threshold: HW <= {threshold} ---")
        for r1 in rounds_of_interest:
            mask = np.array(hw_e[r1]) <= threshold
            count = mask.sum()
            if count == 0:
                print(f"    HW(De{r1})<={threshold}: 0 samples")
                continue
            line = f"    HW(De{r1})<={threshold} ({count} samples) => "
            parts = []
            for r2 in rounds_of_interest:
                if r2 == r1:
                    continue
                cond_mean = np.array(hw_e[r2])[mask].mean()
                uncond_mean = np.array(hw_e[r2]).mean()
                parts.append(f"avg HW(De{r2})={cond_mean:.2f} (uncond={uncond_mean:.2f})")
            print(line + ", ".join(parts))
        print()

    # ── Section 5: Joint HW Distribution for De17 vs De18 ──
    print("="*78)
    print("SECTION 5: Joint HW Distribution — De17 vs De18")
    print("="*78)

    hw17 = np.array(hw_e[17])
    hw18 = np.array(hw_e[18])
    print(f"  Pearson(HW(De17), HW(De18)) = {np.corrcoef(hw17, hw18)[0,1]:+.6f}")
    print()

    # Joint probability of both being low
    for thr in [5, 8, 10]:
        both_low = np.sum((hw17 <= thr) & (hw18 <= thr))
        either17 = np.sum(hw17 <= thr)
        either18 = np.sum(hw18 <= thr)
        p17 = either17 / N_SMALL
        p18 = either18 / N_SMALL
        p_joint = both_low / N_SMALL
        p_indep = p17 * p18
        ratio = p_joint / p_indep if p_indep > 0 else float('inf')
        print(f"  HW<={thr}: P(De17 low)={p17:.6f}, P(De18 low)={p18:.6f}")
        print(f"           P(both low)={p_joint:.6f}, P(indep)={p_indep:.6f}, ratio={ratio:.4f}")
        if ratio > 1:
            print(f"           => POSITIVE correlation (ratio {ratio:.2f}x) — barriers help each other!")
        elif ratio < 1:
            print(f"           => NEGATIVE correlation (ratio {ratio:.2f}x) — barriers fight each other")
        else:
            print(f"           => Independent")
        print()

    # ── Section 6: Per-bit Correlation (top pairs) ──
    print("="*78)
    print("SECTION 6: Per-bit Correlation — De17 vs De18 (top 20 pairs)")
    print("="*78)

    de17_arr = np.array(de_samples[17], dtype=np.uint32)
    de18_arr = np.array(de_samples[18], dtype=np.uint32)

    bit_corrs = []
    for bi in range(32):
        bx = ((de17_arr >> bi) & 1).astype(np.float64)
        for bj in range(32):
            by = ((de18_arr >> bj) & 1).astype(np.float64)
            c = np.corrcoef(bx, by)[0, 1]
            if not np.isnan(c):
                bit_corrs.append((abs(c), c, bi, bj))

    bit_corrs.sort(reverse=True)
    print(f"  {'De17 bit':>10s} {'De18 bit':>10s} {'Correlation':>12s}")
    for rank, (ac, c, bi, bj) in enumerate(bit_corrs[:20]):
        print(f"  {bi:>10d} {bj:>10d} {c:>+12.6f}")

    # Summary stats on bit correlations
    all_c = [c for _, c, _, _ in bit_corrs]
    print(f"\n  All {len(all_c)} bit-pair correlations:")
    print(f"    mean |corr| = {np.mean(np.abs(all_c)):.6f}")
    print(f"    max  |corr| = {np.max(np.abs(all_c)):.6f}")
    print(f"    Pairs with |corr| > 0.05: {sum(1 for c in all_c if abs(c)>0.05)}")
    print(f"    Pairs with |corr| > 0.10: {sum(1 for c in all_c if abs(c)>0.10)}")

    # ── Section 7: Mutual Information (subset of bits for speed) ──
    print("\n" + "="*78)
    print("SECTION 7: Mutual Information (full 32x32 bit-pair)")
    print("="*78)
    print("  Computing MI between barrier pairs...")

    for r1, r2 in [(17,18),(17,19),(17,20),(18,19),(18,20),(19,20)]:
        mi = mutual_info_bits(de_samples[r1][:2000], de_samples[r2][:2000])
        print(f"  MI(De{r1}; De{r2}) = {mi:.4f} bits")

    print("\n  Cross e/a barriers:")
    for r1, r2 in [(17,17),(17,18),(18,18),(18,19)]:
        mi = mutual_info_bits(de_samples[r1][:2000], da_samples[r2][:2000])
        print(f"  MI(De{r1}; Da{r2}) = {mi:.4f} bits")

    # ── Section 8: Large-sample joint zero probability ──
    print("\n" + "="*78)
    print("SECTION 8: Large-Sample Joint Zero Probability (N={})".format(N_LARGE))
    print("="*78)
    print("  Collecting large sample...")

    random.seed(0xC033_8888)
    zero_counts = {r: 0 for r in rounds_of_interest}
    joint_zero = defaultdict(int)  # (r1,r2) -> count both zero
    all4_zero = 0
    low_hw_counts = {r: {5: 0, 8: 0, 10: 0} for r in rounds_of_interest}
    joint_low = defaultdict(lambda: defaultdict(int))  # (r1,r2) -> {thr: count}

    for i in range(N_LARGE):
        msg = list(base_msg)
        msg[14] = random.getrandbits(32)
        des = get_all_de(msg, max_round=21)

        de_vals = {}
        hw_vals = {}
        for r in rounds_of_interest:
            v = des[r]['e']
            de_vals[r] = v
            hw_vals[r] = hw(v)
            if v == 0:
                zero_counts[r] += 1
            for thr in [5, 8, 10]:
                if hw(v) <= thr:
                    low_hw_counts[r][thr] += 1

        # Joint zeros
        for i1, r1 in enumerate(rounds_of_interest):
            for r2 in rounds_of_interest[i1+1:]:
                if de_vals[r1] == 0 and de_vals[r2] == 0:
                    joint_zero[(r1,r2)] += 1
                for thr in [5, 8, 10]:
                    if hw_vals[r1] <= thr and hw_vals[r2] <= thr:
                        joint_low[(r1,r2)][thr] += 1

        if all(de_vals[r] == 0 for r in rounds_of_interest):
            all4_zero += 1

        if (i+1) % 200000 == 0:
            print(f"    ... {i+1}/{N_LARGE} samples processed")

    print()
    print("  Individual barrier zero probabilities:")
    for r in rounds_of_interest:
        p = zero_counts[r] / N_LARGE
        log2p = math.log2(p) if p > 0 else float('-inf')
        print(f"    P(De{r}=0) = {zero_counts[r]}/{N_LARGE} = {p:.2e} (2^{log2p:.1f})")

    print()
    print("  Joint zero probabilities:")
    for i1, r1 in enumerate(rounds_of_interest):
        for r2 in rounds_of_interest[i1+1:]:
            jz = joint_zero[(r1,r2)]
            p_joint = jz / N_LARGE
            p1 = zero_counts[r1] / N_LARGE
            p2 = zero_counts[r2] / N_LARGE
            p_indep = p1 * p2
            ratio = p_joint / p_indep if p_indep > 0 else float('inf')
            print(f"    P(De{r1}=0 AND De{r2}=0) = {jz}/{N_LARGE}", end="")
            if jz > 0:
                print(f" = {p_joint:.2e}, indep={p_indep:.2e}, ratio={ratio:.2f}x")
            else:
                print(f" = 0 (need more samples; indep={p_indep:.2e})")

    print()
    print(f"  P(all four De17..De20 = 0) = {all4_zero}/{N_LARGE}")

    # Joint low-HW probabilities
    print()
    print("  Joint low-HW probabilities (large sample):")
    for thr in [5, 8, 10]:
        print(f"\n  --- HW <= {thr} ---")
        for i1, r1 in enumerate(rounds_of_interest):
            for r2 in rounds_of_interest[i1+1:]:
                jl = joint_low[(r1,r2)][thr]
                p1 = low_hw_counts[r1][thr] / N_LARGE
                p2 = low_hw_counts[r2][thr] / N_LARGE
                p_joint = jl / N_LARGE
                p_indep = p1 * p2
                ratio = p_joint / p_indep if p_indep > 0 else float('inf')
                print(f"    P(HW(De{r1})<={thr} AND HW(De{r2})<={thr}) = {p_joint:.6f}, "
                      f"indep={p_indep:.6f}, ratio={ratio:.4f}")

    # ── Section 9: Conditional boost ──
    print("\n" + "="*78)
    print("SECTION 9: Conditional Zero Probability (from large sample)")
    print("="*78)
    for r1 in rounds_of_interest:
        if zero_counts[r1] == 0:
            print(f"  De{r1}=0: never observed, skip conditional")
            continue
        for r2 in rounds_of_interest:
            if r2 == r1:
                continue
            jz = joint_zero.get((min(r1,r2), max(r1,r2)), 0)
            p_cond = jz / zero_counts[r1] if zero_counts[r1] > 0 else 0
            p_uncond = zero_counts[r2] / N_LARGE
            boost = p_cond / p_uncond if p_uncond > 0 else float('inf')
            print(f"  P(De{r2}=0 | De{r1}=0) = {jz}/{zero_counts[r1]}", end="")
            if zero_counts[r1] > 0 and p_uncond > 0:
                print(f" = {p_cond:.4e}  (unconditional={p_uncond:.4e}, boost={boost:.2f}x)")
            else:
                print()

    # ── Summary ──
    print("\n" + "="*78)
    print("SUMMARY")
    print("="*78)

    # Key correlations
    c_17_18 = np.corrcoef(np.array(hw_e[17]), np.array(hw_e[18]))[0, 1]
    c_17_19 = np.corrcoef(np.array(hw_e[17]), np.array(hw_e[19]))[0, 1]
    c_18_19 = np.corrcoef(np.array(hw_e[18]), np.array(hw_e[19]))[0, 1]
    c_17_20 = np.corrcoef(np.array(hw_e[17]), np.array(hw_e[20]))[0, 1]

    print(f"  Pearson(HW) correlations:")
    print(f"    De17-De18: {c_17_18:+.4f}")
    print(f"    De17-De19: {c_17_19:+.4f}")
    print(f"    De17-De20: {c_17_20:+.4f}")
    print(f"    De18-De19: {c_18_19:+.4f}")
    print()

    if c_17_18 > 0.01:
        print("  FINDING: De17 and De18 are POSITIVELY correlated.")
        print("           Low De17 predicts lower De18 — barriers HELP each other.")
        print("           Nonlinear coupling is FAVORABLE for collision search!")
    elif c_17_18 < -0.01:
        print("  FINDING: De17 and De18 are NEGATIVELY correlated.")
        print("           Low De17 predicts HIGHER De18 — barriers FIGHT each other.")
        print("           Nonlinear coupling makes collision search HARDER.")
    else:
        print("  FINDING: De17 and De18 are essentially UNCORRELATED.")
        print("           Barriers are statistically independent (nonlinear coupling is negligible).")
        print("           Joint probability = product of marginals.")

    # Effective joint difficulty
    p17 = zero_counts[17] / N_LARGE if zero_counts[17] > 0 else None
    p18 = zero_counts[18] / N_LARGE if zero_counts[18] > 0 else None
    if p17 and p18:
        log2_p17 = math.log2(p17)
        log2_p18 = math.log2(p18)
        print(f"\n  Marginal zero rates: P(De17=0) ~ 2^{log2_p17:.1f}, P(De18=0) ~ 2^{log2_p18:.1f}")
        print(f"  If independent: P(both=0) ~ 2^{log2_p17+log2_p18:.1f}")
        jz = joint_zero.get((17,18), 0)
        if jz > 0:
            log2_joint = math.log2(jz / N_LARGE)
            print(f"  Measured joint:  P(both=0) ~ 2^{log2_joint:.1f}")
            diff = log2_joint - (log2_p17 + log2_p18)
            print(f"  Correlation effect: {diff:+.1f} bits")
        else:
            print(f"  Joint zero not observed in {N_LARGE} samples — consistent with independence")

    print("\n" + "="*78)
    print("Done.")

if __name__ == "__main__":
    main()
