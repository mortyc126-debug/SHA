"""
Combined multi-statistic distinguisher for reduced-round SHA-256.
Lean implementation targeting < 90 second runtime.
"""

import numpy as np
from scipy import stats
import time

# === SHA-256 constants and implementation ===
K = [0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
     0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174]
IV = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]

def rotr(x, n):
    return ((x >> n) | (x << (32 - n))) & 0xFFFFFFFF

def sha256_rounds(W, nr):
    w = list(W[:16])
    for i in range(16, max(nr, 16)):
        s0 = rotr(w[i-15], 7) ^ rotr(w[i-15], 18) ^ (w[i-15] >> 3)
        s1 = rotr(w[i-2], 17) ^ rotr(w[i-2], 19) ^ (w[i-2] >> 10)
        w.append((w[i-16] + s0 + w[i-7] + s1) & 0xFFFFFFFF)
    a, b, c, d, e, f, g, h = IV
    for i in range(nr):
        S1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25)
        ch = (e & f) ^ ((~e) & g) & 0xFFFFFFFF
        t1 = (h + S1 + ch + K[i] + w[i]) & 0xFFFFFFFF
        S0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22)
        maj = (a & b) ^ (a & c) ^ (b & c)
        t2 = (S0 + maj) & 0xFFFFFFFF
        h, g, f, e, d, c, b, a = g, f, e, (d + t1) & 0xFFFFFFFF, c, b, a, (t1 + t2) & 0xFFFFFFFF
    return [a, b, c, d, e, f, g, h]

def hw(x):
    """Hamming weight of a 32-bit integer."""
    x = x - ((x >> 1) & 0x55555555)
    x = (x & 0x33333333) + ((x >> 2) & 0x33333333)
    return (((x + (x >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24

def rand_msg():
    """Generate random 16-word message block."""
    return [int(x) for x in np.random.randint(0, 0xFFFFFFFF + 1, 16, dtype=np.uint64)]

# =========================================================
# Phase 1: Four statistical tests
# =========================================================

def test_higher_order_differential(nr, n_subspaces=2000, dim=4):
    """Higher-order differential: XOR-sum over affine subspaces of dim=4 in W[0]."""
    bit0_ones = 0
    for _ in range(n_subspaces):
        base = rand_msg()
        # Pick dim random bit positions in W[0] for the subspace
        bits = np.random.choice(32, dim, replace=False)
        xor_sum = 0
        for mask_val in range(1 << dim):
            w = list(base)
            delta = 0
            for j in range(dim):
                if mask_val & (1 << j):
                    delta ^= (1 << bits[j])
            w[0] = base[0] ^ delta
            out = sha256_rounds(w, nr)
            xor_sum ^= out[4]  # e-register
        bit0_ones += (xor_sum & 1)
    # Under random oracle: bit 0 of XOR-sum is uniform => ~50% ones
    pval = stats.binomtest(bit0_ones, n_subspaces, 0.5).pvalue
    return pval

def test_linear_bias(nr, n_samples=10000):
    """Linear bias: correlation between W[0] bit 31 and output e bit 0."""
    agreements = 0
    for _ in range(n_samples):
        w = rand_msg()
        out = sha256_rounds(w, nr)
        in_bit = (w[0] >> 31) & 1
        out_bit = out[4] & 1  # e bit 0
        if in_bit == out_bit:
            agreements += 1
    pval = stats.binomtest(agreements, n_samples, 0.5).pvalue
    return pval

def test_differential(nr, n_pairs=10000):
    """Differential: ΔW[0]=0x80000000, KS-test on HW of output diff e."""
    hw_vals = np.empty(n_pairs, dtype=np.int32)
    delta = 0x80000000
    for i in range(n_pairs):
        w1 = rand_msg()
        w2 = list(w1)
        w2[0] = w1[0] ^ delta
        out1 = sha256_rounds(w1, nr)
        out2 = sha256_rounds(w2, nr)
        diff_e = out1[4] ^ out2[4]
        hw_vals[i] = hw(diff_e)
    # Under random oracle, HW(diff_e) ~ Binomial(32, 0.5) => mean=16
    # Use KS test against the expected CDF
    expected = stats.binom(32, 0.5)
    ks_stat, pval = stats.kstest(hw_vals, expected.cdf)
    return pval

def test_conditional_bias(nr, n_samples=10000):
    """Conditional bias: among msgs with HW(W[0])>20, chi-squared on e bit 0."""
    bit_counts = [0, 0]  # [zeros, ones]
    for _ in range(n_samples):
        w = rand_msg()
        if hw(w[0]) > 20:
            out = sha256_rounds(w, nr)
            bit_counts[out[4] & 1] += 1
    total = sum(bit_counts)
    if total < 20:
        return 1.0  # not enough samples
    expected_freq = [total / 2.0, total / 2.0]
    chi2_stat = sum((o - e) ** 2 / e for o, e in zip(bit_counts, expected_freq))
    pval = 1.0 - stats.chi2.cdf(chi2_stat, df=1)
    return pval

def fisher_combine(pvals):
    """Fisher's method to combine p-values."""
    pvals_clipped = [max(p, 1e-300) for p in pvals]
    X = -2.0 * sum(np.log(p) for p in pvals_clipped)
    df = 2 * len(pvals_clipped)
    combined_p = 1.0 - stats.chi2.cdf(X, df)
    return combined_p

# =========================================================
# Phase 1 execution
# =========================================================

def run_phase1():
    print("=" * 72)
    print("PHASE 1: Multi-statistic tests for rounds 8-12")
    print("=" * 72)
    # Use reduced sample sizes for speed
    results = {}
    test_names = ["Higher-Order Diff", "Linear Bias", "Differential", "Conditional Bias"]

    for nr in range(8, 13):
        t0 = time.time()
        p1 = test_higher_order_differential(nr, n_subspaces=2000, dim=4)
        p2 = test_linear_bias(nr, n_samples=10000)
        p3 = test_differential(nr, n_pairs=10000)
        p4 = test_conditional_bias(nr, n_samples=10000)
        pvals = [p1, p2, p3, p4]
        combined = fisher_combine(pvals)
        elapsed = time.time() - t0

        results[nr] = {"pvals": pvals, "combined": combined, "time": elapsed}

        print(f"\n--- Round {nr} (elapsed {elapsed:.1f}s) ---")
        for name, p in zip(test_names, pvals):
            sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
            print(f"  {name:20s}: p = {p:.6e} {sig}")
        sig_c = "DISTINGUISHED" if combined < 0.05 else "not distinguished"
        print(f"  {'Combined (Fisher)':20s}: p = {combined:.6e}  => {sig_c}")

    return results

# =========================================================
# Phase 2: Adaptive query budget (rounds 8-10 only)
# =========================================================

def run_adaptive_test(nr, n_samples):
    """Run all 4 tests with given sample size, return combined p-value."""
    # Scale subspace count proportionally
    n_sub = max(200, n_samples // 5)
    p1 = test_higher_order_differential(nr, n_subspaces=n_sub, dim=4)
    p2 = test_linear_bias(nr, n_samples=n_samples)
    p3 = test_differential(nr, n_pairs=n_samples)
    p4 = test_conditional_bias(nr, n_samples=n_samples)
    return fisher_combine([p1, p2, p3, p4]), [p1, p2, p3, p4]

def run_phase2():
    print("\n" + "=" * 72)
    print("PHASE 2: Adaptive query budget (rounds 8-10)")
    print("=" * 72)
    adaptive_results = {}

    for nr in range(8, 11):
        t0 = time.time()
        N = 1000
        best_p = 1.0
        best_N = None
        while N <= 32000:
            elapsed_so_far = time.time() - t0
            if elapsed_so_far > 20:  # safety cap per round
                print(f"  Round {nr}: time cap reached at N={N}")
                break
            combined_p, _ = run_adaptive_test(nr, N)
            print(f"  Round {nr}, N={N:5d}: combined p = {combined_p:.6e}")
            if combined_p < best_p:
                best_p = combined_p
                best_N = N
            if combined_p < 0.01:
                break
            N *= 2

        adaptive_results[nr] = {"best_N": best_N, "best_p": best_p}
        status = f"p < 0.01 at N={best_N}" if best_p < 0.01 else f"best p={best_p:.4e} at N={best_N}"
        print(f"  => Round {nr}: {status}")

    return adaptive_results

# =========================================================
# Phase 3: Summary table
# =========================================================

def run_phase3(phase1_results, phase2_results):
    print("\n" + "=" * 72)
    print("PHASE 3: Summary Table")
    print("=" * 72)
    test_names = ["HO-Diff", "LinBias", "Diff", "CondBias"]

    header = f"{'Round':>5} | {'Best Test':>10} | {'Best p':>12} | {'Combined p':>12} | {'Verdict':>16} | {'Est. Queries':>12}"
    print(header)
    print("-" * len(header))

    for nr in range(8, 13):
        r = phase1_results[nr]
        pvals = r["pvals"]
        best_idx = int(np.argmin(pvals))
        best_test = test_names[best_idx]
        best_p = pvals[best_idx]
        combined = r["combined"]
        verdict = "DISTINGUISHED" if combined < 0.05 else "NOT distinguished"

        # Estimated queries: from phase2 if available, else from phase1
        if nr in phase2_results and phase2_results[nr]["best_p"] < 0.01:
            est_q = phase2_results[nr]["best_N"] * 4  # 4 tests
        else:
            # phase1 used ~10000 per test => ~40000 total, 2000 subspaces * 16 = 32000
            est_q = 10000 * 3 + 2000 * 16  # approximate total hash evaluations

        print(f"{nr:>5} | {best_test:>10} | {best_p:>12.4e} | {combined:>12.4e} | {verdict:>16} | {est_q:>12}")

# =========================================================
# Main
# =========================================================

if __name__ == "__main__":
    np.random.seed(42)
    t_start = time.time()

    print("Combined Multi-Statistic Distinguisher for Reduced-Round SHA-256")
    print(f"Started at {time.strftime('%H:%M:%S')}")
    print()

    phase1_results = run_phase1()

    t_phase1 = time.time() - t_start
    print(f"\n[Phase 1 completed in {t_phase1:.1f}s]")

    # Only run Phase 2 if we have time budget remaining
    if t_phase1 < 60:
        phase2_results = run_phase2()
    else:
        print("\n[Skipping Phase 2 to stay within time budget]")
        phase2_results = {}

    t_phase2 = time.time() - t_start
    print(f"\n[Phase 1+2 completed in {t_phase2:.1f}s]")

    run_phase3(phase1_results, phase2_results)

    total = time.time() - t_start
    print(f"\nTotal runtime: {total:.1f}s")
    if total < 90:
        print("OK: Within 90-second budget.")
    else:
        print("WARNING: Exceeded 90-second budget.")
