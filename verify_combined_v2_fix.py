"""
Fixed verification: use chi-squared test (correct for discrete HW data)
instead of broken KS test on discrete distribution.
"""

import numpy as np
from scipy import stats
import time

# === SHA-256 reduced-round ===
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
    return bin(x).count('1')

def chi2_hw_test(hw_vals, n_total):
    """Chi-squared test: observed HW histogram vs Binomial(32, 0.5).
    Bins HW values into groups with expected count >= 5."""
    expected_dist = stats.binom(32, 0.5)

    # Count observed frequencies for each HW value
    obs_counts = np.bincount(hw_vals, minlength=33)

    # Bin tails so expected >= 5
    # Find range where expected >= 5/n_total fraction
    lo, hi = 0, 32
    while expected_dist.pmf(lo) * n_total < 5:
        lo += 1
    while expected_dist.pmf(hi) * n_total < 5:
        hi -= 1

    # Build binned observed and expected
    obs_binned = []
    exp_binned = []

    # Lower tail: [0, lo)
    obs_binned.append(np.sum(obs_counts[:lo]))
    exp_binned.append(expected_dist.cdf(lo - 1) * n_total)

    # Individual bins [lo, hi]
    for k in range(lo, hi + 1):
        obs_binned.append(obs_counts[k])
        exp_binned.append(expected_dist.pmf(k) * n_total)

    # Upper tail: (hi, 32]
    obs_binned.append(np.sum(obs_counts[hi + 1:]))
    exp_binned.append((1 - expected_dist.cdf(hi)) * n_total)

    obs_binned = np.array(obs_binned, dtype=float)
    exp_binned = np.array(exp_binned, dtype=float)

    # Remove bins with 0 expected
    mask = exp_binned > 0
    obs_binned = obs_binned[mask]
    exp_binned = exp_binned[mask]

    chi2_stat = np.sum((obs_binned - exp_binned) ** 2 / exp_binned)
    df = len(obs_binned) - 1
    pval = 1 - stats.chi2.cdf(chi2_stat, df)
    return pval, chi2_stat, df

# === Test 1: Null hypothesis with CORRECT test ===
def test_null_chi2(n_trials=20, n_pairs=10000):
    """Random oracle with chi-squared test. Should give p ~ Uniform."""
    print("=" * 60)
    print("TEST 1: Null hypothesis (random oracle) — CHI-SQUARED")
    print(f"  {n_trials} trials, {n_pairs} pairs")
    print("=" * 60)

    pvals = []
    for trial in range(n_trials):
        hw_vals = np.array([hw(int(np.random.randint(0, 0xFFFFFFFF + 1, dtype=np.uint64)))
                            for _ in range(n_pairs)])
        # HW of single random 32-bit value ~ Binom(32, 0.5)
        p, chi2, df = chi2_hw_test(hw_vals, n_pairs)
        pvals.append(p)

    pvals = np.array(pvals)
    print(f"  p-values: min={pvals.min():.4f}, median={np.median(pvals):.4f}, max={pvals.max():.4f}")
    print(f"  p < 0.05: {np.sum(pvals < 0.05)}/{n_trials} (expected ~{n_trials*0.05:.0f})")
    # Uniform check
    _, meta_p = stats.kstest(pvals, 'uniform')
    print(f"  Meta-KS (p-values uniform?): p = {meta_p:.4f}")
    print(f"  => {'PASS: null looks correct' if meta_p > 0.01 else 'FAIL: null still broken'}")

# === Test 2: SHA-256 differential with chi-squared ===
def test_sha256_chi2(rounds_list, n_pairs=10000, n_seeds=5):
    """SHA-256 differential test using chi-squared on HW distribution."""
    print("\n" + "=" * 60)
    print("TEST 2: SHA-256 differential — CHI-SQUARED")
    print(f"  Rounds: {rounds_list}, {n_pairs} pairs, {n_seeds} seeds")
    print("=" * 60)

    delta = 0x80000000

    for nr in rounds_list:
        print(f"\n  --- Round {nr} ---")
        pvals = []
        for seed in range(n_seeds):
            np.random.seed(seed * 1000 + 42)
            hw_vals = np.empty(n_pairs, dtype=np.int32)
            for i in range(n_pairs):
                W = [int(x) for x in np.random.randint(0, 0xFFFFFFFF + 1, 16, dtype=np.uint64)]
                W2 = list(W)
                W2[0] = W[0] ^ delta
                out1 = sha256_rounds(W, nr)
                out2 = sha256_rounds(W2, nr)
                diff_e = out1[4] ^ out2[4]
                hw_vals[i] = hw(diff_e)

            p, chi2, df = chi2_hw_test(hw_vals, n_pairs)
            mean_hw = np.mean(hw_vals)
            std_hw = np.std(hw_vals)
            pvals.append(p)
            sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
            print(f"    seed={seed}: p={p:.6f} (chi2={chi2:.2f}, df={df}), mean_HW={mean_hw:.3f}, std={std_hw:.3f} {sig}")

        pvals = np.array(pvals)
        n_sig = np.sum(pvals < 0.05)
        print(f"    Significant: {n_sig}/{n_seeds} seeds (need majority for real signal)")

# === Test 3: Other tests from combined_v2, also verified ===
def test_higher_order(rounds_list, n_subspaces=3000, dim=4, n_seeds=3):
    """Higher-order differential: XOR-sum over affine subspaces."""
    print("\n" + "=" * 60)
    print("TEST 3: Higher-order differential (binomial test)")
    print(f"  Rounds: {rounds_list}, {n_subspaces} subspaces, dim={dim}")
    print("=" * 60)

    for nr in rounds_list:
        print(f"\n  --- Round {nr} ---")
        for seed in range(n_seeds):
            np.random.seed(seed * 1000 + 7)
            bit0_ones = 0
            for _ in range(n_subspaces):
                base = [int(x) for x in np.random.randint(0, 0xFFFFFFFF + 1, 16, dtype=np.uint64)]
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
                    xor_sum ^= out[4]
                bit0_ones += (xor_sum & 1)

            p = stats.binomtest(bit0_ones, n_subspaces, 0.5).pvalue
            sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
            print(f"    seed={seed}: ones={bit0_ones}/{n_subspaces} ({bit0_ones/n_subspaces:.3f}), p={p:.6f} {sig}")

# === Test 4: Linear bias ===
def test_linear_bias(rounds_list, n_samples=20000, n_seeds=3):
    """Linear bias: W[0] bit 31 vs output e bit 0."""
    print("\n" + "=" * 60)
    print("TEST 4: Linear bias (binomial test)")
    print(f"  Rounds: {rounds_list}, {n_samples} samples")
    print("=" * 60)

    for nr in rounds_list:
        print(f"\n  --- Round {nr} ---")
        for seed in range(n_seeds):
            np.random.seed(seed * 1000 + 13)
            agree = 0
            for _ in range(n_samples):
                W = [int(x) for x in np.random.randint(0, 0xFFFFFFFF + 1, 16, dtype=np.uint64)]
                out = sha256_rounds(W, nr)
                if ((W[0] >> 31) & 1) == (out[4] & 1):
                    agree += 1
            p = stats.binomtest(agree, n_samples, 0.5).pvalue
            bias = agree / n_samples - 0.5
            sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
            print(f"    seed={seed}: agree={agree}/{n_samples} (bias={bias:+.5f}), p={p:.6f} {sig}")

if __name__ == "__main__":
    t0 = time.time()

    test_null_chi2(n_trials=20, n_pairs=10000)

    test_sha256_chi2(
        rounds_list=[4, 6, 8, 10, 12],
        n_pairs=10000,
        n_seeds=3
    )

    test_higher_order(
        rounds_list=[4, 6, 8, 10, 12],
        n_subspaces=2000,
        n_seeds=3
    )

    test_linear_bias(
        rounds_list=[4, 6, 8, 10, 12],
        n_samples=20000,
        n_seeds=3
    )

    elapsed = time.time() - t0
    print(f"\n{'='*60}")
    print(f"Total time: {elapsed:.1f}s")
