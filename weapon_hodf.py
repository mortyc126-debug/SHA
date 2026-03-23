#!/usr/bin/env python3
"""
WEAPON: Higher-Order Differential + Filtering (HODF) Combined Attack
=====================================================================

Combines two strong findings for reduced-round SHA-256:
  1. Higher-order differentials through 8 rounds (k=4, advantage ~0.20)
  2. Message pair filtering by intermediate state (5-11x near-collision boost)

Tool 1: Filtered higher-order distinguisher
Tool 2: Adaptive bit selection for affine subspace
Tool 3: Iterated higher-order statistics (chi-squared)
Tool 4: Combined weapon — minimum evaluations for 99% confidence distinguisher
"""

import random
import time
import math
from itertools import combinations
from collections import defaultdict

# ---------------------------------------------------------------------------
# SHA-256 constants and primitives
# ---------------------------------------------------------------------------
MASK = 0xFFFFFFFF

K = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
    0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
    0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
    0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
    0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
    0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
]

IV = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
]

def rotr(x, n):
    return ((x >> n) | (x << (32 - n))) & MASK

def Sig0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sig1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def sig0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sig1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)
def Ch(e, f, g):  return ((e & f) ^ (~e & g)) & MASK
def Maj(a, b, c): return ((a & b) ^ (a & c) ^ (b & c)) & MASK
def add32(*args):
    s = 0
    for a in args:
        s = (s + a) & MASK
    return s
def hw(x): return bin(x & MASK).count('1')


def expand_schedule(W0_15, nr):
    """Expand 16-word message to nr schedule words."""
    W = list(W0_15[:16])
    for i in range(16, nr):
        W.append(add32(sig1(W[i-2]), W[i-7], sig0(W[i-15]), W[i-16]))
    return W


def sha256_t_rounds(msg, t):
    """t-round SHA-256 compression from standard IV. Returns 8-word state."""
    W = expand_schedule(msg, max(t, 16))
    a, b, c, d, e, f, g, h = IV
    for i in range(t):
        T1 = add32(h, Sig1(e), Ch(e, f, g), K[i], W[i])
        T2 = add32(Sig0(a), Maj(a, b, c))
        h, g, f = g, f, e
        e = add32(d, T1)
        d, c, b = c, b, a
        a = add32(T1, T2)
    return [a, b, c, d, e, f, g, h]


def sha256_t_rounds_mid(msg, t_mid):
    """Return intermediate state after t_mid rounds (for filtering)."""
    return sha256_t_rounds(msg, t_mid)


def state_xor(s1, s2):
    """XOR two 8-word states."""
    return [(a ^ b) & MASK for a, b in zip(s1, s2)]


def state_hw(state):
    """Total Hamming weight of 8-word state."""
    return sum(hw(w) for w in state)


def random_msg():
    """Random 16-word message block."""
    return [random.getrandbits(32) for _ in range(16)]


def gen_affine_subspace(base_msg, bit_positions):
    """
    Generate 2^k messages from affine subspace.
    bit_positions: list of k bit positions in W[0] to flip.
    Returns list of 2^k messages.
    """
    k = len(bit_positions)
    msgs = []
    for mask_val in range(1 << k):
        msg = list(base_msg)
        w0 = msg[0]
        for j in range(k):
            if mask_val & (1 << j):
                w0 ^= (1 << bit_positions[j])
        msg[0] = w0 & MASK
        msgs.append(msg)
    return msgs


def xor_sum_hashes(msgs, t):
    """Compute XOR-sum of all hashes in the subspace."""
    acc = [0] * 8
    for m in msgs:
        h = sha256_t_rounds(m, t)
        for i in range(8):
            acc[i] ^= h[i]
    return acc


def xor_sum_is_zero(xsum):
    return all(w == 0 for w in xsum)


def xor_sum_hw(xsum):
    return sum(hw(w) for w in xsum)


# ---------------------------------------------------------------------------
# Tool 1: Filtered Higher-Order Distinguisher
# ---------------------------------------------------------------------------
def tool1_filtered_distinguisher():
    print("=" * 72)
    print("TOOL 1: Filtered Higher-Order Distinguisher")
    print("=" * 72)
    print()
    print("Test: does filtering subspaces by intermediate state similarity")
    print("boost the higher-order differential advantage?")
    print()

    # Reduced sample counts for speed
    k_values = [3, 4, 5]
    t_values = [6, 7, 8, 9, 10]

    # Adaptive sample counts: fewer for larger k (more expensive)
    samples_by_k = {3: 3000, 4: 1500, 5: 600}
    filter_pct = 0.25  # keep top 25% with lowest HW spread at midpoint

    results = {}

    for t in t_values:
        t_mid = t // 2
        for k in k_values:
            n_samples = samples_by_k[k]
            bit_pos = list(range(31, 31 - k, -1))  # high bits of W[0]

            # Collect XOR-sums and filter scores
            xor_sums = []
            filter_scores = []  # lower = more similar intermediate states

            for trial in range(n_samples):
                base = random_msg()
                msgs = gen_affine_subspace(base, bit_pos)

                # Compute XOR-sum of t-round hashes
                xsum = xor_sum_hashes(msgs, t)

                # Compute filter score: max HW of pairwise XOR at midpoint
                mid_states = [sha256_t_rounds_mid(m, t_mid) for m in msgs]
                # Use average pairwise HW difference (sample a few pairs)
                hw_diffs = []
                n_m = len(mid_states)
                # Sample up to 8 pairs for speed
                pair_count = min(8, n_m * (n_m - 1) // 2)
                for p in range(pair_count):
                    i1 = random.randint(0, n_m - 1)
                    i2 = random.randint(0, n_m - 2)
                    if i2 >= i1:
                        i2 += 1
                    diff = state_xor(mid_states[i1], mid_states[i2])
                    hw_diffs.append(state_hw(diff))
                avg_hw = sum(hw_diffs) / max(len(hw_diffs), 1)

                xor_sums.append(xsum)
                filter_scores.append(avg_hw)

            # Split into filtered (low HW spread) and unfiltered
            threshold_idx = int(n_samples * filter_pct)
            sorted_indices = sorted(range(n_samples), key=lambda i: filter_scores[i])
            filtered_idx = set(sorted_indices[:threshold_idx])

            # Measure bias: expected HW of XOR-sum for random = 128 per word
            # Actual for higher-order diff: may deviate (zeros = perfect)
            filt_zero = 0
            filt_total = 0
            filt_hw_sum = 0
            unfilt_zero = 0
            unfilt_total = 0
            unfilt_hw_sum = 0

            for i in range(n_samples):
                xhw = xor_sum_hw(xor_sums[i])
                is_zero = xor_sum_is_zero(xor_sums[i])
                if i in filtered_idx:
                    filt_total += 1
                    filt_hw_sum += xhw
                    if is_zero:
                        filt_zero += 1
                else:
                    unfilt_total += 1
                    unfilt_hw_sum += xhw
                    if is_zero:
                        unfilt_zero += 1

            # Expected HW for random XOR-sum = 128.0
            filt_avg_hw = filt_hw_sum / max(filt_total, 1)
            unfilt_avg_hw = unfilt_hw_sum / max(unfilt_total, 1)
            filt_bias = abs(128.0 - filt_avg_hw)
            unfilt_bias = abs(128.0 - unfilt_avg_hw)

            filt_zero_rate = filt_zero / max(filt_total, 1)
            unfilt_zero_rate = unfilt_zero / max(unfilt_total, 1)

            boost = filt_bias / max(unfilt_bias, 0.01)

            results[(t, k)] = {
                'filt_avg_hw': filt_avg_hw,
                'unfilt_avg_hw': unfilt_avg_hw,
                'filt_bias': filt_bias,
                'unfilt_bias': unfilt_bias,
                'boost': boost,
                'filt_zero_rate': filt_zero_rate,
                'unfilt_zero_rate': unfilt_zero_rate,
            }

    # Print results table
    print(f"{'t':>3} {'k':>3} {'Filt HW':>9} {'Unfilt HW':>10} "
          f"{'Filt Bias':>10} {'Unfilt Bias':>12} {'Boost':>7} "
          f"{'Filt 0-rate':>12} {'Unfilt 0-rate':>14}")
    print("-" * 100)
    for t in t_values:
        for k in k_values:
            r = results[(t, k)]
            print(f"{t:>3} {k:>3} {r['filt_avg_hw']:>9.2f} {r['unfilt_avg_hw']:>10.2f} "
                  f"{r['filt_bias']:>10.3f} {r['unfilt_bias']:>12.3f} {r['boost']:>7.2f}x "
                  f"{r['filt_zero_rate']:>12.6f} {r['unfilt_zero_rate']:>14.6f}")
        print()

    return results


# ---------------------------------------------------------------------------
# Tool 2: Adaptive Bit Selection for Subspace
# ---------------------------------------------------------------------------
def tool2_adaptive_bits():
    print("=" * 72)
    print("TOOL 2: Adaptive Bit Selection for Affine Subspace")
    print("=" * 72)
    print()
    print("Find which k-bit subsets of W[0] give the strongest distinguisher.")
    print()

    t_values = [7, 8, 9]
    results = {}

    for t in t_values:
        print(f"--- t = {t} rounds ---")

        # k=2: test ALL C(32,2) = 496 subsets
        k = 2
        n_trials = 40
        best_bias_k2 = 0
        best_bits_k2 = None
        all_biases_k2 = []

        for bits in combinations(range(32), k):
            bias_sum = 0
            for _ in range(n_trials):
                base = random_msg()
                msgs = gen_affine_subspace(base, list(bits))
                xsum = xor_sum_hashes(msgs, t)
                xhw = xor_sum_hw(xsum)
                bias_sum += abs(128.0 - xhw)
            avg_bias = bias_sum / n_trials
            all_biases_k2.append((bits, avg_bias))
            if avg_bias > best_bias_k2:
                best_bias_k2 = avg_bias
                best_bits_k2 = bits

        all_biases_k2.sort(key=lambda x: -x[1])
        print(f"  k={k}: Tested all {len(all_biases_k2)} subsets")
        print(f"    Top 5 bit subsets (bias from 128):")
        for bits, bias in all_biases_k2[:5]:
            print(f"      bits={str(bits):<15s} avg_bias={bias:.3f}")
        print(f"    Worst 3:")
        for bits, bias in all_biases_k2[-3:]:
            print(f"      bits={str(bits):<15s} avg_bias={bias:.3f}")
        results[(t, k)] = {'best_bits': best_bits_k2, 'best_bias': best_bias_k2,
                           'top5': all_biases_k2[:5]}

        # k=3: sample 300 subsets from C(32,3) = 4960 (exhaustive too slow)
        k = 3
        n_trials = 30
        n_subsets_k3 = 300
        best_bias_k3 = 0
        best_bits_k3 = None
        all_biases_k3 = []

        all_combos_k3 = list(combinations(range(32), k))
        sampled_k3 = random.sample(all_combos_k3, min(n_subsets_k3, len(all_combos_k3)))
        for bits in sampled_k3:
            bias_sum = 0
            for _ in range(n_trials):
                base = random_msg()
                msgs = gen_affine_subspace(base, list(bits))
                xsum = xor_sum_hashes(msgs, t)
                xhw = xor_sum_hw(xsum)
                bias_sum += abs(128.0 - xhw)
            avg_bias = bias_sum / n_trials
            all_biases_k3.append((bits, avg_bias))
            if avg_bias > best_bias_k3:
                best_bias_k3 = avg_bias
                best_bits_k3 = bits

        all_biases_k3.sort(key=lambda x: -x[1])
        print(f"  k={k}: Sampled {len(all_biases_k3)} of {len(all_combos_k3)} subsets")
        print(f"    Top 5 bit subsets:")
        for bits, bias in all_biases_k3[:5]:
            print(f"      bits={str(bits):<20s} avg_bias={bias:.3f}")
        print(f"    Worst 3:")
        for bits, bias in all_biases_k3[-3:]:
            print(f"      bits={str(bits):<20s} avg_bias={bias:.3f}")
        results[(t, k)] = {'best_bits': best_bits_k3, 'best_bias': best_bias_k3,
                           'top5': all_biases_k3[:5]}

        # k=4,5: sample 200 random subsets
        for k in [4, 5]:
            n_trials = 40 if k == 4 else 25
            n_subsets = 200
            best_bias = 0
            best_bits = None
            all_biases = []

            seen = set()
            for _ in range(n_subsets):
                bits = tuple(sorted(random.sample(range(32), k)))
                if bits in seen:
                    continue
                seen.add(bits)
                bias_sum = 0
                for _ in range(n_trials):
                    base = random_msg()
                    msgs = gen_affine_subspace(base, list(bits))
                    xsum = xor_sum_hashes(msgs, t)
                    xhw = xor_sum_hw(xsum)
                    bias_sum += abs(128.0 - xhw)
                avg_bias = bias_sum / n_trials
                all_biases.append((bits, avg_bias))
                if avg_bias > best_bias:
                    best_bias = avg_bias
                    best_bits = bits

            all_biases.sort(key=lambda x: -x[1])
            print(f"  k={k}: Sampled {len(all_biases)} random subsets")
            print(f"    Top 5 bit subsets:")
            for bits, bias in all_biases[:5]:
                print(f"      bits={str(bits):<25s} avg_bias={bias:.3f}")
            results[(t, k)] = {'best_bits': best_bits, 'best_bias': best_bias,
                               'top5': all_biases[:5]}
        print()

    return results


# ---------------------------------------------------------------------------
# Tool 3: Iterated Higher-Order Statistics (Chi-Squared)
# ---------------------------------------------------------------------------
def tool3_iterated_statistics():
    print("=" * 72)
    print("TOOL 3: Iterated Higher-Order Statistics (Chi-Squared)")
    print("=" * 72)
    print()
    print("Collect N XOR-sums, test per-bit distribution for non-uniformity.")
    print()

    k = 4
    bit_pos = [31, 30, 29, 28]  # high bits of W[0]
    t_values = [7, 8, 9, 10]

    # For each t, incrementally collect XOR-sums and check when chi-squared
    # detects non-uniformity (p < 0.01)
    max_N = 3000
    check_interval = 100

    results = {}

    for t in t_values:
        print(f"--- t = {t} rounds, k = {k} ---")

        # Track bit counts across XOR-sum collection
        # We look at all 256 bits of the XOR-sum
        bit_counts = [0] * 256  # count of 1-bits at each position
        total_trials = 0

        detection_N = None  # first N where p < 0.01
        best_chi2 = 0
        best_bit = -1

        for trial_batch in range(max_N // check_interval):
            for _ in range(check_interval):
                base = random_msg()
                msgs = gen_affine_subspace(base, bit_pos)
                xsum = xor_sum_hashes(msgs, t)
                total_trials += 1

                # Update bit counts
                for w_idx in range(8):
                    word = xsum[w_idx]
                    for b in range(32):
                        if word & (1 << b):
                            bit_counts[w_idx * 32 + b] += 1

            # Chi-squared test on each bit position
            # Under null (uniform): each bit is Bernoulli(0.5)
            # Expected count of 1s = total_trials / 2
            # Chi-squared = sum((obs - exp)^2 / exp) per bit
            N = total_trials
            expected = N / 2.0
            if expected < 5:
                continue

            # Find the most biased bit
            max_chi2_this = 0
            max_bit_this = -1
            significant_bits = 0

            for bit_idx in range(256):
                obs_1 = bit_counts[bit_idx]
                obs_0 = N - obs_1
                chi2 = ((obs_1 - expected) ** 2 + (obs_0 - expected) ** 2) / expected
                if chi2 > max_chi2_this:
                    max_chi2_this = chi2
                    max_bit_this = bit_idx
                # chi2 > 6.635 corresponds to p < 0.01 for 1 df
                if chi2 > 6.635:
                    significant_bits += 1

            # Bonferroni correction: for 256 tests, need p < 0.01/256
            # chi2 > 17.37 for individual test to pass Bonferroni
            bonferroni_threshold = 17.37

            if max_chi2_this > best_chi2:
                best_chi2 = max_chi2_this
                best_bit = max_bit_this

            if detection_N is None and max_chi2_this > bonferroni_threshold:
                detection_N = N

            if N % 500 == 0 or trial_batch == max_N // check_interval - 1:
                print(f"  N={N:>5d}: max_chi2={max_chi2_this:>8.2f} "
                      f"(bit {max_bit_this:>3d}) sig_bits={significant_bits:>3d}/256 "
                      f"{'*** DETECTED ***' if max_chi2_this > bonferroni_threshold else ''}")

        if detection_N:
            sha_evals = detection_N * (1 << k)
            print(f"  => Detection at N={detection_N} subspaces "
                  f"({sha_evals} SHA-256 evals, {math.log2(sha_evals):.1f} bits)")
        else:
            print(f"  => NOT detected within N={max_N} subspaces")

        results[t] = {
            'detection_N': detection_N,
            'best_chi2': best_chi2,
            'best_bit': best_bit,
            'total_N': total_trials,
        }
        print()

    return results


# ---------------------------------------------------------------------------
# Tool 4: Combined Weapon
# ---------------------------------------------------------------------------
def tool4_combined_weapon(tool2_results, tool3_results):
    print("=" * 72)
    print("TOOL 4: Combined Weapon — Minimum Evaluations for 99% Confidence")
    print("=" * 72)
    print()
    print("Uses best bits from Tool 2 + filtering from Tool 1 + iterated stats.")
    print()

    t_values = [7, 8, 9, 10]
    k = 4
    filter_pct = 0.25  # keep top 25%
    n_trials = 2000

    results = {}

    for t in t_values:
        print(f"--- t = {t} rounds ---")

        # Get best bits from Tool 2 if available, else use high bits
        best_bits = None
        for kk in [4, 3, 2]:
            key = (t, kk)
            if key in tool2_results and tool2_results[key]['best_bits'] is not None:
                best_bits = list(tool2_results[key]['best_bits'])
                if len(best_bits) < k:
                    # Pad with adjacent high bits
                    extras = [b for b in range(31, -1, -1) if b not in best_bits]
                    best_bits = best_bits + extras[:k - len(best_bits)]
                elif len(best_bits) > k:
                    best_bits = best_bits[:k]
                break
        if best_bits is None:
            best_bits = [31, 30, 29, 28]

        print(f"  Using bit positions: {best_bits}")

        t_mid = t // 2

        # Phase 1: Collect XOR-sums with filter scores
        bit_counts_filt = [0] * 256
        bit_counts_unfilt = [0] * 256
        n_filt = 0
        n_unfilt = 0

        batch_data = []  # (xsum, filter_score)

        for trial in range(n_trials):
            base = random_msg()
            msgs = gen_affine_subspace(base, best_bits)
            xsum = xor_sum_hashes(msgs, t)

            # Filter score: average pairwise HW at midpoint
            mid_states = [sha256_t_rounds_mid(m, t_mid) for m in msgs]
            hw_diffs = []
            n_m = len(mid_states)
            for p in range(min(6, n_m)):
                i1 = p
                i2 = (p + 1) % n_m
                diff = state_xor(mid_states[i1], mid_states[i2])
                hw_diffs.append(state_hw(diff))
            avg_hw = sum(hw_diffs) / max(len(hw_diffs), 1)
            batch_data.append((xsum, avg_hw))

        # Sort by filter score, take top 25%
        batch_data.sort(key=lambda x: x[1])
        threshold_idx = int(len(batch_data) * filter_pct)

        for i, (xsum, _) in enumerate(batch_data):
            target = bit_counts_filt if i < threshold_idx else bit_counts_unfilt
            if i < threshold_idx:
                n_filt += 1
            else:
                n_unfilt += 1
            for w_idx in range(8):
                word = xsum[w_idx]
                for b in range(32):
                    if word & (1 << b):
                        target[w_idx * 32 + b] += 1

        # Chi-squared analysis for filtered set
        exp_filt = n_filt / 2.0
        exp_unfilt = n_unfilt / 2.0

        max_chi2_filt = 0
        max_bit_filt = -1
        sig_bits_filt = 0
        max_chi2_unfilt = 0
        sig_bits_unfilt = 0

        bonferroni = 17.37

        for bit_idx in range(256):
            # Filtered
            obs1 = bit_counts_filt[bit_idx]
            obs0 = n_filt - obs1
            if exp_filt > 0:
                chi2_f = ((obs1 - exp_filt) ** 2 + (obs0 - exp_filt) ** 2) / exp_filt
            else:
                chi2_f = 0
            if chi2_f > max_chi2_filt:
                max_chi2_filt = chi2_f
                max_bit_filt = bit_idx
            if chi2_f > 6.635:
                sig_bits_filt += 1

            # Unfiltered
            obs1u = bit_counts_unfilt[bit_idx]
            obs0u = n_unfilt - obs1u
            if exp_unfilt > 0:
                chi2_u = ((obs1u - exp_unfilt) ** 2 + (obs0u - exp_unfilt) ** 2) / exp_unfilt
            else:
                chi2_u = 0
            if chi2_u > max_chi2_unfilt:
                max_chi2_unfilt = chi2_u
            if chi2_u > 6.635:
                sig_bits_unfilt += 1

        # Estimate minimum evaluations for 99% confidence
        # If we detected bias, extrapolate from chi2 growth rate
        if max_chi2_filt > bonferroni:
            # Detected! Compute evaluations used
            sha_evals = n_filt * (1 << k) + n_trials * (1 << k)  # filtering overhead
            # The filter overhead: we evaluated all n_trials subspaces at midpoint too
            total_evals = n_trials * (1 << k) * 2  # full + midpoint
            status = "DETECTED"
        elif max_chi2_filt > 3.84:
            # Marginal — extrapolate
            # chi2 scales linearly with N, need chi2 > bonferroni
            scale_factor = bonferroni / max_chi2_filt
            est_N = int(n_filt * scale_factor * 1.5)  # safety margin
            total_evals = int(est_N / filter_pct) * (1 << k) * 2
            status = "MARGINAL (extrapolated)"
            sha_evals = total_evals
        else:
            # Not detected — report lower bound
            sha_evals = -1
            total_evals = -1
            status = "NOT DETECTED"

        print(f"  Filtered set: N={n_filt}, max_chi2={max_chi2_filt:.2f} "
              f"(bit {max_bit_filt}), sig_bits={sig_bits_filt}/256")
        print(f"  Unfiltered:   N={n_unfilt}, max_chi2={max_chi2_unfilt:.2f}, "
              f"sig_bits={sig_bits_unfilt}/256")
        if sha_evals > 0:
            print(f"  => {status}: ~{total_evals} SHA-256 evals "
                  f"({math.log2(total_evals):.1f} bits)")
        else:
            print(f"  => {status} within budget")

        # Boost from filtering
        boost = max_chi2_filt / max(max_chi2_unfilt, 0.01)
        print(f"  Filtering chi2 boost: {boost:.2f}x")

        results[t] = {
            'status': status,
            'sha_evals': total_evals,
            'max_chi2_filt': max_chi2_filt,
            'max_chi2_unfilt': max_chi2_unfilt,
            'boost': boost,
            'sig_bits_filt': sig_bits_filt,
            'best_bits': best_bits,
        }
        print()

    # Summary table
    print("=" * 72)
    print("COMBINED WEAPON SUMMARY")
    print("=" * 72)
    print(f"{'Rounds':>7} {'Status':>20} {'SHA-256 Evals':>16} {'Log2':>8} "
          f"{'Chi2 (filt)':>12} {'Boost':>7} {'Sig Bits':>10}")
    print("-" * 82)
    for t in t_values:
        r = results[t]
        evals_str = f"{r['sha_evals']:,}" if r['sha_evals'] > 0 else "N/A"
        log2_str = f"{math.log2(r['sha_evals']):.1f}" if r['sha_evals'] > 0 else "N/A"
        print(f"{t:>7} {r['status']:>20} {evals_str:>16} {log2_str:>8} "
              f"{r['max_chi2_filt']:>12.2f} {r['boost']:>7.2f}x {r['sig_bits_filt']:>10}")

    return results


# ---------------------------------------------------------------------------
# Main — run with timing budget
# ---------------------------------------------------------------------------
def main():
    print("*" * 72)
    print("*  WEAPON: Higher-Order Differential + Filtering (HODF)")
    print("*  Combined SHA-256 Cryptanalysis Attack")
    print("*" * 72)
    print()

    t0 = time.time()

    # Tool 1
    t1_start = time.time()
    tool1_results = tool1_filtered_distinguisher()
    t1_time = time.time() - t1_start
    print(f"[Tool 1 completed in {t1_time:.1f}s]\n")

    # Tool 2 — most expensive, limit scope if running behind
    elapsed = time.time() - t0
    remaining = 110 - elapsed  # target 110s to leave margin
    t2_start = time.time()
    if remaining > 40:
        tool2_results = tool2_adaptive_bits()
    else:
        print("TOOL 2: SKIPPED (time budget exceeded)")
        tool2_results = {}
    t2_time = time.time() - t2_start
    print(f"[Tool 2 completed in {t2_time:.1f}s]\n")

    # Tool 3
    elapsed = time.time() - t0
    remaining = 110 - elapsed
    t3_start = time.time()
    if remaining > 15:
        tool3_results = tool3_iterated_statistics()
    else:
        print("TOOL 3: SKIPPED (time budget exceeded)")
        tool3_results = {}
    t3_time = time.time() - t3_start
    print(f"[Tool 3 completed in {t3_time:.1f}s]\n")

    # Tool 4: Combined weapon
    elapsed = time.time() - t0
    remaining = 110 - elapsed
    t4_start = time.time()
    if remaining > 10:
        tool4_results = tool4_combined_weapon(tool2_results, tool3_results)
    else:
        print("TOOL 4: SKIPPED (time budget exceeded)")
        tool4_results = {}
    t4_time = time.time() - t4_start
    print(f"[Tool 4 completed in {t4_time:.1f}s]\n")

    total = time.time() - t0
    print(f"\nTotal runtime: {total:.1f}s")


if __name__ == "__main__":
    main()
