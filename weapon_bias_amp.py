#!/usr/bin/env python3
"""
SHA-256 Bias Amplification Weapon — Exhaustive condition search, linear and
nonlinear approximation, and attack-complexity estimation for reduced rounds.

Tool 1: Exhaustive condition search on W[0] (mod, byte, bit-pattern, XOR-HW)
Tool 2: Linear approximation search (Matsui Algorithm 1 style)
Tool 3: Nonlinear (quadratic) approximation search
Tool 4: Attack complexity estimation and security curve

Uses real SHA-256 constants and IV throughout.
"""

import numpy as np
import time
import sys
from collections import defaultdict

# ============================================================================
# SHA-256 constants
# ============================================================================
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

MASK32 = 0xFFFFFFFF

# ============================================================================
# Core SHA-256 primitives
# ============================================================================

def rotr(x, n):
    return ((x >> n) | (x << (32 - n))) & MASK32

def shr(x, n):
    return x >> n

def hw(x):
    return bin(x).count('1')

def expand_message_schedule(W0_15, num_rounds):
    W = list(W0_15[:16])
    for i in range(16, num_rounds):
        s0 = rotr(W[i-15], 7) ^ rotr(W[i-15], 18) ^ shr(W[i-15], 3)
        s1 = rotr(W[i-2], 17) ^ rotr(W[i-2], 19) ^ shr(W[i-2], 10)
        W.append((W[i-16] + s0 + W[i-7] + s1) & MASK32)
    return W

def sha256_t_rounds(msg_words, t):
    """Run t rounds of SHA-256 from standard IV. Return 8-word state."""
    W = expand_message_schedule(msg_words, max(t, 16))
    a, b, c, d, e, f, g, h = IV
    for i in range(t):
        S1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25)
        ch = (e & f) ^ ((~e) & MASK32 & g)
        temp1 = (h + S1 + ch + K[i] + W[i]) & MASK32
        S0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22)
        maj = (a & b) ^ (a & c) ^ (b & c)
        temp2 = (S0 + maj) & MASK32
        h = g; g = f; f = e
        e = (d + temp1) & MASK32
        d = c; c = b; b = a
        a = (temp1 + temp2) & MASK32
    return [a, b, c, d, e, f, g, h]

def random_msg(rng):
    return [int(rng.integers(0, 2**32)) for _ in range(16)]

def get_e_bits(state):
    """Extract all 32 bits of e-register (state[4])."""
    e = state[4]
    return [(e >> b) & 1 for b in range(32)]

def get_bit(val, b):
    return (val >> b) & 1

# ============================================================================
# TOOL 1: Exhaustive condition search
# ============================================================================

def tool1_condition_search():
    print("=" * 72)
    print("TOOL 1: Exhaustive Condition Search on W[0]")
    print("=" * 72)
    t0 = time.time()

    rng = np.random.default_rng(42)
    N_TOTAL = 500_000

    # Pre-generate all messages
    print(f"  Generating {N_TOTAL} random messages...")
    all_msgs = []
    for _ in range(N_TOTAL):
        all_msgs.append(random_msg(rng))

    results_by_rounds = {}

    for t_rounds in [4, 6]:
        print(f"\n--- {t_rounds}-round SHA-256 ---")

        # Compute all outputs (e-register bits)
        print(f"  Computing {t_rounds}-round outputs for {N_TOTAL} messages...")
        w0_vals = np.array([m[0] for m in all_msgs], dtype=np.uint64)
        w1_vals = np.array([m[1] for m in all_msgs], dtype=np.uint64)

        # e-register bits for each sample: shape (N, 32)
        e_bits = np.zeros((N_TOTAL, 32), dtype=np.int8)
        for idx in range(N_TOTAL):
            st = sha256_t_rounds(all_msgs[idx], t_rounds)
            e = st[4]
            for b in range(32):
                e_bits[idx, b] = (e >> b) & 1

        # Baseline bias per bit
        baseline_bias = np.abs(e_bits.mean(axis=0) - 0.5)
        best_base_bit = int(np.argmax(baseline_bias))
        print(f"  Baseline max |bias| = {baseline_bias.max():.6f} at e-bit {best_base_bit}")

        # --- Build conditions ---
        conditions = {}

        # W[0] mod 4
        for r in range(4):
            mask = (w0_vals % 4) == r
            conditions[f"W0 mod 4 == {r}"] = mask

        # W[0] & 0xFF == random byte values (50 samples)
        byte_samples = rng.integers(0, 256, size=50)
        for bv in byte_samples:
            mask = (w0_vals & 0xFF) == int(bv)
            conditions[f"W0&0xFF == 0x{int(bv):02x}"] = mask

        # W[0] has specific 2-bit patterns at Sigma1-aligned positions
        # Sigma1 uses rotations by 6, 11, 25 — check 2-bit patterns at those offsets
        for shift in [6, 11, 25]:
            for pat in range(4):
                mask = ((w0_vals >> shift) & 3) == pat
                conditions[f"W0 bits[{shift+1}:{shift}] == {pat:02b}"] = mask

        # W[0] HW ranges
        w0_hw = np.array([hw(int(v)) for v in w0_vals])
        for lo, hi, label in [(0, 8, "low HW 0-8"), (8, 16, "mid-low HW 8-16"),
                               (16, 24, "mid-hi HW 16-24"), (24, 33, "high HW 24-32"),
                               (0, 4, "very low HW 0-4"), (28, 33, "very high HW 28-32")]:
            mask = (w0_hw >= lo) & (w0_hw < hi)
            conditions[f"W0 {label}"] = mask

        # W[0] XOR W[1] has low HW
        xor_hw = np.array([hw(int(w0_vals[i] ^ w1_vals[i])) for i in range(N_TOTAL)])
        for thr in [4, 6, 8, 10]:
            mask = xor_hw <= thr
            conditions[f"W0^W1 HW <= {thr}"] = mask

        # W[0] mod 8, mod 16
        for m in [8, 16]:
            for r in range(m):
                mask = (w0_vals % m) == r
                conditions[f"W0 mod {m} == {r}"] = mask

        # --- Evaluate each condition ---
        print(f"  Testing {len(conditions)} conditions...")
        cond_results = []
        for name, mask in conditions.items():
            count = int(mask.sum())
            if count < 50:
                continue
            sub = e_bits[mask]
            bias_per_bit = np.abs(sub.mean(axis=0) - 0.5)
            best_bit = int(np.argmax(bias_per_bit))
            best_bias = float(bias_per_bit[best_bit])
            amplification = best_bias / max(baseline_bias[best_bit], 1e-9)
            cond_results.append((name, count, best_bit, best_bias, amplification))

        cond_results.sort(key=lambda x: -x[3])

        print(f"\n  Top 15 single conditions (by max |bias| in any e-bit):")
        print(f"  {'Condition':<38} {'N':>7} {'Bit':>4} {'|Bias|':>10} {'Amp':>7}")
        print(f"  {'-'*38} {'-'*7} {'-'*4} {'-'*10} {'-'*7}")
        for name, cnt, bit, bias, amp in cond_results[:15]:
            print(f"  {name:<38} {cnt:>7} {bit:>4} {bias:>10.6f} {amp:>7.2f}x")

        # --- Best pair of conditions (top 20 x top 20) ---
        print(f"\n  Searching best condition pairs (top 20 x top 20)...")
        top_names = [c[0] for c in cond_results[:20]]
        best_pair = None
        best_pair_bias = 0.0

        for i in range(len(top_names)):
            for j in range(i + 1, len(top_names)):
                n1, n2 = top_names[i], top_names[j]
                joint_mask = conditions[n1] & conditions[n2]
                cnt = int(joint_mask.sum())
                if cnt < 30:
                    continue
                sub = e_bits[joint_mask]
                bias_per_bit = np.abs(sub.mean(axis=0) - 0.5)
                best_bit = int(np.argmax(bias_per_bit))
                best_bias = float(bias_per_bit[best_bit])
                if best_bias > best_pair_bias:
                    best_pair_bias = best_bias
                    best_pair = (n1, n2, cnt, best_bit, best_bias)

        if best_pair:
            n1, n2, cnt, bit, bias = best_pair
            amp = bias / max(baseline_bias[bit], 1e-9)
            print(f"  Best pair: [{n1}] AND [{n2}]")
            print(f"    N={cnt}, bit={bit}, |bias|={bias:.6f}, amp={amp:.2f}x")
        else:
            print(f"  No valid pairs found.")

        results_by_rounds[t_rounds] = {
            'best_single': cond_results[0] if cond_results else None,
            'best_pair': best_pair,
            'baseline_max': float(baseline_bias.max()),
        }

    elapsed = time.time() - t0
    print(f"\n  Tool 1 completed in {elapsed:.1f}s")
    return results_by_rounds


# ============================================================================
# TOOL 2: Linear Approximation Search (Matsui's Algorithm 1 for SHA-256)
# ============================================================================

def tool2_linear_approximation():
    print("\n" + "=" * 72)
    print("TOOL 2: Linear Approximation Search (Matsui Algorithm 1)")
    print("=" * 72)
    t0 = time.time()

    rng = np.random.default_rng(123)
    N_MASKS = 10_000      # random mask pairs to test
    N_SAMPLES = 20_000    # samples per mask pair

    # Pre-generate messages for reuse
    print(f"  Pre-generating {N_SAMPLES} messages...")
    msgs = [random_msg(rng) for _ in range(N_SAMPLES)]

    results_by_rounds = {}

    for t_rounds in [4, 5, 6, 7, 8]:
        print(f"\n--- {t_rounds}-round linear search ---")

        # Pre-compute all outputs
        print(f"  Computing {t_rounds}-round outputs...")
        outputs = []
        for m in msgs:
            outputs.append(sha256_t_rounds(m, t_rounds))

        # We focus on single-bit masks for efficiency:
        # input_mask: bit position in W[0] (0..31)
        # output_mask: bit position in e-register = state[4] (0..31)
        # Also try a-register = state[0]
        #
        # Correlation = 2 * Pr[input_bit XOR output_bit = 0] - 1

        best_bias = 0.0
        best_info = None

        # Test all 32x32 = 1024 single-bit pairs for W[0]->e and W[0]->a
        for reg_name, reg_idx in [("e", 4), ("a", 0)]:
            for in_bit in range(32):
                for out_bit in range(32):
                    agree = 0
                    for idx in range(N_SAMPLES):
                        ib = (msgs[idx][0] >> in_bit) & 1
                        ob = (outputs[idx][reg_idx] >> out_bit) & 1
                        agree += (ib ^ ob ^ 1)  # 1 if equal
                    corr = abs(2.0 * agree / N_SAMPLES - 1.0)
                    if corr > best_bias:
                        best_bias = corr
                        best_info = (f"W0[{in_bit}]", f"{reg_name}[{out_bit}]", corr)

        # Also test W[1]->e, W[1]->a
        for reg_name, reg_idx in [("e", 4), ("a", 0)]:
            for in_bit in range(32):
                for out_bit in range(32):
                    agree = 0
                    for idx in range(N_SAMPLES):
                        ib = (msgs[idx][1] >> in_bit) & 1
                        ob = (outputs[idx][reg_idx] >> out_bit) & 1
                        agree += (ib ^ ob ^ 1)
                    corr = abs(2.0 * agree / N_SAMPLES - 1.0)
                    if corr > best_bias:
                        best_bias = corr
                        best_info = (f"W1[{in_bit}]", f"{reg_name}[{out_bit}]", corr)

        # Multi-bit input masks: XOR of 2 input bits -> 1 output bit
        # Sample random pairs
        n_multi = min(N_MASKS, 5000)
        for _ in range(n_multi):
            w_idx1 = int(rng.integers(0, 2))
            b1 = int(rng.integers(0, 32))
            w_idx2 = int(rng.integers(0, 2))
            b2 = int(rng.integers(0, 32))
            reg_idx = int(rng.choice([0, 4]))
            out_bit = int(rng.integers(0, 32))
            reg_name = "a" if reg_idx == 0 else "e"

            agree = 0
            for idx in range(N_SAMPLES):
                ib = ((msgs[idx][w_idx1] >> b1) & 1) ^ ((msgs[idx][w_idx2] >> b2) & 1)
                ob = (outputs[idx][reg_idx] >> out_bit) & 1
                agree += (ib ^ ob ^ 1)
            corr = abs(2.0 * agree / N_SAMPLES - 1.0)
            if corr > best_bias:
                best_bias = corr
                best_info = (f"W{w_idx1}[{b1}]^W{w_idx2}[{b2}]",
                             f"{reg_name}[{out_bit}]", corr)

        epsilon = best_bias / 2.0  # bias eps = (corr)/2 since corr = 2*eps
        print(f"  Best linear approx: {best_info[0]} -> {best_info[1]}")
        print(f"    correlation = {best_info[2]:.6f}, epsilon = {epsilon:.6f}")
        print(f"    Distinguisher queries ~ {1.0/max(epsilon**2, 1e-30):.1f}")

        results_by_rounds[t_rounds] = {
            'input': best_info[0],
            'output': best_info[1],
            'correlation': best_info[2],
            'epsilon': epsilon,
        }

    elapsed = time.time() - t0
    print(f"\n  Tool 2 completed in {elapsed:.1f}s")
    return results_by_rounds


# ============================================================================
# TOOL 3: Nonlinear (Quadratic) Approximation Search
# ============================================================================

def tool3_nonlinear_approximation():
    print("\n" + "=" * 72)
    print("TOOL 3: Nonlinear (Quadratic) Approximation Search")
    print("=" * 72)
    t0 = time.time()

    rng = np.random.default_rng(456)
    N_TRIPLES = 5_000
    N_SAMPLES = 20_000

    # Pre-generate messages
    print(f"  Pre-generating {N_SAMPLES} messages...")
    msgs = [random_msg(rng) for _ in range(N_SAMPLES)]

    results_by_rounds = {}

    for t_rounds in [4, 5, 6]:
        print(f"\n--- {t_rounds}-round quadratic search ---")

        # Pre-compute outputs
        outputs = []
        for m in msgs:
            outputs.append(sha256_t_rounds(m, t_rounds))

        # Extract W[0] bits and output bits as arrays for speed
        w0_bits = np.zeros((N_SAMPLES, 32), dtype=np.int8)
        e_bits = np.zeros((N_SAMPLES, 32), dtype=np.int8)
        a_bits = np.zeros((N_SAMPLES, 32), dtype=np.int8)
        for idx in range(N_SAMPLES):
            w0 = msgs[idx][0]
            for b in range(32):
                w0_bits[idx, b] = (w0 >> b) & 1
            e_bits[idx] = np.array([(outputs[idx][4] >> b) & 1 for b in range(32)], dtype=np.int8)
            a_bits[idx] = np.array([(outputs[idx][0] >> b) & 1 for b in range(32)], dtype=np.int8)

        best_quad_bias = 0.0
        best_quad_info = None
        best_lin_bias = 0.0
        best_lin_info = None

        # First find best linear for comparison
        for in_b in range(32):
            for reg_name, out_arr in [("e", e_bits), ("a", a_bits)]:
                for out_b in range(32):
                    corr = abs(float(np.mean((-1.0) ** (w0_bits[:, in_b] ^ out_arr[:, out_b]))))
                    if corr > best_lin_bias:
                        best_lin_bias = corr
                        best_lin_info = (f"W0[{in_b}]", f"{reg_name}[{out_b}]", corr)

        # Quadratic: (W0[i] AND W0[j]) XOR output_bit k
        # Sample random triples
        tested = 0
        for _ in range(N_TRIPLES):
            i = int(rng.integers(0, 32))
            j = int(rng.integers(0, 32))
            if i == j:
                continue
            k = int(rng.integers(0, 32))
            reg_idx = int(rng.choice([0, 1]))  # 0 = e_bits, 1 = a_bits
            reg_name = "e" if reg_idx == 0 else "a"
            out_arr = e_bits if reg_idx == 0 else a_bits

            quad_input = w0_bits[:, i] & w0_bits[:, j]
            corr = abs(float(np.mean((-1.0) ** (quad_input ^ out_arr[:, k]))))
            tested += 1

            if corr > best_quad_bias:
                best_quad_bias = corr
                best_quad_info = (f"W0[{i}]&W0[{j}]", f"{reg_name}[{k}]", corr)

        # Also try: (W0[i] AND W0[j]) XOR W0[m] -> output bit k
        for _ in range(N_TRIPLES):
            i = int(rng.integers(0, 32))
            j = int(rng.integers(0, 32))
            m = int(rng.integers(0, 32))
            if i == j:
                continue
            k = int(rng.integers(0, 32))
            reg_idx = int(rng.choice([0, 1]))
            reg_name = "e" if reg_idx == 0 else "a"
            out_arr = e_bits if reg_idx == 0 else a_bits

            quad_input = (w0_bits[:, i] & w0_bits[:, j]) ^ w0_bits[:, m]
            corr = abs(float(np.mean((-1.0) ** (quad_input ^ out_arr[:, k]))))

            if corr > best_quad_bias:
                best_quad_bias = corr
                best_quad_info = (f"W0[{i}]&W0[{j}]^W0[{m}]", f"{reg_name}[{k}]", corr)

        quad_eps = best_quad_bias / 2.0
        lin_eps = best_lin_bias / 2.0

        print(f"  Best linear:    {best_lin_info[0]} -> {best_lin_info[1]}, corr={best_lin_info[2]:.6f}, eps={lin_eps:.6f}")
        print(f"  Best quadratic: {best_quad_info[0]} -> {best_quad_info[1]}, corr={best_quad_info[2]:.6f}, eps={quad_eps:.6f}")
        improvement = best_quad_bias / max(best_lin_bias, 1e-30)
        print(f"  Quadratic/Linear ratio: {improvement:.3f}x")

        results_by_rounds[t_rounds] = {
            'best_linear': best_lin_info,
            'best_quadratic': best_quad_info,
            'lin_eps': lin_eps,
            'quad_eps': quad_eps,
            'improvement': improvement,
        }

    elapsed = time.time() - t0
    print(f"\n  Tool 3 completed in {elapsed:.1f}s")
    return results_by_rounds


# ============================================================================
# TOOL 4: Attack Complexity Estimation and Security Curve
# ============================================================================

def tool4_complexity_estimation(tool1_res, tool2_res, tool3_res):
    print("\n" + "=" * 72)
    print("TOOL 4: Attack Complexity Estimation & Security Curve")
    print("=" * 72)

    import math

    # Collect best biases per round count from all tools
    print("\n--- Consolidating best biases across all tools ---\n")

    all_rounds = sorted(set(
        list(tool1_res.keys()) + list(tool2_res.keys()) + list(tool3_res.keys())
    ))

    round_summary = {}

    for t in all_rounds:
        biases = []
        sources = []

        # Tool 1: condition-based bias
        if t in tool1_res:
            r = tool1_res[t]
            if r['best_single']:
                name, cnt, bit, bias, amp = r['best_single']
                biases.append(bias)
                sources.append(f"Cond-single({name}, bit={bit})")
            if r['best_pair']:
                n1, n2, cnt, bit, bias = r['best_pair']
                biases.append(bias)
                sources.append(f"Cond-pair({n1}&{n2}, bit={bit})")

        # Tool 2: linear approximation
        if t in tool2_res:
            r = tool2_res[t]
            biases.append(r['epsilon'])
            sources.append(f"Linear({r['input']}->{r['output']})")

        # Tool 3: quadratic approximation
        if t in tool3_res:
            r = tool3_res[t]
            biases.append(r['quad_eps'])
            sources.append(f"Quadratic({r['best_quadratic'][0]}->{r['best_quadratic'][1]})")
            biases.append(r['lin_eps'])
            sources.append(f"Linear-T3({r['best_linear'][0]}->{r['best_linear'][1]})")

        if biases:
            best_idx = int(np.argmax(biases))
            best_eps = biases[best_idx]
            best_src = sources[best_idx]
        else:
            best_eps = 0.0
            best_src = "none"

        # Queries needed for distinguisher: O(1/eps^2)
        if best_eps > 1e-12:
            queries = 1.0 / (best_eps ** 2)
            log2_queries = math.log2(queries)
        else:
            queries = float('inf')
            log2_queries = float('inf')

        round_summary[t] = {
            'best_eps': best_eps,
            'best_source': best_src,
            'queries': queries,
            'log2_queries': log2_queries,
            'all_biases': list(zip(biases, sources)),
        }

    # --- Print per-round summary ---
    print(f"  {'Rounds':>6} {'Best eps':>12} {'Queries':>14} {'log2(Q)':>10} Source")
    print(f"  {'-'*6} {'-'*12} {'-'*14} {'-'*10} {'-'*40}")
    for t in all_rounds:
        s = round_summary[t]
        if s['log2_queries'] < 1e10:
            print(f"  {t:>6} {s['best_eps']:>12.6f} {s['queries']:>14.0f} {s['log2_queries']:>10.1f} {s['best_source']}")
        else:
            print(f"  {t:>6} {s['best_eps']:>12.6f} {'inf':>14} {'inf':>10} {s['best_source']}")

    # --- Chaining independent biases ---
    print("\n--- Bias chaining analysis (independent output bits) ---\n")
    for t in all_rounds:
        s = round_summary[t]
        all_b = s['all_biases']
        if len(all_b) < 2:
            continue

        # Sort by bias descending
        all_b_sorted = sorted(all_b, key=lambda x: -x[0])
        # Take top biases that come from different output bits/sources
        print(f"  {t}-round: {len(all_b_sorted)} available approximations")
        # Combined advantage if independent: multiply correlations -> product
        # distinguisher advantage = product of correlations (piling-up)
        product = 1.0
        for eps, src in all_b_sorted[:4]:
            corr = 2.0 * eps  # correlation
            product *= corr
            print(f"    corr={corr:.6f} from {src}")
        combined_eps = product / 2.0
        if combined_eps > 1e-30:
            combined_queries = 1.0 / (combined_eps ** 2)
            combined_log2 = math.log2(max(combined_queries, 1))
        else:
            combined_log2 = float('inf')
        print(f"    => Combined eps={combined_eps:.9f}, log2(queries)={combined_log2:.1f}")
        print(f"       (assuming independent channels, piling-up lemma)")

    # --- Security curve ---
    print("\n--- Security Curve ---\n")
    print("  At what round does attack complexity exceed key thresholds?\n")
    thresholds = [(32, "2^32"), (64, "2^64"), (80, "2^80"), (128, "2^128")]

    for thr_log2, thr_label in thresholds:
        exceeded = None
        for t in all_rounds:
            if round_summary[t]['log2_queries'] >= thr_log2:
                exceeded = t
                break
        if exceeded:
            print(f"  Complexity >= {thr_label}: first at {exceeded} rounds "
                  f"(actual log2(Q) = {round_summary[exceeded]['log2_queries']:.1f})")
        else:
            print(f"  Complexity >= {thr_label}: NOT reached within tested rounds ({all_rounds[-1]} max)")

    # --- ASCII security curve plot ---
    print("\n  Security curve (log2 queries vs rounds):\n")
    max_display = 140
    for t in all_rounds:
        lq = round_summary[t]['log2_queries']
        if lq == float('inf'):
            bar_len = max_display
            label = "inf"
        else:
            bar_len = min(int(lq), max_display)
            label = f"{lq:.1f}"
        bar = "#" * max(bar_len // 2, 1)
        print(f"  t={t:>2}: {bar} [{label}]")

    # Markers
    print(f"\n  Legend: each '#' ~ 2 bits of security")
    for thr_log2, thr_label in thresholds:
        pos = thr_log2 // 2
        print(f"  {thr_label} threshold at position {pos}")

    return round_summary


# ============================================================================
# MAIN
# ============================================================================

def main():
    overall_t0 = time.time()
    print("SHA-256 Bias Amplification Weapon")
    print("Exhaustive condition search + Linear/Nonlinear approx + Complexity")
    print(f"{'=' * 72}\n")

    # Tool 1
    tool1_res = tool1_condition_search()

    # Tool 2
    tool2_res = tool2_linear_approximation()

    # Tool 3
    tool3_res = tool3_nonlinear_approximation()

    # Tool 4
    tool4_res = tool4_complexity_estimation(tool1_res, tool2_res, tool3_res)

    total = time.time() - overall_t0
    print(f"\n{'=' * 72}")
    print(f"TOTAL RUNTIME: {total:.1f}s")
    print(f"{'=' * 72}")


if __name__ == "__main__":
    main()
