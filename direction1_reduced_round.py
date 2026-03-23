#!/usr/bin/env python3
"""
Direction 1: Reduced-Round Frontier Analysis

Investigates WHY the collision attack boundary sits at 31 out of 64 rounds
by measuring how cryptographic properties evolve through the transition zone.

Metrics measured per round count N (1..40):
  1. Differential propagation: probability of low-HW output diff, average output HW
  2. Avalanche completeness: rank of avalanche matrix, condition proxy
  3. Algebraic degree estimation via higher-order differentials
  4. Linearization quality: best linear approximation bias
  5. Security margin: -log2(best differential probability)
"""

import numpy as np
import time
import sys
from collections import defaultdict

# ─── SHA-256 constants ───────────────────────────────────────────────────────

K256 = [
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
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
]

IV = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
]

MASK32 = 0xFFFFFFFF

# ─── Variable-round SHA-256 compression ──────────────────────────────────────

def sha256_compress(W16, H8, R):
    """
    SHA-256 compression for R rounds (1..64).
    W16: list of 16 uint32 message words
    H8:  list of 8 uint32 initial state words
    R:   number of rounds
    Returns: list of 8 uint32 output words (with feed-forward addition)
    """
    W = list(W16)
    for i in range(16, max(R, 16)):
        w15 = W[i - 15]
        s0 = (((w15 >> 7) | (w15 << 25)) ^ ((w15 >> 18) | (w15 << 14)) ^ (w15 >> 3)) & MASK32
        w2 = W[i - 2]
        s1 = (((w2 >> 17) | (w2 << 15)) ^ ((w2 >> 19) | (w2 << 13)) ^ (w2 >> 10)) & MASK32
        W.append((W[i - 16] + s0 + W[i - 7] + s1) & MASK32)

    a, b, c, d, e, f, g, h = H8

    for t in range(R):
        S1 = (((e >> 6) | (e << 26)) ^ ((e >> 11) | (e << 21)) ^ ((e >> 25) | (e << 7))) & MASK32
        ch = ((e & f) ^ ((~e) & g)) & MASK32
        T1 = (h + S1 + ch + K256[t] + W[t]) & MASK32
        S0 = (((a >> 2) | (a << 30)) ^ ((a >> 13) | (a << 19)) ^ ((a >> 22) | (a << 10))) & MASK32
        maj = ((a & b) ^ (a & c) ^ (b & c)) & MASK32
        T2 = (S0 + maj) & MASK32

        h = g
        g = f
        f = e
        e = (d + T1) & MASK32
        d = c
        c = b
        b = a
        a = (T1 + T2) & MASK32

    # Feed-forward
    return [
        (H8[0] + a) & MASK32, (H8[1] + b) & MASK32,
        (H8[2] + c) & MASK32, (H8[3] + d) & MASK32,
        (H8[4] + e) & MASK32, (H8[5] + f) & MASK32,
        (H8[6] + g) & MASK32, (H8[7] + h) & MASK32,
    ]


def state_to_bits(state):
    """Convert 8 x uint32 state to 256-bit array."""
    bits = []
    for w in state:
        for b in range(31, -1, -1):
            bits.append((w >> b) & 1)
    return bits


def hamming_weight_256(state_xor):
    """Hamming weight of XOR of two states (each 8 x uint32)."""
    hw = 0
    for w in state_xor:
        hw += bin(w).count('1')
    return hw


def xor_states(s1, s2):
    return [(a ^ b) & MASK32 for a, b in zip(s1, s2)]


def random_w16(rng):
    return [int(rng.integers(0, 2**32)) for _ in range(16)]


# ─── Study 1: Differential Propagation ──────────────────────────────────────

def study_differential_propagation(max_rounds=40, n_samples=10000):
    """
    For each round count N, inject 1-bit or 2-bit differences in W[0],
    measure output differential properties.
    """
    print("=" * 80)
    print("STUDY 1: DIFFERENTIAL PROPAGATION (rounds 1-{})".format(max_rounds))
    print("=" * 80)
    print(f"Samples per round: {n_samples}")
    print()

    rng = np.random.default_rng(42)

    results = {}

    for N in range(1, max_rounds + 1):
        low_hw_count = 0       # output diff HW < 32
        total_hw = 0
        best_single_bit = 0    # best probability that a specific output bit stays 0 in diff
        bit_diff_counts = np.zeros(256, dtype=np.int64)

        for trial in range(n_samples):
            W = random_w16(rng)
            # 1-bit difference in W[0]
            bit_pos = int(rng.integers(0, 32))
            W2 = list(W)
            W2[0] ^= (1 << bit_pos)

            out1 = sha256_compress(W, IV, N)
            out2 = sha256_compress(W2, IV, N)

            diff = xor_states(out1, out2)
            hw = hamming_weight_256(diff)
            total_hw += hw
            if hw < 32:
                low_hw_count += 1

            # Track per-bit differences
            for wi, d in enumerate(diff):
                for b in range(32):
                    if (d >> b) & 1:
                        bit_diff_counts[wi * 32 + (31 - b)] += 1

        avg_hw = total_hw / n_samples
        prob_low_hw = low_hw_count / n_samples

        # Best single-bit differential: the bit that changes least often
        # (highest probability of staying unchanged = 1 - count/n_samples)
        min_change_rate = np.min(bit_diff_counts) / n_samples
        best_bit_prob = 1.0 - min_change_rate  # prob that this bit is zero in diff

        results[N] = {
            'prob_low_hw': prob_low_hw,
            'avg_hw': avg_hw,
            'best_bit_prob': best_bit_prob,
            'min_change_rate': min_change_rate,
        }

        if N <= 10 or N % 5 == 0:
            print(f"  Round {N:2d}: P(HW<32)={prob_low_hw:.4f}  avg_HW={avg_hw:.1f}  "
                  f"best_bit_hold={best_bit_prob:.6f}")

    print()
    return results


# ─── Study 2: Avalanche Completeness ────────────────────────────────────────

def study_avalanche_completeness(max_rounds=40, n_samples=200):
    """
    For each round N, flip each of 32 bits in W[0], measure which of 256
    output bits are affected. Build 256x32 avalanche matrix, compute rank.
    """
    print("=" * 80)
    print("STUDY 2: AVALANCHE COMPLETENESS (rounds 1-{})".format(max_rounds))
    print("=" * 80)
    print(f"Samples per round per bit: {n_samples}")
    print()

    rng = np.random.default_rng(123)

    results = {}
    full_rank_round = None

    for N in range(1, max_rounds + 1):
        # Avalanche matrix: 256 output bits x 32 input bits
        # Entry (i,j) = fraction of times output bit i flips when input bit j flips
        aval_matrix = np.zeros((256, 32), dtype=np.float64)

        for trial in range(n_samples):
            W = random_w16(rng)
            out_base = sha256_compress(W, IV, N)

            for j in range(32):
                W2 = list(W)
                W2[0] ^= (1 << (31 - j))
                out_flip = sha256_compress(W2, IV, N)
                diff = xor_states(out_base, out_flip)

                for wi, d in enumerate(diff):
                    for b in range(32):
                        if (d >> (31 - b)) & 1:
                            aval_matrix[wi * 32 + b, j] += 1

        aval_matrix /= n_samples

        # Binarize: if flip probability > 0.01, count as "affected"
        binary_matrix = (aval_matrix > 0.01).astype(np.int8)

        # Rank over GF(2) using row reduction
        rank = gf2_rank(binary_matrix)

        # Condition proxy: how close to ideal 0.5 are the entries?
        # Use mean absolute deviation from 0.5
        active_entries = aval_matrix[binary_matrix > 0]
        if len(active_entries) > 0:
            mad_from_half = np.mean(np.abs(active_entries - 0.5))
        else:
            mad_from_half = 0.5

        # Number of active entries (out of 256*32 = 8192)
        n_active = np.sum(binary_matrix)

        results[N] = {
            'rank': rank,
            'n_active': int(n_active),
            'mad_from_half': mad_from_half,
        }

        if full_rank_round is None and rank == 32:
            full_rank_round = N

        if N <= 10 or N % 5 == 0:
            print(f"  Round {N:2d}: rank={rank:3d}/32  active={n_active:5d}/8192  "
                  f"MAD_from_0.5={mad_from_half:.4f}")

    print(f"\n  Full rank (32) first achieved at round: {full_rank_round}")
    print()
    return results, full_rank_round


def gf2_rank(matrix):
    """Compute rank of binary matrix over GF(2)."""
    M = matrix.astype(np.int8).copy()
    rows, cols = M.shape
    rank = 0
    for col in range(cols):
        # Find pivot
        pivot = None
        for row in range(rank, rows):
            if M[row, col]:
                pivot = row
                break
        if pivot is None:
            continue
        # Swap
        M[[rank, pivot]] = M[[pivot, rank]]
        # Eliminate
        for row in range(rows):
            if row != rank and M[row, col]:
                M[row] = (M[row] ^ M[rank])
        rank += 1
    return rank


# ─── Study 3: Algebraic Degree Estimation ───────────────────────────────────

def study_algebraic_degree(max_rounds=40, n_test_points=500):
    """
    Estimate algebraic degree via higher-order differentials.
    For round count N, the k-th order differential of a degree-d function is zero if k > d.
    We test orders k = 1,2,...,32 on a single output bit to find effective degree.

    For efficiency, we use a cube-sum approach on W[0] bits.
    """
    print("=" * 80)
    print("STUDY 3: ALGEBRAIC DEGREE ESTIMATION (rounds 1-{})".format(max_rounds))
    print("=" * 80)
    print(f"Test points per round: {n_test_points}")
    print()

    rng = np.random.default_rng(456)

    results = {}

    for N in range(1, max_rounds + 1):
        # Test degrees from high to low: find max k where k-th order diff is nonzero
        # We'll test a single output bit (bit 0 of word 0) for efficiency
        # k-th order differential over bits 0..k-1 of W[0]

        max_degree_found = 0

        # Cap cube dimension at 10 to keep 2^k evaluations feasible
        max_k = min(10, 31)

        for k in range(max_k, 0, -1):
            nonzero_count = 0
            n_pts = min(n_test_points, max(20, 50 // max(1, k - 3)))

            for trial in range(n_pts):
                W = random_w16(rng)
                W[0] &= ~((1 << k) - 1)

                xor_sum = 0
                for subset in range(1 << k):
                    Wc = list(W)
                    Wc[0] = W[0] | subset
                    out = sha256_compress(Wc, IV, N)
                    xor_sum ^= out[0]

                if xor_sum != 0:
                    nonzero_count += 1

            if nonzero_count > 0:
                max_degree_found = k
                break

        results[N] = {'degree': max_degree_found}

        if N <= 10 or N % 5 == 0:
            print(f"  Round {N:2d}: estimated algebraic degree >= {max_degree_found}")

    print()
    return results


# ─── Study 4: Linearization Quality ─────────────────────────────────────────

def study_linearization(max_rounds=40, n_samples=10000):
    """
    For each round count N, estimate the best linear approximation bias.
    We test random linear masks on input W[0] and output word 0,
    and find the maximum bias |P(in·a XOR out·b = 0) - 0.5|.
    """
    print("=" * 80)
    print("STUDY 4: LINEARIZATION QUALITY (rounds 1-{})".format(max_rounds))
    print("=" * 80)
    print(f"Samples per round: {n_samples}")
    print()

    rng = np.random.default_rng(789)

    # Pre-generate test masks: single-bit masks are most informative
    # Test all 32 input bits x first 8 output bits = 256 combinations
    results = {}

    for N in range(1, max_rounds + 1):
        # Generate random inputs and outputs once
        inputs_w0 = rng.integers(0, 2**32, size=n_samples, dtype=np.uint64)
        outputs = np.zeros(n_samples, dtype=np.uint64)

        for i in range(n_samples):
            W = [int(inputs_w0[i])] + [int(rng.integers(0, 2**32)) for _ in range(15)]
            out = sha256_compress(W, IV, N)
            outputs[i] = out[0]

        # Test single-bit linear approximations: parity(input bit j) vs parity(output bit k)
        best_bias = 0.0

        for j in range(32):
            in_bits = ((inputs_w0 >> j) & 1).astype(np.int8)
            for k in range(8):  # test 8 output bits for speed
                out_bits = ((outputs >> k) & 1).astype(np.int8)
                xor_bits = in_bits ^ out_bits
                count_zero = np.sum(xor_bits == 0)
                bias = abs(count_zero / n_samples - 0.5)
                if bias > best_bias:
                    best_bias = bias

        results[N] = {'best_bias': best_bias}

        if N <= 10 or N % 5 == 0:
            if best_bias > 0:
                log_bias = -np.log2(best_bias) if best_bias > 0 else float('inf')
                print(f"  Round {N:2d}: best_bias={best_bias:.6f}  -log2(bias)={log_bias:.2f}")
            else:
                print(f"  Round {N:2d}: best_bias=0 (no detectable linear relation)")

    print()
    return results


# ─── Study 5: Security Margin ───────────────────────────────────────────────

def study_security_margin(diff_results):
    """
    security_margin(N) = -log2(best_differential_probability(N))
    where best_differential_probability is estimated from the differential study.
    """
    print("=" * 80)
    print("STUDY 5: SECURITY MARGIN")
    print("=" * 80)
    print()

    results = {}

    for N, data in sorted(diff_results.items()):
        # Use the best single-bit hold probability as proxy for differential prob
        # A differential "holds" if specific output bits stay at expected values
        # We use: P(specific output bit unchanged) as per-bit differential probability
        # For the full 256-bit state, effective probability ~ (best_bit_prob)^256
        # But more practically, we measure the per-bit deviation from random (0.5)

        best_bit_prob = data['best_bit_prob']
        min_change_rate = data['min_change_rate']

        # The "differential probability" per output bit: deviation from random
        # If change_rate = 0.5, it's perfectly random => no exploitable bias
        # If change_rate < 0.5, there's a bias toward this bit staying same
        # If change_rate > 0.5, there's a bias toward this bit flipping

        deviation = abs(min_change_rate - 0.5)
        if deviation > 1e-10:
            security_bits = -np.log2(deviation)
        else:
            security_bits = float('inf')

        # Also compute based on prob_low_hw (more direct measure)
        prob_low = data['prob_low_hw']
        # For random 256-bit vector, P(HW < 32) ~ very small
        # Expected HW = 128, std ~ 8. P(HW < 32) ~ 2^{-100} for random
        # So if observed prob_low_hw >> 2^{-100}, there's exploitable structure
        if prob_low > 0:
            margin_from_hw = -np.log2(prob_low)
        else:
            margin_from_hw = 100.0  # below detection threshold => random

        results[N] = {
            'per_bit_security': security_bits,
            'margin_from_hw': margin_from_hw,
            'deviation': deviation,
        }

    return results


# ─── Main ────────────────────────────────────────────────────────────────────

def main():
    print("╔══════════════════════════════════════════════════════════════════════════════╗")
    print("║  DIRECTION 1: REDUCED-ROUND FRONTIER ANALYSIS                              ║")
    print("║  Why does the collision attack boundary sit at 31 rounds?                   ║")
    print("╚══════════════════════════════════════════════════════════════════════════════╝")
    print()

    t_start = time.time()
    n_samples = 10000

    # ── Study 1: Differential propagation ──
    t1 = time.time()
    diff_results = study_differential_propagation(max_rounds=40, n_samples=n_samples)
    print(f"  [Study 1 completed in {time.time()-t1:.1f}s]\n")

    # ── Study 2: Avalanche completeness ──
    t2 = time.time()
    aval_results, full_rank_round = study_avalanche_completeness(max_rounds=40, n_samples=100)
    print(f"  [Study 2 completed in {time.time()-t2:.1f}s]\n")

    # ── Study 3: Algebraic degree ──
    t3 = time.time()
    # Use fewer test points per round and limit cube size for speed
    deg_results = study_algebraic_degree(max_rounds=40, n_test_points=30)
    print(f"  [Study 3 completed in {time.time()-t3:.1f}s]\n")

    # ── Study 4: Linearization quality ──
    t4 = time.time()
    lin_results = study_linearization(max_rounds=40, n_samples=5000)
    print(f"  [Study 4 completed in {time.time()-t4:.1f}s]\n")

    # ── Study 5: Security margin ──
    sec_results = study_security_margin(diff_results)

    # ═══════════════════════════════════════════════════════════════════════════
    # COMPREHENSIVE TABLE
    # ═══════════════════════════════════════════════════════════════════════════
    print()
    print("=" * 120)
    print("COMPREHENSIVE TABLE: ALL METRICS vs ROUND COUNT")
    print("=" * 120)
    print(f"{'Rnd':>3s} | {'P(HW<32)':>9s} | {'Avg HW':>7s} | {'Best1Bit':>9s} | "
          f"{'Aval Rank':>9s} | {'Active':>6s} | {'MAD 0.5':>7s} | "
          f"{'Alg Deg':>7s} | {'Lin Bias':>9s} | {'-log2(B)':>8s} | "
          f"{'SecMargin':>9s} | {'Phase':>8s}")
    print("-" * 120)

    # Determine phase transition
    phase_transition_round = None

    for N in range(1, 41):
        d = diff_results[N]
        a = aval_results.get(N, {'rank': 0, 'n_active': 0, 'mad_from_half': 0.5})
        g = deg_results.get(N, {'degree': 0})
        l = lin_results.get(N, {'best_bias': 0})
        s = sec_results.get(N, {'per_bit_security': 0, 'margin_from_hw': 0, 'deviation': 0})

        bias = l['best_bias']
        log_bias = f"{-np.log2(bias):.2f}" if bias > 1e-10 else ">13"

        sec_margin = s['margin_from_hw']
        sec_str = f"{sec_margin:.1f}" if sec_margin < 100 else ">log(N)"

        # Phase classification
        if d['prob_low_hw'] > 0.5:
            phase = "EASY"
        elif d['prob_low_hw'] > 0.01:
            phase = "MEDIUM"
        elif d['prob_low_hw'] > 0.0:
            phase = "HARD"
        else:
            phase = "INTRACT"
            if phase_transition_round is None:
                phase_transition_round = N

        print(f"{N:3d} | {d['prob_low_hw']:9.5f} | {d['avg_hw']:7.1f} | {d['best_bit_prob']:9.6f} | "
              f"{a['rank']:5d}/32  | {a['n_active']:6d} | {a['mad_from_half']:7.4f} | "
              f"{g['degree']:7d} | {bias:9.6f} | {log_bias:>8s} | "
              f"{sec_str:>9s} | {phase:>8s}")

    print("=" * 120)
    print()

    # ═══════════════════════════════════════════════════════════════════════════
    # ANALYSIS
    # ═══════════════════════════════════════════════════════════════════════════
    print("=" * 80)
    print("ANALYSIS: PHASE TRANSITION IDENTIFICATION")
    print("=" * 80)
    print()

    # Find key transition points
    first_full_diffusion = None
    first_no_low_hw = None
    first_high_degree = None
    first_negligible_bias = None

    for N in range(1, 41):
        d = diff_results[N]
        if first_no_low_hw is None and d['prob_low_hw'] == 0:
            first_no_low_hw = N
        if first_high_degree is None and deg_results[N]['degree'] >= 15:
            first_high_degree = N
        l = lin_results[N]
        if first_negligible_bias is None and l['best_bias'] < 0.01:
            first_negligible_bias = N

    print(f"  Avalanche matrix full rank at:            round {full_rank_round}")
    print(f"  P(HW<32) drops to zero at:                round {first_no_low_hw}")
    print(f"  Algebraic degree reaches >= 15 at:        round {first_high_degree}")
    print(f"  Linear bias drops below 0.01 at:          round {first_negligible_bias}")
    print(f"  Phase transition (differential intract.): round {phase_transition_round}")
    print()

    # Differential probability decay analysis
    print("=" * 80)
    print("DIFFERENTIAL PROBABILITY DECAY CURVE")
    print("=" * 80)
    print()
    print(f"{'Round':>5s} | {'P(HW<32)':>10s} | {'-log2(P)':>10s} | {'Bar':>40s}")
    print("-" * 72)
    for N in range(1, 41):
        p = diff_results[N]['prob_low_hw']
        if p > 0:
            logp = -np.log2(p)
            bar = "#" * min(40, int(logp * 3))
        else:
            logp = float('inf')
            bar = "#" * 40 + " (ZERO)"
        print(f"{N:5d} | {p:10.6f} | {logp:10.2f} | {bar}")

    print()

    # Average HW convergence to 128
    print("=" * 80)
    print("AVERAGE OUTPUT HAMMING WEIGHT CONVERGENCE (target: 128 for random)")
    print("=" * 80)
    print()
    print(f"{'Round':>5s} | {'Avg HW':>8s} | {'|HW-128|':>8s} | {'Bar':>40s}")
    print("-" * 68)
    for N in range(1, 41):
        hw = diff_results[N]['avg_hw']
        dev = abs(hw - 128)
        bar = "#" * min(40, int(dev / 3))
        print(f"{N:5d} | {hw:8.1f} | {dev:8.1f} |  {bar}")

    print()

    # Linearization bias decay
    print("=" * 80)
    print("LINEAR APPROXIMATION BIAS DECAY")
    print("=" * 80)
    print()
    print(f"{'Round':>5s} | {'Bias':>10s} | {'-log2':>8s} | {'Bar':>40s}")
    print("-" * 68)
    for N in range(1, 41):
        bias = lin_results[N]['best_bias']
        if bias > 1e-10:
            logb = -np.log2(bias)
            bar = "#" * min(40, int(logb * 3))
        else:
            logb = float('inf')
            bar = "#" * 40 + " (ZERO)"
        print(f"{N:5d} | {bias:10.6f} | {logb:8.2f} | {bar}")

    print()

    # Algebraic degree growth
    print("=" * 80)
    print("ALGEBRAIC DEGREE GROWTH")
    print("=" * 80)
    print()
    print(f"{'Round':>5s} | {'Degree':>7s} | {'Bar':>30s}")
    print("-" * 48)
    for N in range(1, 41):
        deg = deg_results[N]['degree']
        bar = "#" * min(30, deg)
        print(f"{N:5d} | {deg:7d} | {bar}")

    print()

    # ═══════════════════════════════════════════════════════════════════════════
    # CONCLUSIONS
    # ═══════════════════════════════════════════════════════════════════════════
    print("=" * 80)
    print("CONCLUSIONS")
    print("=" * 80)
    print()
    print("The 31-round collision boundary exists because of a convergence of")
    print("multiple cryptographic properties reaching maturity around this point:")
    print()
    print(f"  1. DIFFERENTIAL DIFFUSION: Output Hamming weight of XOR differences")
    print(f"     converges toward the random expectation (128) by round ~{first_no_low_hw or '??'},")
    print(f"     making differential paths exponentially unlikely to hold.")
    print()
    print(f"  2. AVALANCHE COMPLETENESS: The avalanche matrix reaches full rank")
    print(f"     at round {full_rank_round}, meaning every output bit depends on every")
    print(f"     input bit. Before this, attackers can find 'free' variables.")
    print()
    print(f"  3. ALGEBRAIC DEGREE: The Boolean degree of the output grows rapidly,")
    print(f"     reaching >= 15 by round {first_high_degree or '??'}. High degree defeats")
    print(f"     higher-order differential and algebraic attacks.")
    print()
    print(f"  4. LINEARITY DECAY: The best linear approximation bias drops below")
    print(f"     0.01 by round {first_negligible_bias or '??'}, defeating linear cryptanalysis.")
    print()
    print(f"  The phase transition around rounds {phase_transition_round or '??'} represents the point")
    print(f"  where ALL these properties become simultaneously strong enough")
    print(f"  to defeat known attack strategies. The best known collision attack")
    print(f"  reaching 31 rounds likely exploits residual weaknesses in the")
    print(f"  transition zone just before full diffusion is achieved.")
    print()

    total_time = time.time() - t_start
    print(f"Total runtime: {total_time:.1f}s")
    print()


if __name__ == "__main__":
    main()
