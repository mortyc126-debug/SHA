#!/usr/bin/env python3
"""
C20: GF(2) Jacobian Rank Deficiency x Message Schedule

Tests how the GF(2) Jacobian rank deficiency of SHA-256 compression
function grows with the number of rounds.

Background:
- The GF(2) Jacobian of SHA-256 has rank 254-255 (deficiency 1-2) at full 64 rounds
- For a random 256x512 binary matrix, P(rank < 256) ~ 2^{-257} => structural
- Question: how does deficiency grow with rounds?

Experiment:
- For R rounds (1..40), compute the 256x768 GF(2) Jacobian of
  F_R: (W[0..15], H[0..7]) -> H'[0..7]
- Use numerical differentiation: flip each input bit, XOR outputs
- Compute rank over GF(2), deficiency = 256 - rank
- Repeat at N random operating points per R
- Analyze growth pattern of deficiency vs R

OPTIMIZATION: Use message-only Jacobian (256x512) if full is too slow.
"""

import numpy as np
import time
import sys

# SHA-256 constants
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

MASK32 = 0xFFFFFFFF


def sha256_compress(W16, H8, R):
    """
    SHA-256 compression for R rounds.
    W16: list of 16 uint32 message words
    H8: list of 8 uint32 state words
    R: number of rounds (1..64)
    Returns: list of 8 uint32 output state words
    """
    # Message schedule expansion
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

    return [a, b, c, d, e, f, g, h]


def words_to_bits(words):
    """Convert list of uint32 words to bit array (MSB first per word)."""
    bits = []
    for w in words:
        for b in range(31, -1, -1):
            bits.append((w >> b) & 1)
    return bits


def flip_bit_in_words(words, bit_idx):
    """Flip bit bit_idx in the word array. bit_idx is MSB-first overall."""
    words = list(words)
    word_idx = bit_idx // 32
    bit_in_word = 31 - (bit_idx % 32)  # MSB first
    words[word_idx] ^= (1 << bit_in_word)
    return words


def compute_jacobian_message_only(W16, H8, R):
    """
    Compute the 256 x 512 GF(2) Jacobian of R-round SHA-256
    with respect to the 512 message bits only (H fixed).

    J[i, j] = (F(W ^ e_j, H)[i]) XOR (F(W, H)[i])
    """
    base_out = sha256_compress(W16, H8, R)
    base_bits = words_to_bits(base_out)

    J = np.zeros((256, 512), dtype=np.uint8)

    for j in range(512):
        W_flip = flip_bit_in_words(W16, j)
        flip_out = sha256_compress(W_flip, H8, R)
        flip_bits = words_to_bits(flip_out)
        for i in range(256):
            J[i, j] = base_bits[i] ^ flip_bits[i]

    return J


def compute_jacobian_full(W16, H8, R):
    """
    Compute the 256 x 768 GF(2) Jacobian of R-round SHA-256
    with respect to all 768 input bits (512 message + 256 state).
    """
    base_out = sha256_compress(W16, H8, R)
    base_bits = words_to_bits(base_out)

    J = np.zeros((256, 768), dtype=np.uint8)

    # Message bits (columns 0..511)
    for j in range(512):
        W_flip = flip_bit_in_words(W16, j)
        flip_out = sha256_compress(W_flip, H8, R)
        flip_bits = words_to_bits(flip_out)
        for i in range(256):
            J[i, j] = base_bits[i] ^ flip_bits[i]

    # State bits (columns 512..767)
    for j in range(256):
        H_flip = flip_bit_in_words(H8, j)
        flip_out = sha256_compress(W16, H_flip, R)
        flip_bits = words_to_bits(flip_out)
        for i in range(256):
            J[i, 512 + j] = base_bits[i] ^ flip_bits[i]

    return J


def gf2_rank(M):
    """Rank of binary matrix over GF(2) via row reduction."""
    A = M.copy().astype(np.uint8)
    rows, cols = A.shape
    rank = 0
    pivot_col = 0
    for row in range(rows):
        if pivot_col >= cols:
            break
        # Find pivot
        found = False
        for r in range(row, rows):
            if A[r, pivot_col]:
                found = True
                if r != row:
                    A[[row, r]] = A[[r, row]]
                break
        if not found:
            pivot_col += 1
            # retry this row with next column
            continue
        # Eliminate
        for r in range(rows):
            if r != row and A[r, pivot_col]:
                A[r] ^= A[row]
        rank += 1
        pivot_col += 1
    return rank


def gf2_rref_and_pivots(M):
    """
    Row reduce M over GF(2). Returns (RREF matrix, list of pivot columns).
    """
    A = M.copy().astype(np.uint8)
    rows, cols = A.shape
    pivot_cols = []
    pivot_col = 0
    current_row = 0
    for current_row_target in range(rows):
        if pivot_col >= cols:
            break
        # Search for pivot in this column
        found = False
        while pivot_col < cols:
            for r in range(current_row_target, rows):
                if A[r, pivot_col]:
                    found = True
                    if r != current_row_target:
                        A[[current_row_target, r]] = A[[r, current_row_target]]
                    break
            if found:
                break
            pivot_col += 1
        if not found:
            break
        pivot_cols.append(pivot_col)
        # Eliminate
        for r in range(rows):
            if r != current_row_target and A[r, pivot_col]:
                A[r] ^= A[current_row_target]
        pivot_col += 1

    return A, pivot_cols


def gf2_left_nullspace(M):
    """
    Compute the LEFT null space of M over GF(2).
    i.e., vectors v such that v^T M = 0 (equivalently M^T v = 0).
    M is m x n. Left null space has dimension m - rank(M).
    Returns list of m-dimensional vectors.
    """
    # Left null space of M = right null space of M^T
    Mt = M.T.copy().astype(np.uint8)
    n, m = Mt.shape  # Mt is n x m

    # Augment Mt with identity: [Mt | I_n] and row reduce
    # After RREF, rows where Mt part is zero give null space vectors from I part
    # Actually easier: row reduce M^T and find null space

    A, pivot_cols = gf2_rref_and_pivots(Mt)
    n_rows, n_cols = A.shape  # n x m

    # Free columns of Mt correspond to left null space of M
    free_cols = [c for c in range(n_cols) if c not in pivot_cols]

    nullvecs = []
    for fc in free_cols:
        v = np.zeros(n_cols, dtype=np.uint8)
        v[fc] = 1
        for idx, pc in enumerate(pivot_cols):
            if idx < n_rows:
                v[pc] = A[idx, fc]
        nullvecs.append(v)

    return nullvecs


def main():
    print("=" * 72)
    print("C20: GF(2) Jacobian Rank Deficiency x Message Schedule")
    print("=" * 72)
    print()
    print("Testing how GF(2) Jacobian rank deficiency grows with rounds.")
    print("SHA-256 compression: F_R(W, H) -> H'")
    print()

    rng = np.random.RandomState(2024)

    # Timing calibration
    print("Calibrating speed...")
    t0 = time.time()
    W_test = [rng.randint(0, 2**32) for _ in range(16)]
    H_test = [rng.randint(0, 2**32) for _ in range(8)]
    for _ in range(100):
        sha256_compress(W_test, H_test, 20)
    dt = time.time() - t0
    evals_per_sec = 100 / dt
    print(f"  Speed: {evals_per_sec:.0f} SHA evals/sec")

    # Budget: aim for ~5 minutes total
    # Message-only Jacobian: 512 evals per point
    # Full Jacobian: 768 evals per point
    # Use message-only for speed, plus a few full Jacobian checks

    # Phase 1: Message-only Jacobian (256 x 512) across rounds
    # For each R: N points x 512 evals = 512*N evals
    # With 40 round values and N=30: 40 * 30 * 512 = 614,400 evals
    # At ~15000 evals/sec for small R => ~40 seconds

    round_list = list(range(1, 21)) + [24, 28, 32, 36, 40, 48, 56, 64]
    N_points = 20  # operating points per round count

    total_evals_est = len(round_list) * N_points * 512
    est_time = total_evals_est / evals_per_sec
    print(f"  Estimated time for message-only Jacobian: {est_time:.0f}s")
    print(f"  Round values: {len(round_list)}, Points per round: {N_points}")

    if est_time > 500:
        # Reduce
        N_points = max(5, int(N_points * 400 / est_time))
        print(f"  Reduced to N={N_points} points per round to fit budget")

    print()
    print("-" * 72)
    print("Phase 1: Message-only GF(2) Jacobian (256 x 512)")
    print("  Fixing H = random, varying W. Deficiency = 256 - rank(J).")
    print("-" * 72)
    print()
    print(f"{'R':>4s}  {'mean_rank':>10s}  {'mean_def':>10s}  {'min_def':>8s}  {'max_def':>8s}  {'std_def':>8s}")
    print("-" * 60)

    results_msg = {}
    t_start = time.time()

    for R in round_list:
        ranks = []
        for n in range(N_points):
            W = [int(rng.randint(0, 2**32)) for _ in range(16)]
            H = [int(rng.randint(0, 2**32)) for _ in range(8)]

            J = compute_jacobian_message_only(W, H, R)
            r = gf2_rank(J)
            ranks.append(r)

        ranks = np.array(ranks)
        defs = 256 - ranks
        mean_r = np.mean(ranks)
        mean_d = np.mean(defs)
        min_d = np.min(defs)
        max_d = np.max(defs)
        std_d = np.std(defs)

        results_msg[R] = {
            'mean_rank': mean_r, 'mean_def': mean_d,
            'min_def': min_d, 'max_def': max_d, 'std_def': std_d,
            'ranks': ranks.tolist()
        }

        elapsed = time.time() - t_start
        print(f"{R:4d}  {mean_r:10.2f}  {mean_d:10.2f}  {min_d:8d}  {max_d:8d}  {std_d:8.2f}  [{elapsed:.1f}s]")
        sys.stdout.flush()

        # Safety: abort if taking too long
        if elapsed > 360:
            print("  [TIME LIMIT - stopping early]")
            break

    total_phase1 = time.time() - t_start
    print(f"\nPhase 1 completed in {total_phase1:.1f}s")

    # Phase 2: Full Jacobian (256 x 768) at a few key round counts
    print()
    print("-" * 72)
    print("Phase 2: Full GF(2) Jacobian (256 x 768) at selected rounds")
    print("  Includes both message and state bits.")
    print("-" * 72)
    print()

    remaining_time = 480 - total_phase1  # leave 2 min margin from 10 min
    full_rounds = [1, 4, 8, 16, 32, 64]
    N_full = 5
    full_evals = len(full_rounds) * N_full * 768
    est_full = full_evals / evals_per_sec

    if est_full > remaining_time:
        N_full = max(2, int(N_full * remaining_time / est_full))
        print(f"  Reduced to N={N_full} for full Jacobian")

    print(f"{'R':>4s}  {'mean_rank':>10s}  {'mean_def':>10s}  {'min_def':>8s}  {'max_def':>8s}")
    print("-" * 55)

    results_full = {}
    t_start2 = time.time()

    for R in full_rounds:
        ranks = []
        for n in range(N_full):
            W = [int(rng.randint(0, 2**32)) for _ in range(16)]
            H = [int(rng.randint(0, 2**32)) for _ in range(8)]

            J = compute_jacobian_full(W, H, R)
            r = gf2_rank(J)
            ranks.append(r)

        ranks = np.array(ranks)
        defs = 256 - ranks
        mean_r = np.mean(ranks)
        mean_d = np.mean(defs)
        min_d = np.min(defs)
        max_d = np.max(defs)

        results_full[R] = {
            'mean_rank': mean_r, 'mean_def': mean_d,
            'min_def': min_d, 'max_def': max_d,
            'ranks': ranks.tolist()
        }

        elapsed2 = time.time() - t_start2
        print(f"{R:4d}  {mean_r:10.2f}  {mean_d:10.2f}  {min_d:8d}  {max_d:8d}  [{elapsed2:.1f}s]")
        sys.stdout.flush()

        if elapsed2 > remaining_time:
            print("  [TIME LIMIT - stopping early]")
            break

    total_phase2 = time.time() - t_start2
    print(f"\nPhase 2 completed in {total_phase2:.1f}s")

    # Phase 3: Null space analysis at high rounds
    print()
    print("-" * 72)
    print("Phase 3: Null Space Stability Analysis")
    print("  Are the deficient directions FIXED or VARIABLE across operating points?")
    print("-" * 72)
    print()

    # Pick a round count with known deficiency
    test_R_vals = [R for R in results_msg if results_msg[R]['mean_def'] > 0.5]
    if not test_R_vals:
        test_R_vals = [max(results_msg.keys())]

    test_R = max(test_R_vals)  # Use highest round with deficiency
    print(f"Testing null space stability at R={test_R}")

    null_spaces = []
    N_null = min(10, N_points)

    for n in range(N_null):
        W = [int(rng.randint(0, 2**32)) for _ in range(16)]
        H = [int(rng.randint(0, 2**32)) for _ in range(8)]

        J = compute_jacobian_message_only(W, H, test_R)
        rank = gf2_rank(J)
        deficiency = 256 - rank

        if deficiency > 0:
            # Left null space of J = output relations that always hold
            ns = gf2_left_nullspace(J)
            null_spaces.append((n, deficiency, ns))

    if len(null_spaces) >= 2:
        print(f"  Found {len(null_spaces)} points with deficiency > 0")

        # Compare null space vectors across operating points
        # Check if the same output bit patterns appear
        all_ns_sets = []
        for idx, def_val, ns_vecs in null_spaces:
            vec_set = set()
            for v in ns_vecs:
                # Convert to tuple for hashing
                vec_set.add(tuple(v.tolist()))
            all_ns_sets.append(vec_set)
            nz_positions = []
            for v in ns_vecs:
                nz = tuple(np.where(v)[0].tolist())
                nz_positions.append(nz)
            print(f"  Point {idx}: deficiency={def_val}, "
                  f"null vectors: {len(ns_vecs)}, "
                  f"support sizes: {[len(nz) for nz in nz_positions]}")

        # Check overlap between null spaces
        if len(all_ns_sets) >= 2:
            # Check if any vectors are shared
            common = all_ns_sets[0]
            for s in all_ns_sets[1:]:
                common = common & s
            print(f"\n  Common null vectors across ALL points: {len(common)}")
            if len(common) > 0:
                print("  => FIXED null space (global structural weakness)")
                null_stability = "FIXED"
            else:
                # Check pairwise overlap
                n_pairs = 0
                n_overlap = 0
                for i in range(len(all_ns_sets)):
                    for j in range(i+1, len(all_ns_sets)):
                        overlap = len(all_ns_sets[i] & all_ns_sets[j])
                        n_pairs += 1
                        if overlap > 0:
                            n_overlap += 1
                print(f"  Pairwise overlaps: {n_overlap}/{n_pairs} pairs share vectors")
                if n_overlap == 0:
                    print("  => VARIABLE null space (per-message property)")
                    null_stability = "VARIABLE"
                else:
                    print("  => PARTIALLY FIXED null space")
                    null_stability = "PARTIAL"
    else:
        print("  Insufficient deficient points for null space comparison")
        null_stability = "UNKNOWN"

    # Phase 4: Analysis and verdict
    print()
    print("=" * 72)
    print("ANALYSIS")
    print("=" * 72)
    print()

    # Extract deficiency trend
    Rs = sorted(results_msg.keys())
    defs = [results_msg[R]['mean_def'] for R in Rs]

    print("Message-only Jacobian deficiency trend:")
    print(f"  R=1:  deficiency = {results_msg[Rs[0]]['mean_def']:.2f}")
    if len(Rs) > 1:
        mid_idx = len(Rs) // 2
        print(f"  R={Rs[mid_idx]}:  deficiency = {results_msg[Rs[mid_idx]]['mean_def']:.2f}")
    print(f"  R={Rs[-1]}: deficiency = {results_msg[Rs[-1]]['mean_def']:.2f}")

    # Check growth pattern
    early_def = np.mean([results_msg[R]['mean_def'] for R in Rs if R <= 4])
    mid_def = np.mean([results_msg[R]['mean_def'] for R in Rs if 8 <= R <= 16])
    late_def = np.mean([results_msg[R]['mean_def'] for R in Rs if R >= 32])

    print(f"\n  Early (R<=4) avg deficiency:  {early_def:.2f}")
    print(f"  Mid   (8<=R<=16) avg deficiency: {mid_def:.2f}")
    print(f"  Late  (R>=32) avg deficiency:  {late_def:.2f}")

    # Growth analysis
    if late_def > 2 and late_def > early_def + 1:
        if late_def > 2 * mid_def:
            growth = "SUPER-LINEAR"
        else:
            growth = "LINEAR"
    elif late_def > early_def + 0.5:
        growth = "SUB-LINEAR"
    elif late_def > 0.5:
        growth = "CONSTANT"
    else:
        growth = "ZERO"

    print(f"\n  Growth pattern: {growth}")

    # Linear fit on non-zero region
    valid_Rs = [R for R in Rs if results_msg[R]['mean_def'] > 0]
    if len(valid_Rs) >= 3:
        x = np.array(valid_Rs, dtype=float)
        y = np.array([results_msg[R]['mean_def'] for R in valid_Rs])
        slope, intercept = np.polyfit(x, y, 1)
        print(f"  Linear fit: deficiency ~ {slope:.4f} * R + {intercept:.2f}")

        # Log fit
        if np.all(y > 0):
            log_slope, log_intercept = np.polyfit(np.log(x), y, 1)
            print(f"  Log fit: deficiency ~ {log_slope:.4f} * ln(R) + {log_intercept:.2f}")

    # Full Jacobian comparison
    if results_full:
        print(f"\n  Full Jacobian (768 inputs) deficiencies:")
        for R in sorted(results_full.keys()):
            rf = results_full[R]
            print(f"    R={R:3d}: rank={rf['mean_rank']:.1f}, deficiency={rf['mean_def']:.1f}")

    # Null space
    print(f"\n  Null space stability: {null_stability}")

    # VERDICT
    print()
    print("=" * 72)
    print("VERDICT")
    print("=" * 72)

    max_def = max(defs) if defs else 0
    has_deficiency = max_def > 0.5
    grows = late_def > early_def + 0.5 if len(Rs) > 4 else False

    if grows and late_def >= 3:
        verdict = "ALIVE"
        detail = (f"Rank deficiency GROWS with rounds ({growth}).\n"
                  f"  From {early_def:.1f} at low rounds to {late_def:.1f} at high rounds.\n"
                  f"  This accumulation could compound into exploitable weakness.")
    elif has_deficiency and not grows:
        verdict = "ANOMALY"
        detail = (f"Rank deficiency is present (mean={max_def:.1f}) but does NOT grow.\n"
                  f"  Structural but constant — insufficient for attack advantage.\n"
                  f"  Null space: {null_stability}")
    elif has_deficiency and grows:
        verdict = "ALIVE"
        detail = (f"Rank deficiency grows from {early_def:.1f} to {late_def:.1f}.\n"
                  f"  Growth pattern: {growth}.\n"
                  f"  Null space: {null_stability}")
    else:
        verdict = "DEAD"
        detail = (f"No significant rank deficiency observed.\n"
                  f"  Max deficiency = {max_def:.2f} across all rounds.\n"
                  f"  Previous finding may have been an artifact.")

    label = {"ALIVE": "***ALIVE***", "ANOMALY": "**ANOMALY**", "DEAD": "DEAD"}
    print(f"\n  {label[verdict]}: {verdict}")
    print()
    print(f"  {detail}")

    # Random baseline comparison
    print()
    print("-" * 72)
    print("Random baseline: rank of random 256x512 GF(2) matrix")
    print("-" * 72)
    n_rand = 10
    rand_ranks = []
    for _ in range(n_rand):
        M = np.random.randint(0, 2, size=(256, 512), dtype=np.uint8)
        rand_ranks.append(gf2_rank(M))
    print(f"  {n_rand} random matrices: ranks = {rand_ranks}")
    print(f"  Mean rank = {np.mean(rand_ranks):.1f} (expected: 256.0)")
    print(f"  Deficiency = {256 - np.mean(rand_ranks):.1f} (expected: 0.0)")
    all_full = all(r == 256 for r in rand_ranks)
    print(f"  All full rank: {all_full}")

    print()
    total_time = time.time() - t_start
    print(f"Total runtime: {total_time:.1f}s")
    print()
    print(f"RESULT: {verdict}")


if __name__ == "__main__":
    main()
