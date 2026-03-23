#!/usr/bin/env python3
"""
ЗАДАНИЕ 2: Jacobian Rank Cascade — как теряется ранг системы от выхода к входу

Experiments:
A. Jacobian rank as function of backward rounds R
B. Jacobian with schedule constraint (real W[0..15])
C. Conditional rank — Wang chain δhash subspace (PRIORITY)
D. Distinguisher score vs collision distance

Priority: C > A > D > B
"""

import struct
import os
import numpy as np
from collections import defaultdict
import time

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

H0 = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
]


def rotr(x, n):
    return ((x >> n) | (x << (32 - n))) & MASK

def shr(x, n):
    return x >> n

def Sig0(x):
    return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)

def Sig1(x):
    return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)

def sig0(x):
    return rotr(x, 7) ^ rotr(x, 18) ^ shr(x, 3)

def sig1(x):
    return rotr(x, 17) ^ rotr(x, 19) ^ shr(x, 10)

def Ch(e, f, g):
    return (e & f) ^ ((~e) & g) & MASK

def Maj(a, b, c):
    return (a & b) ^ (a & c) ^ (b & c)

def add(*args):
    s = 0
    for x in args:
        s = (s + x) & MASK
    return s

def sub(a, b):
    return (a - b) & MASK


def message_schedule(M_words):
    W = list(M_words[:16])
    for i in range(16, 64):
        W.append(add(sig1(W[i-2]), W[i-7], sig0(W[i-15]), W[i-16]))
    return W


def sha256_compress(W, iv=None):
    """SHA-256 compression: returns (final_state, all_states)."""
    if iv is None:
        state = list(H0)
    else:
        state = list(iv)
    states = [tuple(state)]
    a, b, c, d, e, f, g, h = state
    for r in range(64):
        T1 = add(h, Sig1(e), Ch(e, f, g), K[r], W[r])
        T2 = add(Sig0(a), Maj(a, b, c))
        h, g, f = g, f, e
        e = add(d, T1)
        d, c, b = c, b, a
        a = add(T1, T2)
        states.append((a, b, c, d, e, f, g, h))
    return states[-1], states


def sha256_hash(M_words, iv=None):
    """Full SHA-256 hash of 16-word message, returns 8-word hash."""
    W = message_schedule(M_words)
    if iv is None:
        iv_list = list(H0)
    else:
        iv_list = list(iv)
    final, _ = sha256_compress(W, iv_list)
    return tuple(add(final[i], iv_list[i]) for i in range(8))


def forward_rounds(state_start, W_slice, start_round):
    """Run forward from state_start using W_slice starting at start_round."""
    a, b, c, d, e, f, g, h = state_start
    for i, r in enumerate(range(start_round, start_round + len(W_slice))):
        T1 = add(h, Sig1(e), Ch(e, f, g), K[r], W_slice[i])
        T2 = add(Sig0(a), Maj(a, b, c))
        h, g, f = g, f, e
        e = add(d, T1)
        d, c, b = c, b, a
        a = add(T1, T2)
    return (a, b, c, d, e, f, g, h)


def random_words(n):
    return list(struct.unpack(f'>{n}I', os.urandom(4 * n)))


def hamming_weight(x):
    return bin(x).count('1')


def hw_xor_states(s1, s2):
    return sum(hamming_weight(a ^ b) for a, b in zip(s1, s2))


def state_to_bits(state):
    bits = []
    for word in state:
        for bit in range(31, -1, -1):
            bits.append((word >> bit) & 1)
    return bits


def gf2_rank(matrix):
    """Compute rank of binary matrix over GF(2)."""
    M = np.array(matrix, dtype=np.uint8) % 2
    M = M.copy()
    rows, cols = M.shape
    pivot_row = 0
    for col in range(cols):
        found = -1
        for row in range(pivot_row, rows):
            if M[row, col] == 1:
                found = row
                break
        if found == -1:
            continue
        M[[pivot_row, found]] = M[[found, pivot_row]]
        for row in range(rows):
            if row != pivot_row and M[row, col] == 1:
                M[row] ^= M[pivot_row]
        pivot_row += 1
    return pivot_row


# =============================================================
# Wang chain generation (from methodology П-25/П-26)
# =============================================================
def generate_wang_pair():
    """
    Generate a Wang pair (Wn, Wf) with δe[2..16]=0, P=1.0.
    δW[0] = 0x8000 (bit 15).
    For rounds 1..15: δW[r] = adaptive correction.
    """
    Wn = random_words(16)
    Wf = list(Wn)
    Wf[0] = (Wn[0] ^ 0x8000) & MASK  # XOR difference in W[0]

    # Run round 0 for both
    a_n, b_n, c_n, d_n = H0[0], H0[1], H0[2], H0[3]
    e_n, f_n, g_n, h_n = H0[4], H0[5], H0[6], H0[7]
    a_f, b_f, c_f, d_f = H0[0], H0[1], H0[2], H0[3]
    e_f, f_f, g_f, h_f = H0[4], H0[5], H0[6], H0[7]

    # Round 0 with Wn[0] and Wf[0]
    T1_n = add(h_n, Sig1(e_n), Ch(e_n, f_n, g_n), K[0], Wn[0])
    T2_n = add(Sig0(a_n), Maj(a_n, b_n, c_n))
    T1_f = add(h_f, Sig1(e_f), Ch(e_f, f_f, g_f), K[0], Wf[0])
    T2_f = add(Sig0(a_f), Maj(a_f, b_f, c_f))

    h_n, g_n, f_n = g_n, f_n, e_n
    e_n = add(d_n, T1_n)
    d_n, c_n, b_n = c_n, b_n, a_n
    a_n = add(T1_n, T2_n)

    h_f, g_f, f_f = g_f, f_f, e_f
    e_f = add(d_f, T1_f)
    d_f, c_f, b_f = c_f, b_f, a_f
    a_f = add(T1_f, T2_f)

    # Rounds 1..15: adaptive Wang correction
    for r in range(1, 16):
        # Current deltas (full state difference)
        dd = sub(d_f, d_n)
        dh = sub(h_f, h_n)
        dSig1 = sub(Sig1(e_f), Sig1(e_n))
        dCh = sub(Ch(e_f, f_f, g_f), Ch(e_n, f_n, g_n))

        # δW[r] to make δe[r+1] = 0:
        # e_new_f - e_new_n = (d_f + T1_f) - (d_n + T1_n) = 0
        # => δT1 = -δd
        # δT1 = δh + δSig1 + δCh + δW[r]
        # => δW[r] = -δd - δh - δSig1 - δCh
        dW_r = sub(0, add(dd, dh, dSig1, dCh))

        Wn[r] = Wn[r]  # keep original random
        Wf[r] = add(Wn[r], dW_r)

        # Advance round r for both
        T1_n = add(h_n, Sig1(e_n), Ch(e_n, f_n, g_n), K[r], Wn[r])
        T2_n = add(Sig0(a_n), Maj(a_n, b_n, c_n))
        T1_f = add(h_f, Sig1(e_f), Ch(e_f, f_f, g_f), K[r], Wf[r])
        T2_f = add(Sig0(a_f), Maj(a_f, b_f, c_f))

        h_n, g_n, f_n = g_n, f_n, e_n
        e_n = add(d_n, T1_n)
        d_n, c_n, b_n = c_n, b_n, a_n
        a_n = add(T1_n, T2_n)

        h_f, g_f, f_f = g_f, f_f, e_f
        e_f = add(d_f, T1_f)
        d_f, c_f, b_f = c_f, b_f, a_f
        a_f = add(T1_f, T2_f)

    state_n_16 = (a_n, b_n, c_n, d_n, e_n, f_n, g_n, h_n)
    state_f_16 = (a_f, b_f, c_f, d_f, e_f, f_f, g_f, h_f)

    return Wn, Wf, state_n_16, state_f_16


def verify_wang_pair(Wn, Wf):
    """Verify that δe[2..16]=0 for a Wang pair."""
    Wn_exp = message_schedule(Wn)
    Wf_exp = message_schedule(Wf)
    _, states_n = sha256_compress(Wn_exp)
    _, states_f = sha256_compress(Wf_exp)

    zeros = 0
    for r in range(2, 17):
        de = states_n[r][4] ^ states_f[r][4]  # e-register XOR diff
        if de == 0:
            zeros += 1
    return zeros, states_n, states_f, Wn_exp, Wf_exp


# =============================================================
# EXPERIMENT C: Conditional rank — Wang δhash subspace (PRIORITY)
# =============================================================
def experiment_C(n_wang=1000, n_random=1000):
    print("=" * 70)
    print("EXPERIMENT C: Conditional rank — Wang chain δhash subspace")
    print("  (PRIORITY EXPERIMENT)")
    print("=" * 70)
    t0 = time.time()

    # Part 1: Generate Wang pairs and compute δhash
    wang_dhash_bits = []
    wang_zeros_count = []
    wang_dhash_hw = []
    n_trivial = 0

    for i in range(n_wang):
        Wn, Wf, _, _ = generate_wang_pair()

        # Full SHA-256 hash
        hash_n = sha256_hash(Wn)
        hash_f = sha256_hash(Wf)

        # δhash XOR
        dhash = tuple(a ^ b for a, b in zip(hash_n, hash_f))

        # Rule 7: check δhash ≠ 0
        if all(d == 0 for d in dhash):
            n_trivial += 1
            continue

        hw = sum(hamming_weight(d) for d in dhash)
        wang_dhash_hw.append(hw)
        wang_dhash_bits.append(state_to_bits(dhash))

        # Verify Wang chain
        zeros, _, _, _, _ = verify_wang_pair(Wn, Wf)
        wang_zeros_count.append(zeros)

    print(f"  Wang pairs generated: {n_wang}")
    print(f"  Trivial δhash=0 (collisions!): {n_trivial}")
    print(f"  Valid δhash≠0: {len(wang_dhash_bits)}")
    print(f"  Wang chain δe zeros (expected 15): mean={np.mean(wang_zeros_count):.2f}")
    print(f"  Mean HW(δhash): {np.mean(wang_dhash_hw):.2f} (expected ~128 for random)")
    print(f"  Min HW(δhash): {np.min(wang_dhash_hw)}")
    print(f"  Max HW(δhash): {np.max(wang_dhash_hw)}")

    # GF(2) rank of Wang δhash matrix
    wang_matrix = np.array(wang_dhash_bits, dtype=np.uint8)
    wang_rank = gf2_rank(wang_matrix)
    print(f"\n  *** GF(2) rank of Wang δhash matrix ({wang_matrix.shape[0]}×256): {wang_rank} ***")

    # Part 2: Null hypothesis — random pairs (Rule 12)
    print(f"\n  Null hypothesis: random single-bit-flip pairs")
    random_dhash_bits = []
    random_dhash_hw = []

    for i in range(n_random):
        M = random_words(16)
        # Flip one random bit
        word_idx = np.random.randint(0, 16)
        bit_idx = np.random.randint(0, 32)
        M_prime = list(M)
        M_prime[word_idx] ^= (1 << bit_idx)

        hash_n = sha256_hash(M)
        hash_f = sha256_hash(M_prime)
        dhash = tuple(a ^ b for a, b in zip(hash_n, hash_f))

        if all(d == 0 for d in dhash):
            continue

        hw = sum(hamming_weight(d) for d in dhash)
        random_dhash_hw.append(hw)
        random_dhash_bits.append(state_to_bits(dhash))

    random_matrix = np.array(random_dhash_bits, dtype=np.uint8)
    random_rank = gf2_rank(random_matrix)
    print(f"  Random pairs: {len(random_dhash_bits)}")
    print(f"  Mean HW(δhash): {np.mean(random_dhash_hw):.2f}")
    print(f"  *** GF(2) rank of random δhash matrix ({random_matrix.shape[0]}×256): {random_rank} ***")

    # Part 3: Comparison
    print(f"\n  COMPARISON:")
    print(f"    Wang  rank: {wang_rank}/256")
    print(f"    Random rank: {random_rank}/256")
    if wang_rank < random_rank:
        deficit = random_rank - wang_rank
        print(f"    Rank deficit: {deficit} — Wang δhash lives in a SUBSPACE!")
        print(f"    *** THIS IS A SIGNIFICANT FINDING ***")
    elif wang_rank == random_rank:
        print(f"    No rank deficit — Wang δhash covers same dimension as random")
    else:
        print(f"    Wang rank > random rank — unexpected, investigate")

    # Rule 14: if rank < 256, triple check
    if wang_rank < 256:
        print(f"\n  Rule 14 check: wang_rank={wang_rank} < 256")
        # Re-run with more samples
        more_bits = list(wang_dhash_bits)
        for i in range(500):
            Wn, Wf, _, _ = generate_wang_pair()
            hash_n = sha256_hash(Wn)
            hash_f = sha256_hash(Wf)
            dhash = tuple(a ^ b for a, b in zip(hash_n, hash_f))
            if all(d == 0 for d in dhash):
                continue
            more_bits.append(state_to_bits(dhash))

        extended_matrix = np.array(more_bits, dtype=np.uint8)
        extended_rank = gf2_rank(extended_matrix)
        print(f"  Extended check ({len(more_bits)} Wang pairs): rank = {extended_rank}")
        if extended_rank == wang_rank:
            print(f"  Rank stable at {wang_rank} — confirmed")
        else:
            print(f"  Rank changed {wang_rank} -> {extended_rank} — was sampling artifact")
            wang_rank = extended_rank

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")
    print()
    return wang_rank, random_rank


# =============================================================
# EXPERIMENT A: Jacobian rank as function of backward rounds R
# =============================================================
def experiment_A(n_trials=1000):
    print("=" * 70)
    print("EXPERIMENT A: Jacobian rank vs number of backward rounds R")
    print("=" * 70)
    t0 = time.time()

    R_values = [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 16, 20, 24, 32, 48, 64]

    print(f"\n  {'R':>3} | {'2R':>4} | {'med_rank':>8} | {'P(r<8)':>8} | {'DOF':>6} | {'DOF/R':>6}")
    print(f"  {'---':>3}-+-{'----':>4}-+-{'--------':>8}-+-{'--------':>8}-+-{'------':>6}-+-{'------':>6}")

    results = []

    for R in R_values:
        start_round = 64 - R
        ranks = []

        for trial in range(n_trials):
            # Random state at round 64-R
            state_start = tuple(random_words(8))

            # Two random W sets of length R
            W_a = random_words(R)
            W_b = random_words(R)

            # System: F(W_a, W_b) = state64(state_start, W_a) - state64(state_start, W_b)
            def compute_state64(W_set):
                return forward_rounds(state_start, W_set, start_round)

            def compute_diff(Wa, Wb):
                sa = compute_state64(Wa)
                sb = compute_state64(Wb)
                return tuple(sub(x, y) for x, y in zip(sa, sb))

            base_diff = compute_diff(W_a, W_b)

            # Jacobian: 8 equations × 2R variables
            n_vars = 2 * R
            J = np.zeros((8, n_vars), dtype=np.float64)

            for var in range(n_vars):
                Wa_p = list(W_a)
                Wb_p = list(W_b)
                if var < R:
                    Wa_p[var] = (Wa_p[var] + 1) & MASK
                else:
                    Wb_p[var - R] = (Wb_p[var - R] + 1) & MASK

                new_diff = compute_diff(Wa_p, Wb_p)
                for eq in range(8):
                    delta = sub(new_diff[eq], base_diff[eq])
                    if delta > 0x7FFFFFFF:
                        delta -= 0x100000000
                    J[eq, var] = delta

            rank = np.linalg.matrix_rank(J, tol=1e-6)
            ranks.append(rank)

        med = np.median(ranks)
        p_lt8 = np.mean([r < 8 for r in ranks])
        dof = 2 * R - med
        dof_per_r = dof / R

        print(f"  {R:>3} | {2*R:>4} | {med:>8.0f} | {p_lt8:>8.4f} | {dof:>6.0f} | {dof_per_r:>6.2f}")
        results.append((R, 2*R, med, p_lt8, dof, dof_per_r))

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")
    print()
    return results


# =============================================================
# EXPERIMENT D: Distinguisher score vs collision distance
# =============================================================
def experiment_D(n_random=100000, n_wang=1000):
    print("=" * 70)
    print("EXPERIMENT D: Distinguisher score vs collision distance")
    print("=" * 70)
    t0 = time.time()

    # Ch invariant check: for each output word H[i],
    # check if bits 30,31 of the last Ch output satisfy Ch[b30,b31]=0
    # Simplified: we check bit patterns in H[7] (most affected by carry[63])
    def ch_invariant_score(hash_val):
        """
        Score based on T_CH_INVARIANT: check bit patterns in output words.
        For each H[i], check if bits 30 and 31 are both 0.
        This is a simplified proxy for the carry[63]=0 condition.
        """
        score = 0
        for word in hash_val:
            b30 = (word >> 30) & 1
            b31 = (word >> 31) & 1
            if b30 == 0 and b31 == 0:
                score += 1
        return score

    # Part 1: Random single-bit-flip pairs
    distances_random = []
    scores_random = []

    for i in range(n_random):
        M = random_words(16)
        word_idx = i % 16
        bit_idx = i % 32
        M_prime = list(M)
        M_prime[word_idx] ^= (1 << bit_idx)

        h1 = sha256_hash(M)
        h2 = sha256_hash(M_prime)

        dist = sum(hamming_weight(a ^ b) for a, b in zip(h1, h2))
        score = ch_invariant_score(h1) + ch_invariant_score(h2)

        distances_random.append(dist)
        scores_random.append(score)

    distances_random = np.array(distances_random)
    scores_random = np.array(scores_random)

    # Pearson correlation
    if np.std(distances_random) > 0 and np.std(scores_random) > 0:
        corr_random = np.corrcoef(distances_random, scores_random)[0, 1]
    else:
        corr_random = 0.0

    print(f"  Random pairs ({n_random}):")
    print(f"    Mean collision distance: {np.mean(distances_random):.2f}")
    print(f"    Mean Ch-invariant score: {np.mean(scores_random):.3f}")
    print(f"    Pearson r(score, distance): {corr_random:.6f}")

    # Part 2: Wang pairs
    distances_wang = []
    scores_wang = []

    for i in range(n_wang):
        Wn, Wf, _, _ = generate_wang_pair()
        h1 = sha256_hash(Wn)
        h2 = sha256_hash(Wf)

        dist = sum(hamming_weight(a ^ b) for a, b in zip(h1, h2))
        score = ch_invariant_score(h1) + ch_invariant_score(h2)

        distances_wang.append(dist)
        scores_wang.append(score)

    distances_wang = np.array(distances_wang)
    scores_wang = np.array(scores_wang)

    if np.std(distances_wang) > 0 and np.std(scores_wang) > 0:
        corr_wang = np.corrcoef(distances_wang, scores_wang)[0, 1]
    else:
        corr_wang = 0.0

    print(f"\n  Wang pairs ({n_wang}):")
    print(f"    Mean collision distance: {np.mean(distances_wang):.2f}")
    print(f"    Mean Ch-invariant score: {np.mean(scores_wang):.3f}")
    print(f"    Pearson r(score, distance): {corr_wang:.6f}")

    print(f"\n  Comparison:")
    print(f"    Random distance: {np.mean(distances_random):.2f} ± {np.std(distances_random):.2f}")
    print(f"    Wang distance:   {np.mean(distances_wang):.2f} ± {np.std(distances_wang):.2f}")
    is_closer = np.mean(distances_wang) < np.mean(distances_random)
    print(f"    Wang pairs closer to collision: {'YES' if is_closer else 'NO'}")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")
    print()
    return corr_random, corr_wang


# =============================================================
# EXPERIMENT B: Jacobian with schedule constraint
# =============================================================
def experiment_B(n_samples=1000):
    print("=" * 70)
    print("EXPERIMENT B: Jacobian with schedule constraint (real W[0..15])")
    print("=" * 70)
    t0 = time.time()

    R_values = [4, 8, 16, 32, 48]

    for R_back in R_values:
        start_round = 64 - R_back

        # For each sample: compute Jacobian d(state[64]) / d(W[0..15])
        # using single-bit perturbations (256×512 matrix → SVD rank)
        ranks = []

        for trial in range(n_samples):
            M = random_words(16)
            W = message_schedule(M)
            _, states = sha256_compress(W)
            state64_base = states[64]

            # Jacobian: 8 output words × 16 input words = 8×16 word-level
            J = np.zeros((8, 16), dtype=np.float64)

            for w_idx in range(16):
                M_p = list(M)
                M_p[w_idx] = (M_p[w_idx] + 1) & MASK
                W_p = message_schedule(M_p)
                _, states_p = sha256_compress(W_p)
                state64_p = states_p[64]

                for eq in range(8):
                    delta = sub(state64_p[eq], state64_base[eq])
                    if delta > 0x7FFFFFFF:
                        delta -= 0x100000000
                    J[eq, w_idx] = delta

            rank = np.linalg.matrix_rank(J, tol=1e-6)
            ranks.append(rank)

        med_rank = np.median(ranks)
        rank_deficit = 8 - med_rank
        rank_counts = defaultdict(int)
        for r in ranks:
            rank_counts[r] += 1

        print(f"\n  R_backward={R_back} (rounds {start_round}..63):")
        print(f"    Word-level Jacobian 8×16, rank distribution:")
        for r in sorted(rank_counts.keys()):
            print(f"      rank={r}: {rank_counts[r]} ({100*rank_counts[r]/n_samples:.1f}%)")
        print(f"    Median rank: {med_rank}")
        print(f"    Rank deficit (8-rank): {rank_deficit}")

    # Also do bit-level for R=4
    print(f"\n  Bit-level Jacobian for R_backward=4 (256×512):")
    bit_ranks = []
    for trial in range(min(100, n_samples)):
        M = random_words(16)
        W = message_schedule(M)
        _, states = sha256_compress(W)
        base_bits = state_to_bits(states[64])

        J_bits = np.zeros((256, 512), dtype=np.int8)

        for w_idx in range(16):
            for b in range(32):
                M_p = list(M)
                M_p[w_idx] ^= (1 << b)
                W_p = message_schedule(M_p)
                _, states_p = sha256_compress(W_p)
                p_bits = state_to_bits(states_p[64])

                col = w_idx * 32 + (31 - b)
                for row in range(256):
                    J_bits[row, col] = base_bits[row] ^ p_bits[row]

        rank = gf2_rank(J_bits)
        bit_ranks.append(rank)

    print(f"    GF(2) rank of 256×512 matrix ({len(bit_ranks)} samples):")
    print(f"    Mean: {np.mean(bit_ranks):.1f}, Min: {np.min(bit_ranks)}, Max: {np.max(bit_ranks)}")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")
    print()
    return


# =============================================================
# MAIN
# =============================================================
def main():
    print("=" * 70)
    print("  ЗАДАНИЕ 2: Jacobian Rank Cascade")
    print("  Priority: C > A > D > B")
    print("=" * 70)
    print()

    # Verify Wang chain first
    print("Verification: Wang chain generation...")
    ok = 0
    for _ in range(100):
        Wn, Wf, _, _ = generate_wang_pair()
        zeros, _, _, _, _ = verify_wang_pair(Wn, Wf)
        if zeros == 15:
            ok += 1
    print(f"  Wang chain: {ok}/100 have 15 zeros δe[2..16] (expected: 100)")
    assert ok == 100, f"Wang verification failed: only {ok}/100"
    print()

    # Priority C: Conditional rank
    wang_rank, random_rank = experiment_C(n_wang=1000, n_random=1000)

    # Priority A: Jacobian rank cascade
    results_A = experiment_A(n_trials=1000)

    # Priority D: Distinguisher
    corr_random, corr_wang = experiment_D(n_random=100000, n_wang=1000)

    # Priority B: Schedule constraint
    experiment_B(n_samples=1000)

    # =============================================================
    # SUMMARY
    # =============================================================
    print("=" * 70)
    print("  SUMMARY OF RESULTS")
    print("=" * 70)
    print()
    print(f"  EXPERIMENT C (KEY RESULT):")
    print(f"    Wang δhash GF(2) rank:   {wang_rank}/256")
    print(f"    Random δhash GF(2) rank: {random_rank}/256")
    print(f"    Rank deficit: {random_rank - wang_rank}")
    print()
    print(f"  EXPERIMENT A (Jacobian rank table):")
    print(f"  {'R':>3} | {'2R':>4} | {'med_rank':>8} | {'DOF':>6} | {'DOF/R':>6}")
    for R, twoR, med, p_lt8, dof, dofr in results_A:
        print(f"  {R:>3} | {twoR:>4} | {med:>8.0f} | {dof:>6.0f} | {dofr:>6.2f}")
    print()
    print(f"  EXPERIMENT D:")
    print(f"    Pearson r (random pairs): {corr_random:.6f}")
    print(f"    Pearson r (Wang pairs):   {corr_wang:.6f}")
    print()
    print(f"  Rule 7: δhash≠0 verified for all pairs")
    print(f"  Rule 12: null hypothesis (random) tested alongside Wang")
    print(f"  Rule 1: all numbers reported as computed")
    print()


if __name__ == "__main__":
    main()
