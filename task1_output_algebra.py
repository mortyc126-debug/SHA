#!/usr/bin/env python3
"""
ЗАДАНИЕ 1: Output Algebra — обратная диффузия раундов 60-64
SHA-256 reverse diffusion analysis.

Implements:
1. Full SHA-256 (64 rounds) with state capture
2. T_INVERSE backward pass (state[r+1] -> state[r])
3. Correlation matrix state[60] vs W[0..15] (256x512)
4. Collision funnel: Jacobian rank analysis
5. Backward diffusion speed: HW(Δstate[r])

All results averaged over 10000 samples.
"""

import struct
import os
import numpy as np
from collections import defaultdict

MASK = 0xFFFFFFFF

# SHA-256 constants
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
    """Expand 16-word message to 64-word schedule."""
    W = list(M_words[:16])
    for i in range(16, 64):
        W.append(add(sig1(W[i-2]), W[i-7], sig0(W[i-15]), W[i-16]))
    return W


def sha256_rounds(W, initial_state=None):
    """
    Run SHA-256 round function, return all intermediate states.
    state[0] = initial state (H0 or custom)
    state[r] = state after round r-1 (so state[64] is final)
    """
    if initial_state is None:
        state = list(H0)
    else:
        state = list(initial_state)

    states = [tuple(state)]
    a, b, c, d, e, f, g, h = state

    for r in range(64):
        T1 = add(h, Sig1(e), Ch(e, f, g), K[r], W[r])
        T2 = add(Sig0(a), Maj(a, b, c))
        h = g
        g = f
        f = e
        e = add(d, T1)
        d = c
        c = b
        b = a
        a = add(T1, T2)
        states.append((a, b, c, d, e, f, g, h))

    return states


def inverse_round(state_next, W_r, r):
    """
    T_INVERSE: given state[r+1] and W[r], recover state[r].

    From methodology:
      a=b'; b=c'; c=d'; e=f'; f=g'; g=h'
      T2 = Sig0(b') + Maj(b',c',d')
      T1 = a' - T2  (mod 2^32)
      d  = e' - T1  (mod 2^32)
      h  = T1 - Sig1(f') - Ch(f',g',h') - K[r] - W[r]  (mod 2^32)
    """
    a_next, b_next, c_next, d_next, e_next, f_next, g_next, h_next = state_next

    # Recover shifted registers
    a = b_next
    b = c_next
    c = d_next
    e = f_next
    f = g_next
    g = h_next

    # Compute T2 and T1
    T2 = add(Sig0(b_next), Maj(b_next, c_next, d_next))
    T1 = sub(a_next, T2)

    # Recover d and h
    d = sub(e_next, T1)
    h = sub(sub(sub(sub(T1, Sig1(f_next)), Ch(f_next, g_next, h_next)), K[r]), W_r)

    return (a, b, c, d, e, f, g, h)


def backward_pass(state_end, W, from_round, to_round):
    """Run backward pass from state[from_round] to state[to_round]."""
    state = state_end
    for r in range(from_round - 1, to_round - 1, -1):
        state = inverse_round(state, W[r], r)
    return state


def random_message():
    """Generate random 512-bit message as 16 uint32 words."""
    data = os.urandom(64)
    return list(struct.unpack('>16I', data))


def state_to_bits(state):
    """Convert 8-word state to 256-bit array."""
    bits = []
    for word in state:
        for bit in range(31, -1, -1):
            bits.append((word >> bit) & 1)
    return bits


def words_to_bits(words):
    """Convert list of uint32 words to bit array."""
    bits = []
    for word in words:
        for bit in range(31, -1, -1):
            bits.append((word >> bit) & 1)
    return bits


def hamming_weight_xor(s1, s2):
    """Count differing bits between two states."""
    hw = 0
    for a, b in zip(s1, s2):
        hw += bin(a ^ b).count('1')
    return hw


def gf2_rank(matrix):
    """Compute rank of binary matrix over GF(2)."""
    M = np.array(matrix, dtype=np.int8) % 2
    rows, cols = M.shape
    rank = 0
    used_rows = set()

    for col in range(cols):
        pivot = -1
        for row in range(rows):
            if row not in used_rows and M[row, col] == 1:
                pivot = row
                break
        if pivot == -1:
            continue
        rank += 1
        used_rows.add(pivot)
        for row in range(rows):
            if row != pivot and M[row, col] == 1:
                M[row] = (M[row] ^ M[pivot]) % 2

    return rank


# =============================================================
# EXPERIMENT 1: Round-trip verification (Rule 5)
# =============================================================
def verify_roundtrip(n=1000):
    print("=" * 70)
    print("VERIFICATION: Round-trip identity (forward + backward)")
    print("=" * 70)
    ok = 0
    for _ in range(n):
        M = random_message()
        W = message_schedule(M)
        states = sha256_rounds(W)

        # Backward 64 -> 0
        recovered = states[64]
        for r in range(63, -1, -1):
            recovered = inverse_round(recovered, W[r], r)

        if recovered == states[0]:
            ok += 1

    print(f"  Round-trip test: {ok}/{n} passed")
    assert ok == n, f"FAILED: only {ok}/{n} round-trips matched!"
    print("  STATUS: PASS (identity verified)")
    print()
    return True


# =============================================================
# EXPERIMENT 2: Backward pass state[64] -> state[60] -> state[56]
# =============================================================
def experiment_backward_states(n=10000):
    print("=" * 70)
    print("EXPERIMENT 2: Backward pass state[64]->state[60]->state[56]")
    print("=" * 70)

    # Collect data for correlation analysis
    state60_bits_all = []
    state56_bits_all = []
    w0_15_bits_all = []
    w56_63_all = []
    w60_63_all = []

    for i in range(n):
        M = random_message()
        W = message_schedule(M)
        states = sha256_rounds(W)

        # Backward 4 rounds: state[64] -> state[60]
        state60 = backward_pass(states[64], W, 64, 60)
        # Verify against forward
        assert state60 == states[60], f"Backward 64->60 mismatch at sample {i}"

        # Backward 8 rounds: state[64] -> state[56]
        state56 = backward_pass(states[64], W, 64, 56)
        assert state56 == states[56], f"Backward 64->56 mismatch at sample {i}"

        state60_bits_all.append(state_to_bits(state60))
        state56_bits_all.append(state_to_bits(state56))
        w0_15_bits_all.append(words_to_bits(M))
        w56_63_all.append(W[56:64])
        w60_63_all.append(W[60:64])

    print(f"  Backward pass verified for all {n} samples")
    print(f"  state[60] and state[56] match forward computation: OK")
    print()

    return (np.array(state60_bits_all), np.array(state56_bits_all),
            np.array(w0_15_bits_all), w56_63_all, w60_63_all)


# =============================================================
# EXPERIMENT 3: Correlation matrix 256x512
# =============================================================
def experiment_correlation_matrix(state60_bits, w0_15_bits):
    print("=" * 70)
    print("EXPERIMENT 3: Correlation matrix state[60] vs W[0..15]")
    print("=" * 70)

    n = state60_bits.shape[0]
    n_state_bits = 256
    n_input_bits = 512

    # Compute correlation: convert 0/1 to -1/+1 for proper correlation
    S = 2.0 * state60_bits.astype(np.float64) - 1.0  # (n, 256)
    W = 2.0 * w0_15_bits.astype(np.float64) - 1.0    # (n, 512)

    # Correlation matrix: corr[i,j] = mean(S[:,i] * W[:,j])
    corr = (S.T @ W) / n  # (256, 512)

    # Statistics
    max_corr = np.max(np.abs(corr))
    mean_corr = np.mean(np.abs(corr))
    print(f"  Max |correlation|:  {max_corr:.6f}")
    print(f"  Mean |correlation|: {mean_corr:.6f}")

    # Top-20 largest correlations
    flat = np.abs(corr).flatten()
    top_idx = np.argsort(flat)[-20:][::-1]
    print(f"\n  Top-20 largest |correlations|:")
    print(f"  {'Rank':>4}  {'state[60] bit':>13}  {'W[0..15] bit':>13}  {'corr':>10}")
    for rank, idx in enumerate(top_idx):
        i = idx // 512
        j = idx % 512
        word_i = i // 32
        bit_i = 31 - (i % 32)
        word_j = j // 32
        bit_j = 31 - (j % 32)
        print(f"  {rank+1:>4}  reg{word_i}[{bit_i:>2}] ({i:>3})  W{word_j:>2}[{bit_j:>2}] ({j:>3})  {corr[i,j]:>+10.6f}")

    # GF(2) rank of binarized correlation matrix
    binary_corr = (np.abs(corr) > 0.1).astype(np.int8)
    n_nonzero = np.sum(binary_corr)
    print(f"\n  Binarized matrix (|corr|>0.1): {n_nonzero} nonzero entries")

    if n_nonzero > 0:
        rank = gf2_rank(binary_corr)
        print(f"  GF(2) rank of binarized matrix: {rank}")
    else:
        rank = 0
        print(f"  GF(2) rank: 0 (no entries above threshold)")

    # Weakly dependent bits (corr < 0.01 with ALL input bits)
    max_corr_per_state_bit = np.max(np.abs(corr), axis=1)
    weak_bits = np.where(max_corr_per_state_bit < 0.01)[0]
    moderate_bits = np.where(max_corr_per_state_bit < 0.05)[0]
    print(f"\n  Bits with max|corr|<0.01 (weak dep): {len(weak_bits)}/256")
    print(f"  Bits with max|corr|<0.05 (moderate):  {len(moderate_bits)}/256")

    if len(weak_bits) > 0 and len(weak_bits) <= 30:
        for b in weak_bits:
            word = b // 32
            bit = 31 - (b % 32)
            print(f"    state[60] reg{word}[{bit}]: max|corr| = {max_corr_per_state_bit[b]:.6f}")

    print()
    return corr, rank


# =============================================================
# EXPERIMENT 4: Collision Funnel — Jacobian Rank
# =============================================================
def experiment_collision_funnel(n_trials=10000):
    print("=" * 70)
    print("EXPERIMENT 4: Collision Funnel — free-start, Jacobian rank")
    print("=" * 70)

    jacobian_ranks = []
    dof_list = []

    for trial in range(n_trials):
        # Generate random state[60]
        state60 = tuple(struct.unpack('>8I', os.urandom(32)))

        # Generate two different sets of W[60..63]
        W_a = list(struct.unpack('>4I', os.urandom(16)))
        W_b = list(struct.unpack('>4I', os.urandom(16)))

        # Forward 4 rounds from state[60] with W_a
        state_a = state60
        for r in range(60, 64):
            a, b, c, d, e, f, g, h = state_a
            T1 = add(h, Sig1(e), Ch(e, f, g), K[r], W_a[r - 60])
            T2 = add(Sig0(a), Maj(a, b, c))
            state_a = (add(T1, T2), a, b, c, add(d, T1), e, f, g)
        state64_a = state_a

        # Forward 4 rounds from state[60] with W_b
        state_b = state60
        for r in range(60, 64):
            a, b, c, d, e, f, g, h = state_b
            T1 = add(h, Sig1(e), Ch(e, f, g), K[r], W_b[r - 60])
            T2 = add(Sig0(a), Maj(a, b, c))
            state_b = (add(T1, T2), a, b, c, add(d, T1), e, f, g)
        state64_b = state_b

        # Compute Jacobian numerically: d(state64_a - state64_b) / d(W_a, W_b)
        # System: F(W_a[0..3], W_b[0..3]) = state64_a - state64_b = 0
        # 8 equations (state diff), 8 variables (W_a[0..3] and W_b[0..3])
        # Jacobian via finite differences (bit-level)

        eps = 1  # single-bit perturbation in LSB
        J = np.zeros((8, 8), dtype=np.float64)

        def compute_diff(Wa, Wb):
            """Compute state64(state60, Wa) - state64(state60, Wb) mod 2^32."""
            sa = state60
            for r in range(60, 64):
                a, b, c, d, e, f, g, h = sa
                T1 = add(h, Sig1(e), Ch(e, f, g), K[r], Wa[r - 60])
                T2 = add(Sig0(a), Maj(a, b, c))
                sa = (add(T1, T2), a, b, c, add(d, T1), e, f, g)
            sb = state60
            for r in range(60, 64):
                a, b, c, d, e, f, g, h = sb
                T1 = add(h, Sig1(e), Ch(e, f, g), K[r], Wb[r - 60])
                T2 = add(Sig0(a), Maj(a, b, c))
                sb = (add(T1, T2), a, b, c, add(d, T1), e, f, g)
            return tuple(sub(x, y) for x, y in zip(sa, sb))

        base_diff = compute_diff(W_a, W_b)

        for var in range(8):
            Wa_p = list(W_a)
            Wb_p = list(W_b)
            if var < 4:
                Wa_p[var] = (Wa_p[var] + eps) & MASK
            else:
                Wb_p[var - 4] = (Wb_p[var - 4] + eps) & MASK

            new_diff = compute_diff(Wa_p, Wb_p)
            for eq in range(8):
                # Numerical derivative (mod 2^32 arithmetic)
                delta = sub(new_diff[eq], base_diff[eq])
                # Convert to signed
                if delta > 0x7FFFFFFF:
                    delta -= 0x100000000
                J[eq, var] = delta

        # Rank of Jacobian (over reals)
        rank = np.linalg.matrix_rank(J)
        jacobian_ranks.append(rank)

    rank_counts = defaultdict(int)
    for r in jacobian_ranks:
        rank_counts[r] += 1

    mean_rank = np.mean(jacobian_ranks)
    print(f"  Jacobian rank statistics ({n_trials} trials):")
    for r in sorted(rank_counts.keys()):
        pct = 100.0 * rank_counts[r] / n_trials
        print(f"    rank={r}: {rank_counts[r]:>6} ({pct:>6.2f}%)")
    print(f"  Mean rank: {mean_rank:.3f}")
    print(f"  Degrees of freedom (8 - rank): {8 - mean_rank:.3f}")

    # Analysis
    if mean_rank >= 7.9:
        print("  CONCLUSION: System is (nearly) fully determined.")
        print("    Same state[60] + different W[60..63] -> almost always different state[64]")
        print("    Collision funnel has ~0 degrees of freedom")
    elif mean_rank < 7:
        print(f"  CONCLUSION: System has ~{8 - mean_rank:.1f} degrees of freedom")
        print("    Collision funnel exists with nontrivial dimension")

    print()
    return mean_rank, jacobian_ranks


# =============================================================
# EXPERIMENT 5: Backward Diffusion Speed
# =============================================================
def experiment_backward_diffusion(n=10000):
    print("=" * 70)
    print("EXPERIMENT 5: Backward diffusion — HW(Δstate[r])")
    print("=" * 70)

    # For each round r from 63 down to 60:
    # Fix state[64], vary W[r], measure HW(state[r] XOR state[r]')
    results = {}  # r -> list of HW values

    for r_vary in range(63, 59, -1):
        hw_list = []
        for _ in range(n):
            # Random state[64]
            state64 = tuple(struct.unpack('>8I', os.urandom(32)))

            # Random W for rounds 60..63
            W_base = list(struct.unpack('>4I', os.urandom(16)))
            W_alt = list(W_base)
            # Vary only W[r_vary]: replace with new random value
            new_word = struct.unpack('>I', os.urandom(4))[0]
            W_alt[r_vary - 60] = new_word

            if W_base[r_vary - 60] == new_word:
                # Skip trivial case (Rule 7: ΔW≠0)
                hw_list.append(None)
                continue

            # Backward from state[64] to state[r_vary] with base W
            state_base = state64
            for r in range(63, r_vary - 1, -1):
                state_base = inverse_round(state_base, W_base[r - 60], r)

            # Backward from state[64] to state[r_vary] with alt W
            state_alt = state64
            for r in range(63, r_vary - 1, -1):
                if r == r_vary:
                    state_alt = inverse_round(state_alt, W_alt[r - 60], r)
                else:
                    state_alt = inverse_round(state_alt, W_base[r - 60], r)

            hw = hamming_weight_xor(state_base, state_alt)
            hw_list.append(hw)

        # Filter out None (trivial cases)
        hw_valid = [x for x in hw_list if x is not None]
        results[r_vary] = hw_valid

    print(f"\n  {'Round r':>10}  {'Mean HW(Δstate)':>16}  {'Std':>8}  {'Min':>6}  {'Max':>6}  {'N':>6}")
    print(f"  {'-'*10}  {'-'*16}  {'-'*8}  {'-'*6}  {'-'*6}  {'-'*6}")

    for r in range(63, 59, -1):
        hw = results[r]
        if len(hw) > 0:
            mean_hw = np.mean(hw)
            std_hw = np.std(hw)
            min_hw = np.min(hw)
            max_hw = np.max(hw)
            print(f"  r={r:>7}  {mean_hw:>16.2f}  {std_hw:>8.2f}  {min_hw:>6}  {max_hw:>6}  {len(hw):>6}")

    # Additional: measure cumulative diffusion
    # Vary W[63] only, measure state at r=63,62,61,60
    print(f"\n  Cumulative backward diffusion (varying W[63] only):")
    print(f"  {'Measured at':>12}  {'Mean HW(Δstate)':>16}  {'Std':>8}")
    print(f"  {'-'*12}  {'-'*16}  {'-'*8}")

    for target_r in range(63, 59, -1):
        hw_list = []
        for _ in range(n):
            state64 = tuple(struct.unpack('>8I', os.urandom(32)))
            W_base = list(struct.unpack('>4I', os.urandom(16)))
            W_alt = list(W_base)
            W_alt[3] = struct.unpack('>I', os.urandom(4))[0]  # W[63]

            if W_base[3] == W_alt[3]:
                continue

            # Backward to target_r with base
            state_base = state64
            for r in range(63, target_r - 1, -1):
                state_base = inverse_round(state_base, W_base[r - 60], r)

            # Backward to target_r with alt W[63]
            state_alt = state64
            for r in range(63, target_r - 1, -1):
                w = W_alt[r - 60] if r == 63 else W_base[r - 60]
                state_alt = inverse_round(state_alt, w, r)

            hw_list.append(hamming_weight_xor(state_base, state_alt))

        if hw_list:
            print(f"  state[{target_r}]     {np.mean(hw_list):>16.2f}  {np.std(hw_list):>8.2f}")

    print()
    return results


# =============================================================
# MAIN
# =============================================================
def main():
    print("=" * 70)
    print("  ЗАДАНИЕ 1: Output Algebra — обратная диффузия раундов 60-64")
    print("  SHA-256 Reverse Diffusion Analysis")
    print("  N = 10000 samples")
    print("=" * 70)
    print()

    # Step 0: Verification (Rule 5 — double verification)
    verify_roundtrip(n=1000)

    # Step 2: Backward pass and data collection
    state60_bits, state56_bits, w0_15_bits, w56_63, w60_63 = \
        experiment_backward_states(n=10000)

    # Step 3: Correlation matrix
    corr, gf2_r = experiment_correlation_matrix(state60_bits, w0_15_bits)

    # Step 4: Collision funnel
    mean_rank, jac_ranks = experiment_collision_funnel(n_trials=10000)

    # Step 5: Backward diffusion speed
    diff_results = experiment_backward_diffusion(n=10000)

    # =============================================================
    # SUMMARY
    # =============================================================
    print("=" * 70)
    print("  SUMMARY OF RESULTS")
    print("=" * 70)
    print()
    print("  1. Round-trip verification: PASS (1000/1000)")
    print(f"  2. Correlation matrix 256×512:")
    print(f"     - Max |corr|: {np.max(np.abs(corr)):.6f}")
    print(f"     - Mean |corr|: {np.mean(np.abs(corr)):.6f}")
    print(f"     - GF(2) rank (threshold 0.1): {gf2_r}")
    weak = np.sum(np.max(np.abs(corr), axis=1) < 0.01)
    print(f"     - Weakly dependent bits (<0.01): {weak}/256")
    print(f"  3. Collision funnel Jacobian rank: {np.mean(jac_ranks):.3f} (mean)")
    print(f"     Degrees of freedom: {8 - np.mean(jac_ranks):.3f}")
    print(f"  4. Backward diffusion (HW table):")
    for r in range(63, 59, -1):
        hw = diff_results[r]
        if hw:
            print(f"     r={r}: mean HW = {np.mean(hw):.2f}")
    print()

    # Rule 7 check: verify results are nontrivial
    print("  Rule 7 (artifact check):")
    print(f"    - Correlations near 0 (expected for 60 rounds of mixing): "
          f"{'YES' if np.max(np.abs(corr)) < 0.5 else 'NO - investigate!'}")
    print(f"    - Jacobian rank ~8 (fully determined): "
          f"{'YES' if np.mean(jac_ranks) > 7.5 else 'NO'}")
    print(f"    - HW(Δstate) > 0 (nontrivial diffusion): "
          f"{'YES' if all(np.mean(diff_results[r]) > 0 for r in range(60, 64) if diff_results[r]) else 'NO'}")
    print()
    print("  Rule 1 (honest result): All numbers reported as computed.")
    print()


if __name__ == "__main__":
    main()
