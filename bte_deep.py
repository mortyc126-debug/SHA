"""
BTE DEEP: Push the invariant finding to its limits.

1. Test at 32 and 64 rounds
2. IDENTIFY the invariant vector — what bits are in it?
3. Is it universal (same for all messages) or local?
4. If it persists at 64 rounds → this is REAL structure
"""

import random
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, schedule
)


def compute_sha_round_jacobian(state, r, W):
    """Compute 256×256 Jacobian of one SHA-256 round over GF(2)."""
    a, b, c, d, e, f, g, h = state
    T1 = (h + Sig1(e) + Ch(e, f, g) + K[r] + W[r]) & MASK
    T2 = (Sig0(a) + Maj(a, b, c)) & MASK
    base_out = [(T1 + T2) & MASK, a, b, c, (d + T1) & MASK, e, f, g]

    def to_bits(regs):
        bits = []
        for reg in regs:
            for k in range(32):
                bits.append((reg >> k) & 1)
        return bits

    base_bits = to_bits(base_out)

    J = []
    for reg_idx in range(8):
        for bit_idx in range(32):
            test_state = list(state)
            test_state[reg_idx] ^= (1 << bit_idx)
            ta, tb, tc, td, te, tf, tg, th = test_state
            tT1 = (th + Sig1(te) + Ch(te, tf, tg) + K[r] + W[r]) & MASK
            tT2 = (Sig0(ta) + Maj(ta, tb, tc)) & MASK
            test_out = [(tT1 + tT2) & MASK, ta, tb, tc, (td + tT1) & MASK, te, tf, tg]
            test_bits = to_bits(test_out)
            J.append([base_bits[i] ^ test_bits[i] for i in range(256)])

    return J


def mat_mul_gf2(A, B, n):
    """Multiply two n×n matrices over GF(2)."""
    C = [[0]*n for _ in range(n)]
    for i in range(n):
        for k in range(n):
            if A[i][k]:
                for j in range(n):
                    C[i][j] ^= B[k][j]
    return C


def find_eigenspace(J_product, n):
    """Find eigenvalue-1 subspace: null space of (J^T - I)."""
    JT = [[J_product[i][j] for i in range(n)] for j in range(n)]
    JT_minus_I = [list(JT[i]) for i in range(n)]
    for i in range(n):
        JT_minus_I[i][i] ^= 1

    # Gaussian elimination with back-substitution to find null space
    m = [list(row) for row in JT_minus_I]
    pivot_cols = []
    rank = 0
    for col in range(n):
        pivot = None
        for row in range(rank, n):
            if m[row][col] == 1:
                pivot = row
                break
        if pivot is None:
            continue
        m[rank], m[pivot] = m[pivot], m[rank]
        pivot_cols.append(col)
        for row in range(n):
            if row != rank and m[row][col] == 1:
                for c in range(n):
                    m[row][c] ^= m[rank][c]
        rank += 1

    null_dim = n - rank
    pivot_set = set(pivot_cols)
    free_cols = [c for c in range(n) if c not in pivot_set]

    # Extract null space basis vectors
    null_vectors = []
    for fc in free_cols:
        vec = [0] * n
        vec[fc] = 1
        for r_idx, pc in enumerate(pivot_cols):
            if m[r_idx][fc] == 1:
                vec[pc] = 1
        null_vectors.append(vec)

    return rank, null_dim, null_vectors


def experiment_64_round_invariant():
    """Test invariant at 32 and 64 rounds."""
    print("=" * 80)
    print("BTE DEEP: Invariant subspace at 32 and 64 rounds")
    print("=" * 80)

    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]
    states, W = sha256_round_trace(M)

    J_product = None
    checkpoints = [1, 2, 4, 8, 16, 32, 64]
    results = {}

    for r in range(8, 8 + 64):
        if r >= 72:
            break
        actual_r = r % 64  # Wrap if needed
        if actual_r >= 64:
            break

        J_r = compute_sha_round_jacobian(states[actual_r], actual_r, W)

        if J_product is None:
            J_product = J_r
        else:
            J_product = mat_mul_gf2(J_r, J_product, 256)

        n_rounds = r - 8 + 1
        if n_rounds in checkpoints:
            _, dim, vectors = find_eigenspace(J_product, 256)
            results[n_rounds] = (dim, vectors)
            print(f"  {n_rounds:>2} rounds: eigenvalue-1 dim = {dim}")

    return results


def experiment_identify_invariant():
    """
    WHAT IS the invariant vector?
    For each round count, extract the invariant vector(s) and analyze:
    - Which registers are involved?
    - Which bit positions?
    - Is there a pattern?
    """
    print("\n" + "=" * 80)
    print("BTE DEEP: Identify the invariant vectors")
    print("=" * 80)

    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]
    states, W = sha256_round_trace(M)

    reg_names = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']

    J_product = None
    for r in range(8, 8 + 16):
        J_r = compute_sha_round_jacobian(states[r], r, W)
        if J_product is None:
            J_product = J_r
        else:
            J_product = mat_mul_gf2(J_r, J_product, 256)

    _, dim, vectors = find_eigenspace(J_product, 256)
    print(f"  16 rounds: dim = {dim}")

    for vi, vec in enumerate(vectors[:4]):
        print(f"\n  Invariant vector {vi}:")
        hw = sum(vec)
        print(f"    Total HW: {hw}/256")

        for reg in range(8):
            reg_bits = vec[reg*32:(reg+1)*32]
            reg_hw = sum(reg_bits)
            if reg_hw > 0:
                positions = [k for k in range(32) if reg_bits[k]]
                print(f"    {reg_names[reg]}: HW={reg_hw}, positions={positions}")


def experiment_universality():
    """
    Is the invariant the SAME for all messages?
    If yes → it's a structural property of SHA-256 itself.
    If no → it's local (state-dependent).
    """
    print("\n" + "=" * 80)
    print("BTE DEEP: Universality test — same invariant for all messages?")
    print("=" * 80)

    # For 10 different messages, compute the 16-round invariant
    # Check if they're the same vector or different
    all_vectors = []

    for seed in range(10):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        states, W = sha256_round_trace(M)

        J_product = None
        for r in range(8, 24):
            J_r = compute_sha_round_jacobian(states[r], r, W)
            if J_product is None:
                J_product = J_r
            else:
                J_product = mat_mul_gf2(J_r, J_product, 256)

        _, dim, vectors = find_eigenspace(J_product, 256)
        all_vectors.append((dim, vectors))
        print(f"  Seed {seed}: dim = {dim}, HW of first vector = {sum(vectors[0]) if vectors else 'N/A'}")

    # Compare vectors pairwise
    print(f"\n  Pairwise comparison of invariant vectors:")
    for i in range(min(5, len(all_vectors))):
        for j in range(i+1, min(5, len(all_vectors))):
            if all_vectors[i][1] and all_vectors[j][1]:
                v1 = all_vectors[i][1][0]
                v2 = all_vectors[j][1][0]
                hamming = sum(a ^ b for a, b in zip(v1, v2))
                print(f"    Seeds {i} vs {j}: Hamming distance = {hamming}/256")

    # If Hamming distances are ~128 → vectors are unrelated (local)
    # If Hamming distances are ~0 → same vector (universal)
    # If Hamming distances cluster → partial universality


def experiment_invariant_functional_test():
    """
    THE REAL TEST: Does the invariant actually WORK?

    If v is an eigenvalue-1 vector of J_product (16 rounds):
      v · state[r] should equal v · state[r+16] (mod 2)
    for the SPECIFIC message used to compute J.

    And for OTHER messages: does it still hold approximately?
    """
    print("\n" + "=" * 80)
    print("BTE DEEP: Functional test — does the invariant predict anything?")
    print("=" * 80)

    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]
    states, W = sha256_round_trace(M)

    # Build 16-round Jacobian product starting from round 8
    J_product = None
    for r in range(8, 24):
        J_r = compute_sha_round_jacobian(states[r], r, W)
        if J_product is None:
            J_product = J_r
        else:
            J_product = mat_mul_gf2(J_r, J_product, 256)

    _, dim, vectors = find_eigenspace(J_product, 256)

    if not vectors:
        print("  No invariant vectors found.")
        return

    v = vectors[0]

    def dot_state(vec, state):
        """Compute v · state mod 2."""
        result = 0
        for reg in range(8):
            for k in range(32):
                if vec[reg * 32 + k]:
                    result ^= (state[reg] >> k) & 1
        return result

    # Test on the SAME message
    print(f"  Same message (seed=42):")
    print(f"    v · state[8]  = {dot_state(v, states[8])}")
    print(f"    v · state[24] = {dot_state(v, states[24])}")
    print(f"    v · state[40] = {dot_state(v, states[40])}")
    print(f"    v · state[56] = {dot_state(v, states[56])}")

    match_same = dot_state(v, states[8]) == dot_state(v, states[24])
    print(f"    state[8] vs state[24]: {'MATCH' if match_same else 'DIFFER'}")

    # Test on DIFFERENT messages
    print(f"\n  Different messages (v · state[8] vs v · state[24]):")
    match_count = 0
    N = 500
    for seed in range(N):
        rng_s = random.Random(seed)
        M_s = [rng_s.randint(0, MASK) for _ in range(16)]
        states_s, _ = sha256_round_trace(M_s)
        d8 = dot_state(v, states_s[8])
        d24 = dot_state(v, states_s[24])
        if d8 == d24:
            match_count += 1

    print(f"    Matches: {match_count}/{N} = {match_count/N:.4f} (random = 0.500)")

    # Test at different starting rounds
    print(f"\n  v · state[r] vs v · state[r+16] for various r:")
    for start_r in [0, 4, 8, 16, 32, 48]:
        end_r = start_r + 16
        if end_r > 64:
            break
        match_count = 0
        for seed in range(N):
            rng_s = random.Random(seed)
            M_s = [rng_s.randint(0, MASK) for _ in range(16)]
            states_s, _ = sha256_round_trace(M_s)
            if dot_state(v, states_s[start_r]) == dot_state(v, states_s[end_r]):
                match_count += 1
        print(f"    r={start_r:>2}→{end_r:>2}: {match_count}/{N} = {match_count/N:.4f}")


if __name__ == "__main__":
    results = experiment_64_round_invariant()
    experiment_identify_invariant()
    experiment_universality()
    experiment_invariant_functional_test()
