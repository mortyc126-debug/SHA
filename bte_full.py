"""
BTE FULL: Extend Bi-Temporal Element to model complete SHA-256.

Extensions needed:
1. V = MAJ(T1_bit, T2_bit, carry) instead of AND(x, carry)
2. Two operand streams (T1 and T2) feeding the vertical coupling
3. Multiple additions per round (6 in SHA-256)
4. 8 registers with shift structure

Key question: Does the invariant (dim=1) survive these extensions?
"""

import random
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, schedule
)


def build_full_bte_matrix(n=32):
    """
    Build the COMPLETE BTE evolution matrix for one SHA-256 round.

    One round takes 256 input bits (8 registers × 32 bits)
    and produces 256 output bits.

    The evolution matrix E is 256×256 over GF(2).
    E[i][j] = 1 if output bit i depends on input bit j
    (in the LINEAR part — ignoring carry nonlinearity).

    For the FULL nonlinear case, we compute E by flipping each input
    bit and observing which output bits change.

    Note: this gives the JACOBIAN at a specific point, not a global matrix.
    For a linear function, the Jacobian is constant. For nonlinear SHA-256,
    it varies with the state.
    """
    pass  # Will use numerical approach instead


def compute_sha_round_jacobian(state, r, W):
    """
    Compute 256×256 Jacobian of one SHA-256 round over GF(2).

    J[i][j] = 1 if flipping input bit j changes output bit i.
    Input/output: 256 bits = 8 registers × 32 bits.
    """
    base_state = list(state)

    # Compute base output
    a, b, c, d, e, f, g, h = base_state
    T1 = (h + Sig1(e) + Ch(e, f, g) + K[r] + W[r]) & MASK
    T2 = (Sig0(a) + Maj(a, b, c)) & MASK
    base_out = [(T1 + T2) & MASK, a, b, c, (d + T1) & MASK, e, f, g]

    # Flatten to bits
    def to_bits(regs):
        bits = []
        for reg in regs:
            for k in range(32):
                bits.append((reg >> k) & 1)
        return bits

    base_bits = to_bits(base_out)

    # Build Jacobian: flip each of 256 input bits
    J = []
    for reg_idx in range(8):
        for bit_idx in range(32):
            test_state = list(base_state)
            test_state[reg_idx] ^= (1 << bit_idx)

            ta, tb, tc, td, te, tf, tg, th = test_state
            tT1 = (th + Sig1(te) + Ch(te, tf, tg) + K[r] + W[r]) & MASK
            tT2 = (Sig0(ta) + Maj(ta, tb, tc)) & MASK
            test_out = [(tT1 + tT2) & MASK, ta, tb, tc, (td + tT1) & MASK, te, tf, tg]

            test_bits = to_bits(test_out)
            diff = [base_bits[i] ^ test_bits[i] for i in range(256)]
            J.append(diff)

    return J  # 256 rows (input bits) × 256 columns (output bits)


def experiment_round_jacobian():
    """
    Compute the Jacobian of one SHA-256 round and analyze its structure.

    Key questions:
    1. What is the rank? (Should be 256 — bijection)
    2. Does it have invariant subspaces?
    3. How does it vary with the state?
    """
    print("=" * 80)
    print("FULL BTE: SHA-256 round Jacobian analysis")
    print("=" * 80)

    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]
    states, W = sha256_round_trace(M)

    for r in [0, 8, 32, 63]:
        J = compute_sha_round_jacobian(states[r], r, W)

        # Rank
        rank = gf2_rank_square(J, 256)

        # Row weights (how many output bits does each input bit affect?)
        row_weights = [sum(row) for row in J]
        avg_rw = sum(row_weights) / 256

        # Column weights (how many input bits affect each output bit?)
        col_weights = [sum(J[i][j] for i in range(256)) for j in range(256)]
        avg_cw = sum(col_weights) / 256

        print(f"\n  Round {r}: Jacobian 256×256")
        print(f"    Rank: {rank}/256")
        print(f"    Avg row weight: {avg_rw:.1f} (each input affects ~{avg_rw:.0f} outputs)")
        print(f"    Avg col weight: {avg_cw:.1f} (each output depends on ~{avg_cw:.0f} inputs)")

        # Register-level structure: how much does each input register
        # affect each output register?
        print(f"    Register coupling matrix (HW of 32×32 blocks):")
        reg_names = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
        header = "      " + " ".join(f"{name:>4}" for name in reg_names)
        print(header)
        for out_reg in range(8):
            couplings = []
            for in_reg in range(8):
                total = 0
                for out_bit in range(32):
                    for in_bit in range(32):
                        total += J[in_reg * 32 + in_bit][out_reg * 32 + out_bit]
                couplings.append(total)
            line = f"    {reg_names[out_reg]} " + " ".join(f"{c:>4}" for c in couplings)
            print(line)


def experiment_invariant_subspace():
    """
    Search for INVARIANT SUBSPACES of the round Jacobian.

    An invariant subspace S: J(S) ⊆ S.
    If such S exists → evolution within S is a separate (simpler) system.

    The register shifts create obvious invariant structure:
    b_new = a, c_new = b, etc. → the (b,c,d,f,g,h) subspace is just shifted copies.
    Only a_new and e_new are "computed."

    But are there OTHER invariant subspaces?
    Specifically: are there linear combinations of bits that are preserved?
    """
    print("\n" + "=" * 80)
    print("FULL BTE: Invariant subspace search")
    print("=" * 80)

    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]
    states, W = sha256_round_trace(M)

    # Compute Jacobian at a specific round
    r = 20
    J = compute_sha_round_jacobian(states[r], r, W)

    # J is 256×256 over GF(2). We want to find the invariant subspace:
    # vectors v such that v · J = λ · v for some λ ∈ {0, 1}.
    # Over GF(2), eigenvalues are 0 or 1.
    # Eigenvectors with eigenvalue 1: J^T · v = v → (J^T - I) · v = 0 → null space of (J^T - I).

    # Build J^T
    JT = [[J[i][j] for i in range(256)] for j in range(256)]

    # Compute (J^T - I) = J^T XOR I
    JT_minus_I = []
    for i in range(256):
        row = list(JT[i])
        row[i] ^= 1
        JT_minus_I.append(row)

    # Rank of (J^T - I)
    rank = gf2_rank_square(JT_minus_I, 256)
    null_dim = 256 - rank

    print(f"  Round {r}: rank(J^T - I) = {rank}/256")
    print(f"  Eigenvalue-1 subspace dimension: {null_dim}")

    if null_dim > 0:
        print(f"  *** {null_dim} INVARIANT DIRECTIONS found! ***")

    # Eigenvalue 0: J^T · v = 0 → null space of J^T
    rank_JT = gf2_rank_square(JT, 256)
    null_JT = 256 - rank_JT

    print(f"  rank(J^T) = {rank_JT}/256")
    print(f"  Eigenvalue-0 subspace dimension: {null_JT}")

    # Check if this varies with state
    print(f"\n  Variation across rounds:")
    for r in [0, 8, 16, 32, 63]:
        J = compute_sha_round_jacobian(states[r], r, W)
        JT = [[J[i][j] for i in range(256)] for j in range(256)]
        JT_minus_I = [list(JT[i]) for i in range(256)]
        for i in range(256):
            JT_minus_I[i][i] ^= 1

        rank = gf2_rank_square(JT_minus_I, 256)
        null_dim = 256 - rank
        print(f"    r={r:>2}: eigenvalue-1 dim = {null_dim}")

    # Check multiple messages
    print(f"\n  Variation across messages (round 20):")
    for seed in range(5):
        rng_s = random.Random(seed)
        M_s = [rng_s.randint(0, MASK) for _ in range(16)]
        states_s, W_s = sha256_round_trace(M_s)
        J = compute_sha_round_jacobian(states_s[20], 20, W_s)
        JT = [[J[i][j] for i in range(256)] for j in range(256)]
        JT_minus_I = [list(JT[i]) for i in range(256)]
        for i in range(256):
            JT_minus_I[i][i] ^= 1
        rank = gf2_rank_square(JT_minus_I, 256)
        null_dim = 256 - rank
        print(f"    Seed {seed}: eigenvalue-1 dim = {null_dim}")


def experiment_multi_round_invariant():
    """
    The REAL test: does the eigenvalue-1 subspace persist across MULTIPLE rounds?

    Compute J_1 (round r), J_2 (round r+1). Product J_2 · J_1 = Jacobian of 2 rounds.
    Find eigenvalue-1 subspace of J_2 · J_1. Compare with individual rounds.

    If dim stays constant → persistent invariant → THIS IS THE AXIS.
    If dim drops to 0 → no multi-round invariant → need different approach.
    """
    print("\n" + "=" * 80)
    print("FULL BTE: Multi-round invariant subspace")
    print("=" * 80)

    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]
    states, W = sha256_round_trace(M)

    def mat_mul_gf2(A, B, n):
        """Multiply two n×n matrices over GF(2)."""
        C = [[0]*n for _ in range(n)]
        for i in range(n):
            for k in range(n):
                if A[i][k]:
                    for j in range(n):
                        C[i][j] ^= B[k][j]
        return C

    # Compute product of Jacobians for 1, 2, 4, 8, 16 rounds
    J_product = None
    invariant_dims = []

    for r in range(8, 24):
        J_r = compute_sha_round_jacobian(states[r], r, W)

        if J_product is None:
            J_product = J_r
        else:
            J_product = mat_mul_gf2(J_r, J_product, 256)

        n_rounds = r - 8 + 1
        if n_rounds in [1, 2, 4, 8, 16]:
            # Find eigenvalue-1 subspace
            JT = [[J_product[i][j] for i in range(256)] for j in range(256)]
            JT_minus_I = [list(JT[i]) for i in range(256)]
            for i in range(256):
                JT_minus_I[i][i] ^= 1

            rank = gf2_rank_square(JT_minus_I, 256)
            dim = 256 - rank
            invariant_dims.append((n_rounds, dim))
            print(f"  {n_rounds:>2} rounds: eigenvalue-1 dim = {dim}")

    print(f"\n  If dim persists > 0 → INVARIANT SUBSPACE across rounds!")
    print(f"  If dim → 0 → no persistent invariant")


def gf2_rank_square(matrix, n):
    """Compute rank of n×n binary matrix over GF(2)."""
    m = [list(row[:n]) for row in matrix[:n]]
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
        for row in range(n):
            if row != rank and m[row][col] == 1:
                for c in range(n):
                    m[row][c] ^= m[rank][c]
        rank += 1
    return rank


if __name__ == "__main__":
    experiment_round_jacobian()
    experiment_invariant_subspace()
    experiment_multi_round_invariant()
