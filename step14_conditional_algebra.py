"""
Step 14: Conditional Algebra — Ch as MUX, not multiplication

Problem: Ch(e,f,g) = e·(f⊕g) ⊕ g = degree 2 in GF(2).
Through rounds: degree doubles (Fibonacci) → explosion.

But Ch is conceptually a MUX: "if e then f else g".
A MUX SELECTS — it doesn't multiply degrees.

Approach 1: LINEARIZATION at a known base point.
  ΔCh ≈ (f₀⊕g₀)·Δe ⊕ e₀·Δf ⊕ (1⊕e₀)·Δg
  This is DEGREE 1 in differences, with KNOWN coefficients from base state.
  If we linearize ALL Ch/Maj → entire SHA becomes a LINEAR system.
  Solve by Gaussian elimination, then check/correct quadratic errors.

Approach 2: CONDITIONAL REPRESENTATION.
  Represent f as a DECISION TREE, not a polynomial.
  Ch(e,f,g) adds ONE node to the tree (not doubling the degree).
  Tree SIZE grows linearly per round, not exponentially.

Let's test BOTH on mini-SHA and see which works better.
"""

import numpy as np
import time
from step0_exact_algebra import mini_sha, N, MASK, mobius_transform

N_MSG = 4
N_INPUT = N * N_MSG
N_TOTAL = 1 << N_INPUT


def linearized_collision_solver(M_base, R):
    """
    Approach 1: Linearize Ch and Maj at the base point M_base.

    For the collision system f(M⊕δ) ⊕ f(M) = 0:
    Replace Ch(e⊕Δe, f⊕Δf, g⊕Δg) ⊕ Ch(e,f,g) with its LINEAR part:
      ΔCh ≈ (f⊕g)·Δe ⊕ e·Δf ⊕ ē·Δg

    This makes the ENTIRE system linear in δM.
    Solve by Gaussian elimination.
    Then check which solutions are TRUE collisions.
    """
    print(f"\n{'='*80}")
    print(f"LINEARIZED COLLISION SOLVER — R={R}")
    print(f"{'='*80}")

    a_base, e_base = mini_sha(M_base, R)

    # Build the linearized system:
    # For each output bit, compute the linear approximation of
    # f(M⊕δ) ⊕ f(M) as a function of δ.
    #
    # Method: compute the JACOBIAN — derivative of each output bit
    # with respect to each input bit.

    # Jacobian: J[i][j] = ∂(output_bit_i) / ∂(input_bit_j) over GF(2)
    n_out = 2 * N  # 8 output bits
    J = np.zeros((n_out, N_INPUT), dtype=np.uint8)
    c = np.zeros(n_out, dtype=np.uint8)  # constant term (should be 0 at δ=0)

    for j in range(N_INPUT):
        # Flip input bit j
        delta = 1 << j
        dw = []
        tmp = delta
        for w in range(N_MSG):
            dw.append(tmp & MASK)
            tmp >>= N
        M2 = [(M_base[w] ^ dw[w]) for w in range(N_MSG)]
        a2, e2 = mini_sha(M2, R)

        for i in range(N):
            J[i][j] = ((a_base ^ a2) >> i) & 1
            J[N+i][j] = ((e_base ^ e2) >> i) & 1

    # The linear system: J · δ = 0 (collision means output difference = 0)
    # Solve for kernel of J

    print(f"  Jacobian: {n_out} × {N_INPUT}")
    rank = gf2_rank(J)
    kernel_dim = N_INPUT - rank
    print(f"  Rank: {rank}")
    print(f"  Kernel dimension: {kernel_dim}")
    print(f"  Linear kernel: 2^{kernel_dim} = {2**kernel_dim} candidate δ")

    # Find the kernel basis
    kernel_basis = gf2_kernel(J)
    print(f"  Kernel basis vectors: {len(kernel_basis)}")

    # Enumerate kernel and check TRUE collisions
    t0 = time.time()
    n_kernel = 2 ** kernel_dim
    true_collisions = 0
    false_positives = 0

    if kernel_dim <= 20:  # enumerate if small enough
        for k_idx in range(n_kernel):
            # Build δ from kernel basis
            delta_bits = np.zeros(N_INPUT, dtype=np.uint8)
            for b in range(kernel_dim):
                if (k_idx >> b) & 1:
                    delta_bits ^= kernel_basis[b]

            # Convert to message words
            delta_int = 0
            for j in range(N_INPUT):
                if delta_bits[j]:
                    delta_int |= (1 << j)

            if delta_int == 0:
                continue

            dw = []
            tmp = delta_int
            for w in range(N_MSG):
                dw.append(tmp & MASK)
                tmp >>= N
            M2 = [(M_base[w] ^ dw[w]) for w in range(N_MSG)]
            a2, e2 = mini_sha(M2, R)

            if a2 == a_base and e2 == e_base:
                true_collisions += 1
            else:
                false_positives += 1

        t_solve = time.time() - t0

        # Brute force count for comparison
        brute_collisions = 0
        t0 = time.time()
        for didx in range(1, N_TOTAL):
            dw = []
            tmp = didx
            for w in range(N_MSG):
                dw.append(tmp & MASK)
                tmp >>= N
            M2 = [(M_base[w] ^ dw[w]) for w in range(N_MSG)]
            a2, e2 = mini_sha(M2, R)
            if a2 == a_base and e2 == e_base:
                brute_collisions += 1
        t_brute = time.time() - t0

        print(f"\n  RESULTS:")
        print(f"  Linear kernel size: {n_kernel}")
        print(f"  True collisions in kernel: {true_collisions}")
        print(f"  False positives: {false_positives}")
        print(f"  Precision: {true_collisions/max(n_kernel-1,1)*100:.1f}%")
        print(f"  Total collisions (brute force): {brute_collisions}")
        print(f"  Recall: {true_collisions/max(brute_collisions,1)*100:.1f}%")
        print(f"  Solve time: {t_solve:.3f}s")
        print(f"  Brute force time: {t_brute:.3f}s")
        print(f"  Speedup: {t_brute/max(t_solve,0.001):.1f}×")
        print(f"  Birthday cost: ~{2**N} = {2**N}")

        # KEY METRIC: efficiency
        evals_linear = n_kernel  # enumerate kernel
        print(f"\n  EFFICIENCY:")
        print(f"  Linear: {evals_linear} evaluations → {true_collisions} collisions")
        print(f"  Brute:  {N_TOTAL} evaluations → {brute_collisions} collisions")
        print(f"  Birthday: ~{2**N} evaluations → ~1 collision")
        if true_collisions > 0:
            print(f"  Linear cost per collision: {evals_linear/true_collisions:.0f}")
            print(f"  Brute cost per collision: {N_TOTAL/brute_collisions:.0f}")
            print(f"  Birthday cost per collision: ~{2**N}")

        return true_collisions, brute_collisions, kernel_dim

    else:
        print(f"  Kernel too large to enumerate (2^{kernel_dim})")
        return 0, 0, kernel_dim


def gf2_rank(M):
    """Compute rank of binary matrix over GF(2)."""
    m = M.copy()
    rows, cols = m.shape
    rank = 0
    for col in range(cols):
        pivot = None
        for row in range(rank, rows):
            if m[row][col]:
                pivot = row
                break
        if pivot is None:
            continue
        m[[rank, pivot]] = m[[pivot, rank]]
        for row in range(rows):
            if row != rank and m[row][col]:
                m[row] ^= m[rank]
        rank += 1
    return rank


def gf2_kernel(M):
    """Find kernel basis of binary matrix over GF(2)."""
    m = M.copy()
    rows, cols = m.shape

    # Augment with identity for tracking
    aug = np.zeros((rows, cols + cols), dtype=np.uint8)
    aug[:, :cols] = m

    # Track column operations instead
    # Use null space computation
    # Transpose and find null space of M^T = find left null space
    # Actually, find right null space: Mx = 0

    # RREF
    m = M.copy()
    pivot_cols = []
    rank = 0
    for col in range(cols):
        pivot = None
        for row in range(rank, rows):
            if m[row][col]:
                pivot = row
                break
        if pivot is None:
            continue
        m[[rank, pivot]] = m[[pivot, rank]]
        for row in range(rows):
            if row != rank and m[row][col]:
                m[row] ^= m[rank]
        pivot_cols.append(col)
        rank += 1

    # Free columns
    free_cols = [c for c in range(cols) if c not in pivot_cols]

    # Build kernel basis
    kernel = []
    for fc in free_cols:
        vec = np.zeros(cols, dtype=np.uint8)
        vec[fc] = 1
        for i, pc in enumerate(pivot_cols):
            if m[i][fc]:
                vec[pc] = 1
        kernel.append(vec)

    return kernel


def main():
    import random
    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(N_MSG)]

    for R in [3, 4, 5, 6, 8]:
        linearized_collision_solver(M, R)


if __name__ == "__main__":
    main()
