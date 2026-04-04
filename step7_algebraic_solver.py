"""
Step 7: Algebraic Solver — SOLVE the equations, don't just count solutions

We know: a_new[0] = T1[0] ⊕ T2[0] (degree 2, carry-free).
The collision equation: f(M) ⊕ f(M⊕δ) = 0.

For a degree-d polynomial, f(M⊕δ) ⊕ f(M) is a polynomial in δ (of degree d).
For bit 0 (degree 2): the collision equation is degree 2 in δ.

Method:
1. Compute f(M⊕δ) ⊕ f(M) as a truth table in δ
2. Convert to ANF in δ → get explicit polynomial equations
3. Linearize: x_i·x_j → y_{ij}
4. Solve by Gaussian elimination
5. Recover solutions

This is SOLVING, not brute-forcing.
"""

import numpy as np
from collections import defaultdict
from step0_exact_algebra import mini_sha, N, MASK, mobius_transform

N_MSG = 4
N_INPUT = N * N_MSG  # 16
N_TOTAL = 1 << N_INPUT


def compute_collision_polynomial(M, R, out_bit):
    """
    Compute the collision polynomial: g(δ) = f(M⊕δ)[bit] ⊕ f(M)[bit]
    where f = mini_sha output, as a truth table indexed by δ.
    Returns: truth table of g(δ) for all 2^16 values of δ.
    """
    a_base, e_base = mini_sha(M, R)

    # For output bit: 0..N-1 = a bits, N..2N-1 = e bits
    if out_bit < N:
        base_val = (a_base >> out_bit) & 1
    else:
        base_val = (e_base >> (out_bit - N)) & 1

    tt = np.zeros(N_TOTAL, dtype=np.uint8)

    for delta_idx in range(N_TOTAL):
        delta_words = []
        tmp = delta_idx
        for w in range(N_MSG):
            delta_words.append(tmp & MASK)
            tmp >>= N

        M2 = [(M[w] ^ delta_words[w]) for w in range(N_MSG)]
        a2, e2 = mini_sha(M2, R)

        if out_bit < N:
            val2 = (a2 >> out_bit) & 1
        else:
            val2 = (e2 >> (out_bit - N)) & 1

        tt[delta_idx] = base_val ^ val2

    return tt


def anf_to_equations(anf, n_vars):
    """
    Convert ANF to list of monomials.
    Each monomial = frozenset of variable indices.
    """
    monomials = []
    for idx in range(len(anf)):
        if anf[idx]:
            variables = frozenset(v for v in range(n_vars) if idx & (1 << v))
            monomials.append(variables)
    return monomials


def linearize_and_solve(equations, n_vars, max_degree=2):
    """
    Linearize a system of GF(2) polynomial equations.

    Each equation = list of monomials (frozensets of variable indices).
    Equation = 0 means XOR of all monomials = 0.

    Linearization: each monomial up to max_degree becomes a new variable.
    Solve by Gaussian elimination.

    Returns: kernel dimension and number of solutions.
    """
    # Enumerate all monomials up to degree max_degree
    from itertools import combinations

    mono_to_idx = {}
    idx_to_mono = []

    # Constant (empty set)
    mono_to_idx[frozenset()] = 0
    idx_to_mono.append(frozenset())

    # Degree 1
    for v in range(n_vars):
        m = frozenset([v])
        mono_to_idx[m] = len(idx_to_mono)
        idx_to_mono.append(m)

    # Degree 2
    if max_degree >= 2:
        for v1 in range(n_vars):
            for v2 in range(v1+1, n_vars):
                m = frozenset([v1, v2])
                mono_to_idx[m] = len(idx_to_mono)
                idx_to_mono.append(m)

    # Degree 3 (if needed)
    if max_degree >= 3:
        for combo in combinations(range(n_vars), 3):
            m = frozenset(combo)
            mono_to_idx[m] = len(idx_to_mono)
            idx_to_mono.append(m)

    n_monos = len(idx_to_mono)
    n_eqs = len(equations)

    print(f"    Linearization: {n_vars} vars, degree ≤ {max_degree}, "
          f"{n_monos} linearized variables, {n_eqs} equations")

    # Build matrix
    matrix = np.zeros((n_eqs, n_monos), dtype=np.uint8)

    higher_degree_terms = 0
    for i, eq in enumerate(equations):
        for mono in eq:
            if len(mono) <= max_degree and mono in mono_to_idx:
                matrix[i][mono_to_idx[mono]] = 1
            else:
                higher_degree_terms += 1

    if higher_degree_terms > 0:
        print(f"    WARNING: {higher_degree_terms} terms with degree > {max_degree} IGNORED")

    # Gaussian elimination over GF(2)
    rank = gf2_gaussian_rank(matrix)
    kernel_dim = n_monos - rank

    print(f"    Rank: {rank}/{n_eqs}, Kernel dimension: {kernel_dim}")
    print(f"    Solution space: 2^{kernel_dim} = {2**kernel_dim if kernel_dim < 30 else '2^'+str(kernel_dim)}")

    return rank, kernel_dim, n_monos


def gf2_gaussian_rank(matrix):
    """Compute rank of binary matrix over GF(2)."""
    m = matrix.copy()
    rows, cols = m.shape
    rank = 0
    for col in range(cols):
        # Find pivot
        pivot = None
        for row in range(rank, rows):
            if m[row][col]:
                pivot = row
                break
        if pivot is None:
            continue
        # Swap
        m[[rank, pivot]] = m[[pivot, rank]]
        # Eliminate
        for row in range(rows):
            if row != rank and m[row][col]:
                m[row] ^= m[rank]
        rank += 1
    return rank


def solve_collision_algebraically(R, M=None):
    """
    Extract collision equations and solve them algebraically.
    """
    import random
    if M is None:
        rng = random.Random(42)
        M = [rng.randint(0, MASK) for _ in range(N_MSG)]

    print(f"\n{'='*80}")
    print(f"ALGEBRAIC COLLISION SOLVER — R={R}")
    print(f"Message M = {[hex(w) for w in M]}")
    print(f"{'='*80}")

    # Compute collision polynomial for each output bit
    all_equations = []
    bit_degrees = []

    for out_bit in range(2 * N):
        reg = "a" if out_bit < N else "e"
        bit_idx = out_bit % N

        tt_collision = compute_collision_polynomial(M, R, out_bit)
        anf_collision = mobius_transform(tt_collision, N_INPUT)
        monomials = anf_to_equations(anf_collision, N_INPUT)

        # Get max degree
        max_deg = max(len(m) for m in monomials) if monomials else 0
        n_monos = len(monomials)

        # Check: constant term should be 0 (δ=0 gives collision trivially)
        has_constant = frozenset() in monomials
        if has_constant:
            monomials.remove(frozenset())

        bit_degrees.append(max_deg)
        all_equations.append(monomials)

        print(f"\n  {reg}[{bit_idx}] collision equation: "
              f"degree={max_deg}, {n_monos} monomials, const={has_constant}")

    # ================================================================
    # SOLVE LAYER BY LAYER
    # ================================================================

    print(f"\n{'='*80}")
    print(f"LAYERED ALGEBRAIC SOLVING")
    print(f"{'='*80}")

    # Layer 0: bit 0 of a and e (carry-free, degree 2)
    print(f"\n  --- LAYER 0: a[0] and e[0] (degree {bit_degrees[0]}, {bit_degrees[N]}) ---")
    layer0_eqs = [all_equations[0], all_equations[N]]  # a[0], e[0]
    r0, k0, n0 = linearize_and_solve(layer0_eqs, N_INPUT, max_degree=2)

    # Layer 0+1: add bit 1
    print(f"\n  --- LAYER 0+1: add a[1] and e[1] (degree {bit_degrees[1]}, {bit_degrees[N+1]}) ---")
    layer01_eqs = [all_equations[0], all_equations[N],
                    all_equations[1], all_equations[N+1]]
    r01, k01, n01 = linearize_and_solve(layer01_eqs, N_INPUT, max_degree=2)

    # All equations at degree 2 (ignoring higher terms)
    print(f"\n  --- ALL BITS at degree ≤ 2 (ignoring higher-degree terms) ---")
    r_all2, k_all2, n_all2 = linearize_and_solve(all_equations, N_INPUT, max_degree=2)

    # All equations at degree 3
    print(f"\n  --- ALL BITS at degree ≤ 3 ---")
    r_all3, k_all3, n_all3 = linearize_and_solve(all_equations, N_INPUT, max_degree=3)

    # Verify by brute force
    print(f"\n  --- BRUTE FORCE VERIFICATION ---")
    a_base, e_base = mini_sha(M, R)
    n_collisions = 0
    for delta_idx in range(1, N_TOTAL):
        delta_words = []
        tmp = delta_idx
        for w in range(N_MSG):
            delta_words.append(tmp & MASK)
            tmp >>= N
        M2 = [(M[w] ^ delta_words[w]) for w in range(N_MSG)]
        a2, e2 = mini_sha(M2, R)
        if a_base == a2 and e_base == e2:
            n_collisions += 1

    print(f"    Actual collisions (brute force): {n_collisions}")
    print(f"    Algebraic prediction (layer 0, 2^{k0}): {2**k0}")
    print(f"    Algebraic prediction (all bits, deg≤2, 2^{k_all2}): {2**k_all2}")
    print(f"    Algebraic prediction (all bits, deg≤3, 2^{k_all3}): {2**k_all3}")

    # ================================================================
    # KEY COMPARISON
    # ================================================================
    print(f"\n{'='*80}")
    print(f"ALGEBRAIC SOLVER vs BRUTE FORCE")
    print(f"{'='*80}")
    print(f"  Brute force: enumerate 2^{N_INPUT} = {N_TOTAL} deltas → {n_collisions} solutions")
    print(f"  Birthday:    ~2^{N*2//2} = {2**(N)} evaluations (for {2*N}-bit output)")
    print(f"")
    print(f"  Algebraic (degree ≤ 2 linearization):")
    print(f"    Variables: {n0}, Rank: {r_all2}, Kernel: {k_all2}")
    print(f"    Solution space: 2^{k_all2} = {2**k_all2}")
    print(f"    Includes higher-degree solutions that DON'T satisfy ignored terms")
    print(f"")
    print(f"  Algebraic (degree ≤ 3 linearization):")
    print(f"    Variables: {n_all3}, Rank: {r_all3}, Kernel: {k_all3}")
    print(f"    Solution space: 2^{k_all3} = {2**k_all3}")
    print(f"")
    print(f"  If algebraic kernel ≈ brute force count → linearization captures everything")
    print(f"  If algebraic kernel >> brute force → need higher degree to filter")

    return n_collisions, k0, k_all2, k_all3


def main():
    for R in [4, 6, 8]:
        solve_collision_algebraically(R)


if __name__ == "__main__":
    main()
