"""
Iterative Relinearization for the Q∩T system.

Strategy: after Gaussian elimination on the linear part,
substitute pivot expressions into quadratic equations.
This can reveal NEW linear equations (when quadratic terms
cancel or simplify), compressing the kernel further.

Algorithm:
  1. Gaussian eliminate linear part → pivots + free vars
  2. For each pivot var x_p: extract expression x_p = f(free vars)
  3. Substitute pivots into all quadratic equations
  4. Collect any newly-linear equations → add to linear system
  5. Re-Gaussian-eliminate → possibly more pivots
  6. Repeat until no new linear equations found

This is specialized F4/relinearization for the Q∩T structure.
"""

from qt_solver.gf2 import gf2_gaussian_eliminate, gf2_solve, gf2_rank


def iter_bits(x):
    """Yield positions of set bits."""
    while x:
        b = x & (-x)
        yield b.bit_length() - 1
        x ^= b


class PivotMap:
    """
    Stores pivot variable expressions after Gaussian elimination.

    For each pivot variable p:
      x_p = XOR of (free variables in expr_mask) XOR constant
    """

    def __init__(self, pivots, reduced_rows, num_vars):
        self.num_vars = num_vars
        self.is_pivot = set()
        self.expr = {}  # pivot_col -> (free_var_mask, constant)

        pivot_cols = set(col for _, col in pivots)
        self.is_pivot = pivot_cols
        self.free_vars = sorted(set(range(num_vars)) - pivot_cols)
        self.free_set = set(self.free_vars)

        var_mask = (1 << num_vars) - 1
        const_bit = 1 << num_vars

        for row_idx, col in pivots:
            row = reduced_rows[row_idx]
            # Expression: x_col = other bits XOR constant
            # Remove the pivot bit itself
            expr_mask = (row & var_mask) ^ (1 << col)
            constant = (row >> num_vars) & 1
            self.expr[col] = (expr_mask, constant)

    def is_free(self, v):
        return v in self.free_set

    def get_expr(self, v):
        """Get expression for pivot variable v. Returns (mask, const)."""
        return self.expr[v]

    def substitute_var(self, v):
        """
        Substitute variable v.
        If v is pivot: returns (free_var_mask, constant).
        If v is free: returns (1 << v, 0).
        """
        if v in self.expr:
            return self.expr[v]
        return (1 << v, 0)


def substitute_product(pmap, v1, v2):
    """
    Substitute pivot expressions into product v1 * v2.

    Returns: (linear_mask, quadratic_pairs_set, constant)
    where the substituted product equals:
      constant XOR XOR(linear vars) XOR XOR(products)
    """
    m1, c1 = pmap.substitute_var(v1)
    m2, c2 = pmap.substitute_var(v2)

    # (m1 XOR c1) * (m2 XOR c2) over GF(2)
    # = m1*m2 XOR c1*m2 XOR c2*m1 XOR c1*c2

    constant = c1 & c2
    linear = 0
    quad_pairs = set()

    # c1*m2: if c1=1, add m2 as linear
    if c1:
        linear ^= m2

    # c2*m1: if c2=1, add m1 as linear
    if c2:
        linear ^= m1

    # m1*m2: product of two masks = all pairs (i from m1, j from m2)
    # Collect as set of (min(i,j), max(i,j)) pairs
    # Use XOR for cancellation (GF(2))
    bits1 = list(iter_bits(m1))
    bits2 = list(iter_bits(m2))

    for i in bits1:
        for j in bits2:
            if i == j:
                # x_i * x_i = x_i over GF(2)
                linear ^= (1 << i)
            else:
                pair = (min(i, j), max(i, j))
                quad_pairs ^= {pair}

    return linear, quad_pairs, constant


def substitute_quadratic_eq(pmap, lin_row, quad_pairs, num_vars):
    """
    Substitute all pivot expressions into one quadratic equation.

    Input equation: XOR(linear vars) XOR XOR(products) XOR const = 0

    Returns: (new_linear_mask, new_quad_pairs, new_constant)
    expressed purely in terms of free variables.
    """
    var_mask = (1 << num_vars) - 1
    const_bit = 1 << num_vars

    # Start with the linear part
    orig_linear = lin_row & var_mask
    orig_const = (lin_row >> num_vars) & 1

    # Substitute pivot vars in linear part
    result_linear = 0
    result_const = orig_const
    result_quad = set()

    for v in iter_bits(orig_linear):
        m, c = pmap.substitute_var(v)
        result_linear ^= m
        result_const ^= c

    # Substitute pivot vars in quadratic part
    for v1, v2 in quad_pairs:
        lin, qp, c = substitute_product(pmap, v1, v2)
        result_linear ^= lin
        result_quad ^= qp
        result_const ^= c

    return result_linear, result_quad, result_const


def relinearize_once(linear_rows, quad_equations, num_vars, verbose=False):
    """
    One iteration of relinearization.

    Args:
        linear_rows: list of augmented GF(2) row integers
        quad_equations: list of (lin_row, [(v1,v2),...], const)
        num_vars: total variable count

    Returns:
        new_linear_rows: updated linear rows (may have more)
        remaining_quad: quadratic equations that didn't become linear
        stats: dict with iteration statistics
    """
    # Step 1: Gaussian eliminate linear part
    pivots, reduced, rank = gf2_gaussian_eliminate(list(linear_rows), num_vars)
    pmap = PivotMap(pivots, reduced, num_vars)

    if verbose:
        print(f"    Pivots: {len(pivots)}, Free: {len(pmap.free_vars)}, "
              f"Quad eqs: {len(quad_equations)}")

    # Step 2: Substitute into quadratic equations
    new_linear = []
    remaining_quad = []

    for lin_row, qpairs, const in quad_equations:
        res_lin, res_quad, res_const = substitute_quadratic_eq(
            pmap, lin_row, qpairs, num_vars
        )

        if len(res_quad) == 0:
            # Became linear!
            row = res_lin
            if res_const:
                row |= (1 << num_vars)
            if row & ((1 << num_vars) - 1):  # non-trivial
                new_linear.append(row)
        else:
            # Still quadratic — keep
            row = res_lin
            if res_const:
                row |= (1 << num_vars)
            remaining_quad.append((row, list(res_quad), res_const))

    stats = {
        'rank': rank,
        'free_vars': len(pmap.free_vars),
        'new_linear_from_quad': len(new_linear),
        'remaining_quad': len(remaining_quad),
    }

    # Step 3: Add new linear equations and re-eliminate
    all_linear = list(reduced[:rank]) + new_linear

    return all_linear, remaining_quad, stats


def iterative_relinearize(linear_rows, quad_equations, num_vars,
                           max_iterations=20, verbose=True):
    """
    Full iterative relinearization loop.

    Keeps substituting pivots into quadratics until no new
    linear equations are generated.
    """
    current_linear = list(linear_rows)
    current_quad = list(quad_equations)

    history = []

    for iteration in range(max_iterations):
        new_linear, remaining_quad, stats = relinearize_once(
            current_linear, current_quad, num_vars, verbose=verbose
        )

        history.append(stats)

        if verbose:
            print(f"  Iter {iteration}: rank={stats['rank']}, "
                  f"free={stats['free_vars']}, "
                  f"new_linear={stats['new_linear_from_quad']}, "
                  f"remaining_quad={stats['remaining_quad']}")

        if stats['new_linear_from_quad'] == 0:
            if verbose:
                print(f"  Converged at iteration {iteration}")
            break

        current_linear = new_linear
        current_quad = remaining_quad

    # Final Gaussian elimination
    final_pivots, final_reduced, final_rank = gf2_gaussian_eliminate(
        current_linear, num_vars
    )

    return {
        'final_rank': final_rank,
        'final_kernel_dim': num_vars - final_rank,
        'remaining_quad': len(current_quad),
        'iterations': len(history),
        'history': history,
        'quad_equations': current_quad,
    }


def analyze_relinearization(num_rounds, seed=42, verbose=True):
    """
    Run iterative relinearization on the unified v2 Q∩T system
    and report kernel compression.
    """
    import random
    from qt_solver.sha256_traced import MASK32
    from qt_solver.unified_v2 import build_unified_v2

    rng = random.Random(seed)
    msg = [rng.randint(0, MASK32) for _ in range(16)]

    system, trace = build_unified_v2(msg, num_rounds)
    st = system.stats()
    n = system.vm.num_vars

    if verbose:
        print(f"\n{'='*60}")
        print(f"Relinearization: R={num_rounds}")
        print(f"{'='*60}")
        print(f"Vars: {n}, Linear: {st['num_linear_eq']}, Quad: {st['num_quad_eq']}")

    # Initial rank (linear only)
    init_rank = gf2_rank(list(system.linear_rows), n)
    init_kernel = n - init_rank

    if verbose:
        print(f"Initial: rank={init_rank}, kernel={init_kernel}")
        print(f"\nRelinearization iterations:")

    result = iterative_relinearize(
        system.linear_rows, system.quad_equations, n,
        verbose=verbose
    )

    compression = init_kernel - result['final_kernel_dim']

    if verbose:
        print(f"\nResult:")
        print(f"  Initial kernel: {init_kernel}")
        print(f"  Final kernel:   {result['final_kernel_dim']}")
        print(f"  Compression:    {compression} bits")
        print(f"  Remaining quad: {result['remaining_quad']}")

    return {
        'num_rounds': num_rounds,
        'num_vars': n,
        'initial_rank': init_rank,
        'initial_kernel': init_kernel,
        'final_rank': result['final_rank'],
        'final_kernel': result['final_kernel_dim'],
        'compression': compression,
        'remaining_quad': result['remaining_quad'],
        'iterations': result['iterations'],
    }
