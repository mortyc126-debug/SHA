"""
GF(2) linear algebra using Python arbitrary-precision integers as bit vectors.

Each equation is represented as a big integer where:
  - Bits 0..n-1 are coefficients of variables x_0..x_{n-1}
  - Bit n is the constant (augmented column)

XOR of two equations = Python ^ on their integer representations.
"""


def gf2_gaussian_eliminate(rows, num_vars):
    """
    Gaussian elimination over GF(2) on augmented system.

    Args:
        rows: list of integers (augmented bit vectors, bit num_vars = constant)
        num_vars: number of variables

    Returns:
        pivots: list of (row_index, pivot_column)
        reduced: list of row integers in reduced row echelon form
        rank: rank of the coefficient matrix
    """
    rows = list(rows)  # copy
    n = len(rows)
    pivots = []
    pivot_row = 0

    for col in range(num_vars):
        # Find a row with a 1 in this column
        found = -1
        for i in range(pivot_row, n):
            if (rows[i] >> col) & 1:
                found = i
                break
        if found == -1:
            continue

        # Swap to pivot position
        rows[pivot_row], rows[found] = rows[found], rows[pivot_row]

        # Eliminate column in all other rows
        pivot_val = rows[pivot_row]
        for i in range(n):
            if i != pivot_row and (rows[i] >> col) & 1:
                rows[i] ^= pivot_val

        pivots.append((pivot_row, col))
        pivot_row += 1

    return pivots, rows, pivot_row


def gf2_solve(rows, num_vars):
    """
    Solve GF(2) linear system Ax = b (augmented).

    Args:
        rows: list of augmented row integers
        num_vars: number of variables

    Returns:
        None if inconsistent (no solution)
        (particular, kernel) where:
            particular: integer bit vector (one solution)
            kernel: list of integer bit vectors (null space basis)
    """
    pivots, reduced, rank = gf2_gaussian_eliminate(rows, num_vars)

    # Check consistency: any row with all-zero coefficients but nonzero constant?
    const_bit = 1 << num_vars
    for i in range(rank, len(reduced)):
        if reduced[i] & const_bit:
            return None  # Inconsistent

    pivot_cols = set(col for _, col in pivots)
    free_cols = [c for c in range(num_vars) if c not in pivot_cols]

    # Particular solution: set all free variables to 0
    particular = 0
    for row_idx, col in pivots:
        if (reduced[row_idx] >> num_vars) & 1:
            particular |= (1 << col)

    # Kernel basis: for each free variable, set it to 1
    kernel = []
    for fc in free_cols:
        vec = 1 << fc
        for row_idx, col in pivots:
            if (reduced[row_idx] >> fc) & 1:
                vec |= (1 << col)
        kernel.append(vec)

    return particular, kernel


def gf2_rank(rows, num_vars):
    """Compute rank of a GF(2) matrix."""
    _, _, rank = gf2_gaussian_eliminate(rows, num_vars)
    return rank


def gf2_mat_vec_mul(rows, vec, num_vars):
    """Multiply GF(2) matrix by vector. Returns result vector as integer."""
    result = 0
    for i, row in enumerate(rows):
        # Dot product of row and vec over GF(2)
        bits = row & vec & ((1 << num_vars) - 1)
        parity = bin(bits).count('1') & 1
        result |= (parity << i)
    return result


def popcount(x):
    """Count number of 1 bits."""
    return bin(x).count('1')


def bit_to_list(x, n):
    """Convert bit vector to list of set bit positions."""
    return [i for i in range(n) if (x >> i) & 1]


def list_to_bit(positions):
    """Convert list of positions to bit vector."""
    result = 0
    for p in positions:
        result |= (1 << p)
    return result
