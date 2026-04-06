"""
GF(2) linear algebra using Python integers as bit-vectors.

Each row is a Python int where bit i represents column i.
The rightmost bit (bit 0) is column 0.
For an augmented system Ax=b, bit n (the (n+1)-th bit) is the RHS.
"""


def gf2_rank(rows):
    """Compute rank of a GF(2) matrix (list of ints)."""
    rows = list(rows)
    rank = 0
    used = 0
    for r in rows:
        r ^= (r & used)  # cancel already-used pivots — wrong approach
    # Proper Gaussian elimination
    rows2 = [r for r in rows if r]
    pivots = []
    for r in rows2:
        for p in pivots:
            r = min(r, r ^ p)
        if r:
            pivots.append(r)
    return len(pivots)


def gf2_rank_proper(rows):
    """Proper GF(2) rank via Gaussian elimination."""
    rows = [r for r in rows if r]
    n = 0
    for i in range(len(rows)):
        # find pivot
        pivot_bit = rows[i].bit_length() - 1
        if pivot_bit < 0:
            continue
        # reduce all other rows
        for j in range(len(rows)):
            if j != i and (rows[j] >> pivot_bit) & 1:
                rows[j] ^= rows[i]
        n += 1
    return n


def gf2_gaussian_eliminate(rows, ncols):
    """
    Full Gaussian elimination on GF(2) matrix.
    Returns (echelon_rows, pivot_cols) where echelon_rows is in reduced row echelon form.
    """
    rows = [r for r in rows if r]
    pivot_cols = []

    for col in range(ncols):
        # find row with 1 in this column
        pivot_row = None
        for i in range(len(pivot_cols), len(rows)):
            if (rows[i] >> col) & 1:
                pivot_row = i
                break
        if pivot_row is None:
            continue

        # swap to position
        pos = len(pivot_cols)
        rows[pos], rows[pivot_row] = rows[pivot_row], rows[pos]

        # eliminate column in all other rows
        for i in range(len(rows)):
            if i != pos and (rows[i] >> col) & 1:
                rows[i] ^= rows[pos]

        pivot_cols.append(col)

    return rows[:len(pivot_cols)], pivot_cols


def gf2_solve(rows, ncols):
    """
    Solve GF(2) system. Each row encodes equation: bits 0..ncols-1 are coefficients,
    bit ncols is the RHS (augmented matrix).

    Returns (particular_solution, kernel_basis) or None if inconsistent.
    particular_solution: int (bit-vector of length ncols)
    kernel_basis: list of ints (bit-vectors)
    """
    echelon, pivot_cols = gf2_gaussian_eliminate(list(rows), ncols)

    # Check consistency: any row with all-zero LHS but nonzero RHS
    for row in echelon:
        lhs = row & ((1 << ncols) - 1)
        rhs = (row >> ncols) & 1
        if lhs == 0 and rhs == 1:
            return None  # inconsistent

    # Build particular solution: for each pivot col, read RHS
    particular = 0
    pivot_set = set(pivot_cols)
    for i, col in enumerate(pivot_cols):
        rhs = (echelon[i] >> ncols) & 1
        if rhs:
            particular |= (1 << col)

    # Build kernel basis: one vector per free column
    free_cols = [c for c in range(ncols) if c not in pivot_set]
    kernel = []

    for fc in free_cols:
        vec = 1 << fc  # set the free variable to 1
        # For each pivot, check if it depends on this free col
        for i, pc in enumerate(pivot_cols):
            if (echelon[i] >> fc) & 1:
                vec |= (1 << pc)
        kernel.append(vec)

    return particular, kernel


def gf2_kernel(rows, ncols):
    """Return kernel basis of homogeneous system (no RHS)."""
    # Treat as augmented with RHS=0
    result = gf2_solve(rows, ncols)
    if result is None:
        return []
    return result[1]


def bitvec_to_list(v, n):
    """Convert bit-vector to list of bit values."""
    return [(v >> i) & 1 for i in range(n)]


def list_to_bitvec(bits):
    """Convert list of bit values to bit-vector."""
    v = 0
    for i, b in enumerate(bits):
        if b:
            v |= (1 << i)
    return v


def bitvec_weight(v):
    """Hamming weight of bit-vector."""
    return bin(v).count('1')


def bitvec_dot(a, b):
    """Dot product of two bit-vectors in GF(2)."""
    return bitvec_weight(a & b) & 1
