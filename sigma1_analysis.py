"""
Analysis of SHA-256 Sigma1, Sigma0, sigma0, sigma1 linear operators over GF(2).
Key question: are characteristic polynomials irreducible?
If irreducible -> no invariant subspaces -> SHA-256 optimal
If reducible with small factor -> invariant subspace -> potential weakness
"""

import numpy as np
from functools import reduce

N = 32  # word size

def rotr_matrix(k, n=N):
    """ROTR_k as n×n matrix over GF(2)"""
    M = np.zeros((n, n), dtype=int)
    for i in range(n):
        M[i][(i + k) % n] = 1
    return M

def shr_matrix(k, n=N):
    """SHR_k as n×n matrix over GF(2)"""
    M = np.zeros((n, n), dtype=int)
    for i in range(k, n):
        M[i - k][i] = 1  # bit i goes to position i-k
    # Wait, SHR_k: bit i of input -> bit (i-k) of output, for i >= k
    # Actually: (x >> k)[j] = x[j+k] for j+k < n, else 0
    M2 = np.zeros((n, n), dtype=int)
    for j in range(n):
        if j + k < n:
            M2[j][j + k] = 1
    return M2

def xor_matrices(*matrices):
    """XOR (add mod 2) matrices"""
    result = matrices[0].copy()
    for m in matrices[1:]:
        result = (result + m) % 2
    return result

# SHA-256 linear operators
Sigma1 = xor_matrices(rotr_matrix(6), rotr_matrix(11), shr_matrix(25))
Sigma0 = xor_matrices(rotr_matrix(2), rotr_matrix(13), rotr_matrix(22))
sigma0 = xor_matrices(rotr_matrix(7), rotr_matrix(18), shr_matrix(3))
sigma1 = xor_matrices(rotr_matrix(17), rotr_matrix(19), shr_matrix(10))

def gf2_rank(M):
    """Rank of matrix over GF(2) via Gaussian elimination"""
    m = M.copy() % 2
    rows, cols = m.shape
    rank = 0
    for col in range(cols):
        pivot = None
        for row in range(rank, rows):
            if m[row][col] == 1:
                pivot = row
                break
        if pivot is None:
            continue
        m[[rank, pivot]] = m[[pivot, rank]]
        for row in range(rows):
            if row != rank and m[row][col] == 1:
                m[row] = (m[row] + m[rank]) % 2
        rank += 1
    return rank

def gf2_charpoly(M):
    """
    Characteristic polynomial of M over GF(2).
    Returns list of coefficients [c0, c1, ..., cn] where poly = c0 + c1*x + ... + cn*x^n
    Uses Berkowitz algorithm or direct computation for small n.
    For 32x32: use the fact that charpoly over GF(2) = det(M + x*I) mod 2.
    We compute it via the Faddeev-LeVerrier algorithm adapted for GF(2).
    """
    n = M.shape[0]
    # Over GF(2), we can compute charpoly by finding minimal polynomial
    # or by direct matrix operations.
    # Method: compute M^k for k=0..n, find linear dependency over GF(2)

    # Actually, let's use a direct method:
    # charpoly coefficients via traces of powers (Newton's identities)
    # But over GF(2), trace is just sum of diagonal mod 2

    # Better: use Hessenberg form or Frobenius form
    # For GF(2), let's use the Danilevsky method

    # Simplest for correctness: compute det(xI + M) symbolically
    # For 32x32 over GF(2) this is feasible via Gaussian elimination on polynomial matrix

    # Use Berkowitz algorithm (division-free)
    # Or: compute via factored form

    # Let's use a practical approach: find the minimal polynomial first
    # by computing M^0, M^1, ..., M^n vectors on a random start vector

    # For the FULL characteristic polynomial, use:
    # Build the matrix xI + M over GF(2)[x], compute determinant

    # We'll use the approach of computing the sequence of principal minors
    # Actually, for GF(2), let's just use the companion matrix / rational canonical form

    # Simplest correct method: Frobenius normal form via Krylov
    # Start with random vector v, compute v, Mv, M²v, ... until linear dependence

    import random
    random.seed(42)

    # Find minimal polynomial on random vector (= charpoly if cyclic)
    v = np.array([random.randint(0,1) for _ in range(n)], dtype=int)

    # Krylov sequence
    krylov = [v.copy()]
    current = v.copy()
    for i in range(n):
        current = M.dot(current) % 2
        krylov.append(current.copy())

    # krylov[0], ..., krylov[n]: find first linear dependence
    # Build matrix [v | Mv | ... | M^n v] and find kernel
    K = np.column_stack(krylov) % 2  # n × (n+1)

    # Find the first k such that krylov[k] is in span of krylov[0..k-1]
    # via Gaussian elimination on K^T
    KT = K.T.copy() % 2  # (n+1) × n

    # Row reduce to find dependencies
    pivot_cols = []
    mat = KT.copy()
    rows_m, cols_m = mat.shape
    row_idx = 0
    for col in range(cols_m):
        pivot = None
        for r in range(row_idx, rows_m):
            if mat[r][col] == 1:
                pivot = r
                break
        if pivot is None:
            continue
        mat[[row_idx, pivot]] = mat[[pivot, row_idx]]
        for r in range(rows_m):
            if r != row_idx and mat[r][col] == 1:
                mat[r] = (mat[r] + mat[row_idx]) % 2
        pivot_cols.append(col)
        row_idx += 1

    rank_K = len(pivot_cols)

    # The minimal polynomial has degree = rank of Krylov matrix
    # If rank = n, the minimal polynomial = characteristic polynomial (cyclic vector)

    # Find the dependency: which row of reduced matrix is zero that gives us the relation
    # Actually let's find the minimal poly differently

    # Build the n × (n+1) matrix and find the leftmost dependency
    # [v, Mv, M²v, ..., M^n v] - find coefficients c0,...,cn with Σ ci M^i v = 0

    # Transpose approach: solve for last column in terms of previous
    A = np.column_stack(krylov[:n]) % 2  # n × n
    b = krylov[n] % 2

    # Solve Ax = b over GF(2) (x = coefficients of minimal poly)
    aug = np.column_stack([A, b.reshape(-1,1)]) % 2
    aug_rows, aug_cols = aug.shape

    pivot_r = 0
    pivots = {}
    for col in range(aug_cols - 1):
        found = None
        for r in range(pivot_r, aug_rows):
            if aug[r][col] == 1:
                found = r
                break
        if found is None:
            continue
        aug[[pivot_r, found]] = aug[[found, pivot_r]]
        for r in range(aug_rows):
            if r != pivot_r and aug[r][col] == 1:
                aug[r] = (aug[r] + aug[pivot_r]) % 2
        pivots[col] = pivot_r
        pivot_r += 1

    # Extract solution
    x = np.zeros(n, dtype=int)
    for col, row in pivots.items():
        x[col] = aug[row][-1]

    # Minimal poly: M^n v = Σ x[i] M^i v => M^n + x[n-1]M^{n-1} + ... + x[0]I = 0 on v
    # Coefficients: [x[0], x[1], ..., x[n-1], 1] (monic polynomial of degree n)
    minpoly = list(x) + [1]

    return minpoly

def poly_to_str(coeffs):
    """Pretty print polynomial over GF(2)"""
    terms = []
    for i, c in enumerate(coeffs):
        if c == 1:
            if i == 0:
                terms.append("1")
            elif i == 1:
                terms.append("x")
            else:
                terms.append(f"x^{i}")
    return " + ".join(terms) if terms else "0"

def gf2_poly_gcd(a, b):
    """GCD of two polynomials over GF(2)"""
    while any(x == 1 for x in b):
        # Polynomial division over GF(2)
        while len(a) >= len(b) and any(x == 1 for x in a):
            if a[-1] == 0:
                a = a[:-1]
                continue
            shift = len(a) - len(b)
            for i in range(len(b)):
                a[i + shift] ^= b[i]
            while len(a) > 0 and a[-1] == 0:
                a = a[:-1]
        a, b = b, a
    return a if any(x == 1 for x in a) else [0]

def gf2_poly_mul(a, b):
    """Multiply two polynomials over GF(2)"""
    if not a or not b:
        return [0]
    result = [0] * (len(a) + len(b) - 1)
    for i, ca in enumerate(a):
        if ca:
            for j, cb in enumerate(b):
                if cb:
                    result[i + j] ^= 1
    # Trim trailing zeros
    while len(result) > 1 and result[-1] == 0:
        result = result[:-1]
    return result

def gf2_poly_mod(a, m):
    """a mod m over GF(2)"""
    a = list(a)
    m = list(m)
    while len(a) >= len(m):
        if a[-1] == 0:
            a = a[:-1]
            continue
        shift = len(a) - len(m)
        for i in range(len(m)):
            a[i + shift] ^= m[i]
        while len(a) > 1 and a[-1] == 0:
            a = a[:-1]
    return a

def gf2_poly_powmod(base, exp, mod):
    """base^exp mod mod over GF(2)"""
    result = [1]  # 1
    base = gf2_poly_mod(base, mod)
    while exp > 0:
        if exp & 1:
            result = gf2_poly_mod(gf2_poly_mul(result, base), mod)
        base = gf2_poly_mod(gf2_poly_mul(base, base), mod)
        exp >>= 1
    return result

def is_irreducible_gf2(poly):
    """Test if polynomial is irreducible over GF(2)"""
    n = len(poly) - 1  # degree
    if n <= 0:
        return False
    if poly[0] == 0:  # Must have nonzero constant term
        return False

    # Check: x^(2^k) ≡ x (mod poly) for k = n, and
    # gcd(x^(2^(n/p)) - x, poly) = 1 for each prime p | n

    # First check: x^(2^n) ≡ x mod poly
    x_poly = [0, 1]  # x
    x_pow = gf2_poly_powmod(x_poly, 2**n, poly)
    diff = list(x_pow)
    if len(diff) <= 1:
        diff = diff + [0] * (2 - len(diff))
    if len(diff) > 1:
        diff[1] ^= 1  # subtract x
    else:
        diff.append(1)
    while len(diff) > 1 and diff[-1] == 0:
        diff = diff[:-1]

    if any(d == 1 for d in diff):
        return False  # x^(2^n) ≠ x mod poly

    # Check divisors: for each prime p dividing n, gcd(x^(2^(n/p)) - x, poly) should be 1
    # Find prime factors of n
    primes = set()
    temp = n
    for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]:
        while temp % p == 0:
            primes.add(p)
            temp //= p
    if temp > 1:
        primes.add(temp)

    for p in primes:
        k = n // p
        x_pow_k = gf2_poly_powmod(x_poly, 2**k, poly)
        diff_k = list(x_pow_k)
        if len(diff_k) <= 1:
            diff_k = diff_k + [0] * (2 - len(diff_k))
        if len(diff_k) > 1:
            diff_k[1] ^= 1
        else:
            diff_k.append(1)
        while len(diff_k) > 1 and diff_k[-1] == 0:
            diff_k = diff_k[:-1]

        g = gf2_poly_gcd(list(diff_k), list(poly))
        while len(g) > 1 and g[-1] == 0:
            g = g[:-1]
        if len(g) > 1:  # gcd has degree > 0
            return False

    return True

def factorize_gf2(poly):
    """
    Factor polynomial over GF(2) using Berlekamp or simple trial division.
    Returns list of (factor, multiplicity) pairs.
    """
    n = len(poly) - 1
    if n <= 1:
        return [(poly, 1)]

    if is_irreducible_gf2(poly):
        return [(poly, 1)]

    # Try to find factors by computing gcd with x^(2^k) - x for k = 1, 2, ...
    factors = []
    remaining = list(poly)

    x_poly = [0, 1]

    for k in range(1, n // 2 + 1):
        if len(remaining) - 1 < 2 * k:
            break
        # x^(2^k) - x mod remaining
        x_pow = gf2_poly_powmod(x_poly, 2**k, remaining)
        diff = list(x_pow)
        if len(diff) <= 1:
            diff = diff + [0] * (2 - len(diff))
        if len(diff) > 1:
            diff[1] ^= 1
        else:
            diff.append(1)
        while len(diff) > 1 and diff[-1] == 0:
            diff = diff[:-1]

        g = gf2_poly_gcd(list(diff), list(remaining))
        while len(g) > 1 and g[-1] == 0:
            g = g[:-1]

        if 1 < len(g) - 0 < len(remaining):
            deg_g = len(g) - 1
            if deg_g > 0 and deg_g < len(remaining) - 1:
                # g is a non-trivial factor
                # Extract factors of degree k from g
                sub_factors = factorize_gf2(g)
                factors.extend(sub_factors)
                # Divide remaining by g
                # Polynomial division over GF(2)
                quotient = [0] * (len(remaining) - len(g) + 1)
                temp = list(remaining)
                for i in range(len(quotient) - 1, -1, -1):
                    if len(temp) - 1 >= i + len(g) - 1 and temp[i + len(g) - 1] == 1:
                        quotient[i] = 1
                        for j in range(len(g)):
                            temp[i + j] ^= g[j]
                while len(quotient) > 1 and quotient[-1] == 0:
                    quotient = quotient[:-1]
                remaining = quotient

    if len(remaining) > 1:
        sub = factorize_gf2(remaining)
        factors.extend(sub)
    elif remaining == [1]:
        pass  # unit

    if not factors:
        factors = [(poly, 1)]

    return factors

# ============================================================
# MAIN COMPUTATION
# ============================================================

print("=" * 70)
print("SHA-256 LINEAR OPERATOR ANALYSIS OVER GF(2)")
print("=" * 70)

for name, M in [("Σ₁ (Sigma1)", Sigma1), ("Σ₀ (Sigma0)", Sigma0),
                ("σ₀ (sigma0)", sigma0), ("σ₁ (sigma1)", sigma1)]:
    print(f"\n{'='*50}")
    print(f"Operator: {name}")
    print(f"{'='*50}")

    rank = gf2_rank(M)
    print(f"Rank: {rank}/32")

    if rank < 32:
        print(f"SINGULAR! Kernel dimension: {32 - rank}")
    else:
        print("Full rank (invertible)")

    # Characteristic polynomial
    charpoly = gf2_charpoly(M)
    print(f"Char. poly degree: {len(charpoly) - 1}")
    print(f"Char. poly: {poly_to_str(charpoly)}")

    # Check irreducibility
    irr = is_irreducible_gf2(charpoly)
    print(f"Irreducible over GF(2): {irr}")

    if not irr:
        print("REDUCIBLE! Factoring...")
        factors = factorize_gf2(charpoly)
        print(f"Number of factors: {len(factors)}")
        for i, (f, mult) in enumerate(factors):
            deg = len(f) - 1
            print(f"  Factor {i+1}: degree {deg}, mult {mult}: {poly_to_str(f)}")

        min_deg = min(len(f) - 1 for f, _ in factors)
        print(f"\n  MINIMUM FACTOR DEGREE: {min_deg}")
        if min_deg < 32:
            print(f"  *** INVARIANT SUBSPACE of dimension {min_deg} EXISTS ***")
            print(f"  *** Potential speedup: 2^({32-min_deg}) per barrier ***")
    else:
        print("NO invariant subspaces (operator is maximally mixing)")

# Also check products
print(f"\n{'='*50}")
print("PRODUCT Σ₀ · Σ₁")
print(f"{'='*50}")
product = Sigma0.dot(Sigma1) % 2
rank_p = gf2_rank(product)
print(f"Rank: {rank_p}/32")
charpoly_p = gf2_charpoly(product)
irr_p = is_irreducible_gf2(charpoly_p)
print(f"Char. poly irreducible: {irr_p}")
if not irr_p:
    print("REDUCIBLE! Factoring...")
    factors_p = factorize_gf2(charpoly_p)
    for i, (f, mult) in enumerate(factors_p):
        deg = len(f) - 1
        print(f"  Factor {i+1}: degree {deg}: {poly_to_str(f)}")
    min_deg_p = min(len(f) - 1 for f, _ in factors_p)
    print(f"  MINIMUM FACTOR DEGREE: {min_deg_p}")

# Check the full algebra generated by all 4 operators
print(f"\n{'='*50}")
print("ALGEBRA generated by {Σ₀, Σ₁, σ₀, σ₁}")
print(f"{'='*50}")

# Generate algebra by repeatedly multiplying and adding
basis = []
operators = [np.eye(N, dtype=int), Sigma0, Sigma1, sigma0, sigma1]
queue = list(operators)
seen = set()

def mat_to_key(M):
    return tuple(M.flatten() % 2)

for m in queue:
    seen.add(mat_to_key(m))

while queue and len(basis) < N*N:
    current = queue.pop(0)
    # Check if current is linearly independent of basis
    if basis:
        test_mat = np.vstack([b.flatten() for b in basis] + [current.flatten()]) % 2
        r = gf2_rank(test_mat)
        if r <= len(basis):
            continue  # linearly dependent
    basis.append(current)

    # Generate new elements
    for gen in [Sigma0, Sigma1, sigma0, sigma1]:
        for new_m in [current.dot(gen) % 2, gen.dot(current) % 2, (current + gen) % 2]:
            key = mat_to_key(new_m)
            if key not in seen:
                seen.add(key)
                queue.append(new_m)

    if len(basis) % 100 == 0:
        print(f"  Algebra dimension so far: {len(basis)}")

print(f"Algebra dimension: {len(basis)}")
print(f"Full End(GF(2)^32) dimension: {N*N} = 1024")
if len(basis) == N*N:
    print("ALGEBRA = End(GF(2)^32) — FULL! No invariant subspaces for joint action.")
else:
    print(f"ALGEBRA < End — dimension {len(basis)} < 1024")
    print(f"*** JOINT INVARIANT SUBSPACES MAY EXIST ***")
