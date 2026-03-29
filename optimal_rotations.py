"""
Search for rotation constants giving IRREDUCIBLE characteristic polynomial.
If found: SHA-256 with these constants would have NO fixed vectors = strictly stronger.
If none exist: all choices are equally "weak" = SHA-256 not at fault.

Also: exploit the 2D invariant subspace of Σ₀·Σ₁ product.
"""
import numpy as np
from itertools import combinations

N = 32

def rotr_matrix(k, n=N):
    M = np.zeros((n, n), dtype=int)
    for i in range(n):
        M[i][(i + k) % n] = 1
    return M

def shr_matrix(k, n=N):
    M = np.zeros((n, n), dtype=int)
    for j in range(n):
        if j + k < n:
            M[j][j + k] = 1
    return M

def xor_matrices(*matrices):
    result = matrices[0].copy()
    for m in matrices[1:]:
        result = (result + m) % 2
    return result

def gf2_rank(M):
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
    """Characteristic polynomial via Krylov method"""
    n = M.shape[0]
    import random
    random.seed(42)
    v = np.array([random.randint(0,1) for _ in range(n)], dtype=int)
    krylov = [v.copy()]
    current = v.copy()
    for i in range(n):
        current = M.dot(current) % 2
        krylov.append(current.copy())
    A = np.column_stack(krylov[:n]) % 2
    b = krylov[n] % 2
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
    x = np.zeros(n, dtype=int)
    for col, row in pivots.items():
        x[col] = aug[row][-1]
    return list(x) + [1]

def gf2_poly_mul(a, b):
    if not a or not b:
        return [0]
    result = [0] * (len(a) + len(b) - 1)
    for i, ca in enumerate(a):
        if ca:
            for j, cb in enumerate(b):
                if cb:
                    result[i + j] ^= 1
    while len(result) > 1 and result[-1] == 0:
        result = result[:-1]
    return result

def gf2_poly_mod(a, m):
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
    result = [1]
    base = gf2_poly_mod(base, mod)
    while exp > 0:
        if exp & 1:
            result = gf2_poly_mod(gf2_poly_mul(result, base), mod)
        base = gf2_poly_mod(gf2_poly_mul(base, base), mod)
        exp >>= 1
    return result

def is_irreducible_gf2(poly):
    n = len(poly) - 1
    if n <= 0:
        return False
    if poly[0] == 0:
        return False
    x_poly = [0, 1]
    x_pow = gf2_poly_powmod(x_poly, 2**n, poly)
    diff = list(x_pow)
    if len(diff) <= 1:
        diff = diff + [0] * (2 - len(diff))
    if len(diff) > 1:
        diff[1] ^= 1
    else:
        diff.append(1)
    while len(diff) > 1 and diff[-1] == 0:
        diff = diff[:-1]
    if any(d == 1 for d in diff):
        return False
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
        if len(g) > 1:
            return False
    return True

def gf2_poly_gcd(a, b):
    a = list(a)
    b = list(b)
    while any(x == 1 for x in b):
        while len(a) >= len(b) and any(x == 1 for x in a):
            if a[-1] == 0:
                a = a[:-1]
                continue
            shift = len(a) - len(b)
            for i in range(len(b)):
                a[i + shift] ^= b[i]
            while len(a) > 0 and a[-1] == 0:
                a = a[:-1]
            if not a:
                a = [0]
        a, b = b, a
    return a if any(x == 1 for x in a) else [0]

def fixed_space_dim(M):
    """Dimension of fixed space ker(M+I) over GF(2)"""
    MpI = (M + np.eye(N, dtype=int)) % 2
    return N - gf2_rank(MpI)

def min_factor_degree(charpoly):
    """Find minimum factor degree of charpoly over GF(2)"""
    n = len(charpoly) - 1
    if n <= 1:
        return n
    # Check if (1+x) divides (i.e., poly evaluated at 1 = 0 mod 2)
    if sum(charpoly) % 2 == 0:
        return 1
    # Check small degree factors
    x_poly = [0, 1]
    for k in range(1, n // 2 + 1):
        x_pow = gf2_poly_powmod(x_poly, 2**k, charpoly)
        diff = list(x_pow)
        if len(diff) <= 1:
            diff = diff + [0] * (2 - len(diff))
        if len(diff) > 1:
            diff[1] ^= 1
        else:
            diff.append(1)
        while len(diff) > 1 and diff[-1] == 0:
            diff = diff[:-1]
        g = gf2_poly_gcd(list(diff), list(charpoly))
        while len(g) > 1 and g[-1] == 0:
            g = g[:-1]
        if len(g) > 1:
            return min(k, len(g) - 1)
    return n

# ============================================================
# PART 1: Search for ROTR_a + ROTR_b + SHR_c with irreducible charpoly
# ============================================================
print("=" * 70)
print("SEARCH: ROTR_a ⊕ ROTR_b ⊕ SHR_c with irreducible χ over GF(2)^32")
print("=" * 70)

best_min_deg = 0
best_params = None
irreducible_count = 0
total_count = 0

results_by_min_deg = {}

for a in range(1, 31):
    for b in range(a + 1, 32):
        for c in range(1, 31):
            if c == a or c == b:
                continue
            M = xor_matrices(rotr_matrix(a), rotr_matrix(b), shr_matrix(c))
            r = gf2_rank(M)
            if r < N:
                continue  # singular, skip

            cp = gf2_charpoly(M)
            total_count += 1

            if is_irreducible_gf2(cp):
                irreducible_count += 1
                if irreducible_count <= 5:
                    print(f"  IRREDUCIBLE: a={a}, b={b}, c={c}")
                min_d = N
            else:
                min_d = min_factor_degree(cp)

            if min_d not in results_by_min_deg:
                results_by_min_deg[min_d] = 0
            results_by_min_deg[min_d] += 1

            if min_d > best_min_deg:
                best_min_deg = min_d
                best_params = (a, b, c)

            if total_count % 2000 == 0:
                print(f"  Tested {total_count}... best min_deg so far: {best_min_deg}, irreducible: {irreducible_count}")

print(f"\nTotal tested: {total_count}")
print(f"Irreducible: {irreducible_count} ({100*irreducible_count/total_count:.1f}%)")
print(f"Best min factor degree: {best_min_deg} with params {best_params}")
print(f"\nDistribution of min factor degree:")
for d in sorted(results_by_min_deg.keys()):
    print(f"  min_deg = {d}: {results_by_min_deg[d]} ({100*results_by_min_deg[d]/total_count:.1f}%)")

# SHA-256 Sigma1 for comparison
print(f"\nSHA-256 Σ₁ (6,11,25): min factor degree = {min_factor_degree(gf2_charpoly(xor_matrices(rotr_matrix(6), rotr_matrix(11), shr_matrix(25))))}")

# ============================================================
# PART 2: Pure ROTR_a + ROTR_b + ROTR_c (no SHR) — like Sigma0
# ============================================================
print("\n" + "=" * 70)
print("SEARCH: ROTR_a ⊕ ROTR_b ⊕ ROTR_c (pure rotation, like Σ₀)")
print("=" * 70)

pure_irr = 0
pure_total = 0
pure_best_min = 0
pure_best_params = None

for a in range(1, 30):
    for b in range(a + 1, 31):
        for c_val in range(b + 1, 32):
            M = xor_matrices(rotr_matrix(a), rotr_matrix(b), rotr_matrix(c_val))
            r = gf2_rank(M)
            if r < N:
                continue
            cp = gf2_charpoly(M)
            pure_total += 1
            if is_irreducible_gf2(cp):
                pure_irr += 1
                min_d = N
            else:
                min_d = min_factor_degree(cp)
            if min_d > pure_best_min:
                pure_best_min = min_d
                pure_best_params = (a, b, c_val)

print(f"Total: {pure_total}, Irreducible: {pure_irr}")
print(f"Best min factor: {pure_best_min} with {pure_best_params}")
print(f"SHA-256 Σ₀ (2,13,22): charpoly = (1+x)^32, min factor = 1")

# ============================================================
# PART 3: Invariant subspace of Σ₀·Σ₁ (degree 2 factor)
# ============================================================
print("\n" + "=" * 70)
print("INVARIANT SUBSPACE of Σ₀·Σ₁ (min factor degree 2)")
print("=" * 70)

Sigma0 = xor_matrices(rotr_matrix(2), rotr_matrix(13), rotr_matrix(22))
Sigma1 = xor_matrices(rotr_matrix(6), rotr_matrix(11), shr_matrix(25))
product = Sigma0.dot(Sigma1) % 2

# Factor x^2 + x + 1: eigenvalue is primitive cube root of unity in GF(4)
# Ker of product^2 + product + I (evaluation at roots of x^2+x+1)
test = (product.dot(product) + product + np.eye(N, dtype=int)) % 2
ker_dim = N - gf2_rank(test)
print(f"ker(M² + M + I) dimension: {ker_dim}")

if ker_dim > 0:
    # Find basis of this kernel
    m = test.copy() % 2
    pivot_cols = []
    r = 0
    for col in range(N):
        pivot = None
        for row in range(r, N):
            if m[row][col] == 1:
                pivot = row
                break
        if pivot is None:
            continue
        m[[r, pivot]] = m[[pivot, r]]
        for row in range(N):
            if row != r and m[row][col] == 1:
                m[row] = (m[row] + m[r]) % 2
        pivot_cols.append(col)
        r += 1

    free_cols = [c for c in range(N) if c not in pivot_cols]
    print(f"Free variables: {free_cols[:5]}...")

    for fc in free_cols[:3]:
        v = np.zeros(N, dtype=int)
        v[fc] = 1
        for i in range(len(pivot_cols)-1, -1, -1):
            pc = pivot_cols[i]
            v[pc] = m[i][fc]
        val = sum(v[j] * (1 << j) for j in range(N))
        # Check: (Σ₀·Σ₁)² + (Σ₀·Σ₁) + I annihilates v
        check = test.dot(v) % 2
        print(f"  Invariant vector: 0x{val:08x}, HW={sum(v)}, check={'OK' if sum(check)==0 else 'FAIL'}")

        # Check if this vector is also fixed by Sigma0 or Sigma1 individually
        S0v = Sigma0.dot(v) % 2
        S1v = Sigma1.dot(v) % 2
        print(f"    Σ₀(v)=v: {np.array_equal(S0v, v)}, Σ₁(v)=v: {np.array_equal(S1v, v)}")
