"""Find fixed vectors of SHA-256 linear operators"""
import numpy as np

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

Sigma1 = xor_matrices(rotr_matrix(6), rotr_matrix(11), shr_matrix(25))
Sigma0 = xor_matrices(rotr_matrix(2), rotr_matrix(13), rotr_matrix(22))
sigma0 = xor_matrices(rotr_matrix(7), rotr_matrix(18), shr_matrix(3))
sigma1 = xor_matrices(rotr_matrix(17), rotr_matrix(19), shr_matrix(10))

def find_fixed_space(M, name):
    """Find kernel of (M + I) over GF(2) = fixed vectors of M"""
    MpI = (M + np.eye(N, dtype=int)) % 2

    # Gaussian elimination to find kernel
    m = MpI.copy()
    pivot_cols = []
    row_idx = 0
    for col in range(N):
        pivot = None
        for r in range(row_idx, N):
            if m[r][col] == 1:
                pivot = r
                break
        if pivot is None:
            continue
        m[[row_idx, pivot]] = m[[pivot, row_idx]]
        for r in range(N):
            if r != row_idx and m[r][col] == 1:
                m[r] = (m[r] + m[row_idx]) % 2
        pivot_cols.append(col)
        row_idx += 1

    rank = len(pivot_cols)
    kernel_dim = N - rank

    print(f"\n{name}:")
    print(f"  rank(M+I) = {rank}")
    print(f"  dim(ker(M+I)) = {kernel_dim} = dim(fixed space)")

    if kernel_dim > 0:
        # Find kernel basis
        free_cols = [c for c in range(N) if c not in pivot_cols]

        for fc in free_cols[:5]:  # show up to 5 fixed vectors
            v = np.zeros(N, dtype=int)
            v[fc] = 1
            # Back-substitute
            for i in range(rank-1, -1, -1):
                pc = pivot_cols[i]
                v[pc] = m[i][fc]

            # Verify
            Mv = M.dot(v) % 2
            check = np.array_equal(Mv, v)

            # Convert to hex
            val = sum(v[i] * (1 << i) for i in range(N))
            val_Mv = sum(Mv[i] * (1 << i) for i in range(N))

            print(f"  Fixed vector: 0x{val:08x} (binary: {''.join(str(v[31-i]) for i in range(32))})")
            print(f"  M(v)        : 0x{val_Mv:08x} (check: {'OK' if check else 'FAIL'})")
            print(f"  HW(v) = {sum(v)}")

find_fixed_space(Sigma1, "Σ₁ (Sigma1)")
find_fixed_space(Sigma0, "Σ₀ (Sigma0)")
find_fixed_space(sigma0, "σ₀ (sigma0)")
find_fixed_space(sigma1, "σ₁ (sigma1)")

# Also find common fixed space of ALL four
print("\n" + "="*50)
print("COMMON fixed space of ALL four operators")
print("="*50)

# ker(Σ₁+I) ∩ ker(Σ₀+I) ∩ ker(σ₀+I) ∩ ker(σ₁+I)
# = ker of stacked matrix [(Σ₁+I); (Σ₀+I); (σ₀+I); (σ₁+I)]
stacked = np.vstack([
    (Sigma1 + np.eye(N, dtype=int)) % 2,
    (Sigma0 + np.eye(N, dtype=int)) % 2,
    (sigma0 + np.eye(N, dtype=int)) % 2,
    (sigma1 + np.eye(N, dtype=int)) % 2,
]) % 2

# Find rank of stacked
m = stacked.copy()
rows_m, cols_m = m.shape
pivot_cols = []
row_idx = 0
for col in range(cols_m):
    pivot = None
    for r in range(row_idx, rows_m):
        if m[r][col] == 1:
            pivot = r
            break
    if pivot is None:
        continue
    m[[row_idx, pivot]] = m[[pivot, row_idx]]
    for r in range(rows_m):
        if r != row_idx and m[r][col] == 1:
            m[r] = (m[r] + m[row_idx]) % 2
    pivot_cols.append(col)
    row_idx += 1

rank_s = len(pivot_cols)
common_dim = N - rank_s
print(f"Rank of stacked system: {rank_s}")
print(f"Dimension of COMMON fixed space: {common_dim}")

if common_dim > 0:
    free_cols = [c for c in range(N) if c not in pivot_cols]
    for fc in free_cols:
        v = np.zeros(N, dtype=int)
        v[fc] = 1
        for i in range(rank_s-1, -1, -1):
            pc = pivot_cols[i]
            v[pc] = m[i][fc]
        val = sum(v[i] * (1 << i) for i in range(N))
        print(f"  Common fixed vector: 0x{val:08x}, HW={sum(v)}")

        # Verify all four
        for name, M in [("Σ₁", Sigma1), ("Σ₀", Sigma0), ("σ₀", sigma0), ("σ₁", sigma1)]:
            Mv = M.dot(v) % 2
            check = np.array_equal(Mv, v)
            print(f"    {name}(v) = v: {check}")
elif common_dim == 0:
    print("NO common fixed vector — good for SHA-256 security")
    print("Individual fixed spaces exist but do NOT overlap")
