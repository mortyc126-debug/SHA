"""Check nilpotency degree of Sigma0 + I and invariant subspace structure"""
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

Sigma0 = xor_matrices(rotr_matrix(2), rotr_matrix(13), rotr_matrix(22))
Sigma1 = xor_matrices(rotr_matrix(6), rotr_matrix(11), shr_matrix(25))
I = np.eye(N, dtype=int)

# Nilpotency of Sigma0 + I
S0pI = (Sigma0 + I) % 2

print("Nilpotency degree of (Σ₀ + I):")
power = I.copy()
for k in range(1, 33):
    power = power.dot(S0pI) % 2
    rank = np.linalg.matrix_rank(power.astype(float))
    # Better: GF2 rank
    m = power.copy() % 2
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
        r += 1

    print(f"  (Σ₀+I)^{k}: rank = {r}")
    if r == 0:
        print(f"  *** NILPOTENT OF DEGREE {k} ***")
        break

# Invariant subspace chain (Jordan filtration)
print("\nJordan filtration of Σ₀:")
print("  V_k = ker((Σ₀+I)^k)")
for k in range(1, 33):
    # ker of (S0pI)^k
    power_k = I.copy()
    for _ in range(k):
        power_k = power_k.dot(S0pI) % 2

    m = power_k.copy() % 2
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
        r += 1

    dim_ker = N - r
    print(f"  dim V_{k} = {dim_ker}")
    if dim_ker == N:
        break

# Check: is Sigma1 fixed space contained in any Sigma0 Jordan subspace?
print("\nΣ₁ fixed vectors vs Σ₀ Jordan filtration:")

# Fixed vectors of Sigma1
S1pI = (Sigma1 + I) % 2
m_s1 = S1pI.copy() % 2
r_s1 = 0
pivots_s1 = []
for col in range(N):
    pivot = None
    for row in range(r_s1, N):
        if m_s1[row][col] == 1:
            pivot = row
            break
    if pivot is None:
        continue
    m_s1[[r_s1, pivot]] = m_s1[[pivot, r_s1]]
    for row in range(N):
        if row != r_s1 and m_s1[row][col] == 1:
            m_s1[row] = (m_s1[row] + m_s1[r_s1]) % 2
    pivots_s1.append(col)
    r_s1 += 1

free_cols = [c for c in range(N) if c not in pivots_s1]
fixed_vectors = []
for fc in free_cols:
    v = np.zeros(N, dtype=int)
    v[fc] = 1
    for i in range(r_s1-1, -1, -1):
        pc = pivots_s1[i]
        v[pc] = m_s1[i][fc]
    fixed_vectors.append(v)

for idx, v in enumerate(fixed_vectors):
    val = sum(v[i] * (1 << i) for i in range(N))
    print(f"\n  v_{idx+1} = 0x{val:08x}")

    for k in range(1, 33):
        power_k = I.copy()
        for _ in range(k):
            power_k = power_k.dot(S0pI) % 2
        result = power_k.dot(v) % 2
        if sum(result) == 0:
            print(f"    (Σ₀+I)^{k} · v = 0  →  v ∈ V_{k} (Jordan level {k})")
            break
