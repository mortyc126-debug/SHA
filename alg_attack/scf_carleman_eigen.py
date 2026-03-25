#!/usr/bin/env python3
"""
SCF: Собственные значения линейной части Carleman (матрица A).

Факт: rank(Carleman) = state_dim для ЛЮБОГО R.
Факт: A (линейная часть) имеет полный rank.
Вопрос: какова структура A? Собственные значения, Jordan-форма?

A — это линейная часть одного XOR-SHA раунда над GF(2).
A^R — линейная часть R раундов.
Если A имеет малый минимальный многочлен → A^64 вычислима эффективно.
Если A имеет цикл порядка d → A^d = I → A^64 = A^{64 mod d}.

Для GF(2): собственные значения ∈ GF(2) = {0, 1}.
Rank(A) = dim → det(A) = 1 → A обратима → все собственные значения = 1.
Но Jordan-форма может быть нетривиальной: (A-I)^k = 0 для некоторого k.
"""
import os, sys

def rotr_n(x, r, n):
    return ((x >> r) | (x << (n - r))) & ((1 << n) - 1)

class MiniSHA:
    def __init__(self, n=4):
        self.n = n
        self.mask = (1 << n) - 1
        self.sig0_rots = [max(1, n//8), max(2, n//3), n//2]
        self.sig1_rots = [max(1, n//5), max(2, n//3), max(3, 3*n//4)]

    def Sig0(self, x):
        n = self.n
        return rotr_n(x, self.sig0_rots[0], n) ^ rotr_n(x, self.sig0_rots[1], n) ^ rotr_n(x, self.sig0_rots[2], n)

    def Sig1(self, x):
        n = self.n
        return rotr_n(x, self.sig1_rots[0], n) ^ rotr_n(x, self.sig1_rots[1], n) ^ rotr_n(x, self.sig1_rots[2], n)

    def Ch(self, e, f, g):
        return (e & f) ^ (~e & g) & self.mask

    def Maj(self, a, b, c):
        return (a & b) ^ (a & c) ^ (b & c)

    def xor_round(self, state, W=0, K=0):
        a,b,c,d,e,f,g,h = state
        T1 = h ^ self.Sig1(e) ^ self.Ch(e,f,g) ^ K ^ W
        T2 = self.Sig0(a) ^ self.Maj(a,b,c)
        return [T1^T2, a, b, c, d^T1, e, f, g]

    def state_to_bits(self, state):
        bits = []
        for w in state:
            for i in range(self.n):
                bits.append((w >> i) & 1)
        return bits

    def bits_to_state(self, bits):
        state = []
        for r in range(8):
            val = 0
            for i in range(self.n):
                val |= bits[r * self.n + i] << i
            state.append(val)
        return state


def build_linear_matrix(sha):
    """Build matrix A of the linear part of one XOR-SHA round."""
    n = sha.n
    dim = 8 * n
    f_zero = sha.state_to_bits(sha.xor_round([0]*8))

    A = []
    for out_bit in range(dim):
        row = []
        for in_bit in range(dim):
            test = [0] * dim
            test[in_bit] = 1
            state = sha.bits_to_state(test)
            f_ei = sha.state_to_bits(sha.xor_round(state))
            row.append(f_ei[out_bit] ^ f_zero[out_bit])
        A.append(row)
    return A


def mat_mul_gf2(A, B):
    """Multiply two GF(2) matrices."""
    n = len(A)
    m = len(B[0])
    k = len(B)
    C = [[0]*m for _ in range(n)]
    for i in range(n):
        for j in range(m):
            val = 0
            for l in range(k):
                val ^= A[i][l] & B[l][j]
            C[i][j] = val
    return C


def mat_pow_gf2(A, p):
    """Matrix power A^p over GF(2)."""
    n = len(A)
    result = [[1 if i==j else 0 for j in range(n)] for i in range(n)]  # Identity
    base = [list(row) for row in A]
    while p > 0:
        if p & 1:
            result = mat_mul_gf2(result, base)
        base = mat_mul_gf2(base, base)
        p >>= 1
    return result


def mat_add_gf2(A, B):
    return [[A[i][j] ^ B[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def mat_is_zero(A):
    return all(all(x == 0 for x in row) for row in A)


def mat_is_identity(A):
    n = len(A)
    for i in range(n):
        for j in range(n):
            if A[i][j] != (1 if i==j else 0):
                return False
    return True


def gf2_rank(mat):
    m = [list(r) for r in mat]
    nr, nc = len(m), len(m[0]) if m else 0
    rank = 0
    for col in range(nc):
        pv = -1
        for r in range(rank, nr):
            if m[r][col]: pv = r; break
        if pv == -1: continue
        m[rank], m[pv] = m[pv], m[rank]
        for r in range(nr):
            if r != rank and m[r][col]:
                for c in range(nc): m[r][c] ^= m[rank][c]
        rank += 1
    return rank


# ============================================================
# EXP 1: Order of A — find minimal d such that A^d = I
# ============================================================
def exp1_order(A, max_order=1000):
    n = len(A)
    I = [[1 if i==j else 0 for j in range(n)] for i in range(n)]

    power = [list(row) for row in A]
    for d in range(1, max_order+1):
        if mat_is_identity(power):
            return d
        power = mat_mul_gf2(power, A)
    return None


# ============================================================
# EXP 2: Nilpotency index of (A - I) — Jordan structure
# ============================================================
def exp2_nilpotency(A):
    n = len(A)
    I = [[1 if i==j else 0 for j in range(n)] for i in range(n)]
    N = mat_add_gf2(A, I)  # A - I = A + I over GF(2)

    power = [list(row) for row in N]
    for k in range(1, n+1):
        r = gf2_rank(power)
        if r == 0:
            return k, True  # (A-I)^k = 0 → nilpotent!
        power = mat_mul_gf2(power, N)

    return n, False  # Not nilpotent within n steps


# ============================================================
# EXP 3: Minimal polynomial of A
# ============================================================
def exp3_minimal_poly(A):
    """Find minimal polynomial of A over GF(2).
    Test: A^k + c_{k-1}A^{k-1} + ... + c_0 I = 0
    Using trial: check if A^k is in span of {I, A, ..., A^{k-1}}."""
    n = len(A)

    powers = []
    # Flatten matrices to vectors for rank computation
    I = [[1 if i==j else 0 for j in range(n)] for i in range(n)]
    powers_flat = []

    current = [list(row) for row in I]
    for k in range(n + 1):
        flat = []
        for row in current:
            flat.extend(row)
        powers_flat.append(flat)

        # Check if current power is in span of previous
        if k > 0:
            test_matrix = powers_flat[:k+1]  # rows = flattened I, A, ..., A^k
            r = gf2_rank(test_matrix)
            if r < k + 1:
                return k  # A^k is linearly dependent on {I, A, ..., A^{k-1}}

        current = mat_mul_gf2(current, A)

    return n  # Minimal polynomial has degree n (maximum)


# ============================================================
# EXP 4: What does A^64 look like?
# ============================================================
def exp4_power_structure(A, sha):
    """Compute A^R for various R, compare with actual multi-round SHA."""
    n_bits = len(A)

    print(f"\n  A^R structure:")
    print(f"  {'R':>4} | {'rank(A^R)':>10} | {'A^R=I?':>7} | {'rank(A^R-I)':>12}")
    print("  " + "-"*50)

    for R in [1, 2, 4, 8, 16, 32, 64, 128]:
        AR = mat_pow_gf2(A, R)
        rank_AR = gf2_rank(AR)
        is_I = mat_is_identity(AR)
        n_dim = len(A)
        I_mat = [[1 if i==j else 0 for j in range(n_dim)] for i in range(n_dim)]
        AR_minus_I = mat_add_gf2(AR, I_mat)
        rank_diff = gf2_rank(AR_minus_I)

        marker = " ★ IDENTITY!" if is_I else ""
        print(f"  {R:4d} | {rank_AR:10d} | {'YES' if is_I else 'no':>7} | {rank_diff:12d}{marker}")

        if is_I:
            return R  # A has order R

    return None


# ============================================================
if __name__ == '__main__':
    print("="*70)
    print("SCF: EIGENSTRUCTURE OF CARLEMAN LINEAR PART")
    print("="*70)

    for n in [4, 6, 8]:
        sha = MiniSHA(n)
        A = build_linear_matrix(sha)
        dim = 8 * n

        print(f"\n{'='*60}")
        print(f"n = {n} (state dim = {dim})")
        print(f"{'='*60}")

        # Rank
        rank = gf2_rank(A)
        print(f"  rank(A) = {rank}/{dim}")

        # Order
        print(f"\n  Searching for order of A (A^d = I)...")
        order = exp1_order(A, max_order=500)
        if order:
            print(f"  ★ Order of A = {order}")
            print(f"    A^{order} = I (identity)")
            print(f"    A^64 = A^{64 % order}")
        else:
            print(f"  Order > 500 (not found)")

        # Nilpotency
        print(f"\n  Jordan structure: (A-I) nilpotency index...")
        nil_idx, is_nil = exp2_nilpotency(A)
        if is_nil:
            print(f"  (A-I)^{nil_idx} = 0 → A is UNIPOTENT")
            print(f"    All eigenvalues = 1")
            print(f"    Jordan blocks of size ≤ {nil_idx}")
        else:
            print(f"  (A-I) not nilpotent after {nil_idx} steps")

        # Minimal polynomial
        print(f"\n  Minimal polynomial degree...")
        min_deg = exp3_minimal_poly(A)
        print(f"  deg(min_poly) = {min_deg}")
        print(f"    (vs matrix dim = {dim})")
        if min_deg < dim:
            print(f"    ★ Cayley-Hamilton compression: {dim}/{min_deg} = {dim/min_deg:.1f}x")

        # Power structure
        exp4_power_structure(A, sha)

    print(f"\n{'='*70}")
    print("ИТОГ")
    print(f"{'='*70}")
    print("""
  Если A имеет малый порядок d:
    A^64 = A^{64 mod d} → можно свести 64 раунда к (64 mod d)
    Для коллизии: нужно решить f(x) = f(y) для A^{64 mod d}

  Если A унипотентна ((A-I)^k = 0):
    A^R = I + R(A-I) + C(R,2)(A-I)^2 + ... + C(R,k-1)(A-I)^{k-1}
    Это ПОЛИНОМИАЛЬНАЯ формула для A^R!
    Вместо возведения в степень — полином степени k-1 от R.

  Если min_poly маленький:
    A удовлетворяет p(A) = 0 для полинома степени min_deg.
    A^R можно вычислить как полином степени < min_deg от A.
    """)
