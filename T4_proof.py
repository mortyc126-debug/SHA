"""
ANALYTICAL PROOF OF T4: |{y : rank(J_{C_y}) = k}| = 2·C(n-1, k).

Key insight: the Jacobian of C_y at x=0 has a SPECIFIC STRUCTURE.

C_y(x)[k] = MAJ(x[k-1], y[k-1], C_y(x)[k-1])
           = x[k-1]·y[k-1] + x[k-1]·C[k-1] + y[k-1]·C[k-1]

At x=0: C_y(0) = 0 (all carries zero when first operand = 0).

Jacobian J[k][j] = ∂C_y(x)[k]/∂x[j] evaluated at x=0:

Recursion:
  J[0][j] = 0  (carry into bit 0 = always 0)
  J[k][j] = δ(j, k-1)·y[k-1] + y[k-1]·J[k-1][j]

Solving:
  J[k][j] = y[k-1] · y[k-2] · ... · y[j]  for j < k
  J[k][j] = 0                                for j ≥ k

This is a LOWER TRIANGULAR matrix with:
  - Zeros on and above diagonal
  - Sub-diagonal entry J[k][k-1] = y[k-1]
  - Below that: PRODUCTS of consecutive y-bits

RANK:
  Row k (k≥1) is nonzero iff J[k][k-1] = y[k-1] = 1.
  (If y[k-1]=0, then J[k][j]=0 for all j, since the product includes y[k-1].)

  rank(J) = |{k ∈ {1,...,n-1} : y[k-1] = 1}|
           = HW(y[0..n-2])

  This counts bits 0 through n-2 of y that are 1.
  Bit n-1 of y does NOT affect rank.

DISTRIBUTION:
  For each value of bits y[0..n-2]: rank = HW(y[0..n-2]) = k
  Number of such y[0..n-2] patterns: C(n-1, k)
  Bit y[n-1] can be 0 or 1: factor of 2.

  |{y ∈ {0,...,2^n-1} : rank(J) = k}| = 2 · C(n-1, k)

QED.

This also proves nilpotency of J:
  J is strictly lower triangular → J^n = 0.
  (Any strictly lower triangular n×n matrix is nilpotent with index ≤ n.)
"""

import random, math
from collections import Counter


def carry_jacobian_at_zero(y, n):
    """Compute Jacobian of C_y at x=0 analytically."""
    J = [[0]*n for _ in range(n)]
    for k in range(1, n):
        for j in range(k):
            # J[k][j] = product of y[i] for i = j to k-1
            prod = 1
            for i in range(j, k):
                prod &= (y >> i) & 1
            J[k][j] = prod
    return J


def gf2_rank(matrix, n):
    m = [list(row) for row in matrix]
    rank = 0
    for col in range(n):
        pivot = None
        for row in range(rank, n):
            if m[row][col] == 1:
                pivot = row; break
        if pivot is None: continue
        m[rank], m[pivot] = m[pivot], m[rank]
        for row in range(n):
            if row != rank and m[row][col] == 1:
                for c in range(n): m[row][c] ^= m[rank][c]
        rank += 1
    return rank


def verify_proof():
    print("=" * 80)
    print("PROOF OF T4: Carry Jacobian rank = 2·C(n-1, k)")
    print("=" * 80)

    for n in [4, 6, 8, 10]:
        N = 2**n

        # Method 1: Direct computation of rank
        ranks_direct = []
        for y_val in range(N):
            # Compute C_y Jacobian by flipping bits
            base = 0  # C_y(0) = 0
            J = []
            for j in range(n):
                x_flip = 1 << j
                # Compute C_y(x_flip)
                c = 0; cv = 0
                for k in range(n):
                    cv |= (c << k)
                    xk = (x_flip >> k) & 1
                    yk = (y_val >> k) & 1
                    c = (xk & yk) | (xk & c) | (yk & c)
                J.append([(cv >> k) & 1 for k in range(n)])
            ranks_direct.append(gf2_rank(J, n))

        # Method 2: Analytical formula (rank = HW(y[0..n-2]))
        ranks_formula = []
        for y_val in range(N):
            hw_low = bin(y_val & ((1 << (n-1)) - 1)).count('1')  # HW of bits 0..n-2
            ranks_formula.append(hw_low)

        # Compare
        match = (ranks_direct == ranks_formula)

        # Distribution
        dist_direct = Counter(ranks_direct)
        dist_formula = {k: 2 * math.comb(n-1, k) for k in range(n)}

        dist_match = all(dist_direct.get(k, 0) == dist_formula.get(k, 0) for k in range(n))

        print(f"\n  n={n}:")
        print(f"    Pointwise match (rank = HW(y[0..n-2])): {match}")
        print(f"    Distribution match (2·C(n-1,k)): {dist_match}")

        if n <= 8:
            print(f"    Direct:  {dict(sorted(dist_direct.items()))}")
            print(f"    Formula: {dict(sorted(dist_formula.items()))}")

    print(f"""
    ══════════════════════════════════════════════════════════════
    THEOREM T4 (Carry Jacobian Rank): ANALYTICALLY PROVED
    ══════════════════════════════════════════════════════════════

    THEOREM: For C_y: {{0,1}}^n → {{0,1}}^n (carry operator):
      rank(J_{{C_y}}|_{{x=0}}) = HW(y[0..n-2])

    COROLLARY: |{{y : rank = k}}| = 2 · C(n-1, k)

    PROOF:
      1. J[k][j] = ∏_{{i=j}}^{{k-1}} y[i]  (product of consecutive y-bits)
      2. Row k nonzero ⟺ y[k-1] = 1
      3. rank = count of k with y[k-1]=1 = HW(y[0..n-2])
      4. Bit n-1 doesn't affect rank → factor 2
      5. |{{y: HW(y[0..n-2])=k}}| = C(n-1,k) · 2  ∎

    BONUS: J is strictly lower triangular → J^n = 0.
    This gives ALTERNATIVE PROOF of T3 (nilpotency) at Jacobian level.
    """)


if __name__ == "__main__":
    verify_proof()
