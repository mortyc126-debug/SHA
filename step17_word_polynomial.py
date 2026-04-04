"""
Step 17: Build and solve the WORD-LEVEL polynomial

Step 16 showed: SHA as function Z/2^n → Z/2^n has degree n (not 2^n).
For mini-SHA: degree 4 in each of 4 word variables over Z/16.

Plan:
1. Build the EXPLICIT word-level polynomial using interpolation
2. Solve the collision system (2 word equations in 4 word unknowns)
3. Compare cost with birthday (2^4 = 16)

A degree-d polynomial in k variables over Z/m:
  F(x) = Σ_{|α|≤d} c_α · x^α  mod m
  where α = (α₁,...,α_k) multiindex, x^α = x₁^α₁ · ... · x_k^αk

Number of monomials: C(k+d, d) = C(4+4, 4) = 70.
So F is determined by 70 coefficients per output word.
With 16^4 = 65536 data points and 70 unknowns → MASSIVELY overdetermined.
"""

import numpy as np
from itertools import product as iterproduct
from step0_exact_algebra import mini_sha, N, MASK

N_MSG = 4
N_WORDS = 1 << N  # 16
MODULUS = N_WORDS   # Z/16


def enumerate_monomials(n_vars, max_degree):
    """
    Enumerate all monomials x₁^α₁ · ... · x_k^αk with Σαᵢ ≤ max_degree.
    Returns list of tuples (α₁, ..., α_k).
    """
    monomials = []

    def recurse(var_idx, remaining_degree, current):
        if var_idx == n_vars:
            monomials.append(tuple(current))
            return
        for d in range(remaining_degree + 1):
            current.append(d)
            recurse(var_idx + 1, remaining_degree - d, current)
            current.pop()

    recurse(0, max_degree, [])
    return monomials


def eval_monomial(alpha, x, modulus):
    """Evaluate monomial x^α mod modulus."""
    result = 1
    for i, a in enumerate(alpha):
        result = (result * pow(int(x[i]), int(a), modulus)) % modulus
    return result


def build_word_polynomial(M_base, R, out_idx):
    """
    Build the word-level polynomial for output word out_idx (0=δa, 1=δe).

    F(δW₀, δW₁, δW₂, δW₃) = output_word(M ⊕ δW) ⊕ output_word(M)

    Using XOR differences (δW applied as XOR to message).
    Polynomial over Z/2^n using ADDITIVE differences for the polynomial structure.
    """
    a_base, e_base = mini_sha(M_base, R)
    base_vals = [a_base, e_base]

    # Collect data points: for each δW, compute output difference
    data = {}  # δW tuple → output value
    for dw0 in range(N_WORDS):
        for dw1 in range(N_WORDS):
            for dw2 in range(N_WORDS):
                for dw3 in range(N_WORDS):
                    dw = (dw0, dw1, dw2, dw3)
                    M2 = [M_base[0]^dw0, M_base[1]^dw1, M_base[2]^dw2, M_base[3]^dw3]
                    a2, e2 = mini_sha(M2, R)
                    out = [a2 ^ a_base, e2 ^ e_base]
                    data[dw] = out[out_idx]

    return data


def interpolate_polynomial(data, n_vars, max_degree, modulus):
    """
    Interpolate a polynomial from data points.

    data: dict mapping (x₁,...,x_k) → y, all in Z/modulus
    Returns: coefficients for each monomial.

    Uses least squares over Z/modulus (try all monomial subsets).
    For exact interpolation: solve the Vandermonde system.
    """
    monomials = enumerate_monomials(n_vars, max_degree)
    n_monos = len(monomials)

    # Build Vandermonde matrix: V[i][j] = monomial_j evaluated at point_i
    points = sorted(data.keys())
    n_points = len(points)

    # Use a subset of points for fitting (n_monos points suffice)
    # But we'll use all and check consistency

    # For modular arithmetic, we need modular least squares
    # Simpler approach: use n_monos points, solve exactly

    V = np.zeros((min(n_points, n_monos * 2), n_monos), dtype=np.int64)
    y = np.zeros(min(n_points, n_monos * 2), dtype=np.int64)

    for i in range(min(n_points, n_monos * 2)):
        pt = points[i]
        y[i] = data[pt]
        for j, alpha in enumerate(monomials):
            V[i][j] = eval_monomial(alpha, pt, modulus)

    return V, y, monomials


def solve_modular_system(V, y, modulus):
    """
    Solve V·c ≡ y (mod modulus) for coefficient vector c.
    Uses Gaussian elimination mod modulus.
    Returns coefficients or None if inconsistent.
    """
    n_rows, n_cols = V.shape
    # Augmented matrix
    aug = np.zeros((n_rows, n_cols + 1), dtype=np.int64)
    aug[:, :n_cols] = V % modulus
    aug[:, n_cols] = y % modulus

    rank = 0
    pivot_cols = []
    for col in range(n_cols):
        # Find pivot with invertible leading element
        pivot = None
        for row in range(rank, n_rows):
            if aug[row][col] % modulus != 0:
                # Check if invertible mod modulus
                # For modulus = 2^k: invertible iff odd
                if aug[row][col] % 2 == 1:  # odd = invertible mod 2^k
                    pivot = row
                    break
        if pivot is None:
            # Try non-invertible but nonzero
            for row in range(rank, n_rows):
                if aug[row][col] % modulus != 0:
                    pivot = row
                    break
        if pivot is None:
            continue

        # Swap
        aug[[rank, pivot]] = aug[[pivot, rank]]
        pivot_cols.append(col)

        # Make pivot = 1 if possible
        pv = int(aug[rank][col]) % modulus
        # Find inverse mod modulus (if exists)
        inv = None
        for t in range(modulus):
            if (pv * t) % modulus == 1:
                inv = t
                break

        if inv is not None:
            aug[rank] = (aug[rank] * inv) % modulus

        # Eliminate
        for row in range(n_rows):
            if row != rank and aug[row][col] % modulus != 0:
                factor = int(aug[row][col])
                if inv is not None:
                    aug[row] = (aug[row] - factor * aug[rank]) % modulus
                else:
                    # Can't fully eliminate without inverse
                    pass

        rank += 1

    # Extract solution
    coeffs = np.zeros(n_cols, dtype=np.int64)
    for i, col in enumerate(pivot_cols):
        coeffs[col] = aug[i][n_cols] % modulus

    # Verify
    residual = (V @ coeffs - y) % modulus
    n_exact = np.sum(residual == 0)

    return coeffs, rank, n_exact


def main():
    import random
    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(N_MSG)]

    for R in [4, 6, 8]:
        print(f"\n{'='*80}")
        print(f"WORD-LEVEL POLYNOMIAL — R={R}")
        print(f"M = {[hex(w) for w in M]}")
        print(f"{'='*80}")

        for max_deg in [2, 3, 4]:
            monomials = enumerate_monomials(N_MSG, max_deg)
            print(f"\n  Degree ≤ {max_deg}: {len(monomials)} monomials in {N_MSG} variables")

            for out_idx, out_name in [(0, "δa"), (1, "δe")]:
                data = build_word_polynomial(M, R, out_idx)

                V, y, monos = interpolate_polynomial(data, N_MSG, max_deg, MODULUS)
                coeffs, rank, n_exact = solve_modular_system(V, y, MODULUS)

                # Test on ALL points
                n_correct = 0
                n_total = 0
                for pt in sorted(data.keys()):
                    pred = 0
                    for j, alpha in enumerate(monos):
                        pred = (pred + int(coeffs[j]) * eval_monomial(alpha, pt, MODULUS)) % MODULUS
                    if pred == data[pt]:
                        n_correct += 1
                    n_total += 1

                pct = n_correct / n_total * 100
                print(f"    {out_name}: rank={rank}, fit={n_correct}/{n_total} ({pct:.1f}%)"
                      f"{'  ★ EXACT FIT' if pct > 99.9 else ''}")

                # If exact fit: we have the polynomial!
                if pct > 99.9:
                    # Count nonzero coefficients
                    nonzero = np.sum(coeffs % MODULUS != 0)
                    print(f"         Nonzero coefficients: {nonzero}/{len(monos)}")

                    # Show the polynomial (if small)
                    if nonzero <= 20:
                        terms = []
                        for j, alpha in enumerate(monos):
                            c = int(coeffs[j]) % MODULUS
                            if c != 0:
                                vars_str = "·".join(
                                    f"x{i}^{a}" if a > 1 else f"x{i}"
                                    for i, a in enumerate(alpha) if a > 0
                                )
                                if not vars_str:
                                    vars_str = "1"
                                terms.append(f"{c}·{vars_str}")
                        print(f"         F = {' + '.join(terms)} (mod {MODULUS})")

        # ================================================================
        # SOLVE the word-level system for collision
        # ================================================================
        print(f"\n  COLLISION SOLVING at word level:")

        # Collision: F_a(δW) = 0 AND F_e(δW) = 0
        # If we have the explicit polynomial, we can:
        # 1. Enumerate solutions (brute force on word variables: 16^4 = 65536)
        # 2. Use word-level algebraic methods

        data_a = build_word_polynomial(M, R, 0)
        data_e = build_word_polynomial(M, R, 1)

        # Count collisions
        collisions = []
        for dw in sorted(data_a.keys()):
            if data_a[dw] == 0 and data_e[dw] == 0 and dw != (0,0,0,0):
                collisions.append(dw)

        print(f"  Collisions: {len(collisions)}")
        print(f"  Birthday cost: {N_WORDS} = {N_WORDS}")
        print(f"  Word-brute cost: {N_WORDS**N_MSG} = {N_WORDS**N_MSG}")

        # Key: with degree-d polynomial over Z/m in k vars:
        # Number of solutions ≈ m^(k-2) for 2 equations (generically)
        expected = MODULUS ** (N_MSG - 2)
        print(f"  Expected (generic degree-{N} system): {expected}")
        print(f"  Actual: {len(collisions)}")
        print(f"  Ratio: {len(collisions)/max(expected,1):.3f}")


if __name__ == "__main__":
    main()
