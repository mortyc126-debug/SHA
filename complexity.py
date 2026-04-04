"""
BTE COMPLEXITY THEORY: What does the class structure imply about hardness?

From our theory:
  Layer rank = 2R - 1 (Theorem 1)
  Number of layers to full rank D = ceil(n_msg × n / (2R - 1))
  Each layer = degree-2 GF(2) system (from Ch/Maj)
  Layers connected by carry bridges (also degree 2 per step)

The "effective algebraic degree" of the full system:
  Each layer adds one level of degree-2 composition.
  D layers → degree 2^D.

For SHA-256: D = 4, effective degree = 2^4 = 16.

This is BETWEEN:
  GF(2) thinks: degree 32 (every bit position = degree 32 polynomial)
  Z/2^32 thinks: degree 1 (linear over integers)
  BTE says: degree 16 (4 layers of degree 2)

QUESTION: Does degree 16 give any computational advantage over degree 32?

For Gröbner basis / XL: complexity depends on degree.
  Degree D system in N variables: XL complexity ≈ N^(D+1).
  Degree 32: N^33 → astronomical.
  Degree 16: N^17 → still astronomical but 2^16× better?

But this is NAIVE. The actual question: can we SOLVE the system
layer by layer, spending degree-2 cost per layer?

Layer 0: degree-2 system, 2R-1 = 127 constraints, 512 variables.
  Solution space: 512 - 127 = 385 dimensions.
  Solving degree-2 over GF(2): XL cost ≈ C(385, 2) ≈ 74000.
  That's POLYNOMIAL. Not 2^385.

Layer 1 (given layer 0): degree-2, +127 constraints.
  Solution space: 512 - 254 = 258 dimensions.
  But: layer 1's degree-2 terms involve BOTH layer-1 variables
  AND layer-0 variables (which are now KNOWN from step 0).
  So effectively: degree-2 in 258 unknowns with 127 equations.
  XL cost ≈ C(258, 2) ≈ 33000.

Layer 2: 131 unknowns, 127 equations. XL ≈ C(131, 2) ≈ 8500.
Layer 3: 4 unknowns, 127 equations. XL ≈ C(4, 2) = 6.
Layer 4: 0 unknowns. Done.

BUT: this ignores the FACT that layer-0 has a 385-dim solution SPACE,
not a unique solution. We need to ENUMERATE or search this space.
The 385-dim space has 2^385 points.

So: layered solving = solve layer 0 in 2^385, then each subsequent
layer in polynomial time.

2^385 > 2^256 (brute force) → WORSE, not better.

UNLESS: we can CONSTRAIN layer 0 further using information from
higher layers.

What if we solve ALL layers simultaneously as a degree-2 system?
Total: 512 variables, 512 equations (at full rank), degree 2.
XL for degree-2 in 512 variables: need matrix of size C(512,2) ≈ 131000.
Gaussian elimination on 131000 × 131000 matrix → O(131000^3) ≈ 2^51.

Wait — that's O(2^51)? For SHA-256?

No — this is the XL cost for a RANDOM degree-2 system.
SHA-256's system is NOT random — it's structured.
And the degree-2 approximation ignores carries (which are ALSO degree 2
within our layer framework, but compose to degree 16).

Let me compute more carefully.
"""

import math


def experiment_xl_estimates():
    """
    Estimate XL/Gröbner basis complexity for SHA-256 viewed through BTE.

    The system:
      N = 512 message bits (unknowns)
      M = 512 hash bits (equations, at full rank)
      Degree: 2 within each layer, but layers compose

    For a degree-d system in N variables:
      XL at degree D: matrix size = C(N, D) × C(N, D)
      Need D such that the system becomes overdetermined at degree D.
      For random degree-2 system: D_reg ≈ N/M = 512/512 = 1.
      But degree-2 with D_reg=1 may not be solvable (depends on structure).

    Lazard's theorem: for a generic degree-2 system with M=N:
      XL needs to go to degree D_XL ≈ 2.
      Matrix: C(N+D_XL, D_XL) ≈ C(514, 2) ≈ 132000.
      Gaussian: 132000^ω where ω ≈ 2.37 → 132000^2.37 ≈ 2^41.
    """
    print("=" * 80)
    print("BTE COMPLEXITY: XL/Gröbner estimates for SHA-256")
    print("=" * 80)

    N = 512  # unknowns
    M = 512  # equations (at full rank from 4 layers × 127 + 4)

    # Degree-2 system parameters
    print(f"\n  SHA-256 as degree-2 system (BTE view):")
    print(f"    Unknowns: N = {N}")
    print(f"    Equations: M = {M}")
    print(f"    Degree: d = 2 (within each layer)")
    print(f"    Layers: D = 4 (carry depth)")

    # Number of monomials up to degree d
    def n_monomials(n, d):
        """Number of monomials in n variables up to degree d over GF(2)."""
        total = 0
        for k in range(d + 1):
            total += math.comb(n, k)
        return total

    mono_d2 = n_monomials(N, 2)
    print(f"\n  Monomials up to degree 2: C({N},0)+C({N},1)+C({N},2) = {mono_d2}")

    # XL at degree D: extend each equation to degree D by multiplying with monomials
    for D in [2, 3, 4]:
        n_ext = n_monomials(N, D)
        n_equations = M * n_monomials(N, D - 2)  # multiply M equations by degree D-2 monomials
        print(f"\n  XL at degree {D}:")
        print(f"    Extended monomials: {n_ext}")
        print(f"    Extended equations: {n_equations}")
        print(f"    Matrix: {n_equations} × {n_ext}")
        print(f"    Gaussian: ≈ {n_ext}^2.37 ≈ 2^{2.37 * math.log2(n_ext):.1f}")

    # BUT: SHA-256 is NOT a generic degree-2 system.
    # It has SPECIFIC structure:
    # - 4 layers of 127 equations each
    # - Each layer has degree 2 (Ch/Maj) in ~7 variables per equation
    # - Cross-layer coupling is LINEAR (rotations)
    # - Carry bridges connect layers

    print(f"\n  SHA-256 SPECIFIC structure:")
    print(f"    Each equation involves ~7 quadratic + ~6 linear variables")
    print(f"    Equations are SPARSE (not dense degree-2)")
    print(f"    Cross-layer coupling is LINEAR (reduces effective degree)")

    # Sparse system: XL complexity depends on the number of terms per equation,
    # not the total number of variables.
    # Each SHA-256 equation has ~13 variables (7 quadratic + 6 linear).
    # Quadratic terms: ~C(7,2) = 21 monomials per equation.
    # Total unique monomials across all equations: M × 21 = 10752,
    # but with heavy overlap (same variables appear in many equations).

    # For SPARSE degree-2 systems, specialized algorithms exist:
    # F4/F5 Gröbner: exploit structure.
    # But SHA-256's schedule creates SPECIFIC coupling that may or may not help.

    print(f"\n  Sparse structure per equation:")
    print(f"    Variables per equation: ~13 (7 NL + 6 L)")
    print(f"    Quadratic monomials per equation: ~C(7,2) = 21")
    print(f"    Total if independent: 512 × 21 = 10752")
    print(f"    But with overlap → much less")

    # THE REAL QUESTION: is there a way to exploit the layer structure
    # to solve the system faster than generic degree-2 methods?

    print(f"\n  ═══════════════════════════════════════════════════════════")
    print(f"  LAYERED SOLVING vs MONOLITHIC SOLVING")
    print(f"  ═══════════════════════════════════════════════════════════")
    print(f"")
    print(f"  MONOLITHIC (standard):")
    print(f"    Solve 512 degree-2 equations in 512 variables")
    print(f"    XL degree 3: matrix 2^{2.37 * math.log2(n_monomials(512,3)):.0f}")
    print(f"    This is HUGE. Impractical.")
    print(f"")
    print(f"  LAYERED (BTE approach):")
    print(f"    Layer 0: 127 degree-2 eqs in 512 vars → solution space dim 385")
    print(f"    Layer 1: 127 degree-2 eqs conditioned on layer 0")
    print(f"    Layer 2: 127 degree-2 eqs conditioned on layers 0-1")
    print(f"    Layer 3: 127 degree-2 eqs conditioned on layers 0-2")
    print(f"    Layer 4: 4 eqs → finalize")
    print(f"")
    print(f"  Problem: Layer 0 has 385-dim solution space.")
    print(f"  Can't enumerate 2^385 points.")
    print(f"  Need to carry the solution space SYMBOLICALLY to layer 1.")
    print(f"")
    print(f"  SYMBOLIC LAYERED SOLVING:")
    print(f"    Layer 0: parameterize 385 free variables.")
    print(f"    Layer 1: substitute, get 127 degree-2 eqs in 385 params.")
    print(f"    Solve: 127 eqs in 385 vars, degree 2.")
    print(f"    XL degree 2: matrix C(385,2) ≈ {math.comb(385,2)}")
    print(f"    Gaussian: {math.comb(385,2)}^2.37 ≈ 2^{2.37*math.log2(math.comb(385,2)):.0f}")
    print(f"")
    print(f"    Layer 2: 127 eqs in 258 vars, degree 2.")
    print(f"    XL: C(258,2) ≈ {math.comb(258,2)}")
    print(f"    Gaussian: 2^{2.37*math.log2(math.comb(258,2)):.0f}")
    print(f"")
    print(f"    Layer 3: 127 eqs in 131 vars, degree 2.")
    print(f"    XL: C(131,2) ≈ {math.comb(131,2)}")
    print(f"    Gaussian: 2^{2.37*math.log2(math.comb(131,2)):.0f}")
    print(f"")
    print(f"    Layer 4: 4 eqs in 4 vars. Trivial.")
    print(f"")
    print(f"  Total layered cost: max of per-layer costs = 2^{2.37*math.log2(math.comb(385,2)):.0f}")
    print(f"")
    print(f"  BUT: this is XL on degree-2. For SHA-256 the degree-2")
    print(f"  approximation requires KNOWING the carry values.")
    print(f"  Carry values = function of state = function of M.")
    print(f"  They're NOT known in advance.")
    print(f"")
    print(f"  OMEGA handles this by treating carry as variables.")
    print(f"  This adds ~31 variables per round × 64 rounds = 1984 extra variables.")
    print(f"  Total: 512 + 1984 = 2496 variables, degree 2.")
    print(f"  XL degree 2: C(2496, 2) ≈ {math.comb(2496, 2)}")
    print(f"  Gaussian: 2^{2.37*math.log2(math.comb(2496,2)):.0f}")
    print(f"")

    # COMPARISON TABLE
    print(f"  ═══════════════════════════════════════════════════════════")
    print(f"  COMPARISON")
    print(f"  ═══════════════════════════════════════════════════════════")
    print(f"")

    approaches = [
        ("Brute force", "2^256 (preimage) or 2^128 (collision)"),
        ("OMEGA (carry as vars)", f"XL on 2496 vars, deg 2 → 2^{2.37*math.log2(math.comb(2496,2)):.0f}"),
        ("BTE layered (layer 0)", f"XL on 385 vars, deg 2 → 2^{2.37*math.log2(math.comb(385,2)):.0f}"),
        ("BTE layered (all)", f"Max per-layer ≈ 2^{2.37*math.log2(math.comb(385,2)):.0f}"),
        ("Birthday (collision)", "2^128"),
    ]

    for name, cost in approaches:
        print(f"    {name:>30}: {cost}")

    print(f"\n  NOTE: XL estimates assume the system behaves like RANDOM degree-2.")
    print(f"  SHA-256's system is SPARSE and STRUCTURED → XL may be cheaper or more expensive.")
    print(f"  These are ORDER-OF-MAGNITUDE estimates, not exact.")


if __name__ == "__main__":
    experiment_xl_estimates()
