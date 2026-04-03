"""
XL Linearization for the unified Q∩T system.

Converts quadratic GF(2) system to linear by introducing
auxiliary variables y_ij = x_i * x_j for each product term.

The quadratic equation:
  c ⊕ Σ x_i ⊕ Σ x_i*x_j = 0
becomes linear:
  c ⊕ Σ x_i ⊕ Σ y_ij = 0

Plus multiplication constraints (handled by iterative solving):
  y_ij = x_i * x_j

After linearization, Gaussian elimination on the enlarged system
reveals how much the quadratic constraints compress the kernel.
"""

from qt_solver.gf2 import gf2_gaussian_eliminate, gf2_solve, gf2_rank


class XLSystem:
    """
    Linearized Q∩T system.

    Original variables: x_0 .. x_{n-1}
    Auxiliary variables: y_k for each unique product (x_i * x_j)
    Total variables: n + num_products
    """

    def __init__(self, num_original_vars):
        self.n_orig = num_original_vars
        self.product_map = {}   # (i,j) -> auxiliary var index
        self.n_aux = 0
        self.rows = []          # Augmented GF(2) row vectors

    @property
    def num_vars(self):
        return self.n_orig + self.n_aux

    def _get_product_var(self, v1, v2):
        """Get or create auxiliary variable for x_v1 * x_v2."""
        key = (min(v1, v2), max(v1, v2))
        if key not in self.product_map:
            self.product_map[key] = self.n_orig + self.n_aux
            self.n_aux += 1
        return self.product_map[key]

    def add_linear_row(self, row, n_orig_vars):
        """
        Add a linear equation from the original system.
        The row uses original variable indices [0, n_orig).
        Augmented bit is at position n_orig_vars.
        """
        # Extract constant
        const = (row >> n_orig_vars) & 1
        # Copy variable bits, shift to new layout
        new_row = row & ((1 << n_orig_vars) - 1)
        if const:
            new_row |= (1 << self.num_vars)
        self.rows.append(new_row)

    def add_quadratic_eq(self, linear_row, quad_pairs, n_orig_vars):
        """
        Add a quadratic equation, linearizing the products.

        linear_row: augmented bit vector (original variable space)
        quad_pairs: list of (v1, v2) product pairs
        """
        const = (linear_row >> n_orig_vars) & 1
        new_row = linear_row & ((1 << n_orig_vars) - 1)

        for v1, v2 in quad_pairs:
            aux = self._get_product_var(v1, v2)
            new_row ^= (1 << aux)

        if const:
            new_row |= (1 << self.num_vars)
        self.rows.append(new_row)

    def finalize_rows(self):
        """
        Ensure all rows have consistent width after all products are known.
        Re-set the augmented column position.
        """
        total = self.num_vars
        final_rows = []
        for row in self.rows:
            # Find where the old augmented bit might be
            # We need to move it to position total
            # Old augmented positions could vary — re-extract
            # Actually, we set augmented bit at self.num_vars at add time,
            # but num_vars grows as we add products.
            # Simpler: store rows as (coeff_int, const) and rebuild.
            pass
        # The issue: as we add quadratic equations, num_vars grows,
        # so earlier rows have augmented bit at wrong position.
        # Fix: defer row construction.

    def solve(self):
        """
        Solve the linearized system.
        Returns (rank, kernel_dim, particular, kernel) or None.
        """
        n = self.num_vars
        # Fix augmented column: all rows should have const at bit n
        # Since we add rows with (1 << self.num_vars) and num_vars
        # might have changed, we need to fix this.
        fixed_rows = []
        for row in self.rows:
            # The constant bit could be anywhere above n_orig + n_aux_at_time
            # Extract all bits above the current variable space
            coeff = row & ((1 << n) - 1)
            # Check if any bit above n is set (that's the constant)
            const_bits = row >> n
            const = const_bits & 1
            new_row = coeff
            if const:
                new_row |= (1 << n)
            fixed_rows.append(new_row)

        rank = gf2_rank(fixed_rows, n)
        result = gf2_solve(fixed_rows, n)

        return {
            'rank': rank,
            'kernel_dim': n - rank,
            'num_vars': n,
            'num_orig_vars': self.n_orig,
            'num_aux_vars': self.n_aux,
            'num_equations': len(fixed_rows),
            'solution': result,
        }


def build_xl_from_unified(unified_system):
    """
    Build XL-linearized system from a UnifiedQTSystem.

    Takes all linear and quadratic equations and produces
    a fully linear system in the extended variable space.
    """
    us = unified_system
    n = us.vm.num_vars

    xl = XLSystem(n)

    # Add all linear equations
    for row in us.linear_rows:
        xl.add_linear_row(row, n)

    # Add all quadratic equations (linearized)
    for lin_row, quad_pairs, const in us.quad_equations:
        xl.add_quadratic_eq(lin_row, quad_pairs, n)

    return xl


def analyze_xl(num_rounds, seed=42):
    """
    Build and analyze the XL-linearized unified Q+T system.

    Reports:
    - Original variables vs extended (with auxiliaries)
    - Total equations
    - Rank and kernel dimension
    - Compression ratio vs original Q-only system
    """
    import random
    from qt_solver.sha256_traced import MASK32
    from qt_solver.unified import build_unified

    rng = random.Random(seed)
    msg = [rng.randint(0, MASK32) for _ in range(16)]

    system, trace = build_unified(msg, num_rounds)
    st = system.stats()

    print(f"\n--- XL Analysis: {num_rounds} rounds ---")
    print(f"Unified system: {st['num_vars']} vars, "
          f"{st['num_linear_eq']} linear + {st['num_quad_eq']} quad = "
          f"{st['total_eq']} total")

    xl = build_xl_from_unified(system)
    result = xl.solve()

    print(f"XL system: {result['num_vars']} vars "
          f"({result['num_orig_vars']} orig + {result['num_aux_vars']} aux)")
    print(f"Equations: {result['num_equations']}")
    print(f"Rank: {result['rank']}")
    print(f"Kernel dim: {result['kernel_dim']}")

    if result['solution'] is not None:
        _, kernel = result['solution']
        # Count how many kernel vectors touch only message bits
        msg_only = 0
        for kvec in kernel:
            non_msg = kvec >> 512
            if non_msg == 0:
                msg_only += 1
        print(f"Kernel vectors touching only message bits: {msg_only}/{len(kernel)}")

    return result
