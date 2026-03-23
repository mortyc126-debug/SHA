#!/usr/bin/env python3
"""Step 17a: Find fixed points of Sigma0 and Sigma1.

Nikolic-Biryukov key insight: if Sigma1(x) = x, then
Sigma1(x) - Sigma1(y) = x - y when both are fixed points.
This makes the nonlinear Sigma function preserve modular differences.

We need to find ALL x such that Sigma0(x) = x and Sigma1(x) = x.
"""

MASK32 = 0xFFFFFFFF

def rotr(x, n):
    return ((x >> n) | (x << (32 - n))) & MASK32

def Sigma0(x):
    return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)

def Sigma1(x):
    return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)

def hw(x):
    return bin(x).count('1')

# ============================================================
# Part 1: Brute-force search for fixed points (32-bit space)
# ============================================================
def find_fixed_points():
    """Search all 2^32 values for Sigma fixed points."""
    print("=" * 70)
    print("Part 1: Fixed points of Sigma0 and Sigma1")
    print("=" * 70)

    # Sigma0(x) = ROTR2(x) ^ ROTR13(x) ^ ROTR22(x)
    # Sigma0(x) = x means: ROTR2(x) ^ ROTR13(x) ^ ROTR22(x) = x
    # This is a linear equation over GF(2)!
    # So we can solve it with linear algebra.

    # Build the 32x32 matrix for Sigma0 over GF(2)
    # Sigma0 is a GF(2)-linear function (XOR of rotations)
    # Fixed points: (Sigma0 + I)x = 0, i.e., kernel of (Sigma0 + I)

    import numpy as np

    def build_matrix(func):
        """Build 32x32 GF(2) matrix for a function."""
        M = np.zeros((32, 32), dtype=int)
        for j in range(32):
            val = func(1 << j)
            for i in range(32):
                M[i][j] = (val >> i) & 1
        return M

    def gf2_kernel(M):
        """Find kernel of M over GF(2) using Gaussian elimination."""
        n = M.shape[0]
        m = M.shape[1]
        A = M.copy() % 2

        # Augment for tracking
        pivot_cols = []
        row = 0
        for col in range(m):
            # Find pivot
            found = False
            for r in range(row, n):
                if A[r][col] == 1:
                    A[[row, r]] = A[[r, row]]
                    found = True
                    break
            if not found:
                continue
            pivot_cols.append(col)
            # Eliminate
            for r in range(n):
                if r != row and A[r][col] == 1:
                    A[r] = (A[r] + A[row]) % 2
            row += 1

        # Free variables = columns not in pivot_cols
        free_cols = [c for c in range(m) if c not in pivot_cols]
        rank = len(pivot_cols)

        # Build kernel basis
        kernel_basis = []
        for fc in free_cols:
            vec = np.zeros(m, dtype=int)
            vec[fc] = 1
            for i, pc in enumerate(pivot_cols):
                vec[pc] = A[i][fc]
            kernel_basis.append(vec)

        return kernel_basis, rank

    for name, func in [("Sigma0", Sigma0), ("Sigma1", Sigma1)]:
        M = build_matrix(func)
        # Fixed points: (M + I)x = 0
        MI = (M + np.eye(32, dtype=int)) % 2

        kernel, rank = gf2_kernel(MI)
        dim = len(kernel)

        print(f"\n  {name}:")
        print(f"    Matrix rank of (Sigma + I): {rank}")
        print(f"    Kernel dimension: {dim}")
        print(f"    Number of fixed points: 2^{dim} = {2**dim}")

        if dim > 0 and dim <= 10:
            # Enumerate all fixed points
            fps = []
            for mask in range(2**dim):
                x = 0
                for bit in range(dim):
                    if (mask >> bit) & 1:
                        for j in range(32):
                            x ^= (kernel[bit][j] << j)
                x = int(x) & MASK32
                fps.append(x)
                # Verify
                assert func(x) == x, f"NOT a fixed point: {x:#010x}"

            fps.sort()
            print(f"    Fixed points:")
            for x in fps:
                print(f"      0x{x:08x} (HW={hw(x):2d})")
        elif dim > 10:
            print(f"    Too many fixed points to enumerate (2^{dim})")
            # Show basis vectors
            print(f"    Kernel basis vectors:")
            for i, vec in enumerate(kernel):
                val = 0
                for j in range(32):
                    val |= (int(vec[j]) << j)
                print(f"      v{i}: 0x{val:08x} (HW={hw(val)})")

# ============================================================
# Part 2: Near-fixed-points — Sigma(x) close to x
# ============================================================
def find_near_fixed_points():
    """Analyze additive differential preservation at fixed points."""
    print("\n" + "=" * 70)
    print("Part 2: Additive preservation analysis")
    print("=" * 70)

    import random
    random.seed(42)

    # Known fixed points
    sigma1_fps = [0x00000000, 0x33333333, 0x55555555, 0x66666666,
                  0x99999999, 0xaaaaaaaa, 0xcccccccc, 0xffffffff]
    sigma0_fps = [0x00000000, 0xffffffff]

    for name, func, fps in [("Sigma0", Sigma0, sigma0_fps),
                             ("Sigma1", Sigma1, sigma1_fps)]:
        print(f"\n  {name} additive preservation for delta=1 (random x):")
        good_count = 0
        N = 100000
        for _ in range(N):
            x = random.randint(0, MASK32)
            diff = (func((x + 1) & MASK32) - func(x)) & MASK32
            if diff == 1:
                good_count += 1
        print(f"    P(Sigma(x+1) - Sigma(x) = 1) = {good_count/N:.6f}")

        print(f"\n  {name} additive diff AT fixed points (delta=1):")
        for x in fps:
            d = (func((x + 1) & MASK32) - func(x)) & MASK32
            print(f"    x=0x{x:08x}: Sigma(x+1)-Sigma(x) = 0x{d:08x} (HW={hw(d)})")

# ============================================================
# Part 3: Differential preservation analysis
# ============================================================
def differential_preservation():
    """For a local collision with delta=+1, analyze how Sigma
    preserves or distorts the difference."""
    print("\n" + "=" * 70)
    print("Part 3: Sigma differential preservation for delta=+1")
    print("=" * 70)

    import random
    random.seed(42)

    for name, func in [("Sigma0", Sigma0), ("Sigma1", Sigma1)]:
        # For the local collision, we need:
        # Sigma(x + delta) - Sigma(x) = delta (mod 2^32)
        # This is EXACT (no approximation) when x is a "good" value

        # Probability that Sigma preserves +1 difference
        # depends on carry propagation

        print(f"\n  {name} difference preservation for various deltas:")
        for delta in [1, -1 & MASK32, 2, 0x80000000]:
            good = 0
            N = 200000
            for _ in range(N):
                x = random.randint(0, MASK32)
                actual_diff = (func((x + delta) & MASK32) - func(x)) & MASK32
                if actual_diff == delta:
                    good += 1
            p = good / N
            if p > 0:
                log_p = f"2^{__import__('math').log2(p):.1f}"
            else:
                log_p = "< 2^-17"
            print(f"    delta=0x{delta:08x}: P(preserved) = {p:.6f} ({log_p})")

        # For the Nikolic-Biryukov approach:
        # We don't need Sigma to preserve delta for ALL x.
        # We need to CHOOSE x such that it works.
        # The sufficient condition is related to carry-free addition.

        print(f"\n  {name} carry-free analysis for delta=1:")
        # Adding 1 to x: carry propagates through consecutive 1-bits
        # x + 1 = x XOR (mask of flipped bits)
        # The flipped bits are: bit 0 always, plus consecutive 1-bits
        # Sigma is GF2-linear, so Sigma(x XOR d) = Sigma(x) XOR Sigma(d)
        # But addition is not XOR when carries exist

        # When x has bit 0 = 0: x+1 = x XOR 1, no carry
        # Sigma(x+1) = Sigma(x) XOR Sigma(1)
        # Sigma(x+1) - Sigma(x) = Sigma(1) - 2*(Sigma(x) AND Sigma(1))
        # For this to equal 1, we need specific conditions

        sig_1 = func(1)
        print(f"    {name}(1) = 0x{sig_1:08x} (HW={hw(sig_1)})")

        # When bit 0 of x is 0 (no carry):
        # D_add = Sigma(x XOR 1) - Sigma(x)
        #       = (Sigma(x) XOR Sigma(1)) - Sigma(x)
        #       = Sigma(1) - 2*(Sigma(x) AND Sigma(1))
        # For D_add = 1:
        # Sigma(1) - 2*(Sigma(x) AND Sigma(1)) = 1
        # Let s = Sigma(1), then need:
        # s - 2*(Sigma(x) AND s) = 1
        # => Sigma(x) AND s = (s - 1) / 2

        # This constrains specific bits of Sigma(x)
        target = ((sig_1 - 1) >> 1) & MASK32
        constraint_bits = hw(sig_1)
        print(f"    Constraint: Sigma(x) AND Sigma(1) = 0x{target:08x}")
        print(f"    Number of constrained bits: {constraint_bits}")
        print(f"    Expected P(satisfied | bit0=0) ≈ 2^-{constraint_bits}")

        # Verify empirically
        good = 0
        N = 500000
        for _ in range(N):
            x = random.randint(0, MASK32) & ~1  # Force bit 0 = 0
            actual_diff = (func((x + 1) & MASK32) - func(x)) & MASK32
            if actual_diff == 1:
                good += 1
        print(f"    Empirical P = {good/N:.6f} (2^{__import__('math').log2(good/N) if good else -99:.1f})")


if __name__ == "__main__":
    find_fixed_points()
    find_near_fixed_points()
    differential_preservation()

    print("\n" + "=" * 70)
    print("Step 17a Complete")
    print("=" * 70)
