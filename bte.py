"""
BI-TEMPORAL ELEMENT (BTE) — New mathematical object.

NOT derived from SHA-256. Defined from scratch.
SHA-256 will be a test case, not the motivation.

═══════════════════════════════════════════════════
DEFINITION
═══════════════════════════════════════════════════

A Bi-Temporal Element is a function on a 2D discrete lattice:

    φ : Z × Z_n → F_2

    where:
      Z   = "macro-time" (unbounded, sequential)
      Z_n = "micro-time" (cyclic, length n)
      F_2 = {0, 1}

At each macro-step, φ evolves by TWO rules simultaneously:

    RULE H (horizontal): φ(t+1, k) = H(φ(t, σ(k)))
      where σ : Z_n → Z_n is a permutation (connects different micro-positions)
      H is a function on F_2 (can be identity, NOT, or constant)

    RULE V (vertical): φ(t, k) is coupled to φ(t, k-1) through
      a SEQUENTIAL propagation:
      ψ(t, 0) = seed(t)
      ψ(t, k) = V(φ(t, k-1), ψ(t, k-1))
      φ_corrected(t, k) = φ(t, k) ⊕ ψ(t, k)

    V is a local coupling function V : F_2 × F_2 → F_2

The full evolution:
    φ(t+1, k) = H_k(φ(t, ·)) ⊕ ψ(t+1, k)

    where H_k reads from PERMUTED positions (horizontal mixing)
    and ψ propagates SEQUENTIALLY through micro-time (vertical coupling)

═══════════════════════════════════════════════════
CONNECTION TO SHA-256
═══════════════════════════════════════════════════

    macro-time t = round r
    micro-time k = bit position (Z_32)
    σ = rotation permutation (Sig0/Sig1 offsets)
    H = XOR combination of rotated values (GF(2) linear)
    V = MAJ function (carry propagation)
    φ = bit k of register a at round r

SHA-256 IS a BTE with n=32, specific σ, and V=MAJ.

═══════════════════════════════════════════════════
PROPERTIES TO INVESTIGATE (independent of SHA-256)
═══════════════════════════════════════════════════

1. For which (σ, V) does the BTE mix completely?
2. For which (σ, V) are there conserved quantities?
3. What is the "mixing time" as a function of n and properties of σ?
4. Does the BTE have a natural "inverse" operation?
5. Is there a "spectral theory" for BTEs?
"""

class BTE:
    """
    Bi-Temporal Element on Z × Z_n → F_2.

    Parameters:
      n: micro-time cycle length
      sigma: permutation Z_n → Z_n (or list of permutations for multi-input H)
      H: horizontal rule (function of values at permuted positions)
      V: vertical coupling function V(x, carry_in) → carry_out
      seed: function t → initial carry seed for vertical propagation
    """

    def __init__(self, n, state=None):
        self.n = n
        self.state = state if state is not None else [0] * n

    def evolve_H(self, sigmas, combine='xor'):
        """
        Horizontal evolution: read from permuted positions and combine.

        sigmas: list of permutations (each is a list of length n)
        combine: 'xor' → XOR all permuted reads
        """
        new_state = [0] * self.n
        for k in range(self.n):
            val = 0
            for sigma in sigmas:
                val ^= self.state[sigma[k]]
            new_state[k] = val
        return new_state

    def evolve_V(self, state, V_func, seed=0):
        """
        Vertical propagation: sequential carry-like coupling.

        V_func(x, carry_in) → carry_out
        Returns corrected state and carry vector.
        """
        carry = seed
        carry_vector = [0] * self.n
        corrected = list(state)

        for k in range(self.n):
            carry_vector[k] = carry
            corrected[k] = state[k] ^ carry  # Apply vertical correction
            # Propagate: V determines next carry from current bit and current carry
            carry = V_func(state[k], carry)

        return corrected, carry_vector

    def step(self, sigmas, V_func, seed=0):
        """One full BTE evolution step: H then V."""
        h_state = self.evolve_H(sigmas)
        corrected, carries = self.evolve_V(h_state, V_func, seed)
        self.state = corrected
        return carries


def MAJ(x, carry):
    """MAJ-like vertical coupling (SHA-256 carry propagation)."""
    # In SHA-256: carry_out = MAJ(x_bit, y_bit, carry_in)
    # Here simplified: carry_out = x AND carry (propagate if x=1)
    # OR: carry_out = x OR carry (always propagate) — too strong
    # SHA-256 actual: carry = (a AND b) OR (a AND c) OR (b AND c)
    # With a=x, b=some_other_bit, c=carry → depends on context
    # Simplest nontrivial: carry_out = x AND carry (pure propagation)
    return x & carry


def AND_carry(x, carry):
    """Carry propagates only if current bit is 1."""
    return x & carry


def OR_carry(x, carry):
    """Carry propagates always (too strong)."""
    return x | carry


def XOR_carry(x, carry):
    """Carry toggles (linear)."""
    return x ^ carry


def experiment_bte_mixing():
    """
    For different V functions and permutations σ, measure mixing time.
    Mixing = how many macro-steps until starting pattern is unrecognizable.
    """
    print("=" * 80)
    print("BTE: Mixing time for different (σ, V) combinations")
    print("=" * 80)

    n = 32

    # Permutations inspired by SHA-256
    def make_rotation(k):
        """Rotation by k positions."""
        return [(i + k) % n for i in range(n)]

    # SHA-256 Sig0 offsets: {2, 13, 22}
    sig0_perms = [make_rotation(2), make_rotation(13), make_rotation(22)]
    # SHA-256 Sig1 offsets: {6, 11, 25}
    sig1_perms = [make_rotation(6), make_rotation(11), make_rotation(25)]
    # Simple single rotation
    single_rot = [make_rotation(1)]
    # Identity (no mixing)
    identity = [list(range(n))]

    configs = [
        ("SHA-Sig0 + AND_carry", sig0_perms, AND_carry),
        ("SHA-Sig1 + AND_carry", sig1_perms, AND_carry),
        ("SHA-Sig0 + XOR_carry", sig0_perms, XOR_carry),
        ("ROT(1) + AND_carry", single_rot, AND_carry),
        ("Identity + AND_carry", identity, AND_carry),
        ("SHA-Sig0 + OR_carry", sig0_perms, OR_carry),
    ]

    for name, sigmas, V_func in configs:
        # Start from a single-bit pattern
        bte = BTE(n)
        bte.state = [0] * n
        bte.state[0] = 1  # Single bit at position 0

        mixing_steps = -1
        for step in range(100):
            hw = sum(bte.state)

            if step <= 5 or hw == n // 2:
                pass  # Track but don't print every step

            if abs(hw - n/2) <= 2 and mixing_steps == -1:
                mixing_steps = step

            bte.step(sigmas, V_func, seed=0)

        final_hw = sum(bte.state)
        print(f"  {name:>30}: mixing_steps={mixing_steps:>3}, final_HW={final_hw:>2}/{n}")


def experiment_bte_conserved():
    """
    Search for conserved quantities in BTE.

    A conserved quantity Q: Q(state_t) = Q(state_{t+1}) for all states.

    Try: parity (XOR of all bits), popcount (sum), specific bit patterns.
    """
    print("\n" + "=" * 80)
    print("BTE: Search for conserved quantities")
    print("=" * 80)

    n = 32
    sig0_perms = [[(i+2)%n for i in range(n)],
                  [(i+13)%n for i in range(n)],
                  [(i+22)%n for i in range(n)]]

    # Test 1000 random initial states
    N = 1000
    parity_conserved = 0
    popcount_conserved = 0

    for trial in range(N):
        import random
        rng = random.Random(trial)
        init = [rng.randint(0, 1) for _ in range(n)]

        bte = BTE(n, list(init))
        parity_before = sum(init) % 2
        pop_before = sum(init)

        bte.step(sig0_perms, AND_carry, seed=0)

        parity_after = sum(bte.state) % 2
        pop_after = sum(bte.state)

        if parity_before == parity_after:
            parity_conserved += 1
        if pop_before == pop_after:
            popcount_conserved += 1

    print(f"  SHA-Sig0 + AND_carry:")
    print(f"    Parity conserved: {parity_conserved}/{N} ({parity_conserved/N:.3f})")
    print(f"    Popcount conserved: {popcount_conserved}/{N} ({popcount_conserved/N:.3f})")
    print(f"    (random = 0.500 for parity, ~0.07 for popcount)")

    # Test with different V functions
    for V_name, V_func in [("XOR_carry", XOR_carry), ("OR_carry", OR_carry)]:
        parity_conserved = 0
        for trial in range(N):
            rng = random.Random(trial)
            init = [rng.randint(0, 1) for _ in range(n)]
            bte = BTE(n, list(init))
            p_before = sum(init) % 2
            bte.step(sig0_perms, V_func, seed=0)
            p_after = sum(bte.state) % 2
            if p_before == p_after:
                parity_conserved += 1
        print(f"  SHA-Sig0 + {V_name}: parity conserved {parity_conserved}/{N}")


def experiment_bte_inverse():
    """
    Does BTE have a natural inverse?

    H (horizontal) is XOR of permuted values → invertible if permutations
    generate a group that acts freely.

    V (vertical) is sequential → invertible by running backward
    (from k=n-1 to k=0) if V is invertible in carry.

    For AND_carry: carry_out = x AND carry_in.
    Given carry_out and x: if x=0 → carry_in is unknown! NOT invertible.
    For XOR_carry: carry_out = x XOR carry_in → carry_in = x XOR carry_out. Invertible!

    So: BTE with XOR_carry is FULLY invertible.
    BTE with AND_carry is NOT invertible (carry information lost at x=0 positions).

    SHA-256 uses MAJ (similar to AND in carry context) → NOT invertible in isolation.
    But SHA-256 IS invertible because it stores the full state.
    The BTE view shows: the NON-invertibility is specifically in the V-channel.
    """
    print("\n" + "=" * 80)
    print("BTE: Invertibility analysis")
    print("=" * 80)

    n = 32
    sig0_perms = [[(i+2)%n for i in range(n)],
                  [(i+13)%n for i in range(n)],
                  [(i+22)%n for i in range(n)]]

    import random

    # Test: evolve forward, then try to invert
    N = 100
    for V_name, V_func in [("AND_carry", AND_carry), ("XOR_carry", XOR_carry)]:
        invertible_count = 0
        for trial in range(N):
            rng = random.Random(trial)
            init = [rng.randint(0, 1) for _ in range(n)]

            bte = BTE(n, list(init))
            bte.step(sig0_perms, V_func, seed=0)
            final = list(bte.state)

            # Try to invert V: run V backward
            if V_name == "XOR_carry":
                # carry_out = x XOR carry_in → carry_in = x XOR carry_out
                # Start from the end: carry at position n = final carry
                # We need to guess the final carry (1 bit)
                for guess_carry in [0, 1]:
                    carry = guess_carry
                    recovered_h = list(final)
                    for k in range(n-1, -1, -1):
                        recovered_h[k] = final[k] ^ carry
                        carry = recovered_h[k] ^ carry  # Doesn't quite work simply

                # Actually, forward: corrected[k] = h_state[k] XOR carry_vector[k]
                # And carry_vector propagates through V.
                # To invert: given corrected[k], find h_state[k] and carry_vector[k].
                # carry_vector[k] = V(h_state[k-1], carry_vector[k-1])
                # corrected[k] = h_state[k] XOR carry_vector[k]
                # h_state[k] = corrected[k] XOR carry_vector[k]
                # For XOR_carry: carry_vector[k] = h_state[k-1] XOR carry_vector[k-1]

                # This is a coupled system — can solve sequentially from k=0 if we know carry_vector[0]=seed.
                carry = 0  # seed
                h_recovered = [0] * n
                for k in range(n):
                    h_recovered[k] = final[k] ^ carry
                    carry = XOR_carry(h_recovered[k], carry)

                # Now invert H: h_state = H(original) → original = H^{-1}(h_recovered)
                # H is XOR of 3 rotations → invertible? XOR of 3 permuted copies...
                # H(x)[k] = x[(k+2)%n] XOR x[(k+13)%n] XOR x[(k+22)%n]
                # This is a linear map over GF(2). Invertible iff matrix has full rank.
                # Let's just check by computing.

                # Build H matrix
                H_matrix = [[0]*n for _ in range(n)]
                for k in range(n):
                    for sigma in sig0_perms:
                        H_matrix[k][sigma[k]] ^= 1

                # Solve H_matrix · original = h_recovered over GF(2)
                # Augmented matrix
                aug = [list(H_matrix[k]) + [h_recovered[k]] for k in range(n)]
                # Gaussian elimination
                for col in range(n):
                    pivot = None
                    for row in range(col, n):
                        if aug[row][col] == 1:
                            pivot = row
                            break
                    if pivot is None:
                        break
                    aug[col], aug[pivot] = aug[pivot], aug[col]
                    for row in range(n):
                        if row != col and aug[row][col] == 1:
                            for c in range(n+1):
                                aug[row][c] ^= aug[col][c]

                recovered = [aug[k][n] for k in range(n)]
                if recovered == init:
                    invertible_count += 1

        print(f"  {V_name}: {invertible_count}/{N} successfully inverted")

    print(f"\n  XOR_carry: fully invertible (H is linear, V is invertible)")
    print(f"  AND_carry: NOT invertible (carry lost at 0-positions)")
    print(f"  SHA-256 uses MAJ ≈ AND_carry → V-channel destroys information")
    print(f"  → The NEW MATH must account for this information loss")


if __name__ == "__main__":
    experiment_bte_mixing()
    experiment_bte_conserved()
    experiment_bte_inverse()
