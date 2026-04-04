"""
BTE NONLINEAR: Stop decomposing H and V. Treat H∘V as ONE operation.

The failure of linear invariants tells us: separating H and V loses
the essential structure. SHA-256's power comes from their INTERACTION.

New approach: define a SINGLE operation ★ that captures
what one SHA-256 round does WITHOUT decomposing it.

★ is not addition. Not XOR. Not rotation. It's the ROUND ITSELF
viewed as an atomic operation on some new algebraic structure.

What properties must ★ have?
- ★ is a bijection (round function is invertible)
- ★ depends on parameters (K[r], W[r])
- Composing ★ 64 times gives SHA-256
- ★ has no simple invariants (we proved this)

But ★ DOES have structure:
- It's sparse: each output bit depends on ~5 input bits (measured!)
- The dependency graph has specific topology (measured!)
- The dependency CHANGES with the state (nonlinear)

What if the "algebra" is not about VALUES but about DEPENDENCIES?

Define: for each state s, the DEPENDENCY GRAPH D(s) is:
  D(s)[i][j] = 1 iff output bit i depends on input bit j
  This is exactly the Jacobian J(s).

The composition of two rounds: D(s2) ∘ D(s1) = matrix product over GF(2).
This product captures HOW dependencies propagate.

We already computed this — it's the Jacobian product.
But we only looked at EIGENSPACES.

What about the GRAPH STRUCTURE of the product?
Specifically: after k rounds, what is the DEPENDENCY GRAPH topology?
"""

import random
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, schedule
)


def compute_jacobian(state, r, W):
    """Compute 256x256 Jacobian (which input bits affect which output bits)."""
    a, b, c, d, e, f, g, h = state
    T1 = (h + Sig1(e) + Ch(e, f, g) + K[r] + W[r]) & MASK
    T2 = (Sig0(a) + Maj(a, b, c)) & MASK
    base_out = [(T1 + T2) & MASK, a, b, c, (d + T1) & MASK, e, f, g]

    base_bits = []
    for reg in base_out:
        for k in range(32):
            base_bits.append((reg >> k) & 1)

    J = []
    for reg_idx in range(8):
        for bit_idx in range(32):
            ts = list(state)
            ts[reg_idx] ^= (1 << bit_idx)
            ta, tb, tc, td, te, tf, tg, th = ts
            tT1 = (th + Sig1(te) + Ch(te, tf, tg) + K[r] + W[r]) & MASK
            tT2 = (Sig0(ta) + Maj(ta, tb, tc)) & MASK
            tout = [(tT1 + tT2) & MASK, ta, tb, tc, (td + tT1) & MASK, te, tf, tg]
            tbits = []
            for reg in tout:
                for k in range(32):
                    tbits.append((reg >> k) & 1)
            J.append([base_bits[i] ^ tbits[i] for i in range(256)])
    return J


def mat_mul_gf2(A, B, n):
    C = [[0]*n for _ in range(n)]
    for i in range(n):
        for k in range(n):
            if A[i][k]:
                for j in range(n):
                    C[i][j] ^= B[k][j]
    return C


def experiment_dependency_topology():
    """
    Study the TOPOLOGY of the dependency graph through rounds.

    For the Jacobian product of k rounds:
    - How many nonzero entries? (density)
    - What is the row/column weight distribution?
    - Are there HUBS (bits that influence everything)?
    - Are there ISOLATED bits (influenced by few)?
    """
    print("=" * 80)
    print("DEPENDENCY TOPOLOGY through rounds")
    print("=" * 80)

    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]
    states, W = sha256_round_trace(M)

    J_product = None
    for r in range(64):
        J_r = compute_jacobian(states[r], r, W)
        if J_product is None:
            J_product = J_r
        else:
            J_product = mat_mul_gf2(J_r, J_product, 256)

        n_rounds = r + 1
        if n_rounds in [1, 2, 4, 8, 16, 32, 64]:
            # Density
            total_ones = sum(sum(row) for row in J_product)
            density = total_ones / (256 * 256)

            # Row weights
            row_w = [sum(row) for row in J_product]
            min_rw = min(row_w)
            max_rw = max(row_w)
            avg_rw = sum(row_w) / 256

            # Column weights
            col_w = [sum(J_product[i][j] for i in range(256)) for j in range(256)]
            min_cw = min(col_w)
            max_cw = max(col_w)

            print(f"  {n_rounds:>2} rounds: density={density:.3f}, "
                  f"row_w=[{min_rw},{avg_rw:.0f},{max_rw}], "
                  f"col_w=[{min_cw},{max_cw}]")


def experiment_star_algebra():
    """
    Define ★ as the round function viewed abstractly.

    Instead of tracking bit values, track the TRANSFORMATION ITSELF.

    Two consecutive rounds: ★_r and ★_{r+1}.
    Their composition ★_{r+1} ∘ ★_r is another transformation.

    Question: Is there a PATTERN in how ★ transforms compose?

    Specifically: the Jacobians J_r vary with state, but do they
    share any common structure INDEPENDENT of state?

    Test: For 100 random states at round 20, compute J.
    What fraction of the 256×256 entries are ALWAYS 0 or ALWAYS 1?
    Those entries represent STRUCTURAL (state-independent) dependencies.
    """
    print("\n" + "=" * 80)
    print("★-ALGEBRA: State-independent structure of the round function")
    print("=" * 80)

    N = 200
    r = 20

    # Collect Jacobians for many different states at round 20
    always_one = [[True]*256 for _ in range(256)]  # J[i][j] = 1 for ALL states
    always_zero = [[True]*256 for _ in range(256)]  # J[i][j] = 0 for ALL states

    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        states, W = sha256_round_trace(M)
        J = compute_jacobian(states[r], r, W)

        for i in range(256):
            for j in range(256):
                if J[i][j] == 0:
                    always_one[i][j] = False
                if J[i][j] == 1:
                    always_zero[i][j] = False

    n_always_one = sum(sum(1 for j in range(256) if always_one[i][j]) for i in range(256))
    n_always_zero = sum(sum(1 for j in range(256) if always_zero[i][j]) for i in range(256))
    n_variable = 256*256 - n_always_one - n_always_zero

    print(f"  Round {r}, {N} random states:")
    print(f"    Always 1: {n_always_one}/{256*256} ({n_always_one/65536*100:.1f}%)")
    print(f"    Always 0: {n_always_zero}/{256*256} ({n_always_zero/65536*100:.1f}%)")
    print(f"    Variable: {n_variable}/{256*256} ({n_variable/65536*100:.1f}%)")

    # The always-0 entries = structural zeros (no dependency regardless of state)
    # The always-1 entries = structural dependencies (always depend)
    # The variable entries = carry-dependent (sometimes depend, sometimes not)

    # Show register-level breakdown
    reg_names = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
    print(f"\n  Register-level structural dependencies (always-1 count per 32x32 block):")
    print(f"    in\\out  " + " ".join(f"{n:>4}" for n in reg_names))
    for in_reg in range(8):
        counts = []
        for out_reg in range(8):
            c = 0
            for ib in range(32):
                for ob in range(32):
                    if always_one[in_reg*32+ib][out_reg*32+ob]:
                        c += 1
            counts.append(c)
        print(f"    {reg_names[in_reg]:>4}    " + " ".join(f"{c:>4}" for c in counts))

    print(f"\n  Register-level variable entries (carry-dependent, per 32x32 block):")
    print(f"    in\\out  " + " ".join(f"{n:>4}" for n in reg_names))
    for in_reg in range(8):
        counts = []
        for out_reg in range(8):
            c = 0
            for ib in range(32):
                for ob in range(32):
                    if not always_one[in_reg*32+ib][out_reg*32+ob] and not always_zero[in_reg*32+ib][out_reg*32+ob]:
                        c += 1
            counts.append(c)
        print(f"    {reg_names[in_reg]:>4}    " + " ".join(f"{c:>4}" for c in counts))

    # THE VARIABLE ENTRIES = where carries create state-dependent behavior.
    # THE STRUCTURAL ENTRIES = where the behavior is DETERMINISTIC.
    # The ratio variable/total = the ACTUAL "carry complexity" of SHA-256.

    total = 256 * 256
    print(f"\n  ★-complexity of SHA-256 round function:")
    print(f"    Structural (deterministic): {n_always_one + n_always_zero}/{total} = {(n_always_one+n_always_zero)/total:.1%}")
    print(f"    Carry-dependent (variable): {n_variable}/{total} = {n_variable/total:.1%}")
    print(f"")
    print(f"    This {n_variable/total:.1%} is the EXACT measure of carry's contribution.")
    print(f"    It's the part that no fixed linear algebra can capture.")
    print(f"    The rest ({(n_always_one+n_always_zero)/total:.1%}) is FIXED regardless of carries.")


if __name__ == "__main__":
    experiment_dependency_topology()
    experiment_star_algebra()
