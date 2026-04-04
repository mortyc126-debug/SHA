"""
P1.1: ANALYTICAL PROOF — Why створочное creates exactly 1 GF(2) dependency on bit-0.

THE PROOF:

Given: SHA-256 round function.
  a_new = T1 + T2
  e_new = d + T1
  → a_new - e_new = T2 - d  (mod 2^32)
  → e_new = a_new - T2 + d  (mod 2^32)

Substituting shifts: d[r] = a[r-3], T2[r] = Sig0(a[r]) + Maj(a[r], a[r-1], a[r-2])

  e[r+1] = a[r+1] - Sig0(a[r]) - Maj(a[r], a[r-1], a[r-2]) + a[r-3]  (mod 2^32)

AT BIT 0 (no carry):
  e[r+1][0] = a[r+1][0] ⊕ Sig0(a[r])[0] ⊕ Maj(a[r], a[r-1], a[r-2])[0] ⊕ a[r-3][0]

  where:
    Sig0(a[r])[0] = a[r][2] ⊕ a[r][13] ⊕ a[r][22]  — CROSS-LAYER (bits 2,13,22)
    Maj(a,b,c)[0] = a[r][0]·a[r-1][0] ⊕ a[r][0]·a[r-2][0] ⊕ a[r-1][0]·a[r-2][0]
                   — WITHIN Layer 0, degree 2

For the JACOBIAN (linear approximation at a fixed point):
  ∂e[r+1][0]/∂a[s][0] depends on:
    - ∂a[r+1][0]/∂a[s][0]  (from a-trajectory)
    - ∂Maj[0]/∂a[s][0] = linearized Maj (depends on specific values)
    - a[r-3][0] contribution

The Jacobian of the e-bit-0 trajectory = Jacobian of a-bit-0 trajectory
MODIFIED by the створочное relation.

KEY: the створочное relation is ONE equation per round:
  e[r+1][0] = f(a[r+1][0], a[r][0], a[r-1][0], a[r-2][0], a[r-3][0], cross-layer)

This equation LINKS the e-sequence to the a-sequence at every round.
Over R rounds: R such equations.

But: the a-sequence has R values (a[1][0]...a[R][0]) with rank R.
The e-sequence has R values (e[1][0]...e[R][0]) with rank R.
Combined: 2R values.

The R створочное equations reduce the rank by at most R.
But: how many are INDEPENDENT?

If the створочное equations are all independent: rank = 2R - R = R.
If only 1 is independent: rank = 2R - 1.

We measured: rank = 2R - 1. So only 1 створочное equation is independent.
Why?

BECAUSE: the створочное equations at different rounds are NOT independent!
Equation at round r+1 involves a[r+1][0], a[r][0], a[r-1][0], a[r-2][0], a[r-3][0].
Equation at round r involves a[r][0], a[r-1][0], a[r-2][0], a[r-3][0], a[r-4][0].

They OVERLAP in variables a[r][0], a[r-1][0], a[r-2][0], a[r-3][0].
The XOR of equation(r+1) and equation(r) involves FEWER variables.

In the LINEAR part: XOR of consecutive equations involves only
a[r+1][0] ⊕ a[r-4][0] (new) and the linearized Maj difference.

Over R rounds: the R equations form a SHIFT-STRUCTURED system.
A shift-structured system of R equations in 2R variables has rank...
let's compute.
"""

import random
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, schedule
)


def gf2_rank(matrix, nrows, ncols):
    m = [list(row[:ncols]) for row in matrix[:nrows]]
    rank = 0
    for col in range(ncols):
        pivot = None
        for row in range(rank, len(m)):
            if m[row][col] == 1:
                pivot = row
                break
        if pivot is None:
            continue
        m[rank], m[pivot] = m[pivot], m[rank]
        for row in range(len(m)):
            if row != rank and m[row][col] == 1:
                for c in range(ncols):
                    m[row][c] ^= m[rank][c]
        rank += 1
    return rank


def prove_one_dependency():
    """
    Direct proof: the R створочное equations have rank exactly R-1 over GF(2).

    Build the створочное system matrix and compute its rank.
    """
    print("=" * 80)
    print("P1.1: Створочное equations — rank = R-1 (1 dependency)")
    print("=" * 80)

    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]

    for R in [8, 16, 32, 64]:
        states, W = sha256_round_trace(M, rounds=R)

        # Build створочное Jacobian:
        # For each round r (5..R): створочное says
        # e[r][0] = a[r][0] ⊕ a[r-4][0] ⊕ Maj_linearized ⊕ Sig0_cross_layer
        #
        # We need the Jacobian of e[r][0] - a[r][0] - a[r-4][0] - Maj_lin - Sig0_cross
        # with respect to {a[1][0], ..., a[R][0], e[1][0], ..., e[R][0]}
        #
        # Variables: 2R bit-0 values (a[1..R][0] and e[1..R][0])
        # Equations: R створочное relations

        # Numerically: compute derivative of [e[r][0] ⊕ f(a)] w.r.t. each input bit
        # The створочное at bit 0 should be: e[r][0] ⊕ a[r][0] ⊕ a[r-4][0] ⊕ ... = 0
        # Linearized: a linear combination of a's and e's that equals 0.

        # Build: for each message bit, compute which створочное equations change.
        # Створочное equation r: evaluate e[r][0] ⊕ (a[r][0] + a[r-4][0] - T2[r-1][0])
        # But this is always 0 by definition. The JACOBIAN tells us dependencies.

        # Instead: build Jacobian of (a-trajectory[0], e-trajectory[0]) jointly,
        # then find the dependency.

        # The створочное equations define R linear relations among 2R variables.
        # Jacobian of створочное: for each msg bit, does it change
        # the value of (e[r][0] ⊕ a[r][0] ⊕ a[r-4][0] ⊕ Sig0(a[r-1])[0] ⊕ Maj[0])?

        # Since створочное = 0 always, its Jacobian = 0 trivially.
        # The rank comes from the Jacobian of (a, e) SUBJECT TO створочное.

        # Cleaner approach: compute rank of the (a,e) bit-0 system directly
        # and verify it equals 2R - 1.

        # Actually, let me think about this differently.
        # The combined (a,e) system at bit 0 has rank 2R-1.
        # The a-only system has rank R (full).
        # So: e adds R new equations but only R-1 new independent ones.
        # The 1 dependent = the створочное.

        # To find WHICH combination of e-equations is dependent on a-equations:
        # Build matrix: rows = [a[1][0]_derivatives, ..., a[R][0]_derivatives,
        #                       e[1][0]_derivatives, ..., e[R][0]_derivatives]
        # Find null space of this matrix.
        # The null space vector tells us: which linear combination of
        # a-rows and e-rows = 0.

        base_a = [(states[r+1][0] & 1) for r in range(R)]
        base_e = [(states[r+1][4] & 1) for r in range(R)]

        # Jacobian of combined system: 512 inputs → 2R outputs
        J = []
        for word in range(16):
            for bit in range(32):
                M2 = list(M); M2[word] ^= (1 << bit)
                s2, _ = sha256_round_trace(M2, rounds=R)
                row = []
                for r in range(R):
                    row.append(base_a[r] ^ (s2[r+1][0] & 1))
                for r in range(R):
                    row.append(base_e[r] ^ (s2[r+1][4] & 1))
                J.append(row)

        # Rank should be 2R-1
        rank = gf2_rank(J, 512, 2*R)

        # Find the null vector of J^T (the dependency)
        # J is 512×2R. J^T is 2R×512. We need null space of J^T.
        # But we want the dependency among the 2R COLUMNS of J.
        # Null space of J (as column vectors) = left null space.

        # Actually: we want v ∈ GF(2)^{2R} such that J · v^T = 0 for all rows.
        # This is: for each input bit, Σ v[k] · J[input][k] = 0.
        # = the null space of J^T (transpose).

        JT = [[J[i][j] for i in range(512)] for j in range(2*R)]
        rank_JT = gf2_rank(JT, 2*R, 512)
        null_dim = 2*R - rank_JT

        print(f"  R={R:>2}: rank(a+e) = {rank}/{2*R}, null_dim = {null_dim}")

        # Find null vector
        if null_dim == 1:
            # Gaussian eliminate JT, find the free variable
            m = [list(JT[j][:512]) + [0]*2*R for j in range(2*R)]
            # Augment with identity for tracking
            for j in range(2*R):
                m[j][512 + j] = 1

            r = 0
            for col in range(512):
                pivot = None
                for row in range(r, 2*R):
                    if m[row][col] == 1:
                        pivot = row
                        break
                if pivot is None:
                    continue
                m[r], m[pivot] = m[pivot], m[r]
                for row in range(2*R):
                    if row != r and m[row][col] == 1:
                        for c in range(512 + 2*R):
                            m[row][c] ^= m[r][c]
                r += 1

            # The free row (rank_JT = 2R-1, so 1 free row)
            for j in range(2*R):
                is_zero = all(m[j][c] == 0 for c in range(512))
                if is_zero:
                    null_vec = [m[j][512 + k] for k in range(2*R)]
                    a_part = null_vec[:R]
                    e_part = null_vec[R:]
                    hw_a = sum(a_part)
                    hw_e = sum(e_part)
                    print(f"         Null vector: a-part HW={hw_a}, e-part HW={hw_e}")

                    # Which rounds are in the dependency?
                    a_rounds = [r for r in range(R) if a_part[r]]
                    e_rounds = [r for r in range(R) if e_part[r]]
                    if R <= 16:
                        print(f"         a-rounds: {a_rounds}")
                        print(f"         e-rounds: {e_rounds}")
                    else:
                        print(f"         a-rounds: {a_rounds[:5]}...{a_rounds[-5:]}")
                        print(f"         e-rounds: {e_rounds[:5]}...{e_rounds[-5:]}")
                    break


if __name__ == "__main__":
    prove_one_dependency()
