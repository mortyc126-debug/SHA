"""
COMBINE FACTS: Bit 0 (carry-free) + Schedule (GF(2)-linear) = ?

Fact A: Bit 0 of any addition = pure XOR (no carry-in ever)
Fact B: Schedule W[16..63] = GF(2)-linear function of W[0..15]

Both are exact. Both are GF(2). Both work on all 64 rounds.

ATTEMPT: Track ONLY bit 0 of every register through all 64 rounds.
At each addition, bit 0 = pure XOR. No carry.
Through rotations, bit 0 picks up bits from positions 2, 6, 11, 13, 22, 25.
Those positions DO have carries. But at bit 0 ITSELF — never.

So: define φ₀(r) = bit 0 of a[r].
    φ₀(r+1) = f(bits of state[r] at specific positions) — exact GF(2) degree 2.

The bits at positions 2, 6, etc. ARE contaminated by carry.
But φ₀ ITSELF is a clean GF(2) function of those contaminated bits.

What if we don't care about the carry contamination of the inputs?
What if we just track φ₀ as a BLACK BOX function of M?

φ₀(r) is a specific Boolean function of M (512 bits).
For r=1: degree 2 function of ~15 bits.
For r=2: degree ? function of ~30 bits.
For r=64: degree ? function of ~256 bits.

KEY: at each round, the UPDATE of φ₀ is degree 2.
But the COMPOSITION of 64 degree-2 updates is degree... 2^64?

No! The degree doesn't multiply — it only multiplies through Ch/Maj
(degree 2 per round), but the carry contamination of INPUT bits
means the actual degree of φ₀(r) in terms of M grows differently.

Let's measure the ACTUAL algebraic properties of φ₀ across rounds.
"""

import random
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, schedule
)


def experiment_bit0_composition():
    """
    Track bit 0 of a[r] as a function of M.
    Measure: for each round, what is the ALGEBRAIC DEGREE of φ₀(r)?

    Method: compute the Walsh-Hadamard spectrum or simpler:
    test linearity (degree 1), quadraticity (degree 2), etc.
    """
    print("=" * 80)
    print("FACT COMBINATION: Bit 0 across rounds — algebraic degree")
    print("=" * 80)

    rng = random.Random(42)
    M_base = [rng.randint(0, MASK) for _ in range(16)]

    for R in [1, 2, 3, 4, 6, 8, 16, 32, 64]:
        states_base, _ = sha256_round_trace(M_base, rounds=R)
        base_bit = states_base[R][0] & 1

        # Test degree 1: is f(M XOR e_i XOR e_j) = f(M) XOR f(M XOR e_i) XOR f(M XOR e_j)?
        # If yes for all i,j → degree 1 (affine).
        # If no → degree ≥ 2.

        # Test degree 2: is f(M XOR e_i XOR e_j XOR e_k) =
        #   f(M) XOR f(MXORe_i) XOR f(MXORe_j) XOR f(MXORe_k)
        #   XOR f(MXORe_iXORe_j) XOR f(MXORe_iXORe_k) XOR f(MXORe_jXORe_k)?
        # If yes for all i,j,k → degree ≤ 2.

        # Quick test with random triples
        N_test = 200
        degree_1_violations = 0
        degree_2_violations = 0

        for trial in range(N_test):
            rng2 = random.Random(trial * 1000 + R)
            # Pick random bit positions (word, bit)
            w1, b1 = rng2.randint(0, 3), rng2.randint(0, 31)
            w2, b2 = rng2.randint(0, 3), rng2.randint(0, 31)
            w3, b3 = rng2.randint(0, 3), rng2.randint(0, 31)

            def flip(M, w, b):
                M2 = list(M)
                M2[w] ^= (1 << b)
                return M2

            def f(M):
                s, _ = sha256_round_trace(M, rounds=R)
                return s[R][0] & 1

            f0 = f(M_base)
            f1 = f(flip(M_base, w1, b1))
            f2 = f(flip(M_base, w2, b2))
            f12 = f(flip(flip(M_base, w1, b1), w2, b2))

            # Degree 1 test: f(x+a+b) + f(x+a) + f(x+b) + f(x) = 0 (mod 2)
            d1_check = f0 ^ f1 ^ f2 ^ f12
            if d1_check != 0:
                degree_1_violations += 1

            # Degree 2 test
            f3 = f(flip(M_base, w3, b3))
            f13 = f(flip(flip(M_base, w1, b1), w3, b3))
            f23 = f(flip(flip(M_base, w2, b2), w3, b3))
            f123 = f(flip(flip(flip(M_base, w1, b1), w2, b2), w3, b3))

            d2_check = f0 ^ f1 ^ f2 ^ f3 ^ f12 ^ f13 ^ f23 ^ f123
            if d2_check != 0:
                degree_2_violations += 1

        d1_rate = degree_1_violations / N_test
        d2_rate = degree_2_violations / N_test

        if d1_rate == 0:
            deg = "= 1 (affine)"
        elif d2_rate == 0:
            deg = "= 2 (quadratic)"
        elif d2_rate < 0.1:
            deg = f"≈ 2 ({d2_rate:.1%} violations)"
        elif d2_rate < 0.4:
            deg = f"≈ 3-4 ({d2_rate:.1%} degree-2 violations)"
        else:
            deg = f"≥ 5 ({d2_rate:.1%} = near random)"

        print(f"  R={R:>2}: degree-1 viol={d1_rate:.3f}, degree-2 viol={d2_rate:.3f} → degree {deg}")


def experiment_multibit0_system():
    """
    Instead of one bit, take bit 0 of ALL 8 registers.
    This gives 8 Boolean functions of M at each round.

    Together: 8 × 64 = 512 Boolean functions.
    M has 512 bits.

    If these 512 functions are "structured" (low degree, sparse),
    they define a system that might be solvable.

    Test: what is the GF(2) rank of the system "bit 0 of all 8 regs at all 64 rounds"?
    """
    print("\n" + "=" * 80)
    print("MULTIBIT-0 SYSTEM: 8 bits × 64 rounds = 512 equations")
    print("=" * 80)

    N = 600  # messages
    n_eqs = 8 * 64  # 512 equations

    # Build Jacobian: for each of 512 message bits,
    # which of the 512 bit-0 values change?
    rng = random.Random(42)
    M_base = [rng.randint(0, MASK) for _ in range(16)]
    states_base, _ = sha256_round_trace(M_base)

    base_bits = []
    for r in range(64):
        for reg in range(8):
            base_bits.append(states_base[r+1][reg] & 1)  # bit 0 of reg at round r+1

    # Jacobian: 512 input bits → 512 output bits
    J = []
    for word in range(16):
        for bit in range(32):
            M_flip = list(M_base)
            M_flip[word] ^= (1 << bit)
            states_flip, _ = sha256_round_trace(M_flip)
            flip_bits = []
            for r in range(64):
                for reg in range(8):
                    flip_bits.append(states_flip[r+1][reg] & 1)
            J.append([base_bits[i] ^ flip_bits[i] for i in range(512)])

    # GF(2) rank
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

    rank = gf2_rank(J, 512, 512)
    print(f"  Jacobian of bit-0 system: rank = {rank}/512")
    print(f"  Kernel dimension: {512 - rank}")

    if rank < 512:
        print(f"  *** {512 - rank} KERNEL DIRECTIONS in bit-0 system! ***")
        print(f"  These are message perturbations that DON'T CHANGE any bit 0")
        print(f"  of any register at any round.")
    else:
        print(f"  Full rank: every message bit affects at least one bit-0 value.")

    # Also test: Jacobian of ONLY bit 0 of register a at all 64 rounds = 64 equations
    J_a0 = []
    for word in range(16):
        for bit in range(32):
            M_flip = list(M_base)
            M_flip[word] ^= (1 << bit)
            states_flip, _ = sha256_round_trace(M_flip)
            row = []
            for r in range(64):
                row.append((states_base[r+1][0] & 1) ^ (states_flip[r+1][0] & 1))
            J_a0.append(row)

    rank_a0 = gf2_rank(J_a0, 512, 64)
    print(f"\n  Jacobian of bit-0 of a only: rank = {rank_a0}/64")
    print(f"  Kernel: {512 - rank_a0} directions don't affect bit-0 of a at any round")


if __name__ == "__main__":
    experiment_bit0_composition()
    experiment_multibit0_system()
