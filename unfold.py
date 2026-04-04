"""
ITERATIVE BIT UNFOLDING: Build SHA-256 understanding bit by bit.

bit-0 system: rank=127, kernel=385 — sees 1/4 of message
bit-0+1 system: rank=?, kernel=? — should see more
...
bit-0..31 system: rank=256, kernel=256 — full SHA-256

HOW does the rank grow? Linearly? In steps? With plateaus?
The SHAPE of this growth curve is a fingerprint of SHA-256's structure.
"""

import random
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, schedule
)


def get_bit_system_jacobian(M_base, max_bit, rounds=64):
    """
    Compute Jacobian of the system: bits 0..max_bit of all 8 registers
    at all rounds.

    Output: n_eqs = 8 * rounds * (max_bit + 1) equations
    Input: 512 message bits
    Returns: rank of the Jacobian
    """
    states_base, _ = sha256_round_trace(M_base, rounds=rounds)

    # Base bit values
    base_bits = []
    for r in range(rounds):
        for reg in range(8):
            for b in range(max_bit + 1):
                base_bits.append((states_base[r + 1][reg] >> b) & 1)

    n_eqs = len(base_bits)

    # Jacobian: flip each of 512 input bits
    J = []
    for word in range(16):
        for bit in range(32):
            M_flip = list(M_base)
            M_flip[word] ^= (1 << bit)
            states_flip, _ = sha256_round_trace(M_flip, rounds=rounds)

            flip_bits = []
            for r in range(rounds):
                for reg in range(8):
                    for b in range(max_bit + 1):
                        flip_bits.append((states_flip[r + 1][reg] >> b) & 1)

            J.append([base_bits[i] ^ flip_bits[i] for i in range(n_eqs)])

    return J, n_eqs


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


def experiment_unfolding():
    """Build SHA-256 bit by bit and measure the rank growth."""
    print("=" * 80)
    print("ITERATIVE BIT UNFOLDING: How SHA-256 reveals itself bit by bit")
    print("=" * 80)

    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]

    print(f"\n  Bits included | Equations | Rank/512 | Kernel | New rank per bit")
    print(f"  " + "-" * 65)

    prev_rank = 0
    for max_bit in range(32):
        J, n_eqs = get_bit_system_jacobian(M, max_bit)
        rank = gf2_rank(J, 512, n_eqs)
        kernel = 512 - rank
        new_per_bit = rank - prev_rank

        print(f"  bits 0..{max_bit:>2}     | {n_eqs:>5}     | {rank:>3}/512 | {kernel:>3}    | +{new_per_bit}")

        prev_rank = rank

        if rank >= 512:
            print(f"\n  FULL RANK reached at bit {max_bit}!")
            break


def experiment_rank_per_register():
    """
    Same unfolding but tracking WHICH register contributes rank.

    At each bit level, how much rank comes from a vs e vs others?
    """
    print("\n" + "=" * 80)
    print("RANK PER REGISTER: Which register contributes most at each bit?")
    print("=" * 80)

    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]

    reg_names = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']

    # For a few key bit levels, measure rank from each register separately
    for max_bit in [0, 1, 3, 7, 15, 31]:
        states_base, _ = sha256_round_trace(M)

        ranks_per_reg = []
        for reg in range(8):
            base_bits = []
            for r in range(64):
                for b in range(max_bit + 1):
                    base_bits.append((states_base[r + 1][reg] >> b) & 1)

            n_eqs = len(base_bits)
            J = []
            for word in range(16):
                for bit in range(32):
                    M_flip = list(M)
                    M_flip[word] ^= (1 << bit)
                    states_flip, _ = sha256_round_trace(M_flip)
                    flip_bits = []
                    for r in range(64):
                        for b in range(max_bit + 1):
                            flip_bits.append((states_flip[r + 1][reg] >> b) & 1)
                    J.append([base_bits[i] ^ flip_bits[i] for i in range(n_eqs)])

            rank = gf2_rank(J, 512, n_eqs)
            ranks_per_reg.append(rank)

        reg_str = "  ".join(f"{reg_names[i]}:{ranks_per_reg[i]:>3}" for i in range(8))
        total = sum(ranks_per_reg)
        print(f"  bits 0..{max_bit:>2}: {reg_str}  (sum={total}, overlap={total - min(512, total)})")


if __name__ == "__main__":
    experiment_unfolding()
    experiment_rank_per_register()
