"""
WHY 127? SHA-256 reveals exactly 127 new constraints per bit layer.

512 = 4 × 127 + 4. Why not 128?

Hypotheses:
H1: The створочное число e[r] = a[r] + a[r-4] - T2[r-1] creates
    1 dependency between a and e per round → 1 less independent bit.
    64 rounds but rank of a alone = 64 and e alone = 64.
    Combined: 128 - 1 = 127? Test this.

H2: The schedule creates 1 linear dependency in bit-0 space.
    Test: rank for JUST rounds 0-15 (no schedule constraint).

H3: It's a property of the rotation constants.
    Test: change rotations, see if 127 changes.

H4: 127 = rank of the GF(2) schedule matrix restricted to bit 0.
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


def test_H1_ae_dependency():
    """Is 127 = 128 - 1 because e depends on a?"""
    print("=" * 80)
    print("H1: Is 127 = 128 - 1 because e depends on a?")
    print("=" * 80)

    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]
    states, _ = sha256_round_trace(M)

    # Rank of bit-0 of a ONLY (64 equations)
    base_a = [(states[r+1][0] & 1) for r in range(64)]
    J_a = []
    for word in range(16):
        for bit in range(32):
            M2 = list(M); M2[word] ^= (1 << bit)
            s2, _ = sha256_round_trace(M2)
            J_a.append([(base_a[r] ^ (s2[r+1][0] & 1)) for r in range(64)])
    rank_a = gf2_rank(J_a, 512, 64)

    # Rank of bit-0 of e ONLY (64 equations)
    base_e = [(states[r+1][4] & 1) for r in range(64)]
    J_e = []
    for word in range(16):
        for bit in range(32):
            M2 = list(M); M2[word] ^= (1 << bit)
            s2, _ = sha256_round_trace(M2)
            J_e.append([(base_e[r] ^ (s2[r+1][4] & 1)) for r in range(64)])
    rank_e = gf2_rank(J_e, 512, 64)

    # Rank of bit-0 of a AND e together (128 equations)
    base_ae = base_a + base_e
    J_ae = []
    for word in range(16):
        for bit in range(32):
            M2 = list(M); M2[word] ^= (1 << bit)
            s2, _ = sha256_round_trace(M2)
            row = [(base_a[r] ^ (s2[r+1][0] & 1)) for r in range(64)]
            row += [(base_e[r] ^ (s2[r+1][4] & 1)) for r in range(64)]
            J_ae.append(row)
    rank_ae = gf2_rank(J_ae, 512, 128)

    print(f"  rank(a bit0 only):   {rank_a}/64")
    print(f"  rank(e bit0 only):   {rank_e}/64")
    print(f"  rank(a+e bit0):      {rank_ae}/128")
    print(f"  rank(a) + rank(e):   {rank_a + rank_e}")
    print(f"  Overlap:             {rank_a + rank_e - rank_ae}")
    print(f"")
    if rank_ae == 127:
        print(f"  127 = rank(a+e). The 1 missing = 1 dependency between a and e.")
        print(f"  This is the створочное число: e[r] = f(a[r..r-4]).")
    elif rank_ae == 128:
        print(f"  128 = full. No a-e dependency at bit 0 level.")
        print(f"  127 must come from register shifts (b,c,d,f,g,h).")


def test_H2_schedule():
    """Is 127 from schedule constraint?"""
    print("\n" + "=" * 80)
    print("H2: Does rank change with fewer rounds (before schedule kicks in)?")
    print("=" * 80)

    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]

    for R in [8, 15, 16, 32, 64]:
        states, _ = sha256_round_trace(M, rounds=R)

        base_bits = []
        for r in range(R):
            for reg in range(8):
                base_bits.append((states[r+1][reg] & 1))

        n_eqs = len(base_bits)
        J = []
        for word in range(16):
            for bit in range(32):
                M2 = list(M); M2[word] ^= (1 << bit)
                s2, _ = sha256_round_trace(M2, rounds=R)
                row = []
                for r in range(R):
                    for reg in range(8):
                        row.append((base_bits[r*8+reg] ^ (s2[r+1][reg] & 1)))
                J.append(row)

        rank = gf2_rank(J, 512, n_eqs)
        per_round = rank / R if R > 0 else 0
        print(f"  R={R:>2}: rank={rank:>3}/512, n_eqs={n_eqs:>4}, rank/R={per_round:.2f}")


def test_universality_127():
    """Is 127 the same for all messages?"""
    print("\n" + "=" * 80)
    print("UNIVERSALITY: Is 127 constant across messages?")
    print("=" * 80)

    for seed in range(10):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        states, _ = sha256_round_trace(M)

        base_bits = []
        for r in range(64):
            for reg in range(8):
                base_bits.append((states[r+1][reg] & 1))

        J = []
        for word in range(16):
            for bit in range(32):
                M2 = list(M); M2[word] ^= (1 << bit)
                s2, _ = sha256_round_trace(M2)
                row = []
                for r in range(64):
                    for reg in range(8):
                        row.append((base_bits[r*8+reg] ^ (s2[r+1][reg] & 1)))
                J.append(row)

        rank = gf2_rank(J, 512, 512)
        print(f"  Seed {seed}: rank = {rank}")


if __name__ == "__main__":
    test_H1_ae_dependency()
    test_H2_schedule()
    test_universality_127()
