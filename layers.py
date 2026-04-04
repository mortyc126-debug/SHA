"""
LAYER ANATOMY: What is each of the 4 layers of 127?

Layer 0 (bit 0): rank 127. Carry-free. Pure GF(2).
Layer 1 (bit 1): +127 new. Has carry from bit 0.
Layer 2 (bit 2): +127 new. Has carry from bits 0-1.
Layer 3 (bit 3): +127 new. Has carry from bits 0-2.
Layer 4 (bit 4): +4 remaining. Full rank 512.

Questions:
1. Are the 127 new bits in layer 1 INDEPENDENT of layer 0?
   Or do they share structure?
2. What are the 4 remaining bits? Why exactly 4?
3. Does the 127 per layer decompose further? 127 = ?
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


def experiment_layer_independence():
    """
    For EACH bit position independently, compute rank.
    Are they all 127? Or does it vary?
    """
    print("=" * 80)
    print("LAYER ANATOMY: Rank of each individual bit position")
    print("=" * 80)

    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]
    states, _ = sha256_round_trace(M)

    print(f"\n  Individual bit rank (8 regs × 64 rounds each):")
    ranks = []
    for b in range(32):
        base = []
        for r in range(64):
            for reg in range(8):
                base.append((states[r+1][reg] >> b) & 1)

        J = []
        for word in range(16):
            for bit in range(32):
                M2 = list(M); M2[word] ^= (1 << bit)
                s2, _ = sha256_round_trace(M2)
                row = []
                for r in range(64):
                    for reg in range(8):
                        row.append((base[r*8+reg] ^ (s2[r+1][reg] >> b) & 1))
                J.append(row)

        rank = gf2_rank(J, 512, 512)
        ranks.append(rank)

        if b < 8 or b >= 28 or rank != 127:
            print(f"    bit {b:>2}: rank = {rank}")

    # Pattern analysis
    print(f"\n  Rank distribution:")
    from collections import Counter
    dist = Counter(ranks)
    for rank_val, count in sorted(dist.items()):
        print(f"    rank={rank_val}: {count} bit positions")

    # Is EVERY bit position giving 127?
    if all(r == 127 for r in ranks):
        print(f"\n  ALL 32 bit positions give rank = 127!")
        print(f"  Each bit independently sees the same 127-dim subspace.")
    else:
        unique_ranks = set(ranks)
        print(f"\n  Ranks: {sorted(unique_ranks)}")


def experiment_cumulative_detail():
    """
    Finer-grained cumulative unfolding: measure rank for each bit
    added individually.
    """
    print("\n" + "=" * 80)
    print("CUMULATIVE DETAIL: Rank as each bit is added")
    print("=" * 80)

    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]
    states, _ = sha256_round_trace(M)

    prev_rank = 0
    for max_bit in range(8):  # First 8 bits in detail
        base = []
        for r in range(64):
            for reg in range(8):
                for b in range(max_bit + 1):
                    base.append((states[r+1][reg] >> b) & 1)

        n_eqs = len(base)
        J = []
        for word in range(16):
            for bit in range(32):
                M2 = list(M); M2[word] ^= (1 << bit)
                s2, _ = sha256_round_trace(M2)
                row = []
                for r in range(64):
                    for reg in range(8):
                        for b in range(max_bit + 1):
                            row.append(base[r*8*(max_bit+1)+reg*(max_bit+1)+b] ^
                                      ((s2[r+1][reg] >> b) & 1))
                J.append(row)

        rank = gf2_rank(J, 512, n_eqs)
        delta = rank - prev_rank
        print(f"  bits 0..{max_bit}: rank={rank:>3}, Δ=+{delta:>3}, kernel={512-rank:>3}")
        prev_rank = rank


def experiment_what_are_the_4():
    """
    At bits 0-3: rank = 508. At bits 0-4: rank = 512.
    The 4 missing bits at layer 3 → filled at layer 4.

    What ARE these 4 bits? Which message bits are they?
    Find the 4-dimensional kernel at bits 0-3 level.
    """
    print("\n" + "=" * 80)
    print("THE FINAL 4: What bits complete the rank at layer 4?")
    print("=" * 80)

    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]
    states, _ = sha256_round_trace(M)

    # Build Jacobian for bits 0-3
    base = []
    for r in range(64):
        for reg in range(8):
            for b in range(4):
                base.append((states[r+1][reg] >> b) & 1)

    n_eqs = len(base)
    J = []
    for word in range(16):
        for bit in range(32):
            M2 = list(M); M2[word] ^= (1 << bit)
            s2, _ = sha256_round_trace(M2)
            row = []
            for r in range(64):
                for reg in range(8):
                    for b in range(4):
                        row.append(base[r*32+reg*4+b] ^ ((s2[r+1][reg] >> b) & 1))
            J.append(row)

    # Find the null space (kernel) — these are the 4 directions
    # that DON'T affect bits 0-3 of ANY register at ANY round
    m = [list(row[:n_eqs]) for row in J[:512]]
    pivot_cols = []
    rank = 0
    for col in range(n_eqs):
        pivot = None
        for row in range(rank, 512):
            if m[row][col] == 1:
                pivot = row
                break
        if pivot is None:
            continue
        m[rank], m[pivot] = m[pivot], m[rank]
        pivot_cols.append(col)
        for row in range(512):
            if row != rank and m[row][col] == 1:
                for c in range(n_eqs):
                    m[row][c] ^= m[rank][c]
        rank += 1

    kernel_dim = 512 - rank
    print(f"  Rank of bits 0-3 system: {rank}/512")
    print(f"  Kernel dimension: {kernel_dim}")

    # Find kernel vectors
    pivot_set = set(pivot_cols)
    free_rows = [r for r in range(512) if r not in pivot_set and r < 512]

    # The "free rows" of the input correspond to message bits
    # that don't affect bits 0-3 of any register
    print(f"  Free input bit indices (first 10): {free_rows[:10]}")

    # Translate to (word, bit) positions
    print(f"  As (word, bit) positions:")
    for idx in free_rows[:kernel_dim + 4]:
        word = idx // 32
        bit = idx % 32
        print(f"    M[{word}] bit {bit}")


if __name__ == "__main__":
    experiment_layer_independence()
    experiment_cumulative_detail()
    experiment_what_are_the_4()
