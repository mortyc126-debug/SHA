"""
SYNTHESIS: Reassemble ALL findings into one picture.

BRICKS we have:
  A. SHA-256 = F(a[r..r-7], K[r], W[r]) — one 8th-order recurrence
  B. F = XOR_PART ⊕ CARRY_PART
  C. XOR_PART = exact degree-2 GF(2) polynomial (KNOWN)
  D. CARRY_PART = ~16 bits, in ~8 segments of length ~3.6
  E. Segments rebuilt each round (don't survive rotation)
  F. Bit 0 = always carry-free (pure XOR at every addition)
  G. Rotations {2,6,11,13,22,25} generate Z/32 in 3 steps
  H. Two times: round r and bit k, coupled by rotation
  I. Every scalar observable → random by round 8
  J. Intermediate complexity: 48% carry

SYNTHESIS IDEA:

  OMEGA (from methodology) treats carry as VARIABLES and builds a
  GF(2) system. It works for R≤15, fails at R=16 because
  512 message bits = 512 carry constraints (α-kernel = 0).

  But brick D says: carries are NOT 32 independent constraints per addition.
  They're ~8 SEGMENTS. Within a segment, carry bits are CHAINED
  (each determined by the previous). The independent "choices" are
  just: does each segment generate a carry or not?

  QUESTION: If carries have ~8 independent degrees of freedom per addition
  (not 32), does OMEGA's α-kernel become nonzero at R=16?

  This would mean: OMEGA + segmentation = works beyond R=15.
"""

import random
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, schedule
)


def measure_effective_carry_dof():
    """
    For each addition in SHA-256, measure the EFFECTIVE degrees of freedom
    in the carry pattern.

    Method: For a fixed round, sample N random states.
    For each state, compute the carry pattern of T1+T2.
    The carry pattern is 32 bits. But how many of these bits are
    INDEPENDENT?

    If carries were fully independent: rank of the N×32 carry matrix = 32.
    If carries are segmented: rank < 32 (segments create dependencies).

    Actually, carry bits ARE deterministic given the operands.
    The question is: over random operands, what is the GF(2) rank
    of the carry pattern matrix?

    If rank = 32: carry patterns fill all of GF(2)^32 → no redundancy.
    If rank < 32: there are CONSTRAINTS on carry patterns → redundancy.
    """
    print("=" * 80)
    print("SYNTHESIS TEST 1: Effective carry DOF per addition")
    print("=" * 80)

    N = 500

    for r in [0, 8, 20, 40, 63]:
        carry_matrix = []

        for seed in range(N):
            rng = random.Random(seed)
            M = [rng.randint(0, MASK) for _ in range(16)]
            states, W = sha256_round_trace(M)
            a, b, c, d, e, f, g, h = states[r]

            T1 = (h + Sig1(e) + Ch(e, f, g) + K[r] + W[r]) & MASK
            T2 = (Sig0(a) + Maj(a, b, c)) & MASK

            # Carry pattern
            carry_bits = []
            c_bit = 0
            for k in range(32):
                carry_bits.append(c_bit)
                xk = (T1 >> k) & 1
                yk = (T2 >> k) & 1
                c_bit = (xk & yk) | (xk & c_bit) | (yk & c_bit)

            carry_matrix.append(carry_bits)

        # GF(2) rank
        rank = gf2_rank(carry_matrix, N, 32)
        print(f"  Round {r:>2}: carry GF(2) rank = {rank}/32 (from {N} samples)")

    print(f"\n  If rank < 32 → carry patterns are CONSTRAINED → potential for exploitation")
    print(f"  Bit 0 always = 0 → guaranteed rank ≤ 31")


def measure_all_carries_per_round():
    """
    In one round, there are 6 additions:
      T1: h + Sig1(e) → +Ch → +K → +W  (4 additions)
      T2: Sig0(a) + Maj(a,b,c) (1 addition)
      Final: T1 + T2 (1 addition)

    Total: 6 × 32 = 192 carry bits per round.
    But bit 0 of each = always 0 → 6 redundant → 186 max.

    Over R rounds: 186R carry bits.
    Message freedom: 512 bits.

    If effective carry DOF < 186 per round → OMEGA extends further.

    Specifically: at how many rounds does carry DOF = message freedom (512)?
    OMEGA says: R = 16 (512 = 32 × 16, using ~32 effective per round).
    With segmentation: maybe more rounds needed → OMEGA works further.
    """
    print("\n" + "=" * 80)
    print("SYNTHESIS TEST 2: Total carry DOF across rounds")
    print("=" * 80)

    N = 300

    # For each message, collect ALL carry patterns for ALL additions, ALL rounds
    # Then compute overall GF(2) rank

    for R in [4, 8, 12, 16, 20, 24]:
        all_carry_rows = []

        for seed in range(N):
            rng = random.Random(seed)
            M = [rng.randint(0, MASK) for _ in range(16)]
            states, W = sha256_round_trace(M, rounds=R)

            carry_row = []
            for r in range(R):
                a, b, c, d, e, f, g, h = states[r]

                # Collect carry from T1+T2 only (the dominant addition)
                T1 = (h + Sig1(e) + Ch(e, f, g) + K[r] + W[r]) & MASK
                T2 = (Sig0(a) + Maj(a, b, c)) & MASK

                c_bit = 0
                for k in range(32):
                    carry_row.append(c_bit)
                    xk = (T1 >> k) & 1
                    yk = (T2 >> k) & 1
                    c_bit = (xk & yk) | (xk & c_bit) | (yk & c_bit)

            all_carry_rows.append(carry_row)

        # Total carry bits per message = 32 × R
        n_cols = 32 * R
        rank = gf2_rank(all_carry_rows, N, n_cols)

        print(f"  R={R:>2}: {n_cols} carry bits, GF(2) rank = {rank}/{n_cols}")
        print(f"         Effective DOF = {rank}, message DOF = 512, α_potential = {512 - rank}")


def measure_xor_part_rank():
    """
    SYNTHESIS TEST 3: The XOR part of SHA-256.

    a[r+1] = XOR_PART ⊕ CARRY_PART

    XOR_PART is degree-2 over GF(2). If we could REMOVE the carry part,
    we'd have a pure degree-2 system. OMEGA does exactly this by treating
    carry as variables.

    But what is the GF(2) rank of the XOR_PART mapping M → H?
    If XOR_PART has rank < 256 → there's a kernel → preimages exist in XOR world.
    """
    print("\n" + "=" * 80)
    print("SYNTHESIS TEST 3: GF(2) rank of XOR-only SHA-256 (all + replaced by XOR)")
    print("=" * 80)

    def sha256_xor_only(msg16, rounds):
        """SHA-256 with all additions replaced by XOR."""
        W = list(msg16)
        for i in range(16, 64):
            s1 = rotr(W[i-2], 17) ^ rotr(W[i-2], 19) ^ (W[i-2] >> 10)
            s0 = rotr(W[i-15], 7) ^ rotr(W[i-15], 18) ^ (W[i-15] >> 3)
            W.append(s1 ^ W[i-7] ^ s0 ^ W[i-16])  # XOR instead of +

        a, b, c, d, e, f, g, h = H0
        for r in range(rounds):
            T1 = h ^ Sig1(e) ^ Ch(e, f, g) ^ K[r] ^ W[r]  # XOR instead of +
            T2 = Sig0(a) ^ Maj(a, b, c)  # XOR instead of +
            h, g, f = g, f, e
            e = d ^ T1  # XOR instead of +
            d, c, b = c, b, a
            a = T1 ^ T2  # XOR instead of +

        return [a ^ H0[0], b ^ H0[1], c ^ H0[2], d ^ H0[3],
                e ^ H0[4], f ^ H0[5], g ^ H0[6], h ^ H0[7]]

    # Compute Jacobian of XOR-only SHA over GF(2)
    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]

    for R in [4, 8, 16, 32, 64]:
        base_hash = sha256_xor_only(M, R)
        base_bits = []
        for word in base_hash:
            for b in range(32):
                base_bits.append((word >> b) & 1)

        # Flip each of 512 input bits
        matrix = []
        for word in range(16):
            for bit in range(32):
                M_flip = list(M)
                M_flip[word] ^= (1 << bit)
                flip_hash = sha256_xor_only(M_flip, R)
                row = []
                for i, word_h in enumerate(flip_hash):
                    for b in range(32):
                        row.append(((word_h >> b) & 1) ^ base_bits[i * 32 + b])
                matrix.append(row)

        rank = gf2_rank(matrix, 512, 256)
        null_dim = 512 - rank  # Kernel dimension = free directions

        print(f"  R={R:>2}: XOR-SHA Jacobian rank = {rank}/256, kernel dim = {null_dim}")


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


if __name__ == "__main__":
    measure_effective_carry_dof()
    measure_all_carries_per_round()
    measure_xor_part_rank()
