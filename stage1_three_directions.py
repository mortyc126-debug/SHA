"""
Stage 1, Part 6: Three invented directions — tested concretely.

Direction A: Geometry of the preimage set S_H
Direction B: Deformation of preimages (new analog of differentials)
Direction C: Constraint interaction across rounds

We test on REDUCED SHA-256 (fewer rounds) where we can actually
compute preimage sets and examine their structure.
"""

import random
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, schedule
)


def sha256_reduced(msg16, rounds):
    """Compute SHA-256 reduced to `rounds` rounds, return hash (state + IV)."""
    W = schedule(msg16)
    a, b, c, d, e, f, g, h = H0
    for r in range(rounds):
        T1 = (h + Sig1(e) + Ch(e, f, g) + K[r] + W[r]) & MASK
        T2 = (Sig0(a) + Maj(a, b, c)) & MASK
        h, g, f = g, f, e
        e = (d + T1) & MASK
        d, c, b = c, b, a
        a = (T1 + T2) & MASK
    # Feedforward
    return [
        (a + H0[0]) & MASK, (b + H0[1]) & MASK,
        (c + H0[2]) & MASK, (d + H0[3]) & MASK,
        (e + H0[4]) & MASK, (f + H0[5]) & MASK,
        (g + H0[6]) & MASK, (h + H0[7]) & MASK,
    ]


def direction_A_preimage_geometry():
    """
    DIRECTION A: Geometry of the preimage set.

    For reduced SHA-256 (few rounds), fix most of the message,
    vary 2 words (W[0], W[1]), and map out which (W0, W1) pairs
    give the same hash value for specific output bits.

    This creates a 2D "landscape" where we can see the geometry
    of level sets (preimages of specific hash values).
    """
    print("=" * 80)
    print("DIRECTION A: Geometry of preimage set")
    print("=" * 80)

    rng = random.Random(42)
    base_msg = [rng.randint(0, MASK) for _ in range(16)]

    for R in [2, 4, 8, 16]:
        # Fix W[2..15], vary W[0] and W[1]
        # For each (W0, W1), compute hash bit 0 of register a
        # Map: (W0, W1) → hash_a[bit 0]

        N = 256  # grid size per dimension
        # We'll sample N×N points

        target_hash = None
        collision_count = 0
        hash_to_inputs = {}

        for i in range(N):
            for j in range(N):
                msg = list(base_msg)
                msg[0] = (base_msg[0] + i * 37) & MASK  # Spread out
                msg[1] = (base_msg[1] + j * 73) & MASK

                h = sha256_reduced(msg, R)
                # Take just the low byte of a-register as "hash"
                key = h[0] & 0xFF  # 8-bit hash

                if key not in hash_to_inputs:
                    hash_to_inputs[key] = []
                hash_to_inputs[key].append((i, j))

        # Analyze the distribution
        sizes = [len(v) for v in hash_to_inputs.values()]
        n_keys = len(hash_to_inputs)
        expected_per_key = N * N / 256  # 256 possible 8-bit values

        # Are preimage sets clustered or scattered?
        # Pick the largest preimage set and check if points are "close"
        largest_key = max(hash_to_inputs, key=lambda k: len(hash_to_inputs[k]))
        largest_set = hash_to_inputs[largest_key]

        # Compute average pairwise distance in the (i,j) grid
        if len(largest_set) > 1:
            total_dist = 0
            count = 0
            for a_idx in range(min(50, len(largest_set))):
                for b_idx in range(a_idx + 1, min(50, len(largest_set))):
                    p1 = largest_set[a_idx]
                    p2 = largest_set[b_idx]
                    dist = abs(p1[0] - p2[0]) + abs(p1[1] - p2[1])
                    total_dist += dist
                    count += 1
            avg_dist = total_dist / count if count > 0 else 0
            expected_random_dist = N * 2 / 3  # expected L1 distance for uniform

            print(f"\n  R={R} rounds, {N}×{N} grid, 8-bit hash:")
            print(f"    Distinct hash values: {n_keys}/256")
            print(f"    Largest preimage set: {len(largest_set)} points (expected {expected_per_key:.0f})")
            print(f"    Avg pairwise L1 distance in largest set: {avg_dist:.1f} (random: {expected_random_dist:.0f})")
            print(f"    → {'CLUSTERED' if avg_dist < expected_random_dist * 0.7 else 'SCATTERED (no geometry)'}")


def direction_B_deformation():
    """
    DIRECTION B: How does the preimage set deform when the target changes?

    For reduced SHA-256, find a message M with hash H.
    Then change H by 1 bit → H'. Find a message M' with hash H'.
    Question: How far is M' from M?

    If preimage sets are "close" when hashes are close → exploitable.
    If preimage sets are far apart → no structure.

    We can't find exact preimages for full SHA-256, but we can
    measure the REVERSE: take two close messages, see how close
    their hashes are across many rounds.

    More precisely: for M and M' = M + small perturbation,
    at what RATE does the hash difference grow vs rounds?
    We already measured this with differentials, but now we ask
    a different question: is the GROWTH RATE itself structured?
    """
    print("\n" + "=" * 80)
    print("DIRECTION B: Deformation — growth rate of hash difference")
    print("=" * 80)

    N = 500
    # For each pair (M, M XOR 1bit), compute hash at R rounds
    # Measure: Hamming distance of hashes at each R

    # We already know HW → 128 by round 8.
    # NEW question: Is the APPROACH to 128 the same for all bit positions?
    # I.e., does flipping bit 0 vs bit 15 vs bit 31 give a different CURVE?

    print(f"\n  Growth curves for different input bit positions (N={N}):")
    print(f"  (Average HW of hash difference at each round count)")
    print(f"\n  {'R':>3} | {'bit 0':>6} {'bit 7':>6} {'bit 15':>6} {'bit 24':>6} {'bit 31':>6} | {'W[0]b0':>7} {'W[8]b0':>7} {'W[15]b0':>7}")
    print(f"  " + "-" * 75)

    for R in [1, 2, 3, 4, 5, 6, 8, 12, 16, 32, 64]:
        results = {}

        for flip_desc, flip_word, flip_bit in [
            ('bit 0', 0, 0), ('bit 7', 0, 7), ('bit 15', 0, 15),
            ('bit 24', 0, 24), ('bit 31', 0, 31),
            ('W[0]b0', 0, 0), ('W[8]b0', 8, 0), ('W[15]b0', 15, 0)
        ]:
            total_hw = 0
            for seed in range(N):
                rng = random.Random(seed)
                M = [rng.randint(0, MASK) for _ in range(16)]
                M2 = list(M)
                M2[flip_word] ^= (1 << flip_bit)

                h1 = sha256_reduced(M, R)
                h2 = sha256_reduced(M2, R)
                hw = sum(bin(a ^ b).count('1') for a, b in zip(h1, h2))
                total_hw += hw

            results[flip_desc] = total_hw / N

        r = results
        print(f"  {R:>3} | {r['bit 0']:>6.1f} {r['bit 7']:>6.1f} {r['bit 15']:>6.1f} "
              f"{r['bit 24']:>6.1f} {r['bit 31']:>6.1f} | "
              f"{r['W[0]b0']:>7.1f} {r['W[8]b0']:>7.1f} {r['W[15]b0']:>7.1f}")

    print(f"\n  If all curves converge to 128 at the same rate → isotropic")
    print(f"  If different rates → anisotropy exists → potential structure")


def direction_C_constraints():
    """
    DIRECTION C: Constraint interaction.

    Each round of SHA-256 can be viewed as adding a CONSTRAINT
    on the message M. Specifically:

    Round r produces state[r+1] from state[r] and W[r].
    For a fixed target hash H, round 64 constrains state[63].
    Round 63 constrains state[62]. Etc.

    Going BACKWARD from H, each round constrains the state
    one step further back. After 64 backward steps, we have
    constraints on the initial state = IV + message.

    Key question: Are the 64 constraints INDEPENDENT?
    If yes: each constraint removes ~4 bits of freedom → total = 256 bits → birthday.
    If no: there's REDUNDANCY → fewer effective constraints → cheaper search.

    We can measure this by checking: does knowing the output of round r
    give information about the output of round r+k, BEYOND what
    the intermediate states provide?

    Concretely: for reduced R rounds, compute the RANK of the
    "constraint matrix" relating input bits to output bits.
    """
    print("\n" + "=" * 80)
    print("DIRECTION C: Constraint rank — independence of round constraints")
    print("=" * 80)

    # For different numbers of rounds, compute the effective dimensionality
    # of the hash function viewed as a constraint on message bits.
    #
    # Method: Sample N random messages. For each, compute hash (256 bits).
    # Form N×256 binary matrix (GF(2)). Compute rank.
    # If rank < 256 → constraints are dependent → redundancy exists.

    def hash_bits(msg, R):
        """Return 256-bit hash as list of bits."""
        h = sha256_reduced(msg, R)
        bits = []
        for word in h:
            for b in range(32):
                bits.append((word >> b) & 1)
        return bits

    N = 300

    for R in [1, 2, 4, 8, 16, 32, 64]:
        # Build Jacobian-like matrix:
        # Each row = how hash changes when one input bit flips
        # This is the "derivative" matrix over GF(2)
        rng = random.Random(42)
        M = [rng.randint(0, MASK) for _ in range(16)]

        base_hash = hash_bits(M, R)

        # Flip each of 512 input bits, record which of 256 hash bits change
        matrix = []
        for word in range(16):
            for bit in range(32):
                M_flip = list(M)
                M_flip[word] ^= (1 << bit)
                flip_hash = hash_bits(M_flip, R)

                row = [base_hash[i] ^ flip_hash[i] for i in range(256)]
                matrix.append(row)

        # Compute GF(2) rank of this 512×256 matrix
        rank = gf2_rank_rect(matrix, 512, 256)

        print(f"  R={R:>2} rounds: Jacobian rank = {rank}/256")
        if rank < 256:
            print(f"         → {256 - rank} DEPENDENT constraints! Redundancy exists.")
        else:
            print(f"         → All constraints independent (full rank)")


def gf2_rank_rect(matrix, nrows, ncols):
    """Compute rank of binary matrix over GF(2)."""
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
    direction_A_preimage_geometry()
    direction_B_deformation()
    direction_C_constraints()
