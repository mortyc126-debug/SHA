"""
TRAJECTORY SPACE: The 512-dimensional manifold where BTE theory lives.

Trajectory = {a[1], a[2], ..., a[64]} — the full a-sequence.
(e determined by створочное, other regs by shifts.)

Dimension = 512 (= message bits).
Embedding in: 64 × 32 = 2048 dimensional a-space.

Hash = projection: trajectory → state[64] → H.
Hash dimension = 256.
Preimage per hash: 2^{512-256} = 2^256 trajectories.

BTE GEOMETRY of trajectory space:
  Layer 0 (bit 0): 127-dim submanifold (T1)
  Layer 1 (bit 1): +127-dim, carry-connected to layer 0
  Layer 2: +127-dim
  Layer 3: +127-dim
  Layer 4: +4-dim (completes)

QUESTION: What does trajectory space look like GEOMETRICALLY?
Is it "flat" (linear subspace) or "curved" (nonlinear manifold)?

For LINEAR part (L): trajectory space = 512-dim linear subspace.
For CARRY part (Φ): adds nonlinear curvature.
Combined: 512-dim CURVED manifold in 2048-dim space.

The curvature = measure of how much carry deviates from linearity.
From T2: deficit ~0.022 bits/round = SMALL curvature.
"""

import random
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, schedule
)


def get_trajectory(M, R=64):
    """Get a-sequence trajectory: a[1], a[2], ..., a[R]."""
    states, _ = sha256_round_trace(M, rounds=R)
    return [states[r+1][0] for r in range(R)]


def trajectory_bits(traj, bit):
    """Extract one bit-layer from trajectory."""
    return [(t >> bit) & 1 for t in traj]


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


def experiment_trajectory_geometry():
    """
    Study the SHAPE of the trajectory manifold.

    Generate N trajectories (from random M).
    View them as points in 2048-dim space (64 × 32 bits of a-sequence).
    What is the GF(2) rank? (= linear dimension of the manifold)
    """
    print("=" * 80)
    print("TRAJECTORY GEOMETRY: Shape of the 512-dim manifold")
    print("=" * 80)

    N = 600
    R = 64

    # Full trajectory: 64 × 32 = 2048 bits
    traj_matrix = []
    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        traj = get_trajectory(M, R)
        row = []
        for r in range(R):
            for b in range(32):
                row.append((traj[r] >> b) & 1)
        traj_matrix.append(row)

    rank = gf2_rank(traj_matrix, N, 2048)
    print(f"  Full trajectory (2048 bits): rank = {rank}/{min(N, 2048)}")
    print(f"  Expected: 512 (= message dimension)")

    # Per bit-layer trajectory
    print(f"\n  Per bit-layer trajectory ranks:")
    for bit in [0, 1, 3, 7, 15, 31]:
        bit_matrix = []
        for seed in range(N):
            rng = random.Random(seed)
            M = [rng.randint(0, MASK) for _ in range(16)]
            traj = get_trajectory(M, R)
            row = [(traj[r] >> bit) & 1 for r in range(R)]
            bit_matrix.append(row)

        r = gf2_rank(bit_matrix, N, R)
        print(f"    bit {bit:>2}: rank = {r}/{R}")


def experiment_hash_projection():
    """
    Hash = projection of trajectory to last state.

    state[64] = (a[64], a[63], a[62], a[61], e[64], e[63], e[62], e[61])
    (from register shifts: b=a[-1], c=a[-2], d=a[-3])

    Hash = state[64] + IV (feedforward).

    Which DIRECTIONS in trajectory space map to 0 under hash projection?
    These = the "invisible" directions. Kernel of hash projection.
    """
    print("\n" + "=" * 80)
    print("HASH PROJECTION: Which trajectory directions are invisible?")
    print("=" * 80)

    rng = random.Random(42)
    M_base = [rng.randint(0, MASK) for _ in range(16)]

    # Hash = (a[64]+H0[0], a[63]+H0[1], a[62]+H0[2], a[61]+H0[3],
    #         e[64]+H0[4], e[63]+H0[5], e[62]+H0[6], e[61]+H0[7])
    # From створочное: e[r] = a[r] + a[r-4] - T2[r-1]
    # So hash = function of a[61], a[62], a[63], a[64] (and their e equivalents)
    # = function of a[57..64] (through e dependencies)

    # Hash depends on a[r] for r = 61, 62, 63, 64 (direct)
    # and on a[r-4] for e, which is r = 57, 58, 59, 60 (through створочное)

    # So: hash depends on a[57..64] = 8 words = 256 bits.
    # The other 56 words a[1..56] = INVISIBLE to hash!

    # But: a[57..64] DEPENDS on a[1..56] through the recurrence.
    # So they're not "free" — they're constrained by earlier trajectory.

    # The KERNEL of hash projection: message perturbations that
    # change a[1..56] but NOT a[57..64].

    # Measure: Jacobian of hash w.r.t. message
    base_hash = []
    states_base, _ = sha256_round_trace(M_base)
    for i in range(8):
        for b in range(32):
            base_hash.append(((states_base[64][i] + H0[i]) >> b) & 1)

    J_hash = []
    for word in range(16):
        for bit in range(32):
            M_flip = list(M_base)
            M_flip[word] ^= (1 << bit)
            states_flip, _ = sha256_round_trace(M_flip)
            row = []
            for i in range(8):
                for b in range(32):
                    row.append(base_hash[i*32+b] ^ (((states_flip[64][i] + H0[i]) >> b) & 1))
            J_hash.append(row)

    hash_rank = gf2_rank(J_hash, 512, 256)
    hash_kernel = 512 - hash_rank
    print(f"  Hash Jacobian: 512 → 256, rank = {hash_rank}, kernel = {hash_kernel}")

    # Now: Jacobian of a[r] for each round r
    print(f"\n  Per-round trajectory Jacobian rank (a[r] only, 32 bits):")
    for r in [1, 4, 8, 16, 32, 56, 60, 63, 64]:
        J_r = []
        for word in range(16):
            for bit in range(32):
                M_flip = list(M_base)
                M_flip[word] ^= (1 << bit)
                states_flip, _ = sha256_round_trace(M_flip, rounds=r)
                diff = states_base[r][0] ^ states_flip[r][0]
                row = [(diff >> b) & 1 for b in range(32)]
                J_r.append(row)

        r_rank = gf2_rank(J_r, 512, 32)
        print(f"    a[{r:>2}]: rank = {r_rank}/32")

    # KEY: how many MESSAGE bits affect the LAST 8 rounds (a[57..64])?
    print(f"\n  Message sensitivity to last 8 rounds (a[57..64]):")
    J_last8 = []
    for word in range(16):
        for bit in range(32):
            M_flip = list(M_base)
            M_flip[word] ^= (1 << bit)
            states_flip, _ = sha256_round_trace(M_flip)
            row = []
            for r in range(57, 65):
                diff = states_base[r][0] ^ states_flip[r][0]
                for b in range(32):
                    row.append((diff >> b) & 1)
            J_last8.append(row)

    rank_last8 = gf2_rank(J_last8, 512, 256)
    print(f"    a[57..64] Jacobian rank: {rank_last8}/256")
    print(f"    Kernel (invisible to last 8 rounds): {512 - rank_last8}")


def experiment_trajectory_vs_hash_info():
    """
    Quantify: how much MORE information does trajectory carry than hash?

    Trajectory = 2048 bits, rank 512.
    Hash = 256 bits, rank 256.
    Difference = 512 - 256 = 256 bits of trajectory info lost in hash.

    But: trajectory has STRUCTURE (layers, cocycle).
    Hash doesn't (random).

    Can we RECONSTRUCT trajectory structure from hash + BTE theory?
    I.e., given hash H and BTE constraints (T1-T6), what can we
    infer about the trajectory WITHOUT knowing M?
    """
    print("\n" + "=" * 80)
    print("TRAJECTORY vs HASH: Information gap")
    print("=" * 80)

    # From hash H: we know state[64] = (a[64], a[63], a[62], a[61], e[64]...)
    # From створочное: e[64] = a[64] + a[60] - T2[63]
    #   T2[63] = Sig0(a[63]) + Maj(a[63], a[62], a[61])
    # So: a[60] = e[64] - a[64] + T2[63] = known! (from hash)

    # Similarly: e[63] = a[63] + a[59] - T2[62] → a[59] known.
    # e[62] = a[62] + a[58] - T2[61] → a[58] known.
    # e[61] = a[61] + a[57] - T2[60] → a[57] known... if T2[60] known.
    # T2[60] = Sig0(a[60]) + Maj(a[60], a[59], a[58]) — all known!

    # So from hash: we can recover a[57..64] (8 words = 256 bits) via
    # backward створочное. But to go further (a[56], a[55]...):
    # a[56] = e[60] - ... but e[60] requires a[56] itself!

    # Actually: e[r] = a[r] + a[r-4] - T2[r-1]. Going backward from r=60:
    # We know a[60], a[59], a[58], a[57] (from above).
    # e[60] = state[64][5] = f[64] = e[63] (known from hash).
    # Wait: e[60] = g[63] = f[62] = e[61]. And e[61] is from hash.
    # Actually: state[64] = (a64, a63, a62, a61, e64, e63, e62, e61)
    # where a-branch shifts and e-branch shifts.

    # From hash we know: a[64], a[63], a[62], a[61], e[64], e[63], e[62], e[61].
    # From створочное: e[r] = a[r] + a[r-4] - T2[r-1].
    # e[64] → a[60] = e[64] - a[64] + T2[63]. T2[63] computable from a[63,62,61]. ✓
    # e[63] → a[59] = e[63] - a[63] + T2[62]. T2[62] from a[62,61,60]. ✓
    # e[62] → a[58]. T2[61] from a[61,60,59]. ✓
    # e[61] → a[57]. T2[60] from a[60,59,58]. ✓

    # Now: a[57,58,59,60] known. Can we get a[56]?
    # Need e[60]. e[60] = hash[5] from state[64]? No: state[64][5] = f[64] = e[63].
    # Actually: register at position 5 of state[64] is f = e[-1].
    # state[64] = (a[64], b[64]=a[63], c[64]=a[62], d[64]=a[61],
    #              e[64], f[64]=e[63], g[64]=e[62], h[64]=e[61])
    # So from hash: a[61..64] and e[61..64]. 8 values.
    # From створочное: recover a[57..60]. Now have a[57..64].

    # To get a[56]: use створочное at r=60: e[60] = a[60]+a[56]-T2[59].
    # → a[56] = e[60] - a[60] + T2[59].
    # T2[59] = Sig0(a[59]) + Maj(a[59],a[58],a[57]). KNOWN.
    # But e[60] = ? Not directly from hash. e[60] = g[63] (shift of 3).
    # state[63] is NOT known (only state[64]).

    # To get state[63]: need to INVERT one round from state[64].
    # Round inversion: given state[64] and W[63], recover state[63].
    # But W[63] = schedule(M) — UNKNOWN!

    # So: we can go back 4 rounds from hash (a[57..64]) using створочное,
    # but further requires W[r] which requires M.

    print(f"  From hash (8 words): can recover a[57..64] via створочное.")
    print(f"  = 4 extra rounds back. Total: 8 a-values from hash.")
    print(f"  Further back: requires W[r] = requires M. Blocked.")
    print(f"")
    print(f"  TRAJECTORY INFO AVAILABLE WITHOUT M:")
    print(f"    From hash: a[61..64], e[61..64] = 256 bits (direct)")
    print(f"    From створочное: a[57..60] = 128 bits (inferred)")
    print(f"    Total: a[57..64] = 256 bits trajectory known")
    print(f"    Unknown: a[1..56] = 1792 bits trajectory = REQUIRES M")
    print(f"")
    print(f"  This means: hash reveals 256/2048 = 12.5% of trajectory.")
    print(f"  Створочное adds 128 bits → 384/2048 = 18.8%.")
    print(f"  The other 81.2% is HIDDEN behind the message.")

    # Verify: can we actually recover a[57..60]?
    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]
    states, _ = sha256_round_trace(M)

    # From state[64]:
    a64 = states[64][0]; a63 = states[64][1]; a62 = states[64][2]; a61 = states[64][3]
    e64 = states[64][4]; e63 = states[64][5]; e62 = states[64][6]; e61 = states[64][7]

    # Recover a[60]: e[64] = a[64] + a[60] - T2[63]
    T2_63 = (Sig0(a63) + Maj(a63, a62, a61)) & MASK
    a60_recovered = (e64 - a64 + T2_63) & MASK

    # Recover a[59]: e[63] = a[63] + a[59] - T2[62]
    T2_62 = (Sig0(a62) + Maj(a62, a61, a60_recovered)) & MASK
    a59_recovered = (e63 - a63 + T2_62) & MASK

    # Verify
    a60_actual = states[60][0]
    a59_actual = states[59][0]

    print(f"\n  VERIFICATION:")
    print(f"    a[60] recovered: 0x{a60_recovered:08x}, actual: 0x{a60_actual:08x}, match: {a60_recovered == a60_actual}")
    print(f"    a[59] recovered: 0x{a59_recovered:08x}, actual: 0x{a59_actual:08x}, match: {a59_recovered == a59_actual}")


if __name__ == "__main__":
    experiment_trajectory_geometry()
    experiment_hash_projection()
    experiment_trajectory_vs_hash_info()
