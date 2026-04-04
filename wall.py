"""
BUILDING THE WALL: SHA-256 as a SINGLE system of constraints.

Not round-by-round. All 64 rounds AT ONCE.

The system:
  - 512 unknowns: M[0..15] (message bits)
  - 256 knowns: H[0..7] (hash output)
  - 64 round functions connecting them

In BTE/★ language, this becomes:
  - 512 unknowns: M bits
  - Per round: 93.7% fixed linear constraints + 6.3% P-mask nonlinear
  - Schedule: W[16..63] = GF(2)-linear function of M

Total: 64 × 256 = 16384 intermediate "equations" (state bits at each round)
But most are redundant (state[r+1] is determined by state[r]).

The INDEPENDENT constraints are:
  - 256 final hash equations: state[64] + IV = H
  - 1536 schedule equations: W[16..63] = f(W[0..15])

That's 256 + 1536 = 1792 equations in 512 unknowns.
Massively over-determined → solutions exist only at specific M.

For COLLISION: find M, M' with same H.
  1792 × 2 = 3584 equations, 1024 unknowns (M, M'), 256 shared (same H).

But this is the STANDARD view. What's new in BTE?

NEW: separate the 93.7% (linear) from 6.3% (P-mask).
  Linear part: 93.7% of 16384 = ~15352 linear constraints
  P-mask part: 6.3% = ~1032 nonlinear constraints

If we could solve the linear part first (cheap) and then
handle the P-mask part (small)... but they're coupled.

UNLESS: we can find a basis where they DECOUPLE.
"""

import random
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, schedule
)


def experiment_skeleton_system():
    """
    The 93.7% SKELETON of SHA-256: extract the FIXED linear part.

    For each round, the Jacobian has ~61387 fixed entries (out of 65536).
    These form a FIXED 256×256 matrix over GF(2).

    The product of 64 such matrices = the "skeleton transfer matrix."
    This describes how M maps to H through the LINEAR (carry-free) path.

    If this skeleton has useful structure (rank < 256, specific kernel)
    → it reveals the LINEAR backbone of SHA-256.
    """
    print("=" * 80)
    print("SKELETON: The fixed linear backbone of SHA-256")
    print("=" * 80)

    # The "skeleton" is what SHA-256 would be if ALL P-masks were 0
    # (no carry propagation at all). This is equivalent to XOR-SHA.
    # We already know: XOR-SHA kernel = 256 dim (PHI result).

    # BUT: the skeleton is NOT just XOR-SHA.
    # It's XOR-SHA with specific Ch/Maj contributions.
    # Ch and Maj are degree-2 (nonlinear even without carries).

    # Let me separate differently:
    # TRULY LINEAR part = everything except Ch, Maj, AND carries.
    # = Sig0, Sig1 (rotations), register shifts, K additions, W additions
    # All of these are GF(2)-linear (XOR of known positions).

    # Ch(e,f,g) = ef + (1-e)g = ef + g + eg (degree 2)
    # Maj(a,b,c) = ab + ac + bc (degree 2)

    # So degree-2 part = Ch + Maj = 2 operations per round, each 32 bits.
    # Total degree-2: 64 × 2 × 32 = 4096 "quadratic bit equations" over 64 rounds.

    # But these are the SAME quadratic terms as in OMEGA.
    # OMEGA handles them by treating carry as variables → system stays degree 2.

    # What's OUR new contribution?
    # We know: 93.7% of the Jacobian is STRUCTURALLY FIXED.
    # The 6.3% variable part has a SPECIFIC form: P-mask carry chains.
    # P-mask carry chains = consecutive AND (propagation).

    # So SHA-256 = LINEAR(L) + QUADRATIC(Ch,Maj) + P-MASK_CARRY(sequential AND)

    # L: GF(2)-linear. Dimension: 256 per round. Full rank.
    # Q: GF(2)-quadratic. 64 bits per round (32 Ch + 32 Maj).
    # P: Carry chain. ~16 bits per addition × 6 additions = 96 bits per round.

    # Total: L(big) + Q(small) + P(medium)

    # The KEY: Q and P together = the 6.3%.
    # L alone = 93.7%.

    # In OMEGA: Q is handled (carry as variables makes it degree 2).
    # P is NOT handled (carry variables are treated as independent, but they're chained).

    # Our BTE insight: P has CHAIN STRUCTURE (MAJ sequential).
    # OMEGA ignores this structure → treats 32 carry variables as independent.
    # Reality: the 32 carry variables are CHAINED → only ~8 independent segments.

    # Could this reduce OMEGA's variable count?

    # OMEGA for R rounds: 32R carry variables.
    # With chain structure: ~8R independent segments × ~4 bits each.
    # But the segment boundaries CHANGE with state → can't fix globally.

    # Hmm. Let me compute something concrete instead of theorizing.

    # CONCRETE: For the SKELETON (no carries), what does SHA-256 look like
    # as a system of GF(2) equations?

    # SHA-256 round without carries (XOR-only) + Ch + Maj:
    # a_new = Sig0(a) XOR Maj(a,b,c) XOR h XOR Sig1(e) XOR Ch(e,f,g) XOR K XOR W
    # e_new = d XOR h XOR Sig1(e) XOR Ch(e,f,g) XOR K XOR W

    # These are EXACT (no approximation) for the XOR-only version.
    # Degree = 2 (from Ch and Maj).

    # For the REAL SHA-256: same equations PLUS carry corrections.
    # carry corrections = P-mask dependent.

    # So: real_a_new = xor_a_new XOR carry_correction_a
    #     real_e_new = xor_e_new XOR carry_correction_e

    # carry_correction = function of P-masks (themselves functions of state).

    # OMEGA treats carry_correction as VARIABLES.
    # We now know: carry_correction is determined by P-masks,
    # which have CHAIN structure.

    # Let me count: per round, how many INDEPENDENT carry correction bits?
    # There are 6 additions. Each has 32 carry bits.
    # But carry[0] = 0 always → 31 per addition × 6 = 186 per round.
    # Over 64 rounds: 186 × 64 = 11904 carry variables in OMEGA.
    # But most are determined by state (not free variables).

    # OMEGA Gaussian-eliminates to get effective carry variables.
    # For R=16: effective = 496 (we measured). α-kernel = 16 (trivial).

    # With CHAIN structure: each addition has ~8 segments.
    # At segment boundaries: carry = 0 (deterministic).
    # Within segments: carry propagates (each bit = AND of previous).
    # The "free" information per segment: just whether carry ENTERS the segment.
    # That's 1 bit per segment × 8 segments × 6 additions = 48 bits per round.
    # Over 64 rounds: 48 × 64 = 3072 "segment entry" bits.

    # But these 3072 bits are determined by state (not free).
    # After Gaussian elimination: ??? Need to measure.

    print(f"  SHA-256 decomposition:")
    print(f"    L (linear): Sig0, Sig1, register shifts, K+W = 93.7% of Jacobian")
    print(f"    Q (quadratic): Ch, Maj = degree 2, 64 bits per round")
    print(f"    P (carry): ~186 carry bits per round, CHAINED (not independent)")
    print(f"")
    print(f"  OMEGA: treats all carry bits as independent variables")
    print(f"    → 186 × R variables per R rounds")
    print(f"    → After elimination: 31R effective (we measured: 31×16 = 496)")
    print(f"")
    print(f"  BTE/★ insight: carries are CHAINED")
    print(f"    → ~8 segments per addition, 1 free bit per segment entry")
    print(f"    → 48 × R segment-entry bits per R rounds")
    print(f"    → Potential effective: 48R vs OMEGA's 31R")
    print(f"    → Wait, 48 > 31? That's WORSE, not better.")
    print(f"")
    print(f"  The issue: segment-entry bits are NOT the same as carry variables.")
    print(f"  OMEGA's 31R already accounts for bit-0 = 0.")
    print(f"  Segments add structure but don't REDUCE variable count.")

    # OK so the direct carry-variable approach doesn't help.
    # The chain structure doesn't give fewer variables — it gives
    # DIFFERENT (correlated) variables.

    # Let me try something else: work in P-MASK SPACE directly.

    # P-mask = T1 XOR T2 (32 bits per addition).
    # Over 64 rounds × 6 additions = 384 P-masks = 12288 bits.
    # These are determined by M (512 bits).

    # The mapping M → P-masks is 512 → 12288.
    # What's the RANK of this mapping?

    print(f"\n  Mapping M → P-masks: computing rank...")

    rng = random.Random(42)
    M_base = [rng.randint(0, MASK) for _ in range(16)]

    def get_all_pmasks(M):
        """Get all P-masks as a flat bit vector."""
        states, W = sha256_round_trace(M)
        bits = []
        for r in range(64):
            a, b, c, d, e, f, g, h = states[r]
            sig1e = Sig1(e)
            che = Ch(e, f, g)
            # Just do T1+T2 and d+T1 (the two output additions)
            T1 = (h + sig1e + che + K[r] + W[r]) & MASK
            T2 = (Sig0(a) + Maj(a, b, c)) & MASK
            p_a = T1 ^ T2  # P-mask for a_new
            p_e = d ^ T1    # P-mask for e_new
            for k in range(32):
                bits.append((p_a >> k) & 1)
            for k in range(32):
                bits.append((p_e >> k) & 1)
        return bits  # 64 × 2 × 32 = 4096 bits

    base_pmask_bits = get_all_pmasks(M_base)
    n_pmask_bits = len(base_pmask_bits)

    # Jacobian: 512 input bits → n_pmask_bits output bits
    J_pmask = []
    for word in range(16):
        for bit in range(32):
            M_flip = list(M_base)
            M_flip[word] ^= (1 << bit)
            flip_bits = get_all_pmasks(M_flip)
            J_pmask.append([base_pmask_bits[i] ^ flip_bits[i] for i in range(n_pmask_bits)])

    # Rank
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

    rank = gf2_rank(J_pmask, 512, n_pmask_bits)
    print(f"  M → P-masks Jacobian: {n_pmask_bits} P-mask bits, rank = {rank}/512")
    print(f"  P-mask kernel dimension: {512 - rank}")
    print(f"  (directions in M that don't change ANY P-mask)")


if __name__ == "__main__":
    experiment_skeleton_system()
