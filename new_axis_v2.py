"""
NEW AXIS v2: GPK as PARTIAL INVERSION tool.

Key insight: (g, p) of the addition T1+T2 tells us about T1 and T2 SEPARATELY:
  - g[k]=1 → T1[k]=1 AND T2[k]=1 (both known)
  - kill[k]=1 → T1[k]=0 AND T2[k]=0 (both known)
  - p[k]=1 → T1[k]≠T2[k] (ambiguous: which is 0, which is 1?)

So from (g, p) we KNOW ~16 bits of T1,T2 exactly (g and k positions)
and have ~16 AMBIGUOUS bits (p positions).

If we track this across rounds: each round resolves some ambiguity
from previous rounds (because T1, T2 are functions of the state).

Question: Does the ambiguity SHRINK as we use more rounds of GPK?
If yes → GPK provides a path to inversion.
"""

import random
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, schedule
)

from new_axis import ShaElement, sha_element_round


def experiment_gpk_inversion_info():
    """
    For each round, how many bits of T1 and T2 does GPK reveal?

    At each addition T1+T2:
      g = T1 AND T2 → g[k]=1 means T1[k]=T2[k]=1
      p = T1 XOR T2 → p[k]=1 means T1[k]≠T2[k]
      k = ~g & ~p   → k[k]=1 means T1[k]=T2[k]=0

    Known bits of T1: g-positions (=1) + k-positions (=0) = 32 - HW(p) bits
    Ambiguous bits of T1: HW(p) bits (we know they differ but not which is which)

    For T2: same (symmetric).
    """
    print("=" * 80)
    print("GPK INVERSION: How many bits does GPK reveal per round?")
    print("=" * 80)

    N = 500
    known_per_round = []
    ambig_per_round = []

    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        W = schedule(M)

        state = [ShaElement(h) for h in H0]
        for r in range(64):
            state = sha_element_round(state, r, W[r])
            a_se = state[0]

            # Known bits = 32 - HW(p) = HW(g) + HW(kill)
            known = 32 - bin(a_se.p).count('1')
            ambig = bin(a_se.p).count('1')

            known_per_round.append(known)
            ambig_per_round.append(ambig)

    avg_known = sum(known_per_round) / len(known_per_round)
    avg_ambig = sum(ambig_per_round) / len(ambig_per_round)

    print(f"\n  Per round (avg over {N} messages × 64 rounds):")
    print(f"    Known T1/T2 bits:     {avg_known:.1f}/32")
    print(f"    Ambiguous T1/T2 bits: {avg_ambig:.1f}/32")
    print(f"    Information revealed:  {avg_known:.1f} bits per round")
    print(f"    Over 64 rounds:       {avg_known * 64:.0f} bits total")
    print(f"    Ambiguity over 64 rounds: {avg_ambig * 64:.0f} bits")
    print(f"    Message bits: 512")
    print(f"    Ratio ambiguity/message: {avg_ambig * 64 / 512:.2f}")


def experiment_propagate_resolution():
    """
    At propagate positions: T1[k] ≠ T2[k], but we don't know which is 1.

    However: T2 = Sig0(a) + Maj(a,b,c) depends only on a-branch.
    If we KNOW the a-branch values, we can COMPUTE T2 exactly.
    Then: T1 = a_new - T2 (subtraction), giving T1 exactly.

    But we don't know a-branch values (that's what we're trying to find).

    ALTERNATIVE: At the propagate position k:
      T1[k] XOR T2[k] = 1 (we know this)
      a_new[k] = T1[k] XOR T2[k] XOR carry[k] = 1 XOR carry[k]

    So at propagate positions: a_new[k] = 1 XOR carry[k]
    And carry[k] is determined by bits 0..k-1.

    This means: at propagate positions, a_new[k] is DETERMINED by the carry chain!
    There is NO ambiguity in a_new — it's fully determined.

    The ambiguity is in DECOMPOSING a_new into T1 and T2, not in a_new itself.
    But if we're going FORWARD (computing SHA), a_new is what we need.
    The ambiguity matters only going BACKWARD (inverting).

    For BACKWARD direction: at propagate positions, we know T1 XOR T2 = 1,
    but not T1 and T2 individually. This is 1 bit of ambiguity per propagate position.
    """
    print("\n" + "=" * 80)
    print("PROPAGATE RESOLUTION: What does ambiguity mean for inversion?")
    print("=" * 80)

    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]
    W = schedule(M)

    states, _ = sha256_round_trace(M)

    # For one round: verify that at propagate positions, T1 XOR T2 = 1
    r = 20
    a, b, c, d, e, f, g, h = states[r]
    T1 = (h + Sig1(e) + Ch(e, f, g) + K[r] + W[r]) & MASK
    T2 = (Sig0(a) + Maj(a, b, c)) & MASK
    a_new = states[r+1][0]

    p = T1 ^ T2  # propagate mask
    g_mask = T1 & T2  # generate mask

    # At propagate positions: T1 XOR T2 = 1 (by definition)
    # At generate positions: T1 = T2 = 1
    # At kill positions: T1 = T2 = 0

    print(f"  Round {r}: T1=0x{T1:08x}, T2=0x{T2:08x}")
    print(f"  Generate (both=1): {bin(g_mask).count('1')} positions")
    print(f"  Propagate (differ): {bin(p).count('1')} positions")
    print(f"  Kill (both=0): {bin((~T1 & ~T2) & MASK).count('1')} positions")

    # For backward step: given a_new, g, p, what do we know about state[r]?
    # We know: T1 + T2 = a_new (mod 2^32)
    # We know: T1 XOR T2 = p
    # We know: T1 AND T2 = g
    # From these three: T1 and T2 are FULLY DETERMINED!

    # T1 + T2 = a_new
    # T1 XOR T2 = p
    # Therefore: 2 * (T1 AND T2) = (T1 + T2) - (T1 XOR T2) = a_new - p
    # So: T1 AND T2 = (a_new - p) / 2 (if even)

    # But this is mod 2^32 arithmetic. Let's verify:
    t1_and_t2_calc = ((a_new - p) & MASK)
    # In mod 2^32: (a_new - p) should equal 2*g + carries
    # Actually: a + b = (a XOR b) + 2*(a AND b) is only true without wrapping...
    # In mod 2^32: a + b = (a XOR b) + 2*(a AND b) holds bit by bit with carry:
    #   (a+b)[k] = (a XOR b)[k] XOR carry[k]
    #   carry[k+1] = MAJ(a[k], b[k], carry[k])
    # So: a + b - (a XOR b) = carry_contribution

    carry_contrib = (a_new ^ p) & MASK  # This should equal the carry correction
    # a_new = p XOR carry_vector → carry_vector = a_new XOR p
    carry_vector = a_new ^ p

    print(f"\n  From a_new and p, recover carry vector:")
    print(f"    carry_vector = a_new XOR p = 0x{carry_vector:08x}")
    print(f"    HW(carry) = {bin(carry_vector).count('1')}")

    # From carry_vector and p, we can recover g:
    # At each bit k: if carry_vector[k]=1 and p[k]=1 → carry propagated through
    #                if carry_vector[k]=1 and p[k]=0 → carry was generated (g[k-1]=1) or killed
    # Actually this is more complex. But the point is:

    # Given (a_new, g, p) we can FULLY recover T1 and T2:
    # T1 = g | (p & mask) where mask picks which of T1,T2 has the 1 at propagate positions
    # T2 = g | (p & ~mask)
    # And mask is determined by the carry chain.

    # Let me just verify: can we recover T1 from (a_new, p, g)?
    # T1 = a_new - T2. T2 = Sig0(a) + Maj(a,b,c). We need a,b,c.
    # Going backward: b = a[r+1] (current a, which becomes next b), etc.
    # So we know b,c,d from the NEXT state: b = state[r+1][1] = a[r]... wait.
    # state[r+1] = [a_new, a, b, c, e_new, e, f, g_val]
    # So state[r+1][1] = a = state[r][0]. We KNOW this.

    # Therefore: T2 at round r = Sig0(state[r+1][1]) + Maj(state[r+1][1], state[r+1][2], state[r+1][3])
    a_prev = states[r+1][1]  # = states[r][0]
    b_prev = states[r+1][2]  # = states[r][1]
    c_prev = states[r+1][3]  # = states[r][2]
    T2_recovered = (Sig0(a_prev) + Maj(a_prev, b_prev, c_prev)) & MASK
    T1_recovered = (a_new - T2_recovered) & MASK

    print(f"\n  Recovery of T1, T2 from next state:")
    print(f"    T2_recovered = 0x{T2_recovered:08x} (actual = 0x{T2:08x}) {'OK' if T2_recovered == T2 else 'FAIL'}")
    print(f"    T1_recovered = 0x{T1_recovered:08x} (actual = 0x{T1:08x}) {'OK' if T1_recovered == T1 else 'FAIL'}")

    print(f"\n  CONCLUSION: SHA-256 is trivially invertible given full state.")
    print(f"  GPK doesn't add inversion power — the STATE already inverts.")
    print(f"  The problem is: we don't HAVE the full state. We have only the HASH.")


def experiment_gpk_trajectory_structure():
    """
    THE REAL QUESTION: Is the GPK TRAJECTORY more structured than the VALUE trajectory?

    Define GPK-trajectory: (g[0], p[0], g[1], p[1], ..., g[63], p[63])
    This is 64 × 64 = 4096 bits.

    Value-trajectory: (a[1], a[2], ..., a[64]) = 64 × 32 = 2048 bits.

    Both are determined by M (512 bits). So both live in a 512-dimensional
    subspace of their respective spaces.

    But the GPK-trajectory might have LOWER effective dimension because
    g and p at each round are more constrained than v.

    Specifically: g[k] + p[k] + kill[k] = 32 always (they partition the bits).
    And g = T1 AND T2, p = T1 XOR T2 — both derived from T1, T2.
    So (g, p) at round r has only as much entropy as (T1, T2) at round r.
    And T1, T2 have 64 bits total — but they're computed from 256-bit state.

    The question: across ALL 64 rounds, does the GPK trajectory have
    redundancy (lower rank) that the value trajectory doesn't?
    """
    print("\n" + "=" * 80)
    print("GPK TRAJECTORY: Is it more structured than value trajectory?")
    print("=" * 80)

    N = 600
    R = 64

    # Collect: for each message, the flattened GPK trajectory of a-register
    # (g[r] for r=0..63) = 64 × 32 = 2048 bits

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

    # Value trajectory: just a[r] for r=1..64 → 64×32 = 2048 bits
    v_rows = []
    g_rows = []
    p_rows = []

    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        W = schedule(M)

        state = [ShaElement(h) for h in H0]
        v_row = []
        g_row = []
        p_row = []

        for r in range(R):
            state = sha_element_round(state, r, W[r])
            a_se = state[0]

            for k in range(32):
                v_row.append((a_se.v >> k) & 1)
                g_row.append((a_se.g >> k) & 1)
                p_row.append((a_se.p >> k) & 1)

        v_rows.append(v_row)
        g_rows.append(g_row)
        p_rows.append(p_row)

    ncols = 32 * R  # 2048

    v_rank = gf2_rank(v_rows, N, ncols)
    g_rank = gf2_rank(g_rows, N, ncols)
    p_rank = gf2_rank(p_rows, N, ncols)

    print(f"  Trajectory rank over GF(2) (N={N} messages, {ncols} bits):")
    print(f"    Value trajectory rank: {v_rank}/{ncols}")
    print(f"    Generate trajectory rank: {g_rank}/{ncols}")
    print(f"    Propagate trajectory rank: {p_rank}/{ncols}")
    print(f"")
    print(f"    Both capped by N={N} and message bits=512.")

    if g_rank < v_rank:
        print(f"    *** g-trajectory has LOWER rank ({g_rank} < {v_rank})! ***")
        print(f"    → Generate pattern carries LESS info than values")
        print(f"    → It's a LOSSY projection of the state")
    elif g_rank > v_rank:
        print(f"    *** g-trajectory has HIGHER rank! Impossible if both determined by M. ***")
    else:
        print(f"    Same rank — both saturated by sample size or message dimension.")


if __name__ == "__main__":
    experiment_gpk_inversion_info()
    experiment_propagate_resolution()
    experiment_gpk_trajectory_structure()
