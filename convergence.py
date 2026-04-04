"""
CONVERGENCE: How two trajectories merge to create a collision.

Collision: state1[64] = state2[64] → a1[57..64] = a2[57..64].
But a1[1..56] ≠ a2[1..56] (different M).

The CONVERGENCE POINT: the earliest round where a1[r] = a2[r].
Does it exist? At which round? What forces it?
"""

import random
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, schedule
)


def experiment_convergence_profile():
    """
    For a collision (if we can find one on reduced rounds),
    trace WHEN the trajectories converge.

    For full SHA-256: can't find collisions.
    For reduced rounds (R≤16): Wang chain gives collision pairs.
    Use Wang-like approach to find pairs with matching final state.

    Actually: simpler. We know from створочное that collision ↔ a1[R-7..R] = a2[R-7..R].
    For REDUCED rounds (small R): search by birthday on a[R-3..R] (4 words = 128 bits).

    For R=8: hash = 256 bits. Birthday = 2^128. Too expensive.

    Alternative: use δW[0]=1 cascade from methodology.
    Wang chain gives δe[2..17]=0 for P=1. This means state difference
    is concentrated in a-registers only (δe=0).

    From methodology P-92: Wang chain verified.
    Let's trace the TRAJECTORY CONVERGENCE for Wang pairs.
    """
    print("=" * 80)
    print("CONVERGENCE: Where do collision trajectories merge?")
    print("=" * 80)

    # Instead of finding actual collision (hard), study the STRUCTURE:
    # For two messages differing in δM, what is δa[r] at each round?
    # δa[r] = a1[r] XOR a2[r]

    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]

    # Small perturbation
    M2 = list(M); M2[0] ^= 1

    states1, _ = sha256_round_trace(M)
    states2, _ = sha256_round_trace(M2)

    print(f"\n  δM = W[0] bit 0 flip. Trajectory divergence profile:")
    print(f"  {'r':>3} | {'HW(δa)':>7} {'HW(δe)':>7} | δa[r] first 8 hex")
    print(f"  " + "-" * 60)

    for r in range(65):
        da = states1[r][0] ^ states2[r][0]
        de = states1[r][4] ^ states2[r][4]
        hw_a = bin(da).count('1')
        hw_e = bin(de).count('1')

        if r <= 8 or r >= 57 or r % 16 == 0:
            print(f"  {r:>3} | {hw_a:>7} {hw_e:>7} | 0x{da:08x}")

    # For COLLISION: δa[57..64] = 0 AND δe[57..64] = 0.
    # This requires: all differences to "cancel out" by round 57.
    # How much cancellation is needed?
    total_diff_57 = sum(bin(states1[r][0] ^ states2[r][0]).count('1') for r in range(57, 65))
    total_diff_57 += sum(bin(states1[r][4] ^ states2[r][4]).count('1') for r in range(57, 65))
    print(f"\n  Total HW of δ(a,e) at rounds 57-64: {total_diff_57}")
    print(f"  For collision: this must be EXACTLY 0.")
    print(f"  Probability of random cancellation: 2^-{8*64} (impossibly small)")


def experiment_what_convergence_requires():
    """
    For collision: state1[64] = state2[64].
    This means: a1[64] = a2[64] AND a1[63] = a2[63] AND ... a1[61] = a2[61]
                AND e1[64] = e2[64] AND ... e1[61] = e2[61].

    From recurrence: a[r+1] = F(a[r..r-7], K[r], W[r]).
    a1[64] = a2[64] requires: F(a1[63..56], K[63], W1[63]) = F(a2[63..56], K[63], W2[63]).

    If W1[63] = W2[63] (same schedule word): then a1[56..63] = a2[56..63] is needed.
    If W1[63] ≠ W2[63]: the difference ΔW[63] must COMPENSATE the state difference.

    ΔW[63] comes from schedule: ΔW[63] = sig1(ΔW[61]) + ΔW[56] + sig0(ΔW[48]) + ΔW[47].
    All depend on δM = M1 XOR M2.

    The COMPENSATION equation:
      F(a1[56..63], K[63], W1[63]) = F(a2[56..63], K[63], W2[63])

    In L⊕Φ decomposition:
      L(a1,W1) ⊕ Φ(a1,W1) = L(a2,W2) ⊕ Φ(a2,W2)
      L(δa, δW) = Φ(a1,W1) ⊕ Φ(a2,W2)

    This is the PER-ROUND collision equation.
    For round r: linear part = carry difference.
    """
    print("\n" + "=" * 80)
    print("CONVERGENCE REQUIREMENTS: Per-round collision equation")
    print("=" * 80)

    print(f"""
    For collision at round R:
      ∀ r ∈ {{R-7, ..., R}}: a₁[r] = a₂[r] AND e₁[r] = e₂[r]

    Working backward from round R:
      Round R: a₁[R] = F(a₁[R-1..R-8], K[R-1], W₁[R-1])
             = F(a₂[R-1..R-8], K[R-1], W₂[R-1])

    In BTE decomposition:
      a[r+1] = L_round(state[r], W[r]) ⊕ Φ_round(state[r], W[r])

    For convergence at round r+1:
      L_round(s₁, W₁) ⊕ Φ_round(s₁, W₁) = L_round(s₂, W₂) ⊕ Φ_round(s₂, W₂)

    Rearranging:
      L_round(δs, δW) = Φ_round(s₁, W₁) ⊕ Φ_round(s₂, W₂)

    Left: LINEAR in (δs, δW). Known once δM chosen.
    Right: CARRY DIFFERENCE. Depends on actual s₁, s₂.

    For EACH of rounds R-7 through R: one such equation.
    8 equations, each 32 bits = 256-bit system.
    Unknowns: s₁ (= trajectory of M₁) and δM.

    THIS IS THE COLLISION SYSTEM in BTE language.
    It's a system where:
      - Left side (linear) = known from δM choice
      - Right side (carry diff) = structured by T3-T5
      - 8 equations × 32 bits = 256 constraints
      - Schedule connects δW[r] across rounds

    The CARRY DIFFERENCE Φ(s₁,W₁) ⊕ Φ(s₂,W₂) is the CORE.
    It depends on two trajectories and their carry patterns.
    By T5 (cocycle): carries compose by XOR.
    By T3 (nilpotent): each carry dies in n steps.
    By T4 (binomial rank): carry sensitivity = binomial.

    These properties CONSTRAIN the carry difference.
    A random function would have no such constraints.

    THE QUESTION FOR NEXT STEPS:
    Do T3-T5 constraints on carry difference reduce the
    search space below 2^128 (birthday)?

    Current answer: not obviously (Ψ = random on output).
    But: in TRAJECTORY SPACE (not output space), constraints may be visible.
    """)


def experiment_partial_collision():
    """
    Instead of full collision, look for PARTIAL: a1[64] = a2[64] (just one word).
    This is 32-bit match → birthday in 2^16.
    Then: check if BTE structure makes a2[63] match "more likely" than random.
    """
    print("\n" + "=" * 80)
    print("PARTIAL COLLISION: a1[64] = a2[64] (32-bit match)")
    print("=" * 80)

    N = 100000
    hash_to_msg = {}  # key = a[64], value = message seed

    found_collisions = []
    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        states, _ = sha256_round_trace(M)

        a64 = (states[64][0] + H0[0]) & MASK  # hash word 0

        if a64 in hash_to_msg:
            found_collisions.append((hash_to_msg[a64], seed))
        else:
            hash_to_msg[a64] = seed

    print(f"  {len(found_collisions)} partial collisions in {N} messages")
    print(f"  Expected: ~{N**2 // (2 * 2**32)} (birthday for 32 bits)")

    # For each partial collision: check how many OTHER words also match
    if found_collisions:
        for s1, s2 in found_collisions[:3]:
            rng1 = random.Random(s1)
            M1 = [rng1.randint(0, MASK) for _ in range(16)]
            rng2 = random.Random(s2)
            M2 = [rng2.randint(0, MASK) for _ in range(16)]

            states1, _ = sha256_round_trace(M1)
            states2, _ = sha256_round_trace(M2)

            matching_words = sum(1 for i in range(8) if
                                 (states1[64][i]+H0[i])&MASK == (states2[64][i]+H0[i])&MASK)

            print(f"    Seeds {s1},{s2}: {matching_words}/8 hash words match")

            # Check trajectory convergence: when do a-values start matching?
            first_match = None
            for r in range(64, 0, -1):
                if states1[r][0] == states2[r][0]:
                    first_match = r
                else:
                    break

            print(f"    Trajectory matches from round {first_match+1 if first_match else 'never'} to 64")


if __name__ == "__main__":
    experiment_convergence_profile()
    experiment_what_convergence_requires()
    experiment_partial_collision()
