"""
BTE ATTACK FRAMEWORK: Using theory for structured collision search.

NOT a complete attack. A FRAMEWORK that identifies where BTE gives advantage.

STRATEGY:
  1. Choose δM in 160-dim schedule kernel (δW[53..63]=0)
  2. This ensures: rounds 53-63 receive NO new message difference
  3. State diff at round 53 must CONVERGE to 0 by round 64
  4. Convergence = 256-bit condition on state[53]
  5. State[53] determined by M (forward computation through 53 rounds)

COST ANALYSIS:
  - Naive: search 160-dim kernel for state[53] convergence = 2^160 trials
  - But: 256-bit convergence condition → P(success) ≈ 2^{-256}
  - Total: 2^160 × 2^{-256} = 2^{-96} expected successes. NOT ENOUGH.

  Need: EITHER more freedom (>256 free bits)
        OR convergence probability > 2^{-256} (structure helps)

WHERE BTE HELPS:
  A. Convergence NOT requiring full 256-bit match (partial convergence)
  B. Carry nilpotency reducing effective convergence dimension
  C. Layer structure making some bits converge "for free"
"""

import random
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, schedule
)


def experiment_convergence_dimension():
    """
    When δW=0 for rounds 53-63: how does state diff EVOLVE?
    If state diff at round 53 = δ: what is state diff at round 64?

    Key: without new message injection (δW=0), the diff evolves through
    pure round function. Does it SHRINK, STAY, or GROW?

    From our theory:
    - Thermalized state: HW(δ) ≈ 128 → stays ≈ 128 (no shrinkage)
    - But: if δ starts SMALL → might stay small for a few rounds
    """
    print("=" * 80)
    print("CONVERGENCE DIMENSION: δstate evolution with δW=0")
    print("=" * 80)

    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]

    # Simulate: inject small δ at round 53, propagate with δW=0 through 53-63
    states, W = sha256_round_trace(M)

    # Start with various initial δ sizes at state[53]
    print(f"\n  Injecting δ at round 53, propagating with δW=0:")
    print(f"  {'Init HW(δ)':>12} | {'HW after 1r':>12} {'after 4r':>12} {'after 8r':>12} {'after 11r':>12}")
    print("-" * 70)

    for init_hw in [1, 2, 4, 8, 16, 32, 64, 128]:
        results = []
        N = 100

        for trial in range(N):
            # Create random δ with given HW in register a
            rng2 = random.Random(trial * 100 + init_hw)

            # Random δstate: flip init_hw random bits in a-register
            delta_a = 0
            bits_to_flip = random.sample(range(32), min(init_hw, 32))
            for b in bits_to_flip:
                delta_a ^= (1 << b)

            # Perturbed state at round 53
            state1 = list(states[53])
            state2 = list(states[53])
            state2[0] ^= delta_a  # Flip bits in a-register only

            # Propagate both with SAME W (δW=0) from round 53 to 64
            hw_per_round = []
            for r in range(53, 64):
                a1,b1,c1,d1,e1,f1,g1,h1 = state1
                T1_1 = (h1+Sig1(e1)+Ch(e1,f1,g1)+K[r]+W[r])&MASK
                T2_1 = (Sig0(a1)+Maj(a1,b1,c1))&MASK
                state1 = [(T1_1+T2_1)&MASK, a1, b1, c1, (d1+T1_1)&MASK, e1, f1, g1]

                a2,b2,c2,d2,e2,f2,g2,h2 = state2
                T1_2 = (h2+Sig1(e2)+Ch(e2,f2,g2)+K[r]+W[r])&MASK
                T2_2 = (Sig0(a2)+Maj(a2,b2,c2))&MASK
                state2 = [(T1_2+T2_2)&MASK, a2, b2, c2, (d2+T1_2)&MASK, e2, f2, g2]

                hw = sum(bin(state1[i] ^ state2[i]).count('1') for i in range(8))
                hw_per_round.append(hw)

            results.append(hw_per_round)

        # Average
        avg = [sum(r[i] for r in results) / N for i in range(11)]
        print(f"  {init_hw:>12} | {avg[0]:>12.1f} {avg[3]:>12.1f} {avg[7]:>12.1f} {avg[10]:>12.1f}")

    print(f"\n  If diff SHRINKS with δW=0: convergence is 'cheap' (probability > 2^{{-256}})")
    print(f"  If diff STAYS or GROWS: convergence requires full 256-bit birthday")


def experiment_partial_convergence():
    """
    Full convergence (256 bits) costs 2^256. Too much.
    What about PARTIAL convergence?

    Hash = (a[64]+IV[0], a[63]+IV[1], ..., e[61]+IV[7]).
    If ONLY a[64] matches (32 bits): partial collision.
    Cost: 2^16 (birthday in 32-bit space).

    With schedule kernel (160 free bits):
    P(a[64] matches) = 2^{-32}. Need 2^32 trials. 2^32 < 2^160. FEASIBLE!

    For a[64] AND a[63] match (64 bits):
    P = 2^{-64}. Need 2^64. Still < 2^160. FEASIBLE!

    For full 256-bit match:
    P = 2^{-256}. Need 2^256. > 2^160. NOT feasible with kernel alone.

    QUESTION: can we use MORE of the kernel to match more hash words?
    """
    print("\n" + "=" * 80)
    print("PARTIAL CONVERGENCE: How many hash words can we match?")
    print("=" * 80)

    print(f"""
    Schedule kernel gives 160 free bits in δM.
    Each hash word = 32-bit condition.
    160 / 32 = 5 hash words matchable (in THEORY, if linear).

    But: matching is NONLINEAR (carry-dependent).
    Effective matching capacity < 160 bits.

    Cost analysis:
      1 hash word (32 bits):   need 2^16 birthday.  160 free → OK!
      2 hash words (64 bits):  need 2^32 birthday.  160 free → OK!
      3 hash words (96 bits):  need 2^48 birthday.  160 free → OK!
      4 hash words (128 bits): need 2^64 birthday.  160 free → OK!
      5 hash words (160 bits): need 2^80 birthday.  160 free → MARGINAL
      8 hash words (256 bits): need 2^128 birthday. 160 free → NOT ENOUGH

    For COLLISION (all 8 words, 256 bits):
      Birthday needs 2^128 pairs.
      Each pair uses 1 δM from 160-dim kernel.
      2^128 < 2^160 → FEASIBLE IN KERNEL!

    Wait — 2^128 < 2^160. The kernel is LARGE ENOUGH for birthday!
    Standard birthday in 160-dim kernel:
      Generate 2^128 pairs from kernel. Each = one δM + random M.
      Check full hash match.
      Expected: ~1 collision.

    This is STANDARD BIRTHDAY but in a 160-dim SUBSPACE.
    NOT faster than regular birthday (still 2^128).

    But: what if the kernel CONCENTRATES collisions?
    I.e., P(collision | δM ∈ kernel) > P(collision | δM random)?
    """)

    # Test: for BTE-8, is collision probability HIGHER in schedule kernel?
    # Can't find collision for BTE-8 64-bit hash (need 2^32).
    # But: can compare HASH DISTANCE for kernel δM vs random δM.

    n = 8; n_msg = 8; R = 16; mask = (1 << n) - 1

    def bte8_hash(msg):
        IV = [random.Random(42+i).randint(0,mask) for i in range(8)]
        K_vals = [random.Random(99999+i).randint(0,mask) for i in range(R)]
        W = list(msg[:n_msg])
        for i in range(n_msg, R):
            W.append((W[i-2] ^ W[max(0,i-n_msg)]) & mask)
        a,b,c,d,e,f,g,h = IV
        for r in range(R):
            s1 = ((e>>2)|(e<<6))&mask ^ ((e>>4)|(e<<4))&mask ^ ((e>>6)|(e<<2))&mask
            ch = (e&f)^(~e&g)&mask
            T1=(h+s1+ch+K_vals[r]+W[r])&mask
            s0 = ((a>>1)|(a<<7))&mask ^ ((a>>3)|(a<<5))&mask ^ ((a>>5)|(a<<3))&mask
            mj=(a&b)^(a&c)^(b&c)
            T2=(s0+mj)&mask
            h,g,f=g,f,e; e=(d+T1)&mask; d,c,b=c,b,a; a=(T1+T2)&mask
        return tuple((r+IV[i])&mask for i,r in enumerate([a,b,c,d,e,f,g,h]))

    # Generate random M and compute hash distances for:
    # A: random δM (512-dim)
    # B: kernel δM (40-dim, from schedule kernel for δW[13..15]=0)
    N = 10000
    rng = random.Random(42)

    dist_random = []
    dist_kernel = []

    for trial in range(N):
        M1 = [random.Random(trial).randint(0, mask) for _ in range(n_msg)]

        # Random δM
        M2_rand = [random.Random(trial + N).randint(0, mask) for _ in range(n_msg)]
        H1 = bte8_hash(M1)
        H2r = bte8_hash(M2_rand)
        dr = sum(bin(H1[i] ^ H2r[i]).count('1') for i in range(8))
        dist_random.append(dr)

        # Kernel δM: only change W[0..12], keep W[13..15] schedule-compatible
        # Simplified: change only W[0]
        M2_kern = list(M1)
        M2_kern[0] ^= random.Random(trial + 2*N).randint(1, mask)
        H2k = bte8_hash(M2_kern)
        dk = sum(bin(H1[i] ^ H2k[i]).count('1') for i in range(8))
        dist_kernel.append(dk)

    avg_rand = sum(dist_random) / N
    avg_kern = sum(dist_kernel) / N

    print(f"\n  BTE-8 hash distances:")
    print(f"    Random δM:  avg HW = {avg_rand:.2f}/{n*8} (random = {n*4})")
    print(f"    Kernel δM:  avg HW = {avg_kern:.2f}/{n*8} (random = {n*4})")
    print(f"    Ratio: {avg_kern/avg_rand:.3f}")

    if avg_kern < avg_rand * 0.9:
        print(f"    *** KERNEL pairs are CLOSER! ***")
    else:
        print(f"    No advantage: kernel pairs same distance as random.")


def experiment_attack_cost():
    """
    ATTACK COST ANALYSIS using BTE theory.
    """
    print("\n" + "=" * 80)
    print("ATTACK COST: What BTE theory gives us")
    print("=" * 80)

    print(f"""
    APPROACH 1: Standard birthday in schedule kernel
      Kernel dim: 160 (for SHA-256, δW[53..63]=0)
      Birthday cost: 2^128 (same as standard — kernel is large enough)
      Advantage: NONE (kernel too large, doesn't concentrate collisions)

    APPROACH 2: Layered birthday (BTE layer structure)
      Layer 0: 8-bit partial collision → birthday 2^4
      Layer 1: +8 bits → birthday 2^8 among layer-0 matches
      ...
      Layer 31: full collision → birthday 2^128
      Total: 2^4 × 2^8 × ... = 2^128 (SAME — layers independent F18)
      Advantage: NONE (layers independent at output)

    APPROACH 3: Wang chain + BTE extension
      Wang: δe[2..17]=0 for free (P=1). Barrier at r=17: P=2^{-32}.
      BTE: explains WHY barrier (a/e asymmetry, F45).
      Extension: use 160-dim kernel for rounds 53-63.
      But: gap rounds 17-52 = thermalized = no structure.
      Advantage: NONE beyond Wang's 2^32 per round extension.

    APPROACH 4: Carry algebra exploitation
      T3-T5 describe carry microscopically.
      But: carry structure invisible on output (F31, F38-40).
      All derivatives random at R>20 (T7).
      Advantage: NONE at output level.

    APPROACH 5: Rotation weakness
      SHA-256 rotations: not fastest possible (F55).
      But: R=64 >> R_full=20. Margin 3.2×.
      Even with suboptimal rotations: fully randomized at R=64.
      Advantage: NONE for full SHA-256.

    APPROACH 6: ??? (NEW, from BTE theory)
      What BTE gives that HASN'T been tried:

      A. TRAJECTORY-LEVEL attack: work in 16384-bit trajectory space,
         not 256-bit hash space. Trajectory has STRUCTURE (T1: layers).
         But: trajectory requires knowing M (circular).

      B. MEDIAN ALGEBRA: carry = median. Median has specific algebraic
         properties (T9-T10). Can these be used for EQUATION SOLVING
         (not just description)?

      C. SIMULTANEOUS layer solving: solve ALL 4 layers at once,
         using cross-layer constraints from rotations.
         Unlike OMEGA (all at once): use LAYER STRUCTURE for speedup.

      D. WATERSHED manipulation: engineer δM that shifts watershed
         earlier (round 48 possible, F62). Earlier watershed =
         more rounds of convergence = easier to achieve.
    """)

    print(f"  MOST PROMISING: Approach 6D (watershed manipulation)")
    print(f"    From F62: minimum watershed ≈ 48 (vs standard 57)")
    print(f"    Shifting watershed from 57 to 48 = 9 more convergence rounds")
    print(f"    9 rounds × 32 bits = 288 additional convergence bits")
    print(f"    But: achieving earlier watershed requires structured δM")
    print(f"    This is OPEN — BTE theory points to it but doesn't solve it")


if __name__ == "__main__":
    experiment_convergence_dimension()
    experiment_partial_convergence()
    experiment_attack_cost()
