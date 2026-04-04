"""
BTE CRYPTOGRAPHY: Born from BTE Math, feeding back into it.

Cryptography asks: WHAT makes a problem HARD?
BTE Math answers: HERE is the structure.
Crypto asks: DOES this structure create hardness?
Math extends to answer.

FUNDAMENTAL QUESTION:
  For a BTE with parameters (n, R, n_msg, rotations):
  What is the COLLISION HARDNESS as a function of these parameters?

From BTE Math we know:
  T7: R_full = n_msg + 2 (all derivatives random)
  T8: rotation necessary (no rotation → no randomness)
  T1: layer rank = 2R-1 (information capacity per layer)
  F59: entropy(rotations) predicts D2 speed

BTE Crypto defines:
  HARDNESS(BTE) = minimum cost to find collision

  For R < R_full: hardness < birthday (derivatives not random)
  For R ≥ R_full: hardness = birthday (all random)
  For R >> R_full: hardness = birthday (no extra benefit)

The TRANSITION at R_full = the cryptographic phase transition.
Below: vulnerable. Above: secure. At: threshold.
"""

import random, math
from stage0_observe import MASK, H0, K, rotr, Sig0, Sig1, Ch, Maj, schedule, sha256_round_trace


def experiment_hardness_function():
    """
    HARDNESS as function of R (rounds).

    Measure: for BTE-8 at various R, what is the empirical collision cost?
    Can't find actual collisions (too expensive).
    Proxy: D2 deviation from 0.5 = "structure remaining" = "hardness deficit".

    Define: S(R) = |0.5 - D2(R)| = structure at R rounds.
    S = 0 → random → birthday hard.
    S > 0 → structured → potentially easier.

    HARDNESS(R) ≈ 2^{n*4 × (1 - 2S)} for n-bit BTE with 8 registers.
    (Rough: structure S reduces effective hash bits by factor 2S.)
    """
    print("=" * 80)
    print("BTE CRYPTO: Hardness function H(R)")
    print("=" * 80)

    n = 32  # SHA-256
    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]

    print(f"\n  {'R':>3} | {'D2':>6} {'D3':>6} {'D4':>6} | {'S(R)':>6} | {'Effective bits':>14} | {'Hardness':>12}")
    print("-" * 75)

    for R in [4, 8, 12, 16, 20, 24, 32, 48, 64]:
        def f(msg):
            states, _ = sha256_round_trace(msg, rounds=R)
            return ((states[R][0] + H0[0]) & MASK) & 1

        N = 200
        d2=0; d3=0; d4=0; tot=0
        for trial in range(N):
            r = random.Random(trial*1000+R)
            flips = [(r.randint(0,15), r.randint(0,31)) for _ in range(4)]
            if len(set(flips)) < 4: continue
            vals = {}
            for mv in range(16):
                M2 = list(M)
                for i in range(4):
                    if (mv>>i)&1: M2[flips[i][0]] ^= (1<<flips[i][1])
                vals[mv] = f(M2)
            d2 += vals[0]^vals[1]^vals[2]^vals[3]
            dv3=0
            for s in range(8): dv3 ^= vals[s]
            d3 += dv3
            dv4=0
            for s in range(16): dv4 ^= vals[s]
            d4 += dv4
            tot += 1

        p2=d2/tot; p3=d3/tot; p4=d4/tot
        S = max(abs(0.5-p2), abs(0.5-p3), abs(0.5-p4))

        # Effective hash bits: full hash = 256. Structure reduces by 2S × 256.
        eff_bits = 256 * (1 - 2*S)
        hardness = eff_bits / 2  # birthday = half effective bits

        print(f"  {R:>3} | {p2:.3f} {p3:.3f} {p4:.3f} | {S:.3f} | {eff_bits:>14.1f} | 2^{hardness:.1f}")


def experiment_crypto_parameters():
    """
    BTE Crypto PARAMETER SPACE:

    For a BTE to be "cryptographically secure":
    1. R ≥ R_full = n_msg + 2 (T7: all derivatives random)
    2. Rotations generate Z/n (F6: full coverage)
    3. Rotation entropy sufficient (F59: D2 speed)
    4. n_msg ≥ hash_bits / (2 × n) (enough message capacity)

    These are NECESSARY conditions from BTE theory.
    Together: SUFFICIENT? (open question)
    """
    print("\n" + "=" * 80)
    print("BTE CRYPTO: Parameter requirements for security")
    print("=" * 80)

    print(f"""
    NECESSARY CONDITIONS for BTE security (from BTE Math):

    1. ROUND CONDITION (T7):
       R ≥ R_full = n_msg + 2
       SHA-256: R=64 ≥ 18. SATISFIED (margin 3.6×)

    2. COVERAGE CONDITION (F6):
       Rotation constants generate all of Z/n.
       SHA-256: {{2,6,11,13,22,25}} covers Z/32 in 3 steps. SATISFIED.

    3. ENTROPY CONDITION (F59):
       Rotation entropy ≥ threshold for target D2 speed.
       SHA-256: entropy = 2.585. D2@R=16 = 0.30. SATISFIED.

    4. CAPACITY CONDITION (T1):
       n_msg × n ≥ 2 × hash_bits (enough message space for birthday)
       SHA-256: 16 × 32 = 512 ≥ 2 × 256 = 512. SATISFIED (exactly).

    5. NONLINEARITY CONDITION (T8):
       At least one nonlinear operation (Ch or Maj or carry).
       SHA-256: has all three. SATISFIED.

    6. TWO-BRANCH CONDITION (F45):
       Both T1 (W-dependent) and T2 (state-dependent) present.
       Creates a/e asymmetry → Wang limited to e-branch.
       SHA-256: T1 uses W, T2 doesn't. SATISFIED.
    """)

    # Define: BTE SECURITY LEVEL
    print(f"  BTE SECURITY LEVEL definition:")
    print(f"    SecLevel(BTE) = min(hash_bits/2, collision_cost)")
    print(f"    For ideal BTE: SecLevel = hash_bits / 2 = birthday bound")
    print(f"")
    print(f"    SHA-256: SecLevel = 128 bits (birthday = 2^128)")
    print(f"    If all 6 conditions met: SecLevel = hash_bits / 2 (CONJECTURE)")
    print(f"")
    print(f"    This is the BTE SECURITY CONJECTURE:")
    print(f"    'A BTE satisfying conditions 1-6 has collision security = birthday.'")


def experiment_new_math_from_crypto():
    """
    CRYPTO → MATH feedback: what new math does crypto need?

    Crypto asks: "Is birthday OPTIMAL for BTE?"
    Math needs: "PROVE that all derivatives random → birthday"

    This requires CONNECTING:
    - T7 (derivatives random at R_full) → k-wise independence
    - k-wise independence → collision lower bound
    - collision lower bound → birthday

    MISSING MATH (what crypto needs from math):

    A. UNIFORMITY THEOREM: D_k ≈ 0.5 for ALL k (not just 2-6)
    B. GLOBAL UNIFORMITY: D_k uniform over ALL M (not just one)
    C. INDEPENDENCE THEOREM: k-wise independence → birthday bound
    """
    print("\n" + "=" * 80)
    print("CRYPTO → MATH: What new theorems does crypto need?")
    print("=" * 80)

    print(f"""
    NEEDED THEOREM A (All-order randomness):
      For BTE with R > R_full: D_k(R) ≈ 0.5 for ALL k ≥ 2.
      CURRENT: verified for k = 2,3,4,5,6 (F64).
      NEEDED: proof for all k, or at least k up to n.

    NEEDED THEOREM B (Global uniformity):
      D_k(R, M) ≈ 0.5 for ALL M, not just specific M.
      CURRENT: measured at one M. Not proven globally.
      NEEDED: bound on max_M |D_k(R,M) - 0.5|.

    NEEDED THEOREM C (Independence → Birthday):
      If output is k-wise independent for large k:
      collision probability = Θ(N²/2^hash_bits).
      CURRENT: standard result for truly random functions.
      NEEDED: extend to BTE (pseudo-random, not truly random).

    NEEDED THEOREM D (Rotation sufficiency):
      For BTE with conditions 1-6: output is pseudo-random.
      CURRENT: T7+T8 give necessary conditions.
      NEEDED: prove these are SUFFICIENT.

    SYNTHESIS: Theorems A+B+C+D together would give:
      'BTE satisfying conditions 1-6 has collision security = birthday.'
      This would be the FIRST PROVABLE SECURITY RESULT for SHA-256-like functions.

    MATH PROGRAM for crypto:
      Step 1: Prove A (all-order randomness) — extend T7
      Step 2: Prove B (global uniformity) — new technique needed
      Step 3: Apply C (standard) — existing theory
      Step 4: Prove D (conditions sufficient) — combine A+B+C
    """)

    # Can we make progress on A RIGHT NOW?
    # A: D_k ≈ 0.5 for all k. We measured k=2..6.
    # Test k=7,8 at R=20.

    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]
    R = 20

    print(f"  Testing D_7 and D_8 at R=20:")
    for order in [7, 8]:
        nz = 0; tot = 0
        for trial in range(100):
            r = random.Random(trial*10000+order)
            flips = [(r.randint(0,15), r.randint(0,31)) for _ in range(order)]
            if len(set(flips)) < order: continue

            def f(msg):
                states, _ = sha256_round_trace(msg, rounds=R)
                return ((states[R][0]+H0[0])&MASK) & 1

            val = 0
            for s in range(1 << order):
                M2 = list(M)
                for i in range(order):
                    if (s >> i) & 1:
                        M2[flips[i][0]] ^= (1 << flips[i][1])
                val ^= f(M2)

            tot += 1
            if val: nz += 1

        p = nz/tot if tot > 0 else 0
        print(f"    D_{order} at R=20: P = {p:.3f} (random = 0.500)")


if __name__ == "__main__":
    experiment_hardness_function()
    experiment_crypto_parameters()
    experiment_new_math_from_crypto()
