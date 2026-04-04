"""
FIVE BRANCHES growing from BTE Theory.

Each branch = new mathematical direction.
Plant seeds for each. See what grows.
"""

import random, math
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, schedule
)


def branch_1_design():
    """
    BRANCH 1: Design Theory — optimal BTE construction.

    From T7: R_full = n_msg + 2. From T8: rotation = essential.
    From F59: entropy predicts speed. From F60: parity imbalance helps.

    QUESTION: What is the MINIMUM rounds for a secure BTE?
    R_min = R_full = n_msg + 2.
    For n_msg=16: R_min = 18.
    SHA-256 uses 64 = 3.6× R_min.

    Can we design a BTE with R = R_min that is provably secure?
    Parameters: n=32, n_msg=16, R=18, optimal rotations.

    Optimal rotations: max entropy + max parity imbalance.
    = 6 all-odd rotations uniformly distributed in Z/32.
    = {1, 5, 11, 17, 23, 29} (spacing ≈ 5-6 each).
    """
    print("=" * 80)
    print("BRANCH 1: Design Theory — Minimal Secure BTE")
    print("=" * 80)

    # Optimal rotation set: 6 odd values, max entropy in Z/32
    # Divide Z/32 into regions, pick one odd from each
    optimal = [1, 5, 11, 17, 23, 29]  # odd, spread

    # Compare with SHA-256
    sha_rots = [2, 13, 22, 6, 11, 25]

    def rot_entropy(rots, n=32, n_bins=8):
        bin_size = n // n_bins
        bins = [0] * n_bins
        for r in rots:
            bins[r % n // bin_size] += 1
        total = sum(bins)
        if total == 0: return 0
        ent = 0
        for b in bins:
            if b > 0:
                p = b / total
                ent -= p * math.log2(p)
        return ent

    ent_opt = rot_entropy(optimal)
    ent_sha = rot_entropy(sha_rots)

    n_odd_opt = sum(1 for r in optimal if r % 2 == 1)
    n_odd_sha = sum(1 for r in sha_rots if r % 2 == 1)

    print(f"\n  Optimal rotations: {optimal}")
    print(f"    Entropy: {ent_opt:.3f}, Parity: {n_odd_opt}/6 odd")
    print(f"  SHA-256 rotations: {sha_rots}")
    print(f"    Entropy: {ent_sha:.3f}, Parity: {n_odd_sha}/6 odd")
    print(f"\n  Minimal secure BTE: R = n_msg + 2 = 18 rounds")
    print(f"  SHA-256: R = 64 = 3.6× minimum")
    print(f"  A BTE with R=18 and optimal rotations could be EQUALLY secure")
    print(f"  with 3.6× fewer rounds = 3.6× faster hashing.")


def branch_2_algebra():
    """
    BRANCH 2: Abstract Carry Algebra.

    T3-T5 characterize the carry operator. But carry = just ONE example
    of a "sequential Boolean accumulator":
      c[0] = seed
      c[k] = f(x[k-1], y[k-1], c[k-1])

    For f = MAJ: carry. For other f: different accumulators.
    Define: (n, f)-accumulator = sequential recursion with rule f.

    Do T3-T5 generalize to OTHER f?
    T3 (nilpotent): for which f is the accumulator nilpotent?
    T5 (cocycle): for which f is the accumulator a cocycle?
    """
    print("\n" + "=" * 80)
    print("BRANCH 2: Abstract Carry Algebra — generalize T3-T5")
    print("=" * 80)

    n = 8
    mask = (1 << n) - 1

    # Test different accumulator functions f(x, y, c)
    functions = {
        'MAJ': lambda x, y, c: (x & y) | (x & c) | (y & c),  # carry
        'AND3': lambda x, y, c: x & y & c,
        'OR': lambda x, y, c: x | y | c,
        'XOR': lambda x, y, c: x ^ y ^ c,
        'NAND-carry': lambda x, y, c: 1 - ((x & y) | (x & c) | (y & c)),
        'threshold2': lambda x, y, c: 1 if (x + y + c >= 2) else 0,  # = MAJ
    }

    for fname, f in functions.items():
        # Test nilpotency: apply accumulator repeatedly
        nilpotent = True
        max_depth = 0
        for y in range(0, mask + 1, max(1, mask // 30)):
            for x in range(0, mask + 1, max(1, mask // 30)):
                v = x
                for depth in range(n + 2):
                    # Apply accumulator
                    c = 0
                    new_v = 0
                    for k in range(n):
                        new_v |= (c << k)
                        c = f((v >> k) & 1, (y >> k) & 1, c)
                    v = new_v
                    if v == 0:
                        max_depth = max(max_depth, depth + 1)
                        break
                else:
                    nilpotent = False
                    break
            if not nilpotent:
                break

        # Test cocycle: E(a,b,c) = E(a,b) XOR E(a+b, c)?
        # E(a,b) = accumulator_effect(a,b) = result XOR (a XOR b)
        # where result = sum with accumulator f
        cocycle = True
        for _ in range(100):
            a = random.randint(0, mask)
            b = random.randint(0, mask)
            cc = random.randint(0, mask)

            def acc_add(x, y):
                c = 0
                result = 0
                for k in range(n):
                    xk = (x >> k) & 1
                    yk = (y >> k) & 1
                    rk = xk ^ yk ^ c
                    result |= (rk << k)
                    c = f(xk, yk, c)
                return result & mask

            r_abc = acc_add(acc_add(a, b), cc)
            r_direct = acc_add(a, acc_add(b, cc))

            # Cocycle: associativity → E decomposes
            if r_abc != r_direct:
                cocycle = False
                break

        print(f"  {fname:>15}: nilpotent={nilpotent} (depth≤{max_depth}), "
              f"associative={'YES' if cocycle else 'NO'}")


def branch_3_watershed():
    """
    BRANCH 3: Watershed Geometry.

    Watershed = round 57 for SHA-256 collision (from T_watershed).
    But: for STRUCTURED collisions (not random), watershed could differ.

    Define: watershed(M1, M2) = max{r : a1[r] ≠ a2[r]}.
    For random collision pair: watershed = 56 (with P≈1).

    Question: does WANG-LIKE structured pair have earlier watershed?
    Wang: δe[2..17] = 0. This DOESN'T mean a-convergence, but it
    CONSTRAINS the trajectory. Does it shift the watershed?
    """
    print("\n" + "=" * 80)
    print("BRANCH 3: Watershed Geometry")
    print("=" * 80)

    # Can't find collision for full SHA-256.
    # Theoretical analysis:

    print(f"\n  Watershed for different collision types:")
    print(f"")
    print(f"  RANDOM collision (birthday):")
    print(f"    Paths independent until round 57")
    print(f"    Watershed = 56 (last different round)")
    print(f"    Convergence: sudden at round 57")
    print(f"")
    print(f"  WANG-like collision (structured δe=0):")
    print(f"    e-branch converged for rounds 2-17")
    print(f"    a-branch divergent (HW≈16)")
    print(f"    After round 17: δe≠0 again (schedule adds δW)")
    print(f"    Watershed: still ≈56 (a-branch never converges early)")
    print(f"    Wang doesn't shift watershed — only controls e-branch")
    print(f"")
    print(f"  HYPOTHETICAL a-convergence collision:")
    print(f"    If δa[r]=0 from some round r* to 64:")
    print(f"    Watershed = r* - 1")
    print(f"    Requires: δW[r*..63] all satisfy convergence conditions")
    print(f"    From F44: 160 free bits in schedule kernel for δW[53..63]=0")
    print(f"    → Watershed ≥ 53 (can't go earlier with schedule kernel)")
    print(f"")
    print(f"  MINIMUM POSSIBLE WATERSHED:")
    print(f"    Need 8 consecutive δa=0 for state match (a[r*..r*+7])")
    print(f"    Plus створочне gives δe=0 too")
    print(f"    Schedule must support δW[r*..63]=0")
    print(f"    Schedule rank for k words: 32k bits, rank ≤ 32k")
    print(f"    Free: 512 - 32k = 512-32k bits")
    print(f"    For k=64-r*: free = 512-32(64-r*) = 32r*-1536")
    print(f"    Need free > 0: r* > 48")
    print(f"    → MINIMUM WATERSHED ≈ 48 (from schedule capacity)")


def branch_4_complexity():
    """
    BRANCH 4: Complexity Theory of BTE.

    T8: rotation necessary. T7: R_full = n_msg + 2.
    Together: after R_full rounds, function = random.
    Random function collision = 2^{n*4} birthday (for 8n-bit hash).

    Can we PROVE: for BTE with R > R_full, collision ≥ 2^{n*4}?

    Sketch:
    1. After R_full: ALL derivatives Dk = 0.5 (random) → function
       is k-wise independent for all k ≤ 4 (at least).
    2. k-wise independent function: collision probability = same as random.
    3. Random collision: birthday bound = 2^{hash_bits/2}.
    """
    print("\n" + "=" * 80)
    print("BRANCH 4: Complexity Theory — toward provable bounds")
    print("=" * 80)

    print(f"""
    PROOF STRATEGY for birthday lower bound:

    Step 1 (T7): After R_full = n_msg + 2 rounds:
      All derivatives D2, D3, D4 ≈ 0.5 (random).
      This means: output is 4-wise independent (at least).

    Step 2 (information theory): For a function that is
      k-wise independent, the collision probability is:
      P(collision) ≤ (number of pairs) / (output space)
                    = C(N,2) / 2^hash_bits
      This gives birthday bound: need N ≈ 2^(hash_bits/2).

    Step 3: For SHA-256 with R=64 > R_full=20:
      Output is 4-wise independent → birthday bound holds.
      Collision: 2^128. QED (conditional on Step 1-2 formalization).

    GAPS:
    - Step 1: "D_k ≈ 0.5" needs to be formalized as "k-wise independence".
      D_k = 0.5 for random triples → output is STATISTICALLY k-wise independent.
      But this is measured at ONE point (specific M), not globally.
    - Step 2: k-wise independence → birthday is standard.
      But need UNIFORM k-wise independence (over all M, not just one).
    - Step 3: R=64 > R_full=20, but R_full was measured for D_4 only.
      Higher derivatives D_5, D_6, ... not measured.
      Need: ALL derivatives random, not just D2-D4.

    WHAT'S NEEDED:
    1. Prove D_k ≈ 0.5 for ALL k (not just 2,3,4) at R > R_full.
    2. Prove this holds UNIFORMLY over all M (not just specific M).
    3. Connect to standard complexity-theoretic birthday bound.
    """)


def branch_5_coding():
    """
    BRANCH 5: Coding Theory Connection.

    T1: Layer rank = 2R-1. This is a CODE:
    The trajectory at bit-0 forms a [512, 127]-like code over GF(2)
    (512 message bits → 127-dimensional codeword space).

    Minimum distance of this code = measure of collision difficulty at bit 0.
    If minimum distance = d: need to change ≥ d message bits for any
    bit-0 trajectory change.

    From methodology (T_DMIN_97): minimum distance of SHA-256 hash = 97.
    Is bit-0 layer minimum distance also measurable?
    """
    print("\n" + "=" * 80)
    print("BRANCH 5: Coding Theory — BTE layers as codes")
    print("=" * 80)

    # The bit-0 Jacobian = 512×127 matrix over GF(2) (approximately).
    # This defines a [512, 127] code. Minimum distance?

    # For a random [512, 127] code: d_min ≈ 512 × (1/2 - δ) ≈ 200+.
    # For SHA-256 specific code: measure minimum HW of Jacobian rows.

    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]
    states, _ = sha256_round_trace(M)

    # Jacobian of bit-0 system: 512 input bits → 512 bit-0 output values (8 regs × 64 rounds)
    base_bits = []
    for r in range(64):
        for reg in range(8):
            base_bits.append((states[r+1][reg] >> 0) & 1)

    # Row weights (how many output bits each input bit affects)
    row_weights = []
    for word in range(16):
        for bit in range(32):
            M_flip = list(M)
            M_flip[word] ^= (1 << bit)
            sf, _ = sha256_round_trace(M_flip)
            diff = 0
            for r in range(64):
                for reg in range(8):
                    if base_bits[r*8+reg] != ((sf[r+1][reg] >> 0) & 1):
                        diff += 1
            row_weights.append(diff)

    min_w = min(row_weights)
    max_w = max(row_weights)
    avg_w = sum(row_weights) / len(row_weights)

    print(f"\n  Bit-0 layer Jacobian row weights (≈ minimum distance proxy):")
    print(f"    min = {min_w}, max = {max_w}, avg = {avg_w:.1f}")
    print(f"    (random code [512,127] expected d_min ≈ 200+)")
    print(f"")
    print(f"  Interpretation:")
    print(f"    min weight {min_w} = minimum number of bit-0 trajectory")
    print(f"    values that change when flipping ANY single message bit.")
    print(f"    Higher = more diffusion = stronger mixing.")


if __name__ == "__main__":
    branch_1_design()
    branch_2_algebra()
    branch_3_watershed()
    branch_4_complexity()
    branch_5_coding()
