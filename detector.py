"""
BETTER DETECTOR: The 1.46× signal says correlation EXISTS.
Our tool is too crude to see it fully.

Instead of binary (all-zero / not-all-zero), look at the FULL
8-bit pattern of δH layer-0 and its relationship to layer-1.

δH layer-0 = (δH[0]&1, δH[1]&1, ..., δH[7]&1) = 8-bit vector.
δH layer-1 = (δH[0]>>1&1, ..., δH[7]>>1&1) = 8-bit vector.

There are 256 possible layer-0 patterns and 256 possible layer-1 patterns.
Build the 256×256 transition matrix: P(layer-1 pattern | layer-0 pattern).

If layers independent: every row of the matrix = same distribution.
If layers correlated: rows differ → some layer-0 patterns predict layer-1.
"""

import random
import math
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, schedule
)


def sha256_hash(M):
    states, _ = sha256_round_trace(M)
    return [(states[64][i] + H0[i]) & MASK for i in range(8)]


def experiment_transition_matrix():
    """
    Build the full 256×256 transition matrix between layer-0 and layer-1
    patterns of δH for structured pairs.
    """
    print("=" * 80)
    print("TRANSITION MATRIX: Full layer-0 → layer-1 relationship")
    print("=" * 80)

    N = 200000

    # Count transitions: (layer0_pattern, layer1_pattern) → count
    transitions = {}
    layer0_counts = {}

    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        M2 = list(M); M2[0] ^= 1

        H1 = sha256_hash(M)
        H2 = sha256_hash(M2)
        dH = [H1[i] ^ H2[i] for i in range(8)]

        # Layer-0 pattern: 8 bits
        L0 = 0
        L1 = 0
        for reg in range(8):
            L0 |= ((dH[reg] & 1) << reg)
            L1 |= (((dH[reg] >> 1) & 1) << reg)

        key = (L0, L1)
        transitions[key] = transitions.get(key, 0) + 1
        layer0_counts[L0] = layer0_counts.get(L0, 0) + 1

    # If independent: P(L1 | L0) = P(L1) for all L0.
    # Measure deviation from independence.

    # Global L1 distribution
    L1_global = {}
    for (l0, l1), count in transitions.items():
        L1_global[l1] = L1_global.get(l1, 0) + count

    # For each L0 pattern, compute conditional L1 distribution
    # and measure its KL divergence from global.

    print(f"\n  N = {N}, 256 possible L0 patterns, 256 possible L1 patterns")
    print(f"  Unique L0 patterns observed: {len(layer0_counts)}/256")

    # Compute chi-squared test for independence
    chi2_total = 0
    df = 0
    significant_deviations = 0

    # Focus on L0 patterns with enough samples
    well_sampled = {l0: c for l0, c in layer0_counts.items() if c >= 100}
    print(f"  L0 patterns with ≥100 samples: {len(well_sampled)}")

    # For the top-10 most common L0 patterns:
    top_L0 = sorted(well_sampled.items(), key=lambda x: -x[1])[:10]

    print(f"\n  Top L0 patterns and their most probable L1:")
    for L0_val, L0_count in top_L0:
        # Conditional distribution of L1 given this L0
        cond = {}
        for l1 in range(256):
            cond[l1] = transitions.get((L0_val, l1), 0)

        # Most probable L1
        top_L1 = sorted(cond.items(), key=lambda x: -x[1])[:3]

        # Expected (from global)
        expected_top = sorted(
            [(l1, L1_global.get(l1, 0) * L0_count / N) for l1 in range(256)],
            key=lambda x: -x[1]
        )[:3]

        L0_str = format(L0_val, '08b')
        print(f"    L0={L0_str} (n={L0_count:>5}): "
              f"top L1={top_L1[0][0]:>3}({top_L1[0][1]:>4}) "
              f"{top_L1[1][0]:>3}({top_L1[1][1]:>4}) "
              f"{top_L1[2][0]:>3}({top_L1[2][1]:>4})")

    # KEY TEST: mutual information I(L0; L1)
    # I(L0; L1) = Σ P(l0,l1) log(P(l0,l1) / (P(l0)*P(l1)))
    MI = 0.0
    for (l0, l1), joint_count in transitions.items():
        p_joint = joint_count / N
        p_l0 = layer0_counts.get(l0, 0) / N
        p_l1 = L1_global.get(l1, 0) / N
        if p_joint > 0 and p_l0 > 0 and p_l1 > 0:
            MI += p_joint * math.log2(p_joint / (p_l0 * p_l1))

    print(f"\n  MUTUAL INFORMATION I(L0; L1) = {MI:.6f} bits")
    print(f"  (Independent: 0 bits. Max: 8 bits.)")
    print(f"  This is the EXACT measure of layer correlation.")

    if MI > 0.01:
        print(f"  *** SIGNIFICANT: {MI:.4f} bits of shared info between layers! ***")
    else:
        print(f"  Negligible: layers are effectively independent.")

    # Compare with I(L0; L2), I(L0; L3)
    print(f"\n  Mutual information for other layer pairs:")
    for target_layer in [2, 3, 4, 7]:
        # Recompute transitions for L0 vs L_target
        trans2 = {}
        Ltarget_global = {}
        for seed in range(N):
            rng = random.Random(seed)
            M = [rng.randint(0, MASK) for _ in range(16)]
            M2 = list(M); M2[0] ^= 1
            H1 = sha256_hash(M)
            H2 = sha256_hash(M2)
            dH = [H1[i] ^ H2[i] for i in range(8)]

            L0 = sum(((dH[reg] >> 0) & 1) << reg for reg in range(8))
            Lt = sum(((dH[reg] >> target_layer) & 1) << reg for reg in range(8))

            trans2[(L0, Lt)] = trans2.get((L0, Lt), 0) + 1
            Ltarget_global[Lt] = Ltarget_global.get(Lt, 0) + 1

        MI2 = 0.0
        for (l0, lt), jc in trans2.items():
            pj = jc / N
            pl0 = layer0_counts.get(l0, 0) / N
            plt = Ltarget_global.get(lt, 0) / N
            if pj > 0 and pl0 > 0 and plt > 0:
                MI2 += pj * math.log2(pj / (pl0 * plt))

        print(f"    I(L0; L{target_layer}) = {MI2:.6f} bits")

    # Also: I(L0; L1) for RANDOM pairs (should be ≈ 0)
    print(f"\n  Control: I(L0; L1) for RANDOM pairs (not structured):")
    trans_rand = {}
    L0_rand = {}
    L1_rand_g = {}
    for seed in range(N):
        rng1 = random.Random(seed)
        rng2 = random.Random(seed + N)
        M1 = [rng1.randint(0, MASK) for _ in range(16)]
        M2 = [rng2.randint(0, MASK) for _ in range(16)]
        H1 = sha256_hash(M1)
        H2 = sha256_hash(M2)
        dH = [H1[i] ^ H2[i] for i in range(8)]

        L0 = sum(((dH[reg] >> 0) & 1) << reg for reg in range(8))
        L1 = sum(((dH[reg] >> 1) & 1) << reg for reg in range(8))

        trans_rand[(L0, L1)] = trans_rand.get((L0, L1), 0) + 1
        L0_rand[L0] = L0_rand.get(L0, 0) + 1
        L1_rand_g[L1] = L1_rand_g.get(L1, 0) + 1

    MI_rand = 0.0
    for (l0, l1), jc in trans_rand.items():
        pj = jc / N
        pl0 = L0_rand.get(l0, 0) / N
        pl1 = L1_rand_g.get(l1, 0) / N
        if pj > 0 and pl0 > 0 and pl1 > 0:
            MI_rand += pj * math.log2(pj / (pl0 * pl1))

    print(f"    I(L0; L1) random = {MI_rand:.6f} bits (should be ≈ 0)")
    print(f"    I(L0; L1) structured = {MI:.6f} bits")
    print(f"    Signal = {MI - MI_rand:.6f} bits above noise")


if __name__ == "__main__":
    experiment_transition_matrix()
