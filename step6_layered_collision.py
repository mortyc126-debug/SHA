"""
Step 6: Layered Collision — solving layer by layer

Strategy: find collision by solving bit-by-bit, starting from
the carry-free layer (bit 0) and adding one carry level at a time.

For mini-SHA (n=4, msg=4, 16 input bits, 8 output bits):
  - Full collision: 8 equations, birthday cost ≈ 2^4 = 16
  - Layer 0 (bits 0): 2 equations (a[0], e[0]) — carry-free, degree 2
  - Layer 1 (bits 0-1): 4 equations — one carry level
  - Layer 2 (bits 0-2): 6 equations
  - Layer 3 (bits 0-3): 8 equations (full collision)

Question: can we solve layer 0 CHEAPLY (algebraically),
then extend solutions to deeper layers?

If yes → layered approach beats birthday.
If no → birthday is optimal even with structural knowledge.
"""

import numpy as np
from collections import defaultdict
from step0_exact_algebra import mini_sha, N, MASK

N_MSG = 4
N_INPUT = N * N_MSG  # 16
N_TOTAL = 1 << N_INPUT  # 65536


def find_collisions_by_layer(R, seed_msg=None):
    """
    For a fixed message M, enumerate ALL possible δ (2^16 - 1 nonzero)
    and count how many satisfy collision at each layer level.
    """
    import random
    if seed_msg is None:
        rng = random.Random(42)
        M = [rng.randint(0, MASK) for _ in range(N_MSG)]
    else:
        M = seed_msg

    print(f"\n{'='*80}")
    print(f"LAYERED COLLISION SEARCH — R={R} rounds")
    print(f"Base message M = {[hex(w) for w in M]}")
    print(f"{'='*80}")

    # Compute base hash
    a_base, e_base = mini_sha(M, R)

    # For each nonzero δ, compute hash and check layer by layer
    # δ is applied as XOR to message words
    layer_survivors = {k: 0 for k in range(N + 1)}  # layer k = bits 0..k matched
    collision_deltas = []

    for delta_idx in range(1, N_TOTAL):
        # Decode delta into word-level XOR
        delta_words = []
        tmp = delta_idx
        for w in range(N_MSG):
            delta_words.append(tmp & MASK)
            tmp >>= N

        M2 = [(M[w] ^ delta_words[w]) for w in range(N_MSG)]
        a2, e2 = mini_sha(M2, R)

        # Check layer by layer
        matched_bits = 0
        for bit in range(N):
            a_match = ((a_base >> bit) & 1) == ((a2 >> bit) & 1)
            e_match = ((e_base >> bit) & 1) == ((e2 >> bit) & 1)
            if a_match and e_match:
                matched_bits += 1
            else:
                break

        for k in range(matched_bits + 1):
            layer_survivors[k] += 1

        if matched_bits == N:
            collision_deltas.append((delta_idx, delta_words, M2))

    # Print results
    total = N_TOTAL - 1  # all nonzero deltas

    print(f"\n  LAYER SURVIVAL (of {total} nonzero δ):")
    print(f"  {'Layer':>6} {'Bits matched':>13} {'Survivors':>10} {'Fraction':>10} "
          f"{'Expected (random)':>18} {'Ratio':>7}")

    for k in range(N + 1):
        surv = layer_survivors[k]
        frac = surv / total
        # Random expectation: each bit-pair (a[k], e[k]) filters by 1/4
        # But bit 0 might filter differently
        if k == 0:
            expected = total
        else:
            expected = total / (4 ** k)  # 2 bits per layer (a and e)
        ratio = surv / expected if expected > 0 else float('inf')

        bits_str = f"0..{k-1}" if k > 0 else "(none)"
        print(f"  {k:>6} {bits_str:>13} {surv:>10} {frac:>10.6f} {expected:>18.1f} {ratio:>7.2f}")

    # Collisions found
    print(f"\n  COLLISIONS FOUND: {len(collision_deltas)}")
    if collision_deltas:
        for i, (didx, dwords, m2) in enumerate(collision_deltas[:5]):
            dw_hex = [hex(w) for w in dwords]
            print(f"    δ #{i}: words={dw_hex}")
            print(f"      M  = {[hex(w) for w in M]} → a={hex(a_base)}, e={hex(e_base)}")
            a2, e2 = mini_sha(m2, R)
            print(f"      M' = {[hex(w) for w in m2]} → a={hex(a2)}, e={hex(e2)}")

    return layer_survivors, collision_deltas


def algebraic_layer0_analysis(R, seed_msg=None):
    """
    Analyze the bit-0 layer algebraically.
    Bit 0 = carry-free = degree 2 over GF(2).
    How many solutions does the degree-2 system have?
    """
    import random
    if seed_msg is None:
        rng = random.Random(42)
        M = [rng.randint(0, MASK) for _ in range(N_MSG)]
    else:
        M = seed_msg

    print(f"\n{'='*80}")
    print(f"ALGEBRAIC ANALYSIS OF LAYER 0 — R={R}")
    print(f"{'='*80}")

    a_base, e_base = mini_sha(M, R)

    # Count solutions at each layer for bit 0 only (just a[0] and e[0])
    bit0_solutions = 0
    bit0_only_a = 0
    bit01_solutions = 0
    full_solutions = 0

    # Also track: for bit-0 solutions, how do they distribute across bit 1?
    bit0_sols_list = []

    for delta_idx in range(1, N_TOTAL):
        delta_words = []
        tmp = delta_idx
        for w in range(N_MSG):
            delta_words.append(tmp & MASK)
            tmp >>= N

        M2 = [(M[w] ^ delta_words[w]) for w in range(N_MSG)]
        a2, e2 = mini_sha(M2, R)

        # Bit 0 match
        a0_match = ((a_base ^ a2) & 1) == 0
        e0_match = ((e_base ^ e2) & 1) == 0

        if a0_match:
            bit0_only_a += 1

        if a0_match and e0_match:
            bit0_solutions += 1
            bit0_sols_list.append(delta_idx)

            # Check bit 1
            a1_match = ((a_base ^ a2) >> 1 & 1) == 0
            e1_match = ((e_base ^ e2) >> 1 & 1) == 0
            if a1_match and e1_match:
                bit01_solutions += 1

        # Full match
        if a_base == a2 and e_base == e2:
            full_solutions += 1

    total = N_TOTAL - 1

    print(f"\n  Layer 0 (a[0] match):          {bit0_only_a:>6} / {total} = {bit0_only_a/total*100:.2f}% "
          f"(random: 50%)")
    print(f"  Layer 0 (a[0]+e[0] match):     {bit0_solutions:>6} / {total} = {bit0_solutions/total*100:.2f}% "
          f"(random: 25%)")
    print(f"  Layer 0+1 (bits 0-1 match):    {bit01_solutions:>6} / {total} = {bit01_solutions/total*100:.2f}% "
          f"(random: {100/16:.2f}%)")
    print(f"  Full collision (all bits):      {full_solutions:>6} / {total} = {full_solutions/total*100:.4f}% "
          f"(random: {100/256:.4f}%)")

    # Key metric: filtering power of each layer
    print(f"\n  FILTERING POWER:")
    print(f"  Layer 0 (bit 0): {total} → {bit0_solutions} "
          f"(×{bit0_solutions/total:.4f}, random ×0.25)")
    if bit0_solutions > 0:
        print(f"  Layer 1 (bit 1): {bit0_solutions} → {bit01_solutions} "
              f"(×{bit01_solutions/bit0_solutions:.4f}, random ×0.25)")
    if bit01_solutions > 0:
        print(f"  Full:            {bit01_solutions} → {full_solutions} "
              f"(×{full_solutions/bit01_solutions:.4f}, random ×{1/64:.4f})")

    # Compare filtering power with random expectation
    print(f"\n  LAYERED vs BIRTHDAY:")
    print(f"  Birthday cost for 8-bit output: ~2^4 = 16 evaluations")
    print(f"  Layer 0 solutions: {bit0_solutions} (these are 'free' — degree-2 system)")
    print(f"  If layer 0 has MORE solutions than random → layered is WORSE")
    print(f"  If layer 0 has FEWER solutions → layered might help")

    random_expected_bit0 = total * 0.25  # 2 bits match = 1/4
    ratio = bit0_solutions / random_expected_bit0
    print(f"  Observed/Expected ratio: {ratio:.4f}")
    if ratio > 1.05:
        print(f"  → Layer 0 has MORE solutions than random ({ratio:.2f}×)")
        print(f"    This means bit-0 equations are EASIER than random")
        print(f"    More solutions to filter through → no advantage over birthday")
    elif ratio < 0.95:
        print(f"  → Layer 0 has FEWER solutions than random ({ratio:.2f}×)")
        print(f"    This means bit-0 equations provide EXTRA filtering")
        print(f"    Potential advantage: fewer candidates for deeper layers")
    else:
        print(f"  → Layer 0 ≈ random ({ratio:.2f}×) — no significant advantage")

    return bit0_solutions, bit01_solutions, full_solutions


def multi_seed_analysis(R, n_seeds=20):
    """Run analysis across multiple random messages to get statistics."""
    import random

    print(f"\n{'='*80}")
    print(f"MULTI-SEED ANALYSIS — R={R}, {n_seeds} seeds")
    print(f"{'='*80}")

    bit0_counts = []
    bit01_counts = []
    full_counts = []

    for seed in range(n_seeds):
        rng = random.Random(seed * 137 + 7)
        M = [rng.randint(0, MASK) for _ in range(N_MSG)]

        b0, b01, full = algebraic_layer0_analysis(R, M)
        bit0_counts.append(b0)
        bit01_counts.append(b01)
        full_counts.append(full)

    total = N_TOTAL - 1
    random_bit0 = total * 0.25

    print(f"\n{'='*80}")
    print(f"STATISTICS ACROSS {n_seeds} SEEDS (R={R}):")
    print(f"{'='*80}")
    print(f"  Layer 0 (bit 0 match):")
    print(f"    Mean: {np.mean(bit0_counts):.1f} (random: {random_bit0:.1f})")
    print(f"    Std:  {np.std(bit0_counts):.1f}")
    print(f"    Min:  {min(bit0_counts)}, Max: {max(bit0_counts)}")
    print(f"    Ratio to random: {np.mean(bit0_counts)/random_bit0:.4f}")

    print(f"  Layer 0+1 (bits 0-1 match):")
    random_bit01 = total / 16
    print(f"    Mean: {np.mean(bit01_counts):.1f} (random: {random_bit01:.1f})")
    print(f"    Std:  {np.std(bit01_counts):.1f}")
    print(f"    Ratio: {np.mean(bit01_counts)/random_bit01:.4f}")

    print(f"  Full collision:")
    random_full = total / 256
    print(f"    Mean: {np.mean(full_counts):.1f} (random: {random_full:.1f})")
    print(f"    Found in {sum(1 for f in full_counts if f > 0)}/{n_seeds} seeds")

    # The key question: correlation between layers
    print(f"\n  LAYER CORRELATION:")
    for i in range(n_seeds):
        if bit0_counts[i] > 0:
            filter_rate = bit01_counts[i] / bit0_counts[i]
        else:
            filter_rate = 0
        if i < 5 or full_counts[i] > 0:
            print(f"    Seed {i}: bit0={bit0_counts[i]:>5}, bit01={bit01_counts[i]:>5}, "
                  f"full={full_counts[i]:>3}, filter_01={filter_rate:.4f} (random: 0.0625)")


def main():
    # Start with small rounds to validate
    for R in [4, 6, 8]:
        find_collisions_by_layer(R)

    # Multi-seed statistical analysis at equilibrium
    multi_seed_analysis(R=6, n_seeds=10)


if __name__ == "__main__":
    main()
