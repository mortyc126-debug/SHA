"""
Step 10: Native Algebra — work in SHA's own language

GF(2): carry = 2^n monomials = impossible
Z/2^n: carry = one addition = O(1)

SHA-256 works in Z/2^32. Let's formulate the collision problem
in SHA's NATIVE algebra and see what the equations look like.

Mini-SHA collision at word level:
  - 4 word variables: δW[0..3] ∈ Z/2^4
  - 2 word equations: δa[R] = 0, δe[R] = 0 over Z/2^4
  - Each equation involves: +, ⊕, ROTR, Ch, Maj

Question: what does the WORD-LEVEL collision landscape look like?
And can we navigate it more efficiently than birthday?
"""

import numpy as np
import time
from step0_exact_algebra import mini_sha, N, MASK

N_MSG = 4
N_WORDS = 1 << N  # 16 values per word (4 bits)


def word_collision_landscape(R, M_base):
    """
    Map the COMPLETE collision landscape at word level.

    For each δW combination (4 words × 16 values = 16^4 = 65536),
    compute (δa_R, δe_R). This is the WORD-LEVEL function.

    In Z/2^4: the function is F: (Z/16)^4 → (Z/16)^2.
    Collision = F(δW) = (0, 0).
    """
    a_base, e_base = mini_sha(M_base, R)

    # Map every input to its output
    landscape = {}  # (δW0, δW1, δW2, δW3) → (δa, δe)
    output_count = {}  # (δa, δe) → count

    for dw0 in range(N_WORDS):
        for dw1 in range(N_WORDS):
            for dw2 in range(N_WORDS):
                for dw3 in range(N_WORDS):
                    M2 = [M_base[0]^dw0, M_base[1]^dw1, M_base[2]^dw2, M_base[3]^dw3]
                    a2, e2 = mini_sha(M2, R)
                    da = a2 ^ a_base
                    de = e2 ^ e_base

                    landscape[(dw0,dw1,dw2,dw3)] = (da, de)
                    key = (da, de)
                    output_count[key] = output_count.get(key, 0) + 1

    return landscape, output_count


def analyze_word_structure(R, M_base):
    """Analyze the collision at WORD level."""
    print(f"\n{'='*80}")
    print(f"WORD-LEVEL COLLISION LANDSCAPE — R={R}")
    print(f"M = {[hex(w) for w in M_base]}")
    print(f"{'='*80}")

    t0 = time.time()
    landscape, output_count = word_collision_landscape(R, M_base)
    t_compute = time.time() - t0

    total = N_WORDS ** N_MSG  # 65536
    collisions = output_count.get((0, 0), 0)

    print(f"\n  Domain: (Z/{N_WORDS})^{N_MSG} = {total} points")
    print(f"  Range:  (Z/{N_WORDS})^2 = {N_WORDS**2} possible outputs")
    print(f"  Compute time: {t_compute:.3f}s")
    print(f"  Collisions (δa=δe=0): {collisions}")
    print(f"  Expected (random): {total // (N_WORDS**2)} = {total / N_WORDS**2:.0f}")

    # Output distribution
    print(f"\n  OUTPUT DISTRIBUTION:")
    counts = sorted(output_count.values(), reverse=True)
    print(f"    Max preimages: {counts[0]}")
    print(f"    Min preimages: {counts[-1]}")
    print(f"    Mean: {np.mean(counts):.1f}")
    print(f"    Std: {np.std(counts):.1f}")
    print(f"    Unique outputs: {len(output_count)} / {N_WORDS**2}")

    # Is the output distribution uniform?
    expected = total / (N_WORDS**2)
    chi2 = sum((c - expected)**2 / expected for c in output_count.values())
    print(f"    χ² = {chi2:.1f} (expected for random: ~{N_WORDS**2 - 1})")

    # ================================================================
    # KEY: Structure in the WORD-LEVEL function
    # ================================================================
    print(f"\n  WORD-LEVEL STRUCTURE:")

    # 1. How many words of δ are needed?
    # Fix δW3=0: how many collisions?
    for fixed_word in range(N_MSG):
        count = 0
        for dw in landscape:
            if dw[fixed_word] == 0 and landscape[dw] == (0, 0):
                count += 1
        print(f"    Collisions with δW[{fixed_word}]=0: {count} "
              f"(fraction with δW[{fixed_word}]=0: {count}/{collisions})")

    # 2. Partial collisions: δa=0 only, δe=0 only
    da_zero = sum(1 for (da, de), c in output_count.items() if da == 0 for _ in range(c))
    de_zero = sum(1 for (da, de), c in output_count.items() if de == 0 for _ in range(c))
    print(f"\n    δa=0 (any δe): {da_zero} ({da_zero/total*100:.1f}%)")
    print(f"    δe=0 (any δa): {de_zero} ({de_zero/total*100:.1f}%)")
    print(f"    Both = 0:      {collisions} ({collisions/total*100:.2f}%)")
    print(f"    If independent: {da_zero * de_zero / total:.0f}")
    independent_pred = da_zero * de_zero / total
    print(f"    Ratio actual/independent: {collisions / max(independent_pred, 1):.3f}")

    # 3. The ADDITIVE structure
    # In Z/2^n, δa_R = g(δW) where g involves +, ⊕, Ch, Maj, ROTR
    # Does g have ADDITIVE structure? g(δW1 + δW2) =? g(δW1) + g(δW2)?
    print(f"\n  ADDITIVE LINEARITY TEST (over Z/{N_WORDS}):")
    # Test: g(x+y) = g(x) + g(y) mod 2^n?
    n_test = 1000
    n_linear = 0
    import random
    rng = random.Random(42)
    for _ in range(n_test):
        x = tuple(rng.randint(0, MASK) for _ in range(N_MSG))
        y = tuple(rng.randint(0, MASK) for _ in range(N_MSG))
        xy = tuple((x[i] + y[i]) & MASK for i in range(N_MSG))

        gx = landscape[x]
        gy = landscape[y]
        gxy = landscape[xy]

        gx_plus_gy = ((gx[0] + gy[0]) & MASK, (gx[1] + gy[1]) & MASK)
        if gxy == gx_plus_gy:
            n_linear += 1

    print(f"    g(x+y) = g(x)+g(y) mod {N_WORDS}: {n_linear}/{n_test} "
          f"({n_linear/n_test*100:.1f}%) — {'LINEAR!' if n_linear > n_test*0.95 else 'nonlinear'}")

    # XOR-linearity
    n_xor_linear = 0
    for _ in range(n_test):
        x = tuple(rng.randint(0, MASK) for _ in range(N_MSG))
        y = tuple(rng.randint(0, MASK) for _ in range(N_MSG))
        xy = tuple(x[i] ^ y[i] for i in range(N_MSG))

        gx = landscape[x]
        gy = landscape[y]
        gxy = landscape[xy]

        gx_xor_gy = (gx[0] ^ gy[0], gx[1] ^ gy[1])
        if gxy == gx_xor_gy:
            n_xor_linear += 1

    print(f"    g(x⊕y) = g(x)⊕g(y):        {n_xor_linear}/{n_test} "
          f"({n_xor_linear/n_test*100:.1f}%) — {'XOR-LINEAR!' if n_xor_linear > n_test*0.95 else 'nonlinear'}")

    # ================================================================
    # 4. WORD-LEVEL DIFFERENTIAL: how does δa_R change when we flip ONE word?
    # ================================================================
    print(f"\n  WORD-LEVEL SENSITIVITY:")
    for w in range(N_MSG):
        # For each δ that gives collision, what happens if we change δW[w] by ±1?
        sensitivities = []
        collision_list = [dw for dw in landscape if landscape[dw] == (0, 0)]
        for dw in collision_list[:100]:  # sample
            dw_mod = list(dw)
            dw_mod[w] = (dw_mod[w] + 1) & MASK
            da_new, de_new = landscape[tuple(dw_mod)]
            # How far from collision?
            dist = bin(da_new).count('1') + bin(de_new).count('1')
            sensitivities.append(dist)
        if sensitivities:
            print(f"    δW[{w}] += 1: mean distance = {np.mean(sensitivities):.2f} bits "
                  f"(random: {N:.1f})")

    # ================================================================
    # 5. Can we solve it with WORD-LEVEL Newton's method?
    # ================================================================
    print(f"\n  WORD-LEVEL NEWTON'S METHOD:")
    # Start from random δ, iterate: δ' = δ - J^{-1} * f(δ)
    # J = Jacobian over Z/2^n (word-level)

    # Simple gradient descent in Hamming weight
    n_attempts = 100
    n_success = 0
    total_evals = 0

    for attempt in range(n_attempts):
        # Random start
        dw = [rng.randint(0, MASK) for _ in range(N_MSG)]
        da, de = landscape[tuple(dw)]
        cost = bin(da).count('1') + bin(de).count('1')
        evals = 1

        for step in range(200):
            if cost == 0:
                n_success += 1
                break

            # Try all single-word changes
            best_dw = dw
            best_cost = cost
            for w in range(N_MSG):
                for val in range(N_WORDS):
                    dw_try = list(dw)
                    dw_try[w] = val
                    da_try, de_try = landscape[tuple(dw_try)]
                    c = bin(da_try).count('1') + bin(de_try).count('1')
                    evals += 1
                    if c < best_cost:
                        best_cost = c
                        best_dw = list(dw_try)

            if best_cost >= cost:
                break  # stuck
            dw = best_dw
            cost = best_cost

        total_evals += evals

    avg_evals = total_evals / n_attempts
    print(f"    Success: {n_success}/{n_attempts}")
    print(f"    Avg evaluations: {avg_evals:.0f}")
    print(f"    Birthday cost: ~{N_WORDS} = {N_WORDS}")
    print(f"    Brute force: {total}")
    if n_success > 0:
        print(f"    Greedy speedup vs brute: {total/avg_evals:.1f}×")
        print(f"    Greedy speedup vs birthday: {N_WORDS/avg_evals:.2f}×")

    return collisions, output_count


def main():
    import random
    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(N_MSG)]

    for R in [4, 6, 8]:
        analyze_word_structure(R, M)


if __name__ == "__main__":
    main()
