"""
Step 16: Bit-Free Approaches — can we avoid bits entirely?

Three ideas:
1. Mersenne ring Z/(2^n-1): rotation = multiplication, addition native
2. Masked Word Algebra: Ch = single MUX operation on word triples
3. Word-level orbit analysis: collision = orbit matching in function space

Test: formulate mini-SHA collision with ZERO bit-level operations.
Measure: does the "word count" (not bit count) change with rounds?
"""

import numpy as np
from step0_exact_algebra import mini_sha, N, MASK

N_MSG = 4
N_INPUT = N * N_MSG
N_TOTAL = 1 << N_INPUT
MERSENNE = (1 << N) - 1  # 2^4 - 1 = 15


def to_mersenne(x):
    """Convert x mod 2^n to x mod (2^n - 1)."""
    return x % MERSENNE if MERSENNE > 0 else 0


def mersenne_add(a, b):
    """Addition in Z/(2^n - 1)."""
    s = a + b
    # Carry out = floor((a+b) / 2^n)
    carry_out = s >> N
    result = (s & MASK) + carry_out
    if result >= MERSENNE:
        result -= MERSENNE
    return result


def mersenne_rotr(x, s):
    """Rotation in Mersenne ring = multiplication by 2^{-s}."""
    # ROTR(x, s) = x * 2^{n-s} mod (2^n - 1) for x != 0
    # Actually: ROTR(x, s) in bits = x * 2^{-s} mod (2^n - 1)
    # Since 2^n ≡ 1 mod (2^n-1), 2^{-s} ≡ 2^{n-s} mod (2^n-1)
    if x == 0:
        return 0
    return (x * pow(2, N - s, MERSENNE)) % MERSENNE


def test_mersenne_rotation():
    """Verify: ROTR in bits = multiplication in Mersenne."""
    print(f"\n  MERSENNE ROTATION TEST (n={N}, Mersenne = {MERSENNE}):")

    def rotr_bits(x, s):
        return ((x >> s) | (x << (N - s))) & MASK

    errors = 0
    total = 0
    for x in range(1 << N):  # skip 0? Actually 0 maps to 0 in both
        for s in range(1, N):
            bit_result = rotr_bits(x, s)
            if x == MASK:  # 2^n - 1 maps to 0 in Mersenne
                mersenne_x = 0
            else:
                mersenne_x = x % MERSENNE
            mer_result = mersenne_rotr(mersenne_x, s)

            # Convert back: Mersenne 0 could be 0 or 2^n-1 in bits
            if bit_result % MERSENNE != mer_result:
                errors += 1
            total += 1

    print(f"  Rotation match: {total - errors}/{total} ({(total-errors)/total*100:.1f}%)")
    if errors > 0:
        print(f"  Mismatches: {errors} (likely from x=2^n-1 ↔ 0 ambiguity)")
    return errors == 0


def word_level_collision_analysis(R, M_base):
    """
    Analyze the collision at PURE WORD LEVEL.

    The collision condition: 2 words match (δa=0, δe=0).
    This is 2 equations in word space.

    Question: does the word-level view give us FEWER "variables"
    than the bit-level view?

    At word level: 4 message words → 4 "word variables"
    At bit level: 16 message bits → 16 "bit variables"

    Birthday at word level: ~2^(output_words * word_bits / 2) = 2^(2*4/2) = 2^4
    Same as bit-level birthday.

    But the INTERNAL structure might be different.
    """
    print(f"\n{'='*80}")
    print(f"WORD-LEVEL ORBIT ANALYSIS — R={R}")
    print(f"{'='*80}")

    a_base, e_base = mini_sha(M_base, R)

    # The collision function at WORD level:
    # F: (Z/2^n)^4 → (Z/2^n)^2
    # F(dw0, dw1, dw2, dw3) = (δa_R, δe_R)

    # Word-level "dimension": 4 input words, 2 output words
    # Collision: find (dw0..dw3) with F = (0, 0)

    # Key question: is F a POLYNOMIAL over Z/2^n?
    # If yes: what degree?
    # Polynomial over Z/2^n means: F(x) = Σ a_i x^i mod 2^n

    # Test: is F(dw0, 0, 0, 0) a polynomial in dw0?
    print(f"\n  F(δW[0], 0, 0, 0) — is it polynomial in δW[0] over Z/{MASK+1}?")

    # Compute F for all δW[0] values, keeping δW[1..3] = 0
    f_vals_a = []
    f_vals_e = []
    for dw0 in range(1 << N):
        M2 = [(M_base[0] ^ dw0)] + list(M_base[1:])
        a2, e2 = mini_sha(M2, R)
        f_vals_a.append(a2 ^ a_base)
        f_vals_e.append(e2 ^ e_base)

    # Check: is the sequence f_vals_a a polynomial function mod 2^n?
    # A polynomial of degree d over Z/2^n: Σ_{i=0}^d a_i x^i mod 2^n
    # For d < 2^n: the function values uniquely determine the polynomial.
    # We can use finite differences to find the degree.

    # Finite differences: Δf(x) = f(x+1) - f(x) mod 2^n
    # Δ^k f(x) = 0 for all x iff deg(f) < k

    print(f"\n  Finite differences of δa as function of δW[0]:")
    current = [a ^ a_base for a in [mini_sha([(M_base[0] ^ dw0)] + list(M_base[1:]), R)[0]
               for dw0 in range(1 << N)]]

    for order in range(N + 1):
        max_val = max(current)
        min_val = min(current)
        nonzero = sum(1 for v in current if v != 0)
        all_zero = all(v == 0 for v in current)
        print(f"    Δ^{order}: {'ALL ZERO → degree < ' + str(order) if all_zero else f'nonzero={nonzero}/{len(current)}, max={max_val}'}")
        if all_zero:
            break
        # Next order difference
        new_current = []
        for i in range(len(current) - 1):
            new_current.append((current[i+1] - current[i]) & MASK)
        current = new_current
        if not current:
            break

    # ================================================================
    # ORBIT STRUCTURE
    # ================================================================
    print(f"\n  ORBIT STRUCTURE:")
    print(f"  How many distinct outputs does F produce?")

    outputs = {}
    for didx in range(N_TOTAL):
        dw = []
        tmp = didx
        for w in range(N_MSG):
            dw.append(tmp & MASK)
            tmp >>= N
        M2 = [(M_base[w] ^ dw[w]) for w in range(N_MSG)]
        a2, e2 = mini_sha(M2, R)
        out = (a2 ^ a_base, e2 ^ e_base)
        outputs[out] = outputs.get(out, 0) + 1

    n_outputs = len(outputs)
    collisions = outputs.get((0, 0), 0)
    max_fiber = max(outputs.values())
    min_fiber = min(outputs.values())

    print(f"  Distinct outputs: {n_outputs} / {(1<<N)**2} = {n_outputs/(1<<N)**2*100:.1f}%")
    print(f"  Collisions (0,0): {collisions}")
    print(f"  Max fiber size: {max_fiber}")
    print(f"  Min fiber size: {min_fiber}")
    print(f"  Fiber std: {np.std(list(outputs.values())):.1f}")

    # ================================================================
    # THE KEY BIT-FREE QUESTION
    # ================================================================
    print(f"\n  BIT-FREE FORMULATION:")
    print(f"  Variables: {N_MSG} words ∈ Z/{1<<N}")
    print(f"  Equations: 2 words = 0")
    print(f"  Over Z/{1<<N}: this is 2 equations in 4 unknowns")
    print(f"  Generic solution count: {(1<<N)**(N_MSG-2)} = {(1<<N)**(N_MSG-2)}")
    print(f"  Actual solution count: {collisions}")
    print(f"  Ratio: {collisions / (1<<N)**(N_MSG-2):.3f}")

    # Is the function "polynomial" at word level?
    # If degree d in Z/2^n: can be solved by d-th order methods
    print(f"\n  WORD-LEVEL DEGREE (via finite differences on FULL function):")

    # Test along each axis
    for w_idx in range(N_MSG):
        diffs = []
        for dw_val in range(1 << N):
            dw = [0] * N_MSG
            dw[w_idx] = dw_val
            M2 = [(M_base[i] ^ dw[i]) for i in range(N_MSG)]
            a2, e2 = mini_sha(M2, R)
            diffs.append((a2 ^ a_base))

        # Compute finite difference order
        current = list(diffs)
        deg = N  # default: full
        for order in range(1, N + 1):
            new = [(current[i+1] - current[i]) & MASK for i in range(len(current)-1)]
            if all(v == 0 for v in new):
                deg = order
                break
            current = new

        print(f"    δW[{w_idx}] axis: word-degree = {deg}")


def main():
    import random
    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(N_MSG)]

    test_mersenne_rotation()

    for R in [1, 3, 4, 6, 8]:
        word_level_collision_analysis(R, M)


if __name__ == "__main__":
    main()
