"""
Step 22: Compact Fiber Description — describe, don't search

N2: Can we DESCRIBE {M : SHA(M) = H} in < 2^128 bits?

If the fiber (collision set) has COMPRESSIBLE structure,
then its description is shorter than enumeration.

Test on mini-SHA (n=4, R=8):
- Fiber size: ~256 elements in {0,1}^16
- Brute description: 256 × 16 bits = 4096 bits
- Birthday cost to find one: ~16 evaluations

Can we describe the fiber in < 16 bits?
(16 bits = log₂(birthday cost) × bits_per_element)

Approaches:
A. XOR-span: find minimal generating set under XOR
B. Arithmetic span: generate via word addition
C. Recurrence: f(i+1) = g(f(i)) for some simple g
D. Projection: fiber = preimage of a SIMPLER function
E. Grammar: fiber elements follow a production rule
"""

import numpy as np
from collections import Counter
import random as rng_module

N = 4
MASK = (1 << N) - 1
N_MSG = 4
N_INPUT = N * N_MSG
N_TOTAL = 1 << N_INPUT

def mini_sha(msg, R):
    IV = [0x6, 0xB, 0x3, 0xA, 0x5, 0x9, 0x1, 0xF]
    K = [0x4, 0x2, 0xB, 0x7, 0xA, 0x3, 0xE, 0x5,
         0x9, 0x1, 0xD, 0x6, 0x0, 0x8, 0xC, 0xF]
    a, b, c, d, e, f, g, h = IV[:]
    W = list(msg) + [0] * max(0, R - len(msg))
    for r in range(R):
        def rotr(x, s):
            return ((x >> s) | (x << (N - s))) & MASK
        w_r = W[r] if r < len(W) else 0
        k_r = K[r % len(K)]
        sig1 = rotr(e, 1) ^ rotr(e, 3) ^ (e >> 1)
        ch = (e & f) ^ (~e & g) & MASK
        sig0 = rotr(a, 1) ^ rotr(a, 2) ^ rotr(a, 3)
        maj = (a & b) ^ (a & c) ^ (b & c)
        T1 = (h + sig1 + ch + k_r + w_r) & MASK
        T2 = (sig0 + maj) & MASK
        a_new = (T1 + T2) & MASK
        e_new = (d + T1) & MASK
        h, g, f = g, f, e
        e = e_new
        d, c, b = c, b, a
        a = a_new
    return (a, e)


def compute_fiber(M, R):
    base = mini_sha(M, R)
    return sorted([d for d in range(1, N_TOTAL)
                   if mini_sha([(M[w] ^ ((d >> (w*N)) & MASK)) for w in range(N_MSG)], R) == base])


def approach_A_xor_span(fiber):
    """Find minimal XOR generating set."""
    # Greedy: Gaussian elimination
    basis = []
    reduced = list(fiber)
    for bit in range(N_INPUT):
        # Find element with this bit set
        pivot = None
        for elem in reduced:
            if (elem >> bit) & 1:
                pivot = elem
                break
        if pivot is None:
            continue
        basis.append(pivot)
        reduced = [e ^ pivot if (e >> bit) & 1 else e for e in reduced if e != pivot]

    # Check span
    span = {0}
    for b in basis:
        span = span | {s ^ b for s in span}

    fiber_set = set(fiber) | {0}
    captured = len(span & fiber_set) - 1  # minus zero
    return len(basis), len(span), captured


def approach_B_arithmetic_span(fiber):
    """Check if fiber has additive generators."""
    def word_add(a, b):
        result = 0
        for w in range(N_MSG):
            aw = (a >> (w*N)) & MASK
            bw = (b >> (w*N)) & MASK
            result |= ((aw + bw) & MASK) << (w*N)
        return result

    # Try: generate from first few elements using addition
    if len(fiber) < 2:
        return 0, 0, 0

    generated = set()
    seeds = fiber[:3]

    # Generate all sums/differences of seeds
    for s1 in seeds:
        generated.add(s1)
        for s2 in seeds:
            generated.add(word_add(s1, s2))
            # Also try s1 - s2
            neg_s2 = 0
            for w in range(N_MSG):
                neg_s2 |= ((MASK + 1 - ((s2 >> (w*N)) & MASK)) & MASK) << (w*N)
            generated.add(word_add(s1, neg_s2))

    captured = len(generated & set(fiber))
    return len(seeds), len(generated), captured


def approach_C_recurrence(fiber):
    """Check if fiber elements follow a recurrence: f[i+1] = g(f[i])."""
    if len(fiber) < 3:
        return None, 0

    # Sort fiber and look for patterns
    fiber_set = set(fiber)

    best_g = None
    best_count = 0

    # Try simple recurrences: f[i+1] = f[i] XOR c
    for c in range(1, N_TOTAL):
        chain_len = 0
        current = fiber[0]
        while current in fiber_set and chain_len < len(fiber):
            chain_len += 1
            current = current ^ c
        if chain_len > best_count:
            best_count = chain_len
            best_g = ('xor', c)

    # Try: f[i+1] = word_add(f[i], c)
    for c in range(1, min(N_TOTAL, 1000)):
        chain_len = 0
        current = fiber[0]
        while current in fiber_set and chain_len < len(fiber):
            chain_len += 1
            result = 0
            for w in range(N_MSG):
                result |= (((current >> (w*N)) & MASK) + ((c >> (w*N)) & MASK)) & MASK << (w*N)
            current = result & (N_TOTAL - 1)
        if chain_len > best_count:
            best_count = chain_len
            best_g = ('add', c)

    # Try: f[i+1] = ROTR(f[i], s)  (word-level rotation of entire delta)
    for s in range(1, N_INPUT):
        chain_len = 0
        current = fiber[0]
        visited = set()
        while current in fiber_set and current not in visited and chain_len < len(fiber):
            visited.add(current)
            chain_len += 1
            current = ((current >> s) | (current << (N_INPUT - s))) & (N_TOTAL - 1)
        if chain_len > best_count:
            best_count = chain_len
            best_g = ('rotr', s)

    return best_g, best_count


def approach_D_projection(fiber):
    """Check if fiber = preimage of a simpler function."""
    # If there's a function f: {0,1}^16 → {0,1}^k (k < 8) such that
    # fiber = {x : f(x) = c} for some constant c,
    # then the fiber is a "level set" of a simpler function.

    # Test: is there a LINEAR function L: {0,1}^16 → {0,1}^k
    # that's CONSTANT on the fiber?

    fiber_set = set(fiber)

    # For each possible linear projection (single bit):
    # check if it's constant on the fiber
    constant_bits = 0
    for bit in range(N_INPUT):
        values = set((d >> bit) & 1 for d in fiber)
        if len(values) == 1:
            constant_bits += 1

    # For each pair of bits: check XOR
    constant_xor = 0
    for b1 in range(N_INPUT):
        for b2 in range(b1+1, N_INPUT):
            values = set(((d >> b1) & 1) ^ ((d >> b2) & 1) for d in fiber)
            if len(values) == 1:
                constant_xor += 1

    return constant_bits, constant_xor


def approach_E_word_pattern(fiber):
    """Check if fiber has word-level pattern."""
    # Decompose each delta into 4 words and look for word-level regularities
    word_sets = [set() for _ in range(N_MSG)]
    for d in fiber:
        for w in range(N_MSG):
            word_sets[w].add((d >> (w*N)) & MASK)

    # How many distinct values per word?
    word_sizes = [len(ws) for ws in word_sets]

    # Product of word sizes = upper bound on fiber size if independent
    product = 1
    for s in word_sizes:
        product *= s

    return word_sizes, product


def main():
    rng = rng_module.Random(42)
    M = [rng.randint(0, MASK) for _ in range(N_MSG)]
    R = 8

    print(f"{'='*80}")
    print(f"COMPACT FIBER DESCRIPTION — R={R}")
    print(f"{'='*80}")

    fiber = compute_fiber(M, R)
    print(f"  Fiber size: {len(fiber)}")
    print(f"  Brute description: {len(fiber) * N_INPUT} bits")
    print(f"  Birthday cost: ~{2**N} evaluations")

    # A: XOR span
    print(f"\n  APPROACH A: XOR generators")
    n_gen, span_size, captured = approach_A_xor_span(fiber)
    print(f"    Generators: {n_gen}")
    print(f"    Span size: {span_size}")
    print(f"    Fiber captured: {captured}/{len(fiber)} ({captured/len(fiber)*100:.1f}%)")
    print(f"    Description: {n_gen} × {N_INPUT} = {n_gen * N_INPUT} bits")
    if captured == len(fiber) and span_size == len(fiber) + 1:
        print(f"    ★ FIBER = EXACT XOR-SUBSPACE!")
    else:
        print(f"    Span covers {span_size} elements but fiber has {len(fiber)} → NOT subspace")

    # B: Arithmetic span
    print(f"\n  APPROACH B: Arithmetic generators (word-add)")
    n_seeds, gen_size, captured_b = approach_B_arithmetic_span(fiber)
    print(f"    Seeds: {n_seeds}")
    print(f"    Generated: {gen_size}")
    print(f"    Fiber captured: {captured_b}/{len(fiber)} ({captured_b/len(fiber)*100:.1f}%)")

    # C: Recurrence
    print(f"\n  APPROACH C: Recurrence f[i+1] = g(f[i])")
    best_g, chain_len = approach_C_recurrence(fiber)
    print(f"    Best recurrence: {best_g}")
    print(f"    Chain length: {chain_len}/{len(fiber)} ({chain_len/len(fiber)*100:.1f}%)")
    if chain_len > len(fiber) * 0.5:
        print(f"    ★ >50% of fiber follows this recurrence!")
        print(f"    Description: seed ({N_INPUT} bits) + rule → {chain_len} elements")

    # D: Projection (constant bits)
    print(f"\n  APPROACH D: Constant projections")
    const_bits, const_xor = approach_D_projection(fiber)
    print(f"    Constant individual bits: {const_bits}/{N_INPUT}")
    print(f"    Constant XOR pairs: {const_xor}/{N_INPUT*(N_INPUT-1)//2}")
    if const_bits > 0:
        print(f"    ★ {const_bits} bits are FIXED in all fiber elements!")
        print(f"    Effective freedom: {N_INPUT - const_bits} bits (not {N_INPUT})")

    # E: Word pattern
    print(f"\n  APPROACH E: Word-level decomposition")
    word_sizes, product = approach_E_word_pattern(fiber)
    print(f"    Distinct values per word: {word_sizes}")
    print(f"    Product (if independent): {product}")
    print(f"    Actual fiber: {len(fiber)}")
    print(f"    Ratio actual/product: {len(fiber)/product:.4f}")
    if len(fiber) / product < 0.5:
        print(f"    ★ Words are NOT independent — structure exists!")

    # SUMMARY
    print(f"\n{'='*80}")
    print(f"  DESCRIPTION COST SUMMARY")
    print(f"{'='*80}")
    print(f"  Full enumeration:     {len(fiber) * N_INPUT} bits")
    print(f"  XOR generators:       {n_gen * N_INPUT} bits (captures {captured}/{len(fiber)})")
    print(f"  Recurrence:           {N_INPUT + 16} bits (captures {chain_len}/{len(fiber)})")
    print(f"  Fixed bits:           saves {const_bits} bits/element")
    print(f"  Birthday equivalent:  {N} × log₂(evals) ≈ {4 * N} bits")

    # Is any description SHORTER than birthday?
    birthday_bits = 4 * N  # ~16 bits
    best_description = min(n_gen * N_INPUT, N_INPUT + 16)
    print(f"\n  Best compact description: {best_description} bits")
    print(f"  Birthday cost (bits):     {birthday_bits} bits")
    if best_description < birthday_bits:
        print(f"  ★★★ COMPACT DESCRIPTION BEATS BIRTHDAY!")
    else:
        print(f"  Compact description ≥ birthday. No advantage.")


if __name__ == "__main__":
    main()
