"""
Step 21: Phase Gate — using the OTHER SIDE of the phase transition

Phase transition at bit 5: below = polynomial, above = exponential.
Everyone tries to go FROM poly TO exp. What if we go the OTHER WAY?

Bits 6-31 are "random" (exp complexity). But they CONTAIN information
about bits 0-5 (because they're computed FROM the same state).

Question: can we EXTRACT information about bits 0-5 (easy zone)
from bits 6-31 (hard zone)?

This is INVERSION of the phase transition:
Instead of: solve bits 0-5 (easy) then extend to 6-31 (hard)
Try: READ bits 6-31 (given as output) and DEDUCE bits 0-5

In SHA-256 output: we SEE all 32 bits of each word.
Bits 0-5 of the output are PART OF the output we already have.
But the collision on bits 0-5 depends on intermediate states,
not just the output.

NEW IDEA: use the "random" high bits of INTERMEDIATE states
(which we can infer from the output through backward computation)
to constrain the "structured" low bits.

The high bits KNOW something about the low bits (they were computed
together). The phase transition doesn't destroy the RELATIONSHIP —
it makes the relationship COMPLEX. But if we go BACKWARD from
known high bits, the complexity works FOR us, not against.
"""

import numpy as np
from collections import Counter
from step0_exact_algebra import mini_sha, N, MASK, mobius_transform

N_MSG = 4
N_INPUT = N * N_MSG
N_TOTAL = 1 << N_INPUT


def sha_trace(msg, R):
    """Full trace of mini-SHA."""
    IV = [0x6, 0xB, 0x3, 0xA, 0x5, 0x9, 0x1, 0xF]
    K = [0x4, 0x2, 0xB, 0x7, 0xA, 0x3, 0xE, 0x5,
         0x9, 0x1, 0xD, 0x6, 0x0, 0x8, 0xC, 0xF]
    a, b, c, d, e, f, g, h = IV[:]
    W = list(msg) + [0] * max(0, R - N_MSG)
    trace = [(a, e)]
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
        trace.append((a, e))
    return trace


def test_phase_gate():
    """
    Test: do HIGH bits of output predict LOW bits of collision condition?

    For mini-SHA (n=4):
    - "Low" = bit 0 (carry-free, degree 2)
    - "High" = bits 1-3 (carry-affected, high degree)

    Question: if we know bits 1-3 of the COLLISION DIFFERENCE,
    can we predict bit 0?

    This is asking: does the "hard" part (bits 1-3) contain
    information about the "easy" part (bit 0)?
    """
    import random
    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(N_MSG)]
    R = 8

    print(f"{'='*80}")
    print(f"PHASE GATE TEST — R={R}")
    print(f"{'='*80}")

    base_a, base_e = mini_sha(M, R)

    # For each delta, compute the output difference per bit
    results = []
    for didx in range(1, N_TOTAL):
        dw = []
        tmp = didx
        for w in range(N_MSG):
            dw.append(tmp & MASK)
            tmp >>= N
        M2 = [(M[w] ^ dw[w]) for w in range(N_MSG)]
        a2, e2 = mini_sha(M2, R)
        da = base_a ^ a2
        de = base_e ^ e2
        results.append((didx, da, de))

    # Test 1: If da[bits 1-3] = 0, what's P(da[bit 0] = 0)?
    print(f"\n  TEST 1: P(da[bit 0]=0 | da[bits 1-3]=0)")
    high_zero = [r for r in results if (r[1] >> 1) == 0]  # bits 1-3 of da = 0
    low_zero = [r for r in high_zero if (r[1] & 1) == 0]  # bit 0 also 0

    p_low_given_high = len(low_zero) / max(len(high_zero), 1)
    p_low_uncond = sum(1 for r in results if (r[1] & 1) == 0) / len(results)

    print(f"    P(da[0]=0 | da[1:3]=0) = {p_low_given_high:.4f}")
    print(f"    P(da[0]=0 unconditional) = {p_low_uncond:.4f}")
    print(f"    Lift = {p_low_given_high / max(p_low_uncond, 0.001):.3f}×")

    # Test 2: For ALL possible high-bit patterns, what's the conditional P(low=0)?
    print(f"\n  TEST 2: P(da[0]=0) conditioned on da[1:3] value")
    print(f"  {'da[1:3]':>8} {'Count':>6} {'P(da[0]=0)':>12} {'Lift':>6}")

    for high_val in range(1 << (N-1)):  # 0 to 7 for n=4
        subset = [r for r in results if (r[1] >> 1) == high_val]
        if not subset:
            continue
        p_low = sum(1 for r in subset if (r[1] & 1) == 0) / len(subset)
        lift = p_low / max(p_low_uncond, 0.001)
        marker = " ★" if abs(lift - 1) > 0.1 else ""
        print(f"  {high_val:>8b} {len(subset):>6} {p_low:>12.4f} {lift:>5.2f}×{marker}")

    # Test 3: Same for e-register
    print(f"\n  TEST 3: P(de[0]=0) conditioned on de[1:3] value")
    p_elow_uncond = sum(1 for r in results if (r[2] & 1) == 0) / len(results)
    print(f"  {'de[1:3]':>8} {'Count':>6} {'P(de[0]=0)':>12} {'Lift':>6}")

    for high_val in range(1 << (N-1)):
        subset = [r for r in results if (r[2] >> 1) == high_val]
        if not subset:
            continue
        p_low = sum(1 for r in subset if (r[2] & 1) == 0) / len(subset)
        lift = p_low / max(p_elow_uncond, 0.001)
        marker = " ★" if abs(lift - 1) > 0.1 else ""
        print(f"  {high_val:>8b} {len(subset):>6} {p_low:>12.4f} {lift:>5.2f}×{marker}")

    # Test 4: COMBINED — if BOTH da[1:3]=0 AND de[1:3]=0, what's P(full collision)?
    print(f"\n  TEST 4: Combined gate")
    both_high_zero = [r for r in results if (r[1] >> 1) == 0 and (r[2] >> 1) == 0]
    full_collision = [r for r in results if r[1] == 0 and r[2] == 0]
    full_from_gate = [r for r in both_high_zero if r[1] == 0 and r[2] == 0]

    print(f"    δ with da[1:3]=de[1:3]=0: {len(both_high_zero)}")
    print(f"    Full collisions: {len(full_collision)}")
    print(f"    Full collisions WITHIN gate: {len(full_from_gate)}")

    if both_high_zero:
        precision = len(full_from_gate) / len(both_high_zero)
        recall = len(full_from_gate) / max(len(full_collision), 1)
        print(f"    Precision (collision | gate): {precision:.4f}")
        print(f"    Recall (gate captures): {recall:.4f}")
        print(f"    Gate size: {len(both_high_zero)} (vs brute {N_TOTAL})")
        if len(both_high_zero) > 0 and len(full_from_gate) > 0:
            cost = len(both_high_zero) / len(full_from_gate)
            birthday_cost = N_TOTAL / max(len(full_collision), 1)
            print(f"    Cost per collision (gate): {cost:.1f}")
            print(f"    Cost per collision (brute): {birthday_cost:.1f}")
            print(f"    Speedup: {birthday_cost/cost:.2f}×")

    # Test 5: REVERSE GATE — use output bits 1-3 to find collision on bit 0
    # The idea: filter messages where output high bits match,
    # then check if bit 0 also matches (collision)
    print(f"\n  TEST 5: Reverse gate — filter by output high bits matching")
    print(f"  For each δ, check if H(M)[bits 1-3] = H(M⊕δ)[bits 1-3]")
    print(f"  Then check P(full collision)")

    high_match = [r for r in results if (r[1] & 0xE) == 0 and (r[2] & 0xE) == 0]
    full_in_high = [r for r in high_match if r[1] == 0 and r[2] == 0]

    print(f"    High-bit matches: {len(high_match)}")
    print(f"    Full collisions among them: {len(full_in_high)}")
    if high_match:
        enrichment = (len(full_in_high)/len(high_match)) / (len(full_collision)/len(results))
        print(f"    Enrichment: {enrichment:.2f}× vs random")
        print(f"    Expected enrichment if bits independent: {(1<<(2*(N-1)))/(1<<(2*N))*len(results)/len(high_match):.2f}×")


def main():
    test_phase_gate()


if __name__ == "__main__":
    main()
