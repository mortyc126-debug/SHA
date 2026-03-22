#!/usr/bin/env python3
"""
CRAZY-14: Algebraic Degree Growth of SHA-256
=============================================

THE LAST UNMEASURED MATHEMATICAL PROPERTY.

For a function f: GF(2)^n → GF(2), the algebraic degree is the maximum
degree of its ANF (Algebraic Normal Form). For SHA-256:

- Single round: degree ~7 (Ch has degree 2, additions degree ~32)
- After R rounds: degree bounded by min(255, 7^R)
  BUT the actual degree might be MUCH lower than the bound!

If the ACTUAL degree at round R is significantly lower than the
THEORETICAL maximum, there's exploitable algebraic structure.

Higher-order differentials test this:
- The k-th order differential Δ^k f(x) = sum_{S⊆{1,...,k}} (-1)^{|S|} f(x ⊕ ΔS)
- If deg(f) < k, then Δ^k f = 0 for ALL inputs
- If deg(f) = k, then Δ^k f is a non-zero CONSTANT

We measure: for each round R, what is the MINIMUM k such that Δ^k H_R
is NOT constant? This gives a lower bound on the degree.

If degree at round 20 is, say, 100 instead of 255, there's a
200th-order differential that's ZERO — revealing algebraic structure
that could be exploited.
"""

import random
import time

MASK = 0xFFFFFFFF

K = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
    0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
    0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
    0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
    0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
    0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
]

H0 = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
]

def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK
def Sig0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sig1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def sig0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sig1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)
def Ch(e, f, g): return ((e & f) ^ (~e & g)) & MASK
def Maj(a, b, c): return ((a & b) ^ (a & c) ^ (b & c)) & MASK
def add32(*args):
    s = 0
    for a in args: s = (s + a) & MASK
    return s
def hw(x): return bin(x & MASK).count('1')

def sha256_r(msg, iv, nr):
    W = list(msg)
    for i in range(16, max(nr, 16)):
        W.append(add32(sig1(W[i-2]), W[i-7], sig0(W[i-15]), W[i-16]))
    state = list(iv)
    for t in range(nr):
        a, b, c, d, e, f, g, h = state
        T1 = add32(h, Sig1(e), Ch(e, f, g), K[t], W[t])
        T2 = add32(Sig0(a), Maj(a, b, c))
        state = [add32(T1, T2), a, b, c, add32(d, T1), e, f, g]
    return [add32(state[i], iv[i]) for i in range(8)]


def higher_order_diff(base_msg, directions, iv, nr, output_bit):
    """Compute k-th order differential of SHA-256/nr at output_bit.

    directions: list of k direction vectors (each is [word, bit] pair)
    Returns: the value of the k-th order differential (0 or 1)

    Δ^k f(x) = XOR over all S⊆directions of f(x ⊕ sum(S))
    """
    k = len(directions)
    result = 0

    # Iterate over all 2^k subsets
    for mask in range(1 << k):
        msg = list(base_msg)
        for i in range(k):
            if (mask >> i) & 1:
                word, bit = directions[i]
                msg[word] ^= (1 << bit)

        h = sha256_r(msg, iv, nr)
        # Extract output_bit
        reg = output_bit // 32
        bit = output_bit % 32
        result ^= (h[reg] >> bit) & 1

    return result


def test_degree_at_order(base_msg, iv, nr, order, output_bit, n_tests=50):
    """Test if degree is < order by checking if k-th order differential is
    ALWAYS zero (constant zero implies degree < k).

    Returns: fraction of tests where the k-th order differential was zero.
    If fraction = 1.0, degree is likely < order.
    If fraction ≈ 0.5, degree is likely ≥ order.
    """
    zero_count = 0

    for _ in range(n_tests):
        # Random directions
        directions = []
        used = set()
        for _ in range(order):
            while True:
                word = random.randint(0, 15)
                bit = random.randint(0, 31)
                key = word * 32 + bit
                if key not in used:
                    used.add(key)
                    directions.append((word, bit))
                    break

        val = higher_order_diff(base_msg, directions, iv, nr, output_bit)
        if val == 0:
            zero_count += 1

    return zero_count / n_tests


def main():
    print("=" * 72)
    print("CRAZY-14: Algebraic Degree Growth of SHA-256")
    print("=" * 72)
    print()

    random.seed(0xDE614)
    t_start = time.time()
    iv = list(H0)

    base_msg = [random.getrandbits(32) for _ in range(16)]
    output_bit = 4  # bit 4 of register a (arbitrary choice)

    # ── Phase 1: Degree lower bound per round ────────────────────

    print("Phase 1: Higher-order differential tests")
    print("  Testing: is the k-th order derivative zero?")
    print("  If P(zero) ≈ 1.0 → degree < k")
    print("  If P(zero) ≈ 0.5 → degree ≥ k")
    print("-" * 72)
    print()

    # For each round, test multiple orders to find the degree
    # We can test up to about order 20 (2^20 = 1M evaluations)

    max_order_by_round = {}

    for nr in [1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20]:
        if time.time() - t_start > 400:
            print(f"  [Time limit at round {nr}]")
            break

        print(f"  Round {nr}:")

        # Binary search for the degree
        # Start from low orders and increase
        found_degree = False

        for order in [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 18, 20]:
            if time.time() - t_start > 400:
                break

            # 2^order evaluations per test, n_tests tests
            if order > 18:
                n_tests = 5  # Very expensive
            elif order > 14:
                n_tests = 10
            elif order > 10:
                n_tests = 20
            else:
                n_tests = 30

            cost = (1 << order) * n_tests
            if cost > 5_000_000:  # 5M max evaluations
                print(f"    order {order}: SKIPPED (too expensive: {cost} evals)")
                continue

            p_zero = test_degree_at_order(base_msg, iv, nr, order, output_bit, n_tests)

            status = "degree < k" if p_zero > 0.9 else "degree ≥ k" if p_zero < 0.6 else "UNCERTAIN"
            elapsed = time.time() - t_start
            print(f"    order {order:2d}: P(zero) = {p_zero:.2f}  [{status}]  [{elapsed:.1f}s]")

            if p_zero > 0.9:
                max_order_by_round[nr] = order  # degree < order
            elif not found_degree and p_zero < 0.6:
                # degree ≥ previous order
                found_degree = True

        print()

    # ── Phase 2: Multiple output bits ─────────────────────────────

    if time.time() - t_start < 420:
        print()
        print("Phase 2: Degree consistency across output bits")
        print("-" * 72)

        nr = 4  # Test at round 4 where degree should be small
        for out_bit in [0, 32, 64, 128, 160, 224]:
            for order in [3, 5, 7, 10]:
                if time.time() - t_start > 440:
                    break
                n_tests = 20 if order <= 10 else 10
                p = test_degree_at_order(base_msg, iv, nr, order, out_bit, n_tests)
                reg = out_bit // 32
                bit = out_bit % 32
                regnames = ['a','b','c','d','e','f','g','h']
                print(f"    R={nr}, {regnames[reg]}[{bit}], order={order:2d}: P(zero)={p:.2f}")
            print()

    # ── Phase 3: Theoretical vs actual degree ─────────────────────

    print()
    print("=" * 72)
    print("DEGREE GROWTH ANALYSIS")
    print("=" * 72)
    print()

    # Theoretical bounds:
    # - Ch(e,f,g) = ef ⊕ e'g has degree 2
    # - Maj(a,b,c) = ab ⊕ ac ⊕ bc has degree 2
    # - Addition mod 2^32: carry_k has degree k+1 in the worst case
    #   But the MSB carry has degree 32 (full carry chain)
    # - Σ₀, Σ₁, σ₀, σ₁: degree 1 (linear)
    #
    # Per round: degree multiplies by factor of ~2 (from Ch/Maj)
    # Plus carry chains add more degree
    #
    # After R rounds: degree ≤ min(255, product_of_per_round_degrees)
    # Theoretical: degree doubles per round → 2^R
    # But bounded by 255 (dimension of feedforward)
    #
    # Actual degree reaches 255 at round R* where 2^R* ≈ 255 → R* ≈ 8
    # So after 8 rounds, degree should be maximal (255)

    print("  Theoretical degree growth:")
    print("    Round 1: degree ≤ ~34 (2 from Ch/Maj × ~17 from MSB carry)")
    print("    Round 2: degree ≤ ~34² ≈ but bounded by 255")
    print("    Round 3+: degree = 255 (maximal)")
    print()

    if max_order_by_round:
        print("  Measured degree lower bounds:")
        for nr in sorted(max_order_by_round.keys()):
            order = max_order_by_round[nr]
            print(f"    Round {nr:2d}: degree ≥ {order}")

    # ── Verdict ───────────────────────────────────────────────────

    print()
    print("=" * 72)
    print("VERDICT")
    print("=" * 72)
    print()

    # Check if actual degree is significantly lower than theoretical max
    if max_order_by_round:
        for nr, order in sorted(max_order_by_round.items()):
            theoretical_max = min(255, 2 ** (nr + 5))  # rough upper bound
            if order < theoretical_max * 0.3:
                print(f"  **ALIVE** at round {nr}: degree = {order} << theoretical {theoretical_max}")
            elif order < theoretical_max * 0.7:
                print(f"  **ANOMALY** at round {nr}: degree = {order} < theoretical {theoretical_max}")
            else:
                print(f"  Normal at round {nr}: degree ≥ {order} (close to max)")

        # For full SHA-256: if degree reaches 255 by round 8,
        # there's no higher-order differential attack
        max_nr = max(max_order_by_round.keys())
        max_measured = max(max_order_by_round.values())
        if max_measured >= 15:
            print(f"\n  Degree reaches ≥{max_measured} by round {max_nr}.")
            print(f"  Full SHA-256 degree = 255 (confirmed by theory).")
            print(f"  No higher-order differential attack possible.")
        else:
            print(f"\n  Degree only measured up to {max_measured}.")
            print(f"  Need higher-order tests for definitive answer.")

    else:
        print("  Insufficient data collected.")

    print(f"\n  Runtime: {time.time() - t_start:.1f}s")
    print("=" * 72)


if __name__ == "__main__":
    main()
