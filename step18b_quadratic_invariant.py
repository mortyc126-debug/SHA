#!/usr/bin/env python3
"""Step 18b: Deep analysis of the quadratic semi-invariant.

DISCOVERY: The function Q(state) = (a&e) XOR (b&f) XOR (c&g) XOR (d&h)
is preserved by the SHA-256 round function with P = 0.625, not 0.5!

This is a 25% bias — FAR from random. What does this mean?

Note: this is a per-BIT function (produces 32 bits).
Let's analyze it bit by bit and understand WHY it's preserved.

Algebraic meaning: Q pairs up the (a,b,c,d) and (e,f,g,h) registers
via AND. This is the INNER PRODUCT of the top and bottom halves
of the shift register!

The shift register maps:
  (a,b,c,d) → (a_new, a, b, c)
  (e,f,g,h) → (e_new, e, f, g)

So Q_new = (a_new & e_new) XOR (a & e) XOR (b & f) XOR (c & g)
         = (a_new & e_new) XOR Q_old XOR (d & h)
            - wait, Q_old = (a&e) XOR (b&f) XOR (c&g) XOR (d&h)
         = (a_new & e_new) - (d & h) + Q_old  (XOR arithmetic)

Hmm, not quite. Let me be more careful.

Q_old = (a&e) XOR (b&f) XOR (c&g) XOR (d&h)
Q_new = (a_new & e_new) XOR (b_new & f_new) XOR (c_new & g_new) XOR (d_new & h_new)
      = (a_new & e_new) XOR (a & e) XOR (b & f) XOR (c & g)

So: Q_new XOR Q_old = (a_new & e_new) XOR (a & e) XOR (d & h)  [cancel common terms]
Wait: Q_new = (a_new&e_new) XOR (a&e) XOR (b&f) XOR (c&g)
      Q_old = (a&e) XOR (b&f) XOR (c&g) XOR (d&h)
Q_new XOR Q_old = (a_new&e_new) XOR (d&h)

So Q is preserved iff (a_new & e_new) = (d & h) at each bit!

a_new = T1 + T2
e_new = d + T1

So: (T1+T2) & (d+T1) = d & h (bitwise, modular arithmetic)

This is the condition we need to study!
"""

import random
import math
from collections import defaultdict
MASK32 = 0xFFFFFFFF

def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
def hw(x): return bin(x).count('1')

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
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
]

def sha256_step(state, W_t, K_t):
    a, b, c, d, e, f, g, h = state
    T1 = (h + Sigma1(e) + Ch(e, f, g) + K_t + W_t) & MASK32
    T2 = (Sigma0(a) + Maj(a, b, c)) & MASK32
    return ((T1 + T2) & MASK32, a, b, c, (d + T1) & MASK32, e, f, g)

def Q(state):
    """The quadratic semi-invariant."""
    a, b, c, d, e, f, g, h = state
    return (a & e) ^ (b & f) ^ (c & g) ^ (d & h)

# ============================================================
# Part 1: Per-bit analysis of the invariant
# ============================================================
def per_bit_analysis():
    """Check which bits of Q are most strongly preserved."""
    print("=" * 70)
    print("Part 1: Per-bit preservation of Q = ae XOR bf XOR cg XOR dh")
    print("=" * 70)

    random.seed(42)
    N = 100000

    bit_preserved = [0] * 32

    for _ in range(N):
        s = tuple(random.randint(0, MASK32) for _ in range(8))
        W = random.randint(0, MASK32)
        s_new = sha256_step(s, W, K[0])

        q_old = Q(s)
        q_new = Q(s_new)
        q_diff = q_old ^ q_new  # Bits that changed

        for bit in range(32):
            if not ((q_diff >> bit) & 1):
                bit_preserved[bit] += 1

    print(f"\n  Per-bit P(Q[bit] preserved):")
    for bit in range(32):
        p = bit_preserved[bit] / N
        bias = p - 0.5
        bar = "█" * int(abs(bias) * 200)
        print(f"  bit {bit:2d}: P = {p:.4f} (bias = {bias:+.4f}) {bar}")

    avg_p = sum(bit_preserved) / (32 * N)
    print(f"\n  Average P(bit preserved) = {avg_p:.4f}")
    print(f"  Total preserved bits per round = {avg_p * 32:.1f} / 32")

# ============================================================
# Part 2: WHY is it preserved? Algebraic analysis
# ============================================================
def algebraic_analysis():
    """Q_new XOR Q_old = (a_new & e_new) XOR (d & h).

    a_new = T1 + T2 = (h + Sigma1(e) + Ch(e,f,g) + K + W) + (Sigma0(a) + Maj(a,b,c))
    e_new = d + T1 = d + h + Sigma1(e) + Ch(e,f,g) + K + W

    Let T = T1 = h + Sigma1(e) + Ch(e,f,g) + K + W
    a_new = T + T2
    e_new = d + T

    (a_new & e_new) = (T + T2) & (d + T)

    For Q to be preserved: (T + T2) & (d + T) = d & h

    This is a condition relating T, T2, d, and h.
    T depends on h, e, f, g, K, W.
    T2 depends on a, b, c.

    Key insight: T includes h. So T = h + rest.
    Let R = Sigma1(e) + Ch(e,f,g) + K + W (rest of T1)
    Then T = h + R, a_new = h + R + T2, e_new = d + h + R

    (h + R + T2) & (d + h + R) = d & h

    Let S = h + R. Then:
    (S + T2) & (d + S) = d & h

    The question is: for random S and T2, when does (S+T2) & (d+S) = d & h?
    """
    print("\n" + "=" * 70)
    print("Part 2: Algebraic analysis of preservation condition")
    print("=" * 70)

    random.seed(42)
    N = 500000

    # Test: for random S, T2, d, h: P((S+T2) & (d+S) = d & h) per bit
    preserved_per_bit = [0] * 32
    total_preserved_all = 0

    for _ in range(N):
        S = random.randint(0, MASK32)
        T2 = random.randint(0, MASK32)
        d = random.randint(0, MASK32)
        h = random.randint(0, MASK32)

        lhs = (S + T2) & MASK32
        rhs = (d + S) & MASK32
        prod = lhs & rhs
        target = d & h

        diff = prod ^ target
        if diff == 0:
            total_preserved_all += 1

        for bit in range(32):
            if not ((diff >> bit) & 1):
                preserved_per_bit[bit] += 1

    print(f"\n  For RANDOM (S, T2, d, h):")
    print(f"  P(all 32 bits preserved) = {total_preserved_all/N:.6f}")
    print(f"  Per-bit: avg P = {sum(preserved_per_bit)/(32*N):.4f}")

    # Now compare with ACTUAL SHA-256 (where S, T2, d, h are correlated)
    preserved_sha = [0] * 32
    total_sha = 0

    for _ in range(N):
        s = tuple(random.randint(0, MASK32) for _ in range(8))
        W = random.randint(0, MASK32)
        a, b, c, d, e, f, g, h = s

        R = (Sigma1(e) + Ch(e, f, g) + K[0] + W) & MASK32
        S = (h + R) & MASK32
        T2 = (Sigma0(a) + Maj(a, b, c)) & MASK32

        lhs = (S + T2) & MASK32
        rhs = (d + S) & MASK32
        prod = lhs & rhs
        target = d & h

        diff = prod ^ target
        if diff == 0:
            total_sha += 1

        for bit in range(32):
            if not ((diff >> bit) & 1):
                preserved_sha[bit] += 1

    print(f"\n  For ACTUAL SHA-256 computation:")
    print(f"  P(all 32 bits preserved) = {total_sha/N:.6f}")
    print(f"  Per-bit: avg P = {sum(preserved_sha)/(32*N):.4f}")

    # The difference tells us about the correlation structure
    print(f"\n  Ratio (SHA vs random): {total_sha/max(total_preserved_all,1):.2f}x")
    print(f"  Per-bit ratio: {sum(preserved_sha)/max(sum(preserved_per_bit),1):.4f}x")

# ============================================================
# Part 3: Multi-round accumulation of the invariant
# ============================================================
def multi_round_invariant():
    """Track Q over many rounds. Does the bias accumulate?"""
    print("\n" + "=" * 70)
    print("Part 3: Multi-round Q preservation")
    print("=" * 70)

    random.seed(42)
    N = 100000

    # For each trial, compute Q at round 0 and round t, check preservation
    q_preserved_by_round = defaultdict(int)

    for _ in range(N):
        s = tuple(random.randint(0, MASK32) for _ in range(8))
        q0 = Q(s)

        for t in range(30):
            W = random.randint(0, MASK32)
            s = sha256_step(s, W, K[t % 64])
            qt = Q(s)

            # How many bits of Q are preserved from round 0?
            preserved_bits = 32 - hw(q0 ^ qt)
            q_preserved_by_round[t] += preserved_bits

    print(f"\n  Average preserved Q bits (out of 32) from round 0:")
    for t in range(30):
        avg = q_preserved_by_round[t] / N
        expected = 16.0  # Random expectation
        bias = avg - expected
        bar = "█" * int(abs(bias) * 20) if bias > 0 else "░" * int(abs(bias) * 20)
        print(f"  Round {t:2d}: {avg:.2f} / 32 (bias = {bias:+.2f}) {bar}")

# ============================================================
# Part 4: Search for more quadratic invariants
# ============================================================
def search_more_invariants():
    """Systematically search for other quadratic semi-invariants.

    Form: sum of (reg_i & reg_j) for various pairs (i,j).
    """
    print("\n" + "=" * 70)
    print("Part 4: Systematic search for quadratic semi-invariants")
    print("=" * 70)

    random.seed(42)
    N = 50000
    reg_names = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']

    # Test all pairs (i,j) with i < j
    print(f"\n  Single-pair invariants: P((s[i]&s[j]) XOR bit preserved)")
    results = []
    for i in range(8):
        for j in range(i+1, 8):
            preserved = 0
            for _ in range(N):
                s = tuple(random.randint(0, MASK32) for _ in range(8))
                W = random.randint(0, MASK32)
                s_new = sha256_step(s, W, K[0])

                q_old = s[i] & s[j]
                q_new = s_new[i] & s_new[j]
                preserved += (32 - hw(q_old ^ q_new))

            avg_p = preserved / (32 * N)
            bias = avg_p - 0.75  # For AND of independent bits, P(same) = 3/4
            results.append((abs(bias), i, j, avg_p, bias))

    results.sort(reverse=True)
    print(f"\n  Top biased pairs (expected P = 0.75 for independent AND):")
    for abs_bias, i, j, p, bias in results[:15]:
        marker = " ← SIGNIFICANT" if abs_bias > 0.01 else ""
        print(f"    {reg_names[i]}&{reg_names[j]}: P = {p:.4f} (bias = {bias:+.4f}){marker}")

    # Now test XOR combinations of multiple pairs
    print(f"\n  Multi-pair XOR invariants:")
    combos = [
        ("ae⊕bf⊕cg⊕dh", [(0,4),(1,5),(2,6),(3,7)]),  # Our discovery
        ("ae⊕bf⊕cg", [(0,4),(1,5),(2,6)]),
        ("ae⊕cg", [(0,4),(2,6)]),
        ("bf⊕dh", [(1,5),(3,7)]),
        ("ae⊕dh", [(0,4),(3,7)]),
        ("ab⊕cd⊕ef⊕gh", [(0,1),(2,3),(4,5),(6,7)]),
        ("ac⊕bd⊕eg⊕fh", [(0,2),(1,3),(4,6),(5,7)]),
        ("ad⊕bc⊕eh⊕fg", [(0,3),(1,2),(4,7),(5,6)]),
        ("af⊕bg⊕ch⊕de", [(0,5),(1,6),(2,7),(3,4)]),
        ("ag⊕bh⊕ce⊕df", [(0,6),(1,7),(2,4),(3,5)]),
        ("ah⊕be⊕cf⊕dg", [(0,7),(1,4),(2,5),(3,6)]),
    ]

    for name, pairs in combos:
        preserved = 0
        for _ in range(N):
            s = tuple(random.randint(0, MASK32) for _ in range(8))
            W = random.randint(0, MASK32)
            s_new = sha256_step(s, W, K[0])

            q_old = 0
            q_new = 0
            for i, j in pairs:
                q_old ^= s[i] & s[j]
                q_new ^= s_new[i] & s_new[j]

            preserved += (32 - hw(q_old ^ q_new))

        avg_p = preserved / (32 * N)
        bias = avg_p - 0.5
        marker = " ← SIGNIFICANT" if abs(bias) > 0.05 else ""
        print(f"    {name:25s}: P = {avg_p:.4f} (bias = {bias:+.4f}){marker}")

# ============================================================
# Part 5: Can we exploit the invariant for collision?
# ============================================================
def exploit_analysis():
    """If Q is approximately preserved, what does this tell us about collisions?

    For a collision: H(M) = H(M'), which means state(M) = state(M')
    after feedforward. In particular, Q(state(M)) = Q(state(M')).

    If Q is preserved with bias b > 0, then after n rounds,
    the bias decays as b^n (if independent).

    P(Q preserved after n rounds) ≈ 0.5 + 0.125 * decay_factor

    This gives us a DISTINGUISHER: the output of SHA-256 is not
    perfectly random with respect to Q.
    """
    print("\n" + "=" * 70)
    print("Part 5: Exploitation analysis")
    print("=" * 70)

    random.seed(42)
    N = 200000

    # Check: is Q(IV) → Q(output) biased for full 64-round SHA-256?
    q_match_count = 0

    for _ in range(N):
        W = [random.randint(0, MASK32) for _ in range(16)]
        s = list(tuple(0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
                       0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19))
        q_iv = Q(s)

        # Full 64-round compression
        from step18a_invariants import sha256_step as step_fn
        W_exp = list(W)
        for t in range(16, 64):
            W_exp.append(((rotr(W_exp[t-2],17)^rotr(W_exp[t-2],19)^(W_exp[t-2]>>10))
                         + W_exp[t-7]
                         + (rotr(W_exp[t-15],7)^rotr(W_exp[t-15],18)^(W_exp[t-15]>>3))
                         + W_exp[t-16]) & MASK32)

        for t in range(64):
            s = list(sha256_step(tuple(s), W_exp[t], K[t]))

        # Output hash = IV + state
        h = tuple((s[i] + [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
                           0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19][i]) & MASK32
                  for i in range(8))
        q_out = Q(h)

        # How many bits match?
        match_bits = 32 - hw(q_iv ^ q_out)
        q_match_count += match_bits

    avg_match = q_match_count / N
    expected = 16.0  # Random
    bias_per_bit = avg_match / 32 - 0.5

    print(f"\n  Full 64-round SHA-256:")
    print(f"  Q(IV) vs Q(hash): avg preserved bits = {avg_match:.2f} / 32")
    print(f"  Expected (random): 16.00")
    print(f"  Per-bit bias: {bias_per_bit:+.6f}")
    print(f"  Significant? {'YES' if abs(bias_per_bit) > 0.01 else 'NO (random)'}")

    # Check at reduced rounds
    print(f"\n  Reduced-round Q preservation:")
    for n_rounds in [1, 2, 4, 8, 16, 32, 64]:
        q_match = 0
        N2 = 100000
        for _ in range(N2):
            W = [random.randint(0, MASK32) for _ in range(16)]
            s = list(tuple(0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
                           0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19))
            q_iv = Q(s)

            W_exp = list(W)
            for t in range(16, 64):
                W_exp.append(((rotr(W_exp[t-2],17)^rotr(W_exp[t-2],19)^(W_exp[t-2]>>10))
                             + W_exp[t-7]
                             + (rotr(W_exp[t-15],7)^rotr(W_exp[t-15],18)^(W_exp[t-15]>>3))
                             + W_exp[t-16]) & MASK32)

            for t in range(min(n_rounds, 64)):
                s = list(sha256_step(tuple(s), W_exp[t], K[t]))

            h = tuple((s[i] + [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
                               0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19][i]) & MASK32
                      for i in range(8))
            q_out = Q(h)
            q_match += (32 - hw(q_iv ^ q_out))

        avg = q_match / N2
        bias = avg/32 - 0.5
        print(f"    {n_rounds:2d} rounds: avg preserved = {avg:.2f}/32 "
              f"(bias = {bias:+.6f})")


if __name__ == "__main__":
    per_bit_analysis()
    algebraic_analysis()
    multi_round_invariant()
    search_more_invariants()
    exploit_analysis()

    print("\n" + "=" * 70)
    print("Step 18b Complete")
    print("=" * 70)
