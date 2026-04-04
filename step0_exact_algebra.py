"""
Step 0: Exact Algebra of Mini-SHA

Goal: compute the EXACT algebraic normal form (ANF) of a miniature
SHA-like function. Not statistics — the actual polynomial.

Mini-SHA parameters:
  - Word size: n = 4 bits
  - Message words: 4 (= 16 input bits, 2^16 = 65536 — enumerable)
  - Registers: 8 (a,b,c,d,e,f,g,h) same shift structure as SHA-256
  - Rotations: ROTR(x, 1), ROTR(x, 2), ROTR(x, 3) (adapted for n=4)
  - Ch, Maj: same as SHA-256
  - K constants: fixed small values
  - Schedule: simplified (W[r] = 0 for r >= 4, or simple recurrence)
  - Output: 2n = 8 bits (a, e after R rounds)

This gives f: {0,1}^16 → {0,1}^8.
We compute the full truth table, then extract ANF via Möbius transform.
"""

import itertools
import numpy as np
from collections import Counter

# ============================================================
# Mini-SHA primitives (n=4 bits)
# ============================================================

N = 4  # word size in bits
MASK = (1 << N) - 1  # 0xF for n=4

def rotr(x, r):
    """Rotate right by r positions in n-bit word."""
    return ((x >> r) | (x << (N - r))) & MASK

def ch(e, f, g):
    return (e & f) ^ (~e & g) & MASK

def maj(a, b, c):
    return (a & b) ^ (a & c) ^ (b & c)

def sigma0(x):
    """Σ₀ for mini-SHA: ROTR(1) ⊕ ROTR(2) ⊕ ROTR(3)"""
    return rotr(x, 1) ^ rotr(x, 2) ^ rotr(x, 3)

def sigma1(x):
    """Σ₁ for mini-SHA: ROTR(1) ⊕ ROTR(2) ⊕ ROTR(3) — same for now"""
    # Using different shifts to avoid degeneracy
    return rotr(x, 1) ^ rotr(x, 3) ^ (x >> 1)  # last = SHR, like σ₁

def add(a, b):
    """Addition mod 2^n"""
    return (a + b) & MASK

# ============================================================
# Mini-SHA round function
# ============================================================

# Fixed IV (like SHA-256 H0 but tiny)
IV = [0x6, 0xB, 0x3, 0xA, 0x5, 0x9, 0x1, 0xF]  # 8 registers

# Fixed round constants
K = [0x4, 0x2, 0xB, 0x7, 0xA, 0x3, 0xE, 0x5,
     0x9, 0x1, 0xD, 0x6, 0x0, 0x8, 0xC, 0xF]

def mini_sha(msg_words, R):
    """
    Mini-SHA compress function.
    msg_words: list of 4 words (each n bits)
    R: number of rounds
    Returns: (a, e) after R rounds — 2*n = 8 output bits
    """
    assert len(msg_words) == 4

    # Schedule: W[0..3] = message, W[4..R-1] = 0 (simplest)
    W = list(msg_words) + [0] * max(0, R - 4)

    # State
    a, b, c, d, e, f, g, h = IV[:]

    for r in range(R):
        w_r = W[r] if r < len(W) else 0
        k_r = K[r % len(K)]

        T1 = add(add(add(add(h, sigma1(e)), ch(e, f, g)), k_r), w_r)
        T2 = add(sigma0(a), maj(a, b, c))

        # Shift
        h = g
        g = f
        f = e
        e = add(d, T1)
        d = c
        c = b
        b = a
        a = add(T1, T2)

    return (a, e)

# ============================================================
# Full truth table
# ============================================================

def compute_truth_table(R, n_msg_words=4):
    """
    Compute complete truth table of mini-SHA.
    Input: 16 bits (4 words × 4 bits)
    Output: 8 bits (a, e — 4 bits each)
    """
    n_input = n_msg_words * N  # 16
    n_total = 1 << n_input      # 65536

    # output_bits[i][j] = j-th output bit for input i
    output_bits = np.zeros((n_total, 2 * N), dtype=np.uint8)

    for idx in range(n_total):
        # Decode index into message words
        words = []
        tmp = idx
        for w in range(n_msg_words):
            words.append(tmp & MASK)
            tmp >>= N

        a_out, e_out = mini_sha(words, R)

        # Extract output bits
        for b in range(N):
            output_bits[idx][b] = (a_out >> b) & 1
            output_bits[idx][N + b] = (e_out >> b) & 1

    return output_bits

# ============================================================
# ANF via Möbius transform
# ============================================================

def mobius_transform(truth_table, n_vars):
    """
    Compute ANF from truth table via Möbius transform.
    truth_table: array of 2^n_vars bits (0 or 1)
    Returns: ANF coefficients (same size array)

    anf[S] = 1 iff monomial ∏_{i∈S} x_i appears in the polynomial.
    """
    n = 1 << n_vars
    anf = truth_table.copy()

    for i in range(n_vars):
        step = 1 << i
        for j in range(n):
            if j & step:
                anf[j] ^= anf[j ^ step]

    return anf

def anf_to_monomials(anf, n_vars):
    """Convert ANF coefficient array to list of monomials (as sets of variable indices)."""
    monomials = []
    for idx in range(len(anf)):
        if anf[idx]:
            # Extract which variables are in this monomial
            variables = []
            for v in range(n_vars):
                if idx & (1 << v):
                    variables.append(v)
            monomials.append(tuple(variables))
    return monomials

def analyze_anf(monomials, n_vars):
    """Analyze ANF structure."""
    if not monomials:
        return {"degree": 0, "n_monomials": 0}

    degrees = [len(m) for m in monomials]
    degree_dist = Counter(degrees)

    # Which variables appear in highest-degree monomials
    max_deg = max(degrees)
    top_monomials = [m for m in monomials if len(m) == max_deg]

    # Variable frequency across ALL monomials
    var_freq = Counter()
    for m in monomials:
        for v in m:
            var_freq[v] += 1

    return {
        "degree": max_deg,
        "n_monomials": len(monomials),
        "degree_distribution": dict(degree_dist),
        "top_monomials_count": len(top_monomials),
        "var_frequency": dict(var_freq),
    }

# ============================================================
# Main experiment
# ============================================================

def var_name(idx):
    """Human-readable variable name."""
    word = idx // N
    bit = idx % N
    return f"M[{word}]_{bit}"

def run_experiment():
    n_msg = 4
    n_input = n_msg * N  # 16

    print(f"Mini-SHA: n={N}, msg_words={n_msg}, input_bits={n_input}")
    print(f"Truth table size: 2^{n_input} = {1 << n_input}")
    print()

    for R in [1, 2, 3, 4, 5, 6, 8]:
        print(f"{'='*60}")
        print(f"R = {R} rounds")
        print(f"{'='*60}")

        tt = compute_truth_table(R, n_msg)

        for out_bit in range(2 * N):  # 8 output bits
            reg = "a" if out_bit < N else "e"
            bit_idx = out_bit % N

            # Extract truth table for this output bit
            tt_bit = tt[:, out_bit].copy()

            # Compute ANF
            anf = mobius_transform(tt_bit, n_input)
            monomials = anf_to_monomials(anf, n_input)
            analysis = analyze_anf(monomials, n_input)

            print(f"  {reg}[{bit_idx}]: degree={analysis['degree']}, "
                  f"monomials={analysis['n_monomials']}/{1 << n_input}, "
                  f"degree_dist={analysis['degree_distribution']}")

        # Summary for this R
        all_degrees = []
        total_monomials = []
        for out_bit in range(2 * N):
            tt_bit = tt[:, out_bit].copy()
            anf = mobius_transform(tt_bit, n_input)
            monomials = anf_to_monomials(anf, n_input)
            analysis = analyze_anf(monomials, n_input)
            all_degrees.append(analysis['degree'])
            total_monomials.append(analysis['n_monomials'])

        print(f"\n  SUMMARY R={R}:")
        print(f"    Max degree: {max(all_degrees)}")
        print(f"    Min degree: {min(all_degrees)}")
        print(f"    Avg monomials: {sum(total_monomials)/len(total_monomials):.0f} / {1<<n_input}")
        print(f"    Density: {sum(total_monomials)/len(total_monomials)/(1<<n_input)*100:.1f}%")
        print()

if __name__ == "__main__":
    run_experiment()
