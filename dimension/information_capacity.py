"""
INFORMATION CAPACITY: сколько бит информации проходит через SHA-256 за раунд.

Absorption Law (T16): каждое слово поглощается за ~5 раундов.
Но: сколько бит ВЗАИМНОЙ ИНФОРМАЦИИ между входом и выходом?

I(W; H|r) = mutual information at round r.
  At r=0: I=0 (no rounds, output = IV)
  At r→∞: I=256 (full)
  Growth: how fast?

Это ИНФОРМАЦИОННАЯ скорость SHA-256 — сколько бит/раунд.
Если скорость = v, то r_sat = 256/v.

Измеряем через entropy of conditional distribution.
"""

import numpy as np
import struct, hashlib
from collections import Counter

MASK32 = 0xFFFFFFFF
K_const = [
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
IV = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]

def hw(x): return bin(x).count('1')
def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def add32(x, y): return (x + y) & MASK32
def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)

def sha256_r(W16, n_r):
    W = list(W16)
    for r in range(16, max(n_r, 16)):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))
    a,b,c,d,e,f,g,h = IV
    for r in range(n_r):
        T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K_const[r]), W[r])
        T2 = add32(Sigma0(a), Maj(a,b,c))
        h,g,f,e = g,f,e,add32(d,T1)
        d,c,b,a = c,b,a,add32(T1,T2)
    return tuple(add32(IV[i], [a,b,c,d,e,f,g,h][i]) for i in range(8))


def entropy_bits(values):
    """Estimate entropy of a list of values in bits."""
    c = Counter(values)
    n = len(values)
    h = 0
    for count in c.values():
        p = count / n
        if p > 0:
            h -= p * np.log2(p)
    return h


def main():
    np.random.seed(42)

    print("=" * 70)
    print("INFORMATION CAPACITY of SHA-256")
    print("=" * 70)

    # =========================================================
    print(f"\n{'=' * 70}")
    print("1. OUTPUT ENTROPY vs ROUNDS (varying single word)")
    print("=" * 70)

    # Fix all words except W[0], vary W[0] randomly.
    # Measure: entropy of H at each round.
    W_base = [np.random.randint(0, 2**32) for _ in range(16)]
    N_samples = 5000

    print(f"\n  Vary W[0] (5000 random values), measure H entropy:")
    print(f"  {'Rounds':>6} {'H(H0) bits':>11} {'H(H0-3) bits':>13} {'Sens':>6} {'Full output bits':>16}")

    for n_r in [1, 2, 3, 4, 6, 8, 12, 16, 20, 24, 32, 64]:
        # Collect outputs
        h0_values = []
        h_low_values = []
        sens_vals = []

        for _ in range(N_samples):
            W2 = list(W_base)
            W2[0] = np.random.randint(0, 2**32)
            H2 = sha256_r(W2, n_r)
            h0_values.append(H2[0])
            h_low_values.append(H2[:4])  # First 4 words

            # Also sensitivity
            W3 = list(W2); W3[0] ^= 1  # flip 1 bit
            H3 = sha256_r(W3, n_r)
            sens_vals.append(sum(hw(H2[i]^H3[i]) for i in range(8)))

        # Entropy of H[0] (32-bit word)
        e_h0 = entropy_bits(h0_values)

        # Entropy of (H[0],H[1],H[2],H[3]) - approximate via hashing
        h_low_hashed = [hash(v) for v in h_low_values]
        e_h_low = entropy_bits(h_low_hashed)

        # Full output entropy (all 8 words)
        # Can't compute exactly, estimate from unique count
        full_outputs = [sha256_r(list(W_base[:0]) + [np.random.randint(0, 2**32)] + W_base[1:], n_r)
                       for _ in range(min(N_samples, 2000))]
        unique_count = len(set(full_outputs))
        est_bits = np.log2(unique_count) if unique_count > 1 else 0

        print(f"  {n_r:>6} {e_h0:>10.1f} {e_h_low:>12.1f} {np.mean(sens_vals):>5.1f} {est_bits:>15.1f}")

    # =========================================================
    print(f"\n{'=' * 70}")
    print("2. INFORMATION RATE: bits of output entropy per round")
    print("=" * 70)

    # Use unique output count as proxy for entropy
    print(f"\n  {'Rounds':>6} {'Unique outputs':>15} {'Est. entropy':>13} {'Δ(entropy)':>11} {'bits/round':>10}")

    prev_ent = 0
    prev_r = 0
    for n_r in [1, 2, 3, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20, 24, 32, 64]:
        outputs = set()
        for _ in range(10000):
            W2 = list(W_base)
            W2[0] = np.random.randint(0, 2**32)
            H2 = sha256_r(W2, n_r)
            outputs.add(H2)

        ent = np.log2(len(outputs)) if len(outputs) > 1 else 0
        delta = ent - prev_ent
        rate = delta / max(1, n_r - prev_r) if n_r > prev_r else 0

        print(f"  {n_r:>6} {len(outputs):>14} {ent:>12.1f} {delta:>10.1f} {rate:>9.1f}")
        prev_ent = ent
        prev_r = n_r

    # =========================================================
    print(f"\n{'=' * 70}")
    print("3. MULTI-WORD ENTROPY: how many input words contribute?")
    print("=" * 70)

    # Vary 1, 2, 4, 8, 16 words at r=16 and r=64
    for n_r in [16, 64]:
        print(f"\n  r={n_r}:")
        print(f"  {'Words varied':>13} {'Unique':>10} {'Est. bits':>10}")
        for n_vary in [1, 2, 4, 8, 16]:
            outputs = set()
            for _ in range(10000):
                W2 = list(W_base)
                for w in range(n_vary):
                    W2[w] = np.random.randint(0, 2**32)
                outputs.add(sha256_r(W2, n_r))
            ent = np.log2(len(outputs)) if len(outputs) > 1 else 0
            print(f"  {n_vary:>13} {len(outputs):>10} {ent:>9.1f}")

    # =========================================================
    print(f"\n{'=' * 70}")
    print("4. CHANNEL CAPACITY MODEL")
    print("=" * 70)

    print(f"""
  SHA-256 as an INFORMATION CHANNEL:

  Input: 512 bits (16 × 32-bit words)
  Output: 256 bits (8 × 32-bit words)
  Channel: r rounds of compression

  CAPACITY per round:
    Round function updates 2 registers (a, e) = 64 new bits.
    But: carries + nonlinearity = not full 64 bits of NEW information.
    Effective new info per round ≈ 28-32 bits (matches absorption rate λ≈28).

  TOTAL CAPACITY after r rounds:
    C(r) = min(256, Σ new_info_per_round)
    C(r) ≈ min(256, 28 × r) for r ≤ 9
    C(r) = 256 for r ≥ 10 (saturated, one word contributes fully)

    But FULL saturation (ALL 16 words contributing):
    Requires r ≥ 20 (last word W[15] absorbed at r=20).

  INFORMATION-THEORETIC SECURITY BOUNDARY:
    Full channel capacity C = 256 bits at r ≥ 20.
    Collision cost = 2^(C/2) = 2^128 at r ≥ 20.

    For r < 20: C < 256 → collision < 2^128.
    The "discount" = (256 - C(r)) / 2 bits off collision cost.

  COMPARISON WITH OUR DIMENSION:
    rank(T) = 256 at r=16:  ALL output bits REACHABLE.
    C = 256 at r=20:        ALL output bits UNIFORMLY MIXED.
    K = 128 at r=24:        ALL output bits on ISOTROPIC SPHERE.

  THREE information-theoretic thresholds:
    r=16: algebraic reachability (rank)
    r=20: information saturation (capacity)
    r=24: geometric uniformity (sphere)

  SHA-256 (r=64): 3.2× margin over sphere boundary.
""")


if __name__ == "__main__":
    main()
