"""
ЗАДАЧА A1: Алгебраическая степень SHA-256.

f: GF(2)^512 → GF(2) — один выходной бит.
deg(f) = ? как многочлен Жегалкина (ANF) над GF(2).

Round function: Ch(e,f,g) = ef ⊕ ēg = ef ⊕ g ⊕ eg — degree 2.
                Maj(a,b,c) = ab ⊕ ac ⊕ bc — degree 2.
                ADD mod 2^32: degree ЗАВИСИТ ОТ БИТА.
                  bit 0 of (x+y): x₀ ⊕ y₀ — degree 1.
                  bit 1: x₁ ⊕ y₁ ⊕ x₀y₀ — degree 2.
                  bit k: degree k+1 (carry chain).
                  bit 31: degree 32!

ADD — ГЛАВНЫЙ источник алгебраической степени!
Ch, Maj: degree 2 per round.
ADD: degree UP TO 32 per operation, 7 ADDs per round.

Вопрос: растёт ли степень как 2^r (наивно) или МЕДЛЕННЕЕ?
Если медленнее — есть шанс что deg < 128 при 64 раундах.

МЕТОДЫ ИЗМЕРЕНИЯ:
  1. CUBE ATTACK: order-k differential = 0 iff deg < k.
  2. DIRECT: evaluate ANF coefficient для monomial степени d.
  3. PROBABILISTIC: random higher-order diff, measure zero rate.
"""

import numpy as np
import struct, hashlib
import time

MASK32 = 0xFFFFFFFF

def sha256_bit(W16, out_word, out_bit, n_r=64):
    """Return single output bit of SHA-256."""
    raw = struct.pack('>16I', *W16)
    if n_r == 64:
        H = struct.unpack('>8I', hashlib.sha256(raw).digest())
    else:
        from functools import reduce
        def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
        def add32(x, y): return (x + y) & MASK32
        def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
        def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)
        def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
        def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
        def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
        def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
        K = [0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
             0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
             0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
             0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
             0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
             0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
             0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
             0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2]
        IV = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]
        W = list(W16)
        for r in range(16, max(n_r, 16)):
            W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))
        a,b,c,d,e,f,g,h = IV
        for r in range(n_r):
            T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K[r]), W[r])
            T2 = add32(Sigma0(a), Maj(a,b,c))
            h,g,f,e = g,f,e,add32(d,T1)
            d,c,b,a = c,b,a,add32(T1,T2)
        H = tuple(add32(IV[i], [a,b,c,d,e,f,g,h][i]) for i in range(8))
    return (H[out_word] >> out_bit) & 1


def higher_order_test(n_r, order_k, out_word, out_bit, N_tests=1000):
    """Test if degree < order_k by computing order-k differential."""
    zeros = 0
    for _ in range(N_tests):
        W = [np.random.randint(0, 2**32) for _ in range(16)]

        # Choose k DISTINCT random input bits
        bits = set()
        while len(bits) < order_k:
            bits.add((np.random.randint(0, 16), np.random.randint(0, 32)))
        bits = list(bits)

        # Compute 2^k XOR
        xor_sum = 0
        for mask in range(1 << order_k):
            Wm = list(W)
            for j in range(order_k):
                if (mask >> j) & 1:
                    Wm[bits[j][0]] ^= (1 << bits[j][1])
            xor_sum ^= sha256_bit(Wm, out_word, out_bit, n_r)

        if xor_sum == 0:
            zeros += 1

    return zeros / N_tests


def main():
    np.random.seed(42)

    print("=" * 70)
    print("ЗАДАЧА A1: Алгебраическая степень SHA-256")
    print("=" * 70)

    # ═══════════════════
    print(f"\n{'='*70}")
    print("1. ТЕОРИЯ: degree of ADD mod 2^32")
    print(f"{'='*70}")

    print(f"""
  ADD mod 2^32: z = x + y.
  z₀ = x₀ ⊕ y₀                    (degree 1)
  z₁ = x₁ ⊕ y₁ ⊕ x₀y₀            (degree 2)
  z₂ = x₂ ⊕ y₂ ⊕ (x₁⊕y₁)x₀y₀ ⊕ x₁y₁  (degree 3)
  ...
  z_k = degree k+1 (carry chain adds 1 degree per bit)

  Per round: T1 = h + Σ1(e) + Ch + K + W — 4 ADDs chained.
  Degree: bit 0 of T1 = degree 1. Bit 31 of T1 = degree 32.
  But: T1 feeds into a = T1 + T2 and e = d + T1.
  a bit 31: degree up to 32 + 32 = 64?? No — degree of SUM is
  max(deg(x), deg(y)) + 1 per carry bit.

  ACTUAL ANALYSIS:
  Round 1 input: degree 1 (input bits).
  After Ch: degree 2 (ef ⊕ g ⊕ eg).
  After ADD chain (4 adds): bit 31 → degree up to 2 + 31 = 33?
  But degree doesn't simply add — it's max over carry chain.

  KEY INSIGHT: carry chain degree = degree of OPERANDS + chain length.
  If operands have degree d, and we add k of them:
  Result bit j: degree ≤ d + j (from carry propagation).

  After r rounds with feedback:
  Degree growth: NOT 2^r (that's for composition of degree-2 maps).
  Actual: degree(output bit j, round r) ≈ f(r, j).
""")

    # ═══════════════════
    print(f"{'='*70}")
    print("2. MEASURE: higher-order differential at SINGLE BIT")
    print(f"{'='*70}")

    # For single output bit: test if degree < k
    # If order-k diff = 0 with P=1 → degree < k.
    # If order-k diff = 0 with P=0.5 → degree ≥ k.
    # P(zero | degree ≥ k) = 0.5 exactly.

    print(f"\n  Output bit H[0] bit 0 (LSB, lowest degree expected):")
    print(f"  {'Rounds':>6} {'Order':>6} {'P(zero)':>8} {'Degree<k?':>10}")

    for n_r in [4, 8, 12, 16, 20, 24]:
        for k in [2, 3, 4, 5, 8]:
            if k > 5 and n_r > 16:
                continue  # too slow
            N = 500 if k <= 5 else 200
            p_zero = higher_order_test(n_r, k, 0, 0, N)
            is_below = "YES ★" if p_zero > 0.9 else "maybe" if p_zero > 0.6 else "NO"
            print(f"  {n_r:>6} {k:>6} {p_zero:>8.3f} {is_below:>10}")

    print(f"\n  Output bit H[0] bit 31 (MSB, highest degree expected):")
    for n_r in [4, 8, 12, 16, 20]:
        for k in [2, 3, 4, 5]:
            N = 500
            p_zero = higher_order_test(n_r, k, 0, 31, N)
            is_below = "YES ★" if p_zero > 0.9 else "maybe" if p_zero > 0.6 else "NO"
            print(f"  {n_r:>6} {k:>6} {p_zero:>8.3f} {is_below:>10}")

    # ═══════════════════
    print(f"\n{'='*70}")
    print("3. DEGREE MAP: which (round, bit) combinations have low degree?")
    print(f"{'='*70}")

    # For each output bit position: find the ORDER at which P(zero) drops to 0.5
    print(f"\n  Degree estimation per output bit (H[0]):")
    print(f"  {'Bit':>4} {'r=4':>6} {'r=8':>6} {'r=12':>6} {'r=16':>6} {'r=20':>6}")

    for out_bit in [0, 1, 4, 8, 15, 16, 24, 31]:
        row = f"  {out_bit:>4}"
        for n_r in [4, 8, 12, 16, 20]:
            # Find degree: smallest k where P(zero) < 0.6
            deg = "≥8"
            for k in [2, 3, 4, 5, 6, 7, 8]:
                p = higher_order_test(n_r, k, 0, out_bit, 300)
                if p < 0.6:
                    deg = str(k - 1)
                    break
            row += f" {deg:>6}"
        print(row)

    # ═══════════════════
    print(f"\n{'='*70}")
    print("4. FULL SHA-256 (64 rounds): degree estimate")
    print(f"{'='*70}")

    # At 64 rounds: we expect degree to be VERY high.
    # Test: order-2 through order-8 on H[0] bit 0.
    print(f"\n  H[0] bit 0, 64 rounds:")
    for k in [2, 3, 4, 5]:
        N = 300
        p = higher_order_test(64, k, 0, 0, N)
        print(f"    Order-{k}: P(zero) = {p:.3f} {'→ deg<{}'.format(k) if p > 0.6 else '→ deg≥{}'.format(k)}")

    print(f"\n  H[0] bit 31, 64 rounds:")
    for k in [2, 3, 4, 5]:
        N = 300
        p = higher_order_test(64, k, 0, 31, N)
        print(f"    Order-{k}: P(zero) = {p:.3f} {'→ deg<{}'.format(k) if p > 0.6 else '→ deg≥{}'.format(k)}")

    # ═══════════════════
    print(f"\n{'='*70}")
    print("5. CUBE ATTACK PERSPECTIVE")
    print(f"{'='*70}")

    # Cube attack: fix some input bits ("cube variables"),
    # sum output over all 2^k combinations of cube variables.
    # If sum = constant → found degree-k relation.
    # If sum = linear in other variables → found linear equation.

    # Test: fix all bits except k cube variables in W[15].
    # Sum over 2^k combinations. Is result constant?

    print(f"\n  Cube test: k bits from W[15] as cube variables, r=24:")
    for k in [2, 3, 4, 5]:
        if k > 4:
            N_base = 50
        else:
            N_base = 200

        constant_count = 0
        for _ in range(N_base):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            cube_bits = list(np.random.choice(32, k, replace=False))

            xor_sum = 0
            for mask in range(1 << k):
                Wm = list(W)
                for j in range(k):
                    if (mask >> j) & 1:
                        Wm[15] ^= (1 << cube_bits[j])
                xor_sum ^= sha256_bit(Wm, 0, 0, 24)

            if xor_sum == 0:
                constant_count += 1

        p = constant_count / N_base
        print(f"    k={k}: P(sum=0) = {p:.3f} (0.5 = random, >0.9 = cube found)")

    # ═══════════════════
    print(f"\n{'='*70}")
    print("ВЕРДИКТ")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
