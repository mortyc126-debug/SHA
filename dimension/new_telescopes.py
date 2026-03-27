"""
ТРИ НОВЫХ ТЕЛЕСКОПА: ищем structure где Hamming видит random.

Телескоп 1: МОДУЛЯРНАЯ АРИФМЕТИКА (Z/2^32)
  SHA-256 использует ADD, не XOR. Арифметическая structure может
  сохраняться там где XOR-structure стёрта.
  Метрика: |H1[i] - H2[i]| mod 2^32 вместо HW(H1[i] ⊕ H2[i])

Телескоп 2: АЛГЕБРАИЧЕСКАЯ ГЕОМЕТРИЯ
  State = точка на variety. Ищем АЛГЕБРАИЧЕСКИЕ СООТНОШЕНИЯ
  между регистрами: f(a,b,c,d,e,f,g,h) = 0.

Телескоп 3: 2-АДИЧЕСКИЙ АНАЛИЗ
  v_2(x) = число trailing zeros (2-adic valuation).
  Carry в ADD связан с 2-adic structure.
  Ищем: v_2 инварианты через раунды.
"""

import numpy as np
import struct, hashlib

MASK32 = 0xFFFFFFFF

def hw(x): return bin(x).count('1')
def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def add32(x, y): return (x + y) & MASK32
def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)

K = [
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
]
IV = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]


def sha256_full(W16):
    raw = struct.pack('>16I', *W16)
    return struct.unpack('>8I', hashlib.sha256(raw).digest())


def v2(x):
    """2-adic valuation: number of trailing zeros."""
    if x == 0: return 32
    v = 0
    while (x >> v) & 1 == 0: v += 1
    return v


def arith_dist(a, b):
    """Arithmetic distance in Z/2^32."""
    return min((a - b) & MASK32, (b - a) & MASK32)


def main():
    np.random.seed(42)

    print("=" * 70)
    print("ТРИ НОВЫХ ТЕЛЕСКОПА")
    print("=" * 70)

    # ═══════════════════════════
    print(f"\n{'='*70}")
    print("ТЕЛЕСКОП 1: МОДУЛЯРНАЯ АРИФМЕТИКА")
    print(f"{'='*70}")

    # Arithmetic distance vs Hamming distance
    # For random function: arith_dist ≈ 2^31/2 ≈ 10^9, uncorrelated with Hamming
    # For SHA-256: ADD structure might create correlation

    # Test: δW = +1 (arithmetic) vs δW = XOR bit 0 (same when bit0=0)
    print(f"\n  A1. Arithmetic δH vs Hamming δH (1-bit flip, 10K pairs):")

    hamming_diffs = []
    arith_diffs = []

    for _ in range(10000):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W); W2[0] ^= 1  # XOR flip bit 0
        H1 = sha256_full(W); H2 = sha256_full(W2)

        hamming_diffs.append(sum(hw(H1[i]^H2[i]) for i in range(8)))
        arith_diffs.append(sum(arith_dist(H1[i], H2[i]) for i in range(8)))

    corr = np.corrcoef(hamming_diffs, arith_diffs)[0, 1]
    print(f"    Corr(Hamming, Arith): {corr:.4f}")
    print(f"    Hamming: mean={np.mean(hamming_diffs):.1f}")
    print(f"    Arith: mean={np.mean(arith_diffs):.0f}")

    # Test: does arithmetic distance reveal structure at r=24 that Hamming misses?
    print(f"\n  A2. Arithmetic near-collision search (full SHA-256):")

    # Search for pairs with LOW arithmetic distance
    best_arith = float('inf')
    best_hamming = 256

    for _ in range(50000):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W); W2[15] ^= (1 << 31)
        H1 = sha256_full(W); H2 = sha256_full(W2)

        ad = sum(arith_dist(H1[i], H2[i]) for i in range(8))
        hd = sum(hw(H1[i]^H2[i]) for i in range(8))

        if ad < best_arith: best_arith = ad
        if hd < best_hamming: best_hamming = hd

    # Expected arithmetic distance: 8 × 2^31/3 ≈ 5.7×10^9
    print(f"    Best Hamming: {best_hamming}")
    print(f"    Best Arith: {best_arith:,.0f}")
    print(f"    Expected Arith (random): ~{8 * 2**31 // 3:,}")

    # Key test: is the pair with best arith DIFFERENT from best hamming?
    # If different → arithmetic reveals different structure

    # ═══════════════════════════
    print(f"\n{'='*70}")
    print("ТЕЛЕСКОП 2: АЛГЕБРАИЧЕСКИЕ СООТНОШЕНИЯ")
    print(f"{'='*70}")

    # Look for polynomial relations: f(H[0], H[1], ..., H[7]) = 0
    # Over Z/2^32 or modular arithmetic

    # Test: does H[0] + H[1] + ... + H[7] mod 2^k have bias?
    print(f"\n  B1. Sum of hash words mod small p:")

    for p in [2, 3, 4, 5, 7, 8, 16, 256]:
        counts = {}
        N = 50000
        for _ in range(N):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            H = sha256_full(W)
            s = 0
            for h in H: s = (s + h) & MASK32
            r = s % p
            counts[r] = counts.get(r, 0) + 1

        expected = N / p
        chi2 = sum((counts.get(i, 0) - expected)**2 / expected for i in range(p))
        pval_approx = "BIAS ★" if chi2 > p * 3 else "uniform"
        print(f"    mod {p:>3}: χ²={chi2:>7.1f} (threshold={p*3:>5}) → {pval_approx}")

    # Test: pairwise relations H[i]*H[j] mod p
    print(f"\n  B2. Product relations H[0]*H[4] mod p:")
    for p in [2, 3, 5, 7]:
        counts = {}
        N = 50000
        for _ in range(N):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            H = sha256_full(W)
            r = (H[0] * H[4]) % p
            counts[r] = counts.get(r, 0) + 1

        expected = N / p
        chi2 = sum((counts.get(i, 0) - expected)**2 / expected for i in range(p))
        print(f"    H[0]*H[4] mod {p}: χ²={chi2:.1f} → {'BIAS ★' if chi2 > p*3 else 'uniform'}")

    # Test: does (H[0] - H[4]) mod p have bias? (conservation remnant?)
    print(f"\n  B3. Difference H[0]-H[4] mod p (conservation remnant?):")
    for p in [2, 3, 4, 8, 16, 256]:
        counts = {}
        N = 50000
        for _ in range(N):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            H = sha256_full(W)
            r = ((H[0] - H[4]) & MASK32) % p
            counts[r] = counts.get(r, 0) + 1

        expected = N / p
        chi2 = sum((counts.get(i, 0) - expected)**2 / expected for i in range(p))
        print(f"    (H[0]-H[4]) mod {p:>3}: χ²={chi2:>7.1f} → {'BIAS ★' if chi2 > p*3 else 'uniform'}")

    # ═══════════════════════════
    print(f"\n{'='*70}")
    print("ТЕЛЕСКОП 3: 2-АДИЧЕСКИЙ АНАЛИЗ")
    print(f"{'='*70}")

    # 2-adic valuation v2(x) = trailing zeros
    # For random 32-bit: E[v2] = 1 (geometric distribution)
    # For SHA-256 output: same or different?

    print(f"\n  C1. Distribution of v2(H[i]) for SHA-256 output:")

    v2_counts = {i: 0 for i in range(33)}
    N = 100000

    for _ in range(N):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H = sha256_full(W)
        for h in H:
            v2_counts[v2(h)] += 1

    total = N * 8
    print(f"    {'v2':>3} {'Count':>8} {'Rate':>8} {'Expected':>8} {'Ratio':>6}")
    for v in range(10):
        expected_rate = 1 / (2 ** (v + 1))
        actual_rate = v2_counts[v] / total
        ratio = actual_rate / expected_rate if expected_rate > 0 else 0
        marker = " ★" if abs(ratio - 1) > 0.03 else ""
        print(f"    {v:>3} {v2_counts[v]:>8} {actual_rate:>8.4f} {expected_rate:>8.4f} {ratio:>5.3f}{marker}")

    # Test: v2(H[0] + H[4]) — does conservation leave 2-adic trace?
    print(f"\n  C2. v2(H[0] + H[4]) — conservation remnant?")
    v2_sum_counts = {i: 0 for i in range(33)}
    for _ in range(N):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H = sha256_full(W)
        s = add32(H[0], H[4])
        v2_sum_counts[v2(s)] += 1

    print(f"    {'v2':>3} {'Count':>8} {'Rate':>8} {'Expected':>8} {'Ratio':>6}")
    for v in range(10):
        expected_rate = 1 / (2 ** (v + 1))
        actual_rate = v2_sum_counts[v] / N
        ratio = actual_rate / expected_rate if expected_rate > 0 else 0
        marker = " ★" if abs(ratio - 1) > 0.05 else ""
        print(f"    {v:>3} {v2_sum_counts[v]:>8} {actual_rate:>8.4f} {expected_rate:>8.4f} {ratio:>5.3f}{marker}")

    # Test: v2(δH) when flipping 1 bit — does carry create 2-adic structure?
    print(f"\n  C3. v2(H1[0] ⊕ H2[0]) for 1-bit flip — carry structure?")
    v2_xor_counts = {i: 0 for i in range(33)}
    v2_arith_counts = {i: 0 for i in range(33)}

    for _ in range(N):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W); W2[0] ^= 1
        H1 = sha256_full(W); H2 = sha256_full(W2)
        xor_diff = H1[0] ^ H2[0]
        arith_diff = (H1[0] - H2[0]) & MASK32
        v2_xor_counts[v2(xor_diff)] += 1
        v2_arith_counts[v2(arith_diff)] += 1

    print(f"    {'v2':>3} {'XOR rate':>9} {'ARITH rate':>11} {'Expected':>9} {'XOR ratio':>10} {'ARITH ratio':>12}")
    for v in range(10):
        expected_rate = 1 / (2 ** (v + 1))
        xor_rate = v2_xor_counts[v] / N
        arith_rate = v2_arith_counts[v] / N
        xor_r = xor_rate / expected_rate if expected_rate > 0 else 0
        arith_r = arith_rate / expected_rate if expected_rate > 0 else 0
        xm = " ★" if abs(xor_r - 1) > 0.05 else ""
        am = " ★" if abs(arith_r - 1) > 0.05 else ""
        print(f"    {v:>3} {xor_rate:>9.4f} {arith_rate:>11.4f} {expected_rate:>9.4f} {xor_r:>9.3f}{xm} {arith_r:>11.3f}{am}")

    # ═══════════════════════════
    print(f"\n{'='*70}")
    print("SYNTHESIS")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
