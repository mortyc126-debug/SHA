"""
КВАДРАТИЧНЫЙ ПОИСК: GF2-kernel pairs → near-collision HW=102.

Random pairs: HW(δH) ≈ 128.
GF2-kernel pairs: HW(δH) ≈ 102 (из 1225 pairs).
Compression: 26 бит.

Вопрос: масштабируется ли? Больше pairs → ещё ниже HW?
И: для r=64 то же самое?
"""

import numpy as np
import struct, hashlib

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

def sha256_full(W16):
    return struct.unpack('>8I', hashlib.sha256(struct.pack('>16I', *W16)).digest())


def build_gf2_kernel(W_base):
    """Build GF(2) kernel for full 64-round SHA-256."""
    H_base = sha256_full(W_base)

    T_rows = []
    for word in range(16):
        for bit in range(32):
            W_mod = list(W_base)
            W_mod[word] ^= (1 << bit)
            H_mod = sha256_full(W_mod)
            row = []
            for w in range(8):
                d = H_base[w] ^ H_mod[w]
                for b in range(32): row.append((d >> b) & 1)
            T_rows.append(row)

    T = np.array(T_rows, dtype=np.uint8)

    M = T.T.copy()
    pivots = []; row = 0
    for col in range(512):
        found = False
        for r in range(row, 256):
            if M[r, col] == 1:
                M[[row, r]] = M[[r, row]]; found = True; break
        if not found: continue
        pivots.append(col)
        for r in range(256):
            if r != row and M[r, col] == 1: M[r] ^= M[row]
        row += 1

    free_vars = [c for c in range(512) if c not in pivots]
    return M, pivots, free_vars, H_base


def make_kv(M, pivots, free_vars, fc_idx):
    x = np.zeros(512, dtype=np.uint8)
    x[free_vars[fc_idx]] = 1
    for i in range(len(pivots)-1, -1, -1):
        pc = pivots[i]; val = np.uint8(0)
        for j in range(512):
            if j != pc: val ^= (M[i,j] & x[j])
        x[pc] = val
    return x

def kv_to_dW(x):
    dW = [0]*16
    for word in range(16):
        for bit in range(32):
            if x[word*32+bit]: dW[word] ^= (1<<bit)
    return dW


def main():
    np.random.seed(42)

    print("=" * 70)
    print("КВАДРАТИЧНЫЙ ПОИСК: GF2-kernel pairs для FULL SHA-256")
    print("=" * 70)

    W_base = [np.random.randint(0, 2**32) for _ in range(16)]

    print(f"\n  Building GF(2) kernel for 64-round SHA-256...")
    M, pivots, free_vars, H_base = build_gf2_kernel(W_base)
    print(f"  Kernel dimension: {len(free_vars)}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("1. SINGLE kernel vectors: CE(v)")
    print("=" * 70)

    single_hws = []
    for i in range(min(100, len(free_vars))):
        x = make_kv(M, pivots, free_vars, i)
        dW = kv_to_dW(x)
        W2 = [W_base[k]^dW[k] for k in range(16)]
        H2 = sha256_full(W2)
        dh = sum(hw(H_base[k]^H2[k]) for k in range(8))
        single_hws.append(dh)

    print(f"\n  100 single GF2-kernel vectors → HW(δH):")
    print(f"    Mean: {np.mean(single_hws):.1f}")
    print(f"    Min:  {min(single_hws)}")
    print(f"    Std:  {np.std(single_hws):.1f}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("2. PAIRS of kernel vectors: CE(v1⊕v2)")
    print("=" * 70)

    pair_hws = []
    best_hw = 256

    N_pairs = 5000
    for trial in range(N_pairs):
        i = np.random.randint(0, min(256, len(free_vars)))
        j = np.random.randint(0, min(256, len(free_vars)))
        if i == j: continue

        v1 = make_kv(M, pivots, free_vars, i)
        v2 = make_kv(M, pivots, free_vars, j)
        v12 = v1 ^ v2

        dW = kv_to_dW(v12)
        if all(w == 0 for w in dW): continue

        W2 = [W_base[k]^dW[k] for k in range(16)]
        H2 = sha256_full(W2)
        dh = sum(hw(H_base[k]^H2[k]) for k in range(8))
        pair_hws.append(dh)

        if dh < best_hw:
            best_hw = dh

    print(f"\n  {len(pair_hws)} pairs → HW(δH):")
    print(f"    Mean: {np.mean(pair_hws):.1f}")
    print(f"    Min:  {min(pair_hws)}")
    print(f"    Std:  {np.std(pair_hws):.1f}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("3. TRIPLES of kernel vectors: CE(v1⊕v2⊕v3)")
    print("=" * 70)

    triple_hws = []
    for trial in range(3000):
        ids = np.random.choice(min(256, len(free_vars)), 3, replace=False)
        v = make_kv(M, pivots, free_vars, ids[0])
        for k in range(1, 3):
            v = v ^ make_kv(M, pivots, free_vars, ids[k])

        dW = kv_to_dW(v)
        if all(w == 0 for w in dW): continue

        W2 = [W_base[k]^dW[k] for k in range(16)]
        H2 = sha256_full(W2)
        dh = sum(hw(H_base[k]^H2[k]) for k in range(8))
        triple_hws.append(dh)

    print(f"\n  {len(triple_hws)} triples → HW(δH):")
    print(f"    Mean: {np.mean(triple_hws):.1f}")
    print(f"    Min:  {min(triple_hws)}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("4. RANDOM XOR (K random kernel vectors)")
    print("=" * 70)

    for K in [1, 2, 3, 5, 10, 20, 50]:
        khws = []
        for trial in range(1000):
            ids = np.random.choice(min(256, len(free_vars)), K, replace=False)
            v = np.zeros(512, dtype=np.uint8)
            for k in ids:
                v = v ^ make_kv(M, pivots, free_vars, k)

            dW = kv_to_dW(v)
            if all(w == 0 for w in dW): continue
            W2 = [W_base[k]^dW[k] for k in range(16)]
            H2 = sha256_full(W2)
            dh = sum(hw(H_base[k]^H2[k]) for k in range(8))
            khws.append(dh)

        if khws:
            print(f"    K={K:2d}: mean={np.mean(khws):.1f}, min={min(khws)}, std={np.std(khws):.1f}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("5. СРАВНЕНИЕ: GF2-kernel vs random δW")
    print("=" * 70)

    rand_hws = []
    for trial in range(5000):
        dW = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = [W_base[k]^dW[k] for k in range(16)]
        H2 = sha256_full(W2)
        dh = sum(hw(H_base[k]^H2[k]) for k in range(8))
        rand_hws.append(dh)

    print(f"""
  Random δW (5K):     mean={np.mean(rand_hws):.1f}, min={min(rand_hws)}, std={np.std(rand_hws):.1f}
  GF2-kernel single:  mean={np.mean(single_hws):.1f}, min={min(single_hws)}, std={np.std(single_hws):.1f}
  GF2-kernel pairs:   mean={np.mean(pair_hws):.1f}, min={min(pair_hws)}, std={np.std(pair_hws):.1f}
  GF2-kernel triples: mean={np.mean(triple_hws):.1f}, min={min(triple_hws)}, std={np.std(triple_hws):.1f}

  GF2-kernel advantage (pairs vs random):
    Mean: {np.mean(rand_hws) - np.mean(pair_hws):+.1f} bits
    Min:  {min(rand_hws) - min(pair_hws):+d}
    → {'GF2-kernel pairs BETTER!' if np.mean(pair_hws) < np.mean(rand_hws) - 2 else 'No significant advantage'}
""")


if __name__ == "__main__":
    main()
