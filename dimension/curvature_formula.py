"""
УТОЧНЕНИЕ ФОРМУЛЫ K: кривизна зависит от nodes per round?

Данные:
  SHA-256 (2 nodes, 32-bit): K=128 at r=64
  TinyHash (1 node, 16-bit): K=15 at r=8

Гипотезы:
  H1: K = output_bits / 2 → SHA: 128✓, Tiny: 32✗
  H2: K = nodes × C_bits / 2 → SHA: 2×16=32✗
  H3: K зависит от rounds после saturation
  H4: K = f(nodes, C_bits, excess_rounds)
"""

import numpy as np
import struct, hashlib

MASK16 = 0xFFFF
MASK32 = 0xFFFFFFFF

def hw(x): return bin(x).count('1')
def rotr16(x, n): return ((x >> n) | (x << (16 - n))) & MASK16
def rotr32(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def add16(x, y): return (x + y) & MASK16
def add32(x, y): return (x + y) & MASK32

K_T = [0x428a, 0x7137, 0xb5c0, 0xe9b5, 0x3956, 0x59f1, 0x923f, 0xab1c,
       0xd807, 0x1283, 0x2431, 0x550c, 0x72be, 0x80de, 0x9bdc, 0xc19b]
IV_T4 = [0x6a09, 0xbb67, 0x3c6e, 0xa54f]

K_S = [0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
       0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
       0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
       0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174]
IV_S8 = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
         0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]

def tinyhash(W, n_r):
    a,b,c,d = IV_T4
    for r in range(n_r):
        T1 = add16(add16(add16(d, rotr16(a,5)^rotr16(a,11)^(a>>1)), K_T[r%16]), W[r%len(W)])
        T2 = add16(rotr16(a,2)^rotr16(a,7)^(a>>3), ((a&b)^(a&c)^(b&c))&MASK16)
        d,c,b = c,b,a
        a = add16(T1, T2)
    return tuple(add16(IV_T4[i], [a,b,c,d][i]) for i in range(4))

def sha256_n(W, n_r):
    from dimension.carry_skeleton import Sigma0, Sigma1, Ch, Maj
    a,b,c,d,e,f,g,h = IV_S8
    for r in range(n_r):
        T1 = add32(add32(add32(add32(h, rotr32(e,6)^rotr32(e,11)^rotr32(e,25)),
             (e&f)^(~e&g)&MASK32), K_S[r%16]), W[r%len(W)])
        T2 = add32(rotr32(a,2)^rotr32(a,13)^rotr32(a,22),
             (a&b)^(a&c)^(b&c))
        h,g,f,e = g,f,e,add32(d,T1)
        d,c,b,a = c,b,a,add32(T1,T2)
    return tuple(add32(IV_S8[i], [a,b,c,d,e,f,g,h][i]) for i in range(8))


def measure_K(hash_fn, n_input, n_output, C_bits, n_rounds):
    np.random.seed(42)
    mask = (1 << C_bits) - 1
    W_base = [np.random.randint(0, mask+1) for _ in range(n_input)]
    H_base = hash_fn(W_base, n_rounds)
    total_input_bits = n_input * C_bits

    Ks = []
    for _ in range(500):
        b1, b2 = np.random.choice(total_input_bits, 2, replace=False)
        W1 = list(W_base); W1[b1//C_bits] ^= (1 << (b1%C_bits))
        W2 = list(W_base); W2[b2//C_bits] ^= (1 << (b2%C_bits))
        W12 = list(W_base); W12[b1//C_bits] ^= (1 << (b1%C_bits)); W12[b2//C_bits] ^= (1 << (b2%C_bits))
        H1 = hash_fn(W1, n_rounds)
        H2 = hash_fn(W2, n_rounds)
        H12 = hash_fn(W12, n_rounds)
        nl = sum(hw((H_base[i]^H12[i])^((H_base[i]^H1[i])^(H_base[i]^H2[i]))) for i in range(n_output))
        Ks.append(nl)
    return np.mean(Ks)


def main():
    np.random.seed(42)

    print("=" * 70)
    print("УТОЧНЕНИЕ ФОРМУЛЫ K (кривизна)")
    print("=" * 70)

    # K по раундам для TinyHash
    print(f"\n  TinyHash K по раундам (4-reg, 16-bit, 1 node):")
    for n_r in [4, 6, 8, 10, 12, 16, 24, 32]:
        K = measure_K(tinyhash, 8, 4, 16, n_r)
        print(f"    r={n_r:2d}: K={K:.1f}")

    # K по раундам для SHA-256
    print(f"\n  SHA-256 K по раундам (8-reg, 32-bit, 2 nodes):")
    for n_r in [8, 12, 16, 20, 24, 32, 48, 64]:
        K = measure_K(sha256_n, 16, 8, 32, n_r)
        print(f"    r={n_r:2d}: K={K:.1f}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("ФОРМУЛА K: fit")
    print("=" * 70)

    # Данные для fit:
    # SHA-256 r=64: K=128, output=256, nodes=2
    # SHA-256 r=20: K=128 (sphere)
    # TinyHash r=32: K=? (saturated)
    # TinyHash r=8: K=15

    # K_sat = curvature at saturation
    K_tiny_sat = measure_K(tinyhash, 8, 4, 16, 32)
    K_sha_sat = measure_K(sha256_n, 16, 8, 32, 64)

    print(f"\n  K at saturation:")
    print(f"    TinyHash (4-reg, 16-bit): K_sat = {K_tiny_sat:.1f}")
    print(f"    SHA-256  (8-reg, 32-bit): K_sat = {K_sha_sat:.1f}")
    print(f"    Ratio: {K_sha_sat / max(K_tiny_sat, 0.1):.1f}×")

    # Theory: K_sat = output_bits / 2
    # SHA-256: 256/2 = 128 ✓
    # TinyHash: 64/2 = 32, but measured ≈ K_tiny_sat

    print(f"\n  K_sat vs output_bits/2:")
    print(f"    SHA-256: {K_sha_sat:.0f} vs {256/2:.0f} (ratio {K_sha_sat / 128:.2f})")
    print(f"    TinyHash: {K_tiny_sat:.0f} vs {64/2:.0f} (ratio {K_tiny_sat / 32:.2f})")

    # K depends on n_nodes?
    # SHA-256: 2 nodes → K=128 = output/2
    # TinyHash: 1 node → K=K_tiny_sat ≈ output/4?
    # K = nodes × output_bits / (2 × total_transitions)
    # = nodes × N_reg × C_bits / (2 × (nodes + pipes))

    nodes_sha = 2
    nodes_tiny = 1
    pipes_sha = 6
    pipes_tiny = 3

    K_pred_sha = nodes_sha * 256 / (2 * (nodes_sha + pipes_sha)) * (nodes_sha + pipes_sha)
    # Simpler: K = output_bits × (nodes / total_transitions) × correction

    print(f"""
  CORRECTED K FORMULA:

  K_sat depends on SATURATION DEPTH:
    SHA-256: saturates at r≈20 (well beyond boundary r=16)
      → 4 extra rounds after boundary → K fully develops
    TinyHash: boundary=8, testing at r=8
      → 0 extra rounds → K NOT fully developed!

  Let's check TinyHash at r=32 (well beyond boundary):
    K_sat(TinyHash, r=32) = {K_tiny_sat:.1f}

  If K_sat ≈ output/2 at full saturation:
    TinyHash needs r >> 8 to reach K=32.
    At r=8: K=15 (half of saturation).
    At r=32: K={K_tiny_sat:.1f} ({'fully saturated' if abs(K_tiny_sat - 32) < 5 else 'still growing'}).

  FORMULA: K = output_bits/2 × min(1, (r - boundary) / saturation_depth)
  Where saturation_depth ≈ N_reg (4 for TinyHash, 8 for SHA-256).
""")


if __name__ == "__main__":
    main()
