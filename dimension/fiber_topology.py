"""
НОВОЕ НАПРАВЛЕНИЕ 3: Топология волокон.

Волокно (fiber) = множество всех прообразов одного хеш-значения.
H(x₁) = H(x₂) = ... = H(xₖ) = h → {x₁, x₂, ..., xₖ} = fiber(h).

Вопросы:
  A. Прообразы кластеризуются в пространстве входов?
     (близкие x → один h? или разбросаны случайно?)
  B. Расстояние между прообразами: uniform или structured?
  C. Есть ли "магнитные линии" — направления в которых
     прообразы выстраиваются?

На полном SHA-256: невозможно (2^256 прообразов).
На mini-SHA (8-bit output): можно найти ВСЕ прообразы.
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
]
IV = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]


def sha256_truncated(W16, n_rounds, out_bits):
    """SHA-256 truncated to out_bits."""
    W = list(W16)
    for r in range(16, max(n_rounds, 16)):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))
    a,b,c,d,e,f,g,h = IV
    for r in range(n_rounds):
        T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K[r%16]), W[r])
        T2 = add32(Sigma0(a), Maj(a,b,c))
        h,g,f,e = g,f,e,add32(d,T1)
        d,c,b,a = c,b,a,add32(T1,T2)
    # Truncate: take lowest out_bits from H[0]
    H0 = add32(IV[0], a)
    return H0 & ((1 << out_bits) - 1)


def main():
    np.random.seed(42)

    print("=" * 70)
    print("ТОПОЛОГИЯ ВОЛОКОН: структура прообразов")
    print("=" * 70)

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("A. Mini-SHA: 8-bit output, vary W[0] only (2^16 inputs)")
    print("=" * 70)

    # Fix W[1..15], vary W[0] over 2^16 values
    # Truncate output to 8 bits → 256 possible outputs
    # Each output has ~2^16/256 = 256 preimages

    W_base = [np.random.randint(0, 2**32) for _ in range(16)]
    n_rounds = 8
    out_bits = 8

    # Build fibers
    fibers = {}  # hash_value → list of W[0] values
    N_inputs = 65536  # 2^16

    for i in range(N_inputs):
        W = list(W_base)
        W[0] = i  # vary only W[0] in [0, 65535]
        h = sha256_truncated(W, n_rounds, out_bits)
        if h not in fibers:
            fibers[h] = []
        fibers[h].append(i)

    # Fiber sizes
    sizes = [len(f) for f in fibers.values()]
    print(f"  {N_inputs} inputs → {len(fibers)} distinct outputs")
    print(f"  Fiber sizes: mean={np.mean(sizes):.1f}, std={np.std(sizes):.1f}")
    print(f"    min={min(sizes)}, max={max(sizes)}")

    # Expected for random function: Poisson(λ=256)
    # std = sqrt(256) = 16
    expected_std = np.sqrt(np.mean(sizes))
    print(f"  Expected std (Poisson): {expected_std:.1f}")
    print(f"  Actual std: {np.std(sizes):.1f}")
    print(f"  → {'CLUSTERED (std too high)' if np.std(sizes) > expected_std * 1.5 else 'UNIFORM (matches Poisson)' if np.std(sizes) < expected_std * 1.3 else 'BORDERLINE'}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("B. INTRA-FIBER DISTANCE: how far apart are preimages?")
    print("=" * 70)

    # Pick a fiber with many preimages
    # Measure pairwise Hamming distance between preimages
    biggest_fiber_hash = max(fibers, key=lambda h: len(fibers[h]))
    bf = fibers[biggest_fiber_hash]

    print(f"  Biggest fiber: h={biggest_fiber_hash}, size={len(bf)}")

    # Pairwise distances within fiber
    intra_dists = []
    sample = bf[:200]  # limit for speed
    for i in range(len(sample)):
        for j in range(i+1, len(sample)):
            intra_dists.append(hw(sample[i] ^ sample[j]))

    # Compare with random pairs from ALL inputs
    inter_dists = []
    all_inputs = list(range(N_inputs))
    for _ in range(len(intra_dists)):
        i, j = np.random.choice(N_inputs, 2, replace=False)
        inter_dists.append(hw(i ^ j))

    print(f"\n  Intra-fiber distances (within same hash):")
    print(f"    Mean: {np.mean(intra_dists):.2f}")
    print(f"    Std:  {np.std(intra_dists):.2f}")
    print(f"  Inter-fiber distances (random pairs):")
    print(f"    Mean: {np.mean(inter_dists):.2f}")
    print(f"    Std:  {np.std(inter_dists):.2f}")

    from scipy import stats
    t, p = stats.ttest_ind(intra_dists, inter_dists)
    print(f"  t={t:.2f}, p={p:.2e}")
    print(f"  → {'CLUSTERED! (preimages are closer)' if p < 0.01 and np.mean(intra_dists) < np.mean(inter_dists) else 'UNIFORM (no clustering)'}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("C. FIBER STRUCTURE: are preimages evenly spaced?")
    print("=" * 70)

    # For a fiber: sort preimages, compute gaps
    # If uniform: gaps ~ geometric distribution
    # If clustered: gaps are bimodal (small within cluster, large between)

    mid_fiber_hash = sorted(fibers.keys(), key=lambda h: abs(len(fibers[h]) - 256))[0]
    mf = sorted(fibers[mid_fiber_hash])
    gaps = [mf[i+1] - mf[i] for i in range(len(mf)-1)]

    print(f"  Fiber h={mid_fiber_hash}, size={len(mf)}")
    print(f"  Gaps between sorted preimages:")
    print(f"    Mean: {np.mean(gaps):.1f}")
    print(f"    Std:  {np.std(gaps):.1f}")
    print(f"    Expected (uniform): mean={N_inputs/len(mf):.1f}, std≈{N_inputs/len(mf):.1f}")
    print(f"    Min gap: {min(gaps)}, Max gap: {max(gaps)}")

    # Coefficient of variation
    cv = np.std(gaps) / np.mean(gaps)
    print(f"    CV (std/mean): {cv:.3f}")
    print(f"    Expected CV (exponential): 1.000")
    print(f"    → {'STRUCTURED' if abs(cv - 1.0) > 0.2 else 'RANDOM (matches exponential)'}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("D. SCALE UP: 16-bit output, 2^20 inputs")
    print("=" * 70)

    out_bits = 16
    N_inputs = 1 << 20  # 1M

    fibers16 = {}
    for i in range(N_inputs):
        W = list(W_base)
        W[0] = i & MASK32
        h = sha256_truncated(W, n_rounds, out_bits)
        if h not in fibers16:
            fibers16[h] = []
        fibers16[h].append(i)

    sizes16 = [len(f) for f in fibers16.values()]
    expected_size = N_inputs / (1 << out_bits)
    print(f"  {N_inputs} inputs → {len(fibers16)} distinct outputs (of {1<<out_bits} possible)")
    print(f"  Fiber sizes: mean={np.mean(sizes16):.1f} (expected: {expected_size:.1f})")
    print(f"  Std: {np.std(sizes16):.1f} (expected Poisson: {np.sqrt(expected_size):.1f})")

    # Sample fibers and check distances
    big_fibers = [h for h, f in fibers16.items() if len(f) > 20]
    if big_fibers:
        intra16 = []
        inter16 = []
        for fh in big_fibers[:50]:
            fb = fibers16[fh]
            sample = fb[:30]
            for i in range(len(sample)):
                for j in range(i+1, len(sample)):
                    intra16.append(hw(sample[i] ^ sample[j]))

        for _ in range(len(intra16)):
            i, j = np.random.choice(N_inputs, 2, replace=False)
            inter16.append(hw(i ^ j))

        print(f"\n  Intra-fiber (16-bit): mean={np.mean(intra16):.2f}")
        print(f"  Inter-fiber (random): mean={np.mean(inter16):.2f}")
        t16, p16 = stats.ttest_ind(intra16, inter16)
        print(f"  t={t16:.2f}, p={p16:.2e}")
        print(f"  → {'CLUSTERED!' if p16 < 0.01 and np.mean(intra16) < np.mean(inter16) else 'UNIFORM'}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("E. ROUND DEPENDENCE: do fibers cluster at low rounds?")
    print("=" * 70)

    # At round 1: almost no mixing → fibers should be VERY clustered
    # At round 8: more mixing → less clustered
    # At round 16: full mixing → no clustering

    for nr in [1, 2, 4, 8]:
        fibers_r = {}
        for i in range(N_inputs):
            W = list(W_base)
            W[0] = i & MASK32
            h = sha256_truncated(W, nr, 8)  # 8-bit output
            if h not in fibers_r:
                fibers_r[h] = []
            fibers_r[h].append(i)

        # Sample intra vs inter distances
        big_f = [h for h, f in fibers_r.items() if len(f) > 20]
        if big_f:
            intra_r = []
            for fh in big_f[:20]:
                fb = fibers_r[fh][:30]
                for i in range(len(fb)):
                    for j in range(i+1, len(fb)):
                        intra_r.append(hw(fb[i] ^ fb[j]))

            inter_r = []
            for _ in range(len(intra_r)):
                a, b = np.random.choice(N_inputs, 2, replace=False)
                inter_r.append(hw(a ^ b))

            ratio = np.mean(intra_r) / np.mean(inter_r) if np.mean(inter_r) > 0 else 0
            print(f"  r={nr}: intra={np.mean(intra_r):.2f}, inter={np.mean(inter_r):.2f}, ratio={ratio:.3f}")
        else:
            print(f"  r={nr}: no fibers with >20 preimages")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("РЕЗУЛЬТАТ")
    print("=" * 70)


if __name__ == "__main__":
    main()
