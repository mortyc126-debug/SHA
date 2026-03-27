"""
ПРОЕКТИРУЕМ НОВУЮ ХЕШ-ФУНКЦИЮ ИЗ НАШИХ ТЕОРЕМ.

Наши законы → конкретные дизайн-решения:

  1. λ = C/4.5 → optimal_rounds = input_words + ceil(C/λ) + 2*N
     Для C=32, N=8: 16 + 5 + 16 = 37. Берём 29 (3N+5).

  2. Dual nonlinearity: carry ALONE sufficient.
     → Убираем Ch, Maj. Оставляем ТОЛЬКО ADD mod 2^32.

  3. Pipe conservation lifetime = N/2 = 4.
     → Достаточно 4 раунда для полного pipe cycle.

  4. Sphere onset ≈ input_words + 5.
     → После r=21 все метрики стабилизируются.

ДИЗАЙН: "DimHash-256"
  - 8 регистров × 32 бит = 256 бит output
  - 29 раундов (3×8+5, наш оптимум)
  - ТОЛЬКО carry nonlinearity (no boolean functions)
  - Rotations для диффузии (σ0, σ1, Σ0, Σ1 — сохранены)
  - Schedule: стандартный SHA-256

Потом сравниваем с SHA-256 по ВСЕМ нашим метрикам.
"""

import numpy as np
import struct, hashlib

MASK32 = 0xFFFFFFFF

def hw(x): return bin(x).count('1')
def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def add32(*args):
    r = 0
    for a in args: r = (r + a) & MASK32
    return r
def sub32(x, y): return (x - y) & MASK32

# SHA-256 functions (for comparison)
def Sigma0_sha(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sigma1_sha(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
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


def dimhash_256(W16, n_rounds=29):
    """
    DimHash-256: designed from our dimension's theorems.
    NO boolean nonlinearity (no Ch, no Maj).
    ONLY carry (ADD mod 2^32) + rotations.
    29 rounds (optimal from 3N+5).
    """
    # Schedule: same as SHA-256
    W = list(W16)
    for r in range(16, max(n_rounds, 16)):
        W.append(add32(sigma1(W[r-2]), W[r-7], sigma0(W[r-15]), W[r-16]))

    a,b,c,d,e,f,g,h = IV

    for r in range(n_rounds):
        # T1 = h + Σ1(e) + (e + f + g) + K[r] + W[r]
        # (e+f+g replaces Ch(e,f,g) — carry-only nonlinearity)
        Sig1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25)
        T1 = add32(h, Sig1, add32(e, f, g), K[r], W[r])

        # T2 = Σ0(a) + (a + b + c)
        # (a+b+c replaces Maj(a,b,c) — carry-only nonlinearity)
        Sig0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22)
        T2 = add32(Sig0, add32(a, b, c))

        # Shift register (same as SHA-256)
        h,g,f,e = g,f,e,add32(d,T1)
        d,c,b,a = c,b,a,add32(T1,T2)

    return tuple(add32(IV[i], [a,b,c,d,e,f,g,h][i]) for i in range(8))


def sha256_custom(W16, n_rounds=64):
    """Standard SHA-256 with configurable rounds."""
    W = list(W16)
    for r in range(16, max(n_rounds, 16)):
        W.append(add32(sigma1(W[r-2]), W[r-7], sigma0(W[r-15]), W[r-16]))
    a,b,c,d,e,f,g,h = IV
    for r in range(n_rounds):
        T1 = add32(h, Sigma1_sha(e), Ch(e,f,g), K[r], W[r])
        T2 = add32(Sigma0_sha(a), Maj(a,b,c))
        h,g,f,e = g,f,e,add32(d,T1)
        d,c,b,a = c,b,a,add32(T1,T2)
    return tuple(add32(IV[i], [a,b,c,d,e,f,g,h][i]) for i in range(8))


def measure_metrics(hash_fn, name, n_samples=2000):
    """Measure all our dimension's metrics for a hash function."""
    metrics = {}

    # 1. Avalanche (K metric): HW(δH) for 1-bit flip
    avalanche_hws = []
    for _ in range(n_samples):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H1 = hash_fn(W)
        w = np.random.randint(0, 16)
        b = np.random.randint(0, 32)
        W2 = list(W); W2[w] ^= (1 << b)
        H2 = hash_fn(W2)
        avalanche_hws.append(sum(hw(H1[i]^H2[i]) for i in range(8)))
    metrics['K'] = np.mean(avalanche_hws)
    metrics['K_std'] = np.std(avalanche_hws)

    # 2. Output HW distribution
    out_hws = []
    for _ in range(n_samples):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H = hash_fn(W)
        out_hws.append(sum(hw(h) for h in H))
    metrics['HW_mean'] = np.mean(out_hws)
    metrics['HW_std'] = np.std(out_hws)

    # 3. Bit independence: correlation between output bits
    bit_corr = []
    bit_data = np.zeros((n_samples, 256))
    for trial in range(n_samples):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H = hash_fn(W)
        for w_idx in range(8):
            for b_idx in range(32):
                bit_data[trial, w_idx*32+b_idx] = (H[w_idx] >> b_idx) & 1

    # Sample 100 random bit pairs
    for _ in range(100):
        i, j = np.random.choice(256, 2, replace=False)
        c = abs(np.corrcoef(bit_data[:, i], bit_data[:, j])[0, 1])
        bit_corr.append(c)
    metrics['bit_corr'] = np.mean(bit_corr)

    # 4. Differential uniformity: is δH uniform for different δW?
    dH_per_word = {}
    for w in range(16):
        dHs = []
        for _ in range(200):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            W2 = list(W); W2[w] ^= 1
            H1 = hash_fn(W); H2 = hash_fn(W2)
            dHs.append(sum(hw(H1[i]^H2[i]) for i in range(8)))
        dH_per_word[w] = np.mean(dHs)
    metrics['diff_uniformity'] = np.std(list(dH_per_word.values()))

    # 5. CE matrix rank (GF2)
    W_base = [np.random.randint(0, 2**32) for _ in range(16)]
    H_base = hash_fn(W_base)
    rows = []
    for in_w in range(8):  # sample 256 of 512 input bits
        for in_b in range(32):
            W_mod = list(W_base); W_mod[in_w] ^= (1 << in_b)
            H_mod = hash_fn(W_mod)
            row = []
            for out_w in range(8):
                for out_b in range(32):
                    row.append(((H_base[out_w]^H_mod[out_w]) >> out_b) & 1)
            rows.append(row)
    T = np.array(rows, dtype=np.uint8)
    # GF2 rank
    M = T.copy()
    rank = 0
    for col in range(256):
        found = False
        for r in range(rank, min(len(rows), 256)):
            if M[r, col] == 1:
                M[[rank, r]] = M[[r, rank]]
                found = True
                break
        if not found:
            continue
        for r in range(256):
            if r != rank and M[r, col] == 1:
                M[r] ^= M[rank]
        rank += 1
    metrics['CE_rank'] = rank

    return metrics


def main():
    np.random.seed(42)

    print("=" * 70)
    print("DimHash-256: hash function designed from our theorems")
    print("=" * 70)

    print(f"""
  DESIGN DECISIONS (from our dimension's laws):

  1. Rounds: 29 (from 3N+5 = 3×8+5)
     SHA-256 uses 64 → we save 55% computation.

  2. Nonlinearity: ONLY carry (ADD mod 2^32)
     Ch(e,f,g) → replaced by (e + f + g)
     Maj(a,b,c) → replaced by (a + b + c)
     From dual nonlinearity theorem: carry alone sufficient.

  3. Rotations: same as SHA-256 (Σ0, Σ1, σ0, σ1)
     From our analysis: rotations contribute 1 bit to rank.

  4. Schedule: same as SHA-256
     Our analysis showed schedule is well-designed.

  5. Constants: same K[0..28] and IV
     From anomaly 6: K determines carry profile, proven benign.
""")

    # ═══════════════════
    print(f"{'=' * 70}")
    print("COMPARISON: DimHash-256 vs SHA-256 (64r) vs SHA-256 (29r)")
    print("=" * 70)

    hashes = [
        (lambda W: dimhash_256(W, 29), "DimHash-256 (29r)"),
        (lambda W: sha256_custom(W, 29), "SHA-256 (29r)"),
        (lambda W: sha256_custom(W, 64), "SHA-256 (64r)"),
    ]

    all_metrics = {}
    for hash_fn, name in hashes:
        print(f"\n  Measuring {name}...")
        m = measure_metrics(hash_fn, name, n_samples=2000)
        all_metrics[name] = m

    # Print comparison table
    print(f"\n  {'Metric':<25}", end="")
    for _, name in hashes:
        print(f" {name:>18}", end="")
    print(f" {'Ideal':>10}")

    metrics_to_show = [
        ('K', 'Avalanche (K)', 128.0),
        ('K_std', 'Avalanche std', 8.0),
        ('HW_mean', 'Output HW mean', 128.0),
        ('HW_std', 'Output HW std', 8.0),
        ('bit_corr', 'Bit correlation', 0.0),
        ('diff_uniformity', 'Diff uniformity', 0.0),
        ('CE_rank', 'CE rank (GF2)', 256),
    ]

    for key, label, ideal in metrics_to_show:
        print(f"  {label:<25}", end="")
        for _, name in hashes:
            val = all_metrics[name][key]
            # Color-code: close to ideal = good
            diff = abs(val - ideal) / max(ideal, 1)
            marker = "✓" if diff < 0.05 else "~" if diff < 0.15 else "✗"
            print(f" {val:>15.2f} {marker}", end="")
        print(f" {ideal:>9}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("NEAR-COLLISION RESISTANCE: best δH from 5000 random δW")
    print("=" * 70)

    for hash_fn, name in hashes:
        best = 256
        for _ in range(5000):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            W2 = list(W)
            w = np.random.randint(0, 16)
            b = np.random.randint(0, 32)
            W2[w] ^= (1 << b)
            H1 = hash_fn(W); H2 = hash_fn(W2)
            dH = sum(hw(H1[i]^H2[i]) for i in range(8))
            if dH < best: best = dH
        print(f"  {name}: min HW(δH) = {best}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("CONSERVATION LAW: does DimHash preserve it?")
    print("=" * 70)

    W_test = [np.random.randint(0, 2**32) for _ in range(16)]

    # Compute DimHash states
    W = list(W_test)
    for r in range(16, 29):
        W.append(add32(sigma1(W[r-2]), W[r-7], sigma0(W[r-15]), W[r-16]))
    a,b,c,d,e,f,g,h = IV
    dim_states = [(a,b,c,d,e,f,g,h)]
    for r in range(29):
        Sig1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25)
        T1 = add32(h, Sig1, add32(e, f, g), K[r], W[r])
        Sig0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22)
        T2 = add32(Sig0, add32(a, b, c))
        h,g,f,e = g,f,e,add32(d,T1)
        d,c,b,a = c,b,a,add32(T1,T2)
        dim_states.append((a,b,c,d,e,f,g,h))

    # Check conservation
    violations = 0
    for r in range(28):
        ae_r = add32(dim_states[r][0], dim_states[r][4])
        bf_r1 = add32(dim_states[r+1][1], dim_states[r+1][5])
        if ae_r != bf_r1:
            violations += 1

    print(f"  DimHash-256: conservation violations = {violations}/28")
    print(f"  → {'CONSERVED ✓' if violations == 0 else 'BROKEN ✗'}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("COST ANALYSIS")
    print("=" * 70)

    # Operations per round:
    # SHA-256: 2 rotations (Σ0,Σ1) + 2 boolean (Ch,Maj) + 7 ADD + schedule
    # DimHash: 2 rotations (Σ0,Σ1) + 0 boolean + 9 ADD + schedule (same)
    # Boolean ops: Ch = 3 XOR/AND, Maj = 3 XOR/AND = 6 ops
    # ADD: ~32 full adder gates each

    sha_ops = 64 * (2 + 6 + 7)  # rounds × (rotations + boolean + adds)
    dim_ops = 29 * (2 + 0 + 9)  # rounds × (rotations + 0 + adds)

    print(f"  SHA-256 (64r): {sha_ops} operations")
    print(f"  DimHash (29r): {dim_ops} operations")
    print(f"  Speedup: {sha_ops/dim_ops:.1f}×")
    print(f"  Gate savings: {(1 - dim_ops/sha_ops)*100:.0f}%")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("VERDICT")
    print("=" * 70)

    # Compare all metrics
    dim_m = all_metrics["DimHash-256 (29r)"]
    sha64_m = all_metrics["SHA-256 (64r)"]

    checks = [
        ("Avalanche K", abs(dim_m['K'] - 128) < 2, abs(sha64_m['K'] - 128) < 2),
        ("CE rank", dim_m['CE_rank'] == 256, sha64_m['CE_rank'] == 256),
        ("Bit independence", dim_m['bit_corr'] < 0.05, sha64_m['bit_corr'] < 0.05),
        ("HW distribution", abs(dim_m['HW_mean'] - 128) < 2, abs(sha64_m['HW_mean'] - 128) < 2),
    ]

    all_pass = True
    print(f"\n  {'Check':<25} {'DimHash':>10} {'SHA-256':>10}")
    for name, dim_ok, sha_ok in checks:
        print(f"  {name:<25} {'PASS' if dim_ok else 'FAIL':>10} {'PASS' if sha_ok else 'FAIL':>10}")
        if not dim_ok:
            all_pass = False


if __name__ == "__main__":
    main()
