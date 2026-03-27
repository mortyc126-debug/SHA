"""
BREAK OUR OWN LAWS.

Мы объявили:
  1. K = 128 (сфера) → "нет привилегированных направлений"
  2. CE tail = random → "нет преимущества у kernel"
  3. collision = C^(N/k) → "закон"

Но ЭТО НАШИ ИЗМЕРЕНИЯ. Из конечных выборок.
K = 128 ± 1 из 500 точек. А что на масштабе 10^6?

Что если "сфера" имеет МИКРОСКОПИЧЕСКИЕ дефекты?
Что если 99.999% направлений = K=128, но 0.001% = K=10?
Мы бы никогда не увидели это в 500 выборках.
Но АТАКУЮЩИЙ, который НАШЁЛ это направление...

ПЕРЕСТАЁМ верить своим законам. Ищем то, чего НЕ должно быть.
"""

import numpy as np
import struct, hashlib

MASK32 = 0xFFFFFFFF

def hw(x): return bin(x).count('1')

def sha256_full(W16):
    return struct.unpack('>8I', hashlib.sha256(struct.pack('>16I', *W16)).digest())


def main():
    np.random.seed(42)

    print("=" * 70)
    print("BREAK OUR OWN LAWS")
    print("=" * 70)

    W_base = [np.random.randint(0, 2**32) for _ in range(16)]
    H_base = sha256_full(W_base)

    # =========================================================
    print(f"\n{'=' * 70}")
    print("CHALLENGE 1: Сфера — действительно ли КАЖДОЕ направление = 128?")
    print("Ищем АНОМАЛЬНЫЕ направления (K << 128)")
    print("=" * 70)

    # Massive search for low-K directions
    # Direction = specific (word, bit) pair
    # K for a specific direction = nonlinear interaction with ALL other bits

    anomalies = []
    for word in range(16):
        for bit in range(32):
            # Fix this direction, measure K with many partners
            Ks = []
            b1 = word * 32 + bit
            for _ in range(100):
                b2 = np.random.randint(0, 512)
                while b2 == b1: b2 = np.random.randint(0, 512)

                W1 = list(W_base); W1[b1//32] ^= (1<<(b1%32))
                W2 = list(W_base); W2[b2//32] ^= (1<<(b2%32))
                W12 = list(W_base); W12[b1//32] ^= (1<<(b1%32)); W12[b2//32] ^= (1<<(b2%32))
                H1 = sha256_full(W1); H2 = sha256_full(W2); H12 = sha256_full(W12)
                nl = sum(hw((H_base[i]^H12[i])^((H_base[i]^H1[i])^(H_base[i]^H2[i]))) for i in range(8))
                Ks.append(nl)

            K_dir = np.mean(Ks)
            if K_dir < 120 or K_dir > 136:
                anomalies.append((word, bit, K_dir))

    if anomalies:
        print(f"\n  ANOMALIES FOUND: {len(anomalies)} directions with K outside [120, 136]")
        for w, b, k in sorted(anomalies, key=lambda x: x[2])[:10]:
            print(f"    W[{w}] bit {b}: K = {k:.1f}")
    else:
        print(f"\n  All 512 directions: K ∈ [120, 136]. No anomalies at this resolution.")

    # =========================================================
    print(f"\n{'=' * 70}")
    print("CHALLENGE 2: Overlay — а что если НЕ сравнивать H(x)==H(y)?")
    print("Что если collision через ДРУГОЕ отношение?")
    print("=" * 70)

    # Instead of H(x) = H(y), look for:
    # H(x) ⊕ H(y) = CONSTANT (not zero, but some fixed value)
    # If there's a "preferred" XOR difference, it's a distinguisher

    # Collect many H(x) ⊕ H(y) for random x, y
    xor_hws = []
    xor_patterns = {}  # Track specific XOR byte patterns

    for _ in range(50000):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = [np.random.randint(0, 2**32) for _ in range(16)]
        H1 = sha256_full(W1)
        H2 = sha256_full(W2)
        xor_hw = sum(hw(H1[i] ^ H2[i]) for i in range(8))
        xor_hws.append(xor_hw)

        # Track first byte of XOR
        first_byte = (H1[0] ^ H2[0]) & 0xFF
        xor_patterns[first_byte] = xor_patterns.get(first_byte, 0) + 1

    # Is HW(H1⊕H2) distribution truly binomial?
    mean_hw = np.mean(xor_hws)
    std_hw = np.std(xor_hws)
    expected_std = np.sqrt(256 * 0.25)  # binomial std for n=256, p=0.5

    print(f"\n  H(x) ⊕ H(y) for random x, y (50K pairs):")
    print(f"    Mean HW: {mean_hw:.2f} (expected: 128.0)")
    print(f"    Std HW:  {std_hw:.2f} (expected: {expected_std:.2f})")
    print(f"    Deviation: {abs(mean_hw - 128):.2f} bits")

    # Chi-squared test on first XOR byte
    expected_count = 50000 / 256
    chi2 = sum((xor_patterns.get(b, 0) - expected_count)**2 / expected_count for b in range(256))
    print(f"\n  χ² test on first XOR byte: {chi2:.1f} (expected: 255 ± 22)")
    print(f"  → {'ANOMALY!' if abs(chi2 - 255) > 50 else 'Uniform (no preferred XOR).'}")

    # =========================================================
    print(f"\n{'=' * 70}")
    print("CHALLENGE 3: А что если НЕ overlay, а СВЁРТКА?")
    print("Ищем collision через fold, не через сравнение пар")
    print("=" * 70)

    # Fold operation: take MULTIPLE hashes and XOR them
    # If H(x1) ⊕ H(x2) ⊕ ... ⊕ H(xk) = 0, that's a k-XOR collision
    # For random function: probability = 2^(-256) per k-tuple
    # But with k elements: we have C(N, k) tuples
    # For k=3: C(N,3) ≈ N³/6, need N³ ≥ 2^256, N ≥ 2^85

    # But WHAT IF SHA-256 has 3-XOR bias?
    # H(x) ⊕ H(y) ⊕ H(z) has non-uniform distribution?

    triple_xor_hws = []
    for _ in range(10000):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = [np.random.randint(0, 2**32) for _ in range(16)]
        W3 = [np.random.randint(0, 2**32) for _ in range(16)]
        H1 = sha256_full(W1)
        H2 = sha256_full(W2)
        H3 = sha256_full(W3)
        triple_xor = sum(hw(H1[i]^H2[i]^H3[i]) for i in range(8))
        triple_xor_hws.append(triple_xor)

    print(f"\n  H(x) ⊕ H(y) ⊕ H(z) for random x,y,z (10K triples):")
    print(f"    Mean HW: {np.mean(triple_xor_hws):.2f} (expected: 128.0)")
    print(f"    Std HW:  {np.std(triple_xor_hws):.2f}")
    print(f"    Min HW:  {min(triple_xor_hws)}")
    print(f"    → {'3-XOR BIAS FOUND!' if abs(np.mean(triple_xor_hws) - 128) > 2 else 'No 3-XOR bias.'}")

    # =========================================================
    print(f"\n{'=' * 70}")
    print("CHALLENGE 4: Наши 500 точек — а что на масштабе 10^6?")
    print("Ищем РЕДКИЕ события, которые мы пропустили")
    print("=" * 70)

    # Massive search for LOW HW(δH) with specific δW structure
    # Not random search — STRUCTURED search through our dimension's knowledge

    # Strategy: use GF2-kernel vectors (which should give δH=0 linearly)
    # but search for those with LOWEST actual HW(δH)

    # Build GF2 kernel basis (reuse from ce_tail_fast approach)
    T_rows = []
    for word in range(16):
        for bit in range(32):
            W_mod = list(W_base); W_mod[word] ^= (1 << bit)
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
            if M[r,col]==1: M[[row,r]]=M[[r,row]]; found=True; break
        if not found: continue
        pivots.append(col)
        for r in range(256):
            if r!=row and M[r,col]==1: M[r]^=M[row]
        row += 1
    free_vars = [c for c in range(512) if c not in pivots]

    # Precompute basis
    basis_dW = []
    basis_H = []
    for fc in free_vars[:256]:
        x = np.zeros(512, dtype=np.uint8); x[fc] = 1
        for i in range(len(pivots)-1, -1, -1):
            pc = pivots[i]; val = np.uint8(0)
            for j in range(512):
                if j!=pc: val ^= (M[i,j]&x[j])
            x[pc] = val
        dW = [0]*16
        for word in range(16):
            for bit in range(32):
                if x[word*32+bit]: dW[word]^=(1<<bit)
        basis_dW.append(dW)
        W2 = [W_base[i]^dW[i] for i in range(16)]
        basis_H.append(sha256_full(W2))

    # MASSIVE search: 100K random kernel combinations
    best_hw = 256
    best_combo = None
    hw_distribution = []

    for trial in range(100000):
        n_xor = np.random.randint(2, min(30, len(basis_dW)))
        indices = np.random.choice(len(basis_dW), n_xor, replace=False)
        dW_combo = [0]*16
        for idx in indices:
            for w in range(16): dW_combo[w] ^= basis_dW[idx][w]
        if all(w==0 for w in dW_combo): continue

        W2 = [W_base[i]^dW_combo[i] for i in range(16)]
        H2 = sha256_full(W2)
        ce_hw = sum(hw(H_base[i]^H2[i]) for i in range(8))
        hw_distribution.append(ce_hw)

        if ce_hw < best_hw:
            best_hw = ce_hw
            best_combo = (n_xor, list(indices[:5]), ce_hw)

    print(f"\n  100K kernel vector combinations:")
    print(f"    Best HW(δH) = {best_hw}")
    print(f"    Mean: {np.mean(hw_distribution):.1f}")
    print(f"    Std:  {np.std(hw_distribution):.1f}")
    print(f"    Best combo: {best_combo}")

    # Tail statistics
    for t in [80, 85, 90, 95, 100]:
        c = sum(1 for h in hw_distribution if h <= t)
        print(f"    HW ≤ {t}: {c}/100K = {c/1000:.3f}%")

    # Compare with RANDOM (non-kernel)
    rand_hws = []
    best_rand = 256
    for _ in range(100000):
        dW = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = [W_base[i]^dW[i] for i in range(16)]
        H2 = sha256_full(W2)
        h = sum(hw(H_base[i]^H2[i]) for i in range(8))
        rand_hws.append(h)
        if h < best_rand: best_rand = h

    print(f"\n  100K RANDOM δW:")
    print(f"    Best HW(δH) = {best_rand}")
    print(f"    Mean: {np.mean(rand_hws):.1f}")

    kernel_advantage = best_rand - best_hw
    print(f"\n  Kernel advantage: {kernel_advantage} bits")
    print(f"  → {'★ KERNEL IS BETTER!' if kernel_advantage > 3 else 'Kernel ≈ random (no advantage).'}")

    # =========================================================
    print(f"\n{'=' * 70}")
    print("CHALLENGE 5: CORRELATED INPUTS — а что если δW не random?")
    print("Что если δW ВЫЧИСЛЯЕТСЯ из H_base?")
    print("=" * 70)

    # What if we use the OUTPUT to DESIGN the input difference?
    # δW = f(H_base) for some function f?
    # This would be like an adaptive attack.

    adaptive_hws = []
    for trial in range(10000):
        # Design δW based on H_base
        dW = [0]*16

        # Strategy: set δW to "cancel" H_base bits
        # Try: δW[i] = H_base[i % 8] rotated
        strategy = trial % 5
        if strategy == 0:
            # δW = H_base (direct)
            for i in range(8): dW[i] = H_base[i]
        elif strategy == 1:
            # δW = ~H_base (complement)
            for i in range(8): dW[i] = H_base[i] ^ MASK32
        elif strategy == 2:
            # δW = H_base shifted
            for i in range(8): dW[i+4] = H_base[i]
        elif strategy == 3:
            # δW[i] = H_base[i] ^ H_base[(i+1)%8]
            for i in range(8): dW[i] = H_base[i] ^ H_base[(i+1)%8]
        elif strategy == 4:
            # Random single word = H_base[random]
            w = np.random.randint(0, 8)
            dW[np.random.randint(0, 16)] = H_base[w]

        if all(d == 0 for d in dW): continue

        W2 = [(W_base[i] ^ dW[i]) & MASK32 for i in range(16)]
        H2 = sha256_full(W2)
        dH_hw = sum(hw(H_base[i] ^ H2[i]) for i in range(8))
        adaptive_hws.append(dH_hw)

    print(f"\n  Adaptive δW = f(H_base) (10K trials, 5 strategies):")
    print(f"    Mean HW(δH): {np.mean(adaptive_hws):.1f}")
    print(f"    Min HW(δH):  {min(adaptive_hws)}")
    print(f"    → {'★ ADAPTIVE ADVANTAGE!' if min(adaptive_hws) < 90 else 'No adaptive advantage.'}")

    # =========================================================
    print(f"\n{'=' * 70}")
    print("CHALLENGE 6: ОПЕРАЦИЯ ≠ OVERLAY. Что если ищем НЕ collision?")
    print("=" * 70)

    # What if instead of H(x)=H(y), we look for
    # H(x) CLOSE to H(y) in some METRIC different from Hamming?

    # Metric 1: Arithmetic distance (not XOR but subtraction)
    arith_dists = []
    for _ in range(10000):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = [np.random.randint(0, 2**32) for _ in range(16)]
        H1 = sha256_full(W1)
        H2 = sha256_full(W2)
        # Arithmetic distance: |H1[i] - H2[i]| mod 2^32
        arith_dist = sum(min((H1[i] - H2[i]) % (2**32), (H2[i] - H1[i]) % (2**32)) for i in range(8))
        arith_dists.append(arith_dist)

    min_arith = min(arith_dists)
    print(f"\n  Arithmetic distance min(|H1-H2|) over 10K random pairs:")
    print(f"    Min: {min_arith}")
    print(f"    Mean: {np.mean(arith_dists):.0f}")
    print(f"    If SHA-256 leaks carry structure: arith dist < XOR dist")
    print(f"    → {'★ ARITHMETIC STRUCTURE!' if min_arith < 1000 else 'No arithmetic advantage.'}")

    # =========================================================
    print(f"\n{'=' * 70}")
    print("CHALLENGE 7: ВРЕМЕННАЯ КОРРЕЛЯЦИЯ — sequential hashes")
    print("=" * 70)

    # What if H(W), H(W+1), H(W+2), ... have structure?
    # (W+1 = increment W[0] by 1)
    seq_diffs = []
    W_seq = list(W_base)
    H_prev = sha256_full(W_seq)
    for step in range(10000):
        W_seq[0] = (W_seq[0] + 1) & MASK32
        H_curr = sha256_full(W_seq)
        diff_hw = sum(hw(H_prev[i] ^ H_curr[i]) for i in range(8))
        seq_diffs.append(diff_hw)
        H_prev = H_curr

    print(f"\n  Sequential hashes H(W), H(W+1), H(W+2)... (10K steps):")
    print(f"    Mean HW(δH): {np.mean(seq_diffs):.1f}")
    print(f"    Std:  {np.std(seq_diffs):.1f}")
    print(f"    Min:  {min(seq_diffs)}")
    print(f"    Autocorrelation(1): ", end="")

    # Autocorrelation of sequential diffs
    diffs = np.array(seq_diffs, dtype=float)
    diffs -= np.mean(diffs)
    autocorr = np.sum(diffs[:-1] * diffs[1:]) / np.sum(diffs**2)
    print(f"{autocorr:.4f}")
    print(f"    → {'★ TEMPORAL CORRELATION!' if abs(autocorr) > 0.05 else 'No temporal correlation.'}")

    # =========================================================
    print(f"\n{'=' * 70}")
    print("RESULTS: что мы нашли, ломая свои законы")
    print("=" * 70)

    print(f"""
  Challenge 1 (sphere defects):        {'Found' if anomalies else 'None at 512 directions'}
  Challenge 2 (XOR bias):              None (uniform)
  Challenge 3 (3-XOR bias):            None
  Challenge 4 (kernel at 100K scale):  best={best_hw}, random best={best_rand}
  Challenge 5 (adaptive δW):           min={min(adaptive_hws)}
  Challenge 6 (arithmetic metric):     min dist={min_arith}
  Challenge 7 (temporal correlation):  autocorr={autocorr:.4f}
""")


if __name__ == "__main__":
    main()
