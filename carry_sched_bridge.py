#!/usr/bin/env python3
"""
Carry-Sched Bridge — Stage 1.3: Does schedule nonlinearity penetrate chaos?
============================================================================

carry alone → H[7]: corr = 0.005 (ZERO, proven in G3)
carry_sched = W_real[r] - W_xor[r] — the nonlinear correction.

KEY HYPOTHESIS: carry_sched travels through W[16..63] (OUTSIDE state),
bypassing the chaotic zone entirely. If carry_sched correlates with H[7],
we found the bridge.

Four experiments:
  S1: carry_sched profile — what does it look like?
  S2: carry_sched vs H[7] — the bridge test
  S3: carry_sched decomposition — which rounds carry signal?
  S4: Combined (carry + carry_sched) vs H[7] — synergy test
"""

import numpy as np
from collections import defaultdict, Counter
import time
import sys

MASK = 0xFFFFFFFF

H0 = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
]

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

def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK
def sig0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sig1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)
def Sig0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sig1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)


def full_trace(W16):
    """
    Returns everything: carry profile, carry_sched profile,
    H[0..7], W_real[16..63], W_xor[16..63], SS63.
    """
    # Real schedule
    W_real = list(W16) + [0] * 48
    for i in range(16, 64):
        W_real[i] = (sig1(W_real[i-2]) + W_real[i-7] + sig0(W_real[i-15]) + W_real[i-16]) & MASK

    # XOR schedule
    W_xor = list(W16) + [0] * 48
    for i in range(16, 64):
        W_xor[i] = sig1(W_xor[i-2]) ^ W_xor[i-7] ^ sig0(W_xor[i-15]) ^ W_xor[i-16]

    # carry_sched[r] = W_real[r] - W_xor[r] mod 2^32 for r=16..63
    cs_profile = []
    for r in range(16, 64):
        cs_profile.append((W_real[r] - W_xor[r]) & MASK)

    # SHA-256 compression with carry tracing
    a, b, c, d, e, f, g, h = H0
    carries = []

    for r in range(64):
        raw = h + Sig1(e) + Ch(e, f, g) + K[r] + W_real[r]
        T1 = raw & MASK
        carries.append(1 if raw >= (1 << 32) else 0)

        T2 = (Sig0(a) + Maj(a, b, c)) & MASK
        h, g, f = g, f, e
        e = (d + T1) & MASK
        d, c, b = c, b, a
        a = (T1 + T2) & MASK

    H_out = [(v + iv) & MASK for v, iv in zip([a,b,c,d,e,f,g,h], H0)]

    return {
        'carries': tuple(carries),
        'cs_profile': cs_profile,  # 48 values, carry_sched[16..63]
        'H': H_out,
        'H7': H_out[7],
        'H6': H_out[6],
        'W_real': W_real,
        'W_xor': W_xor,
    }


def hw(x):
    """Hamming weight of 32-bit int."""
    return bin(x & MASK).count('1')


# ============================================================
# S1: carry_sched profile structure
# ============================================================

def experiment_S1(N=15000, seed=60):
    print("=" * 70)
    print("S1: carry_sched Profile Structure")
    print(f"N={N}, seed={seed}")
    print("=" * 70)

    rng = np.random.RandomState(seed)

    # Collect carry_sched statistics per round
    cs_means = np.zeros(48)
    cs_hws = np.zeros(48)
    cs_zeros = np.zeros(48)  # count of cs[r]=0

    # Also collect a compact fingerprint: HW of each cs word
    fingerprints = []

    t0 = time.time()
    for i in range(N):
        W16 = [int(rng.randint(0, 1 << 32)) for _ in range(16)]
        t = full_trace(W16)
        cs = t['cs_profile']

        fp = tuple(hw(cs[j]) for j in range(48))
        fingerprints.append(fp)

        for j in range(48):
            cs_means[j] += cs[j]
            cs_hws[j] += hw(cs[j])
            if cs[j] == 0:
                cs_zeros[j] += 1

        if (i + 1) % 5000 == 0:
            print(f"  {i+1}/{N} ({time.time()-t0:.1f}s)")

    cs_means /= N
    cs_hws /= N
    cs_zeros_frac = cs_zeros / N

    print(f"\nCompleted: {time.time()-t0:.1f}s")

    # Per-round analysis
    print(f"\n{'r':>3} | {'E[cs]/2^32':>11} | {'E[HW(cs)]':>10} | {'P(cs=0)':>8} | {'bias':>6}")
    print("-" * 55)
    for j in range(48):
        r = j + 16
        ratio = cs_means[j] / (1 << 32)
        bias = abs(ratio - 0.5)
        marker = " ***" if bias > 0.02 else ""
        print(f"{r:3d} | {ratio:11.4f} | {cs_hws[j]:10.2f} | {cs_zeros_frac[j]:8.5f} | {bias:6.4f}{marker}")

    # Fingerprint diversity
    unique_fp = len(set(fingerprints))
    print(f"\n--- Fingerprint (HW-vector) Diversity ---")
    print(f"  Unique HW-fingerprints: {unique_fp} / {N}")
    print(f"  log2(unique): {np.log2(max(unique_fp, 1)):.2f}")

    # Pairwise fingerprint distance
    dists = []
    rng2 = np.random.RandomState(seed + 100)
    for _ in range(5000):
        i, j = rng2.randint(0, N), rng2.randint(0, N)
        if i == j: continue
        d = sum(abs(fingerprints[i][k] - fingerprints[j][k]) for k in range(48))
        dists.append(d)
    print(f"  Mean L1 distance between fingerprints: {np.mean(dists):.1f}")

    return cs_means, cs_hws, fingerprints


# ============================================================
# S2: THE BRIDGE TEST — carry_sched vs H[7]
# ============================================================

def experiment_S2(N=20000, seed=61):
    print("\n" + "=" * 70)
    print("S2: THE BRIDGE TEST — carry_sched vs H[7]")
    print(f"N={N}, seed={seed}")
    print("=" * 70)

    rng = np.random.RandomState(seed)

    # Collect data
    data = []
    t0 = time.time()
    for i in range(N):
        W16 = [int(rng.randint(0, 1 << 32)) for _ in range(16)]
        t = full_trace(W16)
        # Store: cs_profile (48 words), H7, H6, carry_profile
        data.append({
            'cs': t['cs_profile'],
            'H7': t['H7'],
            'H6': t['H6'],
            'H': t['H'],
            'carry': t['carries'],
        })
        if (i + 1) % 5000 == 0:
            print(f"  {i+1}/{N} ({time.time()-t0:.1f}s)")

    print(f"Collected: {time.time()-t0:.1f}s")

    # Test A: Per-round carry_sched correlation with H[7] bits
    print(f"\n--- Test A: corr(carry_sched[r], H[7]) per round ---")
    print(f"  Pearson correlation of cs[r] (32-bit value) with each bit of H[7]")

    # For efficiency: correlate cs[r] with full H[7] value
    cs_arrays = np.zeros((48, N))
    h7_array = np.zeros(N)
    h6_array = np.zeros(N)

    for i, d in enumerate(data):
        h7_array[i] = d['H7']
        h6_array[i] = d['H6']
        for j in range(48):
            cs_arrays[j, i] = d['cs'][j]

    best_corr_h7 = 0
    best_r_h7 = -1
    best_corr_h6 = 0
    best_r_h6 = -1

    print(f"\n  {'r':>3} | {'corr(cs,H7)':>12} | {'corr(cs,H6)':>12} | {'corr(HW(cs),H7)':>16}")
    print("  " + "-" * 55)

    significant_rounds = []
    for j in range(48):
        r = j + 16
        # Pearson correlation
        if np.std(cs_arrays[j]) > 0:
            corr_h7 = np.corrcoef(cs_arrays[j], h7_array)[0, 1]
            corr_h6 = np.corrcoef(cs_arrays[j], h6_array)[0, 1]
            hw_cs = np.array([hw(int(cs_arrays[j, i])) for i in range(N)], dtype=float)
            corr_hw_h7 = np.corrcoef(hw_cs, h7_array)[0, 1] if np.std(hw_cs) > 0 else 0.0
        else:
            corr_h7 = corr_h6 = corr_hw_h7 = 0.0

        marker = ""
        if abs(corr_h7) > 0.015:
            marker = " ***"
            significant_rounds.append((r, corr_h7))
        if abs(corr_h7) > abs(best_corr_h7):
            best_corr_h7 = corr_h7
            best_r_h7 = r
        if abs(corr_h6) > abs(best_corr_h6):
            best_corr_h6 = corr_h6
            best_r_h6 = r

        if abs(corr_h7) > 0.005 or abs(corr_h6) > 0.005 or r in [16, 17, 62, 63]:
            print(f"  {r:3d} | {corr_h7:+12.5f} | {corr_h6:+12.5f} | {corr_hw_h7:+16.5f}{marker}")

    print(f"\n  Best corr(cs, H7): {best_corr_h7:+.5f} at r={best_r_h7}")
    print(f"  Best corr(cs, H6): {best_corr_h6:+.5f} at r={best_r_h6}")
    print(f"  Significant rounds (|corr|>0.015): {len(significant_rounds)}")

    # Test B: Bit-level — correlate cs[r] bits with H[7] bits
    print(f"\n--- Test B: Bit-level correlation cs[r][bit] → H[7][bit] ---")
    print(f"  Scanning all 48×32 × 32 = 49,152 pairs...")

    best_bit_corr = 0
    best_bit_info = ""
    n_significant_bits = 0

    for j in range(48):
        r = j + 16
        for cs_bit in range(32):
            cs_bit_arr = np.array([(int(data[i]['cs'][j]) >> cs_bit) & 1 for i in range(N)], dtype=float)
            if np.std(cs_bit_arr) < 0.01:
                continue
            for h_bit in range(32):
                h7_bit_arr = np.array([(data[i]['H7'] >> h_bit) & 1 for i in range(N)], dtype=float)
                if np.std(h7_bit_arr) < 0.01:
                    continue
                c = np.corrcoef(cs_bit_arr, h7_bit_arr)[0, 1]
                if abs(c) > abs(best_bit_corr):
                    best_bit_corr = c
                    best_bit_info = f"cs[{r}][b{cs_bit}] → H7[b{h_bit}]"
                if abs(c) > 0.03:
                    n_significant_bits += 1

    print(f"  Best bit-level corr: {best_bit_corr:+.5f} ({best_bit_info})")
    print(f"  Pairs with |corr|>0.03: {n_significant_bits} / 49152")

    # Test C: Total carry_sched sum vs H[7]
    print(f"\n--- Test C: Total carry_sched sum vs H[7] ---")
    cs_total = np.array([sum(d['cs']) for d in data], dtype=float)
    cs_total_hw = np.array([sum(hw(v) for v in d['cs']) for d in data], dtype=float)

    corr_total = np.corrcoef(cs_total, h7_array)[0, 1] if np.std(cs_total) > 0 else 0
    corr_total_hw = np.corrcoef(cs_total_hw, h7_array)[0, 1] if np.std(cs_total_hw) > 0 else 0

    print(f"  corr(sum(cs), H7) = {corr_total:+.5f}")
    print(f"  corr(sum(HW(cs)), H7) = {corr_total_hw:+.5f}")

    # Test D: Last-round carry_sched — direct influence
    print(f"\n--- Test D: Last rounds carry_sched ---")
    print(f"  cs[63] = W_real[63] - W_xor[63]: direct input to final round")

    cs63_arr = cs_arrays[47]  # index 47 = round 63
    cs62_arr = cs_arrays[46]
    cs61_arr = cs_arrays[45]

    for name, arr, label in [
        ("cs[63]", cs63_arr, "H7"),
        ("cs[62]", cs62_arr, "H7"),
        ("cs[61]", cs61_arr, "H7"),
        ("cs[63]", cs63_arr, "H6"),
        ("cs[62]", cs62_arr, "H6"),
    ]:
        target = h7_array if label == "H7" else h6_array
        c = np.corrcoef(arr, target)[0, 1] if np.std(arr) > 0 else 0
        print(f"  corr({name}, {label}) = {c:+.6f}")

    # Test E: cs[63] quartiles → H[7] bit distribution
    print(f"\n--- Test E: cs[63] quartiles → H[7] bit bias ---")
    cs63_vals = [d['cs'][47] for d in data]
    quartiles = np.percentile(cs63_vals, [25, 50, 75])

    for q_name, lo, hi in [("Q1 (lowest)", 0, quartiles[0]),
                            ("Q4 (highest)", quartiles[2], MASK + 1)]:
        subset = [d for d in data if lo <= d['cs'][47] < hi]
        if len(subset) < 100:
            continue
        for bit in [29, 30, 31]:
            p = np.mean([(d['H7'] >> bit) & 1 for d in subset])
            p_base = np.mean([(d['H7'] >> bit) & 1 for d in data])
            delta = p - p_base
            marker = " ***" if abs(delta) > 0.01 else ""
            print(f"  {q_name} (N={len(subset)}): P(H7[b{bit}]=1) = {p:.4f} (delta={delta:+.4f}){marker}")

    return data, significant_rounds, best_bit_corr, best_bit_info


# ============================================================
# S3: carry_sched decomposition — where does signal live?
# ============================================================

def experiment_S3(N=15000, seed=62):
    print("\n" + "=" * 70)
    print("S3: carry_sched Decomposition — Signal Source")
    print(f"N={N}, seed={seed}")
    print("=" * 70)

    rng = np.random.RandomState(seed)

    # For each round r=16..63: how does cs[r] depend on W[0]?
    # cs[r] = W_real[r] - W_xor[r]
    # W_real[r] uses addition, W_xor[r] uses XOR
    # Their difference = accumulated carries in schedule

    # Measure: Jacobian d(cs[r])/d(W[0][bit]) — does flipping W[0] change cs?
    print(f"\n--- Jacobian d(carry_sched)/d(W[0] bits) ---")

    n_jac = 200
    sensitivity = np.zeros((48, 32))  # sensitivity[r-16][bit]

    t0 = time.time()
    for trial in range(n_jac):
        W16 = [int(rng.randint(0, 1 << 32)) for _ in range(16)]
        t_base = full_trace(W16)
        cs_base = t_base['cs_profile']

        for b in range(32):
            W16_flip = list(W16)
            W16_flip[0] ^= (1 << b)
            t_flip = full_trace(W16_flip)
            cs_flip = t_flip['cs_profile']

            for j in range(48):
                if cs_flip[j] != cs_base[j]:
                    sensitivity[j, b] += 1

    sensitivity /= n_jac
    print(f"Computed: {time.time()-t0:.1f}s")

    # Which rounds are most sensitive to W[0]?
    round_sensitivity = np.mean(sensitivity, axis=1)  # avg over bits
    print(f"\n  Per-round sensitivity to W[0] (fraction of bits that change cs):")

    for j in range(48):
        r = j + 16
        s = round_sensitivity[j]
        if s > 0.01:
            bar = '#' * int(s * 100)
            print(f"    r={r:2d}: {s:.4f} {bar}")

    # Which W[0] bits most affect late carry_sched?
    print(f"\n  Per-bit sensitivity of W[0] on cs[60..63]:")
    late_sens = np.mean(sensitivity[44:48, :], axis=0)  # rounds 60-63
    for b in range(32):
        if late_sens[b] > 0.005:
            bar = '#' * int(late_sens[b] * 200)
            print(f"    bit {b:2d}: {late_sens[b]:.4f} {bar}")

    # Key question: does W[0] influence cs[63] differently than cs[16]?
    print(f"\n  cs[16] sensitivity: {np.mean(sensitivity[0]):.4f} (W[16]=W[0]+... directly)")
    print(f"  cs[63] sensitivity: {np.mean(sensitivity[47]):.4f} (47 schedule steps away)")
    print(f"  Ratio: {np.mean(sensitivity[47])/max(np.mean(sensitivity[0]), 1e-10):.3f}")

    return sensitivity


# ============================================================
# S4: Combined carry + carry_sched vs H[7]
# ============================================================

def experiment_S4(N=15000, seed=63):
    print("\n" + "=" * 70)
    print("S4: Combined (carry + carry_sched) vs H[7]")
    print(f"N={N}, seed={seed}")
    print("=" * 70)

    rng = np.random.RandomState(seed)

    # Build feature vector: [carry_profile(64), cs_hw(48)] = 112 features
    # Test correlation of this combined vector with H[7]

    features = []  # N x 112
    h7_vals = []
    h6_vals = []

    t0 = time.time()
    for i in range(N):
        W16 = [int(rng.randint(0, 1 << 32)) for _ in range(16)]
        t = full_trace(W16)

        feat = list(t['carries'])  # 64 binary features
        feat += [hw(v) for v in t['cs_profile']]  # 48 HW features
        features.append(feat)
        h7_vals.append(t['H7'])
        h6_vals.append(t['H6'])

        if (i + 1) % 5000 == 0:
            print(f"  {i+1}/{N} ({time.time()-t0:.1f}s)")

    X = np.array(features, dtype=float)
    y7 = np.array(h7_vals, dtype=float)
    y6 = np.array(h6_vals, dtype=float)

    print(f"\nCollected: {time.time()-t0:.1f}s")
    print(f"Feature matrix: {X.shape}")

    # Per-feature correlation with H[7]
    print(f"\n--- Per-feature correlation with H[7] ---")

    corrs = []
    for f in range(X.shape[1]):
        if np.std(X[:, f]) > 0.01:
            c = np.corrcoef(X[:, f], y7)[0, 1]
        else:
            c = 0.0
        corrs.append(c)

    # Top features by |corr|
    top_idx = sorted(range(len(corrs)), key=lambda i: abs(corrs[i]), reverse=True)[:20]
    print(f"  Top 20 features by |corr(feature, H7)|:")
    for idx in top_idx:
        if idx < 64:
            name = f"carry[{idx}]"
        else:
            name = f"HW(cs[{idx-64+16}])"
        if abs(corrs[idx]) > 0.001:
            print(f"    {name:20s}: {corrs[idx]:+.5f}")

    # Summary: carry vs carry_sched features
    carry_corrs = [abs(corrs[i]) for i in range(64)]
    cs_corrs = [abs(corrs[i]) for i in range(64, 112)]

    print(f"\n  Mean |corr| from carry features:      {np.mean(carry_corrs):.5f}")
    print(f"  Mean |corr| from carry_sched features: {np.mean(cs_corrs):.5f}")
    print(f"  Max |corr| from carry:                 {max(carry_corrs):.5f}")
    print(f"  Max |corr| from carry_sched:           {max(cs_corrs):.5f}")

    # Per-bit H[7] analysis with best features
    print(f"\n--- Per-bit H[7] correlation with best carry_sched round ---")
    # Find which cs round best predicts each H[7] bit
    for h_bit in [28, 29, 30, 31]:
        h_bit_arr = np.array([(h >> h_bit) & 1 for h in h7_vals], dtype=float)
        best_c = 0
        best_f = ""
        for f in range(112):
            if np.std(X[:, f]) < 0.01:
                continue
            c = np.corrcoef(X[:, f], h_bit_arr)[0, 1]
            if abs(c) > abs(best_c):
                best_c = c
                if f < 64:
                    best_f = f"carry[{f}]"
                else:
                    best_f = f"HW(cs[{f-64+16}])"
        print(f"  H7[b{h_bit}]: best corr={best_c:+.5f} from {best_f}")

    # Regression: can linear combination of all features predict H[7] bit 29?
    print(f"\n--- Linear regression: all 112 features → H[7][b29] ---")
    h7_b29 = np.array([(h >> 29) & 1 for h in h7_vals], dtype=float)

    # Simple: correlation of each feature, then weighted sum
    weights = np.array(corrs)
    # Use only features correlated with bit 29
    b29_corrs = []
    for f in range(112):
        if np.std(X[:, f]) > 0.01:
            c = np.corrcoef(X[:, f], h7_b29)[0, 1]
        else:
            c = 0
        b29_corrs.append(c)

    # Prediction score
    score = X @ np.array(b29_corrs)
    pred = (score > np.median(score)).astype(float)
    accuracy = np.mean(pred == h7_b29)
    print(f"  Accuracy predicting H7[b29]: {accuracy:.4f} (random=0.500)")
    print(f"  Advantage: {accuracy - 0.5:+.4f}")

    return corrs, carry_corrs, cs_corrs


# ============================================================
# MAIN
# ============================================================

if __name__ == '__main__':
    print("Carry-Sched Bridge — Stage 1.3")
    print(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print()

    N1 = 10000
    N2 = 15000
    N3 = 10000
    N4 = 12000

    if len(sys.argv) > 1 and sys.argv[1] == '--fast':
        N1, N2, N3, N4 = 5000, 8000, 5000, 6000

    t_start = time.time()

    cs_means, cs_hws, fingerprints = experiment_S1(N=N1)
    data_s2, sig_rounds, best_corr, best_info = experiment_S2(N=N2)
    sensitivity = experiment_S3(N=N3)
    corrs, carry_corrs, cs_corrs = experiment_S4(N=N4)

    total = time.time() - t_start

    print(f"\n{'='*70}")
    print(f"TOTAL TIME: {total:.1f}s")
    print(f"{'='*70}")

    print(f"""
{'='*70}
SYNTHESIS: Carry-Sched Bridge Test
{'='*70}

QUESTION: Does carry_sched (schedule nonlinearity) correlate with H[7]?

carry alone:      corr ≈ 0.005 (ZERO — proven in G3)
carry_sched:      corr = {best_corr:+.5f} ({best_info})

If carry_sched corr > 0.05: BRIDGE EXISTS
   → schedule nonlinearity carries information through W[16..63]
   → bypasses chaotic zone via schedule path (not state path)
   → new mathematics: carry_sched algebra

If carry_sched corr ≈ 0: NO BRIDGE
   → both internal paths (state AND schedule) are opaque
   → SHA-256 is truly random oracle at output
   → only architectural shortcuts (carry[62,63]) work

Mean |corr| carry features:      {np.mean(carry_corrs):.5f}
Mean |corr| carry_sched features: {np.mean(cs_corrs):.5f}
""")
