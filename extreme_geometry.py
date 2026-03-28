#!/usr/bin/env python3
"""
Extreme Geometry — Stage 2.2: Theory of carry=0 as extreme event
=================================================================

Stage 2.1 showed: GF(2) linear algebra sees nothing — structure
lives in RARE carry=0 events (P≈0.01%). This is the tail of
raw[r] = SS[r] + K[r] + W[r] distribution.

New framework: carry=0 ⟺ raw[r] < 2^32 ⟺ extreme low value.

Experiments:
  E1: RAW DISTRIBUTION — shape of raw[r] for each round
      Is it uniform? Normal? What determines the tail?
  E2: CONDITIONAL GEOMETRY — given carry[63]=0 (via HC),
      what is the distribution of raw[r] for OTHER rounds?
      Are raw values correlated in the extreme?
  E3: JOINT EXTREMES — P(raw[r1] < T1 AND raw[r2] < T2)
      Is there a copula structure? Tail dependence?
  E4: THE BRIDGE — does raw[63] small predict raw[r] small
      better than carry[63]=0 predicts carry[r]=0?
      (Soft threshold vs hard threshold)
"""

import numpy as np
from collections import defaultdict
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

def message_schedule(W16):
    W = list(W16) + [0] * 48
    for i in range(16, 64):
        W[i] = (sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16]) & MASK
    return W

def full_raws(W16):
    """Returns raw[r] for all 64 rounds (before mod 2^32)."""
    W = message_schedule(W16)
    a, b, c, d, e, f, g, h = H0
    raws = []
    for r in range(64):
        raw = h + Sig1(e) + Ch(e, f, g) + K[r] + W[r]
        raws.append(raw)
        T1 = raw & MASK
        T2 = (Sig0(a) + Maj(a, b, c)) & MASK
        h, g, f = g, f, e
        e = (d + T1) & MASK
        d, c, b = c, b, a
        a = (T1 + T2) & MASK
    return raws

def hc_min_raw63(rng, steps=200):
    W0 = int(rng.randint(0, 1 << 32))
    raws = full_raws([W0] + [0]*15)
    best = raws[63]
    for _ in range(steps):
        b = int(rng.randint(0, 32))
        W0_try = W0 ^ (1 << b)
        raws_try = full_raws([W0_try] + [0]*15)
        if raws_try[63] < best:
            W0 = W0_try
            raws = raws_try
            best = raws_try[63]
    return W0, raws


# ============================================================
# E1: Raw Distribution
# ============================================================

def experiment_E1(N=10000, seed=400):
    print("=" * 70)
    print("E1: Raw Distribution — Shape of raw[r] for each round")
    print(f"N={N}, seed={seed}")
    print("=" * 70)

    rng = np.random.RandomState(seed)
    all_raws = np.zeros((N, 64))
    t0 = time.time()

    for i in range(N):
        W16 = [int(rng.randint(0, 1 << 32)) for _ in range(16)]
        raws = full_raws(W16)
        all_raws[i] = raws
        if (i+1) % 5000 == 0:
            print(f"  {i+1}/{N} ({time.time()-t0:.1f}s)")

    print(f"Collected: {time.time()-t0:.1f}s")

    T = float(1 << 32)  # threshold for carry=0

    print(f"\n  {'r':>3} | {'E[raw]/2^32':>12} | {'std/2^32':>9} | {'P(raw<T)':>9} | {'min/T':>7} | {'skewness':>9}")
    print(f"  {'-'*3}-+-{'-'*12}-+-{'-'*9}-+-{'-'*9}-+-{'-'*7}-+-{'-'*9}")

    # Focus on key rounds
    key_rounds = [0,1,4,5,9,10,13,18,19,24,30,33,40,47,48,51,55,59,62,63]
    for r in key_rounds:
        col = all_raws[:, r]
        mean = np.mean(col)
        std = np.std(col)
        p_below = np.mean(col < T)
        mn = np.min(col)
        # Skewness
        if std > 0:
            skew = np.mean(((col - mean)/std)**3)
        else:
            skew = 0
        print(f"  {r:3d} | {mean/T:12.4f} | {std/T:9.4f} | {p_below:9.5f} | {mn/T:7.3f} | {skew:+9.3f}")

    # Key insight: raw[r] = SS[r] + K[r] + W[r]
    # SS ~ sum of 3 random 32-bit values ≈ N(1.5×2^32, σ≈0.87×2^32)
    # K[r] = constant
    # W[r] ~ uniform 32-bit
    # raw ≈ SS + K + W ≈ N(2.5×2^32, σ≈...)
    # Carry=0 ⟺ raw < 2^32 — this is in the EXTREME LEFT TAIL

    mean_all = np.mean(all_raws[:, 63])
    std_all = np.std(all_raws[:, 63])
    z_threshold = (T - mean_all) / std_all
    print(f"\n  For r=63:")
    print(f"    E[raw] = {mean_all/T:.4f} × 2^32")
    print(f"    std    = {std_all/T:.4f} × 2^32")
    print(f"    Threshold (carry=0) = 1.0 × 2^32")
    print(f"    Z-score of threshold: {z_threshold:+.3f}")
    print(f"    Normal approx P(carry=0): {float(_normal_cdf(z_threshold)):.6f}")
    print(f"    Actual P(carry=0): {np.mean(all_raws[:, 63] < T):.6f}")

    return all_raws


def _normal_cdf(z):
    """Approximate normal CDF."""
    from math import erfc, sqrt
    return 0.5 * erfc(-z / sqrt(2))


# ============================================================
# E2: Conditional Geometry — raw[r] given raw[63] small
# ============================================================

def experiment_E2(N=5000, seed=401):
    print("\n" + "=" * 70)
    print("E2: Conditional Geometry — raw[r] | raw[63] in bottom quantile")
    print(f"N={N} HC attempts, seed={seed}")
    print("=" * 70)

    rng = np.random.RandomState(seed)
    T = float(1 << 32)

    # Collect HC-optimized raws
    hc_raws = []
    c0_raws = []

    t0 = time.time()
    for i in range(N):
        W0, raws = hc_min_raw63(rng, steps=200)
        hc_raws.append(raws)
        if raws[63] < T:
            c0_raws.append(raws)
        if (i+1) % 1000 == 0:
            print(f"  {i+1}/{N}, c0={len(c0_raws)} ({time.time()-t0:.1f}s)")

    print(f"\nCollected: {time.time()-t0:.1f}s")
    print(f"  carry[63]=0: {len(c0_raws)} / {N}")

    hc_arr = np.array(hc_raws)

    # Soft analysis: correlation of raw[63] with raw[r] across ALL HC attempts
    print(f"\n--- Soft correlation: corr(raw[63], raw[r]) over HC population ---")
    print(f"  (This uses continuous values, not binary carry)")
    print()
    print(f"  {'r':>3} | {'corr':>8} | {'meaning':>30}")
    print(f"  {'-'*3}-+-{'-'*8}-+-{'-'*30}")

    significant_soft = []
    for r in range(64):
        c = np.corrcoef(hc_arr[:, 63], hc_arr[:, r])[0, 1]
        meaning = ""
        if abs(c) > 0.1:
            meaning = "★★ STRONG"
            significant_soft.append((r, c))
        elif abs(c) > 0.03:
            meaning = "★ moderate"
            significant_soft.append((r, c))

        if abs(c) > 0.02 or r in [0,1,9,10,18,30,47,48,62,63]:
            print(f"  {r:3d} | {c:+8.4f} | {meaning}")

    print(f"\n  Significant (|corr|>0.03): {len(significant_soft)} rounds")

    # KEY: does soft threshold beat hard threshold?
    if c0_raws:
        print(f"\n--- Hard vs Soft threshold ---")
        c0_arr = np.array(c0_raws)

        # For r=0 (island 1): compare
        for r in [0, 1, 9, 10, 30, 47, 48]:
            # Hard: P(carry[r]=0 | carry[63]=0)
            p_hard = np.mean(c0_arr[:, r] < T)
            # Soft: correlation in continuous space
            c_soft = np.corrcoef(hc_arr[:, 63], hc_arr[:, r])[0, 1] if len(hc_arr) > 10 else 0
            # Soft quantile: P(raw[r] in bottom 10% | raw[63] in bottom 1%)
            q1_63 = np.percentile(hc_arr[:, 63], 1)
            mask_q1 = hc_arr[:, 63] < q1_63
            if np.sum(mask_q1) > 5:
                q10_r_cond = np.percentile(hc_arr[mask_q1, r], 10) if np.sum(mask_q1) > 10 else 0
                p_soft_q = np.mean(hc_arr[mask_q1, r] < np.percentile(hc_arr[:, r], 10))
            else:
                p_soft_q = 0

            print(f"  r={r:2d}: P_hard(c[r]=0|c63=0)={p_hard:.3f}, soft_corr={c_soft:+.3f}, P_soft(q10|q1)={p_soft_q:.3f}")

    return hc_arr, c0_raws


# ============================================================
# E3: Joint Extremes — tail dependence
# ============================================================

def experiment_E3(N=15000, seed=402):
    print("\n" + "=" * 70)
    print("E3: Joint Extremes — Tail Dependence of raw values")
    print(f"N={N}, seed={seed}")
    print("=" * 70)

    rng = np.random.RandomState(seed)
    T = float(1 << 32)

    all_raws = np.zeros((N, 64))
    t0 = time.time()
    for i in range(N):
        W16 = [int(rng.randint(0, 1 << 32)) for _ in range(16)]
        all_raws[i] = full_raws(W16)

    print(f"Collected: {time.time()-t0:.1f}s")

    # For each round: compute normalized value (rank within sample / N)
    # Then measure tail dependence: P(U_r1 < q AND U_r2 < q) / q^2
    # If independent: ratio = 1. If tail-dependent: ratio > 1.

    # Rank transform
    ranks = np.zeros_like(all_raws)
    for r in range(64):
        order = np.argsort(all_raws[:, r])
        ranks[order, r] = np.arange(N) / N

    print(f"\n--- Tail Dependence Coefficient λ(r1, r2) ---")
    print(f"  λ = lim P(U_r1 < q, U_r2 < q) / q^2 as q→0")
    print(f"  λ = 1 if independent, λ > 1 if tail-dependent")
    print()

    q_values = [0.10, 0.05, 0.02]
    pairs = [(0, 63), (1, 63), (9, 63), (10, 63), (18, 63),
             (30, 63), (47, 63), (48, 63), (62, 63),
             (0, 9), (0, 47), (9, 47), (47, 62)]

    print(f"  {'pair':>10}", end='')
    for q in q_values:
        print(f" | {'λ(q='+f'{q:.2f}'+')':>12}", end='')
    print(f" | {'interp':>15}")
    print(f"  {'-'*10}", end='')
    for _ in q_values:
        print(f"-+-{'-'*12}", end='')
    print(f"-+-{'-'*15}")

    tail_deps = {}
    for r1, r2 in pairs:
        lambdas = []
        for q in q_values:
            joint = np.mean((ranks[:, r1] < q) & (ranks[:, r2] < q))
            indep = q * q
            lam = joint / max(indep, 1e-10)
            lambdas.append(lam)

        tail_deps[(r1, r2)] = lambdas
        interp = "STRONG TAIL DEP" if lambdas[-1] > 3 else ("moderate" if lambdas[-1] > 1.5 else "independent")
        marker = " ★★★" if lambdas[-1] > 3 else (" ★" if lambdas[-1] > 1.5 else "")

        print(f"  ({r1:2d},{r2:2d})", end='')
        for lam in lambdas:
            print(f" | {lam:12.2f}", end='')
        print(f" | {interp:>15}{marker}")

    # Cross-round raw correlation matrix (key rounds)
    key_r = [0, 1, 9, 10, 18, 30, 47, 48, 62, 63]
    print(f"\n--- Raw correlation matrix (key rounds) ---")
    sub = all_raws[:, key_r]
    corr_mat = np.corrcoef(sub.T)

    print(f"  {'':>5}", end='')
    for r in key_r:
        print(f"  r={r:2d}", end='')
    print()
    for i, r1 in enumerate(key_r):
        print(f"  r={r1:2d}", end='')
        for j, r2 in enumerate(key_r):
            c = corr_mat[i, j]
            marker = '*' if abs(c) > 0.1 and i != j else ' '
            print(f" {c:+5.2f}{marker}", end='')
        print()

    return all_raws, tail_deps


# ============================================================
# E4: Soft Bridge
# ============================================================

def experiment_E4(N=5000, seed=403):
    print("\n" + "=" * 70)
    print("E4: Soft Bridge — Continuous raw[63] vs H[7] bits")
    print(f"N={N} HC, seed={seed}")
    print("=" * 70)

    rng = np.random.RandomState(seed)
    T = float(1 << 32)

    raw63_vals = []
    h7_vals = []
    h6_vals = []

    t0 = time.time()
    for i in range(N):
        W0, raws = hc_min_raw63(rng, steps=200)
        raw63_vals.append(raws[63])

        # Need H[7] — run full SHA
        W = message_schedule([W0] + [0]*15)
        a, b, c, d, e, f, g, h = H0
        for r in range(64):
            raw = h + Sig1(e) + Ch(e, f, g) + K[r] + W[r]
            t1 = raw & MASK
            t2 = (Sig0(a) + Maj(a, b, c)) & MASK
            h, g, f = g, f, e
            e = (d + t1) & MASK
            d, c, b = c, b, a
            a = (t1 + t2) & MASK
        h7 = (h + H0[7]) & MASK
        h6 = (g + H0[6]) & MASK  # g at end = g[63]
        h7_vals.append(h7)
        h6_vals.append(h6)

        if (i+1) % 1000 == 0:
            print(f"  {i+1}/{N} ({time.time()-t0:.1f}s)")

    print(f"Collected: {time.time()-t0:.1f}s")

    raw63 = np.array(raw63_vals)
    h7 = np.array(h7_vals, dtype=float)

    # Overall correlation
    corr_full = np.corrcoef(raw63, h7)[0, 1]
    print(f"\n  corr(raw[63], H[7]): {corr_full:+.5f}")

    # Per-bit correlation: raw[63] vs H[7][bit]
    print(f"\n  Per-bit: corr(raw[63], H[7][bit])")
    for bit in range(32):
        h7_bit = np.array([(int(v) >> bit) & 1 for v in h7_vals], dtype=float)
        c = np.corrcoef(raw63, h7_bit)[0, 1] if np.std(h7_bit) > 0.01 else 0
        if abs(c) > 0.02:
            print(f"    bit {bit:2d}: {c:+.4f}{'  ★' if abs(c) > 0.05 else ''}")

    # Quantile analysis: stratify by raw[63] deciles
    print(f"\n--- Stratification: raw[63] decile → H[7][b29] ---")
    deciles = np.percentile(raw63, np.arange(0, 101, 10))
    for d in range(10):
        lo, hi = deciles[d], deciles[d+1]
        mask = (raw63 >= lo) & (raw63 < hi)
        if np.sum(mask) > 10:
            p_b29 = np.mean([(int(h7_vals[i]) >> 29) & 1 for i in range(N) if mask[i]])
            p_b30 = np.mean([(int(h7_vals[i]) >> 30) & 1 for i in range(N) if mask[i]])
            p_b31 = np.mean([(int(h7_vals[i]) >> 31) & 1 for i in range(N) if mask[i]])
            raw_mean = np.mean(raw63[mask])
            print(f"  D{d}: raw∈[{lo/T:.3f}, {hi/T:.3f}]×2^32, P(b29)={p_b29:.3f}, P(b30)={p_b30:.3f}, P(b31)={p_b31:.3f}")

    # THE KEY: soft threshold (raw[63] < percentile_k) vs hard (carry=0)
    print(f"\n--- Soft vs Hard threshold for H[7] prediction ---")
    for pct in [1, 2, 5, 10, 20, 50]:
        threshold = np.percentile(raw63, pct)
        mask = raw63 < threshold
        n_pass = np.sum(mask)
        if n_pass >= 5:
            p_b29 = np.mean([(int(h7_vals[i]) >> 29) & 1 for i in range(N) if mask[i]])
            phi_b29 = p_b29 - 0.5
            is_c0 = threshold < T
            label = "(carry=0 included)" if is_c0 else "(above carry threshold)"
            print(f"  bottom {pct:2d}%: N={n_pass:4d}, P(b29)={p_b29:.4f}, phi={phi_b29:+.4f} {label}")

    return raw63_vals, h7_vals


# ============================================================
# MAIN
# ============================================================

if __name__ == '__main__':
    print("Extreme Geometry — Stage 2.2")
    print(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print()

    N1, N2, N3, N4 = 8000, 3000, 10000, 3000
    if len(sys.argv) > 1 and sys.argv[1] == '--fast':
        N1, N2, N3, N4 = 5000, 1500, 5000, 1500

    t_start = time.time()
    all_raws = experiment_E1(N=N1)
    hc_raws, c0_raws = experiment_E2(N=N2)
    all_raws_e3, tail_deps = experiment_E3(N=N3)
    raw63, h7 = experiment_E4(N=N4)

    total = time.time() - t_start
    print(f"\n{'='*70}")
    print(f"TOTAL TIME: {total:.1f}s")
    print(f"{'='*70}")

    print(f"""
{'='*70}
SYNTHESIS: Extreme Geometry (Stage 2.2)
{'='*70}

NEW FRAMEWORK: carry=0 as extreme event in raw[r] distribution.

E1: Raw distribution ~ sum of 3-4 random 32-bit values.
    Mean ≈ 1.5 × 2^32, threshold at 1.0 × 2^32.
    Carry=0 lives at Z ≈ -1.5 standard deviations — not THAT extreme.
    But K[r] shifts the mean, making some rounds harder than others.

E2: Conditional geometry — soft correlations in HC population.
    Continuous raw[63] may correlate with raw[r] even when
    binary carry[63] shows no correlation.

E3: Tail dependence — do extreme lows co-occur?
    If λ > 1 for (r, 63) → extreme events cluster across rounds.
    This would be the mathematical basis for the carry-web.

E4: Soft bridge — does continuous raw[63] predict H[7] bits?
    If soft threshold (bottom 10%) shows phi > 0 → the bridge
    works for a LARGER fraction of inputs, not just carry=0.

IMPLICATION: The "new mathematics" may be CONTINUOUS (real-valued
raw distributions and their tail dependence), not DISCRETE
(binary carry profiles and GF(2) algebra).
""")
