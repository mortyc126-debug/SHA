#!/usr/bin/env python3
"""
Multi-Round Amplification — Stage 3.5
=======================================

Stage 2.4: rank62×rank63 → 1.51× amplification for H[7][b31].
Stage 3.4: state coupling lift 2.5× for nearest-neighbor pairs.
          cross-block (4,30) lift 2.1× through state.

Now: chain them. Build a multi-round product score using:
  - Tail block: raw[61], raw[62], raw[63]
  - State coupling: leverage carry-window structure
  - Cross-block: (4→30), (9→10→30), (30→31→...→47→48)

Goal: amplify the 1.51× to 3-5× through multi-round products.

Experiments:
  A1: 3-round product: rank(61)×rank(62)×rank(63) → H[7]
  A2: Cross-block chain score: product across blocks
  A3: Optimal score function: which combination maximizes corr with H?
  A4: Practical distinguisher using the optimal score
"""

import numpy as np
import time
import sys

MASK = 0xFFFFFFFF
T_VAL = float(1 << 32)

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

def get_raws_and_H(W16):
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
    H_out = [(v + iv) & MASK for v, iv in zip([a,b,c,d,e,f,g,h], H0)]
    return raws, H_out


def collect_data(N, seed):
    """Collect raws and H for N random inputs."""
    rng = np.random.RandomState(seed)
    all_raws = np.zeros((N, 64))
    all_H = np.zeros((N, 8))
    t0 = time.time()
    for i in range(N):
        W16 = [int(rng.randint(0, 1 << 32)) for _ in range(16)]
        raws, H_out = get_raws_and_H(W16)
        all_raws[i] = raws
        all_H[i] = H_out
        if (i+1) % 5000 == 0:
            print(f"  {i+1}/{N} ({time.time()-t0:.1f}s)")
    return all_raws, all_H


def rank_normalize(arr):
    """Rank normalize to [0,1]. Small values → small ranks."""
    N = len(arr)
    order = np.argsort(arr)
    ranks = np.zeros(N)
    ranks[order] = np.arange(N) / N
    return ranks


# ============================================================
# A1: 3-round product
# ============================================================

def experiment_A1(all_raws, all_H, N):
    print("=" * 70)
    print("A1: Multi-Round Product Scores → H[7]")
    print(f"N={N}")
    print("=" * 70)

    # Rank normalize each round
    ranks = {}
    for r in range(64):
        ranks[r] = rank_normalize(all_raws[:, r])

    h7_b31 = np.array([(int(all_H[i, 7]) >> 31) & 1 for i in range(N)], dtype=float)
    h7_b29 = np.array([(int(all_H[i, 7]) >> 29) & 1 for i in range(N)], dtype=float)
    h6_b31 = np.array([(int(all_H[i, 6]) >> 31) & 1 for i in range(N)], dtype=float)
    h6_b29 = np.array([(int(all_H[i, 6]) >> 29) & 1 for i in range(N)], dtype=float)

    # Test various product scores
    scores = {}

    # Single rounds
    scores['r63'] = ranks[63]
    scores['r62'] = ranks[62]
    scores['r61'] = ranks[61]

    # 2-round products (tail block)
    scores['r62×r63'] = ranks[62] * ranks[63]
    scores['r61×r63'] = ranks[61] * ranks[63]
    scores['r61×r62'] = ranks[61] * ranks[62]

    # 3-round product (tail block)
    scores['r61×r62×r63'] = ranks[61] * ranks[62] * ranks[63]

    # Cross-block products
    scores['r48×r63'] = ranks[48] * ranks[63]
    scores['r47×r63'] = ranks[47] * ranks[63]
    scores['r30×r63'] = ranks[30] * ranks[63]
    scores['r9×r63'] = ranks[9] * ranks[63]

    # Multi-block chains
    scores['r9×r10×r63'] = ranks[9] * ranks[10] * ranks[63]
    scores['r30×r31×r63'] = ranks[30] * ranks[31] * ranks[63]
    scores['r47×r48×r63'] = ranks[47] * ranks[48] * ranks[63]

    # Full chain: one round per block
    scores['r9×r30×r47×r63'] = ranks[9] * ranks[30] * ranks[47] * ranks[63]
    scores['r0×r9×r30×r47×r63'] = ranks[0] * ranks[9] * ranks[30] * ranks[47] * ranks[63]

    # Sum-based (additive instead of multiplicative)
    scores['sum(r61,r62,r63)'] = ranks[61] + ranks[62] + ranks[63]
    scores['sum(r9,r30,r47,r63)'] = ranks[9] + ranks[30] + ranks[47] + ranks[63]

    # Log-product (= sum of log-ranks, more stable)
    eps = 1e-10
    scores['logprod(r62,r63)'] = np.log(ranks[62] + eps) + np.log(ranks[63] + eps)
    scores['logprod(r61,r62,r63)'] = np.log(ranks[61]+eps) + np.log(ranks[62]+eps) + np.log(ranks[63]+eps)
    scores['logprod(4blocks)'] = np.log(ranks[9]+eps) + np.log(ranks[30]+eps) + np.log(ranks[47]+eps) + np.log(ranks[63]+eps)

    print(f"\n  {'Score':>30} | {'corr(H7)':>9} | {'corr(b31)':>10} | {'corr(b29)':>10} | {'corr(H6[b31])':>14}")
    print(f"  {'-'*30}-+-{'-'*9}-+-{'-'*10}-+-{'-'*10}-+-{'-'*14}")

    results = []
    for name, sc in scores.items():
        c_h7 = np.corrcoef(sc, all_H[:, 7].astype(float))[0, 1]
        c_b31 = np.corrcoef(sc, h7_b31)[0, 1]
        c_b29 = np.corrcoef(sc, h7_b29)[0, 1]
        c_h6b31 = np.corrcoef(sc, h6_b31)[0, 1]
        results.append((name, c_h7, c_b31, c_b29, c_h6b31))
        marker = " ★★" if abs(c_b31) > 0.12 else (" ★" if abs(c_b31) > 0.09 else "")
        print(f"  {name:>30} | {c_h7:+9.4f} | {c_b31:+10.4f} | {c_b29:+10.4f} | {c_h6b31:+14.4f}{marker}")

    # Sort by |corr(b31)|
    results.sort(key=lambda x: -abs(x[2]))
    print(f"\n  TOP-5 by |corr(H7[b31])|:")
    for name, c_h7, c_b31, c_b29, c_h6b31 in results[:5]:
        amplif = abs(c_b31) / abs(results[-1][2]) if abs(results[-1][2]) > 0.001 else 0
        print(f"    {name:>30}: corr(b31)={c_b31:+.4f}")

    # Amplification vs single r63
    base_b31 = abs([r for r in results if r[0] == 'r63'][0][2])
    print(f"\n  Amplification vs single r63 (|corr(b31)|={base_b31:.4f}):")
    for name, c_h7, c_b31, c_b29, c_h6b31 in results[:8]:
        amp = abs(c_b31) / max(base_b31, 0.001)
        print(f"    {name:>30}: {amp:.2f}×")

    return scores, results


# ============================================================
# A2: Optimal linear combination
# ============================================================

def experiment_A2(all_raws, all_H, N):
    print("\n" + "=" * 70)
    print("A2: Optimal Linear Combination of raw[r] → H[7][b31]")
    print(f"N={N}")
    print("=" * 70)

    h7_b31 = np.array([(int(all_H[i, 7]) >> 31) & 1 for i in range(N)], dtype=float)

    # Use last 8 rounds as features
    feature_rounds = [56, 57, 58, 59, 60, 61, 62, 63]
    X = np.column_stack([rank_normalize(all_raws[:, r]) for r in feature_rounds])

    # Per-feature correlation
    print(f"\n  Per-feature corr with H7[b31]:")
    for i, r in enumerate(feature_rounds):
        c = np.corrcoef(X[:, i], h7_b31)[0, 1]
        print(f"    rank(raw[{r}]): {c:+.5f}")

    # Simple linear regression (closed form)
    # β = (X^T X)^-1 X^T y
    X_centered = X - X.mean(axis=0)
    y_centered = h7_b31 - h7_b31.mean()

    try:
        beta = np.linalg.lstsq(X_centered, y_centered, rcond=None)[0]
        y_pred = X_centered @ beta

        corr_opt = np.corrcoef(y_pred, h7_b31)[0, 1]
        print(f"\n  Optimal linear combination (8 features): corr={corr_opt:+.5f}")
        print(f"  Weights:")
        for i, r in enumerate(feature_rounds):
            print(f"    raw[{r}]: {beta[i]:+.5f}")
    except:
        corr_opt = 0
        print(f"  Linear regression failed")

    # Wider: use all carry-window rounds + tail
    wide_rounds = [0, 1, 4, 5, 9, 10, 18, 19, 30, 31, 47, 48, 49, 61, 62, 63]
    X_wide = np.column_stack([rank_normalize(all_raws[:, r]) for r in wide_rounds])
    X_wide_c = X_wide - X_wide.mean(axis=0)

    try:
        beta_wide = np.linalg.lstsq(X_wide_c, y_centered, rcond=None)[0]
        y_pred_wide = X_wide_c @ beta_wide
        corr_wide = np.corrcoef(y_pred_wide, h7_b31)[0, 1]
        print(f"\n  Wide linear combination ({len(wide_rounds)} features): corr={corr_wide:+.5f}")

        # Top weights
        sorted_weights = sorted(zip(wide_rounds, beta_wide), key=lambda x: -abs(x[1]))
        print(f"  Top-5 weights:")
        for r, w in sorted_weights[:5]:
            print(f"    raw[{r}]: {w:+.5f}")

        # Amplification
        base_corr = abs(np.corrcoef(rank_normalize(all_raws[:, 63]), h7_b31)[0, 1])
        print(f"\n  Amplification: {abs(corr_wide)/max(base_corr, 0.001):.2f}× vs single r63")
    except:
        corr_wide = 0

    # Multi-target: optimize for joint (H7[b31], H6[b31], H7[b29], H6[b29])
    print(f"\n  --- Multi-target: best feature for EACH output bit ---")
    targets = {
        'H7[b31]': np.array([(int(all_H[i, 7]) >> 31) & 1 for i in range(N)], dtype=float),
        'H7[b29]': np.array([(int(all_H[i, 7]) >> 29) & 1 for i in range(N)], dtype=float),
        'H6[b31]': np.array([(int(all_H[i, 6]) >> 31) & 1 for i in range(N)], dtype=float),
        'H6[b29]': np.array([(int(all_H[i, 6]) >> 29) & 1 for i in range(N)], dtype=float),
        'H5[b31]': np.array([(int(all_H[i, 5]) >> 31) & 1 for i in range(N)], dtype=float),
    }

    for tname, target in targets.items():
        tc = target - target.mean()
        try:
            beta_t = np.linalg.lstsq(X_wide_c, tc, rcond=None)[0]
            pred_t = X_wide_c @ beta_t
            corr_t = np.corrcoef(pred_t, target)[0, 1]
            best_r = wide_rounds[np.argmax(np.abs(beta_t))]
            print(f"    {tname}: optimal_corr={corr_t:+.5f}, best_round=r{best_r}")
        except:
            print(f"    {tname}: failed")

    return corr_opt, corr_wide


# ============================================================
# A3: Practical distinguisher using optimal score
# ============================================================

def experiment_A3(all_raws, all_H, N):
    print("\n" + "=" * 70)
    print("A3: Practical Distinguisher — Train/Test Split")
    print(f"N={N}")
    print("=" * 70)

    # Split
    N_train = N // 2
    N_test = N - N_train

    wide_rounds = [0, 1, 4, 5, 9, 10, 18, 19, 30, 31, 47, 48, 49, 61, 62, 63]

    # Train
    X_train = np.column_stack([rank_normalize(all_raws[:N_train, r]) for r in wide_rounds])
    h7_b31_train = np.array([(int(all_H[i, 7]) >> 31) & 1 for i in range(N_train)], dtype=float)
    X_train_c = X_train - X_train.mean(axis=0)
    y_train_c = h7_b31_train - h7_b31_train.mean()

    try:
        beta = np.linalg.lstsq(X_train_c, y_train_c, rcond=None)[0]
    except:
        print("  Training failed")
        return

    # Test
    X_test = np.column_stack([rank_normalize(all_raws[N_train:, r]) for r in wide_rounds])
    h7_b31_test = np.array([(int(all_H[i, 7]) >> 31) & 1 for i in range(N_train, N)], dtype=float)
    X_test_c = X_test - X_test.mean(axis=0)

    y_pred_test = X_test_c @ beta
    corr_test = np.corrcoef(y_pred_test, h7_b31_test)[0, 1]

    # Threshold at median
    threshold = np.median(y_pred_test)
    pred_class = (y_pred_test > threshold).astype(float)
    accuracy = np.mean(pred_class == h7_b31_test)
    advantage = accuracy - 0.5

    print(f"\n  Train corr(optimal, H7[b31]): {np.corrcoef(X_train_c @ beta, h7_b31_train)[0,1]:+.5f}")
    print(f"  Test corr(optimal, H7[b31]):  {corr_test:+.5f}")
    print(f"  Test accuracy:                {accuracy:.4f}")
    print(f"  Test advantage:               {advantage:+.4f}")

    # Compare with single r63
    r63_test = rank_normalize(all_raws[N_train:, 63])
    corr_r63_test = np.corrcoef(r63_test, h7_b31_test)[0, 1]
    acc_r63 = np.mean((r63_test < np.median(r63_test)).astype(float) == h7_b31_test)
    print(f"\n  Single r63 test corr:         {corr_r63_test:+.5f}")
    print(f"  Single r63 test accuracy:     {acc_r63:.4f}")
    print(f"  Improvement:                  {(accuracy-acc_r63)/max(acc_r63-0.5, 0.001)*100:+.0f}% over r63")

    # Z-score
    z = advantage * np.sqrt(N_test)
    print(f"  Z-score: {z:.1f}")

    # THIS IS THE KEY: does multi-round amplification ACTUALLY work?
    print(f"\n  ★ MULTI-ROUND AMPLIFICATION RESULT:")
    print(f"    Single raw[63]:     corr={abs(corr_r63_test):.4f}, accuracy={acc_r63:.4f}")
    print(f"    Optimal 16-feature: corr={abs(corr_test):.4f}, accuracy={accuracy:.4f}")
    print(f"    Amplification:      {abs(corr_test)/max(abs(corr_r63_test), 0.001):.2f}×")

    return accuracy, advantage, corr_test


# ============================================================
# MAIN
# ============================================================

if __name__ == '__main__':
    print("Multi-Round Amplification — Stage 3.5")
    print(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print()

    N = 20000
    if len(sys.argv) > 1 and sys.argv[1] == '--fast':
        N = 10000

    t_start = time.time()

    print(f"Collecting data (N={N})...")
    all_raws, all_H = collect_data(N, seed=1100)

    scores, results = experiment_A1(all_raws, all_H, N)
    corr_opt, corr_wide = experiment_A2(all_raws, all_H, N)
    accuracy, advantage, corr_test = experiment_A3(all_raws, all_H, N)

    total = time.time() - t_start
    print(f"\n{'='*70}")
    print(f"TOTAL TIME: {total:.1f}s")
    print(f"{'='*70}")

    print(f"""
{'='*70}
SYNTHESIS: Multi-Round Amplification (Stage 3.5)
{'='*70}

PRODUCT SCORES: rank products across multiple rounds.
OPTIMAL LINEAR: regression on 16 carry-window features.
PRACTICAL TEST: train/test split accuracy.

Single raw[63]:      baseline
Optimal 16-feature:  amplification = {abs(corr_test)/max(abs(np.corrcoef(rank_normalize(all_raws[:N//2, 63]), np.array([(int(all_H[i, 7]) >> 31) & 1 for i in range(N//2, N)], dtype=float))[0,1]), 0.001):.2f}×

KEY QUESTION: Does combining multiple rounds ACTUALLY help?
If amplification > 1.5× → multi-round attack viable.
If amplification ≈ 1.0× → single round is optimal (no new math needed).
""")
