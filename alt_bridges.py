#!/usr/bin/env python3
"""
Alternative Bridges + Schedule Resonance — Stage 3.6/3.7
==========================================================

Stage 3.5 found: H[5][b31] best predicted by r=5 (not tail!).
Stage 3.4: Maj→H[0] bridge doesn't exist. But what about:

  B1: FULL BRIDGE MAP — corr(raw[r], H[w][bit]) for ALL r×w×bit
      Find ALL bridges, not just the known r63→H[7].
      The H[5]←r5 signal hints at UNDISCOVERED bridges.

  B2: SCHEDULE RESONANCE — W[16]=W[0] creates carry resonance 0.951.
      Does W[r] for r=16..63 create raw[r]↔raw[r-16] coupling?
      If yes → schedule creates a 16-round periodic structure.

  B3: REVERSE ENGINEERING — instead of "which raw[r] predicts H[w]",
      ask "which H[w] combination predicts carry[63]=0?"
      This inverts the problem: output → internal state.
"""

import numpy as np
import time
import sys

MASK = 0xFFFFFFFF
T_VAL = float(1 << 32)

H0 = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]

K = [0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
     0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
     0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
     0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
     0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
     0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
     0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
     0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2]

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
        T1 = raw & MASK; T2 = (Sig0(a) + Maj(a, b, c)) & MASK
        h, g, f = g, f, e; e = (d + T1) & MASK
        d, c, b = c, b, a; a = (T1 + T2) & MASK
    H_out = [(v + iv) & MASK for v, iv in zip([a,b,c,d,e,f,g,h], H0)]
    return raws, H_out, W

def rank_normalize(arr):
    N = len(arr)
    order = np.argsort(arr)
    ranks = np.zeros(N)
    ranks[order] = np.arange(N) / N
    return ranks

def collect(N, seed):
    rng = np.random.RandomState(seed)
    raws = np.zeros((N, 64)); H = np.zeros((N, 8)); W64 = np.zeros((N, 64))
    t0 = time.time()
    for i in range(N):
        W16 = [int(rng.randint(0, 1 << 32)) for _ in range(16)]
        r, h, w = get_raws_and_H(W16)
        raws[i] = r; H[i] = h; W64[i] = w
        if (i+1) % 5000 == 0: print(f"  {i+1}/{N} ({time.time()-t0:.1f}s)")
    return raws, H, W64


# ============================================================
# B1: Full Bridge Map — ALL raw[r] → H[w][bit] correlations
# ============================================================

def experiment_B1(raws, H, N):
    print("=" * 70)
    print(f"B1: Full Bridge Map — corr(raw[r], H[w][bit]) for key r,w,bit")
    print(f"N={N}")
    print("=" * 70)

    # For each H word, test top 4 bits against all 64 rounds
    # This finds ALL bridges, not just the known ones

    print(f"\n--- Heatmap: |corr(raw[r], H[w][b31])| for each round r, word w ---")
    print(f"  r\\H ", end='')
    for w in range(8):
        print(f"  H[{w}]  ", end='')
    print()

    bridge_map = {}  # (r, w, bit) → corr
    best_per_word = {}  # w → (r, corr)

    for w in range(8):
        h_b31 = np.array([(int(H[i, w]) >> 31) & 1 for i in range(N)], dtype=float)
        best_r = -1; best_c = 0
        for r in range(64):
            c = np.corrcoef(raws[:, r], h_b31)[0, 1] if np.std(raws[:, r]) > 0 else 0
            bridge_map[(r, w, 31)] = c
            if abs(c) > abs(best_c):
                best_c = c; best_r = r
        best_per_word[w] = (best_r, best_c)

    # Print condensed: only rounds with |corr| > 0.03 for any word
    active_rounds = set()
    for (r, w, b), c in bridge_map.items():
        if abs(c) > 0.03:
            active_rounds.add(r)
    active_rounds = sorted(active_rounds)

    for r in active_rounds:
        print(f"  {r:2d}  ", end='')
        for w in range(8):
            c = bridge_map.get((r, w, 31), 0)
            if abs(c) > 0.08:
                print(f" {c:+.3f}*", end='')
            elif abs(c) > 0.03:
                print(f" {c:+.3f} ", end='')
            else:
                print(f"    ·   ", end='')
        print()

    print(f"\n  Best predictor round for each H word (bit 31):")
    for w in range(8):
        r, c = best_per_word[w]
        print(f"    H[{w}][b31]: best r={r:2d}, corr={c:+.5f}{'  ★' if abs(c) > 0.05 else ''}")

    # Also check bit 29 for comparison
    print(f"\n  Best predictor round for each H word (bit 29):")
    for w in range(8):
        h_b29 = np.array([(int(H[i, w]) >> 29) & 1 for i in range(N)], dtype=float)
        best_r = -1; best_c = 0
        for r in range(64):
            c = np.corrcoef(raws[:, r], h_b29)[0, 1] if np.std(raws[:, r]) > 0 else 0
            if abs(c) > abs(best_c):
                best_c = c; best_r = r
        print(f"    H[{w}][b29]: best r={r:2d}, corr={c:+.5f}{'  ★' if abs(c) > 0.05 else ''}")

    return bridge_map, best_per_word


# ============================================================
# B2: Schedule Resonance — 16-round periodicity
# ============================================================

def experiment_B2(raws, H, W64, N):
    print("\n" + "=" * 70)
    print(f"B2: Schedule Resonance — W[r] ↔ W[r-16] periodicity")
    print(f"N={N}")
    print("=" * 70)

    # W[16] = sig1(W[14]) + W[9] + sig0(W[1]) + W[0]
    # At pad=0: W[1..15]=0 → W[16]=W[0]. Resonance!
    # For random W[0..15]: W[16] depends on W[0,1,9,14]

    # Test: corr(W[r], W[r-16]) for r=16..63
    print(f"\n  --- Schedule word correlation: corr(W[r], W[r-16]) ---")
    for lag in [16, 7, 15, 2]:
        print(f"\n  Lag={lag}:")
        for r in range(lag, min(64, lag+20)):
            c = np.corrcoef(W64[:, r], W64[:, r-lag])[0, 1]
            marker = " ★★" if abs(c) > 0.3 else (" ★" if abs(c) > 0.1 else "")
            if abs(c) > 0.05 or r < lag + 5:
                print(f"    W[{r:2d}] ↔ W[{r-lag:2d}]: corr={c:+.4f}{marker}")

    # Raw resonance: does raw[r] correlate with raw[r-16]?
    print(f"\n  --- Raw resonance: corr(raw[r], raw[r-16]) ---")
    for r in range(16, 64):
        c = np.corrcoef(raws[:, r], raws[:, r-16])[0, 1]
        if abs(c) > 0.03:
            print(f"    raw[{r:2d}] ↔ raw[{r-16:2d}]: corr={c:+.4f}")

    # Key: does W[r] periodicity create raw[r] periodicity?
    # If W[16]=f(W[0]) and raw[16] depends on W[16]+state →
    # raw[16] correlated with raw[0] through W?
    print(f"\n  --- Cross-period raw: corr(raw[r], raw[r+16]) for r=0..47 ---")
    cross_period = []
    for r in range(48):
        c = np.corrcoef(raws[:, r], raws[:, r+16])[0, 1]
        cross_period.append((r, c))
        if abs(c) > 0.02:
            print(f"    raw[{r:2d}] ↔ raw[{r+16:2d}]: corr={c:+.4f}")

    mean_cross = np.mean([abs(c) for _, c in cross_period])
    print(f"\n  Mean |corr| across 16-round lag: {mean_cross:.5f}")
    print(f"  Expected if independent: ~{1/np.sqrt(N):.5f}")

    return cross_period


# ============================================================
# B3: Reverse — H combination → carry[63]=0
# ============================================================

def experiment_B3(raws, H, N):
    print("\n" + "=" * 70)
    print(f"B3: Reverse Bridge — Which H combination predicts carry[63]?")
    print(f"N={N}")
    print("=" * 70)

    carry63 = (raws[:, 63] >= T_VAL).astype(float)
    # carry63=1 for ~99.9%. We want to predict carry63=0 from H.

    # Per-word correlation with carry63 (as continuous value: use raw63)
    raw63_rank = rank_normalize(raws[:, 63])

    print(f"\n  --- corr(H[w], raw63) — which output words track raw[63]? ---")
    h_corrs = []
    for w in range(8):
        c = np.corrcoef(H[:, w], raw63_rank)[0, 1]
        h_corrs.append(c)
        marker = " ★★" if abs(c) > 0.05 else (" ★" if abs(c) > 0.02 else "")
        print(f"    H[{w}]: corr={c:+.5f}{marker}")

    # Per-bit: find the H bits most correlated with raw63
    print(f"\n  --- Top H bits correlated with raw63 ---")
    all_bit_corrs = []
    for w in range(8):
        for bit in range(32):
            h_bit = np.array([(int(H[i, w]) >> bit) & 1 for i in range(N)], dtype=float)
            c = np.corrcoef(h_bit, raw63_rank)[0, 1] if np.std(h_bit) > 0.01 else 0
            all_bit_corrs.append((w, bit, c))

    all_bit_corrs.sort(key=lambda x: -abs(x[2]))
    print(f"  Top-15 H[w][bit] by |corr(raw63)|:")
    for w, bit, c in all_bit_corrs[:15]:
        print(f"    H[{w}][b{bit:2d}]: corr={c:+.5f}")

    # Optimal linear combination of ALL 256 H bits → raw63
    # Too many features. Use top-32 bits instead.
    top_bits = all_bit_corrs[:32]
    X = np.column_stack([
        np.array([(int(H[i, w]) >> bit) & 1 for i in range(N)], dtype=float)
        for w, bit, _ in top_bits
    ])

    # Train/test
    Nt = N // 2
    X_train = X[:Nt] - X[:Nt].mean(axis=0)
    y_train = raw63_rank[:Nt] - raw63_rank[:Nt].mean()
    X_test = X[Nt:] - X[Nt:].mean(axis=0)
    y_test = raw63_rank[Nt:]

    try:
        beta = np.linalg.lstsq(X_train, y_train, rcond=None)[0]
        pred = X_test @ beta
        corr_pred = np.corrcoef(pred, y_test)[0, 1]

        # Can we predict carry63=0 from H?
        # Bottom 1% of pred → should have carry63=0 more often
        q1 = np.percentile(pred, 1)
        q5 = np.percentile(pred, 5)
        mask_q1 = pred < q1
        mask_q5 = pred < q5
        p_c0_q1 = np.mean(raws[Nt:, 63][mask_q1] < T_VAL) if np.sum(mask_q1) > 0 else 0
        p_c0_q5 = np.mean(raws[Nt:, 63][mask_q5] < T_VAL) if np.sum(mask_q5) > 0 else 0
        p_c0_base = np.mean(raws[Nt:, 63] < T_VAL)

        print(f"\n  --- Reverse predictor (32 H bits → raw63) ---")
        print(f"    corr(predicted, actual raw63): {corr_pred:+.5f}")
        print(f"    P(carry63=0 | pred bottom 1%): {p_c0_q1:.5f} (base: {p_c0_base:.5f})")
        print(f"    P(carry63=0 | pred bottom 5%): {p_c0_q5:.5f} (base: {p_c0_base:.5f})")
        if p_c0_base > 0:
            print(f"    Lift (1%): {p_c0_q1/p_c0_base:.1f}×")
            print(f"    Lift (5%): {p_c0_q5/p_c0_base:.1f}×")

        print(f"\n  ★ REVERSE BRIDGE STATUS:")
        if abs(corr_pred) > 0.15:
            print(f"    STRONG reverse bridge: H predicts raw63 with corr={corr_pred:+.3f}")
        elif abs(corr_pred) > 0.05:
            print(f"    MODERATE reverse bridge: corr={corr_pred:+.3f}")
        else:
            print(f"    WEAK/NO reverse bridge: corr={corr_pred:+.3f}")

    except Exception as e:
        print(f"  Regression failed: {e}")
        corr_pred = 0

    return all_bit_corrs, corr_pred


# ============================================================
# MAIN
# ============================================================

if __name__ == '__main__':
    print("Alternative Bridges + Schedule Resonance — Stage 3.6/3.7")
    print(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print()

    N = 15000
    if len(sys.argv) > 1 and sys.argv[1] == '--fast':
        N = 8000

    t_start = time.time()
    print(f"Collecting data (N={N})...")
    raws, H, W64 = collect(N, seed=1200)

    bridge_map, best_per_word = experiment_B1(raws, H, N)
    cross_period = experiment_B2(raws, H, W64, N)
    bit_corrs, corr_reverse = experiment_B3(raws, H, N)

    total = time.time() - t_start
    print(f"\n{'='*70}")
    print(f"TOTAL TIME: {total:.1f}s")
    print(f"{'='*70}")

    print(f"""
{'='*70}
SYNTHESIS: Alternative Bridges + Schedule Resonance (3.6/3.7)
{'='*70}

B1: FULL BRIDGE MAP
  Which raw[r] → H[w] channels exist beyond the known r63→H[7]?
  Best per word: {dict((w, (r, f'{c:+.3f}')) for w, (r, c) in best_per_word.items())}

B2: SCHEDULE RESONANCE
  W[r] ↔ W[r-16]: does 16-round periodicity create raw coupling?
  Mean cross-period |corr|: {np.mean([abs(c) for _, c in cross_period]):.5f}

B3: REVERSE BRIDGE
  Can H predict carry[63]=0?
  Reverse corr: {corr_reverse:+.4f}
  This is the PASSIVE ATTACK path: output → internal state.
""")
