#!/usr/bin/env python3
"""
Deep Nonlinear + State Predictor — Stage 5
=============================================

v7: AUC=0.661 with 51 features (output-only, HC).
Gap to v6.0 (AUC=0.976) = internal state access.

Two paths to close the gap:

  5.1: DEEPER FEATURES — expand to H[4..7], add 3-way products,
       bit-slice features (top byte as integer), and word-level stats.
       Goal: squeeze more from output bits alone.

  5.2: SMARTER HC — instead of minimizing SS[63] blindly,
       use the composite state coupling score to guide HC.
       Minimize sum of raw[r] for coupling rounds {9,10,30,31,47,48,62,63}.

  5.3: COMBINE — v8 distinguisher using deep features + smart HC.
"""

import numpy as np
import time
import sys

MASK = 0xFFFFFFFF
T_VAL = float(1 << 32)

H0 = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
      0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]
K = [0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
     0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
     0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
     0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
     0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
     0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
     0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
     0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2]

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&MASK
def sig0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def sig1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def Sig0(x): return rotr(x,2)^rotr(x,13)^rotr(x,22)
def Sig1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def Ch(e,f,g): return (e&f)^(~e&g)&MASK
def Maj(a,b,c): return (a&b)^(a&c)^(b&c)
def hw(x): return bin(x&MASK).count('1')

def msg_sched(W16):
    W=list(W16)+[0]*48
    for i in range(16,64): W[i]=(sig1(W[i-2])+W[i-7]+sig0(W[i-15])+W[i-16])&MASK
    return W

def sha_full(W16):
    W=msg_sched(W16); a,b,c,d,e,f,g,h=H0; raws=[]
    for r in range(64):
        raw=h+Sig1(e)+Ch(e,f,g)+K[r]+W[r]; raws.append(raw)
        T1=raw&MASK; T2=(Sig0(a)+Maj(a,b,c))&MASK
        h,g,f=g,f,e; e=(d+T1)&MASK; d,c,b=c,b,a; a=(T1+T2)&MASK
    return raws, [(v+iv)&MASK for v,iv in zip([a,b,c,d,e,f,g,h],H0)]


def build_deep_features(H):
    """Build feature vector: deep nonlinear from H[4..7]."""
    feats = []

    # Layer 1: single bits H[4..7], top 8 bits each (32 features)
    for w in [4,5,6,7]:
        for b in range(24, 32):
            feats.append((H[w] >> b) & 1)

    # Layer 2: top-byte as normalized integer (4 features)
    for w in [4,5,6,7]:
        feats.append((H[w] >> 24) / 255.0)

    # Layer 3: HW of top byte (4 features)
    for w in [4,5,6,7]:
        feats.append(hw(H[w] >> 24) / 8.0)

    # Layer 4: adjacent bit products within each word (3 pairs × 4 words = 12)
    for w in [5,6,7]:
        for b in range(29, 32):
            feats.append(((H[w]>>b)&1) * ((H[w]>>(b-1))&1))

    # Layer 5: cross-word products H[w1][b1]×H[w2][b2] (key pairs)
    cross_pairs = [
        (6,31, 7,31), (6,31, 7,30), (6,30, 7,31), (6,30, 7,30),
        (6,31, 6,30), (6,31, 6,29), (6,30, 6,29),
        (7,31, 7,30), (7,31, 7,29), (7,30, 7,29),
        (5,31, 6,31), (5,31, 7,31), (5,30, 6,30),
        (6,31, 5,30), (6,29, 7,30),
    ]
    for w1,b1,w2,b2 in cross_pairs:
        feats.append(((H[w1]>>b1)&1) * ((H[w2]>>b2)&1))

    # Layer 6: triple products (key)
    feats.append(((H[6]>>31)&1)*((H[6]>>30)&1)*((H[7]>>31)&1))
    feats.append(((H[6]>>31)&1)*((H[6]>>29)&1)*((H[7]>>31)&1))
    feats.append(((H[6]>>31)&1)*((H[7]>>30)&1)*((H[7]>>29)&1))
    feats.append(((H[5]>>31)&1)*((H[6]>>31)&1)*((H[7]>>31)&1))

    # Layer 7: XOR features (carry signature pattern)
    feats.append(((H[6]>>31)&1) ^ ((H[7]>>31)&1))  # sign difference
    feats.append(((H[6]>>30)&1) ^ ((H[6]>>31)&1))  # adjacent XOR in H[6]
    feats.append(((H[7]>>30)&1) ^ ((H[7]>>31)&1))

    return np.array(feats, dtype=float)


def hc_smart(rng, steps=200, target='ss63'):
    """Hill climb with different targets."""
    W0 = int(rng.randint(0, 1 << 32))
    W16 = [W0] + [0]*15
    raws, _ = sha_full(W16)

    if target == 'ss63':
        best = raws[63]
    elif target == 'composite':
        # Minimize sum of raw at coupling rounds
        coupling_rounds = [9, 10, 30, 31, 47, 48, 62, 63]
        weights = [0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 2.0, 3.0]  # emphasize late rounds
        best = sum(raws[r]*w for r,w in zip(coupling_rounds, weights))
    elif target == 'tail3':
        best = raws[61] + raws[62] + 2*raws[63]
    else:
        best = raws[63]

    for _ in range(steps):
        b = int(rng.randint(0, 32))
        W0t = W0 ^ (1 << b)
        r2, _ = sha_full([W0t]+[0]*15)
        if target == 'ss63':
            val = r2[63]
        elif target == 'composite':
            val = sum(r2[r]*w for r,w in zip(coupling_rounds, weights))
        elif target == 'tail3':
            val = r2[61] + r2[62] + 2*r2[63]
        else:
            val = r2[63]

        if val < best:
            W0 = W0t; best = val

    raws, H = sha_full([W0]+[0]*15)
    return W0, raws, H


# ============================================================
# MAIN EXPERIMENT
# ============================================================

def run_experiment(N_hc=1500, N_rand=5000, seed=4000):
    print("="*70)
    print("Stage 5: Deep Nonlinear + Smart HC → Distinguisher v8")
    print(f"N_hc={N_hc}, N_rand={N_rand}")
    print("="*70)

    rng = np.random.RandomState(seed)
    n_feats = len(build_deep_features([0]*8))
    print(f"  Deep features: {n_feats}")

    # Collect HC samples with different strategies
    strategies = ['ss63', 'composite', 'tail3']
    hc_data = {}

    for strat in strategies:
        X_hc = []; raw63_hc = []
        t0 = time.time()
        for i in range(N_hc):
            _, raws, H = hc_smart(rng, steps=200, target=strat)
            X_hc.append(build_deep_features(H))
            raw63_hc.append(raws[63])
        hc_data[strat] = {
            'X': np.array(X_hc),
            'raw63': np.array(raw63_hc),
            'p_c0': np.mean(np.array(raw63_hc) < T_VAL),
        }
        print(f"  HC({strat}): P(c63=0)={hc_data[strat]['p_c0']:.4f}, {time.time()-t0:.1f}s")

    # Random oracle baseline
    X_rand = []
    for i in range(N_rand):
        H_rand = [int(rng.randint(0, 1<<32)) for _ in range(8)]
        X_rand.append(build_deep_features(H_rand))
    X_rand = np.array(X_rand)

    # Test each strategy: HC vs Random distinguisher
    print(f"\n  --- Distinguisher comparison ---")
    print(f"  {'Strategy':>12} | {'AUC':>6} | {'Accuracy':>9} | {'Advantage':>10} | {'P(c63=0)':>9}")
    print(f"  {'-'*12}-+-{'-'*6}-+-{'-'*9}-+-{'-'*10}-+-{'-'*9}")

    best_strat = None; best_auc = 0

    for strat in strategies:
        N_use = min(N_hc, N_rand)
        X = np.vstack([hc_data[strat]['X'][:N_use], X_rand[:N_use]])
        y = np.concatenate([np.ones(N_use), np.zeros(N_use)])

        # Shuffle and split
        perm = rng.permutation(2*N_use)
        X = X[perm]; y = y[perm]
        Nt = N_use  # half for each
        X_tr, y_tr = X[:Nt], y[:Nt]
        X_te, y_te = X[Nt:], y[Nt:]

        mu = X_tr.mean(0)
        X_tr_c = X_tr - mu; X_te_c = X_te - mu

        try:
            beta = np.linalg.lstsq(X_tr_c, y_tr - y_tr.mean(), rcond=None)[0]
            score_te = X_te_c @ beta
            pred = score_te > np.median(score_te)
            accuracy = np.mean(pred == y_te)
            advantage = accuracy - 0.5

            # AUC
            pos = score_te[y_te==1]; neg = score_te[y_te==0]
            n_pos = len(pos); n_neg = len(neg)
            if n_pos > 0 and n_neg > 0:
                all_sc = np.concatenate([pos, neg])
                all_lb = np.concatenate([np.ones(n_pos), np.zeros(n_neg)])
                sorted_lb = all_lb[np.argsort(all_sc)[::-1]]
                rank_sum = np.sum(np.where(sorted_lb==1)[0])
                auc = 1.0 - (rank_sum - n_pos*(n_pos-1)/2) / (n_pos*n_neg)
            else:
                auc = 0.5

            p_c0 = hc_data[strat]['p_c0']
            print(f"  {strat:>12} | {auc:6.3f} | {accuracy:9.4f} | {advantage:+10.4f} | {p_c0:9.4f}")

            if auc > best_auc:
                best_auc = auc; best_strat = strat

        except:
            print(f"  {strat:>12} | FAIL")

    # Feature importance for best strategy
    if best_strat:
        print(f"\n  BEST STRATEGY: {best_strat} (AUC={best_auc:.3f})")

        N_use = min(N_hc, N_rand)
        X_best = np.vstack([hc_data[best_strat]['X'][:N_use], X_rand[:N_use]])
        y_best = np.concatenate([np.ones(N_use), np.zeros(N_use)])

        corrs = []
        for k in range(n_feats):
            if np.std(X_best[:,k]) > 0.01:
                c = np.corrcoef(X_best[:,k], y_best)[0,1]
            else:
                c = 0
            corrs.append(c)

        sorted_f = sorted(enumerate(corrs), key=lambda x: -abs(x[1]))
        print(f"\n  Top-10 features for {best_strat}:")
        for idx, c in sorted_f[:10]:
            print(f"    feature {idx:3d}: corr={c:+.4f}")

    # FINAL: compare with v7
    print(f"\n  {'='*50}")
    print(f"  COMPARISON:")
    print(f"    v7 (51 features, ss63 HC):    AUC=0.661")
    print(f"    v8 ({n_feats} features, {best_strat} HC): AUC={best_auc:.3f}")
    if best_auc > 0.661:
        print(f"    IMPROVEMENT: +{best_auc - 0.661:.3f} AUC ({(best_auc-0.661)/0.661*100:+.1f}%)")
    else:
        print(f"    No improvement over v7")

    return best_auc, best_strat


if __name__ == '__main__':
    print(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")

    N_hc, N_rand = 1500, 5000
    if len(sys.argv) > 1 and sys.argv[1] == '--fast':
        N_hc, N_rand = 800, 3000

    auc, strat = run_experiment(N_hc=N_hc, N_rand=N_rand)

    print(f"""
{'='*70}
SYNTHESIS: Stage 5
{'='*70}

Best AUC: {auc:.3f} (strategy: {strat})
v7 baseline: 0.661
v6.0 (methodology): 0.976

Nonlinear deep features + smart HC target.
""")
