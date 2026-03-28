#!/usr/bin/env python3
"""
Distinguisher v7 — Stage 4.5: Nonlinear Multi-Word
====================================================

Stage 4.2 found: H[6][b30]×H[6][b31] gives corr=+0.129 with raw63.
This PRODUCT of bits is stronger than any single bit.

Build a proper distinguisher using nonlinear features from H[5,6,7].
Compare with methodology v6.0 (AUC=0.976 with HC) on the same data.

Key innovation: bit products capture the carry-propagation signature
in the ADDITION H[w] = state + IV. When carry occurs at bit k,
bits k and k+1 become correlated. Products detect this.
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

def hc_minimize_ss63(rng, steps=200):
    """Hill climb to minimize raw[63]."""
    W0=int(rng.randint(0,1<<32))
    W16=[W0]+[0]*15; raws,_=sha_full(W16); best=raws[63]
    for _ in range(steps):
        b=int(rng.randint(0,32)); W0t=W0^(1<<b)
        r2,_=sha_full([W0t]+[0]*15)
        if r2[63]<best: W0=W0t; best=r2[63]
    raws,H=sha_full([W0]+[0]*15)
    return W0,raws,H


def build_features(H_out):
    """Build nonlinear feature vector from H[5,6,7]."""
    feats = []
    # Single bits (24 features: 3 words × 8 bits)
    for w in [5,6,7]:
        for b in range(24, 32):
            feats.append((H_out[w] >> b) & 1)

    # Pairwise products within H[6] (top 4 bits: 6 pairs)
    for b1 in range(28, 32):
        for b2 in range(b1+1, 32):
            feats.append(((H_out[6]>>b1)&1) * ((H_out[6]>>b2)&1))

    # Pairwise products within H[7]
    for b1 in range(28, 32):
        for b2 in range(b1+1, 32):
            feats.append(((H_out[7]>>b1)&1) * ((H_out[7]>>b2)&1))

    # Cross-word products H[6]×H[7] (top 2 bits each: 4 pairs)
    for b6 in [30, 31]:
        for b7 in [30, 31]:
            feats.append(((H_out[6]>>b6)&1) * ((H_out[7]>>b7)&1))

    # Cross-word H[5]×H[6], H[5]×H[7]
    for b5 in [30, 31]:
        for bx in [30, 31]:
            feats.append(((H_out[5]>>b5)&1) * ((H_out[6]>>bx)&1))
            feats.append(((H_out[5]>>b5)&1) * ((H_out[7]>>bx)&1))

    # Triple products H[6][b31]×H[6][b30]×H[7][b31]
    feats.append(((H_out[6]>>31)&1)*((H_out[6]>>30)&1)*((H_out[7]>>31)&1))
    feats.append(((H_out[6]>>31)&1)*((H_out[6]>>29)&1)*((H_out[7]>>30)&1))
    feats.append(((H_out[6]>>31)&1)*((H_out[7]>>30)&1)*((H_out[7]>>29)&1))

    return np.array(feats, dtype=float)


def experiment_v7(N_sha=20000, N_hc=3000, seed=3000):
    print("="*70)
    print("DISTINGUISHER v7: Nonlinear Multi-Word")
    print(f"N_sha={N_sha} (random SHA), N_hc={N_hc} (HC-optimized)")
    print("="*70)

    rng = np.random.RandomState(seed)
    n_feats = len(build_features([0]*8))
    print(f"  Features per sample: {n_feats}")

    # Phase 1: Collect SHA-256 samples (label=1)
    print(f"\n  Collecting {N_sha} random SHA-256 hashes...")
    t0 = time.time()
    X_sha = np.zeros((N_sha, n_feats))
    raw63_sha = np.zeros(N_sha)
    for i in range(N_sha):
        W16=[int(rng.randint(0,1<<32)) for _ in range(16)]
        raws, H = sha_full(W16)
        X_sha[i] = build_features(H)
        raw63_sha[i] = raws[63]
    print(f"  Done: {time.time()-t0:.1f}s")

    # Phase 2: Collect HC-optimized samples (chosen-prefix, label=1)
    print(f"\n  Collecting {N_hc} HC-optimized samples...")
    t0 = time.time()
    X_hc = np.zeros((N_hc, n_feats))
    raw63_hc = np.zeros(N_hc)
    for i in range(N_hc):
        W0, raws, H = hc_minimize_ss63(rng, steps=200)
        X_hc[i] = build_features(H)
        raw63_hc[i] = raws[63]
        if (i+1)%1000==0: print(f"    {i+1}/{N_hc} ({time.time()-t0:.1f}s)")
    print(f"  Done: {time.time()-t0:.1f}s")

    # Phase 3: Generate random oracle samples (label=0)
    print(f"\n  Generating {N_sha} random oracle (uniform H)...")
    X_random = np.zeros((N_sha, n_feats))
    for i in range(N_sha):
        H_rand = [int(rng.randint(0, 1<<32)) for _ in range(8)]
        X_random[i] = build_features(H_rand)

    # Phase 4: Feature importance — SHA vs Random
    print(f"\n  --- Feature correlations with label (SHA=1 vs Random=0) ---")
    X_both = np.vstack([X_sha, X_random])
    y_both = np.concatenate([np.ones(N_sha), np.zeros(N_sha)])

    feat_corrs = []
    for k in range(n_feats):
        if np.std(X_both[:,k]) > 0.01:
            c = np.corrcoef(X_both[:,k], y_both)[0,1]
        else:
            c = 0
        feat_corrs.append(c)

    # Any significant?
    sorted_feats = sorted(enumerate(feat_corrs), key=lambda x: -abs(x[1]))
    print(f"  Top-10 features (SHA vs Random):")
    for idx, c in sorted_feats[:10]:
        print(f"    feature {idx:3d}: corr={c:+.5f}")

    max_sha_random_corr = max(abs(c) for _, c in sorted_feats)
    print(f"\n  Max |corr| SHA vs Random: {max_sha_random_corr:.5f}")
    if max_sha_random_corr < 0.02:
        print(f"  → SHA-256 output INDISTINGUISHABLE from random (no-HC mode)")
        print(f"  → Confirms: passive distinguisher impossible without chosen-prefix")

    # Phase 5: Feature importance — HC vs Random
    print(f"\n  --- Feature correlations with label (HC=1 vs Random=0) ---")
    N_use = min(N_hc, N_sha)
    X_hc_vs_rand = np.vstack([X_hc[:N_use], X_random[:N_use]])
    y_hc_vs_rand = np.concatenate([np.ones(N_use), np.zeros(N_use)])

    feat_corrs_hc = []
    for k in range(n_feats):
        if np.std(X_hc_vs_rand[:,k]) > 0.01:
            c = np.corrcoef(X_hc_vs_rand[:,k], y_hc_vs_rand)[0,1]
        else:
            c = 0
        feat_corrs_hc.append(c)

    sorted_hc = sorted(enumerate(feat_corrs_hc), key=lambda x: -abs(x[1]))
    print(f"  Top-15 features (HC vs Random):")
    for idx, c in sorted_hc[:15]:
        print(f"    feature {idx:3d}: corr={c:+.5f}{'  ★' if abs(c)>0.02 else ''}")

    # Phase 6: Train/Test distinguisher (balanced classes)
    print(f"\n  --- Train/Test Distinguisher (HC vs Random, balanced) ---")
    # Balance: use N_use of each class
    idx_hc = np.arange(N_use)
    idx_rand = np.arange(N_use, 2*N_use)
    rng_split = np.random.RandomState(seed+999)
    rng_split.shuffle(idx_hc); rng_split.shuffle(idx_rand)
    Nt = N_use // 2
    train_idx = np.concatenate([idx_hc[:Nt], idx_rand[:Nt]])
    test_idx = np.concatenate([idx_hc[Nt:], idx_rand[Nt:]])
    X_tr = X_hc_vs_rand[train_idx]
    y_tr = y_hc_vs_rand[train_idx]
    X_te = X_hc_vs_rand[test_idx]
    y_te = y_hc_vs_rand[test_idx]

    # Center
    mu = X_tr.mean(axis=0)
    X_tr_c = X_tr - mu
    X_te_c = X_te - mu

    try:
        beta = np.linalg.lstsq(X_tr_c, y_tr - y_tr.mean(), rcond=None)[0]
        score_tr = X_tr_c @ beta
        score_te = X_te_c @ beta

        # AUC approximation
        from collections import Counter
        pred_te = score_te > np.median(score_te)
        accuracy = np.mean(pred_te == y_te)
        advantage = accuracy - 0.5

        # True positive rate at various thresholds
        pos_scores = score_te[y_te==1]
        neg_scores = score_te[y_te==0]

        # AUC via rank sum
        all_scores = np.concatenate([pos_scores, neg_scores])
        all_labels = np.concatenate([np.ones(len(pos_scores)), np.zeros(len(neg_scores))])
        sorted_idx = np.argsort(all_scores)[::-1]
        sorted_labels = all_labels[sorted_idx]
        n_pos = len(pos_scores)
        n_neg = len(neg_scores)
        rank_sum = np.sum(np.where(sorted_labels==1)[0])
        auc = 1.0 - (rank_sum - n_pos*(n_pos-1)/2) / (n_pos * n_neg)

        z = advantage * np.sqrt(len(y_te))

        print(f"\n  DISTINGUISHER v7 RESULTS:")
        print(f"    Features:    {n_feats} (linear + products + triples)")
        print(f"    AUC:         {auc:.4f}")
        print(f"    Accuracy:    {accuracy:.4f}")
        print(f"    Advantage:   {advantage:+.4f}")
        print(f"    Z-score:     {z:.1f}")
        print(f"    N_test:      {len(y_te)}")

        # Compare with single-feature baseline
        best_single = max(range(n_feats), key=lambda k: abs(feat_corrs_hc[k]))
        pred_single = (X_te[:, best_single] > np.median(X_te[:, best_single]))
        # For HC: higher feature → more likely HC
        if feat_corrs_hc[best_single] < 0:
            pred_single = 1 - pred_single
        acc_single = np.mean(pred_single == y_te)

        print(f"\n  COMPARISON:")
        print(f"    Single best feature:  accuracy={acc_single:.4f} (advantage={acc_single-0.5:+.4f})")
        print(f"    v7 ({n_feats} features):  accuracy={accuracy:.4f} (advantage={advantage:+.4f})")
        print(f"    Improvement:          {(advantage)/(max(acc_single-0.5, 0.001)):.2f}×")

    except Exception as e:
        print(f"  Regression failed: {e}")
        auc = 0.5; accuracy = 0.5; advantage = 0

    # Phase 7: What did v7 learn?
    print(f"\n  --- What v7 learned: top weights ---")
    try:
        top_weights = sorted(enumerate(beta), key=lambda x: -abs(x[1]))[:10]
        for idx, w in top_weights:
            print(f"    feature {idx:3d}: weight={w:+.5f}, corr={feat_corrs_hc[idx]:+.5f}")
    except:
        pass

    return auc, accuracy, advantage


if __name__ == '__main__':
    print("Distinguisher v7 — Stage 4.5")
    print(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print()

    N_sha, N_hc = 15000, 2000
    if len(sys.argv) > 1 and sys.argv[1] == '--fast':
        N_sha, N_hc = 8000, 1000

    auc, acc, adv = experiment_v7(N_sha=N_sha, N_hc=N_hc)

    print(f"""
{'='*70}
FINAL: Distinguisher v7
{'='*70}

AUC:       {auc:.4f}
Accuracy:  {acc:.4f}
Advantage: {adv:+.4f}

Methodology v6.0: AUC=0.976 (with HC + internal state features)
Our v7:          AUC={auc:.4f} (with HC + output-only nonlinear features)

v7 uses ONLY H[5,6,7] bits and their products — no internal state.
If AUC > 0.6 → nonlinear features improve over linear-only.
If AUC > 0.8 → approaching v6.0 without internal access.
""")
