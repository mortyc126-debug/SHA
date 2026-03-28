#!/usr/bin/env python3
"""
Paradigm Shift — Stage 6: Break the 0.66 ceiling
==================================================

Output-only ceiling: AUC=0.66. To break it, need new paradigm.

Three parallel approaches:

  P1: DIFFERENTIAL — pair (M, M'), analyze δH = H(M) XOR H(M')
      Carry-web predicts: δraw[r] has structure when δW small.
      If δH has detectable structure → differential distinguisher.

  P2: MULTI-QUERY — N hashes from related inputs.
      Each gives AUC=0.66 individually. But N combined?
      If errors independent → N hashes → AUC grows as 1-(1-p)^N.
      If correlated → weaker but still improving.

  P3: ALGEBRAIC — exact bit relations through K-invariants.
      T_CH_INVARIANT: Ch62[b31]=0 when carry[63]=0.
      Ch62 = Ch(H[5]-IV[5], H[6]-IV[6], H[7]-IV[7]) — computable!
      Filter: keep only H where Ch62[b31]=0. Among those, P(carry=0) = 4×.
      This is an ALGEBRAIC pre-filter, not statistical.
"""

import numpy as np
import time, sys

MASK = 0xFFFFFFFF; T_VAL = float(1 << 32)
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

def hc_ss63(rng, steps=200):
    W0=int(rng.randint(0,1<<32)); raws,_=sha_full([W0]+[0]*15); best=raws[63]
    for _ in range(steps):
        b=int(rng.randint(0,32)); W0t=W0^(1<<b); r2,_=sha_full([W0t]+[0]*15)
        if r2[63]<best: W0=W0t; best=r2[63]
    raws,H=sha_full([W0]+[0]*15)
    return W0,raws,H


# ============================================================
# P1: DIFFERENTIAL — δH structure for small δW
# ============================================================

def paradigm_P1(N=10000, seed=5000):
    print("="*70)
    print("P1: DIFFERENTIAL — δH structure for related inputs")
    print(f"N={N}")
    print("="*70)

    rng = np.random.RandomState(seed)

    # For each pair (W, W⊕δ): compute δH = H(W) XOR H(W⊕δ)
    # SHA-256 δH should be random (HW≈128 per word).
    # But with HC (small raw63 for W) → maybe δH has structure?

    # Test A: random pairs with 1-bit difference
    print(f"\n  --- Test A: Random pairs, δW[0] = 1 bit ---")
    dh_hws = {w: [] for w in range(8)}

    for i in range(N):
        W16 = [int(rng.randint(0,1<<32)) for _ in range(16)]
        _, H1 = sha_full(W16)
        b = int(rng.randint(0, 32))
        W16_flip = list(W16); W16_flip[0] ^= (1 << b)
        _, H2 = sha_full(W16_flip)
        for w in range(8):
            dh_hws[w].append(hw(H1[w] ^ H2[w]))

    print(f"  E[HW(δH[w])] for random 1-bit flip:")
    for w in range(8):
        mean_hw = np.mean(dh_hws[w])
        print(f"    δH[{w}]: {mean_hw:.2f} (random=16.0)")

    # Test B: HC-optimized pair — W from HC, W' = W ⊕ 1
    print(f"\n  --- Test B: HC pair, W from ss63 HC, W' = W ⊕ 2^b ---")
    N_hc = min(N // 5, 2000)
    dh_hc = {w: [] for w in range(8)}
    raw63_diffs = []

    for i in range(N_hc):
        W0, raws1, H1 = hc_ss63(rng, steps=150)
        b = int(rng.randint(0, 32))
        raws2, H2 = sha_full([(W0 ^ (1<<b))] + [0]*15)
        for w in range(8):
            dh_hc[w].append(hw(H1[w] ^ H2[w]))
        raw63_diffs.append(abs(raws1[63] - raws2[63]))

    print(f"  E[HW(δH[w])] for HC 1-bit flip (N={N_hc}):")
    for w in range(8):
        mean_hw = np.mean(dh_hc[w])
        delta = mean_hw - 16.0
        marker = " ★" if abs(delta) > 0.3 else ""
        print(f"    δH[{w}]: {mean_hw:.2f} (Δ={delta:+.2f}){marker}")

    print(f"  E[|δraw63|]: {np.mean(raw63_diffs)/T_VAL:.4f}×T")

    # Test C: Are δH bits correlated across words?
    print(f"\n  --- Test C: δH cross-word correlation (HC pairs) ---")
    for w1, w2 in [(6,7), (5,6), (5,7), (0,7)]:
        c = np.corrcoef(dh_hc[w1], dh_hc[w2])[0,1]
        marker = " ★" if abs(c) > 0.05 else ""
        print(f"    corr(HW(δH[{w1}]), HW(δH[{w2}])): {c:+.4f}{marker}")

    return dh_hc


# ============================================================
# P2: MULTI-QUERY — N hashes combined
# ============================================================

def paradigm_P2(N_queries=50, N_trials=500, seed=5001):
    print("\n"+"="*70)
    print(f"P2: MULTI-QUERY — combine N={N_queries} HC hashes")
    print(f"N_trials={N_trials}")
    print("="*70)

    rng = np.random.RandomState(seed)

    # For each trial: generate N_queries HC hashes → compute aggregate score
    # Compare with N_queries random hashes

    def score_single(H):
        """v7-style score from single hash."""
        s = 0
        # H[6][b31]: negative = more likely HC
        s -= ((H[6]>>31)&1) * 0.2
        s -= ((H[6]>>30)&1) * 0.1
        s += ((H[6]>>29)&1) * 0.15
        s -= ((H[7]>>31)&1) * 0.15
        s += ((H[7]>>30)&1) * 0.1
        s += ((H[7]>>29)&1) * 0.1
        # Products
        s += ((H[6]>>30)&1)*((H[6]>>31)&1) * 0.13
        s += ((H[6]>>29)&1)*((H[6]>>31)&1) * 0.11
        return s

    # HC trials
    hc_agg_scores = []
    rand_agg_scores = []

    t0 = time.time()
    for trial in range(N_trials):
        # HC aggregate
        hc_scores = []
        for q in range(N_queries):
            _, _, H = hc_ss63(rng, steps=100)
            hc_scores.append(score_single(H))
        hc_agg_scores.append(np.mean(hc_scores))

        # Random aggregate
        rand_scores = []
        for q in range(N_queries):
            H_rand = [int(rng.randint(0,1<<32)) for _ in range(8)]
            rand_scores.append(score_single(H_rand))
        rand_agg_scores.append(np.mean(rand_scores))

        if (trial+1) % 100 == 0:
            print(f"  {trial+1}/{N_trials} ({time.time()-t0:.1f}s)")

    hc_agg = np.array(hc_agg_scores)
    rand_agg = np.array(rand_agg_scores)

    # Distinguisher on aggregate
    all_scores = np.concatenate([hc_agg, rand_agg])
    all_labels = np.concatenate([np.ones(N_trials), np.zeros(N_trials)])
    threshold = np.median(all_scores)
    pred = (all_scores > threshold).astype(float)
    accuracy = np.mean(pred == all_labels)
    advantage = accuracy - 0.5

    # AUC
    pos = hc_agg; neg = rand_agg
    all_sc = np.concatenate([pos, neg])
    all_lb = np.concatenate([np.ones(len(pos)), np.zeros(len(neg))])
    sorted_lb = all_lb[np.argsort(all_sc)[::-1]]
    n_pos=len(pos); n_neg=len(neg)
    rank_sum = np.sum(np.where(sorted_lb==1)[0])
    auc = 1.0 - (rank_sum - n_pos*(n_pos-1)/2) / (n_pos*n_neg)

    z = advantage * np.sqrt(2*N_trials)

    print(f"\n  MULTI-QUERY RESULTS (N_queries={N_queries}):")
    print(f"    AUC:       {auc:.4f}")
    print(f"    Accuracy:  {accuracy:.4f}")
    print(f"    Advantage: {advantage:+.4f}")
    print(f"    Z-score:   {z:.1f}")
    print(f"    E[score|HC]:     {np.mean(hc_agg):+.4f} ± {np.std(hc_agg):.4f}")
    print(f"    E[score|random]: {np.mean(rand_agg):+.4f} ± {np.std(rand_agg):.4f}")

    # Scaling: how does AUC change with N_queries?
    print(f"\n  --- Scaling: AUC vs N_queries ---")
    for nq in [1, 5, 10, 20, 50]:
        if nq > N_queries: break
        hc_sub = np.array([np.mean(sc[:nq]) for sc in
                          [hc_agg_scores]]) if nq == N_queries else None
        # Approximate: std shrinks as 1/sqrt(N)
        d = (np.mean(hc_agg) - np.mean(rand_agg))
        std_hc = np.std(hc_agg) * np.sqrt(N_queries / nq)
        std_rand = np.std(rand_agg) * np.sqrt(N_queries / nq)
        std_avg = (std_hc + std_rand) / 2
        sep = d / max(std_avg, 1e-10)
        # AUC ≈ Φ(d/√2 / std)
        from math import erfc, sqrt
        auc_est = 0.5 * erfc(-sep / sqrt(2))
        print(f"    N={nq:3d}: estimated AUC={auc_est:.3f}, separation={sep:.2f}σ")

    return auc


# ============================================================
# P3: ALGEBRAIC — Ch-filter pre-selection
# ============================================================

def paradigm_P3(N=15000, seed=5002):
    print("\n"+"="*70)
    print("P3: ALGEBRAIC — Ch62 pre-filter from output")
    print(f"N={N}")
    print("="*70)

    rng = np.random.RandomState(seed)

    # For RANDOM inputs: compute Ch62 from H, filter by Ch62[b31]=0
    # Among filtered: is P(H7[b31]=1) shifted?
    # For RANDOM ORACLE: Ch62 computed from random H → P(Ch[b31]=0) = 0.5
    # → P(H7[b31]=1 | Ch[b31]=0) should still be 0.5

    # For SHA-256: Ch62 = Ch(e62, f62, g62) where e62,f62,g62 are state values
    # These are NOT independent uniform — they're correlated through SHA rounds.
    # Ch filter may create different bias for SHA vs random.

    sha_H = []; rand_H = []
    t0 = time.time()
    for i in range(N):
        # SHA-256
        W16 = [int(rng.randint(0,1<<32)) for _ in range(16)]
        _, H = sha_full(W16)
        sha_H.append(H)
        # Random oracle
        rand_H.append([int(rng.randint(0,1<<32)) for _ in range(8)])

    print(f"  Collected: {time.time()-t0:.1f}s")

    def compute_ch62(H):
        e62 = (H[5] - H0[5]) & MASK
        f62 = (H[6] - H0[6]) & MASK
        g62 = (H[7] - H0[7]) & MASK
        return Ch(e62, f62, g62)

    # Filter: Ch62[b31] = 0
    print(f"\n  --- Ch62 filter: keep H where Ch62[b31]=0 ---")

    for label, hashes in [("SHA-256", sha_H), ("Random", rand_H)]:
        n_pass = 0
        h7_b31_pass = []
        h7_b29_pass = []
        h6_b31_pass = []

        for H in hashes:
            ch = compute_ch62(H)
            if not ((ch >> 31) & 1):  # Ch[b31]=0
                n_pass += 1
                h7_b31_pass.append((H[7]>>31)&1)
                h7_b29_pass.append((H[7]>>29)&1)
                h6_b31_pass.append((H[6]>>31)&1)

        p_pass = n_pass / N
        p_h7b31 = np.mean(h7_b31_pass) if h7_b31_pass else 0.5
        p_h7b29 = np.mean(h7_b29_pass) if h7_b29_pass else 0.5
        p_h6b31 = np.mean(h6_b31_pass) if h6_b31_pass else 0.5

        print(f"\n  {label}:")
        print(f"    P(Ch62[b31]=0):        {p_pass:.4f}")
        print(f"    N passing:             {n_pass}")
        print(f"    P(H7[b31]=1 | pass):   {p_h7b31:.4f} (delta from 0.5: {p_h7b31-0.5:+.4f})")
        print(f"    P(H7[b29]=1 | pass):   {p_h7b29:.4f}")
        print(f"    P(H6[b31]=1 | pass):   {p_h6b31:.4f}")

    # Double filter: Ch62[b31]=0 AND Ch62[b30]=0
    print(f"\n  --- Double filter: Ch62[b31,b30]=0 ---")
    for label, hashes in [("SHA-256", sha_H), ("Random", rand_H)]:
        n_pass = 0; h7_b31_pass = []
        for H in hashes:
            ch = compute_ch62(H)
            if not ((ch>>31)&1) and not ((ch>>30)&1):
                n_pass += 1
                h7_b31_pass.append((H[7]>>31)&1)
        p_pass = n_pass / N
        p_h7b31 = np.mean(h7_b31_pass) if h7_b31_pass else 0.5
        print(f"  {label}: P(pass)={p_pass:.4f}, N={n_pass}, P(H7[b31]=1|pass)={p_h7b31:.4f}")

    # KEY: does the filter create DIFFERENT P(H7[b31]) for SHA vs Random?
    print(f"\n  ★ KEY: Differential Ch filter effect (SHA minus Random):")
    for filt_name in ["Ch[b31]=0", "Ch[b31,b30]=0"]:
        sha_vals = []; rand_vals = []
        for H in sha_H:
            ch = compute_ch62(H)
            if filt_name == "Ch[b31]=0":
                passed = not ((ch>>31)&1)
            else:
                passed = not ((ch>>31)&1) and not ((ch>>30)&1)
            if passed: sha_vals.append((H[7]>>31)&1)

        for H in rand_H:
            ch = compute_ch62(H)
            if filt_name == "Ch[b31]=0":
                passed = not ((ch>>31)&1)
            else:
                passed = not ((ch>>31)&1) and not ((ch>>30)&1)
            if passed: rand_vals.append((H[7]>>31)&1)

        if sha_vals and rand_vals:
            delta = np.mean(sha_vals) - np.mean(rand_vals)
            z = delta / np.sqrt(1/(4*len(sha_vals)) + 1/(4*len(rand_vals)))
            print(f"    {filt_name}: delta P(H7[b31]) = {delta:+.5f}, Z={z:.2f}")
            if abs(z) > 2:
                print(f"    ★★ ALGEBRAIC DISTINGUISHER DETECTED (Z={z:.1f})")
            elif abs(z) > 1:
                print(f"    ★ Borderline signal")
            else:
                print(f"    No signal (Ch filter same for SHA and random)")

    return


# ============================================================
# MAIN
# ============================================================

if __name__ == '__main__':
    print("Paradigm Shift — Stage 6")
    print(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")

    N = 8000
    if len(sys.argv) > 1 and sys.argv[1] == '--fast':
        N = 5000

    t_start = time.time()
    dh = paradigm_P1(N=N)
    auc_mq = paradigm_P2(N_queries=30, N_trials=300)
    paradigm_P3(N=N)

    total = time.time() - t_start
    print(f"\n{'='*70}")
    print(f"TOTAL TIME: {total:.1f}s")
    print(f"{'='*70}")

    print(f"""
{'='*70}
SYNTHESIS: Paradigm Shift (Stage 6)
{'='*70}

P1 DIFFERENTIAL: Does δH have structure for HC pairs?
P2 MULTI-QUERY:  N=30 hashes combined → AUC={auc_mq:.3f}
P3 ALGEBRAIC:    Ch62 filter → SHA vs Random difference?

Which paradigm breaks the 0.66 ceiling?
""")
