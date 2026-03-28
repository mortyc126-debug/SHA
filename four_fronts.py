#!/usr/bin/env python3
"""
Four Fronts — Stage 4: All open directions simultaneously
===========================================================

4.1: STATE COUPLING ATTACK — use carry[r_early] to predict carry[63]
4.2: MULTI-WORD TAIL DISTINGUISHER — H[5]↔H[6] corr=+0.201
4.3: CONTINUOUS HENSEL — raw[62] mod 2^k → raw[63] mod 2^k
4.4: MULTI-BLOCK — carry structure propagation across blocks
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
    H_out=[(v+iv)&MASK for v,iv in zip([a,b,c,d,e,f,g,h],H0)]
    return raws,H_out

def sha_compress(IV, W16):
    """SHA-256 compression with custom IV."""
    W=msg_sched(W16); a,b,c,d,e,f,g,h=IV; raws=[]
    for r in range(64):
        raw=h+Sig1(e)+Ch(e,f,g)+K[r]+W[r]; raws.append(raw)
        T1=raw&MASK; T2=(Sig0(a)+Maj(a,b,c))&MASK
        h,g,f=g,f,e; e=(d+T1)&MASK; d,c,b=c,b,a; a=(T1+T2)&MASK
    H_out=[(v+iv)&MASK for v,iv in zip([a,b,c,d,e,f,g,h],IV)]
    return raws,H_out

def rank_norm(a):
    n=len(a); o=np.argsort(a); r=np.zeros(n); r[o]=np.arange(n)/n; return r


# ============================================================
# 4.1: STATE COUPLING ATTACK
# ============================================================

def front_41(N=15000, seed=2000):
    """
    Can we predict raw[63] from EARLY carry events through STATE?
    Key: carry[r_early]=0 → state constrained → propagates forward
    through 50+ rounds of state mixing → affects raw[63]?

    If corr(carry_early, raw63 | W[0] fixed) > 0 → state coupling attack.
    """
    print("="*70)
    print("4.1: STATE COUPLING ATTACK — early carry → raw[63] through state")
    print(f"N={N}")
    print("="*70)

    rng = np.random.RandomState(seed)
    all_raws = np.zeros((N,64)); all_carries = np.zeros((N,64),dtype=np.int8)
    t0=time.time()
    for i in range(N):
        W16=[int(rng.randint(0,1<<32)) for _ in range(16)]
        raws,_=sha_full(W16)
        all_raws[i]=raws; all_carries[i]=[1 if r>=T_VAL else 0 for r in raws]
    print(f"  Collected: {time.time()-t0:.1f}s")

    # For each early round: does carry[r]=0 predict raw[63]?
    raw63=all_raws[:,63]
    raw63_rank=rank_norm(raw63)

    print(f"\n  corr(carry[r]=0 flag, rank(raw63)):")
    print(f"  {'r':>3} | {'corr':>8} | {'E[rank63|c=0]':>14} | {'N(c=0)':>7} | {'signal?':>8}")
    print(f"  {'-'*3}-+-{'-'*8}-+-{'-'*14}-+-{'-'*7}-+-{'-'*8}")

    signals = []
    for r in range(64):
        c_r = 1 - all_carries[:,r]  # 1 if carry=0
        if np.std(c_r) < 0.001: continue
        corr_val = np.corrcoef(c_r, raw63_rank)[0,1]
        n_c0 = np.sum(c_r)
        mean_rank_c0 = np.mean(raw63_rank[c_r==1]) if n_c0>0 else 0.5
        if abs(corr_val) > 0.015 or r in [0,9,30,47,62,63]:
            marker = " ★" if abs(corr_val)>0.02 else ""
            print(f"  {r:3d} | {corr_val:+8.4f} | {mean_rank_c0:14.4f} | {n_c0:7d} | {marker:>8}")
            if abs(corr_val) > 0.02:
                signals.append((r, corr_val, n_c0))

    # COMPOSITE SCORE: sum of carry[r]=0 indicators for best rounds
    if signals:
        best_rounds = [r for r,c,n in signals if abs(c)>0.02][:8]
        composite = np.zeros(N)
        for r in best_rounds:
            composite += (1-all_carries[:,r]) * np.sign(
                np.corrcoef(1-all_carries[:,r], raw63_rank)[0,1])

        corr_composite = np.corrcoef(composite, raw63_rank)[0,1]
        print(f"\n  COMPOSITE SCORE (sum of {len(best_rounds)} carry flags):")
        print(f"    Rounds: {best_rounds}")
        print(f"    corr(composite, rank(raw63)): {corr_composite:+.5f}")

        # Stratify: top/bottom quartile of composite → P(carry63=0)
        q75 = np.percentile(composite, 75)
        mask_hi = composite >= q75
        p_c63_hi = np.mean(all_raws[mask_hi,63] < T_VAL)
        p_c63_lo = np.mean(all_raws[~mask_hi,63] < T_VAL)
        print(f"    P(carry63=0 | composite high): {p_c63_hi:.5f}")
        print(f"    P(carry63=0 | composite low):  {p_c63_lo:.5f}")
        print(f"    Lift: {p_c63_hi/max(p_c63_lo, 1e-10):.1f}×")

    return signals


# ============================================================
# 4.2: MULTI-WORD TAIL DISTINGUISHER
# ============================================================

def front_42(N=15000, seed=2001):
    """
    H[5]↔H[6] corr jumps to +0.201 in raw63 tail (Stage 3.4).
    Build a 3-word distinguisher using H[5,6,7] jointly.
    """
    print("\n"+"="*70)
    print("4.2: MULTI-WORD TAIL DISTINGUISHER — H[5,6,7] joint structure")
    print(f"N={N}")
    print("="*70)

    rng=np.random.RandomState(seed)
    all_H=np.zeros((N,8)); all_raw63=np.zeros(N)
    t0=time.time()
    for i in range(N):
        W16=[int(rng.randint(0,1<<32)) for _ in range(16)]
        raws,H_out=sha_full(W16)
        all_H[i]=H_out; all_raw63[i]=raws[63]
    print(f"  Collected: {time.time()-t0:.1f}s")

    raw63_rank = rank_norm(all_raw63)

    # Build features from H[5,6,7] bits
    features = []
    feat_names = []
    for w in [5,6,7]:
        for b in range(28,32):
            feat = np.array([(int(all_H[i,w])>>b)&1 for i in range(N)], dtype=float)
            features.append(feat)
            feat_names.append(f"H[{w}][b{b}]")

    # Add pairwise products
    for i in range(len(features)):
        for j in range(i+1, len(features)):
            if abs(i-j) <= 4:  # nearby bits only
                features.append(features[i] * features[j])
                feat_names.append(f"{feat_names[i]}×{feat_names[j]}")

    X = np.column_stack(features)

    # Correlation with raw63
    print(f"\n  Features correlated with rank(raw63):")
    corrs = []
    for k in range(X.shape[1]):
        c = np.corrcoef(X[:,k], raw63_rank)[0,1] if np.std(X[:,k])>0.01 else 0
        corrs.append((feat_names[k], c))

    corrs.sort(key=lambda x: -abs(x[1]))
    for name, c in corrs[:15]:
        marker = " ★" if abs(c)>0.05 else ""
        print(f"    {name:>30}: {c:+.5f}{marker}")

    # Optimal multi-word predictor (train/test)
    Nt = N//2
    X_tr = X[:Nt] - X[:Nt].mean(0)
    y_tr = raw63_rank[:Nt] - raw63_rank[:Nt].mean()
    X_te = X[Nt:] - X[Nt:].mean(0)
    y_te = raw63_rank[Nt:]

    try:
        beta = np.linalg.lstsq(X_tr, y_tr, rcond=None)[0]
        pred = X_te @ beta
        corr_multi = np.corrcoef(pred, y_te)[0,1]

        # Compare with single H[6][b31]
        h6b31 = np.array([(int(all_H[i,6])>>31)&1 for i in range(Nt,N)], dtype=float)
        corr_single = np.corrcoef(h6b31, y_te)[0,1]

        print(f"\n  MULTI-WORD DISTINGUISHER:")
        print(f"    Single H[6][b31]:            corr = {corr_single:+.5f}")
        print(f"    Multi-word ({X.shape[1]} features): corr = {corr_multi:+.5f}")
        print(f"    Amplification:               {abs(corr_multi)/max(abs(corr_single),0.001):.2f}×")
    except:
        corr_multi = 0
        print(f"  Regression failed")

    return corr_multi


# ============================================================
# 4.3: CONTINUOUS HENSEL — raw mod 2^k coupling
# ============================================================

def front_43(N=15000, seed=2002):
    """
    Classical Hensel: lift solution mod 2^k → mod 2^{k+1}.
    Continuous version: does raw[63] mod 2^k predict raw[63] mod 2^{k+1}?
    And: does raw[62] mod 2^k help predict raw[63] mod 2^k?
    """
    print("\n"+"="*70)
    print("4.3: CONTINUOUS HENSEL — raw mod 2^k coupling")
    print(f"N={N}")
    print("="*70)

    rng=np.random.RandomState(seed)
    raw62=np.zeros(N); raw63=np.zeros(N)
    t0=time.time()
    for i in range(N):
        W16=[int(rng.randint(0,1<<32)) for _ in range(16)]
        raws,_=sha_full(W16)
        raw62[i]=raws[62]; raw63[i]=raws[63]
    print(f"  Collected: {time.time()-t0:.1f}s")

    # Test 1: raw[63] mod 2^k vs raw[63] mod 2^{k+1}
    print(f"\n  --- Self-lifting: raw[63] mod 2^k → raw[63] mod 2^{{k+1}} ---")
    print(f"  {'k':>3} | {'corr(mod2^k, mod2^{k+1})':>25} | {'MI bits':>8}")
    print(f"  {'-'*3}-+-{'-'*25}-+-{'-'*8}")

    for k in range(1, 20):
        mod_k = (raw63.astype(np.int64)) % (1 << k)
        mod_k1 = (raw63.astype(np.int64)) % (1 << (k+1))
        c = np.corrcoef(mod_k.astype(float), mod_k1.astype(float))[0,1]
        # The top bit of mod_{k+1} is the "new" bit
        new_bit = ((mod_k1.astype(np.int64)) >> k) & 1
        old_val = mod_k.astype(float)
        c_new = np.corrcoef(old_val, new_bit.astype(float))[0,1] if np.std(new_bit)>0.01 else 0
        print(f"  {k:3d} | {c:+25.5f} | {c_new:+8.5f}")

    # Test 2: raw[62] mod 2^k helps predict raw[63] mod 2^k?
    print(f"\n  --- Cross-round: raw[62] mod 2^k → raw[63] mod 2^k ---")
    for k in [4, 8, 12, 16, 20, 24, 28, 32]:
        mod62 = (raw62.astype(np.int64)) % (1 << k)
        mod63 = (raw63.astype(np.int64)) % (1 << k)
        c = np.corrcoef(mod62.astype(float), mod63.astype(float))[0,1]
        print(f"    k={k:2d}: corr(raw62 mod 2^k, raw63 mod 2^k) = {c:+.5f}")

    # Test 3: can we predict bit k of raw[63] from bits 0..k-1?
    print(f"\n  --- Bit prediction: bits 0..k-1 of raw[63] → bit k ---")
    raw63_int = raw63.astype(np.int64)
    for k in [1,2,4,8,16,24,31]:
        low_bits = (raw63_int % (1 << k)).astype(float)
        target_bit = ((raw63_int >> k) & 1).astype(float)
        if np.std(target_bit) > 0.01 and np.std(low_bits) > 0:
            c = np.corrcoef(low_bits, target_bit)[0,1]
        else:
            c = 0
        print(f"    bit {k:2d} from bits 0..{k-1}: corr = {c:+.5f}")

    # KEY: is there 2-adic structure in raw?
    # If bits are independent → no Hensel possible
    # If corr(low_bits, high_bit) > 0 → 2-adic lifting possible
    print(f"\n  ★ 2-ADIC STRUCTURE:")
    total_corr = 0
    for k in range(1, 32):
        low = (raw63_int % (1 << k)).astype(float)
        high = ((raw63_int >> k) & 1).astype(float)
        if np.std(high)>0.01 and np.std(low)>0:
            c = abs(np.corrcoef(low, high)[0,1])
            total_corr += c
    avg_corr = total_corr / 31
    print(f"    Mean |corr(bits 0..k-1, bit k)|: {avg_corr:.5f}")
    print(f"    Expected if independent: ~{1/np.sqrt(N):.5f}")
    if avg_corr > 3/np.sqrt(N):
        print(f"    → 2-ADIC STRUCTURE DETECTED (avg corr {avg_corr/max(1/np.sqrt(N),1e-10):.1f}× noise)")
    else:
        print(f"    → No 2-adic structure (bits independent)")

    return avg_corr


# ============================================================
# 4.4: MULTI-BLOCK — carry structure across blocks
# ============================================================

def front_44(N=8000, seed=2003):
    """
    Block 1: compute SHA → H1 (= IV2).
    Block 2: compute SHA with IV=H1.
    Question: does carry structure of block 1 affect block 2 output?
    """
    print("\n"+"="*70)
    print("4.4: MULTI-BLOCK — carry propagation through IV")
    print(f"N={N}")
    print("="*70)

    rng=np.random.RandomState(seed)

    # Block 1: random message
    # Block 2: random message, IV = H1
    raw63_b1=[]; raw63_b2=[]; h7_b2=[]; h6_b2=[]
    h7_b1=[]; h6_b1=[]

    t0=time.time()
    for i in range(N):
        # Block 1
        W16_b1=[int(rng.randint(0,1<<32)) for _ in range(16)]
        raws1, H1 = sha_full(W16_b1)
        raw63_b1.append(raws1[63])
        h7_b1.append(H1[7]); h6_b1.append(H1[6])

        # Block 2 with IV = H1
        W16_b2=[int(rng.randint(0,1<<32)) for _ in range(16)]
        raws2, H2 = sha_compress(H1, W16_b2)
        raw63_b2.append(raws2[63])
        h7_b2.append(H2[7]); h6_b2.append(H2[6])

    print(f"  Collected: {time.time()-t0:.1f}s")

    r1 = np.array(raw63_b1); r2 = np.array(raw63_b2)
    h7_1 = np.array(h7_b1,dtype=float); h7_2 = np.array(h7_b2,dtype=float)

    # Does raw63 of block 1 affect raw63 of block 2?
    corr_r1_r2 = np.corrcoef(r1, r2)[0,1]
    corr_r1_h72 = np.corrcoef(r1, h7_2)[0,1]

    print(f"\n  --- Block 1 → Block 2 propagation ---")
    print(f"  corr(raw63_b1, raw63_b2): {corr_r1_r2:+.5f}")
    print(f"  corr(raw63_b1, H7_b2):    {corr_r1_h72:+.5f}")

    # H7 of block 1 → H7 of block 2?
    corr_h71_h72 = np.corrcoef(h7_1, h7_2)[0,1]
    print(f"  corr(H7_b1, H7_b2):       {corr_h71_h72:+.5f}")

    # Tail analysis: when raw63_b1 is small, is raw63_b2 affected?
    q5 = np.percentile(r1, 5)
    mask = r1 < q5
    n_tail = np.sum(mask)

    if n_tail > 20:
        mean_r2_tail = np.mean(r2[mask])
        mean_r2_base = np.mean(r2)
        print(f"\n  Tail (raw63_b1 < 5%, N={n_tail}):")
        print(f"    E[raw63_b2 | tail]:    {mean_r2_tail/T_VAL:.4f}×T")
        print(f"    E[raw63_b2 | baseline]: {mean_r2_base/T_VAL:.4f}×T")
        print(f"    Delta: {(mean_r2_tail-mean_r2_base)/T_VAL:+.4f}×T")

        # H7 bit biases in block 2 conditioned on block 1 tail
        for bit in [29,30,31]:
            p_tail = np.mean([(int(h7_b2[i])>>bit)&1 for i in range(N) if mask[i]])
            p_base = np.mean([(int(h7_b2[i])>>bit)&1 for i in range(N)])
            delta = p_tail - p_base
            marker = " ★" if abs(delta)>0.02 else ""
            print(f"    H7_b2[b{bit}]: tail={p_tail:.4f}, base={p_base:.4f}, Δ={delta:+.4f}{marker}")

    # KEY: does structure propagate?
    # Mechanism: small raw63_b1 → constrained state → biased H1 → biased IV2
    # → biased early state in block 2 → ...but 64 rounds of mixing → ?
    print(f"\n  ★ MULTI-BLOCK PROPAGATION:")
    if abs(corr_r1_r2) > 0.02:
        print(f"    SIGNAL: raw63_b1 → raw63_b2 corr={corr_r1_r2:+.4f}")
    elif abs(corr_r1_h72) > 0.02:
        print(f"    WEAK: raw63_b1 → H7_b2 corr={corr_r1_h72:+.4f}")
    else:
        print(f"    NO PROPAGATION: block 2 independent of block 1 carry structure")
        print(f"    64 rounds of mixing destroy any IV bias from block 1")

    return corr_r1_r2, corr_r1_h72


# ============================================================
# MAIN
# ============================================================

if __name__ == '__main__':
    print("Four Fronts — Stage 4: All Open Directions")
    print(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print()

    N = 12000
    if len(sys.argv) > 1 and sys.argv[1] == '--fast':
        N = 6000

    t_start = time.time()

    signals_41 = front_41(N=N)
    corr_42 = front_42(N=N)
    corr_43 = front_43(N=N)
    corr_44_r, corr_44_h = front_44(N=N)

    total = time.time()-t_start
    print(f"\n{'='*70}")
    print(f"TOTAL TIME: {total:.1f}s")
    print(f"{'='*70}")

    print(f"""
{'='*70}
SYNTHESIS: Four Fronts (Stage 4)
{'='*70}

4.1 STATE COUPLING: early carry → raw63 through state?
    Signals found: {len(signals_41) if signals_41 else 0}

4.2 MULTI-WORD: H[5,6,7] joint → raw63 predictor
    Multi-word corr: {corr_42:+.4f}

4.3 CONTINUOUS HENSEL: 2-adic structure in raw[63]?
    Mean bit correlation: {corr_43:.5f}

4.4 MULTI-BLOCK: block 1 carry → block 2 output?
    corr(raw63_b1, raw63_b2): {corr_44_r:+.4f}
    corr(raw63_b1, H7_b2):   {corr_44_h:+.4f}

WHICH FRONTS SHOW PROMISE?
""")
