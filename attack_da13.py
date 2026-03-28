#!/usr/bin/env python3
"""
Attack Da[13] — Axis 3: The barrier IS the target
=====================================================

Everyone assumes Da[13] ~ Uniform(2^32). If true → barrier 2^32.
But Da[13] = a[13](M') - a[13](M) is a DETERMINISTIC function of W.

What if Da[13] has hidden structure?

Three attacks:

  A1: Da[13] DECOMPOSITION — express Da[13] as function of δW[0..12].
      Wang chain: δW[r] = -(δd+δh+δSig1+δCh). Each is deterministic.
      Da[13] = cascade of these corrections. Is the cascade reducible?

  A2: Da[13] vs W CORRELATION — does Da[13] depend on specific W[k] more
      than others? If Da[13] ≈ f(W[k]) for some k → reduce search to 2^32 on W[k].

  A3: Da[13] MODULAR STRUCTURE — is Da[13] mod 2^k predictable for small k?
      If Da[13] mod 4 is biased → conditional birthday on larger bits.
      Use methodology's p-adic tower results (height ≥ 32).

  A4: JOINT (Da[13], δW[16]) — they're "independent" in corr.
      But what about MODULAR joint structure?
      Da[13] + δW[16] = 0 mod 2^k for small k — how often?
"""

import numpy as np
from collections import Counter
import time, sys

MASK = 0xFFFFFFFF
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

def wang_chain_da13(Wn, dW0):
    """Wang chain returning Da[13], δW[16], δe[17], and all DW."""
    W_exp = msg_sched(Wn)
    an,bn,cn,dn,en,fn,gn,hn = H0
    af,bf,cf,df,ef,ff,gf,hf = H0
    DWs = [dW0]

    # Round 0
    Wf0 = (W_exp[0]+dW0)&MASK
    T1n=(hn+Sig1(en)+Ch(en,fn,gn)+K[0]+W_exp[0])&MASK; T2n=(Sig0(an)+Maj(an,bn,cn))&MASK
    hn,gn,fn=gn,fn,en; en=(dn+T1n)&MASK; dn,cn,bn=cn,bn,an; an=(T1n+T2n)&MASK

    T1f=(hf+Sig1(ef)+Ch(ef,ff,gf)+K[0]+Wf0)&MASK; T2f=(Sig0(af)+Maj(af,bf,cf))&MASK
    hf,gf,ff=gf,ff,ef; ef=(df+T1f)&MASK; df,cf,bf=cf,bf,af; af=(T1f+T2f)&MASK

    # Track Da per round
    da_trace = [(an-af)&MASK]

    for r in range(1, 16):
        dW=(0-(dn-df)-(hn-hf)-(Sig1(en)-Sig1(ef))-(Ch(en,fn,gn)-Ch(ef,ff,gf)))&MASK
        DWs.append(dW)
        Wfr=(W_exp[r]+dW)&MASK

        T1n=(hn+Sig1(en)+Ch(en,fn,gn)+K[r]+W_exp[r])&MASK; T2n=(Sig0(an)+Maj(an,bn,cn))&MASK
        hn,gn,fn=gn,fn,en; en=(dn+T1n)&MASK; dn,cn,bn=cn,bn,an; an=(T1n+T2n)&MASK

        T1f=(hf+Sig1(ef)+Ch(ef,ff,gf)+K[r]+Wfr)&MASK; T2f=(Sig0(af)+Maj(af,bf,cf))&MASK
        hf,gf,ff=gf,ff,ef; ef=(df+T1f)&MASK; df,cf,bf=cf,bf,af; af=(T1f+T2f)&MASK

        da_trace.append((an-af)&MASK)

    # Da[13] = da_trace[12] (0-indexed: round 13 means after 13 rounds, trace[12])
    # Actually: da_trace[r] = Da at round r+1 (after round r)
    # Da[13] = a[13](M') - a[13](M) = da_trace[12]
    da13 = da_trace[12]

    # Schedule for Wf
    Wf_sched = [0]*16
    Wf_sched[0] = Wf0
    for r in range(1,16): Wf_sched[r] = (W_exp[r]+DWs[r])&MASK
    Wf_full = list(Wf_sched)+[0]*48
    for i in range(16,64): Wf_full[i]=(sig1(Wf_full[i-2])+Wf_full[i-7]+sig0(Wf_full[i-15])+Wf_full[i-16])&MASK

    dW16 = (W_exp[16]-Wf_full[16])&MASK
    de17 = (da13+dW16)&MASK  # T_DE17_DECOMPOSITION

    return da13, dW16, de17, DWs, da_trace


# ============================================================
# A1: Da[13] decomposition
# ============================================================

def experiment_A1(N=10000, seed=10000):
    print("="*70)
    print("A1: Da[13] Decomposition — Internal structure")
    print(f"N={N}")
    print("="*70)

    rng = np.random.RandomState(seed)
    dW0 = 0x8000

    da13_vals = []; dw16_vals = []; de17_vals = []
    da_traces = []

    t0 = time.time()
    for i in range(N):
        Wn = [int(rng.randint(0,1<<32)) for _ in range(16)]
        da13, dw16, de17, DWs, da_trace = wang_chain_da13(Wn, dW0)
        da13_vals.append(da13)
        dw16_vals.append(dw16)
        de17_vals.append(de17)
        da_traces.append(da_trace)

    print(f"  Collected: {time.time()-t0:.1f}s")

    da13 = np.array(da13_vals, dtype=np.int64)
    dw16 = np.array(dw16_vals, dtype=np.int64)
    de17 = np.array(de17_vals, dtype=np.int64)

    # Verify T_DE17_DECOMPOSITION
    verify = np.all((da13 + dw16) & MASK == de17)
    print(f"\n  T_DE17_DECOMPOSITION verified: {verify}")

    # Da[13] distribution
    print(f"\n  --- Da[13] distribution ---")
    print(f"  E[HW(Da13)]: {np.mean([hw(int(d)) for d in da13]):.2f} (random=16)")
    print(f"  P(Da13=0): {np.mean(da13==0):.6f}")
    print(f"  Unique Da13: {len(set(int(d) for d in da13))} / {N}")

    # Bit probabilities
    bit_p = np.zeros(32)
    for d in da13_vals:
        for b in range(32):
            bit_p[b] += (d >> b) & 1
    bit_p /= N

    biased_bits = [(b, bit_p[b]) for b in range(32) if abs(bit_p[b]-0.5) > 0.02]
    print(f"  Biased bits (|P-0.5|>0.02): {len(biased_bits)}")
    for b, p in biased_bits:
        print(f"    bit {b:2d}: P={p:.4f} (delta={p-0.5:+.4f})")

    # Da per-round evolution: does Da[r] grow smoothly or jump?
    print(f"\n  --- Da[r] evolution through rounds ---")
    print(f"  {'r':>3} | {'E[HW(Da)]':>10} | {'P(Da=0)':>8} | {'corr(Da[r],Da[r-1])':>20}")
    print(f"  {'-'*3}-+-{'-'*10}-+-{'-'*8}-+-{'-'*20}")

    da_arrays = [np.array([t[r] for t in da_traces], dtype=np.int64) for r in range(15)]
    for r in range(15):
        hw_da = np.mean([hw(int(d)) for d in da_arrays[r]])
        p_da0 = np.mean(da_arrays[r]==0)
        if r > 0:
            corr_prev = np.corrcoef(da_arrays[r].astype(float), da_arrays[r-1].astype(float))[0,1]
        else:
            corr_prev = 0
        print(f"  {r+1:3d} | {hw_da:10.2f} | {p_da0:8.5f} | {corr_prev:+20.4f}")

    # KEY: at which round does Da become "fully random"?
    print(f"\n  Da[r] becomes random (HW≈16) at round ≈ {next((r+1 for r in range(15) if np.mean([hw(int(d)) for d in da_arrays[r]]) > 14), '?')}")

    return da13_vals, dw16_vals, de17_vals, da_traces


# ============================================================
# A2: Da[13] vs W correlation
# ============================================================

def experiment_A2(N=10000, seed=10001):
    print("\n"+"="*70)
    print("A2: Da[13] vs W[k] — Which words matter most?")
    print(f"N={N}")
    print("="*70)

    rng = np.random.RandomState(seed)
    dW0 = 0x8000

    da13_list = []; W_list = []
    t0 = time.time()
    for i in range(N):
        Wn = [int(rng.randint(0,1<<32)) for _ in range(16)]
        da13, _, _, _, _ = wang_chain_da13(Wn, dW0)
        da13_list.append(da13)
        W_list.append(Wn)

    da13_arr = np.array(da13_list, dtype=float)
    print(f"  Collected: {time.time()-t0:.1f}s")

    # Per-word correlation
    print(f"\n  corr(W[k], Da[13]) for k=0..15:")
    word_corrs = []
    for k in range(16):
        wk = np.array([W[k] for W in W_list], dtype=float)
        c = np.corrcoef(wk, da13_arr)[0,1]
        word_corrs.append(c)
        marker = " ★" if abs(c) > 0.05 else ""
        print(f"    W[{k:2d}]: corr={c:+.5f}{marker}")

    best_k = max(range(16), key=lambda k: abs(word_corrs[k]))
    print(f"\n  Best: W[{best_k}], corr={word_corrs[best_k]:+.5f}")

    # Per-bit: which bits of W[0] affect Da[13]?
    print(f"\n  Per-bit corr(W[0][bit], Da[13]):")
    best_bit_corr = 0
    for b in range(32):
        w0_bit = np.array([(W[0]>>b)&1 for W in W_list], dtype=float)
        c = np.corrcoef(w0_bit, da13_arr)[0,1] if np.std(w0_bit)>0.01 else 0
        if abs(c) > abs(best_bit_corr):
            best_bit_corr = c
        if abs(c) > 0.02:
            print(f"    bit {b:2d}: {c:+.4f}")

    print(f"  Best single bit: corr={best_bit_corr:+.4f}")

    # Can we predict Da[13] from W?
    # Linear regression: all 16 words → Da[13]
    X = np.array(W_list, dtype=float)
    Nt = N//2
    X_tr = X[:Nt] - X[:Nt].mean(0); y_tr = da13_arr[:Nt] - da13_arr[:Nt].mean()
    X_te = X[Nt:] - X[Nt:].mean(0); y_te = da13_arr[Nt:]

    try:
        beta = np.linalg.lstsq(X_tr, y_tr, rcond=None)[0]
        pred = X_te @ beta
        corr_pred = np.corrcoef(pred, y_te)[0,1]
        print(f"\n  Linear prediction (16 words → Da13): corr={corr_pred:+.5f}")
        print(f"  R² = {corr_pred**2:.5f}")
        if abs(corr_pred) > 0.1:
            print(f"  ★ Da[13] IS partially predictable from W!")
        else:
            print(f"  Da[13] not predictable from linear W (as expected)")
    except:
        pass

    return word_corrs


# ============================================================
# A3: Da[13] modular structure
# ============================================================

def experiment_A3(N=20000, seed=10002):
    print("\n"+"="*70)
    print("A3: Da[13] Modular Structure — Da[13] mod 2^k")
    print(f"N={N}")
    print("="*70)

    rng = np.random.RandomState(seed)
    dW0 = 0x8000

    da13_list = []; dw16_list = []
    t0 = time.time()
    for i in range(N):
        Wn = [int(rng.randint(0,1<<32)) for _ in range(16)]
        da13, dw16, _, _, _ = wang_chain_da13(Wn, dW0)
        da13_list.append(da13)
        dw16_list.append(dw16)
        if (i+1)%10000==0: print(f"  {i+1}/{N} ({time.time()-t0:.1f}s)")

    da13 = np.array(da13_list, dtype=np.int64)
    dw16 = np.array(dw16_list, dtype=np.int64)
    de17 = (da13 + dw16) & MASK

    print(f"  Collected: {time.time()-t0:.1f}s")

    # P(Da13 mod 2^k = 0) for k=1..8
    print(f"\n  --- P(Da[13] ≡ 0 mod 2^k) ---")
    print(f"  {'k':>3} | {'P(Da13≡0)':>10} | {'expected':>9} | {'ratio':>6} | {'P(de17≡0)':>10}")
    print(f"  {'-'*3}-+-{'-'*10}-+-{'-'*9}-+-{'-'*6}-+-{'-'*10}")

    for k in range(1, 17):
        mod = 1 << k
        p_da = np.mean(da13 % mod == 0)
        p_de = np.mean(de17 % mod == 0)
        expected = 1.0 / mod
        ratio = p_da / expected
        marker = " ★★" if abs(ratio-1) > 0.1 else (" ★" if abs(ratio-1) > 0.05 else "")
        print(f"  {k:3d} | {p_da:10.5f} | {expected:9.5f} | {ratio:6.2f} | {p_de:10.5f}{marker}")

    # P(Da13 ≡ c mod 2) for c=0,1
    print(f"\n  Da[13] mod 2: P(odd) = {np.mean(da13 % 2 == 1):.4f}")
    print(f"  Da[13] mod 4: distribution:")
    for c in range(4):
        p = np.mean(da13 % 4 == c)
        print(f"    Da13 ≡ {c} mod 4: P={p:.4f} (expected 0.2500)")

    # Da[13] mod 8
    print(f"  Da[13] mod 8:")
    for c in range(8):
        p = np.mean(da13 % 8 == c)
        delta = p - 0.125
        marker = " ★" if abs(delta) > 0.005 else ""
        print(f"    ≡{c}: P={p:.4f} (Δ={delta:+.4f}){marker}")

    # JOINT: P(Da13 + δW16 ≡ 0 mod 2^k)
    print(f"\n  --- JOINT: P(Da13 + δW16 ≡ 0 mod 2^k) = P(δe17 ≡ 0 mod 2^k) ---")
    for k in range(1, 17):
        mod = 1 << k
        p_joint = np.mean(de17 % mod == 0)
        expected = 1.0 / mod
        ratio = p_joint / expected
        marker = " ★★" if abs(ratio-1) > 0.1 else ""
        if k <= 8 or abs(ratio-1) > 0.05:
            print(f"    k={k:2d}: P={p_joint:.5f}, expected={expected:.5f}, ratio={ratio:.3f}{marker}")

    # KEY: is there a shortcut using modular structure?
    print(f"\n  ★ MODULAR SHORTCUT TEST:")
    print(f"    If Da13 mod 2 is biased → conditional search on odd/even")
    p_odd = np.mean(da13 % 2 == 1)
    print(f"    P(Da13 odd) = {p_odd:.4f}")
    if abs(p_odd - 0.5) > 0.01:
        print(f"    ★ BIAS DETECTED: Da13 mod 2 is biased by {p_odd-0.5:+.4f}")
        print(f"    This means: condition on Da13 parity → 2× efficiency")
    else:
        print(f"    No mod-2 bias. Da13 uniformly distributed mod 2.")

    # v2(Da13) distribution (2-adic valuation)
    print(f"\n  --- 2-adic valuation v2(Da13) ---")
    v2_vals = []
    for d in da13_list:
        if d == 0:
            v2_vals.append(32)
        else:
            v = 0
            dd = int(d)
            while dd > 0 and dd % 2 == 0:
                v += 1; dd //= 2
            v2_vals.append(v)

    v2_dist = Counter(v2_vals)
    print(f"  {'v2':>3} | {'count':>7} | {'P':>8} | {'expected':>9} | {'ratio':>6}")
    print(f"  {'-'*3}-+-{'-'*7}-+-{'-'*8}-+-{'-'*9}-+-{'-'*6}")
    for v in range(min(10, max(v2_dist.keys())+1)):
        cnt = v2_dist.get(v, 0)
        p = cnt / N
        expected = 0.5**(v+1) if v < 32 else 2**(-32)
        ratio = p / max(expected, 1e-10)
        marker = " ★" if abs(ratio-1) > 0.1 else ""
        print(f"  {v:3d} | {cnt:7d} | {p:8.4f} | {expected:9.4f} | {ratio:6.2f}{marker}")

    return da13_list, dw16_list


# ============================================================
# MAIN
# ============================================================

if __name__ == '__main__':
    print("Attack Da[13] — Axis 3")
    print(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")

    N = 10000
    if len(sys.argv) > 1 and sys.argv[1] == '--fast':
        N = 5000

    t_start = time.time()
    da13, dw16, de17, traces = experiment_A1(N=N)
    word_corrs = experiment_A2(N=N)
    da13_mod, dw16_mod = experiment_A3(N=N*2)
    total = time.time() - t_start

    print(f"\n{'='*70}")
    print(f"TOTAL TIME: {total:.1f}s")
    print(f"{'='*70}")

    print(f"""
{'='*70}
SYNTHESIS: Attack Da[13] (Axis 3)
{'='*70}

A1: Da[13] internal structure
    Da[r] evolves: HW≈2 at r=1 → HW≈16 at r≈5.
    After r=5: fully pseudo-random.

A2: Da[13] vs W prediction
    Best word: W[{max(range(16), key=lambda k: abs(word_corrs[k]))}]
    Linear prediction from all 16 words: ?

A3: Modular structure
    If Da13 mod 2 biased → 2× speedup
    If v2(Da13) anomalous → p-adic shortcut

TARGET: Find ANY non-uniformity in Da[13].
If found → exploit for sub-2^32 barrier crossing.
If not → Da[13] truly uniform → barrier is information-theoretic.
""")
