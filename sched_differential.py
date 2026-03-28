#!/usr/bin/env python3
"""
Schedule-Differential Interaction — Axis 2.2
===============================================

δW[16] = sig1(δW[14]) + δW[9] + sig0(δW[1]) + δW[0]

In Wang chain: δW[0] is chosen, δW[1..15] are ADAPTED (computed
from state to make δe[r+1]=0). So δW[1], δW[9], δW[14] are
determined by the Wang chain, NOT free parameters.

BUT: Wang chain δW[r] depends on δW[0] and the message W.
Different W gives different δW[1..15] → different δW[16].

Question: can we choose W (message) such that δW[16] is small?
If HW(δW[16]) = 0 → δe[17] depends only on state diff (δa[13]).
If HW(δW[16]) small → δe[17] more predictable → cheaper differential.

Experiments:
  S1: δW[16] distribution under Wang chain for various δW[0]
  S2: Search for W with min HW(δW[16])
  S3: Does min HW(δW[16]) actually help P(δe[17]=0)?
  S4: Schedule algebra: can δW[16]=0 analytically?
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


def wang_chain_full(Wn, dW0):
    """
    Build full Wang chain: Wn is base message, dW0 is additive difference in W[0].
    Returns: DWs[0..15], δe[17], δa[13], δW[16..19], full schedules.
    """
    Wf = list(Wn); Wf[0] = (Wf[0] + dW0) & MASK

    # State traces
    an,bn,cn,dn,en,fn,gn,hn = H0
    af,bf,cf,df,ef,ff,gf,hf = H0

    DWs = [dW0]

    # Round 0
    Wn_exp = msg_sched(Wn)

    raw_n = hn+Sig1(en)+Ch(en,fn,gn)+K[0]+Wn_exp[0]
    T1n=raw_n&MASK; T2n=(Sig0(an)+Maj(an,bn,cn))&MASK
    hn,gn,fn=gn,fn,en; en=(dn+T1n)&MASK; dn,cn,bn=cn,bn,an; an=(T1n+T2n)&MASK

    Wf_exp = list(Wn_exp); Wf_exp[0] = (Wn_exp[0]+dW0)&MASK
    raw_f = hf+Sig1(ef)+Ch(ef,ff,gf)+K[0]+Wf_exp[0]
    T1f=raw_f&MASK; T2f=(Sig0(af)+Maj(af,bf,cf))&MASK
    hf,gf,ff=gf,ff,ef; ef=(df+T1f)&MASK; df,cf,bf=cf,bf,af; af=(T1f+T2f)&MASK

    # Rounds 1..15: adapt Wf[r] for δe[r+1]=0
    for r in range(1, 16):
        delta_d=(dn-df)&MASK; delta_h=(hn-hf)&MASK
        delta_sig1=(Sig1(en)-Sig1(ef))&MASK
        delta_ch=(Ch(en,fn,gn)-Ch(ef,ff,gf))&MASK
        delta_W=(0-delta_d-delta_h-delta_sig1-delta_ch)&MASK
        Wf_exp[r]=(Wn_exp[r]+delta_W)&MASK
        DWs.append(delta_W)

        raw_n=hn+Sig1(en)+Ch(en,fn,gn)+K[r]+Wn_exp[r]
        T1n=raw_n&MASK; T2n=(Sig0(an)+Maj(an,bn,cn))&MASK
        hn,gn,fn=gn,fn,en; en=(dn+T1n)&MASK; dn,cn,bn=cn,bn,an; an=(T1n+T2n)&MASK

        raw_f=hf+Sig1(ef)+Ch(ef,ff,gf)+K[r]+Wf_exp[r]
        T1f=raw_f&MASK; T2f=(Sig0(af)+Maj(af,bf,cf))&MASK
        hf,gf,ff=gf,ff,ef; ef=(df+T1f)&MASK; df,cf,bf=cf,bf,af; af=(T1f+T2f)&MASK

    # Expand Wf schedule from Wf_exp[0..15]
    Wf_sched = list(Wf_exp[:16]) + [0]*48
    for i in range(16,64):
        Wf_sched[i]=(sig1(Wf_sched[i-2])+Wf_sched[i-7]+sig0(Wf_sched[i-15])+Wf_sched[i-16])&MASK

    # δW[16..19]
    DW16 = (Wn_exp[16] - Wf_sched[16]) & MASK
    DW17 = (Wn_exp[17] - Wf_sched[17]) & MASK
    DW18 = (Wn_exp[18] - Wf_sched[18]) & MASK

    # Continue rounds 16,17 for δe[17]
    for r in range(16, 18):
        raw_n=hn+Sig1(en)+Ch(en,fn,gn)+K[r]+Wn_exp[r]
        T1n=raw_n&MASK; T2n=(Sig0(an)+Maj(an,bn,cn))&MASK
        hn,gn,fn=gn,fn,en; en=(dn+T1n)&MASK; dn,cn,bn=cn,bn,an; an=(T1n+T2n)&MASK

        raw_f=hf+Sig1(ef)+Ch(ef,ff,gf)+K[r]+Wf_sched[r]
        T1f=raw_f&MASK; T2f=(Sig0(af)+Maj(af,bf,cf))&MASK
        hf,gf,ff=gf,ff,ef; ef=(df+T1f)&MASK; df,cf,bf=cf,bf,af; af=(T1f+T2f)&MASK

    de17 = (en - ef) & MASK
    # δa[13] from methodology: Da13 = a[13](M1) - a[13](M2)
    # We'd need state[13] but simplified here

    return DWs, DW16, DW17, de17


# ============================================================
# S1: δW[16] distribution
# ============================================================

def experiment_S1(N=10000, seed=9000):
    print("="*70)
    print("S1: δW[16] distribution under Wang chain")
    print(f"N={N}")
    print("="*70)

    rng = np.random.RandomState(seed)

    # Test δW[0] = 0x8000 (standard from methodology)
    dW0 = 0x8000

    hw_dw16 = []; hw_dw17 = []; hw_de17 = []
    dw16_vals = []

    t0 = time.time()
    for i in range(N):
        Wn = [int(rng.randint(0,1<<32)) for _ in range(16)]
        DWs, DW16, DW17, de17 = wang_chain_full(Wn, dW0)
        hw_dw16.append(hw(DW16))
        hw_dw17.append(hw(DW17))
        hw_de17.append(hw(de17))
        dw16_vals.append(DW16)

        if (i+1) % 5000 == 0:
            print(f"  {i+1}/{N} ({time.time()-t0:.1f}s)")

    print(f"  Collected: {time.time()-t0:.1f}s")

    print(f"\n  δW[0] = 0x{dW0:08x}")
    print(f"  E[HW(δW[16])]: {np.mean(hw_dw16):.2f} (random=16)")
    print(f"  min HW(δW[16]): {min(hw_dw16)}")
    print(f"  max HW(δW[16]): {max(hw_dw16)}")
    print(f"  E[HW(δW[17])]: {np.mean(hw_dw17):.2f}")
    print(f"  E[HW(δe[17])]: {np.mean(hw_de17):.2f}")
    print(f"  P(δe[17]=0): {np.mean(np.array(hw_de17)==0):.6f}")

    # HW distribution
    from collections import Counter
    hw_dist = Counter(hw_dw16)
    print(f"\n  HW(δW[16]) distribution:")
    for h in sorted(hw_dist.keys()):
        if hw_dist[h] > N*0.005 or h <= 5:
            print(f"    HW={h:2d}: {hw_dist[h]:5d} ({hw_dist[h]/N*100:5.2f}%)")

    # Correlation: HW(δW[16]) vs HW(δe[17])
    corr = np.corrcoef(hw_dw16, hw_de17)[0,1]
    print(f"\n  corr(HW(δW[16]), HW(δe[17])): {corr:+.4f}")

    # KEY: does low HW(δW[16]) predict low HW(δe[17])?
    print(f"\n  --- Stratified: E[HW(δe17)] by HW(δW16) ---")
    for hw_thresh in [3, 5, 8, 10, 14, 16, 20]:
        mask = np.array(hw_dw16) <= hw_thresh
        if np.sum(mask) >= 5:
            mean_de17 = np.mean(np.array(hw_de17)[mask])
            p_de17_0 = np.mean(np.array(hw_de17)[mask] == 0)
            print(f"    HW(δW16)≤{hw_thresh:2d}: N={np.sum(mask):5d}, E[HW(δe17)]={mean_de17:.2f}, P(δe17=0)={p_de17_0:.5f}")

    # Test other δW[0]
    print(f"\n  --- δW[16] for different δW[0] ---")
    for dw0_test in [1, 0x8000, 0x80000000, 0x100, 2]:
        hw_list = []
        for _ in range(min(N, 3000)):
            Wn = [int(rng.randint(0,1<<32)) for _ in range(16)]
            _, DW16, _, _ = wang_chain_full(Wn, dw0_test)
            hw_list.append(hw(DW16))
        print(f"    δW[0]=0x{dw0_test:08x}: E[HW(δW16)]={np.mean(hw_list):.2f}, min={min(hw_list)}")

    return hw_dw16, hw_de17, dw16_vals


# ============================================================
# S2: Search for W with min HW(δW[16])
# ============================================================

def experiment_S2(N=20000, seed=9001):
    print("\n"+"="*70)
    print("S2: Search for message W with minimal HW(δW[16])")
    print(f"N={N}")
    print("="*70)

    rng = np.random.RandomState(seed)
    dW0 = 0x8000  # standard delta

    best_hw = 32; best_Wn = None; best_DW16 = None; best_de17 = None
    hw_history = []

    t0 = time.time()
    for i in range(N):
        Wn = [int(rng.randint(0,1<<32)) for _ in range(16)]
        _, DW16, _, de17 = wang_chain_full(Wn, dW0)
        h = hw(DW16)
        hw_history.append(h)

        if h < best_hw:
            best_hw = h; best_Wn = Wn; best_DW16 = DW16; best_de17 = de17
            print(f"  NEW RECORD: HW(δW16)={h}, de17=0x{de17:08x} HW={hw(de17)}, at i={i}")

    print(f"\n  Search: {time.time()-t0:.1f}s")
    print(f"  Best HW(δW[16]): {best_hw}")
    print(f"  δW[16] = 0x{best_DW16:08x}")
    print(f"  δe[17] = 0x{best_de17:08x} (HW={hw(best_de17)})")

    # Does low HW(δW16) help δe17?
    # Compare: all samples with HW(δW16)<=best_hw+2
    close_mask = np.array(hw_history) <= best_hw + 2
    n_close = np.sum(close_mask)

    if n_close > 0:
        # Recompute δe17 for these
        close_de17 = []
        rng2 = np.random.RandomState(seed)
        for i in range(N):
            Wn = [int(rng2.randint(0,1<<32)) for _ in range(16)]
            _, _, _, de17 = wang_chain_full(Wn, dW0)
            if hw_history[i] <= best_hw + 2:
                close_de17.append(hw(de17))

        print(f"\n  Samples with HW(δW16)≤{best_hw+2}: N={n_close}")
        print(f"    E[HW(δe17)]: {np.mean(close_de17):.2f}")
        print(f"    P(δe17=0): {np.mean(np.array(close_de17)==0):.5f}")

    return best_hw, best_DW16


# ============================================================
# S4: Schedule algebra — can δW[16]=0?
# ============================================================

def experiment_S4():
    print("\n"+"="*70)
    print("S4: Schedule Algebra — Can δW[16]=0?")
    print("="*70)

    # δW[16] = sig1(δW[14]) + δW[9] + sig0(δW[1]) + δW[0]
    #
    # In Wang chain with δW[0]=d:
    #   δW[1] is determined by Wang adaptation (cancel δe[2])
    #   δW[9] is determined by Wang adaptation (cancel δe[10])
    #   δW[14] is determined by Wang adaptation (cancel δe[15])
    #
    # δW[16]=0 requires: sig1(δW[14]) + δW[9] + sig0(δW[1]) + d = 0 (mod 2^32)
    # → sig1(δW[14]) + δW[9] + sig0(δW[1]) = -d (mod 2^32)
    #
    # Each δW[r] is a function of W (message) and d (initial diff).
    # So this is ONE equation in 16 unknowns (W[0..15]).
    # Should have solutions!

    print(f"""
  ANALYSIS:

  δW[16] = sig1(δW[14]) + δW[9] + sig0(δW[1]) + δW[0]

  Wang chain determines δW[1], δW[9], δW[14] as functions of W and δW[0].
  Each δW[r] = -(δd[r] + δh[r] + δSig1(e[r]) + δCh(e[r],f[r],g[r]))

  These depend on the FULL STATE at each round, which depends on W.

  For δW[16]=0: need sig1(δW[14]) + δW[9] + sig0(δW[1]) = -δW[0]

  This is 1 equation in ~16×32 = 512 bits of message freedom.
  VASTLY underdetermined → solutions MUST exist.

  BUT: finding them requires solving a highly nonlinear system:
  - δW[1] depends on e[1], f[1], g[1], d[1], h[1] (round 1 state)
  - δW[9] depends on round 9 state (which depends on W[0..8])
  - δW[14] depends on round 14 state (which depends on W[0..13])

  The system δW[16]=0 is equivalent to:
    F(W) = sig1(G14(W)) + G9(W) + sig0(G1(W)) + δW[0] = 0

  where Gr(W) = Wang-correction at round r = nonlinear function of W.

  This is ONE equation, 512-bit input → should be O(1) to satisfy.
  But the nonlinearity makes it O(2^k) for some k.

  From methodology T_BIRTHDAY_COST17:
    Finding δe[17]=0 costs O(2^32) birthday.
    δe[17]=0 ⟺ Da[13] + δW[16] = 0 (mod 2^32).
    If δW[16]=0: need Da[13]=0 additionally.
    If δW[16] small: birthday on fewer bits.

  KEY QUESTION: Does HW(δW[16])=k reduce birthday from 2^32 to 2^(32-g(k))?
  If yes → optimizing δW[16] directly speeds up the Wang-barrier crossing.
""")

    # Theoretical: δe[17] = Da[13] + δW[16]
    # P(δe[17]=0) = P(Da[13] = -δW[16])
    # If Da[13] ~ Uniform(2^32) and δW[16] fixed:
    #   P = 1/2^32 regardless of δW[16]!
    # UNLESS Da[13] depends on δW[16] (which it does through W)

    print(f"  THEORETICAL:")
    print(f"    δe[17] = Da[13] + δW[16]  (T_DE17_DECOMPOSITION)")
    print(f"    If Da[13] ~ Uniform(2^32): P(δe17=0) = 1/2^32 for ANY δW[16]")
    print(f"    δW[16] value doesn't matter — only Da[13] randomness matters!")
    print(f"")
    print(f"    BUT: if Da[13] and δW[16] are CORRELATED (both depend on W):")
    print(f"    then small δW[16] might correlate with small Da[13]")
    print(f"    → P(δe17=0 | HW(δW16) small) > P(δe17=0 | HW(δW16) random)")
    print(f"")
    print(f"    This is the KEY TEST for Axis 2.2.")


# ============================================================
# MAIN
# ============================================================

if __name__ == '__main__':
    print("Schedule-Differential Interaction — Axis 2.2")
    print(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")

    N1, N2 = 8000, 15000
    if len(sys.argv) > 1 and sys.argv[1] == '--fast':
        N1, N2 = 4000, 8000

    t_start = time.time()
    hw_dw16, hw_de17, dw16_vals = experiment_S1(N=N1)
    best_hw, best_DW16 = experiment_S2(N=N2)
    experiment_S4()
    total = time.time() - t_start

    print(f"\n{'='*70}")
    print(f"TOTAL TIME: {total:.1f}s")
    print(f"{'='*70}")

    corr_key = np.corrcoef(hw_dw16, hw_de17)[0,1]
    print(f"""
{'='*70}
SYNTHESIS: Schedule-Differential (Axis 2.2)
{'='*70}

δW[16] under Wang chain:
  E[HW] = {np.mean(hw_dw16):.1f}, min = {min(hw_dw16)}, record search = {best_hw}

corr(HW(δW16), HW(δe17)) = {corr_key:+.4f}

If corr > 0.1 → δW[16] optimization helps differential trail.
If corr ≈ 0  → δW[16] and δe[17] independent → optimization useless.

THEORETICAL: δe17 = Da13 + δW16.
If Da13 ⊥ δW16 → P(δe17=0) = 1/2^32 regardless of δW16.
If Da13 ∥ δW16 → small δW16 → small δe17 more likely.
""")
