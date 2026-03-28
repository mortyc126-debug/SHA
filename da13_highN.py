#!/usr/bin/env python3
"""
High-N Da[13] Confirmation — Axis 3.2
========================================

Stage 3.1 found two borderline signals at N=10K:
  1. Da[13] mod 2 bias: P(odd)=0.5107, Z=2.14σ
  2. δe[17] mod 2^k excess: ratio 2-6× at k≥12 (1-3 events)

N=100K will give Z>6σ if signal is real, or Z<2 if noise.
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

def msg_sched(W16):
    W=list(W16)+[0]*48
    for i in range(16,64): W[i]=(sig1(W[i-2])+W[i-7]+sig0(W[i-15])+W[i-16])&MASK
    return W

def wang_fast(Wn, dW0):
    """Optimized Wang chain returning only Da13, δW16, δe17."""
    W_exp = msg_sched(Wn)
    an,bn,cn,dn,en,fn,gn,hn = H0
    af,bf,cf,df,ef,ff,gf,hf = H0
    DWs_1 = dW0; DWs_9 = 0; DWs_14 = 0

    # Round 0
    T1n=(hn+Sig1(en)+Ch(en,fn,gn)+K[0]+W_exp[0])&MASK; T2n=(Sig0(an)+Maj(an,bn,cn))&MASK
    hn,gn,fn=gn,fn,en; en=(dn+T1n)&MASK; dn,cn,bn=cn,bn,an; an=(T1n+T2n)&MASK

    T1f=(hf+Sig1(ef)+Ch(ef,ff,gf)+K[0]+((W_exp[0]+dW0)&MASK))&MASK; T2f=(Sig0(af)+Maj(af,bf,cf))&MASK
    hf,gf,ff=gf,ff,ef; ef=(df+T1f)&MASK; df,cf,bf=cf,bf,af; af=(T1f+T2f)&MASK

    for r in range(1, 16):
        dW=(0-(dn-df)-(hn-hf)-(Sig1(en)-Sig1(ef))-(Ch(en,fn,gn)-Ch(ef,ff,gf)))&MASK
        if r == 1: DWs_1 = dW
        elif r == 9: DWs_9 = dW
        elif r == 14: DWs_14 = dW
        Wfr=(W_exp[r]+dW)&MASK

        T1n=(hn+Sig1(en)+Ch(en,fn,gn)+K[r]+W_exp[r])&MASK; T2n=(Sig0(an)+Maj(an,bn,cn))&MASK
        hn,gn,fn=gn,fn,en; en=(dn+T1n)&MASK; dn,cn,bn=cn,bn,an; an=(T1n+T2n)&MASK

        T1f=(hf+Sig1(ef)+Ch(ef,ff,gf)+K[r]+Wfr)&MASK; T2f=(Sig0(af)+Maj(af,bf,cf))&MASK
        hf,gf,ff=gf,ff,ef; ef=(df+T1f)&MASK; df,cf,bf=cf,bf,af; af=(T1f+T2f)&MASK

    da13 = (an - af) & MASK  # Da at round 15 = a after round 14...
    # Actually we need Da[13]. Let me track a_n and a_f at round 13.
    # The loop goes r=1..15, so after r=12 we have state at round 13.
    # But I'm computing all 15 rounds. Da[13] = a after round 12 diff.
    # Let me re-do properly.
    return None, None, None  # placeholder

def wang_chain_da13_fast(Wn, dW0):
    """Proper Wang chain returning Da[13], δW[16], δe[17]."""
    W_exp = msg_sched(Wn)
    an,bn,cn,dn,en,fn,gn,hn = H0
    af,bf,cf,df,ef,ff,gf,hf = H0

    # Round 0
    T1n=(hn+Sig1(en)+Ch(en,fn,gn)+K[0]+W_exp[0])&MASK; T2n=(Sig0(an)+Maj(an,bn,cn))&MASK
    hn,gn,fn=gn,fn,en; en=(dn+T1n)&MASK; dn,cn,bn=cn,bn,an; an=(T1n+T2n)&MASK

    T1f=(hf+Sig1(ef)+Ch(ef,ff,gf)+K[0]+((W_exp[0]+dW0)&MASK))&MASK; T2f=(Sig0(af)+Maj(af,bf,cf))&MASK
    hf,gf,ff=gf,ff,ef; ef=(df+T1f)&MASK; df,cf,bf=cf,bf,af; af=(T1f+T2f)&MASK

    Wf_words = [0]*16
    Wf_words[0] = (W_exp[0]+dW0)&MASK

    for r in range(1, 16):
        dW=(0-(dn-df)-(hn-hf)-(Sig1(en)-Sig1(ef))-(Ch(en,fn,gn)-Ch(ef,ff,gf)))&MASK
        Wfr=(W_exp[r]+dW)&MASK
        Wf_words[r] = Wfr

        T1n=(hn+Sig1(en)+Ch(en,fn,gn)+K[r]+W_exp[r])&MASK; T2n=(Sig0(an)+Maj(an,bn,cn))&MASK
        hn,gn,fn=gn,fn,en; en=(dn+T1n)&MASK; dn,cn,bn=cn,bn,an; an=(T1n+T2n)&MASK

        T1f=(hf+Sig1(ef)+Ch(ef,ff,gf)+K[r]+Wfr)&MASK; T2f=(Sig0(af)+Maj(af,bf,cf))&MASK
        hf,gf,ff=gf,ff,ef; ef=(df+T1f)&MASK; df,cf,bf=cf,bf,af; af=(T1f+T2f)&MASK

        # After r=12 (round 13 state): save Da[13]
        if r == 12:
            da13 = (an - af) & MASK

    # Compute δW[16] from schedule
    Wf_sched = list(Wf_words) + [0]*48
    for i in range(16,64):
        Wf_sched[i]=(sig1(Wf_sched[i-2])+Wf_sched[i-7]+sig0(Wf_sched[i-15])+Wf_sched[i-16])&MASK

    dW16 = (W_exp[16] - Wf_sched[16]) & MASK
    de17 = (da13 + dW16) & MASK

    return da13, dW16, de17


def run_highN(N=50000, seed=11000):
    print("="*70)
    print(f"HIGH-N Da[13] CONFIRMATION — N={N}")
    print("="*70)

    rng = np.random.RandomState(seed)
    dW0 = 0x8000

    da13_list = []
    dw16_list = []
    de17_list = []

    t0 = time.time()
    for i in range(N):
        Wn = [int(rng.randint(0,1<<32)) for _ in range(16)]
        da13, dw16, de17 = wang_chain_da13_fast(Wn, dW0)
        da13_list.append(da13)
        dw16_list.append(dw16)
        de17_list.append(de17)

        if (i+1) % 10000 == 0:
            elapsed = time.time()-t0
            rate = (i+1)/elapsed
            print(f"  {i+1}/{N} ({elapsed:.1f}s, {rate:.0f}/s)")

    elapsed = time.time()-t0
    print(f"\n  Completed: {elapsed:.1f}s ({N/elapsed:.0f}/s)")

    da13 = np.array(da13_list, dtype=np.int64)
    dw16 = np.array(dw16_list, dtype=np.int64)
    de17 = np.array(de17_list, dtype=np.int64)

    # Verify
    verify = np.all((da13 + dw16) & MASK == de17)
    print(f"  T_DE17_DECOMPOSITION: {verify}")

    # ===== TEST 1: Da[13] mod 2 bias =====
    print(f"\n  ===== TEST 1: Da[13] mod 2 bias =====")
    p_odd = np.mean(da13 % 2 == 1)
    z_mod2 = (p_odd - 0.5) / (0.5 / np.sqrt(N))
    print(f"  P(Da13 odd) = {p_odd:.5f}")
    print(f"  Z-score: {z_mod2:.2f}")
    if abs(z_mod2) > 3:
        print(f"  ★★★ CONFIRMED at Z={z_mod2:.1f}σ — Da[13] mod 2 IS biased!")
    elif abs(z_mod2) > 2:
        print(f"  ★★ Borderline (Z={z_mod2:.1f}σ)")
    else:
        print(f"  REFUTED — Da[13] mod 2 unbiased (Z={z_mod2:.1f}σ)")

    # ===== TEST 2: Da[13] mod 4 =====
    print(f"\n  ===== TEST 2: Da[13] mod 4 =====")
    for c in range(4):
        p = np.mean(da13 % 4 == c)
        z = (p - 0.25) / np.sqrt(0.25*0.75/N)
        marker = " ★" if abs(z) > 3 else ""
        print(f"  Da13 ≡ {c} mod 4: P={p:.5f} (Z={z:+.2f}){marker}")

    # ===== TEST 3: Da[13] mod 8 =====
    print(f"\n  ===== TEST 3: Da[13] mod 8 =====")
    for c in range(8):
        p = np.mean(da13 % 8 == c)
        z = (p - 0.125) / np.sqrt(0.125*0.875/N)
        marker = " ★" if abs(z) > 3 else ""
        print(f"  Da13 ≡ {c} mod 8: P={p:.5f} (Z={z:+.2f}){marker}")

    # ===== TEST 4: δe[17] mod 2^k excess =====
    print(f"\n  ===== TEST 4: P(δe17 ≡ 0 mod 2^k) =====")
    print(f"  {'k':>3} | {'count':>7} | {'P':>10} | {'expected':>10} | {'ratio':>6} | {'Z':>7}")
    print(f"  {'-'*3}-+-{'-'*7}-+-{'-'*10}-+-{'-'*10}-+-{'-'*6}-+-{'-'*7}")

    for k in range(1, 21):
        mod = 1 << k
        count = int(np.sum(de17 % mod == 0))
        p = count / N
        expected = 1.0 / mod
        n_expected = N * expected
        ratio = p / max(expected, 1e-15)
        # Z-score for binomial
        se = np.sqrt(expected * (1-expected) / N) if expected < 1 else 0.001
        z = (p - expected) / max(se, 1e-10)
        marker = " ★★★" if abs(z)>3 and count>=3 else (" ★★" if abs(z)>2 and count>=2 else "")
        print(f"  {k:3d} | {count:7d} | {p:10.6f} | {expected:10.6f} | {ratio:6.2f} | {z:+7.2f}{marker}")

    # ===== TEST 5: Da[13] bit probabilities =====
    print(f"\n  ===== TEST 5: Da[13] bit probabilities =====")
    biased = []
    for b in range(32):
        p = np.mean((da13 >> b) & 1)
        z = (p - 0.5) / (0.5 / np.sqrt(N))
        if abs(z) > 3:
            biased.append((b, p, z))
            print(f"  bit {b:2d}: P={p:.5f}, Z={z:+.2f} ★★★")

    if not biased:
        print(f"  No biased bits at Z>3σ (N={N})")
    print(f"  Total biased (Z>3): {len(biased)} / 32")

    # ===== TEST 6: v2(Da[13]) =====
    print(f"\n  ===== TEST 6: v2(Da[13]) — 2-adic valuation =====")
    v2_counts = Counter()
    for d in da13_list:
        if d == 0:
            v2_counts[32] += 1; continue
        v = 0; dd = int(d)
        while dd % 2 == 0: v += 1; dd //= 2
        v2_counts[v] += 1

    for v in range(min(12, max(v2_counts.keys())+1)):
        cnt = v2_counts.get(v, 0)
        p = cnt / N
        expected = 0.5**(v+1)
        z = (p - expected) / np.sqrt(expected*(1-expected)/N)
        marker = " ★★★" if abs(z)>3 else ""
        print(f"  v2={v:2d}: count={cnt:6d}, P={p:.5f}, exp={expected:.5f}, Z={z:+.2f}{marker}")

    # ===== TEST 7: δW[0] variation =====
    print(f"\n  ===== TEST 7: Same tests with δW[0]=1 =====")
    rng2 = np.random.RandomState(seed+500)
    da13_alt = []
    N_alt = min(N, 30000)
    for i in range(N_alt):
        Wn = [int(rng2.randint(0,1<<32)) for _ in range(16)]
        da13, _, _ = wang_chain_da13_fast(Wn, 1)  # δW[0]=1
        da13_alt.append(da13)
    da13_a = np.array(da13_alt, dtype=np.int64)
    p_odd_alt = np.mean(da13_a % 2 == 1)
    z_alt = (p_odd_alt - 0.5) / (0.5 / np.sqrt(N_alt))
    print(f"  δW[0]=1: P(Da13 odd) = {p_odd_alt:.5f}, Z={z_alt:.2f}")

    # δW[0]=0x80000000
    da13_b = []
    for i in range(N_alt):
        Wn = [int(rng2.randint(0,1<<32)) for _ in range(16)]
        da13, _, _ = wang_chain_da13_fast(Wn, 0x80000000)
        da13_b.append(da13)
    da13_ba = np.array(da13_b, dtype=np.int64)
    p_odd_b = np.mean(da13_ba % 2 == 1)
    z_b = (p_odd_b - 0.5) / (0.5 / np.sqrt(N_alt))
    print(f"  δW[0]=0x80000000: P(Da13 odd) = {p_odd_b:.5f}, Z={z_b:.2f}")

    return da13, dw16, de17


if __name__ == '__main__':
    print("High-N Da[13] Confirmation — Axis 3.2")
    print(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")

    N = 100000
    if len(sys.argv) > 1 and sys.argv[1] == '--fast':
        N = 50000

    da13, dw16, de17 = run_highN(N=N)
