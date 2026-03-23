#!/usr/bin/env python3
"""
CRAZY-23: Extract linear component of barrier W[14] → De17

CRAZY-19 proved: W[14] has 30-85% bit-to-bit correlation with De17.
This means De17 ≈ L(W[14]) + nonlinear_remainder.

Plan:
1. Build full 32×32 correlation matrix W[14] bits → De17 bits
2. Extract best GF(2)-linear approximation L
3. Compute L(W[14]) for each sample, measure residual R = De17 ⊕ L(W[14])
4. If HW(R) << 16 → linear extraction saves bits!
5. Strategy: solve L(W[14]) = target via GF(2) inverse, then
   brute-force the low-HW residual

This could reduce single-barrier cost from 2^32 to 2^HW(R).
"""

import random
import time

MASK = 0xFFFFFFFF
K = [0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc]
H0 = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&MASK
def Sig0(x): return rotr(x,2)^rotr(x,13)^rotr(x,22)
def Sig1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def sig0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def sig1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def Ch(e,f,g): return ((e&f)^(~e&g))&MASK
def Maj(a,b,c): return ((a&b)^(a&c)^(b&c))&MASK
def add32(*a):
    s=0
    for x in a: s=(s+x)&MASK
    return s
def hw(x): return bin(x&MASK).count('1')

def sha_round(st,w,k):
    a,b,c,d,e,f,g,h=st
    T1=add32(h,Sig1(e),Ch(e,f,g),k,w); T2=add32(Sig0(a),Maj(a,b,c))
    return [add32(T1,T2),a,b,c,add32(d,T1),e,f,g]

def expand_W(W16):
    W=list(W16)
    for i in range(16,20): W.append(add32(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    return W

def wang_chain(msg,iv):
    W=list(msg);Wp=list(W);Wp[0]^=0x80000000
    s=list(iv);sp=list(iv)
    s=sha_round(s,W[0],K[0]);sp=sha_round(sp,Wp[0],K[0])
    for t in range(1,16):
        a,b,c,d,e,f,g,h=s;a2,b2,c2,d2,e2,f2,g2,h2=sp
        tp=add32(h,Sig1(e),Ch(e,f,g),K[t]);tp2=add32(h2,Sig1(e2),Ch(e2,f2,g2),K[t])
        target=add32(d,tp,W[t]);Wp[t]=(target-d2-tp2)&MASK
        s=sha_round(s,W[t],K[t]);sp=sha_round(sp,Wp[t],K[t])
    return Wp

def get_de17(msg,iv):
    Wp=wang_chain(msg,iv);We=expand_W(msg);Wpe=expand_W(Wp)
    s1=list(iv);s2=list(iv)
    for t in range(17): s1=sha_round(s1,We[t],K[t]);s2=sha_round(s2,Wpe[t],K[t])
    return s1[4]^s2[4]


def main():
    random.seed(0xC023)
    t0 = time.time()
    iv = list(H0)

    NUM_MSGS = 10
    N_SAMPLES = 5000

    all_residual_hws = []
    all_ranks = []

    for msg_idx in range(NUM_MSGS):
        base_msg = [random.getrandbits(32) for _ in range(16)]

        # ── Step 1: Build GF(2) Jacobian matrix W[14] → De17 ─────
        # Column j = De17(W14 ^ e_j) XOR De17(W14)
        base_de17 = get_de17(base_msg, iv)
        jacobian_cols = []  # 32 columns, each a 32-bit integer
        for j in range(32):
            msg_f = list(base_msg)
            msg_f[14] ^= (1 << j)
            de17_f = get_de17(msg_f, iv)
            jacobian_cols.append(base_de17 ^ de17_f)

        # ── Step 2: Compute rank of this 32×32 GF(2) matrix ──────
        rows = []
        for r in range(32):
            row = 0
            for c in range(32):
                if (jacobian_cols[c] >> r) & 1:
                    row |= (1 << c)
            rows.append(row)

        # Gaussian elimination
        rows_copy = list(rows)
        rank = 0
        pivot_rows = []
        for col in range(32):
            pivot = -1
            for r in range(rank, 32):
                if (rows_copy[r] >> col) & 1:
                    pivot = r; break
            if pivot == -1: continue
            rows_copy[rank], rows_copy[pivot] = rows_copy[pivot], rows_copy[rank]
            for r in range(32):
                if r != rank and (rows_copy[r] >> col) & 1:
                    rows_copy[r] ^= rows_copy[rank]
            pivot_rows.append(rank)
            rank += 1

        all_ranks.append(rank)

        # ── Step 3: Apply linear approximation to many samples ────
        # L(x) ≈ J * x (GF(2) matrix-vector multiply)
        # Residual R = De17(actual) XOR L(delta_W14)
        residual_hws = []
        for _ in range(N_SAMPLES):
            w14_new = random.getrandbits(32)
            delta_w14 = base_msg[14] ^ w14_new

            # Linear prediction: L(delta_w14) = XOR of jacobian_cols[j] where bit j is set
            lin_pred = 0
            for j in range(32):
                if (delta_w14 >> j) & 1:
                    lin_pred ^= jacobian_cols[j]

            # Actual De17
            msg_new = list(base_msg)
            msg_new[14] = w14_new
            actual_de17 = get_de17(msg_new, iv)

            # Residual
            residual = (actual_de17 ^ base_de17) ^ lin_pred
            residual_hws.append(hw(residual))

        avg_res = sum(residual_hws) / len(residual_hws)
        min_res = min(residual_hws)
        all_residual_hws.append(avg_res)

        if msg_idx < 3 or avg_res < 14:
            print(f"  Msg {msg_idx}: rank={rank}, residual HW: "
                  f"avg={avg_res:.1f}, min={min_res} (random=16)")

    # ── Summary ───────────────────────────────────────────────────
    print(f"\n{'='*72}")
    print("SUMMARY: Linear extraction W[14] → De17")
    print(f"{'='*72}")
    overall_avg = sum(all_residual_hws) / len(all_residual_hws)
    print(f"  GF(2) Jacobian rank: {all_ranks}")
    print(f"  Avg rank: {sum(all_ranks)/len(all_ranks):.1f}")
    print(f"  Residual HW: avg across msgs = {overall_avg:.1f}")
    print(f"  Expected if random: 16.0")
    print(f"  Expected if perfect linear: 0.0")

    savings = 16 - overall_avg
    print(f"\n  Linear approximation saves: {savings:.1f} bits")
    print(f"  Effective barrier cost: 2^{32-savings:.1f} instead of 2^32")

    # ── Step 4: What about SECOND-ORDER (quadratic) correction? ───
    print(f"\n{'='*72}")
    print("Quadratic correction: can we predict the residual?")
    print(f"{'='*72}")

    base_msg = [random.getrandbits(32) for _ in range(16)]
    base_de17 = get_de17(base_msg, iv)
    # Build Jacobian at this point
    jcols = []
    for j in range(32):
        msg_f = list(base_msg); msg_f[14] ^= (1<<j)
        jcols.append(base_de17 ^ get_de17(msg_f, iv))

    # Sample residuals and check if any residual BIT is predictable
    # from PAIRS of input bits (quadratic terms)
    N_Q = 3000
    inputs = []
    residuals = []
    for _ in range(N_Q):
        w14 = random.getrandbits(32)
        dw14 = base_msg[14] ^ w14
        lin = 0
        for j in range(32):
            if (dw14>>j)&1: lin ^= jcols[j]
        msg_n = list(base_msg); msg_n[14] = w14
        actual = get_de17(msg_n, iv)
        res = (actual ^ base_de17) ^ lin
        inputs.append(dw14)
        residuals.append(res)

    # Check: for output bit 0 of residual, does ANY pair of input bits predict it?
    best_quad_corr = 0
    for out_b in [0, 8, 16, 24, 31]:  # sample output bits
        for ib1 in range(0, 32, 4):  # sample input pairs
            for ib2 in range(ib1+1, 32, 4):
                match = 0
                for i in range(N_Q):
                    quad = ((inputs[i]>>ib1)&1) & ((inputs[i]>>ib2)&1)
                    out = (residuals[i]>>out_b)&1
                    if quad == out: match += 1
                c = abs(match/N_Q - 0.5)*2
                if c > best_quad_corr:
                    best_quad_corr = c

    print(f"  Best quadratic correlation: {best_quad_corr:.3f}")
    print(f"  Expected random: ~{2/(N_Q**0.5):.3f}")
    if best_quad_corr > 5/(N_Q**0.5):
        print(f"  ALIVE: quadratic terms have structure!")
    else:
        print(f"  DEAD: residual is random (no quadratic structure)")

    # ── Step 5: Multi-word linear map ─────────────────────────────
    print(f"\n{'='*72}")
    print("Multi-word: Combined W[13]+W[14]+W[15] → De17 linear map")
    print(f"{'='*72}")

    base_msg = [random.getrandbits(32) for _ in range(16)]
    base_de17 = get_de17(base_msg, iv)

    # Build 96-column Jacobian (W[13]: 32 cols, W[14]: 32 cols, W[15]: 32 cols)
    multi_cols = []
    for w in [13, 14, 15]:
        for j in range(32):
            msg_f = list(base_msg); msg_f[w] ^= (1<<j)
            multi_cols.append(base_de17 ^ get_de17(msg_f, iv))

    # Rank of 32×96 matrix
    rows96 = []
    for r in range(32):
        row = 0
        for c in range(96):
            if (multi_cols[c] >> r) & 1:
                row |= (1 << c)
        rows96.append(row)
    rows96c = list(rows96)
    rank96 = 0
    for col in range(96):
        pivot = -1
        for r in range(rank96, 32):
            if (rows96c[r] >> col) & 1: pivot = r; break
        if pivot == -1: continue
        rows96c[rank96], rows96c[pivot] = rows96c[pivot], rows96c[rank96]
        for r in range(32):
            if r != rank96 and (rows96c[r] >> col) & 1: rows96c[r] ^= rows96c[rank96]
        rank96 += 1

    print(f"  Combined Jacobian rank (32×96): {rank96}/32")
    print(f"  Kernel dimension: {96 - rank96}")
    print(f"  Free parameters after linear solve: {96 - rank96}")

    # Residual with combined linear map
    res_hws = []
    for _ in range(2000):
        msg_n = list(base_msg)
        for w in [13,14,15]: msg_n[w] = random.getrandbits(32)
        dw = []
        for w in [13,14,15]:
            for j in range(32):
                dw.append((base_msg[w] ^ msg_n[w] >> j) & 1)
        lin = 0
        for c in range(96):
            idx_w = c // 32  # which word
            idx_b = c % 32   # which bit
            w = [13,14,15][idx_w]
            if ((base_msg[w] ^ msg_n[w]) >> idx_b) & 1:
                lin ^= multi_cols[c]
        actual = get_de17(msg_n, iv)
        res = (actual ^ base_de17) ^ lin
        res_hws.append(hw(res))

    avg_multi = sum(res_hws)/len(res_hws)
    print(f"  Multi-word residual HW: avg={avg_multi:.1f}, min={min(res_hws)}")
    print(f"  Single-word (W[14] only): avg={overall_avg:.1f}")
    print(f"  Improvement: {overall_avg - avg_multi:.1f} bits")

    # ── VERDICT ───────────────────────────────────────────────────
    print(f"\n{'='*72}")
    print("VERDICT")
    print(f"{'='*72}\n")
    print(f"  Linear extraction: saves {savings:.1f} bits per barrier")
    print(f"  Effective single-barrier cost: 2^{32-savings:.0f}")
    if savings > 5:
        print(f"\n  **ALIVE**: Linear component is SIGNIFICANT!")
        print(f"  Strategy: solve linear system, brute-force residual 2^{32-savings:.0f}")
        print(f"  For 4 coupled barriers: 2^{4*(32-savings):.0f} via birthday")
    elif savings > 2:
        print(f"\n  **ANOMALY**: Small but real linear component")
    else:
        print(f"\n  **DEAD**: Linear approximation ≈ useless")

    print(f"\n  Runtime: {time.time()-t0:.1f}s")
    print("=" * 72)


if __name__ == "__main__":
    main()
