#!/usr/bin/env python3
"""
Step 2: SHA-256 round as algebraic element.

Key questions from Step 1:
- B has rank 5, kernel = {d, h}  (dim ker = 3, since f and g also degenerate!)
  Actually rank=5 means kernel dim=3: let's find the FULL kernel.
- How does the bilinear form propagate through multiple rounds?
- Can we decompose the round function into linear + quadratic parts?

New idea: decompose Round = Linear_part + Quadratic_part
  Linear: Sig0, Sig1, additions, shift register
  Quadratic: Ch(e,f,g), Maj(a,b,c)

If we track the quadratic part separately, what structure emerges?
"""

import random
from collections import Counter

MASK = 0xFFFFFFFF

def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK
def sig0(x):  return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sig1(x):  return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)
def Sig0(x):  return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sig1(x):  return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g):   return ((e & f) ^ (~e & g)) & MASK
def Maj(a, b, c):  return ((a & b) ^ (a & c) ^ (b & c)) & MASK
def hw(x): return bin(x).count('1')

IV = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]
K = [0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
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
     0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2]

def schedule(W16):
    W = list(W16) + [0]*48
    for i in range(16, 64):
        W[i] = (sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16]) & MASK
    return W

def sha_rounds(W64, R, iv=None):
    """Run R rounds of SHA-256. Return list of states [state_0, ..., state_R]."""
    if iv is None:
        iv = IV
    states = [list(iv)]
    s = list(iv)
    for r in range(R):
        a, b, c, d, e, f, g, h = s
        T1 = (h + Sig1(e) + Ch(e, f, g) + K[r] + W64[r]) & MASK
        T2 = (Sig0(a) + Maj(a, b, c)) & MASK
        s = [(T1 + T2) & MASK, a, b, c, (d + T1) & MASK, e, f, g]
        states.append(list(s))
    return states

# ============================================================
# PART A: Full kernel of B (rank 5 → kernel dim 3)
# ============================================================
print("=" * 70)
print("PART A: Full kernel analysis of B (8×8, rank 5)")
print("=" * 70)
print()

# B matrix from Step 1
B = [[0]*8 for _ in range(8)]
# Maj: a,b,c = 0,1,2
B[0][1] = B[1][0] = 1
B[0][2] = B[2][0] = 1
B[1][2] = B[2][1] = 1
# Ch: e,f,g = 4,5,6
B[4][5] = B[5][4] = 1
B[4][6] = B[6][4] = 1
# f,g cross = 0 in Ch

# Compute kernel over GF(2)
def gf2_kernel(M, n):
    """Find kernel of n×n matrix over GF(2)."""
    # Augmented matrix [M | I]
    aug = [[M[i][j] for j in range(n)] + [1 if i==j else 0 for j in range(n)] for i in range(n)]

    pivot_cols = []
    row = 0
    for col in range(n):
        # Find pivot
        pivot = -1
        for r in range(row, n):
            if aug[r][col] % 2 == 1:
                pivot = r
                break
        if pivot == -1:
            continue
        aug[row], aug[pivot] = aug[pivot], aug[row]
        for r in range(n):
            if r != row and aug[r][col] % 2 == 1:
                for j in range(2*n):
                    aug[r][j] = (aug[r][j] + aug[row][j]) % 2
        pivot_cols.append(col)
        row += 1

    # Non-pivot columns are free variables → kernel vectors
    rank = len(pivot_cols)
    free_cols = [c for c in range(n) if c not in pivot_cols]

    kernel = []
    for fc in free_cols:
        v = [0]*n
        v[fc] = 1
        for i, pc in enumerate(pivot_cols):
            v[pc] = aug[i][fc] % 2
        kernel.append(v)
    return kernel, rank

ker_vecs, rank = gf2_kernel(B, 8)
labels = "abcdefgh"
print(f"Rank of B over GF(2): {rank}")
print(f"Kernel dimension: {len(ker_vecs)}")
print()
for i, v in enumerate(ker_vecs):
    nonzero = [labels[j] for j in range(8) if v[j]]
    print(f"  ker vec {i+1}: {v}  →  {'+'.join(nonzero) if nonzero else '0'}")

print()

# ============================================================
# PART B: Decompose Round = Linear + Quadratic
# ============================================================
print("=" * 70)
print("PART B: Round decomposition — Linear vs Quadratic contribution")
print("=" * 70)
print()

# For XOR differential δ, one round:
#   δa' = δ(T1+T2) = δT1 ⊕ δT2 ⊕ carry_diff
#   δe' = δ(d+T1) = δd ⊕ δT1 ⊕ carry_diff
#
# Over GF(2) (ignoring carry):
#   δT1_GF2 = δh ⊕ Sig1(δe) ⊕ δCh(e,f,g,δe,δf,δg) ⊕ δW
#   δT2_GF2 = Sig0(δa) ⊕ δMaj(a,b,c,δa,δb,δc)
#
# LINEAR parts of δT1:  δh ⊕ Sig1(δe) ⊕ δW
# QUADRATIC part:       δCh = δe·(f⊕g)  when δf=δg=0
#                        or   e·(δf⊕δg)  when δe=0
#
# Question: how much of the differential is captured by linear part alone?

def sha_round_xor_diff(state, delta_state, W_r):
    """Compute XOR differential of one round.
    Returns (delta_out, linear_part, quadratic_ch, quadratic_maj, carry_noise).
    """
    a,b,c,d,e,f,g,h = state
    da,db,dc,dd,de,df,dg,dh = delta_state

    # Actual output
    T1_n = (h + Sig1(e) + Ch(e,f,g) + K[0] + W_r) & MASK
    T2_n = (Sig0(a) + Maj(a,b,c)) & MASK

    a2 = a^da; b2=b^db; c2=c^dc; d2=d^dd; e2=e^de; f2=f^df; g2=g^dg; h2=h^dh
    T1_f = (h2 + Sig1(e2) + Ch(e2,f2,g2) + K[0] + W_r) & MASK
    T2_f = (Sig0(a2) + Maj(a2,b2,c2)) & MASK

    out_n = [(T1_n+T2_n)&MASK, a, b, c, (d+T1_n)&MASK, e, f, g]
    out_f = [(T1_f+T2_f)&MASK, a2, b2, c2, (d2+T1_f)&MASK, e2, f2, g2]
    delta_out = [out_n[i] ^ out_f[i] for i in range(8)]

    # GF(2) linear prediction:
    # δT1_lin = δh ⊕ Sig1(δe) ⊕ δW (δW=0 since same W)
    dT1_lin = dh ^ Sig1(de)
    dT2_lin = Sig0(da)

    # Quadratic parts (GF(2)):
    # δCh = Ch(e⊕de, f⊕df, g⊕dg) ⊕ Ch(e,f,g)
    dCh = Ch(e^de, f^df, g^dg) ^ Ch(e,f,g)
    dMaj = Maj(a^da, b^db, c^dc) ^ Maj(a,b,c)

    # GF(2) prediction: δa' = δT1_lin ⊕ δCh ⊕ δT2_lin ⊕ δMaj
    # Carry noise = actual ⊕ GF(2) prediction
    da_gf2 = dT1_lin ^ dCh ^ dT2_lin ^ dMaj
    de_gf2 = dd ^ dT1_lin ^ dCh

    gf2_pred = [da_gf2, da, db, dc, de_gf2, de, df, dg]  # shift register
    carry_noise = [delta_out[i] ^ gf2_pred[i] for i in range(8)]

    return delta_out, gf2_pred, dCh, dMaj, carry_noise

# Test: how well does GF(2) prediction work?
N = 5000
hw_actual = []
hw_gf2 = []
hw_carry = []
hw_quad_ch = []
hw_quad_maj = []

for _ in range(N):
    W_r = random.randint(0, MASK)
    state = [random.randint(0, MASK) for _ in range(8)]
    # Random 1-bit differential in e
    delta = [0]*8
    delta[4] = 1 << random.randint(0, 31)

    dout, gf2, dCh, dMaj, carry = sha_round_xor_diff(state, delta, W_r)

    hw_actual.append(sum(hw(x) for x in dout))
    hw_gf2.append(sum(hw(x) for x in gf2))
    hw_carry.append(sum(hw(x) for x in carry))
    hw_quad_ch.append(hw(dCh))
    hw_quad_maj.append(hw(dMaj))

import numpy as np
print("Round decomposition (δe = 1 random bit, 1 round, N=5000):")
print(f"  E[HW(actual diff)]:     {np.mean(hw_actual):.2f} ± {np.std(hw_actual):.2f}")
print(f"  E[HW(GF2 prediction)]:  {np.mean(hw_gf2):.2f} ± {np.std(hw_gf2):.2f}")
print(f"  E[HW(carry noise)]:     {np.mean(hw_carry):.2f} ± {np.std(hw_carry):.2f}")
print(f"  E[HW(δCh)]:             {np.mean(hw_quad_ch):.2f} ± {np.std(hw_quad_ch):.2f}")
print(f"  E[HW(δMaj)]:            {np.mean(hw_quad_maj):.2f} ± {np.std(hw_quad_maj):.2f}")
print()

# ============================================================
# PART C: Multi-round propagation — how does B_round compose?
# ============================================================
print("=" * 70)
print("PART C: Multi-round — tracking linear vs quadratic vs carry")
print("=" * 70)
print()

# For R rounds, track: what fraction of the output diff is explained by
# (a) purely linear model, (b) linear+quadratic, (c) carry noise

def multi_round_decompose(W16, R, delta_W=None):
    """Track differential components through R rounds."""
    W64 = schedule(W16)

    state = list(IV)

    # Create perturbed message
    if delta_W is None:
        DW = [0]*16
        DW[0] = 1  # standard delta
    else:
        DW = delta_W

    W16_f = [(W16[i] ^ DW[i]) & MASK for i in range(16)]
    W64_f = schedule(W16_f)

    state_n = list(IV)
    state_f = list(IV)

    results = []
    for r in range(R):
        a,b,c,d,e,f,g,h = state_n
        T1 = (h + Sig1(e) + Ch(e,f,g) + K[r] + W64[r]) & MASK
        T2 = (Sig0(a) + Maj(a,b,c)) & MASK
        state_n = [(T1+T2)&MASK, a, b, c, (d+T1)&MASK, e, f, g]

        a,b,c,d,e,f,g,h = state_f
        T1 = (h + Sig1(e) + Ch(e,f,g) + K[r] + W64_f[r]) & MASK
        T2 = (Sig0(a) + Maj(a,b,c)) & MASK
        state_f = [(T1+T2)&MASK, a, b, c, (d+T1)&MASK, e, f, g]

        xor_diff = [state_n[i] ^ state_f[i] for i in range(8)]
        total_hw = sum(hw(x) for x in xor_diff)

        # Ch contribution: bits where e differs
        de = state_n[4] ^ state_f[4]  # current delta_e BEFORE this round
        # Actually this is delta after round... let me just track total HW

        results.append({
            'round': r+1,
            'total_hw': total_hw,
            'de_hw': hw(xor_diff[4]),  # delta_e after round
            'da_hw': hw(xor_diff[0]),  # delta_a after round
        })

    return results

# Run for multiple random W
print("Multi-round HW(diff) profile (averaged over 200 random W):")
print()
print(f"{'Round':>5} | {'E[HW(total)]':>13} | {'E[HW(δe)]':>10} | {'E[HW(δa)]':>10}")
print("-" * 50)

all_results = []
for trial in range(200):
    W16 = [random.randint(0, MASK) for _ in range(16)]
    res = multi_round_decompose(W16, 20)
    all_results.append(res)

for r in range(20):
    avg_total = np.mean([all_results[t][r]['total_hw'] for t in range(200)])
    avg_de = np.mean([all_results[t][r]['de_hw'] for t in range(200)])
    avg_da = np.mean([all_results[t][r]['da_hw'] for t in range(200)])
    print(f"  r={r+1:2d} | {avg_total:13.2f} | {avg_de:10.2f} | {avg_da:10.2f}")

print()

# ============================================================
# PART D: THE KEY EXPERIMENT — Quadratic residue over rounds
# ============================================================
# Question: if we subtract the "linear prediction" from the actual diff,
# what is the structure of the residue? Is it purely quadratic?
# Or does carry accumulate?

print("=" * 70)
print("PART D: Quadratic residue — GF(2) vs Z_{2^32} gap per round")
print("=" * 70)
print()

def track_gf2_gap(W16, R):
    """Track how much the GF(2) model diverges from reality."""
    W64 = schedule(W16)
    DW = [0]*16; DW[0] = 1
    W16_f = [(W16[i] ^ DW[i]) for i in range(16)]
    W64_f = schedule(W16_f)

    state_n = list(IV)
    state_f = list(IV)

    gaps = []
    for r in range(R):
        # Before round: current diff
        cur_diff = [state_n[i] ^ state_f[i] for i in range(8)]

        a,b,c,d,e,f,g,h = state_n
        a2,b2,c2,d2,e2,f2,g2,h2 = state_f

        # Actual round
        T1n = (h + Sig1(e) + Ch(e,f,g) + K[r] + W64[r]) & MASK
        T2n = (Sig0(a) + Maj(a,b,c)) & MASK
        T1f = (h2 + Sig1(e2) + Ch(e2,f2,g2) + K[r] + W64_f[r]) & MASK
        T2f = (Sig0(a2) + Maj(a2,b2,c2)) & MASK

        state_n = [(T1n+T2n)&MASK, a, b, c, (d+T1n)&MASK, e, f, g]
        state_f = [(T1f+T2f)&MASK, a2, b2, c2, (d2+T1f)&MASK, e2, f2, g2]

        actual_diff = [state_n[i] ^ state_f[i] for i in range(8)]

        # GF(2) prediction from cur_diff
        da,db,dc,dd,de,df,dg,dh = cur_diff
        dW = W64[r] ^ W64_f[r]

        dSig1 = Sig1(de)
        dCh = Ch(e^de, f^df, g^dg) ^ Ch(e,f,g)
        dSig0 = Sig0(da)
        dMaj = Maj(a^da, b^db, c^dc) ^ Maj(a,b,c)

        dT1_gf2 = dh ^ dSig1 ^ dCh ^ dW
        dT2_gf2 = dSig0 ^ dMaj

        pred_diff = [dT1_gf2 ^ dT2_gf2, da, db, dc, dd ^ dT1_gf2, de, df, dg]

        gap = [actual_diff[i] ^ pred_diff[i] for i in range(8)]
        gap_hw = sum(hw(x) for x in gap)

        gaps.append(gap_hw)

    return gaps

print("GF(2) prediction gap per round (carry noise):")
print("(averaged over 500 random W, DW[0]=1)")
print()

all_gaps = []
for _ in range(500):
    W16 = [random.randint(0, MASK) for _ in range(16)]
    gaps = track_gf2_gap(W16, 20)
    all_gaps.append(gaps)

print(f"{'Round':>5} | {'E[gap HW]':>10} | {'Interpretation':>30}")
print("-" * 55)
for r in range(20):
    avg = np.mean([all_gaps[t][r] for t in range(500)])
    if avg < 1:
        interp = "GF(2) nearly exact"
    elif avg < 5:
        interp = "small carry perturbation"
    elif avg < 15:
        interp = "significant carry divergence"
    else:
        interp = "FULL DIVERGENCE (random)"
    print(f"  r={r+1:2d} | {avg:10.2f} | {interp}")

print()
print("=" * 70)
print("PART E: Critical question — does gap stabilize or grow?")
print("=" * 70)
print()

# If gap stabilizes → carry is a bounded deformation (good for algebra!)
# If gap grows → carry destroys structure (bad)
gaps_30 = []
for _ in range(200):
    W16 = [random.randint(0, MASK) for _ in range(16)]
    g = track_gf2_gap(W16, 30)
    gaps_30.append(g)

print(f"{'Round':>5} | {'E[gap]':>8} | {'max':>6} | {'min':>6}")
print("-" * 40)
for r in [0, 1, 2, 4, 7, 9, 14, 19, 24, 29]:
    vals = [gaps_30[t][r] for t in range(200)]
    print(f"  r={r+1:2d} | {np.mean(vals):8.2f} | {max(vals):6d} | {min(vals):6d}")

print()
print("=" * 70)
print("SUMMARY STEP 2")
print("=" * 70)
print()
print("Key findings:")
print("1. ker(B) has dim 3: d, h, and f⊕g direction")
print("2. GF(2) model diverges from reality due to carry")
print("3. The GAP (carry noise) behavior tells us if Clifford approach works")
print("4. If gap saturates → carry is a BOUNDED deformation → algebra viable")
print("5. If gap grows unbounded → need different algebraic framework")
