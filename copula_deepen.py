#!/usr/bin/env python3
"""
Copula Deepening — Stage 2.5: Three unexplored directions
============================================================

The copula model is complete but found no shortcut.
Three directions that the model HASN'T explored:

  D1: CROSS-WORD COPULA — we tested W[0] only (W[1..15]=0).
      With ALL 16 words active, blocks may MERGE.
      If blocks overlap → isolation breaks → new structure.

  D2: STATE COPULA — raw[r] = SS[r] + K[r] + W[r].
      We modeled raw. But SS[r] = h + Sig1(e) + Ch(e,f,g)
      itself has INTERNAL structure: h, Sig1(e), Ch are NOT
      independent. Their copula may reveal hidden coupling.

  D3: H-COPULA — instead of raw→H (one-directional bridge),
      model the JOINT distribution (raw[61], raw[62], raw[63], H[5..7]).
      The 6-dimensional copula may reveal structure invisible
      in marginal correlations.
"""

import numpy as np
from math import erfc, sqrt
import time
import sys

MASK = 0xFFFFFFFF
T_VAL = float(1 << 32)

H0 = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
]

K = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
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
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
]

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

def full_trace(W16):
    """Returns raws, SS components, H, state at round 62."""
    W = message_schedule(W16)
    a, b, c, d, e, f, g, h = H0
    raws = []
    ss_h = []  # h component of SS
    ss_s1 = []  # Sig1(e) component
    ss_ch = []  # Ch(e,f,g) component
    states = [(a,b,c,d,e,f,g,h)]

    for r in range(64):
        s1_val = Sig1(e)
        ch_val = Ch(e, f, g)
        raw = h + s1_val + ch_val + K[r] + W[r]
        raws.append(raw)
        ss_h.append(h)
        ss_s1.append(s1_val)
        ss_ch.append(ch_val)

        T1 = raw & MASK
        T2 = (Sig0(a) + Maj(a, b, c)) & MASK
        h, g, f = g, f, e
        e = (d + T1) & MASK
        d, c, b = c, b, a
        a = (T1 + T2) & MASK
        states.append((a,b,c,d,e,f,g,h))

    H_out = [(v + iv) & MASK for v, iv in zip([a,b,c,d,e,f,g,h], H0)]
    return raws, ss_h, ss_s1, ss_ch, W, H_out, states


# ============================================================
# D1: Cross-Word Copula — all 16 words active
# ============================================================

def experiment_D1(N=8000, seed=700):
    print("=" * 70)
    print("D1: Cross-Word Copula — All 16 W[k] active")
    print(f"N={N}, seed={seed}")
    print("=" * 70)

    rng = np.random.RandomState(seed)

    # Compare: W[0] only vs all 16 random
    configs = {
        'W0_only': lambda rng: [int(rng.randint(0, 1<<32))] + [0]*15,
        'all_16':  lambda rng: [int(rng.randint(0, 1<<32)) for _ in range(16)],
    }

    for cname, gen_W in configs.items():
        all_raws = np.zeros((N, 64))
        t0 = time.time()
        for i in range(N):
            W16 = gen_W(rng)
            raws, _, _, _, _, _, _ = full_trace(W16)
            all_raws[i] = raws

        # Block correlations
        blocks = {'W1': [0,1], 'W2': [9,10], 'W3': [18,19], 'W4': [30,31], 'W5': [47,48], 'tail': [62,63]}

        within = []
        between = []
        block_names = list(blocks.keys())
        for bn in block_names:
            rs = blocks[bn]
            if len(rs) >= 2:
                c = np.corrcoef(all_raws[:, rs[0]], all_raws[:, rs[1]])[0, 1]
                within.append(c)

        for i, b1 in enumerate(block_names):
            for b2 in block_names[i+1:]:
                c = np.corrcoef(all_raws[:, blocks[b1][0]], all_raws[:, blocks[b2][0]])[0, 1]
                between.append(c)

        mean_within = np.mean(within)
        mean_between = np.mean([abs(b) for b in between])
        isolation = mean_within / max(mean_between, 1e-10)

        # Cross-block for consecutive windows
        cross_w2_w3 = np.corrcoef(all_raws[:, 10], all_raws[:, 18])[0, 1]
        cross_w4_w5 = np.corrcoef(all_raws[:, 31], all_raws[:, 47])[0, 1]

        print(f"\n  {cname}:")
        print(f"    Within-block corr:  {mean_within:+.4f}")
        print(f"    Between-block corr: {mean_between:.4f}")
        print(f"    Isolation ratio:    {isolation:.1f}×")
        print(f"    W2→W3 cross:        {cross_w2_w3:+.4f}")
        print(f"    W4→W5 cross:        {cross_w4_w5:+.4f}")
        print(f"    P(carry[63]=0):     {np.mean(all_raws[:, 63] < T_VAL):.5f}")

    # KEY TEST: does full W[0..15] create cross-block correlations?
    print(f"\n  KEY: Does all_16 break block isolation?")
    print(f"  If isolation drops significantly → blocks merge → new structure")

    return


# ============================================================
# D2: State Copula — Internal structure of SS[r]
# ============================================================

def experiment_D2(N=10000, seed=701):
    print("\n" + "=" * 70)
    print("D2: State Copula — Internal structure of SS = h + Sig1(e) + Ch")
    print(f"N={N}, seed={seed}")
    print("=" * 70)

    rng = np.random.RandomState(seed)

    # Collect SS components for last 4 rounds
    data = {r: {'h': [], 's1': [], 'ch': [], 'ss': [], 'raw': []} for r in [60,61,62,63]}
    h7_vals = []

    t0 = time.time()
    for i in range(N):
        W16 = [int(rng.randint(0, 1 << 32)) for _ in range(16)]
        raws, ss_h, ss_s1, ss_ch, W, H_out, states = full_trace(W16)

        for r in [60,61,62,63]:
            data[r]['h'].append(ss_h[r])
            data[r]['s1'].append(ss_s1[r])
            data[r]['ch'].append(ss_ch[r])
            data[r]['ss'].append(ss_h[r] + ss_s1[r] + ss_ch[r])
            data[r]['raw'].append(raws[r])
        h7_vals.append(H_out[7])

    print(f"Collected: {time.time()-t0:.1f}s")

    # For each round: correlation matrix of (h, Sig1, Ch, W[r])
    for r in [62, 63]:
        h_arr = np.array(data[r]['h'], dtype=float)
        s1_arr = np.array(data[r]['s1'], dtype=float)
        ch_arr = np.array(data[r]['ch'], dtype=float)
        ss_arr = np.array(data[r]['ss'], dtype=float)
        raw_arr = np.array(data[r]['raw'], dtype=float)
        h7_arr = np.array(h7_vals, dtype=float)

        print(f"\n  --- Round {r}: SS = h + Sig1(e) + Ch(e,f,g) ---")

        vars_list = [('h', h_arr), ('Sig1', s1_arr), ('Ch', ch_arr)]
        print(f"    Correlation matrix:")
        print(f"    {'':>6}", end='')
        for name, _ in vars_list:
            print(f" {name:>8}", end='')
        print(f" {'SS':>8} {'raw':>8} {'H7':>8}")

        for name_i, arr_i in vars_list:
            print(f"    {name_i:>6}", end='')
            for name_j, arr_j in vars_list:
                c = np.corrcoef(arr_i, arr_j)[0, 1]
                print(f" {c:+8.3f}", end='')
            c_ss = np.corrcoef(arr_i, ss_arr)[0, 1]
            c_raw = np.corrcoef(arr_i, raw_arr)[0, 1]
            c_h7 = np.corrcoef(arr_i, h7_arr)[0, 1]
            print(f" {c_ss:+8.3f} {c_raw:+8.3f} {c_h7:+8.3f}")

        # Variance decomposition
        var_h = np.var(h_arr)
        var_s1 = np.var(s1_arr)
        var_ch = np.var(ch_arr)
        var_ss = np.var(ss_arr)
        cov_h_s1 = np.cov(h_arr, s1_arr)[0, 1]
        cov_h_ch = np.cov(h_arr, ch_arr)[0, 1]
        cov_s1_ch = np.cov(s1_arr, ch_arr)[0, 1]

        print(f"\n    Variance decomposition of SS[{r}]:")
        print(f"      var(h)/var(SS)    = {var_h/var_ss*100:5.1f}%")
        print(f"      var(Sig1)/var(SS) = {var_s1/var_ss*100:5.1f}%")
        print(f"      var(Ch)/var(SS)   = {var_ch/var_ss*100:5.1f}%")
        print(f"      2cov(h,S1)/var    = {2*cov_h_s1/var_ss*100:+5.1f}%")
        print(f"      2cov(h,Ch)/var    = {2*cov_h_ch/var_ss*100:+5.1f}%")
        print(f"      2cov(S1,Ch)/var   = {2*cov_s1_ch/var_ss*100:+5.1f}%")
        print(f"      Sum               = {(var_h+var_s1+var_ch+2*cov_h_s1+2*cov_h_ch+2*cov_s1_ch)/var_ss*100:5.1f}%")

    # KEY: which SS component best predicts H[7]?
    print(f"\n  --- Which SS component predicts H[7]? ---")
    h7 = np.array(h7_vals, dtype=float)
    for r in [62, 63]:
        for name, key in [('h', 'h'), ('Sig1', 's1'), ('Ch', 'ch'), ('SS', 'ss'), ('raw', 'raw')]:
            arr = np.array(data[r][key], dtype=float)
            c = np.corrcoef(arr, h7)[0, 1]
            # Also per-bit
            h7_b29 = np.array([(int(v) >> 29) & 1 for v in h7_vals], dtype=float)
            h7_b31 = np.array([(int(v) >> 31) & 1 for v in h7_vals], dtype=float)
            c29 = np.corrcoef(arr, h7_b29)[0, 1] if np.std(h7_b29) > 0.01 else 0
            c31 = np.corrcoef(arr, h7_b31)[0, 1] if np.std(h7_b31) > 0.01 else 0
            marker = " ★" if abs(c) > 0.05 else ""
            print(f"    r={r} {name:>4} → H7: {c:+.4f}, → b29: {c29:+.4f}, → b31: {c31:+.4f}{marker}")

    return data


# ============================================================
# D3: H-Copula — Joint (raw[61..63], H[5..7])
# ============================================================

def experiment_D3(N=10000, seed=702):
    print("\n" + "=" * 70)
    print("D3: H-Copula — Joint distribution of (raw[61..63], H[5..7])")
    print(f"N={N}, seed={seed}")
    print("=" * 70)

    rng = np.random.RandomState(seed)

    raw61 = []; raw62 = []; raw63 = []
    h5 = []; h6 = []; h7 = []

    t0 = time.time()
    for i in range(N):
        W16 = [int(rng.randint(0, 1 << 32)) for _ in range(16)]
        raws, _, _, _, _, H_out, _ = full_trace(W16)
        raw61.append(raws[61]); raw62.append(raws[62]); raw63.append(raws[63])
        h5.append(H_out[5]); h6.append(H_out[6]); h7.append(H_out[7])

    print(f"Collected: {time.time()-t0:.1f}s")

    # 6×6 correlation matrix
    names = ['raw61', 'raw62', 'raw63', 'H5', 'H6', 'H7']
    arrays = [np.array(x, dtype=float) for x in [raw61, raw62, raw63, h5, h6, h7]]

    print(f"\n  6×6 Correlation Matrix:")
    print(f"  {'':>7}", end='')
    for n in names:
        print(f" {n:>7}", end='')
    print()
    for i, n_i in enumerate(names):
        print(f"  {n_i:>7}", end='')
        for j, n_j in enumerate(names):
            c = np.corrcoef(arrays[i], arrays[j])[0, 1]
            marker = '*' if abs(c) > 0.05 and i != j else ' '
            print(f" {c:+6.3f}{marker}", end='')
        print()

    # Tail dependence in H-space: when raw[63] small, are H values correlated?
    print(f"\n  --- Tail conditioned: raw[63] in bottom 5% ---")
    q5 = np.percentile(arrays[2], 5)
    mask = arrays[2] < q5
    n_tail = np.sum(mask)
    print(f"  N in tail: {n_tail}")

    if n_tail > 20:
        print(f"\n  H correlation CONDITIONED on raw[63] < 5th percentile:")
        h_names = ['H5', 'H6', 'H7']
        h_arrays = [arrays[3][mask], arrays[4][mask], arrays[5][mask]]
        for i, n_i in enumerate(h_names):
            for j, n_j in enumerate(h_names):
                if j > i:
                    c_cond = np.corrcoef(h_arrays[i], h_arrays[j])[0, 1]
                    c_uncond = np.corrcoef(arrays[3+i], arrays[3+j])[0, 1]
                    print(f"    corr({n_i},{n_j}): uncond={c_uncond:+.4f}, cond={c_cond:+.4f}, delta={c_cond-c_uncond:+.4f}")

    # Per-bit tail analysis: which H bits cluster in the tail?
    print(f"\n  --- H[7] bit probabilities in raw[63] tail vs baseline ---")
    for bit in range(24, 32):
        h7_bit_all = np.array([(int(v) >> bit) & 1 for v in h7], dtype=float)
        h7_bit_tail = np.array([(int(h7[i]) >> bit) & 1 for i in range(N) if mask[i]], dtype=float)
        p_all = np.mean(h7_bit_all)
        p_tail = np.mean(h7_bit_tail)
        delta = p_tail - p_all
        marker = " ★★" if abs(delta) > 0.05 else (" ★" if abs(delta) > 0.02 else "")
        print(f"    H7[b{bit}]: baseline={p_all:.4f}, tail={p_tail:.4f}, delta={delta:+.4f}{marker}")

    # Same for H[6]
    print(f"\n  --- H[6] bit probabilities in raw[63] tail vs baseline ---")
    for bit in range(24, 32):
        h6_bit_all = np.array([(int(v) >> bit) & 1 for v in h6], dtype=float)
        h6_bit_tail = np.array([(int(h6[i]) >> bit) & 1 for i in range(N) if mask[i]], dtype=float)
        p_all = np.mean(h6_bit_all)
        p_tail = np.mean(h6_bit_tail)
        delta = p_tail - p_all
        marker = " ★★" if abs(delta) > 0.05 else (" ★" if abs(delta) > 0.02 else "")
        print(f"    H6[b{bit}]: baseline={p_all:.4f}, tail={p_tail:.4f}, delta={delta:+.4f}{marker}")

    # JOINT tail: raw[62] AND raw[63] both small
    print(f"\n  --- Joint tail: raw[62] < 5% AND raw[63] < 5% ---")
    q5_62 = np.percentile(arrays[1], 5)
    joint_mask = (arrays[1] < q5_62) & (arrays[2] < q5)
    n_joint = np.sum(joint_mask)
    print(f"  N in joint tail: {n_joint} (expected if indep: {N*0.05*0.05:.0f})")

    if n_joint > 5:
        for bit in [29, 30, 31]:
            h7_bit_joint = np.array([(int(h7[i]) >> bit) & 1 for i in range(N) if joint_mask[i]], dtype=float)
            h7_bit_single = np.array([(int(h7[i]) >> bit) & 1 for i in range(N) if mask[i]], dtype=float)
            p_joint = np.mean(h7_bit_joint) if len(h7_bit_joint) > 0 else 0.5
            p_single = np.mean(h7_bit_single)
            p_base = np.mean([(int(v) >> bit) & 1 for v in h7])
            print(f"    H7[b{bit}]: base={p_base:.3f}, single_tail={p_single:.3f}, joint_tail={p_joint:.3f}, amplif={abs(p_joint-0.5)/max(abs(p_single-0.5),0.001):.2f}×")

    return


# ============================================================
# MAIN
# ============================================================

if __name__ == '__main__':
    print("Copula Deepening — Stage 2.5")
    print(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print()

    N1, N2, N3 = 6000, 8000, 8000
    if len(sys.argv) > 1 and sys.argv[1] == '--fast':
        N1, N2, N3 = 3000, 5000, 5000

    t_start = time.time()
    experiment_D1(N=N1)
    data_d2 = experiment_D2(N=N2)
    experiment_D3(N=N3)

    total = time.time() - t_start
    print(f"\n{'='*70}")
    print(f"TOTAL TIME: {total:.1f}s")
    print(f"{'='*70}")

    print(f"""
{'='*70}
SYNTHESIS: Copula Deepening (Stage 2.5)
{'='*70}

D1: Does full W[0..15] break block isolation?
D2: Which SS component (h, Sig1, Ch) carries the bridge?
D3: Joint H-copula — amplification in the joint tail?

These three directions probe DEEPER than the basic copula model.
If any shows new structure → path to stronger distinguisher.
If all three confirm the model → copula theory is complete.
""")
