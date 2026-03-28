#!/usr/bin/env python3
"""
Error Audit — Stage 3.1: Find false negatives and premature closures
=====================================================================

Systematic review of all 15 experiments from Stages 1-2.
For each: re-test the KEY conclusion with higher N or different method.

Known risks:
  R1: Stage 1.2 G3 "carry→H7 corr=0.005" — tested W[0] only, W[1..15]=0.
      All 16 words active could be different.
  R2: Stage 1.3 S2 "carry_sched→H7 corr=0.047" — N=8000, borderline.
      At N=50K this might be significant.
  R3: Stage 1.5 "r=30 signal" closed at N=100 — was it actually there?
  R4: Stage 2.1 "rank=53" — was rank 23 from influence matrix real?
  R5: Stage 2.6 "Ch is tautology" — is the Ch→H CONDITIONAL on carry=0
      a tautology, or does carry=0 create non-tautological structure?
  R6: All experiments used W[1..15]=0. Real SHA uses all 16 words.
      Are we studying a DEGENERATE CASE?
"""

import numpy as np
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
def hw(x): return bin(x & MASK).count('1')

def message_schedule(W16):
    W = list(W16) + [0] * 48
    for i in range(16, 64):
        W[i] = (sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16]) & MASK
    return W

def sha_with_raws(W16):
    W = message_schedule(W16)
    a, b, c, d, e, f, g, h = H0
    carries = []; raws = []
    for r in range(64):
        raw = h + Sig1(e) + Ch(e, f, g) + K[r] + W[r]
        carries.append(1 if raw >= T_VAL else 0)
        raws.append(raw)
        T1 = raw & MASK
        T2 = (Sig0(a) + Maj(a, b, c)) & MASK
        h, g, f = g, f, e
        e = (d + T1) & MASK
        d, c, b = c, b, a
        a = (T1 + T2) & MASK
    H_out = [(v + iv) & MASK for v, iv in zip([a,b,c,d,e,f,g,h], H0)]
    return carries, raws, H_out


# ============================================================
# R1: Was "carry→H7 = 0" tested correctly? (all 16 words)
# ============================================================

def audit_R1(N=15000, seed=900):
    print("=" * 70)
    print("R1 AUDIT: carry→H7 correlation — W[0] only vs all 16 words")
    print(f"N={N}, seed={seed}")
    print("=" * 70)

    rng = np.random.RandomState(seed)

    for config_name, gen_w in [
        ("W[0] only", lambda rng: [int(rng.randint(0,1<<32))] + [0]*15),
        ("all 16 random", lambda rng: [int(rng.randint(0,1<<32)) for _ in range(16)]),
        ("real padding 4B", lambda rng: [int(rng.randint(0,1<<32)), 0x80000000]+[0]*13+[32]),
    ]:
        # Collect carry profiles and H[7]
        carry_sums = []
        h7_vals = []
        h7_b29 = []
        h7_b31 = []
        raw63_vals = []

        t0 = time.time()
        for i in range(N):
            W16 = gen_w(rng)
            carries, raws, H_out = sha_with_raws(W16)
            carry_sums.append(sum(1-c for c in carries))
            h7_vals.append(H_out[7])
            h7_b29.append((H_out[7] >> 29) & 1)
            h7_b31.append((H_out[7] >> 31) & 1)
            raw63_vals.append(raws[63])

        cs = np.array(carry_sums, dtype=float)
        h7 = np.array(h7_vals, dtype=float)
        r63 = np.array(raw63_vals, dtype=float)
        b29 = np.array(h7_b29, dtype=float)
        b31 = np.array(h7_b31, dtype=float)

        corr_cs_h7 = np.corrcoef(cs, h7)[0,1]
        corr_r63_h7 = np.corrcoef(r63, h7)[0,1]
        corr_r63_b31 = np.corrcoef(r63, b31)[0,1]

        # Per-round carry→H7 correlation
        per_round_corr = []
        for r_test in [0, 9, 30, 47, 62, 63]:
            carry_r = []
            for i in range(N):
                W16 = gen_w(np.random.RandomState(seed + i))
                c, _, _ = sha_with_raws(W16)
                carry_r.append(c[r_test])
            # Too slow for full — use sample
            break  # Skip per-round, too expensive

        print(f"\n  {config_name} (N={N}, {time.time()-t0:.1f}s):")
        print(f"    corr(carry_sum, H7):    {corr_cs_h7:+.5f}")
        print(f"    corr(raw63, H7):        {corr_r63_h7:+.5f}")
        print(f"    corr(raw63, H7[b31]):   {corr_r63_b31:+.5f}")
        print(f"    E[carry_sum]:           {np.mean(cs):.2f}")
        print(f"    P(carry63=0):           {np.mean(np.array(raw63_vals) < T_VAL):.5f}")

    return


# ============================================================
# R5: Ch→H CONDITIONAL on carry=0 — tautology or real?
# ============================================================

def audit_R5(N=20000, seed=901):
    print("\n" + "=" * 70)
    print("R5 AUDIT: Is Ch→H conditional on carry=0 a tautology?")
    print(f"N={N}, seed={seed}")
    print("=" * 70)

    rng = np.random.RandomState(seed)

    # The argument: Ch contains g, g becomes H[7], so corr(Ch, H7) is tautological.
    # BUT: is the correlation STRONGER when carry[63]=0?
    # If yes → carry=0 creates ADDITIONAL structure beyond the tautology.
    # If same → pure tautology, carry=0 adds nothing.

    ch_vals = []
    h7_vals = []
    carry63 = []

    t0 = time.time()
    for i in range(N):
        W16 = [int(rng.randint(0, 1 << 32)) for _ in range(16)]
        carries, raws, H_out = sha_with_raws(W16)

        e62 = (H_out[5] - H0[5]) & MASK
        f62 = (H_out[6] - H0[6]) & MASK
        g62 = (H_out[7] - H0[7]) & MASK
        ch62 = Ch(e62, f62, g62)

        ch_vals.append(ch62)
        h7_vals.append(H_out[7])
        carry63.append(carries[63])

    ch_arr = np.array(ch_vals, dtype=float)
    h7_arr = np.array(h7_vals, dtype=float)
    c63_arr = np.array(carry63)

    # Unconditional correlation
    corr_uncond = np.corrcoef(ch_arr, h7_arr)[0, 1]

    # Correlation conditioned on carry[63]=0
    mask_c0 = c63_arr == 0
    n_c0 = np.sum(mask_c0)
    if n_c0 > 10:
        corr_c0 = np.corrcoef(ch_arr[mask_c0], h7_arr[mask_c0])[0, 1]
    else:
        corr_c0 = float('nan')

    # Correlation conditioned on carry[63]=1
    mask_c1 = c63_arr == 1
    corr_c1 = np.corrcoef(ch_arr[mask_c1], h7_arr[mask_c1])[0, 1]

    print(f"\n  Collected: {time.time()-t0:.1f}s, N(c63=0)={n_c0}")
    print(f"\n  corr(Ch62, H7) unconditional: {corr_uncond:+.5f}")
    print(f"  corr(Ch62, H7) | carry63=1:   {corr_c1:+.5f}")
    print(f"  corr(Ch62, H7) | carry63=0:   {corr_c0:+.5f}" if not np.isnan(corr_c0) else "  corr | c63=0: N/A (too few)")

    if not np.isnan(corr_c0):
        print(f"\n  Delta (c0 - c1): {corr_c0 - corr_c1:+.5f}")
        if abs(corr_c0 - corr_c1) > 0.05:
            print(f"  ★ carry=0 CHANGES the Ch-H7 relationship — NOT pure tautology!")
        else:
            print(f"  Carry=0 doesn't change relationship — pure tautological correlation")

    # DEEPER: does Ch→H hold for a RANDOM function with same structure?
    # Construct: random_H[7] = random_g62 + IV[7], Ch = Ch(e62, f62, random_g62)
    # The tautological component: corr(Ch(e,f,g), g+const) for random e,f,g
    print(f"\n  --- Tautological baseline: corr(Ch(e,f,g), g) for random e,f,g ---")
    n_taut = 50000
    rng2 = np.random.RandomState(seed + 100)
    e_rand = rng2.randint(0, 1 << 32, n_taut).astype(np.int64)
    f_rand = rng2.randint(0, 1 << 32, n_taut).astype(np.int64)
    g_rand = rng2.randint(0, 1 << 32, n_taut).astype(np.int64)
    ch_rand = np.array([(int(e) & int(f)) ^ (~int(e) & int(g)) & MASK for e, f, g in zip(e_rand, f_rand, g_rand)], dtype=float)
    g_plus_const = np.array([(int(g) + H0[7]) & MASK for g in g_rand], dtype=float)

    corr_taut = np.corrcoef(ch_rand, g_plus_const)[0, 1]
    print(f"  corr(Ch(e,f,g), g+IV7) for random uniform: {corr_taut:+.5f}")
    print(f"  corr(Ch, H7) in SHA-256:                    {corr_uncond:+.5f}")
    print(f"  Delta (SHA - tautological):                  {corr_uncond - corr_taut:+.5f}")

    if abs(corr_uncond - corr_taut) > 0.02:
        print(f"\n  ★★ SHA-256 creates ADDITIONAL correlation beyond tautological baseline!")
        print(f"  The excess {corr_uncond - corr_taut:+.4f} is a REAL signal, not tautology.")
    else:
        print(f"\n  SHA-256 correlation = tautological baseline. Pure algebraic identity.")

    # Per-bit analysis: Ch[b31]→H7[b31]
    print(f"\n  --- Per-bit tautological analysis ---")
    for bit in [29, 30, 31]:
        # Tautological: P(Ch[bit]=0) → P(H7[bit]=1) for random Ch,g
        ch_bit_rand = np.array([(int(c) >> bit) & 1 for c in ch_rand], dtype=float)
        g_bit_rand = np.array([(int(g) >> bit) & 1 for g in g_plus_const], dtype=float)
        corr_bit_taut = np.corrcoef(ch_bit_rand, g_bit_rand)[0, 1]

        ch_bit_sha = np.array([(c >> bit) & 1 for c in ch_vals], dtype=float)
        h7_bit_sha = np.array([(h >> bit) & 1 for h in h7_vals], dtype=float)
        corr_bit_sha = np.corrcoef(ch_bit_sha, h7_bit_sha)[0, 1] if np.std(ch_bit_sha) > 0.01 else 0

        delta = corr_bit_sha - corr_bit_taut
        marker = " ★★" if abs(delta) > 0.02 else (" ★" if abs(delta) > 0.01 else "")
        print(f"    bit {bit}: taut={corr_bit_taut:+.4f}, SHA={corr_bit_sha:+.4f}, delta={delta:+.4f}{marker}")

    return corr_uncond, corr_taut


# ============================================================
# R6: Degenerate case W[1..15]=0 vs real messages
# ============================================================

def audit_R6(N=10000, seed=902):
    print("\n" + "=" * 70)
    print("R6 AUDIT: Are we studying a degenerate case?")
    print(f"N={N}, seed={seed}")
    print("=" * 70)

    rng = np.random.RandomState(seed)

    configs = {
        'W0_only':     lambda: [int(rng.randint(0,1<<32))] + [0]*15,
        'all_16':      lambda: [int(rng.randint(0,1<<32)) for _ in range(16)],
        'real_4B':     lambda: [int(rng.randint(0,1<<32)), 0x80000000] + [0]*13 + [32],
        'real_16B':    lambda: [int(rng.randint(0,1<<32)) for _ in range(4)] + [0x80000000] + [0]*10 + [128],
        'real_55B':    lambda: [int(rng.randint(0,1<<32)) for _ in range(14)] + [0x80000000, 440],
    }

    print(f"\n  Comparing carry structure across message types:")
    print(f"  {'config':>12} | {'E[cs]':>6} | {'P(c63=0)':>9} | {'E[raw63]/T':>11} | {'std(raw63)/T':>13} | {'corr(r62,r63)':>14}")
    print(f"  {'-'*12}-+-{'-'*6}-+-{'-'*9}-+-{'-'*11}-+-{'-'*13}-+-{'-'*14}")

    for cname, gen_w in configs.items():
        raws_63 = []
        raws_62 = []
        cs_list = []

        for _ in range(N):
            carries, raws, _ = sha_with_raws(gen_w())
            raws_63.append(raws[63])
            raws_62.append(raws[62])
            cs_list.append(sum(1-c for c in carries))

        r63 = np.array(raws_63)
        r62 = np.array(raws_62)
        cs = np.array(cs_list)

        corr_62_63 = np.corrcoef(r62, r63)[0, 1]
        print(f"  {cname:>12} | {np.mean(cs):6.2f} | {np.mean(r63 < T_VAL):9.5f} | {np.mean(r63)/T_VAL:11.4f} | {np.std(r63)/T_VAL:13.4f} | {corr_62_63:+14.4f}")

    # KEY: do the copula model parameters change?
    print(f"\n  KEY: Does σ (std of raw) change across configs?")
    print(f"  If σ is constant → model is universal")
    print(f"  If σ changes → W[1..15]=0 was degenerate")

    # Compare block structure
    print(f"\n  Block isolation comparison:")
    for cname, gen_w in [('W0_only', configs['W0_only']), ('all_16', configs['all_16'])]:
        raws_all = np.zeros((min(N, 5000), 64))
        for i in range(min(N, 5000)):
            _, raws, _ = sha_with_raws(gen_w())
            raws_all[i] = raws

        # Within-block (W5: r=47,48) vs between-block (W2→W5: r=10,47)
        within = np.corrcoef(raws_all[:, 47], raws_all[:, 48])[0, 1]
        between = np.corrcoef(raws_all[:, 10], raws_all[:, 47])[0, 1]
        isolation = abs(within) / max(abs(between), 1e-5)
        print(f"    {cname:>12}: within={within:+.4f}, between={between:+.4f}, isolation={isolation:.1f}×")

    return


# ============================================================
# MAIN
# ============================================================

if __name__ == '__main__':
    print("Error Audit — Stage 3.1")
    print(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print()

    N = 15000
    if len(sys.argv) > 1 and sys.argv[1] == '--fast':
        N = 8000

    t_start = time.time()
    audit_R1(N=N)
    corr_sha, corr_taut = audit_R5(N=N)
    audit_R6(N=N)

    total = time.time() - t_start
    print(f"\n{'='*70}")
    print(f"TOTAL TIME: {total:.1f}s")
    print(f"{'='*70}")

    print(f"""
{'='*70}
AUDIT RESULTS — Stage 3.1
{'='*70}

R1: carry/raw → H7 across configs
  W0_only and all_16 should show SAME or DIFFERENT correlations.
  If different → our W[1..15]=0 experiments missed structure.

R5: Ch→H7 tautological baseline
  Random Ch(e,f,g) vs g: corr = {corr_taut:+.4f}
  SHA-256 Ch62 vs H7:    corr = {corr_sha:+.4f}
  Delta (excess):        {corr_sha - corr_taut:+.4f}
  If delta > 0.02 → SHA-256 adds real signal beyond tautology.

R6: Degenerate case check
  If σ(raw63) changes between W0-only and all-16 → model breaks.
  If block isolation changes → structure is config-dependent.

ERRORS FOUND: [see output above]
""")
