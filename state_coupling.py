#!/usr/bin/env python3
"""
State Coupling + Unexplored Angles — Stage 3.3/3.4
====================================================

Stage 2.0 W3 found: conditioning on W[0] does NOT remove carry
correlations. 52% residual for (carry[1], carry[48]). State creates
independent coupling BEYOND schedule.

Stage 2.5 found: Ch is tautological for H[7]. But what about
the a-BRANCH? H[0..3] come from a,b,c,d — through Maj, Sig0.
Nobody tested whether Maj creates a parallel bridge to H[0..3].

THREE UNEXPLORED DIRECTIONS:

  U1: STATE COUPLING EXPLOITATION
      carry[r1] and carry[r2] are correlated through state even
      after conditioning on W[0]. Can we use this to amplify?
      Specifically: if we know carry[1]=0, does that help predict
      carry[63]=0 through state (not schedule)?

  U2: A-BRANCH BRIDGE — H[0..3] through Maj(a,b,c)
      H[7] = g62 + IV[7] → Ch bridge (tautological)
      H[0] = a63 + IV[0] → Maj bridge? (unexplored!)
      Maj(a,b,c) = (a&b)^(a&c)^(b&c) — does it create signal in H[0]?

  U3: CROSS-WORD DISTINGUISHER — joint H[0..7] structure
      All 8 output words come from state[63] + IV.
      Under carry=0 conditions, state is constrained.
      Do ALL 8 words show joint structure?
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

def message_schedule(W16):
    W = list(W16) + [0] * 48
    for i in range(16, 64):
        W[i] = (sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16]) & MASK
    return W

def full_sha_detailed(W16):
    """Returns carries, raws, H_out, and internal state components."""
    W = message_schedule(W16)
    a, b, c, d, e, f, g, h = H0
    carries = []; raws = []
    # Track T2 components for a-branch analysis
    t2_sig0 = []; t2_maj = []; t2_vals = []
    states = [(a,b,c,d,e,f,g,h)]

    for r in range(64):
        s1 = Sig1(e); ch = Ch(e, f, g)
        raw = h + s1 + ch + K[r] + W[r]
        carries.append(1 if raw >= T_VAL else 0)
        raws.append(raw)

        T1 = raw & MASK
        s0 = Sig0(a); maj = Maj(a, b, c)
        T2 = (s0 + maj) & MASK
        t2_sig0.append(s0); t2_maj.append(maj); t2_vals.append(T2)

        h, g, f = g, f, e
        e = (d + T1) & MASK
        d, c, b = c, b, a
        a = (T1 + T2) & MASK
        states.append((a,b,c,d,e,f,g,h))

    H_out = [(states[64][i] + H0[i]) & MASK for i in range(8)]
    return carries, raws, H_out, states, t2_sig0, t2_maj, t2_vals, W


# ============================================================
# U1: State Coupling — carry[r1] helps predict carry[r2]?
# ============================================================

def experiment_U1(N=15000, seed=1000):
    print("=" * 70)
    print("U1: State Coupling — Does carry[r1]=0 help predict carry[r2]?")
    print(f"N={N}, seed={seed}")
    print("=" * 70)

    rng = np.random.RandomState(seed)

    # Collect carry profiles and raw values
    all_carries = np.zeros((N, 64), dtype=np.int8)
    all_raws = np.zeros((N, 64))

    t0 = time.time()
    for i in range(N):
        W16 = [int(rng.randint(0, 1 << 32)) for _ in range(16)]
        carries, raws, _, _, _, _, _, _ = full_sha_detailed(W16)
        all_carries[i] = carries
        all_raws[i] = raws

    print(f"Collected: {time.time()-t0:.1f}s")

    # For each pair of "easy" rounds (low K), measure conditional P
    easy_rounds = [0, 4, 9, 10, 18, 19, 30, 31, 47, 48, 49]

    print(f"\n--- Conditional carry: P(carry[r2]=0 | carry[r1]=0) ---")
    print(f"  Baseline P(carry[r]=0) for reference:")
    for r in easy_rounds:
        p = np.mean(1 - all_carries[:, r])
        print(f"    r={r:2d}: P={p:.4f}")

    print(f"\n  {'r1→r2':>8} | {'P(r2=0|r1=0)':>13} | {'P(r2=0) base':>13} | {'lift':>6} | {'N(r1=0)':>8}")
    print(f"  {'-'*8}-+-{'-'*13}-+-{'-'*13}-+-{'-'*6}-+-{'-'*8}")

    strong_couplings = []
    for r1 in easy_rounds:
        mask_r1 = all_carries[:, r1] == 0
        n_r1 = np.sum(mask_r1)
        if n_r1 < 10:
            continue
        for r2 in easy_rounds:
            if r2 == r1:
                continue
            p_cond = np.mean(1 - all_carries[mask_r1, r2])
            p_base = np.mean(1 - all_carries[:, r2])
            lift = p_cond / max(p_base, 1e-5)
            if lift > 1.5 and p_cond > 0.01:
                strong_couplings.append((r1, r2, p_cond, p_base, lift, n_r1))
                marker = " ★★" if lift > 3 else " ★"
                print(f"  {r1:2d}→{r2:2d} | {p_cond:13.4f} | {p_base:13.4f} | {lift:6.1f} | {n_r1:8d}{marker}")

    # KEY: Is the coupling through SCHEDULE or STATE?
    # Test: condition on W[0] bin AND carry[r1]=0
    print(f"\n--- Distinguishing schedule vs state coupling ---")
    W0_vals = np.array([int(rng.randint(0, 1 << 32)) for _ in range(N)])
    # Actually need the original W0 values — recompute
    rng2 = np.random.RandomState(seed)
    W0_actual = np.array([int(rng2.randint(0, 1 << 32)) for _ in range(N)])

    # Bin W0 into quartiles
    q25, q75 = np.percentile(W0_actual, [25, 75])

    for r1, r2, p_cond, p_base, lift, n_r1 in strong_couplings[:5]:
        # Within low-W0 quartile
        lo_mask = W0_actual < q25
        r1_and_lo = (all_carries[:, r1] == 0) & lo_mask
        n_joint = np.sum(r1_and_lo)
        if n_joint > 5:
            p_cond_lo = np.mean(1 - all_carries[r1_and_lo, r2])
            # Just r1=0 in lo quartile baseline
            r1_in_lo = np.sum(all_carries[lo_mask, r1] == 0)
            p_base_lo = np.mean(1 - all_carries[lo_mask, r2])
            lift_lo = p_cond_lo / max(p_base_lo, 1e-5)
            print(f"  {r1}→{r2}: full lift={lift:.1f}×, W0_lo lift={lift_lo:.1f}× (N={n_joint})")
            if lift_lo > 1.3:
                print(f"    → STATE coupling persists within W0 bin!")
            else:
                print(f"    → Coupling explained by W0 (schedule only)")

    return strong_couplings


# ============================================================
# U2: A-Branch Bridge — Maj → H[0..3]
# ============================================================

def experiment_U2(N=15000, seed=1001):
    print("\n" + "=" * 70)
    print("U2: A-Branch Bridge — Does Maj create signal in H[0..3]?")
    print(f"N={N}, seed={seed}")
    print("=" * 70)

    rng = np.random.RandomState(seed)

    # H[0] = a63 + IV[0], where a63 = T1[63] + T2[63]
    # T2[63] = Sig0(a[63_input]) + Maj(a[63_input], b[63_input], c[63_input])
    # a[63_input] = a after round 62 = states[63][0]
    # b[63_input] = states[63][1] = a after round 61 = states[62][0]
    # c[63_input] = states[63][2] = a after round 60 = states[61][0]

    # H[1] = b63 + IV[1] = a62 + IV[1] = states[63][0] + IV[1]
    # H[2] = c63 + IV[2] = b62 + IV[2] = a61 + IV[2]
    # H[3] = d63 + IV[3] = c62 + IV[3] = a60 + IV[3]

    # So from H we can extract:
    # a62 = H[1] - IV[1]
    # a61 = H[2] - IV[2]
    # a60 = H[3] - IV[3]
    # Maj(a62, a61, a60) = Maj(H[1]-IV[1], H[2]-IV[2], H[3]-IV[3])

    maj_vals = []
    sig0_vals = []
    h0_vals = []; h1_vals = []; h2_vals = []; h3_vals = []
    raw63_vals = []
    t2_63_vals = []

    t0 = time.time()
    for i in range(N):
        W16 = [int(rng.randint(0, 1 << 32)) for _ in range(16)]
        carries, raws, H_out, states, t2_s0, t2_m, t2_v, W = full_sha_detailed(W16)

        # Maj at round 63 input = Maj(states[63][0], states[63][1], states[63][2])
        a_in = states[63][0]
        b_in = states[63][1]
        c_in = states[63][2]
        maj63 = Maj(a_in, b_in, c_in)
        sig0_63 = Sig0(a_in)

        maj_vals.append(maj63)
        sig0_vals.append(sig0_63)
        t2_63_vals.append(t2_v[63])
        h0_vals.append(H_out[0])
        h1_vals.append(H_out[1])
        h2_vals.append(H_out[2])
        h3_vals.append(H_out[3])
        raw63_vals.append(raws[63])

    print(f"Collected: {time.time()-t0:.1f}s")

    # Can we compute Maj63 from output?
    # Maj63 = Maj(a62, b62, c62) where a62=H[1]-IV[1], b62=a61=H[2]-IV[2], c62=a60=H[3]-IV[3]
    maj_from_H = []
    for i in range(N):
        a62 = (h1_vals[i] - H0[1]) & MASK
        b62 = (h2_vals[i] - H0[2]) & MASK
        c62 = (h3_vals[i] - H0[3]) & MASK
        maj_from_H.append(Maj(a62, b62, c62))

    # Verify
    match = sum(1 for i in range(min(1000, N)) if maj_from_H[i] == maj_vals[i])
    # Actually this computes Maj at round 63 input using states[63] = (a62, b62=a61, c62=a60)
    # From H: H[1]-IV[1]=b63=a62, H[2]-IV[2]=c63=b62=a61, H[3]-IV[3]=d63=c62=a60... wait
    # H[1] = b63 + IV[1], b63 = states[64][1] = a after round 62 = states[63][0]
    # Hmm, states[64][1] = b at end = a from prev round = states[63][0]
    # So H[1]-IV[1] = states[63][0] = a_input_for_round_63 ✓
    # H[2]-IV[2] = states[64][2] = c at end = b from prev = states[63][1] ✓
    # H[3]-IV[3] = states[64][3] = d at end = c from prev = states[63][2] ✓

    print(f"\n  Verification: Maj63 from H matches internal: {match}/1000")

    # Correlation Maj63 → H[0]
    maj_arr = np.array(maj_vals, dtype=float)
    h0_arr = np.array(h0_vals, dtype=float)
    h1_arr = np.array(h1_vals, dtype=float)

    corr_maj_h0 = np.corrcoef(maj_arr, h0_arr)[0, 1]
    corr_maj_h1 = np.corrcoef(maj_arr, h1_arr)[0, 1]

    print(f"\n  --- Maj63 → H correlations ---")
    print(f"  corr(Maj63, H[0]): {corr_maj_h0:+.5f}")
    print(f"  corr(Maj63, H[1]): {corr_maj_h1:+.5f}")

    # T2 → H[0]
    t2_arr = np.array(t2_63_vals, dtype=float)
    sig0_arr = np.array(sig0_vals, dtype=float)
    corr_t2_h0 = np.corrcoef(t2_arr, h0_arr)[0, 1]
    corr_s0_h0 = np.corrcoef(sig0_arr, h0_arr)[0, 1]
    print(f"  corr(T2[63], H[0]):   {corr_t2_h0:+.5f}")
    print(f"  corr(Sig0[63], H[0]): {corr_s0_h0:+.5f}")

    # Compare with Ch → H[7] (the known bridge)
    raw63 = np.array(raw63_vals, dtype=float)
    h7_arr = np.array([0]*N, dtype=float)  # need H[7]
    # Recompute to get H[7]
    rng3 = np.random.RandomState(seed)
    h7_list = []
    for i in range(N):
        W16 = [int(rng3.randint(0, 1 << 32)) for _ in range(16)]
        _, _, H_out, _, _, _, _, _ = full_sha_detailed(W16)
        h7_list.append(H_out[7])
    h7_arr = np.array(h7_list, dtype=float)
    corr_r63_h7 = np.corrcoef(raw63, h7_arr)[0, 1]

    print(f"\n  --- Comparison: a-branch vs e-branch ---")
    print(f"  e-branch: corr(raw63, H[7]) = {corr_r63_h7:+.5f}")
    print(f"  a-branch: corr(Maj63, H[0]) = {corr_maj_h0:+.5f}")
    print(f"  a-branch: corr(T2[63], H[0])= {corr_t2_h0:+.5f}")

    if abs(corr_maj_h0) > 0.05:
        print(f"\n  ★★ MAJ BRIDGE EXISTS! a-branch has signal in H[0]!")
    elif abs(corr_t2_h0) > 0.05:
        print(f"\n  ★★ T2 BRIDGE EXISTS! a-branch has signal through T2!")
    else:
        print(f"\n  No a-branch bridge detected (corr < 0.05)")

    # Tautological baseline for Maj
    print(f"\n  --- Tautological baseline: corr(Maj(a,b,c), a+b+c) ---")
    rng4 = np.random.RandomState(seed + 200)
    a_rand = rng4.randint(0, 1<<32, 50000).astype(np.int64)
    b_rand = rng4.randint(0, 1<<32, 50000).astype(np.int64)
    c_rand = rng4.randint(0, 1<<32, 50000).astype(np.int64)
    maj_rand = np.array([Maj(int(a), int(b), int(c)) for a, b, c in zip(a_rand, b_rand, c_rand)], dtype=float)

    # H[0] = a63 + IV = T1_63 + T2_63 + IV
    # Maj is part of T2. So corr(Maj, H[0]) includes corr(Maj, T2)×corr(T2, H[0])
    # T2 = Sig0 + Maj → corr(Maj, T2) should be ~0.5-0.7
    # Is this tautological? Maj → T2 → a63 → H[0]. Three steps.

    # For random: corr(Maj(a,b,c), T1+Sig0(a)+Maj(a,b,c)+IV) where T1 independent
    t1_rand = rng4.randint(0, 1<<32, 50000).astype(np.int64)
    sig0_rand = np.array([Sig0(int(a)) for a in a_rand], dtype=float)
    h0_rand = np.array([(int(t1) + int(s0) + int(m) + H0[0]) & MASK for t1, s0, m in zip(t1_rand, sig0_rand, maj_rand)], dtype=float)
    corr_taut_maj = np.corrcoef(maj_rand, h0_rand)[0, 1]

    print(f"  Tautological corr(Maj, T1+Sig0+Maj+IV): {corr_taut_maj:+.5f}")
    print(f"  SHA-256 corr(Maj63, H[0]):               {corr_maj_h0:+.5f}")
    print(f"  Delta (excess):                           {corr_maj_h0 - corr_taut_maj:+.5f}")

    return corr_maj_h0, corr_t2_h0


# ============================================================
# U3: Full 8-word output structure
# ============================================================

def experiment_U3(N=10000, seed=1002):
    print("\n" + "=" * 70)
    print("U3: Full 8-Word Output Structure")
    print(f"N={N}, seed={seed}")
    print("=" * 70)

    rng = np.random.RandomState(seed)

    all_H = np.zeros((N, 8), dtype=np.float64)
    all_raw63 = np.zeros(N)

    t0 = time.time()
    for i in range(N):
        W16 = [int(rng.randint(0, 1 << 32)) for _ in range(16)]
        carries, raws, H_out, _, _, _, _, _ = full_sha_detailed(W16)
        all_H[i] = H_out
        all_raw63[i] = raws[63]

    print(f"Collected: {time.time()-t0:.1f}s")

    # 8×8 H correlation matrix (unconditional)
    print(f"\n  --- H[0..7] correlation matrix (unconditional) ---")
    corr_H = np.corrcoef(all_H.T)
    print(f"  {'':>5}", end='')
    for j in range(8):
        print(f" H[{j}]  ", end='')
    print()
    for i in range(8):
        print(f"  H[{i}]", end='')
        for j in range(8):
            c = corr_H[i, j]
            marker = '*' if abs(c) > 0.01 and i != j else ' '
            print(f" {c:+5.3f}{marker}", end='')
        print()

    # Correlation with raw63 (soft bridge)
    print(f"\n  --- corr(raw63, H[i]) for each output word ---")
    for i in range(8):
        c = np.corrcoef(all_raw63, all_H[:, i])[0, 1]
        marker = " ★★" if abs(c) > 0.05 else (" ★" if abs(c) > 0.02 else "")
        print(f"    corr(raw63, H[{i}]): {c:+.5f}{marker}")

    # Tail analysis: bottom 5% of raw63
    q5 = np.percentile(all_raw63, 5)
    mask = all_raw63 < q5
    n_tail = np.sum(mask)

    print(f"\n  --- Tail (raw63 < 5th pct, N={n_tail}): H correlation matrix ---")
    corr_tail = np.corrcoef(all_H[mask].T)
    print(f"  {'':>5}", end='')
    for j in range(8):
        print(f" H[{j}]  ", end='')
    print()
    for i in range(8):
        print(f"  H[{i}]", end='')
        for j in range(8):
            c = corr_tail[i, j]
            c_uncond = corr_H[i, j]
            marker = '*' if abs(c) > 0.05 and i != j else ' '
            print(f" {c:+5.3f}{marker}", end='')
        print()

    # Which H pairs gain the most correlation in the tail?
    print(f"\n  --- Largest correlation INCREASES in tail ---")
    gains = []
    for i in range(8):
        for j in range(i+1, 8):
            delta = corr_tail[i, j] - corr_H[i, j]
            gains.append((i, j, corr_H[i, j], corr_tail[i, j], delta))

    gains.sort(key=lambda x: -abs(x[4]))
    for i, j, c_unc, c_tail, delta in gains[:10]:
        marker = " ★★" if abs(delta) > 0.10 else (" ★" if abs(delta) > 0.05 else "")
        print(f"    H[{i}]↔H[{j}]: uncond={c_unc:+.4f}, tail={c_tail:+.4f}, Δ={delta:+.4f}{marker}")

    # Per-bit sensitivity of each H word to raw63
    print(f"\n  --- Per-word bit sensitivity to raw63 (top bits only) ---")
    for word in range(8):
        best_c = 0
        best_bit = -1
        for bit in range(28, 32):
            h_bit = np.array([(int(all_H[i, word]) >> bit) & 1 for i in range(N)], dtype=float)
            c = np.corrcoef(all_raw63, h_bit)[0, 1] if np.std(h_bit) > 0.01 else 0
            if abs(c) > abs(best_c):
                best_c = c
                best_bit = bit
        marker = " ★" if abs(best_c) > 0.03 else ""
        print(f"    H[{word}]: best bit {best_bit}, corr={best_c:+.4f}{marker}")

    return all_H, corr_H, corr_tail


# ============================================================
# MAIN
# ============================================================

if __name__ == '__main__':
    print("State Coupling + Unexplored Angles — Stage 3.3/3.4")
    print(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print()

    N1, N2, N3 = 12000, 12000, 10000
    if len(sys.argv) > 1 and sys.argv[1] == '--fast':
        N1, N2, N3 = 6000, 6000, 5000

    t_start = time.time()
    couplings = experiment_U1(N=N1)
    corr_maj, corr_t2 = experiment_U2(N=N2)
    all_H, corr_unc, corr_tail = experiment_U3(N=N3)

    total = time.time() - t_start
    print(f"\n{'='*70}")
    print(f"TOTAL TIME: {total:.1f}s")
    print(f"{'='*70}")

    print(f"""
{'='*70}
SYNTHESIS: State Coupling + A-Branch (Stage 3.3/3.4)
{'='*70}

U1: STATE COUPLING
  Strong couplings found: {len(couplings)}
  These are carry[r1]=0 → P(carry[r2]=0) elevated.

U2: A-BRANCH BRIDGE (Maj → H[0])
  corr(Maj63, H[0]) = {corr_maj:+.4f}
  corr(T2[63], H[0]) = {corr_t2:+.4f}
  If > 0.05 → parallel bridge through a-register!

U3: FULL 8-WORD STRUCTURE
  Do ALL output words correlate in the tail?
  If yes → multi-word distinguisher possible.
""")
