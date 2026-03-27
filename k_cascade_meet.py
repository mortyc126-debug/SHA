#!/usr/bin/env python3
"""
K-Cascade Meet — Stage 1.5: Does K-chain cascade backward to meet schedule?
=============================================================================

From Stage 1.4:
  - ALL 64 rounds have K-constraints on SS when carry=0
  - 96 bits of state[62] directly from H
  - corr(Ch62, SS63) = 0.575

THE EXPERIMENT: Use HC to get carry[63]=0, then measure:
  1. P(carry[r]=0 | carry[63]=0) for r=62..0 — backward cascade
  2. If carry[r]=0 cascades back to r~55, we have 8+ rounds of constraints
  3. Forward: schedule sensitivity to W[0] at each round
  4. MEETING POINT: round where backward K-chain AND forward schedule overlap

If meeting point exists → path around chaos without going through it.
"""

import numpy as np
from collections import defaultdict
import time
import sys

MASK = 0xFFFFFFFF

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

def sha256_full_trace(W16):
    W = message_schedule(W16)
    a, b, c, d, e, f, g, h = H0
    states = [(a, b, c, d, e, f, g, h)]
    carries = []
    raws = []
    ss_vals = []  # SS[r] = h + Sig1(e) + Ch(e,f,g) before adding K[r]+W[r]

    for r in range(64):
        ss = h + Sig1(e) + Ch(e, f, g)
        ss_vals.append(ss)
        raw = ss + K[r] + W[r]
        T1 = raw & MASK
        carries.append(1 if raw >= (1 << 32) else 0)
        raws.append(raw)
        T2 = (Sig0(a) + Maj(a, b, c)) & MASK
        h, g, f = g, f, e
        e = (d + T1) & MASK
        d, c, b = c, b, a
        a = (T1 + T2) & MASK
        states.append((a, b, c, d, e, f, g, h))

    H_out = [(states[-1][i] + H0[i]) & MASK for i in range(8)]
    return states, carries, raws, ss_vals, W, H_out


def hc_find_carry0(rng, target_r=63, steps=200):
    """HC to minimize raw[target_r], returns trace with carry[target_r]=0 or best attempt."""
    W0 = int(rng.randint(0, 1 << 32))
    W16 = [W0] + [0] * 15
    states, carries, raws, ss_vals, W, H = sha256_full_trace(W16)
    best_raw = raws[target_r]

    for _ in range(steps):
        b = int(rng.randint(0, 32))
        W0_try = W0 ^ (1 << b)
        W16_try = [W0_try] + [0] * 15
        st2, ca2, ra2, ss2, W2, H2 = sha256_full_trace(W16_try)
        if ra2[target_r] < best_raw:
            W0 = W0_try
            W16, states, carries, raws, ss_vals, W, H = W16_try, st2, ca2, ra2, ss2, W2, H2
            best_raw = ra2[target_r]

    return W16, states, carries, raws, ss_vals, W, H


# ============================================================
# EXPERIMENT M1: Backward Cascade from carry[63]=0
# ============================================================

def exp_backward_cascade(N_c0_target=500, seed=80):
    """
    Collect N_c0_target instances of carry[63]=0 via HC.
    For each, record full carry profile.
    Measure P(carry[r]=0 | carry[63]=0) for all r.
    Compare with baseline P(carry[r]=0).
    """
    print("=" * 70)
    print(f"M1: Backward Cascade — P(carry[r]=0 | carry[63]=0)")
    print(f"Target: {N_c0_target} carry[63]=0 events")
    print("=" * 70)

    rng = np.random.RandomState(seed)

    # Phase 1: collect carry[63]=0 via HC
    c0_carries = []  # carries when carry[63]=0
    c0_ss = []       # SS values when carry[63]=0
    n_attempts = 0
    t0 = time.time()

    while len(c0_carries) < N_c0_target:
        W16, states, carries, raws, ss_vals, W, H = hc_find_carry0(rng, target_r=63, steps=200)
        n_attempts += 1
        if carries[63] == 0:
            c0_carries.append(carries)
            c0_ss.append(ss_vals)

        if n_attempts % 500 == 0:
            print(f"  attempts={n_attempts}, found={len(c0_carries)}, "
                  f"rate=1/{n_attempts/max(len(c0_carries),1):.0f}, {time.time()-t0:.1f}s")

    elapsed = time.time() - t0
    print(f"\nCollected {N_c0_target} carry[63]=0 in {n_attempts} attempts ({elapsed:.1f}s)")
    print(f"Rate: 1/{n_attempts/N_c0_target:.0f}")

    # Phase 2: baseline from random W
    n_baseline = 5000
    base_carry_sum = np.zeros(64)
    for _ in range(n_baseline):
        W16 = [int(rng.randint(0, 1 << 32)) for _ in range(16)]
        _, carries_b, _, _, _, _ = sha256_full_trace(W16)
        for r in range(64):
            base_carry_sum[r] += (1 - carries_b[r])
    base_p = base_carry_sum / n_baseline

    # Phase 3: P(carry[r]=0 | carry[63]=0)
    cond_p = np.zeros(64)
    for carries in c0_carries:
        for r in range(64):
            cond_p[r] += (1 - carries[r])
    cond_p /= N_c0_target

    print(f"\n{'r':>3} | {'P(c=0|c63=0)':>13} | {'P(c=0) base':>12} | {'ratio':>6} | {'lift':>6} | K det bits")
    print("-" * 75)

    # K det bits map
    k_det = {}
    for r in range(64):
        max_ss = ((1 << 32) - K[r]) & MASK
        det = 0
        for b in range(31, -1, -1):
            if max_ss < (1 << (b + 1)):
                det += 1
            else:
                break
        k_det[r] = det

    chain_depth = 0
    for r in range(63, -1, -1):
        ratio = cond_p[r] / max(base_p[r], 1e-10)
        lift = cond_p[r] - base_p[r]
        marker = ""
        if ratio > 2.0 and cond_p[r] > 0.005:
            marker = " ★★★" if ratio > 10 else (" ★★" if ratio > 5 else " ★")
            if r >= 50:
                chain_depth = max(chain_depth, 63 - r + 1)

        print(f"{r:3d} | {cond_p[r]:13.4f} | {base_p[r]:12.4f} | {ratio:6.1f} | {lift:+6.4f} | {k_det[r]}{marker}")

    # Backward chain summary
    print(f"\n--- Backward Chain Summary ---")
    chain_rounds = []
    for r in range(63, -1, -1):
        ratio = cond_p[r] / max(base_p[r], 1e-10)
        if ratio > 1.5 and cond_p[r] > 0.002:
            chain_rounds.append(r)

    print(f"  Rounds with carry[r]=0 elevated (ratio>1.5): {sorted(chain_rounds)}")
    print(f"  Chain from r=63 backward: {[r for r in sorted(chain_rounds, reverse=True) if r >= 45]}")

    # Consecutive chain from r=63
    consec = 0
    for r in range(63, -1, -1):
        ratio = cond_p[r] / max(base_p[r], 1e-10)
        if ratio > 1.3 and cond_p[r] > 0.001:
            consec += 1
        else:
            break
    print(f"  Consecutive chain length from r=63: {consec} rounds")

    return cond_p, base_p, c0_carries, c0_ss


# ============================================================
# EXPERIMENT M2: Forward Schedule Influence Map
# ============================================================

def exp_forward_schedule(N=3000, seed=81):
    """
    For each round r, measure: does W[0] influence state[r]?
    Specifically: corr(W[0], SS[r]) and sensitivity d(SS[r])/d(W[0][bit]).

    This gives the "forward signal strength" at each round.
    """
    print("\n" + "=" * 70)
    print(f"M2: Forward Schedule Influence Map")
    print(f"N={N}")
    print("=" * 70)

    rng = np.random.RandomState(seed)

    # Method: for each W[0], flip bit b, measure |ΔSS[r]|
    n_jac = 300
    sensitivity = np.zeros(64)  # fraction of (W0, bit) pairs that change SS[r]

    t0 = time.time()
    for trial in range(n_jac):
        W0 = int(rng.randint(0, 1 << 32))
        W16 = [W0] + [0] * 15
        _, _, _, ss_base, _, _ = sha256_full_trace(W16)

        for b in range(32):
            W16_flip = [W0 ^ (1 << b)] + [0] * 15
            _, _, _, ss_flip, _, _ = sha256_full_trace(W16_flip)
            for r in range(64):
                if ss_flip[r] != ss_base[r]:
                    sensitivity[r] += 1

    sensitivity /= (n_jac * 32)
    print(f"Computed: {time.time()-t0:.1f}s")

    # Also: direct correlation W[0] → SS[r]
    w0_vals = np.zeros(N)
    ss_matrix = np.zeros((64, N))

    for i in range(N):
        W0 = int(rng.randint(0, 1 << 32))
        W16 = [W0] + [0] * 15
        _, _, _, ss_vals, _, _ = sha256_full_trace(W16)
        w0_vals[i] = W0
        for r in range(64):
            ss_matrix[r, i] = ss_vals[r]

    corr_w0_ss = np.zeros(64)
    for r in range(64):
        if np.std(ss_matrix[r]) > 0:
            corr_w0_ss[r] = np.corrcoef(w0_vals, ss_matrix[r])[0, 1]

    print(f"\n{'r':>3} | {'sensitivity':>11} | {'corr(W0,SS)':>12} | {'signal':>8}")
    print("-" * 50)
    for r in range(64):
        signal = sensitivity[r] * abs(corr_w0_ss[r]) if abs(corr_w0_ss[r]) > 0.01 else 0
        marker = " ★" if signal > 0.01 else ""
        if sensitivity[r] > 0.01 or r < 5 or r > 58:
            print(f"{r:3d} | {sensitivity[r]:11.4f} | {corr_w0_ss[r]:+12.5f} | {signal:8.5f}{marker}")

    return sensitivity, corr_w0_ss


# ============================================================
# EXPERIMENT M3: Meeting Point — Overlap of forward and backward
# ============================================================

def exp_meeting_point(cond_p, base_p, sensitivity, corr_w0_ss):
    """
    THE KEY ANALYSIS: Find rounds where both:
      - Backward signal: P(carry[r]=0 | c63=0) / P(carry[r]=0) > 1.5
      - Forward signal: sensitivity(W[0] → SS[r]) > 0.01

    If such rounds exist → the K-chain connects input to output
    via a path AROUND the chaotic zone.
    """
    print("\n" + "=" * 70)
    print("M3: MEETING POINT — Forward Schedule ∩ Backward K-Chain")
    print("=" * 70)

    print(f"\n{'r':>3} | {'back ratio':>10} | {'fwd sens':>9} | {'fwd corr':>9} | {'product':>8} | status")
    print("-" * 70)

    meeting_rounds = []
    for r in range(64):
        back_ratio = cond_p[r] / max(base_p[r], 1e-10)
        fwd_sens = sensitivity[r]
        fwd_corr = abs(corr_w0_ss[r])
        product = (back_ratio - 1) * fwd_sens  # joint signal strength

        status = ""
        if back_ratio > 1.5 and fwd_sens > 0.1:
            status = "★★★ MEETING"
            meeting_rounds.append(r)
        elif back_ratio > 1.3 and fwd_sens > 0.1:
            status = "★★ NEAR"
            meeting_rounds.append(r)
        elif back_ratio > 1.5 and fwd_sens > 0.01:
            status = "★ PARTIAL"
        elif fwd_sens > 0.5:
            status = "(fwd only)"
        elif back_ratio > 2.0:
            status = "(back only)"

        if status or r < 3 or r > 58:
            print(f"{r:3d} | {back_ratio:10.2f} | {fwd_sens:9.4f} | {fwd_corr:9.5f} | {product:+8.4f} | {status}")

    print(f"\n--- MEETING POINT RESULT ---")
    if meeting_rounds:
        print(f"  Meeting rounds found: {meeting_rounds}")
        print(f"  Range: r={min(meeting_rounds)}..{max(meeting_rounds)}")
        print(f"  >>> PATH AROUND CHAOS EXISTS <<<")
        print(f"  Forward (W[0] → schedule → SS[r]): sensitive at r={min(meeting_rounds)}")
        print(f"  Backward (carry[63]=0 → K-chain → carry[r]=0): elevated at r={max(meeting_rounds)}")
        if max(meeting_rounds) - min(meeting_rounds) > 5:
            print(f"  >>> OVERLAP ZONE: {max(meeting_rounds) - min(meeting_rounds)} rounds <<<")
    else:
        print(f"  No meeting point found.")
        print(f"  Forward signal dies at r ≈ {max((r for r in range(64) if sensitivity[r] > 0.5), default=-1)}")
        back_chain = [r for r in range(64) if cond_p[r] / max(base_p[r], 1e-10) > 1.5]
        print(f"  Backward chain reaches r ≈ {min(back_chain) if back_chain else 63}")
        if back_chain:
            gap = max((r for r in range(64) if sensitivity[r] > 0.5), default=0) - min(back_chain)
            print(f"  GAP: {abs(gap)} rounds of no signal")

    return meeting_rounds


# ============================================================
# EXPERIMENT M4: Double-carry cascade depth
# ============================================================

def exp_double_carry(c0_carries, N_c0):
    """
    Among carry[63]=0 instances, how often is carry[62]=0 too?
    And carry[61]? Build the full backward conditional chain.
    """
    print("\n" + "=" * 70)
    print("M4: Double/Triple Carry Cascade")
    print(f"N={N_c0} carry[63]=0 instances")
    print("=" * 70)

    # Build conditional chain: P(carry[r]=0 | carry[r+1..63]=0)
    # Start with all N_c0 having carry[63]=0
    # Filter for carry[62]=0, then carry[61]=0, etc.

    current_set = list(range(N_c0))  # indices into c0_carries

    print(f"\n  {'condition':>30} | {'N':>6} | {'P(next=0)':>10} | {'filt/total':>11}")
    print(f"  {'-'*30}-+-{'-'*6}-+-{'-'*10}-+-{'-'*11}")
    print(f"  {'carry[63]=0':>30} | {N_c0:6d} | {'—':>10} | {N_c0}/{N_c0}")

    chain_data = [(63, N_c0)]
    for target_r in range(62, max(62 - 15, -1), -1):
        # How many of current_set also have carry[target_r]=0?
        new_set = [i for i in current_set if c0_carries[i][target_r] == 0]
        p = len(new_set) / max(len(current_set), 1)

        condition = f"carry[{target_r}..63]=0"
        print(f"  {condition:>30} | {len(new_set):6d} | {p:10.4f} | {len(new_set)}/{N_c0}")

        chain_data.append((target_r, len(new_set)))

        if len(new_set) < 3:
            print(f"  {'(chain broken — too few events)':>30}")
            break

        current_set = new_set

    # Summary
    print(f"\n  Backward chain depth: {63 - chain_data[-1][0]} rounds from r=63")
    print(f"  Deepest round with >=3 events: r={chain_data[-1][0]}")

    # Conditional P(carry[r]=0 | carry[r+1]=0) for each pair
    print(f"\n--- Per-round conditional cascade P(c[r]=0 | c[r+1]=0) ---")
    for target_r in range(62, max(62 - 12, -1), -1):
        # Among all instances with carry[target_r+1]=0, how many have carry[target_r]=0?
        n_r1_0 = sum(1 for c in c0_carries if c[target_r + 1] == 0)
        n_both = sum(1 for c in c0_carries if c[target_r + 1] == 0 and c[target_r] == 0)
        p_cond = n_both / max(n_r1_0, 1)
        print(f"  P(c[{target_r}]=0 | c[{target_r+1}]=0) = {p_cond:.4f}  (N={n_r1_0})")

    return chain_data


# ============================================================
# MAIN
# ============================================================

if __name__ == '__main__':
    print("K-Cascade Meet — Stage 1.5")
    print(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print()

    N_c0 = 300
    N_fwd = 2000

    if len(sys.argv) > 1 and sys.argv[1] == '--fast':
        N_c0 = 150
        N_fwd = 1000

    t_start = time.time()

    # M1: backward cascade
    cond_p, base_p, c0_carries, c0_ss = exp_backward_cascade(N_c0_target=N_c0)

    # M2: forward schedule
    sensitivity, corr_w0_ss = exp_forward_schedule(N=N_fwd)

    # M3: meeting point
    meeting = exp_meeting_point(cond_p, base_p, sensitivity, corr_w0_ss)

    # M4: double carry depth
    chain_data = exp_double_carry(c0_carries, N_c0)

    total = time.time() - t_start

    print(f"\n{'='*70}")
    print(f"TOTAL TIME: {total:.1f}s")
    print(f"{'='*70}")

    print(f"""
{'='*70}
SYNTHESIS: K-Cascade Meeting Point
{'='*70}

BACKWARD (from carry[63]=0):
  Chain depth: {63 - chain_data[-1][0]} rounds back from r=63
  Deepest: r={chain_data[-1][0]} with {chain_data[-1][1]} events

FORWARD (from W[0] through schedule):
  Sensitivity > 0.5 until r ≈ {max((r for r in range(64) if sensitivity[r] > 0.5), default=-1)}
  Then chaos: sensitivity stays ~1.0 but corr→0

MEETING POINTS: {meeting if meeting else 'NONE FOUND'}

INTERPRETATION:
  If meeting exists → K-chain provides a corridor around chaos
  If gap exists → the corridor is broken, need a different approach
""")
