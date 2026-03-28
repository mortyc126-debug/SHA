#!/usr/bin/env python3
"""
K-Cascade Meeting Point — Stage 1.5
=====================================

THE KEY EXPERIMENT: Does carry[63]=0 create a backward cascade?

If carry[63]=0 → P(carry[62]=0) elevated → P(carry[61]=0) elevated → ...
then K-chain creates a CORRIDOR from output back toward structural zone.

Three tests:
  C1: Backward cascade — P(carry[r]=0 | carry[63]=0) for all r
  C2: K-resonance map — which rounds amplify the backward signal?
  C3: Meeting point — where does backward cascade meet forward structure?
"""

import numpy as np
from collections import defaultdict, Counter
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

def sha256_carries(W16):
    """Returns carry profile (64 ints) and raw values."""
    W = message_schedule(W16)
    a, b, c, d, e, f, g, h = H0
    carries = []
    raws = []
    for r in range(64):
        raw = h + Sig1(e) + Ch(e, f, g) + K[r] + W[r]
        T1 = raw & MASK
        carries.append(1 if raw >= (1 << 32) else 0)
        raws.append(raw)
        T2 = (Sig0(a) + Maj(a, b, c)) & MASK
        h, g, f = g, f, e
        e = (d + T1) & MASK
        d, c, b = c, b, a
        a = (T1 + T2) & MASK
    return carries, raws


def hc_find_carry0(rng, target_r=63, steps=200):
    """Hill climb to minimize raw[target_r], return (W16, carries, raws)."""
    W0 = int(rng.randint(0, 1 << 32))
    W16 = [W0] + [0] * 15
    carries, raws = sha256_carries(W16)
    best = raws[target_r]

    for _ in range(steps):
        b = int(rng.randint(0, 32))
        W16_try = [W0 ^ (1 << b)] + [0] * 15
        c2, r2 = sha256_carries(W16_try)
        if r2[target_r] < best:
            W0 ^= (1 << b)
            W16 = W16_try
            carries, raws = c2, r2
            best = r2[target_r]

    return W16, carries, raws


# ============================================================
# C1: Backward Cascade — conditional carry profile
# ============================================================

def experiment_C1(N=5000, seed=80):
    """
    Collect many carry[63]=0 events via HC.
    Measure P(carry[r]=0 | carry[63]=0) for all r.
    Compare with baseline P(carry[r]=0).
    """
    print("=" * 70)
    print("C1: Backward Cascade — P(carry[r]=0 | carry[63]=0)")
    print(f"N={N} HC attempts, seed={seed}")
    print("=" * 70)

    rng = np.random.RandomState(seed)

    # Phase 1: collect carry[63]=0 events
    c0_profiles = []  # carry profiles where carry[63]=0
    all_profiles = []

    t0 = time.time()
    for i in range(N):
        W16, carries, raws = hc_find_carry0(rng, target_r=63, steps=200)
        all_profiles.append(carries)
        if carries[63] == 0:
            c0_profiles.append(carries)

        if (i + 1) % 1000 == 0:
            print(f"  {i+1}/{N}, carry[63]=0: {len(c0_profiles)} ({time.time()-t0:.1f}s)")

    elapsed = time.time() - t0
    n_c0 = len(c0_profiles)
    print(f"\nCompleted: {elapsed:.1f}s")
    print(f"  carry[63]=0 found: {n_c0} / {N} ({n_c0/N*100:.1f}%)")

    if n_c0 < 10:
        print("  TOO FEW carry=0 events. Increase N or HC steps.")
        return None, None, None

    # Phase 2: baseline from random (no HC)
    rng2 = np.random.RandomState(seed + 500)
    baseline = np.zeros(64)
    n_base = 5000
    for _ in range(n_base):
        W16 = [int(rng2.randint(0, 1 << 32)) for _ in range(16)]
        c, _ = sha256_carries(W16)
        for r in range(64):
            baseline[r] += (1 - c[r])
    baseline /= n_base

    # Phase 3: conditional profile
    cond = np.zeros(64)
    for cp in c0_profiles:
        for r in range(64):
            cond[r] += (1 - cp[r])
    cond /= n_c0

    # Phase 4: HC-all profile (HC but not conditioned on carry[63])
    hc_all = np.zeros(64)
    for cp in all_profiles:
        for r in range(64):
            hc_all[r] += (1 - cp[r])
    hc_all /= len(all_profiles)

    print(f"\n  {'r':>3} | {'P(c=0) base':>12} | {'P(c=0) HC':>10} | {'P(c=0|c63=0)':>14} | {'ratio':>6} | {'cascade?':>9}")
    print(f"  {'-'*3}-+-{'-'*12}-+-{'-'*10}-+-{'-'*14}-+-{'-'*6}-+-{'-'*9}")

    cascade_rounds = []
    for r in range(64):
        ratio = cond[r] / max(baseline[r], 1e-10)
        is_cascade = ratio > 2.0 and cond[r] > 0.01
        if is_cascade:
            cascade_rounds.append(r)
        marker = " ★★★" if ratio > 10 else (" ★★" if ratio > 5 else (" ★" if ratio > 2 else ""))

        if baseline[r] > 0.001 or cond[r] > 0.01 or r >= 55 or r <= 5:
            print(f"  {r:3d} | {baseline[r]:12.5f} | {hc_all[r]:10.5f} | {cond[r]:14.4f} | {ratio:6.1f} | {'YES' if is_cascade else 'no':>9}{marker}")

    print(f"\n  Cascade rounds (ratio > 2×, P > 0.01): {cascade_rounds}")

    # Phase 5: Is there a CHAIN? carry[63]=0 → carry[62]=0 → carry[61]=0 ...
    print(f"\n--- Chain Analysis ---")
    # Count how many carry=0 events we get in sequence from r=63 backward
    chain_lengths = []
    for cp in c0_profiles:
        length = 0
        for r in range(63, -1, -1):
            if cp[r] == 0:
                length += 1
            else:
                break
        chain_lengths.append(length)

    print(f"  Chain from r=63 backward (consecutive carry=0):")
    print(f"    Mean length: {np.mean(chain_lengths):.2f}")
    print(f"    Max length: {max(chain_lengths)}")
    cl_hist = Counter(chain_lengths)
    for length in sorted(cl_hist.keys()):
        pct = cl_hist[length] / n_c0 * 100
        bar = '#' * int(pct)
        print(f"    length {length:2d}: {cl_hist[length]:4d} ({pct:5.1f}%) {bar}")

    # Phase 6: Non-consecutive carry=0 count in r=47..63
    print(f"\n--- Carry=0 count in r=47..63 (W5 window + tail) ---")
    w5_counts = []
    for cp in c0_profiles:
        cnt = sum(1 - cp[r] for r in range(47, 64))
        w5_counts.append(cnt)

    w5_base = []
    for cp in all_profiles[:2000]:
        cnt = sum(1 - cp[r] for r in range(47, 64))
        w5_base.append(cnt)

    print(f"  Conditioned (c63=0): mean={np.mean(w5_counts):.2f}, max={max(w5_counts)}")
    print(f"  HC baseline:         mean={np.mean(w5_base):.2f}, max={max(w5_base)}")
    print(f"  Random baseline:     ~{17 * np.mean(baseline[47:64]):.2f}")

    return cond, baseline, cascade_rounds


# ============================================================
# C2: K-Resonance Map
# ============================================================

def experiment_C2(N=3000, seed=81):
    """
    Which rounds AMPLIFY the backward signal?

    For each target round t=62,61,...,47:
      HC minimize raw[t]
      Measure P(carry[63]=0 | carry[t]=0)

    If P(c63=0|ct=0) >> baseline → round t resonates with 63.
    """
    print("\n" + "=" * 70)
    print("C2: K-Resonance Map — Which rounds amplify backward signal?")
    print(f"N={N} per target, seed={seed}")
    print("=" * 70)

    rng = np.random.RandomState(seed)

    # Baseline P(carry[r]=0) from random
    rng_base = np.random.RandomState(seed + 600)
    baseline = np.zeros(64)
    n_base = 3000
    for _ in range(n_base):
        W16 = [int(rng_base.randint(0, 1 << 32)) for _ in range(16)]
        c, _ = sha256_carries(W16)
        for r in range(64):
            baseline[r] += (1 - c[r])
    baseline /= n_base

    # For each target round, HC to find carry[t]=0, then check carry[63]
    targets = list(range(46, 64))  # r=46..63

    print(f"\n  {'target':>7} | {'c_t=0 found':>12} | {'P(c63=0|ct=0)':>14} | {'P(c63=0) base':>14} | {'lift':>6}")
    print(f"  {'-'*7}-+-{'-'*12}-+-{'-'*14}-+-{'-'*14}-+-{'-'*6}")

    resonance = {}
    n_per_target = N // len(targets)

    t0 = time.time()
    for t in targets:
        n_ct0 = 0
        n_c63_given_ct0 = 0

        for _ in range(n_per_target):
            W16, carries, raws = hc_find_carry0(rng, target_r=t, steps=150)
            if carries[t] == 0:
                n_ct0 += 1
                if carries[63] == 0:
                    n_c63_given_ct0 += 1

        p_c63_given_ct0 = n_c63_given_ct0 / max(n_ct0, 1)
        p_c63_base = baseline[63]
        lift = p_c63_given_ct0 / max(p_c63_base, 1e-10)

        resonance[t] = (n_ct0, p_c63_given_ct0, lift)
        marker = " ★★★" if lift > 10 else (" ★★" if lift > 3 else (" ★" if lift > 1.5 else ""))
        print(f"  r={t:3d} | {n_ct0:12d} | {p_c63_given_ct0:14.4f} | {p_c63_base:14.5f} | {lift:6.1f}{marker}")

    print(f"\n  Completed: {time.time()-t0:.1f}s")

    return resonance


# ============================================================
# C3: Meeting Point — forward structure meets backward cascade
# ============================================================

def experiment_C3(N=3000, seed=82):
    """
    Where does the backward cascade meet the forward structure?

    Forward: W[0] → schedule → W[16..63]. Structure dies at r≈31.
    Backward: carry[63]=0 → constraints on state[62] → ... → ?

    Test: for carry[63]=0 events, measure how much of the carry profile
    in r=16..46 is predictable from W[0] alone.

    If the backward cascade reaches r≈46 (W5 window start),
    and forward schedule structure reaches r≈30, there's a gap of 16 rounds.
    If the gap shrinks — we found the meeting point.
    """
    print("\n" + "=" * 70)
    print("C3: Meeting Point — Forward Structure Meets Backward Cascade")
    print(f"N={N}, seed={seed}")
    print("=" * 70)

    rng = np.random.RandomState(seed)

    # Collect carry[63]=0 events with full data
    c0_data = []  # (W0, carry_profile, raw_profile)
    t0 = time.time()

    for i in range(N):
        W16, carries, raws = hc_find_carry0(rng, target_r=63, steps=200)
        if carries[63] == 0:
            c0_data.append((W16[0], carries, raws))

        if (i + 1) % 1000 == 0:
            print(f"  {i+1}/{N}, c63=0: {len(c0_data)} ({time.time()-t0:.1f}s)")

    print(f"\n  carry[63]=0 found: {len(c0_data)}")

    if len(c0_data) < 10:
        print("  Too few events. Showing what we have.")
        if len(c0_data) == 0:
            return [], {}
        # Continue with what we have

    # For each carry[63]=0 event, compute:
    # 1. Total carry=0 count per "zone"
    # 2. Correlation between zones

    zones = {
        'forward_early': list(range(0, 16)),    # W[0..15] direct
        'forward_late':  list(range(16, 31)),    # schedule nearing wall
        'chaos':         list(range(31, 47)),    # chaotic zone
        'W5_window':     list(range(47, 56)),    # carry window 5
        'tail':          list(range(56, 63)),    # approach to r=63
        'target':        [63],                   # carry[63]=0 (always 0 here)
    }

    zone_scores = {z: [] for z in zones}
    for W0, cp, raws in c0_data:
        for zname, rounds in zones.items():
            score = sum(1 - cp[r] for r in rounds)
            zone_scores[zname].append(score)

    print(f"\n--- Zone Analysis (carry[63]=0 events only) ---")
    print(f"  {'Zone':>16} | {'Rounds':>12} | {'Mean c=0':>9} | {'Max c=0':>8} | {'c=0/total':>10}")
    print(f"  {'-'*16}-+-{'-'*12}-+-{'-'*9}-+-{'-'*8}-+-{'-'*10}")
    for zname in ['forward_early', 'forward_late', 'chaos', 'W5_window', 'tail', 'target']:
        rounds = zones[zname]
        scores = zone_scores[zname]
        print(f"  {zname:>16} | {f'r={min(rounds)}..{max(rounds)}':>12} | {np.mean(scores):9.2f} | {max(scores):8d} | {np.mean(scores)/len(rounds):10.3f}")

    # Correlation between zones
    print(f"\n--- Zone Correlations (carry[63]=0 events) ---")
    zone_names = ['forward_early', 'forward_late', 'chaos', 'W5_window', 'tail']
    for i, z1 in enumerate(zone_names):
        for z2 in zone_names[i+1:]:
            a1 = np.array(zone_scores[z1], dtype=float)
            a2 = np.array(zone_scores[z2], dtype=float)
            if np.std(a1) > 0 and np.std(a2) > 0:
                c = np.corrcoef(a1, a2)[0, 1]
                marker = " ★" if abs(c) > 0.1 else ""
                print(f"  corr({z1}, {z2}) = {c:+.4f}{marker}")

    # KEY TEST: Does W[0] predict carry in tail zone (r=56..63)?
    print(f"\n--- KEY TEST: Does W[0] predict tail carry? ---")
    W0_bits = np.array([[((d[0] >> b) & 1) for b in range(32)] for d in c0_data], dtype=float)
    tail_scores_arr = np.array(zone_scores['tail'], dtype=float)

    best_bit_corr = 0
    best_bit = -1
    for b in range(32):
        if np.std(W0_bits[:, b]) > 0 and np.std(tail_scores_arr) > 0:
            c = np.corrcoef(W0_bits[:, b], tail_scores_arr)[0, 1]
            if abs(c) > abs(best_bit_corr):
                best_bit_corr = c
                best_bit = b

    print(f"  Best corr(W0[bit], tail_score): {best_bit_corr:+.4f} at bit {best_bit}")

    # W[0] → chaos zone
    chaos_arr = np.array(zone_scores['chaos'], dtype=float)
    best_chaos_corr = 0
    for b in range(32):
        if np.std(W0_bits[:, b]) > 0 and np.std(chaos_arr) > 0:
            c = np.corrcoef(W0_bits[:, b], chaos_arr)[0, 1]
            if abs(c) > abs(best_chaos_corr):
                best_chaos_corr = c

    print(f"  Best corr(W0[bit], chaos_score): {best_chaos_corr:+.4f}")

    # W[0] → forward_late
    fl_arr = np.array(zone_scores['forward_late'], dtype=float)
    best_fl_corr = 0
    for b in range(32):
        if np.std(W0_bits[:, b]) > 0 and np.std(fl_arr) > 0:
            c = np.corrcoef(W0_bits[:, b], fl_arr)[0, 1]
            if abs(c) > abs(best_fl_corr):
                best_fl_corr = c

    print(f"  Best corr(W0[bit], forward_late_score): {best_fl_corr:+.4f}")

    # Carry profile heatmap: per-round P(carry=0) conditioned on carry[63]=0
    print(f"\n--- Per-round P(carry=0 | carry[63]=0) for c63=0 events ---")
    cond_p = np.zeros(64)
    for _, cp, _ in c0_data:
        for r in range(64):
            cond_p[r] += (1 - cp[r])
    cond_p /= len(c0_data)

    # Visual map
    print(f"  Carry=0 probability map (conditioned on c63=0):")
    print(f"  {'':>5}", end='')
    for r in range(64):
        if r % 8 == 0:
            print(f" |{r:2d}", end='')
    print()
    print(f"  {'':>5}", end='')
    for r in range(64):
        p = cond_p[r]
        if p > 0.10:
            ch = '█'
        elif p > 0.05:
            ch = '▓'
        elif p > 0.02:
            ch = '▒'
        elif p > 0.005:
            ch = '░'
        else:
            ch = '·'
        if r % 8 == 0:
            print(f" |{ch}", end='')
        else:
            print(ch, end='')
    print()
    print(f"  Legend: █>10% ▓>5% ▒>2% ░>0.5% ·<0.5%")

    return c0_data, zone_scores


# ============================================================
# MAIN
# ============================================================

if __name__ == '__main__':
    print("K-Cascade Meeting Point — Stage 1.5")
    print(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print()

    N1, N2, N3 = 8000, 3000, 6000
    if len(sys.argv) > 1 and sys.argv[1] == '--fast':
        N1, N2, N3 = 4000, 1500, 3000

    t_start = time.time()

    cond, baseline, cascade_rounds = experiment_C1(N=N1, seed=80)
    resonance = experiment_C2(N=N2, seed=81)
    c0_data, zones = experiment_C3(N=N3, seed=82)

    total = time.time() - t_start
    print(f"\n{'='*70}")
    print(f"TOTAL TIME: {total:.1f}s")
    print(f"{'='*70}")

    print(f"""
{'='*70}
SYNTHESIS: K-Cascade Meeting Point
{'='*70}

QUESTION: Does carry[63]=0 create a backward cascade?

If backward cascade reaches r≈46 (W5 start):
  → K-chain corridor: 17 rounds (r=46..63)
  → Forward structure: 30 rounds (r=0..30)
  → Gap: only r=31..46 (15 rounds of chaos)
  → This is SMALLER than the full 29-round chaotic zone

If backward cascade reaches r≈30:
  → MEETING POINT FOUND
  → Forward structure (r=0..30) overlaps backward cascade (r=30..63)
  → Carry-algebra has an algebraic corridor through chaos

Cascade rounds detected: {cascade_rounds if cascade_rounds else 'checking...'}
""")
