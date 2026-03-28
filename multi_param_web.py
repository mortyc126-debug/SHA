#!/usr/bin/env python3
"""
Multi-Parameter Web — Stage 1.7: Bridging the two islands
==========================================================

Stage 1.6 found two-island structure:
  Island 1: W1+W2 (r=0-13)  — strong cascade from carry[63]=0
  Island 2: W5 (r=47-51)    — weak cascade
  GAP: r=18-46              — no signal

With W[1..15]=0, schedule is sparse. Adding active W[k] words
may fill the gap by creating carry=0 at intermediate rounds.

Experiments:
  M1: Which W[k] (k=1..15) most affects carry at gap rounds (r=18-46)?
  M2: Dual-HC — minimize SS[63] using W[0] AND W[k]
  M3: Full web — W[0]+W[k] carry profile vs single W[0]
  M4: Does the gap close? P(carry[r]=0|c63=0) for gap rounds with 2 params
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

def message_schedule(W16):
    W = list(W16) + [0] * 48
    for i in range(16, 64):
        W[i] = (sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16]) & MASK
    return W

def sha256_full_trace(W16):
    """Returns carries (64), raws (64), H_out (8)."""
    W = message_schedule(W16)
    a, b, c, d, e, f, g, h = H0
    carries = []
    raws = []
    for r in range(64):
        raw = h + Sig1(e) + Ch(e, f, g) + K[r] + W[r]
        carries.append(1 if raw >= (1 << 32) else 0)
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
# M1: Which W[k] most affects gap rounds (r=18-46)?
# ============================================================

def experiment_M1(N=300, seed=100):
    """
    For each k=1..15, measure how much W[k] affects carry at gap rounds.
    Sensitivity = P(carry[r] changes when W[k] flipped by 1 bit).
    """
    print("=" * 70)
    print("M1: W[k] Sensitivity at Gap Rounds (r=18-46)")
    print(f"N={N} per k, seed={seed}")
    print("=" * 70)

    rng = np.random.RandomState(seed)
    gap_rounds = list(range(18, 47))

    # For each k, flip each of 32 bits and measure carry change in gap
    print(f"\n  Schedule dependency: W[r] depends on W[k] if...")
    print(f"  W[16]=sig1(W[14])+W[9]+sig0(W[1])+W[0]")
    print(f"  W[17]=sig1(W[15])+W[10]+sig0(W[2])+W[1]")
    print(f"  W[r]=sig1(W[r-2])+W[r-7]+sig0(W[r-15])+W[r-16]")
    print()

    # Analytical: which W[k] enters which W[r] for r=16..63
    print(f"  Direct schedule dependencies (first hop):")
    for k in range(16):
        deps = []
        for r in range(16, 64):
            if r - 2 == k or r - 7 == k or r - 15 == k or r - 16 == k:
                deps.append(r)
        gap_deps = [r for r in deps if 18 <= r <= 46]
        if gap_deps:
            print(f"    W[{k:2d}] → {deps[:8]}{'...' if len(deps)>8 else ''} (gap: {gap_deps})")

    # Empirical sensitivity
    sensitivity = np.zeros((16, 64))  # sensitivity[k][r]

    t0 = time.time()
    for trial in range(N):
        W16_base = [int(rng.randint(0, 1 << 32)) for _ in range(16)]
        c_base, _, _ = sha256_full_trace(W16_base)

        for k in range(16):
            b = int(rng.randint(0, 32))
            W16_flip = list(W16_base)
            W16_flip[k] ^= (1 << b)
            c_flip, _, _ = sha256_full_trace(W16_flip)
            for r in range(64):
                if c_base[r] != c_flip[r]:
                    sensitivity[k][r] += 1

    sensitivity /= N
    print(f"\nComputed: {time.time()-t0:.1f}s")

    # Gap sensitivity per W[k]
    print(f"\n  {'W[k]':>5} | {'gap sens':>9} | {'early':>6} | {'late':>6} | {'best gap r':>11}")
    print(f"  {'-'*5}-+-{'-'*9}-+-{'-'*6}-+-{'-'*6}-+-{'-'*11}")

    gap_scores = []
    for k in range(16):
        gap_s = np.mean(sensitivity[k][18:47])
        early_s = np.mean(sensitivity[k][0:18])
        late_s = np.mean(sensitivity[k][47:64])
        best_gap_r = 18 + np.argmax(sensitivity[k][18:47])
        gap_scores.append((k, gap_s, best_gap_r))
        marker = " ★" if gap_s > 0.05 else ""
        print(f"  W[{k:2d}] | {gap_s:9.4f} | {early_s:6.4f} | {late_s:6.4f} | r={best_gap_r:2d}{marker}")

    gap_scores.sort(key=lambda x: -x[1])
    best_k = gap_scores[0][0]
    print(f"\n  Best W[k] for gap: W[{best_k}] (gap_sens={gap_scores[0][1]:.4f})")
    print(f"  Top-3: {[(f'W[{k}]', f'{s:.4f}') for k, s, _ in gap_scores[:3]]}")

    return sensitivity, best_k


# ============================================================
# M2: Dual-HC — minimize SS[63] using W[0] AND W[k]
# ============================================================

def experiment_M2(best_k, N=2000, seed=101):
    """
    HC with TWO parameters: W[0] and W[best_k].
    Compare carry=0 rate and profile with single-W[0] HC.
    """
    print("\n" + "=" * 70)
    print(f"M2: Dual-HC — W[0] + W[{best_k}] → minimize SS[63]")
    print(f"N={N}, seed={seed}")
    print("=" * 70)

    rng = np.random.RandomState(seed)

    # Single-param HC (W[0] only, W[1..15]=0)
    single_c0 = 0
    single_profiles = []

    # Dual-param HC (W[0] + W[best_k])
    dual_c0 = 0
    dual_profiles = []

    t0 = time.time()
    for i in range(N):
        # --- Single HC ---
        W0 = int(rng.randint(0, 1 << 32))
        W16_s = [W0] + [0] * 15
        _, raws_s, _ = sha256_full_trace(W16_s)
        best_raw = raws_s[63]

        for step in range(200):
            b = int(rng.randint(0, 32))
            W16_try = list(W16_s)
            W16_try[0] ^= (1 << b)
            _, raws_try, _ = sha256_full_trace(W16_try)
            if raws_try[63] < best_raw:
                W16_s = W16_try
                best_raw = raws_try[63]

        c_s, _, _ = sha256_full_trace(W16_s)
        single_profiles.append(c_s)
        if c_s[63] == 0:
            single_c0 += 1

        # --- Dual HC ---
        W0 = int(rng.randint(0, 1 << 32))
        Wk = int(rng.randint(0, 1 << 32))
        W16_d = [0] * 16
        W16_d[0] = W0
        W16_d[best_k] = Wk
        _, raws_d, _ = sha256_full_trace(W16_d)
        best_raw_d = raws_d[63]

        for step in range(200):
            # Alternate flipping W[0] and W[k]
            which = 0 if step % 2 == 0 else best_k
            b = int(rng.randint(0, 32))
            W16_try = list(W16_d)
            W16_try[which] ^= (1 << b)
            _, raws_try, _ = sha256_full_trace(W16_try)
            if raws_try[63] < best_raw_d:
                W16_d = W16_try
                best_raw_d = raws_try[63]

        c_d, _, _ = sha256_full_trace(W16_d)
        dual_profiles.append(c_d)
        if c_d[63] == 0:
            dual_c0 += 1

        if (i + 1) % 500 == 0:
            print(f"  {i+1}/{N}: single c0={single_c0}, dual c0={dual_c0} ({time.time()-t0:.1f}s)")

    print(f"\nCompleted: {time.time()-t0:.1f}s")
    print(f"  Single HC (W[0] only):    carry[63]=0: {single_c0}/{N} ({single_c0/N*100:.2f}%)")
    print(f"  Dual HC (W[0]+W[{best_k}]):  carry[63]=0: {dual_c0}/{N} ({dual_c0/N*100:.2f}%)")
    print(f"  Speedup: {(dual_c0/max(single_c0,1)):.2f}×")

    # Compare profiles: gap rounds
    gap = list(range(18, 47))
    single_gap = np.mean([[1 - c[r] for r in gap] for c in single_profiles], axis=0)
    dual_gap = np.mean([[1 - c[r] for r in gap] for c in dual_profiles], axis=0)

    print(f"\n  P(carry[r]=0) for gap rounds (r=18..46):")
    print(f"  {'r':>3} | {'single':>8} | {'dual':>8} | {'lift':>6}")
    print(f"  {'-'*3}-+-{'-'*8}-+-{'-'*8}-+-{'-'*6}")

    improved_rounds = []
    for i, r in enumerate(gap):
        if single_gap[i] > 0.001 or dual_gap[i] > 0.001:
            lift = dual_gap[i] / max(single_gap[i], 1e-5)
            marker = " ★" if lift > 2.0 and dual_gap[i] > 0.005 else ""
            if marker:
                improved_rounds.append(r)
            print(f"  {r:3d} | {single_gap[i]:8.4f} | {dual_gap[i]:8.4f} | {lift:6.1f}{marker}")

    # Total carry=0 count
    single_total = np.mean([sum(1-c[r] for r in range(64)) for c in single_profiles])
    dual_total = np.mean([sum(1-c[r] for r in range(64)) for c in dual_profiles])
    print(f"\n  Mean total carry=0: single={single_total:.2f}, dual={dual_total:.2f}")

    # Gap carry=0 conditioned on carry[63]=0
    if dual_c0 >= 5:
        print(f"\n  Gap carry=0 | carry[63]=0:")
        dual_c0_gap = []
        for c in dual_profiles:
            if c[63] == 0:
                dual_c0_gap.append(sum(1-c[r] for r in gap))
        print(f"    Dual: mean gap c0 = {np.mean(dual_c0_gap):.2f} (N={len(dual_c0_gap)})")

    if single_c0 >= 5:
        single_c0_gap = []
        for c in single_profiles:
            if c[63] == 0:
                single_c0_gap.append(sum(1-c[r] for r in gap))
        print(f"    Single: mean gap c0 = {np.mean(single_c0_gap):.2f} (N={len(single_c0_gap)})")

    return single_c0, dual_c0, improved_rounds


# ============================================================
# M3: Full web comparison — carry heatmap
# ============================================================

def experiment_M3(best_k, N=1500, seed=102):
    """
    Full carry heatmap: single vs dual vs triple parameter.
    """
    print("\n" + "=" * 70)
    print(f"M3: Carry Heatmap — 1 vs 2 vs 3 parameters")
    print(f"N={N}, seed={seed}")
    print("=" * 70)

    rng = np.random.RandomState(seed)

    # Try different k values: best_k, and also k that enters gap directly
    # From M1 analysis: W[2]→W[17], W[3]→W[18], W[9]→W[16,25], W[11]→W[18,27]
    k_candidates = [best_k]
    # Add k values that directly affect gap rounds via schedule
    for k in [2, 3, 7, 9, 11]:
        if k != best_k and k not in k_candidates:
            k_candidates.append(k)
    k_candidates = k_candidates[:4]

    for k_test in k_candidates:
        # Collect carry profiles with HC on W[0]+W[k_test]
        profiles = []
        n_c0 = 0

        t0 = time.time()
        for i in range(N):
            W0 = int(rng.randint(0, 1 << 32))
            Wk = int(rng.randint(0, 1 << 32))
            W16 = [0] * 16
            W16[0] = W0
            W16[k_test] = Wk
            _, raws, _ = sha256_full_trace(W16)
            best_raw = raws[63]

            for step in range(150):
                which = 0 if step % 2 == 0 else k_test
                b = int(rng.randint(0, 32))
                W16_try = list(W16)
                W16_try[which] ^= (1 << b)
                _, raws_try, _ = sha256_full_trace(W16_try)
                if raws_try[63] < best_raw:
                    W16 = W16_try
                    best_raw = raws_try[63]

            c, _, _ = sha256_full_trace(W16)
            profiles.append(c)
            if c[63] == 0:
                n_c0 += 1

        # Carry heatmap
        p_c0 = np.zeros(64)
        for cp in profiles:
            for r in range(64):
                p_c0[r] += (1 - cp[r])
        p_c0 /= N

        gap_mean = np.mean(p_c0[18:47])
        gap_max = np.max(p_c0[18:47])

        # Visual map
        bar = ""
        for r in range(64):
            p = p_c0[r]
            if p > 0.10: bar += '█'
            elif p > 0.05: bar += '▓'
            elif p > 0.02: bar += '▒'
            elif p > 0.005: bar += '░'
            else: bar += '·'

        print(f"\n  W[0]+W[{k_test:2d}]: c63=0: {n_c0}/{N} ({n_c0/N*100:.1f}%), gap_mean={gap_mean:.4f}, gap_max={gap_max:.4f}")
        print(f"  Map: {bar}")

    # Reference: single W[0]
    profiles_ref = []
    n_c0_ref = 0
    for i in range(N):
        W0 = int(rng.randint(0, 1 << 32))
        W16 = [W0] + [0] * 15
        _, raws, _ = sha256_full_trace(W16)
        best_raw = raws[63]
        for step in range(150):
            b = int(rng.randint(0, 32))
            W16_try = list(W16)
            W16_try[0] ^= (1 << b)
            _, raws_try, _ = sha256_full_trace(W16_try)
            if raws_try[63] < best_raw:
                W16 = W16_try
                best_raw = raws_try[63]
        c, _, _ = sha256_full_trace(W16)
        profiles_ref.append(c)
        if c[63] == 0: n_c0_ref += 1

    p_c0_ref = np.zeros(64)
    for cp in profiles_ref:
        for r in range(64): p_c0_ref[r] += (1 - cp[r])
    p_c0_ref /= N

    bar_ref = ""
    for r in range(64):
        p = p_c0_ref[r]
        if p > 0.10: bar_ref += '█'
        elif p > 0.05: bar_ref += '▓'
        elif p > 0.02: bar_ref += '▒'
        elif p > 0.005: bar_ref += '░'
        else: bar_ref += '·'

    gap_ref = np.mean(p_c0_ref[18:47])
    print(f"\n  W[0] only: c63=0: {n_c0_ref}/{N} ({n_c0_ref/N*100:.1f}%), gap_mean={gap_ref:.4f}")
    print(f"  Map: {bar_ref}")
    print(f"\n  Legend: █>10% ▓>5% ▒>2% ░>0.5% ·<0.5%")
    print(f"  Rounds: 0       8      16      24      32      40      48      56   63")

    return k_candidates


# ============================================================
# MAIN
# ============================================================

if __name__ == '__main__':
    print("Multi-Parameter Web — Stage 1.7")
    print(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print()

    N1, N2, N3 = 200, 1500, 1000
    if len(sys.argv) > 1 and sys.argv[1] == '--fast':
        N1, N2, N3 = 100, 800, 500

    t_start = time.time()

    sensitivity, best_k = experiment_M1(N=N1, seed=100)
    single_c0, dual_c0, improved = experiment_M2(best_k, N=N2, seed=101)
    k_candidates = experiment_M3(best_k, N=N3, seed=102)

    total = time.time() - t_start
    print(f"\n{'='*70}")
    print(f"TOTAL TIME: {total:.1f}s")
    print(f"{'='*70}")

    print(f"""
{'='*70}
SYNTHESIS: Multi-Parameter Web (Stage 1.7)
{'='*70}

QUESTION: Can a second parameter W[{best_k}] bridge the gap (r=18-46)?

Single HC (W[0] only): carry[63]=0 rate = {single_c0}/{N2}
Dual HC (W[0]+W[{best_k}]): carry[63]=0 rate = {dual_c0}/{N2}
Speedup: {dual_c0/max(single_c0,1):.1f}×
Gap rounds improved: {improved}

If dual fills the gap → two-island structure bridged
If not → gap is STRUCTURAL, not parameter-limited

DEEPER QUESTION: Is the gap between islands a fundamental property
of SHA-256 schedule, or an artifact of having too few parameters?
With all 16 words active → gap should close (all schedule words ≠ 0).
But then we lose HC efficiency (512 bits to optimize, not 32).
""")
