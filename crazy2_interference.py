#!/usr/bin/env python3
"""
CRAZY-2: Differential Path Interference (Quantum-Inspired)
===========================================================
Test whether two differential paths can "interfere" — their combined
effect on later SHA-256 rounds partially cancels due to near-linearity.

If differential propagation were perfectly linear:
    De(δ1⊕δ2) = De(δ1) ⊕ De(δ2)  =>  residue R_NL = 0

If perfectly random: HW(R_NL) ≈ 16

We search for structured (δ1, δ2) pairs where the nonlinear residue
stays small even after many rounds.
"""

import random
import time
from collections import defaultdict

# ── SHA-256 constants ──
K = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
]

MASK = 0xFFFFFFFF

def rotr(x, n):
    return ((x >> n) | (x << (32 - n))) & MASK

def shr(x, n):
    return x >> n

def ch(x, y, z):
    return (x & y) ^ (~x & z) & MASK

def maj(x, y, z):
    return (x & y) ^ (x & z) ^ (y & z)

def sigma0(x):
    return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)

def sigma1(x):
    return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)

def hw(x):
    """Hamming weight of 32-bit integer."""
    return bin(x & MASK).count('1')

def sha256_round(state, w, ki):
    """One SHA-256 compression round. state = [a,b,c,d,e,f,g,h]."""
    a, b, c, d, e, f, g, h = state
    t1 = (h + sigma1(e) + ch(e, f, g) + ki + w) & MASK
    t2 = (sigma0(a) + maj(a, b, c)) & MASK
    new_state = [
        (t1 + t2) & MASK,  # a
        a,                   # b
        b,                   # c
        c,                   # d
        (d + t1) & MASK,    # e
        e,                   # f
        f,                   # g
        g,                   # h
    ]
    return new_state

def run_rounds(state, W, num_rounds):
    """Run num_rounds of SHA-256 compression, return list of states per round."""
    s = list(state)
    states = [tuple(s)]
    for r in range(num_rounds):
        s = sha256_round(s, W[r], K[r])
        states.append(tuple(s))
    return states

def state_diff(s1, s2):
    """XOR difference between two 8-word states."""
    return tuple((a ^ b) & MASK for a, b in zip(s1, s2))

def state_hw(diff):
    """Total Hamming weight of state difference."""
    return sum(hw(d) for d in diff)

def rand_state():
    return [random.getrandbits(32) for _ in range(8)]

def rand_w(num_rounds):
    return [random.getrandbits(32) for _ in range(num_rounds)]

# ══════════════════════════════════════════════════════════════════════
# EXPERIMENT
# ══════════════════════════════════════════════════════════════════════

NUM_ROUNDS = 20
N_RANDOM = 5000
N_STRUCTURED = 5000

print("=" * 72)
print("CRAZY-2: Differential Path Interference (Quantum-Inspired)")
print("=" * 72)
print(f"Rounds: {NUM_ROUNDS}, Random trials: {N_RANDOM}, Structured trials: {N_STRUCTURED}")
print()

t0 = time.time()

# ── 1. RANDOM δ PAIRS ──
print("─" * 72)
print("Phase 1: Random δ1, δ2 pairs")
print("─" * 72)

random_residue_hw = defaultdict(list)  # round -> list of HW(R_NL)
best_random = []  # (round, hw, δ1, δ2)

for trial in range(N_RANDOM):
    state = rand_state()
    W = rand_w(NUM_ROUNDS)

    delta1 = random.getrandbits(32)
    delta2 = random.getrandbits(32)
    if delta1 == 0 or delta2 == 0 or delta1 == delta2:
        continue

    delta_combined = delta1 ^ delta2

    # W perturbed at position 0 only
    W_a = list(W); W_a[0] ^= delta1
    W_b = list(W); W_b[0] ^= delta2
    W_c = list(W); W_c[0] ^= delta_combined

    states_orig = run_rounds(state, W, NUM_ROUNDS)
    states_a = run_rounds(state, W_a, NUM_ROUNDS)
    states_b = run_rounds(state, W_b, NUM_ROUNDS)
    states_c = run_rounds(state, W_c, NUM_ROUNDS)

    for r in range(1, NUM_ROUNDS + 1):
        de_a = state_diff(states_orig[r], states_a[r])
        de_b = state_diff(states_orig[r], states_b[r])
        de_c = state_diff(states_orig[r], states_c[r])

        # Nonlinear residue: De_combined ⊕ De_A ⊕ De_B
        residue = tuple((a ^ b ^ c) & MASK for a, b, c in zip(de_a, de_b, de_c))
        h = state_hw(residue)
        random_residue_hw[r].append(h)

        if r >= 10 and h < 8:
            best_random.append((r, h, delta1, delta2, trial))

print(f"  Completed {N_RANDOM} random trials in {time.time()-t0:.1f}s")
print()

# Print mean HW per round
print("  Round | Mean HW(R_NL) | Min HW | % with HW<8 | % with HW<4")
print("  " + "-" * 65)
for r in range(1, NUM_ROUNDS + 1):
    vals = random_residue_hw[r]
    if not vals:
        continue
    mean_hw = sum(vals) / len(vals)
    min_hw = min(vals)
    pct_lt8 = 100.0 * sum(1 for v in vals if v < 8) / len(vals)
    pct_lt4 = 100.0 * sum(1 for v in vals if v < 4) / len(vals)
    marker = " <--" if mean_hw < 14 else ""
    print(f"  {r:5d} | {mean_hw:13.2f} | {min_hw:6d} | {pct_lt8:10.4f}% | {pct_lt4:10.4f}%{marker}")

print()

# ── 2. STRUCTURED δ PAIRS ──
print("─" * 72)
print("Phase 2: Structured δ1, δ2 pairs")
print("─" * 72)

SHA_ROTATIONS = [2, 6, 11, 13, 22, 25]

structured_configs = []

# Config A: δ2 = ROTR(δ1, k) for SHA-256 rotation constants
for k in SHA_ROTATIONS:
    structured_configs.append(("ROTR", k))

# Config B: Adjacent bit flips
for bit in range(31):
    structured_configs.append(("ADJ_BITS", bit))

# Config C: MSB + LSB
structured_configs.append(("MSB_LSB", None))

# Config D: Single-bit δ1 with ROTR δ2
for bit in [0, 1, 15, 16, 30, 31]:
    for k in SHA_ROTATIONS:
        structured_configs.append(("BIT_ROTR", (bit, k)))

structured_residue_hw = defaultdict(lambda: defaultdict(list))  # config_name -> round -> HW list
best_structured = []

t1 = time.time()

trials_per_config = max(1, N_STRUCTURED // len(structured_configs))

for cfg_idx, (cfg_type, cfg_param) in enumerate(structured_configs):
    for trial in range(trials_per_config):
        state = rand_state()
        W = rand_w(NUM_ROUNDS)

        if cfg_type == "ROTR":
            # δ1 random, δ2 = ROTR(δ1, k)
            delta1 = random.getrandbits(32)
            delta2 = rotr(delta1, cfg_param)
            label = f"ROTR-{cfg_param}"
        elif cfg_type == "ADJ_BITS":
            # Adjacent single-bit flips
            bit = cfg_param
            delta1 = 1 << bit
            delta2 = 1 << (bit + 1)
            label = f"ADJ-{bit}"
        elif cfg_type == "MSB_LSB":
            delta1 = 0x80000000
            delta2 = 0x00000001
            label = "MSB+LSB"
        elif cfg_type == "BIT_ROTR":
            bit, k = cfg_param
            delta1 = 1 << bit
            delta2 = rotr(delta1, k)
            label = f"BIT{bit}_ROTR{k}"
        else:
            continue

        if delta1 == 0 or delta2 == 0 or delta1 == delta2:
            continue

        delta_combined = delta1 ^ delta2

        W_a = list(W); W_a[0] ^= delta1
        W_b = list(W); W_b[0] ^= delta2
        W_c = list(W); W_c[0] ^= delta_combined

        states_orig = run_rounds(state, W, NUM_ROUNDS)
        states_a = run_rounds(state, W_a, NUM_ROUNDS)
        states_b = run_rounds(state, W_b, NUM_ROUNDS)
        states_c = run_rounds(state, W_c, NUM_ROUNDS)

        for r in range(1, NUM_ROUNDS + 1):
            de_a = state_diff(states_orig[r], states_a[r])
            de_b = state_diff(states_orig[r], states_b[r])
            de_c = state_diff(states_orig[r], states_c[r])

            residue = tuple((a ^ b ^ c) & MASK for a, b, c in zip(de_a, de_b, de_c))
            h = state_hw(residue)
            structured_residue_hw[label][r].append(h)

            if r >= 10 and h < 8:
                best_structured.append((r, h, delta1, delta2, label))

print(f"  Completed {len(structured_configs)} configs x {trials_per_config} trials in {time.time()-t1:.1f}s")
print()

# Summarize structured results by config type at key rounds
print("  Structured δ residue summary (mean HW at selected rounds):")
print(f"  {'Config':<20s} | R=1  | R=3  | R=5  | R=8  | R=10 | R=13 | R=15 | R=17 | R=20")
print("  " + "-" * 100)

# Group by config category
config_categories = defaultdict(lambda: defaultdict(list))
for label in sorted(structured_residue_hw.keys()):
    # Determine category
    if label.startswith("ROTR"):
        cat = label
    elif label.startswith("ADJ"):
        cat = "ADJ-bits (avg)"
    elif label.startswith("MSB"):
        cat = label
    elif label.startswith("BIT"):
        cat = label.split("_")[0] + " (avg)"
    else:
        cat = label

    for r in range(1, NUM_ROUNDS + 1):
        config_categories[cat][r].extend(structured_residue_hw[label][r])

key_rounds = [1, 3, 5, 8, 10, 13, 15, 17, 20]
for cat in sorted(config_categories.keys()):
    vals = []
    for r in key_rounds:
        data = config_categories[cat][r]
        if data:
            vals.append(f"{sum(data)/len(data):5.1f}")
        else:
            vals.append("  N/A")
    print(f"  {cat:<20s} | {'| '.join(vals)}")

# Add random baseline
vals = []
for r in key_rounds:
    data = random_residue_hw[r]
    if data:
        vals.append(f"{sum(data)/len(data):5.1f}")
    else:
        vals.append("  N/A")
print(f"  {'RANDOM (baseline)':<20s} | {'| '.join(vals)}")

print()

# ── 3. BOOMERANG VARIANT ──
print("─" * 72)
print("Phase 3: Boomerang-style interference")
print("─" * 72)

N_BOOM = 3000
BOOM_SPLIT = 10  # split rounds
boomerang_results = defaultdict(list)  # (fwd_rounds, bwd_rounds) -> HW list

t2 = time.time()

for trial in range(N_BOOM):
    state = rand_state()
    W = rand_w(NUM_ROUNDS)

    delta1 = random.getrandbits(32)
    delta2 = random.getrandbits(32)
    if delta1 == 0 or delta2 == 0:
        continue

    # Forward path: R rounds with original and δ1-perturbed W
    W_fwd = list(W); W_fwd[0] ^= delta1

    for split in [5, 8, 10, 12, 15]:
        if split > NUM_ROUNDS:
            continue

        states_orig = run_rounds(state, W, split)
        states_fwd = run_rounds(state, W_fwd, split)

        # At the split point, inject δ2 into the state difference
        state_at_split_orig = list(states_orig[split])
        state_at_split_fwd = list(states_fwd[split])

        # Apply δ2 to 'e' register of both
        state_at_split_orig_mod = list(state_at_split_orig)
        state_at_split_fwd_mod = list(state_at_split_fwd)
        state_at_split_orig_mod[4] ^= delta2
        state_at_split_fwd_mod[4] ^= delta2

        # Continue forward from split to end
        remaining = NUM_ROUNDS - split
        if remaining <= 0:
            continue

        W_rest = W[split:]

        s1 = run_rounds(state_at_split_orig, W_rest, remaining)
        s2 = run_rounds(state_at_split_fwd, W_rest, remaining)
        s3 = run_rounds(state_at_split_orig_mod, W_rest, remaining)
        s4 = run_rounds(state_at_split_fwd_mod, W_rest, remaining)

        # Boomerang check: diff(s1_end, s2_end) ⊕ diff(s3_end, s4_end) = ?
        d12 = state_diff(s1[remaining], s2[remaining])
        d34 = state_diff(s3[remaining], s4[remaining])
        boom_residue = tuple((a ^ b) & MASK for a, b in zip(d12, d34))
        h = state_hw(boom_residue)
        boomerang_results[split].append(h)

print(f"  Completed {N_BOOM} boomerang trials in {time.time()-t2:.1f}s")
print()
print("  Split Round | Mean HW(boom_residue) | Min HW | % HW<8 | % HW<16")
print("  " + "-" * 65)
for split in sorted(boomerang_results.keys()):
    vals = boomerang_results[split]
    mean_h = sum(vals) / len(vals)
    min_h = min(vals)
    pct8 = 100.0 * sum(1 for v in vals if v < 8) / len(vals)
    pct16 = 100.0 * sum(1 for v in vals if v < 16) / len(vals)
    print(f"  {split:11d} | {mean_h:21.2f} | {min_h:6d} | {pct8:5.2f}% | {pct16:5.2f}%")

print()

# ── 4. DEEP SEARCH for low-residue pairs ──
print("─" * 72)
print("Phase 4: Deep search for low-residue pairs at high rounds")
print("─" * 72)

t3 = time.time()

# Try single-bit δ pairs exhaustively
N_DEEP = 2000
deep_best = []  # (round, hw, δ1, δ2, type)

# Single-bit pairs
for b1 in range(32):
    for b2 in range(b1 + 1, 32):
        delta1 = 1 << b1
        delta2 = 1 << b2
        delta_combined = delta1 ^ delta2

        for trial in range(30):
            state = rand_state()
            W = rand_w(NUM_ROUNDS)

            W_a = list(W); W_a[0] ^= delta1
            W_b = list(W); W_b[0] ^= delta2
            W_c = list(W); W_c[0] ^= delta_combined

            states_orig = run_rounds(state, W, NUM_ROUNDS)
            states_a = run_rounds(state, W_a, NUM_ROUNDS)
            states_b = run_rounds(state, W_b, NUM_ROUNDS)
            states_c = run_rounds(state, W_c, NUM_ROUNDS)

            for r in [10, 13, 15, 17, 20]:
                de_a = state_diff(states_orig[r], states_a[r])
                de_b = state_diff(states_orig[r], states_b[r])
                de_c = state_diff(states_orig[r], states_c[r])
                residue = tuple((a ^ b ^ c) & MASK for a, b, c in zip(de_a, de_b, de_c))
                h = state_hw(residue)
                if h < 8:
                    deep_best.append((r, h, delta1, delta2, f"bits({b1},{b2})"))

print(f"  Single-bit pair search: {32*31//2} pairs x 30 trials in {time.time()-t3:.1f}s")
print(f"  Found {len(deep_best)} instances with HW < 8 at round 10+")

# Show distribution
if deep_best:
    from collections import Counter
    round_counts = Counter()
    for r, h, d1, d2, tp in deep_best:
        round_counts[r] += 1
    print("  Breakdown by round:")
    for r in sorted(round_counts.keys()):
        total_trials = 32 * 31 // 2 * 30
        pct = 100.0 * round_counts[r] / total_trials
        print(f"    Round {r}: {round_counts[r]} instances ({pct:.4f}%)")

    # Best examples
    deep_best.sort(key=lambda x: (x[0], x[1]))  # sort by round desc, then hw
    print("  Best examples (highest round, lowest HW):")
    seen = set()
    count = 0
    for r, h, d1, d2, tp in sorted(deep_best, key=lambda x: (-x[0], x[1])):
        key = (r, tp)
        if key not in seen:
            seen.add(key)
            print(f"    Round {r}: HW={h}, δ1=0x{d1:08x}, δ2=0x{d2:08x} ({tp})")
            count += 1
            if count >= 15:
                break

print()

# ── 5. LINEARITY DECAY ANALYSIS ──
print("─" * 72)
print("Phase 5: Linearity Decay Analysis")
print("─" * 72)

print()
print("  Linearity decay curve (mean HW of nonlinear residue per round):")
print()
print("  Round | Random δ | Expected(random)")
print("  " + "-" * 45)
for r in range(1, NUM_ROUNDS + 1):
    rv = random_residue_hw[r]
    mean_r = sum(rv) / len(rv) if rv else 0
    bar_len = int(mean_r * 2)
    bar = "█" * bar_len
    print(f"  {r:5d} | {mean_r:8.2f} | {16:15d}  {bar}")

print()

# ── 6. STRUCTURED vs RANDOM comparison ──
print("─" * 72)
print("Phase 6: Structured vs Random — statistical comparison")
print("─" * 72)

print()
comparison_rounds = [5, 10, 15, 20]
print(f"  {'Category':<25s}", end="")
for r in comparison_rounds:
    print(f" | R={r:<3d} mean(HW)", end="")
print()
print("  " + "-" * 90)

# Random baseline
print(f"  {'Random δ':<25s}", end="")
for r in comparison_rounds:
    rv = random_residue_hw[r]
    print(f" | {sum(rv)/len(rv):14.2f}", end="")
print()

# Each structured category
for cat in sorted(config_categories.keys()):
    print(f"  {cat:<25s}", end="")
    for r in comparison_rounds:
        data = config_categories[cat][r]
        if data:
            print(f" | {sum(data)/len(data):14.2f}", end="")
        else:
            print(f" | {'N/A':>14s}", end="")
    print()

print()

# ── 7. Check for HW<4 at round 17+ ──
print("─" * 72)
print("Phase 7: Search for HW(R_NL) < 4 at round 17+")
print("─" * 72)

found_17plus = []
for r, h, d1, d2, info in best_random:
    if r >= 17 and h < 4:
        found_17plus.append(("random", r, h, d1, d2, info))

for r, h, d1, d2, label in best_structured:
    if r >= 17 and h < 4:
        found_17plus.append(("structured", r, h, d1, d2, label))

for r, h, d1, d2, tp in deep_best:
    if r >= 17 and h < 4:
        found_17plus.append(("deep", r, h, d1, d2, tp))

print(f"  Total instances with HW < 4 at round 17+: {len(found_17plus)}")
if found_17plus:
    found_17plus.sort(key=lambda x: (x[1], -x[2]))
    for src, r, h, d1, d2, info in found_17plus[:20]:
        print(f"    [{src}] Round {r}: HW={h}, δ1=0x{d1:08x}, δ2=0x{d2:08x} ({info})")
else:
    print("  None found.")

print()

# ── FINAL VERDICT ──
print("=" * 72)
print("VERDICT")
print("=" * 72)

# Check if structured is meaningfully better than random
structured_better = False
random_r15 = random_residue_hw[15]
mean_random_r15 = sum(random_r15) / len(random_r15) if random_r15 else 16

best_struct_mean_r15 = 999
best_struct_cat = ""
for cat, rdata in config_categories.items():
    if rdata[15]:
        m = sum(rdata[15]) / len(rdata[15])
        if m < best_struct_mean_r15:
            best_struct_mean_r15 = m
            best_struct_cat = cat

improvement = mean_random_r15 - best_struct_mean_r15
pct_improvement = 100.0 * improvement / mean_random_r15 if mean_random_r15 > 0 else 0

print(f"  Random baseline at R=15: mean HW = {mean_random_r15:.2f}")
print(f"  Best structured at R=15: mean HW = {best_struct_mean_r15:.2f} ({best_struct_cat})")
print(f"  Improvement: {improvement:.2f} bits ({pct_improvement:.1f}%)")
print()

# Check round 1 linearity
r1_mean = sum(random_residue_hw[1]) / len(random_residue_hw[1])
print(f"  Round 1 residue (measures immediate nonlinearity): mean HW = {r1_mean:.2f}")
print(f"  (0 = perfectly linear, 16 = fully random)")
print()

# Determine alive/dead
if pct_improvement > 10 or len(found_17plus) > 5:
    status = "ALIVE"
    reason = "Structured pairs show measurably lower nonlinear residue"
elif pct_improvement > 3:
    status = "WEAK SIGNAL"
    reason = "Small but detectable difference between structured and random"
elif r1_mean < 8:
    status = "INTERESTING (early rounds)"
    reason = f"Strong linearity at round 1 (HW={r1_mean:.1f}) but decays quickly"
else:
    status = "DEAD"
    reason = "No meaningful difference; nonlinear residue ≈ random at all rounds"

print(f"  Status: {status}")
print(f"  Reason: {reason}")
print()

# Round 1 analysis - the first round IS partially linear
print("  KEY INSIGHT: Round 1 linearity analysis")
print(f"    Mean HW at R=1: {r1_mean:.2f} (if linear: 0, if random: 16)")
r2_mean = sum(random_residue_hw[2]) / len(random_residue_hw[2])
r3_mean = sum(random_residue_hw[3]) / len(random_residue_hw[3])
print(f"    Mean HW at R=2: {r2_mean:.2f}")
print(f"    Mean HW at R=3: {r3_mean:.2f}")
print(f"    Linearity breaks down rapidly (as expected for SHA-256)")

elapsed = time.time() - t0
print()
print(f"  Total runtime: {elapsed:.1f}s")
print("=" * 72)
