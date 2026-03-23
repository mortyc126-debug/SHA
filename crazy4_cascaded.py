#!/usr/bin/env python3
"""
CRAZY-4: Cascaded Navigation
=============================
Instead of optimizing |De_mod| at round 17 only, try CASCADED optimization
across rounds 17->18->19->20.

Previous findings:
  - C9: ~82.6 neutral bits exist for Wang chain
  - C9xC13: Greedy neutral bit selection reduces |De_mod[17]| from 2^30 to 2^25
  - BUT: savings don't persist to round 18

Key idea: Phase-by-phase cascaded optimization:
  Phase 1: Greedy-select neutral bits to minimize |De_mod[17]|
  Phase 2: From Phase 1 solution, try ADDITIONAL flips to minimize |De_mod[18]|
           without worsening |De_mod[17]| by more than 2x
  Phase 3: Same for round 19 (allow 2x degradation of rounds 17-18)
  Phase 4: Same for round 20

Also: Joint GA optimization minimizing sum of log2(|De_mod[t]|) for t=17..20.

Compare: sequential cascade vs joint optimization vs single-round greedy.
"""

import random
import time
import math
import numpy as np
from collections import defaultdict

MASK = 0xFFFFFFFF
MOD = 2**32

# SHA-256 round constants
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
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
]

# SHA-256 initial hash
H0 = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
]

# ── SHA-256 primitives ──────────────────────────────────────────────

def rotr(x, n):
    return ((x >> n) | (x << (32 - n))) & MASK

def Sig0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sig1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def sig0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sig1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)
def Ch(e, f, g):  return ((e & f) ^ ((~e) & g)) & MASK
def Maj(a, b, c): return ((a & b) ^ (a & c) ^ (b & c)) & MASK

def add32(*args):
    s = 0
    for a in args:
        s = (s + a) & MASK
    return s

def sha256_round(state, W_t, K_t):
    a, b, c, d, e, f, g, h = state
    T1 = add32(h, Sig1(e), Ch(e, f, g), K_t, W_t)
    T2 = add32(Sig0(a), Maj(a, b, c))
    new_a = add32(T1, T2)
    new_e = add32(d, T1)
    return [new_a, a, b, c, new_e, e, f, g]

def expand_message(W16):
    W = list(W16)
    for i in range(16, 64):
        W.append(add32(sig1(W[i-2]), W[i-7], sig0(W[i-15]), W[i-16]))
    return W

def run_rounds(init_state, W, num_rounds):
    states = [list(init_state)]
    for t in range(num_rounds):
        s = sha256_round(states[-1], W[t], K[t])
        states.append(s)
    return states

def mod_diff_signed(a, b):
    d = (a - b) & MASK
    if d >= (1 << 31):
        return d - MOD
    return d

def abs_mod_diff(a, b):
    return abs(mod_diff_signed(a, b))

# ── Wang chain ──────────────────────────────────────────────────────

def compute_wang_chain(base_msg, init_state, dw0=0x80000000):
    W = list(base_msg)
    W_prime = list(W)
    W_prime[0] = W[0] ^ dw0

    state = list(init_state)
    state_prime = list(init_state)

    states = [list(state)]
    states_prime = [list(state_prime)]

    state = sha256_round(state, W[0], K[0])
    state_prime = sha256_round(state_prime, W_prime[0], K[0])
    states.append(list(state))
    states_prime.append(list(state_prime))

    for t in range(1, 16):
        a, b, c, d, e, f, g, h = state
        a2, b2, c2, d2, e2, f2, g2, h2 = state_prime

        T1_partial = add32(h, Sig1(e), Ch(e, f, g), K[t])
        T1_partial_prime = add32(h2, Sig1(e2), Ch(e2, f2, g2), K[t])

        target_e = add32(d, T1_partial, W[t])
        W_prime_t = (target_e - d2 - T1_partial_prime) & MASK
        W_prime[t] = W_prime_t

        state = sha256_round(state, W[t], K[t])
        state_prime = sha256_round(state_prime, W_prime_t, K[t])
        states.append(list(state))
        states_prime.append(list(state_prime))

    return W_prime, states, states_prime

def verify_wang_chain(states, states_prime, check_rounds=16):
    for t in range(2, min(check_rounds + 1, len(states))):
        de = (states[t][4] - states_prime[t][4]) & MASK
        if de != 0:
            return False
    return True

# ── Neutral bit detection ───────────────────────────────────────────

def find_neutral_bits(base_msg, init_state):
    W_prime, states, states_prime = compute_wang_chain(base_msg, init_state)
    neutrals = []

    for k in range(1, 16):
        for b in range(32):
            flip = 1 << b
            W_f = list(base_msg)
            W_f[k] ^= flip
            Wp_f = list(W_prime)
            Wp_f[k] ^= flip

            W_exp = expand_message(W_f)
            Wp_exp = expand_message(Wp_f)

            s1 = run_rounds(init_state, W_exp, 16)
            s2 = run_rounds(init_state, Wp_exp, 16)

            ok = True
            for t in range(2, 17):
                de = (s1[t][4] - s2[t][4]) & MASK
                if de != 0:
                    ok = False
                    break

            if ok:
                neutrals.append((k, b))

    return neutrals

# ── Modular differential tracking (extended to round 21) ───────────

def get_de_mod_rounds(base_msg, init_state, W_prime_base, neutrals_to_flip, max_round=21):
    """
    Flip the specified neutral bits in both W and W', expand, run rounds,
    return dict: round_number -> signed De_mod for rounds 17..max_round.
    """
    W = list(base_msg)
    Wp = list(W_prime_base)

    for (k, b) in neutrals_to_flip:
        flip = 1 << b
        W[k] ^= flip
        Wp[k] ^= flip

    W_exp = expand_message(W)
    Wp_exp = expand_message(Wp)

    s1 = run_rounds(init_state, W_exp, max_round)
    s2 = run_rounds(init_state, Wp_exp, max_round)

    de_mod = {}
    for t in range(17, max_round + 1):
        if t < len(s1):
            de_mod[t] = mod_diff_signed(s1[t][4], s2[t][4])
    return de_mod


def safe_log2(val):
    """log2 of absolute value, with floor at 0."""
    av = abs(val)
    if av == 0:
        return 0.0
    return math.log2(av)


# ── Phase-based Cascaded Greedy ─────────────────────────────────────

def cascaded_greedy(base_msg, init_state, W_prime, neutrals, max_steps_per_phase=30):
    """
    Phase 1: minimize |De_mod[17]|
    Phase 2: minimize |De_mod[18]| without worsening |De_mod[17]| by >2x
    Phase 3: minimize |De_mod[19]| without worsening 17,18 by >2x
    Phase 4: minimize |De_mod[20]| without worsening 17,18,19 by >2x

    Returns (final_flips, de_mod_after_each_phase).
    de_mod_after_each_phase: list of dicts {17:val, 18:val, ...} after each phase.
    """
    current_flips = []
    available = list(neutrals)
    phase_results = []

    for phase, target_round in enumerate([17, 18, 19, 20]):
        # Record thresholds: allow 2x degradation of previously optimized rounds
        if phase == 0:
            thresholds = {}  # no constraints in phase 1
        else:
            prev_de = get_de_mod_rounds(base_msg, init_state, W_prime, current_flips)
            thresholds = {}
            for r in range(17, target_round):
                thresholds[r] = abs(prev_de.get(r, 0)) * 2.0 + 1  # +1 to handle zero

        steps = min(len(available), max_steps_per_phase)
        for step in range(steps):
            best_bit = None
            best_target_val = float('inf')

            # Current target value
            cur_de = get_de_mod_rounds(base_msg, init_state, W_prime, current_flips)
            cur_target = abs(cur_de.get(target_round, 0))

            for (k, b) in available:
                test_flips = current_flips + [(k, b)]
                de = get_de_mod_rounds(base_msg, init_state, W_prime, test_flips)

                # Check constraints: don't worsen previous rounds by >2x
                violated = False
                for r, thresh in thresholds.items():
                    if abs(de.get(r, 0)) > thresh:
                        violated = True
                        break
                if violated:
                    continue

                target_val = abs(de.get(target_round, 0))
                if target_val < best_target_val:
                    best_target_val = target_val
                    best_bit = (k, b)

            if best_bit is not None and best_target_val < cur_target:
                current_flips.append(best_bit)
                available.remove(best_bit)
            else:
                break

        # Record results after this phase
        de_after = get_de_mod_rounds(base_msg, init_state, W_prime, current_flips)
        phase_results.append(dict(de_after))

    return current_flips, phase_results


# ── Single-round greedy (optimize r17 only) ─────────────────────────

def single_round_greedy(base_msg, init_state, W_prime, neutrals, max_steps=40):
    """Greedy to minimize |De_mod[17]| only. Returns (flips, de_mod_dict)."""
    current_flips = []
    available = list(neutrals)

    for step in range(min(len(available), max_steps)):
        cur_de = get_de_mod_rounds(base_msg, init_state, W_prime, current_flips)
        cur_17 = abs(cur_de.get(17, 0))
        best_bit = None
        best_val = cur_17

        for (k, b) in available:
            test_flips = current_flips + [(k, b)]
            de = get_de_mod_rounds(base_msg, init_state, W_prime, test_flips)
            val = abs(de.get(17, 0))
            if val < best_val:
                best_val = val
                best_bit = (k, b)

        if best_bit is not None:
            current_flips.append(best_bit)
            available.remove(best_bit)
        else:
            break

    final_de = get_de_mod_rounds(base_msg, init_state, W_prime, current_flips)
    return current_flips, final_de


# ── Joint GA/SA optimization ────────────────────────────────────────

def joint_ga_optimization(base_msg, init_state, W_prime, neutrals,
                          pop_size=40, generations=60, rng=None):
    """
    GA that minimizes sum of log2(|De_mod[t]|) for t=17..20.
    Each individual is a bitmask over neutral bits.
    Returns (best_flips, best_de_mod_dict).
    """
    if rng is None:
        rng = random.Random()

    n = len(neutrals)
    if n == 0:
        return [], {}

    def fitness(mask):
        flips = [neutrals[i] for i in range(n) if mask[i]]
        de = get_de_mod_rounds(base_msg, init_state, W_prime, flips)
        score = 0.0
        for t in range(17, 21):
            score += safe_log2(de.get(t, 0))
        return score

    # Initialize population: mix of empty, random sparse, and single-bit
    population = []
    # Empty
    population.append([0] * n)
    # Single bits (best subset)
    for i in range(min(n, pop_size // 4)):
        ind = [0] * n
        ind[i] = 1
        population.append(ind)
    # Random sparse
    while len(population) < pop_size:
        ind = [0] * n
        density = rng.uniform(0.02, 0.3)
        for i in range(n):
            if rng.random() < density:
                ind[i] = 1
        population.append(ind)

    # Evaluate
    scores = [fitness(ind) for ind in population]

    best_idx = int(np.argmin(scores))
    best_score = scores[best_idx]
    best_ind = list(population[best_idx])

    for gen in range(generations):
        # Tournament selection + crossover + mutation
        new_pop = [list(best_ind)]  # elitism
        new_scores = [best_score]

        while len(new_pop) < pop_size:
            # Tournament
            i1, i2 = rng.sample(range(len(population)), 2)
            p1 = population[i1] if scores[i1] < scores[i2] else population[i2]
            i3, i4 = rng.sample(range(len(population)), 2)
            p2 = population[i3] if scores[i3] < scores[i4] else population[i4]

            # Uniform crossover
            child = [p1[i] if rng.random() < 0.5 else p2[i] for i in range(n)]

            # Mutation
            mut_rate = 1.0 / max(n, 1)
            for i in range(n):
                if rng.random() < mut_rate:
                    child[i] = 1 - child[i]

            sc = fitness(child)
            new_pop.append(child)
            new_scores.append(sc)

        population = new_pop
        scores = new_scores

        cur_best_idx = int(np.argmin(scores))
        if scores[cur_best_idx] < best_score:
            best_score = scores[cur_best_idx]
            best_ind = list(population[cur_best_idx])

    best_flips = [neutrals[i] for i in range(n) if best_ind[i]]
    best_de = get_de_mod_rounds(base_msg, init_state, W_prime, best_flips)
    return best_flips, best_de


# ── Main experiment ─────────────────────────────────────────────────

def run_experiment(num_base_msgs=30, time_limit=480):
    random.seed(0xCAFE)
    rng = random.Random(0xCAFE)
    start_time = time.time()

    print("=" * 76)
    print("CRAZY-4: Cascaded Navigation")
    print("=" * 76)
    print()
    print(f"Parameters: {num_base_msgs} base messages, seed=0xCAFE, time_limit={time_limit}s")
    print()

    # Per-message results
    results = {
        'baseline': [],       # list of dicts {17:abs_de, 18:..., 19:..., 20:...}
        'single_greedy': [],  # single-round greedy (r17 only)
        'cascade_phase1': [],
        'cascade_phase2': [],
        'cascade_phase3': [],
        'cascade_phase4': [],
        'joint_ga': [],
        'neutral_counts': [],
    }

    for msg_idx in range(num_base_msgs):
        elapsed = time.time() - start_time
        if elapsed > time_limit:
            print(f"\n[Time limit reached at message {msg_idx}, stopping]")
            break

        print(f"  Message {msg_idx+1}/{num_base_msgs} (elapsed: {elapsed:.1f}s)")

        base_msg = [rng.getrandbits(32) for _ in range(16)]
        init_state = list(H0)

        # Compute Wang chain
        W_prime, states, states_prime = compute_wang_chain(base_msg, init_state)

        if not verify_wang_chain(states, states_prime, 16):
            print(f"    WARNING: Wang chain failed for message {msg_idx}")
            continue

        # Baseline: no neutral flips
        baseline_de = get_de_mod_rounds(base_msg, init_state, W_prime, [])
        baseline_abs = {t: abs(baseline_de.get(t, 0)) for t in range(17, 21)}
        results['baseline'].append(baseline_abs)

        # Find neutral bits
        neutrals = find_neutral_bits(base_msg, init_state)
        results['neutral_counts'].append(len(neutrals))

        if len(neutrals) == 0:
            print(f"    No neutral bits found, skipping")
            continue

        # Single-round greedy (optimize r17 only)
        sg_flips, sg_de = single_round_greedy(base_msg, init_state, W_prime, neutrals, max_steps=30)
        sg_abs = {t: abs(sg_de.get(t, 0)) for t in range(17, 21)}
        results['single_greedy'].append(sg_abs)

        # Cascaded greedy
        casc_flips, casc_phases = cascaded_greedy(
            base_msg, init_state, W_prime, neutrals, max_steps_per_phase=20
        )
        for pi, phase_key in enumerate(['cascade_phase1', 'cascade_phase2',
                                         'cascade_phase3', 'cascade_phase4']):
            if pi < len(casc_phases):
                pa = {t: abs(casc_phases[pi].get(t, 0)) for t in range(17, 21)}
                results[phase_key].append(pa)

        # Joint GA optimization (smaller for speed)
        elapsed2 = time.time() - start_time
        remaining = time_limit - elapsed2
        # Scale GA effort based on remaining time
        if remaining > 120 and len(neutrals) <= 120:
            ga_flips, ga_de = joint_ga_optimization(
                base_msg, init_state, W_prime, neutrals,
                pop_size=30, generations=40, rng=random.Random(rng.randint(0, 2**32))
            )
            ga_abs = {t: abs(ga_de.get(t, 0)) for t in range(17, 21)}
            results['joint_ga'].append(ga_abs)
        elif remaining > 60:
            ga_flips, ga_de = joint_ga_optimization(
                base_msg, init_state, W_prime, neutrals,
                pop_size=20, generations=20, rng=random.Random(rng.randint(0, 2**32))
            )
            ga_abs = {t: abs(ga_de.get(t, 0)) for t in range(17, 21)}
            results['joint_ga'].append(ga_abs)
        else:
            print(f"    Skipping GA (time constraint)")

    total_time = time.time() - start_time

    # ── Reporting ───────────────────────────────────────────────────
    print()
    print("=" * 76)
    print("RESULTS: CRAZY-4 Cascaded Navigation")
    print("=" * 76)

    def summarize(label, data, rounds=range(17, 21)):
        """Print summary stats for a list of {round: abs_de} dicts."""
        if not data:
            print(f"\n--- {label}: NO DATA ---")
            return {}
        print(f"\n--- {label} (n={len(data)}) ---")
        stats = {}
        for t in rounds:
            vals = [d.get(t, 0) for d in data]
            arr = np.array(vals, dtype=np.float64)
            nonzero = arr[arr > 0]
            if len(nonzero) > 0:
                log2_vals = np.log2(nonzero)
                mean_log2 = np.mean(log2_vals)
                median_log2 = np.median(log2_vals)
            else:
                mean_log2 = 0.0
                median_log2 = 0.0
            stats[t] = {
                'mean': arr.mean(),
                'median': np.median(arr),
                'mean_log2': mean_log2,
                'median_log2': median_log2,
            }
            print(f"  Round {t}: mean |De_mod| = {arr.mean():.3e}  "
                  f"(mean log2 = {mean_log2:.2f}, median log2 = {median_log2:.2f})")
        return stats

    baseline_stats = summarize("Baseline (no optimization)", results['baseline'])
    sg_stats = summarize("Single-round greedy (optimize r17 only)", results['single_greedy'])
    p1_stats = summarize("Cascade Phase 1 (optimize r17)", results['cascade_phase1'])
    p2_stats = summarize("Cascade Phase 2 (optimize r18, protect r17)", results['cascade_phase2'])
    p3_stats = summarize("Cascade Phase 3 (optimize r19, protect r17-18)", results['cascade_phase3'])
    p4_stats = summarize("Cascade Phase 4 (optimize r20, protect r17-19)", results['cascade_phase4'])
    ga_stats = summarize("Joint GA (minimize sum log2 for r17-20)", results['joint_ga'])

    # ── Comparison table ────────────────────────────────────────────
    print()
    print("=" * 76)
    print("COMPARISON TABLE: Mean log2(|De_mod[t]|) by round")
    print("=" * 76)
    print(f"{'Method':<45} {'r17':>8} {'r18':>8} {'r19':>8} {'r20':>8}")
    print("-" * 76)

    all_stats = [
        ("Baseline", baseline_stats),
        ("Single-round greedy (r17)", sg_stats),
        ("Cascade Phase 1 (r17)", p1_stats),
        ("Cascade Phase 2 (r17+r18)", p2_stats),
        ("Cascade Phase 3 (r17+r18+r19)", p3_stats),
        ("Cascade Phase 4 (r17..r20)", p4_stats),
        ("Joint GA (r17..r20)", ga_stats),
    ]

    for label, st in all_stats:
        if st:
            vals = [f"{st[t]['mean_log2']:>8.2f}" if t in st else f"{'N/A':>8}" for t in range(17, 21)]
            print(f"{label:<45} {vals[0]} {vals[1]} {vals[2]} {vals[3]}")
        else:
            print(f"{label:<45} {'N/A':>8} {'N/A':>8} {'N/A':>8} {'N/A':>8}")

    # ── Savings computation ─────────────────────────────────────────
    print()
    print("=" * 76)
    print("BIT SAVINGS vs BASELINE (mean log2 reduction)")
    print("=" * 76)
    print(f"{'Method':<45} {'r17':>8} {'r18':>8} {'r19':>8} {'r20':>8}")
    print("-" * 76)

    if baseline_stats:
        for label, st in all_stats[1:]:  # skip baseline
            if st and baseline_stats:
                savings = []
                for t in range(17, 21):
                    if t in st and t in baseline_stats:
                        s = baseline_stats[t]['mean_log2'] - st[t]['mean_log2']
                        savings.append(f"{s:>8.2f}")
                    else:
                        savings.append(f"{'N/A':>8}")
                print(f"{label:<45} {savings[0]} {savings[1]} {savings[2]} {savings[3]}")

    # ── Persistence analysis ────────────────────────────────────────
    print()
    print("=" * 76)
    print("PERSISTENCE ANALYSIS")
    print("=" * 76)

    # Check: does cascade phase 2 maintain r17 savings while improving r18?
    if p2_stats and baseline_stats and p1_stats:
        r17_savings_p2 = baseline_stats[17]['mean_log2'] - p2_stats[17]['mean_log2']
        r18_savings_p2 = baseline_stats[18]['mean_log2'] - p2_stats[18]['mean_log2']
        r17_savings_p1 = baseline_stats[17]['mean_log2'] - p1_stats[17]['mean_log2']
        print(f"  Phase 1 r17 savings: {r17_savings_p1:.2f} bits")
        print(f"  Phase 2 r17 savings: {r17_savings_p2:.2f} bits (after r18 optimization)")
        print(f"  Phase 2 r18 savings: {r18_savings_p2:.2f} bits")
        print(f"  r17 degradation from Phase 1->2: {r17_savings_p1 - r17_savings_p2:.2f} bits")
        print()

    if sg_stats and baseline_stats:
        sg_r17_savings = baseline_stats[17]['mean_log2'] - sg_stats[17]['mean_log2']
        sg_r18_savings = baseline_stats[18]['mean_log2'] - sg_stats[18]['mean_log2']
        print(f"  Single-round greedy r17 savings: {sg_r17_savings:.2f} bits")
        print(f"  Single-round greedy r18 savings: {sg_r18_savings:.2f} bits (unoptimized)")
        print(f"  -> Does r17 optimization persist to r18? {'YES' if sg_r18_savings > 1 else 'NO'}")
        print()

    # Per-message cascade detail (first 5 messages)
    print()
    print("=" * 76)
    print("PER-MESSAGE DETAIL (first 5 messages)")
    print("=" * 76)

    for i in range(min(5, len(results['baseline']))):
        print(f"\n  Message {i+1}:")
        bl = results['baseline'][i]
        print(f"    Baseline:    r17={safe_log2(bl[17]):6.1f}  r18={safe_log2(bl[18]):6.1f}  "
              f"r19={safe_log2(bl[19]):6.1f}  r20={safe_log2(bl[20]):6.1f}")

        if i < len(results['single_greedy']):
            sg = results['single_greedy'][i]
            print(f"    Greedy(r17): r17={safe_log2(sg[17]):6.1f}  r18={safe_log2(sg[18]):6.1f}  "
                  f"r19={safe_log2(sg[19]):6.1f}  r20={safe_log2(sg[20]):6.1f}")

        for pi, (key, lbl) in enumerate([
            ('cascade_phase1', 'Casc P1'),
            ('cascade_phase2', 'Casc P2'),
            ('cascade_phase3', 'Casc P3'),
            ('cascade_phase4', 'Casc P4'),
        ]):
            if i < len(results[key]):
                cp = results[key][i]
                print(f"    {lbl}:      r17={safe_log2(cp[17]):6.1f}  r18={safe_log2(cp[18]):6.1f}  "
                      f"r19={safe_log2(cp[19]):6.1f}  r20={safe_log2(cp[20]):6.1f}")

        if i < len(results['joint_ga']):
            ga = results['joint_ga'][i]
            print(f"    Joint GA:    r17={safe_log2(ga[17]):6.1f}  r18={safe_log2(ga[18]):6.1f}  "
                  f"r19={safe_log2(ga[19]):6.1f}  r20={safe_log2(ga[20]):6.1f}")

    # ── Neutral bit statistics ──────────────────────────────────────
    if results['neutral_counts']:
        nc = np.array(results['neutral_counts'])
        print(f"\n  Neutral bits per message: mean={nc.mean():.1f}, min={nc.min()}, max={nc.max()}")

    # ── Verdict ─────────────────────────────────────────────────────
    print()
    print("=" * 76)
    print("VERDICT")
    print("=" * 76)

    # Cascaded r18 savings
    cascade_r18_savings = 0.0
    cascade_r17_maintained = False
    if p2_stats and baseline_stats:
        cascade_r18_savings = baseline_stats[18]['mean_log2'] - p2_stats[18]['mean_log2']
        r17_savings_p2 = baseline_stats[17]['mean_log2'] - p2_stats[17]['mean_log2']
        cascade_r17_maintained = r17_savings_p2 >= 2.0  # at least 2 bits at r17

    # Also check joint GA
    ga_r18_savings = 0.0
    if ga_stats and baseline_stats:
        ga_r18_savings = baseline_stats[18]['mean_log2'] - ga_stats[18]['mean_log2']

    best_r18_savings = max(cascade_r18_savings, ga_r18_savings)

    print()
    print(f"  Cascade r18 savings vs baseline: {cascade_r18_savings:.2f} bits")
    print(f"  Joint GA r18 savings vs baseline: {ga_r18_savings:.2f} bits")
    print(f"  Best r18 savings: {best_r18_savings:.2f} bits")
    print(f"  Cascade r17 maintained (>=2 bits): {cascade_r17_maintained}")
    print()

    if best_r18_savings >= 3.0 and cascade_r17_maintained:
        verdict = "ALIVE"
        explanation = ("Cascaded optimization gives >=3 bit savings at round 18 "
                      "while maintaining round 17 savings.")
    elif best_r18_savings >= 1.0:
        verdict = "ANOMALY"
        explanation = ("Some savings at r18 ({:.2f} bits) but less than 3 bits. "
                      "Partial persistence of optimization.".format(best_r18_savings))
    else:
        verdict = "DEAD"
        explanation = ("No persistent savings: r18 returns to near-random "
                      "regardless of r17 optimization.")

    print(f"  >>> VERDICT: {verdict}")
    print(f"  >>> {explanation}")
    print()
    print(f"  Total runtime: {total_time:.1f}s")
    print()


if __name__ == "__main__":
    run_experiment(num_base_msgs=30, time_limit=480)
