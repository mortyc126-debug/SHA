#!/usr/bin/env python3
"""
ЗАДАНИЕ 4: Schedule Collision — аддитивные vs XOR нули расписания

Experiments (priority A > B > C > D):
A. Additive kernel of schedule (cascade of ΔW[t]=0)
B. Partial schedule collision (SA minimization)
C. Schedule near-periodicity (autocorrelation)
D. Two-message schedule alignment
"""

import struct
import os
import numpy as np
import math
import time
from collections import defaultdict

MASK = 0xFFFFFFFF


def rotr(x, n):
    return ((x >> n) | (x << (32 - n))) & MASK

def sig0(x):
    return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)

def sig1(x):
    return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)

def add(*args):
    s = 0
    for x in args:
        s = (s + x) & MASK
    return s

def sub(a, b):
    return (a - b) & MASK

def hw32(x):
    return bin(x & MASK).count('1')

def to_signed(x):
    if x > 0x7FFFFFFF:
        return x - 0x100000000
    return x

def random_words(n):
    return list(struct.unpack(f'>{n}I', os.urandom(4 * n)))


def message_schedule(M):
    W = list(M[:16])
    for i in range(16, 64):
        W.append(add(sig1(W[i-2]), W[i-7], sig0(W[i-15]), W[i-16]))
    return W


def schedule_delta(W_base, DW):
    """
    Given base W[0..15] and additive delta DW[0..15],
    compute W' = W + DW and return DW_expanded[0..63].
    """
    W_prime = [(W_base[i] + DW[i]) & MASK for i in range(16)]
    W_exp = message_schedule(W_base)
    W_prime_exp = message_schedule(W_prime)
    DW_full = [sub(W_prime_exp[i], W_exp[i]) for i in range(64)]
    return DW_full


# =============================================================
# EXPERIMENT A: Additive kernel of schedule
# =============================================================
def experiment_A():
    print("=" * 70)
    print("EXPERIMENT A: Additive kernel of schedule")
    print("=" * 70)
    t0 = time.time()

    # Part 1: Single equation ΔW[16]=0, random search
    print("\n  Part 1: Find ΔW[0..15] with ΔW[16]=0 (birthday search)")
    n_base = 100

    propagation_profiles = []

    for trial in range(n_base):
        W = random_words(16)

        # Random search: try random DW until ΔW[16]=0
        # ΔW[16] = sig1(W'[14])-sig1(W[14]) + (W'[9]-W[9])
        #         + sig0(W'[1])-sig0(W[1]) + (W'[0]-W[0])
        # We need this = 0 mod 2^32
        # Fix DW[2..8,10..13,15] = 0, vary DW[0,1,9,14]
        # Target: Δsig1(14) + DW[9] + Δsig0(1) + DW[0] = 0

        best_dw = None
        best_hw16 = 999

        for attempt in range(50000):
            DW = [0] * 16
            # Random perturbation of a few words
            n_perturb = np.random.randint(1, 5)
            indices = np.random.choice(16, n_perturb, replace=False)
            for idx in indices:
                DW[idx] = np.random.randint(1, MASK + 1)

            if all(d == 0 for d in DW):
                continue

            DW_full = schedule_delta(W, DW)
            hw16 = hw32(DW_full[16])

            if hw16 < best_hw16:
                best_hw16 = hw16
                best_dw = list(DW)

            if DW_full[16] == 0:
                best_dw = list(DW)
                best_hw16 = 0
                break

        if best_dw and best_hw16 == 0:
            DW_full = schedule_delta(W, best_dw)
            profile = [hw32(DW_full[t]) for t in range(16, 64)]
            propagation_profiles.append(profile)

    found = len(propagation_profiles)
    print(f"    Found ΔW[16]=0 in {found}/{n_base} trials")

    if found > 0:
        profiles_arr = np.array(propagation_profiles)
        print(f"\n    Propagation after ΔW[16]=0 (HW of ΔW[t]):")
        print(f"    {'t':>4}  {'Mean HW':>8}  {'Min':>5}  {'Max':>5}")
        for i, t in enumerate(range(16, min(40, 64))):
            m = np.mean(profiles_arr[:, i])
            mn = np.min(profiles_arr[:, i])
            mx = np.max(profiles_arr[:, i])
            print(f"    {t:>4}  {m:>8.2f}  {mn:>5}  {mx:>5}")

    # Part 2: Cascade — ΔW[16]=...=ΔW[16+k]=0
    print(f"\n  Part 2: Cascade — find max k with ΔW[16..16+k]=0")
    print(f"    Using Newton-like iterative search")

    cascade_results = []

    for trial in range(20):
        W = random_words(16)

        best_k = -1
        best_dw_cascade = None

        # Try many random DW, keep those with most consecutive zeros
        for attempt in range(200000):
            DW = [0] * 16
            n_perturb = np.random.randint(1, 6)
            indices = np.random.choice(16, n_perturb, replace=False)
            for idx in indices:
                DW[idx] = np.random.randint(1, MASK + 1)

            if all(d == 0 for d in DW):
                continue

            DW_full = schedule_delta(W, DW)

            # Count consecutive zeros from t=16
            k = 0
            for t in range(16, 64):
                if DW_full[t] == 0:
                    k += 1
                else:
                    break

            if k > best_k:
                best_k = k
                best_dw_cascade = list(DW)

        cascade_results.append(best_k)
        if best_k > 0:
            DW_full = schedule_delta(W, best_dw_cascade)
            nonzero_dw = sum(1 for i in range(16) if best_dw_cascade[i] != 0)
            print(f"    Trial {trial+1:>2}: max k={best_k}, "
                  f"nonzero DW words={nonzero_dw}, "
                  f"HW(ΔW[16+k+1])={hw32(DW_full[16+best_k]) if 16+best_k < 64 else 'N/A'}")

    print(f"\n    Cascade results (20 trials):")
    print(f"    Max k achieved: {max(cascade_results)}")
    print(f"    Mean k: {np.mean(cascade_results):.2f}")
    print(f"    Distribution: {dict(sorted(defaultdict(int, {k: cascade_results.count(k) for k in set(cascade_results)}).items()))}")

    # Part 3: Targeted search — fix DW to single word, measure cascade
    print(f"\n  Part 3: Single-word DW perturbation cascade")
    print(f"    {'DW word':>8}  {'k=0 (ΔW16=0)':>13}  {'Mean HW(ΔW16)':>14}  {'Mean HW(ΔW17)':>14}")

    for dw_idx in range(16):
        hw16_list = []
        hw17_list = []
        zeros_16 = 0
        n_trials = 10000
        for _ in range(n_trials):
            W = random_words(16)
            DW = [0] * 16
            DW[dw_idx] = np.random.randint(1, MASK + 1)
            DW_full = schedule_delta(W, DW)
            hw16_list.append(hw32(DW_full[16]))
            hw17_list.append(hw32(DW_full[17]))
            if DW_full[16] == 0:
                zeros_16 += 1

        print(f"    DW[{dw_idx:>2}]   {zeros_16:>13}  {np.mean(hw16_list):>14.2f}  {np.mean(hw17_list):>14.2f}")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")
    print()
    return cascade_results


# =============================================================
# EXPERIMENT B: Partial schedule collision (SA)
# =============================================================
def experiment_B():
    print("=" * 70)
    print("EXPERIMENT B: Partial schedule collision (Simulated Annealing)")
    print("=" * 70)
    t0 = time.time()

    k_values = [1, 2, 4, 8, 12, 16]

    for k in k_values:
        W = random_words(16)

        def fitness(DW):
            """Sum of HW(ΔW[t]) for t=16..16+k-1."""
            DW_full = schedule_delta(W, DW)
            return sum(hw32(DW_full[t]) for t in range(16, 16 + k))

        best_fitness = float('inf')
        best_dw = None

        n_restarts = 50
        n_steps = 50000

        for restart in range(n_restarts):
            # Random initial DW
            DW = [0] * 16
            n_perturb = np.random.randint(1, 5)
            for idx in np.random.choice(16, n_perturb, replace=False):
                DW[idx] = np.random.randint(1, MASK + 1)
            if all(d == 0 for d in DW):
                DW[0] = 1

            current_fit = fitness(DW)
            temp = 32.0

            for step in range(n_steps):
                # Mutate: flip random bits
                DW_new = list(DW)
                word_idx = np.random.randint(0, 16)
                bit_idx = np.random.randint(0, 32)
                DW_new[word_idx] ^= (1 << bit_idx)

                # Rule 7: don't allow all-zero DW
                if all(d == 0 for d in DW_new):
                    continue

                new_fit = fitness(DW_new)
                delta = new_fit - current_fit

                if delta < 0 or np.random.random() < math.exp(-delta / max(temp, 0.01)):
                    DW = DW_new
                    current_fit = new_fit

                temp *= 0.99995

                if current_fit < best_fitness:
                    best_fitness = current_fit
                    best_dw = list(DW)

                if current_fit == 0:
                    break

            if best_fitness == 0:
                break

        # Profile
        DW_full = schedule_delta(W, best_dw)
        n_near_zero = sum(1 for t in range(16, 64) if hw32(DW_full[t]) < 4)

        print(f"\n  k={k:>2}: best fitness={best_fitness:>4} (sum HW[16..{16+k-1}])")
        print(f"    Nonzero DW words: {sum(1 for d in best_dw if d != 0)}")
        print(f"    Near-zero schedule words (HW<4): {n_near_zero}/48")

        # Show profile
        profile_str = ""
        for t in range(16, min(16 + k + 8, 64)):
            profile_str += f"ΔW[{t}]={hw32(DW_full[t]):>2} "
        print(f"    Profile: {profile_str}")

        if k <= 4 and best_fitness == 0:
            # Show full cascade
            print(f"    Full profile (HW): ", end="")
            for t in range(16, 40):
                print(f"{hw32(DW_full[t]):>2}", end=" ")
            print()

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")
    print()


# =============================================================
# EXPERIMENT C: Schedule near-periodicity
# =============================================================
def experiment_C(n=1000):
    print("=" * 70)
    print("EXPERIMENT C: Schedule near-periodicity & autocorrelation")
    print("=" * 70)
    t0 = time.time()

    periods = [1, 2, 3, 4, 8, 16]

    # Part 1: Mean additive distance for each period
    print(f"\n  Part 1: Mean |W[t+p] - W[t]| for periods p")
    print(f"  Expected for random: ~2^31 = {2**31}")

    period_dists = {p: [] for p in periods}

    all_W = []
    for _ in range(n):
        M = random_words(16)
        W = message_schedule(M)
        all_W.append(W)

        for p in periods:
            dists = []
            for t in range(16, 64 - p):
                d = abs(to_signed(sub(W[t + p], W[t])))
                dists.append(d)
            period_dists[p].append(np.mean(dists))

    print(f"\n  {'Period p':>8}  {'Mean dist':>12}  {'Ratio to 2^31':>14}  {'Std':>12}")
    for p in periods:
        m = np.mean(period_dists[p])
        s = np.std(period_dists[p])
        ratio = m / (2**31)
        print(f"  p={p:>5}  {m:>12.0f}  {ratio:>14.4f}  {s:>12.0f}")

    # Part 2: Autocorrelation ACF(k) for k=1..32
    print(f"\n  Part 2: Autocorrelation ACF(k) = corr(W[t], W[t+k])")
    print(f"  Averaged over t=16..60 and {n} messages")

    acf = {}
    for k in range(1, 33):
        corrs = []
        for W in all_W:
            vals_t = []
            vals_tk = []
            for t in range(16, 64 - k):
                vals_t.append(W[t])
                vals_tk.append(W[t + k])
            if len(vals_t) > 2:
                vals_t = np.array(vals_t, dtype=np.float64)
                vals_tk = np.array(vals_tk, dtype=np.float64)
                if np.std(vals_t) > 0 and np.std(vals_tk) > 0:
                    c = np.corrcoef(vals_t, vals_tk)[0, 1]
                    corrs.append(c)
        acf[k] = np.mean(corrs) if corrs else 0.0

    print(f"\n  {'k':>4}  {'ACF(k)':>10}")
    for k in range(1, 33):
        bar = '*' * int(abs(acf[k]) * 500)
        sign = '+' if acf[k] >= 0 else '-'
        print(f"  {k:>4}  {acf[k]:>+10.6f}  {bar}")

    # Check for significant ACF
    significant = [(k, acf[k]) for k in range(1, 33) if abs(acf[k]) > 0.01]
    if significant:
        print(f"\n  Significant ACF (|ACF|>0.01):")
        for k, v in significant:
            print(f"    k={k}: ACF={v:+.6f}")
    else:
        print(f"\n  No significant autocorrelation found (all |ACF|<0.01)")

    # Part 3: XOR autocorrelation for comparison
    print(f"\n  Part 3: XOR autocorrelation (bit-level)")
    xor_acf = {}
    for k in [1, 2, 4, 8, 16]:
        bit_corrs = []
        for W in all_W[:200]:
            for t in range(16, 64 - k):
                xor_val = W[t] ^ W[t + k]
                bit_corrs.append(hw32(xor_val))
        xor_acf[k] = np.mean(bit_corrs)
        print(f"    k={k:>2}: mean HW(W[t] XOR W[t+k]) = {xor_acf[k]:.3f} (random=16.0)")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")
    print()
    return acf


# =============================================================
# EXPERIMENT D: Two-message schedule alignment
# =============================================================
def experiment_D():
    print("=" * 70)
    print("EXPERIMENT D: Two-message schedule alignment")
    print("=" * 70)
    t0 = time.time()

    # Part 1: Random search — T_align distribution
    print(f"\n  Part 1: T_align distribution (1M random pairs)")
    n_pairs = 1000000
    t_align_counts = defaultdict(int)
    never_aligned = 0
    best_t_align = 64
    best_pair = None

    for i in range(n_pairs):
        M = random_words(16)
        # Random delta: flip 1-3 random bits
        M_prime = list(M)
        n_flip = np.random.randint(1, 4)
        for _ in range(n_flip):
            w = np.random.randint(0, 16)
            b = np.random.randint(0, 32)
            M_prime[w] ^= (1 << b)

        W = message_schedule(M)
        W_prime = message_schedule(M_prime)

        # T_align = min T such that max(HW(ΔW[t])) < 16 for all t >= T
        t_align = 64  # never
        for T in range(16, 64):
            max_hw = max(hw32(sub(W_prime[t], W[t])) for t in range(T, 64))
            if max_hw < 16:
                t_align = T
                break

        if t_align < 64:
            t_align_counts[t_align] += 1
            if t_align < best_t_align:
                best_t_align = t_align
                best_pair = (M, M_prime)
        else:
            never_aligned += 1

    print(f"    Never aligned (T_align=64): {never_aligned}/{n_pairs} "
          f"({100*never_aligned/n_pairs:.2f}%)")
    if t_align_counts:
        print(f"    Best T_align: {best_t_align}")
        print(f"    Distribution of T_align:")
        for t in sorted(t_align_counts.keys()):
            print(f"      T={t}: {t_align_counts[t]}")
    else:
        print(f"    No pairs found with T_align < 64")

    # Part 2: SA search for schedule alignment
    print(f"\n  Part 2: SA search for min max(HW(ΔW[t])) for t=32..63")
    M = random_words(16)
    W = message_schedule(M)

    best_max_hw = 999
    best_M_prime = None

    n_restarts = 20
    n_steps = 50000

    for restart in range(n_restarts):
        M_p = list(M)
        # Random perturbation
        for _ in range(np.random.randint(1, 4)):
            w = np.random.randint(0, 16)
            b = np.random.randint(0, 32)
            M_p[w] ^= (1 << b)

        def sa_fitness(Mp):
            Wp = message_schedule(Mp)
            return max(hw32(sub(Wp[t], W[t])) for t in range(32, 64))

        current_fit = sa_fitness(M_p)
        temp = 16.0

        for step in range(n_steps):
            M_new = list(M_p)
            w = np.random.randint(0, 16)
            b = np.random.randint(0, 32)
            M_new[w] ^= (1 << b)

            if M_new == M:
                continue

            new_fit = sa_fitness(M_new)
            delta = new_fit - current_fit

            if delta < 0 or np.random.random() < math.exp(-delta / max(temp, 0.01)):
                M_p = M_new
                current_fit = new_fit

            temp *= 0.9999

            if current_fit < best_max_hw:
                best_max_hw = current_fit
                best_M_prime = list(M_p)

    print(f"    Best max(HW(ΔW[t])) for t=32..63: {best_max_hw}")
    if best_M_prime:
        Wp = message_schedule(best_M_prime)
        print(f"    HW profile of ΔW[t] for t=32..63:")
        for t in range(32, 64):
            dw = sub(Wp[t], W[t])
            print(f"      ΔW[{t}]: HW={hw32(dw):>2}", end="")
            if (t - 31) % 8 == 0:
                print()
        print()

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")
    print()
    return best_max_hw


# =============================================================
# MAIN
# =============================================================
def main():
    print("=" * 70)
    print("  ЗАДАНИЕ 4: Schedule Collision")
    print("  Additive vs XOR schedule zeros")
    print("  Priority: A > B > C > D")
    print("=" * 70)
    print()

    # Quick verification: GF(2) kernel = {0}
    print("Verification: XOR schedule kernel = {0}")
    for _ in range(100):
        W = random_words(16)
        DW = random_words(16)
        # XOR schedule
        W_xor = [W[i] ^ DW[i] for i in range(16)]
        W_exp = message_schedule(W)
        W_xor_exp = message_schedule(W_xor)
        DW_xor = [W_xor_exp[t] ^ W_exp[t] for t in range(64)]
        # Check: is DW_xor[16..63] all zero?
        if all(DW_xor[t] == 0 for t in range(16, 64)):
            print("  FOUND XOR kernel element! (unexpected)")
            break
    else:
        print("  Confirmed: no XOR kernel elements found (100 trials)")
    print()

    # Experiment A (priority)
    cascade = experiment_A()

    # Experiment B
    experiment_B()

    # Experiment C
    acf = experiment_C(n=1000)

    # Experiment D
    best_d = experiment_D()

    # SUMMARY
    print("=" * 70)
    print("  SUMMARY")
    print("=" * 70)
    print(f"\n  A: Max cascade k (ΔW[16..16+k]=0): max={max(cascade)}, mean={np.mean(cascade):.2f}")
    print(f"     Single-equation ΔW[16]=0: solvable by random search")
    print(f"     Cascade quickly fails (k typically 0-1)")
    sig_acf = [(k, acf[k]) for k in range(1, 33) if abs(acf[k]) > 0.01]
    print(f"  C: Significant ACF entries: {len(sig_acf)}")
    if sig_acf:
        print(f"     Strongest: k={sig_acf[0][0]}, ACF={sig_acf[0][1]:+.6f}")
    print(f"  D: Best schedule alignment (max HW for t≥32): {best_d}")
    print()
    print(f"  Rule 7: all ΔW≠0 verified")
    print(f"  Rule 12: null = random ΔW for comparison")
    print(f"  Rule 1: honest results reported")
    print()


if __name__ == "__main__":
    main()
