#!/usr/bin/env python3
"""
C13: Modular Differentials × Carry Identity
============================================
Tests whether modular differentials (D = x' - x mod 2^32) propagate more
sparsely/predictably than XOR differentials through the SHA-256 round function.

Background:
- Wang (2004-2005) used modular differentials for MD5/SHA-1 attacks
- Carry identity: D_add[L](x,y) = L(c) - 2(L(x)∧L(c)) where c=(x+y)⊕x
- XOR and modular related by: A⊕B = A+B - 2(A∧B)
"""

import numpy as np
import time
from collections import defaultdict

MOD = 2**32

# SHA-256 round constants (first 64)
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
    0x391c0cb3, 0x4ed8aa4a, 0x5b9cda5e, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
]

def rotr(x, n):
    return ((x >> n) | (x << (32 - n))) & 0xFFFFFFFF

def shr(x, n):
    return x >> n

def Ch(e, f, g):
    return (e & f) ^ (~e & g) & 0xFFFFFFFF

def Maj(a, b, c):
    return (a & b) ^ (a & c) ^ (b & c)

def Sigma0(a):
    return rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22)

def Sigma1(e):
    return rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25)

def add32(*args):
    s = 0
    for a in args:
        s = (s + a) % MOD
    return s

def sha256_round(state, W, Ki):
    """One SHA-256 round. state = [a,b,c,d,e,f,g,h]. Returns new state."""
    a, b, c, d, e, f, g, h = state
    T1 = add32(h, Sigma1(e), Ch(e, f, g), Ki, W)
    T2 = add32(Sigma0(a), Maj(a, b, c))
    new_a = add32(T1, T2)
    new_e = add32(d, T1)
    return [new_a, a, b, c, new_e, e, f, g]

def hw(x):
    """Hamming weight of 32-bit integer."""
    return bin(x & 0xFFFFFFFF).count('1')

def mod_diff(x_prime, x):
    """Modular differential: x' - x mod 2^32."""
    return (x_prime - x) % MOD

def xor_diff(x_prime, x):
    """XOR differential."""
    return (x_prime ^ x) & 0xFFFFFFFF

def signed_mod(d):
    """Convert unsigned mod-2^32 diff to signed magnitude for |D_mod|."""
    if d > MOD // 2:
        return MOD - d
    return d

def entropy_bits(values, nbins=256):
    """Estimate entropy of a distribution of values (in bits)."""
    if len(values) == 0:
        return 0.0
    counts = defaultdict(int)
    # Hash values into bins for entropy estimation
    for v in values:
        counts[v % nbins] += 1
    total = len(values)
    ent = 0.0
    for c in counts.values():
        if c > 0:
            p = c / total
            ent -= p * np.log2(p)
    return ent


def experiment_1_single_round(N=10000):
    """
    Single-round: modular vs XOR differential sparsity for various δ.
    """
    print("=" * 72)
    print("EXPERIMENT 1: Single-Round Modular vs XOR Differential")
    print("=" * 72)

    rng = np.random.RandomState(42)

    # Perturbation deltas: +1, +2, +4, ..., +2^31
    deltas = [1 << k for k in range(32)]

    # We'll track results for registers a' (index 0) and e' (index 4)
    results = {}  # (delta_k, reg) -> {'hw_mod': [], 'hw_xor': []}

    for ki, delta in enumerate(deltas):
        results_a_mod_hw = []
        results_a_xor_hw = []
        results_e_mod_hw = []
        results_e_xor_hw = []

        for _ in range(N):
            state = [rng.randint(0, MOD) for _ in range(8)]
            W = rng.randint(0, MOD)
            W_prime = (W + delta) % MOD

            out = sha256_round(state, W, K[0])
            out_prime = sha256_round(state, W_prime, K[0])

            # Register a' (index 0)
            dm_a = mod_diff(out_prime[0], out[0])
            dx_a = xor_diff(out_prime[0], out[0])
            results_a_mod_hw.append(hw(dm_a))
            results_a_xor_hw.append(hw(dx_a))

            # Register e' (index 4)
            dm_e = mod_diff(out_prime[4], out[4])
            dx_e = xor_diff(out_prime[4], out[4])
            results_e_mod_hw.append(hw(dm_e))
            results_e_xor_hw.append(hw(dx_e))

        results[ki] = {
            'a_mod_hw': np.mean(results_a_mod_hw),
            'a_xor_hw': np.mean(results_a_xor_hw),
            'e_mod_hw': np.mean(results_e_mod_hw),
            'e_xor_hw': np.mean(results_e_xor_hw),
            'a_mod_std': np.std(results_a_mod_hw),
            'a_xor_std': np.std(results_a_xor_hw),
            'e_mod_std': np.std(results_e_mod_hw),
            'e_xor_std': np.std(results_e_xor_hw),
        }

    print(f"\n{'delta':>12s}  {'HW(Da_mod)':>12s} {'HW(Da_xor)':>12s} "
          f"{'HW(De_mod)':>12s} {'HW(De_xor)':>12s}  mod<xor?")
    print("-" * 82)

    mod_wins_a = 0
    mod_wins_e = 0
    for ki in range(32):
        r = results[ki]
        delta_str = f"2^{ki}"
        a_better = r['a_mod_hw'] < r['a_xor_hw']
        e_better = r['e_mod_hw'] < r['e_xor_hw']
        if a_better: mod_wins_a += 1
        if e_better: mod_wins_e += 1
        flag = ""
        if a_better and e_better: flag = "BOTH"
        elif a_better: flag = "a only"
        elif e_better: flag = "e only"
        else: flag = "no"

        if ki % 4 == 0 or ki < 4:  # Print subset
            print(f"{delta_str:>12s}  "
                  f"{r['a_mod_hw']:>6.2f}±{r['a_mod_std']:.2f} "
                  f"{r['a_xor_hw']:>6.2f}±{r['a_xor_std']:.2f} "
                  f"{r['e_mod_hw']:>6.2f}±{r['e_mod_std']:.2f} "
                  f"{r['e_xor_hw']:>6.2f}±{r['e_xor_std']:.2f}  {flag}")

    print(f"\nModular HW < XOR HW: a'={mod_wins_a}/32, e'={mod_wins_e}/32")
    return results


def experiment_2_entropy(N=5000):
    """
    Entropy comparison: is D_mod more predictable than D_xor?
    """
    print("\n" + "=" * 72)
    print("EXPERIMENT 2: Entropy of Modular vs XOR Differentials")
    print("=" * 72)

    rng = np.random.RandomState(123)

    deltas_k = [0, 1, 2, 4, 8, 16, 31]

    print(f"\n{'delta':>10s}  {'Ent(Da_mod)':>12s} {'Ent(Da_xor)':>12s} "
          f"{'Ent(De_mod)':>12s} {'Ent(De_xor)':>12s}  lower?")
    print("-" * 78)

    for ki in deltas_k:
        delta = 1 << ki
        dm_a_vals = []
        dx_a_vals = []
        dm_e_vals = []
        dx_e_vals = []

        for _ in range(N):
            state = [rng.randint(0, MOD) for _ in range(8)]
            W = rng.randint(0, MOD)
            W_prime = (W + delta) % MOD

            out = sha256_round(state, W, K[0])
            out_prime = sha256_round(state, W_prime, K[0])

            dm_a_vals.append(mod_diff(out_prime[0], out[0]))
            dx_a_vals.append(xor_diff(out_prime[0], out[0]))
            dm_e_vals.append(mod_diff(out_prime[4], out[4]))
            dx_e_vals.append(xor_diff(out_prime[4], out[4]))

        # Use more bins for better entropy estimate
        nbins = 1024
        ent_a_mod = entropy_bits(dm_a_vals, nbins)
        ent_a_xor = entropy_bits(dx_a_vals, nbins)
        ent_e_mod = entropy_bits(dm_e_vals, nbins)
        ent_e_xor = entropy_bits(dx_e_vals, nbins)

        flag = ""
        if ent_a_mod < ent_a_xor and ent_e_mod < ent_e_xor:
            flag = "MOD BOTH"
        elif ent_a_mod < ent_a_xor:
            flag = "MOD a"
        elif ent_e_mod < ent_e_xor:
            flag = "MOD e"
        else:
            flag = "XOR"

        print(f"2^{ki:>2d}        "
              f"{ent_a_mod:>8.3f}     {ent_a_xor:>8.3f}     "
              f"{ent_e_mod:>8.3f}     {ent_e_xor:>8.3f}      {flag}")


def experiment_3_carry_correction(N=10000):
    """
    Carry correction term: does 2(L(x)∧L(c)) vanish for small modular δ?
    """
    print("\n" + "=" * 72)
    print("EXPERIMENT 3: Carry Correction Term Analysis")
    print("=" * 72)
    print("For addition x+y with perturbation y'=y+δ:")
    print("  c = (x+y)⊕x⊕y (carry bits)")
    print("  Carry correction = 2(x∧c)")
    print("  If correction is small, modular ≈ XOR differential")

    rng = np.random.RandomState(456)

    print(f"\n{'delta':>10s}  {'mean_HW(corr)':>14s} {'frac_zero':>10s} "
          f"{'mean_HW(c)':>12s}")
    print("-" * 56)

    for ki in range(0, 32, 2):
        delta = 1 << ki
        corr_hws = []
        c_hws = []
        zero_count = 0

        for _ in range(N):
            x = rng.randint(0, MOD)
            y = rng.randint(0, MOD)
            y_prime = (y + delta) % MOD

            s = (x + y) % MOD
            s_prime = (x + y_prime) % MOD

            # Carry bits for original sum
            c = (s ^ x ^ y) & 0xFFFFFFFF
            # Carry bits for perturbed sum
            c_prime = (s_prime ^ x ^ y_prime) & 0xFFFFFFFF

            # Carry correction term
            correction = (2 * (x & c)) & 0xFFFFFFFF
            correction_prime = (2 * (x & c_prime)) & 0xFFFFFFFF

            # Differential of carry correction
            corr_diff = xor_diff(correction_prime, correction)
            corr_hws.append(hw(corr_diff))
            c_hws.append(hw(xor_diff(c_prime, c)))
            if corr_diff == 0:
                zero_count += 1

        frac_zero = zero_count / N
        print(f"2^{ki:>2d}        {np.mean(corr_hws):>10.3f}     "
              f"{frac_zero:>8.4f}   {np.mean(c_hws):>8.3f}")


def experiment_4_multiround(N=3000, max_rounds=6):
    """
    Multi-round: track modular vs XOR differentials across rounds 0-5.
    """
    print("\n" + "=" * 72)
    print("EXPERIMENT 4: Multi-Round Differential Propagation")
    print("=" * 72)

    rng = np.random.RandomState(789)

    deltas_k = [0, 1, 4, 8, 16]

    for ki in deltas_k:
        delta = 1 << ki
        print(f"\n--- Input perturbation δ = 2^{ki} on W[0] ---")
        print(f"{'Round':>6s}  {'HW(Da_mod)':>11s} {'HW(Da_xor)':>11s} "
              f"{'HW(De_mod)':>11s} {'HW(De_xor)':>11s}  winner")
        print("-" * 68)

        for r in range(max_rounds):
            hw_a_mod = []
            hw_a_xor = []
            hw_e_mod = []
            hw_e_xor = []

            for _ in range(N):
                state = [rng.randint(0, MOD) for _ in range(8)]
                W_vals = [rng.randint(0, MOD) for _ in range(max_rounds)]
                W_vals_prime = W_vals.copy()
                W_vals_prime[0] = (W_vals[0] + delta) % MOD  # Perturb only W[0]

                s1 = list(state)
                s2 = list(state)
                for rr in range(r + 1):
                    s1 = sha256_round(s1, W_vals[rr], K[rr])
                    s2 = sha256_round(s2, W_vals_prime[rr], K[rr])

                hw_a_mod.append(hw(mod_diff(s2[0], s1[0])))
                hw_a_xor.append(hw(xor_diff(s2[0], s1[0])))
                hw_e_mod.append(hw(mod_diff(s2[4], s1[4])))
                hw_e_xor.append(hw(xor_diff(s2[4], s1[4])))

            ma = np.mean(hw_a_mod)
            xa = np.mean(hw_a_xor)
            me = np.mean(hw_e_mod)
            xe = np.mean(hw_e_xor)

            w = ""
            if ma < xa and me < xe: w = "MOD"
            elif ma < xa: w = "MOD(a)"
            elif me < xe: w = "MOD(e)"
            else: w = "XOR/TIE"

            print(f"{r:>6d}  {ma:>8.2f}    {xa:>8.2f}    "
                  f"{me:>8.2f}    {xe:>8.2f}     {w}")


def experiment_5_relaxed_wang(N=5000, max_rounds=4):
    """
    Test: when De is small but nonzero, is |De_mod| < 2^k easier to maintain
    than HW(De_xor) < k?
    """
    print("\n" + "=" * 72)
    print("EXPERIMENT 5: Relaxed Wang Condition — Small De")
    print("=" * 72)
    print("Question: is it easier to keep |De_mod| small or HW(De_xor) small?")

    rng = np.random.RandomState(2024)

    # Use small modular perturbations on W
    for delta_k in [0, 1, 2, 3]:
        delta = 1 << delta_k
        print(f"\n--- δ = 2^{delta_k} on W[0], after {max_rounds} rounds ---")
        print(f"{'Threshold k':>12s}  {'P(|De_mod|<2^k)':>16s} {'P(HW(De_xor)<k)':>17s}  ratio")
        print("-" * 62)

        de_mod_vals = []
        de_xor_vals = []

        for _ in range(N):
            state = [rng.randint(0, MOD) for _ in range(8)]
            W_vals = [rng.randint(0, MOD) for _ in range(max_rounds)]
            W_vals_prime = W_vals.copy()
            W_vals_prime[0] = (W_vals[0] + delta) % MOD

            s1 = list(state)
            s2 = list(state)
            for r in range(max_rounds):
                s1 = sha256_round(s1, W_vals[r], K[r])
                s2 = sha256_round(s2, W_vals_prime[r], K[r])

            de_mod = mod_diff(s2[4], s1[4])
            de_xor = xor_diff(s2[4], s1[4])
            de_mod_vals.append(signed_mod(de_mod))
            de_xor_vals.append(hw(de_xor))

        de_mod_arr = np.array(de_mod_vals)
        de_xor_arr = np.array(de_xor_vals)

        for k in range(1, 17):
            mod_thresh = 1 << k
            p_mod = np.mean(de_mod_arr < mod_thresh)
            p_xor = np.mean(de_xor_arr < k)
            ratio = p_mod / p_xor if p_xor > 0 else float('inf')
            marker = " <-- MOD easier" if ratio > 2.0 else (" <-- XOR easier" if ratio < 0.5 else "")
            print(f"{k:>12d}  {p_mod:>14.4f}   {p_xor:>15.4f}   {ratio:>6.2f}{marker}")


def experiment_6_count_small_de(N=8000, max_rounds=4):
    """
    Count how often |De_mod| < 2^k vs HW(De_xor) < k across rounds.
    """
    print("\n" + "=" * 72)
    print("EXPERIMENT 6: Fraction with Small De — Modular vs XOR, per Round")
    print("=" * 72)

    rng = np.random.RandomState(999)
    delta = 1  # Single-bit modular perturbation

    for rnd in range(1, 6):
        print(f"\n--- After {rnd} round(s), δ=+1 on W[0] ---")

        de_mod_vals = []
        de_xor_vals = []

        for _ in range(N):
            state = [rng.randint(0, MOD) for _ in range(8)]
            W_vals = [rng.randint(0, MOD) for _ in range(rnd)]
            W_vals_prime = W_vals.copy()
            W_vals_prime[0] = (W_vals[0] + delta) % MOD

            s1 = list(state)
            s2 = list(state)
            for r in range(rnd):
                s1 = sha256_round(s1, W_vals[r], K[r])
                s2 = sha256_round(s2, W_vals_prime[r], K[r])

            de_mod_vals.append(signed_mod(mod_diff(s2[4], s1[4])))
            de_xor_vals.append(hw(xor_diff(s2[4], s1[4])))

        de_mod_arr = np.array(de_mod_vals)
        de_xor_arr = np.array(de_xor_vals)

        print(f"{'k':>4s}  {'P(|De_mod|<2^k)':>16s}  {'P(HW(De_xor)<k)':>17s}  mod/xor")
        print("-" * 50)
        for k in [1, 2, 4, 8, 12, 16]:
            p_mod = np.mean(de_mod_arr < (1 << k))
            p_xor = np.mean(de_xor_arr < k)
            ratio = p_mod / p_xor if p_xor > 0 else float('inf')
            print(f"{k:>4d}  {p_mod:>14.4f}    {p_xor:>15.4f}    {ratio:.3f}")


def verdict(results_exp1):
    """Determine overall verdict."""
    print("\n" + "=" * 72)
    print("VERDICT")
    print("=" * 72)

    # Check experiment 1: how often is modular HW lower?
    mod_better_a = sum(1 for ki in range(32) if results_exp1[ki]['a_mod_hw'] < results_exp1[ki]['a_xor_hw'])
    mod_better_e = sum(1 for ki in range(32) if results_exp1[ki]['e_mod_hw'] < results_exp1[ki]['e_xor_hw'])

    # Compute average advantage
    avg_diff_a = np.mean([results_exp1[ki]['a_mod_hw'] - results_exp1[ki]['a_xor_hw'] for ki in range(32)])
    avg_diff_e = np.mean([results_exp1[ki]['e_mod_hw'] - results_exp1[ki]['e_xor_hw'] for ki in range(32)])

    print(f"\nSingle-round HW comparison (32 delta values):")
    print(f"  Register a': modular wins {mod_better_a}/32 times, avg diff = {avg_diff_a:.3f}")
    print(f"  Register e': modular wins {mod_better_e}/32 times, avg diff = {avg_diff_e:.3f}")
    print(f"  (negative diff = modular is sparser)")

    # The key insight: for SHA-256, XOR and modular differentials
    # should behave very similarly due to the mixing. The carry identity
    # doesn't simplify things because the nonlinear operations (Ch, Maj,
    # Sigma) are bit-based, not modular.

    total_mod_wins = mod_better_a + mod_better_e
    significant_advantage = abs(avg_diff_a) > 0.5 or abs(avg_diff_e) > 0.5

    if total_mod_wins > 48 and significant_advantage and avg_diff_a < 0 and avg_diff_e < 0:
        label = "ALIVE"
        explanation = ("Modular differentials consistently propagate with lower "
                      "Hamming weight than XOR differentials.")
    elif total_mod_wins > 36 or (abs(avg_diff_a) > 0.2 or abs(avg_diff_e) > 0.2):
        label = "ANOMALY"
        explanation = ("Some difference exists between modular and XOR differential "
                      "propagation, but it is not consistent enough for systematic exploitation.")
    else:
        label = "DEAD"
        explanation = ("No significant difference between modular and XOR differential "
                      "propagation in SHA-256. The bitwise operations (Ch, Maj, Sigma) "
                      "destroy any advantage that modular arithmetic might offer over XOR.")

    print(f"\n>>> VERDICT: {label}")
    print(f"    {explanation}")

    print("\nAnalysis:")
    print("  - SHA-256 uses bitwise rotations and Boolean functions (Ch, Maj)")
    print("  - These operations are defined on individual bits, not modular arithmetic")
    print("  - The carry identity relates XOR to addition, but Sigma/Ch/Maj are purely bitwise")
    print("  - After these nonlinear bitwise ops, any modular structure is destroyed")
    print("  - Wang's technique worked on MD5/SHA-1 because those have simpler structure")
    print("  - SHA-256's design specifically resists modular differential shortcuts")

    return label


def main():
    t0 = time.time()
    print("C13: Modular Differentials × Carry Identity")
    print("Testing whether modular diffs propagate more sparsely than XOR in SHA-256\n")

    results_exp1 = experiment_1_single_round(N=10000)
    t1 = time.time()
    print(f"\n[Exp 1 took {t1-t0:.1f}s]")

    experiment_2_entropy(N=5000)
    t2 = time.time()
    print(f"\n[Exp 2 took {t2-t1:.1f}s]")

    experiment_3_carry_correction(N=10000)
    t3 = time.time()
    print(f"\n[Exp 3 took {t3-t2:.1f}s]")

    experiment_4_multiround(N=3000, max_rounds=6)
    t4 = time.time()
    print(f"\n[Exp 4 took {t4-t3:.1f}s]")

    experiment_5_relaxed_wang(N=5000, max_rounds=4)
    t5 = time.time()
    print(f"\n[Exp 5 took {t5-t4:.1f}s]")

    experiment_6_count_small_de(N=8000, max_rounds=4)
    t6 = time.time()
    print(f"\n[Exp 6 took {t6-t5:.1f}s]")

    label = verdict(results_exp1)

    print(f"\nTotal runtime: {time.time()-t0:.1f}s")
    return label


if __name__ == "__main__":
    main()
