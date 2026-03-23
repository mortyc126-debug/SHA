#!/usr/bin/env python3
"""
weapon_cube_v2.py — Advanced Cube Attack on Reduced SHA-256

Pushes cube attack beyond the 6-round key-recovery boundary.
Phase 1: Higher-dimension cubes (dim 6-12) for rounds 7-10
Phase 2: Multi-word cubes (bits spread across W[0] AND W[1])
Phase 3: Optimal bit selection using Sigma1-related positions
Phase 4: Conditional cube attack with favorable bit fixing

Key insight for correct methodology:
  A cube attack recovers a secret bit if the cube sum (superpoly evaluated
  at the secret) is a NON-CONSTANT function of the secret bit. For a single
  output bit, cube_sum is 0 or 1. If the superpoly truly depends on the
  secret, then for a FIXED cube subset and FIXED non-cube/non-secret bits,
  cube_sum(secret=0) != cube_sum(secret=1). We test this across multiple
  random settings of the remaining message words to get statistical power.

  Null hypothesis: superpoly is independent of secret -> P(s0!=s1) = 50%
  True key recovery: superpoly depends on secret -> P(s0!=s1) = 100%
  We require s0!=s1 in ALL of N_CONFIRM trials to declare success.
"""

import random
import time
from collections import defaultdict

# ── SHA-256 constants ──────────────────────────────────────────────
MASK = 0xFFFFFFFF

K = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
]

IV = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
]

# ── Primitives ─────────────────────────────────────────────────────
def rotr(x, n):
    return ((x >> n) | (x << (32 - n))) & MASK

def Sig0(a):
    return rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22)

def Sig1(e):
    return rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25)

def Ch(e, f, g):
    return ((e & f) ^ (~e & g)) & MASK

def Maj(a, b, c):
    return ((a & b) ^ (a & c) ^ (b & c)) & MASK

def add(*args):
    s = 0
    for x in args:
        s = (s + x) & MASK
    return s

def sha256_reduced(W, n_rounds):
    """Run n_rounds of SHA-256 compression, return state [a..h]."""
    a, b, c, d, e, f, g, h = IV
    for i in range(n_rounds):
        S1 = Sig1(e)
        ch = Ch(e, f, g)
        temp1 = add(h, S1, ch, K[i], W[i])
        S0 = Sig0(a)
        maj = Maj(a, b, c)
        temp2 = add(S0, maj)
        h, g, f, e, d, c, b, a = g, f, e, add(d, temp1), c, b, a, add(temp1, temp2)
    return [a, b, c, d, e, f, g, h]


def output_bit(W, n_rounds, bit):
    """Extract single bit from 'a' register after n_rounds."""
    state = sha256_reduced(W, n_rounds)
    return (state[0] >> bit) & 1


# ── Cube attack core ──────────────────────────────────────────────
def cube_sum(W_template, free_positions, n_rounds, out_bit):
    """
    Compute cube sum: XOR of output bit over all 2^dim assignments
    to free_positions. W_template already has secret set.
    """
    dim = len(free_positions)
    total = 0
    for mask in range(1 << dim):
        Wc = list(W_template)
        for idx in range(dim):
            wi, bi = free_positions[idx]
            if mask & (1 << idx):
                Wc[wi] |= (1 << bi)
            else:
                Wc[wi] &= ~(1 << bi)
        total ^= output_bit(Wc, n_rounds, out_bit)
    return total


def test_cube_robust(free_positions, secret_word, secret_bit,
                     n_rounds, out_bit, n_confirm=8):
    """
    Robust cube attack test: for a FIXED cube subset, test across
    n_confirm random settings of all OTHER message words.

    For each random setting:
      - compute cube_sum with secret=0 and secret=1
      - check if they differ

    If secret bit truly appears linearly in the superpoly, ALL n_confirm
    trials should show s0 != s1. Under null hypothesis (random), probability
    of all-differ = (1/2)^n_confirm.

    Returns (n_differ, n_confirm) where n_differ = count of trials with s0!=s1.
    """
    n_differ = 0
    for _ in range(n_confirm):
        W = [random.getrandbits(32) for _ in range(16)]

        # Set secret = 0
        W0 = list(W)
        W0[secret_word] &= ~(1 << secret_bit)
        s0 = cube_sum(W0, free_positions, n_rounds, out_bit)

        # Set secret = 1
        W1 = list(W)
        W1[secret_word] |= (1 << secret_bit)
        s1 = cube_sum(W1, free_positions, n_rounds, out_bit)

        if s0 != s1:
            n_differ += 1

    return n_differ, n_confirm


# ======================================================================
#  PHASE 0: Baseline verification (rounds 5-6 should be ~100%)
# ======================================================================
def phase0(t0, time_budget):
    print("=" * 72)
    print("PHASE 0: Baseline verification")
    print("=" * 72)
    print("  Confirming 100% key-recovery at rounds 5-6 with dim 4-6")
    print("  N_CONFIRM=8 per cube (null prob = 1/256)")
    print()

    secret_word, secret_bit = 1, 0
    out_bit = 0
    n_confirm = 8

    for rnd in [5, 6]:
        for dim in [4, 5, 6]:
            n_full_recovery = 0
            n_tested = 10
            for _ in range(n_tested):
                bits = random.sample(range(32), dim)
                free_pos = [(0, b) for b in bits]
                nd, nc = test_cube_robust(free_pos, secret_word, secret_bit,
                                          rnd, out_bit, n_confirm)
                if nd == nc:
                    n_full_recovery += 1
            print(f"  round={rnd} dim={dim}: {n_full_recovery}/{n_tested} "
                  f"perfect key-recovery (all {n_confirm}/{n_confirm} differ)")

    print()


# ======================================================================
#  PHASE 1: Higher-dimension cubes for rounds 7-10
# ======================================================================
def phase1(t0, time_budget):
    print("=" * 72)
    print("PHASE 1: Higher-dimension cubes — rounds 7-10")
    print("=" * 72)
    print("Free bits in W[0], secret = W[1] bit 0, output = a[0]")
    print("N_CONFIRM=8 per cube (null false-positive prob = 0.39%)")
    print()

    dims = [6, 8, 10, 12]
    rounds_list = [7, 8, 9, 10]
    secret_word, secret_bit = 1, 0
    out_bit = 0
    n_confirm = 8

    # results[rounds][dim] = (perfect_recovery, partial_signal, total)
    results = defaultdict(lambda: defaultdict(lambda: [0, 0, 0]))

    for dim in dims:
        n_subsets = 20 if dim >= 12 else (30 if dim >= 10 else 50)
        for rnd in rounds_list:
            elapsed = time.time() - t0
            if elapsed > time_budget * 0.25:
                n_subsets = min(n_subsets, 10)
            if elapsed > time_budget * 0.32:
                n_subsets = min(n_subsets, 5)
            if elapsed > time_budget * 0.38:
                break

            perfect = 0
            partial = 0
            tested = 0
            for _ in range(n_subsets):
                if time.time() - t0 > time_budget * 0.40:
                    break
                bits = random.sample(range(32), dim)
                free_pos = [(0, b) for b in bits]
                nd, nc = test_cube_robust(free_pos, secret_word, secret_bit,
                                          rnd, out_bit, n_confirm)
                if nd == nc:
                    perfect += 1
                elif nd > nc * 0.75:
                    partial += 1
                tested += 1

            results[rnd][dim] = [perfect, partial, tested]
            if tested > 0:
                pct_p = 100 * perfect / tested
                pct_s = 100 * (perfect + partial) / tested
                print(f"  rounds={rnd:2d}  dim={dim:2d}  "
                      f"perfect: {perfect:3d}/{tested:3d} ({pct_p:5.1f}%)  "
                      f"signal(>75%): {perfect+partial:3d}/{tested:3d} ({pct_s:5.1f}%)")

    # Summary table
    print()
    print("  Summary — perfect key-recovery rate:")
    print(f"  {'':8s}", end="")
    for dim in dims:
        print(f"  dim={dim:2d} ", end="")
    print()
    for rnd in rounds_list:
        print(f"  rnd={rnd:2d}  ", end="")
        for dim in dims:
            p, s, t = results[rnd][dim]
            if t > 0:
                print(f"  {100*p/t:5.1f}%", end="")
            else:
                print(f"     -- ", end="")
        print()

    # Determine failure round
    print()
    for dim in dims:
        last_good = None
        for rnd in rounds_list:
            p, s, t = results[rnd][dim]
            if t > 0 and p / t > 0.1:  # at least 10% perfect recovery
                last_good = rnd
        if last_good is not None:
            p, s, t = results[last_good][dim]
            print(f"  dim={dim:2d}: last round with >10% perfect recovery: "
                  f"round {last_good} ({100*p/t:.0f}%)")
        else:
            print(f"  dim={dim:2d}: no significant key-recovery at rounds 7+")

    return results


# ======================================================================
#  PHASE 2: Multi-word cube (bits in W[0] AND W[1])
# ======================================================================
def phase2(t0, time_budget):
    print()
    print("=" * 72)
    print("PHASE 2: Multi-word cube — 4 bits W[0] + 4 bits W[1]")
    print("=" * 72)
    print("Dim 8 total, secret = W[1] bit 0 (excluded from free bits)")
    print("Compare: 8 bits all in W[0] vs 4+4 split")
    print()

    rounds_list = [7, 8, 9, 10]
    secret_word, secret_bit = 1, 0
    out_bit = 0
    n_confirm = 8
    n_subsets = 20

    results_single = {}
    results_multi = {}

    for rnd in rounds_list:
        if time.time() - t0 > time_budget * 0.58:
            break

        # Single-word: 8 bits in W[0]
        perf_s, tested_s = 0, 0
        for _ in range(n_subsets):
            if time.time() - t0 > time_budget * 0.52:
                break
            bits = random.sample(range(32), 8)
            free_pos = [(0, b) for b in bits]
            nd, nc = test_cube_robust(free_pos, secret_word, secret_bit,
                                      rnd, out_bit, n_confirm)
            if nd == nc:
                perf_s += 1
            tested_s += 1
        results_single[rnd] = (perf_s, tested_s)

        # Multi-word: 4 bits in W[0] + 4 bits in W[1] (avoiding secret bit)
        perf_m, tested_m = 0, 0
        w1_candidates = [b for b in range(32) if b != secret_bit]
        for _ in range(n_subsets):
            if time.time() - t0 > time_budget * 0.56:
                break
            bits_w0 = random.sample(range(32), 4)
            bits_w1 = random.sample(w1_candidates, 4)
            free_pos = [(0, b) for b in bits_w0] + [(1, b) for b in bits_w1]
            nd, nc = test_cube_robust(free_pos, secret_word, secret_bit,
                                      rnd, out_bit, n_confirm)
            if nd == nc:
                perf_m += 1
            tested_m += 1
        results_multi[rnd] = (perf_m, tested_m)

        def pct(s, t):
            return f"{100*s/t:5.1f}%" if t > 0 else "  -- "

        print(f"  rounds={rnd:2d}  single-word: {pct(perf_s, tested_s):>7s}  "
              f"multi-word: {pct(perf_m, tested_m):>7s}")

    print()
    improvement = False
    for rnd in rounds_list:
        if rnd in results_single and rnd in results_multi:
            ss, ts = results_single[rnd]
            sm, tm = results_multi[rnd]
            if ts > 0 and tm > 0:
                if sm / tm > ss / ts + 0.05:
                    improvement = True
    print(f"  Verdict: Multi-word spreading "
          f"{'HELPS' if improvement else 'shows NO significant improvement'}")

    return results_single, results_multi


# ======================================================================
#  PHASE 3: Optimal bit selection (Sigma1-aligned positions)
# ======================================================================
def phase3(t0, time_budget):
    print()
    print("=" * 72)
    print("PHASE 3: Optimal bit selection — Sigma1-aligned cube")
    print("=" * 72)
    print("Sigma1(bit31) -> bits {5, 10, 24}")
    print("Optimal cube: free bits at {5, 10, 24, 31} in W[0]  (dim=4)")
    print("Also test dim=8 with {5,6,10,11,20,24,25,31}")
    print("Compare with random subsets of same dimension")
    print()

    rounds_list = [6, 7, 8, 9, 10]
    secret_word, secret_bit = 1, 0
    out_bit = 0
    n_confirm = 8
    n_subsets = 25

    # Optimal dim-4: bit 31 + Sigma1-related
    optimal_4 = [5, 10, 24, 31]
    # Optimal dim-8: bit 31 + neighbors of Sigma1-related positions
    optimal_8 = [5, 6, 10, 11, 20, 24, 25, 31]

    configs = [
        ("Sig1-opt-d4", [(0, b) for b in optimal_4]),
        ("Sig1-opt-d8", [(0, b) for b in optimal_8]),
    ]

    for label, optimal_pos in configs:
        dim = len(optimal_pos)
        print(f"  [{label}] dim={dim}, bits={[b for _, b in optimal_pos]}")

        for rnd in rounds_list:
            if time.time() - t0 > time_budget * 0.76:
                break

            # Optimal selection (same subset, different random messages)
            perf_o, tested_o = 0, 0
            for _ in range(n_subsets):
                if time.time() - t0 > time_budget * 0.72:
                    break
                nd, nc = test_cube_robust(optimal_pos, secret_word, secret_bit,
                                          rnd, out_bit, n_confirm)
                if nd == nc:
                    perf_o += 1
                tested_o += 1

            # Random subsets of same dimension
            perf_r, tested_r = 0, 0
            for _ in range(n_subsets):
                if time.time() - t0 > time_budget * 0.74:
                    break
                bits = random.sample(range(32), dim)
                free_pos = [(0, b) for b in bits]
                nd, nc = test_cube_robust(free_pos, secret_word, secret_bit,
                                          rnd, out_bit, n_confirm)
                if nd == nc:
                    perf_r += 1
                tested_r += 1

            def pct(s, t):
                return f"{100*s/t:5.1f}%" if t > 0 else "  -- "

            print(f"    round={rnd:2d}  optimal: {pct(perf_o, tested_o):>7s}  "
                  f"random: {pct(perf_r, tested_r):>7s}")

        print()

    return {}


# ======================================================================
#  PHASE 4: Conditional cube attack
# ======================================================================
def phase4(t0, time_budget):
    print()
    print("=" * 72)
    print("PHASE 4: Conditional cube attack")
    print("=" * 72)
    print("Fix W[0] upper byte to favorable patterns, cube over bits 0-7")
    print("Compare conditioned vs unconditioned (dim 8)")
    print()

    rounds_list = [8, 9, 10]
    secret_word, secret_bit = 1, 0
    out_bit = 0
    n_confirm = 8
    n_subsets = 20

    conditions = [
        ("zeros_hi", 0x00000000, 0xFF000000),
        ("ones_hi",  0xFF000000, 0xFF000000),
        ("alt_hi",   0xAA000000, 0xFF000000),
        ("none",     0x00000000, 0x00000000),
    ]

    free_pos = [(0, b) for b in range(8)]  # bits 0-7 of W[0]

    results = {}

    for cond_name, cond_val, cond_mask in conditions:
        results[cond_name] = {}
        for rnd in rounds_list:
            if time.time() - t0 > time_budget * 0.95:
                break

            perf, tested = 0, 0
            for _ in range(n_subsets):
                if time.time() - t0 > time_budget * 0.93:
                    break

                # For conditioned test: generate random base, apply condition,
                # then test with robust method
                n_differ = 0
                for __ in range(n_confirm):
                    W = [random.getrandbits(32) for _ in range(16)]
                    W[0] = (W[0] & ~cond_mask) | (cond_val & cond_mask)

                    W0 = list(W)
                    W0[secret_word] &= ~(1 << secret_bit)
                    s0 = cube_sum(W0, free_pos, rnd, out_bit)

                    W1 = list(W)
                    W1[secret_word] |= (1 << secret_bit)
                    s1 = cube_sum(W1, free_pos, rnd, out_bit)

                    if s0 != s1:
                        n_differ += 1

                if n_differ == n_confirm:
                    perf += 1
                tested += 1

            results[cond_name][rnd] = (perf, tested)

    # Print results
    print(f"  {'Condition':<12s}", end="")
    for rnd in rounds_list:
        print(f"  round={rnd}", end="")
    print()
    print("  " + "-" * 50)

    for cond_name, _, _ in conditions:
        print(f"  {cond_name:<12s}", end="")
        for rnd in rounds_list:
            if rnd in results.get(cond_name, {}):
                s, t = results[cond_name][rnd]
                if t > 0:
                    print(f"  {100*s/t:6.1f}%", end="")
                else:
                    print(f"      -- ", end="")
            else:
                print(f"      -- ", end="")
        print()

    print()
    cond_helps = False
    for rnd in rounds_list:
        none_s, none_t = results.get("none", {}).get(rnd, (0, 0))
        if none_t == 0:
            continue
        none_rate = none_s / none_t
        for cn in ["zeros_hi", "ones_hi", "alt_hi"]:
            cs, ct = results.get(cn, {}).get(rnd, (0, 0))
            if ct > 0 and cs / ct > none_rate + 0.1:
                cond_helps = True
    print(f"  Verdict: Conditioning "
          f"{'IMPROVES cube attack reach' if cond_helps else 'shows NO significant improvement'}")

    return results


# ======================================================================
#  MAIN
# ======================================================================
def main():
    random.seed(0xC0BEA772)
    t0 = time.time()
    TIME_BUDGET = 115  # seconds

    print("*" * 72)
    print("  WEAPON: Advanced Cube Attack on Reduced SHA-256 (v2)")
    print("  Baseline: 100% key-recovery through 6 rounds at dim 4-6")
    print("  Goal: Push cube attack boundary further")
    print("  Methodology: N_CONFIRM=8 trials per cube subset")
    print("    (false positive prob under null = (1/2)^8 = 0.39%)")
    print("*" * 72)
    print()

    # Phase 0: baseline
    phase0(t0, TIME_BUDGET)

    # Phase 1
    p1_results = phase1(t0, TIME_BUDGET)

    # Phase 2
    p2_single, p2_multi = phase2(t0, TIME_BUDGET)

    # Phase 3
    phase3(t0, TIME_BUDGET)

    # Phase 4
    p4_results = phase4(t0, TIME_BUDGET)

    # ── Final Summary ─────────────────────────────────────────────
    elapsed = time.time() - t0
    print()
    print("=" * 72)
    print("FINAL SUMMARY")
    print("=" * 72)

    # Determine max round with significant perfect recovery
    max_round_alive = 6  # baseline
    for rnd in [7, 8, 9, 10]:
        for dim in [6, 8, 10, 12]:
            p, s, t = p1_results.get(rnd, {}).get(dim, [0, 0, 0])
            if t > 0 and p / t > 0.05:  # >5% perfect recovery
                max_round_alive = max(max_round_alive, rnd)
        if rnd in p2_multi:
            pm, tm = p2_multi[rnd]
            if tm > 0 and pm / tm > 0.05:
                max_round_alive = max(max_round_alive, rnd)
        if rnd in p2_single:
            ps, ts = p2_single[rnd]
            if ts > 0 and ps / ts > 0.05:
                max_round_alive = max(max_round_alive, rnd)

    print(f"\n  Cube attack key-recovery boundary: round {max_round_alive}")
    print(f"  (>5% perfect key-recovery with robust 8-trial confirmation)")

    # First round with zero signal
    first_dead = None
    for rnd in [7, 8, 9, 10]:
        all_zero = True
        for dim in [6, 8, 10, 12]:
            p, s, t = p1_results.get(rnd, {}).get(dim, [0, 0, 0])
            if t > 0 and p > 0:
                all_zero = False
        if rnd in p2_multi:
            pm, tm = p2_multi[rnd]
            if tm > 0 and pm > 0:
                all_zero = False
        if rnd in p2_single:
            ps, ts = p2_single[rnd]
            if ts > 0 and ps > 0:
                all_zero = False
        if all_zero and first_dead is None:
            first_dead = rnd

    if first_dead:
        print(f"  Complete failure (0% recovery) starts at round {first_dead}")
    else:
        print(f"  Some recovery signal persists through round 10")

    # Multi-word verdict
    print()
    any_multi_better = False
    for rnd in [7, 8, 9, 10]:
        if rnd in p2_single and rnd in p2_multi:
            ss, ts = p2_single[rnd]
            sm, tm = p2_multi[rnd]
            if ts > 0 and tm > 0 and sm / tm > ss / ts + 0.03:
                any_multi_better = True
    print(f"  Multi-word cubes help? {'YES' if any_multi_better else 'NO'}")

    # Conditional verdict
    cond_helps = False
    for rnd in [8, 9, 10]:
        none_s, none_t = p4_results.get("none", {}).get(rnd, (0, 0))
        if none_t == 0:
            continue
        none_rate = none_s / none_t
        for cn in ["zeros_hi", "ones_hi", "alt_hi"]:
            cs, ct = p4_results.get(cn, {}).get(rnd, (0, 0))
            if ct > 0 and cs / ct > none_rate + 0.05:
                cond_helps = True
    print(f"  Conditioning helps? {'YES' if cond_helps else 'NO'}")

    print(f"\n  Runtime: {elapsed:.1f}s")
    print("=" * 72)


if __name__ == "__main__":
    main()
