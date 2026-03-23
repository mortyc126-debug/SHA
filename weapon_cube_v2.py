#!/usr/bin/env python3
"""
weapon_cube_v2.py — Advanced Cube Attack on Reduced SHA-256

Pushes cube attack beyond the 6-round key-recovery boundary.
Phase 1: Higher-dimension cubes (dim 6-12) for rounds 7-10
Phase 2: Multi-word cubes (bits spread across W[0] AND W[1])
Phase 3: Optimal bit selection using Sigma1-related positions
Phase 4: Conditional cube attack with favorable bit fixing

Methodology:
  The cube attack recovers secret bit s if the cube sum (XOR over all 2^dim
  assignments to free positions) depends on s. For a single output bit:
    cube_sum(s=0) XOR cube_sum(s=1) = 1  if superpoly captures s
    cube_sum(s=0) XOR cube_sum(s=1) = 0  if superpoly is independent of s

  For a FIXED cube subset and FIXED remaining message words ("key"),
  this is deterministic. Across random keys, if the superpoly truly
  depends on s, the difference is 1 for 100% of keys. If independent,
  the difference is random with ~50% chance of being 1.

  We test N_KEYS random key settings per cube subset, and report the
  fraction where cube_sum differs. True recovery = 100%. Null = ~50%.
  We declare a cube subset "works" if ALL N_KEYS trials show difference.
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
    to free_positions. W_template already has secret/key set.
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


def test_cube_single_key(W_key, free_positions, secret_word, secret_bit,
                         n_rounds, out_bit):
    """
    For a single key setting, test if cube sum differs between secret=0/1.
    Returns 1 if they differ, 0 if same.
    """
    W0 = list(W_key)
    W0[secret_word] &= ~(1 << secret_bit)
    s0 = cube_sum(W0, free_positions, n_rounds, out_bit)

    W1 = list(W_key)
    W1[secret_word] |= (1 << secret_bit)
    s1 = cube_sum(W1, free_positions, n_rounds, out_bit)

    return 1 if s0 != s1 else 0


def test_cube_multi_key(free_positions, secret_word, secret_bit,
                        n_rounds, out_bit, n_keys=10):
    """
    Test cube across n_keys random key settings.
    Returns fraction of keys where cube sum differs (s0 != s1).
    True cube attack: fraction = 1.0
    Null (random): fraction ~ 0.5
    """
    n_differ = 0
    for _ in range(n_keys):
        W = [random.getrandbits(32) for _ in range(16)]
        n_differ += test_cube_single_key(W, free_positions, secret_word,
                                         secret_bit, n_rounds, out_bit)
    return n_differ / n_keys


# ======================================================================
#  PHASE 0: Baseline verification
# ======================================================================
def phase0(t0, time_budget):
    print("=" * 72)
    print("PHASE 0: Baseline verification")
    print("=" * 72)
    print("  Testing rounds 3-7, dim 4-6, N_KEYS=20 per cube subset")
    print("  Null hypothesis: 50% differ rate")
    print("  True cube attack: 100% differ rate")
    print()

    secret_word, secret_bit = 1, 0
    out_bit = 0
    n_keys = 20
    n_subsets = 15

    for rnd in [3, 4, 5, 6, 7]:
        for dim in [4, 6]:
            if time.time() - t0 > time_budget * 0.08:
                break
            rates = []
            n_perfect = 0
            for _ in range(n_subsets):
                bits = random.sample(range(32), dim)
                free_pos = [(0, b) for b in bits]
                rate = test_cube_multi_key(free_pos, secret_word, secret_bit,
                                           rnd, out_bit, n_keys)
                rates.append(rate)
                if rate == 1.0:
                    n_perfect += 1

            avg_rate = sum(rates) / len(rates)
            print(f"  round={rnd} dim={dim}: avg differ rate = {avg_rate:.3f}  "
                  f"perfect(100%) = {n_perfect}/{n_subsets}  "
                  f"{'<-- WORKING' if avg_rate > 0.8 else '(~random)' if avg_rate < 0.6 else '(weak signal)'}")

    print()


# ======================================================================
#  PHASE 1: Higher-dimension cubes for rounds 7-10
# ======================================================================
def phase1(t0, time_budget):
    print("=" * 72)
    print("PHASE 1: Higher-dimension cubes — rounds 7-10")
    print("=" * 72)
    print("Free bits in W[0], secret = W[1] bit 0, output = a[0]")
    print("N_KEYS=20 per cube subset, N_SUBSETS varies by dim")
    print()

    dims = [6, 8, 10, 12]
    rounds_list = [7, 8, 9, 10]
    secret_word, secret_bit = 1, 0
    out_bit = 0
    n_keys = 20

    results = defaultdict(lambda: defaultdict(lambda: (0.0, 0, 0)))

    for dim in dims:
        n_subsets = 20 if dim >= 12 else (25 if dim >= 10 else 50)
        for rnd in rounds_list:
            elapsed = time.time() - t0
            if elapsed > time_budget * 0.35:
                n_subsets = min(n_subsets, 8)
            if elapsed > time_budget * 0.42:
                break

            rates = []
            n_perfect = 0
            tested = 0
            for _ in range(n_subsets):
                if time.time() - t0 > time_budget * 0.44:
                    break
                bits = random.sample(range(32), dim)
                free_pos = [(0, b) for b in bits]
                rate = test_cube_multi_key(free_pos, secret_word, secret_bit,
                                           rnd, out_bit, n_keys)
                rates.append(rate)
                if rate == 1.0:
                    n_perfect += 1
                tested += 1

            if tested > 0:
                avg = sum(rates) / len(rates)
                results[rnd][dim] = (avg, n_perfect, tested)
                signal = "WORKING" if avg > 0.8 else "weak" if avg > 0.6 else "~null"
                print(f"  rounds={rnd:2d}  dim={dim:2d}  avg_differ={avg:.3f}  "
                      f"perfect={n_perfect:3d}/{tested:3d}  [{signal}]")

    # Summary table
    print()
    print("  Average differ rate (null=0.50, perfect=1.00):")
    print(f"  {'':8s}", end="")
    for dim in dims:
        print(f"  dim={dim:2d} ", end="")
    print()
    for rnd in rounds_list:
        print(f"  rnd={rnd:2d}  ", end="")
        for dim in dims:
            avg, pf, t = results[rnd][dim]
            if t > 0:
                print(f"   {avg:.3f}", end="")
            else:
                print(f"     -- ", end="")
        print()

    # Determine boundary
    print()
    for dim in dims:
        best_rnd = None
        for rnd in rounds_list:
            avg, pf, t = results[rnd][dim]
            if t > 0 and avg > 0.7:
                best_rnd = rnd
        if best_rnd:
            avg, pf, t = results[best_rnd][dim]
            print(f"  dim={dim:2d}: signal through round {best_rnd} "
                  f"(avg_differ={avg:.3f})")
        else:
            # Check if any round has even weak signal
            for rnd in rounds_list:
                avg, pf, t = results[rnd][dim]
                if t > 0 and avg > 0.55:
                    print(f"  dim={dim:2d}: WEAK signal at round {rnd} "
                          f"(avg_differ={avg:.3f})")
                    break
            else:
                print(f"  dim={dim:2d}: no signal above null at rounds 7+")

    return results


# ======================================================================
#  PHASE 2: Multi-word cube
# ======================================================================
def phase2(t0, time_budget):
    print()
    print("=" * 72)
    print("PHASE 2: Multi-word cube — 4 bits W[0] + 4 bits W[1]")
    print("=" * 72)
    print("Dim 8 total, secret = W[1] bit 0 (excluded from free bits)")
    print()

    rounds_list = [7, 8, 9, 10]
    secret_word, secret_bit = 1, 0
    out_bit = 0
    n_keys = 20
    n_subsets = 20

    results_single = {}
    results_multi = {}

    for rnd in rounds_list:
        if time.time() - t0 > time_budget * 0.62:
            break

        # Single-word: 8 bits in W[0]
        rates_s = []
        for _ in range(n_subsets):
            if time.time() - t0 > time_budget * 0.56:
                break
            bits = random.sample(range(32), 8)
            free_pos = [(0, b) for b in bits]
            rate = test_cube_multi_key(free_pos, secret_word, secret_bit,
                                       rnd, out_bit, n_keys)
            rates_s.append(rate)
        results_single[rnd] = rates_s

        # Multi-word: 4+4
        rates_m = []
        w1_cands = [b for b in range(32) if b != secret_bit]
        for _ in range(n_subsets):
            if time.time() - t0 > time_budget * 0.60:
                break
            bits_w0 = random.sample(range(32), 4)
            bits_w1 = random.sample(w1_cands, 4)
            free_pos = [(0, b) for b in bits_w0] + [(1, b) for b in bits_w1]
            rate = test_cube_multi_key(free_pos, secret_word, secret_bit,
                                       rnd, out_bit, n_keys)
            rates_m.append(rate)
        results_multi[rnd] = rates_m

        avg_s = sum(rates_s) / len(rates_s) if rates_s else 0
        avg_m = sum(rates_m) / len(rates_m) if rates_m else 0
        print(f"  rounds={rnd:2d}  single-word avg={avg_s:.3f}  "
              f"multi-word avg={avg_m:.3f}  "
              f"{'MULTI BETTER' if avg_m > avg_s + 0.05 else ''}")

    print()
    improvement = False
    for rnd in rounds_list:
        rs = results_single.get(rnd, [])
        rm = results_multi.get(rnd, [])
        if rs and rm:
            if sum(rm) / len(rm) > sum(rs) / len(rs) + 0.05:
                improvement = True
    print(f"  Verdict: Multi-word spreading "
          f"{'HELPS' if improvement else 'shows NO significant improvement'}")

    return results_single, results_multi


# ======================================================================
#  PHASE 3: Optimal bit selection
# ======================================================================
def phase3(t0, time_budget):
    print()
    print("=" * 72)
    print("PHASE 3: Sigma1-aligned cube vs random")
    print("=" * 72)
    print("Optimal dim-4: free bits {5, 10, 24, 31}")
    print("Optimal dim-8: free bits {5, 6, 10, 11, 20, 24, 25, 31}")
    print()

    rounds_list = [5, 6, 7, 8, 9, 10]
    secret_word, secret_bit = 1, 0
    out_bit = 0
    n_keys = 20
    n_subsets = 20

    configs = [
        ("Sig1-d4", [(0, b) for b in [5, 10, 24, 31]], 4),
        ("Sig1-d8", [(0, b) for b in [5, 6, 10, 11, 20, 24, 25, 31]], 8),
    ]

    for label, optimal_pos, dim in configs:
        print(f"  [{label}]")
        for rnd in rounds_list:
            if time.time() - t0 > time_budget * 0.80:
                break

            # Optimal (same subset each time, different keys)
            rate_opt = test_cube_multi_key(optimal_pos, secret_word, secret_bit,
                                           rnd, out_bit, n_keys * 2)

            # Random subsets
            rates_rnd = []
            for _ in range(n_subsets):
                if time.time() - t0 > time_budget * 0.78:
                    break
                bits = random.sample(range(32), dim)
                free_pos = [(0, b) for b in bits]
                r = test_cube_multi_key(free_pos, secret_word, secret_bit,
                                        rnd, out_bit, n_keys)
                rates_rnd.append(r)

            avg_rnd = sum(rates_rnd) / len(rates_rnd) if rates_rnd else 0
            marker = ""
            if rate_opt > avg_rnd + 0.1:
                marker = " <-- OPTIMAL BETTER"
            elif avg_rnd > rate_opt + 0.1:
                marker = " <-- RANDOM BETTER"
            print(f"    round={rnd:2d}  optimal={rate_opt:.3f}  "
                  f"random_avg={avg_rnd:.3f}{marker}")
        print()


# ======================================================================
#  PHASE 4: Conditional cube attack
# ======================================================================
def phase4(t0, time_budget):
    print()
    print("=" * 72)
    print("PHASE 4: Conditional cube attack")
    print("=" * 72)
    print("Fix upper byte of W[0], cube over bits 0-7 (dim 8)")
    print()

    rounds_list = [8, 9, 10]
    secret_word, secret_bit = 1, 0
    out_bit = 0
    n_keys = 20
    n_subsets = 15

    conditions = [
        ("zeros_hi", 0x00000000, 0xFF000000),
        ("ones_hi",  0xFF000000, 0xFF000000),
        ("alt_hi",   0xAA000000, 0xFF000000),
        ("none",     0x00000000, 0x00000000),
    ]

    free_pos = [(0, b) for b in range(8)]

    results = {}

    for cond_name, cond_val, cond_mask in conditions:
        results[cond_name] = {}
        for rnd in rounds_list:
            if time.time() - t0 > time_budget * 0.95:
                break

            rates = []
            for _ in range(n_subsets):
                if time.time() - t0 > time_budget * 0.93:
                    break
                # Test with conditioned keys
                n_differ = 0
                for __ in range(n_keys):
                    W = [random.getrandbits(32) for _ in range(16)]
                    W[0] = (W[0] & ~cond_mask) | (cond_val & cond_mask)
                    n_differ += test_cube_single_key(W, free_pos, secret_word,
                                                     secret_bit, rnd, out_bit)
                rates.append(n_differ / n_keys)

            if rates:
                avg = sum(rates) / len(rates)
                results[cond_name][rnd] = avg
            else:
                results[cond_name][rnd] = 0.0

    # Print
    print(f"  {'Condition':<12s}", end="")
    for rnd in rounds_list:
        print(f"  round={rnd}", end="")
    print()
    print("  " + "-" * 50)

    for cond_name, _, _ in conditions:
        print(f"  {cond_name:<12s}", end="")
        for rnd in rounds_list:
            val = results.get(cond_name, {}).get(rnd, None)
            if val is not None:
                print(f"   {val:.3f} ", end="")
            else:
                print(f"     --  ", end="")
        print()

    print()
    cond_helps = False
    for rnd in rounds_list:
        none_rate = results.get("none", {}).get(rnd, 0.5)
        for cn in ["zeros_hi", "ones_hi", "alt_hi"]:
            cr = results.get(cn, {}).get(rnd, 0.5)
            if cr > none_rate + 0.08:
                cond_helps = True
    print(f"  Verdict: Conditioning "
          f"{'IMPROVES cube attack' if cond_helps else 'shows NO significant improvement'}")

    return results


# ======================================================================
#  MAIN
# ======================================================================
def main():
    random.seed(0xC0BEA772)
    t0 = time.time()
    TIME_BUDGET = 115

    print("*" * 72)
    print("  WEAPON: Advanced Cube Attack on Reduced SHA-256 (v2)")
    print("  Baseline: 100% key-recovery through 6 rounds at dim 4-6")
    print("  Goal: Push cube attack boundary further")
    print("*" * 72)
    print()
    print("  Methodology:")
    print("    For each cube subset, test N_KEYS random key settings.")
    print("    Compute cube_sum(secret=0) XOR cube_sum(secret=1).")
    print("    True key recovery: differ rate = 1.000 (all keys)")
    print("    Null (no signal):  differ rate ~ 0.500 (random)")
    print()

    phase0(t0, TIME_BUDGET)
    p1 = phase1(t0, TIME_BUDGET)
    p2s, p2m = phase2(t0, TIME_BUDGET)
    phase3(t0, TIME_BUDGET)
    p4 = phase4(t0, TIME_BUDGET)

    # ── Final Summary ─────────────────────────────────────────────
    elapsed = time.time() - t0
    print()
    print("=" * 72)
    print("FINAL SUMMARY")
    print("=" * 72)

    # Find max round with avg_differ > 0.7
    max_round = 6
    for rnd in [7, 8, 9, 10]:
        for dim in [6, 8, 10, 12]:
            avg, pf, t = p1.get(rnd, {}).get(dim, (0.0, 0, 0))
            if t > 0 and avg > 0.7:
                max_round = max(max_round, rnd)
        for d in [p2s, p2m]:
            rs = d.get(rnd, [])
            if rs and sum(rs) / len(rs) > 0.7:
                max_round = max(max_round, rnd)

    print(f"\n  Key-recovery boundary (avg differ > 0.70): round {max_round}")

    # Find max round with any signal > 0.55
    max_signal = 6
    for rnd in [7, 8, 9, 10]:
        for dim in [6, 8, 10, 12]:
            avg, pf, t = p1.get(rnd, {}).get(dim, (0.0, 0, 0))
            if t > 0 and avg > 0.55:
                max_signal = max(max_signal, rnd)

    print(f"  Weak signal boundary (avg differ > 0.55): round {max_signal}")

    # Breakdown
    first_dead = None
    for rnd in [7, 8, 9, 10]:
        all_null = True
        for dim in [6, 8, 10, 12]:
            avg, pf, t = p1.get(rnd, {}).get(dim, (0.0, 0, 0))
            if t > 0 and avg > 0.55:
                all_null = False
        if all_null and first_dead is None:
            first_dead = rnd
    if first_dead:
        print(f"  Cube attack fully dead at round {first_dead}")
    else:
        print(f"  Some signal persists through round 10")

    # Multi-word
    multi_helps = False
    for rnd in [7, 8, 9, 10]:
        rs = p2s.get(rnd, [])
        rm = p2m.get(rnd, [])
        if rs and rm:
            as_ = sum(rs) / len(rs)
            am = sum(rm) / len(rm)
            if am > as_ + 0.05:
                multi_helps = True
    print(f"\n  Multi-word cubes help? {'YES' if multi_helps else 'NO'}")

    # Conditioning
    cond_helps = False
    for rnd in [8, 9, 10]:
        none_r = p4.get("none", {}).get(rnd, 0.5)
        for cn in ["zeros_hi", "ones_hi", "alt_hi"]:
            cr = p4.get(cn, {}).get(rnd, 0.5)
            if cr > none_r + 0.05:
                cond_helps = True
    print(f"  Conditioning helps? {'YES' if cond_helps else 'NO'}")

    print(f"\n  Runtime: {elapsed:.1f}s")
    print("=" * 72)


if __name__ == "__main__":
    main()
