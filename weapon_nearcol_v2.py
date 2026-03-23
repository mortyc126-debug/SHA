#!/usr/bin/env python3
"""
SHA-256 Near-Collision Finding Machine v2
=========================================
Combines guided differential search, filtered message conditions, multi-path
exploration, and birthday-enhanced collision finding for reduced-round SHA-256.

Phase 1: Guided near-collision search (local collision path, 200k trials)
Phase 2: Filtered search (HW/XOR conditions on W[0], 100k trials)
Phase 3: Multi-path search (4 differential paths, 100k trials each)
Phase 4: Birthday-enhanced search (sorted partial-hash matching)

Target runtime: under 180 seconds.
"""

import time
import random
from collections import defaultdict

# ============================================================================
# SHA-256 constants (real, first 16 round constants)
# ============================================================================
K = [0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
     0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
     0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
     0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174]

IV = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]

M = 0xFFFFFFFF

# ============================================================================
# Primitives — inlined for speed where possible
# ============================================================================

def hw(x):
    return bin(x).count('1')

def sha256_rounds(W, nr):
    """Run nr rounds of SHA-256 compression from standard IV."""
    a, b, c, d, e, f, g, h = 0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, \
                               0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
    for i in range(nr):
        S1 = (((e >> 6) | (e << 26)) ^ ((e >> 11) | (e << 21)) ^ ((e >> 25) | (e << 7))) & M
        ch = (e & f) ^ ((~e & M) & g)
        t1 = (h + S1 + ch + K[i] + W[i]) & M
        S0 = (((a >> 2) | (a << 30)) ^ ((a >> 13) | (a << 19)) ^ ((a >> 22) | (a << 10))) & M
        maj = (a & b) ^ (a & c) ^ (b & c)
        t2 = (S0 + maj) & M
        h = g; g = f; f = e; e = (d + t1) & M
        d = c; c = b; b = a; a = (t1 + t2) & M
    return (a, b, c, d, e, f, g, h)

def state_hamming(s1, s2):
    """Total Hamming distance across all 8 registers."""
    t = 0
    for i in range(8):
        t += bin(s1[i] ^ s2[i]).count('1')
    return t

def per_register_hw(s1, s2):
    return [bin(s1[i] ^ s2[i]).count('1') for i in range(8)]

def rand_msg(rng, n=16):
    r = rng.randint
    return [r(0, M) for _ in range(n)]

def make_high_hw_word(rng, min_hw=21):
    """Generate a 32-bit word with Hamming weight > min_hw."""
    while True:
        w = rng.randint(0, M)
        if bin(w).count('1') > min_hw:
            return w

def make_close_word(rng, base, max_hw=8):
    """Generate a word that has HW(base ^ word) <= max_hw by flipping few bits."""
    nbits = rng.randint(1, max_hw)
    positions = rng.sample(range(32), nbits)
    flip = 0
    for p in positions:
        flip |= (1 << p)
    return base ^ flip

# ============================================================================
# Phase 1: Guided near-collision search with local collision path
# ============================================================================

def phase1():
    print("=" * 72)
    print("PHASE 1: Guided Near-Collision Search (Local Collision Path)")
    print("  DeltaW[0]=0x80000000, DeltaW[1]=0x821001c0, DeltaW[2..15]=0")
    print("  200,000 random base messages per round count")
    print("=" * 72)
    t0 = time.time()

    delta = [0x80000000, 0x821001c0] + [0] * 14
    N = 200000
    rng = random.Random(12345)

    for nr in [4, 5, 6, 7, 8]:
        best_hw_val = 256
        best_per_reg = None
        buckets = {10: 0, 20: 0, 30: 0}

        for _ in range(N):
            W = rand_msg(rng)
            W2 = [W[j] ^ delta[j] for j in range(16)]
            h1 = sha256_rounds(W, nr)
            h2 = sha256_rounds(W2, nr)
            d = state_hamming(h1, h2)

            if d <= 30:
                buckets[30] += 1
                if d <= 20:
                    buckets[20] += 1
                    if d <= 10:
                        buckets[10] += 1

            if d < best_hw_val:
                best_hw_val = d
                best_per_reg = per_register_hw(h1, h2)

        print(f"\n  {nr}-round SHA-256:")
        print(f"    Best near-collision HW = {best_hw_val}")
        print(f"    Per-register HW: {best_per_reg}")
        close_regs = [i for i, v in enumerate(best_per_reg) if v <= 3]
        if close_regs:
            reg_names = "abcdefgh"
            print(f"    Close registers (HW<=3): {[reg_names[i] for i in close_regs]}")
        print(f"    Distribution (of {N}):")
        for thr in [10, 20, 30]:
            pct = buckets[thr] / N * 100
            print(f"      HW <= {thr:2d}: {buckets[thr]:7d}  ({pct:.3f}%)")

    elapsed = time.time() - t0
    print(f"\n  Phase 1 time: {elapsed:.1f}s")
    return elapsed

# ============================================================================
# Phase 2: Filtered search — condition on W[0] properties
# ============================================================================

def phase2():
    print("\n" + "=" * 72)
    print("PHASE 2: Filtered Search (Message Conditions)")
    print("  Filter: HW(W[0]) > 20 AND HW(W[0] XOR W[1]) <= 8")
    print("  100,000 filtered + 100,000 unfiltered pairs per round count")
    print("=" * 72)
    t0 = time.time()

    delta = [0x80000000, 0x821001c0] + [0] * 14
    N = 100000

    for nr in [4, 5, 6, 7, 8]:
        # --- Unfiltered run ---
        rng_u = random.Random(54321 + nr)
        best_hw_unf = 256
        best_per_reg_unf = None
        buckets_unf = {10: 0, 20: 0, 30: 0}

        for _ in range(N):
            W = rand_msg(rng_u)
            W2 = [W[j] ^ delta[j] for j in range(16)]
            h1 = sha256_rounds(W, nr)
            h2 = sha256_rounds(W2, nr)
            d = state_hamming(h1, h2)
            if d <= 30:
                buckets_unf[30] += 1
                if d <= 20:
                    buckets_unf[20] += 1
                    if d <= 10:
                        buckets_unf[10] += 1
            if d < best_hw_unf:
                best_hw_unf = d
                best_per_reg_unf = per_register_hw(h1, h2)

        # --- Filtered run: directly construct satisfying messages ---
        rng_f = random.Random(98765 + nr)
        best_hw_flt = 256
        best_per_reg_flt = None
        buckets_flt = {10: 0, 20: 0, 30: 0}

        for _ in range(N):
            # Construct W[0] with high HW and W[1] close to W[0]
            w0 = make_high_hw_word(rng_f, min_hw=20)
            w1 = make_close_word(rng_f, w0, max_hw=8)
            W = [w0, w1] + [rng_f.randint(0, M) for _ in range(14)]
            W2 = [W[j] ^ delta[j] for j in range(16)]
            h1 = sha256_rounds(W, nr)
            h2 = sha256_rounds(W2, nr)
            d = state_hamming(h1, h2)
            if d <= 30:
                buckets_flt[30] += 1
                if d <= 20:
                    buckets_flt[20] += 1
                    if d <= 10:
                        buckets_flt[10] += 1
            if d < best_hw_flt:
                best_hw_flt = d
                best_per_reg_flt = per_register_hw(h1, h2)

        print(f"\n  {nr}-round SHA-256:")
        print(f"    Unfiltered best HW = {best_hw_unf}  per-reg: {best_per_reg_unf}")
        print(f"    Filtered   best HW = {best_hw_flt}  per-reg: {best_per_reg_flt}")
        improvement = best_hw_unf - best_hw_flt
        tag = "BETTER" if improvement > 0 else ("SAME" if improvement == 0 else "worse")
        print(f"    Improvement: {improvement:+d} bits [{tag}]")
        print(f"    Distribution comparison (HW <= threshold):")
        for thr in [10, 20, 30]:
            pu = buckets_unf[thr] / N * 100
            pf = buckets_flt[thr] / N * 100
            ratio = pf / pu if pu > 0 else float('inf')
            print(f"      HW<={thr:2d}: unfiltered={buckets_unf[thr]:6d} ({pu:.3f}%)  "
                  f"filtered={buckets_flt[thr]:6d} ({pf:.3f}%)  ratio={ratio:.2f}x")

    elapsed = time.time() - t0
    print(f"\n  Phase 2 time: {elapsed:.1f}s")
    return elapsed

# ============================================================================
# Phase 3: Multi-path search
# ============================================================================

def hill_climb_dw2(rng, nr, base_dw0, base_dw1, n_seeds=200, steps=50):
    """Hill-climb to find a good DW[2] given fixed DW[0], DW[1].
    Uses a small evaluation sample for speed."""
    sample_size = 50  # small for speed
    msgs = [rand_msg(rng) for _ in range(sample_size)]

    best_dw2_overall = 0
    best_score_overall = 999.0

    for seed_i in range(n_seeds):
        rng_local = random.Random(seed_i * 1000 + nr)
        dw2 = rng_local.randint(0, M)

        def evaluate(dw2_val):
            delta_l = [base_dw0, base_dw1, dw2_val] + [0] * 13
            total_d = 0
            for W in msgs:
                W2 = [W[j] ^ delta_l[j] for j in range(16)]
                h1 = sha256_rounds(W, nr)
                h2 = sha256_rounds(W2, nr)
                total_d += state_hamming(h1, h2)
            return total_d / sample_size

        score = evaluate(dw2)

        for _ in range(steps):
            bit = rng_local.randint(0, 31)
            candidate = dw2 ^ (1 << bit)
            new_score = evaluate(candidate)
            if new_score < score:
                score = new_score
                dw2 = candidate

        if score < best_score_overall:
            best_score_overall = score
            best_dw2_overall = dw2

    return best_dw2_overall, best_score_overall

def phase3():
    print("\n" + "=" * 72)
    print("PHASE 3: Multi-Path Search")
    print("  Path A: DW[0]=0x80000000, DW[1]=0x821001c0")
    print("  Path B: DW[0]=0x80000000, DW[1]=0xfe100040")
    print("  Path C: DW[0]=0x00000024, DW[1]=0 (2-bit, no correction)")
    print("  Path D: DW[0]=0x80000000, DW[1]=0, DW[2]=hill-climbed")
    print("  100,000 trials each at 6 and 8 rounds")
    print("=" * 72)
    t0 = time.time()

    N = 100000

    # Define paths A, B, C
    paths_base = {
        'A': [0x80000000, 0x821001c0] + [0] * 14,
        'B': [0x80000000, 0xfe100040] + [0] * 14,
        'C': [0x00000024, 0x00000000] + [0] * 14,
    }

    # Path D: hill-climb DW[2] for each round count
    print("\n  Hill-climbing Path D (DW[2]) — 200 seeds x 50 steps...")
    path_d_dw2 = {}
    for nr in [6, 8]:
        hc_rng = random.Random(99999 + nr)
        dw2, score = hill_climb_dw2(hc_rng, nr, 0x80000000, 0, n_seeds=200, steps=50)
        path_d_dw2[nr] = dw2
        print(f"    {nr}-round: best DW[2] = 0x{dw2:08x} (mean HW={score:.1f})")

    for nr in [6, 8]:
        paths = dict(paths_base)
        paths['D'] = [0x80000000, 0x00000000, path_d_dw2[nr]] + [0] * 13

        print(f"\n  --- {nr}-round results ({N} trials per path) ---")

        for pname in ['A', 'B', 'C', 'D']:
            delta = paths[pname]
            rng = random.Random(77777 + ord(pname[0]) + nr * 100)
            best_hw_val = 256
            best_per_reg = None
            buckets = {10: 0, 20: 0, 30: 0}
            hw_sum = 0

            for _ in range(N):
                W = rand_msg(rng)
                W2 = [W[j] ^ delta[j] for j in range(16)]
                h1 = sha256_rounds(W, nr)
                h2 = sha256_rounds(W2, nr)
                d = state_hamming(h1, h2)
                hw_sum += d
                if d <= 30:
                    buckets[30] += 1
                    if d <= 20:
                        buckets[20] += 1
                        if d <= 10:
                            buckets[10] += 1
                if d < best_hw_val:
                    best_hw_val = d
                    best_per_reg = per_register_hw(h1, h2)

            mean_hw = hw_sum / N
            delta_nonzero = [(j, delta[j]) for j in range(16) if delta[j] != 0]
            delta_str = ", ".join(f"DW[{j}]=0x{v:08x}" for j, v in delta_nonzero)
            if not delta_str:
                delta_str = "all zero"
            print(f"    Path {pname} [{delta_str}]:")
            print(f"      Best HW = {best_hw_val}, Mean HW = {mean_hw:.1f}")
            print(f"      Per-register: {best_per_reg}")
            print(f"      HW<=10: {buckets[10]}  HW<=20: {buckets[20]}  HW<=30: {buckets[30]}")

    elapsed = time.time() - t0
    print(f"\n  Phase 3 time: {elapsed:.1f}s")
    return elapsed

# ============================================================================
# Phase 4: Birthday-enhanced search
# ============================================================================

def phase4():
    print("\n" + "=" * 72)
    print("PHASE 4: Birthday-Enhanced Search")
    print("  Generate N hashes, sort by partial hash, find closest pairs")
    print("  Compare with random birthday bound")
    print("=" * 72)
    t0 = time.time()

    N = 50000
    delta = [0x80000000, 0x821001c0] + [0] * 14

    for nr in [6, 8]:
        print(f"\n  --- {nr}-round SHA-256 (N={N}) ---")

        rng = random.Random(314159 + nr)

        # Generate messages and compute hash pairs
        hashes_orig = []
        hashes_flip = []

        for _ in range(N):
            W = rand_msg(rng)
            h1 = sha256_rounds(W, nr)
            W2 = [W[j] ^ delta[j] for j in range(16)]
            h2 = sha256_rounds(W2, nr)
            hashes_orig.append(h1)
            hashes_flip.append(h2)

        # --- Method A: Standard differential (best pair from direct diff) ---
        best_hw_std = 256
        best_idx_std = 0
        for i in range(N):
            d = state_hamming(hashes_orig[i], hashes_flip[i])
            if d < best_hw_std:
                best_hw_std = d
                best_idx_std = i

        print(f"\n    [A] Standard differential search ({N} pairs):")
        print(f"      Best near-collision HW = {best_hw_std}")
        print(f"      Per-register: {per_register_hw(hashes_orig[best_idx_std], hashes_flip[best_idx_std])}")

        # --- Method B: Birthday on register 'e' (index 4) upper 16 bits ---
        print(f"\n    [B] Birthday on upper 16 bits of register 'e':")

        partial_bits = 16
        partial_shift = 32 - partial_bits

        # Orig-orig birthday
        keyed = sorted((hashes_orig[i][4] >> partial_shift, i) for i in range(N))
        matches_oo = 0
        best_hw_oo = 256
        for i in range(len(keyed) - 1):
            if keyed[i][0] == keyed[i + 1][0]:
                matches_oo += 1
                ia, ib = keyed[i][1], keyed[i + 1][1]
                d = state_hamming(hashes_orig[ia], hashes_orig[ib])
                if d < best_hw_oo:
                    best_hw_oo = d

        expected_birthday = N * (N - 1) / (2 * (1 << partial_bits))
        print(f"      Orig-orig partial matches: {matches_oo} (expected random: {expected_birthday:.0f})")
        ratio_oo = matches_oo / expected_birthday if expected_birthday > 0 else 0
        print(f"      Ratio: {ratio_oo:.2f}x")
        if best_hw_oo < 256:
            print(f"      Best orig-orig near-collision HW = {best_hw_oo}")

        # Flip-flip birthday
        keyed_f = sorted((hashes_flip[i][4] >> partial_shift, i) for i in range(N))
        matches_ff = 0
        best_hw_ff = 256
        for i in range(len(keyed_f) - 1):
            if keyed_f[i][0] == keyed_f[i + 1][0]:
                matches_ff += 1
                ia, ib = keyed_f[i][1], keyed_f[i + 1][1]
                d = state_hamming(hashes_flip[ia], hashes_flip[ib])
                if d < best_hw_ff:
                    best_hw_ff = d

        print(f"      Flip-flip partial matches: {matches_ff} (expected: {expected_birthday:.0f})")
        if best_hw_ff < 256:
            print(f"      Best flip-flip near-collision HW = {best_hw_ff}")

        # --- Method C: Cross-birthday (orig[i] vs flip[j]) ---
        print(f"\n    [C] Cross-birthday (orig[i] vs flip[j]) on upper 16 bits of 'e':")

        flip_by_key = defaultdict(list)
        for i in range(N):
            pk = hashes_flip[i][4] >> partial_shift
            flip_by_key[pk].append(i)

        cross_matches = 0
        best_cross_hw = 256
        best_cross_pair = None
        for i in range(N):
            pk = hashes_orig[i][4] >> partial_shift
            if pk in flip_by_key:
                for j in flip_by_key[pk]:
                    if i == j:
                        continue
                    cross_matches += 1
                    d = state_hamming(hashes_orig[i], hashes_flip[j])
                    if d < best_cross_hw:
                        best_cross_hw = d
                        best_cross_pair = (i, j)

        expected_cross = N * (N - 1) / (1 << partial_bits)
        print(f"      Cross-matches found: {cross_matches}")
        print(f"      Expected random: {expected_cross:.0f}")
        ratio_cross = cross_matches / expected_cross if expected_cross > 0 else 0
        print(f"      Ratio (observed/expected): {ratio_cross:.2f}x")
        if best_cross_pair:
            i, j = best_cross_pair
            print(f"      Best cross near-collision HW = {best_cross_hw}")
            print(f"      Per-register: {per_register_hw(hashes_orig[i], hashes_flip[j])}")

        # --- Method D: Full birthday — closest pair among all 2N hashes ---
        print(f"\n    [D] Full birthday (closest pair among all 2N={2*N} hashes):")

        # Use upper 8 bits of registers a and e as 16-bit sort key
        all_entries = []
        for i in range(N):
            key_o = ((hashes_orig[i][0] >> 24) << 8) | (hashes_orig[i][4] >> 24)
            all_entries.append((key_o, 'O', i))
            key_f = ((hashes_flip[i][0] >> 24) << 8) | (hashes_flip[i][4] >> 24)
            all_entries.append((key_f, 'F', i))

        all_entries.sort()

        matches_16bit = 0
        best_bday_hw = 256
        best_bday_desc = ""
        best_bday_per_reg = None
        for i in range(len(all_entries) - 1):
            if all_entries[i][0] == all_entries[i + 1][0]:
                matches_16bit += 1
                t1, idx1 = all_entries[i][1], all_entries[i][2]
                t2, idx2 = all_entries[i + 1][1], all_entries[i + 1][2]
                h1 = hashes_orig[idx1] if t1 == 'O' else hashes_flip[idx1]
                h2 = hashes_orig[idx2] if t2 == 'O' else hashes_flip[idx2]
                d = state_hamming(h1, h2)
                if d < best_bday_hw:
                    best_bday_hw = d
                    best_bday_desc = f"{t1}[{idx1}] vs {t2}[{idx2}]"
                    best_bday_per_reg = per_register_hw(h1, h2)

        total_2n = 2 * N
        expected_16 = total_2n * (total_2n - 1) / (2 * (1 << 16))
        print(f"      16-bit partial matches: {matches_16bit} (expected: {expected_16:.0f})")
        if best_bday_hw < 256:
            print(f"      Best near-collision HW = {best_bday_hw} ({best_bday_desc})")
            print(f"      Per-register: {best_bday_per_reg}")
        print(f"      Improvement over standard diff: {best_hw_std - best_bday_hw:+d} bits")

    elapsed = time.time() - t0
    print(f"\n  Phase 4 time: {elapsed:.1f}s")
    return elapsed

# ============================================================================
# Main
# ============================================================================

def main():
    print("+" + "=" * 70 + "+")
    print("|    SHA-256 Near-Collision Finding Machine v2                       |")
    print("|    Combining differential, filtered, multi-path, birthday          |")
    print("+" + "=" * 70 + "+")
    print()

    total_start = time.time()

    t1 = phase1()
    t2 = phase2()
    t3 = phase3()
    t4 = phase4()

    total = time.time() - total_start

    print("\n" + "=" * 72)
    print("TIMING SUMMARY")
    print("=" * 72)
    print(f"  Phase 1 (Guided search):    {t1:6.1f}s")
    print(f"  Phase 2 (Filtered search):  {t2:6.1f}s")
    print(f"  Phase 3 (Multi-path):       {t3:6.1f}s")
    print(f"  Phase 4 (Birthday):         {t4:6.1f}s")
    print(f"  TOTAL:                      {total:6.1f}s")
    if total < 180:
        print(f"  Within 180s budget with {180 - total:.1f}s to spare.")
    else:
        print(f"  WARNING: Exceeded 180s budget by {total - 180:.1f}s.")

if __name__ == "__main__":
    main()
