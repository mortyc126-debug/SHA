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
import struct
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
# Primitives
# ============================================================================

def rotr(x, n):
    return ((x >> n) | (x << (32 - n))) & M

def hw(x):
    return bin(x).count('1')

def sha256_rounds(W, nr):
    """Run nr rounds of SHA-256 compression from standard IV.
    W must have at least nr words. Returns 8-word state tuple."""
    a, b, c, d, e, f, g, h = IV
    for i in range(nr):
        S1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25)
        ch = (e & f) ^ ((~e & M) & g)
        t1 = (h + S1 + ch + K[i] + W[i]) & M
        S0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22)
        maj = (a & b) ^ (a & c) ^ (b & c)
        t2 = (S0 + maj) & M
        h = g; g = f; f = e; e = (d + t1) & M
        d = c; c = b; b = a; a = (t1 + t2) & M
    return (a, b, c, d, e, f, g, h)

def state_hamming(s1, s2):
    """Total Hamming distance across all 8 registers."""
    total = 0
    for x, y in zip(s1, s2):
        total += hw(x ^ y)
    return total

def per_register_hw(s1, s2):
    """Return list of per-register Hamming weights."""
    return [hw(x ^ y) for x, y in zip(s1, s2)]

def rand_word(rng):
    return rng.randint(0, M)

def rand_msg(rng, n=16):
    return [rng.randint(0, M) for _ in range(n)]

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
        best_msg = None
        best_per_reg = None
        buckets = {10: 0, 20: 0, 30: 0}

        for _ in range(N):
            W = rand_msg(rng)
            W2 = [W[j] ^ delta[j] for j in range(16)]
            h1 = sha256_rounds(W, nr)
            h2 = sha256_rounds(W2, nr)
            d = state_hamming(h1, h2)

            for thr in [10, 20, 30]:
                if d <= thr:
                    buckets[thr] += 1

            if d < best_hw_val:
                best_hw_val = d
                best_msg = W[:]
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
            expected_random = N * sum(1 for _ in [0]) * 0  # placeholder
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
    print("  100,000 filtered pairs per round count")
    print("=" * 72)
    t0 = time.time()

    delta = [0x80000000, 0x821001c0] + [0] * 14
    N = 100000
    rng = random.Random(54321)

    for nr in [4, 5, 6, 7, 8]:
        best_hw_unf = 256
        best_hw_flt = 256
        best_per_reg_unf = None
        best_per_reg_flt = None
        count_unf = 0
        count_flt = 0
        total_generated = 0
        buckets_unf = {10: 0, 20: 0, 30: 0}
        buckets_flt = {10: 0, 20: 0, 30: 0}

        while count_flt < N or count_unf < N:
            W = rand_msg(rng)
            total_generated += 1
            W2 = [W[j] ^ delta[j] for j in range(16)]
            h1 = sha256_rounds(W, nr)
            h2 = sha256_rounds(W2, nr)
            d = state_hamming(h1, h2)

            # Unfiltered tracking
            if count_unf < N:
                count_unf += 1
                for thr in [10, 20, 30]:
                    if d <= thr:
                        buckets_unf[thr] += 1
                if d < best_hw_unf:
                    best_hw_unf = d
                    best_per_reg_unf = per_register_hw(h1, h2)

            # Filtered tracking
            hw_w0 = hw(W[0])
            hw_xor = hw(W[0] ^ W[1])
            if hw_w0 > 20 and hw_xor <= 8:
                if count_flt < N:
                    count_flt += 1
                    for thr in [10, 20, 30]:
                        if d <= thr:
                            buckets_flt[thr] += 1
                    if d < best_hw_flt:
                        best_hw_flt = d
                        best_per_reg_flt = per_register_hw(h1, h2)

        print(f"\n  {nr}-round SHA-256:")
        print(f"    Generated {total_generated} messages to get {N} filtered")
        filt_rate = N / total_generated * 100
        print(f"    Filter acceptance rate: {filt_rate:.1f}%")
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

def hill_climb_dw2(rng, nr, base_dw0, base_dw1, steps=50):
    """Hill-climb to find a good DW[2] given fixed DW[0], DW[1].
    Minimize mean Hamming distance over a small sample."""
    sample_size = 200
    msgs = [rand_msg(rng) for _ in range(sample_size)]

    best_dw2 = rng.randint(0, M)
    delta = [base_dw0, base_dw1, best_dw2] + [0] * 13

    def evaluate(dw2_val):
        delta[2] = dw2_val
        total_d = 0
        for W in msgs:
            W2 = [W[j] ^ delta[j] for j in range(16)]
            h1 = sha256_rounds(W, nr)
            h2 = sha256_rounds(W2, nr)
            total_d += state_hamming(h1, h2)
        return total_d / sample_size

    best_score = evaluate(best_dw2)

    for _ in range(steps):
        # Flip a random bit
        bit = rng.randint(0, 31)
        candidate = best_dw2 ^ (1 << bit)
        score = evaluate(candidate)
        if score < best_score:
            best_score = score
            best_dw2 = candidate

    return best_dw2, best_score

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
    paths = {
        'A': [0x80000000, 0x821001c0] + [0] * 14,
        'B': [0x80000000, 0xfe100040] + [0] * 14,
        'C': [0x00000024, 0x00000000] + [0] * 14,
    }

    # Path D: hill-climb DW[2] for each round count
    print("\n  Hill-climbing Path D (DW[2])...")
    path_d_dw2 = {}
    for nr in [6, 8]:
        hc_rng = random.Random(99999 + nr)
        best_dw2_overall = 0
        best_score_overall = 999.0
        # 200 seeds x 50 steps
        seeds_per_batch = 200
        steps_per_seed = 50
        for seed_i in range(seeds_per_batch):
            hc_rng_local = random.Random(seed_i * 1000 + nr)
            dw2, score = hill_climb_dw2(hc_rng_local, nr, 0x80000000, 0, steps=steps_per_seed)
            if score < best_score_overall:
                best_score_overall = score
                best_dw2_overall = dw2
        path_d_dw2[nr] = best_dw2_overall
        print(f"    {nr}-round: best DW[2] = 0x{best_dw2_overall:08x} (mean HW={best_score_overall:.1f})")

    for nr in [6, 8]:
        paths_with_d = dict(paths)
        paths_with_d['D'] = [0x80000000, 0x00000000, path_d_dw2[nr]] + [0] * 13

        print(f"\n  --- {nr}-round results ({N} trials per path) ---")

        for pname in ['A', 'B', 'C', 'D']:
            delta = paths_with_d[pname]
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
                for thr in [10, 20, 30]:
                    if d <= thr:
                        buckets[thr] += 1
                if d < best_hw_val:
                    best_hw_val = d
                    best_per_reg = per_register_hw(h1, h2)

            mean_hw = hw_sum / N
            delta_hex = [f"0x{d:08x}" for d in delta[:3] if d != 0]
            if not delta_hex:
                delta_hex = ["0x00000000"]
            print(f"    Path {pname} [{','.join(delta_hex)}]:")
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

    delta_w0 = 0x80000000  # bit 31 in W[0]
    N = 50000

    for nr in [6, 8]:
        print(f"\n  --- {nr}-round SHA-256 (N={N}) ---")

        rng = random.Random(314159 + nr)

        # Generate messages and compute hashes for msg and msg^delta
        msgs = []
        hashes_orig = []
        hashes_flip = []

        delta = [0x80000000, 0x821001c0] + [0] * 14

        for _ in range(N):
            W = rand_msg(rng)
            msgs.append(W)
            h1 = sha256_rounds(W, nr)
            W2 = [W[j] ^ delta[j] for j in range(16)]
            h2 = sha256_rounds(W2, nr)
            hashes_orig.append(h1)
            hashes_flip.append(h2)

        # --- Method A: Standard differential (best pair) ---
        best_hw_std = 256
        best_idx_std = -1
        for i in range(N):
            d = state_hamming(hashes_orig[i], hashes_flip[i])
            if d < best_hw_std:
                best_hw_std = d
                best_idx_std = i

        print(f"    Standard differential search (N pairs):")
        print(f"      Best near-collision HW = {best_hw_std}")
        print(f"      Per-register: {per_register_hw(hashes_orig[best_idx_std], hashes_flip[best_idx_std])}")

        # --- Method B: Birthday on register 'e' (index 4) upper 16 bits ---
        # For each message, extract upper 16 bits of register e from original hash
        # Sort and find collisions on those 16 bits
        # Then among those collisions, check full Hamming distance between
        # their flipped hashes

        print(f"\n    Birthday search on upper 16 bits of register 'e':")

        # Build (partial_key, index) and sort
        partial_bits = 16
        partial_mask = ((1 << partial_bits) - 1) << (32 - partial_bits)

        keyed_orig = [(hashes_orig[i][4] & partial_mask, i) for i in range(N)]
        keyed_orig.sort()

        # Find matching partial keys in original hashes
        birthday_matches_orig = 0
        birthday_pairs_orig = []
        for i in range(len(keyed_orig) - 1):
            if keyed_orig[i][0] == keyed_orig[i + 1][0]:
                birthday_matches_orig += 1
                idx_a = keyed_orig[i][1]
                idx_b = keyed_orig[i + 1][1]
                birthday_pairs_orig.append((idx_a, idx_b))

        expected_random = N * (N - 1) / (2 * (1 << partial_bits))
        print(f"      Partial matches in original hashes: {birthday_matches_orig}")
        print(f"      Expected random birthday: {expected_random:.1f}")
        ratio = birthday_matches_orig / expected_random if expected_random > 0 else 0
        print(f"      Ratio (observed/expected): {ratio:.2f}x")

        # --- Method C: Birthday on flipped hashes ---
        keyed_flip = [(hashes_flip[i][4] & partial_mask, i) for i in range(N)]
        keyed_flip.sort()

        birthday_matches_flip = 0
        for i in range(len(keyed_flip) - 1):
            if keyed_flip[i][0] == keyed_flip[i + 1][0]:
                birthday_matches_flip += 1

        print(f"      Partial matches in flipped hashes: {birthday_matches_flip}")

        # --- Method D: Cross-birthday: match orig[i] with flip[j] ---
        # This finds near-collisions between H(m_i) and H(m_j ^ delta)
        # Meaning: two DIFFERENT messages whose hashes nearly collide
        print(f"\n    Cross-birthday (orig[i] vs flip[j]) on upper 16 bits of 'e':")

        # Build dict of partial keys for flipped
        flip_dict = defaultdict(list)
        for i in range(N):
            pk = hashes_flip[i][4] & partial_mask
            flip_dict[pk].append(i)

        cross_matches = 0
        best_cross_hw = 256
        best_cross_pair = None
        for i in range(N):
            pk = hashes_orig[i][4] & partial_mask
            if pk in flip_dict:
                for j in flip_dict[pk]:
                    if i == j:
                        continue
                    cross_matches += 1
                    d = state_hamming(hashes_orig[i], hashes_flip[j])
                    if d < best_cross_hw:
                        best_cross_hw = d
                        best_cross_pair = (i, j)

        expected_cross = N * (N - 1) / (1 << partial_bits)
        print(f"      Cross-matches found: {cross_matches}")
        print(f"      Expected random: {expected_cross:.1f}")
        ratio_cross = cross_matches / expected_cross if expected_cross > 0 else 0
        print(f"      Ratio (observed/expected): {ratio_cross:.2f}x")
        if best_cross_pair:
            i, j = best_cross_pair
            print(f"      Best cross near-collision HW = {best_cross_hw}")
            print(f"      Per-register: {per_register_hw(hashes_orig[i], hashes_flip[j])}")

        # --- Method E: Full birthday — find closest pair among ALL hashes ---
        # Combine orig and flip hashes, sort by concatenated partial keys
        # across multiple registers
        print(f"\n    Full birthday (closest pair among all 2N hashes):")

        # Use upper 8 bits of registers a,e as 16-bit key
        all_hashes = []
        for i in range(N):
            all_hashes.append(('O', i, hashes_orig[i]))
            all_hashes.append(('F', i, hashes_flip[i]))

        # Sort by combined partial key (upper 8 bits of a + upper 8 bits of e)
        def sort_key(entry):
            h = entry[2]
            return ((h[0] >> 24) << 8) | (h[4] >> 24)

        all_hashes.sort(key=sort_key)

        best_bday_hw = 256
        best_bday_desc = ""
        matches_16bit = 0
        for i in range(len(all_hashes) - 1):
            if sort_key(all_hashes[i]) == sort_key(all_hashes[i + 1]):
                matches_16bit += 1
                h1 = all_hashes[i][2]
                h2 = all_hashes[i + 1][2]
                d = state_hamming(h1, h2)
                if d < best_bday_hw:
                    best_bday_hw = d
                    t1, idx1 = all_hashes[i][0], all_hashes[i][1]
                    t2, idx2 = all_hashes[i + 1][0], all_hashes[i + 1][1]
                    best_bday_desc = f"{t1}[{idx1}] vs {t2}[{idx2}]"
                    best_bday_per_reg = per_register_hw(h1, h2)

        expected_16 = (2 * N) * (2 * N - 1) / (2 * (1 << 16))
        print(f"      16-bit partial matches: {matches_16bit}")
        print(f"      Expected random: {expected_16:.1f}")
        if best_bday_hw < 256:
            print(f"      Best near-collision HW = {best_bday_hw} ({best_bday_desc})")
            print(f"      Per-register: {best_bday_per_reg}")
        print(f"      Improvement over standard: {best_hw_std - best_bday_hw:+d} bits")

    elapsed = time.time() - t0
    print(f"\n  Phase 4 time: {elapsed:.1f}s")
    return elapsed

# ============================================================================
# Main
# ============================================================================

def main():
    print("╔══════════════════════════════════════════════════════════════════════╗")
    print("║        SHA-256 Near-Collision Finding Machine v2                    ║")
    print("║        Combining differential, filtered, multi-path, birthday       ║")
    print("╚══════════════════════════════════════════════════════════════════════╝")
    print()

    total_start = time.time()

    t1 = phase1()
    t2 = phase2()
    t3 = phase3()
    t4 = phase4()

    total = time.time() - total_start

    print("\n" + "=" * 72)
    print("SUMMARY")
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
