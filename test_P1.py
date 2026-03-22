#!/usr/bin/env python3
"""
TEST P1: SHA-256/20 collision at cost 2^7?

From LAWS.md T3: Cost_real(20) = 2^{7×(20-19)²} = 2^7 = 128 operations.

The theory says: with cascade shift to R=19 base and late differential,
round 20 is the first "paid" barrier, costing only 2^7.

Test plan:
1. For diff in W[12] (best Eff.DOF from E4): try 2^15 random messages
2. For diff in W[9] (most neutrals): try 2^15 random messages
3. For diff in W[15] (lowest first barrier): try 2^15 random messages
4. Measure actual HW at 20 rounds
5. Does ANY reach HW=0?

Also test P2: SHA-256/20 with d=12 at cost 2^15.
"""

import random
import time

MASK = 0xFFFFFFFF

K = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
    0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
]

H0 = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
]

def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK
def Sig0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sig1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def sig0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sig1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)
def Ch(e, f, g): return ((e & f) ^ (~e & g)) & MASK
def Maj(a, b, c): return ((a & b) ^ (a & c) ^ (b & c)) & MASK
def add32(*args):
    s = 0
    for a in args: s = (s + a) & MASK
    return s
def hw(x): return bin(x & MASK).count('1')

def sha256_20(msg, iv):
    W = list(msg)
    for i in range(16, 20):
        W.append(add32(sig1(W[i-2]), W[i-7], sig0(W[i-15]), W[i-16]))
    state = list(iv)
    for t in range(20):
        a, b, c, d, e, f, g, h = state
        T1 = add32(h, Sig1(e), Ch(e, f, g), K[t], W[t])
        T2 = add32(Sig0(a), Maj(a, b, c))
        state = [add32(T1, T2), a, b, c, add32(d, T1), e, f, g]
    return [add32(state[i], iv[i]) for i in range(8)]


def test_delta(diff_word, diff_bit, n_trials, label):
    iv = list(H0)
    delta_val = 1 << diff_bit
    best_hw = 999
    best_msg = None
    best_regs = None

    t0 = time.time()
    for trial in range(n_trials):
        msg = [random.getrandbits(32) for _ in range(16)]
        msg2 = list(msg)
        msg2[diff_word] ^= delta_val

        h1 = sha256_20(msg, iv)
        h2 = sha256_20(msg2, iv)

        total = sum(hw(h1[i] ^ h2[i]) for i in range(8))

        if total < best_hw:
            best_hw = total
            best_msg = msg
            best_regs = [hw(h1[i] ^ h2[i]) for i in range(8)]
            if total <= 20:
                print(f"  [{label}] trial {trial}: HW={total} regs={best_regs}")
            if total == 0:
                print(f"\n  *** COLLISION FOUND! ***")
                print(f"  M  = {[hex(w) for w in msg]}")
                print(f"  M' = {[hex(w) for w in msg2]}")
                print(f"  H  = {[hex(h) for h in h1]}")
                return 0, msg, best_regs

        if trial % 100000 == 0 and trial > 0:
            elapsed = time.time() - t0
            print(f"  [{label}] {trial}/{n_trials} best={best_hw} [{elapsed:.0f}s]")

    elapsed = time.time() - t0
    return best_hw, best_msg, best_regs


def main():
    print("=" * 72)
    print("TEST P1: SHA-256/20 collision search")
    print("=" * 72)
    print()

    random.seed(0x1AAA)
    t_start = time.time()

    N = 500_000  # 2^19 trials per delta

    # ── Test 1: W[15] single-bit deltas (lowest first barrier) ────

    print(f"Test 1: W[15] bit 28 (CRAZY-15 best), {N} trials at R=20")
    print("-" * 72)
    hw1, msg1, regs1 = test_delta(15, 28, N, "W15b28")
    print(f"  Result: best HW = {hw1}, regs = {regs1}")
    print()

    # ── Test 2: W[12] single-bit deltas (best Eff.DOF) ───────────

    if time.time() - t_start < 120:
        print(f"Test 2: W[12] bit 28 (best Eff.DOF from E4), {N} trials")
        print("-" * 72)
        hw2, msg2, regs2 = test_delta(12, 28, N, "W12b28")
        print(f"  Result: best HW = {hw2}, regs = {regs2}")
        print()

    # ── Test 3: W[9] (most neutrals) ─────────────────────────────

    if time.time() - t_start < 240:
        print(f"Test 3: W[9] bit 28 (most neutrals from Z2), {N} trials")
        print("-" * 72)
        hw3, msg3, regs3 = test_delta(9, 28, N, "W9b28")
        print(f"  Result: best HW = {hw3}, regs = {regs3}")
        print()

    # ── Test 4: Scan ALL single-bit deltas in W[10..15] ──────────

    if time.time() - t_start < 300:
        print(f"Test 4: Best single-bit delta scan (W[10..15], all 32 bits)")
        print("-" * 72)
        overall_best = 999
        overall_config = None
        n_scan = 50_000  # per delta

        for word in range(10, 16):
            for bit in [0, 7, 15, 23, 28, 31]:  # sample 6 bits per word
                if time.time() - t_start > 360:
                    break
                hw_s, _, regs_s = test_delta(word, bit, n_scan, f"W{word}b{bit}")
                if hw_s < overall_best:
                    overall_best = hw_s
                    overall_config = (word, bit, regs_s)

        if overall_config:
            print(f"  Best delta: W[{overall_config[0]}] bit {overall_config[1]}")
            print(f"  Best HW: {overall_best}, regs = {overall_config[2]}")
        print()

    # ── Verdict ───────────────────────────────────────────────────

    print("=" * 72)
    print("VERDICT on P1")
    print("=" * 72)
    print()

    results = []
    if 'hw1' in dir(): results.append(("W[15]b28", hw1))
    if 'hw2' in dir(): results.append(("W[12]b28", hw2))
    if 'hw3' in dir(): results.append(("W[9]b28", hw3))

    for label, hw_val in results:
        print(f"  {label}: best HW = {hw_val}")

    best_overall = min(r[1] for r in results) if results else 999

    if best_overall == 0:
        print(f"\n  P1 CONFIRMED: SHA-256/20 collision FOUND!")
    elif best_overall < 20:
        print(f"\n  P1 PARTIAL: near-collision HW={best_overall}.")
        print(f"  Theory predicted 2^7 cost. Actual: > 2^19 and counting.")
        print(f"  Theory OVERESTIMATES — conditions aren't fully independent.")
    else:
        print(f"\n  P1 REFUTED at this budget: HW={best_overall} >> 0.")
        print(f"  Theory's 2^7 prediction was too optimistic.")
        print(f"  Law Z10 (condition growth) needs recalibration.")

    print(f"\n  Runtime: {time.time() - t_start:.1f}s")
    print("=" * 72)


if __name__ == "__main__":
    main()
