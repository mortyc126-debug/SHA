#!/usr/bin/env python3
"""
Optimal differential path search for SHA-256.
Studies A-D: single-bit, 2-bit, local collision, and best path verification.
Target runtime: under 90 seconds.
"""

import struct
import time
import random
from collections import defaultdict

# ---------- SHA-256 constants & primitives ----------

K = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f118f7, 0x923f82a4, 0xab1c5ed5,
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

M32 = 0xFFFFFFFF

def rotr(x, n):
    return ((x >> n) | (x << (32 - n))) & M32

def sigma1(e):
    return rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25)

def sigma0(a):
    return rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22)

def ch(e, f, g):
    return (e & f) ^ ((~e) & g) & M32

def maj(a, b, c):
    return (a & b) ^ (a & c) ^ (b & c)

def hw(x):
    """Hamming weight via popcount."""
    return bin(x).count('1')

# ---------- SHA-256 partial rounds ----------

# IV for SHA-256
IV = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
]

def sha256_rounds(W_list, n_rounds):
    """Run n_rounds of SHA-256 compression from IV with given W words.
    W_list must have at least n_rounds entries.
    Returns state [a,b,c,d,e,f,g,h]."""
    a, b, c, d, e, f, g, h = IV
    for i in range(n_rounds):
        s1 = sigma1(e)
        ch_val = ch(e, f, g)
        temp1 = (h + s1 + ch_val + K[i] + W_list[i]) & M32
        s0 = sigma0(a)
        maj_val = maj(a, b, c)
        temp2 = (s0 + maj_val) & M32
        h, g, f, e, d, c, b, a = g, f, e, (d + temp1) & M32, c, b, a, (temp1 + temp2) & M32
    return [a, b, c, d, e, f, g, h]

def state_diff_hw(s1, s2):
    """Total HW of XOR diff between two states."""
    return sum(hw(a ^ b) for a, b in zip(s1, s2))

def state_xor(s1, s2):
    """XOR diff of two states."""
    return [a ^ b for a, b in zip(s1, s2)]

# ---------- Fast random word generation ----------

_rand = random.Random(42)
_randint = _rand.randint

def rand_w():
    return _randint(0, M32)

def rand_w_list(n):
    return [_randint(0, M32) for _ in range(n)]

# ============================================================
# STUDY A: Single-bit differential simulation (all 32 positions)
# ============================================================

def study_a():
    print("=" * 65)
    print("STUDY A: Single-bit differential simulation (32 positions)")
    print("=" * 65)

    N = 30000
    n_rounds_list = [4, 6]
    max_rounds = max(n_rounds_list)

    results = {}  # bit_pos -> {round -> {mean_hw, p_hw_e_le3, p_hw_e_le5}}

    for bit in range(32):
        delta = 1 << bit
        results[bit] = {}

        for nr in n_rounds_list:
            total_hw = 0
            count_e_le3 = 0
            count_e_le5 = 0

            for _ in range(N):
                W = rand_w_list(max_rounds)
                W2 = W[:]
                W2[0] ^= delta

                s1 = sha256_rounds(W, nr)
                s2 = sha256_rounds(W2, nr)

                diff = state_xor(s1, s2)
                total = sum(hw(d) for d in diff)
                hw_e = hw(diff[4])  # e is index 4 in [a,b,c,d,e,f,g,h]

                total_hw += total
                if hw_e <= 3:
                    count_e_le3 += 1
                if hw_e <= 5:
                    count_e_le5 += 1

            results[bit][nr] = {
                'mean_hw': total_hw / N,
                'p_e_le3': count_e_le3 / N,
                'p_e_le5': count_e_le5 / N,
            }

    # Rank by P(HW_e <= 3) at round 4
    ranking = sorted(range(32), key=lambda b: -results[b][4]['p_e_le3'])

    print(f"\nTrials per position: {N}")
    print(f"\nTop 10 bit positions ranked by P(HW_e <= 3) at round 4:")
    print(f"{'Bit':>4} {'MeanHW_r4':>10} {'P(e<=3)_r4':>12} {'P(e<=5)_r4':>12} {'MeanHW_r6':>10} {'P(e<=3)_r6':>12}")
    print("-" * 65)
    for bit in ranking[:10]:
        r4 = results[bit][4]
        r6 = results[bit][6]
        print(f"{bit:>4} {r4['mean_hw']:>10.2f} {r4['p_e_le3']:>12.6f} {r4['p_e_le5']:>12.6f} {r6['mean_hw']:>10.2f} {r6['p_e_le3']:>12.6f}")

    print(f"\nAll 32 positions at round 4:")
    print(f"{'Bit':>4} {'MeanHW':>8} {'P(e<=3)':>10} {'P(e<=5)':>10}")
    print("-" * 36)
    for bit in ranking:
        r4 = results[bit][4]
        print(f"{bit:>4} {r4['mean_hw']:>8.2f} {r4['p_e_le3']:>10.6f} {r4['p_e_le5']:>10.6f}")

    return ranking, results

# ============================================================
# STUDY B: Best 2-bit differences (top from C(32,2)=496)
# ============================================================

def study_b():
    print("\n" + "=" * 65)
    print("STUDY B: Best 2-bit W[0] differences (496 pairs, 4 rounds)")
    print("=" * 65)

    N = 10000
    nr = 4
    max_rounds = 6  # need at least nr words

    scores = []

    for i in range(32):
        for j in range(i + 1, 32):
            delta = (1 << i) | (1 << j)
            count_e_le5 = 0
            count_e_le3 = 0
            total_hw = 0

            for _ in range(N):
                W = rand_w_list(max_rounds)
                W2 = W[:]
                W2[0] ^= delta

                s1 = sha256_rounds(W, nr)
                s2 = sha256_rounds(W2, nr)

                diff = state_xor(s1, s2)
                hw_e = hw(diff[4])
                total_hw += sum(hw(d) for d in diff)
                if hw_e <= 3:
                    count_e_le3 += 1
                if hw_e <= 5:
                    count_e_le5 += 1

            scores.append((count_e_le5 / N, count_e_le3 / N, total_hw / N, i, j))

    scores.sort(reverse=True)

    print(f"\nTrials per pair: {N}")
    print(f"\nTop 10 two-bit differences by P(HW_e_r4 <= 5):")
    print(f"{'Bits':>10} {'Delta(hex)':>12} {'P(e<=5)':>10} {'P(e<=3)':>10} {'MeanHW':>8}")
    print("-" * 55)
    for p5, p3, mhw, i, j in scores[:10]:
        delta = (1 << i) | (1 << j)
        print(f"  ({i:>2},{j:>2}) {delta:>12s} {p5:>10.6f} {p3:>10.6f} {mhw:>8.2f}".replace(
            f"{delta:>12s}", f"0x{delta:08x}"))

    # Fix the formatting
    print(f"\nTop 10 (re-formatted):")
    print(f"{'Bits':>10} {'Delta(hex)':>12} {'P(e<=5)':>10} {'P(e<=3)':>10} {'MeanHW':>8}")
    print("-" * 55)
    for p5, p3, mhw, i, j in scores[:10]:
        delta = (1 << i) | (1 << j)
        print(f"  ({i:>2},{j:>2})   0x{delta:08x} {p5:>10.6f} {p3:>10.6f} {mhw:>8.2f}")

    return scores[:50]

# ============================================================
# STUDY C: Local collision via hill climbing
# ============================================================

def study_c():
    print("\n" + "=" * 65)
    print("STUDY C: Local collision via hill climbing")
    print("         DeltaW[0] = bit 31, searching best DeltaW[1]")
    print("=" * 65)

    delta_w0 = 1 << 31
    n_seeds = 1000
    n_steps = 100
    eval_trials = 200  # trials per evaluation for speed

    def evaluate(dw0, dw1, nr=2, trials=200):
        """Average total state diff HW after nr rounds with given W diffs."""
        total = 0
        for _ in range(trials):
            W = rand_w_list(nr)
            W2 = W[:]
            W2[0] ^= dw0
            W2[1] ^= dw1
            s1 = sha256_rounds(W, nr)
            s2 = sha256_rounds(W2, nr)
            total += state_diff_hw(s1, s2)
        return total / trials

    best_dw1 = 0
    best_score = evaluate(delta_w0, 0)

    print(f"\nBaseline (DeltaW[1]=0): mean state diff HW = {best_score:.2f}")
    print(f"Running {n_seeds} random seeds x {n_steps} hill-climb steps...")

    for seed_i in range(n_seeds):
        # Random starting DeltaW[1]
        dw1 = rand_w()
        score = evaluate(delta_w0, dw1, trials=eval_trials)

        for step in range(n_steps):
            # Flip a random bit
            flip_bit = _randint(0, 31)
            candidate = dw1 ^ (1 << flip_bit)
            cand_score = evaluate(delta_w0, candidate, trials=eval_trials)

            if cand_score < score:
                dw1 = candidate
                score = cand_score

        if score < best_score:
            best_score = score
            best_dw1 = dw1

    # Verify best with more trials
    verified_score = evaluate(delta_w0, best_dw1, nr=2, trials=5000)

    print(f"\nBest DeltaW[1] found: 0x{best_dw1:08x} (HW={hw(best_dw1)})")
    print(f"  Hill-climb score:  {best_score:.2f}")
    print(f"  Verified (5000t):  {verified_score:.2f}")

    # Show the state diff distribution for best
    diffs_detail = []
    for _ in range(5000):
        W = rand_w_list(2)
        W2 = W[:]
        W2[0] ^= delta_w0
        W2[1] ^= best_dw1
        s1 = sha256_rounds(W, 2)
        s2 = sha256_rounds(W2, 2)
        d = state_xor(s1, s2)
        diffs_detail.append([hw(x) for x in d])

    labels = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
    print(f"\n  Mean HW per register after 2 rounds (5000 trials):")
    for idx, lbl in enumerate(labels):
        vals = [d[idx] for d in diffs_detail]
        print(f"    {lbl}: mean={sum(vals)/len(vals):.2f}, "
              f"P(<=3)={sum(1 for v in vals if v<=3)/len(vals):.4f}, "
              f"P(=0)={sum(1 for v in vals if v==0)/len(vals):.4f}")

    # Also try extending to 3 rounds
    score_3r = evaluate(delta_w0, best_dw1, nr=3, trials=2000)
    print(f"\n  Extended to 3 rounds: mean state diff HW = {score_3r:.2f}")

    return best_dw1, verified_score

# ============================================================
# STUDY D: Best path summary & verification
# ============================================================

def study_d(ranking_a, results_a, top_2bit, best_dw1_c, score_c):
    print("\n" + "=" * 65)
    print("STUDY D: Best path summary & verification (50000 trials)")
    print("=" * 65)

    N = 50000

    # Best single-bit from Study A
    best_1bit = ranking_a[0]
    delta_1 = 1 << best_1bit

    # Best 2-bit from Study B
    _, _, _, bi, bj = top_2bit[0]
    delta_2 = (1 << bi) | (1 << bj)

    # Best local collision from Study C
    delta_w0_c = 1 << 31

    candidates = [
        (f"1-bit W[0] bit {best_1bit}", [(0, delta_1)]),
        (f"2-bit W[0] bits ({bi},{bj})", [(0, delta_2)]),
        (f"Local collision: dW0=bit31, dW1=0x{best_dw1_c:08x}", [(0, delta_w0_c), (1, best_dw1_c)]),
    ]

    for desc, diffs in candidates:
        print(f"\n--- {desc} ---")
        max_r = 6
        n_w = max(idx for idx, _ in diffs) + 1
        n_w = max(n_w, max_r)

        for nr in [2, 3, 4, 5, 6]:
            total_hw = 0
            count_e_le3 = 0
            count_e_le5 = 0
            count_total_le10 = 0

            for _ in range(N):
                W = rand_w_list(n_w)
                W2 = W[:]
                for idx, dv in diffs:
                    W2[idx] ^= dv

                s1 = sha256_rounds(W, nr)
                s2 = sha256_rounds(W2, nr)

                diff = state_xor(s1, s2)
                hw_total = sum(hw(d) for d in diff)
                hw_e = hw(diff[4])

                total_hw += hw_total
                if hw_e <= 3:
                    count_e_le3 += 1
                if hw_e <= 5:
                    count_e_le5 += 1
                if hw_total <= 10:
                    count_total_le10 += 1

            mean_hw = total_hw / N
            print(f"  Round {nr}: MeanHW={mean_hw:>7.2f}  "
                  f"P(e<=3)={count_e_le3/N:.6f}  "
                  f"P(e<=5)={count_e_le5/N:.6f}  "
                  f"P(totalHW<=10)={count_total_le10/N:.6f}")

    # Final summary
    print("\n" + "=" * 65)
    print("FINAL SUMMARY")
    print("=" * 65)
    print(f"\nBest single-bit position: bit {best_1bit}")
    r4 = results_a[best_1bit][4]
    print(f"  Round 4: P(HW_e<=3)={r4['p_e_le3']:.6f}, P(HW_e<=5)={r4['p_e_le5']:.6f}")

    print(f"\nBest 2-bit difference: bits ({bi},{bj}), delta=0x{delta_2:08x}")
    p5, p3, mhw, _, _ = top_2bit[0]
    print(f"  Round 4: P(HW_e<=5)={p5:.6f}, P(HW_e<=3)={p3:.6f}, MeanHW={mhw:.2f}")

    print(f"\nBest local collision (2-step):")
    print(f"  DeltaW[0]=0x{delta_w0_c:08x} (bit 31)")
    print(f"  DeltaW[1]=0x{best_dw1_c:08x} (HW={hw(best_dw1_c)})")
    print(f"  Mean state diff HW after 2 rounds: {score_c:.2f}")

    print(f"\nKey insight: SHA-256 diffusion is extremely fast.")
    print(f"  After 4 rounds, even the best single-bit diff yields")
    print(f"  mean state HW ~ {results_a[best_1bit][4]['mean_hw']:.1f} (out of 256 bits).")
    print(f"  Practical differential attacks on full SHA-256 remain infeasible.")

# ============================================================
# MAIN
# ============================================================

def main():
    t0 = time.time()

    print("SHA-256 Differential Path Search")
    print(f"Started at {time.strftime('%H:%M:%S')}")
    print()

    # Study A
    ta = time.time()
    ranking_a, results_a = study_a()
    print(f"\n[Study A completed in {time.time()-ta:.1f}s]")

    # Study B
    tb = time.time()
    top_2bit = study_b()
    print(f"\n[Study B completed in {time.time()-tb:.1f}s]")

    # Study C
    tc = time.time()
    best_dw1, score_c = study_c()
    print(f"\n[Study C completed in {time.time()-tc:.1f}s]")

    # Study D
    td = time.time()
    study_d(ranking_a, results_a, top_2bit, best_dw1, score_c)
    print(f"\n[Study D completed in {time.time()-td:.1f}s]")

    total = time.time() - t0
    print(f"\n{'='*65}")
    print(f"Total runtime: {total:.1f}s")
    print(f"{'='*65}")

if __name__ == "__main__":
    main()
