#!/usr/bin/env python3
"""
ЗАДАНИЕ 13, Часть 3: Mutual Information — информационный поток SHA-256
Measure MI between each input bit and specific output bits.
"""
import random, time, math

MASK = 0xFFFFFFFF

K = [
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2,
]
H0 = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]

def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK
def _sig0(x):   return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def _sig1(x):   return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)
def Sig0(x):    return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sig1(x):    return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return ((e & f) ^ ((~e & MASK) & g)) & MASK
def Maj(a, b, c): return ((a & b) ^ (a & c) ^ (b & c)) & MASK
def add(*args):
    s = 0
    for x in args: s = (s + x) & MASK
    return s
def hw(x): return bin(x & MASK).count('1')

def schedule(M):
    W = list(M[:16])
    for i in range(16, 64):
        W.append(add(_sig1(W[i-2]), W[i-7], _sig0(W[i-15]), W[i-16]))
    return W

def sha256_hash(W):
    s = list(H0)
    Ws = schedule(W)
    for r in range(64):
        T1 = add(s[7], Sig1(s[4]), Ch(s[4], s[5], s[6]), K[r], Ws[r])
        T2 = add(Sig0(s[0]), Maj(s[0], s[1], s[2]))
        s = [add(T1, T2), s[0], s[1], s[2], add(s[3], T1), s[4], s[5], s[6]]
    return [add(H0[i], s[i]) for i in range(8)]

def all_states(W):
    s = list(H0)
    Ws = schedule(W)
    states = [tuple(s)]
    for r in range(64):
        T1 = add(s[7], Sig1(s[4]), Ch(s[4], s[5], s[6]), K[r], Ws[r])
        T2 = add(Sig0(s[0]), Maj(s[0], s[1], s[2]))
        s = [add(T1, T2), s[0], s[1], s[2], add(s[3], T1), s[4], s[5], s[6]]
        states.append(tuple(s))
    return states

def compute_mi(counts, n_total):
    """Compute MI from 2x2 contingency table.
    counts = [n00, n01, n10, n11] where first index=input bit, second=output bit."""
    n00, n01, n10, n11 = counts
    mi = 0.0
    for w in range(2):
        for h in range(2):
            p_wh = counts[w * 2 + h] / n_total
            p_w = (counts[w * 2 + 0] + counts[w * 2 + 1]) / n_total
            p_h = (counts[0 * 2 + h] + counts[1 * 2 + h]) / n_total
            if p_wh > 0 and p_w > 0 and p_h > 0:
                mi += p_wh * math.log2(p_wh / (p_w * p_h))
    return mi

def measure_mi_for_output_bit(out_word, out_bit, n_samples=10000):
    """Measure MI between every input bit and one output bit.
    Returns list of (MI, word_idx, bit_idx) for all 512 input bits."""
    random.seed(12345)

    # Pre-generate random messages and their hashes
    messages = []
    hashes = []
    for _ in range(n_samples):
        W = [random.getrandbits(32) for _ in range(16)]
        messages.append(W)
        hashes.append(sha256_hash(W))

    results = []

    for w_idx in range(16):
        for b_idx in range(32):
            # Count (input_bit, output_bit) pairs
            counts = [0, 0, 0, 0]  # [n00, n01, n10, n11]
            for s in range(n_samples):
                input_bit = (messages[s][w_idx] >> b_idx) & 1
                output_bit = (hashes[s][out_word] >> out_bit) & 1
                counts[input_bit * 2 + output_bit] += 1

            mi = compute_mi(counts, n_samples)
            results.append((mi, w_idx, b_idx))

    return results

def measure_mi_midstate(mid_round, mid_word, mid_bit, out_word, out_bit, n_samples=10000):
    """Measure MI between a midstate bit and an output bit."""
    random.seed(12345)

    counts = [0, 0, 0, 0]
    for _ in range(n_samples):
        W = [random.getrandbits(32) for _ in range(16)]
        states = all_states(W)
        H = [add(H0[i], states[64][i]) for i in range(8)]

        mid_val = (states[mid_round][mid_word] >> mid_bit) & 1
        out_val = (H[out_word] >> out_bit) & 1
        counts[mid_val * 2 + out_val] += 1

    return compute_mi(counts, n_samples)

def main():
    print("=" * 70)
    print("ЗАДАНИЕ 13, Часть 3: MUTUAL INFORMATION")
    print("=" * 70)

    # ── Three output bits to analyze ──
    output_bits = [
        (0, 0,  "H[0][0]  (LSB of first hash word)"),
        (0, 31, "H[0][31] (MSB of first hash word)"),
        (7, 0,  "H[7][0]  (LSB of last hash word)"),
    ]

    N_SAMPLES = 10000

    all_results = {}

    for out_word, out_bit, label in output_bits:
        print(f"\n{'─'*70}")
        print(f"MI analysis for output bit: {label}")
        print(f"N_samples = {N_SAMPLES}")
        print(f"{'─'*70}")

        t0 = time.time()
        results = measure_mi_for_output_bit(out_word, out_bit, N_SAMPLES)
        elapsed = time.time() - t0
        print(f"  Computed in {elapsed:.1f}s")

        # Sort by MI descending
        results.sort(reverse=True)
        all_results[(out_word, out_bit)] = results

        max_mi = results[0][0]
        print(f"\n  Max MI = {max_mi:.6f}")

        print(f"\n  Top-10 input bits by MI:")
        print(f"  {'Rank':>4} | {'Input bit':>15} | {'MI':>12}")
        print(f"  {'─'*4}-+-{'─'*15}-+-{'─'*12}")
        for rank in range(10):
            mi, w_idx, b_idx = results[rank]
            print(f"  {rank+1:4d} | W[{w_idx:2d}][{b_idx:2d}]      | {mi:12.6f}")

        # Statistics
        all_mi = [r[0] for r in results]
        mean_mi = sum(all_mi) / len(all_mi)
        print(f"\n  Mean MI across all 512 input bits: {mean_mi:.6f}")
        print(f"  Median MI: {sorted(all_mi)[256]:.6f}")
        n_above_001 = sum(1 for m in all_mi if m > 0.001)
        n_above_0005 = sum(1 for m in all_mi if m > 0.0005)
        print(f"  Input bits with MI > 0.001: {n_above_001}")
        print(f"  Input bits with MI > 0.0005: {n_above_0005}")

    # ── Compare top-MI bits across output bits ──
    print(f"\n{'─'*70}")
    print("Cross-comparison: are top-MI input bits the same for different outputs?")
    print(f"{'─'*70}")

    for i, (ow1, ob1, l1) in enumerate(output_bits):
        for j, (ow2, ob2, l2) in enumerate(output_bits):
            if j <= i: continue
            top10_1 = set((r[1], r[2]) for r in all_results[(ow1, ob1)][:10])
            top10_2 = set((r[1], r[2]) for r in all_results[(ow2, ob2)][:10])
            overlap = top10_1 & top10_2
            print(f"  {l1} vs {l2}: {len(overlap)}/10 overlap")
            if overlap:
                print(f"    Common bits: {sorted(overlap)}")

    # ── Verdict ──
    print(f"\n{'─'*70}")
    print("VERDICT:")
    print(f"{'─'*70}")

    max_mi_overall = max(all_results[(ow, ob)][0][0] for ow, ob, _ in output_bits)

    if max_mi_overall > 0.001:
        print(f"  Max MI = {max_mi_overall:.6f} > 0.001")
        print(f"  → Information channel EXISTS")
        print(f"  → RULE 14: Triple-checking required")

        # Triple check with different seeds
        print(f"\n  TRIPLE CHECK:")
        check_results = []
        for seed in [999, 7777, 31415]:
            random.seed(seed)
            # Just check the top input bit for H[0][0]
            top_w, top_b = all_results[(0, 0)][0][1], all_results[(0, 0)][0][2]
            counts = [0, 0, 0, 0]
            for _ in range(N_SAMPLES):
                W = [random.getrandbits(32) for _ in range(16)]
                H = sha256_hash(W)
                inp = (W[top_w] >> top_b) & 1
                out = (H[0] >> 0) & 1
                counts[inp * 2 + out] += 1
            mi = compute_mi(counts, N_SAMPLES)
            check_results.append(mi)
            print(f"    Seed {seed}: MI(W[{top_w}][{top_b}], H[0][0]) = {mi:.6f}")

        mean_check = sum(check_results) / len(check_results)
        print(f"    Mean across checks: {mean_check:.6f}")
        if mean_check > 0.001:
            print(f"    → CONFIRMED: information channel real")
        else:
            print(f"    → REFUTED: initial result was noise")
    else:
        print(f"  Max MI = {max_mi_overall:.6f} < 0.001")
        print(f"  → SHA-256 fully destroys MI over 64 rounds")
        print(f"  → Information flow does NOT provide new attack surface")

    # ── BONUS: MI between midstate and output ──
    print(f"\n{'─'*70}")
    print("BONUS: MI(state[r][0][0], H[0][0]) for various rounds r")
    print(f"{'─'*70}")

    N_MID = 5000  # fewer samples for speed

    for r in [0, 8, 16, 24, 32, 40, 48, 56, 64]:
        random.seed(12345)
        counts = [0, 0, 0, 0]
        for _ in range(N_MID):
            W = [random.getrandbits(32) for _ in range(16)]
            states = all_states(W)
            H = [add(H0[i], states[64][i]) for i in range(8)]

            if r == 64:
                mid_val = (H[0] >> 0) & 1  # final hash = state[64] + H0
            else:
                mid_val = (states[r][0] >> 0) & 1
            out_val = (H[0] >> 0) & 1
            counts[mid_val * 2 + out_val] += 1

        mi = compute_mi(counts, N_MID)
        marker = ""
        if r == 0: marker = " (IV, constant → MI=0)"
        if r == 64: marker = " (output = self → MI=1.0)"
        print(f"  Round {r:2d}: MI = {mi:.6f}{marker}")

    print(f"\n{'='*70}")
    print("Part 3 COMPLETE")
    print(f"{'='*70}")

if __name__ == "__main__":
    main()
