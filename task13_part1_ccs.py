#!/usr/bin/env python3
"""
ЗАДАНИЕ 13, Часть 1: CCS Pattern Clusters
Generate 5000 Wang pairs (birthday on W[0]),
extract 896-bit carry fingerprints (rounds 17-20),
cluster by Hamming distance, check if similar fingerprints → similar hashes.
"""
import random, time, math
from collections import defaultdict

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
def rw(): return random.getrandbits(32)

def schedule(M):
    W = list(M[:16])
    for i in range(16, 64):
        W.append(add(_sig1(W[i-2]), W[i-7], _sig0(W[i-15]), W[i-16]))
    return W

def wang_chain(Wn):
    Wf = list(Wn); Wf[0] = Wn[0] ^ 0x8000
    sn = list(H0); sf = list(H0)
    T1n = add(sn[7], Sig1(sn[4]), Ch(sn[4], sn[5], sn[6]), K[0], Wn[0])
    T2n = add(Sig0(sn[0]), Maj(sn[0], sn[1], sn[2]))
    T1f = add(sf[7], Sig1(sf[4]), Ch(sf[4], sf[5], sf[6]), K[0], Wf[0])
    T2f = add(Sig0(sf[0]), Maj(sf[0], sf[1], sf[2]))
    sn = [add(T1n, T2n), sn[0], sn[1], sn[2], add(sn[3], T1n), sn[4], sn[5], sn[6]]
    sf = [add(T1f, T2f), sf[0], sf[1], sf[2], add(sf[3], T1f), sf[4], sf[5], sf[6]]
    for r in range(1, 16):
        dd = (sf[3] - sn[3]) & MASK
        dh = (sf[7] - sn[7]) & MASK
        dS = (Sig1(sf[4]) - Sig1(sn[4])) & MASK
        dC = (Ch(sf[4], sf[5], sf[6]) - Ch(sn[4], sn[5], sn[6])) & MASK
        dWr = (-(dd + dh + dS + dC)) & MASK
        Wf[r] = add(Wn[r], dWr)
        T1n = add(sn[7], Sig1(sn[4]), Ch(sn[4], sn[5], sn[6]), K[r], Wn[r])
        T2n = add(Sig0(sn[0]), Maj(sn[0], sn[1], sn[2]))
        T1f = add(sf[7], Sig1(sf[4]), Ch(sf[4], sf[5], sf[6]), K[r], Wf[r])
        T2f = add(Sig0(sf[0]), Maj(sf[0], sf[1], sf[2]))
        sn = [add(T1n, T2n), sn[0], sn[1], sn[2], add(sn[3], T1n), sn[4], sn[5], sn[6]]
        sf = [add(T1f, T2f), sf[0], sf[1], sf[2], add(sf[3], T1f), sf[4], sf[5], sf[6]]
    return Wf

def sha256_hash(W):
    s = list(H0)
    Ws = schedule(W)
    for r in range(64):
        T1 = add(s[7], Sig1(s[4]), Ch(s[4], s[5], s[6]), K[r], Ws[r])
        T2 = add(Sig0(s[0]), Maj(s[0], s[1], s[2]))
        s = [add(T1, T2), s[0], s[1], s[2], add(s[3], T1), s[4], s[5], s[6]]
    return [add(H0[i], s[i]) for i in range(8)]

def compute_carries(x, y):
    """Carry vector for x + y mod 2^32."""
    carries = 0; c = 0
    for i in range(32):
        xi = (x >> i) & 1; yi = (y >> i) & 1
        c = (xi & yi) | (xi & c) | (yi & c)
        if c: carries |= (1 << i)
    return carries

def extract_carry_fingerprint(Wn, Wf, target_rounds=(17, 18, 19, 20)):
    """Extract carry_xor for 7 sub-additions on each target round.
    Returns list of 32-bit XOR carry values (7 per round, 4 rounds = 28 values = 896 bits)."""
    Wn_s = schedule(Wn); Wf_s = schedule(Wf)
    sn = list(H0); sf = list(H0)

    fingerprint = []  # list of 32-bit carry_xor values

    for r in range(max(target_rounds) + 1):
        if r in target_rounds:
            # 7 sub-additions in T1 and T2 computation:
            # T1 = h + Sig1(e) + Ch(e,f,g) + K[r] + W[r]
            #   sub1: h + Sig1(e)
            #   sub2: (h+Sig1(e)) + Ch(e,f,g)
            #   sub3: (...) + K[r]
            #   sub4: (...) + W[r]  → T1
            # T2 = Sig0(a) + Maj(a,b,c)
            #   sub5: Sig0(a) + Maj(a,b,c) → T2
            # new_a = T1 + T2
            #   sub6: T1 + T2
            # new_e = d + T1
            #   sub7: d + T1

            s1n = Sig1(sn[4]); s1f = Sig1(sf[4])
            chn = Ch(sn[4], sn[5], sn[6]); chf = Ch(sf[4], sf[5], sf[6])
            s0n = Sig0(sn[0]); s0f = Sig0(sf[0])
            mjn = Maj(sn[0], sn[1], sn[2]); mjf = Maj(sf[0], sf[1], sf[2])

            # sub1: h + Sig1(e)
            c1n = compute_carries(sn[7], s1n); c1f = compute_carries(sf[7], s1f)
            sum1n = add(sn[7], s1n); sum1f = add(sf[7], s1f)

            # sub2: sum1 + Ch
            c2n = compute_carries(sum1n, chn); c2f = compute_carries(sum1f, chf)
            sum2n = add(sum1n, chn); sum2f = add(sum1f, chf)

            # sub3: sum2 + K[r]
            c3n = compute_carries(sum2n, K[r]); c3f = compute_carries(sum2f, K[r])
            sum3n = add(sum2n, K[r]); sum3f = add(sum2f, K[r])

            # sub4: sum3 + W[r] → T1
            c4n = compute_carries(sum3n, Wn_s[r]); c4f = compute_carries(sum3f, Wf_s[r])
            T1n = add(sum3n, Wn_s[r]); T1f = add(sum3f, Wf_s[r])

            # sub5: Sig0(a) + Maj → T2
            c5n = compute_carries(s0n, mjn); c5f = compute_carries(s0f, mjf)
            T2n = add(s0n, mjn); T2f = add(s0f, mjf)

            # sub6: T1 + T2 → new_a
            c6n = compute_carries(T1n, T2n); c6f = compute_carries(T1f, T2f)

            # sub7: d + T1 → new_e
            c7n = compute_carries(sn[3], T1n); c7f = compute_carries(sf[3], T1f)

            carry_xors = [c1n ^ c1f, c2n ^ c2f, c3n ^ c3f, c4n ^ c4f,
                          c5n ^ c5f, c6n ^ c6f, c7n ^ c7f]
            fingerprint.extend(carry_xors)

        # Advance state
        T1n = add(sn[7], Sig1(sn[4]), Ch(sn[4], sn[5], sn[6]), K[r], Wn_s[r])
        T2n = add(Sig0(sn[0]), Maj(sn[0], sn[1], sn[2]))
        T1f = add(sf[7], Sig1(sf[4]), Ch(sf[4], sf[5], sf[6]), K[r], Wf_s[r])
        T2f = add(Sig0(sf[0]), Maj(sf[0], sf[1], sf[2]))
        sn = [add(T1n, T2n), sn[0], sn[1], sn[2], add(sn[3], T1n), sn[4], sn[5], sn[6]]
        sf = [add(T1f, T2f), sf[0], sf[1], sf[2], add(sf[3], T1f), sf[4], sf[5], sf[6]]

    return fingerprint  # 28 x 32-bit values

def fingerprint_hamming(fp1, fp2):
    """Hamming distance between two fingerprints (list of 32-bit values)."""
    d = 0
    for a, b in zip(fp1, fp2):
        d += hw(a ^ b)
    return d

def main():
    random.seed(42)
    print("=" * 70)
    print("ЗАДАНИЕ 13, Часть 1: CCS PATTERN CLUSTERS")
    print("=" * 70)

    # ── Generate 5000 Wang pairs ──
    N_PAIRS = 5000
    print(f"\nGenerating {N_PAIRS} Wang pairs...")
    t0 = time.time()

    pairs = []  # list of (Wn, Wf, hash_n, hash_f, hw_delta_hash)
    fingerprints = []  # list of carry fingerprints

    for i in range(N_PAIRS):
        Wn = [rw() for _ in range(16)]
        Wf = wang_chain(Wn)
        Hn = sha256_hash(schedule(Wn)[:16])  # Need raw message
        Hf = sha256_hash(schedule(Wf)[:16])

        # Actually we should hash the message words directly
        Hn = sha256_hash(Wn)
        Hf = sha256_hash(Wf)

        hw_dh = sum(hw(Hn[j] ^ Hf[j]) for j in range(8))
        pairs.append((Wn, Wf, Hn, Hf, hw_dh))

        fp = extract_carry_fingerprint(Wn, Wf)
        fingerprints.append(fp)

        if (i + 1) % 1000 == 0:
            print(f"  Generated {i+1}/{N_PAIRS} ({time.time()-t0:.1f}s)")

    elapsed = time.time() - t0
    print(f"  Done in {elapsed:.1f}s")

    # ── Stats on HW(δhash) ──
    all_hw = [p[4] for p in pairs]
    mean_hw = sum(all_hw) / len(all_hw)
    print(f"\nOverall mean HW(δhash) = {mean_hw:.2f} (expected ~128)")

    # ── Find 20 closest fingerprint pairs ──
    print(f"\nFinding 20 closest fingerprint pairs (by Hamming distance)...")
    t0 = time.time()

    # Since N=5000, pairwise = 12.5M comparisons. Too many for pure Python.
    # Subsample: compare each pair to 500 random others, track top-20 minimum.

    # Actually let's be smarter: use a heap of size 20
    import heapq

    # We need minimum distances. Use a max-heap of size 20 (negate distances).
    top20 = []  # max-heap: (-distance, i, j)

    # For 5000 pairs, full pairwise = 12.5M. Each comparison is ~28 XOR+popcount.
    # In Python this might take a while. Let's do first 2000 pairs (2M comparisons).
    N_COMPARE = min(N_PAIRS, 2000)
    print(f"  Comparing {N_COMPARE} pairs ({N_COMPARE*(N_COMPARE-1)//2} comparisons)...")

    for i in range(N_COMPARE):
        for j in range(i + 1, N_COMPARE):
            d = fingerprint_hamming(fingerprints[i], fingerprints[j])
            if len(top20) < 20:
                heapq.heappush(top20, (-d, i, j))
            elif d < -top20[0][0]:
                heapq.heapreplace(top20, (-d, i, j))

        if (i + 1) % 500 == 0:
            print(f"    {i+1}/{N_COMPARE} ({time.time()-t0:.1f}s)")

    elapsed = time.time() - t0
    print(f"  Done in {elapsed:.1f}s")

    # Sort by distance (ascending)
    top20_sorted = sorted([(-d, i, j) for d, i, j in top20])

    print(f"\n{'─'*70}")
    print(f"Top-20 closest carry fingerprint pairs:")
    print(f"{'─'*70}")
    print(f"{'Rank':>4} | {'Pair':>12} | {'FP dist':>8} | {'HW(δhash_A)':>12} | {'HW(δhash_B)':>12}")
    print(f"{'─'*4}-+-{'─'*12}-+-{'─'*8}-+-{'─'*12}-+-{'─'*12}")

    close_hw_list = []
    for rank, (dist, i, j) in enumerate(top20_sorted):
        hw_i = pairs[i][4]
        hw_j = pairs[j][4]
        close_hw_list.append(hw_i)
        close_hw_list.append(hw_j)
        print(f"{rank+1:4d} | ({i:4d},{j:4d}) | {dist:8d} | {hw_i:12d} | {hw_j:12d}")

    mean_close = sum(close_hw_list) / len(close_hw_list)
    print(f"\nMean HW(δhash) for close-fingerprint pairs: {mean_close:.2f}")
    print(f"Overall mean HW(δhash): {mean_hw:.2f}")
    print(f"Difference: {mean_close - mean_hw:.2f}")

    # ── Null control: 20 random pairs, their fingerprint distances and HW ──
    print(f"\n{'─'*70}")
    print("NULL CONTROL: 20 random pairs")
    print(f"{'─'*70}")

    rand_indices = random.sample(range(N_COMPARE), 40)
    rand_pairs = [(rand_indices[2*k], rand_indices[2*k+1]) for k in range(20)]

    null_fp_dists = []
    null_hw_list = []
    for i, j in rand_pairs:
        d = fingerprint_hamming(fingerprints[i], fingerprints[j])
        null_fp_dists.append(d)
        null_hw_list.append(pairs[i][4])
        null_hw_list.append(pairs[j][4])

    print(f"Mean fingerprint distance (random): {sum(null_fp_dists)/len(null_fp_dists):.1f}")
    print(f"Mean fingerprint distance (top-20 close): {sum(d for d,_,_ in top20_sorted)/len(top20_sorted):.1f}")
    print(f"Mean HW(δhash) (random pairs): {sum(null_hw_list)/len(null_hw_list):.2f}")
    print(f"Mean HW(δhash) (close pairs): {mean_close:.2f}")

    # ── Correlation: fingerprint distance vs HW(δhash) ──
    # Sample 500 random pairs and compute correlation
    print(f"\n{'─'*70}")
    print("Correlation: fingerprint distance vs |HW(δhash_i) - HW(δhash_j)|")
    print(f"{'─'*70}")

    n_corr = 500
    fp_dists = []
    hw_diffs = []
    for _ in range(n_corr):
        i = random.randint(0, N_COMPARE - 1)
        j = random.randint(0, N_COMPARE - 1)
        while i == j:
            j = random.randint(0, N_COMPARE - 1)
        d = fingerprint_hamming(fingerprints[i], fingerprints[j])
        hw_diff = abs(pairs[i][4] - pairs[j][4])
        fp_dists.append(d)
        hw_diffs.append(hw_diff)

    mean_d = sum(fp_dists) / n_corr
    mean_h = sum(hw_diffs) / n_corr
    cov = sum((fp_dists[k] - mean_d) * (hw_diffs[k] - mean_h) for k in range(n_corr)) / n_corr
    std_d = math.sqrt(sum((x - mean_d)**2 for x in fp_dists) / n_corr)
    std_h = math.sqrt(sum((x - mean_h)**2 for x in hw_diffs) / n_corr)
    corr = cov / (std_d * std_h) if std_d > 0 and std_h > 0 else 0

    print(f"Pearson correlation: {corr:.4f}")
    print(f"  (positive = similar fingerprints → similar hash diffs)")

    if abs(corr) > 0.1:
        print(f"  → SIGNIFICANT correlation detected!")
    else:
        print(f"  → No significant correlation (|r| < 0.1)")

    print(f"\n{'='*70}")
    print("Part 1 COMPLETE")
    print(f"{'='*70}")

if __name__ == "__main__":
    main()
