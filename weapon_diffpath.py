#!/usr/bin/env python3
"""
Optimal differential path search for SHA-256.
Studies A-D: single-bit, 2-bit, local collision, and best path verification.
Target runtime: under 90 seconds.
"""

import time
import random

# ---------- SHA-256 constants ----------

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

M = 0xFFFFFFFF
IV = (0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19)

K0, K1, K2, K3, K4, K5 = K[0], K[1], K[2], K[3], K[4], K[5]
IV_a, IV_b, IV_c, IV_d, IV_e, IV_f, IV_g, IV_h = IV

def hw(x):
    return bin(x).count('1')

def _rounds(W, nr):
    """SHA-256 rounds. Returns (a,b,c,d,e,f,g,h)."""
    a, b, c, d, e, f, g, h = IV
    for i in range(nr):
        s1 = (((e >> 6) | (e << 26)) ^ ((e >> 11) | (e << 21)) ^ ((e >> 25) | (e << 7))) & M
        cv = (e & f) ^ ((e ^ M) & g)
        t1 = (h + s1 + cv + K[i] + W[i]) & M
        s0 = (((a >> 2) | (a << 30)) ^ ((a >> 13) | (a << 19)) ^ ((a >> 22) | (a << 10))) & M
        mv = (a & b) ^ (a & c) ^ (b & c)
        t2 = (s0 + mv) & M
        h = g; g = f; f = e; e = (d + t1) & M
        d = c; c = b; b = a; a = (t1 + t2) & M
    return (a, b, c, d, e, f, g, h)

_rng = random.Random(42)

def _make_w_flat(count, length):
    """Flat list of random words."""
    r = _rng.randint
    return [r(0, M) for _ in range(count * length)]


# ============================================================
# Core: run N trial pairs, W[0] xored by delta, measure stats
# Returns (mean_total_hw, p_e_le3, p_e_le5)
# ============================================================

def _scan_delta_w0(delta, nr, N, bulk, wlen):
    """Given pre-generated bulk random words, test delta on W[0]."""
    e_le3 = 0
    e_le5 = 0
    thw_sum = 0
    pc = bin  # local ref

    for t in range(N):
        off = t * wlen
        W0 = bulk[off]; W1 = bulk[off+1]
        W2_ = bulk[off+2] if wlen > 2 else 0
        W3_ = bulk[off+3] if wlen > 3 else 0
        W4_ = bulk[off+4] if wlen > 4 else 0
        W5_ = bulk[off+5] if wlen > 5 else 0

        # Original path
        a1=IV_a; b1=IV_b; c1=IV_c; d1=IV_d; e1=IV_e; f1=IV_f; g1=IV_g; h1=IV_h
        # Round 0
        s1v = (((e1>>6)|(e1<<26))^((e1>>11)|(e1<<21))^((e1>>25)|(e1<<7)))&M
        cv = (e1&f1)^((e1^M)&g1)
        t1 = (h1+s1v+cv+K0+W0)&M
        s0v = (((a1>>2)|(a1<<30))^((a1>>13)|(a1<<19))^((a1>>22)|(a1<<10)))&M
        mv = (a1&b1)^(a1&c1)^(b1&c1)
        t2 = (s0v+mv)&M
        h1=g1;g1=f1;f1=e1;e1=(d1+t1)&M;d1=c1;c1=b1;b1=a1;a1=(t1+t2)&M
        if nr > 1:
            s1v=(((e1>>6)|(e1<<26))^((e1>>11)|(e1<<21))^((e1>>25)|(e1<<7)))&M
            cv=(e1&f1)^((e1^M)&g1);t1=(h1+s1v+cv+K1+W1)&M
            s0v=(((a1>>2)|(a1<<30))^((a1>>13)|(a1<<19))^((a1>>22)|(a1<<10)))&M
            mv=(a1&b1)^(a1&c1)^(b1&c1);t2=(s0v+mv)&M
            h1=g1;g1=f1;f1=e1;e1=(d1+t1)&M;d1=c1;c1=b1;b1=a1;a1=(t1+t2)&M
        if nr > 2:
            s1v=(((e1>>6)|(e1<<26))^((e1>>11)|(e1<<21))^((e1>>25)|(e1<<7)))&M
            cv=(e1&f1)^((e1^M)&g1);t1=(h1+s1v+cv+K2+W2_)&M
            s0v=(((a1>>2)|(a1<<30))^((a1>>13)|(a1<<19))^((a1>>22)|(a1<<10)))&M
            mv=(a1&b1)^(a1&c1)^(b1&c1);t2=(s0v+mv)&M
            h1=g1;g1=f1;f1=e1;e1=(d1+t1)&M;d1=c1;c1=b1;b1=a1;a1=(t1+t2)&M
        if nr > 3:
            s1v=(((e1>>6)|(e1<<26))^((e1>>11)|(e1<<21))^((e1>>25)|(e1<<7)))&M
            cv=(e1&f1)^((e1^M)&g1);t1=(h1+s1v+cv+K3+W3_)&M
            s0v=(((a1>>2)|(a1<<30))^((a1>>13)|(a1<<19))^((a1>>22)|(a1<<10)))&M
            mv=(a1&b1)^(a1&c1)^(b1&c1);t2=(s0v+mv)&M
            h1=g1;g1=f1;f1=e1;e1=(d1+t1)&M;d1=c1;c1=b1;b1=a1;a1=(t1+t2)&M
        if nr > 4:
            s1v=(((e1>>6)|(e1<<26))^((e1>>11)|(e1<<21))^((e1>>25)|(e1<<7)))&M
            cv=(e1&f1)^((e1^M)&g1);t1=(h1+s1v+cv+K4+W4_)&M
            s0v=(((a1>>2)|(a1<<30))^((a1>>13)|(a1<<19))^((a1>>22)|(a1<<10)))&M
            mv=(a1&b1)^(a1&c1)^(b1&c1);t2=(s0v+mv)&M
            h1=g1;g1=f1;f1=e1;e1=(d1+t1)&M;d1=c1;c1=b1;b1=a1;a1=(t1+t2)&M
        if nr > 5:
            s1v=(((e1>>6)|(e1<<26))^((e1>>11)|(e1<<21))^((e1>>25)|(e1<<7)))&M
            cv=(e1&f1)^((e1^M)&g1);t1=(h1+s1v+cv+K5+W5_)&M
            s0v=(((a1>>2)|(a1<<30))^((a1>>13)|(a1<<19))^((a1>>22)|(a1<<10)))&M
            mv=(a1&b1)^(a1&c1)^(b1&c1);t2=(s0v+mv)&M
            h1=g1;g1=f1;f1=e1;e1=(d1+t1)&M;d1=c1;c1=b1;b1=a1;a1=(t1+t2)&M

        # Modified path - only W[0] differs
        W0m = W0 ^ delta
        a2=IV_a; b2=IV_b; c2=IV_c; d2=IV_d; e2=IV_e; f2=IV_f; g2=IV_g; h2=IV_h
        s1v = (((e2>>6)|(e2<<26))^((e2>>11)|(e2<<21))^((e2>>25)|(e2<<7)))&M
        cv = (e2&f2)^((e2^M)&g2)
        t1 = (h2+s1v+cv+K0+W0m)&M
        s0v = (((a2>>2)|(a2<<30))^((a2>>13)|(a2<<19))^((a2>>22)|(a2<<10)))&M
        mv = (a2&b2)^(a2&c2)^(b2&c2)
        t2 = (s0v+mv)&M
        h2=g2;g2=f2;f2=e2;e2=(d2+t1)&M;d2=c2;c2=b2;b2=a2;a2=(t1+t2)&M
        if nr > 1:
            s1v=(((e2>>6)|(e2<<26))^((e2>>11)|(e2<<21))^((e2>>25)|(e2<<7)))&M
            cv=(e2&f2)^((e2^M)&g2);t1=(h2+s1v+cv+K1+W1)&M
            s0v=(((a2>>2)|(a2<<30))^((a2>>13)|(a2<<19))^((a2>>22)|(a2<<10)))&M
            mv=(a2&b2)^(a2&c2)^(b2&c2);t2=(s0v+mv)&M
            h2=g2;g2=f2;f2=e2;e2=(d2+t1)&M;d2=c2;c2=b2;b2=a2;a2=(t1+t2)&M
        if nr > 2:
            s1v=(((e2>>6)|(e2<<26))^((e2>>11)|(e2<<21))^((e2>>25)|(e2<<7)))&M
            cv=(e2&f2)^((e2^M)&g2);t1=(h2+s1v+cv+K2+W2_)&M
            s0v=(((a2>>2)|(a2<<30))^((a2>>13)|(a2<<19))^((a2>>22)|(a2<<10)))&M
            mv=(a2&b2)^(a2&c2)^(b2&c2);t2=(s0v+mv)&M
            h2=g2;g2=f2;f2=e2;e2=(d2+t1)&M;d2=c2;c2=b2;b2=a2;a2=(t1+t2)&M
        if nr > 3:
            s1v=(((e2>>6)|(e2<<26))^((e2>>11)|(e2<<21))^((e2>>25)|(e2<<7)))&M
            cv=(e2&f2)^((e2^M)&g2);t1=(h2+s1v+cv+K3+W3_)&M
            s0v=(((a2>>2)|(a2<<30))^((a2>>13)|(a2<<19))^((a2>>22)|(a2<<10)))&M
            mv=(a2&b2)^(a2&c2)^(b2&c2);t2=(s0v+mv)&M
            h2=g2;g2=f2;f2=e2;e2=(d2+t1)&M;d2=c2;c2=b2;b2=a2;a2=(t1+t2)&M
        if nr > 4:
            s1v=(((e2>>6)|(e2<<26))^((e2>>11)|(e2<<21))^((e2>>25)|(e2<<7)))&M
            cv=(e2&f2)^((e2^M)&g2);t1=(h2+s1v+cv+K4+W4_)&M
            s0v=(((a2>>2)|(a2<<30))^((a2>>13)|(a2<<19))^((a2>>22)|(a2<<10)))&M
            mv=(a2&b2)^(a2&c2)^(b2&c2);t2=(s0v+mv)&M
            h2=g2;g2=f2;f2=e2;e2=(d2+t1)&M;d2=c2;c2=b2;b2=a2;a2=(t1+t2)&M
        if nr > 5:
            s1v=(((e2>>6)|(e2<<26))^((e2>>11)|(e2<<21))^((e2>>25)|(e2<<7)))&M
            cv=(e2&f2)^((e2^M)&g2);t1=(h2+s1v+cv+K5+W5_)&M
            s0v=(((a2>>2)|(a2<<30))^((a2>>13)|(a2<<19))^((a2>>22)|(a2<<10)))&M
            mv=(a2&b2)^(a2&c2)^(b2&c2);t2=(s0v+mv)&M
            h2=g2;g2=f2;f2=e2;e2=(d2+t1)&M;d2=c2;c2=b2;b2=a2;a2=(t1+t2)&M

        he = pc(e1^e2).count('1')
        thw = (pc(a1^a2).count('1')+pc(b1^b2).count('1')+
               pc(c1^c2).count('1')+pc(d1^d2).count('1')+
               he+pc(f1^f2).count('1')+
               pc(g1^g2).count('1')+pc(h1^h2).count('1'))
        thw_sum += thw
        if he <= 3: e_le3 += 1
        if he <= 5: e_le5 += 1

    return (thw_sum / N, e_le3 / N, e_le5 / N)


# ============================================================
# STUDY A
# ============================================================

def study_a():
    print("=" * 65)
    print("STUDY A: Single-bit differential simulation (32 positions)")
    print("=" * 65)

    N = 30000
    results = {}

    for nr in (4, 6):
        wlen = nr
        bulk = _make_w_flat(N, wlen)
        for bit in range(32):
            mhw, p3, p5 = _scan_delta_w0(1 << bit, nr, N, bulk, wlen)
            results[(bit, nr)] = (mhw, p3, p5)

    ranking = sorted(range(32),
                     key=lambda b: (-results[(b,4)][1], -results[(b,4)][2], results[(b,4)][0]))

    print(f"\nTrials per position: {N}")
    print(f"\nTop 10 bit positions ranked by P(HW_e <= 3) at round 4:")
    print(f"{'Bit':>4} {'MeanHW_r4':>10} {'P(e<=3)_r4':>12} {'P(e<=5)_r4':>12} "
          f"{'MeanHW_r6':>10} {'P(e<=3)_r6':>12}")
    print("-" * 70)
    for bit in ranking[:10]:
        mh4, p34, p54 = results[(bit, 4)]
        mh6, p36, p56 = results[(bit, 6)]
        print(f"{bit:>4} {mh4:>10.2f} {p34:>12.6f} {p54:>12.6f} "
              f"{mh6:>10.2f} {p36:>12.6f}")

    print(f"\nAll 32 positions at round 4 (sorted):")
    print(f"{'Bit':>4} {'MeanHW':>8} {'P(e<=3)':>10} {'P(e<=5)':>10}")
    print("-" * 36)
    for bit in ranking:
        mh, p3, p5 = results[(bit, 4)]
        print(f"{bit:>4} {mh:>8.2f} {p3:>10.6f} {p5:>10.6f}")

    return ranking, results


# ============================================================
# STUDY B: Best 2-bit differences
# ============================================================

def study_b():
    print("\n" + "=" * 65)
    print("STUDY B: Best 2-bit W[0] differences (496 pairs, 4 rounds)")
    print("=" * 65)

    # Phase 1: coarse scan with 3000 trials
    N1 = 3000
    nr = 4; wlen = nr
    bulk = _make_w_flat(N1, wlen)

    scores = []
    for i in range(32):
        for j in range(i + 1, 32):
            delta = (1 << i) | (1 << j)
            mhw, p3, p5 = _scan_delta_w0(delta, nr, N1, bulk, wlen)
            scores.append((p5, p3, mhw, i, j))

    scores.sort(reverse=True)

    # Phase 2: refine top 50 with 10000 trials
    N2 = 10000
    bulk2 = _make_w_flat(N2, wlen)
    refined = []
    for _, _, _, i, j in scores[:50]:
        delta = (1 << i) | (1 << j)
        mhw, p3, p5 = _scan_delta_w0(delta, nr, N2, bulk2, wlen)
        refined.append((p5, p3, mhw, i, j))

    refined.sort(reverse=True)

    print(f"\nCoarse scan: {N1} trials/pair, refined top 50 with {N2} trials")
    print(f"\nTop 10 two-bit differences by P(HW_e_r4 <= 5):")
    print(f"{'Bits':>10} {'Delta(hex)':>12} {'P(e<=5)':>10} {'P(e<=3)':>10} {'MeanHW':>8}")
    print("-" * 55)
    for p5, p3, mhw, i, j in refined[:10]:
        delta = (1 << i) | (1 << j)
        print(f"  ({i:>2},{j:>2})   0x{delta:08x} {p5:>10.6f} {p3:>10.6f} {mhw:>8.2f}")

    return refined[:50]


# ============================================================
# STUDY C: Local collision via hill climbing
# ============================================================

def study_c():
    print("\n" + "=" * 65)
    print("STUDY C: Local collision via hill climbing")
    print("         DeltaW[0] = bit 31, searching best DeltaW[1]")
    print("=" * 65)

    delta_w0 = 1 << 31
    n_seeds = 300
    n_steps = 40
    eval_N = 100

    rng = random.Random(123)

    def evaluate(dw1, nr=2, N=100):
        total = 0
        for _ in range(N):
            W = [rng.randint(0, M), rng.randint(0, M)]
            if nr > 2:
                W += [rng.randint(0, M) for _ in range(nr - 2)]
            W2 = list(W)
            W2[0] ^= delta_w0
            W2[1] ^= dw1
            s1 = _rounds(W, nr)
            s2 = _rounds(W2, nr)
            for r in range(8):
                total += bin(s1[r] ^ s2[r]).count('1')
        return total / N

    best_dw1 = 0
    best_score = evaluate(0, N=500)

    print(f"\nBaseline (DeltaW[1]=0): mean state diff HW = {best_score:.2f}")
    print(f"Running {n_seeds} seeds x {n_steps} hill-climb steps...")

    for _ in range(n_seeds):
        dw1 = rng.randint(0, M)
        score = evaluate(dw1, N=eval_N)

        for __ in range(n_steps):
            flip = rng.randint(0, 31)
            cand = dw1 ^ (1 << flip)
            cs = evaluate(cand, N=eval_N)
            if cs < score:
                dw1 = cand
                score = cs

        if score < best_score:
            best_score = score
            best_dw1 = dw1

    # Verify
    verified = evaluate(best_dw1, N=5000)

    print(f"\nBest DeltaW[1] found: 0x{best_dw1:08x} (HW={hw(best_dw1)})")
    print(f"  Hill-climb score:  {best_score:.2f}")
    print(f"  Verified (5000t):  {verified:.2f}")

    # Per-register breakdown
    reg_hw = [0]*8
    p_le3 = [0]*8
    p_zero = [0]*8
    Nv = 5000
    for _ in range(Nv):
        W = [rng.randint(0, M), rng.randint(0, M)]
        W2 = list(W); W2[0] ^= delta_w0; W2[1] ^= best_dw1
        s1 = _rounds(W, 2); s2 = _rounds(W2, 2)
        for r in range(8):
            h = bin(s1[r]^s2[r]).count('1')
            reg_hw[r] += h
            if h <= 3: p_le3[r] += 1
            if h == 0: p_zero[r] += 1

    labels = ['a','b','c','d','e','f','g','h']
    print(f"\n  Mean HW per register after 2 rounds ({Nv} trials):")
    for idx, lbl in enumerate(labels):
        print(f"    {lbl}: mean={reg_hw[idx]/Nv:.2f}, "
              f"P(<=3)={p_le3[idx]/Nv:.4f}, "
              f"P(=0)={p_zero[idx]/Nv:.4f}")

    score_3r = evaluate(best_dw1, nr=3, N=2000)
    print(f"\n  Extended to 3 rounds: mean state diff HW = {score_3r:.2f}")

    return best_dw1, verified


# ============================================================
# STUDY D: Best path summary & verification
# ============================================================

def study_d(ranking_a, results_a, top_2bit, best_dw1_c, score_c):
    print("\n" + "=" * 65)
    print("STUDY D: Best path summary & verification (50000 trials)")
    print("=" * 65)

    N = 50000
    best_1bit = ranking_a[0]
    delta_1 = 1 << best_1bit
    _, _, _, bi, bj = top_2bit[0]
    delta_2 = (1 << bi) | (1 << bj)
    delta_w0_c = 1 << 31

    candidates = [
        (f"1-bit W[0] bit {best_1bit}", [(0, delta_1)]),
        (f"2-bit W[0] bits ({bi},{bj})", [(0, delta_2)]),
        (f"Local collision: dW0=bit31, dW1=0x{best_dw1_c:08x}",
         [(0, delta_w0_c), (1, best_dw1_c)]),
    ]

    for desc, diffs in candidates:
        print(f"\n--- {desc} ---")
        max_idx = max(idx for idx, _ in diffs) + 1
        for nr in [2, 3, 4, 5, 6]:
            wlen = max(max_idx, nr)
            bulk = _make_w_flat(N, wlen)

            thw_sum = 0; e_le3 = 0; e_le5 = 0; t_le10 = 0
            for t in range(N):
                off = t * wlen
                W = bulk[off:off+wlen]
                W2 = list(W)
                for idx, dv in diffs:
                    W2[idx] ^= dv
                s1 = _rounds(W, nr)
                s2 = _rounds(W2, nr)
                thw = 0
                for r in range(8):
                    thw += bin(s1[r]^s2[r]).count('1')
                he = bin(s1[4]^s2[4]).count('1')
                thw_sum += thw
                if he <= 3: e_le3 += 1
                if he <= 5: e_le5 += 1
                if thw <= 10: t_le10 += 1

            print(f"  Round {nr}: MeanHW={thw_sum/N:>7.2f}  "
                  f"P(e<=3)={e_le3/N:.6f}  P(e<=5)={e_le5/N:.6f}  "
                  f"P(totalHW<=10)={t_le10/N:.6f}")

    # Final summary
    print("\n" + "=" * 65)
    print("FINAL SUMMARY")
    print("=" * 65)
    mh4, p34, p54 = results_a[(best_1bit, 4)]
    print(f"\nBest single-bit position: bit {best_1bit}")
    print(f"  Round 4: P(HW_e<=3)={p34:.6f}, P(HW_e<=5)={p54:.6f}, MeanHW={mh4:.1f}")

    print(f"\nBest 2-bit difference: bits ({bi},{bj}), delta=0x{delta_2:08x}")
    p5, p3, mhw, _, _ = top_2bit[0]
    print(f"  Round 4: P(HW_e<=5)={p5:.6f}, P(HW_e<=3)={p3:.6f}, MeanHW={mhw:.2f}")

    print(f"\nBest local collision (2-step):")
    print(f"  DeltaW[0]=0x{delta_w0_c:08x} (bit 31)")
    print(f"  DeltaW[1]=0x{best_dw1_c:08x} (HW={hw(best_dw1_c)})")
    print(f"  Mean state diff HW after 2 rounds: {score_c:.2f}")

    print(f"\nKey insight: SHA-256 diffusion is extremely fast.")
    print(f"  After 4 rounds, even the best single-bit diff yields")
    print(f"  mean state HW ~ {mh4:.1f} (out of 256 bits).")
    print(f"  Practical differential attacks on full SHA-256 remain infeasible.")


# ============================================================
# MAIN
# ============================================================

def main():
    t0 = time.time()
    print("SHA-256 Differential Path Search")
    print(f"Started at {time.strftime('%H:%M:%S')}\n")

    ta = time.time()
    ranking_a, results_a = study_a()
    print(f"\n[Study A completed in {time.time()-ta:.1f}s]")

    tb = time.time()
    top_2bit = study_b()
    print(f"\n[Study B completed in {time.time()-tb:.1f}s]")

    tc = time.time()
    best_dw1, score_c = study_c()
    print(f"\n[Study C completed in {time.time()-tc:.1f}s]")

    td = time.time()
    study_d(ranking_a, results_a, top_2bit, best_dw1, score_c)
    print(f"\n[Study D completed in {time.time()-td:.1f}s]")

    total = time.time() - t0
    print(f"\n{'='*65}")
    print(f"Total runtime: {total:.1f}s")
    print(f"{'='*65}")

if __name__ == "__main__":
    main()
