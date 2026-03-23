#!/usr/bin/env python3
"""
weapon_localcol_v2.py — 3-step and 4-step local collision search for SHA-256.

Extends the best 2-step finding (DW0=bit31, DW1=0x821001c0, HW=14.23, P(HW_e<=3)=49%)
to 3-step, 4-step, and greedy sequential 8-step local collision paths.

Phases:
  1. 3-step local collision search (fix DW0,DW1, hill-climb DW2; joint DW1+DW2)
  2. 4-step local collision (extend best 3-step, hill-climb DW3)
  3. Greedy sequential extension (DW0..DW7)
  4. Verification with 50000 trials

SHA-256 round function:
  S1 = ROTR(e,6) ^ ROTR(e,11) ^ ROTR(e,25)
  ch = (e & f) ^ (~e & g)
  temp1 = (h + S1 + ch + K[i] + W[i]) mod 2^32
  S0 = ROTR(a,2) ^ ROTR(a,13) ^ ROTR(a,22)
  maj = (a & b) ^ (a & c) ^ (b & c)
  temp2 = (S0 + maj) mod 2^32
  h,g,f,e,d,c,b,a = g,f,e,(d+temp1),c,b,a,(temp1+temp2)
"""

import time
import random

M = 0xFFFFFFFF

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
    0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
]

IV = (0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19)

def hw(x):
    return bin(x).count('1')

def make_pool(n_trials, n_words, seed=42):
    rng = random.Random(seed)
    ri = rng.randint
    return [ri(0, M) for _ in range(n_trials * n_words)]

# ── Core scoring: run both messages through nr rounds, return mean total XOR HW ──

def score_deltas(deltas, nr, pool, nt, nw):
    """Mean total XOR HW after nr rounds. deltas applied to W[0..nd-1]."""
    nd = len(deltas)
    total = 0
    iv = IV
    for t in range(nt):
        off = t * nw
        a, b, c, d, e, f, g, h = iv
        a2, b2, c2, d2, e2, f2, g2, h2 = iv
        for i in range(nr):
            wi = pool[off + i]
            s1 = (((e >> 6) | (e << 26)) ^ ((e >> 11) | (e << 21)) ^ ((e >> 25) | (e << 7))) & M
            ch = (e & f) ^ ((e ^ M) & g)
            t1 = (h + s1 + ch + K[i] + wi) & M
            s0 = (((a >> 2) | (a << 30)) ^ ((a >> 13) | (a << 19)) ^ ((a >> 22) | (a << 10))) & M
            mj = (a & b) ^ (a & c) ^ (b & c)
            t2 = (s0 + mj) & M
            h = g; g = f; f = e; e = (d + t1) & M
            d = c; c = b; b = a; a = (t1 + t2) & M

            wi2 = wi ^ deltas[i] if i < nd else wi
            s1 = (((e2 >> 6) | (e2 << 26)) ^ ((e2 >> 11) | (e2 << 21)) ^ ((e2 >> 25) | (e2 << 7))) & M
            ch = (e2 & f2) ^ ((e2 ^ M) & g2)
            t1 = (h2 + s1 + ch + K[i] + wi2) & M
            s0 = (((a2 >> 2) | (a2 << 30)) ^ ((a2 >> 13) | (a2 << 19)) ^ ((a2 >> 22) | (a2 << 10))) & M
            mj = (a2 & b2) ^ (a2 & c2) ^ (b2 & c2)
            t2 = (s0 + mj) & M
            h2 = g2; g2 = f2; f2 = e2; e2 = (d2 + t1) & M
            d2 = c2; c2 = b2; b2 = a2; a2 = (t1 + t2) & M

        total += (hw(a^a2) + hw(b^b2) + hw(c^c2) + hw(d^d2) +
                  hw(e^e2) + hw(f^f2) + hw(g^g2) + hw(h^h2))
    return total / nt

def full_eval(deltas, nr, pool, nt, nw):
    """Returns (mean_thw, mean_ehw, p_e_le3, p_e_le5, p_thw_le20)."""
    nd = len(deltas)
    thw_sum = 0; ehw_sum = 0; e_le3 = 0; e_le5 = 0; t_le20 = 0
    iv = IV
    for t in range(nt):
        off = t * nw
        a, b, c, d, e, f, g, h = iv
        a2, b2, c2, d2, e2, f2, g2, h2 = iv
        for i in range(nr):
            wi = pool[off + i]
            s1 = (((e >> 6) | (e << 26)) ^ ((e >> 11) | (e << 21)) ^ ((e >> 25) | (e << 7))) & M
            ch = (e & f) ^ ((e ^ M) & g)
            t1 = (h + s1 + ch + K[i] + wi) & M
            s0 = (((a >> 2) | (a << 30)) ^ ((a >> 13) | (a << 19)) ^ ((a >> 22) | (a << 10))) & M
            mj = (a & b) ^ (a & c) ^ (b & c)
            t2 = (s0 + mj) & M
            h = g; g = f; f = e; e = (d + t1) & M
            d = c; c = b; b = a; a = (t1 + t2) & M

            wi2 = wi ^ deltas[i] if i < nd else wi
            s1 = (((e2 >> 6) | (e2 << 26)) ^ ((e2 >> 11) | (e2 << 21)) ^ ((e2 >> 25) | (e2 << 7))) & M
            ch = (e2 & f2) ^ ((e2 ^ M) & g2)
            t1 = (h2 + s1 + ch + K[i] + wi2) & M
            s0 = (((a2 >> 2) | (a2 << 30)) ^ ((a2 >> 13) | (a2 << 19)) ^ ((a2 >> 22) | (a2 << 10))) & M
            mj = (a2 & b2) ^ (a2 & c2) ^ (b2 & c2)
            t2 = (s0 + mj) & M
            h2 = g2; g2 = f2; f2 = e2; e2 = (d2 + t1) & M
            d2 = c2; c2 = b2; b2 = a2; a2 = (t1 + t2) & M

        thw = hw(a^a2) + hw(b^b2) + hw(c^c2) + hw(d^d2) + hw(e^e2) + hw(f^f2) + hw(g^g2) + hw(h^h2)
        ehw = hw(e ^ e2)
        thw_sum += thw; ehw_sum += ehw
        if ehw <= 3: e_le3 += 1
        if ehw <= 5: e_le5 += 1
        if thw <= 20: t_le20 += 1
    n = float(nt)
    return (thw_sum/n, ehw_sum/n, e_le3/n, e_le5/n, t_le20/n)

# ── Hill climbing ─────────────────────────────────────────────────────────

def hill_climb_single(deltas, idx, nr, n_steps, pool, nt, nw, seed=0):
    rng = random.Random(seed)
    ds = list(deltas)
    best_score = score_deltas(ds, nr, pool, nt, nw)
    best_val = ds[idx]
    ri = rng.randint
    for _ in range(n_steps):
        bit = 1 << ri(0, 31)
        ds[idx] ^= bit
        s = score_deltas(ds, nr, pool, nt, nw)
        if s < best_score:
            best_score = s
            best_val = ds[idx]
        else:
            ds[idx] ^= bit
    ds[idx] = best_val
    return ds, best_score

def hill_climb_pair(deltas, idx1, idx2, nr, n_steps, pool, nt, nw, seed=0):
    rng = random.Random(seed)
    ds = list(deltas)
    best_score = score_deltas(ds, nr, pool, nt, nw)
    ri = rng.randint
    for step in range(n_steps):
        idx = idx1 if step % 2 == 0 else idx2
        bit = 1 << ri(0, 31)
        ds[idx] ^= bit
        s = score_deltas(ds, nr, pool, nt, nw)
        if s < best_score:
            best_score = s
        else:
            ds[idx] ^= bit
    return ds, best_score

# ── Full report ───────────────────────────────────────────────────────────

def full_report(deltas, label, n_trials):
    nw = 16
    pool = make_pool(n_trials, nw, seed=77777)
    print(f"\n{'='*72}")
    print(f"  {label}")
    print(f"  Deltas: {['0x%08x' % d for d in deltas]}")
    print(f"  Verification: {n_trials} trials")
    print(f"{'='*72}")
    max_round = min(len(deltas) + 4, 8)
    print(f"  {'Rnd':>4}  {'MeanTotalHW':>12}  {'MeanE_HW':>9}  {'P(e<=3)':>8}  {'P(e<=5)':>8}  {'P(THW<=20)':>10}")
    print(f"  {'-'*4}  {'-'*12}  {'-'*9}  {'-'*8}  {'-'*8}  {'-'*10}")
    for r in range(1, max_round + 1):
        thw, ehw, pe3, pe5, pt20 = full_eval(deltas, r, pool, n_trials, nw)
        print(f"  {r:>4d}  {thw:>12.2f}  {ehw:>9.2f}  {pe3:>8.4f}  {pe5:>8.4f}  {pt20:>10.4f}")


# ══════════════════════════════════════════════════════════════════════════

def main():
    t0 = time.time()
    best_paths = {}

    DW0 = 0x80000000
    DW1_base = 0x821001c0

    # HC parameters: 50 trials per eval, 16 words per trial
    NT = 50; NW = 16
    pool_hc = make_pool(NT, NW, seed=42)
    pool_check = make_pool(3000, NW, seed=54321)

    # ──────────────────────────────────────────────────────────────────────
    # PHASE 1a: 3-step (fix DW0, DW1, hill-climb DW2)
    # target: minimize state diff after round 3
    # ──────────────────────────────────────────────────────────────────────
    print("=" * 72)
    print("PHASE 1a: 3-step local collision (fix DW0, DW1, hill-climb DW2)")
    print("  DW0=0x80000000, DW1=0x821001c0, optimize DW2")
    print("  1000 seeds, 200 hill-climb steps, 50 trials/eval")
    print("=" * 72)

    best_score_1a = 999.0
    best_ds_1a = None
    for seed in range(1000):
        rng_init = random.Random(seed)
        dw2 = rng_init.randint(0, M)
        ds = [DW0, DW1_base, dw2]
        ds, score = hill_climb_single(ds, 2, 3, 200, pool_hc, NT, NW, seed)
        if score < best_score_1a:
            best_score_1a = score
            best_ds_1a = list(ds)
            print(f"  seed {seed:4d}: DW2=0x{ds[2]:08x}  score={score:.2f}")

    print(f"\n  BEST Phase 1a: DW2=0x{best_ds_1a[2]:08x}  mean_total_hw={best_score_1a:.2f}")
    best_paths['phase1a'] = list(best_ds_1a)
    print(f"  [elapsed: {time.time()-t0:.1f}s]")

    # ──────────────────────────────────────────────────────────────────────
    # PHASE 1b: Joint DW1+DW2 optimization
    # ──────────────────────────────────────────────────────────────────────
    print(f"\n{'='*72}")
    print("PHASE 1b: 3-step (DW0=bit31, jointly optimize DW1+DW2)")
    print("  500 seeds, alternating optimization, 200 steps")
    print("=" * 72)

    best_score_1b = 999.0
    best_ds_1b = None
    for seed in range(500):
        rng_init = random.Random(seed + 10000)
        dw1 = rng_init.randint(0, M)
        dw2 = rng_init.randint(0, M)
        ds = [DW0, dw1, dw2]
        ds, score = hill_climb_pair(ds, 1, 2, 3, 200, pool_hc, NT, NW, seed + 10000)
        if score < best_score_1b:
            best_score_1b = score
            best_ds_1b = list(ds)
            print(f"  seed {seed:4d}: DW1=0x{ds[1]:08x} DW2=0x{ds[2]:08x}  score={score:.2f}")

    print(f"\n  BEST Phase 1b: DW1=0x{best_ds_1b[1]:08x} DW2=0x{best_ds_1b[2]:08x}  score={best_score_1b:.2f}")
    best_paths['phase1b'] = list(best_ds_1b)
    print(f"  [elapsed: {time.time()-t0:.1f}s]")

    # Pick best 3-step
    if best_score_1a <= best_score_1b:
        best_3step = list(best_ds_1a)
        best_3score = best_score_1a
        src = "1a"
    else:
        best_3step = list(best_ds_1b)
        best_3score = best_score_1b
        src = "1b"
    print(f"\n  >>> Best 3-step from Phase {src}: {['0x%08x' % d for d in best_3step]}  hw={best_3score:.2f}")

    # Intermediate check
    print("\n  3-step intermediate check (3000 trials):")
    for r in range(1, 7):
        thw, ehw, pe3, pe5, pt20 = full_eval(best_3step, r, pool_check, 3000, NW)
        print(f"    Round {r}: thw={thw:.2f}  ehw={ehw:.2f}  P(e<=3)={pe3:.4f}  P(THW<=20)={pt20:.4f}")

    # ──────────────────────────────────────────────────────────────────────
    # PHASE 2: 4-step local collision
    # ──────────────────────────────────────────────────────────────────────
    elapsed = time.time() - t0
    print(f"\n{'='*72}")
    print(f"PHASE 2: 4-step local collision (extend best 3-step, hill-climb DW3)")
    print(f"  500 seeds, 100 steps, 50 trials/eval")
    print(f"  [elapsed: {elapsed:.1f}s]")
    print("=" * 72)

    best_score_2 = 999.0
    best_ds_2 = None
    for seed in range(500):
        rng_init = random.Random(seed + 20000)
        dw3 = rng_init.randint(0, M)
        ds = best_3step + [dw3]
        ds, score = hill_climb_single(ds, 3, 4, 100, pool_hc, NT, NW, seed + 20000)
        if score < best_score_2:
            best_score_2 = score
            best_ds_2 = list(ds)
            print(f"  seed {seed:4d}: DW3=0x{ds[3]:08x}  score={score:.2f}")

    print(f"\n  BEST Phase 2: DW3=0x{best_ds_2[3]:08x}  mean_total_hw={best_score_2:.2f}")
    best_paths['phase2'] = list(best_ds_2)

    print("\n  4-step intermediate check (3000 trials):")
    for r in range(1, 7):
        thw, ehw, pe3, pe5, pt20 = full_eval(best_ds_2, r, pool_check, 3000, NW)
        print(f"    Round {r}: thw={thw:.2f}  ehw={ehw:.2f}  P(e<=3)={pe3:.4f}  P(THW<=20)={pt20:.4f}")
    print(f"  [elapsed: {time.time()-t0:.1f}s]")

    # ──────────────────────────────────────────────────────────────────────
    # PHASE 3: Greedy sequential extension
    # ──────────────────────────────────────────────────────────────────────
    print(f"\n{'='*72}")
    print("PHASE 3: Greedy sequential extension (DW0=bit31, DW1..DW7)")
    print("  For each step t: 300 seeds, 100 steps, 50 trials/eval (200 for check)")
    print("=" * 72)

    greedy_ds = [DW0]
    print(f"  Step 0: DW0=0x{DW0:08x} (fixed)")

    for t in range(1, 8):
        elapsed = time.time() - t0
        # Adaptive budget
        if elapsed > 130:
            ns, nst = 80, 40
        elif elapsed > 100:
            ns, nst = 150, 60
        else:
            ns, nst = 300, 100

        target_r = t + 1
        best_score_t = 999.0
        best_dw_t = 0
        for seed in range(ns):
            rng_init = random.Random(seed + 30000 + t * 1000)
            dw = rng_init.randint(0, M)
            ds = greedy_ds + [dw]
            ds, score = hill_climb_single(ds, t, target_r, nst, pool_hc, NT, NW,
                                          seed + 30000 + t * 1000)
            if score < best_score_t:
                best_score_t = score
                best_dw_t = ds[t]

        greedy_ds.append(best_dw_t)
        thw, ehw, pe3, pe5, pt20 = full_eval(greedy_ds, target_r, pool_check, 3000, NW)
        print(f"  Step {t}: DW{t}=0x{best_dw_t:08x}  round={target_r}  "
              f"thw={thw:.2f}  ehw={ehw:.2f}  P(e<=3)={pe3:.4f}  [{time.time()-t0:.1f}s]")

    best_paths['greedy'] = list(greedy_ds)
    print(f"\n  Greedy path: {['0x%08x' % d for d in greedy_ds]}")
    print(f"  [elapsed: {time.time()-t0:.1f}s]")

    # ──────────────────────────────────────────────────────────────────────
    # PHASE 4: Verification and comparison
    # ──────────────────────────────────────────────────────────────────────
    elapsed = time.time() - t0
    remaining = 175 - elapsed
    if remaining > 70:
        verify_n = 50000
    elif remaining > 40:
        verify_n = 20000
    elif remaining > 15:
        verify_n = 10000
    else:
        verify_n = 3000

    print(f"\n{'='*72}")
    print(f"PHASE 4: Verification and comparison ({verify_n} trials)")
    print(f"  [elapsed: {elapsed:.1f}s, remaining budget: {remaining:.1f}s]")
    print("=" * 72)

    full_report([DW0, DW1_base], "BASELINE: 2-step (DW0=bit31, DW1=0x821001c0)", verify_n)

    full_report(best_3step, "BEST 3-step local collision", verify_n)

    full_report(best_ds_2, "BEST 4-step local collision", verify_n)

    full_report(greedy_ds, "GREEDY 8-step sequential extension", verify_n)

    # Phase 1b if different from winner
    if best_paths.get('phase1b') and best_paths['phase1b'] != best_3step:
        vn2 = min(verify_n, 20000)
        full_report(best_paths['phase1b'], "Phase 1b: joint DW1+DW2 optimization", vn2)

    # ── Summary ──
    elapsed = time.time() - t0
    print(f"\n{'='*72}")
    print("SUMMARY")
    print(f"{'='*72}")
    print(f"  2-step baseline:  {['0x%08x' % d for d in [DW0, DW1_base]]}")
    print(f"  Best 3-step:      {['0x%08x' % d for d in best_3step]}  (from Phase {src})")
    print(f"  Best 4-step:      {['0x%08x' % d for d in best_ds_2]}")
    print(f"  Greedy 8-step:    {['0x%08x' % d for d in greedy_ds]}")
    print(f"\n  Total runtime: {elapsed:.1f}s")

if __name__ == '__main__':
    main()
