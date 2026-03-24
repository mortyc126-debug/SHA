#!/usr/bin/env python3
"""
ЗАДАНИЕ 12: Precision Carry Control — δ(d+T1) alignment on round 18

Experiments A-D: v2(S18) analysis, carry propagation depth,
cost-benefit of carry alignment.

Methodology rules:
  Rule 12: null = v2(S18)=0 (random, no alignment)
  Rule 14: if correlation v2(S18)↔HW(δe[25]) > 0.1 → triple check
"""
import random, sys, time, math
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

# ── Primitives ──
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
def sub(a, b): return (a - b) & MASK
def hw(x): return bin(x & MASK).count('1')
def rw(): return random.getrandbits(32)

def schedule(M):
    W = list(M[:16])
    for i in range(16, 64):
        W.append(add(_sig1(W[i-2]), W[i-7], _sig0(W[i-15]), W[i-16]))
    return W

def all_states(W, iv=None):
    """Return states[0..64], each is (a,b,c,d,e,f,g,h)."""
    if iv is None: iv = H0
    s = list(iv)
    states = [tuple(s)]
    for r in range(64):
        T1 = add(s[7], Sig1(s[4]), Ch(s[4], s[5], s[6]), K[r], W[r])
        T2 = add(Sig0(s[0]), Maj(s[0], s[1], s[2]))
        s = [add(T1, T2), s[0], s[1], s[2], add(s[3], T1), s[4], s[5], s[6]]
        states.append(tuple(s))
    return states

def sha256_hash(W):
    """Full SHA-256 hash (single block), returns 8 words."""
    st = all_states(W)
    return [add(H0[i], st[64][i]) for i in range(8)]

# ── Wang chain ──
def wang_chain(Wn):
    """Given Wn[0..15], produce Wf[0..15] with δe[2..17]=0.
    Returns (Wf, states_n_at_16, states_f_at_16) for intermediate access."""
    Wf = list(Wn)
    Wf[0] = Wn[0] ^ 0x8000
    sn = list(H0); sf = list(H0)
    # Round 0
    T1n = add(sn[7], Sig1(sn[4]), Ch(sn[4], sn[5], sn[6]), K[0], Wn[0])
    T2n = add(Sig0(sn[0]), Maj(sn[0], sn[1], sn[2]))
    T1f = add(sf[7], Sig1(sf[4]), Ch(sf[4], sf[5], sf[6]), K[0], Wf[0])
    T2f = add(Sig0(sf[0]), Maj(sf[0], sf[1], sf[2]))
    sn = [add(T1n, T2n), sn[0], sn[1], sn[2], add(sn[3], T1n), sn[4], sn[5], sn[6]]
    sf = [add(T1f, T2f), sf[0], sf[1], sf[2], add(sf[3], T1f), sf[4], sf[5], sf[6]]
    # Rounds 1..15
    for r in range(1, 16):
        dd = sub(sf[3], sn[3])
        dh = sub(sf[7], sn[7])
        dS = sub(Sig1(sf[4]), Sig1(sn[4]))
        dC = sub(Ch(sf[4], sf[5], sf[6]), Ch(sn[4], sn[5], sn[6]))
        dWr = (-(dd + dh + dS + dC)) & MASK
        Wf[r] = add(Wn[r], dWr)
        T1n = add(sn[7], Sig1(sn[4]), Ch(sn[4], sn[5], sn[6]), K[r], Wn[r])
        T2n = add(Sig0(sn[0]), Maj(sn[0], sn[1], sn[2]))
        T1f = add(sf[7], Sig1(sf[4]), Ch(sf[4], sf[5], sf[6]), K[r], Wf[r])
        T2f = add(Sig0(sf[0]), Maj(sf[0], sf[1], sf[2]))
        sn = [add(T1n, T2n), sn[0], sn[1], sn[2], add(sn[3], T1n), sn[4], sn[5], sn[6]]
        sf = [add(T1f, T2f), sf[0], sf[1], sf[2], add(sf[3], T1f), sf[4], sf[5], sf[6]]
    return Wf, tuple(sn), tuple(sf)

def wang_de17(Wn):
    """Compute δe[17] for given Wn. Returns δe[17] as uint32."""
    Wf, sn16, sf16 = wang_chain(Wn)
    Wn_s = schedule(Wn)
    Wf_s = schedule(Wf)
    # Round 16: no message adaptation
    T1n = add(sn16[7], Sig1(sn16[4]), Ch(sn16[4], sn16[5], sn16[6]), K[16], Wn_s[16])
    T1f = add(sf16[7], Sig1(sf16[4]), Ch(sf16[4], sf16[5], sf16[6]), K[16], Wf_s[16])
    en17 = add(sn16[3], T1n)
    ef17 = add(sf16[3], T1f)
    return sub(ef17, en17)

# ── 2-adic valuation ──
def v2(x):
    """2-adic valuation of x. v2(0) = 32."""
    x = x & MASK
    if x == 0: return 32
    v = 0
    while (x & 1) == 0:
        v += 1
        x >>= 1
    return v

# ── Birthday search for δe[17]=0 pairs ──
def find_wang_pairs(target_count, batch_size=500000, verbose=True):
    """Find Wang pairs with δe[2..17]=0 via birthday on W[0].
    For each random Wn[1..15], sweep W[0] through random values."""
    pairs = []
    attempts = 0
    t0 = time.time()

    while len(pairs) < target_count:
        # Random background
        Wn_base = [0] + [rw() for _ in range(15)]
        # Try random W[0] values
        for _ in range(batch_size):
            Wn_base[0] = rw()
            de17 = wang_de17(Wn_base)
            attempts += 1
            if de17 == 0:
                Wn = list(Wn_base)
                Wf, _, _ = wang_chain(Wn)
                pairs.append((Wn, Wf))
                if verbose and len(pairs) % 10 == 0:
                    elapsed = time.time() - t0
                    rate = attempts / elapsed if elapsed > 0 else 0
                    print(f"  Found {len(pairs)}/{target_count} pairs "
                          f"({attempts:.2e} attempts, {rate:.0f}/s, {elapsed:.1f}s)")
                break  # new background for next search

    elapsed = time.time() - t0
    if verbose:
        print(f"  Birthday search complete: {len(pairs)} pairs in {elapsed:.1f}s "
              f"({attempts:.2e} attempts)")
    return pairs

# ── Collect per-pair data ──
def analyze_pair(Wn, Wf):
    """For a Wang pair, compute all metrics needed for experiments A-D."""
    Wn_s = schedule(Wn)
    Wf_s = schedule(Wf)
    sn_states = all_states(Wn_s)
    sf_states = all_states(Wf_s)

    # δe[r] for r=0..64
    de = [sub(sf_states[r][4], sn_states[r][4]) for r in range(65)]
    # δa[r] for r=0..64
    da = [sub(sf_states[r][0], sn_states[r][0]) for r in range(65)]

    # δa[14] = a_f[14] - a_n[14]
    da14 = da[14]
    # δW[17] from expanded schedule
    dW17 = sub(Wf_s[17], Wn_s[17])
    # S18 = δa[14] + δW[17] mod 2^32
    S18 = add(da14, dW17)
    v2_S18 = v2(S18)

    # Hash
    hn = sha256_hash(Wn_s)
    hf = sha256_hash(Wf_s)
    dhash = 0
    hw_hash = 0
    for i in range(8):
        d = hn[i] ^ hf[i]
        hw_hash += hw(d)

    # HW of δe for rounds 18..30
    hw_de = {}
    for r in range(18, 31):
        hw_de[r] = hw(de[r])
    hw_de[64] = hw(de[64])

    # Low-bit check: HW of lower k bits of δe[18] where k = v2_S18
    de18 = de[18]
    low_k_hw = 0
    if v2_S18 > 0 and v2_S18 <= 32:
        low_mask = (1 << min(v2_S18, 32)) - 1
        low_k_hw = hw(de18 & low_mask)

    return {
        'da14': da14, 'dW17': dW17, 'S18': S18, 'v2_S18': v2_S18,
        'de18': de18, 'hw_de18': hw(de18),
        'hw_de': hw_de,  # r -> hw(δe[r])
        'hw_hash': hw_hash,
        'low_k_hw': low_k_hw,
        'de': de,
    }


# ═══════════════════════════════════════════════════════════════
# EXPERIMENT A: v2(S18) vs HW(δe[18]) and HW(δhash)
# ═══════════════════════════════════════════════════════════════
def experiment_A(data_list):
    print("=" * 80)
    print("EXPERIMENT A: v2(S18) vs HW(δe[18]) — Carry Alignment Effect")
    print("=" * 80)
    print(f"Total pairs: {len(data_list)}")

    # Group by v2(S18)
    groups = defaultdict(list)
    for d in data_list:
        groups[d['v2_S18']].append(d)

    print(f"\n{'v2(S18)':>8} | {'count':>6} | {'mean HW(δe18)':>14} | {'std':>6} "
          f"| {'mean HW(δe19)':>14} | {'mean HW(δhash)':>15} | {'mean low_k_hw':>14}")
    print("-" * 100)

    for v in sorted(groups.keys()):
        g = groups[v]
        n = len(g)
        hw18 = [d['hw_de18'] for d in g]
        hw19 = [d['hw_de'].get(19, 0) for d in g]
        hwh = [d['hw_hash'] for d in g]
        lkh = [d['low_k_hw'] for d in g]

        mean18 = sum(hw18) / n
        std18 = (sum((x - mean18)**2 for x in hw18) / max(n-1, 1)) ** 0.5
        mean19 = sum(hw19) / n
        meanh = sum(hwh) / n
        meanlk = sum(lkh) / n

        print(f"{v:>8} | {n:>6} | {mean18:>14.2f} | {std18:>6.2f} "
              f"| {mean19:>14.2f} | {meanh:>15.2f} | {meanlk:>14.2f}")

    # Theoretical check: HW(δe[18][0..k-1]) should be ~0 when v2(S18)≥k
    print(f"\n--- Theoretical check: HW of lower k bits of δe[18] ---")
    print("If v2(S18) = k, then δe[18] bits 0..k-1 should be carry-free.")
    print("low_k_hw = HW(δe[18] & ((1<<k)-1)) should be ≈ 0 if carry alignment works.\n")

    for v in sorted(groups.keys()):
        g = groups[v]
        if len(g) < 3: continue
        lkh = [d['low_k_hw'] for d in g]
        mean_lkh = sum(lkh) / len(g)
        zero_count = sum(1 for x in lkh if x == 0)
        print(f"  v2={v:>2}: mean low_k_hw = {mean_lkh:.3f}, "
              f"zero_fraction = {zero_count}/{len(g)} = {zero_count/len(g):.3f}")

    # Null comparison (Rule 12)
    null_group = groups.get(0, [])
    if null_group:
        null_hw18 = sum(d['hw_de18'] for d in null_group) / len(null_group)
        null_hwh = sum(d['hw_hash'] for d in null_group) / len(null_group)
        print(f"\n  Null (v2=0): mean HW(δe18) = {null_hw18:.2f}, mean HW(δhash) = {null_hwh:.2f}")
        high_v2 = [d for d in data_list if d['v2_S18'] >= 4]
        if high_v2:
            hv_hw18 = sum(d['hw_de18'] for d in high_v2) / len(high_v2)
            hv_hwh = sum(d['hw_hash'] for d in high_v2) / len(high_v2)
            print(f"  v2≥4 ({len(high_v2)} pairs): mean HW(δe18) = {hv_hw18:.2f}, "
                  f"mean HW(δhash) = {hv_hwh:.2f}")
            print(f"  Δ HW(δe18) = {hv_hw18 - null_hw18:.2f}, Δ HW(δhash) = {hv_hwh - null_hwh:.2f}")


# ═══════════════════════════════════════════════════════════════
# EXPERIMENT C: Carry propagation depth — correlation analysis
# ═══════════════════════════════════════════════════════════════
def experiment_C(data_list):
    print("\n" + "=" * 80)
    print("EXPERIMENT C: Carry Propagation Depth — v2(S18) correlation with δe[r]")
    print("=" * 80)
    print(f"Total pairs: {len(data_list)}")

    # Pearson correlation between v2(S18) and HW(δe[r])
    v2_vals = [d['v2_S18'] for d in data_list]
    mean_v2 = sum(v2_vals) / len(v2_vals)
    var_v2 = sum((x - mean_v2)**2 for x in v2_vals) / len(v2_vals)
    std_v2 = var_v2 ** 0.5

    print(f"\nPearson correlation: corr(v2(S18), HW(δe[r]))")
    print(f"  Negative correlation = higher v2 → lower HW (good!)")
    print(f"  Rule 14: |corr| > 0.1 at r≥25 → triple check\n")
    print(f"{'round':>6} | {'corr':>8} | {'mean HW':>8} | {'std HW':>7} | {'significance':>12}")
    print("-" * 55)

    significant_rounds = []
    for r in range(18, 31):
        hw_vals = [d['hw_de'].get(r, hw(d['de'][r])) for d in data_list]
        mean_hw = sum(hw_vals) / len(hw_vals)
        var_hw = sum((x - mean_hw)**2 for x in hw_vals) / len(hw_vals)
        std_hw = var_hw ** 0.5

        if std_v2 > 0 and std_hw > 0:
            cov = sum((v2_vals[i] - mean_v2) * (hw_vals[i] - mean_hw)
                      for i in range(len(data_list))) / len(data_list)
            corr = cov / (std_v2 * std_hw)
        else:
            corr = 0.0

        # Significance: |corr| * sqrt(n) as t-statistic approximation
        n = len(data_list)
        t_stat = abs(corr) * (n - 2)**0.5 / max((1 - corr**2)**0.5, 1e-10)
        sig = "***" if t_stat > 3.29 else "**" if t_stat > 2.58 else "*" if t_stat > 1.96 else ""

        print(f"  r={r:>2} | {corr:>+8.4f} | {mean_hw:>8.2f} | {std_hw:>7.2f} | {sig:>12}")

        if abs(corr) > 0.1 and r >= 25:
            significant_rounds.append((r, corr))

    # Rule 14: triple check for significant late-round correlations
    if significant_rounds:
        print(f"\n  *** Rule 14 triggered: significant correlation at late rounds ***")
        for r, corr in significant_rounds:
            print(f"    r={r}: corr={corr:+.4f} — needs triple verification")

    # Additional: correlation with HW(δhash)
    hwh_vals = [d['hw_hash'] for d in data_list]
    mean_hwh = sum(hwh_vals) / len(hwh_vals)
    var_hwh = sum((x - mean_hwh)**2 for x in hwh_vals) / len(hwh_vals)
    std_hwh = var_hwh ** 0.5
    if std_v2 > 0 and std_hwh > 0:
        cov = sum((v2_vals[i] - mean_v2) * (hwh_vals[i] - mean_hwh)
                  for i in range(len(data_list))) / len(data_list)
        corr_hash = cov / (std_v2 * std_hwh)
    else:
        corr_hash = 0.0
    print(f"\n  corr(v2(S18), HW(δhash)) = {corr_hash:+.4f}")

    # Grouped means for visualization
    print(f"\n--- Grouped: mean HW(δe[r]) by v2(S18) bucket ---")
    buckets = [(0, 0, "v2=0"), (1, 1, "v2=1"), (2, 3, "v2=2-3"),
               (4, 7, "v2=4-7"), (8, 31, "v2=8-31"), (32, 32, "v2=32")]

    header = f"{'bucket':>10}"
    for r in range(18, 26):
        header += f" | {'r='+str(r):>7}"
    header += f" | {'hash':>7}"
    print(header)
    print("-" * (12 + 10 * 9))

    for lo, hi, label in buckets:
        grp = [d for d in data_list if lo <= d['v2_S18'] <= hi]
        if not grp: continue
        row = f"{label:>10}"
        for r in range(18, 26):
            m = sum(hw(d['de'][r]) for d in grp) / len(grp)
            row += f" | {m:>7.2f}"
        mh = sum(d['hw_hash'] for d in grp) / len(grp)
        row += f" | {mh:>7.2f}"
        row += f"  (n={len(grp)})"
        print(row)


# ═══════════════════════════════════════════════════════════════
# EXPERIMENT D: Cost-benefit table for k=1,2,4,8,12,16
# ═══════════════════════════════════════════════════════════════
def experiment_D(data_list):
    print("\n" + "=" * 80)
    print("EXPERIMENT D: Optimal k — Cost vs Propagation Benefit")
    print("=" * 80)
    print(f"Total pairs: {len(data_list)}")

    # Check if S18 distribution is uniform
    print("\n--- Distribution check: P(v2(S18) >= k) vs 2^{-k} ---")
    n = len(data_list)
    for k in [1, 2, 4, 8, 12, 16]:
        count = sum(1 for d in data_list if d['v2_S18'] >= k)
        p_obs = count / n
        p_exp = 2**(-k)
        ratio = p_obs / p_exp if p_exp > 0 else float('inf')
        print(f"  k={k:>2}: observed P={p_obs:.6f} ({count}/{n}), "
              f"expected 2^{{-{k}}}={p_exp:.6f}, ratio={ratio:.3f}")

    # Cost-benefit table
    print(f"\n--- Cost-benefit table ---")
    print(f"  Cost = birthday pairs needed ≈ 2^(32+k) if P(v2≥k) = 2^{-k}")
    print(f"  Benefit = reduction in HW(δe[18]), HW(δe[22]), HW(δhash)\n")

    # Null group
    null = [d for d in data_list if d['v2_S18'] == 0]
    if not null:
        null = data_list  # fallback
    null_hw18 = sum(d['hw_de18'] for d in null) / len(null)
    null_hw22 = sum(hw(d['de'][22]) for d in null) / len(null)
    null_hwh = sum(d['hw_hash'] for d in null) / len(null)

    print(f"{'k':>4} | {'cost':>12} | {'count':>6} | {'HW(δe18)':>10} | {'HW(δe22)':>10} "
          f"| {'HW(δhash)':>10} | {'Δe18':>7} | {'Δe22':>7} | {'Δhash':>7}")
    print("-" * 100)
    print(f"{'null':>4} | {'baseline':>12} | {len(null):>6} | {null_hw18:>10.2f} | {null_hw22:>10.2f} "
          f"| {null_hwh:>10.2f} | {'—':>7} | {'—':>7} | {'—':>7}")

    for k in [1, 2, 4, 8, 12, 16]:
        grp = [d for d in data_list if d['v2_S18'] >= k]
        if not grp:
            print(f"{k:>4} | {'2^'+str(32+k):>12} | {0:>6} | {'n/a':>10} | {'n/a':>10} "
                  f"| {'n/a':>10} | {'n/a':>7} | {'n/a':>7} | {'n/a':>7}")
            continue
        m18 = sum(d['hw_de18'] for d in grp) / len(grp)
        m22 = sum(hw(d['de'][22]) for d in grp) / len(grp)
        mh = sum(d['hw_hash'] for d in grp) / len(grp)
        d18 = m18 - null_hw18
        d22 = m22 - null_hw22
        dh = mh - null_hwh
        print(f"{k:>4} | {'2^'+str(32+k):>12} | {len(grp):>6} | {m18:>10.2f} | {m22:>10.2f} "
              f"| {mh:>10.2f} | {d18:>+7.2f} | {d22:>+7.2f} | {dh:>+7.2f}")


# ═══════════════════════════════════════════════════════════════
# EXPERIMENT B: Two-level birthday — δe[17]=0 AND v2(S18)≥8
# ═══════════════════════════════════════════════════════════════
def experiment_B(data_list):
    print("\n" + "=" * 80)
    print("EXPERIMENT B: Two-Level Birthday — δe[17]=0 AND v2(S18) ≥ k")
    print("=" * 80)
    print(f"Total pairs with δe[17]=0: {len(data_list)}")

    for k_thresh in [4, 8, 12, 16]:
        filtered = [d for d in data_list if d['v2_S18'] >= k_thresh]
        expected = len(data_list) / (2**k_thresh)
        print(f"\n--- v2(S18) ≥ {k_thresh}: found {len(filtered)}, expected ≈ {expected:.1f} ---")

        if not filtered:
            print("  No pairs found.")
            continue

        hw18 = [d['hw_de18'] for d in filtered]
        hw19 = [d['hw_de'].get(19, 0) for d in filtered]
        hwh = [d['hw_hash'] for d in filtered]

        all_hw18 = [d['hw_de18'] for d in data_list]
        all_hw19 = [d['hw_de'].get(19, 0) for d in data_list]
        all_hwh = [d['hw_hash'] for d in data_list]

        print(f"  HW(δe[18]): filtered mean = {sum(hw18)/len(hw18):.2f}, "
              f"baseline mean = {sum(all_hw18)/len(all_hw18):.2f}")
        print(f"  HW(δe[19]): filtered mean = {sum(hw19)/len(hw19):.2f}, "
              f"baseline mean = {sum(all_hw19)/len(all_hw19):.2f}")
        print(f"  HW(δhash):  filtered mean = {sum(hwh)/len(hwh):.2f}, "
              f"baseline mean = {sum(all_hwh)/len(all_hwh):.2f}")

        # Cascade viability assessment
        if len(filtered) >= 5:
            improvement_e18 = sum(all_hw18)/len(all_hw18) - sum(hw18)/len(hw18)
            improvement_hash = sum(all_hwh)/len(all_hwh) - sum(hwh)/len(hwh)
            print(f"  Improvement: ΔHW(δe18) = {improvement_e18:+.2f} bits, "
                  f"ΔHW(δhash) = {improvement_hash:+.2f} bits")
            if improvement_e18 > 2:
                print(f"  → Cascade potential: YES (significant δe18 reduction)")
            else:
                print(f"  → Cascade potential: marginal or none")


# ═══════════════════════════════════════════════════════════════
# PAIR LOADING from C generator output
# ═══════════════════════════════════════════════════════════════
def parse_c_pairs(filename):
    """Parse PAIR/Wn/Wf output from C generators."""
    pairs = []
    with open(filename) as f:
        lines = [l.strip() for l in f if l.strip() and not l.startswith('#')]
    i = 0
    seen = set()
    while i < len(lines):
        if lines[i] == 'PAIR':
            if i+2 < len(lines):
                wn_hex = lines[i+1].replace('Wn =', '').strip().split()
                wf_hex = lines[i+2].replace('Wf =', '').strip().split()
                if len(wn_hex) == 16 and len(wf_hex) == 16:
                    Wn = [int(x, 16) for x in wn_hex]
                    Wf = [int(x, 16) for x in wf_hex]
                    key = tuple(Wn)
                    if key not in seen:
                        seen.add(key)
                        pairs.append((Wn, Wf))
                i += 3
            else:
                i += 1
        else:
            i += 1
    return pairs


# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════
def main():
    import os, glob as glob_mod

    print("=" * 80)
    print(f"ЗАДАНИЕ 12: Precision Carry Control — δ(d+T1) alignment on round 18")
    print("=" * 80)

    # Step 1: Load pairs from C output OR generate
    pairs = []

    # Try loading from file
    pair_file = 'task12_pairs.txt'
    if os.path.exists(pair_file) and os.path.getsize(pair_file) > 0:
        print(f"\n[Step 1] Loading pairs from {pair_file}...")
        pairs = parse_c_pairs(pair_file)
        print(f"  Loaded {len(pairs)} unique pairs")

    if not pairs:
        N_PAIRS = 10 if len(sys.argv) <= 1 else int(sys.argv[1])
        print(f"\n[Step 1] No pre-generated pairs. Generating {N_PAIRS} via birthday search...")
        t0 = time.time()
        pairs = find_wang_pairs(N_PAIRS, batch_size=500000, verbose=True)
        print(f"  Search time: {time.time()-t0:.1f}s")

    print(f"  Total pairs for analysis: {len(pairs)}")

    # Step 2: Analyze all pairs
    print(f"\n[Step 2] Analyzing {len(pairs)} pairs...")
    t0 = time.time()
    data_list = []
    for i, (Wn, Wf) in enumerate(pairs):
        d = analyze_pair(Wn, Wf)
        data_list.append(d)
        if (i+1) % 1000 == 0:
            print(f"  Analyzed {i+1}/{len(pairs)}...")
    t_analyze = time.time() - t0
    print(f"  Analysis time: {t_analyze:.1f}s")

    # Verify: all δe[17] should be 0
    de17_nonzero = sum(1 for d in data_list if d['de'][17] != 0)
    print(f"  Verification: δe[17]≠0 count = {de17_nonzero} (should be 0)")

    # Important theoretical note
    print(f"\n  NOTE: δe[18] = S18 = δa[14] + δW[17] exactly (when δe[17]=0)")
    print(f"        So HW(δe[18]) = HW(S18) by definition.")
    print(f"        The INTERESTING question is propagation to rounds 19+.")

    # v2(S18) distribution summary
    v2_dist = defaultdict(int)
    for d in data_list:
        v2_dist[d['v2_S18']] += 1
    print(f"\n  v2(S18) distribution:")
    for v in sorted(v2_dist.keys()):
        print(f"    v2={v}: {v2_dist[v]} ({100*v2_dist[v]/len(data_list):.2f}%)")

    # Step 3: Run experiments (priority A > C > D > B)
    print()
    experiment_A(data_list)
    experiment_C(data_list)
    experiment_D(data_list)
    experiment_B(data_list)

    # Summary / Conclusions
    print("\n" + "=" * 80)
    print("THEORETICAL ANALYSIS (supplements empirical data)")
    print("=" * 80)
    print("""
When δe[2..17]=0 (Wang chain + birthday on round 17):

1. δe[18] = S18 = δa[14] + δW[17] EXACTLY (mod 2^32)
   - This is because δh[17]=δe[14]=0, δSig1(e[17])=0, δCh(e[17],f[17],g[17])=0
   - So δT1[17] = δW[17] only, and δe[18] = δd[17] + δT1[17] = δa[14] + δW[17]

2. v2(S18) = k means low k bits of δe[18] are 0
   - Expected HW(δe[18]) when v2(S18)=k: ≈ (32-k)/2 = 16 - k/2
   - This is a TRIVIAL consequence, not a SHA-256 property

3. The KEY question: does v2(S18) affect rounds 19+?
   - δe[19] depends on δSig1(e[18]) and δCh(e[18],f[18],g[18])
   - These are NONLINEAR in δe[18], so empirical measurement is essential
   - If v2(S18)=k, the low k bits of e_f[18] and e_n[18] AGREE on carries
   - This means Sig1(e_f[18]) and Sig1(e_n[18]) partially agree (via ROTR6,11,25)
   - The partial agreement propagates to δe[19] through the carry chain

4. Cascade feasibility:
   - If corr(v2(S18), HW(δe[r])) > 0 for r >> 18: carry alignment propagates
   - Cost of v2(S18)≥k: birthday 2^{32+k} (if S18 uniform mod 2^k)
   - Each level "freezes" k carry bits → need ceil(32/k) levels for full collision
   - Naive estimate: 6 levels × 2^8 = 2^48 (for k=8)
   - BUT this requires carry alignment to PROPAGATE, which SHA-256 resists
""")

    print("=" * 80)
    print("DONE")
    print("=" * 80)

if __name__ == '__main__':
    main()
