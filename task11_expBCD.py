#!/usr/bin/env python3
"""
ЗАДАНИЕ 11, Experiments B, C, D:
  B: δ-persistent bits analysis
  C: Alternating cascade vs Wang carry agreement
  D: Round-by-round healing analysis
"""
import json, os, sys, re, struct, math, random
from collections import Counter

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
    for i in range(16, 64): W.append(add(_sig1(W[i-2]), W[i-7], _sig0(W[i-15]), W[i-16]))
    return W

def all_states(W, iv=None):
    if iv is None: iv = H0
    s = list(iv)
    states = [tuple(s)]
    for r in range(64):
        T1 = add(s[7], Sig1(s[4]), Ch(s[4], s[5], s[6]), K[r], W[r])
        T2 = add(Sig0(s[0]), Maj(s[0], s[1], s[2]))
        s = [add(T1, T2), s[0], s[1], s[2], add(s[3], T1), s[4], s[5], s[6]]
        states.append(tuple(s))
    return states

def wang_chain(Wn, iv=None):
    if iv is None: iv = H0
    Wf = list(Wn); Wf[0] = Wn[0] ^ 0x8000
    sn = list(iv); sf = list(iv)
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
        dWr = (-(dd + dh + dS + dC) & MASK) & MASK
        Wf[r] = add(Wn[r], dWr)
        T1n = add(sn[7], Sig1(sn[4]), Ch(sn[4], sn[5], sn[6]), K[r], Wn[r])
        T2n = add(Sig0(sn[0]), Maj(sn[0], sn[1], sn[2]))
        T1f = add(sf[7], Sig1(sf[4]), Ch(sf[4], sf[5], sf[6]), K[r], Wf[r])
        T2f = add(Sig0(sf[0]), Maj(sf[0], sf[1], sf[2]))
        sn = [add(T1n, T2n), sn[0], sn[1], sn[2], add(sn[3], T1n), sn[4], sn[5], sn[6]]
        sf = [add(T1f, T2f), sf[0], sf[1], sf[2], add(sf[3], T1f), sf[4], sf[5], sf[6]]
    return Wf

def compute_carries(x, y):
    carries = 0; c = 0
    for i in range(32):
        xi = (x >> i) & 1; yi = (y >> i) & 1
        c = (xi & yi) | (xi & c) | (yi & c)
        if c: carries |= (1 << i)
    return carries

def carry_agreement_pair(Wn, Wf):
    """Return carry_agreement(r) for r=0..63."""
    Wn_s = schedule(Wn); Wf_s = schedule(Wf)
    sn = list(H0); sf = list(H0)
    ca_list = []
    for r in range(64):
        # compute round intermediates for both
        s1n = Sig1(sn[4]); s1f = Sig1(sf[4])
        chn = Ch(sn[4], sn[5], sn[6]); chf = Ch(sf[4], sf[5], sf[6])
        s0n = Sig0(sn[0]); s0f = Sig0(sf[0])
        mjn = Maj(sn[0], sn[1], sn[2]); mjf = Maj(sf[0], sf[1], sf[2])

        sum1n = add(sn[7], s1n); sum1f = add(sf[7], s1f)
        sum2n = add(sum1n, chn); sum2f = add(sum1f, chf)
        sum3n = add(sum2n, K[r]); sum3f = add(sum2f, K[r])
        sum4n = add(sum3n, Wn_s[r]); sum4f = add(sum3f, Wf_s[r])
        sum5n = add(s0n, mjn); sum5f = add(s0f, mjf)

        cn = [compute_carries(sn[7], s1n), compute_carries(sum1n, chn),
              compute_carries(sum2n, K[r]), compute_carries(sum3n, Wn_s[r]),
              compute_carries(s0n, mjn)]
        cf = [compute_carries(sf[7], s1f), compute_carries(sum1f, chf),
              compute_carries(sum2f, K[r]), compute_carries(sum3f, Wf_s[r]),
              compute_carries(s0f, mjf)]

        agree = sum(32 - hw(cn[i] ^ cf[i]) for i in range(5))
        ca_list.append(agree / 160.0)

        T1n = sum4n; T1f = sum4f; T2n = sum5n; T2f = sum5f
        sn = [add(T1n, T2n), sn[0], sn[1], sn[2], add(sn[3], T1n), sn[4], sn[5], sn[6]]
        sf = [add(T1f, T2f), sf[0], sf[1], sf[2], add(sf[3], T1f), sf[4], sf[5], sf[6]]

    return ca_list

def parse_pairs(filename):
    pairs = []
    with open(filename) as f:
        content = f.read()
    wn_lines = re.findall(r'Wn = ([0-9a-f ]+)', content)
    wf_lines = re.findall(r'Wf = ([0-9a-f ]+)', content)
    for wn_str, wf_str in zip(wn_lines, wf_lines):
        Wn = [int(x, 16) for x in wn_str.strip().split()]
        Wf = [int(x, 16) for x in wf_str.strip().split()]
        if len(Wn) == 16 and len(Wf) == 16:
            pairs.append((Wn, Wf))
    return pairs


# ═══════════════════════════════════════════════════════════════
# EXPERIMENT B: δ-persistent bits
# ═══════════════════════════════════════════════════════════════
def experiment_B(pairs):
    print("=" * 80)
    print("EXPERIMENT B: δ-persistent bits analysis")
    print("=" * 80)

    n_pairs = len(pairs)
    # persistence[pair][word][bit] = count of rounds 17..64 where δ=0
    R_START, R_END = 17, 64  # states[17]..states[64] = 48 rounds
    N_ROUNDS = R_END - R_START + 1  # 48

    all_persistence = []  # list of 8x32 arrays

    for pi, (Wn, Wf) in enumerate(pairs):
        Wn_s = schedule(Wn); Wf_s = schedule(Wf)
        sn_states = all_states(Wn_s)
        sf_states = all_states(Wf_s)

        pers = [[0]*32 for _ in range(8)]
        for r in range(R_START, R_END + 1):
            for w in range(8):
                d = sn_states[r][w] ^ sf_states[r][w]
                for b in range(32):
                    if ((d >> b) & 1) == 0:
                        pers[w][b] += 1
        all_persistence.append(pers)

    # 1. Per-pair: bits with persistence > 0.5 and > 0.7
    print(f"\n--- Per-pair persistence counts ---")
    print(f"(out of 256 state bits, {N_ROUNDS} rounds checked)")
    for pi in range(n_pairs):
        count_05 = sum(1 for w in range(8) for b in range(32)
                       if all_persistence[pi][w][b] / N_ROUNDS > 0.5)
        count_06 = sum(1 for w in range(8) for b in range(32)
                       if all_persistence[pi][w][b] / N_ROUNDS > 0.6)
        count_07 = sum(1 for w in range(8) for b in range(32)
                       if all_persistence[pi][w][b] / N_ROUNDS > 0.7)
        print(f"  Pair {pi}: pers>0.5: {count_05}, pers>0.6: {count_06}, pers>0.7: {count_07}")

    # 2. Bits with persistence > 0.7 in any pair
    print(f"\n--- Bits with persistence > 0.7 (any pair) ---")
    reg_names = ['a','b','c','d','e','f','g','h']
    high_pers_bits = []
    for pi in range(n_pairs):
        for w in range(8):
            for b in range(32):
                p = all_persistence[pi][w][b] / N_ROUNDS
                if p > 0.7:
                    high_pers_bits.append((pi, w, b, p))
    if high_pers_bits:
        print(f"  Found {len(high_pers_bits)} (pair,word,bit) with persistence > 0.7:")
        for pi, w, b, p in sorted(high_pers_bits, key=lambda x: -x[3])[:30]:
            print(f"    Pair {pi}: reg={reg_names[w]} bit={b:2d} persistence={p:.3f} ({all_persistence[pi][w][b]}/{N_ROUNDS})")
    else:
        print("  None found.")

    # 3. Cross-persistence: for each bit, how many pairs have persistence > 0.6?
    print(f"\n--- Cross-persistence (bits persistent across multiple pairs) ---")
    cross = [[0]*32 for _ in range(8)]
    for pi in range(n_pairs):
        for w in range(8):
            for b in range(32):
                if all_persistence[pi][w][b] / N_ROUNDS > 0.6:
                    cross[w][b] += 1

    # Distribution
    cross_dist = Counter()
    for w in range(8):
        for b in range(32):
            cross_dist[cross[w][b]] += 1
    print(f"  Cross-persistence distribution (how many of 7 pairs have pers>0.6):")
    for k in sorted(cross_dist.keys()):
        print(f"    cross={k}: {cross_dist[k]} bits")

    # Bits with cross >= 5 (persistent in majority of pairs)
    structural_bits = []
    for w in range(8):
        for b in range(32):
            if cross[w][b] >= 5:
                structural_bits.append((w, b, cross[w][b]))
    if structural_bits:
        print(f"\n  Bits with cross_persistence >= 5:")
        for w, b, c in sorted(structural_bits, key=lambda x: (-x[2], x[0], x[1])):
            print(f"    reg={reg_names[w]} bit={b:2d}: persistent in {c}/7 pairs")
    else:
        print(f"\n  No bits with cross_persistence >= 5.")

    # Bits with cross = 7 (ALL pairs)
    all7 = [(w, b) for w in range(8) for b in range(32) if cross[w][b] == 7]
    if all7:
        print(f"\n  *** STRUCTURAL: {len(all7)} bits persistent in ALL 7 pairs: ***")
        for w, b in all7:
            print(f"    reg={reg_names[w]} bit={b:2d}")
    else:
        print(f"\n  No bits persistent in all 7 pairs.")

    # 4. Register and position breakdown
    print(f"\n--- Breakdown of bits with cross_persistence >= 3 ---")
    # By register
    reg_counts = [0]*8
    pos_counts = [0]*32
    for w in range(8):
        for b in range(32):
            if cross[w][b] >= 3:
                reg_counts[w] += 1
                pos_counts[b] += 1
    print("  By register:")
    for w in range(8):
        print(f"    {reg_names[w]}: {reg_counts[w]}/32 bits")
    print("  By bit position (high bits = 31..24, mid=23..8, low=7..0):")
    low = sum(pos_counts[b] for b in range(8))
    mid = sum(pos_counts[b] for b in range(8, 24))
    high = sum(pos_counts[b] for b in range(24, 32))
    print(f"    Low (0-7): {low} across all regs")
    print(f"    Mid (8-23): {mid}")
    print(f"    High (24-31): {high}")

    # 5. Random baseline: how many bits have persistence > 0.6 for random pairs?
    print(f"\n--- Random baseline (100 random Wang pairs) ---")
    random_counts = []
    for _ in range(100):
        Wn = [rw() for _ in range(16)]
        Wf = wang_chain(Wn)
        sn = all_states(schedule(Wn))
        sf = all_states(schedule(Wf))
        count = 0
        for w in range(8):
            for b in range(32):
                zero_count = 0
                for r in range(R_START, R_END + 1):
                    if ((sn[r][w] ^ sf[r][w]) >> b) & 1 == 0:
                        zero_count += 1
                if zero_count / N_ROUNDS > 0.6:
                    count += 1
        random_counts.append(count)
    avg_r = sum(random_counts) / len(random_counts)
    print(f"  Mean bits with pers>0.6 for random Wang: {avg_r:.1f}")
    print(f"  Our 7 pairs: {', '.join(str(sum(1 for w in range(8) for b in range(32) if all_persistence[pi][w][b]/N_ROUNDS > 0.6)) for pi in range(n_pairs))}")


# ═══════════════════════════════════════════════════════════════
# EXPERIMENT C: Alternating cascade vs Wang
# ═══════════════════════════════════════════════════════════════
def alternating_chain(Wn, iv=None):
    """Build alternating chain: even rounds δe=0, odd rounds δa=0.
    Returns Wf[0..15]."""
    if iv is None: iv = H0
    Wf = list(Wn); Wf[0] = Wn[0] ^ 0x8000
    sn = list(iv); sf = list(iv)

    for r in range(16):
        if r == 0:
            # Standard: just apply δW[0]
            T1n = add(sn[7], Sig1(sn[4]), Ch(sn[4], sn[5], sn[6]), K[0], Wn[0])
            T2n = add(Sig0(sn[0]), Maj(sn[0], sn[1], sn[2]))
            T1f = add(sf[7], Sig1(sf[4]), Ch(sf[4], sf[5], sf[6]), K[0], Wf[0])
            T2f = add(Sig0(sf[0]), Maj(sf[0], sf[1], sf[2]))
            sn = [add(T1n, T2n), sn[0], sn[1], sn[2], add(sn[3], T1n), sn[4], sn[5], sn[6]]
            sf = [add(T1f, T2f), sf[0], sf[1], sf[2], add(sf[3], T1f), sf[4], sf[5], sf[6]]
            continue

        if r % 2 == 0:
            # Even rounds: target δe[r+1] = 0
            # new_e = d + T1, so δ(new_e) = δd + δT1 = 0
            # δT1 = δh + δ(Sig1(e)) + δ(Ch(e,f,g)) + δW[r]
            # => δW[r] = -(δd + δh + δSig1 + δCh)
            dd = (sf[3] - sn[3]) & MASK
            dh = (sf[7] - sn[7]) & MASK
            dS = (Sig1(sf[4]) - Sig1(sn[4])) & MASK
            dC = (Ch(sf[4], sf[5], sf[6]) - Ch(sn[4], sn[5], sn[6])) & MASK
            dWr = (-(dd + dh + dS + dC)) & MASK
        else:
            # Odd rounds: target δa[r+1] = 0
            # new_a = T1 + T2, so δ(new_a) = δT1 + δT2 = 0
            # δT1 = δh + δSig1 + δCh + δW[r]
            # δT2 = δSig0 + δMaj
            # => δW[r] = -(δh + δSig1 + δCh + δSig0 + δMaj)
            dh = (sf[7] - sn[7]) & MASK
            dS1 = (Sig1(sf[4]) - Sig1(sn[4])) & MASK
            dC = (Ch(sf[4], sf[5], sf[6]) - Ch(sn[4], sn[5], sn[6])) & MASK
            dS0 = (Sig0(sf[0]) - Sig0(sn[0])) & MASK
            dM = (Maj(sf[0], sf[1], sf[2]) - Maj(sn[0], sn[1], sn[2])) & MASK
            dWr = (-(dh + dS1 + dC + dS0 + dM)) & MASK

        Wf[r] = add(Wn[r], dWr)
        T1n = add(sn[7], Sig1(sn[4]), Ch(sn[4], sn[5], sn[6]), K[r], Wn[r])
        T2n = add(Sig0(sn[0]), Maj(sn[0], sn[1], sn[2]))
        T1f = add(sf[7], Sig1(sf[4]), Ch(sf[4], sf[5], sf[6]), K[r], Wf[r])
        T2f = add(Sig0(sf[0]), Maj(sf[0], sf[1], sf[2]))
        sn = [add(T1n, T2n), sn[0], sn[1], sn[2], add(sn[3], T1n), sn[4], sn[5], sn[6]]
        sf = [add(T1f, T2f), sf[0], sf[1], sf[2], add(sf[3], T1f), sf[4], sf[5], sf[6]]

    return Wf

def experiment_C():
    print("\n" + "=" * 80)
    print("EXPERIMENT C: Alternating cascade vs Wang carry agreement")
    print("=" * 80)

    N = 100
    wang_ca = [0.0] * 64
    alt_ca = [0.0] * 64

    # Also verify alternating chain correctness
    alt_correct = 0
    alt_total = 0

    for trial in range(N):
        Wn = [rw() for _ in range(16)]

        # Wang chain
        Wf_wang = wang_chain(Wn)
        ca_w = carry_agreement_pair(Wn, Wf_wang)
        for r in range(64): wang_ca[r] += ca_w[r]

        # Alternating chain
        Wf_alt = alternating_chain(Wn)
        ca_a = carry_agreement_pair(Wn, Wf_alt)
        for r in range(64): alt_ca[r] += ca_a[r]

        # Verify alternating chain: check δe on even rounds, δa on odd rounds
        sn = all_states(schedule(Wn))
        sf = all_states(schedule(Wf_alt))
        for r in range(2, 17):
            alt_total += 1
            if r % 2 == 0:
                de = (sn[r][4] ^ sf[r][4])
                if de == 0: alt_correct += 1
            else:
                da = (sn[r][0] ^ sf[r][0])
                if da == 0: alt_correct += 1

    for r in range(64):
        wang_ca[r] /= N
        alt_ca[r] /= N

    print(f"\n  Alternating chain verification: {alt_correct}/{alt_total} "
          f"({alt_correct/alt_total*100:.1f}%) correct δ=0")

    print(f"\n  {'r':>3} | {'Wang CA':>8} | {'Alt CA':>8} | {'Diff(Alt-Wang)':>14}")
    print("  " + "-" * 45)
    for r in range(28):
        diff = alt_ca[r] - wang_ca[r]
        marker = " *" if abs(diff) > 0.02 else ""
        print(f"  {r:3d} | {wang_ca[r]:.4f} | {alt_ca[r]:.4f} | {diff:+.4f}{marker}")

    # Summary ranges
    for label, r1, r2 in [("1..16", 1, 16), ("17..25", 17, 25), ("18..40", 18, 40)]:
        w_mean = sum(wang_ca[r] for r in range(r1, r2+1)) / (r2-r1+1)
        a_mean = sum(alt_ca[r] for r in range(r1, r2+1)) / (r2-r1+1)
        print(f"\n  Mean carry_agreement({label}): Wang={w_mean:.4f}, Alt={a_mean:.4f}, "
              f"Diff={a_mean-w_mean:+.4f}")


# ═══════════════════════════════════════════════════════════════
# EXPERIMENT D: Round-by-round healing analysis
# ═══════════════════════════════════════════════════════════════
def experiment_D(pairs):
    print("\n" + "=" * 80)
    print("EXPERIMENT D: Round-by-round healing analysis")
    print("=" * 80)

    n_pairs = len(pairs)
    all_new = [[0]*64 for _ in range(n_pairs)]
    all_healed = [[0]*64 for _ in range(n_pairs)]
    all_ratio = [[0.0]*64 for _ in range(n_pairs)]

    for pi, (Wn, Wf) in enumerate(pairs):
        sn = all_states(schedule(Wn))
        sf = all_states(schedule(Wf))
        for r in range(64):
            new_d = 0; healed_d = 0
            for w in range(8):
                d_cur = sn[r][w] ^ sf[r][w]
                d_nxt = sn[r+1][w] ^ sf[r+1][w]
                new_d += hw((~d_cur & MASK) & d_nxt)
                healed_d += hw(d_cur & (~d_nxt & MASK))
            all_new[pi][r] = new_d
            all_healed[pi][r] = healed_d
            all_ratio[pi][r] = healed_d / max(new_d, 1)

    # 1. Per-round data
    print(f"\n--- Per-round new_diffs and healed_diffs (mean ± std) ---")
    print(f"  {'r':>3} | {'new':>6} | {'healed':>8} | {'ratio':>6} | "
          f"{'heal>1.5 count':>14}")
    print("  " + "-" * 50)

    healing_rounds_strict = []  # ratio > 1.5 in ALL pairs
    for r in range(64):
        new_vals = [all_new[pi][r] for pi in range(n_pairs)]
        heal_vals = [all_healed[pi][r] for pi in range(n_pairs)]
        ratio_vals = [all_ratio[pi][r] for pi in range(n_pairs)]
        mn = sum(new_vals) / n_pairs
        mh = sum(heal_vals) / n_pairs
        mr = sum(ratio_vals) / n_pairs
        count_15 = sum(1 for pi in range(n_pairs) if all_ratio[pi][r] > 1.5)
        print(f"  {r:3d} | {mn:6.1f} | {mh:8.1f} | {mr:6.3f} | {count_15:>14}")

        if all(all_ratio[pi][r] > 1.5 for pi in range(n_pairs)):
            healing_rounds_strict.append(r)

    # 2. Healing rounds (ratio > 1.5 in ALL pairs)
    print(f"\n--- Healing rounds (ratio > 1.5 in ALL 7 pairs) ---")
    if healing_rounds_strict:
        print(f"  Rounds: {healing_rounds_strict}")
    else:
        print(f"  None found with ratio > 1.5 in ALL pairs.")

    # Relaxed: ratio > 1.0 in majority (>=5)
    healing_majority = []
    for r in range(64):
        count = sum(1 for pi in range(n_pairs) if all_healed[pi][r] > all_new[pi][r])
        if count >= 5:
            healing_majority.append((r, count))
    print(f"\n--- Rounds where healed > new in ≥5/7 pairs ---")
    if healing_majority:
        for r, c in healing_majority:
            vals_n = [all_new[pi][r] for pi in range(n_pairs)]
            vals_h = [all_healed[pi][r] for pi in range(n_pairs)]
            print(f"  r={r:2d}: {c}/7 pairs, mean_new={sum(vals_n)/n_pairs:.1f}, "
                  f"mean_healed={sum(vals_h)/n_pairs:.1f}")
    else:
        print("  None found.")

    # 3. For rounds with mean(healed) > mean(new): what bits are healing?
    print(f"\n--- Bit analysis for top healing rounds ---")
    # Find rounds with significant mean healing excess
    top_heal_rounds = []
    for r in range(64):
        mn = sum(all_new[pi][r] for pi in range(n_pairs)) / n_pairs
        mh = sum(all_healed[pi][r] for pi in range(n_pairs)) / n_pairs
        if mh > mn + 1:
            top_heal_rounds.append((r, mh - mn))
    top_heal_rounds.sort(key=lambda x: -x[1])

    reg_names = ['a','b','c','d','e','f','g','h']
    for r, excess in top_heal_rounds[:5]:
        print(f"\n  Round {r} (healing excess = {excess:.1f}):")
        # For each bit position, how many pairs have it healing?
        heal_by_pos = [[0]*32 for _ in range(8)]
        for pi, (Wn, Wf) in enumerate(pairs):
            sn = all_states(schedule(Wn))
            sf = all_states(schedule(Wf))
            for w in range(8):
                d_cur = sn[r][w] ^ sf[r][w]
                d_nxt = sn[r+1][w] ^ sf[r+1][w]
                healed = d_cur & (~d_nxt & MASK)
                for b in range(32):
                    if (healed >> b) & 1:
                        heal_by_pos[w][b] += 1
        # Show positions healed in >=3 pairs
        common = []
        for w in range(8):
            for b in range(32):
                if heal_by_pos[w][b] >= 3:
                    common.append((w, b, heal_by_pos[w][b]))
        if common:
            common.sort(key=lambda x: -x[2])
            for w, b, c in common[:10]:
                print(f"    reg={reg_names[w]} bit={b:2d}: heals in {c}/7 pairs")

    # 4. Carry agreement on healing rounds vs non-healing
    print(f"\n--- Carry agreement: healing rounds vs non-healing ---")
    heal_set = set(r for r, _ in top_heal_rounds[:10])
    non_heal_set = set(range(20, 60)) - heal_set

    all_ca_heal = []
    all_ca_nonheal = []
    for pi, (Wn, Wf) in enumerate(pairs):
        ca = carry_agreement_pair(Wn, Wf)
        for r in heal_set:
            if r < 64: all_ca_heal.append(ca[r])
        for r in non_heal_set:
            if r < 64: all_ca_nonheal.append(ca[r])

    if all_ca_heal and all_ca_nonheal:
        mh = sum(all_ca_heal) / len(all_ca_heal)
        mn = sum(all_ca_nonheal) / len(all_ca_nonheal)
        print(f"  Mean carry_agreement on healing rounds: {mh:.4f}")
        print(f"  Mean carry_agreement on non-healing rounds: {mn:.4f}")
        print(f"  Difference: {mh-mn:+.4f}")

    # 5. Connection to K[r] or W[r]?
    print(f"\n--- K[r] at healing rounds ---")
    for r, excess in top_heal_rounds[:10]:
        print(f"  r={r:2d}: K[r]=0x{K[r]:08x} (HW={hw(K[r]):2d}), excess={excess:.1f}")

    # 6. Random baseline: are specific healing rounds structural?
    print(f"\n--- Structural check: healing rounds in 100 random Wang pairs ---")
    heal_freq = [0] * 64
    for _ in range(100):
        Wn = [rw() for _ in range(16)]
        Wf = wang_chain(Wn)
        sn = all_states(schedule(Wn))
        sf = all_states(schedule(Wf))
        for r in range(64):
            new_d = 0; healed_d = 0
            for w in range(8):
                d_cur = sn[r][w] ^ sf[r][w]
                d_nxt = sn[r+1][w] ^ sf[r+1][w]
                new_d += hw((~d_cur & MASK) & d_nxt)
                healed_d += hw(d_cur & (~d_nxt & MASK))
            if healed_d > new_d:
                heal_freq[r] += 1

    print(f"  Frequency of healed>new in 100 random Wang pairs:")
    print(f"  {'r':>3} | {'freq':>5} | structural?")
    print("  " + "-" * 35)
    for r in range(64):
        structural = "YES" if heal_freq[r] > 60 else ("maybe" if heal_freq[r] > 45 else "no")
        if heal_freq[r] > 40 or r in heal_set:
            marker = " <-- healing round" if r in heal_set else ""
            print(f"  {r:3d} | {heal_freq[r]:5d} | {structural}{marker}")


# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════
def main():
    print("=" * 80)
    print("  ЗАДАНИЕ 11: Follow the Carry — Two Anomalies from Tomography")
    print("  Experiments B, C, D")
    print("=" * 80)

    # Load pairs
    search_output = None
    for f in ["task8_results.txt", "/tmp/task8_results.txt"]:
        if os.path.exists(f):
            search_output = f
            break

    if not search_output:
        print("ERROR: Cannot find task8_results.txt")
        sys.exit(1)

    pairs = parse_pairs(search_output)
    print(f"\nLoaded {len(pairs)} pairs from {search_output}")

    experiment_B(pairs)
    experiment_C()
    experiment_D(pairs)

    print("\n" + "=" * 80)
    print("  ALL EXPERIMENTS COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
