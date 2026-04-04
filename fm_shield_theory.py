#!/usr/bin/env python3
"""
Flow Mathematics v0.5 — Shield Theory

Goal: Classify when and where △ DIES (doesn't propagate).
Build the INVERSE of Theorem 1'' — not "how fast △ spreads"
but "how slowly can we MAKE it spread".

Key concept: Shield potential S(r) = probability that a △ entering
round r gets absorbed (doesn't reach round r+1).

If S(r) is high → △ dies at round r → differential path survives.
"""

import random
from collections import defaultdict

MASK32 = 0xFFFFFFFF

K = [0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2]
IV = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&MASK32
def sig0(x): return rotr(x,2)^rotr(x,13)^rotr(x,22)
def sig1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def ssig0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def ssig1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def ch(e,f,g): return (e&f)^(~e&g)&MASK32
def maj(a,b,c): return (a&b)^(a&c)^(b&c)
def hw(x): return bin(x&MASK32).count('1')

def gpk_char(a_bit, b_bit):
    if a_bit == 1 and b_bit == 1: return 'G'
    if a_bit == 0 and b_bit == 0: return 'K'
    return 'P'


# ================================================================
# PART 1: Per-bit shield analysis
#
# When does △ on bit k get ABSORBED?
#
# In SHA-256 round: e_new = d + T1
# carry[k] = MAJ(d[k], T1[k], carry[k-1])
# e_new[k] = d[k] XOR T1[k] XOR carry[k-1]
#
# △ on bit k of e_new means: for two messages M1, M2,
# e_new[k] differs. This happens when d[k] XOR T1[k] XOR carry[k-1]
# gives different results.
#
# Shield = when the difference CANCELS:
# d1[k] XOR T1_1[k] XOR c1[k-1] = d2[k] XOR T1_2[k] XOR c2[k-1]
# even though some inputs differ.
#
# Three shield mechanisms:
# S1: d[k] shields — d1[k]=d2[k] (shift register preserves, no diff)
# S2: T1[k] shields — T1 differences cancel internally
# S3: carry shields — carry[k-1] absorbs the difference
# ================================================================

print("=" * 70)
print("PART 1: Per-bit shield mechanisms in SHA-256")
print("=" * 70)

# Measure: for each bit k, when two Wang pairs have δe=0 but δa≠0,
# HOW does the δa difference get absorbed in the e-computation?
# Specifically: at which bits does carry ABSORB vs PROPAGATE the diff?

def full_trace_pair(W1, W2, R):
    """Run SHA-256 on both messages, return per-round per-bit details."""
    def expand(M):
        W = list(M) + [0]*(64-len(M))
        for i in range(16,64):
            W[i] = (ssig1(W[i-2])+W[i-7]+ssig0(W[i-15])+W[i-16])&MASK32
        return W

    Wexp1 = expand(W1)
    Wexp2 = expand(W2)

    a1,b1,c1,d1,e1,f1,g1,h1 = IV
    a2,b2,c2,d2,e2,f2,g2,h2 = IV

    rounds = []
    for r in range(R):
        # Message 1
        T1_1 = (h1+sig1(e1)+ch(e1,f1,g1)+K[r]+Wexp1[r])&MASK32
        T2_1 = (sig0(a1)+maj(a1,b1,c1))&MASK32

        # Message 2
        T1_2 = (h2+sig1(e2)+ch(e2,f2,g2)+K[r]+Wexp2[r])&MASK32
        T2_2 = (sig0(a2)+maj(a2,b2,c2))&MASK32

        # Per-bit analysis of e_new = d + T1
        shields = {'K_absorb': 0, 'G_absorb': 0, 'P_pass': 0, 'P_block': 0}
        bit_detail = []

        for k in range(32):
            d1k = (d1 >> k) & 1
            d2k = (d2 >> k) & 1
            t1k_1 = (T1_1 >> k) & 1
            t1k_2 = (T1_2 >> k) & 1

            gpk1 = gpk_char(d1k, t1k_1)
            gpk2 = gpk_char(d2k, t1k_2)

            # Does the input to this bit differ between M1 and M2?
            input_differs = (d1k != d2k) or (t1k_1 != t1k_2)

            # Result: does e_new[k] differ?
            e_new_1 = ((d1 + T1_1) >> k) & 1
            e_new_2 = ((d2 + T1_2) >> k) & 1
            output_differs = (e_new_1 != e_new_2)

            if input_differs and not output_differs:
                # SHIELD EVENT: input differed but output same
                if gpk1 == gpk2:
                    if gpk1 == 'K': shields['K_absorb'] += 1
                    elif gpk1 == 'G': shields['G_absorb'] += 1
                    else: shields['P_block'] += 1
                else:
                    shields['P_block'] += 1  # GPK changed but result same

            elif input_differs and output_differs:
                shields['P_pass'] += 1  # difference passed through

            bit_detail.append({
                'gpk1': gpk1, 'gpk2': gpk2,
                'in_diff': input_differs, 'out_diff': output_differs
            })

        # State update
        h1,g1,f1 = g1,f1,e1; e1=(d1+T1_1)&MASK32; d1,c1,b1=c1,b1,a1; a1=(T1_1+T2_1)&MASK32
        h2,g2,f2 = g2,f2,e2; e2=(d2+T1_2)&MASK32; d2,c2,b2=c2,b2,a2; a2=(T1_2+T2_2)&MASK32

        da_xor = a1 ^ a2
        de_xor = e1 ^ e2

        rounds.append({
            'r': r,
            'shields': shields,
            'hw_da': hw(da_xor),
            'hw_de': hw(de_xor),
            'T1_diff': hw(T1_1 ^ T1_2),
            'd_diff': hw(d1 ^ d2) if r > 0 else 0,
            'bit_detail': bit_detail,
        })

    return rounds


def wang_cascade(W_base, dW0, R):
    W1 = list(W_base)
    DW = [0]*16; DW[0] = dW0
    # Simple trace for cascade
    def quick_sha(M, rounds):
        W = list(M)+[0]*(64-len(M))
        for i in range(16,64): W[i]=(ssig1(W[i-2])+W[i-7]+ssig0(W[i-15])+W[i-16])&MASK32
        a,b,c,d,e,f,g,h = IV
        for r in range(rounds):
            T1=(h+sig1(e)+ch(e,f,g)+K[r]+W[r])&MASK32
            T2=(sig0(a)+maj(a,b,c))&MASK32
            h,g,f=g,f,e; e=(d+T1)&MASK32; d,c,b=c,b,a; a=(T1+T2)&MASK32
        return [a,b,c,d,e,f,g,h]

    for step in range(min(R-1,15)):
        wi = step+1
        W2 = [(W1[i]+DW[i])&MASK32 for i in range(16)]
        s1 = quick_sha(W1, step+2)
        s2 = quick_sha(W2, step+2)
        de = (s2[4]-s1[4])&MASK32
        DW[wi] = (-de)&MASK32

    W2 = [(W1[i]+DW[i])&MASK32 for i in range(16)]
    return W1, W2, DW


# Run analysis on many Wang pairs
random.seed(42)
N = 500
R = 16

# Aggregate shield statistics per round
shield_stats = defaultdict(lambda: defaultdict(list))

for trial in range(N):
    W_base = [random.randint(0, MASK32) for _ in range(16)]
    W1, W2, DW = wang_cascade(W_base, 1, R)
    rounds = full_trace_pair(W1, W2, R)

    for rd in rounds:
        r = rd['r']
        for stype in ['K_absorb', 'G_absorb', 'P_pass', 'P_block']:
            shield_stats[r][stype].append(rd['shields'][stype])
        shield_stats[r]['hw_da'].append(rd['hw_da'])
        shield_stats[r]['hw_de'].append(rd['hw_de'])
        shield_stats[r]['T1_diff'].append(rd['T1_diff'])
        shield_stats[r]['d_diff'].append(rd['d_diff'])

print(f"\n  N={N} Wang pairs, R={R}")
print(f"\n  {'r':>2} | {'HW(δa)':>6} {'HW(δe)':>6} | {'d_diff':>6} {'T1_diff':>7} | {'K_abs':>5} {'G_abs':>5} {'P_blk':>5} {'P_pass':>6} | {'S_pot':>5}")

for r in range(R):
    da = sum(shield_stats[r]['hw_da'])/N
    de = sum(shield_stats[r]['hw_de'])/N
    dd = sum(shield_stats[r]['d_diff'])/N
    dt1 = sum(shield_stats[r]['T1_diff'])/N
    k_abs = sum(shield_stats[r]['K_absorb'])/N
    g_abs = sum(shield_stats[r]['G_absorb'])/N
    p_blk = sum(shield_stats[r]['P_block'])/N
    p_pass = sum(shield_stats[r]['P_pass'])/N

    # Shield potential = fraction of input-diffs that get absorbed
    total_in_diff = k_abs + g_abs + p_blk + p_pass
    s_pot = (k_abs + g_abs + p_blk) / total_in_diff if total_in_diff > 0 else 1.0

    print(f"  {r:>2} | {da:>6.1f} {de:>6.1f} | {dd:>6.1f} {dt1:>7.1f} | {k_abs:>5.1f} {g_abs:>5.1f} {p_blk:>5.1f} {p_pass:>6.1f} | {s_pot:>5.2f}")


# ================================================================
# PART 2: Shield potential as function of δa (the "shadow")
#
# Key insight from GPK analysis:
# δe = 0 (Wang cascade), but δa ≠ 0.
# The δa difference enters e-computation through d = a[r-3].
# At round r: d_diff = δa[r-3] (from 3 rounds ago).
#
# Shield potential should DEPEND on HW(δa):
# Low HW(δa) → few input diffs → easy to shield → S high
# High HW(δa) → many input diffs → hard to shield → S low
# ================================================================

print()
print("=" * 70)
print("PART 2: Shield potential vs HW(δa)")
print("=" * 70)

# Group rounds by HW(δa) and measure shield potential
hw_da_buckets = defaultdict(lambda: {'shielded': 0, 'total': 0})

random.seed(123)
N2 = 1000

for trial in range(N2):
    W_base = [random.randint(0, MASK32) for _ in range(16)]
    W1, W2, DW = wang_cascade(W_base, 1, R)
    rounds = full_trace_pair(W1, W2, R)

    for rd in rounds:
        hw_da = rd['hw_da']
        bucket = (hw_da // 4) * 4  # group by 4
        shielded = rd['shields']['K_absorb'] + rd['shields']['G_absorb'] + rd['shields']['P_block']
        total = shielded + rd['shields']['P_pass']

        hw_da_buckets[bucket]['shielded'] += shielded
        hw_da_buckets[bucket]['total'] += total

print(f"\n  N={N2} pairs × {R} rounds")
print(f"\n  {'HW(δa)':>8} | {'S_potential':>11} | {'samples':>7}")

for bucket in sorted(hw_da_buckets.keys()):
    d = hw_da_buckets[bucket]
    if d['total'] > 0:
        s_pot = d['shielded'] / d['total']
        print(f"  {bucket:>4}-{bucket+3:<3} | {s_pot:>11.3f} | {d['total']:>7}")


# ================================================================
# PART 3: Critical question — do LOW-HW(δa) rounds have higher S?
# If yes → seek differential paths that keep δa small.
# This connects to Wang: his paths try to keep δa small.
# ================================================================

print()
print("=" * 70)
print("PART 3: Shield potential by round (with δa context)")
print("=" * 70)

# For each round, show S_pot AND the HW(δa) at that round
print(f"\n  {'r':>2} | {'E[HW(δa)]':>10} | {'S_pot':>5} | {'S_pot explanation':>40}")

for r in range(R):
    da = sum(shield_stats[r]['hw_da'])/N
    k_abs = sum(shield_stats[r]['K_absorb'])/N
    g_abs = sum(shield_stats[r]['G_absorb'])/N
    p_blk = sum(shield_stats[r]['P_block'])/N
    p_pass = sum(shield_stats[r]['P_pass'])/N
    total = k_abs + g_abs + p_blk + p_pass
    s_pot = (k_abs + g_abs + p_blk) / total if total > 0 else 1.0

    if total < 0.1:
        expl = "no diffs (δd=0, δT1=0)"
    elif s_pot > 0.95:
        expl = "TOTAL SHIELD (Wang cascade active)"
    elif s_pot > 0.7:
        expl = f"strong shield (few δd bits: {sum(shield_stats[r]['d_diff'])/N:.1f})"
    elif s_pot > 0.5:
        expl = f"partial shield"
    else:
        expl = f"weak shield (δa saturated at ~16)"

    print(f"  {r:>2} | {da:>10.1f} | {s_pot:>5.2f} | {expl}")


# ================================================================
# PART 4: Bit-position analysis — which bits shield most?
#
# Carry propagates from bit 0 upward. Lower bits have shorter
# carry chains → more predictable → better shields?
# ================================================================

print()
print("=" * 70)
print("PART 4: Per-bit shield rate (which bit positions shield best?)")
print("=" * 70)

# Collect per-bit shield events across all rounds 5-15 (where δa≠0)
bit_shield = [0] * 32
bit_total = [0] * 32

random.seed(456)
N3 = 500

for trial in range(N3):
    W_base = [random.randint(0, MASK32) for _ in range(16)]
    W1, W2, DW = wang_cascade(W_base, 1, 16)
    rounds = full_trace_pair(W1, W2, 16)

    for rd in rounds:
        r = rd['r']
        if r < 5 or r > 15:
            continue

        for k in range(32):
            bd = rd['bit_detail'][k]
            if bd['in_diff']:
                bit_total[k] += 1
                if not bd['out_diff']:
                    bit_shield[k] += 1

print(f"\n  N={N3} pairs, rounds 5-15")
print(f"\n  {'bit':>3} | {'shield_rate':>11} | {'total_diffs':>11} | note")

for k in range(32):
    if bit_total[k] > 0:
        rate = bit_shield[k] / bit_total[k]
        note = ""
        if k == 0: note = "← LSB (carry-free)"
        elif k <= 4: note = "← low bits (short carry)"
        elif k >= 28: note = "← high bits (long carry)"

        # Highlight anomalies
        if rate > 0.6: note += " ★ HIGH SHIELD"
        elif rate < 0.3: note += " ✗ LOW SHIELD"

        print(f"  {k:>3} | {rate:>11.3f} | {bit_total[k]:>11} | {note}")


# ================================================================
# PART 5: THE KEY QUESTION
# Can we MAXIMIZE shield potential by choosing δW carefully?
#
# Wang uses δW[0] = 1 (bit 0). Is there a better choice?
# Test: different single-bit δW[0] values, measure total shield.
# ================================================================

print()
print("=" * 70)
print("PART 5: Shield potential vs choice of δW[0]")
print("=" * 70)

N4 = 200
R_test = 16

results = []

for bit in range(32):
    dW0 = 1 << bit
    total_shield = 0
    total_diff = 0

    random.seed(789)
    for trial in range(N4):
        W_base = [random.randint(0, MASK32) for _ in range(16)]
        W1, W2, DW = wang_cascade(W_base, dW0, R_test)
        rounds = full_trace_pair(W1, W2, R_test)

        for rd in rounds:
            s = rd['shields']
            shielded = s['K_absorb'] + s['G_absorb'] + s['P_block']
            total = shielded + s['P_pass']
            total_shield += shielded
            total_diff += total

    s_pot = total_shield / total_diff if total_diff > 0 else 0
    results.append((bit, s_pot, total_diff))

# Sort by shield potential (best first)
results.sort(key=lambda x: -x[1])

print(f"\n  N={N4} pairs × {R_test} rounds per δW[0] choice")
print(f"\n  {'rank':>4} | {'bit':>3} | {'S_potential':>11} | note")

for rank, (bit, s_pot, total) in enumerate(results):
    note = ""
    if bit == 0: note = "← standard (Wang)"
    if rank == 0: note += " ★ BEST"
    if rank < 5: note += " (top 5)"

    print(f"  {rank+1:>4} | {bit:>3} | {s_pot:>11.4f} | {note}")
    if rank >= 9 and bit != 0:
        continue  # show top 10 + bit 0

# Find where bit 0 ranks
bit0_rank = next(i for i, (b, _, _) in enumerate(results) if b == 0) + 1
print(f"\n  Bit 0 (Wang standard) ranks #{bit0_rank} out of 32")
print(f"  Best: bit {results[0][0]} (S={results[0][1]:.4f})")
print(f"  Worst: bit {results[-1][0]} (S={results[-1][1]:.4f})")
print(f"  Range: {results[0][1] - results[-1][1]:.4f}")
