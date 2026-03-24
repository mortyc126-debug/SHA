#!/usr/bin/env python3
"""
ЗАДАНИЕ 10: Complete SHA-256 Tomography — Hypothesis-Free Observation
Full 6-layer analysis of 7 verified 17-round Wang pairs.
"""
import json, struct, os, sys, re, math
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

# ── Primitives ──
def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK
def sig0(x):    return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sig1(x):    return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)
def Sig0(x):    return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sig1(x):    return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return ((e & f) ^ ((~e & MASK) & g)) & MASK
def Maj(a, b, c): return ((a & b) ^ (a & c) ^ (b & c)) & MASK
def add(*args):
    s = 0
    for x in args: s = (s + x) & MASK
    return s
def hw(x): return bin(x & MASK).count('1')
def sub(a, b): return (a - b) & MASK

def schedule(M):
    W = list(M[:16])
    for i in range(16, 64):
        W.append(add(sig1(W[i-2]), W[i-7], sig0(W[i-15]), W[i-16]))
    return W

# ── Carry computation ──
def compute_carries(x, y):
    """Compute 32-bit carry vector for x + y.
    carry[i] = carry INTO bit position i+1 (i.e., carry OUT of bit i).
    Returns 32-bit integer where bit i = carry out of position i."""
    carries = 0
    c = 0  # carry into bit 0 = 0
    for i in range(32):
        xi = (x >> i) & 1
        yi = (y >> i) & 1
        c = (xi & yi) | (xi & c) | (yi & c)
        if c:
            carries |= (1 << i)
    return carries

# ── Full state computation with intermediates ──
def full_round_data(W, iv=None):
    """Compute all states + all intermediates for 64 rounds."""
    if iv is None: iv = H0
    s = list(iv)
    states = [tuple(s)]  # states[0..64]
    rounds = []  # rounds[0..63], each is a dict of intermediates

    for r in range(64):
        a, b, c, d, e, f, g, h = s

        # Function outputs
        sig1_e = Sig1(e)
        ch_efg = Ch(e, f, g)
        sig0_a = Sig0(a)
        maj_abc = Maj(a, b, c)

        # T1 sub-additions: h + Sig1(e), then + Ch, then + K[r], then + W[r]
        sum1 = add(h, sig1_e)
        sum2 = add(sum1, ch_efg)
        sum3 = add(sum2, K[r])
        sum4 = add(sum3, W[r])  # = T1

        # T2 sub-addition: Sig0(a) + Maj(a,b,c)
        sum5 = add(sig0_a, maj_abc)  # = T2

        T1 = sum4
        T2 = sum5

        # Carries for each sub-addition
        carry1 = compute_carries(h, sig1_e)
        carry2 = compute_carries(sum1, ch_efg)
        carry3 = compute_carries(sum2, K[r])
        carry4 = compute_carries(sum3, W[r])
        carry5 = compute_carries(sig0_a, maj_abc)

        # New state
        new_a = add(T1, T2)
        new_e = add(d, T1)
        carry_a = compute_carries(T1, T2)
        carry_e = compute_carries(d, T1)

        s = [new_a, a, b, c, new_e, e, f, g]
        states.append(tuple(s))

        rounds.append({
            'T1': T1, 'T2': T2,
            'Sig1': sig1_e, 'Ch': ch_efg, 'Sig0': sig0_a, 'Maj': maj_abc,
            'sum1': sum1, 'sum2': sum2, 'sum3': sum3, 'sum4': sum4, 'sum5': sum5,
            'carry': [carry1, carry2, carry3, carry4, carry5],
            'carry_a': carry_a, 'carry_e': carry_e,
        })

    return states, rounds


def tomography_pair(Wn, Wf, pair_idx):
    """Full 6-layer tomography for one pair."""
    Wn_sched = schedule(Wn)
    Wf_sched = schedule(Wf)

    states_n, rounds_n = full_round_data(Wn_sched)
    states_f, rounds_f = full_round_data(Wf_sched)

    result = {
        'pair_idx': pair_idx,
        'Wn': [f'{w:08x}' for w in Wn],
        'Wf': [f'{w:08x}' for w in Wf],
    }

    # ── LAYER 1: Full state dump ──
    layer1 = []
    for r in range(65):
        sn = states_n[r]
        sf = states_f[r]
        xor_diff = tuple(sn[i] ^ sf[i] for i in range(8))
        add_diff = tuple(sub(sf[i], sn[i]) for i in range(8))
        layer1.append({
            'r': r,
            'state_n': [f'{x:08x}' for x in sn],
            'state_f': [f'{x:08x}' for x in sf],
            'xor_diff': [f'{x:08x}' for x in xor_diff],
            'add_diff': [f'{x:08x}' for x in add_diff],
            'hw_xor': [hw(x) for x in xor_diff],
        })
    result['layer1'] = layer1

    # ── LAYER 2: Round intermediates ──
    layer2 = []
    for r in range(64):
        rn = rounds_n[r]
        rf = rounds_f[r]
        entry = {
            'r': r,
            'T1_n': rn['T1'], 'T1_f': rf['T1'],
            'T2_n': rn['T2'], 'T2_f': rf['T2'],
            'dT1_xor': rn['T1'] ^ rf['T1'],
            'dT2_xor': rn['T2'] ^ rf['T2'],
            'dT1_add': sub(rf['T1'], rn['T1']),
            'dT2_add': sub(rf['T2'], rn['T2']),
            'Sig1_n': rn['Sig1'], 'Sig1_f': rf['Sig1'],
            'Ch_n': rn['Ch'], 'Ch_f': rf['Ch'],
            'Sig0_n': rn['Sig0'], 'Sig0_f': rf['Sig0'],
            'Maj_n': rn['Maj'], 'Maj_f': rf['Maj'],
            'dSig1': rn['Sig1'] ^ rf['Sig1'],
            'dCh': rn['Ch'] ^ rf['Ch'],
            'dSig0': rn['Sig0'] ^ rf['Sig0'],
            'dMaj': rn['Maj'] ^ rf['Maj'],
        }
        layer2.append(entry)
    result['layer2'] = layer2

    # ── LAYER 3: Carry analysis ──
    layer3 = []
    for r in range(64):
        rn = rounds_n[r]
        rf = rounds_f[r]
        carry_data = []
        for si in range(5):
            cn = rn['carry'][si]
            cf = rf['carry'][si]
            dc = cn ^ cf
            carry_data.append({
                'carry_n': f'{cn:08x}', 'carry_f': f'{cf:08x}',
                'dcarry': f'{dc:08x}', 'hw_dcarry': hw(dc),
            })
        layer3.append({'r': r, 'carries': carry_data})
    result['layer3'] = layer3

    # ── LAYER 4: Message schedule ──
    layer4 = []
    for t in range(64):
        wn = Wn_sched[t]
        wf = Wf_sched[t]
        entry = {
            't': t,
            'W_n': f'{wn:08x}', 'W_f': f'{wf:08x}',
            'dW_xor': f'{(wn ^ wf):08x}',
            'dW_add': f'{sub(wf, wn):08x}',
            'hw_dW': hw(wn ^ wf),
        }
        if t >= 16:
            s0n = sig0(Wn_sched[t - 15])
            s0f = sig0(Wf_sched[t - 15])
            s1n = sig1(Wn_sched[t - 2])
            s1f = sig1(Wf_sched[t - 2])
            entry['sig0_n'] = f'{s0n:08x}'
            entry['sig0_f'] = f'{s0f:08x}'
            entry['dsig0'] = f'{(s0n ^ s0f):08x}'
            entry['sig1_n'] = f'{s1n:08x}'
            entry['sig1_f'] = f'{s1f:08x}'
            entry['dsig1'] = f'{(s1n ^ s1f):08x}'
        layer4.append(entry)
    result['layer4'] = layer4

    # ── LAYER 5: Aggregated statistics ──
    layer5 = []
    for r in range(65):
        sn = states_n[r]
        sf = states_f[r]
        xor_state = [sn[i] ^ sf[i] for i in range(8)]
        hw_state = sum(hw(x) for x in xor_state)
        hw_e = hw(xor_state[4])
        hw_a = hw(xor_state[0])

        entry = {
            'r': r,
            'hw_dstate': hw_state,
            'hw_de': hw_e,
            'hw_da': hw_a,
        }

        if r < 64:
            rn = rounds_n[r]
            rf = rounds_f[r]
            entry['hw_dT1'] = hw(rn['T1'] ^ rf['T1'])
            entry['hw_dT2'] = hw(rn['T2'] ^ rf['T2'])
            entry['hw_dCh'] = hw(rn['Ch'] ^ rf['Ch'])
            entry['hw_dMaj'] = hw(rn['Maj'] ^ rf['Maj'])

            # sum_carry_diff
            sum_carry_hw = sum(hw(rn['carry'][i] ^ rf['carry'][i]) for i in range(5))
            entry['sum_carry_diff'] = sum_carry_hw

            # carry_agreement: fraction of 160 carry bits where carry_n == carry_f
            total_agree = sum(32 - hw(rn['carry'][i] ^ rf['carry'][i]) for i in range(5))
            entry['carry_agreement'] = total_agree / 160.0

            # function_agreement: fraction of 64 bits where Ch_n==Ch_f AND Maj_n==Maj_f
            dCh = rn['Ch'] ^ rf['Ch']
            dMaj = rn['Maj'] ^ rf['Maj']
            # bit is agreeing if BOTH Ch and Maj agree at that position
            func_disagree = dCh | dMaj  # bit=1 if EITHER differs
            func_agree_count = 32 - hw(func_disagree)
            entry['function_agreement'] = func_agree_count / 32.0

            # linear_residual: HW(δT1 XOR (δh XOR δSig1 XOR δCh XOR δW))
            dh = sn[7] ^ sf[7]  # h register at start of this round = state[r][7]
            dSig1 = rn['Sig1'] ^ rf['Sig1']
            dW = Wn_sched[r] ^ Wf_sched[r]
            linear_approx = dh ^ dSig1 ^ dCh ^ dW
            dT1 = rn['T1'] ^ rf['T1']
            entry['linear_residual'] = hw(dT1 ^ linear_approx)

        layer5.append(entry)
    result['layer5'] = layer5

    # ── LAYER 6: Transition patterns ──
    layer6 = []
    for r in range(64):
        sn_cur = states_n[r]
        sf_cur = states_f[r]
        sn_nxt = states_n[r + 1]
        sf_nxt = states_f[r + 1]

        # Flatten 8×32=256 bit delta vectors
        new_diffs = 0
        healed_diffs = 0
        bit_flip_count = 0

        for w in range(8):
            d_cur = sn_cur[w] ^ sf_cur[w]
            d_nxt = sn_nxt[w] ^ sf_nxt[w]
            # bits that changed their delta status
            changed = d_cur ^ d_nxt
            bit_flip_count += hw(changed)
            # new diffs: was 0, now 1
            new_d = (~d_cur & MASK) & d_nxt
            new_diffs += hw(new_d)
            # healed diffs: was 1, now 0
            healed_d = d_cur & (~d_nxt & MASK)
            healed_diffs += hw(healed_d)

        stability = 256 - new_diffs - healed_diffs

        # bit_flip_map as 8 x 32-bit words
        bfm = []
        for w in range(8):
            d_cur = sn_cur[w] ^ sf_cur[w]
            d_nxt = sn_nxt[w] ^ sf_nxt[w]
            bfm.append(d_cur ^ d_nxt)

        layer6.append({
            'r': r,
            'new_diffs': new_diffs,
            'healed_diffs': healed_diffs,
            'stability': stability,
            'bit_flip_count': bit_flip_count,
            'bit_flip_map': [f'{x:08x}' for x in bfm],
        })
    result['layer6'] = layer6

    # Final hash
    hn = tuple(add(states_n[64][i], H0[i]) for i in range(8))
    hf = tuple(add(states_f[64][i], H0[i]) for i in range(8))
    dhash = tuple(hn[i] ^ hf[i] for i in range(8))
    result['hash_n'] = [f'{x:08x}' for x in hn]
    result['hash_f'] = [f'{x:08x}' for x in hf]
    result['dhash'] = [f'{x:08x}' for x in dhash]
    result['hw_dhash'] = [hw(x) for x in dhash]
    result['hw_dhash_total'] = sum(hw(x) for x in dhash)

    return result


def parse_pairs(filename):
    """Parse pairs from task8_results.txt."""
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


def aggregate_tables(all_results):
    """Compute aggregate Tables A, B, C and anomalies across all pairs."""
    n_pairs = len(all_results)
    output = []

    # ══════════════════════════════════════════════════════════════════
    # TABLE B: Transition patterns (PRIORITY)
    # ══════════════════════════════════════════════════════════════════
    output.append("=" * 100)
    output.append("TABLE B: Transition Patterns (mean ± std over {} pairs)".format(n_pairs))
    output.append("=" * 100)
    output.append(f"{'r':>3} | {'new_diffs':>12} | {'healed_diffs':>14} | {'stability':>12} | {'ratio n/h':>10}")
    output.append("-" * 60)

    table_b_data = []
    for r in range(64):
        new_vals = [res['layer6'][r]['new_diffs'] for res in all_results]
        heal_vals = [res['layer6'][r]['healed_diffs'] for res in all_results]
        stab_vals = [res['layer6'][r]['stability'] for res in all_results]

        mn = sum(new_vals) / n_pairs
        mh = sum(heal_vals) / n_pairs
        ms = sum(stab_vals) / n_pairs
        sn = math.sqrt(sum((x - mn)**2 for x in new_vals) / n_pairs) if n_pairs > 1 else 0
        sh = math.sqrt(sum((x - mh)**2 for x in heal_vals) / n_pairs) if n_pairs > 1 else 0
        ss = math.sqrt(sum((x - ms)**2 for x in stab_vals) / n_pairs) if n_pairs > 1 else 0

        ratio = f"{mn/mh:.3f}" if mh > 0 else "inf"
        output.append(f"{r:3d} | {mn:6.1f}±{sn:4.1f} | {mh:7.1f}±{sh:4.1f} | {ms:6.1f}±{ss:4.1f} | {ratio:>10}")
        table_b_data.append({'r': r, 'mean_new': mn, 'std_new': sn,
                             'mean_healed': mh, 'std_healed': sh,
                             'mean_stability': ms, 'std_stability': ss,
                             'ratio': mn / mh if mh > 0 else float('inf')})

    # ══════════════════════════════════════════════════════════════════
    # TABLE A: HW profile
    # ══════════════════════════════════════════════════════════════════
    output.append("")
    output.append("=" * 140)
    output.append("TABLE A: HW Profile (mean ± std over {} pairs)".format(n_pairs))
    output.append("=" * 140)
    output.append(f"{'r':>3} | {'HW(dstate)':>12} | {'HW(de)':>8} | {'HW(da)':>8} | "
                  f"{'HW(dT1)':>9} | {'HW(dT2)':>9} | {'carry_agr':>10} | "
                  f"{'func_agr':>9} | {'lin_resid':>10} | {'Σcarry_d':>9}")
    output.append("-" * 140)

    for r in range(65):
        vals_state = [res['layer5'][r]['hw_dstate'] for res in all_results]
        vals_e = [res['layer5'][r]['hw_de'] for res in all_results]
        vals_a = [res['layer5'][r]['hw_da'] for res in all_results]

        ms = sum(vals_state) / n_pairs
        me = sum(vals_e) / n_pairs
        ma = sum(vals_a) / n_pairs
        ss = math.sqrt(sum((x - ms)**2 for x in vals_state) / n_pairs)
        se = math.sqrt(sum((x - me)**2 for x in vals_e) / n_pairs)
        sa = math.sqrt(sum((x - ma)**2 for x in vals_a) / n_pairs)

        if r < 64:
            vals_t1 = [res['layer5'][r]['hw_dT1'] for res in all_results]
            vals_t2 = [res['layer5'][r]['hw_dT2'] for res in all_results]
            vals_ca = [res['layer5'][r]['carry_agreement'] for res in all_results]
            vals_fa = [res['layer5'][r]['function_agreement'] for res in all_results]
            vals_lr = [res['layer5'][r]['linear_residual'] for res in all_results]
            vals_cd = [res['layer5'][r]['sum_carry_diff'] for res in all_results]

            mt1 = sum(vals_t1) / n_pairs
            mt2 = sum(vals_t2) / n_pairs
            mca = sum(vals_ca) / n_pairs
            mfa = sum(vals_fa) / n_pairs
            mlr = sum(vals_lr) / n_pairs
            mcd = sum(vals_cd) / n_pairs

            st1 = math.sqrt(sum((x - mt1)**2 for x in vals_t1) / n_pairs)
            st2 = math.sqrt(sum((x - mt2)**2 for x in vals_t2) / n_pairs)
            sca = math.sqrt(sum((x - mca)**2 for x in vals_ca) / n_pairs)
            sfa = math.sqrt(sum((x - mfa)**2 for x in vals_fa) / n_pairs)
            slr = math.sqrt(sum((x - mlr)**2 for x in vals_lr) / n_pairs)
            scd = math.sqrt(sum((x - mcd)**2 for x in vals_cd) / n_pairs)

            output.append(
                f"{r:3d} | {ms:5.1f}±{ss:4.1f} | {me:3.1f}±{se:3.1f} | {ma:3.1f}±{sa:3.1f} | "
                f"{mt1:4.1f}±{st1:3.1f} | {mt2:4.1f}±{st2:3.1f} | "
                f"{mca:.3f}±{sca:.3f} | {mfa:.3f}±{sfa:.3f} | "
                f"{mlr:4.1f}±{slr:3.1f} | {mcd:4.1f}±{scd:3.1f}"
            )
        else:
            output.append(
                f"{r:3d} | {ms:5.1f}±{ss:4.1f} | {me:3.1f}±{se:3.1f} | {ma:3.1f}±{sa:3.1f} | "
                f"{'':>9} | {'':>9} | {'':>10} | {'':>9} | {'':>10} | {'':>9}"
            )

    # ══════════════════════════════════════════════════════════════════
    # TABLE C: Per-word HW(δH[i]) for final hash
    # ══════════════════════════════════════════════════════════════════
    output.append("")
    output.append("=" * 90)
    output.append("TABLE C: Per-word HW(δHash[i]) — 8 words × {} pairs".format(n_pairs))
    output.append("=" * 90)
    header = f"{'Pair':>5} |"
    for i in range(8):
        header += f" {'H'+str(i):>5}"
    header += f" | {'Total':>6}"
    output.append(header)
    output.append("-" * 70)

    all_hw_per_word = [[] for _ in range(8)]
    for res in all_results:
        line = f"{res['pair_idx']:5d} |"
        for i in range(8):
            line += f" {res['hw_dhash'][i]:5d}"
            all_hw_per_word[i].append(res['hw_dhash'][i])
        line += f" | {res['hw_dhash_total']:6d}"
        output.append(line)

    # Means
    line = f"{'Mean':>5} |"
    for i in range(8):
        m = sum(all_hw_per_word[i]) / n_pairs
        line += f" {m:5.1f}"
    total_mean = sum(sum(all_hw_per_word[i]) for i in range(8)) / n_pairs
    line += f" | {total_mean:6.1f}"
    output.append(line)

    # ══════════════════════════════════════════════════════════════════
    # ANOMALIES
    # ══════════════════════════════════════════════════════════════════
    output.append("")
    output.append("=" * 100)
    output.append("ANOMALIES")
    output.append("=" * 100)

    # 1. Rounds where linear_residual < 10
    output.append("\n--- Rounds where mean(linear_residual) < 10 (carry nearly invisible) ---")
    found = False
    for r in range(64):
        vals = [res['layer5'][r]['linear_residual'] for res in all_results]
        m = sum(vals) / n_pairs
        if m < 10:
            output.append(f"  r={r:2d}: mean_linear_residual={m:.1f}, values={vals}")
            found = True
    if not found:
        output.append("  None found.")

    # 2. Rounds where healed_diffs > new_diffs (δ DECREASING)
    output.append("\n--- Rounds where mean(healed_diffs) > mean(new_diffs) (δ DECREASING) ---")
    found = False
    for r in range(64):
        new_vals = [res['layer6'][r]['new_diffs'] for res in all_results]
        heal_vals = [res['layer6'][r]['healed_diffs'] for res in all_results]
        mn = sum(new_vals) / n_pairs
        mh = sum(heal_vals) / n_pairs
        if mh > mn:
            output.append(f"  r={r:2d}: mean_new={mn:.1f}, mean_healed={mh:.1f}, "
                          f"diff={mh-mn:.1f}")
            found = True
    if not found:
        output.append("  None found.")

    # 3. State bits that are δ=0 for more than 50% of rounds 17..64
    output.append("\n--- State bits with δ=0 for >50% of rounds 17..64 ---")
    # For each of 256 bits (8 words × 32 bits), count how many rounds r=17..64 have δ=0
    rounds_range = range(17, 65)  # states[17] to states[64]
    n_rounds = len(list(rounds_range))
    threshold = n_rounds * 0.5

    for pair_res in all_results:
        stable_bits = []
        for w in range(8):
            for b in range(32):
                count_zero = 0
                for r in rounds_range:
                    sn = int(pair_res['layer1'][r]['state_n'][w], 16)
                    sf = int(pair_res['layer1'][r]['state_f'][w], 16)
                    if ((sn ^ sf) >> b) & 1 == 0:
                        count_zero += 1
                if count_zero > threshold:
                    stable_bits.append((w, b, count_zero))
        if stable_bits:
            output.append(f"  Pair {pair_res['pair_idx']}: {len(stable_bits)} bits with δ=0 >50% of rounds 17..64")
            # Show top 10
            stable_bits.sort(key=lambda x: -x[2])
            for w, b, cnt in stable_bits[:10]:
                output.append(f"    word={w} bit={b:2d}: δ=0 in {cnt}/{n_rounds} rounds ({cnt/n_rounds*100:.0f}%)")
        else:
            output.append(f"  Pair {pair_res['pair_idx']}: no such bits found")

    # 4. Bit flip map patterns repeating between pairs
    output.append("\n--- Bit-flip-map patterns repeating between pairs ---")
    # For each round, check if any word of the bit_flip_map is identical across ≥2 pairs
    repeated_patterns = []
    for r in range(64):
        for w in range(8):
            patterns = {}
            for res in all_results:
                val = res['layer6'][r]['bit_flip_map'][w]
                if val != '00000000':
                    if val not in patterns:
                        patterns[val] = []
                    patterns[val].append(res['pair_idx'])
            for val, pair_indices in patterns.items():
                if len(pair_indices) >= 2:
                    repeated_patterns.append((r, w, val, pair_indices))

    if repeated_patterns:
        output.append(f"  Found {len(repeated_patterns)} repeated non-zero patterns:")
        for r, w, val, pids in repeated_patterns[:30]:
            output.append(f"    r={r:2d} word={w}: 0x{val} appears in pairs {pids}")
    else:
        output.append("  No exact repeated patterns found.")

    # 5. Additional: Message schedule HW(δW) profile
    output.append("\n--- Message schedule δW profile (mean HW across pairs) ---")
    output.append(f"{'t':>3} | {'mean HW(dW)':>12}")
    output.append("-" * 20)
    for t in range(64):
        vals = [int(res['layer4'][t]['hw_dW']) for res in all_results]
        m = sum(vals) / n_pairs
        s = math.sqrt(sum((x - m)**2 for x in vals) / n_pairs)
        marker = ""
        if m < 5:
            marker = " <<<< LOW"
        elif m == 0:
            marker = " <<<< ZERO"
        output.append(f"{t:3d} | {m:5.1f}±{s:4.1f}{marker}")

    return "\n".join(output)


def main():
    print("=" * 70)
    print("  ЗАДАНИЕ 10: Complete SHA-256 Tomography")
    print("  Hypothesis-Free Observation of 7 Wang Pairs")
    print("=" * 70)

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
    if len(pairs) < 7:
        print(f"WARNING: Expected 7 pairs, got {len(pairs)}")

    # Run tomography on each pair
    all_results = []
    for i, (Wn, Wf) in enumerate(pairs):
        print(f"\nProcessing pair {i}...")
        result = tomography_pair(Wn, Wf, i)
        all_results.append(result)
        print(f"  HW(δhash) = {result['hw_dhash_total']}")

    # Save full per-pair data as JSON
    print("\nSaving full per-pair data...")
    for res in all_results:
        fname = f"task10_pair{res['pair_idx']}.json"
        with open(fname, 'w') as f:
            json.dump(res, f, indent=1)
        print(f"  Saved {fname}")

    # Generate aggregate tables
    print("\nGenerating aggregate tables...")
    tables = aggregate_tables(all_results)

    with open("task10_tables.txt", 'w') as f:
        f.write(tables)
    print(f"Saved task10_tables.txt")

    # Also print tables
    print("\n" + tables)

    print("\n" + "=" * 70)
    print("  TOMOGRAPHY COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()
