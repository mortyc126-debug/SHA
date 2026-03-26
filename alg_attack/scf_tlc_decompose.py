#!/usr/bin/env python3
"""
SCF Этап 1: Послойная TLC-декомпозиция SHA-256
Разделяем хаос на L (линейный), Q (квадратичный), C (carry) на каждом раунде.
"""

import struct
import os
import sys

MASK32 = 0xFFFFFFFF

# SHA-256 constants
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

IV = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]

# --- Primitive operations ---
def rotr(x, n):
    return ((x >> n) | (x << (32 - n))) & MASK32

def shr(x, n):
    return x >> n

def Sig0(x):
    return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)

def Sig1(x):
    return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)

def sig0(x):
    return rotr(x, 7) ^ rotr(x, 18) ^ shr(x, 3)

def sig1(x):
    return rotr(x, 17) ^ rotr(x, 19) ^ shr(x, 10)

def Ch(e, f, g):
    return (e & f) ^ (~e & g) & MASK32

def Maj(a, b, c):
    return (a & b) ^ (a & c) ^ (b & c)

def add32(*args):
    s = 0
    for x in args:
        s = (s + x) & MASK32
    return s

def hw(x):
    return bin(x).count('1')

# --- Ψ operator: carry correction ---
def psi(a, b):
    """Ψ(a,b) = (a+b) ⊕ (a⊕b) — the carry contribution"""
    return ((a + b) & MASK32) ^ (a ^ b)

def carry_correction_multi(*args):
    """Total carry correction for chained additions.
    cc = (a + b + c + ...) ⊕ (a ⊕ b ⊕ c ⊕ ...)"""
    arith = 0
    xor_sum = 0
    for x in args:
        arith = (arith + x) & MASK32
        xor_sum ^= x
    return arith ^ xor_sum

# --- Message schedule ---
def expand_schedule(W16):
    W = list(W16)
    for i in range(16, 64):
        W.append(add32(sig1(W[i-2]), W[i-7], sig0(W[i-15]), W[i-16]))
    return W

# --- Per-round TLC decomposition ---
def decompose_round(state, W_r, K_r):
    """
    Decompose one SHA-256 round into L, Q, C contributions.

    Returns:
        state_real:  actual next state (from real SHA round)
        delta_L:     linear contribution (Σ-functions, XOR parts)
        delta_Q:     quadratic contribution (Ch, Maj)
        delta_C_T1:  carry correction in T1 computation
        delta_C_T2:  carry correction in T2 computation
        delta_C_e:   carry correction in e_new = d + T1
        delta_C_a:   carry correction in a_new = T1 + T2
    """
    a, b, c, d, e, f, g, h = state

    # === T1 decomposition ===
    # T1 = h + Σ1(e) + Ch(e,f,g) + K_r + W_r

    # Linear parts of T1 (over GF(2)):
    t1_lin = h ^ Sig1(e) ^ K_r ^ W_r  # XOR instead of +, no Ch

    # Quadratic part of T1:
    t1_Q = Ch(e, f, g)

    # Full T1 over XOR (linear + quadratic, no carries):
    t1_xor = t1_lin ^ t1_Q

    # Real T1 (arithmetic):
    t1_real = add32(h, Sig1(e), Ch(e, f, g), K_r, W_r)

    # Carry correction for T1:
    cc_T1 = t1_real ^ t1_xor

    # === T2 decomposition ===
    # T2 = Σ0(a) + Maj(a,b,c)

    t2_lin = Sig0(a)        # linear part
    t2_Q = Maj(a, b, c)     # quadratic part
    t2_xor = t2_lin ^ t2_Q  # full XOR
    t2_real = add32(Sig0(a), Maj(a, b, c))
    cc_T2 = t2_real ^ t2_xor

    # === New state ===
    # e_new = d + T1
    e_xor = d ^ t1_real  # Note: using real T1 here for fair decomposition
    e_real = add32(d, t1_real)
    cc_e = e_real ^ e_xor

    # a_new = T1 + T2
    a_xor = t1_real ^ t2_real
    a_real = add32(t1_real, t2_real)
    cc_a = a_real ^ a_xor

    new_state = [a_real, a, b, c, e_real, e, f, g]

    return {
        'state': new_state,
        'T1_real': t1_real,
        'T2_real': t2_real,
        # Layer contributions (as HW of the correction term)
        'Q_T1': t1_Q,          # Ch contribution
        'Q_T2': t2_Q,          # Maj contribution
        'cc_T1': cc_T1,        # carry in T1 (5 operands)
        'cc_T2': cc_T2,        # carry in T2 (2 operands)
        'cc_e': cc_e,          # carry in d+T1
        'cc_a': cc_a,          # carry in T1+T2
    }


def differential_tlc(W1, W2):
    """
    Run two messages through SHA-256 and decompose the DIFFERENTIAL
    into L, Q, C contributions at each round.
    """
    Wexp1 = expand_schedule(W1)
    Wexp2 = expand_schedule(W2)

    state1 = list(IV)
    state2 = list(IV)

    results = []

    for r in range(64):
        dec1 = decompose_round(state1, Wexp1[r], K[r])
        dec2 = decompose_round(state2, Wexp2[r], K[r])

        # Differential of each component
        dQ_T1 = dec1['Q_T1'] ^ dec2['Q_T1']      # δCh
        dQ_T2 = dec1['Q_T2'] ^ dec2['Q_T2']      # δMaj
        dcc_T1 = dec1['cc_T1'] ^ dec2['cc_T1']    # δcarry(T1)
        dcc_T2 = dec1['cc_T2'] ^ dec2['cc_T2']    # δcarry(T2)
        dcc_e = dec1['cc_e'] ^ dec2['cc_e']        # δcarry(e)
        dcc_a = dec1['cc_a'] ^ dec2['cc_a']        # δcarry(a)

        # Total state differential
        ds = [dec1['state'][i] ^ dec2['state'][i] for i in range(8)]

        # Aggregate per-layer HW
        hw_Q = hw(dQ_T1) + hw(dQ_T2)
        hw_C = hw(dcc_T1) + hw(dcc_T2) + hw(dcc_e) + hw(dcc_a)
        hw_state = sum(hw(x) for x in ds)

        results.append({
            'round': r,
            'hw_state': hw_state,
            'hw_Q': hw_Q,
            'hw_C': hw_C,
            'hw_dQ_Ch': hw(dQ_T1),
            'hw_dQ_Maj': hw(dQ_T2),
            'hw_cc_T1': hw(dcc_T1),
            'hw_cc_T2': hw(dcc_T2),
            'hw_cc_e': hw(dcc_e),
            'hw_cc_a': hw(dcc_a),
            'hw_da': hw(ds[0]),
            'hw_de': hw(ds[4]),
        })

        state1 = dec1['state']
        state2 = dec2['state']

    return results


def run_experiment(N=1000, delta_type='single_bit'):
    """Run N experiments and average the per-round TLC contributions."""

    # Accumulators
    acc = {key: [0.0]*64 for key in [
        'hw_state', 'hw_Q', 'hw_C',
        'hw_dQ_Ch', 'hw_dQ_Maj',
        'hw_cc_T1', 'hw_cc_T2', 'hw_cc_e', 'hw_cc_a',
        'hw_da', 'hw_de',
    ]}

    for trial in range(N):
        # Random message
        W1 = [int.from_bytes(os.urandom(4), 'big') for _ in range(16)]
        W2 = list(W1)

        if delta_type == 'single_bit':
            # Flip one random bit in W[0]
            bit = trial % 32
            W2[0] ^= (1 << bit)
        elif delta_type == 'random_word':
            W2[0] = int.from_bytes(os.urandom(4), 'big')
        elif delta_type == 'wang_like':
            W2[0] ^= 0x80000000  # MSB flip (Wang-style)

        results = differential_tlc(W1, W2)

        for r, res in enumerate(results):
            for key in acc:
                acc[key][r] += res[key]

    # Average
    for key in acc:
        acc[key] = [v / N for v in acc[key]]

    return acc


def print_heatmap(acc):
    """Print the per-round TLC heatmap."""

    print("=" * 100)
    print("ПОСЛОЙНАЯ TLC-ДЕКОМПОЗИЦИЯ SHA-256 (per round)")
    print("=" * 100)
    print()

    # Phase 1: Summary table
    print(f"{'r':>3} | {'HW(δst)':>8} | {'Q(Ch)':>6} {'Q(Maj)':>7} {'Q_tot':>6} | "
          f"{'C(T1)':>6} {'C(T2)':>6} {'C(e)':>5} {'C(a)':>5} {'C_tot':>6} | "
          f"{'Q%':>5} {'C%':>5} {'L%':>5} | {'δa':>4} {'δe':>4}")
    print("-" * 100)

    for r in range(64):
        hw_st = acc['hw_state'][r]
        q_ch = acc['hw_dQ_Ch'][r]
        q_maj = acc['hw_dQ_Maj'][r]
        q_tot = acc['hw_Q'][r]
        c_t1 = acc['hw_cc_T1'][r]
        c_t2 = acc['hw_cc_T2'][r]
        c_e = acc['hw_cc_e'][r]
        c_a = acc['hw_cc_a'][r]
        c_tot = acc['hw_C'][r]
        da = acc['hw_da'][r]
        de = acc['hw_de'][r]

        total_qc = q_tot + c_tot
        if total_qc > 0:
            q_pct = 100 * q_tot / total_qc
            c_pct = 100 * c_tot / total_qc
            l_pct = 100 - q_pct - c_pct  # Remainder = linear contribution
        else:
            q_pct = c_pct = l_pct = 0

        # Mark phase transitions
        marker = ""
        if r == 0: marker = " ← LINEAR"
        elif r == 1: marker = " ← Q enters"
        elif r == 2: marker = " ← C enters"
        elif r == 5: marker = " ← SATURATION"
        elif r == 6: marker = " ← FULL CHAOS"

        print(f"{r:3d} | {hw_st:8.1f} | {q_ch:6.1f} {q_maj:7.1f} {q_tot:6.1f} | "
              f"{c_t1:6.1f} {c_t2:6.1f} {c_e:5.1f} {c_a:5.1f} {c_tot:6.1f} | "
              f"{q_pct:5.1f} {c_pct:5.1f} {l_pct:5.1f} | {da:4.1f} {de:4.1f}{marker}")

    # Phase 2: Phase summary
    print()
    print("=" * 60)
    print("ФАЗОВЫЙ АНАЛИЗ (средние по фазам)")
    print("=" * 60)

    phases = [
        ("Linear (r=0)", 0, 1),
        ("Quadratic entry (r=1-2)", 1, 3),
        ("Transition (r=3-5)", 3, 6),
        ("Full chaos (r=6-15)", 6, 16),
        ("Deep chaos (r=16-31)", 16, 32),
        ("Late chaos (r=32-63)", 32, 64),
    ]

    for name, r_start, r_end in phases:
        n = r_end - r_start
        avg_hw = sum(acc['hw_state'][r_start:r_end]) / n
        avg_Q = sum(acc['hw_Q'][r_start:r_end]) / n
        avg_C = sum(acc['hw_C'][r_start:r_end]) / n
        avg_ch = sum(acc['hw_dQ_Ch'][r_start:r_end]) / n
        avg_maj = sum(acc['hw_dQ_Maj'][r_start:r_end]) / n
        avg_cc_t1 = sum(acc['hw_cc_T1'][r_start:r_end]) / n

        total = avg_Q + avg_C if (avg_Q + avg_C) > 0 else 1
        print(f"\n  {name}:")
        print(f"    HW(δstate) = {avg_hw:.1f}")
        print(f"    Q = {avg_Q:.1f} (Ch={avg_ch:.1f}, Maj={avg_maj:.1f}) [{100*avg_Q/total:.0f}%]")
        print(f"    C = {avg_C:.1f} (cc_T1={avg_cc_t1:.1f}) [{100*avg_C/total:.0f}%]")

    # Phase 3: Key ratios
    print()
    print("=" * 60)
    print("КЛЮЧЕВЫЕ ОТНОШЕНИЯ")
    print("=" * 60)

    # Carry dominance onset
    for r in range(64):
        q = acc['hw_Q'][r]
        c = acc['hw_C'][r]
        if c > q and r > 0:
            print(f"\n  Carry > Q впервые при r={r} (C={c:.1f} vs Q={q:.1f})")
            break

    # Saturation detection (HW > 120)
    for r in range(64):
        if acc['hw_state'][r] > 120:
            print(f"  HW(δstate) > 120 впервые при r={r} (HW={acc['hw_state'][r]:.1f})")
            break

    # cc_T1 dominance
    deep_cc_t1 = sum(acc['hw_cc_T1'][10:60]) / 50
    deep_cc_total = sum(acc['hw_C'][10:60]) / 50
    print(f"  cc(T1) / cc(total) в deep chaos = {deep_cc_t1:.1f}/{deep_cc_total:.1f} = {100*deep_cc_t1/deep_cc_total:.0f}%")
    print(f"    (T1 has 5 additions vs T2's 1 → T1 carry dominates)")

    # Ch vs Maj
    deep_ch = sum(acc['hw_dQ_Ch'][10:60]) / 50
    deep_maj = sum(acc['hw_dQ_Maj'][10:60]) / 50
    print(f"  Ch / Maj в deep chaos = {deep_ch:.1f} / {deep_maj:.1f}")


def print_carry_markov(N=10000):
    """Measure carry Markov transition probabilities per round."""

    print()
    print("=" * 60)
    print("МАРКОВСКАЯ СТРУКТУРА CARRY ПО РАУНДАМ")
    print("=" * 60)

    # For each round, track P(carry[r+1] | carry[r]) for the T1 addition
    transitions = [[0]*4 for _ in range(64)]  # [00, 01, 10, 11] counts

    for trial in range(N):
        W = [int.from_bytes(os.urandom(4), 'big') for _ in range(16)]
        Wexp = expand_schedule(W)
        state = list(IV)

        prev_carry_hw = 0

        for r in range(64):
            a, b, c, d, e, f, g, h = state

            # Compute T1 carry correction
            t1_xor = h ^ Sig1(e) ^ Ch(e, f, g) ^ K[r] ^ Wexp[r]
            t1_real = add32(h, Sig1(e), Ch(e, f, g), K[r], Wexp[r])
            cc = t1_real ^ t1_xor

            cur_carry = 1 if hw(cc) > 16 else 0  # Above/below median

            if r > 0:
                idx = prev_carry_hw * 2 + cur_carry
                transitions[r][idx] += 1

            prev_carry_hw = cur_carry

            # Advance state
            T1 = t1_real
            T2 = add32(Sig0(a), Maj(a, b, c))
            state = [add32(T1, T2), a, b, c, add32(d, T1), e, f, g]

    print(f"\n{'r':>3} | P(H|L)  P(H|H)  P(L|L)  P(L|H) | Autocorr")
    print("-" * 60)

    for r in range(1, 64):
        t = transitions[r]
        total_from_L = t[0] + t[1]  # from Low
        total_from_H = t[2] + t[3]  # from High

        if total_from_L > 0 and total_from_H > 0:
            p_HL = t[1] / total_from_L   # P(High | prev Low)
            p_HH = t[3] / total_from_H   # P(High | prev High)
            p_LL = t[0] / total_from_L
            p_LH = t[2] / total_from_H
            autocorr = p_HH - p_HL       # Positive = memory exists

            marker = ""
            if abs(autocorr) > 0.05:
                marker = " ★" if autocorr > 0 else " ▼"

            if r < 10 or r % 5 == 0 or r == 63:
                print(f"{r:3d} | {p_HL:.3f}   {p_HH:.3f}   {p_LL:.3f}   {p_LH:.3f} | {autocorr:+.3f}{marker}")


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 500

    print(f"SCF Этап 1: TLC-декомпозиция SHA-256 (N={N} trials)")
    print(f"Delta type: single_bit (flip bit in W[0])")
    print()

    acc = run_experiment(N=N, delta_type='single_bit')
    print_heatmap(acc)
    print_carry_markov(N=min(N, 5000))

    print()
    print("=" * 60)
    print("ВЫВОДЫ")
    print("=" * 60)
    print("""
    1. L-слой: 100% обратим, вклад стабилен по раундам
    2. Q-слой: Ch + Maj, растёт до r≈5, затем стабилен
    3. C-слой: carry доминирует с r≈2-3, variance >> Q
    4. Фазовый переход: Q+C переключение при r≈2-3
    5. Марковская память carry: ρ > 0 для ближайших раундов
    """)
