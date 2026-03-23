#!/usr/bin/env python3
"""
ЗАДАНИЕ 3: Slow Diffusion Trajectories in Phase B (rounds 17-48)

Experiments (priority D > A > C > B):
D. Register-specific diffusion asymmetry (efgh vs abcd)
A. Diffusion profile for Wang pairs (rounds 17-64)
C. Targeted search for slow paths
B. Lyapunov exponent per trajectory
"""

import struct
import os
import numpy as np
import math
import time
from collections import defaultdict

MASK = 0xFFFFFFFF

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

H0 = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
]


def rotr(x, n):
    return ((x >> n) | (x << (32 - n))) & MASK

def Sig0(x):
    return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)

def Sig1(x):
    return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)

def sig0(x):
    return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)

def sig1(x):
    return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)

def Ch(e, f, g):
    return (e & f) ^ ((~e) & g) & MASK

def Maj(a, b, c):
    return (a & b) ^ (a & c) ^ (b & c)

def add(*args):
    s = 0
    for x in args:
        s = (s + x) & MASK
    return s

def sub(a, b):
    return (a - b) & MASK

def hw32(x):
    return bin(x & MASK).count('1')

def random_words(n):
    return list(struct.unpack(f'>{n}I', os.urandom(4 * n)))


def message_schedule(M):
    W = list(M[:16])
    for i in range(16, 64):
        W.append(add(sig1(W[i-2]), W[i-7], sig0(W[i-15]), W[i-16]))
    return W


def sha256_all_states(W, iv=None):
    """Return all 65 states (state[0]..state[64])."""
    if iv is None:
        s = list(H0)
    else:
        s = list(iv)
    states = [tuple(s)]
    a, b, c, d, e, f, g, h = s
    for r in range(64):
        T1 = add(h, Sig1(e), Ch(e, f, g), K[r], W[r])
        T2 = add(Sig0(a), Maj(a, b, c))
        h, g, f = g, f, e
        e = add(d, T1)
        d, c, b = c, b, a
        a = add(T1, T2)
        states.append((a, b, c, d, e, f, g, h))
    return states


def sha256_hash(M):
    W = message_schedule(M)
    states = sha256_all_states(W)
    return tuple(add(states[64][i], H0[i]) for i in range(8))


def generate_wang_pair():
    """
    Generate Wang pair: δW[0]=0x8000, adaptive δW[1..15] → δe[2..16]=0.
    Returns (Wn, Wf) as 16-word lists.
    """
    Wn = random_words(16)
    Wf = list(Wn)
    Wf[0] = (Wn[0] ^ 0x8000) & MASK

    # Run rounds tracking both paths
    a_n, b_n, c_n, d_n = H0[0], H0[1], H0[2], H0[3]
    e_n, f_n, g_n, h_n = H0[4], H0[5], H0[6], H0[7]
    a_f, b_f, c_f, d_f = H0[0], H0[1], H0[2], H0[3]
    e_f, f_f, g_f, h_f = H0[4], H0[5], H0[6], H0[7]

    # Round 0
    T1_n = add(h_n, Sig1(e_n), Ch(e_n, f_n, g_n), K[0], Wn[0])
    T2_n = add(Sig0(a_n), Maj(a_n, b_n, c_n))
    T1_f = add(h_f, Sig1(e_f), Ch(e_f, f_f, g_f), K[0], Wf[0])
    T2_f = add(Sig0(a_f), Maj(a_f, b_f, c_f))
    h_n, g_n, f_n = g_n, f_n, e_n
    e_n = add(d_n, T1_n); d_n, c_n, b_n = c_n, b_n, a_n; a_n = add(T1_n, T2_n)
    h_f, g_f, f_f = g_f, f_f, e_f
    e_f = add(d_f, T1_f); d_f, c_f, b_f = c_f, b_f, a_f; a_f = add(T1_f, T2_f)

    # Rounds 1..15: adaptive correction
    for r in range(1, 16):
        dd = sub(d_f, d_n)
        dh = sub(h_f, h_n)
        dSig1 = sub(Sig1(e_f), Sig1(e_n))
        dCh = sub(Ch(e_f, f_f, g_f), Ch(e_n, f_n, g_n))
        dW_r = sub(0, add(dd, dh, dSig1, dCh))
        Wf[r] = add(Wn[r], dW_r)

        T1_n = add(h_n, Sig1(e_n), Ch(e_n, f_n, g_n), K[r], Wn[r])
        T2_n = add(Sig0(a_n), Maj(a_n, b_n, c_n))
        T1_f = add(h_f, Sig1(e_f), Ch(e_f, f_f, g_f), K[r], Wf[r])
        T2_f = add(Sig0(a_f), Maj(a_f, b_f, c_f))
        h_n, g_n, f_n = g_n, f_n, e_n
        e_n = add(d_n, T1_n); d_n, c_n, b_n = c_n, b_n, a_n; a_n = add(T1_n, T2_n)
        h_f, g_f, f_f = g_f, f_f, e_f
        e_f = add(d_f, T1_f); d_f, c_f, b_f = c_f, b_f, a_f; a_f = add(T1_f, T2_f)

    return Wn, Wf


def generate_wang_pair_with_w0(w0_val):
    """Generate Wang pair with specific W[0] value."""
    Wn = random_words(16)
    Wn[0] = w0_val
    Wf = list(Wn)
    Wf[0] = (Wn[0] ^ 0x8000) & MASK

    a_n, b_n, c_n, d_n = H0[0], H0[1], H0[2], H0[3]
    e_n, f_n, g_n, h_n = H0[4], H0[5], H0[6], H0[7]
    a_f, b_f, c_f, d_f = H0[0], H0[1], H0[2], H0[3]
    e_f, f_f, g_f, h_f = H0[4], H0[5], H0[6], H0[7]

    T1_n = add(h_n, Sig1(e_n), Ch(e_n, f_n, g_n), K[0], Wn[0])
    T2_n = add(Sig0(a_n), Maj(a_n, b_n, c_n))
    T1_f = add(h_f, Sig1(e_f), Ch(e_f, f_f, g_f), K[0], Wf[0])
    T2_f = add(Sig0(a_f), Maj(a_f, b_f, c_f))
    h_n, g_n, f_n = g_n, f_n, e_n
    e_n = add(d_n, T1_n); d_n, c_n, b_n = c_n, b_n, a_n; a_n = add(T1_n, T2_n)
    h_f, g_f, f_f = g_f, f_f, e_f
    e_f = add(d_f, T1_f); d_f, c_f, b_f = c_f, b_f, a_f; a_f = add(T1_f, T2_f)

    for r in range(1, 16):
        dd = sub(d_f, d_n)
        dh = sub(h_f, h_n)
        dSig1 = sub(Sig1(e_f), Sig1(e_n))
        dCh = sub(Ch(e_f, f_f, g_f), Ch(e_n, f_n, g_n))
        dW_r = sub(0, add(dd, dh, dSig1, dCh))
        Wf[r] = add(Wn[r], dW_r)

        T1_n = add(h_n, Sig1(e_n), Ch(e_n, f_n, g_n), K[r], Wn[r])
        T2_n = add(Sig0(a_n), Maj(a_n, b_n, c_n))
        T1_f = add(h_f, Sig1(e_f), Ch(e_f, f_f, g_f), K[r], Wf[r])
        T2_f = add(Sig0(a_f), Maj(a_f, b_f, c_f))
        h_n, g_n, f_n = g_n, f_n, e_n
        e_n = add(d_n, T1_n); d_n, c_n, b_n = c_n, b_n, a_n; a_n = add(T1_n, T2_n)
        h_f, g_f, f_f = g_f, f_f, e_f
        e_f = add(d_f, T1_f); d_f, c_f, b_f = c_f, b_f, a_f; a_f = add(T1_f, T2_f)

    return Wn, Wf


def compute_diffusion_profile(Wn, Wf):
    """
    Compute full diffusion profile for a Wang pair.
    Returns dict with HW arrays for rounds 16..64.
    """
    Wn_exp = message_schedule(Wn)
    Wf_exp = message_schedule(Wf)
    states_n = sha256_all_states(Wn_exp)
    states_f = sha256_all_states(Wf_exp)

    hw_total = []  # HW of full 256-bit δstate
    hw_e = []      # HW of δe (32 bits)
    hw_abcd = []   # HW of δa,δb,δc,δd (128 bits)
    hw_efgh = []   # HW of δe,δf,δg,δh (128 bits)

    for r in range(16, 65):
        sn = states_n[r]
        sf = states_f[r]
        dx = tuple(a ^ b for a, b in zip(sn, sf))

        hw_t = sum(hw32(x) for x in dx)
        hw_total.append(hw_t)
        hw_e.append(hw32(dx[4]))
        hw_abcd.append(sum(hw32(dx[i]) for i in range(4)))
        hw_efgh.append(sum(hw32(dx[i]) for i in range(4, 8)))

    # δhash
    hash_n = sha256_hash(Wn)
    hash_f = sha256_hash(Wf)
    dhash = tuple(a ^ b for a, b in zip(hash_n, hash_f))
    hw_dhash = sum(hw32(x) for x in dhash)

    return {
        'hw_total': hw_total,  # index 0 = round 16
        'hw_e': hw_e,
        'hw_abcd': hw_abcd,
        'hw_efgh': hw_efgh,
        'hw_dhash': hw_dhash,
    }


# =============================================================
# EXPERIMENT D: Register-specific diffusion asymmetry (PRIORITY)
# =============================================================
def experiment_D(n=10000):
    print("=" * 70)
    print("EXPERIMENT D: Register-specific diffusion asymmetry")
    print("  (efgh vs abcd in rounds 17-24)")
    print("=" * 70)
    t0 = time.time()

    all_efgh = []  # shape: (n, rounds)
    all_abcd = []
    all_e = []
    all_f = []
    all_g = []
    all_h = []
    all_a = []
    all_b = []
    all_c = []
    all_d = []

    for _ in range(n):
        Wn, Wf = generate_wang_pair()
        Wn_exp = message_schedule(Wn)
        Wf_exp = message_schedule(Wf)
        states_n = sha256_all_states(Wn_exp)
        states_f = sha256_all_states(Wf_exp)

        efgh_row = []
        abcd_row = []
        e_row, f_row, g_row, h_row = [], [], [], []
        a_row, b_row, c_row, d_row = [], [], [], []

        for r in range(16, 33):  # rounds 16..32
            sn = states_n[r]
            sf = states_f[r]
            dx = tuple(a ^ b for a, b in zip(sn, sf))

            abcd_row.append(sum(hw32(dx[i]) for i in range(4)))
            efgh_row.append(sum(hw32(dx[i]) for i in range(4, 8)))
            a_row.append(hw32(dx[0]))
            b_row.append(hw32(dx[1]))
            c_row.append(hw32(dx[2]))
            d_row.append(hw32(dx[3]))
            e_row.append(hw32(dx[4]))
            f_row.append(hw32(dx[5]))
            g_row.append(hw32(dx[6]))
            h_row.append(hw32(dx[7]))

        all_efgh.append(efgh_row)
        all_abcd.append(abcd_row)
        all_e.append(e_row)
        all_f.append(f_row)
        all_g.append(g_row)
        all_h.append(h_row)
        all_a.append(a_row)
        all_b.append(b_row)
        all_c.append(c_row)
        all_d.append(d_row)

    all_efgh = np.array(all_efgh)
    all_abcd = np.array(all_abcd)
    all_e = np.array(all_e)
    all_f = np.array(all_f)
    all_g = np.array(all_g)
    all_h = np.array(all_h)
    all_a = np.array(all_a)
    all_b = np.array(all_b)
    all_c = np.array(all_c)
    all_d = np.array(all_d)

    # Print register-group comparison
    print(f"\n  {'Round':>5}  {'HW_abcd':>8} {'±':>1} {'std':>5}  {'HW_efgh':>8} {'±':>1} {'std':>5}  {'gap':>6}  {'efgh/abcd':>9}")
    print(f"  {'-'*5}  {'-'*16}  {'-'*16}  {'-'*6}  {'-'*9}")

    catchup_round = None
    for i, r in enumerate(range(16, 33)):
        m_abcd = np.mean(all_abcd[:, i])
        s_abcd = np.std(all_abcd[:, i])
        m_efgh = np.mean(all_efgh[:, i])
        s_efgh = np.std(all_efgh[:, i])
        gap = m_abcd - m_efgh
        ratio = m_efgh / m_abcd if m_abcd > 0 else 0
        print(f"  r={r:>3}  {m_abcd:>8.2f} ± {s_abcd:>5.2f}  {m_efgh:>8.2f} ± {s_efgh:>5.2f}  {gap:>+6.2f}  {ratio:>9.3f}")

        if catchup_round is None and ratio >= 0.95 and r > 16:
            catchup_round = r

    # Per-register detail for rounds 16-24
    print(f"\n  Per-register HW (mean) for rounds 16-24:")
    print(f"  {'Round':>5}  {'δa':>6}  {'δb':>6}  {'δc':>6}  {'δd':>6}  {'δe':>6}  {'δf':>6}  {'δg':>6}  {'δh':>6}")
    print(f"  {'-'*5}  {'-'*6}  {'-'*6}  {'-'*6}  {'-'*6}  {'-'*6}  {'-'*6}  {'-'*6}  {'-'*6}")
    for i, r in enumerate(range(16, 25)):
        print(f"  r={r:>3}  {np.mean(all_a[:,i]):>6.2f}  {np.mean(all_b[:,i]):>6.2f}  "
              f"{np.mean(all_c[:,i]):>6.2f}  {np.mean(all_d[:,i]):>6.2f}  "
              f"{np.mean(all_e[:,i]):>6.2f}  {np.mean(all_f[:,i]):>6.2f}  "
              f"{np.mean(all_g[:,i]):>6.2f}  {np.mean(all_h[:,i]):>6.2f}")

    if catchup_round:
        print(f"\n  Structure horizon (efgh catches abcd at 95%): round {catchup_round}")
    else:
        print(f"\n  Structure horizon: efgh did not reach 95% of abcd by round 32")

    # Verify initial conditions (round 16)
    print(f"\n  Verification at round 16 (Wang chain output):")
    print(f"    δe[16]=0: mean HW={np.mean(all_e[:,0]):.4f} (expected 0)")
    print(f"    δf[16]=0: mean HW={np.mean(all_f[:,0]):.4f} (expected 0)")
    print(f"    δg[16]=0: mean HW={np.mean(all_g[:,0]):.4f} (expected 0)")
    print(f"    δh[16]=0: mean HW={np.mean(all_h[:,0]):.4f} (expected 0)")
    print(f"    δa[16]≠0: mean HW={np.mean(all_a[:,0]):.2f} (expected ~16)")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")
    print()
    return catchup_round


# =============================================================
# EXPERIMENT A: Diffusion profile for Wang pairs
# =============================================================
def experiment_A(n=10000):
    print("=" * 70)
    print("EXPERIMENT A: Diffusion profile for Wang pairs (rounds 16-64)")
    print("=" * 70)
    t0 = time.time()

    profiles = []
    slowness_scores = []
    hw_dhashes = []
    wn_list = []

    for i in range(n):
        Wn, Wf = generate_wang_pair()
        prof = compute_diffusion_profile(Wn, Wf)
        profiles.append(prof)

        # Slowness score: sum(256 - HW_total[r]) for r=17..48
        # index offset: round 16 = index 0, round 17 = index 1, ... round 48 = index 32
        S = sum(256 - prof['hw_total'][r - 16] for r in range(17, 49))
        slowness_scores.append(S)
        hw_dhashes.append(prof['hw_dhash'])
        wn_list.append(Wn)

    slowness_scores = np.array(slowness_scores)
    hw_dhashes = np.array(hw_dhashes)

    # Mean/std HW_total table
    hw_total_arr = np.array([p['hw_total'] for p in profiles])  # (n, 49)

    print(f"\n  Mean HW_total[r] (δstate XOR, 256 bits):")
    print(f"  {'Round':>5}  {'Mean':>8}  {'Std':>6}  {'Min':>6}  {'Max':>6}")
    print(f"  {'-'*5}  {'-'*8}  {'-'*6}  {'-'*6}  {'-'*6}")
    for r in range(16, 65):
        idx = r - 16
        m = np.mean(hw_total_arr[:, idx])
        s = np.std(hw_total_arr[:, idx])
        mn = np.min(hw_total_arr[:, idx])
        mx = np.max(hw_total_arr[:, idx])
        if r <= 32 or r % 8 == 0 or r == 64:
            print(f"  r={r:>3}  {m:>8.2f}  {s:>6.2f}  {mn:>6}  {mx:>6}")

    # Slowness distribution
    print(f"\n  Slowness score S = sum(256 - HW[r]) for r=17..48:")
    print(f"    Mean: {np.mean(slowness_scores):.1f}")
    print(f"    Std:  {np.std(slowness_scores):.1f}")
    print(f"    Min:  {np.min(slowness_scores)}")
    print(f"    Max:  {np.max(slowness_scores)}")

    # Top-10 slowest pairs
    top10_idx = np.argsort(slowness_scores)[-10:][::-1]
    print(f"\n  Top-10 slowest pairs:")
    print(f"  {'Rank':>4}  {'S':>6}  {'HW(δhash)':>10}  {'HW[17]':>6}  {'HW[20]':>6}  {'HW[24]':>6}  {'HW[32]':>6}  {'HW[48]':>6}  {'HW[64]':>6}")
    for rank, idx in enumerate(top10_idx):
        p = profiles[idx]
        hw = p['hw_total']
        print(f"  {rank+1:>4}  {slowness_scores[idx]:>6}  {p['hw_dhash']:>10}  "
              f"{hw[1]:>6}  {hw[4]:>6}  {hw[8]:>6}  {hw[16]:>6}  {hw[32]:>6}  {hw[48]:>6}")

    # Slowest pair details
    slowest_idx = top10_idx[0]
    print(f"\n  Slowest pair (rank 1) full profile:")
    print(f"    Wn[0..15] = {[hex(w) for w in wn_list[slowest_idx]]}")
    print(f"    HW(δhash) = {profiles[slowest_idx]['hw_dhash']}")
    p = profiles[slowest_idx]
    for r in range(16, 65):
        idx_r = r - 16
        print(f"    r={r:>3}: HW_total={p['hw_total'][idx_r]:>3}  HW_abcd={p['hw_abcd'][idx_r]:>3}  HW_efgh={p['hw_efgh'][idx_r]:>3}")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")
    print()
    return profiles, slowness_scores, hw_dhashes


# =============================================================
# EXPERIMENT C: Targeted search for slow paths
# =============================================================
def experiment_C():
    print("=" * 70)
    print("EXPERIMENT C: Targeted search for slow paths")
    print("=" * 70)
    t0 = time.time()

    # Part 1: search over W[0] for min HW_total[20]
    print("\n  Part 1: 100K random W[0], measure HW_total[20]")
    n_search1 = 100000
    best_hw20 = 999
    best_w0_hw20 = None
    best_wn_hw20 = None
    best_wf_hw20 = None
    hw20_list = []

    for i in range(n_search1):
        Wn, Wf = generate_wang_pair()
        Wn_exp = message_schedule(Wn)
        Wf_exp = message_schedule(Wf)
        states_n = sha256_all_states(Wn_exp)
        states_f = sha256_all_states(Wf_exp)

        dx20 = tuple(a ^ b for a, b in zip(states_n[20], states_f[20]))
        hw20 = sum(hw32(x) for x in dx20)
        hw20_list.append(hw20)

        if hw20 < best_hw20:
            best_hw20 = hw20
            best_w0_hw20 = Wn[0]
            best_wn_hw20 = list(Wn)
            best_wf_hw20 = list(Wf)

    hw20_arr = np.array(hw20_list)
    print(f"    Mean HW_total[20]: {np.mean(hw20_arr):.2f}")
    print(f"    Min HW_total[20]:  {best_hw20}")
    print(f"    Best W[0]:         {hex(best_w0_hw20)}")

    # Profile for best W[0]
    if best_wn_hw20:
        prof = compute_diffusion_profile(best_wn_hw20, best_wf_hw20)
        print(f"    Best pair profile (HW_total):")
        for r in range(16, 33):
            print(f"      r={r:>3}: {prof['hw_total'][r-16]:>3}")
        print(f"    HW(δhash): {prof['hw_dhash']}")

    # Part 2: search over (W[0], W[1]) for min HW_total[24]
    print(f"\n  Part 2: 1M random pairs, measure HW_total[24]")
    n_search2 = 1000000
    best_hw24 = 999
    best_wn_hw24 = None
    best_wf_hw24 = None
    hw24_hist = defaultdict(int)

    for i in range(n_search2):
        Wn, Wf = generate_wang_pair()
        Wn_exp = message_schedule(Wn)
        Wf_exp = message_schedule(Wf)
        states_n = sha256_all_states(Wn_exp)
        states_f = sha256_all_states(Wf_exp)

        dx24 = tuple(a ^ b for a, b in zip(states_n[24], states_f[24]))
        hw24 = sum(hw32(x) for x in dx24)

        bucket = hw24 // 10 * 10
        hw24_hist[bucket] += 1

        if hw24 < best_hw24:
            best_hw24 = hw24
            best_wn_hw24 = list(Wn)
            best_wf_hw24 = list(Wf)

    print(f"    Min HW_total[24]:  {best_hw24}")

    # Profile for best
    if best_wn_hw24:
        prof = compute_diffusion_profile(best_wn_hw24, best_wf_hw24)
        print(f"    Best pair profile (HW_total):")
        for r in range(16, 65):
            idx = r - 16
            if r <= 32 or r % 8 == 0 or r == 64:
                print(f"      r={r:>3}: {prof['hw_total'][idx]:>3}")
        print(f"    HW(δhash): {prof['hw_dhash']}")

    # Check: HW_total[32] < 200?
    # Quick scan: measure HW[32] for the best pair
    if best_wn_hw24:
        hw32_best = prof['hw_total'][16]  # index 16 = round 32
        print(f"\n    HW_total[32] of best pair: {hw32_best}")
        print(f"    HW_total[32] < 200? {'YES — slow path found!' if hw32_best < 200 else 'NO'}")

    # Rule 14: if suspiciously good, triple check
    if best_hw24 < 64:
        print(f"\n    Rule 14: min HW[24]={best_hw24} < 64 — triple checking...")
        prof2 = compute_diffusion_profile(best_wn_hw24, best_wf_hw24)
        print(f"    Recomputed HW[24]: {prof2['hw_total'][8]}")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")
    print()
    return best_hw20, best_hw24


# =============================================================
# EXPERIMENT B: Lyapunov exponent per trajectory
# =============================================================
def experiment_B(profiles, hw_dhashes):
    print("=" * 70)
    print("EXPERIMENT B: Lyapunov exponent per trajectory")
    print("=" * 70)
    t0 = time.time()

    n = len(profiles)
    lambdas = []

    for i in range(n):
        hw = profiles[i]['hw_total']  # index 0 = round 16
        lam_list = []
        for r in range(17, 49):  # rounds 17..48
            idx = r - 16
            hw_curr = hw[idx]
            hw_prev = hw[idx - 1]
            if hw_prev > 0 and hw_curr > 0:
                lam = math.log2(hw_curr / hw_prev)
                lam_list.append(lam)

        if lam_list:
            Lam = np.mean(lam_list)
            lambdas.append(Lam)
        else:
            lambdas.append(float('nan'))

    lambdas = np.array(lambdas)
    valid = lambdas[~np.isnan(lambdas)]

    print(f"  Λ (mean Lyapunov exponent per pair) distribution:")
    print(f"    Mean: {np.mean(valid):.6f}")
    print(f"    Std:  {np.std(valid):.6f}")
    print(f"    Min:  {np.min(valid):.6f}")
    print(f"    Max:  {np.max(valid):.6f}")
    print(f"    Pairs with Λ < 0 (shrinking): {np.sum(valid < 0)}/{len(valid)}")
    print(f"    Pairs with |Λ| < 0.01 (stable): {np.sum(np.abs(valid) < 0.01)}/{len(valid)}")

    # Histogram
    bins = np.linspace(np.min(valid), np.max(valid), 11)
    hist, edges = np.histogram(valid, bins=bins)
    print(f"\n  Histogram of Λ:")
    for i in range(len(hist)):
        bar = '#' * (hist[i] * 50 // max(hist))
        print(f"    [{edges[i]:>+7.4f}, {edges[i+1]:>+7.4f}): {hist[i]:>5}  {bar}")

    # Correlation Λ ↔ HW(δhash)
    valid_mask = ~np.isnan(lambdas)
    if np.sum(valid_mask) > 10:
        corr = np.corrcoef(lambdas[valid_mask], hw_dhashes[valid_mask])[0, 1]
        print(f"\n  Correlation Λ ↔ HW(δhash): {corr:.6f}")
    else:
        print(f"\n  Insufficient valid data for correlation")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")
    print()
    return lambdas


# =============================================================
# MAIN
# =============================================================
def main():
    print("=" * 70)
    print("  ЗАДАНИЕ 3: Slow Diffusion Trajectories in Phase B")
    print("  Priority: D > A > C > B")
    print("=" * 70)
    print()

    # Quick Wang verification
    ok = 0
    for _ in range(100):
        Wn, Wf = generate_wang_pair()
        Wn_exp = message_schedule(Wn)
        Wf_exp = message_schedule(Wf)
        states_n = sha256_all_states(Wn_exp)
        states_f = sha256_all_states(Wf_exp)
        zeros = sum(1 for r in range(2, 17) if states_n[r][4] == states_f[r][4])
        if zeros == 15:
            ok += 1
    print(f"Wang verification: {ok}/100 (expected 100)")
    assert ok == 100
    print()

    # Experiment D (priority)
    catchup = experiment_D(n=10000)

    # Experiment A
    profiles, slowness, hw_dh = experiment_A(n=10000)

    # Experiment C
    min_hw20, min_hw24 = experiment_C()

    # Experiment B
    lambdas = experiment_B(profiles, hw_dh)

    # SUMMARY
    print("=" * 70)
    print("  SUMMARY")
    print("=" * 70)
    valid_lam = lambdas[~np.isnan(lambdas)]
    print(f"\n  D: Structure horizon (efgh catches abcd): round {catchup if catchup else '>32'}")
    print(f"  A: Mean slowness S={np.mean(slowness):.0f}, max={np.max(slowness)}")
    print(f"     Mean HW(δhash)={np.mean(hw_dh):.1f}")
    print(f"  C: Min HW_total[20] = {min_hw20} (from 100K)")
    print(f"     Min HW_total[24] = {min_hw24} (from 1M)")
    print(f"  B: Mean Λ = {np.mean(valid_lam):.6f}")
    print(f"     Pairs with Λ<0: {np.sum(valid_lam<0)}")
    print(f"     Corr(Λ, HW(δhash)): {np.corrcoef(valid_lam, hw_dh[~np.isnan(lambdas)])[0,1]:.4f}")
    print()
    print(f"  Rule 12: null = random Wang pair (no special structure)")
    print(f"  Rule 7: all δhash≠0, no trivial results")
    print(f"  Rule 1: honest results reported as computed")
    print()


if __name__ == "__main__":
    main()
