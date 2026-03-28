#!/usr/bin/env python3
"""
Ch-Bridge Exploitation — Stage 2.6
=====================================

Stage 2.5 discovered: Ch(e62,f62,g62) is THE bridge function.
corr(Ch62, H7) = -0.201 — stronger than raw63 (-0.097).

The algebraic reason: g[62] = H[7] - IV[7]. Ch contains g directly.

Stage 2.6 exploits this:

  A1: INVERSE BRIDGE — given H, compute Ch62 WITHOUT knowing input.
      Ch62 = Ch(H[5]-IV[5], H[6]-IV[6], H[7]-IV[7])
      This is a PASSIVE distinguisher (no chosen-prefix needed)!

  A2: Ch62 INVARIANT VERIFICATION — T_CH_INVARIANT says Ch62[b31,b30]=0
      when carry[63]=0. Can we detect this from output alone?

  A3: MULTI-WORD DISTINGUISHER — use Ch62 + correlations H5,H6,H7
      in the tail to build a stronger passive distinguisher.

  A4: THEORETICAL LIMIT — what is the maximum possible distinguisher
      advantage using only the Ch-bridge?
"""

import numpy as np
from math import erfc, sqrt
import time
import sys

MASK = 0xFFFFFFFF
T_VAL = float(1 << 32)

H0 = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
]

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
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
]

def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK
def sig0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sig1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)
def Sig0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sig1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
def hw(x): return bin(x & MASK).count('1')

def message_schedule(W16):
    W = list(W16) + [0] * 48
    for i in range(16, 64):
        W[i] = (sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16]) & MASK
    return W

def full_sha(W16):
    W = message_schedule(W16)
    a, b, c, d, e, f, g, h = H0
    raws = []
    for r in range(64):
        raw = h + Sig1(e) + Ch(e, f, g) + K[r] + W[r]
        raws.append(raw)
        T1 = raw & MASK
        T2 = (Sig0(a) + Maj(a, b, c)) & MASK
        h, g, f = g, f, e
        e = (d + T1) & MASK
        d, c, b = c, b, a
        a = (T1 + T2) & MASK
    H_out = [(v + iv) & MASK for v, iv in zip([a,b,c,d,e,f,g,h], H0)]
    return raws, H_out


def ch_from_output(H_out):
    """Compute Ch(e62, f62, g62) from output H alone.

    State at end of round 63: (a63, b63, c63, d63, e63, f63, g63, h63)
    H[i] = state63[i] + IV[i]

    Register shifts: f63=e62, g63=f62=e61, h63=g62=f61=e60
    So: e62 = f63 = H[5] - IV[5]
        f62 = g63 = H[6] - IV[6]   (actually f62=e61, g63=f62)
        g62 = h63 = H[7] - IV[7]   (actually g62=f61=e60, h63=g62)

    Wait — let me be precise about the shift register:
    After round r, the state update does:
      h_new = g_old, g_new = f_old, f_new = e_old
    So after round 63:
      h63 = g62, g63 = f62, f63 = e62

    H[7] = h63 + IV[7] → h63 = H[7] - IV[7] = g62
    H[6] = g63 + IV[6] → g63 = H[6] - IV[6] = f62
    H[5] = f63 + IV[5] → f63 = H[5] - IV[5] = e62

    Ch(e62, f62, g62) = Ch(H[5]-IV[5], H[6]-IV[6], H[7]-IV[7])
    """
    e62 = (H_out[5] - H0[5]) & MASK
    f62 = (H_out[6] - H0[6]) & MASK
    g62 = (H_out[7] - H0[7]) & MASK
    return Ch(e62, f62, g62)


def sig1_from_output(H_out):
    """Compute Sig1(e62) from output H alone."""
    e62 = (H_out[5] - H0[5]) & MASK
    return Sig1(e62)


# ============================================================
# A1: Inverse Bridge — Passive Distinguisher
# ============================================================

def experiment_A1(N=20000, seed=800):
    print("=" * 70)
    print("A1: INVERSE BRIDGE — Compute Ch62 from output H alone")
    print(f"N={N}, seed={seed}")
    print("=" * 70)

    rng = np.random.RandomState(seed)

    # Verify: Ch62 computed from output matches Ch62 computed internally
    print(f"\n--- Verification: Ch_from_H == Ch_internal ---")
    n_verify = 1000
    match = 0
    for i in range(n_verify):
        W16 = [int(rng.randint(0, 1 << 32)) for _ in range(16)]
        # Internal computation
        W = message_schedule(W16)
        a, b, c, d, e, f, g, h = H0
        for r in range(64):
            raw = h + Sig1(e) + Ch(e, f, g) + K[r] + W[r]
            T1 = raw & MASK
            T2 = (Sig0(a) + Maj(a, b, c)) & MASK
            if r == 61:
                ch62_internal = Ch(e, f, g)  # Ch computed at START of round 62
                # Wait — Ch at round 62 uses state BEFORE round 62 update
                # Actually: Ch in round r uses current e,f,g = state after round r-1
                # So Ch in round 62 = Ch(e[62], f[62], g[62]) where state[62] is after 62 rounds
            h, g, f = g, f, e
            e = (d + T1) & MASK
            d, c, b = c, b, a
            a = (T1 + T2) & MASK

        H_out = [(v + iv) & MASK for v, iv in zip([a,b,c,d,e,f,g,h], H0)]

        # Note: round 62 uses state[62] which is state AFTER 62 rounds
        # The h,g,f,e in the loop at r=62 are state[62] values
        # But we reassign h,g,f = g,f,e AFTER computing raw
        # So at r=62: e,f,g used in Ch ARE state[62].e, state[62].f, state[62].g
        # After round 63: f63 = e62, g63 = f62, h63 = g62
        # H[5] = f63 + IV[5], H[6] = g63 + IV[6], H[7] = h63 + IV[7]

        ch62_from_H = ch_from_output(H_out)

        # Internal: we need Ch at the START of round 62 = Ch(state[62].e, .f, .g)
        # But our loop computes Ch at start of each round using current e,f,g
        # At r=62, current e,f,g = state after round 61 = state[62] in our indexing
        # Hmm, let me recompute properly

    # Redo with explicit state tracking
    match = 0
    for i in range(n_verify):
        W16 = [int(rng.randint(0, 1 << 32)) for _ in range(16)]
        W = message_schedule(W16)
        a, b, c, d, e, f, g, h = H0
        states = [(a,b,c,d,e,f,g,h)]
        for r in range(64):
            raw = h + Sig1(e) + Ch(e, f, g) + K[r] + W[r]
            T1 = raw & MASK
            T2 = (Sig0(a) + Maj(a, b, c)) & MASK
            h, g, f = g, f, e
            e = (d + T1) & MASK
            d, c, b = c, b, a
            a = (T1 + T2) & MASK
            states.append((a,b,c,d,e,f,g,h))

        H_out = [(states[64][i] + H0[i]) & MASK for i in range(8)]

        # State BEFORE round 63 (= state after round 62 = states[62])
        # Ch in round 63 uses states[63].e, .f, .g
        # Wait: in the loop, at r=63, the e,f,g at START of iteration
        # are states[63] (result of 63 rounds, 0-indexed as states[63])
        # But SHA indexes: state[0]=IV, state[r+1]=after round r
        # So states[63] = after 62 rounds = before round 63
        # Ch used in round 63 = Ch(states[63][4], states[63][5], states[63][6])

        # Actually let's think simpler:
        # After all 64 rounds: final state = states[64]
        # H[5] = states[64][5] + IV[5] = f_final + IV[5]
        # f_final = states[64][5]
        # In the round loop: after round 63, f was set to e (from round 62's state)
        # states[64][5] = f after round 63 = e before round 63 = states[63][4]
        # So H[5] - IV[5] = states[63][4] = e before round 63

        # Similarly: g_final = states[64][6] = f before round 63 = states[63][5]
        # h_final = states[64][7] = g before round 63 = states[63][6]

        e_before_63 = states[63][4]
        f_before_63 = states[63][5]
        g_before_63 = states[63][6]
        ch63_internal = Ch(e_before_63, f_before_63, g_before_63)

        ch63_from_H = ch_from_output(H_out)

        if ch63_internal == ch63_from_H:
            match += 1

    print(f"  Ch63 match (internal vs from-H): {match}/{n_verify} ({match/n_verify*100:.1f}%)")

    if match < n_verify * 0.99:
        # Try Ch62 instead
        match62 = 0
        for i in range(n_verify):
            W16 = [int(rng.randint(0, 1 << 32)) for _ in range(16)]
            W = message_schedule(W16)
            a, b, c, d, e, f, g, h = H0
            states = [(a,b,c,d,e,f,g,h)]
            for r in range(64):
                raw = h + Sig1(e) + Ch(e, f, g) + K[r] + W[r]
                T1 = raw & MASK
                T2 = (Sig0(a) + Maj(a, b, c)) & MASK
                h, g, f = g, f, e
                e = (d + T1) & MASK
                d, c, b = c, b, a
                a = (T1 + T2) & MASK
                states.append((a,b,c,d,e,f,g,h))

            H_out = [(states[64][i] + H0[i]) & MASK for i in range(8)]

            # Ch used in round 62: Ch(states[62][4], states[62][5], states[62][6])
            ch62_int = Ch(states[62][4], states[62][5], states[62][6])

            # From H: e_before_62 = states[62][4]
            # states[63][5] = f after round 62 = e before round 62 = states[62][4]
            # states[64][6] = g after round 63 = f after round 63 = ...
            # This gets complicated. Let's just use the formula and check.

            # The key: Ch_from_output computes Ch(e62_fromH, f62_fromH, g62_fromH)
            # where e62_fromH = H[5]-IV[5], f62_fromH = H[6]-IV[6], g62_fromH = H[7]-IV[7]
            # These correspond to states[63][4], states[63][5], states[63][6]
            # So ch_from_output gives Ch for round 63, not round 62!

            ch63_int = Ch(states[63][4], states[63][5], states[63][6])
            ch63_fromH = ch_from_output(H_out)
            if ch63_int == ch63_fromH:
                match62 += 1

        print(f"  Ch (round 63 input) match: {match62}/{n_verify} ({match62/n_verify*100:.1f}%)")

    # Now the actual passive distinguisher
    print(f"\n--- PASSIVE DISTINGUISHER: Ch from output only ---")

    ch_vals = []
    ch_hw = []
    h7_b29 = []
    h7_b31 = []
    h6_b29 = []
    h6_b31 = []
    carry63 = []

    t0 = time.time()
    for i in range(N):
        W16 = [int(rng.randint(0, 1 << 32)) for _ in range(16)]
        raws, H_out = full_sha(W16)

        ch = ch_from_output(H_out)
        ch_vals.append(ch)
        ch_hw.append(hw(ch))
        h7_b29.append((H_out[7] >> 29) & 1)
        h7_b31.append((H_out[7] >> 31) & 1)
        h6_b29.append((H_out[6] >> 29) & 1)
        h6_b31.append((H_out[6] >> 31) & 1)
        carry63.append(1 if raws[63] >= T_VAL else 0)

        if (i+1) % 10000 == 0:
            print(f"  {i+1}/{N} ({time.time()-t0:.1f}s)")

    ch_arr = np.array(ch_vals, dtype=float)
    ch_hw_arr = np.array(ch_hw, dtype=float)

    # T_CH_INVARIANT: when carry[63]=0, Ch[b31]=Ch[b30]=0
    # Can we detect carry[63]=0 from Ch alone?
    print(f"\n  Ch distribution:")
    print(f"    E[HW(Ch)]: {np.mean(ch_hw_arr):.2f}")
    print(f"    P(Ch[b31]=0): {np.mean([1 - ((c >> 31) & 1) for c in ch_vals]):.4f}")
    print(f"    P(Ch[b30]=0): {np.mean([1 - ((c >> 30) & 1) for c in ch_vals]):.4f}")

    # Correlation Ch → H bits
    print(f"\n  Correlation Ch → H[7] bits (PASSIVE, no input knowledge):")
    for bit in [28, 29, 30, 31]:
        h7_bit = np.array([(H >> bit) & 1 for H in [ch_vals[i] for i in range(N)]], dtype=float)
        # Wait, we want corr(Ch, H7[bit])
        h7_bit_arr = np.array(h7_b29 if bit == 29 else (h7_b31 if bit == 31 else
                              [((ch_vals[i]) >> bit) & 1 for i in range(N)]), dtype=float)
        # Actually let me be explicit
    for bit in [28, 29, 30, 31]:
        target = np.array([((H_out_dummy >> bit) & 1) for H_out_dummy in [0]*N], dtype=float)
        # Recompute properly
        pass

    # Recompute H values
    all_H7 = np.array([0]*N, dtype=float)
    all_H6 = np.array([0]*N, dtype=float)
    # Already have h7_b29, h7_b31 etc as lists

    print(f"\n  corr(Ch_value, H7_bit):")
    for name, arr in [('H7[b29]', h7_b29), ('H7[b31]', h7_b31), ('H6[b29]', h6_b29), ('H6[b31]', h6_b31)]:
        target = np.array(arr, dtype=float)
        c = np.corrcoef(ch_arr, target)[0, 1]
        c_hw = np.corrcoef(ch_hw_arr, target)[0, 1]
        print(f"    corr(Ch, {name}): {c:+.5f},  corr(HW(Ch), {name}): {c_hw:+.5f}")

    # Per-bit Ch → H7 correlation
    print(f"\n  Per-bit: corr(Ch[bit_i], H7[b31]):")
    h7_b31_arr = np.array(h7_b31, dtype=float)
    best_ch_bit = -1
    best_ch_corr = 0
    for bit in range(32):
        ch_bit = np.array([(c >> bit) & 1 for c in ch_vals], dtype=float)
        if np.std(ch_bit) > 0.01:
            c = np.corrcoef(ch_bit, h7_b31_arr)[0, 1]
            if abs(c) > abs(best_ch_corr):
                best_ch_corr = c
                best_ch_bit = bit
            if abs(c) > 0.02:
                print(f"    Ch[b{bit:2d}] → H7[b31]: {c:+.4f}")

    print(f"    Best: Ch[b{best_ch_bit}] → H7[b31]: {best_ch_corr:+.4f}")

    # KEY: Ch bit 31 as PASSIVE filter
    print(f"\n--- KEY: Ch[b31]=0 as passive carry-indicator ---")
    ch_b31_0 = [i for i in range(N) if not ((ch_vals[i] >> 31) & 1)]
    ch_b31_1 = [i for i in range(N) if ((ch_vals[i] >> 31) & 1)]

    p_h7b31_given_ch0 = np.mean([h7_b31[i] for i in ch_b31_0])
    p_h7b31_given_ch1 = np.mean([h7_b31[i] for i in ch_b31_1])
    p_h7b29_given_ch0 = np.mean([h7_b29[i] for i in ch_b31_0])
    p_h7b29_given_ch1 = np.mean([h7_b29[i] for i in ch_b31_1])

    print(f"  N(Ch[b31]=0): {len(ch_b31_0)} ({len(ch_b31_0)/N*100:.1f}%)")
    print(f"  N(Ch[b31]=1): {len(ch_b31_1)} ({len(ch_b31_1)/N*100:.1f}%)")
    print(f"  P(H7[b31]=1 | Ch[b31]=0): {p_h7b31_given_ch0:.4f}")
    print(f"  P(H7[b31]=1 | Ch[b31]=1): {p_h7b31_given_ch1:.4f}")
    print(f"  Phi H7[b31]: {p_h7b31_given_ch0 - p_h7b31_given_ch1:+.4f}")
    print(f"  P(H7[b29]=1 | Ch[b31]=0): {p_h7b29_given_ch0:.4f}")
    print(f"  P(H7[b29]=1 | Ch[b31]=1): {p_h7b29_given_ch1:.4f}")
    print(f"  Phi H7[b29]: {p_h7b29_given_ch0 - p_h7b29_given_ch1:+.4f}")

    # ALSO: Ch[b31]=0 AND Ch[b30]=0 (double filter)
    ch_top2_0 = [i for i in range(N) if not ((ch_vals[i] >> 31) & 1) and not ((ch_vals[i] >> 30) & 1)]
    print(f"\n  Ch[b31,b30]=0: N={len(ch_top2_0)} ({len(ch_top2_0)/N*100:.1f}%)")
    if ch_top2_0:
        p31 = np.mean([h7_b31[i] for i in ch_top2_0])
        p29 = np.mean([h7_b29[i] for i in ch_top2_0])
        print(f"    P(H7[b31]=1): {p31:.4f} (baseline: {np.mean(h7_b31):.4f}, phi={p31-np.mean(h7_b31):+.4f})")
        print(f"    P(H7[b29]=1): {p29:.4f} (baseline: {np.mean(h7_b29):.4f}, phi={p29-np.mean(h7_b29):+.4f})")

    # PASSIVE DISTINGUISHER FORMULA
    print(f"\n  ★ PASSIVE DISTINGUISHER (no chosen-prefix!):")
    print(f"    Given only H[5,6,7]:")
    print(f"    1. e62 = H[5] - IV[5]")
    print(f"    2. f62 = H[6] - IV[6]")
    print(f"    3. g62 = H[7] - IV[7]")
    print(f"    4. Ch62 = (e62 & f62) ^ (~e62 & g62)")
    print(f"    5. If Ch62[b31]=0 AND Ch62[b30]=0 → lean 'SHA-256'")
    print(f"    Cost: 3 subtractions + 1 Ch computation = O(1)")

    # Compare with random oracle
    print(f"\n  For RANDOM oracle: H[5,6,7] uniform → e62,f62,g62 uniform")
    print(f"  P(Ch[b31]=0 | random) = 0.5000 exactly")
    print(f"  P(Ch[b31]=0 | SHA-256) = {np.mean([1-((c>>31)&1) for c in ch_vals]):.4f}")
    print(f"  These should be EQUAL for random inputs (no HC).")
    print(f"  Difference would mean SHA-256 ≠ random oracle PASSIVELY.")

    return ch_vals, ch_hw


# ============================================================
# MAIN
# ============================================================

if __name__ == '__main__':
    print("Ch-Bridge Exploitation — Stage 2.6")
    print(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print()

    N = 20000
    if len(sys.argv) > 1 and sys.argv[1] == '--fast':
        N = 10000

    t_start = time.time()
    ch_vals, ch_hw = experiment_A1(N=N)

    total = time.time() - t_start
    print(f"\n{'='*70}")
    print(f"TOTAL TIME: {total:.1f}s")
    print(f"{'='*70}")
