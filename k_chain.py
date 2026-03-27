#!/usr/bin/env python3
"""
K-Chain — Stage 1.4: Algebraic invariants through K[r] constants
=================================================================

Three paths around the chaotic zone (r=31..59):

PATH 1: K-CHAIN — deterministic invariants via K[r] magnitude
  Each K[r] is large enough to force certain bits of intermediate
  registers to 0 when carry=0. This creates a CHAIN of constraints
  from r=63 backward. How far back does it reach?

PATH 2: MULTI-BLOCK — carry structure of block 1 affects IV of block 2
  If carry[63]=0 in block 1, IV2 is biased. Does block 2 inherit structure?

PATH 3: BACKWARD — from H to state[62], skipping chaotic zone entirely
  H[7] = g[62] + IV[7]. If we know H[7], we know g[62].
  Can we chain backwards from g[62] to find structure in earlier rounds?
"""

import numpy as np
from collections import defaultdict, Counter
import time
import sys

MASK = 0xFFFFFFFF

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

def sha256_rounds(W16, num_rounds=64):
    """Returns list of states (a,b,c,d,e,f,g,h) for rounds 0..num_rounds."""
    W = message_schedule(W16)
    a, b, c, d, e, f, g, h = H0
    states = [(a, b, c, d, e, f, g, h)]
    carries = []
    raws = []

    for r in range(min(num_rounds, 64)):
        raw = h + Sig1(e) + Ch(e, f, g) + K[r] + W[r]
        T1 = raw & MASK
        carries.append(1 if raw >= (1 << 32) else 0)
        raws.append(raw)
        T2 = (Sig0(a) + Maj(a, b, c)) & MASK
        h, g, f = g, f, e
        e = (d + T1) & MASK
        d, c, b = c, b, a
        a = (T1 + T2) & MASK
        states.append((a, b, c, d, e, f, g, h))

    H_out = [(states[-1][i] + H0[i]) & MASK for i in range(8)]
    return states, carries, raws, W, H_out


def hill_climb_ss63(rng, steps=200):
    """Find W[0] with small SS[63] via hill climbing. Returns (W16, states, carries)."""
    W0 = int(rng.randint(0, 1 << 32))
    W16 = [W0] + [0] * 15
    states, carries, raws, W, H = sha256_rounds(W16)
    best_ss = raws[63] if len(raws) > 63 else float('inf')

    for _ in range(steps):
        b = int(rng.randint(0, 32))
        W16_try = [W0 ^ (1 << b)] + [0] * 15
        st2, ca2, ra2, W2, H2 = sha256_rounds(W16_try)
        ss2 = ra2[63] if len(ra2) > 63 else float('inf')
        if ss2 < best_ss:
            W0 = W0 ^ (1 << b)
            W16 = W16_try
            states, carries, raws, W, H = st2, ca2, ra2, W2, H2
            best_ss = ss2

    return W16, states, carries, raws, W, H, best_ss


# ============================================================
# PATH 1: K-CHAIN — Backward invariant chain through K[r]
# ============================================================

def path1_k_chain(N=3000, seed=70):
    """
    For each round r=63..0, determine:
    1. What is the maximum SS[r] compatible with carry[r]=0?
       max_SS = 2^32 - K[r] - 1 (since SS + K + W < 2^32, and W >= 0)
       Actually: SS[r] + K[r] + W[r] < 2^32, so SS[r] < 2^32 - K[r] (when W[r]=0)
    2. Which bits of Ch[r-1], g[r-1], e[r-1] are forced to 0?
    3. Does carry[r]=0 force carry[r-1]=0 with any probability?
    """
    print("=" * 70)
    print("PATH 1: K-CHAIN — Backward Invariant Chain")
    print(f"N={N}, seed={seed}")
    print("=" * 70)

    # Analytical: max SS[r] when carry[r]=0 and W[r]=0
    print(f"\n--- Analytical: Max SS[r] when carry[r]=0 ---")
    print(f"  carry[r]=0 ⟺ h[r-1] + Sig1(e[r-1]) + Ch(e,f,g)[r-1] + K[r] + W[r] < 2^32")
    print(f"  With W[r]=0: SS[r] < 2^32 - K[r]")
    print(f"  SS[r] = h[r-1] + Sig1(e[r-1]) + Ch(e,f,g)[r-1]")
    print()

    print(f"  {'r':>3} | {'K[r]':>12} | {'K[r]/2^32':>10} | {'max_SS':>12} | {'max_SS/2^32':>12} | {'det bits':>9}")
    print(f"  {'-'*3}-+-{'-'*12}-+-{'-'*10}-+-{'-'*12}-+-{'-'*12}-+-{'-'*9}")

    k_chain_info = []
    for r in range(64):
        max_ss = ((1 << 32) - K[r]) & MASK  # When W[r]=0
        k_ratio = K[r] / (1 << 32)
        ss_ratio = max_ss / (1 << 32)

        # How many top bits of SS are forced to 0?
        # If max_ss < 2^31: bit 31 = 0 (deterministic)
        # If max_ss < 2^30: bits 31,30 = 0
        det_bits = 0
        for b in range(31, -1, -1):
            if max_ss < (1 << (b + 1)):
                det_bits += 1
            else:
                break

        k_chain_info.append((r, K[r], max_ss, det_bits))

        marker = " ★" if det_bits >= 1 else ""
        if det_bits >= 1 or r >= 58 or r <= 5:
            print(f"  {r:3d} | 0x{K[r]:08x} | {k_ratio:10.4f} | 0x{max_ss:08x} | {ss_ratio:10.4f} | {det_bits:9d}{marker}")

    # Rounds with deterministic bits
    det_rounds = [(r, K_r, max_ss, db) for r, K_r, max_ss, db in k_chain_info if db >= 1]
    print(f"\n  Rounds with deterministic bit constraints (carry=0 → SS top bits = 0):")
    print(f"  Total: {len(det_rounds)} / 64")
    for r, K_r, max_ss, db in det_rounds:
        bits_str = ','.join([f'b{31-i}' for i in range(db)])
        print(f"    r={r:2d}: K=0x{K_r:08x}, {db} det bits ({bits_str}), max_SS=0x{max_ss:08x}")

    # Can we chain? If carry[r]=0 forces SS[r] small, does that force carry[r-1]=0?
    print(f"\n--- Backward Chain Test ---")
    print(f"  Question: carry[r]=0 → P(carry[r-1]=0) > baseline?")

    rng = np.random.RandomState(seed)

    # Use HC to find carry=0 events, then check backward chain
    n_c0_found = 0
    backward_chain = defaultdict(lambda: [0, 0])  # [count_c0, total]

    t0 = time.time()
    for i in range(N):
        W16, states, carries, raws, W, H, best_ss = hill_climb_ss63(rng, steps=150)

        if carries[63] == 0:
            n_c0_found += 1
            # Check backward: which previous rounds also have carry=0?
            for r in range(63):
                backward_chain[r][1] += 1
                if carries[r] == 0:
                    backward_chain[r][0] += 1

        if (i + 1) % 1000 == 0:
            print(f"  {i+1}/{N}, carry[63]=0: {n_c0_found} ({time.time()-t0:.1f}s)")

    print(f"\n  carry[63]=0 found: {n_c0_found} / {N}")

    if n_c0_found > 10:
        # Baseline P(carry[r]=0) from all samples
        all_c0 = np.zeros(64)
        for i in range(min(N, 2000)):
            W16 = [int(rng.randint(0, 1 << 32)) for _ in range(16)]
            _, carries_all, _, _, _ = sha256_rounds(W16)
            for r in range(64):
                all_c0[r] += (1 - carries_all[r])
        all_c0 /= min(N, 2000)

        print(f"\n  {'r':>3} | {'P(c[r]=0|c63=0)':>16} | {'P(c[r]=0) base':>15} | {'ratio':>6} | {'chain?':>7}")
        print(f"  {'-'*3}-+-{'-'*16}-+-{'-'*15}-+-{'-'*6}-+-{'-'*7}")

        chain_rounds = []
        for r in range(64):
            if backward_chain[r][1] > 0:
                p_cond = backward_chain[r][0] / backward_chain[r][1]
                p_base = all_c0[r]
                ratio = p_cond / max(p_base, 1e-10)
                is_chain = ratio > 2.0 and backward_chain[r][0] >= 3
                if is_chain:
                    chain_rounds.append(r)
                marker = " ★★★" if ratio > 5 else (" ★★" if ratio > 3 else (" ★" if ratio > 2 else ""))
                if p_base > 0.001 or p_cond > 0.01 or r >= 55:
                    print(f"  {r:3d} | {p_cond:16.4f} | {p_base:15.4f} | {ratio:6.1f} | {'YES' if is_chain else 'no':>7}{marker}")

        print(f"\n  Chain rounds (ratio > 2×): {chain_rounds}")
        print(f"  Chain length from r=63: {len([r for r in chain_rounds if r >= 55])}")

    return k_chain_info, det_rounds


# ============================================================
# PATH 2: MULTI-BLOCK — Block 1 carry structure → Block 2
# ============================================================

def path2_multiblock(N=2000, seed=71):
    print("\n" + "=" * 70)
    print("PATH 2: MULTI-BLOCK — Carry structure across block boundary")
    print(f"N={N}, seed={seed}")
    print("=" * 70)

    rng = np.random.RandomState(seed)

    # Block 1: HC to find carry[63]=0, compute H1 = new IV for block 2
    # Block 2: random message, compute with IV = H1
    # Question: does carry[63]=0 in block 1 affect block 2 carry structure?

    bias_b2_carry = np.zeros(64)  # P(carry[r]=0 in block 2 | c63_b1=0)
    base_b2_carry = np.zeros(64)  # P(carry[r]=0 in block 2 | random)
    n_c0 = 0
    n_base = 0

    # Also measure H2[7] bit 29 bias
    h2_b29_c0 = []
    h2_b29_base = []

    t0 = time.time()
    for i in range(N):
        # Block 1 with HC
        W16_b1, states1, carries1, raws1, W1, H1, ss = hill_climb_ss63(rng, steps=150)

        # Block 2: random message, IV = H1
        W16_b2 = [int(rng.randint(0, 1 << 32)) for _ in range(16)]
        W_b2 = message_schedule(W16_b2)

        # Run SHA-256 block 2 with IV = H1
        a, b, c, d, e, f, g, h = H1
        carries2 = []
        for r in range(64):
            raw = h + Sig1(e) + Ch(e, f, g) + K[r] + W_b2[r]
            T1 = raw & MASK
            carries2.append(1 if raw >= (1 << 32) else 0)
            T2 = (Sig0(a) + Maj(a, b, c)) & MASK
            h, g, f = g, f, e
            e = (d + T1) & MASK
            d, c, b = c, b, a
            a = (T1 + T2) & MASK
        H2 = [(v + H1[i]) & MASK for i, v in enumerate([a,b,c,d,e,f,g,h])]

        if carries1[63] == 0:
            n_c0 += 1
            for r in range(64):
                bias_b2_carry[r] += (1 - carries2[r])
            h2_b29_c0.append((H2[7] >> 29) & 1)
        else:
            n_base += 1
            for r in range(64):
                base_b2_carry[r] += (1 - carries2[r])
            h2_b29_base.append((H2[7] >> 29) & 1)

        if (i + 1) % 500 == 0:
            print(f"  {i+1}/{N}, c63=0: {n_c0} ({time.time()-t0:.1f}s)")

    print(f"\n  Block 1 carry[63]=0: {n_c0} / {N}")

    if n_c0 > 5 and n_base > 5:
        bias_b2_carry /= max(n_c0, 1)
        base_b2_carry /= max(n_base, 1)

        print(f"\n--- Block 2 carry profile: conditioned on block 1 carry[63]=0 ---")
        print(f"  {'r':>3} | {'P(c2=0|c1_63=0)':>16} | {'P(c2=0|random)':>15} | {'delta':>8}")
        print(f"  {'-'*3}-+-{'-'*16}-+-{'-'*15}-+-{'-'*8}")

        significant = []
        for r in range(64):
            delta = bias_b2_carry[r] - base_b2_carry[r]
            marker = " ***" if abs(delta) > 0.02 else ""
            if abs(delta) > 0.005 or r in [0, 1, 62, 63]:
                print(f"  {r:3d} | {bias_b2_carry[r]:16.4f} | {base_b2_carry[r]:15.4f} | {delta:+8.4f}{marker}")
            if abs(delta) > 0.02:
                significant.append((r, delta))

        print(f"\n  Significant rounds (|delta|>0.02): {len(significant)}")

        # H2[7] bit 29 bias
        if h2_b29_c0:
            p_c0 = np.mean(h2_b29_c0)
            p_base = np.mean(h2_b29_base) if h2_b29_base else 0.5
            print(f"\n  H2[7][b29] | c1_63=0: P={p_c0:.4f} (N={len(h2_b29_c0)})")
            print(f"  H2[7][b29] | random:  P={p_base:.4f} (N={len(h2_b29_base)})")
            print(f"  Delta: {p_c0 - p_base:+.4f}")

    return n_c0


# ============================================================
# PATH 3: BACKWARD — From H[7] to state[62]
# ============================================================

def path3_backward(N=5000, seed=72):
    print("\n" + "=" * 70)
    print("PATH 3: BACKWARD — From H[7] to internal state")
    print(f"N={N}, seed={seed}")
    print("=" * 70)

    rng = np.random.RandomState(seed)

    # H[7] = g[62] + IV[7]
    # H[6] = f[62] + IV[6] = e[61] + IV[6]
    # H[5] = e[62] + IV[5]
    # So knowing H gives us direct access to state at round 62

    # Question: are there correlations between state[62] components
    # that are NOT present in random state?

    print(f"\n--- Analytical: H → state[62] ---")
    print(f"  H[7] = h[63] + IV[7], h[63] = g[62]  →  g[62] = H[7] - IV[7]")
    print(f"  H[6] = g[63] + IV[6], g[63] = f[62]  →  f[62] = H[6] - IV[6]")
    print(f"  H[5] = f[63] + IV[5], f[63] = e[62]  →  e[62] = H[5] - IV[5]")
    print(f"  H[4] = e[63] + IV[4]                  →  e[63] = H[4] - IV[4]")
    print(f"  e[63] = d[62] + T1[63]")

    # Collect state[62] and check correlations
    g62_vals = []
    f62_vals = []
    e62_vals = []
    ch62_vals = []  # Ch(e[62], f[62], g[62])
    ss63_vals = []
    carry63_vals = []

    t0 = time.time()
    for i in range(N):
        W16 = [int(rng.randint(0, 1 << 32)) for _ in range(16)]
        states, carries, raws, W, H = sha256_rounds(W16)

        # state[62] = states[62]
        a62, b62, c62, d62, e62, f62, g62, h62 = states[62]

        g62_vals.append(g62)
        f62_vals.append(f62)
        e62_vals.append(e62)
        ch62_vals.append(Ch(e62, f62, g62))
        ss63_vals.append(h62 + Sig1(e62) + Ch(e62, f62, g62))
        carry63_vals.append(carries[63] if len(carries) > 63 else 1)

    g62 = np.array(g62_vals, dtype=float)
    f62 = np.array(f62_vals, dtype=float)
    e62 = np.array(e62_vals, dtype=float)
    ch62 = np.array(ch62_vals, dtype=float)
    ss63 = np.array(ss63_vals, dtype=float)
    c63 = np.array(carry63_vals, dtype=float)

    # Correlations between state[62] components
    print(f"\n--- Correlations between state[62] components ---")
    pairs = [('g62', g62), ('f62', f62), ('e62', e62), ('Ch62', ch62), ('SS63', ss63)]
    print(f"  {'':>6}", end='')
    for name, _ in pairs:
        print(f" | {name:>8}", end='')
    print()
    for name_i, arr_i in pairs:
        print(f"  {name_i:>6}", end='')
        for name_j, arr_j in pairs:
            if np.std(arr_i) > 0 and np.std(arr_j) > 0:
                c = np.corrcoef(arr_i, arr_j)[0, 1]
            else:
                c = 0
            print(f" | {c:+8.4f}", end='')
        print()

    # Key test: do any state[62] components predict carry[63]?
    print(f"\n--- state[62] → carry[63] prediction ---")
    for name, arr in pairs:
        if np.std(arr) > 0 and np.std(c63) > 0:
            c = np.corrcoef(arr, c63)[0, 1]
            print(f"  corr({name}, carry[63]) = {c:+.5f}")

    # Bit-level: which bits of g[62] predict carry[63]?
    print(f"\n--- g[62] bit → carry[63] ---")
    best_g62_bit = -1
    best_g62_corr = 0
    for b in range(32):
        g62_bit = np.array([(int(v) >> b) & 1 for v in g62_vals], dtype=float)
        if np.std(g62_bit) > 0:
            c = np.corrcoef(g62_bit, c63)[0, 1]
            if abs(c) > abs(best_g62_corr):
                best_g62_corr = c
                best_g62_bit = b
            if abs(c) > 0.02:
                print(f"    g62[b{b:2d}]: corr = {c:+.5f}")

    print(f"  Best: g62[b{best_g62_bit}] corr = {best_g62_corr:+.5f}")

    # SS63 distribution
    print(f"\n--- SS63 distribution ---")
    print(f"  E[SS63] = {np.mean(ss63):.0f} ({np.mean(ss63)/(1<<32):.4f} × 2^32)")
    print(f"  std[SS63] = {np.std(ss63):.0f}")
    print(f"  Threshold for carry=0: SS63 < {(1<<32) - K[63]} = {(1<<32) - K[63]}")
    print(f"  P(SS63 < threshold) = {np.mean(ss63 < ((1<<32) - K[63])):.5f}")

    # How many state[62] bits are determined by H[7,6,5]?
    print(f"\n--- Bits of state[62] accessible from H ---")
    print(f"  g[62] = H[7] - IV[7]: 32 bits (full)")
    print(f"  f[62] = H[6] - IV[6]: 32 bits (full)")
    print(f"  e[62] = H[5] - IV[5]: 32 bits (full)")
    print(f"  d[62] = e[63] - T1[63] = (H[4]-IV[4]) - T1[63]: needs T1[63]")
    print(f"  Total directly accessible: 96 bits / 256")
    print(f"  T1[63] = h[62] + Sig1(e[62]) + Ch(e62,f62,g62) + K[63] + W[63]")
    print(f"  All components known from H except h[62] and W[63]")

    return g62_vals, f62_vals, e62_vals


# ============================================================
# MAIN
# ============================================================

if __name__ == '__main__':
    print("K-Chain — Stage 1.4: Three Paths Around Chaos")
    print(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print()

    N1, N2, N3 = 2000, 1500, 5000
    if len(sys.argv) > 1 and sys.argv[1] == '--fast':
        N1, N2, N3 = 1000, 800, 3000

    t_start = time.time()
    k_info, det_rounds = path1_k_chain(N=N1)
    n_c0 = path2_multiblock(N=N2)
    g62, f62, e62 = path3_backward(N=N3)
    total = time.time() - t_start

    print(f"\n{'='*70}")
    print(f"TOTAL TIME: {total:.1f}s")
    print(f"{'='*70}")

    print(f"""
{'='*70}
SYNTHESIS: Three Paths Around the Chaotic Zone
{'='*70}

PATH 1 (K-CHAIN):
  {len(det_rounds)} rounds have K[r] large enough to force SS bits to 0
  Chain from r=63 backward: carry[63]=0 → constraints on state[62]
  T_CH_INVARIANT (methodology): Ch[62][b31,b30]=0 is one such constraint

PATH 2 (MULTI-BLOCK):
  carry[63]=0 in block 1 → IV2 biased → block 2 structure?
  Result depends on N (need many carry=0 events)

PATH 3 (BACKWARD):
  96 bits of state[62] directly computable from H[5,6,7]
  SS63 = h[62] + Sig1(e[62]) + Ch(e62,f62,g62)
  All components of SS63 accessible from H except h[62]

KEY INSIGHT: The chaotic zone blocks FORWARD information flow.
But BACKWARD (from output) gives direct access to state[62].
And K-chain gives DETERMINISTIC constraints at state[62].
The meeting point is state[62] — not through chaos, but around it.
""")
