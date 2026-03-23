#!/usr/bin/env python3
"""
CRAZY 42-45: Four structural analysis experiments on SHA-256 Wang chain.

CRAZY-42: Fault sensitivity (simulated) — flip bits in state at round 14
CRAZY-43: Higher-order differential — algebraic degree of W[14] -> De17[0]
CRAZY-44: State-split MITM — which state halves affect T1 at round 14
CRAZY-45: Correlation immunity — single/pair/triple bit correlations

All self-contained. Uses standard SHA-256 internals + Wang chain.
"""

import random
import math
import time
from collections import defaultdict
from itertools import combinations

MASK = 0xFFFFFFFF

K = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
    0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
]

H0 = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]

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
    return ((e & f) ^ (~e & g)) & MASK

def Maj(a, b, c):
    return ((a & b) ^ (a & c) ^ (b & c)) & MASK

def add32(*args):
    s = 0
    for x in args:
        s = (s + x) & MASK
    return s

def sha_round(st, w, k):
    a, b, c, d, e, f, g, h = st
    T1 = add32(h, Sig1(e), Ch(e, f, g), k, w)
    T2 = add32(Sig0(a), Maj(a, b, c))
    return [add32(T1, T2), a, b, c, add32(d, T1), e, f, g]

def expand_W(W16, num):
    W = list(W16)
    for i in range(16, num):
        W.append(add32(sig1(W[i-2]), W[i-7], sig0(W[i-15]), W[i-16]))
    return W

def rand32():
    return random.getrandbits(32)

def rand_msg16():
    return [rand32() for _ in range(16)]

def popcount(x):
    return bin(x).count('1')

# ── Wang chain: DW[0]=0x80000000, force De=0 for rounds 1-16 ──

def wang_chain_states(msg, max_round=18):
    """Run Wang chain and return (states_orig, states_prime, W_orig, W_prime)
    where states[t] is state AFTER round t (0-indexed), states[0] = IV."""
    iv = list(H0)
    W = list(msg)
    Wp = list(W)
    Wp[0] ^= 0x80000000

    s = list(iv)
    sp = list(iv)
    s = sha_round(s, W[0], K[0])
    sp = sha_round(sp, Wp[0], K[0])

    for t in range(1, 16):
        a, b, c, d, e, f, g, h = s
        a2, b2, c2, d2, e2, f2, g2, h2 = sp
        tp = add32(h, Sig1(e), Ch(e, f, g), K[t])
        tp2 = add32(h2, Sig1(e2), Ch(e2, f2, g2), K[t])
        target = add32(d, tp, W[t])
        Wp[t] = (target - d2 - tp2) & MASK
        s = sha_round(s, W[t], K[t])
        sp = sha_round(sp, Wp[t], K[t])

    We = expand_W(W, max_round)
    Wpe = expand_W(Wp, max_round)

    # Full run collecting all states
    s1 = list(iv); s2 = list(iv)
    states1 = [list(iv)]
    states2 = [list(iv)]
    for t in range(max_round):
        s1 = sha_round(s1, We[t], K[t])
        s2 = sha_round(s2, Wpe[t], K[t])
        states1.append(list(s1))
        states2.append(list(s2))
    return states1, states2, We, Wpe

def get_De17(msg):
    """Return De17 = e_orig[17] ^ e_prime[17]."""
    st1, st2, _, _ = wang_chain_states(msg, 18)
    return st1[17][4] ^ st2[17][4]

def get_De17_from_W14(w14_val):
    """Compute De17 with a random base message where W[14] is set to w14_val."""
    msg = rand_msg16()
    msg[14] = w14_val
    return get_De17(msg)

# ============================================================
# CRAZY-42: Fault Sensitivity (Simulated)
# ============================================================
def crazy42():
    print("=" * 72)
    print("CRAZY-42: Fault Sensitivity Analysis (Simulated)")
    print("=" * 72)
    print()
    print("Flip bit 0 of each register at round 14, measure De17 change.")
    print("1K messages, 8 registers x 1 bit flip each.")
    print()

    N = 1000
    reg_names = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']

    # fault_response[reg] = list of HW(De17_faulted XOR De17_original)
    fault_hw = {r: [] for r in range(8)}
    # sensitivity matrix: 8 regs x 32 De17 bits — count how often flipping reg changes De17 bit
    sens_count = [[0]*32 for _ in range(8)]

    for trial in range(N):
        msg = rand_msg16()
        st1, st2, We, Wpe = wang_chain_states(msg, 18)

        # Original De17
        De17_orig = st1[17][4] ^ st2[17][4]

        # State at round 14 (after 14 rounds)
        state14_orig = list(st1[14])
        state14_prime = list(st2[14])

        for reg in range(8):
            # Flip bit 0 of register 'reg' in the ORIGINAL state at round 14
            faulted_s1 = list(state14_orig)
            faulted_s1[reg] ^= 1  # flip bit 0

            # Re-run rounds 14-17 from faulted state
            s_f = list(faulted_s1)
            s_p = list(state14_prime)
            for t in range(14, 17):
                s_f = sha_round(s_f, We[t], K[t])
                s_p = sha_round(s_p, Wpe[t], K[t])

            De17_faulted = s_f[4] ^ s_p[4]
            diff = De17_faulted ^ De17_orig
            hw = popcount(diff)
            fault_hw[reg].append(hw)

            for b in range(32):
                if (diff >> b) & 1:
                    sens_count[reg][b] += 1

    print("  Register   Mean HW(De17 change)   Std Dev")
    print("  " + "-" * 50)
    means = []
    for reg in range(8):
        vals = fault_hw[reg]
        mean = sum(vals) / len(vals)
        var = sum((v - mean)**2 for v in vals) / len(vals)
        std = var ** 0.5
        means.append(mean)
        print(f"  {reg_names[reg]:>8s}   {mean:18.4f}   {std:8.4f}")

    most = reg_names[means.index(max(means))]
    least = reg_names[means.index(min(means))]
    print()
    print(f"  MOST disruptive flip:  register '{most}' (mean HW = {max(means):.4f})")
    print(f"  LEAST disruptive flip: register '{least}' (mean HW = {min(means):.4f})")

    # Sensitivity matrix rank (binary: threshold at N/2 = significant)
    threshold = N * 0.1  # 10% threshold for "affects"
    binary_matrix = []
    for reg in range(8):
        row = [1 if sens_count[reg][b] > threshold else 0 for b in range(32)]
        binary_matrix.append(row)

    # Compute rank via Gaussian elimination over GF(2)
    mat = [list(row) for row in binary_matrix]
    rows = len(mat)
    cols = len(mat[0])
    rank = 0
    for col in range(cols):
        pivot = None
        for r in range(rank, rows):
            if mat[r][col]:
                pivot = r
                break
        if pivot is None:
            continue
        mat[rank], mat[pivot] = mat[pivot], mat[rank]
        for r in range(rows):
            if r != rank and mat[r][col]:
                mat[r] = [(mat[r][c] ^ mat[rank][c]) for c in range(cols)]
        rank += 1

    print(f"\n  Binary sensitivity matrix (8x32) rank over GF(2): {rank}")
    print(f"  (max possible = 8; indicates linear independence of fault effects)")

    # Show the binary matrix compactly
    print("\n  Sensitivity heatmap (reg vs De17 bits 0-31, '1' = affected):")
    print(f"  {'Reg':>4s}  {'Bits 0-31':32s}")
    for reg in range(8):
        row_str = ''.join(str(binary_matrix[reg][b]) for b in range(32))
        print(f"  {reg_names[reg]:>4s}  {row_str}")

    print()
    print("  INTERPRETATION:")
    print(f"  - Registers e,f,g,h feed directly into T1 -> expect high sensitivity")
    print(f"  - Registers a,b,c,d affect next round only via shift -> lower sensitivity")
    print(f"  - Rank {rank}/8: {'full rank — all faults linearly independent' if rank == 8 else 'deficient — some fault effects are linearly dependent'}")
    print()


# ============================================================
# CRAZY-43: Higher-Order Differential
# ============================================================
def crazy43():
    print("=" * 72)
    print("CRAZY-43: Higher-Order Differential Analysis")
    print("=" * 72)
    print()
    print("f: W[14] -> De17[bit 0]. Test algebraic degree via D^k f.")
    print("For k=1..10, sample 20 random k-dim affine subspaces.")
    print()

    # Use a fixed base message, vary only W[14]
    base_msg = rand_msg16()

    def f(w14):
        """Evaluate f: W[14] -> De17[bit 0]."""
        msg = list(base_msg)
        msg[14] = w14 & MASK
        de17 = get_De17(msg)
        return de17 & 1

    max_k = 10
    n_subspaces = 20

    print(f"  {'k':>3s}  {'#zero':>6s}  {'#total':>7s}  {'frac_zero':>10s}  {'degree < k?':>12s}")
    print("  " + "-" * 50)

    algebraic_degree = None
    for k in range(1, max_k + 1):
        zero_count = 0
        for _ in range(n_subspaces):
            # Random affine subspace: base point + k linearly independent directions
            base = rand32()
            directions = []
            for _ in range(k):
                d = rand32()
                while d == 0:
                    d = rand32()
                directions.append(d)

            # Evaluate f over all 2^k points in the subspace
            xor_sum = 0
            for mask in range(1 << k):
                point = base
                for bit_idx in range(k):
                    if (mask >> bit_idx) & 1:
                        point ^= directions[bit_idx]
                point &= MASK
                xor_sum ^= f(point)

            if xor_sum == 0:
                zero_count += 1

        frac = zero_count / n_subspaces
        is_zero = "YES" if zero_count == n_subspaces else "no"
        print(f"  {k:3d}  {zero_count:6d}  {n_subspaces:7d}  {frac:10.4f}  {is_zero:>12s}")

        if zero_count == n_subspaces and algebraic_degree is None:
            algebraic_degree = k - 1

    print()
    if algebraic_degree is not None:
        print(f"  Algebraic degree upper bound: {algebraic_degree}")
        print(f"  (D^{algebraic_degree+1} f = 0 for all tested subspaces)")
    else:
        print(f"  Algebraic degree >= {max_k} (no k found where D^k f is always 0)")
        print(f"  Expected: SHA-256 is effectively full-degree (32) in W[14]")

    print()
    print("  INTERPRETATION:")
    print("  - Random Boolean function has degree n-1 in n variables")
    print("  - SHA-256 round function involves additions and nonlinear ops")
    print("  - If degree < 32, there may be exploitable algebraic structure")
    print("  - Fraction of zeros near 0.5 for all k -> behaves like random function")
    print()


# ============================================================
# CRAZY-44: State-Split MITM
# ============================================================
def crazy44():
    print("=" * 72)
    print("CRAZY-44: State-Split Meet-in-the-Middle Analysis")
    print("=" * 72)
    print()
    print("At round 14: T1 = h + Sig1(e) + Ch(e,f,g) + K[14] + W[14]")
    print("e_new = d + T1. Which half {a,b,c,d} vs {e,f,g,h} affects T1?")
    print()

    N = 1000
    reg_names_left = ['a', 'b', 'c', 'd']
    reg_names_right = ['e', 'f', 'g', 'h']

    # For each register bit (bit 0), measure how much it affects T1
    # We do this by flipping bit 0 in the register and measuring T1 change

    print("  PART 1: Direct T1 sensitivity (flip bit 0 of each register)")
    print()

    t1_sensitivity = {}

    for reg_idx in range(8):
        reg_name = (['a','b','c','d','e','f','g','h'])[reg_idx]
        diff_count = 0
        for _ in range(N):
            state = [rand32() for _ in range(8)]
            W14 = rand32()
            a, b, c, d, e, f, g, h = state

            T1_orig = add32(h, Sig1(e), Ch(e, f, g), K[14], W14)

            faulted = list(state)
            faulted[reg_idx] ^= 1
            a2, b2, c2, d2, e2, f2, g2, h2 = faulted
            T1_faulted = add32(h2, Sig1(e2), Ch(e2, f2, g2), K[14], W14)

            if T1_orig != T1_faulted:
                diff_count += 1

        frac = diff_count / N
        t1_sensitivity[reg_name] = frac
        print(f"    {reg_name}: T1 changes in {frac*100:.1f}% of cases")

    left_bits = sum(1 for r in reg_names_left if t1_sensitivity[r] > 0.01)
    right_bits = sum(1 for r in reg_names_right if t1_sensitivity[r] > 0.01)

    print()
    print(f"  Left half {{a,b,c,d}} bits affecting T1:  {left_bits}")
    print(f"  Right half {{e,f,g,h}} bits affecting T1: {right_bits}")
    print()

    # Part 2: e_new = d + T1 — d is from left half
    print("  PART 2: e_new = d + T1 sensitivity")
    print()

    enew_sensitivity = {}
    for reg_idx in range(8):
        reg_name = (['a','b','c','d','e','f','g','h'])[reg_idx]
        diff_count = 0
        for _ in range(N):
            state = [rand32() for _ in range(8)]
            W14 = rand32()
            a, b, c, d, e, f, g, h = state

            T1 = add32(h, Sig1(e), Ch(e, f, g), K[14], W14)
            e_new_orig = add32(d, T1)

            faulted = list(state)
            faulted[reg_idx] ^= 1
            a2, b2, c2, d2, e2, f2, g2, h2 = faulted
            T1_f = add32(h2, Sig1(e2), Ch(e2, f2, g2), K[14], W14)
            e_new_faulted = add32(d2, T1_f)

            if e_new_orig != e_new_faulted:
                diff_count += 1

        frac = diff_count / N
        enew_sensitivity[reg_name] = frac
        print(f"    {reg_name}: e_new changes in {frac*100:.1f}% of cases")

    left_enew = sum(1 for r in reg_names_left if enew_sensitivity[r] > 0.01)
    right_enew = sum(1 for r in reg_names_right if enew_sensitivity[r] > 0.01)

    print()
    print(f"  Left {{a,b,c,d}} affecting e_new:  {left_enew} registers")
    print(f"  Right {{e,f,g,h}} affecting e_new: {right_enew} registers")

    mitm_quality = min(left_enew, right_enew)
    print(f"\n  MITM split quality = min(left, right) = {mitm_quality}")

    # Part 3: Multi-round propagation (partial correlation through rounds 14-17)
    print("\n  PART 3: Multi-round MITM — correlation of left/right halves to De17")
    print()

    left_corr_hw = []
    right_corr_hw = []

    for _ in range(N):
        msg = rand_msg16()
        st1, st2, We, Wpe = wang_chain_states(msg, 18)
        De17_orig = st1[17][4] ^ st2[17][4]

        state14 = list(st1[14])
        state14p = list(st2[14])

        # Flip all left half bits (register d, bit 0 — most direct)
        faulted_left = list(state14)
        faulted_left[3] ^= 1  # d
        s_f = list(faulted_left)
        s_p = list(state14p)
        for t in range(14, 17):
            s_f = sha_round(s_f, We[t], K[t])
            s_p = sha_round(s_p, Wpe[t], K[t])
        De17_left = s_f[4] ^ s_p[4]
        left_corr_hw.append(popcount(De17_left ^ De17_orig))

        # Flip right half (register e, bit 0)
        faulted_right = list(state14)
        faulted_right[4] ^= 1  # e
        s_f = list(faulted_right)
        s_p = list(state14p)
        for t in range(14, 17):
            s_f = sha_round(s_f, We[t], K[t])
            s_p = sha_round(s_p, Wpe[t], K[t])
        De17_right = s_f[4] ^ s_p[4]
        right_corr_hw.append(popcount(De17_right ^ De17_orig))

    mean_left = sum(left_corr_hw) / len(left_corr_hw)
    mean_right = sum(right_corr_hw) / len(right_corr_hw)

    print(f"    Flip d[0] -> mean HW(De17 change): {mean_left:.4f}")
    print(f"    Flip e[0] -> mean HW(De17 change): {mean_right:.4f}")

    print()
    print("  INTERPRETATION:")
    print("  - T1 depends DIRECTLY on {e,f,g,h} (right half) and W[14]")
    print("  - T1 does NOT depend directly on {a,b,c} (only d enters via e_new = d + T1)")
    print(f"  - MITM quality = {mitm_quality}: {'good split' if mitm_quality >= 2 else 'poor split — one side dominates'}")
    print("  - For attack: right half controls T1, left half only shifts via d")
    print()


# ============================================================
# CRAZY-45: Correlation Immunity
# ============================================================
def crazy45():
    print("=" * 72)
    print("CRAZY-45: Correlation Immunity of W[14] -> De17[bit 0]")
    print("=" * 72)
    print()
    print("Measure correlation of input bit subsets with output bit.")
    print()

    base_msg = rand_msg16()
    N_single = 50000
    N_pair = 50000
    N_triple = 50000

    # Pre-compute: generate N samples of (W14, De17_bit0)
    print(f"  Generating {N_single} samples...")
    t0 = time.time()
    w14_vals = [rand32() for _ in range(N_single)]
    de17_bit0 = []
    for w14 in w14_vals:
        msg = list(base_msg)
        msg[14] = w14
        de17 = get_De17(msg)
        de17_bit0.append(de17 & 1)

    elapsed = time.time() - t0
    print(f"  Generated in {elapsed:.1f}s")

    # Output bias
    ones = sum(de17_bit0)
    bias = ones / N_single - 0.5
    print(f"\n  Output bias: P(De17[0]=1) = {ones/N_single:.6f}, bias = {bias:+.6f}")

    # ── Single-bit correlations ──
    print(f"\n  SINGLE-BIT CORRELATIONS (32 bits of W[14]):")
    single_corr = []
    for i in range(32):
        # corr = E[(-1)^(x_i XOR f(x))] = (agree - disagree) / N
        agree = 0
        for idx in range(N_single):
            xi = (w14_vals[idx] >> i) & 1
            agree += 1 if xi == de17_bit0[idx] else -1
        corr = agree / N_single
        single_corr.append(corr)

    # Show top 5 and bottom 5
    indexed = sorted(enumerate(single_corr), key=lambda x: abs(x[1]), reverse=True)
    print(f"    Top 5 by |correlation|:")
    for bit, c in indexed[:5]:
        print(f"      bit {bit:2d}: corr = {c:+.6f}")
    print(f"    Bottom 5 by |correlation|:")
    for bit, c in indexed[-5:]:
        print(f"      bit {bit:2d}: corr = {c:+.6f}")

    max_single = max(abs(c) for c in single_corr)
    print(f"    Max |single-bit corr| = {max_single:.6f}")

    ci_order = 0
    CI_THRESHOLD = 0.02  # threshold for "zero" correlation

    all_single_zero = all(abs(c) < CI_THRESHOLD for c in single_corr)
    if all_single_zero:
        ci_order = 1
        print(f"    All single-bit correlations < {CI_THRESHOLD} -> CI order >= 1")
    else:
        print(f"    Some single-bit correlations >= {CI_THRESHOLD} -> CI order = 0")

    # ── Pair correlations ──
    print(f"\n  PAIR CORRELATIONS (all C(32,2) = 496 pairs):")
    pair_corrs = []
    for i, j in combinations(range(32), 2):
        agree = 0
        for idx in range(N_pair):
            xij = ((w14_vals[idx] >> i) & 1) ^ ((w14_vals[idx] >> j) & 1)
            agree += 1 if xij == de17_bit0[idx] else -1
        corr = agree / N_pair
        pair_corrs.append((i, j, corr))

    pair_corrs_sorted = sorted(pair_corrs, key=lambda x: abs(x[2]), reverse=True)
    print(f"    Top 5 pairs by |correlation|:")
    for i, j, c in pair_corrs_sorted[:5]:
        print(f"      bits ({i:2d},{j:2d}): corr = {c:+.6f}")

    max_pair = max(abs(c) for _, _, c in pair_corrs)
    print(f"    Max |pair corr| = {max_pair:.6f}")

    if ci_order >= 1:
        all_pair_zero = all(abs(c) < CI_THRESHOLD for _, _, c in pair_corrs)
        if all_pair_zero:
            ci_order = 2
            print(f"    All pair correlations < {CI_THRESHOLD} -> CI order >= 2")
        else:
            print(f"    Some pair correlations >= {CI_THRESHOLD} -> CI order = 1")

    # ── Triple correlations (sample 1000 random triples) ──
    print(f"\n  TRIPLE CORRELATIONS (1000 random triples):")
    triple_corrs = []
    for _ in range(1000):
        i, j, k = random.sample(range(32), 3)
        agree = 0
        for idx in range(N_triple):
            xijk = ((w14_vals[idx] >> i) & 1) ^ ((w14_vals[idx] >> j) & 1) ^ ((w14_vals[idx] >> k) & 1)
            agree += 1 if xijk == de17_bit0[idx] else -1
        corr = agree / N_triple
        triple_corrs.append((i, j, k, corr))

    triple_sorted = sorted(triple_corrs, key=lambda x: abs(x[3]), reverse=True)
    print(f"    Top 5 triples by |correlation|:")
    for i, j, k, c in triple_sorted[:5]:
        print(f"      bits ({i:2d},{j:2d},{k:2d}): corr = {c:+.6f}")

    max_triple = max(abs(c) for _, _, _, c in triple_corrs)
    print(f"    Max |triple corr| = {max_triple:.6f}")

    if ci_order >= 2:
        all_triple_zero = all(abs(c) < CI_THRESHOLD for _, _, _, c in triple_corrs)
        if all_triple_zero:
            ci_order = 3
            print(f"    All triple correlations < {CI_THRESHOLD} -> CI order >= 3")
        else:
            print(f"    Some triple correlations >= {CI_THRESHOLD} -> CI order = 2")

    print()
    print(f"  CORRELATION IMMUNITY ORDER: {ci_order}")
    print()
    print("  INTERPRETATION:")
    print("  - CI order 0: some single input bits are correlated with output")
    print("  - CI order >= 1: no single bit predicts output (ideal for crypto)")
    print("  - Higher CI order -> harder to mount correlation attacks")
    print(f"  - For SHA-256, expect CI order 0 (additions create bit correlations)")
    print(f"  - Max correlations: single={max_single:.4f}, pair={max_pair:.4f}, triple={max_triple:.4f}")
    print()


# ============================================================
# Main
# ============================================================
if __name__ == "__main__":
    random.seed(42)
    t_start = time.time()

    crazy42()
    t1 = time.time()
    print(f"[CRAZY-42 completed in {t1 - t_start:.1f}s]\n")

    crazy43()
    t2 = time.time()
    print(f"[CRAZY-43 completed in {t2 - t1:.1f}s]\n")

    crazy44()
    t3 = time.time()
    print(f"[CRAZY-44 completed in {t3 - t2:.1f}s]\n")

    crazy45()
    t4 = time.time()
    print(f"[CRAZY-45 completed in {t4 - t3:.1f}s]\n")

    print("=" * 72)
    print(f"ALL CRAZY 42-45 COMPLETE — total time: {t4 - t_start:.1f}s")
    print("=" * 72)
