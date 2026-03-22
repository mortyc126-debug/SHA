"""
crazy25_stochastic_linearization.py
Stochastic linearization of SHA-256 barrier De17(W[14])

Idea: Instead of linearizing at ONE point (residual ~15.7 bits),
use 2^K base points. For each test input, find nearest base point
(Hamming distance) and use its precomputed Jacobian.

If the barrier function is "locally linear" in small Hamming
neighborhoods, the residual should drop as K increases.
"""

import random
import time
import struct

# ─── SHA-256 primitives ───
MASK = 0xFFFFFFFF
K_SHA = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
]
H0 = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]

def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK
def Sig0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sig1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def sig0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sig1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)
def Ch(e, f, g): return ((e & f) ^ (~e & g)) & MASK
def Maj(a, b, c): return ((a & b) ^ (a & c) ^ (b & c)) & MASK

def add32(*a):
    s = 0
    for x in a:
        s = (s + x) & MASK
    return s

def hw(x): return bin(x & MASK).count('1')

def sha_round(st, w, k):
    a, b, c, d, e, f, g, h = st
    T1 = add32(h, Sig1(e), Ch(e, f, g), k, w)
    T2 = add32(Sig0(a), Maj(a, b, c))
    return [add32(T1, T2), a, b, c, add32(d, T1), e, f, g]

def expand_W(W16):
    W = list(W16)
    for i in range(16, 20):
        W.append(add32(sig1(W[i-2]), W[i-7], sig0(W[i-15]), W[i-16]))
    return W

def get_de17(msg, iv=None):
    """Compute De17 = e17 XOR e17' where msg'[0] = msg[0] ^ 0x80000000
    and msg'[1..15] are adjusted to cancel differences through round 15."""
    if iv is None:
        iv = list(H0)
    W = list(msg)
    Wp = list(W)
    Wp[0] ^= 0x80000000
    s = list(iv)
    sp = list(iv)
    s = sha_round(s, W[0], K_SHA[0])
    sp = sha_round(sp, Wp[0], K_SHA[0])
    for t in range(1, 16):
        a, b, c, d, e, f, g, h = s
        a2, b2, c2, d2, e2, f2, g2, h2 = sp
        tp = add32(h, Sig1(e), Ch(e, f, g), K_SHA[t])
        tp2 = add32(h2, Sig1(e2), Ch(e2, f2, g2), K_SHA[t])
        target = add32(d, tp, W[t])
        Wp[t] = (target - d2 - tp2) & MASK
        s = sha_round(s, W[t], K_SHA[t])
        sp = sha_round(sp, Wp[t], K_SHA[t])
    We = expand_W(msg)
    Wpe = expand_W(Wp)
    s1 = list(iv)
    s2 = list(iv)
    for t in range(17):
        s1 = sha_round(s1, We[t], K_SHA[t])
        s2 = sha_round(s2, Wpe[t], K_SHA[t])
    return s1[4] ^ s2[4]


def make_msg_with_w14(base_msg, w14_val):
    """Return a copy of base_msg with W[14] replaced."""
    m = list(base_msg)
    m[14] = w14_val & MASK
    return m


def compute_jacobian(base_msg, w14_center):
    """Compute 32x32 GF(2) Jacobian of De17 w.r.t. W[14] at w14_center.
    Returns (de17_center, jacobian_columns) where jacobian_columns[j] is
    the XOR response of De17 to flipping bit j of W[14]."""
    msg_c = make_msg_with_w14(base_msg, w14_center)
    de17_c = get_de17(msg_c)
    jac = [0] * 32
    for j in range(32):
        w14_flip = w14_center ^ (1 << j)
        msg_f = make_msg_with_w14(base_msg, w14_flip)
        de17_f = get_de17(msg_f)
        jac[j] = de17_c ^ de17_f  # response column
    return de17_c, jac


def linear_predict(de17_center, jacobian, w14_center, w14_target):
    """Predict De17(w14_target) using linear approximation around w14_center.
    pred = de17_center XOR (sum of jacobian columns for differing bits)."""
    delta = w14_center ^ w14_target
    pred = de17_center
    while delta:
        j = (delta & -delta).bit_length() - 1
        pred ^= jacobian[j]
        delta &= delta - 1
    return pred


def hamming_dist(a, b):
    return hw(a ^ b)


def main():
    random.seed(0xC025)

    # Fixed base message (16 words)
    base_msg = [random.getrandbits(32) for _ in range(16)]
    print("Base message W[0..3]:", [hex(w) for w in base_msg[:4]])
    print()

    # ─── Single-point baseline ───
    print("=" * 70)
    print("SINGLE-POINT BASELINE (linearization at one random point)")
    print("=" * 70)
    w14_single = random.getrandbits(32)
    de17_single, jac_single = compute_jacobian(base_msg, w14_single)

    N_TEST = 10000
    test_w14s = [random.getrandbits(32) for _ in range(N_TEST)]

    # Compute true De17 for all test points
    print(f"Computing true De17 for {N_TEST} test points...")
    t0 = time.time()
    true_de17 = []
    for w in test_w14s:
        msg = make_msg_with_w14(base_msg, w)
        true_de17.append(get_de17(msg))
    t_true = time.time() - t0
    print(f"  Done in {t_true:.1f}s")

    # Single-point residuals
    residuals_single = []
    for i, w in enumerate(test_w14s):
        pred = linear_predict(de17_single, jac_single, w14_single, w)
        residuals_single.append(hw(pred ^ true_de17[i]))
    avg_single = sum(residuals_single) / len(residuals_single)
    min_single = min(residuals_single)
    print(f"  Single-point avg residual HW: {avg_single:.2f}")
    print(f"  Single-point min residual HW: {min_single}")
    print()

    # ─── Multi-point stochastic linearization ───
    print("=" * 70)
    print("STOCHASTIC LINEARIZATION (multiple base points)")
    print("=" * 70)

    results = []

    # Try K = 8, 10, 12, 14 (but time-gate K>12)
    for K in [8, 10, 12, 14]:
        num_bp = 1 << K
        print(f"\n--- K={K}, base points = 2^{K} = {num_bp} ---")

        # Time estimate from previous iteration
        if results and K >= 14:
            # Extrapolate: each K+2 is 4x slower
            est_time = results[-1][4] * 4
            if est_time > 600:
                print(f"  Estimated time: {est_time:.0f}s > 600s, SKIPPING (will extrapolate)")
                results.append((K, num_bp, None, None, None))
                continue

        t0 = time.time()

        # Generate base points (random W[14] values)
        base_points = [random.getrandbits(32) for _ in range(num_bp)]

        # Compute Jacobians at all base points
        print(f"  Computing {num_bp} Jacobians ({num_bp * 32} De17 calls)...")
        base_data = []  # (w14, de17, jacobian)
        for idx, bp in enumerate(base_points):
            de17_bp, jac_bp = compute_jacobian(base_msg, bp)
            base_data.append((bp, de17_bp, jac_bp))
            if (idx + 1) % max(1, num_bp // 4) == 0:
                elapsed = time.time() - t0
                print(f"    {idx+1}/{num_bp} done ({elapsed:.1f}s)")

        t_jac = time.time() - t0
        print(f"  Jacobians computed in {t_jac:.1f}s")

        # For each test point, find nearest base point and predict
        print(f"  Evaluating {N_TEST} test points...")
        residuals = []
        hdist_used = []
        for i, w in enumerate(test_w14s):
            # Find nearest base point
            best_dist = 33
            best_idx = 0
            for j, (bp, _, _) in enumerate(base_data):
                d = hamming_dist(w, bp)
                if d < best_dist:
                    best_dist = d
                    best_idx = j
                    if d == 0:
                        break
            bp_w14, bp_de17, bp_jac = base_data[best_idx]
            pred = linear_predict(bp_de17, bp_jac, bp_w14, w)
            res = hw(pred ^ true_de17[i])
            residuals.append(res)
            hdist_used.append(best_dist)

        t_total = time.time() - t0
        avg_res = sum(residuals) / len(residuals)
        min_res = min(residuals)
        avg_hdist = sum(hdist_used) / len(hdist_used)
        min_hdist = min(hdist_used)

        print(f"  Avg nearest Hamming dist: {avg_hdist:.2f}, min: {min_hdist}")
        print(f"  Avg residual HW: {avg_res:.2f}")
        print(f"  Min residual HW: {min_res}")
        print(f"  Total time: {t_total:.1f}s")

        results.append((K, num_bp, avg_res, min_res, t_total))

    # ─── Summary table ───
    print("\n" + "=" * 70)
    print("SUMMARY TABLE")
    print("=" * 70)
    print(f"{'K':>4} {'#bases':>8} {'avg_resid':>10} {'min_resid':>10} {'time(s)':>10} {'eff_barrier':>14}")
    print("-" * 60)
    print(f"{'0':>4} {'1':>8} {avg_single:>10.2f} {min_single:>10} {'--':>10} {'2^' + f'{avg_single:.1f}':>14}")

    for K, num_bp, avg_res, min_res, t_total in results:
        if avg_res is not None:
            eff = f"2^{avg_res:.1f}"
            print(f"{K:>4} {num_bp:>8} {avg_res:>10.2f} {min_res:>10} {t_total:>10.1f} {eff:>14}")
        else:
            # Extrapolate from trend
            # Use last two valid points for log-linear extrapolation
            valid = [(r[0], r[2]) for r in results if r[2] is not None]
            if len(valid) >= 2:
                k1, r1 = valid[-2]
                k2, r2 = valid[-1]
                slope = (r2 - r1) / (k2 - k1)
                extrap = r2 + slope * (K - k2)
                eff = f"~2^{extrap:.1f}"
                print(f"{K:>4} {num_bp:>8} {'~'+f'{extrap:.2f}':>10} {'(extrap)':>10} {'SKIP':>10} {eff:>14}")
            else:
                print(f"{K:>4} {num_bp:>8} {'SKIPPED':>10} {'':>10} {'':>10}")

    # ─── Analysis ───
    print("\n" + "=" * 70)
    print("ANALYSIS")
    print("=" * 70)

    valid = [(r[0], r[2]) for r in results if r[2] is not None]
    if len(valid) >= 2:
        k_first, r_first = valid[0]
        k_last, r_last = valid[-1]
        improvement = avg_single - r_last
        improvement_per_K = improvement / (k_last - k_first) if k_last != k_first else 0

        print(f"Single-point residual:       {avg_single:.2f} bits")
        print(f"Best multi-point residual:   {r_last:.2f} bits (K={k_last}, {1<<k_last} bases)")
        print(f"Total improvement:           {improvement:.2f} bits")
        print(f"Improvement per doubling K:  {improvement_per_K * 2:.2f} bits")
        print()

        # How many base points to reach residual = 0?
        if improvement > 0 and improvement_per_K > 0:
            k_needed = k_first + avg_single / improvement_per_K
            print(f"Extrapolated K for residual=0: K~{k_needed:.1f} => 2^{k_needed:.1f} base points")
            print(f"That means precomputing 2^{k_needed:.1f} Jacobians = 2^{k_needed+5:.1f} De17 calls")
        else:
            print("No meaningful improvement observed.")

        print()
        print("CONCLUSION:")
        if improvement < 1.0:
            print("  Stochastic linearization provides NEGLIGIBLE improvement (<1 bit).")
            print("  The De17 barrier is NOT locally linear even in small Hamming neighborhoods.")
            print("  The nonlinearity is INTRINSIC, not an artifact of linearization point choice.")
            print("  Effective barrier remains ~2^16 per word — confirming SHA-256's resistance.")
        elif improvement < 3.0:
            print("  Stochastic linearization provides MODEST improvement (1-3 bits).")
            print("  Some local linearity exists, but the barrier reduction is insufficient")
            print("  to meaningfully accelerate collision search.")
        else:
            print("  Stochastic linearization provides SIGNIFICANT improvement (>3 bits).")
            print("  This suggests exploitable local linearity in the De17 barrier.")

    print("\nDone.")


if __name__ == "__main__":
    main()
