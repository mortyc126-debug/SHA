#!/usr/bin/env python3
"""
CRAZY-8: Barrier Function Correlation Analysis
===============================================
Critical question: Are De[17], De[18], ..., De[24] INDEPENDENT
as functions of neutral bits?

If they're correlated, the EFFECTIVE dimension of the barrier space
is less than 8×32=256 bits. Birthday on a k-dimensional space costs
2^(k/2), so correlation directly reduces attack cost.

We measure:
1. Pairwise correlation between De[t] and De[t'] over neutral bit choices
2. Joint entropy of (De[17], De[18]) — is it 64 bits or less?
3. Conditional entropy: H(De[18] | De[17]=small) — does low De[17]
   predict anything about De[18]?
4. The EFFECTIVE dimension of the barrier manifold

This determines whether birthday on barriers could be cheaper than 2^128.
"""

import random
import time
import math
from collections import Counter

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

def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK
def Sig0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sig1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def sig0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sig1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)
def Ch(e, f, g): return ((e & f) ^ (~e & g)) & MASK
def Maj(a, b, c): return ((a & b) ^ (a & c) ^ (b & c)) & MASK
def add32(*args):
    s = 0
    for a in args: s = (s + a) & MASK
    return s
def hw(x): return bin(x & MASK).count('1')

def sha256_round(state, W_t, K_t):
    a, b, c, d, e, f, g, h = state
    T1 = add32(h, Sig1(e), Ch(e, f, g), K_t, W_t)
    T2 = add32(Sig0(a), Maj(a, b, c))
    return [add32(T1, T2), a, b, c, add32(d, T1), e, f, g]

def expand_W(W16):
    W = list(W16)
    for i in range(16, 64):
        W.append(add32(sig1(W[i-2]), W[i-7], sig0(W[i-15]), W[i-16]))
    return W

def run_rounds(init_state, W, num_rounds):
    states = [list(init_state)]
    for t in range(num_rounds):
        s = sha256_round(states[-1], W[t], K[t])
        states.append(s)
    return states

def compute_wang_chain(base_msg, init_state, dw0=0x80000000):
    W = list(base_msg)
    W_prime = list(W)
    W_prime[0] = W[0] ^ dw0

    state = list(init_state)
    state_prime = list(init_state)

    state = sha256_round(state, W[0], K[0])
    state_prime = sha256_round(state_prime, W_prime[0], K[0])

    for t in range(1, 16):
        a, b, c, d, e, f, g, h = state
        a2, b2, c2, d2, e2, f2, g2, h2 = state_prime
        T1_partial = add32(h, Sig1(e), Ch(e, f, g), K[t])
        T1_partial_prime = add32(h2, Sig1(e2), Ch(e2, f2, g2), K[t])
        target_e = add32(d, T1_partial, W[t])
        W_prime[t] = (target_e - d2 - T1_partial_prime) & MASK

        state = sha256_round(state, W[t], K[t])
        state_prime = sha256_round(state_prime, W_prime[t], K[t])

    return W_prime

def find_neutral_bits(base_msg, init_state):
    W_prime = compute_wang_chain(base_msg, init_state)
    neutrals = []
    for k in range(1, 16):
        for b in range(32):
            flip = 1 << b
            W_f = list(base_msg)
            W_f[k] ^= flip
            Wp_f = list(W_prime)
            Wp_f[k] ^= flip

            W_exp = expand_W(W_f)
            Wp_exp = expand_W(Wp_f)
            s1 = run_rounds(init_state, W_exp, 16)
            s2 = run_rounds(init_state, Wp_exp, 16)

            ok = True
            for t in range(2, 17):
                if (s1[t][4] - s2[t][4]) & MASK != 0:
                    ok = False
                    break
            if ok:
                neutrals.append((k, b))
    return neutrals


def sample_barriers(base_msg, init_state, neutrals, n_samples, max_round=21):
    """Sample De[17..max_round] for random neutral bit subsets."""
    results = []
    for _ in range(n_samples):
        W = list(base_msg)
        W_prime = compute_wang_chain(base_msg, init_state)
        # Flip random subset of neutral bits
        for k, b in neutrals:
            if random.random() < 0.5:
                flip = 1 << b
                W[k] ^= flip
                W_prime[k] ^= flip

        # Recompute Wang chain
        W_prime = compute_wang_chain(W, init_state)
        W_exp = expand_W(W)
        Wp_exp = expand_W(W_prime)

        s1 = run_rounds(init_state, W_exp, max_round)
        s2 = run_rounds(init_state, Wp_exp, max_round)

        de_vals = {}
        for t in range(17, min(max_round + 1, len(s1))):
            de_vals[t] = (s1[t][4] ^ s2[t][4]) & MASK
        results.append(de_vals)

    return results


def entropy_from_samples(values):
    """Estimate entropy from a list of values."""
    n = len(values)
    if n == 0:
        return 0
    counts = Counter(values)
    ent = 0
    for c in counts.values():
        p = c / n
        if p > 0:
            ent -= p * math.log2(p)
    return ent


def main():
    print("=" * 72)
    print("CRAZY-8: Barrier Function Correlation Analysis")
    print("=" * 72)
    print()

    random.seed(0xC0DE8)
    t_start = time.time()
    init_state = list(H0)

    NUM_MSGS = 5
    N_SAMPLES = 2000  # per message

    all_correlations = []
    all_joint_entropies = []
    all_conditional = []
    all_de17_hw = []
    all_de18_hw = []
    all_bit_correlations = {(t1, t2): [] for t1 in range(17, 21) for t2 in range(t1+1, 21)}

    for msg_idx in range(NUM_MSGS):
        elapsed = time.time() - t_start
        if elapsed > 400:
            print(f"[Time limit at message {msg_idx}]")
            break

        print(f"Message {msg_idx+1}/{NUM_MSGS} (elapsed: {elapsed:.1f}s)")
        base_msg = [random.getrandbits(32) for _ in range(16)]

        neutrals = find_neutral_bits(base_msg, init_state)
        print(f"  Found {len(neutrals)} neutral bits")

        if len(neutrals) == 0:
            continue

        # Sample barrier values
        samples = sample_barriers(base_msg, init_state, neutrals, N_SAMPLES, max_round=21)

        # ── Analysis 1: Per-round entropy ──────────────────────────
        for t in range(17, 21):
            vals = [s[t] for s in samples if t in s]
            ent = entropy_from_samples(vals)
            hw_vals = [hw(v) for v in vals]
            avg_hw = sum(hw_vals) / len(hw_vals) if hw_vals else 0
            if t == 17:
                all_de17_hw.extend(hw_vals)
            elif t == 18:
                all_de18_hw.extend(hw_vals)
            print(f"  De[{t}]: entropy={ent:.1f} bits (of {math.log2(len(vals)):.1f} possible), "
                  f"avg HW={avg_hw:.1f}, unique={len(set(vals))}/{len(vals)}")

        # ── Analysis 2: Pairwise bit-level correlation ─────────────
        for t1 in range(17, 21):
            for t2 in range(t1+1, 21):
                vals1 = [s.get(t1, 0) for s in samples]
                vals2 = [s.get(t2, 0) for s in samples]

                # Bit-level correlation
                corr_sum = 0
                for bit in range(32):
                    mask = 1 << bit
                    bits1 = [(v >> bit) & 1 for v in vals1]
                    bits2 = [(v >> bit) & 1 for v in vals2]
                    # Pearson correlation
                    n = len(bits1)
                    mean1 = sum(bits1) / n
                    mean2 = sum(bits2) / n
                    cov = sum((b1-mean1)*(b2-mean2) for b1, b2 in zip(bits1, bits2)) / n
                    std1 = (sum((b-mean1)**2 for b in bits1) / n) ** 0.5
                    std2 = (sum((b-mean2)**2 for b in bits2) / n) ** 0.5
                    if std1 > 0 and std2 > 0:
                        r = cov / (std1 * std2)
                        corr_sum += abs(r)
                    # else: constant bit, corr = 0

                avg_corr = corr_sum / 32
                all_bit_correlations[(t1, t2)].append(avg_corr)

        # ── Analysis 3: Joint entropy of (De[17], De[18]) ─────────
        joint_vals = [(s.get(17, 0), s.get(18, 0)) for s in samples]
        joint_ent = entropy_from_samples(joint_vals)
        marginal17 = entropy_from_samples([v[0] for v in joint_vals])
        marginal18 = entropy_from_samples([v[1] for v in joint_vals])
        mutual_info = marginal17 + marginal18 - joint_ent
        all_joint_entropies.append((joint_ent, marginal17, marginal18, mutual_info))

        print(f"  Joint H(De17,De18)={joint_ent:.1f}, H(De17)={marginal17:.1f}, "
              f"H(De18)={marginal18:.1f}, MI={mutual_info:.2f}")

        # ── Analysis 4: Conditional — when De17 is small ───────────
        small_de17 = [(s.get(17, 0), s.get(18, 0)) for s in samples if hw(s.get(17, 0)) <= 8]
        if small_de17:
            cond_hw18 = sum(hw(v[1]) for v in small_de17) / len(small_de17)
            uncond_hw18 = sum(hw(s.get(18, 0)) for s in samples) / len(samples)
            all_conditional.append((cond_hw18, uncond_hw18, len(small_de17)))
            print(f"  Conditional: when HW(De17)≤8 ({len(small_de17)} samples): "
                  f"avg HW(De18)={cond_hw18:.1f} vs unconditional {uncond_hw18:.1f}")

        # ── Analysis 5: HW correlation ────────────────────────────
        hw17 = [hw(s.get(17, 0)) for s in samples]
        hw18 = [hw(s.get(18, 0)) for s in samples]
        n = len(hw17)
        mean17 = sum(hw17) / n
        mean18 = sum(hw18) / n
        cov = sum((a-mean17)*(b-mean18) for a, b in zip(hw17, hw18)) / n
        std17 = (sum((a-mean17)**2 for a in hw17) / n) ** 0.5
        std18 = (sum((b-mean18)**2 for b in hw18) / n) ** 0.5
        if std17 > 0 and std18 > 0:
            hw_corr = cov / (std17 * std18)
        else:
            hw_corr = 0
        all_correlations.append(hw_corr)
        print(f"  HW correlation r(HW(De17), HW(De18)) = {hw_corr:.4f}")
        print()

    # ── Global Summary ────────────────────────────────────────────
    print("=" * 72)
    print("GLOBAL SUMMARY")
    print("=" * 72)
    print()

    # Bit-level correlations
    print("Pairwise bit-level correlations (avg |r| across 32 bits):")
    for (t1, t2), corrs in all_bit_correlations.items():
        if corrs:
            avg = sum(corrs) / len(corrs)
            print(f"  De[{t1}] vs De[{t2}]: avg |r| = {avg:.4f}")

    print()
    # HW correlations
    if all_correlations:
        avg_hw_corr = sum(all_correlations) / len(all_correlations)
        print(f"HW correlation r(HW(De17), HW(De18)): mean = {avg_hw_corr:.4f}")

    # Joint entropy
    if all_joint_entropies:
        print()
        print("Joint entropy analysis:")
        for i, (j, m1, m2, mi) in enumerate(all_joint_entropies):
            print(f"  Msg {i+1}: Joint={j:.1f}, Sum(marginals)={m1+m2:.1f}, MI={mi:.2f}")
        avg_mi = sum(mi for _, _, _, mi in all_joint_entropies) / len(all_joint_entropies)
        print(f"  Average mutual information: {avg_mi:.2f} bits")
        print(f"  (0 = independent, positive = correlated)")

    # Conditional
    if all_conditional:
        print()
        print("Conditional analysis (HW(De18) | HW(De17) ≤ 8):")
        for i, (cond, uncond, n) in enumerate(all_conditional):
            print(f"  Msg {i+1}: conditional={cond:.1f}, unconditional={uncond:.1f} "
                  f"(Δ={cond-uncond:+.1f}, n={n})")
        avg_delta = sum(c-u for c, u, _ in all_conditional) / len(all_conditional)
        print(f"  Average Δ(HW): {avg_delta:+.2f}")

    # ── Verdict ───────────────────────────────────────────────────
    print()
    print("=" * 72)
    print("VERDICT")
    print("=" * 72)
    print()

    if all_correlations:
        avg_corr = abs(sum(all_correlations) / len(all_correlations))
        if all_joint_entropies:
            avg_mi = sum(mi for _, _, _, mi in all_joint_entropies) / len(all_joint_entropies)
        else:
            avg_mi = 0

        if avg_mi > 5:
            print(f"  **ALIVE** — Significant mutual information ({avg_mi:.1f} bits)!")
            print(f"  Barriers are CORRELATED — effective dimension is reduced")
            print(f"  Birthday cost could be 2^{{k/2}} where k < 256")
        elif avg_mi > 1:
            print(f"  **ANOMALY** — Some mutual information ({avg_mi:.1f} bits)")
            print(f"  Slight correlation, but unlikely to help significantly")
        else:
            print(f"  **DEAD** — Negligible mutual information ({avg_mi:.2f} bits)")
            print(f"  Barriers are effectively INDEPENDENT")
            print(f"  Birthday on barriers = birthday on independent 32-bit values")
            print(f"  No shortcut from correlation")

    if all_conditional:
        avg_delta = sum(c-u for c, u, _ in all_conditional) / len(all_conditional)
        print()
        if abs(avg_delta) > 2:
            print(f"  Conditional effect: HW(De18) changes by {avg_delta:+.1f} when De17 is small")
            print(f"  This means barriers interact — could be exploited!")
        else:
            print(f"  No conditional effect: Δ(HW) = {avg_delta:+.2f}")
            print(f"  Low De17 does NOT predict low De18")

    print()
    print(f"  Runtime: {time.time() - t_start:.1f}s")
    print("=" * 72)


if __name__ == "__main__":
    main()
