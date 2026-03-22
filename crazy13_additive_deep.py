#!/usr/bin/env python3
"""
CRAZY-13: Deep Additive Differential Analysis
==============================================

CRAZY-12 found: at R=18, additive diff HW=39 vs XOR HW=54 (15 bits better!)
At R=20-24: additive still ~7-9 bits better.

This is a NEW finding. ALL our previous experiments used XOR differentials.
If the additive landscape has genuine structure, it could shift the attack.

Key questions:
1. Is the additive advantage UNIVERSAL (across messages) or overfitting?
2. What's the expected random min for additive diffs? (different distribution!)
3. Does the advantage grow with GA budget?
4. Can we combine additive + Wang chain + neutral bits?
5. Phase transition for additive: later than round 19?

CRUCIAL INSIGHT: SHA-256 uses MODULAR additions as its primary operation.
Additive differentials are "native" to the algorithm. XOR diffs fight the
addition structure; additive diffs flow WITH it.
"""

import random
import time
import math

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

def sha256_r(msg, iv, nr):
    W = list(msg)
    for i in range(16, max(nr, 16)):
        W.append(add32(sig1(W[i-2]), W[i-7], sig0(W[i-15]), W[i-16]))
    state = list(iv)
    for t in range(nr):
        a, b, c, d, e, f, g, h = state
        T1 = add32(h, Sig1(e), Ch(e, f, g), K[t], W[t])
        T2 = add32(Sig0(a), Maj(a, b, c))
        state = [add32(T1, T2), a, b, c, add32(d, T1), e, f, g]
    return [add32(state[i], iv[i]) for i in range(8)]


def add_diff_metric(msg, delta, iv, nr):
    """Additive differential: M' = M + delta (mod 2^32 per word).
    Returns total 'distance' = sum of min(hw(d), hw(-d)) for each register."""
    msg2 = [(msg[i] + delta[i]) & MASK for i in range(16)]
    h1 = sha256_r(msg, iv, nr)
    h2 = sha256_r(msg2, iv, nr)
    total = 0
    per_reg = []
    for i in range(8):
        d = (h1[i] - h2[i]) & MASK
        m = min(hw(d), hw((-d) & MASK))
        total += m
        per_reg.append(m)
    return total, per_reg


def xor_diff_metric(msg, delta, iv, nr):
    """XOR differential: M' = M ^ delta.
    Returns total HW and per-register."""
    msg2 = [(msg[i] ^ delta[i]) & MASK for i in range(16)]
    h1 = sha256_r(msg, iv, nr)
    h2 = sha256_r(msg2, iv, nr)
    total = 0
    per_reg = []
    for i in range(8):
        d = h1[i] ^ h2[i]
        m = hw(d)
        total += m
        per_reg.append(m)
    return total, per_reg


def ga_search_add(base_msg, iv, nr, pop_size=80, n_gen=200):
    """GA optimizing additive differential."""
    def fitness(delta):
        t, _ = add_diff_metric(base_msg, delta, iv, nr)
        return -t

    pop = []
    for _ in range(pop_size):
        d = [0] * 16
        # Low additive diff: small values
        n_words = random.randint(1, 4)
        for _ in range(n_words):
            w = random.randint(0, 15)
            sign = random.choice([1, -1])
            mag = 1 << random.randint(0, 31)
            d[w] = (d[w] + sign * mag) & MASK
        pop.append(d)

    fits = [fitness(d) for d in pop]
    best_f = max(fits)
    best_d = pop[fits.index(best_f)]

    for g in range(n_gen):
        new_pop = [best_d]
        while len(new_pop) < pop_size:
            idxs = random.sample(range(pop_size), 5)
            p1 = pop[max(idxs, key=lambda i: fits[i])]
            idxs = random.sample(range(pop_size), 5)
            p2 = pop[max(idxs, key=lambda i: fits[i])]
            # Crossover: per-word choice
            child = [(p1[w] if random.random() < 0.5 else p2[w]) for w in range(16)]
            # Mutation: add/subtract small value to random word
            for _ in range(random.randint(1, 3)):
                w = random.randint(0, 15)
                mag = 1 << random.randint(0, 31)
                if random.random() < 0.5:
                    child[w] = (child[w] + mag) & MASK
                else:
                    child[w] = (child[w] - mag) & MASK
            # Ensure non-zero
            if all(x == 0 for x in child):
                child[0] = 1
            new_pop.append(child)
        pop = new_pop
        fits = [fitness(d) for d in pop]
        gf = max(fits)
        if gf > best_f:
            best_f = gf
            best_d = pop[fits.index(gf)]

    return -best_f, best_d


def ga_search_xor(base_msg, iv, nr, pop_size=80, n_gen=200):
    """GA optimizing XOR differential."""
    def fitness(delta):
        t, _ = xor_diff_metric(base_msg, delta, iv, nr)
        return -t

    pop = []
    for _ in range(pop_size):
        d = [0] * 16
        for _ in range(random.randint(1, 5)):
            w = random.randint(0, 15)
            b = random.randint(0, 31)
            d[w] ^= (1 << b)
        pop.append(d)

    fits = [fitness(d) for d in pop]
    best_f = max(fits)
    best_d = pop[fits.index(best_f)]

    for g in range(n_gen):
        new_pop = [best_d]
        while len(new_pop) < pop_size:
            idxs = random.sample(range(pop_size), 5)
            p1 = pop[max(idxs, key=lambda i: fits[i])]
            idxs = random.sample(range(pop_size), 5)
            p2 = pop[max(idxs, key=lambda i: fits[i])]
            child = [p1[w] if random.random() < 0.5 else p2[w] for w in range(16)]
            for _ in range(random.randint(1, 2)):
                w = random.randint(0, 15)
                b = random.randint(0, 31)
                child[w] ^= (1 << b)
            if all(x == 0 for x in child):
                child[0] = 1
            new_pop.append(child)
        pop = new_pop
        fits = [fitness(d) for d in pop]
        gf = max(fits)
        if gf > best_f:
            best_f = gf
            best_d = pop[fits.index(gf)]

    return -best_f, best_d


def main():
    print("=" * 72)
    print("CRAZY-13: Deep Additive Differential Analysis")
    print("=" * 72)
    print()

    random.seed(0xADD13)
    t_start = time.time()
    iv = list(H0)

    NUM_MSGS = 10
    base_msgs = [[random.getrandbits(32) for _ in range(16)] for _ in range(NUM_MSGS)]

    # ── Phase 1: Random baseline — distribution of additive vs XOR ──

    print("Phase 1: Random differential distribution (no GA)")
    print("-" * 72)

    for nr in [18, 20, 24, 64]:
        xor_hws = []
        add_hws = []
        for _ in range(2000):
            delta = [random.getrandbits(32) for _ in range(16)]
            msg = base_msgs[0]
            xh, _ = xor_diff_metric(msg, delta, iv, nr)
            ah, _ = add_diff_metric(msg, delta, iv, nr)
            xor_hws.append(xh)
            add_hws.append(ah)

        xor_avg = sum(xor_hws) / len(xor_hws)
        add_avg = sum(add_hws) / len(add_hws)
        xor_min = min(xor_hws)
        add_min = min(add_hws)
        print(f"  R={nr:2d}: XOR avg={xor_avg:.1f} min={xor_min}, "
              f"ADD avg={add_avg:.1f} min={add_min}, "
              f"Δavg={add_avg-xor_avg:+.1f}, Δmin={add_min-xor_min:+d}")

    # ── Phase 2: GA comparison across multiple messages ───────────

    print()
    print("Phase 2: GA comparison (XOR vs Additive, multiple messages)")
    print("-" * 72)

    for nr in [17, 18, 19, 20, 24]:
        if time.time() - t_start > 300:
            print(f"  [Time limit at R={nr}]")
            break

        xor_results = []
        add_results = []

        for msg_idx in range(min(5, NUM_MSGS)):
            base_msg = base_msgs[msg_idx]
            xor_hw, _ = ga_search_xor(base_msg, iv, nr, pop_size=60, n_gen=150)
            add_hw, _ = ga_search_add(base_msg, iv, nr, pop_size=60, n_gen=150)
            xor_results.append(xor_hw)
            add_results.append(add_hw)

        xor_avg = sum(xor_results) / len(xor_results)
        add_avg = sum(add_results) / len(add_results)
        xor_best = min(xor_results)
        add_best = min(add_results)

        elapsed = time.time() - t_start
        print(f"  R={nr:2d}: XOR avg={xor_avg:.1f} best={xor_best}, "
              f"ADD avg={add_avg:.1f} best={add_best}, "
              f"advantage={xor_avg-add_avg:+.1f} [{elapsed:.0f}s]")

    # ── Phase 3: Single-value additive census ─────────────────────

    if time.time() - t_start < 340:
        print()
        print("Phase 3: Additive single-value census (ΔW[k] = ±1, ±2, ±2^j)")
        print("-" * 72)

        base_msg = base_msgs[0]
        nr = 20

        best_add_single = 999
        best_add_config = None

        for word in range(16):
            for val in [1, -1, 2, -2, 4, -4, 0x80000000, 0x7FFFFFFF]:
                delta = [0] * 16
                delta[word] = val & MASK
                total, _ = add_diff_metric(base_msg, delta, iv, nr)
                if total < best_add_single:
                    best_add_single = total
                    best_add_config = (word, val)

        print(f"  Best single-value additive diff at R=20:")
        print(f"    W[{best_add_config[0]}] + {best_add_config[1] if best_add_config[1] < 0x80000000 else best_add_config[1] - 0x100000000}")
        print(f"    Additive metric: {best_add_single}/256")
        print()

        # Compare: single-bit XOR
        best_xor_single = 999
        best_xor_config = None
        for word in range(16):
            for bit in range(32):
                delta = [0] * 16
                delta[word] = 1 << bit
                total, _ = xor_diff_metric(base_msg, delta, iv, nr)
                if total < best_xor_single:
                    best_xor_single = total
                    best_xor_config = (word, bit)

        print(f"  Best single-bit XOR diff at R=20:")
        print(f"    W[{best_xor_config[0]}] bit {best_xor_config[1]}")
        print(f"    XOR metric: {best_xor_single}/256")

    # ── Phase 4: Additive phase transition ────────────────────────

    if time.time() - t_start < 380:
        print()
        print("Phase 4: Additive differential phase transition")
        print("-" * 72)

        # Expected min for additive metric
        # Additive metric uses min(hw(d), hw(-d)) per register
        # For random d: E[min(hw(d), hw(-d))] ≈ 12.5 for 32-bit
        # (because if hw(d) > 16, we use hw(-d) which is ≈ 16 - (hw(d)-16))
        # Actually min(hw(d), hw(32-hw(d)+ correction))
        # For random d: E[metric per reg] ≈ 12.8 (empirical from Phase 1)
        # E[total] ≈ 8 × 12.8 = 102

        # Expected min: harder to compute analytically
        # From Phase 1, random min at R=64 was about 73-78

        base_msg = base_msgs[0]
        n_evals = 60 * 150

        for nr in range(16, 25):
            if time.time() - t_start > 420:
                break

            add_hw, add_delta = ga_search_add(base_msg, iv, nr, pop_size=60, n_gen=150)
            xor_hw, xor_delta = ga_search_xor(base_msg, iv, nr, pop_size=60, n_gen=150)

            # Per-register for additive
            _, add_regs = add_diff_metric(base_msg, add_delta, iv, nr)
            _, xor_regs = xor_diff_metric(base_msg, xor_delta, iv, nr)

            add_zeros = sum(1 for r in add_regs if r == 0)
            xor_zeros = sum(1 for r in xor_regs if r == 0)

            print(f"  R={nr:2d}: ADD={add_hw:3d} (zeros={add_zeros}), "
                  f"XOR={xor_hw:3d} (zeros={xor_zeros}), "
                  f"ADD advantage={xor_hw-add_hw:+d}")

    # ── Phase 5: Hybrid — XOR+Additive ────────────────────────────

    if time.time() - t_start < 450:
        print()
        print("Phase 5: Hybrid — start with additive, refine with XOR")
        print("-" * 72)

        for nr in [18, 20, 24]:
            if time.time() - t_start > 470:
                break

            base_msg = base_msgs[0]

            # Step 1: Find best additive delta
            add_hw, add_delta = ga_search_add(base_msg, iv, nr, pop_size=80, n_gen=200)

            # Step 2: Convert to XOR: M' = M + delta, so XOR_delta = M ⊕ (M + delta)
            msg2_add = [(base_msg[i] + add_delta[i]) & MASK for i in range(16)]
            xor_equiv = [(base_msg[i] ^ msg2_add[i]) & MASK for i in range(16)]

            # Step 3: Measure this XOR delta
            xor_of_add, _ = xor_diff_metric(base_msg, xor_equiv, iv, nr)

            # Step 4: Use the additive-derived XOR delta as seed for XOR GA
            # (skip for now, just compare)

            # Also: direct XOR GA
            xor_hw, _ = ga_search_xor(base_msg, iv, nr, pop_size=80, n_gen=200)

            print(f"  R={nr:2d}: Direct ADD={add_hw}, "
                  f"ADD→XOR={xor_of_add}, Direct XOR={xor_hw}")

    # ── Verdict ───────────────────────────────────────────────────

    print()
    print("=" * 72)
    print("VERDICT")
    print("=" * 72)
    print()

    print("  KEY INSIGHT: Additive vs XOR differential metrics measure")
    print("  DIFFERENT things. A collision requires H(M) = H(M'), which is")
    print("  the SAME in both metrics (zero = zero regardless of representation).")
    print()
    print("  The question is: does the additive metric provide a better")
    print("  GRADIENT toward zero than the XOR metric?")
    print()
    print("  If additive consistently finds lower metric values, it means")
    print("  the additive landscape is smoother and easier to navigate.")
    print("  This could mean additive-based GA converges to collisions faster.")

    print(f"\n  Runtime: {time.time() - t_start:.1f}s")
    print("=" * 72)


if __name__ == "__main__":
    main()
