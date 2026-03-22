#!/usr/bin/env python3
"""
CRAZY-12: Phase Transition & Additive Differentials
====================================================

CRAZY-11 found: 16 rounds → HW=7, 20 rounds → HW=95.
Where EXACTLY does the phase transition happen?

Also tests additive (modular) differentials as an alternative metric.

Key experiments:
A) Round-by-round GA: rounds 16..24, one at a time, with high resolution
B) Compare XOR vs additive differential landscapes
C) Rotational differential test
D) Multi-message phase transition (is the transition universal?)
"""

import random
import time

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
    """SHA-256 compression reduced to nr rounds."""
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


def xor_diff_hw(msg, delta, iv, nr):
    """HW of XOR hash difference."""
    msg2 = [(msg[i] ^ delta[i]) & MASK for i in range(16)]
    h1 = sha256_r(msg, iv, nr)
    h2 = sha256_r(msg2, iv, nr)
    return sum(hw(h1[i] ^ h2[i]) for i in range(8))


def add_diff_hw(msg, delta, iv, nr):
    """HW of additive (modular) hash difference."""
    msg2 = [(msg[i] + delta[i]) & MASK for i in range(16)]
    h1 = sha256_r(msg, iv, nr)
    h2 = sha256_r(msg2, iv, nr)
    total = 0
    for i in range(8):
        d = (h1[i] - h2[i]) & MASK
        # "HW" of modular diff = min(hw(d), hw(-d))  (two's complement)
        total += min(hw(d), hw((-d) & MASK))
    return total


def ga_search(fitness_fn, pop_size=80, n_gen=200, msg_init='low_weight'):
    """Generic GA search for optimal ΔM. Returns best fitness and ΔM."""
    pop = []
    for _ in range(pop_size):
        d = [0] * 16
        if msg_init == 'low_weight':
            for _ in range(random.randint(1, 5)):
                w = random.randint(0, 15)
                b = random.randint(0, 31)
                d[w] ^= (1 << b)
        else:
            d = [random.getrandbits(32) for _ in range(16)]
        pop.append(d)

    fits = [fitness_fn(d) for d in pop]
    best_f = max(fits)
    best_d = pop[fits.index(best_f)]

    for g in range(n_gen):
        new_pop = [best_d]
        while len(new_pop) < pop_size:
            idxs = random.sample(range(pop_size), min(5, pop_size))
            p1 = pop[max(idxs, key=lambda i: fits[i])]
            idxs = random.sample(range(pop_size), min(5, pop_size))
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
        fits = [fitness_fn(d) for d in pop]
        gf = max(fits)
        if gf > best_f:
            best_f = gf
            best_d = pop[fits.index(gf)]

    return -best_f, best_d


def expected_min(n_bits, n_trials):
    """Expected minimum of n_trials samples from Binomial(n_bits, 0.5)."""
    import math
    mu = n_bits / 2
    sigma = (n_bits / 4) ** 0.5
    # Approximate: min ≈ μ - σ * sqrt(2 * ln(N))
    return mu - sigma * (2 * math.log(n_trials)) ** 0.5


def main():
    print("=" * 72)
    print("CRAZY-12: Phase Transition & Additive Differentials")
    print("=" * 72)
    print()

    random.seed(0xF12E)
    t_start = time.time()
    iv = list(H0)

    # ── Phase A: High-resolution phase transition ─────────────────

    print("Phase A: Round-by-round GA (rounds 14..24)")
    print("-" * 72)
    print()

    NUM_MSGS = 3
    base_msgs = [[random.getrandbits(32) for _ in range(16)] for _ in range(NUM_MSGS)]

    n_evals = 80 * 200  # GA evaluations
    exp_random_min = expected_min(256, n_evals)

    print(f"  GA budget: {n_evals} evaluations per run")
    print(f"  Expected random min (no structure): {exp_random_min:.1f}/256")
    print()
    print(f"  {'Rounds':>6} {'Best HW':>8} {'Avg best':>9} {'Random min':>10} {'Δ':>6} {'Per-reg':>50}")
    print(f"  {'':>6} {'(best)':>8} {'(3 msgs)':>9} {'expected':>10} {'':>6} {'':>50}")
    print("  " + "-" * 90)

    transition_data = {}

    for nr in range(14, 25):
        if time.time() - t_start > 300:
            print(f"  [Time limit at round {nr}]")
            break

        best_hws = []
        best_regs = None
        best_overall = 999

        for msg_idx, base_msg in enumerate(base_msgs):
            fitness_fn = lambda d, m=base_msg, r=nr: -xor_diff_hw(m, d, iv, r)
            best_hw, best_delta = ga_search(fitness_fn, pop_size=80, n_gen=200)

            best_hws.append(best_hw)
            if best_hw < best_overall:
                best_overall = best_hw
                # Per-register
                msg2 = [(base_msg[i] ^ best_delta[i]) & MASK for i in range(16)]
                h1 = sha256_r(base_msg, iv, nr)
                h2 = sha256_r(msg2, iv, nr)
                best_regs = [hw(h1[i] ^ h2[i]) for i in range(8)]

        avg_best = sum(best_hws) / len(best_hws)
        delta = avg_best - exp_random_min
        zeros = sum(1 for r in best_regs if r == 0) if best_regs else 0
        transition_data[nr] = (best_overall, avg_best, delta)

        elapsed = time.time() - t_start
        print(f"  R={nr:2d}:  {best_overall:6d}   {avg_best:8.1f}  {exp_random_min:8.1f}  "
              f"{delta:+5.1f}  {best_regs}  zeros={zeros}  [{elapsed:.0f}s]")

    # ── Phase B: Additive differential comparison ─────────────────

    if time.time() - t_start < 350:
        print()
        print("Phase B: Additive (modular) vs XOR differentials")
        print("-" * 72)

        base_msg = base_msgs[0]

        for nr in [16, 17, 18, 20, 24, 32, 64]:
            if time.time() - t_start > 400:
                break

            # XOR diff GA
            fitness_xor = lambda d, m=base_msg, r=nr: -xor_diff_hw(m, d, iv, r)
            xor_hw, _ = ga_search(fitness_xor, pop_size=60, n_gen=150)

            # Additive diff GA
            fitness_add = lambda d, m=base_msg, r=nr: -add_diff_hw(m, d, iv, r)
            add_hw, _ = ga_search(fitness_add, pop_size=60, n_gen=150)

            print(f"  R={nr:2d}: XOR best HW = {xor_hw:3d}, "
                  f"Additive best HW = {add_hw:3d}, "
                  f"advantage = {xor_hw - add_hw:+d}")

    # ── Phase C: Rotational differential ──────────────────────────

    if time.time() - t_start < 420:
        print()
        print("Phase C: Rotational differential test")
        print("-" * 72)

        base_msg = base_msgs[0]

        for rot_k in [1, 2, 3, 4, 7, 8, 15, 16]:
            if time.time() - t_start > 450:
                break

            # Rotate each word by k bits
            msg_rot = [rotr(w, rot_k) for w in base_msg]

            # Compare hashes
            h1 = sha256_r(base_msg, iv, 64)
            h2 = sha256_r(msg_rot, iv, 64)
            h1_rot = [rotr(h, rot_k) for h in h1]

            # XOR difference between rotated hash and hash of rotated message
            rot_diff = sum(hw(h1_rot[i] ^ h2[i]) for i in range(8))
            # Direct hash difference
            direct_diff = sum(hw(h1[i] ^ h2[i]) for i in range(8))

            print(f"  ROT-{rot_k:2d}: rotational deviation = {rot_diff}/256, "
                  f"direct diff = {direct_diff}/256")

    # ── Phase D: Multi-message transition universality ────────────

    if time.time() - t_start < 460:
        print()
        print("Phase D: Is the phase transition universal?")
        print("-" * 72)

        for nr in [16, 17, 18]:
            if time.time() - t_start > 480:
                break

            hws_per_msg = []
            for base_msg in base_msgs:
                fitness_fn = lambda d, m=base_msg, r=nr: -xor_diff_hw(m, d, iv, r)
                best_hw, _ = ga_search(fitness_fn, pop_size=60, n_gen=100)
                hws_per_msg.append(best_hw)

            print(f"  R={nr}: per-message HWs = {hws_per_msg}, "
                  f"avg = {sum(hws_per_msg)/len(hws_per_msg):.1f}")

    # ── Verdict ───────────────────────────────────────────────────

    print()
    print("=" * 72)
    print("PHASE TRANSITION ANALYSIS")
    print("=" * 72)
    print()

    if transition_data:
        # Find where Δ crosses zero (structure → random)
        structured = [(nr, d) for nr, (_, _, d) in transition_data.items() if d < -10]
        random_like = [(nr, d) for nr, (_, _, d) in transition_data.items() if abs(d) < 10]

        if structured:
            last_structured = max(nr for nr, _ in structured)
            print(f"  Last round with significant structure: R={last_structured}")
            print(f"    (GA finds solutions {-transition_data[last_structured][2]:.0f} bits better than random)")
        if random_like:
            first_random = min(nr for nr, _ in random_like)
            print(f"  First round that looks random: R={first_random}")

        print()
        print("  Phase transition summary:")
        for nr in sorted(transition_data.keys()):
            best, avg, delta = transition_data[nr]
            status = "STRUCTURED" if delta < -20 else "TRANSITION" if delta < -5 else "RANDOM"
            bar_len = max(0, int(-delta / 2))
            bar = "#" * bar_len
            print(f"    R={nr:2d}: Δ={delta:+6.1f}  [{status:>10}]  {bar}")

    print()
    print("=" * 72)
    print("VERDICT")
    print("=" * 72)
    print()

    if transition_data:
        for nr in sorted(transition_data.keys()):
            _, _, delta = transition_data[nr]
            if delta >= -5:
                print(f"  Phase transition completes at round {nr}.")
                print(f"  For R < {nr}: exploitable structure (GA beats random)")
                print(f"  For R ≥ {nr}: indistinguishable from random oracle")
                print(f"  SHA-256 security margin: {64 - nr} rounds beyond transition")
                break

    print(f"\n  Runtime: {time.time() - t_start:.1f}s")
    print("=" * 72)


if __name__ == "__main__":
    main()
