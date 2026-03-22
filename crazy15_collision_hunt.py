#!/usr/bin/env python3
"""
CRAZY-15: SHA-256/17 Collision Hunt
====================================

CRAZY-12 found near-collision at 17 rounds: HW=14 with per-reg [7,1,0,0,5,1,0,0].
That's only 14 bits from a REAL collision!

Strategy:
1. Large GA (pop=300, gen=1000) to find best near-collision
2. Local hill climbing from best GA result
3. Multiple restarts with different base messages
4. If HW gets below 8: exhaustive neighborhood search

If we find HW=0, we have a SHA-256/17 collision — a real result!
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

NR = 17  # Target: 17 rounds

def sha256_17(msg, iv):
    W = list(msg)
    # Need W[16] for round 16
    W.append(add32(sig1(W[14]), W[9], sig0(W[1]), W[0]))
    state = list(iv)
    for t in range(NR):
        a, b, c, d, e, f, g, h = state
        T1 = add32(h, Sig1(e), Ch(e, f, g), K[t], W[t])
        T2 = add32(Sig0(a), Maj(a, b, c))
        state = [add32(T1, T2), a, b, c, add32(d, T1), e, f, g]
    return [add32(state[i], iv[i]) for i in range(8)]


def hash_diff(msg, delta, iv):
    """Returns (total_hw, per_reg_hw, hash_xor)"""
    msg2 = [(msg[i] ^ delta[i]) & MASK for i in range(16)]
    h1 = sha256_17(msg, iv)
    h2 = sha256_17(msg2, iv)
    per_reg = []
    total = 0
    xor_vals = []
    for i in range(8):
        d = h1[i] ^ h2[i]
        h = hw(d)
        per_reg.append(h)
        total += h
        xor_vals.append(d)
    return total, per_reg, xor_vals


def main():
    print("=" * 72)
    print("CRAZY-15: SHA-256/17 Collision Hunt")
    print("=" * 72)
    print()

    t_start = time.time()
    iv = list(H0)

    # ── Phase 1: Large GA with multiple restarts ──────────────────

    POP_SIZE = 200
    N_GEN = 500
    N_RESTARTS = 50

    overall_best_hw = 999
    overall_best_msg = None
    overall_best_delta = None
    overall_best_regs = None

    print(f"Phase 1: GA search ({N_RESTARTS} restarts × {POP_SIZE} pop × {N_GEN} gen)")
    print(f"  Total evaluations: ~{N_RESTARTS * POP_SIZE * N_GEN:,}")
    print("-" * 72)

    for restart in range(N_RESTARTS):
        if time.time() - t_start > 420:
            print(f"\n  [Time limit at restart {restart}]")
            break

        random.seed(restart * 7919 + 12345)
        base_msg = [random.getrandbits(32) for _ in range(16)]

        # Initialize population
        pop = []
        for _ in range(POP_SIZE):
            d = [0] * 16
            for _ in range(random.randint(1, 6)):
                w = random.randint(0, 15)
                b = random.randint(0, 31)
                d[w] ^= (1 << b)
            pop.append(d)

        def fitness(delta):
            t, _, _ = hash_diff(base_msg, delta, iv)
            return -t

        fits = [fitness(d) for d in pop]
        best_f = max(fits)
        best_d = pop[fits.index(best_f)]

        for gen in range(N_GEN):
            new_pop = [best_d]
            while len(new_pop) < POP_SIZE:
                idxs = random.sample(range(POP_SIZE), 5)
                p1 = pop[max(idxs, key=lambda i: fits[i])]
                idxs = random.sample(range(POP_SIZE), 5)
                p2 = pop[max(idxs, key=lambda i: fits[i])]
                child = [p1[w] if random.random() < 0.5 else p2[w] for w in range(16)]
                n_mut = random.choices([1, 2, 3], weights=[0.6, 0.3, 0.1])[0]
                for _ in range(n_mut):
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

        this_hw = -best_f
        if this_hw < overall_best_hw:
            overall_best_hw = this_hw
            overall_best_msg = base_msg
            overall_best_delta = best_d
            _, overall_best_regs, _ = hash_diff(base_msg, best_d, iv)

        if restart % 5 == 0 or this_hw < 15:
            elapsed = time.time() - t_start
            print(f"  Restart {restart:3d}: HW={this_hw:3d}  "
                  f"Overall best: {overall_best_hw}  [{elapsed:.0f}s]")

    # ── Phase 2: Local hill climbing from best ────────────────────

    if overall_best_msg and time.time() - t_start < 440:
        print()
        print(f"Phase 2: Hill climbing from HW={overall_best_hw}")
        print(f"  Per-register: {overall_best_regs}")
        print("-" * 72)

        base_msg = overall_best_msg
        current_delta = list(overall_best_delta)
        current_hw = overall_best_hw
        improved = True
        iterations = 0

        while improved and time.time() - t_start < 470:
            improved = False
            iterations += 1

            # Try flipping each bit
            for w in range(16):
                for b in range(32):
                    trial = list(current_delta)
                    trial[w] ^= (1 << b)
                    if all(x == 0 for x in trial):
                        continue
                    t, regs, _ = hash_diff(base_msg, trial, iv)
                    if t < current_hw:
                        current_hw = t
                        current_delta = trial
                        improved = True
                        print(f"    Iter {iterations}: HW={current_hw}, regs={regs}")
                        if current_hw == 0:
                            break
                if current_hw == 0:
                    break

        print(f"  Hill climb: HW={current_hw} after {iterations} passes")

    # ── Phase 3: Neighborhood search if HW is low ────────────────

    if overall_best_hw <= 10 and time.time() - t_start < 480:
        print()
        print(f"Phase 3: Exhaustive 2-bit neighborhood (HW={overall_best_hw})")
        print("-" * 72)

        base_msg = overall_best_msg
        best_delta = overall_best_delta

        # Try all pairs of bit flips
        found = False
        checked = 0
        for w1 in range(16):
            for b1 in range(32):
                for w2 in range(w1, 16):
                    start_b2 = b1+1 if w2 == w1 else 0
                    for b2 in range(start_b2, 32):
                        trial = list(best_delta)
                        trial[w1] ^= (1 << b1)
                        trial[w2] ^= (1 << b2)
                        if all(x == 0 for x in trial):
                            continue
                        t, regs, _ = hash_diff(base_msg, trial, iv)
                        checked += 1
                        if t < overall_best_hw:
                            overall_best_hw = t
                            best_delta = trial
                            overall_best_regs = regs
                            print(f"    Found HW={t}: {regs}")
                            if t == 0:
                                found = True
                                break
                    if found: break
                if found: break
            if found: break

        print(f"  Checked {checked} 2-bit neighbors")

    # ── Final Report ──────────────────────────────────────────────

    print()
    print("=" * 72)
    print("RESULT")
    print("=" * 72)
    print()
    print(f"  Best near-collision: HW = {overall_best_hw} / 256")
    print(f"  Per-register: {overall_best_regs}")
    if overall_best_msg:
        msg2 = [(overall_best_msg[i] ^ overall_best_delta[i]) & MASK for i in range(16)]
        print(f"  M:  {['%08x'%w for w in overall_best_msg]}")
        print(f"  M': {['%08x'%w for w in msg2]}")
        print(f"  ΔM: {['%08x'%d for d in overall_best_delta]}")

    if overall_best_hw == 0:
        print()
        print("  *** SHA-256/17 COLLISION FOUND! ***")
        h1 = sha256_17(overall_best_msg, iv)
        h2 = sha256_17(msg2, iv)
        print(f"  H(M)  = {['%08x'%h for h in h1]}")
        print(f"  H(M') = {['%08x'%h for h in h2]}")
    elif overall_best_hw <= 5:
        print(f"\n  VERY CLOSE! Only {overall_best_hw} bits from collision.")
        print(f"  With 2^{overall_best_hw} more work, a collision is likely reachable.")
    else:
        print(f"\n  Near-collision but {overall_best_hw} bits still differ.")
        print(f"  Estimated work to find collision: ~2^{overall_best_hw}")

    print(f"\n  Runtime: {time.time() - t_start:.1f}s")
    print("=" * 72)


if __name__ == "__main__":
    main()
