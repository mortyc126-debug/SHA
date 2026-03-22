#!/usr/bin/env python3
"""
CRAZY-11: Direct Differential Search with GA
=============================================

COMPLETELY NEW ANGLE: forget Wang chains, forget neutral bits.
Directly evolve ΔM (512-bit message differential) to minimize
HW(H(M) ⊕ H(M⊕ΔM)).

If the fitness landscape has gradient in the FULL ΔM space,
we can find near-collisions without any structural assumptions.

Sub-experiments:
A) Single-bit census: which bit flip gives lowest HW?
B) Direct GA in full ΔM space
C) GA starting from low-weight differentials
D) Multi-message: does the optimal ΔM depend on M?
E) Partial collision: can we zero out specific registers?

The KEY question: is the landscape random (HW ≈ 128 for everything)
or does it have structure (some ΔM give HW << 128)?
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

def sha256_compress(msg, iv):
    W = list(msg)
    for i in range(16, 64):
        W.append(add32(sig1(W[i-2]), W[i-7], sig0(W[i-15]), W[i-16]))
    state = list(iv)
    for t in range(64):
        a, b, c, d, e, f, g, h = state
        T1 = add32(h, Sig1(e), Ch(e, f, g), K[t], W[t])
        T2 = add32(Sig0(a), Maj(a, b, c))
        state = [add32(T1, T2), a, b, c, add32(d, T1), e, f, g]
    return [add32(state[i], iv[i]) for i in range(8)]


def hash_diff_hw(msg, delta_msg, iv):
    """Compute HW(H(M) ⊕ H(M⊕ΔM))"""
    msg2 = [(msg[i] ^ delta_msg[i]) & MASK for i in range(16)]
    h1 = sha256_compress(msg, iv)
    h2 = sha256_compress(msg2, iv)
    return sum(hw(h1[i] ^ h2[i]) for i in range(8))


def hash_diff_per_reg(msg, delta_msg, iv):
    """Compute per-register HW of hash difference"""
    msg2 = [(msg[i] ^ delta_msg[i]) & MASK for i in range(16)]
    h1 = sha256_compress(msg, iv)
    h2 = sha256_compress(msg2, iv)
    return [hw(h1[i] ^ h2[i]) for i in range(8)]


def main():
    print("=" * 72)
    print("CRAZY-11: Direct Differential Search with GA")
    print("=" * 72)
    print()

    random.seed(0xDA11)
    t_start = time.time()
    iv = list(H0)

    # ── Phase A: Single-bit census ────────────────────────────────

    print("Phase A: Single-bit differential census (512 bits)")
    print("-" * 72)

    NUM_MSGS = 20
    base_msgs = [[random.getrandbits(32) for _ in range(16)] for _ in range(NUM_MSGS)]

    bit_scores = []
    for word in range(16):
        for bit in range(32):
            delta = [0] * 16
            delta[word] = 1 << bit
            total_hw = 0
            for msg in base_msgs:
                total_hw += hash_diff_hw(msg, delta, iv)
            avg = total_hw / NUM_MSGS
            bit_scores.append((word, bit, avg))

    bit_scores.sort(key=lambda x: x[2])
    print(f"  Top 10 easiest single-bit differentials:")
    for i, (w, b, avg) in enumerate(bit_scores[:10]):
        print(f"    #{i+1}: W[{w}] bit {b:2d}  →  avg HW = {avg:.1f}")
    print(f"  Bottom 3:")
    for w, b, avg in bit_scores[-3:]:
        print(f"         W[{w}] bit {b:2d}  →  avg HW = {avg:.1f}")
    print(f"  Overall range: {bit_scores[0][2]:.1f} — {bit_scores[-1][2]:.1f}")
    print(f"  Expected for random: 128.0")
    print()

    # ── Phase B: Direct GA in full ΔM space ──────────────────────

    print("Phase B: GA search for optimal ΔM (fixed message)")
    print("-" * 72)

    base_msg = base_msgs[0]
    base_hash = sha256_compress(base_msg, iv)

    POP_SIZE = 100
    N_GEN = 200
    N_BITS = 512  # Total bits in ΔM

    # Initialize population: low-weight differentials
    population = []
    for _ in range(POP_SIZE):
        delta = [0] * 16
        # Flip 1-5 random bits
        n_flips = random.randint(1, 5)
        for _ in range(n_flips):
            w = random.randint(0, 15)
            b = random.randint(0, 31)
            delta[w] ^= (1 << b)
        population.append(delta)

    def fitness(delta):
        return -hash_diff_hw(base_msg, delta, iv)  # maximize (minimize HW)

    fitnesses = [fitness(d) for d in population]
    best_f = max(fitnesses)
    best_delta = population[fitnesses.index(best_f)]

    history = [(-best_f,)]  # track best HW over generations

    for gen in range(N_GEN):
        if time.time() - t_start > 180:
            break

        new_pop = [best_delta]  # elitism

        while len(new_pop) < POP_SIZE:
            # Tournament selection
            idxs = random.sample(range(POP_SIZE), 5)
            p1 = population[max(idxs, key=lambda i: fitnesses[i])]
            idxs = random.sample(range(POP_SIZE), 5)
            p2 = population[max(idxs, key=lambda i: fitnesses[i])]

            # Crossover: uniform per-word
            child = []
            for w in range(16):
                if random.random() < 0.5:
                    child.append(p1[w])
                else:
                    child.append(p2[w])

            # Mutation: flip 1-3 random bits
            n_mut = random.choices([1, 2, 3], weights=[0.7, 0.2, 0.1])[0]
            for _ in range(n_mut):
                w = random.randint(0, 15)
                b = random.randint(0, 31)
                child[w] ^= (1 << b)

            # Ensure non-zero
            if all(x == 0 for x in child):
                child[0] = 1

            new_pop.append(child)

        population = new_pop
        fitnesses = [fitness(d) for d in population]
        gen_best = max(fitnesses)
        if gen_best > best_f:
            best_f = gen_best
            best_delta = population[fitnesses.index(gen_best)]
        history.append((-best_f,))

        if gen % 50 == 0:
            print(f"  Gen {gen:3d}: best HW = {-best_f}, "
                  f"pop avg = {-sum(fitnesses)/len(fitnesses):.1f}")

    print(f"  Final: best HW = {-best_f}")
    print(f"  Best ΔM: {['%08x'%d for d in best_delta]}")
    print(f"  ΔM Hamming weight: {sum(hw(d) for d in best_delta)}")
    print()

    # Per-register breakdown
    reg_hw = hash_diff_per_reg(base_msg, best_delta, iv)
    print(f"  Per-register: {reg_hw}")
    print(f"  Registers with HW≤5: {sum(1 for h in reg_hw if h <= 5)}")
    print()

    # ── Phase C: GA with different base messages ──────────────────

    print("Phase C: Is the optimal ΔM message-dependent?")
    print("-" * 72)

    # Test best ΔM on other messages
    other_hws = []
    for msg in base_msgs[1:10]:
        h = hash_diff_hw(msg, best_delta, iv)
        other_hws.append(h)
    print(f"  Best ΔM on other messages: {other_hws}")
    print(f"  Average: {sum(other_hws)/len(other_hws):.1f}")
    print(f"  (vs {-best_f} on original message)")
    print()

    # ── Phase D: GA per register (partial collision) ──────────────

    if time.time() - t_start < 240:
        print("Phase D: Can we zero ONE register? (partial collision)")
        print("-" * 72)

        for target_reg in range(8):
            if time.time() - t_start > 300:
                break

            # Fitness = -HW of target register only
            def fitness_reg(delta, reg=target_reg):
                msg2 = [(base_msg[i] ^ delta[i]) & MASK for i in range(16)]
                h1 = sha256_compress(base_msg, iv)
                h2 = sha256_compress(msg2, iv)
                return -hw(h1[reg] ^ h2[reg])

            # Quick GA
            pop = []
            for _ in range(80):
                d = [0] * 16
                for _ in range(random.randint(1, 8)):
                    w = random.randint(0, 15)
                    b = random.randint(0, 31)
                    d[w] ^= (1 << b)
                pop.append(d)

            fits = [fitness_reg(d) for d in pop]
            best_reg_f = max(fits)
            best_reg_d = pop[fits.index(best_reg_f)]

            for g in range(100):
                new_pop = [best_reg_d]
                while len(new_pop) < 80:
                    idxs = random.sample(range(80), 5)
                    p1 = pop[max(idxs, key=lambda i: fits[i])]
                    idxs = random.sample(range(80), 5)
                    p2 = pop[max(idxs, key=lambda i: fits[i])]
                    child = [p1[w] if random.random() < 0.5 else p2[w] for w in range(16)]
                    for _ in range(random.randint(1, 2)):
                        w = random.randint(0, 15)
                        b = random.randint(0, 31)
                        child[w] ^= (1 << b)
                    if all(x == 0 for x in child): child[0] = 1
                    new_pop.append(child)
                pop = new_pop
                fits = [fitness_reg(d) for d in pop]
                gf = max(fits)
                if gf > best_reg_f:
                    best_reg_f = gf
                    best_reg_d = pop[fits.index(gf)]

            regnames = ['a','b','c','d','e','f','g','h']
            print(f"  Register {regnames[target_reg]} (#{target_reg}): "
                  f"best HW = {-best_reg_f}/32, "
                  f"expected random: 16")

    # ── Phase E: Multi-message GA (message-independent ΔM) ───────

    if time.time() - t_start < 350:
        print()
        print("Phase E: Message-independent ΔM search")
        print("-" * 72)

        # Fitness = average HW over multiple messages
        eval_msgs = base_msgs[:5]

        def fitness_multi(delta):
            total = 0
            for msg in eval_msgs:
                total += hash_diff_hw(msg, delta, iv)
            return -total / len(eval_msgs)

        pop = []
        for _ in range(60):
            d = [0] * 16
            for _ in range(random.randint(1, 5)):
                w = random.randint(0, 15)
                b = random.randint(0, 31)
                d[w] ^= (1 << b)
            pop.append(d)

        fits = [fitness_multi(d) for d in pop]
        best_mf = max(fits)
        best_md = pop[fits.index(best_mf)]

        for g in range(150):
            if time.time() - t_start > 400:
                break
            new_pop = [best_md]
            while len(new_pop) < 60:
                idxs = random.sample(range(60), 5)
                p1 = pop[max(idxs, key=lambda i: fits[i])]
                idxs = random.sample(range(60), 5)
                p2 = pop[max(idxs, key=lambda i: fits[i])]
                child = [p1[w] if random.random() < 0.5 else p2[w] for w in range(16)]
                for _ in range(random.randint(1, 2)):
                    w = random.randint(0, 15)
                    b = random.randint(0, 31)
                    child[w] ^= (1 << b)
                if all(x == 0 for x in child): child[0] = 1
                new_pop.append(child)
            pop = new_pop
            fits = [fitness_multi(d) for d in pop]
            gf = max(fits)
            if gf > best_mf:
                best_mf = gf
                best_md = pop[fits.index(gf)]

            if g % 50 == 0:
                print(f"  Gen {g:3d}: best avg HW = {-best_mf:.1f}")

        print(f"  Final: best avg HW = {-best_mf:.1f}")
        # Test on held-out messages
        holdout_hws = [hash_diff_hw(msg, best_md, iv) for msg in base_msgs[5:]]
        if holdout_hws:
            print(f"  Held-out messages: avg HW = {sum(holdout_hws)/len(holdout_hws):.1f}")

    # ── Phase F: Reduced rounds with direct GA ────────────────────

    if time.time() - t_start < 420:
        print()
        print("Phase F: Direct GA on reduced-round SHA-256")
        print("-" * 72)

        def sha256_compress_r(msg, iv, nr):
            W = list(msg)
            for i in range(16, 64):
                W.append(add32(sig1(W[i-2]), W[i-7], sig0(W[i-15]), W[i-16]))
            state = list(iv)
            for t in range(nr):
                a, b, c, d, e, f, g, h = state
                T1 = add32(h, Sig1(e), Ch(e, f, g), K[t], W[t])
                T2 = add32(Sig0(a), Maj(a, b, c))
                state = [add32(T1, T2), a, b, c, add32(d, T1), e, f, g]
            return [add32(state[i], iv[i]) for i in range(8)]

        for nr in [16, 20, 24, 28, 32]:
            if time.time() - t_start > 450:
                break

            def fitness_r(delta, n=nr):
                m2 = [(base_msg[i] ^ delta[i]) & MASK for i in range(16)]
                h1 = sha256_compress_r(base_msg, iv, n)
                h2 = sha256_compress_r(m2, iv, n)
                return -sum(hw(h1[i] ^ h2[i]) for i in range(8))

            pop_r = []
            for _ in range(80):
                d = [0] * 16
                for _ in range(random.randint(1, 5)):
                    w = random.randint(0, 15)
                    b = random.randint(0, 31)
                    d[w] ^= (1 << b)
                pop_r.append(d)

            fits_r = [fitness_r(d) for d in pop_r]
            best_rf = max(fits_r)
            best_rd = pop_r[fits_r.index(best_rf)]

            for g in range(200):
                if time.time() - t_start > 450:
                    break
                new_pop = [best_rd]
                while len(new_pop) < 80:
                    idxs = random.sample(range(80), 5)
                    p1 = pop_r[max(idxs, key=lambda i: fits_r[i])]
                    idxs = random.sample(range(80), 5)
                    p2 = pop_r[max(idxs, key=lambda i: fits_r[i])]
                    child = [p1[w] if random.random() < 0.5 else p2[w] for w in range(16)]
                    for _ in range(random.randint(1, 2)):
                        w = random.randint(0, 15)
                        b = random.randint(0, 31)
                        child[w] ^= (1 << b)
                    if all(x == 0 for x in child): child[0] = 1
                    new_pop.append(child)
                pop_r = new_pop
                fits_r = [fitness_r(d) for d in pop_r]
                gf = max(fits_r)
                if gf > best_rf:
                    best_rf = gf
                    best_rd = pop_r[fits_r.index(gf)]

            reg_hw = []
            m2 = [(base_msg[i] ^ best_rd[i]) & MASK for i in range(16)]
            h1 = sha256_compress_r(base_msg, iv, nr)
            h2 = sha256_compress_r(m2, iv, nr)
            for i in range(8):
                reg_hw.append(hw(h1[i] ^ h2[i]))

            print(f"  {nr:2d} rounds: best total HW = {-best_rf:3d}/256, "
                  f"per-reg: {reg_hw}, zeros: {sum(1 for h in reg_hw if h == 0)}")

    # ── Verdict ───────────────────────────────────────────────────

    print()
    print("=" * 72)
    print("VERDICT")
    print("=" * 72)
    print()

    best_single = bit_scores[0]
    print(f"  Best single-bit: W[{best_single[0]}] bit {best_single[1]}, "
          f"avg HW = {best_single[2]:.1f}")
    print(f"  Best GA (fixed msg): HW = {-best_f}")
    if other_hws:
        print(f"  Same ΔM on other msgs: avg HW = {sum(other_hws)/len(other_hws):.1f}")

    if -best_f < 110:
        print(f"\n  **ALIVE** — GA finds ΔM with HW = {-best_f} < 110!")
        print(f"  Landscape has STRONG gradient in ΔM space")
    elif -best_f < 120:
        print(f"\n  **ANOMALY** — GA finds HW = {-best_f}, slightly below random")
    else:
        print(f"\n  **DEAD** — GA finds HW = {-best_f} ≈ 128 (random)")
        print(f"  No gradient in ΔM space")

    print(f"\n  Runtime: {time.time() - t_start:.1f}s")
    print("=" * 72)


if __name__ == "__main__":
    main()
