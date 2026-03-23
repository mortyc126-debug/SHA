#!/usr/bin/env python3
"""
CRAZY-6: Hybrid GA+Birthday Complexity Estimation
==================================================
Estimate the REAL complexity of a hybrid attack combining all known techniques:

Phase 1: Measure optimized barrier costs for rounds 17-24 using neutral bits + mini-GA
Phase 2: Track full 8-register state differential after optimization
Phase 3: Compute total attack complexity budget

Combines: Wang chain (r1-16) + free words (r17-20) + neutral bits + GA + birthday
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

def rotr(x, n):
    return ((x >> n) | (x << (32 - n))) & MASK

def Sig0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sig1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def sig0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sig1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)
def Ch(e, f, g):  return ((e & f) ^ (~e & g)) & MASK
def Maj(a, b, c): return ((a & b) ^ (a & c) ^ (b & c)) & MASK

def add32(*args):
    s = 0
    for a in args:
        s = (s + a) & MASK
    return s

def hw(x):
    return bin(x & MASK).count('1')

def expand_W(W16):
    W = list(W16)
    for i in range(16, 64):
        W.append(add32(sig1(W[i-2]), W[i-7], sig0(W[i-15]), W[i-16]))
    return W

def sha256_round(state, W_t, K_t):
    a, b, c, d, e, f, g, h = state
    T1 = add32(h, Sig1(e), Ch(e, f, g), K_t, W_t)
    T2 = add32(Sig0(a), Maj(a, b, c))
    return [add32(T1, T2), a, b, c, add32(d, T1), e, f, g]

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

# ── Fitness: minimize total HW of De at target rounds ──────────────

def evaluate_config(base_msg, init_state, neutrals, flip_indices, max_round=25):
    """Flip specified neutral bits, run SHA-256, return per-round De and full state diffs."""
    W = list(base_msg)
    W_prime = compute_wang_chain(base_msg, init_state)
    for idx in flip_indices:
        k, b = neutrals[idx]
        flip = 1 << b
        W[k] ^= flip
        W_prime[k] ^= flip

    # Recompute wang chain for modified message
    W_prime = compute_wang_chain(W, init_state)
    W_exp = expand_W(W)
    Wp_exp = expand_W(W_prime)

    s1 = run_rounds(init_state, W_exp, max_round)
    s2 = run_rounds(init_state, Wp_exp, max_round)

    de_per_round = {}
    state_diff_hw = {}
    for t in range(17, min(max_round + 1, len(s1))):
        de = (s1[t][4] ^ s2[t][4]) & MASK
        de_per_round[t] = de
        total_hw = sum(hw((s1[t][r] ^ s2[t][r]) & MASK) for r in range(8))
        state_diff_hw[t] = total_hw

    return de_per_round, state_diff_hw

def mini_ga(base_msg, init_state, neutrals, target_rounds, pop_size=50, n_gen=80):
    """Mini-GA optimizing sum(HW(De[t])) for target_rounds."""
    n = len(neutrals)
    if n == 0:
        return set(), float('inf')

    # Initialize population
    population = []
    for _ in range(pop_size):
        ind = set()
        for i in range(n):
            if random.random() < 0.5:
                ind.add(i)
        population.append(ind)

    def fitness(ind):
        de_map, _ = evaluate_config(base_msg, init_state, neutrals, ind, max_round=max(target_rounds) + 1)
        total = 0
        for t in target_rounds:
            if t in de_map:
                total += hw(de_map[t])
        return -total  # maximize (minimize HW)

    fitnesses = [fitness(ind) for ind in population]
    best_f = max(fitnesses)
    best_ind = population[fitnesses.index(best_f)]
    mut_rate = 1.0 / n

    for gen in range(n_gen):
        new_pop = [best_ind]  # elitism
        while len(new_pop) < pop_size:
            # tournament
            idxs = random.sample(range(pop_size), min(5, pop_size))
            p1 = population[max(idxs, key=lambda i: fitnesses[i])]
            idxs = random.sample(range(pop_size), min(5, pop_size))
            p2 = population[max(idxs, key=lambda i: fitnesses[i])]
            # crossover
            child = set()
            for bit in p1 | p2:
                if bit in p1 and bit in p2:
                    child.add(bit)
                elif random.random() < 0.5:
                    child.add(bit)
            # mutation
            for i in range(n):
                if random.random() < mut_rate:
                    child.symmetric_difference_update({i})
            new_pop.append(child)

        population = new_pop
        fitnesses = [fitness(ind) for ind in population]
        gen_best = max(fitnesses)
        if gen_best > best_f:
            best_f = gen_best
            best_ind = population[fitnesses.index(gen_best)]

    return best_ind, best_f


def main():
    print("=" * 72)
    print("CRAZY-6: Hybrid GA+Birthday Complexity Estimation")
    print("=" * 72)
    print()

    random.seed(0xFACE)
    t_start = time.time()

    NUM_MSGS = 20  # reduced for time
    init_state = list(H0)

    # ── Phase 1: Per-round barrier cost with optimization ──────────

    print("Phase 1: Measuring optimized barrier costs (rounds 17-24)")
    print("-" * 72)

    # Track per-round statistics
    round_stats = {t: {
        'baseline_de_hw': [],
        'optimized_de_hw': [],
        'baseline_state_hw': [],
        'optimized_state_hw': [],
        'de_zero_count': 0,
    } for t in range(17, 25)}

    for msg_idx in range(NUM_MSGS):
        elapsed = time.time() - t_start
        if elapsed > 420:  # 7 min limit
            print(f"\n[Time limit at message {msg_idx}]")
            break

        if msg_idx % 5 == 0:
            print(f"  Message {msg_idx}/{NUM_MSGS} (elapsed: {elapsed:.1f}s)")

        base_msg = [random.getrandbits(32) for _ in range(16)]

        # Baseline: no optimization
        de_base, state_hw_base = evaluate_config(base_msg, init_state, [], set(), max_round=25)
        for t in range(17, 25):
            if t in de_base:
                round_stats[t]['baseline_de_hw'].append(hw(de_base[t]))
                round_stats[t]['baseline_state_hw'].append(state_hw_base[t])

        # Find neutral bits
        neutrals = find_neutral_bits(base_msg, init_state)
        if len(neutrals) == 0:
            continue

        # Mini-GA optimization for each target round
        for target_round in range(17, 25):
            if time.time() - t_start > 420:
                break

            best_ind, best_f = mini_ga(
                base_msg, init_state, neutrals,
                target_rounds=[target_round],
                pop_size=30, n_gen=50
            )

            de_opt, state_hw_opt = evaluate_config(
                base_msg, init_state, neutrals, best_ind,
                max_round=target_round + 1
            )

            if target_round in de_opt:
                de_hw = hw(de_opt[target_round])
                round_stats[target_round]['optimized_de_hw'].append(de_hw)
                round_stats[target_round]['optimized_state_hw'].append(state_hw_opt[target_round])
                if de_opt[target_round] == 0:
                    round_stats[target_round]['de_zero_count'] += 1

    # ── Phase 1 Report ─────────────────────────────────────────────

    print()
    print("=" * 72)
    print("Phase 1 Results: Per-Round Barrier Cost")
    print("=" * 72)
    print()
    print(f"{'Round':>6} {'Baseline':>10} {'Optimized':>10} {'Savings':>8} "
          f"{'Base State':>11} {'Opt State':>10} {'De=0':>6}")
    print(f"{'':>6} {'HW(De)':>10} {'HW(De)':>10} {'(bits)':>8} "
          f"{'HW(all)':>11} {'HW(all)':>10} {'count':>6}")
    print("-" * 72)

    barrier_costs = {}
    optimized_state_hws = {}

    for t in range(17, 25):
        stats = round_stats[t]
        if stats['baseline_de_hw'] and stats['optimized_de_hw']:
            base_avg = sum(stats['baseline_de_hw']) / len(stats['baseline_de_hw'])
            opt_avg = sum(stats['optimized_de_hw']) / len(stats['optimized_de_hw'])
            savings = base_avg - opt_avg
            base_state = sum(stats['baseline_state_hw']) / len(stats['baseline_state_hw'])
            opt_state = sum(stats['optimized_state_hw']) / len(stats['optimized_state_hw'])

            barrier_costs[t] = opt_avg
            optimized_state_hws[t] = opt_state

            print(f"  r{t:2d}  {base_avg:10.1f} {opt_avg:10.1f} {savings:8.1f} "
                  f"{base_state:11.1f} {opt_state:10.1f} {stats['de_zero_count']:6d}")
        else:
            print(f"  r{t:2d}  (insufficient data)")

    # ── Phase 2: Multi-round joint optimization ────────────────────

    print()
    print("=" * 72)
    print("Phase 2: Joint Multi-Round GA Optimization (rounds 17-20)")
    print("=" * 72)

    joint_de_hw = {t: [] for t in range(17, 21)}
    joint_state_hw = {t: [] for t in range(17, 21)}

    for msg_idx in range(min(10, NUM_MSGS)):
        if time.time() - t_start > 420:
            break

        base_msg = [random.getrandbits(32) for _ in range(16)]
        neutrals = find_neutral_bits(base_msg, init_state)
        if len(neutrals) == 0:
            continue

        best_ind, best_f = mini_ga(
            base_msg, init_state, neutrals,
            target_rounds=[17, 18, 19, 20],
            pop_size=40, n_gen=80
        )

        de_opt, state_hw_opt = evaluate_config(
            base_msg, init_state, neutrals, best_ind,
            max_round=21
        )

        for t in range(17, 21):
            if t in de_opt:
                joint_de_hw[t].append(hw(de_opt[t]))
                joint_state_hw[t].append(state_hw_opt[t])

    print()
    print(f"{'Round':>6} {'Joint HW(De)':>13} {'Joint State HW':>15}")
    print("-" * 40)
    for t in range(17, 21):
        if joint_de_hw[t]:
            avg_de = sum(joint_de_hw[t]) / len(joint_de_hw[t])
            avg_st = sum(joint_state_hw[t]) / len(joint_state_hw[t])
            print(f"  r{t:2d}  {avg_de:13.1f} {avg_st:15.1f}")

    # ── Phase 3: Total Complexity Budget ───────────────────────────

    print()
    print("=" * 72)
    print("Phase 3: Total Attack Complexity Budget")
    print("=" * 72)
    print()

    # Model:
    # Rounds 1-16: Wang chain, cost O(1) = 2^0
    # Rounds 17-20: free words available, cost per barrier = 2^(opt_HW_De)
    #   (finding De[t]=0 requires trying ~2^HW candidates)
    # Rounds 21-24: neutral bits only (no free words), cost = 2^(opt_HW_De)
    # Rounds 25-64: no control, birthday on full state diff
    #   cost = 2^(state_HW / 2) per round via birthday
    #   but actually we need ALL registers to match at round 64,
    #   so it's more like 2^128 for the remaining block

    total_log2 = 0
    print("  Optimistic estimate (each barrier independent):")
    print()

    # Rounds 1-16
    print(f"    Rounds  1-16: Wang chain           cost = 2^0")

    # Rounds 17-24 with optimization
    for t in range(17, 25):
        if t in barrier_costs:
            cost = barrier_costs[t]
            # Cost to find De[t]=0 ≈ 2^HW(De) (birthday on 32-bit register)
            # But we measure OPTIMIZED HW, so actual cost ≈ 2^(opt_HW)
            # More precisely: P(De=0) ≈ 2^(-HW) since each bit independently ~1/2
            # But De is a function of W, not independent bits
            # Conservative: cost = 2^32 per barrier (generic)
            # Optimistic: cost = 2^(32 - savings) where savings from GA
            baseline_hw = 16.0  # expected HW of random 32-bit
            savings_bits = baseline_hw - cost
            barrier_log2 = max(32 - savings_bits, 16)  # at least 2^16
            total_log2 += barrier_log2
            print(f"    Round  {t:2d}:   GA-optimized barrier  cost = 2^{barrier_log2:.1f} "
                  f"(HW: {cost:.1f}, saved {savings_bits:.1f} bits)")
        else:
            total_log2 += 32
            print(f"    Round  {t:2d}:   unoptimized barrier   cost = 2^32")

    # Rounds 25-64: 40 barriers, no optimization possible
    unopt_barriers = 40
    unopt_cost = unopt_barriers * 32
    total_log2 += unopt_cost
    print(f"    Rounds 25-64: {unopt_barriers} unoptimized barriers   cost = {unopt_barriers} × 2^32 = 2^{unopt_cost}")

    print()
    print(f"    TOTAL (sum of barriers, optimistic): 2^{total_log2:.1f}")
    print()

    # Birthday bound comparison
    # The above treats each barrier as independent — but they're NOT.
    # The REAL cost is min(product of barriers, birthday bound)
    # Birthday on 256-bit state = 2^128
    # Birthday on 32-bit e register per round = 2^16
    # But barriers are sequential, not parallel

    # Realistic model: Wang chain + free words cover r1-20
    # After r20: 44 barriers × 2^32 each = 2^(44*32) = 2^1408 (sequential)
    # vs birthday attack: 2^128 (parallel on full hash)
    # Birthday wins massively

    print("  Realistic comparison:")
    print()
    print(f"    Optimistic sequential: 2^{total_log2:.0f}")
    print(f"    Birthday bound:       2^128")
    print(f"    Brute force:          2^256")
    print()

    # What the optimization ACTUALLY buys:
    # The savings from neutral bits + GA only apply to rounds 17-24
    # Savings: sum(savings_bits) over 8 rounds
    total_savings = 0
    for t in range(17, 25):
        if t in barrier_costs:
            total_savings += 16.0 - barrier_costs[t]

    print(f"    Total savings from optimization (rounds 17-24): {total_savings:.1f} bits")
    print(f"    Effective complexity: 2^(128 - {total_savings:.1f}) = 2^{128 - total_savings:.1f}")
    print(f"    (assuming optimization reduces birthday bound proportionally)")
    print()

    # ── Verdict ────────────────────────────────────────────────────

    effective = 128 - total_savings
    print("=" * 72)
    print("VERDICT")
    print("=" * 72)
    print()

    if effective < 120:
        verdict = "ALIVE"
        print(f"  **ALIVE** — Effective complexity 2^{effective:.1f} < 2^120")
        print(f"  GA+neutral bits save {total_savings:.1f} bits overall")
    elif effective < 126:
        verdict = "ANOMALY"
        print(f"  **ANOMALY** — Effective complexity 2^{effective:.1f} (savings: {total_savings:.1f} bits)")
        print(f"  Some reduction but not enough for practical impact")
    else:
        verdict = "DEAD"
        print(f"  **DEAD** — Effective complexity 2^{effective:.1f} ≈ 2^128")
        print(f"  Total savings only {total_savings:.1f} bits — negligible vs birthday bound")

    print()
    print("  Key insight: Even with ALL optimization techniques combined,")
    print("  the attack reduces to birthday bound 2^128 because:")
    print("  1. Optimization only helps rounds 17-24 (8 out of 64)")
    print("  2. Rounds 25-64 have no free variables to optimize")
    print("  3. Sequential barrier approach (2^1408) loses to birthday (2^128)")
    print("  4. Birthday attack doesn't benefit from per-round optimization")
    print()

    elapsed = time.time() - t_start
    print(f"  Total runtime: {elapsed:.1f}s")
    print("=" * 72)


if __name__ == "__main__":
    main()
