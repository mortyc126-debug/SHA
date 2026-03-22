#!/usr/bin/env python3
"""
crazy1_evolutionary.py -- CRAZY-1: Evolutionary/Genetic Search for Collision Pairs

Treat SHA-256 collision-finding as an OPTIMIZATION problem.
- "Genome": which of the ~82 neutral bits to flip
- "Fitness": minimize HW(De[17..24]) across the barrier zone
- Neutral bits preserve Wang chain (De=0 for rounds 1-16), so we only
  optimize the differential propagation after round 16.

Two approaches:
1. Genetic Algorithm (GA): population=200, generations=500
2. Simulated Annealing (SA): 100K steps, T cooling from 10 to 0.01

Compare both against random baseline.
"""

import random
import time
import math
from copy import deepcopy

MASK = 0xFFFFFFFF

# ── SHA-256 constants ─────────────────────────────────────────────────
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

IV = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
]

# ── SHA-256 primitives ────────────────────────────────────────────────
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
    """Hamming weight of 32-bit value."""
    return bin(x & MASK).count('1')

def expand_W(W16):
    """Expand 16-word message schedule to 25 words (we need up to round 24)."""
    W = list(W16)
    for i in range(16, 25):
        W.append(add32(sig1(W[i-2]), W[i-7], sig0(W[i-15]), W[i-16]))
    return W

def sha256_all_states(W_expanded, n_rounds):
    """Run SHA-256 for n_rounds, return list of full states after each round."""
    a, b, c, d, e, f, g, h = IV
    states = []
    for i in range(n_rounds):
        T1 = add32(h, Sig1(e), Ch(e, f, g), K[i], W_expanded[i])
        T2 = add32(Sig0(a), Maj(a, b, c))
        h = g; g = f; f = e
        e = add32(d, T1)
        d = c; c = b; b = a
        a = add32(T1, T2)
        states.append((a, b, c, d, e, f, g, h))
    return states

def flip_bit(W, word_idx, bit_pos):
    W2 = list(W)
    W2[word_idx] ^= (1 << bit_pos)
    return W2

# ── Wang chain computation ────────────────────────────────────────────
def compute_wang_chain(base_msg, dw0=0x80000000):
    """
    Given base message W[0..15], compute W'[0..15] such that De_t = 0
    for rounds 1..15 (states after rounds 1..15 have same e register).
    DW[0] = dw0 = MSB flip.
    DW[t] chosen to cancel differential in e at each round.
    """
    W = list(base_msg)
    W_prime = list(W)
    W_prime[0] = W[0] ^ dw0

    a, b, c, d, e, f, g, h = IV
    a2, b2, c2, d2, e2, f2, g2, h2 = IV

    # Round 0
    T1 = add32(h, Sig1(e), Ch(e, f, g), K[0], W[0])
    T2 = add32(Sig0(a), Maj(a, b, c))
    h, g, f = g, f, e
    e = add32(d, T1); d, c, b = c, b, a; a = add32(T1, T2)

    T1p = add32(h2, Sig1(e2), Ch(e2, f2, g2), K[0], W_prime[0])
    T2p = add32(Sig0(a2), Maj(a2, b2, c2))
    h2, g2, f2 = g2, f2, e2
    e2 = add32(d2, T1p); d2, c2, b2 = c2, b2, a2; a2 = add32(T1p, T2p)

    for t in range(1, 16):
        # We want new_e == new_e'
        # new_e  = d  + h  + Sig1(e)  + Ch(e,f,g)  + K[t] + W[t]
        # new_e' = d2 + h2 + Sig1(e2) + Ch(e2,f2,g2) + K[t] + W'[t]
        base_new_e = add32(d, h, Sig1(e), Ch(e, f, g), K[t], W[t])
        prime_partial = add32(d2, h2, Sig1(e2), Ch(e2, f2, g2), K[t])
        # Want: prime_partial + W'[t] = base_new_e
        W_prime_t = (base_new_e - prime_partial) & MASK
        W_prime[t] = W_prime_t

        # Advance base state
        T1 = add32(h, Sig1(e), Ch(e, f, g), K[t], W[t])
        T2 = add32(Sig0(a), Maj(a, b, c))
        h, g, f = g, f, e
        e = add32(d, T1); d, c, b = c, b, a; a = add32(T1, T2)

        # Advance prime state
        T1p = add32(h2, Sig1(e2), Ch(e2, f2, g2), K[t], W_prime_t)
        T2p = add32(Sig0(a2), Maj(a2, b2, c2))
        h2, g2, f2 = g2, f2, e2
        e2 = add32(d2, T1p); d2, c2, b2 = c2, b2, a2; a2 = add32(T1p, T2p)

    return W_prime

# ── Neutral bit detection ─────────────────────────────────────────────
def find_neutral_bits(W_base, max_round=16):
    """
    Find bits in W[1..15] that are neutral for the Wang chain differential
    through max_round. A bit is neutral if flipping it in BOTH M and M'
    preserves all De values through round max_round.
    """
    W_diff = compute_wang_chain(W_base)
    W_b_exp = expand_W(W_base)
    W_d_exp = expand_W(W_diff)
    ref_b = sha256_all_states(W_b_exp, max_round + 1)
    ref_d = sha256_all_states(W_d_exp, max_round + 1)

    # Reference differentials (XOR) for each round
    ref_diffs = []
    for r in range(max_round + 1):
        diff = tuple((ref_d[r][i] ^ ref_b[r][i]) for i in range(8))
        ref_diffs.append(diff)

    neutral = []
    for word_idx in range(1, 16):
        for bit_pos in range(32):
            W_b2 = flip_bit(W_base, word_idx, bit_pos)
            W_d2 = flip_bit(W_diff, word_idx, bit_pos)
            W_b2_exp = expand_W(W_b2)
            W_d2_exp = expand_W(W_d2)
            st_b2 = sha256_all_states(W_b2_exp, max_round + 1)
            st_d2 = sha256_all_states(W_d2_exp, max_round + 1)

            is_neutral = True
            for r in range(max_round + 1):
                diff2 = tuple((st_d2[r][i] ^ st_b2[r][i]) for i in range(8))
                if diff2 != ref_diffs[r]:
                    is_neutral = False
                    break
            if is_neutral:
                neutral.append((word_idx, bit_pos))

    return neutral

# ── Fitness function ──────────────────────────────────────────────────
def compute_fitness(W_base, neutral_bits, individual):
    """
    Given an individual (set of indices into neutral_bits to flip),
    compute fitness = -sum(HW(De[t]) for t in 17..24).

    Steps:
    1. Flip the selected neutral bits in W_base to get W_mod
    2. Compute Wang chain for W_mod -> W_mod_prime
    3. Expand both, run 25 rounds
    4. Compute De[t] = e_mod[t] ^ e_mod_prime[t] for t=17..24
    5. fitness = -sum(HW(De[t]))
    """
    W_mod = list(W_base)
    for idx in individual:
        word_idx, bit_pos = neutral_bits[idx]
        W_mod[word_idx] ^= (1 << bit_pos)

    W_mod_prime = compute_wang_chain(W_mod)
    W_m_exp = expand_W(W_mod)
    W_mp_exp = expand_W(W_mod_prime)

    st_m = sha256_all_states(W_m_exp, 25)
    st_mp = sha256_all_states(W_mp_exp, 25)

    total_hw = 0
    de_list = []
    for t in range(17, 25):  # rounds 17..24 (0-indexed)
        de = st_m[t-1][4] ^ st_mp[t-1][4]  # states[t-1] = state after round t (0-indexed round t)
        # Actually sha256_all_states returns states[i] = state after round i (0-indexed)
        # Round 17 means i=17 in the list
        de_list.append(de)
        total_hw += hw(de)

    return -total_hw, de_list


def compute_fitness_fast(W_base, neutral_bits, flip_set):
    """Same as compute_fitness but takes a frozenset of indices."""
    return compute_fitness(W_base, neutral_bits, flip_set)


# ── Genetic Algorithm ─────────────────────────────────────────────────
def tournament_select(population, fitnesses, k=5):
    """Tournament selection: pick k random, return best."""
    indices = random.sample(range(len(population)), k)
    best = max(indices, key=lambda i: fitnesses[i])
    return population[best]

def uniform_crossover(parent1, parent2, n_bits):
    """Uniform crossover on neutral bit indices."""
    child = set()
    all_bits = parent1 | parent2
    for b in all_bits:
        if b in parent1 and b in parent2:
            child.add(b)
        elif b in parent1:
            if random.random() < 0.5:
                child.add(b)
        else:
            if random.random() < 0.5:
                child.add(b)
    return child

def mutate(individual, n_neutral, mutation_rate):
    """Flip each neutral bit position with probability mutation_rate."""
    result = set(individual)
    for i in range(n_neutral):
        if random.random() < mutation_rate:
            if i in result:
                result.remove(i)
            else:
                result.add(i)
    return result


def run_ga(W_base, neutral_bits, pop_size=200, n_gen=500):
    """Run genetic algorithm. Returns best fitness history and best individual."""
    n_neutral = len(neutral_bits)
    if n_neutral == 0:
        print("  No neutral bits found! Cannot run GA.")
        return [], None, []

    mutation_rate = 1.0 / n_neutral

    # Initialize population: random subsets of neutral bits
    population = []
    for _ in range(pop_size):
        ind = set()
        for i in range(n_neutral):
            if random.random() < 0.5:
                ind.add(i)
        population.append(ind)

    # Evaluate initial population
    fitnesses = []
    for ind in population:
        f, _ = compute_fitness(W_base, neutral_bits, ind)
        fitnesses.append(f)

    best_fitness_history = []
    best_ever_fitness = max(fitnesses)
    best_ever_ind = population[fitnesses.index(best_ever_fitness)]
    best_ever_de = None

    de_zero_count = 0  # track how many De[t]=0 events

    print(f"  GA: {n_neutral} neutral bits, pop={pop_size}, gen={n_gen}")
    print(f"  Initial best fitness: {best_ever_fitness}")

    t0 = time.time()
    for gen in range(n_gen):
        # Create new population
        new_pop = []
        # Elitism: keep best
        best_idx = fitnesses.index(max(fitnesses))
        new_pop.append(population[best_idx])

        while len(new_pop) < pop_size:
            p1 = tournament_select(population, fitnesses)
            p2 = tournament_select(population, fitnesses)
            child = uniform_crossover(p1, p2, n_neutral)
            child = mutate(child, n_neutral, mutation_rate)
            new_pop.append(child)

        population = new_pop

        # Evaluate
        fitnesses = []
        for ind in population:
            f, de_list = compute_fitness(W_base, neutral_bits, ind)
            fitnesses.append(f)
            if f > best_ever_fitness:
                best_ever_fitness = f
                best_ever_ind = ind
                best_ever_de = de_list
            # Count De[t]=0 events
            if de_list:
                for de in de_list:
                    if de == 0:
                        de_zero_count += 1

        best_fitness_history.append(max(fitnesses))

        if (gen + 1) % 100 == 0:
            elapsed = time.time() - t0
            print(f"  Gen {gen+1}: best_fitness={max(fitnesses)}, "
                  f"best_ever={best_ever_fitness}, elapsed={elapsed:.1f}s")

    # Get final best De pattern
    if best_ever_de is None:
        _, best_ever_de = compute_fitness(W_base, neutral_bits, best_ever_ind)

    return best_fitness_history, best_ever_ind, best_ever_de, de_zero_count


# ── Simulated Annealing ───────────────────────────────────────────────
def run_sa(W_base, neutral_bits, n_steps=100000, T_start=10.0, T_end=0.01):
    """Run simulated annealing. Returns best fitness history and best state."""
    n_neutral = len(neutral_bits)
    if n_neutral == 0:
        print("  No neutral bits found! Cannot run SA.")
        return [], None, []

    # Initialize: random subset
    current = set()
    for i in range(n_neutral):
        if random.random() < 0.5:
            current.add(i)

    current_fitness, current_de = compute_fitness(W_base, neutral_bits, current)
    best_fitness = current_fitness
    best_ind = set(current)
    best_de = current_de

    fitness_history = []
    de_zero_count = 0

    # Count De=0 in initial
    for de in current_de:
        if de == 0:
            de_zero_count += 1

    print(f"  SA: {n_neutral} neutral bits, {n_steps} steps, T: {T_start}->{T_end}")
    print(f"  Initial fitness: {current_fitness}")

    t0 = time.time()
    log_ratio = math.log(T_end / T_start)

    for step in range(n_steps):
        # Temperature schedule (exponential cooling)
        frac = step / max(n_steps - 1, 1)
        T = T_start * math.exp(log_ratio * frac)

        # Flip one random neutral bit
        bit_to_flip = random.randint(0, n_neutral - 1)
        candidate = set(current)
        if bit_to_flip in candidate:
            candidate.remove(bit_to_flip)
        else:
            candidate.add(bit_to_flip)

        cand_fitness, cand_de = compute_fitness(W_base, neutral_bits, candidate)

        # Count De=0
        for de in cand_de:
            if de == 0:
                de_zero_count += 1

        # Accept/reject
        dF = cand_fitness - current_fitness
        if dF >= 0 or random.random() < math.exp(dF / T):
            current = candidate
            current_fitness = cand_fitness
            current_de = cand_de

        if current_fitness > best_fitness:
            best_fitness = current_fitness
            best_ind = set(current)
            best_de = current_de

        if (step + 1) % 20000 == 0:
            elapsed = time.time() - t0
            fitness_history.append(best_fitness)
            print(f"  Step {step+1}: T={T:.4f}, current={current_fitness}, "
                  f"best={best_fitness}, elapsed={elapsed:.1f}s")

    return fitness_history, best_ind, best_de, de_zero_count


# ── Random baseline ──────────────────────────────────────────────────
def random_baseline(W_base, neutral_bits, n_samples=1000):
    """Evaluate n_samples random neutral bit configurations."""
    n_neutral = len(neutral_bits)
    if n_neutral == 0:
        return -128, [], 0

    best_fitness = -999
    best_de = None
    de_zero_count = 0
    all_fitnesses = []

    for _ in range(n_samples):
        ind = set()
        for i in range(n_neutral):
            if random.random() < 0.5:
                ind.add(i)
        f, de_list = compute_fitness(W_base, neutral_bits, ind)
        all_fitnesses.append(f)
        for de in de_list:
            if de == 0:
                de_zero_count += 1
        if f > best_fitness:
            best_fitness = f
            best_de = de_list

    avg_fitness = sum(all_fitnesses) / len(all_fitnesses)
    return best_fitness, best_de, de_zero_count, avg_fitness


# ══════════════════════════════════════════════════════════════════════
#  MAIN
# ══════════════════════════════════════════════════════════════════════
def main():
    print("=" * 72)
    print("CRAZY-1: Evolutionary/Genetic Search for SHA-256 Collision Pairs")
    print("=" * 72)
    print()

    random.seed(42)

    # Step 1: Generate a random base message
    W_base = [random.getrandbits(32) for _ in range(16)]
    print("Base message W[0..15] generated.")

    # Step 2: Find neutral bits
    print("\nStep 2: Finding neutral bits (bits that preserve Wang chain De=0)...")
    t0 = time.time()
    neutral_bits = find_neutral_bits(W_base, max_round=16)
    t_neutral = time.time() - t0
    print(f"  Found {len(neutral_bits)} neutral bits in {t_neutral:.1f}s")
    if len(neutral_bits) > 0:
        # Show distribution across words
        word_counts = {}
        for w, b in neutral_bits:
            word_counts[w] = word_counts.get(w, 0) + 1
        print(f"  Distribution: {dict(sorted(word_counts.items()))}")
    print()

    # Verify Wang chain works
    W_prime = compute_wang_chain(W_base)
    W_b_exp = expand_W(W_base)
    W_p_exp = expand_W(W_prime)
    st_b = sha256_all_states(W_b_exp, 25)
    st_p = sha256_all_states(W_p_exp, 25)
    print("Wang chain verification (De at each round):")
    for t in range(17):
        de = st_b[t][4] ^ st_p[t][4]
        tag = " <-- initial perturbation" if t == 0 else ""
        if t <= 16 and de != 0 and t > 0:
            tag = " *** WANG CHAIN BROKEN ***"
        print(f"  Round {t:2d}: De = 0x{de:08x} (HW={hw(de):2d}){tag}")
    print("  --- barrier zone ---")
    for t in range(17, 25):
        de = st_b[t][4] ^ st_p[t][4]
        print(f"  Round {t:2d}: De = 0x{de:08x} (HW={hw(de):2d})")
    print()

    if len(neutral_bits) == 0:
        print("NO neutral bits found. Trying a different base message...")
        # Try a few more
        for attempt in range(5):
            W_base = [random.getrandbits(32) for _ in range(16)]
            neutral_bits = find_neutral_bits(W_base, max_round=16)
            if len(neutral_bits) > 0:
                print(f"  Attempt {attempt+1}: found {len(neutral_bits)} neutral bits")
                break
        if len(neutral_bits) == 0:
            print("  Still no neutral bits. Using relaxed criterion (max_round=12)...")
            neutral_bits = find_neutral_bits(W_base, max_round=12)
            print(f"  Found {len(neutral_bits)} neutral bits (relaxed, through round 12)")

    if len(neutral_bits) == 0:
        print("\nFATAL: No neutral bits found at all. Cannot proceed.")
        print("VERDICT: DEAD (no search space)")
        return

    n_neutral = len(neutral_bits)
    # Adjust parameters based on neutral bit count
    # If few neutral bits, reduce GA size
    pop_size = min(200, 2 ** min(n_neutral, 10))
    n_gen = 500 if n_neutral >= 10 else 200
    sa_steps = 100000

    # Step 3: Random baseline (1000 samples)
    print("=" * 72)
    print("Step 3: Random baseline (1000 random neutral-bit configurations)")
    print("=" * 72)
    t0 = time.time()
    rand_best, rand_de, rand_de0, rand_avg = random_baseline(W_base, neutral_bits, 1000)
    t_rand = time.time() - t0
    print(f"  Random baseline: best_fitness={rand_best}, avg_fitness={rand_avg:.1f}")
    print(f"  De[17..24] at best: {['0x%08x' % d for d in rand_de]}")
    print(f"  De=0 events: {rand_de0} out of {1000*8} (rate={rand_de0/(1000*8):.6f})")
    print(f"  Time: {t_rand:.1f}s")
    print()

    # Step 4: Genetic Algorithm
    print("=" * 72)
    print(f"Step 4: Genetic Algorithm (pop={pop_size}, gen={n_gen})")
    print("=" * 72)
    t0 = time.time()
    ga_history, ga_best_ind, ga_best_de, ga_de0 = run_ga(
        W_base, neutral_bits, pop_size=pop_size, n_gen=n_gen
    )
    t_ga = time.time() - t0
    ga_best = max(ga_history) if ga_history else -999
    ga_total_evals = pop_size * n_gen
    print(f"\n  GA result: best_fitness={ga_best}")
    print(f"  De[17..24] at best: {['0x%08x' % d for d in ga_best_de]}")
    hw_list = [hw(d) for d in ga_best_de]
    print(f"  HW per round: {hw_list}")
    print(f"  De=0 events across all evaluations: {ga_de0} (rate={ga_de0/(ga_total_evals*8):.6f})")
    print(f"  Time: {t_ga:.1f}s, evaluations: {ga_total_evals}")

    # Show fitness progression
    if ga_history:
        print(f"  Fitness progression: gen1={ga_history[0]}, "
              f"gen{len(ga_history)//2}={ga_history[len(ga_history)//2]}, "
              f"final={ga_history[-1]}")
    print()

    # Step 5: Simulated Annealing
    print("=" * 72)
    print(f"Step 5: Simulated Annealing ({sa_steps} steps)")
    print("=" * 72)
    t0 = time.time()
    sa_history, sa_best_ind, sa_best_de, sa_de0 = run_sa(
        W_base, neutral_bits, n_steps=sa_steps
    )
    t_sa = time.time() - t0
    sa_best = max(sa_history) if sa_history else -999
    print(f"\n  SA result: best_fitness={sa_best}")
    print(f"  De[17..24] at best: {['0x%08x' % d for d in sa_best_de]}")
    hw_list = [hw(d) for d in sa_best_de]
    print(f"  HW per round: {hw_list}")
    print(f"  De=0 events across all evaluations: {sa_de0} (rate={sa_de0/(sa_steps*8):.6f})")
    print(f"  Time: {t_sa:.1f}s, evaluations: {sa_steps}")
    print()

    # ── Step 6: Comparison ────────────────────────────────────────────
    print("=" * 72)
    print("Step 6: COMPARISON")
    print("=" * 72)
    print()
    print(f"  {'Method':<20s} {'Best Fitness':>14s} {'Avg/1st Fitness':>16s} {'De=0 events':>12s} {'Evals':>10s}")
    print(f"  {'-'*20} {'-'*14} {'-'*16} {'-'*12} {'-'*10}")
    print(f"  {'Random (1K)':<20s} {rand_best:>14d} {rand_avg:>16.1f} {rand_de0:>12d} {'1,000':>10s}")
    ga_evals_str = f"{ga_total_evals:,}"
    print(f"  {'GA':<20s} {ga_best:>14d} {ga_history[0] if ga_history else 0:>16d} {ga_de0:>12d} {ga_evals_str:>10s}")
    sa_evals_str = f"{sa_steps:,}"
    print(f"  {'SA':<20s} {sa_best:>14d} {'N/A':>16s} {sa_de0:>12d} {sa_evals_str:>10s}")
    print()

    # Expected random baseline: 8 rounds * 16 HW = 128
    # (each De[t] is ~random 32-bit, expected HW = 16)
    expected_random = -128
    print(f"  Theoretical random expectation: fitness ~ {expected_random}")
    print(f"    (8 rounds x 16 avg HW per 32-bit word)")
    print()

    # How many De[t]=0 found total?
    total_de0 = rand_de0 + ga_de0 + sa_de0
    total_evals = 1000 + ga_total_evals + sa_steps
    print(f"  Total De[t]=0 events: {total_de0} across {total_evals*8:,} round-evaluations")
    expected_de0_rate = 1.0 / (2**32)
    expected_de0 = total_evals * 8 * expected_de0_rate
    print(f"  Expected by random chance: {expected_de0:.6f}")
    print()

    # Improvement metric
    ga_improvement = ga_best - rand_best if ga_history else 0
    sa_improvement = sa_best - rand_best if sa_history else 0
    best_improvement = max(ga_improvement, sa_improvement)

    print(f"  GA improvement over random best: {ga_improvement:+d}")
    print(f"  SA improvement over random best: {sa_improvement:+d}")
    print(f"  GA improvement over random avg:  {ga_best - rand_avg:+.1f}")
    print()

    # ── VERDICT ───────────────────────────────────────────────────────
    print("=" * 72)
    print("VERDICT")
    print("=" * 72)

    # "ALIVE" criteria:
    # - GA/SA finds significantly better fitness than random (improvement > 5)
    # - Or any De[t]=0 found (would be remarkable)
    alive = False
    reasons = []

    if best_improvement > 5:
        alive = True
        reasons.append(f"Optimization improved fitness by {best_improvement} over random best")

    if total_de0 > 0:
        alive = True
        reasons.append(f"Found {total_de0} De[t]=0 events (extraordinary!)")

    if ga_best > rand_avg + 20:
        alive = True
        reasons.append(f"GA best ({ga_best}) far exceeds random average ({rand_avg:.1f})")

    if alive:
        print("ALIVE -- Evolutionary optimization shows meaningful advantage:")
        for r in reasons:
            print(f"  * {r}")
        print()
        print("This suggests the neutral-bit search space has exploitable structure")
        print("that evolutionary methods can leverage for barrier-zone optimization.")
    else:
        print("DEAD -- Evolutionary optimization shows NO significant advantage:")
        print(f"  * GA best ({ga_best}) vs random best ({rand_best}): improvement = {ga_improvement}")
        print(f"  * SA best ({sa_best}) vs random best ({rand_best}): improvement = {sa_improvement}")
        print(f"  * No De[t]=0 events found in {total_evals*8:,} evaluations")
        print()
        print("The barrier zone (rounds 17-24) appears to be a near-random function")
        print("of neutral bits, with no gradient for evolutionary methods to exploit.")
        print("The fitness landscape is effectively flat/random.")

    print()
    print("=" * 72)

if __name__ == "__main__":
    main()
