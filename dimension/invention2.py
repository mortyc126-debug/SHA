"""
ТРИ НОВЫХ ИЗОБРЕТЕНИЯ.

1. CORRELATED δ: выбрать δW[14], δW[15] так чтобы
   σ1(δW[14]) ⊕ σ1(δW[15]) в schedule ГАСИЛИ друг друга.

2. LOW-BIT COLLISION: H = IV + state (ADD mod 2^32).
   ADD preserves low bits: если state1 ≡ state2 mod 2^k,
   то H1 ≡ H2 mod 2^k. Collision в k LSBs = ЛЕГЧЕ.

3. EVOLUTIONARY: популяция (M, δM) пар, скрещивание + мутация,
   fitness = 256 - HW(δH). Evolve toward collision.
"""

import numpy as np
import struct, hashlib
import time

MASK32 = 0xFFFFFFFF

def hw(x): return bin(x).count('1')
def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def add32(x, y): return (x + y) & MASK32
def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)

def sha256_full(W16):
    raw = struct.pack('>16I', *W16)
    return struct.unpack('>8I', hashlib.sha256(raw).digest())


def main():
    np.random.seed(42)

    print("=" * 70)
    print("ТРИ НОВЫХ ИЗОБРЕТЕНИЯ")
    print("=" * 70)

    # ═══════════════════
    print(f"\n{'='*70}")
    print("1. CORRELATED δ: cancel schedule propagation")
    print(f"{'='*70}")

    # δW[16] = σ1(δW[14])  — from W[14] change
    # δW[17] = σ1(δW[15])  — from W[15] change
    # δW[19] = σ1(δW[17])  — second generation from W[15]
    # δW[21] = σ1(δW[19])  — third generation

    # What if δW[14] is chosen so σ1(δW[14]) CANCELS with
    # something from δW[15] chain?
    # σ1(δW[14]) ⊕ σ1(σ1(δW[15])) = 0? (cancel at W[19])
    # Need: σ1(δW[14]) = σ1(σ1(δW[15]))

    # Search: for each δW[15], find δW[14] that cancels at W[19]
    print(f"\n  Finding δW[14] that cancels δW[15]'s effect at W[19]:")

    best_cancel = None
    min_total_sched_hw = 999

    for b15 in range(32):
        dW15 = 1 << b15
        # δW[17] = σ1(dW15)
        dW17 = sigma1(dW15)
        # δW[19] from W[15] chain: σ1(dW17)
        dW19_from_15 = sigma1(dW17)

        # We want: δW[19] from W[14] to cancel this
        # δW[16] = σ1(dW14), δW[18] = σ1(δW[16]) = σ1²(dW14)
        # But δW[19] depends on δW[17] (from W[15]) AND δW[12] (=0)...
        # Actually: δW[19] = σ1(δW[17]) + δW[12] + σ0(δW[4]) + δW[3]
        # If ONLY W[14,15] change: δW[19] = σ1(δW[17]) = σ1(σ1(dW15))

        # For W[14] to affect W[19]:
        # δW[19] = σ1(δW[17]) where δW[17] gets contribution from...
        # Actually δW[17] = σ1(δW[15]) only. W[14] affects W[16], not W[17].
        # W[14] → W[16] = σ1(W[14]), W[18] = σ1(W[16]), W[20] = σ1(W[18])
        # W[15] → W[17] = σ1(W[15]), W[19] = σ1(W[17]), W[21] = σ1(W[19])

        # They're on DIFFERENT schedule tracks (even vs odd)!
        # W[14] → W[16,18,20,...] (even)
        # W[15] → W[17,19,21,...] (odd)
        # They NEVER directly cancel in schedule!

        # BUT: both affect the STATE at every round.
        # At round 16: δW[16] enters state (from W[14]).
        # At round 17: δW[17] enters state (from W[15]).
        # Can they cancel IN THE STATE?

        for b14 in range(32):
            dW14 = 1 << b14
            dW16 = sigma1(dW14)
            dW17_15 = sigma1(dW15)
            # Total schedule HW for first 8 expanded words
            total = hw(dW16) + hw(dW17_15) + hw(sigma1(dW16)) + hw(sigma1(dW17_15))
            if total < min_total_sched_hw:
                min_total_sched_hw = total
                best_cancel = (b14, b15, total)

    print(f"  Best correlated pair: W[14]bit{best_cancel[0]}, W[15]bit{best_cancel[1]}")
    print(f"  Total schedule HW (4 words): {best_cancel[2]}")

    # Test this pair vs independent
    b14, b15 = best_cancel[0], best_cancel[1]
    N = 20000
    corr_best = 256
    indep_best = 256

    for _ in range(N):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        # Correlated
        Wc = list(W); Wc[14] ^= (1 << b14); Wc[15] ^= (1 << b15)
        H1 = sha256_full(W); Hc = sha256_full(Wc)
        dc = sum(hw(H1[i]^Hc[i]) for i in range(8))
        if dc < corr_best: corr_best = dc

        # Independent (only W[15])
        Wi = list(W); Wi[15] ^= (1 << b15)
        Hi = sha256_full(Wi)
        di = sum(hw(H1[i]^Hi[i]) for i in range(8))
        if di < indep_best: indep_best = di

    print(f"\n  Full SHA-256 (20K search):")
    print(f"    Correlated (W[14]+W[15]): best δH = {corr_best}")
    print(f"    Independent (W[15] only): best δH = {indep_best}")
    print(f"    Advantage: {indep_best - corr_best:+d}")

    # ═══════════════════
    print(f"\n{'='*70}")
    print("2. LOW-BIT COLLISION: match k LSBs of each hash word")
    print(f"{'='*70}")

    # H[i] = (IV[i] + state[i]) mod 2^32
    # Low bits: H1[i] mod 2^k = H2[i] mod 2^k
    # iff (IV[i]+s1[i]) mod 2^k = (IV[i]+s2[i]) mod 2^k
    # iff s1[i] mod 2^k = s2[i] mod 2^k
    # iff (s1[i]-s2[i]) mod 2^k = 0
    # iff δstate[i] mod 2^k = 0 (low k bits of difference = 0)

    # Birthday on k-bit match: 2^(4k) for 8 words with k bits each.
    # k=1: match bit 0 of all 8 words: 2^4 = 16 trials
    # k=2: 2^8 = 256 trials
    # k=4: 2^16 = 65K trials
    # k=8: 2^32 ≈ 4 billion

    print(f"\n  Low-bit matching: find M1,M2 with H1[i] ≡ H2[i] mod 2^k for all i")

    for k in [1, 2, 3, 4]:
        mask_k = (1 << k) - 1
        found = 0
        N_search = 100000

        # Birthday approach: hash table on low bits
        low_table = {}
        for trial in range(N_search):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            H = sha256_full(W)
            key = tuple(H[i] & mask_k for i in range(8))

            if key in low_table:
                W_prev = low_table[key]
                H_prev = sha256_full(W_prev)
                # Verify
                match = all((H[i] & mask_k) == (H_prev[i] & mask_k) for i in range(8))
                if match:
                    full_dH = sum(hw(H[i]^H_prev[i]) for i in range(8))
                    if found < 3:
                        print(f"    k={k}: MATCH at trial {trial}, full δH = {full_dH}")
                    found += 1
            else:
                low_table[key] = W

        expected = N_search**2 / (2 ** (8*k + 1))
        print(f"    k={k}: {found} matches in {N_search} (expected ≈{expected:.0f})")

    # ═══════════════════
    print(f"\n{'='*70}")
    print("3. EVOLUTIONARY SEARCH")
    print(f"{'='*70}")

    # Population of messages. Fitness = 256 - HW(δH).
    # Crossover: combine words from two parents.
    # Mutation: flip random bit.

    POP_SIZE = 200
    GENERATIONS = 500
    ELITE = 20
    δ_MASK = 1 << 31  # fixed δ in W[15]

    # Initialize population
    population = []
    for _ in range(POP_SIZE):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W); W2[15] ^= δ_MASK
        H1 = sha256_full(W); H2 = sha256_full(W2)
        fitness = 256 - sum(hw(H1[i]^H2[i]) for i in range(8))
        population.append((fitness, W))

    best_ever = max(population, key=lambda x: x[0])
    history = [best_ever[0]]

    t0 = time.time()
    for gen in range(GENERATIONS):
        # Sort by fitness
        population.sort(key=lambda x: -x[0])

        # Elite
        new_pop = population[:ELITE]

        # Breed
        while len(new_pop) < POP_SIZE:
            # Tournament selection
            p1 = max(np.random.choice(len(population), 3, replace=False),
                     key=lambda i: population[i][0])
            p2 = max(np.random.choice(len(population), 3, replace=False),
                     key=lambda i: population[i][0])

            # Crossover: mix words from both parents
            child_W = []
            for i in range(16):
                if np.random.random() < 0.5:
                    child_W.append(population[p1][1][i])
                else:
                    child_W.append(population[p2][1][i])

            # Mutation: flip 1-3 random bits (NOT in W[15])
            for _ in range(np.random.randint(1, 4)):
                w = np.random.randint(0, 15)  # not 15
                b = np.random.randint(0, 32)
                child_W[w] ^= (1 << b)

            # Evaluate
            W2 = list(child_W); W2[15] ^= δ_MASK
            H1 = sha256_full(child_W); H2 = sha256_full(W2)
            fit = 256 - sum(hw(H1[i]^H2[i]) for i in range(8))
            new_pop.append((fit, child_W))

        population = new_pop

        current_best = max(population, key=lambda x: x[0])
        if current_best[0] > best_ever[0]:
            best_ever = current_best
        history.append(best_ever[0])

    elapsed = time.time() - t0

    print(f"  Evolutionary search: {POP_SIZE} pop × {GENERATIONS} gen = {POP_SIZE*GENERATIONS} evals")
    print(f"  Time: {elapsed:.1f}s")
    print(f"  Best fitness: {best_ever[0]} (= δH = {256 - best_ever[0]})")
    print(f"  History: gen 0={history[0]}, gen 50={history[50]}, "
          f"gen 200={history[200]}, gen {GENERATIONS}={history[-1]}")

    # Compare with random search (same number of evaluations)
    best_random = 0
    for _ in range(POP_SIZE * GENERATIONS):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W); W2[15] ^= δ_MASK
        H1 = sha256_full(W); H2 = sha256_full(W2)
        fit = 256 - sum(hw(H1[i]^H2[i]) for i in range(8))
        if fit > best_random: best_random = fit

    print(f"  Random (same evals): best fitness = {best_random} (δH = {256-best_random})")
    print(f"  Evolution advantage: {best_ever[0] - best_random:+d} fitness")

    # ═══════════════════
    print(f"\n{'='*70}")
    print("ИТОГ")
    print(f"{'='*70}")

    print(f"""
  1. CORRELATED δ: W[14]+W[15] schedule on DIFFERENT tracks (even/odd).
     Cannot cancel in schedule. In state: advantage = {indep_best - corr_best:+d} bits.

  2. LOW-BIT COLLISION: birthday on k LSBs WORKS!
     k=1: instant (match bit 0 of all words)
     k=2: ~256 trials
     k=4: ~65K trials
     BUT: matching k LSBs ≠ collision. Remaining 32-k bits = random.
     To reach full collision: need k=32 → 2^128 (birthday bound).

  3. EVOLUTIONARY SEARCH: fitness {best_ever[0]} vs random {best_random}.
     {'Evolution BETTER!' if best_ever[0] > best_random else 'Evolution ≈ random.'}
     SHA-256 landscape has NO gradient — fitness landscape is FLAT.
     Evolution = random walk on flat landscape = no advantage.

  FUNDAMENTAL: SHA-256 output is a FLAT landscape in EVERY metric.
  No gradient, no shortcuts, no reducible structure.
  This is EXACTLY what makes it a good hash function.
""")


if __name__ == "__main__":
    main()
