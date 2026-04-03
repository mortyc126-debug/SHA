"""
V7 Final Assault on R=6: three parallel strategies.

Dir 1: Brute search — V3 with 10K steps + random restarts
Dir 2: Population evolution — crossover stays in Q-variety
Dir 3: Multi-product reduction — exploit product structure

Uses V3 (best for R=5) as foundation.
"""

import random
import time
import math

from qt_solver.sha256_traced import (
    MASK32, IV256, sha256_compress, sha256_compress_traced,
    get_all_carry_chains, get_bit,
)
from qt_solver.gf2 import gf2_solve, gf2_rank
from qt_solver.qt_system import build_qt_system
from qt_solver.v3_solver import _iv_constant_value, build_v3_system


def _extract_msg(bv, vm):
    msg = []
    for w in range(16):
        word = 0
        for b in range(32):
            if (bv >> vm.w(w, b)) & 1:
                word |= (1 << b)
        msg.append(word)
    return msg


def _evaluate(bv, vm, n, ref_carries, num_rounds, remaining_quad):
    m = _extract_msg(bv, vm)
    t = sha256_compress_traced(m, num_rounds)
    c = get_all_carry_chains(t)
    cmm = sum(bin(ref_carries[i] ^ c[i]).count('1') for i in range(len(ref_carries)))
    assignment = {v: (bv >> v) & 1 for v in range(n)}
    qv = 0
    for lin_terms, quad_terms, const in remaining_quad:
        val = const
        for v in lin_terms:
            val ^= assignment.get(v, 0)
        for v1, v2 in quad_terms:
            val ^= (assignment.get(v1, 0) & assignment.get(v2, 0))
        if val != 0:
            qv += 1
    return cmm + qv * 5, cmm, qv


# ═══════════ Dir 1: Heavy Search with Restarts ═══════════

def dir1_heavy_search(num_rounds=6, max_restarts=20, steps_per=2000,
                       seed=42, verbose=True):
    rng = random.Random(seed)

    for restart in range(max_restarts):
        msg = [rng.randint(0, MASK32) for _ in range(16)]
        lr, rq, system, trace, info = build_v3_system(msg, num_rounds)
        vm = system.vm
        n = info['n']
        target = trace.hash_words
        ref_carries = get_all_carry_chains(trace)

        result = gf2_solve(lr, n)
        if result is None:
            continue
        particular, kernel = result

        current = particular
        cost, cmm, qviol = _evaluate(current, vm, n, ref_carries, num_rounds, rq)
        best_cost, best_bv = cost, current
        temperature = 100.0

        for step in range(steps_per):
            # SA with aggressive exploration
            trials = min(80, len(kernel))
            indices = rng.sample(range(len(kernel)), trials)

            for ki in indices:
                kvec = kernel[ki]
                if kvec & ((1 << 512) - 1) == 0:
                    continue
                cand = current ^ kvec
                c, cm, qv = _evaluate(cand, vm, n, ref_carries, num_rounds, rq)
                if c < cost:
                    cost, cmm, qviol = c, cm, qv
                    current = cand
                    if c < best_cost:
                        best_cost = c
                        best_bv = cand
                elif temperature > 0.1:
                    delta = c - cost
                    if delta < temperature and rng.random() < math.exp(-delta / temperature):
                        cost, cmm, qviol = c, cm, qv
                        current = cand

            temperature *= 0.998

            if best_cost == 0:
                found_msg = _extract_msg(best_bv, vm)
                found_hash = sha256_compress(found_msg, num_rounds)
                if found_hash == target and found_msg != msg:
                    if verbose:
                        print(f"  Dir1: PREIMAGE at restart={restart}, step={step}")
                    return {'success': True, 'restart': restart, 'step': step}

        if verbose and restart % 5 == 0:
            print(f"  Dir1: restart={restart}, best_cost={best_cost}")

    if verbose:
        print(f"  Dir1: no solution in {max_restarts} restarts")
    return {'success': False, 'best_cost': best_cost}


# ═══════════ Dir 2: Population Evolution ═══════════

def dir2_population(num_rounds=6, pop_size=200, generations=500,
                     seed=42, verbose=True):
    rng = random.Random(seed)
    msg = [rng.randint(0, MASK32) for _ in range(16)]
    lr, rq, system, trace, info = build_v3_system(msg, num_rounds)
    vm = system.vm
    n = info['n']
    target = trace.hash_words
    ref_carries = get_all_carry_chains(trace)

    result = gf2_solve(lr, n)
    if result is None:
        return {'success': False}
    particular, kernel = result
    kd = len(kernel)

    # Initialize population
    pop = []
    for _ in range(pop_size):
        p = particular
        for kvec in kernel:
            if rng.getrandbits(1):
                p ^= kvec
        c, cm, qv = _evaluate(p, vm, n, ref_carries, num_rounds, rq)
        pop.append((c, p))

    pop.sort()
    best_ever = pop[0][0]

    for gen in range(generations):
        # Selection: top 25%
        elite = pop[:pop_size // 4]

        # Offspring
        children = []
        while len(children) < pop_size - len(elite):
            p1 = elite[rng.randint(0, len(elite) - 1)][1]
            p2 = elite[rng.randint(0, len(elite) - 1)][1]

            # Crossover: p1 ⊕ p2 ⊕ particular (stays in variety)
            child = p1 ^ p2 ^ particular

            # Mutation: flip 1-5 random kernel vectors
            for _ in range(rng.randint(1, 5)):
                child ^= kernel[rng.randint(0, kd - 1)]

            c, cm, qv = _evaluate(child, vm, n, ref_carries, num_rounds, rq)
            children.append((c, child))

            if c == 0:
                found_msg = _extract_msg(child, vm)
                found_hash = sha256_compress(found_msg, num_rounds)
                if found_hash == target and found_msg != msg:
                    if verbose:
                        print(f"  Dir2: PREIMAGE at gen={gen}")
                    return {'success': True, 'generation': gen}

        pop = elite + children
        pop.sort()

        if pop[0][0] < best_ever:
            best_ever = pop[0][0]

        if verbose and gen % 100 == 0:
            print(f"  Dir2: gen={gen}, best={pop[0][0]}, "
                  f"mean={sum(c for c,_ in pop)/len(pop):.0f}")

    if verbose:
        print(f"  Dir2: best={best_ever} after {generations} gens")
    return {'success': False, 'best_cost': best_ever}


# ═══════════ Dir 3: Multi-product reduction ═══════════

def dir3_reduce_multiproduct(num_rounds=6, seed=42, verbose=True):
    """
    Analyze multi-product equations: can some products be pre-determined?

    For each multi-product eq with k products:
      L(x) ⊕ p1 ⊕ p2 ⊕ ... ⊕ pk = 0

    If we can determine (k-1) products, the last one becomes linear.
    Strategy: find products shared between equations → fixing one
    satisfies constraints in multiple equations.
    """
    rng = random.Random(seed)
    msg = [rng.randint(0, MASK32) for _ in range(16)]
    lr, rq, system, trace, info = build_v3_system(msg, num_rounds)
    n = info['n']

    if verbose:
        print(f"\n  Dir3: Multi-product analysis for R={num_rounds}")
        print(f"  {len(rq)} multi-product equations")

    # Collect all products across all multi-product eqs
    from collections import Counter, defaultdict
    product_freq = Counter()
    eq_products = []
    for lin_terms, quad_terms, const in rq:
        prods = list(quad_terms)
        eq_products.append(prods)
        for p in prods:
            product_freq[p] += 1

    total_products = sum(len(p) for p in eq_products)
    unique_products = len(product_freq)
    shared = sum(1 for p, c in product_freq.items() if c > 1)

    if verbose:
        print(f"  Total products: {total_products}")
        print(f"  Unique products: {unique_products}")
        print(f"  Shared (in >1 eq): {shared}")
        prod_per_eq = [len(p) for p in eq_products]
        print(f"  Products/eq: min={min(prod_per_eq)}, max={max(prod_per_eq)}, "
              f"mean={sum(prod_per_eq)/len(prod_per_eq):.1f}")

        freq_dist = Counter(product_freq.values())
        print(f"  Frequency distribution: {dict(freq_dist)}")

    # Key metric: if we could fix ALL shared products,
    # how many equations become single-product?
    if shared > 0:
        would_simplify = 0
        for prods in eq_products:
            non_shared = sum(1 for p in prods if product_freq[p] == 1)
            if non_shared <= 1:
                would_simplify += 1
        if verbose:
            print(f"  Eqs that would become single-product: {would_simplify}/{len(rq)}")

    return {
        'total_products': total_products,
        'unique_products': unique_products,
        'shared_products': shared,
    }


# ═══════════ Run all three ═══════════

def run_all_r6(seed=42):
    print(f"\n{'═'*60}")
    print(f"V7 FINAL ASSAULT ON R=6")
    print(f"{'═'*60}")

    print("\n--- Dir 3: Multi-product analysis ---")
    dir3_reduce_multiproduct(6, seed)

    print("\n--- Dir 2: Population evolution ---")
    r2 = dir2_population(6, pop_size=150, generations=300, seed=seed)

    print("\n--- Dir 1: Heavy search (20 restarts × 2K steps) ---")
    r1 = dir1_heavy_search(6, max_restarts=20, steps_per=2000, seed=seed)

    return r1, r2
