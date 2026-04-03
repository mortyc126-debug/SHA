"""
Improved Q∩T solvers exploiting the key insight:
the ORIGINAL message lies in the Q-variety with T-mismatch = 0.

Instead of random sampling (mismatch ~48-270 bits), we start
from the known-good point and navigate the variety intelligently.

Strategy 1: Greedy Kernel Walk
  From M (mismatch=0), try each kernel vector k_i.
  M' = M ⊕ k_i is still in Q-variety.
  Keep moves that minimize T-mismatch.

Strategy 2: Neutral Bit Analysis
  Find kernel vectors that PRESERVE carry consistency.
  These are the "neutral bits" from methodology_v20.md.

Strategy 3: Population Search
  Crossover of two Q-variety points stays in the variety.
  Evolve a population toward T-consistency.
"""

import random
import time
from collections import defaultdict

from qt_solver.sha256_traced import (
    MASK32, sha256_compress_traced, sha256_compress,
    get_all_carry_chains,
)
from qt_solver.gf2 import gf2_solve, gf2_rank
from qt_solver.qt_system import QTSystem, VarMap, build_qt_system


def _extract_msg(bitvec, vm):
    """Extract 16 message words from a bit vector."""
    msg = []
    for w in range(16):
        word = 0
        for b in range(32):
            if (bitvec >> vm.w(w, b)) & 1:
                word |= (1 << b)
        msg.append(word)
    return msg


def _msg_to_bitvec(msg, vm):
    """Convert message words to bit vector."""
    bv = 0
    for w in range(16):
        for b in range(32):
            if (msg[w] >> b) & 1:
                bv |= (1 << vm.w(w, b))
    return bv


def _compute_mismatch(msg, ref_carries, num_rounds):
    """Compute carry mismatch between msg's actual carries and reference."""
    trace = sha256_compress_traced(msg, num_rounds)
    actual = get_all_carry_chains(trace)
    total = 0
    for i in range(len(ref_carries)):
        total += bin(ref_carries[i] ^ actual[i]).count('1')
    return total, trace


def _compute_combined_cost(bitvec, vm, ref_carries, num_rounds, quad_eqs):
    """
    Combined fitness: carry mismatch + quadratic violations.

    IMPORTANT: quad equations are evaluated using the Q-PREDICTED state
    (from bitvec), NOT the actual SHA trace. This ensures that cost=0
    truly means all Q equations are satisfied, guaranteeing hash match.
    """
    msg = _extract_msg(bitvec, vm)
    carry_mm, trace = _compute_mismatch(msg, ref_carries, num_rounds)

    # Build assignment from bitvec (Q-predicted values)
    assignment = {}
    for v in range(vm.num_vars):
        assignment[v] = (bitvec >> v) & 1

    quad_viol = sum(1 for eq in quad_eqs if eq.evaluate(assignment) != 0)
    return carry_mm + quad_viol * 5, carry_mm, quad_viol


# ════════════════════════════════════════════════════════════
# Strategy 1: Greedy Kernel Walk
# ════════════════════════════════════════════════════════════

def greedy_kernel_walk(num_rounds, max_steps=500, seed=42, verbose=True):
    """
    Start from M with mismatch=0, greedily walk the Q-variety
    by flipping kernel vectors to minimize T-mismatch for a
    DIFFERENT target hash.

    For collision: find M' in Q-variety(M) where SHA(M') = SHA(M) but M'≠M.
    For preimage: find M' where SHA(M') = target.

    Key insight: at each step, we try ALL kernel vectors and pick the
    best one. This is O(kernel_dim) per step but each step is cheap.
    """
    rng = random.Random(seed)

    if verbose:
        print(f"\n{'='*60}")
        print(f"Greedy Kernel Walk: R={num_rounds}")
        print(f"{'='*60}")

    # Generate base message and target
    base_msg = [rng.randint(0, MASK32) for _ in range(16)]
    system, trace = build_qt_system(base_msg, num_rounds)
    target = trace.hash_words
    ref_carries = get_all_carry_chains(trace)
    vm = system.vm

    # Solve linear system
    rows = system.to_linear_rows()
    result = gf2_solve(rows, vm.num_vars)
    if result is None:
        print("ERROR: Linear system inconsistent")
        return {}

    particular, kernel = result

    if verbose:
        print(f"Kernel dimension: {len(kernel)}")

    # Verify base message is a solution
    base_bv = _msg_to_bitvec(base_msg, vm)

    # Current point = particular solution (may differ from base_msg in state vars)
    # We need the full solution that matches base_msg
    # Since base_msg generated the carries, it should satisfy Q
    # Let's start from particular and find the kernel combo matching base_msg

    # Actually, start from particular and walk
    current = particular
    current_msg = _extract_msg(current, vm)
    current_mismatch, _ = _compute_mismatch(current_msg, ref_carries, num_rounds)

    if verbose:
        print(f"Starting mismatch: {current_mismatch}")

    best_mismatch = current_mismatch
    best_msg = current_msg
    history = [current_mismatch]

    # Pre-compute message-only parts of kernel vectors for fast evaluation
    msg_kernel = []
    for kvec in kernel:
        msg_part = kvec & ((1 << 512) - 1)
        msg_kernel.append(msg_part)

    t0 = time.time()
    no_improve_count = 0

    for step in range(max_steps):
        # Try each kernel vector
        best_next = current_mismatch
        best_kvec_idx = -1

        for ki, kvec in enumerate(kernel):
            # Quick check: does this vector change message bits?
            if msg_kernel[ki] == 0:
                continue  # Only changes state/schedule vars, not message

            candidate = current ^ kvec
            candidate_msg = _extract_msg(candidate, vm)
            mm, _ = _compute_mismatch(candidate_msg, ref_carries, num_rounds)

            if mm < best_next:
                best_next = mm
                best_kvec_idx = ki

        if best_kvec_idx >= 0 and best_next < current_mismatch:
            current ^= kernel[best_kvec_idx]
            current_msg = _extract_msg(current, vm)
            current_mismatch = best_next
            no_improve_count = 0

            if current_mismatch < best_mismatch:
                best_mismatch = current_mismatch
                best_msg = list(current_msg)
        else:
            no_improve_count += 1
            # Random perturbation: flip random kernel vectors
            for _ in range(rng.randint(1, min(5, len(kernel)))):
                ki = rng.randint(0, len(kernel) - 1)
                current ^= kernel[ki]
            current_msg = _extract_msg(current, vm)
            current_mismatch, _ = _compute_mismatch(current_msg, ref_carries, num_rounds)

        history.append(current_mismatch)

        if verbose and step % 50 == 0:
            elapsed = time.time() - t0
            print(f"  step {step}: current={current_mismatch}, "
                  f"best={best_mismatch}, no_improve={no_improve_count}, "
                  f"{elapsed:.1f}s")

        if best_mismatch == 0:
            # Check hash
            h = sha256_compress(best_msg, num_rounds)
            if h == target and best_msg != base_msg:
                if verbose:
                    print(f"  *** COLLISION FOUND at step {step}! ***")
                return {'success': True, 'msg': best_msg, 'steps': step}
            elif h == target:
                if verbose:
                    print(f"  step {step}: mismatch=0 but same message")

    elapsed = time.time() - t0
    if verbose:
        print(f"\nResult: best_mismatch={best_mismatch} ({elapsed:.1f}s)")
        print(f"History (first 20): {history[:20]}")

    return {
        'best_mismatch': best_mismatch,
        'history': history,
        'steps': max_steps,
    }


# ════════════════════════════════════════════════════════════
# Strategy 2: Neutral Bit Analysis
# ════════════════════════════════════════════════════════════

def neutral_bit_analysis(num_rounds, num_bases=5, seed=42, verbose=True):
    """
    Find kernel vectors that preserve carry consistency.

    For each kernel vector k_i, measure T-mismatch(M ⊕ k_i).
    Vectors with mismatch=0 are "carry-neutral" — they can be
    flipped without disrupting the carry chains.

    These correspond to the "neutral bits" from methodology_v20.md
    (P-99: 103/512 neutral positions, |Sol_Wang| >= 2^104).
    """
    rng = random.Random(seed)

    if verbose:
        print(f"\n{'='*60}")
        print(f"Neutral Bit Analysis: R={num_rounds}")
        print(f"{'='*60}")

    all_results = []

    for base_idx in range(num_bases):
        base_msg = [rng.randint(0, MASK32) for _ in range(16)]
        system, trace = build_qt_system(base_msg, num_rounds)
        ref_carries = get_all_carry_chains(trace)
        vm = system.vm

        rows = system.to_linear_rows()
        result = gf2_solve(rows, vm.num_vars)
        if result is None:
            continue

        particular, kernel = result

        # Test each kernel vector
        mismatches = []
        neutral_count = 0
        low_mismatch_count = 0

        for ki, kvec in enumerate(kernel):
            msg_part = kvec & ((1 << 512) - 1)
            if msg_part == 0:
                # Doesn't change message — automatically neutral
                mismatches.append(0)
                neutral_count += 1
                continue

            candidate = particular ^ kvec
            candidate_msg = _extract_msg(candidate, vm)
            mm, _ = _compute_mismatch(candidate_msg, ref_carries, num_rounds)
            mismatches.append(mm)

            if mm == 0:
                neutral_count += 1
            if mm <= 5:
                low_mismatch_count += 1

        all_results.append({
            'kernel_dim': len(kernel),
            'neutral': neutral_count,
            'low_mismatch': low_mismatch_count,
            'mismatches': mismatches,
        })

        if verbose:
            msg_changing = sum(1 for k in kernel if k & ((1 << 512) - 1) != 0)
            print(f"\n  Base {base_idx}: kernel={len(kernel)}, "
                  f"msg-changing={msg_changing}")
            print(f"  Neutral (mm=0): {neutral_count}/{len(kernel)}")
            print(f"  Low mismatch (mm<=5): {low_mismatch_count}/{len(kernel)}")
            if mismatches:
                nonzero = [m for m in mismatches if m > 0]
                if nonzero:
                    print(f"  Non-neutral: min={min(nonzero)}, "
                          f"mean={sum(nonzero)/len(nonzero):.1f}, "
                          f"max={max(nonzero)}")

                # Distribution
                buckets = defaultdict(int)
                for m in mismatches:
                    buckets[m] += 1
                top = sorted(buckets.items(), key=lambda x: -x[1])[:8]
                print(f"  Top mismatch values: "
                      f"{', '.join(f'{k}:{v}' for k,v in top)}")

    return all_results


# ════════════════════════════════════════════════════════════
# Strategy 3: Population Search on Q-Variety
# ════════════════════════════════════════════════════════════

def population_search(num_rounds, pop_size=100, generations=200,
                      seed=42, verbose=True):
    """
    Evolutionary search on Q-variety.

    Key property: if p1, p2 are in affine variety x0 + span{K},
    then p1 ⊕ p2 ⊕ x0 is also in the variety (crossover).

    - Population: points in Q-variety
    - Fitness: -T_mismatch (lower is better)
    - Crossover: p1 ⊕ p2 ⊕ x0 (stays in variety)
    - Mutation: flip random kernel vectors
    - Selection: tournament
    """
    rng = random.Random(seed)

    if verbose:
        print(f"\n{'='*60}")
        print(f"Population Search: R={num_rounds}")
        print(f"{'='*60}")

    base_msg = [rng.randint(0, MASK32) for _ in range(16)]
    system, trace = build_qt_system(base_msg, num_rounds)
    ref_carries = get_all_carry_chains(trace)
    target = trace.hash_words
    vm = system.vm

    rows = system.to_linear_rows()
    result = gf2_solve(rows, vm.num_vars)
    if result is None:
        return {}

    particular, kernel = result
    kernel_dim = len(kernel)

    if verbose:
        print(f"Kernel dim: {kernel_dim}, Pop: {pop_size}, Gens: {generations}")

    # Initialize population
    population = []
    for _ in range(pop_size):
        point = particular
        for kvec in kernel:
            if rng.getrandbits(1):
                point ^= kvec
        msg = _extract_msg(point, vm)
        mm, _ = _compute_mismatch(msg, ref_carries, num_rounds)
        population.append((point, mm))

    t0 = time.time()
    best_ever = min(mm for _, mm in population)
    history = [best_ever]

    for gen in range(generations):
        # Sort by fitness
        population.sort(key=lambda x: x[1])
        best_gen = population[0][1]

        if best_gen < best_ever:
            best_ever = best_gen

        if best_ever == 0:
            if verbose:
                print(f"  Gen {gen}: mismatch=0 FOUND!")
            break

        # Selection: keep top half
        survivors = population[:pop_size // 2]

        # Generate offspring
        offspring = []
        while len(offspring) < pop_size - len(survivors):
            # Tournament selection
            p1 = survivors[rng.randint(0, len(survivors) - 1)][0]
            p2 = survivors[rng.randint(0, len(survivors) - 1)][0]

            # Crossover: p1 ⊕ p2 ⊕ particular (stays in variety)
            child = p1 ^ p2 ^ particular

            # Mutation: flip 1-3 random kernel vectors
            for _ in range(rng.randint(1, 3)):
                child ^= kernel[rng.randint(0, kernel_dim - 1)]

            child_msg = _extract_msg(child, vm)
            mm, _ = _compute_mismatch(child_msg, ref_carries, num_rounds)
            offspring.append((child, mm))

        population = survivors + offspring
        history.append(best_ever)

        if verbose and gen % 50 == 0:
            fitnesses = [mm for _, mm in population]
            print(f"  Gen {gen}: best={best_gen}, best_ever={best_ever}, "
                  f"mean={sum(fitnesses)/len(fitnesses):.1f}, "
                  f"min_pop={min(fitnesses)}")

    elapsed = time.time() - t0

    if verbose:
        print(f"\nResult: best_ever={best_ever} ({elapsed:.1f}s)")
        print(f"Improvement: random_baseline→{history[0]}, final→{best_ever}")

    return {
        'best_mismatch': best_ever,
        'history': history,
        'elapsed': elapsed,
    }


# ════════════════════════════════════════════════════════════
# Combined benchmark
# ════════════════════════════════════════════════════════════

def greedy_combined_walk(num_rounds, max_steps=300, seed=42, verbose=True):
    """
    Greedy walk minimizing COMBINED cost:
      cost = carry_mismatch + 5 * quad_violations

    This ensures both T-consistency (carry) and Q-consistency (Ch/Maj).
    When cost=0: carries match AND all quadratic eqs satisfied
    → SHA(M) = target GUARANTEED.
    """
    rng = random.Random(seed)

    if verbose:
        print(f"\n{'='*60}")
        print(f"Greedy Combined Walk: R={num_rounds}")
        print(f"{'='*60}")

    base_msg = [rng.randint(0, MASK32) for _ in range(16)]
    system, trace = build_qt_system(base_msg, num_rounds)
    target = trace.hash_words
    ref_carries = get_all_carry_chains(trace)
    vm = system.vm

    rows = system.to_linear_rows()
    result = gf2_solve(rows, vm.num_vars)
    if result is None:
        return {}
    particular, kernel = result
    quad_eqs = system.q_quadratic

    if verbose:
        print(f"Kernel: {len(kernel)}, Quad eqs: {len(quad_eqs)}")

    current = particular
    cost, cmm, qviol = _compute_combined_cost(
        current, vm, ref_carries, num_rounds, quad_eqs
    )
    best_cost, best_cmm, best_qviol = cost, cmm, qviol
    best_bv = current

    if verbose:
        print(f"Start: cost={cost} (carry={cmm}, quad_viol={qviol})")

    history = [(cmm, qviol)]
    t0 = time.time()

    for step in range(max_steps):
        improved = False
        for ki, kvec in enumerate(kernel):
            if kvec & ((1 << 512) - 1) == 0:
                continue
            candidate = current ^ kvec
            c, cm, qv = _compute_combined_cost(
                candidate, vm, ref_carries, num_rounds, quad_eqs
            )
            if c < cost:
                cost, cmm, qviol = c, cm, qv
                current = candidate
                improved = True
                if c < best_cost:
                    best_cost, best_cmm, best_qviol = c, cm, qv
                    best_bv = current
                break  # First improvement (faster than best improvement)

        if not improved:
            # Random perturbation
            for _ in range(rng.randint(1, 5)):
                current ^= kernel[rng.randint(0, len(kernel)-1)]
            cost, cmm, qviol = _compute_combined_cost(
                current, vm, ref_carries, num_rounds, quad_eqs
            )

        history.append((cmm, qviol))

        if verbose and step % 30 == 0:
            print(f"  step {step}: cost={cost} (carry={cmm}, quad={qviol}), "
                  f"best=({best_cmm},{best_qviol}), {time.time()-t0:.1f}s")

        if best_cost == 0:
            # VERIFY
            final_msg = _extract_msg(best_bv, vm)
            final_hash = sha256_compress(final_msg, num_rounds)
            if verbose:
                match = final_hash == target
                same = final_msg == base_msg
                print(f"\n  cost=0 at step {step}!")
                print(f"  Hash match: {match}, Same msg: {same}")
                if match and not same:
                    print(f"  *** SECOND PREIMAGE FOUND! ***")
            return {
                'success': final_hash == target,
                'best_cost': 0,
                'steps': step,
                'msg': final_msg,
                'history': history,
            }

    elapsed = time.time() - t0
    if verbose:
        print(f"\nBest: carry={best_cmm}, quad_viol={best_qviol} ({elapsed:.1f}s)")
    return {
        'success': False,
        'best_cost': best_cost,
        'best_carry_mm': best_cmm,
        'best_quad_viol': best_qviol,
        'history': history,
    }


def benchmark_all(num_rounds=2, seed=42):
    """Run all strategies and compare."""
    print(f"\n{'═'*60}")
    print(f"BENCHMARK: R={num_rounds}")
    print(f"{'═'*60}")

    # Baseline: random sampling
    print("\n[Baseline: Random Sampling]")
    rng = random.Random(seed)
    base_msg = [rng.randint(0, MASK32) for _ in range(16)]
    system, trace = build_qt_system(base_msg, num_rounds)
    ref_carries = get_all_carry_chains(trace)
    vm = system.vm
    rows = system.to_linear_rows()
    result = gf2_solve(rows, vm.num_vars)
    if result:
        particular, kernel = result
        mismatches = []
        for _ in range(500):
            sample = particular
            for kvec in kernel:
                if rng.getrandbits(1):
                    sample ^= kvec
            msg = _extract_msg(sample, vm)
            mm, _ = _compute_mismatch(msg, ref_carries, num_rounds)
            mismatches.append(mm)
        print(f"  500 samples: min={min(mismatches)}, "
              f"mean={sum(mismatches)/len(mismatches):.1f}, "
              f"max={max(mismatches)}")

    # Strategy 1
    greedy_kernel_walk(num_rounds, max_steps=200, seed=seed)

    # Strategy 2
    neutral_bit_analysis(num_rounds, num_bases=3, seed=seed)

    # Strategy 3
    population_search(num_rounds, pop_size=80, generations=150, seed=seed)
