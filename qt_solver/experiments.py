"""
Experiments for the Q∩T solver framework.

Validates correctness and measures structural properties
of the Q∩T decomposition at different round counts.
"""

import random
import time

from qt_solver.sha256_traced import (
    sha256_compress, sha256_compress_traced,
    get_all_carry_chains, get_carry_out_vector,
    count_additions, MASK32, IV256,
)
from qt_solver.gf2 import gf2_rank, gf2_solve, popcount
from qt_solver.qt_system import QTSystem, VarMap, build_qt_system, GF2Equation


def verify_sha256_traced(num_tests=100, num_rounds=64):
    """Verify that traced SHA-256 produces correct hashes."""
    rng = random.Random(42)
    errors = 0
    for _ in range(num_tests):
        msg = [rng.randint(0, MASK32) for _ in range(16)]
        h_std = sha256_compress(msg, num_rounds)
        trace = sha256_compress_traced(msg, num_rounds)
        if h_std != trace.hash_words:
            errors += 1
    return errors, num_tests


def experiment_carry_statistics(num_rounds=64, num_samples=100, seed=42):
    """
    Measure carry statistics for SHA-256.

    From methodology_v20.md:
      - 7 additions per round × 64 rounds = 448 compression additions
      - 3 additions per schedule word × 48 words = 144 schedule additions
      - 8 final IV additions
      - Total: 600 additions, 592 carry chains
      - GF(2) rank of carry-out: 589/592
    """
    rng = random.Random(seed)
    print(f"\n{'='*60}")
    print(f"Carry Statistics: {num_rounds} rounds, {num_samples} samples")
    print(f"{'='*60}")

    carry_out_counts = []
    total_carry_bits = []

    for _ in range(num_samples):
        msg = [rng.randint(0, MASK32) for _ in range(16)]
        trace = sha256_compress_traced(msg, num_rounds)

        carry_outs = get_carry_out_vector(trace)
        n_sched, n_round, n_final = count_additions(trace)

        carry_out_counts.append(sum(carry_outs))
        total_carry_bits.append(len(carry_outs))

    total = total_carry_bits[0]
    mean_outs = sum(carry_out_counts) / len(carry_out_counts)
    print(f"Total additions tracked: {total}")
    print(f"  Schedule: {n_sched}, Round: {n_round}, Final: {n_final}")
    print(f"Mean carry-out = 1: {mean_outs:.1f}/{total} ({100*mean_outs/total:.1f}%)")
    print(f"  (Expected for random: 50%)")
    return total


def experiment_system_structure(rounds_list=None, seed=42):
    """
    Analyze Q∩T system structure at different round counts.

    Measures:
      - Number of variables and equations
      - Rank of linear subsystem
      - Kernel dimension (= degrees of freedom)
      - Number of quadratic equations and terms
    """
    if rounds_list is None:
        rounds_list = [1, 2, 3, 4, 5, 6, 8, 10, 12, 16]

    rng = random.Random(seed)

    print(f"\n{'='*60}")
    print(f"Q∩T System Structure Analysis")
    print(f"{'='*60}")
    print(f"{'R':>3} | {'Vars':>6} | {'LinEq':>6} | {'QuadEq':>6} | "
          f"{'Rank':>6} | {'Kernel':>6} | {'QuadTerms':>9}")
    print(f"{'-'*3}-+-{'-'*6}-+-{'-'*6}-+-{'-'*6}-+-{'-'*6}-+-{'-'*6}-+-{'-'*9}")

    results = []
    for R in rounds_list:
        msg = [rng.randint(0, MASK32) for _ in range(16)]
        system, trace = build_qt_system(msg, R)
        st = system.stats()

        # Compute rank of linear part
        rows = system.to_linear_rows()
        rank = gf2_rank(rows, st['num_vars'])
        kernel_dim = st['num_vars'] - rank

        quad_terms = st['quad_terms_total']

        print(f"{R:>3} | {st['num_vars']:>6} | {st['num_linear_eq']:>6} | "
              f"{st['num_quadratic_eq']:>6} | {rank:>6} | {kernel_dim:>6} | "
              f"{quad_terms:>9}")

        results.append({
            'rounds': R,
            'vars': st['num_vars'],
            'linear_eq': st['num_linear_eq'],
            'quadratic_eq': st['num_quadratic_eq'],
            'rank': rank,
            'kernel_dim': kernel_dim,
            'quad_terms': quad_terms,
        })

    return results


def experiment_q_variety_sampling(num_rounds=4, num_samples=5, num_variety=50,
                                  seed=42):
    """
    Sample points from the Q-variety and measure T-consistency.

    For each random message M:
    1. Build Q system with M's carries
    2. Solve linear part → get Q-variety
    3. Sample from Q-variety
    4. For each sample, compute actual carries and measure T-mismatch

    This measures the "distance" between Q-variety and T-feasible set.
    """
    rng = random.Random(seed)

    print(f"\n{'='*60}")
    print(f"Q-Variety Sampling: {num_rounds} rounds")
    print(f"{'='*60}")

    all_mismatches = []

    for s in range(num_samples):
        msg = [rng.randint(0, MASK32) for _ in range(16)]
        system, trace = build_qt_system(msg, num_rounds)
        target = trace.hash_words
        carries = get_all_carry_chains(trace)

        rows = system.to_linear_rows()
        vm = system.vm
        lin_result = gf2_solve(rows, vm.num_vars)

        if lin_result is None:
            print(f"  Sample {s}: linear system INCONSISTENT (shouldn't happen)")
            continue

        particular, kernel = lin_result

        # The original message should satisfy the system (mismatch = 0)
        # Sample from variety
        mismatches = []
        hash_matches = 0
        for v in range(num_variety):
            sample_bits = particular
            for kvec in kernel:
                if rng.random() < 0.5:
                    sample_bits ^= kvec

            # Extract message
            sample_msg = []
            for w in range(16):
                word = 0
                for b in range(32):
                    if (sample_bits >> vm.w(w, b)) & 1:
                        word |= (1 << b)
                sample_msg.append(word)

            # Check T-consistency
            sample_trace = sha256_compress_traced(sample_msg, num_rounds)
            sample_carries = get_all_carry_chains(sample_trace)

            mismatch = sum(
                bin(carries[i] ^ sample_carries[i]).count('1')
                for i in range(len(carries))
            )
            mismatches.append(mismatch)

            # Check hash match
            if sample_trace.hash_words == target:
                hash_matches += 1

        all_mismatches.extend(mismatches)

        if mismatches:
            print(f"  Sample {s}: kernel_dim={len(kernel)}, "
                  f"mismatch min={min(mismatches)} max={max(mismatches)} "
                  f"mean={sum(mismatches)/len(mismatches):.1f}, "
                  f"hash_matches={hash_matches}/{num_variety}")

    if all_mismatches:
        zeros = sum(1 for m in all_mismatches if m == 0)
        print(f"\nOverall: {zeros}/{len(all_mismatches)} T-consistent "
              f"({100*zeros/len(all_mismatches):.2f}%)")
        print(f"  Mean mismatch: {sum(all_mismatches)/len(all_mismatches):.1f} bits")
        if zeros > 0:
            print(f"  *** T-consistent solutions found in Q-variety! ***")

    return all_mismatches


def experiment_quadratic_structure(num_rounds=4, seed=42):
    """
    Analyze the quadratic equations in detail.

    For each quadratic equation, report:
    - Which variables are involved
    - Number of quadratic terms
    - Sparsity pattern
    """
    rng = random.Random(seed)
    msg = [rng.randint(0, MASK32) for _ in range(16)]
    system, trace = build_qt_system(msg, num_rounds)
    vm = system.vm

    print(f"\n{'='*60}")
    print(f"Quadratic Equation Structure: {num_rounds} rounds")
    print(f"{'='*60}")
    print(f"Total quadratic equations: {len(system.q_quadratic)}")

    if not system.q_quadratic:
        print("  (No quadratic equations — system is purely linear)")
        return

    # Analyze quadratic terms
    term_counts = [len(eq.quadratic) for eq in system.q_quadratic]
    var_counts = [len(eq.linear) + 2*len(eq.quadratic) for eq in system.q_quadratic]

    print(f"Quadratic terms per equation: min={min(term_counts)}, "
          f"max={max(term_counts)}, mean={sum(term_counts)/len(term_counts):.1f}")
    print(f"Variables per equation: min={min(var_counts)}, "
          f"max={max(var_counts)}, mean={sum(var_counts)/len(var_counts):.1f}")

    # Show first few equations
    print(f"\nFirst 5 quadratic equations:")
    for i, eq in enumerate(system.q_quadratic[:5]):
        lin_vars = sorted(eq.linear)
        quad_vars = sorted(eq.quadratic)
        lin_desc = [vm.describe(v) for v in lin_vars[:3]]
        quad_desc = [f"{vm.describe(v1)}*{vm.describe(v2)}" for v1, v2 in quad_vars[:3]]
        print(f"  [{i}] const={eq.constant}, "
              f"linear({len(lin_vars)})={lin_desc}{'...' if len(lin_vars)>3 else ''}, "
              f"quad({len(quad_vars)})={quad_desc}{'...' if len(quad_vars)>3 else ''}")


def experiment_cgaw_convergence(rounds_list=None, max_iter=50, seed=42):
    """
    Test CGAW convergence at different round counts.
    """
    if rounds_list is None:
        rounds_list = [1, 2, 3, 4]

    from qt_solver.solvers import solve_cgaw

    print(f"\n{'='*60}")
    print(f"CGAW Convergence Experiment")
    print(f"{'='*60}")

    for R in rounds_list:
        print(f"\n--- {R} rounds ---")
        result = solve_cgaw(R, max_iter=max_iter, seed=seed, verbose=True)
        if result.carry_mismatches:
            valid = [x for x in result.carry_mismatches if x >= 0]
            if valid:
                print(f"  Mismatch trajectory (first 10): "
                      f"{valid[:10]}")


def experiment_verify_q_system(num_rounds=4, seed=42):
    """
    Verify that the original message satisfies the Q system.

    The message used to generate the carries should be a solution
    of the Q system (since carries are self-consistent).
    """
    rng = random.Random(seed)

    print(f"\n{'='*60}")
    print(f"Q System Verification: {num_rounds} rounds")
    print(f"{'='*60}")

    num_tests = 10
    all_ok = True

    for t in range(num_tests):
        msg = [rng.randint(0, MASK32) for _ in range(16)]
        trace = sha256_compress_traced(msg, num_rounds)
        carries = get_all_carry_chains(trace)
        target = trace.hash_words

        vm = VarMap(num_rounds)
        system = QTSystem(num_rounds, carries, vm)
        system.build(target)

        # Build assignment from actual message and states
        assignment = {}

        # Message bits
        for w in range(16):
            for b in range(32):
                assignment[vm.w(w, b)] = (msg[w] >> b) & 1

        # Schedule bits
        for w in range(16, num_rounds):
            for b in range(32):
                assignment[vm.w(w, b)] = (trace.schedule[w] >> b) & 1

        # State bits
        for r in range(1, num_rounds + 1):
            state = trace.states[r]
            for reg in range(8):
                for b in range(32):
                    assignment[vm.s(r, reg, b)] = (state[reg] >> b) & 1

        # Verify all linear equations
        lin_violations = 0
        for eq in system.q_linear:
            if eq.evaluate(assignment) != 0:
                lin_violations += 1

        # Verify all quadratic equations
        quad_violations = 0
        for eq in system.q_quadratic:
            if eq.evaluate(assignment) != 0:
                quad_violations += 1

        status = "OK" if (lin_violations == 0 and quad_violations == 0) else "FAIL"
        if status == "FAIL":
            all_ok = False
        print(f"  Test {t}: {status} "
              f"(linear: {lin_violations}/{len(system.q_linear)}, "
              f"quad: {quad_violations}/{len(system.q_quadratic)})")

    print(f"\nVerification: {'ALL PASSED' if all_ok else 'FAILURES DETECTED'}")
    return all_ok
