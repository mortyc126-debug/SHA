#!/usr/bin/env python3
"""
Q∩T Algebra Solver — Main Entry Point

Global solver for the hybrid SHA-256 system:
  Q: 256 quadratic GF(2) equations (SHA-256 with fixed carry)
  T: 14,336 threshold equations (carry self-consistency)
  512 variables (message bits)

Based on methodology_v20.md Section 216 — the sole remaining
promising direction for SHA-256 cryptanalysis.

Usage:
  python run.py                  # Run all experiments
  python run.py verify           # Verify system correctness
  python run.py structure        # Analyze system structure
  python run.py cgaw             # Test CGAW solver
  python run.py sampling         # Q-variety sampling experiment
  python run.py full             # Full experiment suite
"""

import sys
import time


def run_verification():
    """Step 1: Verify correctness of all components."""
    from qt_solver.experiments import verify_sha256_traced, experiment_verify_q_system

    print("=" * 60)
    print("PHASE 1: VERIFICATION")
    print("=" * 60)

    # 1a. SHA-256 traced vs standard
    print("\n--- SHA-256 Traced vs Standard ---")
    for R in [1, 4, 16, 64]:
        errors, total = verify_sha256_traced(100, R)
        status = "PASS" if errors == 0 else f"FAIL ({errors} errors)"
        print(f"  {R:>2} rounds: {status} ({total} tests)")

    # 1b. Q system verification
    print("\n--- Q System Self-Consistency ---")
    for R in [1, 2, 3, 4, 8]:
        print(f"\n  {R} rounds:")
        ok = experiment_verify_q_system(R)
        if not ok:
            print(f"  *** VERIFICATION FAILED at {R} rounds ***")
            return False

    print("\n✓ All verifications passed")
    return True


def run_structure_analysis():
    """Step 2: Analyze Q∩T system structure."""
    from qt_solver.experiments import (
        experiment_carry_statistics,
        experiment_system_structure,
        experiment_quadratic_structure,
    )

    print("\n" + "=" * 60)
    print("PHASE 2: STRUCTURAL ANALYSIS")
    print("=" * 60)

    experiment_carry_statistics(num_rounds=64, num_samples=50)
    experiment_system_structure(rounds_list=[1, 2, 3, 4, 5, 6, 8, 10, 12, 16])
    experiment_quadratic_structure(num_rounds=4)


def run_sampling():
    """Step 3: Q-variety sampling experiment."""
    from qt_solver.experiments import experiment_q_variety_sampling

    print("\n" + "=" * 60)
    print("PHASE 3: Q-VARIETY SAMPLING")
    print("=" * 60)

    for R in [1, 2, 3, 4]:
        experiment_q_variety_sampling(num_rounds=R, num_samples=5, num_variety=100)


def run_cgaw():
    """Step 4: CGAW solver experiment."""
    from qt_solver.experiments import experiment_cgaw_convergence

    print("\n" + "=" * 60)
    print("PHASE 4: CGAW SOLVER")
    print("=" * 60)

    experiment_cgaw_convergence(rounds_list=[1, 2, 3, 4], max_iter=50)


def run_solvers():
    """Step 5: Solver comparison."""
    from qt_solver.solvers import solve_cgaw, solve_carry_partition, analyze_linear_relaxation

    print("\n" + "=" * 60)
    print("PHASE 5: SOLVER COMPARISON")
    print("=" * 60)

    for R in [1, 2, 3, 4]:
        print(f"\n{'─'*40}")
        print(f"  {R} rounds")
        print(f"{'─'*40}")

        # Linear relaxation analysis
        analyze_linear_relaxation(R, num_samples=5)

        # CGAW
        print(f"\n  CGAW:")
        r1 = solve_cgaw(R, max_iter=100, seed=42)
        print(f"    Result: {r1}")

        # Carry Partition Walk
        print(f"\n  Carry Partition Walk:")
        r2 = solve_carry_partition(R, max_iter=100, seed=42)
        print(f"    Result: {r2}")


def main():
    args = sys.argv[1:] if len(sys.argv) > 1 else ['full']
    cmd = args[0].lower()

    t0 = time.time()

    if cmd == 'verify':
        run_verification()
    elif cmd == 'structure':
        run_structure_analysis()
    elif cmd == 'sampling':
        run_sampling()
    elif cmd == 'cgaw':
        run_cgaw()
    elif cmd == 'solvers':
        run_solvers()
    elif cmd == 'full':
        print("╔══════════════════════════════════════════════════════════╗")
        print("║   Q∩T Algebra Solver for SHA-256 — Full Experiment     ║")
        print("║   Based on methodology_v20.md Section 216              ║")
        print("╚══════════════════════════════════════════════════════════╝")

        if not run_verification():
            print("\n*** Verification failed, aborting ***")
            sys.exit(1)
        run_structure_analysis()
        run_sampling()
        run_cgaw()

    elapsed = time.time() - t0
    print(f"\n{'='*60}")
    print(f"Total elapsed: {elapsed:.1f}s")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
