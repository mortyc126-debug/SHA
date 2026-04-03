"""
Solver strategies for the Q∩T system.

Key constraint (methodology_v20.md):
  - NOT layer-by-layer (reduces to forward/reverse computation)
  - NOT metric-based (thermostat kills all metrics by round 30)
  - NOT single-algebra (F39: can't linearize + and XOR simultaneously)
  - MUST handle Q and T simultaneously

Implemented strategies:
  1. CGAW (Carry-Guided Algebraic Walk) - main novel algorithm
  2. Partial carry enumeration - for small systems
  3. Linear relaxation analysis - structural insight
"""

import random
import time

from qt_solver.sha256_traced import (
    sha256_compress_traced, sha256_compress,
    get_all_carry_chains, MASK32,
)
from qt_solver.gf2 import (
    gf2_gaussian_eliminate, gf2_solve, gf2_rank, popcount,
)
from qt_solver.qt_system import QTSystem, VarMap, build_qt_system


# ─── Solver Results ───

class SolverResult:
    """Result from a solver run."""
    def __init__(self):
        self.success = False
        self.message = None       # Solution message words (if found)
        self.iterations = 0
        self.carry_mismatches = [] # History of T-violation counts
        self.elapsed_sec = 0.0
        self.extra = {}           # Solver-specific data

    def __repr__(self):
        status = "SUCCESS" if self.success else "FAILED"
        return (f"SolverResult({status}, iter={self.iterations}, "
                f"elapsed={self.elapsed_sec:.3f}s)")


# ─── Strategy 1: CGAW (Carry-Guided Algebraic Walk) ───

def solve_cgaw(num_rounds, target_hash=None, max_iter=100, seed=None, verbose=True):
    """
    Carry-Guided Algebraic Walk.

    Algorithm:
    1. Start with random message M₀
    2. Compute SHA(M₀) → get carry vector C₀ (self-consistent by construction)
    3. If target_hash is None, set target = SHA(M₀) (preimage of known hash)
    4. Build Q system with carry C₀ fixed
    5. Solve Q's linear part via Gaussian elimination
    6. Sample a solution M* from the Q-variety (random free variables)
    7. Compute actual carry C(M*)
    8. If C(M*) == C₀: SUCCESS (found valid preimage)
    9. Else: C₀ ← C(M*), goto step 4

    The key insight: Q-solving provides structure-aware jumps in message space,
    potentially much better than random bit flips. Each iteration re-linearizes
    around the current carry assignment, creating a global (not layer-by-layer)
    exploration of the Q∩T intersection.
    """
    rng = random.Random(seed)
    result = SolverResult()
    t0 = time.time()

    # Step 1: Random initial message
    msg = [rng.randint(0, MASK32) for _ in range(16)]

    # Step 2: Compute trace
    trace = sha256_compress_traced(msg, num_rounds)
    current_carries = get_all_carry_chains(trace)

    # Step 3: Set target
    if target_hash is None:
        target_hash = trace.hash_words
        if verbose:
            print(f"[CGAW] Target = SHA-{num_rounds}(random M)")

    if verbose:
        print(f"[CGAW] System: {num_rounds} rounds, target_hash[0]=0x{target_hash[0]:08x}")
        print(f"[CGAW] Starting CGAW with max_iter={max_iter}")

    for iteration in range(max_iter):
        # Step 4: Build Q system with current carries
        vm = VarMap(num_rounds)
        system = QTSystem(num_rounds, current_carries, vm)
        system.build(target_hash)

        # Step 5: Solve linear part
        linear_rows = system.to_linear_rows()
        lin_result = gf2_solve(linear_rows, vm.num_vars)

        if lin_result is None:
            # Linear system inconsistent with this carry assignment
            # Perturb and retry
            msg = [rng.randint(0, MASK32) for _ in range(16)]
            trace = sha256_compress_traced(msg, num_rounds)
            current_carries = get_all_carry_chains(trace)
            result.carry_mismatches.append(-1)  # -1 = inconsistent
            if verbose and iteration % 10 == 0:
                print(f"  iter {iteration}: linear system inconsistent, restarting")
            continue

        particular, kernel = lin_result

        # Step 6: Sample from Q-variety
        # particular + random combination of kernel vectors
        candidate_bits = particular
        for kvec in kernel:
            if rng.random() < 0.5:
                candidate_bits ^= kvec

        # Extract message words from candidate
        candidate_msg = []
        for w in range(16):
            word = 0
            for b in range(32):
                var_idx = vm.w(w, b)
                if (candidate_bits >> var_idx) & 1:
                    word |= (1 << b)
            candidate_msg.append(word)

        # Step 7: Compute actual carries for candidate
        candidate_trace = sha256_compress_traced(candidate_msg, num_rounds)
        candidate_carries = get_all_carry_chains(candidate_trace)

        # Step 8: Check carry consistency
        mismatch_bits = 0
        for i in range(len(current_carries)):
            diff = current_carries[i] ^ candidate_carries[i]
            mismatch_bits += bin(diff).count('1')

        result.carry_mismatches.append(mismatch_bits)

        if mismatch_bits == 0:
            # Check actual hash
            actual_hash = candidate_trace.hash_words
            if actual_hash == target_hash:
                result.success = True
                result.message = candidate_msg
                result.iterations = iteration + 1
                result.elapsed_sec = time.time() - t0
                if verbose:
                    print(f"  iter {iteration}: SUCCESS! Carry-consistent solution found.")
                return result

        # Step 9: Update carries
        current_carries = candidate_carries

        if verbose and iteration % 10 == 0:
            print(f"  iter {iteration}: carry_mismatch_bits={mismatch_bits}, "
                  f"kernel_dim={len(kernel)}")

    result.iterations = max_iter
    result.elapsed_sec = time.time() - t0
    result.extra['final_kernel_dim'] = len(kernel) if lin_result else 0
    if verbose:
        print(f"[CGAW] Finished {max_iter} iterations without convergence")
        if result.carry_mismatches:
            valid = [x for x in result.carry_mismatches if x >= 0]
            if valid:
                print(f"  min mismatch: {min(valid)}, "
                      f"max: {max(valid)}, mean: {sum(valid)/len(valid):.1f}")
    return result


# ─── Strategy 2: Brute-force Q∩T for tiny systems ───

def solve_bruteforce(num_rounds, target_hash, max_attempts=2**20, seed=None,
                     verbose=True):
    """
    Brute-force search: random messages, check hash match.
    Baseline for comparison.
    """
    rng = random.Random(seed)
    result = SolverResult()
    t0 = time.time()

    for i in range(max_attempts):
        msg = [rng.randint(0, MASK32) for _ in range(16)]
        h = sha256_compress(msg, num_rounds)
        if h == target_hash:
            result.success = True
            result.message = msg
            result.iterations = i + 1
            result.elapsed_sec = time.time() - t0
            if verbose:
                print(f"[BF] Found preimage at attempt {i+1}")
            return result

    result.iterations = max_attempts
    result.elapsed_sec = time.time() - t0
    if verbose:
        print(f"[BF] No solution in {max_attempts} attempts ({result.elapsed_sec:.1f}s)")
    return result


# ─── Strategy 3: Linear Relaxation Analysis ───

def analyze_linear_relaxation(num_rounds, num_samples=10, seed=None, verbose=True):
    """
    Analyze the linear part of the Q system.

    For multiple random messages:
    1. Build Q system with their carries
    2. Compute rank of linear part
    3. Measure kernel dimension
    4. Check how many Q-variety samples satisfy T

    Returns statistics about the Q∩T structure.
    """
    rng = random.Random(seed)
    stats = {
        'ranks': [],
        'kernel_dims': [],
        'q_variety_sizes': [],
        'carry_mismatches': [],
        'quad_eq_counts': [],
    }

    for sample in range(num_samples):
        msg = [rng.randint(0, MASK32) for _ in range(16)]
        trace = sha256_compress_traced(msg, num_rounds)
        carries = get_all_carry_chains(trace)
        target = trace.hash_words

        vm = VarMap(num_rounds)
        system = QTSystem(num_rounds, carries, vm)
        system.build(target)

        # Analyze linear part
        rows = system.to_linear_rows()
        rank = gf2_rank(rows, vm.num_vars)
        kernel_dim = vm.num_vars - rank

        stats['ranks'].append(rank)
        stats['kernel_dims'].append(kernel_dim)
        stats['quad_eq_counts'].append(len(system.q_quadratic))

        # Sample from Q-variety and check T
        lin_result = gf2_solve(rows, vm.num_vars)
        if lin_result is not None:
            particular, kernel = lin_result

            # Test a few random samples from the variety
            mismatches = []
            for _ in range(min(20, max(1, 2**min(len(kernel), 10)))):
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

                # Check carry consistency
                sample_trace = sha256_compress_traced(sample_msg, num_rounds)
                sample_carries = get_all_carry_chains(sample_trace)

                mismatch = 0
                for i in range(len(carries)):
                    mismatch += bin(carries[i] ^ sample_carries[i]).count('1')
                mismatches.append(mismatch)

            stats['carry_mismatches'].append(mismatches)

    if verbose:
        print(f"\n=== Linear Relaxation Analysis: {num_rounds} rounds ===")
        print(f"Variables: {vm.num_vars}")
        print(f"  Message:  512")
        print(f"  Schedule: {vm.num_sched_words * 32}")
        print(f"  State:    {num_rounds * 256}")
        print(f"Linear equations: {len(system.q_linear)}")
        print(f"Quadratic equations: {len(system.q_quadratic)}")
        print(f"Linear rank: {min(stats['ranks'])} - {max(stats['ranks'])} "
              f"(mean {sum(stats['ranks'])/len(stats['ranks']):.1f})")
        print(f"Kernel dim: {min(stats['kernel_dims'])} - {max(stats['kernel_dims'])} "
              f"(mean {sum(stats['kernel_dims'])/len(stats['kernel_dims']):.1f})")

        all_mm = []
        for mm_list in stats['carry_mismatches']:
            all_mm.extend(mm_list)
        if all_mm:
            print(f"Carry mismatches from Q-variety samples:")
            print(f"  min={min(all_mm)}, max={max(all_mm)}, "
                  f"mean={sum(all_mm)/len(all_mm):.1f}")
            zeros = sum(1 for m in all_mm if m == 0)
            print(f"  T-consistent samples: {zeros}/{len(all_mm)} "
                  f"({100*zeros/len(all_mm):.1f}%)")

    return stats


# ─── Strategy 4: Carry Partition Walk ───

def solve_carry_partition(num_rounds, target_hash=None, max_iter=200,
                          partition_size=32, seed=None, verbose=True):
    """
    Carry Partition Walk: fix most carries, enumerate a small subset.

    Instead of fixing ALL carries (which creates self-consistency issues),
    partition carries into 'anchored' (fixed from current best) and
    'free' (enumerated). For each free-carry assignment, solve the
    linear system and check global T consistency.

    This exploits the block-diagonal structure of SHA-256 carries
    (documented as 5 carry windows with 42-111x isolation between blocks).
    """
    rng = random.Random(seed)
    result = SolverResult()
    t0 = time.time()

    # Start with random message
    msg = [rng.randint(0, MASK32) for _ in range(16)]
    trace = sha256_compress_traced(msg, num_rounds)
    best_carries = get_all_carry_chains(trace)

    if target_hash is None:
        target_hash = trace.hash_words

    num_carry_chains = len(best_carries)
    best_mismatch = 0  # Starting carries are self-consistent

    if verbose:
        print(f"[CPW] Carry Partition Walk: {num_rounds} rounds, "
              f"{num_carry_chains} carry chains")

    for iteration in range(max_iter):
        # Select a random subset of carry chains to perturb
        chain_indices = list(range(num_carry_chains))
        rng.shuffle(chain_indices)
        perturb_indices = chain_indices[:partition_size]

        # Create perturbed carry vector
        trial_carries = list(best_carries)
        for idx in perturb_indices:
            # Flip random bits in the carry chain
            num_flips = rng.randint(1, 3)
            for _ in range(num_flips):
                bit = rng.randint(0, 31)
                trial_carries[idx] ^= (1 << bit)

        # Build and solve Q system with trial carries
        vm = VarMap(num_rounds)
        system = QTSystem(num_rounds, trial_carries, vm)
        system.build(target_hash)

        rows = system.to_linear_rows()
        lin_result = gf2_solve(rows, vm.num_vars)

        if lin_result is None:
            continue  # Inconsistent, try another perturbation

        particular, kernel = lin_result

        # Sample from variety
        candidate_bits = particular
        for kvec in kernel:
            if rng.random() < 0.5:
                candidate_bits ^= kvec

        # Extract message
        candidate_msg = []
        for w in range(16):
            word = 0
            for b in range(32):
                if (candidate_bits >> vm.w(w, b)) & 1:
                    word |= (1 << b)
            candidate_msg.append(word)

        # Check T-consistency
        candidate_trace = sha256_compress_traced(candidate_msg, num_rounds)
        actual_carries = get_all_carry_chains(candidate_trace)

        mismatch = sum(bin(trial_carries[i] ^ actual_carries[i]).count('1')
                       for i in range(num_carry_chains))

        result.carry_mismatches.append(mismatch)

        # Check for solution
        if candidate_trace.hash_words == target_hash:
            result.success = True
            result.message = candidate_msg
            result.iterations = iteration + 1
            result.elapsed_sec = time.time() - t0
            if verbose:
                print(f"  iter {iteration}: SUCCESS!")
            return result

        # Update best if actual carries are self-consistent
        # (they always are for actual SHA computation)
        best_carries = actual_carries

        if verbose and iteration % 20 == 0:
            print(f"  iter {iteration}: mismatch={mismatch}, kernel_dim={len(kernel)}")

    result.iterations = max_iter
    result.elapsed_sec = time.time() - t0
    if verbose:
        print(f"[CPW] Finished {max_iter} iterations")
    return result
