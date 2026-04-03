"""
Three parallel improvement directions for the Q∩T solver.

Direction 1: Degree-2 Elevation
  Multiply remaining quadratic equations by their own variables
  to generate degree-3 → linearize → more constraints.

Direction 2: Structural Decomposition
  Analyze sparsity of remaining quadratics.
  Check if system decomposes into independent sub-problems.

Direction 3: Exploit Fully-Linear R=2
  R=2 after relinearization has 0 remaining quadratics!
  946 free variables, any assignment solves Q.
  Sample from Q-variety and check T-consistency.
"""

import random
import time
from collections import defaultdict

from qt_solver.sha256_traced import (
    MASK32, sha256_compress_traced, get_all_carry_chains,
)
from qt_solver.gf2 import gf2_gaussian_eliminate, gf2_solve, gf2_rank
from qt_solver.unified_v2 import build_unified_v2
from qt_solver.relinearize import (
    iterative_relinearize, relinearize_once, PivotMap, iter_bits,
)


# ════════════════════════════════════════════════════════════
# Direction 1: Degree-2 Elevation
# ════════════════════════════════════════════════════════════

def degree2_elevation(num_rounds, seed=42, verbose=True):
    """
    Multiply remaining quadratic equations by variables
    to generate new constraints.

    For each quad eq E with products involving vars {v1..vk}:
      For each vi in {v1..vk}:
        Compute E * vi → degree-3 equation
        Linearize: x_i*x_j*x_k → fresh aux variable z_ijk
        This produces a new linear equation in extended space.

    The key: if two quad equations share variables, their
    degree-2 extensions may create redundant constraints that
    reveal hidden linear relationships.
    """
    rng = random.Random(seed)
    msg = [rng.randint(0, MASK32) for _ in range(16)]
    system, trace = build_unified_v2(msg, num_rounds)
    n = system.vm.num_vars

    if verbose:
        print(f"\n{'='*60}")
        print(f"Direction 1: Degree-2 Elevation (R={num_rounds})")
        print(f"{'='*60}")

    # Step 1: Relinearize to get remaining quad in free-var space
    relin = iterative_relinearize(
        system.linear_rows, system.quad_equations, n, verbose=False
    )

    remaining = relin['quad_equations']
    init_kernel = relin['final_kernel_dim']

    if not remaining:
        if verbose:
            print(f"No remaining quadratics — system is fully linear!")
            print(f"Kernel dim: {init_kernel}")
        return {'compression': 0, 'final_kernel': init_kernel}

    if verbose:
        print(f"After relinearization: kernel={init_kernel}, "
              f"remaining_quad={len(remaining)}")

    # Step 2: Collect all variables appearing in quadratic terms
    quad_vars = set()
    for lin_row, qpairs, const in remaining:
        for v1, v2 in qpairs:
            quad_vars.add(v1)
            quad_vars.add(v2)

    if verbose:
        print(f"Variables in quadratic terms: {len(quad_vars)}")

    # Step 3: For each quad eq × each of its quad variables,
    # generate degree-3 equation (linearized)
    # Product map for degree-2 (quadratic) auxiliary vars
    product2_map = {}  # (i,j) → aux_id
    # Product map for degree-3 (cubic) auxiliary vars
    product3_map = {}  # (i,j,k) → aux_id
    next_aux = n

    def get_prod2(a, b):
        nonlocal next_aux
        key = (min(a,b), max(a,b))
        if key not in product2_map:
            product2_map[key] = next_aux
            next_aux += 1
        return product2_map[key]

    def get_prod3(a, b, c):
        nonlocal next_aux
        key = tuple(sorted([a, b, c]))
        # Handle x*x*y = x*y (GF(2))
        if key[0] == key[1]:
            return get_prod2(key[0], key[2]) if key[0] != key[2] else key[0]
        if key[1] == key[2]:
            return get_prod2(key[0], key[1]) if key[0] != key[1] else key[0]
        if key not in product3_map:
            product3_map[key] = next_aux
            next_aux += 1
        return product3_map[key]

    new_rows = []

    for lin_row, qpairs, const in remaining:
        # Get variables in this equation's quadratic terms
        eq_qvars = set()
        for v1, v2 in qpairs:
            eq_qvars.add(v1)
            eq_qvars.add(v2)

        for mult_var in eq_qvars:
            # Multiply entire equation by mult_var
            # E * mult_var = 0
            # E = const ⊕ Σ(lin_vars) ⊕ Σ(v1*v2)
            # E * m = const*m ⊕ Σ(lin_var * m) ⊕ Σ(v1*v2*m)

            new_row = 0

            # const * mult_var
            if const:
                new_row ^= (1 << mult_var)

            # linear_var * mult_var → quadratic aux
            var_mask = (1 << n) - 1
            lin_bits = lin_row & var_mask
            for v in iter_bits(lin_bits):
                if v == mult_var:
                    # v * v = v
                    new_row ^= (1 << v)
                else:
                    aux = get_prod2(v, mult_var)
                    new_row ^= (1 << aux)

            # v1 * v2 * mult_var → cubic aux (or simpler)
            for v1, v2 in qpairs:
                if mult_var == v1:
                    # v1 * v2 * v1 = v1 * v2 (GF2: x²=x)
                    aux = get_prod2(v1, v2)
                    new_row ^= (1 << aux)
                elif mult_var == v2:
                    aux = get_prod2(v1, v2)
                    new_row ^= (1 << aux)
                else:
                    aux = get_prod3(v1, v2, mult_var)
                    if isinstance(aux, int) and aux < n:
                        # Collapsed to single var
                        new_row ^= (1 << aux)
                    else:
                        new_row ^= (1 << aux)

            new_rows.append(new_row)

    # Also add the original quad equations (linearized at degree 2)
    for lin_row, qpairs, const in remaining:
        row = lin_row & ((1 << n) - 1)
        for v1, v2 in qpairs:
            aux = get_prod2(v1, v2)
            row ^= (1 << aux)
        if const:
            row ^= (1 << next_aux)  # This is wrong, should be augmented col
        # Actually need augmented column at position = total_vars
        # We'll fix this after knowing total vars
        new_rows.append(row | (const << next_aux) if const else row)

    total_vars = next_aux

    if verbose:
        print(f"Degree-2 aux vars (products): {len(product2_map)}")
        print(f"Degree-3 aux vars (triples): {len(product3_map)}")
        print(f"Total extended vars: {total_vars}")
        print(f"New equations from elevation: {len(new_rows)}")

    # Fix augmented column position
    fixed_rows = []
    for row in new_rows:
        # Extract potential constant bits above variable space
        coeff = row & ((1 << total_vars) - 1)
        # Check if any high bit is set (as constant marker)
        high = row >> total_vars
        const_bit = high & 1
        fixed_row = coeff | (const_bit << total_vars)
        fixed_rows.append(fixed_row)

    # Gaussian eliminate
    rank = gf2_rank(fixed_rows, total_vars)
    final_kernel = total_vars - rank

    if verbose:
        print(f"\nAfter degree-2 Gaussian elimination:")
        print(f"  Rank: {rank}")
        print(f"  Kernel: {final_kernel}")
        print(f"  (vs initial extended kernel ≈ {init_kernel + len(product2_map) + len(product3_map)})")
        effective_compression = (init_kernel + len(product2_map) + len(product3_map)) - final_kernel
        print(f"  Effective compression from degree-2: {effective_compression} bits")

    return {
        'total_vars': total_vars,
        'rank': rank,
        'final_kernel': final_kernel,
        'num_prod2': len(product2_map),
        'num_prod3': len(product3_map),
        'new_equations': len(new_rows),
    }


# ════════════════════════════════════════════════════════════
# Direction 2: Structural Decomposition
# ════════════════════════════════════════════════════════════

def analyze_quad_structure(num_rounds, seed=42, verbose=True):
    """
    Analyze sparsity and structure of remaining quadratic equations.

    Checks:
    - Variable count per equation (sparsity)
    - Variable overlap between equations
    - Connected components in the variable-equation graph
    - Bit-position clustering
    """
    rng = random.Random(seed)
    msg = [rng.randint(0, MASK32) for _ in range(16)]
    system, trace = build_unified_v2(msg, num_rounds)
    n = system.vm.num_vars

    if verbose:
        print(f"\n{'='*60}")
        print(f"Direction 2: Structural Decomposition (R={num_rounds})")
        print(f"{'='*60}")

    relin = iterative_relinearize(
        system.linear_rows, system.quad_equations, n, verbose=False
    )
    remaining = relin['quad_equations']

    if not remaining:
        if verbose:
            print("No remaining quadratics!")
        return {}

    # Analyze each equation
    eq_stats = []
    all_vars = set()
    var_to_eqs = defaultdict(set)

    for i, (lin_row, qpairs, const) in enumerate(remaining):
        var_mask = (1 << n) - 1
        lin_vars = set(iter_bits(lin_row & var_mask))
        quad_vars = set()
        for v1, v2 in qpairs:
            quad_vars.add(v1)
            quad_vars.add(v2)

        all_eq_vars = lin_vars | quad_vars
        for v in all_eq_vars:
            var_to_eqs[v].add(i)
        all_vars |= all_eq_vars

        eq_stats.append({
            'lin_count': len(lin_vars),
            'quad_pairs': len(qpairs),
            'total_vars': len(all_eq_vars),
            'quad_only_vars': len(quad_vars - lin_vars),
        })

    if verbose:
        lin_counts = [s['lin_count'] for s in eq_stats]
        qp_counts = [s['quad_pairs'] for s in eq_stats]
        tv_counts = [s['total_vars'] for s in eq_stats]

        print(f"Remaining quadratic equations: {len(remaining)}")
        print(f"Total unique variables involved: {len(all_vars)}")
        print(f"\nPer-equation statistics:")
        print(f"  Linear terms: min={min(lin_counts)}, max={max(lin_counts)}, "
              f"mean={sum(lin_counts)/len(lin_counts):.1f}")
        print(f"  Quad pairs:   min={min(qp_counts)}, max={max(qp_counts)}, "
              f"mean={sum(qp_counts)/len(qp_counts):.1f}")
        print(f"  Total vars:   min={min(tv_counts)}, max={max(tv_counts)}, "
              f"mean={sum(tv_counts)/len(tv_counts):.1f}")

    # Connected components via union-find on equations
    parent = list(range(len(remaining)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        a, b = find(a), find(b)
        if a != b:
            parent[a] = b

    # Two equations are connected if they share a variable
    for v, eqs in var_to_eqs.items():
        eqs = list(eqs)
        for i in range(1, len(eqs)):
            union(eqs[0], eqs[i])

    components = defaultdict(list)
    for i in range(len(remaining)):
        components[find(i)].append(i)

    comp_sizes = sorted([len(c) for c in components.values()], reverse=True)

    if verbose:
        print(f"\nConnected components: {len(components)}")
        print(f"  Sizes: {comp_sizes[:10]}{'...' if len(comp_sizes) > 10 else ''}")
        if len(components) > 1:
            print(f"  *** System DECOMPOSES into {len(components)} independent parts! ***")
        else:
            print(f"  System is fully connected (no decomposition)")

    # Variable sharing density
    shared_counts = [len(eqs) for eqs in var_to_eqs.values()]
    if verbose and shared_counts:
        print(f"\nVariable sharing:")
        print(f"  Vars appearing in 1 eq: {sum(1 for c in shared_counts if c == 1)}")
        print(f"  Vars appearing in 2 eq: {sum(1 for c in shared_counts if c == 2)}")
        print(f"  Vars appearing in 3+ eq: {sum(1 for c in shared_counts if c >= 3)}")
        print(f"  Max sharing: {max(shared_counts)}")

    return {
        'num_quad': len(remaining),
        'total_vars': len(all_vars),
        'components': len(components),
        'component_sizes': comp_sizes,
        'eq_stats': eq_stats,
    }


# ════════════════════════════════════════════════════════════
# Direction 3: Exploit Fully-Linear R=2
# ════════════════════════════════════════════════════════════

def solve_r2_linear(num_samples=1000, seed=42, verbose=True):
    """
    R=2 is fully linear after relinearization (0 remaining quad).
    946 free variables — any assignment solves Q.

    Strategy:
    1. Build unified system for random M, target=SHA(M)
    2. Relinearize → get pivot map
    3. Sample from Q-variety (random free vars)
    4. Extract message from sample
    5. Check T-consistency (carry self-consistency)

    If the original M is a Q-solution with mismatch=0,
    samples "near" M in the Q-variety might also have low mismatch.
    """
    rng = random.Random(seed)

    if verbose:
        print(f"\n{'='*60}")
        print(f"Direction 3: Fully-Linear R=2 Exploitation")
        print(f"{'='*60}")

    msg = [rng.randint(0, MASK32) for _ in range(16)]
    system, trace = build_unified_v2(msg, 2)
    n = system.vm.num_vars
    target = trace.hash_words
    true_carries = get_all_carry_chains(trace)

    # Relinearize
    # First get the linear rows after relinearization
    pivots_data, reduced, rank = gf2_gaussian_eliminate(
        list(system.linear_rows), n
    )
    pmap = PivotMap(pivots_data, reduced, n)

    # Add quadratic equations → they become linear after pivot substitution
    all_rows = list(system.linear_rows)
    # Process quad equations through relinearization
    relin = iterative_relinearize(
        system.linear_rows, system.quad_equations, n, verbose=False
    )

    # Solve the fully linear system
    final_pivots, final_reduced, final_rank = gf2_gaussian_eliminate(
        list(system.linear_rows) +
        [r for r in _extract_new_linear(system, n)],
        n
    )

    result = gf2_solve(
        list(system.linear_rows) +
        [r for r in _extract_new_linear(system, n)],
        n
    )

    if result is None:
        if verbose:
            print("ERROR: System inconsistent!")
        return {}

    particular, kernel = result
    if verbose:
        print(f"Kernel dimension: {len(kernel)}")
        print(f"Sampling {num_samples} points from Q-variety...")

    # Verify original message is a solution
    orig_assignment = _msg_to_bitvec(msg, system)
    # Check it's in the variety
    orig_in_variety = _check_linear_system(
        list(system.linear_rows) + list(_extract_new_linear(system, n)),
        orig_assignment, n
    )

    if verbose:
        print(f"Original message in Q-variety: {orig_in_variety}")

    # Sample and check T-consistency
    mismatches = []
    hash_matches = 0
    min_mismatch = float('inf')
    best_msg = None

    t0 = time.time()
    for i in range(num_samples):
        # Random point in Q-variety
        sample = particular
        for kvec in kernel:
            if rng.getrandbits(1):
                sample ^= kvec

        # Extract message words
        sample_msg = _bitvec_to_msg(sample, system.vm)

        # Check T-consistency
        sample_trace = sha256_compress_traced(sample_msg, 2)
        sample_carries = get_all_carry_chains(sample_trace)

        mismatch = sum(
            bin(true_carries[j] ^ sample_carries[j]).count('1')
            for j in range(len(true_carries))
        )
        mismatches.append(mismatch)

        if sample_trace.hash_words == target:
            hash_matches += 1

        if mismatch < min_mismatch:
            min_mismatch = mismatch
            best_msg = sample_msg

    elapsed = time.time() - t0

    if verbose:
        print(f"\nResults ({elapsed:.2f}s):")
        print(f"  Hash matches: {hash_matches}/{num_samples}")
        print(f"  Carry mismatches:")
        print(f"    min={min(mismatches)}, max={max(mismatches)}, "
              f"mean={sum(mismatches)/len(mismatches):.1f}")
        print(f"    P(mismatch=0): {sum(1 for m in mismatches if m==0)}/{num_samples}")

        # Distribution
        buckets = defaultdict(int)
        for m in mismatches:
            buckets[m // 10 * 10] += 1
        print(f"  Distribution (by 10s):")
        for k in sorted(buckets.keys())[:8]:
            bar = '#' * (buckets[k] * 40 // num_samples)
            print(f"    [{k:>3}-{k+9:>3}]: {buckets[k]:>4} {bar}")

    return {
        'kernel_dim': len(kernel),
        'num_samples': num_samples,
        'hash_matches': hash_matches,
        'min_mismatch': min(mismatches),
        'mean_mismatch': sum(mismatches) / len(mismatches),
        'mismatches': mismatches,
    }


def solve_r2_targeted(num_outer=50, num_inner=200, seed=42, verbose=True):
    """
    Targeted R=2 solver: use DIFFERENT base messages to get
    different carry assignments, then sample from each Q-variety.

    The idea: different M give different carry vectors C(M).
    For each C(M), the Q-variety is different.
    We search for an M where the Q-variety has T-consistent points.
    """
    rng = random.Random(seed)

    if verbose:
        print(f"\n{'='*60}")
        print(f"Direction 3b: Targeted R=2 (multi-base)")
        print(f"{'='*60}")

    best_global = float('inf')

    for outer in range(num_outer):
        msg = [rng.randint(0, MASK32) for _ in range(16)]
        system, trace = build_unified_v2(msg, 2)
        n = system.vm.num_vars
        target = trace.hash_words

        all_rows = list(system.linear_rows) + list(_extract_new_linear(system, n))
        result = gf2_solve(all_rows, n)
        if result is None:
            continue

        particular, kernel = result

        best_local = float('inf')
        for inner in range(num_inner):
            sample = particular
            for kvec in kernel:
                if rng.getrandbits(1):
                    sample ^= kvec

            sample_msg = _bitvec_to_msg(sample, system.vm)
            sample_trace = sha256_compress_traced(sample_msg, 2)
            sample_carries = get_all_carry_chains(sample_trace)
            true_carries = get_all_carry_chains(trace)

            mismatch = sum(
                bin(true_carries[j] ^ sample_carries[j]).count('1')
                for j in range(len(true_carries))
            )

            if mismatch < best_local:
                best_local = mismatch

            if sample_trace.hash_words == target:
                if verbose:
                    print(f"  *** PREIMAGE FOUND at outer={outer}, inner={inner}! ***")
                return {'success': True, 'msg': sample_msg, 'iterations': outer * num_inner + inner}

        if best_local < best_global:
            best_global = best_local

        if verbose and outer % 10 == 0:
            print(f"  outer={outer}: best_local={best_local}, best_global={best_global}")

    if verbose:
        print(f"  Best mismatch across all attempts: {best_global}")
    return {'success': False, 'best_mismatch': best_global}


# ─── Helpers ───

def _extract_new_linear(system, n):
    """Extract linearized quadratic equations via pivot substitution."""
    pivots, reduced, rank = gf2_gaussian_eliminate(list(system.linear_rows), n)
    pmap = PivotMap(pivots, reduced, n)

    new_rows = []
    for lin_row, qpairs, const in system.quad_equations:
        from qt_solver.relinearize import substitute_quadratic_eq
        res_lin, res_quad, res_const = substitute_quadratic_eq(
            pmap, lin_row, qpairs, n
        )
        if len(res_quad) == 0:
            row = res_lin
            if res_const:
                row |= (1 << n)
            new_rows.append(row)
    return new_rows


def _msg_to_bitvec(msg_words, system):
    """Convert message words to bit vector."""
    bv = 0
    for w in range(16):
        for b in range(32):
            if (msg_words[w] >> b) & 1:
                bv |= (1 << system.vm.w(w, b))
    return bv


def _bitvec_to_msg(bv, vm):
    """Extract message words from bit vector."""
    msg = []
    for w in range(16):
        word = 0
        for b in range(32):
            if (bv >> vm.w(w, b)) & 1:
                word |= (1 << b)
        msg.append(word)
    return msg


def _check_linear_system(rows, assignment, n):
    """Check if assignment satisfies linear system."""
    for row in rows:
        coeff = row & ((1 << n) - 1)
        const = (row >> n) & 1
        val = bin(coeff & assignment).count('1') & 1
        if val != const:
            return False
    return True
