"""
Component-based solver for the decomposed Q∩T quadratic system.

After relinearization, the remaining quadratic system decomposes
into 29 independent connected components. This module:

1. Extracts and analyzes each component
2. Solves small components exhaustively (brute force on their vars)
3. Propagates solutions back to narrow the main linear system
4. Measures total kernel compression from component solving
"""

import random
import time
from collections import defaultdict
from itertools import product as cartesian

from qt_solver.sha256_traced import (
    MASK32, sha256_compress_traced, get_all_carry_chains,
)
from qt_solver.gf2 import gf2_gaussian_eliminate, gf2_solve, gf2_rank
from qt_solver.unified_v2 import build_unified_v2
from qt_solver.relinearize import (
    iterative_relinearize, PivotMap, iter_bits,
    relinearize_once, substitute_quadratic_eq,
)


class Component:
    """One independent component of the quadratic system."""

    def __init__(self, comp_id, eq_indices):
        self.id = comp_id
        self.eq_indices = eq_indices
        self.equations = []       # (lin_mask, quad_pairs, const)
        self.variables = set()    # all free vars in this component
        self.quad_vars = set()    # vars appearing in quadratic terms
        self.solutions = None     # list of valid assignments (if solved)
        self.num_solutions = None

    def __repr__(self):
        return (f"Component({self.id}: {len(self.equations)} eqs, "
                f"{len(self.variables)} vars, {len(self.quad_vars)} quad_vars)")


def extract_components(num_rounds, seed=42, verbose=True):
    """
    Extract independent components from the remaining quadratic system.

    Returns:
        components: list of Component objects
        context: dict with system info needed for propagation
    """
    rng = random.Random(seed)
    msg = [rng.randint(0, MASK32) for _ in range(16)]
    system, trace = build_unified_v2(msg, num_rounds)
    n = system.vm.num_vars

    # Relinearize
    relin = iterative_relinearize(
        system.linear_rows, system.quad_equations, n, verbose=False
    )
    remaining = relin['quad_equations']
    final_rank = relin['final_rank']

    if not remaining:
        if verbose:
            print(f"R={num_rounds}: No remaining quadratics")
        return [], {'system': system, 'trace': trace, 'n': n,
                    'final_rank': final_rank, 'relin': relin}

    # Build variable-equation adjacency
    eq_vars = []
    for i, (lin_row, qpairs, const) in enumerate(remaining):
        var_mask = (1 << n) - 1
        lin_vars = set(iter_bits(lin_row & var_mask))
        qv = set()
        for v1, v2 in qpairs:
            qv.add(v1)
            qv.add(v2)
        eq_vars.append((lin_vars | qv, qv))

    # Union-find for connected components
    parent = list(range(len(remaining)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    var_to_eqs = defaultdict(list)
    for i, (all_v, _) in enumerate(eq_vars):
        for v in all_v:
            var_to_eqs[v].append(i)

    for v, eqs in var_to_eqs.items():
        for j in range(1, len(eqs)):
            union(eqs[0], eqs[j])

    # Group equations by component
    comp_groups = defaultdict(list)
    for i in range(len(remaining)):
        comp_groups[find(i)].append(i)

    # Build Component objects
    components = []
    for comp_id, (root, indices) in enumerate(
            sorted(comp_groups.items(), key=lambda x: -len(x[1]))):
        comp = Component(comp_id, indices)
        for idx in indices:
            comp.equations.append(remaining[idx])
            all_v, qv = eq_vars[idx]
            comp.variables |= all_v
            comp.quad_vars |= qv
        components.append(comp)

    if verbose:
        print(f"\n{'='*60}")
        print(f"Component Analysis: R={num_rounds}")
        print(f"{'='*60}")
        print(f"Total remaining quad: {len(remaining)}")
        print(f"Components: {len(components)}")
        print(f"\n{'ID':>3} | {'Eqs':>4} | {'Vars':>5} | {'QVars':>5} | Status")
        print(f"{'-'*3}-+-{'-'*4}-+-{'-'*5}-+-{'-'*5}-+-------")
        for c in components:
            solvable = "solvable" if len(c.quad_vars) <= 20 else "large"
            print(f"{c.id:>3} | {len(c.equations):>4} | "
                  f"{len(c.variables):>5} | {len(c.quad_vars):>5} | {solvable}")

    context = {
        'system': system, 'trace': trace, 'n': n,
        'final_rank': final_rank, 'relin': relin,
        'msg': msg,
    }
    return components, context


def solve_component_exhaustive(comp, n, max_vars=24, verbose=False):
    """
    Solve a component by exhaustive enumeration of its quadratic variables.

    For each assignment of quad_vars:
      - Substitute into equations
      - Check if all equations are satisfied
      - If yes: record as valid assignment

    Returns number of solutions found.
    """
    qvars = sorted(comp.quad_vars)
    if len(qvars) > max_vars:
        if verbose:
            print(f"  Component {comp.id}: {len(qvars)} quad vars > {max_vars}, skipping")
        return None

    # For each equation, precompute the structure
    # Each eq: linear_mask XOR quad_pairs XOR const = 0
    # When we fix quad_vars, the quad terms become constants,
    # and the linear terms involving quad_vars also become constants.
    # The remaining linear terms (non-quad vars) must be satisfiable.

    # Actually, we need to check: for a given assignment of quad_vars,
    # do the equations REDUCE to consistent linear constraints on the
    # remaining (non-quad) variables?

    # Simpler approach for small components:
    # Enumerate all assignments of ALL variables in the component.
    all_vars = sorted(comp.variables)

    if len(all_vars) > max_vars:
        # Too many — just enumerate quad vars
        return _solve_by_quad_enum(comp, qvars, n, verbose)

    # Enumerate all assignments of all vars
    solutions = []
    for bits in range(1 << len(all_vars)):
        assignment = {}
        for i, v in enumerate(all_vars):
            assignment[v] = (bits >> i) & 1

        # Check all equations
        satisfied = True
        for lin_row, qpairs, const in comp.equations:
            var_mask = (1 << n) - 1
            val = const
            for v in iter_bits(lin_row & var_mask):
                val ^= assignment.get(v, 0)
            for v1, v2 in qpairs:
                val ^= (assignment.get(v1, 0) & assignment.get(v2, 0))
            if val != 0:
                satisfied = False
                break

        if satisfied:
            solutions.append(dict(assignment))

    comp.solutions = solutions
    comp.num_solutions = len(solutions)

    if verbose:
        print(f"  Component {comp.id}: {len(all_vars)} vars, "
              f"{len(comp.equations)} eqs → {len(solutions)} solutions")

    return len(solutions)


def _solve_by_quad_enum(comp, qvars, n, verbose):
    """Enumerate quad var assignments, check linear consistency."""
    solutions = []
    non_quad_vars = sorted(comp.variables - comp.quad_vars)

    for bits in range(1 << len(qvars)):
        q_assign = {}
        for i, v in enumerate(qvars):
            q_assign[v] = (bits >> i) & 1

        # Substitute quad assignment into equations
        # Each equation becomes linear in non_quad_vars
        linear_rows = []
        consistent = True

        for lin_row, qpairs, const in comp.equations:
            var_mask = (1 << n) - 1
            new_const = const

            # Evaluate quad terms with assignment
            for v1, v2 in qpairs:
                new_const ^= (q_assign.get(v1, 0) & q_assign.get(v2, 0))

            # Evaluate linear terms that are quad vars
            new_lin = 0
            for v in iter_bits(lin_row & var_mask):
                if v in q_assign:
                    new_const ^= q_assign[v]
                else:
                    new_lin ^= (1 << v)

            # This is now: new_lin XOR new_const = 0
            if new_lin == 0:
                if new_const != 0:
                    consistent = False
                    break
                # Trivially satisfied
            else:
                row = new_lin | (new_const << n)
                linear_rows.append(row)

        if not consistent:
            continue

        # Check if the linear system on non_quad_vars is consistent
        if not linear_rows:
            # All equations trivially satisfied
            solutions.append(dict(q_assign))
            continue

        result = gf2_solve(linear_rows, n)
        if result is not None:
            # Consistent! Record the quad assignment
            solutions.append(dict(q_assign))

    comp.solutions = solutions
    comp.num_solutions = len(solutions)

    if verbose:
        print(f"  Component {comp.id}: {len(qvars)} qvars, "
              f"{len(non_quad_vars)} non-qvars → "
              f"{len(solutions)}/{1 << len(qvars)} quad assignments valid")

    return len(solutions)


def solve_all_components(components, n, max_vars=20, verbose=True):
    """
    Solve all small components, report statistics.
    """
    if verbose:
        print(f"\n{'='*60}")
        print(f"Solving Components (max_vars={max_vars})")
        print(f"{'='*60}")

    t0 = time.time()
    solved = 0
    total_constraints_bits = 0

    for comp in components:
        if len(comp.quad_vars) <= max_vars:
            nsol = solve_component_exhaustive(comp, n, max_vars, verbose)
            if nsol is not None:
                solved += 1
                # Information gained: log2(solutions) vs log2(search space)
                search_space = 1 << len(comp.quad_vars)
                if nsol > 0:
                    import math
                    bits_gained = len(comp.quad_vars) - math.log2(nsol) if nsol > 0 else len(comp.quad_vars)
                    total_constraints_bits += bits_gained
                else:
                    if verbose:
                        print(f"    Component {comp.id}: NO SOLUTIONS (contradiction!)")

    elapsed = time.time() - t0

    if verbose:
        print(f"\nSummary:")
        print(f"  Solved: {solved}/{len(components)} components")
        print(f"  Total constraint bits gained: {total_constraints_bits:.1f}")
        print(f"  Elapsed: {elapsed:.2f}s")

        # Detail per solved component
        for comp in components:
            if comp.num_solutions is not None:
                import math
                qv = len(comp.quad_vars)
                ns = comp.num_solutions
                bits = qv - (math.log2(ns) if ns > 0 else 0)
                print(f"    C{comp.id}: {qv} qvars, {ns} solutions, "
                      f"{bits:.1f} bits constrained, "
                      f"P(valid)={ns/(1<<qv):.4f}")

    return solved, total_constraints_bits


def propagate_solutions(components, context, verbose=True):
    """
    Use solved component solutions to add linear constraints
    to the main system, compressing the kernel further.

    For each solved component with solutions:
      - If exactly 1 solution: all vars are fixed → add as linear eqs
      - If k solutions: the solution set forms an affine subspace
        → extract linear constraints on variables
    """
    n = context['n']
    system = context['system']

    if verbose:
        print(f"\n{'='*60}")
        print(f"Propagating Component Solutions")
        print(f"{'='*60}")

    new_linear_rows = []

    for comp in components:
        if comp.solutions is None or not comp.solutions:
            continue

        sols = comp.solutions
        if len(sols) == 1:
            # Unique solution: fix all variables
            for v, val in sols[0].items():
                row = (1 << v) | (val << n)
                new_linear_rows.append(row)
            if verbose:
                print(f"  C{comp.id}: unique solution → {len(sols[0])} vars fixed")

        elif len(sols) >= 2:
            # Multiple solutions: find linear constraints
            # A linear constraint L(x)=c holds iff it's the same for ALL solutions
            vars_list = sorted(comp.variables)

            # Check each variable: is it constant across all solutions?
            fixed_vars = {}
            for v in vars_list:
                vals = set(s.get(v, 0) for s in sols)
                if len(vals) == 1:
                    fixed_vars[v] = vals.pop()

            for v, val in fixed_vars.items():
                row = (1 << v) | (val << n)
                new_linear_rows.append(row)

            # Check pairwise XOR constraints: x_i XOR x_j = c
            pair_constraints = 0
            for i, v1 in enumerate(vars_list):
                for v2 in vars_list[i+1:]:
                    xor_vals = set(
                        (s.get(v1, 0) ^ s.get(v2, 0)) for s in sols
                    )
                    if len(xor_vals) == 1:
                        c = xor_vals.pop()
                        row = (1 << v1) ^ (1 << v2)
                        if c:
                            row |= (1 << n)
                        new_linear_rows.append(row)
                        pair_constraints += 1

            if verbose:
                print(f"  C{comp.id}: {len(sols)} solutions → "
                      f"{len(fixed_vars)} fixed vars, "
                      f"{pair_constraints} pair constraints")

    if not new_linear_rows:
        if verbose:
            print("  No new constraints to propagate")
        return 0

    # Add to main system and measure compression
    all_rows = list(system.linear_rows) + new_linear_rows

    # Include linearized quad equations
    from qt_solver.advanced import _extract_new_linear
    all_rows += list(_extract_new_linear(system, n))

    old_rank = context['final_rank']
    new_rank = gf2_rank(all_rows, n)
    compression = new_rank - old_rank

    if verbose:
        print(f"\n  New constraints added: {len(new_linear_rows)}")
        print(f"  Old rank: {old_rank}, New rank: {new_rank}")
        print(f"  Kernel compression: {compression} bits")
        print(f"  Old kernel: {n - old_rank}, New kernel: {n - new_rank}")

    return compression


def full_component_pipeline(num_rounds=3, seed=42, verbose=True):
    """
    Full pipeline: extract → solve → propagate → measure.
    """
    if verbose:
        print(f"\n{'='*60}")
        print(f"FULL COMPONENT PIPELINE: R={num_rounds}")
        print(f"{'='*60}")

    t0 = time.time()

    # Extract
    components, context = extract_components(num_rounds, seed, verbose)

    if not components:
        return {'compression': 0}

    # Solve
    solved, bits_gained = solve_all_components(
        components, context['n'], max_vars=20, verbose=verbose
    )

    # Propagate
    compression = propagate_solutions(components, context, verbose)

    elapsed = time.time() - t0

    if verbose:
        print(f"\n{'─'*60}")
        print(f"Pipeline complete: {elapsed:.2f}s")
        print(f"Components solved: {solved}/{len(components)}")
        print(f"Bits gained from components: {bits_gained:.1f}")
        print(f"Kernel compression: {compression}")
        print(f"{'─'*60}")

    return {
        'num_components': len(components),
        'solved': solved,
        'bits_gained': bits_gained,
        'compression': compression,
        'elapsed': elapsed,
    }
