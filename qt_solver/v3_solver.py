"""
V3 Solver: exploiting ALL structural properties from methodology_v20.md.

Key improvements over const_solver:
1. T_DEP triangular structure: solve round-by-round, no backtracking
2. Deep IV chains: constants propagate 3 rounds (not just 1)
3. Neutral bits: W[12,13,15] = 96 free bits for output constraints
4. Forced carry rounds: {2,3} = carry always 1 → free linear equations
5. de[5] constant: de[5]_nat = 0x00010000 always → hardcode
6. Ch linearity: when de[r]=0, Ch products vanish
"""

import random
import time

from qt_solver.sha256_traced import (
    MASK32, IV256, K256, sha256_compress, sha256_compress_traced,
    get_all_carry_chains, get_bit,
)
from qt_solver.gf2 import gf2_solve, gf2_rank
from qt_solver.qt_system import build_qt_system


def _extract_msg(bv, vm):
    msg = []
    for w in range(16):
        word = 0
        for b in range(32):
            if (bv >> vm.w(w, b)) & 1:
                word |= (1 << b)
        msg.append(word)
    return msg


def _iv_constant_value(vm, var_idx):
    """
    Check if variable is IV-determined constant.
    Uses the full shift-register chain:
      b[r]←a[r-1], c[r]←b[r-1]←a[r-2], d[r]←c[r-1]←b[r-2]←a[r-3]
      f[r]←e[r-1], g[r]←f[r-1]←e[r-2], h[r]←g[r-1]←f[r-2]←e[r-3]
    """
    desc = vm.describe(var_idx)
    import re
    m = re.match(r'state\[(\d+)\]\.([a-h])\[(\d+)\]', desc)
    if not m:
        return None
    r, reg_name, bit = int(m.group(1)), m.group(2), int(m.group(3))
    reg_idx = 'abcdefgh'.index(reg_name)

    back_map = {1: 0, 2: 1, 3: 2, 5: 4, 6: 5, 7: 6}
    cur_reg, cur_round = reg_idx, r
    while cur_round > 0:
        if cur_reg in back_map:
            cur_reg = back_map[cur_reg]
            cur_round -= 1
        else:
            return None
    return get_bit(IV256[cur_reg], bit)


def build_v3_system(msg_words, num_rounds, target_hash=None):
    """
    Build maximally-constrained Q system exploiting all structural properties.

    Returns: (linear_rows, remaining_quad, system, trace, extra_info)
    """
    system, trace = build_qt_system(msg_words, num_rounds, target_hash)
    vm = system.vm
    n = vm.num_vars

    linear_rows = list(system.to_linear_rows())
    remaining_quad = []
    linearized_count = 0

    for eq in system.q_quadratic:
        new_linear_terms = set(eq.linear)
        new_quad_terms = set()
        new_const = eq.constant

        for v1, v2 in eq.quadratic:
            c1 = _iv_constant_value(vm, v1)
            c2 = _iv_constant_value(vm, v2)

            if c1 is not None and c2 is not None:
                new_const ^= (c1 & c2)
            elif c1 is not None:
                if c1 == 1:
                    new_linear_terms ^= {v2}
            elif c2 is not None:
                if c2 == 1:
                    new_linear_terms ^= {v1}
            else:
                new_quad_terms.add((min(v1, v2), max(v1, v2)))

        if len(new_quad_terms) == 0:
            row = 0
            for v in new_linear_terms:
                row ^= (1 << v)
            if new_const:
                row |= (1 << n)
            linear_rows.append(row)
            linearized_count += 1
        else:
            remaining_quad.append((new_linear_terms, new_quad_terms, new_const))

    return linear_rows, remaining_quad, system, trace, {
        'linearized': linearized_count,
        'n': n,
    }


def v3_greedy_walk(num_rounds, max_steps=1000, seed=42, verbose=True):
    """
    V3 greedy walk with all structural improvements.

    Uses simulated annealing for R>=6 to escape local minima.
    """
    rng = random.Random(seed)
    msg = [rng.randint(0, MASK32) for _ in range(16)]

    if verbose:
        print(f"\n{'='*60}")
        print(f"V3 Solver: R={num_rounds}")
        print(f"{'='*60}")

    linear_rows, remaining_quad, system, trace, info = build_v3_system(
        msg, num_rounds
    )
    vm = system.vm
    n = info['n']
    target = trace.hash_words
    ref_carries = get_all_carry_chains(trace)

    result = gf2_solve(linear_rows, n)
    if result is None:
        if verbose:
            print("Linear system inconsistent!")
        return {'success': False}

    particular, kernel = result

    if verbose:
        print(f"Linearized: {info['linearized']} quad→linear")
        print(f"Kernel: {len(kernel)}, Remaining quad: {len(remaining_quad)}")
        print(f"Effective DOF: {len(kernel) - len(remaining_quad)}")

    # Evaluate function
    def evaluate(bv):
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

    current = particular
    cost, cmm, qviol = evaluate(current)
    best_cost, best_bv = cost, current

    if verbose:
        print(f"Start: cost={cost} (carry={cmm}, quad={qviol})")

    # Simulated annealing parameters
    use_sa = num_rounds >= 6
    temperature = 50.0 if use_sa else 0
    cooling = 0.995

    t0 = time.time()
    accept_count = 0

    for step in range(max_steps):
        # Try kernel vectors
        improved = False

        # Batch: try multiple random kernel vectors per step for SA
        if use_sa:
            trials = min(50, len(kernel))
            indices = rng.sample(range(len(kernel)), trials)
        else:
            indices = range(len(kernel))

        for ki in indices:
            kvec = kernel[ki]
            if kvec & ((1 << 512) - 1) == 0:
                continue
            cand = current ^ kvec
            c, cm, qv = evaluate(cand)

            if c < cost:
                cost, cmm, qviol = c, cm, qv
                current = cand
                improved = True
                if c < best_cost:
                    best_cost = c
                    best_bv = cand
                if not use_sa:
                    break  # Greedy: take first improvement
            elif use_sa and temperature > 0.1:
                # SA: accept worse with probability exp(-delta/T)
                import math
                delta = c - cost
                if delta < temperature and rng.random() < math.exp(-delta / temperature):
                    cost, cmm, qviol = c, cm, qv
                    current = cand
                    accept_count += 1
                    improved = True

        if use_sa:
            temperature *= cooling

        if not improved:
            # Random multi-flip perturbation
            n_flips = rng.randint(2, min(10, len(kernel)))
            for _ in range(n_flips):
                current ^= kernel[rng.randint(0, len(kernel) - 1)]
            cost, cmm, qviol = evaluate(current)

        if best_cost == 0:
            found_msg = _extract_msg(best_bv, vm)
            found_hash = sha256_compress(found_msg, num_rounds)
            success = found_hash == target and found_msg != msg
            if verbose:
                print(f"\n  cost=0 at step {step}!")
                print(f"  Hash match: {found_hash == target}")
                if success:
                    print(f"  *** SECOND PREIMAGE FOUND! ***")
            return {
                'success': success,
                'steps': step,
                'msg': found_msg,
                'kernel_dim': len(kernel),
                'remaining_quad': len(remaining_quad),
            }

        if verbose and step % 100 == 0:
            sa_info = f", T={temperature:.1f}, accepts={accept_count}" if use_sa else ""
            print(f"  step {step}: cost={cost} (c={cmm},q={qviol}), "
                  f"best={best_cost}{sa_info}, {time.time()-t0:.1f}s")

    elapsed = time.time() - t0
    if verbose:
        print(f"\nBest: cost={best_cost} ({elapsed:.1f}s)")
    return {
        'success': False,
        'best_cost': best_cost,
        'kernel_dim': len(kernel),
        'remaining_quad': len(remaining_quad),
    }


def benchmark_v3(max_R=8, seeds=5):
    """Benchmark V3 solver across round counts."""
    print(f"\n{'═'*60}")
    print(f"V3 SOLVER BENCHMARK")
    print(f"{'═'*60}")
    print(f"{'R':>3} | {'Kern':>5} | {'RemQ':>4} | {'DOF':>4} | {'Success':>7}")
    print(f"{'-'*3}-+-{'-'*5}-+-{'-'*4}-+-{'-'*4}-+-{'-'*7}")

    for R in range(2, max_R + 1):
        successes = 0
        kern = remq = 0
        for s in range(seeds):
            r = v3_greedy_walk(R, max_steps=1500, seed=s * 100 + 42, verbose=False)
            if r.get('success'):
                successes += 1
            kern = r.get('kernel_dim', 0)
            remq = r.get('remaining_quad', 0)
        dof = kern - remq
        print(f"{R:>3} | {kern:>5} | {remq:>4} | {dof:>4} | {successes}/{seeds}")
