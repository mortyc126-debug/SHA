"""
V4 Wang-Constraint Solver.

Key insight from methodology_v20.md:
  When delta_e[r] = 0: Ch(e,f,g) becomes LINEAR → 2 quad products vanish per bit.
  The Wang chain achieves delta_e[2..16] = 0 with P=1.0.

Strategy: ADD delta_e = 0 as explicit linear constraints.
  - For the Q system built from message M with carries C(M):
    state[r] variables represent the ACTUAL SHA state.
  - If we constrain e[r] to equal a SPECIFIC value (from the base message),
    then Ch products involving e[r] become constants.
  - This converts Ch-quadratic equations to linear.

Effect: each round with constrained e loses 32 quad eqs (Ch part).
  Maj still contributes 32 quad eqs.
  Net: 32 quad/round instead of 64 → DOF = kern - 32*(R-1) instead of kern - 64*(R-1).

For R=6: DOF = 512 - 32*5 = 352 (vs 256 without).
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
    """Check if variable is IV-determined constant via shift register chain."""
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


def build_wang_constrained_system(msg_words, num_rounds, target_hash=None,
                                   constrain_e_rounds=None):
    """
    Build Q system with Wang-style e-register constraints.

    For each round r in constrain_e_rounds:
      Add constraints: e[r][b] = base_e[r][b] for all b=0..31
      This fixes e[r] to the value from the base message trace.

    When e[r] is fixed (constant), Ch(e[r],f[r],g[r]) products become:
      e[r][b]*f[r][b] → const * f[r][b] = linear (or 0)
      e[r][b]*g[r][b] → const * g[r][b] = linear (or 0)

    This eliminates Ch-quadratic terms for those rounds!
    """
    system, trace = build_qt_system(msg_words, num_rounds, target_hash)
    vm = system.vm
    n = vm.num_vars

    if constrain_e_rounds is None:
        # Constrain e for rounds where it helps most: rounds 2+ (after IV)
        # Round 0: state is IV (already constant)
        # Round 1+: e depends on message
        constrain_e_rounds = list(range(1, num_rounds))

    # Get actual e-register values from trace
    e_values = {}
    for r in constrain_e_rounds:
        if r < len(trace.states):
            e_values[r] = trace.states[r][4]  # reg 4 = e

    # Start with standard linear rows
    linear_rows = list(system.to_linear_rows())

    # Add e-constraint equations: e[r][b] = known_value
    e_constraints_added = 0
    for r in constrain_e_rounds:
        if r not in e_values:
            continue
        e_val = e_values[r]
        for b in range(32):
            var = vm.s(r, 4, b)  # state[r].e[b]
            bit_val = get_bit(e_val, b)
            row = (1 << var)
            if bit_val:
                row |= (1 << n)
            linear_rows.append(row)
            e_constraints_added += 1

    # Now process quadratic equations with IV constants AND e-constraints
    remaining_quad = []
    linearized = 0

    for eq in system.q_quadratic:
        new_lin = set(eq.linear)
        new_quad = set()
        new_const = eq.constant

        for v1, v2 in eq.quadratic:
            # Check IV constants
            c1 = _iv_constant_value(vm, v1)
            c2 = _iv_constant_value(vm, v2)

            # Check e-constraints
            if c1 is None:
                d1 = vm.describe(v1)
                import re
                m1 = re.match(r'state\[(\d+)\]\.e\[(\d+)\]', d1)
                if m1:
                    r1 = int(m1.group(1))
                    b1 = int(m1.group(2))
                    if r1 in e_values:
                        c1 = get_bit(e_values[r1], b1)

            if c2 is None:
                d2 = vm.describe(v2)
                import re
                m2 = re.match(r'state\[(\d+)\]\.e\[(\d+)\]', d2)
                if m2:
                    r2 = int(m2.group(1))
                    b2 = int(m2.group(2))
                    if r2 in e_values:
                        c2 = get_bit(e_values[r2], b2)

            # Substitute known values
            if c1 is not None and c2 is not None:
                new_const ^= (c1 & c2)
            elif c1 is not None:
                if c1 == 1:
                    new_lin ^= {v2}
            elif c2 is not None:
                if c2 == 1:
                    new_lin ^= {v1}
            else:
                new_quad.add((min(v1, v2), max(v1, v2)))

        if len(new_quad) == 0:
            row = 0
            for v in new_lin:
                row ^= (1 << v)
            if new_const:
                row |= (1 << n)
            linear_rows.append(row)
            linearized += 1
        else:
            remaining_quad.append((new_lin, new_quad, new_const))

    return linear_rows, remaining_quad, system, trace, {
        'n': n,
        'e_constraints': e_constraints_added,
        'linearized': linearized,
        'constrained_rounds': constrain_e_rounds,
    }


def v4_wang_walk(num_rounds, max_steps=1500, seed=42, verbose=True):
    """
    V4 solver: greedy walk on Wang-constrained Q-variety.

    The e-constraints fix e-register values, linearizing Ch.
    The walk then only needs to satisfy Maj quadratics + carry consistency.
    """
    rng = random.Random(seed)
    msg = [rng.randint(0, MASK32) for _ in range(16)]

    if verbose:
        print(f"\n{'='*60}")
        print(f"V4 Wang-Constraint Solver: R={num_rounds}")
        print(f"{'='*60}")

    linear_rows, remaining_quad, system, trace, info = build_wang_constrained_system(
        msg, num_rounds
    )
    vm = system.vm
    n = info['n']
    target = trace.hash_words
    ref_carries = get_all_carry_chains(trace)

    if verbose:
        print(f"E-constraints added: {info['e_constraints']}")
        print(f"Quad linearized: {info['linearized']}")
        print(f"Remaining quad: {len(remaining_quad)}")

    result = gf2_solve(linear_rows, n)
    if result is None:
        if verbose:
            print("Linear system inconsistent!")
        return {'success': False, 'kernel_dim': 0, 'remaining_quad': len(remaining_quad)}

    particular, kernel = result

    if verbose:
        print(f"Kernel: {len(kernel)}")
        print(f"Effective DOF: {len(kernel) - len(remaining_quad)}")

    # Greedy + SA walk
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

    use_sa = num_rounds >= 6
    temperature = 80.0 if use_sa else 0
    cooling = 0.997

    t0 = time.time()
    for step in range(max_steps):
        improved = False

        if use_sa:
            trials = min(60, len(kernel))
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
                    break
            elif use_sa and temperature > 0.1:
                delta = c - cost
                if delta < temperature and rng.random() < math.exp(-delta / temperature):
                    cost, cmm, qviol = c, cm, qv
                    current = cand
                    improved = True

        if use_sa:
            temperature *= cooling

        if not improved:
            n_flips = rng.randint(2, min(10, len(kernel)))
            for _ in range(n_flips):
                current ^= kernel[rng.randint(0, len(kernel) - 1)]
            cost, cmm, qviol = evaluate(current)

        if best_cost == 0:
            found_msg = _extract_msg(best_bv, vm)
            found_hash = sha256_compress(found_msg, num_rounds)
            success = found_hash == target and found_msg != msg
            if verbose:
                print(f"\n  cost=0 at step {step}! Hash match: {found_hash == target}")
                if success:
                    print(f"  *** SECOND PREIMAGE FOUND! ***")
            return {
                'success': success,
                'steps': step,
                'kernel_dim': len(kernel),
                'remaining_quad': len(remaining_quad),
            }

        if verbose and step % 200 == 0:
            sa_info = f", T={temperature:.1f}" if use_sa else ""
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


def benchmark_v4(max_R=8, seeds=5):
    """Benchmark V4 solver."""
    print(f"\n{'═'*60}")
    print(f"V4 WANG-CONSTRAINT BENCHMARK")
    print(f"{'═'*60}")
    print(f"{'R':>3} | {'Kern':>5} | {'RemQ':>4} | {'DOF':>4} | {'Success':>7}")
    print(f"{'-'*3}-+-{'-'*5}-+-{'-'*4}-+-{'-'*4}-+-{'-'*7}")

    for R in range(2, max_R + 1):
        successes = 0
        kern = remq = 0
        for s in range(seeds):
            r = v4_wang_walk(R, max_steps=1500, seed=s * 100 + 42, verbose=False)
            if r.get('success'):
                successes += 1
            kern = r.get('kernel_dim', 0)
            remq = r.get('remaining_quad', 0)
        dof = kern - remq
        print(f"{R:>3} | {kern:>5} | {remq:>4} | {dof:>4} | {successes}/{seeds}")
