"""
Constant-Aware Q∩T Solver.

Key discovery: 96+ quadratic products are IV-determined CONSTANTS.
By substituting these constants, quadratic equations become linear,
expanding the effective kernel and pushing the solvability boundary.

Round 1: b[1]=IV.a, c[1]=IV.b, f[1]=IV.e, g[1]=IV.f → all products constant
Round 2: c[2]=b[1]=IV.a, g[2]=f[1]=IV.e → partially constant

This turns ~50% of quad equations into linear for first 2 rounds.
"""

import random
import time

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


def _iv_constant_value(var_desc):
    """
    If this variable is an IV-determined constant, return its value.
    Otherwise return None.

    state[1].b[k] = IV[0] bit k  (b[1] = a[0] = IV.a)
    state[1].c[k] = IV[1] bit k  (c[1] = b[0] = IV.b)
    state[1].f[k] = IV[4] bit k  (f[1] = e[0] = IV.e)
    state[1].g[k] = IV[5] bit k  (g[1] = f[0] = IV.f)
    state[1].h[k] = IV[6] bit k  (h[1] = g[0] = IV.g)
    state[1].d[k] = IV[2] bit k  (d[1] = c[0] = IV.c)

    state[2].c[k] = state[1].b[k] = IV[0] bit k
    state[2].d[k] = state[1].c[k] = IV[1] bit k
    state[2].g[k] = state[1].f[k] = IV[4] bit k
    state[2].h[k] = state[1].g[k] = IV[5] bit k

    state[3].d[k] = state[2].c[k] = IV[0] bit k
    state[3].h[k] = state[2].g[k] = IV[4] bit k
    """
    # Parse "state[R].reg[bit]"
    import re
    m = re.match(r'state\[(\d+)\]\.([a-h])\[(\d+)\]', var_desc)
    if not m:
        return None

    r, reg, bit = int(m.group(1)), m.group(2), int(m.group(3))
    reg_idx = 'abcdefgh'.index(reg)

    # Shift register: b[r]=a[r-1], c[r]=b[r-1]=a[r-2], d[r]=c[r-1]=a[r-3]
    #                 f[r]=e[r-1], g[r]=f[r-1]=e[r-2], h[r]=g[r-1]=e[r-3]
    # Shift register traces back to IV:
    # b[r]=a[r-1], c[r]=b[r-1]=a[r-2], d[r]=c[r-1]=a[r-3]
    # f[r]=e[r-1], g[r]=f[r-1]=e[r-2], h[r]=g[r-1]=e[r-3]
    # If the source round is 0: the value = IV[source_reg][bit]

    # Shift register: b[r]=a[r-1], c[r]=b[r-1], d[r]=c[r-1]
    #                  f[r]=e[r-1], g[r]=f[r-1], h[r]=g[r-1]
    # Tracing back: reg at round r → reg' at round r-1 → ... → round 0
    # b[r] ← a[r-1], c[r] ← b[r-1] ← a[r-2], d[r] ← c[r-1] ← b[r-2] ← a[r-3]
    # f[r] ← e[r-1], g[r] ← f[r-1] ← e[r-2], h[r] ← g[r-1] ← f[r-2] ← e[r-3]
    # At round 0: register i = IV[i]

    # Trace back to round 0
    cur_reg = reg_idx
    cur_round = r
    while cur_round > 0:
        # One step back through shift register
        back_map = {1: 0, 2: 1, 3: 2, 5: 4, 6: 5, 7: 6}  # b←a, c←b, d←c, f←e, g←f, h←g
        if cur_reg in back_map:
            cur_reg = back_map[cur_reg]
            cur_round -= 1
        else:
            return None  # a[r] and e[r] depend on message, not just shifts
    return get_bit(IV256[cur_reg], bit)
        # For shift chain: base_reg at src_round
        # b at src_round = a at src_round-1, etc.
        # This recurses — only round 0 gives constants
        # Simplified: only if src_round == 0
    return None


def build_constant_aware_system(msg_words, num_rounds, target_hash=None):
    """
    Build Q system with constant products pre-substituted as linear equations.

    For each quadratic equation with products involving IV-constant registers:
    - If BOTH factors of a product are IV-constants: substitute constant value
    - If ONE factor is IV-constant: substitute, reducing product to linear or 0
    - The resulting equation may become fully linear
    """
    system, trace = build_qt_system(msg_words, num_rounds, target_hash)
    vm = system.vm
    n = vm.num_vars

    linear_rows = system.to_linear_rows()
    new_linear = []
    remaining_quad = []

    for eq in system.q_quadratic:
        # Process each quadratic term
        new_linear_terms = set(eq.linear)
        new_quad_terms = set()
        new_const = eq.constant

        for v1, v2 in eq.quadratic:
            d1 = vm.describe(v1)
            d2 = vm.describe(v2)
            c1 = _iv_constant_value(d1)
            c2 = _iv_constant_value(d2)

            if c1 is not None and c2 is not None:
                # Both constant: product is a constant
                new_const ^= (c1 & c2)
            elif c1 is not None:
                # v1 is constant: product = c1 * v2
                if c1 == 1:
                    new_linear_terms ^= {v2}
                # if c1 == 0: product = 0, nothing to add
            elif c2 is not None:
                # v2 is constant: product = v1 * c2
                if c2 == 1:
                    new_linear_terms ^= {v1}
            else:
                # Both variable: keep quadratic
                new_quad_terms.add((min(v1, v2), max(v1, v2)))

        if len(new_quad_terms) == 0:
            # Became fully linear!
            row = 0
            for v in new_linear_terms:
                row ^= (1 << v)
            if new_const:
                row |= (1 << n)
            new_linear.append(row)
        else:
            remaining_quad.append((new_linear_terms, new_quad_terms, new_const))

    all_linear = linear_rows + new_linear
    return all_linear, remaining_quad, system, trace


def const_aware_greedy_walk(num_rounds, max_steps=800, seed=42, verbose=True):
    """
    Greedy walk on constant-aware Q-variety.

    With constant substitution, the linear system is larger
    → smaller kernel → but fewer quad constraints to satisfy.
    """
    rng = random.Random(seed)
    msg = [rng.randint(0, MASK32) for _ in range(16)]

    if verbose:
        print(f"\n{'='*60}")
        print(f"Constant-Aware Greedy Walk: R={num_rounds}")
        print(f"{'='*60}")

    all_linear, remaining_quad, system, trace = build_constant_aware_system(
        msg, num_rounds
    )
    vm = system.vm
    n = vm.num_vars
    target = trace.hash_words
    ref_carries = get_all_carry_chains(trace)

    if verbose:
        orig_quad = len(system.q_quadratic)
        new_lin = len(all_linear) - len(system.to_linear_rows())
        print(f"Original: {orig_quad} quad eqs")
        print(f"After constant substitution: {new_lin} became linear, "
              f"{len(remaining_quad)} remain quadratic")

    # Solve enlarged linear system
    result = gf2_solve(all_linear, n)
    if result is None:
        if verbose:
            print("Linear system inconsistent!")
        return {'success': False}

    particular, kernel = result
    if verbose:
        print(f"Kernel dim: {len(kernel)} (was 480 without constants)")
        print(f"Remaining quad: {len(remaining_quad)}")
        print(f"Effective DOF: {len(kernel)} - {len(remaining_quad)} = "
              f"{len(kernel) - len(remaining_quad)}")

    # Greedy walk minimizing carry_mm + quad_violations
    current = particular

    def evaluate(bv):
        m = _extract_msg(bv, vm)
        t = sha256_compress_traced(m, num_rounds)
        c = get_all_carry_chains(t)
        cmm = sum(bin(ref_carries[i] ^ c[i]).count('1') for i in range(len(ref_carries)))

        # Quad violations using Q-predicted state
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

    cost, cmm, qviol = evaluate(current)
    best_cost, best_bv = cost, current

    if verbose:
        print(f"Start: cost={cost} (carry={cmm}, quad={qviol})")

    t0 = time.time()
    for step in range(max_steps):
        improved = False
        for kvec in kernel:
            if kvec & ((1 << 512) - 1) == 0:
                continue
            c, cm, qv = evaluate(current ^ kvec)
            if c < cost:
                cost, cmm, qviol = c, cm, qv
                current ^= kvec
                improved = True
                if c < best_cost:
                    best_cost = c
                    best_bv = current
                break

        if not improved:
            for _ in range(rng.randint(1, 5)):
                current ^= kernel[rng.randint(0, len(kernel) - 1)]
            cost, cmm, qviol = evaluate(current)

        if best_cost == 0:
            found_msg = _extract_msg(best_bv, vm)
            found_hash = sha256_compress(found_msg, num_rounds)
            if verbose:
                print(f"\n  cost=0 at step {step}!")
                print(f"  Hash match: {found_hash == target}")
                if found_hash == target and found_msg != msg:
                    print(f"  *** SECOND PREIMAGE FOUND! ***")
            return {
                'success': found_hash == target,
                'steps': step,
                'msg': found_msg,
                'kernel_dim': len(kernel),
                'remaining_quad': len(remaining_quad),
            }

        if verbose and step % 100 == 0:
            print(f"  step {step}: cost={cost} (c={cmm},q={qviol}), "
                  f"best={best_cost}, {time.time()-t0:.1f}s")

    elapsed = time.time() - t0
    if verbose:
        print(f"\nBest: cost={best_cost} ({elapsed:.1f}s)")
    return {
        'success': False,
        'best_cost': best_cost,
        'kernel_dim': len(kernel),
        'remaining_quad': len(remaining_quad),
    }


def benchmark_const_aware(max_R=8, seeds=5):
    """Compare constant-aware vs standard solver."""
    print(f"\n{'═'*60}")
    print(f"BENCHMARK: Constant-Aware vs Standard")
    print(f"{'═'*60}")
    print(f"{'R':>3} | {'Kern':>5} | {'RemQ':>4} | {'DOF':>4} | {'Success':>7}")
    print(f"{'-'*3}-+-{'-'*5}-+-{'-'*4}-+-{'-'*4}-+-{'-'*7}")

    for R in range(2, max_R + 1):
        successes = 0
        kern = 0
        remq = 0
        for s in range(seeds):
            r = const_aware_greedy_walk(R, max_steps=800, seed=s*100+42, verbose=False)
            if r.get('success'):
                successes += 1
            kern = r.get('kernel_dim', 0)
            remq = r.get('remaining_quad', 0)

        dof = kern - remq
        print(f"{R:>3} | {kern:>5} | {remq:>4} | {dof:>4} | {successes}/{seeds}")
