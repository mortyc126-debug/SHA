"""
Wang+Q∩T Hybrid Solver.

Injects Wang chain's e,f,g values as constants → Ch linearized.
Only Maj quadratics remain. α-kernel > 0 at R≥8.
"""

import random, time, re, math
from qt_solver.sha256_traced import (
    MASK32, IV256, get_bit, sha256_compress, sha256_compress_traced,
    get_all_carry_chains,
)
from qt_solver.gf2 import gf2_solve, gf2_rank
from qt_solver.qt_system import build_qt_system
from qt_solver.v3_solver import _iv_constant_value
from qt_solver.v10_from_M import _encode_M, _extract_msg


def _wang_constant(vm, var_idx, e_values):
    """Check if var is known from IV chain OR Wang e-chain."""
    # IV constant?
    c = _iv_constant_value(vm, var_idx)
    if c is not None:
        return c

    d = vm.describe(var_idx)

    # e[r] directly
    m = re.match(r'state\[(\d+)\]\.e\[(\d+)\]', d)
    if m:
        r, b = int(m.group(1)), int(m.group(2))
        return e_values.get((r, b))

    # f[r] = e[r-1]
    m = re.match(r'state\[(\d+)\]\.f\[(\d+)\]', d)
    if m:
        r, b = int(m.group(1)), int(m.group(2))
        return e_values.get((r - 1, b))

    # g[r] = e[r-2] (= f[r-1] = e[r-2])
    m = re.match(r'state\[(\d+)\]\.g\[(\d+)\]', d)
    if m:
        r, b = int(m.group(1)), int(m.group(2))
        return e_values.get((r - 2, b))

    # h[r] = e[r-3] — DO NOT inject as constant!
    # Over-constraining h kills α-kernel (tested: kernel→0).
    # Leave h as variable to preserve algebraic freedom.

    return None


def build_hybrid_system(msg, R):
    """Build Q∩T with Wang e,f,g,h constants injected."""
    system, trace = build_qt_system(msg, R)
    vm = system.vm
    n = vm.num_vars

    # Extract e-values from trace (Wang's known constants)
    e_values = {}
    for r in range(1, R + 1):
        for b in range(32):
            e_values[(r, b)] = get_bit(trace.states[r][4], b)

    linear_rows = list(system.to_linear_rows())
    remaining_quad = []

    for eq in system.q_quadratic:
        new_lin = set(eq.linear)
        new_quad = set()
        new_const = eq.constant

        for v1, v2 in eq.quadratic:
            c1 = _wang_constant(vm, v1, e_values)
            c2 = _wang_constant(vm, v2, e_values)

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

        if not new_quad:
            row = 0
            for v in new_lin:
                row ^= (1 << v)
            if new_const:
                row |= (1 << n)
            linear_rows.append(row)
        else:
            remaining_quad.append((new_lin, new_quad, new_const))

    return linear_rows, remaining_quad, system, trace, n


def hybrid_solve(R, max_steps=3000, seed=42, verbose=True):
    """Wang+Q∩T hybrid greedy solver."""
    rng = random.Random(seed)
    msg = [rng.randint(0, MASK32) for _ in range(16)]

    if verbose:
        print(f"\n{'='*60}")
        print(f"Wang+Q∩T Hybrid: R={R}")
        print(f"{'='*60}")

    linear_rows, remaining_quad, system, trace, n = build_hybrid_system(msg, R)
    vm = system.vm
    target = trace.hash_words
    ref_carries = get_all_carry_chains(trace)

    result = gf2_solve(linear_rows, n)
    if not result:
        if verbose:
            print("Inconsistent!")
        return {'success': False}
    particular, kernel = result

    rel_mask = 0
    for w in range(R):
        for b in range(32):
            rel_mask |= (1 << vm.w(w, b))
    filtered = [k for k in kernel if (k & rel_mask) != 0]

    if verbose:
        print(f"Remaining quad (Maj only): {len(remaining_quad)}")
        print(f"Filtered kernel: {len(filtered)}")
        print(f"Effective DOF: {len(filtered) - len(remaining_quad)}")

    if not filtered:
        if verbose:
            print("No filtered kernel vectors!")
        return {'success': False}

    M_bv = _encode_M(msg, trace, vm, R)

    def evaluate(bv):
        m = _extract_msg(bv, vm)
        t = sha256_compress_traced(m, R)
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
            if val:
                qv += 1
        return cmm + qv * 5, cmm, qv

    # Start from best 1-hop AWAY from M (force M'≠M)
    t0 = time.time()
    best_cost = 99999
    best_bv = None

    for kvec in filtered:
        cand = M_bv ^ kvec
        cand_msg = _extract_msg(cand, vm)
        if cand_msg == msg:
            continue  # Skip trivial (back to M)
        c, cm, qv = evaluate(cand)
        if c < best_cost:
            best_cost = c
            best_bv = cand

    if verbose:
        print(f"Best non-M 1-hop: cost={best_cost} ({time.time()-t0:.1f}s)")

    if best_cost == 0 and best_bv is not None:
        found = _extract_msg(best_bv, vm)
        h = sha256_compress(found, R)
        if h == target and found != msg:
            if verbose:
                print(f"*** R={R} SECOND PREIMAGE (1-hop)! ***")
            return {'success': True, 'hops': 1}

    # Phase 2: greedy+SA walk from best point
    current = best_bv
    cost, cmm, qviol = evaluate(current)
    temperature = 80.0
    cooling = 0.998

    for step in range(max_steps):
        improved = False
        trials = min(len(filtered), 80)
        indices = rng.sample(range(len(filtered)), trials)

        for ki in indices:
            kvec = filtered[ki]
            cand = current ^ kvec
            c, cm, qv = evaluate(cand)
            if c < cost:
                cost, cmm, qviol = c, cm, qv
                current = cand
                improved = True
                if c < best_cost:
                    best_cost = c
                    best_bv = cand
            elif temperature > 0.1:
                delta = c - cost
                if delta < temperature and rng.random() < math.exp(-delta / temperature):
                    cost, cmm, qviol = c, cm, qv
                    current = cand
                    improved = True

        temperature *= cooling

        if not improved:
            nf = rng.randint(1, max(2, min(6, len(filtered))))
            for _ in range(nf):
                current ^= filtered[rng.randint(0, len(filtered) - 1)]
            cost, cmm, qviol = evaluate(current)

        if best_cost == 0:
            found = _extract_msg(best_bv, vm)
            if found == msg:
                # Returned to M trivially — skip and keep searching
                best_cost = 1  # Reset to force continued search
                continue
            h = sha256_compress(found, R)
            success = h == target
            if verbose:
                print(f"\n  cost=0 at step {step}! Hash: {h == target}, diff_msg=True")
                if success:
                    diffs = sum(1 for a, b in zip(msg, found) if a != b)
                    print(f"  *** R={R} SECOND PREIMAGE! {diffs}/16 words differ ***")
            return {'success': success, 'steps': step, 'best_cost': 0}

        if verbose and step % 500 == 0:
            print(f"  step {step}: cost={cost} (c={cmm},q={qviol}), "
                  f"best={best_cost}, T={temperature:.1f}, {time.time()-t0:.1f}s")

    if verbose:
        print(f"\nBest: cost={best_cost}")
    return {'success': False, 'best_cost': best_cost}


def benchmark_hybrid(max_R=12, seeds=3):
    print(f"\n{'═'*60}")
    print(f"WANG+Q∩T HYBRID BENCHMARK")
    print(f"{'═'*60}")
    for R in range(2, max_R + 1):
        successes = 0
        best = 9999
        for s in range(seeds):
            r = hybrid_solve(R, max_steps=2000, seed=s * 100 + 42, verbose=False)
            if r.get('success'):
                successes += 1
            bc = r.get('best_cost', 9999)
            if bc < best:
                best = bc
        print(f"R={R:>2}: {successes}/{seeds}, best_cost={best}")
