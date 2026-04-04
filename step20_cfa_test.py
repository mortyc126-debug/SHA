"""
Step 20: Test the Collision Flow Algebra (CFA)

THE KEY TEST: Is the collision flow COLDER than random at the output?

T[r] = entropy of state_r projected from collision set / max_entropy
If T[R] < 1 → collisions are CLUSTERED → exploitable structure
If T[R] = 1 → birthday optimal
"""

import numpy as np
from collections import Counter
import math

N = 4
MASK = (1 << N) - 1
N_MSG = 4
N_INPUT = N * N_MSG
N_TOTAL = 1 << N_INPUT

IV = [0x6, 0xB, 0x3, 0xA, 0x5, 0x9, 0x1, 0xF]
K = [0x4, 0x2, 0xB, 0x7, 0xA, 0x3, 0xE, 0x5,
     0x9, 0x1, 0xD, 0x6, 0x0, 0x8, 0xC, 0xF]


def rotr(x, s):
    return ((x >> s) | (x << (N - s))) & MASK


def sha_full_trace(msg, R):
    """Return state at EVERY round."""
    a, b, c, d, e, f, g, h = IV[:]
    W = list(msg) + [0] * max(0, R - N_MSG)
    trace = [(a, b, c, d, e, f, g, h)]

    for r in range(R):
        w_r = W[r] if r < len(W) else 0
        k_r = K[r % len(K)]
        sig1 = rotr(e, 1) ^ rotr(e, 3) ^ (e >> 1)
        ch_val = (e & f) ^ (~e & g) & MASK
        sig0 = rotr(a, 1) ^ rotr(a, 2) ^ rotr(a, 3)
        maj_val = (a & b) ^ (a & c) ^ (b & c)
        T1 = (h + sig1 + ch_val + k_r + w_r) & MASK
        T2 = (sig0 + maj_val) & MASK
        a_new = (T1 + T2) & MASK
        e_new = (d + T1) & MASK
        h, g, f = g, f, e
        e = e_new
        d, c, b = c, b, a
        a = a_new
        trace.append((a, b, c, d, e, f, g, h))

    return trace


def entropy(counts, total):
    """Shannon entropy in bits."""
    H = 0
    for c in counts.values():
        if c > 0:
            p = c / total
            H -= p * math.log2(p)
    return H


def test_cfa(M_base, R):
    """Test Collision Flow Algebra predictions."""
    print(f"\n{'='*80}")
    print(f"COLLISION FLOW ALGEBRA TEST — R={R}")
    print(f"{'='*80}")

    # Compute base trace
    base_trace = sha_full_trace(M_base, R)
    base_output = (base_trace[R][0], base_trace[R][4])  # (a, e) at final round

    # Find ALL collisions and their traces
    collision_traces = []
    for didx in range(1, N_TOTAL):
        dw = []
        tmp = didx
        for w in range(N_MSG):
            dw.append(tmp & MASK)
            tmp >>= N
        M2 = [(M_base[w] ^ dw[w]) for w in range(N_MSG)]
        trace2 = sha_full_trace(M2, R)
        output2 = (trace2[R][0], trace2[R][4])
        if output2 == base_output:
            collision_traces.append(trace2)

    n_col = len(collision_traces)
    print(f"  Collisions: {n_col}")

    if n_col == 0:
        print(f"  No collisions found!")
        return

    # Generate RANDOM subset of same size (control)
    import random
    rng = random.Random(42)
    random_traces = []
    random_indices = rng.sample(range(1, N_TOTAL), min(n_col, N_TOTAL - 1))
    for didx in random_indices:
        dw = []
        tmp = didx
        for w in range(N_MSG):
            dw.append(tmp & MASK)
            tmp >>= N
        M2 = [(M_base[w] ^ dw[w]) for w in range(N_MSG)]
        random_traces.append(sha_full_trace(M2, R))

    # Measure TEMPERATURE at each round
    max_entropy_ae = 2 * N  # 8 bits for (a, e)
    max_entropy_full = 8 * N  # 32 bits for full state

    print(f"\n  TEMPERATURE PROFILE (a,e projection, max entropy = {max_entropy_ae} bits):")
    print(f"  {'Round':>6} {'T_collision':>12} {'T_random':>10} {'Ratio':>7} {'Signal':>8}")

    temperatures_col = []
    temperatures_rnd = []

    for r in range(R + 1):
        # Collision flow: distribution of (a,e) at round r
        col_states = Counter()
        for trace in collision_traces:
            ae = (trace[r][0], trace[r][4])  # (a, e)
            col_states[ae] += 1

        rnd_states = Counter()
        for trace in random_traces:
            ae = (trace[r][0], trace[r][4])
            rnd_states[ae] += 1

        H_col = entropy(col_states, n_col)
        H_rnd = entropy(rnd_states, len(random_traces))

        T_col = H_col / max_entropy_ae
        T_rnd = H_rnd / max_entropy_ae

        temperatures_col.append(T_col)
        temperatures_rnd.append(T_rnd)

        ratio = T_col / max(T_rnd, 0.001)
        signal = "★" if ratio < 0.95 else ("★★" if ratio < 0.8 else "")

        print(f"  R={r:>3}  {T_col:>11.4f}  {T_rnd:>9.4f}  {ratio:>6.3f}  {signal:>7}")

    # KEY RESULT: Temperature at output
    print(f"\n  ★ OUTPUT TEMPERATURE: T_col={temperatures_col[R]:.4f}, T_rnd={temperatures_rnd[R]:.4f}")
    print(f"  ★ Ratio: {temperatures_col[R]/max(temperatures_rnd[R], 0.001):.4f}")

    if temperatures_col[R] < temperatures_rnd[R] * 0.95:
        print(f"  ★★★ COLLISION FLOW IS COLDER! Structure exists at output!")
    else:
        print(f"  Collision flow ≈ random at output.")

    # CORRELATION between rounds
    print(f"\n  INTER-ROUND CORRELATION (within collision flow):")
    print(f"  {'R1':>4} {'R2':>4} {'Corr_col':>10} {'Corr_rnd':>10} {'Excess':>8}")

    for r1 in [0, 2, 4]:
        for r2 in [R-2, R-1, R]:
            if r1 >= r2:
                continue

            # Correlation: P(same (a,e) at r2 | same (a,e) at r1)
            # For collision set: track (state_r1, state_r2) pairs
            col_pairs = Counter()
            for trace in collision_traces:
                ae1 = (trace[r1][0], trace[r1][4])
                ae2 = (trace[r2][0], trace[r2][4])
                col_pairs[(ae1, ae2)] += 1

            rnd_pairs = Counter()
            for trace in random_traces:
                ae1 = (trace[r1][0], trace[r1][4])
                ae2 = (trace[r2][0], trace[r2][4])
                rnd_pairs[(ae1, ae2)] += 1

            # Mutual information I(state_r1; state_r2)
            H_joint_col = entropy(col_pairs, n_col)
            H_r1_col = entropy(
                Counter(trace[r1][0] << N | trace[r1][4] for trace in collision_traces),
                n_col)
            H_r2_col = entropy(
                Counter(trace[r2][0] << N | trace[r2][4] for trace in collision_traces),
                n_col)
            MI_col = H_r1_col + H_r2_col - H_joint_col

            H_joint_rnd = entropy(rnd_pairs, len(random_traces))
            H_r1_rnd = entropy(
                Counter(trace[r1][0] << N | trace[r1][4] for trace in random_traces),
                len(random_traces))
            H_r2_rnd = entropy(
                Counter(trace[r2][0] << N | trace[r2][4] for trace in random_traces),
                len(random_traces))
            MI_rnd = H_r1_rnd + H_r2_rnd - H_joint_rnd

            excess = MI_col - MI_rnd
            signal = " ★" if excess > 0.1 else ""
            print(f"  {r1:>4} {r2:>4} {MI_col:>10.4f} {MI_rnd:>10.4f} {excess:>+7.4f}{signal}")

    # DENSITY at key rounds
    print(f"\n  DENSITY PROFILE (how many distinct states does collision flow visit?):")
    print(f"  {'Round':>6} {'Distinct_col':>13} {'Distinct_rnd':>13} {'Max_possible':>13}")

    for r in [0, 2, 4, R//2, R-2, R-1, R]:
        if r > R:
            continue
        col_distinct = len(set((trace[r][0], trace[r][4]) for trace in collision_traces))
        rnd_distinct = len(set((trace[r][0], trace[r][4]) for trace in random_traces))
        max_possible = (1 << N) ** 2  # 256

        print(f"  R={r:>3}  {col_distinct:>12}  {rnd_distinct:>12}  {max_possible:>12}")


def main():
    import random
    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(N_MSG)]

    for R in [4, 6, 8]:
        test_cfa(M, R)


if __name__ == "__main__":
    main()
