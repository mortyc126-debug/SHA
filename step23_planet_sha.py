"""
Planet SHA-256 — First Observations

No numbers. No algebra. Only PATHS.

A path = sequence of states through 64 layers.
Collision = two paths arriving at the same endpoint.

Questions for the planet:
1. How do paths FLOW? (traffic at each layer)
2. Are there FUNNELS? (states that attract many paths)
3. When paths MEET, what happens? (convergence dynamics)
4. Can we predict WHERE paths will end, from WHERE they are mid-way?
"""

import numpy as np
from collections import Counter, defaultdict

N = 4
MASK = (1 << N) - 1
N_MSG = 4
N_TOTAL = 1 << (N * N_MSG)

IV = [0x6, 0xB, 0x3, 0xA, 0x5, 0x9, 0x1, 0xF]
K = [0x4, 0x2, 0xB, 0x7, 0xA, 0x3, 0xE, 0x5,
     0x9, 0x1, 0xD, 0x6, 0x0, 0x8, 0xC, 0xF]


def rotr(x, s):
    return ((x >> s) | (x << (N - s))) & MASK


def full_path(msg, R):
    """Compute the FULL PATH of a message through R layers."""
    a, b, c, d, e, f, g, h = IV[:]
    W = list(msg) + [0] * max(0, R - len(msg))
    path = [(a, e)]  # Track (a, e) as the "position" at each layer

    for r in range(R):
        w_r = W[r] if r < len(W) else 0
        k_r = K[r % len(K)]
        sig1 = rotr(e, 1) ^ rotr(e, 3) ^ (e >> 1)
        ch = (e & f) ^ (~e & g) & MASK
        sig0 = rotr(a, 1) ^ rotr(a, 2) ^ rotr(a, 3)
        maj = (a & b) ^ (a & c) ^ (b & c)
        T1 = (h + sig1 + ch + k_r + w_r) & MASK
        T2 = (sig0 + maj) & MASK
        a_new = (T1 + T2) & MASK
        e_new = (d + T1) & MASK
        h, g, f = g, f, e
        e = e_new
        d, c, b = c, b, a
        a = a_new
        path.append((a, e))

    return path


def observe_planet(R=8):
    """First observations on Planet SHA-256."""

    print(f"{'='*70}")
    print(f"PLANET SHA-256 — First Observations")
    print(f"{'='*70}")
    print(f"Paths: {N_TOTAL}, Layers: {R}, States per layer: {(1<<N)**2} = {(1<<N)**2}")

    # Compute ALL paths
    all_paths = {}
    for midx in range(N_TOTAL):
        msg = []
        tmp = midx
        for w in range(N_MSG):
            msg.append(tmp & MASK)
            tmp >>= N
        all_paths[midx] = full_path(msg, R)

    # ================================================================
    # OBSERVATION 1: Traffic at each layer
    # How many paths visit each state?
    # ================================================================
    print(f"\n--- OBSERVATION 1: Traffic Flow ---")
    print(f"{'Layer':>6} {'Distinct':>9} {'Max traffic':>12} {'Min traffic':>12} {'Funnel ratio':>13}")

    for r in range(R + 1):
        traffic = Counter()
        for midx, path in all_paths.items():
            traffic[path[r]] += 1

        distinct = len(traffic)
        max_t = max(traffic.values())
        min_t = min(traffic.values())
        # Funnel ratio = max_traffic / average_traffic
        avg_t = N_TOTAL / distinct
        funnel = max_t / avg_t

        print(f"  {r:>4}  {distinct:>8}  {max_t:>11}  {min_t:>11}  {funnel:>12.2f}×")

    # ================================================================
    # OBSERVATION 2: Funnels — states that attract MANY paths
    # ================================================================
    print(f"\n--- OBSERVATION 2: Biggest Funnels ---")

    for r in [0, R//2, R-1, R]:
        traffic = Counter()
        for midx, path in all_paths.items():
            traffic[path[r]] += 1

        top5 = traffic.most_common(5)
        print(f"  Layer {r}: top states by traffic:")
        for state, count in top5:
            pct = count / N_TOTAL * 100
            print(f"    state={state} → {count} paths ({pct:.1f}%)")

    # ================================================================
    # OBSERVATION 3: Convergence dynamics
    # When two paths share a state at layer r, what happens next?
    # ================================================================
    print(f"\n--- OBSERVATION 3: Path Convergence ---")
    print(f"If two paths share a state at layer r, do they share at layer r+1?")

    for r in range(R):
        # Group paths by state at layer r
        groups = defaultdict(list)
        for midx, path in all_paths.items():
            groups[path[r]].append(midx)

        # For each group (paths sharing state at r):
        # Check if they also share state at r+1
        total_pairs = 0
        same_next = 0
        diff_next = 0

        for state, members in groups.items():
            if len(members) < 2:
                continue
            for i in range(min(len(members), 50)):
                for j in range(i+1, min(len(members), 50)):
                    total_pairs += 1
                    m1, m2 = members[i], members[j]
                    if all_paths[m1][r+1] == all_paths[m2][r+1]:
                        same_next += 1
                    else:
                        diff_next += 1

        if total_pairs > 0:
            p_stay = same_next / total_pairs
            print(f"  Layer {r}→{r+1}: {total_pairs} pairs share state at {r}. "
                  f"Stay together: {same_next} ({p_stay*100:.1f}%)")

    # ================================================================
    # OBSERVATION 4: Endpoint prediction from midpoint
    # If I know a path's state at layer r, can I predict its endpoint?
    # ================================================================
    print(f"\n--- OBSERVATION 4: Predicting the Endpoint ---")
    print(f"Knowing state at layer r, how many possible endpoints?")

    for r in range(R + 1):
        groups = defaultdict(set)
        for midx, path in all_paths.items():
            groups[path[r]].add(path[R])

        endpoint_counts = [len(ends) for ends in groups.values()]
        avg_endpoints = np.mean(endpoint_counts)
        max_endpoints = max(endpoint_counts)
        min_endpoints = min(endpoint_counts)

        print(f"  From layer {r}: avg {avg_endpoints:.1f} endpoints, "
              f"min {min_endpoints}, max {max_endpoints}")

    # ================================================================
    # OBSERVATION 5: Collision geography
    # Paths that collide (same endpoint) — where do they SEPARATE?
    # ================================================================
    print(f"\n--- OBSERVATION 5: Collision Geography ---")
    print(f"Collision pairs: where do they first DIVERGE?")

    # Group by endpoint
    endpoint_groups = defaultdict(list)
    for midx, path in all_paths.items():
        endpoint_groups[path[R]].append(midx)

    # For each collision group: find earliest divergence
    diverge_counts = Counter()
    total_col_pairs = 0

    for endpoint, members in endpoint_groups.items():
        if len(members) < 2:
            continue
        for i in range(min(len(members), 20)):
            for j in range(i+1, min(len(members), 20)):
                total_col_pairs += 1
                m1, m2 = members[i], members[j]
                # Find first layer where paths differ
                for r in range(R + 1):
                    if all_paths[m1][r] != all_paths[m2][r]:
                        diverge_counts[r] += 1
                        break

    if total_col_pairs > 0:
        print(f"  Total collision pairs sampled: {total_col_pairs}")
        for r in sorted(diverge_counts.keys()):
            pct = diverge_counts[r] / total_col_pairs * 100
            print(f"  First diverge at layer {r}: {diverge_counts[r]} pairs ({pct:.1f}%)")

    # ================================================================
    # OBSERVATION 6: Last convergence — when do collision paths MEET?
    # ================================================================
    print(f"\n--- OBSERVATION 6: Last Convergence ---")
    print(f"Collision pairs: when do they LAST share a state (before endpoint)?")

    converge_counts = Counter()

    for endpoint, members in endpoint_groups.items():
        if len(members) < 2:
            continue
        for i in range(min(len(members), 20)):
            for j in range(i+1, min(len(members), 20)):
                m1, m2 = members[i], members[j]
                # Find LAST layer before R where paths share state
                last_shared = -1
                for r in range(R):
                    if all_paths[m1][r] == all_paths[m2][r]:
                        last_shared = r
                converge_counts[last_shared] += 1

    if converge_counts:
        print(f"  Last shared state before endpoint:")
        for r in sorted(converge_counts.keys()):
            pct = converge_counts[r] / total_col_pairs * 100
            label = "(never shared)" if r == -1 else f"layer {r}"
            print(f"  {label:>20}: {converge_counts[r]} pairs ({pct:.1f}%)")


if __name__ == "__main__":
    observe_planet(R=8)
