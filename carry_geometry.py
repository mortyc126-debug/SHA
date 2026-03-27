#!/usr/bin/env python3
"""
Carry Geometry — Stage 1.2: Topological structure of carry-manifold
====================================================================
Three experiments:

  G1: Connectivity graph of carry profiles
      - Build adjacency via single-bit flips of W[0]
      - Connected components, diameter, hub structure

  G2: Cluster structure via carry-profile similarity
      - Hamming-based clustering
      - Carry "islands" — groups of profiles close in Hamming distance

  G3: THE KEY TEST — carry clusters vs H[7] clusters
      - Do inputs with similar carry profiles produce similar H[7]?
      - If yes → structure penetrates chaotic zone
      - If no  → carry-geometry is internal only

Based on Stage 1 findings: carry lives in ~2^10, rank 10, no XOR closure.
"""

import numpy as np
from collections import defaultdict, Counter
import time
import sys

# ============================================================
# SHA-256 CORE (minimal, from sha256_carry_algebra.py)
# ============================================================

MASK = 0xFFFFFFFF

H0 = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
]

K = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
    0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
    0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
    0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
    0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
    0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
]

def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK
def sig0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sig1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)
def Sig0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sig1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)

def message_schedule(W16):
    W = list(W16) + [0] * 48
    for i in range(16, 64):
        W[i] = (sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16]) & MASK
    return W

def sha256_full(W16):
    """Returns (carry_profile_tuple, H[7], H[6], full_H)"""
    W = message_schedule(W16)
    a, b, c, d, e, f, g, h = H0
    carries = []

    for r in range(64):
        raw = h + Sig1(e) + Ch(e, f, g) + K[r] + W[r]
        T1 = raw & MASK
        carry = 1 if raw >= (1 << 32) else 0
        carries.append(carry)

        T2 = (Sig0(a) + Maj(a, b, c)) & MASK
        h, g, f = g, f, e
        e = (d + T1) & MASK
        d, c, b = c, b, a
        a = (T1 + T2) & MASK

    H_out = [(v + iv) & MASK for v, iv in zip([a,b,c,d,e,f,g,h], H0)]
    return tuple(carries), H_out[7], H_out[6], H_out


# ============================================================
# Helper: compact carry profile to only variable rounds
# ============================================================

def find_variable_rounds(N=5000, seed=42):
    """Identify rounds where carry is not fixed (P between 0.005 and 0.995)."""
    rng = np.random.RandomState(seed)
    counts = np.zeros(64)
    for _ in range(N):
        W16 = [int(rng.randint(0, 1 << 32)) for _ in range(16)]
        cp, _, _, _ = sha256_full(W16)
        for r in range(64):
            counts[r] += cp[r]
    freq = counts / N
    variable = [r for r in range(64) if 0.005 < freq[r] < 0.995]
    return variable, freq


# ============================================================
# G1: Connectivity Graph
# ============================================================

def experiment_G1(N=8000, seed=50):
    """
    Build connectivity graph of carry profiles via single-bit W[0] flips.

    For each W[0], compute carry profile. Then flip each of 32 bits
    and check if the new profile is different (edge) or same (neutral bit).

    Measures:
      - Number of distinct neighbors per node
      - Connected components (via union-find on sampled graph)
      - Hub structure (nodes with many neighbors)
    """
    print("=" * 70)
    print("G1: Carry Profile Connectivity Graph")
    print(f"N={N}, seed={seed}")
    print("=" * 70)

    rng = np.random.RandomState(seed)

    # Phase 1: collect profiles and their neighbors
    profile_to_id = {}
    id_to_profile = {}
    edges = defaultdict(set)  # id -> set of neighbor ids
    neutral_counts = []  # how many bits are neutral per W[0]
    w0_to_id = {}
    next_id = 0

    t0 = time.time()
    for i in range(N):
        W0 = int(rng.randint(0, 1 << 32))
        W16 = [W0] + [0] * 15

        cp, h7, h6, _ = sha256_full(W16)

        if cp not in profile_to_id:
            profile_to_id[cp] = next_id
            id_to_profile[next_id] = cp
            next_id += 1
        my_id = profile_to_id[cp]
        w0_to_id[W0] = my_id

        n_neutral = 0
        for b in range(32):
            W16_flip = [W0 ^ (1 << b)] + [0] * 15
            cp_flip, _, _, _ = sha256_full(W16_flip)

            if cp_flip == cp:
                n_neutral += 1
            else:
                if cp_flip not in profile_to_id:
                    profile_to_id[cp_flip] = next_id
                    id_to_profile[next_id] = cp_flip
                    next_id += 1
                nb_id = profile_to_id[cp_flip]
                edges[my_id].add(nb_id)
                edges[nb_id].add(my_id)

        neutral_counts.append(n_neutral)

        if (i + 1) % 2000 == 0:
            elapsed = time.time() - t0
            print(f"  {i+1}/{N} profiles, {next_id} unique, {elapsed:.1f}s")

    elapsed = time.time() - t0
    n_nodes = next_id
    n_edges = sum(len(v) for v in edges.values()) // 2

    print(f"\nCompleted: {elapsed:.1f}s")
    print(f"  Nodes (unique profiles): {n_nodes}")
    print(f"  Edges (1-bit flip connections): {n_edges}")
    print(f"  Avg degree: {2*n_edges/max(n_nodes,1):.2f}")

    # Neutral bits stats
    print(f"\n--- Neutral Bits ---")
    print(f"  Mean neutral bits per W[0]: {np.mean(neutral_counts):.2f} / 32")
    print(f"  Median: {np.median(neutral_counts):.0f}")
    print(f"  Max: {max(neutral_counts)}")
    nc = Counter(neutral_counts)
    for k in sorted(nc.keys()):
        if nc[k] > N * 0.01:
            print(f"    {k} neutral: {nc[k]} ({nc[k]/N*100:.1f}%)")

    # Connected components via BFS
    visited = set()
    components = []
    for node in range(n_nodes):
        if node in visited:
            continue
        # BFS
        component = set()
        queue = [node]
        while queue:
            v = queue.pop()
            if v in visited:
                continue
            visited.add(v)
            component.add(v)
            for nb in edges.get(v, []):
                if nb not in visited:
                    queue.append(nb)
        components.append(component)

    components.sort(key=len, reverse=True)
    print(f"\n--- Connected Components ---")
    print(f"  Total components: {len(components)}")
    for i, comp in enumerate(components[:10]):
        print(f"    Component {i}: {len(comp)} nodes ({len(comp)/n_nodes*100:.1f}%)")
    if len(components) > 10:
        isolated = sum(1 for c in components if len(c) == 1)
        print(f"    ... {len(components)-10} more (isolated: {isolated})")

    # Degree distribution
    degrees = [len(edges.get(i, set())) for i in range(n_nodes)]
    print(f"\n--- Degree Distribution ---")
    print(f"  Mean degree: {np.mean(degrees):.2f}")
    print(f"  Max degree: {max(degrees)} (hub)")
    print(f"  Isolated (degree 0): {sum(1 for d in degrees if d == 0)}")
    deg_hist = Counter(degrees)
    for d in sorted(deg_hist.keys())[:15]:
        bar = '#' * min(50, deg_hist[d] * 50 // max(deg_hist.values()))
        print(f"    degree {d:3d}: {deg_hist[d]:5d} {bar}")

    return profile_to_id, edges, components, neutral_counts, n_nodes


# ============================================================
# G2: Cluster Structure
# ============================================================

def experiment_G2(N=10000, seed=51):
    """
    Cluster structure of carry profiles in Hamming space.

    Measures:
      - Pairwise Hamming distance distribution
      - Clusters: groups of profiles within distance <=2
      - "Islands": connected components at distance threshold
      - Carry-sum vs cluster membership
    """
    print("\n" + "=" * 70)
    print("G2: Carry Profile Cluster Structure")
    print(f"N={N}, seed={seed}")
    print("=" * 70)

    rng = np.random.RandomState(seed)

    # Collect profiles with their W[0] and H[7]
    data = []  # list of (W0, carry_profile, H7, H6, carry_sum)
    t0 = time.time()

    for i in range(N):
        W0 = int(rng.randint(0, 1 << 32))
        W16 = [W0] + [0] * 15
        cp, h7, h6, _ = sha256_full(W16)
        cs = sum(cp)
        data.append((W0, cp, h7, h6, cs))

        if (i + 1) % 5000 == 0:
            print(f"  {i+1}/{N} ({time.time()-t0:.1f}s)")

    elapsed = time.time() - t0
    print(f"Completed: {elapsed:.1f}s")

    # Unique profiles
    unique_profiles = set(d[1] for d in data)
    print(f"\n--- Profile Space ---")
    print(f"  Unique profiles: {len(unique_profiles)}")
    print(f"  log2(unique): {np.log2(max(len(unique_profiles), 1)):.2f}")

    # Group by profile
    profile_groups = defaultdict(list)
    for W0, cp, h7, h6, cs in data:
        profile_groups[cp].append((W0, h7, h6))

    group_sizes = [len(v) for v in profile_groups.values()]
    print(f"  Avg group size: {np.mean(group_sizes):.2f}")
    print(f"  Max group size: {max(group_sizes)}")
    print(f"  Singletons: {sum(1 for s in group_sizes if s == 1)}")

    # Pairwise Hamming distances (sample)
    profiles_list = list(unique_profiles)
    n_dist = min(len(profiles_list), 2000)
    sample_profiles = profiles_list[:n_dist]

    print(f"\n--- Pairwise Hamming Distances (sample {n_dist}) ---")
    distances = []
    for i in range(min(5000, n_dist * (n_dist - 1) // 2)):
        a = i % n_dist
        b = (i * 7 + 13) % n_dist
        if a == b:
            continue
        d = sum(sample_profiles[a][r] ^ sample_profiles[b][r] for r in range(64))
        distances.append(d)

    if distances:
        print(f"  Mean Hamming distance: {np.mean(distances):.2f}")
        print(f"  Min: {min(distances)}")
        print(f"  Max: {max(distances)}")
        # Distribution
        dist_hist = Counter(distances)
        for d in sorted(dist_hist.keys())[:10]:
            print(f"    distance {d}: {dist_hist[d]} ({dist_hist[d]/len(distances)*100:.1f}%)")

    # Clusters: profiles within Hamming distance <=2
    print(f"\n--- Cluster Analysis (Hamming threshold = 2) ---")
    # For efficiency, only cluster the first 500 profiles
    cluster_n = min(500, len(profiles_list))
    cluster_profiles = profiles_list[:cluster_n]

    # Union-Find for clustering
    parent = list(range(cluster_n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py

    for i in range(cluster_n):
        for j in range(i + 1, cluster_n):
            d = sum(cluster_profiles[i][r] ^ cluster_profiles[j][r] for r in range(64))
            if d <= 2:
                union(i, j)

    cluster_ids = defaultdict(list)
    for i in range(cluster_n):
        cluster_ids[find(i)].append(i)

    cluster_sizes = sorted([len(v) for v in cluster_ids.values()], reverse=True)
    print(f"  Profiles clustered: {cluster_n}")
    print(f"  Clusters (d<=2): {len(cluster_ids)}")
    print(f"  Largest cluster: {cluster_sizes[0] if cluster_sizes else 0}")
    print(f"  Singletons: {sum(1 for s in cluster_sizes if s == 1)}")

    return data, profile_groups, unique_profiles


# ============================================================
# G3: THE KEY TEST — Carry Clusters vs H[7]
# ============================================================

def experiment_G3(N=15000, seed=52):
    """
    THE CRITICAL EXPERIMENT: Do similar carry profiles → similar H[7]?

    If carry-geometry correlates with output geometry,
    we have found structure that penetrates the chaotic zone.

    Tests:
      1. For pairs with IDENTICAL carry profiles: how close are their H[7]?
      2. For pairs with Hamming distance 1 in carry: how close are their H[7]?
      3. For random pairs: what's the baseline H[7] distance?
      4. Correlation: Hamming(carry_A, carry_B) vs Hamming(H7_A, H7_B)
    """
    print("\n" + "=" * 70)
    print("G3: THE KEY TEST — Carry Profile vs H[7] Correlation")
    print(f"N={N}, seed={seed}")
    print("=" * 70)

    rng = np.random.RandomState(seed)

    # Collect data
    data = []
    t0 = time.time()
    for i in range(N):
        W0 = int(rng.randint(0, 1 << 32))
        W16 = [W0] + [0] * 15
        cp, h7, h6, H_full = sha256_full(W16)
        data.append((W0, cp, h7, h6, H_full))

        if (i + 1) % 5000 == 0:
            print(f"  {i+1}/{N} ({time.time()-t0:.1f}s)")

    elapsed = time.time() - t0
    print(f"Collected: {elapsed:.1f}s")

    # Group by carry profile
    profile_groups = defaultdict(list)
    for W0, cp, h7, h6, H_full in data:
        profile_groups[cp].append((W0, h7, h6, H_full))

    multi_groups = {k: v for k, v in profile_groups.items() if len(v) >= 2}
    print(f"\n  Unique profiles: {len(profile_groups)}")
    print(f"  Profiles with >=2 entries: {len(multi_groups)}")
    total_same_pairs = sum(len(v) * (len(v) - 1) // 2 for v in multi_groups.values())
    print(f"  Total same-profile pairs: {total_same_pairs}")

    # Test 1: H[7] distance for IDENTICAL carry profiles
    print(f"\n--- Test 1: Identical Carry → H[7] Distance ---")
    same_h7_dists = []
    same_h7_exact = 0
    same_h6_dists = []
    same_full_dists = []

    for cp, members in multi_groups.items():
        for i in range(len(members)):
            for j in range(i + 1, min(len(members), i + 20)):  # limit pairs per group
                h7_a, h7_b = members[i][1], members[j][1]
                h6_a, h6_b = members[i][2], members[j][2]
                H_a, H_b = members[i][3], members[j][3]

                d_h7 = bin(h7_a ^ h7_b).count('1')
                d_h6 = bin(h6_a ^ h6_b).count('1')
                d_full = sum(bin(H_a[w] ^ H_b[w]).count('1') for w in range(8))

                same_h7_dists.append(d_h7)
                same_h6_dists.append(d_h6)
                same_full_dists.append(d_full)
                if h7_a == h7_b:
                    same_h7_exact += 1

    if same_h7_dists:
        print(f"  Pairs measured: {len(same_h7_dists)}")
        print(f"  H[7] HW distance:  mean={np.mean(same_h7_dists):.2f}, min={min(same_h7_dists)}, max={max(same_h7_dists)}")
        print(f"  H[6] HW distance:  mean={np.mean(same_h6_dists):.2f}")
        print(f"  Full H distance:   mean={np.mean(same_full_dists):.2f}")
        print(f"  H[7] exact match:  {same_h7_exact} ({same_h7_exact/len(same_h7_dists)*100:.2f}%)")
    else:
        print(f"  No pairs with identical carry profiles found.")

    # Test 2: Random baseline
    print(f"\n--- Test 2: Random Baseline ---")
    random_h7_dists = []
    random_full_dists = []
    n_random = min(10000, N * (N - 1) // 2)

    rng2 = np.random.RandomState(seed + 200)
    for _ in range(n_random):
        i = rng2.randint(0, len(data))
        j = rng2.randint(0, len(data))
        if i == j:
            continue
        h7_a, h7_b = data[i][2], data[j][2]
        H_a, H_b = data[i][4], data[j][4]
        d_h7 = bin(h7_a ^ h7_b).count('1')
        d_full = sum(bin(H_a[w] ^ H_b[w]).count('1') for w in range(8))
        random_h7_dists.append(d_h7)
        random_full_dists.append(d_full)

    print(f"  Random H[7] distance: mean={np.mean(random_h7_dists):.2f}")
    print(f"  Random Full H distance: mean={np.mean(random_full_dists):.2f}")

    # The key comparison
    if same_h7_dists:
        delta_h7 = np.mean(random_h7_dists) - np.mean(same_h7_dists)
        ratio_h7 = np.mean(same_h7_dists) / np.mean(random_h7_dists)
        print(f"\n  *** KEY RESULT ***")
        print(f"  Same-carry H[7] distance:   {np.mean(same_h7_dists):.3f}")
        print(f"  Random H[7] distance:       {np.mean(random_h7_dists):.3f}")
        print(f"  Delta:                       {delta_h7:+.3f}")
        print(f"  Ratio (same/random):         {ratio_h7:.4f}")
        if ratio_h7 < 0.95:
            print(f"  >>> SIGNAL DETECTED: same carry profiles produce closer H[7] <<<")
        elif ratio_h7 > 1.05:
            print(f"  >>> ANTI-SIGNAL: same carry profiles produce MORE distant H[7] <<<")
        else:
            print(f"  >>> NO SIGNAL: carry profile does not predict H[7] proximity <<<")

    # Test 3: Carry Hamming distance vs H[7] Hamming distance correlation
    print(f"\n--- Test 3: Carry Distance vs H[7] Distance Correlation ---")
    carry_dists_sample = []
    h7_dists_sample = []
    n_corr = min(20000, N * (N - 1) // 2)

    for _ in range(n_corr):
        i = rng2.randint(0, len(data))
        j = rng2.randint(0, len(data))
        if i == j:
            continue
        cp_a, cp_b = data[i][1], data[j][1]
        h7_a, h7_b = data[i][2], data[j][2]

        d_carry = sum(cp_a[r] ^ cp_b[r] for r in range(64))
        d_h7 = bin(h7_a ^ h7_b).count('1')

        carry_dists_sample.append(d_carry)
        h7_dists_sample.append(d_h7)

    carry_arr = np.array(carry_dists_sample, dtype=float)
    h7_arr = np.array(h7_dists_sample, dtype=float)

    if len(carry_arr) > 10 and np.std(carry_arr) > 0 and np.std(h7_arr) > 0:
        corr = np.corrcoef(carry_arr, h7_arr)[0, 1]
    else:
        corr = 0.0

    print(f"  Pairs: {len(carry_dists_sample)}")
    print(f"  corr(Hamming_carry, Hamming_H7) = {corr:+.5f}")
    print(f"  Interpretation:")
    if abs(corr) > 0.05:
        print(f"    >>> CORRELATION DETECTED ({corr:+.3f}): carry geometry affects output <<<")
    elif abs(corr) > 0.01:
        print(f"    Weak signal ({corr:+.3f}): borderline, needs larger N")
    else:
        print(f"    No correlation: carry and H[7] geometries are independent")

    # Test 4: Stratified analysis — group by carry distance, measure H[7] distance
    print(f"\n--- Test 4: H[7] Distance Stratified by Carry Distance ---")
    strata = defaultdict(list)
    for cd, hd in zip(carry_dists_sample, h7_dists_sample):
        strata[cd].append(hd)

    print(f"  {'carry_dist':>11} | {'N_pairs':>8} | {'mean H7 dist':>13} | {'delta vs 16':>12}")
    print(f"  {'-'*11}-+-{'-'*8}-+-{'-'*13}-+-{'-'*12}")
    baseline_h7 = np.mean(h7_arr)
    for cd in sorted(strata.keys()):
        if len(strata[cd]) >= 20:
            mean_hd = np.mean(strata[cd])
            delta = mean_hd - baseline_h7
            marker = " ***" if abs(delta) > 0.3 else ""
            print(f"  {cd:11d} | {len(strata[cd]):8d} | {mean_hd:13.3f} | {delta:+12.3f}{marker}")

    # Test 5: Same carry profile → H[7] XOR analysis
    print(f"\n--- Test 5: H[7] XOR for Same Carry Profiles ---")
    if same_h7_dists:
        xor_h7_values = []
        for cp, members in multi_groups.items():
            for i in range(len(members)):
                for j in range(i+1, min(len(members), i+10)):
                    xor_val = members[i][1] ^ members[j][1]
                    xor_h7_values.append(xor_val)

        if xor_h7_values:
            # Are there any collisions (XOR = 0)?
            n_collisions = sum(1 for x in xor_h7_values if x == 0)
            print(f"  H[7] XOR pairs: {len(xor_h7_values)}")
            print(f"  H[7] exact collisions: {n_collisions}")

            # Distribution of HW(XOR)
            hw_xor = [bin(x).count('1') for x in xor_h7_values]
            print(f"  Mean HW(H7_a XOR H7_b): {np.mean(hw_xor):.2f} (random=16.0)")

            # Most common H[7] XOR values (if any structure)
            xor_counts = Counter(xor_h7_values)
            top_xor = xor_counts.most_common(5)
            print(f"  Top-5 H[7] XOR values:")
            for val, cnt in top_xor:
                print(f"    0x{val:08x} (HW={bin(val).count('1'):2d}): {cnt}")

    return data, profile_groups, multi_groups


# ============================================================
# MAIN
# ============================================================

if __name__ == '__main__':
    print("Carry Geometry — Stage 1.2")
    print(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print()

    N1 = 5000
    N2 = 8000
    N3 = 12000

    if len(sys.argv) > 1 and sys.argv[1] == '--fast':
        N1, N2, N3 = 2000, 3000, 5000

    t_start = time.time()

    profile_to_id, edges, components, neutral_counts, n_nodes = experiment_G1(N=N1)
    data_g2, profile_groups_g2, unique_profiles_g2 = experiment_G2(N=N2)
    data_g3, profile_groups_g3, multi_groups_g3 = experiment_G3(N=N3)

    total = time.time() - t_start

    print(f"\n{'='*70}")
    print(f"TOTAL TIME: {total:.1f}s")
    print(f"{'='*70}")

    print(f"""
{'='*70}
SYNTHESIS: Carry Geometry Stage 1.2
{'='*70}

THE KEY QUESTION: Does carry-geometry predict output geometry?

If corr(Hamming_carry, Hamming_H7) > 0.05 → carry algebra has
a "window" through the chaotic zone (rounds 31-59).

If corr ≈ 0 → carry structure is internal only, invisible from output.
The chaotic zone is truly opaque.

Either result is valuable:
  Signal → path to new mathematics
  No signal → carry-algebra alone is insufficient,
              need to combine with schedule dual algebra
""")
