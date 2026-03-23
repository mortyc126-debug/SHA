#!/usr/bin/env python3
"""
Step 26a: FUNCTIONAL GRAPH OF SHA-256
═══════════════════════════════════════

Random mapping f: [N] → [N] has known statistics (Flajolet-Odlyzko 1990):
  - Expected tail length (rho):  √(πN/8)
  - Expected cycle length (lambda): √(πN/8)
  - Expected rho+lambda:           √(πN/2)
  - Number of components:          (1/2)·ln(N)
  - Number of cyclic nodes:        √(πN/2)
  - Largest component:             0.7582·N
  - Largest tree:                  0.4834·N

QUESTION: Does mini-SHA match random mapping statistics?
If YES → no structural shortcut from graph theory.
If NO  → deviation reveals exploitable structure.
"""

import struct, hashlib, time, math, random
from collections import Counter

MASK32 = 0xFFFFFFFF

# ─── Mini SHA-256 compression (reduced to n-bit output) ─────────────

def rotr(x, n, bits=32):
    return ((x >> n) | (x << (bits - n))) & ((1 << bits) - 1)

def sha256_compress_rounds(msg_words, rounds=64):
    """Full SHA-256 compression with configurable rounds."""
    K = [
        0x428a2f98, 0x71374491, 0xb5c0f6e2, 0xe9b5dba5,
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
    IV = [
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
        0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
    ]

    W = list(msg_words) + [0] * (64 - len(msg_words))
    for i in range(16, 64):
        s0 = rotr(W[i-15], 7) ^ rotr(W[i-15], 18) ^ (W[i-15] >> 3)
        s1 = rotr(W[i-2], 17) ^ rotr(W[i-2], 19) ^ (W[i-2] >> 10)
        W[i] = (W[i-16] + s0 + W[i-7] + s1) & MASK32

    a, b, c, d, e, f, g, h = IV
    for i in range(min(rounds, 64)):
        S1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25)
        ch = (e & f) ^ (~e & g) & MASK32
        temp1 = (h + S1 + ch + K[i] + W[i]) & MASK32
        S0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22)
        maj = (a & b) ^ (a & c) ^ (b & c)
        temp2 = (S0 + maj) & MASK32

        h = g; g = f; f = e
        e = (d + temp1) & MASK32
        d = c; c = b; b = a
        a = (temp1 + temp2) & MASK32

    return [(IV[i] + v) & MASK32 for i, v in enumerate([a, b, c, d, e, f, g, h])]


def mini_sha_func(x, bits=16, rounds=64):
    """Map n-bit → n-bit using SHA-256 compression truncated to n bits."""
    # Pack x into message words
    W = [x & MASK32] + [0] * 15
    state = sha256_compress_rounds(W, rounds)
    # Truncate to n bits
    return state[0] & ((1 << bits) - 1)


# ─── Experiment 1: Functional Graph Statistics ───────────────────────

def measure_rho_lambda(f, N, samples=1000):
    """Measure tail (rho) and cycle (lambda) lengths via Floyd's algorithm."""
    rhos = []
    lambdas = []

    for _ in range(samples):
        x0 = random.randint(0, N - 1)

        # Floyd's cycle detection
        tortoise = f(x0)
        hare = f(f(x0))
        while tortoise != hare:
            tortoise = f(tortoise)
            hare = f(f(hare))

        # Find mu (tail length)
        mu = 0
        tortoise = x0
        while tortoise != hare:
            tortoise = f(tortoise)
            hare = f(hare)
            mu += 1

        # Find lambda (cycle length)
        lam = 1
        hare = f(tortoise)
        while tortoise != hare:
            hare = f(hare)
            lam += 1

        rhos.append(mu)
        lambdas.append(lam)

    return rhos, lambdas


def full_graph_analysis(f, N):
    """Complete functional graph analysis for small N."""
    # Build the full graph
    succ = [f(x) for x in range(N)]

    # Find all cycles
    visited = [False] * N
    in_cycle = [False] * N

    for start in range(N):
        if visited[start]:
            continue
        path = []
        x = start
        while not visited[x]:
            visited[x] = True
            path.append(x)
            x = succ[x]
        # x is now either in a cycle or already processed
        if x in path:
            idx = path.index(x)
            for node in path[idx:]:
                in_cycle[node] = True

    num_cyclic = sum(in_cycle)

    # Find components
    # Build reverse graph
    rev = [[] for _ in range(N)]
    for x in range(N):
        rev[succ[x]].append(x)

    # Each cyclic node is root of a tree
    # Count components = number of cycles
    cycle_visited = [False] * N
    num_components = 0
    component_sizes = []

    for x in range(N):
        if in_cycle[x] and not cycle_visited[x]:
            num_components += 1
            # BFS to find full component
            comp_size = 0
            # First trace the cycle
            cycle_nodes = [x]
            y = succ[x]
            while y != x:
                cycle_nodes.append(y)
                y = succ[y]

            # BFS from all cycle nodes
            queue = list(cycle_nodes)
            comp_visited = set(queue)
            while queue:
                node = queue.pop(0)
                comp_size += 1
                for pred in rev[node]:
                    if pred not in comp_visited:
                        comp_visited.add(pred)
                        queue.append(pred)

            for node in comp_visited:
                if in_cycle[node]:
                    cycle_visited[node] = True

            component_sizes.append(comp_size)

    # Image size (number of distinct values in range)
    image = len(set(succ))

    # Tail lengths for all nodes
    tail_lengths = [0] * N
    for x in range(N):
        if in_cycle[x]:
            tail_lengths[x] = 0
        else:
            length = 0
            y = x
            while not in_cycle[y]:
                y = succ[y]
                length += 1
            tail_lengths[x] = length

    return {
        'num_cyclic': num_cyclic,
        'num_components': num_components,
        'component_sizes': sorted(component_sizes, reverse=True),
        'image_size': image,
        'max_tail': max(tail_lengths),
        'avg_tail': sum(tail_lengths) / N,
        'tail_distribution': tail_lengths
    }


def random_func_expected(N):
    """Flajolet-Odlyzko expected values for random mapping [N]→[N]."""
    sqrtN = math.sqrt(N)
    return {
        'rho': math.sqrt(math.pi * N / 8),
        'lambda': math.sqrt(math.pi * N / 8),
        'rho_plus_lambda': math.sqrt(math.pi * N / 2),
        'num_components': 0.5 * math.log(N),
        'num_cyclic': math.sqrt(math.pi * N / 2),
        'largest_component': 0.7582 * N,
        'image_size': (1 - 1/math.e) * N,  # ≈ 0.6321·N
    }


# ─── Experiment 2: Compare SHA vs truly random ───────────────────────

def run_comparison(bits_list=[12, 14, 16], rounds_list=[8, 16, 32, 64]):
    """Compare mini-SHA functional graph vs random mapping predictions."""

    print("FUNCTIONAL GRAPH: mini-SHA vs RANDOM MAPPING")
    print("=" * 70)
    print()
    print("Flajolet-Odlyzko (1990) predictions for random f: [N]→[N]:")
    print("  E[ρ] = E[λ] = √(πN/8)")
    print("  E[components] = (1/2)·ln(N)")
    print("  E[cyclic nodes] = √(πN/2)")
    print("  E[image size] = (1-1/e)·N ≈ 0.632·N")
    print()

    for bits in bits_list:
        N = 1 << bits
        expected = random_func_expected(N)

        print(f"{'='*70}")
        print(f"  N = 2^{bits} = {N}")
        print(f"{'='*70}")
        print()
        print(f"  {'Metric':<22} {'Random theory':>14} ", end="")
        for R in rounds_list:
            print(f"  {'R='+str(R):>8}", end="")
        print()
        print(f"  {'-'*22} {'-'*14} ", end="")
        for _ in rounds_list:
            print(f"  {'-'*8}", end="")
        print()

        results = {}
        for R in rounds_list:
            t0 = time.time()
            f = lambda x, r=R, b=bits: mini_sha_func(x, bits=b, rounds=r)
            stats = full_graph_analysis(f, N)
            results[R] = stats
            elapsed = time.time() - t0
            if bits >= 16:
                print(f"    [R={R}: {elapsed:.1f}s]", flush=True)

        # Print comparison
        metrics = [
            ('Cyclic nodes', 'num_cyclic', expected['num_cyclic']),
            ('Components', 'num_components', expected['num_components']),
            ('Largest comp', lambda s: s['component_sizes'][0] if s['component_sizes'] else 0, expected['largest_component']),
            ('Image size', 'image_size', expected['image_size']),
            ('Max tail', 'max_tail', expected['rho'] * 3),  # rough max ≈ 3·E[rho]
            ('Avg tail (non-cyc)', 'avg_tail', expected['rho']),
        ]

        for name, key, exp_val in metrics:
            print(f"  {name:<22} {exp_val:>14.1f} ", end="")
            for R in rounds_list:
                if callable(key):
                    val = key(results[R])
                else:
                    val = results[R][key]
                # Color: deviation from expected
                ratio = val / exp_val if exp_val > 0 else float('inf')
                flag = " " if 0.5 < ratio < 2.0 else "*"
                print(f"  {val:>7.1f}{flag}", end="")
            print()
        print()


# ─── Experiment 3: Birthday collision in mini-SHA ────────────────────

def birthday_collision_test(bits=16, rounds=64, trials=100):
    """Measure actual collision probability vs birthday bound."""
    N = 1 << bits
    expected_birthday = math.sqrt(math.pi * N / 2)  # ≈ 1.177·√N

    print(f"BIRTHDAY COLLISION TEST (n={bits} bits, R={rounds})")
    print(f"  Expected steps to collision (birthday): {expected_birthday:.1f}")
    print(f"  (= {math.log2(expected_birthday):.1f} bits of work)")
    print()

    steps_list = []
    for trial in range(trials):
        seen = {}
        x = random.randint(0, N - 1)
        for step in range(N * 4):  # safety limit
            h = mini_sha_func(x, bits=bits, rounds=rounds)
            if h in seen and seen[h] != x:
                steps_list.append(step)
                break
            seen[h] = x
            x = random.randint(0, N - 1)

    if steps_list:
        avg = sum(steps_list) / len(steps_list)
        med = sorted(steps_list)[len(steps_list) // 2]
        print(f"  Avg steps to collision:    {avg:.1f}")
        print(f"  Median steps:              {med:.1f}")
        print(f"  Ratio (avg/expected):      {avg/expected_birthday:.3f}")
        print(f"  → {'MATCHES random' if 0.8 < avg/expected_birthday < 1.3 else 'DEVIATES from random'}")
    print()


# ─── Experiment 4: Collision structure — multi-collision ─────────────

def multicollision_test(bits=16, rounds=64):
    """
    Joux (2004): For random oracle, k-collision in k·√(N) work.
    For Merkle-Damgård: k-collision in k·2^(n/2) = trivial!

    But for COMPRESSION FUNCTION (fixed IV, single block):
    k-collision should cost ≈ N^((k-1)/k) for random function.
    """
    N = 1 << bits
    print(f"MULTI-COLLISION TEST (n={bits} bits, R={rounds})")
    print(f"  Building image histogram...")

    # Count how many preimages each output has
    image_count = Counter()
    for x in range(N):
        h = mini_sha_func(x, bits=bits, rounds=rounds)
        image_count[h] += 1

    # Distribution
    count_dist = Counter(image_count.values())

    print(f"  Image size: {len(image_count)}/{N} ({len(image_count)/N:.4f})")
    print(f"  Expected (random): {(1-1/math.e):.4f}")
    print()
    print(f"  {'k-collision':>12} {'Count':>8} {'Expected (Poisson)':>20}")
    print(f"  {'-'*12} {'-'*8} {'-'*20}")

    for k in range(0, 8):
        actual = count_dist.get(k, 0)
        # Poisson: P(X=k) = e^{-1}/k! for random mapping
        expected_poisson = N * math.exp(-1) / math.factorial(k)
        print(f"  {k:>12} {actual:>8} {expected_poisson:>20.1f}")

    max_coll = max(image_count.values())
    print(f"\n  Max collision size: {max_coll}")
    print(f"  Expected max (random): ~ln(N)/ln(ln(N)) ≈ {math.log(N)/math.log(math.log(N)):.1f}")
    print()


# ─── Experiment 5: Transition zone — partial randomness ─────────────

def transition_zone_test(bits=16):
    """
    Measure HOW FAST mini-SHA approaches random mapping statistics.
    The transition zone (rounds 1-8) might have exploitable structure.
    """
    N = 1 << bits
    expected = random_func_expected(N)

    print(f"TRANSITION ZONE: approach to randomness (n={bits})")
    print(f"{'='*70}")
    print()
    print(f"  {'Rounds':>6}  {'Image%':>8}  {'Cyclic':>8}  {'Comps':>8}  "
          f"{'Largest%':>10}  {'Score':>8}")
    print(f"  {'-'*6}  {'-'*8}  {'-'*8}  {'-'*8}  {'-'*10}  {'-'*8}")

    for R in [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 16, 24, 32, 64]:
        f = lambda x, r=R, b=bits: mini_sha_func(x, bits=b, rounds=r)
        stats = full_graph_analysis(f, N)

        image_pct = stats['image_size'] / N
        cyclic = stats['num_cyclic']
        comps = stats['num_components']
        largest_pct = stats['component_sizes'][0] / N if stats['component_sizes'] else 0

        # "Randomness score": how close to random (0=structured, 1=random)
        scores = []
        if expected['image_size'] > 0:
            scores.append(1 - abs(image_pct - 0.6321) / 0.6321)
        if expected['num_cyclic'] > 0:
            scores.append(1 - min(abs(cyclic - expected['num_cyclic']) / expected['num_cyclic'], 1))
        if expected['num_components'] > 0:
            scores.append(1 - min(abs(comps - expected['num_components']) / expected['num_components'], 1))

        score = sum(scores) / len(scores) if scores else 0

        flag = "  ← STRUCTURED" if score < 0.7 else ""
        print(f"  {R:>6}  {image_pct:>8.4f}  {cyclic:>8}  {comps:>8}  "
              f"{largest_pct:>10.4f}  {score:>8.3f}{flag}")

    print()
    print(f"  Random mapping:  image=0.6321  cyclic={expected['num_cyclic']:.0f}  "
          f"comps={expected['num_components']:.1f}  largest=0.758")
    print()


# ─── Experiment 6: Can we exploit the rho structure? ─────────────────

def pollard_rho_collision(bits=16, rounds=64, trials=20):
    """
    Pollard's rho: find collision using O(1) memory, O(√N) time.
    This is the BEST generic attack on random functions.

    For SHA-256 (256-bit): O(2^128) time, O(1) memory.
    Question: does mini-SHA give collisions faster than √N?
    """
    N = 1 << bits
    expected = math.sqrt(math.pi * N / 2)

    print(f"POLLARD'S RHO COLLISION ({bits}-bit, R={rounds})")
    print(f"  Expected: √(πN/2) = {expected:.1f} steps")
    print()

    steps_list = []
    for trial in range(trials):
        x0 = random.randint(0, N - 1)

        # Floyd's: find collision point
        tortoise = mini_sha_func(x0, bits, rounds)
        hare = mini_sha_func(mini_sha_func(x0, bits, rounds), bits, rounds)
        steps = 1

        while tortoise != hare:
            tortoise = mini_sha_func(tortoise, bits, rounds)
            hare = mini_sha_func(mini_sha_func(hare, bits, rounds), bits, rounds)
            steps += 1

        # Find the actual collision
        tortoise = x0
        while mini_sha_func(tortoise, bits, rounds) != mini_sha_func(hare, bits, rounds):
            tortoise = mini_sha_func(tortoise, bits, rounds)
            hare = mini_sha_func(hare, bits, rounds)
            steps += 1

        if tortoise != hare:
            steps_list.append(steps)

    if steps_list:
        avg = sum(steps_list) / len(steps_list)
        print(f"  Avg steps: {avg:.1f}  (expected {expected:.1f})")
        print(f"  Ratio:     {avg/expected:.3f}")
        print(f"  → {'MATCHES random' if 0.7 < avg/expected < 1.5 else 'DEVIATION!'}")
    else:
        print(f"  No collisions found (all fixed points?)")
    print()


# ─── MAIN ────────────────────────────────────────────────────────────

if __name__ == '__main__':
    print("STEP 26a: FUNCTIONAL GRAPH OF SHA-256")
    print("Does SHA-256 match random mapping statistics?")
    print("=" * 70)
    print()

    t_start = time.time()

    # Exp 1: Full graph analysis (small N)
    print("EXPERIMENT 1: Full graph analysis")
    print("-" * 40)
    run_comparison(bits_list=[12, 14], rounds_list=[4, 8, 16, 64])

    # Exp 2: Birthday collision
    print("EXPERIMENT 2: Birthday collision timing")
    print("-" * 40)
    birthday_collision_test(bits=14, rounds=64, trials=200)
    birthday_collision_test(bits=14, rounds=8, trials=200)
    birthday_collision_test(bits=14, rounds=4, trials=200)

    # Exp 3: Multi-collision distribution
    print("EXPERIMENT 3: Multi-collision (Poisson test)")
    print("-" * 40)
    multicollision_test(bits=14, rounds=64)
    multicollision_test(bits=14, rounds=4)

    # Exp 4: Transition zone
    print("EXPERIMENT 4: Transition zone")
    print("-" * 40)
    transition_zone_test(bits=14)

    # Exp 5: Pollard's rho
    print("EXPERIMENT 5: Pollard's rho collision")
    print("-" * 40)
    pollard_rho_collision(bits=14, rounds=64, trials=50)
    pollard_rho_collision(bits=14, rounds=4, trials=50)

    elapsed = time.time() - t_start

    print("=" * 70)
    print(f"TOTAL TIME: {elapsed:.1f}s")
    print()
    print("IMPLICATIONS FOR FULL SHA-256 (256-bit):")
    print("  If SHA matches random mapping → birthday bound 2^128 is TIGHT")
    print("  Pollard rho: 2^128 time, O(1) memory")
    print("  van Oorschot-Wiener: 2^128/M time with M processors")
    print("  Wagner k-tree: needs INDEPENDENT lists (SHA gives only 1)")
    print("  → No generic sub-birthday attack possible")
