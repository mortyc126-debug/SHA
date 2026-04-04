"""
WEAK BTE: Construct a BTE instance where collision is CHEAPER than birthday.

SHA-256 has:
  - n=32, R=64, rotations {2,6,11,13,22,25}
  - 4 pure layers, birthday = 2^128

Can we make a BTE with FEWER pure layers?
If pure_layers = 2 → effective hash complexity = 2 × (2R-1) ≈ 4R
For n=32, R=64: 4×64 = 256 bits effective... same as hash.

Or: can we make a BTE where layers are NOT independent?
If layers CORRELATE → birthday in smaller space → cheaper.

Strategy: choose rotation constants that MINIMIZE mixing.
  - SHA-256: coverage in 3 steps (optimal)
  - Bad choice: coverage in many steps (poor mixing)
  - Worst: no coverage (some bits never reached)

If rotation constants don't generate Z/n → some bits are NEVER
connected to others → hash has STRUCTURAL weakness.
"""

import random
import time


def rotr_n(x, k, n):
    mask = (1 << n) - 1
    return ((x >> k) | (x << (n - k))) & mask


def ch_fn(e, f, g, mask):
    return (e & f) ^ (~e & g) & mask

def maj_fn(a, b, c, mask):
    return (a & b) ^ (a & c) ^ (b & c)


def bte_hash(msg, n, R, sig0_rots, sig1_rots):
    """Generic BTE hash."""
    mask = (1 << n) - 1
    n_msg = 16

    IV = [random.Random(42 + i).randint(0, mask) for i in range(8)]
    K = [random.Random(99999 + i).randint(0, mask) for i in range(R)]

    W = list(msg[:n_msg])
    for i in range(n_msg, R):
        W.append((W[i-2] ^ W[i-7 if i >= 7 else 0] ^ W[i-n_msg if i >= n_msg else 0]) & mask)

    def sig0(x):
        r = 0
        for rot in sig0_rots:
            r ^= rotr_n(x, rot, n)
        return r

    def sig1(x):
        r = 0
        for rot in sig1_rots:
            r ^= rotr_n(x, rot, n)
        return r

    a, b, c, d, e, f, g, h = IV
    for r in range(R):
        T1 = (h + sig1(e) + ch_fn(e, f, g, mask) + K[r] + W[r]) & mask
        T2 = (sig0(a) + maj_fn(a, b, c, mask)) & mask
        h, g, f = g, f, e
        e = (d + T1) & mask
        d, c, b = c, b, a
        a = (T1 + T2) & mask

    return (a, b, c, d, e, f, g, h)


def coverage_steps(n, rots):
    """How many steps to reach all n bit positions from position 0."""
    reached = {0}
    for step in range(n):
        new = set()
        for k in reached:
            for r in rots:
                new.add((k + r) % n)
        if new <= reached:
            return -1  # Never covers
        reached = reached | new
        if len(reached) == n:
            return step + 1
    return -1


def find_weak_rotations(n):
    """Find rotation sets that DON'T generate all of Z/n."""
    print(f"  Searching for non-generating rotations (n={n}):")

    weak = []
    for r1 in range(1, n):
        for r2 in range(r1+1, n):
            rots = [r1, r2]
            reached = {0}
            for _ in range(n):
                new = set()
                for k in reached:
                    for r in rots:
                        new.add((k + r) % n)
                if new <= reached:
                    break
                reached = reached | new
            if len(reached) < n:
                weak.append((r1, r2, len(reached)))

    if weak:
        for r1, r2, cover in weak[:10]:
            print(f"    Rots {{{r1},{r2}}}: covers {cover}/{n} positions")
    else:
        print(f"    No non-generating pairs found for n={n}")

    return weak


def experiment_coverage_vs_collision():
    """
    For n=8: test BTE instances with different coverage properties.
    Measure collision hardness (via hash entropy/distribution).
    """
    print("=" * 80)
    print("WEAK BTE: Coverage vs collision hardness")
    print("=" * 80)

    n = 8
    R = 16

    # Find weak rotations
    weak = find_weak_rotations(n)

    # Strong rotations (full coverage)
    strong = [(1, 3, 5, 2, 4, 6)]  # SHA-like, covers all 8

    # Test hash distribution for each
    configs = []

    # Add some weak configs
    for r1, r2, cover in weak[:5]:
        configs.append((f"Weak {{{r1},{r2}}} cover={cover}", [r1], [r2]))

    # Add strong
    configs.append(("Strong {1,3,5,2,4,6}", [1, 3, 5], [2, 4, 6]))
    configs.append(("Medium {1,4,2,5}", [1, 4], [2, 5]))

    N = 50000
    for name, s0, s1 in configs:
        # Generate N random hashes
        hash_set = set()
        for seed in range(N):
            rng = random.Random(seed)
            M = tuple(rng.randint(0, (1 << n) - 1) for _ in range(16))
            H = bte_hash(M, n, R, s0, s1)
            hash_set.add(H)

        unique = len(hash_set)
        collision_rate = 1 - unique / N  # Fraction of hashes that collided
        # Expected for random: 1 - (1 - 1/2^(n*8))^N ≈ N / 2^64 for n=8
        # ≈ 50000 / 2^64 ≈ 0

        # Better measure: number of collisions
        from collections import Counter
        hash_list = []
        for seed in range(N):
            rng = random.Random(seed)
            M = tuple(rng.randint(0, (1 << n) - 1) for _ in range(16))
            H = bte_hash(M, n, R, s0, s1)
            hash_list.append(H)

        counts = Counter(hash_list)
        n_collisions = sum(c - 1 for c in counts.values() if c > 1)
        max_collision = max(counts.values())

        # Hash entropy
        n_unique = len(counts)
        if n_unique > 1:
            probs = [c/N for c in counts.values()]
            entropy = -sum(p * (p if p > 0 else 1) for p in probs)  # rough
        else:
            entropy = 0

        print(f"\n  {name}:")
        print(f"    Unique hashes: {n_unique}/{N}")
        print(f"    Collisions: {n_collisions}")
        print(f"    Max same hash: {max_collision}")
        print(f"    Effective hash bits: ~{n_unique.bit_length() if n_unique > 0 else 0}")

        if n_collisions > N * 0.1:
            print(f"    *** WEAK: {n_collisions/N:.1%} collision rate! ***")
        elif n_collisions > 0:
            print(f"    Moderate: some collisions ({n_collisions})")
        else:
            print(f"    Strong: no collisions in {N} samples")


if __name__ == "__main__":
    experiment_coverage_vs_collision()
