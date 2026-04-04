"""
BTE-8 ATTACK: Can our theory beat birthday on a small BTE?

BTE-8: n=8, R=16, hash=64 bits. Birthday collision = 2^32.

Our theory says:
  - Layer rank = 2R-1 = 31 (per trajectory layer)
  - Hash rank per layer = 8 (limited by output)
  - 4 layers needed for full coverage
  - Each layer = degree-2 GF(2) system
  - System is SPARSE (each eq involves ~13 vars)

Can we use layer structure to find collision CHEAPER than 2^32?

Approach: Generate many messages. For each, compute layer-0 hash
(bit-0 of each register). Group by layer-0 hash (2^8 = 256 groups).
Within each group: messages agree on layer 0.
Check layer 1: how many also agree? Layer 2? Full hash?

If layer agreement PROPAGATES → we filter efficiently → fewer full-hash checks.
We showed MI=0 for SHA-256... but BTE-8 might be different (weaker).
"""

import random
import time


MASK8 = 0xFF


def rotr8(x, k):
    return ((x >> k) | (x << (8 - k))) & MASK8


def bte8_hash(msg, R=16):
    """Compute BTE-8 hash."""
    n = 8
    mask = MASK8
    n_msg = 16

    IV = [0x6a, 0xbb, 0x3c, 0xa5, 0x51, 0x9b, 0x1f, 0x5b]
    K = [(i * 37 + 113) & mask for i in range(R)]

    W = list(msg[:n_msg])
    for i in range(n_msg, R):
        W.append((W[i-2] ^ W[i-7] ^ W[i-n_msg]) & mask)

    def sig0(x): return rotr8(x, 1) ^ rotr8(x, 3) ^ rotr8(x, 5)
    def sig1(x): return rotr8(x, 2) ^ rotr8(x, 4) ^ rotr8(x, 6)
    def ch(e, f, g): return (e & f) ^ (~e & g) & mask
    def maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)

    a, b, c, d, e, f, g, h = IV
    for r in range(R):
        T1 = (h + sig1(e) + ch(e, f, g) + K[r] + W[r]) & mask
        T2 = (sig0(a) + maj(a, b, c)) & mask
        h, g, f = g, f, e
        e = (d + T1) & mask
        d, c, b = c, b, a
        a = (T1 + T2) & mask

    return (a, b, c, d, e, f, g, h)


def experiment_birthday_baseline():
    """Standard birthday attack on BTE-8. Measure actual cost."""
    print("=" * 80)
    print("BASELINE: Birthday collision on BTE-8 (hash=64 bits)")
    print("=" * 80)

    t0 = time.time()
    seen = {}
    attempts = 0

    for seed in range(10_000_000):
        M = [(seed * 37 + i * 71 + seed // 256 * 13) & MASK8 for i in range(16)]
        H = bte8_hash(M)
        attempts += 1

        if H in seen:
            t1 = time.time()
            M_other = seen[H]
            if M != M_other:
                print(f"  COLLISION found at attempt {attempts} ({t1-t0:.2f}s)")
                print(f"    M1 = {M[:4]}...")
                print(f"    M2 = {M_other[:4]}...")
                print(f"    H  = {H}")
                print(f"  Expected birthday: 2^32 ≈ {2**32}")
                print(f"  Actual: {attempts} ≈ 2^{attempts.bit_length()-1}")
                return attempts
        else:
            seen[H] = M

        if attempts % 1_000_000 == 0:
            print(f"    {attempts/1e6:.0f}M attempts ({time.time()-t0:.1f}s)...")

    print(f"  No collision in {attempts} attempts")
    return attempts


def experiment_layered_sieve():
    """
    Layered sieve: group by layer-0, then filter by layer-1, etc.

    Birthday expects 2^32 attempts for 64-bit hash.
    Layer-0 birthday expects 2^4 for 8-bit layer-0.
    If we can cheaply filter layer-0 matches to also match layer-1...

    Approach:
    1. Generate N messages, compute full hash.
    2. Group by layer-0 (bit-0 of each reg = 8-bit key). ~N/256 per group.
    3. Within each group: check layer-1 (bit-1 of each reg).
    4. Matches → check full hash.

    Cost: N hash computations + comparison cost.
    If layers independent: need N ≈ 2^32 anyway (same as birthday).
    If layers correlated for BTE-8 (weaker than SHA-256): need N < 2^32.
    """
    print("\n" + "=" * 80)
    print("LAYERED SIEVE: Can layer structure reduce collision cost?")
    print("=" * 80)

    N = 2_000_000  # Try 2M messages
    t0 = time.time()

    # Store messages by layer-0 hash
    layer0_groups = {}
    full_hashes = {}
    collision_found = False

    for seed in range(N):
        M = [(seed * 37 + i * 71 + seed // 256 * 13) & MASK8 for i in range(16)]
        H = bte8_hash(M)

        # Check full collision first
        if H in full_hashes:
            M_other = full_hashes[H]
            if M != M_other:
                t1 = time.time()
                print(f"  FULL COLLISION at attempt {seed+1} ({t1-t0:.2f}s)")
                collision_found = True
                break
        full_hashes[H] = M

        # Layer-0 key
        l0 = tuple((H[reg] & 1) for reg in range(8))
        if l0 not in layer0_groups:
            layer0_groups[l0] = []
        layer0_groups[l0].append((M, H))

    t1 = time.time()

    if not collision_found:
        print(f"  No full collision in {N} attempts ({t1-t0:.2f}s)")

    # Analyze layer-0 groups
    group_sizes = [len(v) for v in layer0_groups.values()]
    n_groups = len(layer0_groups)

    print(f"\n  Layer-0 groups: {n_groups}/256")
    print(f"  Avg group size: {sum(group_sizes)/n_groups:.0f} (expected {N//256})")
    print(f"  Max group size: {max(group_sizes)}")

    # Within each group: how many layer-1 matches?
    total_l0_pairs = 0
    total_l1_matches = 0
    total_full_matches = 0

    for l0_key, members in layer0_groups.items():
        if len(members) < 2:
            continue

        for i in range(min(50, len(members))):
            for j in range(i+1, min(50, len(members))):
                Mi, Hi = members[i]
                Mj, Hj = members[j]
                total_l0_pairs += 1

                # Layer-1 match?
                l1_match = all(((Hi[r] >> 1) & 1) == ((Hj[r] >> 1) & 1) for r in range(8))
                if l1_match:
                    total_l1_matches += 1

                    # Full match?
                    if Hi == Hj:
                        total_full_matches += 1

    print(f"\n  Within layer-0 groups:")
    print(f"    L0 pairs examined: {total_l0_pairs}")
    print(f"    L1 also matches: {total_l1_matches}")
    print(f"    Full hash matches: {total_full_matches}")

    if total_l0_pairs > 0:
        p_l1 = total_l1_matches / total_l0_pairs
        expected_random = 1 / 256  # 2^-8 for 8-bit layer
        ratio = p_l1 / expected_random if expected_random > 0 else 0
        print(f"    P(L1 match | L0 match) = {p_l1:.6f} (random = {expected_random:.6f})")
        print(f"    Ratio: {ratio:.3f}×")

        if ratio > 1.5:
            print(f"    *** LAYERS CORRELATED IN BTE-8! ***")
            saving = ratio
            effective_cost = N / saving
            print(f"    Effective cost: {N} / {saving:.1f} ≈ {effective_cost:.0f} ≈ 2^{effective_cost.bit_length()}")
        else:
            print(f"    Layers independent in BTE-8 too.")


def experiment_algebraic_collision():
    """
    ALGEBRAIC approach: use the degree-2 structure.

    For BTE-8 collision: find M1 ≠ M2 with H(M1) = H(M2).
    Equivalently: find δM ≠ 0 with H(M) XOR H(M XOR δM) = 0.

    This is 64 equations in 128 bits (δM) with degree... what?
    If we fix M and vary δM: δH(δM) = H(M) XOR H(M XOR δM).
    δH is a function of 128 bits (δM). Its degree in δM?

    At bit-0 level: we showed degree ≤ 2 (from Ch/Maj).
    But δH includes CARRY effects which raise degree.

    For BTE-8: effective degree through our layers = degree 2^D where D ≈ 3-4.
    So degree ≈ 2^3 = 8 or 2^4 = 16.

    XL on degree-8 system in 128 variables:
    C(128, 8) ≈ 2^34. Hmm, worse than birthday (2^32).

    But: the system is SPARSE. And our layered structure might help.
    """
    print("\n" + "=" * 80)
    print("ALGEBRAIC: Degree of δH as function of δM for BTE-8")
    print("=" * 80)

    rng = random.Random(42)
    M = [rng.randint(0, MASK8) for _ in range(16)]
    H_base = bte8_hash(M)

    # Measure degree of δH[0] bit 0 as function of δM
    # Test degree 1, 2, 3
    N = 300
    deg1_violations = 0
    deg2_violations = 0
    total = 0

    for trial in range(N):
        rng2 = random.Random(trial * 1000)
        w1, b1 = rng2.randint(0, 15), rng2.randint(0, 7)
        w2, b2 = rng2.randint(0, 15), rng2.randint(0, 7)
        w3, b3 = rng2.randint(0, 15), rng2.randint(0, 7)
        if (w1, b1) == (w2, b2) or (w1, b1) == (w3, b3) or (w2, b2) == (w3, b3):
            continue

        def dH0b0(flips):
            M2 = list(M)
            for w, b in flips:
                M2[w] ^= (1 << b)
            H2 = bte8_hash(M2)
            return (H_base[0] ^ H2[0]) & 1

        f0 = 0  # dH for no flip
        f1 = dH0b0([(w1, b1)])
        f2 = dH0b0([(w2, b2)])
        f12 = dH0b0([(w1, b1), (w2, b2)])

        d1 = f0 ^ f1 ^ f2 ^ f12
        if d1:
            deg1_violations += 1

        f3 = dH0b0([(w3, b3)])
        f13 = dH0b0([(w1, b1), (w3, b3)])
        f23 = dH0b0([(w2, b2), (w3, b3)])
        f123 = dH0b0([(w1, b1), (w2, b2), (w3, b3)])

        d2 = f0 ^ f1 ^ f2 ^ f3 ^ f12 ^ f13 ^ f23 ^ f123
        if d2:
            deg2_violations += 1

        total += 1

    print(f"  δH[0] bit 0 as function of δM (N={total}):")
    print(f"    Degree-1 violations: {deg1_violations}/{total} = {deg1_violations/total:.3f}")
    print(f"    Degree-2 violations: {deg2_violations}/{total} = {deg2_violations/total:.3f}")

    if deg1_violations == 0:
        print(f"    → AFFINE (degree 1)")
    elif deg2_violations == 0:
        print(f"    → QUADRATIC (degree 2)")
    elif deg2_violations < total * 0.1:
        print(f"    → Nearly degree 2 ({deg2_violations/total:.1%} violations)")
    else:
        print(f"    → HIGH DEGREE (≥3)")


if __name__ == "__main__":
    experiment_algebraic_collision()
    experiment_layered_sieve()
    experiment_birthday_baseline()
