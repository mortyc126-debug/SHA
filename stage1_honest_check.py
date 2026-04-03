"""
Stage 1, Part 7: Verify clustering + honest check.

Two things to settle:
1. Is clustering real or artifact of linear grid sampling?
2. Does it survive at 64 rounds?

Then: regardless of result, THINK about what to invent next.
"""

import random
from stage0_observe import (
    sha256_round_trace, MASK, H0, K, schedule,
    rotr, Sig0, Sig1, Ch, Maj
)


def sha256_reduced(msg16, rounds):
    W = schedule(msg16)
    a, b, c, d, e, f, g, h = H0
    for r in range(rounds):
        T1 = (h + Sig1(e) + Ch(e, f, g) + K[r] + W[r]) & MASK
        T2 = (Sig0(a) + Maj(a, b, c)) & MASK
        h, g, f = g, f, e
        e = (d + T1) & MASK
        d, c, b = c, b, a
        a = (T1 + T2) & MASK
    return [(a + H0[0]) & MASK, (b + H0[1]) & MASK,
            (c + H0[2]) & MASK, (d + H0[3]) & MASK,
            (e + H0[4]) & MASK, (f + H0[5]) & MASK,
            (g + H0[6]) & MASK, (h + H0[7]) & MASK]


def honest_clustering_check():
    """
    Check clustering with RANDOM sampling (not grid).
    Compare L1 distances within same-hash groups vs random groups.
    Test at R=4, 16, 64.
    """
    print("=" * 80)
    print("HONEST CHECK: Is preimage clustering real or sampling artifact?")
    print("=" * 80)

    for R in [4, 16, 64]:
        N = 50000  # Many random messages

        rng = random.Random(42)
        base_msg = [rng.randint(0, MASK) for _ in range(16)]

        # Only vary W[0] (32-bit), keep rest fixed
        # Hash down to 8 bits (low byte of register a)
        hash_to_w0 = {}

        for trial in range(N):
            w0 = random.Random(trial).randint(0, MASK)
            msg = list(base_msg)
            msg[0] = w0
            h = sha256_reduced(msg, R)
            key = h[0] & 0xFF

            if key not in hash_to_w0:
                hash_to_w0[key] = []
            hash_to_w0[key].append(w0)

        # Pick groups with >= 10 members
        groups = [v for v in hash_to_w0.values() if len(v) >= 10]

        if not groups:
            print(f"\n  R={R}: Not enough groups with >=10 members")
            continue

        # For each group: compute average |w0_i - w0_j| mod 2^32 (additive distance)
        # Compare with random pairs
        within_dists = []
        for group in groups[:50]:
            for i in range(min(20, len(group))):
                for j in range(i+1, min(20, len(group))):
                    d = abs(group[i] - group[j])
                    if d > (1 << 31):
                        d = (1 << 32) - d
                    within_dists.append(d)

        # Random pairs (any two random W0 values)
        random_dists = []
        all_w0s = []
        for group in groups:
            all_w0s.extend(group[:5])
        for i in range(min(5000, len(all_w0s))):
            for j in range(i+1, min(i+5, len(all_w0s))):
                d = abs(all_w0s[i] - all_w0s[j])
                if d > (1 << 31):
                    d = (1 << 32) - d
                random_dists.append(d)

        avg_within = sum(within_dists) / len(within_dists) if within_dists else 0
        avg_random = sum(random_dists) / len(random_dists) if random_dists else 0

        # Also: XOR distance
        within_xor = []
        for group in groups[:50]:
            for i in range(min(20, len(group))):
                for j in range(i+1, min(20, len(group))):
                    within_xor.append(bin(group[i] ^ group[j]).count('1'))

        random_xor = []
        for i in range(min(5000, len(all_w0s))):
            for j in range(i+1, min(i+5, len(all_w0s))):
                random_xor.append(bin(all_w0s[i] ^ all_w0s[j]).count('1'))

        avg_within_xor = sum(within_xor) / len(within_xor) if within_xor else 0
        avg_random_xor = sum(random_xor) / len(random_xor) if random_xor else 0

        ratio_add = avg_within / avg_random if avg_random > 0 else 0
        ratio_xor = avg_within_xor / avg_random_xor if avg_random_xor > 0 else 0

        print(f"\n  R={R} rounds, {N} random messages, varying W[0] only:")
        print(f"    Groups with ≥10 members: {len(groups)}")
        print(f"    Additive distance: within={avg_within:.0f} random={avg_random:.0f} ratio={ratio_add:.3f}")
        print(f"    XOR distance:      within={avg_within_xor:.1f}  random={avg_random_xor:.1f}  ratio={ratio_xor:.3f}")
        print(f"    → {'CLUSTERED (real signal!)' if ratio_add < 0.9 else 'NOT CLUSTERED (artifact)'}")


if __name__ == "__main__":
    honest_clustering_check()
