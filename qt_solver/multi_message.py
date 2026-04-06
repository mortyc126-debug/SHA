"""
Direction 2: Multi-message OMEGA — birthday search among R=15 solutions.

Concept:
  At R=15, α-kernel = 32. One message M generates 2^32 messages with same R=15 hash.
  Extending to R=16: these 2^32 messages produce DIFFERENT R=16 hashes.
  Question: how much variation? Is there structure to exploit?

  If R=16 hash variation spans only ~32 bits (one round's contribution),
  birthday search in 32-bit space needs ~2^16 attempts → O(2^16) collision at R=16.

  If variation is full 256 bits (random) → birthday needs 2^128 → no advantage.

Analysis steps:
  1. Generate α-kernel at R=15
  2. Sample messages from kernel (enumerate or random linear combos)
  3. Compute their R=16 hashes
  4. Measure hash variation: how many bits actually change?
  5. If variation is restricted → birthday attack is feasible
"""

import random
from qt_solver.sha256_traced import MASK32, get_bit, sha256_compress
from qt_solver.omega_solve import omega_linearize_and_solve
from qt_solver.gf2 import bitvec_weight


def analyze_r16_variation(msg=None, seed=42, n_samples=256, verbose=True):
    """
    Generate R=15 preimage family, extend to R=16, measure variation.
    """
    if msg is None:
        rng = random.Random(seed)
        msg = [rng.randint(0, MASK32) for _ in range(16)]

    if verbose:
        print("=" * 60)
        print("Multi-message OMEGA: R=15 → R=16 variation analysis")
        print("=" * 60)

    # Get α-kernel at R=15
    res15 = omega_linearize_and_solve(15, msg, verbose=False)
    kernel = res15['kernel']
    alpha_dim = len([k for k in kernel
                     if sha256_compress(_apply_delta(msg, k), 15) == res15['hash']])

    if verbose:
        print(f"  R=15 α-kernel: {res15['alpha_dim']} vectors")
        print(f"  Verified preimage vectors: {alpha_dim}")

    # Filter to actual preimage vectors
    good_kernel = [k for k in kernel
                   if sha256_compress(_apply_delta(msg, k), 15) == res15['hash']]

    if not good_kernel:
        print("  ERROR: no valid preimage vectors found")
        return None

    # Verify R=15 hash is same for all
    ref_hash_15 = sha256_compress(msg, 15)

    # Sample messages from kernel and compute R=16 hashes
    rng = random.Random(seed + 1)
    r16_hashes = []
    r16_hashes_set = set()

    # Base message R=16 hash
    base_hash_16 = sha256_compress(msg, 16)
    r16_hashes.append(tuple(base_hash_16))
    r16_hashes_set.add(tuple(base_hash_16))

    if verbose:
        print(f"\n  Base R=16 hash: {[hex(w) for w in base_hash_16]}")
        print(f"  Sampling {n_samples} kernel combinations...")

    collisions_16 = 0

    for sample in range(n_samples):
        # Random linear combination of kernel vectors
        delta = 0
        for kv in good_kernel:
            if rng.randint(0, 1):
                delta ^= kv

        if delta == 0:
            continue

        msg2 = _apply_delta(msg, delta)

        # Verify R=15 still matches
        h15 = sha256_compress(msg2, 15)
        if h15 != ref_hash_15:
            continue  # not a valid preimage (linearization error)

        # Compute R=16 hash
        h16 = sha256_compress(msg2, 16)
        h16_tuple = tuple(h16)
        r16_hashes.append(h16_tuple)

        if h16_tuple in r16_hashes_set:
            collisions_16 += 1
            if verbose:
                print(f"  ★ R=16 COLLISION found at sample {sample}!")
        r16_hashes_set.add(h16_tuple)

    if verbose:
        print(f"\n  Results:")
        print(f"  Valid R=15 preimage samples: {len(r16_hashes)}")
        print(f"  Distinct R=16 hashes: {len(r16_hashes_set)}")
        print(f"  R=16 collisions: {collisions_16}")

    # Analyze hash variation
    if len(r16_hashes) >= 2:
        # Compute pairwise Hamming distances
        dists = []
        for i in range(min(100, len(r16_hashes))):
            for j in range(i+1, min(100, len(r16_hashes))):
                d = sum(bin(r16_hashes[i][w] ^ r16_hashes[j][w]).count('1') for w in range(8))
                dists.append(d)

        avg_dist = sum(dists) / len(dists) if dists else 0
        min_dist = min(dists) if dists else 0
        max_dist = max(dists) if dists else 0

        if verbose:
            print(f"\n  R=16 hash distances (pairwise Hamming):")
            print(f"    Mean: {avg_dist:.1f} / 256")
            print(f"    Min:  {min_dist} / 256")
            print(f"    Max:  {max_dist} / 256")
            print(f"    (128.0 = random, <128 = structure)")

        # Per-word analysis
        print(f"\n  Per-word hash variation:")
        for w in range(8):
            word_vals = set()
            for h in r16_hashes:
                word_vals.add(h[w])
            # Count distinct bits that change
            all_xor = 0
            for h in r16_hashes:
                all_xor |= (h[w] ^ r16_hashes[0][w])
            changing_bits = bin(all_xor).count('1')
            print(f"    H[{w}]: {len(word_vals):4d} distinct values, {changing_bits}/32 bits vary")

        return {
            'n_samples': len(r16_hashes),
            'n_distinct': len(r16_hashes_set),
            'collisions': collisions_16,
            'avg_dist': avg_dist,
            'min_dist': min_dist,
            'max_dist': max_dist,
        }

    return None


def birthday_search_r16(msg=None, seed=42, max_attempts=10000, verbose=True):
    """
    Birthday search for R=16 collision among R=15 preimage family.

    Strategy:
      1. Get R=15 kernel (32 dimensions → 2^32 messages)
      2. Sample random messages from family
      3. Hash each at R=16
      4. Look for collision (birthday paradox)

    If R=16 hashes are truly random 256-bit: need ~2^128 → hopeless
    If R=16 hashes constrained to k-bit subspace: need ~2^(k/2) → possible if k ≤ 64
    """
    if msg is None:
        rng = random.Random(seed)
        msg = [rng.randint(0, MASK32) for _ in range(16)]

    res15 = omega_linearize_and_solve(15, msg, verbose=False)
    kernel = res15['kernel']
    good_kernel = [k for k in kernel
                   if sha256_compress(_apply_delta(msg, k), 15) == res15['hash']]

    if verbose:
        print(f"\nBirthday search R=16")
        print(f"  R=15 preimage kernel: {len(good_kernel)} vectors")
        print(f"  Family size: 2^{len(good_kernel)}")
        print(f"  Max attempts: {max_attempts}")

    rng = random.Random(seed + 999)
    seen = {}  # hash_tuple -> msg
    collisions = []

    for attempt in range(max_attempts):
        # Random kernel combination
        delta = 0
        for kv in good_kernel:
            if rng.randint(0, 1):
                delta ^= kv

        msg2 = _apply_delta(msg, delta)

        # Verify R=15
        h15 = sha256_compress(msg2, 15)
        if h15 != sha256_compress(msg, 15):
            continue

        h16 = tuple(sha256_compress(msg2, 16))

        if h16 in seen:
            msg_prev = seen[h16]
            if msg_prev != msg2:
                collisions.append((msg_prev, msg2))
                if verbose:
                    print(f"  ★ COLLISION at attempt {attempt}!")
                    print(f"    M1: {[hex(w) for w in msg_prev[:4]]}...")
                    print(f"    M2: {[hex(w) for w in msg2[:4]]}...")
                    d = sum(1 for i in range(16) if msg_prev[i] != msg2[i])
                    print(f"    Differ in {d} words")
        seen[h16] = msg2

    if verbose:
        print(f"\n  Attempted: {max_attempts}")
        print(f"  Distinct hashes: {len(seen)}")
        print(f"  Collisions found: {len(collisions)}")
        if not collisions:
            print(f"  (Expected: 0 — hash variation is 256-bit random)")

    return collisions


def _apply_delta(msg, delta):
    """Apply bit-vector delta to message words."""
    msg2 = list(msg)
    for wi in range(16):
        word_bits = 0
        for bi in range(32):
            if (delta >> (wi * 32 + bi)) & 1:
                word_bits |= (1 << bi)
        msg2[wi] ^= word_bits
    return msg2


if __name__ == '__main__':
    analyze_r16_variation(n_samples=500, verbose=True)
    birthday_search_r16(max_attempts=5000, verbose=True)
