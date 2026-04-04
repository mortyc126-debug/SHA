"""
COLLISION через BIT-LAYER DECOMPOSITION.

Collision: M ≠ M', SHA(M) = SHA(M').

В наших терминах: hash H = function of ALL 4 layers.
Два M дают одинаковый H ↔ все 256 бит hash совпадают.
256 бит hash = 32 бит × 8 регистров, каждый бит в СВОЁМ слое.

Hash bit k ∈ Layer(k). Для совпадения hash: ВСЕ 32 слоя совпадают.

Но: слои 5-31 ОПРЕДЕЛЕНЫ слоями 0-4 (через ротации).
Значит: коллизия ↔ слои 0-4 дают одинаковый hash.

Слой 0: 127 constraints. Коллизия в слое 0: два M дают одинаковый
         bit-0 hash (8 бит = bits 0 of 8 hash words).
Слой 1: +127 constraints. Adds bit-1 hash match (8 more bits).
...
Слой 4: +4 constraints. Completes.

ВОПРОС: Можно ли найти коллизию ПОСЛОЙНО?
  Step 0: Find M, M' with same bit-0 of hash (8-bit match).
          Cost: birthday in 2^8 → need ~2^4 pairs.
          But constrained by layer-0 structure (127 constraints).
  Step 1: Among those pairs, find one matching bit-1 too.
          ...

Let's compute: what's the birthday cost in each layer?
"""

import random
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, schedule
)


def sha256_hash(M):
    states, _ = sha256_round_trace(M)
    return [(states[64][i] + H0[i]) & MASK for i in range(8)]


def experiment_layered_collision():
    """
    For each bit position k, find pairs with matching hash at bit k.
    Measure: how many pairs needed? Does it factor across layers?
    """
    print("=" * 80)
    print("LAYERED COLLISION: Birthday cost per bit-layer")
    print("=" * 80)

    N = 50000  # messages to try

    # Compute hashes for N random messages
    hashes = []
    messages = []
    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        H = sha256_hash(M)
        hashes.append(H)
        messages.append(M)

    # For each "layer mask" (bits 0, bits 0-1, bits 0-2, ..., bits 0-7),
    # find collisions: pairs with identical masked hash.
    for n_bits in [1, 2, 3, 4, 8]:
        mask = (1 << n_bits) - 1  # e.g., 1, 3, 7, 15, 255

        # Hash projection: take bits 0..n_bits-1 of each of 8 hash words
        # This gives n_bits × 8 bits = n_bits × 8 bit hash projection.
        projected = {}
        collisions = 0
        first_collision = -1

        for i in range(N):
            key = tuple((hashes[i][reg] & mask) for reg in range(8))
            if key in projected:
                collisions += 1
                if first_collision == -1:
                    first_collision = i
            else:
                projected[key] = i

        proj_bits = n_bits * 8
        expected_birthday = 2 ** (proj_bits / 2)

        print(f"  Bits 0..{n_bits-1} ({proj_bits}-bit projection): "
              f"{collisions} collisions in {N}, "
              f"first at message #{first_collision}, "
              f"birthday expects 2^{proj_bits/2:.0f} = {expected_birthday:.0f}")

    # KEY: For 1-bit layer (8-bit projection): birthday = 2^4 = 16.
    # For 2-bit layer (16-bit): birthday = 2^8 = 256.
    # For 4-bit layer (32-bit): birthday = 2^16 = 65536.
    # For 8-bit layer (64-bit): birthday = 2^32.
    # For full 32-bit (256-bit): birthday = 2^128.

    # LAYERED approach: find collision in layer 0, then extend to layer 1, etc.
    # Layer 0: 8-bit projection → birthday 2^4. CHEAP.
    # Among layer-0 collisions: fraction matching layer 1 too?


def experiment_layer_conditional():
    """
    CRITICAL: Among pairs matching in layer 0 (bit-0 of all 8 regs),
    what fraction ALSO matches in layer 1 (bit-1)?

    If fraction = 2^{-8} (random): layers are independent for collision.
    Total cost = 2^4 × 2^8 × 2^8 × 2^8 × ... → birthday same as before.

    If fraction > 2^{-8}: layers are CORRELATED for collision → potential speedup.
    """
    print("\n" + "=" * 80)
    print("LAYER CONDITIONAL: Among bit-0 collisions, how many match bit-1?")
    print("=" * 80)

    N = 200000

    # Build hash table for bit-0 projection (8-bit)
    bit0_groups = {}  # key = (h[0]&1, h[1]&1, ..., h[7]&1)
    all_hashes = []

    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        H = sha256_hash(M)
        all_hashes.append(H)

        key0 = tuple((H[reg] & 1) for reg in range(8))
        if key0 not in bit0_groups:
            bit0_groups[key0] = []
        bit0_groups[key0].append(seed)

    # Find groups with many members (bit-0 collisions)
    large_groups = [(k, v) for k, v in bit0_groups.items() if len(v) >= 20]
    print(f"  Bit-0 groups with ≥20 members: {len(large_groups)}")
    print(f"  (Expected: 256 groups × {N}/256 = {N//256} per group)")

    if not large_groups:
        print("  Not enough for analysis")
        return

    # Among bit-0 collision pairs: what fraction also matches bit-1?
    bit1_matches = 0
    bit1_total = 0
    bit01_matches = 0

    for key0, members in large_groups[:50]:
        for i in range(min(30, len(members))):
            for j in range(i+1, min(30, len(members))):
                idx_i = members[i]
                idx_j = members[j]

                # These already match on bit 0 (by construction)
                # Do they also match on bit 1?
                match_bit1 = all((all_hashes[idx_i][reg] >> 1) & 1 ==
                                 (all_hashes[idx_j][reg] >> 1) & 1
                                 for reg in range(8))

                # Do they match on bits 0 AND 1 together?
                match_bit01 = all(all_hashes[idx_i][reg] & 3 ==
                                  all_hashes[idx_j][reg] & 3
                                  for reg in range(8))

                bit1_total += 1
                if match_bit1:
                    bit1_matches += 1
                if match_bit01:
                    bit01_matches += 1

    if bit1_total > 0:
        p_bit1 = bit1_matches / bit1_total
        p_bit01 = bit01_matches / bit1_total
        expected_random = 2 ** (-8)  # 8 independent bits

        print(f"\n  Among bit-0 collision pairs:")
        print(f"    P(also match bit-1): {p_bit1:.6f} (random = {expected_random:.6f})")
        print(f"    P(match bits 0+1):   {p_bit01:.6f} (random = {2**(-8):.6f})")
        print(f"    Ratio: {p_bit1/expected_random:.2f}×")

        if p_bit1 > expected_random * 1.5:
            print(f"    *** CORRELATED! Bit-0 match predicts bit-1 match! ***")
            print(f"    This would mean: layered collision is CHEAPER than birthday.")
        elif p_bit1 < expected_random * 0.7:
            print(f"    *** ANTI-CORRELATED! Bit-0 match makes bit-1 harder! ***")
        else:
            print(f"    Independent: bit-0 match gives no advantage for bit-1.")
            print(f"    Layered collision cost = product of per-layer costs = birthday.")


if __name__ == "__main__":
    experiment_layered_collision()
    experiment_layer_conditional()
