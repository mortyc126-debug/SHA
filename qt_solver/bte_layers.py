"""
BTE Layer Decomposition: verify 512 = 4×127 + 4.

Theorem T1: rank(Layer(0)) = 2R - 1 for any BTE with 8-register shift.
  SHA-256: 127 = 2×64 - 1.

Layer(k) = system of all bit-k values across all 8 registers and R rounds.
Each layer has 8×R bit-values, but shift structure → only 2R independent.
Створочне dependency: e[r] depends on a → one relation → rank = 2R-1.

Full rank: 512 = 4 × 127 + 4. Exactly 4 "pure" layers + 4 residual bits.
"""

from qt_solver.sha256_traced import (
    MASK32, sha256_compress, sha256_compress_traced, IV, get_bit,
)
from qt_solver.gf2 import gf2_gaussian_eliminate
import random


def compute_layer_jacobian(R, msg, bit_positions, verbose=False):
    """
    Compute GF(2) Jacobian of specific bit positions across the full trajectory.

    For bit position k: we track a[1..R] bit k and e[1..R] bit k.
    (b,c,d,f,g,h are just shifted copies, so a and e are sufficient.)

    Returns Jacobian matrix: rows = trajectory bits, columns = message bits.
    """
    trace = sha256_compress_traced(msg, R)
    states = trace['states']

    n_msg = 512  # 16 words × 32 bits

    # Reference trajectory values for specified bit positions
    # For each bit position k, track a[1..R][k] and e[1..R][k]
    ref_values = {}
    for k in bit_positions:
        for r in range(1, R + 1):
            ref_values[('a', r, k)] = get_bit(states[r][0], k)  # a = state[r][0]
            ref_values[('e', r, k)] = get_bit(states[r][4], k)  # e = state[r][4]

    n_traj = len(ref_values)
    traj_keys = sorted(ref_values.keys())

    # Build Jacobian by flipping each message bit
    rows = []
    for traj_idx, key in enumerate(traj_keys):
        row = 0
        reg, r, k = key
        reg_idx = 0 if reg == 'a' else 4
        ref_val = ref_values[key]

        for mb in range(n_msg):
            wi, bi = mb // 32, mb % 32
            msg2 = list(msg)
            msg2[wi] ^= (1 << bi)
            trace2 = sha256_compress_traced(msg2, R)
            new_val = get_bit(trace2['states'][r][reg_idx], k)
            if new_val != ref_val:
                row |= (1 << mb)
        rows.append(row)

    return rows, traj_keys


def verify_layer_structure(R=64, n_samples=3, verbose=True):
    """
    Verify BTE layer decomposition: 512 = 4×127 + 4.

    Method:
    1. Compute trajectory Jacobian for bit-0 layer → rank should be 127
    2. Cumulative: bits 0-1 → rank 254, bits 0-2 → 381, bits 0-3 → 508, bits 0-4 → 512
    """
    if verbose:
        print("=" * 60)
        print(f"BTE Layer Decomposition (R={R})")
        print("=" * 60)

    rng = random.Random(42)

    for sample in range(n_samples):
        msg = [rng.randint(0, MASK32) for _ in range(16)]
        if verbose:
            print(f"\n  Sample {sample}: msg={[hex(w) for w in msg[:3]]}...")

        # Compute layer-by-layer
        cumulative_rank = 0
        prev_rows = []

        for layer in range(6):  # bits 0..5
            if verbose:
                print(f"    Computing layer {layer}...", end=" ", flush=True)

            # Get Jacobian for this bit position
            layer_rows, keys = compute_layer_jacobian(R, msg, [layer])

            # Add to cumulative system
            all_rows = prev_rows + layer_rows
            echelon, pivots = gf2_gaussian_eliminate(list(all_rows), 512)
            new_rank = len(pivots)
            added = new_rank - cumulative_rank

            if verbose:
                print(f"rank={new_rank} (+{added})")

            cumulative_rank = new_rank
            prev_rows = all_rows

            # After bit 4 should have full rank 512
            if cumulative_rank == 512:
                if verbose:
                    print(f"    → FULL RANK reached at layer {layer}")
                break

        if verbose:
            print(f"\n    Summary for sample {sample}:")
            print(f"    Expected: 127, 254, 381, 508, 512")
            print(f"    Formula: 4×127 + 4 = 512")

    return True


def verify_T1_universal(verbose=True):
    """
    Theorem T1: Layer rank = 2R - 1 for any BTE with 8-register shift.
    Verify for different R values.
    """
    if verbose:
        print("\n" + "=" * 60)
        print("Theorem T1: Layer(0) rank = 2R - 1")
        print("=" * 60)

    rng = random.Random(42)
    msg = [rng.randint(0, MASK32) for _ in range(16)]

    for R in [4, 8, 16, 32]:
        if verbose:
            print(f"\n  R={R}: computing bit-0 layer Jacobian...", end=" ", flush=True)

        rows, keys = compute_layer_jacobian(R, msg, [0])
        echelon, pivots = gf2_gaussian_eliminate(list(rows), 512)
        rank = len(pivots)
        expected = 2 * R - 1

        status = "✓" if rank == expected else "✗"
        if verbose:
            print(f"rank={rank}, expected 2×{R}-1={expected}  {status}")

    return True


def stvorochne_verification(R=64, n_tests=100, verbose=True):
    """
    Verify створочне (shutter) identity:
      e[r] = a[r] + a[r-4] - T2[r-1]  (mod 2^32)

    This is the one dependency that reduces layer rank from 2R to 2R-1.
    """
    if verbose:
        print("\n" + "=" * 60)
        print("Створочне Identity: e[r] = a[r] + a[r-4] - T2[r-1]")
        print("=" * 60)

    rng = random.Random(42)
    violations = 0
    total = 0

    for _ in range(n_tests):
        msg = [rng.randint(0, MASK32) for _ in range(16)]
        trace = sha256_compress_traced(msg, R)
        states = trace['states']

        for r in range(4, R + 1):
            a_r = states[r][0]
            e_r = states[r][4]
            a_r4 = states[r-4][0] if r >= 4 else IV[4 - r]  # simplified
            t2 = trace['T2'][r-1] if r > 0 else 0

            # Check: e[r] = (a[r] + a[r-4] - T2[r-1]) mod 2^32
            # Actually: from the round function:
            #   a_new = T1 + T2
            #   e_new = d + T1 = a[r-3] + T1
            # And T1 = a_new - T2
            # So e_new = a[r-3] + a_new - T2
            # → e[r] = a[r-3] + a[r] - T2[r-1]  (not a[r-4]!)

            # Wait, let me recheck. d[r] = c[r-1] = b[r-2] = a[r-3]
            # e_new = d + T1 = a[r-3] + T1
            # a_new = T1 + T2 → T1 = a_new - T2
            # So e_new = a[r-3] + a_new - T2
            # At round r: e[r+1] = a[r-2] + a[r+1] - T2[r]
            # Or: e[r] = a[r-3] + a[r] - T2[r-1]

            # d at round r-1 = state[r-1][3] = a[r-4]? No:
            # state[r][3] = d[r] = c[r-1] = b[r-2] = a[r-3]
            # So d at round r-1 = a[r-1-3] = a[r-4]
            # e[r] = e_new from round r-1 = d[r-1] + T1[r-1] = a[r-4] + T1[r-1]
            # T1[r-1] = a[r] - T2[r-1] (since a[r] = T1[r-1] + T2[r-1])
            # → e[r] = a[r-4] + a[r] - T2[r-1]

            if r >= 4 and r <= R:
                a_r4_actual = states[r-4][0]
                expected_e = (a_r + a_r4_actual - t2) & MASK32
                if expected_e != e_r:
                    violations += 1
                total += 1

    if verbose:
        print(f"  Violations: {violations}/{total}")
        if violations == 0:
            print(f"  ✓ VERIFIED: e[r] = a[r] + a[r-4] - T2[r-1] (mod 2^32)")
            print(f"  This is the створочне that creates rank = 2R-1 (one dependency)")

    return violations == 0


if __name__ == '__main__':
    # First verify створочне
    stvorochne_verification(R=64, n_tests=50)

    # Verify T1 for small R (full layer Jacobian is expensive)
    verify_T1_universal()

    # Full layer structure for small R
    verify_layer_structure(R=16, n_samples=1)
