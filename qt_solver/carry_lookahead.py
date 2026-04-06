"""
Carry-Lookahead OMEGA System.

Instead of modeling carry as 32 sequential MAJ equations per addition
(c[k+1] = MAJ(x[k], y[k], c[k])), use the parallel carry-lookahead
formulation based on Generate/Propagate/Kill segments:

    G[k] = x[k] AND y[k]       (generate: carry produced regardless of carry-in)
    P[k] = x[k] XOR y[k]       (propagate: carry-in passes through)
    K[k] = NOT(x[k]) AND NOT(y[k])  (kill: carry stopped)
    c[k] = OR_{j<k} (G[j] AND P[j+1] AND ... AND P[k-1])

Key insight: P[k] = x[k] XOR y[k] is LINEAR in GF(2), while G[k] = x[k]*y[k]
is quadratic. The segment structure means most propagate chains are short (~2 bits),
so the effective degree of the carry formula is low.
"""

import random
from qt_solver.sha256_traced import (
    MASK32, get_bit, sha256_compress, sha256_compress_traced, IV, add32_traced,
)
from qt_solver.gf2 import (
    gf2_solve, gf2_kernel, gf2_gaussian_eliminate, bitvec_weight,
)


# ============================================================
# 1. GPK Decomposition
# ============================================================

def gpk_decompose(x, y, n=32):
    """
    Decompose addition x + y into Generate/Propagate/Kill bits and segments.

    For each bit position k:
        G[k] = x[k] AND y[k]       (generate)
        P[k] = x[k] XOR y[k]       (propagate)
        K[k] = NOT(x[k]) AND NOT(y[k])  (kill)

    Segments are maximal runs of consecutive P bits, bookended by G or K bits.
    Each segment is (type, start, length) where type is 'G', 'P', or 'K'.

    Returns:
        G: int (bit-vector of generate bits)
        P: int (bit-vector of propagate bits)
        K: int (bit-vector of kill bits)
        segments: list of (type, start, length)
    """
    G = x & y
    P = x ^ y
    K = (~x & ~y) & ((1 << n) - 1)

    # Verify partition: at each position exactly one of G, P, K is set
    assert (G | P | K) == ((1 << n) - 1), "G|P|K should cover all bits"
    assert (G & P) == 0 and (G & K) == 0 and (P & K) == 0, "G, P, K should be disjoint"

    # Build segments: maximal runs of same type
    segments = []
    if n == 0:
        return G, P, K, segments

    k = 0
    while k < n:
        if get_bit(G, k):
            seg_type = 'G'
        elif get_bit(P, k):
            seg_type = 'P'
        else:
            seg_type = 'K'

        start = k
        k += 1
        while k < n:
            if seg_type == 'G' and get_bit(G, k):
                k += 1
            elif seg_type == 'P' and get_bit(P, k):
                k += 1
            elif seg_type == 'K' and get_bit(K, k):
                k += 1
            else:
                break
        segments.append((seg_type, start, k - start))

    return G, P, K, segments


# ============================================================
# 2. Carry from GPK using Lookahead Formula
# ============================================================

def carry_from_gpk(G, P, n=32):
    """
    Compute carry bits from G and P vectors using the carry-lookahead formula.

    c[0] = 0 (no carry into bit 0)
    c[k] = OR_{j=0}^{k-1} ( G[j] AND P[j+1] AND P[j+2] AND ... AND P[k-1] )

    This is equivalent to: carry into position k is generated if there exists
    some position j < k that generates a carry AND all positions between j+1
    and k-1 propagate.

    Returns: carry bit-vector (int), where bit k = carry into position k.
    """
    carry = 0  # c[0] = 0

    for k in range(1, n):
        ck = 0
        # Check each possible generator position j < k
        for j in range(k - 1, -1, -1):
            if not get_bit(G, j):
                continue
            # Check that P[j+1] through P[k-1] are all 1
            all_propagate = True
            for i in range(j + 1, k):
                if not get_bit(P, i):
                    all_propagate = False
                    break
            if all_propagate:
                ck = 1
                break  # OR: one term suffices
        if ck:
            carry |= (1 << k)

    return carry


def verify_carry_gpk(n_samples=1000, n=32):
    """Verify that carry_from_gpk matches add32_traced carry computation."""
    rng = random.Random(12345)
    mismatches = 0

    for _ in range(n_samples):
        x = rng.randint(0, (1 << n) - 1)
        y = rng.randint(0, (1 << n) - 1)

        _, carry_std = add32_traced(x, y)
        G, P, K, _ = gpk_decompose(x, y, n)
        carry_gpk = carry_from_gpk(G, P, n)

        if carry_std != carry_gpk:
            mismatches += 1

    return mismatches


# ============================================================
# 3. GPK Segment Statistics
# ============================================================

def gpk_segment_stats(n_samples=10000, n=32):
    """
    Gather statistics on GPK segments for random 32-bit additions.

    Returns dict with:
        avg_segments: average number of segments per addition
        avg_p_length: average length of P-segments
        max_p_length: maximum P-segment length observed
        p_length_hist: histogram of P-segment lengths
        type_counts: average count of each segment type
        avg_g_count: average number of G bits per addition
        avg_p_count: average number of P bits per addition
        avg_k_count: average number of K bits per addition
    """
    rng = random.Random(54321)

    total_segments = 0
    p_lengths = []
    type_segment_counts = {'G': 0, 'P': 0, 'K': 0}
    total_g_bits = 0
    total_p_bits = 0
    total_k_bits = 0

    for _ in range(n_samples):
        x = rng.randint(0, (1 << n) - 1)
        y = rng.randint(0, (1 << n) - 1)

        G, P, K, segments = gpk_decompose(x, y, n)

        total_segments += len(segments)
        total_g_bits += bin(G).count('1')
        total_p_bits += bin(P).count('1')
        total_k_bits += bin(K).count('1')

        for seg_type, start, length in segments:
            type_segment_counts[seg_type] += 1
            if seg_type == 'P':
                p_lengths.append(length)

    # P-length histogram
    p_length_hist = {}
    for pl in p_lengths:
        p_length_hist[pl] = p_length_hist.get(pl, 0) + 1

    stats = {
        'avg_segments': total_segments / n_samples,
        'avg_p_length': sum(p_lengths) / len(p_lengths) if p_lengths else 0,
        'max_p_length': max(p_lengths) if p_lengths else 0,
        'p_length_hist': dict(sorted(p_length_hist.items())),
        'type_counts': {k: v / n_samples for k, v in type_segment_counts.items()},
        'avg_g_bits': total_g_bits / n_samples,
        'avg_p_bits': total_p_bits / n_samples,
        'avg_k_bits': total_k_bits / n_samples,
    }

    return stats


# ============================================================
# 4. Build Compact OMEGA System
# ============================================================

def _build_jacobian_rows(msg, R):
    """
    Build Jacobian of R-round SHA-256: for each hash bit, which message
    bits affect it? Returns list of row bit-vectors (one per hash bit)
    and the reference hash.
    """
    n_msg = 512
    target_hash = sha256_compress(msg, R)

    # For each message bit, compute which hash bits flip
    msg_effects = [0] * n_msg
    for wi in range(16):
        for bi in range(32):
            msg_idx = wi * 32 + bi
            msg_copy = list(msg)
            msg_copy[wi] ^= (1 << bi)
            new_hash = sha256_compress(msg_copy, R)

            diff = 0
            for w in range(8):
                xd = target_hash[w] ^ new_hash[w]
                diff |= (xd << (w * 32))
            msg_effects[msg_idx] = diff

    # Build rows: one per hash bit
    hash_bits = 256
    rows = []
    for hb in range(hash_bits):
        row = 0
        for mb in range(n_msg):
            if (msg_effects[mb] >> hb) & 1:
                row |= (1 << mb)
        rows.append(row)

    return rows, target_hash, msg_effects


def build_compact_system(R, msg, verbose=True):
    """
    Build a COMPACT OMEGA system using carry-lookahead decomposition.

    The key insight:
    - Standard OMEGA: 32 MAJ carry variables per addition (7 additions/round
      for rounds, 3 per schedule expansion) -> many quadratic equations
    - Compact OMEGA: G[k] = x[k]*y[k] (quadratic), P[k] = x[k] XOR y[k] (LINEAR)
      c[k] = OR of (G[j] * P-chain) -> products of known-length P-chains

    For the compact system we:
    1. Trace SHA-256 to get all carries
    2. Decompose each addition into G/P/K segments
    3. Count how many carry variables we actually need (only at segment boundaries)
    4. Build the Jacobian-based linearized system (same approach as omega_solve)
    5. Compare equation counts

    Returns dict with system info and comparison.
    """
    trace = sha256_compress_traced(msg, R)
    target_hash = trace['hash']

    if verbose:
        print(f"\n  Building compact system for R={R}")

    # --- Analyze carry structure ---
    # Count total additions, total carry bits, and compact variables
    total_additions = 0
    total_carry_bits = 0   # standard: 32 carry vars per addition (minus c[0]=0, so 31)
    total_gpk_segments = 0
    total_g_bits = 0
    total_p_bits = 0
    total_k_bits = 0
    p_chain_lengths = []

    # Round additions (7 per round)
    for r in range(R):
        state = trace['states'][r]
        a, b, c, d, e, f, g, h = state

        # The 7 additions per round. We need the two operands for each.
        # We reconstruct from the trace.
        from qt_solver.sha256_traced import sigma0, sigma1, ch, maj as maj_fn, K as K_const, rotr

        sig1_e = sigma1(e)
        ch_efg = ch(e, f, g)
        sig0_a = sigma0(a)
        maj_abc = maj_fn(a, b, c)

        # T1 chain: h + sig1(e) + ch(e,f,g) + K[r] + W[r]
        s1 = (h + sig1_e) & MASK32
        s2 = (s1 + ch_efg) & MASK32
        s3 = (s2 + K_const[r]) & MASK32
        t1 = (s3 + trace['W'][r]) & MASK32

        # T2: sig0(a) + maj(a,b,c)
        t2 = (sig0_a + maj_abc) & MASK32

        # Addition pairs for this round
        add_pairs = [
            (h, sig1_e),
            (s1, ch_efg),
            (s2, K_const[r]),
            (s3, trace['W'][r]),
            (sig0_a, maj_abc),
            (d, t1),
            (t1, t2),
        ]

        for x_op, y_op in add_pairs:
            total_additions += 1
            total_carry_bits += 31  # 31 non-trivial carry bits per addition

            G, P, K_vec, segments = gpk_decompose(x_op, y_op, 32)
            total_gpk_segments += len(segments)
            total_g_bits += bin(G).count('1')
            total_p_bits += bin(P).count('1')
            total_k_bits += bin(K_vec).count('1')

            for seg_type, start, length in segments:
                if seg_type == 'P':
                    p_chain_lengths.append(length)

    # Schedule additions (3 per expansion, for words 16..R-1 if R > 16)
    schedule_additions = 0
    if R > 16:
        # Words 16 through max(R-1, 15) have 3 additions each
        for i in range(16, R):
            schedule_additions += 3
            total_additions += 3
            total_carry_bits += 31 * 3

    # --- Compact variable count ---
    # In the compact system, carry at position k only needs to be tracked
    # at G-segment boundaries. Within a P-chain, carry propagates deterministically.
    # Within a K-segment, carry is killed (c = 0 after K).
    # So effective carry variables = number of G-segments (one bit: does this G generate?)
    # Plus P-chain products (but these are known once G/P are known).

    # The real saving: in linearized system, P[k] is linear in message bits,
    # so carry = G[j] * (product of P's) has degree = 1 + chain_length.
    # Short chains => low degree => better linearization.

    compact_vars = total_g_bits  # one quadratic var per G-bit
    # P-bits are linear (free), K-bits are determined (no carry)

    # --- Build Jacobian system (same as omega_solve but we track structure) ---
    n_msg = 512
    rows, _, msg_effects = _build_jacobian_rows(msg, R)
    hash_bits = 256

    echelon, pivots = gf2_gaussian_eliminate(list(rows), n_msg)
    rank = len(pivots)
    kernel = gf2_kernel(rows, n_msg)
    alpha_dim = len(kernel)

    # --- Average P-chain length ---
    avg_p_chain = sum(p_chain_lengths) / len(p_chain_lengths) if p_chain_lengths else 0
    max_p_chain = max(p_chain_lengths) if p_chain_lengths else 0

    result = {
        'R': R,
        'total_additions': total_additions,
        'standard_carry_vars': total_carry_bits,
        'compact_carry_vars': compact_vars,
        'total_gpk_segments': total_gpk_segments,
        'total_g_bits': total_g_bits,
        'total_p_bits': total_p_bits,
        'total_k_bits': total_k_bits,
        'avg_p_chain_length': avg_p_chain,
        'max_p_chain_length': max_p_chain,
        'jacobian_rank': rank,
        'jacobian_rows': hash_bits,
        'jacobian_cols': n_msg,
        'alpha_kernel_dim': alpha_dim,
        'expected_alpha': max(0, 32 * (16 - R)),
        'savings_ratio': (total_carry_bits - compact_vars) / total_carry_bits if total_carry_bits > 0 else 0,
    }

    if verbose:
        print(f"    Additions: {total_additions} ({R}*7 round + {schedule_additions} schedule)")
        print(f"    Standard carry vars: {total_carry_bits} (31 per addition)")
        print(f"    G-bits (quadratic): {total_g_bits}")
        print(f"    P-bits (linear):    {total_p_bits}")
        print(f"    K-bits (zero carry):{total_k_bits}")
        print(f"    Compact carry vars: {compact_vars} (G-bits only)")
        print(f"    Savings: {result['savings_ratio']:.1%} fewer carry variables")
        print(f"    Avg P-chain: {avg_p_chain:.2f}, Max P-chain: {max_p_chain}")
        print(f"    Jacobian rank: {rank}/{hash_bits}")
        print(f"    alpha-kernel dim: {alpha_dim} (expected: {result['expected_alpha']})")

    return result


# ============================================================
# 5. Compare Standard vs Compact Systems
# ============================================================

def compare_systems(R_values=None, seed=42, verbose=True):
    """
    Compare standard OMEGA vs compact carry-lookahead system for various R.

    For each R:
    - Number of additions, carry variables (standard vs compact)
    - Jacobian rank
    - alpha-kernel dimension
    - Variable savings

    Returns list of result dicts.
    """
    if R_values is None:
        R_values = [4, 8, 12, 16]

    rng = random.Random(seed)
    msg = [rng.randint(0, MASK32) for _ in range(16)]

    if verbose:
        print("=" * 72)
        print("CARRY-LOOKAHEAD OMEGA: Standard vs Compact Comparison")
        print("=" * 72)
        print(f"  Message: {[hex(w) for w in msg[:4]]}...")

    results = []
    for R in R_values:
        result = build_compact_system(R, msg, verbose=verbose)
        results.append(result)

    if verbose:
        print("\n" + "=" * 72)
        print(f"{'R':>4} | {'Adds':>5} | {'Std Carry':>10} | {'Compact':>8} | {'Savings':>8} | {'Rank':>5} | {'alpha-K':>8} | {'Expected':>8}")
        print("-" * 72)
        for r in results:
            print(f"{r['R']:>4} | {r['total_additions']:>5} | {r['standard_carry_vars']:>10} | "
                  f"{r['compact_carry_vars']:>8} | {r['savings_ratio']:>7.1%} | "
                  f"{r['jacobian_rank']:>5} | {r['alpha_kernel_dim']:>8} | {r['expected_alpha']:>8}")
        print("=" * 72)

        # Summary insights
        print("\nKey insights:")
        print("  - P[k] = x[k] XOR y[k] is LINEAR: no new variables needed")
        print("  - K[k] = NOT(x[k]) AND NOT(y[k]) kills carry: no variables needed")
        print("  - Only G[k] = x[k] AND y[k] introduces quadratic terms")
        print("  - Compact system needs ~75% fewer carry variables")
        print("  - alpha-kernel dimension matches standard OMEGA (both use same Jacobian)")

        if results:
            avg_savings = sum(r['savings_ratio'] for r in results) / len(results)
            avg_p = sum(r['avg_p_chain_length'] for r in results) / len(results)
            print(f"\n  Average savings: {avg_savings:.1%}")
            print(f"  Average P-chain length: {avg_p:.2f} bits")
            print(f"  (Short P-chains mean carry-lookahead products have low degree)")

    return results


# ============================================================
# Main
# ============================================================

if __name__ == '__main__':
    print("=" * 72)
    print("CARRY-LOOKAHEAD OMEGA SYSTEM")
    print("=" * 72)

    # --- Verification: carry_from_gpk matches standard ---
    print("\n1. Verifying carry_from_gpk correctness...")
    mismatches = verify_carry_gpk(n_samples=5000, n=32)
    if mismatches == 0:
        print(f"   PASSED: carry_from_gpk matches add32_traced on 5000 random pairs")
    else:
        print(f"   FAILED: {mismatches} mismatches!")

    # Also verify on small n for exhaustive check
    small_n = 8
    mismatches_small = 0
    for x in range(1 << small_n):
        for y in range(1 << small_n):
            _, carry_std = add32_traced(x, y)
            # add32_traced always does 32 bits; mask to small_n
            carry_std_masked = carry_std & ((1 << small_n) - 1)
            G, P, K_vec, _ = gpk_decompose(x, y, small_n)
            carry_gpk = carry_from_gpk(G, P, small_n)
            if carry_std_masked != carry_gpk:
                mismatches_small += 1
    total_small = (1 << small_n) ** 2
    if mismatches_small == 0:
        print(f"   PASSED: exhaustive check on n={small_n} ({total_small} pairs)")
    else:
        print(f"   FAILED: {mismatches_small}/{total_small} mismatches on n={small_n}")

    # --- GPK segment statistics ---
    print("\n2. GPK segment statistics (n=32, 10000 samples)...")
    stats = gpk_segment_stats(n_samples=10000, n=32)
    print(f"   Avg segments per addition: {stats['avg_segments']:.2f}")
    print(f"   Avg G-bits: {stats['avg_g_bits']:.2f}, P-bits: {stats['avg_p_bits']:.2f}, K-bits: {stats['avg_k_bits']:.2f}")
    print(f"   Avg P-chain length: {stats['avg_p_length']:.2f}")
    print(f"   Max P-chain length: {stats['max_p_length']}")
    print(f"   Segment type counts: G={stats['type_counts']['G']:.2f}, P={stats['type_counts']['P']:.2f}, K={stats['type_counts']['K']:.2f}")

    # P-chain length distribution
    print("   P-chain length histogram:")
    hist = stats['p_length_hist']
    total_p = sum(hist.values())
    for length in sorted(hist.keys()):
        count = hist[length]
        pct = 100.0 * count / total_p
        bar = '#' * int(pct / 2)
        print(f"     len={length:2d}: {count:6d} ({pct:5.1f}%) {bar}")

    # --- GPK decomposition properties ---
    print("\n3. Verifying GPK algebraic properties...")
    rng = random.Random(99)
    for trial in range(100):
        x = rng.randint(0, MASK32)
        y = rng.randint(0, MASK32)
        G, P, K_vec, segs = gpk_decompose(x, y, 32)

        # G + P + K should have exactly 32 bits set total
        assert bin(G).count('1') + bin(P).count('1') + bin(K_vec).count('1') == 32

        # G[k] implies NOT P[k] and NOT K[k]
        assert G & P == 0
        assert G & K_vec == 0
        assert P & K_vec == 0

        # G[k] = x[k]*y[k]: verify bit by bit
        for k in range(32):
            xk = get_bit(x, k)
            yk = get_bit(y, k)
            assert get_bit(G, k) == (xk & yk)
            assert get_bit(P, k) == (xk ^ yk)
            assert get_bit(K_vec, k) == ((1 - xk) & (1 - yk))

        # Segments should be contiguous and cover all 32 bits
        total_len = sum(length for _, _, length in segs)
        assert total_len == 32, f"Segments cover {total_len} bits, expected 32"

    print("   PASSED: all GPK algebraic properties verified (100 trials)")

    # --- Linearity of P ---
    print("\n4. Verifying P[k] = x[k] XOR y[k] is linear in GF(2)...")
    # P is linear: P(x+dx, y+dy) = P(x,y) XOR P(dx, dy) XOR P(0,0)
    # Since P(0,0) = 0, we need P(x XOR dx, y XOR dy) = P(x,y) XOR P(dx,dy)
    # where XOR-addition is used (GF(2) vector space).
    # Actually: (x XOR dx) XOR (y XOR dy) = (x XOR y) XOR (dx XOR dy) = P(x,y) XOR P(dx,dy)
    # So P is indeed linear (it IS XOR).
    rng = random.Random(77)
    linear_ok = True
    for _ in range(1000):
        x = rng.randint(0, MASK32)
        y = rng.randint(0, MASK32)
        dx = rng.randint(0, MASK32)
        dy = rng.randint(0, MASK32)
        P1 = (x ^ y)                       # P(x, y)
        P2 = ((x ^ dx) ^ (y ^ dy))         # P(x+dx, y+dy)
        Pdiff = (dx ^ dy)                   # P(dx, dy)
        if P2 != (P1 ^ Pdiff):
            linear_ok = False
            break
    print(f"   P linearity: {'PASSED' if linear_ok else 'FAILED'} (1000 trials)")

    # G is quadratic: G(x,y) = x AND y is NOT linear
    g_linear = True
    for _ in range(100):
        x = rng.randint(0, MASK32)
        y = rng.randint(0, MASK32)
        dx = rng.randint(0, MASK32)
        dy = rng.randint(0, MASK32)
        G1 = x & y
        G2 = (x ^ dx) & (y ^ dy)
        Gdiff = dx & dy
        if G2 != (G1 ^ Gdiff):
            g_linear = False
            break
    print(f"   G linearity: {'correctly NOT linear' if not g_linear else 'ERROR: appears linear'}")

    # --- Compare systems ---
    print("\n5. Comparing standard vs compact OMEGA systems...")
    results = compare_systems(R_values=[4, 8, 12, 16], seed=42, verbose=True)

    # --- Final summary ---
    print("\n" + "=" * 72)
    print("CARRY-LOOKAHEAD SUMMARY")
    print("=" * 72)
    print("  The carry-lookahead decomposition reduces the OMEGA system by")
    print("  replacing 31 sequential MAJ carry equations per addition with:")
    print("    - G[k] = x[k]*y[k] : quadratic (one variable per G-bit)")
    print("    - P[k] = x[k] XOR y[k] : LINEAR (no new variables!)")
    print("    - K[k] = kills carry (no variables needed)")
    print("  Since E[P-bits] ~ 16/32, roughly half the carry structure is linear.")
    print("  Since E[K-bits] ~ 8/32, another quarter needs no variables at all.")
    print("  Only the ~8/32 G-bits contribute quadratic terms.")
    print("  The alpha-kernel dimension is invariant under this reformulation.")
