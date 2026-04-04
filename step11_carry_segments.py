"""
Step 11: Carry-Segment Algebra — keep carry as STRUCTURE, not bits

Problem: GF(2) expands carry into 2^n monomials. Z/2^n hides carry entirely.
We need something in between.

Key observation from Step 5: carry chain = MAJ recurrence.
  carry[0] = 0
  carry[k] = MAJ(a[k-1], b[k-1], carry[k-1])

This is a RECURSIVE structure. The recursion creates SEGMENTS:
consecutive bits where carry propagates (G), generates (G), or kills (K).

From BTE: mean segment length = 3.64, geometric distribution.

IDEA: Instead of expanding carry into GF(2) monomials,
represent it as SEGMENTS. Each segment is a compact object:
  - Generate segment: carry = 1 regardless of input carry
  - Propagate segment: carry_out = carry_in
  - Kill segment: carry = 0 regardless

For addition a + b, the GPK classification at each bit:
  G(k): a[k]=1 AND b[k]=1 → carry_out = 1
  P(k): a[k] XOR b[k] = 1 → carry_out = carry_in
  K(k): a[k]=0 AND b[k]=0 → carry_out = 0

The CARRY at bit k = determined by the rightmost G or K below k.

THIS is the "native" representation: carry is not a degree-2^k polynomial,
it's a SCAN OPERATION on a GPK string.

Let's formulate the collision in GPK language and see if it stays structured.
"""

import numpy as np
from step0_exact_algebra import mini_sha, N, MASK

N_MSG = 4
N_INPUT = N * N_MSG
N_TOTAL = 1 << N_INPUT


def gpk_classify(a, b, n=N):
    """Classify each bit position as G, P, or K for addition a + b."""
    result = []
    for k in range(n):
        ak = (a >> k) & 1
        bk = (b >> k) & 1
        if ak == 1 and bk == 1:
            result.append('G')
        elif ak == 0 and bk == 0:
            result.append('K')
        else:
            result.append('P')
    return result


def carry_from_gpk(gpk):
    """Compute carry bits from GPK string."""
    carries = [0]  # carry into bit 0 = 0
    for k in range(len(gpk)):
        if gpk[k] == 'G':
            carries.append(1)
        elif gpk[k] == 'K':
            carries.append(0)
        else:  # P
            carries.append(carries[-1])
    return carries[1:]  # carry OUT of each position


def gpk_segments(gpk):
    """Break GPK string into segments of consecutive same type."""
    if not gpk:
        return []
    segments = []
    current = gpk[0]
    length = 1
    start = 0
    for k in range(1, len(gpk)):
        if gpk[k] == current:
            length += 1
        else:
            segments.append((current, start, length))
            current = gpk[k]
            length = 1
            start = k
    segments.append((current, start, length))
    return segments


def analyze_gpk_collision():
    """
    For the collision problem: when f(M) = f(M⊕δ),
    what do the GPK patterns of the intermediate additions look like?
    Is there GPK-level structure that survives more rounds than GF(2)?
    """
    import random
    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(N_MSG)]

    print(f"{'='*80}")
    print(f"GPK ANALYSIS OF COLLISION LANDSCAPE")
    print(f"M = {[hex(w) for w in M]}")
    print(f"{'='*80}")

    # For each round, during the computation of T1+T2:
    # The GPK pattern of (T1, T2) determines the carry.
    # For a collision pair (M, M⊕δ):
    #   GPK_1 = gpk(T1(M), T2(M))
    #   GPK_2 = gpk(T1(M⊕δ), T2(M⊕δ))
    # Question: how similar are GPK_1 and GPK_2 for collision δ?

    for R in [4, 6, 8]:
        print(f"\n{'='*80}")
        print(f"R = {R}")
        print(f"{'='*80}")

        # Run base computation, track T1 and T2 at each round
        def run_with_intermediates(msg, R):
            a, b, c, d, e, f, g, h = 0x6, 0xB, 0x3, 0xA, 0x5, 0x9, 0x1, 0xF
            W = list(msg) + [0] * max(0, R - len(msg))
            K = [0x4, 0x2, 0xB, 0x7, 0xA, 0x3, 0xE, 0x5,
                 0x9, 0x1, 0xD, 0x6, 0x0, 0x8, 0xC, 0xF]
            intermediates = []  # (T1, T2) at each round

            for r in range(R):
                w_r = W[r] if r < len(W) else 0
                k_r = K[r % len(K)]
                def rotr(x, s):
                    return ((x >> s) | (x << (N - s))) & MASK
                sig1 = rotr(e, 1) ^ rotr(e, 3) ^ (e >> 1)
                ch_val = (e & f) ^ (~e & g) & MASK
                sig0 = rotr(a, 1) ^ rotr(a, 2) ^ rotr(a, 3)
                maj_val = (a & b) ^ (a & c) ^ (b & c)
                T1 = (h + sig1 + ch_val + k_r + w_r) & MASK
                T2 = (sig0 + maj_val) & MASK
                intermediates.append((T1, T2))
                a_new = (T1 + T2) & MASK
                e_new = (d + T1) & MASK
                h, g, f = g, f, e
                e = e_new
                d, c, b = c, b, a
                a = a_new
            return (a, e), intermediates

        base_out, base_inter = run_with_intermediates(M, R)

        # Find collision δs
        collisions = []
        for didx in range(1, N_TOTAL):
            dw = []
            tmp = didx
            for w in range(N_MSG):
                dw.append(tmp & MASK)
                tmp >>= N
            M2 = [(M[w] ^ dw[w]) for w in range(N_MSG)]
            out2, inter2 = run_with_intermediates(M2, R)
            if out2 == base_out:
                collisions.append((didx, inter2))

        print(f"  Collisions found: {len(collisions)}")

        if not collisions:
            continue

        # Analyze GPK patterns for collisions vs random
        print(f"\n  GPK PATTERN COMPARISON (per round):")
        print(f"  {'Round':>6} {'GPK_match%':>11} {'Seg_match%':>11} "
              f"{'P_count_base':>13} {'P_count_coll':>13} {'δP':>5}")

        for r in range(R):
            T1_base, T2_base = base_inter[r]
            gpk_base = gpk_classify(T1_base, T2_base)

            gpk_matches = 0
            total_bits = 0
            p_count_base = gpk_base.count('P')
            p_counts_coll = []
            seg_matches = 0
            total_segs = 0

            for didx, inter2 in collisions[:50]:  # sample
                T1_c, T2_c = inter2[r]
                gpk_coll = gpk_classify(T1_c, T2_c)

                # Bit-level GPK match
                for k in range(N):
                    if gpk_base[k] == gpk_coll[k]:
                        gpk_matches += 1
                    total_bits += 1

                p_counts_coll.append(gpk_coll.count('P'))

                # Segment-level match
                segs_base = gpk_segments(gpk_base)
                segs_coll = gpk_segments(gpk_coll)
                for sb in segs_base:
                    if sb in segs_coll:
                        seg_matches += 1
                    total_segs += 1

            gpk_pct = gpk_matches / max(total_bits, 1) * 100
            seg_pct = seg_matches / max(total_segs, 1) * 100
            avg_p_coll = np.mean(p_counts_coll) if p_counts_coll else 0
            delta_p = avg_p_coll - p_count_base

            print(f"  R={r:>3}  {gpk_pct:>10.1f}% {seg_pct:>10.1f}% "
                  f"{p_count_base:>12} {avg_p_coll:>12.1f} {delta_p:>+5.2f}")

        # Compare with RANDOM non-collision δs
        print(f"\n  GPK MATCH: collision δ vs random δ")
        print(f"  {'Round':>6} {'Collision':>10} {'Random':>10} {'Excess':>8}")

        for r in range(R):
            T1_base, T2_base = base_inter[r]
            gpk_base = gpk_classify(T1_base, T2_base)

            # Collision GPK match
            coll_match = 0
            coll_total = 0
            for didx, inter2 in collisions[:50]:
                T1_c, T2_c = inter2[r]
                gpk_coll = gpk_classify(T1_c, T2_c)
                for k in range(N):
                    if gpk_base[k] == gpk_coll[k]:
                        coll_match += 1
                    coll_total += 1
            coll_pct = coll_match / max(coll_total, 1) * 100

            # Random GPK match
            rand_match = 0
            rand_total = 0
            for _ in range(50):
                didx = rng.randint(1, N_TOTAL - 1)
                dw = []
                tmp = didx
                for w in range(N_MSG):
                    dw.append(tmp & MASK)
                    tmp >>= N
                M2 = [(M[w] ^ dw[w]) for w in range(N_MSG)]
                _, inter2 = run_with_intermediates(M2, R)
                T1_r, T2_r = inter2[r]
                gpk_rand = gpk_classify(T1_r, T2_r)
                for k in range(N):
                    if gpk_base[k] == gpk_rand[k]:
                        rand_match += 1
                    rand_total += 1
            rand_pct = rand_match / max(rand_total, 1) * 100

            excess = coll_pct - rand_pct
            marker = " ★" if abs(excess) > 5 else ""
            print(f"  R={r:>3}  {coll_pct:>9.1f}% {rand_pct:>9.1f}% {excess:>+7.1f}%{marker}")

    # ================================================================
    # KEY TEST: GPK as compact representation
    # ================================================================
    print(f"\n{'='*80}")
    print(f"GPK COMPRESSION: Can we represent collision in GPK space?")
    print(f"{'='*80}")

    R = 6
    base_out, base_inter = run_with_intermediates(M, R)

    # For each collision δ, compute the GPK SIGNATURE:
    # the sequence of GPK strings across all rounds
    # Is this signature more structured than the bit-level representation?

    gpk_signatures = {}  # signature → count
    for didx in range(1, N_TOTAL):
        dw = []
        tmp = didx
        for w in range(N_MSG):
            dw.append(tmp & MASK)
            tmp >>= N
        M2 = [(M[w] ^ dw[w]) for w in range(N_MSG)]
        out2, inter2 = run_with_intermediates(M2, R)

        if out2 != base_out:
            continue

        # GPK signature = tuple of GPK strings for all rounds
        sig = []
        for r in range(R):
            T1_c, T2_c = inter2[r]
            gpk = tuple(gpk_classify(T1_c, T2_c))
            sig.append(gpk)
        sig = tuple(sig)
        gpk_signatures[sig] = gpk_signatures.get(sig, 0) + 1

    print(f"\n  R={R}: {len(collisions)} collisions, {len(gpk_signatures)} unique GPK signatures")
    print(f"  Compression ratio: {len(collisions)}/{len(gpk_signatures)} = "
          f"{len(collisions)/max(len(gpk_signatures),1):.2f}")
    print(f"  If every collision has unique GPK signature → no GPK-level degeneracy")
    print(f"  If many share a signature → GPK captures equivalence classes")

    # Top GPK signatures
    top = sorted(gpk_signatures.items(), key=lambda x: -x[1])[:5]
    for sig, count in top:
        print(f"    Count={count}: GPK[R-1]={''.join(sig[-1])}")


if __name__ == "__main__":
    analyze_gpk_collision()
