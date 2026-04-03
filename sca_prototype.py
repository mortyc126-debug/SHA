"""
SEGMENTED CARRY ALGEBRA (SCA) — First prototype.

The idea: a 32-bit word in SHA-256 is not a monolith (GF(2)^32)
and not an atom (Z/2^32). It's a CHAIN OF SEGMENTS separated
by carry-breaks.

A carry-break at position k means: carry[k] = 0 regardless of
what happens at positions 0..k-1. This makes bits k..next_break
INDEPENDENT of bits 0..k-1.

The SCA representation of a 32-bit addition x + y:
  1. Compute G[i] = x[i] AND y[i], P[i] = x[i] XOR y[i], K[i] = NOT(G[i] OR P[i])
  2. Find segments: consecutive runs of P's terminated by G or K
  3. Within each segment: carry propagation is LOCAL
  4. Between segments: NO carry coupling

This gives a FACTORED representation of the addition.

Key question: Does this factoring SURVIVE through SHA-256 rounds?
Or do rotations destroy the segment boundaries?
"""

import random
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, schedule
)


class SCAWord:
    """
    Segmented Carry Algebra representation of a 32-bit word.

    A word is represented by its VALUE plus the SEGMENT BOUNDARIES
    from the most recent addition that produced it.

    Segment boundary at position k means: carry chain broke here.
    Bits within a segment are carry-coupled.
    Bits in different segments are carry-independent.
    """

    def __init__(self, value, segments=None):
        self.value = value & MASK
        # segments = list of (start, length) tuples
        if segments is None:
            # No addition history → treat as single segment
            self.segments = [(0, 32)]
        else:
            self.segments = segments

    @staticmethod
    def from_addition(x_val, y_val):
        """Create SCAWord from x + y, computing segment boundaries."""
        result = (x_val + y_val) & MASK

        # Compute G, P, K per bit
        segments = []
        seg_start = 0
        carry = 0

        for k in range(32):
            xk = (x_val >> k) & 1
            yk = (y_val >> k) & 1
            g = xk & yk  # generate
            p = xk ^ yk  # propagate

            new_carry = g | (p & carry)

            # A segment BREAKS when carry becomes 0 after being nonzero
            # OR when we hit a Kill (G=0, P=0) with no incoming carry
            if not new_carry and carry:
                # Carry chain just ended
                segments.append((seg_start, k + 1 - seg_start))
                seg_start = k + 1
            elif not new_carry and not carry and not g and not p:
                # Kill position with no carry — also a boundary
                if k > seg_start:
                    segments.append((seg_start, k - seg_start))
                    seg_start = k

            carry = new_carry

        # Final segment
        if seg_start < 32:
            segments.append((seg_start, 32 - seg_start))

        return SCAWord(result, segments)

    def n_segments(self):
        return len(self.segments)

    def max_segment_len(self):
        return max(l for _, l in self.segments) if self.segments else 0

    def segment_lengths(self):
        return [l for _, l in self.segments]

    def __repr__(self):
        segs = ','.join(f'{s}:{l}' for s, l in self.segments)
        return f'SCA(0x{self.value:08x}, [{segs}])'


def experiment_sca_basic():
    """
    Basic SCA: compute segment statistics for SHA-256 additions.
    """
    print("=" * 80)
    print("SCA PROTOTYPE: Segment statistics across SHA-256 rounds")
    print("=" * 80)

    N = 200
    all_n_segs = []
    all_max_lens = []

    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        states, W = sha256_round_trace(M)

        for r in range(64):
            a, b, c, d, e, f, g, h = states[r]

            # The final addition: T1 + T2 → a_new
            T1 = (h + Sig1(e) + Ch(e, f, g) + K[r] + W[r]) & MASK
            T2 = (Sig0(a) + Maj(a, b, c)) & MASK

            sca = SCAWord.from_addition(T1, T2)
            all_n_segs.append(sca.n_segments())
            all_max_lens.append(sca.max_segment_len())

    avg_segs = sum(all_n_segs) / len(all_n_segs)
    avg_max = sum(all_max_lens) / len(all_max_lens)

    print(f"\n  Across {N} messages × 64 rounds:")
    print(f"  Average segments per word: {avg_segs:.2f}")
    print(f"  Average max segment length: {avg_max:.2f}")

    # Distribution of segment counts
    from collections import Counter
    seg_dist = Counter(all_n_segs)
    print(f"\n  Distribution of segment count per word:")
    for n in sorted(seg_dist.keys()):
        pct = seg_dist[n] / len(all_n_segs) * 100
        print(f"    {n:>2} segments: {pct:>5.1f}%")


def experiment_sca_through_rounds():
    """
    KEY TEST: Do segments SURVIVE through rounds?

    After round r produces a_new via T1+T2, the segments are defined.
    In round r+1, a_new becomes one of the inputs to new additions
    (through Sig0, Maj for T2, or through the e-branch for T1).

    Sig0(a) = ROTR(a,2) XOR ROTR(a,13) XOR ROTR(a,22)

    This ROTATES the segment boundaries. If a had a break at position k,
    ROTR(a,2) has a break at position (k-2) mod 32.

    XOR of three rotated copies: breaks survive only if they coincide
    in all three rotations. A break at position k in a survives in Sig0(a)
    only if k-2, k-13, AND k-22 are all break positions (mod 32).

    This is RESTRICTIVE. Most breaks will NOT survive.

    Let's measure: how many of the segment boundaries in a_new
    are INHERITED from the segment boundaries of a (via Sig0)?
    """
    print("\n" + "=" * 80)
    print("KEY TEST: Do SCA segments survive through Sig0 rotation?")
    print("=" * 80)

    N = 500
    survival_rates = []

    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        states, W = sha256_round_trace(M)

        for r in range(8, 60):
            a = states[r][0]
            b = states[r][1]
            c = states[r][2]

            # Segments of the addition Sig0(a) + Maj(a,b,c)
            s0 = Sig0(a)
            maj = Maj(a, b, c)
            sca_T2 = SCAWord.from_addition(s0, maj)

            # Segments of a itself (from the addition that created it)
            a_prev_r = r - 1
            if a_prev_r >= 0:
                ap = states[a_prev_r]
                T1p = (ap[7] + Sig1(ap[4]) + Ch(ap[4], ap[5], ap[6]) + K[a_prev_r] + W[a_prev_r]) & MASK
                T2p = (Sig0(ap[0]) + Maj(ap[0], ap[1], ap[2])) & MASK
                sca_a = SCAWord.from_addition(T1p, T2p)

                # Break positions in a
                breaks_a = set()
                for start, length in sca_a.segments:
                    breaks_a.add(start)

                # Break positions in T2 = Sig0(a) + Maj(a,b,c)
                breaks_T2 = set()
                for start, length in sca_T2.segments:
                    breaks_T2.add(start)

                # How many T2 breaks were also breaks in a?
                if breaks_a and breaks_T2:
                    inherited = len(breaks_T2 & breaks_a)
                    rate = inherited / len(breaks_T2)
                    survival_rates.append(rate)

    avg_survival = sum(survival_rates) / len(survival_rates) if survival_rates else 0
    print(f"\n  Fraction of T2 segment breaks that were also breaks in a:")
    print(f"  Average: {avg_survival:.4f}")

    # What about random coincidence? If breaks are random at ~25% density:
    # P(random match) = 0.25
    random_expected = 0.25
    print(f"  Random expectation: {random_expected:.4f}")
    print(f"  → {'INHERITED (real structure!)' if avg_survival > random_expected * 1.5 else 'RANDOM (no inheritance)'}")


def experiment_sca_composition():
    """
    The real question: Can we compose SCA across rounds?

    For a FULL preimage attack, we'd need:
    1. Start with target H (known)
    2. Work backward: what segment structure could produce H?
    3. Each backward step: segment structure constrains the previous state
    4. After 64 backward steps: constraints on M

    If each step preserves O(s) segment information (s ≈ 8 segments),
    total information = 64 × log2(configurations per segment).

    Let's measure: how many distinct segment PATTERNS appear?
    """
    print("\n" + "=" * 80)
    print("SCA COMPOSITION: Distinct segment patterns across rounds")
    print("=" * 80)

    N = 1000

    # For each round, collect the segment pattern (just lengths, not positions)
    pattern_counts = {}
    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        states, W = sha256_round_trace(M)

        for r in [0, 4, 16, 32, 63]:
            a = states[r+1][0]
            # Get segments from the addition that produced a
            ap = states[r]
            T1 = (ap[7] + Sig1(ap[4]) + Ch(ap[4], ap[5], ap[6]) + K[r] + W[r]) & MASK
            T2 = (Sig0(ap[0]) + Maj(ap[0], ap[1], ap[2])) & MASK
            sca = SCAWord.from_addition(T1, T2)

            pattern = tuple(sca.segment_lengths())
            key = (r, pattern)
            pattern_counts[key] = pattern_counts.get(key, 0) + 1

    # How many unique patterns per round?
    for r in [0, 4, 16, 32, 63]:
        patterns_at_r = {k: v for k, v in pattern_counts.items() if k[0] == r}
        n_unique = len(patterns_at_r)
        total = sum(patterns_at_r.values())
        # Most common pattern
        if patterns_at_r:
            most_common = max(patterns_at_r.items(), key=lambda x: x[1])
            print(f"  Round {r:>2}: {n_unique} unique segment patterns out of {total}")
            print(f"           Most common: {most_common[0][1]} ({most_common[1]}/{total} = {most_common[1]/total:.1%})")


def experiment_sca_information():
    """
    KEY QUESTION: How much information do segments carry about M?

    If segments carry information → they could guide a search.
    If segments are independent of M → they're useless.

    Test: for two different M with same segment pattern at round 64,
    are they "closer" than random M pairs?
    """
    print("\n" + "=" * 80)
    print("SCA INFORMATION: Do segment patterns carry info about M?")
    print("=" * 80)

    N = 5000
    pattern_to_messages = {}

    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        states, W = sha256_round_trace(M)

        # Segment pattern at final round
        ap = states[63]
        T1 = (ap[7] + Sig1(ap[4]) + Ch(ap[4], ap[5], ap[6]) + K[63] + W[63]) & MASK
        T2 = (Sig0(ap[0]) + Maj(ap[0], ap[1], ap[2])) & MASK
        sca = SCAWord.from_addition(T1, T2)
        pattern = tuple(sca.segment_lengths())

        if pattern not in pattern_to_messages:
            pattern_to_messages[pattern] = []
        pattern_to_messages[pattern].append(M[0])  # Just track W[0]

    # Find patterns with multiple messages
    multi = {p: ms for p, ms in pattern_to_messages.items() if len(ms) >= 5}
    print(f"  Patterns with ≥5 messages: {len(multi)}")

    if multi:
        # For each such pattern, compute avg |W0_i - W0_j| within group
        within_dists = []
        for pattern, w0s in list(multi.items())[:20]:
            for i in range(min(10, len(w0s))):
                for j in range(i+1, min(10, len(w0s))):
                    d = abs(w0s[i] - w0s[j])
                    if d > (1 << 31):
                        d = (1 << 32) - d
                    within_dists.append(d)

        # Random distances
        all_w0 = []
        for ms in pattern_to_messages.values():
            all_w0.extend(ms[:3])
        random_dists = []
        for i in range(min(3000, len(all_w0))):
            j = (i + 7) % len(all_w0)
            d = abs(all_w0[i] - all_w0[j])
            if d > (1 << 31):
                d = (1 << 32) - d
            random_dists.append(d)

        avg_within = sum(within_dists) / len(within_dists) if within_dists else 0
        avg_random = sum(random_dists) / len(random_dists) if random_dists else 1

        ratio = avg_within / avg_random if avg_random > 0 else 0
        print(f"  Avg distance within same pattern: {avg_within:.0f}")
        print(f"  Avg random distance: {avg_random:.0f}")
        print(f"  Ratio: {ratio:.4f}")
        print(f"  → {'CARRIES INFO' if ratio < 0.9 else 'NO INFO (patterns independent of M)'}")
    else:
        print(f"  Not enough multi-message patterns for analysis")


if __name__ == "__main__":
    experiment_sca_basic()
    experiment_sca_through_rounds()
    experiment_sca_composition()
    experiment_sca_information()
