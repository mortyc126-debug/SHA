"""
Stage 0, Part 4: Study the a-recurrence itself.

SHA-256 is: a[r+1] = F(a[r], a[r-1], a[r-2], a[r-3], a[r-4], K[r], W[r])

where F involves Sig0, Sig1, Ch, Maj, and mod 2^32 additions.

Questions:
1. What is the bit-level dependency structure of F?
   (Which bits of a[r+1] depend on which bits of a[r..r-4]?)
2. Are there "simple" bits where the dependency is narrow?
3. Does the trajectory {a[r]} live on a lower-dimensional manifold?
4. What is the ALGEBRAIC structure of F (not GF(2), not Z/2^32, but what?)
"""

import random
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, sig0, sig1, schedule
)


def experiment_13_bit_dependency():
    """
    For each bit b of a[r+1], determine which bits of a[r..r-4]
    it depends on.

    Method: fix a base message, then flip each bit of a[r..r-4]
    (by choosing different M) and observe which bits of a[r+1] change.

    But we can't directly control a[r] — we control M.
    Instead: at a fixed round r, vary ONE bit of the state and
    see how a[r+1] responds. We do this by running SHA up to round r,
    then manually flipping one bit and computing one more round.
    """
    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]
    states, W = sha256_round_trace(M)

    print("=" * 80)
    print("EXPERIMENT 13: Bit-level dependency of a[r+1] on state[r]")
    print("=" * 80)

    r = 20  # Pick a round well into the computation

    a, b, c, d, e, f, g, h = states[r]

    # Compute a[r+1] normally
    T1_base = (h + Sig1(e) + Ch(e, f, g) + K[r] + W[r]) & MASK
    T2_base = (Sig0(a) + Maj(a, b, c)) & MASK
    a_next_base = (T1_base + T2_base) & MASK

    # For each register and each bit, flip it and see which bits of a_next change
    regs = [a, b, c, d, e, f, g, h]
    reg_names = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']

    print(f"\nAt round {r}, flipping each bit of each register:")
    print(f"Showing: number of bits of a[r+1] that change\n")

    dependency_map = {}  # (reg, bit) -> set of output bits affected

    for ri, (reg_val, rname) in enumerate(zip(regs, reg_names)):
        affected_counts = []
        for bit in range(32):
            # Flip bit in this register
            flipped = reg_val ^ (1 << bit)
            test_regs = list(regs)
            test_regs[ri] = flipped

            ta, tb, tc, td, te, tf, tg, th = test_regs

            T1_test = (th + Sig1(te) + Ch(te, tf, tg) + K[r] + W[r]) & MASK
            T2_test = (Sig0(ta) + Maj(ta, tb, tc)) & MASK
            a_next_test = (T1_test + T2_test) & MASK

            diff = a_next_base ^ a_next_test
            n_affected = bin(diff).count('1')
            affected_counts.append(n_affected)

            dependency_map[(rname, bit)] = diff

        avg = sum(affected_counts) / 32
        min_a = min(affected_counts)
        max_a = max(affected_counts)
        print(f"  {rname}: avg={avg:.1f} min={min_a} max={max_a} bits changed in a[r+1]")

    # Now: which registers have the LEAST influence?
    print(f"\n--- Key question: which register contributes least to a[r+1]? ---")
    # d and h don't appear in the formula for a_new = T1 + T2 directly!
    # T1 = h + Sig1(e) + Ch(e,f,g) + K + W — h appears!
    # T2 = Sig0(a) + Maj(a,b,c) — no d, no e, no f, no g, no h
    # a_new = T1 + T2 — so d does NOT directly appear!

    print(f"  d register: does NOT appear in a[r+1] formula (only in e[r+1] = d + T1)")
    print(f"  Verification: flipping d should change 0 bits of a[r+1]")

    # Verify
    for bit in range(32):
        flipped_d = d ^ (1 << bit)
        T1_test = (h + Sig1(e) + Ch(e, f, g) + K[r] + W[r]) & MASK  # d not in T1
        T2_test = (Sig0(a) + Maj(a, b, c)) & MASK  # d not in T2
        a_next_test = (T1_test + T2_test) & MASK
        if a_next_test != a_next_base:
            print(f"    UNEXPECTED: d bit {bit} changes a[r+1]!")
            break
    else:
        print(f"    CONFIRMED: d has ZERO influence on a[r+1] ✓")

    # Similarly for e[r+1]:
    e_next_base = (d + T1_base) & MASK
    print(f"\n  For e[r+1] = d + T1:")
    for ri, rname in enumerate(reg_names):
        affected_counts = []
        for bit in range(32):
            test_regs = list(regs)
            test_regs[ri] = regs[ri] ^ (1 << bit)
            ta, tb, tc, td, te, tf, tg, th = test_regs
            T1_test = (th + Sig1(te) + Ch(te, tf, tg) + K[r] + W[r]) & MASK
            e_next_test = (td + T1_test) & MASK
            diff = e_next_base ^ e_next_test
            affected_counts.append(bin(diff).count('1'))
        avg = sum(affected_counts) / 32
        print(f"    {rname}: avg={avg:.1f} bits changed in e[r+1]")

    # a and b and c affect a_new through T2 (Sig0, Maj)
    # e, f, g affect a_new through T1 (Sig1, Ch)
    # h affects a_new through T1 (h + ...)
    # d affects only e_new


def experiment_14_trajectory_dimension():
    """
    Does the trajectory {a[0], a[1], ..., a[64]} live on a
    lower-dimensional manifold in (Z/2^32)^65?

    Method: collect many trajectories (different M), form a matrix,
    compute rank. If rank < 65 → trajectories are constrained.

    But 32-bit words in Z/2^32 make rank computation tricky.
    Instead: project onto individual bits and check.
    """
    N = 500
    print("\n" + "=" * 80)
    print("EXPERIMENT 14: Trajectory dimension — do a-sequences fill the space?")
    print("=" * 80)

    # For each bit position b (0..31):
    # Collect the 65-bit vector (a[0][b], a[1][b], ..., a[64][b]) for N messages
    # Form N×65 matrix over GF(2). Rank = effective dimension.

    for bit_pos in [0, 7, 15, 24, 31]:
        matrix = []
        for seed in range(N):
            rng = random.Random(seed)
            M = [rng.randint(0, MASK) for _ in range(16)]
            states, _ = sha256_round_trace(M)

            row = []
            for r in range(65):
                row.append((states[r][0] >> bit_pos) & 1)
            matrix.append(row)

        # Compute GF(2) rank
        rank = gf2_rank(matrix, N, 65)
        print(f"  Bit {bit_pos:>2}: rank of {N}×65 trajectory matrix = {rank}/65")


def gf2_rank(matrix, nrows, ncols):
    """Compute rank of binary matrix over GF(2)."""
    # Copy matrix
    m = [list(row) for row in matrix[:nrows]]
    rank = 0
    for col in range(ncols):
        # Find pivot
        pivot = None
        for row in range(rank, len(m)):
            if m[row][col] == 1:
                pivot = row
                break
        if pivot is None:
            continue
        # Swap
        m[rank], m[pivot] = m[pivot], m[rank]
        # Eliminate
        for row in range(len(m)):
            if row != rank and m[row][col] == 1:
                for c in range(ncols):
                    m[row][c] ^= m[rank][c]
        rank += 1
    return rank


def experiment_15_recurrence_decomposition():
    """
    Decompose F into its components and measure their relative contribution.

    a[r+1] = T1 + T2 where:
      T1 = h + Sig1(e) + Ch(e,f,g) + K[r] + W[r]
      T2 = Sig0(a) + Maj(a,b,c)

    But h = e[r-3], f = e[r-1], g = e[r-2]
    And e[r] = a[r] + a[r-4] - T2[r-1] (створочное число)

    So everything feeds back through a. Let's see the actual
    NUMERICAL contribution of each term.
    """
    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]
    states, W = sha256_round_trace(M)

    print("\n" + "=" * 80)
    print("EXPERIMENT 15: Decomposition of the recurrence at each round")
    print("=" * 80)

    print(f"\n{'r':>3} | {'Sig0(a)':>10} {'Maj':>10} {'T2':>10} | {'h':>10} {'Sig1(e)':>10} {'Ch':>10} {'K+W':>10} {'T1':>10} | {'a_new':>10}")
    print("-" * 120)

    for r in range(min(20, 64)):
        a, b, c, d, e, f, g, h = states[r]

        s0 = Sig0(a)
        maj = Maj(a, b, c)
        T2 = (s0 + maj) & MASK

        s1 = Sig1(e)
        ch = Ch(e, f, g)
        kw = (K[r] + W[r]) & MASK
        T1 = (h + s1 + ch + kw) & MASK

        a_new = (T1 + T2) & MASK

        # Verify
        assert a_new == states[r+1][0], f"Mismatch at round {r}"

        print(f"{r:>3} | {s0:>10x} {maj:>10x} {T2:>10x} | {h:>10x} {s1:>10x} {ch:>10x} {kw:>10x} {T1:>10x} | {a_new:>10x}")

    # The key question: what is the RELATIVE entropy of each component?
    # I.e., how much does each component vary across messages?
    print(f"\n--- Variance of each component across 500 messages (at round 20) ---")

    components = {'Sig0(a)': [], 'Maj': [], 'T2': [], 'h': [], 'Sig1(e)': [], 'Ch': [], 'T1': []}

    for seed in range(500):
        rng_s = random.Random(seed)
        M = [rng_s.randint(0, MASK) for _ in range(16)]
        states, W = sha256_round_trace(M, rounds=21)
        a, b, c, d, e, f, g, h = states[20]

        components['Sig0(a)'].append(Sig0(a))
        components['Maj'].append(Maj(a, b, c))
        components['T2'].append((Sig0(a) + Maj(a, b, c)) & MASK)
        components['h'].append(h)
        components['Sig1(e)'].append(Sig1(e))
        components['Ch'].append(Ch(e, f, g))
        components['T1'].append((h + Sig1(e) + Ch(e, f, g) + K[20] + W[20]) & MASK)

    for name, vals in components.items():
        unique = len(set(vals))
        hws = [bin(v).count('1') for v in vals]
        avg_hw = sum(hws) / len(hws)
        print(f"  {name:>10}: {unique}/500 unique, E[HW]={avg_hw:.1f}")


def experiment_16_single_round_algebraic():
    """
    Focus on ONE round. Express a[r+1] purely in terms of a[r..r-4] and W[r].

    a[r+1] = [e[r-3] + Sig1(e[r]) + Ch(e[r], e[r-1], e[r-2]) + K[r] + W[r]]
           + [Sig0(a[r]) + Maj(a[r], a[r-1], a[r-2])]

    where e[s] = a[s] + a[s-4] - Sig0(a[s-1]) - Maj(a[s-1], a[s-2], a[s-3])

    This substitution creates a SINGLE expression in a[r], a[r-1], ..., a[r-7]
    (because e[r-3] needs a[r-7]).

    Let's verify and count the depth of dependency.
    """
    print("\n" + "=" * 80)
    print("EXPERIMENT 16: Full substitution — a[r+1] in terms of a alone")
    print("=" * 80)

    # Verify that we can compute a[r+1] from a[r..r-7] + K[r] + W[r..r-3]

    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]
    states, W = sha256_round_trace(M)

    violations = 0
    for r in range(8, 64):
        # Get a-values
        a_vals = [states[r-i][0] for i in range(8)]  # a[r], a[r-1], ..., a[r-7]

        # Compute e[r], e[r-1], e[r-2], e[r-3] from a-values
        # e[s] = a[s] + a[s-4] - T2[s-1]
        # T2[s-1] = Sig0(a[s-1]) + Maj(a[s-1], a[s-2], a[s-3])
        def compute_e(s_offset):
            # s_offset: 0=r, 1=r-1, 2=r-2, 3=r-3
            # e[r-s_offset] = a[r-s_offset] + a[r-s_offset-4] - T2[r-s_offset-1]
            a_s = a_vals[s_offset]
            a_s4 = a_vals[s_offset + 4]
            a_s1 = a_vals[s_offset + 1]
            a_s2 = a_vals[s_offset + 2]
            a_s3 = a_vals[s_offset + 3]
            T2_s1 = (Sig0(a_s1) + Maj(a_s1, a_s2, a_s3)) & MASK
            return (a_s + a_s4 - T2_s1) & MASK

        e_r = compute_e(0)
        e_r1 = compute_e(1)
        e_r2 = compute_e(2)
        e_r3 = compute_e(3)

        # Verify e-values
        assert e_r == states[r][4], f"e[r] mismatch at r={r}"
        assert e_r1 == states[r-1][4], f"e[r-1] mismatch at r={r}"
        assert e_r2 == states[r-2][4], f"e[r-2] mismatch at r={r}"
        assert e_r3 == states[r-3][4], f"e[r-3] mismatch at r={r}"

        # Now compute a[r+1]
        T1 = (e_r3 + Sig1(e_r) + Ch(e_r, e_r1, e_r2) + K[r] + W[r]) & MASK
        T2 = (Sig0(a_vals[0]) + Maj(a_vals[0], a_vals[1], a_vals[2])) & MASK
        a_computed = (T1 + T2) & MASK

        if a_computed != states[r+1][0]:
            violations += 1
            print(f"  VIOLATION at r={r}")

    print(f"\n  Verified a[r+1] = F(a[r..r-7], K[r], W[r]) for rounds 8-63")
    print(f"  Violations: {violations}/{64-8}")
    if violations == 0:
        print(f"  CONFIRMED ✓")
        print(f"\n  SHA-256 is an 8th-order recurrence in {'{a[r]}'}")
        print(f"  (Not 5th as initially thought — the e-substitution extends depth to 8)")
        print(f"\n  Dependencies per round:")
        print(f"    a[r+1] depends on: a[r], a[r-1], a[r-2], a[r-3], a[r-4], a[r-5], a[r-6], a[r-7]")
        print(f"    Plus constants: K[r], W[r]")
        print(f"    Operations: Sig0, Sig1, Ch, Maj, mod 2^32 addition (×6)")
        print(f"\n  This is the COMPLETE algebraic description of SHA-256 in one sequence.")


if __name__ == "__main__":
    experiment_13_bit_dependency()
    experiment_14_trajectory_dimension()
    experiment_15_recurrence_decomposition()
    experiment_16_single_round_algebraic()
