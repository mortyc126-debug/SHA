"""
Step 15: Conditional Propagation Algebra (CPA)

Core idea: Ch(e,f,g) is degree-2 ONLY when 'e' is unknown.
When e is KNOWN (e.g., δe=0 from collision constraint):
  Ch(e, f, g) with known δe → δCh determined in O(1) (degree 0!)

Strategy: propagate BACKWARD from collision condition (δstate_R=0).
At each round, the known bits from the collision constraint
may RESOLVE Ch/Maj constraints for free.

Every resolved Ch = one degree-2 constraint REMOVED.
Count: how many Ch/Maj constraints survive after propagation?
If few survive → system is solvable faster than birthday.

This is a NEW ALGEBRA because it tracks (value, knowledge_state) pairs
and uses knowledge to reduce degree at each operation.
"""

import numpy as np
from step0_exact_algebra import mini_sha, N, MASK

N_MSG = 4
N_INPUT = N * N_MSG
N_TOTAL = 1 << N_INPUT

IV = [0x6, 0xB, 0x3, 0xA, 0x5, 0x9, 0x1, 0xF]
K_CONST = [0x4, 0x2, 0xB, 0x7, 0xA, 0x3, 0xE, 0x5,
           0x9, 0x1, 0xD, 0x6, 0x0, 0x8, 0xC, 0xF]


def rotr(x, s):
    return ((x >> s) | (x << (N - s))) & MASK


class CPAState:
    """
    CPA State: tracks which bits of δ-state are KNOWN at each round.

    A bit can be:
      - KNOWN (value determined by constraints) → Ch/Maj with this bit = FREE
      - UNKNOWN (free variable) → Ch/Maj with this bit = degree-2 constraint

    The more bits are known, the fewer degree-2 constraints survive.
    """

    def __init__(self):
        # δ-state for 8 registers, N bits each
        # None = unknown, 0/1 = known value
        self.regs = [[None]*N for _ in range(8)]  # a,b,c,d,e,f,g,h
        self.n_constraints = 0  # degree-2 constraints accumulated
        self.n_free_resolves = 0  # Ch/Maj resolved for free (known bits)

    def known_count(self):
        count = 0
        for reg in self.regs:
            for bit in reg:
                if bit is not None:
                    count += 1
        return count

    def unknown_count(self):
        return 8 * N - self.known_count()

    def set_all_zero(self):
        """Set all δ-state bits to 0 (collision condition)."""
        for i in range(8):
            for j in range(N):
                self.regs[i][j] = 0

    def set_all_unknown(self):
        """Set all δ-state bits to unknown."""
        for i in range(8):
            for j in range(N):
                self.regs[i][j] = None

    def copy(self):
        s = CPAState()
        for i in range(8):
            for j in range(N):
                s.regs[i][j] = self.regs[i][j]
        s.n_constraints = self.n_constraints
        s.n_free_resolves = self.n_free_resolves
        return s


def backward_propagate(M_base, R):
    """
    Propagate collision constraint BACKWARD through R rounds.

    Start: δstate[R] = 0 (all 8*N bits known = 0)
    At each round (backward): determine what constraints exist on δstate[r-1]

    For each Ch/Maj: check if the controlling bit is KNOWN.
    If yes → resolve for FREE (no degree-2 constraint).
    If no → add degree-2 constraint.

    Count total degree-2 constraints after full backward propagation.
    """
    print(f"\n{'='*80}")
    print(f"CPA BACKWARD PROPAGATION — R={R}")
    print(f"M = {[hex(w) for w in M_base]}")
    print(f"{'='*80}")

    # Run forward computation to get base state at each round
    a, b, c, d, e, f, g, h = IV[:]
    W = list(M_base) + [0] * max(0, R - N_MSG)

    states = [(a, b, c, d, e, f, g, h)]
    T1_vals = []
    T2_vals = []

    for r in range(R):
        w_r = W[r] if r < len(W) else 0
        k_r = K_CONST[r % len(K_CONST)]
        sig1 = rotr(e, 1) ^ rotr(e, 3) ^ (e >> 1)
        ch_val = (e & f) ^ (~e & g) & MASK
        sig0 = rotr(a, 1) ^ rotr(a, 2) ^ rotr(a, 3)
        maj_val = (a & b) ^ (a & c) ^ (b & c)
        T1 = (h + sig1 + ch_val + k_r + w_r) & MASK
        T2 = (sig0 + maj_val) & MASK
        T1_vals.append(T1)
        T2_vals.append(T2)
        a_new = (T1 + T2) & MASK
        e_new = (d + T1) & MASK
        h, g, f = g, f, e
        e = e_new
        d, c, b = c, b, a
        a = a_new
        states.append((a, b, c, d, e, f, g, h))

    # Start backward propagation
    # State at round R: all δ = 0 (collision)
    cpa = CPAState()
    cpa.set_all_zero()

    total_ch_ops = 0
    total_maj_ops = 0
    free_ch = 0
    free_maj = 0
    constrained_ch = 0
    constrained_maj = 0

    print(f"\n  {'Round':>6} {'Known':>6} {'Unknown':>8} {'Ch_free':>8} {'Ch_cost':>8} "
          f"{'Maj_free':>9} {'Maj_cost':>9} {'Constraints':>12}")

    for r in range(R-1, -1, -1):
        # At round r: we know δstate[r+1].
        # Need to determine what we know about δstate[r].

        # Current known δstate[r+1]
        da_next = cpa.regs[0]  # δa at round r+1
        db_next = cpa.regs[1]  # δb = δa[r] (shift register)
        dc_next = cpa.regs[2]
        dd_next = cpa.regs[3]
        de_next = cpa.regs[4]  # δe at round r+1
        df_next = cpa.regs[5]  # δf = δe[r] (shift register)
        dg_next = cpa.regs[6]
        dh_next = cpa.regs[7]

        # From shift register: δb[r+1] = δa[r], δc[r+1] = δb[r], etc.
        # So: δa[r] = δb[r+1] (KNOWN if δb[r+1] is known)
        # δe[r] = δf[r+1] (KNOWN if δf[r+1] is known)

        da_r = list(db_next)  # δa[r] = δb[r+1]
        db_r = list(dc_next)  # δb[r] = δc[r+1]
        dc_r = list(dd_next)  # δc[r] = δd[r+1]
        de_r = list(df_next)  # δe[r] = δf[r+1]
        df_r = list(dg_next)  # δf[r] = δg[r+1]
        dg_r = list(dh_next)  # δg[r] = δh[r+1]

        # δe[r+1] = δd[r] + δT1[r] (mod 2^n)
        # → δT1[r] = δe[r+1] - δd[r]
        # δa[r+1] = δT1[r] + δT2[r]
        # → δT2[r] = δa[r+1] - δT1[r]

        # δd[r] = δc[r-1] = ... (from shift register of previous round)
        # For now: δd[r] known if δc[r+1] propagated...
        # Actually: δd[r] = δd_next... wait.
        # δd[r+1] = δc[r]. δd[r] = δc[r-1]. But δd at round r is what becomes d_next.
        # Hmm, δd[r] = δc[r] of the PREVIOUS state... this is just the shift register going backward.

        # The shift register backward:
        # state[r] = (a[r], b[r]=a[r-1], c[r]=a[r-2], d[r]=a[r-3],
        #             e[r], f[r]=e[r-1], g[r]=e[r-2], h[r]=e[r-3])

        # So δd[r] = δa[r-3]. For backward prop from r+1 to r:
        # We now know δa[r] from shift register.
        # δd[r] is δa[r-3] which might be known from earlier backward steps.

        # For this analysis: just track which δ-bits are known/unknown
        dd_r = [None] * N  # δd[r] = δa[r-3], may be unknown
        dh_r = [None] * N  # δh[r] = δe[r-3], may be unknown

        # If we've propagated back far enough, δa[r-3] might be known
        # For first 3 backward rounds: δd and δh unknown (not yet propagated)
        # After 3 backward rounds: δd = δa from 3 rounds ahead (known)

        if r + 3 < R:
            # δd[r] = δa[r] shifted 3 rounds... actually from the forward direction,
            # d[r] = c[r-1] = b[r-2] = a[r-3].
            # In backward: we know δa[r] = δb[r+1]. We need δd[r] = δa[r-3].
            # If r >= 3: δa[r-3] was computed in a PREVIOUS backward step.
            # For now: mark as unknown until we track it
            pass

        # Count Ch/Maj operations and whether they're free
        base = states[r]
        base_e = base[4]

        # Ch(e, f, g): controlling bit = e
        for bit in range(N):
            total_ch_ops += 1
            e_bit_known = de_r[bit] is not None
            if e_bit_known:
                free_ch += 1
                cpa.n_free_resolves += 1
            else:
                constrained_ch += 1
                cpa.n_constraints += 1

        # Maj(a, b, c): no single controlling bit, but if 2 of 3 known and equal → free
        for bit in range(N):
            total_maj_ops += 1
            a_known = da_r[bit] is not None
            b_known = db_r[bit] is not None
            c_known = dc_r[bit] is not None

            n_known = sum([a_known, b_known, c_known])
            if n_known >= 2:
                free_maj += 1
                cpa.n_free_resolves += 1
            else:
                constrained_maj += 1
                cpa.n_constraints += 1

        # Update CPA state for round r
        new_cpa = CPAState()
        new_cpa.regs[0] = da_r  # a[r]
        new_cpa.regs[1] = db_r
        new_cpa.regs[2] = dc_r
        new_cpa.regs[3] = dd_r
        new_cpa.regs[4] = de_r
        new_cpa.regs[5] = df_r
        new_cpa.regs[6] = dg_r
        new_cpa.regs[7] = dh_r
        new_cpa.n_constraints = cpa.n_constraints
        new_cpa.n_free_resolves = cpa.n_free_resolves

        known = new_cpa.known_count()
        unknown = new_cpa.unknown_count()

        print(f"  R={r:>4} {known:>6} {unknown:>8} {free_ch:>8} {constrained_ch:>8} "
              f"{free_maj:>9} {constrained_maj:>9} {cpa.n_constraints:>12}")

        cpa = new_cpa

    # Final summary
    print(f"\n{'='*80}")
    print(f"  CPA SUMMARY")
    print(f"{'='*80}")
    print(f"  Total Ch operations:  {total_ch_ops}")
    print(f"  Free Ch (known e):    {free_ch} ({free_ch/max(total_ch_ops,1)*100:.1f}%)")
    print(f"  Constrained Ch:       {constrained_ch} ({constrained_ch/max(total_ch_ops,1)*100:.1f}%)")
    print(f"  Total Maj operations: {total_maj_ops}")
    print(f"  Free Maj (≥2 known):  {free_maj} ({free_maj/max(total_maj_ops,1)*100:.1f}%)")
    print(f"  Constrained Maj:      {constrained_maj} ({constrained_maj/max(total_maj_ops,1)*100:.1f}%)")
    print(f"  Total degree-2 constraints: {cpa.n_constraints}")
    print(f"  Total free resolves:  {cpa.n_free_resolves}")
    print(f"  Constraint ratio: {cpa.n_constraints}/{total_ch_ops+total_maj_ops} "
          f"= {cpa.n_constraints/(total_ch_ops+total_maj_ops)*100:.1f}%")

    # Cost estimate
    # Each surviving degree-2 constraint contributes ~1 bit of information
    # Total search space: 2^(n_constraints) if constraints are independent
    # Birthday cost: 2^(N) = 2^4 = 16 for mini-SHA
    print(f"\n  COST ESTIMATE:")
    print(f"  CPA search space: ~2^{cpa.n_constraints} (if constraints independent)")
    print(f"  Birthday cost:    2^{N} = {2**N}")
    print(f"  Brute force:      2^{N_INPUT} = {N_TOTAL}")
    if cpa.n_constraints < N:
        print(f"  ★ CPA potentially BEATS birthday! ({cpa.n_constraints} < {N})")
    elif cpa.n_constraints < N_INPUT:
        print(f"  CPA between birthday and brute ({N} < {cpa.n_constraints} < {N_INPUT})")
    else:
        print(f"  CPA = brute force ({cpa.n_constraints} ≥ {N_INPUT})")

    return cpa.n_constraints, free_ch + free_maj


def main():
    import random
    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(N_MSG)]

    for R in [3, 4, 6, 8]:
        backward_propagate(M, R)


if __name__ == "__main__":
    main()
