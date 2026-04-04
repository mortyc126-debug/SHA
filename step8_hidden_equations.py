"""
Step 8: Hidden Equations — the intermediate rounds as constraints

Step 7 showed: 8 final equations in 16 unknowns of degree 15 = hopeless.

But SHA computes through INTERMEDIATE STATES. Each round adds:
  - Carry bits (new variables)
  - MAJ/Ch equations (degree-2 constraints on carries)
  - Round function equations (linking rounds)

If we DON'T eliminate intermediates, the system stays degree-2
but has MORE equations. Question: does the ratio improve?

For 2n→n function (collision):
  Final only: n equations, 2n variables → ratio 0.5 → birthday
  With intermediates: ~R×K equations, ~R×K variables → ratio ≈ 1.0?

Let's count EXACTLY for our mini-SHA and see if intermediates help.
"""

import numpy as np
from step0_exact_algebra import N, MASK

N_MSG = 4
N_INPUT = N * N_MSG


def count_system_size(R):
    """
    Count variables and equations in the FULL GF(2) system for mini-SHA collision.

    Variables:
      - 16 δM bits (message difference)
      - Per round: carry bits from each addition

    Equations:
      - Per round: defining each output bit as function of inputs
      - Final: δstate[R] = 0 (8 equations)
    """

    # Per round operations for mini-SHA:
    # T1 = h + sig1(e) + ch(e,f,g) + K + W
    #   = ((h + sig1(e)) + ch(e,f,g)) + K) + W)
    #   4 additions, each with N-1 carry bits
    #   Plus Ch: N AND operations (degree 2)
    #
    # T2 = sig0(a) + maj(a,b,c)
    #   1 addition, N-1 carry bits
    #   Plus Maj: N AND operations (degree 2)
    #
    # a_new = T1 + T2
    #   1 addition, N-1 carry bits
    #
    # e_new = d + T1
    #   1 addition, N-1 carry bits
    #
    # Total additions per round: 4 + 1 + 1 + 1 = 7
    # Carry bits per addition: N-1 = 3
    # Total carry bits per round: 7 * 3 = 21

    adds_per_round = 7  # h+sig1, +ch, +K, +W, sig0+maj, T1+T2, d+T1
    carry_per_add = N - 1  # 3 for n=4
    carry_per_round = adds_per_round * carry_per_add

    # New state bits per round: a_new (N bits) + e_new (N bits)
    # b,c,d,f,g,h are just shifted copies — no new info
    state_new_per_round = 2 * N  # 8

    # Equations per round:
    # Each output bit of each addition: N equations per addition
    # But some are linear (XOR) — don't add variables
    # The degree-2 equations come from:
    #   - Carry bits: carry[k] = MAJ(a[k-1], b[k-1], carry[k-1]) → 1 eq per carry
    #   - Ch bits: ch[k] = e[k]·f[k] ⊕ NOT(e[k])·g[k] → 1 eq per bit
    #   - Maj bits: maj[k] = a[k]·b[k] ⊕ a[k]·c[k] ⊕ b[k]·c[k] → 1 eq per bit
    #   - Sum bits: sum[k] = a[k] ⊕ b[k] ⊕ carry[k] → linear (free)

    # Degree-2 equations per round:
    carry_equations = carry_per_round  # 21
    ch_equations = N  # 4
    maj_equations = N  # 4
    deg2_equations_per_round = carry_equations + ch_equations + maj_equations  # 29

    # Linear equations per round (sum bits, rotations, shifts):
    sum_equations = adds_per_round * N  # 7 * 4 = 28
    rotation_equations = 2 * N  # sig0, sig1 (each N bits, but linear)
    linear_per_round = sum_equations + rotation_equations  # 36

    # Variables per round:
    vars_per_round = carry_per_round + state_new_per_round  # 21 + 8 = 29
    # (Ch and Maj outputs can be expressed directly — not separate vars)

    # Total for R rounds:
    msg_vars = N_INPUT  # 16
    intermediate_vars = R * vars_per_round
    total_vars = msg_vars + intermediate_vars

    total_deg2_eqs = R * deg2_equations_per_round
    total_linear_eqs = R * linear_per_round
    output_eqs = 2 * N  # δa[R]=0, δe[R]=0 → 8 equations
    total_eqs = total_deg2_eqs + total_linear_eqs + output_eqs

    return {
        'msg_vars': msg_vars,
        'intermediate_vars': intermediate_vars,
        'total_vars': total_vars,
        'deg2_eqs': total_deg2_eqs,
        'linear_eqs': total_linear_eqs,
        'output_eqs': output_eqs,
        'total_eqs': total_eqs,
        'vars_per_round': vars_per_round,
        'eqs_per_round': deg2_equations_per_round + linear_per_round,
    }


def analyze():
    print(f"HIDDEN EQUATIONS ANALYSIS")
    print(f"Mini-SHA: n={N}, msg_words={N_MSG}, input_bits={N_INPUT}")
    print(f"{'='*80}\n")

    # Compare: final-only vs full system
    print(f"{'COMPARISON: Final-only vs Full system':^80}")
    print(f"{'='*80}")
    print(f"\n  {'Round':>6} │ {'Final only':^25} │ {'Full system':^35}")
    print(f"  {'':>6} │ {'Eqs':>6} {'Vars':>6} {'Ratio':>7} │ {'Eqs':>8} {'Vars':>8} {'Ratio':>7} {'Deg2':>6}")
    print(f"  {'─'*6}─┼─{'─'*25}─┼─{'─'*35}")

    for R in range(1, 9):
        info = count_system_size(R)

        # Final-only system
        final_eqs = 2 * N  # 8
        final_vars = N_INPUT  # 16
        final_ratio = final_eqs / final_vars

        # Full system
        full_eqs = info['total_eqs']
        full_vars = info['total_vars']
        full_ratio = full_eqs / full_vars
        deg2 = info['deg2_eqs']

        print(f"  R={R:>3} │ {final_eqs:>6} {final_vars:>6} {final_ratio:>7.3f} │ "
              f"{full_eqs:>8} {full_vars:>8} {full_ratio:>7.3f} {deg2:>6}")

    # Detailed breakdown
    R = 8
    info = count_system_size(R)
    print(f"\n{'='*80}")
    print(f"DETAILED BREAKDOWN for R={R}")
    print(f"{'='*80}")
    print(f"  Message variables:       {info['msg_vars']}")
    print(f"  Intermediate variables:  {info['intermediate_vars']} ({info['vars_per_round']} per round × {R})")
    print(f"  TOTAL variables:         {info['total_vars']}")
    print(f"")
    print(f"  Degree-2 equations:      {info['deg2_eqs']}")
    print(f"  Linear equations:        {info['linear_eqs']}")
    print(f"  Output equations:        {info['output_eqs']}")
    print(f"  TOTAL equations:         {info['total_eqs']}")
    print(f"")
    print(f"  Ratio (eqs/vars):        {info['total_eqs']/info['total_vars']:.3f}")
    print(f"  Excess equations:        {info['total_eqs'] - info['total_vars']}")

    # After Gaussian elimination on linear equations
    # Linear equations can eliminate many intermediate variables
    # What remains: degree-2 equations on remaining variables
    remaining_vars = info['total_vars'] - info['linear_eqs']  # optimistic
    remaining_eqs = info['deg2_eqs'] + info['output_eqs']

    print(f"\n  After eliminating linear equations (optimistic):")
    print(f"  Remaining variables:     {remaining_vars}")
    print(f"  Remaining equations:     {remaining_eqs} (all degree 2)")
    print(f"  Ratio:                   {remaining_eqs/max(remaining_vars,1):.3f}")

    # For SHA-256
    print(f"\n{'='*80}")
    print(f"SCALING TO SHA-256 (n=32, R=64)")
    print(f"{'='*80}")

    n_sha = 32
    R_sha = 64

    # SHA-256 has the same structure but n=32
    adds_sha = 7
    carry_sha = n_sha - 1  # 31
    carry_per_round_sha = adds_sha * carry_sha  # 217
    state_new_sha = 2 * n_sha  # 64

    vars_per_round_sha = carry_per_round_sha + state_new_sha  # 281
    deg2_per_round_sha = carry_per_round_sha + 2 * n_sha  # 217 + 64 = 281
    linear_per_round_sha = adds_sha * n_sha + 2 * n_sha  # 224 + 64 = 288

    msg_vars_sha = 512
    total_vars_sha = msg_vars_sha + R_sha * vars_per_round_sha
    total_eqs_sha = R_sha * (deg2_per_round_sha + linear_per_round_sha) + 256

    print(f"  Message variables:       {msg_vars_sha}")
    print(f"  Variables per round:     {vars_per_round_sha}")
    print(f"  Total variables:         {total_vars_sha}")
    print(f"  Total equations:         {total_eqs_sha}")
    print(f"  Ratio:                   {total_eqs_sha/total_vars_sha:.3f}")
    print(f"  Excess:                  {total_eqs_sha - total_vars_sha}")

    remaining_sha = total_vars_sha - R_sha * linear_per_round_sha
    remaining_eqs_sha = R_sha * deg2_per_round_sha + 256
    print(f"\n  After eliminating linear equations:")
    print(f"  Remaining variables:     {remaining_sha}")
    print(f"  Remaining deg-2 eqs:     {remaining_eqs_sha}")
    print(f"  Ratio:                   {remaining_eqs_sha/max(remaining_sha,1):.3f}")

    # The critical number: what XL degree is needed?
    print(f"\n{'='*80}")
    print(f"XL ANALYSIS: What degree elevation is needed?")
    print(f"{'='*80}")

    print(f"\n  For n variables, m degree-2 equations:")
    print(f"  XL at degree D succeeds when m × C(n, D-2) > C(n, D)")
    print(f"  i.e., when m > C(n,D) / C(n,D-2) = (D-1)×D / ((n-D+1)(n-D+2))")
    print(f"")

    # Mini-SHA after linear elimination
    n_mini = max(remaining_vars, 1)
    m_mini = remaining_eqs
    print(f"  Mini-SHA (n={n_mini}, m={m_mini}):")
    for D in range(3, min(n_mini+1, 10)):
        from math import comb
        xl_eqs = m_mini * comb(n_mini, D-2) if D-2 <= n_mini else 0
        xl_monos = comb(n_mini, D) if D <= n_mini else 0
        if xl_monos > 0:
            success = "YES" if xl_eqs > xl_monos else "no"
            print(f"    D={D}: XL equations={xl_eqs}, monomials={xl_monos}, ratio={xl_eqs/xl_monos:.2f} → {success}")

    # SHA-256
    n_real = max(remaining_sha, 1)
    m_real = remaining_eqs_sha
    print(f"\n  SHA-256 (n={n_real}, m={m_real}):")
    for D in range(3, 8):
        from math import comb
        xl_eqs = m_real * comb(n_real, D-2) if D-2 <= n_real else 0
        xl_monos = comb(n_real, D) if D <= n_real else 0
        if xl_monos > 0 and xl_eqs > 0:
            ratio = xl_eqs / xl_monos
            success = "YES" if ratio > 1 else "no"
            print(f"    D={D}: ratio={ratio:.6f} → {success}")
            if ratio > 1:
                print(f"    ★ XL succeeds at degree {D}!")
                print(f"    Cost: O(C(n,D)^ω) = O({comb(n_real,D)}^2.37)")
                break


if __name__ == "__main__":
    analyze()
