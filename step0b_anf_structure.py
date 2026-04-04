"""
Step 0b: Is there STRUCTURE inside the 50% ANF?

A random Boolean function on 16 variables has:
  - C(16, k) / 2 monomials of degree k on average
  - Total: 2^16 / 2 = 32768 monomials

We compare mini-SHA's ANF with the random expectation at each degree.
If they match → ANF is truly random (no exploitable structure).
If they differ → there's algebraic structure surviving the rounds.

Also: check WHICH VARIABLES appear most/least in monomials.
A random function uses all variables equally.
SHA has structure (shift register, schedule) — does it survive?
"""

import numpy as np
from math import comb
from collections import Counter
from step0_exact_algebra import (
    mini_sha, compute_truth_table, mobius_transform,
    anf_to_monomials, N, MASK
)

def random_anf_expected(n_vars):
    """Expected number of monomials at each degree for random Boolean function."""
    return {k: comb(n_vars, k) / 2.0 for k in range(n_vars + 1)}

def variable_participation(monomials, n_vars):
    """How often each variable appears in monomials."""
    counts = [0] * n_vars
    for m in monomials:
        for v in m:
            counts[v] += 1
    return counts

def pair_participation(monomials, n_vars):
    """How often each PAIR of variables appears together in a monomial."""
    pair_counts = {}
    for m in monomials:
        for i in range(len(m)):
            for j in range(i+1, len(m)):
                pair = (m[i], m[j])
                pair_counts[pair] = pair_counts.get(pair, 0) + 1
    return pair_counts

def var_name(idx):
    word = idx // N
    bit = idx % N
    return f"W{word}b{bit}"

def analyze_structure():
    n_msg = 4
    n_input = n_msg * N  # 16
    expected = random_anf_expected(n_input)

    for R in [1, 3, 5, 6, 8]:
        print(f"\n{'='*70}")
        print(f"R = {R} rounds")
        print(f"{'='*70}")

        tt = compute_truth_table(R, n_msg)

        # Analyze output bit a[0] (representative)
        for out_idx, out_name in [(0, "a[0]"), (N-1, "a[3]"), (N, "e[0]"), (2*N-1, "e[3]")]:
            tt_bit = tt[:, out_idx].copy()
            anf = mobius_transform(tt_bit, n_input)
            monomials = anf_to_monomials(anf, n_input)

            # Degree distribution: observed vs expected random
            degree_dist = Counter(len(m) for m in monomials)
            total_mono = len(monomials)

            print(f"\n  {out_name}: {total_mono} monomials")
            print(f"  {'Deg':>4} {'Obs':>7} {'Expected':>9} {'Ratio':>7} {'Excess':>8}")
            print(f"  {'-'*40}")

            total_excess = 0
            for k in range(n_input + 1):
                obs = degree_dist.get(k, 0)
                exp = expected.get(k, 0)
                if exp > 0:
                    ratio = obs / exp
                    excess = obs - exp
                else:
                    ratio = float('inf') if obs > 0 else 0
                    excess = obs
                if obs > 0 or exp > 1:
                    print(f"  {k:>4} {obs:>7} {exp:>9.1f} {ratio:>7.3f} {excess:>+8.1f}")
                total_excess += abs(excess)

            print(f"  Total |excess| = {total_excess:.1f}")

            # Variable participation
            vp = variable_participation(monomials, n_input)
            expected_vp = total_mono * (n_input // 2) / n_input if total_mono > 0 else 0
            # For random ANF: each var appears in ~half the monomials
            # Specifically: E[participation] = total_monomials * avg_degree / n_vars

            if total_mono > 10:
                avg_deg = sum(len(m) for m in monomials) / total_mono
                expected_per_var = total_mono * avg_deg / n_input

                vp_sorted = sorted(enumerate(vp), key=lambda x: -x[1])
                print(f"\n  Variable participation (expected ≈ {expected_per_var:.0f} per var):")
                print(f"  {'Var':>8} {'Count':>7} {'Ratio':>7}")
                for idx, count in vp_sorted:
                    print(f"  {var_name(idx):>8} {count:>7} {count/expected_per_var:>7.3f}")

    # ============================================================
    # Key test: correlation between ANF coefficients at R=8
    # ============================================================
    print(f"\n{'='*70}")
    print(f"CORRELATION TEST: Are ANF coefficients of different output bits independent?")
    print(f"{'='*70}")

    R = 8
    tt = compute_truth_table(R, n_msg)

    # Get ANF for all 8 output bits
    anfs = []
    for out_idx in range(2 * N):
        tt_bit = tt[:, out_idx].copy()
        anf = mobius_transform(tt_bit, n_input)
        anfs.append(anf)

    # Correlation between ANF coefficient vectors
    # For random independent functions: correlation ≈ 0
    print(f"\n  Correlation matrix (ANF coeff vectors):")
    n_total = 1 << n_input
    corr_matrix = np.zeros((2*N, 2*N))
    for i in range(2*N):
        for j in range(2*N):
            # Correlation = (number of agreeing positions - disagreeing) / total
            agree = np.sum(anfs[i] == anfs[j])
            corr_matrix[i][j] = (2*agree - n_total) / n_total

    labels = [f"a{i}" for i in range(N)] + [f"e{i}" for i in range(N)]
    print(f"  {'':>5}", end="")
    for l in labels:
        print(f" {l:>6}", end="")
    print()
    for i, li in enumerate(labels):
        print(f"  {li:>5}", end="")
        for j in range(len(labels)):
            v = corr_matrix[i][j]
            if i == j:
                print(f"  {'1.00':>5}", end="")
            elif abs(v) > 0.02:
                print(f" {v:>6.3f}", end="")
            else:
                print(f"   {'~0':>4}", end="")
        print()

    # ============================================================
    # Test: Shared monomials between output bits
    # ============================================================
    print(f"\n  Shared monomial analysis (R={R}):")
    for i in range(2*N):
        mono_i = set(np.where(anfs[i] == 1)[0])
        for j in range(i+1, 2*N):
            mono_j = set(np.where(anfs[j] == 1)[0])
            shared = len(mono_i & mono_j)
            expected_shared = len(mono_i) * len(mono_j) / n_total
            if abs(shared - expected_shared) > 100:
                print(f"  {labels[i]}-{labels[j]}: shared={shared}, expected={expected_shared:.0f}, excess={shared-expected_shared:+.0f}")

    # Expected shared for two random functions: p1*p2*N where p1,p2 ≈ 0.5
    print(f"\n  (Only pairs with |excess| > 100 shown)")
    print(f"  For random: expected shared ≈ {n_total/4:.0f} (= N/4)")

    # ============================================================
    # KEY: Monomial pattern by message word
    # ============================================================
    print(f"\n{'='*70}")
    print(f"WORD-LEVEL STRUCTURE (R={R})")
    print(f"{'='*70}")

    # For each monomial, classify by which message words it involves
    tt_bit = tt[:, 0].copy()  # a[0]
    anf = mobius_transform(tt_bit, n_input)
    monomials = anf_to_monomials(anf, n_input)

    word_patterns = Counter()
    for m in monomials:
        words_involved = frozenset(v // N for v in m)
        word_patterns[tuple(sorted(words_involved))] += 1

    print(f"\n  a[0] monomials by message-word pattern:")
    print(f"  {'Pattern':>20} {'Count':>7} {'Frac':>7}")
    for pattern, count in sorted(word_patterns.items(), key=lambda x: -x[1])[:20]:
        wnames = ",".join(f"W{w}" for w in pattern)
        if not wnames:
            wnames = "(const)"
        print(f"  {wnames:>20} {count:>7} {count/len(monomials)*100:>6.1f}%")


if __name__ == "__main__":
    analyze_structure()
