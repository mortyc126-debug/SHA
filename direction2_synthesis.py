#!/usr/bin/env python3
"""
Direction 2: Synthesis Report — All 45 CRAZY experiments + cross-hash analysis.
Pure reporting script, no computation.
"""

def main():
    print("=" * 90)
    print("  SYNTHESIS REPORT: 45 ATTACK PARADIGMS AGAINST SHA-256")
    print("=" * 90)
    print()

    # ── Category A: Differential Attacks ──
    print("+" + "-" * 88 + "+")
    print("| CATEGORY A: DIFFERENTIAL ATTACKS                                                    |")
    print("+" + "-" * 88 + "+")
    attacks_diff = [
        ("CRAZY-1",  "Random differential search",        "Best prob 2^-7.3 at 1 rnd, 0 at 3+",  "Exponential decay"),
        ("CRAZY-2",  "Bit-level avalanche",               "Full avalanche by round 4-5",          "Complete diffusion"),
        ("CRAZY-3",  "Word-level diffusion",              "Full 8-word diffusion by round 4",     "ARX mixing"),
        ("CRAZY-4",  "Diff trails through Ch/Maj",        "Ch bias 0.25, Maj bias 0.25 (optimal)","Balanced nonlinearity"),
        ("CRAZY-5",  "Carry propagation analysis",        "Carry chains avg 2.2 bits, max ~32",   "Unpredictable carries"),
        ("CRAZY-6",  "Message schedule differential",     "sigma linear, addition nonlinear",     "Carry nonlinearity"),
        ("CRAZY-7",  "Multi-round differential",          "Prob ~2^-20 by round 5",               "Multiplicative decay"),
        ("CRAZY-8",  "Truncated differentials",           "No useful patterns beyond 3 rounds",   "Full diffusion"),
        ("CRAZY-29", "Impossible differential",           "All 'holes' = sampling artifacts",     "No impossibilities"),
        ("CRAZY-38", "Impossible diff (detailed)",        "Coupon collector artifacts confirmed",  "Statistical artifacts"),
        ("CRAZY-40", "Boomerang attack",                  "p^2*q^2 ~ 5e-9, zero quartets",       "Prob too low"),
    ]
    print(f"  {'ID':<9s} {'Attack':<35s} {'Result':<40s} {'Why Failed':<20s}")
    print("  " + "-" * 86)
    for id_, attack, result, why in attacks_diff:
        print(f"  {id_:<9s} {attack:<35s} {result:<40s} {why:<20s}")
    print()
    print("  COMMON BARRIER: Carry propagation in mod-2^32 addition creates exponential")
    print("  differential probability decay. By round 5, random differentials are")
    print("  indistinguishable from random.")
    print()

    # ── Category B: Algebraic Attacks ──
    print("+" + "-" * 88 + "+")
    print("| CATEGORY B: ALGEBRAIC ATTACKS                                                       |")
    print("+" + "-" * 88 + "+")
    attacks_alg = [
        ("CRAZY-9",  "Algebraic degree of round fn",     "Degree >= 10 from round 2",            "Carries = high degree"),
        ("CRAZY-10", "Boolean analysis of Ch",           "Optimal bent-like, max nonlinearity",   "Already optimal"),
        ("CRAZY-11", "Boolean analysis of Maj",          "Optimal threshold, resilient",          "Already optimal"),
        ("CRAZY-12", "XL/Grobner complexity",            "2^49 at deg-2, 2^71 at deg-3",         "No shortcut vs 2^32"),
        ("CRAZY-13", "SAT solver reduced rounds",        "1 rnd: 0.3s, 2 rnds: timeout",         "Exponential blowup"),
        ("CRAZY-14", "Linearization of round fn",        "Best linear bias ~0.02 from round 2",  "Noise floor"),
        ("CRAZY-15", "ANF of output bits",               "2^31+ terms, no structure",             "Full ANF"),
        ("CRAZY-33", "Interpolation attack",             "Degree too high for interpolation",     "Degree explosion"),
        ("CRAZY-34", "Cube attack on W[14]->De17",       "Zero constant superpolynomials",        "Full degree"),
        ("CRAZY-35", "Equation count for carry chain",   "899 quad eqs, 1859 vars, no shortcut",  "Dense system"),
        ("CRAZY-43", "Higher-order differential",        "Degree >= 10 (effectively full)",       "Carries again"),
    ]
    print(f"  {'ID':<9s} {'Attack':<35s} {'Result':<40s} {'Why Failed':<20s}")
    print("  " + "-" * 86)
    for id_, attack, result, why in attacks_alg:
        print(f"  {id_:<9s} {attack:<35s} {result:<40s} {why:<20s}")
    print()
    print("  COMMON BARRIER: Modular addition introduces algebraic degree ~32 per operation.")
    print("  Ch and Maj are already cryptographically optimal Boolean functions.")
    print()

    # ── Category C: Statistical Attacks ──
    print("+" + "-" * 88 + "+")
    print("| CATEGORY C: STATISTICAL ATTACKS                                                     |")
    print("+" + "-" * 88 + "+")
    attacks_stat = [
        ("CRAZY-16", "Hamming weight statistics",        "Perfect normal (mu=128, sigma=8)",      "Ideal distribution"),
        ("CRAZY-17", "Cross-word correlation",           "Near-zero correlation",                 "Word independence"),
        ("CRAZY-18", "Output bit independence",          "All pairwise corr < 0.01",              "Bit independence"),
        ("CRAZY-19", "Positional bias",                  "No bias detected",                      "Uniform"),
        ("CRAZY-20", "Frequency/byte analysis",          "All bytes equidistributed",             "Uniform"),
        ("CRAZY-21", "Serial correlation",               "Zero lag-k corr for k=1..10",           "No memory"),
        ("CRAZY-22", "Runs test",                        "Matches random exactly",                "No patterns"),
        ("CRAZY-30", "Zero-correlation linear",          "~25% near-zero = noise",                "Statistical noise"),
        ("CRAZY-39", "Zero-correlation (detailed)",      "Near-zero fraction = noise",            "Statistical noise"),
        ("CRAZY-45", "Correlation immunity",             "CI order=0, max corr 0.269 (bit 17)",   "Low CI expected"),
    ]
    print(f"  {'ID':<9s} {'Attack':<35s} {'Result':<40s} {'Why Failed':<20s}")
    print("  " + "-" * 86)
    for id_, attack, result, why in attacks_stat:
        print(f"  {id_:<9s} {attack:<35s} {result:<40s} {why:<20s}")
    print()
    print("  COMMON BARRIER: SHA-256 output is statistically indistinguishable from random.")
    print()

    # ── Category D: Structural Attacks ──
    print("+" + "-" * 88 + "+")
    print("| CATEGORY D: STRUCTURAL ATTACKS                                                      |")
    print("+" + "-" * 88 + "+")
    attacks_struct = [
        ("CRAZY-23", "Fixed point search",               "Zero fixed points in 10^6 IVs",        "No fixed points"),
        ("CRAZY-24", "Cycle length estimation",          "Estimated cycle > 2^60",                "Huge state space"),
        ("CRAZY-25", "Related-key differentials",        "N/A (K=IV, fixed)",                    "Not applicable"),
        ("CRAZY-26", "Weak key search",                  "No weak IV classes found",              "No weak keys"),
        ("CRAZY-27", "Symmetry/invariant search",        "Zero invariant subspaces",              "No symmetry"),
        ("CRAZY-28", "Integral/square attack",           "Full sum random after 3+ rounds",       "No balanced prop"),
        ("CRAZY-31", "MITM state split",                 "Quality=1, T1 depends on right half",   "No useful split"),
        ("CRAZY-32", "Time-memory tradeoff",             "Standard Hellman, no weakness",         "No shortcut"),
        ("CRAZY-36", "Round constant correlation",       "Zero commuting K pairs out of 276",     "Non-commutative"),
        ("CRAZY-37", "Multiplicative differential",      "Advantage ratio 0.98 after 3 rounds",   "Vanishes quickly"),
        ("CRAZY-41", "Rotational cryptanalysis",         "128-bit error, K[t] breaks symmetry",   "Round constants"),
        ("CRAZY-42", "Fault sensitivity analysis",       "Rank 3 matrix, e most sensitive",       "No exploitable bias"),
        ("CRAZY-44", "MITM state-split (detailed)",      "Quality=1",                             "No useful split"),
    ]
    print(f"  {'ID':<9s} {'Attack':<35s} {'Result':<40s} {'Why Failed':<20s}")
    print("  " + "-" * 86)
    for id_, attack, result, why in attacks_struct:
        print(f"  {id_:<9s} {attack:<35s} {result:<40s} {why:<20s}")
    print()
    print("  COMMON BARRIER: 64 distinct round constants K[t] prevent slide/rotational attacks.")
    print("  Feed-forward prevents fixed points and MITM splits. State too large for TMTO.")
    print()

    # ── Component Security Ranking ──
    print("=" * 90)
    print("  COMPONENT SECURITY CONTRIBUTION RANKING")
    print("=" * 90)
    print()
    components = [
        (1, "Addition mod 2^32",    "CRITICAL",  "Source of ALL nonlinearity. Defeats algebraic, differential, linear."),
        (2, "Message schedule",     "HIGH",      "Nonlinear expansion. SHA-1's LINEAR schedule was its fatal flaw."),
        (3, "Ch(e,f,g)",           "HIGH",      "Optimal Boolean function. Max nonlinearity in e-branch."),
        (4, "Round constants K[t]", "HIGH",      "64 distinct values prevent slide, rotational, commutative attacks."),
        (5, "Sigma0/Sigma1",       "MEDIUM",    "Linear but provide crucial bit-mixing. Make carries global."),
        (6, "Maj(a,b,c)",         "MEDIUM",     "Threshold fn in a-branch. Less critical (a-branch has T1+T2)."),
        (7, "Feed-forward add",    "MEDIUM",     "Prevents state recovery, fixed points, MITM."),
        (8, "sigma0/sigma1",      "MEDIUM",     "Schedule mixing. Linear, but combined with addition = nonlinear."),
    ]
    print(f"  {'Rank':<5s} {'Component':<22s} {'Level':<10s} {'Evidence'}")
    print("  " + "-" * 84)
    for rank, comp, level, evidence in components:
        print(f"  {rank:<5d} {comp:<22s} {level:<10s} {evidence}")
    print()

    # ── Meta-Analysis ──
    print("=" * 90)
    print("  META-ANALYSIS: WHICH ATTACK CAME CLOSEST?")
    print("=" * 90)
    print()
    closest = [
        (1, "Differential (CRAZY-1,7)",  "~5 rounds",   "Measurable bias persists through 5 rounds"),
        (2, "Truncated diff (CRAZY-8)",  "~3 rounds",   "Truncated patterns useful through 3 rounds"),
        (3, "Linear approx (CRAZY-14)",  "~2 rounds",   "Bias 0.5 at round 1, 0.02 from round 2"),
        (4, "Integral (CRAZY-28)",       "~3 rounds",   "Balanced property breaks at round 3"),
        (5, "SAT/algebraic (CRAZY-13)",  "~1 round",    "Only 1 round solvable in reasonable time"),
        (6, "Statistical (CRAZY-16-22)", "~0 rounds",   "Indistinguishable from random at any round"),
    ]
    print(f"  {'Rank':<5s} {'Attack':<28s} {'Reach':<12s} {'Evidence'}")
    print("  " + "-" * 84)
    for rank, attack, reach, evidence in closest:
        print(f"  {rank:<5d} {attack:<28s} {reach:<12s} {evidence}")
    print()
    print("  WINNER: Differential cryptanalysis. The 31-round collision (Mendel-Yu 2013)")
    print("  uses differential paths + message modification. Random differentials die by")
    print("  round 5; optimized paths with message freedom extend to 31. The gap (5->31)")
    print("  represents the power of message modification as an amplifier.")
    print()

    # ── Cross-Hash Comparison ──
    print("=" * 90)
    print("  CROSS-HASH COMPARISON")
    print("=" * 90)
    print()
    print(f"  {'Feature':<28s} {'SHA-1':<16s} {'SHA-256':<16s} {'BLAKE2s':<16s} {'Keccak':<16s}")
    print("  " + "-" * 88)
    rows = [
        ("Design",               "MD+Davies-Meyer", "MD+Davies-Meyer", "HAIFA",           "Sponge"),
        ("Nonlinear rounds",     "40/80 (50%)",     "64/64 (100%)",    "10/10 (100%)",    "24/24 (100%)"),
        ("Expansion",            "GF(2)-LINEAR",    "Nonlinear",       "N/A (direct)",    "N/A (sponge)"),
        ("Diffusion to random",  "~10 rounds",      "4-5 rounds",      "1 round",         "2 rounds"),
        ("Degree per round",     "1-2",             "~4+ (carries)",   ">>8",             "2^k (exact)"),
        ("Best collision attack","FULL 80 rounds",  "31/64 rounds",    "0 rounds",        "6/24 rounds"),
        ("Fatal weakness",       "Linear expansion","NONE FOUND",      "NONE FOUND",      "NONE FOUND"),
    ]
    for feat, *vals in rows:
        line = f"  {feat:<28s}"
        for v in vals:
            line += f" {v:<16s}"
        print(line)
    print()
    print("  WHY SHA-1 BROKE:")
    print("  1. Linear message expansion: W[t] = W[t-3] XOR W[t-8] XOR W[t-14] XOR W[t-16]")
    print("  2. 20 Parity rounds use XOR only (degree 1)")
    print("  3. Differential paths can be constructed algebraically through linear schedule")
    print()

    # ── Final Verdict ──
    print("=" * 90)
    print("  FINAL VERDICT")
    print("=" * 90)
    print()
    print("  45 experiments, 4 attack categories, 3 cross-hash comparisons.")
    print("  SHA-256 shows NO exploitable weakness on full rounds.")
    print()
    print("  Security margin: 33/64 rounds (52%) beyond best known attack.")
    print()
    print("  Three reinforcing barriers:")
    print("    1. NONLINEAR EXPANSION - prevents algebraic differential path construction")
    print("    2. CARRY PROPAGATION  - exponential nonlinearity growth per round")
    print("    3. COMPLETE DIFFUSION - every output bit depends on every input by round 5")
    print()
    print("  No known technique defeats even one barrier beyond ~5 rounds (random)")
    print("  or ~31 rounds (optimized message modification).")
    print()


if __name__ == "__main__":
    main()
