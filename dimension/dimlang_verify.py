"""
DimLang VERIFICATION: отделяем реальные открытия от статистических артефактов.

Подозрительные результаты:
  - "MITM advantage 28 bits" → может быть неправильная формула expected
  - "Gradient 33 bits за 1 шаг" → реально?
  - "Schedule autocorrelation 0.228" → реально?
  - "Min dist 91" → ожидаемо для birthday?

Тестируем на RANDOM FUNCTION (не SHA-256) для сравнения.
"""

import numpy as np
from scipy import stats
from dimlang import *
from dimlang import hw

np.random.seed(42)


def random_hash(msg_raw):
    """'Random hash': SHA-256-like but replace with actual random mapping."""
    # Use numpy hash as proxy for random function
    h = np.random.RandomState(hash(tuple(msg_raw)) % (2**32)).randint(0, 2**32, 8)
    return State(h.tolist())


def verify_min_dist_expectation():
    """Правильная формула для expected min distance."""
    print(f"{'=' * 70}")
    print("1. CORRECT EXPECTED MIN DISTANCE")
    print("=" * 70)

    # For N pairs, each pair has HW(XOR) ~ Binomial(256, 0.5)
    # mean = 128, std = 8
    # Min of N samples from this distribution:
    # Expected min ≈ 128 - 8 * Φ^(-1)(1 - 1/N)

    from scipy.stats import norm

    print(f"\n  Theoretical expected min distance for N random 256-bit pairs:")
    print(f"  {'N':>10} {'Expected min':>13} {'Our formula':>12}")
    for N in [100, 1000, 10000, 40000, 100000, 1000000]:
        # Correct: using order statistics
        p = 1.0 / N
        expected_min = 128 - 8 * norm.ppf(1 - p)
        our_formula = 128 - np.log2(N) / 2
        print(f"  {N:>10} {expected_min:>12.1f} {our_formula:>11.1f}")

    print(f"\n  Our formula was WRONG by ~20-30 bits!")
    print(f"  Correct expected for 40K pairs: {128 - 8 * norm.ppf(1 - 1/40000):.1f}")
    print(f"  We got: 91-97 → WITHIN EXPECTED RANGE!")


def verify_gradient():
    """Is the 33-bit gradient improvement REAL or statistical artifact?"""
    print(f"\n{'=' * 70}")
    print("2. VERIFY GRADIENT")
    print("=" * 70)

    # Test: for RANDOM function, does single-bit flip also give 33-bit improvement?
    # Compare SHA-256 vs random
    print(f"\n  SHA-256 single-bit gradient (10 trials):")
    sha_improvements = []
    for _ in range(10):
        m1 = Message()
        m2 = Message()
        p1 = Path(m1)
        p2 = Path(m2)
        d_base = overlay(p1, p2).weight

        best_imp = 0
        for word in range(16):
            for bit in range(32):
                m_test = m1.copy()
                m_test[word] = m_test[word].val ^ (1 << bit)
                d_test = overlay(Path(m_test), p2).weight
                imp = d_base - d_test
                if imp > best_imp: best_imp = imp
        sha_improvements.append(best_imp)

    print(f"    Mean best improvement: {np.mean(sha_improvements):.1f}")
    print(f"    Std: {np.std(sha_improvements):.1f}")
    print(f"    Min: {min(sha_improvements)}, Max: {max(sha_improvements)}")

    # Expected: for 512 trials of Binomial(256,0.5) compared to baseline,
    # what's the expected max improvement?
    # Each trial: new HW ~ Binomial(256, 0.5), improvement = d_base - new_HW
    # Max of 512 such: E[max] ≈ d_base - (128 - 8*Φ^(-1)(1-1/512))

    from scipy.stats import norm
    d_base_typical = 128
    expected_max_imp = 8 * norm.ppf(1 - 1.0/512)
    print(f"\n    Expected max improvement (512 random trials): {expected_max_imp:.1f}")
    print(f"    → SHA-256 gradient = {'ANOMALOUS!' if np.mean(sha_improvements) > expected_max_imp + 5 else 'EXPECTED (not special).'}")


def verify_mitm():
    """Is MITM advantage real?"""
    print(f"\n{'=' * 70}")
    print("3. VERIFY MITM ADVANTAGE")
    print("=" * 70)

    from scipy.stats import norm

    # SHA-256 MITM
    print(f"\n  SHA-256 MITM (200×200 at r=32):")
    sha_min_dists = []
    for trial in range(5):
        forward = []
        for _ in range(200):
            p = Path(Message())
            forward.append(p.at(32))

        target_path = Path(Message())
        backward = []
        for _ in range(200):
            s = State(target_path.states[64].regs)
            for r in range(63, 31, -1):
                s = invert_round(s, Word(np.random.randint(0, 2**32)), r)
            backward.append(s)

        min_d = 256
        for i in range(200):
            for j in range(200):
                d = forward[i].delta(backward[j]).weight
                if d < min_d: min_d = d
        sha_min_dists.append(min_d)

    print(f"    5 trials: {sha_min_dists}")
    print(f"    Mean: {np.mean(sha_min_dists):.1f}")

    # PURE RANDOM: same exercise with random 256-bit vectors
    print(f"\n  RANDOM (same N=200×200 = 40K pairs):")
    rand_min_dists = []
    for trial in range(5):
        # 200 random 256-bit vectors (forward)
        fw = [np.random.randint(0, 2, 256) for _ in range(200)]
        # 200 random 256-bit vectors (backward)
        bw = [np.random.randint(0, 2, 256) for _ in range(200)]

        min_d = 256
        for i in range(200):
            for j in range(200):
                d = np.sum(fw[i] != bw[j])
                if d < min_d: min_d = d
        rand_min_dists.append(min_d)

    print(f"    5 trials: {rand_min_dists}")
    print(f"    Mean: {np.mean(rand_min_dists):.1f}")

    expected = 128 - 8 * norm.ppf(1 - 1.0/40000)
    print(f"\n    Expected min for 40K pairs: {expected:.1f}")
    print(f"    SHA-256 mean: {np.mean(sha_min_dists):.1f}")
    print(f"    Random mean:  {np.mean(rand_min_dists):.1f}")
    advantage = np.mean(rand_min_dists) - np.mean(sha_min_dists)
    print(f"    SHA-256 advantage over random: {advantage:.1f} bits")
    print(f"    → {'REAL ADVANTAGE!' if advantage > 5 else 'NO ADVANTAGE (statistical noise).'}")


def verify_schedule_autocorrelation():
    """Is schedule autocorrelation real and exploitable?"""
    print(f"\n{'=' * 70}")
    print("4. VERIFY SCHEDULE AUTOCORRELATION")
    print("=" * 70)

    # Many messages, measure autocorrelation of expanded schedule
    autocorrs_by_lag = {1: [], 2: [], 7: [], 15: [], 16: []}

    for _ in range(100):
        m = Message()
        W = m.expand()
        vals = np.array([w.val / 2**32 for w in W])  # normalize to [0,1]
        vals -= np.mean(vals)
        norm_sq = np.sum(vals**2)
        if norm_sq == 0: continue

        for lag in autocorrs_by_lag:
            if lag < 64:
                c = np.sum(vals[:64-lag] * vals[lag:64]) / norm_sq
                autocorrs_by_lag[lag].append(c)

    print(f"\n  Schedule autocorrelation (100 random messages):")
    for lag in sorted(autocorrs_by_lag):
        vals = autocorrs_by_lag[lag]
        if vals:
            mean_c = np.mean(vals)
            std_c = np.std(vals)
            significant = abs(mean_c) > 2 * std_c / np.sqrt(len(vals))
            print(f"    lag={lag:2d}: mean={mean_c:+.4f} ± {std_c:.4f}  {'★ SIGNIFICANT' if significant else ''}")

    # Compare with random (no schedule structure)
    print(f"\n  RANDOM comparison (no schedule structure):")
    rand_autocorrs = {1: [], 16: []}
    for _ in range(100):
        vals = np.random.randn(64)
        norm_sq = np.sum(vals**2)
        for lag in rand_autocorrs:
            c = np.sum(vals[:64-lag] * vals[lag:64]) / norm_sq
            rand_autocorrs[lag].append(c)

    for lag in rand_autocorrs:
        vals = rand_autocorrs[lag]
        print(f"    lag={lag:2d}: mean={np.mean(vals):+.4f} ± {np.std(vals):.4f}")


def verify_repair_advantage():
    """Is repair advantage real?"""
    print(f"\n{'=' * 70}")
    print("5. VERIFY REPAIR ADVANTAGE")
    print("=" * 70)

    # Repair: random δW[0], measure HW(δH)
    # Compare: random δW (all 16 words), measure HW(δH)
    # Both use 1 call to SHA-256

    print(f"\n  Repair (strategy 'a', 10K trials):")
    repair_hws = []
    for _ in range(10000):
        m = Message()
        dw0 = np.random.randint(1, 2**32)
        _, dH = dual_repair(m, dw0, 'a')
        repair_hws.append(dH)

    print(f"    Mean: {np.mean(repair_hws):.1f}")
    print(f"    Std:  {np.std(repair_hws):.1f}")
    print(f"    Min:  {min(repair_hws)}")

    print(f"\n  Random δW (single bit, 10K trials):")
    random_hws = []
    for _ in range(10000):
        m = Message()
        p1 = Path(m)
        m2 = m.copy()
        w = np.random.randint(0, 16)
        b = np.random.randint(0, 32)
        m2[w] = m2[w].val ^ (1 << b)
        p2 = Path(m2)
        random_hws.append(overlay(p1, p2).weight)

    print(f"    Mean: {np.mean(random_hws):.1f}")
    print(f"    Std:  {np.std(random_hws):.1f}")
    print(f"    Min:  {min(random_hws)}")

    print(f"\n  Random δW (random 16 words, 10K trials):")
    random_full_hws = []
    for _ in range(10000):
        m1 = Message()
        m2 = Message()
        p1 = Path(m1)
        p2 = Path(m2)
        random_full_hws.append(overlay(p1, p2).weight)

    print(f"    Mean: {np.mean(random_full_hws):.1f}")
    print(f"    Std:  {np.std(random_full_hws):.1f}")
    print(f"    Min:  {min(random_full_hws)}")

    # Statistical test
    t_stat, p_val = stats.ttest_ind(repair_hws, random_hws)
    print(f"\n  T-test repair vs random single-bit: t={t_stat:.2f}, p={p_val:.6f}")
    print(f"  → {'SIGNIFICANT' if p_val < 0.01 else 'NOT SIGNIFICANT'}")

    # Key: repair uses structured δW (only W[0] differs + compensation in W[1..15])
    # vs random uses single-bit δW (all compensation through natural hash)
    # If repair ≈ random → repair adds no value over just random 1-bit flip


def main():
    verify_min_dist_expectation()
    verify_gradient()
    verify_mitm()
    verify_schedule_autocorrelation()
    verify_repair_advantage()

    print(f"\n{'=' * 70}")
    print("TRUTH TABLE: what's real, what's artifact")
    print("=" * 70)
    print(f"""
  Claimed finding               Real?    Evidence
  ─────────────────────────────  ──────   ────────
  MITM 28-bit advantage         ???      Need to compare with random
  Gradient 33 bits / step       ???      Need to compare with random
  Schedule autocorrelation      ???      Need statistical significance
  Repair advantage              ???      Need comparison with random
  Min dist = 91                 ???      Expected from birthday stats

  The HONEST answer requires COMPARING with random functions.
  If SHA-256 gives the SAME results as random → no advantage.
  If SHA-256 gives DIFFERENT results → real structure to exploit.
""")


if __name__ == "__main__":
    main()
