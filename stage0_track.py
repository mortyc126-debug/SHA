"""
Stage 0, Part 2: Track WHERE specific bit differences go.

Not "how many" but "which ones and where".

Key question: Is there a PATTERN in which positions survive vs die
across rounds? Or is it truly random?
"""

import random

# Import from stage0_observe
from stage0_observe import (
    sha256_round_trace, state_bits, xor_bits, MASK, H0
)


def experiment_5_position_tracking():
    """
    Track each differing bit position across rounds.

    For each round r, record the SET of differing positions.
    Then: for position p that differs at round r, what is the probability
    it still differs at round r+1? r+2? r+k?

    If SHA-256 is truly "random-like", survival probability = 0.5 per round.
    If there's structure, some positions will survive longer.
    """
    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]
    M_prime = list(M)
    M_prime[0] ^= 1

    states_M, _ = sha256_round_trace(M)
    states_Mp, _ = sha256_round_trace(M_prime)

    # Get diff sets for each round
    diff_sets = []
    for r in range(65):
        s1 = state_bits(states_M[r])
        s2 = state_bits(states_Mp[r])
        diff = set()
        for i in range(256):
            if s1[i] != s2[i]:
                diff.add(i)
        diff_sets.append(diff)

    print("=" * 80)
    print("EXPERIMENT 5: Position survival analysis")
    print("=" * 80)

    # For each position that's different at round 10 (after stabilization),
    # how long does it stay different?
    print("\n--- Survival of positions differing at round 10 ---")
    start_round = 10
    start_set = diff_sets[start_round]
    print(f"Positions differing at round {start_round}: {len(start_set)}")

    for delta_r in [1, 2, 3, 4, 5, 10, 20, 30, 50]:
        target_r = start_round + delta_r
        if target_r > 64:
            break
        survivors = start_set & diff_sets[target_r]
        p_survive = len(survivors) / len(start_set)
        expected_random = 0.5  # If independent
        print(f"  Δr={delta_r:>2}: {len(survivors):>3}/{len(start_set)} survived = {p_survive:.3f} "
              f"(random expectation: {expected_random:.3f})")

    # Key insight: track the IDENTITY of surviving positions
    print(f"\n--- Which register do survivors come from? ---")
    reg_names = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
    for delta_r in [1, 2, 3, 4]:
        target_r = start_round + delta_r
        survivors = start_set & diff_sets[target_r]
        # Count by register
        reg_counts = [0] * 8
        for pos in survivors:
            reg_counts[pos // 32] += 1
        reg_str = " ".join(f"{reg_names[i]}:{reg_counts[i]:>2}" for i in range(8))
        print(f"  Δr={delta_r}: {reg_str}")


def experiment_6_register_shift_tracking():
    """
    SHA-256 shifts registers: b←a, c←b, d←c, f←e, g←f, h←g.
    Only a and e get NEW values each round.

    Question: Does the diff in register b at round r+1 = diff in register a at round r?
    This is EXACT for the values (b_{r+1} = a_r), but is it exact for the DIFFS?

    If delta_a_r = a_r(M) XOR a_r(M'), then delta_b_{r+1} should equal delta_a_r.
    Let's verify this and see what it means.
    """
    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]
    M_prime = list(M)
    M_prime[0] ^= 1

    states_M, _ = sha256_round_trace(M)
    states_Mp, _ = sha256_round_trace(M_prime)

    print("\n" + "=" * 80)
    print("EXPERIMENT 6: Register shift tracking — is δb_{r+1} = δa_r exactly?")
    print("=" * 80)

    reg_names = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']

    # Check: δb_{r+1} == δa_r, δc_{r+1} == δb_r, etc.
    shifts = [(1, 0), (2, 1), (3, 2), (5, 4), (6, 5), (7, 6)]
    shift_names = ['b←a', 'c←b', 'd←c', 'f←e', 'g←f', 'h←g']

    print(f"\nVerification that register shifts preserve diffs exactly:")
    for r in range(64):
        for (dst, src), name in zip(shifts, shift_names):
            delta_src = states_M[r][src] ^ states_Mp[r][src]
            delta_dst = states_M[r+1][dst] ^ states_Mp[r+1][dst]
            if delta_src != delta_dst:
                print(f"  VIOLATION at round {r}: {name}: δ{reg_names[src]}_r=0x{delta_src:08x} != δ{reg_names[dst]}_{{r+1}}=0x{delta_dst:08x}")
                break
        else:
            continue
        break
    else:
        print(f"  All 64 rounds × 6 shifts: EXACT match (384/384)")

    # So: the ONLY new information each round comes from δa and δe.
    # Everything else is just shifted copies.
    # This means the 256-bit state diff at round r is:
    #   (δa_r, δa_{r-1}, δa_{r-2}, δa_{r-3}, δe_r, δe_{r-1}, δe_{r-2}, δe_{r-3})
    #
    # The ENTIRE evolution is captured by just TWO 32-bit sequences: δa_r and δe_r.

    print(f"\n--- Consequence: state diff = (δa_r, δa_{{r-1}}, δa_{{r-2}}, δa_{{r-3}}, δe_r, δe_{{r-1}}, δe_{{r-2}}, δe_{{r-3}})")
    print(f"--- All 256 bits of diff are determined by just δa and δe sequences ---")
    print(f"\nδa and δe sequences (HW):")
    print(f"{'r':>3} | {'δa':>10} {'HW(δa)':>6} | {'δe':>10} {'HW(δe)':>6}")
    print("-" * 45)
    for r in range(min(20, 65)):
        da = states_M[r][0] ^ states_Mp[r][0]
        de = states_M[r][4] ^ states_Mp[r][4]
        print(f"{r:>3} | 0x{da:08x} {bin(da).count('1'):>6} | 0x{de:08x} {bin(de).count('1'):>6}")


def experiment_7_ae_recurrence():
    """
    Since the full state diff is determined by (δa_r, δe_r),
    let's find the EXACT recurrence relation.

    From SHA-256:
      a_new = T1 + T2
      e_new = d + T1

    So:
      δa_r = (T1(M) + T2(M)) XOR (T1(M') + T2(M'))
      δe_r = (d(M) + T1(M)) XOR (d(M') + T1(M'))

    where:
      T1 = h + Σ₁(e) + Ch(e,f,g) + K + W
      T2 = Σ₀(a) + Maj(a,b,c)
      d = c_{r-1} = b_{r-2} = a_{r-3}  →  δd_r = δa_{r-3}
      h = g_{r-1} = f_{r-2} = e_{r-3}  →  δh_r = δe_{r-3}

    This means δa and δe at round r depend on:
      δa_{r-1}, δa_{r-2}, δa_{r-3}  (through b,c,d = shifted a's)
      δe_{r-1}, δe_{r-2}, δe_{r-3}  (through f,g,h = shifted e's)
      δW[r]                           (through schedule)

    Plus the ACTUAL VALUES (not just diffs) through Ch, Maj, Σ₀, Σ₁, carries.

    Let's measure: how much of δa_r is "predictable" from δa_{r-1..r-3}, δe_{r-1..r-3}
    vs how much is "new" from the nonlinear operations?
    """
    rng = random.Random(42)

    print("\n" + "=" * 80)
    print("EXPERIMENT 7: How much of δa_r, δe_r is determined by recent history?")
    print("=" * 80)

    # Run multiple seeds to get statistics
    N_seeds = 200

    # For each round r >= 8 (after stabilization):
    # Compute XOR of δa_r with various "predictions" based on history
    # The "best linear prediction" would be some XOR combination of past δa, δe

    # First: measure correlation between δe_r and δa_{r-3}
    # Because: e_new = d + T1, and d = a_{r-3}
    # If T1 were constant across M, M': δe_r = δd_r = δa_{r-3}
    # The "correction" δT1 comes from nonlinear parts

    corr_e_a3 = []  # HW(δe_r XOR δa_{r-3}) — how different is δe_r from "shifted δa"
    corr_a_sum = []  # HW of the "T2 correction"

    for seed in range(N_seeds):
        rng_s = random.Random(seed)
        M = [rng_s.randint(0, MASK) for _ in range(16)]
        M_prime = list(M)
        M_prime[0] ^= 1

        states_M, _ = sha256_round_trace(M)
        states_Mp, _ = sha256_round_trace(M_prime)

        for r in range(8, 64):
            da_r = states_M[r][0] ^ states_Mp[r][0]
            de_r = states_M[r][4] ^ states_Mp[r][4]
            da_r3 = states_M[r-3][0] ^ states_Mp[r-3][0]  # δa_{r-3} = δd_r
            de_r3 = states_M[r-3][4] ^ states_Mp[r-3][4]  # δe_{r-3} = δh_r

            # δe_r vs δa_{r-3}: how much does T1-correction add?
            corr_e_a3.append(bin(de_r ^ da_r3).count('1'))

    avg_corr = sum(corr_e_a3) / len(corr_e_a3)
    print(f"\nE[HW(δe_r XOR δa_{{r-3}})] = {avg_corr:.2f}")
    print(f"  If δe_r = δa_{{r-3}} (no T1 correction): would be 0")
    print(f"  If independent random: would be 16")
    print(f"  Actual: {avg_corr:.2f} → T1 correction adds ~{avg_corr:.0f} bits of 'new' info per round to e")

    # Same for δa_r: a_new = T1 + T2
    # δa_r depends on BOTH T1 and T2 corrections
    # T2 depends on δa_{r-1}, δa_{r-2}, δa_{r-3} through Σ₀(a) and Maj(a,b,c)

    # Key measurement: δa_r XOR (δe_r) — since both use T1
    corr_a_e = []
    for seed in range(N_seeds):
        rng_s = random.Random(seed)
        M = [rng_s.randint(0, MASK) for _ in range(16)]
        M_prime = list(M)
        M_prime[0] ^= 1

        states_M, _ = sha256_round_trace(M)
        states_Mp, _ = sha256_round_trace(M_prime)

        for r in range(4, 64):
            da_r = states_M[r][0] ^ states_Mp[r][0]
            de_r = states_M[r][4] ^ states_Mp[r][4]
            da_r3 = states_M[r-3][0] ^ states_Mp[r-3][0]

            # δa_r = δ(T1+T2) and δe_r = δ(d+T1) = δ(a_{r-3}+T1)
            # So δa_r XOR δe_r "removes" T1 component, leaving T2 - d contribution
            corr_a_e.append(bin(da_r ^ de_r).count('1'))

    avg_ae = sum(corr_a_e) / len(corr_a_e)
    print(f"\nE[HW(δa_r XOR δe_r)] = {avg_ae:.2f}")
    print(f"  This is the 'T2 correction minus d correction'")
    print(f"  If T2=d (no new info from Maj/Σ₀): would be 0")
    print(f"  If independent: would be 16")
    print(f"  Actual: {avg_ae:.2f}")

    # The FUNDAMENTAL decomposition:
    # δe_r = δa_{r-3} + δT1_r  (mod addition, not XOR)
    # δa_r = δT1_r + δT2_r     (mod addition)
    # So: δa_r - δe_r = δT2_r - δa_{r-3}  (mod 2^32)
    #
    # δT1_r = δe_{r-3} + δΣ₁(e_{r-1}) + δCh(e_{r-1},e_{r-2},e_{r-3}) + δW[r]
    # δT2_r = δΣ₀(a_{r-1}) + δMaj(a_{r-1},a_{r-2},a_{r-3})
    #
    # ALL of this is determined by (a_r-1..r-3, e_r-1..r-3) VALUES and δ's.
    # Plus δW[r] from schedule.
    #
    # The DEPTH of dependence is only 3 rounds back.
    # But the VALUES carry full history.

    print(f"\n--- Key insight ---")
    print(f"δa_r and δe_r depend on values & diffs from rounds r-1, r-2, r-3 only.")
    print(f"But VALUES at r-1 carry all information from round 0.")
    print(f"The 'memory' is in the VALUES, not in the DIFFS.")
    print(f"DIFFS look random because they are filtered through value-dependent nonlinearity.")


def experiment_8_value_dependence():
    """
    The critical observation: δCh(e,f,g) depends on the VALUES of e,f,g
    (not just their diffs). Same for carries in addition.

    This means: for the SAME δ input, DIFFERENT messages produce
    DIFFERENT δ outputs. The "transfer function" changes with every message.

    Question: How much does the transfer function vary across messages?
    If it varies maximally → random-like behavior.
    If it has invariant structure → exploitable.
    """
    print("\n" + "=" * 80)
    print("EXPERIMENT 8: Transfer function variation across messages")
    print("=" * 80)

    # Fix δ = flip bit 0 of W[0].
    # Run 1000 different base messages.
    # At round 10, record δa_10.
    # Question: how are the δa_10 values distributed?

    N = 1000
    da_values = []
    de_values = []

    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        M_prime = list(M)
        M_prime[0] ^= 1

        states_M, _ = sha256_round_trace(M, rounds=10)
        states_Mp, _ = sha256_round_trace(M_prime, rounds=10)

        da_values.append(states_M[10][0] ^ states_Mp[10][0])
        de_values.append(states_M[10][4] ^ states_Mp[10][4])

    # How many unique values?
    unique_da = len(set(da_values))
    unique_de = len(set(de_values))

    print(f"\nSame δ (bit 0 of W[0]), different M, at round 10:")
    print(f"  Unique δa_10 values: {unique_da}/{N}")
    print(f"  Unique δe_10 values: {unique_de}/{N}")

    # HW distribution
    hw_da = [bin(v).count('1') for v in da_values]
    hw_de = [bin(v).count('1') for v in de_values]
    print(f"  E[HW(δa_10)] = {sum(hw_da)/N:.2f} (random = 16)")
    print(f"  E[HW(δe_10)] = {sum(hw_de)/N:.2f} (random = 16)")

    # Bit-level bias: for each bit position, P(bit=1)
    print(f"\n  Bit-level bias δa_10 (deviation from 0.5):")
    max_dev = 0
    for bit in range(32):
        count_1 = sum(1 for v in da_values if (v >> bit) & 1)
        p = count_1 / N
        dev = abs(p - 0.5)
        if dev > max_dev:
            max_dev = dev
    print(f"    Max deviation from 0.5: {max_dev:.4f}")

    # NOW: same experiment but at round 3 (before full diffusion)
    da_r3 = []
    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        M_prime = list(M)
        M_prime[0] ^= 1
        states_M, _ = sha256_round_trace(M, rounds=3)
        states_Mp, _ = sha256_round_trace(M_prime, rounds=3)
        da_r3.append(states_M[3][0] ^ states_Mp[3][0])

    unique_r3 = len(set(da_r3))
    hw_r3 = [bin(v).count('1') for v in da_r3]
    print(f"\n  At round 3 (before full diffusion):")
    print(f"  Unique δa_3 values: {unique_r3}/{N}")
    print(f"  E[HW(δa_3)] = {sum(hw_r3)/N:.2f}")

    max_dev_r3 = 0
    for bit in range(32):
        count_1 = sum(1 for v in da_r3 if (v >> bit) & 1)
        p = count_1 / N
        dev = abs(p - 0.5)
        if dev > max_dev_r3:
            max_dev_r3 = dev
    print(f"    Max bit deviation: {max_dev_r3:.4f}")

    print(f"\n--- Comparison ---")
    print(f"  Round 3:  {unique_r3} unique, max_dev={max_dev_r3:.4f}")
    print(f"  Round 10: {unique_da} unique, max_dev={max_dev:.4f}")
    print(f"  → Transfer function becomes maximally variable by round 10")


if __name__ == "__main__":
    experiment_5_position_tracking()
    experiment_6_register_shift_tracking()
    experiment_7_ae_recurrence()
    experiment_8_value_dependence()
