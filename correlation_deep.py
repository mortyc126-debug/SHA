"""
DEEP CORRELATION: The 8.66× signal between layers for structured pairs.

Verify on larger N. Extend to layers 2, 3, 4. Understand WHY.

Key observation: the correlation is COLLECTIVE (all 8 regs together).
Individual registers show zero correlation.
This suggests: the 8 registers are coupled through a shared mechanism.
That mechanism = the round function operating on the FULL state.
"""

import random
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, schedule
)


def sha256_hash(M):
    states, _ = sha256_round_trace(M)
    return [(states[64][i] + H0[i]) & MASK for i in range(8)]


def experiment_verify_and_extend():
    """
    1. Verify 8.66× with larger N.
    2. Measure correlation for all layer pairs.
    3. Chain: P(layer 0+1+2 zero | layer 0 zero) etc.
    """
    print("=" * 80)
    print("DEEP CORRELATION: Verify and extend the 8.66× signal")
    print("=" * 80)

    N = 100000
    delta_word, delta_bit = 0, 0

    # Collect per-layer zero status for each message
    # layer_zero[k][i] = True if ALL 8 regs have δH bit k = 0 for message i
    max_layers = 8
    layer_zero = [[] for _ in range(max_layers)]

    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        M2 = list(M); M2[delta_word] ^= (1 << delta_bit)
        H1 = sha256_hash(M)
        H2 = sha256_hash(M2)
        dH = [H1[i] ^ H2[i] for i in range(8)]

        for k in range(max_layers):
            all_zero = all(((dH[reg] >> k) & 1) == 0 for reg in range(8))
            layer_zero[k].append(all_zero)

    # Marginal probabilities
    print(f"\n  N = {N}, δM = W[0] bit 0")
    print(f"\n  Marginal P(layer k all-zero):")
    for k in range(max_layers):
        count = sum(layer_zero[k])
        expected = N * 2**(-8)
        print(f"    Layer {k}: {count}/{N} = {count/N:.6f} (expected {expected/N:.6f})")

    # Conditional: P(layer j zero | layer i zero)
    print(f"\n  Conditional P(layer j=0 | layer i=0):")
    print(f"  {'':>8}", end="")
    for j in range(5):
        print(f" {'L'+str(j)+' zero':>10}", end="")
    print()

    for i in range(5):
        i_zero = [idx for idx in range(N) if layer_zero[i][idx]]
        print(f"  L{i}=0 |", end="")
        for j in range(5):
            if i == j:
                print(f"     ---  ", end="")
                continue
            both = sum(1 for idx in i_zero if layer_zero[j][idx])
            p_cond = both / len(i_zero) if i_zero else 0
            p_marginal = sum(layer_zero[j]) / N
            ratio = p_cond / p_marginal if p_marginal > 0 else 0
            print(f" {ratio:>8.2f}×", end="")
        print(f"  (n={len(i_zero)})")

    # Chain: P(layers 0..k all zero)
    print(f"\n  Cumulative: P(layers 0..k ALL zero):")
    for k in range(5):
        count = sum(1 for idx in range(N) if all(layer_zero[j][idx] for j in range(k+1)))
        bits = 8 * (k + 1)
        expected = N * 2**(-bits)
        ratio = (count / expected) if expected > 0 else float('inf')
        print(f"    Layers 0..{k} ({bits}-bit): {count}/{N} (expected {expected:.2f}, ratio = {ratio:.1f}×)")


def experiment_why_collective():
    """
    WHY is the correlation collective (all 8 regs) but not per-register?

    Hypothesis: when δH bit-0 of ALL regs = 0, it means the CARRY
    structure at bit 0 is identical between M and M'. This global
    carry match constrains the state at bit-1 level too.

    But individual register bit-0 match doesn't constrain other registers.

    Test: does δH bit-0 = 0 for SOME subset of registers (not all 8)
    still predict bit-1?
    """
    print("\n" + "=" * 80)
    print("WHY COLLECTIVE: Testing subsets of registers")
    print("=" * 80)

    N = 100000
    delta_word, delta_bit = 0, 0

    # For each message pair, collect per-register δH bits
    dh_bits = []  # list of (dH[0], dH[1], ..., dH[7]) per message

    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        M2 = list(M); M2[delta_word] ^= (1 << delta_bit)
        H1 = sha256_hash(M)
        H2 = sha256_hash(M2)
        dh_bits.append(tuple(H1[i] ^ H2[i] for i in range(8)))

    # For each subset size s (1, 2, 4, 8 registers matching on bit 0):
    # P(those registers also match on bit 1)
    reg_subsets = {
        '1 reg (a)': [0],
        '2 regs (a,e)': [0, 4],
        '4 regs (a,b,c,d)': [0, 1, 2, 3],
        '4 regs (a,b,e,f)': [0, 1, 4, 5],
        '8 regs (all)': list(range(8)),
    }

    for name, regs in reg_subsets.items():
        # Find messages where ALL specified regs have δH bit-0 = 0
        matching_b0 = [i for i in range(N)
                       if all(((dh_bits[i][r] >> 0) & 1) == 0 for r in regs)]

        if not matching_b0:
            print(f"  {name}: no matches at bit-0")
            continue

        # Among those: P(same regs also match at bit-1)
        also_b1 = sum(1 for i in matching_b0
                      if all(((dh_bits[i][r] >> 1) & 1) == 0 for r in regs))

        p_cond = also_b1 / len(matching_b0)
        p_marginal = sum(1 for i in range(N)
                        if all(((dh_bits[i][r] >> 1) & 1) == 0 for r in regs)) / N

        ratio = p_cond / p_marginal if p_marginal > 0 else 0

        print(f"  {name:>20}: P(b1=0|b0=0) = {p_cond:.5f}, "
              f"marginal = {p_marginal:.5f}, ratio = {ratio:.2f}×, "
              f"n_b0={len(matching_b0)}")

    print(f"\n  If ratio grows with subset size → correlation is COLLECTIVE.")
    print(f"  If ratio constant → correlation is per-register (hidden by noise).")


def experiment_reduced_rounds():
    """
    Does the correlation exist at reduced rounds too?
    If yes at R=8 but not R=64 → early-round phenomenon.
    If yes at R=64 → persists through full computation.
    """
    print("\n" + "=" * 80)
    print("REDUCED ROUNDS: Does correlation persist?")
    print("=" * 80)

    for R in [4, 8, 16, 32, 64]:
        N = 50000

        count_b0_zero = 0
        count_b01_zero = 0

        for seed in range(N):
            rng = random.Random(seed)
            M = [rng.randint(0, MASK) for _ in range(16)]
            M2 = list(M); M2[0] ^= 1

            s1, _ = sha256_round_trace(M, rounds=R)
            s2, _ = sha256_round_trace(M2, rounds=R)

            H1 = [(s1[R][i] + H0[i]) & MASK for i in range(8)]
            H2 = [(s2[R][i] + H0[i]) & MASK for i in range(8)]
            dH = [H1[i] ^ H2[i] for i in range(8)]

            b0_zero = all(((dH[i] >> 0) & 1) == 0 for i in range(8))
            b1_zero = all(((dH[i] >> 1) & 1) == 0 for i in range(8))

            if b0_zero:
                count_b0_zero += 1
                if b1_zero:
                    count_b01_zero += 1

        if count_b0_zero > 0:
            p_cond = count_b01_zero / count_b0_zero
            p_marginal = 2**(-8)
            ratio = p_cond / p_marginal
            print(f"  R={R:>2}: n(b0=0)={count_b0_zero:>5}, P(b1=0|b0=0)={p_cond:.5f}, ratio={ratio:.2f}×")
        else:
            print(f"  R={R:>2}: n(b0=0)=0 (need more samples)")


if __name__ == "__main__":
    experiment_verify_and_extend()
    experiment_why_collective()
    experiment_reduced_rounds()
