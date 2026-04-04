"""
STRUCTURED COLLISION: layers may be correlated for FIXED δM.

Previous test: random pairs → layers independent. This is BIRTHDAY = 2^128.

New test: FIXED δM (e.g., flip bit 0 of W[0]).
For M and M' = M ⊕ δM: are layers correlated?

If yes → structured pairs can exploit layer correlation → better than birthday.
"""

import random
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, schedule
)


def sha256_hash(M):
    states, _ = sha256_round_trace(M)
    return [(states[64][i] + H0[i]) & MASK for i in range(8)]


def experiment_structured_layers():
    """
    Fix δM = flip bit 0 of W[0] (single-bit XOR difference).
    For many M, compute H(M) and H(M ⊕ δM).
    δH = H(M) ⊕ H(M ⊕ δM).

    Question: are bit-0 of δH and bit-1 of δH CORRELATED?
    If yes → layers are NOT independent under structured differences.
    """
    print("=" * 80)
    print("STRUCTURED PAIRS: Are layers correlated for fixed δM?")
    print("=" * 80)

    N = 10000

    # δM = flip bit 0 of W[0]
    delta_word, delta_bit = 0, 0

    dh_bit0 = []  # bit 0 of δH[0..7] for each M
    dh_bit1 = []  # bit 1 of δH[0..7] for each M
    dh_bit2 = []
    dh_bit3 = []

    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        M2 = list(M)
        M2[delta_word] ^= (1 << delta_bit)

        H1 = sha256_hash(M)
        H2 = sha256_hash(M2)

        dH = [H1[i] ^ H2[i] for i in range(8)]

        # Extract per-layer bits of δH
        b0 = tuple((dH[i] >> 0) & 1 for i in range(8))
        b1 = tuple((dH[i] >> 1) & 1 for i in range(8))
        b2 = tuple((dH[i] >> 2) & 1 for i in range(8))
        b3 = tuple((dH[i] >> 3) & 1 for i in range(8))

        dh_bit0.append(b0)
        dh_bit1.append(b1)
        dh_bit2.append(b2)
        dh_bit3.append(b3)

    # Correlation: P(δH bit-1 = 0 | δH bit-0 = 0)
    # "bit-k = 0" means all 8 registers have δH[reg][k] = 0
    bit0_zero = [i for i in range(N) if all(b == 0 for b in dh_bit0[i])]
    bit1_zero = [i for i in range(N) if all(b == 0 for b in dh_bit1[i])]
    bit0_and_1_zero = [i for i in range(N) if i in set(bit0_zero) and i in set(bit1_zero)]

    print(f"\n  δM = flip bit 0 of W[0], N = {N}")
    print(f"  P(δH bit-0 all zero): {len(bit0_zero)}/{N} = {len(bit0_zero)/N:.6f}")
    print(f"  P(δH bit-1 all zero): {len(bit1_zero)}/{N} = {len(bit1_zero)/N:.6f}")
    print(f"  P(both zero):         {len(bit0_and_1_zero)}/{N} = {len(bit0_and_1_zero)/N:.6f}")

    if len(bit0_zero) > 0:
        p_cond = len(bit0_and_1_zero) / len(bit0_zero)
        p_marginal = len(bit1_zero) / N
        print(f"  P(bit-1=0 | bit-0=0): {p_cond:.6f}")
        print(f"  P(bit-1=0) marginal:  {p_marginal:.6f}")
        print(f"  Ratio: {p_cond/p_marginal:.3f}×" if p_marginal > 0 else "  N/A")
    else:
        print(f"  No bit-0 zero cases found (expected: ~{N * 2**-8:.0f})")

    # Wider test: just look at INDIVIDUAL register correlations
    # For δH[0] (first hash word): are bit 0 and bit 1 of δH[0] correlated?
    print(f"\n  Per-register correlation of δH bits:")
    for reg in range(1):  # Just register 0 for clarity
        b0_vals = [(dh_bit0[i][reg]) for i in range(N)]
        b1_vals = [(dh_bit1[i][reg]) for i in range(N)]
        b2_vals = [(dh_bit2[i][reg]) for i in range(N)]
        b3_vals = [(dh_bit3[i][reg]) for i in range(N)]

        # P(b1=1 | b0=1) vs P(b1=1 | b0=0)
        b0_is_1 = [i for i in range(N) if b0_vals[i] == 1]
        b0_is_0 = [i for i in range(N) if b0_vals[i] == 0]

        p_b1_given_b0_1 = sum(b1_vals[i] for i in b0_is_1) / len(b0_is_1) if b0_is_1 else 0
        p_b1_given_b0_0 = sum(b1_vals[i] for i in b0_is_0) / len(b0_is_0) if b0_is_0 else 0

        print(f"    H[{reg}]: P(δb1=1|δb0=1)={p_b1_given_b0_1:.4f}, P(δb1=1|δb0=0)={p_b1_given_b0_0:.4f}")
        print(f"           Difference: {abs(p_b1_given_b0_1 - p_b1_given_b0_0):.4f}")

        # Same for b0 vs b2, b0 vs b3
        p_b2_given_b0_1 = sum(b2_vals[i] for i in b0_is_1) / len(b0_is_1) if b0_is_1 else 0
        p_b2_given_b0_0 = sum(b2_vals[i] for i in b0_is_0) / len(b0_is_0) if b0_is_0 else 0
        p_b3_given_b0_1 = sum(b3_vals[i] for i in b0_is_1) / len(b0_is_1) if b0_is_1 else 0
        p_b3_given_b0_0 = sum(b3_vals[i] for i in b0_is_0) / len(b0_is_0) if b0_is_0 else 0

        print(f"           P(δb2=1|δb0=1)={p_b2_given_b0_1:.4f}, P(δb2=1|δb0=0)={p_b2_given_b0_0:.4f}, diff={abs(p_b2_given_b0_1-p_b2_given_b0_0):.4f}")
        print(f"           P(δb3=1|δb0=1)={p_b3_given_b0_1:.4f}, P(δb3=1|δb0=0)={p_b3_given_b0_0:.4f}, diff={abs(p_b3_given_b0_1-p_b3_given_b0_0):.4f}")

    # DIFFERENT δM: try several
    print(f"\n  Testing multiple δM values:")
    for dw, db in [(0, 0), (0, 15), (0, 31), (7, 0), (15, 0)]:
        b0b1_agree = 0
        total = 0
        for seed in range(5000):
            rng = random.Random(seed + dw*100000 + db*1000)
            M = [rng.randint(0, MASK) for _ in range(16)]
            M2 = list(M); M2[dw] ^= (1 << db)
            H1 = sha256_hash(M)
            H2 = sha256_hash(M2)
            dH0 = H1[0] ^ H2[0]
            # Does bit 0 of δH predict bit 1?
            b0 = dH0 & 1
            b1 = (dH0 >> 1) & 1
            if b0 == b1:
                b0b1_agree += 1
            total += 1

        rate = b0b1_agree / total
        print(f"    δM = W[{dw}] bit {db}: P(δb0 = δb1) = {rate:.4f} (random = 0.5000)")


if __name__ == "__main__":
    experiment_structured_layers()
