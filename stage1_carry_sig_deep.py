"""
Candidate 2 showed max bit bias = 0.5000 in carry signature.
This is HUGE — some bit of the rolling carry XOR is DETERMINISTIC.

Let's find which bit and why.
"""

import random
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, schedule
)


def get_carry_pattern(x, y):
    """Return 32-bit carry-into pattern for x + y."""
    carry = 0
    pattern = 0
    for i in range(32):
        xi = (x >> i) & 1
        yi = (y >> i) & 1
        if carry:
            pattern |= (1 << i)
        carry = (xi & yi) | (xi & carry) | (yi & carry)
    return pattern


def deep_carry_signature():
    N = 1000

    print("=" * 80)
    print("DEEP ANALYSIS: Which bits of carry signature are biased?")
    print("=" * 80)

    # Collect per-bit statistics
    bit_ones = [0] * 32

    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        states, W = sha256_round_trace(M)

        cs = 0
        for r in range(64):
            a, b, c, d, e, f, g, h = states[r]
            T1 = (h + Sig1(e) + Ch(e, f, g) + K[r] + W[r]) & MASK
            T2 = (Sig0(a) + Maj(a, b, c)) & MASK
            carry_pattern = get_carry_pattern(T1, T2)
            cs ^= carry_pattern

        for bit in range(32):
            if (cs >> bit) & 1:
                bit_ones[bit] += 1

    print(f"\nPer-bit P(CS[bit]=1) for rolling XOR of 64 carry patterns:")
    for bit in range(32):
        p = bit_ones[bit] / N
        marker = " *** BIASED ***" if abs(p - 0.5) > 0.1 else ""
        print(f"  bit {bit:>2}: P={p:.3f}{marker}")

    # Find the most biased bit
    most_biased = max(range(32), key=lambda b: abs(bit_ones[b]/N - 0.5))
    p_biased = bit_ones[most_biased] / N
    print(f"\nMost biased: bit {most_biased}, P={p_biased:.3f}")

    # Now: is this bias from the carry signature, or from the ADDITION itself?
    # Let's also check: XOR of T1+T2 RESULTS (not just carries)
    print(f"\n--- Control: XOR of a[r+1] values (not carry patterns) ---")
    bit_ones_control = [0] * 32
    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        states, _ = sha256_round_trace(M)
        xor_all = 0
        for r in range(1, 65):
            xor_all ^= states[r][0]
        for bit in range(32):
            if (xor_all >> bit) & 1:
                bit_ones_control[bit] += 1

    for bit in range(32):
        p = bit_ones_control[bit] / N
        if abs(p - 0.5) > 0.1:
            print(f"  bit {bit}: P={p:.3f} *** BIASED ***")

    # Let's check: is the bias in carry signature JUST bit 0?
    # Bit 0 of carry = always 0 (no carry into bit 0).
    # XOR of 64 zeros = 0. So bit 0 of CS = always 0.
    print(f"\n--- Bit 0 of carry pattern ---")
    print(f"  Carry into bit 0 = always 0 (by definition)")
    print(f"  XOR of 64 zeros = 0")
    print(f"  So CS[0] = 0 always → P(CS[0]=1) = {bit_ones[0]/N:.3f}")
    print(f"  This is TRIVIAL — not a real signal.")

    # Remove bit 0 and recheck
    print(f"\n--- After removing trivial bit 0: ---")
    real_biased = []
    for bit in range(1, 32):
        p = bit_ones[bit] / N
        if abs(p - 0.5) > 0.05:
            real_biased.append((bit, p))
            print(f"  bit {bit}: P={p:.3f}")

    if not real_biased:
        print(f"  No non-trivial biased bits found.")
        print(f"  The 0.5000 max bias was from bit 0 (trivial: carry-in to bit 0 = always 0)")

    # Also check: what about PER-ADDITION carry patterns, not just T1+T2?
    # There are ~6 additions per round. Let's check each.
    print(f"\n--- Per-addition carry analysis (round 20, 500 messages) ---")
    N2 = 500
    additions = ['h+Sig1(e)', '+Ch', '+K_W', 'Sig0+Maj', 'd+T1', 'T1+T2']

    for add_idx, add_name in enumerate(additions):
        bit0_count = 0
        for seed in range(N2):
            rng = random.Random(seed)
            M = [rng.randint(0, MASK) for _ in range(16)]
            states, W = sha256_round_trace(M)
            a, b, c, d, e, f, g, h = states[20]

            s1e = Sig1(e)
            che = Ch(e, f, g)
            kw = (K[20] + W[20]) & MASK
            s0a = Sig0(a)
            maja = Maj(a, b, c)

            if add_idx == 0:
                cp = get_carry_pattern(h, s1e)
            elif add_idx == 1:
                cp = get_carry_pattern((h + s1e) & MASK, che)
            elif add_idx == 2:
                cp = get_carry_pattern((h + s1e + che) & MASK, kw)
            elif add_idx == 3:
                cp = get_carry_pattern(s0a, maja)
            elif add_idx == 4:
                T1 = (h + s1e + che + kw) & MASK
                cp = get_carry_pattern(d, T1)
            else:
                T1 = (h + s1e + che + kw) & MASK
                T2 = (s0a + maja) & MASK
                cp = get_carry_pattern(T1, T2)

            # Count carry at bit 1 (first non-trivial carry)
            if (cp >> 1) & 1:
                bit0_count += 1

        print(f"  {add_name:>12}: P(carry into bit 1) = {bit0_count/N2:.3f}")


if __name__ == "__main__":
    deep_carry_signature()
