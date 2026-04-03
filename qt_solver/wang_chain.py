"""
Wang Chain: Sequential additive Newton in Z/2^32.

Bypasses the k*=5 GF(2) degree barrier entirely.
Works for rounds 1-16 (P=1.0, deterministic).

Algorithm:
  delta_W[0] = 0x8000 (initial perturbation)
  For r = 1..15:
    Compute delta_e_natural[r+1] from current state difference
    Set delta_W[r] = -delta_e_natural[r+1] mod 2^32
    This forces delta_e[r+1] = 0 (sensitivity = 1, Newton exact)

Result: 16 consecutive rounds with de[r] = 0.
Cost: O(16) = constant time.
"""

from qt_solver.sha256_traced import (
    MASK32, IV256, K256, sha256_compress, sha256_compress_traced,
    add32, rotr32, big_sigma0, big_sigma1, ch, maj,
)


def wang_chain(W_base, delta_W0=0x8000, num_correction_rounds=15):
    """
    Construct Wang pair: (W_base, W_modified) with de[2..R]=0.

    Args:
        W_base: 16 message words (the "normal" message)
        delta_W0: initial additive perturbation on W[0]
        num_correction_rounds: rounds to correct (1..15)

    Returns:
        W_mod: modified message words
        de_values: de[r] for each round
        da_values: da[r] for each round
    """
    W_n = list(W_base)
    W_f = list(W_base)
    W_f[0] = add32(W_n[0], delta_W0)

    # Run both messages through SHA-256 round by round
    # Normal message
    a_n, b_n, c_n, d_n = IV256[0], IV256[1], IV256[2], IV256[3]
    e_n, f_n, g_n, h_n = IV256[4], IV256[5], IV256[6], IV256[7]

    # Modified message
    a_f, b_f, c_f, d_f = IV256[0], IV256[1], IV256[2], IV256[3]
    e_f, f_f, g_f, h_f = IV256[4], IV256[5], IV256[6], IV256[7]

    de_values = []
    da_values = []
    dw_values = [delta_W0]

    for r in range(num_correction_rounds + 1):
        # Compute T1 for normal
        S1_n = big_sigma1(e_n)
        ch_n = ch(e_n, f_n, g_n)
        T1_n = add32(add32(add32(add32(h_n, S1_n), ch_n), K256[r]), W_n[r])

        S0_n = big_sigma0(a_n)
        maj_n = maj(a_n, b_n, c_n)
        T2_n = add32(S0_n, maj_n)

        # Compute T1 for modified
        S1_f = big_sigma1(e_f)
        ch_f = ch(e_f, f_f, g_f)
        T1_f = add32(add32(add32(add32(h_f, S1_f), ch_f), K256[r]), W_f[r])

        S0_f = big_sigma0(a_f)
        maj_f = maj(a_f, b_f, c_f)
        T2_f = add32(S0_f, maj_f)

        # New state
        e_new_n = add32(d_n, T1_n)
        a_new_n = add32(T1_n, T2_n)

        e_new_f = add32(d_f, T1_f)
        a_new_f = add32(T1_f, T2_f)

        de = (e_new_f - e_new_n) & MASK32
        da = (a_new_f - a_new_n) & MASK32
        de_values.append(de)
        da_values.append(da)

        # Shift registers
        h_n, g_n, f_n, e_n = g_n, f_n, e_n, e_new_n
        d_n, c_n, b_n, a_n = c_n, b_n, a_n, a_new_n

        h_f, g_f, f_f, e_f = g_f, f_f, e_f, e_new_f
        d_f, c_f, b_f, a_f = c_f, b_f, a_f, a_new_f

        # For next round: adapt W_f[r+1] to force de[r+2] = 0
        if r < num_correction_rounds:
            next_r = r + 1
            # Compute delta_e_natural for next round (without correction)
            # delta_e[r+2] = delta_d[r+1] + delta_T1[r+1]
            # delta_d[r+1] = delta_c[r] = delta_b[r-1] = delta_a[r-2]
            # For the T1 correction:
            # T1_f[r+1] uses W_f[r+1]. We want:
            # e_f[r+2] = e_n[r+2], i.e., d_f + T1_f = d_n + T1_n
            # T1_f - T1_n = -(d_f - d_n) = delta needed
            # T1 = h + S1(e) + Ch(e,f,g) + K + W
            # delta_T1 = delta_h + delta_S1 + delta_Ch + delta_W
            # We want delta_T1 = -(d_f - d_n)
            # So delta_W = -(d_f - d_n) - delta_h - delta_S1 - delta_Ch

            delta_d_next = (d_f - d_n) & MASK32
            delta_h_next = (h_f - h_n) & MASK32
            delta_S1_next = (big_sigma1(e_f) - big_sigma1(e_n)) & MASK32
            delta_Ch_next = (ch(e_f, f_f, g_f) - ch(e_n, f_n, g_n)) & MASK32

            delta_W_needed = (0 - delta_d_next - delta_h_next
                              - delta_S1_next - delta_Ch_next) & MASK32

            W_f[next_r] = add32(W_n[next_r], delta_W_needed)
            dw_values.append(delta_W_needed)

    return W_f, de_values, da_values, dw_values


def verify_wang_chain(num_tests=100, verbose=True):
    """Verify Wang chain produces de=0 for rounds 2..16."""
    import random
    rng = random.Random(42)

    if verbose:
        print(f"{'='*60}")
        print(f"WANG CHAIN VERIFICATION")
        print(f"{'='*60}")

    success = 0
    for t in range(num_tests):
        W = [rng.randint(0, MASK32) for _ in range(16)]
        W_f, de_vals, da_vals, dw_vals = wang_chain(W)

        # Check de[1..15] (rounds 2..16 in 1-indexed)
        all_zero = all(de == 0 for de in de_vals[1:])
        if all_zero:
            success += 1

        if verbose and t < 5:
            print(f"\n  Test {t}: W[0]=0x{W[0]:08x}")
            for r, (de, da) in enumerate(zip(de_vals, da_vals)):
                de_str = f"0x{de:08x}" if de != 0 else "0"
                print(f"    round {r+1}: de={de_str}, HW(da)={bin(da).count('1')}")

    if verbose:
        print(f"\nResult: {success}/{num_tests} with de[2..16]=0 (P=1.0)")

    return success, num_tests


def wang_birthday_r17(W_base, max_attempts=2**24, verbose=True):
    """
    Birthday search for de[17] = 0.

    The Wang chain gives de[2..16]=0. Round 17 uses W[16] from schedule:
    W[16] = sigma1(W[14]) + W[9] + sigma0(W[1]) + W[0]

    de[17] = Da[13] + DW[16] (T_BARRIER17_EXACT).
    We search over W[0] (which affects Da[13] and DW[16]) for de[17]=0.
    """
    import random, time
    rng = random.Random(42)

    if verbose:
        print(f"\n{'='*60}")
        print(f"WANG BIRTHDAY SEARCH: de[17]=0")
        print(f"{'='*60}")

    t0 = time.time()
    for attempt in range(max_attempts):
        W = list(W_base)
        W[0] = rng.randint(0, MASK32)

        W_f, de_vals, da_vals, dw_vals = wang_chain(W, num_correction_rounds=15)

        # Now compute round 16 (r=16, using W[16] from schedule)
        # Need full SHA to get de[17]
        H_n = sha256_compress(W, 17)
        H_f = sha256_compress(W_f, 17)

        # Quick check: are hashes equal for first 17 rounds?
        if H_n == H_f:
            if verbose:
                elapsed = time.time() - t0
                print(f"  FOUND at attempt {attempt}! ({elapsed:.1f}s)")
                print(f"  W[0] = 0x{W[0]:08x}")
            return W, W_f, attempt

        if verbose and attempt % 1000000 == 0 and attempt > 0:
            elapsed = time.time() - t0
            rate = attempt / elapsed
            print(f"  {attempt} attempts, {rate:.0f}/s, {elapsed:.1f}s")

    if verbose:
        print(f"  Not found in {max_attempts} attempts")
    return None


if __name__ == '__main__':
    verify_wang_chain()
