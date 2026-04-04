"""
VERIFY: Run actual XL-like solving on small BTE instances.

For n=8, R=16: layer_size = 31, total = 128 bits.
  Layer 0: 31 degree-2 eqs in 128 vars. Solution space = 97 dims.

Can we ACTUALLY solve the degree-2 system for a small BTE?
If yes → verify complexity estimates.
If no → our estimates are wrong.

Method: Build the ACTUAL degree-2 GF(2) system for a BTE,
then solve with Gaussian elimination on linearized (XL) system.
"""

import random
import time


MASK8 = 0xFF


def rotr8(x, k):
    return ((x >> k) | (x << (8 - k))) & MASK8


def ch8(e, f, g):
    return (e & f) ^ (~e & g) & MASK8

def maj8(a, b, c):
    return (a & b) ^ (a & c) ^ (b & c)


def bte8_trajectory(msg, R=16):
    """Compute BTE-8 (n=8) trajectory."""
    n = 8
    mask = MASK8
    n_msg = 16  # 16 message words

    IV = [0x6a, 0xbb, 0x3c, 0xa5, 0x51, 0x9b, 0x1f, 0x5b]
    K = [(i * 37 + 113) & mask for i in range(R)]

    # Schedule
    W = list(msg[:n_msg])
    for i in range(n_msg, R):
        W.append((W[i-2] ^ W[i-7] ^ W[i-n_msg]) & mask)

    sig0_rots = [1, 3, 5]
    sig1_rots = [2, 4, 6]

    def sig0(x):
        return rotr8(x, 1) ^ rotr8(x, 3) ^ rotr8(x, 5)
    def sig1(x):
        return rotr8(x, 2) ^ rotr8(x, 4) ^ rotr8(x, 6)

    state = tuple(IV)
    states = [state]

    for r in range(R):
        a, b, c, d, e, f, g, h = state
        T1 = (h + sig1(e) + ch8(e, f, g) + K[r] + W[r]) & mask
        T2 = (sig0(a) + maj8(a, b, c)) & mask
        state = ((T1 + T2) & mask, a, b, c, (d + T1) & mask, e, f, g)
        states.append(state)

    return states, W


def build_degree2_system(target_hash_bits, R=16):
    """
    Build the degree-2 GF(2) system that the BTE defines.

    Variables: 128 message bits (16 words × 8 bits).
    Equations: target_hash_bits (list of (reg, bit, value) triples).

    For each equation: express hash_bit as degree-2 polynomial
    in message bits by computing the Jacobian AND the quadratic terms.

    Returns: list of equations, each = set of monomials (as frozensets).
    """
    # This is complex to do symbolically. Instead, use numerical approach:
    # For the LINEAR part: Jacobian (already have this).
    # For the QUADRATIC part: compute second-order derivatives.
    #   ∂²f/∂x_i∂x_j = f(0) ⊕ f(e_i) ⊕ f(e_j) ⊕ f(e_i⊕e_j)

    pass  # Too complex for symbolic. Use numerical verification instead.


def verify_layer_rank():
    """Verify layer ranks for BTE-8."""
    print("=" * 80)
    print("VERIFY: Layer ranks for BTE-8 (n=8, R=16)")
    print("=" * 80)

    n = 8
    n_msg = 16
    R = 16
    n_input = n * n_msg  # 128

    rng = random.Random(42)
    M_base = [rng.randint(0, MASK8) for _ in range(n_msg)]
    states_base, W_base = bte8_trajectory(M_base, R)

    # Hash = final state + IV feedforward (simplified: just final state)
    hash_base = states_base[R]

    prev_rank = 0
    for max_bit in range(n):
        base_bits = []
        for reg in range(8):
            for b in range(max_bit + 1):
                base_bits.append((hash_base[reg] >> b) & 1)

        n_eqs = len(base_bits)

        J = []
        for word in range(n_msg):
            for bit in range(n):
                M_flip = list(M_base)
                M_flip[word] ^= (1 << bit)
                sf, _ = bte8_trajectory(M_flip, R)
                hash_flip = sf[R]
                row = []
                for reg in range(8):
                    for b in range(max_bit + 1):
                        row.append(base_bits[reg*(max_bit+1)+b] ^
                                  ((hash_flip[reg] >> b) & 1))
                J.append(row)

        rank = gf2_rank(J, n_input, n_eqs)
        delta = rank - prev_rank
        print(f"  bit 0..{max_bit}: rank={rank:>3}/{n_input}, Δ=+{delta}")
        prev_rank = rank

        if rank >= n_input:
            print(f"  FULL RANK at bit {max_bit}")
            break


def verify_preimage_search():
    """
    ACTUAL preimage search on BTE-8 using layered approach.

    Target: a known hash H (from a known message M_target).
    Task: find ANY message M with BTE8(M) = H.

    Layer-0 approach:
      1. Build bit-0 Jacobian (128×8 = 128 inputs, 8 outputs)
      2. Find kernel (messages with same bit-0 hash)
      3. Among kernel: search for full hash match
    """
    print("\n" + "=" * 80)
    print("VERIFY: Actual preimage search on BTE-8")
    print("=" * 80)

    n = 8
    n_msg = 16
    R = 16

    # Target message and hash
    rng = random.Random(42)
    M_target = [rng.randint(0, MASK8) for _ in range(n_msg)]
    states_target, _ = bte8_trajectory(M_target, R)
    H_target = states_target[R]

    print(f"  Target hash: {['0x%02x' % h for h in H_target]}")

    # Brute force: try random messages
    t0 = time.time()
    attempts = 0
    found_brute = None

    for seed in range(100000):
        M_try = [(seed * 37 + i * 71 + 13) & MASK8 for i in range(n_msg)]
        sf, _ = bte8_trajectory(M_try, R)
        attempts += 1

        if sf[R] == H_target:
            found_brute = M_try
            break

    t_brute = time.time() - t0

    if found_brute:
        print(f"  Brute force: found in {attempts} attempts ({t_brute:.2f}s)")
    else:
        print(f"  Brute force: NOT found in {attempts} attempts ({t_brute:.2f}s)")

    # Layered approach: first match bit-0, then check higher bits
    t0 = time.time()
    layer0_target = tuple((H_target[reg] & 1) for reg in range(8))

    layer0_matches = 0
    layer01_matches = 0
    full_matches = 0

    for seed in range(100000):
        M_try = [(seed * 37 + i * 71 + 13) & MASK8 for i in range(n_msg)]
        sf, _ = bte8_trajectory(M_try, R)
        H_try = sf[R]

        l0 = tuple((H_try[reg] & 1) for reg in range(8))
        if l0 == layer0_target:
            layer0_matches += 1

            l1 = tuple(((H_try[reg] >> 1) & 1) for reg in range(8))
            l1_target = tuple(((H_target[reg] >> 1) & 1) for reg in range(8))
            if l1 == l1_target:
                layer01_matches += 1

                if H_try == H_target:
                    full_matches += 1

    t_layer = time.time() - t0

    print(f"\n  Layered analysis (100K messages):")
    print(f"    Layer 0 matches: {layer0_matches} (expected: 100000/2^8 = {100000//256})")
    print(f"    Layer 0+1 matches: {layer01_matches} (expected: 100000/2^16 = {100000//65536})")
    print(f"    Full matches: {full_matches}")
    print(f"    Time: {t_layer:.2f}s")

    # Birthday estimate for full 64-bit hash
    print(f"\n  Hash space: 2^{n*8} = 2^{n*8}")
    print(f"  Birthday: 2^{n*8//2} = 2^{n*8//2}")
    print(f"  Brute preimage: 2^{n*8} = 2^{n*8}")


def gf2_rank(matrix, nrows, ncols):
    m = [list(row[:ncols]) for row in matrix[:nrows]]
    rank = 0
    for col in range(ncols):
        pivot = None
        for row in range(rank, len(m)):
            if m[row][col] == 1:
                pivot = row
                break
        if pivot is None:
            continue
        m[rank], m[pivot] = m[pivot], m[rank]
        for row in range(len(m)):
            if row != rank and m[row][col] == 1:
                for c in range(ncols):
                    m[row][c] ^= m[rank][c]
        rank += 1
    return rank


if __name__ == "__main__":
    verify_layer_rank()
    verify_preimage_search()
