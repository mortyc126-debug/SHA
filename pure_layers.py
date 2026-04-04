"""
PURE LAYERS: What determines the number of "full" layers?

We have: n=8→3, n=16→8, n=32→4 pure layers.
But rotations were DIFFERENT for each n.

Test: fix n=8, vary rotations. Does pure_layers change?
If YES → it's rotation-dependent.
If NO → it's n-dependent.
"""

import random


def rotr_n(x, k, n):
    mask = (1 << n) - 1
    return ((x >> k) | (x << (n - k))) & mask


def ch_fn(e, f, g, mask):
    return (e & f) ^ (~e & g) & mask

def maj_fn(a, b, c, mask):
    return (a & b) ^ (a & c) ^ (b & c)


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


def measure_pure_layers(n, sig_rots, R=None):
    """Count pure layers for given BTE config."""
    if R is None:
        R = 2 * n
    mask = (1 << n) - 1
    n_msg = 2 * n

    rng_iv = random.Random(12345)
    IV = [rng_iv.randint(0, mask) for _ in range(8)]
    K = [random.Random(i + 99999).randint(0, mask) for i in range(R)]

    def sig0(x):
        r = 0
        for rot in sig_rots[:len(sig_rots)//2] or sig_rots[:2]:
            r ^= rotr_n(x, rot, n)
        return r

    def sig1(x):
        r = 0
        for rot in sig_rots[len(sig_rots)//2:] or sig_rots[2:]:
            r ^= rotr_n(x, rot, n)
        return r

    def compute_trajectory(msg):
        W = list(msg[:n_msg])
        for i in range(n_msg, R):
            W.append((W[i-2] ^ W[i-7 if i>=7 else 0] ^ W[i-n_msg]) & mask)

        state = tuple(IV)
        states = [state]
        for r in range(R):
            a, b, c, d, e, f, g, h = state
            T1 = (h + sig1(e) + ch_fn(e, f, g, mask) + K[r] + W[r]) & mask
            T2 = (sig0(a) + maj_fn(a, b, c, mask)) & mask
            state = ((T1 + T2) & mask, a, b, c, (d + T1) & mask, e, f, g)
            states.append(state)
        return states

    rng = random.Random(42)
    M_base = [rng.randint(0, mask) for _ in range(n_msg)]
    states_base = compute_trajectory(M_base)
    n_input = n * n_msg

    expected_layer = 2 * R - 1
    pure_count = 0
    prev_rank = 0

    for max_bit in range(n):
        base_bits = []
        for r in range(R):
            for reg in range(8):
                for b in range(max_bit + 1):
                    base_bits.append((states_base[r+1][reg] >> b) & 1)

        n_eqs = len(base_bits)
        J = []
        for word in range(n_msg):
            for bit in range(n):
                M_flip = list(M_base)
                M_flip[word] ^= (1 << bit)
                sf = compute_trajectory(M_flip)
                row = []
                for r in range(R):
                    for reg in range(8):
                        for b in range(max_bit + 1):
                            row.append(base_bits[r*8*(max_bit+1)+reg*(max_bit+1)+b] ^
                                      ((sf[r+1][reg] >> b) & 1))
                J.append(row)

        rank = gf2_rank(J, n_input, n_eqs)
        delta = rank - prev_rank

        if delta == expected_layer:
            pure_count += 1
        else:
            break  # First non-pure layer

        prev_rank = rank
        if rank >= n_input:
            break

    return pure_count


def experiment_vary_rotations():
    """Fix n=8, vary rotation sets. Measure pure layers."""
    print("=" * 80)
    print("PURE LAYERS: Vary rotations at fixed n=8")
    print("=" * 80)

    n = 8
    # Different rotation sets (all for n=8)
    configs = [
        ("SHA-like {1,3,5,2,4,6}", [1,3,5,2,4,6]),
        ("Alt {1,2,4,3,5,7}", [1,2,4,3,5,7]),
        ("Min {1,4,2,5}", [1,4,2,5]),
        ("Wide {1,2,3,4,5,6}", [1,2,3,4,5,6]),
        ("Narrow {1,2,5,6}", [1,2,5,6]),
        ("Single {3,5}", [3,5]),
        ("Coprime {1,3,2,5}", [1,3,2,5]),
        ("Powers {1,2,4,3,5,7}", [1,2,4,3,5,7]),
    ]

    for name, rots in configs:
        pure = measure_pure_layers(n, rots)
        print(f"  {name:>30}: {pure} pure layers")


def experiment_vary_n():
    """Fix rotation RATIOS, vary n. Measure pure layers."""
    print("\n" + "=" * 80)
    print("PURE LAYERS: Vary n with proportional rotations")
    print("=" * 80)

    for n in [4, 6, 8, 10, 12, 16]:
        # Proportional rotations: ~{n/16, n/4, n/2-3} for sig0, ~{n/4-1, n/3, n-3} for sig1
        r1 = max(1, n//8)
        r2 = max(2, n//4)
        r3 = max(3, n//2 - 1)
        r4 = max(1, n//4 - 1)
        r5 = max(2, n//3)
        r6 = max(3, n - 3)
        rots = list(set([r1, r2, r3, r4, r5, r6]))  # Remove duplicates
        if len(rots) < 4:
            rots = [1, n//4, n//2, n-1]

        pure = measure_pure_layers(n, rots)
        expected_layer = 2 * (2*n) - 1  # R = 2n
        total_bits = n * 2 * n
        print(f"  n={n:>2} (R={2*n:>2}): {pure} pure layers, "
              f"total={total_bits}, layer_size={expected_layer}, "
              f"pure×layer={pure*expected_layer}")


if __name__ == "__main__":
    experiment_vary_rotations()
    experiment_vary_n()
