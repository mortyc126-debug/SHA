"""
BTE CLASS THEORY: Study the class of ALL BTE, not just SHA-256.

SHA-256 = BTE(n=32, σ=SHA-rotations, V=MAJ, schedule=SHA-schedule).

What if we vary these parameters and find UNIVERSAL properties?

Parameters to vary:
  n: bit-width (8, 16, 32, 64)
  σ: rotation constants
  V: carry function (MAJ, AND, XOR, OR)
  schedule: constraint linking rounds

For each configuration: measure the same structural constants
we found for SHA-256 (4 layers, 127 per layer, etc.).

If these constants are UNIVERSAL across configurations → they're
properties of the BTE CLASS, not just SHA-256.
If they VARY → they depend on parameters → we can find the FORMULA.
"""

import random


MASK8 = 0xFF
MASK16 = 0xFFFF
MASK32 = 0xFFFFFFFF


def rotr_n(x, k, n):
    """Rotate right by k in n-bit word."""
    mask = (1 << n) - 1
    return ((x >> k) | (x << (n - k))) & mask


def make_bte_round(n, sig0_rots, sig1_rots, ch_fn, maj_fn, K_val, W_val):
    """
    One round of a generic BTE hash function.

    State: (a, b, c, d, e, f, g, h) — 8 registers of n bits each.
    Returns new state.
    """
    mask = (1 << n) - 1

    def sig0(x):
        result = 0
        for r in sig0_rots:
            result ^= rotr_n(x, r, n)
        return result

    def sig1(x):
        result = 0
        for r in sig1_rots:
            result ^= rotr_n(x, r, n)
        return result

    def round_fn(state):
        a, b, c, d, e, f, g, h = state
        T1 = (h + sig1(e) + ch_fn(e, f, g, mask) + K_val + W_val) & mask
        T2 = (sig0(a) + maj_fn(a, b, c, mask)) & mask
        a_new = (T1 + T2) & mask
        e_new = (d + T1) & mask
        return (a_new, a, b, c, e_new, e, f, g)

    return round_fn


def ch_standard(e, f, g, mask):
    return (e & f) ^ (~e & g) & mask

def maj_standard(a, b, c, mask):
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


def measure_layer_structure(n, sig0_rots, sig1_rots, R=None):
    """
    For a BTE with given parameters, measure the bit-layer unfolding.
    Returns: list of cumulative ranks for bits 0, 0-1, 0-2, ..., 0-(n-1).
    """
    if R is None:
        R = 2 * n  # enough rounds

    mask = (1 << n) - 1
    n_msg = 2 * n  # message words (like 16 for n=32)
    state_regs = 8

    # IV (random but fixed)
    rng_iv = random.Random(12345)
    IV = [rng_iv.randint(0, mask) for _ in range(state_regs)]

    # K constants
    K = [random.Random(i + 99999).randint(0, mask) for i in range(R)]

    def compute_trajectory(msg):
        """Compute full state trajectory for a message."""
        # Simple schedule: W[i] = msg[i] for i < n_msg, then linear extension
        W = list(msg[:n_msg])
        for i in range(n_msg, R):
            # Simple linear schedule (XOR of earlier words)
            val = W[i - 2] ^ W[i - 7 if i >= 7 else 0] ^ W[i - n_msg]
            W.append(val & mask)

        state = tuple(IV)
        states = [state]

        for r in range(R):
            round_fn = make_bte_round(n, sig0_rots, sig1_rots,
                                       ch_standard, maj_standard, K[r], W[r])
            state = round_fn(state)
            states.append(state)

        return states

    # Base message
    rng = random.Random(42)
    M_base = [rng.randint(0, mask) for _ in range(n_msg)]
    states_base = compute_trajectory(M_base)

    # For each bit level, compute rank
    n_input = n * n_msg  # total message bits
    ranks = []

    for max_bit in range(n):
        # Base bits for this level
        base_bits = []
        for r in range(R):
            for reg in range(state_regs):
                for b in range(max_bit + 1):
                    base_bits.append((states_base[r + 1][reg] >> b) & 1)

        n_eqs = len(base_bits)

        # Jacobian
        J = []
        for word in range(n_msg):
            for bit in range(n):
                M_flip = list(M_base)
                M_flip[word] ^= (1 << bit)
                states_flip = compute_trajectory(M_flip)
                row = []
                for r in range(R):
                    for reg in range(state_regs):
                        for b in range(max_bit + 1):
                            row.append(base_bits[r * state_regs * (max_bit + 1) +
                                                  reg * (max_bit + 1) + b] ^
                                       ((states_flip[r + 1][reg] >> b) & 1))
                J.append(row)

        rank = gf2_rank(J, n_input, n_eqs)
        ranks.append(rank)

        if rank >= n_input:
            break

    return ranks


def experiment_class_theory():
    """
    Vary BTE parameters and measure layer structure.
    """
    print("=" * 80)
    print("BTE CLASS THEORY: Layer structure across different BTE instances")
    print("=" * 80)

    # Configuration: (name, n, sig0_rots, sig1_rots)
    configs = [
        # SHA-256 (scaled to n=8 for speed)
        ("SHA-like n=8", 8, [1, 3, 5], [2, 4, 6]),
        # Different rotations
        ("Alt-rot n=8", 8, [1, 2, 4], [3, 5, 7]),
        # Minimal rotations (just 2)
        ("Min-rot n=8", 8, [1, 4], [2, 5]),
        # Single rotation
        ("Single-rot n=8", 8, [3], [5]),
        # SHA-256 actual (n=16 for medium speed)
        ("SHA-like n=16", 16, [1, 7, 11], [3, 8, 13]),
    ]

    for name, n, s0, s1 in configs:
        R = min(2 * n, 32)  # Cap rounds for speed
        print(f"\n  {name} (n={n}, R={R}):")
        print(f"    Sig0 rotations: {s0}")
        print(f"    Sig1 rotations: {s1}")

        ranks = measure_layer_structure(n, s0, s1, R=R)

        n_input = n * 2 * n  # n_msg = 2n
        prev = 0
        for i, rank in enumerate(ranks):
            delta = rank - prev
            print(f"    bit 0..{i}: rank = {rank}/{n_input}, Δ = +{delta}")
            prev = rank
            if rank >= n_input:
                print(f"    FULL RANK at bit {i}")
                break

        # Compute the "structural constants"
        if len(ranks) >= 2:
            layer_size = ranks[0]  # First layer rank
            n_layers = 0
            for r in ranks:
                if r >= n_input:
                    break
                n_layers += 1

            # Formula check: rank = layer_size × (layer_index + 1)?
            predicted = [(i + 1) * layer_size for i in range(n_layers)]
            actual = ranks[:n_layers]
            match = all(abs(p - a) <= 2 for p, a in zip(predicted, actual))

            print(f"    Layer size = {layer_size}")
            print(f"    Layers to full rank = {n_layers}")
            print(f"    {n_input} = {n_layers} × {layer_size} + {n_input - n_layers * layer_size}")
            print(f"    Linear pattern: {'YES' if match else 'NO'}")

            # The KEY formula: layer_size = 2R - 1?
            expected_127 = 2 * R - 1
            print(f"    2R - 1 = {expected_127}, actual layer_size = {layer_size}")
            print(f"    Match: {'YES' if layer_size == expected_127 else 'NO'}")


if __name__ == "__main__":
    experiment_class_theory()
