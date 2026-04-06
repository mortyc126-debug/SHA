"""
Gröbner basis / XL test on mini-BTE systems.

Key question: is D_reg (degree of regularity) lower for structured
BTE systems than for random quadratic GF(2) systems?
If yes → Gröbner can exploit structure → potential sub-birthday attack.
"""

import random
from qt_solver.sha256_traced import get_bit
from qt_solver.gf2 import gf2_gaussian_eliminate


def mini_bte_round(state, w_r, k_r, n_bits, rots0, rots1):
    """One round of mini BTE."""
    mask = (1 << n_bits) - 1

    def rotr(x, r):
        return ((x >> r) | (x << (n_bits - r))) & mask

    def sig0(x):
        v = 0
        for r in rots0:
            v ^= rotr(x, r % n_bits)
        return v

    def sig1(x):
        v = 0
        for r in rots1:
            v ^= rotr(x, r % n_bits)
        return v

    a, b, c, d, e, f, g, h = state
    t1 = (h + sig1(e) + ((e & f) ^ (~e & g) & mask) + k_r + w_r) & mask
    t2 = (sig0(a) + ((a & b) ^ (a & c) ^ (b & c))) & mask
    return ((t1 + t2) & mask, a, b, c, (d + t1) & mask, e, f, g)


def mini_bte_hash_bits(msg_bits, n_bits=4, n_msg_words=4, R=8):
    """
    Mini BTE hash returning specific output bits.
    msg_bits: list of 0/1 values, length = n_msg_words * n_bits.
    Returns: list of 0/1 values (hash bits).
    """
    mask = (1 << n_bits) - 1
    # Reconstruct words from bits
    msg = []
    for w in range(n_msg_words):
        val = 0
        for b in range(n_bits):
            if msg_bits[w * n_bits + b]:
                val |= (1 << b)
        msg.append(val)

    # Simple schedule
    w_list = list(msg)
    while len(w_list) < R:
        i = len(w_list)
        w_list.append((w_list[i-2] ^ w_list[max(0, i-n_msg_words)]) & mask)

    iv = [(0x6a >> (8 - n_bits)) & mask] * 8
    state = tuple(iv)
    k_vals = [(i * 0x9e + 0x37) & mask for i in range(R)]

    for r in range(R):
        state = mini_bte_round(state, w_list[r], k_vals[r], n_bits, [1], [2])

    # Output hash bits
    hash_bits = []
    for reg in range(8):
        val = (state[reg] + iv[reg]) & mask
        for b in range(n_bits):
            hash_bits.append((val >> b) & 1)

    return hash_bits


def build_gf2_quadratic_system(n_bits=3, n_msg_words=3, R=4):
    """
    Build the quadratic GF(2) system for mini BTE preimage.

    Returns (polys, n_vars, n_eqs) where polys is list of sets of frozensets.
    Each poly = set of monomials. Monomial = frozenset of var indices.
    """
    n_vars = n_msg_words * n_bits
    n_hash = 8 * n_bits  # 8 registers × n_bits

    rng = random.Random(42)
    msg = [rng.randint(0, (1 << n_bits) - 1) for _ in range(n_msg_words)]
    msg_bits = []
    for w in range(n_msg_words):
        for b in range(n_bits):
            msg_bits.append((msg[w] >> b) & 1)

    target = mini_bte_hash_bits(msg_bits, n_bits, n_msg_words, R)

    # Build system by evaluating all 2^n_vars inputs (brute force ANF)
    # For each hash bit h: compute ANF via Möbius transform
    polys = []

    # Truth table for each hash bit
    for hb in range(min(n_hash, n_vars)):  # only need n_vars equations
        tt = [0] * (1 << n_vars)
        for x in range(1 << n_vars):
            bits = [(x >> i) & 1 for i in range(n_vars)]
            h = mini_bte_hash_bits(bits, n_bits, n_msg_words, R)
            tt[x] = h[hb] ^ target[hb]  # Want h[hb] = target[hb]

        # Möbius transform → ANF
        anf = list(tt)
        for i in range(n_vars):
            for x in range(1 << n_vars):
                if (x >> i) & 1:
                    anf[x] ^= anf[x ^ (1 << i)]

        # Extract monomials
        poly = set()
        for x in range(1 << n_vars):
            if anf[x]:
                mono = frozenset(i for i in range(n_vars) if (x >> i) & 1)
                poly.add(mono)

        polys.append(poly)

    return polys, n_vars, len(polys)


def xl_analyze(polys, n_vars, max_degree=4):
    """
    XL analysis: measure rank growth with degree.
    At each degree d, multiply equations by all monomials of degree ≤ d-2,
    linearize, and compute rank.
    """
    results = {}

    # Collect all monomials up to max_degree
    def gen_monomials(n, max_d):
        """Generate all monomials up to degree max_d."""
        monos = [frozenset()]  # constant
        for d in range(1, max_d + 1):
            from itertools import combinations
            for combo in combinations(range(n), d):
                monos.append(frozenset(combo))
        return monos

    all_monos = gen_monomials(n_vars, max_degree)
    mono_to_idx = {m: i for i, m in enumerate(all_monos)}
    n_monos = len(all_monos)

    for target_deg in range(2, max_degree + 1):
        # Multiplier monomials: degree ≤ target_deg - 2
        multipliers = [m for m in all_monos if len(m) <= target_deg - 2]

        rows = []
        for poly in polys:
            max_poly_deg = max((len(m) for m in poly), default=0)
            for mult in multipliers:
                if max_poly_deg + len(mult) > target_deg:
                    continue

                # Multiply poly by mult
                new_poly = set()
                for mono in poly:
                    product = mono.symmetric_difference(mult)  # XOR in GF(2)
                    if len(product) <= max_degree:
                        if product in new_poly:
                            new_poly.remove(product)
                        else:
                            new_poly.add(product)

                # Convert to row vector
                row = 0
                for mono in new_poly:
                    if mono in mono_to_idx:
                        row |= (1 << mono_to_idx[mono])
                if row:
                    rows.append(row)

        if rows:
            ech, pivots = gf2_gaussian_eliminate(list(rows), n_monos)
            rank = len(pivots)
        else:
            rank = 0

        results[target_deg] = {
            'rank': rank,
            'n_monos': n_monos,
            'n_rows': len(rows),
            'n_multipliers': len(multipliers),
        }

    return results


def random_quadratic_system(n_vars, n_eqs, seed=99):
    """Generate random quadratic GF(2) system."""
    rng = random.Random(seed)
    polys = []
    for _ in range(n_eqs):
        poly = set()
        # Random quadratic terms
        for i in range(n_vars):
            for j in range(i + 1, n_vars):
                if rng.randint(0, 1):
                    poly.add(frozenset([i, j]))
        # Random linear terms
        for i in range(n_vars):
            if rng.randint(0, 1):
                mono = frozenset([i])
                if mono in poly:
                    poly.remove(mono)
                else:
                    poly.add(mono)
        # Random constant
        if rng.randint(0, 1):
            mono = frozenset()
            if mono in poly:
                poly.remove(mono)
            else:
                poly.add(mono)
        polys.append(poly)
    return polys


def measure_d_reg():
    """Compare D_reg for BTE vs random systems."""
    print("=" * 65)
    print("Gröbner / XL: D_reg Comparison (BTE vs Random)")
    print("=" * 65)

    configs = [
        (3, 3, 4),  # n_bits=3, n_msg=3, R=4 → 9 vars
        (3, 3, 6),  # R=6
        (3, 4, 4),  # 12 vars
        (4, 3, 4),  # 12 vars
    ]

    for n_bits, n_msg, R in configs:
        n_vars = n_bits * n_msg
        print(f"\n--- n_bits={n_bits}, n_msg={n_msg}, R={R} (n_vars={n_vars}) ---")

        # BTE system
        print(f"  Building BTE system...", end=" ", flush=True)
        polys, nv, ne = build_gf2_quadratic_system(n_bits, n_msg, R)
        print(f"{ne} equations")

        # Analyze degree distribution
        degrees = []
        for p in polys:
            d = max((len(m) for m in p), default=0)
            degrees.append(d)
        print(f"  Max poly degree: {max(degrees)}, avg: {sum(degrees)/len(degrees):.1f}")

        # XL analysis
        print(f"  XL analysis (BTE):")
        bte_xl = xl_analyze(polys, nv, max_degree=min(4, nv))
        for d, info in sorted(bte_xl.items()):
            solved = "SOLVED" if info['rank'] >= info['n_monos'] else ""
            print(f"    deg={d}: rank={info['rank']}/{info['n_monos']} "
                  f"({info['n_rows']} rows) {solved}")

        # Random system of same size
        rand_polys = random_quadratic_system(nv, ne)
        print(f"  XL analysis (RANDOM):")
        rand_xl = xl_analyze(rand_polys, nv, max_degree=min(4, nv))
        for d, info in sorted(rand_xl.items()):
            solved = "SOLVED" if info['rank'] >= info['n_monos'] else ""
            print(f"    deg={d}: rank={info['rank']}/{info['n_monos']} "
                  f"({info['n_rows']} rows) {solved}")

        # Compare
        for d in sorted(bte_xl.keys()):
            if d in rand_xl:
                bte_r = bte_xl[d]['rank']
                rand_r = rand_xl[d]['rank']
                diff = bte_r - rand_r
                print(f"    deg={d}: BTE-Random = {diff:+d} "
                      f"({'BTE denser' if diff > 0 else 'Random denser' if diff < 0 else 'equal'})")


if __name__ == '__main__':
    measure_d_reg()
    print("\n" + "=" * 65)
    print("CONCLUSION")
    print("=" * 65)
    print("""
  If BTE rank grows FASTER than random at low degree:
    → BTE has LOWER D_reg → Gröbner exploits structure
    → Potential sub-birthday attack
  If BTE rank ≈ random:
    → No structural advantage for Gröbner
    → Birthday optimal confirmed algebraically
    """)
