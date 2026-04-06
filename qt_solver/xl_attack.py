"""
Direction 1: XL (eXtended Linearization) attack on R=16 SHA-256.

Idea: The Jacobian (degree-1) gives α-kernel = 0 at R=16.
But maybe degree-2 analysis reveals redundancy in the quadratic carry system.

XL approach:
1. Compute Hessian: H[i][j][k] = ∂²hash_k / (∂msg_i · ∂msg_j) in GF(2)
2. Build extended system: for each pair (i,j) of message bits, treat
   δm_i * δm_j as a new variable → linearize degree-2 relationships
3. Look for contradictions or rank deficiency in extended system

Key question: are the ~4680 quadratic carry equations at R=16
truly independent, or does XL find hidden redundancy?
"""

import random
from qt_solver.sha256_traced import MASK32, get_bit, sha256_compress
from qt_solver.gf2 import gf2_gaussian_eliminate


def compute_hessian_sample(R, msg, hash_bit_indices=None, msg_bit_indices=None):
    """
    Compute Hessian (second derivative) of SHA-256 hash bits wrt message bits.

    H[hb][mi][mj] = f(m) ⊕ f(m⊕ei) ⊕ f(m⊕ej) ⊕ f(m⊕ei⊕ej)

    This is the GF(2) second derivative. If H = 0, the function is affine
    in those directions. Non-zero H entries indicate quadratic interactions.

    To keep computation tractable, we sample specific hash/msg bit subsets.
    """
    if hash_bit_indices is None:
        hash_bit_indices = list(range(256))
    if msg_bit_indices is None:
        msg_bit_indices = list(range(512))

    base_hash = sha256_compress(msg, R)
    n_hb = len(hash_bit_indices)
    n_mb = len(msg_bit_indices)

    # Precompute single-flip hashes
    single_flip = {}
    for mi in msg_bit_indices:
        wi, bi = mi // 32, mi % 32
        m2 = list(msg)
        m2[wi] ^= (1 << bi)
        single_flip[mi] = sha256_compress(m2, R)

    # For each pair (mi, mj), compute double-flip hash
    hessian_nonzero = 0
    hessian_total = 0
    hessian_rows = {}  # (mi, mj) -> which hash bits have H=1

    # Sample pairs to keep tractable
    # Full would be C(512,2) ≈ 131K pairs × 256 hash bits
    pairs = []
    for i in range(len(msg_bit_indices)):
        for j in range(i+1, len(msg_bit_indices)):
            pairs.append((msg_bit_indices[i], msg_bit_indices[j]))

    print(f"  Computing Hessian for {len(pairs)} pairs...")

    for mi, mj in pairs:
        wi_i, bi_i = mi // 32, mi % 32
        wi_j, bi_j = mj // 32, mj % 32

        # Double flip
        m2 = list(msg)
        m2[wi_i] ^= (1 << bi_i)
        m2[wi_j] ^= (1 << bi_j)
        double_hash = sha256_compress(m2, R)

        # H = f(m) ^ f(m^ei) ^ f(m^ej) ^ f(m^ei^ej)
        for hb in hash_bit_indices:
            hw, hbi = hb // 32, hb % 32
            h_val = (get_bit(base_hash[hw], hbi) ^
                     get_bit(single_flip[mi][hw], hbi) ^
                     get_bit(single_flip[mj][hw], hbi) ^
                     get_bit(double_hash[hw], hbi))
            if h_val:
                hessian_nonzero += 1
            hessian_total += 1

    frac = hessian_nonzero / hessian_total if hessian_total else 0
    print(f"  Hessian nonzero: {hessian_nonzero}/{hessian_total} = {frac:.4f}")
    print(f"  (0.0 = affine, 0.5 = random)")

    return {
        'nonzero': hessian_nonzero,
        'total': hessian_total,
        'fraction': frac,
    }


def xl_degree2_attack(R, msg=None, seed=42, n_msg_bits=32, verbose=True):
    """
    XL degree-2 attack on R-round SHA-256.

    Reduced version: use first n_msg_bits message bits (instead of all 512)
    to keep computation tractable.

    1. Compute Jacobian J (256 × n) and Hessian H (256 × n × n)
    2. Build XL system: each row is a hash-bit equation in degree-2 variables
       y_i = δm_i (original), z_{ij} = δm_i * δm_j (products)
    3. Total XL variables: n + C(n,2)
    4. Find rank of XL system — if rank < variables, there's redundancy
    """
    if msg is None:
        rng = random.Random(seed)
        msg = [rng.randint(0, MASK32) for _ in range(16)]

    msg_indices = list(range(n_msg_bits))  # first n_msg_bits
    n = len(msg_indices)
    n_products = n * (n - 1) // 2
    n_xl_vars = n + n_products

    if verbose:
        print(f"\nXL Degree-2 Attack on R={R}")
        print(f"  Message bits used: {n}")
        print(f"  Product variables: {n_products}")
        print(f"  Total XL variables: {n_xl_vars}")

    base_hash = sha256_compress(msg, R)

    # Compute single-flip effects
    single_effects = {}  # mi -> 256-bit hash diff
    for mi in msg_indices:
        wi, bi = mi // 32, mi % 32
        m2 = list(msg)
        m2[wi] ^= (1 << bi)
        h2 = sha256_compress(m2, R)
        diff = 0
        for w in range(8):
            for b in range(32):
                if get_bit(base_hash[w], b) != get_bit(h2[w], b):
                    diff |= (1 << (w * 32 + b))
        single_effects[mi] = diff

    # Compute Hessian entries for pairs
    pair_index = {}
    idx = n  # products start after linear variables
    for i in range(n):
        for j in range(i+1, n):
            pair_index[(msg_indices[i], msg_indices[j])] = idx
            idx += 1

    # Build XL rows: for each hash bit, one equation
    xl_rows = []
    for hb in range(256):
        row = 0  # bit-vector over n_xl_vars

        # Linear part: J[hb][mi]
        for i, mi in enumerate(msg_indices):
            if (single_effects[mi] >> hb) & 1:
                row |= (1 << i)

        # Quadratic part: H[hb][mi][mj]
        for i in range(n):
            mi = msg_indices[i]
            wi_i, bi_i = mi // 32, mi % 32
            for j in range(i+1, n):
                mj = msg_indices[j]
                wi_j, bi_j = mj // 32, mj % 32

                # Double flip
                m2 = list(msg)
                m2[wi_i] ^= (1 << bi_i)
                m2[wi_j] ^= (1 << bi_j)
                h2 = sha256_compress(m2, R)

                hw, hbi = hb // 32, hb % 32
                h_val = (get_bit(base_hash[hw], hbi) ^
                         get_bit(single_effects[mi] >> hb, 0) ^  # wrong
                         get_bit(single_effects[mj] >> hb, 0))  # wrong
                # Fix: compute properly
                h_base = get_bit(base_hash[hw], hbi)
                h_ei = get_bit(sha256_compress(
                    [msg[w] ^ ((1 << bi_i) if w == wi_i else 0) for w in range(16)], R)[hw], hbi)
                h_ej = get_bit(sha256_compress(
                    [msg[w] ^ ((1 << bi_j) if w == wi_j else 0) for w in range(16)], R)[hw], hbi)
                h_eij = get_bit(h2[hw], hbi)

                hessian_bit = h_base ^ h_ei ^ h_ej ^ h_eij

                if hessian_bit:
                    col = pair_index[(mi, mj)]
                    row |= (1 << col)

        xl_rows.append(row)

    # Analyze XL system rank
    echelon, pivots = gf2_gaussian_eliminate(list(xl_rows), n_xl_vars)
    rank = len(pivots)
    kernel_dim = n_xl_vars - rank

    if verbose:
        print(f"  XL system: {len(xl_rows)} equations × {n_xl_vars} variables")
        print(f"  XL rank: {rank}")
        print(f"  XL kernel: {kernel_dim}")
        print(f"  Linear-only rank would be: {min(256, n)}")

        # Check how many pivot columns are in product range
        linear_pivots = sum(1 for p in pivots if p < n)
        product_pivots = sum(1 for p in pivots if p >= n)
        print(f"  Pivots in linear vars: {linear_pivots}")
        print(f"  Pivots in product vars: {product_pivots}")

    return {
        'R': R,
        'n_msg_bits': n,
        'n_xl_vars': n_xl_vars,
        'rank': rank,
        'kernel_dim': kernel_dim,
        'linear_pivots': linear_pivots,
        'product_pivots': product_pivots,
    }


def xl_redundancy_search(R=16, verbose=True):
    """
    Search for XL redundancy at R=16 across different message bit subsets.
    """
    rng = random.Random(42)
    msg = [rng.randint(0, MASK32) for _ in range(16)]

    if verbose:
        print("=" * 60)
        print(f"XL Redundancy Search at R={R}")
        print("=" * 60)

    # Start small, increase
    for n_bits in [8, 16, 24, 32]:
        res = xl_degree2_attack(R, msg=msg, n_msg_bits=n_bits, verbose=verbose)

        # Check: is there more info in degree-2 than degree-1?
        # If XL rank > linear rank, degree-2 carries new constraints
        # If XL kernel < expected, there's structure to exploit
        expected_kernel = n_bits + n_bits*(n_bits-1)//2 - 256
        if expected_kernel < 0:
            expected_kernel = 0

        if verbose:
            print(f"  Expected kernel (if independent): {max(0, res['n_xl_vars'] - 256)}")
            actual_excess = res['kernel_dim'] - max(0, res['n_xl_vars'] - 256)
            print(f"  Excess kernel (redundancy signal): {actual_excess}")
            print()


if __name__ == '__main__':
    # First check Hessian profile at different R
    rng = random.Random(42)
    msg = [rng.randint(0, MASK32) for _ in range(16)]

    print("=" * 60)
    print("Hessian Profile")
    print("=" * 60)
    for R in [4, 8, 12, 16]:
        print(f"\nR={R}:")
        compute_hessian_sample(R, msg,
                               hash_bit_indices=list(range(32)),  # first word
                               msg_bit_indices=list(range(32)))   # first word

    print("\n")
    xl_redundancy_search(R=16)
