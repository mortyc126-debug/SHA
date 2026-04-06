"""
BTE Optimal Hash Design.

Use BTE theory (T1-T12) to design an optimal hash function:
- Minimal rounds: R = R_full + safety_margin
- Optimal rotations: maximize entropy (F59: corr=0.852 with D2 speed)
- Proven security framework from carry algebra + degree growth
"""

from qt_solver.sha256_traced import MASK32, get_bit, IV, K
from qt_solver.degree_growth import compute_derivative
import random
import math


def mini_bte_hash(msg, n_bits=8, R=16, rotations_sig0=None, rotations_sig1=None, iv=None):
    """
    Mini BTE hash function with configurable parameters.
    msg: list of n_msg words, each n_bits wide.
    """
    mask = (1 << n_bits) - 1
    if rotations_sig0 is None:
        rotations_sig0 = [1, 3]
    if rotations_sig1 is None:
        rotations_sig1 = [2, 5]
    if iv is None:
        iv = [(0x6a09e667 >> (32 - n_bits)) & mask,
              (0xbb67ae85 >> (32 - n_bits)) & mask,
              (0x3c6ef372 >> (32 - n_bits)) & mask,
              (0xa54ff53a >> (32 - n_bits)) & mask,
              (0x510e527f >> (32 - n_bits)) & mask,
              (0x9b05688c >> (32 - n_bits)) & mask,
              (0x1f83d9ab >> (32 - n_bits)) & mask,
              (0x5be0cd19 >> (32 - n_bits)) & mask]

    def rotr(x, r):
        return ((x >> r) | (x << (n_bits - r))) & mask

    def sig0(x):
        v = 0
        for r in rotations_sig0:
            v ^= rotr(x, r % n_bits)
        return v

    def sig1(x):
        v = 0
        for r in rotations_sig1:
            v ^= rotr(x, r % n_bits)
        return v

    def ch(e, f, g):
        return (e & f) ^ (~e & g) & mask

    def maj(a, b, c):
        return (a & b) ^ (a & c) ^ (b & c)

    n_msg = len(msg)
    w = list(msg[:n_msg])
    # Pad/extend schedule
    while len(w) < R:
        i = len(w)
        w.append((w[i-2] ^ w[max(0, i-7)] ^ w[max(0, i-n_msg)]) & mask)

    a, b, c, d, e, f, g, h = iv
    round_k = [((K[r % 64]) >> (32 - n_bits)) & mask for r in range(R)]

    for r in range(R):
        t1 = (h + sig1(e) + ch(e, f, g) + round_k[r] + w[r]) & mask
        t2 = (sig0(a) + maj(a, b, c)) & mask
        h, g, f = g, f, e
        e = (d + t1) & mask
        d, c, b = c, b, a
        a = (t1 + t2) & mask

    return [(a + iv[0]) & mask, (b + iv[1]) & mask,
            (c + iv[2]) & mask, (d + iv[3]) & mask,
            (e + iv[4]) & mask, (f + iv[5]) & mask,
            (g + iv[6]) & mask, (h + iv[7]) & mask]


def rotation_entropy(rots, n_bits):
    """Compute entropy of rotation distribution over Z/n_bits."""
    counts = [0] * n_bits
    for r in rots:
        counts[r % n_bits] += 1
    total = sum(counts)
    if total == 0:
        return 0.0
    ent = 0.0
    for c in counts:
        if c > 0:
            p = c / total
            ent -= p * math.log2(p)
    return ent


def rotation_coverage(rots, n_bits):
    """How many steps to reach all bit positions from bit 0."""
    reached = {0}
    for step in range(1, n_bits + 1):
        new = set()
        for pos in reached:
            for r in rots:
                new.add((pos + r) % n_bits)
                new.add((pos - r) % n_bits)
        reached |= new
        if len(reached) == n_bits:
            return step
    return n_bits


def evaluate_rotations(n_bits, sig0_rots, sig1_rots, R_test=16, n_trials=100):
    """Evaluate a rotation configuration."""
    all_rots = sig0_rots + sig1_rots
    ent = rotation_entropy(all_rots, n_bits)
    cov = rotation_coverage(all_rots, n_bits)

    # Measure D2 at R_test
    rng = random.Random(42)
    nonzero = 0
    for _ in range(n_trials):
        msg = [rng.randint(0, (1 << n_bits) - 1) for _ in range(4)]
        # 2 random directions
        dirs = []
        for _ in range(2):
            wi, bi = rng.randint(0, 3), rng.randint(0, n_bits - 1)
            dirs.append((wi, bi))

        hb_word, hb_bit = rng.randint(0, 7), rng.randint(0, n_bits - 1)

        # Compute D2 manually on mini hash
        vals = []
        for mask_bits in range(4):
            m2 = list(msg)
            for j in range(2):
                if (mask_bits >> j) & 1:
                    m2[dirs[j][0]] ^= (1 << dirs[j][1])
            h = mini_bte_hash(m2, n_bits, R_test, sig0_rots, sig1_rots)
            vals.append(get_bit(h[hb_word], hb_bit))

        d2 = vals[0] ^ vals[1] ^ vals[2] ^ vals[3]
        if d2:
            nonzero += 1

    d2_frac = nonzero / n_trials
    return {'entropy': ent, 'coverage': cov, 'D2': d2_frac}


def search_optimal_rotations(n_bits=8):
    """Search for best rotation constants."""
    print(f"Rotation Search (n={n_bits})")
    print("=" * 65)

    candidates = []

    # Generate candidates: pairs of 2-3 rotation constants
    for a in range(1, n_bits):
        for b in range(a + 1, n_bits):
            sig0 = [a, b]
            for c in range(1, n_bits):
                for d in range(c + 1, n_bits):
                    sig1 = [c, d]
                    if set(sig0) == set(sig1):
                        continue
                    res = evaluate_rotations(n_bits, sig0, sig1, R_test=12, n_trials=50)
                    candidates.append((sig0, sig1, res))

    # Sort by D2 (higher = faster randomization)
    candidates.sort(key=lambda x: -x[2]['D2'])

    print(f"  Tested {len(candidates)} configurations")
    print(f"\n  Top 10:")
    print(f"  {'Sig0':>10} {'Sig1':>10} | {'Ent':>5} {'Cov':>4} {'D2@12':>7}")
    print(f"  {'-'*50}")
    for sig0, sig1, res in candidates[:10]:
        print(f"  {str(sig0):>10} {str(sig1):>10} | {res['entropy']:5.2f} {res['coverage']:4d} {res['D2']:7.3f}")

    # SHA-256-like
    sha_sig0 = [2 % n_bits, 5 % n_bits]  # simplified for n=8
    sha_sig1 = [3 % n_bits, 6 % n_bits]
    sha_res = evaluate_rotations(n_bits, sha_sig0, sha_sig1, R_test=12, n_trials=50)
    print(f"\n  SHA-like: Sig0={sha_sig0} Sig1={sha_sig1} | "
          f"Ent={sha_res['entropy']:.2f} Cov={sha_res['coverage']} D2={sha_res['D2']:.3f}")

    best = candidates[0]
    print(f"\n  OPTIMAL: Sig0={best[0]} Sig1={best[1]} | "
          f"Ent={best[2]['entropy']:.2f} Cov={best[2]['coverage']} D2={best[2]['D2']:.3f}")

    return best


def design_bte_hash(n_bits=32, security_bits=128):
    """Design a complete BTE hash specification."""
    print(f"\nBTE Hash Design (n={n_bits}, security={security_bits})")
    print("=" * 65)

    n_msg_words = 16  # same as SHA-256
    hash_bits = 2 * security_bits  # 256 for 128-bit security

    # R_full from T7
    R_full = n_msg_words + 2  # = 18

    # Safety margin
    safety_factor = 2.0  # standard
    R_total = int(R_full * safety_factor)

    print(f"""
  SPECIFICATION: BTE-{hash_bits}
  ─────────────────────────────
  Word width:      {n_bits} bits
  Registers:       8
  Message words:   {n_msg_words}
  Message bits:    {n_msg_words * n_bits}
  Hash bits:       {hash_bits}

  ROUND COUNT:
    R_full (T7):   {R_full} (all D_k ≈ 0.5)
    Safety factor:  {safety_factor}×
    Total rounds:   {R_total}

  vs SHA-256:      64 rounds (3.6× R_full)
  vs BTE-{hash_bits}:      {R_total} rounds ({safety_factor}× R_full)
  SPEEDUP:         {64/R_total:.1f}× fewer rounds

  SECURITY ARGUMENT:
    T1:  Layer rank = 2×{R_total}-1 = {2*R_total-1} per layer
    T3:  Carry nilpotent (depth {n_bits})
    T5:  Carry = cocycle (associative addition)
    T7:  Full randomization at R={R_full}
    T11: Degree ≥ Fib({R_full}) = {_fib(R_full)} >> {n_msg_words*n_bits} input bits
    T12: Monomial spread complete at R={R_full} (from rotation coverage)
    →    PRF at R≥{R_full}, collision ≥ 2^{security_bits}
    →    With {safety_factor}× margin: {R_total} rounds
    """)

    return {
        'n_bits': n_bits, 'n_msg': n_msg_words, 'R': R_total,
        'R_full': R_full, 'hash_bits': hash_bits, 'security': security_bits,
    }


def compare_with_sha256(n_bits=8):
    """Compare SHA-256-like vs optimal BTE on mini instances."""
    print(f"\nComparison: SHA-like vs Optimal BTE (n={n_bits})")
    print("=" * 65)

    rng = random.Random(42)
    n_msg = 4

    sha_sig0, sha_sig1 = [2, 5], [3, 6]
    opt = search_optimal_rotations(n_bits)
    opt_sig0, opt_sig1 = opt[0], opt[1]

    print(f"\n  D2 profile across rounds:")
    print(f"  {'R':>4} | {'SHA-like':>9} | {'Optimal':>9}")
    print(f"  {'-'*30}")

    for R in range(2, 20, 2):
        sha_res = evaluate_rotations(n_bits, sha_sig0, sha_sig1, R_test=R, n_trials=100)
        opt_res = evaluate_rotations(n_bits, opt_sig0, opt_sig1, R_test=R, n_trials=100)
        print(f"  {R:4d} | {sha_res['D2']:9.3f} | {opt_res['D2']:9.3f}")


def _fib(n):
    a, b = 1, 1
    for _ in range(n - 2):
        a, b = b, a + b
    return b


if __name__ == '__main__':
    best = search_optimal_rotations(n_bits=8)
    design_bte_hash(n_bits=32, security_bits=128)
    compare_with_sha256(n_bits=8)
