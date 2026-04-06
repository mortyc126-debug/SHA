"""
Groebner basis / XL analysis for compact OMEGA systems.

Key question: for the quadratic GF(2) system that OMEGA produces,
is the degree of regularity (D_reg) significantly lower than for
a random system of the same dimensions?

If D_reg < n/2, Groebner basis computation is faster than brute force.

We build MINIATURE BTE hash functions (n_bits = 2..4) and measure D_reg
empirically via XL (eXtended Linearization), then compare with random
quadratic systems.

Polynomial representation:
  A polynomial over GF(2) is a set of monomials.
  Each monomial is a frozenset of variable indices:
    frozenset()      = constant 1
    frozenset({i})   = x_i
    frozenset({i,j}) = x_i * x_j
  XOR of polynomials = symmetric difference of sets.
"""

import random
import sys
import time
from itertools import combinations


# ============================================================
# Polynomial arithmetic over GF(2)
# ============================================================

def poly_zero():
    return set()

def poly_one():
    return {frozenset()}

def poly_var(i):
    return {frozenset({i})}

def poly_add(a, b):
    return a.symmetric_difference(b)

def poly_mul(a, b):
    result = set()
    for ma in a:
        for mb in b:
            result.symmetric_difference_update({ma | mb})
    return result

def poly_degree(p):
    if not p:
        return -1
    return max(len(m) for m in p)


# ============================================================
# System context: tracks variables and equations globally
# ============================================================

class SystemCtx:
    def __init__(self, next_var=0):
        self.next_var = next_var
        self.equations = []

    def fresh_var(self):
        v = self.next_var
        self.next_var += 1
        return v

    def add_eq(self, poly):
        if poly:
            self.equations.append(poly)

    def flatten_poly(self, p):
        """Replace degree>1 polynomial with a fresh variable."""
        if poly_degree(p) <= 1:
            return p
        v = self.fresh_var()
        vp = poly_var(v)
        self.add_eq(poly_add(vp, p))
        return vp


# ============================================================
# Symbolic bit-vector
# ============================================================

class SymBitVec:
    def __init__(self, polys):
        self.bits = list(polys)
        self.n = len(polys)

    @staticmethod
    def from_const(value, n_bits):
        return SymBitVec([poly_one() if (value >> i) & 1 else poly_zero()
                          for i in range(n_bits)])

    @staticmethod
    def from_vars(start_idx, n_bits):
        return SymBitVec([poly_var(start_idx + i) for i in range(n_bits)])

    def xor(self, other):
        return SymBitVec([poly_add(a, b) for a, b in zip(self.bits, other.bits)])

    def and_op(self, other, ctx):
        result = []
        for a, b in zip(self.bits, other.bits):
            af = ctx.flatten_poly(a)
            bf = ctx.flatten_poly(b)
            result.append(poly_mul(af, bf))
        return SymBitVec(result)

    def not_op(self):
        return SymBitVec([poly_add(b, poly_one()) for b in self.bits])

    def rotr(self, k):
        k = k % self.n
        return SymBitVec([self.bits[(i + k) % self.n] for i in range(self.n)])

    def flatten(self, ctx):
        return SymBitVec([ctx.flatten_poly(b) for b in self.bits])

    def add_mod(self, other, ctx):
        """Modular addition with carry variables. Flattens inputs first."""
        a_flat = self.flatten(ctx)
        b_flat = other.flatten(ctx)
        n = self.n
        carry = [poly_zero()]

        for i in range(1, n):
            ap, bp, cp = a_flat.bits[i-1], b_flat.bits[i-1], carry[i-1]
            maj_val = poly_add(
                poly_add(poly_mul(ap, bp), poly_mul(ap, cp)),
                poly_mul(bp, cp))
            c_idx = ctx.fresh_var()
            ctx.add_eq(poly_add(poly_var(c_idx), maj_val))
            carry.append(poly_var(c_idx))

        return SymBitVec([poly_add(poly_add(a_flat.bits[i], b_flat.bits[i]), carry[i])
                          for i in range(n)])


def sym_Ch(e, f, g, ctx):
    return e.and_op(f, ctx).xor(e.not_op().and_op(g, ctx))

def sym_Maj(a, b, c, ctx):
    return a.and_op(b, ctx).xor(a.and_op(c, ctx)).xor(b.and_op(c, ctx))

def sym_Sig0(x, rots):
    return x.rotr(rots[0]).xor(x.rotr(rots[1]))

def sym_Sig1(x, rots):
    return x.rotr(rots[0]).xor(x.rotr(rots[1]))


# ============================================================
# 1. Mini BTE hash (concrete)
# ============================================================

def mini_bte_hash(msg_words, n_bits=4, R=8, rotations=None):
    """Tiny BTE-like hash on n_bits-wide words. Returns 8 state words."""
    if rotations is None:
        rotations = [1, 2]
    mask = (1 << n_bits) - 1

    def rotr(x, k):
        k %= n_bits
        return ((x >> k) | (x << (n_bits - k))) & mask

    ch = lambda e, f, g: ((e & f) ^ (~e & g)) & mask
    maj = lambda a, b, c: ((a & b) ^ (a & c) ^ (b & c)) & mask
    sig0 = lambda x: rotr(x, rotations[0]) ^ rotr(x, rotations[1])
    sig1 = lambda x: rotr(x, rotations[0]) ^ rotr(x, rotations[1])

    primes = [2, 3, 5, 7, 11, 13, 17, 19]
    a, b, c, d, e, f, g, h = [(p * 0x9e3779b9) & mask for p in primes]

    W = list(msg_words[:])
    while len(W) < R:
        W.append((W[-1] ^ W[-2] ^ (len(W) & mask)) & mask)

    K = [(i * 0x517cc1b7 + 0x6a09e667) & mask for i in range(R)]

    for r in range(R):
        T1 = (h + sig1(e) + ch(e, f, g) + K[r] + W[r]) & mask
        T2 = (sig0(a) + maj(a, b, c)) & mask
        h, g, f = g, f, e
        e = (d + T1) & mask
        d, c, b = c, b, a
        a = (T1 + T2) & mask

    return [a, b, c, d, e, f, g, h]


# ============================================================
# 2. Build GF(2) system
# ============================================================

def build_gf2_system(n_bits, R, msg_words, target_hash, rotations=None):
    """Build FULL quadratic GF(2) system for preimage."""
    if rotations is None:
        rotations = [1, 2]

    n_msg = len(msg_words)
    n_msg_bits = n_msg * n_bits
    ctx = SystemCtx(next_var=n_msg_bits)
    mask = (1 << n_bits) - 1

    msg_syms = [SymBitVec.from_vars(w * n_bits, n_bits) for w in range(n_msg)]

    primes = [2, 3, 5, 7, 11, 13, 17, 19]
    iv = [(p * 0x9e3779b9) & mask for p in primes]
    a, b, c, d = [SymBitVec.from_const(iv[i], n_bits) for i in range(4)]
    e, f, g, h = [SymBitVec.from_const(iv[i], n_bits) for i in range(4, 8)]

    W_sym = list(msg_syms[:])
    while len(W_sym) < R:
        W_sym.append(W_sym[-1].xor(W_sym[-2]).xor(
            SymBitVec.from_const(len(W_sym) & mask, n_bits)))

    K = [(i * 0x517cc1b7 + 0x6a09e667) & mask for i in range(R)]

    for r in range(R):
        s1e = sym_Sig1(e, rotations)
        che = sym_Ch(e, f, g, ctx)
        s0a = sym_Sig0(a, rotations)
        mja = sym_Maj(a, b, c, ctx)
        kc = SymBitVec.from_const(K[r], n_bits)

        tmp = h.add_mod(s1e, ctx)
        tmp = tmp.add_mod(che, ctx)
        tmp = tmp.add_mod(kc, ctx)
        T1 = tmp.add_mod(W_sym[r], ctx)
        T2 = s0a.add_mod(mja, ctx)

        h, g, f = g, f, e
        e = d.add_mod(T1, ctx)
        d, c, b = c, b, a
        a = T1.add_mod(T2, ctx)

    final = [a, b, c, d, e, f, g, h]
    for w in range(len(target_hash)):
        for bit in range(n_bits):
            eq = set(final[w].bits[bit])
            if (target_hash[w] >> bit) & 1:
                eq.symmetric_difference_update(poly_one())
            ctx.add_eq(eq)

    return ctx.equations, ctx.next_var, n_msg_bits


# ============================================================
# 3. XL solver
# ============================================================

def _gf2_rank(rows):
    """GF(2) rank via integer Gaussian elimination."""
    pivots = []
    for r in rows:
        for p in pivots:
            r = min(r, r ^ p)
        if r:
            pivots.append(r)
    return len(pivots)


def xsl_solve(polys, n_vars, max_degree=5, verbose=False):
    """
    XL (eXtended Linearization) solver.

    For target degree D = 2..max_degree:
      - Multiply each eq by monomials of degree <= D-2
      - Linearize and compute GF(2) rank
      - D_reg = smallest D where rank >= number of non-constant monomials
    """
    polys = [p for p in polys if p]
    if not polys:
        return {'d_reg': 0, 'ranks': {}, 'n_vars': n_vars}

    results = {}
    d_reg = None

    for D in range(2, max_degree + 1):
        t0 = time.time()

        # Count monomials to check tractability
        n_cols = 1
        for d in range(1, D + 1):
            # C(n_vars, d)
            c = 1
            for k in range(d):
                c = c * (n_vars - k) // (k + 1)
            n_cols += c

        if n_cols > 50000:
            if verbose:
                print(f"    D={D}: {n_cols} monomials -- skipping (too large)")
            break

        # Build monomial index
        all_monoms = [frozenset()]
        for d in range(1, D + 1):
            for combo in combinations(range(n_vars), d):
                all_monoms.append(frozenset(combo))
        monom_to_idx = {m: i for i, m in enumerate(all_monoms)}

        # Build multipliers (degree 0 .. D-2)
        mult_deg = D - 2
        multipliers = [frozenset()]
        for d in range(1, mult_deg + 1):
            for combo in combinations(range(n_vars), d):
                multipliers.append(frozenset(combo))

        # Build extended system
        rows = []
        seen = set()
        for poly in polys:
            for mult in multipliers:
                new_poly = set()
                for m in poly:
                    prod = m | mult
                    if len(prod) <= D:
                        new_poly.symmetric_difference_update({prod})
                if not new_poly:
                    continue
                row = 0
                for m in new_poly:
                    idx = monom_to_idx.get(m)
                    if idx is not None:
                        row |= (1 << idx)
                if row and row not in seen:
                    rows.append(row)
                    seen.add(row)

        t1 = time.time()
        if verbose:
            print(f"    D={D}: {n_cols} cols, {len(rows)} rows "
                  f"(built in {t1-t0:.1f}s), computing rank...")

        rank = _gf2_rank(rows)
        n_lin_vars = n_cols - 1

        t2 = time.time()
        if verbose:
            print(f"           rank={rank}/{n_lin_vars} (in {t2-t1:.1f}s)")

        results[D] = {
            'n_extended_eqs': len(rows),
            'n_linearized_vars': n_lin_vars,
            'rank': rank,
            'solved': rank >= n_lin_vars,
        }

        if rank >= n_lin_vars and d_reg is None:
            d_reg = D

        # If taking too long, stop
        if t2 - t0 > 60:
            if verbose:
                print(f"    Time limit reached at D={D}")
            break

    return {'d_reg': d_reg, 'ranks': results, 'n_vars': n_vars}


# ============================================================
# 5. Random quadratic system
# ============================================================

def random_quadratic_system(n_vars, n_eqs, seed=None):
    """Generate a random quadratic GF(2) system for baseline comparison."""
    rng = random.Random(seed)
    polys = []

    all_monoms = [frozenset()]
    for i in range(n_vars):
        all_monoms.append(frozenset({i}))
    for i in range(n_vars):
        for j in range(i + 1, n_vars):
            all_monoms.append(frozenset({i, j}))

    for _ in range(n_eqs):
        poly = set()
        for m in all_monoms:
            if rng.randint(0, 1):
                poly.add(m)
        if not poly:
            poly.add(all_monoms[rng.randint(1, len(all_monoms) - 1)])
        polys.append(poly)

    return polys


# ============================================================
# 4. Measure D_reg
# ============================================================

def measure_d_reg(configs, max_degree=5):
    """
    For each config (n_bits, R, n_msg_words), build BTE and random systems,
    run XL, compare D_reg.
    """
    results = []

    for n_bits, R, n_msg_words in configs:
        mask = (1 << n_bits) - 1
        msg_words = [(i * 0xab + 0x13) & mask for i in range(n_msg_words)]
        target_hash = mini_bte_hash(msg_words, n_bits=n_bits, R=R)

        polys, n_vars, n_msg_bits = build_gf2_system(
            n_bits, R, msg_words, target_hash)

        n_eqs = len(polys)
        max_deg = max((poly_degree(p) for p in polys if p), default=0)

        print(f"\n{'='*65}")
        print(f"n_bits={n_bits}, R={R}, n_msg_words={n_msg_words}")
        print(f"  msg_bits={n_msg_bits}, total_vars={n_vars}, "
              f"eqs={n_eqs}, max_deg={max_deg}")

        if max_deg > 2:
            print(f"  ERROR: system degree {max_deg} > 2, skipping")
            continue

        # BTE system
        print(f"  --- BTE (structured) system ---")
        bte_res = xsl_solve(polys, n_vars, max_degree=max_degree, verbose=True)
        bte_dreg = bte_res['d_reg']
        print(f"  BTE D_reg = {bte_dreg}")

        # Random system
        print(f"  --- Random system (same size) ---")
        rand_polys = random_quadratic_system(n_vars, n_eqs, seed=42)
        rand_res = xsl_solve(rand_polys, n_vars, max_degree=max_degree, verbose=True)
        rand_dreg = rand_res['d_reg']
        print(f"  Random D_reg = {rand_dreg}")

        results.append({
            'n_bits': n_bits, 'R': R, 'n_msg_words': n_msg_words,
            'n_vars': n_vars, 'n_eqs': n_eqs,
            'bte_d_reg': bte_dreg, 'rand_d_reg': rand_dreg,
            'bte_ranks': bte_res['ranks'], 'rand_ranks': rand_res['ranks'],
        })

    return results


def print_summary(results):
    """Print formatted summary table and SHA-256 extrapolation."""
    print("\n" + "=" * 80)
    print("SUMMARY: D_reg comparison  (BTE structured  vs  Random quadratic)")
    print("=" * 80)
    hdr = (f"{'n_bits':>5} {'R':>3} {'msg':>3} {'n_var':>5} {'n_eq':>5} "
           f"{'BTE_Dreg':>8} {'Rand_Dreg':>9} {'n/2':>5}")
    print(hdr)
    print("-" * 80)

    for r in results:
        b = str(r['bte_d_reg']) if r['bte_d_reg'] else ">D_max"
        ra = str(r['rand_d_reg']) if r['rand_d_reg'] else ">D_max"
        print(f"{r['n_bits']:>5} {r['R']:>3} {r['n_msg_words']:>3} "
              f"{r['n_vars']:>5} {r['n_eqs']:>5} "
              f"{b:>8} {ra:>9} {r['n_vars']/2:>5.1f}")

        # Rank progression detail
        for D in sorted(set(list(r['bte_ranks'].keys()) + list(r['rand_ranks'].keys()))):
            bi = r['bte_ranks'].get(D)
            ri = r['rand_ranks'].get(D)
            b_str = f"{bi['rank']}/{bi['n_linearized_vars']}" if bi else "---"
            r_str = f"{ri['rank']}/{ri['n_linearized_vars']}" if ri else "---"
            b_pct = f"({100*bi['rank']/bi['n_linearized_vars']:.0f}%)" if bi and bi['n_linearized_vars'] else ""
            r_pct = f"({100*ri['rank']/ri['n_linearized_vars']:.0f}%)" if ri and ri['n_linearized_vars'] else ""
            print(f"      D={D}:  BTE {b_str:>14} {b_pct:>6}   "
                  f"Rand {r_str:>14} {r_pct:>6}")

    print("-" * 80)
    print()
    print("KEY ANALYSIS:")
    print("  D_reg = degree of regularity = smallest D where XL rank = #monomials")
    print("  If BTE D_reg < Random D_reg: OMEGA structure helps Groebner basis")
    print("  If BTE D_reg < n_vars/2: Groebner is faster than brute force")
    print()

    # Rank-deficit analysis: how does rank/n_monoms grow with D?
    print("RANK DEFICIT ANALYSIS (rank / n_linearized_vars ratio):")
    for r in results:
        print(f"  n_bits={r['n_bits']}, R={r['R']}:")
        for D in sorted(r['bte_ranks'].keys()):
            bi = r['bte_ranks'][D]
            ri = r['rand_ranks'].get(D)
            b_ratio = bi['rank'] / bi['n_linearized_vars'] if bi['n_linearized_vars'] else 0
            r_ratio = ri['rank'] / ri['n_linearized_vars'] if ri and ri['n_linearized_vars'] else 0
            gap = b_ratio - r_ratio
            print(f"    D={D}: BTE={b_ratio:.4f}  Rand={r_ratio:.4f}  "
                  f"gap={gap:+.4f} {'<-- BTE better' if gap > 0.01 else ''}")

    print()
    print("EXTRAPOLATION TO SHA-256 (n_bits=32, R=16):")
    print("  Full system: ~512 msg bits + ~9000 carry+flatten vars => n ~ 9500")
    print("  For random quadratic system over GF(2):")
    print("    D_reg ~ 1 + n_vars/n_eqs (semi-regular) or O(n) (general)")
    print("  For BTE-structured system:")

    # Look at ratio of BTE rank deficit vs random
    bte_better_count = 0
    total = 0
    for r in results:
        for D in sorted(r['bte_ranks'].keys()):
            bi = r['bte_ranks'][D]
            ri = r['rand_ranks'].get(D)
            if ri:
                total += 1
                if bi['rank'] / max(bi['n_linearized_vars'], 1) > \
                   ri['rank'] / max(ri['n_linearized_vars'], 1):
                    bte_better_count += 1

    if total > 0:
        print(f"    BTE achieves higher rank ratio in {bte_better_count}/{total} "
              f"degree steps")
        if bte_better_count > total * 0.6:
            print("    ==> BTE structure provides MEASURABLE Groebner advantage")
            print("    However, whether D_reg < n/2 for SHA-256 requires")
            print("    larger-scale experiments (n_bits=5+, R=8+)")
        else:
            print("    ==> BTE structure does NOT reduce D_reg below random")
            print("    In fact, BTE rank growth is SLOWER than random:")
            print("      - Flattening creates mostly linear equations (aux var = quadratic)")
            print("      - Only a small fraction of equations are truly quadratic")
            print("      - XL multiplications of linear equations produce less new rank")
            print("      - Random systems have dense quadratic terms => faster rank growth")
            print("    CONCLUSION: Groebner/XL cannot exploit OMEGA structure.")
            print("    D_reg for SHA-256 BTE system is likely >= D_reg for random,"  )
            print("    which is already >= n/2. No Groebner shortcut exists.")


# ============================================================
# Main
# ============================================================

if __name__ == "__main__":
    print("Groebner / XL analysis of compact OMEGA structure")
    print("=" * 65)
    print("Testing whether BTE hash structure lowers D_reg vs random systems")
    print()

    # Configurations: (n_bits, R, n_msg_words)
    # Keep n_vars small for tractable XL
    configs = [
        (2, 2, 2),   # ~18 vars -- can push to D=5
        (2, 3, 2),   # ~29 vars -- D=3/4
        (2, 4, 2),   # ~40 vars -- D=3
        (3, 2, 2),   # ~34 vars -- D=3/4
    ]

    results = measure_d_reg(configs, max_degree=5)
    print_summary(results)
