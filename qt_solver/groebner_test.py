"""
Groebner basis / XL analysis for compact OMEGA systems.

Key question: for the quadratic GF(2) system that OMEGA produces,
is the degree of regularity (D_reg) significantly lower than for
a random system of the same dimensions?

If D_reg < n/2, Groebner basis computation is faster than brute force.

We build MINIATURE BTE hash functions (n_bits = 3..5) and measure D_reg
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
from itertools import combinations

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False


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
            prod = ma | mb
            result.symmetric_difference_update({prod})
    return result

def poly_degree(p):
    if not p:
        return -1
    return max(len(m) for m in p)


# ============================================================
# System context: tracks variables and equations globally
# ============================================================

class SystemCtx:
    """Global context for building a quadratic GF(2) system."""
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
        """If p has degree > 1, introduce a fresh variable v and add
        equation v + p = 0. Return the variable polynomial."""
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
        bits = []
        for i in range(n_bits):
            bits.append(poly_one() if (value >> i) & 1 else poly_zero())
        return SymBitVec(bits)

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
        one = poly_one()
        return SymBitVec([poly_add(b, one) for b in self.bits])

    def rotr(self, k):
        n = self.n
        k = k % n
        return SymBitVec([self.bits[(i + k) % n] for i in range(n)])

    def flatten(self, ctx):
        return SymBitVec([ctx.flatten_poly(b) for b in self.bits])

    def add_mod(self, other, ctx):
        """Modular addition with carry variables. Inputs flattened first."""
        a_flat = self.flatten(ctx)
        b_flat = other.flatten(ctx)
        n = self.n
        carry_bits = [poly_zero()]

        for i in range(1, n):
            ap = a_flat.bits[i - 1]
            bp = b_flat.bits[i - 1]
            cp = carry_bits[i - 1]
            # MAJ(a,b,c) = ab + ac + bc; all inputs degree <= 1
            maj_val = poly_add(
                poly_add(poly_mul(ap, bp), poly_mul(ap, cp)),
                poly_mul(bp, cp)
            )
            c_idx = ctx.fresh_var()
            ctx.add_eq(poly_add(poly_var(c_idx), maj_val))
            carry_bits.append(poly_var(c_idx))

        sum_bits = [
            poly_add(poly_add(a_flat.bits[i], b_flat.bits[i]), carry_bits[i])
            for i in range(n)
        ]
        return SymBitVec(sum_bits)


def sym_Ch(e, f, g, ctx):
    return e.and_op(f, ctx).xor(e.not_op().and_op(g, ctx))

def sym_Maj(a, b, c, ctx):
    return a.and_op(b, ctx).xor(a.and_op(c, ctx)).xor(b.and_op(c, ctx))

def sym_Sig0(x, rotations):
    return x.rotr(rotations[0]).xor(x.rotr(rotations[1]))

def sym_Sig1(x, rotations):
    return x.rotr(rotations[0]).xor(x.rotr(rotations[1]))


# ============================================================
# 1. Mini BTE hash (concrete)
# ============================================================

def mini_bte_hash(msg_words, n_bits=4, R=8, rotations=None):
    """Tiny BTE-like hash on n_bits-wide words. Returns 8 state words."""
    if rotations is None:
        rotations = [1, 2]
    mask = (1 << n_bits) - 1

    def rotr(x, k):
        k = k % n_bits
        return ((x >> k) | (x << (n_bits - k))) & mask

    def ch(e, f, g):
        return ((e & f) ^ (~e & g)) & mask

    def maj(a, b, c):
        return ((a & b) ^ (a & c) ^ (b & c)) & mask

    def sig0(x):
        return rotr(x, rotations[0]) ^ rotr(x, rotations[1])

    def sig1(x):
        return rotr(x, rotations[0]) ^ rotr(x, rotations[1])

    primes = [2, 3, 5, 7, 11, 13, 17, 19]
    state = [(p * 0x9e3779b9) & mask for p in primes]
    a, b, c, d, e, f, g, h = state

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
    """Build FULL quadratic GF(2) system for preimage. Returns (polys, n_vars, n_msg_bits)."""
    if rotations is None:
        rotations = [1, 2]

    n_msg = len(msg_words)
    n_msg_bits = n_msg * n_bits
    ctx = SystemCtx(next_var=n_msg_bits)

    msg_syms = [SymBitVec.from_vars(w * n_bits, n_bits) for w in range(n_msg)]

    mask = (1 << n_bits) - 1
    primes = [2, 3, 5, 7, 11, 13, 17, 19]
    init_state = [(p * 0x9e3779b9) & mask for p in primes]

    a = SymBitVec.from_const(init_state[0], n_bits)
    b = SymBitVec.from_const(init_state[1], n_bits)
    c = SymBitVec.from_const(init_state[2], n_bits)
    d = SymBitVec.from_const(init_state[3], n_bits)
    e = SymBitVec.from_const(init_state[4], n_bits)
    f = SymBitVec.from_const(init_state[5], n_bits)
    g = SymBitVec.from_const(init_state[6], n_bits)
    h = SymBitVec.from_const(init_state[7], n_bits)

    W_sym = list(msg_syms[:])
    while len(W_sym) < R:
        W_sym.append(W_sym[-1].xor(W_sym[-2]).xor(
            SymBitVec.from_const(len(W_sym) & mask, n_bits)))

    K = [(i * 0x517cc1b7 + 0x6a09e667) & mask for i in range(R)]

    for r in range(R):
        sig1_e = sym_Sig1(e, rotations)
        ch_efg = sym_Ch(e, f, g, ctx)
        sig0_a = sym_Sig0(a, rotations)
        maj_abc = sym_Maj(a, b, c, ctx)
        k_const = SymBitVec.from_const(K[r], n_bits)

        tmp = h.add_mod(sig1_e, ctx)
        tmp = tmp.add_mod(ch_efg, ctx)
        tmp = tmp.add_mod(k_const, ctx)
        T1 = tmp.add_mod(W_sym[r], ctx)

        T2 = sig0_a.add_mod(maj_abc, ctx)

        h, g, f = g, f, e
        e = d.add_mod(T1, ctx)
        d, c, b = c, b, a
        a = T1.add_mod(T2, ctx)

    final_state = [a, b, c, d, e, f, g, h]
    for w in range(len(target_hash)):
        for bit in range(n_bits):
            target_bit = (target_hash[w] >> bit) & 1
            eq = set(final_state[w].bits[bit])
            if target_bit:
                eq.symmetric_difference_update(poly_one())
            ctx.add_eq(eq)

    return ctx.equations, ctx.next_var, n_msg_bits


# ============================================================
# 3. XL solver (numpy-accelerated)
# ============================================================

def _build_monom_index(n_vars, max_deg):
    """Build monomial list and lookup, up to degree max_deg."""
    monoms = [frozenset()]
    for d in range(1, max_deg + 1):
        for combo in combinations(range(n_vars), d):
            monoms.append(frozenset(combo))
    idx = {m: i for i, m in enumerate(monoms)}
    return monoms, idx


def _gf2_rank_numpy(rows_list, n_cols):
    """GF(2) rank using numpy uint8 matrix and row reduction."""
    if not rows_list:
        return 0
    n_rows = len(rows_list)
    # Pack into numpy array
    mat = np.zeros((n_rows, n_cols), dtype=np.uint8)
    for i, row in enumerate(rows_list):
        for j in range(n_cols):
            if (row >> j) & 1:
                mat[i, j] = 1

    # Gaussian elimination over GF(2)
    pivot_row = 0
    for col in range(n_cols):
        # Find pivot
        found = -1
        for r in range(pivot_row, n_rows):
            if mat[r, col]:
                found = r
                break
        if found < 0:
            continue
        # Swap
        if found != pivot_row:
            mat[[pivot_row, found]] = mat[[found, pivot_row]]
        # Eliminate
        mask = mat[:, col].astype(bool)
        mask[pivot_row] = False
        mat[mask] ^= mat[pivot_row]
        pivot_row += 1

    return pivot_row


def _gf2_rank_int(rows):
    """GF(2) rank via integer Gaussian elimination."""
    pivots = []
    for r in rows:
        for p in pivots:
            r = min(r, r ^ p)
        if r:
            pivots.append(r)
    return len(pivots)


def xsl_solve(polys, n_vars, max_degree=4, verbose=False):
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
        mult_deg = D - 2

        # Build multiplier monomials
        multipliers = [frozenset()]
        for d in range(1, mult_deg + 1):
            for combo in combinations(range(n_vars), d):
                multipliers.append(frozenset(combo))

        # Build column index (all monomials up to degree D)
        all_monoms, monom_to_idx = _build_monom_index(n_vars, D)
        n_cols = len(all_monoms)

        if verbose:
            print(f"    D={D}: {len(multipliers)} multipliers, "
                  f"{n_cols} monomials, generating rows...")

        # Build extended rows
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
                    if m in monom_to_idx:
                        row |= (1 << monom_to_idx[m])
                if row and row not in seen:
                    rows.append(row)
                    seen.add(row)

        # Compute rank
        if HAS_NUMPY and n_cols > 64:
            rank = _gf2_rank_numpy(rows, n_cols)
        else:
            rank = _gf2_rank_int(rows)

        n_lin_vars = n_cols - 1  # exclude constant monomial

        results[D] = {
            'n_extended_eqs': len(rows),
            'n_linearized_vars': n_lin_vars,
            'rank': rank,
            'solved': rank >= n_lin_vars,
        }

        if rank >= n_lin_vars and d_reg is None:
            d_reg = D

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
        if poly:
            polys.append(poly)
        else:
            polys.append({all_monoms[rng.randint(1, len(all_monoms) - 1)]})

    return polys


# ============================================================
# 4. Measure D_reg
# ============================================================

def measure_d_reg(n_bits_list=None, R_list=None, max_degree=4, n_msg_words=2):
    """
    For each (n_bits, R), build mini BTE system and random system,
    run XL, compare D_reg.
    """
    if n_bits_list is None:
        n_bits_list = [3, 4]
    if R_list is None:
        R_list = [2, 3, 4]

    results = []

    for n_bits in n_bits_list:
        for R in R_list:
            mask = (1 << n_bits) - 1
            msg_words = [(i * 0xab + 0x13) & mask for i in range(n_msg_words)]

            target_hash = mini_bte_hash(msg_words, n_bits=n_bits, R=R)

            print(f"\n{'='*60}")
            print(f"n_bits={n_bits}, R={R}, n_msg_words={n_msg_words}")

            polys, n_vars, n_msg_bits = build_gf2_system(
                n_bits, R, msg_words, target_hash)

            n_eqs = len(polys)
            max_deg_actual = max((poly_degree(p) for p in polys if p), default=0)

            print(f"  Message bits: {n_msg_bits}")
            print(f"  Total vars (msg + carry + flatten): {n_vars}")
            print(f"  Equations: {n_eqs}")
            print(f"  Max equation degree: {max_deg_actual}")

            if max_deg_actual > 2:
                print(f"  WARNING: system has degree {max_deg_actual}, expected <= 2")
                continue

            # Estimate tractability: degree-3 has C(n,3) + C(n,2) + n + 1 columns
            n3_cols = 1 + n_vars + n_vars*(n_vars-1)//2
            if max_degree >= 3:
                n3_cols += n_vars*(n_vars-1)*(n_vars-2)//6
            print(f"  Linearized vars at D=3: ~{n3_cols}")
            if n3_cols > 500000:
                print(f"  SKIPPING: too many monomials for tractable XL")
                continue

            # XL on BTE system
            print(f"  Running XL on BTE system...")
            bte_result = xsl_solve(polys, n_vars, max_degree=max_degree, verbose=True)

            bte_dreg = bte_result['d_reg']
            print(f"  BTE D_reg: {bte_dreg}")
            for D, info in sorted(bte_result['ranks'].items()):
                print(f"    D={D}: rank={info['rank']}/{info['n_linearized_vars']} "
                      f"({info['n_extended_eqs']} ext eqs)"
                      f"{'  ** SOLVED **' if info['solved'] else ''}")

            # Random system of same dimensions
            print(f"  Running XL on random system ({n_vars} vars, {n_eqs} eqs)...")
            rand_polys = random_quadratic_system(n_vars, n_eqs, seed=42)
            rand_result = xsl_solve(rand_polys, n_vars, max_degree=max_degree, verbose=True)

            rand_dreg = rand_result['d_reg']
            print(f"  Random D_reg: {rand_dreg}")
            for D, info in sorted(rand_result['ranks'].items()):
                print(f"    D={D}: rank={info['rank']}/{info['n_linearized_vars']} "
                      f"({info['n_extended_eqs']} ext eqs)"
                      f"{'  ** SOLVED **' if info['solved'] else ''}")

            results.append({
                'n_bits': n_bits,
                'R': R,
                'n_vars': n_vars,
                'n_eqs': n_eqs,
                'bte_d_reg': bte_dreg,
                'rand_d_reg': rand_dreg,
                'n_over_2': n_vars / 2,
            })

    return results


def print_summary(results):
    """Print formatted summary table and SHA-256 extrapolation."""
    print("\n" + "=" * 80)
    print("SUMMARY: D_reg comparison (BTE structured vs Random)")
    print("=" * 80)
    print(f"{'n_bits':>6} {'R':>4} {'n_vars':>7} {'n_eqs':>6} "
          f"{'BTE D_reg':>10} {'Rand D_reg':>11} {'n/2':>6} {'Advantage':>10}")
    print("-" * 80)

    for r in results:
        bte_s = str(r['bte_d_reg']) if r['bte_d_reg'] else ">max_D"
        rand_s = str(r['rand_d_reg']) if r['rand_d_reg'] else ">max_D"

        if r['bte_d_reg'] and r['rand_d_reg']:
            adv_s = f"{r['rand_d_reg'] - r['bte_d_reg']:+d}"
        else:
            adv_s = "N/A"

        print(f"{r['n_bits']:>6} {r['R']:>4} {r['n_vars']:>7} {r['n_eqs']:>6} "
              f"{bte_s:>10} {rand_s:>11} {r['n_over_2']:>6.1f} {adv_s:>10}")

    print("-" * 80)
    print()
    print("Interpretation:")
    print("  BTE D_reg < Random D_reg  =>  structure helps Groebner basis")
    print("  BTE D_reg < n/2           =>  Groebner faster than brute force")
    print()

    print("Extrapolation to SHA-256 (n_bits=32, R=16):")
    print("  Full SHA-256: ~512 msg bits + ~9000+ carry/flatten vars")
    print("  n_vars ~ 9500, n/2 ~ 4750")

    bte_ratios = [r['bte_d_reg'] / r['n_vars']
                  for r in results if r['bte_d_reg']]
    rand_ratios = [r['rand_d_reg'] / r['n_vars']
                   for r in results if r['rand_d_reg']]

    if bte_ratios:
        avg_bte = sum(bte_ratios) / len(bte_ratios)
        print(f"  Avg BTE D_reg/n_vars: {avg_bte:.4f}  =>  SHA-256 est: ~{avg_bte*9500:.0f}")
    if rand_ratios:
        avg_rand = sum(rand_ratios) / len(rand_ratios)
        print(f"  Avg Rand D_reg/n_vars: {avg_rand:.4f}  =>  SHA-256 est: ~{avg_rand*9500:.0f}")

    if bte_ratios:
        est = avg_bte * 9500
        if est < 4750:
            print("\n  ** BTE structure MAY allow Groebner advantage over brute force **")
        else:
            print("\n  BTE structure does NOT reduce D_reg below n/2.")
            print("  Groebner basis unlikely to beat brute force for SHA-256.")


# ============================================================
# Main
# ============================================================

if __name__ == "__main__":
    print("Groebner / XL analysis of compact OMEGA structure")
    print("Comparing BTE-structured systems vs random quadratic systems")

    # Use small parameters: 2 message words keeps n_vars tractable
    # n_bits=3, 2 msg words => 6 msg bits, ~30-50 total vars depending on R
    # n_bits=4, 2 msg words => 8 msg bits, ~40-70 total vars
    results = measure_d_reg(
        n_bits_list=[3, 4],
        R_list=[2, 3, 4],
        max_degree=4,
        n_msg_words=2,
    )

    print_summary(results)
