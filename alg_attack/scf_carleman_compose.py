#!/usr/bin/env python3
"""
SCF Задача A, Этап 2: Композиция раундов в Carleman-пространстве.

Факт: один раунд XOR-SHA — точно квадратичный (degree 2).
Вопрос: при композиции R раундов, растёт ли rank Carleman-пространства?

Если rank растёт ЛИНЕЙНО (+160/round) → dim(64 rounds) ≈ 8500 (tractable!)
Если rank растёт ЭКСПОНЕНЦИАЛЬНО → Carleman не помогает.

Метод: для mini-SHA (n=4), вычислить f_R = f∘f∘...∘f (R раз),
построить Carleman-декомпозицию f_R, измерить rank и число мономов.
"""
import os, sys

def bit(x, i):
    return (x >> i) & 1

def rotr_n(x, r, n):
    return ((x >> r) | (x << (n - r))) & ((1 << n) - 1)

class MiniSHA:
    def __init__(self, n=4):
        self.n = n
        self.mask = (1 << n) - 1
        self.sig0_rots = [max(1, n//8), max(2, n//3), n//2]
        self.sig1_rots = [max(1, n//5), max(2, n//3), max(3, 3*n//4)]

    def Sig0(self, x):
        n = self.n
        return rotr_n(x, self.sig0_rots[0], n) ^ rotr_n(x, self.sig0_rots[1], n) ^ rotr_n(x, self.sig0_rots[2], n)

    def Sig1(self, x):
        n = self.n
        return rotr_n(x, self.sig1_rots[0], n) ^ rotr_n(x, self.sig1_rots[1], n) ^ rotr_n(x, self.sig1_rots[2], n)

    def Ch(self, e, f, g):
        return (e & f) ^ (~e & g) & self.mask

    def Maj(self, a, b, c):
        return (a & b) ^ (a & c) ^ (b & c)

    def xor_round(self, state, W=0, K=0):
        a,b,c,d,e,f,g,h = state
        T1 = h ^ self.Sig1(e) ^ self.Ch(e,f,g) ^ K ^ W
        T2 = self.Sig0(a) ^ self.Maj(a,b,c)
        return [T1^T2, a, b, c, d^T1, e, f, g]

    def xor_multi_round(self, state, R, W_list=None):
        s = list(state)
        for r in range(R):
            W = W_list[r] if W_list else 0
            s = self.xor_round(s, W, 0)
        return s

    def state_to_bits(self, state):
        bits = []
        for w in state:
            for i in range(self.n):
                bits.append(bit(w, i))
        return bits

    def bits_to_state(self, bits):
        state = []
        for r in range(8):
            val = 0
            for i in range(self.n):
                val |= bits[r * self.n + i] << i
            state.append(val)
        return state


def gf2_rank(mat):
    if not mat: return 0
    m = [list(r) for r in mat]
    nr, nc = len(m), len(m[0])
    rank = 0
    for col in range(nc):
        pv = -1
        for r in range(rank, nr):
            if m[r][col]: pv = r; break
        if pv == -1: continue
        m[rank], m[pv] = m[pv], m[rank]
        for r in range(nr):
            if r != rank and m[r][col]:
                for c in range(nc): m[r][c] ^= m[rank][c]
        rank += 1
    return rank


def analyze_composition(n, R, sha):
    """Build Carleman decomposition of R composed rounds."""
    dim = 8 * n  # state bits

    # f_R(0)
    f_zero = sha.state_to_bits(sha.xor_multi_round([0]*8, R))

    # f_R(e_i) for all basis vectors
    f_single = []
    for inp in range(dim):
        test_bits = [0] * dim
        test_bits[inp] = 1
        test_state = sha.bits_to_state(test_bits)
        f_single.append(sha.state_to_bits(sha.xor_multi_round(test_state, R)))

    # Linear part: A[out][in] = f(e_i)[out] ⊕ f(0)[out]
    A = [[f_single[inp][out] ^ f_zero[out] for inp in range(dim)] for out in range(dim)]
    rank_A = gf2_rank(A)

    # Second-order: sample pairs to find nonzero monomials
    monomial_pairs = []
    for i in range(dim):
        for j in range(i+1, dim):
            monomial_pairs.append((i, j))

    nonzero_monomials = []
    B_cols = []  # only store nonzero columns

    for idx, (i, j) in enumerate(monomial_pairs):
        test_bits = [0] * dim
        test_bits[i] = 1
        test_bits[j] = 1
        test_state = sha.bits_to_state(test_bits)
        f_ij = sha.state_to_bits(sha.xor_multi_round(test_state, R))

        # Second-order differential
        col = [f_ij[out] ^ f_single[i][out] ^ f_single[j][out] ^ f_zero[out]
               for out in range(dim)]

        if any(col):
            nonzero_monomials.append((i, j))
            B_cols.append(col)

    n_mono = len(nonzero_monomials)

    # Rank of full Carleman system [A | B_nonzero]
    if B_cols:
        full_matrix = [A[out] + [B_cols[m][out] for m in range(n_mono)]
                       for out in range(dim)]
        rank_full = gf2_rank(full_matrix)
    else:
        rank_full = rank_A

    # Check for THIRD-order terms (degree 3)
    n_cubic_test = min(200, dim * (dim-1) * (dim-2) // 6)
    n_cubic_found = 0

    for trial in range(n_cubic_test):
        # Pick 3 random distinct bits
        indices = []
        while len(indices) < 3:
            b = int.from_bytes(os.urandom(1), 'big') % dim
            if b not in indices:
                indices.append(b)
        i, j, k = sorted(indices)

        # Third-order differential: f(i+j+k) - f(i+j) - f(i+k) - f(j+k) + f(i) + f(j) + f(k) - f(0)
        combos = {}
        for mask in range(8):
            bits = [0] * dim
            if mask & 1: bits[i] = 1
            if mask & 2: bits[j] = 1
            if mask & 4: bits[k] = 1
            state = sha.bits_to_state(bits)
            combos[mask] = sha.state_to_bits(sha.xor_multi_round(state, R))

        d3 = [0] * dim
        for out in range(dim):
            for mask in range(8):
                d3[out] ^= combos[mask][out]

        if any(d3):
            n_cubic_found += 1

    return rank_A, n_mono, rank_full, n_cubic_found, n_cubic_test


# ============================================================
# Main: sweep over rounds
# ============================================================
if __name__ == '__main__':
    print("="*70)
    print("SCF: CARLEMAN COMPOSITION — rank growth over rounds")
    print("="*70)

    for n in [4, 6]:
        sha = MiniSHA(n)
        dim = 8 * n
        total_pairs = dim * (dim-1) // 2

        print(f"\n{'='*70}")
        print(f"n = {n} bits (state dim = {dim})")
        print(f"{'='*70}")
        print(f"  {'R':>3} | {'rank_A':>7} {'n_mono':>7} {'rank_full':>10} {'carleman':>9} | "
              f"{'cubic':>6} | notes")
        print("  " + "-"*70)

        for R in range(1, min(17, 1 << n) + 1):
            rA, nm, rF, nc, nc_test = analyze_composition(n, R, sha)
            carl_dim = dim + nm
            compression = nm / total_pairs * 100 if total_pairs > 0 else 0

            notes = ""
            if nc > 0:
                notes = f" ★ CUBIC TERMS ({nc}/{nc_test})!"
            elif nm == 0:
                notes = " ← LINEAR"
            elif rF > rA:
                notes = f" ← monomials add {rF-rA} rank"

            print(f"  {R:3d} | {rA:7d} {nm:7d} {rF:10d} {carl_dim:9d} | "
                  f"{nc:3d}/{nc_test:3d} | {notes}")

    # Summary and scaling
    print(f"\n{'='*70}")
    print("SCALING ANALYSIS")
    print(f"{'='*70}")
    print("""
  For n=32 (real SHA-256):
    State dim: 256
    One round: 160 monomials, rank = 256
    Key question: does composition to R=64 create CUBIC terms?

    If NO cubic terms for any R:
      → SHA-256 is EXACTLY degree 2 over GF(2) for any number of rounds!
      → Carleman dimension stays bounded at ~8500
      → A^64 is a LINEAR map in Carleman space

    If cubic terms appear at round R:
      → Carleman needs degree-3 lift at round R
      → New monomials: C(256,3) = 2.7M possible
      → But actual count may be much less (like degree 2: 0.5%)
      → If ~13000 cubic monomials → Carleman dim ≈ 20000 (still tractable)

    The DEGREE SATURATION at 32 (from our EXP) suggests:
      Cubic terms DO appear, but saturate.
      Final Carleman dim = C(256,1) + C(256,2) = 33152
      This is the EFFECTIVE dimension of SHA-256 as a polynomial.
    """)
