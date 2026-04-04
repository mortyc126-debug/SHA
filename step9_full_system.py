"""
Step 9: Build and Solve the FULL GF(2) System

Stop counting, start solving.

Build the complete collision system for mini-SHA:
  - Keep ALL intermediate variables (carries, state bits)
  - Write ALL equations (degree-2 and linear)
  - Solve by: Gaussian elimination → back-substitution → enumerate kernel

Compare solve cost vs brute force (2^16) and birthday (2^4).

If the system ACTUALLY solves faster → we have an algorithm.
If not → we see EXACTLY where the bottleneck is.
"""

import numpy as np
import time
from step0_exact_algebra import mini_sha, N, MASK

N_MSG = 4
N_INPUT = N * N_MSG  # 16
N_TOTAL = 1 << N_INPUT


class GF2System:
    """
    System of GF(2) equations with variables.
    Supports linear and degree-2 (quadratic) equations.

    For quadratic: linearize by introducing y_{ij} = x_i * x_j.
    """

    def __init__(self):
        self.var_names = []      # variable name → index
        self.var_index = {}      # name → index
        self.n_vars = 0
        self.equations = []      # list of (set_of_var_indices, constant)
        # For quadratic: product variables
        self.product_vars = {}   # (i,j) → var index (i < j)

    def add_var(self, name):
        """Add a variable, return its index."""
        if name in self.var_index:
            return self.var_index[name]
        idx = self.n_vars
        self.var_names.append(name)
        self.var_index[name] = idx
        self.n_vars += 1
        return idx

    def get_product_var(self, i, j):
        """Get or create linearized variable for x_i * x_j."""
        if i > j:
            i, j = j, i
        if (i, j) not in self.product_vars:
            name = f"({self.var_names[i]}*{self.var_names[j]})"
            idx = self.add_var(name)
            self.product_vars[(i, j)] = idx
        return self.product_vars[(i, j)]

    def add_equation(self, terms, constant=0):
        """
        Add equation: XOR of terms = constant (mod 2).
        terms: list of variable indices (linear) or tuples (i,j) for products.
        """
        var_set = set()
        for t in terms:
            if isinstance(t, tuple):
                # Product term → linearize
                i, j = t
                if i == j:
                    # x_i * x_i = x_i over GF(2)
                    var_set ^= {i}
                else:
                    prod_var = self.get_product_var(i, j)
                    var_set ^= {prod_var}
            else:
                var_set ^= {t}
        self.equations.append((frozenset(var_set), constant))

    def to_matrix(self):
        """Convert to augmented matrix [A | b] over GF(2)."""
        n_eqs = len(self.equations)
        # Matrix: n_eqs × (n_vars + 1), last column = constants
        M = np.zeros((n_eqs, self.n_vars + 1), dtype=np.uint8)
        for i, (vars_set, const) in enumerate(self.equations):
            for v in vars_set:
                M[i][v] = 1
            M[i][self.n_vars] = const
        return M

    def solve(self):
        """
        Gaussian elimination over GF(2).
        Returns: rank, kernel_dimension, and whether system is consistent.
        """
        M = self.to_matrix()
        n_rows, n_cols = M.shape
        n_vars = n_cols - 1

        # Forward elimination
        pivot_col = [None] * n_rows
        rank = 0

        for col in range(n_vars):
            # Find pivot
            pivot_row = None
            for row in range(rank, n_rows):
                if M[row][col]:
                    pivot_row = row
                    break
            if pivot_row is None:
                continue

            # Swap
            M[[rank, pivot_row]] = M[[pivot_row, rank]]
            pivot_col[rank] = col

            # Eliminate
            for row in range(n_rows):
                if row != rank and M[row][col]:
                    M[row] ^= M[rank]

            rank += 1

        # Check consistency
        consistent = True
        for row in range(rank, n_rows):
            if M[row][n_vars]:  # 0 = 1 → inconsistent
                consistent = False
                break

        kernel_dim = n_vars - rank

        return rank, kernel_dim, consistent, M[:rank]


def build_collision_system(M_base, R):
    """
    Build the FULL GF(2) system for finding collision δ.

    For each round r: state(M⊕δ) vs state(M).
    The difference propagation creates equations in δ-variables.
    """
    sys = GF2System()

    # Message difference variables
    dm = []
    for w in range(N_MSG):
        word_vars = []
        for b in range(N):
            idx = sys.add_var(f"dM[{w}][{b}]")
            word_vars.append(idx)
        dm.append(word_vars)

    # We need to track the ACTUAL computation to get state-dependent coefficients
    # Run the base computation to get all intermediate values
    states = []  # states[r] = (a,b,c,d,e,f,g,h) as integers

    a, b, c, d, e, f, g, h = 0x6, 0xB, 0x3, 0xA, 0x5, 0x9, 0x1, 0xF
    W = list(M_base) + [0] * max(0, R - N_MSG)
    K = [0x4, 0x2, 0xB, 0x7, 0xA, 0x3, 0xE, 0x5,
         0x9, 0x1, 0xD, 0x6, 0x0, 0x8, 0xC, 0xF]

    states.append((a, b, c, d, e, f, g, h))

    for r in range(R):
        w_r = W[r] if r < len(W) else 0
        k_r = K[r % len(K)]

        def rotr(x, s):
            return ((x >> s) | (x << (N - s))) & MASK

        sig1 = rotr(e, 1) ^ rotr(e, 3) ^ (e >> 1)
        ch_val = (e & f) ^ (~e & g) & MASK
        sig0 = rotr(a, 1) ^ rotr(a, 2) ^ rotr(a, 3)
        maj_val = (a & b) ^ (a & c) ^ (b & c)

        T1_val = (h + sig1 + ch_val + k_r + w_r) & MASK
        T2_val = (sig0 + maj_val) & MASK
        a_new = (T1_val + T2_val) & MASK
        e_new = (d + T1_val) & MASK

        h, g, f = g, f, e
        e = e_new
        d, c, b = c, b, a
        a = a_new

        states.append((a, b, c, d, e, f, g, h))

    # Now build the difference system
    # δstate[0] = 0 (same IV)
    # For each round: propagate differences

    # Current δstate variables (initialized to 0 = no variable)
    # da[bit], de[bit] etc. — but b=a_prev, so db = da_prev
    da = [None] * N  # None = known zero
    de = [None] * N

    # Track all δ-variables through rounds
    total_degree2_eqs = 0
    total_linear_eqs = 0

    for r in range(R):
        # Get base state at round r
        a_val, b_val, c_val, d_val, e_val, f_val, g_val, h_val = states[r]
        w_r = W[r] if r < len(W) else 0
        k_r = K[r % len(K)]

        # δW[r] = message difference for word r (or 0 if r >= N_MSG)
        if r < N_MSG:
            dw = dm[r]  # list of N variable indices
        else:
            dw = [None] * N  # zero

        # δh = δe from 3 rounds ago (shift register)
        # For r=0: δh = 0 (no difference in IV)
        # General: δh = δg_prev = δf_prev2 = δe_prev3
        dh = de if r == 0 else prev_de_3  # Will track this

        # For simplicity in this first version:
        # Track δa and δe through rounds using shift register
        # da_prev[k] = da[k] from round r-k for b,c,d
        # de_prev[k] = de[k] from round r-k for f,g,h

        # Actually, let me use a simpler approach:
        # maintain the full 8-register δ-state

        pass  # Complex... let me use truth tables instead

    # SIMPLER APPROACH: Use truth tables to build the system
    # For each output bit, the collision polynomial gives us the equation
    # But we want INTERMEDIATE equations too

    return sys  # placeholder


def direct_solve_with_truth_tables(M_base, R):
    """
    More direct approach: compute the collision polynomial for EACH
    intermediate state bit, not just the output. This gives us
    intermediate equations "for free".

    Each intermediate bit at round r is a polynomial in δM.
    The ROUND FUNCTION gives degree-2 constraints between consecutive rounds.
    """
    from step0_exact_algebra import mobius_transform as mobius16

    print(f"\n{'='*80}")
    print(f"FULL SYSTEM SOLVE — R={R}, M={[hex(w) for w in M_base]}")
    print(f"{'='*80}")

    # Step 1: Compute truth tables for ALL intermediate δ-bits
    # For each δ (0..2^16-1), run both M and M⊕δ through R rounds,
    # record δ at each intermediate point

    a_base, e_base = mini_sha(M_base, R)

    # Brute force: find all collisions (ground truth)
    t0 = time.time()
    collisions = []
    for didx in range(1, N_TOTAL):
        dw = []
        tmp = didx
        for w in range(N_MSG):
            dw.append(tmp & MASK)
            tmp >>= N
        M2 = [(M_base[w] ^ dw[w]) for w in range(N_MSG)]
        a2, e2 = mini_sha(M2, R)
        if a2 == a_base and e2 == e_base:
            collisions.append(didx)
    t_brute = time.time() - t0
    print(f"\n  Brute force: {len(collisions)} collisions in {t_brute:.3f}s")

    # Step 2: For each INTERMEDIATE round, compute δa and δe as truth tables in δ
    # Then extract degree-2 part and see if it constrains the system

    # Compute all intermediate states for M_base
    def run_rounds(msg, R):
        a, b, c, d, e, f, g, h = 0x6, 0xB, 0x3, 0xA, 0x5, 0x9, 0x1, 0xF
        W = list(msg) + [0] * max(0, R - len(msg))
        K = [0x4, 0x2, 0xB, 0x7, 0xA, 0x3, 0xE, 0x5,
             0x9, 0x1, 0xD, 0x6, 0x0, 0x8, 0xC, 0xF]
        history = [(a, e)]
        for r in range(R):
            def rotr(x, s):
                return ((x >> s) | (x << (N - s))) & MASK
            w_r = W[r] if r < len(W) else 0
            k_r = K[r % len(K)]
            sig1 = rotr(e, 1) ^ rotr(e, 3) ^ (e >> 1)
            ch_val = (e & f) ^ (~e & g) & MASK
            sig0 = rotr(a, 1) ^ rotr(a, 2) ^ rotr(a, 3)
            maj_val = (a & b) ^ (a & c) ^ (b & c)
            T1 = (h + sig1 + ch_val + k_r + w_r) & MASK
            T2 = (sig0 + maj_val) & MASK
            a_new = (T1 + T2) & MASK
            e_new = (d + T1) & MASK
            h, g, f = g, f, e
            e = e_new
            d, c, b = c, b, a
            a = a_new
            history.append((a, e))
        return history

    base_history = run_rounds(M_base, R)

    # For each δ, compute difference at each round
    # Build truth table of δa[r][bit] and δe[r][bit] as functions of δ
    intermediate_tts = {}  # (round, 'a'/'e', bit) → truth table

    for r in range(R + 1):
        for bit in range(N):
            intermediate_tts[(r, 'a', bit)] = np.zeros(N_TOTAL, dtype=np.uint8)
            intermediate_tts[(r, 'e', bit)] = np.zeros(N_TOTAL, dtype=np.uint8)

    for didx in range(N_TOTAL):
        dw = []
        tmp = didx
        for w in range(N_MSG):
            dw.append(tmp & MASK)
            tmp >>= N
        M2 = [(M_base[w] ^ dw[w]) for w in range(N_MSG)]
        hist2 = run_rounds(M2, R)

        for r in range(R + 1):
            da = base_history[r][0] ^ hist2[r][0]
            de = base_history[r][1] ^ hist2[r][1]
            for bit in range(N):
                intermediate_tts[(r, 'a', bit)][didx] = (da >> bit) & 1
                intermediate_tts[(r, 'e', bit)][didx] = (de >> bit) & 1

    # Step 3: For each intermediate δ-bit, compute ANF and get its degree
    print(f"\n  INTERMEDIATE POLYNOMIAL DEGREES:")
    print(f"  {'Round':>6} {'da degrees':>30} {'de degrees':>30}")

    for r in range(R + 1):
        da_degs = []
        de_degs = []
        for bit in range(N):
            anf_a = mobius16(intermediate_tts[(r, 'a', bit)], N_INPUT)
            anf_e = mobius16(intermediate_tts[(r, 'e', bit)], N_INPUT)
            deg_a = max((bin(m).count('1') for m in range(N_TOTAL) if anf_a[m]), default=0)
            deg_e = max((bin(m).count('1') for m in range(N_TOTAL) if anf_e[m]), default=0)
            da_degs.append(deg_a)
            de_degs.append(deg_e)
        print(f"  R={r:>3}  a:{da_degs}  e:{de_degs}")

    # Step 4: BUILD LINEAR SYSTEM from degree-2 truncation of ALL rounds
    print(f"\n  BUILDING LINEARIZED SYSTEM (degree ≤ 2, all rounds)...")

    # Enumerate all degree-≤2 monomials
    mono_list = [0]  # constant
    for i in range(N_INPUT):
        mono_list.append(1 << i)
    for i in range(N_INPUT):
        for j in range(i+1, N_INPUT):
            mono_list.append((1 << i) | (1 << j))

    n_monos = len(mono_list)
    mono_to_col = {m: i for i, m in enumerate(mono_list)}

    # For each intermediate δ-bit at each round: extract degree-≤2 ANF
    rows = []
    for r in range(1, R + 1):  # skip r=0 (all zero)
        for reg in ['a', 'e']:
            for bit in range(N):
                anf = mobius16(intermediate_tts[(r, reg, bit)], N_INPUT)
                # Extract degree-≤2 monomials
                row = np.zeros(n_monos + 1, dtype=np.uint8)  # +1 for constant
                has_higher = False
                for m_idx in range(N_TOTAL):
                    if anf[m_idx]:
                        deg = bin(m_idx).count('1')
                        if deg <= 2 and m_idx in mono_to_col:
                            row[mono_to_col[m_idx]] = 1
                        elif deg > 2:
                            has_higher = True
                        elif m_idx == 0:
                            row[0] = 1

                # For OUTPUT round: equation is δ-bit = 0
                if r == R:
                    # This is the actual constraint
                    rows.append(row)
                else:
                    # Intermediate: provides structure but isn't a "target = 0" constraint
                    # We can use it to define auxiliary variables
                    rows.append(row)

    matrix = np.array(rows, dtype=np.uint8)
    n_eqs_total = matrix.shape[0]

    print(f"    Total equations (degree ≤ 2): {n_eqs_total}")
    print(f"    Linearized variables: {n_monos}")

    # Gaussian elimination
    t0 = time.time()
    M = matrix.copy()
    rank = 0
    for col in range(n_monos):
        pivot = None
        for row in range(rank, n_eqs_total):
            if M[row][col]:
                pivot = row
                break
        if pivot is None:
            continue
        M[[rank, pivot]] = M[[pivot, rank]]
        for row in range(n_eqs_total):
            if row != rank and M[row][col]:
                M[row] ^= M[rank]
        rank += 1
    t_gauss = time.time() - t0

    kernel = n_monos - rank
    print(f"    Rank: {rank}")
    print(f"    Kernel: {kernel}")
    print(f"    Gauss time: {t_gauss:.3f}s")

    # How many solutions?
    if kernel < 30:
        predicted = 2 ** kernel
    else:
        predicted = f"2^{kernel}"

    print(f"\n{'='*80}")
    print(f"  RESULTS COMPARISON")
    print(f"{'='*80}")
    print(f"  Actual collisions:           {len(collisions)}")
    print(f"  Degree-2 linearized kernel:  2^{kernel} = {predicted}")
    print(f"  Brute force cost:            {N_TOTAL} evaluations ({t_brute:.3f}s)")
    print(f"  Gauss cost:                  {n_monos}×{n_eqs_total} matrix ({t_gauss:.3f}s)")
    print(f"  Birthday cost:               ~{2**N} evaluations")

    # The key question: does the degree-2 system ACTUALLY find the solutions?
    if kernel < 20:
        print(f"\n  ENUMERATING kernel solutions (2^{kernel} = {2**kernel})...")
        # Find actual solutions in the kernel
        # For now: just report the numbers
        pass

    return len(collisions), rank, kernel


def main():
    import random
    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(N_MSG)]

    for R in [3, 4, 6, 8]:
        direct_solve_with_truth_tables(M, R)


if __name__ == "__main__":
    main()
