"""
Q∩T System Builder for SHA-256.

Constructs the hybrid system:
  Q: Quadratic GF(2) equations (SHA-256 with fixed carry → degree 2)
  T: Carry self-consistency constraints (MAJ recursion)

Key insight (methodology_v20.md, Section 212):
  When carry-out vector c is FIXED, addition becomes XOR + constant.
  SHA-256 with fixed carry = quadratic (degree 2 from Ch and Maj only).

System architecture with intermediate state variables:
  - Variables: message bits (512) + state bits per round (256 each)
  - Each round transition equation is at most degree 2
  - Linear part: ~85% of equations (shifts, XOR, additions with fixed carry)
  - Quadratic part: ~15% (Ch and Maj terms, 64 per round)
"""

from qt_solver.sha256_traced import (
    K256, IV256, MASK32,
    sha256_compress_traced, get_bit,
    big_sigma0_bit_sources, big_sigma1_bit_sources,
    sigma0_bit_sources, sigma1_bit_sources,
)


# ─── Variable Indexing ───

class VarMap:
    """
    Maps symbolic variables to flat integer indices.

    Layout:
      [0, 512)              : message bits W[0..15] (free variables)
      [512, 512 + S*32)     : schedule bits W[16..R-1] (defined by schedule)
      [512+S*32, ...)       : state bits state[1..R] per round
                               (state[0] = IV = constants, no variables)

    State register order: a=0, b=1, c=2, d=3, e=4, f=5, g=6, h=7
    """

    def __init__(self, num_rounds):
        self.R = num_rounds
        self._idx = 0

        # Message variables: W[0..15], 32 bits each
        self.msg_base = 0
        self._idx = 512

        # Schedule variables: W[16..R-1] (only if R > 16)
        self.sched_base = self._idx
        self.num_sched_words = max(0, num_rounds - 16)
        self._idx += self.num_sched_words * 32

        # State variables: state[1..R], 8 registers × 32 bits each
        self.state_base = self._idx
        self._idx += num_rounds * 256  # R rounds of state

        self.num_vars = self._idx

    def w(self, word_idx, bit):
        """Variable index for W[word_idx][bit]."""
        assert 0 <= bit < 32
        if word_idx < 16:
            return self.msg_base + word_idx * 32 + bit
        else:
            assert word_idx < self.R, f"W[{word_idx}] out of range for {self.R} rounds"
            return self.sched_base + (word_idx - 16) * 32 + bit

    def s(self, round_idx, reg, bit):
        """Variable index for state[round_idx][reg][bit].
        round_idx: 1..R (state[0] = IV, no variables)
        reg: 0=a, 1=b, 2=c, 3=d, 4=e, 5=f, 6=g, 7=h
        """
        assert 1 <= round_idx <= self.R
        assert 0 <= reg < 8
        assert 0 <= bit < 32
        return self.state_base + (round_idx - 1) * 256 + reg * 32 + bit

    def describe(self, var_idx):
        """Human-readable description of a variable."""
        if var_idx < 512:
            w, b = divmod(var_idx, 32)
            return f"W[{w}][{b}]"
        elif var_idx < self.sched_base + self.num_sched_words * 32:
            off = var_idx - self.sched_base
            w, b = divmod(off, 32)
            return f"W[{w + 16}][{b}]"
        else:
            off = var_idx - self.state_base
            r, rem = divmod(off, 256)
            reg, b = divmod(rem, 32)
            reg_names = 'abcdefgh'
            return f"state[{r + 1}].{reg_names[reg]}[{b}]"


# ─── Equation Representation ───

class GF2Equation:
    """
    A GF(2) equation: constant ⊕ Σ(linear vars) ⊕ Σ(quadratic products) = 0

    Linear part stored as a set of variable indices (XOR = symmetric difference).
    Quadratic part stored as a set of (i,j) tuples with i < j.
    """
    __slots__ = ('constant', 'linear', 'quadratic')

    def __init__(self, constant=0, linear=None, quadratic=None):
        self.constant = constant & 1
        self.linear = set(linear) if linear else set()
        self.quadratic = set(quadratic) if quadratic else set()

    def is_linear(self):
        return len(self.quadratic) == 0

    def xor(self, other):
        """XOR two equations (addition in GF(2))."""
        return GF2Equation(
            constant=self.constant ^ other.constant,
            linear=self.linear ^ other.linear,
            quadratic=self.quadratic ^ other.quadratic,
        )

    def add_var(self, v):
        """Toggle variable v in linear part."""
        self.linear ^= {v}

    def add_product(self, v1, v2):
        """Toggle product v1*v2 in quadratic part."""
        if v1 == v2:
            # x*x = x over GF(2)
            self.linear ^= {v1}
            return
        if v1 > v2:
            v1, v2 = v2, v1
        self.quadratic ^= {(v1, v2)}

    def add_const(self, c):
        """Toggle constant."""
        self.constant ^= (c & 1)

    def clone(self):
        return GF2Equation(self.constant, set(self.linear), set(self.quadratic))

    def evaluate(self, assignment):
        """
        Evaluate equation given a variable assignment (dict var_id -> 0/1).
        Returns 0 if satisfied, 1 if violated.
        """
        val = self.constant
        for v in self.linear:
            val ^= assignment.get(v, 0)
        for v1, v2 in self.quadratic:
            val ^= (assignment.get(v1, 0) & assignment.get(v2, 0))
        return val

    def __repr__(self):
        parts = []
        if self.constant:
            parts.append('1')
        for v in sorted(self.linear):
            parts.append(f'x{v}')
        for v1, v2 in sorted(self.quadratic):
            parts.append(f'x{v1}*x{v2}')
        return ' ⊕ '.join(parts) + ' = 0' if parts else '0 = 0'


# ─── Q∩T System ───

class QTSystem:
    """
    The Q∩T system for SHA-256.

    Q equations: SHA-256 round transitions with fixed carry (quadratic GF(2))
    T equations: carry self-consistency (carry[i] = MAJ(a[i], b[i], carry[i-1]))

    With intermediate state variables, each equation is at most degree 2.
    """

    def __init__(self, num_rounds, carry_chains, vm=None):
        """
        Args:
            num_rounds: number of SHA-256 rounds
            carry_chains: list of 32-bit carry chain integers
                          (from sha256_compress_traced)
            vm: VarMap instance (created if None)
        """
        self.R = num_rounds
        self.carry_chains = carry_chains
        self.vm = vm or VarMap(num_rounds)

        self.q_linear = []     # Linear Q equations
        self.q_quadratic = []  # Quadratic Q equations
        self.t_equations = []  # T (carry self-consistency) equations

        self._carry_idx = 0  # Index into carry_chains list

    def _next_carry_chain(self):
        """Get next carry chain from the list."""
        c = self.carry_chains[self._carry_idx]
        self._carry_idx += 1
        return c

    def build(self, target_hash=None):
        """
        Build the complete Q∩T system.

        Args:
            target_hash: 8 uint32 words (for preimage); None for structure analysis
        """
        self._carry_idx = 0
        self.q_linear = []
        self.q_quadratic = []
        self.t_equations = []

        # 1. Schedule equations (for W[16..R-1])
        self._build_schedule_equations()

        # 2. Round transition equations
        self._build_round_equations()

        # 3. Output equations (if target given)
        if target_hash is not None:
            self._build_output_equations(target_hash)

    def _get_state_bits(self, round_idx, reg):
        """
        Get variable indices for state[round_idx][reg][0..31].
        For round 0: returns None (use IV constants).
        For round > 0: returns list of variable indices.
        """
        if round_idx == 0:
            return None  # IV constants
        return [self.vm.s(round_idx, reg, b) for b in range(32)]

    def _get_state_val(self, round_idx, reg):
        """Get the constant value for state[0] (IV)."""
        if round_idx != 0:
            raise ValueError("Only round 0 has constant state")
        return IV256[reg]  # reg 0..7 maps to a..h from IV

    def _build_schedule_equations(self):
        """
        Build schedule expansion equations for W[16..R-1].
        W[i] = W[i-16] + sigma0(W[i-15]) + W[i-7] + sigma1(W[i-2])

        With fixed carry, each addition is linear:
        result[b] = op1[b] ⊕ op2[b] ⊕ carry_into[b]
        """
        for i in range(16, self.R):
            # 3 additions in schedule expansion
            # add1: W[i-16] + sigma0(W[i-15])
            # add2: add1_result + W[i-7]
            # add3: add2_result + sigma1(W[i-2])
            # Final result = W[i]

            # We track intermediate results symbolically
            # Since all additions are linear (carry fixed), we can compose them

            # For each output bit b of W[i]:
            for b in range(32):
                eq = GF2Equation()

                # The output variable: W[i][b]
                eq.add_var(self.vm.w(i, b))

                # Trace through 3 additions:
                # add3_result[b] = add2_result[b] ⊕ sigma1(W[i-2])[b] ⊕ carry3[b]
                # add2_result[b] = add1_result[b] ⊕ W[i-7][b] ⊕ carry2[b]
                # add1_result[b] = W[i-16][b] ⊕ sigma0(W[i-15])[b] ⊕ carry1[b]

                # Get carry chains for this schedule word
                sched_idx = i - 16
                carries_for_word = self.carry_chains[sched_idx * 3: sched_idx * 3 + 3]

                c1 = carries_for_word[0]
                c2 = carries_for_word[1]
                c3 = carries_for_word[2]

                # Combined constant from carries
                carry_const = get_bit(c1, b) ^ get_bit(c2, b) ^ get_bit(c3, b)
                eq.add_const(carry_const)

                # W[i-16][b]
                eq.add_var(self.vm.w(i - 16, b))

                # sigma0(W[i-15])[b] = XOR of specific bits of W[i-15]
                for src_bit in sigma0_bit_sources(b):
                    eq.add_var(self.vm.w(i - 15, src_bit))

                # W[i-7][b]
                eq.add_var(self.vm.w(i - 7, b))

                # sigma1(W[i-2])[b] = XOR of specific bits of W[i-2]
                for src_bit in sigma1_bit_sources(b):
                    eq.add_var(self.vm.w(i - 2, src_bit))

                self.q_linear.append(eq)

    def _build_round_equations(self):
        """
        Build round transition equations for rounds 0..R-1.

        For each round r:
          state[r+1] is determined from state[r] and W[r].
          With fixed carry:
            - Shift register copies (b=a, c=b, d=c, f=e, g=f, h=g): LINEAR
            - e_new = d + T1: LINEAR (carry fixed) but T1 contains Ch → QUADRATIC
            - a_new = T1 + T2: LINEAR (carry fixed) but T1 has Ch, T2 has Maj → QUADRATIC
        """
        # Carry index offset: skip schedule carries
        num_sched_carries = max(0, self.R - 16) * 3
        carry_base = num_sched_carries

        for r in range(self.R):
            rc_base = carry_base + r * 7  # 7 carry chains per round

            # Get carry chains for this round's 7 additions
            c_T1_1 = self.carry_chains[rc_base + 0]  # h + Sigma1(e)
            c_T1_2 = self.carry_chains[rc_base + 1]  # + Ch(e,f,g)
            c_T1_3 = self.carry_chains[rc_base + 2]  # + K[r]
            c_T1_4 = self.carry_chains[rc_base + 3]  # + W[r]
            c_T2   = self.carry_chains[rc_base + 4]  # Sigma0(a) + Maj
            c_e    = self.carry_chains[rc_base + 5]  # d + T1
            c_a    = self.carry_chains[rc_base + 6]  # T1 + T2

            # ── Shift register equations (LINEAR) ──
            # b[r+1] = a[r], c[r+1] = b[r], d[r+1] = c[r]
            # f[r+1] = e[r], g[r+1] = f[r], h[r+1] = g[r]
            shift_pairs = [
                (1, 0),  # b[r+1] = a[r]
                (2, 1),  # c[r+1] = b[r]
                (3, 2),  # d[r+1] = c[r]
                (5, 4),  # f[r+1] = e[r]
                (6, 5),  # g[r+1] = f[r]
                (7, 6),  # h[r+1] = g[r]
            ]
            for dst_reg, src_reg in shift_pairs:
                for b in range(32):
                    eq = GF2Equation()
                    # state[r+1][dst_reg][b] = state[r][src_reg][b]
                    eq.add_var(self.vm.s(r + 1, dst_reg, b))
                    if r == 0:
                        # state[0] is IV constant
                        eq.add_const(get_bit(IV256[src_reg], b))
                    else:
                        eq.add_var(self.vm.s(r, src_reg, b))
                    self.q_linear.append(eq)

            # ── e_new and a_new equations (QUADRATIC from Ch and Maj) ──
            for b in range(32):
                self._build_e_equation(r, b, c_T1_1, c_T1_2, c_T1_3, c_T1_4, c_e)
                self._build_a_equation(r, b, c_T1_1, c_T1_2, c_T1_3, c_T1_4,
                                       c_T2, c_a)

    def _build_e_equation(self, r, b, c1, c2, c3, c4, c_e):
        """
        Build equation for e[r+1][b].
        e_new = d + T1, where T1 = h + Sigma1(e) + Ch(e,f,g) + K[r] + W[r]

        With fixed carry:
        e_new[b] = d[b] ⊕ T1[b] ⊕ carry_e[b]
        T1[b] = h[b] ⊕ Sigma1(e)[b] ⊕ Ch(e,f,g)[b] ⊕ K[r][b] ⊕ W[r][b]
                 ⊕ carry_T1_combined[b]

        Ch(e,f,g)[b] = e[b]*f[b] ⊕ e[b]*g[b] ⊕ g[b]  (QUADRATIC)
        """
        eq = GF2Equation()

        # LHS: e[r+1][b]
        eq.add_var(self.vm.s(r + 1, 4, b))  # reg 4 = e

        # Combined carry constant from all 5 additions
        carry_const = get_bit(c1, b) ^ get_bit(c2, b) ^ get_bit(c3, b) \
                      ^ get_bit(c4, b) ^ get_bit(c_e, b)
        eq.add_const(carry_const)

        # K[r][b] constant
        eq.add_const(get_bit(K256[r], b))

        # Helper to add state variable or IV constant
        def add_state_var(reg, bit):
            if r == 0:
                eq.add_const(get_bit(IV256[reg], bit))
            else:
                eq.add_var(self.vm.s(r, reg, bit))

        # d[r][b]
        add_state_var(3, b)  # reg 3 = d

        # h[r][b]
        add_state_var(7, b)  # reg 7 = h

        # Sigma1(e[r])[b] = XOR of e[r] bits
        for src_bit in big_sigma1_bit_sources(b):
            add_state_var(4, src_bit)  # reg 4 = e

        # W[r][b]
        eq.add_var(self.vm.w(r, b))

        # Ch(e,f,g)[b] = e[b]*f[b] ⊕ e[b]*g[b] ⊕ g[b]
        if r == 0:
            # state[0] is constant → Ch is constant
            e_val = get_bit(IV256[4], b)
            f_val = get_bit(IV256[5], b)
            g_val = get_bit(IV256[6], b)
            ch_val = (e_val & f_val) ^ (e_val & g_val) ^ g_val
            eq.add_const(ch_val)
            self.q_linear.append(eq)
        else:
            # g[r][b] (linear term of Ch)
            eq.add_var(self.vm.s(r, 6, b))

            # e[b]*f[b] (quadratic)
            eq.add_product(self.vm.s(r, 4, b), self.vm.s(r, 5, b))

            # e[b]*g[b] (quadratic)
            eq.add_product(self.vm.s(r, 4, b), self.vm.s(r, 6, b))

            self.q_quadratic.append(eq)

    def _build_a_equation(self, r, b, c1, c2, c3, c4, c_T2, c_a):
        """
        Build equation for a[r+1][b].
        a_new = T1 + T2
        T1 = h + Sigma1(e) + Ch(e,f,g) + K[r] + W[r]
        T2 = Sigma0(a) + Maj(a,b,c)

        Ch contributes quadratic terms, Maj contributes quadratic terms.
        """
        eq = GF2Equation()

        # LHS: a[r+1][b]
        eq.add_var(self.vm.s(r + 1, 0, b))  # reg 0 = a

        # Combined carry from T1 (4 adds) + T2 (1 add) + a_new (1 add)
        carry_const = get_bit(c1, b) ^ get_bit(c2, b) ^ get_bit(c3, b) \
                      ^ get_bit(c4, b) ^ get_bit(c_T2, b) ^ get_bit(c_a, b)
        eq.add_const(carry_const)

        # K[r][b]
        eq.add_const(get_bit(K256[r], b))

        def add_state_var(reg, bit):
            if r == 0:
                eq.add_const(get_bit(IV256[reg], bit))
            else:
                eq.add_var(self.vm.s(r, reg, bit))

        # T1 linear terms: h[b], Sigma1(e)[b], W[r][b]
        add_state_var(7, b)
        for src_bit in big_sigma1_bit_sources(b):
            add_state_var(4, src_bit)
        eq.add_var(self.vm.w(r, b))

        # T2 linear term: Sigma0(a)[b]
        for src_bit in big_sigma0_bit_sources(b):
            add_state_var(0, src_bit)

        if r == 0:
            # All state is IV → Ch and Maj are constants
            e_val = get_bit(IV256[4], b)
            f_val = get_bit(IV256[5], b)
            g_val = get_bit(IV256[6], b)
            ch_val = (e_val & f_val) ^ (e_val & g_val) ^ g_val
            eq.add_const(ch_val)

            a_val = get_bit(IV256[0], b)
            b_val = get_bit(IV256[1], b)
            c_val = get_bit(IV256[2], b)
            maj_val = (a_val & b_val) ^ (a_val & c_val) ^ (b_val & c_val)
            eq.add_const(maj_val)

            self.q_linear.append(eq)
        else:
            # Ch(e,f,g)[b] = e*f ⊕ e*g ⊕ g
            eq.add_var(self.vm.s(r, 6, b))  # g[b] linear
            eq.add_product(self.vm.s(r, 4, b), self.vm.s(r, 5, b))  # e*f
            eq.add_product(self.vm.s(r, 4, b), self.vm.s(r, 6, b))  # e*g

            # Maj(a,b,c)[b] = a*b ⊕ a*c ⊕ b*c
            eq.add_product(self.vm.s(r, 0, b), self.vm.s(r, 1, b))  # a*b
            eq.add_product(self.vm.s(r, 0, b), self.vm.s(r, 2, b))  # a*c
            eq.add_product(self.vm.s(r, 1, b), self.vm.s(r, 2, b))  # b*c

            self.q_quadratic.append(eq)

    def _build_output_equations(self, target_hash):
        """
        Build output equations: H[i] = state[R][i] + IV[i] = target[i]
        → state[R][i] = target[i] - IV[i] mod 2^32

        With fixed carry on the final addition:
        state[R][reg][b] ⊕ IV[reg][b] ⊕ carry_final[reg][b] = target[b]
        """
        num_sched_carries = max(0, self.R - 16) * 3
        num_round_carries = self.R * 7
        final_base = num_sched_carries + num_round_carries

        for reg in range(8):
            carry_chain = self.carry_chains[final_base + reg]
            for b in range(32):
                eq = GF2Equation()

                # state[R][reg][b]
                eq.add_var(self.vm.s(self.R, reg, b))

                # ⊕ IV[reg][b] (constant)
                eq.add_const(get_bit(IV256[reg], b))

                # ⊕ carry[b] (constant, from fixed carry)
                eq.add_const(get_bit(carry_chain, b))

                # = target[reg][b]
                eq.add_const(get_bit(target_hash[reg], b))

                self.q_linear.append(eq)

    def build_t_constraints(self, msg_words, trace):
        """
        Build T (carry self-consistency) constraints.

        For each addition's carry chain, verify:
          carry[0] = 0 (always)
          carry[i] = MAJ(a[i], b[i], carry[i-1]) for i = 1..31

        These are quadratic GF(2) equations:
          carry[i] = a[i]*b[i] ⊕ a[i]*carry[i-1] ⊕ b[i]*carry[i-1]

        Returns: number of violated T constraints
        """
        violations = 0
        all_chains = []

        for sc in trace.sched_carries:
            for chain, _ in sc:
                all_chains.append(chain)
        for rc in trace.round_carries:
            for chain, _ in rc:
                all_chains.append(chain)
        for chain, _ in trace.final_carries:
            all_chains.append(chain)

        assumed_chains = self.carry_chains

        for idx in range(min(len(all_chains), len(assumed_chains))):
            if all_chains[idx] != assumed_chains[idx]:
                # Count differing bits
                diff = all_chains[idx] ^ assumed_chains[idx]
                violations += bin(diff).count('1')

        return violations

    def stats(self):
        """Return system statistics."""
        return {
            'num_rounds': self.R,
            'num_vars': self.vm.num_vars,
            'num_msg_vars': 512,
            'num_sched_vars': self.vm.num_sched_words * 32,
            'num_state_vars': self.R * 256,
            'num_linear_eq': len(self.q_linear),
            'num_quadratic_eq': len(self.q_quadratic),
            'total_eq': len(self.q_linear) + len(self.q_quadratic),
            'quad_terms_total': sum(len(eq.quadratic) for eq in self.q_quadratic),
        }

    def to_linear_rows(self):
        """
        Convert linear equations to GF(2) augmented row format
        for Gaussian elimination.

        Returns list of integers (bit vectors).
        """
        rows = []
        for eq in self.q_linear:
            row = 0
            for v in eq.linear:
                row |= (1 << v)
            if eq.constant:
                row |= (1 << self.vm.num_vars)  # augmented column
            rows.append(row)
        return rows


def build_qt_system(msg_words, num_rounds=64, target_hash=None):
    """
    Convenience function: compute SHA-256 trace and build Q∩T system.

    Args:
        msg_words: 16 uint32 message words
        num_rounds: number of rounds
        target_hash: target hash for preimage (None = use computed hash)

    Returns:
        (system, trace) tuple
    """
    from qt_solver.sha256_traced import sha256_compress_traced, get_all_carry_chains

    trace = sha256_compress_traced(msg_words, num_rounds)

    if target_hash is None:
        target_hash = trace.hash_words

    carry_chains = get_all_carry_chains(trace)

    system = QTSystem(num_rounds, carry_chains)
    system.build(target_hash)

    return system, trace
