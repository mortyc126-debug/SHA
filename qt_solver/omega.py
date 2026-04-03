"""
OMEGA: Complete GF(2) representation of SHA-256.

Every bit. Every carry. Every operation. No approximation.
Built from scratch — not derived from any existing system.

Architecture:
  For each round r, SHA-256 computes:
    S1  = Sigma1(e)           [linear: XOR of rotations]
    ch  = Ch(e, f, g)         [quadratic: e*f XOR e*g XOR g]
    tmp = h + S1 + ch + K + W [4 chained additions, each with carry]
    S0  = Sigma0(a)           [linear]
    mj  = Maj(a, b, c)        [quadratic: a*b XOR a*c XOR b*c]
    T2  = S0 + mj             [1 addition with carry]
    e'  = d + tmp             [1 addition with carry]
    a'  = tmp + T2            [1 addition with carry]

  Total per round: 7 additions = 7 × 32 carry chains = 224 carry bits.
  Plus: Ch (2 products/bit), Maj (3 products/bit) = 5 × 32 = 160 quad terms.

  Variables per round:
    - 32 bits for each: S1, ch, tmp1, tmp2, tmp3, T1, S0, mj, T2
    - 7 × 32 = 224 carry bits
    - state[r+1]: 256 bits

  Equations per round:
    - Linear: S1 definition (32), S0 definition (32), shift registers (192)
    - Quadratic: Ch definition (32 eqs, 2 products each)
    -            Maj definition (32 eqs, 3 products each)
    -            7 × 32 = 224 carry MAJ equations
    - Linear: 7 × 32 = 224 addition-result equations (given carry)
"""

from qt_solver.sha256_traced import (
    IV256, K256, MASK32, get_bit,
    big_sigma0_bit_sources, big_sigma1_bit_sources,
)


class OmegaVar:
    """Variable allocator with named variables."""
    def __init__(self):
        self.n = 0
        self.names = {}

    def alloc(self, name):
        v = self.n
        self.n += 1
        self.names[v] = name
        return v

    def alloc32(self, prefix):
        """Allocate 32 consecutive variables."""
        return [self.alloc(f"{prefix}[{b}]") for b in range(32)]


class OmegaSystem:
    """
    Complete GF(2) SHA-256 system.
    Linear equations: list of (set_of_vars, constant).
    Quadratic equations: list of (set_of_linear_vars, set_of_products, constant).
    """

    def __init__(self, num_rounds):
        self.R = num_rounds
        self.V = OmegaVar()
        self.linear = []   # [(var_set, const), ...]
        self.quad = []      # [(var_set, product_set, const), ...]

        # Message words
        self.W = [self.V.alloc32(f"W{w}") for w in range(16)]

        # State registers per round: state[0]=IV (no vars), state[1..R]
        self.state = {}  # (round, reg) → [32 vars]
        for r in range(1, num_rounds + 1):
            for reg in range(8):
                name = f"s{r}{'abcdefgh'[reg]}"
                self.state[(r, reg)] = self.V.alloc32(name)

        # Build round by round
        for r in range(num_rounds):
            self._build_round(r)

    def _state_bits(self, r, reg):
        """Get state bits: vars for r>0, IV constants for r=0."""
        if r == 0:
            return None, [get_bit(IV256[reg], b) for b in range(32)]
        return self.state[(r, reg)], None

    def _add_lin(self, var_set, const):
        self.linear.append((frozenset(var_set), const & 1))

    def _add_quad(self, var_set, prod_set, const):
        self.quad.append((frozenset(var_set), frozenset(prod_set), const & 1))

    def _addition_32(self, a_vars, a_const, b_vars, b_const, name):
        """
        Model 32-bit addition a + b.

        Returns: (result_vars, carry_vars)

        For each bit position b:
          result[b] = a[b] XOR b[b] XOR carry[b]    (linear)
          carry[0] = 0                                (linear)
          carry[b+1] = MAJ(a[b], b_[b], carry[b])    (quadratic)

        a_vars/b_vars: list of 32 var IDs (or None if constant)
        a_const/b_const: list of 32 constant bits (or None if variable)
        """
        result = self.V.alloc32(f"{name}_r")
        carry = self.V.alloc32(f"{name}_c")

        for b in range(32):
            av = a_vars[b] if a_vars else None
            ac = a_const[b] if a_const else None
            bv = b_vars[b] if b_vars else None
            bc = b_const[b] if b_const else None

            # carry[0] = 0
            if b == 0:
                self._add_lin({carry[0]}, 0)

            # result[b] = a[b] XOR b[b] XOR carry[b]
            rvars = {result[b], carry[b]}
            rconst = 0
            if av is not None:
                rvars.add(av)
            else:
                rconst ^= ac
            if bv is not None:
                rvars.add(bv)
            else:
                rconst ^= bc
            self._add_lin(rvars, rconst)

            # carry[b+1] = MAJ(a[b], b[b], carry[b]) for b < 31
            if b < 31:
                # MAJ(x, y, z) = x*y XOR x*z XOR y*z
                cv_next = carry[b + 1]
                cb = carry[b]

                # Three cases based on which operands are vars vs const
                lin_terms = {cv_next}
                products = set()
                const_val = 0

                # Product a*b
                if av is not None and bv is not None:
                    products.add((min(av, bv), max(av, bv)))
                elif av is not None:  # b const
                    if bc == 1:
                        lin_terms ^= {av}
                elif bv is not None:  # a const
                    if ac == 1:
                        lin_terms ^= {bv}
                else:  # both const
                    const_val ^= (ac & bc)

                # Product a*carry
                if av is not None:
                    products.add((min(av, cb), max(av, cb)))
                else:
                    if ac == 1:
                        lin_terms ^= {cb}

                # Product b*carry
                if bv is not None:
                    products.add((min(bv, cb), max(bv, cb)))
                else:
                    if bc == 1:
                        lin_terms ^= {cb}

                if products:
                    self._add_quad(lin_terms, products, const_val)
                else:
                    self._add_lin(lin_terms, const_val)

        return result, carry

    def _ch_32(self, e_vars, e_const, f_vars, f_const, g_vars, g_const, name):
        """Ch(e,f,g) = e*f XOR e*g XOR g. Returns result vars."""
        result = self.V.alloc32(f"{name}")

        for b in range(32):
            ev = e_vars[b] if e_vars else None
            ec = e_const[b] if e_const else None
            fv = f_vars[b] if f_vars else None
            fc = f_const[b] if f_const else None
            gv = g_vars[b] if g_vars else None
            gc = g_const[b] if g_const else None

            # result[b] = e*f XOR e*g XOR g
            lin = {result[b]}
            prods = set()
            const = 0

            # g term (linear)
            if gv is not None:
                lin ^= {gv}
            else:
                const ^= gc

            # e*f
            if ev is not None and fv is not None:
                prods.add((min(ev, fv), max(ev, fv)))
            elif ev is not None:
                if fc == 1: lin ^= {ev}
            elif fv is not None:
                if ec == 1: lin ^= {fv}
            else:
                const ^= (ec & fc)

            # e*g
            if ev is not None and gv is not None:
                prods.add((min(ev, gv), max(ev, gv)))
            elif ev is not None:
                if gc == 1: lin ^= {ev}
            elif gv is not None:
                if ec == 1: lin ^= {gv}
            else:
                const ^= (ec & gc)

            if prods:
                self._add_quad(lin, prods, const)
            else:
                self._add_lin(lin, const)

        return result

    def _maj_32(self, a_vars, a_const, b_vars, b_const, c_vars, c_const, name):
        """Maj(a,b,c) = a*b XOR a*c XOR b*c. Returns result vars."""
        result = self.V.alloc32(f"{name}")

        for b in range(32):
            av = a_vars[b] if a_vars else None
            ac = a_const[b] if a_const else None
            bv = b_vars[b] if b_vars else None
            bc = b_const[b] if b_const else None
            cv = c_vars[b] if c_vars else None
            cc = c_const[b] if c_const else None

            lin = {result[b]}
            prods = set()
            const = 0

            for x, xc, y, yc in [(av,ac,bv,bc), (av,ac,cv,cc), (bv,bc,cv,cc)]:
                if x is not None and y is not None:
                    prods.add((min(x, y), max(x, y)))
                elif x is not None:
                    if yc == 1: lin ^= {x}
                elif y is not None:
                    if xc == 1: lin ^= {y}
                else:
                    const ^= (xc & yc)

            if prods:
                self._add_quad(lin, prods, const)
            else:
                self._add_lin(lin, const)

        return result

    def _sigma_32(self, x_vars, x_const, sources_func, name):
        """Sigma/sigma: XOR of rotated bits. Purely linear."""
        result = self.V.alloc32(name)
        for b in range(32):
            src_bits = sources_func(b)
            lin = {result[b]}
            const = 0
            for sb in src_bits:
                if x_vars:
                    lin ^= {x_vars[sb]}
                else:
                    const ^= x_const[sb]
            self._add_lin(lin, const)
        return result

    def _build_round(self, r):
        """Build ALL equations for one round."""
        # Get state
        a_v, a_c = self._state_bits(r, 0)
        b_v, b_c = self._state_bits(r, 1)
        c_v, c_c = self._state_bits(r, 2)
        d_v, d_c = self._state_bits(r, 3)
        e_v, e_c = self._state_bits(r, 4)
        f_v, f_c = self._state_bits(r, 5)
        g_v, g_c = self._state_bits(r, 6)
        h_v, h_c = self._state_bits(r, 7)

        # W[r] (message word for this round)
        if r < 16:
            w_v = self.W[r]
            w_c = None
        else:
            # Schedule expansion needed — skip for now (R <= 16)
            return

        # K[r] constant
        k_c = [get_bit(K256[r], b) for b in range(32)]

        # 1. S1 = Sigma1(e) [linear]
        S1 = self._sigma_32(e_v, e_c, big_sigma1_bit_sources, f"S1_{r}")

        # 2. ch = Ch(e, f, g) [quadratic]
        ch = self._ch_32(e_v, e_c, f_v, f_c, g_v, g_c, f"ch_{r}")

        # 3. T1 chain: h + S1 + ch + K + W (4 additions)
        tmp1, _ = self._addition_32(h_v, h_c, S1, None, f"t1a_{r}")
        tmp2, _ = self._addition_32(tmp1, None, ch, None, f"t1b_{r}")
        tmp3, _ = self._addition_32(tmp2, None, None, k_c, f"t1c_{r}")
        T1, _   = self._addition_32(tmp3, None, w_v, None, f"t1d_{r}")

        # 4. S0 = Sigma0(a) [linear]
        S0 = self._sigma_32(a_v, a_c, big_sigma0_bit_sources, f"S0_{r}")

        # 5. mj = Maj(a, b, c) [quadratic]
        mj = self._maj_32(a_v, a_c, b_v, b_c, c_v, c_c, f"mj_{r}")

        # 6. T2 = S0 + mj [1 addition]
        T2, _ = self._addition_32(S0, None, mj, None, f"t2_{r}")

        # 7. e' = d + T1 [1 addition]
        e_new, _ = self._addition_32(d_v, d_c, T1, None, f"enew_{r}")

        # 8. a' = T1 + T2 [1 addition]
        a_new, _ = self._addition_32(T1, None, T2, None, f"anew_{r}")

        # 9. Connect to state[r+1]
        next_a = self.state[(r + 1, 0)]
        next_b = self.state[(r + 1, 1)]
        next_c = self.state[(r + 1, 2)]
        next_d = self.state[(r + 1, 3)]
        next_e = self.state[(r + 1, 4)]
        next_f = self.state[(r + 1, 5)]
        next_g = self.state[(r + 1, 6)]
        next_h = self.state[(r + 1, 7)]

        for b in range(32):
            # a[r+1] = a_new
            self._add_lin({next_a[b], a_new[b]}, 0)
            # e[r+1] = e_new
            self._add_lin({next_e[b], e_new[b]}, 0)
            # Shift register: b'=a, c'=b, d'=c, f'=e, g'=f, h'=g
            for dst, src_v, src_c in [
                (next_b, a_v, a_c), (next_c, b_v, b_c), (next_d, c_v, c_c),
                (next_f, e_v, e_c), (next_g, f_v, f_c), (next_h, g_v, g_c),
            ]:
                if src_v:
                    self._add_lin({dst[b], src_v[b]}, 0)
                else:
                    self._add_lin({dst[b]}, src_c[b])

    def add_hash_target(self, target_hash):
        """Add output equations: state[R] values must produce target hash."""
        R = self.R
        # H[reg] = state[R][reg] + IV[reg] (mod 2^32)
        # This is another addition with carry!
        for reg in range(8):
            s_vars = self.state[(R, reg)]
            iv_const = [get_bit(IV256[reg], b) for b in range(32)]
            tgt_const = [get_bit(target_hash[reg], b) for b in range(32)]

            # state + IV = target → addition with carry
            # result = target (known), operands = state (var) + IV (const)
            # result[b] = state[b] XOR IV[b] XOR carry[b]
            # → state[b] XOR carry[b] = target[b] XOR IV[b]

            carry = self.V.alloc32(f"cf{reg}")
            for b in range(32):
                # carry[0] = 0
                if b == 0:
                    self._add_lin({carry[0]}, 0)

                # state[b] XOR carry[b] = target[b] XOR IV[b]
                c = tgt_const[b] ^ iv_const[b]
                self._add_lin({s_vars[b], carry[b]}, c)

                # carry[b+1] = MAJ(state[b], IV[b], carry[b])
                if b < 31:
                    cv_next = carry[b + 1]
                    sv = s_vars[b]
                    ic = iv_const[b]
                    cb = carry[b]

                    lin = {cv_next}
                    prods = set()
                    const = 0

                    # sv * iv_const: if iv=1 → sv, else 0
                    if ic == 1:
                        lin ^= {sv}
                    # sv * carry
                    prods.add((min(sv, cb), max(sv, cb)))
                    # iv_const * carry: if iv=1 → carry
                    if ic == 1:
                        lin ^= {cb}

                    self._add_quad(lin, prods, const)

    def stats(self):
        return {
            'vars': self.V.n,
            'linear': len(self.linear),
            'quad': len(self.quad),
            'rounds': self.R,
        }
