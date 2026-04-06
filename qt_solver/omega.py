"""
OMEGA System: SHA-256 as GF(2) equations with carry as VARIABLES.

Key innovation over Q∩T: carry bits are unknowns, not constants.
This gives α-kernel = max(0, 32*(16-R)) for R-round SHA-256.

Variables:
  - Message bits: W[w][b] for w=0..15, b=0..31 (512 vars)
  - State bits: registers at each round
  - Carry bits: for every addition at every bit position

Equations:
  - Linear: XOR operations, rotations, register shifts
  - Quadratic: carry = MAJ(x, y, c_in), Ch, Maj functions
"""

from qt_solver.sha256_traced import (
    K, IV, MASK32, get_bit, rotr,
    sigma0, sigma1, ssigma0, ssigma1,
    ch, maj, sha256_compress_traced,
)


class VarAllocator:
    """Allocate variable indices for GF(2) system."""

    def __init__(self):
        self.n = 0
        self.names = {}

    def alloc(self, name):
        idx = self.n
        self.n += 1
        self.names[idx] = name
        return idx

    def alloc_word(self, prefix):
        """Allocate 32 variables for a word, return list of 32 indices."""
        return [self.alloc(f"{prefix}[{b}]") for b in range(32)]


class OmegaSystem:
    """
    Build GF(2) equation system for R-round SHA-256.
    Carries are VARIABLES (not constants).

    Equations stored as:
      linear_eqs: list of (row, rhs) where row is a set of variable indices
                  meaning XOR of those variables = rhs (0 or 1)
      quadratic_eqs: list of (var_a, var_b, var_c, val) meaning
                     var_a AND var_b XOR var_c = val
                     (represents carry = MAJ equation)
    """

    def __init__(self, R):
        self.R = R
        self.V = VarAllocator()
        self.linear_eqs = []    # (set_of_vars, rhs)
        self.quadratic_eqs = [] # for carry MAJ equations

        # Allocate message variables
        self.W = []  # W[word][bit] = var index
        for w in range(16):
            self.W.append(self.V.alloc_word(f"W{w}"))

        # Allocate schedule words W[16..R-1] as derived (for R > 16 they
        # are determined by W[0..15] via linear recurrence, so we add
        # new vars + constraints)
        self.W_extended = list(self.W)  # will grow
        self.schedule_carries = {}

        # Build schedule
        self._build_schedule()

        # Allocate state variables and build round equations
        # State: only track a and e (b=a[-1], c=a[-2], etc.)
        self.a = []  # a[round][bit]
        self.e = []  # e[round][bit]

        # IV state — these are constants, not variables
        # We represent them by adding equations: var = constant_bit
        self.a.append(self.V.alloc_word("a0"))
        self.e.append(self.V.alloc_word("e0"))

        # Fix IV values
        iv_a, iv_b, iv_c, iv_d = IV[0], IV[1], IV[2], IV[3]
        iv_e, iv_f, iv_g, iv_h = IV[4], IV[5], IV[6], IV[7]

        # a[0] = IV[0], but we also need b[0]=IV[1], c[0]=IV[2], d[0]=IV[3]
        # and e[0]=IV[4], f[0]=IV[5], g[0]=IV[6], h[0]=IV[7]
        # Instead of separate b,c,d,f,g,h variables, we track history:
        # b[r] = a[r-1], c[r] = a[r-2], d[r] = a[r-3]
        # f[r] = e[r-1], g[r] = e[r-2], h[r] = e[r-3]

        # So we need a[-3..0] and e[-3..0] for the IV
        # a[0] = IV[0], a[-1] = IV[1], a[-2] = IV[2], a[-3] = IV[3]
        # e[0] = IV[4], e[-1] = IV[5], e[-2] = IV[6], e[-3] = IV[7]

        # Store IV history as constant arrays (no variables needed)
        self.iv_a_hist = [IV[3], IV[2], IV[1], IV[0]]  # a[-3], a[-2], a[-1], a[0]
        self.iv_e_hist = [IV[7], IV[6], IV[5], IV[4]]  # e[-3], e[-2], e[-1], e[0]

        # Fix a[0] and e[0] to IV values
        for b in range(32):
            self._fix_var(self.a[0][b], get_bit(IV[0], b))
            self._fix_var(self.e[0][b], get_bit(IV[4], b))

        # Build round equations
        self.round_carries = {}
        self._build_rounds()

    def _fix_var(self, var, val):
        """Add equation: var = val (constant)."""
        self.linear_eqs.append(({var}, val))

    def _add_xor_eq(self, var_set, rhs):
        """Add equation: XOR of vars in var_set = rhs."""
        self.linear_eqs.append((set(var_set), rhs))

    def _get_a(self, r):
        """Get a[r] word (list of 32 var indices). r can be negative (IV)."""
        if r <= 0:
            # Return constant bits from IV history
            idx = r + 3  # a[-3] -> idx 0, a[0] -> idx 3
            if 0 <= idx <= 3:
                val = self.iv_a_hist[idx]
                return val  # return integer, not var list
        return self.a[r]  # return var list

    def _get_e(self, r):
        """Get e[r] word. r can be negative (IV)."""
        if r <= 0:
            idx = r + 3
            if 0 <= idx <= 3:
                val = self.iv_e_hist[idx]
                return val
        return self.e[r]

    def _is_const(self, word):
        """Check if word is a constant (int) vs variable list."""
        return isinstance(word, int)

    def _get_bit_var_or_const(self, word, bit):
        """
        Get (is_const, value_or_var) for bit position of a word.
        If word is constant int: returns (True, bit_value)
        If word is var list: returns (False, var_index)
        """
        if isinstance(word, int):
            return (True, get_bit(word, bit))
        else:
            return (False, word[bit])

    def _build_schedule(self):
        """Build message schedule W[16..] from W[0..15] using σ₀, σ₁."""
        needed = max(self.R, 16)
        for i in range(16, needed):
            # W[i] = σ₁(W[i-2]) + W[i-7] + σ₀(W[i-15]) + W[i-16]
            # σ₀ and σ₁ are GF(2)-linear (rotations + shift)
            # + is mod 2^32, needs carry variables

            # Allocate result word and carry variables for 3 additions
            wi = self.V.alloc_word(f"W{i}")
            c1 = self.V.alloc_word(f"Wc{i}_1")  # carry for add 1
            c2 = self.V.alloc_word(f"Wc{i}_2")  # carry for add 2
            c3 = self.V.alloc_word(f"Wc{i}_3")  # carry for add 3

            self.schedule_carries[i] = (c1, c2, c3)

            # Intermediate sums
            s1 = self.V.alloc_word(f"Ws{i}_1")  # σ₁(W[i-2]) + W[i-7]
            s2 = self.V.alloc_word(f"Ws{i}_2")  # s1 + σ₀(W[i-15])

            # σ₁(W[i-2]): ROTR(17) ^ ROTR(19) ^ SHR(10)
            # σ₀(W[i-15]): ROTR(7) ^ ROTR(18) ^ SHR(3)
            # These are linear in GF(2) — we express them as XOR of rotated bits

            w_im2 = self.W_extended[i-2]
            w_im7 = self.W_extended[i-7]
            w_im15 = self.W_extended[i-15]
            w_im16 = self.W_extended[i-16]

            # Build addition equations for each bit position
            for b in range(32):
                # Compute σ₁(W[i-2]) at bit b:
                # ROTR(x,17)[b] = x[(b+17)%32]
                # ROTR(x,19)[b] = x[(b+19)%32]
                # SHR(x,10)[b] = x[b+10] if b+10 < 32 else 0
                sig1_vars = set()
                sig1_vars.add(w_im2[(b+17) % 32])
                sig1_vars.add(w_im2[(b+19) % 32])
                if b + 10 < 32:
                    sig1_vars.add(w_im2[b+10])
                # If a var appears twice in XOR, it cancels
                # We handle via symmetric difference (set XOR)

                # σ₀(W[i-15]) at bit b:
                sig0_vars = set()
                sig0_vars.add(w_im15[(b+7) % 32])
                sig0_vars.add(w_im15[(b+18) % 32])
                if b + 3 < 32:
                    sig0_vars.add(w_im15[b+3])

                # Addition 1: s1 = σ₁(W[i-2]) + W[i-7]
                # s1[b] = σ₁[b] XOR W[i-7][b] XOR c1[b]
                # c1[0] = 0
                # c1[b] = MAJ(σ₁[b-1], W[i-7][b-1], c1[b-1]) for b>0
                # This is hard because σ₁ is itself a XOR of bits...
                # We need auxiliary variables for the σ values

            # Actually, the proper way: allocate σ₁ and σ₀ as
            # auxiliary word variables constrained by linear equations,
            # then do additions on them.

            sig1_word = self.V.alloc_word(f"sig1_W{i}")
            sig0_word = self.V.alloc_word(f"sig0_W{i}")

            for b in range(32):
                # σ₁(W[i-2])[b] = W[i-2][(b+17)%32] ^ W[i-2][(b+19)%32] ^ (W[i-2][b+10] if b+10<32 else 0)
                xor_set = {sig1_word[b]}  # we want sig1_word[b] XOR ... = 0
                xor_set.add(w_im2[(b+17) % 32])
                xor_set.add(w_im2[(b+19) % 32])
                if b + 10 < 32:
                    xor_set.add(w_im2[b + 10])
                # Handle duplicates by symmetric difference
                self._add_xor_eq(self._sym_diff(xor_set), 0)

                # σ₀(W[i-15])[b]
                xor_set2 = {sig0_word[b]}
                xor_set2.add(w_im15[(b+7) % 32])
                xor_set2.add(w_im15[(b+18) % 32])
                if b + 3 < 32:
                    xor_set2.add(w_im15[b + 3])
                self._add_xor_eq(self._sym_diff(xor_set2), 0)

            # Now 3 additions:
            # s1 = sig1_word + w_im7
            self._add_addition(sig1_word, w_im7, s1, c1)
            # s2 = s1 + sig0_word
            self._add_addition(s1, sig0_word, s2, c2)
            # wi = s2 + w_im16
            self._add_addition(s2, w_im16, wi, c3)

            self.W_extended.append(wi)

    def _sym_diff(self, s):
        """Compute symmetric difference for XOR (handle duplicate cancellation)."""
        # When building XOR equations, if a variable appears twice it cancels.
        # We need proper multiset → set reduction.
        # Since we use sets, duplicates are already removed.
        # But the issue is when we ADD the same var twice to a set — it stays once.
        # We need to track parity:
        counts = {}
        for v in s:
            counts[v] = counts.get(v, 0) + 1
        return {v for v, c in counts.items() if c % 2 == 1}

    def _add_addition(self, x_vars, y_vars, result_vars, carry_vars):
        """
        Add equations for 32-bit addition: result = x + y.
        x_vars, y_vars, result_vars, carry_vars: lists of 32 var indices each.

        Equations per bit b:
          result[b] = x[b] XOR y[b] XOR carry[b]       (linear)
          carry[0] = 0                                   (constant)
          carry[b+1] = MAJ(x[b], y[b], carry[b])        (quadratic)
        """
        # carry[0] = 0
        self._fix_var(carry_vars[0], 0)

        for b in range(32):
            # result[b] = x[b] ^ y[b] ^ carry[b]
            self._add_xor_eq({result_vars[b], x_vars[b], y_vars[b], carry_vars[b]}, 0)

            # carry[b+1] = MAJ(x[b], y[b], carry[b])
            # MAJ(a,b,c) = ab ^ ac ^ bc
            # carry_out = x[b]*y[b] ^ x[b]*carry[b] ^ y[b]*carry[b]
            if b < 31:
                self.quadratic_eqs.append({
                    'type': 'carry_maj',
                    'result': carry_vars[b+1],
                    'inputs': (x_vars[b], y_vars[b], carry_vars[b]),
                })

    def _add_addition_const(self, x_vars, const_val, result_vars, carry_vars):
        """
        Add equations for result = x + constant.
        Constant bits are known, so carry becomes:
          carry[b+1] = x[b]*const[b] ^ x[b]*carry[b] ^ const[b]*carry[b]
        When const[b]=0: carry[b+1] = 0 (if carry[b]=0 initially, it stays 0... no)
        Actually: carry[b+1] = MAJ(x[b], const[b], carry[b])
        With const[b] known, this simplifies.
        """
        self._fix_var(carry_vars[0], 0)

        for b in range(32):
            cb = get_bit(const_val, b)
            # result[b] = x[b] ^ const[b] ^ carry[b]
            self._add_xor_eq({result_vars[b], x_vars[b], carry_vars[b]}, cb)

            if b < 31:
                # carry[b+1] = MAJ(x[b], const_bit, carry[b])
                # If const_bit = 0: MAJ(x,0,c) = x*c (AND)
                # If const_bit = 1: MAJ(x,1,c) = x + c - x*c = x^c^(x*c)... no
                # MAJ(x,1,c) = x*1 ^ x*c ^ 1*c = x ^ xc ^ c
                # Actually MAJ(a,b,c) = ab^ac^bc
                # MAJ(x,1,c) = x^xc^c  ... hmm
                # Let's just use the general quadratic form but mark const
                self.quadratic_eqs.append({
                    'type': 'carry_maj_const',
                    'result': carry_vars[b+1],
                    'var_inputs': (x_vars[b], carry_vars[b]),
                    'const_bit': cb,
                })

    def _build_rounds(self):
        """Build round function equations for R rounds."""
        for r in range(self.R):
            # Get register values for this round
            # a[r], b[r]=a[r-1], c[r]=a[r-2], d[r]=a[r-3]
            # e[r], f[r]=e[r-1], g[r]=e[r-2], h[r]=e[r-3]
            a_r = self._get_a(r) if r > 0 else self.a[0]
            b_r = self._get_a(r-1) if r >= 1 else self._get_a(r-1)
            c_r = self._get_a(r-2) if r >= 2 else self._get_a(r-2)
            d_r = self._get_a(r-3) if r >= 3 else self._get_a(r-3)
            e_r = self._get_e(r) if r > 0 else self.e[0]
            f_r = self._get_e(r-1) if r >= 1 else self._get_e(r-1)
            g_r = self._get_e(r-2) if r >= 2 else self._get_e(r-2)
            h_r = self._get_e(r-3) if r >= 3 else self._get_e(r-3)

            # For round r:
            # T1 = h + Σ₁(e) + Ch(e,f,g) + K[r] + W[r]
            # T2 = Σ₀(a) + Maj(a,b,c)
            # a[r+1] = T1 + T2
            # e[r+1] = d + T1

            # Allocate new state variables
            a_new = self.V.alloc_word(f"a{r+1}")
            e_new = self.V.alloc_word(f"e{r+1}")
            self.a.append(a_new)
            self.e.append(e_new)

            # Allocate intermediate and carry variables
            # Σ₁(e) — linear, needs aux var
            sig1_e = self.V.alloc_word(f"Sig1_e_r{r}")
            # Σ₀(a) — linear, needs aux var
            sig0_a = self.V.alloc_word(f"Sig0_a_r{r}")
            # Ch(e,f,g) — quadratic, needs aux var
            ch_efg = self.V.alloc_word(f"Ch_r{r}")
            # Maj(a,b,c) — quadratic, needs aux var
            maj_abc = self.V.alloc_word(f"Maj_r{r}")

            # T1 intermediates
            t1_s1 = self.V.alloc_word(f"T1s1_r{r}")  # h + Σ₁(e)
            t1_s2 = self.V.alloc_word(f"T1s2_r{r}")  # + Ch
            t1_s3 = self.V.alloc_word(f"T1s3_r{r}")  # + K[r]
            t1 = self.V.alloc_word(f"T1_r{r}")        # + W[r]
            t2 = self.V.alloc_word(f"T2_r{r}")        # Σ₀(a) + Maj

            # Carry variables for each addition (7 per round)
            c0 = self.V.alloc_word(f"c_r{r}_0")  # h + Σ₁(e)
            c1 = self.V.alloc_word(f"c_r{r}_1")  # + Ch
            c2 = self.V.alloc_word(f"c_r{r}_2")  # + K
            c3 = self.V.alloc_word(f"c_r{r}_3")  # + W
            c4 = self.V.alloc_word(f"c_r{r}_4")  # Σ₀ + Maj
            c5 = self.V.alloc_word(f"c_r{r}_5")  # d + T1
            c6 = self.V.alloc_word(f"c_r{r}_6")  # T1 + T2

            self.round_carries[r] = (c0, c1, c2, c3, c4, c5, c6)

            # === Build Σ₁(e) equations (linear) ===
            # Σ₁(e)[b] = e[(b+6)%32] ^ e[(b+11)%32] ^ e[(b+25)%32]
            for b in range(32):
                eq_vars = {sig1_e[b]}
                self._add_rot_xor(eq_vars, e_r, [(b+6)%32, (b+11)%32, (b+25)%32])
                # The equation is: sig1_e[b] ^ e_rot1 ^ e_rot2 ^ e_rot3 = 0

            # === Build Σ₀(a) equations (linear) ===
            for b in range(32):
                eq_vars = {sig0_a[b]}
                self._add_rot_xor(eq_vars, a_r, [(b+2)%32, (b+13)%32, (b+22)%32])

            # === Build Ch(e,f,g) equations (quadratic) ===
            # Ch(e,f,g)[b] = e[b]*f[b] ^ (1^e[b])*g[b] = e[b]*f[b] ^ g[b] ^ e[b]*g[b]
            for b in range(32):
                self.quadratic_eqs.append({
                    'type': 'ch',
                    'result': ch_efg[b],
                    'e': self._var_or_const(e_r, b),
                    'f': self._var_or_const(f_r, b),
                    'g': self._var_or_const(g_r, b),
                })

            # === Build Maj(a,b,c) equations (quadratic) ===
            # Maj(a,b,c)[b] = a[b]*b[b] ^ a[b]*c[b] ^ b[b]*c[b]
            for b in range(32):
                self.quadratic_eqs.append({
                    'type': 'maj',
                    'result': maj_abc[b],
                    'a': self._var_or_const(a_r, b),
                    'b': self._var_or_const(b_r, b),
                    'c': self._var_or_const(c_r, b),
                })

            # === Build T1 additions ===
            if self._is_const(h_r):
                self._add_addition_const(sig1_e, h_r, t1_s1, c0)
            else:
                self._add_addition(h_r, sig1_e, t1_s1, c0)

            self._add_addition(t1_s1, ch_efg, t1_s2, c1)
            self._add_addition_const(t1_s2, K[r], t1_s3, c2)

            w_r = self.W_extended[r] if r < len(self.W_extended) else self.W_extended[r]
            self._add_addition(t1_s3, w_r, t1, c3)

            # === Build T2 ===
            self._add_addition(sig0_a, maj_abc, t2, c4)

            # === e_new = d + T1 ===
            if self._is_const(d_r):
                self._add_addition_const(t1, d_r, e_new, c5)
            else:
                self._add_addition(d_r, t1, e_new, c5)

            # === a_new = T1 + T2 ===
            self._add_addition(t1, t2, a_new, c6)

    def _var_or_const(self, word, bit):
        """Return (is_const, val_or_var) for a word's bit."""
        if isinstance(word, int):
            return ('const', get_bit(word, bit))
        return ('var', word[bit])

    def _add_rot_xor(self, eq_vars, src_word, rot_positions):
        """
        Add rotation-XOR equation.
        eq_vars starts with {result_var}.
        src_word is either int (constant) or list of var indices.
        rot_positions: list of bit indices to XOR from src_word.
        Result: result_var = XOR of src_word[pos] for pos in rot_positions.
        """
        rhs = 0
        if isinstance(src_word, int):
            # All constant
            for pos in rot_positions:
                rhs ^= get_bit(src_word, pos)
            self._add_xor_eq(eq_vars, rhs)
        else:
            for pos in rot_positions:
                if src_word[pos] in eq_vars:
                    eq_vars.remove(src_word[pos])  # XOR cancels
                else:
                    eq_vars.add(src_word[pos])
            self._add_xor_eq(eq_vars, 0)

    def get_msg_var_mask(self):
        """Return bitmask of all message variables (W[0..15])."""
        mask = 0
        for w in range(16):
            for b in range(32):
                mask |= (1 << self.W[w][b])
        return mask

    def summary(self):
        """Print system summary."""
        print(f"OMEGA System for R={self.R} rounds")
        print(f"  Variables: {self.V.n}")
        print(f"  Linear equations: {len(self.linear_eqs)}")
        print(f"  Quadratic equations: {len(self.quadratic_eqs)}")
        print(f"  Message variables: 512 (W[0..15])")
        print(f"  Expected α-kernel: {max(0, 32*(16-self.R))}")
