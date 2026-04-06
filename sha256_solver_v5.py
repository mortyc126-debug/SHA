#!/usr/bin/env python3
"""
SHA-256 Q∩T Solver v5 — Best of all + Constant Propagation.

Combines:
1. ANF Ch encoding (Dir 1 winner: 2.1× on 64R)
2. CMS + native XOR (v4: 2× on 64R)
3. Constant propagation: words that are fully known → no SAT vars
4. Backward pinning (v3)

Key idea: "lazy words" — either a concrete int OR a list of SAT vars.
If both operands are concrete → compute directly, 0 SAT overhead.
"""

from pycryptosat import Solver as CMS
import struct, time, hashlib
from multiprocessing import Process, Queue

MOD = 2**32
IV = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
      0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]
K = [0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
     0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
     0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
     0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
     0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
     0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
     0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
     0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2]

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&0xFFFFFFFF
def fS0(x): return rotr(x,2)^rotr(x,13)^rotr(x,22)
def fS1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def fCh(e,f,g): return ((e&f)^((~e)&g))&0xFFFFFFFF
def fMaj(a,b,c): return (a&b)^(a&c)^(b&c)
def fs0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def fs1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def add32(*a):
    s=0
    for x in a: s=(s+x)%MOD
    return s

def is_const(w):
    """Check if word is a concrete constant (int)."""
    return isinstance(w, int)

# ============================================================
# Layer 1: SAT primitives (bit-level)
# ============================================================

class SATCore:
    def __init__(self):
        self.nv = 0
        self.cnf = []
        self.xor = []
        self.n_const_ops = 0  # count of folded operations

    def var(self):
        self.nv += 1; return self.nv

    def word_vars(self, n=32):
        return [self.var() for _ in range(n)]

    def fix_word(self, w, val):
        for i in range(32):
            self.cnf.append([w[i] if (val>>i)&1 else -w[i]])

    def const_word(self, val):
        w = self.word_vars(); self.fix_word(w, val); return w

    # Bit-level gates
    def xor2(self, a, b):
        c = self.var(); self.xor.append(([a,b,c], False)); return c

    def xor3(self, a, b, c):
        d = self.var(); self.xor.append(([a,b,c,d], False)); return d

    def and2(self, a, b):
        c = self.var(); self.cnf += [[-a,-b,c],[a,-c],[b,-c]]; return c

    def or2(self, a, b):
        c = self.var(); self.cnf += [[a,b,-c],[-a,c],[-b,c]]; return c


# ============================================================
# Layer 2: Word-level ops with constant propagation
# ============================================================

class WordOps:
    def __init__(self, core):
        self.c = core

    def make_const(self, val):
        """Create a concrete constant (no SAT vars)."""
        self.c.n_const_ops += 1
        return val & 0xFFFFFFFF

    def make_var_word(self, n=32):
        """Create a variable word (SAT vars)."""
        return self.c.word_vars(n)

    def fix(self, w, val):
        """Pin a variable word to a value."""
        if is_const(w):
            return  # already const, nothing to do
        self.c.fix_word(w, val)

    def emit_const(self, val):
        """Convert constant to SAT word (when needed for constraints)."""
        return self.c.const_word(val & 0xFFFFFFFF)

    # --- Rotation/shift (zero cost for both const and var) ---

    def wrotr(self, a, n):
        if is_const(a): return rotr(a, n)
        return [a[(i+n)%32] for i in range(32)]

    def wshr(self, a, n):
        if is_const(a): return (a >> n) & 0xFFFFFFFF
        r = []
        for i in range(32):
            if i+n < 32: r.append(a[i+n])
            else:
                z = self.c.var(); self.c.cnf.append([-z]); r.append(z)
        return r

    # --- XOR (free for constants) ---

    def wxor3(self, a, b, c):
        if is_const(a) and is_const(b) and is_const(c):
            self.c.n_const_ops += 1
            return (a ^ b ^ c) & 0xFFFFFFFF
        # Ensure all are var words
        if is_const(a): a = self.emit_const(a)
        if is_const(b): b = self.emit_const(b)
        if is_const(c): c = self.emit_const(c)
        return [self.c.xor3(a[i],b[i],c[i]) for i in range(32)]

    # --- Sigma functions ---

    def sig0(self, x):
        return self.wxor3(self.wrotr(x,2), self.wrotr(x,13), self.wrotr(x,22))

    def sig1(self, x):
        return self.wxor3(self.wrotr(x,6), self.wrotr(x,11), self.wrotr(x,25))

    def ss0(self, x):
        return self.wxor3(self.wrotr(x,7), self.wrotr(x,18), self.wshr(x,3))

    def ss1(self, x):
        return self.wxor3(self.wrotr(x,17), self.wrotr(x,19), self.wshr(x,10))

    # --- Ch/Maj (ANF encoding, quadratic) ---

    def ch_bit(self, e, f, g):
        """Ch = ef ⊕ eg ⊕ g (ANF, no NOT needed)."""
        ef = self.c.and2(e, f)
        eg = self.c.and2(e, g)
        return self.c.xor3(ef, eg, g)

    def maj_bit(self, a, b, c):
        """Maj = ab ⊕ ac ⊕ bc."""
        ab = self.c.and2(a, b)
        ac = self.c.and2(a, c)
        bc = self.c.and2(b, c)
        return self.c.xor3(ab, ac, bc)

    def wch(self, e, f, g):
        if is_const(e) and is_const(f) and is_const(g):
            self.c.n_const_ops += 1
            return fCh(e, f, g)
        if is_const(e): e = self.emit_const(e)
        if is_const(f): f = self.emit_const(f)
        if is_const(g): g = self.emit_const(g)
        return [self.ch_bit(e[i],f[i],g[i]) for i in range(32)]

    def wmaj(self, a, b, c):
        if is_const(a) and is_const(b) and is_const(c):
            self.c.n_const_ops += 1
            return fMaj(a, b, c)
        if is_const(a): a = self.emit_const(a)
        if is_const(b): b = self.emit_const(b)
        if is_const(c): c = self.emit_const(c)
        return [self.maj_bit(a[i],b[i],c[i]) for i in range(32)]

    # --- Addition mod 2^32 (key: constant folding) ---

    def wadd(self, a, b, use_opt_const=False):
        """Add two words. If both const → direct compute."""
        if is_const(a) and is_const(b):
            self.c.n_const_ops += 1
            return (a + b) & 0xFFFFFFFF

        # Emit constants as SAT words, use standard var+var adder
        # (bit-level wadd_const OR/XNOR structure hurts CDCL)
        if is_const(b): b = self.emit_const(b)
        if is_const(a): a = self.emit_const(a)

        return self._wadd_var(a, b)

    def _wadd_var(self, a, b):
        """Standard addition: both operands are variable SAT words."""
        r = []; cy = self.c.var(); self.c.cnf.append([-cy])
        for i in range(32):
            r.append(self.c.xor3(a[i], b[i], cy))
            if i < 31:
                cy = self.maj_bit(a[i], b[i], cy)
        return r

    def _wadd_const(self, var_w, const_val):
        """Optimized: add variable word + known constant.

        When const bit = 0: carry[k] = var[k] AND carry[k-1]  (simpler!)
        When const bit = 1: carry[k] = var[k] OR carry[k-1]   (simpler!)
        Sum bit: var[k] XOR const_bit XOR carry[k-1]
          = var[k] XOR carry[k-1]           when const_bit=0
          = NOT(var[k] XOR carry[k-1])      when const_bit=1
        """
        r = []
        cy = self.c.var()
        self.c.cnf.append([-cy])  # carry[-1] = 0

        for i in range(32):
            cbit = (const_val >> i) & 1

            if cbit == 0:
                # sum = var[i] XOR carry
                r.append(self.c.xor2(var_w[i], cy))
                # carry_out = var[i] AND carry
                if i < 31:
                    cy = self.c.and2(var_w[i], cy)
            else:
                # sum = NOT(var[i] XOR carry) = var[i] XNOR carry
                xv = self.c.xor2(var_w[i], cy)
                s = self.c.var()
                self.c.xor.append(([xv, s], True))  # s = NOT xv
                r.append(s)
                # carry_out = var[i] OR carry
                if i < 31:
                    cy = self.c.or2(var_w[i], cy)

        return r

    def waddm(self, ws):
        """Add multiple words with constant folding.
        Pre-sum all constants, then add to variable part."""
        const_sum = 0
        var_words = []
        for w in ws:
            if is_const(w):
                const_sum = (const_sum + w) & 0xFFFFFFFF
            else:
                var_words.append(w)

        if not var_words:
            self.c.n_const_ops += 1
            return const_sum  # all constants!

        # Add variable words together
        result = var_words[0]
        for vw in var_words[1:]:
            result = self._wadd_var(result, vw)

        # Add constant sum
        if const_sum != 0:
            result = self._wadd_const(result, const_sum)

        return result

    def weq(self, a, b):
        """Assert a == b via XOR clauses."""
        if is_const(a) and is_const(b):
            assert a == b, f"Constant mismatch: {a} != {b}"
            return
        if is_const(a):
            self.fix(b, a); return
        if is_const(b):
            self.fix(a, b); return
        for i in range(32):
            self.c.xor.append(([a[i], b[i]], False))

    def constrain_ascii(self, bv):
        self.c.cnf.append([-bv[7]])
        self.c.cnf.append([bv[5], bv[6]])


# ============================================================
# Layer 3: SHA-256 solver
# ============================================================

def backward_chain(target_H):
    ak = {}; ek = {}
    for i in range(4):
        ak[64-i] = (target_H[i] - IV[i]) % MOD
        ek[64-i] = (target_H[4+i] - IV[4+i]) % MOD
    for step in range(4):
        r = 63 - step
        T2 = add32(fS0(ak[r]), fMaj(ak[r], ak[r-1], ak[r-2]))
        T1 = (ak[r+1] - T2) % MOD
        ak[r-3] = (ek[r+1] - T1) % MOD
    return ak, ek


def solve_v5(target_hex, msg_len, n_rounds=64, use_backward=True, timeout=180):
    target_H = [int(target_hex[i*8:(i+1)*8], 16) for i in range(8)]

    t0 = time.time()
    core = SATCore()
    w = WordOps(core)

    # Backward chain
    ak, ek = backward_chain(target_H) if use_backward and n_rounds == 64 else ({}, {})

    # Message bytes (variable)
    msg_vars = []
    for i in range(msg_len):
        bv = [core.var() for _ in range(8)]
        msg_vars.append(bv)
        w.constrain_ascii(bv)

    # Build W[0..15] — mix of variable and constant!
    W = []
    for wi in range(16):
        word_bytes = []  # (byte_value_or_vars, is_const)
        all_const = True
        byte_vals = []

        for bp in range(4):
            gb = wi * 4 + bp
            if gb < msg_len:
                all_const = False
            elif gb == msg_len:
                byte_vals.append((bp, 0x80))
            elif wi == 15:
                lbv = (msg_len * 8 >> (8 * (3 - bp))) & 0xFF
                byte_vals.append((bp, lbv))
            else:
                byte_vals.append((bp, 0x00))

        if all_const:
            # Entire word is constant!
            val = 0
            for bp, bval in byte_vals:
                val |= (bval << (8 * (3 - bp)))
            W.append(val)  # concrete int, no SAT vars
        else:
            # Mix: some bytes are message vars, some are constant
            # Build word from bytes (big-endian)
            word_bits = []
            for bp in range(4):
                gb = wi * 4 + bp
                if gb < msg_len:
                    # Variable byte (message), stored big-endian
                    byte_bits = msg_vars[gb]  # LSB-first bits
                else:
                    # Constant byte
                    bval = 0
                    if gb == msg_len:
                        bval = 0x80
                    elif wi == 15:
                        bval = (msg_len * 8 >> (8 * (3 - bp))) & 0xFF
                    cbits = []
                    for bi in range(8):
                        v = core.var()
                        core.cnf.append([v if (bval >> bi) & 1 else -v])
                        cbits.append(v)
                    byte_bits = cbits
                word_bits = byte_bits + word_bits  # prepend (big-endian)
            W.append(word_bits)

    # Schedule W[16..n_rounds-1] with constant propagation
    for i in range(16, n_rounds):
        W.append(w.waddm([w.ss1(W[i-2]), W[i-7], w.ss0(W[i-15]), W[i-16]]))

    # Count const schedule words
    n_const_w = sum(1 for i in range(16, len(W)) if is_const(W[i]))

    # Forward rounds with constant propagation
    sa, sb, sc, sd = IV[0], IV[1], IV[2], IV[3]  # concrete constants!
    se, sf, sg, sh = IV[4], IV[5], IV[6], IV[7]

    a_v = {0: sa}; e_v = {0: se}
    n_const_rounds = 0  # rounds where T2 is fully constant

    for r in range(n_rounds):
        # T1 = h + Σ₁(e) + Ch(e,f,g) + K[r] + W[r]
        sig1_e = w.sig1(se)
        ch_efg = w.wch(se, sf, sg)
        T1 = w.waddm([sh, sig1_e, ch_efg, K[r], W[r]])

        # T2 = Σ₀(a) + Maj(a,b,c)
        sig0_a = w.sig0(sa)
        maj_abc = w.wmaj(sa, sb, sc)
        T2 = w.wadd(sig0_a, maj_abc)

        if is_const(T2):
            n_const_rounds += 1

        na = w.wadd(T1, T2)
        ne = w.wadd(sd, T1)

        sh, sg, sf = sg, sf, se; se = ne
        sd, sc, sb = sc, sb, sa; sa = na
        a_v[r+1] = na; e_v[r+1] = ne

    # Fix output hash
    for i, (sv, iv) in enumerate(zip([sa,sb,sc,sd,se,sf,sg,sh], IV)):
        h_word = w.wadd(sv, iv)  # iv is const → wadd_const
        if is_const(h_word):
            # If final state+IV is somehow const (impossible for real case), just check
            assert h_word == target_H[i], f"Hash mismatch at word {i}"
        else:
            w.fix(h_word, target_H[i])

    # Pin backward chain
    n_pinned = 0
    if use_backward and n_rounds == 64:
        for r in sorted(ak):
            if r in a_v and 57 <= r <= 64 and not is_const(a_v[r]):
                w.fix(a_v[r], ak[r]); n_pinned += 1
        for r in sorted(ek):
            if r in e_v and 61 <= r <= 64 and not is_const(e_v[r]):
                w.fix(e_v[r], ek[r]); n_pinned += 1

    build_time = time.time() - t0

    # Solve
    t0 = time.time()
    solver = CMS(threads=1)
    for c in core.cnf:
        solver.add_clause(c)
    for lits, rhs in core.xor:
        solver.add_xor_clause(lits, rhs)

    sat, model = solver.solve()
    solve_time = time.time() - t0

    result = None
    if sat and model:
        rb = bytearray()
        for bv in msg_vars:
            val = 0
            for i in range(8):
                if model[bv[i]]: val |= (1 << i)
            rb.append(val)
        result = bytes(rb)

    return {
        'msg': result, 'time': solve_time, 'build': build_time,
        'vars': core.nv, 'cnf': len(core.cnf), 'xor': len(core.xor),
        'pinned': n_pinned, 'const_ops': core.n_const_ops,
        'const_w': n_const_w, 'const_rounds': n_const_rounds,
    }


# ============================================================
# BENCHMARK
# ============================================================

def get_reduced_hash(msg_bytes, n_rounds):
    m = bytearray(msg_bytes); ml = len(msg_bytes)*8; m.append(0x80)
    while len(m)%64!=56: m.append(0)
    m += struct.pack('>Q', ml)
    Wt = list(struct.unpack('>16I', m[:64]))
    for i in range(16,64): Wt.append(add32(fs1(Wt[i-2]),Wt[i-7],fs0(Wt[i-15]),Wt[i-16]))
    a,b,c,d,e,f,g,h = IV
    for r in range(n_rounds):
        T1=add32(h,fS1(e),fCh(e,f,g),K[r],Wt[r])
        T2=add32(fS0(a),fMaj(a,b,c))
        h,g,f=g,f,e; e=add32(d,T1); d,c,b=c,b,a; a=add32(T1,T2)
    fake_H = [add32(x,iv) for x,iv in zip([a,b,c,d,e,f,g,h],IV)]
    return ''.join(f'{w:08x}' for w in fake_H)


print("="*60)
print("SHA-256 SOLVER v5 — COMBINED BEST + CONST PROPAGATION")
print("="*60)

# Test 1: Verify on "Hi"
print("\n=== Test 1: Verify correctness ===")
target = hashlib.sha256(b"Hi").hexdigest()
r = solve_v5(target, 2, 64)
if r['msg']:
    check = hashlib.sha256(r['msg']).hexdigest() == target
    try: ms = r['msg'].decode()
    except: ms = r['msg'].hex()
    print(f"  FOUND: '{ms}' hash={'ok' if check else 'FAIL'}")
print(f"  Time: {r['time']:.2f}s (build: {r['build']:.2f}s)")
print(f"  Vars: {r['vars']} | CNF: {r['cnf']} | XOR: {r['xor']}")
print(f"  Const folded: {r['const_ops']} ops, {r['const_w']} sched words, {r['const_rounds']} T2-const rounds")
print(f"  Pinned: {r['pinned']}")

# Compare with v4
print(f"\n  vs v2 baseline (~96s): {96/r['time']:.2f}×")
print(f"  vs v4 Q∩T ANF (~41s): {41/r['time']:.2f}×")

# Test 2: Round scaling
print(f"\n=== Test 2: Round scaling 'Hi' ===")
print(f"  {'NR':>3s} {'time':>8s} {'vars':>7s} {'cnf':>7s} {'xor':>7s} {'folded':>7s}")
print(f"  {'─'*3} {'─'*8} {'─'*7} {'─'*7} {'─'*7} {'─'*7}")

for NR in [8, 16, 20, 24, 28, 32, 48, 64]:
    tgt = get_reduced_hash(b"Hi", NR) if NR < 64 else target
    bk = (NR == 64)
    r = solve_v5(tgt, 2, NR, use_backward=bk)
    ok = "ok" if r['msg'] else "FAIL"
    print(f"  {NR:3d} {r['time']:8.2f} {r['vars']:7d} {r['cnf']:7d} {r['xor']:7d} {r['const_ops']:7d} [{ok}]")

# Test 3: 1-byte (easy)
print(f"\n=== Test 3: 'A' 1-byte 64R ===")
tgt = hashlib.sha256(b"A").hexdigest()
r = solve_v5(tgt, 1, 64)
print(f"  Time: {r['time']:.2f}s, Vars: {r['vars']}, Folded: {r['const_ops']}")

# Test 4: 3-byte short rounds only
print(f"\n=== Test 4: 3-byte 'Abc' short rounds ===")
for NR in [8, 16]:
    tgt = get_reduced_hash(b"Abc", NR)
    r = solve_v5(tgt, 3, NR, use_backward=False, timeout=180)
    ok = "ok" if r['msg'] else "TIMEOUT/FAIL"
    print(f"  {NR}R: {r['time']:.2f}s, {r['vars']}v, {r['const_ops']} folded [{ok}]")

print("\n=== DONE ===")
