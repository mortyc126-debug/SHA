#!/usr/bin/env python3
"""
SHA-256 Q∩T Hybrid Solver v4.

Key insight from methodology (section 212-216):
  SHA-256 = Q (quadratic GF(2), degree 2) ∩ T (threshold/carry)
  
  v3 SAT solver treats EVERYTHING as T (CNF clauses).
  v4 separates Q and T:
    - XOR/rotation/Sigma = encoded as native XOR clauses (CryptoMiniSat)
    - Ch/Maj = quadratic GF(2), encoded as XOR of ANDs
    - Carry = threshold, encoded as standard CNF
    
  CryptoMiniSat does Gaussian elimination on XOR clauses internally.
  This lets the solver USE the algebraic (Q) structure, not just the
  threshold (T) structure.
"""

from pycryptosat import Solver as CMS
from pysat.solvers import Glucose4
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

# ============================================================
# Q∩T SAT builder — separates XOR (Q) from AND/carry (T)
# ============================================================

class QT_SAT:
    """Hybrid Q∩T encoder.
    
    Q-part (XOR clauses): rotations, Sigma, sigma, XOR operations
    T-part (CNF clauses): carry chains, Ch AND gates, Maj AND gates
    
    CryptoMiniSat handles both, with Gaussian elimination on XOR.
    """
    def __init__(self):
        self.nv = 0  # CMS is 0-indexed internally but 1-indexed in API
        self.cnf_clauses = []
        self.xor_clauses = []  # (lits, rhs) where XOR of lits = rhs
        
    def var(self):
        self.nv += 1
        return self.nv
    
    def word(self, n=32):
        return [self.var() for _ in range(n)]
    
    def fix(self, w, val):
        """Fix word to constant value."""
        for i in range(len(w)):
            if (val >> i) & 1:
                self.cnf_clauses.append([w[i]])
            else:
                self.cnf_clauses.append([-w[i]])
    
    def const(self, val, n=32):
        w = self.word(n)
        self.fix(w, val)
        return w
    
    # ---- Q-part: XOR operations (native XOR clauses) ----
    
    def xor2(self, a, b):
        """c = a XOR b. Encoded as XOR clause."""
        c = self.var()
        # a XOR b XOR c = 0  (i.e., c = a XOR b)
        self.xor_clauses.append(([a, b, c], False))
        return c
    
    def xor3(self, a, b, c):
        """d = a XOR b XOR c. Encoded as single XOR clause."""
        d = self.var()
        # a XOR b XOR c XOR d = 0
        self.xor_clauses.append(([a, b, c, d], False))
        return d
    
    def xor2_cnf(self, a, b):
        """c = a XOR b via CNF (for comparison)."""
        c = self.var()
        self.cnf_clauses += [[-a,-b,-c],[-a,b,c],[a,-b,c],[a,b,-c]]
        return c
        
    # ---- T-part: AND operations (CNF clauses) ----
    
    def and2(self, a, b):
        """c = a AND b."""
        c = self.var()
        self.cnf_clauses += [[-a,-b,c],[a,-c],[b,-c]]
        return c
    
    def not1(self, a):
        """c = NOT a. Via XOR: a XOR c = 1."""
        c = self.var()
        self.xor_clauses.append(([a, c], True))  # a XOR c = 1
        return c
    
    # ---- Hybrid: Ch and Maj (Q: quadratic = AND + XOR) ----
    
    def ch(self, e, f, g):
        """Ch(e,f,g) = (e AND f) XOR (NOT_e AND g).
        Degree 2 over GF(2). Encoded as AND (T) + XOR (Q)."""
        ef = self.and2(e, f)
        ne = self.not1(e)
        neg = self.and2(ne, g)
        return self.xor2(ef, neg)
    
    def maj(self, a, b, c):
        """Maj(a,b,c) = (a AND b) XOR (a AND c) XOR (b AND c).
        Degree 2 over GF(2). Three ANDs (T) + two XORs (Q)."""
        ab = self.and2(a, b)
        ac = self.and2(a, c)
        bc = self.and2(b, c)
        # ab XOR ac XOR bc
        return self.xor3(ab, ac, bc)
    
    # ---- Word-level operations ----
    
    def wxor3(self, a, b, c):
        """Word XOR3 — all Q (native XOR)."""
        return [self.xor3(a[i], b[i], c[i]) for i in range(32)]
    
    def wrotr(self, a, n):
        """Rotation — pure rewiring, no clauses."""
        return [a[(i+n)%32] for i in range(32)]
    
    def wshr(self, a, n):
        """Shift right — rewiring + zero constants."""
        r = []
        for i in range(32):
            if i+n < 32:
                r.append(a[i+n])
            else:
                z = self.var()
                self.cnf_clauses.append([-z])  # z = 0
                r.append(z)
        return r
    
    def wadd(self, a, b):
        """Word addition mod 2^32.
        Q-part: sum[k] = a[k] XOR b[k] XOR carry[k-1]
        T-part: carry[k] = MAJ(a[k], b[k], carry[k-1])
        """
        r = []
        cy = self.var()
        self.cnf_clauses.append([-cy])  # carry[-1] = 0
        for i in range(32):
            # Q: sum bit = XOR of three values
            r.append(self.xor3(a[i], b[i], cy))
            # T: carry = MAJ (threshold function)
            if i < 31:
                cy = self.maj(a[i], b[i], cy)
        return r
    
    def waddm(self, ws):
        r = ws[0]
        for w in ws[1:]:
            r = self.wadd(r, w)
        return r
    
    def sig0(self, x):
        return self.wxor3(self.wrotr(x,2), self.wrotr(x,13), self.wrotr(x,22))
    
    def sig1(self, x):
        return self.wxor3(self.wrotr(x,6), self.wrotr(x,11), self.wrotr(x,25))
    
    def ss0(self, x):
        return self.wxor3(self.wrotr(x,7), self.wrotr(x,18), self.wshr(x,3))
    
    def ss1(self, x):
        return self.wxor3(self.wrotr(x,17), self.wrotr(x,19), self.wshr(x,10))
    
    def wch(self, e, f, g):
        return [self.ch(e[i], f[i], g[i]) for i in range(32)]
    
    def wmaj(self, a, b, c):
        return [self.maj(a[i], b[i], c[i]) for i in range(32)]
    
    def weq(self, a, b):
        """a == b, enforced via XOR: a[i] XOR b[i] = 0."""
        for i in range(32):
            self.xor_clauses.append(([a[i], b[i]], False))
    
    def constrain_ascii_byte(self, byte_vars):
        """ASCII printable: bit7=0, at least one of bit5/bit6 set."""
        self.cnf_clauses.append([-byte_vars[7]])
        self.cnf_clauses.append([byte_vars[5], byte_vars[6]])
    
    def extract(self, model, w):
        v = 0
        for i in range(len(w)):
            if model[w[i]-1]:  # CMS returns True/False
                v |= (1 << i)
        return v


def backward_chain_from_H(target_H):
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


# ============================================================
# Solver function
# ============================================================

def solve_qt(target_hash_hex, msg_len, charset='ascii', n_rounds=64,
             use_xor=True, use_backward=True, timeout_sec=180):
    """Q∩T hybrid solver.
    
    use_xor=True: CryptoMiniSat with native XOR (Q∩T)
    use_xor=False: Glucose4 with all-CNF encoding (pure T, like v2)
    """
    target_H = [int(target_hash_hex[i*8:(i+1)*8], 16) for i in range(8)]
    
    t0 = time.time()
    s = QT_SAT()
    
    # Backward chain
    ak = {}
    if use_backward and n_rounds == 64:
        ak, _ = backward_chain_from_H(target_H)
    
    # Message bytes
    msg_bytes_vars = []
    for i in range(msg_len):
        byte_v = [s.var() for _ in range(8)]
        msg_bytes_vars.append(byte_v)
        if charset == 'ascii':
            s.constrain_ascii_byte(byte_v)
    
    # Build W[0..15]
    W_vars = []
    for wi in range(16):
        w = []
        for byte_pos in range(4):
            global_byte = wi * 4 + byte_pos
            if global_byte < msg_len:
                w = msg_bytes_vars[global_byte] + w
            elif global_byte == msg_len:
                w = list(s.const(0x80, 8)) + w
            elif wi == 15:
                len_bits = msg_len * 8
                len_byte_val = (len_bits >> (8 * (3 - byte_pos))) & 0xFF
                w = list(s.const(len_byte_val, 8)) + w
            else:
                w = list(s.const(0, 8)) + w
        W_vars.append(w)
    
    # Schedule
    for i in range(16, n_rounds):
        W_vars.append(s.waddm([s.ss1(W_vars[i-2]), W_vars[i-7],
                                s.ss0(W_vars[i-15]), W_vars[i-16]]))
    
    # Forward rounds
    sa=s.const(IV[0]); sb=s.const(IV[1]); sc=s.const(IV[2]); sd=s.const(IV[3])
    se=s.const(IV[4]); sf=s.const(IV[5]); sg=s.const(IV[6]); sh=s.const(IV[7])
    
    a_vars = {0: sa}; e_vars = {0: se}
    
    for r in range(n_rounds):
        kw = s.const(K[r])
        T1 = s.waddm([sh, s.sig1(se), s.wch(se,sf,sg), kw, W_vars[r]])
        T2 = s.waddm([s.sig0(sa), s.wmaj(sa,sb,sc)])
        na = s.wadd(T1, T2)
        ne = s.wadd(sd, T1)
        sh,sg,sf = sg,sf,se; se = ne
        sd,sc,sb = sc,sb,sa; sa = na
        a_vars[r+1] = na; e_vars[r+1] = ne
    
    # Fix output hash
    for i, (sv, iv_val) in enumerate(zip([sa,sb,sc,sd,se,sf,sg,sh], IV)):
        h_word = s.wadd(sv, s.const(iv_val))
        s.fix(h_word, target_H[i])
    
    # Backward pinning
    n_pinned = 0
    if use_backward and n_rounds == 64:
        for r in sorted(ak.keys()):
            if r in a_vars and 57 <= r <= 64:
                s.fix(a_vars[r], ak[r])
                n_pinned += 1
    
    build_time = time.time() - t0
    n_xor = len(s.xor_clauses)
    n_cnf = len(s.cnf_clauses)
    
    # Solve
    t0 = time.time()
    
    if use_xor:
        # CryptoMiniSat: native XOR + CNF
        solver = CMS(threads=1)
        for c in s.cnf_clauses:
            solver.add_clause(c)
        for lits, rhs in s.xor_clauses:
            solver.add_xor_clause(lits, rhs)
        
        sat, model = solver.solve()
    else:
        # Glucose4: convert XOR to CNF
        # Each XOR of n vars = 2^(n-1) clauses
        all_clauses = list(s.cnf_clauses)
        for lits, rhs in s.xor_clauses:
            if len(lits) == 2:
                a, b = lits
                if rhs:  # a XOR b = 1
                    all_clauses += [[a,b],[-a,-b]]
                else:    # a XOR b = 0
                    all_clauses += [[-a,b],[a,-b]]
            elif len(lits) == 3:
                a, b, c = lits
                if rhs:  # a XOR b XOR c = 1
                    all_clauses += [[a,b,c],[-a,-b,c],[-a,b,-c],[a,-b,-c]]
                else:    # a XOR b XOR c = 0
                    all_clauses += [[-a,-b,-c],[-a,b,c],[a,-b,c],[a,b,-c]]
            elif len(lits) == 4:
                a, b, c, d = lits
                # Generate 2^(n-1)=8 clauses forbidding wrong-parity assignments
                xor_cnf = []
                for mask in range(16):
                    bits = [(mask>>i)&1 for i in range(4)]
                    parity = sum(bits) % 2
                    target_parity = 1 if rhs else 0
                    if parity != target_parity:
                        # Forbid this assignment
                        clause = []
                        for i, (lit, bit) in enumerate(zip([a,b,c,d], bits)):
                            clause.append(-lit if bit else lit)
                        xor_cnf.append(clause)
                all_clauses += xor_cnf
            else:
                raise ValueError(f"XOR clause with {len(lits)} literals not supported in CNF mode")
        
        solver = Glucose4()
        for c in all_clauses:
            solver.add_clause(c)
        
        sat = solver.solve()
        if sat:
            raw_model = solver.get_model()
            model = [None] + [lit > 0 for lit in raw_model]  # 1-indexed
        solver.delete()
    
    solve_time = time.time() - t0
    
    result_msg = None
    if sat:
        if use_xor:
            # CMS model is 0-indexed list of True/False, but our vars are 1-indexed
            result_bytes = bytearray()
            for byte_vars in msg_bytes_vars:
                val = 0
                for i in range(8):
                    if model[byte_vars[i]]:  # CMS: model[var] directly (0-indexed internally, but var is 1-based and model is 0-based shifted)
                        val |= (1 << i)
                result_bytes.append(val)
            result_msg = bytes(result_bytes)
        else:
            result_bytes = bytearray()
            for byte_vars in msg_bytes_vars:
                val = 0
                for i in range(8):
                    if model[byte_vars[i]]:
                        val |= (1 << i)
                result_bytes.append(val)
            result_msg = bytes(result_bytes)
    
    return result_msg, solve_time, s.nv, n_cnf, n_xor, n_pinned, build_time


# ============================================================
# Quick test: verify CMS XOR clause interface
# ============================================================

print("=== Quick CMS XOR test ===")
s = CMS(threads=1)
# x1 XOR x2 = 1 (they must differ)
s.add_xor_clause([1, 2], True)
# x1 = 1
s.add_clause([1])
sat, model = s.solve()
print(f"  x1 XOR x2 = 1, x1=1: sat={sat}, x1={model[1]}, x2={model[2]}")
assert sat and model[1] == True and model[2] == False
print("  OK")

# ============================================================
# Benchmark: v4 Q∩T (CMS+XOR) vs v2-style (Glucose4, all CNF)
# ============================================================

print("\n" + "="*60)
print("SHA-256 Q∩T HYBRID SOLVER v4 — BENCHMARK")
print("="*60)

target = hashlib.sha256(b"Hi").hexdigest()

print(f"\nTarget: SHA-256('Hi') = {target[:16]}...")

# Test 1: CMS with XOR (Q∩T)
print("\n--- v4 Q∩T (CryptoMiniSat + native XOR) ---")
r, t, nv, nc, nx, np, bt = solve_qt(target, 2, 'ascii', 64, use_xor=True, use_backward=True)
if r:
    check = hashlib.sha256(r).hexdigest() == target
    try: msg_s = r.decode()
    except: msg_s = r.hex()
    print(f"  FOUND: '{msg_s}' hash={'ok' if check else 'FAIL'}")
print(f"  Time: {t:.2f}s (build: {bt:.2f}s)")
print(f"  {nv} vars, {nc} CNF clauses, {nx} XOR clauses, {np} pinned")
t_qt = t

# Test 2: CMS without XOR (all CNF, but still CMS solver)  
print("\n--- v4 CNF-only (CryptoMiniSat, XOR→CNF converted) ---")
r2, t2, nv2, nc2, nx2, np2, bt2 = solve_qt(target, 2, 'ascii', 64, use_xor=False, use_backward=True)
if r2:
    check = hashlib.sha256(r2).hexdigest() == target
    try: msg_s = r2.decode()
    except: msg_s = r2.hex()
    print(f"  FOUND: '{msg_s}' hash={'ok' if check else 'FAIL'}")
print(f"  Time: {t2:.2f}s (build: {bt2:.2f}s)")
# Can't easily count total CNF after conversion without modifying, skip

# Test 3: No backward pinning
print("\n--- v4 Q∩T without backward pinning ---")
r3, t3, *_ = solve_qt(target, 2, 'ascii', 64, use_xor=True, use_backward=False)
if r3:
    check = hashlib.sha256(r3).hexdigest() == target
    try: msg_s = r3.decode()
    except: msg_s = r3.hex()
    print(f"  FOUND: '{msg_s}' hash={'ok' if check else 'FAIL'}")
print(f"  Time: {t3:.2f}s")

print(f"\n=== SUMMARY for 'Hi' (2B, 64R) ===")
print(f"  v2 baseline (Glucose4):     ~96s")
print(f"  v3 backward+stvorochne:     ~84s (1.14x)")
print(f"  v4 Q∩T (CMS+XOR+bkwd):     {t_qt:.1f}s ({96/t_qt:.2f}x)")
print(f"  v4 CNF (CMS, no XOR):      {t2:.1f}s ({96/t2:.2f}x)")
print(f"  v4 Q∩T (no bkwd):          {t3:.1f}s ({96/t3:.2f}x)")

# Test 4: Round scaling
print(f"\n=== Round scaling — 'Hi' (2B) ===")
print(f"  {'NR':>3s} {'Q∩T':>8s} {'CNF':>8s} {'ratio':>6s}")
print(f"  {'─'*3} {'─'*8} {'─'*8} {'─'*6}")

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
    fake_H = [add32(x,iv) for x,iv in zip([a,b,c,d,e,f,g,h], IV)]
    return ''.join(f'{w:08x}' for w in fake_H)

for NR in [8, 16, 20, 24, 28, 32, 48, 64]:
    tgt = get_reduced_hash(b"Hi", NR) if NR < 64 else hashlib.sha256(b"Hi").hexdigest()
    bk = (NR == 64)
    
    rq, tq, *_ = solve_qt(tgt, 2, 'ascii', NR, use_xor=True, use_backward=bk, timeout_sec=180)
    rc, tc, *_ = solve_qt(tgt, 2, 'ascii', NR, use_xor=False, use_backward=bk, timeout_sec=180)
    
    sp = f"{tc/tq:.2f}x" if tq > 0.001 else "inf"
    rq_s = "ok" if rq else "FAIL"
    print(f"  {NR:3d} {tq:8.2f} {tc:8.2f} {sp:>6s} [{rq_s}]")

# Test 5: 1-byte message (very easy — should be instant)
print(f"\n=== 1-byte 'A' (8 bits, 64R) ===")
tgt = hashlib.sha256(b"A").hexdigest()
rq, tq, *_ = solve_qt(tgt, 1, 'ascii', 64, use_xor=True, use_backward=True)
rc, tc, *_ = solve_qt(tgt, 1, 'ascii', 64, use_xor=False, use_backward=True)
print(f"  Q∩T: {tq:.2f}s, CNF: {tc:.2f}s, ratio: {tc/tq:.2f}x")

print("\n=== DONE ===")
