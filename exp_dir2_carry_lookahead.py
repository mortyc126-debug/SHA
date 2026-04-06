#!/usr/bin/env python3
"""
Direction 2: Carry-lookahead SAT encoding.

Current: ripple carry — c[k] = MAJ(a[k], b[k], c[k-1])
  Depth = 32 (sequential chain, CDCL has O(n) propagation depth)

New: parallel prefix carry (Brent-Kung style)
  g[k] = a[k] AND b[k]  (generate)
  p[k] = a[k] XOR b[k]  (propagate)
  c[k] = g[k] OR (p[k] AND c[k-1])

  Combine: (G,P) for ranges using prefix computation
  (G_ij, P_ij) = (G_ik | (P_ik & G_{k+1,j}), P_ik & P_{k+1,j})

  Depth = O(log n) = 5 for 32 bits

Hypothesis: shorter carry chains help CDCL by enabling faster unit propagation.
"""

from pycryptosat import Solver as CMS
import struct, time, hashlib

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
def fMaj(a,b,c): return (a&b)^(a&c)^(b&c)
def add32(*a):
    s=0
    for x in a: s=(s+x)%MOD
    return s

class QT_CLA:
    """Q∩T with carry-lookahead adder."""
    def __init__(self, adder='cla'):
        self.nv = 0; self.cnf = []; self.xor = []
        self.adder = adder  # 'ripple' or 'cla'

    def var(self):
        self.nv += 1; return self.nv
    def word(self, n=32): return [self.var() for _ in range(n)]
    def fix(self, w, val):
        for i in range(len(w)): self.cnf.append([w[i] if (val>>i)&1 else -w[i]])
    def const(self, val, n=32):
        w = self.word(n); self.fix(w, val); return w
    def xor2(self, a, b):
        c = self.var(); self.xor.append(([a,b,c], False)); return c
    def xor3(self, a, b, c):
        d = self.var(); self.xor.append(([a,b,c,d], False)); return d
    def and2(self, a, b):
        c = self.var(); self.cnf += [[-a,-b,c],[a,-c],[b,-c]]; return c
    def or2(self, a, b):
        c = self.var(); self.cnf += [[a,b,-c],[-a,c],[-b,c]]; return c
    def not1(self, a):
        c = self.var(); self.xor.append(([a,c], True)); return c

    def ch(self, e, f, g):
        ef = self.and2(e, f); eg = self.and2(e, g)
        return self.xor3(ef, eg, g)
    def maj(self, a, b, c):
        ab = self.and2(a,b); ac = self.and2(a,c); bc = self.and2(b,c)
        return self.xor3(ab, ac, bc)

    def _combine_gp(self, g_hi, p_hi, g_lo, p_lo):
        """Combine (G,P) pairs: G_new = G_hi | (P_hi & G_lo), P_new = P_hi & P_lo"""
        p_hi_and_g_lo = self.and2(p_hi, g_lo)
        g_new = self.or2(g_hi, p_hi_and_g_lo)
        p_new = self.and2(p_hi, p_lo)
        return g_new, p_new

    def wadd_cla(self, a, b):
        """32-bit addition with carry-lookahead (parallel prefix).

        Step 1: compute g[k], p[k] for each bit
        Step 2: parallel prefix to get G[0:k], P[0:k]
        Step 3: carry[k] = G[0:k], sum[k] = p[k] XOR carry[k-1]
        """
        n = 32
        # Step 1: generate and propagate for each bit
        g = [self.and2(a[i], b[i]) for i in range(n)]
        p = [self.xor2(a[i], b[i]) for i in range(n)]

        # Step 2: Brent-Kung parallel prefix
        # We need carry[k] = G[0:k] where G[0:k] means "generate from bit 0 to k"
        # Use Kogge-Stone for simplicity (each level doubles the span)

        # Build prefix G, P arrays
        # Level 0: G[i] = g[i], P[i] = p[i]
        G = list(g)
        P = list(p)

        # Kogge-Stone: log2(32) = 5 levels
        span = 1
        while span < n:
            G_new = list(G)
            P_new = list(P)
            for i in range(span, n):
                G_new[i], P_new[i] = self._combine_gp(G[i], P[i], G[i-span], P[i-span])
            G = G_new
            P = P_new
            span *= 2

        # Step 3: carry[k] = G[k] (prefix generate from bit 0)
        # carry[-1] = 0, so carry[0] = g[0], carry[1] = G[0:1], etc.
        # Actually G[k] after prefix = carry OUT of position k
        # sum[k] = p[k] XOR carry[k-1]

        result = []
        for i in range(n):
            if i == 0:
                # sum[0] = a[0] XOR b[0] XOR 0 = p[0]
                result.append(p[i])
            else:
                # sum[i] = p[i] XOR carry[i-1] = p[i] XOR G[i-1]
                result.append(self.xor2(p[i], G[i-1]))

        return result

    def wadd_ripple(self, a, b):
        """Standard ripple carry (for comparison)."""
        r = []; cy = self.var(); self.cnf.append([-cy])
        for i in range(32):
            r.append(self.xor3(a[i], b[i], cy))
            if i < 31: cy = self.maj(a[i], b[i], cy)
        return r

    def wadd(self, a, b):
        if self.adder == 'cla':
            return self.wadd_cla(a, b)
        else:
            return self.wadd_ripple(a, b)

    def waddm(self, ws):
        r = ws[0]
        for w in ws[1:]: r = self.wadd(r, w)
        return r

    def wxor3(self, a, b, c): return [self.xor3(a[i],b[i],c[i]) for i in range(32)]
    def wrotr(self, a, n): return [a[(i+n)%32] for i in range(32)]
    def wshr(self, a, n):
        r = []
        for i in range(32):
            if i+n<32: r.append(a[i+n])
            else: z = self.var(); self.cnf.append([-z]); r.append(z)
        return r
    def sig0(self, x): return self.wxor3(self.wrotr(x,2),self.wrotr(x,13),self.wrotr(x,22))
    def sig1(self, x): return self.wxor3(self.wrotr(x,6),self.wrotr(x,11),self.wrotr(x,25))
    def ss0(self, x): return self.wxor3(self.wrotr(x,7),self.wrotr(x,18),self.wshr(x,3))
    def ss1(self, x): return self.wxor3(self.wrotr(x,17),self.wrotr(x,19),self.wshr(x,10))
    def wch(self, e, f, g): return [self.ch(e[i],f[i],g[i]) for i in range(32)]
    def wmaj(self, a, b, c): return [self.maj(a[i],b[i],c[i]) for i in range(32)]
    def weq(self, a, b):
        for i in range(32): self.xor.append(([a[i],b[i]], False))
    def constrain_ascii(self, bv):
        self.cnf.append([-bv[7]]); self.cnf.append([bv[5], bv[6]])


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


def solve_cla(target_hex, msg_len, adder='cla', n_rounds=64):
    target_H = [int(target_hex[i*8:(i+1)*8], 16) for i in range(8)]
    t0 = time.time()
    s = QT_CLA(adder=adder)

    ak, ek = backward_chain(target_H) if n_rounds == 64 else ({}, {})

    msg_vars = []
    for i in range(msg_len):
        bv = [s.var() for _ in range(8)]; msg_vars.append(bv)
        s.constrain_ascii(bv)

    W = []
    for wi in range(16):
        w = []
        for bp in range(4):
            gb = wi*4+bp
            if gb < msg_len: w = msg_vars[gb] + w
            elif gb == msg_len: w = list(s.const(0x80,8)) + w
            elif wi == 15:
                lbv = (msg_len*8>>(8*(3-bp)))&0xFF; w = list(s.const(lbv,8)) + w
            else: w = list(s.const(0,8)) + w
        W.append(w)
    for i in range(16, n_rounds):
        W.append(s.waddm([s.ss1(W[i-2]),W[i-7],s.ss0(W[i-15]),W[i-16]]))

    sa,sb,sc,sd = s.const(IV[0]),s.const(IV[1]),s.const(IV[2]),s.const(IV[3])
    se,sf,sg,sh = s.const(IV[4]),s.const(IV[5]),s.const(IV[6]),s.const(IV[7])
    a_v = {0:sa}; e_v = {0:se}

    for r in range(n_rounds):
        kw = s.const(K[r])
        T1 = s.waddm([sh, s.sig1(se), s.wch(se,sf,sg), kw, W[r]])
        T2 = s.waddm([s.sig0(sa), s.wmaj(sa,sb,sc)])
        na = s.wadd(T1, T2); ne = s.wadd(sd, T1)
        sh,sg,sf = sg,sf,se; se = ne
        sd,sc,sb = sc,sb,sa; sa = na
        a_v[r+1] = na; e_v[r+1] = ne

    for i,(sv,iv) in enumerate(zip([sa,sb,sc,sd,se,sf,sg,sh],IV)):
        hw = s.wadd(sv, s.const(iv)); s.fix(hw, target_H[i])

    if n_rounds == 64:
        for r in sorted(ak):
            if r in a_v and 57<=r<=64: s.fix(a_v[r], ak[r])
        for r in sorted(ek):
            if r in e_v and 61<=r<=64: s.fix(e_v[r], ek[r])

    build_time = time.time()-t0

    t0 = time.time()
    solver = CMS(threads=1)
    for c in s.cnf: solver.add_clause(c)
    for lits, rhs in s.xor: solver.add_xor_clause(lits, rhs)
    sat, model = solver.solve()
    solve_time = time.time()-t0

    result = None
    if sat and model:
        rb = bytearray()
        for bv in msg_vars:
            val = 0
            for i in range(8):
                if model[bv[i]]: val |= (1<<i)
            rb.append(val)
        result = bytes(rb)

    return result, solve_time, build_time, s.nv, len(s.cnf), len(s.xor)


# ============================================================
print("="*60)
print("DIR 2: Carry-Lookahead vs Ripple Carry")
print("="*60)

target = hashlib.sha256(b"Hi").hexdigest()

def get_reduced_hash(msg_bytes, n_rounds):
    m = bytearray(msg_bytes); ml = len(msg_bytes)*8; m.append(0x80)
    while len(m)%64!=56: m.append(0)
    m += struct.pack('>Q', ml)
    Wt = list(struct.unpack('>16I', m[:64]))
    for i in range(16,64):
        Wt.append(add32(((Wt[i-2]>>17|Wt[i-2]<<15)^(Wt[i-2]>>19|Wt[i-2]<<13)^(Wt[i-2]>>10))&0xFFFFFFFF,
                        Wt[i-7],
                        ((Wt[i-15]>>7|Wt[i-15]<<25)^(Wt[i-15]>>18|Wt[i-15]<<14)^(Wt[i-15]>>3))&0xFFFFFFFF,
                        Wt[i-16]))
    a,b,c,d,e,f,g,h = IV
    for r in range(n_rounds):
        T1=add32(h,((e>>6|e<<26)^(e>>11|e<<21)^(e>>25|e<<7))&0xFFFFFFFF,((e&f)^((~e)&g))&0xFFFFFFFF,K[r],Wt[r])
        T2=add32(((a>>2|a<<30)^(a>>13|a<<19)^(a>>22|a<<10))&0xFFFFFFFF,(a&b)^(a&c)^(b&c))
        h,g,f=g,f,e; e=add32(d,T1); d,c,b=c,b,a; a=add32(T1,T2)
    fake_H = [add32(x,iv) for x,iv in zip([a,b,c,d,e,f,g,h],IV)]
    return ''.join(f'{w:08x}' for w in fake_H)


print(f"\n{'Adder':>8s} {'NR':>3s} {'Time':>8s} {'Vars':>8s} {'CNF':>8s} {'XOR':>8s} {'Result'}")
print(f"{'─'*8} {'─'*3} {'─'*8} {'─'*8} {'─'*8} {'─'*8} {'─'*10}")

for NR in [8, 16, 24, 32, 48, 64]:
    tgt = get_reduced_hash(b"Hi", NR) if NR < 64 else target
    for adder in ['ripple', 'cla']:
        r, t, bt, nv, nc, nx = solve_cla(tgt, 2, adder=adder, n_rounds=NR)
        ok = ""
        if r:
            try: ok = f"'{r.decode()}'"
            except: ok = r.hex()[:8]
            if NR == 64:
                ok += " " + ("hash_ok" if hashlib.sha256(r).hexdigest() == tgt else "FAIL")
        else:
            ok = "FAIL"
        print(f"{adder:>8s} {NR:3d} {t:8.2f} {nv:8d} {nc:8d} {nx:8d} {ok}")

print("\nDone.")
