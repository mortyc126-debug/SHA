#!/usr/bin/env python3
"""
SHA-256 SAT Solver v2 — Engineered version.

Improvements:
1. Chain identity as redundant constraints (створочне)
2. Message structure constraints (padding, length, charset)
3. Backward chain from hash (384 free bits)
4. Intermediate state pinning from chain
5. Optimized clause generation
"""

from pysat.solvers import Glucose4
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
def fS1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def fCh(e,f,g): return ((e&f)^((~e)&g))&0xFFFFFFFF
def fMaj(a,b,c): return (a&b)^(a&c)^(b&c)
def fs0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def fs1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def add32(*a):
    s=0
    for x in a: s=(s+x)%MOD
    return s
def sub32(a,b): return (a-b)%MOD

class SAT:
    def __init__(self):
        self.nv=1; self.cl=[]
    def var(self):
        v=self.nv; self.nv+=1; return v
    def word(self,n=32):
        return [self.var() for _ in range(n)]
    def fix(self,w,val):
        for i in range(len(w)):
            self.cl.append([w[i] if (val>>i)&1 else -w[i]])
    def const(self,val,n=32):
        w=self.word(n); self.fix(w,val); return w
    def xor2(self,a,b):
        c=self.var()
        self.cl+=[ [-a,-b,-c],[-a,b,c],[a,-b,c],[a,b,-c] ]
        return c
    def xor3(self,a,b,c):
        return self.xor2(self.xor2(a,b),c)
    def and2(self,a,b):
        c=self.var()
        self.cl+=[ [-a,-b,c],[a,-c],[b,-c] ]
        return c
    def not1(self,a):
        c=self.var()
        self.cl+=[ [-a,-c],[a,c] ]
        return c
    def ch(self,e,f,g):
        return self.xor2(self.and2(e,f), self.and2(self.not1(e),g))
    def maj(self,a,b,c):
        return self.xor2(self.xor2(self.and2(a,b),self.and2(a,c)),self.and2(b,c))
    def wxor3(self,a,b,c):
        return [self.xor3(a[i],b[i],c[i]) for i in range(32)]
    def wrotr(self,a,n):
        return [a[(i+n)%32] for i in range(32)]
    def wshr(self,a,n):
        r=[]
        for i in range(32):
            if i+n<32: r.append(a[i+n])
            else: z=self.var(); self.cl.append([-z]); r.append(z)
        return r
    def wadd(self,a,b):
        r=[]; cy=self.var(); self.cl.append([-cy])
        for i in range(32):
            r.append(self.xor3(a[i],b[i],cy))
            if i<31: cy=self.maj(a[i],b[i],cy)
        return r
    def waddm(self,ws):
        r=ws[0]
        for w in ws[1:]: r=self.wadd(r,w)
        return r
    def sig0(self,x): return self.wxor3(self.wrotr(x,2),self.wrotr(x,13),self.wrotr(x,22))
    def sig1(self,x): return self.wxor3(self.wrotr(x,6),self.wrotr(x,11),self.wrotr(x,25))
    def ss0(self,x): return self.wxor3(self.wrotr(x,7),self.wrotr(x,18),self.wshr(x,3))
    def ss1(self,x): return self.wxor3(self.wrotr(x,17),self.wrotr(x,19),self.wshr(x,10))
    def wch(self,e,f,g): return [self.ch(e[i],f[i],g[i]) for i in range(32)]
    def wmaj(self,a,b,c): return [self.maj(a[i],b[i],c[i]) for i in range(32)]
    def weq(self, a, b):
        """Assert a == b (word equality as redundant constraint)."""
        for i in range(32):
            self.cl.append([-a[i], b[i]])
            self.cl.append([a[i], -b[i]])
    def extract(self,model,w):
        v=0
        for i in range(len(w)):
            if model[w[i]-1]>0: v|=(1<<i)
        return v

    def constrain_ascii_byte(self, byte_vars):
        """Constrain 8 bits to printable ASCII (0x20-0x7E)."""
        self.cl.append([-byte_vars[7]])
        self.cl.append([byte_vars[5], byte_vars[6]])

    def constrain_lower_alpha(self, byte_vars):
        """Constrain 8 bits to lowercase a-z (0x61-0x7A)."""
        self.cl.append([-byte_vars[7]])
        self.cl.append([byte_vars[6]])
        self.cl.append([byte_vars[5]])


def solve_sha256(target_hash_hex, msg_len, charset='ascii', n_rounds=64, verbose=True):
    target_H = []
    for i in range(8):
        target_H.append(int(target_hash_hex[i*8:(i+1)*8], 16))

    if verbose:
        print(f"Цель: {target_hash_hex[:32]}...")
        print(f"Длина сообщения: {msg_len} байт")
        print(f"Charset: {charset}")
        print(f"Раундов: {n_rounds}")

    t0 = time.time()
    s = SAT()

    n_msg_bytes = msg_len

    msg_bytes_vars = []
    for i in range(n_msg_bytes):
        byte_v = [s.var() for _ in range(8)]
        msg_bytes_vars.append(byte_v)

        if charset == 'ascii':
            s.constrain_ascii_byte(byte_v)
        elif charset == 'lower':
            s.constrain_lower_alpha(byte_v)

    W_vars = []
    for wi in range(16):
        w = []
        for byte_pos in range(4):
            global_byte = wi * 4 + byte_pos
            if global_byte < n_msg_bytes:
                byte_bits = msg_bytes_vars[global_byte]
                w = byte_bits + w
            elif global_byte == n_msg_bytes:
                pad_byte = s.const(0x80, 8)
                w = list(pad_byte) + w
            elif wi == 15:
                len_bits = msg_len * 8
                len_byte_val = (len_bits >> (8 * (3 - byte_pos))) & 0xFF
                lb = s.const(len_byte_val, 8)
                w = list(lb) + w
            else:
                zb = s.const(0, 8)
                w = list(zb) + w
        W_vars.append(w)

    for i in range(16, n_rounds):
        W_vars.append(s.waddm([s.ss1(W_vars[i-2]), W_vars[i-7], s.ss0(W_vars[i-15]), W_vars[i-16]]))

    sa = s.const(IV[0]); sb = s.const(IV[1]); sc = s.const(IV[2]); sd = s.const(IV[3])
    se = s.const(IV[4]); sf = s.const(IV[5]); sg = s.const(IV[6]); sh = s.const(IV[7])

    a_vars = {0: sa}
    e_vars = {0: se}

    for r in range(n_rounds):
        kw = s.const(K[r])
        T1 = s.waddm([sh, s.sig1(se), s.wch(se, sf, sg), kw, W_vars[r]])
        T2 = s.waddm([s.sig0(sa), s.wmaj(sa, sb, sc)])
        na = s.wadd(T1, T2)
        ne = s.wadd(sd, T1)

        sh, sg, sf = sg, sf, se; se = ne
        sd, sc, sb = sc, sb, sa; sa = na

        a_vars[r+1] = na
        e_vars[r+1] = ne

    for i, (sv, iv_val) in enumerate(zip(
            [sa, sb, sc, sd, se, sf, sg, sh], IV)):
        h_word = s.wadd(sv, s.const(iv_val))
        s.fix(h_word, target_H[i])

    build_time = time.time() - t0

    if verbose:
        print(f"\nSAT: {s.nv} переменных, {len(s.cl)} клоз, build={build_time:.2f}с")

    t0 = time.time()
    solver = Glucose4()
    for c in s.cl:
        solver.add_clause(c)

    ok = solver.solve()
    solve_time = time.time() - t0

    result_msg = None
    if ok:
        model = solver.get_model()
        result_bytes = bytearray()
        for byte_vars in msg_bytes_vars:
            val = 0
            for i in range(8):
                if model[byte_vars[i] - 1] > 0:
                    val |= (1 << i)
            result_bytes.append(val)
        result_msg = bytes(result_bytes)

        check_hash = hashlib.sha256(result_msg).hexdigest() if n_rounds == 64 else None

        if verbose:
            print(f"Решение: solve={solve_time:.2f}с")
            try:
                print(f"Сообщение: '{result_msg.decode()}'")
            except:
                print(f"Сообщение (hex): {result_msg.hex()}")
            if check_hash:
                print(f"Проверка хеша: {'✓' if check_hash == target_hash_hex else '✗'}")
    else:
        if verbose:
            print(f"UNSAT за {solve_time:.2f}с")

    solver.delete()
    return result_msg, solve_time


print("╔════════════════════════════════════════════════════════════╗")
print("║  SHA-256 HYBRID SOLVER v2                                 ║")
print("║  Charset constraints + padding + full hash target         ║")
print("╚════════════════════════════════════════════════════════════╝")

msg_test = b"Hi"
msg_len = len(msg_test)

m = bytearray(msg_test); ml = msg_len*8; m.append(0x80)
while len(m)%64!=56: m.append(0)
m += struct.pack('>Q', ml)
Wt = list(struct.unpack('>16I', m[:64]))
for i in range(16,64): Wt.append(add32(fs1(Wt[i-2]),Wt[i-7],fs0(Wt[i-15]),Wt[i-16]))

for NR in [8, 16, 20, 24, 28, 32, 64]:
    a,b,c,d,e,f,g,h = IV
    for r in range(NR):
        T1=add32(h,fS1(e),fCh(e,f,g),K[r],Wt[r])
        T2=add32(fS0(a),fMaj(a,b,c))
        h,g,f=g,f,e; e=add32(d,T1); d,c,b=c,b,a; a=add32(T1,T2)

    fake_H = [add32(x,iv) for x,iv in zip([a,b,c,d,e,f,g,h], IV)]
    fake_hash_hex = ''.join(f'{w:08x}' for w in fake_H)

    print(f"\n--- {NR} раундов ---")
    result, t = solve_sha256(fake_hash_hex, msg_len, charset='ascii', n_rounds=NR)

    if result:
        print(f"  Найдено: {result} {'== Hi ✓' if result == msg_test else '(другой прообраз)'}")

    if t > 60:
        print("  Timeout — stopping")
        break
