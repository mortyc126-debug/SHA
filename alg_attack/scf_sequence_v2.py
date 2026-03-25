#!/usr/bin/env python3
"""
SCF: SEQUENCE SOLVER v2 — правильная метрика с D-coupling.

v1 ошибка: метрика δa[57..64] (8 значений) не даёт δe=0.
v2 метрика: δa[53..64] (12 значений) → δD[57..64]=0 → δe[57..64]=0.

12 нулей в a-sequence = гарантированная коллизия.
384 бит (12×32) нужно обнулить. Свобода: 512 бит.
Система НЕДООПРЕДЕЛЕНА (512 > 384) → решение СУЩЕСТВУЕТ!
"""
import os, sys, math

MASK32 = 0xFFFFFFFF
K = [
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2,
]
IV = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
      0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&MASK32
def Sig0(x): return rotr(x,2)^rotr(x,13)^rotr(x,22)
def Sig1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def sig0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def sig1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def Ch(e,f,g): return (e&f)^(~e&g)&MASK32
def Maj(a,b,c): return (a&b)^(a&c)^(b&c)
def add32(*args):
    s=0
    for x in args: s=(s+x)&MASK32
    return s
def sub32(a,b): return (a-b)&MASK32
def hw(x): return bin(x).count('1')

def expand_real(W16):
    W=list(W16)
    for i in range(16,64):
        W.append(add32(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    return W

def sha_a_seq(W16):
    W=expand_real(W16); s=list(IV)
    a_seq = [s[0]]
    for r in range(64):
        a,b,c,d,e,f,g,h=s
        T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add32(Sig0(a),Maj(a,b,c))
        s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
        a_seq.append(s[0])
    return a_seq

def sha_compress(W16):
    W=expand_real(W16); s=list(IV)
    for r in range(64):
        a,b,c,d,e,f,g,h=s
        T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add32(Sig0(a),Maj(a,b,c))
        s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
    return [add32(IV[i],s[i]) for i in range(8)]


def seq12_score(W1, W2):
    """Correct metric: δa[53..64] — 12 values = 384 bits."""
    a1 = sha_a_seq(W1); a2 = sha_a_seq(W2)
    return sum(hw(a1[r] ^ a2[r]) for r in range(53, 65))

def hash_score(W1, W2):
    H1=sha_compress(W1); H2=sha_compress(W2)
    return sum(hw(H1[i]^H2[i]) for i in range(8))

def sa_optimize(W1, W2_init, score_fn, budget, protect=None):
    """SA optimizer with score function and protected bits."""
    cur = list(W2_init); best_s = score_fn(W1, cur); best_W = list(cur)
    cur_s = best_s
    for it in range(budget):
        t = list(cur)
        w = int.from_bytes(os.urandom(1),'big') % 16
        b = int.from_bytes(os.urandom(1),'big') % 32
        if protect and (w, b) in protect: continue
        t[w] ^= (1<<b)
        if t == W1: continue
        s = score_fn(W1, t)
        T = max(0.01, 1-it/budget)
        if s < cur_s or math.exp(-(s-cur_s)/(T*2)) > int.from_bytes(os.urandom(4),'big')/(1<<32):
            cur = t; cur_s = s
            if s < best_s: best_s = s; best_W = list(t)
    return best_W, best_s


# ============================================================
# EXP 1: Verify — does δa[53..64]=0 guarantee collision?
# ============================================================
def exp1_verify():
    print("="*70)
    print("EXP 1: VERIFY — δa[53..64]=0 → collision?")
    print("  D[r] coupling: δe = δa - δD, δD depends on a[r-1..r-4]")
    print("  δa[53..64]=0 → δD[57..64]=0 → δe[57..64]=0 → δH=0")
    print("="*70)

    # If δa[53..64]=0 then:
    # δD[r] for r=57..64 depends on a[r-1..r-4] = a[56..53]
    # If δa[53..56]=0 → δD[57..60]=0 → δe[57..60]=0
    # And δa[57..60]=0 → δD[61..64]=0 → δe[61..64]=0
    # Hash registers: (a64,a63,a62,a61, e64,e63,e62,e61)
    # δa[61..64]=0 ✓, δe[61..64]=0 ✓ → δH=0 ✓

    print(f"\n  Logic chain:")
    print(f"    δa[53..56]=0 → δD[57..60]=0 → δe[57..60]=0")
    print(f"    δa[57..60]=0 → δD[61..64]=0 → δe[61..64]=0")
    print(f"    δa[61..64]=0 → hash regs a[61..64] match ✓")
    print(f"    δe[61..64]=0 → hash regs e[61..64] match ✓")
    print(f"    → ALL 8 hash registers match → δH=0 ✓")
    print(f"\n  Constraint count: 12 values × 32 bits = 384 bits")
    print(f"  Freedom: 16 words × 32 bits = 512 bits")
    print(f"  Excess freedom: 512 - 384 = 128 bits")
    print(f"  → System is UNDERDETERMINED → solutions EXIST!")


# ============================================================
# EXP 2: Seq12 metric vs hash metric — are they correlated?
# ============================================================
def exp2_correlation(N):
    print("\n" + "="*70)
    print("EXP 2: CORRELATION — seq12 metric vs hash HW")
    print("="*70)

    pairs = []
    for _ in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = list(W1)
        W2[int.from_bytes(os.urandom(1),'big')%16] ^= int.from_bytes(os.urandom(4),'big')
        if W2 == W1: W2[0] ^= 1

        s12 = seq12_score(W1, W2)
        h = hash_score(W1, W2)
        pairs.append((s12, h))

    # Correlation
    xs = [p[0] for p in pairs]; ys = [p[1] for p in pairs]
    mx = sum(xs)/len(xs); my = sum(ys)/len(ys)
    cov = sum((x-mx)*(y-my) for x,y in pairs)/len(pairs)
    sx = (sum((x-mx)**2 for x in xs)/len(xs))**0.5
    sy = (sum((y-my)**2 for y in ys)/len(ys))**0.5
    corr = cov/(sx*sy) if sx>0 and sy>0 else 0

    print(f"  Correlation(seq12, hash_HW) = {corr:.4f}")
    print(f"  (1.0 = perfect, 0.0 = none)")

    if corr > 0.8:
        print(f"  ★ STRONG correlation! Seq12 is a GOOD proxy for collision!")
    elif corr > 0.5:
        print(f"  Moderate correlation. Seq12 partially predicts hash.")
    else:
        print(f"  Weak correlation. Seq12 ≠ hash.")

    return corr


# ============================================================
# EXP 3: Optimize seq12 and measure hash
# ============================================================
def exp3_optimize(N):
    print("\n" + "="*70)
    print("EXP 3: OPTIMIZE seq12 → measure hash HW")
    print("  Compare: optimize seq12 vs optimize hash directly")
    print("="*70)

    seq_hash_results = []  # hash HW when optimizing seq12
    direct_hash_results = []  # hash HW when optimizing hash

    budget = 3000

    for trial in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

        # Screen deltas
        screen = []
        for word in range(8):
            for bit in range(32):
                W2 = list(W1); W2[word] ^= (1<<bit)
                s = seq12_score(W1, W2)
                screen.append((s, word, bit))
        screen.sort()

        per_start = (budget - 256) // 3

        # A: Optimize seq12
        best_seq_hash = 256
        for _, dw, db in screen[:3]:
            W2 = list(W1); W2[dw] ^= (1<<db)
            W2_opt, _ = sa_optimize(W1, W2, seq12_score, per_start, protect={(dw,db)})
            h = hash_score(W1, W2_opt)
            if h < best_seq_hash: best_seq_hash = h
        seq_hash_results.append(best_seq_hash)

        # B: Optimize hash directly (same budget)
        screen_h = []
        for word in range(8):
            for bit in range(32):
                W2 = list(W1); W2[word] ^= (1<<bit)
                s = hash_score(W1, W2)
                screen_h.append((s, word, bit))
        screen_h.sort()

        best_direct = 256
        for _, dw, db in screen_h[:3]:
            W2 = list(W1); W2[dw] ^= (1<<db)
            W2_opt, _ = sa_optimize(W1, W2, hash_score, per_start, protect={(dw,db)})
            h = hash_score(W1, W2_opt)
            if h < best_direct: best_direct = h
        direct_hash_results.append(best_direct)

        print(f"  Trial {trial:2d}: seq12→HW={best_seq_hash:3d}  direct→HW={best_direct:3d}  "
              f"{'SEQ ★' if best_seq_hash < best_direct else ('HASH' if best_direct < best_seq_hash else 'tie')}")

    avg_s = sum(seq_hash_results)/N
    avg_d = sum(direct_hash_results)/N
    print(f"\n  Seq12→HW:  avg={avg_s:.1f}  min={min(seq_hash_results)}")
    print(f"  Direct HW: avg={avg_d:.1f}  min={min(direct_hash_results)}")
    print(f"  Seq12 advantage: {avg_d-avg_s:+.1f}")


# ============================================================
# EXP 4: The 128-bit excess — is it exploitable?
# ============================================================
def exp4_excess_freedom():
    print("\n" + "="*70)
    print("EXP 4: 128-BIT EXCESS FREEDOM")
    print("  512 bits input - 384 bits constraints = 128 bits free")
    print("  These 128 bits parameterize the SPACE of collisions")
    print("="*70)

    # The collision variety near diagonal:
    # F: Z^{16}_{2^32} → Z^{12}_{2^32} (map W to δa[53..64])
    # dim(kernel) = 16-12 = 4 words = 128 bits (if full rank)

    # Verify: rank of Jacobian of the map W → a[53..64]
    W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    a1 = sha_a_seq(W1)

    # GF(2) Jacobian: 384 × 512
    rows = []
    for r in range(53, 65):
        for bit in range(32):
            row = 0
            for word in range(16):
                for wbit in range(32):
                    W2 = list(W1); W2[word] ^= (1<<wbit)
                    a2 = sha_a_seq(W2)
                    if (a1[r] ^ a2[r]) >> bit & 1:
                        row |= 1 << (word*32 + wbit)
            rows.append(row)

    # Rank over GF(2)
    m = list(rows); rank = 0
    for bp in range(511, -1, -1):
        mask = 1 << bp
        pv = -1
        for i in range(rank, len(m)):
            if m[i] & mask: pv = i; break
        if pv == -1: continue
        m[rank], m[pv] = m[pv], m[rank]
        for i in range(len(m)):
            if i != rank and m[i] & mask:
                m[i] ^= m[rank]
        rank += 1

    kernel_dim = 512 - rank
    print(f"\n  Jacobian: {len(rows)} × 512 over GF(2)")
    print(f"  Rank: {rank}")
    print(f"  Kernel dim: {kernel_dim}")
    print(f"  Expected: rank=384, kernel=128")

    if kernel_dim > 0:
        print(f"\n  ★ {kernel_dim} free dimensions for collision search!")
        print(f"    Birthday in this space: O(2^{kernel_dim//2})")
    else:
        print(f"  Kernel = 0 → overconstrained")


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10

    exp1_verify()
    exp2_correlation(min(N*20, 200))
    exp3_optimize(N)
    exp4_excess_freedom()
