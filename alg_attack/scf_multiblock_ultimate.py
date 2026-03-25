#!/usr/bin/env python3
"""
SCF: MULTIBLOCK ULTIMATE — два прогона:
  Run A: SA-only multiblock (baseline)
  Run B: наши инструменты (adaptive δ + structure + biased targeting)

Сравниваем напрямую при ОДИНАКОВОМ бюджете.
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

def sha_compress(W16, iv=None):
    if iv is None: iv=IV
    W=expand_real(W16); s=list(iv)
    for r in range(64):
        a,b,c,d,e,f,g,h=s
        T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add32(Sig0(a),Maj(a,b,c))
        s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
    return [add32(iv[i],s[i]) for i in range(8)]

def two_block_hash(M1, M2):
    """Two-block hash: H = compress(compress(IV, M1), M2)."""
    H1 = sha_compress(M1)
    return sha_compress(M2, iv=H1)

def two_block_diff(M1, M1p, M2, M2p):
    """HW difference between two 2-block messages."""
    H_a = two_block_hash(M1, M2)
    H_b = two_block_hash(M1p, M2p)
    return sum(hw(H_a[i]^H_b[i]) for i in range(8))


# ============================================================
# RUN A: SA-only multiblock (blind baseline)
# ============================================================
def run_a_sa_multiblock(budget=3000):
    """Pure SA on two-block collision.
    Strategy: fix M1/M1' (random delta), optimize M2/M2'."""

    M1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    M1p = list(M1); M1p[0] ^= 0x80000000

    # Same M2 for both (simplest multiblock)
    M2 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    M2p = list(M2)  # Same M2 — different IV does the work

    H1 = sha_compress(M1); H1p = sha_compress(M1p)

    def score(m2):
        Ha = sha_compress(m2, iv=H1)
        Hb = sha_compress(m2, iv=H1p)
        return sum(hw(Ha[i]^Hb[i]) for i in range(8))

    best = score(M2); cur = list(M2)

    for it in range(budget):
        t = list(cur)
        w = int.from_bytes(os.urandom(1),'big') % 16
        b = int.from_bytes(os.urandom(1),'big') % 32
        t[w] ^= (1<<b)
        s = score(t)
        T = max(0.01, 1-it/budget)
        if s < best or math.exp(-(s-best)/(T*2)) > int.from_bytes(os.urandom(4),'big')/(1<<32):
            cur = t
            if s < best: best = s

    return best


# ============================================================
# RUN B: Our tools multiblock (structure-aware)
# ============================================================
def run_b_structured_multiblock(budget=3000):
    """Our tools:
    1. Screen δ for Block 1 (adaptive)
    2. Structure-solve Block 2 (exact inversion targeting)
    3. SA polish
    """
    # Phase 1: Screen Block 1 δ (50 evals)
    M1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    H1_base = sha_compress(M1)

    best_delta = (256, 0, 0)
    for word in range(4):  # Screen first 4 words
        for bit in range(32):
            M1p = list(M1); M1p[word] ^= (1<<bit)
            H1p = sha_compress(M1p)
            dH1 = sum(hw(H1_base[i]^H1p[i]) for i in range(8))
            if dH1 < best_delta[0]:
                best_delta = (dH1, word, bit)

    _, d_word, d_bit = best_delta
    M1p = list(M1); M1p[d_word] ^= (1<<d_bit)
    H1 = sha_compress(M1); H1p = sha_compress(M1p)
    evals = 128

    # Phase 2: Structure-solve Block 2
    M2 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

    def score(m2):
        Ha = sha_compress(m2, iv=H1)
        Hb = sha_compress(m2, iv=H1p)
        return sum(hw(Ha[i]^Hb[i]) for i in range(8))

    best = score(M2); cur = list(M2)
    evals += 1

    # Structure step: for each word, try exact-inversion-style corrections
    # Target: make compress(H1, M2) ≈ compress(H1p, M2)
    for struct_iter in range(10):
        improved = False
        for w in range(16):
            # Try ±1 additive correction
            for delta in [1, -1, 3, -3, 7, -7, 0xFF, -0xFF]:
                t = list(cur)
                t[w] = (cur[w] + delta) & MASK32
                s = score(t); evals += 1
                if s < best:
                    best = s; cur = list(t); improved = True
                    break
            if evals > budget // 2: break
        if not improved or evals > budget // 2: break

    # Phase 3: SA polish with remaining budget
    sa_budget = budget - evals
    for it in range(max(sa_budget, 0)):
        t = list(cur)
        w = int.from_bytes(os.urandom(1),'big') % 16
        b = int.from_bytes(os.urandom(1),'big') % 32
        t[w] ^= (1<<b)
        s = score(t)
        T = max(0.01, 1-it/max(sa_budget,1))
        if s < best or math.exp(-(s-best)/(T*2)) > int.from_bytes(os.urandom(4),'big')/(1<<32):
            cur = t
            if s < best: best = s

    return best


# ============================================================
# RUN C: Multi-start structured (try multiple Block 1 deltas)
# ============================================================
def run_c_multistart(budget=3000):
    """Try 5 different Block 1 pairs, pick best Block 2."""
    M1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    per_start = budget // 5

    global_best = 256

    # Screen top-5 Block1 deltas
    screen = []
    for word in range(8):
        for bit in range(32):
            M1p = list(M1); M1p[word] ^= (1<<bit)
            H1 = sha_compress(M1); H1p = sha_compress(M1p)
            dH = sum(hw(H1[i]^H1p[i]) for i in range(8))
            screen.append((dH, word, bit, H1, H1p))
    screen.sort()

    for dH, word, bit, H1, H1p in screen[:5]:
        M2 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

        def score(m2):
            Ha = sha_compress(m2, iv=H1)
            Hb = sha_compress(m2, iv=H1p)
            return sum(hw(Ha[i]^Hb[i]) for i in range(8))

        best = score(M2); cur = list(M2)

        for it in range(per_start):
            t = list(cur)
            w = int.from_bytes(os.urandom(1),'big') % 16
            b = int.from_bytes(os.urandom(1),'big') % 32
            t[w] ^= (1<<b)
            s = score(t)
            T = max(0.01, 1-it/per_start)
            if s < best or math.exp(-(s-best)/(T*2)) > int.from_bytes(os.urandom(4),'big')/(1<<32):
                cur = t
                if s < best: best = s

        if best < global_best:
            global_best = best

    return global_best


# ============================================================
if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    budget = 3000

    print("="*70)
    print(f"MULTIBLOCK SHOOTOUT (budget={budget}, N={N})")
    print("="*70)

    a_results = []; b_results = []; c_results = []

    for trial in range(N):
        a = run_a_sa_multiblock(budget)
        b = run_b_structured_multiblock(budget)
        c = run_c_multistart(budget)

        a_results.append(a); b_results.append(b); c_results.append(c)

        best_label = "A" if a <= b and a <= c else ("B" if b <= c else "C")
        print(f"  Trial {trial:2d}: A(SA)={a:3d}  B(struct)={b:3d}  C(multi-start)={c:3d}  → {best_label}")

    print(f"\n{'='*70}")
    print(f"SUMMARY")
    print(f"  A (SA-only):       avg={sum(a_results)/N:.1f}  min={min(a_results)}")
    print(f"  B (structured):    avg={sum(b_results)/N:.1f}  min={min(b_results)}")
    print(f"  C (multi-start):   avg={sum(c_results)/N:.1f}  min={min(c_results)}")
    print(f"")

    a_wins = sum(1 for i in range(N) if a_results[i] <= b_results[i] and a_results[i] <= c_results[i])
    b_wins = sum(1 for i in range(N) if b_results[i] < a_results[i] and b_results[i] <= c_results[i])
    c_wins = sum(1 for i in range(N) if c_results[i] < a_results[i] and c_results[i] < b_results[i])
    print(f"  Wins: A={a_wins}  B={b_wins}  C={c_wins}")

    overall_best = min(min(a_results), min(b_results), min(c_results))
    print(f"\n  OVERALL BEST: HW = {overall_best}")
    print(f"  Single-block best: HW ≈ 95")
    print(f"  Multi-block improvement: {95 - overall_best:+d}")
