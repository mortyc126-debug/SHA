#!/usr/bin/env python3
"""
SCF: FINAL PUSH — C+B combo + 3-block + increased budget.

Лучшие техники объединены:
1. Multi-start δ screening (C) — лучший avg
2. Structure corrections (B) — лучший min
3. 3-block — ещё 512 бит свободы
4. Increased budget — 10K evals для финального рекорда
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


def sa_optimize(score_fn, W_init, budget, W_ref=None):
    """Generic SA optimizer. score_fn(W) → int to minimize."""
    cur = list(W_init); best_s = score_fn(cur); best_W = list(cur); cur_s = best_s
    for it in range(budget):
        t = list(cur)
        w = int.from_bytes(os.urandom(1),'big') % 16
        b = int.from_bytes(os.urandom(1),'big') % 32
        t[w] ^= (1<<b)
        if W_ref and t == W_ref: continue
        s = score_fn(t)
        T = max(0.01, 1-it/budget)
        if s < cur_s or math.exp(-(s-cur_s)/(T*2)) > int.from_bytes(os.urandom(4),'big')/(1<<32):
            cur = t; cur_s = s
            if s < best_s: best_s = s; best_W = list(t)
    return best_W, best_s


# ============================================================
# COMBO C+B: Multi-start screen + structure + SA
# ============================================================
def combo_cb_2block(budget=5000):
    """Best 2-block: screen Block1 δ + structure Block2 + SA."""
    M1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

    # Screen Block1 deltas (256 evals for words 0-7)
    screen = []
    for word in range(8):
        for bit in range(32):
            M1p = list(M1); M1p[word] ^= (1<<bit)
            H1 = sha_compress(M1); H1p = sha_compress(M1p)
            dH = sum(hw(H1[i]^H1p[i]) for i in range(8))
            screen.append((dH, word, bit))
    screen.sort()

    global_best = 256
    per_start = (budget - 256) // 5

    for dH, word, bit in screen[:5]:
        M1p = list(M1); M1p[word] ^= (1<<bit)
        H1 = sha_compress(M1); H1p = sha_compress(M1p)

        def score(m2):
            return sum(hw(sha_compress(m2,iv=H1)[i] ^ sha_compress(m2,iv=H1p)[i]) for i in range(8))

        M2 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

        # Structure phase: additive corrections
        cur = list(M2); best = score(cur)
        for _ in range(min(per_start//10, 20)):
            improved = False
            for w in range(16):
                for d in [1,-1,3,-3,0xF,-0xF,0xFF]:
                    t = list(cur); t[w] = (cur[w]+d)&MASK32
                    s = score(t)
                    if s < best: best=s; cur=list(t); improved=True; break
                if improved: break
            if not improved: break

        # SA polish
        _, sa_best = sa_optimize(score, cur, per_start - per_start//10)

        if sa_best < global_best:
            global_best = sa_best

    return global_best


# ============================================================
# 3-BLOCK: Three blocks of freedom
# ============================================================
def three_block(budget=5000):
    """3-block collision attempt.
    M1/M1' → H1/H1' (block 1 near-collision)
    M2 shared → H2_a = compress(H1,M2), H2_b = compress(H1',M2)
    M3 shared → final_a = compress(H2_a,M3), final_b = compress(H2_b,M3)

    Optimize M2 AND M3 to minimize final diff.
    """
    M1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

    # Screen Block1 δ
    best_delta = (256, 0, 0)
    for word in range(4):
        for bit in range(32):
            M1p = list(M1); M1p[word] ^= (1<<bit)
            H1 = sha_compress(M1); H1p = sha_compress(M1p)
            dH = sum(hw(H1[i]^H1p[i]) for i in range(8))
            if dH < best_delta[0]:
                best_delta = (dH, word, bit)

    _, dw, db = best_delta
    M1p = list(M1); M1p[dw] ^= (1<<db)
    H1 = sha_compress(M1); H1p = sha_compress(M1p)

    M2 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    M3 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

    def score_3block(m2, m3):
        H2a = sha_compress(m2, iv=H1); H2b = sha_compress(m2, iv=H1p)
        H3a = sha_compress(m3, iv=H2a); H3b = sha_compress(m3, iv=H2b)
        return sum(hw(H3a[i]^H3b[i]) for i in range(8))

    best = score_3block(M2, M3)

    # Alternate: optimize M2, then M3, repeat
    half = budget // 4

    for phase in range(4):
        # Optimize M2
        def score_m2(m2): return score_3block(m2, M3)
        M2, s = sa_optimize(score_m2, M2, half)
        if s < best: best = s

        # Optimize M3
        H2a = sha_compress(M2, iv=H1); H2b = sha_compress(M2, iv=H1p)
        def score_m3(m3):
            return sum(hw(sha_compress(m3,iv=H2a)[i]^sha_compress(m3,iv=H2b)[i]) for i in range(8))
        M3, s = sa_optimize(score_m3, M3, half)
        if s < best: best = s

    return best


# ============================================================
# MAIN
# ============================================================
if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    budget = 5000

    print("="*70)
    print(f"SCF: FINAL PUSH (budget={budget}, N={N})")
    print("="*70)

    cb2_results = []; b3_results = []

    for trial in range(N):
        cb2 = combo_cb_2block(budget)
        b3 = three_block(budget)
        cb2_results.append(cb2); b3_results.append(b3)

        best = min(cb2, b3)
        label = "2B" if cb2 <= b3 else "3B"
        star = " ★" if best < 90 else ""
        print(f"  Trial {trial:2d}: 2-block(C+B)={cb2:3d}  3-block={b3:3d}  → {label}{star}")

    print(f"\n{'='*70}")
    print(f"FINAL RESULTS")
    print(f"  2-block (C+B): avg={sum(cb2_results)/N:.1f}  min={min(cb2_results)}")
    print(f"  3-block:       avg={sum(b3_results)/N:.1f}  min={min(b3_results)}")
    overall = min(min(cb2_results), min(b3_results))
    print(f"\n  ★ OVERALL BEST: HW = {overall}")
    print(f"  Session record: 128 → {overall} = {128-overall} bits improvement")
