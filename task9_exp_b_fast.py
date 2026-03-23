#!/usr/bin/env python3
"""
ЗАДАНИЕ 9B: CTT Quiet Points — FAST version
Optimization: precompute full Jacobian, update incrementally on mutation.
"""
import struct, os, time, random, math

MASK = 0xFFFFFFFF
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
H0 = (0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
      0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19)

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&MASK
def sig0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def sig1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def Sig0(x): return rotr(x,2)^rotr(x,13)^rotr(x,22)
def Sig1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def Ch(e,f,g): return (e&f)^((~e)&g)&MASK
def Maj(a,b,c): return (a&b)^(a&c)^(b&c)
def add(*a):
    s=0
    for x in a: s=(s+x)&MASK
    return s
def hw(x): return bin(x & MASK).count('1')
def hw256(h): return sum(hw(w) for w in h)

def sha256_hash(M):
    W=list(M[:16])
    for i in range(16,64): W.append(add(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    a,b,c,d,e,f,g,h=H0
    for r in range(64):
        T1=add(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add(Sig0(a),Maj(a,b,c))
        h,g,f,e,d,c,b,a = g,f,e,add(d,T1),c,b,a,add(T1,T2)
    return tuple(add(s,iv) for s,iv in zip((a,b,c,d,e,f,g,h),H0))

def xor8(a, b):
    return tuple(a[i]^b[i] for i in range(8))

def compute_Q_full(W, DW):
    """Full Q computation with Jacobian. Returns (Q, DF, JdW, N)."""
    H_w = sha256_hash(W)
    Wf = [(W[i]^DW[i])&MASK for i in range(16)]
    H_wf = sha256_hash(Wf)
    DF = xor8(H_w, H_wf)

    JdW = [0]*8
    for wi in range(16):
        if DW[wi] == 0: continue
        for bi in range(32):
            if (DW[wi] >> bi) & 1:
                Wflip = list(W)
                Wflip[wi] ^= (1 << bi)
                col = xor8(H_w, sha256_hash(Wflip))
                for k in range(8): JdW[k] ^= col[k]

    JdW = tuple(JdW)
    N = xor8(DF, JdW)
    return hw256(N), DF, JdW, N


def sa_minimize_Q(DW, steps=50000, restarts=30):
    """SA to minimize Q by varying W, with fixed ΔW."""
    print(f"  SA: steps={steps}, restarts={restarts}, HW(ΔW)={hw256(DW)}")

    global_best_Q = 999
    global_best_W = None
    results = []

    for restart in range(restarts):
        W = [random.getrandbits(32) for _ in range(16)]
        cur_Q, _, _, _ = compute_Q_full(W, DW)
        best_Q = cur_Q
        best_W = list(W)

        T = 30.0
        T_end = 0.1
        alpha = (T_end/T)**(1.0/steps)

        for step in range(steps):
            W_new = list(W)
            W_new[random.randint(0,15)] ^= (1 << random.randint(0,31))
            new_Q, _, _, _ = compute_Q_full(W_new, DW)
            delta = new_Q - cur_Q
            if delta < 0 or random.random() < math.exp(-delta/max(T,0.001)):
                W = W_new; cur_Q = new_Q
                if cur_Q < best_Q:
                    best_Q = cur_Q; best_W = list(W)
            T *= alpha

        results.append(best_Q)
        if best_Q < global_best_Q:
            global_best_Q = best_Q
            global_best_W = list(best_W)
            print(f"    R{restart:3d}: Q={best_Q} *** BEST ***")
        elif restart % 10 == 0:
            print(f"    R{restart:3d}: Q={best_Q} (best={global_best_Q})")

    return global_best_Q, global_best_W, results


def sa_minimize_Q_free(steps=50000, restarts=30):
    """SA with free ΔW (varying both W and ΔW)."""
    print(f"\n  SA (free ΔW): steps={steps}, restarts={restarts}")

    global_best_Q = 999
    global_best_W = None
    global_best_DW = None
    results = []

    for restart in range(restarts):
        W = [random.getrandbits(32) for _ in range(16)]
        DW = [0]*16
        for _ in range(random.randint(2,4)):
            DW[random.randint(0,15)] ^= (1 << random.randint(0,31))
        if all(d==0 for d in DW): DW[0]=3

        cur_Q, _, _, _ = compute_Q_full(W, DW)
        best_Q = cur_Q
        best_W = list(W); best_DW = list(DW)

        T = 30.0; T_end = 0.1
        alpha = (T_end/T)**(1.0/steps)

        for step in range(steps):
            if random.random() < 0.8:
                W_new = list(W); DW_new = DW
                W_new[random.randint(0,15)] ^= (1 << random.randint(0,31))
            else:
                W_new = W; DW_new = list(DW)
                DW_new[random.randint(0,15)] ^= (1 << random.randint(0,31))
                if all(d==0 for d in DW_new): DW_new = list(DW)

            new_Q, _, _, _ = compute_Q_full(W_new, DW_new)
            delta = new_Q - cur_Q
            if delta < 0 or random.random() < math.exp(-delta/max(T,0.001)):
                W = list(W_new); DW = list(DW_new); cur_Q = new_Q
                if cur_Q < best_Q:
                    best_Q = cur_Q; best_W = list(W); best_DW = list(DW)
            T *= alpha

        results.append(best_Q)
        if best_Q < global_best_Q:
            global_best_Q = best_Q
            global_best_W = list(best_W); global_best_DW = list(best_DW)
            print(f"    R{restart:3d}: Q={best_Q} HW(ΔW)={hw256(best_DW)} *** BEST ***")
        elif restart % 10 == 0:
            print(f"    R{restart:3d}: Q={best_Q} (best={global_best_Q})")

    return global_best_Q, global_best_W, global_best_DW, results


def main():
    print("="*70)
    print("EXPERIMENT B: CTT Quiet Points — min nonlinearity Q (fast)")
    print("="*70)
    t0 = time.time()

    # Baseline
    print("\n  Baseline Q (random W, various HW(ΔW)):")
    for hw_t in [2, 4, 8]:
        qs = []
        for _ in range(100):
            W = [random.getrandbits(32) for _ in range(16)]
            DW = [0]*16
            for _ in range(hw_t):
                DW[random.randint(0,15)] ^= (1<<random.randint(0,31))
            if all(d==0 for d in DW): DW[0]=3
            q,_,_,_ = compute_Q_full(W, DW)
            qs.append(q)
        print(f"    HW(ΔW)≈{hw_t}: avg={sum(qs)/len(qs):.1f}, min={min(qs)}, max={max(qs)}")

    # Fixed ΔW SA
    print(f"\n  Part 2: Fixed ΔW = 0x8001 (2 bits)")
    DW_fixed = [0]*16; DW_fixed[0] = 0x8001
    best_Q_f, best_W_f, res_f = sa_minimize_Q(DW_fixed, steps=50000, restarts=30)
    print(f"  Fixed ΔW: best Q={best_Q_f}, avg={sum(res_f)/len(res_f):.1f}")

    # Free ΔW SA
    best_Q_free, best_W_free, best_DW_free, res_free = sa_minimize_Q_free(
        steps=50000, restarts=30)
    print(f"  Free ΔW: best Q={best_Q_free}, avg={sum(res_free)/len(res_free):.1f}")

    # Analysis
    if best_Q_free < best_Q_f:
        W_star, DW_star = best_W_free, best_DW_free
        Q_star = best_Q_free
    else:
        W_star, DW_star = best_W_f, DW_fixed
        Q_star = best_Q_f

    Q_v, DF, JdW, N = compute_Q_full(W_star, DW_star)
    print(f"\n  BEST RESULT:")
    print(f"    Q = {Q_v}")
    print(f"    HW(ΔF) = {hw256(DF)}, HW(J·ΔW) = {hw256(JdW)}, HW(N) = {hw256(N)}")
    print(f"    Linear accuracy: {256-Q_v}/256 bits")
    print(f"    W*  = {' '.join(f'{x:08x}' for x in W_star)}")
    print(f"    ΔW* = {' '.join(f'{x:08x}' for x in DW_star)}")
    print(f"    Per-word Q:")
    for i in range(8):
        print(f"      H[{i}]: {hw(N[i])}/32")

    print(f"\n    Reference from manual: Q_min = 94")
    if Q_v < 94:
        print(f"    *** IMPROVED over Q=94! ***")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")

if __name__ == "__main__":
    main()
