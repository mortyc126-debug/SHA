#!/usr/bin/env python3
"""
ЗАДАНИЕ 9B: CTT Quiet Points — minimize nonlinearity Q(W, ΔW)

Q = HW(N) where N = ΔF ⊕ J·ΔW
ΔF = SHA256(W⊕ΔW) ⊕ SHA256(W)
J·ΔW = XOR of Jacobian columns for set bits in ΔW
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
def rw(n): return [random.getrandbits(32) for _ in range(n)]
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

def xor256(a, b):
    return tuple(a[i]^b[i] for i in range(8))

def compute_Q(W, DW):
    """Compute Q(W, ΔW) = HW(N) where N = ΔF ⊕ J·ΔW."""
    # ΔF = SHA256(W⊕ΔW) ⊕ SHA256(W)
    Wf = [(W[i] ^ DW[i]) & MASK for i in range(16)]
    H_w = sha256_hash(W)
    H_wf = sha256_hash(Wf)
    DF = xor256(H_w, H_wf)

    # J·ΔW = XOR of columns where ΔW has bit set
    JdW = [0]*8
    for word_idx in range(16):
        if DW[word_idx] == 0:
            continue
        for bit_idx in range(32):
            if (DW[word_idx] >> bit_idx) & 1:
                # Column for this input bit: SHA256(W ⊕ e_i) ⊕ SHA256(W)
                W_flip = list(W)
                W_flip[word_idx] ^= (1 << bit_idx)
                H_flip = sha256_hash(W_flip)
                col = xor256(H_w, H_flip)
                for k in range(8):
                    JdW[k] ^= col[k]

    JdW = tuple(JdW)
    N = xor256(DF, JdW)
    Q = hw256(N)
    return Q, DF, JdW, N


def compute_Q_fast(W, DW, H_w=None):
    """Faster Q computation — reuses H_w if provided."""
    if H_w is None:
        H_w = sha256_hash(W)
    Wf = [(W[i] ^ DW[i]) & MASK for i in range(16)]
    H_wf = sha256_hash(Wf)
    DF = xor256(H_w, H_wf)

    JdW = [0]*8
    for word_idx in range(16):
        if DW[word_idx] == 0:
            continue
        for bit_idx in range(32):
            if (DW[word_idx] >> bit_idx) & 1:
                W_flip = list(W)
                W_flip[word_idx] ^= (1 << bit_idx)
                H_flip = sha256_hash(W_flip)
                col = xor256(H_w, H_flip)
                for k in range(8):
                    JdW[k] ^= col[k]

    N = tuple(DF[k] ^ JdW[k] for k in range(8))
    return hw256(N), H_w


def sa_minimize_Q_fixed_DW(DW, steps=200000, restarts=30):
    """SA to minimize Q with fixed ΔW, varying W."""
    print(f"  SA minimizing Q with fixed ΔW (HW={hw256(DW)})")
    print(f"  Steps={steps}, restarts={restarts}")

    global_best_Q = 999
    global_best_W = None
    results = []

    for restart in range(restarts):
        W = rw(16)
        H_w = sha256_hash(W)
        cur_Q, H_w = compute_Q_fast(W, DW, H_w)
        best_Q = cur_Q
        best_W = list(W)

        T = 30.0
        T_end = 0.01
        alpha = (T_end / T) ** (1.0 / steps)

        for step in range(steps):
            # Mutate: flip random bit in W
            widx = random.randint(0, 15)
            bidx = random.randint(0, 31)
            W_new = list(W)
            W_new[widx] ^= (1 << bidx)

            new_Q, H_new = compute_Q_fast(W_new, DW)

            delta = new_Q - cur_Q
            if delta < 0 or random.random() < math.exp(-delta / T):
                W = W_new
                H_w = H_new
                cur_Q = new_Q
                if cur_Q < best_Q:
                    best_Q = cur_Q
                    best_W = list(W)

            T *= alpha

        results.append(best_Q)
        if best_Q < global_best_Q:
            global_best_Q = best_Q
            global_best_W = list(best_W)
            print(f"    Restart {restart}: Q={best_Q} *** NEW BEST ***")
        elif restart % 10 == 0:
            print(f"    Restart {restart}: Q={best_Q} (global best={global_best_Q})")

    return global_best_Q, global_best_W, results


def sa_minimize_Q_free_DW(steps=200000, restarts=30):
    """SA to minimize Q with free ΔW (at least 1 bit set)."""
    print(f"\n  SA minimizing Q with FREE ΔW")
    print(f"  Steps={steps}, restarts={restarts}")

    global_best_Q = 999
    global_best_W = None
    global_best_DW = None
    results = []

    for restart in range(restarts):
        W = rw(16)
        # Random ΔW with low HW (1-3 bits)
        DW = [0]*16
        n_bits = random.randint(1, 3)
        for _ in range(n_bits):
            DW[random.randint(0,15)] ^= (1 << random.randint(0,31))
        if all(d == 0 for d in DW):
            DW[0] = 1

        cur_Q, H_w = compute_Q_fast(W, DW)
        best_Q = cur_Q
        best_W = list(W)
        best_DW = list(DW)

        T = 30.0
        T_end = 0.01
        alpha = (T_end / T) ** (1.0 / steps)

        for step in range(steps):
            # 80% mutate W, 20% mutate DW
            if random.random() < 0.8:
                widx = random.randint(0, 15)
                bidx = random.randint(0, 31)
                W_new = list(W)
                W_new[widx] ^= (1 << bidx)
                DW_new = DW
            else:
                W_new = W
                DW_new = list(DW)
                widx = random.randint(0, 15)
                bidx = random.randint(0, 31)
                DW_new[widx] ^= (1 << bidx)
                if all(d == 0 for d in DW_new):
                    DW_new = list(DW)  # reject: DW must be nonzero

            new_Q, _ = compute_Q_fast(W_new, DW_new)

            delta = new_Q - cur_Q
            if delta < 0 or random.random() < math.exp(-delta / T):
                W = list(W_new)
                DW = list(DW_new)
                cur_Q = new_Q
                if cur_Q < best_Q:
                    best_Q = cur_Q
                    best_W = list(W)
                    best_DW = list(DW)

            T *= alpha

        results.append(best_Q)
        if best_Q < global_best_Q:
            global_best_Q = best_Q
            global_best_W = list(best_W)
            global_best_DW = list(best_DW)
            print(f"    Restart {restart}: Q={best_Q} HW(ΔW)={hw256(best_DW)} *** NEW BEST ***")
        elif restart % 10 == 0:
            print(f"    Restart {restart}: Q={best_Q} (global best={global_best_Q})")

    return global_best_Q, global_best_W, global_best_DW, results


def main():
    print("="*70)
    print("EXPERIMENT B: CTT Quiet Points — minimize nonlinearity Q")
    print("="*70)
    t0 = time.time()

    # Part 1: Baseline — Q for random (W, ΔW)
    # NOTE: For single-bit ΔW, Q=0 trivially (J column = exact derivative).
    # Nonlinearity only manifests with multi-bit ΔW.
    print("\n  Part 1: Baseline Q for random (W, ΔW) with multi-bit ΔW")
    for hw_target in [2, 4, 8, 16, 32]:
        baseline_Qs = []
        for _ in range(200):
            W = rw(16)
            DW = [0]*16
            for _ in range(hw_target):
                DW[random.randint(0,15)] ^= (1 << random.randint(0,31))
            if all(d == 0 for d in DW): DW[0] = 3
            Q, _, _, _ = compute_Q(W, DW)
            baseline_Qs.append(Q)
        avg_baseline = sum(baseline_Qs) / len(baseline_Qs)
        print(f"    HW(ΔW)≈{hw_target:2d}: avg Q = {avg_baseline:.1f}, "
              f"min = {min(baseline_Qs)}, max = {max(baseline_Qs)}")

    # Part 2: SA with fixed ΔW = [0x8000, 0, ..., 0] (Wang-like single word diff)
    print("\n  Part 2: SA with fixed ΔW = [0x8000, 0, ..., 0]")
    DW_fixed = [0]*16
    DW_fixed[0] = 0x8000  # 1-bit diff in bit 15 of word 0
    # For 1-bit ΔW, Q=0 trivially. Use 2-bit instead.
    print("    (1-bit ΔW gives Q=0 trivially — using 2-bit ΔW)")
    DW_fixed[0] = 0x8001  # 2 bits set → nonlinearity possible
    best_Q_fixed, best_W_fixed, results_fixed = sa_minimize_Q_fixed_DW(
        DW_fixed, steps=200000, restarts=30)

    print(f"\n  Fixed ΔW results:")
    print(f"    Best Q = {best_Q_fixed}")
    print(f"    Distribution: min={min(results_fixed)}, "
          f"max={max(results_fixed)}, avg={sum(results_fixed)/len(results_fixed):.1f}")

    # Part 3: SA with free ΔW
    print("\n  Part 3: SA with free ΔW")
    best_Q_free, best_W_free, best_DW_free, results_free = sa_minimize_Q_free_DW(
        steps=200000, restarts=30)

    print(f"\n  Free ΔW results:")
    print(f"    Best Q = {best_Q_free}")
    print(f"    Best ΔW HW = {hw256(best_DW_free)}")
    print(f"    Distribution: min={min(results_free)}, "
          f"max={max(results_free)}, avg={sum(results_free)/len(results_free):.1f}")

    # Part 4: Analysis of best point
    print("\n  Part 4: Analysis of best (W*, ΔW*)")
    if best_Q_free < best_Q_fixed:
        W_star, DW_star = best_W_free, best_DW_free
        print(f"    Best overall: Q={best_Q_free} (free ΔW)")
    else:
        W_star, DW_star = best_W_fixed, DW_fixed
        print(f"    Best overall: Q={best_Q_fixed} (fixed ΔW)")

    Q_star, DF_star, JdW_star, N_star = compute_Q(W_star, DW_star)
    print(f"    Q = {Q_star}")
    print(f"    HW(ΔF) = {hw256(DF_star)} (total difference)")
    print(f"    HW(J·ΔW) = {hw256(JdW_star)} (linear prediction)")
    print(f"    HW(N) = {hw256(N_star)} (nonlinear residual = Q)")
    print(f"    Linear accuracy: {256 - Q_star}/256 bits predicted correctly")

    print(f"\n    W*  = {' '.join(f'{x:08x}' for x in W_star)}")
    print(f"    ΔW* = {' '.join(f'{x:08x}' for x in DW_star)}")
    print(f"    ΔF  = {' '.join(f'{x:08x}' for x in DF_star)}")
    print(f"    J·ΔW= {' '.join(f'{x:08x}' for x in JdW_star)}")
    print(f"    N   = {' '.join(f'{x:08x}' for x in N_star)}")

    # Per-word nonlinearity
    print(f"\n    Per-word nonlinearity (Q per word):")
    for i in range(8):
        print(f"      H[{i}]: HW(N[{i}])={hw(N_star[i])}/32")

    elapsed = time.time() - t0
    print(f"\n  Total time: {elapsed:.1f}s")

    # Summary
    print(f"\n  SUMMARY:")
    print(f"    Baseline Q (random): {avg_baseline:.1f}")
    print(f"    Best Q (fixed ΔW=[1,0..]): {best_Q_fixed}")
    print(f"    Best Q (free ΔW): {best_Q_free}")
    print(f"    Reference from manual: Q_min=94 (SA)")
    if min(best_Q_fixed, best_Q_free) < 94:
        print(f"    *** IMPROVED over manual's Q=94! ***")
    else:
        print(f"    Consistent with manual's result")


if __name__ == "__main__":
    main()
