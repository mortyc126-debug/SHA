"""
NO ROTATIONS SHA-256: rank(CE)=255 → ker(CE)=1 → ONE collision!

SHA-256 без Σ₀, Σ₁ (rotations = identity): rank(CE)=255.
ker(CE) = 1-dimensional → ОДНА collision direction.

Это collision на MODIFIED SHA-256 (no rotations).
Не стандартная SHA-256, но показывает роль ротаций.
"""

import numpy as np
import struct

MASK32 = 0xFFFFFFFF
K_const = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
]
IV = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]

def hw(x): return bin(x).count('1')
def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def add32(x, y): return (x + y) & MASK32
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)

def sha256_no_rot(W16, n_rounds=16):
    """SHA-256 with Σ₀=Σ₁=identity (no rotations)."""
    a,b,c,d,e,f,g,h = IV
    for r in range(n_rounds):
        T1 = add32(add32(add32(add32(h, e), Ch(e,f,g)), K_const[r%16]), W16[r%16])
        T2 = add32(a, Maj(a,b,c))  # Σ₀(a) = a
        h,g,f,e = g,f,e,add32(d,T1)
        d,c,b,a = c,b,a,add32(T1,T2)
    return tuple(add32(IV[i], [a,b,c,d,e,f,g,h][i]) for i in range(8))


def main():
    np.random.seed(42)
    print("=" * 70)
    print("NO ROTATIONS: rank(CE)=255 → FINDING THE ONE COLLISION")
    print("=" * 70)

    W_base = [np.random.randint(0, 2**32) for _ in range(16)]
    H_base = sha256_no_rot(W_base)

    # T matrix
    T_rows = []
    for word in range(16):
        for bit in range(32):
            W_mod = list(W_base)
            W_mod[word] ^= (1 << bit)
            H_mod = sha256_no_rot(W_mod)
            row = []
            for w in range(8):
                d = H_base[w] ^ H_mod[w]
                for b in range(32): row.append((d >> b) & 1)
            T_rows.append(row)

    T = np.array(T_rows, dtype=np.uint8)
    print(f"\n  rank(T) = {np.linalg.matrix_rank(T.astype(float))}")

    # GF2 kernel
    M = T.T.copy()
    pivots = []; row = 0
    for col in range(512):
        found = False
        for r in range(row, 256):
            if M[r,col]==1: M[[row,r]]=M[[r,row]]; found=True; break
        if not found: continue
        pivots.append(col)
        for r in range(256):
            if r!=row and M[r,col]==1: M[r]^=M[row]
        row += 1

    free_vars = [c for c in range(512) if c not in pivots]
    print(f"  GF2 kernel dim = {len(free_vars)}")

    # CE matrix
    basis_vecs = []
    basis_ces = []
    for fc in free_vars[:256]:
        x = np.zeros(512, dtype=np.uint8); x[fc] = 1
        for i in range(len(pivots)-1, -1, -1):
            pc = pivots[i]; val = np.uint8(0)
            for j in range(512):
                if j!=pc: val ^= (M[i,j]&x[j])
            x[pc] = val
        dW = [0]*16
        for word in range(16):
            for bit in range(32):
                if x[word*32+bit]: dW[word]^=(1<<bit)
        W2 = [W_base[i]^dW[i] for i in range(16)]
        H2 = sha256_no_rot(W2)
        ce = []
        for w in range(8):
            d = H_base[w]^H2[w]
            for b in range(32): ce.append((d>>b)&1)
        basis_vecs.append((x, dW))
        basis_ces.append(ce)

    CE = np.array(basis_ces[:256], dtype=np.uint8)
    rank_CE = np.linalg.matrix_rank(CE.astype(float))
    print(f"  rank(CE) = {rank_CE}")
    print(f"  ker(CE) = {256 - rank_CE}")

    if rank_CE < 256:
        # Find CE kernel
        M_ce = CE.T.copy()
        ce_pivots = []; ce_row = 0
        for col in range(256):
            found = False
            for r in range(ce_row, 256):
                if M_ce[r,col]==1: M_ce[[ce_row,r]]=M_ce[[r,ce_row]]; found=True; break
            if not found: continue
            ce_pivots.append(col)
            for r in range(256):
                if r!=ce_row and M_ce[r,col]==1: M_ce[r]^=M_ce[ce_row]
            ce_row += 1

        ce_free = [c for c in range(256) if c not in ce_pivots]
        print(f"  CE kernel free vars: {len(ce_free)}")

        for fc in ce_free:
            alpha = np.zeros(256, dtype=np.uint8); alpha[fc] = 1
            for i in range(len(ce_pivots)-1, -1, -1):
                pc = ce_pivots[i]; val = np.uint8(0)
                for j in range(256):
                    if j!=pc: val ^= (M_ce[i,j]&alpha[j])
                alpha[pc] = val

            combo_dW = [0]*16
            for i in range(256):
                if alpha[i]:
                    for w in range(16): combo_dW[w] ^= basis_vecs[i][1][w]

            if all(w==0 for w in combo_dW): continue

            W2 = [W_base[i]^combo_dW[i] for i in range(16)]
            H1 = sha256_no_rot(W_base)
            H2 = sha256_no_rot(W2)
            dh = sum(hw(H1[i]^H2[i]) for i in range(8))

            active_words = sum(1 for w in combo_dW if w != 0)
            total_hw = sum(hw(w) for w in combo_dW)

            print(f"\n  CE-kernel vector:")
            print(f"    δW: {active_words}/16 words active, HW={total_hw}")
            print(f"    HW(δH) = {dh}")
            if dh == 0:
                print(f"\n    ★★★ COLLISION ON NO-ROTATION SHA-256! ★★★")
                print(f"    H1 = {[hex(w) for w in H1]}")
                print(f"    H2 = {[hex(w) for w in H2]}")
                print(f"    H1 == H2: {H1 == H2}")

                # Which words differ?
                print(f"    δW nonzero words: ", end="")
                for w in range(16):
                    if combo_dW[w] != 0:
                        print(f"W[{w}]={hex(combo_dW[w])} ", end="")
                print()
            else:
                print(f"    Not exact collision (carry residual = {dh})")

    # Stability check
    print(f"\n  Stability: rank(CE) for 10 random W_base:")
    for trial in range(10):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H = sha256_no_rot(W)
        T_r = []
        for word in range(16):
            for bit in range(32):
                W_m = list(W); W_m[word] ^= (1<<bit)
                H_m = sha256_no_rot(W_m)
                row = []
                for w in range(8):
                    d = H[w]^H_m[w]
                    for b in range(32): row.append((d>>b)&1)
                T_r.append(row)
        T_r = np.array(T_r, dtype=np.uint8)
        rT = np.linalg.matrix_rank(T_r.astype(float))
        # Quick CE check
        if rT == 256:
            print(f"    Trial {trial}: rank(T)={rT} (CE=check needed)")
        else:
            print(f"    Trial {trial}: rank(T)={rT}")


if __name__ == "__main__":
    main()
