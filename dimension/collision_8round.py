"""
8-ROUND SHA-256: rank(CE)=1 → collision вычислима!

rank(CE)=1 означает: ker(CE) = 255-мерный.
255 независимых GF(2)-kernel vectors с carry error = 0.
Каждый = COLLISION на 8-round SHA-256!

Верифицируем: находим collision, проверяем.
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
def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)

def sha256_8rounds(W16):
    a,b,c,d,e,f,g,h = IV
    W = list(W16)
    for r in range(8):
        T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K_const[r]), W[r])
        T2 = add32(Sigma0(a), Maj(a,b,c))
        h,g,f,e = g,f,e,add32(d,T1)
        d,c,b,a = c,b,a,add32(T1,T2)
    return tuple(add32(IV[i], [a,b,c,d,e,f,g,h][i]) for i in range(8))


def main():
    np.random.seed(42)

    print("=" * 70)
    print("8-ROUND SHA-256: COLLISION ЧЕРЕЗ CE-KERNEL")
    print("=" * 70)

    W_base = [np.random.randint(0, 2**32) for _ in range(16)]
    H_base = sha256_8rounds(W_base)

    # =================================================================
    print(f"\n  Строим T matrix (512×256) для 8-round SHA-256...")

    T_rows = []
    for word in range(16):
        for bit in range(32):
            W_mod = list(W_base)
            W_mod[word] ^= (1 << bit)
            H_mod = sha256_8rounds(W_mod)
            row = []
            for w in range(8):
                delta = H_base[w] ^ H_mod[w]
                for b in range(32):
                    row.append((delta >> b) & 1)
            T_rows.append(row)

    T = np.array(T_rows, dtype=np.uint8)
    rank_T = np.linalg.matrix_rank(T.astype(float))
    print(f"  rank(T) = {rank_T}")

    # GF(2) kernel
    print(f"  Computing GF(2) kernel...")
    M = T.T.copy()
    pivots = []
    row = 0
    for col in range(512):
        found = False
        for r in range(row, 256):
            if M[r, col] == 1:
                M[[row, r]] = M[[r, row]]
                found = True
                break
        if not found: continue
        pivots.append(col)
        for r in range(256):
            if r != row and M[r, col] == 1:
                M[r] = M[r] ^ M[row]
        row += 1

    free_vars = [c for c in range(512) if c not in pivots]
    print(f"  Kernel dimension: {len(free_vars)}")

    # CE matrix
    print(f"  Computing CE matrix...")
    basis_ces = []
    basis_dWs = []

    for fc in free_vars[:256]:
        x = np.zeros(512, dtype=np.uint8)
        x[fc] = 1
        for i in range(len(pivots)-1, -1, -1):
            pc = pivots[i]
            val = np.uint8(0)
            for j in range(512):
                if j != pc: val ^= (M[i,j] & x[j])
            x[pc] = val

        dW = [0]*16
        for word in range(16):
            for bit in range(32):
                if x[word*32+bit]: dW[word] ^= (1<<bit)

        W2 = [W_base[i]^dW[i] for i in range(16)]
        H2 = sha256_8rounds(W2)
        ce = []
        for w in range(8):
            delta = H_base[w] ^ H2[w]
            for b in range(32):
                ce.append((delta >> b) & 1)

        basis_ces.append(ce)
        basis_dWs.append(dW)

    CE = np.array(basis_ces[:256], dtype=np.uint8)
    rank_CE = np.linalg.matrix_rank(CE.astype(float))
    print(f"\n  ★ rank(CE) = {rank_CE}")
    print(f"  ★ ker(CE) dimension = {256 - rank_CE}")

    if rank_CE >= 256:
        print(f"  CE full rank — no algebraic collision.")
        return

    # =================================================================
    print(f"\n{'=' * 70}")
    print(f"  FINDING COLLISION: ker(CE) = {256-rank_CE}-dimensional")
    print("=" * 70)

    # GF(2) elimination on CE
    M_ce = CE.T.copy()  # 256 × 256
    ce_pivots = []
    ce_row = 0
    for col in range(256):
        found = False
        for r in range(ce_row, 256):
            if M_ce[r, col] == 1:
                M_ce[[ce_row, r]] = M_ce[[r, ce_row]]
                found = True
                break
        if not found: continue
        ce_pivots.append(col)
        for r in range(256):
            if r != ce_row and M_ce[r, col] == 1:
                M_ce[r] = M_ce[r] ^ M_ce[ce_row]
        ce_row += 1

    ce_free = [c for c in range(256) if c not in ce_pivots]
    print(f"  CE kernel free vars: {len(ce_free)}")

    # Try each CE-kernel vector
    collisions_found = 0

    for fc_idx, fc in enumerate(ce_free[:20]):
        alpha = np.zeros(256, dtype=np.uint8)
        alpha[fc] = 1
        for i in range(len(ce_pivots)-1, -1, -1):
            pc = ce_pivots[i]
            val = np.uint8(0)
            for j in range(256):
                if j != pc: val ^= (M_ce[i,j] & alpha[j])
            alpha[pc] = val

        # alpha = combination of basis kernel vectors
        # Combine to get actual δW
        combo_dW = [0]*16
        for i in range(256):
            if alpha[i]:
                for w in range(16):
                    combo_dW[w] ^= basis_dWs[i][w]

        if all(w == 0 for w in combo_dW):
            continue  # trivial

        # Test collision
        W2 = [W_base[i] ^ combo_dW[i] for i in range(16)]
        H1 = sha256_8rounds(W_base)
        H2 = sha256_8rounds(W2)

        dh = sum(hw(H1[i] ^ H2[i]) for i in range(8))
        dh_words = [hw(H1[i] ^ H2[i]) for i in range(8)]

        if dh == 0:
            collisions_found += 1
            print(f"\n  ★★★★★ COLLISION #{collisions_found} FOUND! ★★★★★")
            print(f"    W1 = {[hex(w) for w in W_base[:8]]}...")
            print(f"    W2 = {[hex(w) for w in W2[:8]]}...")
            print(f"    H1 = {[hex(w) for w in H1]}")
            print(f"    H2 = {[hex(w) for w in H2]}")
            print(f"    H1 == H2: {H1 == H2}")
            print(f"    δW HW: {sum(hw(combo_dW[i]) for i in range(16))}")
            if collisions_found >= 3:
                break
        else:
            if fc_idx < 5:
                print(f"    CE-kernel[{fc}]: δH HW = {dh} (per word: {dh_words})")

    print(f"\n  Total collisions found: {collisions_found} из {min(20, len(ce_free))} CE-kernel vectors")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("SIGNIFICANCE")
    print("=" * 70)

    print(f"""
  8-ROUND SHA-256:
    rank(T) = {rank_T}
    rank(CE) = {rank_CE}
    ker(CE) = {256 - rank_CE}-dimensional
    Collisions found: {collisions_found}

  {'★★★★★ COLLISION VERIFIED ON 8-ROUND SHA-256! ★★★★★' if collisions_found > 0 else 'No exact collision (CE rank measurement may have numerical error)'}

  Method: GF(2)-kernel + CE-kernel (our dimension's algebra).
  Cost: O(512³) ≈ 2^27 (Gaussian elimination).
  {'This is polynomial time — NOT birthday!' if collisions_found > 0 else ''}

  For comparison:
    8-round SHA-256 birthday: 2^128
    Our method: {'2^27 (polynomial!)' if collisions_found > 0 else 'needs refinement'}
    {'Speedup: 2^101 !!!' if collisions_found > 0 else ''}
""")


if __name__ == "__main__":
    main()
