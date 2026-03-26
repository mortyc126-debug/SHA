"""
RANK CARRY ERROR: главный вопрос.

Carry error ≈ ЛИНЕЙНЫЙ оператор CE на GF(2)-kernel (разница 0.6 из 128).

CE: GF(2)^256 (kernel) → GF(2)^256 (carry error)
rank(CE) = ???

Если rank(CE) < 256: ker(CE) непуст → REAL COLLISION вычислима!
Если rank(CE) = 256: ker(CE) = {0} → нет бесплатной collision.

Строим матрицу CE из 256 базисных carry errors и вычисляем rank.
"""

import numpy as np
import struct, hashlib

MASK32 = 0xFFFFFFFF

def sha256_words(W16):
    return struct.unpack('>8I', hashlib.sha256(struct.pack('>16I', *W16)).digest())

def hw(x): return bin(x).count('1')


def main():
    np.random.seed(42)

    print("=" * 70)
    print("RANK CARRY ERROR: collision за polynomial time?")
    print("=" * 70)

    # Build T matrix
    W_base = [np.random.randint(0, 2**32) for _ in range(16)]
    H_base = sha256_words(W_base)

    print(f"\n  Строим матрицу T (512×256)...")
    T_rows = []
    for word in range(16):
        for bit in range(32):
            W_mod = list(W_base)
            W_mod[word] ^= (1 << bit)
            H_mod = sha256_words(W_mod)
            row = []
            for w in range(8):
                delta = H_base[w] ^ H_mod[w]
                for b in range(32):
                    row.append((delta >> b) & 1)
            T_rows.append(row)

    T = np.array(T_rows, dtype=np.uint8)  # 512 × 256

    # GF(2) kernel basis
    print(f"  GF(2) elimination...")
    M = T.T.copy()  # 256 × 512
    pivots = []
    row = 0
    for col in range(512):
        found = False
        for r in range(row, 256):
            if M[r, col] == 1:
                M[[row, r]] = M[[r, row]]
                found = True
                break
        if not found:
            continue
        pivots.append(col)
        for r in range(256):
            if r != row and M[r, col] == 1:
                M[r] = M[r] ^ M[row]
        row += 1

    free_vars = [c for c in range(512) if c not in pivots]
    print(f"  Kernel dimension: {len(free_vars)}")

    # Precompute all 256 basis kernel vectors
    print(f"  Precomputing {len(free_vars)} kernel basis vectors...")
    basis_vectors = []
    for fc_idx, fc in enumerate(free_vars):
        x = np.zeros(512, dtype=np.uint8)
        x[fc] = 1
        for i in range(len(pivots)-1, -1, -1):
            pc = pivots[i]
            val = np.uint8(0)
            for j in range(512):
                if j != pc:
                    val ^= (M[i, j] & x[j])
            x[pc] = val
        basis_vectors.append(x)

    # Compute carry error for each basis vector
    print(f"  Computing carry errors for {len(basis_vectors)} basis vectors...")
    carry_error_matrix = []  # 256 rows × 256 columns (CE basis)

    for x in basis_vectors:
        dW = [0] * 16
        for word in range(16):
            for bit in range(32):
                if x[word * 32 + bit]:
                    dW[word] ^= (1 << bit)

        W2 = [W_base[i] ^ dW[i] for i in range(16)]
        H2 = sha256_words(W2)

        # Carry error = actual δH (which should be 0 in GF(2) but isn't in Z/2^32)
        ce_bits = []
        for w in range(8):
            delta = H_base[w] ^ H2[w]
            for b in range(32):
                ce_bits.append((delta >> b) & 1)

        carry_error_matrix.append(ce_bits)

    CE = np.array(carry_error_matrix, dtype=np.uint8)  # 256 × 256

    print(f"\n  Carry error matrix CE: {CE.shape}")

    # RANK of CE
    rank_CE = np.linalg.matrix_rank(CE.astype(float))
    print(f"\n  ★★★ RANK(CE) = {rank_CE} ★★★")
    print(f"  Kernel dimension: 256")
    print(f"  If rank < 256: ker(CE) dimension = {256 - rank_CE}")

    if rank_CE < 256:
        print(f"\n  !!!  ker(CE) = {256 - rank_CE}-dimensional  !!!")
        print(f"  → {256 - rank_CE} independent REAL COLLISIONS exist!")
        print(f"  → Computable by Gaussian elimination = O(256³) = 2^24")
        print(f"  → SHA-256 BROKEN by carry-linear algebra???")

        # Find ker(CE) vector
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
            if not found:
                continue
            ce_pivots.append(col)
            for r in range(256):
                if r != ce_row and M_ce[r, col] == 1:
                    M_ce[r] = M_ce[r] ^ M_ce[ce_row]
            ce_row += 1

        ce_free = [c for c in range(256) if c not in ce_pivots]
        print(f"  CE kernel free vars: {len(ce_free)}")

        if ce_free:
            # Construct kernel vector of CE
            alpha = np.zeros(256, dtype=np.uint8)
            alpha[ce_free[0]] = 1
            for i in range(len(ce_pivots)-1, -1, -1):
                pc = ce_pivots[i]
                val = np.uint8(0)
                for j in range(256):
                    if j != pc:
                        val ^= (M_ce[i, j] & alpha[j])
                alpha[pc] = val

            # alpha defines: XOR of which basis kernel vectors
            # gives carry error = 0
            combo_x = np.zeros(512, dtype=np.uint8)
            for i in range(256):
                if alpha[i]:
                    combo_x = combo_x ^ basis_vectors[i]

            dW = [0] * 16
            for word in range(16):
                for bit in range(32):
                    if combo_x[word * 32 + bit]:
                        dW[word] ^= (1 << bit)

            if all(w == 0 for w in dW):
                print(f"  Trivial solution (δW=0). Trying next...")
                if len(ce_free) > 1:
                    alpha2 = np.zeros(256, dtype=np.uint8)
                    alpha2[ce_free[1]] = 1
                    for i in range(len(ce_pivots)-1, -1, -1):
                        pc = ce_pivots[i]
                        val = np.uint8(0)
                        for j in range(256):
                            if j != pc:
                                val ^= (M_ce[i, j] & alpha2[j])
                        alpha2[pc] = val

                    combo_x2 = np.zeros(512, dtype=np.uint8)
                    for i in range(256):
                        if alpha2[i]:
                            combo_x2 = combo_x2 ^ basis_vectors[i]

                    dW = [0] * 16
                    for word in range(16):
                        for bit in range(32):
                            if combo_x2[word * 32 + bit]:
                                dW[word] ^= (1 << bit)

            # Verify
            W2 = [W_base[i] ^ dW[i] for i in range(16)]
            H1 = sha256_words(W_base)
            H2 = sha256_words(W2)
            real_dh = sum(hw(H1[i] ^ H2[i]) for i in range(8))
            dh_words = [hw(H1[i] ^ H2[i]) for i in range(8)]

            print(f"\n  CE-kernel vector → δW:")
            print(f"    δW support: {sum(1 for w in dW if w != 0)}/16 nonzero")
            print(f"    δW HW total: {sum(hw(w) for w in dW)}")
            print(f"    Real SHA-256 δH: HW = {real_dh}")
            print(f"    Per word: {dh_words}")

            if real_dh == 0:
                print(f"\n  ★★★★★ COLLISION FOUND!!! ★★★★★")
                print(f"    W1 = {[hex(w) for w in W_base]}")
                print(f"    W2 = {[hex(w) for w in W2]}")
                print(f"    H1 = {[hex(w) for w in H1]}")
                print(f"    H2 = {[hex(w) for w in H2]}")
            else:
                print(f"\n  CE-kernel vector ≠ real collision (HW={real_dh})")
                print(f"  CE approximation error: carry IS almost linear but not exact.")

    else:
        print(f"\n  rank(CE) = 256 = FULL RANK")
        print(f"  ker(CE) = {{0}} → no free collision from linear algebra")
        print(f"  Carry error operator is INVERTIBLE on GF(2)-kernel")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("ВЕРДИКТ")
    print("=" * 70)

    print(f"""
  GF(2)-kernel of T: dimension 256
  Carry error operator CE: 256 × 256 matrix
  Rank(CE) = {rank_CE}

  {'CE IS SINGULAR! ker(CE) = ' + str(256 - rank_CE) + '-dimensional' if rank_CE < 256 else 'CE is full rank. No algebraic shortcut.'}

  Carry linearity test: actual vs predicted diff = 0.6/128 (99.5% linear!)
  But rank still = {rank_CE} → {'degenerate' if rank_CE < 256 else 'full rank despite linearity'}.

  {'NEXT: explore CE kernel for actual collision!' if rank_CE < 256 else 'NEXT: CE is approximation. Real CE nonlinear. Rank measurement includes noise.'}
""")


if __name__ == "__main__":
    main()
