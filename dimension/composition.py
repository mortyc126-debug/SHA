"""
COMPOSITION ТКАНЕЙ: SHA-256 ∘ SHA-256.

F² = F ∘ F: вход → SHA-256 → SHA-256 → выход.
В нашем измерении: двойная ткань = F ⊘ F (свёртка ткани с собой).

Вопросы:
  1. rank(CE) для F²: 256 (как F) или другой?
  2. Кривизна K для F²: 128 (как F) или выше/ниже?
  3. Если rank(CE_F²) < rank(CE_F): composition ОСЛАБЛЯЕТ!
  4. Security boundary F² vs F: такой же r=16?
"""

import numpy as np
import struct, hashlib

def sha256_words(W16):
    return struct.unpack('>8I', hashlib.sha256(struct.pack('>16I', *W16)).digest())

def sha256_double(W16):
    """F² = SHA-256(SHA-256(W16)) with padding."""
    H1 = sha256_words(W16)
    # H1 = 8 words. Pad to 16 words for second SHA-256.
    W2 = list(H1) + [0x80000000, 0, 0, 0, 0, 0, 0, 0x00000100]
    return sha256_words(W2)

def hw(x): return bin(x).count('1')


def compute_ranks(hash_fn, W_base):
    H_base = hash_fn(W_base)
    T_rows = []
    for word in range(16):
        for bit in range(32):
            W_mod = list(W_base)
            W_mod[word] ^= (1 << bit)
            H_mod = hash_fn(W_mod)
            row = []
            for w in range(8):
                d = H_base[w] ^ H_mod[w]
                for b in range(32): row.append((d >> b) & 1)
            T_rows.append(row)

    T = np.array(T_rows, dtype=np.uint8)
    rank_T = np.linalg.matrix_rank(T.astype(float))

    if rank_T < 256:
        return rank_T, None, None

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
    ces = []
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
        H2 = hash_fn(W2)
        ce = []
        for w in range(8):
            d = H_base[w]^H2[w]
            for b in range(32): ce.append((d>>b)&1)
        ces.append(ce)

    CE = np.array(ces[:256], dtype=np.uint8)
    rank_CE = np.linalg.matrix_rank(CE.astype(float))

    # Curvature
    Ks = []
    for _ in range(200):
        b1,b2 = np.random.choice(512, 2, replace=False)
        W1=list(W_base); W1[b1//32]^=(1<<(b1%32))
        W2=list(W_base); W2[b2//32]^=(1<<(b2%32))
        W12=list(W_base); W12[b1//32]^=(1<<(b1%32)); W12[b2//32]^=(1<<(b2%32))
        F1=hash_fn(W1); F2=hash_fn(W2); F12=hash_fn(W12)
        nonlin = sum(hw((H_base[i]^F12[i])^((H_base[i]^F1[i])^(H_base[i]^F2[i]))) for i in range(8))
        Ks.append(nonlin)

    return rank_T, rank_CE, np.mean(Ks)


def main():
    np.random.seed(42)

    print("=" * 70)
    print("COMPOSITION: SHA-256 ∘ SHA-256 (F²)")
    print("=" * 70)

    W_base = [np.random.randint(0, 2**32) for _ in range(16)]

    # Single SHA-256
    print(f"\n  Computing F (single SHA-256)...")
    rT_1, rCE_1, K_1 = compute_ranks(sha256_words, W_base)
    print(f"  F:  rank(T)={rT_1}, rank(CE)={rCE_1}, K={K_1:.1f}")

    # Double SHA-256
    print(f"\n  Computing F² (double SHA-256)...")
    rT_2, rCE_2, K_2 = compute_ranks(sha256_double, W_base)
    print(f"  F²: rank(T)={rT_2}, rank(CE)={rCE_2}, K={K_2:.1f}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("STABILITY: F² rank(CE) for multiple W_base")
    print("=" * 70)

    ranks_2 = []
    for trial in range(5):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        _, rCE, _ = compute_ranks(sha256_double, W)
        ranks_2.append(rCE)
        print(f"  Trial {trial}: rank(CE_F²) = {rCE}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("COMPARISON")
    print("=" * 70)

    print(f"""
  {'':>15} {'rank(T)':>8} {'rank(CE)':>9} {'K(curv)':>8}
  {'F (SHA-256)':>15} {rT_1:>7} {str(rCE_1):>9} {K_1:>7.1f}
  {'F² (double)':>15} {rT_2:>7} {str(rCE_2):>9} {K_2:>7.1f}

  Composition {'PRESERVES' if rCE_2 == rCE_1 else 'CHANGES'} rank(CE).
  Composition {'PRESERVES' if abs(K_2 - K_1) < 10 else 'CHANGES'} curvature K.

  В нашем измерении:
    F² = свёртка F с собой = "утолщённая" ткань (128 слоёв).
    {'Утолщение НЕ ОСЛАБЛЯЕТ (rank preserved).' if rCE_2 == 256 else 'Утолщение ОСЛАБЛЯЕТ! (rank decreased).'}
    {'SHA-256 СТАБИЛЬНА под composition.' if rCE_2 == 256 else 'SHA-256 НЕСТАБИЛЬНА!'}
""")


if __name__ == "__main__":
    main()
