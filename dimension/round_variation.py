"""
ТКАНЬ С РАЗНЫМ ЧИСЛОМ СЛОЁВ: как геометрия зависит от глубины?

SHA-256 = 64 слоя → K=128, g=128, полная сфера.
А если 8, 16, 32 слоёв? Когда ткань становится "сферой"?

Это покажет: СКОЛЬКО слоёв нужно для максимальной кривизны.
Из Этапа 2.1: диффузия 90% за 5 раундов.
Из gap anatomy: amplification 4 раунда.
Ожидание: K достигает 128 за ~8-16 слоёв.
"""

import numpy as np
import struct, hashlib

MASK32 = 0xFFFFFFFF
K_const = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
    0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
    0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
    0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
    0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
    0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
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
def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)

def sha256_n_rounds(W16, n_rounds):
    """SHA-256 с N раундами."""
    W = list(W16)
    for r in range(16, max(n_rounds, 16)):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))
    a,b,c,d,e,f,g,h = IV
    for r in range(n_rounds):
        T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K_const[r % 64]), W[r])
        T2 = add32(Sigma0(a), Maj(a,b,c))
        h,g,f,e = g,f,e,add32(d,T1)
        d,c,b,a = c,b,a,add32(T1,T2)
    return tuple(add32(IV[i], [a,b,c,d,e,f,g,h][i]) for i in range(8))


def main():
    np.random.seed(42)

    print("=" * 70)
    print("ТКАНЬ С РАЗНЫМ ЧИСЛОМ СЛОЁВ")
    print("=" * 70)

    # =================================================================
    print(f"\n{'=' * 70}")
    print("1. КРИВИЗНА K vs число раундов")
    print("=" * 70)

    print(f"\n  {'Rounds':>6} {'K (curv)':>9} {'g (stretch)':>12} {'rank(T)':>8} {'Sphere?':>8}")

    for n_rounds in [1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 32, 48, 64]:
        # Кривизна
        Ks = []
        gs = []
        for _ in range(500):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            H_base = sha256_n_rounds(W, n_rounds)

            b1 = np.random.randint(0, 512)
            b2 = np.random.randint(0, 512)
            if b1 == b2: continue

            W1 = list(W); W1[b1//32] ^= (1<<(b1%32))
            W2 = list(W); W2[b2//32] ^= (1<<(b2%32))
            W12 = list(W); W12[b1//32] ^= (1<<(b1%32)); W12[b2//32] ^= (1<<(b2%32))

            H1 = sha256_n_rounds(W1, n_rounds)
            H2 = sha256_n_rounds(W2, n_rounds)
            H12 = sha256_n_rounds(W12, n_rounds)

            F1 = tuple(H_base[i] ^ H1[i] for i in range(8))
            F2 = tuple(H_base[i] ^ H2[i] for i in range(8))
            F12 = tuple(H_base[i] ^ H12[i] for i in range(8))
            nonlin = sum(hw(F12[i] ^ (F1[i]^F2[i])) for i in range(8))
            Ks.append(nonlin)

            # Stretch
            gs.append(sum(hw(F1[i]) for i in range(8)))

        K_mean = np.mean(Ks)
        g_mean = np.mean(gs)

        # rank(T) — approximate from stretch distribution
        # If g < 128: not fully diffused → rank < 256
        est_rank = min(256, int(g_mean * 2))

        sphere = "★" if K_mean > 120 and g_mean > 120 else ""

        print(f"  {n_rounds:6d} {K_mean:8.1f} {g_mean:11.1f} {est_rank:>7} {sphere:>8}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("2. ПЕРЕХОД К СФЕРЕ: когда K достигает максимума?")
    print("=" * 70)

    # Найдём точку перехода: K > 0.9 × K_max
    transition = None
    for n_r in range(1, 65):
        Ks = []
        for _ in range(200):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            H = sha256_n_rounds(W, n_r)
            b1, b2 = np.random.choice(512, 2, replace=False)
            W1 = list(W); W1[b1//32] ^= (1<<(b1%32))
            W2 = list(W); W2[b2//32] ^= (1<<(b2%32))
            W12 = list(W); W12[b1//32] ^= (1<<(b1%32)); W12[b2//32] ^= (1<<(b2%32))
            H_b = sha256_n_rounds(W, n_r)
            H1 = sha256_n_rounds(W1, n_r)
            H2 = sha256_n_rounds(W2, n_r)
            H12 = sha256_n_rounds(W12, n_r)
            nonlin = sum(hw((H_b[i]^H12[i]) ^ ((H_b[i]^H1[i])^(H_b[i]^H2[i]))) for i in range(8))
            Ks.append(nonlin)

        K_mean = np.mean(Ks)
        if K_mean > 115 and transition is None:
            transition = n_r
            print(f"  ПЕРЕХОД на r={n_r}: K={K_mean:.1f} (>115)")
            break

    if transition:
        print(f"\n  ★ SHA-256 становится СФЕРОЙ на раунде {transition}.")
        print(f"  Раунды 1..{transition-1}: кривизна нарастает.")
        print(f"  Раунды {transition}..64: постоянная кривизна (сфера).")
        print(f"  {64 - transition} из 64 раундов = ИЗБЫТОЧНЫ для геометрии.")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("3. CE RANK vs число раундов")
    print("=" * 70)

    # Для каких раундов rank(CE) < 256? (= algebraic weakness)
    for n_r in [4, 6, 8, 10, 12, 16, 20, 32, 64]:
        W_base = [np.random.randint(0, 2**32) for _ in range(16)]
        H_base = sha256_n_rounds(W_base, n_r)

        # T matrix
        T_rows = []
        for word in range(16):
            for bit in range(32):
                W_mod = list(W_base)
                W_mod[word] ^= (1 << bit)
                H_mod = sha256_n_rounds(W_mod, n_r)
                row = []
                for w in range(8):
                    delta = H_base[w] ^ H_mod[w]
                    for b in range(32):
                        row.append((delta >> b) & 1)
                T_rows.append(row)

        T = np.array(T_rows, dtype=np.uint8)
        rank_T = np.linalg.matrix_rank(T.astype(float))

        if rank_T < 256:
            print(f"  {n_r:3d} rounds: rank(T) = {rank_T} ← RANK DEFICIENT!")
        else:
            # Compute CE rank
            M = T.T.copy()
            pivots = []
            row = 0
            for col in range(512):
                found = False
                for r in range(row, min(256, M.shape[0])):
                    if M[r, col] == 1:
                        M[[row, r]] = M[[r, row]]
                        found = True
                        break
                if not found: continue
                pivots.append(col)
                for r in range(M.shape[0]):
                    if r != row and M[r, col] == 1:
                        M[r] = M[r] ^ M[row]
                row += 1

            free_vars = [c for c in range(512) if c not in pivots]

            # CE from basis
            ces = []
            for fc in free_vars[:min(256, len(free_vars))]:
                x = np.zeros(512, dtype=np.uint8)
                x[fc] = 1
                for i in range(len(pivots)-1, -1, -1):
                    pc = pivots[i]
                    val = np.uint8(0)
                    for j in range(512):
                        if j != pc:
                            val ^= (M[i,j] & x[j])
                    x[pc] = val
                dW = [0]*16
                for word in range(16):
                    for bit in range(32):
                        if x[word*32+bit]:
                            dW[word] ^= (1<<bit)
                W2 = [W_base[i]^dW[i] for i in range(16)]
                H2 = sha256_n_rounds(W2, n_r)
                ce = []
                for w in range(8):
                    delta = H_base[w] ^ H2[w]
                    for b in range(32):
                        ce.append((delta >> b) & 1)
                ces.append(ce)

            if ces:
                CE = np.array(ces[:256], dtype=np.uint8)
                rank_CE = np.linalg.matrix_rank(CE.astype(float))
                marker = " ★★★" if rank_CE < 200 else " ★★" if rank_CE < 240 else " ★" if rank_CE < 256 else ""
                print(f"  {n_r:3d} rounds: rank(T)={rank_T}, rank(CE)={rank_CE}{marker}")
            else:
                print(f"  {n_r:3d} rounds: rank(T)={rank_T}, CE=N/A (no kernel)")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("4. ГЕОМЕТРИЯ vs БЕЗОПАСНОСТЬ")
    print("=" * 70)

    print(f"""
  В нашем измерении:
    Геометрия (кривизна K) достигает максимума на раунде ~{transition or '?'}.
    Алгебра (rank CE) = 256 начиная с ~??? раундов.
    Безопасность (2^128) определяется rank(CE)=256 + K=128.

  Если rank(CE) < 256 для малого числа раундов:
    → REDUCED-ROUND SHA-256 имеет алгебраическую СЛАБОСТЬ!
    → GF(2)-kernel collision СУЩЕСТВУЕТ и вычислима!
    → Но только для reduced rounds, не для full 64.

  SHA-256 с 64 раундами: геометрически = идеальная сфера.
  С < {transition or '?'} раундами: геометрически = "мягкая" (K < 128).
  Мягкость = потенциальная уязвимость.
""")


if __name__ == "__main__":
    main()
