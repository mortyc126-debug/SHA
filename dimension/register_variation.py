"""
ВАРИАЦИЯ РЕГИСТРОВ: 4 vs 8 регистров.

SHA-256 standard: 8 регистров (a,b,c,d,e,f,g,h) → 256 бит выход.
"SHA-128" (toy): 4 регистра (a,b,c,d) → 128 бит выход.

В нашем измерении:
  8 regs: 16 вход → 8 выход → C^4 = 2^128
  4 regs: 16 вход → 4 выход → C^2 = 2^64 ← ДЕШЕВЛЕ!

Вопросы:
  1. rank(CE) для 4-reg: 128 (full) или меньше?
  2. Security boundary для 4-reg: r=???
  3. Кривизна K: та же?
  4. Формула rank(CE) = f(n_regs, n_rounds)?
"""

import numpy as np

MASK32 = 0xFFFFFFFF
K_const = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
]

def hw(x): return bin(x).count('1')
def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def add32(x, y): return (x + y) & MASK32
def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)

IV_4 = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a]
IV_8 = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
        0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]

def sha_4reg(W16, n_rounds=16):
    """4-register variant: only (a,b,c,d), single node."""
    a,b,c,d = IV_4
    for r in range(n_rounds):
        # T1 = d + Σ₁(a) + K + W (simplified, no Ch)
        T1 = add32(add32(add32(d, Sigma1(a)), K_const[r%16]), W16[r%16])
        T2 = add32(Sigma0(a), Maj(a,b,c))
        d,c,b = c,b,a
        a = add32(T1, T2)
    return tuple(add32(IV_4[i], [a,b,c,d][i]) for i in range(4))

def sha_8reg(W16, n_rounds=16):
    """Standard 8-register SHA-256."""
    a,b,c,d,e,f,g,h = IV_8
    for r in range(n_rounds):
        T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K_const[r%16]), W16[r%16])
        T2 = add32(Sigma0(a), Maj(a,b,c))
        h,g,f,e = g,f,e,add32(d,T1)
        d,c,b,a = c,b,a,add32(T1,T2)
    return tuple(add32(IV_8[i], [a,b,c,d,e,f,g,h][i]) for i in range(8))


def compute_ce_rank(hash_fn, W_base, n_output):
    H_base = hash_fn(W_base)
    n_out_bits = n_output * 32
    T_rows = []
    for word in range(16):
        for bit in range(32):
            W_mod = list(W_base)
            W_mod[word] ^= (1 << bit)
            H_mod = hash_fn(W_mod)
            row = []
            for w in range(n_output):
                d = H_base[w] ^ H_mod[w]
                for b in range(32): row.append((d >> b) & 1)
            T_rows.append(row)

    T = np.array(T_rows, dtype=np.uint8)  # 512 × n_out_bits
    rank_T = np.linalg.matrix_rank(T.astype(float))

    if rank_T < n_out_bits:
        return rank_T, None

    M = T.T.copy()
    pivots = []; row = 0
    for col in range(512):
        found = False
        for r in range(row, n_out_bits):
            if M[r,col]==1: M[[row,r]]=M[[r,row]]; found=True; break
        if not found: continue
        pivots.append(col)
        for r in range(n_out_bits):
            if r!=row and M[r,col]==1: M[r]^=M[row]
        row += 1

    free_vars = [c for c in range(512) if c not in pivots]
    ces = []
    for fc in free_vars[:min(n_out_bits, len(free_vars))]:
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
        for w in range(n_output):
            d = H_base[w]^H2[w]
            for b in range(32): ce.append((d>>b)&1)
        ces.append(ce)

    if len(ces) < n_out_bits:
        CE = np.array(ces, dtype=np.uint8)
    else:
        CE = np.array(ces[:n_out_bits], dtype=np.uint8)
    rank_CE = np.linalg.matrix_rank(CE.astype(float))
    return rank_T, rank_CE


def main():
    np.random.seed(42)
    W_base = [np.random.randint(0, 2**32) for _ in range(16)]

    print("=" * 70)
    print("ВАРИАЦИЯ РЕГИСТРОВ: 4-reg vs 8-reg")
    print("=" * 70)

    # =================================================================
    print(f"\n{'=' * 70}")
    print("1. rank(CE) по раундам: 4-reg vs 8-reg")
    print("=" * 70)

    print(f"\n  {'Rounds':>6} {'4-reg T':>8} {'4-reg CE':>9} {'8-reg T':>8} {'8-reg CE':>9}")

    for n_r in [4, 6, 8, 10, 12, 14, 16]:
        rT4, rCE4 = compute_ce_rank(lambda W: sha_4reg(W, n_r), W_base, 4)
        rT8, rCE8 = compute_ce_rank(lambda W: sha_8reg(W, n_r), W_base, 8)

        print(f"  {n_r:6d} {rT4:>7} {str(rCE4) if rCE4 is not None else '—':>9} "
              f"{rT8:>7} {str(rCE8) if rCE8 is not None else '—':>9}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("2. SECURITY BOUNDARY: когда rank(CE) = full?")
    print("=" * 70)

    for name, hash_fn, n_out in [("4-reg", lambda W, nr=16: sha_4reg(W, nr), 4),
                                   ("8-reg", lambda W, nr=16: sha_8reg(W, nr), 8)]:
        full_rank = n_out * 32
        boundary = None
        for n_r in range(4, 20):
            fn = lambda W, _nr=n_r: sha_4reg(W, _nr) if n_out == 4 else sha_8reg(W, _nr)
            rT, rCE = compute_ce_rank(fn, W_base, n_out)
            if rCE == full_rank:
                boundary = n_r
                break

        print(f"  {name}: full rank={full_rank}, secure at r={boundary}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("3. COLLISION НА 4-REG (r < boundary)")
    print("=" * 70)

    # 4-reg с r=12: rank(CE) < 128?
    for n_r in [8, 10, 12, 14, 15, 16]:
        fn = lambda W, _nr=n_r: sha_4reg(W, _nr)
        H_base = fn(W_base)
        rT, rCE = compute_ce_rank(fn, W_base, 4)

        if rCE is not None and rCE < 128:
            # Find collision
            T_rows = []
            for word in range(16):
                for bit in range(32):
                    W_mod = list(W_base)
                    W_mod[word] ^= (1<<bit)
                    H_mod = fn(W_mod)
                    row = []
                    for w in range(4):
                        d = H_base[w]^H_mod[w]
                        for b in range(32): row.append((d>>b)&1)
                    T_rows.append(row)

            T = np.array(T_rows, dtype=np.uint8)
            # Quick collision check: flip each bit
            colls = 0
            for word in range(16):
                for bit in range(32):
                    W2 = list(W_base); W2[word]^=(1<<bit)
                    if fn(W2) == H_base:
                        colls += 1

            print(f"  4-reg r={n_r}: rank(CE)={rCE}, ker={128-rCE if rCE else '?'}, 1-bit collisions={colls}")
        else:
            print(f"  4-reg r={n_r}: rank(CE)={rCE} {'← SECURE' if rCE==128 else ''}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("4. ФОРМУЛА В НАШЕМ ИЗМЕРЕНИИ")
    print("=" * 70)

    print(f"""
  ОБОБЩЁННАЯ ФОРМУЛА:

  SHA с N_reg регистрами, C = 2^32 per register:
    Выход: N_reg позиций × C values = C^N_reg possibilities
    Вход: 16 позиций × C = C^16
    Избыток: C^16 / C^N_reg = C^(16-N_reg)

  Collision cost (overlay):
    C^(N_reg/2)

  rank(CE) для r ≥ N_reg: N_reg × 32 = full rank
  Security boundary: r = max(N_reg, 16) (all words must be active)

  Для SHA-256 (N_reg=8): C^4 = 2^128
  Для SHA-128 (N_reg=4): C^2 = 2^64
  Для SHA-512 (N_reg=8, C=2^64): (2^64)^4 = 2^256

  В нашем измерении: стоимость = C^(N_reg/2).
  Определяется ДВУМЯ числами: ёмкость C и число регистров N_reg.
  ВСЁ ОСТАЛЬНОЕ (rounds, schedule, Ch, Maj) = defense-in-depth.
""")


if __name__ == "__main__":
    main()
