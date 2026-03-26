#!/usr/bin/env python3
"""
SCF: АЛГЕБРАИЧЕСКАЯ ГЕОМЕТРИЯ SHA-256.

SHA(W) = H — это отображение f: Z_{2^32}^{16} → Z_{2^32}^{8}.
Коллизия: f(W1) = f(W2), W1 ≠ W2.
Эквивалентно: F(W1,W2) = f(W1) - f(W2) = 0, (W1,W2) ∉ diagonal.

F определяет МНОГООБРАЗИЕ V = {(W1,W2) : F=0} в Z^{32}_{2^32}.
dim(V) = 32 - 8 = 24 (ожидаемая). В битах: 24×32 = 768 бит.

Но V содержит ДИАГОНАЛЬ D = {(W,W)} как тривиальную компоненту.
dim(D) = 16 (слов) = 512 бит.

Коллизия = точка на V\D (вне диагонали).

Вопросы:
1. Какова РЕАЛЬНАЯ размерность V рядом с D? (≥ 24?)
2. Есть ли ДРУГИЕ компоненты V кроме D?
3. Как V пересекает гиперплоскости H[reg]=0?

Измеряем через: локальная размерность (rank Якобиана F вблизи D).
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


# ============================================================
# EXP 1: Dimension of collision variety near diagonal
# ============================================================
def exp1_local_dimension(N):
    """Near the diagonal (W1≈W2), compute rank of Jacobian of F(W1,W2)=SHA(W1)-SHA(W2).
    F has 8×32=256 components and 32×32=1024 variables (W1+W2).
    rank(J_F) = number of independent constraints.
    Local dimension of V = 1024 - rank(J_F)."""

    print("="*70)
    print("EXP 1: LOCAL DIMENSION OF COLLISION VARIETY")
    print("  F(W1,W2) = SHA(W1) - SHA(W2)")
    print("  J_F ∈ M(256 × 1024) over GF(2)")
    print("  dim(V) = 1024 - rank(J_F)")
    print("="*70)

    for trial in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        # W2 near W1 (small perturbation)
        W2 = list(W1); W2[0] ^= 1

        H1 = sha_compress(W1); H2 = sha_compress(W2)

        # Jacobian: for each input bit (512 for W1 + 512 for W2),
        # measure which output bits of F change.
        # F = H1 ⊕ H2 (XOR for GF(2) approximation)

        rows = []  # Each row = 1024-bit vector (effect on 256 output bits)

        # Perturbations of W1
        for word in range(16):
            for bit in range(32):
                W1t = list(W1); W1t[word] ^= (1<<bit)
                H1t = sha_compress(W1t)
                # ΔF = (H1t ⊕ H2) ⊕ (H1 ⊕ H2) = H1t ⊕ H1
                dF = [H1t[i] ^ H1[i] for i in range(8)]
                row = 0
                for i in range(8):
                    row |= dF[i] << (i*32)
                rows.append(row)

        # Perturbations of W2
        for word in range(16):
            for bit in range(32):
                W2t = list(W2); W2t[word] ^= (1<<bit)
                H2t = sha_compress(W2t)
                dF = [H2t[i] ^ H2[i] for i in range(8)]
                row = 0
                for i in range(8):
                    row |= dF[i] << (i*32)
                rows.append(row)

        # Rank computation over GF(2) using integer representation
        m = list(rows)
        rank = 0
        for bit_pos in range(255, -1, -1):
            mask = 1 << bit_pos
            pivot = -1
            for i in range(rank, len(m)):
                if m[i] & mask:
                    pivot = i; break
            if pivot == -1: continue
            m[rank], m[pivot] = m[pivot], m[rank]
            for i in range(len(m)):
                if i != rank and m[i] & mask:
                    m[i] ^= m[rank]
            rank += 1

        local_dim = 1024 - rank
        print(f"  Trial {trial}: rank(J_F) = {rank}, local dim(V) = {local_dim}")
        print(f"    Expected: rank=256 (full) → dim=768")
        print(f"    Diagonal dim: 512")
        print(f"    Extra dim: {local_dim - 512}")

    return rank


# ============================================================
# EXP 2: Homotopy — continuous path from trivial to nontrivial
# ============================================================
def exp2_homotopy(N):
    """Idea: parameterize W2(t) = W1 + t·δW for t ∈ {0, 1, ..., 2^32-1}.
    At t=0: W2=W1 (trivial collision).
    At t=1: W2=W1+δW (near-collision HW≈128).

    Question: is there a SMOOTH path in t where HW(F(t)) stays low?
    If yes → follow the path from t=0 to large t with HW=0."""

    print("\n" + "="*70)
    print("EXP 2: HOMOTOPY PATH — track HW along t·δW")
    print("="*70)

    W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    H1 = sha_compress(W1)

    # Direction δW
    dW = [0]*16; dW[0] = 1  # Simplest: add 1 to W[0]

    print(f"  δW = [1, 0, 0, ..., 0] (add 1 to W[0])")
    print(f"  {'t':>10} | HW(δH) | notes")
    print("  " + "-"*40)

    prev_hw = 0
    for exp in range(33):
        t = 1 << exp if exp < 32 else MASK32
        W2 = list(W1)
        W2[0] = add32(W1[0], t)
        H2 = sha_compress(W2)
        d = sum(hw(H1[i]^H2[i]) for i in range(8))

        marker = ""
        if d == 0: marker = " ★ ZERO!"
        elif d < prev_hw: marker = " ↓"
        prev_hw = d

        if exp < 10 or exp % 4 == 0 or d < 100:
            print(f"  {t:10d} | {d:6d} | 2^{exp}{marker}")

    # Try many small t values
    print(f"\n  Fine scan t=1..1000:")
    min_hw = 256; min_t = 0
    for t in range(1, 1001):
        W2 = list(W1); W2[0] = add32(W1[0], t)
        H2 = sha_compress(W2)
        d = sum(hw(H1[i]^H2[i]) for i in range(8))
        if d < min_hw:
            min_hw = d; min_t = t

    print(f"  Best: t={min_t}, HW={min_hw}")

    # Multiple directions
    print(f"\n  Multiple directions (t=1):")
    for word in range(16):
        W2 = list(W1); W2[word] = add32(W1[word], 1)
        H2 = sha_compress(W2)
        d = sum(hw(H1[i]^H2[i]) for i in range(8))
        print(f"    δW[{word:2d}]=+1: HW={d}")


# ============================================================
# EXP 3: Variety intersection — how many H[reg]=0 can we have?
# ============================================================
def exp3_reg_zeroing(N):
    """For random pairs, count how many registers have H1[i]=H2[i].
    This measures the intersection of V with hyperplanes."""

    print("\n" + "="*70)
    print("EXP 3: REGISTER ZEROING — how many regs can match?")
    print("="*70)

    max_zeros = [0] * (N)

    for trial in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        H1 = sha_compress(W1)

        best = 0
        # SA to maximize matching registers
        W2 = list(W1); W2[0] ^= 0x80000000
        cur = list(W2)

        for it in range(2000):
            t = list(cur)
            w = int.from_bytes(os.urandom(1),'big') % 16
            b = int.from_bytes(os.urandom(1),'big') % 32
            t[w] ^= (1<<b)
            if t == W1: continue

            H2 = sha_compress(t)
            n_zero_regs = sum(1 for i in range(8) if H1[i] == H2[i])
            n_total_hw = sum(hw(H1[i]^H2[i]) for i in range(8))

            # Optimize for max zero registers (secondary: min total HW)
            cur_H2 = sha_compress(cur)
            cur_zeros = sum(1 for i in range(8) if H1[i] == cur_H2[i])
            cur_hw = sum(hw(H1[i]^cur_H2[i]) for i in range(8))

            score_new = n_zero_regs * 1000 - n_total_hw
            score_cur = cur_zeros * 1000 - cur_hw

            T = max(0.01, 1-it/2000)
            if score_new > score_cur or \
               math.exp((score_new-score_cur)/(T*50)) > int.from_bytes(os.urandom(4),'big')/(1<<32):
                cur = t
                if n_zero_regs > best:
                    best = n_zero_regs

        max_zeros[trial] = best
        if best >= 1:
            # Show which regs are zero
            H2 = sha_compress(cur)
            zero_regs = [i for i in range(8) if H1[i] == H2[i]]
            per_reg = [hw(H1[i]^H2[i]) for i in range(8)]
            print(f"  Trial {trial}: {best} zero regs {zero_regs}, per_reg={per_reg}")

    from collections import Counter
    dist = Counter(max_zeros)
    print(f"\n  Distribution of max zero registers:")
    for k in sorted(dist): print(f"    {k} regs: {dist[k]} trials")
    print(f"  Expected random: P(k regs match) = C(8,k) × 2^{{-32k}}")


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    exp1_local_dimension(min(N, 3))
    exp2_homotopy(1)
    exp3_reg_zeroing(min(N*2, 10))
