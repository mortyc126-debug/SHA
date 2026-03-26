#!/usr/bin/env python3
"""
SCF: GLOBAL DIRECTION FINDER — ΔW что улучшает ВСЕ 8 регистров.

Проблема: плотный Якобиан → коррекция одного reg портит другие.
Решение: найти ΔW где J·ΔW < 0 по ВСЕМ 8 компонентам.

Метод: НЕ случайный поиск. Используем структуру J.
J — матрица 8×16 над Z_{2^32}. Найти ΔW ∈ Z_{2^32}^{16}
такой что J·ΔW ≈ −dH (вектор текущих ошибок).

Это ЗАДАЧА НАИМЕНЬШИХ КВАДРАТОВ над Z_{2^32}:
  minimize ||J·ΔW + dH||^2

Но Z_{2^32} не поле → нет стандартного LS.
Наш инструмент: решаем ПОБИТОВО, от LSB к MSB (Hensel-style).
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

def hash_hw_pair(W1, W2):
    H1=sha_compress(W1); H2=sha_compress(W2)
    return sum(hw(H1[i]^H2[i]) for i in range(8))

def hash_per_reg(W1, W2):
    H1=sha_compress(W1); H2=sha_compress(W2)
    return [hw(H1[i]^H2[i]) for i in range(8)]


# ============================================================
# Tool 1: MULTI-WORD CORRECTION — solve for ΔW as a vector
# ============================================================
def multi_word_correct(W1, W2):
    """Find additive correction ΔW[0..15] that improves ALL registers.

    For each register reg: compute sensitivity to each word.
    Then find ΔW that makes the SIGNED improvement positive for all regs.

    Key insight: we solve 8 equations simultaneously,
    using 16 words as unknowns. System is UNDERDETERMINED (8<16).
    The 8 extra degrees of freedom let us satisfy ALL constraints."""

    H1 = sha_compress(W1); H2 = sha_compress(W2)
    dH = [sub32(H2[i], H1[i]) for i in range(8)]

    # Compute full 8×16 Jacobian
    J = [[0]*16 for _ in range(8)]
    for w in range(16):
        W2t = list(W2); W2t[w] = add32(W2[w], 1)
        H2t = sha_compress(W2t)
        for reg in range(8):
            J[reg][w] = sub32(sub32(H2t[reg], H1[reg]), dH[reg])

    # Strategy: solve for 8 registers using 8 words (pick best 8 words)
    # Use words with most odd (invertible) sensitivities
    word_scores = []
    for w in range(16):
        n_odd = sum(1 for reg in range(8) if J[reg][w] & 1)
        word_scores.append((n_odd, w))
    word_scores.sort(reverse=True)
    active_words = [w for _, w in word_scores[:8]]

    # Build 8×8 subsystem: J_sub · ΔW_sub = -dH
    J_sub = [[J[reg][w] for w in active_words] for reg in range(8)]

    # Solve mod 2 first (GF(2))
    # Then lift to mod 4, 8, ..., 2^32 (Hensel)

    # GF(2) solve
    mat_gf2 = [[(J_sub[reg][j] & 1) for j in range(8)] + [(dH[reg] & 1)] for reg in range(8)]

    # Gaussian elimination mod 2
    sol_gf2 = [0]*8
    m = [list(r) for r in mat_gf2]
    pivots = []; ri = 0
    for col in range(8):
        pv = -1
        for r in range(ri, 8):
            if m[r][col]: pv=r; break
        if pv == -1: continue
        m[ri],m[pv] = m[pv],m[ri]
        for r in range(8):
            if r!=ri and m[r][col]:
                for c in range(9): m[r][c] ^= m[ri][c]
        pivots.append(col); ri+=1

    rank = len(pivots)
    consistent = all(m[r][8]==0 for r in range(rank, 8))

    if consistent:
        for i, pc in enumerate(pivots):
            sol_gf2[pc] = m[i][8]

    # Apply GF(2) correction
    corrections = [0]*16
    for j, w in enumerate(active_words):
        corrections[w] = sol_gf2[j]

    W2_corrected = list(W2)
    for w in range(16):
        if corrections[w]:
            W2_corrected[w] = add32(W2[w], corrections[w])

    if W2_corrected == W1:
        return W2, hash_hw_pair(W1, W2), False

    new_hw = hash_hw_pair(W1, W2_corrected)
    old_hw = hash_hw_pair(W1, W2)

    return W2_corrected, new_hw, new_hw < old_hw


# ============================================================
# Tool 2: ITERATIVE MULTI-WORD with recomputation
# ============================================================
def iterative_multi_word(W1, W2_init, max_iter=30):
    W2 = list(W2_init)
    best_hw = hash_hw_pair(W1, W2)
    best_W2 = list(W2)

    for it in range(max_iter):
        W2_new, new_hw, improved = multi_word_correct(W1, W2)

        if improved and W2_new != W1:
            W2 = W2_new
            if new_hw < best_hw:
                best_hw = new_hw
                best_W2 = list(W2)
        else:
            # Perturbation to escape
            w = it % 16
            W2[w] ^= (1 << (it % 32))
            if W2 == W1: W2[w] ^= 2

    return best_W2, best_hw


# ============================================================
# Tool 3: COMBINED — multi-word + SA from result
# ============================================================
def combined_solver(W1, delta_w0=0x80000000, budget=3000):
    W2 = list(W1); W2[0] ^= delta_w0
    H1 = sha_compress(W1)

    # Phase 1: Screen δ (128 evals)
    screen = []
    for word in range(4):
        for bit in range(32):
            W2t = list(W1); W2t[word] ^= (1<<bit)
            screen.append((hash_hw_pair(W1, W2t), word, bit))
    screen.sort()

    global_best = 256; global_W2 = None

    # Phase 2: For top-3 δ, run multi-word + SA
    per_start = (budget - 128) // 3

    for _, dw, db in screen[:3]:
        W2 = list(W1); W2[dw] ^= (1<<db)

        # Multi-word iterative (50 evals)
        W2, mw_hw = iterative_multi_word(W1, W2, max_iter=10)

        # SA polish (remaining budget)
        cur = list(W2); best = mw_hw
        for it in range(per_start - 20):
            t = list(cur)
            w = int.from_bytes(os.urandom(1),'big') % 16
            b = int.from_bytes(os.urandom(1),'big') % 32
            t[w] ^= (1<<b)
            if t == W1: continue
            s = hash_hw_pair(W1, t)
            T = max(0.01, 1-it/max(per_start-20,1))
            if s < best or math.exp(-(s-best)/(T*2)) > int.from_bytes(os.urandom(4),'big')/(1<<32):
                cur = t
                if s < best: best = s

        if best < global_best:
            global_best = best
            global_W2 = list(cur)

    return global_W2, global_best


# ============================================================
if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10

    print("="*70)
    print("SCF: GLOBAL DIRECTION FINDER")
    print("="*70)

    combined_results = []
    sa_results = []

    for trial in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

        W2_c, hw_c = combined_solver(W1, budget=3000)
        combined_results.append(hw_c)

        # SA baseline
        W2_sa = list(W1); W2_sa[0] ^= 0x80000000
        best_sa = hash_hw_pair(W1, W2_sa); cur = list(W2_sa)
        for it in range(3000):
            t = list(cur)
            w = int.from_bytes(os.urandom(1),'big') % 16
            b = int.from_bytes(os.urandom(1),'big') % 32
            t[w] ^= (1<<b)
            if t == W1: continue
            s = hash_hw_pair(W1, t)
            T = max(0.01,1-it/3000)
            if s<best_sa or math.exp(-(s-best_sa)/(T*2))>int.from_bytes(os.urandom(4),'big')/(1<<32):
                cur=t
                if s<best_sa: best_sa=s
        sa_results.append(best_sa)

        winner = "OURS" if hw_c < best_sa else ("SA" if best_sa < hw_c else "tie")
        print(f"  Trial {trial:2d}: Ours={hw_c:3d}  SA={best_sa:3d}  → {winner}")

    avg_c = sum(combined_results)/N
    avg_sa = sum(sa_results)/N
    wins = sum(1 for i in range(N) if combined_results[i] < sa_results[i])

    print(f"\n{'='*70}")
    print(f"  Ours:  avg={avg_c:.1f}  min={min(combined_results)}")
    print(f"  SA:    avg={avg_sa:.1f}  min={min(sa_results)}")
    print(f"  Wins:  Ours={wins}  SA={N-wins}")
    print(f"  Advantage: {avg_sa-avg_c:+.1f}")
