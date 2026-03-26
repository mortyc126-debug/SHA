#!/usr/bin/env python3
"""
SCF: Масштабирование defect — поиск W_base с максимальным defect.

Факты:
- Defect=1-3 найден при R=56 (kernel dim=256)
- 63% баз имеют defect≥1
- Вопрос: существуют ли базы с defect >> 3?

Стратегия:
1. Массовый поиск W_base (1000+)
2. Custom IV (multiblock) — другой carry landscape
3. Комбинация null-vector + Wang для amplification
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
def hw(x): return bin(x).count('1')

def expand_xor(W16):
    W=list(W16)
    for i in range(16,64):
        W.append(sig1(W[i-2])^W[i-7]^sig0(W[i-15])^W[i-16])
    return W

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

def bits_to_words(bits):
    words=[]
    for w in range(16):
        val=0
        for b in range(32):
            if bits[w*32+b]: val|=(1<<b)
        words.append(val)
    return words

def state_to_bits(state):
    bits=[]
    for w in state:
        for b in range(32): bits.append((w>>b)&1)
    return bits

def extract_kernel_basis(R_start):
    n_in=512; M=[]
    for ow in range(R_start,64):
        for ob in range(32):
            row=[0]*n_in
            for iw in range(16):
                for ib in range(32):
                    W=[0]*16; W[iw]=1<<ib
                    We=expand_xor(W)
                    if (We[ow]>>ob)&1: row[iw*32+ib]=1
            M.append(row)
    m=[list(r) for r in M]; nr,nc=len(m),n_in
    pivots=[]; ri=0
    for col in range(nc):
        pv=-1
        for r in range(ri,nr):
            if m[r][col]==1: pv=r; break
        if pv==-1: continue
        m[ri],m[pv]=m[pv],m[ri]
        for r in range(nr):
            if r!=ri and m[r][col]==1:
                for c in range(nc): m[r][c]^=m[ri][c]
        pivots.append(col); ri+=1
    free=[c for c in range(nc) if c not in pivots]
    basis=[]
    for fc in free:
        v=[0]*nc; v[fc]=1
        for i,pc in enumerate(pivots):
            if m[i][fc]==1: v[pc]=1
        basis.append(v)
    return basis

def fast_rank(rows):
    """Fast GF(2) rank of list of 256-bit vectors represented as ints."""
    m = list(rows)
    rank = 0
    for bit in range(255, -1, -1):
        mask = 1 << bit
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
    return rank

def state_to_int(state):
    """Convert 8×32-bit state to single 256-bit integer."""
    val = 0
    for i, w in enumerate(state):
        val |= w << (i * 32)
    return val


# ============================================================
# EXP 1: Mass search for high-defect W_base
# ============================================================
def exp1_mass_search(basis, N_bases):
    print("="*70)
    print(f"EXP 1: MASS SEARCH FOR HIGH DEFECT (N={N_bases})")
    print("="*70)

    n = len(basis)
    # Precompute basis words
    basis_words = [bits_to_words(v) for v in basis]

    from collections import Counter
    defect_dist = Counter()
    best_defect = 0
    best_W = None

    for bi in range(N_bases):
        W_base = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        H_base = sha_compress(W_base)
        H_int = state_to_int(H_base)

        # Compute δH as integers for fast rank
        dH_ints = []
        for bw in basis_words:
            W_mod = [(W_base[j]^bw[j]) for j in range(16)]
            H_mod = sha_compress(W_mod)
            dH_ints.append(state_to_int(H_base) ^ state_to_int(H_mod))

        rank = fast_rank(dH_ints)
        defect = n - rank
        defect_dist[defect] += 1

        if defect > best_defect:
            best_defect = defect
            best_W = list(W_base)
            print(f"  [{bi:5d}] NEW RECORD: defect = {defect}")

        if (bi+1) % (N_bases//5) == 0:
            print(f"  [{bi+1:5d}] scanned, best so far: {best_defect}")

    print(f"\n  Distribution: {dict(sorted(defect_dist.items()))}")
    print(f"  Best defect: {best_defect}")

    return best_defect, best_W, defect_dist


# ============================================================
# EXP 2: Custom IV search
# ============================================================
def exp2_custom_iv(basis, N_ivs):
    print("\n" + "="*70)
    print(f"EXP 2: CUSTOM IV SEARCH (N={N_ivs})")
    print("  Does different IV give higher defect?")
    print("="*70)

    n = len(basis)
    basis_words = [bits_to_words(v) for v in basis]

    from collections import Counter
    defect_dist = Counter()
    best_defect = 0

    for ii in range(N_ivs):
        custom_iv = [int.from_bytes(os.urandom(4),'big') for _ in range(8)]
        W_base = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        H_base = sha_compress(W_base, iv=custom_iv)

        dH_ints = []
        for bw in basis_words:
            W_mod = [(W_base[j]^bw[j]) for j in range(16)]
            H_mod = sha_compress(W_mod, iv=custom_iv)
            dH_ints.append(state_to_int(H_base) ^ state_to_int(H_mod))

        rank = fast_rank(dH_ints)
        defect = n - rank
        defect_dist[defect] += 1

        if defect > best_defect:
            best_defect = defect
            print(f"  [{ii:5d}] NEW RECORD: defect = {defect} (custom IV)")

    print(f"\n  Distribution: {dict(sorted(defect_dist.items()))}")
    print(f"  Best defect with custom IV: {best_defect}")

    return best_defect, defect_dist


# ============================================================
# EXP 3: Null vector → SA → best near-collision across many bases
# ============================================================
def exp3_best_near_collision(basis, N_bases, N_sa):
    print("\n" + "="*70)
    print(f"EXP 3: BEST NEAR-COLLISION VIA NULL+SA (N_bases={N_bases})")
    print("="*70)

    n = len(basis)
    basis_words = [bits_to_words(v) for v in basis]
    global_best = 256

    for bi in range(N_bases):
        W_base = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        H_base = sha_compress(W_base)

        def score_bits(bits):
            dW = bits_to_words(bits)
            W_mod = [(W_base[j]^dW[j]) for j in range(16)]
            H_mod = sha_compress(W_mod)
            return sum(hw(H_base[j]^H_mod[j]) for j in range(8))

        # Start from random kernel element
        cur_bits = list(basis[bi % n])
        cur_score = score_bits(cur_bits)
        best_local = cur_score

        for it in range(N_sa):
            trial = list(cur_bits)
            idx = int.from_bytes(os.urandom(2),'big') % n
            for j in range(512): trial[j] ^= basis[idx][j]
            if all(b==0 for b in trial): continue

            sc = score_bits(trial)
            T = max(0.01, 1.0 - it/N_sa)
            if sc < cur_score or \
               math.exp(-(sc-cur_score)/(T*2)) > int.from_bytes(os.urandom(4),'big')/(1<<32):
                cur_bits = trial
                cur_score = sc
            if cur_score < best_local:
                best_local = cur_score

        if best_local < global_best:
            global_best = best_local
            print(f"  [{bi:3d}] New best near-collision: HW = {global_best}")

    print(f"\n  GLOBAL BEST NEAR-COLLISION (kernel-constrained SA): HW = {global_best}")

    # Compare with unconstrained SA
    print(f"\n  Comparison: unconstrained SA")
    unc_best = 256
    for bi in range(N_bases):
        W_base = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        H_base = sha_compress(W_base)
        W2 = list(W_base)
        W2[0] ^= int.from_bytes(os.urandom(4),'big')
        cur_score = sum(hw(H_base[j]^sha_compress(W2)[j]) for j in range(8))

        for it in range(N_sa):
            trial_W2 = list(W2)
            w = int.from_bytes(os.urandom(1),'big') % 16
            b = int.from_bytes(os.urandom(1),'big') % 32
            trial_W2[w] ^= (1 << b)
            H2 = sha_compress(trial_W2)
            sc = sum(hw(H_base[j]^H2[j]) for j in range(8))
            T = max(0.01, 1.0 - it/N_sa)
            if sc < cur_score or \
               math.exp(-(sc-cur_score)/(T*2)) > int.from_bytes(os.urandom(4),'big')/(1<<32):
                W2 = trial_W2
                cur_score = sc
        if cur_score < unc_best:
            unc_best = cur_score

    print(f"  Unconstrained SA best: HW = {unc_best}")
    print(f"  Kernel advantage: {unc_best - global_best:+d} bits")


# ============================================================
if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 500

    print("="*70)
    print("SCF: DEFECT SCALING — массовый поиск")
    print("="*70)

    print("\nExtracting kernel R=56...")
    basis = extract_kernel_basis(56)
    print(f"  dim = {len(basis)}")

    best_d, best_W, dist1 = exp1_mass_search(basis, N_bases=N)
    best_d_iv, dist2 = exp2_custom_iv(basis, N_ivs=min(N, 300))
    exp3_best_near_collision(basis, N_bases=min(N//10+1, 30), N_sa=min(N*20, 5000))

    print("\n" + "="*70)
    print("ФИНАЛЬНЫЙ ИТОГ")
    print("="*70)
    print(f"""
  Standard IV: max defect = {best_d} (across {N} bases)
  Custom IV:   max defect = {best_d_iv} (across {min(N,300)} IVs)

  Если defect растёт с N → структурная слабость масштабируется
  Если defect стабилен ~3 → фундаментальный предел

  Defect = d даёт birthday O(2^{{(256-d)/2}})
  Defect 3 → O(2^{{126.5}}) — на 1.5 бита лучше standard
  Defect 10 → O(2^{{123}}) — на 5 бит лучше
  Defect 128 → O(2^{{64}}) — ПРОРЫВ
    """)
