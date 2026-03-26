#!/usr/bin/env python3
"""
SCF: NEWTON В ПРОСТРАНСТВЕ ЯДРА — детерминированная атака.

Факт: SA/random даёт HW(δe)=5-7. Это случайный поиск в 2^128.
Факт: Wang chain использует Newton (sensitivity=1) и находит ТОЧНЫЙ ноль.

Идея: вычислить Якобиан J = d(δe[r])/d(kernel_coeff).
J — матрица 32×128 над GF(2) (для одного раунда).
Решить: J·c = target → c = kernel coefficients.

Если rank(J) = 32 → система решаема → δe[r]=0 ТОЧНО!
Если rank(J) < 32 → ядро J → дополнительная свобода.
"""
import os, sys

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

def expand_real(W16):
    W=list(W16)
    for i in range(16,64):
        W.append(add32(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    return W

def expand_xor(W16):
    W=list(W16)
    for i in range(16,64):
        W.append(sig1(W[i-2])^W[i-7]^sig0(W[i-15])^W[i-16])
    return W

def sha_compress_full(W16):
    W=expand_real(W16); s=list(IV)
    for r in range(64):
        a,b,c,d,e,f,g,h=s
        T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add32(Sig0(a),Maj(a,b,c))
        s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
    return s

def sha_state_at_r(W16, R):
    W=expand_real(W16); s=list(IV)
    for r in range(R):
        a,b,c,d,e,f,g,h=s
        T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add32(Sig0(a),Maj(a,b,c))
        s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
    return s

def bits_to_words(bits):
    words=[]
    for w in range(16):
        val=0
        for b in range(32):
            if bits[w*32+b]: val|=(1<<b)
        words.append(val)
    return words

def extract_kernel_basis(R_start=52):
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

def gf2_solve(M, target):
    """Solve M·x = target over GF(2). M is n×m, target is n-bit.
    Returns (solution, rank) or (None, rank) if no solution."""
    n = len(M)
    m = len(M[0]) if M else 0
    # Augmented matrix [M | target]
    aug = [list(M[i]) + [(target>>i)&1] for i in range(n)]

    pivots = []
    ri = 0
    for col in range(m):
        pv = -1
        for r in range(ri, n):
            if aug[r][col]: pv=r; break
        if pv == -1: continue
        aug[ri], aug[pv] = aug[pv], aug[ri]
        for r in range(n):
            if r != ri and aug[r][col]:
                for c in range(m+1): aug[r][c] ^= aug[ri][c]
        pivots.append(col)
        ri += 1

    rank = len(pivots)

    # Check consistency
    for r in range(rank, n):
        if aug[r][m]:
            return None, rank  # Inconsistent

    # Extract solution (set free vars to 0)
    sol = [0] * m
    for i, pc in enumerate(pivots):
        sol[pc] = aug[i][m]

    return sol, rank


# ============================================================
# EXP 1: Jacobian of kernel → δe[r] map
# ============================================================
def exp1_jacobian(basis, W_base, target_round):
    """Compute GF(2) Jacobian: d(δe[r]) / d(kernel_coeff_i).

    For each basis vector i: flip it, measure change in δe[target_round].
    J[bit][i] = 1 if flipping basis[i] flips bit of δe.
    """
    n_basis = len(basis)
    basis_words = [bits_to_words(v) for v in basis]

    # Compute δe at target_round for "zero kernel" (δW = basis[0])
    ref_bits = list(basis[0])
    ref_words = bits_to_words(ref_bits)
    W_ref = [(W_base[j] ^ ref_words[j]) for j in range(16)]

    state_base = sha_state_at_r(W_base, target_round + 1)
    state_ref = sha_state_at_r(W_ref, target_round + 1)
    de_ref = state_base[4] ^ state_ref[4]  # δe at target_round+1

    # Jacobian: for each basis vector, XOR it with ref and measure δe change
    J = []  # 32 rows (bits of δe) × n_basis cols
    for bit in range(32):
        row = []
        for i in range(n_basis):
            # Flip basis[i]
            test_bits = list(ref_bits)
            for j in range(512):
                test_bits[j] ^= basis[i][j]

            test_words = bits_to_words(test_bits)
            W_test = [(W_base[j] ^ test_words[j]) for j in range(16)]
            state_test = sha_state_at_r(W_test, target_round + 1)
            de_test = state_base[4] ^ state_test[4]

            # Change in bit of δe
            row.append(((de_ref >> bit) ^ (de_test >> bit)) & 1)
        J.append(row)

    return J, de_ref


def exp1_run(basis):
    print("="*70)
    print("EXP 1: JACOBIAN RANK — kernel → δe[r]")
    print("="*70)

    W_base = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    n_basis = len(basis)

    print(f"  Kernel dim: {n_basis}")
    print(f"  {'r':>3} | {'rank(J)':>8} | {'solvable':>9} | notes")
    print("  " + "-"*45)

    for target_r in range(17, 52, 2):
        J, de_ref = exp1_jacobian(basis, W_base, target_r)

        # Rank
        from copy import deepcopy
        m = deepcopy(J)
        rank = 0
        for col in range(n_basis):
            pv = -1
            for r in range(rank, 32):
                if m[r][col]: pv=r; break
            if pv==-1: continue
            m[rank],m[pv] = m[pv],m[rank]
            for r in range(32):
                if r!=rank and m[r][col]:
                    for c in range(n_basis): m[r][c]^=m[rank][c]
            rank+=1

        # Can we solve J·c = de_ref (to zero δe)?
        sol, _ = gf2_solve(J, de_ref)
        solvable = "YES" if sol is not None else "no"

        marker = " ★★★" if sol is not None else (" ★" if rank >= 30 else "")
        print(f"  {target_r:3d} | {rank:8d} | {solvable:>9} |{marker}")

        if sol is not None:
            # Verify!
            combined = list(basis[0])  # Start with ref
            for i in range(n_basis):
                if sol[i]:
                    for j in range(512):
                        combined[j] ^= basis[i][j]

            dW = bits_to_words(combined)
            W_test = [(W_base[j]^dW[j]) for j in range(16)]
            state_test = sha_state_at_r(W_test, target_r + 1)
            state_base_r = sha_state_at_r(W_base, target_r + 1)
            de_actual = state_base_r[4] ^ state_test[4]

            print(f"        Verification: HW(δe) = {hw(de_actual)}")
            if de_actual == 0:
                print(f"        ★★★ EXACT ZERO! Newton solved it!")

    return W_base


# ============================================================
# EXP 2: Newton iteration — correct for nonlinearity
# ============================================================
def exp2_newton(basis, W_base, target_r, max_iter=20):
    """
    Newton's method in kernel space:
    1. Start with current kernel coefficient vector c
    2. Compute δe[r] = f(c)
    3. Compute Jacobian J = df/dc
    4. Solve J·Δc = δe[r]
    5. Update c ← c ⊕ Δc
    6. Repeat until δe[r] = 0
    """
    print(f"\n  Newton iteration for r={target_r}:")

    n_basis = len(basis)

    # Start: c = [0, 0, ..., 0] (use first basis vector as δW)
    current_c = [0] * n_basis
    current_c[0] = 1  # Start with basis[0]

    for iteration in range(max_iter):
        # Build current δW from coefficients
        current_bits = [0] * 512
        for i in range(n_basis):
            if current_c[i]:
                for j in range(512):
                    current_bits[j] ^= basis[i][j]

        if all(b==0 for b in current_bits):
            current_c[1] = 1  # Avoid zero
            continue

        dW = bits_to_words(current_bits)
        W_mod = [(W_base[j]^dW[j]) for j in range(16)]

        state_base = sha_state_at_r(W_base, target_r + 1)
        state_mod = sha_state_at_r(W_mod, target_r + 1)
        de = state_base[4] ^ state_mod[4]

        if de == 0:
            print(f"    Iter {iteration}: ★★★ CONVERGED! δe = 0")
            return current_c, True

        print(f"    Iter {iteration}: HW(δe) = {hw(de)}")

        # Compute Jacobian at current point
        J = []
        for bit in range(32):
            row = []
            for i in range(n_basis):
                # Flip coefficient i
                test_c = list(current_c)
                test_c[i] ^= 1

                test_bits = [0]*512
                for k in range(n_basis):
                    if test_c[k]:
                        for j in range(512):
                            test_bits[j] ^= basis[k][j]
                if all(b==0 for b in test_bits):
                    row.append(0); continue

                test_dW = bits_to_words(test_bits)
                W_t = [(W_base[j]^test_dW[j]) for j in range(16)]
                st_t = sha_state_at_r(W_t, target_r + 1)
                de_t = state_base[4] ^ st_t[4]

                row.append(((de >> bit) ^ (de_t >> bit)) & 1)
            J.append(row)

        # Solve J·Δc = de
        sol, rank = gf2_solve(J, de)

        if sol is None:
            print(f"    Iter {iteration}: J rank={rank}, system inconsistent")
            # Try random perturbation
            idx = int.from_bytes(os.urandom(2),'big') % n_basis
            current_c[idx] ^= 1
            continue

        # Update c
        for i in range(n_basis):
            current_c[i] ^= sol[i]

    return current_c, False


# ============================================================
# EXP 3: Newton for HASH (not individual round)
# ============================================================
def exp3_newton_hash(basis, max_iter=30):
    print("\n" + "="*70)
    print("EXP 3: NEWTON FOR HASH — minimize HW(δH) directly")
    print("="*70)

    n_basis = len(basis)
    W_base = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

    current_c = [0]*n_basis
    current_c[0] = 1

    for iteration in range(max_iter):
        current_bits = [0]*512
        for i in range(n_basis):
            if current_c[i]:
                for j in range(512): current_bits[j] ^= basis[i][j]
        if all(b==0 for b in current_bits):
            current_c[iteration % n_basis] = 1; continue

        dW = bits_to_words(current_bits)
        W_mod = [(W_base[j]^dW[j]) for j in range(16)]
        H_base = sha_compress_full(W_base)
        H_mod = sha_compress_full(W_mod)
        dH = [H_base[i]^H_mod[i] for i in range(8)]
        dH_hw = sum(hw(h) for h in dH)

        if dH_hw == 0:
            print(f"  Iter {iteration}: ★★★ COLLISION FOUND!")
            return current_c, True

        # Pick worst register, solve for it
        worst_reg = max(range(8), key=lambda i: hw(dH[i]))
        target = dH[worst_reg]

        # Jacobian for this register
        J = []
        for bit in range(32):
            row = []
            for i in range(n_basis):
                test_c = list(current_c); test_c[i] ^= 1
                tb = [0]*512
                for k in range(n_basis):
                    if test_c[k]:
                        for j in range(512): tb[j] ^= basis[k][j]
                if all(b==0 for b in tb):
                    row.append(0); continue
                tdW = bits_to_words(tb)
                Wt = [(W_base[j]^tdW[j]) for j in range(16)]
                Ht = sha_compress_full(Wt)
                dt = H_base[worst_reg] ^ Ht[worst_reg]
                row.append(((target>>bit)^(dt>>bit))&1)
            J.append(row)

        sol, rank = gf2_solve(J, target)

        if sol is not None:
            for i in range(n_basis): current_c[i] ^= sol[i]
            print(f"  Iter {iteration}: HW(δH)={dH_hw}, target reg={worst_reg}, J rank={rank}, Newton step applied")
        else:
            # Random perturbation fallback
            idx = int.from_bytes(os.urandom(2),'big') % n_basis
            current_c[idx] ^= 1
            print(f"  Iter {iteration}: HW(δH)={dH_hw}, J rank={rank}, random perturb")

    # Final score
    current_bits = [0]*512
    for i in range(n_basis):
        if current_c[i]:
            for j in range(512): current_bits[j] ^= basis[i][j]
    dW = bits_to_words(current_bits)
    W_mod = [(W_base[j]^dW[j]) for j in range(16)]
    H_base = sha_compress_full(W_base)
    H_mod = sha_compress_full(W_mod)
    final_hw = sum(hw(H_base[i]^H_mod[i]) for i in range(8))
    print(f"\n  Final HW(δH) = {final_hw}")

    return current_c, False


if __name__ == '__main__':
    print("="*70)
    print("SCF: NEWTON В ПРОСТРАНСТВЕ ЯДРА")
    print("="*70)

    print("\nExtracting kernel R=52...")
    basis = extract_kernel_basis(52)
    print(f"  dim = {len(basis)}")

    W_base = exp1_run(basis)

    # Try Newton on a few gap rounds
    for tr in [20, 30, 40]:
        c, success = exp2_newton(basis, W_base, tr, max_iter=10)
        if success:
            print(f"  ★ Newton CONVERGED at r={tr}!")

    exp3_newton_hash(basis, max_iter=20)
