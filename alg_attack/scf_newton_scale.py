#!/usr/bin/env python3
"""
SCF: МАСШТАБИРОВАНИЕ NEWTON — 3 пути к 576 dims.

Факт: Newton решает КАЖДЫЙ раунд gap при rank(J)=32.
Факт: нужно 18 раундов × 32 бит = 576 dims, есть 128.
Зазор: 576/128 = 4.5×

Путь A: Больший kernel (R=56: 256 dims → 8 раундов)
Путь B: Одновременно δe И δa (2 узла × 32 бит = 64/round → 128/2=2 раунда?)
Путь C: Несколько раундов Newton — решаем 4 самых важных
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

def sha_state_at_r(W16, R):
    W=expand_real(W16); s=list(IV)
    for r in range(R):
        a,b,c,d,e,f,g,h=s
        T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add32(Sig0(a),Maj(a,b,c))
        s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
    return s

def sha_compress_full(W16):
    return sha_state_at_r(W16, 64)

def bits_to_words(bits):
    words=[]
    for w in range(16):
        val=0
        for b in range(32):
            if bits[w*32+b]: val|=(1<<b)
        words.append(val)
    return words

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

def gf2_solve(M, target_bits):
    """Solve M·x = target over GF(2). M is list of rows, target is list of bits."""
    n = len(M); m = len(M[0]) if M else 0
    aug = [list(M[i]) + [target_bits[i]] for i in range(n)]
    pivots=[]; ri=0
    for col in range(m):
        pv=-1
        for r in range(ri,n):
            if aug[r][col]: pv=r; break
        if pv==-1: continue
        aug[ri],aug[pv]=aug[pv],aug[ri]
        for r in range(n):
            if r!=ri and aug[r][col]:
                for c in range(m+1): aug[r][c]^=aug[ri][c]
        pivots.append(col); ri+=1
    rank=len(pivots)
    for r in range(rank,n):
        if aug[r][m]: return None, rank
    sol=[0]*m
    for i,pc in enumerate(pivots): sol[pc]=aug[i][m]
    return sol, rank


def compute_jacobian(basis, W_base, target_rounds, target_regs):
    """Compute GF(2) Jacobian for multiple rounds and registers.
    target_regs: list of register indices (0=a, 4=e) per round.
    Returns matrix with 32*len(target_rounds) rows × len(basis) cols."""
    n_basis = len(basis)

    # Reference point: basis[0]
    ref_bits = list(basis[0])
    ref_words = bits_to_words(ref_bits)
    W_ref = [(W_base[j]^ref_words[j]) for j in range(16)]

    # Compute reference δ at all target rounds
    ref_deltas = {}
    for tr, reg in zip(target_rounds, target_regs):
        s_base = sha_state_at_r(W_base, tr+1)
        s_ref = sha_state_at_r(W_ref, tr+1)
        ref_deltas[(tr,reg)] = s_base[reg] ^ s_ref[reg]

    J = []
    target_vec = []

    for tr, reg in zip(target_rounds, target_regs):
        delta_ref = ref_deltas[(tr,reg)]
        for bit in range(32):
            target_vec.append((delta_ref >> bit) & 1)
            row = []
            for i in range(n_basis):
                test_bits = list(ref_bits)
                for j in range(512): test_bits[j] ^= basis[i][j]
                tw = bits_to_words(test_bits)
                Wt = [(W_base[j]^tw[j]) for j in range(16)]
                st = sha_state_at_r(Wt, tr+1)
                dt = st[reg] ^ sha_state_at_r(W_base, tr+1)[reg]
                row.append(((delta_ref>>bit) ^ (dt>>bit)) & 1)
            J.append(row)

    return J, target_vec


# ============================================================
# PATH A: Larger kernel (R=56, dim=256) — more simultaneous rounds
# ============================================================
def path_a(W_base):
    print("="*70)
    print("PATH A: LARGER KERNEL (R=56, dim=256)")
    print("  256 dims / 32 bits = 8 simultaneous rounds possible")
    print("="*70)

    basis = extract_kernel_basis(56)
    n = len(basis)
    print(f"  Kernel dim: {n}")

    for n_rounds in [1, 2, 4, 6, 8]:
        # Pick evenly spaced rounds in gap
        step = max(1, 34 // n_rounds)
        rounds = [17 + i*step for i in range(n_rounds) if 17+i*step < 52]
        rounds = rounds[:n_rounds]
        regs = [4]*len(rounds)  # δe for all

        total_constraints = len(rounds) * 32

        J, target = compute_jacobian(basis, W_base, rounds, regs)
        sol, rank = gf2_solve(J, target)

        solvable = sol is not None
        print(f"\n  {len(rounds)} rounds {rounds}: {total_constraints} constraints, rank={rank}")
        print(f"    Solvable: {'YES ★' if solvable else 'no'}")

        if solvable:
            # Verify
            combined = list(basis[0])
            for i in range(n):
                if sol[i]:
                    for j in range(512): combined[j] ^= basis[i][j]
            dW = bits_to_words(combined)
            W_mod = [(W_base[j]^dW[j]) for j in range(16)]

            for tr in rounds:
                s1 = sha_state_at_r(W_base, tr+1)
                s2 = sha_state_at_r(W_mod, tr+1)
                de = hw(s1[4]^s2[4])
                da = hw(s1[0]^s2[0])
                print(f"    r={tr}: HW(δe)={de}, HW(δa)={da}")


# ============================================================
# PATH B: Both δe AND δa — use all 64 bits per round
# ============================================================
def path_b(W_base):
    print("\n" + "="*70)
    print("PATH B: BOTH δe AND δa — 64 bits/round")
    print("  128 dims / 64 bits = 2 rounds of FULL state control")
    print("="*70)

    basis = extract_kernel_basis(52)
    n = len(basis)

    for n_rounds in [1, 2, 3, 4]:
        step = max(1, 34 // n_rounds)
        rounds = [17 + i*step for i in range(n_rounds) if 17+i*step < 52][:n_rounds]

        # Target both a (reg 0) and e (reg 4)
        all_rounds = []
        all_regs = []
        for r in rounds:
            all_rounds.extend([r, r])
            all_regs.extend([0, 4])  # a and e

        total_constraints = len(rounds) * 64
        J, target = compute_jacobian(basis, W_base, all_rounds, all_regs)
        sol, rank = gf2_solve(J, target)

        print(f"\n  {len(rounds)} rounds, both a+e: {total_constraints} constraints, rank={rank}")
        print(f"    Solvable: {'YES ★' if sol is not None else 'no'}")

        if sol is not None:
            combined = list(basis[0])
            for i in range(n):
                if sol[i]:
                    for j in range(512): combined[j] ^= basis[i][j]
            dW = bits_to_words(combined)
            W_mod = [(W_base[j]^dW[j]) for j in range(16)]
            for tr in rounds:
                s1 = sha_state_at_r(W_base, tr+1)
                s2 = sha_state_at_r(W_mod, tr+1)
                de = hw(s1[4]^s2[4])
                da = hw(s1[0]^s2[0])
                print(f"    r={tr}: HW(δe)={de}, HW(δa)={da}")


# ============================================================
# PATH C: Strategic 4 rounds — maximize coverage with shift reg
# ============================================================
def path_c(W_base):
    print("\n" + "="*70)
    print("PATH C: STRATEGIC 4 ROUNDS — maximize shift register coverage")
    print("  4 rounds × 32 bits = 128 constraints = exactly kernel dim")
    print("  δe[r]=0 → δf[r+1]=0 → δg[r+2]=0 → δh[r+3]=0")
    print("  So 4 rounds with spacing 4 gives 4×4 = 16 covered rounds!")
    print("="*70)

    basis = extract_kernel_basis(52)
    n = len(basis)

    # Best spacing: every 4 rounds (to use shift register propagation)
    configs = [
        [20, 24, 28, 32],       # Early gap
        [24, 28, 32, 36],       # Mid gap
        [28, 32, 36, 40],       # Mid-late
        [32, 36, 40, 44],       # Late gap
        [20, 28, 36, 44],       # Spread out
        [17, 25, 33, 41],       # Max spread
    ]

    for rounds in configs:
        regs = [4]*4  # δe only
        J, target = compute_jacobian(basis, W_base, rounds, regs)
        sol, rank = gf2_solve(J, target)

        print(f"\n  Rounds {rounds}: rank={rank}, solvable={'YES' if sol is not None else 'no'}")

        if sol is not None:
            combined = list(basis[0])
            for i in range(n):
                if sol[i]:
                    for j in range(512): combined[j] ^= basis[i][j]
            dW = bits_to_words(combined)
            W_mod = [(W_base[j]^dW[j]) for j in range(16)]

            # Show FULL gap profile
            zeros_e = 0; zeros_any = 0
            for r in range(17, 52):
                s1 = sha_state_at_r(W_base, r+1)
                s2 = sha_state_at_r(W_mod, r+1)
                de = hw(s1[4]^s2[4])
                da = hw(s1[0]^s2[0])
                if de == 0: zeros_e += 1
                if de == 0 or da == 0: zeros_any += 1
                if de == 0 or da == 0:
                    regs_zero = []
                    for i,name in enumerate("abcdefgh"):
                        if s1[i]==s2[i]: regs_zero.append(name)
                    print(f"    r={r}: δe={de:2d} δa={da:2d} zeros={','.join(regs_zero)}")

            print(f"    Total δe=0 rounds: {zeros_e}/35")
            print(f"    Total any-zero rounds: {zeros_any}/35")

            # Hash diff
            H1 = [add32(IV[i], sha_state_at_r(W_base,64)[i]) for i in range(8)]
            H2 = [add32(IV[i], sha_state_at_r(W_mod,64)[i]) for i in range(8)]
            hash_hw = sum(hw(H1[i]^H2[i]) for i in range(8))
            print(f"    Hash diff: HW = {hash_hw}")


if __name__ == '__main__':
    print("="*70)
    print("SCF: NEWTON SCALING — 3 PATHS")
    print("="*70)

    W_base = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

    path_a(W_base)
    path_b(W_base)
    path_c(W_base)
