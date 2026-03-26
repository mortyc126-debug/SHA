#!/usr/bin/env python3
"""
SCF: NEWTON НАД Z — с carry, не GF(2).

Проблема: GF(2) Newton находит ядро, но все решения тривиальны (δW=0).
Причина: GF(2) не видит carry. Carry — 48.9% дифференциала.

Решение: Newton над Z_{2^32}. Якобиан = числовая производная.
J[i][j] = (f(x + e_j)[i] - f(x)[i]) mod 2^32 — АРИФМЕТИЧЕСКИЙ.

Для одного раунда: δe = δd + δT1.
δT1 = δh + δΣ1(e) + δCh + δW[r] (с carry!)

Newton над Z: решаем J_Z · Δc ≡ -f(c) (mod 2^32).
Это ЛИНЕЙНАЯ система mod 2^32 — решаема через Hensel/CRT!
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
def sub32(a,b): return (a-b)&MASK32
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


# ============================================================
# Z-Newton: arithmetic Jacobian and solving
# ============================================================

def compute_z_jacobian_single_round(basis, W_base, target_r, target_reg=4):
    """Compute ARITHMETIC Jacobian: ∂(δ_reg[r]) / ∂(basis_coeff_i).

    For each basis vector i: compute δ_reg when using basis[i] as δW.
    J_Z[i] = δ_reg(W_base ⊕ basis_words[i]) — a 32-bit value.

    This is a 1×n_basis matrix over Z_{2^32} (one target word).
    For n target words: n×n_basis matrix.
    """
    n_basis = len(basis)
    basis_words = [bits_to_words(v) for v in basis]

    s_base = sha_state_at_r(W_base, target_r + 1)

    J_col = []  # Each entry = δ_reg as 32-bit integer
    for i in range(n_basis):
        W_mod = [(W_base[j] ^ basis_words[i][j]) for j in range(16)]
        s_mod = sha_state_at_r(W_mod, target_r + 1)
        delta = sub32(s_mod[target_reg], s_base[target_reg])
        J_col.append(delta)

    return J_col, s_base[target_reg]


def newton_z_single(basis, W_base, target_r, max_iter=30):
    """Newton over Z_{2^32} for a single round.

    Goal: find kernel coefficients c such that
    state_base[r][reg] = state_mod[r][reg] (additive collision on one word).

    Method: p-adic Newton (Hensel lifting).
    Solve mod 2, then mod 4, then mod 8, ..., mod 2^32.
    """
    n_basis = len(basis)
    basis_words = [bits_to_words(v) for v in basis]

    print(f"\n  Newton-Z for r={target_r}:")

    # Current coefficient vector (over Z/2^k)
    c = [0] * n_basis

    for k in range(1, 33):
        mod = 1 << k
        mask = mod - 1

        # Compute current δe mod 2^k
        bits = [0]*512
        for i in range(n_basis):
            if c[i]:
                for j in range(512): bits[j] ^= basis[i][j]

        if all(b==0 for b in bits):
            # Need to start with something nontrivial
            c[0] = 1
            for j in range(512): bits[j] ^= basis[0][j]

        dW = bits_to_words(bits)
        W_mod = [(W_base[j] ^ dW[j]) for j in range(16)]
        s_base = sha_state_at_r(W_base, target_r + 1)
        s_mod = sha_state_at_r(W_mod, target_r + 1)

        de = sub32(s_mod[4], s_base[4]) & mask

        if de == 0:
            if k <= 5 or k == 32:
                print(f"    k={k:2d} (mod 2^{k}): δe ≡ 0 ✓")
            continue

        # Need to correct: find Δc such that J·Δc ≡ -de (mod 2)
        # At level k: we only flip the k-th bit of some coefficient
        # Try each basis vector: does flipping it reduce de mod 2^k?

        best_flip = -1
        best_de = de

        for i in range(n_basis):
            test_c = list(c)
            test_c[i] ^= (1 << (k-1))  # Flip bit k-1 of coeff i

            test_bits = [0]*512
            for ii in range(n_basis):
                if test_c[ii] & 1:  # Only GF(2) part matters for XOR
                    for j in range(512): test_bits[j] ^= basis[ii][j]

            # Hmm, c[i] can be multi-bit but basis is GF(2)...
            # Actually: kernel coefficients are GF(2) (0 or 1).
            # We can't do p-adic Newton on GF(2) coefficients.
            break

        if k <= 5 or k == 32:
            print(f"    k={k:2d} (mod 2^{k}): δe ≡ {de} (mod {mod}), HW={hw(de)}")

        # GF(2) coefficients mean we can't do Hensel lifting on c.
        # Instead: try all 2^min(n_basis, 20) combinations for small kernel.
        break

    return c


# ============================================================
# Better approach: ADDITIVE Newton on the MESSAGE directly
# ============================================================
def newton_additive(W_base, target_r, max_iter=50):
    """Newton over Z: perturb W[0..15] ADDITIVELY (not XOR).

    δW_add = W_mod - W_base (mod 2^32) — additive perturbation.
    Jacobian: J[reg][word] = ∂(state_reg) / ∂(W_word) mod 2^32.

    This is REAL Newton, not GF(2)!
    """
    print(f"\n  Additive Newton for r={target_r}:")

    W_current = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    s_target = sha_state_at_r(W_base, target_r + 1)

    for iteration in range(max_iter):
        s_current = sha_state_at_r(W_current, target_r + 1)

        # Residual: what we need to zero
        de = sub32(s_current[4], s_target[4])
        da = sub32(s_current[0], s_target[0])

        hw_de = hw(de)
        hw_da = hw(da)

        if de == 0 and da == 0:
            print(f"    Iter {iteration}: ★★★ CONVERGED! δe=0 AND δa=0")
            # Check if W_current ≠ W_base
            if W_current != W_base:
                print(f"    W_current ≠ W_base: NONTRIVIAL!")
                H1 = sha_compress_full(W_base)
                H2 = sha_compress_full(W_current)
                hash_hw = sum(hw(H1[i]^H2[i]) for i in range(8))
                print(f"    Hash diff: HW = {hash_hw}")
            return W_current, True

        if de == 0:
            if iteration < 5 or iteration % 10 == 0:
                print(f"    Iter {iteration}: δe=0 ✓, δa HW={hw_da}")
        elif iteration < 5 or iteration % 10 == 0:
            print(f"    Iter {iteration}: δe HW={hw_de}, δa HW={hw_da}")

        # Compute Jacobian: ∂(e[r+1]) / ∂(W[word])
        # Approximate: flip each W word by +1, measure Δe
        best_word = -1
        best_residual_hw = 32

        for word in range(16):
            W_test = list(W_current)
            # Try additive perturbation that cancels de
            # If J[word] ≈ sensitivity, then ΔW = -de / J[word]

            # Measure sensitivity: (e(W+e_word) - e(W)) mod 2^32
            W_test[word] = add32(W_current[word], 1)
            s_test = sha_state_at_r(W_test, target_r + 1)
            sensitivity = sub32(s_test[4], s_current[4])

            if sensitivity == 0:
                continue

            # Newton correction: ΔW = -de * sensitivity^{-1} mod 2^32
            # Finding modular inverse
            try:
                # sensitivity must be odd for inverse to exist
                if sensitivity & 1:
                    inv = pow(sensitivity, -1, 1 << 32)
                    correction = (-(de) * inv) & MASK32

                    W_trial = list(W_current)
                    W_trial[word] = add32(W_current[word], correction)
                    s_trial = sha_state_at_r(W_trial, target_r + 1)
                    de_trial = sub32(s_trial[4], s_target[4])

                    if hw(de_trial) < best_residual_hw:
                        best_residual_hw = hw(de_trial)
                        best_word = word
                        best_correction = correction
                        best_W = list(W_trial)
            except (ValueError, ZeroDivisionError):
                pass

        if best_word >= 0:
            W_current = best_W
        else:
            # Random perturbation
            word = int.from_bytes(os.urandom(1),'big') % 16
            W_current[word] ^= int.from_bytes(os.urandom(4),'big')

    s_final = sha_state_at_r(W_current, target_r + 1)
    de_final = sub32(s_final[4], s_target[4])
    print(f"    Final: δe HW={hw(de_final)}")
    return W_current, False


# ============================================================
# EXP: Newton-Z for single round, then chain
# ============================================================
if __name__ == '__main__':
    print("="*70)
    print("SCF: NEWTON НАД Z_{2^32}")
    print("="*70)

    W_base = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

    # Single round Newton-Z
    for tr in [20, 30, 40]:
        W_result, success = newton_additive(W_base, tr, max_iter=40)
        if success:
            n_diff = sum(1 for i in range(16) if W_base[i] != W_result[i])
            print(f"  r={tr}: converged! {n_diff} words differ")

    # Chain: solve r=20, then use result to solve r=21, etc.
    print("\n" + "="*70)
    print("CHAIN: sequential Newton-Z through gap")
    print("="*70)

    W_current = list(W_base)
    W_target = list(W_base)  # We want state_current = state_target at each round

    # Actually: we want to find W2 such that SHA(W2) = SHA(W_base)
    # Newton-Z on the HASH directly
    print("\n  Newton-Z on FULL HASH:")
    W2 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    H_target = sha_compress_full(W_base)

    for iteration in range(100):
        H_current = sha_compress_full(W2)
        dH = [sub32(H_current[i], H_target[i]) for i in range(8)]
        total_hw = sum(hw(d) for d in dH)

        if total_hw == 0:
            print(f"  Iter {iteration}: ★★★ HASH COLLISION!")
            print(f"  W_base: {[hex(w) for w in W_base[:4]]}...")
            print(f"  W2:     {[hex(w) for w in W2[:4]]}...")
            break

        if iteration < 10 or iteration % 20 == 0:
            print(f"  Iter {iteration}: HW(δH) = {total_hw}")

        # Pick register with largest diff, Newton-correct one word
        worst = max(range(8), key=lambda i: hw(dH[i]))
        target_val = dH[worst]

        best_word = -1
        best_hw = 32

        for word in range(16):
            W_test = list(W2)
            W_test[word] = add32(W2[word], 1)
            H_test = sha_compress_full(W_test)
            sens = sub32(H_test[worst], H_current[worst])

            if sens & 1:  # Invertible
                try:
                    inv = pow(sens, -1, 1 << 32)
                    corr = (-(target_val) * inv) & MASK32
                    W_trial = list(W2)
                    W_trial[word] = add32(W2[word], corr)
                    H_trial = sha_compress_full(W_trial)
                    new_hw = hw(sub32(H_trial[worst], H_target[worst]))
                    if new_hw < best_hw:
                        best_hw = new_hw
                        best_word = word
                        best_W = list(W_trial)
                except:
                    pass

        if best_word >= 0:
            W2 = best_W
        else:
            W2[iteration % 16] ^= int.from_bytes(os.urandom(4),'big')
