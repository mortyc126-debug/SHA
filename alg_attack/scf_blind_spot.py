#!/usr/bin/env python3
"""
SCF: BLIND SPOT EXPLOIT — O(2^16) collision через kernel(∂H/∂IV).

ОТКРЫТИЕ: rank(∂H/∂IV) = 254-255 для конкретного M2.
Kernel = 1-2 направления δIV где hash НЕ МЕНЯЕТСЯ.

СТРАТЕГИЯ:
1. Выбрать M2 (фиксировать)
2. Найти kernel(∂H/∂IV) — 1-2 вектора δIV* (32-64 бит target)
3. Block 1: найти M1, M1' такие что δH1 = δIV* (birthday O(2^16))
4. Тогда: hash(M1||M2) = hash(M1'||M2) → COLLISION!

Стоимость: O(2^16) для 1-dim kernel, O(2^32) для 2-dim.
vs birthday O(2^128). Speedup: 2^{96}-2^{112}!

НО: kernel GF(2) ≠ kernel Z. Нужно проверить:
реально ли δIV из GF(2) kernel даёт δH=0 в РЕАЛЬНОЙ SHA?
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

def sha_hash(W16, iv):
    W=expand_real(W16); s=list(iv)
    for r in range(64):
        a,b,c,d,e,f,g,h=s
        T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add32(Sig0(a),Maj(a,b,c))
        s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
    return [add32(iv[i],s[i]) for i in range(8)]


def find_iv_kernel(M2, iv_base):
    """Find GF(2) kernel of ∂H/∂IV for fixed M2 at iv_base."""
    H_base = sha_hash(M2, iv_base)

    # Build 256×256 GF(2) Jacobian
    rows = []  # Each row = 256-bit integer
    for out_reg in range(8):
        for out_bit in range(32):
            row = 0
            for in_reg in range(8):
                for in_bit in range(32):
                    iv_flip = list(iv_base)
                    iv_flip[in_reg] ^= (1 << in_bit)
                    H_flip = sha_hash(M2, iv_flip)
                    if (H_base[out_reg] ^ H_flip[out_reg]) >> out_bit & 1:
                        row |= 1 << (in_reg*32 + in_bit)
            rows.append(row)

    # RREF to find kernel
    m = list(rows)
    pivots = []; ri = 0
    for col in range(256):
        mask = 1 << col
        pv = -1
        for r in range(ri, 256):
            if m[r] & mask: pv=r; break
        if pv == -1: continue
        m[ri],m[pv] = m[pv],m[ri]
        for r in range(256):
            if r != ri and m[r] & mask:
                m[r] ^= m[ri]
        pivots.append(col); ri += 1

    rank = len(pivots)
    free_cols = [c for c in range(256) if c not in pivots]

    kernel_vecs = []
    for fc in free_cols:
        vec = 0
        vec |= (1 << fc)
        for i, pc in enumerate(pivots):
            if m[i] & (1 << fc):
                vec |= (1 << pc)
        kernel_vecs.append(vec)

    return kernel_vecs, rank


def int_to_iv(val):
    """256-bit int → 8×32-bit IV."""
    return [(val >> (i*32)) & MASK32 for i in range(8)]

def iv_to_int(iv):
    val = 0
    for i in range(8):
        val |= iv[i] << (i*32)
    return val


# ============================================================
# EXP 1: Find kernel, verify it ACTUALLY gives δH≈0
# ============================================================
def exp1_kernel_verify(N):
    print("="*70)
    print("EXP 1: FIND IV-KERNEL AND VERIFY")
    print("="*70)

    for trial in range(N):
        M2 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        iv_base = [int.from_bytes(os.urandom(4),'big') for _ in range(8)]

        kernel_vecs, rank = find_iv_kernel(M2, iv_base)
        print(f"\n  Trial {trial}: rank={rank}, kernel_dim={len(kernel_vecs)}")

        for ki, kv in enumerate(kernel_vecs[:3]):
            # Verify: does δIV = kernel_vec give δH=0?
            δiv = int_to_iv(kv)
            iv_modified = [iv_base[i] ^ δiv[i] for i in range(8)]

            H1 = sha_hash(M2, iv_base)
            H2 = sha_hash(M2, iv_modified)

            δH_hw = sum(hw(H1[i]^H2[i]) for i in range(8))
            δIV_hw = sum(hw(δiv[i]) for i in range(8))

            marker = " ★★★ COLLISION!" if δH_hw == 0 else (" ★ near!" if δH_hw < 10 else "")
            print(f"    kernel[{ki}]: HW(δIV)={δIV_hw}, HW(δH)={δH_hw}{marker}")

            if δH_hw == 0 and δIV_hw > 0:
                print(f"    ★★★ REAL IV-COLLISION FOUND!")
                print(f"    δIV = {[hex(d) for d in δiv]}")
                return kv, M2, iv_base

    return None, None, None


# ============================================================
# EXP 2: If kernel gives δH≈small, can Block 1 hit the target?
# ============================================================
def exp2_block1_targeting(kernel_vec, M2, iv_target_base, N_search):
    """Block 1 must produce δH1 = δIV (from kernel).
    Search: find M1, M1' with sha(IV, M1) ⊕ sha(IV, M1') = δIV."""

    if kernel_vec is None:
        print("\n  No kernel vector to target. Using random search.")
        # Use a random near-zero kernel-ish direction
        M2 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        iv_base = [int.from_bytes(os.urandom(4),'big') for _ in range(8)]
        kernel_vecs, rank = find_iv_kernel(M2, iv_base)
        if not kernel_vecs:
            print("  No kernel found. Abort.")
            return
        kernel_vec = kernel_vecs[0]
        iv_target_base = iv_base

    δIV_target = int_to_iv(kernel_vec)
    target_hw = sum(hw(d) for d in δIV_target)

    print(f"\n" + "="*70)
    print(f"EXP 2: BLOCK 1 TARGETING")
    print(f"  Target: δH1 = δIV* (kernel vector, HW={target_hw})")
    print(f"  Method: birthday — find M1,M1' with sha(IV,M1)⊕sha(IV,M1')=δIV*")
    print(f"="*70)

    # Birthday: store hash(IV, M1), look for hash(IV, M1') = hash + δIV*
    # Truncated birthday: match on first register

    H_target_xor = δIV_target  # We want H(M1) ⊕ H(M1') = this

    seen = {}  # H[0] → (M1, full_hash)

    for trial in range(N_search):
        M1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        H1 = sha_hash(M1, IV)

        # We want another M1' with H1' = H1 ⊕ δIV_target
        # Equivalent: H1'[0] = H1[0] ⊕ δIV_target[0]
        H_target_0 = H1[0] ^ δIV_target[0]

        if H_target_0 in seen:
            M1_other, H_other = seen[H_target_0]
            # Check full match
            full_match = all(H1[i] ^ H_other[i] == δIV_target[i] for i in range(8))
            if full_match:
                print(f"  ★★★ BLOCK 1 MATCH at trial {trial}!")
                print(f"    M1  = {[hex(m) for m in M1[:4]]}...")
                print(f"    M1' = {[hex(m) for m in M1_other[:4]]}...")

                # Verify 2-block collision
                H1_a = sha_hash(M1, IV)
                H1_b = sha_hash(M1_other, IV)
                H2_a = sha_hash(M2, H1_a)
                H2_b = sha_hash(M2, H1_b)
                δH_final = sum(hw(H2_a[i]^H2_b[i]) for i in range(8))

                print(f"    δH_block1 = {sum(hw(H1_a[i]^H1_b[i]) for i in range(8))} (should = {target_hw})")
                print(f"    δH_final = {δH_final}")

                if δH_final == 0:
                    print(f"    ★★★★★ TWO-BLOCK COLLISION FOUND! ★★★★★")
                else:
                    print(f"    GF(2) kernel ≠ Z kernel. δH_final = {δH_final}")
                return

        # Also store this hash
        seen[H1[0]] = (M1, H1)

        if trial % 50000 == 0 and trial > 0:
            print(f"    {trial} trials, {len(seen)} stored...")

    # Birthday probability estimate
    print(f"\n  No match in {N_search} trials ({len(seen)} stored)")
    print(f"  Birthday expects match at ~2^{16} = {1<<16} for 32-bit target")
    if N_search < (1<<16):
        print(f"  Need more trials! (tried {N_search}, need ~{1<<16})")


# ============================================================
# EXP 3: Mass kernel search — find M2 with largest kernel
# ============================================================
def exp3_mass_kernel(N):
    print(f"\n" + "="*70)
    print(f"EXP 3: MASS KERNEL SEARCH — find M2 with biggest kernel")
    print(f"="*70)

    from collections import Counter
    kernel_dist = Counter()
    best_kernel = 0
    best_M2 = None

    for trial in range(N):
        M2 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        iv = [int.from_bytes(os.urandom(4),'big') for _ in range(8)]

        _, rank = find_iv_kernel(M2, iv)
        kdim = 256 - rank
        kernel_dist[kdim] += 1

        if kdim > best_kernel:
            best_kernel = kdim
            best_M2 = (list(M2), list(iv))
            print(f"  [{trial}] NEW BEST: kernel_dim = {kdim}")

    print(f"\n  Distribution: {dict(sorted(kernel_dist.items()))}")
    print(f"  Best kernel: {best_kernel}")

    return best_M2


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    kv, M2, iv = exp1_kernel_verify(N)
    exp3_mass_kernel(min(N*2, 10))
    exp2_block1_targeting(kv, M2, iv, N_search=min(N*20000, 100000))
