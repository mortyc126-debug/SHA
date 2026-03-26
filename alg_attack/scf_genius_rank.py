#!/usr/bin/env python3
"""
SCF: GENIUS MOVE — Output rank of kernel image.

The kernel of Schedule_GF2[52..63] has dim=128.
Each kernel element δW produces a hash difference δH = SHA(W⊕δW) ⊕ SHA(W).
Question: what is rank of the map kernel → δH?

If rank = 128 → output is 128-dim → birthday on 2^128 elements costs 2^64
If rank < 128 → output is compressed → birthday even cheaper
If rank = 256 → impossible (128 inputs can't span 256), so max is 128

KEY INSIGHT: 2^128 kernel elements live in a 128-dim input subspace.
Their OUTPUTS live in at most 128-dim subspace of the 256-dim hash space.
This means 128 bits of hash are DETERMINED by the other 128.
Birthday in 128-dim output costs O(2^64), not O(2^128)!

BUT: this only helps if we ACTUALLY work within the kernel coset.
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

def sha_compress(W16, iv=None):
    if iv is None: iv = IV
    W=expand_real(W16); s=list(iv)
    for r in range(64):
        a,b,c,d,e,f,g,h=s
        T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add32(Sig0(a),Maj(a,b,c))
        s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
    return [add32(iv[i],s[i]) for i in range(8)]

def bits_to_words(bits):
    words = []
    for w in range(16):
        val = 0
        for b in range(32):
            if bits[w*32+b]: val |= (1<<b)
        words.append(val)
    return words

def state_to_bits(state):
    bits = []
    for w in state:
        for b in range(32):
            bits.append((w>>b)&1)
    return bits

def extract_kernel_basis(R_start=52):
    n_in_bits = 512
    M = []
    for out_word in range(R_start, 64):
        for out_bit in range(32):
            row = [0]*n_in_bits
            for in_word in range(16):
                for in_bit in range(32):
                    W_test = [0]*16
                    W_test[in_word] = 1<<in_bit
                    Wexp = expand_xor(W_test)
                    if (Wexp[out_word]>>out_bit)&1:
                        row[in_word*32+in_bit] = 1
            M.append(row)
    m = [list(row) for row in M]
    n_rows, n_cols = len(m), n_in_bits
    pivot_cols = []
    row_idx = 0
    for col in range(n_cols):
        pivot = -1
        for r in range(row_idx, n_rows):
            if m[r][col]==1: pivot=r; break
        if pivot==-1: continue
        m[row_idx],m[pivot] = m[pivot],m[row_idx]
        for r in range(n_rows):
            if r!=row_idx and m[r][col]==1:
                for c in range(n_cols): m[r][c]^=m[row_idx][c]
        pivot_cols.append(col)
        row_idx += 1
    free_cols = [c for c in range(n_cols) if c not in pivot_cols]
    kernel_basis = []
    for fc in free_cols:
        vec = [0]*n_cols; vec[fc]=1
        for i,pc in enumerate(pivot_cols):
            if m[i][fc]==1: vec[pc]=1
        kernel_basis.append(vec)
    return kernel_basis

def gf2_rank(matrix):
    if not matrix: return 0
    m = [list(row) for row in matrix]
    nrows,ncols = len(m),len(m[0])
    rank = 0
    for col in range(ncols):
        pivot = -1
        for r in range(rank, nrows):
            if m[r][col]==1: pivot=r; break
        if pivot==-1: continue
        m[rank],m[pivot] = m[pivot],m[rank]
        for r in range(nrows):
            if r!=rank and m[r][col]==1:
                for c in range(ncols): m[r][c]^=m[rank][c]
        rank += 1
    return rank


# ============================================================
# THE KEY EXPERIMENT: Rank of kernel → δhash map
# ============================================================
def measure_output_rank(basis, N_bases=5):
    print("="*70)
    print("THE KEY EXPERIMENT: RANK OF KERNEL → δHASH")
    print("  kernel has dim=128. Map: δW → δH = SHA(W⊕δW) ⊕ SHA(W)")
    print("  If rank(image) < 128 → output is COMPRESSED → cheaper birthday")
    print("="*70)

    for base_idx in range(N_bases):
        W_base = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        H_base = sha_compress(W_base)

        # Compute δH for each basis vector
        output_rows = []
        for i, vec in enumerate(basis):
            dW = bits_to_words(vec)
            W_mod = [(W_base[j]^dW[j]) for j in range(16)]
            H_mod = sha_compress(W_mod)
            dH = [H_base[j]^H_mod[j] for j in range(8)]
            output_rows.append(state_to_bits(dH))

        # Rank of output matrix (128 rows × 256 cols)
        r = gf2_rank(output_rows)

        print(f"\n  Base {base_idx}: rank(kernel → δH) = {r} / {len(basis)}")

        if r < len(basis):
            defect = len(basis) - r
            print(f"  ★★★ RANK DEFECT = {defect}!")
            print(f"  Output image has dimension {r}, not {len(basis)}")
            print(f"  {defect} bits of hash are DEPENDENT on the rest!")
            print(f"  Birthday cost: O(2^{r//2}) instead of O(2^{128//2})=O(2^64)")
        else:
            print(f"  Full rank: image = {r}-dimensional (no compression)")

    # Also check with random (non-kernel) inputs for comparison
    print(f"\n  Comparison: random 128 δW vectors (not from kernel)")
    for base_idx in range(2):
        W_base = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        H_base = sha_compress(W_base)

        output_rows = []
        for i in range(128):
            dW = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
            W_mod = [(W_base[j]^dW[j]) for j in range(16)]
            H_mod = sha_compress(W_mod)
            dH = [H_base[j]^H_mod[j] for j in range(8)]
            output_rows.append(state_to_bits(dH))

        r = gf2_rank(output_rows)
        print(f"  Random base {base_idx}: rank = {r} / 128")


# ============================================================
# Check: different kernel sizes
# ============================================================
def measure_multiple_kernels():
    print("\n" + "="*70)
    print("RANK VS KERNEL SIZE")
    print("="*70)

    for R_start in [60, 56, 52]:
        basis = extract_kernel_basis(R_start)
        kdim = len(basis)

        W_base = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        H_base = sha_compress(W_base)

        output_rows = []
        for vec in basis:
            dW = bits_to_words(vec)
            W_mod = [(W_base[j]^dW[j]) for j in range(16)]
            H_mod = sha_compress(W_mod)
            dH = [H_base[j]^H_mod[j] for j in range(8)]
            output_rows.append(state_to_bits(dH))

        r = gf2_rank(output_rows)
        defect = kdim - r

        print(f"  R={R_start}: kernel_dim={kdim}, output_rank={r}, defect={defect}")

        if defect > 0:
            # This would be huge!
            birthday_cost = r // 2
            print(f"  ★★★ Birthday in kernel coset: O(2^{birthday_cost})")
        else:
            print(f"  Full rank → birthday = O(2^{kdim//2})")


# ============================================================
# Birthday within kernel coset
# ============================================================
def birthday_in_kernel(basis, N_trials):
    print("\n" + "="*70)
    print("BIRTHDAY WITHIN KERNEL COSET")
    print(f"  Coset size = 2^{len(basis)}. Expected collisions: ~0.5")
    print(f"  (We test with {N_trials} samples — partial birthday)")
    print("="*70)

    W_base = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

    # Generate random kernel-coset elements and hash them
    hashes = {}
    best_collision_hw = 256
    n_basis = len(basis)

    for trial in range(N_trials):
        # Random kernel element = XOR of random subset of basis
        bits = [0]*512
        mask = int.from_bytes(os.urandom(16),'big')  # 128-bit random
        for i in range(n_basis):
            if (mask >> i) & 1:
                for j in range(512):
                    bits[j] ^= basis[i][j]

        if all(b==0 for b in bits):
            continue

        dW = bits_to_words(bits)
        W_mod = [(W_base[j]^dW[j]) for j in range(16)]
        H = sha_compress(W_mod)
        H_tuple = tuple(H)

        # Check for collision with previous
        # Truncated hash for partial collision detection
        for trunc_bits in [32, 64, 96, 128]:
            mask_val = (1 << trunc_bits) - 1
            H_trunc = tuple(h & mask_val for h in H)

            # Store and check
            key = (trunc_bits, H_trunc)
            if key not in hashes:
                hashes[key] = (trial, dW)
            elif hashes[key][1] != dW:
                # Partial collision found!
                old_trial, old_dW = hashes[key]
                # Measure full collision distance
                W_old = [(W_base[j]^old_dW[j]) for j in range(16)]
                H_old = sha_compress(W_old)
                full_hw = sum(hw(H[i]^H_old[i]) for i in range(8))

                if full_hw < best_collision_hw:
                    best_collision_hw = full_hw
                    if trunc_bits <= 64:
                        print(f"  Partial collision ({trunc_bits}-bit): "
                              f"trial {old_trial} vs {trial}, full HW={full_hw}")

    print(f"\n  Best near-collision in kernel coset ({N_trials} trials): HW={best_collision_hw}")
    print(f"  Expected for {N_trials} random: HW≈{256 - int(2*math.log2(N_trials+1))}")

    # Compare with non-kernel random pairs
    best_random = 256
    for trial in range(N_trials):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        H1 = sha_compress(W1); H2 = sha_compress(W2)
        d = sum(hw(H1[i]^H2[i]) for i in range(8))
        if d < best_random: best_random = d

    print(f"  Best random pair: HW={best_random}")
    print(f"  Kernel advantage: {best_random - best_collision_hw:+d} bits")


import math

if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 200

    print("="*70)
    print("SCF: GENIUS MOVE — OUTPUT RANK OF KERNEL IMAGE")
    print("="*70)

    print("\nExtracting kernel (R=52)...")
    basis = extract_kernel_basis(52)
    print(f"  dim = {len(basis)}")

    measure_output_rank(basis, N_bases=min(N//20+1, 5))
    measure_multiple_kernels()
    birthday_in_kernel(basis, N_trials=min(N*50, 20000))

    print("\n" + "="*70)
    print("ФИНАЛЬНЫЙ АНАЛИЗ")
    print("="*70)
    print(f"""
  Schedule kernel: dim={len(basis)} для δW_xor[52..63]=0.
  Map kernel → δhash: rank измерен выше.

  Если rank = 128 (полный):
    Birthday в coset: O(2^64) time, O(2^64) memory
    НО: это СТАНДАРТНЫЙ birthday на 2^128 элементах!
    Не даёт улучшения над общим birthday O(2^128).

  Если rank < 128 (сжатие):
    δhash лежит в подпространстве dim < 128
    → больше коллизий в coset → дешевле birthday!
    Это было бы НАСТОЯЩИМ прорывом.

  Значение для SHA-256:
    Если rank < 128 → schedule kernel создаёт СТРУКТУРНУЮ слабость
    Если rank = 128 → schedule kernel — лишь красивая математика,
                       не дающая криптографического преимущества.
    """)
