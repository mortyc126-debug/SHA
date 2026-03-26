#!/usr/bin/env python3
"""
SCF: ПРОБИВАЕМ ABSORBING STATE — Wang-in-the-Gap

Shadow Algebra: absorbing state выходит ТОЛЬКО при δW=cancel.
Schedule kernel: 128 dims свободы в δW[0..15].
Каждый dim → разный δW_real на каждом раунде gap.

Стратегия: найти kernel element где δW_real[r] = -δT1_rest[r]
→ δe[r+1]=0 → shift register даёт ещё 3 нуля → ПРОБОЙ.

128 dims / 32 bits = 4 одновременных условия → 4 пробоя в gap.
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

def sha_states(W16):
    W=expand_real(W16); s=list(IV); states=[list(s)]
    for r in range(64):
        a,b,c,d,e,f,g,h=s
        T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add32(Sig0(a),Maj(a,b,c))
        s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
        states.append(list(s))
    return states, W

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


# ============================================================
# EXP 1: For each gap round, find kernel element that zeros δe
# ============================================================
def exp1_wang_in_gap(N_search):
    print("="*70)
    print("EXP 1: WANG-IN-THE-GAP — zero δe at specific gap rounds")
    print("="*70)

    basis = extract_kernel_basis(52)
    n_basis = len(basis)
    print(f"  Kernel dim: {n_basis}")

    # Precompute basis words
    basis_words = [bits_to_words(v) for v in basis]

    W_base = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    states_base, Wexp_base = sha_states(W_base)

    print(f"\n  For each gap round r: find kernel δW that minimizes HW(δe[r+1])")
    print(f"  {'r':>3} | {'best HW(δe)':>12} | {'tries':>6} | status")
    print("  " + "-"*50)

    for target_r in range(20, 52, 4):
        best_hw = 32
        best_dW = None

        for trial in range(N_search):
            # Random kernel element
            bits = [0]*512
            mask = int.from_bytes(os.urandom(16),'big')
            for i in range(n_basis):
                if (mask>>i)&1:
                    for j in range(512):
                        bits[j] ^= basis[i][j]
            if all(b==0 for b in bits): continue

            dW = bits_to_words(bits)
            W_mod = [(W_base[j]^dW[j]) for j in range(16)]
            states_mod, _ = sha_states(W_mod)

            # Check δe at target_r+1
            de = states_base[target_r+1][4] ^ states_mod[target_r+1][4]
            h = hw(de)
            if h < best_hw:
                best_hw = h
                best_dW = dW

        marker = " ★★★ ZERO!" if best_hw==0 else (" ★★" if best_hw<=3 else (" ★" if best_hw<=8 else ""))
        print(f"  {target_r:3d} | {best_hw:12d} | {N_search:6d} |{marker}")

    return basis


# ============================================================
# EXP 2: Simultaneous zeros — multiple gap rounds
# ============================================================
def exp2_simultaneous(basis, N_search):
    print("\n" + "="*70)
    print("EXP 2: SIMULTANEOUS ZEROS — multiple gap rounds at once")
    print("="*70)

    basis_words = [bits_to_words(v) for v in basis]
    n_basis = len(basis)

    W_base = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    states_base, _ = sha_states(W_base)

    # Score: total HW of δe at target rounds
    target_rounds = [24, 32, 40, 48]  # 4 targets, spaced 8 apart

    def score(dW):
        W_mod = [(W_base[j]^dW[j]) for j in range(16)]
        states_mod, _ = sha_states(W_mod)
        total = 0
        for tr in target_rounds:
            total += hw(states_base[tr][4] ^ states_mod[tr][4])
        return total

    # SA search over kernel
    import math
    cur_bits = list(basis[0])
    cur_words = bits_to_words(cur_bits)
    cur_score = score(cur_words)
    best_score = cur_score
    best_words = list(cur_words)

    for it in range(N_search):
        trial_bits = list(cur_bits)
        idx = int.from_bytes(os.urandom(2),'big') % n_basis
        for j in range(512): trial_bits[j] ^= basis[idx][j]
        if all(b==0 for b in trial_bits): continue

        tw = bits_to_words(trial_bits)
        sc = score(tw)

        T = max(0.01, 1.0-it/N_search)
        if sc < cur_score or \
           math.exp(-(sc-cur_score)/(T*3)) > int.from_bytes(os.urandom(4),'big')/(1<<32):
            cur_bits = trial_bits
            cur_words = tw
            cur_score = sc
        if cur_score < best_score:
            best_score = cur_score
            best_words = list(cur_words)

    # Analyze best
    W_mod = [(W_base[j]^best_words[j]) for j in range(16)]
    states_mod, _ = sha_states(W_mod)

    print(f"\n  Target rounds: {target_rounds}")
    print(f"  Best total HW(δe) = {best_score} (random ≈ {len(target_rounds)*16})")
    for tr in target_rounds:
        de = hw(states_base[tr][4] ^ states_mod[tr][4])
        da = hw(states_base[tr][0] ^ states_mod[tr][0])
        print(f"    r={tr}: HW(δe)={de:2d}, HW(δa)={da:2d}")

    # Full profile — count ALL zero-register instances
    print(f"\n  Full gap profile (zeros in any register):")
    n_zeros_total = 0
    for r in range(17, 53):
        zeros = []
        regs = "abcdefgh"
        for i in range(8):
            if states_base[r][i] == states_mod[r][i]:
                zeros.append(regs[i])
                n_zeros_total += 1
        if zeros:
            print(f"    r={r}: δ=0 in registers: {','.join(zeros)}")

    print(f"\n  Total zero-register instances in gap: {n_zeros_total}")
    print(f"  (Expected random: ~{35*8/2**32:.6f} per register ≈ 0)")


# ============================================================
# EXP 3: Best possible — minimize TOTAL differential in gap
# ============================================================
def exp3_minimize_gap(basis, N_search):
    print("\n" + "="*70)
    print("EXP 3: MINIMIZE TOTAL δstate IN THE GAP")
    print("="*70)

    n_basis = len(basis)
    W_base = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    states_base, _ = sha_states(W_base)

    def gap_score(dW):
        W_mod = [(W_base[j]^dW[j]) for j in range(16)]
        states_mod, _ = sha_states(W_mod)
        return sum(sum(hw(states_base[r][i]^states_mod[r][i]) for i in range(8))
                   for r in range(17, 52))

    # Also measure final hash diff for comparison
    def hash_score(dW):
        W_mod = [(W_base[j]^dW[j]) for j in range(16)]
        states_mod, _ = sha_states(W_mod)
        return sum(hw(add32(IV[i],states_base[64][i]) ^ add32(IV[i],states_mod[64][i]))
                   for i in range(8))

    import math
    cur_bits = list(basis[0])
    cur_score = gap_score(bits_to_words(cur_bits))
    best_score = cur_score
    best_bits = list(cur_bits)
    best_hash = hash_score(bits_to_words(cur_bits))

    for it in range(N_search):
        trial = list(cur_bits)
        idx = int.from_bytes(os.urandom(2),'big') % n_basis
        for j in range(512): trial[j] ^= basis[idx][j]
        if all(b==0 for b in trial): continue

        sc = gap_score(bits_to_words(trial))
        T = max(0.01, 1.0-it/N_search)
        if sc < cur_score or \
           math.exp(-(sc-cur_score)/(T*10)) > int.from_bytes(os.urandom(4),'big')/(1<<32):
            cur_bits = trial; cur_score = sc
        if cur_score < best_score:
            best_score = cur_score
            best_bits = list(cur_bits)
            best_hash = hash_score(bits_to_words(cur_bits))

    random_gap = 35 * 128  # 35 rounds × 128 bits diff
    print(f"  Best gap total HW: {best_score} (random ≈ {random_gap})")
    print(f"  Reduction: {(random_gap-best_score)/random_gap*100:.1f}%")
    print(f"  Corresponding hash diff: HW = {best_hash}")

    # Compare: random (non-kernel) δW
    best_random = 999999
    best_rand_hash = 256
    for _ in range(min(N_search, 1000)):
        dW = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        sc = gap_score(dW)
        if sc < best_random:
            best_random = sc
            best_rand_hash = hash_score(dW)

    print(f"\n  Comparison (random δW, {min(N_search,1000)} tries):")
    print(f"    Best gap HW: {best_random}")
    print(f"    Best hash diff: {best_rand_hash}")
    print(f"    Kernel advantage (gap): {best_random - best_score:+d}")
    print(f"    Kernel advantage (hash): {best_rand_hash - best_hash:+d}")


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 500

    basis = exp1_wang_in_gap(min(N*10, 5000))
    exp2_simultaneous(basis, min(N*20, 10000))
    exp3_minimize_gap(basis, min(N*20, 10000))
