#!/usr/bin/env python3
"""
SCF: Newton multi-round — FIXED (no trivial solutions).

Баг: предыдущий код использовал basis[0] как ref, и solve
мог отменить ref → δW=0. Фикс: solve без ref, проверка δW≠0.

Вопрос: сколько раундов gap можно ОДНОВРЕМЕННО обнулить
с НЕТРИВИАЛЬНЫМ δW из ядра?
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


def compute_jacobian_direct(basis, W_base, target_rounds, target_reg=4):
    """Compute Jacobian: for each basis vector, what is δ(register) at each target round?
    NO reference point — directly measure each basis vector's effect.
    Returns J (n_constraints × n_basis) and target (n_constraints bits)."""
    n_basis = len(basis)
    basis_words = [bits_to_words(v) for v in basis]

    # For each basis vector i, compute SHA(W_base ⊕ basis_words[i])
    # and measure δ(reg) at each target round
    state_base_cache = {}
    for tr in target_rounds:
        state_base_cache[tr] = sha_state_at_r(W_base, tr+1)

    # Column i of J = bits of δ(reg)[target_rounds] for basis vector i
    # Row j of J = how bit j depends on each basis vector
    # Target = the actual δ values (what we want to zero)

    # First: compute δ for each basis vector individually
    delta_per_basis = []  # [basis_idx] → list of (round, δ_value)
    for i in range(n_basis):
        W_mod = [(W_base[j] ^ basis_words[i][j]) for j in range(16)]
        deltas = []
        for tr in target_rounds:
            s_mod = sha_state_at_r(W_mod, tr+1)
            d = state_base_cache[tr][target_reg] ^ s_mod[target_reg]
            deltas.append(d)
        delta_per_basis.append(deltas)

    # Build J: each row = one output bit, each col = one basis vector
    # J[row][col] = XOR sensitivity
    # For GF(2) linearization: J[bit][i] = δ_i[bit]
    # (first-order approximation: flipping basis[i] flips these output bits)
    n_rounds = len(target_rounds)
    J = []
    target_bits = []

    for r_idx in range(n_rounds):
        for bit in range(32):
            row = []
            for i in range(n_basis):
                row.append((delta_per_basis[i][r_idx] >> bit) & 1)
            J.append(row)
            target_bits.append(0)  # We want δ = 0, so target = current δ

    # Target: for the "zero" kernel element (δW=0), δ=0.
    # But we need NONTRIVIAL solution. The system J·c = 0 with c≠0.
    # This is finding the KERNEL of J (not solving J·c = target).

    return J, delta_per_basis


def gf2_kernel(M):
    """Find kernel of matrix M over GF(2). Returns list of kernel vectors."""
    if not M: return []
    n = len(M); m = len(M[0])
    aug = [list(row) for row in M]
    pivots = []; ri = 0
    for col in range(m):
        pv = -1
        for r in range(ri, n):
            if aug[r][col]: pv=r; break
        if pv == -1: continue
        aug[ri], aug[pv] = aug[pv], aug[ri]
        for r in range(n):
            if r != ri and aug[r][col]:
                for c in range(m): aug[r][c] ^= aug[ri][c]
        pivots.append(col); ri += 1

    rank = len(pivots)
    free = [c for c in range(m) if c not in pivots]

    ker = []
    for fc in free:
        v = [0]*m; v[fc] = 1
        for i, pc in enumerate(pivots):
            if i < len(aug) and aug[i][fc]:
                v[pc] = 1
        ker.append(v)

    return ker, rank


# ============================================================
# EXP 1: For N simultaneous rounds, find kernel of Jacobian
# ============================================================
def exp1_simultaneous(basis, W_base):
    print("="*70)
    print("EXP 1: SIMULTANEOUS NEWTON — kernel of Jacobian")
    print("  Find c ≠ 0 in schedule kernel s.t. δe[r]=0 at multiple rounds")
    print("="*70)

    n_basis = len(basis)

    for n_rounds in [1, 2, 3, 4, 6, 8]:
        step = max(1, 34 // n_rounds)
        rounds = [17 + i*step for i in range(n_rounds) if 17+i*step < 52][:n_rounds]

        J, delta_per_basis = compute_jacobian_direct(basis, W_base, rounds)

        ker, rank = gf2_kernel(J)
        ker_dim = len(ker)

        print(f"\n  {n_rounds} rounds {rounds}:")
        print(f"    J size: {len(J)}×{n_basis}, rank={rank}, kernel dim={ker_dim}")

        if ker_dim > 0:
            # Find nontrivial kernel element (c ≠ 0 AND resulting δW ≠ 0)
            best_hw_hash = 256
            best_dW = None

            for ki in range(min(ker_dim, 50)):
                # Try each kernel vector and combinations
                for mask in range(1, min(1 << min(ker_dim, 10), 1024)):
                    c = [0] * n_basis
                    for k in range(min(ker_dim, 10)):
                        if (mask >> k) & 1:
                            for i in range(n_basis):
                                c[i] ^= ker[k][i]

                    # Build δW
                    bits = [0]*512
                    for i in range(n_basis):
                        if c[i]:
                            for j in range(512): bits[j] ^= basis[i][j]

                    dW = bits_to_words(bits)
                    if all(w == 0 for w in dW):
                        continue  # Skip trivial

                    W_mod = [(W_base[j]^dW[j]) for j in range(16)]

                    # Verify δe at target rounds
                    all_zero = True
                    for tr in rounds:
                        s1 = sha_state_at_r(W_base, tr+1)
                        s2 = sha_state_at_r(W_mod, tr+1)
                        if s1[4] != s2[4]:
                            all_zero = False
                            break

                    if all_zero:
                        # Compute hash diff
                        H1 = sha_compress_full(W_base)
                        H2 = sha_compress_full(W_mod)
                        hash_hw = sum(hw(H1[i]^H2[i]) for i in range(8))
                        n_diff = sum(1 for i in range(16) if dW[i] != 0)

                        if hash_hw < best_hw_hash:
                            best_hw_hash = hash_hw
                            best_dW = dW

                        if hash_hw < 128:
                            print(f"    ★ VERIFIED: δe=0 at {n_rounds} rounds, "
                                  f"δW≠0 ({n_diff} words), hash HW={hash_hw}")
                            break

                    if all_zero:
                        break

                if best_dW and best_hw_hash < 128:
                    break

            if best_dW is None:
                print(f"    No nontrivial solution found in kernel")
            elif best_hw_hash >= 128:
                print(f"    Nontrivial solutions exist but hash HW={best_hw_hash}")
        else:
            print(f"    Kernel empty — system overdetermined, no solution")


# ============================================================
# EXP 2: Sweep W_base — different bases may have different ranks
# ============================================================
def exp2_sweep_bases(basis, N_bases):
    print("\n" + "="*70)
    print("EXP 2: SWEEP BASES — rank varies by W_base?")
    print("="*70)

    n_basis = len(basis)
    rounds_4 = [17, 25, 33, 41]

    rank_dist = {}
    best_ker_dim = 0

    for bi in range(N_bases):
        W_base = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        J, _ = compute_jacobian_direct(basis, W_base, rounds_4)
        _, rank = gf2_kernel(J)
        ker_dim = n_basis - rank

        rank_dist[rank] = rank_dist.get(rank, 0) + 1

        if ker_dim > best_ker_dim:
            best_ker_dim = ker_dim
            print(f"  [{bi}] rank={rank}, kernel={ker_dim} ★ NEW BEST")

    print(f"\n  Rank distribution (4 rounds, {N_bases} bases):")
    for r in sorted(rank_dist):
        ker = n_basis - r
        print(f"    rank={r}, kernel={ker}: {rank_dist[r]} bases")

    print(f"  Best kernel dim: {best_ker_dim}")
    print(f"  (Need kernel > 0 for nontrivial solution)")


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 50

    print("="*70)
    print("SCF: NEWTON MULTI-ROUND — FIXED")
    print("="*70)

    print("\nKernel R=52 (dim=128):")
    basis52 = extract_kernel_basis(52)

    W_base = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    exp1_simultaneous(basis52, W_base)
    exp2_sweep_bases(basis52, min(N, 30))

    print("\n\nKernel R=56 (dim=256):")
    basis56 = extract_kernel_basis(56)
    exp1_simultaneous(basis56, W_base)
    exp2_sweep_bases(basis56, min(N, 30))
