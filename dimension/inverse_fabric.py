"""
ОБРАТНАЯ ТКАНЬ R⁻¹: самостоятельный объект.

R: state + W → state' (прямая)
R⁻¹: state' + W → state (обратная, 100% exact)

В нашем измерении: R⁻¹ = ткань, прочитанная СПРАВА НАЛЕВО.
Та же структура, но другое НАПРАВЛЕНИЕ.

Вопросы:
  1. GF(2)-kernel R⁻¹ = GF(2)-kernel R? (или другой?)
  2. rank(CE⁻¹) = rank(CE)? (carry error при обратном проходе)
  3. A-repair работает на R⁻¹? (backward a-repair)
  4. Backward SHA-256 = другая "хеш-функция" с другими свойствами?

Идея: если rank(CE⁻¹) < rank(CE) = 256, то ОБРАТНАЯ ткань СЛАБЕЕ.
Collision в обратной ткани → collision в прямой (через инверсию).
"""

import numpy as np
import struct, hashlib

MASK32 = 0xFFFFFFFF
K_const = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
    0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
    0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
    0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
    0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
    0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
]
IV = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]

def hw(x): return bin(x).count('1')
def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def add32(x, y): return (x + y) & MASK32
def sub32(x, y): return (x - y) & MASK32
def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)

def R(state, W_r, r_idx):
    a,b,c,d,e,f,g,h = state
    T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K_const[r_idx]), W_r)
    T2 = add32(Sigma0(a), Maj(a,b,c))
    return (add32(T1,T2),a,b,c,add32(d,T1),e,f,g)

def R_inv(sn, W_r, r_idx):
    ap,bp,cp,dp,ep,fp,gp,hp = sn
    a,b,c,e,f,g = bp,cp,dp,fp,gp,hp
    T2 = add32(Sigma0(a), Maj(a,b,c))
    T1 = sub32(ap, T2)
    d = sub32(ep, T1)
    h = sub32(sub32(sub32(sub32(T1, Sigma1(e)), Ch(e,f,g)), K_const[r_idx]), W_r)
    return (a,b,c,d,e,f,g,h)

def expand_schedule(W16):
    W = list(W16)
    for r in range(16, 64):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))
    return W

def sha256_words(W16):
    return struct.unpack('>8I', hashlib.sha256(struct.pack('>16I', *W16)).digest())


def backward_sha256(state64, W_full):
    """SHA-256 backward: state[64] → state[0] using R⁻¹."""
    s = state64
    for r in range(63, -1, -1):
        s = R_inv(s, W_full[r], r)
    return s


def main():
    np.random.seed(42)

    print("=" * 70)
    print("ОБРАТНАЯ ТКАНЬ R⁻¹: самостоятельный объект")
    print("=" * 70)

    # =================================================================
    print(f"\n{'=' * 70}")
    print("1. BACKWARD HASH: F⁻¹(state64, W) → state0")
    print("=" * 70)

    # Forward: SHA-256(W) = H = state[64] + IV
    # state[64] = H - IV
    # Backward: R⁻¹ chain from state[64] to state[0]
    # If state[0] = IV → collision check.

    # Backward "hash": given target state[64], compute state[0]
    # This is a DIFFERENT function: B(target, W) = state[0]

    W_base = [np.random.randint(0, 2**32) for _ in range(16)]
    W_full = expand_schedule(W_base)

    # Forward
    s = tuple(IV)
    for r in range(64):
        s = R(s, W_full[r], r)
    state64 = s
    H = tuple(add32(IV[i], state64[i]) for i in range(8))

    # Backward from state64
    s_back = backward_sha256(state64, W_full)

    print(f"  Forward: IV → state64 → H")
    print(f"  Backward: state64 → state0")
    print(f"  state0 == IV: {s_back == tuple(IV)}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("2. GF(2)-KERNEL обратной ткани")
    print("=" * 70)

    # Backward "hash": B(state64, W) = state0
    # Фиксируем state64, варьируем W → разные state0.
    # T⁻¹ matrix: 512 input bits (W) → 256 output bits (state0)

    # Используем ОДНО state64 (from random W_base)
    T_inv_rows = []
    for word in range(16):
        for bit in range(32):
            W_mod = list(W_base)
            W_mod[word] ^= (1 << bit)
            W_mod_full = expand_schedule(W_mod)

            s0_mod = backward_sha256(state64, W_mod_full)

            row = []
            for w in range(8):
                delta = s_back[w] ^ s0_mod[w]
                for b in range(32):
                    row.append((delta >> b) & 1)
            T_inv_rows.append(row)

    T_inv = np.array(T_inv_rows, dtype=np.uint8)  # 512 × 256
    rank_T_inv = np.linalg.matrix_rank(T_inv.astype(float))

    print(f"\n  Backward T⁻¹ matrix: 512 × 256")
    print(f"  rank(T⁻¹) = {rank_T_inv}")
    print(f"  Kernel dim = {512 - rank_T_inv}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("3. CARRY ERROR RANK обратной ткани")
    print("=" * 70)

    if rank_T_inv == 256:
        # GF(2) kernel basis
        M = T_inv.T.copy()
        pivots = []
        row = 0
        for col in range(512):
            found = False
            for r in range(row, 256):
                if M[r, col] == 1:
                    M[[row, r]] = M[[r, row]]
                    found = True
                    break
            if not found:
                continue
            pivots.append(col)
            for r in range(256):
                if r != row and M[r, col] == 1:
                    M[r] = M[r] ^ M[row]
            row += 1

        free_vars = [c for c in range(512) if c not in pivots]
        print(f"  GF(2) kernel: {len(free_vars)} free vars")

        # Compute CE⁻¹ matrix
        basis_ces = []
        for fc in free_vars[:256]:
            x = np.zeros(512, dtype=np.uint8)
            x[fc] = 1
            for i in range(len(pivots)-1, -1, -1):
                pc = pivots[i]
                val = np.uint8(0)
                for j in range(512):
                    if j != pc:
                        val ^= (M[i, j] & x[j])
                x[pc] = val

            dW = [0]*16
            for word in range(16):
                for bit in range(32):
                    if x[word*32+bit]:
                        dW[word] ^= (1 << bit)

            W_mod_full = expand_schedule(dW)  # Wrong! Need W_base ^ dW
            # Actually: we need backward with W_base modified
            W_mod = [W_base[i] ^ dW[i] for i in range(16)]
            W_mod_full = expand_schedule(W_mod)
            s0_mod = backward_sha256(state64, W_mod_full)

            ce_bits = []
            for w in range(8):
                delta = s_back[w] ^ s0_mod[w]
                for b in range(32):
                    ce_bits.append((delta >> b) & 1)
            basis_ces.append(ce_bits)

        CE_inv = np.array(basis_ces[:256], dtype=np.uint8)
        rank_CE_inv = np.linalg.matrix_rank(CE_inv.astype(float))

        print(f"\n  ★ rank(CE⁻¹) = {rank_CE_inv}")
        print(f"  ★ rank(CE forward) = 256")

        if rank_CE_inv < 256:
            print(f"\n  !!! CE⁻¹ HAS LOWER RANK !!!")
            print(f"  ker(CE⁻¹) = {256 - rank_CE_inv}-dimensional")
            print(f"  → BACKWARD COLLISION easier than forward!")
        else:
            print(f"  CE⁻¹ = full rank (same as forward)")

        # CE carry error HW distribution
        ce_hws = [sum(hw(int.from_bytes(bytes(ce[w*8:(w+1)*8]), 'big')) for w in range(32)) for ce in basis_ces[:50]]
        # Simpler:
        ce_hws = []
        for ce in basis_ces[:100]:
            h = sum(ce)
            ce_hws.append(h)

        print(f"\n  CE⁻¹ carry error HW: mean={np.mean(ce_hws):.1f} (forward CE: 128)")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("4. АСИММЕТРИЯ R vs R⁻¹: side-by-side")
    print("=" * 70)

    # Forward: rank(T)=256, rank(CE)=256
    # Backward: rank(T⁻¹)=?, rank(CE⁻¹)=?

    # Also: linearity test for backward
    lin_diffs = []
    for _ in range(200):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        W_full = expand_schedule(W)

        s = tuple(IV)
        for r in range(64):
            s = R(s, W_full[r], r)
        s64 = s

        s0_base = backward_sha256(s64, W_full)

        w1, b1 = np.random.randint(0,16), np.random.randint(0,32)
        w2, b2 = np.random.randint(0,16), np.random.randint(0,32)

        dW1 = [0]*16; dW1[w1] ^= (1<<b1)
        dW2 = [0]*16; dW2[w2] ^= (1<<b2)
        dW12 = [dW1[i]^dW2[i] for i in range(16)]

        s0_1 = backward_sha256(s64, expand_schedule([W[i]^dW1[i] for i in range(16)]))
        s0_2 = backward_sha256(s64, expand_schedule([W[i]^dW2[i] for i in range(16)]))
        s0_12 = backward_sha256(s64, expand_schedule([W[i]^dW12[i] for i in range(16)]))

        T_1 = tuple(s0_base[i] ^ s0_1[i] for i in range(8))
        T_2 = tuple(s0_base[i] ^ s0_2[i] for i in range(8))
        T_12 = tuple(s0_base[i] ^ s0_12[i] for i in range(8))
        T_sum = tuple(T_1[i] ^ T_2[i] for i in range(8))

        diff = sum(hw(T_12[i] ^ T_sum[i]) for i in range(8))
        lin_diffs.append(diff)

    print(f"""
  FORWARD (R):
    rank(T) = 256
    rank(CE) = 256
    Linearity: |T(a⊕b) - T(a)⊕T(b)| = 128 (fully nonlinear)

  BACKWARD (R⁻¹):
    rank(T⁻¹) = {rank_T_inv}
    rank(CE⁻¹) = {rank_CE_inv if rank_T_inv == 256 else 'N/A'}
    Linearity: |T⁻¹(a⊕b) - T⁻¹(a)⊕T⁻¹(b)| = {np.mean(lin_diffs):.1f}

  ASYMMETRY: {'IDENTICAL' if rank_T_inv == 256 and rank_CE_inv == 256 else 'DIFFERENT!'}
""")


if __name__ == "__main__":
    main()
