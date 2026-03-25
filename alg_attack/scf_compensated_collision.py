#!/usr/bin/env python3
"""
SCF: COMPENSATED COLLISION — δstate ≠ 0, но δIV + δstate = 0.

Стандарт: collision = δstate[64] = 0 (256 бит constraints)
Новое:    collision = δIV + δstate[64] = 0 (256 бит, но ДРУГИЕ!)

В 2-block: δIV_block2 = δH_block1 (мы КОНТРОЛИРУЕМ).
           δstate_block2 зависит от δIV И от M2.

Если δIV = −δstate → δH = 0 (collision!)
Но δstate при δIV≠0 ≠ δstate при δIV=0.

Ключевой вопрос: какова размерность constraints в ЭТОЙ задаче?
Старое: 12 значений a-seq × 32 бит = 384 бит constraints
Новое:  8 hash registers × 32 бит = 256 бит (МЕНЬШЕ!)

512 − 256 = 256 свободных бит → O(2^{128})?
НО: multiblock добавляет свободу block1 (512 бит).
Total freedom: 1024 бит. Constraints: 256 бит.
Kernel: 1024 − 256 = 768 → O(2^{384})?

Нет — constraints НЕЛИНЕЙНЫ. Реальный kernel может быть меньше.
Измерим.
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

def sha_raw_state(W16, iv):
    """Return raw state[64] (before MD addition)."""
    W=expand_real(W16); s=list(iv)
    for r in range(64):
        a,b,c,d,e,f,g,h=s
        T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add32(Sig0(a),Maj(a,b,c))
        s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
    return s

def sha_hash(W16, iv):
    """Return hash = iv + state[64]."""
    s = sha_raw_state(W16, iv)
    return [add32(iv[i], s[i]) for i in range(8)]

def two_block_hash(M1, M2):
    H1 = sha_hash(M1, IV)
    return sha_hash(M2, H1)


# ============================================================
# EXP 1: Compensated collision — δstate ≠ 0 but δIV = −δstate
# ============================================================
def exp1_compensated(N):
    print("="*70)
    print("EXP 1: COMPENSATED COLLISION")
    print("  2-block: H = compress(compress(IV, M1), M2)")
    print("  Pair: (M1,M2) vs (M1',M2') where M1≠M1', M2=M2'")
    print("  Block 1: creates δH1 (near-collision)")
    print("  Block 2: SAME M2, different IV → δH2 = δIV + δstate2")
    print("  Want: δH2 = 0, i.e., δstate2 = −δIV")
    print("  This does NOT require δstate2=0!")
    print("="*70)

    # The key insight: in block 2, δIV ≠ 0.
    # state2 depends on IV nonlinearly.
    # δH2 = δIV + δstate2(IV, M2)
    # = δIV + [state2(IV+δIV, M2) - state2(IV, M2)]
    # ≈ δIV + J_IV · δIV (linear approx)
    # = (I + J_IV) · δIV
    # For δH2 = 0: (I + J_IV) · δIV = 0

    # This is a LINEAR system in δIV!
    # If rank(I + J_IV) < 256 → nontrivial δIV exists → collision!

    print(f"\n  Measuring rank(I + J_IV) over GF(2):")

    for trial in range(min(N, 5)):
        M2 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        iv_base = [int.from_bytes(os.urandom(4),'big') for _ in range(8)]

        H_base = sha_hash(M2, iv_base)

        # Jacobian J_IV: ∂H/∂IV (GF(2))
        # (I + J_IV)[bit_out][bit_in] = ∂(H[out] ⊕ IV[out]) / ∂IV[in]
        # = ∂state[out] / ∂IV[in]
        # H = IV + state → δH = δIV + δstate
        # So I+J_IV = I + ∂state/∂IV = ∂H/∂IV

        rows = []
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

        # Rank of ∂H/∂IV
        m = list(rows); rank = 0
        for bp in range(255, -1, -1):
            mask = 1 << bp
            pv = -1
            for i in range(rank, len(m)):
                if m[i] & mask: pv=i; break
            if pv == -1: continue
            m[rank],m[pv] = m[pv],m[rank]
            for i in range(len(m)):
                if i!=rank and m[i]&mask: m[i]^=m[rank]
            rank += 1

        kernel = 256 - rank
        print(f"    Trial {trial}: rank(∂H/∂IV) = {rank}, kernel = {kernel}")

        if kernel > 0:
            print(f"    ★★★ KERNEL > 0! {kernel} directions of IV-compensation!")


# ============================================================
# EXP 2: Actual 2-block compensated collision search
# ============================================================
def exp2_compensated_search(N):
    print("\n" + "="*70)
    print("EXP 2: COMPENSATED 2-BLOCK SEARCH")
    print("  Block 1: M1/M1' → H1/H1' with δH1 controlled")
    print("  Block 2: M2/M2' → want δ(H1+state2) = 0")
    print("  Strategy: optimize (M1, M2) jointly")
    print("="*70)

    # Strategy: optimize M1 delta to create δH1 that Block 2 can compensate
    # Block 2 with different IV: δH2 = δIV + δstate2
    # If J_IV has specific structure → choose δIV that self-cancels

    best_overall = 256

    for trial in range(N):
        M1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        M2 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

        # Screen M1 delta
        best_delta = (256, 0, 0)
        for word in range(4):
            for bit in range(32):
                M1p = list(M1); M1p[word] ^= (1<<bit)
                H1a = sha_hash(M1, IV); H1b = sha_hash(M1p, IV)
                # 2-block hash
                H2a = sha_hash(M2, H1a); H2b = sha_hash(M2, H1b)
                d = sum(hw(H2a[i]^H2b[i]) for i in range(8))
                if d < best_delta[0]:
                    best_delta = (d, word, bit)

        raw_hw, dw, db = best_delta
        M1p = list(M1); M1p[dw] ^= (1<<db)

        # Now optimize M2 with SA
        H1a = sha_hash(M1, IV); H1b = sha_hash(M1p, IV)

        def score(m2):
            Ha = sha_hash(m2, H1a); Hb = sha_hash(m2, H1b)
            return sum(hw(Ha[i]^Hb[i]) for i in range(8))

        best = score(M2); cur = list(M2)
        for it in range(1000):
            t = list(cur)
            w = int.from_bytes(os.urandom(1),'big') % 16
            b = int.from_bytes(os.urandom(1),'big') % 32
            t[w] ^= (1<<b)
            s = score(t)
            T = max(0.01, 1-it/1000)
            if s < best or math.exp(-(s-best)/(T*2)) > int.from_bytes(os.urandom(4),'big')/(1<<32):
                cur = t
                if s < best: best = s

        if best < best_overall:
            best_overall = best

        if trial < 5 or best < 90:
            print(f"  Trial {trial}: δH1_hw={sum(hw(H1a[i]^H1b[i]) for i in range(8))}, "
                  f"raw_2block={raw_hw}, optimized={best}")

    print(f"\n  Best compensated 2-block: HW = {best_overall}")


# ============================================================
# EXP 3: Dimension of COMPENSATED collision variety
# ============================================================
def exp3_compensated_dimension():
    print("\n" + "="*70)
    print("EXP 3: DIMENSION OF COMPENSATED COLLISION VARIETY")
    print("  Standard: 512 freedom - 384 constraints = 128 kernel")
    print("  Compensated: HOW MANY constraints?")
    print("="*70)

    # In compensated 2-block:
    # Variables: M1[0..15] + M2[0..15] = 1024 bits
    # Function: F(M1,M2) = sha(sha(IV,M1), M2)
    # Collision: F(M1,M2) = F(M1',M2') where (M1,M2)≠(M1',M2')
    # This is a 1024→256 map. Expected kernel: 1024-256 = 768.

    # But for the SPECIFIC case M2=M2' (same block 2):
    # Variables: M1[0..15] = 512 bits
    # δIV = sha(IV,M1) ⊕ sha(IV,M1')
    # Output: sha(δIV+IV, M2) ⊕ sha(IV, M2) = 0
    # This is 512 → 256. Kernel: 512-256 = 256? Or less?

    M2 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    M1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

    # Jacobian of the COMPOSED map:
    # F(M1) = sha(sha(IV,M1), M2)
    # ∂F/∂M1 = (∂sha/∂IV)|_{IV=H1} · (∂H1/∂M1)
    # Size: 256 × 512. Rank ≤ 256.

    H1_base = sha_hash(M1, IV)
    F_base = sha_hash(M2, H1_base)

    rows = []
    for out_reg in range(8):
        for out_bit in range(32):
            row = 0
            for word in range(16):
                for bit in range(32):
                    M1p = list(M1); M1p[word] ^= (1<<bit)
                    H1p = sha_hash(M1p, IV)
                    Fp = sha_hash(M2, H1p)
                    if (F_base[out_reg] ^ Fp[out_reg]) >> out_bit & 1:
                        row |= 1 << (word*32 + bit)
            rows.append(row)

    m = list(rows); rank = 0
    for bp in range(511, -1, -1):
        mask = 1 << bp
        pv = -1
        for i in range(rank, 256):
            if m[i] & mask: pv=i; break
        if pv == -1: continue
        m[rank],m[pv] = m[pv],m[rank]
        for i in range(256):
            if i!=rank and m[i]&mask: m[i]^=m[rank]
        rank += 1

    kernel = 512 - rank
    print(f"  Jacobian ∂F_2block/∂M1 (M2 fixed):")
    print(f"    Size: 256 × 512")
    print(f"    Rank: {rank}")
    print(f"    Kernel: {kernel}")
    print(f"    Birthday in kernel: O(2^{kernel//2})")

    # Compare with standard single-block
    print(f"\n  Standard single-block: kernel=128 → O(2^64)")
    print(f"  Compensated 2-block (M2 fixed): kernel={kernel} → O(2^{kernel//2})")

    if kernel > 128:
        print(f"  ★ LARGER KERNEL! {kernel-128} extra dimensions!")
    elif kernel == 128:
        print(f"  Same kernel — compensation doesn't help for fixed M2")
    else:
        print(f"  Smaller kernel — something lost")


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10

    exp1_compensated(N)
    exp3_compensated_dimension()
    exp2_compensated_search(min(N, 10))
