#!/usr/bin/env python3
"""
SCF: Integral + Nonlinear Invariant Attack

НОВОЕ НАПРАВЛЕНИЕ: атакуем не компоненты, а СТРУКТУРУ СВЯЗЕЙ.

1. INTEGRAL ATTACK (сатурация):
   Берём 2^k сообщений где k бит = все значения.
   XOR-сумма выходов после R раундов.
   Если XOR-sum имеет нулевые биты → отличитель.

2. NONLINEAR INVARIANT:
   Ищем g: state → bit, такую что g(round(state)) = g(state).
   Если найдём → collision в подпространстве dim/2.

3. PARTIAL INTEGRAL:
   Не полная сатурация, а отдельные регистры.
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

def sha_R_rounds(W16, R, iv=None):
    if iv is None: iv=IV
    W=expand_real(W16); s=list(iv)
    for r in range(R):
        a,b,c,d,e,f,g,h=s
        T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add32(Sig0(a),Maj(a,b,c))
        s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
    return s


# ============================================================
# EXP 1: INTEGRAL DISTINGUISHER
# ============================================================
def exp1_integral(N_contexts=10):
    """
    For k active bits in W[0], compute XOR-sum of states after R rounds.
    If any output bit of XOR-sum = 0 for ALL contexts → integral property.
    Random function: each bit of XOR-sum is 0 with P=0.5.
    """
    print("="*70)
    print("EXP 1: INTEGRAL (SATURATION) DISTINGUISHER")
    print("  Set 2^k messages with k bits of W[0] taking all values")
    print("  XOR-sum of states after R rounds")
    print("="*70)

    for k in [8, 10]:
        n_set = 1 << k

        for R in [1, 2, 3, 4, 5, 6, 8, 16, 32, 64]:
            # Count how many output bits are ALWAYS zero across contexts
            always_zero = [True] * 256  # 8 registers × 32 bits

            for ctx in range(N_contexts):
                # Random context: fix W[1..15], vary W[0] low k bits
                W_base = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

                xor_sum = [0] * 8
                for v in range(n_set):
                    W = list(W_base)
                    W[0] = (W[0] & ~((1 << k) - 1)) | v  # Set low k bits
                    state = sha_R_rounds(W, R)
                    for i in range(8):
                        xor_sum[i] ^= state[i]

                # Check which bits are zero
                for reg in range(8):
                    for bit in range(32):
                        if (xor_sum[reg] >> bit) & 1:
                            always_zero[reg * 32 + bit] = False

            n_zero = sum(always_zero)
            n_expected = 0  # Random: P(always 0 in N contexts) = 2^{-N}

            marker = ""
            if n_zero > 0:
                marker = f" ★★★ {n_zero} BALANCED BITS!"

            if R <= 8 or n_zero > 0:
                print(f"  k={k:2d}, R={R:2d}: {n_zero:3d}/256 always-zero bits{marker}")


# ============================================================
# EXP 2: NONLINEAR INVARIANT SEARCH
# ============================================================
def exp2_nonlinear_invariant(N):
    """
    Search for g(state) ∈ GF(2) such that g(round(state)) = g(state).

    Test candidates:
    - Single bit: g = state[reg][bit]
    - XOR of 2 bits: g = state[r1][b1] ⊕ state[r2][b2]
    - Parity of register: g = ⊕_i state[reg][i]
    - Custom: a[0]⊕e[0], parity(a)⊕parity(e), etc.
    """
    print("\n" + "="*70)
    print("EXP 2: NONLINEAR INVARIANT SEARCH")
    print("  Find g(state) preserved by one round of SHA")
    print("="*70)

    # Generate random states and their round images
    pairs = []
    for _ in range(N):
        W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        state = [int.from_bytes(os.urandom(4),'big') for _ in range(8)]
        next_state = sha_R_rounds(W, 1, iv=state)
        pairs.append((state, next_state, W))

    # Test type 1: single bit invariant
    print(f"\n  Type 1: Single bit g = state[reg][bit]")
    for reg in range(8):
        for bit in [0, 1, 15, 16, 31]:
            preserved = sum(1 for s, ns, _ in pairs
                          if ((s[reg]>>bit)&1) == ((ns[reg]>>bit)&1))
            if abs(preserved/N - 0.5) > 0.1:
                print(f"    g=state[{reg}][{bit}]: preserved {preserved}/{N} = {preserved/N:.3f}")

    # Test type 2: parity of register
    print(f"\n  Type 2: Parity g = ⊕ state[reg][0..31]")
    for reg in range(8):
        preserved = sum(1 for s, ns, _ in pairs
                       if hw(s[reg]) % 2 == hw(ns[reg]) % 2)
        bias = abs(preserved/N - 0.5)
        marker = " ★" if bias > 0.05 else ""
        print(f"    parity(reg[{reg}]): {preserved}/{N} = {preserved/N:.3f}{marker}")

    # Test type 3: XOR of two register parities
    print(f"\n  Type 3: g = parity(reg_i) ⊕ parity(reg_j)")
    best_bias = 0
    best_pair = None
    for r1 in range(8):
        for r2 in range(r1+1, 8):
            preserved = sum(1 for s, ns, _ in pairs
                           if (hw(s[r1])^hw(s[r2])) % 2 == (hw(ns[r1])^hw(ns[r2])) % 2)
            bias = abs(preserved/N - 0.5)
            if bias > best_bias:
                best_bias = bias
                best_pair = (r1, r2)
    if best_pair:
        print(f"    Best: parity(reg[{best_pair[0]}])⊕parity(reg[{best_pair[1]}])")
        print(f"    Bias: {best_bias:.4f} ({'★ SIGNIFICANT' if best_bias > 0.05 else 'noise'})")

    # Test type 4: LSB combinations
    print(f"\n  Type 4: g = state[r1][0] ⊕ state[r2][0] (LSB)")
    for r1 in range(8):
        for r2 in range(r1+1, 8):
            preserved = sum(1 for s, ns, _ in pairs
                           if ((s[r1]^s[r2])&1) == ((ns[r1]^ns[r2])&1))
            bias = abs(preserved/N - 0.5)
            if bias > 0.03:
                print(f"    LSB({r1})⊕LSB({r2}): {preserved/N:.3f} (bias={bias:.3f})")

    # Test type 5: shift register invariant
    # b=a[-1] means: after round, b_new = a_old
    # So g = a[bit] should equal g' = b[bit] after one round
    print(f"\n  Type 5: Shift register consistency")
    for bit in [0, 1, 15, 31]:
        # g(state) = a[bit], g(next) = b[bit] (since b_new = a_old)
        match = sum(1 for s, ns, _ in pairs
                   if ((s[0]>>bit)&1) == ((ns[1]>>bit)&1))
        print(f"    a[{bit}] → b[{bit}]: {match}/{N} = {match/N:.4f} (expect 1.0)")


# ============================================================
# EXP 3: MULTI-ROUND INTEGRAL — deeper analysis
# ============================================================
def exp3_deep_integral(N_ctx=20):
    """
    More detailed integral: track which SPECIFIC output bits are balanced.
    Also test: XOR-SHA vs Real SHA (does carry break integral property?)
    """
    print("\n" + "="*70)
    print("EXP 3: DEEP INTEGRAL — which bits stay balanced longest?")
    print("="*70)

    k = 8  # 2^8 = 256 elements in the set

    # For each round R, count balanced bits across N contexts
    print(f"  k={k} active bits in W[0]")
    print(f"  {'R':>3} | {'balanced':>9} | {'per_register':>50}")
    print("  " + "-"*70)

    for R in range(1, 17):
        balanced_count = [0] * 256

        for ctx in range(N_ctx):
            W_base = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
            xor_sum = [0] * 8

            for v in range(1 << k):
                W = list(W_base)
                W[0] = (W[0] & ~((1 << k) - 1)) | v
                state = sha_R_rounds(W, R)
                for i in range(8):
                    xor_sum[i] ^= state[i]

            for reg in range(8):
                for bit in range(32):
                    if (xor_sum[reg] >> bit) & 1 == 0:
                        balanced_count[reg*32 + bit] += 1

        # A bit is "balanced" if zero in ALL contexts
        n_always = sum(1 for c in balanced_count if c == N_ctx)

        # Per-register breakdown
        per_reg = []
        for reg in range(8):
            reg_always = sum(1 for b in range(32)
                           if balanced_count[reg*32+b] == N_ctx)
            per_reg.append(reg_always)

        regs = "abcdefgh"
        reg_str = " ".join(f"{regs[i]}={per_reg[i]:2d}" for i in range(8))

        marker = " ★" if n_always > 0 else ""
        print(f"  {R:3d} | {n_always:9d} | {reg_str}{marker}")


# ============================================================
# EXP 4: INVARIANT SUBSPACE — does round preserve any subspace?
# ============================================================
def exp4_invariant_subspace(N):
    """
    Test: for various subspaces V, is round(V) ⊆ V?
    Subspace candidates: fix some registers to zero.
    """
    print("\n" + "="*70)
    print("EXP 4: INVARIANT SUBSPACE SEARCH")
    print("  Does round(V) ⊆ V for any nontrivial V?")
    print("="*70)

    # Test: states where specific registers = 0
    print(f"\n  Test: fix register(s) to 0, check if stays 0 after round")

    for fixed_regs in [[3], [7], [3,7], [1,2,3], [5,6,7], [1,2,3,5,6,7]]:
        preserved = 0
        total = 0

        for _ in range(N):
            W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
            state = [int.from_bytes(os.urandom(4),'big') for _ in range(8)]
            for r in fixed_regs:
                state[r] = 0

            next_state = sha_R_rounds(W, 1, iv=state)
            total += 1

            if all(next_state[r] == 0 for r in fixed_regs):
                preserved += 1

        regs_str = ",".join(str(r) for r in fixed_regs)
        print(f"    Fix regs [{regs_str}]=0: preserved {preserved}/{total}")

    # Test: states where a = e (mirror symmetry)
    print(f"\n  Test: a=e symmetry")
    preserved = 0
    for _ in range(N):
        W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        val = int.from_bytes(os.urandom(4),'big')
        state = [val,0,0,0,val,0,0,0]  # a=e, rest=0
        ns = sha_R_rounds(W, 1, iv=state)
        if ns[0] == ns[4]:  # a_new == e_new?
            preserved += 1
    print(f"    a=e preserved: {preserved}/{N}")

    # Test: a⊕e = constant
    print(f"\n  Test: a⊕e = const after round?")
    for const in [0, 0xFFFFFFFF]:
        count = 0
        for _ in range(N):
            W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
            state = [int.from_bytes(os.urandom(4),'big') for _ in range(8)]
            state[4] = state[0] ^ const  # e = a ⊕ const
            ns = sha_R_rounds(W, 1, iv=state)
            if ns[0] ^ ns[4] == const:
                count += 1
        print(f"    a⊕e=0x{const:08x}: preserved {count}/{N}")


# ============================================================
if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 500

    exp1_integral(N_contexts=min(N//10+1, 15))
    exp3_deep_integral(N_ctx=min(N//20+1, 15))
    exp2_nonlinear_invariant(min(N, 2000))
    exp4_invariant_subspace(min(N, 1000))

    print("\n" + "="*70)
    print("ИТОГ")
    print("="*70)
    print("""
  INTEGRAL: если balanced bits > 0 при R > 5 → отличитель!
  INVARIANT: если bias > 0.05 → неслучайная структура!
  SUBSPACE: если preserved > 0 → инвариантное подпространство!

  Любой из трёх → потенциальный shortcut для collision.
    """)
