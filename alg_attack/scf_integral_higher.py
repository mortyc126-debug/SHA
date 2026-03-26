#!/usr/bin/env python3
"""
SCF: Higher-Order Integral — расширяем 5-раундовый отличитель.

Факт: k=8 даёт 14 balanced бит при R=5, 0 при R=6.
Вопрос: k=16,20,32 — больше раундов?

Теория: integral degree d нужен 2^d элементов.
SHA round = degree 2. Через R раундов: degree 2^R.
Для полной сатурации degree-d integral: нужно R > log2(d).
k=8: d=8, log2(8)=3, но мы видим 5 раундов (shift register bonus +2).
k=32: d=32, log2(32)=5, ожидаем ~7 раундов?

Плюс: integral по НЕСКОЛЬКИМ словам W[0..m] одновременно.
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

def sha_R(W16, R):
    W=expand_real(W16); s=list(IV)
    for r in range(R):
        a,b,c,d,e,f,g,h=s
        T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add32(Sig0(a),Maj(a,b,c))
        s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
    return s


# ============================================================
# EXP 1: Vary k (active bits in W[0])
# ============================================================
def exp1_vary_k(N_ctx):
    print("="*70)
    print("EXP 1: INTEGRAL vs k (active bits in W[0])")
    print("="*70)

    # k=8: 256 elements. k=16: 65536. k=20: 1M. k=32: 4G (too big).
    # For large k: use RANDOM SAMPLING with XOR-accumulation
    # Property: for truly balanced bits, XOR of random 2^k subset = 0
    # We approximate by XOR of N_sample random elements

    for k in [8, 12, 16, 20]:
        n_set = 1 << k
        if n_set > 2000000:
            print(f"\n  k={k}: 2^{k}={n_set} too large, using sampling")
            continue

        print(f"\n  k={k} (set size = 2^{k} = {n_set}):")
        print(f"  {'R':>3} | balanced bits (per register)")
        print("  " + "-"*60)

        for R in range(1, 12):
            always_zero_count = [0]*256

            for ctx in range(N_ctx):
                W_base = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
                xor_sum = [0]*8

                for v in range(n_set):
                    W = list(W_base)
                    W[0] = (W[0] & ~((1<<k)-1)) | v
                    state = sha_R(W, R)
                    for i in range(8):
                        xor_sum[i] ^= state[i]

                for reg in range(8):
                    for bit in range(32):
                        if (xor_sum[reg]>>bit)&1 == 0:
                            always_zero_count[reg*32+bit] += 1

            n_always = sum(1 for c in always_zero_count if c == N_ctx)
            per_reg = []
            for reg in range(8):
                per_reg.append(sum(1 for b in range(32)
                                  if always_zero_count[reg*32+b]==N_ctx))

            regs = "abcdefgh"
            reg_str = " ".join(f"{regs[i]}={per_reg[i]:2d}" for i in range(8))
            marker = " ★" if n_always > 0 else ""
            print(f"  {R:3d} | {n_always:3d}/256  {reg_str}{marker}")

            if n_always == 0:
                break  # No more balanced bits


# ============================================================
# EXP 2: Multi-word integral (vary bits across W[0] AND W[1])
# ============================================================
def exp2_multi_word(N_ctx):
    print("\n" + "="*70)
    print("EXP 2: MULTI-WORD INTEGRAL")
    print("  Active bits split across W[0] and W[1]")
    print("="*70)

    # k bits in W[0], k bits in W[1] → 2^{2k} elements
    for k_per_word in [4, 6, 8]:
        n_set = (1 << k_per_word) ** 2
        print(f"\n  {k_per_word} bits/word × 2 words (set = 2^{2*k_per_word} = {n_set}):")
        print(f"  {'R':>3} | balanced")
        print("  " + "-"*30)

        mask = (1 << k_per_word) - 1

        for R in range(1, 12):
            always_zero = [True]*256

            for ctx in range(N_ctx):
                W_base = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
                xor_sum = [0]*8

                for v0 in range(1 << k_per_word):
                    for v1 in range(1 << k_per_word):
                        W = list(W_base)
                        W[0] = (W[0] & ~mask) | v0
                        W[1] = (W[1] & ~mask) | v1
                        state = sha_R(W, R)
                        for i in range(8):
                            xor_sum[i] ^= state[i]

                for idx in range(256):
                    reg, bit = idx//32, idx%32
                    if (xor_sum[reg]>>bit)&1:
                        always_zero[idx] = False

            n_always = sum(always_zero)
            marker = " ★" if n_always > 0 else ""
            print(f"  {R:3d} | {n_always:3d}/256{marker}")

            if n_always == 0:
                break


# ============================================================
# EXP 3: Which SPECIFIC bits stay balanced longest?
# ============================================================
def exp3_survivor_bits(N_ctx):
    print("\n" + "="*70)
    print("EXP 3: SURVIVOR ANALYSIS — which bits balanced at R=5?")
    print("="*70)

    k = 8
    R = 5

    bit_balanced = [0]*256

    for ctx in range(N_ctx * 3):
        W_base = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        xor_sum = [0]*8

        for v in range(1 << k):
            W = list(W_base)
            W[0] = (W[0] & ~((1<<k)-1)) | v
            state = sha_R(W, R)
            for i in range(8):
                xor_sum[i] ^= state[i]

        for reg in range(8):
            for bit in range(32):
                if (xor_sum[reg]>>bit)&1 == 0:
                    bit_balanced[reg*32+bit] += 1

    total = N_ctx * 3
    print(f"\n  R={R}, k={k}, {total} contexts")
    print(f"  Bits balanced in ≥90% of contexts:")

    regs = "abcdefgh"
    for reg in range(8):
        bits_90 = []
        for bit in range(32):
            if bit_balanced[reg*32+bit] >= total * 0.9:
                bits_90.append(bit)
        if bits_90:
            print(f"    {regs[reg]}: bits {bits_90}")

    # Which register stays balanced longest on AVERAGE?
    print(f"\n  Average fraction balanced at R={R}:")
    for reg in range(8):
        avg = sum(bit_balanced[reg*32+b] for b in range(32)) / (32 * total)
        bar = "█" * int(avg * 40)
        print(f"    {regs[reg]}: {avg:.3f} {bar}")


# ============================================================
# EXP 4: Integral from the END (backward integral)
# ============================================================
def exp4_backward_integral(N_ctx):
    print("\n" + "="*70)
    print("EXP 4: BACKWARD INTEGRAL")
    print("  Fix output state, vary input → integral from the END")
    print("  Use free-start: vary IV bits instead of W")
    print("="*70)

    k = 8

    for R in range(1, 10):
        always_zero = [True]*256

        for ctx in range(N_ctx):
            W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
            iv_base = [int.from_bytes(os.urandom(4),'big') for _ in range(8)]

            xor_sum = [0]*8
            Wexp = expand_real(W)

            for v in range(1 << k):
                iv = list(iv_base)
                iv[0] = (iv[0] & ~((1<<k)-1)) | v  # Vary k bits of IV[0]

                s = list(iv)
                for r in range(R):
                    a,b,c,d,e,f,g,h = s
                    T1=add32(h,Sig1(e),Ch(e,f,g),K[r],Wexp[r])
                    T2=add32(Sig0(a),Maj(a,b,c))
                    s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]

                for i in range(8):
                    xor_sum[i] ^= s[i]

            for idx in range(256):
                reg, bit = idx//32, idx%32
                if (xor_sum[reg]>>bit)&1:
                    always_zero[idx] = False

        n_always = sum(always_zero)
        marker = " ★" if n_always > 0 else ""
        print(f"  R={R:2d}: {n_always:3d}/256 balanced{marker}")

        if n_always == 0 and R > 3:
            break


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10

    exp1_vary_k(N_ctx=N)
    exp2_multi_word(N_ctx=N)
    exp3_survivor_bits(N_ctx=N)
    exp4_backward_integral(N_ctx=N)

    print("\n" + "="*70)
    print("ИТОГ: HIGHER-ORDER INTEGRAL")
    print("="*70)
