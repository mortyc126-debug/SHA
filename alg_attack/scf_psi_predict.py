#!/usr/bin/env python3
"""
SCF: Ψ-PREDICTION — копаем в 1.4 бита предсказуемости.

EXP4 показал: avg HW(δΨ_total) = 126.6 при perturbation W[word]+=1.
Random = 128. Зазор 1.4 бита = Ψ_total НЕ полностью random.

Вопросы:
1. КАКИЕ биты Ψ_total более предсказуемы?
2. Зависит ли предсказуемость от СЛОВА (word)?
3. Зависит ли от НАПРАВЛЕНИЯ (бит W vs additive)?
4. Можно ли усилить предсказуемость через выбор δW?
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

def sha_hash(W16):
    W=expand_real(W16); s=list(IV)
    for r in range(64):
        a,b,c,d,e,f,g,h=s
        T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add32(Sig0(a),Maj(a,b,c))
        s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
    return [add32(IV[i],s[i]) for i in range(8)]

def sha_xor(W16):
    W=list(W16)
    for i in range(16,64):
        W.append(sig1(W[i-2])^W[i-7]^sig0(W[i-15])^W[i-16])
    s=list(IV)
    for r in range(64):
        a,b,c,d,e,f,g,h=s
        T1=h^Sig1(e)^Ch(e,f,g)^K[r]^W[r]
        T2=Sig0(a)^Maj(a,b,c)
        s=[T1^T2,a,b,c,d^T1,e,f,g]
    return [IV[i]^s[i] for i in range(8)]

def psi_total(W):
    """Ψ_total(W) = SHA_real(W) ⊕ SHA_xor(W)."""
    Hr = sha_hash(W); Hx = sha_xor(W)
    return [Hr[i]^Hx[i] for i in range(8)]


# ============================================================
# EXP 1: Per-bit Ψ predictability
# ============================================================
def exp1_per_bit(N):
    print("="*70)
    print("EXP 1: PER-BIT Ψ STABILITY")
    print("  For random W: how often is Ψ_total[reg][bit] = 1?")
    print("  If P ≠ 0.5 → that bit is BIASED → predictable!")
    print("="*70)

    bit_counts = [[0]*32 for _ in range(8)]

    for _ in range(N):
        W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        P = psi_total(W)
        for reg in range(8):
            for bit in range(32):
                if (P[reg] >> bit) & 1:
                    bit_counts[reg][bit] += 1

    # Find biased bits
    print(f"\n  Bits with |P-0.5| > 0.02 (N={N}):")
    biased = []
    for reg in range(8):
        for bit in range(32):
            p = bit_counts[reg][bit] / N
            if abs(p - 0.5) > 0.02:
                biased.append((abs(p-0.5), reg, bit, p))

    biased.sort(reverse=True)
    if biased:
        for bias, reg, bit, p in biased[:15]:
            print(f"    Ψ[{reg}][{bit:2d}]: P={p:.4f} (bias={bias:.4f})")
    else:
        print(f"    None found at threshold 0.02")

    # Average bias
    all_biases = [abs(bit_counts[r][b]/N - 0.5) for r in range(8) for b in range(32)]
    avg_bias = sum(all_biases)/len(all_biases)
    max_bias = max(all_biases)
    print(f"\n  Average |bias|: {avg_bias:.4f}")
    print(f"  Max |bias|: {max_bias:.4f}")
    print(f"  Expected random (N={N}): ~{0.5/N**0.5:.4f}")


# ============================================================
# EXP 2: δΨ predictability — how Ψ changes with perturbation
# ============================================================
def exp2_dpsi_stability(N):
    print("\n" + "="*70)
    print("EXP 2: δΨ STABILITY across W_base")
    print("  Fix δW, vary W_base: how stable is δΨ_total?")
    print("  High stability → δΨ depends mainly on δW, not W_base")
    print("="*70)

    for word in [0, 4, 8, 12]:
        # Fix δ: W[word] += 1
        dpsi_samples = []

        for _ in range(N):
            W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
            W2 = list(W1); W2[word] = add32(W1[word], 1)

            P1 = psi_total(W1); P2 = psi_total(W2)
            dpsi = [P1[i] ^ P2[i] for i in range(8)]
            dpsi_samples.append(dpsi)

        # Per-bit consistency: for each bit, how often is δΨ[bit]=1?
        bit_freq = [[0]*32 for _ in range(8)]
        for dpsi in dpsi_samples:
            for reg in range(8):
                for bit in range(32):
                    if (dpsi[reg]>>bit)&1:
                        bit_freq[reg][bit] += 1

        # Predictability = |freq/N - 0.5|. Higher = more predictable.
        predictabilities = []
        for reg in range(8):
            for bit in range(32):
                p = bit_freq[reg][bit] / N
                predictabilities.append(abs(p - 0.5))

        avg_pred = sum(predictabilities)/len(predictabilities)
        max_pred = max(predictabilities)
        n_biased = sum(1 for p in predictabilities if p > 0.05)

        print(f"\n  δW = W[{word}]+=1:")
        print(f"    Avg predictability: {avg_pred:.4f}")
        print(f"    Max predictability: {max_pred:.4f}")
        print(f"    Bits with |bias|>0.05: {n_biased}/256")
        print(f"    Expected random: ~{0.5/N**0.5:.4f}")

        if max_pred > 0.1:
            # Show most predictable bits
            bits_sorted = sorted([(abs(bit_freq[r][b]/N-0.5), r, b, bit_freq[r][b]/N)
                                  for r in range(8) for b in range(32)], reverse=True)
            print(f"    Top 5 predictable δΨ bits:")
            for pred, reg, bit, p in bits_sorted[:5]:
                print(f"      δΨ[{reg}][{bit:2d}]: P(=1)={p:.3f} (pred={pred:.3f})")


# ============================================================
# EXP 3: Can we AMPLIFY predictability through δW choice?
# ============================================================
def exp3_amplify(N):
    print("\n" + "="*70)
    print("EXP 3: AMPLIFY PREDICTABILITY")
    print("  Screen: which δW gives most predictable δΨ?")
    print("="*70)

    W_base = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    P_base = psi_total(W_base)

    best_pred = 0
    best_desc = ""

    # Test many δW types
    results = []
    for word in range(16):
        for delta_type in ['xor_bit0', 'xor_bit31', 'add_1', 'add_minus1']:
            W2 = list(W_base)
            if delta_type == 'xor_bit0': W2[word] ^= 1
            elif delta_type == 'xor_bit31': W2[word] ^= 0x80000000
            elif delta_type == 'add_1': W2[word] = add32(W_base[word], 1)
            elif delta_type == 'add_minus1': W2[word] = add32(W_base[word], MASK32)

            # Measure δΨ across many W_base samples
            dpsi_hw_list = []
            for _ in range(min(N, 50)):
                Wb = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
                Wb2 = list(Wb)
                if delta_type == 'xor_bit0': Wb2[word] ^= 1
                elif delta_type == 'xor_bit31': Wb2[word] ^= 0x80000000
                elif delta_type == 'add_1': Wb2[word] = add32(Wb[word], 1)
                elif delta_type == 'add_minus1': Wb2[word] = add32(Wb[word], MASK32)

                P1 = psi_total(Wb); P2 = psi_total(Wb2)
                dpsi_hw_list.append(sum(hw(P1[i]^P2[i]) for i in range(8)))

            avg_hw = sum(dpsi_hw_list)/len(dpsi_hw_list)
            deviation = 128 - avg_hw  # Positive = more predictable

            results.append((deviation, word, delta_type, avg_hw))

    results.sort(reverse=True)

    print(f"\n  Top 10 most predictable δW (highest deviation from 128):")
    for dev, word, dtype, avg in results[:10]:
        marker = " ★" if dev > 3 else ""
        print(f"    W[{word:2d}] {dtype:12s}: avg HW(δΨ)={avg:.1f}, deviation={dev:+.1f}{marker}")

    print(f"\n  Bottom 3 (least predictable):")
    for dev, word, dtype, avg in results[-3:]:
        print(f"    W[{word:2d}] {dtype:12s}: avg HW(δΨ)={avg:.1f}, deviation={dev:+.1f}")

    # Key metric: how many bits of Ψ are predictable?
    top = results[0]
    print(f"\n  BEST: W[{top[1]}] {top[2]}: {top[3]:.1f} bits change (vs 128 random)")
    print(f"  Predictable bits: {128 - top[3]:.1f}")
    print(f"  Exploitable? {128 - top[3]:.1f} > 2 → maybe!")


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 200

    exp1_per_bit(N)
    exp2_dpsi_stability(min(N, 100))
    exp3_amplify(min(N, 100))
