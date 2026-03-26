#!/usr/bin/env python3
"""
SCF: INFORMATION-THEORETIC CHANNEL MODEL — SHA как зашумлённый канал.

SHA(W) = SHA_xor(W) ⊕ Ψ(W)
       = signal     ⊕ noise

Signal: SHA_xor(W) — ВЫЧИСЛИМ (линейная функция W)
Noise:  Ψ(W) — random 256-bit (доказано: H=255.6)

Collision = два сообщения с одинаковым выходом.
В channel coding: два codewords дают одинаковый received signal
после прохождения через noisy channel.

Channel: input = δW (512 bit), output = δH (256 bit)
Noise:   δΨ (256 bit, ~random)

δH = δSHA_xor ⊕ δΨ

Collision: δH = 0 → δSHA_xor = δΨ.

Channel capacity: C = max_{p(δW)} I(δW; δH)
Если C > 0: можно кодировать информацию через SHA → collision дешевле.
Если C = 0: SHA = чистый шум → birthday оптимален.

Измерим I(δW; δH) — mutual information.
"""
import os, sys, math
from collections import Counter

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

def binary_entropy(p):
    if p <= 0 or p >= 1: return 0
    return -p*math.log2(p) - (1-p)*math.log2(1-p)


# ============================================================
# EXP 1: Channel model — SHA as BSC (Binary Symmetric Channel)
# ============================================================
def exp1_channel_model(N):
    print("="*70)
    print("EXP 1: SHA AS BINARY SYMMETRIC CHANNEL")
    print("  Input:  δW (message difference)")
    print("  Signal: δSHA_xor (computable)")
    print("  Noise:  δΨ (carry correction)")
    print("  Output: δH = δSHA_xor ⊕ δΨ")
    print("="*70)

    # For a BSC with crossover probability p:
    # Capacity C = 1 - H(p) bits per channel use

    # Measure: for each bit position, what's P(δH[bit] ≠ δSHA_xor[bit])?
    # This = P(δΨ[bit] = 1) = crossover probability

    crossover_per_bit = [[0]*32 for _ in range(8)]
    total = 0

    for _ in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = list(W1); W2[0] ^= 0x80000000

        H1 = sha_hash(W1); H2 = sha_hash(W2)
        X1 = sha_xor(W1); X2 = sha_xor(W2)

        dH = [H1[i]^H2[i] for i in range(8)]
        dX = [X1[i]^X2[i] for i in range(8)]
        dPsi = [dH[i]^dX[i] for i in range(8)]  # noise

        total += 1
        for reg in range(8):
            for bit in range(32):
                if (dPsi[reg]>>bit)&1:
                    crossover_per_bit[reg][bit] += 1

    # Per-bit crossover probability
    p_cross = [[crossover_per_bit[r][b]/total for b in range(32)] for r in range(8)]

    # Per-bit capacity
    total_capacity = 0
    bit_capacities = []
    for reg in range(8):
        for bit in range(32):
            p = p_cross[reg][bit]
            c = 1 - binary_entropy(p)
            total_capacity += c
            bit_capacities.append((c, reg, bit, p))

    bit_capacities.sort(reverse=True)

    print(f"\n  Per-bit crossover probability (P(δΨ[bit]=1)):")
    print(f"  Expected: 0.5 for pure noise")

    # Summary
    avg_p = sum(p_cross[r][b] for r in range(8) for b in range(32))/256
    print(f"\n  Average crossover: P = {avg_p:.4f} (ideal noise: 0.500)")
    print(f"  Total channel capacity: {total_capacity:.4f} bits / 256 uses")
    print(f"  Per-bit average capacity: {total_capacity/256:.6f}")

    if total_capacity > 1:
        print(f"\n  ★ Nonzero capacity! {total_capacity:.2f} bits of information leak through!")
    else:
        print(f"\n  Near-zero capacity — SHA is essentially a random oracle")

    # Top capacity bits
    print(f"\n  Top 10 highest-capacity bits:")
    for c, reg, bit, p in bit_capacities[:10]:
        print(f"    H[{reg}][{bit:2d}]: crossover={p:.4f}, capacity={c:.6f}")

    return total_capacity


# ============================================================
# EXP 2: Mutual information I(δW; δH) per bit of input
# ============================================================
def exp2_mutual_info(N):
    print("\n" + "="*70)
    print("EXP 2: MUTUAL INFORMATION I(δW_bit; δH)")
    print("  How much does ONE input bit tell about the hash?")
    print("="*70)

    # For each input bit: measure I(δW[word][bit]; δH)
    # I = H(δH) - H(δH | δW[bit])
    # Since δH has 256 bits, measure per output bit

    W_base = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    H_base = sha_hash(W_base)

    # Collect δH for each input bit flip
    results = []

    for word in range(16):
        for in_bit in [0, 15, 31]:  # Sample 3 bits per word
            # Collect δH for this input bit flip across N bases
            dH_samples = []
            for _ in range(N):
                W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
                W2 = list(W1); W2[word] ^= (1 << in_bit)
                H1 = sha_hash(W1); H2 = sha_hash(W2)
                dH = tuple(H1[i]^H2[i] for i in range(8))
                dH_samples.append(dH)

            # Per output bit: measure P(δH[out_bit]=1) for this input flip
            # If P ≠ 0.5 → mutual information > 0
            total_mi = 0
            for out_reg in range(8):
                for out_bit in range(32):
                    p1 = sum(1 for d in dH_samples if (d[out_reg]>>out_bit)&1)/N
                    h = binary_entropy(p1)
                    mi = 1.0 - h  # Mutual info for this bit pair
                    total_mi += mi

            results.append((total_mi, word, in_bit))

    results.sort(reverse=True)
    print(f"\n  Total MI per input bit (sum over 256 output bits):")
    print(f"  {'input':>12} | MI (bits)")
    print("  " + "-"*30)
    for mi, word, bit in results[:10]:
        print(f"  W[{word:2d}][{bit:2d}]  | {mi:.4f}")
    for mi, word, bit in results[-3:]:
        print(f"  W[{word:2d}][{bit:2d}]  | {mi:.4f}")

    avg_mi = sum(r[0] for r in results)/len(results)
    print(f"\n  Average MI per input bit: {avg_mi:.4f}")
    print(f"  Total MI (all 512 input bits): ≈{avg_mi*512:.1f}")
    print(f"  Expected for random function: ~0")

    return avg_mi


# ============================================================
# EXP 3: Shannon limit — collision cost from channel capacity
# ============================================================
def exp3_shannon_limit(capacity):
    print("\n" + "="*70)
    print("EXP 3: SHANNON LIMIT FOR COLLISION")
    print("="*70)

    print(f"""
  Channel model:
    Input: δW ∈ {{0,1}}^512
    Output: δH ∈ {{0,1}}^256
    Noise: δΨ ∈ {{0,1}}^256 (carry correction)

    Channel capacity per bit: C ≈ {capacity/256:.6f}
    Total capacity: {capacity:.4f} bits per 256-bit block

  Collision = output δH = 0.
  This is a ZERO-ERROR coding problem:
    Find two inputs (W1, W2) that produce the SAME output.

  For a BSC with capacity C:
    Zero-error capacity C₀ ≤ C.
    For crossover p ≈ 0.5: C ≈ 0, C₀ = 0.

  Shannon limit for collision search:
    If C₀ = 0: no encoding can guarantee collision.
    Collision requires EXHAUSTIVE search → birthday O(2^{{128}}).

  If C₀ > 0 (even tiny):
    Can encode δW to "steer" through noise.
    Collision cost = O(2^{{(256-C_total)/2}}).
    With C_total = {capacity:.2f}: cost = O(2^{{{(256-capacity)/2:.1f}}}).

  CONCLUSION:
    Measured C = {capacity:.2f} bits total.
    This gives collision cost = 2^{{{(256-capacity)/2:.1f}}} (vs 2^128 birthday).
    Advantage: {128-(256-capacity)/2:.2f} bits.
    """)

    if capacity < 1:
        print(f"  ★ Capacity < 1 bit → advantage < 0.5 bits → NEGLIGIBLE")
        print(f"  SHA-256 is information-theoretically at birthday bound.")
    else:
        print(f"  ★ Capacity = {capacity:.1f} bits → advantage = {128-(256-capacity)/2:.1f} bits")


# ============================================================
# EXP 4: Is noise INPUT-DEPENDENT? (key question)
# ============================================================
def exp4_noise_independence(N):
    print("\n" + "="*70)
    print("EXP 4: IS NOISE INPUT-DEPENDENT?")
    print("  If δΨ depends on δW → channel is NOT memoryless BSC")
    print("  → potentially exploitable via input design")
    print("="*70)

    # For two DIFFERENT δW: do they produce CORRELATED δΨ?
    W_base = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

    correlations = []

    for _ in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

        # δW type A: flip W[0] MSB
        W2a = list(W1); W2a[0] ^= 0x80000000
        Ha1 = sha_hash(W1); Ha2 = sha_hash(W2a)
        Xa1 = sha_xor(W1); Xa2 = sha_xor(W2a)
        dPsi_a = tuple((Ha1[i]^Ha2[i])^(Xa1[i]^Xa2[i]) for i in range(8))

        # δW type B: flip W[0] LSB
        W2b = list(W1); W2b[0] ^= 1
        Hb1 = sha_hash(W1); Hb2 = sha_hash(W2b)
        Xb1 = sha_xor(W1); Xb2 = sha_xor(W2b)
        dPsi_b = tuple((Hb1[i]^Hb2[i])^(Xb1[i]^Xb2[i]) for i in range(8))

        # Correlation: HW(δΨ_a ⊕ δΨ_b)
        # If independent: HW ≈ 128. If correlated: HW < 128.
        xor_hw = sum(hw(dPsi_a[i]^dPsi_b[i]) for i in range(8))
        correlations.append(xor_hw)

    avg_corr = sum(correlations)/len(correlations)
    print(f"  HW(δΨ_MSB ⊕ δΨ_LSB): avg={avg_corr:.1f} (128=independent)")
    print(f"  Deviation: {128-avg_corr:+.1f}")

    if abs(128 - avg_corr) > 2:
        print(f"  ★ Noise IS input-dependent! δΨ for different δW are correlated!")
    else:
        print(f"  Noise is input-INDEPENDENT. SHA = memoryless BSC. Birthday optimal.")


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 500

    C = exp1_channel_model(N)
    mi = exp2_mutual_info(min(N, 200))
    exp3_shannon_limit(C)
    exp4_noise_independence(min(N, 300))
