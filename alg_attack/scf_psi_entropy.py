#!/usr/bin/env python3
"""
SCF: Ψ-ENTROPY — точная информационная мера предсказуемости δΨ.

Вопрос: сколько бит РЕАЛЬНОЙ энтропии в δΨ_total?
Если H(δΨ) < 256 → предсказуемость > 0 → advantage.

advantage = 256 - H(δΨ) бит.
Collision cost = 2^{(128 + H(δΨ)/2 - 128)} = 2^{H(δΨ)/2}?
No: collision cost = 2^{H(δΨ)/2} (birthday на entropy space).

Если H(δΨ) = 250 → birthday = 2^{125} (3 бита advantage).
Если H(δΨ) = 240 → birthday = 2^{120} (8 бит advantage).

Измеряем H(δΨ) точно через per-bit entropy.
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
# EXP 1: Per-bit entropy of δΨ (for fixed δW, varying W_base)
# ============================================================
def exp1_dpsi_entropy(N, delta_word=0, delta_val=0x80000000):
    print("="*70)
    print(f"EXP 1: PER-BIT ENTROPY OF δΨ")
    print(f"  δW: W[{delta_word}] ^= 0x{delta_val:08x}")
    print(f"  N = {N} random W_base messages")
    print("="*70)

    bit_ones = [[0]*32 for _ in range(8)]

    for _ in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = list(W1); W2[delta_word] ^= delta_val

        Hr1 = sha_hash(W1); Hr2 = sha_hash(W2)
        Hx1 = sha_xor(W1); Hx2 = sha_xor(W2)

        P1 = [Hr1[i]^Hx1[i] for i in range(8)]
        P2 = [Hr2[i]^Hx2[i] for i in range(8)]
        dP = [P1[i]^P2[i] for i in range(8)]

        for reg in range(8):
            for bit in range(32):
                if (dP[reg]>>bit)&1:
                    bit_ones[reg][bit] += 1

    # Per-bit entropy
    total_entropy = 0
    low_entropy_bits = []

    for reg in range(8):
        for bit in range(32):
            p = bit_ones[reg][bit] / N
            h = binary_entropy(p)
            total_entropy += h
            if h < 0.99:
                low_entropy_bits.append((h, reg, bit, p))

    low_entropy_bits.sort()

    print(f"\n  Total per-bit entropy: {total_entropy:.2f} / 256 bits")
    print(f"  Entropy deficit: {256 - total_entropy:.2f} bits")
    print(f"  Effective security: 2^{total_entropy/2:.1f} (vs 2^128 ideal)")

    if low_entropy_bits:
        print(f"\n  Low-entropy bits (H < 0.99):")
        for h, reg, bit, p in low_entropy_bits[:15]:
            deficit = 1 - h
            print(f"    δΨ[{reg}][{bit:2d}]: P={p:.3f}, H={h:.4f}, deficit={deficit:.4f}")

    return total_entropy


# ============================================================
# EXP 2: Entropy for MULTIPLE δW choices
# ============================================================
def exp2_entropy_by_delta(N):
    print("\n" + "="*70)
    print("EXP 2: ENTROPY BY δW CHOICE")
    print("="*70)

    results = []

    deltas = [
        (0, 0x80000000, "W[0] MSB"),
        (0, 0x00000001, "W[0] LSB"),
        (0, 0x00000003, "W[0] 2bit"),
        (1, 0x80000000, "W[1] MSB"),
        (4, 0x80000000, "W[4] MSB"),
        (7, 0x80000000, "W[7] MSB"),
        (15, 0x80000000, "W[15] MSB"),
    ]

    for word, val, name in deltas:
        bit_ones = [[0]*32 for _ in range(8)]

        for _ in range(N):
            W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
            W2 = list(W1); W2[word] ^= val

            P1 = [sha_hash(W1)[i]^sha_xor(W1)[i] for i in range(8)]
            P2 = [sha_hash(W2)[i]^sha_xor(W2)[i] for i in range(8)]
            dP = [P1[i]^P2[i] for i in range(8)]

            for reg in range(8):
                for bit in range(32):
                    if (dP[reg]>>bit)&1:
                        bit_ones[reg][bit] += 1

        total_h = sum(binary_entropy(bit_ones[r][b]/N)
                      for r in range(8) for b in range(32))
        deficit = 256 - total_h
        results.append((deficit, name, total_h))

    results.sort(reverse=True)

    print(f"\n  {'δW':>15} | {'H(δΨ)':>8} | {'deficit':>8} | {'security':>10}")
    print("  " + "-"*50)
    for deficit, name, h in results:
        print(f"  {name:>15} | {h:8.2f} | {deficit:8.2f} | 2^{h/2:.1f}")

    best = results[0]
    print(f"\n  Best: {best[1]} with deficit={best[0]:.2f} bits")
    print(f"  → Security = 2^{best[2]/2:.1f} (vs 2^128 ideal)")
    print(f"  → Advantage: {128 - best[2]/2:.1f} bits below birthday")


# ============================================================
# EXP 3: Pair correlations in δΨ
# ============================================================
def exp3_pair_correlations(N):
    print("\n" + "="*70)
    print("EXP 3: PAIR CORRELATIONS IN δΨ")
    print("  Are δΨ bits correlated? If yes → entropy lower than per-bit estimate")
    print("="*70)

    delta_word, delta_val = 0, 0x80000000

    # Sample δΨ vectors
    samples = []
    for _ in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = list(W1); W2[delta_word] ^= delta_val
        P1 = [sha_hash(W1)[i]^sha_xor(W1)[i] for i in range(8)]
        P2 = [sha_hash(W2)[i]^sha_xor(W2)[i] for i in range(8)]
        dP = tuple(P1[i]^P2[i] for i in range(8))
        samples.append(dP)

    # Check pair correlations between bits of different registers
    n_correlated = 0
    n_tested = 0

    for r1 in range(8):
        for r2 in range(r1+1, 8):
            for b in [0, 15, 31]:
                agree = sum(1 for s in samples
                           if ((s[r1]>>b)&1) == ((s[r2]>>b)&1))
                p_agree = agree / N
                n_tested += 1
                if abs(p_agree - 0.5) > 0.05:
                    n_correlated += 1

    print(f"  Cross-register correlations (same bit position):")
    print(f"    Tested: {n_tested}, Correlated (|P-0.5|>0.05): {n_correlated}")

    # HW distribution — is it normal?
    hw_list = [sum(hw(s[i]) for i in range(8)) for s in samples]
    avg_hw = sum(hw_list)/N
    std_hw = (sum((h-avg_hw)**2 for h in hw_list)/N)**0.5

    print(f"\n  HW(δΨ) distribution:")
    print(f"    avg={avg_hw:.1f}, std={std_hw:.1f}")
    print(f"    Expected: avg=128, std=8 (for 256 independent bits)")
    print(f"    Deviation from normal: avg={abs(avg_hw-128):.1f}, std_diff={abs(std_hw-8):.1f}")

    if abs(std_hw - 8) > 1:
        effective_indep = (std_hw/0.5)**2  # σ² = np(1-p), n=effective independent bits
        print(f"    Effective independent bits: {effective_indep:.0f} (vs 256)")
        print(f"    → Entropy ≈ {effective_indep:.0f} bits")
        print(f"    → Security ≈ 2^{effective_indep/2:.0f}")


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 500

    h1 = exp1_dpsi_entropy(N, delta_word=0, delta_val=0x80000000)
    exp2_entropy_by_delta(min(N, 300))
    exp3_pair_correlations(min(N, 300))

    print(f"\n{'='*70}")
    print(f"ИТОГ: Ψ-ЭНТРОПИЯ")
    print(f"{'='*70}")
    print(f"  H(δΨ_total) ≈ {h1:.1f} бит (из 256 возможных)")
    print(f"  Deficit = {256-h1:.1f} бит")
    print(f"  SHA-256 collision security = 2^{h1/2:.1f}")
    print(f"  Theoretical maximum = 2^128")
    print(f"  OUR RESULT: 2^{h1/2:.1f} ← точный пол безопасности")
