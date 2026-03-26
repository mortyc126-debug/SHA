#!/usr/bin/env python3
"""
SCF: MULTIBLOCK + COMPENSATION — два блока, две степени свободы.

Идея 1 (Multiblock):
  Block 1: (M1, M1') → H1, H1' с ΔH1 = H1⊕H1'
  Block 2: compress(H1, M2) vs compress(H1', M2)
    Разные IV! ΔIV = ΔH1 ≈ 100 бит.
    Нужно: compress(H1, M2) = compress(H1', M2)
    → collision с ΔIV≠0 (semi-free-start)

Идея 2 (Compensation):
  Найти M2, M2' такие что ΔH2 = −ΔH1 (мод Z).
  Тогда H_final1 = H1 + state1 = H1' + state2 → collision!
  Нужно: ΔH2 точно компенсирует ΔH1.

Идея 3 (Birthday на near-collisions):
  Сгенерировать МНОГО пар с HW(ΔH)≈100.
  Найти две пары где ΔH1 = ΔH2.
  Тогда M1||M2 и M1'||M2' — collision!
  Стоимость: birthday на 2^100 space → O(2^50)?
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

def sha_compress(W16, iv=None):
    if iv is None: iv = IV
    W=expand_real(W16); s=list(iv)
    for r in range(64):
        a,b,c,d,e,f,g,h=s
        T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add32(Sig0(a),Maj(a,b,c))
        s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
    return [add32(iv[i],s[i]) for i in range(8)]

def hash_diff_iv(W1, W2, iv):
    H1 = sha_compress(W1, iv); H2 = sha_compress(W2, iv)
    return sum(hw(H1[i]^H2[i]) for i in range(8))


# ============================================================
# EXP 1: Multiblock — how much does custom IV help?
# ============================================================
def exp1_custom_iv():
    print("="*70)
    print("EXP 1: CUSTOM IV — does different IV improve near-collision?")
    print("="*70)

    results_standard = []
    results_custom = []

    for trial in range(20):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = list(W1); W2[0] ^= 0x80000000

        # Standard IV
        h_std = hash_diff_iv(W1, W2, IV)
        results_standard.append(h_std)

        # Random custom IV (as if from block 1)
        custom_iv = [int.from_bytes(os.urandom(4),'big') for _ in range(8)]
        h_custom = hash_diff_iv(W1, W2, custom_iv)
        results_custom.append(h_custom)

    avg_std = sum(results_standard)/len(results_standard)
    avg_cust = sum(results_custom)/len(results_custom)
    print(f"  Standard IV: avg={avg_std:.1f}, min={min(results_standard)}")
    print(f"  Custom IV:   avg={avg_cust:.1f}, min={min(results_custom)}")


# ============================================================
# EXP 2: Near-collision PROFILE — where are the nonzero bits?
# ============================================================
def exp2_profile(N):
    print("\n" + "="*70)
    print("EXP 2: NEAR-COLLISION ΔH PROFILE")
    print("  Which registers carry most of the HW≈100?")
    print("="*70)

    per_reg = [0.0]*8
    per_bit = [[0]*32 for _ in range(8)]

    for _ in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = list(W1); W2[0] ^= 0x80000000

        # SA to get near-collision
        H1 = sha_compress(W1)
        best = 256; cur = list(W2)
        for it in range(500):
            t = list(cur)
            w = int.from_bytes(os.urandom(1),'big') % 16
            b = int.from_bytes(os.urandom(1),'big') % 32
            t[w] ^= (1<<b)
            if t == W1: continue
            H2 = sha_compress(t)
            s = sum(hw(H1[i]^H2[i]) for i in range(8))
            T = max(0.01, 1-it/500)
            if s < best or math.exp(-(s-best)/(T*2)) > int.from_bytes(os.urandom(4),'big')/(1<<32):
                cur = t
                if s < best: best = s; best_H2 = H2

        # Profile the ΔH
        for reg in range(8):
            dh = H1[reg] ^ best_H2[reg]
            per_reg[reg] += hw(dh)
            for bit in range(32):
                per_bit[reg][bit] += (dh >> bit) & 1

    regs = "abcdefgh"
    print(f"\n  Per-register avg HW of ΔH:")
    for reg in range(8):
        avg = per_reg[reg] / N
        bar = "█" * int(avg * 2)
        print(f"    H[{regs[reg]}]: {avg:.1f} {bar}")

    # Any biased bits?
    biased = []
    for reg in range(8):
        for bit in range(32):
            p = per_bit[reg][bit] / N
            if abs(p - 0.5) > 0.08:
                biased.append((abs(p-0.5), reg, bit, p))
    biased.sort(reverse=True)
    if biased:
        print(f"\n  Biased bits (|P-0.5| > 0.08):")
        for bias, reg, bit, p in biased[:10]:
            print(f"    H[{regs[reg]}][{bit}]: P={p:.3f} (bias={bias:.3f})")
    else:
        print(f"\n  No biased bits found (all ≈ 0.5)")


# ============================================================
# EXP 3: Birthday on near-collisions — can two ΔH match?
# ============================================================
def exp3_birthday_nc(N):
    print("\n" + "="*70)
    print("EXP 3: BIRTHDAY ON NEAR-COLLISIONS")
    print("  Generate many ΔH, find two that MATCH on most bits")
    print("  If ΔH1 = ΔH2 → two-block collision!")
    print("="*70)

    # Generate near-collisions
    deltas = []
    for trial in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = list(W1)
        # Random single-bit delta
        w = int.from_bytes(os.urandom(1),'big') % 16
        b = int.from_bytes(os.urandom(1),'big') % 32
        W2[w] ^= (1 << b)

        H1 = sha_compress(W1); H2 = sha_compress(W2)
        dH = tuple(H1[i] ^ H2[i] for i in range(8))
        deltas.append((dH, W1, W2))

    # Find closest pair (minimum XOR distance)
    best_dist = 256
    best_pair = None

    # Check all pairs (O(N^2))
    for i in range(len(deltas)):
        for j in range(i+1, len(deltas)):
            dH_i = deltas[i][0]
            dH_j = deltas[j][0]
            dist = sum(hw(dH_i[r] ^ dH_j[r]) for r in range(8))
            if dist < best_dist:
                best_dist = dist
                best_pair = (i, j)

    print(f"  Generated {N} near-collisions")
    print(f"  Closest pair: HW(ΔH1 ⊕ ΔH2) = {best_dist}")
    print(f"  (Need 0 for 2-block collision)")
    print(f"  Expected for N={N} random 256-bit vectors: ~{max(0, 256-int(2*math.log2(max(N,2))))}")

    if best_pair:
        i, j = best_pair
        dH_i = deltas[i][0]; dH_j = deltas[j][0]
        print(f"\n  Per-register distance:")
        for r in range(8):
            d = hw(dH_i[r] ^ dH_j[r])
            print(f"    H[{r}]: {d} bits differ")

    # Truncated birthday: match on first k bits
    for k_bits in [32, 64, 96, 128]:
        mask = (1 << min(k_bits, 32)) - 1
        seen = {}
        found = 0
        for idx, (dH, _, _) in enumerate(deltas):
            # Truncate to k_bits
            trunc = tuple(dH[r] & mask for r in range(min(k_bits//32 + 1, 8)))
            if trunc in seen:
                found += 1
                if found == 1:
                    print(f"\n  Truncated match at {k_bits} bits! (pair {seen[trunc]} and {idx})")
                    other = seen[trunc]
                    full_dist = sum(hw(deltas[idx][0][r] ^ deltas[other][0][r]) for r in range(8))
                    print(f"    Full distance: {full_dist} bits")
            else:
                seen[trunc] = idx

        if found == 0:
            print(f"  No {k_bits}-bit truncated match in {N} pairs")


# ============================================================
# EXP 4: TWO-BLOCK COLLISION ATTEMPT
# ============================================================
def exp4_two_block(N):
    print("\n" + "="*70)
    print("EXP 4: TWO-BLOCK COLLISION ATTEMPT")
    print("  Block 1: M1, M1' → H1, H1'")
    print("  Block 2: same M2 → compress(H1,M2) vs compress(H1',M2)")
    print("  Need: ΔH2_output = 0 starting from ΔIV = ΔH1")
    print("="*70)

    best_two_block = 256

    for trial in range(N):
        # Block 1: create a near-collision pair
        M1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        M1p = list(M1)
        M1p[int.from_bytes(os.urandom(1),'big')%16] ^= (1 << (int.from_bytes(os.urandom(1),'big')%32))

        H1 = sha_compress(M1)
        H1p = sha_compress(M1p)

        dH1 = sum(hw(H1[i]^H1p[i]) for i in range(8))

        # Block 2: try to make compress(H1, M2) = compress(H1p, M2)
        # This means: same M2, different IVs, want same output
        M2 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

        out_a = sha_compress(M2, iv=H1)
        out_b = sha_compress(M2, iv=H1p)
        raw_dist = sum(hw(out_a[i]^out_b[i]) for i in range(8))

        # SA on M2 to minimize dist
        best_m2 = raw_dist; cur_M2 = list(M2)
        for it in range(500):
            t = list(cur_M2)
            w = int.from_bytes(os.urandom(1),'big') % 16
            b = int.from_bytes(os.urandom(1),'big') % 32
            t[w] ^= (1<<b)
            oa = sha_compress(t, iv=H1)
            ob = sha_compress(t, iv=H1p)
            s = sum(hw(oa[i]^ob[i]) for i in range(8))
            T = max(0.01, 1-it/500)
            if s < best_m2 or math.exp(-(s-best_m2)/(T*2)) > int.from_bytes(os.urandom(4),'big')/(1<<32):
                cur_M2 = t
                if s < best_m2: best_m2 = s

        if best_m2 < best_two_block:
            best_two_block = best_m2
            if trial < 10 or best_m2 < 100:
                print(f"  Trial {trial}: ΔH1={dH1:3d}, Block2 raw={raw_dist:3d}, optimized={best_m2:3d}")

    print(f"\n  Best two-block: HW = {best_two_block}")
    print(f"  Single-block baseline: HW ≈ 100")
    print(f"  Two-block advantage: {100 - best_two_block:+d}")


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 50

    exp1_custom_iv()
    exp2_profile(min(N, 30))
    exp3_birthday_nc(min(N*2, 200))
    exp4_two_block(min(N, 50))
