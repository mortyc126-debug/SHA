#!/usr/bin/env python3
"""
Noise source isolation.

From methodology T8: rotation is the ONLY essential noise source.
But T8 measured D2 (second derivative). We measure word-level structure.

Hypothesis: Σ₀/Σ₁ are the main POSITIONAL mixers.
Without them: structure survives longer → more rounds visible.

Experiment: SHA-256 variants with specific operations disabled.
Measure: at which round does δa mod 8 become uniform?

Variants:
V1: Full SHA-256 (baseline)
V2: No Σ₀, no Σ₁ (replace with identity: sig0(x)=x, sig1(x)=x)
V3: No Ch, no Maj (replace with XOR: ch=e^f^g, maj=a^b^c)
V4: No carry (replace + with XOR)
V5: No Σ₁ only (keep Σ₀)
V6: No Σ₀ only (keep Σ₁)
"""

import random, math
from collections import defaultdict

MASK32 = 0xFFFFFFFF
K=[0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2]
IV=[0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]
def rotr(x,n):return((x>>n)|(x<<(32-n)))&MASK32
def ssig0(x):return rotr(x,7)^rotr(x,18)^(x>>3)
def ssig1(x):return rotr(x,17)^rotr(x,19)^(x>>10)
def hw(x):return bin(x&MASK32).count('1')

def sha_variant(M, variant="full"):
    """SHA-256 with selectively disabled components."""
    W=list(M)+[0]*(64-len(M))
    for i in range(16,64):W[i]=(ssig1(W[i-2])+W[i-7]+ssig0(W[i-15])+W[i-16])&MASK32
    a,b,c,d,e,f,g,h=IV
    states=[(a,b,c,d,e,f,g,h)]

    for r in range(64):
        # Σ₁(e) — round rotation
        if variant in ("no_sigma", "no_sig1"):
            s1e = e  # identity
        else:
            s1e = rotr(e,6)^rotr(e,11)^rotr(e,25)

        # Ch(e,f,g)
        if variant == "no_ch_maj":
            chv = e ^ f ^ g  # linear substitute
        else:
            chv = (e&f)^(~e&g)&MASK32

        # T1 = h + Σ₁(e) + Ch(e,f,g) + K[r] + W[r]
        if variant == "no_carry":
            T1 = h ^ s1e ^ chv ^ K[r] ^ W[r]  # XOR instead of ADD
        else:
            T1 = (h + s1e + chv + K[r] + W[r]) & MASK32

        # Σ₀(a)
        if variant in ("no_sigma", "no_sig0"):
            s0a = a  # identity
        else:
            s0a = rotr(a,2)^rotr(a,13)^rotr(a,22)

        # Maj(a,b,c)
        if variant == "no_ch_maj":
            majv = a ^ b ^ c  # linear substitute
        else:
            majv = (a&b)^(a&c)^(b&c)

        # T2 = Σ₀(a) + Maj(a,b,c)
        if variant == "no_carry":
            T2 = s0a ^ majv
        else:
            T2 = (s0a + majv) & MASK32

        # e_new = d + T1
        if variant == "no_carry":
            e_new = d ^ T1
        else:
            e_new = (d + T1) & MASK32

        # a_new = T1 + T2
        if variant == "no_carry":
            a_new = T1 ^ T2
        else:
            a_new = (T1 + T2) & MASK32

        h,g,f=g,f,e; e=e_new; d,c,b=c,b,a; a=a_new
        states.append((a,b,c,d,e,f,g,h))

    return states

N = 3000
random.seed(42)

variants = [
    ("full",       "Full SHA-256"),
    ("no_sigma",   "No Σ₀ and Σ₁ (identity)"),
    ("no_sig1",    "No Σ₁ only"),
    ("no_sig0",    "No Σ₀ only"),
    ("no_ch_maj",  "No Ch/Maj (XOR substitute)"),
    ("no_carry",   "No carry (XOR instead of ADD)"),
]

# ================================================================
# For each variant: measure δa mod 8 uniformity per round
# Structure survival = last round where mod 8 is non-uniform
# ================================================================

print("=" * 70)
print("Noise source isolation: structure survival per variant")
print("=" * 70)
print(f"\n  N={N}, δM = W[0] bit 0 flip")

for var_name, var_label in variants:
    # Collect δa mod 8 at each round
    mod8_chi2 = []

    for r in range(65):
        counts = [0]*8
        random.seed(42)

        for _ in range(N):
            M1 = [random.randint(0,MASK32) for _ in range(16)]
            M2 = list(M1); M2[0] ^= 1

            ss1 = sha_variant(M1, var_name)
            ss2 = sha_variant(M2, var_name)

            da = (ss2[r][0] - ss1[r][0]) & MASK32
            counts[da % 8] += 1

        expected = N / 8
        chi2 = sum((counts[k] - expected)**2 / expected for k in range(8))
        chi2_norm = chi2 / 7  # normalized by DOF
        mod8_chi2.append(chi2_norm)

    # Find last round where structure exists (chi2_norm > 3)
    last_structured = 0
    for r in range(65):
        if mod8_chi2[r] > 3:
            last_structured = r

    # Also measure HW(δstate) convergence
    hw_at_8 = 0
    hw_at_16 = 0
    hw_at_64 = 0
    random.seed(42)
    for _ in range(min(N, 500)):
        M1 = [random.randint(0,MASK32) for _ in range(16)]
        M2 = list(M1); M2[0] ^= 1
        ss1 = sha_variant(M1, var_name)
        ss2 = sha_variant(M2, var_name)
        hw_at_8 += sum(hw(ss1[8][i]^ss2[8][i]) for i in range(8)) / 500
        hw_at_16 += sum(hw(ss1[16][i]^ss2[16][i]) for i in range(8)) / 500
        hw_at_64 += sum(hw(ss1[64][i]^ss2[64][i]) for i in range(8)) / 500

    print(f"\n  {var_label}:")
    print(f"    Structure survives until round: {last_structured}")
    print(f"    HW(δstate): r=8:{hw_at_8:.0f}  r=16:{hw_at_16:.0f}  r=64:{hw_at_64:.0f}")
    print(f"    mod8 χ² by round: r=1:{mod8_chi2[1]:.0f} r=2:{mod8_chi2[2]:.0f} r=3:{mod8_chi2[3]:.0f} r=4:{mod8_chi2[4]:.0f} r=5:{mod8_chi2[5]:.0f} r=8:{mod8_chi2[8]:.0f} r=16:{mod8_chi2[min(16,64)]:.1f}")


print()
print("=" * 70)
print("SUMMARY")
print("=" * 70)
print("""
  Structure survival (rounds until δa mod 8 = uniform):

  Full SHA-256:           → rounds until uniform (baseline)
  No Σ₀/Σ₁ (no rotation): → if MUCH longer = Σ is main noise source
  No Ch/Maj (linear):     → if MUCH longer = Ch/Maj is noise source
  No carry (XOR only):    → if MUCH longer = carry is noise source

  Also: HW(δstate) convergence speed shows TOTAL mixing rate.
  HW→128 fast = good mixing. HW stays low = poor mixing.
""")
