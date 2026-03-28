#!/usr/bin/env python3
"""
Axis 5: Beyond Degree — Zero-Sum, Bit-Slice, Rotational
==========================================================

Three genuinely new angles:

  Z1: ZERO-SUM — for subspace V of dim k, is Σ_{x∈V} f(x) = 0?
      Degree n-1 at n≤28 gives (n-1)-th order zero-sum automatically.
      But what about LOWER order? k-th order zero-sum for k << n?
      If k=20 works → distinguisher at 2^20 for full 32-bit SHA-256.

  Z2: BIT-SLICE — treat each output bit independently.
      SHA-256 bit j of word w = function of 512 input bits.
      Carries couple adjacent bits. What if we model without carries?
      The carry-free approximation is LINEAR over GF(2).
      How good is this approximation for specific bits?

  Z3: NEAR-ROTATIONAL — Rot^k(SHA(M)) vs SHA(Rot^k(M)).
      T_ROTATIONAL_NEGATIVE says exact equality impossible.
      But how CLOSE do they get? HW(XOR) = ? Average vs minimum.
"""

import numpy as np
import time, sys

MASK = 0xFFFFFFFF
H0 = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
      0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]
K = [0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
     0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
     0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
     0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
     0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
     0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
     0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
     0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2]

def rotr32(x,n): return ((x>>n)|(x<<(32-n)))&MASK
def sig0(x): return rotr32(x,7)^rotr32(x,18)^(x>>3)
def sig1(x): return rotr32(x,17)^rotr32(x,19)^(x>>10)
def Sig0(x): return rotr32(x,2)^rotr32(x,13)^rotr32(x,22)
def Sig1(x): return rotr32(x,6)^rotr32(x,11)^rotr32(x,25)
def Ch(e,f,g): return (e&f)^(~e&g)&MASK
def Maj(a,b,c): return (a&b)^(a&c)^(b&c)
def hw(x): return bin(x&MASK).count('1')

def msg_sched(W16):
    W=list(W16)+[0]*48
    for i in range(16,64): W[i]=(sig1(W[i-2])+W[i-7]+sig0(W[i-15])+W[i-16])&MASK
    return W

def sha256(W16):
    W=msg_sched(W16); a,b,c,d,e,f,g,h=H0
    for r in range(64):
        T1=(h+Sig1(e)+Ch(e,f,g)+K[r]+W[r])&MASK; T2=(Sig0(a)+Maj(a,b,c))&MASK
        h,g,f=g,f,e; e=(d+T1)&MASK; d,c,b=c,b,a; a=(T1+T2)&MASK
    return [(v+iv)&MASK for v,iv in zip([a,b,c,d,e,f,g,h],H0)]


# ============================================================
# Z1: ZERO-SUM at low order
# ============================================================

def experiment_Z1(seed=15000):
    print("="*70)
    print("Z1: ZERO-SUM — Low-order zero-sum properties")
    print("="*70)

    # k-th order zero-sum: pick k directions e_{i1},...,e_{ik},
    # compute XOR over 2^k cube corners:
    # S = XOR_{S⊆{i1,...,ik}} f(x ⊕ (S subset bits))
    # If deg(f) < k → S = 0 always (zero-sum).
    # If deg(f) = k → S = top coefficient (may be 0 or 1).
    # If deg(f) > k → S depends on x.

    # We know degree = 32 at n=32 (full). So k=32 gives S=1 (non-zero).
    # But: does S=0 happen for SPECIFIC (x, directions) at k < 32?

    rng = np.random.RandomState(seed)

    # Test: random k-dimensional cubes in W[0] space
    print(f"\n  --- k-th order differential sums for H[7][b23] at 64 rounds ---")
    print(f"  For each k: pick random x and k random bit directions.")
    print(f"  Compute S = XOR over 2^k cube corners of H[7][b23].")
    print(f"  P(S=0) for random f of degree d: 0.5 if k<d, ~0 if k=d.")
    print()

    N_trials = 200
    w, b = 7, 23

    print(f"  {'k':>3} | {'P(S=0)':>8} | {'expected':>9} | {'Z':>6} | {'note':>20}")
    print(f"  {'-'*3}-+-{'-'*8}-+-{'-'*9}-+-{'-'*6}-+-{'-'*20}")

    for k in range(1, 22):
        n_zero = 0
        n_tests = min(N_trials, 500 if k <= 10 else 100)

        for _ in range(n_tests):
            x = int(rng.randint(0, 1 << 32))
            # Pick k random distinct bit positions
            bits = rng.choice(32, size=k, replace=False)

            # Compute XOR over 2^k cube corners
            xor_sum = 0
            for mask_val in range(1 << k):
                point = x
                for j in range(k):
                    if (mask_val >> j) & 1:
                        point ^= (1 << bits[j])
                H = sha256([point & MASK] + [0]*15)
                xor_sum ^= (H[w] >> b) & 1

            if xor_sum == 0:
                n_zero += 1

        p_zero = n_zero / n_tests
        # For degree-d function: if k<d, P(S=0)=0.5; if k=d, P(S=0)≈0; if k>d, P=0
        expected = 0.5  # since degree=32 >> k for all tested k
        z = (p_zero - 0.5) / (0.5 / np.sqrt(n_tests)) if n_tests > 0 else 0

        note = ""
        if p_zero > 0.6:
            note = "★ BIASED toward 0"
        elif p_zero < 0.4:
            note = "★ BIASED toward 1"

        print(f"  {k:3d} | {p_zero:8.3f} | {expected:9.3f} | {z:+6.2f} | {note}")

        if k >= 15 and abs(z) < 1:
            # Early stop if clearly random
            pass

    return


# ============================================================
# Z2: BIT-SLICE — carry-free linear approximation
# ============================================================

def experiment_Z2(N=5000, seed=15001):
    print(f"\n{'='*70}")
    print("Z2: BIT-SLICE — Carry-free (GF2) approximation accuracy")
    print(f"N={N}")
    print("="*70)

    rng = np.random.RandomState(seed)

    # SHA-256 without carries: replace all + with XOR.
    # This is the GF(2) linearization of SHA-256.
    # How close is SHA_GF2(M) to SHA_real(M)?

    def sha256_xor(W16):
        """SHA-256 with XOR instead of addition (carry-free)."""
        W = list(W16)+[0]*48
        for i in range(16,64):
            W[i] = sig1(W[i-2]) ^ W[i-7] ^ sig0(W[i-15]) ^ W[i-16]

        a,b,c,d,e,f,g,h = H0
        for r in range(64):
            T1 = h ^ Sig1(e) ^ Ch(e,f,g) ^ K[r] ^ W[r]
            T2 = Sig0(a) ^ Maj(a,b,c)
            h,g,f = g,f,e; e = d^T1; d,c,b = c,b,a; a = T1^T2
        return [(v^iv)&MASK for v,iv in zip([a,b,c,d,e,f,g,h],H0)]

    # Compare SHA_real vs SHA_XOR
    hw_diffs = np.zeros(8)
    bit_agree = np.zeros((8, 32))  # per-bit agreement rate

    t0 = time.time()
    for i in range(N):
        W16 = [int(rng.randint(0,1<<32)) for _ in range(16)]
        H_real = sha256(W16)
        H_xor = sha256_xor(W16)

        for w in range(8):
            diff = H_real[w] ^ H_xor[w]
            hw_diffs[w] += hw(diff)
            for b in range(32):
                if not ((diff >> b) & 1):
                    bit_agree[w, b] += 1

    hw_diffs /= N
    bit_agree /= N
    print(f"  Collected: {time.time()-t0:.1f}s")

    print(f"\n  --- HW(SHA_real XOR SHA_xor) per word ---")
    for w in range(8):
        print(f"    H[{w}]: E[HW(diff)] = {hw_diffs[w]:.2f} / 32 (random=16)")

    # Per-bit: which bits agree most?
    print(f"\n  --- Best-agreement bits (real ≈ XOR most often) ---")
    all_bits = []
    for w in range(8):
        for b in range(32):
            all_bits.append((w, b, bit_agree[w, b]))

    all_bits.sort(key=lambda x: -x[2])
    print(f"  Top-10 most predictable bits (carry-free model):")
    for w, b, ag in all_bits[:10]:
        marker = " ★" if ag > 0.55 else ""
        print(f"    H[{w}][b{b:2d}]: agreement = {ag:.4f}{marker}")

    print(f"\n  Bottom-10 (worst agreement):")
    for w, b, ag in all_bits[-10:]:
        print(f"    H[{w}][b{b:2d}]: agreement = {ag:.4f}")

    # Overall: how good is the GF2 approximation?
    mean_agree = np.mean(bit_agree)
    print(f"\n  Mean per-bit agreement: {mean_agree:.4f} (random=0.500)")
    print(f"  This means: carry-free model predicts {mean_agree*100:.1f}% of bits correctly.")
    if mean_agree > 0.52:
        print(f"  ★ GF2 approximation has measurable signal!")
    else:
        print(f"  GF2 approximation ≈ random. Carries destroy all structure.")

    return hw_diffs, bit_agree


# ============================================================
# Z3: NEAR-ROTATIONAL
# ============================================================

def experiment_Z3(N=5000, seed=15002):
    print(f"\n{'='*70}")
    print("Z3: NEAR-ROTATIONAL — HW(Rot(SHA(M)) XOR SHA(Rot(M)))")
    print(f"N={N}")
    print("="*70)

    rng = np.random.RandomState(seed)

    def rot_words(W16, k):
        """Rotate each 32-bit word by k positions."""
        return [rotr32(w, k) for w in W16]

    # For each rotation k=1..31: compute HW(Rot^k(SHA(M)) XOR SHA(Rot^k(M)))
    print(f"\n  {'k':>3} | {'E[HW_total]':>12} | {'min_HW':>8} | {'P(HW<200)':>10} | {'note':>15}")
    print(f"  {'-'*3}-+-{'-'*12}-+-{'-'*8}-+-{'-'*10}-+-{'-'*15}")

    best_k = -1; best_min = 256

    for k in [1, 2, 3, 4, 7, 8, 11, 13, 15, 16, 25]:
        hw_list = []
        for _ in range(N):
            W16 = [int(rng.randint(0,1<<32)) for _ in range(16)]
            H = sha256(W16)
            H_rot = [rotr32(h, k) for h in H]

            W16_rot = rot_words(W16, k)
            H_of_rot = sha256(W16_rot)

            total_hw = sum(hw(H_rot[w] ^ H_of_rot[w]) for w in range(8))
            hw_list.append(total_hw)

        mean_hw = np.mean(hw_list)
        min_hw = min(hw_list)
        p_low = np.mean(np.array(hw_list) < 200)

        note = ""
        if min_hw < best_min:
            best_min = min_hw; best_k = k
            note = "★ best"
        if mean_hw < 127:
            note += " ANOMALY"

        print(f"  {k:3d} | {mean_hw:12.2f} | {min_hw:8d} | {p_low:10.4f} | {note}")

    print(f"\n  Best k={best_k}: min HW = {best_min} (random=128)")
    if best_min < 100:
        print(f"  ★★★ NEAR-ROTATIONAL ANOMALY: min HW significantly below 128!")
    else:
        print(f"  No rotational anomaly. min HW ≈ random minimum.")

    return


# ============================================================
# MAIN
# ============================================================

if __name__ == '__main__':
    print("Axis 5: Beyond Degree — Zero-Sum, Bit-Slice, Rotational")
    print(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")

    N = 3000
    if len(sys.argv) > 1 and sys.argv[1] == '--fast':
        N = 1500

    t_start = time.time()
    experiment_Z1()
    experiment_Z2(N=N)
    experiment_Z3(N=N)

    total = time.time() - t_start
    print(f"\n{'='*70}")
    print(f"TOTAL TIME: {total:.1f}s")
    print(f"{'='*70}")
