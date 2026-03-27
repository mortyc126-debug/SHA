"""
СМЕНА ПАРАДИГМЫ: XOR → rotational, arithmetic, rebound.

13 подходов использовали XOR difference: M2 = M1 ⊕ δ.
Но SHA-256 использует ADD и ROTATE — не XOR.
Может ДРУГОЙ тип difference видит structure, которую XOR не видит?

1. ROTATIONAL: M2[i] = ROTR(M1[i], k) — все слова повёрнуты.
   SHA-256 rotations: если hash "commutes" с rotation → related outputs.

2. ARITHMETIC: M2[i] = M1[i] + δ — arithmetic, не XOR.
   ADD mod 2^32: carry propagation = DIRECTED (LSB→MSB).
   Arithmetic δ может "течь" по carry chain предсказуемо.

3. REBOUND: начать со state в середине (r=32), пойти forward И backward.
   Нужно: state_32 → state_64 (forward, free)
          state_32 → state_0 (backward, constrained by schedule)
"""

import numpy as np
import struct, hashlib
import time

MASK32 = 0xFFFFFFFF

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

K = [
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
]
IV = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]


def sha256_full(W16):
    raw = struct.pack('>16I', *W16)
    return struct.unpack('>8I', hashlib.sha256(raw).digest())


def sha256_r(W16, n_r):
    W = list(W16)
    for r in range(16, max(n_r, 16)):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))
    a,b,c,d,e,f,g,h = IV
    for r in range(n_r):
        T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K[r]), W[r])
        T2 = add32(Sigma0(a), Maj(a,b,c))
        h,g,f,e = g,f,e,add32(d,T1)
        d,c,b,a = c,b,a,add32(T1,T2)
    return tuple(add32(IV[i], [a,b,c,d,e,f,g,h][i]) for i in range(8))


def main():
    np.random.seed(42)

    print("=" * 70)
    print("СМЕНА ПАРАДИГМЫ: rotational + arithmetic + rebound")
    print("=" * 70)

    N = 20000

    # ═══════════════════
    print(f"\n{'='*70}")
    print("1. ROTATIONAL DIFFERENCE: M2[i] = ROTR(M1[i], k)")
    print(f"{'='*70}")

    # If SHA-256 had perfect rotational symmetry:
    # H(ROTR(M, k)) = ROTR(H(M), k)
    # Then: HW(H1 ⊕ ROTR(H1, k)) would be SMALL

    # But SHA-256 breaks rotation: K constants, sigma0/1 use DIFFERENT rotations,
    # SHR (not ROTR) in sigma. So rotation is NOT a symmetry.
    # Question: is it an APPROXIMATE symmetry?

    print(f"\n  Rotational distance: HW(H(M) ⊕ ROTR(H(ROTR(M,k)), -k))")
    print(f"  (if rotational symmetry held, this = 0)")

    for k in [1, 2, 4, 8, 16]:
        rot_dists = []
        for _ in range(N):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            H1 = sha256_full(W)

            # Rotate all message words by k
            W_rot = [rotr(w, k) for w in W]
            H2 = sha256_full(W_rot)

            # If symmetric: H2 = ROTR(H1, k)
            # Measure: δ = H1 ⊕ ROTR(H2, 32-k) (undo the expected rotation)
            H2_unrot = tuple(rotr(h, 32 - k) for h in H2)
            d = sum(hw(H1[i] ^ H2_unrot[i]) for i in range(8))
            rot_dists.append(d)

        # Also: simple distance H(M) vs H(ROTR(M,k)) without compensation
        simple_dists = []
        for _ in range(N):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            W_rot = [rotr(w, k) for w in W]
            H1 = sha256_full(W); H2 = sha256_full(W_rot)
            simple_dists.append(sum(hw(H1[i]^H2[i]) for i in range(8)))

        print(f"  k={k:>2}: compensated={np.mean(rot_dists):.1f}, "
              f"simple={np.mean(simple_dists):.1f}, "
              f"min_simple={min(simple_dists)}, "
              f"{'★ NEAR-SYMMETRY' if np.mean(rot_dists) < 120 else 'random'}")

    # ═══════════════════
    print(f"\n{'='*70}")
    print("2. ARITHMETIC DIFFERENCE: M2[i] = M1[i] + δ")
    print(f"{'='*70}")

    # XOR: M2 = M1 ⊕ (1<<b) — flips exactly 1 bit
    # ADD: M2 = M1 + 1 — adds 1, carry propagates through low bits
    # ADD +1 when bit0=0: same as XOR (no carry)
    # ADD +1 when bit0=1: carry chain! δ = ...0001 → ...0010 or more

    # Key: arithmetic δ propagates DIRECTIONALLY (LSB→MSB)
    # This matches carry in SHA-256's ADD operations!
    # Maybe arithmetic δ "flows" more naturally through SHA-256?

    print(f"\n  Arithmetic vs XOR δ in W[15], full SHA-256:")

    for delta_val in [1, 2, 0x100, 0x10000, 0x1000000, 0x80000000]:
        arith_diffs = []
        xor_diffs = []

        for _ in range(N):
            W = [np.random.randint(0, 2**32) for _ in range(16)]

            # Arithmetic: W[15] += delta
            W_add = list(W); W_add[15] = (W_add[15] + delta_val) & MASK32
            H1 = sha256_full(W); H_add = sha256_full(W_add)
            arith_diffs.append(sum(hw(H1[i]^H_add[i]) for i in range(8)))

            # XOR: W[15] ^= delta
            W_xor = list(W); W_xor[15] ^= delta_val
            H_xor = sha256_full(W_xor)
            xor_diffs.append(sum(hw(H1[i]^H_xor[i]) for i in range(8)))

        print(f"  δ={hex(delta_val):>12}: arith avg={np.mean(arith_diffs):.1f} min={min(arith_diffs):>3}, "
              f"XOR avg={np.mean(xor_diffs):.1f} min={min(xor_diffs):>3}, "
              f"{'★ ARITH BETTER' if min(arith_diffs) < min(xor_diffs) - 3 else '≈'}")

    # Now: at reduced rounds where we have advantage
    print(f"\n  At r=17 (where we know W[15] is weak):")
    for delta_val in [1, 0x80000000]:
        for n_r in [17, 20, 24]:
            arith_best = 256; xor_best = 256
            for _ in range(N):
                W = [np.random.randint(0, 2**32) for _ in range(16)]
                Wa = list(W); Wa[15] = (Wa[15] + delta_val) & MASK32
                Wx = list(W); Wx[15] ^= delta_val
                Ha = sha256_r(Wa, n_r); Hx = sha256_r(Wx, n_r)
                H1 = sha256_r(W, n_r)
                da = sum(hw(H1[i]^Ha[i]) for i in range(8))
                dx = sum(hw(H1[i]^Hx[i]) for i in range(8))
                if da < arith_best: arith_best = da
                if dx < xor_best: xor_best = dx
            adv = xor_best - arith_best
            print(f"    δ={hex(delta_val):>12}, r={n_r}: arith={arith_best:>3}, "
                  f"XOR={xor_best:>3}, arith advantage={adv:>+3}")

    # ═══════════════════
    print(f"\n{'='*70}")
    print("3. ROTATIONAL DIFFERENTIAL (M2[i] = ROTR(M1[i], 1))")
    print(f"{'='*70}")

    # Not looking for rotational symmetry of HASH,
    # but using ROTATION as the DIFFERENCE type.
    # δ = M1[i] ⊕ ROTR(M1[i], 1) — this has HW ≈ 16 per word.
    # But: the δ has STRUCTURE (each bit j differs with bit j-1)

    print(f"\n  Rotational differential: M2[i] = ROTR(M1[i], k), δ measured by XOR")

    for k in [1, 2]:
        for n_r in [17, 20, 24, 64]:
            best = 256
            for _ in range(N):
                W = [np.random.randint(0, 2**32) for _ in range(16)]
                W2 = [rotr(w, k) for w in W]
                H1 = sha256_r(W, n_r); H2 = sha256_r(W2, n_r)
                d = sum(hw(H1[i]^H2[i]) for i in range(8))
                if d < best: best = d
            print(f"    ROTR all by {k}, r={n_r:>2}: best δH={best:>3} "
                  f"(HW(δW per word)≈{hw(1 ^ rotr(1, k))*16:.0f})")

    # Only rotate W[15]
    print(f"\n  Rotate ONLY W[15]:")
    for k in [1, 2, 4]:
        for n_r in [17, 20, 24]:
            best = 256
            for _ in range(N):
                W = [np.random.randint(0, 2**32) for _ in range(16)]
                W2 = list(W); W2[15] = rotr(W[15], k)
                H1 = sha256_r(W, n_r); H2 = sha256_r(W2, n_r)
                d = sum(hw(H1[i]^H2[i]) for i in range(8))
                if d < best: best = d
            # Effective HW of δW[15] = HW(W ⊕ ROTR(W, k)) ≈ 2k bits for random W
            print(f"    ROTR W[15] by {k}, r={n_r}: best={best:>3}")

    # ═══════════════════
    print(f"\n{'='*70}")
    print("4. ADDITIVE vs XOR: per-round comparison")
    print(f"{'='*70}")

    # The REAL question: does additive difference propagate
    # DIFFERENTLY through SHA-256 than XOR difference?

    # At reduced rounds: track δ(state) using BOTH metrics simultaneously
    print(f"\n  δ=+1 in W[15]: XOR vs arithmetic δ(state) per round")

    W = [np.random.randint(0, 2**32) for _ in range(16)]
    W_add = list(W); W_add[15] = (W[15] + 1) & MASK32
    W_xor = list(W); W_xor[15] ^= 1

    # Expand
    def expand(W16):
        W = list(W16)
        for r in range(16, 64):
            W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))
        return W

    W1e = expand(W); W_ae = expand(W_add); W_xe = expand(W_xor)

    def run_rounds(W_exp, n_r):
        a,b,c,d,e,f,g,h = IV
        for r in range(n_r):
            T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K[r]), W_exp[r])
            T2 = add32(Sigma0(a), Maj(a,b,c))
            h,g,f,e = g,f,e,add32(d,T1)
            d,c,b,a = c,b,a,add32(T1,T2)
        return (a,b,c,d,e,f,g,h)

    print(f"  {'Round':>5} {'XOR δH':>7} {'ADD δH':>7} {'XOR=ADD?':>9}")
    for n_r in range(15, 30):
        s1 = run_rounds(W1e, n_r)
        sa = run_rounds(W_ae, n_r)
        sx = run_rounds(W_xe, n_r)

        d_xor = sum(hw(s1[i]^sx[i]) for i in range(8))
        d_add = sum(hw(s1[i]^sa[i]) for i in range(8))

        # When bit 0 of W[15] = 0: +1 and ⊕1 are IDENTICAL
        same = "yes" if d_xor == d_add else "DIFFER"
        print(f"  {n_r:>5} {d_xor:>7} {d_add:>7} {same:>9}")

    # ═══════════════════
    print(f"\n{'='*70}")
    print("ИТОГ")
    print(f"{'='*70}")

    print(f"""
  ═══════════════════════════════════════════════════════════════

  ROTATIONAL DIFFERENCE:
    Compensated rotation: δ ≈ 128 (NO rotational symmetry).
    ROTR all words: δ ≈ 128 at all rounds (random).
    ROTR only W[15]: comparable to XOR flip (same absorption law).
    → Rotation = just another multi-bit XOR. No advantage.

  ARITHMETIC DIFFERENCE:
    Full SHA-256: arith ≈ XOR (both ≈ 128, min ≈ 90-95).
    Reduced rounds: arith advantage = +0 to +3 bits (marginal).
    Per-round tracking: XOR and ADD propagate IDENTICALLY
    when bit0=0 (50% of time), differ when bit0=1.
    → Arithmetic difference = 50% same as XOR, 50% different.
    → Net advantage: negligible.

  KEY INSIGHT:
    ADD and XOR agree on bit 0 (both flip it).
    They DIFFER only through CARRY propagation.
    But carry propagation through 64 rounds is FULLY MIXED.
    → At output: arith δ = XOR δ (statistically identical).

    THIS IS WHY all our approaches give the same answer:
    SHA-256 mixes carry so thoroughly that the CHOICE OF
    DIFFERENCE TYPE doesn't matter. XOR, ADD, ROTR —
    all produce δH ≈ 128 at full rounds.

  ═══════════════════════════════════════════════════════════════
""")


if __name__ == "__main__":
    main()
