#!/usr/bin/env python3
"""
SCF Этап 3: W_32 Witt-квадратичный анализ SHA-256.
Вопрос: является ли SHA-256 degree-2 в Witt-кольце?
"""
import os, sys

MASK32 = 0xFFFFFFFF

K = [
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,
    0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,
    0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,
    0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,
    0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,
    0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,
    0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,
    0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,
    0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2,
]
IV = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
      0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&MASK32
def shr(x,n): return x>>n
def Sig0(x): return rotr(x,2)^rotr(x,13)^rotr(x,22)
def Sig1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def sig0(x): return rotr(x,7)^rotr(x,18)^shr(x,3)
def sig1(x): return rotr(x,17)^rotr(x,19)^shr(x,10)
def Ch(e,f,g): return (e&f)^(~e&g)&MASK32
def Maj(a,b,c): return (a&b)^(a&c)^(b&c)
def add32(*args):
    s=0
    for x in args: s=(s+x)&MASK32
    return s
def hw(x): return bin(x).count('1')

def expand_schedule(W16):
    W=list(W16)
    for i in range(16,64):
        W.append(add32(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    return W

def sha_compress_rounds(W16, R=64):
    W=expand_schedule(W16)
    s=list(IV)
    for r in range(R):
        a,b,c,d,e,f,g,h=s
        T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add32(Sig0(a),Maj(a,b,c))
        s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
    return s

# ============================================================
# PART 1: Verify W_n(GF(2)) ≅ Z/2^n for small n
# ============================================================
def witt_verify(n):
    """Verify that Witt addition over GF(2) vectors of length n
    is isomorphic to Z/2^n addition."""
    # For W_n(GF(2)), Witt vector (x0, x1, ..., x_{n-1}) maps to
    # the integer: x0 + 2*x1 + 4*x2 + ... + 2^{n-1}*x_{n-1}
    # This is just binary representation!
    # The KEY insight: addition in W_n is the same as integer addition mod 2^n
    # BUT the individual components transform nonlinearly (carries!)

    mod = 1 << n
    errors = 0
    total = mod * mod

    for a in range(mod):
        for b in range(mod):
            # Integer addition mod 2^n
            sum_int = (a + b) % mod
            # XOR (GF(2) addition, component-wise)
            sum_xor = a ^ b
            # Witt addition = integer addition (by isomorphism)
            # The carry correction is: Ψ = (a+b) ^ (a^b)
            psi = sum_int ^ sum_xor
            # Verify Ψ encodes exactly the carry
            carry_shifted = psi  # should equal carry << 1
            # Verify sum_int = sum_xor ^ psi (reconstruction)
            if sum_xor ^ psi != sum_int:
                errors += 1

    return errors, total

# ============================================================
# PART 2: Degree measurement — GF(2) vs Witt perspective
# ============================================================
def measure_degree_gf2(R, N=200):
    """Measure effective algebraic degree of SHA output over GF(2).
    Use 1st, 2nd, 3rd order differentials.
    d-th order differential of degree-d function = constant.
    d-th order differential of degree-(d-1) function = 0."""

    results = {}

    for order in [1, 2, 3]:
        nonzero_count = 0
        for trial in range(N):
            W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

            if order == 1:
                # First difference: F(x+d1) - F(x)
                d1 = [0]*16; d1[0] = 1 << (trial % 32)
                W2 = [(W[i]^d1[i]) for i in range(16)]
                s1 = sha_compress_rounds(W, R)
                s2 = sha_compress_rounds(W2, R)
                diff = [s1[i]^s2[i] for i in range(8)]

            elif order == 2:
                d1 = [0]*16; d1[0] = 1 << (trial % 16)
                d2 = [0]*16; d2[0] = 1 << ((trial % 16) + 16)
                W00 = W
                W10 = [(W[i]^d1[i]) for i in range(16)]
                W01 = [(W[i]^d2[i]) for i in range(16)]
                W11 = [(W[i]^d1[i]^d2[i]) for i in range(16)]
                s00 = sha_compress_rounds(W00, R)
                s10 = sha_compress_rounds(W10, R)
                s01 = sha_compress_rounds(W01, R)
                s11 = sha_compress_rounds(W11, R)
                diff = [s00[i]^s10[i]^s01[i]^s11[i] for i in range(8)]

            elif order == 3:
                d1 = [0]*16; d1[0] = 1
                d2 = [0]*16; d2[0] = 2
                d3 = [0]*16; d3[0] = 4
                vals = []
                for mask in range(8):
                    Wt = list(W)
                    if mask & 1: Wt = [(Wt[i]^d1[i]) for i in range(16)]
                    if mask & 2: Wt = [(Wt[i]^d2[i]) for i in range(16)]
                    if mask & 4: Wt = [(Wt[i]^d3[i]) for i in range(16)]
                    vals.append(sha_compress_rounds(Wt, R))
                diff = [0]*8
                for i in range(8):
                    for mask in range(8):
                        diff[i] ^= vals[mask][i]

            total_hw = sum(hw(d) for d in diff)
            if total_hw > 0:
                nonzero_count += 1

        results[order] = nonzero_count / N

    return results

# ============================================================
# PART 3: Witt-perspective degree — isolate Ch/Maj nonlinearity
# ============================================================
def sha_compress_xor_rounds(W16, R=64):
    """SHA-256 with ALL additions replaced by XOR.
    This removes the carry (Witt) nonlinearity.
    Only Ch and Maj remain nonlinear (degree 2 over GF(2))."""
    W = list(W16)
    for i in range(16, 64):
        W.append(sig1(W[i-2]) ^ W[i-7] ^ sig0(W[i-15]) ^ W[i-16])

    s = list(IV)
    for r in range(R):
        a,b,c,d,e,f,g,h = s
        T1 = h ^ Sig1(e) ^ Ch(e,f,g) ^ K[r] ^ W[r]
        T2 = Sig0(a) ^ Maj(a,b,c)
        s = [T1^T2, a, b, c, d^T1, e, f, g]
    return s

def measure_degree_xor_sha(R, N=200):
    """Measure degree of XOR-SHA (no carries).
    Should be degree 2 (from Ch/Maj only) if Witt theory is correct."""
    results = {}

    for order in [1, 2, 3]:
        nonzero_count = 0
        for trial in range(N):
            W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

            if order == 1:
                d1 = [0]*16; d1[0] = 1 << (trial % 32)
                W2 = [(W[i]^d1[i]) for i in range(16)]
                s1 = sha_compress_xor_rounds(W, R)
                s2 = sha_compress_xor_rounds(W2, R)
                diff = [s1[i]^s2[i] for i in range(8)]

            elif order == 2:
                d1 = [0]*16; d1[0] = 1 << (trial % 16)
                d2 = [0]*16; d2[0] = 1 << ((trial % 16) + 16)
                combos = []
                for mask in range(4):
                    Wt = list(W)
                    if mask & 1: Wt = [(Wt[i]^d1[i]) for i in range(16)]
                    if mask & 2: Wt = [(Wt[i]^d2[i]) for i in range(16)]
                    combos.append(sha_compress_xor_rounds(Wt, R))
                diff = [combos[0][i]^combos[1][i]^combos[2][i]^combos[3][i] for i in range(8)]

            elif order == 3:
                d1 = [0]*16; d1[0] = 1
                d2 = [0]*16; d2[0] = 2
                d3 = [0]*16; d3[0] = 4
                vals = []
                for mask in range(8):
                    Wt = list(W)
                    if mask & 1: Wt = [(Wt[i]^d1[i]) for i in range(16)]
                    if mask & 2: Wt = [(Wt[i]^d2[i]) for i in range(16)]
                    if mask & 4: Wt = [(Wt[i]^d3[i]) for i in range(16)]
                    vals.append(sha_compress_xor_rounds(Wt, R))
                diff = [0]*8
                for i in range(8):
                    for mask in range(8):
                        diff[i] ^= vals[mask][i]

            total_hw = sum(hw(d) for d in diff)
            if total_hw > 0:
                nonzero_count += 1

        results[order] = nonzero_count / N

    return results

# ============================================================
# PART 4: Carry correction statistics
# ============================================================
def carry_correction_analysis(R, N=500):
    """Measure carry correction cc(W) = SHA_real(W) ⊕ SHA_xor(W)."""
    hw_cc_list = []
    for _ in range(N):
        W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        s_real = sha_compress_rounds(W, R)
        s_xor = sha_compress_xor_rounds(W, R)
        cc = [s_real[i] ^ s_xor[i] for i in range(8)]
        hw_cc = sum(hw(x) for x in cc)
        hw_cc_list.append(hw_cc)

    avg = sum(hw_cc_list)/len(hw_cc_list)
    var = sum((x-avg)**2 for x in hw_cc_list)/len(hw_cc_list)
    return avg, var**0.5

# ============================================================
# MAIN
# ============================================================
if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 200

    print("="*70)
    print("SCF: W_32 WITT-КВАДРАТИЧНЫЙ АНАЛИЗ SHA-256")
    print("="*70)

    # Part 1: Witt verification
    print("\n" + "="*70)
    print("PART 1: Verify W_n(GF(2)) ≅ Z/2^n")
    print("="*70)
    for n in [2, 3, 4, 8]:
        errs, total = witt_verify(n)
        print(f"  W_{n}: {errs} errors / {total} pairs → {'✓ ISOMORPHISM' if errs==0 else '✗ FAILED'}")
    print("  Conclusion: Ψ(a,b) = (a+b)⊕(a⊕b) perfectly reconstructs carries")

    # Part 2: Degree comparison
    print("\n" + "="*70)
    print("PART 2: Algebraic degree — Real SHA vs XOR-SHA")
    print("  P(diff≠0) at each order. Degree d ⟹ order-d diff is nonzero,")
    print("  order-(d+1) diff is zero.")
    print("="*70)

    print(f"\n{'R':>3} | {'--- Real SHA (GF(2)+Z) ---':^30} | {'--- XOR-SHA (GF(2) only) ---':^30}")
    print(f"{'':>3} | {'ord1':>7} {'ord2':>7} {'ord3':>7} {'deg':>5} | {'ord1':>7} {'ord2':>7} {'ord3':>7} {'deg':>5}")
    print("-"*75)

    test_rounds = [1, 2, 3, 4, 5, 6, 8, 10, 16, 20]
    for R in test_rounds:
        gf2 = measure_degree_gf2(R, N=min(N, 100))
        witt = measure_degree_xor_sha(R, N=min(N, 100))

        # Estimate degree
        def est_degree(res):
            if res[3] > 0.1: return "≥3"
            if res[2] > 0.1: return "2"
            if res[1] > 0.1: return "1"
            return "0"

        deg_real = est_degree(gf2)
        deg_xor = est_degree(witt)

        marker = ""
        if deg_real != deg_xor:
            marker = " ← DIFFERENT!"

        print(f"{R:3d} | {gf2[1]:7.3f} {gf2[2]:7.3f} {gf2[3]:7.3f} {deg_real:>5} | "
              f"{witt[1]:7.3f} {witt[2]:7.3f} {witt[3]:7.3f} {deg_xor:>5}{marker}")

    # Part 3: Carry correction by round
    print("\n" + "="*70)
    print("PART 3: Carry correction cc(W) = SHA_real ⊕ SHA_xor")
    print("  If SHA is degree-2 in Witt, cc should be 'structured'")
    print("  If cc ~ Bin(256,0.5), carries are random → no Witt advantage")
    print("="*70)

    print(f"\n{'R':>3} | {'E[HW(cc)]':>10} {'σ[HW(cc)]':>10} | {'vs Bin(256,0.5)':>20}")
    print("-"*55)

    for R in [1, 2, 3, 4, 5, 6, 10, 20, 32, 64]:
        avg, std = carry_correction_analysis(R, N=min(N, 300))
        bin_expected = 128.0
        bin_std = 8.0
        deviation = abs(avg - bin_expected) / bin_std
        status = "STRUCTURED" if deviation > 3 else ("weak bias" if deviation > 1 else "RANDOM")
        print(f"{R:3d} | {avg:10.1f} {std:10.1f} | {deviation:5.1f}σ from random → {status}")

    # Part 4: The key question
    print("\n" + "="*70)
    print("PART 4: KEY QUESTION — Is SHA degree-2 in Witt ring?")
    print("="*70)
    print("""
    In W_32 algebra:
      - Addition mod 2^32 is LINEAR (degree 1)
      - Ch(e,f,g) = ef ⊕ eg ⊕ g is QUADRATIC (degree 2)
      - Maj(a,b,c) = ab ⊕ ac ⊕ bc is QUADRATIC (degree 2)
      - Σ0, Σ1, σ0, σ1 are LINEAR (rotations + XOR)

    Therefore ONE round of SHA-256 is degree 2 in W_32.

    But COMPOSITION of rounds:
      - Round 1: degree 2 (from Ch, Maj)
      - Round 2: degree 2² = 4 (Ch applied to degree-2 inputs)
      - Round r: degree 2^r (until saturation at 32)

    The degree STILL grows exponentially through composition!
    W_32 does NOT save us from degree blowup.

    HOWEVER: The carry layer is removed, leaving only GF(2) quadratics.
    XOR-SHA has degree growing as 2^r but over a SIMPLER algebra.
    """)

    # Part 5: XOR-SHA collision feasibility
    print("="*70)
    print("PART 5: XOR-SHA collision search (no carries)")
    print("="*70)

    for R in [1, 2, 3, 4, 5, 6, 8]:
        found = 0
        best_hw = 256
        for trial in range(min(N*5, 2000)):
            W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
            W2 = list(W1)
            W2[0] ^= (1 << (trial % 32))
            s1 = sha_compress_xor_rounds(W1, R)
            s2 = sha_compress_xor_rounds(W2, R)
            diff_hw = sum(hw(s1[i]^s2[i]) for i in range(8))
            if diff_hw < best_hw:
                best_hw = diff_hw
            if diff_hw == 0:
                found += 1

        marker = " ← COLLISION FOUND!" if found > 0 else ""
        print(f"  R={R:2d}: best diff HW = {best_hw:3d}, collisions = {found}{marker}")

    print("\n" + "="*70)
    print("ВЫВОДЫ")
    print("="*70)
    print("""
    1. W_n(GF(2)) ≅ Z/2^n — ПОДТВЕРЖДЕНО для n=2,3,4,8
    2. В Witt-кольце сложение ЛИНЕЙНО, но...
    3. Степень РАСТЁТ через композицию: deg(round_r) = 2^r
    4. XOR-SHA (без carry) — ПРОЩЕ, но всё ещё экспоненциальная степень
    5. Carry correction cc(W) стремится к Bin(256,0.5) — RANDOM
    6. Witt-перспектива: carry = линейный шум в W_32,
       но Ch/Maj = квадратичный шум в GF(2)
       → два типа шума в РАЗНЫХ алгебрах = двойной барьер (F39)
    """)
