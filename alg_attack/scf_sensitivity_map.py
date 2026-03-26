#!/usr/bin/env python3
"""
SCF: КАРТА ЧУВСТВИТЕЛЬНОСТИ — где в gap sensitivity ≈ 1?

Wang работает при sensitivity = 1 (δe линейно зависит от δW).
В gap sensitivity — случайное число.
Вопрос: для каких (раунд, слово, бит) sensitivity близка к 1?

sensitivity(r, word, bit) = (e[r+1](W + 2^bit·e_word) - e[r+1](W)) / 2^bit

Если sens = 1: Newton сходится за 1 шаг для этого параметра.
Если sens = odd: Newton сходится, но с большой коррекцией.
Если sens = even: Newton не сходится (необратимо mod 2^32).
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
def sub32(a,b): return (a-b)&MASK32
def hw(x): return bin(x).count('1')

def expand_real(W16):
    W=list(W16)
    for i in range(16,64):
        W.append(add32(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    return W

def sha_state_at_r(W16, R):
    W=expand_real(W16); s=list(IV)
    for r in range(R):
        a,b,c,d,e,f,g,h=s
        T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add32(Sig0(a),Maj(a,b,c))
        s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
    return s

def v2(x):
    """2-adic valuation of x."""
    if x == 0: return 32
    v = 0
    while (x >> v) & 1 == 0: v += 1
    return v


# ============================================================
# EXP 1: Sensitivity map — sens(r, word) for δe
# ============================================================
def exp1_sensitivity_map(N_bases):
    print("="*70)
    print("EXP 1: SENSITIVITY MAP — sens(r, word) = δe / δW[word]")
    print("  sens=1: Wang-like. sens=odd: Newton OK. sens=even: bad.")
    print("="*70)

    # For each (round, word), measure sensitivity across N_bases
    sens_eq_1 = [[0]*16 for _ in range(64)]
    sens_odd = [[0]*16 for _ in range(64)]
    sens_v2_sum = [[0]*16 for _ in range(64)]

    for trial in range(N_bases):
        W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

        for target_r in range(17, 52):
            s_base = sha_state_at_r(W, target_r + 1)

            for word in range(16):
                W_test = list(W)
                W_test[word] = add32(W[word], 1)
                s_test = sha_state_at_r(W_test, target_r + 1)
                sens = sub32(s_test[4], s_base[4])

                if sens == 1:
                    sens_eq_1[target_r][word] += 1
                if sens & 1:
                    sens_odd[target_r][word] += 1
                sens_v2_sum[target_r][word] += v2(sens)

    # Summary: which (round, word) pairs have sens=1 most often?
    print(f"\n  P(sensitivity = 1) for δe — {N_bases} bases")
    print(f"  {'r':>3} | best_word P(s=1) | avg P(s=odd) | avg v2(s)")
    print("  " + "-"*55)

    any_nonzero = False
    for r in range(17, 52):
        best_w = max(range(16), key=lambda w: sens_eq_1[r][w])
        p_eq1 = sens_eq_1[r][best_w] / N_bases
        avg_odd = sum(sens_odd[r][w] for w in range(16)) / (16 * N_bases)
        avg_v2 = sum(sens_v2_sum[r][w] for w in range(16)) / (16 * N_bases)

        marker = ""
        if p_eq1 > 0:
            marker = f" ★ W[{best_w}]"
            any_nonzero = True

        if r < 22 or r > 47 or p_eq1 > 0 or r % 5 == 0:
            print(f"  {r:3d} | {p_eq1:16.4f} | {avg_odd:12.4f} | {avg_v2:9.2f}{marker}")

    if not any_nonzero:
        print(f"\n  ★ P(sens=1) = 0 for ALL (round, word) pairs in gap!")
        print(f"    This confirms: Wang-like Newton is IMPOSSIBLE in the gap.")


# ============================================================
# EXP 2: v2(sensitivity) distribution — how "invertible" is it?
# ============================================================
def exp2_v2_distribution(N_bases):
    print("\n" + "="*70)
    print("EXP 2: v2(SENSITIVITY) DISTRIBUTION IN THE GAP")
    print("  v2=0 (odd): Newton step exists (invertible mod 2^32)")
    print("  v2>0 (even): Newton step doesn't exist for full precision")
    print("  v2=32 (zero): no sensitivity at all")
    print("="*70)

    v2_counts = [0] * 33
    total = 0

    for trial in range(N_bases):
        W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        for r in range(17, 52):
            s_base = sha_state_at_r(W, r+1)
            for word in range(16):
                W_test = list(W)
                W_test[word] = add32(W[word], 1)
                s_test = sha_state_at_r(W_test, r+1)
                sens = sub32(s_test[4], s_base[4])
                v = v2(sens)
                v2_counts[v] += 1
                total += 1

    print(f"\n  v2 | count | fraction | meaning")
    print("  " + "-"*50)
    for v in range(min(10, 33)):
        if v2_counts[v] > 0:
            frac = v2_counts[v] / total
            meaning = "ODD (invertible)" if v==0 else f"divisible by 2^{v}"
            print(f"  {v:2d} | {v2_counts[v]:6d} | {frac:.4f}   | {meaning}")

    p_odd = v2_counts[0] / total
    print(f"\n  P(sensitivity is odd) = {p_odd:.4f}")
    print(f"  Expected (random): 0.5000")
    print(f"  → {'Normal (≈ random)' if abs(p_odd-0.5) < 0.05 else '★ BIASED!'}")


# ============================================================
# EXP 3: ADDITIVE sensitivity — δW = specific value, not +1
# ============================================================
def exp3_additive_target(N_bases):
    print("\n" + "="*70)
    print("EXP 3: TARGETED SENSITIVITY")
    print("  For gap round r: find δW[word] that makes δe[r+1]=0")
    print("  Method: δW = -δe_partial / sensitivity (Newton step)")
    print("  Check: does the Newton step ACTUALLY zero δe?")
    print("="*70)

    success_count = 0
    near_count = 0
    total = 0

    for trial in range(N_bases):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = list(W1)
        W2[0] ^= 0x80000000  # Create a difference

        for r in range(20, 45, 5):
            s1 = sha_state_at_r(W1, r+1)
            s2 = sha_state_at_r(W2, r+1)
            de = sub32(s2[4], s1[4])
            if de == 0: continue

            total += 1

            # Try to correct via each word
            for word in range(16):
                # Sensitivity of W2[word]
                W2_test = list(W2)
                W2_test[word] = add32(W2[word], 1)
                s2_test = sha_state_at_r(W2_test, r+1)
                sens = sub32(s2_test[4], s2[4])

                if sens & 1 == 0: continue  # Not invertible

                inv = pow(sens, -1, 1 << 32)
                correction = ((-de) * inv) & MASK32

                W2_corrected = list(W2)
                W2_corrected[word] = add32(W2[word], correction)
                s2_corr = sha_state_at_r(W2_corrected, r+1)
                de_corr = sub32(s2_corr[4], s1[4])

                if de_corr == 0:
                    success_count += 1
                    if trial < 3:
                        print(f"  r={r}, W[{word}]: δe ZEROED! ★")
                    break
                elif hw(de_corr) <= 5:
                    near_count += 1
                    break

    print(f"\n  Total gap corrections attempted: {total}")
    print(f"  Exact zeros (δe=0): {success_count} ({success_count/max(total,1)*100:.1f}%)")
    print(f"  Near zeros (HW≤5): {near_count} ({near_count/max(total,1)*100:.1f}%)")
    print(f"  Expected if Newton exact: 100%")
    print(f"  Expected if random: ~0%")

    if success_count > 0:
        print(f"\n  ★ Newton-Z WORKS for {success_count/max(total,1)*100:.1f}% of gap rounds!")
        print(f"    (Wang works for 100% of forward rounds)")
        print(f"    Difference: carry + Ch/Maj add nonlinear remainder")


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 50

    exp1_sensitivity_map(min(N, 20))
    exp2_v2_distribution(min(N, 30))
    exp3_additive_target(min(N, 50))

    print("\n" + "="*70)
    print("ИТОГ: КАРТА ЧУВСТВИТЕЛЬНОСТИ")
    print("="*70)
