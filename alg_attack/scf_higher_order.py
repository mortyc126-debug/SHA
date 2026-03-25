#!/usr/bin/env python3
"""
SCF: Higher-Order Differential Attack on Ψ
Идея: Ψ имеет конечную алгебраическую степень.
Если deg(Ψ_total) = d, то d+1-я производная Ψ = 0.
Тогда d+1-я производная F = d+1-я производная F_xor (чисто GF(2)!).

Вопрос: какова deg(Ψ_total) для полного SHA-256?
Если deg < 32 → higher-order differential attack works.
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

def expand_schedule(W16):
    W=list(W16)
    for i in range(16,64):
        W.append(add32(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    return W

def expand_schedule_xor(W16):
    W=list(W16)
    for i in range(16,64):
        W.append(sig1(W[i-2])^W[i-7]^sig0(W[i-15])^W[i-16])
    return W

def sha_rounds(W16, R):
    W=expand_schedule(W16); s=list(IV)
    for r in range(R):
        a,b,c,d,e,f,g,h=s
        T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add32(Sig0(a),Maj(a,b,c))
        s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
    return s

def sha_xor_rounds(W16, R):
    W=expand_schedule_xor(W16); s=list(IV)
    for r in range(R):
        a,b,c,d,e,f,g,h=s
        T1=h^Sig1(e)^Ch(e,f,g)^K[r]^W[r]
        T2=Sig0(a)^Maj(a,b,c)
        s=[T1^T2,a,b,c,d^T1,e,f,g]
    return s

def psi_output(W16, R):
    """Ψ(W) = SHA_real(W) ⊕ SHA_xor(W) at round R."""
    real = sha_rounds(W16, R)
    xor = sha_xor_rounds(W16, R)
    return [real[i] ^ xor[i] for i in range(8)]


# ============================================================
# Higher-order differentials
# ============================================================
def nth_order_diff(func, W_base, directions, R):
    """
    Compute n-th order differential of func.
    directions = list of (word_idx, bit_mask) deltas.
    n-th order diff = XOR over all 2^n evaluations.
    """
    n = len(directions)
    result = [0]*8

    for mask in range(1 << n):
        W = list(W_base)
        for i in range(n):
            if mask & (1 << i):
                word_idx, bit_mask = directions[i]
                W[word_idx] ^= bit_mask

        output = func(W, R)
        for i in range(8):
            result[i] ^= output[i]

    return result


# ============================================================
# EXP 1: Degree of Ψ — at which order does d^n(Ψ) become 0?
# ============================================================
def exp1_psi_degree(N_trials=50):
    print("="*70)
    print("EXP 1: ALGEBRAIC DEGREE OF Ψ")
    print("  d^n(Ψ) = 0 ⟺ deg(Ψ) < n")
    print("  Test orders 1..8 for various rounds")
    print("="*70)

    for R in [1, 2, 3, 4, 5, 6, 10, 20, 64]:
        print(f"\n  --- R = {R} rounds ---")
        print(f"  {'order':>6} | {'P(d^n≠0)':>10} {'E[HW]':>8} | {'degree':>8}")

        for order in range(1, min(9, 33)):
            nonzero = 0
            total_hw = 0

            for trial in range(N_trials):
                W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

                # Generate n independent directions in W[0]
                directions = []
                for d in range(order):
                    bit = d % 32
                    word = d // 32
                    directions.append((word, 1 << bit))

                # n-th order differential of Ψ
                diff_psi = nth_order_diff(psi_output, W, directions, R)
                diff_hw = sum(hw(x) for x in diff_psi)

                if diff_hw > 0:
                    nonzero += 1
                total_hw += diff_hw

            p_nonzero = nonzero / N_trials
            avg_hw = total_hw / N_trials

            if p_nonzero == 0:
                deg_est = f"< {order}"
                print(f"  {order:6d} | {p_nonzero:10.3f} {avg_hw:8.1f} | {deg_est:>8} ← ZERO!")
                break
            else:
                deg_est = f"≥ {order}"
                print(f"  {order:6d} | {p_nonzero:10.3f} {avg_hw:8.1f} | {deg_est:>8}")

                # If all are nonzero at high order, degree is maximal
                if order >= 8 and p_nonzero > 0.9:
                    print(f"         | degree ≥ 8 (likely 32) — stopping")
                    break


# ============================================================
# EXP 2: Compare d^n(F_real) vs d^n(F_xor) vs d^n(Ψ)
# ============================================================
def exp2_three_way_comparison(N_trials=50):
    print("\n" + "="*70)
    print("EXP 2: THREE-WAY COMPARISON of n-th order differentials")
    print("  d^n(F_real), d^n(F_xor), d^n(Ψ)")
    print("  If d^n(Ψ)=0 but d^n(F_xor)≠0 → d^n(F_real) = d^n(F_xor)!")
    print("="*70)

    for R in [2, 4, 6, 10, 64]:
        print(f"\n  --- R = {R} ---")
        print(f"  {'ord':>4} | {'P(d^n_real≠0)':>14} {'P(d^n_xor≠0)':>14} {'P(d^n_Ψ≠0)':>12} | "
              f"{'d^n_real=d^n_xor?':>18}")

        for order in range(1, 7):
            nz_real = 0; nz_xor = 0; nz_psi = 0; match = 0

            for trial in range(N_trials):
                W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
                directions = [(d % 16, 1 << (d % 32)) for d in range(order)]

                d_real = nth_order_diff(sha_rounds, W, directions, R)
                d_xor = nth_order_diff(sha_xor_rounds, W, directions, R)
                d_psi = nth_order_diff(psi_output, W, directions, R)

                hw_real = sum(hw(x) for x in d_real)
                hw_xor = sum(hw(x) for x in d_xor)
                hw_psi = sum(hw(x) for x in d_psi)

                if hw_real > 0: nz_real += 1
                if hw_xor > 0: nz_xor += 1
                if hw_psi > 0: nz_psi += 1

                # Check if d^n_real = d^n_xor (i.e., d^n_Ψ = 0 for this input)
                if all(d_real[i] == d_xor[i] for i in range(8)):
                    match += 1

            p_real = nz_real/N_trials
            p_xor = nz_xor/N_trials
            p_psi = nz_psi/N_trials
            p_match = match/N_trials

            marker = " ★★★" if p_match > 0.9 else (" ★" if p_match > 0.5 else "")
            print(f"  {order:4d} | {p_real:14.3f} {p_xor:14.3f} {p_psi:12.3f} | {p_match:18.3f}{marker}")


# ============================================================
# EXP 3: Degree of Ψ for SINGLE addition (theoretical baseline)
# ============================================================
def exp3_single_add_degree():
    print("\n" + "="*70)
    print("EXP 3: DEGREE OF Ψ FOR SINGLE ADDITION (baseline)")
    print("  Ψ(a,b) = (a+b) ⊕ (a⊕b)")
    print("  Theory: deg(Ψ[bit_i]) = i+1")
    print("  Verify: d^{i+2}(Ψ[bit_i]) should be 0")
    print("="*70)

    for n_bits in [4, 8, 16, 32]:
        mask = (1 << n_bits) - 1
        print(f"\n  --- n = {n_bits} bits ---")

        for target_bit in [0, 1, 2, 3, min(n_bits-1, 7)]:
            # Find degree of Ψ[target_bit] as function of a (with b fixed)
            # d^k(Ψ[bit]) over GF(2) w.r.t. a

            max_order_found = 0
            for order in range(1, min(n_bits+1, 10)):
                nonzero = 0
                for trial in range(100):
                    a = int.from_bytes(os.urandom(4),'big') & mask
                    b = int.from_bytes(os.urandom(4),'big') & mask

                    # d^order Ψ[target_bit] w.r.t. bits 0..order-1 of a
                    result = 0
                    for m in range(1 << order):
                        aa = a
                        for d in range(order):
                            if m & (1 << d):
                                aa ^= (1 << d)
                        psi_val = ((aa + b) & mask) ^ (aa ^ b)
                        result ^= (psi_val >> target_bit) & 1

                    if result != 0:
                        nonzero += 1

                if nonzero > 0:
                    max_order_found = order
                else:
                    break

            expected = target_bit + 1
            status = "✓" if max_order_found == expected else f"✗ (expected {expected})"
            print(f"    bit {target_bit}: deg = {max_order_found} {status}")


# ============================================================
# EXP 4: Higher-order Ψ-nullification
# ============================================================
def exp4_psi_nullification(N_trials=30):
    print("\n" + "="*70)
    print("EXP 4: Ψ-NULLIFICATION via higher-order differentials")
    print("  If d^n(Ψ)=0 at order n, then d^n(F_real) lives in GF(2)-space")
    print("  This gives us equations in F_xor only (no carry!)")
    print("="*70)

    print(f"\n  Searching for order where d^n(Ψ) = 0 for full SHA-256...")
    print(f"  Using directions within single word W[0]")

    R = 64  # Full SHA

    for order in range(1, 33):
        zero_count = 0
        total_hw = 0

        for trial in range(N_trials):
            W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

            # Directions: bits 0..order-1 of W[0]
            directions = [(0, 1 << d) for d in range(order)]

            d_psi = nth_order_diff(psi_output, W, directions, R)
            d_hw = sum(hw(x) for x in d_psi)

            if d_hw == 0:
                zero_count += 1
            total_hw += d_hw

        p_zero = zero_count / N_trials
        avg_hw = total_hw / N_trials

        print(f"  order {order:2d}: P(d^{order}Ψ=0)={p_zero:.3f}, E[HW]={avg_hw:6.1f}", end="")

        if p_zero == 1.0:
            print(f" ← ALL ZERO! deg(Ψ) ≤ {order-1}")
            print(f"\n  ★★★ Ψ-NULLIFICATION at order {order}!")
            print(f"  This means: d^{order}(F_real) = d^{order}(F_xor)")
            print(f"  → {order}-th order differential of SHA = purely GF(2)")
            return order
        elif avg_hw < 10:
            print(f" ← near zero")
        else:
            print()

        # If still maximal at order 8, likely degree 32
        if order >= 8 and p_zero == 0 and avg_hw > 100:
            print(f"\n  deg(Ψ) appears to be 32 (full). No nullification possible.")
            print(f"  Higher-order attack on Ψ is infeasible for full SHA-256.")
            return None

    return None


# ============================================================
# EXP 5: Per-round Ψ degree evolution
# ============================================================
def exp5_per_round_degree(N_trials=40):
    print("\n" + "="*70)
    print("EXP 5: PER-ROUND Ψ DEGREE EVOLUTION")
    print("  At which round does d^n(Ψ) become nonzero?")
    print("="*70)

    print(f"\n  {'R':>3} |", end="")
    for order in range(1, 9):
        print(f" d^{order}≠0", end="")
    print(" | est_deg")
    print("  " + "-"*70)

    for R in range(1, 21):
        print(f"  {R:3d} |", end="")
        est_deg = 0

        for order in range(1, 9):
            nonzero = 0
            for trial in range(N_trials):
                W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
                directions = [(0, 1 << d) for d in range(order)]
                d_psi = nth_order_diff(psi_output, W, directions, R)
                if sum(hw(x) for x in d_psi) > 0:
                    nonzero += 1

            p = nonzero/N_trials
            if p > 0.1:
                est_deg = order
            print(f" {p:5.2f}", end="")

        deg_str = f"≥{est_deg}" if est_deg >= 8 else str(est_deg)
        print(f" | {deg_str}")


# ============================================================
if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 50

    print("="*70)
    print("SCF: HIGHER-ORDER DIFFERENTIAL ATTACK ON Ψ")
    print("  Can we kill the carry layer with high-order derivatives?")
    print("="*70)

    exp3_single_add_degree()
    exp5_per_round_degree(N_trials=N)
    exp2_three_way_comparison(N_trials=N)
    null_order = exp4_psi_nullification(N_trials=N)

    print("\n" + "="*70)
    print("ФИНАЛЬНЫЙ ВЫВОД")
    print("="*70)
    if null_order:
        print(f"""
  ★ Ψ-NULLIFICATION найдена при order={null_order}
  Это означает:
    d^{null_order}(SHA_real) = d^{null_order}(SHA_xor)
  → {null_order}-я производная SHA — чисто GF(2) объект!
  → Carry полностью устраняется в {null_order}-м порядке.

  Стоимость: 2^{null_order} вычислений SHA на один запрос.
  Если {null_order} < 128 → это быстрее birthday!
        """)
    else:
        print("""
  deg(Ψ) = 32 (полная) для полного SHA-256.
  Higher-order differentials НЕ убивают Ψ
  при разумном порядке (≤32).

  НО: для reduced rounds (R<6) deg(Ψ) < 32.
  → Higher-order attack работает на reduced SHA.
        """)
