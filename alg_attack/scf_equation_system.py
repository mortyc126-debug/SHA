#!/usr/bin/env python3
"""
SCF: SHA-256 КАК СИСТЕМА УРАВНЕНИЙ — не оптимизация, а решение.

Стена HW=95 — предел ОПТИМИЗАЦИИ. Нужно РЕШЕНИЕ.

Факт: Newton rank=32 для каждого раунда. Система ИМЕЕТ решение
в линейном приближении. Проблема: нелинейный остаток (carry).

Новый подход: ПОСЛЕДОВАТЕЛЬНАЯ ЭЛИМИНАЦИЯ.
Раунд за раундом, от конца к началу:
1. На раунде 63: решаем δe[64]=0 → получаем уравнение на W[0..15]
2. На раунде 62: решаем δe[63]=0 → ещё одно уравнение
3. ...
4. На раунде 0: решаем δe[1]=0 → последнее уравнение

64 уравнения на 16 неизвестных (W[0..15]).
Система переопределена (64 > 16), но уравнения ЗАВИСИМЫ
(через shift register и schedule).

Вопрос: каков РЕАЛЬНЫЙ ранг этой системы?
Если ранг < 16 → есть свободные переменные → решение существует!
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

def sha_compress(W16, iv=None):
    if iv is None: iv=IV
    W=expand_real(W16); s=list(iv)
    for r in range(64):
        a,b,c,d,e,f,g,h=s
        T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add32(Sig0(a),Maj(a,b,c))
        s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
    return [add32(iv[i],s[i]) for i in range(8)]


# ============================================================
# Tool: FULL JACOBIAN — d(hash)/d(W[word]) over Z_{2^32}
# ============================================================
def full_z_jacobian(W1, W2):
    """Compute 8×16 Jacobian: ∂H_diff[reg] / ∂W2[word] over Z_{2^32}.
    Each entry = sensitivity (32-bit integer)."""
    H1 = sha_compress(W1)
    H2_base = sha_compress(W2)
    dH_base = [sub32(H2_base[i], H1[i]) for i in range(8)]

    J = [[0]*16 for _ in range(8)]
    for w in range(16):
        W2_test = list(W2)
        W2_test[w] = add32(W2[w], 1)
        H2_test = sha_compress(W2_test)
        for reg in range(8):
            J[reg][w] = sub32(sub32(H2_test[reg], H1[reg]), dH_base[reg])
    return J, dH_base


# ============================================================
# EXP 1: Rank of the full Z-Jacobian
# ============================================================
def exp1_z_rank(N):
    print("="*70)
    print("EXP 1: RANK OF Z-JACOBIAN (hash sensitivity)")
    print("  J[reg][word] = ∂(H2[reg]-H1[reg]) / ∂(W2[word])")
    print("  8 equations (hash registers) × 16 unknowns (message words)")
    print("="*70)

    for trial in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = list(W1); W2[0] ^= 0x80000000

        J, dH = full_z_jacobian(W1, W2)

        # Count invertible entries
        n_odd = sum(1 for reg in range(8) for w in range(16) if J[reg][w] & 1)
        n_zero = sum(1 for reg in range(8) for w in range(16) if J[reg][w] == 0)

        # Can we solve for one register?
        # H_diff[reg] = 0 → need Σ J[reg][w] * ΔW[w] = -dH[reg]
        # This is ONE equation in 16 unknowns → always solvable if any J[reg][w] is odd

        solvable_regs = 0
        for reg in range(8):
            has_odd = any(J[reg][w] & 1 for w in range(16))
            if has_odd: solvable_regs += 1

        if trial < 3:
            print(f"\n  Trial {trial}: dH_hw={sum(hw(d) for d in dH)}")
            print(f"    Odd sensitivities: {n_odd}/128")
            print(f"    Zero sensitivities: {n_zero}/128")
            print(f"    Solvable registers: {solvable_regs}/8")

            # Show J structure
            print(f"    J v2-valuation map (0=odd=invertible):")
            for reg in range(8):
                vals = []
                for w in range(16):
                    v = 0
                    x = J[reg][w]
                    if x == 0:
                        v = 32
                    else:
                        while (x >> v) & 1 == 0: v += 1
                    vals.append(v)
                print(f"      H[{reg}]: " + " ".join(f"{v:2d}" for v in vals))


# ============================================================
# EXP 2: Sequential elimination — solve register by register
# ============================================================
def exp2_sequential_elimination(N):
    print("\n" + "="*70)
    print("EXP 2: SEQUENTIAL ELIMINATION")
    print("  Solve H[0]=0 using W[0], then H[1]=0 using W[1], etc.")
    print("  8 registers, 16 words → 8 free words after solving")
    print("="*70)

    for trial in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = list(W1); W2[0] ^= 0x80000000

        H1 = sha_compress(W1)
        init_hw = sum(hw(sha_compress(W2)[i] ^ H1[i]) for i in range(8))

        # Solve register by register
        for target_reg in range(8):
            H2 = sha_compress(W2)
            dH = sub32(H2[target_reg], H1[target_reg])
            if dH == 0:
                continue

            # Find best word to correct this register
            best_word = -1; best_residual = hw(dH)
            for w in range(16):
                W2_test = list(W2)
                W2_test[w] = add32(W2[w], 1)
                H2_test = sha_compress(W2_test)
                sens = sub32(sub32(H2_test[target_reg], H1[target_reg]), dH)

                if sens & 1:  # Invertible
                    inv = pow(sens, -1, 1<<32)
                    corr = (sub32(0, dH) * inv) & MASK32
                    W2_c = list(W2); W2_c[w] = add32(W2[w], corr)
                    if W2_c == W1: continue
                    H2_c = sha_compress(W2_c)
                    resid = hw(sub32(H2_c[target_reg], H1[target_reg]))
                    if resid < best_residual:
                        best_residual = resid
                        best_word = w
                        best_W2 = list(W2_c)

            if best_word >= 0:
                W2 = best_W2

        H2_final = sha_compress(W2)
        final_hw = sum(hw(H2_final[i]^H1[i]) for i in range(8))
        n_diff = sum(1 for i in range(16) if W1[i] != W2[i])

        per_reg = [hw(H2_final[i]^H1[i]) for i in range(8)]

        if trial < 5:
            print(f"\n  Trial {trial}: {init_hw} → {final_hw} ({n_diff} words changed)")
            print(f"    Per-register: " + " ".join(f"H{i}={per_reg[i]}" for i in range(8)))


# ============================================================
# EXP 3: Iterative elimination — solve, recheck, resolve
# ============================================================
def exp3_iterative(N):
    print("\n" + "="*70)
    print("EXP 3: ITERATIVE ELIMINATION — solve, recheck, repeat")
    print("="*70)

    results = []

    for trial in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = list(W1); W2[0] ^= 0x80000000
        H1 = sha_compress(W1)

        best_hw = sum(hw(sha_compress(W2)[i]^H1[i]) for i in range(8))

        for outer in range(20):
            H2 = sha_compress(W2)

            # Find worst register
            reg_hw = [hw(H2[i]^H1[i]) for i in range(8)]
            target = max(range(8), key=lambda i: reg_hw[i])
            dH = sub32(H2[target], H1[target])
            if dH == 0: continue

            # Try all 16 words, pick one that best fixes target WITHOUT hurting others
            best_word = -1; best_total = sum(reg_hw)

            for w in range(16):
                W2_test = list(W2); W2_test[w] = add32(W2[w], 1)
                H2_test = sha_compress(W2_test)
                sens = sub32(sub32(H2_test[target], H1[target]), dH)
                if sens == 0 or sens & 1 == 0: continue

                inv = pow(sens, -1, 1<<32)
                corr = (sub32(0, dH) * inv) & MASK32
                W2_c = list(W2); W2_c[w] = add32(W2[w], corr)
                if W2_c == W1: continue

                H2_c = sha_compress(W2_c)
                total = sum(hw(H2_c[i]^H1[i]) for i in range(8))
                if total < best_total:
                    best_total = total
                    best_word = w
                    best_W2 = list(W2_c)

            if best_word >= 0 and best_total < sum(reg_hw):
                W2 = best_W2
                if best_total < best_hw:
                    best_hw = best_total
            else:
                break

        results.append(best_hw)
        if trial < 5:
            H2 = sha_compress(W2)
            per_reg = [hw(H2[i]^H1[i]) for i in range(8)]
            n_diff = sum(1 for i in range(16) if W1[i]!=W2[i])
            print(f"  Trial {trial}: HW={best_hw}, {n_diff} words, "
                  f"regs={per_reg}")

    avg = sum(results)/len(results)
    print(f"\n  Summary: avg={avg:.1f}, min={min(results)}, max={max(results)}")
    return results


# ============================================================
if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10

    print("="*70)
    print("SCF: SHA КАК СИСТЕМА УРАВНЕНИЙ")
    print("="*70)

    exp1_z_rank(3)
    exp2_sequential_elimination(5)
    results = exp3_iterative(N)

    print(f"\n{'='*70}")
    print(f"ИТОГ")
    print(f"  Iterative elimination: avg={sum(results)/len(results):.1f} min={min(results)}")
    print(f"  SA baseline:          avg≈101 min≈95")
    print(f"  Advantage: {101-sum(results)/len(results):+.1f}")
