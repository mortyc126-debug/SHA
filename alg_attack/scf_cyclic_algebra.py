#!/usr/bin/env python3
"""
SCF: ЦИКЛИЧЕСКАЯ АЛГЕБРА R = GF(2)[x]/(x^32 + 1)

Проблема: Σ-rotations разрушают p-адическую башню.
Решение: работать в алгебре, где ROTR — линейная операция.

В R = GF(2)[x]/(x^32 + 1):
  - 32-bit слово a = a_0 + a_1·x + ... + a_31·x^31
  - ROTR(a, r) = x^r · a  (умножение на x^r — ЛИНЕЙНО!)
  - XOR = сложение в R (линейно)
  - AND = ?? (нужно исследовать)

SHA-256 в этом кольце:
  - Σ0(a) = x^2·a + x^13·a + x^22·a = (x^2 + x^13 + x^22)·a  ЛИНЕЙНО!
  - Σ1(e) = x^6·e + x^11·e + x^25·e = (x^6 + x^11 + x^25)·e  ЛИНЕЙНО!
  - Ch(e,f,g) = e·f + ē·g  (AND = умножение в R? нет, побитовое)
  - Carry(a+b) = ??? в R

Ключевой вопрос: что такое побитовый AND и carry в циклической алгебре?
"""
import os, sys

# ============================================================
# GF(2)[x]/(x^32 + 1) arithmetic
# ============================================================

def poly_from_word(w):
    """32-bit word → polynomial representation (just the int, ops differ)."""
    return w  # Same bit representation, different interpretation

def poly_mul(a, b):
    """Multiply two polynomials in GF(2)[x]/(x^32 + 1).
    This is NOT bitwise AND — it's polynomial multiplication mod x^32+1."""
    result = 0
    for i in range(32):
        if (b >> i) & 1:
            # a * x^i mod (x^32 + 1) = rotate a left by i
            rotated = ((a << i) | (a >> (32 - i))) & 0xFFFFFFFF
            result ^= rotated
    return result

def poly_add(a, b):
    """Add in GF(2)[x]/(x^32+1) = XOR."""
    return a ^ b

def rotr_as_poly(a, r):
    """ROTR(a, r) = x^r * a in the ring."""
    return poly_mul(a, 1 << r)

def rotr_direct(x, n):
    return ((x >> n) | (x << (32 - n))) & 0xFFFFFFFF

def hw(x):
    return bin(x).count('1')


# ============================================================
# EXP 1: Verify ROTR = polynomial multiplication
# ============================================================
def exp1_verify_rotr(N):
    print("="*70)
    print("EXP 1: VERIFY ROTR(a,r) = x^r · a in R")
    print("="*70)

    errors = 0
    for _ in range(N):
        a = int.from_bytes(os.urandom(4), 'big')
        for r in [2, 6, 11, 13, 22, 25]:
            via_rotr = rotr_direct(a, r)
            via_poly = poly_mul(a, 1 << r)
            if via_rotr != via_poly:
                errors += 1

    print(f"  {errors} errors in {N*6} tests → {'✓ VERIFIED' if errors==0 else '✗ FAILED'}")


# ============================================================
# EXP 2: Σ functions as linear operators in R
# ============================================================
def exp2_sigma_linear(N):
    print("\n" + "="*70)
    print("EXP 2: Σ FUNCTIONS AS LINEAR OPERATORS IN R")
    print("  Σ0(a) = (x^2 + x^13 + x^22) · a")
    print("  Σ1(e) = (x^6 + x^11 + x^25) · e")
    print("="*70)

    # The operator polynomial
    S0_poly = (1 << 2) ^ (1 << 13) ^ (1 << 22)  # x^2 + x^13 + x^22
    S1_poly = (1 << 6) ^ (1 << 11) ^ (1 << 25)  # x^6 + x^11 + x^25

    print(f"  S0 operator: 0x{S0_poly:08x} (HW={hw(S0_poly)})")
    print(f"  S1 operator: 0x{S1_poly:08x} (HW={hw(S1_poly)})")

    # Verify Σ0(a) = S0 * a
    errors_0 = 0
    errors_1 = 0
    for _ in range(N):
        a = int.from_bytes(os.urandom(4), 'big')

        sig0_direct = rotr_direct(a, 2) ^ rotr_direct(a, 13) ^ rotr_direct(a, 22)
        sig0_poly = poly_mul(S0_poly, a)
        if sig0_direct != sig0_poly:
            errors_0 += 1

        sig1_direct = rotr_direct(a, 6) ^ rotr_direct(a, 11) ^ rotr_direct(a, 25)
        sig1_poly = poly_mul(S1_poly, a)
        if sig1_direct != sig1_poly:
            errors_1 += 1

    print(f"\n  Σ0: {errors_0} errors / {N} → {'✓' if errors_0==0 else '✗'}")
    print(f"  Σ1: {errors_1} errors / {N} → {'✓' if errors_1==0 else '✗'}")

    # Key property: are S0, S1 invertible in R?
    # S0 is invertible iff gcd(S0_poly, x^32+1) = 1 over GF(2)
    # x^32 + 1 = (x+1)^32 over GF(2)
    # S0 = x^2 + x^13 + x^22, evaluate at x=1: 1+1+1 = 1 ≠ 0
    # So gcd(S0, (x+1)^32) = 1 → S0 is invertible!

    s0_at_1 = hw(S0_poly) % 2
    s1_at_1 = hw(S1_poly) % 2
    print(f"\n  S0(1) = {s0_at_1} → {'invertible in R!' if s0_at_1 == 1 else 'NOT invertible'}")
    print(f"  S1(1) = {s1_at_1} → {'invertible in R!' if s1_at_1 == 1 else 'NOT invertible'}")

    # Compute order of S0 in R* (multiplicative group)
    current = S0_poly
    for order in range(1, 1025):
        if current == 1:
            print(f"  Order of S0 in R*: {order}")
            break
        current = poly_mul(current, S0_poly)
    else:
        print(f"  Order of S0 > 1024")

    current = S1_poly
    for order in range(1, 1025):
        if current == 1:
            print(f"  Order of S1 in R*: {order}")
            break
        current = poly_mul(current, S1_poly)
    else:
        print(f"  Order of S1 > 1024")

    return S0_poly, S1_poly


# ============================================================
# EXP 3: AND operation in R — the key question
# ============================================================
def exp3_and_in_R(N):
    print("\n" + "="*70)
    print("EXP 3: BITWISE AND IN CYCLIC ALGEBRA R")
    print("  AND is NOT polynomial multiplication in R!")
    print("  Poly mul: (Σ a_i x^i)(Σ b_j x^j) = Σ_k (Σ_{i+j≡k} a_i b_j) x^k")
    print("  Bitwise AND: (a AND b)_k = a_k · b_k")
    print("  These differ: poly mul MIXES bit positions, AND does not.")
    print("="*70)

    # Measure: how different is poly_mul from AND?
    diffs = []
    for _ in range(N):
        a = int.from_bytes(os.urandom(4), 'big')
        b = int.from_bytes(os.urandom(4), 'big')

        and_result = a & b
        poly_result = poly_mul(a, b)
        diff = hw(and_result ^ poly_result)
        diffs.append(diff)

    avg_diff = sum(diffs) / len(diffs)
    print(f"\n  HW(AND(a,b) ⊕ poly_mul(a,b)): avg = {avg_diff:.1f} (expected ~16 if random)")

    # Can AND be expressed as a polynomial operation in R?
    # AND(a,b)[k] = a[k] * b[k] (componentwise)
    # This is the Hadamard (Schur) product, NOT the ring product.
    # In the Fourier domain: Hadamard product ↔ convolution.

    # Walsh-Hadamard transform of a 32-bit vector:
    # â[k] = Σ_j (-1)^{bit(k,j)} a[j]
    # In GF(2): the Walsh transform is the same as poly evaluation at roots of x^32+1

    # Key: x^32 + 1 = (x+1)^32 over GF(2), so it has ONE root: x=1.
    # This means R is NOT a product of fields — it's a LOCAL ring!
    # R ≅ GF(2)[x]/(x+1)^32 — a local ring with maximal ideal (x+1)

    print(f"\n  Ring structure of R = GF(2)[x]/(x^32+1):")
    print(f"  x^32 + 1 = (x+1)^32 over GF(2)")
    print(f"  R is a LOCAL RING with maximal ideal m = (x+1)")
    print(f"  R/m ≅ GF(2)")
    print(f"  m^32 = 0 → nilpotent!")

    # The nilradical (x+1) has depth 32.
    # Let y = x+1 (the nilpotent element). Then y^32 = 0.
    # R = GF(2)[y]/(y^32) — a truncated polynomial ring!

    # In terms of y: ROTR(a, r) = (y+1)^r · a
    # (y+1)^r = Σ C(r,k) y^k (binomial expansion)

    print(f"\n  Change of variable: y = x + 1, y^32 = 0")
    print(f"  R = GF(2)[y]/(y^32) — truncated polynomial ring!")
    print(f"  ROTR(a, r) = (1+y)^r · a = Σ C(r,k) y^k · a")
    print(f"\n  Key insight: in the y-basis, ROTR is a BINOMIAL operator!")

    # Compute (1+y)^r for small r
    print(f"\n  (1+y)^r expansion:")
    for r in [2, 6, 11, 13, 22, 25]:
        coeff = 0
        for k in range(32):
            # C(r, k) mod 2 = product of bits (Lucas' theorem)
            c = 1
            rr, kk = r, k
            while kk > 0:
                if (kk & 1) > (rr & 1):
                    c = 0
                    break
                rr >>= 1
                kk >>= 1
            if c:
                coeff |= (1 << k)
        print(f"    (1+y)^{r:2d} = 0x{coeff:08x} (HW={hw(coeff)})")


# ============================================================
# EXP 4: SHA round in y-basis (nilpotent representation)
# ============================================================
def exp4_sha_in_y_basis(N):
    print("\n" + "="*70)
    print("EXP 4: SHA-256 ROUND IN y-BASIS")
    print("  y = x+1, y^32 = 0. R = GF(2)[y]/(y^32)")
    print("  In y-basis: y^k represents 'k-th derivative' of the word")
    print("="*70)

    # Convert word to y-basis
    # a(x) = Σ a_i x^i = Σ a_i (1+y)^i = Σ_k (Σ_i a_i C(i,k)) y^k
    def to_y_basis(word):
        """Convert from x-basis (standard) to y-basis."""
        result = 0
        for k in range(32):
            coeff = 0
            for i in range(32):
                if (word >> i) & 1:
                    # C(i, k) mod 2 via Lucas
                    c = 1
                    ii, kk = i, k
                    while kk > 0:
                        if (kk & 1) > (ii & 1):
                            c = 0; break
                        ii >>= 1; kk >>= 1
                    coeff ^= c
            result |= (coeff << k)
        return result

    def from_y_basis(y_word):
        """Convert from y-basis back to x-basis."""
        # Same transform! (involutory for this specific case)
        return to_y_basis(y_word)

    # Verify round-trip
    errors = 0
    for _ in range(min(N, 100)):
        w = int.from_bytes(os.urandom(4), 'big')
        w2 = from_y_basis(to_y_basis(w))
        if w != w2:
            errors += 1
    print(f"  Round-trip verification: {errors} errors → {'✓' if errors==0 else '✗'}")

    # What does y-basis look like?
    print(f"\n  y-basis representation of sample words:")
    for w in [0x00000001, 0x80000000, 0xFFFFFFFF, 0x6a09e667]:
        yb = to_y_basis(w)
        print(f"    x: 0x{w:08x} (HW={hw(w):2d}) → y: 0x{yb:08x} (HW={hw(yb):2d})")

    # Key: in y-basis, ROTR(a,r) = (1+y)^r · a
    # y^k · a shifts the y-coefficients up by k (truncating at 32)
    # So (1+y)^r · a is a LINEAR combination of {a, y·a, y^2·a, ...}

    # What about Σ0?
    # Σ0(a) = ROTR(a,2) + ROTR(a,13) + ROTR(a,22) (over GF(2))
    # = ((1+y)^2 + (1+y)^13 + (1+y)^22) · a
    # The operator S0_y = (1+y)^2 + (1+y)^13 + (1+y)^22

    def binomial_expand(r):
        """Compute (1+y)^r in GF(2)[y]/(y^32)."""
        result = 0
        for k in range(32):
            c = 1
            rr, kk = r, k
            while kk > 0:
                if (kk & 1) > (rr & 1):
                    c = 0; break
                rr >>= 1; kk >>= 1
            if c:
                result |= (1 << k)
        return result

    S0_y = binomial_expand(2) ^ binomial_expand(13) ^ binomial_expand(22)
    S1_y = binomial_expand(6) ^ binomial_expand(11) ^ binomial_expand(25)

    print(f"\n  Σ-operators in y-basis:")
    print(f"    S0_y = 0x{S0_y:08x} (HW={hw(S0_y)})")
    print(f"    S1_y = 0x{S1_y:08x} (HW={hw(S1_y)})")

    # CRITICAL: what is the y-VALUATION of S0_y?
    # v(S0_y) = min k such that y^k appears
    # This determines: how much "derivative information" is preserved
    v_S0 = 0
    for k in range(32):
        if (S0_y >> k) & 1:
            v_S0 = k; break

    v_S1 = 0
    for k in range(32):
        if (S1_y >> k) & 1:
            v_S1 = k; break

    print(f"    v(S0_y) = {v_S0} (lowest y-power)")
    print(f"    v(S1_y) = {v_S1} (lowest y-power)")

    if v_S0 > 0 or v_S1 > 0:
        print(f"\n  ★ Σ-operator has valuation > 0!")
        print(f"    This means: Σ ANNIHILATES low-order y-components")
        print(f"    The first {max(v_S0,v_S1)} y-coefficients are LOST")
    else:
        print(f"\n  Σ-operators preserve the constant term (y^0)")
        print(f"  No information is lost at the lowest level")

    # What about the CARRY operator in y-basis?
    print(f"\n  Carry (addition mod 2^32) in y-basis:")
    print(f"  a +_Z b in x-basis → ? in y-basis")

    # Measure: carry in y-basis
    y_carry_hw = []
    for _ in range(min(N, 500)):
        a = int.from_bytes(os.urandom(4), 'big')
        b = int.from_bytes(os.urandom(4), 'big')

        sum_real = (a + b) & 0xFFFFFFFF
        sum_xor = a ^ b

        # Carry in x-basis
        carry_x = sum_real ^ sum_xor

        # Carry in y-basis
        carry_y = to_y_basis(carry_x)
        y_carry_hw.append(hw(carry_y))

    avg_y = sum(y_carry_hw) / len(y_carry_hw)
    print(f"  HW(carry in y-basis): avg = {avg_y:.1f}")
    print(f"  HW(carry in x-basis): avg = ~16")

    # y-valuation of carry
    y_val_counts = [0] * 33
    for _ in range(min(N, 1000)):
        a = int.from_bytes(os.urandom(4), 'big')
        b = int.from_bytes(os.urandom(4), 'big')
        carry_x = ((a+b) & 0xFFFFFFFF) ^ (a^b)
        carry_y = to_y_basis(carry_x)
        if carry_y == 0:
            y_val_counts[32] += 1
        else:
            for k in range(32):
                if (carry_y >> k) & 1:
                    y_val_counts[k] += 1
                    break

    print(f"\n  y-valuation distribution of carry:")
    for v in range(6):
        print(f"    v = {v}: {y_val_counts[v]} ({y_val_counts[v]/min(N,1000)*100:.1f}%)")
    print(f"    v = 32 (zero carry): {y_val_counts[32]}")

    return S0_y, S1_y


# ============================================================
if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 500

    print("="*70)
    print("SCF: CYCLIC ALGEBRA R = GF(2)[x]/(x^32+1)")
    print("  Where ROTR = polynomial multiplication (LINEAR!)")
    print("="*70)

    exp1_verify_rotr(min(N, 500))
    S0, S1 = exp2_sigma_linear(min(N, 500))
    exp3_and_in_R(min(N, 500))
    S0_y, S1_y = exp4_sha_in_y_basis(N)

    print("\n" + "="*70)
    print("ИТОГ: ЦИКЛИЧЕСКАЯ АЛГЕБРА")
    print("="*70)
    print(f"""
  R = GF(2)[x]/(x^32+1) = GF(2)[y]/(y^32) где y = x+1.

  Свойства:
  - ROTR(a, r) = (1+y)^r · a — ЛИНЕЙНО (✓)
  - XOR = сложение — ЛИНЕЙНО (✓)
  - Σ0 = (1+y)^2 + (1+y)^13 + (1+y)^22 — ЛИНЕЙНО (✓)
  - AND ≠ poly_mul — РАЗЛИЧНЫЕ операции
  - Carry в y-базисе: HW ≈ {sum(y_carry_hw)/len(y_carry_hw):.0f}

  y-базис интерпретация:
  - y^0 коэффициент = "среднее" (сумма всех бит mod 2)
  - y^k коэффициент = "k-я производная" (k-й порядок разностей)
  - Низкие y-степени = грубая информация
  - Высокие y-степени = тонкая информация

  Σ-операторы в y-базисе:
  - v(S0_y) = {0 if S0_y & 1 else next(k for k in range(32) if (S0_y>>k)&1)}
  - v(S1_y) = {0 if S1_y & 1 else next(k for k in range(32) if (S1_y>>k)&1)}

  ЕСЛИ v(S) > 0: Σ уничтожает информацию на нижних уровнях
  → y-фильтрация = "правильная" p-адическая башня для SHA
  → в отличие от обычной mod 2^k башни, УСТОЙЧИВА к ротациям!
    """)
