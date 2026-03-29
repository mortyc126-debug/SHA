"""
VERIFICATION: Is Π_{k=1}^∞ (1 + (3/4)^k) = e² · √3 ?

Also: derive new formulas from the connection.
Explore q-Pochhammer identities for q = 3/4.
"""
import math

# ============================================================
# PART 1: High-precision numerical computation
# ============================================================
print("=" * 70)
print("PART 1: Numerical verification of Π(1+(3/4)^k) = e²√3")
print("=" * 70)

q = 3/4
product = 1.0
log_product = 0.0

for k in range(1, 1000):
    term = 1 + q**k
    product *= term
    log_product += math.log(term)

print(f"  Π_{{k=1}}^{{999}} (1 + (3/4)^k) = {product:.15f}")
print(f"  ln(Π) = {log_product:.15f}")
print(f"")
print(f"  e² · √3 = {math.e**2 * math.sqrt(3):.15f}")
print(f"  ln(e²√3) = 2 + ln(√3) = {2 + math.log(math.sqrt(3)):.15f}")
print(f"")
print(f"  Difference: {product - math.e**2 * math.sqrt(3):.15e}")
print(f"  Relative error: {abs(product - math.e**2*math.sqrt(3)) / product:.15e}")

# ============================================================
# PART 2: More precise series computation
# ============================================================
print("\n" + "=" * 70)
print("PART 2: Series computation of ln(Π)")
print("=" * 70)

# ln(Π) = Σ ln(1 + q^k) = Σ_k Σ_j (-1)^{j+1} q^{jk}/j
# = Σ_j (-1)^{j+1}/j · q^j/(1-q^j)

series_sum = 0.0
for j in range(1, 200):
    term = ((-1)**(j+1)) / j * q**j / (1 - q**j)
    series_sum += term

print(f"  Series sum (200 terms): {series_sum:.15f}")
print(f"  2 + ln(√3) = {2 + 0.5*math.log(3):.15f}")
print(f"  Difference: {series_sum - 2 - 0.5*math.log(3):.15e}")

# ============================================================
# PART 3: Check term by term
# ============================================================
print("\n" + "=" * 70)
print("PART 3: Decomposition of the series")
print("=" * 70)

# Term j=1: q/(1-q) = (3/4)/(1/4) = 3
# Term j=2: -q²/(2(1-q²)) = -(9/16)/(2·(1-9/16)) = -(9/16)/(2·7/16) = -9/14
# Term j=3: q³/(3(1-q³)) = (27/64)/(3·(1-27/64)) = (27/64)/(3·37/64) = 27/111 = 9/37

print("  Individual terms (-1)^{j+1}/j · q^j/(1-q^j):")
partial = 0.0
for j in range(1, 20):
    term = ((-1)**(j+1)) / j * q**j / (1 - q**j)
    partial += term
    print(f"    j={j:2d}: {term:+.10f}  partial={partial:.10f}")

print(f"\n  Target: 2 + ln(√3) = {2 + 0.5*math.log(3):.10f}")

# ============================================================
# PART 4: Separate into known parts
# ============================================================
print("\n" + "=" * 70)
print("PART 4: Algebraic decomposition attempt")
print("=" * 70)

# Σ (-1)^{j+1}/j · q^j/(1-q^j) = Σ (-1)^{j+1}/j · Σ_{m=1}^∞ q^{jm}
# = Σ_m Σ_j (-1)^{j+1}/j · q^{jm}
# = Σ_m ln(1 + q^m)  ... circular!

# Alternative: use partial fractions
# q^j/(1-q^j) = q^j + q^{2j} + q^{3j} + ...
# So Σ_j (-1)^{j+1}/j · q^j/(1-q^j) = Σ_j (-1)^{j+1}/j · Σ_m q^{jm}
# = Σ_m Σ_j (-1)^{j+1}/j · (q^m)^j
# = Σ_m ln(1 + q^m)
# ... which is circular (= original sum)

# Try Jacobi triple product or theta function identity
# Π_{k=1}^∞ (1 + q^k) = Π (1+q^k) = related to θ functions

# Jacobi triple product: Π_{n=1}^∞ (1-q^{2n})(1+q^{2n-1}z)(1+q^{2n-1}/z) = Σ q^{n²} z^n
# At z=1: Π (1-q^{2n})(1+q^{2n-1})² = θ₃(q)

# For our product: Π (1+q^k) = Π (1+q^{2k}) · Π (1+q^{2k-1})
# = [Π (1+q^{2k})] · [Π (1+q^{2k-1})]

# Known identity: Π_{k=1}^∞ (1+q^k) = Π_{k=1}^∞ 1/(1-q^{2k-1})
# Because (1-q^{2k}) = (1-q^k)(1+q^k), so Π(1+q^k) = Π(1-q^{2k})/Π(1-q^k)

# Π(1-q^{2k})/Π(1-q^k) where the product is k=1..∞
# = Π_{k=1}^∞ (1-q^{2k}) / Π_{k=1}^∞ (1-q^k)
# Numerator = (q²;q²)_∞, Denominator = (q;q)_∞

# So: Π(1+q^k) = (q²;q²)_∞ / (q;q)_∞

print("  Using identity: Π(1+q^k) = (q²;q²)_∞ / (q;q)_∞")

# Compute (q;q)_∞ = Π(1-q^k)
euler_q = 1.0
for k in range(1, 1000):
    euler_q *= (1 - q**k)
print(f"  (q;q)_∞ = Π(1-(3/4)^k) = {euler_q:.15f}")

# Compute (q²;q²)_∞ = Π(1-q^{2k})
q2 = q**2  # = 9/16
euler_q2 = 1.0
for k in range(1, 1000):
    euler_q2 *= (1 - q2**k)
print(f"  (q²;q²)_∞ = Π(1-(9/16)^k) = {euler_q2:.15f}")

ratio = euler_q2 / euler_q
print(f"  Ratio = (q²;q²)/(q;q) = {ratio:.15f}")
print(f"  Our product = {product:.15f}")
print(f"  Match: {abs(ratio - product) < 1e-10}")

# Now: is (q²;q²)/(q;q) = e²√3 for q=3/4?
print(f"\n  Is (q²;q²)/(q;q) = e²√3 ?")
print(f"  e²√3 = {math.e**2 * math.sqrt(3):.15f}")
print(f"  Ratio = {ratio:.15f}")
print(f"  EQUAL? {abs(ratio - math.e**2*math.sqrt(3)) < 1e-6}")

# ============================================================
# PART 5: Try other closed forms
# ============================================================
print("\n" + "=" * 70)
print("PART 5: Search for closed form")
print("=" * 70)

target = product
print(f"  Target value: {target:.15f}")

# Try various combinations
candidates = {
    "e²√3": math.e**2 * math.sqrt(3),
    "e²·3^{1/2}": math.e**2 * 3**0.5,
    "4π/√3": 4*math.pi/math.sqrt(3),
    "π²/√3": math.pi**2/math.sqrt(3),
    "e^{5/2}": math.e**2.5,
    "e²·√3·(1+ε)": math.e**2 * math.sqrt(3),  # same as first
    "3·e": 3*math.e,
    "4·π": 4*math.pi,
    "e²+e": math.e**2 + math.e,
    "2π²/3": 2*math.pi**2/3,
    "7π/√(6)": 7*math.pi/math.sqrt(6),
    "√(2)·e²": math.sqrt(2)*math.e**2,
    "12·ln(3)": 12*math.log(3),
    "e^{2+ln√3}": math.exp(2+math.log(math.sqrt(3))),
    "e^{(5+ln3)/2}": math.exp((5+math.log(3))/2),
}

print(f"  {'Expression':>25} | {'Value':>18} | {'Difference':>15} | {'Rel Error':>12}")
print("  " + "-" * 80)
for name, val in sorted(candidates.items(), key=lambda x: abs(x[1]-target)):
    diff = val - target
    rel = abs(diff/target)
    marker = " ***" if rel < 1e-6 else " **" if rel < 1e-3 else ""
    print(f"  {name:>25} | {val:>18.12f} | {diff:>+15.2e} | {rel:>12.2e}{marker}")

# ============================================================
# PART 6: Test for OTHER q values
# ============================================================
print("\n" + "=" * 70)
print("PART 6: Π(1+q^k) for various q, check for e^a·b^c patterns")
print("=" * 70)

for q_val, q_name in [(1/2, "1/2"), (1/3, "1/3"), (2/3, "2/3"), (3/4, "3/4"), (1/4, "1/4"), (3/8, "3/8")]:
    prod = 1.0
    for k in range(1, 2000):
        prod *= (1 + q_val**k)
    ln_prod = math.log(prod)
    print(f"  q={q_name}: Π = {prod:.10f}, ln = {ln_prod:.10f}")

    # Check if ln(Π) = a + b·ln(c) for small integers
    # ln(Π) for q=1/2: should be related to partition function
    for a in range(0, 5):
        for b_num in range(1, 5):
            for b_den in range(1, 4):
                for c in [2, 3, 5, 7]:
                    test = a + (b_num/b_den)*math.log(c)
                    if abs(test - ln_prod) < 0.001:
                        print(f"    MATCH: ln(Π) ≈ {a} + ({b_num}/{b_den})·ln({c}) = {test:.6f} (err={test-ln_prod:.6f})")
