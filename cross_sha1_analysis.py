#!/usr/bin/env python3
"""
Cross-analysis: SHA-1 vs SHA-256 structural comparison.
Why SHA-1 fell and SHA-256 stands.
"""

import random
import numpy as np
from itertools import product

MASK = 0xFFFFFFFF

print("=" * 72)
print("CROSS-ANALYSIS: SHA-1 vs SHA-256 — WHY SHA-1 FELL")
print("=" * 72)

# ─── SHA-1 primitives ───────────────────────────────────────────────────

def rotl(x, n):
    return ((x << n) | (x >> (32 - n))) & MASK

def rotr(x, n):
    return ((x >> n) | (x << (32 - n))) & MASK

def sha1_round(state, w, t):
    a, b, c, d, e = state
    if t < 20:
        f = (b & c) | (~b & d) & MASK   # Ch
        k = 0x5A827999
    elif t < 40:
        f = b ^ c ^ d                    # Parity
        k = 0x6ED9EBA1
    elif t < 60:
        f = (b & c) | (b & d) | (c & d) & MASK  # Maj
        k = 0x8F1BBCDC
    else:
        f = b ^ c ^ d                    # Parity
        k = 0xCA62C1D6
    temp = (rotl(a, 5) + f + e + k + w) & MASK
    return [temp, a, rotl(b, 30), c, d]

def expand_sha1(W16):
    W = list(W16)
    for i in range(16, 80):
        W.append(rotl(W[i-3] ^ W[i-8] ^ W[i-14] ^ W[i-16], 1))
    return W

# ─── SHA-256 primitives ─────────────────────────────────────────────────

def sha256_Sigma0(x):
    return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)

def sha256_Sigma1(x):
    return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)

def sha256_sigma0(x):
    return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)

def sha256_sigma1(x):
    return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)

def sha256_Ch(e, f, g):
    return (e & f) ^ (~e & g) & MASK

def sha256_Maj(a, b, c):
    return (a & b) ^ (a & c) ^ (b & c)

SHA256_K = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
    0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
    0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
    0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
    0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
    0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
]

def expand_sha256(W16):
    W = list(W16)
    for i in range(16, 64):
        W.append((sha256_sigma1(W[i-2]) + W[i-7] +
                   sha256_sigma0(W[i-15]) + W[i-16]) & MASK)
    return W

def sha256_round(state, w, t):
    a, b, c, d, e, f, g, h = state
    T1 = (h + sha256_Sigma1(e) + sha256_Ch(e, f, g) + SHA256_K[t] + w) & MASK
    T2 = (sha256_Sigma0(a) + sha256_Maj(a, b, c)) & MASK
    return [(T1 + T2) & MASK, a, b, c, (d + T1) & MASK, e, f, g]


# ═══════════════════════════════════════════════════════════════════════
# ANALYSIS 1: BILINEAR FORM RANK OF ROUND FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════

print("\n" + "=" * 72)
print("ANALYSIS 1: BILINEAR FORM RANK OF BOOLEAN FUNCTIONS")
print("=" * 72)

def compute_bilinear_rank(func, n_vars):
    """
    Given a Boolean function f: {0,1}^n -> {0,1}, compute the rank of its
    bilinear form. The bilinear form matrix B[i][j] is the coefficient of
    x_i * x_j in the ANF (algebraic normal form).
    """
    # Compute ANF using Moebius transform
    N = 1 << n_vars
    anf = [0] * N
    # Evaluate function at all points
    for x in range(N):
        bits = [(x >> i) & 1 for i in range(n_vars)]
        anf[x] = func(*bits)
    # Moebius transform to get ANF coefficients
    for i in range(n_vars):
        for x in range(N):
            if x & (1 << i):
                anf[x] ^= anf[x ^ (1 << i)]

    # Extract bilinear form matrix (degree-2 terms)
    B = np.zeros((n_vars, n_vars), dtype=int)
    for i in range(n_vars):
        for j in range(i + 1, n_vars):
            idx = (1 << i) | (1 << j)
            B[i][j] = anf[idx]
            B[j][i] = anf[idx]

    # Rank over GF(2)
    return gf2_rank(B), B

def gf2_rank(M):
    """Compute rank of matrix over GF(2)."""
    m = M.copy() % 2
    rows, cols = m.shape
    rank = 0
    for col in range(cols):
        pivot = None
        for row in range(rank, rows):
            if m[row, col] == 1:
                pivot = row
                break
        if pivot is None:
            continue
        m[[rank, pivot]] = m[[pivot, rank]]
        for row in range(rows):
            if row != rank and m[row, col] == 1:
                m[row] = (m[row] + m[rank]) % 2
        rank += 1
    return rank

# SHA-1 round functions (single-bit versions)
def Ch_bit(b, c, d):
    return (b & c) ^ ((1 ^ b) & d)

def Parity_bit(b, c, d):
    return b ^ c ^ d

def Maj_bit(b, c, d):
    return (b & c) ^ (b & d) ^ (c & d)

# SHA-256 also uses Ch and Maj (same single-bit functions)

print("\nComputing bilinear form rank for each Boolean function (3 variables):")

for name, func in [("Ch(b,c,d)", Ch_bit), ("Parity(b,c,d)", Parity_bit), ("Maj(b,c,d)", Maj_bit)]:
    rank, B = compute_bilinear_rank(func, 3)
    print(f"\n  {name}:")
    print(f"    Bilinear matrix B =")
    for row in B:
        print(f"      {list(row)}")
    print(f"    Rank over GF(2) = {rank}")

    # Also compute algebraic degree
    N = 1 << 3
    anf = [0] * N
    for x in range(N):
        bits = [(x >> i) & 1 for i in range(3)]
        anf[x] = func(*bits)
    for i in range(3):
        for x in range(N):
            if x & (1 << i):
                anf[x] ^= anf[x ^ (1 << i)]
    deg = 0
    for x in range(N):
        if anf[x]:
            deg = max(deg, bin(x).count('1'))
    print(f"    Algebraic degree = {deg}")
    nl = 0
    for x in range(N):
        bits = [(x >> i) & 1 for i in range(3)]
        nl_check = func(*bits)
    # Compute nonlinearity properly via Walsh spectrum
    max_walsh = 0
    for a in range(N):
        s = 0
        for x in range(N):
            bits = [(x >> i) & 1 for i in range(3)]
            fx = func(*bits)
            lin = 0
            for i in range(3):
                lin ^= ((a >> i) & 1) & ((x >> i) & 1)
            s += (-1) ** (fx ^ lin)
        max_walsh = max(max_walsh, abs(s))
    nl = (N - max_walsh) // 2
    print(f"    Nonlinearity = {nl}")

print("\n  KEY FINDING:")
print("  - Ch and Maj have bilinear rank 2 (quadratic, nonlinear)")
print("  - Parity has bilinear rank 0 (purely linear!)")
print("  - SHA-1 uses Parity in rounds 20-39 and 60-79 (50% of rounds)")
print("  - SHA-256 uses Ch+Maj in ALL 64 rounds (100% nonlinear)")


# ═══════════════════════════════════════════════════════════════════════
# ANALYSIS 2: MESSAGE EXPANSION — LINEAR vs NONLINEAR
# ═══════════════════════════════════════════════════════════════════════

print("\n" + "=" * 72)
print("ANALYSIS 2: MESSAGE EXPANSION LINEARITY")
print("=" * 72)

print("\n--- SHA-1 expansion: XOR + rotation (GF(2)-linear) ---")
print("Testing: if DW[0] = 1, all other DW[i] = 0 for i=1..15")
print("Then DW[t] for t>=16 is DETERMINISTIC (no random input needed)")

# SHA-1 difference propagation in message schedule
DW_sha1 = [0] * 80
DW_sha1[0] = 1  # single-bit difference in W[0]
for i in range(16, 80):
    DW_sha1[i] = rotl(DW_sha1[i-3] ^ DW_sha1[i-8] ^ DW_sha1[i-14] ^ DW_sha1[i-16], 1)

print(f"\n  SHA-1 DW propagation (DW[0]=1):")
for t in range(16, 31):
    hw = bin(DW_sha1[t]).count('1')
    print(f"    DW[{t:2d}] = 0x{DW_sha1[t]:08x}  (HW={hw:2d})")

print(f"\n  Verification: SHA-1 expansion is LINEAR over GF(2)")
print(f"  Testing linearity: expand(A XOR B) == expand(A) XOR expand(B)")

random.seed(42)
A = [random.getrandbits(32) for _ in range(16)]
B = [random.getrandbits(32) for _ in range(16)]
AxB = [(a ^ b) & MASK for a, b in zip(A, B)]

EA = expand_sha1(A)
EB = expand_sha1(B)
EAxB = expand_sha1(AxB)

linear_check = all(EAxB[i] == (EA[i] ^ EB[i]) for i in range(80))
print(f"  Result: expand(A^B) == expand(A)^expand(B) → {linear_check}")

print(f"\n--- SHA-256 expansion: sigma + ADDITION (nonlinear) ---")
print(f"Testing: same DW[0]=1 difference, but now probabilistic")

# SHA-256: need to test with actual values because addition is nonlinear
N_TESTS = 2000
print(f"\n  Sampling {N_TESTS} random base messages, injecting DW[0]=1...")

hw_stats = {t: [] for t in range(16, 31)}
for _ in range(N_TESTS):
    W = [random.getrandbits(32) for _ in range(16)]
    W_prime = list(W)
    W_prime[0] ^= 1  # inject difference

    EW = expand_sha256(W)
    EW_prime = expand_sha256(W_prime)

    for t in range(16, 31):
        dw = (EW[t] ^ EW_prime[t]) & MASK
        hw_stats[t].append(bin(dw).count('1'))

print(f"\n  SHA-256 DW propagation (DW[0]=1, averaged over {N_TESTS} samples):")
for t in range(16, 31):
    hws = hw_stats[t]
    mean_hw = np.mean(hws)
    std_hw = np.std(hws)
    unique_vals = len(set(hws))
    print(f"    DW[{t:2d}]: mean_HW={mean_hw:5.1f} +/- {std_hw:4.1f}  "
          f"(unique HW values: {unique_vals})")

print(f"\n  SHA-256 expansion linearity test:")
EA256 = expand_sha256(A[:16])
EB256 = expand_sha256(B[:16])
EAxB256 = expand_sha256(AxB[:16])
linear256 = all(EAxB256[i] == (EA256[i] ^ EB256[i]) for i in range(64))
print(f"  expand(A^B) == expand(A)^expand(B) → {linear256}")

first_fail = None
for i in range(64):
    if EAxB256[i] != (EA256[i] ^ EB256[i]):
        first_fail = i
        break
if first_fail is not None:
    print(f"  First failure at index {first_fail} (expansion becomes nonlinear at W[{first_fail}])")


# ═══════════════════════════════════════════════════════════════════════
# ANALYSIS 3: DIFFERENTIAL PROPAGATION SPEED
# ═══════════════════════════════════════════════════════════════════════

print("\n" + "=" * 72)
print("ANALYSIS 3: DIFFERENTIAL PROPAGATION SPEED")
print("=" * 72)

N_SAMPLES = 5000
N_ROUNDS = 10

print(f"\nInjecting De[0] = +1 (flip bit 0 of state word e), tracking {N_ROUNDS} rounds")
print(f"Samples: {N_SAMPLES}")

# SHA-1 differential propagation
print(f"\n--- SHA-1 differential propagation ---")
sha1_hw_per_round = []
for r in range(N_ROUNDS):
    hws = []
    for _ in range(N_SAMPLES):
        # Random state and message word
        state = [random.getrandbits(32) for _ in range(5)]
        ws = [random.getrandbits(32) for _ in range(N_ROUNDS)]

        # Perturbed state: flip bit 0 of e (index 4)
        state_p = list(state)
        state_p[4] ^= 1

        # Run r+1 rounds
        s1 = list(state)
        s2 = list(state_p)
        for t in range(r + 1):
            s1 = sha1_round(s1, ws[t], t)  # Use round index for function selection
            s2 = sha1_round(s2, ws[t], t)

        # Compute difference Hamming weight
        diff_hw = sum(bin((s1[i] ^ s2[i]) & MASK).count('1') for i in range(5))
        hws.append(diff_hw)

    mean_hw = np.mean(hws)
    sha1_hw_per_round.append(mean_hw)
    print(f"  Round {r+1:2d}: mean diff HW = {mean_hw:6.2f} / 160 bits "
          f"({mean_hw/160*100:5.1f}%)")

# SHA-256 differential propagation
print(f"\n--- SHA-256 differential propagation ---")
sha256_hw_per_round = []
for r in range(N_ROUNDS):
    hws = []
    for _ in range(N_SAMPLES):
        state = [random.getrandbits(32) for _ in range(8)]
        ws = [random.getrandbits(32) for _ in range(N_ROUNDS)]

        # Flip bit 0 of e (index 4)
        state_p = list(state)
        state_p[4] ^= 1

        s1 = list(state)
        s2 = list(state_p)
        for t in range(r + 1):
            s1 = sha256_round(s1, ws[t], t)
            s2 = sha256_round(s2, ws[t], t)

        diff_hw = sum(bin((s1[i] ^ s2[i]) & MASK).count('1') for i in range(8))
        hws.append(diff_hw)

    mean_hw = np.mean(hws)
    sha256_hw_per_round.append(mean_hw)
    print(f"  Round {r+1:2d}: mean diff HW = {mean_hw:6.2f} / 256 bits "
          f"({mean_hw/256*100:5.1f}%)")

# Normalized comparison
print(f"\n--- Normalized diffusion comparison (% of state affected) ---")
print(f"  {'Round':>5s}  {'SHA-1':>8s}  {'SHA-256':>8s}")
for r in range(N_ROUNDS):
    s1_pct = sha1_hw_per_round[r] / 160 * 100
    s256_pct = sha256_hw_per_round[r] / 256 * 100
    print(f"  {r+1:5d}  {s1_pct:7.1f}%  {s256_pct:7.1f}%")

# Check: how many rounds to reach >40% diffusion?
sha1_40 = next((r+1 for r, hw in enumerate(sha1_hw_per_round) if hw/160 > 0.4), ">10")
sha256_40 = next((r+1 for r, hw in enumerate(sha256_hw_per_round) if hw/256 > 0.4), ">10")
print(f"\n  Rounds to reach >40% diffusion: SHA-1={sha1_40}, SHA-256={sha256_40}")


# ═══════════════════════════════════════════════════════════════════════
# ANALYSIS 4: STRUCTURAL COMPARISON SUMMARY
# ═══════════════════════════════════════════════════════════════════════

print("\n" + "=" * 72)
print("ANALYSIS 4: WHY SHA-1 FELL — STRUCTURAL COMPARISON")
print("=" * 72)

# Compute some derived metrics
sha1_nonlinear_pct = 40 / 80 * 100
sha256_nonlinear_pct = 64 / 64 * 100

print(f"""
┌─────────────────────────┬──────────────────────┬──────────────────────┐
│ Feature                 │ SHA-1                │ SHA-256              │
├─────────────────────────┼──────────────────────┼──────────────────────┤
│ Total rounds            │ 80                   │ 64                   │
│ State width             │ 160 bits (5 words)   │ 256 bits (8 words)   │
│ Nonlinear rounds        │ 40/80 (Ch+Maj only)  │ 64/64 (Ch+Maj)      │
│ Linear rounds           │ 40/80 (Parity)       │ 0/64                │
│ Nonlinear round %       │ {sha1_nonlinear_pct:.0f}%                  │ {sha256_nonlinear_pct:.0f}%                 │
│ Message expansion       │ GF(2)-linear         │ Nonlinear (mod-add)  │
│ Expansion linearity     │ {linear_check!s:20s} │ {linear256!s:20s} │
│ Bilinear rank (Parity)  │ 0                    │ N/A (not used)       │
│ Bilinear rank (Ch)      │ 2                    │ 2                    │
│ Bilinear rank (Maj)     │ 2                    │ 2                    │
│ Best known collision    │ Full 80 rounds       │ 31 rounds            │
│                         │ (SHAttered, 2017)    │ (no full collision)  │
│ Diff propagation in exp │ Deterministic (p=1)  │ Probabilistic (p<1)  │
│ Rounds to 40% diffusion │ {str(sha1_40):20s} │ {str(sha256_40):20s} │
└─────────────────────────┴──────────────────────┴──────────────────────┘
""")

print("KEY STRUCTURAL WEAKNESSES OF SHA-1:")
print("=" * 50)
print("""
1. LINEAR MESSAGE EXPANSION
   SHA-1's expansion W[i] = ROTL(W[i-3]^W[i-8]^W[i-14]^W[i-16], 1)
   is completely linear over GF(2). This means:
   - Differential paths through expansion are DETERMINISTIC
   - Attacker can compute exact difference masks for all W[t]
   - No probabilistic carry propagation to fight against

   SHA-256 uses sigma0/sigma1 with MODULAR ADDITION, making
   expansion nonlinear. Difference propagation becomes probabilistic.

2. PARITY ROUNDS (40 out of 80 = 50%)
   f(b,c,d) = b XOR c XOR d has:
   - Algebraic degree 1 (linear!)
   - Bilinear form rank 0
   - Nonlinearity 0

   These rounds provide ZERO resistance to differential attacks.
   They merely permute and XOR — a linear operation.

   SHA-256 uses Ch+Maj in ALL rounds — always nonlinear.

3. COMBINED EFFECT: LINEAR EXPANSION + LINEAR ROUNDS
   An attacker controlling the message can:
   - Choose input differences that propagate predictably through expansion
   - Navigate through Parity rounds with probability 1
   - Only needs to handle Ch/Maj rounds (40 out of 80)

   This is exactly what Wang et al. (2005) and subsequent attacks exploited.

4. STATE WIDTH
   160 bits vs 256 bits means SHA-1 has a smaller "playground" for
   differences, making collision search fundamentally easier even
   before considering the structural weaknesses above.
""")

print("=" * 72)
print("CONCLUSION")
print("=" * 72)
print("""
SHA-1's fall was not an accident — it was a structural inevitability.
The combination of:
  (a) GF(2)-linear message expansion (deterministic diff propagation)
  (b) 50% linear rounds (zero nonlinear resistance)
  (c) Smaller state width (160 vs 256 bits)

created a hash function where sophisticated differential attacks could
systematically navigate the nonlinear barriers. SHA-256 fixes ALL of
these issues: nonlinear expansion, 100% nonlinear rounds, wider state.

The bilinear rank analysis makes this precise:
  - SHA-1 Parity rounds: rank 0 (linear) — trivially navigable
  - SHA-256 all rounds: rank 2+ (quadratic) — each round is a barrier
""")
