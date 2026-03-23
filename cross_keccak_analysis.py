#!/usr/bin/env python3
"""
Cross-analysis: Keccak-f[200] vs SHA-256 nonlinear structure.
Small samples, fast execution (<2 min).
"""

import random
import time

random.seed(42)

# ─── Keccak-f[200] Implementation ───

def keccak_f200(state_bytes):
    """Keccak-f[200] permutation on 25 bytes (5x5 lane width 8)."""
    A = [[0]*5 for _ in range(5)]
    for x in range(5):
        for y in range(5):
            A[x][y] = state_bytes[5*y + x]

    RC = [0x01, 0x82, 0x8a, 0x00, 0x8b, 0x01, 0x81, 0x09,
          0x8a, 0x88, 0x09, 0x0a, 0x8b, 0x8b, 0x89, 0x03, 0x02, 0x80]

    for rnd in range(18):
        # Theta
        C = [0]*5
        for x in range(5):
            C[x] = A[x][0] ^ A[x][1] ^ A[x][2] ^ A[x][3] ^ A[x][4]
        D = [0]*5
        for x in range(5):
            D[x] = C[(x-1)%5] ^ (((C[(x+1)%5] << 1) | (C[(x+1)%5] >> 7)) & 0xFF)
        for x in range(5):
            for y in range(5):
                A[x][y] ^= D[x]

        # Rho + Pi
        B = [[0]*5 for _ in range(5)]
        rho_offsets = [
            [0,1,62,28,27],
            [36,44,6,55,20],
            [3,10,43,25,39],
            [41,45,15,21,8],
            [18,2,61,56,14]
        ]
        for x in range(5):
            for y in range(5):
                rot = rho_offsets[x][y] % 8
                B[y][(2*x+3*y)%5] = ((A[x][y] << rot) | (A[x][y] >> (8-rot))) & 0xFF

        # Chi
        for x in range(5):
            for y in range(5):
                A[x][y] = B[x][y] ^ ((~B[(x+1)%5][y] & B[(x+2)%5][y]) & 0xFF)

        # Iota
        A[0][0] ^= RC[rnd]

    result = [0]*25
    for x in range(5):
        for y in range(5):
            result[5*y + x] = A[x][y]
    return result


def keccak_f200_rounds(state_bytes, num_rounds):
    """Keccak-f[200] with limited rounds."""
    A = [[0]*5 for _ in range(5)]
    for x in range(5):
        for y in range(5):
            A[x][y] = state_bytes[5*y + x]

    RC = [0x01, 0x82, 0x8a, 0x00, 0x8b, 0x01, 0x81, 0x09,
          0x8a, 0x88, 0x09, 0x0a, 0x8b, 0x8b, 0x89, 0x03, 0x02, 0x80]

    for rnd in range(num_rounds):
        C = [0]*5
        for x in range(5):
            C[x] = A[x][0] ^ A[x][1] ^ A[x][2] ^ A[x][3] ^ A[x][4]
        D = [0]*5
        for x in range(5):
            D[x] = C[(x-1)%5] ^ (((C[(x+1)%5] << 1) | (C[(x+1)%5] >> 7)) & 0xFF)
        for x in range(5):
            for y in range(5):
                A[x][y] ^= D[x]

        B = [[0]*5 for _ in range(5)]
        rho_offsets = [
            [0,1,62,28,27],
            [36,44,6,55,20],
            [3,10,43,25,39],
            [41,45,15,21,8],
            [18,2,61,56,14]
        ]
        for x in range(5):
            for y in range(5):
                rot = rho_offsets[x][y] % 8
                B[y][(2*x+3*y)%5] = ((A[x][y] << rot) | (A[x][y] >> (8-rot))) & 0xFF

        for x in range(5):
            for y in range(5):
                A[x][y] = B[x][y] ^ ((~B[(x+1)%5][y] & B[(x+2)%5][y]) & 0xFF)

        A[0][0] ^= RC[rnd]

    result = [0]*25
    for x in range(5):
        for y in range(5):
            result[5*y + x] = A[x][y]
    return result


def state_to_bits(state_bytes):
    """Convert 25 bytes to 200-bit list."""
    bits = []
    for b in state_bytes:
        for i in range(8):
            bits.append((b >> i) & 1)
    return bits


def bits_to_state(bits):
    """Convert 200-bit list to 25 bytes."""
    state = [0]*25
    for i in range(200):
        if bits[i]:
            state[i // 8] |= (1 << (i % 8))
    return state


# ═══════════════════════════════════════════════════════════
print("=" * 70)
print("CROSS-ANALYSIS: Keccak-f[200] vs SHA-256 Nonlinear Structure")
print("=" * 70)

# ─── Analysis 1: Chi nonlinearity characterization ───
print("\n" + "─" * 70)
print("ANALYSIS 1: Chi Nonlinearity Characterization")
print("─" * 70)

print("\nChi operation: a[x] ^= (~a[x+1]) & a[x+2]")
print("Applied to each row of 5 lanes in the state.")
print()

# Demonstrate Chi on a single 5-element row
print("Chi truth table for 3 consecutive bits (x, x+1, x+2) -> output bit x:")
print("  x  x+1  x+2  |  x ^ (~x+1 & x+2)")
print("  " + "-" * 38)
for a in range(2):
    for b in range(2):
        for c in range(2):
            out = a ^ ((~b & 1) & c)
            print(f"  {a}   {b}    {c}   |  {out}")

# Count nonlinear monomials
print("\nAlgebraic Normal Form (ANF) of Chi output bit:")
print("  out_x = x XOR ((NOT x+1) AND x+2)")
print("        = x XOR (x+2 XOR (x+1 * x+2))")
print("        = x XOR x+2 XOR x+1*x+2")
print("  -> Degree 2 (exactly one AND = one degree-2 monomial)")
print()

# Compare with SHA-256 nonlinear operations
print("SHA-256 Ch(e,f,g) = (e AND f) XOR (NOT e AND g)")
print("  = ef XOR g XOR eg  -> Degree 2")
print()
print("SHA-256 Maj(a,b,c) = (a AND b) XOR (a AND c) XOR (b AND c)")
print("  = ab XOR ac XOR bc -> Degree 2")
print()
print("SHA-256 addition carries: a+b carry = a AND b -> Degree 2")
print("  But multi-bit addition chains carries, reaching higher degree quickly")
print()

# AND gate count per round
print("AND gates per round:")
print(f"  Keccak-f[200]: 5 rows x 5 lanes x 1 AND = 25 AND gates")
print(f"  Keccak-f[1600]: 5 rows x 5 lanes x 64 bits x 1 AND = 320 AND gates (per lane bit)")
# SHA-256: Ch has 2 ANDs per bit * 32 bits = 64, Maj has 3 ANDs * 32 = 96,
# plus additions (multiple carries)
print(f"  SHA-256: Ch(64) + Maj(96) + ~32 carry chains = ~192+ AND gates")

# Empirical: measure nonlinearity via Walsh-Hadamard on small Chi
print("\nEmpirical: Chi nonlinearity on 5-bit row")
# Chi on a 5-bit row
def chi_row(row_5bits):
    out = [0]*5
    for x in range(5):
        out[x] = row_5bits[x] ^ ((~row_5bits[(x+1)%5] & 1) & row_5bits[(x+2)%5])
    return out

# Compute nonlinearity of each output bit as boolean function of 5 inputs
for out_bit in range(5):
    # Build truth table
    tt = []
    for inp in range(32):
        bits = [(inp >> i) & 1 for i in range(5)]
        out = chi_row(bits)
        tt.append(out[out_bit])

    # Walsh-Hadamard transform
    n = 5
    wh = [0] * 32
    for w in range(32):
        s = 0
        for x in range(32):
            lin = 0
            for i in range(n):
                lin ^= ((w >> i) & 1) & ((x >> i) & 1)
            s += (-1)**tt[x] * (-1)**lin
        wh[w] = s

    max_wh = max(abs(v) for v in wh)
    nonlinearity = (32 - max_wh) // 2
    print(f"  Output bit {out_bit}: nonlinearity = {nonlinearity} / 12 (max for 5-bit fn)")


# ─── Analysis 2: Algebraic degree growth ───
print("\n" + "─" * 70)
print("ANALYSIS 2: Algebraic Degree Growth Through Rounds")
print("─" * 70)

print("\nTheoretical degree bound after k rounds of Keccak:")
print("  Round 1: degree <= 2")
print("  Round 2: degree <= 2*2 = 4")
print("  Round 3: degree <= 2*4 = 8")
print("  Round k: degree <= min(2^k, n-1) where n = state size")
print()

# Empirical: use higher-order differential to estimate degree
# D^k f = 0 iff degree(f) < k
# For a random function of degree d, D^d f is constant (nonzero), D^{d+1} f = 0

print("Empirical: Higher-order differentials on Keccak-f[200]")
print("Testing: for each round count, find max order k where D^k != 0")
print("(Using 200 random subspace samples per round count, output bit 0)")
print()

def higher_order_diff(func, dim, num_samples=200):
    """
    Test if D^dim of func is zero.
    Pick a random base point and dim random directions.
    Compute the alternating XOR sum over the 2^dim cube.
    Returns fraction of samples where D^dim != 0.
    """
    nonzero_count = 0
    for _ in range(num_samples):
        base = [random.randint(0, 255) for _ in range(25)]
        # Random directions (each is a 200-bit vector packed in 25 bytes)
        directions = []
        for _ in range(dim):
            d = [random.randint(0, 255) for _ in range(25)]
            directions.append(d)

        # Evaluate func at all 2^dim corners of the cube
        xor_sum = 0
        for mask in range(1 << dim):
            point = list(base)
            for j in range(dim):
                if (mask >> j) & 1:
                    point = [point[b] ^ directions[j][b] for b in range(25)]
            out = func(point)
            # Extract bit 0
            xor_sum ^= (out[0] & 1)

        if xor_sum != 0:
            nonzero_count += 1

    return nonzero_count / num_samples


print(f"{'Rounds':<8} {'D^1':>8} {'D^2':>8} {'D^3':>8} {'D^4':>8} {'D^5':>8} {'D^6':>8}")
print("-" * 56)

for num_rounds in range(1, 7):
    t0 = time.time()
    func = lambda s, nr=num_rounds: keccak_f200_rounds(s, nr)
    results = []
    for dim in range(1, 7):
        if dim > 2**num_rounds:
            # Theoretically zero above degree bound
            results.append(0.0)
        else:
            frac = higher_order_diff(func, dim, num_samples=150)
            results.append(frac)
    elapsed = time.time() - t0
    print(f"{num_rounds:<8} {results[0]:>8.3f} {results[1]:>8.3f} {results[2]:>8.3f} "
          f"{results[3]:>8.3f} {results[4]:>8.3f} {results[5]:>8.3f}  ({elapsed:.1f}s)")

print()
print("Interpretation: D^k nonzero fraction > 0 means degree >= k")
print("  After 1 round: degree = 2 (D^2 nonzero, D^3 should drop)")
print("  After 2 rounds: degree = 4 (D^4 nonzero, D^5 should drop)")
print("  Higher rounds: degree saturates toward 199")


# ─── Analysis 3: Diffusion speed ───
print("\n" + "─" * 70)
print("ANALYSIS 3: Diffusion Speed (Avalanche)")
print("─" * 70)

print("\nFlip 1 input bit, measure Hamming distance of output after R rounds")
print("Ideal: 100 bits changed (50% of 200)")
print()

NUM_SAMPLES = 500

print(f"{'Rounds':<8} {'Mean HD':>10} {'% of 200':>10} {'Min':>6} {'Max':>6} {'StdDev':>8}")
print("-" * 52)

for num_rounds in range(1, 19):
    hds = []
    for _ in range(NUM_SAMPLES):
        state = [random.randint(0, 255) for _ in range(25)]
        # Flip one random bit
        flip_byte = random.randint(0, 24)
        flip_bit = random.randint(0, 7)
        state2 = list(state)
        state2[flip_byte] ^= (1 << flip_bit)

        out1 = keccak_f200_rounds(state, num_rounds)
        out2 = keccak_f200_rounds(state2, num_rounds)

        # Hamming distance
        hd = sum(bin(a ^ b).count('1') for a, b in zip(out1, out2))
        hds.append(hd)

    mean_hd = sum(hds) / len(hds)
    min_hd = min(hds)
    max_hd = max(hds)
    std_hd = (sum((h - mean_hd)**2 for h in hds) / len(hds)) ** 0.5
    print(f"{num_rounds:<8} {mean_hd:>10.2f} {mean_hd/200*100:>9.1f}% {min_hd:>6} {max_hd:>6} {std_hd:>8.2f}")

    # Stop early if we've converged and done enough rounds
    if num_rounds >= 6 and abs(mean_hd - 100) < 2:
        print(f"  ... converged, skipping remaining rounds")
        break

print()
print("Note: Keccak's Theta step provides immediate column mixing,")
print("spreading each bit to ~11 positions in a single round.")
print("SHA-256 relies on slower shift-register-style diffusion.")


# ─── Analysis 4: Summary Comparison Table ───
print("\n" + "─" * 70)
print("ANALYSIS 4: Summary Comparison Table")
print("─" * 70)
print()

headers = ["Feature", "SHA-256", "Keccak-f[200]", "Keccak-f[1600] (SHA-3)"]
rows = [
    ["State bits", "256", "200", "1600"],
    ["Construction", "Merkle-Damgard", "Sponge", "Sponge"],
    ["Nonlinear ops", "Ch+Maj+carries", "Chi only", "Chi only"],
    ["Degree per round", "~4+ (carries)", "2 (exact)", "2 (exact)"],
    ["Total rounds", "64", "18", "24"],
    ["AND gates/round", "~192+", "25", "320"],
    ["Diffusion mechanism", "Sigma shifts", "Theta col parity", "Theta col parity"],
    ["Diffusion speed", "Slow (shift reg)", "Fast (~2 rounds)", "Fast (~2 rounds)"],
    ["Rounds to full diff.", "~4-5 rounds", "~2-3 rounds", "~2-3 rounds"],
    ["Design philosophy", "Many weak rounds", "Few strong rounds", "Few strong rounds"],
    ["Best collision atk", "31 of 64 rounds", "—", "6 of 24 rounds"],
    ["Key insight", "ARX complexity", "Wide state + Chi", "Wide state + Chi"],
]

col_widths = [22, 18, 18, 24]
header_line = "| " + " | ".join(h.ljust(w) for h, w in zip(headers, col_widths)) + " |"
sep_line = "|" + "|".join("-" * (w + 2) for w in col_widths) + "|"

print(header_line)
print(sep_line)
for row in rows:
    print("| " + " | ".join(val.ljust(w) for val, w in zip(row, col_widths)) + " |")


# ─── Key Structural Insights ───
print("\n" + "─" * 70)
print("KEY STRUCTURAL INSIGHTS")
print("─" * 70)
print("""
1. NONLINEARITY SOURCE:
   - SHA-256: Multiple nonlinear operations (Ch, Maj, modular addition carries)
     create algebraic complexity through operation diversity.
   - Keccak: Single nonlinear operation (Chi) with degree exactly 2.
     Security comes from ITERATING this through wide-state permutation.

2. DEGREE GROWTH STRATEGY:
   - SHA-256: Degree grows fast per round (~4+ due to carry chains),
     compensated by relatively small state (256 bits).
   - Keccak: Degree doubles each round (2^k), but the huge state (1600 bits)
     means even degree-2 Chi applied across 320 AND gates per round
     creates massive algebraic complexity.

3. DIFFUSION PHILOSOPHY:
   - SHA-256: Slow, local mixing via shifts/rotations. Needs many rounds
     for full diffusion. Each round touches only part of the state.
   - Keccak: Theta provides IMMEDIATE global mixing (column parity).
     Full diffusion in 2-3 rounds. Each round touches ALL state bits.

4. SECURITY MARGIN:
   - SHA-256: 64 rounds, best attack on 31 rounds -> ~2x margin
   - Keccak/SHA-3: 24 rounds, best attack on 6 rounds -> ~4x margin
   - Keccak's wider margin comes from each round being "stronger"

5. WHY BOTH WORK:
   - SHA-256: Many weak rounds * operation diversity = security
   - Keccak: Few strong rounds * wide state * simple nonlinearity = security
   - Both achieve the fundamental goal: make algebraic degree approach
     state size, making linearization attacks infeasible.
""")

print("=" * 70)
print("Analysis complete.")
print("=" * 70)
