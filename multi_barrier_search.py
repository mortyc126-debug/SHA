"""
MULTI-BARRIER NEAR-COLLISION SEARCH

Extend the 8-bit result (4+4 on r=17,18) to more rounds.

Key insight: W[j] controls δW[j+2] through schedule (σ₁ term).
So W[15]→δW[17], W[14]→δW[16]... wait, δW[16] depends on W[14] via σ₁(W[14]).

Schedule: W[i] = W[i-16] + σ₀(W[i-15]) + W[i-7] + σ₁(W[i-2])
So W[i-2] controls W[i] through σ₁.

For δW[r] (r≥16): the σ₁ term uses W[r-2].
Wang chain fixes W[1..15] adaptively. But WHICH words are truly free?

Analysis:
- W[0]: controls Da₁₃ → δe[17]. USED for barrier 1.
- W[15]: neutral for δe[17], controls δW[17] → δe[18]. FREE for barrier 2.
- W[14]: controls δW[16] via σ₁(W[14]). But δW[16] also depends on W[0].
  W[14] is USED by Wang chain (for δe[15]=0). NOT free.
- W[13]: used by Wang for δe[14]=0. NOT free.

So only W[0] and W[15] are "controllable" for barrier optimization.
But W[15] only affects ODD schedule words (17, 19, ...).

NEW IDEA: Instead of 2-phase (W[0] then W[15]),
do JOINT optimization of (W[0], W[15]) minimizing combined HW across rounds 17-20.
"""
import random
import time

M32 = 0xFFFFFFFF

def ror(x, n, bits=32):
    return ((x >> n) | (x << (bits - n))) & ((1 << bits) - 1)

K = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0xfc19dc68, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
]
IV = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]

def compute_de_rounds(W, max_round=20):
    """Compute δe[r] for r=17..max_round with Wang chain, return list of (δe, δa) pairs"""
    W1 = list(W[:16])
    W2 = list(W[:16])
    W2[0] ^= 0x80000000

    a1,b1,c1,d1,e1,f1,g1,h1 = list(IV)
    a2,b2,c2,d2,e2,f2,g2,h2 = list(IV)

    # Wang chain rounds 0-15
    for r in range(16):
        S1_1=(ror(e1,6)^ror(e1,11)^(e1>>25))&M32
        ch1=((e1&f1)^(~e1&g1))&M32
        t1_1=(h1+S1_1+ch1+K[r]+W1[r])&M32
        S0_1=(ror(a1,2)^ror(a1,13)^ror(a1,22))&M32
        maj1=((a1&b1)^(a1&c1)^(b1&c1))&M32
        t2_1=(S0_1+maj1)&M32
        e1_new=(d1+t1_1)&M32; a1_new=(t1_1+t2_1)&M32

        if r >= 1:
            S1_2=(ror(e2,6)^ror(e2,11)^(e2>>25))&M32
            ch2=((e2&f2)^(~e2&g2))&M32
            W2[r]=(e1_new-d2-h2-S1_2-ch2-K[r])&M32

        S1_2=(ror(e2,6)^ror(e2,11)^(e2>>25))&M32
        ch2=((e2&f2)^(~e2&g2))&M32
        t1_2=(h2+S1_2+ch2+K[r]+W2[r])&M32
        S0_2=(ror(a2,2)^ror(a2,13)^ror(a2,22))&M32
        maj2=((a2&b2)^(a2&c2)^(b2&c2))&M32
        t2_2=(S0_2+maj2)&M32
        e2_new=(d2+t1_2)&M32; a2_new=(t1_2+t2_2)&M32

        h1=g1;g1=f1;f1=e1;e1=e1_new;d1=c1;c1=b1;b1=a1;a1=a1_new
        h2=g2;g2=f2;f2=e2;e2=e2_new;d2=c2;c2=b2;b2=a2;a2=a2_new

    # Extend schedules
    for msg in [W1, W2]:
        for i in range(16, max_round):
            if i >= len(msg):
                s0=(ror(msg[i-15],7)^ror(msg[i-15],18)^(msg[i-15]>>3))&M32
                s1=(ror(msg[i-2],17)^ror(msg[i-2],19)^(msg[i-2]>>10))&M32
                msg.append((msg[i-16]+s0+msg[i-7]+s1)&M32)

    # Continue rounds 16..max_round-1, tracking δe and δa
    deltas = []
    for r in range(16, max_round):
        S1_1=(ror(e1,6)^ror(e1,11)^(e1>>25))&M32
        ch1=((e1&f1)^(~e1&g1))&M32
        t1_1=(h1+S1_1+ch1+K[r]+W1[r])&M32
        S0_1=(ror(a1,2)^ror(a1,13)^ror(a1,22))&M32
        maj1=((a1&b1)^(a1&c1)^(b1&c1))&M32
        t2_1=(S0_1+maj1)&M32
        e1_new=(d1+t1_1)&M32; a1_new=(t1_1+t2_1)&M32

        S1_2=(ror(e2,6)^ror(e2,11)^(e2>>25))&M32
        ch2=((e2&f2)^(~e2&g2))&M32
        t1_2=(h2+S1_2+ch2+K[r]+W2[r])&M32
        S0_2=(ror(a2,2)^ror(a2,13)^ror(a2,22))&M32
        maj2=((a2&b2)^(a2&c2)^(b2&c2))&M32
        t2_2=(S0_2+maj2)&M32
        e2_new=(d2+t1_2)&M32; a2_new=(t1_2+t2_2)&M32

        de = (e1_new - e2_new) & M32
        da = (a1_new - a2_new) & M32
        deltas.append((de, da, r+1))

        h1=g1;g1=f1;f1=e1;e1=e1_new;d1=c1;c1=b1;b1=a1;a1=a1_new
        h2=g2;g2=f2;f2=e2;e2=e2_new;d2=c2;c2=b2;b2=a2;a2=a2_new

    return deltas

def hw(x):
    return bin(x).count('1')

random.seed(42)

# ============================================================
# SEARCH 1: Joint (W[0], W[15]) optimization, budget 2M
# ============================================================
print("=" * 70)
print("SEARCH 1: Joint (W[0], W[15]) optimization")
print("Budget: 2M evaluations")
print("=" * 70)

W_rest = [random.randint(0, M32) for _ in range(16)]

best_combined = [256] * 4  # best HW(δe[17..20])
best_total = 256
best_params = None
best_deltas = None

t0 = time.time()
N_BUDGET = 2000000

for attempt in range(N_BUDGET):
    W_rest[0] = random.randint(0, M32)
    W_rest[15] = random.randint(0, M32)

    deltas = compute_de_rounds(W_rest, max_round=21)
    hws = [hw(d[0]) for d in deltas]  # HW of δe for rounds 17-20

    total = sum(hws[:4])  # Sum of first 4 rounds (17-20)

    if total < best_total:
        best_total = total
        best_params = (W_rest[0], W_rest[15])
        best_deltas = [(d[2], hw(d[0]), hw(d[1])) for d in deltas]
        elapsed = time.time() - t0
        print(f"  [{attempt+1:>7d} | {elapsed:5.1f}s] NEW BEST total={total}: " +
              " ".join(f"r{d[2]}:de={hw(d[0]):2d},da={hw(d[1]):2d}" for d in deltas[:5]))

    if attempt % 500000 == 499999:
        print(f"  ... {attempt+1} attempts, best={best_total}, {time.time()-t0:.0f}s")

print(f"\nFINAL: best combined HW(δe[17..20]) = {best_total}")
print(f"  W[0]=0x{best_params[0]:08x}, W[15]=0x{best_params[1]:08x}")
for r, hwe, hwa in best_deltas[:5]:
    print(f"  r={r}: HW(δe)={hwe}, HW(δa)={hwa}")

# ============================================================
# SEARCH 2: Two-phase with larger budget
# ============================================================
print("\n" + "=" * 70)
print("SEARCH 2: Phase 1 = best W[0] for δe[17], Phase 2 = best W[15] for δe[18]")
print("=" * 70)

# Phase 1: Find best W[0] for smallest HW(δe[17])
W_rest2 = [random.randint(0, M32) for _ in range(16)]
best_hw17 = 32
best_w0 = None

t0 = time.time()
for attempt in range(1000000):
    W_rest2[0] = random.randint(0, M32)
    deltas = compute_de_rounds(W_rest2, max_round=17)
    h17 = hw(deltas[0][0])
    if h17 < best_hw17:
        best_hw17 = h17
        best_w0 = W_rest2[0]
        if h17 <= 4:
            print(f"  Phase 1 [{attempt+1}]: HW(δe[17])={h17}, W[0]=0x{best_w0:08x}")

print(f"  Phase 1 result: best HW(δe[17])={best_hw17} ({time.time()-t0:.1f}s)")

# Phase 2: Fix W[0], find best W[15]
W_rest2[0] = best_w0
best_hw18 = 32
best_w15 = None

t0 = time.time()
for attempt in range(1000000):
    W_rest2[15] = random.randint(0, M32)
    deltas = compute_de_rounds(W_rest2, max_round=18)
    h18 = hw(deltas[1][0])  # δe[18]
    if h18 < best_hw18:
        best_hw18 = h18
        best_w15 = W_rest2[15]

print(f"  Phase 2 result: best HW(δe[18])={best_hw18} ({time.time()-t0:.1f}s)")
print(f"  Combined: HW(δe[17])={best_hw17} + HW(δe[18])={best_hw18} = {best_hw17+best_hw18}")

# Phase 2b: also check rounds 19, 20 for this (W[0], W[15])
W_rest2[0] = best_w0
W_rest2[15] = best_w15
deltas_best = compute_de_rounds(W_rest2, max_round=25)
print(f"\n  Full δe profile for best pair:")
for d in deltas_best:
    print(f"    r={d[2]:2d}: HW(δe)={hw(d[0]):2d}, HW(δa)={hw(d[1]):2d}")

# ============================================================
# SEARCH 3: Hill-climbing on (W[0], W[15]) pair
# ============================================================
print("\n" + "=" * 70)
print("SEARCH 3: Hill-climbing on (W[0], W[15])")
print("Minimize combined HW(δe[17]) + HW(δe[18])")
print("=" * 70)

W_hc = [random.randint(0, M32) for _ in range(16)]

# Initial
deltas = compute_de_rounds(W_hc, max_round=18)
current_score = hw(deltas[0][0]) + hw(deltas[1][0])
best_score = current_score
best_hc = (W_hc[0], W_hc[15])

t0 = time.time()
n_improvements = 0

for iteration in range(500000):
    # Randomly flip bits in W[0] or W[15]
    W_trial = list(W_hc)
    target = random.choice([0, 15])
    bit = random.randint(0, 31)
    W_trial[target] ^= (1 << bit)

    deltas = compute_de_rounds(W_trial, max_round=18)
    score = hw(deltas[0][0]) + hw(deltas[1][0])

    if score < current_score:
        current_score = score
        W_hc = W_trial
        n_improvements += 1
        if score < best_score:
            best_score = score
            best_hc = (W_hc[0], W_hc[15])
            d17 = deltas[0][0]
            d18 = deltas[1][0]
            print(f"  [{iteration+1:>6d}] score={score}: HW(δe17)={hw(d17)}, HW(δe18)={hw(d18)}")

    elif score == current_score:
        # Accept with 50% probability (simulated annealing-like)
        if random.random() < 0.3:
            current_score = score
            W_hc = W_trial

    # Restart if stuck
    if iteration % 50000 == 49999:
        if current_score > best_score + 2:
            W_hc[0] = best_hc[0] ^ random.randint(0, M32)
            W_hc[15] = best_hc[1] ^ random.randint(0, M32)
            deltas = compute_de_rounds(W_hc, max_round=18)
            current_score = hw(deltas[0][0]) + hw(deltas[1][0])

print(f"\n  Hill-climbing result: best score={best_score}")
print(f"  ({n_improvements} improvements in 500K iterations, {time.time()-t0:.1f}s)")

# Show full profile for best hill-climbing result
W_hc[0], W_hc[15] = best_hc
deltas_hc = compute_de_rounds(W_hc, max_round=25)
print(f"  Full δe profile:")
for d in deltas_hc:
    print(f"    r={d[2]:2d}: HW(δe)={hw(d[0]):2d}, HW(δa)={hw(d[1]):2d}")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY: Near-collision records")
print("=" * 70)
print(f"  FINAL_REPORT single-word r=17:     HW(ΔH) = 11")
print(f"  Joint (W[0],W[15]) random search:  Σ HW(δe[17..20]) = {best_total}")
print(f"  Two-phase search:                  HW(δe17)+HW(δe18) = {best_hw17}+{best_hw18} = {best_hw17+best_hw18}")
print(f"  Hill-climbing:                     HW(δe17)+HW(δe18) = {best_score}")
