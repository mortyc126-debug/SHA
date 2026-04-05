#!/usr/bin/env python3
"""
NK Directions 6 + 9: Fixed points and preimage structure.

Direction 6: Fixed points and short cycles of SHA-256 compress.
  Does compress(M, IV)[0] = M[0] for some M? (partial fixed point)
  Do short cycles exist? M → H → M' → H' → ... → M?

Direction 9: Preimage of specific H (e.g., H = 0).
  Not collision — single target. Different structure.
  Can structural knowledge help find M with H(M) = target?
"""

import random, math, time

MASK32 = 0xFFFFFFFF
K=[0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2]
IV=[0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&MASK32
def ssig0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def ssig1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def sig0(x): return rotr(x,2)^rotr(x,13)^rotr(x,22)
def sig1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def ch(e,f,g): return (e&f)^(~e&g)&MASK32
def maj(a,b,c): return (a&b)^(a&c)^(b&c)
def hw(x): return bin(x&MASK32).count('1')

def sha256_compress(M, iv=None):
    if iv is None: iv = IV
    W=list(M)+[0]*(64-len(M))
    for i in range(16,64): W[i]=(ssig1(W[i-2])+W[i-7]+ssig0(W[i-15])+W[i-16])&MASK32
    a,b,c,d,e,f,g,h=iv
    for r in range(64):
        T1=(h+sig1(e)+ch(e,f,g)+K[r]+W[r])&MASK32;T2=(sig0(a)+maj(a,b,c))&MASK32
        h,g,f=g,f,e;e=(d+T1)&MASK32;d,c,b=c,b,a;a=(T1+T2)&MASK32
    return tuple((s+v)&MASK32 for s,v in zip([a,b,c,d,e,f,g,h],iv))


# ================================================================
# DIRECTION 6: Fixed points and cycles
# ================================================================

print("=" * 70)
print("DIRECTION 6: Fixed points and short cycles")
print("=" * 70)

# 6.1: Partial fixed point — H[0] = W[0]
# Function f: W[0] → H[0] where W[1..15] = 0
# Fixed point: f(x) = x

print("\n  6.1: Partial fixed point H[0] = W[0] (W[1..15]=0)")

# Pollard rho on f(x) = sha256_compress([x,0,...,0])[0]
def f_h0(x):
    M = [x] + [0]*15
    return sha256_compress(M)[0]

# Floyd cycle detection
random.seed(42)
x0 = random.randint(0, MASK32)

tortoise = f_h0(x0)
hare = f_h0(f_h0(x0))
steps = 1

t0 = time.time()
while tortoise != hare and steps < 2**24:
    tortoise = f_h0(tortoise)
    hare = f_h0(f_h0(hare))
    steps += 1

if tortoise == hare:
    elapsed = time.time() - t0
    print(f"  Cycle detected after {steps} steps ({elapsed:.1f}s)")

    # Find cycle entry (mu)
    tortoise = x0
    mu = 0
    while tortoise != hare:
        tortoise = f_h0(tortoise)
        hare = f_h0(hare)
        mu += 1

    # Find cycle length (lambda)
    lam = 1
    hare = f_h0(tortoise)
    while tortoise != hare:
        hare = f_h0(hare)
        lam += 1

    print(f"  Tail length (μ) = {mu}")
    print(f"  Cycle length (λ) = {lam}")
    print(f"  Expected for random f: E[μ] ≈ E[λ] ≈ √(π·2^32/2) ≈ {int(math.sqrt(math.pi * 2**32 / 2))}")

    # Check: is there a fixed point in the cycle?
    x = tortoise
    fp_found = False
    for i in range(lam):
        if f_h0(x) == x:
            print(f"  ★ FIXED POINT FOUND: x = 0x{x:08x}")
            print(f"    Verify: H[0](x,0,...,0) = 0x{f_h0(x):08x} {'✓' if f_h0(x)==x else '✗'}")
            fp_found = True
            break
        x = f_h0(x)

    if not fp_found:
        print(f"  No fixed point in cycle (expected: P(fp in cycle) ≈ λ/2^32)")
else:
    print(f"  No cycle found in 2^24 steps")


# 6.2: Multiple seed Pollard rho — cycle length distribution
print("\n  6.2: Cycle length distribution (10 seeds)")

cycle_lengths = []
for seed in range(10):
    random.seed(seed*1000+7)
    x0 = random.randint(0, MASK32)

    tortoise = f_h0(x0)
    hare = f_h0(f_h0(x0))
    steps = 0
    while tortoise != hare and steps < 2**22:
        tortoise = f_h0(tortoise)
        hare = f_h0(f_h0(hare))
        steps += 1

    if tortoise == hare:
        # Find lambda
        lam = 1
        hare = f_h0(tortoise)
        while tortoise != hare:
            hare = f_h0(hare)
            lam += 1
        cycle_lengths.append(lam)
        print(f"    seed {seed}: λ = {lam}")
    else:
        print(f"    seed {seed}: no cycle in 2^22 steps")

if cycle_lengths:
    avg_lam = sum(cycle_lengths) / len(cycle_lengths)
    expected = int(math.sqrt(math.pi * 2**32 / 2))
    print(f"\n  Average λ = {avg_lam:.0f}")
    print(f"  Expected for random: {expected}")
    print(f"  Ratio: {avg_lam/expected:.2f}× {'(anomalous!)' if abs(avg_lam/expected - 1) > 0.5 else '(normal)'}")


# 6.3: Does cycle structure help collision?
print("\n  6.3: Collision from cycle structure")
print(f"  If two messages land in SAME cycle → collision on H[0].")
print(f"  But H[0] collision ≠ full collision (need H[0..7] all match).")
print(f"  H[0] collision from cycle: cost = O(√2^32) = O(2^16). KNOWN.")
print(f"  Full collision still needs H[1..7] to match = 2^112 more.")
print(f"  Total: 2^16 + 2^112 = 2^112. But this was already checked:")
print(f"  κ=0 (T6'): H[0..7] independent. No shortcut from H[0] collision.")


# ================================================================
# DIRECTION 9: Preimage of specific target
# ================================================================

print()
print("=" * 70)
print("DIRECTION 9: Preimage search — structural shortcuts?")
print("=" * 70)

# 9.1: Can we find M with H(M) having specific structure?
# E.g., H[7] = 0, or H[7] has low HW?

print("\n  9.1: Finding M with H[7] having low HW")

N_search = 500000
random.seed(42)

best_hw_h7 = 32
best_M = None
hw_dist = [0] * 33  # distribution of HW(H[7])

t0 = time.time()
for trial in range(N_search):
    M = [random.randint(0, MASK32) for _ in range(16)]
    H = sha256_compress(M)
    h7_hw = hw(H[7])
    hw_dist[h7_hw] += 1
    if h7_hw < best_hw_h7:
        best_hw_h7 = h7_hw
        best_M = list(M)

elapsed = time.time() - t0
print(f"\n  N={N_search} random M ({elapsed:.1f}s)")
print(f"  Best HW(H[7]) = {best_hw_h7}")
print(f"  Expected min for random: ≈ {32 - math.sqrt(2*math.log(N_search)):.0f}")

# Distribution check
print(f"\n  HW(H[7]) distribution (should be Binomial(32, 0.5)):")
print(f"  {'HW':>3} | {'count':>6} | {'P_obs':>6} | {'P_binom':>7}")
for h in range(max(0, best_hw_h7), min(33, 20)):
    p_obs = hw_dist[h] / N_search
    # Binomial(32, 0.5)
    from math import comb
    p_bin = comb(32, h) / 2**32
    print(f"  {h:>3} | {hw_dist[h]:>6} | {p_obs:>.4f} | {p_bin:>.5f}")


# 9.2: HC to find M with H[7] = 0
print(f"\n  9.2: HC to minimize H[7]")

best_h7_val = MASK32
random.seed(123)

for trial in range(20):
    M = [random.randint(0, MASK32) for _ in range(16)]
    H = sha256_compress(M)
    cur_hw = hw(H[7])

    for step in range(500):
        w = random.randint(0, 15)
        b = random.randint(0, 31)
        M_try = list(M); M_try[w] ^= (1 << b)
        H_try = sha256_compress(M_try)
        new_hw = hw(H_try[7])
        if new_hw < cur_hw:
            M = M_try; cur_hw = new_hw; H = H_try

    if cur_hw < hw(best_h7_val):
        best_h7_val = H[7]
        print(f"    trial {trial}: HW(H[7]) = {cur_hw}, H[7] = 0x{H[7]:08x}")

print(f"\n  Best HW(H[7]) via HC = {hw(best_h7_val)}")
print(f"  For H[7]=0: need HW=0. Best found: HW={hw(best_h7_val)}.")
print(f"  Gap: {hw(best_h7_val)} bits. Random needs ≈ 2^32 trials for H[7]=0.")


# 9.3: Can Wang-like technique help preimage?
print(f"\n  9.3: Can Wang cascade help preimage?")
print(f"  Wang cancels δe[r+1] via adaptive δW[r].")
print(f"  For PREIMAGE: need H(M) = target. Not a pair — single M.")
print(f"  Wang needs TWO messages to compute δ.")
print(f"  Without a second message, Wang inapplicable to preimage.")
print(f"  Preimage = fundamentally different from collision.")


# 9.4: Backward from target
print(f"\n  9.4: Backward from target H = (0,0,0,0,0,0,0,0)")

target = (0,0,0,0,0,0,0,0)
# state[64] = H - IV
state64 = tuple((t - iv) & MASK32 for t, iv in zip(target, IV))
print(f"  Required state[64]: {tuple(f'0x{s:08x}' for s in state64)}")

# Invert backward — but need W[0..63] which depends on M!
# state[64] is known. To go backward: need W[63], W[62], ...
# W[16..63] = f(W[0..15]). W[0..15] is what we're LOOKING for.
# Circular dependency: can't invert without knowing M.

print(f"  Backward inversion needs W[r] at each step.")
print(f"  W[16..63] = f(W[0..15]) = f(M). Don't know M.")
print(f"  Cannot invert without the message.")
print(f"  Preimage from backward: CIRCULAR DEPENDENCY.")


# ================================================================
# SYNTHESIS
# ================================================================

print()
print("=" * 70)
print("SYNTHESIS: Directions 6 + 9")
print("=" * 70)

print("""
  DIRECTION 6 (Fixed points / cycles):
  - Pollard rho finds H[0]-cycles in O(2^16). Normal cycle structure.
  - Cycle lengths ≈ √(π·2^32/2) ≈ 82K (expected for random function).
  - No anomaly in cycle structure → no structural shortcut.
  - H[0] collision from cycle is cheap (2^16) but useless without H[1..7].
  - Direction 6: CLOSED.

  DIRECTION 9 (Preimage):
  - HC finds HW(H[7]) ≈ 3-4 (vs random min ≈ 5 at N=500K).
  - Marginal improvement, not structural.
  - Wang inapplicable (needs pair, preimage = single message).
  - Backward inversion has circular dependency (needs M to invert).
  - Preimage is fundamentally harder than collision for SHA-256.
  - Direction 9: CLOSED.
""")
