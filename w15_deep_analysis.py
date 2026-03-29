"""
W[15] DEEP ANALYSIS: The only truly neutral word.

Da₁₆ does NOT depend on W[15]. But:
- W[15] enters round 15 (T1[15])
- W[15] enters schedule: W[17] = σ₁(W[15]) + W[10] + σ₀(W[2]) + W[1]
  W[30] depends on W[15] via σ₁
  W[32] depends on W[17] which depends on W[15]

Questions:
1. What does W[15] control in the schedule?
2. Can W[15] freedom be used to control δW[16+]?
3. Does W[15] affect δe[17] through SCHEDULE (not through Da)?
4. Multi-target: does W[15] create multiple viable δW[16] values?
"""
import random

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

def full_schedule(W):
    """Compute full 64-word schedule from 16 input words"""
    w = list(W[:16])
    for i in range(16, 64):
        s0 = (ror(w[i-15], 7) ^ ror(w[i-15], 18) ^ (w[i-15] >> 3)) & M32
        s1 = (ror(w[i-2], 17) ^ ror(w[i-2], 19) ^ (w[i-2] >> 10)) & M32
        w.append((w[i-16] + s0 + w[i-7] + s1) & M32)
    return w

def wang_chain_full(W, dW0=0x80000000):
    """Wang chain: returns (W1, W2, states1, states2) for full analysis"""
    W1 = list(W)
    W2 = list(W)
    W2[0] ^= dW0

    a1,b1,c1,d1,e1,f1,g1,h1 = list(IV)
    a2,b2,c2,d2,e2,f2,g2,h2 = list(IV)

    states1 = [(a1,b1,c1,d1,e1,f1,g1,h1)]
    states2 = [(a2,b2,c2,d2,e2,f2,g2,h2)]

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
            needed=(e1_new-d2-h2-S1_2-ch2-K[r])&M32
            W2[r]=needed

        S1_2=(ror(e2,6)^ror(e2,11)^(e2>>25))&M32
        ch2=((e2&f2)^(~e2&g2))&M32
        t1_2=(h2+S1_2+ch2+K[r]+W2[r])&M32
        S0_2=(ror(a2,2)^ror(a2,13)^ror(a2,22))&M32
        maj2=((a2&b2)^(a2&c2)^(b2&c2))&M32
        t2_2=(S0_2+maj2)&M32
        e2_new=(d2+t1_2)&M32; a2_new=(t1_2+t2_2)&M32

        h1=g1;g1=f1;f1=e1;e1=e1_new;d1=c1;c1=b1;b1=a1;a1=a1_new
        h2=g2;g2=f2;f2=e2;e2=e2_new;d2=c2;c2=b2;b2=a2;a2=a2_new

        states1.append((a1,b1,c1,d1,e1,f1,g1,h1))
        states2.append((a2,b2,c2,d2,e2,f2,g2,h2))

    return W1, W2, states1, states2

random.seed(42)

# ============================================================
# EXP 1: How W[15] affects the schedule difference δW[16..63]
# ============================================================
print("=" * 70)
print("EXP 1: Schedule difference δW[r] as function of W[15]")
print("=" * 70)

W_base = [random.randint(0, M32) for _ in range(16)]
W1_base, W2_base, _, _ = wang_chain_full(W_base)

sched1_base = full_schedule(W1_base)
sched2_base = full_schedule(W2_base)
dW_base = [(s1 - s2) & M32 for s1, s2 in zip(sched1_base, sched2_base)]

# Now vary W[15] and see how δW[16..] changes
print("\n  Varying W[15]: which schedule words change?")
n_test = 1000
dW_changed = [0] * 64

for _ in range(n_test):
    W_test = list(W_base)
    W_test[15] = random.randint(0, M32)
    W1_t, W2_t, _, _ = wang_chain_full(W_test)
    sched1_t = full_schedule(W1_t)
    sched2_t = full_schedule(W2_t)
    dW_t = [(s1 - s2) & M32 for s1, s2 in zip(sched1_t, sched2_t)]

    for r in range(64):
        if dW_t[r] != dW_base[r]:
            dW_changed[r] += 1

print(f"  {'Round':>5} | {'Changed %':>10} | {'Note':>20}")
print("  " + "-" * 45)
for r in range(64):
    pct = 100 * dW_changed[r] / n_test
    note = ""
    if r < 16 and r != 15:
        note = "Wang-fixed"
    elif r == 15:
        note = "W[15] direct"
    elif pct > 99:
        note = "FULLY depends"
    elif pct < 1:
        note = "INDEPENDENT"
    if r <= 20 or r % 8 == 0 or pct < 1 or (r < 25 and pct < 100):
        print(f"  {r:>5} | {pct:>9.1f}% | {note:>20}")

# ============================================================
# EXP 2: Does W[15] affect δe[17]?
# ============================================================
print("\n" + "=" * 70)
print("EXP 2: Does W[15] affect δe[17]?")
print("=" * 70)

# δe[17] = Da₁₆ + δW[16]
# Da₁₆ doesn't depend on W[15] (confirmed in EXP 3 of parallel_experiments)
# δW[16] = schedule(W1)[16] - schedule(W2)[16]
# W[16] = W[0] + σ₀(W[1]) + W[9] + σ₁(W[14])
# W[15] NOT in this formula!

# But W2[r] depends on Wang corrections which depend on state which...
# Actually: Wang corrections for r=1..14 don't involve W[15]
# Wang correction for r=15: W2[15] chosen to zero δe[16]

# W2[15] = needed_val. When we change W1[15], W2[15] changes (to maintain δe[16]=0)
# But δW[15] = W1[15] - W2[15] changes. Does δW[16] depend on δW[15]?

# W[16] = W[0] + σ₀(W[1]) + W[9] + σ₁(W[14])
# δW[16] = δW[0] + σ₀(δW[1]) + δW[9] + σ₁(δW[14])
# W[15] NOT in the formula for W[16]!
# So δW[16] does NOT depend on W[15] NOR on δW[15].

# Verify:
dW16_values = set()
de17_values = set()
n_test = 10000
for _ in range(n_test):
    W_test = list(W_base)
    W_test[15] = random.randint(0, M32)
    W1_t, W2_t, s1, s2 = wang_chain_full(W_test)
    sched1 = full_schedule(W1_t)
    sched2 = full_schedule(W2_t)
    dW16 = (sched1[16] - sched2[16]) & M32
    dW16_values.add(dW16)

    # Compute δe[17]
    da16 = (s1[16][0] - s2[16][0]) & M32  # Da at round 16
    de17 = (da16 + dW16) & M32  # Simplified
    # Actually: δe[17] = δd[16] + δT1[16]
    # δd[16] = δa[13]
    # δT1[16] depends on δe[16]=0, δh[16]=δe[13]=0, δW[16]
    # So δT1[16] = δW[16] (since all state diffs in e-branch = 0)
    # δe[17] = δa[13] + δW[16] = Da₁₃ + δW[16]
    de17_values.add(de17)

print(f"  Unique δW[16] values: {len(dW16_values)} (out of {n_test})")
print(f"  Unique δe[17] values: {len(de17_values)} (out of {n_test})")
print(f"  δW[16] = constant when varying W[15]: {'YES' if len(dW16_values)==1 else 'NO'}")
print(f"  δe[17] = constant when varying W[15]: {'YES' if len(de17_values)==1 else 'NO'}")

if len(dW16_values) == 1:
    print(f"  δW[16] = 0x{list(dW16_values)[0]:08x}")
if len(de17_values) == 1:
    print(f"  δe[17] = 0x{list(de17_values)[0]:08x}")

# ============================================================
# EXP 3: W[15] controls δW[17] — can we exploit?
# ============================================================
print("\n" + "=" * 70)
print("EXP 3: W[15] control over δW[17] and beyond")
print("=" * 70)

# W[17] = W[1] + σ₀(W[2]) + W[10] + σ₁(W[15])
# δW[17] = δW[1] + σ₀(δW[2]) + δW[10] + σ₁(δW[15])
# W[15] enters through σ₁(δW[15])!

# When we vary W1[15] (with Wang correction adjusting W2[15]):
# δW[15] = W1[15] - W2[15] changes
# → δW[17] changes via σ₁(δW[15])

# How much freedom do we have in δW[17]?
dW17_values = set()
dW15_values = set()
for _ in range(n_test):
    W_test = list(W_base)
    W_test[15] = random.randint(0, M32)
    W1_t, W2_t, _, _ = wang_chain_full(W_test)
    dW15 = (W1_t[15] - W2_t[15]) & M32
    dW15_values.add(dW15)

    sched1 = full_schedule(W1_t)
    sched2 = full_schedule(W2_t)
    dW17 = (sched1[17] - sched2[17]) & M32
    dW17_values.add(dW17)

print(f"  Unique δW[15] values: {len(dW15_values)} (should ≈ {n_test})")
print(f"  Unique δW[17] values: {len(dW17_values)}")
print(f"  δW[17] coverage: {100*len(dW17_values)/n_test:.1f}%")

# ============================================================
# EXP 4: Can W[15] zero δe[18]?
# ============================================================
print("\n" + "=" * 70)
print("EXP 4: Can W[15] control δe[18]?")
print("=" * 70)

# δe[18] = δd[17] + δT1[17]
# δd[17] = δa[14]  (from Wang chain, determined by W[0..14])
# δT1[17] depends on δW[17] (which W[15] controls!) and δstate[17]
# δstate[17] = (δa[17], δa[16], δa[15], δa[14], δe[17], 0, 0, 0)
# If δe[17] ≠ 0: δstate[17] has both a-branch and e-branch ≠ 0

# We want: δe[18] = 0
# This requires: δd[17] + δT1[17] = 0
# δd[17] = δa[14] (fixed, doesn't depend on W[15])
# δT1[17] = δh[17] + Σ₁(δe[17]) + δCh[17] + δW[17]
# = δe[14] + Σ₁(δe[17]) + Ch_diff + δW[17]
# δe[14] = 0 (Wang chain), so δh[17] = δe[14] = 0
# δCh[17] depends on δe[17], δf[17]=δe[16]=0, δg[17]=δe[15]=0
# δCh = δe[17] & (f[17] ^ g[17]) (from T_CH_LINEAR)

# So: δT1[17] = Σ₁(δe[17]) + δe[17]·(f[17]⊕g[17]) + δW[17]
# And: δe[18] = δa[14] + Σ₁(δe[17]) + δe[17]·mask + δW[17] = 0
# We control δW[17] through W[15]. Can we solve?
# δW[17] = -δa[14] - Σ₁(δe[17]) - δe[17]·mask

# The question: can W[15] give us the NEEDED δW[17]?

# First: what δW[17] do we need?
n_test_small = 1000
n_found = 0
min_de18 = M32

for trial in range(n_test_small):
    W_test = list(W_base)
    W_test[0] = random.randint(0, M32)

    # Try all W[15] values (subsample)
    best_de18 = M32
    for w15_trial in range(256):
        W_test[15] = random.randint(0, M32)
        W1_t, W2_t, s1, s2 = wang_chain_full(W_test)

        # Compute δe[17] and δe[18] by running one more round
        sched1 = full_schedule(W1_t)
        sched2 = full_schedule(W2_t)

        # State at round 16
        a1,b1,c1,d1,e1,f1,g1,h1 = s1[16]
        a2,b2,c2,d2,e2,f2,g2,h2 = s2[16]

        # Round 16
        S1_1=(ror(e1,6)^ror(e1,11)^(e1>>25))&M32
        ch1=((e1&f1)^(~e1&g1))&M32
        t1_1=(h1+S1_1+ch1+K[16]+sched1[16])&M32
        S0_1=(ror(a1,2)^ror(a1,13)^ror(a1,22))&M32
        maj1=((a1&b1)^(a1&c1)^(b1&c1))&M32
        t2_1=(S0_1+maj1)&M32

        S1_2=(ror(e2,6)^ror(e2,11)^(e2>>25))&M32
        ch2=((e2&f2)^(~e2&g2))&M32
        t1_2=(h2+S1_2+ch2+K[16]+sched2[16])&M32
        S0_2=(ror(a2,2)^ror(a2,13)^ror(a2,22))&M32
        maj2=((a2&b2)^(a2&c2)^(b2&c2))&M32
        t2_2=(S0_2+maj2)&M32

        e1_17=(d1+t1_1)&M32; a1_17=(t1_1+t2_1)&M32
        e2_17=(d2+t1_2)&M32; a2_17=(t1_2+t2_2)&M32

        f1_17=e1; g1_17=f1; h1_17=g1; b1_17=a1; c1_17=b1; d1_17=c1
        f2_17=e2; g2_17=f2; h2_17=g2; b2_17=a2; c2_17=b2; d2_17=c2

        de17 = (e1_17 - e2_17) & M32

        # Round 17
        S1_1r=(ror(e1_17,6)^ror(e1_17,11)^(e1_17>>25))&M32
        ch1r=((e1_17&f1_17)^(~e1_17&g1_17))&M32
        t1_1r=(h1_17+S1_1r+ch1r+K[17]+sched1[17])&M32

        S1_2r=(ror(e2_17,6)^ror(e2_17,11)^(e2_17>>25))&M32
        ch2r=((e2_17&f2_17)^(~e2_17&g2_17))&M32
        t1_2r=(h2_17+S1_2r+ch2r+K[17]+sched2[17])&M32

        e1_18=(d1_17+t1_1r)&M32
        e2_18=(d2_17+t1_2r)&M32

        de18 = (e1_18 - e2_18) & M32
        hw_de18 = bin(de18).count('1')

        if hw_de18 < best_de18:
            best_de18 = hw_de18

    if best_de18 == 0:
        n_found += 1
    min_de18 = min(min_de18, best_de18)

    if trial < 10:
        print(f"  W[0]=trial {trial}: best HW(δe[18]) = {best_de18} (from 256 W[15] trials)")

print(f"\n  Summary ({n_test_small} random W[0], 256 W[15] each):")
print(f"  δe[18]=0 found: {n_found}/{n_test_small}")
print(f"  Minimum HW(δe[18]): {min_de18}")
print(f"  Expected by random (256 trials, 32-bit): P(HW≤{min_de18}) ≈ ...")

# ============================================================
# EXP 5: Birthday on δe[17] using W[0] freedom
# ============================================================
print("\n" + "=" * 70)
print("EXP 5: Birthday search for δe[17]=0")
print("=" * 70)

# For each W[0], compute δe[17]. Find W[0] with δe[17] = 0.
# This is single-target preimage on 32 bits → expected 2^32 trials.
# But we have W[15] as ADDITIONAL freedom. Can we do 2D search?

# δe[17] = Da₁₃(W[0]) + δW[16](W[0..14])
# W[15] doesn't affect δe[17]. So W[15] freedom is useless for δe[17]=0.

# But W[15] CAN affect δe[18]:
# δe[18] depends on δW[17] which depends on W[15].
# So: find W[0] with δe[17]=0 (cost 2^32), THEN use W[15] for δe[18]=0.

# If W[15] gives enough freedom to zero δe[18]:
# Total cost = 2^32 (for δe[17]) + 2^32 (for δe[18] via W[15])
# = 2^33 for TWO barriers instead of 2^64!

# Check: does W[15] have 2^32 distinct effects on δe[18]?

print("  Testing δW[17] range from W[15]...")
W_fixed = [random.randint(0, M32) for _ in range(16)]
dW17_set = set()
for i in range(min(100000, 1 << 17)):
    W_fixed[15] = random.randint(0, M32)
    W1, W2, _, _ = wang_chain_full(W_fixed)
    s1 = full_schedule(W1)
    s2 = full_schedule(W2)
    dW17 = (s1[17] - s2[17]) & M32
    dW17_set.add(dW17)

print(f"  Unique δW[17] from {min(100000, 1<<17)} W[15] trials: {len(dW17_set)}")
print(f"  Coverage: {100*len(dW17_set)/min(100000, 1<<17):.1f}%")

if len(dW17_set) == min(100000, 1 << 17):
    print(f"  δW[17] appears to be a BIJECTION of W[15] — full 2^32 range!")
    print(f"  *** W[15] can independently target δW[17] ***")
    print(f"  *** Two barriers (δe[17]=0, δe[18]=0) for cost 2^32 + 2^32 = 2^33 ***")
else:
    print(f"  δW[17] NOT bijective — limited range")
