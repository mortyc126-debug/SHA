"""
TWO-BARRIER ATTACK VERIFICATION

Theory: W[0] controls δe[17], W[15] controls δe[18] INDEPENDENTLY.
Cost: 2^32 + 2^32 = 2^33 for TWO barriers (vs 2^64 standard).

This experiment VERIFIES the independence and tests the attack on small scale.
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

def compute_de17_de18(W):
    """Given 16-word message, compute Wang chain and return (δe[17], δe[18], W2)"""
    W1 = list(W)
    W2 = list(W)
    W2[0] ^= 0x80000000

    a1,b1,c1,d1,e1,f1,g1,h1 = list(IV)
    a2,b2,c2,d2,e2,f2,g2,h2 = list(IV)

    for r in range(18):
        # Schedule
        if r >= 16:
            for rr in range(len(W1), r+1):
                s0 = (ror(W1[rr-15],7)^ror(W1[rr-15],18)^(W1[rr-15]>>3))&M32
                s1 = (ror(W1[rr-2],17)^ror(W1[rr-2],19)^(W1[rr-2]>>10))&M32
                W1.append((W1[rr-16]+s0+W1[rr-7]+s1)&M32)
            for rr in range(len(W2), r+1):
                s0 = (ror(W2[rr-15],7)^ror(W2[rr-15],18)^(W2[rr-15]>>3))&M32
                s1 = (ror(W2[rr-2],17)^ror(W2[rr-2],19)^(W2[rr-2]>>10))&M32
                W2.append((W2[rr-16]+s0+W2[rr-7]+s1)&M32)

        # Round r for msg1
        S1_1=(ror(e1,6)^ror(e1,11)^(e1>>25))&M32
        ch1=((e1&f1)^(~e1&g1))&M32
        t1_1=(h1+S1_1+ch1+K[r]+W1[r])&M32
        S0_1=(ror(a1,2)^ror(a1,13)^ror(a1,22))&M32
        maj1=((a1&b1)^(a1&c1)^(b1&c1))&M32
        t2_1=(S0_1+maj1)&M32
        e1_new=(d1+t1_1)&M32; a1_new=(t1_1+t2_1)&M32

        # Wang correction for msg2 (rounds 1..15)
        if 1 <= r <= 15:
            S1_2=(ror(e2,6)^ror(e2,11)^(e2>>25))&M32
            ch2=((e2&f2)^(~e2&g2))&M32
            needed=(e1_new-d2-h2-S1_2-ch2-K[r])&M32
            W2[r]=needed

        # Round r for msg2
        S1_2=(ror(e2,6)^ror(e2,11)^(e2>>25))&M32
        ch2=((e2&f2)^(~e2&g2))&M32
        t1_2=(h2+S1_2+ch2+K[r]+W2[r])&M32
        S0_2=(ror(a2,2)^ror(a2,13)^ror(a2,22))&M32
        maj2=((a2&b2)^(a2&c2)^(b2&c2))&M32
        t2_2=(S0_2+maj2)&M32
        e2_new=(d2+t1_2)&M32; a2_new=(t1_2+t2_2)&M32

        h1=g1;g1=f1;f1=e1;e1=e1_new;d1=c1;c1=b1;b1=a1;a1=a1_new
        h2=g2;g2=f2;f2=e2;e2=e2_new;d2=c2;c2=b2;b2=a2;a2=a2_new

    de17 = (e1_new - e2_new) & M32  # Actually this is de[18]...

    # We need to track de[17] and de[18] separately
    # Let me redo with explicit tracking
    return None

def compute_deltas(W):
    """Compute δe[r] for r=17,18 with Wang chain on W[0..15]"""
    W1 = list(W[:16])
    W2 = list(W[:16])
    W2[0] ^= 0x80000000

    states1 = []; states2 = []
    a1,b1,c1,d1,e1,f1,g1,h1 = list(IV)
    a2,b2,c2,d2,e2,f2,g2,h2 = list(IV)

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

    # Now rounds 16 and 17 (using schedule)
    # Extend schedule
    for msg in [W1, W2]:
        for i in range(16, 18):
            s0=(ror(msg[i-15],7)^ror(msg[i-15],18)^(msg[i-15]>>3))&M32
            s1=(ror(msg[i-2],17)^ror(msg[i-2],19)^(msg[i-2]>>10))&M32
            msg.append((msg[i-16]+s0+msg[i-7]+s1)&M32)

    # Round 16
    for label, (a,b,c,d,e,f,g,h,ww) in [
        ("1", (a1,b1,c1,d1,e1,f1,g1,h1,W1)),
        ("2", (a2,b2,c2,d2,e2,f2,g2,h2,W2))]:

        S1_r=(ror(e,6)^ror(e,11)^(e>>25))&M32
        ch_r=((e&f)^(~e&g))&M32
        t1_r=(h+S1_r+ch_r+K[16]+ww[16])&M32
        S0_r=(ror(a,2)^ror(a,13)^ror(a,22))&M32
        maj_r=((a&b)^(a&c)^(b&c))&M32
        t2_r=(S0_r+maj_r)&M32
        e_new=(d+t1_r)&M32; a_new=(t1_r+t2_r)&M32

        if label == "1":
            h1_16=g;g1_16=f;f1_16=e;e1_16=e_new;d1_16=c;c1_16=b;b1_16=a;a1_16=a_new
        else:
            h2_16=g;g2_16=f;f2_16=e;e2_16=e_new;d2_16=c;c2_16=b;b2_16=a;a2_16=a_new

    de17 = (e1_16 - e2_16) & M32

    # Round 17
    for label, (a,b,c,d,e,f,g,h,ww) in [
        ("1", (a1_16,b1_16,c1_16,d1_16,e1_16,f1_16,g1_16,h1_16,W1)),
        ("2", (a2_16,b2_16,c2_16,d2_16,e2_16,f2_16,g2_16,h2_16,W2))]:

        S1_r=(ror(e,6)^ror(e,11)^(e>>25))&M32
        ch_r=((e&f)^(~e&g))&M32
        t1_r=(h+S1_r+ch_r+K[17]+ww[17])&M32
        S0_r=(ror(a,2)^ror(a,13)^ror(a,22))&M32
        maj_r=((a&b)^(a&c)^(b&c))&M32
        t2_r=(S0_r+maj_r)&M32
        e_new=(d+t1_r)&M32

        if label == "1":
            e1_17 = e_new
        else:
            e2_17 = e_new

    de18 = (e1_17 - e2_17) & M32

    return de17, de18

random.seed(42)

# ============================================================
# TEST 1: Verify independence of δe[17] from W[15]
# ============================================================
print("=" * 70)
print("TEST 1: Independence verification")
print("=" * 70)

W_base = [random.randint(0, M32) for _ in range(16)]

de17_set = set()
de18_set = set()
for _ in range(1000):
    W_test = list(W_base)
    W_test[15] = random.randint(0, M32)
    d17, d18 = compute_deltas(W_test)
    de17_set.add(d17)
    de18_set.add(d18)

print(f"  Varying W[15] (W[0] fixed):")
print(f"    Unique δe[17]: {len(de17_set)} (should be 1 if independent)")
print(f"    Unique δe[18]: {len(de18_set)} (should be many if W[15] controls)")

# ============================================================
# TEST 2: Verify independence of δe[18] from W[0] GIVEN δe[17]=fixed
# ============================================================
print("\n" + "=" * 70)
print("TEST 2: Is δe[18] decomposable as f(W[0]) + g(W[15])?")
print("=" * 70)

# If δe[18] = f(W[0]) + g(W[15]) mod 2^32:
# Then δe[18](W0_a, W15_a) + δe[18](W0_b, W15_b)
#    = δe[18](W0_a, W15_b) + δe[18](W0_b, W15_a)  (additivity test)

n_pass = 0
n_test = 5000
for _ in range(n_test):
    W0_a = random.randint(0, M32)
    W0_b = random.randint(0, M32)
    W15_a = random.randint(0, M32)
    W15_b = random.randint(0, M32)

    W_rest = [random.randint(0, M32) for _ in range(16)]

    W_aa = list(W_rest); W_aa[0] = W0_a; W_aa[15] = W15_a
    W_bb = list(W_rest); W_bb[0] = W0_b; W_bb[15] = W15_b
    W_ab = list(W_rest); W_ab[0] = W0_a; W_ab[15] = W15_b
    W_ba = list(W_rest); W_ba[0] = W0_b; W_ba[15] = W15_a

    _, d18_aa = compute_deltas(W_aa)
    _, d18_bb = compute_deltas(W_bb)
    _, d18_ab = compute_deltas(W_ab)
    _, d18_ba = compute_deltas(W_ba)

    # Additivity: d18(a,a) + d18(b,b) = d18(a,b) + d18(b,a) mod 2^32
    lhs = (d18_aa + d18_bb) & M32
    rhs = (d18_ab + d18_ba) & M32

    if lhs == rhs:
        n_pass += 1

print(f"  Additivity test: {n_pass}/{n_test} pass ({100*n_pass/n_test:.1f}%)")
print(f"  (100% = separable, 0% = fully entangled)")

# ============================================================
# TEST 3: Small-scale two-barrier search
# ============================================================
print("\n" + "=" * 70)
print("TEST 3: Two-barrier search (small scale)")
print("=" * 70)

# Phase 1: Find W[0] with δe[17] = 0
# Phase 2: Find W[15] with δe[18] = 0

W_rest = [random.randint(0, M32) for _ in range(16)]

print("  Phase 1: Searching for δe[17]=0 by varying W[0]...")
t0 = time.time()
found_w0 = None
for attempt in range(1 << 20):  # Try up to 1M
    W_rest[0] = random.randint(0, M32)
    d17, _ = compute_deltas(W_rest)
    if d17 == 0:
        found_w0 = W_rest[0]
        print(f"  FOUND! W[0]=0x{found_w0:08x} after {attempt+1} attempts ({time.time()-t0:.1f}s)")
        break
    if attempt % 100000 == 99999:
        print(f"    {attempt+1} attempts... (expected: ~2^32)")

if found_w0 is None:
    print(f"  Not found in 1M attempts (expected: need ~2^32)")
    print(f"  Testing two-phase with PARTIAL match instead...")

    # Find W[0] with smallest HW(δe[17])
    best_hw17 = 32
    best_w0 = None
    for attempt in range(100000):
        W_rest[0] = random.randint(0, M32)
        d17, _ = compute_deltas(W_rest)
        h = bin(d17).count('1')
        if h < best_hw17:
            best_hw17 = h
            best_w0 = W_rest[0]

    print(f"  Best partial: HW(δe[17])={best_hw17} with W[0]=0x{best_w0:08x}")
    W_rest[0] = best_w0

    # Phase 2: Fix this W[0], find W[15] with smallest HW(δe[18])
    print(f"  Phase 2: Fix W[0]=0x{best_w0:08x}, vary W[15]...")
    best_hw18 = 32
    for attempt in range(100000):
        W_rest[15] = random.randint(0, M32)
        d17, d18 = compute_deltas(W_rest)
        h17 = bin(d17).count('1')
        h18 = bin(d18).count('1')
        total = h17 + h18
        if h17 + h18 < best_hw17 + best_hw18:
            best_hw18 = h18
            print(f"    NEW BEST: HW(δe[17])={h17} + HW(δe[18])={h18} = {total}")

    print(f"\n  RESULT: Best combined HW = {best_hw17} + {best_hw18} = {best_hw17 + best_hw18}")
    print(f"  (Standard single-word delta at r=17: best HW(ΔH)=11)")
    print(f"  (Our two-barrier: combined HW on δe[17]+δe[18] = {best_hw17 + best_hw18})")
else:
    # Phase 2: Find W[15] with δe[18] = 0
    W_rest[0] = found_w0
    print(f"  Phase 2: Fix W[0]=0x{found_w0:08x}, search W[15] for δe[18]=0...")
    found_w15 = None
    for attempt in range(1 << 20):
        W_rest[15] = random.randint(0, M32)
        d17, d18 = compute_deltas(W_rest)
        assert d17 == 0, f"δe[17] should be 0 but is {d17}"
        if d18 == 0:
            found_w15 = W_rest[15]
            print(f"  FOUND! W[15]=0x{found_w15:08x} after {attempt+1} attempts")
            print(f"  δe[17]=0 AND δe[18]=0 SIMULTANEOUSLY!")
            break
        if attempt % 100000 == 99999:
            print(f"    {attempt+1} attempts...")

    if found_w15 is None:
        print(f"  Not found in 1M attempts")

# ============================================================
# TEST 4: Correlation between δe[17] and δe[18]
# ============================================================
print("\n" + "=" * 70)
print("TEST 4: Correlation between δe[17] and δe[18]")
print("=" * 70)

corr_sum = 0
n_corr = 10000
for _ in range(n_corr):
    W_test = [random.randint(0, M32) for _ in range(16)]
    d17, d18 = compute_deltas(W_test)
    # Correlation proxy: does small δe[17] predict small δe[18]?
    h17 = bin(d17).count('1')
    h18 = bin(d18).count('1')
    corr_sum += (h17 - 16) * (h18 - 16)

corr = corr_sum / n_corr / (8 * 8)  # normalize by std²
print(f"  Correlation(HW(δe[17]), HW(δe[18])): {corr:.4f}")
print(f"  (0 = independent, ±1 = fully correlated)")
