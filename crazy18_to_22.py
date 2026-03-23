#!/usr/bin/env python3
"""
CRAZY-18 to CRAZY-22: Five parallel insane ideas on verified foundation.

CRAZY-18: Carry correction size in the 128-bit barrier system
  → If linearized 128×128 system solution is CLOSE to real, cost << 2^64

CRAZY-19: Can we LEARN the barrier function with a simple model?
  → Train linear/quadratic model: W[12..15] → De17
  → If R² > 0.5, model-guided search beats random

CRAZY-20: Recurrent inversion — don't seek De17=0, seek De17=TARGET
  → If we can hit ANY target cheaply, chain: De17=T1, De18=T2 where T1,T2 cancel

CRAZY-21: De-profile collisions — don't need De=0, need De(M1)=De(M2)
  → Birthday on De VALUES, not on De=0. Cost: 2^16 per round!

CRAZY-22: 3-block attack with IV control
  → Block 1: create favorable H1 (via Joux multicollision)
  → Block 2: Wang chain with H1 as IV (more neutral bits?)
  → Block 3: finish collision
"""

import random
import time
import math

MASK = 0xFFFFFFFF
K = [0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2]
H0 = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&MASK
def Sig0(x): return rotr(x,2)^rotr(x,13)^rotr(x,22)
def Sig1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def sig0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def sig1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def Ch(e,f,g): return ((e&f)^(~e&g))&MASK
def Maj(a,b,c): return ((a&b)^(a&c)^(b&c))&MASK
def add32(*a):
    s=0
    for x in a: s=(s+x)&MASK
    return s
def hw(x): return bin(x&MASK).count('1')

def sha_round(st, w, k):
    a,b,c,d,e,f,g,h = st
    T1 = add32(h,Sig1(e),Ch(e,f,g),k,w)
    T2 = add32(Sig0(a),Maj(a,b,c))
    return [add32(T1,T2),a,b,c,add32(d,T1),e,f,g]

def expand_W(W16):
    W = list(W16)
    for i in range(16,64):
        W.append(add32(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    return W

def wang_chain(msg, iv):
    W = list(msg); Wp = list(W); Wp[0] ^= 0x80000000
    s = list(iv); sp = list(iv)
    s = sha_round(s,W[0],K[0]); sp = sha_round(sp,Wp[0],K[0])
    for t in range(1,16):
        a,b,c,d,e,f,g,h = s; a2,b2,c2,d2,e2,f2,g2,h2 = sp
        tp = add32(h,Sig1(e),Ch(e,f,g),K[t])
        tp2 = add32(h2,Sig1(e2),Ch(e2,f2,g2),K[t])
        target = add32(d,tp,W[t])
        Wp[t] = (target - d2 - tp2) & MASK
        s = sha_round(s,W[t],K[t]); sp = sha_round(sp,Wp[t],K[t])
    return Wp

def get_de_multi(msg, iv, rounds):
    Wp = wang_chain(msg, iv)
    We = expand_W(msg); Wpe = expand_W(Wp)
    s1 = list(iv); s2 = list(iv)
    results = {}
    for t in range(max(rounds)):
        s1 = sha_round(s1, We[t], K[t])
        s2 = sha_round(s2, Wpe[t], K[t])
        if t+1 in rounds:
            results[t+1] = s1[4] ^ s2[4]
    return results

def compress(msg, iv, nr):
    W = expand_W(msg)
    s = list(iv)
    for t in range(nr): s = sha_round(s, W[t], K[t])
    return [add32(s[i], iv[i]) for i in range(8)]

random.seed(0xA115)
t0 = time.time()
iv = list(H0)

# ══════════════════════════════════════════════════════════════════
print("=" * 72)
print("CRAZY-18: Carry correction size in 128-bit barrier system")
print("=" * 72)
# Linearize: replace + with ^ in Wang chain + barrier computation
# Measure: HW(real_De17 XOR linear_De17) = carry perturbation

def wang_chain_linear(msg, iv):
    """Wang chain with + replaced by ^"""
    W = list(msg); Wp = list(W); Wp[0] ^= 0x80000000
    s = list(iv); sp = list(iv)
    # Round 0
    a,b,c,d,e,f,g,h = s
    T1 = h^Sig1(e)^Ch(e,f,g)^K[0]^W[0]; T2 = Sig0(a)^Maj(a,b,c)
    s = [T1^T2,a,b,c,d^T1,e,f,g]
    a2,b2,c2,d2,e2,f2,g2,h2 = sp
    T1p = h2^Sig1(e2)^Ch(e2,f2,g2)^K[0]^Wp[0]; T2p = Sig0(a2)^Maj(a2,b2,c2)
    sp = [T1p^T2p,a2,b2,c2,d2^T1p,e2,f2,g2]
    for t in range(1,16):
        a,b,c,d,e,f,g,h = s; a2,b2,c2,d2,e2,f2,g2,h2 = sp
        tp = h^Sig1(e)^Ch(e,f,g)^K[t]; tp2 = h2^Sig1(e2)^Ch(e2,f2,g2)^K[t]
        target = d^tp^W[t]
        Wp[t] = target^d2^tp2
        T1 = tp^W[t]; T2 = Sig0(a)^Maj(a,b,c)
        s = [T1^T2,a,b,c,d^T1,e,f,g]
        T1p = tp2^Wp[t]; T2p = Sig0(a2)^Maj(a2,b2,c2)
        sp = [T1p^T2p,a2,b2,c2,d2^T1p,e2,f2,g2]
    return Wp

N = 1000
carry_hws = []
for _ in range(N):
    msg = [random.getrandbits(32) for _ in range(16)]
    # Real De17
    de_real = get_de_multi(msg, iv, [17])
    # Linear De17 (using linearized Wang chain)
    Wp_lin = wang_chain_linear(msg, iv)
    We_lin = list(msg)
    for i in range(16,18):
        We_lin.append(sig1(We_lin[i-2])^We_lin[i-7]^sig0(We_lin[i-15])^We_lin[i-16])
    Wpe_lin = list(Wp_lin)
    for i in range(16,18):
        Wpe_lin.append(sig1(Wpe_lin[i-2])^Wpe_lin[i-7]^sig0(Wpe_lin[i-15])^Wpe_lin[i-16])
    s1 = list(iv); s2 = list(iv)
    for t in range(17):
        a,b,c,d,e,f,g,h = s1
        T1 = h^Sig1(e)^Ch(e,f,g)^K[t]^We_lin[t]; T2 = Sig0(a)^Maj(a,b,c)
        s1 = [T1^T2,a,b,c,d^T1,e,f,g]
        a,b,c,d,e,f,g,h = s2
        T1 = h^Sig1(e)^Ch(e,f,g)^K[t]^Wpe_lin[t]; T2 = Sig0(a)^Maj(a,b,c)
        s2 = [T1^T2,a,b,c,d^T1,e,f,g]
    de_lin = s1[4] ^ s2[4]
    # Carry perturbation = difference between real and linear
    carry_perturb = hw(de_real[17] ^ de_lin)
    carry_hws.append(carry_perturb)

avg_cp = sum(carry_hws)/len(carry_hws)
min_cp = min(carry_hws)
print(f"  Carry perturbation HW(real XOR linear): avg={avg_cp:.1f}, min={min_cp}")
print(f"  Expected if random: 16.0")
if avg_cp < 12:
    print(f"  ALIVE: carry correction is SMALL ({avg_cp:.1f} < 16)!")
    print(f"  Linear solution + local search could work!")
else:
    print(f"  DEAD: carry correction ≈ random ({avg_cp:.1f} ≈ 16)")


# ══════════════════════════════════════════════════════════════════
print(f"\n{'='*72}")
print("CRAZY-19: Can we LEARN the barrier? (linear model)")
print("=" * 72)
# Build dataset: (W[14] bits) → De17 bits, fit linear model over GF(2)
# If any De17 bit is predictable from W[14] bits, model-guided search helps

base_msg = [random.getrandbits(32) for _ in range(16)]
N_data = 2000
# Collect data
X = []  # W[14] values
Y = []  # De17 values
for _ in range(N_data):
    msg = list(base_msg)
    msg[14] = random.getrandbits(32)
    de17 = get_de_multi(msg, iv, [17])[17]
    X.append(msg[14])
    Y.append(de17)

# For each output bit of De17, measure correlation with each input bit of W[14]
max_corr = 0
best_pair = None
for out_bit in range(32):
    for in_bit in range(32):
        match = sum(1 for i in range(N_data) if ((X[i]>>in_bit)&1) == ((Y[i]>>out_bit)&1))
        corr = abs(match/N_data - 0.5) * 2
        if corr > max_corr:
            max_corr = corr
            best_pair = (in_bit, out_bit)

print(f"  Max bit-to-bit correlation: {max_corr:.4f}")
print(f"  Best pair: W[14] bit {best_pair[0]} → De17 bit {best_pair[1]}")
print(f"  Expected random: ~{1/math.sqrt(N_data):.4f}")
if max_corr > 3/math.sqrt(N_data):
    print(f"  ALIVE: significant correlation found!")
else:
    print(f"  DEAD: correlation ≈ noise. Barrier is unpredictable.")


# ══════════════════════════════════════════════════════════════════
print(f"\n{'='*72}")
print("CRAZY-21: De-PROFILE collisions (De1=De2 instead of De=0)")
print("=" * 72)
# KEY IDEA: We don't need De17=0. We need De17(M1) = De17(M2).
# This is a BIRTHDAY on De17 values: cost 2^16 (not 2^32)!
# Two different messages with same De17 → differential cancels at round 17.

# But wait — for a collision, we need H(M1)=H(M2), not just De17 match.
# If De17(M1)=De17(M2), does that help?
# It means: the state DIFFERENCE at round 17 is the same for both pairs.
# This is NOT the same as a collision... unless we can chain it.

# Test: collect De17 values, find birthday matches
N_birthday = 50000
de17_dict = {}
collisions = 0
for trial in range(N_birthday):
    msg = [random.getrandbits(32) for _ in range(16)]
    de17 = get_de_multi(msg, iv, [17])[17]
    if de17 in de17_dict:
        collisions += 1
    else:
        de17_dict[de17] = msg

expected = N_birthday**2 / (2 * 2**32)
print(f"  De17 birthday matches in {N_birthday} messages: {collisions}")
print(f"  Expected: {expected:.1f}")
print(f"  (This is birthday on 32 bits, cost ≈ 2^16)")

# Now: for matched De17, check De18, De19...
if collisions > 0:
    # Find actual matches
    de17_dict2 = {}
    match_pairs = []
    for trial in range(N_birthday):
        msg = [random.getrandbits(32) for _ in range(16)]
        de17 = get_de_multi(msg, iv, [17])[17]
        if de17 in de17_dict2:
            match_pairs.append((de17_dict2[de17], msg))
        de17_dict2[de17] = msg

    print(f"  Found {len(match_pairs)} De17-matching pairs")
    for i, (m1, m2) in enumerate(match_pairs[:5]):
        de1 = get_de_multi(m1, iv, [17,18,19,20])
        de2 = get_de_multi(m2, iv, [17,18,19,20])
        print(f"    Pair {i}: De17 match={de1[17]==de2[17]}, "
              f"De18 hw diff={hw(de1[18]^de2[18])}")


# ══════════════════════════════════════════════════════════════════
print(f"\n{'='*72}")
print("CRAZY-22: Multi-block with controlled IV")
print("=" * 72)
# In multi-block SHA-256, H1 = compress(H0, M1).
# If we choose M1 wisely, H1 can have specific properties
# that make block 2 easier to attack.
# Question: does the NUMBER OF NEUTRAL BITS depend on IV?

N_iv = 10
neutral_counts = []
for _ in range(N_iv):
    # Random IV (simulating H1 from block 1)
    rand_iv = [random.getrandbits(32) for _ in range(8)]
    msg = [random.getrandbits(32) for _ in range(16)]
    Wp = wang_chain(msg, rand_iv)

    count = 0
    for k in range(1, 16):
        for b in range(32):
            Wf = list(msg); Wf[k] ^= (1<<b)
            Wpf = list(Wp); Wpf[k] ^= (1<<b)
            s1 = list(rand_iv); s2 = list(rand_iv)
            s1 = sha_round(s1,Wf[0],K[0]); s2 = sha_round(s2,Wpf[0],K[0])
            ok = True
            for t in range(1,16):
                s1 = sha_round(s1,Wf[t],K[t]); s2 = sha_round(s2,Wpf[t],K[t])
                if t >= 2 and (s1[4]-s2[4])&MASK != 0:
                    ok = False; break
            if ok: count += 1
    neutral_counts.append(count)
    print(f"  IV #{_}: {count} neutral bits")

avg_n = sum(neutral_counts)/len(neutral_counts)
std_n = (sum((x-avg_n)**2 for x in neutral_counts)/len(neutral_counts))**0.5
print(f"\n  Average: {avg_n:.1f} ± {std_n:.1f}")
print(f"  Standard IV: 82.6")
print(f"  Range: {min(neutral_counts)} — {max(neutral_counts)}")

if max(neutral_counts) > 120:
    print(f"  ALIVE: some IVs give {max(neutral_counts)} neutrals (vs 82.6 standard)!")
    print(f"  Multi-block: optimize IV in block 1 for max neutrals in block 2!")
elif avg_n > 90:
    print(f"  ANOMALY: random IVs give more neutrals on average")
else:
    print(f"  Neutral count ≈ independent of IV")


# ══════════════════════════════════════════════════════════════════
print(f"\n{'='*72}")
print("CRAZY-20: Recurrent inversion — hit TARGET instead of 0")
print("=" * 72)
# Instead of De17=0 (cost 2^32), what if we seek De17=TARGET
# where TARGET is chosen to make De18 easier?
# Specifically: find TARGET such that De18(given De17=TARGET) is more likely 0.
# This requires understanding the CONDITIONAL distribution P(De18 | De17).

# We already know barriers are independent (CRAZY-8), so this shouldn't help.
# But let's verify on the FREE WORD system where barriers are COUPLED.
# Maybe P(De18=0 | De17=specific_value) varies?

base_msg = [random.getrandbits(32) for _ in range(16)]
N_sample = 20000
de17_de18_pairs = []
for _ in range(N_sample):
    msg = list(base_msg)
    msg[14] = random.getrandbits(32)
    msg[15] = random.getrandbits(32)
    des = get_de_multi(msg, iv, [17, 18])
    de17_de18_pairs.append((des[17], des[18]))

# Bucket De17 by low bits, check if De18 distribution changes
for mask_bits in [1, 2, 4, 8]:
    mask = (1 << mask_bits) - 1
    buckets = {}
    for d17, d18 in de17_de18_pairs:
        key = d17 & mask
        if key not in buckets:
            buckets[key] = []
        buckets[key].append(hw(d18))

    # Check if mean HW(De18) varies across buckets
    means = [(k, sum(v)/len(v), len(v)) for k, v in buckets.items() if len(v) > 10]
    if means:
        min_mean = min(m for _,m,_ in means)
        max_mean = max(m for _,m,_ in means)
        print(f"  Bucket by De17 low {mask_bits} bits: "
              f"mean HW(De18) range = [{min_mean:.2f}, {max_mean:.2f}] "
              f"(Δ={max_mean-min_mean:.2f})")


# ══════════════════════════════════════════════════════════════════
print(f"\n{'='*72}")
print("GRAND SUMMARY")
print("=" * 72)
print(f"  Runtime: {time.time()-t0:.1f}s")
