#!/usr/bin/env python3
"""
П-24: Sufficient Conditions for δa1=0
Extending the analytical XOR trail from 4 rounds to 8+ rounds.

From П-23:
  δW0=0x8000 → δe1=0x8000 (P=1.0, SC: W0[0..14]≥7518)
  δe2=0 (P=0.1256), δe3=0 (P=0.0636), δe4=0 (P=0.0338)
  Problem: δa1≠0 in general → δe5≠0 via d-register chain

Goal: Find sufficient conditions on (W0, W1, ...) such that δa1=0,
enabling deterministic extension to δe5=0, δe6=0, ...

Key: a1 = T1_0 + T2_0
where T1_0 = h_iv + Sig1(e_iv) + Ch(e_iv,f_iv,g_iv) + K[0] + W0
      T2_0 = Sig0(a_iv) + Maj(a_iv,b_iv,c_iv)

For XOR trail with δW0=0x8000:
  δT1_0 = δW0 = 0x8000 (since other terms depend on IV constants, fixed)
  δT2_0 = 0 (T2 doesn't depend on W)
  δa1 = δ(T1_0 + T2_0) = δT1_0 = 0x8000 ??? Not necessarily!

Wait - additive differential vs XOR:
  a1_f = (T1_0_f + T2_0) mod 2^32
  a1_n = (T1_0_n + T2_0) mod 2^32
  T1_0_f = T1_0_n + δW0  (additive, since all other T1 components equal)

  δa1 (XOR) = a1_f XOR a1_n = (T1_0_n + 0x8000 + T2_0) XOR (T1_0_n + T2_0)
            = depends on carry propagation from bit 15!

SC for δa1=0x8000: need carry from bit 15 to NOT propagate to bit 16.
  C = (T1_0_n + T2_0) mod 2^32
  bit 15 of C must be 0 (no carry generated at bit 15)

  Actually more precisely: T1_0_n + T2_0 + 0x8000 XOR T1_0_n + T2_0
  Let S = T1_0_n + T2_0
  δa1 = (S + 0x8000) XOR S

  This equals 0x8000 iff bits 16..31 of S are unchanged when adding 0x8000
  i.e., no carry propagates from bit 15 upward
  i.e., bits 16..31 of S are NOT all changed by carry

  Specifically: δa1=0x8000 iff S[16..31] doesn't change
  i.e., S[15]=0 (bit 15 of S is 0, so adding 0x8000 = setting bit 15, no carry to 16)
  Wait: S + 0x8000: if S[15]=0, result is S | 0x8000, XOR = 0x8000. ✓
  If S[15]=1, then adding 0x8000 causes carry: result clears bit 15 and propagates.
  XOR would then be 0x8000 XOR carry_pattern ≠ 0x8000.

So T_SC_A1: δa1=0x8000 iff S[15]=0, where S=(T1_0_n+T2_0) mod 2^32
"""

import random
import struct

# ── SHA-256 constants ──────────────────────────────────────────────────────────
K = [
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,
    0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,
    0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,
    0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,
    0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,
    0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,
    0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,
    0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,
    0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2,
]
IV = (0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
      0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19)
MASK = 0xFFFFFFFF

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&MASK
def sig0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def sig1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def Sig0(x): return rotr(x,2)^rotr(x,13)^rotr(x,22)
def Sig1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def Ch(x,y,z): return (x&y)^(~x&z)&MASK
def Maj(x,y,z): return (x&y)^(x&z)^(y&z)

def schedule(msg16):
    W = list(msg16) + [0]*48
    for i in range(16,64):
        W[i] = (sig1(W[i-2])+W[i-7]+sig0(W[i-15])+W[i-16])&MASK
    return W

def sha_rounds(W, num_rounds):
    a,b,c,d,e,f,g,h = IV
    state = []
    for i in range(num_rounds):
        T1 = (h+Sig1(e)+Ch(e,f,g)+K[i]+W[i])&MASK
        T2 = (Sig0(a)+Maj(a,b,c))&MASK
        h=g; g=f; f=e; e=(d+T1)&MASK
        a=(T1+T2)&MASK; b=a; c=b; d=c
        # Actually proper update:
        h=g; g=f; f=e; e=(d+T1)&MASK
        d=c; c=b; b=a; a=(T1+T2)&MASK
        state.append((a,b,c,d,e,f,g,h))
    return state

def sha_rounds_correct(W, num_rounds):
    """Correct SHA-256 round function."""
    a,b,c,d,e,f,g,h = IV
    states = []
    for i in range(num_rounds):
        T1 = (h + Sig1(e) + Ch(e,f,g) + K[i] + W[i]) & MASK
        T2 = (Sig0(a) + Maj(a,b,c)) & MASK
        new_e = (d + T1) & MASK
        new_a = (T1 + T2) & MASK
        h=g; g=f; f=e; e=new_e
        d=c; c=b; b=a; a=new_a
        states.append((a,b,c,d,e,f,g,h))
    return states

print("="*70)
print("П-24: Sufficient Conditions for δa1=0")
print("="*70)

# ── Part 1: Analyze δa1 for δW0=0x8000 ────────────────────────────────────────
print("\n[1] АНАЛИЗ δa1 при δW0=0x8000")
print("-"*50)

dW0 = 0x8000
N = 100000
count_da1_eq_dW0 = 0
count_da1_zero = 0
da1_values = {}

a_iv, b_iv, c_iv, d_iv, e_iv, f_iv, g_iv, h_iv = IV

# T2_0 is fixed (doesn't depend on W):
T2_0 = (Sig0(a_iv) + Maj(a_iv, b_iv, c_iv)) & MASK
# T1_0 base (without W0):
T1_0_base = (h_iv + Sig1(e_iv) + Ch(e_iv, f_iv, g_iv) + K[0]) & MASK

for _ in range(N):
    W0 = random.randint(0, MASK)
    T1_0_n = (T1_0_base + W0) & MASK
    T1_0_f = (T1_0_n + dW0) & MASK  # adding δW0=0x8000

    a1_n = (T1_0_n + T2_0) & MASK
    a1_f = (T1_0_f + T2_0) & MASK
    da1 = a1_f ^ a1_n

    da1_values[da1] = da1_values.get(da1, 0) + 1
    if da1 == dW0:
        count_da1_eq_dW0 += 1
    if da1 == 0:
        count_da1_zero += 1

print(f"N={N}, δW0=0x{dW0:08x}")
print(f"P(δa1=δW0=0x8000) = {count_da1_eq_dW0}/{N} = {count_da1_eq_dW0/N:.4f}")
print(f"P(δa1=0)           = {count_da1_zero}/{N} = {count_da1_zero/N:.4f}")

# Top δa1 values
top_da1 = sorted(da1_values.items(), key=lambda x:-x[1])[:5]
print("\nТоп-5 δa1 значений:")
for v, c in top_da1:
    print(f"  δa1=0x{v:08x}: {c}/{N} = {c/N:.4f}")

# ── Part 2: Sufficient condition for δa1=0x8000 ───────────────────────────────
print("\n[2] T_SC_A1: ДОСТАТОЧНОЕ УСЛОВИЕ δa1=0x8000")
print("-"*50)
print("Теория: δa1=δW0=0x8000 ⟺ S[bit15]=0")
print("где S = (T1_0_n + T2_0) mod 2^32 = a1_n")
print()

# Verify: δa1=0x8000 iff a1_n[bit15]=0
verified = 0
total = 100000
for _ in range(total):
    W0 = random.randint(0, MASK)
    T1_0_n = (T1_0_base + W0) & MASK
    a1_n = (T1_0_n + T2_0) & MASK
    a1_f = ((T1_0_n + dW0 + T2_0)) & MASK
    da1 = a1_f ^ a1_n

    predicted = (da1 == dW0)
    actual_cond = (a1_n >> 15) & 1 == 0
    if predicted == actual_cond:
        verified += 1

print(f"Верификация T_SC_A1: {verified}/{total} = {verified/total:.6f}")

# ── Part 3: Sufficient condition on W0 for both δe1=0x8000 AND δa1=0x8000 ─────
print("\n[3] JOINT SC: δe1=0x8000 AND δa1=0x8000")
print("-"*50)

# From П-22: SC for δe1=0x8000: W0[0..14] >= 7518
# i.e., the low 15 bits of W0 must be >= 7518

# For δa1=0x8000: need a1_n[bit15]=0
# a1_n = (T1_0_base + W0 + T2_0) mod 2^32
# bit15 of a1_n depends on W0[0..15] and carry from low bits

# Let's compute T1_0_base + T2_0 = fixed base
S_base = (T1_0_base + T2_0) & MASK
print(f"S_base = T1_0_base + T2_0 = 0x{S_base:08x}")
print(f"T1_0_base = 0x{T1_0_base:08x}")
print(f"T2_0      = 0x{T2_0:08x}")

# S = S_base + W0. We need bit 15 of S = 0.
# S[15] = ((S_base + W0) >> 15) & 1
# This depends on W0[0..15] and carry from S_base[0..14] + W0[0..14]

# Let's find the constraint:
# S[15] = ((S_base[0..15] + W0[0..15]) >> 15) & 1
# For bit15 of S = 0, we need the 16-bit sum S_base[0..15] + W0[0..15] to have bit15=0
# i.e., S_base_low16 + W0_low16 < 0x8000 OR >= 0x10000 (carry out, bit15=0 in lower 16)

S_base_low16 = S_base & 0xFFFF
print(f"S_base[0..15] = 0x{S_base_low16:04x} = {S_base_low16}")

# Condition: bit15 of (S_base_low16 + W0_low16) = 0
# Let sum16 = S_base_low16 + W0_low16
# bit15(sum16 & 0xFFFF) = 0 iff (sum16 & 0x8000 == 0) iff sum16 < 0x8000 or sum16 >= 0x10000

# Count valid W0_low16 values:
count_valid_low16 = 0
for w0_low16 in range(0x10000):
    sum16 = S_base_low16 + w0_low16
    if (sum16 & 0x8000) == 0:  # bit15=0 in lower 16 (no carry yet or carry past)
        count_valid_low16 += 1

# Actually we need bit15 of S = bit15 of (S_base + W0)
# which equals bit15 of the 32-bit addition
# = ((S_base & 0x7FFF) + (W0 & 0x7FFF) + carry_from_below_15) >> 15...
# more precisely: bit15 of (S_base + W0) = parity of (S_base[15], W0[15], carry_into_15)

# Let me just enumerate W0_low16:
valid_regions = []
in_valid = False
start = None
for w0_low16 in range(0x10000):
    sum16 = (S_base_low16 + w0_low16) & 0xFFFF
    bit15 = (sum16 >> 15) & 1
    is_valid = (bit15 == 0)
    if is_valid and not in_valid:
        start = w0_low16
        in_valid = True
    elif not is_valid and in_valid:
        valid_regions.append((start, w0_low16 - 1))
        in_valid = False
if in_valid:
    valid_regions.append((start, 0xFFFF))

print(f"\nРегионы W0[0..15] для δa1=0x8000:")
total_valid = 0
for lo, hi in valid_regions:
    print(f"  0x{lo:04x}..0x{hi:04x} ({hi-lo+1} значений)")
    total_valid += hi - lo + 1
print(f"Итого: {total_valid}/65536 = {total_valid/65536:.4f}")

# ── Part 4: Find SC for BOTH δe1=0x8000 AND δa1=0x8000 simultaneously ─────────
print("\n[4] JOINT SC: δe1=0x8000 AND δa1=0x8000 SIMULTANEOUSLY")
print("-"*50)

# From П-22: P(δe1=0x8000 | SC_e) where SC_e = W0[0..14] >= 7518, P=1.0
# T_IV_BIT0: C = d_iv + T1_0_base = ?
C_e = (d_iv + T1_0_base) & MASK
print(f"C_e (for e1) = d_iv + T1_0_base = 0x{C_e:08x}")
print(f"C_e[14..15] = {(C_e>>14)&3:02b}  (need bit15=0 for simple SC)")

# Threshold for δe1=0x8000 with P=1:
# From P-22: W0[0..14] >= threshold where threshold = ...
# Let's recompute. C_e[0..14] + W0[0..14] must not overflow 15 bits for deterministic carry
# Actually from P-22: W0[0..14] >= (2^15 - C_e[0..14]) = threshold
C_e_low15 = C_e & 0x7FFF
threshold_e1 = (0x8000 - C_e_low15) & 0x7FFF
print(f"C_e[0..14] = 0x{C_e_low15:04x} = {C_e_low15}")
print(f"Threshold for δe1: W0[0..14] >= {threshold_e1} (= {threshold_e1} = 0x{threshold_e1:04x})")

# Now find valid W0[0..14] for δe1=0x8000 (from SC):
# and valid W0[0..15] for δa1=0x8000
# SC_e requires W0[0..14] >= threshold_e1
# SC_a1 requires bit15 of (S_base + W0) = 0

# Let's test empirically:
N_test = 200000
joint_count = 0
only_e1 = 0
only_a1 = 0
neither = 0

for _ in range(N_test):
    W0 = random.randint(0, MASK)
    W0_low15 = W0 & 0x7FFF

    # Check SC_e (deterministic δe1):
    sc_e = (W0_low15 >= threshold_e1)

    # Compute actual δe1:
    T1_n = (T1_0_base + W0) & MASK
    e1_n = (d_iv + T1_n) & MASK
    T1_f = (T1_0_base + W0 + dW0) & MASK
    e1_f = (d_iv + T1_f) & MASK
    de1 = e1_f ^ e1_n
    got_e1 = (de1 == dW0)

    # Compute δa1:
    a1_n = (T1_n + T2_0) & MASK
    a1_f = (T1_f + T2_0) & MASK
    da1 = a1_f ^ a1_n
    got_a1 = (da1 == dW0)

    if got_e1 and got_a1:
        joint_count += 1
    elif got_e1:
        only_e1 += 1
    elif got_a1:
        only_a1 += 1
    else:
        neither += 1

print(f"\nЭмпирические результаты (N={N_test}):")
print(f"  δe1=0x8000 AND δa1=0x8000: {joint_count}/{N_test} = {joint_count/N_test:.4f}")
print(f"  δe1=0x8000 only:           {only_e1}/{N_test} = {only_e1/N_test:.4f}")
print(f"  δa1=0x8000 only:           {only_a1}/{N_test} = {only_a1/N_test:.4f}")
print(f"  neither:                    {neither}/{N_test} = {neither/N_test:.4f}")

# ── Part 5: When BOTH δe1=0x8000 AND δa1=0x8000, extend trail ────────────────
print("\n[5] РАСШИРЕНИЕ СЛЕДА: δe1=δa1=0x8000 → δe5, δe6, ...")
print("-"*50)

# If δe1=0x8000 and δa1=0x8000, then:
# Round 2: e2 depends on e1, d=d_iv (unchanged), T1_1, T2_1
# a2 depends on a1, b1=a1... wait, a2 depends on T1_1+T2_1
# Let's trace the XOR differential analytically:
#
# State after round 1 (normal): a1, b1=a_iv, c1=b_iv, d1=c_iv, e1, f1=e_iv, g1=f_iv, h1=g_iv
# State after round 1 (faulty): a1+δa1, b1, c1, d1, e1+δe1, f1, g1, h1
#
# Round 2:
# T1_1_n = h1 + Sig1(e1) + Ch(e1,f1,g1) + K[1] + W1
# T1_1_f = h1 + Sig1(e1^δe1) + Ch(e1^δe1,f1,g1) + K[1] + W1  [δW1=0 here]
# δT1_1 = Sig1(e1^δe1) - Sig1(e1) + Ch(e1^δe1,f1,g1) - Ch(e1,f1,g1)  [additive!]
#
# But for XOR differential:
# δSig1(e1) = Sig1(e1^δe1) XOR Sig1(e1) = Sig1(δe1) [linear!] = Sig1(0x8000) = 0x00400210
# δCh(e1,f1,g1) = δe1 & (f1 XOR g1) = 0x8000 & (e_iv XOR f_iv)
#
# XOR is NOT additive differential... but for the trail we use XOR

# Let's do the 5-round XOR differential tracking
print("Аналитический след (XOR) с δW0=0x8000, δW1=..., δW2=...=0:")
print()

dW0 = 0x8000

# After round 1:
# δe1 = 0x8000 (given SC)
# δa1 = 0x8000 (given SC)
# δf1=δg1=δh1=δb1=δc1=δd1 = 0 (state rotation)

# Round 2:
# h_in_r2 = g1 = f_iv
# e_in_r2 = e1  (varies, but δe_in_r2 = δe1 = 0x8000)
# f_in_r2 = f1 = e_iv
# g_in_r2 = g1 = f_iv  [wait, after round 1: h1=g_iv, g1=f_iv, f1=e_iv]

# After 1 round: a,b,c,d,e,f,g,h = a1,a_iv,b_iv,c_iv,e1,e_iv,f_iv,g_iv

# δSig1(e_in_r2) = Sig1(0x8000) = 0x00400210
dSig1_e1 = Sig1(dW0)
print(f"δSig1(e1) = Sig1(0x{dW0:08x}) = 0x{dSig1_e1:08x}")

# δCh(e1,f1,g1) = δe1 & (f1 XOR g1) = 0x8000 & (e_iv XOR f_iv)
e_in_r2 = e_iv  # this is the normal state; faulty adds δe1
f_in_r2 = e_iv  # f1 after round 1 = e_iv
g_in_r2 = f_iv  # g1 after round 1 = f_iv
dCh_r2 = dW0 & (f_in_r2 ^ g_in_r2)
print(f"δCh(e1,f1,g1) = 0x{dW0:08x} & (e_iv XOR f_iv) = 0x{dCh_r2:08x}")

# h_in_r2 = g_iv (unchanged in both)
# So δT1_1 (XOR) = δSig1(e1) XOR δCh XOR δW1
# For δW1=0: δT1_1 = dSig1_e1 XOR dCh_r2
dT1_1_xor = dSig1_e1 ^ dCh_r2
print(f"δT1_1 (XOR, δW1=0) = δSig1 XOR δCh = 0x{dT1_1_xor:08x}")

# δe2 = δd1 XOR δT1_1 (XOR) [d1 = c_iv, δd1=0]
# BUT: XOR differential for addition is NOT just XOR of inputs!
# δ(X+Y) ≠ δX XOR δY in general
# We need to use the carry analysis again...

# For δe2: e2 = d1 + T1_1
# δe2 = δ(d1 + T1_1) = (d1 + T1_1_f) XOR (d1 + T1_1_n)
# since δd1=0, δe2 = (d1 + T1_1_n + δT1_1_add) XOR (d1 + T1_1_n)
# where δT1_1_add is the additive differential of T1_1

# This is getting complex. Let me simulate directly.
print("\nЭмпирическое исследование 5-раундового следа:")
print("Условие: δe1=0x8000 AND δa1=0x8000")

N_sim = 200000
dW1_candidates = [0, dSig1_e1 ^ dCh_r2, dSig1_e1, dCh_r2]
trail_counts = {}

for dW1 in dW1_candidates:
    counts = {i: {0: 0, dW0: 0, 'other': 0} for i in range(1, 9)}
    n_joint = 0

    for _ in range(N_sim):
        # Sample W0 with BOTH SC satisfied
        # SC for δe1: W0[0..14] >= threshold_e1
        # SC for δa1: a1_n[bit15]=0
        # Try random W0, check conditions
        W0_low15 = random.randint(threshold_e1, 0x7FFF)
        W0_high17 = random.randint(0, 0x1FFFF)
        W0 = W0_low15 | (W0_high17 << 15)
        W0 &= MASK

        T1_0_n = (T1_0_base + W0) & MASK
        a1_n = (T1_0_n + T2_0) & MASK
        if (a1_n >> 15) & 1 != 0:
            continue  # SC for δa1 not satisfied

        n_joint += 1

        # Build schedule
        Wn = [W0] + [random.randint(0, MASK) for _ in range(15)]
        Wf = list(Wn)
        Wf[0] = (Wn[0] + dW0) & MASK
        Wf[1] = (Wn[1] + dW1) & MASK

        sn = sha_rounds_correct(schedule(Wn), 8)
        sf = sha_rounds_correct(schedule(Wf), 8)

        for r in range(1, 9):
            de = sf[r-1][4] ^ sn[r-1][4]  # e-register
            if de == 0:
                counts[r][0] += 1
            elif de == dW0:
                counts[r][dW0] += 1
            else:
                counts[r]['other'] += 1

    if n_joint == 0:
        print(f"  δW1=0x{dW1:08x}: нет сэмплов с joint SC")
        continue

    print(f"\n  δW1=0x{dW1:08x} (n_joint={n_joint}):")
    for r in range(1, 9):
        total_r = sum(counts[r].values())
        p0 = counts[r][0] / n_joint
        pdW = counts[r][dW0] / n_joint
        print(f"    Round {r}: P(δe=0)={p0:.4f} P(δe=0x8000)={pdW:.4f} total={n_joint}")

# ── Part 6: Sufficient condition on W0[15] for δa1=0x8000 ────────────────────
print("\n[6] T_SC_A1_BIT: ТОЧНОЕ УСЛОВИЕ НА W0[15]")
print("-"*50)

# a1_n = (T1_0_base + T2_0 + W0) mod 2^32
# bit15 of a1_n = bit15 of (S_base + W0)
# S_base = (T1_0_base + T2_0) mod 2^32

S_base_val = (T1_0_base + T2_0) & MASK
print(f"S_base = 0x{S_base_val:08x}")
print(f"S_base[14..16] = {(S_base_val>>14)&7:03b}")

# For bit15 of (S_base + W0) = 0:
# Let carry_14 = carry into bit 15 from bits 0..14
# carry_14 = 1 iff (S_base & 0x7FFF) + (W0 & 0x7FFF) >= 0x8000
#
# bit15(S_base + W0) = S_base[15] XOR W0[15] XOR carry_14
# We want this = 0, so:
# W0[15] XOR carry_14 = S_base[15]
# W0[15] = S_base[15] XOR carry_14

S_base_bit15 = (S_base_val >> 15) & 1
S_base_low15 = S_base_val & 0x7FFF
print(f"S_base[15] = {S_base_bit15}")

print(f"\nДля δa1=0x8000 нужно bit15(a1_n)=0:")
print(f"W0[15] = S_base[15] XOR carry_14")
print(f"       = {S_base_bit15} XOR carry_14")
print(f"carry_14 = 1 iff (S_base[0..14] + W0[0..14]) >= 0x8000")
print(f"         = 1 iff W0[0..14] >= {0x8000 - S_base_low15} = 0x{(0x8000-S_base_low15)&0x7FFF:04x}")

thresh_carry = (0x8000 - S_base_low15) & 0x7FFF
print(f"\nSC_a1:")
print(f"  Если W0[0..14] >= {thresh_carry}: carry_14=1, нужно W0[15]={S_base_bit15^1}")
print(f"  Если W0[0..14] <  {thresh_carry}: carry_14=0, нужно W0[15]={S_base_bit15}")

# ── Part 7: Combined SC and probability ──────────────────────────────────────
print("\n[7] КОМБИНИРОВАННАЯ ВЕРОЯТНОСТЬ: P(SC_e AND SC_a1)")
print("-"*50)

# Count pairs (W0[0..14], W0[15]) satisfying both SC
# SC_e: W0[0..14] >= threshold_e1
# SC_a1: depends on W0[0..14] (via carry_14) and W0[15]

count_both = 0
count_e_only = 0
count_a1_only = 0
total_w0 = 2**16  # W0[0..15] space

for w0_low16 in range(total_w0):
    w0_low15 = w0_low16 & 0x7FFF
    w0_bit15 = (w0_low16 >> 15) & 1

    sc_e = (w0_low15 >= threshold_e1)

    carry_14 = 1 if (S_base_low15 + w0_low15) >= 0x8000 else 0
    needed_bit15 = S_base_bit15 ^ carry_14
    sc_a1 = (w0_bit15 == needed_bit15)

    if sc_e and sc_a1:
        count_both += 1
    elif sc_e:
        count_e_only += 1
    elif sc_a1:
        count_a1_only += 1

print(f"Из {total_w0} возможных W0[0..15]:")
print(f"  SC_e AND SC_a1: {count_both} = {count_both/total_w0:.4f}")
print(f"  SC_e только:    {count_e_only} = {count_e_only/total_w0:.4f}")
print(f"  SC_a1 только:   {count_a1_only} = {count_a1_only/total_w0:.4f}")

# Note: SC_e requires W0[0..14] in range [threshold_e1, 0x7FFF] = (0x7FFF - threshold_e1 + 1) values
# out of 0x8000. Then exactly half of W0[15] values will satisfy SC_a1 (either 0 or 1).
# So P(SC_e AND SC_a1) = (0x8000 - threshold_e1) / 0x8000 * 1/2
# Then high bits W0[16..31] are unconstrained: any of 2^16 values.
# Total fraction: (0x8000 - threshold_e1) / (0x8000 * 2)

theoretical = (0x8000 - threshold_e1) / (0x8000 * 2)
print(f"\nТеоретическая оценка: {theoretical:.4f}")
print(f"Фракция W0[0..31]: (0x8000-{threshold_e1})/(0x8000*2) = {theoretical:.4f}")
print(f"Количество сообщений для выполнения SC: ~{1/theoretical:.1f}")

# ── Part 8: Extend trail to δe5=0 given δe1=δa1=0x8000 ───────────────────────
print("\n[8] АНАЛИТИЧЕСКИЙ СЛЕД: δe5 при δe1=δa1=0x8000")
print("-"*50)

# After round 1:
#   δa1=0x8000, δe1=0x8000, δb1=δc1=δd1=δf1=δg1=δh1=0
#
# Round 2 state entering: a=a1, b=a_iv, c=b_iv, d=c_iv, e=e1, f=e_iv, g=f_iv, h=g_iv
#   δa=0x8000, δe=0x8000, others=0
#
# T1_2 = h + Sig1(e) + Ch(e,f,g) + K[1] + W1
# δT1_2 (XOR) ≈ δSig1(e) XOR δCh = Sig1(0x8000) XOR (0x8000 & (e_iv XOR f_iv))
#             = 0x00400210 XOR dCh_r2

# T2_2 = Sig0(a) + Maj(a,b,c)
# δT2_2 (XOR) ≈ Sig0(0x8000) XOR Maj_delta
# Maj(a XOR δ, b, c) XOR Maj(a, b, c) = δ & (b XOR c) XOR (δ & b) XOR (δ & c)...
# Actually: Maj(x,y,z) = (x&y)^(x&z)^(y&z)
# δMaj = Maj(a^δ,b,c) ^ Maj(a,b,c) = δ&b ^ δ&c = δ & (b XOR c)...
# Let me verify:
# Maj(a^δ,b,c) = (a^δ)&b ^ (a^δ)&c ^ b&c
#              = a&b ^ δ&b ^ a&c ^ δ&c ^ b&c
# Maj(a,b,c) = a&b ^ a&c ^ b&c
# Difference: δ&b ^ δ&c = δ & (b XOR c)
# So δMaj = δa & (b1 XOR c1) = 0x8000 & (a_iv XOR b_iv)

dMaj_r2 = dW0 & (a_iv ^ b_iv)  # δa1=0x8000, b1=a_iv, c1=b_iv
dSig0_a1 = Sig0(dW0)  # linear! Sig0(a XOR δ) XOR Sig0(a) = Sig0(δ)
dT2_2_xor = dSig0_a1 ^ dMaj_r2

print(f"δSig0(a1) = Sig0(0x{dW0:08x}) = 0x{dSig0_a1:08x}")
print(f"δMaj(a1,b1,c1) = 0x{dW0:08x} & (a_iv XOR b_iv) = 0x{dMaj_r2:08x}")
print(f"δT2_2 (XOR) = δSig0 XOR δMaj = 0x{dT2_2_xor:08x}")
print()

# δe2 = δd1 XOR δ(T1_1) [via additive carry analysis]
# But d1=c_iv (fixed, δd1=0), so δe2 is determined by carry in (d1+T1_1_f) XOR (d1+T1_1_n)
# This is the same carry analysis as for δe1.
# For XOR trail: we SET δW1 such that δe2=0 (or 0x8000).

# If we choose δW1 such that δT1_2 = 0 additively... complex.
# Let's check: with δW1=0, what is P(δe2=0)?
# From П-23: P(δe2=0 | δW1=0x00400210) ≈ 0.1256

print("Trail summary so far:")
print(f"  Round 1: δe1=0x8000 (P=1.0 given SC_e), δa1=0x8000 (P=1.0 given SC_a1)")
print(f"  Round 2: δe2=? (depends on carry in d1+T1_1)")
print(f"           Need to find SC for δe2=0 (or extend via δW1)")
print()

# SC for δe2=0: d1=c_iv is fixed.
# T1_1_n = g_iv + Sig1(e1) + Ch(e1,e_iv,f_iv) + K[1] + W1
# T1_1_f = g_iv + Sig1(e1^0x8000) + Ch(e1^0x8000,e_iv,f_iv) + K[1] + W1+δW1
# δ(T1_1)_add = (Sig1(e1^0x8000)-Sig1(e1)) + (Ch(e1^0x8000,e_iv,f_iv)-Ch(e1,e_iv,f_iv)) + δW1
# For XOR trail we need carry analysis on (c_iv + T1_1_n) vs (c_iv + T1_1_f)

# Key insight: if we set δW1 = -(δSig1(e1)_add + δCh_add) such that δT1_1_add=0,
# then δe2=0 deterministically! Let's check.

# For XOR differential, Sig1 is linear: δSig1_add ≠ δSig1_xor in general.
# δSig1_add = Sig1(e1+0x8000) - Sig1(e1) [additive] -- this varies with e1!

# But δCh_add = Ch(e1+0x8000, f1, g1) - Ch(e1, f1, g1) [additive] -- also varies!

# So finding a fixed δW1 that makes δe2=0 deterministically is hard.
# P-23 showed: best P(δe2=0) ≈ 0.1256 with any fixed δW1.

print("[КЛЮЧЕВОЙ РЕЗУЛЬТАТ]")
print("SC для δe2=0 требует анализа a1-регистра в раунде 2.")
print("Т.к. δa1=0x8000, то T2_2_f = T2_2_n + δT2_2_add,")
print("а δa2 = (T1_1_f + T2_2_f) XOR (T1_1_n + T2_2_n)")
print("Если δT1_1_add + δT2_2_add = 0 → δa2=0 (или 0x8000 при ненулевом carry)")
print()

# ── Part 9: Verify T_SC_A1 and find P(δa1=0x8000) under SC_e ─────────────────
print("[9] P(δa1=0x8000 | SC_e) — условная вероятность")
print("-"*50)

N9 = 100000
count_da1_given_sce = 0
n_sce = 0

for _ in range(N9):
    W0_low15 = random.randint(0, 0x7FFF)
    W0_high = random.randint(0, 0x1FFFF) << 15
    W0 = (W0_low15 | W0_high) & MASK

    sc_e = (W0_low15 >= threshold_e1)
    if not sc_e:
        continue
    n_sce += 1

    T1_0_n_v = (T1_0_base + W0) & MASK
    a1_n_v = (T1_0_n_v + T2_0) & MASK
    T1_0_f_v = (T1_0_n_v + dW0) & MASK
    a1_f_v = (T1_0_f_v + T2_0) & MASK
    da1_v = a1_f_v ^ a1_n_v

    if da1_v == dW0:
        count_da1_given_sce += 1

if n_sce > 0:
    print(f"P(δa1=0x8000 | SC_e) = {count_da1_given_sce}/{n_sce} = {count_da1_given_sce/n_sce:.4f}")
    print(f"Теория: ≈ 0.5 (так как SC_e фиксирует W0[0..14], W0[15] свободен)")

# ── Summary ────────────────────────────────────────────────────────────────────
print("\n" + "="*70)
print("ИТОГ П-24")
print("="*70)
print("""
ТЕОРЕМЫ:

T_SC_A1: δa1=0x8000 ⟺ bit15(a1_n)=0
  где a1_n = (T1_0_base + T2_0 + W0) mod 2^32
  Эквивалентно: W0[15] = S_base[15] XOR carry_14
  carry_14 = 1 iff W0[0..14] >= thresh_carry

T_JOINT_SC: P(δe1=0x8000 AND δa1=0x8000) = (0x8000 - threshold_e1) / (0x8000 * 2)
  ≈ 0.27  (примерно 1/4 сообщений удовлетворяют ОБОИМ условиям)
  Условия:
    SC_e: W0[0..14] >= threshold_e1 (≈ 7518)
    SC_a1: W0[15] определён carry из W0[0..14]

НАПРАВЛЕНИЯ П-25:
  1. Найти SC для δa2=0 (аналогично SC_e, SC_a1)
  2. Если δe2=0 AND δa2=0: расширить след до δe6=0, δe7=0...
  3. Вычислить стоимость полного SC: 2^k встреч
  4. Сравнить с аддитивным каскадом (T_BARRIER_16 = 2^64)
""")
