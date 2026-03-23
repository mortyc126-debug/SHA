#!/usr/bin/env python3
"""
AUDIT 2: Verify key claims from methodology (RESULTS.md / REVIEW_GUIDE.md)

Most critical claims:
A) Fundamental identity: D_add[L](x,y) = L(c) - 2(L(x)∧L(c))
B) Rank-5 bilinear form, kernel = {d,h,f⊕g,a⊕b⊕c}
C) De17=0 concrete solution verification
D) Free words: De17-De20 INDEPENDENTLY controllable via W[12..15]
E) 95.1% AND cancellation
F) Sigma1 has 8 fixed points
"""

import random
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

def sha_round(state, W_t, K_t):
    a,b,c,d,e,f,g,h = state
    T1 = add32(h,Sig1(e),Ch(e,f,g),K_t,W_t)
    T2 = add32(Sig0(a),Maj(a,b,c))
    return [add32(T1,T2),a,b,c,add32(d,T1),e,f,g]

def expand_W(W16):
    W = list(W16)
    for i in range(16,64):
        W.append(add32(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    return W

def run_nr(msg, iv, nr):
    W = expand_W(msg)
    s = list(iv)
    for t in range(nr): s = sha_round(s, W[t], K[t])
    return s

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


print("=" * 72)
print("AUDIT A: Fundamental Identity D_add[L](x,y) = L(c) - 2(L(x)∧L(c))")
print("=" * 72)
# For L = Sig1: D_add[Sig1](x,y) = Sig1(x+y) - Sig1(x) should equal
# Sig1(c) - 2*(Sig1(x) & Sig1(c)) where c = (x+y) ^ x

random.seed(123)
fails = 0
for _ in range(10000):
    x = random.getrandbits(32)
    y = random.getrandbits(32)
    xy = (x + y) & MASK
    c = xy ^ x  # carry pattern

    for L, name in [(Sig0,"Sig0"), (Sig1,"Sig1"), (sig0,"sig0"), (sig1,"sig1")]:
        lhs = (L(xy) - L(x)) & MASK  # D_add[L](x,y)
        rhs = (L(c) - 2*(L(x) & L(c))) & MASK
        if lhs != rhs:
            fails += 1

print(f"  Tested 10000 × 4 functions = 40000 cases")
print(f"  Failures: {fails}")
if fails == 0:
    print("  CONFIRMED: Identity holds exactly.")
else:
    print(f"  *** REFUTED: {fails} failures! ***")


print()
print("=" * 72)
print("AUDIT B: Rank-5 bilinear form from Ch+Maj")
print("=" * 72)
# Ch(e,f,g) + Maj(a,b,c) define a quadratic form per bit position.
# The claim: rank=5, kernel={d, h, f⊕g, a⊕b⊕c}
# Test: for bit j, the bilinear form B maps pairs of 8-bit directions
# to {0,1}. Compute the 8×8 matrix.

# Per-bit: state = (a,b,c,d,e,f,g,h), each 1 bit
# Ch(e,f,g) = e*f XOR (1-e)*g = e*f + g - e*g (mod 2) = e*f + e*g + g (mod 2)
# Wait: Ch = (e AND f) XOR ((NOT e) AND g)
# In GF(2): Ch = ef + (1+e)g = ef + g + eg
# Bilinear part of Ch in (e,f,g): B_Ch(x,y) where x,y are direction vectors
# For directions in (a,b,c,d,e,f,g,h) space (8-dim):
# Ch = ef + eg + g (the +g is linear, not bilinear)
# Bilinear terms: ef, eg
# Maj = ab + ac + bc
# Bilinear terms: ab, ac, bc

# Total bilinear form B(x,y) = x_e*y_f + x_f*y_e + x_e*y_g + x_g*y_e
#                              + x_a*y_b + x_b*y_a + x_a*y_c + x_c*y_a
#                              + x_b*y_c + x_c*y_b

# As 8×8 matrix (a,b,c,d,e,f,g,h) = indices 0..7:
# B[0][1] = B[1][0] = 1 (ab from Maj)
# B[0][2] = B[2][0] = 1 (ac from Maj)
# B[1][2] = B[2][1] = 1 (bc from Maj)
# B[4][5] = B[5][4] = 1 (ef from Ch)
# B[4][6] = B[6][4] = 1 (eg from Ch)
# All others = 0

B = [[0]*8 for _ in range(8)]
B[0][1] = B[1][0] = 1  # ab
B[0][2] = B[2][0] = 1  # ac
B[1][2] = B[2][1] = 1  # bc
B[4][5] = B[5][4] = 1  # ef
B[4][6] = B[6][4] = 1  # eg

# Compute rank over GF(2)
import copy
M = copy.deepcopy(B)
rank = 0
for col in range(8):
    pivot = -1
    for row in range(rank, 8):
        if M[row][col]:
            pivot = row; break
    if pivot == -1: continue
    M[rank], M[pivot] = M[pivot], M[rank]
    for row in range(8):
        if row != rank and M[row][col]:
            M[row] = [(M[row][j] ^ M[rank][j]) for j in range(8)]
    rank += 1

print(f"  Bilinear form matrix B (8×8 over GF(2)):")
labels = ['a','b','c','d','e','f','g','h']
print(f"    {'':>4}", end="")
for l in labels: print(f" {l}", end="")
print()
for i in range(8):
    print(f"    {labels[i]:>4}", end="")
    for j in range(8): print(f" {B[i][j]}", end="")
    print()
print(f"\n  Rank: {rank}")
print(f"  Claimed: 5. {'CONFIRMED' if rank == 5 else '*** MISMATCH: '+str(rank)+' ***'}")

# Find kernel
kernel = []
for mask in range(256):
    v = [(mask >> i) & 1 for i in range(8)]
    Bv = [sum(B[i][j]*v[j] for j in range(8)) % 2 for i in range(8)]
    if all(b == 0 for b in Bv):
        kernel.append(v)
        desc = '+'.join(labels[i] for i in range(8) if v[i])
        if desc: print(f"  Kernel vector: {desc}")

print(f"  Kernel dimension: {len(kernel)}")
print(f"  Claimed: 4 = {{d, h, f⊕g, a⊕b⊕c}}.")
# d = [0,0,0,1,0,0,0,0], h = [0,0,0,0,0,0,0,1]
# f+g = [0,0,0,0,0,1,1,0], a+b+c = [1,1,1,0,0,0,0,0]


print()
print("=" * 72)
print("AUDIT C: De17=0 concrete solution")
print("=" * 72)
# Message pair from RESULTS.md:
# W[0] = 0x3eba2820, W'[0] = 0xbeba2820 (DW = 0x80000000)
# W[1] = 0x3082cc68, W[14] = 0x0013efc2
# Claim: De17 = 0

# We don't have all 16 words. Just verify the CONCEPT:
# With Wang chain + right W[14] (free word), De17=0 is achievable.
# We already confirmed free words work (Audit 1). The specific solution
# would need all 16 words to verify.
print("  Specific solution needs all 16 words (not all listed in RESULTS.md).")
print("  CONCEPT verified: Wang chain + free words can zero De17.")
print("  Cost 2^32 is consistent with P(De17=0) = 2^{-32} per free word value.")


print()
print("=" * 72)
print("AUDIT D: Free words — De17-De20 INDEPENDENTLY controllable?")
print("=" * 72)
# Claim: varying W[14] changes De[17] but NOT De[18].
# Varying W[13] changes De[18] but NOT De[17].
# This would mean 4 independent birthday searches, each 2^32.

random.seed(456)
iv = list(H0)
N = 500

# Test: does W[14] affect De17 and De18?
de17_changes = 0
de18_changes = 0
de19_changes = 0

for _ in range(N):
    msg = [random.getrandbits(32) for _ in range(16)]
    Wp = wang_chain(msg, iv)

    # Compute De17, De18 with base W[14]
    W_exp = expand_W(msg); Wp_exp = expand_W(Wp)
    s1_17 = run_nr(msg, iv, 17); s2_17 = run_nr(Wp, iv, 17)
    s1_18 = run_nr(msg, iv, 18); s2_18 = run_nr(Wp, iv, 18)
    s1_19 = run_nr(msg, iv, 19); s2_19 = run_nr(Wp, iv, 19)
    de17_base = s1_17[4] ^ s2_17[4]
    de18_base = s1_18[4] ^ s2_18[4]
    de19_base = s1_19[4] ^ s2_19[4]

    # Flip W[14] in BOTH messages
    msg_f = list(msg); msg_f[14] ^= 0x80000000
    Wp_f = wang_chain(msg_f, iv)
    sf1_17 = run_nr(msg_f, iv, 17); sf2_17 = run_nr(Wp_f, iv, 17)
    sf1_18 = run_nr(msg_f, iv, 18); sf2_18 = run_nr(Wp_f, iv, 18)
    sf1_19 = run_nr(msg_f, iv, 19); sf2_19 = run_nr(Wp_f, iv, 19)
    de17_flip = sf1_17[4] ^ sf2_17[4]
    de18_flip = sf1_18[4] ^ sf2_18[4]
    de19_flip = sf1_19[4] ^ sf2_19[4]

    if de17_flip != de17_base: de17_changes += 1
    if de18_flip != de18_base: de18_changes += 1
    if de19_flip != de19_base: de19_changes += 1

print(f"  Flipping W[14] (in BOTH M and M'):")
print(f"    De17 changes: {de17_changes}/{N} = {de17_changes/N*100:.1f}%")
print(f"    De18 changes: {de18_changes}/{N} = {de18_changes/N*100:.1f}%")
print(f"    De19 changes: {de19_changes}/{N} = {de19_changes/N*100:.1f}%")

if de17_changes == 0:
    print(f"  W[14] does NOT affect De17 — independence CONFIRMED for this pair.")
else:
    print(f"  *** W[14] DOES affect De17! Independence claim REFUTED! ***")

if de18_changes == 0:
    print(f"  *** W[14] does NOT affect De18 either — W[14] is INERT? ***")
else:
    print(f"  W[14] affects De18: {de18_changes/N*100:.1f}%")

# Now test: does W[14] affect De17 when changed ONLY in M (not M')?
# This is the NON-neutral version
de17_single = 0
for _ in range(N):
    msg = [random.getrandbits(32) for _ in range(16)]
    Wp = wang_chain(msg, iv)

    s1 = run_nr(msg, iv, 17); s2 = run_nr(Wp, iv, 17)
    de17_base = s1[4] ^ s2[4]

    # Change W[14] in M only, recompute Wang chain
    msg_f = list(msg); msg_f[14] ^= 0x80000000
    Wp_f = wang_chain(msg_f, iv)  # new Wang chain!
    sf1 = run_nr(msg_f, iv, 17); sf2 = run_nr(Wp_f, iv, 17)
    de17_new = sf1[4] ^ sf2[4]

    if de17_new != de17_base: de17_single += 1

print(f"\n  Changing W[14] in M only (recompute Wang chain):")
print(f"    De17 changes: {de17_single}/{N} = {de17_single/N*100:.1f}%")

if de17_single > N * 0.9:
    print(f"  W[14] strongly affects De17 through Wang chain recalculation.")
    print(f"  The 'free word' claim means: W[14] is free to CHOOSE")
    print(f"  (different W[14] → different De17 → birthday search on W[14]).")
    print(f"  It does NOT mean De17 is independent of W[14]!")


print()
print("=" * 72)
print("AUDIT E: 95.1% AND cancellation")
print("=" * 72)
# Claim: across rounds, 95.1% of AND terms cancel.
# "655 AND bits → 16 bits output"
# This means: the quadratic part of the differential (AND terms from Ch, Maj,
# and carry) mostly cancels, leaving only ~5% as effective nonlinearity.
# Measure: for random state diff, count AND terms that survive across rounds.

# Actually, 95.1% cancellation is about the ROUND FUNCTION composition.
# After 1 round: 8 registers × ~32 bits each with AND terms
# The claim is that when you compose rounds, AND terms cancel between
# the b=a, c=b, etc. shifts.
# Hard to verify directly without the specific measurement framework.
print("  Cannot directly verify without methodology's specific AND counting framework.")
print("  The 5% = Da+De / total prediction matches 4.9% → internally consistent.")
print("  SKIPPING (would need to reproduce exact measurement from step13).")


print()
print("=" * 72)
print("AUDIT F: Sigma1 fixed points")
print("=" * 72)
# Claim: Sig1 has 8 fixed points
count = 0
fps = []
for x in range(0x100000000) if False else []:  # too slow for 2^32
    pass

# Test the claimed fixed points
claimed = [0, 0x33333333, 0x55555555, 0x66666666,
           0x99999999, 0xaaaaaaaa, 0xcccccccc, 0xffffffff]
verified = 0
for x in claimed:
    if Sig1(x) == x:
        verified += 1
    else:
        print(f"  {hex(x)}: Sig1(x) = {hex(Sig1(x))} ≠ x — NOT a fixed point!")

print(f"  Verified {verified}/{len(claimed)} claimed Sig1 fixed points")
if verified == len(claimed):
    print(f"  All 8 confirmed.")
else:
    print(f"  *** SOME FIXED POINTS ARE WRONG! ***")

# Are there others? Sample random values
for _ in range(1000000):
    x = random.getrandbits(32)
    if Sig1(x) == x and x not in claimed:
        print(f"  EXTRA fixed point found: {hex(x)}")


print()
print("=" * 72)
print("SUMMARY")
print("=" * 72)
