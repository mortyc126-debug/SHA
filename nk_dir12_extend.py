#!/usr/bin/env python3
"""
NK Direction 12: Reduced-round collisions → extend.

Plan:
1. Find ACTUAL collision on SHA-256/16 via brute force (birthday on 16 rounds)
2. Analyze the collision pair's δstate[16] structure
3. Ask: can this δstate[16] be "extended" to round 17, 18, ... ?
4. What's the cost of extension per round?

If extension cost < 32 bits/round → cheaper than Wang barrier.
"""

import random, time
from collections import defaultdict

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

def sha_r(M, R):
    """SHA-256 truncated to R rounds, returns (state_after_R, H_R)."""
    W=list(M)+[0]*(64-len(M))
    for i in range(16,64): W[i]=(ssig1(W[i-2])+W[i-7]+ssig0(W[i-15])+W[i-16])&MASK32
    a,b,c,d,e,f,g,h=IV
    for r in range(R):
        T1=(h+sig1(e)+ch(e,f,g)+K[r]+W[r])&MASK32;T2=(sig0(a)+maj(a,b,c))&MASK32
        h,g,f=g,f,e;e=(d+T1)&MASK32;d,c,b=c,b,a;a=(T1+T2)&MASK32
    state = (a,b,c,d,e,f,g,h)
    H = tuple((s+v)&MASK32 for s,v in zip(state, IV))
    return state, H, W

def state_diff_hw(s1, s2):
    return sum(hw(a^b) for a,b in zip(s1,s2))


# ================================================================
# STEP 1: Find collision on SHA-256/R for small R
# Birthday on H[0] (32 bits): need ~2^16 messages.
# Full collision (256 bits): need ~2^128 — too expensive.
# STATE collision (256 bits at round R): need ~2^128 — too expensive.
#
# Instead: find H[0]-collision at round R, then check δstate.
# H[0]-collision: two M with same H[0] after R rounds.
# δH[0] = 0, but δH[1..7] ≠ 0 in general.
# ================================================================

print("=" * 70)
print("STEP 1: Find H[0]-collisions on reduced SHA-256")
print("=" * 70)

for R in [8, 16, 24, 32]:
    random.seed(R * 100 + 42)
    N_SEARCH = 300000  # birthday on 32 bits needs ~2^16 = 65K

    ht = {}  # H[0] → (M, state, H_full)
    collisions = []
    t0 = time.time()

    for trial in range(N_SEARCH):
        M = tuple(random.randint(0, MASK32) for _ in range(16))
        state, H, W = sha_r(M, R)

        key = H[0]
        if key in ht:
            M_prev, state_prev, H_prev = ht[key]
            if M != M_prev:
                collisions.append((M_prev, M, state_prev, state, H_prev, H))
                if len(collisions) >= 5:
                    break
        else:
            ht[key] = (M, state, H)

    elapsed = time.time() - t0
    print(f"\n  R={R}: {len(collisions)} H[0]-collisions in {trial+1} trials ({elapsed:.1f}s)")

    if collisions:
        # Analyze first collision
        M1, M2, s1, s2, H1, H2 = collisions[0]

        print(f"    δH: ", end="")
        for i in range(8):
            d = H1[i] ^ H2[i]
            print(f"H[{i}]={'0' if d==0 else hw(d)}", end=" ")
        print()

        # How many H-words match?
        h_match = sum(1 for i in range(8) if H1[i] == H2[i])
        print(f"    H-words matching: {h_match}/8")

        # δstate analysis
        ds_hw = state_diff_hw(s1, s2)
        print(f"    HW(δstate[{R}]) = {ds_hw}")

        # Per-register
        reg = ['a','b','c','d','e','f','g','h']
        for i in range(8):
            d = s1[i] ^ s2[i]
            if d == 0:
                print(f"    δ{reg[i]}[{R}] = 0 ★")
            else:
                print(f"    δ{reg[i]}[{R}] = 0x{d:08x} (HW={hw(d)})")

        # Now: extend this collision to R+1, R+2, ...
        # Run both messages through MORE rounds and check
        print(f"\n    Extension analysis:")
        for R_ext in range(R+1, min(R+9, 65)):
            _, H1e, _ = sha_r(M1, R_ext)
            _, H2e, _ = sha_r(M2, R_ext)

            h0_match = (H1e[0] == H2e[0])
            h_match_ext = sum(1 for i in range(8) if H1e[i] == H2e[i])
            hw_dH = sum(hw(H1e[i]^H2e[i]) for i in range(8))

            marker = "★ H[0] still matches!" if h0_match else ""
            print(f"      R={R_ext}: H-words matching={h_match_ext}/8, HW(δH)={hw_dH} {marker}")


# ================================================================
# STEP 2: Can we find MULTI-WORD collisions on reduced rounds?
# H[0] ∧ H[1] collision (64 bits): need ~2^32 — expensive but doable
# Check: at what R does birthday on H[0]∧H[1] become feasible?
# ================================================================

print()
print("=" * 70)
print("STEP 2: Multi-word partial collisions")
print("=" * 70)

# Birthday on H[0] only (32 bits)
for R in [8, 16]:
    random.seed(R * 200 + 7)
    N = 300000

    # Count: how many H-words match in the BEST collision?
    ht = {}
    best_match = 0
    best_pair = None

    for trial in range(N):
        M = tuple(random.randint(0, MASK32) for _ in range(16))
        _, H, _ = sha_r(M, R)

        key = H[0]
        if key in ht:
            M_prev, H_prev = ht[key]
            if M != M_prev:
                match = sum(1 for i in range(8) if H[i] == H_prev[i])
                if match > best_match:
                    best_match = match
                    best_pair = (M_prev, M, H_prev, H)
        else:
            ht[key] = (M, H)

    print(f"\n  R={R}: best H[0]-collision has {best_match}/8 matching words")
    if best_pair:
        M1, M2, H1, H2 = best_pair
        for i in range(8):
            d = H1[i] ^ H2[i]
            status = "MATCH ★" if d == 0 else f"HW={hw(d)}"
            print(f"    H[{i}]: {status}")


# ================================================================
# STEP 3: Structure of δstate at collision point
# For H[0]-collisions: what does δstate look like?
# Is it structured (few nonzero registers) or random?
# ================================================================

print()
print("=" * 70)
print("STEP 3: δstate structure at H[0]-collision point")
print("=" * 70)

R = 16
random.seed(42)
N = 300000

ht = {}
ds_collection = []

for trial in range(N):
    M = tuple(random.randint(0, MASK32) for _ in range(16))
    s, H, _ = sha_r(M, R)

    key = H[0]
    if key in ht:
        M_prev, s_prev, H_prev = ht[key]
        if M != M_prev:
            ds = tuple(s[i] ^ s_prev[i] for i in range(8))
            ds_hw = sum(hw(d) for d in ds)
            n_zero = sum(1 for d in ds if d == 0)
            ds_collection.append((ds_hw, n_zero, ds))
    else:
        ht[key] = (M, s, H)

print(f"\n  R={R}, {len(ds_collection)} H[0]-collisions found")

if ds_collection:
    # Distribution of HW(δstate) and zero registers
    avg_hw = sum(x[0] for x in ds_collection) / len(ds_collection)
    avg_zeros = sum(x[1] for x in ds_collection) / len(ds_collection)

    print(f"  E[HW(δstate)] = {avg_hw:.1f} (expect ~128 for random)")
    print(f"  E[zero registers] = {avg_zeros:.2f} (expect ~0 for random)")

    # Distribution of zero register count
    zero_dist = defaultdict(int)
    for _, nz, _ in ds_collection:
        zero_dist[nz] += 1

    print(f"\n  Zero register distribution:")
    for nz in sorted(zero_dist.keys()):
        pct = zero_dist[nz] / len(ds_collection) * 100
        print(f"    {nz} zero regs: {zero_dist[nz]} ({pct:.1f}%)")

    # Best cases (most zero registers)
    ds_collection.sort(key=lambda x: (-x[1], x[0]))
    print(f"\n  Top 5 (most zero registers):")
    for ds_hw, n_zero, ds in ds_collection[:5]:
        reg = ['a','b','c','d','e','f','g','h']
        zeros = [reg[i] for i in range(8) if ds[i] == 0]
        print(f"    {n_zero} zeros ({','.join(zeros)}), HW={ds_hw}")


# ================================================================
# STEP 4: Extension probability
# Given H[0]-collision at round R: P(H[0] still matches at R+1)?
# This is: P(δH[0] = 0 at R+1 | δH[0] = 0 at R).
# If > 2^{-32}: there's CORRELATION → extension is cheaper.
# ================================================================

print()
print("=" * 70)
print("STEP 4: Extension probability (key test)")
print("=" * 70)

for R in [8, 16]:
    random.seed(R * 300 + 13)
    N = 300000

    ht = {}
    n_coll_R = 0
    n_still_at_R1 = 0
    n_still_at_R2 = 0

    for trial in range(N):
        M = tuple(random.randint(0, MASK32) for _ in range(16))
        _, H_R, _ = sha_r(M, R)

        key = H_R[0]
        if key in ht:
            M_prev = ht[key]
            if M != M_prev:
                n_coll_R += 1

                # Check R+1
                _, H_R1, _ = sha_r(M, R+1)
                _, H_R1_prev, _ = sha_r(M_prev, R+1)
                if H_R1[0] == H_R1_prev[0]:
                    n_still_at_R1 += 1

                # Check R+2
                _, H_R2, _ = sha_r(M, R+2)
                _, H_R2_prev, _ = sha_r(M_prev, R+2)
                if H_R2[0] == H_R2_prev[0]:
                    n_still_at_R2 += 1
        else:
            ht[key] = M

    print(f"\n  R={R}: {n_coll_R} H[0]-collisions at R")
    if n_coll_R > 0:
        p_r1 = n_still_at_R1 / n_coll_R
        p_r2 = n_still_at_R2 / n_coll_R
        expected = 1 / 2**32
        ratio_r1 = p_r1 / expected if expected > 0 else 0

        print(f"    Still H[0]-collision at R+1: {n_still_at_R1}/{n_coll_R} = {p_r1:.6f}")
        print(f"    Still H[0]-collision at R+2: {n_still_at_R2}/{n_coll_R} = {p_r2:.6f}")
        print(f"    Expected if independent: {expected:.2e}")
        print(f"    Ratio (R+1): {ratio_r1:.1f}×")

        if p_r1 > 0:
            print(f"    ★ NON-ZERO extension probability!")
            print(f"    Cost of 1-round extension: ~{1/p_r1:.0f} trials (vs 2^32 = {2**32:.0f} random)")
        else:
            print(f"    Extension probability = 0 → no correlation")


# ================================================================
# SYNTHESIS
# ================================================================

print()
print("=" * 70)
print("SYNTHESIS: Direction 12")
print("=" * 70)

print("""
  Reduced-round H[0]-collisions found easily (birthday 2^16).

  But extending from R to R+1:
  P(H[0] still matches at R+1 | matches at R) = ???

  If = 0/N: collision at R and collision at R+1 are INDEPENDENT.
  Extension = new birthday per round = 2^16 per round.
  Total for 64 rounds: 64 × 2^16 = 2^22 (if sequential).
  But need SAME pair for all rounds → back to 2^128.

  If > 0: some correlation → extension cheaper than new birthday.
""")
