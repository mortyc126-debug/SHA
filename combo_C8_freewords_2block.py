#!/usr/bin/env python3
"""
Combo C8: Free Words × 2-Block Messages

Question: Does a 2-block message give more freedom for differential control?
Block 1: Wang chain to round 20 → DH1 (state diff after 64 rounds + feed-forward)
Block 2: starts with DH1 as initial differential, 16 new message words

Key test: How damaged is DH1? Can Wang chain work in block 2?
"""

import random
MASK = 0xFFFFFFFF

def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK
def Sig0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sig1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def sig0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sig1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)
def Ch(e, f, g): return ((e & f) ^ (~e & g)) & MASK
def Maj(a, b, c): return ((a & b) ^ (a & c) ^ (b & c)) & MASK
def add(*args):
    s = 0
    for a in args: s = (s + a) & MASK
    return s
def hw(x): return bin(x).count('1')

IV = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]
K = [0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
     0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
     0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
     0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
     0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
     0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
     0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
     0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2]

def expand_schedule(W16):
    W = list(W16)
    for t in range(16, 64):
        W.append(add(sig1(W[t-2]), W[t-7], sig0(W[t-15]), W[t-16]))
    return W

def compress(state, W16):
    W = expand_schedule(W16)
    a,b,c,d,e,f,g,h = state
    for t in range(64):
        T1 = add(h, Sig1(e), Ch(e,f,g), K[t], W[t])
        T2 = add(Sig0(a), Maj(a,b,c))
        h,g,f,e,d,c,b,a = g,f,e,add(d,T1),c,b,a,add(T1,T2)
    return [add(state[i], [a,b,c,d,e,f,g,h][i]) for i in range(8)]

def sha256_round_diff(state, state_p, W, Wp, rounds):
    """Run rounds, return state diffs at each round."""
    a,b,c,d,e,f,g,h = state
    ap,bp,cp,dp,ep,fp,gp,hp = state_p
    diffs = []
    for t in range(rounds):
        T1 = add(h, Sig1(e), Ch(e,f,g), K[t], W[t])
        T2 = add(Sig0(a), Maj(a,b,c))
        T1p = add(hp, Sig1(ep), Ch(ep,fp,gp), K[t], Wp[t])
        T2p = add(Sig0(ap), Maj(ap,bp,cp))
        h,g,f,e,d,c,b,a = g,f,e,add(d,T1),c,b,a,add(T1,T2)
        hp,gp,fp,ep,dp,cp,bp,ap = gp,fp,ep,add(dp,T1p),cp,bp,ap,add(T1p,T2p)
        De = e ^ ep
        Da = a ^ ap
        diffs.append((Da, De, [ai^bi for ai,bi in zip([a,b,c,d,e,f,g,h],[ap,bp,cp,dp,ep,fp,gp,hp])]))
    return diffs, [a,b,c,d,e,f,g,h], [ap,bp,cp,dp,ep,fp,gp,hp]

random.seed(0xC8B0)
print("=" * 72)
print("COMBO C8: Free Words × 2-Block Messages")
print("=" * 72)

# ============================================================
# PART 1: How damaged is DH1 after block 1?
# ============================================================
print("\nPART 1: State differential after block 1 (64 rounds)")
print("-" * 50)

N = 500
dh1_hws = []
dh1_per_reg = [[] for _ in range(8)]

for trial in range(N):
    W = [random.randint(0, MASK) for _ in range(16)]
    Wp = list(W)
    Wp[0] = (Wp[0] ^ 0x80000000) & MASK  # DW[0] = MSB flip

    H1 = compress(IV, W)
    H1p = compress(IV, Wp)

    total_hw = sum(hw(H1[i] ^ H1p[i]) for i in range(8))
    dh1_hws.append(total_hw)
    for i in range(8):
        dh1_per_reg[i].append(hw(H1[i] ^ H1p[i]))

avg_hw = sum(dh1_hws) / N
reg_names = ['a','b','c','d','e','f','g','h']
print(f"  Mean HW(DH1): {avg_hw:.1f} / 256 (random = 128)")
print(f"  Min HW(DH1):  {min(dh1_hws)}")
print(f"  Max HW(DH1):  {max(dh1_hws)}")
print(f"  Per register:")
for i in range(8):
    avg = sum(dh1_per_reg[i]) / N
    mn = min(dh1_per_reg[i])
    print(f"    {reg_names[i]}: mean={avg:.1f}/32, min={mn}")

# ============================================================
# PART 2: Can Wang chain work in block 2 with DH1 ≠ 0?
# ============================================================
print(f"\nPART 2: Wang chain feasibility in block 2")
print("-" * 50)

wang_rounds_b2 = []

for trial in range(min(N, 200)):
    W = [random.randint(0, MASK) for _ in range(16)]
    Wp = list(W)
    Wp[0] = (Wp[0] ^ 0x80000000) & MASK

    H1 = compress(IV, W)
    H1p = compress(IV, Wp)

    # Block 2: random message
    W2 = [random.randint(0, MASK) for _ in range(16)]
    W2_sched = expand_schedule(W2)

    # Try Wang chain in block 2: choose DW2[t] to force De=0
    # Start with DH1 as initial state differential
    state2 = list(H1)
    state2p = list(H1p)

    a,b,c,d,e,f,g,h = state2
    ap,bp,cp,dp,ep,fp,gp,hp = state2p

    rounds_ok = 0
    W2p_sched = list(W2_sched)

    for t in range(64):
        # Compute T1 components WITHOUT W contribution
        T1_noW = add(h, Sig1(e), Ch(e,f,g), K[t])
        T1p_noW = add(hp, Sig1(ep), Ch(ep,fp,gp), K[t])
        T2 = add(Sig0(a), Maj(a,b,c))
        T2p = add(Sig0(ap), Maj(ap,bp,cp))

        # For De' = 0: need d + T1 = dp + T1p (mod 2^32)
        # T1 = T1_noW + W[t], T1p = T1p_noW + W'[t]
        # d + T1_noW + W[t] = dp + T1p_noW + W'[t]
        # W'[t] = W[t] + (d + T1_noW) - (dp + T1p_noW)

        target_e = add(d, T1_noW, W2_sched[t])  # desired e' = ep'
        # W'[t] such that dp + T1p_noW + W'[t] = target_e
        needed_Wp = (target_e - add(dp, T1p_noW)) & MASK

        if t < 16:
            W2p_sched[t] = needed_Wp  # Free to choose in block 2
        # For t >= 16, W'[t] is determined by schedule — check if it matches

        # Compute actual round
        T1 = add(T1_noW, W2_sched[t])
        T1p = add(T1p_noW, W2p_sched[t] if t < 16 else W2_sched[t])

        new_e = add(d, T1)
        new_ep = add(dp, T1p)
        new_a = add(T1, T2)
        new_ap = add(T1p, T2p)

        De = new_e ^ new_ep

        if De == 0:
            rounds_ok += 1
        else:
            break

        h,g,f,e,d,c,b,a = g,f,e,new_e,c,b,a,new_a
        hp,gp,fp,ep,dp,cp,bp,ap = gp,fp,ep,new_ep,cp,bp,ap,new_ap

    wang_rounds_b2.append(rounds_ok)

avg_rounds = sum(wang_rounds_b2) / len(wang_rounds_b2)
max_rounds = max(wang_rounds_b2)
print(f"  Wang chain in block 2 (with DH1 ≠ 0):")
print(f"  Mean rounds De=0: {avg_rounds:.1f}")
print(f"  Max rounds De=0:  {max_rounds}")
print(f"  Distribution:")
from collections import Counter
dist = Counter(wang_rounds_b2)
for r in sorted(dist.keys()):
    print(f"    {r} rounds: {dist[r]} ({100*dist[r]/len(wang_rounds_b2):.1f}%)")

# ============================================================
# PART 3: Free words in block 2
# ============================================================
print(f"\nPART 3: Free word analysis in block 2")
print("-" * 50)
print("  Block 2 free words are ALWAYS W'[0..15] (all 16 words)")
print("  because we can choose DW'[t] freely for t=0..15.")
print("  This is identical to block 1: 16 words × 32 bits = 512 bits freedom.")
print()
print("  BUT: the initial state DH1 ≠ 0 means:")
print("  - Da[0], De[0] etc. are NON-ZERO entering block 2")
print("  - Wang chain must FIRST cancel these initial diffs")
print("  - This costs some of the 16 free words")

# ============================================================
# PART 4: Compare total freedom
# ============================================================
print(f"\nPART 4: Total freedom comparison")
print("-" * 50)

# Block 2 with DH1=0 (ideal case)
wang_b2_zeroinit = []
for trial in range(200):
    W2 = [random.randint(0, MASK) for _ in range(16)]
    W2_sched = expand_schedule(W2)
    # Same IV for both → DH1 = 0
    state2 = list(IV)
    state2p = list(IV)

    a,b,c,d,e,f,g,h = state2
    ap,bp,cp,dp,ep,fp,gp,hp = state2p
    # No diff → Wang chain trivially holds forever (De always 0 with DW=0)
    # Not interesting.
    wang_b2_zeroinit.append(64)

# The real question: 2-block where block1 carries the differential
# and block2 tries to cancel it
print(f"  1-block attack:")
print(f"    Free words: W[12-15] = 128 bits")
print(f"    Wang chain: rounds 1-16 free, 17-20 via free words")
print(f"    Barrier: round 21+ (44 rounds, no freedom)")
print()
print(f"  2-block attack:")
print(f"    Block 1: standard (rounds 1-20 controlled)")
print(f"    After block 1: DH1 ≈ {avg_hw:.0f}/256 HW (random)")
print(f"    Block 2: 16 new free words, BUT start with DH1≠0")
print(f"    Wang chain in block 2: {avg_rounds:.1f} rounds mean, {max_rounds} max")
print()

# ============================================================
# PART 5: What if block 1 does NOT carry differential?
# ============================================================
print(f"\nPART 5: Alternative — differential only in block 2")
print("-" * 50)
print("  If DW only in block 2: DH1 = 0 (blocks identical)")
print("  Block 2 = fresh single-block attack with standard IV")
print("  → Equivalent to 1-block attack. No advantage.")
print()

# ============================================================
# PART 6: Partial differential from block 1
# ============================================================
print(f"\nPART 6: Low-weight DH1 search")
print("-" * 50)

# Among our trials, find messages where DH1 has lowest HW
best_idx = min(range(N), key=lambda i: dh1_hws[i])
print(f"  Best DH1 found: HW = {dh1_hws[best_idx]}")
print(f"  Per register at best:")
for i in range(8):
    print(f"    {reg_names[i]}: HW = {dh1_per_reg[i][best_idx]}")

# How many have DH1 with any register at 0?
zero_reg_count = 0
for trial in range(N):
    if any(dh1_per_reg[i][trial] == 0 for i in range(8)):
        zero_reg_count += 1
print(f"\n  Trials with any DH1 register = 0: {zero_reg_count}/{N}")
print(f"  Expected (per register ~2^-32): ~0")

# ============================================================
# VERDICT
# ============================================================
print()
print("=" * 72)
print("VERDICT")
print("=" * 72)

if avg_rounds >= 10:
    verdict = "ALIVE"
    print(f"  {verdict}: Wang chain extends {avg_rounds:.1f} rounds in block 2!")
elif max_rounds >= 16:
    verdict = "ANOMALY"
    print(f"  {verdict}: Some block 2 Wang chains reach {max_rounds} rounds")
elif avg_hw < 100:
    verdict = "ANOMALY"
    print(f"  {verdict}: DH1 below random ({avg_hw:.1f}/256)")
else:
    verdict = "DEAD"
    print(f"  {verdict}: DH1 is random ({avg_hw:.1f}/256), Wang chain in block 2")
    print(f"  averages only {avg_rounds:.1f} rounds. 2-block provides no advantage")
    print(f"  over 1-block attack.")

print()
print("KEY INSIGHT:")
if avg_hw > 120:
    print("  After 64 rounds, the differential from DW[0]=0x80000000 is")
    print("  fully diffused. Block 2 starts with a RANDOM-looking state")
    print("  differential. Wang chain cannot effectively cancel a random")
    print("  256-bit state diff using only 16 message words.")
    print()
    print("  The feed-forward H1[i] = state64[i] + IV[i] preserves the")
    print("  full differential from block 1. No cancellation occurs.")
