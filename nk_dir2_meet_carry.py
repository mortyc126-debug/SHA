#!/usr/bin/env python3
"""
NK Direction 2: Wang forward + backward inversion, meet via carry.

NOT classical MITM (prohibited by N6 — schedule couples both halves).
Instead: can we combine FORWARD structure (Wang, δe=0 for 16 rounds)
with BACKWARD structure (schedule kernel, δW[52..63]=0 for 12 rounds)?

The gap: rounds 17-51 (35 rounds of "chaos").

Concrete questions:
Q1: If Wang gives δe[2..16]=0 AND schedule kernel gives δW[52..63]=0,
    what's the minimum δstate in the gap (rounds 17-51)?
Q2: Does the backward kernel constraint help the forward Wang at all?
Q3: Can we build M that satisfies BOTH constraints simultaneously?
"""

import random, math
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

def expand_schedule(M):
    W=list(M)+[0]*(64-len(M))
    for i in range(16,64): W[i]=(ssig1(W[i-2])+W[i-7]+ssig0(W[i-15])+W[i-16])&MASK32
    return W

def sha_full_trace(M):
    W = expand_schedule(M)
    a,b,c,d,e,f,g,h = IV
    states = [(a,b,c,d,e,f,g,h)]
    for r in range(64):
        T1=(h+sig1(e)+ch(e,f,g)+K[r]+W[r])&MASK32
        T2=(sig0(a)+maj(a,b,c))&MASK32
        h,g,f=g,f,e; e=(d+T1)&MASK32; d,c,b=c,b,a; a=(T1+T2)&MASK32
        states.append((a,b,c,d,e,f,g,h))
    H = tuple((states[64][i]+IV[i])&MASK32 for i in range(8))
    return states, W, H

def wang_cascade(W_base, dW0, R):
    W1=list(W_base); DW=[0]*16; DW[0]=dW0
    for step in range(min(R-1,15)):
        wi=step+1
        W2=[(W1[i]+DW[i])&MASK32 for i in range(16)]
        s1,_,_=sha_full_trace(W1); s2,_,_=sha_full_trace(W2)
        # Only need state up to step+2
        de=(s2[step+2][4]-s1[step+2][4])&MASK32
        DW[wi]=(-de)&MASK32
    W2=[(W1[i]+DW[i])&MASK32 for i in range(16)]
    return W1, W2, DW

def state_xor_hw(s1, s2):
    return sum(hw(a^b) for a,b in zip(s1,s2))


# ================================================================
# Q1: Wang + schedule kernel overlap
#
# Wang uses DW[0..15] to build cascade.
# Schedule kernel at R=52 has dim=128 → 128 DOF in DW[0..15]
# that don't affect W[52..63].
#
# Can Wang cascade STAY WITHIN the kernel?
# I.e., are the Wang corrections DW[1..15] compatible with
# the constraint δW[52..63]=0?
# ================================================================

print("=" * 70)
print("Q1: Are Wang corrections compatible with schedule kernel?")
print("=" * 70)

N = 200
random.seed(42)

# For each Wang pair, check: does δW[52..63] = 0?
dw_52_63_hw = []

for trial in range(N):
    M_base = [random.randint(0,MASK32) for _ in range(16)]
    W1, W2, DW = wang_cascade(M_base, 1, 16)

    Wexp1 = expand_schedule(W1)
    Wexp2 = expand_schedule(W2)

    # Check δW[52..63]
    dw_late = 0
    for r in range(52, 64):
        dw_late += hw(Wexp1[r] ^ Wexp2[r])
    dw_52_63_hw.append(dw_late)

avg_hw = sum(dw_52_63_hw) / N
min_hw = min(dw_52_63_hw)
max_hw = max(dw_52_63_hw)

print(f"\n  N={N} Wang pairs (dW0=1)")
print(f"  HW(δW[52..63]): avg={avg_hw:.1f}, min={min_hw}, max={max_hw}")
print(f"  Expected for random: 12 words × 16 bits = 192")

if avg_hw > 150:
    print(f"\n  Wang corrections are NOT in schedule kernel.")
    print(f"  δW[52..63] ≈ random. Wang and kernel are independent.")
else:
    print(f"\n  Some overlap detected!")


# ================================================================
# Q2: What if we CONSTRAIN Wang to the kernel?
#
# The kernel at R=52 is a 128-dim subspace of input space (512 bits).
# Standard Wang uses ALL 512 bits freely.
# Kernel-constrained Wang: only 128 bits free → can only cancel
# δe for 128/32 = 4 rounds (not 15).
#
# Actually: Wang needs 1 word per round (32 bits per δe cancellation).
# With 128 DOF: can cancel 4 rounds of δe.
# But which 4? Early rounds (less carry noise) or late (near the gap)?
# ================================================================

print()
print("=" * 70)
print("Q2: Kernel-constrained Wang cascade")
print("=" * 70)

# First: what does the kernel look like?
# Kernel(52) = {δM : δW[52..63] = 0 in XOR model}
# This is a linear constraint on δM = δW[0..15].

# Find kernel basis by Gaussian elimination
def find_kernel_basis(R_start):
    """Find basis of δW[0..15] such that δW[R_start..63]=0 (XOR model)."""
    # Build constraint matrix: each row = one output bit constraint
    # We want: for output rounds R_start..63, all bits zero
    cols = []
    for in_w in range(16):
        for in_b in range(32):
            M = [0]*16; M[in_w] = 1 << in_b
            W = expand_schedule(M)
            col = []
            for out_r in range(R_start, 64):
                for out_b in range(32):
                    col.append((W[out_r]>>out_b)&1)
            cols.append(col)

    nrows = len(cols[0])
    ncols = 512
    # Transpose for elimination
    mat = [[cols[c][r] for c in range(ncols)] for r in range(nrows)]

    # Gaussian elimination, track pivots
    pivot_cols = []
    pivot_row = 0
    for col in range(ncols):
        found = False
        for row in range(pivot_row, nrows):
            if mat[row][col] == 1:
                mat[pivot_row], mat[row] = mat[row], mat[pivot_row]
                found = True
                break
        if not found:
            continue
        pivot_cols.append(col)
        for row in range(nrows):
            if row != pivot_row and mat[row][col] == 1:
                for c in range(ncols):
                    mat[row][c] ^= mat[pivot_row][c]
        pivot_row += 1

    rank = len(pivot_cols)
    free_cols = [c for c in range(ncols) if c not in pivot_cols]
    return rank, free_cols, len(free_cols)

rank52, free_cols_52, kdim52 = find_kernel_basis(52)
print(f"  Kernel(52): rank={rank52}, kernel_dim={kdim52}")
print(f"  Free columns (first 20): {free_cols_52[:20]}...")

# Which MESSAGE WORDS do free columns correspond to?
free_words = set()
for c in free_cols_52:
    word = c // 32
    free_words.add(word)
print(f"  Free columns span words: {sorted(free_words)}")

# Can we do Wang cascade using ONLY free columns?
# Wang needs to set DW[r] for r=1..15 to cancel δe.
# If word r is NOT in free_words → cannot modify it → cannot cancel δe[r+1].

print(f"\n  Wang needs to modify: W[1..15]")
print(f"  Kernel allows modifying: words {sorted(free_words)}")
wang_words = set(range(1, 16))
usable = wang_words & free_words
blocked = wang_words - free_words
print(f"  Usable for Wang: {sorted(usable)} ({len(usable)} words)")
print(f"  Blocked: {sorted(blocked)} ({len(blocked)} words)")
print(f"  → Can cancel δe at rounds: {sorted(w+1 for w in usable)}")


# ================================================================
# Q3: δstate profile in the gap (rounds 17-51)
# for standard Wang pairs — baseline measurement
# ================================================================

print()
print("=" * 70)
print("Q3: δstate through the gap (rounds 17-51)")
print("=" * 70)

N3 = 300
random.seed(123)

hw_by_round = defaultdict(list)

for trial in range(N3):
    M_base = [random.randint(0,MASK32) for _ in range(16)]
    W1, W2, DW = wang_cascade(M_base, 1, 16)

    s1, _, _ = sha_full_trace(W1)
    s2, _, _ = sha_full_trace(W2)

    for r in range(64):
        hw_state = state_xor_hw(s1[r+1], s2[r+1])
        hw_by_round[r].append(hw_state)

print(f"\n  N={N3} Wang pairs, full 64 rounds")
print(f"\n  {'r':>3} | {'E[HW(δstate)]':>14} | {'min':>4} {'max':>4} | zone")

for r in range(64):
    avg = sum(hw_by_round[r])/N3
    mn = min(hw_by_round[r])
    mx = max(hw_by_round[r])

    if r == 0: zone = "INIT"
    elif r <= 16: zone = "WANG (δe=0)"
    elif r <= 20: zone = "BARRIER"
    elif r <= 51: zone = "GAP"
    elif r <= 63: zone = "LATE"
    else: zone = ""

    # Only print interesting rounds
    if r <= 2 or (r >= 15 and r <= 22) or r % 8 == 0 or r >= 60:
        print(f"  {r:>3} | {avg:>14.1f} | {mn:>4} {mx:>4} | {zone}")


# ================================================================
# Q4: Can backward from r=64 give any constraint?
#
# Backward: if we FIX H (target), we can compute state[64] = H - IV.
# Then invert rounds backward: state[63], state[62], ...
# Each backward round needs W[r] — which is fixed by schedule.
#
# For collision: H(M1) = H(M2) → state[64](M1) + IV = state[64](M2) + IV
# → state[64](M1) = state[64](M2).
#
# If we fix M1 and vary M2: state[64](M2) must equal state[64](M1).
# We can run backward from state[64](M1) using W(M2).
# But W(M2) ≠ W(M1) (different messages!) → backward gives different state[0].
#
# The question: at what round do backward-propagated states diverge?
# If they stay close for many rounds → "backward kernel" is large.
# ================================================================

print()
print("=" * 70)
print("Q4: Backward divergence — how far back does state stay close?")
print("=" * 70)

def backward_round(state_after, W_r, K_r):
    """Invert one SHA-256 round. state_after = (a,b,c,d,e,f,g,h) after round."""
    a,b,c,d,e,f,g,h = state_after
    # Undo shift: a_before = b, b_before = c, c_before = d
    # e_before = f, f_before = g, g_before = h
    a_before = b
    b_before = c
    c_before = d
    e_before = f
    f_before = g
    g_before = h

    # T2 = sig0(a_before) + maj(a_before, b_before, c_before)
    T2 = (sig0(a_before) + maj(a_before, b_before, c_before)) & MASK32
    # a = T1 + T2 → T1 = a - T2
    T1 = (a - T2) & MASK32
    # e = d_before + T1 → d_before = e - T1
    d_before = (e - T1) & MASK32
    # T1 = h_before + sig1(e_before) + ch(e_before,f_before,g_before) + K + W
    # → h_before = T1 - sig1(e_before) - ch(e_before,f_before,g_before) - K - W
    h_before = (T1 - sig1(e_before) - ch(e_before,f_before,g_before) - K_r - W_r) & MASK32

    return (a_before, b_before, c_before, d_before, e_before, f_before, g_before, h_before)

# Verify backward round
random.seed(42)
M_test = [random.randint(0,MASK32) for _ in range(16)]
states_fwd, W_test, H_test = sha_full_trace(M_test)

# Run backward from state[64]
state_back = states_fwd[64]
backward_ok = True
for r in range(63, -1, -1):
    state_back = backward_round(state_back, W_test[r], K[r])
    if state_back != states_fwd[r]:
        backward_ok = False
        print(f"  Backward mismatch at round {r}!")
        break

if backward_ok:
    print(f"  Backward inversion verified: 64 rounds, exact. ✓")

# Now: take M1, fix state[64](M1), apply backward with W(M2).
# Measure divergence.
print(f"\n  Backward divergence for different M2 (from state[64](M1)):")

N4 = 200
random.seed(42)
M1 = [random.randint(0,MASK32) for _ in range(16)]
s1_states, W1_exp, H1 = sha_full_trace(M1)
target_state64 = s1_states[64]

back_hw = defaultdict(list)

for trial in range(N4):
    M2 = [random.randint(0,MASK32) for _ in range(16)]
    W2_exp = expand_schedule(M2)

    # Backward from target_state64 using M2's schedule
    state = target_state64
    for r in range(63, -1, -1):
        state = backward_round(state, W2_exp[r], K[r])
        hw_diff = state_xor_hw(state, s1_states[r])
        back_hw[r].append(hw_diff)

print(f"\n  {'r':>3} | {'E[HW(δstate)]':>14} | {'min':>4} | note")
for r in [63, 62, 61, 60, 56, 52, 48, 40, 32, 16, 0]:
    avg = sum(back_hw[r])/N4
    mn = min(back_hw[r])
    note = ""
    if avg < 10: note = "★ LOW divergence"
    elif avg < 64: note = "partial divergence"
    else: note = "full divergence"

    print(f"  {r:>3} | {avg:>14.1f} | {mn:>4} | {note}")

# ================================================================
# SYNTHESIS
# ================================================================

print()
print("=" * 70)
print("SYNTHESIS: Forward + Backward")
print("=" * 70)

print("""
  FORWARD (Wang):
    Rounds 1-16: δe = 0 (free via Wang cascade)
    Round 17+: δe ≠ 0, HW(δstate) → 128

  BACKWARD (inversion with different schedule):
    Round 64: δstate = 0 (by target)
    Round 63: depends on δW[63] — one different schedule word
    Round 62: depends on δW[62..63]
    ...

  KEY QUESTION: at what round does backward δstate reach 128?
  If it stays low for many rounds → backward "buys" those rounds free.
  Then: forward covers 1-16, backward covers X-63.
  Gap = 17 to X. If X is small → gap is small → cheaper.

  But: W(M2) ≈ random relative to W(M1) → backward divergence
  should be INSTANT (one round = full diffusion, just like forward).
""")
