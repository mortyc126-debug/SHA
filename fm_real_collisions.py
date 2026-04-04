#!/usr/bin/env python3
"""
Flow Mathematics meets real collisions.

Step 1: Find actual collisions on reduced-round SHA-256 (8 rounds)
Step 2: Analyze them in FM language — where △ lives, where it dies,
        what the GPK/carry structure looks like at each round.
Step 3: Look for patterns that scale.
"""

import random
import time
from collections import defaultdict

MASK32 = 0xFFFFFFFF

K = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
    0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
    0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
    0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
    0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
    0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
]

IV = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
]

def rotr(x, n):
    return ((x >> n) | (x << (32 - n))) & MASK32

def sig0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def sig1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def ssig0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def ssig1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)
def ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
def hw(x): return bin(x & MASK32).count('1')

def sha256_round_trace(msg_words, rounds=8):
    """Full trace: state after each round + T1, T2, carry info."""
    W = list(msg_words) + [0] * (64 - len(msg_words))
    for i in range(16, 64):
        W[i] = (ssig1(W[i-2]) + W[i-7] + ssig0(W[i-15]) + W[i-16]) & MASK32

    a, b, c, d, e, f, g, h = IV
    trace = []
    trace.append({'r': -1, 'state': [a, b, c, d, e, f, g, h]})

    for r in range(rounds):
        T1_parts = {
            'h': h,
            'sig1_e': sig1(e),
            'ch_efg': ch(e, f, g),
            'K': K[r],
            'W': W[r],
        }
        T1 = (h + sig1(e) + ch(e, f, g) + K[r] + W[r]) & MASK32
        T2_parts = {
            'sig0_a': sig0(a),
            'maj_abc': maj(a, b, c),
        }
        T2 = (sig0(a) + maj(a, b, c)) & MASK32

        # Carry analysis for e_new = d + T1
        e_new = (d + T1) & MASK32
        carry_e = ((d + T1) >> 32) & 1  # overflow carry
        # Bit-level carry for d + T1
        carry_bits_e = compute_carry_bits(d, T1)

        # Carry for a_new = T1 + T2
        a_new = (T1 + T2) & MASK32
        carry_bits_a = compute_carry_bits(T1, T2)

        h, g, f = g, f, e
        e = e_new
        d, c, b = c, b, a
        a = a_new

        trace.append({
            'r': r,
            'state': [a, b, c, d, e, f, g, h],
            'T1': T1, 'T2': T2,
            'T1_parts': T1_parts,
            'T2_parts': T2_parts,
            'carry_e': carry_bits_e,
            'carry_a': carry_bits_a,
            'W': W[r],
        })

    H = [(trace[-1]['state'][i] + IV[i]) & MASK32 for i in range(8)]
    return trace, H

def compute_carry_bits(a, b):
    """Compute carry bit at each position for a + b."""
    carries = 0
    c = 0
    for k in range(32):
        ak = (a >> k) & 1
        bk = (b >> k) & 1
        c = (ak & bk) | (ak & c) | (bk & c)
        carries |= (c << k)
    return carries

def gpk_string(a, b):
    """Compute GPK string for addition a + b."""
    result = []
    for k in range(32):
        ak = (a >> k) & 1
        bk = (b >> k) & 1
        if ak == 1 and bk == 1:
            result.append('G')
        elif ak == 0 and bk == 0:
            result.append('K')
        else:
            result.append('P')
    return result

def diff_bits(x, y):
    """Return list of bit positions where x and y differ."""
    d = x ^ y
    return [k for k in range(32) if (d >> k) & 1]


# ================================================================
# STEP 1: Find collisions on 8-round SHA-256
# Birthday on 256 bits needs 2^128 — too expensive.
# But: birthday on single OUTPUT WORD needs 2^16.
# Strategy: find pairs where ALL 8 output words match.
# For 8 rounds: reduced diffusion → easier than full SHA-256.
#
# Actually for 8 rounds the output is still 256 bits.
# Birthday on 256 bits: need ~2^128 pairs.
# We can't find full collisions by brute force.
#
# Alternative: use Wang-style adaptive cascade.
# From the methodology: adaptive ΔW gives δe[2..R]=0 for free.
# Then birthday on the remaining constraint.
# ================================================================

print("=" * 70)
print("STEP 1: Finding collision material on reduced SHA-256")
print("=" * 70)

def wang_cascade(W_base, dW0=1, rounds=8):
    """
    Wang-style adaptive cascade.
    Returns (W1, W2) where W2 = W1 + DW with δe[2..rounds]=0.
    DW[0] = dW0, DW[r] adapted for r=1..rounds-2.
    """
    W1 = list(W_base)
    DW = [0] * 16
    DW[0] = dW0

    # Build W2 adaptively
    for step in range(rounds - 1):  # steps 0..rounds-2
        wi = step + 1  # word index to adapt (W[1], W[2], ...)
        if wi >= 16:
            break

        # Current W2
        W2 = [(W1[i] + DW[i]) & MASK32 for i in range(16)]

        # Compute traces up to round step+2
        t1, _ = sha256_round_trace(W1, rounds=step + 2)
        t2, _ = sha256_round_trace(W2, rounds=step + 2)

        # δe at round step+2
        de = (t2[-1]['state'][4] - t1[-1]['state'][4]) & MASK32

        # Adapt DW[wi] to cancel δe
        DW[wi] = (-de) & MASK32

    W2 = [(W1[i] + DW[i]) & MASK32 for i in range(16)]
    return W1, W2, DW

# Find pairs with δe[2..8]=0, then check full hash
print("\nSearching for 8-round near-collisions via Wang cascade + birthday...")

R = 8
found_pairs = []
N_SEARCH = 500000
t0 = time.time()

# Birthday on remaining hash difference
hash_table = {}  # H_tuple -> (W1, W2, DW)

random.seed(42)
for trial in range(N_SEARCH):
    W_base = [random.randint(0, MASK32) for _ in range(16)]
    W1, W2, DW = wang_cascade(W_base, dW0=1, rounds=R)

    # Compute hashes
    _, H1 = sha256_round_trace(W1, rounds=R)
    _, H2 = sha256_round_trace(W2, rounds=R)

    # Check δe traces
    t1, _ = sha256_round_trace(W1, rounds=R)
    t2, _ = sha256_round_trace(W2, rounds=R)

    # Hash difference
    dH = tuple((H1[i] ^ H2[i]) for i in range(8))
    hw_total = sum(hw(d) for d in dH)

    if hw_total == 0:
        found_pairs.append((W1, W2, DW, H1, H2))
        print(f"  ★ FULL COLLISION FOUND at trial {trial}!")
        break

    # Track best
    if trial < 10 or trial % 100000 == 0:
        de_trace = []
        for r_idx in range(1, R + 1):
            de = (t2[r_idx]['state'][4] - t1[r_idx]['state'][4]) & MASK32
            de_trace.append(de)

    if trial == 0:
        print(f"  Trial 0: HW(δH)={hw_total}, δe trace: {['0' if x==0 else f'{hw(x)}' for x in de_trace]}")

elapsed = time.time() - t0
print(f"\n  Searched {N_SEARCH} trials in {elapsed:.1f}s")
print(f"  Full collisions found: {len(found_pairs)}")

# ================================================================
# STEP 2: Analyze near-collisions in FM language
# Even without full collision, analyze the STRUCTURE of δstate
# ================================================================

print()
print("=" * 70)
print("STEP 2: Analyze differential structure in FM language")
print("=" * 70)

# Take a sample pair and trace everything
random.seed(42)
W_base = [random.randint(0, MASK32) for _ in range(16)]

for R in [8, 16]:
    print(f"\n  ─── R = {R} rounds ───")
    W1, W2, DW = wang_cascade(W_base, dW0=1, rounds=R)
    trace1, H1 = sha256_round_trace(W1, rounds=R)
    trace2, H2 = sha256_round_trace(W2, rounds=R)

    print(f"\n  δW: {[f'{d:#010x}' if d != 0 else '0' for d in DW[:R]]}")

    print(f"\n  {'r':>3} | {'HW(δa)':>6} {'HW(δe)':>6} | {'δe':>10} | {'HW(δstate)':>10} | carry_e GPK | shields")

    for r_idx in range(R + 1):
        s1 = trace1[r_idx]['state']
        s2 = trace2[r_idx]['state']

        da = (s2[0] - s1[0]) & MASK32
        de = (s2[4] - s1[4]) & MASK32

        # XOR differences for each register
        diffs = [s1[i] ^ s2[i] for i in range(8)]
        hw_state = sum(hw(d) for d in diffs)

        r = trace1[r_idx]['r']

        # GPK and carry analysis for this round
        gpk_info = ""
        shield_info = ""

        if r_idx > 0 and r_idx <= R:
            t1r = trace1[r_idx]
            t2r = trace2[r_idx]

            # GPK for e_new = d + T1
            d1 = trace1[r_idx - 1]['state'][3] if r_idx > 0 else IV[3]
            d2 = trace2[r_idx - 1]['state'][3] if r_idx > 0 else IV[3]

            # Count G, P, K in carry of this round
            if r_idx > 0:
                gpk1 = gpk_string(d1, t1r['T1'])
                gpk2 = gpk_string(d2, t2r['T1'])

                # How many GPK symbols differ between the two traces?
                gpk_diff = sum(1 for a, b in zip(gpk1, gpk2) if a != b)
                gpk_same = 32 - gpk_diff

                g1 = gpk1.count('G')
                p1 = gpk1.count('P')
                k1 = gpk1.count('K')

                gpk_info = f"G{g1:2d} P{p1:2d} K{k1:2d} Δ{gpk_diff:2d}"

                # Shield analysis: where does △ die?
                # △ dies when both traces have same GPK AND carry resolves identically
                # Shield = position where carry differs in GPK but result is same
                carry1 = t1r['carry_e']
                carry2 = t2r['carry_e']
                carry_diff = carry1 ^ carry2
                shields = hw(carry_diff)  # positions where carry differs
                shield_info = f"carry_diff={shields:2d}"

        r_label = f"IV" if r == -1 else f"{r:2d}"
        de_str = f"0x{de:08x}" if de != 0 else "    0     "
        print(f"  {r_label:>3} | {hw(da):>6} {hw(de):>6} | {de_str} | {hw_state:>10} | {gpk_info:>20} | {shield_info}")

    # Final hash difference
    dH = [H1[i] ^ H2[i] for i in range(8)]
    hw_hash = sum(hw(d) for d in dH)
    print(f"\n  Hash: HW(δH) = {hw_hash}")
    for i in range(8):
        if dH[i] != 0:
            print(f"    H[{i}]: δ = 0x{dH[i]:08x} (HW={hw(dH[i])})")

# ================================================================
# STEP 3: GPK pattern analysis across many pairs
# ================================================================

print()
print("=" * 70)
print("STEP 3: GPK statistics across 1000 Wang pairs (R=8)")
print("=" * 70)

R = 8
N_PAIRS = 1000
random.seed(123)

# Collect: for each round, how many GPK symbols differ between traces?
gpk_diff_by_round = defaultdict(list)
carry_diff_by_round = defaultdict(list)
de_zero_count = [0] * (R + 1)
shield_events = defaultdict(int)  # round -> count of shield events

for _ in range(N_PAIRS):
    W_base = [random.randint(0, MASK32) for _ in range(16)]
    W1, W2, DW = wang_cascade(W_base, dW0=1, rounds=R)
    trace1, H1 = sha256_round_trace(W1, rounds=R)
    trace2, H2 = sha256_round_trace(W2, rounds=R)

    for r_idx in range(1, R + 1):
        t1r = trace1[r_idx]
        t2r = trace2[r_idx]

        de = (t2r['state'][4] - t1r['state'][4]) & MASK32
        if de == 0:
            de_zero_count[r_idx] += 1

        # GPK diff for e_new computation
        d1_prev = trace1[r_idx - 1]['state'][3]
        d2_prev = trace2[r_idx - 1]['state'][3]

        gpk1 = gpk_string(d1_prev, t1r['T1'])
        gpk2 = gpk_string(d2_prev, t2r['T1'])
        gpk_diff = sum(1 for a, b in zip(gpk1, gpk2) if a != b)
        gpk_diff_by_round[r_idx].append(gpk_diff)

        carry1 = t1r['carry_e']
        carry2 = t2r['carry_e']
        carry_diff = hw(carry1 ^ carry2)
        carry_diff_by_round[r_idx].append(carry_diff)

        # Shield: positions where GPK differs but carry result is same
        for k in range(32):
            g1, g2 = gpk1[k], gpk2[k]
            c1 = (carry1 >> k) & 1
            c2 = (carry2 >> k) & 1
            if g1 != g2 and c1 == c2:
                shield_events[r_idx] += 1

print(f"\n  {'Round':>5} | {'δe=0':>5} | {'E[GPK_diff]':>11} | {'E[carry_diff]':>13} | {'shields/pair':>12}")
for r in range(1, R + 1):
    e_gpk = sum(gpk_diff_by_round[r]) / len(gpk_diff_by_round[r])
    e_carry = sum(carry_diff_by_round[r]) / len(carry_diff_by_round[r])
    shields_per = shield_events[r] / N_PAIRS
    print(f"  {r:>5} | {de_zero_count[r]:>5} | {e_gpk:>11.2f} | {e_carry:>13.2f} | {shields_per:>12.2f}")

print(f"\n  Wang cascade: δe=0 from round 2 to {R}: {all(de_zero_count[r] == N_PAIRS for r in range(2, R+1))}")

# ================================================================
# STEP 4: What happens at the BARRIER (round where δe stops being 0)
# ================================================================

print()
print("=" * 70)
print("STEP 4: Barrier analysis — what kills δe=0?")
print("=" * 70)

# For R rounds, δe[2..R]=0. At round R+1 (if we extend), δe ≠ 0.
# Analyze round R+1 in detail.

R_cascade = 8
R_extend = R_cascade + 2  # look at 2 rounds past barrier

random.seed(42)
W_base = [random.randint(0, MASK32) for _ in range(16)]
W1, W2, DW = wang_cascade(W_base, dW0=1, rounds=R_cascade)

# Extend computation to R_extend rounds (DW[R_cascade..15] = 0)
trace1_ext, H1_ext = sha256_round_trace(W1, rounds=R_extend)
trace2_ext, H2_ext = sha256_round_trace(W2, rounds=R_extend)

print(f"\n  Wang cascade built for R={R_cascade}, extending to R={R_extend}")
print(f"\n  Round-by-round at barrier:")

for r_idx in range(R_cascade - 1, R_extend + 1):
    s1 = trace1_ext[r_idx]['state']
    s2 = trace2_ext[r_idx]['state']
    r = trace1_ext[r_idx]['r']

    da = (s2[0] - s1[0]) & MASK32
    de = (s2[4] - s1[4]) & MASK32
    da_xor = s1[0] ^ s2[0]
    de_xor = s1[4] ^ s2[4]

    # Which specific bits differ?
    da_bits = diff_bits(s1[0], s2[0])
    de_bits = diff_bits(s1[4], s2[4])

    print(f"\n  Round {r}:")
    print(f"    δa = 0x{da:08x} (HW={hw(da_xor):2d}) XOR=0x{da_xor:08x}")
    print(f"    δe = 0x{de:08x} (HW={hw(de_xor):2d}) XOR=0x{de_xor:08x}")

    if r_idx > 0 and r_idx <= R_extend:
        t1r = trace1_ext[r_idx]
        t2r = trace2_ext[r_idx]

        # T1 difference
        dT1 = (t2r['T1'] - t1r['T1']) & MASK32
        dT2 = (t2r['T2'] - t1r['T2']) & MASK32
        print(f"    δT1 = 0x{dT1:08x} (HW={hw(t1r['T1'] ^ t2r['T1']):2d})")
        print(f"    δT2 = 0x{dT2:08x} (HW={hw(t1r['T2'] ^ t2r['T2']):2d})")

        # GPK for both e and a computations
        d1_prev = trace1_ext[r_idx - 1]['state'][3]
        d2_prev = trace2_ext[r_idx - 1]['state'][3]

        gpk1_e = gpk_string(d1_prev, t1r['T1'])
        gpk2_e = gpk_string(d2_prev, t2r['T1'])
        gpk_diff = sum(1 for a, b in zip(gpk1_e, gpk2_e) if a != b)

        # Count G/P/K in each
        for label, gpk in [("M1", gpk1_e), ("M2", gpk2_e)]:
            g = gpk.count('G')
            p = gpk.count('P')
            k = gpk.count('K')
            print(f"    GPK(d+T1) {label}: G={g:2d} P={p:2d} K={k:2d}")
        print(f"    GPK diff: {gpk_diff}/32 positions")

        if de == 0 and r >= 2:
            print(f"    ★ δe=0: CASCADE ACTIVE — Wang correction working")
        elif de != 0 and r > R_cascade:
            print(f"    ✗ δe≠0: BARRIER HIT — no more free W[r] to correct")
            print(f"    δe bits: {de_bits[:10]}{'...' if len(de_bits)>10 else ''}")

print()
print("=" * 70)
print("DONE. Next: analyze with SAT-found collisions for deeper structure.")
print("=" * 70)
