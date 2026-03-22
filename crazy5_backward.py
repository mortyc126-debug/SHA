#!/usr/bin/env python3
"""
CRAZY-5: Meet-in-the-Middle Backward Analysis

Key idea: Work BACKWARD from the collision condition.
For a collision, H(M) = H(M'), meaning final state differences must cancel
when added to IV.

Analysis:
1. Forward (Wang chain): De=0 through rounds 1-16, costs ~2^32/round after that
2. Backward: start from De[64]=Da[64]=0 (collision), go backward —
   how far can we maintain low differentials?
3. The "gap" between forward and backward chains gives true complexity.
"""

import random
import numpy as np
import time

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
def sub(a, b): return (a - b) & MASK
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

def sha256_compress(W16, init_state=None):
    """Full 64-round compression, return final state."""
    W = expand_schedule(W16)
    state = list(init_state) if init_state else list(IV)
    a, b, c, d, e, f, g, h = state
    for t in range(64):
        T1 = add(h, Sig1(e), Ch(e, f, g), K[t], W[t])
        T2 = add(Sig0(a), Maj(a, b, c))
        h, g, f, e, d, c, b, a = g, f, e, add(d, T1), c, b, a, add(T1, T2)
    return [a, b, c, d, e, f, g, h]

def sha256_rounds(W16, R, init_state=None):
    """Run R rounds, return state after round R."""
    W = expand_schedule(W16)
    state = list(init_state) if init_state else list(IV)
    a, b, c, d, e, f, g, h = state
    for t in range(R):
        T1 = add(h, Sig1(e), Ch(e, f, g), K[t], W[t])
        T2 = add(Sig0(a), Maj(a, b, c))
        h, g, f, e, d, c, b, a = g, f, e, add(d, T1), c, b, a, add(T1, T2)
    return [a, b, c, d, e, f, g, h]

def sha256_rounds_from(state, W_expanded, start_round, end_round):
    """Run rounds start_round..end_round-1 given expanded W."""
    a, b, c, d, e, f, g, h = state
    for t in range(start_round, end_round):
        T1 = add(h, Sig1(e), Ch(e, f, g), K[t], W_expanded[t])
        T2 = add(Sig0(a), Maj(a, b, c))
        h, g, f, e, d, c, b, a = g, f, e, add(d, T1), c, b, a, add(T1, T2)
    return [a, b, c, d, e, f, g, h]

def inverse_round(state_after, W_t, K_t):
    """
    Invert one SHA-256 round.

    After round t, state is: [new_a, a_prev, b_prev, c_prev, new_e, e_prev, f_prev, g_prev]
    where new_a = T1 + T2, new_e = d_prev + T1
    and T1 = h_prev + Sig1(e_prev) + Ch(e_prev,f_prev,g_prev) + K[t] + W[t]
        T2 = Sig0(a_prev) + Maj(a_prev,b_prev,c_prev)

    Returns state before round t: [a_prev, b_prev, c_prev, d_prev, e_prev, f_prev, g_prev, h_prev]
    """
    new_a = state_after[0]
    a_prev = state_after[1]
    b_prev = state_after[2]
    c_prev = state_after[3]
    new_e = state_after[4]
    e_prev = state_after[5]
    f_prev = state_after[6]
    g_prev = state_after[7]

    # d_prev: new_e = d_prev + T1, and T1 = new_e - d_prev... but we need T1 first
    # T2 is computable: T2 = Sig0(a_prev) + Maj(a_prev, b_prev, c_prev)
    T2 = add(Sig0(a_prev), Maj(a_prev, b_prev, c_prev))
    # new_a = T1 + T2, so T1 = new_a - T2
    T1 = sub(new_a, T2)
    # d_prev = new_e - T1
    d_prev = sub(new_e, T1)
    # h_prev = T1 - Sig1(e_prev) - Ch(e_prev,f_prev,g_prev) - K_t - W_t
    h_prev = sub(sub(sub(sub(T1, Sig1(e_prev)), Ch(e_prev, f_prev, g_prev)), K_t), W_t)

    return [a_prev, b_prev, c_prev, d_prev, e_prev, f_prev, g_prev, h_prev]

def verify_inverse():
    """Verify that inverse_round correctly inverts forward_round."""
    rng = random.Random(12345)
    for _ in range(100):
        W16 = [rng.randint(0, MASK) for _ in range(16)]
        W = expand_schedule(W16)
        state = list(IV)
        a, b, c, d, e, f, g, h = state
        for t in range(10):
            T1 = add(h, Sig1(e), Ch(e, f, g), K[t], W[t])
            T2 = add(Sig0(a), Maj(a, b, c))
            h, g, f, e, d, c, b, a = g, f, e, add(d, T1), c, b, a, add(T1, T2)
        state_after = [a, b, c, d, e, f, g, h]
        recovered = inverse_round(state_after, W[9], K[9])
        # Re-run 9 rounds to get expected state
        state9 = sha256_rounds(W16, 9)
        assert recovered == state9, f"Inverse mismatch at verification"
    print("  Inverse round verification: PASSED (100 tests)")

# ============================================================================
random.seed(0xBEEF)
np.random.seed(0xBEEF)
start_time = time.time()

print("=" * 72)
print("CRAZY-5: Meet-in-the-Middle Backward Analysis")
print("=" * 72)

# Verify inverse round is correct
verify_inverse()

NUM_MESSAGES = 20

def gen_messages():
    msgs = []
    for _ in range(NUM_MESSAGES):
        W16 = [random.randint(0, MASK) for _ in range(16)]
        msgs.append(W16)
    return msgs

base_messages = gen_messages()

# ============================================================================
# PART 1: Forward Chain — Wang-style De=0 baseline
# ============================================================================
print("\n" + "=" * 72)
print("PART 1: Forward Chain — Wang-style De=0 baseline")
print("=" * 72)
print("For each message, apply DW[0] = single-bit diff, measure De[t] through rounds")
print(f"Using {NUM_MESSAGES} base messages\n")

forward_max_de0 = []
for mi, W16 in enumerate(base_messages):
    W = expand_schedule(W16)
    W16p = list(W16)
    W16p[0] ^= (1 << 15)  # single bit diff in W[0]
    Wp = expand_schedule(W16p)

    state = list(IV)
    statep = list(IV)
    a, b, c, d, e, f, g, h = state
    ap, bp, cp, dp, ep, fp, gp, hp = statep

    last_de0 = 0
    de_trace = []
    da_trace = []
    for t in range(64):
        T1 = add(h, Sig1(e), Ch(e, f, g), K[t], W[t])
        T2 = add(Sig0(a), Maj(a, b, c))
        h, g, f, e, d, c, b, a = g, f, e, add(d, T1), c, b, a, add(T1, T2)

        T1p = add(hp, Sig1(ep), Ch(ep, fp, gp), K[t], Wp[t])
        T2p = add(Sig0(ap), Maj(ap, bp, cp))
        hp, gp, fp, ep, dp, cp, bp, ap = gp, fp, ep, add(dp, T1p), cp, bp, ap, add(T1p, T2p)

        de = hw(e ^ ep)
        da = hw(a ^ ap)
        de_trace.append(de)
        da_trace.append(da)
        if de == 0:
            last_de0 = t + 1

    forward_max_de0.append(last_de0)
    if mi < 3:
        print(f"  Msg {mi}: De=0 through round {last_de0}")
        print(f"    De trace (r1-20): {de_trace[:20]}")
        print(f"    Da trace (r1-20): {da_trace[:20]}")

avg_fwd = np.mean(forward_max_de0)
print(f"\n  Average last De=0 round (forward): {avg_fwd:.1f}")
print(f"  Min: {min(forward_max_de0)}, Max: {max(forward_max_de0)}")
print(f"  Note: with Wang-style CHOSEN DW, this extends to round ~16.")
print(f"  Here we use a simple single-bit DW[0] (non-optimized).")

# ============================================================================
# PART 2: Backward Chain — from DH=0 going backward
# ============================================================================
print("\n" + "=" * 72)
print("PART 2: Backward Chain — from De[64]=Da[64]=0 going backward")
print("=" * 72)
print("Start from final state (round 64) with zero differential.")
print("Go backward choosing DW[t] to minimize HW(De) and HW(Da).\n")

backward_results = []

for mi, W16 in enumerate(base_messages):
    if time.time() - start_time > 420:  # 7 min guard
        print(f"  [Time limit approaching, stopping at message {mi}]")
        break

    W = expand_schedule(W16)

    # Compute all intermediate states going forward
    states = [list(IV)]
    a, b, c, d, e, f, g, h = list(IV)
    for t in range(64):
        T1 = add(h, Sig1(e), Ch(e, f, g), K[t], W[t])
        T2 = add(Sig0(a), Maj(a, b, c))
        h, g, f, e, d, c, b, a = g, f, e, add(d, T1), c, b, a, add(T1, T2)
        states.append([a, b, c, d, e, f, g, h])

    # Backward chain: start from final state = states[64]
    # The "primed" path will have states'[64] = states[64] (zero diff at end)
    # Going backward: at each round t (from 63 down to 0), we choose DW[t]
    # to try to keep De and Da small.
    #
    # State after round t: [a,b,c,d,e,f,g,h]
    # Inverse gives state before round t given W[t].
    # For the primed path, W'[t] = W[t] + DW[t].
    # We try different DW[t] values to minimize differential.

    # Start: state'[64] = states[64] (zero diff)
    state_prime = list(states[64])
    DW_back = [0] * 64

    de_backward = [0]  # De at round 64 = 0
    da_backward = [0]  # Da at round 64 = 0
    hw_total_backward = [0]

    best_backward_low_de = 0  # how many rounds back we keep De <= threshold

    for t in range(63, -1, -1):
        # Invert round t for the original path: should give states[t]
        # (This is just verification; we already have states[t])

        # For the primed path, we try DW[t] values.
        # inverse_round needs state_after, W_t, K_t
        # state_prime is the state AFTER round t for the primed path.
        # We want to choose W'[t] = W[t] + DW[t] s.t. the state BEFORE round t
        # (for primed path) is close to states[t].

        # Try DW[t] = 0 first (no message diff at this round)
        best_de = 999
        best_da = 999
        best_total = 999
        best_dw = 0
        best_prev_prime = None

        # Strategy: try DW=0 and several single-bit DW values
        candidates = [0]
        # Add single-bit candidates
        for bit in range(32):
            candidates.append(1 << bit)

        for dw in candidates:
            Wp_t = (W[t] + dw) & MASK
            prev_prime = inverse_round(state_prime, Wp_t, K[t])
            de = hw(prev_prime[4] ^ states[t][4])  # e register
            da = hw(prev_prime[0] ^ states[t][0])  # a register
            total = sum(hw(prev_prime[i] ^ states[t][i]) for i in range(8))

            if total < best_total or (total == best_total and de < best_de):
                best_de = de
                best_da = da
                best_total = total
                best_dw = dw
                best_prev_prime = prev_prime

        DW_back[t] = best_dw
        state_prime = best_prev_prime
        de_backward.append(best_de)
        da_backward.append(best_da)
        hw_total_backward.append(best_total)

    # Reverse so index 0 = round 0
    de_backward = de_backward[::-1]
    da_backward = da_backward[::-1]
    hw_total_backward = hw_total_backward[::-1]

    # Find how far back from round 64 we maintain low De
    backward_low_de_rounds = 0
    for r in range(64, -1, -1):
        if de_backward[r] <= 8:
            backward_low_de_rounds = 64 - r
        else:
            break

    backward_results.append({
        'msg_idx': mi,
        'de_trace': de_backward,
        'da_trace': da_backward,
        'total_trace': hw_total_backward,
        'backward_low_de': backward_low_de_rounds,
        'DW_back': DW_back,
    })

    if mi < 5:
        print(f"  Msg {mi}: Backward low-De chain: {backward_low_de_rounds} rounds")
        print(f"    De (r60-64): {de_backward[60:65]}")
        print(f"    De (r56-60): {de_backward[56:61]}")
        print(f"    De (r48-56): {de_backward[48:57]}")
        print(f"    Da (r60-64): {da_backward[60:65]}")
        print(f"    Total HW (r56-64): {hw_total_backward[56:65]}")
        nz_dw = sum(1 for d in DW_back if d != 0)
        print(f"    Non-zero DW choices: {nz_dw}/64")

# ============================================================================
# PART 3: Enhanced Backward Chain — multi-bit DW search
# ============================================================================
print("\n" + "=" * 72)
print("PART 3: Enhanced Backward Chain — broader DW search")
print("=" * 72)
print("Try more DW candidates per round (random low-HW values).\n")

enhanced_results = []

for mi, W16 in enumerate(base_messages[:10]):
    if time.time() - start_time > 360:
        print(f"  [Time limit approaching, stopping at message {mi}]")
        break

    W = expand_schedule(W16)

    states = [list(IV)]
    a, b, c, d, e, f, g, h = list(IV)
    for t in range(64):
        T1 = add(h, Sig1(e), Ch(e, f, g), K[t], W[t])
        T2 = add(Sig0(a), Maj(a, b, c))
        h, g, f, e, d, c, b, a = g, f, e, add(d, T1), c, b, a, add(T1, T2)
        states.append([a, b, c, d, e, f, g, h])

    state_prime = list(states[64])
    de_backward = [0]
    da_backward = [0]
    hw_total_backward = [0]

    rng = random.Random(mi * 1000 + 0xBEEF)

    for t in range(63, -1, -1):
        best_total = 999
        best_de = 999
        best_prev_prime = None

        # Candidates: 0, single bits, random 2-bit, random low-HW
        candidates = [0]
        for bit in range(32):
            candidates.append(1 << bit)
        # 2-bit diffs
        for _ in range(64):
            b1 = rng.randint(0, 31)
            b2 = rng.randint(0, 31)
            candidates.append((1 << b1) | (1 << b2))
        # Random low HW (3-4 bits)
        for _ in range(64):
            v = 0
            for _ in range(rng.randint(1, 4)):
                v |= (1 << rng.randint(0, 31))
            candidates.append(v)
        # Also try subtraction (negative DW)
        for bit in range(32):
            candidates.append((-1 << bit) & MASK)

        for dw in candidates:
            Wp_t = (W[t] + dw) & MASK
            prev_prime = inverse_round(state_prime, Wp_t, K[t])
            de = hw(prev_prime[4] ^ states[t][4])
            da = hw(prev_prime[0] ^ states[t][0])
            total = sum(hw(prev_prime[i] ^ states[t][i]) for i in range(8))

            if total < best_total or (total == best_total and de < best_de):
                best_de = de
                best_da = da
                best_total = total
                best_prev_prime = prev_prime

        state_prime = best_prev_prime
        de_backward.append(best_de)
        da_backward.append(best_da)
        hw_total_backward.append(best_total)

    de_backward = de_backward[::-1]
    da_backward = da_backward[::-1]
    hw_total_backward = hw_total_backward[::-1]

    backward_low_de_rounds = 0
    for r in range(64, -1, -1):
        if de_backward[r] <= 8:
            backward_low_de_rounds = 64 - r
        else:
            break

    enhanced_results.append({
        'backward_low_de': backward_low_de_rounds,
        'de_trace': de_backward,
        'total_trace': hw_total_backward,
    })

    if mi < 5:
        print(f"  Msg {mi}: Enhanced backward low-De chain: {backward_low_de_rounds} rounds")
        print(f"    De (r56-64): {de_backward[56:65]}")
        print(f"    De (r48-56): {de_backward[48:57]}")
        print(f"    Total HW (r48-64): {hw_total_backward[48:65]}")

# ============================================================================
# PART 4: Statistical Reachability — random mid-state → collision at round 64
# ============================================================================
print("\n" + "=" * 72)
print("PART 4: Statistical Reachability")
print("=" * 72)
print("At round R, perturb one bit of state. Measure De at round 64.\n")

N_samples = 2000
test_start_rounds = [56, 48, 40, 32, 24, 16]

print(f"{'Start R':>8s}  {'Avg De[64]':>10s}  {'De[64]=0 frac':>14s}  {'Avg HW total':>13s}")
print("-" * 50)

for R in test_start_rounds:
    if time.time() - start_time > 400:
        print(f"  [Time limit]")
        break

    de64_list = []
    zero_de64 = 0
    total_hw_list = []

    for _ in range(N_samples):
        W16 = [random.randint(0, MASK) for _ in range(16)]
        W = expand_schedule(W16)

        state_R = sha256_rounds(W16, R)
        # Perturb one bit
        state_R_p = list(state_R)
        reg = random.randint(0, 7)
        bit = random.randint(0, 31)
        state_R_p[reg] ^= (1 << bit)

        # Run both to round 64
        final1 = sha256_rounds_from(state_R, W, R, 64)
        final2 = sha256_rounds_from(state_R_p, W, R, 64)

        de = hw(final1[4] ^ final2[4])
        total = sum(hw(final1[i] ^ final2[i]) for i in range(8))
        de64_list.append(de)
        total_hw_list.append(total)
        if de == 0:
            zero_de64 += 1

    avg_de = np.mean(de64_list)
    frac = zero_de64 / N_samples
    avg_total = np.mean(total_hw_list)
    print(f"{R:8d}  {avg_de:10.2f}  {frac:14.6f}  {avg_total:13.2f}")

# ============================================================================
# PART 5: Gap Analysis
# ============================================================================
print("\n" + "=" * 72)
print("PART 5: Gap Analysis & Complexity Estimate")
print("=" * 72)

# Forward reach
avg_forward = np.mean(forward_max_de0)

# Backward reach (basic)
basic_back = [r['backward_low_de'] for r in backward_results]
avg_basic_back = np.mean(basic_back) if basic_back else 0

# Backward reach (enhanced)
enhanced_back = [r['backward_low_de'] for r in enhanced_results]
avg_enhanced_back = np.mean(enhanced_back) if enhanced_back else 0

print(f"\n  Forward chain (simple 1-bit DW[0]):  De=0 through ~round {avg_forward:.1f}")
print(f"  Forward chain (Wang-optimized):      De=0 through ~round 16 (literature)")
print(f"  Backward chain (basic, 33 DW cands): {avg_basic_back:.1f} rounds from end")
print(f"  Backward chain (enhanced, ~225 DW):   {avg_enhanced_back:.1f} rounds from end")

gap_basic = 64 - 16 - avg_basic_back  # Using Wang's round 16 as forward
gap_enhanced = 64 - 16 - avg_enhanced_back

print(f"\n  Gap (Wang fwd + basic bkwd):    {gap_basic:.1f} rounds")
print(f"  Gap (Wang fwd + enhanced bkwd): {gap_enhanced:.1f} rounds")

# Each uncovered round has ~32-bit barrier (from CRAZY-2/previous experiments)
complexity_basic = gap_basic * 32 if gap_basic > 0 else 0
complexity_enhanced = gap_enhanced * 32 if gap_enhanced > 0 else 0

print(f"\n  Complexity estimate (basic):    ~2^{complexity_basic:.0f}")
print(f"  Complexity estimate (enhanced): ~2^{complexity_enhanced:.0f}")
print(f"  Birthday bound:                 2^128")
print(f"  Brute force:                    2^256")

# Detailed backward trace summary
print("\n  Backward De trace summary (enhanced, avg over messages):")
if enhanced_results:
    num_msgs = len(enhanced_results)
    for r in [64, 63, 62, 61, 60, 58, 56, 52, 48, 40, 32]:
        avg_de = np.mean([res['de_trace'][r] for res in enhanced_results])
        avg_tot = np.mean([res['total_trace'][r] for res in enhanced_results])
        print(f"    Round {r:2d}: avg De={avg_de:5.1f}, avg total HW={avg_tot:5.1f}")

# ============================================================================
# VERDICT
# ============================================================================
print("\n" + "=" * 72)
print("VERDICT")
print("=" * 72)

best_backward = max(avg_basic_back, avg_enhanced_back)
gap = 64 - 16 - best_backward

if best_backward >= 5:
    verdict = "ALIVE"
    print(f"  ALIVE: Backward chain extends {best_backward:.1f} rounds with low De")
    print(f"  Gap = {gap:.1f} rounds (forward 16 + backward {best_backward:.1f} = {16+best_backward:.1f} of 64)")
elif best_backward >= 2:
    verdict = "ANOMALY"
    print(f"  ANOMALY: Backward chain extends {best_backward:.1f} rounds")
    print(f"  Gap ~{gap:.1f} rounds")
else:
    verdict = "DEAD"
    print(f"  DEAD: Backward chain extends only {best_backward:.1f} rounds")
    print(f"  Gap >= {gap:.1f} rounds — no improvement over brute force approach")

print(f"\n  Interpretation:")
print(f"  - Forward (Wang): 16 rounds of De=0 (established in literature)")
print(f"  - Backward: {best_backward:.1f} rounds of low De from collision end")
print(f"  - Gap: {gap:.1f} rounds where BOTH forward and backward have high De")
print(f"  - Each gap round costs ~2^32 (one unconstrained addition)")
print(f"  - Estimated complexity: ~2^{gap*32:.0f}")
if gap * 32 >= 128:
    print(f"  - This EXCEEDS birthday bound (2^128) — no shortcut from backward analysis")
else:
    print(f"  - This is BELOW birthday bound — backward chain provides useful structure")

elapsed = time.time() - start_time
print(f"\n  Total runtime: {elapsed:.1f}s")
print("=" * 72)
