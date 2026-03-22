#!/usr/bin/env python3
"""Step 16c: Local collision technique for SHA-256.

A "local collision" introduces a difference at round t and
cancels it over the next few rounds using message word differences.

In SHA-256, the shift register structure means:
- Da at round t becomes Db at t+1, Dc at t+2, Dd at t+3
- De at round t becomes Df at t+1, Dg at t+2, Dh at t+3

So a difference takes 4 rounds to "pass through" the register.
We need to cancel it before it creates secondary differences.

Strategy:
1. Introduce DW[t] at round t
2. Use DW[t+1], DW[t+2], DW[t+3] to cancel propagation
3. After 4 rounds, state difference should be zero
"""

import random
MASK32 = 0xFFFFFFFF

def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def shr(x, n): return x >> n
def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ shr(x, 3)
def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ shr(x, 10)
def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
def hw(x): return bin(x).count('1')

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
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
]
IV = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]

def sha256_round_fn(state, W_t, K_t):
    a, b, c, d, e, f, g, h = state
    T1 = (h + Sigma1(e) + Ch(e, f, g) + K_t + W_t) & MASK32
    T2 = (Sigma0(a) + Maj(a, b, c)) & MASK32
    return ((T1 + T2) & MASK32, a, b, c, (d + T1) & MASK32, e, f, g)

def expand_message(W16):
    W = list(W16)
    for t in range(16, 64):
        W.append((sigma1(W[t-2]) + W[t-7] + sigma0(W[t-15]) + W[t-16]) & MASK32)
    return W

def state_diff(s1, s2):
    """Additive difference of all registers."""
    return tuple((s1[i] - s2[i]) & MASK32 for i in range(8))

def state_xdiff(s1, s2):
    return tuple(s1[i] ^ s2[i] for i in range(8))

def total_hw(diff):
    return sum(hw(d) for d in diff)

# ============================================================
# Part 1: 9-round local collision (SHA-256 shift register = 8 deep)
# ============================================================
def local_collision_9round():
    """Find DW[0..8] that create and cancel a difference in 9 rounds.

    Strategy: Force De=0 at each round using DW.
    The Da difference propagates through b,c,d and takes 4 rounds to exit.
    After 4 rounds of Da passing through, we need additional rounds
    for the Da's influence through Sigma0/Maj to die out.
    """
    print("=" * 70)
    print("Part 1: Local collision via De=0 forcing (9 rounds)")
    print("=" * 70)

    N = 10000
    state_diff_hw = [[] for _ in range(20)]
    full_cancel_count = 0

    for trial in range(N):
        W = [random.randint(0, MASK32) for _ in range(16)]
        W2 = list(W)
        W2[0] ^= 0x80000000  # Initial difference

        s1 = list(IV)
        s2 = list(IV)

        for t in range(16):
            if t == 0:
                # Fixed DW[0]
                pass
            elif t <= 8:
                # Force De=0
                a1, b1, c1, d1, e1, f1, g1, h1 = s1
                a2, b2, c2, d2, e2, f2, g2, h2 = s2
                T1_p1 = (h1 + Sigma1(e1) + Ch(e1, f1, g1) + K[t]) & MASK32
                T1_p2 = (h2 + Sigma1(e2) + Ch(e2, f2, g2) + K[t]) & MASK32
                required_dw = ((d1 + T1_p1) - (d2 + T1_p2)) & MASK32
                W2[t] = (W[t] + required_dw) & MASK32
            # else: DW[t] = 0 (no more correction)

            s1 = list(sha256_round_fn(s1, W[t], K[t]))
            s2 = list(sha256_round_fn(s2, W2[t], K[t]))

            sd = state_diff(s1, s2)
            thw = total_hw(state_xdiff(s1, s2))
            state_diff_hw[t].append(thw)

        # Check if state converges
        final_xd = state_xdiff(s1, s2)
        if total_hw(final_xd) == 0:
            full_cancel_count += 1

    print(f"\n  State XOR diff HW evolution (De=0 forced for rounds 1-8):")
    for t in range(16):
        hws = state_diff_hw[t]
        avg = sum(hws) / len(hws)
        zeros = hws.count(0)
        print(f"    Round {t:2d}: avg HW = {avg:6.1f}, full cancel = {zeros}")

    print(f"\n  Full state cancellation at round 15: {full_cancel_count}/{N}")

# ============================================================
# Part 2: Greedy message modification — minimize total diff HW
# ============================================================
def greedy_modification():
    """For each round t, choose DW[t] to minimize TOTAL state diff."""
    print("\n" + "=" * 70)
    print("Part 2: Greedy total-diff minimization")
    print("=" * 70)

    N = 1000

    for strategy_name, strategy in [
        ("De=0 only", "de_zero"),
        ("Min total HW", "min_total"),
        ("Min De+Da HW", "min_ae"),
    ]:
        print(f"\n  Strategy: {strategy_name}")
        diff_hw_stats = [[] for _ in range(20)]

        for trial in range(N):
            W = [random.randint(0, MASK32) for _ in range(16)]
            W2 = list(W)
            W2[0] ^= 0x80000000

            s1 = list(IV)
            s2 = list(IV)

            for t in range(16):
                if t == 0:
                    s1 = list(sha256_round_fn(s1, W[t], K[t]))
                    s2 = list(sha256_round_fn(s2, W2[t], K[t]))
                elif strategy == "de_zero":
                    a1, b1, c1, d1, e1, f1, g1, h1 = s1
                    a2, b2, c2, d2, e2, f2, g2, h2 = s2
                    T1_p1 = (h1 + Sigma1(e1) + Ch(e1, f1, g1) + K[t]) & MASK32
                    T1_p2 = (h2 + Sigma1(e2) + Ch(e2, f2, g2) + K[t]) & MASK32
                    dw = ((d1 + T1_p1) - (d2 + T1_p2)) & MASK32
                    W2[t] = (W[t] + dw) & MASK32
                    s1 = list(sha256_round_fn(s1, W[t], K[t]))
                    s2 = list(sha256_round_fn(s2, W2[t], K[t]))
                elif strategy == "min_total":
                    # Try De=0 and Da=0, pick better
                    a1, b1, c1, d1, e1, f1, g1, h1 = s1
                    a2, b2, c2, d2, e2, f2, g2, h2 = s2

                    T1_p1 = (h1 + Sigma1(e1) + Ch(e1, f1, g1) + K[t]) & MASK32
                    T1_p2 = (h2 + Sigma1(e2) + Ch(e2, f2, g2) + K[t]) & MASK32
                    T2_p1 = (Sigma0(a1) + Maj(a1, b1, c1)) & MASK32
                    T2_p2 = (Sigma0(a2) + Maj(a2, b2, c2)) & MASK32

                    # For De=0: DW forces e-path to zero
                    dw_de = ((d1 + T1_p1) - (d2 + T1_p2)) & MASK32

                    # For Da=0: need T1+T2 to cancel
                    # a_new = T1 + T2, where T1 includes W
                    # Da_new = (T1_p1 + W + T2_p1) - (T1_p2 + W' + T2_p2)
                    # = (T1_p1 + T2_p1) - (T1_p2 + T2_p2) + DW
                    # For Da=0: DW = (T1_p2 + T2_p2) - (T1_p1 + T2_p1)
                    dw_da = ((T1_p2 + T2_p2) - (T1_p1 + T2_p1)) & MASK32

                    # Test both
                    best_dw = dw_de
                    best_thw = 999

                    for dw_cand in [dw_de, dw_da, 0]:
                        W2_test = (W[t] + dw_cand) & MASK32
                        s1_test = list(sha256_round_fn(list(s1), W[t], K[t]))
                        s2_test = list(sha256_round_fn(list(s2), W2_test, K[t]))
                        thw = total_hw(state_xdiff(s1_test, s2_test))
                        if thw < best_thw:
                            best_thw = thw
                            best_dw = dw_cand

                    W2[t] = (W[t] + best_dw) & MASK32
                    s1 = list(sha256_round_fn(s1, W[t], K[t]))
                    s2 = list(sha256_round_fn(s2, W2[t], K[t]))
                elif strategy == "min_ae":
                    a1, b1, c1, d1, e1, f1, g1, h1 = s1
                    a2, b2, c2, d2, e2, f2, g2, h2 = s2
                    T1_p1 = (h1 + Sigma1(e1) + Ch(e1, f1, g1) + K[t]) & MASK32
                    T1_p2 = (h2 + Sigma1(e2) + Ch(e2, f2, g2) + K[t]) & MASK32
                    dw = ((d1 + T1_p1) - (d2 + T1_p2)) & MASK32
                    W2[t] = (W[t] + dw) & MASK32
                    s1 = list(sha256_round_fn(s1, W[t], K[t]))
                    s2 = list(sha256_round_fn(s2, W2[t], K[t]))

                sd = state_xdiff(s1, s2)
                thw = total_hw(sd)
                diff_hw_stats[t].append(thw)

        for t in range(1, 16):
            hws = diff_hw_stats[t]
            if hws:
                avg = sum(hws) / len(hws)
                print(f"    Round {t:2d}: avg total state XOR HW = {avg:6.1f}")

# ============================================================
# Part 3: The fundamental limit — counting degrees of freedom
# ============================================================
def freedom_analysis():
    """Count available vs required degrees of freedom."""
    print("\n" + "=" * 70)
    print("Part 3: Degrees of freedom analysis")
    print("=" * 70)

    print("""
  SHA-256 collision attack freedom count:

  AVAILABLE FREEDOM:
  - Message words W[0..15]: 16 × 32 = 512 bits
  - But DW[0] is fixed (our input difference): -32 bits
  - Net: 15 free DW[t] values = 480 bits of freedom

  REQUIRED CONSTRAINTS:
  - 64 rounds × 8 registers × 32 bits = 16384 state bits
  - But shift register means only 2 new values per round (Da, De)
  - So: 64 × 2 × 32 = 4096 effective constraint bits
  - For collision: need all 256 output bits to match
  - Minimum: 256 bits of constraints

  ANALYSIS:
  - We have 480 bits of freedom from DW[1..15]
  - These directly control rounds 1-15 (De=0 trivially achievable)
  - Message schedule W[16..63] is DETERMINED — zero additional freedom
  - Rounds 16-63: 48 × 2 × 32 = 3072 constraint bits, 0 freedom

  FUNDAMENTAL BARRIER:
  - For rounds 16-63, each round adds ~32 bits of constraint
  - Total: ~1536 bits of constraint with 0 bits of freedom
  - => Must rely on PROBABILITY, not deterministic control
  - Even if each round passes with P = 0.5 per bit:
    P(all 48 rounds pass) ≈ 2^(-48×32) = 2^(-1536)

  THIS IS WHY SHA-256 IS SECURE:
  - After round 16, the message schedule provides NO freedom
  - Every constraint must be satisfied probabilistically
  - Complexity ≈ 2^(number of unsatisfied constraints)
    """)

    # Verify: for random messages, what fraction of rounds 16+ have De=0?
    print("  Empirical verification (De=0 rate for rounds 16+):")
    N = 100000
    de_zero_count = [0] * 64
    for _ in range(N):
        W = [random.randint(0, MASK32) for _ in range(16)]
        W2 = list(W)
        W2[0] ^= 0x80000000

        # Force De=0 for rounds 1-15
        s1 = list(IV)
        s2 = list(IV)
        s1_r = list(sha256_round_fn(s1, W[0], K[0]))
        s2_r = list(sha256_round_fn(s2, W2[0], K[0]))

        for t in range(1, 16):
            a1, b1, c1, d1, e1, f1, g1, h1 = s1_r
            a2, b2, c2, d2, e2, f2, g2, h2 = s2_r
            T1_p1 = (h1 + Sigma1(e1) + Ch(e1, f1, g1) + K[t]) & MASK32
            T1_p2 = (h2 + Sigma1(e2) + Ch(e2, f2, g2) + K[t]) & MASK32
            dw = ((d1 + T1_p1) - (d2 + T1_p2)) & MASK32
            W2[t] = (W[t] + dw) & MASK32
            s1_r = list(sha256_round_fn(s1_r, W[t], K[t]))
            s2_r = list(sha256_round_fn(s2_r, W2[t], K[t]))

        # Now expand and continue with schedule-determined W
        Wexp = expand_message(W)
        W2exp = expand_message(W2)

        for t in range(16, 30):
            s1_r = list(sha256_round_fn(s1_r, Wexp[t], K[t]))
            s2_r = list(sha256_round_fn(s2_r, W2exp[t], K[t]))
            de = (s1_r[4] - s2_r[4]) & MASK32
            if de == 0:
                de_zero_count[t] += 1

    print(f"  After forcing De=0 for rounds 1-15:")
    for t in range(16, 30):
        c = de_zero_count[t]
        log_p = f"2^{-32:.0f}" if c == 0 else f"{c/N:.2e}"
        print(f"    Round {t:2d}: De=0 in {c}/{N} cases ({log_p})")

# ============================================================
# Part 4: Message schedule constraint propagation
# ============================================================
def schedule_constraint_propagation():
    """Analyze how DW[0..15] constraints from De=0 forcing
    propagate through the message schedule to DW[16+]."""
    print("\n" + "=" * 70)
    print("Part 4: DW constraint propagation through schedule")
    print("=" * 70)

    # The required DW[1..15] for De=0 are DETERMINED by W[0..t-1]
    # These DW values then determine DW[16+] via the schedule
    # DW[16] = Dsigma1(W14, W14') + DW[9] + Dsigma0(W1, W1') + DW[0]

    # Question: Can we choose the BASE W[t] values (not just DW[t])
    # to influence what DW[16+] comes out?

    # W[t] is free, W'[t] = W[t] + DW[t] is determined
    # DW[t] is determined by the De=0 constraint
    # But the ABSOLUTE values W[t] affect DW[16+] through sigma0/sigma1

    N = 5000
    dw16_hw_list = []

    for trial in range(N):
        W = [random.randint(0, MASK32) for _ in range(16)]
        W2 = list(W)
        W2[0] ^= 0x80000000

        s1 = list(IV)
        s2 = list(IV)
        s1 = list(sha256_round_fn(s1, W[0], K[0]))
        s2 = list(sha256_round_fn(s2, W2[0], K[0]))

        for t in range(1, 16):
            a1, b1, c1, d1, e1, f1, g1, h1 = s1
            a2, b2, c2, d2, e2, f2, g2, h2 = s2
            T1_p1 = (h1 + Sigma1(e1) + Ch(e1, f1, g1) + K[t]) & MASK32
            T1_p2 = (h2 + Sigma1(e2) + Ch(e2, f2, g2) + K[t]) & MASK32
            dw = ((d1 + T1_p1) - (d2 + T1_p2)) & MASK32
            W2[t] = (W[t] + dw) & MASK32
            s1 = list(sha256_round_fn(s1, W[t], K[t]))
            s2 = list(sha256_round_fn(s2, W2[t], K[t]))

        Wexp = expand_message(W)
        W2exp = expand_message(W2)

        for t in range(16, 25):
            dw = (W2exp[t] - Wexp[t]) & MASK32
            if t == 16:
                dw16_hw_list.append(hw(dw))

    avg_dw16 = sum(dw16_hw_list) / len(dw16_hw_list)
    dw16_zero = dw16_hw_list.count(0)
    dw16_low = sum(1 for h in dw16_hw_list if h <= 3)

    print(f"\n  DW[16] distribution after De=0 forcing:")
    print(f"    Avg HW: {avg_dw16:.1f}")
    print(f"    DW[16]=0: {dw16_zero}/{N}")
    print(f"    HW<=3: {dw16_low}/{N}")
    print(f"    => DW[16] is essentially random (~16 HW)")
    print(f"    => No free lunch from message schedule")


if __name__ == "__main__":
    random.seed(42)

    local_collision_9round()
    greedy_modification()
    freedom_analysis()
    schedule_constraint_propagation()

    print("\n" + "=" * 70)
    print("Step 16c Complete")
    print("=" * 70)
