#!/usr/bin/env python3
"""Step 17e: Corrected NB local collision with proper sign convention.

Bug fix: Da = a' - a (perturbed minus original) = +1, not -1.
s2 is the perturbed state (W2[0] = W[0]+1), so Da = s2[0] - s1[0].

Re-derive the local collision path with correct signs.
"""

import random
import math
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

def sha256_step(state, W_t, K_t):
    a, b, c, d, e, f, g, h = state
    T1 = (h + Sigma1(e) + Ch(e, f, g) + K_t + W_t) & MASK32
    T2 = (Sigma0(a) + Maj(a, b, c)) & MASK32
    return ((T1 + T2) & MASK32, a, b, c, (d + T1) & MASK32, e, f, g)

def expand_message(W16):
    W = list(W16)
    for t in range(16, 64):
        W.append((sigma1(W[t-2]) + W[t-7] + sigma0(W[t-15]) + W[t-16]) & MASK32)
    return W

def diff(s_orig, s_pert):
    """Compute additive differential: perturbed - original."""
    return tuple((s_pert[i] - s_orig[i]) & MASK32 for i in range(8))

def diff_hw(s1, s2):
    return sum(hw((s2[i] - s1[i]) & MASK32) for i in range(8))

# ============================================================
# Part 1: Verify round 0 differential
# ============================================================
def verify_round0():
    """Confirm Da=+1, De=+1 after round 0 with DW[0]=+1."""
    print("=" * 70)
    print("Part 1: Verify round 0 differential (DW[0]=+1)")
    print("=" * 70)

    random.seed(42)
    N = 10000
    all_ok = True

    for _ in range(N):
        W = [random.randint(0, MASK32) for _ in range(16)]
        W2 = list(W)
        W2[0] = (W[0] + 1) & MASK32

        s1 = list(sha256_step(IV, W[0], K[0]))
        s2 = list(sha256_step(IV, W2[0], K[0]))

        d = diff(s1, s2)  # perturbed - original
        if d[0] != 1 or d[4] != 1:
            all_ok = False
            break
        # All other diffs should be 0
        for i in [1, 2, 3, 5, 6, 7]:
            if d[i] != 0:
                all_ok = False
                break

    print(f"  Da=+1, De=+1 after round 0: {'ALWAYS ✓' if all_ok else 'FAILED'}")
    # Show one example
    W = [0x12345678] + [0]*15
    W2 = [0x12345679] + [0]*15
    s1 = list(sha256_step(IV, W[0], K[0]))
    s2 = list(sha256_step(IV, W2[0], K[0]))
    d = diff(s1, s2)
    reg_names = ['Da', 'Db', 'Dc', 'Dd', 'De', 'Df', 'Dg', 'Dh']
    print("  Example: " + ", ".join(f"{reg_names[i]}={d[i]}" for i in range(8)))

# ============================================================
# Part 2: Trace the ideal local collision (corrected signs)
# ============================================================
def trace_ideal_lc():
    """Trace the ideal local collision with corrected convention.

    After round 0: Da=+1, De=+1, others=0
    Round 1 (DW[1] = 0):
      State diffs: Db=Da_old=+1, Dc=0, Dd=0
                   Df=De_old=+1, Dg=0, Dh=0
      T1 = h + Sigma1(e) + Ch(e,f,g) + K + W
      DT1 = Dh + DSigma1(e, De=+1) + DCh(De=+1, Df=+1, Dg=0) + DW

      DSigma1: depends on specific e value
      DCh: Ch(e+1, f+1, g) - Ch(e, f, g) — not simple!

    We can't solve this analytically without knowing state values.
    Let's simulate and find what the actual (Da,De) are after each round.
    """
    print("\n" + "=" * 70)
    print("Part 2: Actual differential path per round (DW[0]=+1, DW[1..]=0)")
    print("=" * 70)

    random.seed(42)
    N = 100000

    # Track (Da, De) at each round
    da_counts = {}  # round -> {val: count}
    de_counts = {}

    for _ in range(N):
        W = [random.randint(0, MASK32) for _ in range(16)]
        W2 = list(W)
        W2[0] = (W[0] + 1) & MASK32

        s1 = list(IV)
        s2 = list(IV)

        for t in range(13):
            wt = W[t]
            wt2 = W2[t] if t == 0 else W[t]  # Only DW[0] differs

            s1 = list(sha256_step(s1, wt, K[t]))
            s2 = list(sha256_step(s2, wt2, K[t]))

            d = diff(s1, s2)
            da, de = d[0], d[4]

            if t not in da_counts:
                da_counts[t] = {}
                de_counts[t] = {}
            da_counts[t][da] = da_counts[t].get(da, 0) + 1
            de_counts[t][de] = de_counts[t].get(de, 0) + 1

    print(f"\n  Most common (Da, De) after each round:")
    print(f"  {'Round':>5} | {'Top Da':>12} {'P':>7} | {'Top De':>12} {'P':>7} | {'Da=0':>6} | {'De=0':>6}")
    print("  " + "-" * 65)
    for t in range(13):
        top_da = max(da_counts[t], key=da_counts[t].get)
        top_de = max(de_counts[t], key=de_counts[t].get)
        p_da = da_counts[t][top_da] / N
        p_de = de_counts[t][top_de] / N
        da_zero = da_counts[t].get(0, 0) / N
        de_zero = de_counts[t].get(0, 0) / N
        print(f"  {t:5d} | 0x{top_da:08x} {p_da:6.4f} | 0x{top_de:08x} {p_de:6.4f} | "
              f"{da_zero:6.4f} | {de_zero:6.4f}")

# ============================================================
# Part 3: NB-style — find DW[1..3,8] to close local collision
# ============================================================
def nb_local_collision_correct():
    """Correct NB local collision: DW at steps 0,1,2,3,8.

    After round 0: Da=+1, De=+1
    We want all state diffs = 0 after round 8.

    The shift register means:
    After round 0: Da=1, De=1
    After round 1: Db=1 (from Da), Df=1 (from De), Da=?, De=?
    After round 2: Bc=1, Dg=1, Db=Da1, Df=De1, Da=?, De=?
    After round 3: Dd=1, Dh=1, Dc=Da1, Dg=De1, Db=Da2, Df=De2, Da=?, De=?
    After round 4: Dd=Da1, Dh=De1, Dc=Da2, Dg=De2, Db=Da3, Df=De3, Da=?, De=?

    For state to be zero after round 8:
    Need Da_8 = 0, De_8 = 0 AND all shifted values = 0
    Shifted: Db_8=Da_7, Dc_8=Da_6, Dd_8=Da_5, Df_8=De_7, Dg_8=De_6, Dh_8=De_5
    So need Da_5=Da_6=Da_7=0 and De_5=De_6=De_7=0 AND Da_8=De_8=0

    Actually, need ALL state to be zero means:
    Da_t = 0 for t=5,6,7,8 and De_t = 0 for t=5,6,7,8.

    We control DW[0..3] and DW[8]. That's 5 control points.
    Constraints: Da=0 for t=5..8 and De=0 for t=5..8 = 8 constraints.
    But we only have 4 control points after DW[0] (DW[1..3,8]).

    Each DW controls De at the next step (via T1).
    Da is then determined. So we have 4 equations for De_t:
    - DW[1] → De_1
    - DW[2] → De_2
    - DW[3] → De_3
    - DW[8] → De_8

    For De=0 at t=5,6,7: no DW correction available — these must be zero naturally!
    This is the probabilistic part.
    """
    print("\n" + "=" * 70)
    print("Part 3: Corrected NB local collision search")
    print("=" * 70)

    random.seed(42)
    N = 2000000

    # Strategy: force De=0 at rounds 1,2,3 and 8.
    # Check if Da and De are also zero at rounds 4,5,6,7.

    best_diff_r8 = 999
    best_diff_r4 = 999
    collisions = 0
    de_zero_4567 = 0
    da_zero_4567 = 0

    for trial in range(N):
        W = [random.randint(0, MASK32) for _ in range(16)]
        W2 = list(W)
        W2[0] = (W[0] + 1) & MASK32

        s1 = list(IV)
        s2 = list(IV)

        # Round 0: fixed DW=+1
        s1 = list(sha256_step(s1, W[0], K[0]))
        s2 = list(sha256_step(s2, W2[0], K[0]))

        # Rounds 1-3: force De=0
        for t in range(1, 4):
            # Force De'=0: d' + T1' = d + T1
            # T1' = T1_partial' + W2[t], T1 = T1_partial + W[t]
            # d' + T1_partial' + W2[t] = d + T1_partial + W[t]
            # W2[t] = W[t] + (d + T1_partial) - (d' + T1_partial')
            T1_p1 = (s1[7] + Sigma1(s1[4]) + Ch(s1[4], s1[5], s1[6]) + K[t]) & MASK32
            T1_p2 = (s2[7] + Sigma1(s2[4]) + Ch(s2[4], s2[5], s2[6]) + K[t]) & MASK32
            dw = ((s1[3] + T1_p1) - (s2[3] + T1_p2)) & MASK32
            W2[t] = (W[t] + dw) & MASK32

            s1 = list(sha256_step(s1, W[t], K[t]))
            s2 = list(sha256_step(s2, W2[t], K[t]))

        # Rounds 4-7: no DW correction (natural propagation)
        all_zero_47 = True
        for t in range(4, 8):
            s1 = list(sha256_step(s1, W[t], K[t]))
            s2 = list(sha256_step(s2, W[t], K[t]))

            d = diff(s1, s2)
            if d[0] != 0 or d[4] != 0:
                all_zero_47 = False

        if all_zero_47:
            da_zero_4567 += 1

        # Check state diff after round 7
        d7 = diff(s1, s2)
        thw_7 = sum(hw(x) for x in d7)
        if thw_7 < best_diff_r4:
            best_diff_r4 = thw_7

        # Round 8: force De=0
        t = 8
        T1_p1 = (s1[7] + Sigma1(s1[4]) + Ch(s1[4], s1[5], s1[6]) + K[t]) & MASK32
        T1_p2 = (s2[7] + Sigma1(s2[4]) + Ch(s2[4], s2[5], s2[6]) + K[t]) & MASK32
        dw = ((s1[3] + T1_p1) - (s2[3] + T1_p2)) & MASK32
        W2[t] = (W[t] + dw) & MASK32

        s1 = list(sha256_step(s1, W[t], K[t]))
        s2 = list(sha256_step(s2, W2[t], K[t]))

        d8 = diff(s1, s2)
        thw_8 = sum(hw(x) for x in d8)

        if thw_8 < best_diff_r8:
            best_diff_r8 = thw_8
            if thw_8 <= 20:
                print(f"  [trial {trial}] After round 8: state diff HW = {thw_8}")
                print(f"    Diffs: {[f'0x{x:08x}' for x in d8]}")

        if thw_8 == 0:
            collisions += 1
            if collisions <= 3:
                print(f"  *** LOCAL COLLISION CLOSED at trial {trial}! ***")
                print(f"  W  = {[f'0x{w:08x}' for w in W[:9]]}")
                print(f"  W' = {[f'0x{w:08x}' for w in W2[:9]]}")

        if trial % 500000 == 0 and trial > 0:
            print(f"  ... {trial/1e6:.1f}M trials, best r7 diff={best_diff_r4}, "
                  f"best r8 diff={best_diff_r8}, closed={collisions}")

    p_zero = da_zero_4567 / N if N > 0 else 0
    print(f"\n  Results after {N} trials:")
    print(f"    P(Da,De=0 at rounds 4-7): {da_zero_4567}/{N} ({p_zero:.2e})")
    print(f"    Best state diff after round 7 (before correction): HW = {best_diff_r4}")
    print(f"    Best state diff after round 8 (with correction): HW = {best_diff_r8}")
    print(f"    Full closures (collision): {collisions}")

# ============================================================
# Part 4: What if we also correct Da? (force Da=0 + De=0)
# ============================================================
def force_da_de_zero():
    """At correction steps (1,2,3,8), force BOTH Da=0 and De=0.

    Only possible when T2_diff = d_diff (see step16c).
    Track how often this compatibility condition holds.
    """
    print("\n" + "=" * 70)
    print("Part 4: Da=De=0 compatibility at correction steps")
    print("=" * 70)

    random.seed(42)
    N = 1000000

    compat_all = 0  # All 4 correction steps compatible
    compat_count = [0, 0, 0, 0]  # Per correction step

    for trial in range(N):
        W = [random.randint(0, MASK32) for _ in range(16)]
        W2 = list(W)
        W2[0] = (W[0] + 1) & MASK32

        s1 = list(IV)
        s2 = list(IV)
        s1 = list(sha256_step(s1, W[0], K[0]))
        s2 = list(sha256_step(s2, W2[0], K[0]))

        all_compat = True
        for ci, t in enumerate([1, 2, 3, 8]):
            # First, advance to step t if needed (for t=8, do steps 4-7)
            if ci == 3:  # Step 8, need to do 4-7 first
                for tt in range(4, 8):
                    s1 = list(sha256_step(s1, W[tt], K[tt]))
                    s2 = list(sha256_step(s2, W[tt], K[tt]))  # W2[tt] = W[tt]

            T2_1 = (Sigma0(s1[0]) + Maj(s1[0], s1[1], s1[2])) & MASK32
            T2_2 = (Sigma0(s2[0]) + Maj(s2[0], s2[1], s2[2])) & MASK32
            T2_diff = (T2_1 - T2_2) & MASK32
            d_diff = (s1[3] - s2[3]) & MASK32

            compat = (T2_diff == d_diff)
            if compat:
                compat_count[ci] += 1
            else:
                all_compat = False

            # Actually do the step with De=0 forcing
            T1_p1 = (s1[7] + Sigma1(s1[4]) + Ch(s1[4], s1[5], s1[6]) + K[t]) & MASK32
            T1_p2 = (s2[7] + Sigma1(s2[4]) + Ch(s2[4], s2[5], s2[6]) + K[t]) & MASK32
            dw = ((s1[3] + T1_p1) - (s2[3] + T1_p2)) & MASK32
            W2[t] = (W[t] + dw) & MASK32
            s1 = list(sha256_step(s1, W[t], K[t]))
            s2 = list(sha256_step(s2, W2[t], K[t]))

        if all_compat:
            compat_all += 1

    print(f"\n  P(T2_diff = d_diff) at each correction step:")
    for ci, t in enumerate([1, 2, 3, 8]):
        p = compat_count[ci] / N
        log_p = f"2^{math.log2(p):.1f}" if p > 0 else "0"
        print(f"    Step {t}: P = {p:.6f} ({log_p})")

    p_all = compat_all / N
    log_p_all = f"2^{math.log2(p_all):.1f}" if p_all > 0 else "0"
    print(f"    All 4 steps: P = {p_all:.6f} ({log_p_all})")
    print(f"    => ~{compat_all} in {N} have dual Da=De=0 compatibility")

if __name__ == "__main__":
    verify_round0()
    trace_ideal_lc()
    nb_local_collision_correct()
    force_da_de_zero()

    print("\n" + "=" * 70)
    print("Step 17e Complete")
    print("=" * 70)
