#!/usr/bin/env python3
"""Step 17b (revised): Nikolic-Biryukov local collision for reduced SHA-256.

The actual NB approach:
1. Local collision at steps i to i+8
2. DW nonzero only at positions i, i+1, i+2, i+3, i+8
3. All state diffs return to zero after step i+8
4. Steps before i and after i+8 have zero diff → collision

Key insight: DW at each correction step is COMPUTED to cancel BOTH
Da and De simultaneously. Since Da-De = T2_diff-d_diff (independent
of DW), this is only possible when T2_diff-d_diff has the right value.
This is where sufficient conditions come in.

For a 21-step collision: local collision at steps 3-11.
Steps 0-2: no diff (3 free steps)
Steps 3-11: local collision
Steps 12-20: no diff
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

def state_equal(s1, s2):
    return all(s1[i] == s2[i] for i in range(8))

def state_diff_hw(s1, s2):
    return sum(hw((s1[i] - s2[i]) & MASK32) for i in range(8))

# ============================================================
# Core: Try local collision and check if it closes
# ============================================================
def try_local_collision(W, start, n_total_rounds):
    """Try a local collision starting at step `start`.

    DW[start] = +1. For steps start+1 to start+3, compute DW to
    force De=0. Then DW=0 for steps start+4 to start+7.
    At step start+8, compute DW to force De=0 (final correction).

    Return: (success, state_diff_hw, W2)
    where success means ALL state diffs are zero after step start+8.
    """
    W2 = list(W)
    W2[start] = (W[start] + 1) & MASK32  # DW[start] = +1

    s1 = list(IV)
    s2 = list(IV)

    # Steps before start: identical
    for t in range(start):
        s1 = list(sha256_step(s1, W[t], K[t]))
        s2 = list(sha256_step(s2, W[t], K[t]))

    # Step start: inject difference
    s1 = list(sha256_step(s1, W[start], K[start]))
    s2 = list(sha256_step(s2, W2[start], K[start]))

    # Steps start+1 to start+3: force De=0 via DW
    for t in range(start + 1, min(start + 4, 16)):
        a1, b1, c1, d1, e1, f1, g1, h1 = s1
        a2, b2, c2, d2, e2, f2, g2, h2 = s2
        T1_p1 = (h1 + Sigma1(e1) + Ch(e1, f1, g1) + K[t]) & MASK32
        T1_p2 = (h2 + Sigma1(e2) + Ch(e2, f2, g2) + K[t]) & MASK32
        # Force De_new = 0: need d1 + T1_p1 + W = d2 + T1_p2 + W'
        # DW = (d1 + T1_p1) - (d2 + T1_p2)
        dw = ((d1 + T1_p1) - (d2 + T1_p2)) & MASK32
        W2[t] = (W[t] + dw) & MASK32

        s1 = list(sha256_step(s1, W[t], K[t]))
        s2 = list(sha256_step(s2, W2[t], K[t]))

    # Steps start+4 to start+7: no correction (DW = 0)
    for t in range(start + 4, min(start + 8, 16)):
        s1 = list(sha256_step(s1, W[t], K[t]))
        s2 = list(sha256_step(s2, W2[t], K[t]))

    # Step start+8: final correction (force De=0)
    t = start + 8
    if t < 16:
        a1, b1, c1, d1, e1, f1, g1, h1 = s1
        a2, b2, c2, d2, e2, f2, g2, h2 = s2
        T1_p1 = (h1 + Sigma1(e1) + Ch(e1, f1, g1) + K[t]) & MASK32
        T1_p2 = (h2 + Sigma1(e2) + Ch(e2, f2, g2) + K[t]) & MASK32
        dw = ((d1 + T1_p1) - (d2 + T1_p2)) & MASK32
        W2[t] = (W[t] + dw) & MASK32
        s1 = list(sha256_step(s1, W[t], K[t]))
        s2 = list(sha256_step(s2, W2[t], K[t]))

    # Check: is the state equal now?
    success = state_equal(s1, s2)
    diff = state_diff_hw(s1, s2)

    if not success:
        return False, diff, W2, s1, s2

    # Continue remaining rounds (should stay equal)
    for t in range(start + 9, n_total_rounds):
        wt = W[t] if t < 16 else expand_message(W)[t]
        wt2 = W2[t] if t < 16 else expand_message(W2)[t]
        s1 = list(sha256_step(s1, wt, K[t]))
        s2 = list(sha256_step(s2, wt2, K[t]))

    # Check final hash
    h1 = tuple((IV[i] + s1[i]) & MASK32 for i in range(8))
    h2 = tuple((IV[i] + s2[i]) & MASK32 for i in range(8))
    hash_match = all(h1[i] == h2[i] for i in range(8))

    return hash_match, 0, W2, h1, h2

# ============================================================
# Part 1: Statistics for local collision closure probability
# ============================================================
def closure_probability():
    """Estimate P(local collision closes) for different start positions."""
    print("=" * 70)
    print("Part 1: Local collision closure probability")
    print("=" * 70)

    random.seed(42)
    N = 500000

    for start in range(8):
        success_count = 0
        best_diff = 999
        diff_sum = 0

        for trial in range(N):
            W = [random.randint(0, MASK32) for _ in range(16)]
            ok, diff, W2, _, _ = try_local_collision(W, start, start + 9)
            if ok:
                success_count += 1
            if diff < best_diff:
                best_diff = diff
            diff_sum += diff

        avg_diff = diff_sum / N
        p = success_count / N
        log_p = f"2^{__import__('math').log2(p):.1f}" if p > 0 else f"< 2^{-__import__('math').log2(N):.0f}"
        print(f"  Start={start}: P(close) = {p:.6f} ({log_p}), "
              f"avg diff HW = {avg_diff:.1f}, best = {best_diff}")

# ============================================================
# Part 2: Search for actual collision with message modification
# ============================================================
def search_collision():
    """Search for reduced-round collision using local collision + modification."""
    print("\n" + "=" * 70)
    print("Part 2: Collision search with local collision")
    print("=" * 70)

    random.seed(42)

    for start in [0, 2, 4, 6]:
        n_rounds = start + 12  # Local collision + 3 rounds after
        print(f"\n  === Start={start}, {n_rounds}-round collision ===")

        N = 2000000
        found = 0

        for trial in range(N):
            W = [random.randint(0, MASK32) for _ in range(16)]
            ok, diff, W2, h1, h2 = try_local_collision(W, start, n_rounds)

            if ok:
                found += 1
                if found <= 3:
                    print(f"    COLLISION at trial {trial}!")
                    print(f"    W[{start}]=0x{W[start]:08x}, W'[{start}]=0x{W2[start]:08x}")
                    # Show DW pattern
                    dws = [(W2[i] - W[i]) & MASK32 for i in range(16)]
                    nonzero = [(i, dws[i]) for i in range(16) if dws[i] != 0]
                    print(f"    DW nonzero at: {[(i, f'0x{d:08x}') for i, d in nonzero]}")
                    if isinstance(h1, tuple):
                        print(f"    H  = 0x{''.join(f'{x:08x}' for x in h1)}")
                        print(f"    H' = 0x{''.join(f'{x:08x}' for x in h2)}")

        p = found / N
        if p > 0:
            log_p = f"2^{__import__('math').log2(p):.1f}"
        else:
            log_p = f"< 2^{-__import__('math').log2(N):.0f}"
        print(f"    Found {found}/{N} collisions (P ≈ {log_p})")

# ============================================================
# Part 3: Force Da=0 too — what DW is needed?
# ============================================================
def analyze_da_residual():
    """After forcing De=0 at correction steps, analyze what Da looks like.
    If Da is also small, the local collision might close naturally."""
    print("\n" + "=" * 70)
    print("Part 3: Da residual analysis after De=0 forcing")
    print("=" * 70)

    random.seed(42)
    N = 100000
    start = 3

    # Track Da at each step of the local collision
    da_stats = {t: [] for t in range(start, start + 9)}

    for trial in range(N):
        W = [random.randint(0, MASK32) for _ in range(16)]
        W2 = list(W)
        W2[start] = (W[start] + 1) & MASK32

        s1 = list(IV)
        s2 = list(IV)

        for t in range(start):
            s1 = list(sha256_step(s1, W[t], K[t]))
            s2 = list(sha256_step(s2, W[t], K[t]))

        for t in range(start, start + 9):
            if t == start:
                s1 = list(sha256_step(s1, W[t], K[t]))
                s2 = list(sha256_step(s2, W2[t], K[t]))
            elif t <= start + 3 or t == start + 8:
                a1, b1, c1, d1, e1, f1, g1, h1 = s1
                a2, b2, c2, d2, e2, f2, g2, h2 = s2
                T1_p1 = (h1 + Sigma1(e1) + Ch(e1, f1, g1) + K[t]) & MASK32
                T1_p2 = (h2 + Sigma1(e2) + Ch(e2, f2, g2) + K[t]) & MASK32
                dw = ((d1 + T1_p1) - (d2 + T1_p2)) & MASK32
                W2[t] = (W[t] + dw) & MASK32
                s1 = list(sha256_step(s1, W[t], K[t]))
                s2 = list(sha256_step(s2, W2[t], K[t]))
            else:
                s1 = list(sha256_step(s1, W[t], K[t]))
                s2 = list(sha256_step(s2, W2[t], K[t]))

            da = (s1[0] - s2[0]) & MASK32
            de = (s1[4] - s2[4]) & MASK32
            da_stats[t].append((hw(da), hw(de)))

    print(f"\n  Local collision at steps {start}-{start+8} (De=0 forced at {start+1}-{start+3},{start+8}):")
    print(f"  {'Step':>5} | {'Avg Da HW':>10} | {'Avg De HW':>10} | {'Da=0':>5} | {'De=0':>5}")
    print("  " + "-" * 50)
    for t in range(start, start + 9):
        da_hws = [x[0] for x in da_stats[t]]
        de_hws = [x[1] for x in da_stats[t]]
        da_avg = sum(da_hws) / N
        de_avg = sum(de_hws) / N
        da_zero = da_hws.count(0)
        de_zero = de_hws.count(0)
        print(f"  {t:5d} | {da_avg:10.2f} | {de_avg:10.2f} | {da_zero:5d} | {de_zero:5d}")

    # After step start+8, check full state diff
    print(f"\n  After step {start+8}, full state equality check:")
    eq_count = 0
    for trial in range(N):
        W = [random.randint(0, MASK32) for _ in range(16)]
        ok, diff, _, _, _ = try_local_collision(W, start, start + 9)
        if ok:
            eq_count += 1

    p = eq_count / N
    log_p = f"2^{__import__('math').log2(p):.1f}" if p > 0 else f"< 2^{-__import__('math').log2(N):.0f}"
    print(f"    P(full state equal) = {eq_count}/{N} ({log_p})")

# ============================================================
# Part 4: Alternative — force BOTH Da=0 and De=0
# ============================================================
def force_both_da_de():
    """Can we find DW values that force BOTH Da=0 AND De=0?

    Da_new = T1 + T2, De_new = d + T1
    DW affects T1 only.
    Da_new - De_new = T2 - d (independent of DW)

    So Da_new = 0 AND De_new = 0 requires T2_diff = d_diff.
    This is a constraint on the state, not controllable via DW.

    Let's measure P(T2_diff = d_diff) at each step.
    """
    print("\n" + "=" * 70)
    print("Part 4: P(Da=De can both be zeroed) at each correction step")
    print("=" * 70)

    random.seed(42)
    N = 500000
    start = 3

    # For each correction step, check if T2_diff = d_diff
    compat_stats = {}

    for trial in range(N):
        W = [random.randint(0, MASK32) for _ in range(16)]
        W2 = list(W)
        W2[start] = (W[start] + 1) & MASK32

        s1 = list(IV)
        s2 = list(IV)

        for t in range(start):
            s1 = list(sha256_step(s1, W[t], K[t]))
            s2 = list(sha256_step(s2, W[t], K[t]))

        for t in range(start, start + 9):
            a1, b1, c1, d1, e1, f1, g1, h1 = s1
            a2, b2, c2, d2, e2, f2, g2, h2 = s2

            T2_1 = (Sigma0(a1) + Maj(a1, b1, c1)) & MASK32
            T2_2 = (Sigma0(a2) + Maj(a2, b2, c2)) & MASK32
            T2_diff = (T2_1 - T2_2) & MASK32
            d_diff = (d1 - d2) & MASK32

            # Da=De=0 possible iff T2_diff = d_diff
            compat = (T2_diff == d_diff)

            if t not in compat_stats:
                compat_stats[t] = 0
            if compat:
                compat_stats[t] += 1

            # Actually step forward
            if t == start:
                s1 = list(sha256_step(s1, W[t], K[t]))
                s2 = list(sha256_step(s2, W2[t], K[t]))
            elif t <= start + 3 or t == start + 8:
                T1_p1 = (h1 + Sigma1(e1) + Ch(e1, f1, g1) + K[t]) & MASK32
                T1_p2 = (h2 + Sigma1(e2) + Ch(e2, f2, g2) + K[t]) & MASK32
                dw = ((d1 + T1_p1) - (d2 + T1_p2)) & MASK32
                W2[t] = (W[t] + dw) & MASK32
                s1 = list(sha256_step(s1, W[t], K[t]))
                s2 = list(sha256_step(s2, W2[t], K[t]))
            else:
                s1 = list(sha256_step(s1, W[t], K[t]))
                s2 = list(sha256_step(s2, W2[t], K[t]))

    print(f"\n  P(T2_diff = d_diff) at each step of local collision (start={start}):")
    for t in sorted(compat_stats):
        p = compat_stats[t] / N
        log_p = f"2^{__import__('math').log2(p):.1f}" if p > 0 else "0"
        print(f"    Step {t}: P = {p:.6f} ({log_p})")


if __name__ == "__main__":
    closure_probability()
    search_collision()
    analyze_da_residual()
    force_both_da_de()

    print("\n" + "=" * 70)
    print("Step 17b (revised) Complete")
    print("=" * 70)
