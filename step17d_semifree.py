#!/usr/bin/env python3
"""Step 17d: Semi-free-start collision for reduced SHA-256.

Semi-free-start: we choose the IV (chaining value).
This gives us 256 extra bits of freedom.

Strategy:
1. Choose IV so that e is a Sigma1 fixed point (ensures D[Sigma1]=De)
2. Choose IV so that a is a Sigma0 fixed point
3. With these, the differential propagates predictably
4. Use DW to cancel differences

Alternative: brute-force search for short reduced rounds.
For 9 rounds, birthday attack on internal state might work.
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

# ============================================================
# Part 1: Semi-free-start with fixed-point IV
# ============================================================
def semifree_fixedpoint():
    """Choose IV with e=0xffffffff (Sigma1 fixed point)
    and a=0xffffffff (Sigma0 fixed point).
    Then Sigma functions preserve +1 differences.
    """
    print("=" * 70)
    print("Part 1: Semi-free-start with Sigma fixed-point IV")
    print("=" * 70)

    # Set a=e=0xffffffff for fixed point behavior
    # Other registers: choose to satisfy Ch/Maj conditions
    # Ch condition: f = g (then DCh = 0 when De changes)
    # Maj condition: b = c (then DMaj = 0 when Da changes)

    # Ideal IV: a = 0xffffffff, b = c (equal), d = anything,
    #           e = 0xffffffff, f = g (equal), h = anything

    N = 2000000
    found = 0

    for trial in range(N):
        # Construct semi-free IV with conditions
        fg_val = random.randint(0, MASK32)
        bc_val = random.randint(0, MASK32)
        d_val = random.randint(0, MASK32)
        h_val = random.randint(0, MASK32)

        iv = [0xffffffff, bc_val, bc_val, d_val,
              0xffffffff, fg_val, fg_val, h_val]

        W = [random.randint(0, MASK32) for _ in range(16)]
        W2 = list(W)
        W2[0] = (W[0] + 1) & MASK32  # DW[0] = +1

        # Run 1 round
        s1 = list(sha256_step(iv, W[0], K[0]))
        s2 = list(sha256_step(iv, W2[0], K[0]))

        # Check Da and De
        da = (s1[0] - s2[0]) & MASK32
        de = (s1[4] - s2[4]) & MASK32

        if da == 1 and de == 1:
            # Perfect! The +1 difference preserved through Sigma and Ch/Maj
            # Now we need the next rounds to also work

            # Check: does f = g hold in new state? (for next Ch condition)
            # New state: f = e_old = 0xffffffff, g = f_old = fg_val
            # f != g unless fg_val = 0xffffffff
            # Similarly: b = a_old = 0xffffffff, c = b_old = bc_val
            # b != c unless bc_val = 0xffffffff

            # For next round conditions to hold, we need DIFFERENT conditions
            # This is the problem: conditions shift each round

            found += 1
            if found <= 5:
                print(f"  Trial {trial}: Da=+1, De=+1 after round 0 ✓")
                # Check round 1
                # At round 1: e = s1[4], a = s1[0]
                # For Sigma1 to preserve: need e to be fixed point
                e1_is_fp = (Sigma1(s1[4]) == s1[4])
                a1_is_fp = (Sigma0(s1[0]) == s1[0])
                # f = g? b = c?
                fg_eq = (s1[5] == s1[6])
                bc_eq = (s1[1] == s1[2])
                print(f"    Round 1 state: e_fp={e1_is_fp}, a_fp={a1_is_fp}, "
                      f"f=g:{fg_eq}, b=c:{bc_eq}")

    print(f"\n  Da=+1 AND De=+1 after round 0: {found}/{N} "
          f"(P={found/N:.6f}, expected ~1/4 from f=g, b=c)")

# ============================================================
# Part 2: Round-by-round condition satisfaction with fixed IV
# ============================================================
def round_by_round_conditions():
    """For the all-0xff / f=g, b=c IV, trace how conditions
    propagate and what probability each round has."""
    print("\n" + "=" * 70)
    print("Part 2: Condition satisfaction probability per round")
    print("=" * 70)

    N = 500000

    # Track how many rounds maintain Da=+1, De=+1
    max_rounds_list = []

    for trial in range(N):
        fg_val = random.randint(0, MASK32)
        bc_val = random.randint(0, MASK32)
        d_val = random.randint(0, MASK32)
        h_val = random.randint(0, MASK32)

        iv = [0xffffffff, bc_val, bc_val, d_val,
              0xffffffff, fg_val, fg_val, h_val]

        W = [random.randint(0, MASK32) for _ in range(16)]
        W2 = list(W)
        W2[0] = (W[0] + 1) & MASK32

        s1 = list(iv)
        s2 = list(iv)

        max_ok = 0
        for t in range(16):
            s1 = list(sha256_step(s1, W[t], K[t]))
            s2 = list(sha256_step(s2, W2[t] if t == 0 else W[t], K[t]))

            da = (s1[0] - s2[0]) & MASK32
            de = (s1[4] - s2[4]) & MASK32

            if da == 1 and de == 1:
                max_ok = t + 1
            else:
                break

        max_rounds_list.append(max_ok)

    print(f"\n  Rounds where Da=+1, De=+1 maintained (no DW correction):")
    for r in range(8):
        count = sum(1 for x in max_rounds_list if x >= r + 1)
        p = count / N
        log_p = f"2^{math.log2(p):.1f}" if p > 0 else "<2^-19"
        print(f"    >= {r+1} rounds: {count}/{N} ({log_p})")

# ============================================================
# Part 3: Direct reduced-round collision via birthday
# ============================================================
def birthday_reduced():
    """For very short rounds (2-5), try birthday attack."""
    print("\n" + "=" * 70)
    print("Part 3: Birthday attack on very short reduced SHA-256")
    print("=" * 70)

    for n_rounds in [2, 3, 4, 5]:
        print(f"\n  === {n_rounds}-round semi-free-start birthday ===")

        # Compute reduced-round hash for many random (IV, M) pairs
        # and look for collisions via birthday

        seen = {}
        N = min(2**20, 1048576)  # 1M trials
        collisions = 0

        for trial in range(N):
            # Random IV and message
            iv = [random.randint(0, MASK32) for _ in range(8)]
            W = [random.randint(0, MASK32) for _ in range(16)]

            s = list(iv)
            for t in range(n_rounds):
                s = list(sha256_step(s, W[t], K[t]))

            # Hash = IV + state
            h = tuple((iv[i] + s[i]) & MASK32 for i in range(8))

            if h in seen:
                collisions += 1
                if collisions <= 2:
                    iv2, W2 = seen[h]
                    print(f"    COLLISION at trial {trial}!")
                    print(f"    IV1 = {[f'0x{x:08x}' for x in iv[:4]]}")
                    print(f"    IV2 = {[f'0x{x:08x}' for x in iv2[:4]]}")
                    print(f"    H   = {[f'0x{x:08x}' for x in h[:4]]}")
            else:
                seen[h] = (iv, W)

        p_est = N * (N-1) / 2 / 2**256  # Expected for random
        print(f"    {n_rounds} rounds: {collisions} collisions in {N} trials "
              f"(random expected: {p_est:.2e})")

# ============================================================
# Part 4: Fixed-IV reduced-round collision with DW structure
# ============================================================
def fixed_iv_reduced_collision():
    """Search for actual reduced-round collision with standard IV.

    Use DW[0] = +1 and free message words for correction.
    For each trial, compute both hashes and check if they match.
    """
    print("\n" + "=" * 70)
    print("Part 4: Fixed-IV reduced-round collision (DW[0]=+1)")
    print("=" * 70)

    # For t rounds with standard IV:
    # We get 15 × 32 = 480 bits of freedom (W[1..15] free)
    # Output: 256 bits
    # Expected birthday complexity: 2^128

    # But with message modification:
    # For each round, we can FIX one 32-bit state register via DW
    # So effectively we reduce the number of independent constraints

    # Let's try: force BOTH Da=0 AND De=0 at rounds where possible
    # Da=0 AND De=0 requires T2_diff = d_diff (condition on state)
    # When this holds, we set DW to zero both.
    # When it doesn't, we force De=0 only.

    for n_rounds in [9, 13, 17]:
        print(f"\n  === {n_rounds}-round collision (standard IV, DW[0]=+1) ===")

        N = 2000000
        best_hw = 999
        collisions = 0

        for trial in range(N):
            W = [random.randint(0, MASK32) for _ in range(16)]
            W2 = list(W)
            W2[0] = (W[0] + 1) & MASK32

            s1 = list(IV)
            s2 = list(IV)

            # Round 0: fixed DW
            s1 = list(sha256_step(s1, W[0], K[0]))
            s2 = list(sha256_step(s2, W2[0], K[0]))

            # Rounds 1-15: try to force BOTH Da=0 and De=0
            for t in range(1, min(n_rounds, 16)):
                a1, b1, c1, d1, e1, f1, g1, h1 = s1
                a2, b2, c2, d2, e2, f2, g2, h2 = s2

                T1_p1 = (h1 + Sigma1(e1) + Ch(e1, f1, g1) + K[t]) & MASK32
                T1_p2 = (h2 + Sigma1(e2) + Ch(e2, f2, g2) + K[t]) & MASK32
                T2_1 = (Sigma0(a1) + Maj(a1, b1, c1)) & MASK32
                T2_2 = (Sigma0(a2) + Maj(a2, b2, c2)) & MASK32

                T2_diff = (T2_1 - T2_2) & MASK32
                d_diff = (d1 - d2) & MASK32

                if T2_diff == d_diff:
                    # Can zero both Da and De!
                    # De=0: DT1 = -d_diff → DW = DT1 - (h_diff + Sig1_diff + Ch_diff)
                    dw = ((d1 + T1_p1) - (d2 + T1_p2)) & MASK32
                    W2[t] = (W[t] + dw) & MASK32
                else:
                    # Can only zero De, Da will be nonzero
                    dw = ((d1 + T1_p1) - (d2 + T1_p2)) & MASK32
                    W2[t] = (W[t] + dw) & MASK32

                s1 = list(sha256_step(s1, W[t], K[t]))
                s2 = list(sha256_step(s2, W2[t], K[t]))

            # For rounds >= 16: use schedule
            if n_rounds > 16:
                Wexp = expand_message(W)
                W2exp = expand_message(W2)
                for t in range(16, n_rounds):
                    s1 = list(sha256_step(s1, Wexp[t], K[t]))
                    s2 = list(sha256_step(s2, W2exp[t], K[t]))

            # Check output
            h1 = tuple((IV[i] + s1[i]) & MASK32 for i in range(8))
            h2 = tuple((IV[i] + s2[i]) & MASK32 for i in range(8))
            thw = sum(hw(h1[i] ^ h2[i]) for i in range(8))

            if thw < best_hw:
                best_hw = thw
            if thw == 0:
                collisions += 1
                if collisions <= 2:
                    print(f"    COLLISION at trial {trial}!")
                    print(f"    W  = {[f'0x{w:08x}' for w in W[:4]]} ...")
                    print(f"    W' = {[f'0x{w:08x}' for w in W2[:4]]} ...")

        print(f"    Best HW={best_hw}, collisions={collisions}/{N}")


if __name__ == "__main__":
    random.seed(42)

    semifree_fixedpoint()
    round_by_round_conditions()
    birthday_reduced()
    fixed_iv_reduced_collision()

    print("\n" + "=" * 70)
    print("Step 17d Complete")
    print("=" * 70)
