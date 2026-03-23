#!/usr/bin/env python3
"""Step 17b: 9-step local collision for SHA-256.

Based on Nikolic-Biryukov (FSE 2008):
- Introduce modular difference +1 at step i via DW[i]
- Cancel it over steps i+1 through i+8
- Total: differences in DW at positions i, i+1, i+2, i+3, i+8

The state differential path for a local collision starting at step i:
  Step i:   Da=+1, De=+1 (from DW[i]=+1)
  Step i+1: Da=?, De=?   (Db=+1 from previous Da)
  Step i+2: Da=?, De=?   (Dc=+1, Df=+1)
  Step i+3: Da=?, De=?   (Dd=+1, Dg=+1)
  Step i+4: Da=0, De=0   (Dh=+1 from original De, but DW cancels)

Actually, let me derive the exact path by simulation.
"""

import random
from collections import defaultdict
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

def compress_full(W64, iv=None):
    if iv is None: iv = IV
    s = list(iv)
    for t in range(64):
        s = list(sha256_round_fn(s, W64[t], K[t]))
    return tuple((iv[i] + s[i]) & MASK32 for i in range(8))

# ============================================================
# Part 1: Derive the ideal differential path for local collision
# ============================================================
def derive_ideal_path():
    """Derive the differential path assuming Sigma preserves +1."""
    print("=" * 70)
    print("Part 1: Ideal 9-step local collision differential path")
    print("=" * 70)

    # SHA-256 round function:
    # T1 = h + Sigma1(e) + Ch(e,f,g) + K + W
    # T2 = Sigma0(a) + Maj(a,b,c)
    # a_new = T1 + T2
    # e_new = d + T1
    # b_new = a, c_new = b, d_new = c
    # f_new = e, g_new = f, h_new = g

    # State = (a, b, c, d, e, f, g, h)

    # Ideal assumptions (Nikolic-Biryukov):
    # 1. Sigma0 and Sigma1 preserve modular differences
    #    i.e., D[Sigma(x)] = Sigma(x + Dx) - Sigma(x) = Dx
    # 2. Ch and Maj have zero differential when Df=Dg=0 or Db=Dc=0
    #    DCh = 0 when only De changes AND f=g at affected bits
    #    DMaj = 0 when only Da changes AND b=c at affected bits

    # Under these ideal assumptions, let's trace a DW[0] = +1:

    print("\n  Tracing ideal differential path (DW[0] = +1):")
    print("  Assumptions: D[Sigma] = Dx, DCh=0 when f=g, DMaj=0 when b=c\n")

    # Initial state: Da=Db=Dc=Dd=De=Df=Dg=Dh = 0
    # Labels: registers at step t BEFORE the round function
    # After step 0:
    #   DT1 = Dh + D[Sigma1(e)] + DCh + DK + DW = 0 + 0 + 0 + 0 + 1 = 1
    #   DT2 = D[Sigma0(a)] + DMaj = 0 + 0 = 0
    #   Da_new = DT1 + DT2 = 1
    #   De_new = Dd + DT1 = 0 + 1 = 1
    # State diff after step 0: (1, 0, 0, 0, 1, 0, 0, 0)

    # After step 1 (with DW[1] to be determined):
    #   Dh_1 = Dg_0 = 0
    #   De_1 = 1 → D[Sigma1] = 1 (ideal)
    #   Df_1 = De_0 = 1 → DCh depends on De_1 and Df_1
    #   Wait — De_1 is the old De which is now at register e
    #   Let me be more careful with register naming.

    # Let's use arrays. State after round t: S[t] = (a_t, b_t, c_t, d_t, e_t, f_t, g_t, h_t)
    # S[t+1][1] = S[t][0]  (b = prev_a)
    # S[t+1][2] = S[t][1]  (c = prev_b)
    # S[t+1][3] = S[t][2]  (d = prev_c)
    # S[t+1][5] = S[t][4]  (f = prev_e)
    # S[t+1][6] = S[t][5]  (g = prev_f)
    # S[t+1][7] = S[t][6]  (h = prev_g)

    reg_names = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']

    # Track additive diffs
    D = [0] * 8  # Current differential
    DW_needed = {}

    for step in range(10):
        # Shift register part
        Dh = D[6]  # h_new = g_old
        Dg = D[5]  # g_new = f_old
        Df = D[4]  # f_new = e_old
        Dd = D[2]  # d_new = c_old
        Dc = D[1]  # c_new = b_old
        Db = D[0]  # b_new = a_old

        # For T1: uses h, e, f, g from CURRENT state (before shift)
        # DT1 = D[h] + D[Sigma1(e)] + DCh(e,f,g) + DW
        # Ideal: D[Sigma1(e)] = D[e] when De small and Sigma preserves
        # DCh: need to determine case by case

        # Ch(e,f,g) = ef XOR (NOT e)g
        # When De != 0 and Df = Dg:
        #   If Df = Dg = 0: DCh = De * (f XOR g) at each bit → needs f=g for DCh=0
        # When only Df != 0: DCh = e * Df → needs e=0 at affected bits
        # When only Dg != 0: DCh = (NOT e) * Dg → needs e=1 at affected bits

        # For simplicity, assume ideal: DCh = 0 when only De changes (f=g condition)
        # and DCh = 0 in other cases too (sufficient conditions)
        # DMaj = 0 similarly (b=c condition when Da != 0)

        # Under ideal assumptions:
        D_Sigma1_e = D[4]  # D[e]
        D_Sigma0_a = D[0]  # D[a]
        DCh = 0  # ideal
        DMaj = 0  # ideal

        # What DW do we need?
        # We want Da_new = 0 and De_new = 0 after the local collision closes
        # But during the collision, we CHOOSE DW to control the path

        # The round function computes:
        # DT1 = D[h] + D_Sigma1_e + DCh + DW
        # DT2 = D_Sigma0_a + DMaj
        # Da_new = DT1 + DT2
        # De_new = D[d] + DT1

        # Given: we want a SPECIFIC (Da_new, De_new) target
        # From De_new: DT1 = De_new - D[d]
        # From Da_new: DT2 = Da_new - DT1
        # But DT2 = D_Sigma0_a + DMaj is DETERMINED by state
        # So Da_new = DT1 + DT2 is determined once DT1 is chosen
        # We can only control DT1 via DW:
        # DW = DT1 - D[h] - D_Sigma1_e - DCh

        # Strategy: force De_new to specific target
        if step <= 3:
            # During perturbation, accept the natural differential
            DW_step = 0
            if step == 0:
                DW_step = 1  # Inject the difference
            elif step == 1:
                # Cancel De: want De_new = 0
                # DT1 = De_new - D[d] = 0 - Dd
                # DW = DT1 - D[h] - D_Sigma1_e - DCh
                target_De = 0
                DT1 = (target_De - D[2]) & MASK32  # Dd from current state
                DW_step = (DT1 - D[7] - D_Sigma1_e - DCh) & MASK32
            elif step == 2:
                target_De = 0
                DT1 = (target_De - D[2]) & MASK32
                DW_step = (DT1 - D[7] - D_Sigma1_e - DCh) & MASK32
            elif step == 3:
                target_De = 0
                DT1 = (target_De - D[2]) & MASK32
                DW_step = (DT1 - D[7] - D_Sigma1_e - DCh) & MASK32
            else:
                DW_step = 0

            DW_needed[step] = DW_step

            DT1 = (D[7] + D_Sigma1_e + DCh + DW_step) & MASK32
            DT2 = (D_Sigma0_a + DMaj) & MASK32
            Da_new = (DT1 + DT2) & MASK32
            De_new = (D[2] + DT1) & MASK32  # d = c from 2 steps ago? No, D[2] is current c
            # Wait — d_new = c_old, but De_new = d_old + DT1
            De_new = (D[3] + DT1) & MASK32

            # Correct: De_new = d_current + T1, where d_current = D[3]
            Da_new = (DT1 + DT2) & MASK32
        else:
            # After perturbation: DW = 0 (or specific correction)
            DW_step = 0

            # Check if we need a correction at step 8
            DT1 = (D[7] + D_Sigma1_e + DCh + DW_step) & MASK32
            DT2 = (D_Sigma0_a + DMaj) & MASK32
            Da_new = (DT1 + DT2) & MASK32
            De_new = (D[3] + DT1) & MASK32

            DW_needed[step] = DW_step

        # Update state diffs
        D = [Da_new, Db, Dc, Dd, De_new, Df, Dg, Dh]

        diff_str = " ".join(f"D{reg_names[i]}={D[i]:+d}" if D[i] != 0 and D[i] < 0x80000000
                           else f"D{reg_names[i]}=0" if D[i] == 0
                           else f"D{reg_names[i]}=0x{D[i]:08x}"
                           for i in range(8))
        dw_str = f"DW={DW_step:+d}" if DW_step < 0x80000000 else f"DW=0x{DW_step:08x}"
        print(f"  Step {step}: {dw_str:16s} | {diff_str}")

    return DW_needed

# ============================================================
# Part 2: Empirical search — 21-step collision attempt
# ============================================================
def search_21step_collision():
    """Search for 21-step reduced SHA-256 collision.

    Nikolic-Biryukov approach:
    - Local collision starts at step 3 (first 3 steps = free)
    - DW[3] = +1, corrections at DW[4], DW[5], DW[6], DW[11]
    - Total: 21 steps with collision

    We use a simpler approach: try random messages and check if
    the reduced-round compression produces a collision.
    """
    print("\n" + "=" * 70)
    print("Part 2: Searching for reduced-round collisions")
    print("=" * 70)

    # For a t-round reduced collision:
    # H = compress_t_rounds(IV, M) = compress_t_rounds(IV, M')
    # where DM[i] = specific pattern

    # Start simple: DW[0] = 1 (additive), check reduced rounds
    for n_rounds in [9, 13, 17, 21]:
        print(f"\n  === {n_rounds}-round collision search (DW[0]=+1) ===")
        N = 500000

        best_hw = 999
        best_msg = None

        for trial in range(N):
            W = [random.randint(0, MASK32) for _ in range(16)]
            W2 = list(W)
            W2[0] = (W[0] + 1) & MASK32  # Additive +1

            Wexp = expand_message(W)
            W2exp = expand_message(W2)

            s1 = list(IV)
            s2 = list(IV)
            for t in range(n_rounds):
                s1 = list(sha256_round_fn(s1, Wexp[t], K[t]))
                s2 = list(sha256_round_fn(s2, W2exp[t], K[t]))

            # Check hash difference (with feedforward)
            h1 = tuple((IV[i] + s1[i]) & MASK32 for i in range(8))
            h2 = tuple((IV[i] + s2[i]) & MASK32 for i in range(8))

            total_diff_hw = sum(hw(h1[i] ^ h2[i]) for i in range(8))

            if total_diff_hw < best_hw:
                best_hw = total_diff_hw
                best_msg = W
                if total_diff_hw == 0:
                    print(f"    COLLISION FOUND at trial {trial}!")
                    print(f"    W = {[f'0x{w:08x}' for w in W]}")
                    break

        print(f"    Best: total hash diff HW = {best_hw} (from {N} trials)")

    # Now try with DW at later positions (Nikolic-Biryukov style)
    # Local collision starting at step 3: DW[3] = +1
    print(f"\n  === DW[3]=+1 (local collision at step 3) ===")
    for n_rounds in [13, 17, 21]:
        N = 500000
        best_hw = 999

        for trial in range(N):
            W = [random.randint(0, MASK32) for _ in range(16)]
            W2 = list(W)
            W2[3] = (W[3] + 1) & MASK32

            Wexp = expand_message(W)
            W2exp = expand_message(W2)

            s1 = list(IV)
            s2 = list(IV)
            for t in range(n_rounds):
                s1 = list(sha256_round_fn(s1, Wexp[t], K[t]))
                s2 = list(sha256_round_fn(s2, W2exp[t], K[t]))

            h1 = tuple((IV[i] + s1[i]) & MASK32 for i in range(8))
            h2 = tuple((IV[i] + s2[i]) & MASK32 for i in range(8))
            total_diff_hw = sum(hw(h1[i] ^ h2[i]) for i in range(8))

            if total_diff_hw < best_hw:
                best_hw = total_diff_hw
                if total_diff_hw == 0:
                    print(f"    {n_rounds}-round COLLISION at trial {trial}!")
                    break

        print(f"    {n_rounds} rounds, DW[3]=+1: best hash diff HW = {best_hw}")

# ============================================================
# Part 3: Multi-word DW pattern from Nikolic-Biryukov
# ============================================================
def nikolic_biryukov_pattern():
    """Try the actual NB pattern: DW at positions i, i+1, i+2, i+3, i+8."""
    print("\n" + "=" * 70)
    print("Part 3: Nikolic-Biryukov multi-word DW pattern")
    print("=" * 70)

    # The NB local collision at step i uses DW at i, i+1, i+2, i+3, i+8
    # The DW values are determined by the state (sufficient conditions)
    # For now, let's try: DW[i]=+1, and determine DW[i+1..i+3,i+8]
    # by forcing De=0 at each intermediate step

    for start in [0, 2, 4, 6, 8]:
        print(f"\n  Local collision starting at step {start}:")

        N = 200000
        best_hw = {9: 999, 13: 999, 17: 999, 21: 999}

        for trial in range(N):
            W = [random.randint(0, MASK32) for _ in range(16)]
            W2 = list(W)

            # Run steps before the local collision (identical)
            s1 = list(IV)
            s2 = list(IV)
            Wexp = expand_message(W)

            for t in range(start):
                s1 = list(sha256_round_fn(s1, Wexp[t], K[t]))
                s2 = list(sha256_round_fn(s2, Wexp[t], K[t]))

            # Step `start`: inject DW = +1
            W2[start] = (W[start] + 1) & MASK32

            # Steps start+1 to start+3: force De=0
            for t in range(start, min(start + 9, 16)):
                if t == start:
                    s1 = list(sha256_round_fn(s1, W[t], K[t]))
                    s2 = list(sha256_round_fn(s2, W2[t], K[t]))
                elif t <= start + 3 and t < 16:
                    # Force De=0 by adjusting DW
                    a1, b1, c1, d1, e1, f1, g1, h1 = s1
                    a2, b2, c2, d2, e2, f2, g2, h2 = s2
                    T1_p1 = (h1 + Sigma1(e1) + Ch(e1, f1, g1) + K[t]) & MASK32
                    T1_p2 = (h2 + Sigma1(e2) + Ch(e2, f2, g2) + K[t]) & MASK32
                    dw = ((d1 + T1_p1) - (d2 + T1_p2)) & MASK32
                    W2[t] = (W[t] + dw) & MASK32
                    s1 = list(sha256_round_fn(s1, W[t], K[t]))
                    s2 = list(sha256_round_fn(s2, W2[t], K[t]))
                elif t == start + 8 and t < 16:
                    # Correction at step i+8
                    a1, b1, c1, d1, e1, f1, g1, h1 = s1
                    a2, b2, c2, d2, e2, f2, g2, h2 = s2
                    T1_p1 = (h1 + Sigma1(e1) + Ch(e1, f1, g1) + K[t]) & MASK32
                    T1_p2 = (h2 + Sigma1(e2) + Ch(e2, f2, g2) + K[t]) & MASK32
                    dw = ((d1 + T1_p1) - (d2 + T1_p2)) & MASK32
                    W2[t] = (W[t] + dw) & MASK32
                    s1 = list(sha256_round_fn(s1, W[t], K[t]))
                    s2 = list(sha256_round_fn(s2, W2[t], K[t]))
                else:
                    s1 = list(sha256_round_fn(s1, W[t], K[t]))
                    s2 = list(sha256_round_fn(s2, W2[t], K[t]))

            # Continue with expanded schedule
            W2exp = expand_message(W2)

            for t in range(min(start + 9, 16), 21):
                wt1 = Wexp[t] if t < len(Wexp) else 0
                wt2 = W2exp[t] if t < len(W2exp) else 0
                s1 = list(sha256_round_fn(s1, wt1, K[t]))
                s2 = list(sha256_round_fn(s2, wt2, K[t]))

                round_num = t + 1
                if round_num in best_hw:
                    h1 = tuple((IV[j] + s1[j]) & MASK32 for j in range(8))
                    h2 = tuple((IV[j] + s2[j]) & MASK32 for j in range(8))
                    thw = sum(hw(h1[j] ^ h2[j]) for j in range(8))
                    if thw < best_hw[round_num]:
                        best_hw[round_num] = thw
                        if thw == 0:
                            print(f"    {round_num}-round COLLISION at trial {trial}!")

        for nr in sorted(best_hw):
            print(f"    {nr} rounds: best hash diff HW = {best_hw[nr]}")


# ============================================================
# Part 4: De=0 forcing for all modified rounds
# ============================================================
def full_de_forcing():
    """Force De=0 at ALL modifiable rounds (0-15), then measure
    how far the collision extends."""
    print("\n" + "=" * 70)
    print("Part 4: Full De=0 forcing + reduced-round collision check")
    print("=" * 70)

    N = 1000000

    # Track minimum hash diff for each number of rounds
    best_hw_by_rounds = {}
    collision_found = {}

    for trial in range(N):
        W = [random.randint(0, MASK32) for _ in range(16)]
        W2 = list(W)
        W2[0] = (W[0] + 1) & MASK32  # DW[0] = +1

        s1 = list(IV)
        s2 = list(IV)

        # Round 0: fixed DW
        s1 = list(sha256_round_fn(s1, W[0], K[0]))
        s2 = list(sha256_round_fn(s2, W2[0], K[0]))

        # Rounds 1-15: force De=0
        for t in range(1, 16):
            a1, b1, c1, d1, e1, f1, g1, h1 = s1
            a2, b2, c2, d2, e2, f2, g2, h2 = s2
            T1_p1 = (h1 + Sigma1(e1) + Ch(e1, f1, g1) + K[t]) & MASK32
            T1_p2 = (h2 + Sigma1(e2) + Ch(e2, f2, g2) + K[t]) & MASK32
            dw = ((d1 + T1_p1) - (d2 + T1_p2)) & MASK32
            W2[t] = (W[t] + dw) & MASK32
            s1 = list(sha256_round_fn(s1, W[t], K[t]))
            s2 = list(sha256_round_fn(s2, W2[t], K[t]))

        # Expand messages
        Wexp = expand_message(W)
        W2exp = expand_message(W2)

        # Continue rounds 16+
        for t in range(16, 21):
            s1 = list(sha256_round_fn(s1, Wexp[t], K[t]))
            s2 = list(sha256_round_fn(s2, W2exp[t], K[t]))

            nr = t + 1
            h1 = tuple((IV[j] + s1[j]) & MASK32 for j in range(8))
            h2 = tuple((IV[j] + s2[j]) & MASK32 for j in range(8))
            thw = sum(hw(h1[j] ^ h2[j]) for j in range(8))

            if nr not in best_hw_by_rounds or thw < best_hw_by_rounds[nr]:
                best_hw_by_rounds[nr] = thw
                if thw == 0 and nr not in collision_found:
                    collision_found[nr] = trial
                    print(f"    {nr}-round COLLISION at trial {trial}!")
                    print(f"    W = {[f'0x{w:08x}' for w in W[:8]]}")
                    print(f"    W2= {[f'0x{w:08x}' for w in W2[:8]]}")

        if trial % 200000 == 0 and trial > 0:
            print(f"    ... {trial} trials done, best: " +
                  ", ".join(f"{nr}r={best_hw_by_rounds.get(nr, '?')}"
                           for nr in [17, 18, 19, 20, 21]))

    print(f"\n  Results after {N} trials (DW[0]=+1, De=0 forced for rounds 1-15):")
    for nr in sorted(best_hw_by_rounds):
        coll_str = f" *** COLLISION at trial {collision_found[nr]}!" if nr in collision_found else ""
        print(f"    {nr} rounds: best hash diff HW = {best_hw_by_rounds[nr]}{coll_str}")


if __name__ == "__main__":
    random.seed(42)

    derive_ideal_path()
    search_21step_collision()
    nikolic_biryukov_pattern()
    full_de_forcing()

    print("\n" + "=" * 70)
    print("Step 17b Complete")
    print("=" * 70)
