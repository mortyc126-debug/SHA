#!/usr/bin/env python3
"""Step 19a: Algebraic dependencies between consecutive barriers.

KEY QUESTION: Are De17=0, De18=0, De19=0 ... independent equations,
or do they share algebraic structure that can be exploited?

What we know:
1. De_new = Dd + T1 = c_old + h + Sigma1(e) + Ch(e,f,g) + K + W
2. Da_new - De_new = T2 - d = Sigma0(a) + Maj(a,b,c) - c  (conservation law)
3. Shift register: De_{t+1} depends on state at round t, which depends on De_t

The barriers are NOT independent — they're RECURSIVELY coupled.
De_{t+1} is a function of the state AFTER De_t has been set.

Research plan:
A) Express De_{t+1} as a function of De_t and other variables
B) Find the "transfer function" between consecutive barriers
C) Look for algebraic cancellations in the transfer function
D) Test if solving De_t=0 constrains De_{t+1} in a useful way
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
# Part 1: Express De_{t+1} in terms of state at round t
# ============================================================
def barrier_transfer_function():
    """Derive: given De_t, what is De_{t+1}?

    State at round t (after computing De_t):
      (a_t, b_t, c_t, d_t, e_t, f_t, g_t, h_t)
    For the perturbed path:
      (a_t', b_t', c_t', d_t', e_t', f_t', g_t', h_t')

    De_t = e_t' - e_t.

    At round t+1:
      e_{t+1} = d_t + T1_{t+1}  where T1 = h_t + Sigma1(e_t) + Ch(e_t,f_t,g_t) + K_{t+1} + W_{t+1}
      e_{t+1}' = d_t' + T1'_{t+1}

    De_{t+1} = (d_t' - d_t) + (T1' - T1)
             = Dd_t + DT1_{t+1}

    Now: Dd_t = c_{t-1}' - c_{t-1} = Dc_{t-1} = Db_{t-2} = Da_{t-3}
    (from the shift register: d_t = c_{t-1} = b_{t-2} = a_{t-3})

    And: DT1 = Dh_t + D[Sigma1(e_t)] + D[Ch(e_t,f_t,g_t)] + DW_{t+1}

    So: De_{t+1} = Da_{t-3} + Dh_t + D[Sigma1] + D[Ch] + DW_{t+1}

    Using conservation: Da_t = De_t + DT2_t - Dd_t
    And shift register backwards:
      Da_{t-3} = De_{t-3} + DT2_{t-3} - Dd_{t-3}

    This creates a RECURRENCE relating De values across rounds!
    """
    print("=" * 70)
    print("Part 1: Barrier transfer function")
    print("=" * 70)

    print("""
  De_{t+1} = Dd_t + Dh_t + D[Sigma1(e_t)] + D[Ch(e_t,f_t,g_t)] + DW_{t+1}

  Where:
    Dd_t = Da_{t-3}      (shift register, 3 rounds delay)
    Dh_t = De_{t-3}      (shift register, 3 rounds delay)
    D[Sigma1] = Sigma1(e_t') - Sigma1(e_t)  (depends on De_t and carry structure)
    D[Ch] = Ch(e_t',f_t',g_t') - Ch(e_t,f_t,g_t)  (depends on De_t, Df_t, Dg_t)

  Substituting shift register:
    Df_t = De_{t-1}
    Dg_t = De_{t-2}

  So: De_{t+1} = F(Da_{t-3}, De_{t-3}, De_t, De_{t-1}, De_{t-2}, DW_{t+1}, state)

  And from conservation law: Da_{t-3} = De_{t-3} + [Sigma0+Maj diff at t-3] - Dc_{t-3}

  This is a RECURRENCE of depth 4 in De!
  De_{t+1} depends on De_t, De_{t-1}, De_{t-2}, De_{t-3}
  """)

    # Verify empirically: measure correlation between consecutive De values
    print("  Empirical correlation between consecutive barriers:")
    random.seed(42)
    N = 200000
    DW0 = 0x80000000

    de_pairs = defaultdict(list)

    for _ in range(N):
        W = [random.randint(0, MASK32) for _ in range(16)]
        W2 = list(W)
        W2[0] ^= DW0

        Wexp = expand_message(W)
        W2exp = expand_message(W2)

        s1 = list(IV)
        s2 = list(IV)
        de_history = []

        for t in range(25):
            s1 = list(sha256_step(s1, Wexp[t], K[t]))
            s2 = list(sha256_step(s2, W2exp[t], K[t]))
            de = (s2[4] - s1[4]) & MASK32
            de_history.append(de)

        for t in range(16, 23):
            de_pairs[t].append((hw(de_history[t]), hw(de_history[t+1])))

    print(f"\n  Correlation(HW(De_t), HW(De_{t+1})):")
    import numpy as np
    for t in range(16, 23):
        pairs = de_pairs[t]
        x = [p[0] for p in pairs]
        y = [p[1] for p in pairs]
        corr = np.corrcoef(x, y)[0, 1]
        print(f"    rounds {t}-{t+1}: r = {corr:.4f}")

# ============================================================
# Part 2: Conditional barrier probability
# ============================================================
def conditional_barriers():
    """If De_t = 0, what is P(De_{t+1} = 0)?

    If barriers were independent: P = 2^-32.
    If they're coupled: P might be higher or lower.

    We need massive sampling because P is tiny.
    Alternative: measure P(HW(De_{t+1}) | De_t has specific HW).
    """
    print("\n" + "=" * 70)
    print("Part 2: Conditional barrier probabilities")
    print("=" * 70)

    random.seed(42)
    N = 5000000  # 5M trials
    DW0 = 0x80000000

    # Track: when De_t has low HW, what HW does De_{t+1} have?
    conditional_hw = defaultdict(lambda: defaultdict(list))

    for _ in range(N):
        W = [random.randint(0, MASK32) for _ in range(16)]
        W2 = list(W)
        W2[0] ^= DW0

        Wexp = expand_message(W)
        W2exp = expand_message(W2)

        s1 = list(IV)
        s2 = list(IV)

        de_prev_hw = None
        for t in range(25):
            s1 = list(sha256_step(s1, Wexp[t], K[t]))
            s2 = list(sha256_step(s2, W2exp[t], K[t]))
            de = (s2[4] - s1[4]) & MASK32
            de_hw = hw(de)

            if t >= 16 and de_prev_hw is not None:
                if de_prev_hw <= 5:
                    conditional_hw[t][de_prev_hw].append(de_hw)

            de_prev_hw = de_hw

    print(f"\n  HW(De_{t+1}) | HW(De_t) = k  (for rounds 16-23):")
    for t in range(17, 24):
        for k in range(6):
            data = conditional_hw[t].get(k, [])
            if len(data) > 10:
                avg = sum(data) / len(data)
                low_count = sum(1 for x in data if x <= 5)
                print(f"    Round {t}, given HW(De_{t-1})={k}: "
                      f"avg HW(De_t) = {avg:.1f}, "
                      f"P(HW≤5) = {low_count}/{len(data)} "
                      f"(n={len(data)})")

# ============================================================
# Part 3: The conservation law as a coupling mechanism
# ============================================================
def conservation_coupling():
    """The conservation law Da-De = T2-d creates a CONSTRAINT
    between the a-path and e-path barriers.

    If De_t = 0, then Da_t = DT2_t - Dd_t.
    This Da_t becomes Dd_{t+3} via the shift register.
    Which feeds into De_{t+4} = Dd_{t+3} + DT1_{t+4}.

    So: De_t = 0 → Da_t = specific value → Dd_{t+3} = Da_t →
        De_{t+4} = Da_t + DT1_{t+4}

    The conservation law creates a 4-ROUND ECHO:
    Setting De_t = 0 constrains De_{t+4} through the Da path!

    Question: can we use this echo constructively?
    If we solve De_17 = 0, does it help with De_21?
    """
    print("\n" + "=" * 70)
    print("Part 3: Conservation law 4-round echo")
    print("=" * 70)

    random.seed(42)
    N = 200000
    DW0 = 0x80000000

    # When De_17 is low, track De_21 (the 4-round echo)
    echo_data = defaultdict(list)

    for _ in range(N):
        W = [random.randint(0, MASK32) for _ in range(16)]
        W2 = list(W)
        W2[0] ^= DW0

        Wexp = expand_message(W)
        W2exp = expand_message(W2)

        s1 = list(IV)
        s2 = list(IV)
        de_history = []
        da_history = []

        for t in range(25):
            s1 = list(sha256_step(s1, Wexp[t], K[t]))
            s2 = list(sha256_step(s2, W2exp[t], K[t]))
            de = (s2[4] - s1[4]) & MASK32
            da = (s2[0] - s1[0]) & MASK32
            de_history.append(de)
            da_history.append(da)

        # Track echo: De_t → De_{t+4}
        for t in range(16, 21):
            hw_de_t = hw(de_history[t])
            hw_de_t4 = hw(de_history[t+4]) if t+4 < len(de_history) else None
            hw_da_t = hw(da_history[t])
            if hw_de_t4 is not None:
                echo_data[t].append((hw_de_t, hw_de_t4, hw_da_t))

    print(f"\n  4-round echo: correlation between De_t and De_{{t+4}}:")
    import numpy as np
    for t in range(16, 21):
        data = echo_data[t]
        de_t = [d[0] for d in data]
        de_t4 = [d[1] for d in data]
        da_t = [d[2] for d in data]
        corr_de = np.corrcoef(de_t, de_t4)[0, 1]
        corr_da = np.corrcoef(de_t, da_t)[0, 1]
        print(f"    De_{t} → De_{t+4}: r = {corr_de:.4f}")
        print(f"    De_{t} → Da_{t}:  r = {corr_da:.4f}")

    # Critical: when De_t is low, is De_{t+4} also lower than random?
    print(f"\n  Conditional: avg HW(De_{{t+4}}) when HW(De_t) is low:")
    for t in range(16, 21):
        data = echo_data[t]
        for threshold in [5, 8, 10]:
            low_de_t = [(d[1], d[2]) for d in data if d[0] <= threshold]
            if low_de_t:
                avg_de_t4 = sum(d[0] for d in low_de_t) / len(low_de_t)
                avg_da_t = sum(d[1] for d in low_de_t) / len(low_de_t)
                print(f"    t={t}, HW(De_t)≤{threshold}: "
                      f"avg HW(De_{{t+4}})={avg_de_t4:.1f}, "
                      f"avg HW(Da_t)={avg_da_t:.1f} "
                      f"(n={len(low_de_t)})")

# ============================================================
# Part 4: Joint barrier equation — algebraic form
# ============================================================
def joint_barrier_equation():
    """Express (De_t, De_{t+1}) as a joint equation.

    De_t = e_t' - e_t = 0     (barrier at t)
    De_{t+1} = e_{t+1}' - e_{t+1} = 0  (barrier at t+1)

    De_{t+1} = Dd_t + DT1_{t+1}
             = Da_{t-3} + [Dh_t + DSigma1(e_t) + DCh_t + DW_{t+1}]

    If De_t = 0, then: e_t' = e_t (same e register)
    This simplifies DSigma1 and DCh at round t+1!

    DSigma1 at t+1: uses e_{t+1} and e_{t+1}'.
    But e_{t+1} = d_t + T1, e_{t+1}' = d_t' + T1'.
    De_{t+1} = Dd_t + DT1 = what we're trying to solve.

    Wait, DSigma1 at round t+1 uses the state at round t+1's inputs:
    e at round t+1 is e_t (from shift), and similarly e' is e_t'.
    Actually no: in the round function at step t+1:
    - e_input = e_t (from previous state)
    - f_input = f_t = e_{t-1}
    - g_input = g_t = e_{t-2}

    So when De_t = 0: the e input at round t+1 is the SAME for both.
    DSigma1(e_t) = 0!
    And DCh depends on De_t=0, Df_t=De_{t-1}, Dg_t=De_{t-2}.

    THIS IS HUGE: When De_t = 0, the Sigma1 differential at round t+1 VANISHES!
    """
    print("\n" + "=" * 70)
    print("Part 4: Joint barrier simplification")
    print("=" * 70)

    print("""
  THEOREM: When De_t = 0:
    1. DSigma1 at round t+1 = 0  (because e_t' = e_t)
    2. DCh at round t+1 depends only on De_{t-1} and De_{t-2}
       (via f_t = e_{t-1} and g_t = e_{t-2})
    3. Dh_t = Dg_{t-1} = Df_{t-2} = De_{t-3}

  So: De_{t+1} = Da_{t-3} + De_{t-3} + DCh(De_{t-1}, De_{t-2}) + DW_{t+1}

  If additionally De_{t-1} = De_{t-2} = De_{t-3} = 0 (Wang chain):
    DCh = 0, De_{t-3} = 0
    De_{t+1} = Da_{t-3} + DW_{t+1}

  And from conservation: Da_{t-3} = DT2_{t-3}  (when De_{t-3}=0 and Dd_{t-3}=0)

  So: De_{t+1} = DT2_{t-3} + DW_{t+1}

  DT2_{t-3} is DETERMINED by the state at round t-3.
  DW_{t+1} is determined by the message schedule.

  CONCLUSION: Under the Wang chain condition (De=0 for rounds 1..t),
  each new barrier is a SINGLE 32-bit equation:
    DW_{t+1} = -DT2_{t-3}

  This is exactly what we found before — each barrier costs 2^32.
  But NOW we can see WHY and look for structure in DT2.
  """)

    # Verify: compute DT2 at various rounds
    random.seed(42)
    N = 100000

    # Under Wang chain (De=0 forced), compute DT2 at each round
    dt2_hw_stats = defaultdict(list)

    for _ in range(N):
        W = [random.randint(0, MASK32) for _ in range(16)]
        W2 = list(W)
        W2[0] ^= 0x80000000

        s1 = list(IV)
        s2 = list(IV)

        # Round 0
        s1 = list(sha256_step(s1, W[0], K[0]))
        s2 = list(sha256_step(s2, W2[0], K[0]))

        # Rounds 1-15: force De=0
        for t in range(1, 16):
            T1_p1 = (s1[7] + Sigma1(s1[4]) + Ch(s1[4], s1[5], s1[6]) + K[t]) & MASK32
            T1_p2 = (s2[7] + Sigma1(s2[4]) + Ch(s2[4], s2[5], s2[6]) + K[t]) & MASK32
            dw = ((s1[3] + T1_p1) - (s2[3] + T1_p2)) & MASK32
            W2[t] = (W[t] + dw) & MASK32
            s1 = list(sha256_step(s1, W[t], K[t]))
            s2 = list(sha256_step(s2, W2[t], K[t]))

        # Compute DT2 at each round (under Wang chain)
        # DT2 = [Sigma0(a') + Maj(a',b',c')] - [Sigma0(a) + Maj(a,b,c)]
        T2_1 = (Sigma0(s1[0]) + Maj(s1[0], s1[1], s1[2])) & MASK32
        T2_2 = (Sigma0(s2[0]) + Maj(s2[0], s2[1], s2[2])) & MASK32
        DT2 = (T2_2 - T2_1) & MASK32
        dt2_hw_stats['round15'].append(hw(DT2))

    avg = sum(dt2_hw_stats['round15']) / N
    zeros = dt2_hw_stats['round15'].count(0)
    print(f"  DT2 at round 15 (under Wang chain):")
    print(f"    Avg HW = {avg:.1f}")
    print(f"    DT2 = 0: {zeros}/{N}")
    print(f"    => The barrier equation is DW[t+1] = -DT2")
    print(f"    => P(DW[t+1] = -DT2) depends on message schedule freedom")

# ============================================================
# Part 5: DT2 structure — is it exploitable?
# ============================================================
def dt2_structure():
    """DT2 = D[Sigma0(a)] + D[Maj(a,b,c)].

    Under Wang chain, Da is accumulated and large (~16 HW).
    But what about the ALGEBRAIC structure of DT2?

    DT2 depends on (a, a', b, b', c, c') = (a, a+Da, b, b+Db, c, c+Dc)
    where Da, Db, Dc are the accumulated differences.

    Key: Db = Da from 1 round ago, Dc = Da from 2 rounds ago.
    So all of Da, Db, Dc are COUPLED through the conservation law.

    Is there a choice of initial message that makes DT2 small or structured?
    """
    print("\n" + "=" * 70)
    print("Part 5: DT2 structure analysis")
    print("=" * 70)

    random.seed(42)
    N = 200000

    # Under Wang chain, track Da and DT2 evolution
    da_hw_by_round = defaultdict(list)
    dt2_hw_by_round = defaultdict(list)

    # Also track: what fraction of DT2 comes from DSigma0 vs DMaj?
    dsig0_hw_by_round = defaultdict(list)
    dmaj_hw_by_round = defaultdict(list)

    for _ in range(N):
        W = [random.randint(0, MASK32) for _ in range(16)]
        W2 = list(W)
        W2[0] ^= 0x80000000

        s1 = list(IV)
        s2 = list(IV)

        s1 = list(sha256_step(s1, W[0], K[0]))
        s2 = list(sha256_step(s2, W2[0], K[0]))

        for t in range(1, 16):
            # Record Da, DT2 components BEFORE forcing De=0
            Da = (s2[0] - s1[0]) & MASK32
            DSigma0 = (Sigma0(s2[0]) - Sigma0(s1[0])) & MASK32
            DMaj = (Maj(s2[0], s2[1], s2[2]) - Maj(s1[0], s1[1], s1[2])) & MASK32
            DT2 = (DSigma0 + DMaj) & MASK32

            da_hw_by_round[t].append(hw(Da))
            dt2_hw_by_round[t].append(hw(DT2))
            dsig0_hw_by_round[t].append(hw(DSigma0))
            dmaj_hw_by_round[t].append(hw(DMaj))

            # Force De=0
            T1_p1 = (s1[7] + Sigma1(s1[4]) + Ch(s1[4], s1[5], s1[6]) + K[t]) & MASK32
            T1_p2 = (s2[7] + Sigma1(s2[4]) + Ch(s2[4], s2[5], s2[6]) + K[t]) & MASK32
            dw = ((s1[3] + T1_p1) - (s2[3] + T1_p2)) & MASK32
            W2[t] = (W[t] + dw) & MASK32
            s1 = list(sha256_step(s1, W[t], K[t]))
            s2 = list(sha256_step(s2, W2[t], K[t]))

    print(f"\n  Under Wang chain (De=0 forced), differential evolution:")
    print(f"  {'Round':>5} | {'Da HW':>6} | {'DT2 HW':>7} | {'DSig0 HW':>8} | {'DMaj HW':>8} | {'DT2=0':>6}")
    print("  " + "-" * 55)
    for t in range(1, 16):
        da_avg = sum(da_hw_by_round[t]) / N
        dt2_avg = sum(dt2_hw_by_round[t]) / N
        ds0_avg = sum(dsig0_hw_by_round[t]) / N
        dm_avg = sum(dmaj_hw_by_round[t]) / N
        dt2_zeros = dt2_hw_by_round[t].count(0)
        print(f"  {t:5d} | {da_avg:6.1f} | {dt2_avg:7.1f} | {ds0_avg:8.1f} | {dm_avg:8.1f} | {dt2_zeros:6d}")

    # KEY: is DT2 correlated across rounds?
    print(f"\n  Correlation of DT2 HW across consecutive rounds:")
    import numpy as np
    for t in range(2, 15):
        corr = np.corrcoef(dt2_hw_by_round[t], dt2_hw_by_round[t+1])[0, 1]
        print(f"    DT2_{t} → DT2_{t+1}: r = {corr:.4f}")


if __name__ == "__main__":
    barrier_transfer_function()
    conditional_barriers()
    conservation_coupling()
    joint_barrier_equation()
    dt2_structure()

    print("\n" + "=" * 70)
    print("Step 19a Complete")
    print("=" * 70)
