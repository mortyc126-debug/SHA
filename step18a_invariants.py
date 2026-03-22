#!/usr/bin/env python3
"""Step 18a: Invariant hunting for SHA-256.

What our data whispers but we haven't listened to:

1. 97.6% AND cancellation = the round function is "almost" preserving
   some structure. What structure?

2. Kernel {d, h, f⊕g, a⊕b⊕c} = hidden symmetry group.
   These combinations are "invisible" to the nonlinear part.

3. Da - De = T2 - d = Sigma0(a) + Maj(a,b,c) - c
   This is a CONSERVATION LAW — independent of the message word.

Research questions:
A) Are there functions I(state) preserved by the round function?
B) What algebraic structure does the conservation law reveal?
C) Can we find approximate invariants that decay slowly?
D) Does the kernel generate a group action that commutes with rounds?
"""

import random
import numpy as np
from collections import defaultdict
MASK32 = 0xFFFFFFFF

def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
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

def sha256_step(state, W_t, K_t):
    a, b, c, d, e, f, g, h = state
    T1 = (h + Sigma1(e) + Ch(e, f, g) + K_t + W_t) & MASK32
    T2 = (Sigma0(a) + Maj(a, b, c)) & MASK32
    return ((T1 + T2) & MASK32, a, b, c, (d + T1) & MASK32, e, f, g)

# ============================================================
# Part 1: The conservation law Da-De = T2-d
# ============================================================
def conservation_law():
    """Explore the conservation law: Da_new - De_new = T2_diff - d_diff.

    This is EXACT and holds for ANY DW value.
    What does it mean algebraically?

    Let S = (a,b,c,d,e,f,g,h) be the state, S' the perturbed state.
    After one round:
      a_new = T1 + T2,  a'_new = T1' + T2'
      e_new = d + T1,   e'_new = d' + T1'

    Da_new = (T1' + T2') - (T1 + T2) = DT1 + DT2
    De_new = (d' + T1') - (d + T1) = Dd + DT1

    Da_new - De_new = DT2 - Dd

    DT2 = [Sigma0(a') + Maj(a',b',c')] - [Sigma0(a) + Maj(a,b,c)]

    This depends ONLY on (a,b,c) and (a',b',c') — NOT on W!

    So the quantity Q = Da - De satisfies:
    Q_new = DT2 - Dd
    which is determined by the state, not the message.

    This means: the EVOLUTION of Q is message-independent.
    Q is a "semi-invariant" — its dynamics are decoupled from W.
    """
    print("=" * 70)
    print("Part 1: Conservation law Q = Da - De")
    print("=" * 70)

    random.seed(42)
    N = 100000

    # Verify the conservation law
    violations = 0
    for _ in range(N):
        a = random.randint(0, MASK32)
        b = random.randint(0, MASK32)
        c = random.randint(0, MASK32)
        d = random.randint(0, MASK32)
        e = random.randint(0, MASK32)
        f = random.randint(0, MASK32)
        g = random.randint(0, MASK32)
        h = random.randint(0, MASK32)
        W = random.randint(0, MASK32)
        DW = random.randint(0, MASK32)

        s = (a, b, c, d, e, f, g, h)
        s2 = (a, b, c, d, e, f, g, h)  # Same start, different W
        s_new = sha256_step(s, W, K[0])
        s2_new = sha256_step(s2, (W + DW) & MASK32, K[0])

        Da_new = (s2_new[0] - s_new[0]) & MASK32
        De_new = (s2_new[4] - s_new[4]) & MASK32
        Q_new = (Da_new - De_new) & MASK32

        # Predicted: Q = DT2 - Dd = DT2 (since Dd = 0 when same state)
        T2 = (Sigma0(a) + Maj(a, b, c)) & MASK32
        # With same state, DT2 = 0 and Dd = 0, so Q = 0
        if Q_new != 0:
            violations += 1

    print(f"  Same state, different W: Q=0 always? violations={violations}/{N}")

    # Now with DIFFERENT states (differential scenario)
    print(f"\n  Different states (Da=1, De=1, others=0):")
    for _ in range(5):
        a = random.randint(0, MASK32)
        b = random.randint(0, MASK32)
        c = random.randint(0, MASK32)
        d = random.randint(0, MASK32)
        e = random.randint(0, MASK32)
        f = random.randint(0, MASK32)
        g = random.randint(0, MASK32)
        h = random.randint(0, MASK32)
        W = random.randint(0, MASK32)

        s = (a, b, c, d, e, f, g, h)
        s2 = ((a+1)&MASK32, b, c, d, (e+1)&MASK32, f, g, h)

        s_new = sha256_step(s, W, K[0])
        s2_new = sha256_step(s2, W, K[0])  # SAME W (no DW)

        Da = (s2_new[0] - s_new[0]) & MASK32
        De = (s2_new[4] - s_new[4]) & MASK32
        Q = (Da - De) & MASK32

        T2_1 = (Sigma0(a) + Maj(a, b, c)) & MASK32
        T2_2 = (Sigma0((a+1)&MASK32) + Maj((a+1)&MASK32, b, c)) & MASK32
        DT2 = (T2_2 - T2_1) & MASK32
        Dd = 0  # d is same

        predicted_Q = (DT2 - Dd) & MASK32

        print(f"    Da=0x{Da:08x}, De=0x{De:08x}, Q=0x{Q:08x}, "
              f"DT2=0x{DT2:08x}, predicted=0x{predicted_Q:08x}, "
              f"match={'✓' if Q == predicted_Q else '✗'}")

# ============================================================
# Part 2: Multi-round Q evolution (message-independent!)
# ============================================================
def q_evolution():
    """Track Q = Da - De over multiple rounds.
    Since Q's evolution is message-independent, it only depends on
    the INITIAL differential and the STATE values.

    Key question: does Q have attractors, cycles, or decay patterns?
    """
    print("\n" + "=" * 70)
    print("Part 2: Multi-round Q evolution")
    print("=" * 70)

    random.seed(42)
    N = 50000

    # Start with Da=1, De=1 (Q_0 = 0)
    # Track Q over rounds with different random messages
    q_hw_by_round = defaultdict(list)
    q_zero_by_round = defaultdict(int)

    for trial in range(N):
        a = random.randint(0, MASK32)
        b = random.randint(0, MASK32)
        c = random.randint(0, MASK32)
        d = random.randint(0, MASK32)
        e = random.randint(0, MASK32)
        f = random.randint(0, MASK32)
        g = random.randint(0, MASK32)
        h = random.randint(0, MASK32)

        s1 = (a, b, c, d, e, f, g, h)
        s2 = ((a+1)&MASK32, b, c, d, (e+1)&MASK32, f, g, h)

        for t in range(20):
            W = random.randint(0, MASK32)
            s1 = sha256_step(s1, W, K[t])
            s2 = sha256_step(s2, W, K[t])  # SAME W

            Da = (s2[0] - s1[0]) & MASK32
            De = (s2[4] - s1[4]) & MASK32
            Q = (Da - De) & MASK32

            q_hw_by_round[t].append(hw(Q))
            if Q == 0:
                q_zero_by_round[t] += 1

    print(f"\n  Q = Da - De evolution (Da0=1, De0=1, Q0=0, same W):")
    print(f"  {'Round':>5} | {'Avg Q HW':>9} | {'Q=0 count':>10} | {'P(Q=0)':>10}")
    print("  " + "-" * 45)
    for t in range(20):
        hws = q_hw_by_round[t]
        avg = sum(hws) / len(hws)
        zeros = q_zero_by_round[t]
        p = zeros / N
        print(f"  {t:5d} | {avg:9.2f} | {zeros:10d} | {p:10.6f}")

# ============================================================
# Part 3: Linear invariants modulo 2^32
# ============================================================
def linear_invariants():
    """Search for linear combinations of state registers that are
    preserved (modulo 2^32) by the round function.

    I(state) = sum_i c_i * state[i]  (mod 2^32)
    where c_i are constants.

    For I to be invariant: I(F(state, W)) = I(state) for all W.
    Since W appears linearly in T1, and T1 affects both a and e:
      I(new) = c_a*(T1+T2) + c_b*a + c_c*b + c_d*c + c_e*(d+T1) + c_f*e + c_g*f + c_h*g
    For this to equal I(old) = c_a*a + c_b*b + ... + c_h*h for all W:
      The W-dependent part must vanish: (c_a + c_e)*T1 = (c_a + c_e)*(h + Sigma1(e) + Ch(e,f,g) + K + W)
      For this to cancel for all W: c_a + c_e = 0, i.e., c_e = -c_a
    """
    print("\n" + "=" * 70)
    print("Part 3: Linear invariants (mod 2^32)")
    print("=" * 70)

    # For I = c_a*a + c_b*b + ... + c_h*h to be invariant:
    # Constraint 1: c_a + c_e = 0 (W-independence)
    # Remaining: c_a*(T2) + c_b*a + c_c*b + c_d*c + c_e*d + c_f*e + c_g*f + c_h*g
    #          = c_a*a + c_b*b + c_c*b + ... (original)
    #
    # Shift register: b_new = a, c_new = b, d_new = c, f_new = e, g_new = f, h_new = g
    # So: c_b*a_old + c_c*b_old + c_d*c_old + c_f*e_old + c_g*f_old + c_h*g_old
    #   should match parts of the original.
    #
    # From shift: b_new = a → coefficient of a in new = c_b
    #             c_new = b → coefficient of b in new = c_c
    #             d_new = c → coefficient of c in new = c_d
    #             f_new = e → coefficient of e in new = c_f
    #             g_new = f → coefficient of f in new = c_g
    #             h_new = g → coefficient of g in new = c_h
    #
    # For invariance (ignoring nonlinear T2):
    # c_b = c_a (from a → b shift)
    # c_c = c_b → c_c = c_a
    # c_d = c_c → c_d = c_a
    # c_f = c_e = -c_a
    # c_g = c_f = -c_a
    # c_h = c_g = -c_a
    #
    # So: I = c_a*(a + b + c + d - e - f - g - h)
    #
    # Check: I_new = c_a*[(T1+T2) + a + b + c - (d+T1) - e - f - g]
    #       = c_a*[T2 + a + b + c - d - e - f - g]
    #       = c_a*[(a+b+c+d-e-f-g-h) + T2 - d - d + h]  ... getting complicated

    # Let me verify empirically
    print("\n  Testing I = a + b + c + d - e - f - g - h:")
    random.seed(42)
    N = 10000

    for name, coeffs in [
        ("a+b+c+d-e-f-g-h", [1,1,1,1,-1,-1,-1,-1]),
        ("a+b+c+d+e+f+g+h", [1,1,1,1,1,1,1,1]),
        ("a-e", [1,0,0,0,-1,0,0,0]),
        ("a+b+c-e-f-g", [1,1,1,0,-1,-1,-1,0]),
        ("a+b+c+d", [1,1,1,1,0,0,0,0]),
        ("e+f+g+h", [0,0,0,0,1,1,1,1]),
        ("a-b+c-d+e-f+g-h", [1,-1,1,-1,1,-1,1,-1]),
    ]:
        diffs = []
        for _ in range(N):
            s = tuple(random.randint(0, MASK32) for _ in range(8))
            W = random.randint(0, MASK32)
            s_new = sha256_step(s, W, K[0])

            I_old = sum(c * s[i] for i, c in enumerate(coeffs))
            I_new = sum(c * s_new[i] for i, c in enumerate(coeffs))
            diff = (I_new - I_old) & MASK32
            diffs.append(diff)

        # Is the diff constant (doesn't depend on W)?
        unique_diffs = len(set(diffs))
        zero_count = diffs.count(0)
        print(f"  {name:25s}: unique diffs = {unique_diffs}, zeros = {zero_count}/{N}")

# ============================================================
# Part 4: The kernel as a symmetry group
# ============================================================
def kernel_symmetry():
    """The kernel {d, h, f⊕g, a⊕b⊕c} means these combinations
    don't participate in the quadratic part.

    Question: if we perturb along kernel directions, what happens?
    Does the perturbation propagate differently?

    Perturbation along d: change d only.
    Perturbation along h: change h only.
    Perturbation along f⊕g: flip same bits in f and g.
    Perturbation along a⊕b⊕c: flip same bits in a, b, and c.
    """
    print("\n" + "=" * 70)
    print("Part 4: Kernel direction perturbations")
    print("=" * 70)

    random.seed(42)
    N = 50000

    perturbations = {
        "Dd=1 (kernel)": lambda s: (s[0],s[1],s[2],(s[3]+1)&MASK32,s[4],s[5],s[6],s[7]),
        "Dh=1 (kernel)": lambda s: (s[0],s[1],s[2],s[3],s[4],s[5],s[6],(s[7]+1)&MASK32),
        "Df=Dg=XOR1 (kernel)": lambda s: (s[0],s[1],s[2],s[3],s[4],s[5]^1,s[6]^1,s[7]),
        "Da=Db=Dc=XOR1 (kernel)": lambda s: (s[0]^1,s[1]^1,s[2]^1,s[3],s[4],s[5],s[6],s[7]),
        "Da=1 (NOT kernel)": lambda s: ((s[0]+1)&MASK32,s[1],s[2],s[3],s[4],s[5],s[6],s[7]),
        "De=1 (NOT kernel)": lambda s: (s[0],s[1],s[2],s[3],(s[4]+1)&MASK32,s[5],s[6],s[7]),
    }

    for name, perturb in perturbations.items():
        diff_hw_by_round = []
        for _ in range(N):
            s = tuple(random.randint(0, MASK32) for _ in range(8))
            s2 = perturb(s)

            for t in range(12):
                W = random.randint(0, MASK32)
                s = sha256_step(s, W, K[t])
                s2 = sha256_step(s2, W, K[t])

            total_hw = sum(hw((s2[i] - s[i]) & MASK32) for i in range(8))
            diff_hw_by_round.append(total_hw)

        avg = sum(diff_hw_by_round) / N
        zeros = diff_hw_by_round.count(0)
        print(f"  {name:30s}: avg diff HW after 12 rounds = {avg:6.1f}, zeros = {zeros}")

# ============================================================
# Part 5: Sum conservation — does a+b+c+d+e+f+g+h have structure?
# ============================================================
def sum_conservation():
    """The total sum S = a+b+c+d+e+f+g+h changes by T2 + T1 per round.

    S_new = (T1+T2) + a + b + c + (d+T1) + e + f + g
          = 2*T1 + T2 + (a+b+c+d+e+f+g)
          = 2*T1 + T2 + S - h

    So: DS = S_new - S = 2*T1 + T2 - h

    Note: T1 = h + Sigma1(e) + Ch(e,f,g) + K + W
    So: 2*T1 = 2h + 2*Sigma1(e) + 2*Ch(e,f,g) + 2K + 2W
    DS = 2h + 2*Sigma1(e) + 2*Ch(e,f,g) + 2K + 2W + T2 - h
       = h + 2*Sigma1(e) + 2*Ch(e,f,g) + 2K + 2W + Sigma0(a) + Maj(a,b,c)

    This depends on W → not invariant.

    But what about S mod 2? Or S mod some other modulus?
    """
    print("\n" + "=" * 70)
    print("Part 5: Sum and parity conservation")
    print("=" * 70)

    random.seed(42)
    N = 50000

    # Check various modular sums
    for mod in [2, 4, 8, 16]:
        for name, func in [
            ("sum(all)", lambda s: sum(s)),
            ("a+b+c+d", lambda s: s[0]+s[1]+s[2]+s[3]),
            ("e+f+g+h", lambda s: s[4]+s[5]+s[6]+s[7]),
            ("a-e+b-f+c-g+d-h", lambda s: s[0]-s[4]+s[1]-s[5]+s[2]-s[6]+s[3]-s[7]),
            ("a*e XOR b*f XOR c*g XOR d*h",
             lambda s: (s[0]&s[4])^(s[1]&s[5])^(s[2]&s[6])^(s[3]&s[7])),
        ]:
            preserved = 0
            for _ in range(N):
                s = tuple(random.randint(0, MASK32) for _ in range(8))
                W = random.randint(0, MASK32)
                s_new = sha256_step(s, W, K[0])

                val_old = func(s) % mod
                val_new = func(s_new) % mod
                if val_old == val_new:
                    preserved += 1

            p = preserved / N
            expected = 1/mod  # Random expectation
            if abs(p - expected) > 0.01:
                marker = " ← ANOMALY!" if abs(p - expected) > 0.05 else " ← slight"
            else:
                marker = ""
            print(f"  mod {mod:2d} {name:35s}: P(preserved) = {p:.4f} "
                  f"(random = {expected:.4f}){marker}")

# ============================================================
# Part 6: XOR parity invariants (GF(2))
# ============================================================
def xor_invariants():
    """Search for XOR (GF2) linear invariants.

    I(state) = XOR of selected bits from registers.
    For invariance: I(F(state, W)) = I(state) for all state, W.

    Since W appears linearly (modulo carries), XOR invariants are rare.
    But maybe some exist due to the specific Sigma/Ch/Maj structure.
    """
    print("\n" + "=" * 70)
    print("Part 6: GF(2) linear invariant search")
    print("=" * 70)

    random.seed(42)
    N = 10000

    # For each pair of bit positions (reg_i, bit_j) and (reg_k, bit_l),
    # check if XOR is invariant.

    # First: single-bit invariants (any single bit preserved?)
    print("\n  Single-bit invariants:")
    single_bit_scores = {}
    for reg in range(8):
        for bit in range(32):
            preserved = 0
            for _ in range(N):
                s = tuple(random.randint(0, MASK32) for _ in range(8))
                W = random.randint(0, MASK32)
                s_new = sha256_step(s, W, K[0])

                b_old = (s[reg] >> bit) & 1
                b_new = (s_new[reg] >> bit) & 1
                if b_old == b_new:
                    preserved += 1

            score = abs(preserved/N - 0.5)
            if score > 0.02:
                single_bit_scores[(reg, bit)] = preserved/N

    if single_bit_scores:
        top = sorted(single_bit_scores.items(), key=lambda x: -abs(x[1]-0.5))[:10]
        for (reg, bit), p in top:
            print(f"    reg[{reg}] bit {bit:2d}: P(preserved) = {p:.4f}")
    else:
        print("    No single-bit invariants found (all within 0.02 of 0.5)")

    # XOR pairs: bit_i XOR bit_j
    print("\n  Two-bit XOR invariants (searching for anomalies):")
    best_pairs = []
    reg_names = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']

    # Sample: check specific interesting pairs
    interesting_pairs = [
        # Kernel-inspired
        ((5, 0), (6, 0)),  # f[0] XOR g[0]
        ((0, 0), (1, 0)),  # a[0] XOR b[0]
        ((0, 0), (4, 0)),  # a[0] XOR e[0]
        ((3, 0), (7, 0)),  # d[0] XOR h[0]
        # Sigma-inspired
        ((4, 6), (4, 11)),  # e[6] XOR e[11]
        ((0, 2), (0, 13)),  # a[2] XOR a[13]
    ]

    for (r1, b1), (r2, b2) in interesting_pairs:
        preserved = 0
        for _ in range(N):
            s = tuple(random.randint(0, MASK32) for _ in range(8))
            W = random.randint(0, MASK32)
            s_new = sha256_step(s, W, K[0])

            xor_old = ((s[r1] >> b1) ^ (s[r2] >> b2)) & 1
            xor_new = ((s_new[r1] >> b1) ^ (s_new[r2] >> b2)) & 1
            if xor_old == xor_new:
                preserved += 1

        p = preserved / N
        bias = abs(p - 0.5)
        marker = " ← BIAS" if bias > 0.02 else ""
        print(f"    {reg_names[r1]}[{b1:2d}] XOR {reg_names[r2]}[{b2:2d}]: "
              f"P={p:.4f}{marker}")

# ============================================================
# Part 7: Deep structure — correlation between rounds
# ============================================================
def inter_round_correlation():
    """Look for correlations between state values across rounds.

    If there's a hidden invariant, states at different rounds should
    be correlated in ways not explained by the shift register alone.
    """
    print("\n" + "=" * 70)
    print("Part 7: Inter-round correlations")
    print("=" * 70)

    random.seed(42)
    N = 100000

    # Track specific quantities across rounds
    print("\n  a[t] XOR e[t] parity correlation between rounds:")
    corr_matrix = np.zeros((8, 8))

    parities = []
    for _ in range(N):
        s = tuple(random.randint(0, MASK32) for _ in range(8))
        round_parities = []
        for t in range(8):
            W = random.randint(0, MASK32)
            s = sha256_step(s, W, K[t])
            # Parity of a XOR e
            p = hw(s[0] ^ s[4]) % 2
            round_parities.append(p)
        parities.append(round_parities)

    parities = np.array(parities)
    for i in range(8):
        for j in range(8):
            corr_matrix[i][j] = np.corrcoef(parities[:, i], parities[:, j])[0, 1]

    print("  Correlation matrix (a XOR e parity between rounds):")
    print("       " + " ".join(f"  r{j}" for j in range(8)))
    for i in range(8):
        row = " ".join(f"{corr_matrix[i][j]:+.2f}" for j in range(8))
        print(f"  r{i}:  {row}")


if __name__ == "__main__":
    conservation_law()
    q_evolution()
    linear_invariants()
    kernel_symmetry()
    sum_conservation()
    xor_invariants()
    inter_round_correlation()

    print("\n" + "=" * 70)
    print("Step 18a Complete")
    print("=" * 70)
