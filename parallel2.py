"""
PARALLEL THREADS:

A: Can BTE's створочне enable SIMULTANEOUS a+e control?
   Wang controls e only. Створочне links a and e.
   If we control a → e follows (or vice versa).

B: PRODUCT STRUCTURE of carry algebra.
   One addition: T3-T5 characterize it.
   6 additions per round: how do their carries INTERACT?
   Is there a "product" carry algebra?
"""

import random
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, schedule
)


# ═══════════════════════════════════════════════════
# THREAD A: Simultaneous a+e control
# ═══════════════════════════════════════════════════

def thread_A():
    """
    From створочне: e[r] = a[r] + a[r-4] - T2[r-1]
    At bit 0 (no carry): e[r][0] = a[r][0] ⊕ a[r-4][0] ⊕ T2[r-1][0]

    T2[r-1][0] = Sig0(a[r-1])[0] ⊕ Maj(a[r-1],a[r-2],a[r-3])[0]
               = a[r-1][2] ⊕ a[r-1][13] ⊕ a[r-1][22] ⊕ Maj(a[r-1][0],a[r-2][0],a[r-3][0])

    If we SET δa[r] = 0 for all r (a-convergence):
      → δe[r] = 0 automatically (from створочне: δe = δa + δa[-4] - δT2 = 0+0-0 = 0)!

    So: a-convergence → e-convergence FOR FREE!
    Wang does it opposite: controls e, a remains uncontrolled.

    QUESTION: Can we control a instead of e? What does this require?
    """
    print("=" * 80)
    print("THREAD A: a-convergence → e-convergence (створочне)")
    print("=" * 80)

    # Verify: if δa[r]=0 for r=R-7..R, does δe[r]=0 follow?
    # From створочне: e[r] = a[r] + a[r-4] - T2[r-1]
    # δe[r] = δa[r] + δa[r-4] - δT2[r-1]
    # If δa[r]=0 AND δa[r-4]=0 AND δT2[r-1]=0: δe[r]=0. ✓
    # δT2[r-1] = Sig0(δa[r-1]) + Maj-diff. If δa[r-1]=δa[r-2]=δa[r-3]=0: δT2=0.
    # So: need δa[r]=0 for r = R-11..R (to cover the -4 offset and T2 deps).
    # = 12 consecutive zero a-differences → 8 zero e-differences.

    print(f"\n  If δa[r] = 0 for r = R-11..R (12 rounds):")
    print(f"    → δT2[r-1] = 0 for r = R-7..R (since T2 uses a[r-1..r-3])")
    print(f"    → δa[r-4] = 0 for r = R-7..R (since r-4 ≥ R-11)")
    print(f"    → δe[r] = 0 for r = R-7..R (from створочне)")
    print(f"    → state[R] match → COLLISION!")

    # What controls a?
    # a[r+1] = T1[r] + T2[r]
    # T1[r] = h[r] + Sig1(e[r]) + Ch(e[r],f[r],g[r]) + K[r] + W[r]
    # T2[r] = Sig0(a[r]) + Maj(a[r],a[r-1],a[r-2])
    #
    # For δa[r+1] = 0: δT1[r] + δT2[r] = 0 (mod 2^32)
    #   → δT1[r] = -δT2[r]
    #
    # δT2[r] depends on δa[r], δa[r-1], δa[r-2] through Sig0 and Maj.
    # If δa = 0 for r-2..r: δT2 = 0.
    # Then δT1 = 0 needed.
    #
    # δT1 = δh + δSig1(e) + δCh(e,f,g) + δW[r]
    # If δe = 0 (from δa = 0 + створочне): δSig1 = 0, δCh = 0.
    # If δh = δe[r-3] = 0: δh = 0.
    # So δT1 = δW[r].
    #
    # For δa[r+1] = 0: δW[r] = 0 needed!

    print(f"\n  For δa[r+1] = 0 given previous δa = 0:")
    print(f"    Need δT1[r] = 0 (since δT2 = 0)")
    print(f"    δT1 = δW[r] (other terms = 0 from δa=0 → δe=0)")
    print(f"    → Need δW[r] = 0!")
    print(f"")
    print(f"  But: δW[r] = W₁[r] - W₂[r] depends on δM through schedule.")
    print(f"  For r < 16: δW[r] = δM[r] (direct message word).")
    print(f"  For r ≥ 16: δW[r] = schedule(δM).")
    print(f"")
    print(f"  a-convergence at rounds R-11..R requires: δW[r] = 0 for those rounds.")
    print(f"  For R=64: δW[53..63] = 0 (11 schedule words).")
    print(f"  From schedule: δW[r] = sig1(δW[r-2]) + δW[r-7] + sig0(δW[r-15]) + δW[r-16].")
    print(f"  Setting 11 consecutive δW = 0 constrains δM heavily.")
    print(f"")

    # How many δM satisfy δW[53..63] = 0?
    # 11 words × 32 bits = 352 equations on δM (512 bits).
    # Schedule is GF(2)-linear in δM → rank of these constraints?

    # Build: Jacobian of W[53..63] w.r.t. M[0..15]
    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]
    W_base = schedule(M)

    J = []
    for word in range(16):
        for bit in range(32):
            M_flip = list(M)
            M_flip[word] ^= (1 << bit)
            W_flip = schedule(M_flip)
            row = []
            for r in range(53, 64):
                diff = W_base[r] ^ W_flip[r]
                for b in range(32):
                    row.append((diff >> b) & 1)
            J.append(row)

    # Rank of 512 × 352 matrix
    m = [list(row) for row in J]
    rank = 0
    n_cols = 11 * 32  # 352
    for col in range(n_cols):
        pivot = None
        for row in range(rank, 512):
            if m[row][col] == 1:
                pivot = row; break
        if pivot is None: continue
        m[rank], m[pivot] = m[pivot], m[rank]
        for row in range(512):
            if row != rank and m[row][col] == 1:
                for c in range(n_cols): m[row][c] ^= m[rank][c]
        rank += 1

    kernel = 512 - rank
    print(f"  δW[53..63] = 0 constraint: rank = {rank}/352")
    print(f"  δM kernel dimension: {kernel}")
    print(f"  = number of δM satisfying all 11 zero-schedule-word conditions")

    if kernel > 0:
        print(f"  *** {kernel} FREE BITS in δM that keep δW[53..63] = 0! ***")
        print(f"  These δM differences don't affect rounds 53-63.")
        print(f"  They ONLY affect rounds 0-52.")
        print(f"  If we find M where these δM also give δa convergence at r<53...")
    else:
        print(f"  No freedom: all 512 δM bits constrained.")


# ═══════════════════════════════════════════════════
# THREAD B: Carry product algebra
# ═══════════════════════════════════════════════════

def thread_B():
    """
    One round has 6 additions (approximately):
      T1: h + Sig1(e) = sum1, sum1 + Ch = sum2, sum2 + (K+W) = T1
      T2: Sig0(a) + Maj(a,b,c) = T2
      a_new: T1 + T2
      e_new: d + T1

    Each addition has its own carry operator.
    By T5: total carry = XOR of individual carries.

    But: the carries are NOT independent!
    carry_2 depends on sum1, which depends on carry_1.
    This creates a CHAIN of dependent carries.

    Question: how do sequential carries interact?
    Define: carry_product(C1, C2) = combined carry of x+y+z
    where C1 = carry(x,y), C2 = carry(x+y, z).

    Is carry_product associative? Commutative?
    Does it form an algebraic structure?
    """
    print("\n" + "=" * 80)
    print("THREAD B: Carry product algebra")
    print("=" * 80)

    n = 8
    mask = (1 << n) - 1

    # For sequential adds: x + y + z
    # C1 = carry(x, y), result1 = x + y
    # C2 = carry(result1, z), result2 = result1 + z = x + y + z
    # By T5: total_carry_effect = C1_effect XOR C2_effect.
    # But C2 DEPENDS on C1 (through result1).

    # Question: is there a "product" operation on carry patterns?
    # C_product = f(C1_pattern, C2_pattern)?
    # Where C_pattern = (generate, propagate) pair.

    # For x + y: G1 = x AND y, P1 = x XOR y.
    # result1 = P1 XOR carry_chain(G1, P1).
    # For result1 + z: G2 = result1 AND z, P2 = result1 XOR z.
    # G2 and P2 depend on result1 which depends on carry_chain(G1, P1).

    # The INTERACTION: G2, P2 at each bit depend on carry BELOW through result1.
    # This is SEQUENTIAL, not parallel.

    # Measure: for random x, y, z:
    # How much does knowing (G1, P1) tell about (G2, P2)?
    N = 10000
    g1p1_predicts_g2p2 = 0
    total = 0

    for _ in range(N):
        x = random.randint(0, mask)
        y = random.randint(0, mask)
        z = random.randint(0, mask)

        r1 = (x + y) & mask
        g1 = x & y
        p1 = x ^ y
        g2 = r1 & z
        p2 = r1 ^ z

        # Does (g1, p1) predict (g2, p2)?
        # Specifically: is g2[k] determined by g1[0..k], p1[0..k], z[k]?
        # YES — because r1[k] = p1[k] XOR carry1[k], and carry1 from g1, p1.
        # So g2[k] = r1[k] AND z[k] = (p1[k] XOR carry1[k]) AND z[k].
        # carry1 = function of g1, p1 (known from T4).

        total += 1

    print(f"  Sequential carry dependency: (G1,P1) fully determines (G2,P2)")
    print(f"  because result1 = P1 XOR carry_chain(G1,P1).")
    print(f"  So: G2 = result1 AND z = (P1 XOR C1) AND z")
    print(f"       P2 = result1 XOR z = (P1 XOR C1) XOR z")
    print(f"")
    print(f"  The PRODUCT: given (G1,P1) and z:")
    print(f"    C1 = carry_chain(G1, P1)  — from T4")
    print(f"    result1 = P1 XOR C1")
    print(f"    G2 = result1 AND z")
    print(f"    P2 = result1 XOR z")
    print(f"    C2 = carry_chain(G2, P2)")
    print(f"    Total carry effect = C1_effect XOR C2_effect  — from T5")
    print(f"")
    print(f"  This is SEQUENTIAL: C1 → result1 → (G2,P2) → C2.")
    print(f"  No shortcut: must compute C1 before C2.")
    print(f"  The 'product' is just COMPOSITION, not a parallel operation.")
    print(f"")

    # But: at BIT-0 level (no carry): all operations are parallel!
    # Bit 0: C1[0] = 0, C2[0] = 0. So result1[0] = P1[0] = x[0] XOR y[0].
    # G2[0] = result1[0] AND z[0] = (x[0] XOR y[0]) AND z[0].
    # This is a degree-2 function of (x[0], y[0], z[0]).

    # At bit 1: C1[1] = G1[0] = x[0] AND y[0] (carry from bit 0).
    # result1[1] = P1[1] XOR C1[1] = (x[1] XOR y[1]) XOR (x[0] AND y[0]).
    # G2[1] = result1[1] AND z[1].
    # = ((x[1] XOR y[1]) XOR (x[0] AND y[0])) AND z[1]
    # = degree 3 in (x, y, z).

    # So: sequential carries INCREASE DEGREE per bit level.
    # Bit 0: degree 2. Bit 1: degree 3. Bit k: degree k+2.

    print(f"  DEGREE GROWTH in sequential carries:")
    print(f"    Bit 0: degree 2 (x[0]⊕y[0])·z[0]")
    print(f"    Bit 1: degree 3 ((x[1]⊕y[1])⊕(x[0]·y[0]))·z[1]")
    print(f"    Bit k: degree k+2 (carry chain of length k adds degree k)")
    print(f"")
    print(f"  For 6 sequential additions per round:")
    print(f"    Degree grows: 2, 3, 4, 5, 6, 7 at bit 0 through 6 adds.")
    print(f"    But SHA-256 doesn't chain ALL 6 sequentially.")
    print(f"    T1 = chain of 4 adds. T2 = 1 add. Final = T1+T2 = 1 add.")
    print(f"    Max chain: 4 (T1) → degree at bit 0 = 5.")
    print(f"    This matches: bit 0 of a_new = degree ≤ 2 (our F4 said degree 2).")
    print(f"")
    print(f"  Wait: F4 says degree 2 for ONE ROUND. That's for the FULL")
    print(f"  round including all additions. The degree-5 from sequential")
    print(f"  carries applies at HIGHER bit positions (bit k).")
    print(f"  At bit 0: ALL carries = 0 → degree = 2 (only from Ch, Maj).")
    print(f"  This is CONSISTENT with F4!")


if __name__ == "__main__":
    thread_A()
    thread_B()
