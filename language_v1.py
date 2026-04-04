"""
LANGUAGE v1: Formal definitions of SHA-256 in our new language.

We don't say "almost independent" or "approximately affine".
We define EXACTLY what happens and give it names.

═══════════════════════════════════════════════════
DEFINITIONS
═══════════════════════════════════════════════════

DEFINITION 1: BIT-LAYER
  For a 32-bit SHA-256 computation, Layer(k) is the system of ALL
  equations involving bit position k of ALL registers at ALL rounds.
  Layer(k) has 8 × 64 = 512 equations (one per register per round).

DEFINITION 2: LAYER RANK
  R(k) = GF(2) rank of Layer(k)'s Jacobian.
  Measured: R(0) = 127, R(k) = 127 or 128 for all k.

DEFINITION 3: СТВОРОЧНОЕ ЧИСЛО (formal)
  The створочное число is the UNIQUE GF(2)-linear dependency between
  the a-register and e-register bit-0 trajectories.
  It reduces rank from 128 to 127 at certain bit positions.

DEFINITION 4: CARRY BRIDGE
  CarryBridge(k) = the carry bit from position k-1 to position k
  in any addition. It is a function of Layer(0)..Layer(k-1) values.
  CarryBridge(0) = 0 always (no carry into bit 0).
  CarryBridge(k) = MAJ(x[k-1], y[k-1], CarryBridge(k-1)).

DEFINITION 5: CONDITIONAL LAYER
  CondLayer(k | 0..k-1) = Layer(k) with all CarryBridges replaced by constants
  (since layers 0..k-1 determine them).
  CondLayer is AFFINE over GF(2) by construction.

DEFINITION 6: LAYER STACK
  SHA-256 = Stack(Layer(0), Layer(1), ..., Layer(31))
  where each Layer(k) is connected to Layer(k-1) by CarryBridge(k).

  The Stack has:
    Total rank = 512
    Layer ranks = [127, 127, 127, 127, 4, 0, 0, ...]
    (cumulative: 127, 254, 381, 508, 512, 512, ...)

═══════════════════════════════════════════════════
Now let's VERIFY these definitions hold exactly, not approximately.
═══════════════════════════════════════════════════
"""

import random
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, schedule
)


def gf2_rank(matrix, nrows, ncols):
    m = [list(row[:ncols]) for row in matrix[:nrows]]
    rank = 0
    for col in range(ncols):
        pivot = None
        for row in range(rank, len(m)):
            if m[row][col] == 1:
                pivot = row
                break
        if pivot is None:
            continue
        m[rank], m[pivot] = m[pivot], m[rank]
        for row in range(len(m)):
            if row != rank and m[row][col] == 1:
                for c in range(ncols):
                    m[row][c] ^= m[rank][c]
        rank += 1
    return rank


def verify_conditional_affinity():
    """
    EXACT TEST: Is CondLayer(1 | 0) truly affine?

    Method: Fix a message M. Compute bit-0 values of all registers at all rounds.
    These determine all CarryBridge(1) values.
    Now: evaluate bit-1 system with these fixed carries.
    Test DEGREE of the resulting system:
      Pick 3 random message perturbations (flipping only bit-1 positions).
      Check: f(a⊕b⊕c) = f(a) ⊕ f(b) ⊕ f(c) ⊕ f(0)?
      If yes for all triples → affine (degree 1).
    """
    print("=" * 80)
    print("EXACT VERIFICATION: Is CondLayer(1|0) truly affine?")
    print("=" * 80)

    # The idea: if we fix all bit-0 values, the carry into bit-1 of each
    # addition is FIXED. Then bit-1 of the result = x[1] XOR y[1] XOR carry(fixed).
    # This is AFFINE by construction (XOR with constant).

    # But SHA-256 has MULTIPLE additions chained: T1 = h + Sig1(e) + Ch + K + W.
    # The intermediate sums' bit-0 values affect carry into bit-1 of later additions.
    # Are ALL these carries also fixed when we fix the original state's bit-0 values?

    # YES: fixing bit-0 of state → fixing bit-0 of all intermediate values
    # (because bit-0 of sum = x[0] XOR y[0] XOR carry_in[0] = x[0] XOR y[0] XOR 0
    #  = pure XOR, ALSO determined by bit-0 of operands).

    # So: ALL bit-0 values everywhere in the computation are determined by
    # the bit-0 values of the initial state and message words.
    # And ALL carry-into-bit-1 values are determined by these bit-0 values.
    # Therefore: CondLayer(1|0) IS EXACTLY AFFINE.

    # Let's verify numerically: the degree-1 test.
    N = 500
    rng = random.Random(42)
    M_base = [rng.randint(0, MASK) for _ in range(16)]

    # We need to perturb ONLY bit 1 of message words (keeping bit 0 fixed).
    # Then bit-1 of output should be affine in bit-1 of input.

    violations = 0
    total = 0

    for trial in range(N):
        rng2 = random.Random(trial * 100)

        # Three random perturbations at bit position 1
        words_a = [rng2.randint(0, 15) for _ in range(3)]  # which word to flip
        # Flip bit 1 of those words
        M_a = list(M_base); M_a[words_a[0]] ^= 2  # bit 1 = value 2
        M_b = list(M_base); M_b[words_a[1]] ^= 2
        M_ab = list(M_base); M_ab[words_a[0]] ^= 2; M_ab[words_a[1]] ^= 2

        s0, _ = sha256_round_trace(M_base)
        sa, _ = sha256_round_trace(M_a)
        sb, _ = sha256_round_trace(M_b)
        sab, _ = sha256_round_trace(M_ab)

        # Check affinity at bit 1 of register a at round 64:
        # f(0) ⊕ f(a) ⊕ f(b) ⊕ f(a⊕b) should = 0 for affine f
        for reg in range(8):
            f0 = (s0[64][reg] >> 1) & 1
            fa = (sa[64][reg] >> 1) & 1
            fb = (sb[64][reg] >> 1) & 1
            fab = (sab[64][reg] >> 1) & 1

            check = f0 ^ fa ^ fb ^ fab
            total += 1
            if check != 0:
                violations += 1

    print(f"\n  Affinity test: flip ONLY bit 1 of message words")
    print(f"  f(0)⊕f(a)⊕f(b)⊕f(a⊕b) = 0 for affine f")
    print(f"  Violations: {violations}/{total}")

    if violations == 0:
        print(f"  ✓ CondLayer(1) IS EXACTLY AFFINE when only bit-1 inputs vary!")
    else:
        print(f"  ✗ NOT affine. Violations = {violations/total:.1%}")
        print(f"  But wait: we're flipping bit 1 of M, which ALSO changes bit-0")
        print(f"  through carry chains in the schedule... Need to be more careful.")

    # MORE CAREFUL: bit 1 of M[i] (= 2) doesn't change bit 0 of M[i].
    # But the schedule: W[16] = sig1(W[14]) + W[9] + sig0(W[1]) + W[0]
    # Changing bit 1 of W[0] → changes bit 1 of W[16], but also
    # sig0(W[1]) involves rotations that MIX bits.
    # sig0(x) = ROTR(x,7) XOR ROTR(x,18) XOR (x>>3).
    # Changing bit 1 of W[1] → changes bits 1-7=26, 1-18=15, and 1>>3=no effect.
    # So bit 1 of W[1] affects bits 15 and 26 of sig0(W[1])!
    # These are NOT bit 1 → our "only bit 1" assumption is VIOLATED by the schedule.

    print(f"\n  ISSUE: Flipping bit 1 of M[i] affects OTHER bit positions")
    print(f"  through the schedule's sig0/sig1 rotations.")
    print(f"  True layer isolation requires fixing ALL non-bit-1 effects.")


def verify_within_one_round():
    """
    Simpler test: within ONE ROUND (not 64), is the layer structure exact?

    For one round: a_new = T1 + T2.
    bit 1 of a_new = T1[1] XOR T2[1] XOR carry_from_bit_0.
    carry_from_bit_0 = T1[0] AND T2[0] = determined by bit-0 layer.

    If we fix the state (known) and vary only W[r] bit 1:
    Does a_new bit 1 change linearly with W[r] bit 1?
    """
    print("\n" + "=" * 80)
    print("SINGLE ROUND: Is bit-1 exactly affine in W[r] bit 1?")
    print("=" * 80)

    N = 1000
    affine_ok = 0
    total = 0

    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        states, W = sha256_round_trace(M, rounds=1)
        state = states[0]  # IV

        # Round 0: vary only bit 1 of W[0]
        for test_bit_val in [0, 1]:
            W_test = list(W)
            W_test[0] = (W[0] & ~2) | (test_bit_val << 1)  # Set bit 1

            a, b, c, d, e, f, g, h = state
            T1 = (h + Sig1(e) + Ch(e, f, g) + K[0] + W_test[0]) & MASK
            T2 = (Sig0(a) + Maj(a, b, c)) & MASK
            a_new = (T1 + T2) & MASK

            # a_new bit 1 should be affine in W[0] bit 1
            # (all other inputs fixed)

        # Full affinity test: flip bit 1 of two different inputs
        a, b, c, d, e, f, g, h = state

        def round_bit1(w0_val):
            T1 = (h + Sig1(e) + Ch(e, f, g) + K[0] + w0_val) & MASK
            T2 = (Sig0(a) + Maj(a, b, c)) & MASK
            return ((T1 + T2) >> 1) & 1

        w0 = W[0]
        f00 = round_bit1(w0)
        f10 = round_bit1(w0 ^ 2)  # flip bit 1

        # Also flip bit 1 of... we only have W[0] as free input at round 0.
        # Within one round, there's only one W input.
        # So affinity = f(w0 ^ 2) = f(w0) ^ constant.
        # This is just: does flipping bit 1 of W[0] always flip the same
        # bits of a_new?

        # Measure over many states: is the bit-1 effect of W[0] bit 1 always the same?
        # f(w0 ^ 2) XOR f(w0) should be constant across states... NO!
        # It depends on carry_from_bit_0, which depends on T1[0] AND T2[0],
        # which depends on W[0] bit 0 AND state.

        # So the effect of W[0] bit 1 on a_new bit 1 is:
        # Δ(a_new[1]) = Δ(T1[1]) XOR Δ(carry[1])
        # Δ(T1[1]) = Δ(W[0][1]) = 1 (direct)
        # Δ(carry[1]) = depends on whether T1[0] AND T2[0] changes
        #             = T1 with flipped W has different T1[0]?
        #             → T1[0] = h[0] XOR Sig1(e)[0] XOR Ch[0] XOR K[0] XOR W[0][0]
        #             → flipping W[0] bit 1 does NOT change W[0] bit 0
        #             → T1[0] does NOT change
        #             → carry[1] does NOT change
        #             → Δ(a_new[1]) = 1 always!

        # So within one round, flipping only bit 1 of W always flips bit 1 of a_new!
        # (Because bit 1 perturbation doesn't affect bit 0, hence carry is unchanged.)

        if f00 ^ f10 == 1:  # Always flips
            affine_ok += 1
        total += 1

    print(f"  Flipping bit 1 of W[0] always flips bit 1 of a_new: {affine_ok}/{total}")

    # Now: what about flipping bit 1 of the STATE? (not W)
    # Flipping a[1] → affects T2 through Sig0(a) and Maj(a,b,c).
    # Sig0(a)[1] = a[3] XOR a[14] XOR a[23] — does NOT include a[1]!
    # Maj(a,b,c)[1] = MAJ(a[1], b[1], c[1]) — DOES include a[1].
    # So flipping a[1] flips Maj[1] iff b[1] != c[1].
    # And carry into Sig0+Maj bit 1 depends on bit 0.

    # This is where it gets interesting: the effect of a[1] on a_new[1]
    # depends on b[1] and c[1] (through Maj). This is the DEGREE-2 part
    # from Ch and Maj, NOT from carry.

    # So: CondLayer(1|0) is affine W.R.T. carry, but STILL degree 2 from Ch/Maj!

    print(f"\n  IMPORTANT CORRECTION:")
    print(f"  Carry from bit 0 → bit 1: AFFINE (carry is constant given bit-0).")
    print(f"  BUT: Ch and Maj are DEGREE 2 at EVERY bit position.")
    print(f"  Ch(e,f,g)[k] = e[k]·f[k] ⊕ (1-e[k])·g[k] — degree 2 in e[k],f[k],g[k].")
    print(f"  This degree-2 nonlinearity is NOT from carry — it's from Boolean functions.")
    print(f"  It exists at EVERY bit position, including bit 0.")
    print(f"")
    print(f"  So: CondLayer(k|0..k-1) removes carry nonlinearity")
    print(f"  but PRESERVES Ch/Maj nonlinearity (degree 2).")
    print(f"  CondLayer is DEGREE 2, not degree 1 (not affine).")
    print(f"")
    print(f"  The Jacobian showed rank(layer k) = 127 because it's LINEAR approximation.")
    print(f"  The actual layer is QUADRATIC (degree 2 from Ch/Maj).")
    print(f"  But: quadratic with 127 independent constraints is still STRUCTURED.")
    print(f"  And: OMEGA already handles degree-2 systems (carry as variables).")
    print(f"")
    print(f"  REVISED PICTURE:")
    print(f"    SHA-256 = 4 layers, each QUADRATIC (degree 2), connected by carry.")
    print(f"    Carry between layers = the ONLY degree > 2 contribution.")
    print(f"    Within each layer: degree exactly 2 (from Ch/Maj).")
    print(f"    Total: degree 2 × 4 layers × carry connections = degree 2^4 = 16?")
    print(f"    No — layers are SEQUENTIAL, not parallel. Degree doesn't multiply.")
    print(f"    Effective: 4 coupled degree-2 systems.")


if __name__ == "__main__":
    verify_conditional_affinity()
    verify_within_one_round()
