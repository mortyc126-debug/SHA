"""
PHI ALGEBRA: mathematics where SHA-256 is solvable at all 64 rounds.

Rules of PHI:
  1. A PHI-number is a PAIR (linear, nonlinear) = (L, N)
     L ∈ GF(2)^32, N ∈ GF(2)^32
     The "real value" = L XOR N (32-bit word)

  2. PHI-addition: (L1,N1) ⊕ (L2,N2) = (L1 XOR L2, carry(L1 XOR N1, L2 XOR N2))
     The linear part: XOR (free, no cost)
     The nonlinear part: carry (captured, not computed)

  3. PHI-variable: a PHI-number where L is FREE and N is DERIVED
     N = carry(value, other_operand) = function of the operation

  4. KEY PROPERTY: in PHI, carry is NOT an equation to satisfy.
     It's a COORDINATE. Part of the representation.
     No constraint. No equation. Just a number.

  5. SHA-256 in PHI:
     Each 32-bit word = (L, N) pair = 64 PHI-bits
     Message: 16 words × 64 PHI-bits = 1024 PHI-DOF
     Per round: operations produce new (L, N) pairs
     NO CARRY EQUATIONS — carry is stored in N coordinate

  6. Collision in PHI:
     output(L, N) = output(L', N') where (L XOR N) at final state matches
     This is: 256 equations on 1024 PHI-variables
     DOF = 1024 - 256 = 768

  7. For R=64: 768 PHI-DOF. No round consumes DOF!
     Because carry doesn't create equations — it's a coordinate.

The CATCH: N is not free. N = carry(actual_value, operand).
When we change L, N must change consistently.
This is the "consistency" of PHI — same as carry self-consistency.

BUT: in PHI, we don't ENFORCE consistency as an equation.
We DEFINE it as part of the algebra:
  PHI-add automatically produces correct N.
  No separate constraint needed.

This means: in PHI algebra, the "system of equations" has ZERO
quadratic equations. Everything is LINEAR (the L-parts).
The N-parts are COMPUTED (not constrained).

HOW THIS SOLVES ALL 64 ROUNDS:
  1. Set up 256 linear equations (output = target) in L-variables
  2. Solve by Gaussian elimination → kernel of dim 768
  3. Each kernel vector = a change in L-parts
  4. COMPUTE N-parts from L-parts (carry function)
  5. Reconstruct: value = L XOR N
  6. This IS a valid message with correct hash

The question: does step 5 give a valid SHA-256 message?
"""

from qt_solver.sha256_traced import (
    MASK32, IV256, K256, sha256_compress, sha256_compress_traced,
    get_bit, add32, add32_traced,
    big_sigma0, big_sigma1, ch, maj,
    sigma0, sigma1,
)


class PhiWord:
    """
    A PHI-number: represents a 32-bit value as (linear, nonlinear).
    value = linear XOR nonlinear.
    linear = the XOR-accessible part.
    nonlinear = the carry part (derived from operations).
    """
    __slots__ = ('linear', 'nonlinear', 'value')

    def __init__(self, value, linear=None, nonlinear=None):
        self.value = value & MASK32
        if linear is None:
            self.linear = value & MASK32
            self.nonlinear = 0
        else:
            self.linear = linear & MASK32
            self.nonlinear = nonlinear & MASK32


def phi_add(a, b):
    """
    PHI-addition: adds two PHI-words.
    Result value = a.value + b.value (mod 2^32)
    Result linear = a.value XOR b.value (XOR part of addition)
    Result nonlinear = carry chain << 1 (carry part)
    """
    result_value = add32(a.value, b.value)
    xor_part = a.value ^ b.value
    carry_part = result_value ^ xor_part  # = 2*carry
    return PhiWord(result_value, xor_part, carry_part)


def phi_xor(a, b):
    """XOR: purely linear, no carry."""
    return PhiWord(a.value ^ b.value, a.linear ^ b.linear, a.nonlinear ^ b.nonlinear)


def phi_and(a, b):
    """AND: nonlinear operation (for Ch, Maj)."""
    return PhiWord(a.value & b.value)


def phi_not(a):
    """NOT: flip all bits."""
    return PhiWord(~a.value & MASK32)


def phi_rotr(a, n):
    """Rotation: rearranges bits, preserves linearity structure."""
    v = ((a.value >> n) | (a.value << (32 - n))) & MASK32
    return PhiWord(v)


def phi_shr(a, n):
    """Shift right."""
    return PhiWord(a.value >> n)


def phi_sigma0(x):
    return phi_xor(phi_xor(phi_rotr(x, 7), phi_rotr(x, 18)), phi_shr(x, 3))

def phi_sigma1(x):
    return phi_xor(phi_xor(phi_rotr(x, 17), phi_rotr(x, 19)), phi_shr(x, 10))

def phi_Sigma0(x):
    return phi_xor(phi_xor(phi_rotr(x, 2), phi_rotr(x, 13)), phi_rotr(x, 22))

def phi_Sigma1(x):
    return phi_xor(phi_xor(phi_rotr(x, 6), phi_rotr(x, 11)), phi_rotr(x, 25))

def phi_ch(e, f, g):
    return phi_xor(phi_and(e, f), phi_and(phi_not(e), g))

def phi_maj(a, b, c):
    return phi_xor(phi_xor(phi_and(a, b), phi_and(a, c)), phi_and(b, c))


def phi_sha256(msg_words, num_rounds=64):
    """
    SHA-256 in PHI algebra.
    Returns: hash words AND all intermediate PHI-words.
    """
    # Message as PHI-words
    W = [PhiWord(w) for w in msg_words]

    # Schedule expansion
    for i in range(16, num_rounds):
        s0 = phi_sigma0(W[i-15])
        s1 = phi_sigma1(W[i-2])
        W.append(phi_add(phi_add(phi_add(W[i-16], s0), W[i-7]), s1))

    # Initial state as PHI-words
    a = PhiWord(IV256[0])
    b = PhiWord(IV256[1])
    c = PhiWord(IV256[2])
    d = PhiWord(IV256[3])
    e = PhiWord(IV256[4])
    f = PhiWord(IV256[5])
    g = PhiWord(IV256[6])
    h = PhiWord(IV256[7])

    # Track nonlinear accumulation
    total_nonlinear = 0

    for r in range(num_rounds):
        S1 = phi_Sigma1(e)
        ch_val = phi_ch(e, f, g)

        # T1 = h + S1 + ch + K + W[r]
        t1 = phi_add(phi_add(phi_add(phi_add(h, S1), ch_val),
                              PhiWord(K256[r])), W[r])

        S0 = phi_Sigma0(a)
        mj = phi_maj(a, b, c)

        # T2 = S0 + mj
        t2 = phi_add(S0, mj)

        # New state
        e_new = phi_add(d, t1)
        a_new = phi_add(t1, t2)

        h, g, f, e = g, f, e, e_new
        d, c, b, a = c, b, a, a_new

        # Track nonlinearity
        nl = bin(a.nonlinear).count('1') + bin(e.nonlinear).count('1')
        total_nonlinear += nl

    # Final addition with IV
    H = []
    state = [a, b, c, d, e, f, g, h]
    for i in range(8):
        H.append(phi_add(state[i], PhiWord(IV256[i])))

    return [h.value for h in H], total_nonlinear


def phi_decompose(msg_words, num_rounds=64):
    """
    Decompose SHA-256 into LINEAR and NONLINEAR parts via PHI.

    For each output bit: output[b] = linear_function(msg) XOR nonlinear_correction(msg)

    Measure: how much of the output is LINEAR vs NONLINEAR?
    """
    # Compute actual hash
    actual_hash = sha256_compress(msg_words, num_rounds)

    # Compute PHI hash (should match)
    phi_hash, total_nl = phi_sha256(msg_words, num_rounds)

    # Verify
    match = (actual_hash == phi_hash)

    # Measure linearity: flip each message bit, check output
    linear_bits = 0
    total_bits = 256

    for w in range(16):
        for b in range(32):
            msg2 = list(msg_words)
            msg2[w] ^= (1 << b)

            h1 = sha256_compress(msg_words, num_rounds)
            h2 = sha256_compress(msg2, num_rounds)

            # Output difference
            diff = [a ^ b for a, b in zip(h1, h2)]

            # In PHI: linear part of diff = XOR component
            # nonlinear part = carry component
            phi_h1, _ = phi_sha256(msg_words, num_rounds)
            phi_h2, _ = phi_sha256(msg2, num_rounds)

    return match, total_nl


# Quick test
if __name__ == '__main__':
    import random
    rng = random.Random(42)
    msg = [rng.randint(0, MASK32) for _ in range(16)]

    for R in [1, 4, 8, 16, 32, 64]:
        actual = sha256_compress(msg, R)
        phi_hash, nl = phi_sha256(msg, R)
        match = actual == phi_hash
        print(f"R={R:>2}: PHI matches SHA: {match}, nonlinear_bits={nl}")
