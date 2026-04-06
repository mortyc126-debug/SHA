"""
SHA-256 compression function with full carry tracing.

Tracks every intermediate value and every carry bit at every addition.
Used for building and verifying the OMEGA system.
"""

MASK32 = 0xFFFFFFFF

# SHA-256 round constants
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
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
]

# SHA-256 initial hash values
IV = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
]


def get_bit(x, k):
    """Get bit k of integer x."""
    return (x >> k) & 1


def rotr(x, n):
    """32-bit right rotation."""
    return ((x >> n) | (x << (32 - n))) & MASK32


def sigma0(x):
    """Σ₀(x) = ROTR(x,2) ⊕ ROTR(x,13) ⊕ ROTR(x,22)"""
    return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)


def sigma1(x):
    """Σ₁(x) = ROTR(x,6) ⊕ ROTR(x,11) ⊕ ROTR(x,25)"""
    return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)


def ssigma0(x):
    """σ₀(x) = ROTR(x,7) ⊕ ROTR(x,18) ⊕ (x >> 3) — schedule"""
    return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)


def ssigma1(x):
    """σ₁(x) = ROTR(x,17) ⊕ ROTR(x,19) ⊕ (x >> 10) — schedule"""
    return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)


def ch(e, f, g):
    """Ch(e,f,g) = (e ∧ f) ⊕ (¬e ∧ g)"""
    return (e & f) ^ (~e & g) & MASK32


def maj(a, b, c):
    """Maj(a,b,c) = (a ∧ b) ⊕ (a ∧ c) ⊕ (b ∧ c)"""
    return (a & b) ^ (a & c) ^ (b & c)


def add32_traced(x, y):
    """
    32-bit addition with carry tracing.
    Returns (result, carry_bits) where carry_bits is a 32-bit int,
    carry_bits[k] = carry INTO position k.
    carry_bits[0] = 0 always.
    """
    result = (x + y) & MASK32
    carry = 0
    c = 0
    for k in range(32):
        if k > 0:
            carry |= (c << k)
        s = get_bit(x, k) + get_bit(y, k) + c
        c = s >> 1  # carry out
    return result, carry


def sha256_compress(msg, rounds=64, iv=None):
    """
    SHA-256 compression. msg = list of 16 uint32.
    Returns list of 8 uint32 (hash).
    """
    if iv is None:
        iv = IV
    w = list(msg[:16])
    for i in range(16, max(rounds, 16)):
        w.append((ssigma1(w[i-2]) + w[i-7] + ssigma0(w[i-15]) + w[i-16]) & MASK32)

    a, b, c, d, e, f, g, h = iv

    for r in range(rounds):
        t1 = (h + sigma1(e) + ch(e, f, g) + K[r] + w[r]) & MASK32
        t2 = (sigma0(a) + maj(a, b, c)) & MASK32
        h = g
        g = f
        f = e
        e = (d + t1) & MASK32
        d = c
        c = b
        b = a
        a = (t1 + t2) & MASK32

    return [(a + iv[0]) & MASK32, (b + iv[1]) & MASK32,
            (c + iv[2]) & MASK32, (d + iv[3]) & MASK32,
            (e + iv[4]) & MASK32, (f + iv[5]) & MASK32,
            (g + iv[6]) & MASK32, (h + iv[7]) & MASK32]


def sha256_compress_traced(msg, rounds=64, iv=None):
    """
    SHA-256 compression with FULL trace of every intermediate.

    Returns dict with:
      'hash': list of 8 uint32
      'W': list of schedule words (len >= rounds)
      'states': list of (a,b,c,d,e,f,g,h) at each round (len = rounds+1, [0] = IV)
      'T1': list of T1 values per round
      'T2': list of T2 values per round
      'carries': dict mapping (round, addition_index) -> carry_bits (32-bit int)
          addition_index:
            0: h + Σ₁(e)
            1: prev + Ch(e,f,g)
            2: prev + K[r]
            3: prev + W[r]          -> T1 complete
            4: Σ₀(a) + Maj(a,b,c)  -> T2 complete
            5: d + T1               -> e_new
            6: T1 + T2              -> a_new
    """
    if iv is None:
        iv = list(IV)
    else:
        iv = list(iv)

    # message schedule
    w = list(msg[:16])
    schedule_carries = {}
    for i in range(16, max(rounds, 16)):
        # W[i] = σ₁(W[i-2]) + W[i-7] + σ₀(W[i-15]) + W[i-16]
        # 3 additions
        v1 = ssigma1(w[i-2])
        v2 = w[i-7]
        v3 = ssigma0(w[i-15])
        v4 = w[i-16]
        s1, c1 = add32_traced(v1, v2)
        s2, c2 = add32_traced(s1, v3)
        s3, c3 = add32_traced(s2, v4)
        w.append(s3)
        schedule_carries[i] = (c1, c2, c3)

    a, b, c, d, e, f, g, h = iv
    states = [(a, b, c, d, e, f, g, h)]
    T1_list = []
    T2_list = []
    carries = {}

    for r in range(rounds):
        # T1 = h + Σ₁(e) + Ch(e,f,g) + K[r] + W[r]
        s1, c0 = add32_traced(h, sigma1(e))
        s2, c1 = add32_traced(s1, ch(e, f, g))
        s3, c2 = add32_traced(s2, K[r])
        t1, c3 = add32_traced(s3, w[r])

        # T2 = Σ₀(a) + Maj(a,b,c)
        t2, c4 = add32_traced(sigma0(a), maj(a, b, c))

        # e_new = d + T1
        e_new, c5 = add32_traced(d, t1)

        # a_new = T1 + T2
        a_new, c6 = add32_traced(t1, t2)

        carries[(r, 0)] = c0
        carries[(r, 1)] = c1
        carries[(r, 2)] = c2
        carries[(r, 3)] = c3
        carries[(r, 4)] = c4
        carries[(r, 5)] = c5
        carries[(r, 6)] = c6

        T1_list.append(t1)
        T2_list.append(t2)

        h = g
        g = f
        f = e
        e = e_new
        d = c
        c = b
        b = a
        a = a_new

        states.append((a, b, c, d, e, f, g, h))

    # final addition
    final_carries = {}
    hash_out = []
    final_state = [a, b, c, d, e, f, g, h]
    for i in range(8):
        val, fc = add32_traced(final_state[i], iv[i])
        hash_out.append(val)
        final_carries[i] = fc

    return {
        'hash': hash_out,
        'W': w,
        'states': states,
        'T1': T1_list,
        'T2': T2_list,
        'carries': carries,
        'schedule_carries': schedule_carries,
        'final_carries': final_carries,
    }
