"""
SHA-256 implementation with full bit-level carry tracking.

Every modular addition is decomposed into XOR + carry chain,
enabling the Q∩T decomposition:
  result[i] = a[i] ⊕ b[i] ⊕ carry[i-1]     (linear: Q part)
  carry[i]  = MAJ(a[i], b[i], carry[i-1])    (quadratic: T part)
"""

MASK32 = 0xFFFFFFFF

# SHA-256 round constants
K256 = [
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

# Initial hash values
IV256 = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
]


def rotr32(x, n):
    """Right rotate 32-bit word."""
    return ((x >> n) | (x << (32 - n))) & MASK32


def sigma0(x):
    """Small sigma_0 (schedule)."""
    return rotr32(x, 7) ^ rotr32(x, 18) ^ (x >> 3)


def sigma1(x):
    """Small sigma_1 (schedule)."""
    return rotr32(x, 17) ^ rotr32(x, 19) ^ (x >> 10)


def big_sigma0(x):
    """Capital Sigma_0 (compression)."""
    return rotr32(x, 2) ^ rotr32(x, 13) ^ rotr32(x, 22)


def big_sigma1(x):
    """Capital Sigma_1 (compression)."""
    return rotr32(x, 6) ^ rotr32(x, 11) ^ rotr32(x, 25)


def ch(e, f, g):
    """Choice function: Ch(e,f,g) = (e AND f) XOR (NOT e AND g)."""
    return ((e & f) ^ (~e & g)) & MASK32


def maj(a, b, c):
    """Majority function: Maj(a,b,c) = (a AND b) XOR (a AND c) XOR (b AND c)."""
    return (a & b) ^ (a & c) ^ (b & c)


def add32(a, b):
    """Standard mod 2^32 addition."""
    return (a + b) & MASK32


def add32_traced(a, b):
    """
    Add two 32-bit numbers with full carry tracking.

    Returns:
        result: (a + b) mod 2^32
        carry_chain: 32-bit int where bit i = carry INTO position i
                     (bit 0 is always 0: no carry into LSB)
        carry_out: carry out of bit 31 (overflow bit)
    """
    result = (a + b) & MASK32
    carry_chain = 0
    c = 0
    for i in range(32):
        # carry INTO position i
        carry_chain |= (c << i)
        ai = (a >> i) & 1
        bi = (b >> i) & 1
        # carry OUT of position i = MAJ(a[i], b[i], carry_in)
        c = (ai & bi) | (ai & c) | (bi & c)
    return result, carry_chain, c


def get_bit(word, bit):
    """Extract bit from 32-bit word."""
    return (word >> bit) & 1


def set_bit(word, bit, val):
    """Set a bit in a 32-bit word."""
    if val:
        return word | (1 << bit)
    return word & ~(1 << bit)


# ─── Bit-level mappings for GF(2)-linear operations ───

def sigma0_bit_sources(out_bit):
    """For sigma0: which input bits XOR to produce out_bit.
    sigma0(x) = ROTR(x,7) ^ ROTR(x,18) ^ SHR(x,3)
    Returns list of input bit indices (XOR them together).
    """
    sources = []
    sources.append((out_bit + 7) % 32)   # ROTR 7
    sources.append((out_bit + 18) % 32)  # ROTR 18
    if out_bit + 3 < 32:                  # SHR 3
        sources.append(out_bit + 3)
    return sources


def sigma1_bit_sources(out_bit):
    """For sigma1: which input bits XOR to produce out_bit.
    sigma1(x) = ROTR(x,17) ^ ROTR(x,19) ^ SHR(x,10)
    """
    sources = []
    sources.append((out_bit + 17) % 32)
    sources.append((out_bit + 19) % 32)
    if out_bit + 10 < 32:
        sources.append(out_bit + 10)
    return sources


def big_sigma0_bit_sources(out_bit):
    """For Sigma0: ROTR(x,2) ^ ROTR(x,13) ^ ROTR(x,22)."""
    return [
        (out_bit + 2) % 32,
        (out_bit + 13) % 32,
        (out_bit + 22) % 32,
    ]


def big_sigma1_bit_sources(out_bit):
    """For Sigma1: ROTR(x,6) ^ ROTR(x,11) ^ ROTR(x,25)."""
    return [
        (out_bit + 6) % 32,
        (out_bit + 11) % 32,
        (out_bit + 25) % 32,
    ]


# ─── Full SHA-256 compression with carry tracking ───

class SHA256Trace:
    """
    SHA-256 compression function result with full carry tracking.

    Attributes:
        hash_words: 8 output hash words
        schedule: W[0..63] expanded schedule
        states: state[0..R] where state[r] = (a,b,c,d,e,f,g,h)
        round_carries: per-round carry info for 7 additions
        sched_carries: schedule expansion carry info
        final_carries: final IV-addition carry info
    """

    def __init__(self):
        self.hash_words = []
        self.schedule = []
        self.states = []
        self.round_carries = []
        self.sched_carries = []
        self.final_carries = []
        self.num_rounds = 0


def sha256_compress_traced(msg_words, num_rounds=64, iv=None):
    """
    SHA-256 compression with full carry tracking.

    Args:
        msg_words: list of 16 uint32 message words
        num_rounds: number of compression rounds (1-64)
        iv: initial hash values (default: standard IV)

    Returns:
        SHA256Trace with all carries, states, and schedule
    """
    if iv is None:
        iv = list(IV256)
    assert len(msg_words) == 16

    trace = SHA256Trace()
    trace.num_rounds = num_rounds

    # ── Expand schedule ──
    w = list(msg_words)
    for i in range(16, max(num_rounds, 16)):
        s0 = sigma0(w[i - 15])
        s1 = sigma1(w[i - 2])
        t1, c1, co1 = add32_traced(w[i - 16], s0)
        t2, c2, co2 = add32_traced(t1, w[i - 7])
        t3, c3, co3 = add32_traced(t2, s1)
        w.append(t3)
        trace.sched_carries.append([(c1, co1), (c2, co2), (c3, co3)])
    trace.schedule = w

    # ── Compression rounds ──
    a, b, c, d, e, f, g, h = iv
    trace.states.append((a, b, c, d, e, f, g, h))

    for r in range(num_rounds):
        S1 = big_sigma1(e)
        ch_val = ch(e, f, g)
        S0 = big_sigma0(a)
        maj_val = maj(a, b, c)

        # T1 = h + Sigma1(e) + Ch(e,f,g) + K[r] + W[r]  (4 additions)
        t1_1, c1, co1 = add32_traced(h, S1)
        t1_2, c2, co2 = add32_traced(t1_1, ch_val)
        t1_3, c3, co3 = add32_traced(t1_2, K256[r])
        T1, c4, co4 = add32_traced(t1_3, w[r])

        # T2 = Sigma0(a) + Maj(a,b,c)  (1 addition)
        T2, c5, co5 = add32_traced(S0, maj_val)

        # e_new = d + T1  (1 addition)
        e_new, c6, co6 = add32_traced(d, T1)

        # a_new = T1 + T2  (1 addition)
        a_new, c7, co7 = add32_traced(T1, T2)

        trace.round_carries.append([
            (c1, co1),  # T1 add1: h + S1
            (c2, co2),  # T1 add2: + Ch
            (c3, co3),  # T1 add3: + K[r]
            (c4, co4),  # T1 add4: + W[r]
            (c5, co5),  # T2: S0 + Maj
            (c6, co6),  # e_new: d + T1
            (c7, co7),  # a_new: T1 + T2
        ])

        h, g, f, e = g, f, e, e_new
        d, c, b, a = c, b, a, a_new
        trace.states.append((a, b, c, d, e, f, g, h))

    # ── Final addition with IV ──
    final_state = [a, b, c, d, e, f, g, h]
    for i in range(8):
        result, fc, fco = add32_traced(final_state[i], iv[i])
        trace.hash_words.append(result)
        trace.final_carries.append((fc, fco))

    return trace


def get_all_carry_chains(trace):
    """
    Extract ALL carry chains as a flat list of 32-bit ints.
    Each entry is the carry-into chain for one addition.

    Order: schedule carries, round carries, final carries.
    """
    chains = []
    for sc in trace.sched_carries:
        for chain, _ in sc:
            chains.append(chain)
    for rc in trace.round_carries:
        for chain, _ in rc:
            chains.append(chain)
    for chain, _ in trace.final_carries:
        chains.append(chain)
    return chains


def get_carry_out_vector(trace):
    """Extract carry-out bits (one per addition) as a list of 0/1."""
    outs = []
    for sc in trace.sched_carries:
        for _, out in sc:
            outs.append(out)
    for rc in trace.round_carries:
        for _, out in rc:
            outs.append(out)
    for _, out in trace.final_carries:
        outs.append(out)
    return outs


def count_additions(trace):
    """Count total additions tracked."""
    n_sched = len(trace.sched_carries) * 3
    n_round = len(trace.round_carries) * 7
    n_final = len(trace.final_carries)
    return n_sched, n_round, n_final


# ─── Standard (non-traced) SHA-256 for verification ───

def sha256_compress(msg_words, num_rounds=64, iv=None):
    """Standard SHA-256 compression (no tracing, faster)."""
    if iv is None:
        iv = list(IV256)

    w = list(msg_words)
    for i in range(16, max(num_rounds, 16)):
        w.append(add32(add32(add32(sigma0(w[i-15]), w[i-16]),
                              w[i-7]), sigma1(w[i-2])))

    a, b, c, d, e, f, g, h = iv
    for r in range(num_rounds):
        T1 = add32(add32(add32(add32(h, big_sigma1(e)),
                                ch(e, f, g)), K256[r]), w[r])
        T2 = add32(big_sigma0(a), maj(a, b, c))
        h, g, f, e = g, f, e, add32(d, T1)
        d, c, b, a = c, b, a, add32(T1, T2)

    return [add32(x, y) for x, y in zip([a, b, c, d, e, f, g, h], iv)]
