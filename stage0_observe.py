"""
Stage 0: Full bitwise observation of SHA-256.

Goal: See what SHA-256 ACTUALLY does to information, round by round.
No metrics, no aggregates — raw bit-level data.

We trace two messages M and M' = M xor delta through all 64 rounds,
recording the COMPLETE 256-bit state at each step, and examining
where bits from the input "go" — not statistically, but deterministically.
"""

import random
import struct

# SHA-256 constants
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

H0 = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
]

MASK = 0xFFFFFFFF

def rotr(x, n):
    return ((x >> n) | (x << (32 - n))) & MASK

def sig0(x):
    return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)

def sig1(x):
    return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)

def Sig0(x):
    return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)

def Sig1(x):
    return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)

def Ch(e, f, g):
    return (e & f) ^ (~e & g) & MASK

def Maj(a, b, c):
    return (a & b) ^ (a & c) ^ (b & c)

def schedule(W16):
    W = list(W16)
    for i in range(16, 64):
        W.append((sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16]) & MASK)
    return W

def sha256_round_trace(msg16, rounds=64):
    """Run SHA-256 and return full state after each round."""
    W = schedule(msg16)
    a, b, c, d, e, f, g, h = H0

    states = []
    states.append((a, b, c, d, e, f, g, h))  # round 0 = IV

    for r in range(rounds):
        T1 = (h + Sig1(e) + Ch(e, f, g) + K[r] + W[r]) & MASK
        T2 = (Sig0(a) + Maj(a, b, c)) & MASK

        h = g
        g = f
        f = e
        e = (d + T1) & MASK
        d = c
        c = b
        b = a
        a = (T1 + T2) & MASK

        states.append((a, b, c, d, e, f, g, h))

    return states, W


def bits(word):
    """Return list of 32 bits (LSB first)."""
    return [(word >> i) & 1 for i in range(32)]

def state_bits(state):
    """Return 256 bits of state (a[0]..a[31], b[0]..b[31], ..., h[0]..h[31])."""
    result = []
    for word in state:
        result.extend(bits(word))
    return result

def xor_bits(b1, b2):
    return [a ^ b for a, b in zip(b1, b2)]

def hw(bitlist):
    return sum(bitlist)


def experiment_1_full_trace():
    """
    Experiment 1: Full bitwise trace of delta propagation.

    For M and M' = M xor (1 bit), trace the COMPLETE difference pattern
    through all 64 rounds. Look at raw data, not statistics.
    """
    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]

    # Flip bit 0 of W[0]
    M_prime = list(M)
    M_prime[0] ^= 1

    states_M, W_M = sha256_round_trace(M)
    states_Mp, W_Mp = sha256_round_trace(M_prime)

    print("=" * 80)
    print("EXPERIMENT 1: Full delta trace, M' = M xor bit0(W[0])")
    print("=" * 80)

    # For each round, show the FULL 256-bit XOR difference
    print(f"\n{'r':>3} | {'HW':>4} | a_diff       b_diff       c_diff       d_diff       e_diff       f_diff       g_diff       h_diff")
    print("-" * 130)

    for r in range(65):
        s1 = states_M[r]
        s2 = states_Mp[r]
        diffs = [s1[i] ^ s2[i] for i in range(8)]
        hws = [bin(d).count('1') for d in diffs]
        total_hw = sum(hws)

        diff_strs = [f"{d:08x}({hws[i]:>2})" for i, d in enumerate(diffs)]
        print(f"{r:>3} | {total_hw:>4} | {' '.join(diff_strs)}")

    return states_M, states_Mp, W_M, W_Mp


def experiment_2_bit_fate():
    """
    Experiment 2: Trace where ONE specific input bit ends up.

    For each of 32 bits in W[0], flip that bit and record which
    output bits change at round 64. This gives a 32x256 "fate matrix".
    """
    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]

    states_M, _ = sha256_round_trace(M)
    base_final = state_bits(states_M[64])

    print("\n" + "=" * 80)
    print("EXPERIMENT 2: Bit fate matrix — where does each input bit end up?")
    print("=" * 80)

    fate_matrix = []  # 32 rows (input bits) x 256 cols (output bits)

    for bit_pos in range(32):
        M_prime = list(M)
        M_prime[0] ^= (1 << bit_pos)
        states_Mp, _ = sha256_round_trace(M_prime)
        final_p = state_bits(states_Mp[64])

        diff = xor_bits(base_final, final_p)
        fate_matrix.append(diff)

        affected = sum(diff)
        # Which registers are affected?
        reg_hw = []
        for reg in range(8):
            reg_hw.append(sum(diff[reg*32:(reg+1)*32]))

        print(f"W[0] bit {bit_pos:>2} -> {affected:>3} output bits changed | "
              f"a:{reg_hw[0]:>2} b:{reg_hw[1]:>2} c:{reg_hw[2]:>2} d:{reg_hw[3]:>2} "
              f"e:{reg_hw[4]:>2} f:{reg_hw[5]:>2} g:{reg_hw[6]:>2} h:{reg_hw[7]:>2}")

    # Now: is there ANY structure in this matrix?
    # Check: are any two input bits producing IDENTICAL output patterns?
    print(f"\nChecking for identical fate rows...")
    for i in range(32):
        for j in range(i+1, 32):
            if fate_matrix[i] == fate_matrix[j]:
                print(f"  IDENTICAL: bit {i} and bit {j}")

    # Check: how many output bits are affected by ALL 32 input bits?
    all_affected = [1] * 256
    for bit_pos in range(32):
        for out_bit in range(256):
            if fate_matrix[bit_pos][out_bit] == 0:
                all_affected[out_bit] = 0
    print(f"Output bits affected by ALL 32 input bits: {sum(all_affected)}/256")

    # Check: how many output bits are affected by NO input bit?
    any_affected = [0] * 256
    for bit_pos in range(32):
        for out_bit in range(256):
            if fate_matrix[bit_pos][out_bit] == 1:
                any_affected[out_bit] = 1
    print(f"Output bits affected by at least one input bit: {sum(any_affected)}/256")

    # Check: XOR of all fate rows (linearity test)
    xor_all = [0] * 256
    for bit_pos in range(32):
        for out_bit in range(256):
            xor_all[out_bit] ^= fate_matrix[bit_pos][out_bit]
    print(f"XOR of all 32 fate rows: HW = {sum(xor_all)}/256")

    return fate_matrix


def experiment_3_carry_trace():
    """
    Experiment 3: Trace carry bits explicitly through each addition.

    For each round, record which bit positions generate carries
    in each of the additions (T1 components, T2 components).
    This is the "hidden" information that no existing approach captures fully.
    """
    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]
    W = schedule(M)

    a, b, c, d, e, f, g, h = H0

    print("\n" + "=" * 80)
    print("EXPERIMENT 3: Carry trace — where carries happen in each round")
    print("=" * 80)

    def add_with_carry(x, y):
        """Return (sum, carry_bits) where carry_bits[i] = carry INTO position i."""
        result = (x + y) & MASK
        carries = []
        c = 0
        for i in range(32):
            xi = (x >> i) & 1
            yi = (y >> i) & 1
            carries.append(c)
            c = (xi & yi) | (xi & c) | (yi & c)
        return result, carries

    total_carry_bits_per_round = []

    for r in range(64):
        # T1 = h + Sig1(e) + Ch(e,f,g) + K[r] + W[r]
        # Break down each addition
        s1 = Sig1(e)
        ch = Ch(e, f, g)

        sum1, c1 = add_with_carry(h, s1)
        sum2, c2 = add_with_carry(sum1, ch)
        sum3, c3 = add_with_carry(sum2, K[r])
        T1, c4 = add_with_carry(sum3, W[r])

        # T2 = Sig0(a) + Maj(a,b,c)
        s0 = Sig0(a)
        maj = Maj(a, b, c)
        T2, c5 = add_with_carry(s0, maj)

        # e_new = d + T1
        e_new, c6 = add_with_carry(d, T1)

        # a_new = T1 + T2
        a_new, c7 = add_with_carry(T1, T2)

        total_carries = sum(c1) + sum(c2) + sum(c3) + sum(c4) + sum(c5) + sum(c6) + sum(c7)
        total_carry_bits_per_round.append(total_carries)

        if r < 8 or r >= 60:
            print(f"r={r:>2}: carries per addition: "
                  f"h+S1={sum(c1):>2} +Ch={sum(c2):>2} +K={sum(c3):>2} +W={sum(c4):>2} | "
                  f"S0+Maj={sum(c5):>2} | d+T1={sum(c6):>2} T1+T2={sum(c7):>2} | "
                  f"total={total_carries}")
        elif r == 8:
            print("...")

        # Update state
        h = g; g = f; f = e; e = e_new
        d = c; c = b; b = a; a = a_new

    print(f"\nTotal carry bits over 64 rounds: {sum(total_carry_bits_per_round)}")
    print(f"Average per round: {sum(total_carry_bits_per_round)/64:.1f}")
    print(f"Max carries in a round: {max(total_carry_bits_per_round)} (round {total_carry_bits_per_round.index(max(total_carry_bits_per_round))})")
    print(f"Min carries in a round: {min(total_carry_bits_per_round)} (round {total_carry_bits_per_round.index(min(total_carry_bits_per_round))})")


def experiment_4_round_by_round_info():
    """
    Experiment 4: How much "new information" does each round create?

    Measure: for two messages differing in 1 bit, at each round,
    how many NEW bit positions become different that weren't different before?
    And how many positions that WERE different become same again?

    This shows the FLOW of the difference — not just its size (HW).
    """
    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]
    M_prime = list(M)
    M_prime[0] ^= 1

    states_M, _ = sha256_round_trace(M)
    states_Mp, _ = sha256_round_trace(M_prime)

    print("\n" + "=" * 80)
    print("EXPERIMENT 4: Information flow — new diffs, healed diffs per round")
    print("=" * 80)

    prev_diff = set()

    print(f"{'r':>3} | {'HW':>4} | {'new':>4} | {'healed':>4} | {'kept':>4} | {'turnover%':>9}")
    print("-" * 55)

    for r in range(65):
        s1_bits = state_bits(states_M[r])
        s2_bits = state_bits(states_Mp[r])

        curr_diff = set()
        for i in range(256):
            if s1_bits[i] != s2_bits[i]:
                curr_diff.add(i)

        new_diffs = curr_diff - prev_diff
        healed = prev_diff - curr_diff
        kept = curr_diff & prev_diff

        if len(curr_diff) > 0:
            turnover = (len(new_diffs) + len(healed)) / max(len(curr_diff), len(prev_diff)) * 100 if max(len(curr_diff), len(prev_diff)) > 0 else 0
        else:
            turnover = 0

        print(f"{r:>3} | {len(curr_diff):>4} | {len(new_diffs):>4} | {len(healed):>4} | {len(kept):>4} | {turnover:>8.1f}%")

        prev_diff = curr_diff


if __name__ == "__main__":
    experiment_1_full_trace()
    experiment_2_bit_fate()
    experiment_3_carry_trace()
    experiment_4_round_by_round_info()
