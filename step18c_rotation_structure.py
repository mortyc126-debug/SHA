#!/usr/bin/env python3
"""Step 18c: SHA-256-specific structure in rotation constants.

SHA-256 uses specific rotation constants:
  Sigma0: ROTR(2) XOR ROTR(13) XOR ROTR(22)
  Sigma1: ROTR(6) XOR ROTR(11) XOR ROTR(25)
  sigma0: ROTR(7) XOR ROTR(18) XOR SHR(3)
  sigma1: ROTR(17) XOR ROTR(19) XOR SHR(10)

These create specific GF(2)-linear maps. Their interaction with
modular addition creates the nonlinearity. But the LINEAR part
has structure determined by these constants.

Research directions:
1. What algebraic relationships exist between these specific constants?
2. Do the Sigma maps commute, or have specific composition structure?
3. What is the spectral structure (eigenvalues over GF(2))?
4. Do certain bit positions "resonate" across rounds?
5. What happens at the message schedule level — is there a short cycle?
"""

import numpy as np
import random
MASK32 = 0xFFFFFFFF

def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def hw(x): return bin(x).count('1')

# ============================================================
# Part 1: GF(2) matrix analysis of Sigma functions
# ============================================================
def gf2_matrix_analysis():
    """Build and analyze the GF(2) matrices for all four Sigma-like functions."""
    print("=" * 70)
    print("Part 1: GF(2) matrix spectral analysis")
    print("=" * 70)

    def build_gf2_matrix(func):
        M = np.zeros((32, 32), dtype=int)
        for j in range(32):
            val = func(1 << j)
            for i in range(32):
                M[i][j] = (val >> i) & 1
        return M

    Sigma0_func = lambda x: rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
    Sigma1_func = lambda x: rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
    sigma0_func = lambda x: rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
    sigma1_func = lambda x: rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)

    matrices = {
        'Sigma0': build_gf2_matrix(Sigma0_func),
        'Sigma1': build_gf2_matrix(Sigma1_func),
        'sigma0': build_gf2_matrix(sigma0_func),
        'sigma1': build_gf2_matrix(sigma1_func),
    }

    for name, M in matrices.items():
        rank = np.linalg.matrix_rank(M.astype(float))
        # Check if M^k = I for some k (order of M as GF(2) matrix)
        current = M.copy()
        order = None
        I = np.eye(32, dtype=int)
        for k in range(1, 100):
            current = np.mod(current @ M, 2)
            if np.array_equal(current, I):
                order = k + 1
                break

        # Compute M + I kernel (fixed points)
        MI = np.mod(M + I, 2)
        fp_rank = np.linalg.matrix_rank(MI.astype(float))
        fp_dim = 32 - fp_rank

        print(f"\n  {name}:")
        print(f"    Rank: {rank}")
        print(f"    Order (M^k = I): {order if order else '>100'}")
        print(f"    Fixed point space dim: {fp_dim}")

    # Check compositions
    print("\n  Composition analysis:")
    S0 = matrices['Sigma0']
    S1 = matrices['Sigma1']
    s0 = matrices['sigma0']
    s1 = matrices['sigma1']

    # Does Sigma0 * Sigma1 = Sigma1 * Sigma0?
    S0S1 = np.mod(S0 @ S1, 2)
    S1S0 = np.mod(S1 @ S0, 2)
    commute = np.array_equal(S0S1, S1S0)
    print(f"    Sigma0 * Sigma1 = Sigma1 * Sigma0? {commute}")

    s0s1 = np.mod(s0 @ s1, 2)
    s1s0 = np.mod(s1 @ s0, 2)
    commute2 = np.array_equal(s0s1, s1s0)
    print(f"    sigma0 * sigma1 = sigma1 * sigma0? {commute2}")

    # Cross-type
    S0s0 = np.mod(S0 @ s0, 2)
    s0S0 = np.mod(s0 @ S0, 2)
    print(f"    Sigma0 * sigma0 = sigma0 * Sigma0? {np.array_equal(S0s0, s0S0)}")

    # Sigma0 + Sigma1
    S0pS1 = np.mod(S0 + S1, 2)
    rank_sum = np.linalg.matrix_rank(S0pS1.astype(float))
    print(f"    Rank of (Sigma0 + Sigma1): {rank_sum}")

    # Sigma0 * Sigma1 + I
    S0S1pI = np.mod(S0S1 + np.eye(32, dtype=int), 2)
    rank_prod = np.linalg.matrix_rank(S0S1pI.astype(float))
    fp_prod = 32 - rank_prod
    print(f"    Fixed points of Sigma0*Sigma1: dim {fp_prod}")

    return matrices

# ============================================================
# Part 2: Message schedule as linear recurrence
# ============================================================
def schedule_recurrence():
    """Analyze the message schedule W[t] = σ1(W[t-2]) + W[t-7] + σ0(W[t-15]) + W[t-16]

    Over GF(2) (ignoring carries), this is a LINEAR recurrence.
    The characteristic polynomial determines its period and structure.
    """
    print("\n" + "=" * 70)
    print("Part 2: Message schedule as GF(2) linear recurrence")
    print("=" * 70)

    # Build the 512x512 state transition matrix over GF(2)
    # State = (W[t], W[t-1], ..., W[t-15]) = 16 x 32 = 512 bits
    # W[t+1] = sigma1(W[t-1]) XOR W[t-6] XOR sigma0(W[t-14]) XOR W[t-15]

    def build_sigma_matrix(func):
        M = np.zeros((32, 32), dtype=int)
        for j in range(32):
            val = func(1 << j)
            for i in range(32):
                M[i][j] = (val >> i) & 1
        return M

    sigma0_mat = build_sigma_matrix(lambda x: rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3))
    sigma1_mat = build_sigma_matrix(lambda x: rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10))

    # Build 512x512 state transition matrix
    # State vector: [W[t-15], W[t-14], ..., W[t-1], W[t]]
    # indexed as: block 0 = W[t-15], block 15 = W[t]
    # New W[t+1] = sigma1(W[t-1]) + W[t-6] + sigma0(W[t-14]) + W[t-15]
    # = sigma1(block 14) + block 9 + sigma0(block 1) + block 0

    n = 512  # 16 * 32
    T = np.zeros((n, n), dtype=int)

    # Shift: block i -> block i-1 (W[t-k] becomes W[t-k-1] in next state)
    # Actually: new state[0..14] = old state[1..15]
    for i in range(15):
        for b in range(32):
            T[i*32 + b][(i+1)*32 + b] = 1  # block i_new = block (i+1)_old

    # Block 15 (new W[t+1]):
    # = sigma1(old block 14) + old block 9 + sigma0(old block 1) + old block 0
    for b_out in range(32):
        # sigma1(block 14): old_block_14 = state[14*32 : 15*32]
        for b_in in range(32):
            if sigma1_mat[b_out][b_in]:
                T[15*32 + b_out][14*32 + b_in] ^= 1

        # block 9: old_block_9 = state[9*32 : 10*32]
        T[15*32 + b_out][9*32 + b_out] ^= 1

        # sigma0(block 1): old_block_1 = state[1*32 : 2*32]
        for b_in in range(32):
            if sigma0_mat[b_out][b_in]:
                T[15*32 + b_out][1*32 + b_in] ^= 1

        # block 0: old_block_0 = state[0 : 32]
        T[15*32 + b_out][0*32 + b_out] ^= 1

    rank_T = np.linalg.matrix_rank(T.astype(float))
    print(f"\n  State transition matrix T (512x512):")
    print(f"    Rank: {rank_T}")

    # Check period: T^k = I?
    # For 512x512 GF(2) matrix, period could be up to 2^512 - 1
    # Let's check small periods
    current = T.copy()
    I512 = np.eye(n, dtype=int)
    print(f"\n  Checking small periods of T (T^k = I):")
    for k in [32, 64, 128, 256, 512]:
        # Compute T^k by repeated squaring
        result = np.eye(n, dtype=int)
        base = T.copy()
        exp = k
        while exp > 0:
            if exp % 2 == 1:
                result = np.mod(result @ base, 2)
            base = np.mod(base @ base, 2)
            exp //= 2
        is_identity = np.array_equal(result, I512)
        print(f"    T^{k} = I? {is_identity}")

    # Fixed points of T: kernel of T + I
    TpI = np.mod(T + I512, 2)
    fp_rank = np.linalg.matrix_rank(TpI.astype(float))
    fp_dim = n - fp_rank
    print(f"\n  Fixed points of T: dim {fp_dim}")
    print(f"  => {2**fp_dim if fp_dim <= 20 else '2^'+str(fp_dim)} fixed message schedules")

    # T^2 fixed points (period-2 cycles)
    T2 = np.mod(T @ T, 2)
    T2pI = np.mod(T2 + I512, 2)
    fp2_rank = np.linalg.matrix_rank(T2pI.astype(float))
    fp2_dim = n - fp2_rank
    print(f"  Period-2 cycles: dim {fp2_dim} (includes fixed points)")

    return T

# ============================================================
# Part 3: Interaction between Sigma and message schedule
# ============================================================
def sigma_schedule_interaction():
    """The round function uses Sigma0(a) and Sigma1(e),
    while the schedule uses sigma0 and sigma1.

    Do the rotation constants have algebraic relationships?
    """
    print("\n" + "=" * 70)
    print("Part 3: Rotation constant relationships")
    print("=" * 70)

    # SHA-256 rotation constants:
    # Sigma0: 2, 13, 22
    # Sigma1: 6, 11, 25
    # sigma0: 7, 18, 3 (SHR)
    # sigma1: 17, 19, 10 (SHR)

    constants = {
        'Sigma0': (2, 13, 22),
        'Sigma1': (6, 11, 25),
        'sigma0': (7, 18, 3),
        'sigma1': (17, 19, 10),
    }

    print("\n  Rotation constants:")
    for name, (r1, r2, r3) in constants.items():
        print(f"    {name}: ({r1}, {r2}, {r3})")
        print(f"      Sum: {r1+r2+r3}")
        print(f"      Pairwise sums: {r1+r2}, {r1+r3}, {r2+r3}")
        print(f"      Differences: {r2-r1}, {r3-r2}, {r3-r1}")
        print(f"      Mod 32: ({r1%32}, {r2%32}, {r3%32})")

    print("\n  Cross-function relationships:")
    # Sigma0 and Sigma1
    s0 = constants['Sigma0']
    s1 = constants['Sigma1']
    print(f"    Sigma0 + Sigma1 pairwise: ", end="")
    pairs = []
    for r0 in s0:
        for r1 in s1:
            s = (r0 + r1) % 32
            pairs.append(s)
            print(f"{r0}+{r1}={s} ", end="")
    print()

    # Check if any sum/difference creates a cycle
    print(f"\n    Differences Sigma0[i] - Sigma1[j] mod 32:")
    for r0 in s0:
        for r1 in s1:
            d = (r0 - r1) % 32
            print(f"      {r0} - {r1} = {d} mod 32")

    # Sigma0 rotation 2: what happens if we compose with itself?
    # ROTR(2) applied 16 times = ROTR(32) = identity
    print(f"\n  Periodicity of individual rotations:")
    for r in [2, 6, 7, 11, 13, 17, 18, 19, 22, 25]:
        period = 32 // __import__('math').gcd(r, 32)
        print(f"    ROTR({r:2d}): period = {period}")

    # The GCD structure
    print(f"\n  GCD analysis:")
    all_rotations = [2, 6, 7, 11, 13, 17, 18, 19, 22, 25]
    from math import gcd
    g = all_rotations[0]
    for r in all_rotations[1:]:
        g = gcd(g, r)
    print(f"    GCD of all rotation constants: {g}")
    print(f"    GCD(Sigma0): {gcd(gcd(2,13),22)}")
    print(f"    GCD(Sigma1): {gcd(gcd(6,11),25)}")
    print(f"    GCD(sigma0): {gcd(gcd(7,18),3)}")
    print(f"    GCD(sigma1): {gcd(gcd(17,19),10)}")

# ============================================================
# Part 4: Looking for resonance — do specific bit patterns recur?
# ============================================================
def resonance_search():
    """Track specific bit positions through multiple rounds.

    The rotation constants create specific "pathways" for bits.
    A bit at position p in register a after Sigma0 ends up at
    positions (p-2)%32, (p-13)%32, (p-22)%32 (XOR'd).

    After many rounds, do certain bit positions "resonate" —
    i.e., return to their starting position more often than random?
    """
    print("\n" + "=" * 70)
    print("Part 4: Bit position resonance through rounds")
    print("=" * 70)

    # Track where a single bit goes through the linear part of SHA-256
    # Focus on the (a -> Sigma0(a) -> Maj -> T2 -> a_new) path
    # and (e -> Sigma1(e) -> Ch -> T1 -> e_new) path

    # For Sigma0(a), bit p → bits at (p+2)%32, (p+13)%32, (p+22)%32
    # (rotations are right-rotations, so ROTR(2) sends bit p to (p+2)%32)

    # Track bit position evolution
    def sigma0_bits(p):
        """Which bits does input bit p of Sigma0 contribute to?"""
        return {(p+2)%32, (p+13)%32, (p+22)%32}

    def sigma1_bits(p):
        return {(p+6)%32, (p+11)%32, (p+25)%32}

    # For each starting bit, track how many rounds until it "returns"
    print("\n  Sigma0 bit pathway (linear XOR model):")
    for start_bit in [0, 1, 7, 15, 31]:
        bits = {start_bit}
        for r in range(10):
            new_bits = set()
            for b in bits:
                new_bits ^= sigma0_bits(b)  # XOR set (symmetric difference)
            bits = new_bits
            if start_bit in bits:
                print(f"    bit {start_bit}: returns after {r+1} rounds (|support|={len(bits)})")
                break
        else:
            print(f"    bit {start_bit}: no return in 10 rounds (|support|={len(bits)})")

    print("\n  Sigma1 bit pathway:")
    for start_bit in [0, 1, 7, 15, 31]:
        bits = {start_bit}
        for r in range(10):
            new_bits = set()
            for b in bits:
                new_bits ^= sigma1_bits(b)
            bits = new_bits
            if start_bit in bits:
                print(f"    bit {start_bit}: returns after {r+1} rounds (|support|={len(bits)})")
                break
        else:
            print(f"    bit {start_bit}: no return in 10 rounds (|support|={len(bits)})")

    # How quickly do bits spread?
    print("\n  Bit spreading rate:")
    for name, func in [("Sigma0", sigma0_bits), ("Sigma1", sigma1_bits)]:
        bits = {0}
        print(f"    {name}: ", end="")
        for r in range(8):
            new_bits = set()
            for b in bits:
                new_bits ^= func(b)
            bits = new_bits
            print(f"|round{r+1}|={len(bits)} ", end="")
        print()

# ============================================================
# Part 5: Does the specific K schedule create structure?
# ============================================================
def k_schedule_analysis():
    """The round constants K[0..63] are derived from cube roots of primes.
    Do they have algebraic properties that create structure?"""
    print("\n" + "=" * 70)
    print("Part 5: Round constant K analysis")
    print("=" * 70)

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

    # XOR consecutive K values
    print("\n  K[t] XOR K[t+1] HW:")
    xor_hws = []
    for t in range(63):
        xh = hw(K[t] ^ K[t+1])
        xor_hws.append(xh)
    avg_xor_hw = sum(xor_hws) / len(xor_hws)
    print(f"    Average: {avg_xor_hw:.1f} (random expectation: 16.0)")

    # K[t] + K[t+1] mod 2^32 patterns
    print(f"\n  K[t] + K[t+1] mod 2^32 HW:")
    sum_hws = []
    for t in range(63):
        sh = hw((K[t] + K[t+1]) & MASK32)
        sum_hws.append(sh)
    print(f"    Average: {sum(sum_hws)/len(sum_hws):.1f}")

    # Are there K values where K[t] = Sigma0(K[t']) or K[t] = Sigma1(K[t'])?
    Sigma0_func = lambda x: rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
    Sigma1_func = lambda x: rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)

    print(f"\n  K[t] that are Sigma images of other K values:")
    for t1 in range(64):
        for t2 in range(64):
            if t1 != t2:
                if Sigma0_func(K[t2]) == K[t1]:
                    print(f"    K[{t1}] = Sigma0(K[{t2}])")
                if Sigma1_func(K[t2]) == K[t1]:
                    print(f"    K[{t1}] = Sigma1(K[{t2}])")

    # Linear dependencies mod 2 (XOR)
    print(f"\n  XOR linear dependencies among K[0..63]:")
    deps_found = 0
    for t1 in range(64):
        for t2 in range(t1+1, 64):
            for t3 in range(t2+1, 64):
                if K[t1] ^ K[t2] ^ K[t3] == 0:
                    print(f"    K[{t1}] XOR K[{t2}] XOR K[{t3}] = 0")
                    deps_found += 1
    if deps_found == 0:
        print(f"    None found (K values are XOR-independent in triples)")


if __name__ == "__main__":
    matrices = gf2_matrix_analysis()
    T = schedule_recurrence()
    sigma_schedule_interaction()
    resonance_search()
    k_schedule_analysis()

    print("\n" + "=" * 70)
    print("Step 18c Complete")
    print("=" * 70)
