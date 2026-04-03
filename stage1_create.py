"""
Stage 1, Part 5: CREATE candidate mathematical objects from scratch.

No analogies. No existing frameworks. Just define objects and test them.

Each candidate is:
1. A new OBJECT defined on SHA-256 states
2. A RULE for how it transforms across rounds
3. A TEST: does it carry non-trivial information through all 64 rounds?

"Non-trivial" = NOT random. If the object at round 64 is predictable
from round 0 better than random, it works.

Candidates:
1. PARITY CHAIN: track parity (XOR) of specific bit groups across rounds
2. CARRY SIGNATURE: the pattern of which additions generate carries
3. RESIDUE ORBIT: state mod small number (3, 5, 7) — non-binary arithmetic
4. COUPLING TENSOR: how a-bits and e-bits interact at each round
5. TRAJECTORY HASH: a rolling hash of the a-sequence itself
"""

import random
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, schedule
)


def candidate_1_parity_chain():
    """
    CANDIDATE 1: Parity Chain

    Define: P[r] = XOR of all 256 state bits at round r.
    This is a SINGLE BIT that summarizes the entire state.

    P[r+1] depends on P[r] through the round function.
    Since XOR is linear, and the only nonlinear parts are Ch, Maj, carries:
    P[r+1] = P[r] XOR (parity of new_a) XOR (parity of new_e) XOR (parity of lost d,h)

    Question: Is P[r] predictable? Does it carry information about M?
    """
    print("=" * 80)
    print("CANDIDATE 1: Parity Chain — XOR of all state bits")
    print("=" * 80)

    N = 1000
    # For each message, compute P[0], P[32], P[64]
    p0_vals = []
    p64_vals = []

    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        states, _ = sha256_round_trace(M)

        def parity(state):
            p = 0
            for word in state:
                p ^= bin(word).count('1') % 2
            return p

        p0_vals.append(parity(states[0]))
        p64_vals.append(parity(states[64]))

    # P[0] is always the same (IV is fixed)
    print(f"  P[0] = {p0_vals[0]} (constant, from IV)")

    # Is P[64] biased?
    p64_ones = sum(p64_vals)
    print(f"  P[64]: {p64_ones}/{N} ones = {p64_ones/N:.3f} (random = 0.500)")

    # Is P[64] correlated with any input bit?
    rng = random.Random(42)
    M_base = [rng.randint(0, MASK) for _ in range(16)]
    s_base, _ = sha256_round_trace(M_base)

    def parity(state):
        p = 0
        for word in state:
            p ^= bin(word).count('1') % 2
        return p

    p_base = parity(s_base[64])
    flips = 0
    for word in range(16):
        for bit in range(32):
            M_flip = list(M_base)
            M_flip[word] ^= (1 << bit)
            s_flip, _ = sha256_round_trace(M_flip)
            if parity(s_flip[64]) != p_base:
                flips += 1

    print(f"  P[64] flips when flipping input bit: {flips}/512 = {flips/512:.3f}")
    print(f"  (random = 0.500)")

    bias = abs(p64_ones/N - 0.5)
    sensitivity = abs(flips/512 - 0.5)
    print(f"  Bias: {bias:.4f}, Sensitivity deviation: {sensitivity:.4f}")
    print(f"  → {'INTERESTING' if bias > 0.05 or sensitivity > 0.05 else 'Random (no signal)'}")


def candidate_2_carry_signature():
    """
    CANDIDATE 2: Carry Signature

    For each round, record WHERE carries happen in the final addition (T1+T2).
    This gives a 32-bit "carry signature" per round.

    The carry signature is deterministic (given state and W).
    Question: Does the SEQUENCE of carry signatures over 64 rounds
    have structure that depends on M in a non-random way?

    Specifically: define CS = XOR of all 64 carry signatures.
    Is CS biased? Does it depend on M non-randomly?
    """
    print("\n" + "=" * 80)
    print("CANDIDATE 2: Carry Signature — rolling XOR of carry patterns")
    print("=" * 80)

    def get_carry_pattern(x, y):
        """Return 32-bit carry pattern for x + y."""
        carry = 0
        pattern = 0
        for i in range(32):
            xi = (x >> i) & 1
            yi = (y >> i) & 1
            if carry:
                pattern |= (1 << i)
            carry = (xi & yi) | (xi & carry) | (yi & carry)
        return pattern

    N = 500
    cs_values = []

    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        states, W = sha256_round_trace(M)

        cs = 0  # Rolling XOR of carry signatures
        for r in range(64):
            a, b, c, d, e, f, g, h = states[r]
            T1 = (h + Sig1(e) + Ch(e, f, g) + K[r] + W[r]) & MASK
            T2 = (Sig0(a) + Maj(a, b, c)) & MASK
            carry_pattern = get_carry_pattern(T1, T2)
            cs ^= carry_pattern

        cs_values.append(cs)

    unique = len(set(cs_values))
    hws = [bin(v).count('1') for v in cs_values]
    avg_hw = sum(hws) / N

    print(f"  Unique CS values: {unique}/{N}")
    print(f"  E[HW(CS)] = {avg_hw:.2f} (random = 16)")

    # Bit-level bias
    max_bias = 0
    for bit in range(32):
        ones = sum(1 for v in cs_values if (v >> bit) & 1)
        bias = abs(ones/N - 0.5)
        if bias > max_bias:
            max_bias = bias

    print(f"  Max bit bias: {max_bias:.4f} (significant if > 0.05)")
    print(f"  → {'INTERESTING' if max_bias > 0.05 else 'Random (no signal)'}")


def candidate_3_residue_orbit():
    """
    CANDIDATE 3: Residue Orbit

    SHA-256 works in Z/2^32. What if we project to Z/p for small prime p?
    This is a completely different algebraic world.

    Define: R_p[r] = a[r] mod p, for p = 3, 5, 7.

    In Z/p:
      - XOR becomes a strange operation (not + or ×)
      - Addition mod 2^32 becomes addition mod p (with wrap-around error)
      - Rotations become... something else

    Key: (a + b) mod p = ((a mod p) + (b mod p)) mod p  — EXACT.
    So addition LIFTS perfectly to Z/p.

    But: XOR does NOT lift to Z/p. (a XOR b) mod p ≠ f(a mod p, b mod p).

    So Z/p sees the additive structure PERFECTLY and ignores XOR.
    This is the OPPOSITE of GF(2) which sees XOR perfectly and ignores addition.

    Question: Does R_p carry information through 64 rounds?
    """
    print("\n" + "=" * 80)
    print("CANDIDATE 3: Residue Orbit — a[r] mod p for small primes")
    print("=" * 80)

    for p in [3, 5, 7, 13, 31]:
        N = 1000
        # For each message, compute (a[64] mod p)
        residues = []
        for seed in range(N):
            rng = random.Random(seed)
            M = [rng.randint(0, MASK) for _ in range(16)]
            states, _ = sha256_round_trace(M)
            residues.append(states[64][0] % p)

        # Count distribution
        counts = [0] * p
        for r in residues:
            counts[r] += 1

        expected = N / p
        max_dev = max(abs(c - expected) / expected for c in counts)
        chi2 = sum((c - expected)**2 / expected for c in counts)

        # Also: is residue at round 64 correlated with residue at round 0?
        # (round 0 = IV, so a[0] mod p is always the same)
        a0_mod_p = H0[0] % p

        print(f"  p={p:>2}: a[64] mod {p} distribution: {counts}")
        print(f"         max deviation: {max_dev:.3f}, χ²={chi2:.1f} (critical ~{2*p:.0f})")
        print(f"         a[0] mod {p} = {a0_mod_p} (constant)")
        print(f"         → {'BIASED' if chi2 > 3*p else 'Uniform (no signal)'}")


def candidate_4_coupling_measure():
    """
    CANDIDATE 4: Coupling Measure

    At each round, a_new and e_new share T1. Define:
      C[r] = a[r+1] XOR e[r+1] (their XOR difference)
      C[r] = (T1 + T2) XOR (d + T1)
            = T2 XOR d XOR carry_diff

    Since T2 depends only on a-branch and d = a[r-3], this is a
    function of the a-branch only. The e-branch appears only through carries.

    Question: Does C[r] mod 2 (just the LSB) have a pattern?
    C[r] mod 2 = T2[0] XOR d[0] (no carry at bit 0!)
               = (Sig0(a)[0] XOR Maj(a,b,c)[0]) XOR a[r-3][0]

    This is a PURE GF(2) function of a-branch bit 0 values!
    Track C[r] mod 2 across all 64 rounds.
    """
    print("\n" + "=" * 80)
    print("CANDIDATE 4: Coupling LSB — C[r] = (a[r+1] XOR e[r+1]) mod 2")
    print("=" * 80)

    N = 500
    # For each message, get the 64-bit string of C[r] mod 2
    all_chains = []
    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        states, _ = sha256_round_trace(M)

        chain = 0
        for r in range(64):
            c_bit = ((states[r+1][0] ^ states[r+1][4]) & 1)
            chain |= (c_bit << r)
        all_chains.append(chain)

    # How many unique chains?
    unique = len(set(all_chains))
    print(f"  Unique 64-bit chains: {unique}/{N}")

    # Average HW of chain (how many rounds have C[r] mod 2 = 1?)
    avg_hw = sum(bin(c).count('1') for c in all_chains) / N
    print(f"  E[HW(chain)] = {avg_hw:.2f} (random = 32)")

    # Per-round bias
    print(f"  Per-round P(C[r] mod 2 = 1):")
    biased_rounds = 0
    for r in range(64):
        ones = sum(1 for c in all_chains if (c >> r) & 1)
        p = ones / N
        if abs(p - 0.5) > 0.05:
            biased_rounds += 1
            if r < 10 or abs(p - 0.5) > 0.1:
                print(f"    r={r:>2}: P={p:.3f} {'*** BIASED ***' if abs(p-0.5) > 0.1 else '* mild *'}")

    print(f"  Total biased rounds (|P-0.5|>0.05): {biased_rounds}/64")

    # Autocorrelation of chain
    print(f"  Chain autocorrelation (lag 1-8):")
    for lag in range(1, 9):
        corr_count = 0
        total = 0
        for c in all_chains:
            for r in range(64 - lag):
                b1 = (c >> r) & 1
                b2 = (c >> (r + lag)) & 1
                if b1 == b2:
                    corr_count += 1
                total += 1
        corr = corr_count / total
        print(f"    lag {lag}: agreement = {corr:.4f} (random = 0.5000)")


def candidate_5_mod_trajectory():
    """
    CANDIDATE 5: Modular Trajectory Distance

    For two different messages M1, M2, define NOT their XOR difference,
    but their MODULAR distance:
      D[r] = |a[r](M1) - a[r](M2)| mod 2^32

    This is the additive distance, not XOR distance.
    XOR distance = Hamming weight of XOR. Additive distance = integer distance.

    These are VERY different. XOR treats each bit independently.
    Additive distance respects the carry structure.

    Question: Does additive distance evolve differently from XOR distance?
    Does it preserve information longer?
    """
    print("\n" + "=" * 80)
    print("CANDIDATE 5: Additive distance between trajectories")
    print("=" * 80)

    N = 500
    xor_dists = {r: [] for r in range(0, 65, 4)}
    add_dists = {r: [] for r in range(0, 65, 4)}

    for seed in range(N):
        rng = random.Random(seed)
        M1 = [rng.randint(0, MASK) for _ in range(16)]
        M2 = list(M1)
        M2[0] ^= 1  # 1-bit difference

        s1, _ = sha256_round_trace(M1)
        s2, _ = sha256_round_trace(M2)

        for r in range(0, 65, 4):
            xor_d = bin(s1[r][0] ^ s2[r][0]).count('1')
            add_d = (s1[r][0] - s2[r][0]) & MASK
            # Normalize additive distance to [0, 2^31]
            if add_d > (1 << 31):
                add_d = (1 << 32) - add_d
            # Convert to "bits" for comparison: log2(add_d + 1)
            import math
            add_bits = math.log2(add_d + 1) if add_d > 0 else 0

            xor_dists[r].append(xor_d)
            add_dists[r].append(add_bits)

    print(f"{'r':>3} | {'E[HW(XOR)]':>10} {'E[log2(ADD)]':>12} | {'σ(XOR)':>6} {'σ(ADD)':>6}")
    print("-" * 55)

    import math
    for r in sorted(xor_dists.keys()):
        avg_xor = sum(xor_dists[r]) / N
        avg_add = sum(add_dists[r]) / N
        std_xor = (sum((x - avg_xor)**2 for x in xor_dists[r]) / N) ** 0.5
        std_add = (sum((x - avg_add)**2 for x in add_dists[r]) / N) ** 0.5
        print(f"{r:>3} | {avg_xor:>10.2f} {avg_add:>12.2f} | {std_xor:>6.2f} {std_add:>6.2f}")

    print(f"\n  If ADD distance behaves differently from XOR → additive world has structure")
    print(f"  If both converge to same → no additional information in additive distance")


if __name__ == "__main__":
    candidate_1_parity_chain()
    candidate_2_carry_signature()
    candidate_3_residue_orbit()
    candidate_4_coupling_measure()
    candidate_5_mod_trajectory()
