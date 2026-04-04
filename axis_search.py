"""
AXIS SEARCH: What to build mathematics AROUND?

We need ONE central object. Testing candidates:

AXIS A: The CARRY TRAJECTORY as primary object
  Not "state at round r" but "all carries across all rounds"
  = 64 rounds × ~7 additions × 32 bits ≈ 14K bits
  These 14K bits are FULLY determined by M (512 bits).
  The mapping M → carry_trajectory is a 512→14K function.
  What kind of function is it? How compressible?

AXIS B: The ROUND FUNCTION as operator on function spaces
  Not F(x) for specific x, but F as an operator.
  What are the EIGENVALUES / EIGENSPACES of F (mod 2)?
  An eigenspace of F = a subspace preserved by the round function.
  If such a subspace exists → invariant through all rounds.

AXIS C: MODULAR DECOMPOSITION at word level
  SHA-256 uses 32-bit words. What if we decompose each word
  not into bits (GF(2)) or keep it whole (Z/2^32),
  but split into NIBBLES (4-bit chunks)?
  4 bits = one hex digit. 8 nibbles per word.
  Within a nibble, carry chain ≤ 3 bits (degree ≤ 4).
  Between nibbles: only the carry-out bit connects them.
  This is a NATURAL granularity (carry segments average ~4 bits).

AXIS D: The W-SCHEDULE as the only external input
  The recurrence a[r+1] = F(a[r..r-7], K[r], W[r])
  has K[r] as known constants and W[r] as the only message-dependent part.
  For r < 16: W[r] = M[r] (direct message word).
  For r ≥ 16: W[r] = schedule(M) — deterministic from M.

  What if we think of SHA-256 as a FILTER driven by W[r]?
  The state is the filter's internal state.
  The output H is the final state.
  Finding a preimage = finding the INPUT SIGNAL W[0..63] that produces H.
  But W is constrained: W[16..63] = f(W[0..15]).

AXIS E: The JOINT a+e value (створочное число)
  We proved: e[r] = a[r] + a[r-4] - T2[r-1]
  So a[r] - e[r] = T2[r-1] - a[r-4]
  Define: ψ[r] = a[r] - e[r] mod 2^32
  ψ[r] = Sig0(a[r-1]) + Maj(a[r-1], a[r-2], a[r-3]) - a[r-4]
  This depends only on a[r-1..r-4] — NO message, NO constants.
  ψ is a PURE structural function of the a-trajectory.
"""

import random
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, schedule
)


def axis_C_nibble():
    """
    AXIS C: Nibble decomposition.

    Split each 32-bit word into 8 nibbles (4 bits each).
    Within a nibble, carry is at most 3 levels deep.
    Between nibbles: only 1 carry bit connects them.

    SHA-256 viewed as operations on nibble-vectors:
    - XOR: parallel per nibble (no cross-nibble coupling)
    - AND: parallel per nibble
    - Rotation: SHUFFLES nibbles (but also shifts within nibbles!)
    - Addition: carry CHAINS through nibbles sequentially

    ROTR(x, 2): shifts all bits by 2. Nibble boundaries at 0,4,8,12,...
    A 2-bit rotation takes bits {0,1,2,3} → {2,3,4,5}. This CROSSES nibble boundary!

    So rotation doesn't respect nibble boundaries. BUT: it crosses them
    in a KNOWN, FIXED way. We can track exactly how.
    """
    print("=" * 80)
    print("AXIS C: Nibble decomposition — carry coupling between nibbles")
    print("=" * 80)

    # For each addition in SHA-256, the carry between nibbles is just 1 bit.
    # So the inter-nibble coupling is MINIMAL: one bit per boundary.
    # 7 boundaries per word × ~6 additions per round = ~42 inter-nibble carry bits.

    # The INTRA-nibble structure is degree ≤ 4 (carry chain ≤ 3).

    # Question: If we know the inter-nibble carries, how does the problem decompose?
    # With inter-nibble carries fixed: each nibble is an INDEPENDENT 4-bit system.
    # 8 nibbles × 4 bits = 32 bits per word, but only 4 bits per nibble → 2^4 = 16 values.

    # For R rounds: 42 × R inter-nibble carry bits to "guess".
    # At R=64: 42 × 64 = 2688 bits.
    # But message = 512 bits → only 512 bits of freedom.
    # So most inter-nibble carries are DETERMINED.

    # What's the EFFECTIVE freedom in inter-nibble carries?
    N = 800
    R = 16

    all_inter = []
    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        states, W = sha256_round_trace(M, rounds=R)

        inter_carries = []
        for r in range(R):
            a, b, c, d, e, f, g, h = states[r]
            T1 = (h + Sig1(e) + Ch(e, f, g) + K[r] + W[r]) & MASK
            T2 = (Sig0(a) + Maj(a, b, c)) & MASK

            # Inter-nibble carries for T1+T2: carry out of positions 3, 7, 11, 15, 19, 23, 27
            c_bit = 0
            for k in range(32):
                xk = (T1 >> k) & 1
                yk = (T2 >> k) & 1
                c_bit = (xk & yk) | (xk & c_bit) | (yk & c_bit)
                if k in [3, 7, 11, 15, 19, 23, 27]:
                    inter_carries.append(c_bit)

        all_inter.append(inter_carries)

    # GF(2) rank of inter-nibble carry matrix
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

    ncols = 7 * R  # 7 inter-nibble boundaries per round
    rank = gf2_rank(all_inter, N, ncols)
    print(f"  R={R}: {ncols} inter-nibble carry bits, GF(2) rank = {rank}/{ncols}")
    print(f"  Free after fixing inter-nibble: 512 - {rank} = {512 - rank}")
    print(f"  With inter-nibble carries known: problem decomposes into 8 independent 4-bit subsystems")


def axis_D_filter():
    """
    AXIS D: SHA-256 as a driven filter.

    a[r+1] = F(a[r..r-7]) + W[r]  (simplifying)

    W[r] is the "driving signal". For r < 16: W[r] = M[r].
    For r ≥ 16: W[r] = schedule(M[0..15]).

    The TRANSFER FUNCTION H(z) of this filter:
    In Z-transform: A(z) = F(z) · A(z) + W(z) → A(z) = W(z) / (1 - F(z))

    But this is for LINEAR systems. SHA-256 is nonlinear.

    Still: the W signal enters the system ADDITIVELY (mod 2^32).
    This means: a_new - a_new_with_W=0 = contribution of W[r] to output.

    How much of the output H is "directly from W" vs "from state evolution"?
    """
    print("\n" + "=" * 80)
    print("AXIS D: SHA-256 as driven filter — W contribution to output")
    print("=" * 80)

    N = 500
    # Compute SHA-256 with real W and with W=0 (but keep K)
    # The difference = W's contribution

    hw_diffs = []
    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]

        # Real SHA
        states_real, _ = sha256_round_trace(M)

        # SHA with W[r] = 0 for all r (but K[r] remains)
        M_zero = [0] * 16
        states_zero, _ = sha256_round_trace(M_zero)

        # Difference at output
        diff = 0
        for i in range(8):
            diff += bin(states_real[64][i] ^ states_zero[64][i]).count('1')
        hw_diffs.append(diff)

    avg_diff = sum(hw_diffs) / N
    print(f"  E[HW(SHA(M) XOR SHA(0))] at round 64 = {avg_diff:.1f}/256")
    print(f"  (If W contributed nothing: would be 0. If everything: 128)")

    # Per round: how quickly does the W-signal spread?
    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]
    M_zero = [0] * 16

    states_M, _ = sha256_round_trace(M)
    states_0, _ = sha256_round_trace(M_zero)

    print(f"\n  W-signal spread per round (HW of state difference):")
    for r in [1, 2, 4, 8, 16, 32, 64]:
        diff = sum(bin(states_M[r][i] ^ states_0[r][i]).count('1') for i in range(8))
        print(f"    r={r:>2}: HW = {diff}/256")


def axis_E_psi():
    """
    AXIS E: ψ[r] = a[r] - e[r] mod 2^32.

    From створочное число: ψ[r] = Sig0(a[r-1]) + Maj(a[r-1], a[r-2], a[r-3]) - a[r-4]
    = T2[r-1] - a[r-4]

    ψ depends only on a-trajectory, NOT on M or K or W.
    It's a PURE structural function.

    ψ[r+1] - ψ[r] = ?
    ψ[r+1] = T2[r] - a[r-3]
    ψ[r]   = T2[r-1] - a[r-4]
    Δψ     = T2[r] - T2[r-1] - a[r-3] + a[r-4]
            = (T2[r] - T2[r-1]) - (a[r-3] - a[r-4])
            = ΔT2[r] - Δa[r-4]

    where ΔT2[r] = Sig0(a[r]) + Maj(a[r], a[r-1], a[r-2]) - Sig0(a[r-1]) - Maj(a[r-1], a[r-2], a[r-3])

    Question: Does ψ have any regularity? Any pattern that persists?
    """
    print("\n" + "=" * 80)
    print("AXIS E: ψ[r] = a[r] - e[r] — pure structural function")
    print("=" * 80)

    N = 200
    psi_hws = [[] for _ in range(65)]

    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        states, _ = sha256_round_trace(M)

        for r in range(65):
            a_val = states[r][0]
            e_val = states[r][4]
            psi = (a_val - e_val) & MASK
            psi_hws[r].append(bin(psi).count('1'))

    print(f"  E[HW(ψ[r])] across rounds (N={N}):")
    for r in [0, 1, 2, 4, 8, 16, 32, 64]:
        avg = sum(psi_hws[r]) / N
        std = (sum((x - avg)**2 for x in psi_hws[r]) / N) ** 0.5
        print(f"    r={r:>2}: E[HW]={avg:.2f}, σ={std:.2f}")

    # Is ψ "simpler" than a or e individually?
    # Check: bit-level bias of ψ
    print(f"\n  Bit-level bias of ψ[16] across {N} messages:")
    max_bias = 0
    for bit in range(32):
        ones = 0
        for seed in range(N):
            rng = random.Random(seed)
            M = [rng.randint(0, MASK) for _ in range(16)]
            states, _ = sha256_round_trace(M, rounds=16)
            psi = (states[16][0] - states[16][4]) & MASK
            if (psi >> bit) & 1:
                ones += 1
        bias = abs(ones/N - 0.5)
        if bias > max_bias:
            max_bias = bias

    print(f"  Max bit bias of ψ[16]: {max_bias:.4f} (random < 0.05)")

    # ψ at round 0 is FIXED (from IV)
    psi_0 = (H0[0] - H0[4]) & MASK
    print(f"\n  ψ[0] = a_iv - e_iv = 0x{psi_0:08x} (CONSTANT)")
    print(f"  ψ[0] HW = {bin(psi_0).count('1')}")

    # Does ψ return close to its initial value? (periodic behavior?)
    print(f"\n  HW(ψ[r] XOR ψ[0]) — distance from initial ψ:")
    for seed in range(3):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        states, _ = sha256_round_trace(M)
        dists = []
        for r in range(65):
            psi_r = (states[r][0] - states[r][4]) & MASK
            dists.append(bin(psi_r ^ psi_0).count('1'))
        print(f"    Seed {seed}: {[dists[r] for r in [0,1,4,8,16,32,64]]}")


if __name__ == "__main__":
    axis_C_nibble()
    axis_D_filter()
    axis_E_psi()
