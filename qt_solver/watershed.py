"""
Direction 4: Watershed manipulation analysis.

From BTE theory (SESSION_RESULTS.md):
  - Watershed round ≈ 57 (where collision paths converge)
  - Minimum watershed ≈ 48 (from schedule capacity: 32r*-1536 > 0 → r* > 48)
  - Before watershed: paths INDEPENDENT (P(match) ≈ 0)
  - At watershed: paths CONVERGE (sudden junction)
  - After watershed: paths IDENTICAL

Watershed = the boundary between "two rivers" flowing into "one lake".

Questions:
  1. Can we push watershed earlier than round 57?
  2. What is the actual convergence behavior for different message differences?
  3. Does the schedule kernel (160 free bits for δW[53..63]=0) help?
  4. Is there structure in which message differences converge faster?

Approach: Wang chain + backward extension analysis.
"""

from qt_solver.sha256_traced import (
    MASK32, sha256_compress, sha256_compress_traced, IV,
    sigma0, sigma1, ch, maj, get_bit,
)
import random


def find_convergence_round(msg1, msg2, R=64):
    """
    Find the round where two messages' states start matching.
    Returns the earliest round r where state[r:] are identical.
    """
    t1 = sha256_compress_traced(msg1, R)
    t2 = sha256_compress_traced(msg2, R)

    # Find last round where states differ
    last_diff = -1
    for r in range(R + 1):
        if t1['states'][r] != t2['states'][r]:
            last_diff = r

    if last_diff == -1:
        return 0  # same from the start (identical messages or same hash)

    # Convergence round = first round after last_diff where states match
    conv_round = last_diff + 1 if last_diff < R else None

    return {
        'last_diff_round': last_diff,
        'convergence_round': conv_round,
        'total_rounds': R,
    }


def analyze_state_convergence(msg1, msg2, R=64):
    """
    Detailed analysis of how two messages' states converge.
    """
    t1 = sha256_compress_traced(msg1, R)
    t2 = sha256_compress_traced(msg2, R)

    print(f"State convergence analysis (R={R}):")
    print(f"  M1: {[hex(w) for w in msg1[:4]]}...")
    print(f"  M2: {[hex(w) for w in msg2[:4]]}...")
    print()

    # Per-register difference at each round
    reg_names = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
    for r in range(R + 1):
        s1 = t1['states'][r]
        s2 = t2['states'][r]
        diffs = []
        total_hw = 0
        for i in range(8):
            hw = bin(s1[i] ^ s2[i]).count('1')
            total_hw += hw
            if hw > 0:
                diffs.append(f"{reg_names[i]}:{hw}")

        if total_hw > 0 and (r < 5 or r > R - 10 or total_hw < 20):
            print(f"  Round {r:2d}: diff={total_hw:3d} bits  [{', '.join(diffs)}]")
        elif total_hw == 0:
            print(f"  Round {r:2d}: MATCH (converged)")

    # Schedule difference
    print(f"\n  Schedule difference:")
    w1 = t1['W']
    w2 = t2['W']
    zero_from = None
    for i in range(min(len(w1), len(w2))):
        hw = bin(w1[i] ^ w2[i]).count('1')
        if hw == 0 and zero_from is None:
            zero_from = i
        elif hw > 0:
            zero_from = None
        if i < R and (i < 5 or i > R - 5):
            print(f"    W[{i:2d}]: diff={hw:2d} bits")

    if zero_from is not None:
        print(f"    δW=0 from round {zero_from}")


def schedule_kernel_analysis():
    """
    Analyze the schedule kernel: which message differences make δW[r..63]=0?

    If δW[r..63]=0 for some r, then rounds r..63 see identical W values,
    which means state convergence is purely driven by the state difference
    (no new divergence from schedule).
    """
    print("=" * 60)
    print("Schedule Kernel Analysis")
    print("=" * 60)

    # Message schedule is LINEAR over GF(2) (rotations + XOR)
    # But NONLINEAR over Z/2^32 (additions create carry)
    # Over GF(2): W[i] = σ₁(W[i-2]) ⊕ W[i-7] ⊕ σ₀(W[i-15]) ⊕ W[i-16]
    # XOR-linear schedule has large kernel.
    # But actual schedule uses + (additions with carry).

    # For the ACTUAL schedule: compute how many message differences
    # lead to δW[r]=0 for late rounds.

    rng = random.Random(42)
    msg = [rng.randint(0, MASK32) for _ in range(16)]

    # Check: for single-bit message differences, when does δW become 0?
    print("\nSingle-bit message changes → δW propagation:")
    for wi in range(16):
        for bi in [0, 15, 31]:  # sample bits
            msg2 = list(msg)
            msg2[wi] ^= (1 << bi)

            w1 = _expand_schedule(msg)
            w2 = _expand_schedule(msg2)

            last_nonzero = -1
            for r in range(64):
                if w1[r] != w2[r]:
                    last_nonzero = r

            if bi == 0:
                print(f"  W[{wi}] bit {bi}: δW nonzero through round {last_nonzero}")


def _expand_schedule(msg):
    """Expand message schedule to 64 words."""
    from qt_solver.sha256_traced import ssigma0, ssigma1
    w = list(msg[:16])
    for i in range(16, 64):
        w.append((ssigma1(w[i-2]) + w[i-7] + ssigma0(w[i-15]) + w[i-16]) & MASK32)
    return w


def wang_chain_analysis(R=64, seeds=5):
    """
    Wang chain: sequential additive Newton in Z/2^32.
    Can ensure de[2..16]=0 → partial collision H[4..7].

    Analyze how far Wang chain can extend convergence.
    """
    print("=" * 60)
    print("Wang Chain — Partial Collision Analysis")
    print("=" * 60)

    rng = random.Random(42)

    for seed in range(seeds):
        msg = [rng.randint(0, MASK32) for _ in range(16)]
        h = sha256_compress(msg, R)

        # Wang chain: perturb W[0] to kill e-difference at rounds 2..16
        # e_new = d + T1. If we can make T1_diff = -d_diff, then e_new_diff = 0.
        # T1 depends on W[r], so adjusting W[r] can control T1.

        # Simple test: flip 1 bit in W[0], observe e-convergence
        msg2 = list(msg)
        msg2[0] ^= 1  # flip bit 0 of W[0]

        t1 = sha256_compress_traced(msg, R)
        t2 = sha256_compress_traced(msg2, R)

        # Track e-register difference
        e_diffs = []
        a_diffs = []
        for r in range(min(R+1, len(t1['states']))):
            s1 = t1['states'][r]
            s2 = t2['states'][r]
            e_hw = bin(s1[4] ^ s2[4]).count('1')  # e = index 4
            a_hw = bin(s1[0] ^ s2[0]).count('1')  # a = index 0
            e_diffs.append(e_hw)
            a_diffs.append(a_hw)

        if seed == 0:
            print(f"\nSeed {seed}: flip W[0] bit 0")
            print(f"  e-branch convergence:")
            for r in range(min(20, len(e_diffs))):
                bar_e = '#' * min(e_diffs[r], 40)
                bar_a = '#' * min(a_diffs[r], 40)
                print(f"    R={r:2d}: δe={e_diffs[r]:2d} {bar_e}")
            print(f"  a-branch (harder to control):")
            for r in range(min(20, len(a_diffs))):
                bar_a = '#' * min(a_diffs[r], 40)
                print(f"    R={r:2d}: δa={a_diffs[r]:2d} {bar_a}")


def backward_extension_analysis():
    """
    створочне backward: from matching hash at round R, extend backward.

    If hash matches at round R: a[R], e[R] (and shifted copies) are known.
    Can we recover a[R-1], e[R-1] from this?

    a[R] = T1 + T2 = known → T1 = known - T2
    e[R] = d[R-1] + T1 = known → T1 = e[R] - d[R-1] = e[R] - a[R-4]

    With SHA-256 register structure, backward is partially deterministic.
    """
    print("=" * 60)
    print("Backward Extension from Matching Hash")
    print("=" * 60)

    rng = random.Random(42)
    msg = [rng.randint(0, MASK32) for _ in range(16)]

    R = 64
    trace = sha256_compress_traced(msg, R)
    states = trace['states']

    # From final state, try backward computation
    # state[R] = (a[R], a[R-1], a[R-2], a[R-3], e[R], e[R-1], e[R-2], e[R-3])
    # We know state[R] → know a[R], a[R-1], a[R-2], a[R-3], e[R], e[R-1], e[R-2], e[R-3]
    #
    # From e[R] = a[R-3-1] + T1[R-1]:
    #   T1[R-1] = e[R] - a[R-4]  ... but a[R-4] is NOT in state[R]!
    #   a[R-4] = d[R-1] = state[R-1][3] — we don't know state[R-1]

    # Actually: state[R] gives us a[R..R-3] and e[R..R-3].
    # But e[R-3] = state[R][7] = h[R], which we know.
    # And a[R-3] = state[R][3] = d[R], which we know.
    #
    # The створочне: e[r] = a[r] + a[r-4] - T2[r-1]
    # From state[R]: we know a[R..R-3] and e[R..R-3].
    # So: e[R] = a[R] + a[R-4] - T2[R-1]
    # We know e[R] and a[R], but NOT a[R-4] or T2[R-1].
    #
    # However: a[R-4] is NOT in state[R] for general R.
    # For backward from hash: we need to go step by step.

    # Key insight from F62: minimum watershed = round 48
    # Because schedule kernel (δW[53..63]=0) has 160 free bits.
    # To push watershed to round r*, need: state at r* to match,
    # which requires 32*r* - 1536 > 0 → r* > 48.

    print(f"\n  From BTE theory:")
    print(f"    Watershed minimum: round 48 (schedule capacity constraint)")
    print(f"    Default watershed: round 57")
    print(f"    Schedule kernel δW[53..63]=0: 160 free message bits")
    print(f"")
    print(f"  To push watershed to round r*:")
    print(f"    Need: state convergence at round r*")
    print(f"    Need: δW[r*..63] = 0 (schedule agrees)")
    print(f"    Schedule constraint: 32*(64-r*) equations on 512 msg bits")
    print(f"    Free bits = 512 - 32*(64-r*) = 32*r* - 1536")
    print(f"")
    for r_star in [48, 50, 52, 54, 56, 57]:
        free = 32 * r_star - 1536
        print(f"    r*={r_star}: free bits = {free}, cost ≈ 2^{max(0, 256-free)}")


if __name__ == '__main__':
    schedule_kernel_analysis()
    print()
    wang_chain_analysis()
    print()
    backward_extension_analysis()
