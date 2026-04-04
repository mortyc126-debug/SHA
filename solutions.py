"""
SECOND THEOREM CANDIDATE: Number of solutions per layer.

For Layer 0 (bit-0 system): 127 constraints on 512 message bits.
Linear theory: solution space = 2^{512-127} = 2^{385}.
But Layer 0 is NOT linear — it's QUADRATIC (Ch, Maj).

Question: how many ACTUAL solutions does the bit-0 system have?
If = 2^385 → quadratic terms are "invisible" (same as linear)
If < 2^385 → quadratic terms REDUCE solutions → system harder
If > 2^385 → impossible for additional constraints

We can test this on small BTE where we can enumerate.

For BTE-8 (n=8, R=16):
  Layer 0: rank 31, so linear solution space = 2^{128-31} = 2^97.
  Bit-0 of hash = 8 bits. For a FIXED target: 2^{128-8} = 2^120 messages.
  Among those: how many satisfy the full bit-0 TRAJECTORY?

Wait — for hash-only we have 8 equations (not 31).
For trajectory we have 31 equations.
We can't enumerate 2^128 or 2^97. Too big even for n=8.

Need SMALLER: n=4, R=8. Total bits = n × n_msg = 4 × 8 = 32.
Layer 0 rank = 2R-1 = 15. Solution space = 2^{32-15} = 2^17 = 131072.
Small enough to enumerate!
"""

import random
import time


def rotr_n(x, k, n):
    mask = (1 << n) - 1
    return ((x >> k) | (x << (n - k))) & mask


def bte_trajectory(msg, n, n_msg, R, sig0_rots, sig1_rots):
    """Compute full trajectory of BTE."""
    mask = (1 << n) - 1

    IV = [random.Random(42 + i).randint(0, mask) for i in range(8)]
    K = [random.Random(99999 + i).randint(0, mask) for i in range(R)]

    W = list(msg[:n_msg])
    for i in range(n_msg, R):
        idx = max(0, i - 7) if i >= 7 else 0
        W.append((W[i-2] ^ W[idx] ^ W[max(0, i-n_msg)]) & mask)

    def sig0(x):
        r = 0
        for rot in sig0_rots:
            r ^= rotr_n(x, rot, n)
        return r

    def sig1(x):
        r = 0
        for rot in sig1_rots:
            r ^= rotr_n(x, rot, n)
        return r

    def ch(e, f, g):
        return (e & f) ^ (~e & g) & mask

    def maj(a, b, c):
        return (a & b) ^ (a & c) ^ (b & c)

    state = tuple(IV)
    states = [state]

    for r in range(R):
        a, b, c, d, e, f, g, h = state
        T1 = (h + sig1(e) + ch(e, f, g) + K[r] + W[r]) & mask
        T2 = (sig0(a) + maj(a, b, c)) & mask
        state = ((T1 + T2) & mask, a, b, c, (d + T1) & mask, e, f, g)
        states.append(state)

    return states


def get_bit0_trajectory(states, R):
    """Extract bit 0 of registers a and e at each round."""
    bits = []
    for r in range(1, R + 1):
        bits.append(states[r][0] & 1)  # a[r] bit 0
        bits.append(states[r][4] & 1)  # e[r] bit 0
    return tuple(bits)


def get_bit0_hash(states, R):
    """Extract bit 0 of final hash (8 registers)."""
    return tuple((states[R][reg] & 1) for reg in range(8))


def experiment_enumerate_solutions():
    """
    For BTE-4 (n=4): enumerate ALL 2^32 messages.
    For each: compute bit-0 trajectory.
    Count: how many messages give each bit-0 trajectory?

    Expected from linear theory: 2^{32-15} = 2^17 = 131072 per trajectory.
    If actual ≈ 131072 → quadratic terms don't reduce solutions.
    If actual < 131072 → quadratic terms matter → STRUCTURE.
    """
    print("=" * 80)
    print("ENUMERATE: BTE-4 solution counts per bit-0 trajectory")
    print("=" * 80)

    n = 4
    n_msg = 8  # 8 message words
    R = 8
    sig0_rots = [1, 2]
    sig1_rots = [1, 3]

    total_bits = n * n_msg  # 32
    expected_traj_rank = 2 * R - 1  # 15
    expected_sol_count = 2 ** (total_bits - expected_traj_rank)  # 2^17

    print(f"  BTE-4: n={n}, n_msg={n_msg}, R={R}")
    print(f"  Total message bits: {total_bits}")
    print(f"  Expected trajectory rank: {expected_traj_rank}")
    print(f"  Expected solutions per trajectory: 2^{total_bits - expected_traj_rank} = {expected_sol_count}")

    # Enumerate all 2^32 messages
    total_msgs = 2 ** total_bits
    print(f"  Enumerating {total_msgs} messages...")

    t0 = time.time()

    # Group messages by bit-0 trajectory
    traj_counts = {}

    for msg_int in range(total_msgs):
        # Decode integer to message words
        msg = []
        tmp = msg_int
        for _ in range(n_msg):
            msg.append(tmp & ((1 << n) - 1))
            tmp >>= n

        states = bte_trajectory(msg, n, n_msg, R, sig0_rots, sig1_rots)
        traj = get_bit0_trajectory(states, R)

        if traj not in traj_counts:
            traj_counts[traj] = 0
        traj_counts[traj] += 1

    t1 = time.time()

    n_trajs = len(traj_counts)
    counts = list(traj_counts.values())
    min_c = min(counts)
    max_c = max(counts)
    avg_c = sum(counts) / len(counts)

    print(f"\n  Enumeration complete ({t1-t0:.1f}s)")
    print(f"  Distinct bit-0 trajectories: {n_trajs}")
    print(f"  Expected: 2^{expected_traj_rank} = {2**expected_traj_rank}")
    print(f"  Solutions per trajectory: min={min_c}, max={max_c}, avg={avg_c:.1f}")
    print(f"  Expected (linear): {expected_sol_count}")

    if min_c == max_c:
        print(f"  *** ALL trajectories have EXACTLY {min_c} solutions ***")
        print(f"  This means: quadratic terms DON'T reduce solution count.")
        print(f"  The system is EFFECTIVELY LINEAR in terms of solution count.")
    else:
        print(f"  Trajectories have DIFFERENT solution counts!")
        print(f"  Ratio max/min = {max_c/min_c:.2f}")

        # Distribution
        from collections import Counter
        count_dist = Counter(counts)
        print(f"  Distribution of solution counts:")
        for c, freq in sorted(count_dist.items()):
            print(f"    {c} solutions: {freq} trajectories")

    # Also: group by bit-0 HASH (just final state bit 0)
    print(f"\n  --- By HASH bit-0 (8 bits) ---")
    hash_counts = {}
    for msg_int in range(min(total_msgs, 1000000)):
        msg = []
        tmp = msg_int
        for _ in range(n_msg):
            msg.append(tmp & ((1 << n) - 1))
            tmp >>= n

        states = bte_trajectory(msg, n, n_msg, R, sig0_rots, sig1_rots)
        h = get_bit0_hash(states, R)

        if h not in hash_counts:
            hash_counts[h] = 0
        hash_counts[h] += 1

    n_hashes = len(hash_counts)
    h_counts = list(hash_counts.values())
    h_min = min(h_counts)
    h_max = max(h_counts)
    h_avg = sum(h_counts) / len(h_counts)
    expected_per_hash = total_msgs / (2 ** 8)

    print(f"  Distinct hash bit-0 values: {n_hashes}/256")
    print(f"  Solutions per hash: min={h_min}, max={h_max}, avg={h_avg:.0f}")
    print(f"  Expected: {expected_per_hash:.0f}")

    if h_min == h_max:
        print(f"  *** PERFECTLY BALANCED: every hash has {h_min} preimages ***")
    else:
        ratio = h_max / h_min if h_min > 0 else float('inf')
        print(f"  Ratio max/min = {ratio:.3f}")
        if ratio < 1.1:
            print(f"  Nearly balanced (< 10% variation)")
        else:
            print(f"  *** IMBALANCED: some hashes have {ratio:.1f}× more preimages ***")


if __name__ == "__main__":
    experiment_enumerate_solutions()
