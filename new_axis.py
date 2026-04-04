"""
NEW AXIS: Built from scratch. Not GF(2), not Z/2^32, not lattice.

Core idea: the fundamental object is NOT the value and NOT the difference.
It's the VALUE + HOW IT WAS PRODUCED.

Two values can be numerically identical but have different "production history."
In SHA-256, this history = which carries happened, which rotations applied.

Define: a SHA-ELEMENT is a triple (v, g, p) where:
  v = 32-bit value (the number itself)
  g = 32-bit "generate mask" (which positions generated carry last time)
  p = 32-bit "propagate mask" (which positions propagated carry last time)

When we add two SHA-elements:
  (v1, g1, p1) + (v2, g2, p2) = (v1+v2, new_g, new_p)
  where new_g[k] = v1[k] AND v2[k]  (generate at this addition)
        new_p[k] = v1[k] XOR v2[k]  (propagate at this addition)

When we XOR two SHA-elements:
  (v1, g1, p1) XOR (v2, g2, p2) = (v1^v2, 0, v1^v2)
  XOR generates no carries (g=0) and everything "propagates" (p = result)

When we rotate:
  ROT((v, g, p), k) = (ROT(v,k), ROT(g,k), ROT(p,k))
  Rotation preserves the GPK structure.

The GPK triple carries MORE information than v alone.
Specifically: g tells you WHERE the nonlinearity concentrated.
              p tells you WHERE the result is "carry-transparent."
              ~g & ~p (= kill) tells you WHERE carries were blocked.

This is our new algebra. Let's test if it reveals anything.
"""

import random
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, schedule
)


class ShaElement:
    """A value with its GPK (generate-propagate-kill) production history."""

    __slots__ = ['v', 'g', 'p']

    def __init__(self, v, g=0, p=None):
        self.v = v & MASK
        self.g = g & MASK
        self.p = (p if p is not None else v) & MASK  # default: fresh value

    def __add__(self, other):
        """Mod 2^32 addition with GPK tracking."""
        result = (self.v + other.v) & MASK
        new_g = self.v & other.v  # bit positions where both are 1
        new_p = self.v ^ other.v  # bit positions where exactly one is 1
        return ShaElement(result, new_g, new_p)

    def __xor__(self, other):
        """XOR with GPK tracking."""
        result = self.v ^ other.v
        new_g = 0  # XOR never generates carry
        new_p = result  # XOR result = propagate pattern
        return ShaElement(result, new_g, new_p)

    def __and__(self, other):
        """AND with GPK tracking."""
        result = self.v & other.v
        return ShaElement(result, result, 0)  # AND = all generate, no propagate

    def __invert__(self):
        """NOT."""
        return ShaElement(~self.v & MASK, ~self.g & MASK, ~self.p & MASK)

    def rotr(self, k):
        """Rotate right by k. Preserves GPK structure."""
        return ShaElement(rotr(self.v, k), rotr(self.g, k), rotr(self.p, k))

    @property
    def kill(self):
        """Kill mask: positions where carry is blocked."""
        return (~self.g & ~self.p) & MASK

    def carry_chain(self):
        """Compute actual carry chain from g and p."""
        c = 0
        carry_word = 0
        for k in range(32):
            carry_word |= (c << k)
            gk = (self.g >> k) & 1
            pk = (self.p >> k) & 1
            c = gk | (pk & c)
        return carry_word

    def __repr__(self):
        return f'SE(v=0x{self.v:08x}, g={bin(self.g).count("1"):>2}g, p={bin(self.p).count("1"):>2}p, k={bin(self.kill).count("1"):>2}k)'


def sha_element_round(state_se, r, W_val):
    """One round of SHA-256 using ShaElements."""
    a, b, c, d, e, f, g, h = state_se

    # Sig1(e) = ROTR(e,6) XOR ROTR(e,11) XOR ROTR(e,25)
    sig1_e = e.rotr(6) ^ e.rotr(11) ^ e.rotr(25)

    # Ch(e,f,g) = (e AND f) XOR (NOT e AND g)
    ch = (e & f) ^ (~e & g)

    # T1 = h + Sig1(e) + Ch(e,f,g) + K[r] + W[r]
    K_se = ShaElement(K[r])
    W_se = ShaElement(W_val)
    T1 = h + sig1_e + ch + K_se + W_se

    # Sig0(a) = ROTR(a,2) XOR ROTR(a,13) XOR ROTR(a,22)
    sig0_a = a.rotr(2) ^ a.rotr(13) ^ a.rotr(22)

    # Maj(a,b,c) = (a AND b) XOR (a AND c) XOR (b AND c)
    maj = (a & b) ^ (a & c) ^ (b & c)

    # T2 = Sig0(a) + Maj(a,b,c)
    T2 = sig0_a + maj

    # New state
    a_new = T1 + T2
    e_new = d + T1
    return [a_new, a, b, c, e_new, e, f, g]


def experiment_sha_element():
    """Run SHA-256 with ShaElement tracking and examine GPK evolution."""
    print("=" * 80)
    print("NEW AXIS: SHA-Element (value + GPK production history)")
    print("=" * 80)

    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]
    W = schedule(M)

    # Initialize state with fresh ShaElements
    state = [ShaElement(h) for h in H0]

    # Verify: values match standard SHA-256
    states_check, _ = sha256_round_trace(M)

    print(f"\n  Tracking GPK through 64 rounds:")
    print(f"  {'r':>3} | {'a_val':>10} | {'g_bits':>6} {'p_bits':>6} {'k_bits':>6} | {'check':>5} | {'carry_HW':>9}")
    print(f"  " + "-" * 65)

    for r in range(64):
        state = sha_element_round(state, r, W[r])

        a_se = state[0]
        a_check = states_check[r+1][0]

        # Verify value matches
        match = "OK" if a_se.v == a_check else "FAIL"

        g_bits = bin(a_se.g).count('1')
        p_bits = bin(a_se.p).count('1')
        k_bits = bin(a_se.kill).count('1')
        carry_hw = bin(a_se.carry_chain()).count('1')

        if r < 8 or r >= 60 or r % 16 == 0:
            print(f"  {r:>3} | 0x{a_se.v:08x} | {g_bits:>6} {p_bits:>6} {k_bits:>6} | {match:>5} | {carry_hw:>9}")

    # Final state
    print(f"\n  Final state (round 64):")
    reg_names = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
    for i, (name, se) in enumerate(zip(reg_names, state)):
        print(f"    {name}: {se}")


def experiment_gpk_statistics():
    """
    Statistics of GPK across many messages.

    KEY QUESTION: Is GPK at round 64 more structured than the value alone?
    If HW(g), HW(p), HW(k) have DIFFERENT distributions than random
    → GPK carries structural information.
    """
    print("\n" + "=" * 80)
    print("GPK STATISTICS: Does GPK carry extra info beyond value?")
    print("=" * 80)

    N = 500
    g_hws = []
    p_hws = []
    k_hws = []
    v_hws = []

    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        W = schedule(M)

        state = [ShaElement(h) for h in H0]
        for r in range(64):
            state = sha_element_round(state, r, W[r])

        a_final = state[0]
        g_hws.append(bin(a_final.g).count('1'))
        p_hws.append(bin(a_final.p).count('1'))
        k_hws.append(bin(a_final.kill).count('1'))
        v_hws.append(bin(a_final.v).count('1'))

    avg_g = sum(g_hws) / N
    avg_p = sum(p_hws) / N
    avg_k = sum(k_hws) / N
    avg_v = sum(v_hws) / N

    std_g = (sum((x-avg_g)**2 for x in g_hws)/N)**0.5
    std_p = (sum((x-avg_p)**2 for x in p_hws)/N)**0.5
    std_k = (sum((x-avg_k)**2 for x in k_hws)/N)**0.5
    std_v = (sum((x-avg_v)**2 for x in v_hws)/N)**0.5

    print(f"  E[HW] at round 64 (N={N}):")
    print(f"    value v: {avg_v:.2f} ± {std_v:.2f}  (random = 16.00 ± 2.83)")
    print(f"    gen   g: {avg_g:.2f} ± {std_g:.2f}  (random = 8.00 ± 2.45)")
    print(f"    prop  p: {avg_p:.2f} ± {std_p:.2f}  (random = 16.00 ± 2.83)")
    print(f"    kill  k: {avg_k:.2f} ± {std_k:.2f}  (random = 8.00 ± 2.45)")

    # For random value v: g = v1 AND v2 where v = v1 + v2
    # E[HW(g)] = E[HW(v1 AND v2)] = 8 for random v1, v2
    # E[HW(p)] = E[HW(v1 XOR v2)] = 16 for random v1, v2
    # E[HW(k)] = 32 - 8 - 16 = 8

    # If actual matches random expectations → GPK carries no extra info.
    # If actual differs → GPK has structure!

    g_dev = abs(avg_g - 8.0)
    p_dev = abs(avg_p - 16.0)
    k_dev = abs(avg_k - 8.0)

    print(f"\n  Deviation from random:")
    print(f"    g: {g_dev:.2f} ({'SIGNAL' if g_dev > 0.5 else 'random'})")
    print(f"    p: {p_dev:.2f} ({'SIGNAL' if p_dev > 0.5 else 'random'})")
    print(f"    k: {k_dev:.2f} ({'SIGNAL' if k_dev > 0.5 else 'random'})")


def experiment_gpk_correlation():
    """
    Does GPK at round r correlate with GPK at round r+1?
    (Beyond what value alone shows)

    If GPK-GPK correlation > value-value correlation → GPK carries
    information about the COMPUTATION that value alone doesn't.
    """
    print("\n" + "=" * 80)
    print("GPK CORRELATION: Does GPK track computation beyond values?")
    print("=" * 80)

    N = 300
    # Correlation of g[r] with g[r+1] (HW correlation)
    gg_corr = []
    vv_corr = []

    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        W = schedule(M)

        state = [ShaElement(h) for h in H0]
        prev_g = None
        prev_v = None

        for r in range(64):
            state = sha_element_round(state, r, W[r])
            curr_g = bin(state[0].g).count('1')
            curr_v = bin(state[0].v).count('1')

            if prev_g is not None and r >= 8:
                gg_corr.append(curr_g * prev_g)
                vv_corr.append(curr_v * prev_v)

            prev_g = curr_g
            prev_v = curr_v

    # Pearson correlation
    def pearson(xy_list, x_mean, y_mean):
        n = len(xy_list) // 2  # rough
        return sum(xy_list) / len(xy_list)

    avg_gg = sum(gg_corr) / len(gg_corr)
    avg_vv = sum(vv_corr) / len(vv_corr)

    # For uncorrelated: E[HW*HW] = E[HW]^2 = 16^2 = 256 (value) or 8^2 = 64 (g)
    print(f"  E[HW(g[r]) × HW(g[r+1])] = {avg_gg:.2f}  (uncorrelated = 64)")
    print(f"  E[HW(v[r]) × HW(v[r+1])] = {avg_vv:.2f}  (uncorrelated = 256)")

    g_excess = avg_gg - 64
    v_excess = avg_vv - 256
    print(f"  g correlation excess: {g_excess:.2f}")
    print(f"  v correlation excess: {v_excess:.2f}")
    print(f"  → {'GPK more correlated than values!' if abs(g_excess) > abs(v_excess) else 'No extra correlation in GPK'}")


def experiment_gpk_two_messages():
    """
    THE REAL TEST: For two messages with same hash value,
    do their GPK patterns differ from two messages with different hash?

    If same-hash pairs have more similar GPK → GPK "sees" the collision.
    If GPK is identical for same-hash → GPK is determined by hash alone.
    If GPK is random regardless → GPK carries no collision info.

    We can't find same-hash pairs for 64 rounds.
    But for REDUCED rounds (where we can find collisions via birthday),
    we can test this.
    """
    print("\n" + "=" * 80)
    print("GPK TWO-MESSAGE TEST: Does GPK distinguish collision paths?")
    print("=" * 80)

    # For R=1 round, find pairs with same hash (easy: just 2^32 space)
    # Actually: for reduced SHA with R rounds, fix W[1..15] and vary W[0].
    # Hash = state[R] + IV. Project to 8 bits.
    # Find collision in the 8-bit hash, then compare GPK.

    R = 4
    rng = random.Random(42)
    base_msg = [rng.randint(0, MASK) for _ in range(16)]

    hash_to_gpk = {}
    N = 10000

    for trial in range(N):
        msg = list(base_msg)
        msg[0] = random.Random(trial).randint(0, MASK)
        W = schedule(msg)

        state = [ShaElement(h) for h in H0]
        for r in range(R):
            state = sha_element_round(state, r, W[r])

        hash_val = (state[0].v + H0[0]) & 0xFF  # 8-bit hash
        gpk_info = (state[0].g, state[0].p)  # GPK of final a

        if hash_val not in hash_to_gpk:
            hash_to_gpk[hash_val] = []
        hash_to_gpk[hash_val].append(gpk_info)

    # For same-hash pairs: compare GPK
    same_hash_gpk_diffs = []
    for hash_val, entries in hash_to_gpk.items():
        if len(entries) < 2:
            continue
        for i in range(min(5, len(entries))):
            for j in range(i+1, min(5, len(entries))):
                g_diff = bin(entries[i][0] ^ entries[j][0]).count('1')
                p_diff = bin(entries[i][1] ^ entries[j][1]).count('1')
                same_hash_gpk_diffs.append((g_diff, p_diff))

    # For random pairs: compare GPK
    all_entries = []
    for entries in hash_to_gpk.values():
        all_entries.extend(entries[:3])

    random_gpk_diffs = []
    for i in range(min(1000, len(all_entries))):
        j = (i + 7) % len(all_entries)
        g_diff = bin(all_entries[i][0] ^ all_entries[j][0]).count('1')
        p_diff = bin(all_entries[i][1] ^ all_entries[j][1]).count('1')
        random_gpk_diffs.append((g_diff, p_diff))

    if same_hash_gpk_diffs:
        avg_g_same = sum(d[0] for d in same_hash_gpk_diffs) / len(same_hash_gpk_diffs)
        avg_p_same = sum(d[1] for d in same_hash_gpk_diffs) / len(same_hash_gpk_diffs)
        avg_g_rand = sum(d[0] for d in random_gpk_diffs) / len(random_gpk_diffs)
        avg_p_rand = sum(d[1] for d in random_gpk_diffs) / len(random_gpk_diffs)

        print(f"  R={R} rounds, 8-bit hash collisions:")
        print(f"    Same-hash GPK distance: g={avg_g_same:.2f}, p={avg_p_same:.2f}")
        print(f"    Random GPK distance:    g={avg_g_rand:.2f}, p={avg_p_rand:.2f}")
        print(f"    Ratio (same/random): g={avg_g_same/avg_g_rand:.3f}, p={avg_p_same/avg_p_rand:.3f}")
        print(f"    → {'GPK DISTINGUISHES collision paths!' if avg_g_same/avg_g_rand < 0.9 else 'GPK same as random for collisions'}")


if __name__ == "__main__":
    experiment_sha_element()
    experiment_gpk_statistics()
    experiment_gpk_correlation()
    experiment_gpk_two_messages()
