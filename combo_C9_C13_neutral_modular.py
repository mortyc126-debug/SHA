#!/usr/bin/env python3
"""
C9×C13 (Layer 2): Neutral Bits × Modular Differentials
=======================================================
Combines two Layer 1 results:
  C9:  82.6 neutral bits exist for Wang chain (bits that don't break De_xor=0)
  C13: Modular differentials are deterministic for 1 round (HW=1, entropy=0)

Key idea: Can we SELECT neutral bits so that the MODULAR differential stays
small for MORE than 1 round past the Wang chain boundary (round 16)?

Standard neutral bits preserve XOR differential. But maybe some subset
preserves MODULAR differential better. If we can keep |De_mod| < 2^k for
2-3 extra rounds by choosing the right neutral bits, we get a hybrid advantage.

Experiment:
  1. SHA-256 with differential tracking (XOR and modular)
  2. Wang chain: DW[0] = 0x80000000, force De=0 for rounds 1-16
  3. Find neutral bits (bits in W[1..15] preserving De_xor=0 through round 16)
  4. For each neutral bit, track De_mod at rounds 17-20
  5. Rank neutral bits by |De_mod[17]|
  6. Test pairs of top-ranked neutral bits for constructive/destructive interference
  7. Greedy search: flip neutral bits one at a time to minimize |De_mod[17..20]|
  8. Report statistics across 50 random base messages
"""

import random
import time
import numpy as np
from collections import defaultdict

MASK = 0xFFFFFFFF
MOD = 2**32

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

# SHA-256 initial hash
H0 = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
]

# ── SHA-256 primitives ──────────────────────────────────────────────

def rotr(x, n):
    return ((x >> n) | (x << (32 - n))) & MASK

def Sig0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sig1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def sig0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sig1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)
def Ch(e, f, g):  return ((e & f) ^ ((~e) & g)) & MASK
def Maj(a, b, c): return ((a & b) ^ (a & c) ^ (b & c)) & MASK

def add32(*args):
    s = 0
    for a in args:
        s = (s + a) & MASK
    return s

def sha256_round(state, W_t, K_t):
    a, b, c, d, e, f, g, h = state
    T1 = add32(h, Sig1(e), Ch(e, f, g), K_t, W_t)
    T2 = add32(Sig0(a), Maj(a, b, c))
    new_a = add32(T1, T2)
    new_e = add32(d, T1)
    return [new_a, a, b, c, new_e, e, f, g]

def expand_message(W16):
    """Expand 16-word message to 64 words using SHA-256 schedule."""
    W = list(W16)
    for i in range(16, 64):
        W.append(add32(sig1(W[i-2]), W[i-7], sig0(W[i-15]), W[i-16]))
    return W

def run_rounds(init_state, W, num_rounds):
    """Run num_rounds of SHA-256. Returns list of states [init, after_r0, after_r1, ...]."""
    states = [list(init_state)]
    for t in range(num_rounds):
        s = sha256_round(states[-1], W[t], K[t])
        states.append(s)
    return states

def mod_diff_signed(a, b):
    """Signed modular difference: a - b mod 2^32, mapped to [-2^31, 2^31)."""
    d = (a - b) & MASK
    if d >= (1 << 31):
        return d - MOD
    return d

def abs_mod_diff(a, b):
    """Absolute modular difference."""
    return abs(mod_diff_signed(a, b))

# ── Wang chain ──────────────────────────────────────────────────────

def compute_wang_chain(base_msg, init_state, dw0=0x80000000):
    """
    Build Wang chain: choose DW[t] so De_t = 0 for rounds 1..15.
    DW[0] = dw0 (MSB flip). Returns (W_prime, states, states_prime).
    """
    W = list(base_msg)
    W_prime = list(W)
    W_prime[0] = W[0] ^ dw0

    state = list(init_state)
    state_prime = list(init_state)

    states = [list(state)]
    states_prime = [list(state_prime)]

    # Round 0
    state = sha256_round(state, W[0], K[0])
    state_prime = sha256_round(state_prime, W_prime[0], K[0])
    states.append(list(state))
    states_prime.append(list(state_prime))

    # Rounds 1..15: choose W'[t] to force new_e' = new_e
    for t in range(1, 16):
        a, b, c, d, e, f, g, h = state
        a2, b2, c2, d2, e2, f2, g2, h2 = state_prime

        T1_partial = add32(h, Sig1(e), Ch(e, f, g), K[t])
        T1_partial_prime = add32(h2, Sig1(e2), Ch(e2, f2, g2), K[t])

        target_e = add32(d, T1_partial, W[t])
        # new_e' = d2 + T1_partial_prime + W'[t] = target_e
        W_prime_t = (target_e - d2 - T1_partial_prime) & MASK
        W_prime[t] = W_prime_t

        state = sha256_round(state, W[t], K[t])
        state_prime = sha256_round(state_prime, W_prime_t, K[t])
        states.append(list(state))
        states_prime.append(list(state_prime))

    return W_prime, states, states_prime

def verify_wang_chain(states, states_prime, check_rounds=16):
    """Verify De=0 for states[2..check_rounds]. Returns True if all zero."""
    for t in range(2, min(check_rounds + 1, len(states))):
        de = (states[t][4] - states_prime[t][4]) & MASK
        if de != 0:
            return False
    return True

# ── Neutral bit detection ───────────────────────────────────────────

def find_neutral_bits(base_msg, init_state):
    """
    Find bits in W[1..15] that are neutral for the Wang chain.
    A bit (k, b) is neutral if flipping it in BOTH W and W' preserves De=0
    for rounds 1..15.
    Returns list of (word_idx, bit_idx) that are fully neutral.
    """
    W_prime, states, states_prime = compute_wang_chain(base_msg, init_state)
    neutrals = []

    for k in range(1, 16):
        for b in range(32):
            flip = 1 << b
            W_f = list(base_msg)
            W_f[k] ^= flip
            Wp_f = list(W_prime)
            Wp_f[k] ^= flip

            # Expand and run
            W_exp = expand_message(W_f)
            Wp_exp = expand_message(Wp_f)

            s1 = run_rounds(init_state, W_exp, 16)
            s2 = run_rounds(init_state, Wp_exp, 16)

            # Check De=0 for states[2..16]
            ok = True
            for t in range(2, 17):
                de = (s1[t][4] - s2[t][4]) & MASK
                if de != 0:
                    ok = False
                    break

            if ok:
                neutrals.append((k, b))

    return neutrals

# ── Modular differential tracking ──────────────────────────────────

def track_modular_diff(base_msg, init_state, W_prime_base, neutrals_to_flip):
    """
    Given base message and its Wang-chain partner W_prime_base,
    flip the specified neutral bits in both W and W', expand to 64 words,
    run 21 rounds, and return De_mod (signed) at rounds 17-20.
    """
    W = list(base_msg)
    Wp = list(W_prime_base)

    for (k, b) in neutrals_to_flip:
        flip = 1 << b
        W[k] ^= flip
        Wp[k] ^= flip

    W_exp = expand_message(W)
    Wp_exp = expand_message(Wp)

    s1 = run_rounds(init_state, W_exp, 21)
    s2 = run_rounds(init_state, Wp_exp, 21)

    de_mod = []
    for t in range(17, 22):  # rounds 17, 18, 19, 20, 21 (indices 17-21 in states)
        if t < len(s1):
            de = mod_diff_signed(s1[t][4], s2[t][4])
            de_mod.append(de)
        else:
            de_mod.append(None)

    return de_mod

# ── Main experiment ─────────────────────────────────────────────────

def run_experiment(num_base_msgs=50, top_n=30):
    random.seed(0xDEAD)
    start_time = time.time()

    print("=" * 72)
    print("C9×C13 (Layer 2): Neutral Bits × Modular Differentials")
    print("=" * 72)
    print()
    print(f"Parameters: {num_base_msgs} base messages, top {top_n} neutral bits for pairs")
    print()

    # Collect statistics across all base messages
    all_single_de17 = []           # |De_mod[17]| for each (msg, neutral_bit)
    all_greedy_best_17 = []        # best |De_mod[17]| from greedy search per msg
    all_greedy_best_18 = []
    all_base_de17 = []             # |De_mod[17]| with NO neutral flips
    neutral_bit_scores = defaultdict(list)  # (k,b) -> list of |De_mod[17]| values
    neutral_bit_helpful = defaultdict(int)  # (k,b) -> count of times it helped

    pair_constructive = 0
    pair_destructive = 0
    pair_neutral_count = 0
    pair_total = 0

    greedy_flips_history = []  # how many flips greedy used per message

    for msg_idx in range(num_base_msgs):
        elapsed = time.time() - start_time
        if elapsed > 280:
            print(f"\n[Time limit approaching at message {msg_idx}, stopping early]")
            break

        if msg_idx % 10 == 0:
            print(f"  Message {msg_idx}/{num_base_msgs} (elapsed: {elapsed:.1f}s)")

        base_msg = [random.getrandbits(32) for _ in range(16)]
        init_state = list(H0)

        # Compute Wang chain
        W_prime, states, states_prime = compute_wang_chain(base_msg, init_state)

        # Verify Wang chain
        if not verify_wang_chain(states, states_prime, 16):
            print(f"    WARNING: Wang chain failed for message {msg_idx}")
            continue

        # Baseline: De_mod with no neutral flips at rounds 17-20
        W_exp = expand_message(base_msg)
        Wp_exp = expand_message(W_prime)
        s1_base = run_rounds(init_state, W_exp, 21)
        s2_base = run_rounds(init_state, Wp_exp, 21)

        base_de_mod = []
        for t in range(17, 22):
            if t < len(s1_base):
                base_de_mod.append(mod_diff_signed(s1_base[t][4], s2_base[t][4]))
        all_base_de17.append(abs(base_de_mod[0]))

        # Find neutral bits for this message
        neutrals = find_neutral_bits(base_msg, init_state)

        if len(neutrals) == 0:
            continue

        # ── Step 4: Track De_mod for each single neutral bit flip ──
        single_scores = []  # (abs_de17, (k, b), de_mod_list)
        for (k, b) in neutrals:
            de_mod = track_modular_diff(base_msg, init_state, W_prime, [(k, b)])
            abs_de17 = abs(de_mod[0])
            single_scores.append((abs_de17, (k, b), de_mod))
            all_single_de17.append(abs_de17)
            neutral_bit_scores[(k, b)].append(abs_de17)

        # ── Step 5: Rank by |De_mod[17]| ──
        single_scores.sort(key=lambda x: x[0])

        # Track which bits are "helpful" (reduce |De_mod[17]| vs baseline)
        for abs_de17, (k, b), _ in single_scores:
            if abs_de17 < abs(base_de_mod[0]):
                neutral_bit_helpful[(k, b)] += 1

        # ── Step 6: Pairwise test of top-ranked neutral bits ──
        top_bits = [s[1] for s in single_scores[:top_n]]
        top_scores_dict = {s[1]: s[0] for s in single_scores[:top_n]}

        n_pairs_this = min(len(top_bits) * (len(top_bits) - 1) // 2, 200)
        pair_count = 0
        for i in range(len(top_bits)):
            if pair_count >= n_pairs_this:
                break
            for j in range(i + 1, len(top_bits)):
                if pair_count >= n_pairs_this:
                    break
                b1, b2 = top_bits[i], top_bits[j]
                de_mod_pair = track_modular_diff(base_msg, init_state, W_prime, [b1, b2])
                abs_de17_pair = abs(de_mod_pair[0])

                # Compare pair to best individual
                best_individual = min(top_scores_dict[b1], top_scores_dict[b2])
                worst_individual = max(top_scores_dict[b1], top_scores_dict[b2])

                pair_total += 1
                if abs_de17_pair < best_individual:
                    pair_constructive += 1
                elif abs_de17_pair > worst_individual:
                    pair_destructive += 1
                else:
                    pair_neutral_count += 1
                pair_count += 1

        # ── Step 7: Greedy search ──
        # Start with base message, greedily flip neutral bits to minimize |De_mod[17]|
        current_flips = []
        current_de17 = abs(base_de_mod[0])
        available = list(neutrals)

        greedy_steps = 0
        max_greedy_steps = min(len(available), 40)

        for step in range(max_greedy_steps):
            best_bit = None
            best_de17 = current_de17

            for (k, b) in available:
                test_flips = current_flips + [(k, b)]
                de_mod = track_modular_diff(base_msg, init_state, W_prime, test_flips)
                abs_de17 = abs(de_mod[0])
                if abs_de17 < best_de17:
                    best_de17 = abs_de17
                    best_bit = (k, b)

            if best_bit is not None:
                current_flips.append(best_bit)
                available.remove(best_bit)
                current_de17 = best_de17
                greedy_steps += 1
            else:
                break

        all_greedy_best_17.append(current_de17)
        greedy_flips_history.append(greedy_steps)

        # Also track greedy effect on round 18
        if current_flips:
            de_mod_final = track_modular_diff(base_msg, init_state, W_prime, current_flips)
            all_greedy_best_18.append(abs(de_mod_final[1]) if de_mod_final[1] is not None else None)

    # ── Reporting ───────────────────────────────────────────────────

    total_time = time.time() - start_time
    n_msgs = len(all_base_de17)

    print()
    print("=" * 72)
    print("RESULTS")
    print("=" * 72)

    # 1. Distribution of |De_mod[17]| for single neutral flips
    print()
    print("--- Distribution of |De_mod[17]| across single neutral bit flips ---")
    if all_single_de17:
        arr = np.array(all_single_de17, dtype=np.float64)
        print(f"  Total single-bit experiments: {len(arr)}")
        print(f"  Mean |De_mod[17]|:   {arr.mean():.2e}")
        print(f"  Median |De_mod[17]|: {np.median(arr):.2e}")
        print(f"  Min |De_mod[17]|:    {arr.min():.2e}")
        print(f"  Max |De_mod[17]|:    {arr.max():.2e}")
        print(f"  Std |De_mod[17]|:    {arr.std():.2e}")

        # Distribution in powers of 2
        print()
        print("  Distribution by magnitude:")
        for exp in range(0, 33, 4):
            lo = 2**exp if exp > 0 else 0
            hi = 2**(exp + 4)
            count = np.sum((arr >= lo) & (arr < hi))
            pct = 100.0 * count / len(arr)
            if count > 0:
                print(f"    [{lo:>12} , {hi:>12}): {count:>6} ({pct:5.1f}%)")

    # 2. Baseline vs single-flip comparison
    print()
    print("--- Baseline |De_mod[17]| (no neutral flips) ---")
    if all_base_de17:
        base_arr = np.array(all_base_de17, dtype=np.float64)
        print(f"  Mean baseline:  {base_arr.mean():.2e}")
        print(f"  Median baseline: {np.median(base_arr):.2e}")

    # 3. Greedy search results
    print()
    print("--- Greedy search: best |De_mod[17]| ---")
    if all_greedy_best_17:
        garr = np.array(all_greedy_best_17, dtype=np.float64)
        print(f"  Messages tested: {len(garr)}")
        print(f"  Mean greedy best:   {garr.mean():.2e}")
        print(f"  Median greedy best: {np.median(garr):.2e}")
        print(f"  Min greedy best:    {garr.min():.2e}")
        print(f"  Max greedy best:    {garr.max():.2e}")
        print(f"  Mean greedy flips:  {np.mean(greedy_flips_history):.1f}")

        # Compare to random expectation (2^31 ~ 2.15e9)
        random_expect = 2**31
        print()
        print(f"  Random expectation (2^31): {random_expect:.2e}")
        print(f"  Greedy mean / random:      {garr.mean() / random_expect:.4f}")
        print(f"  Greedy min / random:        {garr.min() / random_expect:.6f}")
        beat_random = np.sum(garr < random_expect)
        print(f"  Messages beating random:   {beat_random}/{len(garr)} ({100*beat_random/len(garr):.1f}%)")

    # 4. Greedy effect on round 18
    print()
    print("--- Greedy search: effect on |De_mod[18]| ---")
    valid_18 = [x for x in all_greedy_best_18 if x is not None]
    if valid_18:
        a18 = np.array(valid_18, dtype=np.float64)
        print(f"  Mean |De_mod[18]| after greedy: {a18.mean():.2e}")
        print(f"  Min |De_mod[18]| after greedy:  {a18.min():.2e}")

    # 5. Pairwise interference
    print()
    print("--- Pairwise interference of top neutral bits ---")
    if pair_total > 0:
        print(f"  Total pairs tested: {pair_total}")
        print(f"  Constructive (pair < best individual): {pair_constructive} ({100*pair_constructive/pair_total:.1f}%)")
        print(f"  Destructive (pair > worst individual):  {pair_destructive} ({100*pair_destructive/pair_total:.1f}%)")
        print(f"  Neutral (in between):                   {pair_neutral_count} ({100*pair_neutral_count/pair_total:.1f}%)")

    # 6. Consistently helpful neutral bits
    print()
    print("--- Neutral bits consistently helpful across messages ---")
    print("  (bits that reduce |De_mod[17]| below baseline in many messages)")
    helpful_sorted = sorted(neutral_bit_helpful.items(), key=lambda x: -x[1])
    shown = 0
    for (k, b), count in helpful_sorted[:15]:
        if count >= 2:
            avg_score = np.mean(neutral_bit_scores[(k, b)])
            print(f"    W[{k:>2}] bit {b:>2}: helpful in {count:>3}/{n_msgs} msgs, "
                  f"avg |De_mod[17]| = {avg_score:.2e}")
            shown += 1
    if shown == 0:
        print("    No bits found helpful in 2+ messages.")

    # 7. Top 10 neutral bits by average |De_mod[17]|
    print()
    print("--- Top 10 neutral bits with smallest avg |De_mod[17]| ---")
    bit_avg = {}
    for (k, b), scores in neutral_bit_scores.items():
        if len(scores) >= 2:
            bit_avg[(k, b)] = np.mean(scores)
    top10 = sorted(bit_avg.items(), key=lambda x: x[1])[:10]
    for (k, b), avg in top10:
        n_obs = len(neutral_bit_scores[(k, b)])
        print(f"    W[{k:>2}] bit {b:>2}: avg |De_mod[17]| = {avg:.2e} (n={n_obs})")

    print()
    print("=" * 72)
    print("INTERPRETATION")
    print("=" * 72)

    if all_greedy_best_17 and all_base_de17:
        garr = np.array(all_greedy_best_17, dtype=np.float64)
        base_arr = np.array(all_base_de17, dtype=np.float64)
        ratio = garr.mean() / base_arr.mean() if base_arr.mean() > 0 else float('inf')
        random_ratio = garr.mean() / (2**31)

        print()
        print(f"  Greedy mean / baseline mean: {ratio:.4f}")
        print(f"  Greedy mean / 2^31:          {random_ratio:.4f}")
        print()

        if ratio < 0.5:
            print("  ** STRONG SIGNAL: Greedy neutral-bit selection significantly reduces")
            print("     modular differential at round 17 compared to baseline.")
            print("     The C9×C13 combination yields exploitable structure.")
        elif ratio < 0.9:
            print("  * MODERATE SIGNAL: Some reduction in |De_mod[17]| from greedy search.")
            print("    Neutral bits have differential modular-differential sensitivity.")
        else:
            print("  WEAK/NO SIGNAL: Greedy search does not significantly reduce |De_mod[17]|.")
            print("  Neutral bits appear to affect modular differential roughly uniformly.")

        if pair_total > 0:
            c_ratio = pair_constructive / pair_total
            print()
            if c_ratio > 0.5:
                print(f"  Pairs are mostly constructive ({c_ratio:.1%}): combinatorial gains possible.")
            elif c_ratio > 0.3:
                print(f"  Mixed pair interference ({c_ratio:.1%} constructive): selective pairing needed.")
            else:
                print(f"  Pairs mostly destructive/neutral ({c_ratio:.1%} constructive): limited pair benefit.")

    print()
    print(f"  Total runtime: {total_time:.1f}s")
    print()


if __name__ == "__main__":
    run_experiment(num_base_msgs=50, top_n=30)
