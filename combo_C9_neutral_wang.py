#!/usr/bin/env python3
"""
combo_C9_neutral_wang.py -- Wang Chain x Neutral Bits

Tests which message bits are NEUTRAL for the Wang chain property in SHA-256.

Wang chain: choose DW[t] to force De_t = 0 for rounds t=1..16.
Neutral bits (De Canniere-Rechberger 2006): bits of the message that can be
flipped without affecting certain round outputs.

The combination question: which bits of W[0..15] are neutral for the Wang
chain property? Flipping them doesn't break De=0 for rounds 1-16.

These neutral bits provide ADDITIONAL freedom for controlling barriers at
rounds 17+, without losing the Wang chain.
"""

import random
import numpy as np
from collections import defaultdict

MASK = 0xFFFFFFFF

# ── SHA-256 constants ─────────────────────────────────────────────────
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
H0 = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
]

# ── SHA-256 primitives ────────────────────────────────────────────────

def rotr(x, n):
    return ((x >> n) | (x << (32 - n))) & MASK

def Sig0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sig1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g):  return ((e & f) ^ (~e & g)) & MASK
def Maj(a, b, c): return ((a & b) ^ (a & c) ^ (b & c)) & MASK
def add32(*args):
    s = 0
    for a in args:
        s = (s + a) & MASK
    return s


def sha256_round(state, W_t, K_t):
    """One SHA-256 compression round. Returns new state [a,b,c,d,e,f,g,h]."""
    a, b, c, d, e, f, g, h = state
    T1 = add32(h, Sig1(e), Ch(e, f, g), K_t, W_t)
    T2 = add32(Sig0(a), Maj(a, b, c))
    new_a = add32(T1, T2)
    new_e = add32(d, T1)
    return [new_a, a, b, c, new_e, e, f, g]


def run_rounds(init_state, W, num_rounds):
    """Run num_rounds of SHA-256, return list of states (state[0]=init, state[t+1]=after round t)."""
    states = [list(init_state)]
    for t in range(num_rounds):
        s = sha256_round(states[-1], W[t], K[t])
        states.append(s)
    return states


# ── Wang Chain: compute DW[1..15] to force De_t = 0 ──────────────────

def compute_wang_dw(base_msg, init_state, dw0=0x80000000):
    """
    Given a base message W[0..15] and initial state, compute DW[0..15]
    such that De_t = 0 for t=1..16.

    DW[0] = dw0 (the input difference, MSB flip).
    DW[t] for t=1..15 is chosen to cancel the difference in T1 at round t,
    so that De_t = d_t + T1_t has zero difference.

    Wang chain condition: De_t = e_t - e'_t = 0 for all rounds 1..16.
    Since new_e = d + T1, and De = Dd + DT1, we need DT1 = -Dd.
    T1 = h + Sig1(e) + Ch(e,f,g) + K + W, so DT1 includes DW contribution.
    If De=0 for all prior rounds, then Df=Dg=Dh=0 (since e shifts to f,g,h).
    So DT1 = DSig1(e) + DCh(e,f,g) + DW when De_prev=0 => DSig1=0, DCh=0.
    Wait - that's only true after enough rounds. Let me think more carefully.

    Actually the Wang chain works by choosing DW[t] to cancel whatever
    difference arises in e at each round. Let me implement it directly.
    """
    W = list(base_msg)
    DW = [0] * 16
    DW[0] = dw0

    # Create the perturbed message
    W_prime = list(W)
    W_prime[0] = (W[0] ^ dw0)  # XOR difference for W[0]

    # Run the base message forward
    state = list(init_state)
    state_prime = list(init_state)

    states = [list(state)]
    states_prime = [list(state_prime)]

    # Round 0: compute with base W[0] and W'[0]
    state = sha256_round(state, W[0], K[0])
    state_prime = sha256_round(state_prime, W_prime[0], K[0])
    states.append(list(state))
    states_prime.append(list(state_prime))

    # For rounds 1..15, choose DW[t] to force e' = e (De = 0)
    for t in range(1, 16):
        # Current state after round t-1
        a, b, c, d, e, f, g, h = state
        a2, b2, c2, d2, e2, f2, g2, h2 = state_prime

        # In round t:
        # new_e  = d  + h  + Sig1(e)  + Ch(e,f,g)  + K[t] + W[t]
        # new_e' = d2 + h2 + Sig1(e2) + Ch(e2,f2,g2) + K[t] + W'[t]
        # We want new_e = new_e', so:
        # W'[t] = W[t] + (d + h + Sig1(e) + Ch(e,f,g)) - (d2 + h2 + Sig1(e2) + Ch(e2,f2,g2))

        T1_partial = add32(h, Sig1(e), Ch(e, f, g), K[t])
        T1_partial_prime = add32(h2, Sig1(e2), Ch(e2, f2, g2), K[t])

        # new_e = d + T1_partial + W[t]
        # new_e' = d2 + T1_partial_prime + W'[t]
        # Set new_e = new_e':
        # d + T1_partial + W[t] = d2 + T1_partial_prime + W'[t]
        # W'[t] = W[t] + (d + T1_partial) - (d2 + T1_partial_prime)

        base_sum = add32(d, T1_partial, W[t])
        prime_need = base_sum  # we want new_e' = new_e = base_sum
        # new_e' = d2 + T1_partial_prime + W'[t] = base_sum
        # W'[t] = base_sum - d2 - T1_partial_prime
        W_prime_t = (base_sum - d2 - T1_partial_prime) & MASK

        DW[t] = (W_prime_t ^ W[t]) if False else (W_prime_t - W[t]) & MASK
        # Store as XOR difference for consistency
        DW[t] = W_prime_t ^ W[t]

        W_prime[t] = W_prime_t

        # Advance both states
        state = sha256_round(state, W[t], K[t])
        state_prime = sha256_round(state_prime, W_prime_t, K[t])
        states.append(list(state))
        states_prime.append(list(state_prime))

    return DW, W_prime, states, states_prime


def verify_wang_chain(states, states_prime, num_rounds=16):
    """Check De_t = 0 for states[2..num_rounds]. Returns list of De values.
    State[1] is intentionally nonzero (initial perturbation)."""
    de_values = []
    for t in range(2, min(num_rounds + 1, len(states))):
        # e is at index 4 in state
        de = (states[t][4] - states_prime[t][4]) & MASK
        de_values.append(de)
    return de_values


# ── Neutral bit detection ─────────────────────────────────────────────

def test_neutral_bits_full(base_msg, init_state):
    """
    For each bit position b in each word W[k] (k=1..15, b=0..31):
    Flip bit b in W[k] in BOTH M and M' (preserving DW[k]).
    Check how many rounds maintain De_t = 0.

    Returns dict: (k, b) -> number of rounds with De=0 (max 15)

    Note: State[1] (after round 0) has De = 0x80000000 by design (the initial
    perturbation). The Wang chain enforces De=0 for states[2..16] (after rounds
    1..15). We check those 15 states.
    """
    DW, W_prime, states_base, states_prime_base = compute_wang_chain(base_msg, init_state)

    results = {}

    for k in range(1, 16):
        for b in range(32):
            flip = 1 << b

            # Flip bit b in W[k] in BOTH messages
            W_flipped = list(base_msg)
            W_flipped[k] ^= flip

            W_prime_flipped = list(W_prime)
            W_prime_flipped[k] ^= flip

            # Recompute all rounds
            s1 = run_rounds(init_state, W_flipped, 16)
            s2 = run_rounds(init_state, W_prime_flipped, 16)

            # Check De_t for states[2..16] (the Wang chain region)
            # State[1] is intentionally nonzero (De = DW[0] propagation)
            rounds_ok = 0
            for t in range(2, 17):
                de = (s1[t][4] - s2[t][4]) & MASK
                if de == 0:
                    rounds_ok += 1
                else:
                    break  # partial neutrality ends here

            results[(k, b)] = rounds_ok

    return results


def compute_wang_chain(base_msg, init_state):
    """Wrapper that returns DW, W_prime, states, states_prime."""
    return compute_wang_dw(base_msg, init_state)


def test_weak_neutral_bits(base_msg, init_state):
    """
    For each bit position b in each word W[k] (k=1..15, b=0..31):
    Flip bit b in W[k] in M ONLY (changing DW[k]).
    Then recompute the Wang chain from that point to see if De=0 can still
    be maintained (by re-choosing subsequent DW values).

    Returns dict: (k, b) -> True if weakly neutral (chain can be readjusted)
    """
    results = {}

    for k in range(1, 16):
        for b in range(32):
            flip = 1 << b

            # Flip bit b in W[k] in base message only
            W_modified = list(base_msg)
            W_modified[k] ^= flip

            # Recompute Wang chain from scratch with modified message
            try:
                DW_new, W_prime_new, states_new, states_prime_new = compute_wang_dw(
                    W_modified, init_state
                )

                # Verify the chain holds
                de_vals = verify_wang_chain(states_new, states_prime_new, 16)
                all_zero = all(d == 0 for d in de_vals) and len(de_vals) == 15
                results[(k, b)] = all_zero
            except Exception:
                results[(k, b)] = False

    return results


# ── Main experiment ───────────────────────────────────────────────────

def run_experiment(num_trials=100):
    random.seed(42)

    # Accumulators
    full_neutral_counts = np.zeros((15, 32), dtype=int)    # [k-1][b], fully neutral
    partial_depth_counts = np.zeros((15, 32), dtype=float)  # average depth
    weak_neutral_counts = np.zeros((15, 32), dtype=int)     # weakly neutral

    # Per-word summaries
    word_full_counts = np.zeros(15, dtype=int)   # total fully-neutral bits across trials
    word_weak_counts = np.zeros(15, dtype=int)

    # Distribution of partial neutrality depths
    depth_distribution = defaultdict(int)  # depth -> count

    print(f"Running {num_trials} trials...")
    print()

    for trial in range(num_trials):
        if trial % 10 == 0:
            print(f"  Trial {trial}/{num_trials}...")

        # Random base message
        base_msg = [random.getrandbits(32) for _ in range(16)]
        init_state = list(H0)

        # Test full neutrality (flip in both M and M')
        full_results = test_neutral_bits_full(base_msg, init_state)

        # Test weak neutrality (flip in M only, readjust chain)
        weak_results = test_weak_neutral_bits(base_msg, init_state)

        for k in range(1, 16):
            for b in range(32):
                depth = full_results[(k, b)]
                partial_depth_counts[k-1][b] += depth
                depth_distribution[depth] += 1

                if depth == 15:  # all 15 Wang chain rounds preserved
                    full_neutral_counts[k-1][b] += 1
                    word_full_counts[k-1] += 1

                if weak_results[(k, b)]:
                    weak_neutral_counts[k-1][b] += 1
                    word_weak_counts[k-1] += 1

    # ── Report results ────────────────────────────────────────────────
    print()
    print("=" * 78)
    print("COMBO C9: Wang Chain x Neutral Bits")
    print("=" * 78)
    print()

    # Table: word W[k] x neutral_count
    print("─" * 78)
    print("Table 1: Neutral bit counts per word (across {} trials)".format(num_trials))
    print("─" * 78)
    print(f"{'Word':>8s} {'Fully Neutral':>14s} {'Avg Full/Trial':>15s} "
          f"{'Weakly Neutral':>15s} {'Avg Weak/Trial':>15s}")
    print(f"{'':>8s} {'(total hits)':>14s} {'(out of 32)':>15s} "
          f"{'(total hits)':>15s} {'(out of 32)':>15s}")
    print("─" * 78)

    total_full = 0
    total_weak = 0
    for k in range(1, 16):
        fc = word_full_counts[k-1]
        wc = word_weak_counts[k-1]
        total_full += fc
        total_weak += wc
        print(f"  W[{k:2d}]  {fc:14d} {fc/num_trials:15.1f} "
              f"{wc:15d} {wc/num_trials:15.1f}")
    print("─" * 78)
    print(f"  TOTAL  {total_full:14d} {total_full/num_trials:15.1f} "
          f"{total_weak:15d} {total_weak/num_trials:15.1f}")
    print()

    # Average partial neutrality depth per word
    print("─" * 78)
    print("Table 2: Average partial neutrality depth per word (max 15)")
    print("─" * 78)
    for k in range(1, 16):
        avg_depth = partial_depth_counts[k-1].sum() / (32 * num_trials)
        print(f"  W[{k:2d}]:  avg depth = {avg_depth:.2f} / 15")
    print()

    # Distribution of partial neutrality depths
    print("─" * 78)
    print("Table 3: Distribution of partial neutrality depths")
    print("─" * 78)
    total_bits = 15 * 32 * num_trials
    print(f"{'Depth':>8s} {'Count':>10s} {'Fraction':>10s}")
    for d in range(16):
        c = depth_distribution[d]
        print(f"{d:8d} {c:10d} {c/total_bits:10.4f}")
    print()

    # Comparison: early vs late words
    print("─" * 78)
    print("Table 4: Early words (W[1..4]) vs Late words (W[12..15])")
    print("─" * 78)
    early_full = sum(word_full_counts[k-1] for k in range(1, 5))
    late_full = sum(word_full_counts[k-1] for k in range(12, 16))
    early_weak = sum(word_weak_counts[k-1] for k in range(1, 5))
    late_weak = sum(word_weak_counts[k-1] for k in range(12, 16))

    early_avg_depth = sum(partial_depth_counts[k-1].sum() for k in range(1, 5)) / (4 * 32 * num_trials)
    late_avg_depth = sum(partial_depth_counts[k-1].sum() for k in range(12, 16)) / (4 * 32 * num_trials)

    print(f"  Early W[1..4]:   fully neutral = {early_full:6d} ({early_full/num_trials:.1f}/trial), "
          f"weakly neutral = {early_weak:6d} ({early_weak/num_trials:.1f}/trial), "
          f"avg depth = {early_avg_depth:.2f}")
    print(f"  Late  W[12..15]: fully neutral = {late_full:6d} ({late_full/num_trials:.1f}/trial), "
          f"weakly neutral = {late_weak:6d} ({late_weak/num_trials:.1f}/trial), "
          f"avg depth = {late_avg_depth:.2f}")
    print()

    # Most consistently neutral bit positions
    print("─" * 78)
    print("Table 5: Top 20 most consistently fully-neutral bit positions")
    print("─" * 78)
    bit_scores = []
    for k in range(1, 16):
        for b in range(32):
            bit_scores.append((full_neutral_counts[k-1][b], k, b))
    bit_scores.sort(reverse=True)

    print(f"{'Rank':>6s} {'Word':>6s} {'Bit':>5s} {'Neutral in N trials':>20s} {'Rate':>8s}")
    for rank, (count, k, b) in enumerate(bit_scores[:20], 1):
        print(f"{rank:6d} W[{k:2d}] {b:5d} {count:20d} {count/num_trials:8.2%}")
    print()

    # ── Verdict ───────────────────────────────────────────────────────
    # Count average fully neutral bits per trial
    avg_full_per_trial = total_full / num_trials
    avg_weak_per_trial = total_weak / num_trials

    print("=" * 78)
    print("VERDICT")
    print("=" * 78)
    print(f"  Average fully neutral bits per trial:   {avg_full_per_trial:.1f}")
    print(f"  Average weakly neutral bits per trial:  {avg_weak_per_trial:.1f}")
    print()

    if avg_full_per_trial >= 10:
        verdict = "ALIVE"
        reason = (f"Found {avg_full_per_trial:.1f} fully neutral bits on average, "
                  "providing substantial extra freedom for controlling barriers at rounds 17+.")
    elif avg_full_per_trial >= 1:
        verdict = "ANOMALY"
        reason = (f"Found {avg_full_per_trial:.1f} fully neutral bits on average. "
                  "Some freedom exists but limited.")
    else:
        verdict = "DEAD"
        reason = "No fully neutral bits found. Wang chain is rigid with no extra freedom."

    # Also consider weak neutrality for context
    if verdict == "DEAD" and avg_weak_per_trial >= 10:
        reason += (f"\n  However, {avg_weak_per_trial:.1f} weakly neutral bits exist "
                   "(chain can be readjusted).")

    print(f"  Verdict: ** {verdict} **")
    print(f"  Reason:  {reason}")
    print()
    print("  Criteria:")
    print("    ALIVE   if >= 10 fully neutral bits (extra freedom for barriers)")
    print("    ANOMALY if 1-9 fully neutral bits")
    print("    DEAD    if 0 fully neutral bits")
    print("=" * 78)


if __name__ == "__main__":
    run_experiment(num_trials=100)
