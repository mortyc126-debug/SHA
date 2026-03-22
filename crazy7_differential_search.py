#!/usr/bin/env python3
"""
CRAZY-7: Starting Differential Optimization
============================================
We've been using ONE fixed differential: ΔW[0] = 0x80000000 (MSB flip).
But SHA-256 has 2^32 possible single-word differentials and 2^512 multi-word.

Key questions:
1. Does the CHOICE of starting differential affect neutral bit count?
2. Does it affect barrier cost (HW(De) at rounds 17+)?
3. Can we find a differential with dramatically more neutral bits or lower barriers?
4. What about multi-word differentials?

This is a CRITICAL unexplored direction — if the optimal differential gives
200 neutral bits instead of 82, or barriers with HW=3 instead of 7,
the attack landscape changes fundamentally.
"""

import random
import time

MASK = 0xFFFFFFFF

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

def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK
def Sig0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sig1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def sig0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sig1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)
def Ch(e, f, g): return ((e & f) ^ (~e & g)) & MASK
def Maj(a, b, c): return ((a & b) ^ (a & c) ^ (b & c)) & MASK
def add32(*args):
    s = 0
    for a in args: s = (s + a) & MASK
    return s
def hw(x): return bin(x & MASK).count('1')

def sha256_round(state, W_t, K_t):
    a, b, c, d, e, f, g, h = state
    T1 = add32(h, Sig1(e), Ch(e, f, g), K_t, W_t)
    T2 = add32(Sig0(a), Maj(a, b, c))
    return [add32(T1, T2), a, b, c, add32(d, T1), e, f, g]

def expand_W(W16):
    W = list(W16)
    for i in range(16, 64):
        W.append(add32(sig1(W[i-2]), W[i-7], sig0(W[i-15]), W[i-16]))
    return W

def run_rounds(init_state, W, num_rounds):
    states = [list(init_state)]
    for t in range(num_rounds):
        s = sha256_round(states[-1], W[t], K[t])
        states.append(s)
    return states

def compute_wang_chain_general(base_msg, init_state, diff_word, diff_val):
    """Wang chain with arbitrary differential: W'[diff_word] = W[diff_word] ^ diff_val.
    Adjusts W'[t] for t != diff_word to keep De[1..16] = 0."""
    W = list(base_msg)
    W_prime = list(W)
    W_prime[diff_word] = W[diff_word] ^ diff_val

    state = list(init_state)
    state_prime = list(init_state)

    # Process rounds 0..15
    for t in range(16):
        if t == diff_word:
            # This round has the actual difference
            state = sha256_round(state, W[t], K[t])
            state_prime = sha256_round(state_prime, W_prime[t], K[t])
        elif t < diff_word:
            # Before difference: W[t] = W'[t], same processing
            state = sha256_round(state, W[t], K[t])
            state_prime = sha256_round(state_prime, W[t], K[t])
        else:
            # After difference: adjust W'[t] to cancel De[t]
            # We want the 'e' register to match: e_{t+1} = e'_{t+1}
            # e_{t+1} = d_t + T1_t = d_t + h_t + Sig1(e_t) + Ch(e_t,f_t,g_t) + K_t + W_t
            a, b, c, d, e, f, g, h = state
            a2, b2, c2, d2, e2, f2, g2, h2 = state_prime

            T1_partial = add32(h, Sig1(e), Ch(e, f, g), K[t])
            T1_partial_prime = add32(h2, Sig1(e2), Ch(e2, f2, g2), K[t])

            target_e = add32(d, T1_partial, W[t])  # desired e_{t+1}
            # e'_{t+1} = d2 + T1_partial_prime + W'[t]
            # We want e'_{t+1} = target_e
            W_prime[t] = (target_e - d2 - T1_partial_prime) & MASK

            state = sha256_round(state, W[t], K[t])
            state_prime = sha256_round(state_prime, W_prime[t], K[t])

    return W_prime


def count_neutral_bits(base_msg, init_state, diff_word, diff_val, max_check=20):
    """Count neutral bits for a given starting differential.
    A bit (k, b) is neutral if flipping W[k] bit b AND W'[k] bit b
    preserves the Wang chain property De[diff_word+2..16] = 0.

    This matches the original semantics: flip in BOTH messages, don't recompute."""
    W_prime = compute_wang_chain_general(base_msg, init_state, diff_word, diff_val)

    neutrals = []
    for k in range(16):
        if k == diff_word:
            continue
        for b in range(32):
            flip = 1 << b
            W_f = list(base_msg)
            W_f[k] ^= flip
            Wp_f = list(W_prime)
            Wp_f[k] ^= flip  # Same flip in BOTH messages

            W_exp = expand_W(W_f)
            Wp_exp = expand_W(Wp_f)
            s1 = run_rounds(init_state, W_exp, 16)
            s2 = run_rounds(init_state, Wp_exp, 16)

            ok = True
            for t in range(diff_word + 2, 17):
                if (s1[t][4] - s2[t][4]) & MASK != 0:
                    ok = False
                    break
            if ok:
                neutrals.append((k, b))
    return neutrals


def measure_barriers(base_msg, init_state, diff_word, diff_val, max_round=21):
    """Measure De at rounds 17-20 for a given differential."""
    W = list(base_msg)
    W_prime = compute_wang_chain_general(W, init_state, diff_word, diff_val)

    W_exp = expand_W(W)
    Wp_exp = expand_W(W_prime)

    s1 = run_rounds(init_state, W_exp, max_round)
    s2 = run_rounds(init_state, Wp_exp, max_round)

    de_hw = {}
    state_hw = {}
    for t in range(diff_word + 1, min(max_round + 1, len(s1))):
        de = (s1[t][4] ^ s2[t][4]) & MASK
        de_hw[t] = hw(de)
        total = sum(hw((s1[t][r] ^ s2[t][r]) & MASK) for r in range(8))
        state_hw[t] = total

    return de_hw, state_hw


def main():
    print("=" * 72)
    print("CRAZY-7: Starting Differential Optimization")
    print("=" * 72)
    print()

    random.seed(0xDEAD7)
    t_start = time.time()
    init_state = list(H0)

    # ── Phase 1: Single-bit differentials in W[0] ─────────────────
    print("Phase 1: All 32 single-bit differentials in W[0]")
    print("-" * 72)

    NUM_MSGS = 5  # average over 5 messages for stability
    base_msgs = [[random.getrandbits(32) for _ in range(16)] for _ in range(NUM_MSGS)]

    results_phase1 = []
    for bit in range(32):
        diff_val = 1 << bit
        total_neutrals = 0
        total_de17 = 0
        total_de18 = 0
        count = 0

        for base_msg in base_msgs:
            neutrals = count_neutral_bits(base_msg, init_state, 0, diff_val)
            de_hw, _ = measure_barriers(base_msg, init_state, 0, diff_val)
            total_neutrals += len(neutrals)
            total_de17 += de_hw.get(17, 16)
            total_de18 += de_hw.get(18, 16)
            count += 1

        avg_n = total_neutrals / count
        avg_de17 = total_de17 / count
        avg_de18 = total_de18 / count
        results_phase1.append((bit, avg_n, avg_de17, avg_de18))

        elapsed = time.time() - t_start
        if bit % 8 == 7:
            print(f"  Bit {bit:2d}: neutrals={avg_n:.1f}, HW(De17)={avg_de17:.1f}, "
                  f"HW(De18)={avg_de18:.1f}  [{elapsed:.1f}s]")

    print()
    print("Summary — Single-bit differentials in W[0]:")
    print(f"{'Bit':>4} {'Neutrals':>9} {'HW(De17)':>9} {'HW(De18)':>9} {'Score':>7}")
    print("-" * 42)
    # Score = neutrals - 2*(de17 + de18)
    results_phase1.sort(key=lambda x: -x[1] + 2*(x[2] + x[3]))
    for bit, n, d17, d18 in results_phase1:
        score = n - 2*(d17 + d18)
        print(f"  {bit:2d}  {n:9.1f} {d17:9.1f} {d18:9.1f} {score:7.1f}")

    best_bit = max(results_phase1, key=lambda x: x[1])
    print(f"\nBest neutral count: bit {best_bit[0]} with {best_bit[1]:.1f} neutrals")
    baseline_bit = [r for r in results_phase1 if r[0] == 31][0]
    print(f"Baseline (bit 31/MSB): {baseline_bit[1]:.1f} neutrals, "
          f"HW(De17)={baseline_bit[2]:.1f}")

    if time.time() - t_start > 300:
        print("[Time limit reached, skipping further phases]")
        return

    # ── Phase 2: Differential in different words ──────────────────
    print()
    print("=" * 72)
    print("Phase 2: MSB differential in different words W[0]..W[15]")
    print("-" * 72)

    diff_val = 0x80000000  # MSB flip
    results_phase2 = []

    for word in range(16):
        if time.time() - t_start > 300:
            print(f"  [Time limit at word {word}]")
            break

        total_neutrals = 0
        total_de = {t: 0 for t in range(17, 21)}
        count = 0

        for base_msg in base_msgs[:3]:  # fewer messages for speed
            neutrals = count_neutral_bits(base_msg, init_state, word, diff_val)
            de_hw, _ = measure_barriers(base_msg, init_state, word, diff_val, max_round=21)
            total_neutrals += len(neutrals)
            for t in range(17, 21):
                total_de[t] += de_hw.get(t, 16)
            count += 1

        avg_n = total_neutrals / count
        avg_de = {t: total_de[t] / count for t in range(17, 21)}
        results_phase2.append((word, avg_n, avg_de))

        elapsed = time.time() - t_start
        print(f"  W[{word:2d}]: neutrals={avg_n:.1f}, De17={avg_de[17]:.1f}, "
              f"De18={avg_de[18]:.1f}, De19={avg_de[19]:.1f}, De20={avg_de[20]:.1f}  [{elapsed:.1f}s]")

    print()
    print("Key insight: which word gives best neutral bits?")
    if results_phase2:
        best_word = max(results_phase2, key=lambda x: x[1])
        print(f"  Best: W[{best_word[0]}] with {best_word[1]:.1f} neutrals")
        print(f"  Worst: W[{min(results_phase2, key=lambda x: x[1])[0]}] with "
              f"{min(results_phase2, key=lambda x: x[1])[1]:.1f} neutrals")

    if time.time() - t_start > 300:
        print("[Time limit reached, skipping Phase 3]")
        return

    # ── Phase 3: Multi-word differentials ─────────────────────────
    print()
    print("=" * 72)
    print("Phase 3: Two-word differentials (W[0] + W[k]) with MSB flips")
    print("-" * 72)

    # For 2-word diff, we modify Wang chain to handle two diff words
    # This is more complex — let's test a few combinations

    results_phase3 = []
    for word2 in range(1, 8):
        if time.time() - t_start > 400:
            break

        # Differential: W[0] ^= MSB, W[word2] ^= MSB
        # Wang chain must handle both
        total_neutrals = 0
        total_de17 = 0
        count = 0

        for base_msg in base_msgs[:2]:
            # Create two-word differential manually
            W = list(base_msg)
            W_prime = list(W)
            W_prime[0] ^= 0x80000000

            # Run rounds 0..word2-1 with Wang chain (adjusting for W[0] diff)
            state = list(init_state)
            state_prime = list(init_state)

            for t in range(16):
                if t == 0:
                    state = sha256_round(state, W[t], K[t])
                    state_prime = sha256_round(state_prime, W_prime[t], K[t])
                elif t == word2:
                    # Add second differential here
                    W_prime[t] = W[t] ^ 0x80000000
                    # But also adjust for Wang chain
                    a, b, c, d, e, f, g, h = state
                    a2, b2, c2, d2, e2, f2, g2, h2 = state_prime
                    # Don't adjust — let the second diff propagate
                    state = sha256_round(state, W[t], K[t])
                    state_prime = sha256_round(state_prime, W_prime[t], K[t])
                else:
                    # Standard Wang chain adjustment
                    a, b, c, d, e, f, g, h = state
                    a2, b2, c2, d2, e2, f2, g2, h2 = state_prime
                    T1_p = add32(h, Sig1(e), Ch(e, f, g), K[t])
                    T1_p2 = add32(h2, Sig1(e2), Ch(e2, f2, g2), K[t])
                    target_e = add32(d, T1_p, W[t])
                    W_prime[t] = (target_e - d2 - T1_p2) & MASK
                    state = sha256_round(state, W[t], K[t])
                    state_prime = sha256_round(state_prime, W_prime[t], K[t])

            # Check De at round 17
            W_exp = expand_W(W)
            Wp_exp = expand_W(W_prime)
            s1 = run_rounds(init_state, W_exp, 18)
            s2 = run_rounds(init_state, Wp_exp, 18)

            de17 = hw((s1[17][4] ^ s2[17][4]) & MASK)
            total_de17 += de17

            # Count neutrals (simplified — just check a subset)
            n_count = 0
            for k in range(16):
                if k == 0 or k == word2:
                    continue
                for b in [0, 7, 15, 23, 31]:  # sample 5 bits per word
                    flip = 1 << b
                    W_f = list(W)
                    W_f[k] ^= flip
                    W_exp_f = expand_W(W_f)
                    s_f = run_rounds(init_state, W_exp_f, 17)

                    Wp_f = list(W_prime)
                    Wp_f[k] ^= flip
                    Wp_exp_f = expand_W(Wp_f)
                    sp_f = run_rounds(init_state, Wp_exp_f, 17)

                    ok = True
                    for t in range(2, 17):
                        if (s_f[t][4] ^ sp_f[t][4]) & MASK != (s1[t][4] ^ s2[t][4]) & MASK:
                            ok = False
                            break
                    if ok:
                        n_count += 1

            total_neutrals += n_count
            count += 1

        if count > 0:
            avg_n = total_neutrals / count * (32/5)  # extrapolate from 5 bits sample
            avg_de17 = total_de17 / count
            results_phase3.append((word2, avg_n, avg_de17))
            print(f"  W[0]+W[{word2}]: est. neutrals~{avg_n:.0f}, HW(De17)={avg_de17:.1f}")

    # ── Phase 4: Low-HW differentials ─────────────────────────────
    if time.time() - t_start < 400:
        print()
        print("=" * 72)
        print("Phase 4: Low Hamming weight differentials in W[0]")
        print("-" * 72)

        # Test 2-bit differentials (C(32,2) = 496, but sample)
        results_phase4 = []
        tested = 0
        for b1 in range(0, 32, 4):
            for b2 in range(b1+1, 32, 4):
                if time.time() - t_start > 450:
                    break
                diff_val = (1 << b1) | (1 << b2)
                total_n = 0
                total_de17 = 0
                count = 0
                for base_msg in base_msgs[:2]:
                    neutrals = count_neutral_bits(base_msg, init_state, 0, diff_val)
                    de_hw, _ = measure_barriers(base_msg, init_state, 0, diff_val)
                    total_n += len(neutrals)
                    total_de17 += de_hw.get(17, 16)
                    count += 1
                avg_n = total_n / count
                avg_de17 = total_de17 / count
                results_phase4.append((b1, b2, avg_n, avg_de17))
                tested += 1

        if results_phase4:
            best = max(results_phase4, key=lambda x: x[2])
            worst = min(results_phase4, key=lambda x: x[2])
            print(f"  Tested {tested} 2-bit differentials")
            print(f"  Best:  bits ({best[0]},{best[1]}), neutrals={best[2]:.1f}, De17={best[3]:.1f}")
            print(f"  Worst: bits ({worst[0]},{worst[1]}), neutrals={worst[2]:.1f}, De17={worst[3]:.1f}")
            avg_all = sum(r[2] for r in results_phase4) / len(results_phase4)
            print(f"  Average neutrals: {avg_all:.1f}")
            print(f"  vs 1-bit baseline (bit 31): {baseline_bit[1]:.1f}")

    # ── Final Summary ─────────────────────────────────────────────
    print()
    print("=" * 72)
    print("VERDICT")
    print("=" * 72)
    print()

    # Compare best single-bit with baseline
    best_single = max(results_phase1, key=lambda x: x[1])
    print(f"  Best single-bit differential: bit {best_single[0]}")
    print(f"    Neutrals: {best_single[1]:.1f} (vs baseline bit 31: {baseline_bit[1]:.1f})")
    print(f"    Improvement: {best_single[1] - baseline_bit[1]:.1f} bits")
    print()

    if results_phase2:
        best_word_r = max(results_phase2, key=lambda x: x[1])
        worst_word_r = min(results_phase2, key=lambda x: x[1])
        print(f"  Best differential word: W[{best_word_r[0]}] ({best_word_r[1]:.1f} neutrals)")
        print(f"  Worst differential word: W[{worst_word_r[0]}] ({worst_word_r[1]:.1f} neutrals)")
        spread = best_word_r[1] - worst_word_r[1]
        print(f"  Spread: {spread:.1f} bits")
        print()

    # Is the choice of differential a game-changer?
    if best_single[1] > baseline_bit[1] * 1.5:
        print("  **ALIVE** — Differential choice MATTERS significantly!")
        print(f"  Optimal differential gives {best_single[1]/baseline_bit[1]:.1f}x more neutral bits")
    elif best_single[1] > baseline_bit[1] * 1.2:
        print("  **ANOMALY** — Some improvement from differential choice")
        print(f"  But only {(best_single[1]/baseline_bit[1]-1)*100:.0f}% more neutral bits")
    else:
        print("  **DEAD** — Differential choice does NOT matter much")
        print(f"  All single-bit differentials give similar neutral bit counts")
        print(f"  Range: {min(r[1] for r in results_phase1):.1f} — {max(r[1] for r in results_phase1):.1f}")

    print()
    print(f"  Runtime: {time.time() - t_start:.1f}s")
    print("=" * 72)


if __name__ == "__main__":
    main()
