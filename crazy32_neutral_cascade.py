#!/usr/bin/env python3
"""
crazy32_neutral_cascade.py — Neutral bit cascading experiment for SHA-256 collision research.

Tests whether neutral bits (preserving Wang chain through round 16) form a
linear subspace by checking pairs, triples, and larger k-subsets of
simultaneously flipped bits.

Neutral bit: flipping bit (word, pos) in BOTH M and M' preserves
all De values through round 16. ~82 such bits exist per message.

Key question: if we flip MULTIPLE neutral bits at once, do the De values
still remain preserved? If yes → 2^82 free combinations (linear subspace).
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

IV = [
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

def expand_W(W16, n=25):
    W = list(W16)
    for i in range(16, n):
        W.append(add32(sig1(W[i-2]), W[i-7], sig0(W[i-15]), W[i-16]))
    return W

def sha256_states(W_expanded, n_rounds):
    """Run SHA-256 for n_rounds, return list of full 8-tuple states after each round."""
    a, b, c, d, e, f, g, h = IV
    states = []
    for i in range(n_rounds):
        T1 = add32(h, Sig1(e), Ch(e, f, g), K[i], W_expanded[i])
        T2 = add32(Sig0(a), Maj(a, b, c))
        h = g; g = f; f = e
        e = add32(d, T1)
        d = c; c = b; b = a
        a = add32(T1, T2)
        states.append((a, b, c, d, e, f, g, h))
    return states

def compute_wang_chain(base_msg, dw0=0x80000000):
    """Compute twin message M' such that De_t = 0 for rounds 0..15."""
    W = list(base_msg)
    Wp = list(W)
    Wp[0] = W[0] ^ dw0

    a, b, c, d, e, f, g, h = IV
    a2, b2, c2, d2, e2, f2, g2, h2 = IV

    # Round 0
    T1 = add32(h, Sig1(e), Ch(e, f, g), K[0], W[0])
    T2 = add32(Sig0(a), Maj(a, b, c))
    h, g, f = g, f, e
    e = add32(d, T1); d, c, b = c, b, a; a = add32(T1, T2)

    T1p = add32(h2, Sig1(e2), Ch(e2, f2, g2), K[0], Wp[0])
    T2p = add32(Sig0(a2), Maj(a2, b2, c2))
    h2, g2, f2 = g2, f2, e2
    e2 = add32(d2, T1p); d2, c2, b2 = c2, b2, a2; a2 = add32(T1p, T2p)

    for t in range(1, 16):
        base_new_e = add32(d, h, Sig1(e), Ch(e, f, g), K[t], W[t])
        prime_partial = add32(d2, h2, Sig1(e2), Ch(e2, f2, g2), K[t])
        Wp[t] = (base_new_e - prime_partial) & MASK

        T1 = add32(h, Sig1(e), Ch(e, f, g), K[t], W[t])
        T2 = add32(Sig0(a), Maj(a, b, c))
        h, g, f = g, f, e
        e = add32(d, T1); d, c, b = c, b, a; a = add32(T1, T2)

        T1p = add32(h2, Sig1(e2), Ch(e2, f2, g2), K[t], Wp[t])
        T2p = add32(Sig0(a2), Maj(a2, b2, c2))
        h2, g2, f2 = g2, f2, e2
        e2 = add32(d2, T1p); d2, c2, b2 = c2, b2, a2; a2 = add32(T1p, T2p)

    return Wp

def find_neutral_bits(W_base, max_round=13):
    """Find bits in W[1..15] that are neutral for Wang chain through max_round.
    A bit is neutral if flipping it in BOTH M and M' preserves all 8-register
    differentials through max_round.

    Note: max_round=13 gives ~71 neutral bits, max_round=12 gives ~101.
    The '82.6' figure is the average across seeds at an intermediate threshold.
    """
    Wp = compute_wang_chain(W_base)
    We = expand_W(W_base, max_round + 1)
    Wpe = expand_W(Wp, max_round + 1)
    ref_b = sha256_states(We, max_round + 1)
    ref_d = sha256_states(Wpe, max_round + 1)

    ref_diffs = []
    for r in range(max_round + 1):
        diff = tuple(ref_d[r][i] ^ ref_b[r][i] for i in range(8))
        ref_diffs.append(diff)

    neutral = []
    for word_idx in range(1, 16):
        for bit_pos in range(32):
            W2 = list(W_base)
            W2[word_idx] ^= (1 << bit_pos)
            Wp2 = list(Wp)
            Wp2[word_idx] ^= (1 << bit_pos)

            We2 = expand_W(W2, max_round + 1)
            Wpe2 = expand_W(Wp2, max_round + 1)
            st_b2 = sha256_states(We2, max_round + 1)
            st_d2 = sha256_states(Wpe2, max_round + 1)

            is_neutral = True
            for r in range(max_round + 1):
                diff2 = tuple(st_d2[r][i] ^ st_b2[r][i] for i in range(8))
                if diff2 != ref_diffs[r]:
                    is_neutral = False
                    break
            if is_neutral:
                neutral.append((word_idx, bit_pos))

    return neutral


def apply_flips(W_base, Wp_base, bit_list):
    """Flip the given (word, bit) positions in both M and M'."""
    W = list(W_base)
    Wp = list(Wp_base)
    for (word_idx, bit_pos) in bit_list:
        W[word_idx] ^= (1 << bit_pos)
        Wp[word_idx] ^= (1 << bit_pos)
    return W, Wp


def compute_de_rounds(W, Wp, max_round=20):
    """Compute De_t (XOR of e registers) for rounds 0..max_round-1."""
    We = expand_W(W, max_round)
    Wpe = expand_W(Wp, max_round)
    st = sha256_states(We, max_round)
    stp = sha256_states(Wpe, max_round)
    des = []
    for r in range(max_round):
        des.append(st[r][4] ^ stp[r][4])
    return des


def check_neutral_subset(W_base, Wp_base, ref_diffs, subset, max_round=14):
    """Check if flipping all bits in subset preserves differentials through max_round."""
    W, Wp = apply_flips(W_base, Wp_base, subset)
    We = expand_W(W, max_round)
    Wpe = expand_W(Wp, max_round)
    st_b = sha256_states(We, max_round)
    st_d = sha256_states(Wpe, max_round)

    for r in range(max_round):
        diff = tuple(st_d[r][i] ^ st_b[r][i] for i in range(8))
        if diff != ref_diffs[r]:
            return False, r  # fails at round r
    return True, max_round


def compute_De_at_round(W_base, Wp_base, subset, target_round):
    """Compute De (e XOR) at a specific round after flipping subset bits."""
    W, Wp = apply_flips(W_base, Wp_base, subset)
    We = expand_W(W, target_round + 1)
    Wpe = expand_W(Wp, target_round + 1)
    st = sha256_states(We, target_round + 1)
    stp = sha256_states(Wpe, target_round + 1)
    return st[target_round][4] ^ stp[target_round][4]


def main():
    random.seed(0xC032)
    print("=" * 72)
    print("NEUTRAL BIT CASCADING EXPERIMENT")
    print("Seed: 0xC032 | Wang chain with DW[0] = MSB flip")
    print("Neutral bits: preserve full 8-register differential through round 13")
    print("Cascade test: do k simultaneous flips still preserve the differential?")
    print("=" * 72)

    # Step 1: Generate base message and Wang twin
    msg = [random.getrandbits(32) for _ in range(16)]
    Wp = compute_wang_chain(msg)
    print(f"\nBase message W[0..3]: {[hex(w) for w in msg[:4]]}")

    NEUTRAL_ROUND = 13  # neutral bits preserve differential through this round

    # Reference differentials through round 20
    max_check = 21
    We = expand_W(msg, max_check)
    Wpe = expand_W(Wp, max_check)
    ref_st = sha256_states(We, max_check)
    ref_stp = sha256_states(Wpe, max_check)
    ref_diffs = []
    for r in range(max_check):
        diff = tuple(ref_st[r][i] ^ ref_stp[r][i] for i in range(8))
        ref_diffs.append(diff)

    # Show De values
    print("\nReference De values (no neutral bit flips):")
    for r in range(NEUTRAL_ROUND, 21):
        de = ref_st[r][4] ^ ref_stp[r][4]
        print(f"  De{r:2d} = {hex(de)}, hw = {hw(de)}")

    # Step 2: Find neutral bits
    print("\n" + "-" * 72)
    print(f"STEP 2: Finding neutral bits (preserving differential through round {NEUTRAL_ROUND})")
    print("-" * 72)

    t0 = time.time()
    neutral_bits = find_neutral_bits(msg, max_round=NEUTRAL_ROUND)
    elapsed = time.time() - t0
    N = len(neutral_bits)
    print(f"Found {N} neutral bits in W[1..15] (out of 480) in {elapsed:.2f}s")

    by_word = {}
    for (w, b) in neutral_bits:
        by_word.setdefault(w, []).append(b)
    for w in sorted(by_word.keys()):
        bits = sorted(by_word[w])
        print(f"  W[{w:2d}]: {len(bits):2d} neutral bits -> positions {bits}")

    if N < 2:
        print("Too few neutral bits. Aborting.")
        return

    # Step 3: Test ALL PAIRS — do they preserve differential through NEUTRAL_ROUND?
    print("\n" + "-" * 72)
    print(f"STEP 3: Testing ALL pairs of {N} neutral bits ({N*(N-1)//2} pairs)")
    print(f"  Check: differential preserved through round {NEUTRAL_ROUND}")
    print("-" * 72)

    t0 = time.time()
    total_pairs = 0
    compatible_pairs = 0
    interfering_pairs = []

    for i in range(N):
        for j in range(i+1, N):
            ok, fail_r = check_neutral_subset(msg, Wp, ref_diffs,
                                               [neutral_bits[i], neutral_bits[j]],
                                               max_round=NEUTRAL_ROUND + 1)
            total_pairs += 1
            if ok:
                compatible_pairs += 1
            else:
                interfering_pairs.append((neutral_bits[i], neutral_bits[j], fail_r))

    elapsed = time.time() - t0
    frac_pairs = compatible_pairs / total_pairs if total_pairs > 0 else 0
    print(f"Compatible pairs: {compatible_pairs} / {total_pairs} = {frac_pairs*100:.2f}%")
    print(f"Interfering pairs: {len(interfering_pairs)} (elapsed {elapsed:.2f}s)")
    if interfering_pairs:
        print("First 15 interfering pairs (bit_A, bit_B, first_failing_round):")
        for a, b, fr in interfering_pairs[:15]:
            print(f"    W[{a[0]}].{a[1]:2d} + W[{b[0]}].{b[1]:2d}  fails at round {fr}")

    # Step 4: Test RANDOM TRIPLES
    print("\n" + "-" * 72)
    print("STEP 4: Testing random TRIPLES of neutral bits (sample 10000)")
    print("-" * 72)

    t0 = time.time()
    n_triple_samples = min(10000, N*(N-1)*(N-2)//6)
    triple_ok = 0
    for _ in range(n_triple_samples):
        subset = random.sample(neutral_bits, 3)
        ok, _ = check_neutral_subset(msg, Wp, ref_diffs, subset,
                                      max_round=NEUTRAL_ROUND + 1)
        if ok:
            triple_ok += 1
    elapsed = time.time() - t0
    frac_triple = triple_ok / n_triple_samples
    print(f"Triples preserving diff through r{NEUTRAL_ROUND}: {triple_ok} / {n_triple_samples} = {frac_triple*100:.2f}% ({elapsed:.2f}s)")

    # Step 5: Test RANDOM k-SUBSETS
    print("\n" + "-" * 72)
    print(f"STEP 5: Random k-subsets — fraction preserving differential through r{NEUTRAL_ROUND}")
    print("-" * 72)

    all_k_frac = {2: frac_pairs, 3: frac_triple}
    n_samples = 1000
    k_list = [5, 8, 10, 15, 20, 25, 30, 40, 50, 60]

    for k in k_list:
        if k > N:
            print(f"  k={k:3d}: SKIP (only {N} neutral bits)")
            continue
        t0 = time.time()
        ok = 0
        for _ in range(n_samples):
            subset = random.sample(neutral_bits, k)
            passed, _ = check_neutral_subset(msg, Wp, ref_diffs, subset,
                                              max_round=NEUTRAL_ROUND + 1)
            if passed:
                ok += 1
        elapsed = time.time() - t0
        frac = ok / n_samples
        all_k_frac[k] = frac
        print(f"  k={k:3d}: {ok}/{n_samples} = {frac*100:.1f}% preserving  ({elapsed:.2f}s)")

    # Step 6: Cascade effect on De at later rounds (17, 18, 19, 20)
    # Key question: neutral bits preserve diff through r13, but what happens to De17?
    print("\n" + "-" * 72)
    print("STEP 6: Cascade effect on De at later rounds")
    print("Neutral bits preserve diff through r13 — what happens at r17, r18, r19?")
    print("-" * 72)

    for target_r in [14, 15, 16, 17, 18, 19, 20]:
        label = f"De{target_r}"
        de_ref = ref_st[target_r][4] ^ ref_stp[target_r][4]
        print(f"\n  --- {label} (ref hw={hw(de_ref)}) ---")

        # Singles
        single_preserved = 0
        single_hw_diffs = []
        for nb in neutral_bits:
            de_new = compute_De_at_round(msg, Wp, [nb], target_r)
            if de_new == de_ref:
                single_preserved += 1
            else:
                single_hw_diffs.append(hw(de_new ^ de_ref))
        print(f"  Singles: {single_preserved}/{N} = {single_preserved*100/N:.1f}% preserved", end="")
        if single_hw_diffs:
            print(f"  (avg hw_change={sum(single_hw_diffs)/len(single_hw_diffs):.1f})")
        else:
            print()

        # k-subsets — sample for each k
        for k in [2, 5, 10, 20, 40]:
            if k > N:
                continue
            ns = 500
            ok = 0
            hw_diffs = []
            for _ in range(ns):
                subset = random.sample(neutral_bits, k)
                de_new = compute_De_at_round(msg, Wp, subset, target_r)
                if de_new == de_ref:
                    ok += 1
                else:
                    hw_diffs.append(hw(de_new ^ de_ref))
            frac = ok / ns
            extra = ""
            if hw_diffs:
                extra = f"  (avg hw_change={sum(hw_diffs)/len(hw_diffs):.1f})"
            print(f"  k={k:3d}: {ok}/{ns} = {frac*100:.1f}%{extra}")

    # Step 7: Damage analysis for subspace structure
    print("\n" + "-" * 72)
    print("STEP 7: Damage analysis at r13 — hw of diff change when subsets break neutrality")
    print("-" * 72)

    for k in [5, 10, 20, 40]:
        if k > N:
            continue
        hws = []
        for _ in range(1000):
            subset = random.sample(neutral_bits, k)
            ok, fail_r = check_neutral_subset(msg, Wp, ref_diffs, subset,
                                               max_round=NEUTRAL_ROUND + 1)
            if not ok:
                W2, Wp2 = apply_flips(msg, Wp, subset)
                We2 = expand_W(W2, fail_r + 1)
                Wpe2 = expand_W(Wp2, fail_r + 1)
                st2 = sha256_states(We2, fail_r + 1)
                stp2 = sha256_states(Wpe2, fail_r + 1)
                for reg in range(8):
                    diff_new = st2[fail_r][reg] ^ stp2[fail_r][reg]
                    diff_ref = ref_diffs[fail_r][reg]
                    if diff_new != diff_ref:
                        hws.append(hw(diff_new ^ diff_ref))
                        break
        n_broken = len(hws)
        if hws:
            print(f"  k={k:3d}: {n_broken}/1000 broken, avg hw(diff_change) = {sum(hws)/len(hws):.2f}, max = {max(hws)}")
        else:
            print(f"  k={k:3d}: ALL preserved (0/1000 broken)")

    # ── VERDICT ──
    print("\n" + "=" * 72)
    print("VERDICT")
    print("=" * 72)

    # Effective dimension
    eff_dim = 0
    for k in sorted(all_k_frac.keys()):
        if all_k_frac[k] >= 0.5:
            eff_dim = k

    print(f"\nNeutral bits found (through r{NEUTRAL_ROUND}): {N}")
    print(f"Pair compatibility: {frac_pairs*100:.1f}%")
    print(f"Triple compatibility: {frac_triple*100:.1f}%")

    print(f"\nFraction preserving differential through r{NEUTRAL_ROUND} by k:")
    for k in sorted(all_k_frac.keys()):
        bar = "#" * int(all_k_frac[k] * 50)
        print(f"  k={k:3d}: {all_k_frac[k]*100:5.1f}% |{bar}")

    print(f"\nEffective neutral subspace dimension (>50% preserved): k = {eff_dim}")

    if frac_pairs > 0.95:
        print("\nCONCLUSION: Neutral bits are HIGHLY COMPATIBLE in pairs.")
        if frac_triple > 0.90:
            print("Triples also highly compatible — strong evidence for near-linear subspace.")
        elif frac_triple > 0.50:
            print("Triples partially compatible — some higher-order interactions.")
        else:
            print("But triples show major degradation — pairwise linearity does not extend.")
    elif frac_pairs > 0.50:
        print("\nCONCLUSION: Neutral bits are PARTIALLY COMPATIBLE. Interactions reduce effective space.")
    else:
        print("\nCONCLUSION: Neutral bits are STRONGLY INTERACTING. NOT a linear subspace.")

    if eff_dim >= 20:
        print(f"Effective dimension {eff_dim} is LARGE — cascaded flips remain powerful.")
    elif eff_dim >= 5:
        print(f"Effective dimension {eff_dim} is MODERATE — some cascading works.")
    else:
        print(f"Effective dimension {eff_dim} is SMALL — cascading breaks down quickly.")

    print()


if __name__ == "__main__":
    main()
