#!/usr/bin/env python3
"""
combo_C9_C12_neutral_cube.py
LAYER 2 synergy test: C9 (Neutral Bits) × C12 (Cube Attack)

Key idea: C12 alone tested cube attack on De17 using only W[14] bits,
finding degree >= 7 (DEAD). But C9 found ~82 neutral bits spread across
W[1..15]. Using these NEUTRAL bits as cube variables might reveal
lower-degree algebraic structure in De17 because:
  - Different neutral bits enter at different rounds
  - Multi-word variation may create cancellations
  - Cross-word neutral cubes might yield simpler superpolynomials

A bit is NEUTRAL if flipping it in both M and M' preserves ALL
differential values De_r for rounds 0..16. This means the bit doesn't
disturb the differential characteristic at all through round 16.
"""

import struct
import random
import time
from collections import defaultdict

# ── SHA-256 constants ──────────────────────────────────────────────
K = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786,
]

IV = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
]

def rr(x, n):
    return ((x >> n) | (x << (32 - n))) & 0xFFFFFFFF

def shr(x, n):
    return x >> n

def add32(*args):
    s = 0
    for a in args:
        s = (s + a) & 0xFFFFFFFF
    return s

def Ch(e, f, g):
    return (e & f) ^ ((~e) & g) & 0xFFFFFFFF

def Maj(a, b, c):
    return (a & b) ^ (a & c) ^ (b & c)

def Sigma0(a):
    return rr(a, 2) ^ rr(a, 13) ^ rr(a, 22)

def Sigma1(e):
    return rr(e, 6) ^ rr(e, 11) ^ rr(e, 25)

def sigma0(x):
    return rr(x, 7) ^ rr(x, 18) ^ shr(x, 3)

def sigma1(x):
    return rr(x, 17) ^ rr(x, 19) ^ shr(x, 10)

def expand_W(W16):
    """Expand 16-word message to 18 words."""
    W = list(W16)
    for i in range(16, 18):
        W.append(add32(sigma1(W[i-2]), W[i-7], sigma0(W[i-15]), W[i-16]))
    return W

def sha256_rounds(W, n_rounds):
    """Run SHA-256 for n_rounds, return state [a,b,c,d,e,f,g,h]."""
    a, b, c, d, e, f, g, h = IV
    for i in range(n_rounds):
        T1 = add32(h, Sigma1(e), Ch(e, f, g), K[i], W[i])
        T2 = add32(Sigma0(a), Maj(a, b, c))
        h = g
        g = f
        f = e
        e = add32(d, T1)
        d = c
        c = b
        b = a
        a = add32(T1, T2)
    return [a, b, c, d, e, f, g, h]

def sha256_all_states(W, n_rounds):
    """Run SHA-256 for n_rounds, return list of (e, a) after each round."""
    a, b, c, d, e, f, g, h = IV
    states = []
    for i in range(n_rounds):
        T1 = add32(h, Sigma1(e), Ch(e, f, g), K[i], W[i])
        T2 = add32(Sigma0(a), Maj(a, b, c))
        h = g
        g = f
        f = e
        e = add32(d, T1)
        d = c
        c = b
        b = a
        a = add32(T1, T2)
        states.append((a, b, c, d, e, f, g, h))
    return states

def get_e(state):
    return state[4]

def flip_bit(W, word_idx, bit_pos):
    W2 = list(W)
    W2[word_idx] ^= (1 << bit_pos)
    return W2

def find_neutral_bits(W_base, W_diff, max_round=16):
    """
    Find bits in W[1..15] that are neutral for the differential characteristic
    through round max_round.

    A bit (word, pos) is NEUTRAL if flipping it in BOTH W_base and W_diff
    does NOT change any of the state differences (Da, Db, ..., Dh) for
    rounds 0..max_round.

    We check: for each round r in [0..max_round], the full 8-register
    state difference must be identical before and after the bit flip.
    """
    # Compute reference differential states
    W_b_exp = expand_W(W_base)
    W_d_exp = expand_W(W_diff)
    ref_states_b = sha256_all_states(W_b_exp, max_round + 1)
    ref_states_d = sha256_all_states(W_d_exp, max_round + 1)

    # Reference differentials for each round
    ref_diffs = []
    for r in range(max_round + 1):
        diff = tuple(ref_states_d[r][i] ^ ref_states_b[r][i] for i in range(8))
        ref_diffs.append(diff)

    neutral = []
    for word_idx in range(1, 16):
        for bit_pos in range(32):
            W_b2 = flip_bit(W_base, word_idx, bit_pos)
            W_d2 = flip_bit(W_diff, word_idx, bit_pos)
            W_b2_exp = expand_W(W_b2)
            W_d2_exp = expand_W(W_d2)

            states_b2 = sha256_all_states(W_b2_exp, max_round + 1)
            states_d2 = sha256_all_states(W_d2_exp, max_round + 1)

            is_neutral = True
            for r in range(max_round + 1):
                diff2 = tuple(states_d2[r][i] ^ states_b2[r][i] for i in range(8))
                if diff2 != ref_diffs[r]:
                    is_neutral = False
                    break

            if is_neutral:
                neutral.append((word_idx, bit_pos))

    return neutral

def compute_De17_from_msgs(W_base, W_diff):
    """Compute De at round 17 (after 18 rounds, 0-indexed)."""
    W_b = expand_W(W_base)
    W_d = expand_W(W_diff)
    st_b = sha256_rounds(W_b, 18)
    st_d = sha256_rounds(W_d, 18)
    return get_e(st_d) ^ get_e(st_b)


def superpoly_eval(W_base, W_diff, cube_positions, free_bits, output_bit, n_eval=50):
    """
    Evaluate the superpoly of De17[output_bit] summed over cube_positions,
    at n_eval random settings of the free_bits.

    cube_positions: list of (word, bit) to sum over (cube variables)
    free_bits: list of (word, bit) to randomize (superpoly variables)
    output_bit: which bit of De17 to extract

    Returns list of superpoly values (0 or 1).
    """
    results = []
    n_cube = len(cube_positions)

    for _ in range(n_eval):
        # Randomize free bits
        W_rb = list(W_base)
        W_rd = list(W_diff)
        for (w, b) in free_bits:
            if random.randint(0, 1):
                W_rb = flip_bit(W_rb, w, b)
                W_rd = flip_bit(W_rd, w, b)

        # Sum over cube
        cube_sum = 0
        for mask in range(1 << n_cube):
            W_cb = list(W_rb)
            W_cd = list(W_rd)
            for idx in range(n_cube):
                if mask & (1 << idx):
                    w, bpos = cube_positions[idx]
                    W_cb = flip_bit(W_cb, w, bpos)
                    W_cd = flip_bit(W_cd, w, bpos)

            De17 = compute_De17_from_msgs(W_cb, W_cd)
            cube_sum ^= (De17 >> output_bit) & 1

        results.append(cube_sum)

    return results


def check_constant(values):
    return len(set(values)) == 1


def blr_linearity_test(W_base, W_diff, cube_positions, free_bits,
                       output_bit, n_tests=15):
    """
    BLR linearity test on the superpoly as a function of free_bits.
    f(x) + f(y) == f(x XOR y) + f(0) for random x, y.
    """
    n_free = len(free_bits)
    if n_free == 0:
        return 1.0

    n_cube = len(cube_positions)

    def eval_sp(settings):
        W_rb = list(W_base)
        W_rd = list(W_diff)
        for i, (w, b) in enumerate(free_bits):
            if settings[i]:
                W_rb = flip_bit(W_rb, w, b)
                W_rd = flip_bit(W_rd, w, b)

        cube_sum = 0
        for mask in range(1 << n_cube):
            W_cb = list(W_rb)
            W_cd = list(W_rd)
            for idx in range(n_cube):
                if mask & (1 << idx):
                    w, bpos = cube_positions[idx]
                    W_cb = flip_bit(W_cb, w, bpos)
                    W_cd = flip_bit(W_cd, w, bpos)
            De17 = compute_De17_from_msgs(W_cb, W_cd)
            cube_sum ^= (De17 >> output_bit) & 1
        return cube_sum

    passes = 0
    zero = [0] * n_free
    f0 = eval_sp(zero)

    for _ in range(n_tests):
        x = [random.randint(0, 1) for _ in range(n_free)]
        y = [random.randint(0, 1) for _ in range(n_free)]
        xy = [a ^ b for a, b in zip(x, y)]

        fx = eval_sp(x)
        fy = eval_sp(y)
        fxy = eval_sp(xy)

        if (fx ^ fy) == (fxy ^ f0):
            passes += 1

    return passes / n_tests


def select_cross_word(neutral_bits, k):
    """Select k neutral bits from DIFFERENT words."""
    by_word = defaultdict(list)
    for (w, b) in neutral_bits:
        by_word[w].append((w, b))
    words = list(by_word.keys())
    if len(words) < k:
        return None
    chosen_words = random.sample(words, k)
    return [random.choice(by_word[w]) for w in chosen_words]


def select_same_word(neutral_bits, k):
    """Select k neutral bits from the SAME word."""
    by_word = defaultdict(list)
    for (w, b) in neutral_bits:
        by_word[w].append((w, b))
    eligible = [(w, bits) for w, bits in by_word.items() if len(bits) >= k]
    if not eligible:
        return None
    w, bits = random.choice(eligible)
    return random.sample(bits, k)


def select_non_neutral(neutral_set, k):
    """Select k NON-neutral bit positions from W[1..15]."""
    non_neutral = []
    for w in range(1, 16):
        for b in range(32):
            if (w, b) not in neutral_set:
                non_neutral.append((w, b))
    if len(non_neutral) < k:
        return None
    return random.sample(non_neutral, k)


def main():
    random.seed(42)
    t0 = time.time()

    print("=" * 72)
    print("LAYER 2: C9+C12 Neutral Bits x Cube Attack Synergy Test")
    print("=" * 72)

    # ── Generate random base message ──
    W_base = [random.getrandbits(32) for _ in range(16)]
    W_diff = list(W_base)
    W_diff[0] ^= 0x80000000  # DW[0] = MSB flip (Wang differential)

    print("\n[1] Finding neutral bits (preserve full state diff, rounds 0-16)...")
    t1 = time.time()
    neutral_bits = find_neutral_bits(W_base, W_diff, max_round=16)
    t_neutral = time.time() - t1
    print(f"    Found {len(neutral_bits)} neutral bits in {t_neutral:.1f}s")

    by_word = defaultdict(int)
    for (w, b) in neutral_bits:
        by_word[w] += 1
    print("    Distribution by word:")
    for w in sorted(by_word.keys()):
        print(f"      W[{w:2d}]: {by_word[w]:2d} neutral bits")

    n_words_with_neutral = len(by_word)
    print(f"    Words with neutral bits: {n_words_with_neutral}")

    neutral_set = set(neutral_bits)

    # If very few neutral bits found, also try a relaxed definition:
    # only check De (e-register diff), not full state diff
    if len(neutral_bits) < 20:
        print("\n    Few full-state neutral bits found. Trying e-register-only neutral...")
        neutral_bits_e = []
        W_b_exp = expand_W(W_base)
        W_d_exp = expand_W(W_diff)
        ref_De = []
        for r in range(17):
            st_b = sha256_rounds(W_b_exp, r + 1)
            st_d = sha256_rounds(W_d_exp, r + 1)
            ref_De.append(get_e(st_d) ^ get_e(st_b))

        for word_idx in range(1, 16):
            for bit_pos in range(32):
                W_b2 = flip_bit(W_base, word_idx, bit_pos)
                W_d2 = flip_bit(W_diff, word_idx, bit_pos)
                W_b2_exp = expand_W(W_b2)
                W_d2_exp = expand_W(W_d2)

                ok = True
                for r in range(17):
                    st_b2 = sha256_rounds(W_b2_exp, r + 1)
                    st_d2 = sha256_rounds(W_d2_exp, r + 1)
                    De2 = get_e(st_d2) ^ get_e(st_b2)
                    if De2 != ref_De[r]:
                        ok = False
                        break
                if ok:
                    neutral_bits_e.append((word_idx, bit_pos))

        print(f"    Found {len(neutral_bits_e)} e-register neutral bits")
        by_word_e = defaultdict(int)
        for (w, b) in neutral_bits_e:
            by_word_e[w] += 1
        for w in sorted(by_word_e.keys()):
            print(f"      W[{w:2d}]: {by_word_e[w]:2d} e-neutral bits")

        if len(neutral_bits_e) > len(neutral_bits):
            print("    Using e-register neutral bits for cube experiment.")
            neutral_bits = neutral_bits_e
            neutral_set = set(neutral_bits)
            by_word = by_word_e
            n_words_with_neutral = len(by_word)

    # If still very few, try progressively shorter round ranges
    if len(neutral_bits) < 20:
        for max_r in [12, 10, 8, 6]:
            print(f"\n    Still few bits. Trying max_round={max_r}...")
            nb = find_neutral_bits(W_base, W_diff, max_round=max_r)
            print(f"    Found {len(nb)} neutral bits (rounds 0-{max_r})")
            bw = defaultdict(int)
            for (w, b) in nb:
                bw[w] += 1
            for w in sorted(bw.keys()):
                print(f"      W[{w:2d}]: {bw[w]:2d} bits")
            if len(nb) >= 20:
                print(f"    Using round-{max_r} neutral bits for cube experiment.")
                neutral_bits = nb
                neutral_set = set(nb)
                by_word = bw
                n_words_with_neutral = len(bw)
                break

    # Verify De17 is nonzero
    De17_base = compute_De17_from_msgs(W_base, W_diff)
    print(f"\n    De17 (base msg) = 0x{De17_base:08x}")
    print(f"    Total neutral bits for experiment: {len(neutral_bits)}")

    if len(neutral_bits) < 2:
        print("\n    ERROR: Not enough neutral bits found. Trying multiple random messages...")
        best_nb = neutral_bits
        for attempt in range(20):
            W_b2 = [random.getrandbits(32) for _ in range(16)]
            W_d2 = list(W_b2)
            W_d2[0] ^= 0x80000000
            nb2 = find_neutral_bits(W_b2, W_d2, max_round=10)
            if len(nb2) > len(best_nb):
                best_nb = nb2
                W_base = W_b2
                W_diff = W_d2
                print(f"    Attempt {attempt}: found {len(nb2)} neutral bits (best so far)")
            if len(best_nb) >= 20:
                break

        neutral_bits = best_nb
        neutral_set = set(neutral_bits)
        by_word = defaultdict(int)
        for (w, b) in neutral_bits:
            by_word[w] += 1
        n_words_with_neutral = len(by_word)
        De17_base = compute_De17_from_msgs(W_base, W_diff)
        print(f"    Final: {len(neutral_bits)} neutral bits, De17=0x{De17_base:08x}")

    # ── Cube attack experiment ──
    print("\n" + "=" * 72)
    print("[2] Cube attack: neutral bits vs non-neutral bits as cube variables")
    print("=" * 72)
    print("    Testing output bits j = 0, 8, 16, 24 of De17")
    print("    Cube dimensions k = 2, 3, 4, 5, 6")
    print("    20 random cube selections per (k, type)")
    print("    N=50 evaluations per superpoly check")

    output_bits = [0, 8, 16, 24]
    cube_dims = [2, 3, 4, 5, 6]
    N_CUBES = 20
    N_EVAL = 50
    N_LIN = 15

    results = {
        'cross_neutral': defaultdict(lambda: {'constant': 0, 'linear': 0.0, 'total': 0}),
        'same_neutral': defaultdict(lambda: {'constant': 0, 'linear': 0.0, 'total': 0}),
        'non_neutral': defaultdict(lambda: {'constant': 0, 'linear': 0.0, 'total': 0}),
    }

    for j in output_bits:
        print(f"\n  Output bit j={j}:")
        for k in cube_dims:
            cross_c, cross_l, cross_n = 0, 0.0, 0
            same_c, same_l, same_n = 0, 0.0, 0
            non_c, non_l, non_n = 0, 0.0, 0

            for trial in range(N_CUBES):
                elapsed = time.time() - t0
                if elapsed > 270:
                    break

                # ── Cross-word neutral cubes ──
                cube = select_cross_word(neutral_bits, k)
                if cube is not None:
                    cube_set = set(cube)
                    remaining = [b for b in neutral_bits if b not in cube_set]
                    if len(remaining) > 20:
                        remaining = random.sample(remaining, 20)

                    vals = superpoly_eval(W_base, W_diff, cube, remaining, j, N_EVAL)
                    if check_constant(vals):
                        cross_c += 1
                    if k <= 5 and len(remaining) > 0:
                        lr = blr_linearity_test(W_base, W_diff, cube,
                                                remaining[:8], j, N_LIN)
                        cross_l += lr
                    else:
                        cross_l += 0.5
                    cross_n += 1

                # ── Same-word neutral cubes ──
                cube = select_same_word(neutral_bits, k)
                if cube is not None:
                    cube_set = set(cube)
                    remaining = [b for b in neutral_bits if b not in cube_set]
                    if len(remaining) > 20:
                        remaining = random.sample(remaining, 20)

                    vals = superpoly_eval(W_base, W_diff, cube, remaining, j, N_EVAL)
                    if check_constant(vals):
                        same_c += 1
                    if k <= 5 and len(remaining) > 0:
                        lr = blr_linearity_test(W_base, W_diff, cube,
                                                remaining[:8], j, N_LIN)
                        same_l += lr
                    else:
                        same_l += 0.5
                    same_n += 1

                # ── Non-neutral control ──
                cube = select_non_neutral(neutral_set, k)
                if cube is not None:
                    # Use neutral bits as free variables for the superpoly
                    remaining = list(neutral_bits)
                    if len(remaining) > 20:
                        remaining = random.sample(remaining, 20)

                    vals = superpoly_eval(W_base, W_diff, cube, remaining, j, N_EVAL)
                    if check_constant(vals):
                        non_c += 1
                    if k <= 5 and len(remaining) > 0:
                        lr = blr_linearity_test(W_base, W_diff, cube,
                                                remaining[:8], j, N_LIN)
                        non_l += lr
                    else:
                        non_l += 0.5
                    non_n += 1

            def fmt(c, l, n):
                if n == 0:
                    return "N/A (0 trials)"
                return f"const={c/n*100:5.1f}%  linear={l/n*100:5.1f}%  (n={n})"

            print(f"    k={k}: cross-word: {fmt(cross_c, cross_l, cross_n)}")
            print(f"         same-word:  {fmt(same_c, same_l, same_n)}")
            print(f"         non-neutr:  {fmt(non_c, non_l, non_n)}")

            for key, c, l, n in [('cross_neutral', cross_c, cross_l, cross_n),
                                  ('same_neutral', same_c, same_l, same_n),
                                  ('non_neutral', non_c, non_l, non_n)]:
                results[key][k]['constant'] += c
                results[key][k]['linear'] += l
                results[key][k]['total'] += n

            if time.time() - t0 > 270:
                print("    [time limit approaching]")
                break
        if time.time() - t0 > 270:
            break

    # ── Summary ──
    print("\n" + "=" * 72)
    print("SUMMARY")
    print("=" * 72)
    print(f"\nNeutral bits found: {len(neutral_bits)}")
    print(f"Words with neutral bits: {sorted(by_word.keys())}")

    print(f"\n{'Type':<25s} {'k':>3s} {'n':>5s} {'ConstR%':>8s} {'LinearR%':>9s}")
    print("-" * 55)
    for label, key in [("Cross-word neutral", "cross_neutral"),
                        ("Same-word neutral", "same_neutral"),
                        ("Non-neutral (ctrl)", "non_neutral")]:
        for k in cube_dims:
            r = results[key][k]
            if r['total'] > 0:
                cr = r['constant'] / r['total'] * 100
                lr = r['linear'] / r['total'] * 100
                print(f"  {label:<23s} {k:>3d} {r['total']:>5d} {cr:>7.1f}% {lr:>8.1f}%")

    def aggregate(key):
        tc, tl, tt = 0, 0.0, 0
        for k in cube_dims:
            r = results[key][k]
            tc += r['constant']
            tl += r['linear']
            tt += r['total']
        if tt == 0:
            return 0, 0, 0
        return tc / tt * 100, tl / tt * 100, tt

    cross_cr, cross_lr, cross_n = aggregate('cross_neutral')
    same_cr, same_lr, same_n = aggregate('same_neutral')
    non_cr, non_lr, non_n = aggregate('non_neutral')

    print(f"\nAggregate:")
    print(f"  Cross-word neutral: const={cross_cr:.1f}%  linear={cross_lr:.1f}%  (n={cross_n})")
    print(f"  Same-word neutral:  const={same_cr:.1f}%  linear={same_lr:.1f}%  (n={same_n})")
    print(f"  Non-neutral ctrl:   const={non_cr:.1f}%  linear={non_lr:.1f}%  (n={non_n})")

    # ── Verdict ──
    print("\n" + "=" * 72)
    print("VERDICT")
    print("=" * 72)

    cross_adv_c = cross_cr - non_cr if non_n > 0 else cross_cr
    cross_adv_l = cross_lr - non_lr if non_n > 0 else cross_lr
    same_adv_c = same_cr - non_cr if non_n > 0 else same_cr
    same_adv_l = same_lr - non_lr if non_n > 0 else same_lr

    print(f"\n  Cross-word neutral vs non-neutral control:")
    print(f"    Constant rate advantage: {cross_adv_c:+.1f}%")
    print(f"    Linearity rate advantage: {cross_adv_l:+.1f}%")
    print(f"\n  Same-word neutral vs non-neutral control:")
    print(f"    Constant rate advantage: {same_adv_c:+.1f}%")
    print(f"    Linearity rate advantage: {same_adv_l:+.1f}%")
    print(f"\n  Cross-word vs same-word neutral:")
    print(f"    Constant rate diff: {cross_cr - same_cr:+.1f}%")
    print(f"    Linearity rate diff: {cross_lr - same_lr:+.1f}%")

    # Check total experiments run
    total_n = cross_n + same_n + non_n

    if total_n == 0:
        verdict = "DEAD"
        reason = ("Could not find enough neutral bits to run the experiment. "
                  "The synergy between neutral bits and cube attack cannot be "
                  "established because neutrality is too rare for random messages.")
    elif cross_adv_c > 10 or cross_adv_l > 10 or same_adv_c > 10 or same_adv_l > 10:
        verdict = "ALIVE"
        reason = ("Neutral-bit cubes produce measurably lower-degree superpolys "
                  "than non-neutral cubes. The C9 neutral bits provide genuinely "
                  "useful cube variables for algebraic analysis of De17.")
    else:
        verdict = "DEAD"
        reason = ("Neutral-bit cubes show no significant advantage over "
                  "non-neutral cubes. De17 remains high-degree regardless of "
                  "whether cube variables preserve the differential through "
                  "round 16. Neutrality does not simplify algebraic structure "
                  "at round 17.")

    print(f"\n  >>> C9+C12 Synergy: {verdict} <<<")
    print(f"  Reason: {reason}")

    elapsed = time.time() - t0
    print(f"\n  Total runtime: {elapsed:.1f}s")
    print("=" * 72)


if __name__ == "__main__":
    main()
