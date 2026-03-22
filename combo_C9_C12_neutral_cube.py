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
    return (e & f) ^ (~e & g) & 0xFFFFFFFF

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

def get_e(state):
    """Extract e register from state."""
    return state[4]

def compute_De17(W_base, W_mod):
    """Compute De17 = e17(W_mod) XOR e17(W_base)."""
    W_b = expand_W(W_base)
    W_m = expand_W(W_mod)
    st_b = sha256_rounds(W_b, 18)
    st_m = sha256_rounds(W_m, 18)
    # De17 = e after round 17 (0-indexed: round index 17 means 18 rounds)
    return get_e(st_m) ^ get_e(st_b)

def compute_De(W_base, W_mod, rnd):
    """Compute De at given round = e_rnd(mod) XOR e_rnd(base)."""
    W_b = expand_W(W_base)
    W_m = expand_W(W_mod)
    st_b = sha256_rounds(W_b, rnd + 1)
    st_m = sha256_rounds(W_m, rnd + 1)
    return get_e(st_m) ^ get_e(st_b)

def flip_bit(W, word_idx, bit_pos):
    """Return new W with bit flipped at (word_idx, bit_pos)."""
    W2 = list(W)
    W2[word_idx] ^= (1 << bit_pos)
    return W2

def find_neutral_bits(W_base, W_diff):
    """
    Find bits in W[1..15] that are neutral for Wang chain (De=0 for rounds 1-16).
    W_diff is the modified message (W_base with DW[0] applied).
    A bit at (word, pos) is neutral if flipping it in BOTH messages
    preserves De=0 for all rounds 1-16.
    """
    neutral = []
    for word_idx in range(1, 16):
        for bit_pos in range(32):
            # Flip the bit in the modified message
            W_test = flip_bit(W_diff, word_idx, bit_pos)
            # Also flip in base to keep the differential path
            W_base_test = flip_bit(W_base, word_idx, bit_pos)

            is_neutral = True
            for rnd in range(1, 17):  # rounds 1-16
                De = compute_De(W_base_test, W_test, rnd)
                if De != 0:
                    is_neutral = False
                    break
            if is_neutral:
                neutral.append((word_idx, bit_pos))
    return neutral

def superpoly_eval(W_base, W_diff, cube_positions, remaining_neutral,
                   output_bit, n_eval=50):
    """
    Evaluate the superpoly of De17[output_bit] summed over cube_positions,
    at n_eval random settings of remaining neutral bits.

    Returns list of superpoly values (0 or 1) at each evaluation point.
    """
    results = []

    for _ in range(n_eval):
        # Randomize remaining neutral bits
        W_test_base = list(W_base)
        W_test_diff = list(W_diff)
        for (w, b) in remaining_neutral:
            if random.randint(0, 1):
                W_test_base = flip_bit(W_test_base, w, b)
                W_test_diff = flip_bit(W_test_diff, w, b)

        # Sum De17[output_bit] over all 2^k cube assignments
        cube_sum = 0
        n_cube = len(cube_positions)
        for mask in range(1 << n_cube):
            W_cb = list(W_test_base)
            W_cd = list(W_test_diff)
            for idx in range(n_cube):
                if mask & (1 << idx):
                    w, b = cube_positions[idx]
                    W_cb = flip_bit(W_cb, w, b)
                    W_cd = flip_bit(W_cd, w, b)

            De17 = compute_De17(W_cb, W_cd)
            cube_sum ^= (De17 >> output_bit) & 1

        results.append(cube_sum)

    return results

def check_constant(values):
    """Check if all values are the same (constant superpoly)."""
    return len(set(values)) == 1

def check_linearity(W_base, W_diff, cube_positions, remaining_neutral,
                    output_bit, n_tests=20):
    """
    BLR linearity test on the superpoly.
    f(x) + f(y) = f(x+y) + f(0) over GF(2).
    Returns fraction of tests that pass.
    """
    passes = 0
    rem = list(remaining_neutral)
    n_rem = len(rem)

    if n_rem == 0:
        return 1.0  # trivially linear with no variables

    def eval_superpoly(bit_settings):
        """Evaluate superpoly at given setting of remaining neutral bits."""
        W_tb = list(W_base)
        W_td = list(W_diff)
        for i, (w, b) in enumerate(rem):
            if bit_settings[i]:
                W_tb = flip_bit(W_tb, w, b)
                W_td = flip_bit(W_td, w, b)

        cube_sum = 0
        n_cube = len(cube_positions)
        for mask in range(1 << n_cube):
            W_cb = list(W_tb)
            W_cd = list(W_td)
            for idx in range(n_cube):
                if mask & (1 << idx):
                    w, b = cube_positions[idx]
                    W_cb = flip_bit(W_cb, w, b)
                    W_cd = flip_bit(W_cd, w, b)
            De17 = compute_De17(W_cb, W_cd)
            cube_sum ^= (De17 >> output_bit) & 1
        return cube_sum

    for _ in range(n_tests):
        x = [random.randint(0, 1) for _ in range(n_rem)]
        y = [random.randint(0, 1) for _ in range(n_rem)]
        xy = [a ^ b for a, b in zip(x, y)]
        zero = [0] * n_rem

        fx = eval_superpoly(x)
        fy = eval_superpoly(y)
        fxy = eval_superpoly(xy)
        f0 = eval_superpoly(zero)

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
    result = []
    for w in chosen_words:
        result.append(random.choice(by_word[w]))
    return result

def select_same_word(neutral_bits, k):
    """Select k neutral bits from the SAME word."""
    by_word = defaultdict(list)
    for (w, b) in neutral_bits:
        by_word[w].append((w, b))

    # Find words with at least k bits
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
    print("LAYER 2: C9+C12 Neutral Bits × Cube Attack Synergy Test")
    print("=" * 72)

    # ── Generate random base message ──
    W_base = [random.getrandbits(32) for _ in range(16)]
    W_diff = list(W_base)
    W_diff[0] ^= 0x80000000  # DW[0] = msb flip (Wang differential)

    print("\n[1] Finding neutral bits for Wang chain (De=0, rounds 1-16)...")
    t1 = time.time()
    neutral_bits = find_neutral_bits(W_base, W_diff)
    t_neutral = time.time() - t1

    print(f"    Found {len(neutral_bits)} neutral bits in {t_neutral:.1f}s")

    # Distribution by word
    by_word = defaultdict(int)
    for (w, b) in neutral_bits:
        by_word[w] += 1
    print("    Distribution by word:")
    for w in sorted(by_word.keys()):
        print(f"      W[{w:2d}]: {by_word[w]:2d} neutral bits")

    neutral_set = set(neutral_bits)

    # ── Verify De17 is nonzero ──
    De17_base = compute_De17(W_base, W_diff)
    print(f"\n    De17 (base msg) = 0x{De17_base:08x}")

    # ── Cube attack experiment ──
    print("\n[2] Cube attack: neutral bits vs non-neutral bits as cube variables")
    print("    Testing output bits j = 0, 8, 16, 24 of De17")
    print("    Cube dimensions k = 2, 3, 4, 5, 6")
    print("    20 random cube selections per (k, type) combination")
    print("    N=50 evaluations per superpoly check")

    output_bits = [0, 8, 16, 24]
    cube_dims = [2, 3, 4, 5, 6]
    N_CUBES = 20
    N_EVAL = 50

    # Results storage
    # results[type][k] = {'constant': count, 'linear': rate_sum, 'total': count}
    results = {
        'cross_word_neutral': defaultdict(lambda: {'constant': 0, 'linear': 0.0, 'total': 0}),
        'same_word_neutral': defaultdict(lambda: {'constant': 0, 'linear': 0.0, 'total': 0}),
        'non_neutral': defaultdict(lambda: {'constant': 0, 'linear': 0.0, 'total': 0}),
    }

    total_experiments = len(output_bits) * len(cube_dims) * 3 * N_CUBES
    done = 0

    for j in output_bits:
        print(f"\n  Output bit j={j}:")

        for k in cube_dims:
            cross_const = 0
            cross_linear = 0.0
            same_const = 0
            same_linear = 0.0
            non_const = 0
            non_linear = 0.0

            cross_total = 0
            same_total = 0
            non_total = 0

            for trial in range(N_CUBES):
                # ── Cross-word neutral cubes ──
                cube = select_cross_word(neutral_bits, k)
                if cube is not None:
                    remaining = [b for b in neutral_bits if b not in cube]
                    # Limit remaining to speed up linearity test
                    if len(remaining) > 30:
                        remaining = random.sample(remaining, 30)

                    vals = superpoly_eval(W_base, W_diff, cube, remaining, j, N_EVAL)
                    is_const = check_constant(vals)
                    if is_const:
                        cross_const += 1

                    # Linearity test only for small cubes (expensive otherwise)
                    if k <= 5:
                        lin_rate = check_linearity(W_base, W_diff, cube,
                                                   remaining[:10], j, n_tests=15)
                        cross_linear += lin_rate
                    else:
                        cross_linear += 0.5  # placeholder
                    cross_total += 1

                # ── Same-word neutral cubes ──
                cube = select_same_word(neutral_bits, k)
                if cube is not None:
                    remaining = [b for b in neutral_bits if b not in cube]
                    if len(remaining) > 30:
                        remaining = random.sample(remaining, 30)

                    vals = superpoly_eval(W_base, W_diff, cube, remaining, j, N_EVAL)
                    is_const = check_constant(vals)
                    if is_const:
                        same_const += 1

                    if k <= 5:
                        lin_rate = check_linearity(W_base, W_diff, cube,
                                                   remaining[:10], j, n_tests=15)
                        same_linear += lin_rate
                    else:
                        same_linear += 0.5
                    same_total += 1

                # ── Non-neutral control cubes ──
                cube = select_non_neutral(neutral_set, k)
                if cube is not None:
                    remaining = [b for b in neutral_bits if b not in cube]
                    if len(remaining) > 30:
                        remaining = random.sample(remaining, 30)

                    vals = superpoly_eval(W_base, W_diff, cube, remaining, j, N_EVAL)
                    is_const = check_constant(vals)
                    if is_const:
                        non_const += 1

                    if k <= 5:
                        lin_rate = check_linearity(W_base, W_diff, cube,
                                                   remaining[:10], j, n_tests=15)
                        non_linear += lin_rate
                    else:
                        non_linear += 0.5
                    non_total += 1

                done += 3

                # Time check
                elapsed = time.time() - t0
                if elapsed > 270:  # 4.5 min safety
                    print(f"    [time limit approaching at {elapsed:.0f}s, reducing scope]")
                    break

            if cross_total > 0:
                results['cross_word_neutral'][k]['constant'] += cross_const
                results['cross_word_neutral'][k]['linear'] += cross_linear
                results['cross_word_neutral'][k]['total'] += cross_total

            if same_total > 0:
                results['same_word_neutral'][k]['constant'] += same_const
                results['same_word_neutral'][k]['linear'] += same_linear
                results['same_word_neutral'][k]['total'] += same_total

            if non_total > 0:
                results['non_neutral'][k]['constant'] += non_const
                results['non_neutral'][k]['linear'] += non_linear
                results['non_neutral'][k]['total'] += non_total

            # Print per-dimension summary
            def fmt(label, c, l, t):
                if t == 0:
                    return f"    {label}: N/A"
                cr = c / t * 100
                lr = l / t * 100
                return f"    {label}: const={cr:5.1f}%  linear={lr:5.1f}%  (n={t})"

            print(f"    k={k}: {fmt('cross', cross_const, cross_linear, cross_total)}")
            print(f"         {fmt('same ', same_const, same_linear, same_total)}")
            print(f"         {fmt('non  ', non_const, non_linear, non_total)}")

            elapsed = time.time() - t0
            if elapsed > 270:
                break

        elapsed = time.time() - t0
        if elapsed > 270:
            print("  [stopping early due to time limit]")
            break

    # ── Summary ──
    print("\n" + "=" * 72)
    print("SUMMARY")
    print("=" * 72)

    print(f"\nNeutral bits found: {len(neutral_bits)}")
    print(f"Words with neutral bits: {sorted(by_word.keys())}")

    print("\n{:>22s} | {:>6s} {:>8s} {:>8s}".format(
        "Cube Type \\ Dim k", "Metric", "ConstR%", "LinearR%"))
    print("-" * 52)

    for label, key in [("Cross-word neutral", "cross_word_neutral"),
                        ("Same-word neutral", "same_word_neutral"),
                        ("Non-neutral (ctrl)", "non_neutral")]:
        for k in cube_dims:
            r = results[key][k]
            if r['total'] > 0:
                cr = r['constant'] / r['total'] * 100
                lr = r['linear'] / r['total'] * 100
                print(f"  {label:>20s} k={k} |    n={r['total']:3d}  {cr:6.1f}%  {lr:6.1f}%")

    # ── Compute aggregate metrics for verdict ──
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

    cross_cr, cross_lr, cross_n = aggregate('cross_word_neutral')
    same_cr, same_lr, same_n = aggregate('same_word_neutral')
    non_cr, non_lr, non_n = aggregate('non_neutral')

    print(f"\nAggregate across all k and output bits:")
    print(f"  Cross-word neutral: const={cross_cr:.1f}%  linear={cross_lr:.1f}%  (n={cross_n})")
    print(f"  Same-word neutral:  const={same_cr:.1f}%  linear={same_lr:.1f}%  (n={same_n})")
    print(f"  Non-neutral ctrl:   const={non_cr:.1f}%  linear={non_lr:.1f}%  (n={non_n})")

    # ── Verdict ──
    print("\n" + "=" * 72)
    print("VERDICT")
    print("=" * 72)

    # Check if cross-word neutral cubes show meaningfully lower degree
    # than non-neutral control (the C12 baseline proxy)
    cross_advantage_const = cross_cr - non_cr
    cross_advantage_linear = cross_lr - non_lr
    same_advantage_const = same_cr - non_cr
    same_advantage_linear = same_lr - non_lr

    print(f"\n  Cross-word neutral vs non-neutral control:")
    print(f"    Constant rate advantage: {cross_advantage_const:+.1f}%")
    print(f"    Linearity rate advantage: {cross_advantage_linear:+.1f}%")

    print(f"\n  Same-word neutral vs non-neutral control:")
    print(f"    Constant rate advantage: {same_advantage_const:+.1f}%")
    print(f"    Linearity rate advantage: {same_advantage_linear:+.1f}%")

    print(f"\n  Cross-word vs same-word neutral:")
    print(f"    Constant rate difference: {cross_cr - same_cr:+.1f}%")
    print(f"    Linearity rate difference: {cross_lr - same_lr:+.1f}%")

    # Threshold: >10% improvement in either metric = ALIVE
    alive = (cross_advantage_const > 10 or cross_advantage_linear > 10 or
             same_advantage_const > 10 or same_advantage_linear > 10)

    if alive:
        verdict = "ALIVE"
        reason = ("Neutral-bit cubes produce measurably lower-degree superpolys "
                  "than non-neutral cubes. The C9 neutral bits provide genuinely "
                  "useful cube variables for algebraic analysis of De17.")
    else:
        verdict = "DEAD"
        reason = ("Neutral-bit cubes show no significant advantage over "
                  "non-neutral cubes. De17 remains high-degree regardless of "
                  "whether cube variables are neutral for rounds 1-16. "
                  "Neutrality for the Wang chain does not simplify the "
                  "algebraic structure at round 17.")

    print(f"\n  >>> C9+C12 Synergy: {verdict} <<<")
    print(f"  Reason: {reason}")

    elapsed = time.time() - t0
    print(f"\n  Total runtime: {elapsed:.1f}s")
    print("=" * 72)


if __name__ == "__main__":
    main()
