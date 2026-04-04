"""
Step 4: Prediction Levels — measuring our understanding of the "coin"

The coin: survive(m) = ANF(T1)[m] ⊕ ANF(T2)[m] ⊕ ANF(carry)[m]

We measure prediction accuracy at increasing levels of model sophistication:

  Level 0: Random guess                    → 50% (baseline)
  Level 1: Use degree of m                 → ?
  Level 2: Use previous round presence     → ?
  Level 3: Know T1 only                    → ?
  Level 4: Know T1 and T2                  → ? (bit 0: should be 100%)
  Level 5: Know all three                  → 100% (trivial)

The gap between levels shows WHERE the understanding resides.
Building the model = closing these gaps.
"""

import numpy as np
from collections import defaultdict
from step1_monomial_genealogy import (
    N, MASK, N_MSG, N_INPUT, N_TOTAL, IV, K,
    make_constant, make_input_bit, make_const_word,
    tt_xor, tt_and, tt_not,
    mobius_transform, monomial_set,
    symbolic_rotr, symbolic_shr, symbolic_xor_word,
    symbolic_ch, symbolic_maj, symbolic_sigma0, symbolic_sigma1,
    symbolic_add, symbolic_add_no_carry
)


def run_all_rounds(R_max):
    """Run R_max rounds, collect ANFs of all components at each round."""
    state_words = []
    for reg in range(8):
        word = []
        for bit in range(N):
            val = (IV[reg] >> bit) & 1
            word.append(make_constant(val))
        state_words.append(word)

    W_msg = []
    for w in range(N_MSG):
        word = []
        for bit in range(N):
            var_idx = w * N + bit
            word.append(make_input_bit(var_idx))
        W_msg.append(word)

    a, b, c, d, e, f, g, h = state_words

    rounds_data = []

    for r in range(R_max):
        w_r = W_msg[r] if r < len(W_msg) else [make_constant(0)] * N
        k_r = make_const_word(K[r % len(K)])

        sig1_e = symbolic_sigma1(e)
        ch_efg = symbolic_ch(e, f, g)
        sig0_a = symbolic_sigma0(a)
        maj_abc = symbolic_maj(a, b, c)

        t1_1, _ = symbolic_add(h, sig1_e)
        t1_2, _ = symbolic_add(t1_1, ch_efg)
        t1_3, _ = symbolic_add(t1_2, k_r)
        T1, _ = symbolic_add(t1_3, w_r)
        T2, _ = symbolic_add(sig0_a, maj_abc)
        a_new, carry_a = symbolic_add(T1, T2)
        e_new, carry_e = symbolic_add(d, T1)

        # Store ANFs for each bit
        round_info = {'T1': [], 'T2': [], 'carry': [], 'a_new': [], 'a_prev': [],
                       'h': [], 'sig1': [], 'ch': [], 'sig0': [], 'maj': []}

        for bit in range(N):
            round_info['T1'].append(mobius_transform(T1[bit]))
            round_info['T2'].append(mobius_transform(T2[bit]))
            round_info['a_new'].append(mobius_transform(a_new[bit]))
            round_info['a_prev'].append(mobius_transform(a[bit]))
            round_info['h'].append(mobius_transform(h[bit]))
            round_info['sig1'].append(mobius_transform(sig1_e[bit]))
            round_info['ch'].append(mobius_transform(ch_efg[bit]))
            round_info['sig0'].append(mobius_transform(sig0_a[bit]))
            round_info['maj'].append(mobius_transform(maj_abc[bit]))

            if bit > 0:
                round_info['carry'].append(mobius_transform(carry_a[bit - 1]))
            else:
                round_info['carry'].append(np.zeros(N_TOTAL, dtype=np.uint8))

        rounds_data.append(round_info)

        # Update state
        h = [g[k].copy() for k in range(N)]
        g = [f[k].copy() for k in range(N)]
        f = [e[k].copy() for k in range(N)]
        e = e_new
        d = [c[k].copy() for k in range(N)]
        c = [b[k].copy() for k in range(N)]
        b = [a[k].copy() for k in range(N)]
        a = a_new

    return rounds_data


def prediction_accuracy(predicted, actual):
    """Fraction of correct predictions."""
    return np.mean(predicted == actual)


def measure_prediction_levels(rounds_data, R, bit):
    """Measure prediction accuracy at all levels for round R, bit position."""

    if R >= len(rounds_data):
        return

    data = rounds_data[R]
    actual = data['a_new'][bit]  # ground truth: which monomials are in a_new

    n_monomials = N_TOTAL
    n_present = np.sum(actual)

    # Exclude constant monomial (index 0) for cleaner statistics
    mask = np.ones(N_TOTAL, dtype=bool)
    # mask[0] = False  # keep all for simplicity

    print(f"\n  Round {R}, a_new[{bit}]: {n_present}/{n_monomials} monomials present "
          f"({n_present/n_monomials*100:.1f}%)")

    # ==== Level 0: Random guess ====
    # Best random: predict majority class
    if n_present > n_monomials / 2:
        pred_0 = np.ones(n_monomials, dtype=np.uint8)
    else:
        pred_0 = np.zeros(n_monomials, dtype=np.uint8)
    acc_0 = prediction_accuracy(pred_0, actual)
    print(f"  Level 0 (majority guess):     {acc_0*100:>6.2f}%")

    # ==== Level 1: Use degree ====
    # For each degree, predict the majority class at that degree
    degrees = np.array([bin(m).count('1') for m in range(n_monomials)])
    pred_1 = np.zeros(n_monomials, dtype=np.uint8)
    for deg in range(N_INPUT + 1):
        deg_mask = degrees == deg
        if np.sum(deg_mask) == 0:
            continue
        deg_present = np.sum(actual[deg_mask])
        deg_total = np.sum(deg_mask)
        # Predict majority class for this degree
        if deg_present > deg_total / 2:
            pred_1[deg_mask] = 1
    acc_1 = prediction_accuracy(pred_1, actual)
    print(f"  Level 1 (degree):             {acc_1*100:>6.2f}%")

    # ==== Level 2: Previous round presence ====
    if R > 0:
        prev_a = rounds_data[R-1]['a_new'][bit]
        # Predict: monomial survives if it was present last round
        acc_2a = prediction_accuracy(prev_a, actual)
        # Also try: predict present if NOT in previous (complement)
        acc_2b = prediction_accuracy(1 - prev_a, actual)
        acc_2 = max(acc_2a, acc_2b)
        print(f"  Level 2 (prev round same):    {acc_2a*100:>6.2f}%")
        print(f"  Level 2 (prev round invert):  {acc_2b*100:>6.2f}%")
    else:
        print(f"  Level 2: N/A (first round)")

    # ==== Level 3a: Know only T1 ====
    # If we know T1[m] but not T2 or carry:
    # a_new[m] = T1[m] ⊕ T2[m] ⊕ carry[m]
    # Best prediction: a_new[m] = T1[m] (assume T2⊕carry = 0 most likely)
    t1 = data['T1'][bit]
    acc_3a_same = prediction_accuracy(t1, actual)
    acc_3a_inv = prediction_accuracy(1 - t1, actual)
    acc_3a = max(acc_3a_same, acc_3a_inv)
    print(f"  Level 3a (T1 = answer):       {acc_3a_same*100:>6.2f}%")
    print(f"  Level 3a (T1 inverted):       {acc_3a_inv*100:>6.2f}%")

    # ==== Level 3b: Know only T2 ====
    t2 = data['T2'][bit]
    acc_3b_same = prediction_accuracy(t2, actual)
    acc_3b_inv = prediction_accuracy(1 - t2, actual)
    acc_3b = max(acc_3b_same, acc_3b_inv)
    print(f"  Level 3b (T2 = answer):       {acc_3b_same*100:>6.2f}%")
    print(f"  Level 3b (T2 inverted):       {acc_3b_inv*100:>6.2f}%")

    # ==== Level 3c: Know individual components ====
    # h, sig1, ch, sig0, maj separately
    for name in ['h', 'sig1', 'ch', 'sig0', 'maj']:
        comp = data[name][bit]
        acc_same = prediction_accuracy(comp, actual)
        acc_inv = prediction_accuracy(1 - comp, actual)
        acc = max(acc_same, acc_inv)
        better = "same" if acc_same >= acc_inv else "inv"
        print(f"  Level 3c ({name:>4} {better}):       {acc*100:>6.2f}%")

    # ==== Level 4a: Know T1 and T2 (not carry) ====
    # a_new = T1 ⊕ T2 ⊕ carry
    # If we know T1, T2: predict a_new = T1 ⊕ T2 (assume carry = 0)
    pred_4a = t1 ^ t2
    acc_4a = prediction_accuracy(pred_4a, actual)
    print(f"  Level 4a (T1⊕T2, no carry):  {acc_4a*100:>6.2f}%")

    # ==== Level 4b: Know T1 and carry (not T2) ====
    carry = data['carry'][bit]
    pred_4b = t1 ^ carry
    acc_4b = prediction_accuracy(pred_4b, actual)
    print(f"  Level 4b (T1⊕carry, no T2):  {acc_4b*100:>6.2f}%")

    # ==== Level 4c: Know T2 and carry (not T1) ====
    pred_4c = t2 ^ carry
    acc_4c = prediction_accuracy(pred_4c, actual)
    print(f"  Level 4c (T2⊕carry, no T1):  {acc_4c*100:>6.2f}%")

    # ==== Level 5: Know all three ====
    pred_5 = t1 ^ t2 ^ carry
    acc_5 = prediction_accuracy(pred_5, actual)
    print(f"  Level 5  (T1⊕T2⊕carry):     {acc_5*100:>6.2f}%  {'✓' if acc_5 > 0.999 else '✗'}")

    # ==== SUMMARY ====
    print(f"\n  PREDICTION LADDER:")
    print(f"  {'Level':>25} {'Accuracy':>9} {'Gap to next':>12}")
    levels = [
        ("0: Random", acc_0),
        ("1: Degree", acc_1),
        ("3a: T1 alone", acc_3a),
        ("3b: T2 alone", acc_3b),
        ("4a: T1⊕T2", acc_4a),
        ("5: T1⊕T2⊕carry", acc_5),
    ]
    for i, (name, acc) in enumerate(levels):
        if i < len(levels) - 1:
            gap = levels[i+1][1] - acc
            print(f"  {name:>25} {acc*100:>8.2f}%  {gap*100:>+11.2f}%")
        else:
            print(f"  {name:>25} {acc*100:>8.2f}%  (complete)")

    # What fraction of understanding is in each component?
    total_gap = acc_5 - acc_0
    if total_gap > 0:
        print(f"\n  UNDERSTANDING DECOMPOSITION (total gap = {total_gap*100:.2f}%):")
        print(f"    Degree alone:     {(acc_1-acc_0)/total_gap*100:>6.1f}% of understanding")
        print(f"    T1 alone:         {(acc_3a-acc_0)/total_gap*100:>6.1f}%")
        print(f"    T2 alone:         {(acc_3b-acc_0)/total_gap*100:>6.1f}%")
        print(f"    T1⊕T2 combined:   {(acc_4a-acc_0)/total_gap*100:>6.1f}%")
        print(f"    +carry (complete): {(acc_5-acc_0)/total_gap*100:>6.1f}%")


def main():
    print("Computing all rounds...")
    rounds_data = run_all_rounds(R_max=7)
    print("Done.\n")

    print("="*80)
    print("PREDICTION LEVELS: How well can we predict the coin?")
    print("="*80)

    for R in [1, 2, 3, 4, 5, 6]:
        print(f"\n{'='*80}")
        print(f"ROUND {R}")
        print(f"{'='*80}")
        for bit in [0, N-1]:
            measure_prediction_levels(rounds_data, R, bit)


if __name__ == "__main__":
    main()
