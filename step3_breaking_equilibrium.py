"""
Step 3: Breaking the Equilibrium

The "random" in SHA = three equal votes: T1, T2, carry ≈ 25% each.
Question: is this balance UNIFORM across all monomials?
Or do certain CLASSES of monomials have skewed votes?

If monomials of degree d, or involving variable set S, have
a different T1/T2/carry balance — that's exploitable structure.

We check:
  1. Vote balance by DEGREE of monomial
  2. Vote balance by which MESSAGE WORDS are involved
  3. Vote balance for bit 0 (carry-free, only 2 votes)
  4. Cross-bit vote correlation (carry chain links bits)
"""

import numpy as np
from collections import Counter, defaultdict
from step1_monomial_genealogy import (
    N, MASK, N_MSG, N_INPUT, N_TOTAL, IV, K,
    make_constant, make_input_bit, make_const_word,
    tt_xor, tt_and, tt_not,
    mobius_transform, monomial_set,
    symbolic_rotr, symbolic_shr, symbolic_xor_word,
    symbolic_ch, symbolic_maj, symbolic_sigma0, symbolic_sigma1,
    symbolic_add, symbolic_add_no_carry
)


def compute_round_state(R_target):
    """Run R_target rounds, return T1, T2, carry, a_new for that round."""
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

    for r in range(R_target + 1):
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

        if r == R_target:
            return T1, T2, carry_a, a_new, e_new, carry_e

        h_n = [g[k].copy() for k in range(N)]
        g_n = [f[k].copy() for k in range(N)]
        f_n = [e[k].copy() for k in range(N)]
        d_n = [c[k].copy() for k in range(N)]
        c_n = [b[k].copy() for k in range(N)]
        b_n = [a[k].copy() for k in range(N)]
        e = e_new
        a = a_new
        h, g, f, d, c, b = h_n, g_n, f_n, d_n, c_n, b_n


def degree_of(monomial_idx):
    return bin(monomial_idx).count('1')

def words_of(monomial_idx):
    """Which message words does this monomial involve?"""
    words = set()
    for v in range(N_INPUT):
        if monomial_idx & (1 << v):
            words.add(v // N)
    return frozenset(words)


def analyze_vote_by_degree(T1, T2, carry_a, bit, R):
    """Check if vote balance differs by monomial degree."""
    anf_T1 = mobius_transform(T1[bit])
    anf_T2 = mobius_transform(T2[bit])
    if bit > 0:
        anf_carry = mobius_transform(carry_a[bit - 1])
    else:
        anf_carry = np.zeros(N_TOTAL, dtype=np.uint8)

    # For each monomial, determine its vote pattern
    # vote = (in_T1, in_T2, in_carry) ∈ {0,1}^3
    # Monomial survives iff vote has odd weight (1 or 3 ones)

    vote_by_degree = defaultdict(lambda: Counter())

    for m in range(N_TOTAL):
        deg = degree_of(m)
        if deg == 0 and m != 0:
            continue
        v = (int(anf_T1[m]), int(anf_T2[m]), int(anf_carry[m]))
        vote_by_degree[deg][v] += 1

    print(f"\n  a_new[{bit}] R={R}: Vote pattern by degree")
    print(f"  {'Deg':>4} {'Total':>6} | {'(0,0,0)':>8} {'(1,0,0)':>8} {'(0,1,0)':>8} {'(0,0,1)':>8} "
          f"{'(1,1,0)':>8} {'(1,0,1)':>8} {'(0,1,1)':>8} {'(1,1,1)':>8} | {'%surv':>6}")

    for deg in sorted(vote_by_degree.keys()):
        if deg > N_INPUT:
            continue
        vc = vote_by_degree[deg]
        total = sum(vc.values())
        if total == 0:
            continue

        # Survivors = those with odd number of 1s in vote
        survived = (vc[(1,0,0)] + vc[(0,1,0)] + vc[(0,0,1)] + vc[(1,1,1)])
        pct = survived / total * 100

        row = f"  {deg:>4} {total:>6} |"
        for v in [(0,0,0), (1,0,0), (0,1,0), (0,0,1), (1,1,0), (1,0,1), (0,1,1), (1,1,1)]:
            count = vc[v]
            row += f" {count/total*100:>7.1f}%"
        row += f" | {pct:>5.1f}%"
        print(row)

    return vote_by_degree


def analyze_vote_by_words(T1, T2, carry_a, bit, R):
    """Check if vote balance differs by which message words are involved."""
    anf_T1 = mobius_transform(T1[bit])
    anf_T2 = mobius_transform(T2[bit])
    if bit > 0:
        anf_carry = mobius_transform(carry_a[bit - 1])
    else:
        anf_carry = np.zeros(N_TOTAL, dtype=np.uint8)

    vote_by_words = defaultdict(lambda: Counter())

    for m in range(1, N_TOTAL):  # skip constant
        ws = words_of(m)
        v = (int(anf_T1[m]), int(anf_T2[m]), int(anf_carry[m]))
        vote_by_words[ws][v] += 1

    print(f"\n  a_new[{bit}] R={R}: Vote pattern by message words involved")
    print(f"  {'Words':>15} {'Total':>6} | {'OnlyT1':>7} {'OnlyT2':>7} {'OnlyC':>7} {'All3':>7} | {'%surv':>6} {'%T1':>5} {'%T2':>5} {'%C':>5}")

    for ws in sorted(vote_by_words.keys(), key=lambda x: (len(x), x)):
        vc = vote_by_words[ws]
        total = sum(vc.values())
        if total < 10:
            continue

        survived = vc[(1,0,0)] + vc[(0,1,0)] + vc[(0,0,1)] + vc[(1,1,1)]
        only_t1 = vc[(1,0,0)]
        only_t2 = vc[(0,1,0)]
        only_c = vc[(0,0,1)]
        all3 = vc[(1,1,1)]

        pct_surv = survived / total * 100
        pct_t1 = only_t1 / max(survived, 1) * 100
        pct_t2 = only_t2 / max(survived, 1) * 100
        pct_c = only_c / max(survived, 1) * 100

        wname = ",".join(f"W{w}" for w in sorted(ws))
        print(f"  {wname:>15} {total:>6} | {only_t1:>7} {only_t2:>7} {only_c:>7} {all3:>7} "
              f"| {pct_surv:>5.1f}% {pct_t1:>4.0f}% {pct_t2:>4.0f}% {pct_c:>4.0f}%")


def analyze_cross_bit_correlation(T1, T2, carry_a, R):
    """Check if survival of the same monomial across different bits is correlated."""
    print(f"\n  Cross-bit survival correlation (R={R}):")
    print(f"  Are votes for the same monomial correlated across bit positions?")

    # For each pair of bits, check: P(both survive) vs P(survive)^2
    survival = {}
    for bit in range(N):
        anf_T1 = mobius_transform(T1[bit])
        anf_T2 = mobius_transform(T2[bit])
        if bit > 0:
            anf_carry = mobius_transform(carry_a[bit - 1])
        else:
            anf_carry = np.zeros(N_TOTAL, dtype=np.uint8)

        # survive[m] = 1 iff monomial m survives in a_new[bit]
        surv = (anf_T1 ^ anf_T2 ^ anf_carry)
        survival[bit] = surv

    print(f"  {'Bit_i':>6} {'Bit_j':>6} {'P(i)':>7} {'P(j)':>7} {'P(i∧j)':>8} {'P(i)*P(j)':>10} {'Excess':>8}")

    for i in range(N):
        for j in range(i+1, N):
            p_i = np.mean(survival[i][1:])  # skip constant monomial
            p_j = np.mean(survival[j][1:])
            p_ij = np.mean(survival[i][1:] & survival[j][1:])
            p_indep = p_i * p_j
            excess = p_ij - p_indep
            print(f"  {i:>6} {j:>6} {p_i:>7.4f} {p_j:>7.4f} {p_ij:>8.4f} {p_indep:>10.4f} {excess:>+8.4f}")


def main():
    for R in [2, 4, 5]:
        print(f"\n{'='*80}")
        print(f"ROUND {R}")
        print(f"{'='*80}")

        T1, T2, carry_a, a_new, e_new, carry_e = compute_round_state(R)

        # 1. Vote by degree
        for bit in [0, N-1]:
            analyze_vote_by_degree(T1, T2, carry_a, bit, R)

        # 2. Vote by words (only for equilibrium rounds)
        if R >= 4:
            analyze_vote_by_words(T1, T2, carry_a, 0, R)
            analyze_vote_by_words(T1, T2, carry_a, N-1, R)

        # 3. Cross-bit correlation
        analyze_cross_bit_correlation(T1, T2, carry_a, R)


if __name__ == "__main__":
    main()
