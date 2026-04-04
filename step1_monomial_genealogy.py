"""
Step 1: Monomial Genealogy — trace the BIRTH of every monomial

Instead of looking at the final ANF (which looks "random"), we track
the ANF of every intermediate bit at every round.

For each monomial in the final output, we can answer:
  - At which round was it BORN?
  - Which OPERATION created it? (Ch, Maj, carry, rotation, XOR)
  - From which PARENT monomials?

This reveals the HISTORY hidden inside the "random" polynomial.

Method: Symbolic truth-table computation.
  - Each bit = truth table of 2^16 entries (one per input)
  - Operations on truth tables: XOR = bitwise XOR, AND = bitwise AND
  - Addition mod 2^n: bit-by-bit with carry chain
  - After each round: convert truth tables to ANF via Möbius transform
"""

import numpy as np
from collections import Counter, defaultdict

# ============================================================
# Parameters (same mini-SHA as step0)
# ============================================================
N = 4           # word size
MASK = (1 << N) - 1
N_MSG = 4       # message words
N_INPUT = N_MSG * N  # 16 input bits
N_TOTAL = 1 << N_INPUT  # 65536

# IV and K constants
IV = [0x6, 0xB, 0x3, 0xA, 0x5, 0x9, 0x1, 0xF]
K = [0x4, 0x2, 0xB, 0x7, 0xA, 0x3, 0xE, 0x5,
     0x9, 0x1, 0xD, 0x6, 0x0, 0x8, 0xC, 0xF]

# ============================================================
# Symbolic truth table representation
# ============================================================

def make_constant(val):
    """Truth table for a constant bit (0 or 1)."""
    if val:
        return np.ones(N_TOTAL, dtype=np.uint8)
    else:
        return np.zeros(N_TOTAL, dtype=np.uint8)

def make_input_bit(var_idx):
    """Truth table for input variable x_{var_idx}."""
    tt = np.zeros(N_TOTAL, dtype=np.uint8)
    for i in range(N_TOTAL):
        tt[i] = (i >> var_idx) & 1
    return tt

def tt_xor(a, b):
    return a ^ b

def tt_and(a, b):
    return a & b

def tt_not(a):
    return 1 - a

def tt_or(a, b):
    return a | b

# ============================================================
# Möbius transform (same as step0)
# ============================================================

def mobius_transform(truth_table):
    """Compute ANF from truth table."""
    anf = truth_table.copy()
    for i in range(N_INPUT):
        step = 1 << i
        for j in range(N_TOTAL):
            if j & step:
                anf[j] ^= anf[j ^ step]
    return anf

def count_monomials_by_degree(anf):
    """Count monomials at each degree."""
    counts = defaultdict(int)
    for idx in range(N_TOTAL):
        if anf[idx]:
            deg = bin(idx).count('1')
            counts[deg] += 1
    return dict(counts)

def monomial_set(anf):
    """Return set of monomial indices present in ANF."""
    return set(np.where(anf == 1)[0])

# ============================================================
# Build symbolic state (truth tables for all bits)
# ============================================================

def init_state():
    """Initialize state with IV constants and message input variables."""
    # State = list of 8 words, each word = list of N truth tables (one per bit)
    state = []
    for reg in range(8):
        word = []
        for bit in range(N):
            val = (IV[reg] >> bit) & 1
            word.append(make_constant(val))
        state.append(word)

    # Message words: W[0..3], each bit = input variable
    W = []
    for w in range(N_MSG):
        word = []
        for bit in range(N):
            var_idx = w * N + bit
            word.append(make_input_bit(var_idx))
        W.append(word)

    return state, W

def symbolic_rotr(word, r):
    """Rotate right: bit k of result = bit (k+r) mod N of input."""
    return [word[(k + r) % N] for k in range(N)]

def symbolic_shr(word, r):
    """Shift right: bit k of result = bit (k+r) of input if k+r < N, else 0."""
    result = []
    for k in range(N):
        if k + r < N:
            result.append(word[k + r].copy())
        else:
            result.append(make_constant(0))
    return result

def symbolic_xor_word(a, b):
    return [tt_xor(a[k], b[k]) for k in range(N)]

def symbolic_ch(e, f, g):
    """Ch(e,f,g) = (e AND f) XOR (NOT e AND g)"""
    return [tt_xor(tt_and(e[k], f[k]), tt_and(tt_not(e[k]), g[k])) for k in range(N)]

def symbolic_maj(a, b, c):
    """Maj(a,b,c) = (a AND b) XOR (a AND c) XOR (b AND c)"""
    return [tt_xor(tt_xor(tt_and(a[k], b[k]), tt_and(a[k], c[k])), tt_and(b[k], c[k])) for k in range(N)]

def symbolic_sigma0(x):
    return symbolic_xor_word(symbolic_xor_word(symbolic_rotr(x, 1), symbolic_rotr(x, 2)), symbolic_rotr(x, 3))

def symbolic_sigma1(x):
    return symbolic_xor_word(symbolic_xor_word(symbolic_rotr(x, 1), symbolic_rotr(x, 3)), symbolic_shr(x, 1))

def symbolic_add(a, b):
    """Addition mod 2^N with carry tracking. Returns (sum_word, carry_bits)."""
    result = []
    carries = []
    carry = make_constant(0)
    for k in range(N):
        # sum[k] = a[k] XOR b[k] XOR carry
        s = tt_xor(tt_xor(a[k], b[k]), carry)
        result.append(s)
        # carry[k+1] = MAJ(a[k], b[k], carry)
        new_carry = tt_xor(tt_xor(tt_and(a[k], b[k]), tt_and(a[k], carry)), tt_and(b[k], carry))
        carries.append(new_carry)
        carry = new_carry
    return result, carries

def symbolic_add_no_carry(a, b):
    """Addition, return just the sum."""
    s, _ = symbolic_add(a, b)
    return s

def make_const_word(val):
    """Make a word of constant truth tables."""
    return [make_constant((val >> k) & 1) for k in range(N)]

# ============================================================
# Main: Run rounds and track monomial births
# ============================================================

def run_genealogy(R_max=8):
    state, W_msg = init_state()
    a, b, c, d, e, f, g, h = state

    # Track ANF at each round for output bits (a and e)
    # round_anfs[r][(reg, bit)] = set of monomial indices
    round_anfs = {}

    # Record ANF of initial state
    round_anfs[-1] = {}
    for bit in range(N):
        anf_a = mobius_transform(a[bit])
        anf_e = mobius_transform(e[bit])
        round_anfs[-1][('a', bit)] = monomial_set(anf_a)
        round_anfs[-1][('e', bit)] = monomial_set(anf_e)

    # Track intermediate ANFs to identify birth sources
    births_by_round = {}  # round -> {bit: {source: count}}

    for r in range(R_max):
        w_r = W_msg[r] if r < len(W_msg) else [make_constant(0)] * N
        k_r = make_const_word(K[r % len(K)])

        # Compute intermediates symbolically
        sig1_e = symbolic_sigma1(e)
        ch_efg = symbolic_ch(e, f, g)
        sig0_a = symbolic_sigma0(a)
        maj_abc = symbolic_maj(a, b, c)

        # T1 = h + sig1(e) + ch(e,f,g) + K[r] + W[r]
        t1_1, c1 = symbolic_add(h, sig1_e)
        t1_2, c2 = symbolic_add(t1_1, ch_efg)
        t1_3, c3 = symbolic_add(t1_2, k_r)
        T1, c4 = symbolic_add(t1_3, w_r)

        # T2 = sig0(a) + maj(a,b,c)
        T2, c5 = symbolic_add(sig0_a, maj_abc)

        # New state
        h_new = [g[k].copy() for k in range(N)]
        g_new = [f[k].copy() for k in range(N)]
        f_new = [e[k].copy() for k in range(N)]
        e_new = symbolic_add_no_carry(d, T1)
        d_new = [c[k].copy() for k in range(N)]
        c_new = [b[k].copy() for k in range(N)]
        b_new = [a[k].copy() for k in range(N)]
        a_new = symbolic_add_no_carry(T1, T2)

        # Compute ANFs and track births
        round_anfs[r] = {}
        births = {}

        for bit in range(N):
            # a_new ANF
            anf_a = mobius_transform(a_new[bit])
            new_monos_a = monomial_set(anf_a)
            round_anfs[r][('a', bit)] = new_monos_a

            # e_new ANF
            anf_e = mobius_transform(e_new[bit])
            new_monos_e = monomial_set(anf_e)
            round_anfs[r][('e', bit)] = new_monos_e

            # Identify NEW monomials (born this round)
            # Compare with previous round's a and e ANFs
            if r == 0:
                # All non-constant monomials are new
                prev_a = round_anfs[-1].get(('a', bit), set())
                prev_e = round_anfs[-1].get(('e', bit), set())
            else:
                prev_a = round_anfs[r-1].get(('a', bit), set())
                prev_e = round_anfs[r-1].get(('e', bit), set())

            born_a = new_monos_a - prev_a
            died_a = prev_a - new_monos_a
            born_e = new_monos_e - prev_e
            died_e = prev_e - new_monos_e

            births[('a', bit)] = {'born': len(born_a), 'died': len(died_a),
                                  'total': len(new_monos_a)}
            births[('e', bit)] = {'born': len(born_e), 'died': len(died_e),
                                  'total': len(new_monos_e)}

            # Identify WHERE new monomials come from
            # Check: which intermediates contain these new monomials?
            if bit == 0 and len(born_a) > 0 and len(born_a) < 500:
                anf_ch = mobius_transform(ch_efg[bit])
                anf_maj = mobius_transform(maj_abc[bit])
                anf_sig1 = mobius_transform(sig1_e[bit])
                anf_sig0 = mobius_transform(sig0_a[bit])

                ch_monos = monomial_set(anf_ch)
                maj_monos = monomial_set(anf_maj)

                # For carry births: monomials in T1 that aren't in any
                # XOR-only component
                anf_T1 = mobius_transform(T1[bit])
                t1_monos = monomial_set(anf_T1)

                # Carry contributions for each addition
                carry_monos = set()
                for carries in [c1, c2, c3, c4, c5]:
                    if bit > 0:  # carry affects bits 1+
                        anf_c = mobius_transform(carries[bit-1])
                        carry_monos |= monomial_set(anf_c)

                births[('a', bit)]['from_ch'] = len(born_a & ch_monos)
                births[('a', bit)]['from_maj'] = len(born_a & maj_monos)
                births[('a', bit)]['from_carry'] = len(born_a & carry_monos)

        births_by_round[r] = births

        # Update state
        a, b, c, d, e, f, g, h = a_new, b_new, c_new, d_new, e_new, f_new, g_new, h_new

    return round_anfs, births_by_round

def print_results(round_anfs, births_by_round, R_max=8):
    print(f"MONOMIAL GENEALOGY: Mini-SHA (n={N}, msg={N_MSG}, input={N_INPUT})")
    print(f"{'='*70}")

    # Birth/death table
    print(f"\n{'ROUND':>6} {'BIT':>5} {'TOTAL':>7} {'BORN':>7} {'DIED':>7} {'NET':>7} {'%NEW':>6}")
    print(f"{'-'*50}")

    for r in range(R_max):
        births = births_by_round[r]
        for reg in ['a', 'e']:
            for bit in [0, N-1]:  # Just bit 0 and MSB
                b = births[(reg, bit)]
                net = b['born'] - b['died']
                pct_new = b['born'] / max(b['total'], 1) * 100
                print(f"  R={r:>2}  {reg}[{bit}]  {b['total']:>6}  {b['born']:>6}  {b['died']:>6}  {net:>+6}  {pct_new:>5.1f}%")
        print()

    # Monomial degree evolution for a[0]
    print(f"\n{'='*70}")
    print(f"DEGREE EVOLUTION: a[0] through rounds")
    print(f"{'='*70}")
    print(f"{'Round':>6} {'Deg0':>5} {'Deg1':>5} {'Deg2':>5} {'Deg3':>5} "
          f"{'Deg4':>5} {'Deg5':>5} {'Deg6':>5} {'Deg7':>5} {'Deg8+':>6} {'MaxDeg':>7}")

    for r in range(-1, R_max):
        monos = round_anfs[r].get(('a', 0), set())
        if not monos:
            continue
        deg_counts = Counter()
        max_deg = 0
        for m in monos:
            deg = bin(m).count('1')
            deg_counts[deg] += 1
            max_deg = max(max_deg, deg)

        row = f"  {r:>4}"
        for d in range(8):
            row += f" {deg_counts.get(d, 0):>5}"
        higher = sum(v for k, v in deg_counts.items() if k >= 8)
        row += f" {higher:>6}"
        row += f" {max_deg:>7}"
        print(row)

    # Birth sources for a[0] (when trackable)
    print(f"\n{'='*70}")
    print(f"BIRTH SOURCES: a[0] — which operation created new monomials")
    print(f"{'='*70}")
    for r in range(R_max):
        b = births_by_round[r].get(('a', 0), {})
        if 'from_ch' in b:
            total_born = b['born']
            from_ch = b.get('from_ch', 0)
            from_maj = b.get('from_maj', 0)
            from_carry = b.get('from_carry', 0)
            other = total_born - from_ch - from_maj - from_carry
            if total_born > 0:
                print(f"  R={r}: born={total_born}, "
                      f"Ch={from_ch}({from_ch/total_born*100:.0f}%), "
                      f"Maj={from_maj}({from_maj/total_born*100:.0f}%), "
                      f"carry={from_carry}({from_carry/total_born*100:.0f}%), "
                      f"mixed={other}({other/total_born*100:.0f}%)")

    # The key question: MONOMIAL SURVIVAL
    print(f"\n{'='*70}")
    print(f"MONOMIAL SURVIVAL: how many monomials from round r survive to round R_max-1")
    print(f"{'='*70}")

    final_a0 = round_anfs[R_max-1][('a', 0)]
    final_e0 = round_anfs[R_max-1][('e', 0)]

    for r in range(-1, R_max-1):
        monos_a = round_anfs[r].get(('a', 0), set())
        monos_e = round_anfs[r].get(('e', 0), set())

        if monos_a:
            survived_a = len(monos_a & final_a0)
            pct_a = survived_a / len(monos_a) * 100 if monos_a else 0
        else:
            survived_a = 0
            pct_a = 0

        if monos_e:
            survived_e = len(monos_e & final_e0)
            pct_e = survived_e / len(monos_e) * 100 if monos_e else 0
        else:
            survived_e = 0
            pct_e = 0

        print(f"  Round {r:>2}: a[0] {survived_a:>6}/{len(monos_a):>6} survived ({pct_a:>5.1f}%), "
              f"e[0] {survived_e:>6}/{len(monos_e):>6} survived ({pct_e:>5.1f}%)")

    # CHURN: what fraction of monomials is replaced each round
    print(f"\n{'='*70}")
    print(f"CHURN RATE: fraction of monomials replaced each round")
    print(f"{'='*70}")
    print(f"{'Round':>6} {'a[0] churn':>12} {'a[3] churn':>12} {'e[0] churn':>12} {'e[3] churn':>12}")
    for r in range(R_max):
        row = f"  R={r:>2}"
        for reg, bit in [('a',0), ('a',N-1), ('e',0), ('e',N-1)]:
            b = births_by_round[r][(reg, bit)]
            total = max(b['total'], 1)
            churn = (b['born'] + b['died']) / (2 * total)
            row += f" {churn:>11.3f}"
        print(row)


if __name__ == "__main__":
    round_anfs, births = run_genealogy(R_max=8)
    print_results(round_anfs, births, R_max=8)
