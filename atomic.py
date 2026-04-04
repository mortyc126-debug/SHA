"""
ATOMIC LEVEL: Track SHA-256 bit by bit, carry by carry.

"Random" = we don't see the pattern. Not "no pattern."
SHA-256 is DETERMINISTIC. Every bit has a CAUSE.

For ONE specific message M:
  Every carry bit at every addition at every round = DETERMINED.
  The question: what determines it? And is there a PATTERN
  in the determination that we've been averaging away?

Instead of measuring statistics over many M:
  Follow ONE M through all 64 rounds.
  At each step: WHY is this bit 0 or 1?
  The chain of causes = the STRUCTURE of this specific computation.
"""

import random
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, schedule
)


def trace_one_round_atomic(state, r, W_r):
    """
    Trace ONE round at ATOMIC level.
    Return: every intermediate value and every carry bit.
    """
    a, b, c, d, e, f, g, h = state
    mask = MASK

    # Step 1: Sig1(e)
    sig1_e = Sig1(e)

    # Step 2: Ch(e, f, g)
    ch = Ch(e, f, g)

    # Step 3: T1 = h + Sig1(e) + Ch(e,f,g) + K[r] + W[r]
    # Decompose each addition:
    def add_atomic(x, y, name):
        result = (x + y) & mask
        xor_part = x ^ y
        carry_effect = result ^ xor_part

        # Carry chain detail
        carries = []
        generates = []
        propagates = []
        c = 0
        for k in range(32):
            xk = (x >> k) & 1
            yk = (y >> k) & 1
            carries.append(c)
            generates.append(xk & yk)
            propagates.append(xk ^ yk)
            c = (xk & yk) | (xk & c) | (yk & c)

        return {
            'name': name,
            'x': x, 'y': y,
            'result': result,
            'xor': xor_part,
            'carry_effect': carry_effect,
            'carries': carries,
            'generates': generates,
            'propagates': propagates,
            'n_carries': sum(carries),
            'n_generates': sum(generates),
            'n_propagates': sum(propagates),
        }

    add1 = add_atomic(h, sig1_e, 'h + Sig1(e)')
    add2 = add_atomic(add1['result'], ch, '+ Ch')
    add3 = add_atomic(add2['result'], K[r], '+ K[r]')
    add4 = add_atomic(add3['result'], W_r, '+ W[r]')
    T1 = add4['result']

    # Step 4: Sig0(a), Maj(a,b,c)
    sig0_a = Sig0(a)
    maj = Maj(a, b, c)

    add5 = add_atomic(sig0_a, maj, 'Sig0(a) + Maj')
    T2 = add5['result']

    # Step 5: a_new = T1 + T2
    add6 = add_atomic(T1, T2, 'T1 + T2')
    a_new = add6['result']

    # Step 6: e_new = d + T1
    add7 = add_atomic(d, T1, 'd + T1')
    e_new = add7['result']

    return {
        'adds': [add1, add2, add3, add4, add5, add6, add7],
        'a_new': a_new, 'e_new': e_new,
        'T1': T1, 'T2': T2,
        'sig1_e': sig1_e, 'ch': ch,
        'sig0_a': sig0_a, 'maj': maj,
    }


def experiment_atomic_trace():
    """Trace one message through SHA-256 at atomic level."""
    print("=" * 80)
    print("ATOMIC TRACE: One message, every carry, every bit")
    print("=" * 80)

    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]
    states, W = sha256_round_trace(M)

    # Trace rounds 0, 1, 2 in full detail
    for r in range(3):
        trace = trace_one_round_atomic(states[r], r, W[r])

        print(f"\n  === ROUND {r} ===")
        print(f"  State: a=0x{states[r][0]:08x} e=0x{states[r][4]:08x}")

        total_carries = 0
        total_generates = 0
        total_propagates = 0

        for add in trace['adds']:
            print(f"    {add['name']:>15}: carries={add['n_carries']:>2}, "
                  f"gen={add['n_generates']:>2}, prop={add['n_propagates']:>2}, "
                  f"kill={32-add['n_generates']-add['n_propagates']:>2}")
            total_carries += add['n_carries']
            total_generates += add['n_generates']
            total_propagates += add['n_propagates']

        print(f"    {'TOTAL':>15}: carries={total_carries}, gen={total_generates}, prop={total_propagates}")

    # Now: the KEY question. For the carry chain in each addition,
    # what determines WHETHER a carry happens?
    # carry[k] = MAJ(x[k-1], y[k-1], carry[k-1])
    # = 1 iff at least 2 of {x[k-1], y[k-1], carry[k-1]} = 1.
    #
    # The CAUSE chain:
    # carry[0] = 0 (axiom)
    # carry[1] = x[0] AND y[0] (determined by bit 0 of operands)
    # carry[2] = MAJ(x[1], y[1], carry[1])
    #          = MAJ(x[1], y[1], x[0] AND y[0])
    # Each carry = function of ALL LOWER BITS of both operands.

    # But x and y are themselves RESULTS of previous operations.
    # x = h (register value from 3 rounds ago)
    # y = Sig1(e) (rotated + XORed current e)

    # So: carry[k] at round r depends on:
    #   h[0..k-1] (= e[r-3][0..k-1])
    #   Sig1(e[r])[0..k-1] = bits of e[r] at rotated positions

    # This creates a DEPENDENCY WEB: each carry bit at round r
    # depends on specific bits of state at previous rounds.


def experiment_carry_web():
    """
    Build the CARRY DEPENDENCY WEB for one specific computation.

    For one message M: at each round r, each carry bit C[r][add][k]
    depends on specific bits of the state.

    Track: for the FINAL carry (a_new = T1 + T2, bit 31):
    What is its complete dependency tree?
    """
    print("\n" + "=" * 80)
    print("CARRY WEB: Dependency tree of one specific carry bit")
    print("=" * 80)

    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]
    states, W = sha256_round_trace(M)

    # Focus: a_new at round 0, bit 15 (mid-range carry).
    # a_new = T1 + T2.
    # carry[15] of this addition depends on T1[0..14] and T2[0..14].

    r = 0
    trace = trace_one_round_atomic(states[r], r, W[r])

    # The final addition: T1 + T2
    add_final = trace['adds'][5]  # T1 + T2
    T1 = trace['T1']
    T2 = trace['T2']

    print(f"\n  Round 0, a_new = T1 + T2:")
    print(f"  T1 = 0x{T1:08x}")
    print(f"  T2 = 0x{T2:08x}")
    print(f"  a_new = 0x{add_final['result']:08x}")

    # Carry chain for T1 + T2:
    carries = add_final['carries']
    gens = add_final['generates']
    props = add_final['propagates']

    print(f"\n  Carry chain (T1+T2):")
    print(f"  bit | T1  T2 | G   P   K | carry | a_new")
    print(f"  " + "-" * 50)

    for k in range(32):
        t1k = (T1 >> k) & 1
        t2k = (T2 >> k) & 1
        gk = gens[k]
        pk = props[k]
        kk = 1 - gk - pk  # kill
        ck = carries[k]
        ak = (add_final['result'] >> k) & 1

        # WHY is carry[k] what it is?
        if k == 0:
            why = "axiom: carry[0]=0"
        elif ck == 1:
            if gens[k-1]:
                why = f"GENERATE at bit {k-1}"
            elif props[k-1] and carries[k-1]:
                why = f"PROPAGATE from bit {k-1}"
            else:
                why = f"MAJ vote at bit {k-1}"
        else:
            if not gens[k-1] and not props[k-1]:
                why = f"KILL at bit {k-1}"
            elif props[k-1] and not carries[k-1]:
                why = f"PROPAGATE blocked (carry[{k-1}]=0)"
            else:
                why = f"MAJ vote at bit {k-1}"

        if k < 16:  # Show first 16 bits
            print(f"  {k:>3} |  {t1k}   {t2k} |  {gk}   {pk}   {kk} |   {ck}   | {ak}   ← {why}")

    # SEGMENT analysis
    segments = []
    seg_start = 0
    seg_type = 'K' if not gens[0] and not props[0] else ('G' if gens[0] else 'P')

    for k in range(1, 32):
        curr = 'K' if not gens[k] and not props[k] else ('G' if gens[k] else 'P')
        if curr != seg_type or curr == 'G':
            segments.append((seg_start, k - seg_start, seg_type))
            seg_start = k
            seg_type = curr
    segments.append((seg_start, 32 - seg_start, seg_type))

    print(f"\n  Carry segments: (start, length, type)")
    for start, length, stype in segments:
        print(f"    [{start}:{start+length}] len={length} type={stype}")

    print(f"\n  Segment types: G=generate(carry starts), P=propagate(carry flows), K=kill(carry stops)")
    print(f"  Total: {len(segments)} segments")
    print(f"  This EXACT pattern is the 'atomic fingerprint' of this specific computation.")

    # HOW STABLE is this pattern? Does changing M[0] bit 0 change the segments?
    print(f"\n  --- Perturbation: flip M[0] bit 0 ---")
    M2 = list(M); M2[0] ^= 1
    states2, W2 = sha256_round_trace(M2)
    trace2 = trace_one_round_atomic(states2[0], 0, W2[0])
    add_final2 = trace2['adds'][5]

    same_carry = sum(1 for k in range(32) if add_final['carries'][k] == add_final2['carries'][k])
    same_gpk = sum(1 for k in range(32)
                    if add_final['generates'][k] == add_final2['generates'][k]
                    and add_final['propagates'][k] == add_final2['propagates'][k])

    print(f"  Same carry bits: {same_carry}/32")
    print(f"  Same GPK pattern: {same_gpk}/32")
    print(f"  Changed carries: {32-same_carry}/32")

    # Round 0: W[0] changes → T1 changes → carry pattern changes.
    # How MUCH changes depends on WHERE the carry chain is affected.

    # Bit 0 flip in W[0] → T1 bit 0 flips → carry[1] may flip
    # → propagation through P-segments → cascade
    print(f"\n  Cascade: which carry bits flipped?")
    for k in range(32):
        if add_final['carries'][k] != add_final2['carries'][k]:
            print(f"    carry[{k}] flipped ({add_final['carries'][k]} → {add_final2['carries'][k]})")


if __name__ == "__main__":
    experiment_atomic_trace()
    experiment_carry_web()
