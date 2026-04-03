"""
Stage 0, Part 3: Search for AUTO-RELATIONS in SHA-256.

Differentials compare f(x) vs f(x+δ) — TWO messages.
Auto-relations compare state_r vs state_{r+k} — ONE message.

Question: Are there functions g(state_r, state_{r+k}) that equal
a CONSTANT (or a simple function of K,W) for ALL messages M?

If such relations exist, they pass through all 64 rounds by definition,
because they hold for ANY message at ANY round.

Known auto-relation: b_{r+1} = a_r (register shift). This is trivial.
We're looking for NON-TRIVIAL ones.
"""

import random
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, sig0, sig1, schedule, bits
)


def experiment_9_xor_relations():
    """
    Search for XOR relations between registers at different rounds.

    For each pair (reg_i at round r, reg_j at round r+k),
    compute reg_i[r] XOR reg_j[r+k] over many messages.
    If the result is CONSTANT → we found an auto-relation.
    If HW is far from 16 → partial relation (bias).
    """
    N = 500
    print("=" * 80)
    print("EXPERIMENT 9: XOR auto-relations between registers at different rounds")
    print("=" * 80)

    reg_names = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']

    # Collect states for many messages
    all_states = []
    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        states, _ = sha256_round_trace(M)
        all_states.append(states)

    # Test: for round pairs (r, r+k), registers (i, j),
    # is states[r][i] XOR states[r+k][j] constant across messages?
    print("\nSearching for constant XOR relations (r, r+k, reg_i, reg_j)...")
    print("(Showing only cases where all N messages give same value)\n")

    found = 0
    # Check small k first (nearby rounds)
    for k in range(1, 8):
        for r in range(4, 20):  # skip IV rounds, check mid-range
            if r + k > 64:
                continue
            for i in range(8):
                for j in range(8):
                    # Skip trivial shifts
                    if k <= 3 and (
                        (i == 0 and j == k) or  # a→b→c→d
                        (i == 4 and j == 4+k and 4+k < 8)  # e→f→g→h
                    ):
                        continue

                    vals = set()
                    for msg_idx in range(N):
                        v = all_states[msg_idx][r][i] ^ all_states[msg_idx][r+k][j]
                        vals.add(v)
                        if len(vals) > 1:
                            break

                    if len(vals) == 1:
                        found += 1
                        val = list(vals)[0]
                        print(f"  CONSTANT: {reg_names[i]}[{r}] XOR {reg_names[j]}[{r+k}] = 0x{val:08x} (k={k})")

    print(f"\nTotal non-trivial constant XOR relations found: {found}")


def experiment_10_add_relations():
    """
    Search for ADDITIVE relations: reg_i[r] + reg_j[r+k] mod 2^32.

    SHA-256 uses mod-2^32 addition, not XOR. So additive relations
    might exist where XOR ones don't.
    """
    N = 500
    print("\n" + "=" * 80)
    print("EXPERIMENT 10: Additive auto-relations")
    print("=" * 80)

    reg_names = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']

    all_states = []
    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        states, _ = sha256_round_trace(M)
        all_states.append(states)

    # Known: a_{r+1} = T1_r + T2_r, e_{r+1} = d_r + T1_r
    # So: a_{r+1} - e_{r+1} = T2_r - d_r (mod 2^32)
    # And: T2_r = Sig0(a_r) + Maj(a_r, b_r, c_r)
    # And: d_r = a_{r-3}
    # So: a_{r+1} - e_{r+1} = Sig0(a_r) + Maj(a_r, a_{r-1}, a_{r-2}) - a_{r-3}

    # This is a relation involving only a-branch values!
    # Let's verify it and see if it simplifies.

    print("\nVerifying: a[r+1] - e[r+1] = Sig0(a[r]) + Maj(a[r],b[r],c[r]) - d[r]")
    violations = 0
    for msg_idx in range(N):
        states = all_states[msg_idx]
        for r in range(1, 63):
            a_r1 = states[r+1][0]
            e_r1 = states[r+1][4]
            a_r = states[r][0]
            b_r = states[r][1]
            c_r = states[r][2]
            d_r = states[r][3]

            lhs = (a_r1 - e_r1) & MASK
            rhs = (Sig0(a_r) + Maj(a_r, b_r, c_r) - d_r) & MASK

            if lhs != rhs:
                violations += 1
                break
        else:
            continue
        break

    print(f"  Violations: {violations}/{N*62}")
    if violations == 0:
        print(f"  CONFIRMED: a[r+1] - e[r+1] = T2[r] - d[r] for ALL messages, ALL rounds ✓")

    # This means: a[r+1] - e[r+1] depends ONLY on (a_r, a_{r-1}, a_{r-2}, a_{r-3})
    # (because b=a_{r-1}, c=a_{r-2}, d=a_{r-3})
    # This is a 4th-order recurrence in the a-sequence alone!

    print(f"\n  Rewriting with shifts: a[r+1] - e[r+1] = Sig0(a[r]) + Maj(a[r], a[r-1], a[r-2]) - a[r-3]")
    print(f"  This is a relation involving ONLY the a-sequence!")

    # Now: what about e[r+1]?
    # e[r+1] = d[r] + T1[r] = a[r-3] + h[r] + Sig1(e[r]) + Ch(e[r],f[r],g[r]) + K[r] + W[r]
    # = a[r-3] + e[r-3] + Sig1(e[r]) + Ch(e[r], e[r-1], e[r-2]) + K[r] + W[r]
    #
    # This involves BOTH a and e sequences plus K and W.
    # So the e-recurrence is: e[r+1] = a[r-3] + e[r-3] + Sig1(e[r]) + Ch(e[r],e[r-1],e[r-2]) + K[r] + W[r]

    print(f"\n  e[r+1] = a[r-3] + e[r-3] + Sig1(e[r]) + Ch(e[r],e[r-1],e[r-2]) + K[r] + W[r]")

    # Verify this
    violations2 = 0
    for msg_idx in range(min(100, N)):
        states = all_states[msg_idx]
        _, W = sha256_round_trace([random.Random(msg_idx).randint(0, MASK) for _ in range(16)])
        # Recompute with same seed
        rng = random.Random(msg_idx)
        M = [rng.randint(0, MASK) for _ in range(16)]
        states, W = sha256_round_trace(M)

        for r in range(3, 63):
            e_r1 = states[r+1][4]
            a_r3 = states[r-3][0] if r >= 3 else H0[0]
            e_r3 = states[r-3][4] if r >= 3 else H0[4]
            e_r = states[r][4]
            e_r1_prev = states[r-1][4] if r >= 1 else H0[4]
            e_r2 = states[r-2][4] if r >= 2 else H0[4]

            # e[r+1] should equal a[r-3] + e[r-3] + Sig1(e[r]) + Ch(e[r],e[r-1],e[r-2]) + K[r] + W[r]
            computed = (a_r3 + e_r3 + Sig1(e_r) + Ch(e_r, e_r1_prev, e_r2) + K[r] + W[r]) & MASK

            if e_r1 != computed:
                violations2 += 1

    print(f"\n  Verification: {violations2} violations out of {100*60}")
    if violations2 == 0:
        print(f"  CONFIRMED: e[r+1] = a[r-3] + e[r-3] + Sig1(e[r]) + Ch(e[r],e[r-1],e[r-2]) + K[r] + W[r] ✓")


def experiment_11_combined_invariant():
    """
    We now have TWO recurrences:
      (A): a[r+1] - e[r+1] = Sig0(a[r]) + Maj(a[r], a[r-1], a[r-2]) - a[r-3]
      (E): e[r+1] = a[r-3] + e[r-3] + Sig1(e[r]) + Ch(e[r], e[r-1], e[r-2]) + K[r] + W[r]

    From (A): a[r+1] = e[r+1] + Sig0(a[r]) + Maj(a[r], a[r-1], a[r-2]) - a[r-3]

    Substituting (E) into (A):
      a[r+1] = [a[r-3] + e[r-3] + Sig1(e[r]) + Ch(e[r],e[r-1],e[r-2]) + K[r] + W[r]]
             + Sig0(a[r]) + Maj(a[r], a[r-1], a[r-2]) - a[r-3]

      a[r+1] = e[r-3] + Sig1(e[r]) + Ch(e[r],e[r-1],e[r-2]) + K[r] + W[r]
             + Sig0(a[r]) + Maj(a[r], a[r-1], a[r-2])

    This is the FULL recurrence for a alone (given e history and constants).

    Now: the PAIR (a[r], e[r]) forms a coupled dynamical system on (Z/2^32)^2
    with memory depth 3.

    Question: Does the system have any CONSERVED QUANTITY?
    I.e., a function C(a[r], a[r-1], a[r-2], e[r], e[r-1], e[r-2]) that is
    constant (or depends only on K, W) across ALL rounds?
    """
    N = 200
    print("\n" + "=" * 80)
    print("EXPERIMENT 11: Search for conserved quantities in (a,e) dynamics")
    print("=" * 80)

    # Try simple candidates:
    # C1: a[r] + e[r] (mod 2^32)
    # C2: a[r] XOR e[r]
    # C3: a[r] + e[r] + a[r-1] + e[r-1] (running sum)
    # C4: a[r] - e[r] - Sig0(a[r-1]) (from relation A)

    all_states = []
    all_W = []
    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        states, W = sha256_round_trace(M)
        all_states.append(states)
        all_W.append(W)

    candidates = {
        'a+e': lambda s, r: (s[r][0] + s[r][4]) & MASK,
        'a^e': lambda s, r: s[r][0] ^ s[r][4],
        'a-e': lambda s, r: (s[r][0] - s[r][4]) & MASK,
        'a+e+b+f': lambda s, r: (s[r][0] + s[r][4] + s[r][1] + s[r][5]) & MASK,
        'a^e^b^f': lambda s, r: s[r][0] ^ s[r][4] ^ s[r][1] ^ s[r][5],
        'sum_all8': lambda s, r: sum(s[r][i] for i in range(8)) & MASK,
        'xor_all8': lambda s, r: s[r][0] ^ s[r][1] ^ s[r][2] ^ s[r][3] ^ s[r][4] ^ s[r][5] ^ s[r][6] ^ s[r][7],
    }

    print(f"\nFor each candidate C, check if C[r+1] - C[r] is constant across messages:")
    print(f"(Testing rounds 4-60, {N} messages)\n")

    for name, func in candidates.items():
        # Compute C[r] for all messages and rounds
        # Check: is C[r+1] - C[r] the same for all messages?
        deltas_by_round = {}
        for r in range(4, 60):
            vals = set()
            for msg_idx in range(N):
                c_r1 = func(all_states[msg_idx], r+1)
                c_r = func(all_states[msg_idx], r)
                delta = (c_r1 - c_r) & MASK
                vals.add(delta)
                if len(vals) > 1:
                    break
            deltas_by_round[r] = vals

        # How many rounds have constant delta?
        const_rounds = sum(1 for v in deltas_by_round.values() if len(v) == 1)
        print(f"  {name:>12}: {const_rounds}/56 rounds with constant Δ")

    # Now try something deeper: subtract the KNOWN contributions
    # From (A): a[r+1] - e[r+1] - Sig0(a[r]) - Maj(a[r],b[r],c[r]) + d[r] = 0
    # This is always 0 — trivially. Let's find something less trivial.

    # Try: Σ_{r=0}^{R} (a[r] + e[r]) — cumulative sum
    print(f"\n--- Cumulative sum a+e: does it have structure? ---")
    for seed in [0, 1, 2]:
        states = all_states[seed]
        W = all_W[seed]
        cum = 0
        print(f"  Seed {seed}:")
        for r in range(20):
            cum = (cum + states[r][0] + states[r][4]) & MASK
            kw = (K[r] + W[r]) & MASK if r < 64 else 0
            adjusted = (cum - sum(K[i] + W[i] for i in range(r+1))) & MASK
            print(f"    r={r:>2}: cum(a+e)=0x{cum:08x}  cum-ΣKW=0x{adjusted:08x}")


def experiment_12_створка():
    """
    From the methodology (section 202): створочное число
    a[r-4] = e[r] - a[r] + T2[r-1]
    Verified 100%, 640K tests.

    Let's verify this and explore what it means structurally.
    This is a KNOWN auto-relation from the methodology.
    """
    N = 500
    print("\n" + "=" * 80)
    print("EXPERIMENT 12: Verify створочное число (section 202 of methodology)")
    print("=" * 80)
    print("Relation: a[r-4] = e[r] - a[r] + T2[r-1]")

    violations = 0
    total = 0
    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        states, W = sha256_round_trace(M)

        for r in range(5, 64):
            a_r = states[r][0]
            e_r = states[r][4]
            a_r4 = states[r-4][0]
            # T2[r-1] = Sig0(a[r-1]) + Maj(a[r-1], b[r-1], c[r-1])
            a_prev = states[r-1][0]
            b_prev = states[r-1][1]
            c_prev = states[r-1][2]
            T2_prev = (Sig0(a_prev) + Maj(a_prev, b_prev, c_prev)) & MASK

            lhs = a_r4
            rhs = (e_r - a_r + T2_prev) & MASK

            if lhs != rhs:
                violations += 1
            total += 1

    print(f"  Violations: {violations}/{total}")
    if violations == 0:
        print(f"  CONFIRMED: a[r-4] = e[r] - a[r] + T2[r-1] for ALL messages, ALL rounds ✓")
        print(f"\n  This means: e[r] = a[r] + a[r-4] - T2[r-1]")
        print(f"  Or: e[r] - a[r] = a[r-4] - T2[r-1]")
        print(f"  The DIFFERENCE e-a at round r is determined by a[r-4] and T2[r-1]")
        print(f"  T2 depends only on a-branch: T2 = Sig0(a) + Maj(a,b,c) = Sig0(a[r-1]) + Maj(a[r-1],a[r-2],a[r-3])")
        print(f"\n  CONSEQUENCE: e[r] is FULLY determined by a[r], a[r-1], a[r-2], a[r-3], a[r-4]")
        print(f"  The e-sequence carries NO independent information beyond the a-sequence!")
        print(f"  SHA-256 state evolution is a SINGLE 5th-order recurrence in a alone.")


if __name__ == "__main__":
    experiment_9_xor_relations()
    experiment_10_add_relations()
    experiment_11_combined_invariant()
    experiment_12_створка()
