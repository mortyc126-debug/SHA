"""
BTE THEORY: Develop the Bi-Temporal Element as a rigorous framework.

1. PROVE parity conservation for AND_carry analytically
2. Test parity in FULL SHA-256 (not simplified BTE)
3. Search for OTHER invariants of BTE
4. Quantify information destruction by V-channel
"""

import random
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, schedule
)
from bte import BTE, AND_carry, XOR_carry, OR_carry


def prove_parity_conservation():
    """
    THEOREM: BTE with V=AND_carry conserves parity.

    PROOF ATTEMPT:

    Let state = [x_0, x_1, ..., x_{n-1}]
    H-step: h_k = XOR of x_{σ_i(k)} for permutations σ_i
    V-step: corrected_k = h_k XOR carry_k
            carry_0 = 0 (seed)
            carry_k = h_{k-1} AND carry_{k-1}

    Parity of corrected = Σ (h_k XOR carry_k) mod 2
                        = Σ h_k  XOR  Σ carry_k  (mod 2)

    Part 1: Σ h_k mod 2
      H is XOR of permuted copies. For each permutation σ:
        Σ_k x_{σ(k)} = Σ_k x_k  (σ is a bijection)
      So: Σ_k h_k = Σ_k [XOR of x_{σ_i(k)} for i=1..m]
                   = XOR_{i=1}^{m} [Σ_k x_{σ_i(k)}]
                   = XOR_{i=1}^{m} [Σ_k x_k]
                   = m * (Σ x_k) mod 2

      For SHA-256 Sig0/Sig1: m=3 (three rotations XORed).
      3 * (Σ x_k) mod 2 = Σ x_k mod 2 (since 3 is odd).
      So: parity(H-output) = parity(input). H PRESERVES parity.

    Part 2: Σ carry_k mod 2
      carry_0 = 0
      carry_k = h_{k-1} AND carry_{k-1}

      AND_carry: carry can only be 1 if BOTH h_{k-1}=1 AND carry_{k-1}=1.
      Once carry becomes 0, it stays 0 forever.
      So carry = [0, 0, ..., 0, 1, 1, ..., 1, 0, 0, ..., 0]
      (a single run of 1s, possibly empty)

      Wait — that's not right. carry_k = h_{k-1} AND carry_{k-1}.
      carry_0 = 0, so carry_1 = h_0 AND 0 = 0.
      carry_2 = h_1 AND carry_1 = h_1 AND 0 = 0.
      ...
      ALL carries are 0 when seed = 0!

      So with seed=0: V-step does NOTHING. corrected = h_state.
      Parity is conserved trivially because V has no effect.

      Let me re-examine with seed=1...
    """
    print("=" * 80)
    print("THEOREM: Parity conservation in BTE")
    print("=" * 80)

    n = 32

    # Test with seed=0: carries should all be 0
    bte = BTE(n, [random.randint(0,1) for _ in range(n)])
    sig0_perms = [[(i+2)%n for i in range(n)],
                  [(i+13)%n for i in range(n)],
                  [(i+22)%n for i in range(n)]]

    h_state = bte.evolve_H(sig0_perms)
    corrected, carries = bte.evolve_V(h_state, AND_carry, seed=0)

    print(f"  With seed=0: carry vector = {carries[:10]}...")
    print(f"  All zeros: {all(c == 0 for c in carries)}")
    print(f"  → V does nothing with seed=0!")
    print(f"  → Parity conservation is TRIVIAL (H preserves parity, V is identity)")

    # With seed=1: carries propagate until h_{k-1}=0
    corrected1, carries1 = bte.evolve_V(h_state, AND_carry, seed=1)
    print(f"\n  With seed=1: carry vector = {carries1[:10]}...")
    print(f"  Sum of carries: {sum(carries1)}")

    # Parity of carries with seed=1
    carry_parity = sum(carries1) % 2
    print(f"  Parity of carry vector: {carry_parity}")

    # Test many seeds and initial states
    N = 1000
    parity_ok = {0: 0, 1: 0}
    for seed_val in [0, 1]:
        count = 0
        for trial in range(N):
            rng = random.Random(trial)
            init = [rng.randint(0, 1) for _ in range(n)]
            bte = BTE(n, list(init))
            p_before = sum(init) % 2
            bte.step(sig0_perms, AND_carry, seed=seed_val)
            p_after = sum(bte.state) % 2
            if p_before == p_after:
                count += 1
        parity_ok[seed_val] = count

    print(f"\n  Parity conservation test (N={N}):")
    print(f"    seed=0: {parity_ok[0]}/{N}")
    print(f"    seed=1: {parity_ok[1]}/{N}")

    # The initial test showed 1000/1000 with seed=0.
    # That's because V does nothing with AND_carry and seed=0!
    # The "conservation law" was trivial.

    # BUT: SHA-256 doesn't use seed=0. The carry seed comes from
    # the GENERATE bits of the addition. Let me reconsider.

    print(f"\n  REVISION: With AND_carry and seed=0, V is identity.")
    print(f"  The 1000/1000 parity conservation was TRIVIAL.")
    print(f"  Need to test with REALISTIC carry seeds (from SHA-256).")


def test_sha256_parity():
    """
    Does FULL SHA-256 conserve parity of any register/combination?

    Test: parity(a[r]) vs parity(a[r+1]) across many messages.
    """
    print("\n" + "=" * 80)
    print("FULL SHA-256: Parity conservation tests")
    print("=" * 80)

    N = 2000
    tests = {
        'parity(a)': lambda s: sum((s[0] >> k) & 1 for k in range(32)) % 2,
        'parity(e)': lambda s: sum((s[4] >> k) & 1 for k in range(32)) % 2,
        'parity(a XOR e)': lambda s: sum(((s[0] ^ s[4]) >> k) & 1 for k in range(32)) % 2,
        'parity(all 8 regs)': lambda s: sum(sum((s[i] >> k) & 1 for k in range(32)) for i in range(8)) % 2,
        'parity(a+e mod2)': lambda s: ((s[0] + s[4]) & MASK) & 1,
        'bit0(a)': lambda s: s[0] & 1,
        'bit0(a) XOR bit0(e)': lambda s: (s[0] & 1) ^ (s[4] & 1),
    }

    for test_name, test_fn in tests.items():
        # For each round pair (r, r+1): does test_fn(state[r]) predict test_fn(state[r+1])?
        conserved_count = 0
        total = 0

        for seed in range(N):
            rng = random.Random(seed)
            M = [rng.randint(0, MASK) for _ in range(16)]
            states, _ = sha256_round_trace(M)

            for r in range(8, 60):  # Skip early rounds (not stabilized)
                val_r = test_fn(states[r])
                val_r1 = test_fn(states[r+1])
                if val_r == val_r1:
                    conserved_count += 1
                total += 1

        rate = conserved_count / total
        print(f"  {test_name:>25}: P(same) = {rate:.4f} (random = 0.5000)")


def experiment_bte_general_invariants():
    """
    Search for invariants of BTE systematically.

    An invariant I(state) satisfies: I(state) = I(evolved_state) for all states.

    For BTE on F_2^n, possible invariants:
    - Linear: I(x) = Σ c_k x_k mod 2 for some mask c
    - Quadratic: I(x) = Σ c_{ij} x_i x_j mod 2
    - Pattern: I(x) = 1 iff specific bits match specific pattern

    Let's search for LINEAR invariants: c such that c · x = c · Evolve(x) mod 2.
    This means: c is in the LEFT null space of (Evolve_matrix - Identity).
    """
    print("\n" + "=" * 80)
    print("BTE: Systematic search for linear invariants")
    print("=" * 80)

    n = 32
    sig0_perms = [[(i+2)%n for i in range(n)],
                  [(i+13)%n for i in range(n)],
                  [(i+22)%n for i in range(n)]]

    # Build the evolution matrix over GF(2)
    # For each basis vector e_k, compute Evolve(e_k)
    # This gives the columns of the evolution matrix

    for V_name, V_func, seed_val in [
        ("AND_carry seed=0", AND_carry, 0),
        ("AND_carry seed=1", AND_carry, 1),
        ("XOR_carry seed=0", XOR_carry, 0),
    ]:
        evolve_matrix = []  # n × n matrix over GF(2)

        for basis_k in range(n):
            init = [0] * n
            init[basis_k] = 1
            bte = BTE(n, list(init))
            bte.step(sig0_perms, V_func, seed=seed_val)
            evolve_matrix.append(list(bte.state))

        # Transpose to get column representation
        # evolve_matrix[k] = image of e_k

        # Compute (E - I) mod 2
        diff_matrix = []
        for k in range(n):
            row = [evolve_matrix[k][j] ^ (1 if j == k else 0) for j in range(n)]
            diff_matrix.append(row)

        # Find left null space of diff_matrix: vectors c with c · (E-I) = 0
        # Equivalently: null space of (E-I)^T
        # Transpose diff_matrix
        diff_T = [[diff_matrix[i][j] for i in range(n)] for j in range(n)]

        # Gaussian elimination to find rank
        m = [list(row) for row in diff_T]
        rank = 0
        for col in range(n):
            pivot = None
            for row in range(rank, n):
                if m[row][col] == 1:
                    pivot = row
                    break
            if pivot is None:
                continue
            m[rank], m[pivot] = m[pivot], m[rank]
            for row in range(n):
                if row != rank and m[row][col] == 1:
                    for c in range(n):
                        m[row][c] ^= m[rank][c]
            rank += 1

        null_dim = n - rank
        print(f"  {V_name:>25}: rank(E-I) = {rank}/{n}, invariant dimension = {null_dim}")

        if null_dim > 0:
            # Find the invariant vectors
            # Free variables
            pivot_rows = {}
            for r in range(rank):
                for c in range(n):
                    if m[r][c] == 1:
                        pivot_rows[c] = r
                        break

            free_cols = [c for c in range(n) if c not in pivot_rows]
            print(f"    Free variables: {free_cols}")

            for fc in free_cols[:3]:  # Show first 3
                inv = [0] * n
                inv[fc] = 1
                for pc, pr in pivot_rows.items():
                    if m[pr][fc] == 1:
                        inv[pc] = 1
                hw = sum(inv)
                print(f"    Invariant: {inv} (HW={hw})")


def experiment_info_destruction():
    """
    Quantify information destroyed by V-channel in one BTE step.

    Method: count how many distinct states map to the same output.
    For a map f: F_2^n → F_2^n:
      - If bijective: 0 info lost (each output has 1 preimage)
      - If k-to-1: log2(k) bits lost per step

    For BTE: Evolve is a map F_2^n → F_2^n.
    Count: how many pairs (x, y) with x ≠ y have Evolve(x) = Evolve(y)?
    """
    print("\n" + "=" * 80)
    print("BTE: Information destruction quantification")
    print("=" * 80)

    # For small n, we can enumerate all 2^n states
    for n in [8, 12, 16]:
        sig0_perms = [[(i+2)%n for i in range(n)],
                      [(i+5)%n for i in range(n)],  # Adapted for smaller n
                      [(i+int(n*0.7))%n for i in range(n)]]

        for V_name, V_func, seed_val in [
            ("AND seed=0", AND_carry, 0),
            ("XOR seed=0", XOR_carry, 0),
        ]:
            output_count = {}
            total_states = 2 ** n

            for state_int in range(total_states):
                init = [(state_int >> k) & 1 for k in range(n)]
                bte = BTE(n, list(init))
                bte.step(sig0_perms, V_func, seed=seed_val)
                out = sum(bte.state[k] << k for k in range(n))
                output_count[out] = output_count.get(out, 0) + 1

            distinct_outputs = len(output_count)
            max_preimages = max(output_count.values())
            avg_preimages = total_states / distinct_outputs
            info_lost = n - (distinct_outputs.bit_length() - 1)  # Rough estimate

            print(f"  n={n:>2}, {V_name:>12}: {distinct_outputs}/{total_states} distinct outputs, "
                  f"max_preimage={max_preimages}, avg={avg_preimages:.1f}")

            if n <= 8:
                # Show preimage size distribution
                size_dist = {}
                for count in output_count.values():
                    size_dist[count] = size_dist.get(count, 0) + 1
                print(f"    Preimage sizes: {dict(sorted(size_dist.items()))}")


if __name__ == "__main__":
    prove_parity_conservation()
    test_sha256_parity()
    experiment_bte_general_invariants()
    experiment_info_destruction()
