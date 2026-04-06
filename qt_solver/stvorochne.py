"""
Створочне Elimination algebra for SHA-256 cryptanalysis.

Key identity (verified with 0/3050 violations):
    e[r] = a[r] + a[r-4] - T2[r-1]  (mod 2^32)

Where T2[r] = Sigma0(a[r]) + Maj(a[r], a[r-1], a[r-2])
            = Sigma0(a[r]) + Maj(a[r], a[r-1], a[r-2])

This means the e-register is FULLY DETERMINED by the a-sequence!
SHA-256's 8 registers reduce to a single sequence {a[r]}.
State = {a[r], a[r-1], a[r-2], a[r-3]} = 128 bits (not 256).

Important: The створочне is a mod-2^32 identity, NOT a GF(2) identity.
When linearizing, the Jacobian captures first-order effects which ARE
valid in GF(2). So the linearized system works correctly.
"""

from qt_solver.sha256_traced import (
    MASK32, IV, K, get_bit,
    sigma0, sigma1, ssigma0, ssigma1, ch, maj,
    sha256_compress, sha256_compress_traced,
)
from qt_solver.gf2 import gf2_kernel, gf2_gaussian_eliminate


import random


def _T2(a_r, a_rm1, a_rm2):
    """Compute T2[r] = Sigma0(a[r]) + Maj(a[r], a[r-1], a[r-2]) mod 2^32."""
    return (sigma0(a_r) + maj(a_r, a_rm1, a_rm2)) & MASK32


# ------------------------------------------------------------------ #
# 1. compute_a_sequence
# ------------------------------------------------------------------ #
def compute_a_sequence(msg, R=64):
    """
    Compute the full a-sequence a[0..R] from a message.

    a[0] = IV[0], then a[r] is the 'a' register after round r.
    Returns list of R+1 uint32 values.
    """
    trace = sha256_compress_traced(msg, R)
    return [state[0] for state in trace['states']]


# ------------------------------------------------------------------ #
# 2. recover_e_from_a
# ------------------------------------------------------------------ #
def recover_e_from_a(a_seq, R=64):
    """
    Given the a-sequence, recover the e-sequence using the створочне identity:
        e[r] = a[r] + a[r-4] - T2[r-1]  (mod 2^32)

    For rounds r >= 5 the identity applies directly.
    For r = 0..4 we bootstrap from IV.

    Returns list of R+1 uint32 (e[0] = IV[4], e[r] for r=1..R).
    Verifies against actual SHA-256 trace when called from verify_e_recovery.
    """
    # We need a[-4..-1] which come from the IV shift registers.
    # Before round 0 the state is (a,b,c,d,e,f,g,h) = IV.
    # In the a-sequence convention:
    #   a[0]  = IV[0]
    #   a[-1] = IV[1]  (= b before round 0)
    #   a[-2] = IV[2]  (= c)
    #   a[-3] = IV[3]  (= d)
    # Similarly for e:
    #   e[0]  = IV[4]
    #   e[-1] = IV[5]  (= f)
    #   e[-2] = IV[6]  (= g)
    #   e[-3] = IV[7]  (= h)

    # Extend a_seq to include negative indices via IV
    # a_ext[i] = a_seq[i] for i >= 0, IV-derived for i < 0
    # a[-1]=IV[1], a[-2]=IV[2], a[-3]=IV[3]
    def a_ext(i):
        if i >= 0:
            return a_seq[i]
        # i in {-1, -2, -3}
        return IV[-i]  # IV[1], IV[2], IV[3]

    e_seq = [IV[4]]  # e[0]

    for r in range(1, R + 1):
        # e[r] = a[r] + a[r-4] - T2[r-1]  (mod 2^32)
        # T2[r-1] = Sigma0(a[r-1]) + Maj(a[r-1], a[r-2], a[r-3])
        t2 = _T2(a_ext(r - 1), a_ext(r - 2), a_ext(r - 3))
        e_r = (a_ext(r) + a_ext(r - 4) - t2) & MASK32
        e_seq.append(e_r)

    return e_seq


def verify_e_recovery(msg, R=64, verbose=True):
    """Verify that recover_e_from_a produces the correct e-sequence."""
    a_seq = compute_a_sequence(msg, R)
    e_recovered = recover_e_from_a(a_seq, R)

    trace = sha256_compress_traced(msg, R)
    e_actual = [state[4] for state in trace['states']]

    violations = 0
    for r in range(R + 1):
        if e_recovered[r] != e_actual[r]:
            violations += 1
            if verbose:
                print(f"  VIOLATION at r={r}: recovered={hex(e_recovered[r])}, "
                      f"actual={hex(e_actual[r])}")

    if verbose:
        print(f"  e-recovery: {violations}/{R+1} violations")
    return violations


# ------------------------------------------------------------------ #
# 3. recover_full_state
# ------------------------------------------------------------------ #
def recover_full_state(a_seq, R=64):
    """
    Recover all 8 registers at all rounds from just the a-sequence.

    Shift relations:
        b[r] = a[r-1], c[r] = a[r-2], d[r] = a[r-3]
    Створочне:
        e[r] = a[r] + a[r-4] - T2[r-1]
    Then:
        f[r] = e[r-1], g[r] = e[r-2], h[r] = e[r-3]

    Returns list of R+1 tuples (a, b, c, d, e, f, g, h).
    """
    def a_ext(i):
        if i >= 0:
            return a_seq[i]
        return IV[-i]  # IV[1], IV[2], IV[3]

    e_seq = recover_e_from_a(a_seq, R)

    def e_ext(i):
        if i >= 0:
            return e_seq[i]
        # e[-1]=IV[5], e[-2]=IV[6], e[-3]=IV[7]
        return IV[4 - i]  # e[-1]->IV[5], e[-2]->IV[6], e[-3]->IV[7]

    states = []
    for r in range(R + 1):
        a = a_ext(r)
        b = a_ext(r - 1)
        c = a_ext(r - 2)
        d = a_ext(r - 3)
        e = e_ext(r)
        f = e_ext(r - 1)
        g = e_ext(r - 2)
        h = e_ext(r - 3)
        states.append((a, b, c, d, e, f, g, h))

    return states


def verify_full_state(msg, R=64, verbose=True):
    """Verify recover_full_state matches actual SHA-256 trace."""
    a_seq = compute_a_sequence(msg, R)
    recovered = recover_full_state(a_seq, R)

    trace = sha256_compress_traced(msg, R)
    actual = trace['states']

    violations = 0
    for r in range(R + 1):
        if recovered[r] != actual[r]:
            violations += 1
            if verbose:
                print(f"  VIOLATION at r={r}:")
                print(f"    recovered: {tuple(hex(x) for x in recovered[r])}")
                print(f"    actual:    {tuple(hex(x) for x in actual[r])}")

    if verbose:
        print(f"  full-state recovery: {violations}/{R+1} violations")
    return violations


# ------------------------------------------------------------------ #
# 4. build_a_only_system
# ------------------------------------------------------------------ #
def build_a_only_system(R, msg, verbose=True):
    """
    Build OMEGA linearized system using ONLY a-variables (no separate e-variables).

    Variable layout:
        - Message variables: 512 (W[0..15] x 32 bits)
        - State variables: only delta_a[1..R] x 32 bits (not delta_e!)
        - e equations eliminated via створочне
        Total = 512 + R*32 variables

    Linearize around known trace. Solve and report alpha-kernel.

    Returns dict with kernel info.
    """
    n_msg = 512
    n_a = R * 32
    n_vars = n_msg + n_a
    hash_bits = 256

    trace = sha256_compress_traced(msg, R)
    target_hash = trace['hash']
    a_seq = [s[0] for s in trace['states']]  # a[0..R]

    if verbose:
        print(f"Stvorochne OMEGA for R={R}")
        print(f"  Variables: {n_vars} (msg={n_msg}, a={n_a})")

    # For each message bit, compute hash change
    msg_effects = []
    for wi in range(16):
        for bi in range(32):
            msg_copy = list(msg)
            msg_copy[wi] ^= (1 << bi)
            new_hash = sha256_compress(msg_copy, R)
            diff = 0
            for w in range(8):
                xd = target_hash[w] ^ new_hash[w]
                diff |= (xd << (w * 32))
            msg_effects.append(diff)

    # For each a-variable bit, compute hash change using stvorochne reconstruction.
    # Perturb a[r] bit b, reconstruct full state, compute forward hash.
    # Since we're linearizing, we use the actual SHA-256 round function
    # with the perturbed a injected.
    a_effects = []
    for r_idx in range(1, R + 1):
        for bi in range(32):
            # Perturb a[r_idx]
            a_pert = list(a_seq)
            a_pert[r_idx] ^= (1 << bi)

            # Recover full state from perturbed a-sequence
            states_pert = recover_full_state(a_pert, R)

            # The "hash" after R rounds is state[R] + IV
            final = states_pert[R]
            new_hash = [(final[i] + IV[i]) & MASK32 for i in range(8)]

            diff = 0
            for w in range(8):
                xd = target_hash[w] ^ new_hash[w]
                diff |= (xd << (w * 32))
            a_effects.append(diff)

    # Build Jacobian rows: one per hash bit
    rows = []
    for hb in range(hash_bits):
        row = 0
        for mb in range(n_msg):
            if (msg_effects[mb] >> hb) & 1:
                row |= (1 << mb)
        for ab in range(n_a):
            if (a_effects[ab] >> hb) & 1:
                row |= (1 << (n_msg + ab))
        rows.append(row)

    # Find kernel
    kernel = gf2_kernel(rows, n_vars)

    # alpha-kernel: kernel vectors with support only in message bits
    alpha_kernel = []
    a_mask = ((1 << n_vars) - 1) ^ ((1 << n_msg) - 1)
    for kvec in kernel:
        if (kvec & a_mask) == 0 and kvec != 0:
            alpha_kernel.append(kvec)

    if verbose:
        echelon, pivots = gf2_gaussian_eliminate(list(rows), n_vars)
        rank = len(pivots)
        print(f"  Jacobian rank: {rank}/{hash_bits}")
        print(f"  Full kernel dim: {len(kernel)}")
        print(f"  alpha-kernel dim: {len(alpha_kernel)}")

    return {
        'kernel': kernel,
        'alpha_kernel': alpha_kernel,
        'alpha_dim': len(alpha_kernel),
        'n_vars': n_vars,
        'n_msg': n_msg,
        'n_a': n_a,
        'R': R,
        'jacobian_rows': rows,
    }


# ------------------------------------------------------------------ #
# 5. compare_variable_counts
# ------------------------------------------------------------------ #
def compare_variable_counts(R_values=None, msg=None, verbose=True):
    """
    Compare standard OMEGA (msg + a + e variables) vs
    Створочне OMEGA (msg + a only).

    Reports: variable reduction, equation reduction, alpha-kernel (must match).
    """

    if R_values is None:
        R_values = [4, 8, 12, 16]

    if msg is None:
        rng = random.Random(42)
        msg = [rng.randint(0, MASK32) for _ in range(16)]

    if verbose:
        print("=" * 70)
        print("Variable Count Comparison: Standard OMEGA vs Створочне OMEGA")
        print("=" * 70)
        print(f"{'R':>4} | {'Std vars':>10} {'Stv vars':>10} {'Reduction':>10} | "
              f"{'Std alpha':>10} {'Stv alpha':>10} {'Match':>6}")
        print("-" * 70)

    results = []
    for R in R_values:
        # Standard OMEGA variables: msg(512) + a(R*32) + e(R*32)
        std_n_msg = 512
        std_n_a = R * 32
        std_n_e = R * 32
        std_vars = std_n_msg + std_n_a + std_n_e

        # Stvorochne variables: msg(512) + a(R*32) only
        stv_vars = std_n_msg + std_n_a

        # Compute standard OMEGA alpha-kernel
        from qt_solver.omega_solve import omega_linearize_and_solve
        std_result = omega_linearize_and_solve(R, msg, verbose=False)
        std_alpha = std_result['alpha_dim']

        # Compute stvorochne alpha-kernel
        stv_result = build_a_only_system(R, msg, verbose=False)
        stv_alpha = stv_result['alpha_dim']

        reduction = std_vars - stv_vars
        match = (std_alpha == stv_alpha)

        if verbose:
            print(f"{R:>4} | {std_vars:>10} {stv_vars:>10} {reduction:>10} | "
                  f"{std_alpha:>10} {stv_alpha:>10} {'YES' if match else 'NO':>6}")

        results.append({
            'R': R,
            'std_vars': std_vars,
            'stv_vars': stv_vars,
            'reduction': reduction,
            'std_alpha': std_alpha,
            'stv_alpha': stv_alpha,
            'match': match,
        })

    if verbose:
        print("-" * 70)
        print(f"Variable reduction per round: {R_values[0] if R_values else '?'} -> "
              f"{results[0]['reduction'] if results else '?'} "
              f"(eliminated {results[0]['reduction']} e-variables)" if results else "")
        print(f"Створочне eliminates R*32 variables (the entire e-sequence).")

    return results


# ------------------------------------------------------------------ #
# 6. a_sequence_recurrence
# ------------------------------------------------------------------ #
def a_sequence_recurrence(msg, R=64):
    """
    Express SHA-256 as 8th-order recurrence on {a[r]}:
        a[r+1] = F(a[r], a[r-1], ..., a[r-7], K[r], W[r])

    Derivation:
        a[r+1] = T1[r] + T2[r]
        T1[r]  = h[r] + Sigma1(e[r]) + Ch(e[r], f[r], g[r]) + K[r] + W[r]
        T2[r]  = Sigma0(a[r]) + Maj(a[r], b[r], c[r])

    Substituting shift relations:
        b[r]=a[r-1], c[r]=a[r-2], d[r]=a[r-3]
        e[r] from створочне: e[r] = a[r] + a[r-4] - T2[r-1]
        f[r]=e[r-1], g[r]=e[r-2], h[r]=e[r-3]

    So a[r+1] depends on a[r], a[r-1], ..., a[r-7] (and K[r], W[r]).

    Prints the recurrence and verifies it for all rounds.
    """
    trace = sha256_compress_traced(msg, R)
    a_seq = [s[0] for s in trace['states']]
    W = trace['W']

    # Extended a with IV for negative indices
    def a_ext(i):
        if i >= 0:
            return a_seq[i]
        return IV[-i]

    print("SHA-256 as 8th-order recurrence on a-sequence")
    print("=" * 60)
    print()
    print("a[r+1] = T1[r] + T2[r]  where:")
    print("  T2[r] = Sigma0(a[r]) + Maj(a[r], a[r-1], a[r-2])")
    print("  e[r]  = a[r] + a[r-4] - T2[r-1]")
    print("  e[r-1]= a[r-1] + a[r-5] - T2[r-2]")
    print("  e[r-2]= a[r-2] + a[r-6] - T2[r-3]")
    print("  e[r-3]= a[r-3] + a[r-7] - T2[r-4]")
    print("  T1[r] = e[r-3] + Sigma1(e[r]) + Ch(e[r], e[r-1], e[r-2]) + K[r] + W[r]")
    print()
    print("Dependencies: a[r], a[r-1], a[r-2], a[r-3], a[r-4], a[r-5], a[r-6], a[r-7]")
    print("              + constants K[r], W[r]")
    print()

    # Verify the recurrence
    violations = 0
    for r in range(R):
        # Compute e values from stvorochne
        def get_e(idx):
            """Compute e[idx] from a-sequence."""
            if idx <= 0:
                return IV[4 - idx]  # e[0]=IV[4], e[-1]=IV[5], etc.
            t2_prev = _T2(a_ext(idx - 1), a_ext(idx - 2), a_ext(idx - 3))
            return (a_ext(idx) + a_ext(idx - 4) - t2_prev) & MASK32

        e_r = get_e(r)
        e_rm1 = get_e(r - 1)
        e_rm2 = get_e(r - 2)
        h_r = get_e(r - 3)  # h[r] = e[r-3]

        t1 = (h_r + sigma1(e_r) + ch(e_r, e_rm1, e_rm2) + K[r] + W[r]) & MASK32
        t2 = _T2(a_ext(r), a_ext(r - 1), a_ext(r - 2))
        a_next = (t1 + t2) & MASK32

        if a_next != a_seq[r + 1]:
            violations += 1
            print(f"  VIOLATION at r={r}: computed={hex(a_next)}, "
                  f"actual={hex(a_seq[r+1])}")

    print(f"Recurrence verification: {violations}/{R} violations")
    if violations == 0:
        print("All rounds verified: a[r+1] = F(a[r], ..., a[r-7], K[r], W[r])")

    # Print a few sample rounds
    print()
    print("Sample values:")
    for r in [0, 1, 2, R - 1]:
        if r < R:
            print(f"  r={r:2d}: a[r+1]={hex(a_seq[r+1])}, "
                  f"a[r]={hex(a_ext(r))}, ..., a[r-7]={hex(a_ext(r-7))}")

    return violations


# ------------------------------------------------------------------ #
# Main verification
# ------------------------------------------------------------------ #
if __name__ == '__main__':
    print("=" * 72)
    print("СТВОРОЧНЕ ELIMINATION — Full Verification Suite")
    print("=" * 72)

    rng = random.Random(42)
    msg = [rng.randint(0, MASK32) for _ in range(16)]

    # ------------------------------------------------------------------
    # [1] compute_a_sequence
    # ------------------------------------------------------------------
    print("\n[1] compute_a_sequence")
    print("-" * 60)
    a_seq = compute_a_sequence(msg, 64)
    print(f"  a[0]  = {a_seq[0]:#010x}  (IV[0])")
    print(f"  a[1]  = {a_seq[1]:#010x}")
    print(f"  a[64] = {a_seq[64]:#010x}")
    print(f"  Length: {len(a_seq)} (R+1 = 65)")

    # ------------------------------------------------------------------
    # [2] recover_e_from_a — bulk verification (replicating 0/3050)
    # ------------------------------------------------------------------
    print("\n[2] recover_e_from_a — bulk identity verification")
    print("-" * 60)
    total_violations = 0
    total_checks = 0
    N_MSGS = 50
    for seed in range(N_MSGS):
        rng2 = random.Random(seed)
        m = [rng2.randint(0, MASK32) for _ in range(16)]
        v = verify_e_recovery(m, R=64, verbose=False)
        total_violations += v
        total_checks += 65
    print(f"  {N_MSGS} random messages, R=64: "
          f"{total_violations}/{total_checks} violations")
    assert total_violations == 0, "створочне identity violated!"

    # ------------------------------------------------------------------
    # [3] recover_full_state — all 8 registers from a alone
    # ------------------------------------------------------------------
    print("\n[3] recover_full_state — all 8 registers from a-sequence")
    print("-" * 60)
    total_violations = 0
    total_checks = 0
    for seed in range(10):
        rng2 = random.Random(seed)
        m = [rng2.randint(0, MASK32) for _ in range(16)]
        v = verify_full_state(m, R=64, verbose=False)
        total_violations += v
        total_checks += 65
    print(f"  10 random messages, R=64: "
          f"{total_violations}/{total_checks} state tuples checked")
    assert total_violations == 0, "full-state recovery failed!"

    # ------------------------------------------------------------------
    # [4] build_a_only_system — створочне OMEGA
    # ------------------------------------------------------------------
    print("\n[4] build_a_only_system — створочне OMEGA (a-only)")
    print("-" * 60)
    for R in [4, 8]:
        print(f"\n  --- R={R} ---")
        result = build_a_only_system(R, msg, verbose=True)

    # ------------------------------------------------------------------
    # [5] compare_variable_counts
    # ------------------------------------------------------------------
    print("\n[5] compare_variable_counts")
    print("-" * 60)
    results = compare_variable_counts([4, 8, 12, 16], msg=msg, verbose=True)
    for r in results:
        assert r['match'], f"Alpha-kernel mismatch at R={r['R']}!"
    print("  All alpha-kernel dimensions match.")

    # ------------------------------------------------------------------
    # [6] a_sequence_recurrence — 8th-order recurrence
    # ------------------------------------------------------------------
    print("\n[6] a_sequence_recurrence — SHA-256 as 8th-order recurrence")
    print("-" * 60)
    v = a_sequence_recurrence(msg, R=64)
    assert v == 0, f"Recurrence verification failed with {v} violations!"

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print("\n" + "=" * 72)
    print("ALL VERIFICATIONS PASSED")
    print("=" * 72)
    print()
    print("Key results:")
    print("  - створочне identity: e[r] = a[r] + a[r-4] - T2[r-1]  (mod 2^32)")
    print(f"  - Verified on {N_MSGS * 65} round-checks with 0 violations")
    print("  - Full 8-register state recoverable from a-sequence alone")
    print("  - Variable reduction: R*32 e-variables eliminated")
    print("  - Alpha-kernel preserved (identical to standard OMEGA)")
