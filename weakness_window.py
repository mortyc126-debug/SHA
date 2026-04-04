"""
WEAKNESS WINDOW: SHA-256 at R ≤ 12 has structured Hessian.

At R=12: P(H[j][k]=1) ≈ 0.20. 80% of input pairs are affine.
This means: 80% of the degree-2 system is actually degree-1.

Can we IDENTIFY which pairs are affine and exploit this?

If we know WHICH pairs (j,k) have H=0 (affine) vs H=1 (quadratic):
we can partition the system into linear and quadratic parts.
Linear part: solvable in O(n^2).
Quadratic part: only 20% of terms → much sparser → cheaper XL.
"""

import random
from stage0_observe import sha256_round_trace, MASK, H0, K

def sha256_hash_R(M, R):
    states, _ = sha256_round_trace(M, rounds=R)
    return [(states[R][i] + H0[i]) & MASK for i in range(8)]

def experiment_affine_pairs():
    """
    For R=12: compute the full Hessian of one output bit.
    Identify which input pairs (j,k) are affine (H=0) vs quadratic (H=1).
    """
    print("=" * 80)
    print("WEAKNESS WINDOW: Affine structure at R=12")
    print("=" * 80)

    R = 12
    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]

    def f(msg):
        """Bit 0 of hash register a at R rounds."""
        return sha256_hash_R(msg, R)[0] & 1

    f0 = f(M)

    # Compute FULL Hessian for first 64 input bits (2 words × 32 bits)
    # H[j][k] for j,k in first 64 bits
    n_test = 64  # First 64 input bits (W[0] and W[1])

    hessian = [[0]*n_test for _ in range(n_test)]
    n_nonzero = 0

    for j in range(n_test):
        wj, bj = j // 32, j % 32
        Mj = list(M); Mj[wj] ^= (1 << bj)
        fj = f(Mj)

        for k in range(j+1, n_test):
            wk, bk = k // 32, k % 32
            Mk = list(M); Mk[wk] ^= (1 << bk)
            Mjk = list(M); Mjk[wj] ^= (1 << bj); Mjk[wk] ^= (1 << bk)

            fk = f(Mk)
            fjk = f(Mjk)

            h = f0 ^ fj ^ fk ^ fjk
            hessian[j][k] = h
            hessian[k][j] = h
            if h:
                n_nonzero += 1

    total_pairs = n_test * (n_test - 1) // 2
    p_nonzero = n_nonzero / total_pairs

    print(f"\n  R={R}, first {n_test} input bits:")
    print(f"    Nonzero Hessian entries: {n_nonzero}/{total_pairs} = {p_nonzero:.3f}")
    print(f"    Affine pairs: {total_pairs - n_nonzero}/{total_pairs} = {1-p_nonzero:.3f}")

    # Structure of nonzero entries: are they clustered or scattered?
    # Count per row: how many other bits each input bit interacts with nonlinearly
    row_counts = [sum(hessian[j][k] for k in range(n_test) if k != j) for j in range(n_test)]
    avg_row = sum(row_counts) / n_test
    min_row = min(row_counts)
    max_row = max(row_counts)

    print(f"    Per-bit nonlinear interactions: min={min_row}, max={max_row}, avg={avg_row:.1f}")

    # Are W[0] bits more nonlinear than W[1] bits? (W[0] enters first)
    w0_avg = sum(row_counts[:32]) / 32
    w1_avg = sum(row_counts[32:64]) / 32
    print(f"    W[0] bits avg NL interactions: {w0_avg:.1f}")
    print(f"    W[1] bits avg NL interactions: {w1_avg:.1f}")

    # Compare with R=8 and R=16
    for R_test in [8, 16, 64]:
        def f_test(msg):
            return sha256_hash_R(msg, R_test)[0] & 1

        f0_t = f_test(M)
        nz = 0
        n_pairs = 0

        for j in range(n_test):
            wj, bj = j//32, j%32
            Mj = list(M); Mj[wj] ^= (1<<bj)
            fj = f_test(Mj)

            for k in range(j+1, min(j+20, n_test)):  # Sample nearby pairs
                wk, bk = k//32, k%32
                Mk = list(M); Mk[wk] ^= (1<<bk)
                Mjk = list(M); Mjk[wj] ^= (1<<bj); Mjk[wk] ^= (1<<bk)
                h = f0_t ^ fj ^ f_test(Mk) ^ f_test(Mjk)
                n_pairs += 1
                if h: nz += 1

        print(f"    R={R_test}: P(H=1) ≈ {nz/n_pairs:.3f} (sampled)")


def experiment_is_hessian_stable():
    """
    CRITICAL: Is the Hessian (which pairs are affine) the SAME for all messages?
    If yes → structural property, exploitable.
    If no → state-dependent, not directly exploitable.
    """
    print("\n" + "=" * 80)
    print("STABILITY: Is the Hessian pattern the same for all messages?")
    print("=" * 80)

    R = 12
    n_test = 32  # Just W[0] bits for speed

    # Compute Hessian for 10 different messages
    hessians = []
    for seed in range(10):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]

        def f(msg):
            return sha256_hash_R(msg, R)[0] & 1

        f0 = f(M)
        H = [[0]*n_test for _ in range(n_test)]

        for j in range(n_test):
            Mj = list(M); Mj[0] ^= (1<<j)
            fj = f(Mj)
            for k in range(j+1, n_test):
                Mk = list(M); Mk[0] ^= (1<<k)
                Mjk = list(M); Mjk[0] ^= (1<<j); Mjk[0] ^= (1<<k)
                h = f0 ^ fj ^ f(Mk) ^ f(Mjk)
                H[j][k] = h

        hessians.append(H)

    # Compare pairwise: Hamming distance between Hessian patterns
    print(f"\n  Pairwise Hessian Hamming distances (R={R}, {n_test} bits):")
    dists = []
    for i in range(10):
        for j_idx in range(i+1, 10):
            d = 0
            total = 0
            for j in range(n_test):
                for k in range(j+1, n_test):
                    if hessians[i][j][k] != hessians[j_idx][j][k]:
                        d += 1
                    total += 1
            dists.append(d/total)

    avg_dist = sum(dists)/len(dists)
    print(f"    Average fractional distance: {avg_dist:.4f}")
    print(f"    (0 = identical, 0.5 = random)")

    if avg_dist < 0.1:
        print(f"    *** STABLE: Hessian pattern nearly same for all messages! ***")
        print(f"    → Exploitable structural weakness at R={R}.")
    elif avg_dist < 0.3:
        print(f"    Partially stable: some structure shared across messages.")
    else:
        print(f"    Unstable: Hessian depends on message. State-dependent.")


if __name__ == "__main__":
    experiment_affine_pairs()
    experiment_is_hessian_stable()
