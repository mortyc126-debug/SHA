#!/usr/bin/env python3
"""
CRAZY-31: Linear trail analysis of De17/De18/De19 barrier bits.

Question: Do linear combinations of De_t output bits vanish for ALL W[14]?
If so, those are FREE equations — reducing the effective barrier rank.

Method:
1. Fix base message, sample many W[14] values, compute De_t for each.
2. Build GF(2) matrix of output bits, find null space.
3. Check pairwise and triple XOR relations.
4. Cross-validate across multiple base messages.
"""

import random
import numpy as np
from itertools import combinations

MASK = 0xFFFFFFFF
K = [0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
     0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
     0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa]
H0 = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&MASK
def Sig0(x): return rotr(x,2)^rotr(x,13)^rotr(x,22)
def Sig1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def sig0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def sig1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def Ch(e,f,g): return ((e&f)^(~e&g))&MASK
def Maj(a,b,c): return ((a&b)^(a&c)^(b&c))&MASK
def add32(*a):
    s=0
    for x in a: s=(s+x)&MASK
    return s

def sha_round(st,w,k):
    a,b,c,d,e,f,g,h=st
    T1=add32(h,Sig1(e),Ch(e,f,g),k,w)
    T2=add32(Sig0(a),Maj(a,b,c))
    return [add32(T1,T2),a,b,c,add32(d,T1),e,f,g]

def expand_W(W16, num=22):
    W=list(W16)
    for i in range(16, num):
        W.append(add32(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    return W

def compute_state_up_to(W16, num_rounds):
    """Compute SHA-256 state after num_rounds rounds."""
    W = expand_W(W16, num=max(num_rounds, 16))
    st = list(H0)
    for i in range(num_rounds):
        st = sha_round(st, W[i], K[i])
    return st

def compute_De(base_msg, w14_val, target_round):
    """
    Compute De_t = e_t(W[14]=w14_val) - e_t(W[14]=0) mod 2^32
    where e_t is the 'e' register after round t.

    In SHA-256 state [a,b,c,d,e,f,g,h], e is index 4.
    """
    msg0 = list(base_msg)
    msg0[14] = 0

    msg1 = list(base_msg)
    msg1[14] = w14_val

    st0 = compute_state_up_to(msg0, target_round)
    st1 = compute_state_up_to(msg1, target_round)

    # De = difference in e register
    e0 = st0[4]
    e1 = st1[4]
    De = (e1 - e0) & MASK
    return De

def gf2_rank(matrix):
    """Compute rank of binary matrix over GF(2). matrix is numpy uint8 array."""
    M = matrix.copy()
    rows, cols = M.shape
    rank = 0
    for col in range(cols):
        # Find pivot
        pivot = None
        for row in range(rank, rows):
            if M[row, col] == 1:
                pivot = row
                break
        if pivot is None:
            continue
        # Swap
        M[[rank, pivot]] = M[[pivot, rank]]
        # Eliminate
        for row in range(rows):
            if row != rank and M[row, col] == 1:
                M[row] = M[row] ^ M[rank]
        rank += 1
    return rank

def gf2_nullspace(matrix):
    """Find null space of binary matrix over GF(2). Returns list of null vectors."""
    M = matrix.copy()
    rows, cols = M.shape
    # Augment with identity
    aug = np.zeros((rows, cols + rows), dtype=np.uint8)
    aug[:, :cols] = M
    for i in range(rows):
        aug[i, cols + i] = 1

    # Row reduce
    pivot_col = [None] * rows
    r = 0
    for col in range(cols):
        pivot = None
        for row in range(r, rows):
            if aug[row, col] == 1:
                pivot = row
                break
        if pivot is None:
            continue
        aug[[r, pivot]] = aug[[pivot, r]]
        pivot_col[r] = col
        for row in range(rows):
            if row != r and aug[row, col] == 1:
                aug[row] = aug[row] ^ aug[r]
        r += 1

    # Null vectors: rows from r onwards of the identity part give null space
    null_vecs = []
    for i in range(r, rows):
        vec = aug[i, cols:]
        # Verify it's in null space
        if np.all((matrix.T @ vec.reshape(-1, 1)) % 2 == 0):
            null_vecs.append(vec)
    return null_vecs

def bits_to_vec(x, nbits=32):
    """Extract bits of x into a vector."""
    return np.array([(x >> k) & 1 for k in range(nbits)], dtype=np.uint8)

def analyze_barrier(base_msg, target_round, n_samples=10000, rng=None):
    """Analyze linear structure of De_t bits."""
    if rng is None:
        rng = random.Random(0xC031)

    round_name = f"De{target_round}"

    # Sample W[14] values and compute De
    w14_samples = [rng.randint(0, MASK) for _ in range(n_samples)]
    de_values = []
    for w14 in w14_samples:
        de = compute_De(base_msg, w14, target_round)
        de_values.append(de)

    # Build matrix: n_samples x 32
    M = np.zeros((n_samples, 32), dtype=np.uint8)
    for i, de in enumerate(de_values):
        M[i] = bits_to_vec(de)

    # 1. Rank
    rank = gf2_rank(M)
    nullity = 32 - rank

    print(f"\n{'='*70}")
    print(f"  {round_name} ANALYSIS (target round {target_round})")
    print(f"{'='*70}")
    print(f"  Samples: {n_samples}")
    print(f"  GF(2) rank of {n_samples}x32 matrix: {rank}")
    print(f"  Nullity (32 - rank): {nullity}")

    if nullity > 0:
        print(f"  *** FOUND {nullity} LINEAR RELATION(S)! ***")
    else:
        print(f"  All 32 bits are GF(2)-independent over samples.")

    # 2. Find null space of M^T (32 x n_samples)
    # We want vectors v in {0,1}^32 such that M @ v = 0 (mod 2) for all rows
    # Equivalently: v^T M^T = 0, i.e., null space of M (treating rows as constraints on v)
    # Actually: we want v such that for each sample i, sum_k v_k * M[i,k] = 0 mod 2
    # This means M v = 0 (mod 2), where M is n_samples x 32
    # So we need null space of M (viewed as linear map from R^32 to R^n_samples over GF(2))

    # Compute rank of M^T (32 x n_samples) - same as rank of M
    # Null space of M: vectors v in GF(2)^32 with Mv = 0
    # Use row reduction on M^T to find these

    # Simpler: row-reduce M itself
    M2 = M.copy()
    rows, cols = M2.shape
    pivots = []
    r = 0
    for col in range(cols):
        pivot = None
        for row in range(r, rows):
            if M2[row, col] == 1:
                pivot = row
                break
        if pivot is None:
            continue
        M2[[r, pivot]] = M2[[pivot, r]]
        pivots.append(col)
        for row in range(rows):
            if row != r and M2[row, col] == 1:
                M2[row] ^= M2[r]
        r += 1

    free_cols = [c for c in range(32) if c not in pivots]

    if free_cols:
        print(f"\n  Null space basis vectors (bit masks):")
        null_vectors = []
        for fc in free_cols:
            # Build null vector: set free variable to 1, solve for pivot variables
            v = np.zeros(32, dtype=np.uint8)
            v[fc] = 1
            # From reduced form, for each pivot row i with pivot column pivots[i]:
            # M2[i, pivots[i]] = 1, and M2[i, fc] gives the coefficient
            # Actually need to re-reduce properly
            # Use the RREF form
            pass

        # Better approach: directly find null space
        # Row reduce the transpose (32 x n_samples) with identity augmentation
        Mt = M.T.copy()  # 32 x n_samples
        aug = np.zeros((32, n_samples + 32), dtype=np.uint8)
        aug[:, :n_samples] = Mt
        for i in range(32):
            aug[i, n_samples + i] = 1

        r = 0
        for col in range(n_samples):
            pivot = None
            for row in range(r, 32):
                if aug[row, col] == 1:
                    pivot = row
                    break
            if pivot is None:
                continue
            aug[[r, pivot]] = aug[[pivot, r]]
            for row in range(32):
                if row != r and aug[row, col] == 1:
                    aug[row] ^= aug[r]
            r += 1
            if r >= 32:
                break

        # Rows r..31 of the augmented identity part are null space of M^T
        # But we want null space of M, not M^T
        # Null space of M: {v : Mv = 0}
        # Let's just use the direct method

    # Direct null space computation for M (n x 32 over GF(2))
    # Row reduce M, find free variables, back-substitute
    M3 = M.copy()
    pivot_row = {}
    r = 0
    for col in range(32):
        pivot = None
        for row in range(r, n_samples):
            if M3[row, col] == 1:
                pivot = row
                break
        if pivot is None:
            continue
        M3[[r, pivot]] = M3[[pivot, r]]
        pivot_row[col] = r
        for row in range(n_samples):
            if row != r and M3[row, col] == 1:
                M3[row] ^= M3[r]
        r += 1

    free_vars = sorted(set(range(32)) - set(pivot_row.keys()))
    null_vecs = []

    for fv in free_vars:
        v = np.zeros(32, dtype=np.uint8)
        v[fv] = 1
        for pcol, prow in pivot_row.items():
            if M3[prow, fv] == 1:
                v[pcol] = 1
        null_vecs.append(v)

    if null_vecs:
        print(f"\n  Null space basis ({len(null_vecs)} vectors):")
        for idx, v in enumerate(null_vecs):
            bits_set = [k for k in range(32) if v[k]]
            mask_val = sum(1 << k for k in bits_set)
            print(f"    v{idx}: bits {bits_set} (mask 0x{mask_val:08x})")

            # Verify: check that XOR of these bits is always 0
            violations = 0
            for de in de_values:
                parity = 0
                for b in bits_set:
                    parity ^= (de >> b) & 1
                if parity != 0:
                    violations += 1
            print(f"         Verification: {violations}/{n_samples} violations")

    # 3. Pairwise XOR analysis
    print(f"\n  Pairwise XOR analysis (always 0 = identical, always 1 = complement):")
    pair_found = False
    for j in range(32):
        for k in range(j+1, 32):
            xor_sum = 0
            for de in de_values:
                xor_sum += ((de >> j) ^ (de >> k)) & 1
            if xor_sum == 0:
                print(f"    Bits {j} and {k}: ALWAYS EQUAL")
                pair_found = True
            elif xor_sum == n_samples:
                print(f"    Bits {j} and {k}: ALWAYS COMPLEMENT")
                pair_found = True
    if not pair_found:
        print(f"    No always-equal or always-complement pairs found.")

    # 4. Triple XOR (only if rank < 32 or for completeness, check a subset)
    triple_count = 0
    if nullity > 0 or rank <= 30:
        print(f"\n  Triple XOR analysis:")
        for j in range(32):
            for k in range(j+1, 32):
                for l in range(k+1, 32):
                    xor_sum = 0
                    for de in de_values:
                        xor_sum += ((de >> j) ^ (de >> k) ^ (de >> l)) & 1
                    if xor_sum == 0 or xor_sum == n_samples:
                        val = 0 if xor_sum == 0 else 1
                        print(f"    Bits {j},{k},{l}: ALWAYS {val}")
                        triple_count += 1
        if triple_count == 0:
            print(f"    No always-zero/one triples found.")
    else:
        print(f"\n  Triple XOR: skipped (rank=32, no relations expected).")

    # 5. Bit entropy analysis (how biased is each bit?)
    print(f"\n  Per-bit bias (fraction of 1s):")
    biases = []
    for k in range(32):
        count1 = sum(((de >> k) & 1) for de in de_values)
        bias = count1 / n_samples
        biases.append(bias)

    min_bias = min(biases)
    max_bias = max(biases)
    most_biased_bit = min(range(32), key=lambda k: abs(biases[k] - 0.5))
    least_uniform = max(range(32), key=lambda k: abs(biases[k] - 0.5))

    print(f"    Most uniform bit: bit {most_biased_bit} (bias={biases[most_biased_bit]:.4f})")
    print(f"    Most biased bit:  bit {least_uniform} (bias={biases[least_uniform]:.4f})")
    print(f"    Bias range: [{min_bias:.4f}, {max_bias:.4f}]")

    # Summary stats
    always_zero_bits = [k for k in range(32) if biases[k] == 0.0]
    always_one_bits = [k for k in range(32) if biases[k] == 1.0]
    if always_zero_bits:
        print(f"    ALWAYS-ZERO bits: {always_zero_bits}")
    if always_one_bits:
        print(f"    ALWAYS-ONE bits: {always_one_bits}")

    return rank, nullity, null_vecs, de_values, biases


def cross_validate(target_round, n_messages=5, n_samples=5000):
    """Check if null-space vectors are universal across different messages."""
    print(f"\n{'='*70}")
    print(f"  CROSS-VALIDATION: De{target_round} across {n_messages} base messages")
    print(f"{'='*70}")

    all_null_vecs = []
    all_ranks = []

    seeds = [0xC031, 0xDEAD, 0xBEEF, 0x1337, 0xCAFE]

    for msg_idx in range(n_messages):
        seed = seeds[msg_idx]
        rng_msg = random.Random(seed)
        base_msg = [rng_msg.randint(0, MASK) for _ in range(16)]
        base_msg[14] = 0  # Will be varied

        rng_sample = random.Random(seed * 7 + 1)

        rank, nullity, nvecs, _, _ = analyze_barrier(
            base_msg, target_round, n_samples=n_samples, rng=rng_sample
        )
        all_ranks.append(rank)
        all_null_vecs.append(nvecs)

    print(f"\n{'='*70}")
    print(f"  CROSS-VALIDATION SUMMARY: De{target_round}")
    print(f"{'='*70}")
    print(f"  Ranks across messages: {all_ranks}")

    if all(r == 32 for r in all_ranks):
        print(f"  All messages show rank 32 — no universal linear relations.")
    else:
        # Find common null vectors
        if all_null_vecs[0]:
            print(f"  Checking for UNIVERSAL null vectors...")
            for v_idx, v in enumerate(all_null_vecs[0]):
                bits = tuple(k for k in range(32) if v[k])
                mask = sum(1 << k for k in bits)
                universal = True
                for msg_idx in range(1, n_messages):
                    # Check if this vector is in the null space of other messages
                    found = False
                    for v2 in all_null_vecs[msg_idx]:
                        if np.array_equal(v, v2):
                            found = True
                            break
                    if not found:
                        universal = False
                        break
                status = "UNIVERSAL" if universal else "message-specific"
                print(f"    v{v_idx} bits={bits} mask=0x{mask:08x}: {status}")

    return all_ranks, all_null_vecs


def main():
    print("=" * 70)
    print("  CRAZY-31: Linear Trail Analysis of SHA-256 Barrier Bits")
    print("  Finding free equations in De17/De18/De19")
    print("=" * 70)

    # Base message from seed 0xC031
    rng_base = random.Random(0xC031)
    base_msg = [rng_base.randint(0, MASK) for _ in range(16)]
    base_msg[14] = 0  # Will be varied

    print(f"\n  Base message (W[14] set to 0 as reference):")
    for i in range(0, 16, 4):
        print(f"    W[{i:2d}..{i+3:2d}]: " + " ".join(f"0x{base_msg[j]:08x}" for j in range(i, i+4)))

    # ---- De17 Analysis ----
    rng17 = random.Random(0xC031_17)
    rank17, null17, nvecs17, de17_vals, bias17 = analyze_barrier(
        base_msg, 17, n_samples=10000, rng=rng17
    )

    # ---- De18 Analysis ----
    rng18 = random.Random(0xC031_18)
    rank18, null18, nvecs18, de18_vals, bias18 = analyze_barrier(
        base_msg, 18, n_samples=10000, rng=rng18
    )

    # ---- De19 Analysis ----
    rng19 = random.Random(0xC031_19)
    rank19, null19, nvecs19, de19_vals, bias19 = analyze_barrier(
        base_msg, 19, n_samples=10000, rng=rng19
    )

    # ---- Comparison ----
    print(f"\n{'='*70}")
    print(f"  COMPARISON: Ranks and Nullities")
    print(f"{'='*70}")
    print(f"  De17: rank={rank17}, nullity={null17}")
    print(f"  De18: rank={rank18}, nullity={null18}")
    print(f"  De19: rank={rank19}, nullity={null19}")

    if null17 > 0 or null18 > 0 or null19 > 0:
        print(f"\n  *** FREE EQUATIONS FOUND! ***")
        total_free = null17 + null18 + null19
        print(f"  Total free equations across De17-De19: {total_free}")
        print(f"  Effective barrier reduction: 96 bits -> {96 - total_free} bits")
    else:
        print(f"\n  No linear relations found in any barrier.")
        print(f"  All barriers have full rank 32 — the 96-bit barrier is solid.")

    # ---- Cross-validation ----
    print(f"\n\n{'#'*70}")
    print(f"  CROSS-VALIDATION PHASE")
    print(f"{'#'*70}")

    for t in [17, 18, 19]:
        cross_validate(t, n_messages=5, n_samples=5000)

    # ---- Final verdict ----
    print(f"\n\n{'#'*70}")
    print(f"  FINAL VERDICT")
    print(f"{'#'*70}")
    print(f"""
  De17 rank: {rank17}/32  (nullity {null17})
  De18 rank: {rank18}/32  (nullity {null18})
  De19 rank: {rank19}/32  (nullity {null19})

  If all ranks are 32:
    -> Each barrier round contributes a FULL 32-bit constraint.
    -> No free linear equations exist.
    -> The coupled 96-bit barrier (De17+De18+De19) cannot be reduced
       by linear algebra over GF(2).
    -> The rank-128 coupled system from CRAZY-22.4 is NOT reducible
       via output-bit linear dependencies.

  If any rank < 32:
    -> That barrier has {32 - min(rank17,rank18,rank19)} free equation(s).
    -> These reduce the effective search space.
    -> Cross-validation determines if this is universal or message-specific.
""")


if __name__ == "__main__":
    main()
