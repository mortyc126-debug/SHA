#!/usr/bin/env python3
"""
Carry Basis — Stage 2.1: The 23-dimensional manifold
======================================================

Stage 2.0 found rank(carry-web) = 23 over GF(2).
This means 23 carry bits form a BASIS, and 41 others are determined.

Stage 2.1:
  B1: Find the 23 basis rounds via GF(2) Gaussian elimination
  B2: Express all 64 carries as functions of the 23 basis carries
  B3: Test: does the basis predict carry[63]? (bridge through algebra)
  B4: Stability: is the same basis valid across different W[0..15]?
"""

import numpy as np
from collections import defaultdict, Counter
import time
import sys

MASK = 0xFFFFFFFF

H0 = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
]

K = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
    0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
    0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
    0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
    0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
    0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
]

def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK
def sig0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sig1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)
def Sig0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sig1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)

def message_schedule(W16):
    W = list(W16) + [0] * 48
    for i in range(16, 64):
        W[i] = (sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16]) & MASK
    return W

def get_carries(W16):
    W = message_schedule(W16)
    a, b, c, d, e, f, g, h = H0
    carries = []
    for r in range(64):
        raw = h + Sig1(e) + Ch(e, f, g) + K[r] + W[r]
        carries.append(1 if raw >= (1 << 32) else 0)
        T1 = raw & MASK
        T2 = (Sig0(a) + Maj(a, b, c)) & MASK
        h, g, f = g, f, e
        e = (d + T1) & MASK
        d, c, b = c, b, a
        a = (T1 + T2) & MASK
    return carries


def gf2_rref(M):
    """GF(2) reduced row echelon form. Returns (rref, pivot_cols, rank)."""
    A = M.copy().astype(np.int8)
    rows, cols = A.shape
    pivot_cols = []
    rank = 0
    for col in range(cols):
        pivot = None
        for row in range(rank, rows):
            if A[row, col] == 1:
                pivot = row
                break
        if pivot is None:
            continue
        A[[rank, pivot]] = A[[pivot, rank]]
        for row in range(rows):
            if row != rank and A[row, col] == 1:
                A[row] = (A[row] + A[rank]) % 2
        pivot_cols.append(col)
        rank += 1
    return A, pivot_cols, rank


# ============================================================
# B1: Find the 23 basis rounds
# ============================================================

def experiment_B1(N=10000, seed=300):
    """
    Collect N carry profiles, build 64-column matrix, find GF(2) basis.
    The basis rounds = pivot columns of RREF.
    """
    print("=" * 70)
    print("B1: Find the 23 Basis Rounds of Carry-Web")
    print(f"N={N}, seed={seed}")
    print("=" * 70)

    rng = np.random.RandomState(seed)

    # Collect carry profiles
    profiles = np.zeros((N, 64), dtype=np.int8)
    t0 = time.time()
    for i in range(N):
        W16 = [int(rng.randint(0, 1 << 32)) for _ in range(16)]
        c = get_carries(W16)
        profiles[i] = c
        if (i + 1) % 5000 == 0:
            print(f"  {i+1}/{N} ({time.time()-t0:.1f}s)")

    print(f"Collected: {time.time()-t0:.1f}s")

    # GF(2) RREF of the profile matrix (N × 64)
    # Rows = samples, Columns = rounds
    # RREF on M gives pivot COLUMNS = basis rounds
    M = profiles  # N × 64
    rref, pivots, rank = gf2_rref(M)

    print(f"\n--- GF(2) Basis ---")
    print(f"  Rank: {rank}")
    print(f"  Basis rounds (pivot columns of profile^T): {pivots}")
    print(f"  Dependent rounds: {[r for r in range(64) if r not in pivots]}")

    # Classify basis rounds by zone
    zones = {'island1': [], 'gap_early': [], 'gap_deep': [], 'island2': [], 'tail': []}
    for r in pivots:
        if r <= 13: zones['island1'].append(r)
        elif r <= 23: zones['gap_early'].append(r)
        elif r <= 46: zones['gap_deep'].append(r)
        elif r <= 55: zones['island2'].append(r)
        else: zones['tail'].append(r)

    print(f"\n  Basis by zone:")
    for zname, rounds in zones.items():
        print(f"    {zname:15s}: {len(rounds)} rounds: {rounds}")

    # Per-round: P(carry=0)
    p_c0 = np.mean(1 - profiles, axis=0)
    print(f"\n  Basis rounds P(carry=0):")
    for r in pivots:
        k_ratio = K[r] / (1 << 32)
        print(f"    r={r:2d}: P(c=0)={p_c0[r]:.4f}, K/2^32={k_ratio:.3f}")

    # Fixed rounds (always carry=1)
    fixed = [r for r in range(64) if p_c0[r] < 0.0005]
    print(f"\n  Fixed carry=1 rounds (P(c=0)<0.0005): {len(fixed)}")
    print(f"    {fixed}")

    # Non-basis non-fixed rounds: these are linearly dependent
    dependent_nonfixed = [r for r in range(64) if r not in pivots and r not in fixed]
    print(f"  Dependent non-fixed rounds: {len(dependent_nonfixed)}")
    print(f"    {dependent_nonfixed}")

    return profiles, pivots, rank, rref, fixed


# ============================================================
# B2: Express dependent carries as GF(2) functions of basis
# ============================================================

def experiment_B2(profiles, pivots, rank, N_test=5000, seed=301):
    """
    For each non-basis round, find GF(2) linear combination of basis
    rounds that predicts it. Measure prediction accuracy.
    """
    print("\n" + "=" * 70)
    print("B2: Express Carries as Functions of Basis")
    print(f"seed={seed}")
    print("=" * 70)

    N = profiles.shape[0]
    n_basis = len(pivots)

    # Build basis matrix: N × n_basis
    B = profiles[:, pivots]  # N × n_basis

    # For each non-basis round, solve B @ x = profiles[:, r] over GF(2)
    # Using least-squares over GF(2) (majority vote per equation)

    non_basis = [r for r in range(64) if r not in pivots]

    # Split train/test
    rng = np.random.RandomState(seed)
    perm = rng.permutation(N)
    train = perm[:N//2]
    test = perm[N//2:]

    B_train = B[train]
    B_test = B[test]

    results = {}
    print(f"\n  {'r':>3} | {'accuracy':>9} | {'formula (basis rounds XORed)':>40} | {'zone':>10}")
    print(f"  {'-'*3}-+-{'-'*9}-+-{'-'*40}-+-{'-'*10}")

    perfect_rounds = []
    good_rounds = []

    for r in non_basis:
        target_train = profiles[train, r]
        target_test = profiles[test, r]

        # Check if this round is always 1 (fixed)
        if np.std(target_train) < 0.01:
            results[r] = ('fixed', 1.0, [])
            continue

        # GF(2) regression: try all subsets up to size 4 (brute force for small basis)
        # For efficiency: use GF(2) Gaussian elimination
        # Augment: [B_train | target_train]
        aug = np.hstack([B_train, target_train.reshape(-1, 1)])  # N_train × (n_basis+1)
        aug_T = aug.T  # (n_basis+1) × N_train

        rref_aug, pivots_aug, rank_aug = gf2_rref(aug_T)

        # If the target column is a pivot → not in span → imperfect
        # If not a pivot → in span → express as combination
        target_is_pivot = n_basis in pivots_aug

        if not target_is_pivot:
            # Find which basis rounds contribute
            # The target column (index n_basis) is expressed via other pivots
            # Read off from RREF: row where col n_basis has a 1
            contributors = []
            for i, pc in enumerate(pivots_aug):
                if pc < n_basis and rref_aug[i, n_basis] == 1:
                    contributors.append(pivots[pc])

            # Predict on test set
            pred = np.zeros(len(test), dtype=np.int8)
            for c_r in contributors:
                c_idx = pivots.index(c_r)
                pred = (pred + B_test[:, c_idx]) % 2

            accuracy = np.mean(pred == target_test)
        else:
            contributors = []
            accuracy = np.mean(target_test)  # predict majority class

        zone = "island1" if r <= 13 else ("gap" if r <= 46 else ("island2" if r <= 55 else "tail"))
        formula = ' ⊕ '.join([f'c[{c}]' for c in contributors]) if contributors else '(no GF2 formula)'

        if accuracy > 0.99:
            perfect_rounds.append(r)
        elif accuracy > 0.90:
            good_rounds.append(r)

        marker = " ★★★" if accuracy > 0.99 else (" ★★" if accuracy > 0.95 else (" ★" if accuracy > 0.90 else ""))
        print(f"  {r:3d} | {accuracy:9.4f} | {formula:>40s} | {zone:>10}{marker}")

        results[r] = (formula, accuracy, contributors)

    print(f"\n  Perfect prediction (>99%): {len(perfect_rounds)} rounds: {perfect_rounds}")
    print(f"  Good prediction (>90%):    {len(good_rounds)} rounds: {good_rounds}")
    print(f"  Imperfect (<90%):          {64 - len(pivots) - len(perfect_rounds) - len(good_rounds)} rounds")

    return results


# ============================================================
# B3: Does the basis predict carry[63]?
# ============================================================

def experiment_B3(profiles, pivots, N_test=3000, seed=302):
    """
    KEY TEST: Can carry[63] be predicted from the 23 basis carries?
    If yes → algebraic bridge through the carry-web manifold.
    """
    print("\n" + "=" * 70)
    print("B3: Can Basis Predict carry[63]? (The Bridge Test)")
    print(f"seed={seed}")
    print("=" * 70)

    N = profiles.shape[0]

    # Is carry[63] in the basis?
    if 63 in pivots:
        print(f"\n  carry[63] IS a basis round — it's independent, not predictable from others.")
        print(f"  This means: carry[63] adds 1 bit of information beyond the other basis carries.")
        print(f"  The bridge does NOT exist through GF(2) linear algebra alone.")
        return 'basis_member'

    # carry[63] is dependent — it CAN be expressed via basis
    print(f"\n  carry[63] is a DEPENDENT round — it's in the span of basis carries!")

    # Find the GF(2) formula
    B = profiles[:, pivots]
    target = profiles[:, 63]

    if np.std(target) < 0.01:
        print(f"  carry[63] is FIXED (always {int(np.mean(target))})")
        return 'fixed'

    # GF(2) solve: find which basis rounds XOR to give carry[63]
    # Augmented system
    aug_T = np.vstack([B.T, target.reshape(1, -1)])
    rref_aug, pivots_aug, rank_aug = gf2_rref(aug_T)

    n_basis = len(pivots)
    if n_basis in pivots_aug:
        print(f"  carry[63] NOT in GF(2) span — this contradicts rank analysis!")
        print(f"  Possible: nonlinear dependence only.")
        return 'nonlinear'

    # Read formula
    contributors = []
    for i, pc in enumerate(pivots_aug):
        if pc < n_basis and rref_aug[i, n_basis] == 1:
            contributors.append(pivots[pc])

    formula = ' ⊕ '.join([f'c[{r}]' for r in contributors])
    print(f"\n  FORMULA: carry[63] = {formula}")
    print(f"  Number of terms: {len(contributors)}")

    # Verify on data
    pred = np.zeros(N, dtype=np.int8)
    for r in contributors:
        idx = pivots.index(r)
        pred = (pred + B[:, idx]) % 2

    accuracy = np.mean(pred == target)
    n_correct = np.sum(pred == target)
    print(f"\n  Verification: {n_correct}/{N} correct ({accuracy*100:.2f}%)")

    if accuracy > 0.99:
        print(f"\n  ★★★★★ ALGEBRAIC BRIDGE FOUND ★★★★★")
        print(f"  carry[63] is a GF(2) linear function of {len(contributors)} basis carries!")
        print(f"  This means: knowing carries at rounds {contributors}")
        print(f"  DETERMINES carry[63] with {accuracy*100:.1f}% accuracy.")
        print(f"  The carry-web manifold provides an algebraic path to carry[63].")
    elif accuracy > 0.90:
        print(f"\n  ★★★ Strong but imperfect bridge (nonlinear residual)")
    else:
        print(f"\n  Weak bridge — GF(2) formula insufficient")

    # Which zones do the contributors come from?
    z_counts = {'island1': 0, 'gap': 0, 'island2': 0, 'tail': 0}
    for r in contributors:
        if r <= 13: z_counts['island1'] += 1
        elif r <= 46: z_counts['gap'] += 1
        elif r <= 55: z_counts['island2'] += 1
        else: z_counts['tail'] += 1

    print(f"\n  Contributors by zone: {z_counts}")
    if z_counts['gap'] > 0:
        print(f"  ★ GAP ROUNDS CONTRIBUTE — the formula bridges through the gap!")
    if z_counts['island1'] > 0 and z_counts['island2'] > 0:
        print(f"  ★ BOTH ISLANDS CONTRIBUTE — the formula spans the full range!")

    return contributors


# ============================================================
# B4: Basis stability across different W[0..15]
# ============================================================

def experiment_B4(N_per_config=3000, n_configs=10, seed=303):
    """
    Is the basis stable? Run GF(2) RREF on different random W[0..15] sets.
    If same 23 rounds → basis is UNIVERSAL (structural property of SHA-256).
    If different → basis is W-dependent (less useful).
    """
    print("\n" + "=" * 70)
    print("B4: Basis Stability — Is the same basis valid for all W?")
    print(f"N_per_config={N_per_config}, n_configs={n_configs}, seed={seed}")
    print("=" * 70)

    rng = np.random.RandomState(seed)
    all_bases = []
    all_ranks = []

    t0 = time.time()
    for config in range(n_configs):
        profiles = np.zeros((N_per_config, 64), dtype=np.int8)
        for i in range(N_per_config):
            W16 = [int(rng.randint(0, 1 << 32)) for _ in range(16)]
            profiles[i] = get_carries(W16)

        M = profiles.T
        _, pivots, rank = gf2_rref(M)
        all_bases.append(set(pivots))
        all_ranks.append(rank)

    print(f"Computed: {time.time()-t0:.1f}s")

    # Stability analysis
    print(f"\n  Rank distribution: {Counter(all_ranks)}")

    # Intersection of all bases
    common = all_bases[0]
    for b in all_bases[1:]:
        common = common & b
    union = all_bases[0]
    for b in all_bases[1:]:
        union = union | b

    print(f"\n  Common basis rounds (in ALL configs): {sorted(common)} ({len(common)} rounds)")
    print(f"  Union (in ANY config): {sorted(union)} ({len(union)} rounds)")
    print(f"  Stability ratio: {len(common)}/{len(union)} = {len(common)/max(len(union),1)*100:.0f}%")

    # Per-round: how often is it in the basis?
    round_freq = Counter()
    for b in all_bases:
        for r in b:
            round_freq[r] += 1

    print(f"\n  Per-round basis membership ({n_configs} configs):")
    stable_rounds = []
    for r in sorted(round_freq.keys()):
        freq = round_freq[r]
        bar = '#' * (freq * 5)
        stability = "STABLE" if freq == n_configs else ("semi" if freq > n_configs * 0.5 else "unstable")
        if freq == n_configs:
            stable_rounds.append(r)
        print(f"    r={r:2d}: {freq:2d}/{n_configs} {bar} {stability}")

    print(f"\n  Universally stable basis rounds: {stable_rounds}")

    # Jaccard similarity between bases
    jaccards = []
    for i in range(n_configs):
        for j in range(i+1, n_configs):
            inter = len(all_bases[i] & all_bases[j])
            union_ij = len(all_bases[i] | all_bases[j])
            jaccards.append(inter / max(union_ij, 1))
    print(f"  Mean Jaccard similarity: {np.mean(jaccards):.3f}")

    return all_bases, stable_rounds


# ============================================================
# MAIN
# ============================================================

if __name__ == '__main__':
    print("Carry Basis — Stage 2.1: The 23-Dimensional Manifold")
    print(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print()

    N1 = 8000
    if len(sys.argv) > 1 and sys.argv[1] == '--fast':
        N1 = 4000

    t_start = time.time()

    profiles, pivots, rank, rref, fixed = experiment_B1(N=N1)
    results = experiment_B2(profiles, pivots, rank)
    bridge = experiment_B3(profiles, pivots)
    all_bases, stable = experiment_B4(N_per_config=3000, n_configs=10)

    total = time.time() - t_start
    print(f"\n{'='*70}")
    print(f"TOTAL TIME: {total:.1f}s")
    print(f"{'='*70}")

    print(f"""
{'='*70}
SYNTHESIS: Carry Basis (Stage 2.1)
{'='*70}

THE 23-DIMENSIONAL MANIFOLD:
  Basis rounds: {pivots}
  Rank: {rank}
  Fixed rounds: {len(fixed)}

CARRY[63] BRIDGE:
  {'FOUND — carry[63] expressible via basis' if isinstance(bridge, list) else bridge}

BASIS STABILITY:
  Universally stable rounds: {stable}
  Jaccard similarity: {np.mean([len(a&b)/max(len(a|b),1) for a in all_bases for b in all_bases if a!=b]):.3f}

IMPLICATION:
  If carry[63] = XOR of basis carries → the manifold provides
  an algebraic formula connecting early rounds to the final carry.
  This is NOT through state (chaotic zone blocks that).
  This is NOT through schedule alone (W3 showed state coupling).
  This IS through the GF(2) structure of the carry-web itself.
  A genuinely new algebraic object.
""")
