#!/usr/bin/env python3
"""
CRAZY-9/10: Linearized SHA-256 Kernel & Carry Perturbation
===========================================================

COMPLETELY NEW ATTACK ANGLE.

Previous 18 experiments all worked in DIFFERENTIAL space (Wang chain, neutral bits).
This experiment works in ALGEBRAIC space:

1. Replace ALL modular additions (+) with XOR (⊕) → linearized SHA-256
2. This gives a GF(2)-linear map L: {0,1}^512 → {0,1}^256
3. Compute kernel(L) — the set of message differences that give zero output
   in the LINEAR approximation
4. For kernel elements ΔM, compute REAL SHA-256:
   H(M) ⊕ H(M ⊕ ΔM) should be 0 in linear world, but carries perturb it
5. Measure the "carry perturbation": HW(H(M) ⊕ H(M ⊕ ΔM))
6. If carry perturbation is LOW (say 30-50 bits instead of 128),
   we have a NEW attack vector: start from linear solutions, correct with GA

This is genuinely novel — nobody has systematically studied the kernel
of linearized SHA-256 as a collision starting point.
"""

import random
import time

MASK = 0xFFFFFFFF

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
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
]

H0 = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
]

def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK
def Sig0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sig1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def sig0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sig1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)
def Ch(e, f, g): return ((e & f) ^ (~e & g)) & MASK
def Maj(a, b, c): return ((a & b) ^ (a & c) ^ (b & c)) & MASK
def add32(*args):
    s = 0
    for a in args: s = (s + a) & MASK
    return s
def xor32(*args):
    s = 0
    for a in args: s ^= a
    return s & MASK
def hw(x): return bin(x & MASK).count('1')

# ── Standard SHA-256 ──────────────────────────────────────────────

def expand_W(W16):
    W = list(W16)
    for i in range(16, 64):
        W.append(add32(sig1(W[i-2]), W[i-7], sig0(W[i-15]), W[i-16]))
    return W

def sha256_compress(msg, iv):
    W = expand_W(msg)
    state = list(iv)
    for t in range(64):
        a, b, c, d, e, f, g, h = state
        T1 = add32(h, Sig1(e), Ch(e, f, g), K[t], W[t])
        T2 = add32(Sig0(a), Maj(a, b, c))
        state = [add32(T1, T2), a, b, c, add32(d, T1), e, f, g]
    return [add32(state[i], iv[i]) for i in range(8)]

# ── Linearized SHA-256 (all + → ⊕) ──────────────────────────────

def expand_W_linear(W16):
    """Message schedule with + replaced by ⊕"""
    W = list(W16)
    for i in range(16, 64):
        W.append(sig1(W[i-2]) ^ W[i-7] ^ sig0(W[i-15]) ^ W[i-16])
    return W

def sha256_compress_linear(msg, iv):
    """SHA-256 compression with all + replaced by ⊕"""
    W = expand_W_linear(msg)
    state = list(iv)
    for t in range(64):
        a, b, c, d, e, f, g, h = state
        T1 = xor32(h, Sig1(e), Ch(e, f, g), K[t], W[t])
        T2 = xor32(Sig0(a), Maj(a, b, c))
        state = [T1 ^ T2, a, b, c, d ^ T1, e, f, g]
    return [state[i] ^ iv[i] for i in range(8)]


# ── GF(2) Matrix for linearized SHA-256 ─────────────────────────

def build_linearized_matrix(iv, num_rounds=64):
    """Build the GF(2) matrix of linearized SHA-256.

    The matrix maps 512 message bits → 256 hash bits.
    Column i = linearized_SHA256(e_i) where e_i has only bit i set.

    Since the linearized system is NOT truly linear (Ch, Maj are nonlinear),
    we use the DIFFERENTIAL: L(M) = lin_SHA256(M) ⊕ lin_SHA256(0)
    This gives the linear part of the function.
    """
    # Compute lin_SHA256(0) as the "constant" term
    zero_msg = [0] * 16

    # For the linearized function, we need to be careful:
    # Ch(e,f,g) = ef ⊕ ~eg = ef ⊕ e'g (where e' = ~e)
    # This is BILINEAR, not linear. So the linearization depends on
    # the operating point.
    #
    # True linearization: take the Jacobian at a specific point.
    # J[i,j] = ∂(output bit i) / ∂(input bit j)
    #
    # We compute this numerically: flip each input bit, observe output change.

    # Base output (no perturbation)
    base_output = sha256_compress_linear(zero_msg, iv)
    base_hash = 0
    for i, val in enumerate(base_output):
        base_hash |= (val << (i * 32))

    # Build matrix: 256 rows × 512 columns
    # Each column = effect of flipping one message bit
    matrix = []  # list of 256-bit integers (columns)

    for word in range(16):
        for bit in range(32):
            msg = [0] * 16
            msg[word] = 1 << bit

            output = sha256_compress_linear(msg, iv)
            out_hash = 0
            for i, val in enumerate(output):
                out_hash |= (val << (i * 32))

            diff = base_hash ^ out_hash  # 256-bit difference
            matrix.append(diff)

    return matrix


def build_real_jacobian(base_msg, iv, num_rounds=64):
    """Build the Jacobian of REAL SHA-256 at a specific operating point.
    J[i,j] = ∂H[i]/∂M[j] computed numerically (flip each bit, observe change)."""

    base_output = sha256_compress(base_msg, iv)
    base_hash = 0
    for i, val in enumerate(base_output):
        base_hash |= (val << (i * 32))

    matrix = []
    for word in range(16):
        for bit in range(32):
            msg = list(base_msg)
            msg[word] ^= (1 << bit)

            output = sha256_compress(msg, iv)
            out_hash = 0
            for i, val in enumerate(output):
                out_hash |= (val << (i * 32))

            diff = base_hash ^ out_hash
            matrix.append(diff)

    return matrix


def gf2_kernel(matrix, n_cols=512, n_rows=256):
    """Compute kernel of GF(2) matrix (list of n_cols 256-bit column vectors).
    Returns list of 512-bit kernel vectors."""

    # Transpose: work with rows for Gaussian elimination
    # Original: n_cols column vectors, each n_rows bits
    # We want to find v ∈ {0,1}^n_cols such that sum(v[i] * col[i]) = 0

    # Augment each column with its index for tracking
    # Do elimination on the TRANSPOSE matrix (n_rows × n_cols)

    # Build row-based representation of transpose
    rows = []
    for r in range(n_rows):
        row = 0
        for c in range(n_cols):
            if (matrix[c] >> r) & 1:
                row |= (1 << c)
        rows.append(row)

    # Gaussian elimination (reduce rows)
    pivot_cols = []
    for r in range(n_rows):
        # Find pivot
        pivot = -1
        for c in range(n_cols):
            if c in pivot_cols:
                continue
            if (rows[r] >> c) & 1:
                pivot = c
                break
        if pivot == -1:
            continue  # zero row
        pivot_cols.append(pivot)

        # Eliminate other rows
        for r2 in range(n_rows):
            if r2 != r and (rows[r2] >> pivot) & 1:
                rows[r2] ^= rows[r]

    # Free variables = columns NOT in pivot_cols
    free_cols = [c for c in range(n_cols) if c not in pivot_cols]
    rank = len(pivot_cols)

    # Build kernel basis
    # For each free variable, set it to 1 and solve for pivot variables
    kernel_basis = []
    for fc in free_cols:
        vec = 1 << fc  # set the free variable
        # For each pivot, check if it depends on this free variable
        for i, pc in enumerate(pivot_cols):
            if (rows[i] >> fc) & 1:
                vec |= (1 << pc)
        kernel_basis.append(vec)

    return kernel_basis, rank


def vec_to_msg(vec):
    """Convert 512-bit vector to 16 message words."""
    msg = []
    for w in range(16):
        val = 0
        for b in range(32):
            if (vec >> (w * 32 + b)) & 1:
                val |= (1 << b)
        msg.append(val)
    return msg


def main():
    print("=" * 72)
    print("CRAZY-9/10: Linearized SHA-256 & Carry Perturbation")
    print("=" * 72)
    print()

    random.seed(0x11AE)
    t_start = time.time()
    iv = list(H0)

    # ── Phase 1: Build linearized matrix ──────────────────────────

    print("Phase 1: Building GF(2) matrix of linearized SHA-256")
    print("-" * 72)

    lin_matrix = build_linearized_matrix(iv)
    print(f"  Matrix: 256 × 512 over GF(2)")
    print(f"  Built in {time.time() - t_start:.2f}s")

    # ── Phase 2: Compute kernel ───────────────────────────────────

    print()
    print("Phase 2: Computing kernel of linearized map")
    print("-" * 72)

    kernel_basis, rank = gf2_kernel(lin_matrix)
    kernel_dim = len(kernel_basis)

    print(f"  Rank: {rank}")
    print(f"  Kernel dimension: {kernel_dim}")
    print(f"  Expected: rank ≈ 256, kernel ≈ 256")
    print()

    # ── Phase 3: Test kernel elements against REAL SHA-256 ────────

    print("Phase 3: Carry perturbation — kernel elements vs real SHA-256")
    print("-" * 72)
    print()

    NUM_TESTS = 200
    NUM_BASE_MSGS = 10

    # Test random kernel elements
    perturbation_hws = []
    min_hw = 256
    best_delta = None
    best_msg = None

    for msg_idx in range(NUM_BASE_MSGS):
        base_msg = [random.getrandbits(32) for _ in range(16)]
        base_hash = sha256_compress(base_msg, iv)

        for trial in range(NUM_TESTS // NUM_BASE_MSGS):
            # Random kernel element = random XOR combination of basis vectors
            delta_vec = 0
            for kb in kernel_basis:
                if random.random() < 0.5:
                    delta_vec ^= kb

            if delta_vec == 0:
                continue

            delta_msg = vec_to_msg(delta_vec)

            # Apply perturbation
            perturbed_msg = [(base_msg[i] ^ delta_msg[i]) & MASK for i in range(16)]

            # Compute real hash
            perturbed_hash = sha256_compress(perturbed_msg, iv)

            # Measure carry perturbation
            hash_diff = 0
            for i in range(8):
                hash_diff += hw(base_hash[i] ^ perturbed_hash[i])

            perturbation_hws.append(hash_diff)

            if hash_diff < min_hw:
                min_hw = hash_diff
                best_delta = delta_msg
                best_msg = base_msg

    if perturbation_hws:
        avg_hw = sum(perturbation_hws) / len(perturbation_hws)
        min_hw_all = min(perturbation_hws)
        max_hw_all = max(perturbation_hws)

        # Distribution
        buckets = {}
        for h in perturbation_hws:
            b = (h // 10) * 10
            buckets[b] = buckets.get(b, 0) + 1

        print(f"  Tested {len(perturbation_hws)} random kernel elements")
        print(f"  Average carry perturbation: {avg_hw:.1f} / 256 bits")
        print(f"  Min: {min_hw_all}, Max: {max_hw_all}")
        print(f"  Expected for random: 128 / 256")
        print()
        print("  Distribution:")
        for b in sorted(buckets.keys()):
            bar = '#' * (buckets[b] * 40 // max(buckets.values()))
            print(f"    {b:3d}-{b+9:3d}: {buckets[b]:4d} {bar}")

    # ── Phase 4: Compare with random (non-kernel) perturbations ───

    print()
    print("Phase 4: Control — random (non-kernel) perturbations")
    print("-" * 72)

    random_hws = []
    for _ in range(NUM_TESTS):
        base_msg = [random.getrandbits(32) for _ in range(16)]
        delta_msg = [random.getrandbits(32) for _ in range(16)]
        perturbed_msg = [(base_msg[i] ^ delta_msg[i]) & MASK for i in range(16)]

        h1 = sha256_compress(base_msg, iv)
        h2 = sha256_compress(perturbed_msg, iv)

        diff = sum(hw(h1[i] ^ h2[i]) for i in range(8))
        random_hws.append(diff)

    if random_hws:
        avg_random = sum(random_hws) / len(random_hws)
        print(f"  Average random perturbation: {avg_random:.1f} / 256 bits")
        print(f"  Min: {min(random_hws)}, Max: {max(random_hws)}")

    # ── Phase 5: Real Jacobian kernel ─────────────────────────────

    print()
    print("Phase 5: Real (nonlinear) Jacobian kernel analysis")
    print("-" * 72)

    # The real Jacobian at a specific point might have lower rank
    # than the linearized version (because carries add nonlinear terms)
    real_ranks = []
    real_perturbations = []

    for msg_idx in range(5):
        base_msg = [random.getrandbits(32) for _ in range(16)]
        real_matrix = build_real_jacobian(base_msg, iv)
        _, real_rank = gf2_kernel(real_matrix)
        real_ranks.append(real_rank)

        # Test real kernel elements
        real_kernel, _ = gf2_kernel(real_matrix)
        if real_kernel:
            base_hash = sha256_compress(base_msg, iv)
            for trial in range(min(20, len(real_kernel))):
                delta_vec = real_kernel[trial]
                delta_msg = vec_to_msg(delta_vec)
                perturbed_msg = [(base_msg[i] ^ delta_msg[i]) & MASK for i in range(16)]
                perturbed_hash = sha256_compress(perturbed_msg, iv)
                diff = sum(hw(base_hash[i] ^ perturbed_hash[i]) for i in range(8))
                real_perturbations.append(diff)

        elapsed = time.time() - t_start
        print(f"  Msg {msg_idx+1}: real Jacobian rank = {real_rank} [{elapsed:.1f}s]")

    if real_perturbations:
        avg_real = sum(real_perturbations) / len(real_perturbations)
        print(f"\n  Real Jacobian kernel perturbation: avg = {avg_real:.1f}, "
              f"min = {min(real_perturbations)}, max = {max(real_perturbations)}")

    # ── Phase 6: Reduced rounds analysis ──────────────────────────

    print()
    print("Phase 6: Carry perturbation vs number of rounds")
    print("-" * 72)

    def sha256_compress_rounds(msg, iv, num_rounds):
        W = expand_W(msg)
        state = list(iv)
        for t in range(num_rounds):
            a, b, c, d, e, f, g, h = state
            T1 = add32(h, Sig1(e), Ch(e, f, g), K[t], W[t])
            T2 = add32(Sig0(a), Maj(a, b, c))
            state = [add32(T1, T2), a, b, c, add32(d, T1), e, f, g]
        return [add32(state[i], iv[i]) for i in range(8)]

    def sha256_compress_linear_rounds(msg, iv, num_rounds):
        W = expand_W_linear(msg)
        state = list(iv)
        for t in range(num_rounds):
            a, b, c, d, e, f, g, h = state
            T1 = xor32(h, Sig1(e), Ch(e, f, g), K[t], W[t])
            T2 = xor32(Sig0(a), Maj(a, b, c))
            state = [T1 ^ T2, a, b, c, d ^ T1, e, f, g]
        return [state[i] ^ iv[i] for i in range(8)]

    for num_rounds in [1, 2, 4, 8, 16, 24, 32, 48, 64]:
        if time.time() - t_start > 300:
            break

        # Build linearized matrix for this many rounds
        zero_msg = [0] * 16
        base_lin = sha256_compress_linear_rounds(zero_msg, iv, num_rounds)
        base_lin_hash = 0
        for i, val in enumerate(base_lin):
            base_lin_hash |= (val << (i * 32))

        matrix_r = []
        for word in range(16):
            for bit in range(32):
                msg = [0] * 16
                msg[word] = 1 << bit
                out = sha256_compress_linear_rounds(msg, iv, num_rounds)
                h = 0
                for i, val in enumerate(out):
                    h |= (val << (i * 32))
                matrix_r.append(base_lin_hash ^ h)

        kernel_r, rank_r = gf2_kernel(matrix_r)

        # Test perturbation
        perturb_samples = []
        for trial in range(min(50, max(1, len(kernel_r)))):
            # Random kernel element
            delta_vec = 0
            for kb in kernel_r:
                if random.random() < 0.5:
                    delta_vec ^= kb
            if delta_vec == 0:
                continue

            base_msg = [random.getrandbits(32) for _ in range(16)]
            delta_msg = vec_to_msg(delta_vec)
            perturbed_msg = [(base_msg[i] ^ delta_msg[i]) & MASK for i in range(16)]

            h1 = sha256_compress_rounds(base_msg, iv, num_rounds)
            h2 = sha256_compress_rounds(perturbed_msg, iv, num_rounds)
            diff = sum(hw(h1[i] ^ h2[i]) for i in range(8))
            perturb_samples.append(diff)

        if perturb_samples:
            avg_p = sum(perturb_samples) / len(perturb_samples)
            min_p = min(perturb_samples)
        else:
            avg_p = -1
            min_p = -1

        print(f"  Rounds {num_rounds:2d}: rank={rank_r:3d}, kernel_dim={len(kernel_r):3d}, "
              f"carry_perturbation: avg={avg_p:.1f}, min={min_p}")

    # ── Verdict ───────────────────────────────────────────────────

    print()
    print("=" * 72)
    print("VERDICT")
    print("=" * 72)
    print()

    if perturbation_hws:
        avg_hw = sum(perturbation_hws) / len(perturbation_hws)
        if avg_hw < 100:
            print(f"  **ALIVE** — Carry perturbation {avg_hw:.1f} < 128 (random)!")
            print(f"  Linearized kernel gives BETTER-THAN-RANDOM starting points")
            print(f"  Potential for correction-based attack")
        elif avg_hw < 120:
            print(f"  **ANOMALY** — Carry perturbation {avg_hw:.1f}, slightly below random 128")
        else:
            print(f"  **DEAD** — Carry perturbation {avg_hw:.1f} ≈ 128 (random)")
            print(f"  Linearized kernel provides no advantage")

    if real_perturbations:
        avg_real = sum(real_perturbations) / len(real_perturbations)
        if avg_real < avg_hw - 10:
            print(f"\n  Real Jacobian kernel: {avg_real:.1f} — BETTER than linearized!")
        else:
            print(f"\n  Real Jacobian kernel: {avg_real:.1f} — similar to linearized")

    print()
    print(f"  Runtime: {time.time() - t_start:.1f}s")
    print("=" * 72)


if __name__ == "__main__":
    main()
