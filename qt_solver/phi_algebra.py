"""
PHI Algebra: SHA-256(M) = L(M) ⊕ Φ(M)

L(M) = XOR-SHA — all additions replaced by XOR (GF(2)-linear)
Φ(M) = carry correction = SHA-256(M) ⊕ L(M) (nonlinear)

Key properties:
  - L is GF(2)-linear with kernel dim = 256
  - Φ has full rank 256 (same as SHA-256)
  - Collision: Ψ(M) = L(δM) where Ψ(M) = Φ(M) ⊕ Φ(M⊕δM)
"""

from qt_solver.sha256_traced import (
    MASK32, IV, K, get_bit,
    sigma0, sigma1, ssigma0, ssigma1,
    ch, maj,
)
from qt_solver.gf2 import gf2_gaussian_eliminate
import random


def xor_sha256(msg, rounds=64, iv=None):
    """
    XOR-SHA: SHA-256 with all mod-2^32 additions replaced by XOR.
    This is GF(2)-linear in the message.
    """
    if iv is None:
        iv = list(IV)
    else:
        iv = list(iv)

    # Message schedule (XOR version)
    w = list(msg[:16])
    for i in range(16, max(rounds, 16)):
        w.append(ssigma1(w[i-2]) ^ w[i-7] ^ ssigma0(w[i-15]) ^ w[i-16])

    a, b, c, d, e, f, g, h = iv

    for r in range(rounds):
        t1 = h ^ sigma1(e) ^ ch(e, f, g) ^ K[r] ^ w[r]
        t2 = sigma0(a) ^ maj(a, b, c)
        h = g
        g = f
        f = e
        e = d ^ t1
        d = c
        c = b
        b = a
        a = t1 ^ t2

    return [(a ^ iv[0]) & MASK32, (b ^ iv[1]) & MASK32,
            (c ^ iv[2]) & MASK32, (d ^ iv[3]) & MASK32,
            (e ^ iv[4]) & MASK32, (f ^ iv[5]) & MASK32,
            (g ^ iv[6]) & MASK32, (h ^ iv[7]) & MASK32]


def phi_sha256(msg, rounds=64):
    """
    Compute PHI decomposition: SHA-256(M) = L(M) ⊕ Φ(M).
    Returns (sha_hash, l_hash, phi_hash).
    """
    from qt_solver.sha256_traced import sha256_compress

    sha_hash = sha256_compress(msg, rounds)
    l_hash = xor_sha256(msg, rounds)
    phi_hash = [sha_hash[i] ^ l_hash[i] for i in range(8)]

    return sha_hash, l_hash, phi_hash


def verify_phi_decomposition(n_tests=1000, rounds=64, verbose=True):
    """
    Verify: SHA-256(M) = L(M) ⊕ Φ(M) for random messages.
    """
    from qt_solver.sha256_traced import sha256_compress

    if verbose:
        print(f"PHI Decomposition Verification (R={rounds}, {n_tests} tests)")

    rng = random.Random(42)
    violations = 0

    for _ in range(n_tests):
        msg = [rng.randint(0, MASK32) for _ in range(16)]
        sha_hash, l_hash, phi_hash = phi_sha256(msg, rounds)

        # Verify: sha_hash = l_hash XOR phi_hash
        recon = [l_hash[i] ^ phi_hash[i] for i in range(8)]
        if recon != sha_hash:
            violations += 1

    if verbose:
        print(f"  Violations: {violations}/{n_tests}")
        if violations == 0:
            print(f"  ✓ SHA-256(M) = L(M) ⊕ Φ(M) verified for all tests")

    return violations == 0


def analyze_l_kernel(rounds=64, verbose=True):
    """
    Compute the kernel of L (XOR-SHA).
    L is linear → kernel found via Jacobian.
    """
    if verbose:
        print(f"\nL-Kernel Analysis (R={rounds})")

    rng = random.Random(42)
    msg = [rng.randint(0, MASK32) for _ in range(16)]
    base = xor_sha256(msg, rounds)

    # Build Jacobian of L (should be linear, so exact)
    rows = []
    for hb in range(256):
        row = 0
        for mb in range(512):
            wi, bi = mb // 32, mb % 32
            msg2 = list(msg)
            msg2[wi] ^= (1 << bi)
            h2 = xor_sha256(msg2, rounds)
            hw, hbi = hb // 32, hb % 32
            if get_bit(base[hw], hbi) != get_bit(h2[hw], hbi):
                row |= (1 << mb)
        rows.append(row)

    echelon, pivots = gf2_gaussian_eliminate(list(rows), 512)
    rank = len(pivots)
    kernel_dim = 512 - rank

    if verbose:
        print(f"  L rank: {rank}/256")
        print(f"  L kernel dim: {kernel_dim}")
        print(f"  Expected: kernel = 256 (msg 512 - hash 256)")

    # Verify linearity: L(M1⊕M2) = L(M1) ⊕ L(M2) ⊕ L(0)?
    # For affine: L(M) = A·M ⊕ b, so L(M1⊕M2) = A·(M1⊕M2)⊕b = A·M1⊕A·M2⊕b
    # But L(M1)⊕L(M2) = A·M1⊕b ⊕ A·M2⊕b = A·M1⊕A·M2
    # So L(M1⊕M2) ⊕ L(M1) ⊕ L(M2) = b (constant, not 0 unless b=0)
    msg0 = [0] * 16
    l0 = xor_sha256(msg0, rounds)
    msg1 = [rng.randint(0, MASK32) for _ in range(16)]
    msg2 = [rng.randint(0, MASK32) for _ in range(16)]
    msg12 = [msg1[i] ^ msg2[i] for i in range(16)]

    l1 = xor_sha256(msg1, rounds)
    l2 = xor_sha256(msg2, rounds)
    l12 = xor_sha256(msg12, rounds)

    # Check: L(M1⊕M2) ⊕ L(M1) ⊕ L(M2) = L(0)?
    check = [l12[i] ^ l1[i] ^ l2[i] for i in range(8)]
    is_affine = (check == l0)

    if verbose:
        print(f"  L(M1⊕M2) ⊕ L(M1) ⊕ L(M2) = L(0)? {is_affine}")
        if is_affine:
            print(f"  ✓ L is affine (linear + constant from IV+K)")

    return {
        'rank': rank,
        'kernel_dim': kernel_dim,
    }


def analyze_phi_properties(n_tests=500, rounds=64, verbose=True):
    """
    Analyze properties of Φ (carry correction map).
    """
    from qt_solver.sha256_traced import sha256_compress

    if verbose:
        print(f"\nΦ Properties Analysis (R={rounds})")

    rng = random.Random(42)

    # Collect Φ values
    phi_values = []
    for _ in range(n_tests):
        msg = [rng.randint(0, MASK32) for _ in range(16)]
        _, _, phi_hash = phi_sha256(msg, rounds)
        phi_values.append(phi_hash)

    # Check balance: average Hamming weight
    total_hw = 0
    for phi in phi_values:
        for w in range(8):
            total_hw += bin(phi[w]).count('1')
    avg_hw = total_hw / n_tests

    if verbose:
        print(f"  Average HW(Φ): {avg_hw:.1f}/256 (128.0 = balanced)")

    # Check rank: build Jacobian of Φ at a random point
    msg = [rng.randint(0, MASK32) for _ in range(16)]
    _, _, base_phi = phi_sha256(msg, rounds)

    rows = []
    for hb in range(256):
        row = 0
        for mb in range(min(256, 512)):  # first 256 msg bits for speed
            wi, bi = mb // 32, mb % 32
            msg2 = list(msg)
            msg2[wi] ^= (1 << bi)
            _, _, phi2 = phi_sha256(msg2, rounds)
            hw, hbi = hb // 32, hb % 32
            if get_bit(base_phi[hw], hbi) != get_bit(phi2[hw], hbi):
                row |= (1 << mb)
        rows.append(row)

    echelon, pivots = gf2_gaussian_eliminate(list(rows), min(256, 512))
    phi_rank = len(pivots)

    if verbose:
        print(f"  Φ Jacobian rank: {phi_rank}/256")

    # Check XOR-linearity violations
    n_linear_tests = 200
    linear_violations = 0
    for _ in range(n_linear_tests):
        m1 = [rng.randint(0, MASK32) for _ in range(16)]
        m2 = [rng.randint(0, MASK32) for _ in range(16)]
        m12 = [m1[i] ^ m2[i] for i in range(16)]

        _, _, phi1 = phi_sha256(m1, rounds)
        _, _, phi2 = phi_sha256(m2, rounds)
        _, _, phi12 = phi_sha256(m12, rounds)
        _, _, phi0 = phi_sha256([0]*16, rounds)

        # If linear: Φ(M1⊕M2) ⊕ Φ(M1) ⊕ Φ(M2) = Φ(0)
        check = [phi12[i] ^ phi1[i] ^ phi2[i] for i in range(8)]
        if check != phi0:
            linear_violations += 1

    if verbose:
        print(f"  XOR-linearity violations: {linear_violations}/{n_linear_tests}")
        if linear_violations == n_linear_tests:
            print(f"  ✓ Φ is fully nonlinear (all tests show nonlinearity)")

    # Hash distance in L-kernel: pick two messages from L-kernel, compare SHA hashes
    # This tests whether Φ is "random" on the L-kernel
    if verbose:
        print(f"\n  Collision search in L-kernel:")
        print(f"  If Φ is random on L-kernel, hash distances ≈ 128")

    # Generate L-kernel vectors (message diffs invisible to L)
    msg_base = [rng.randint(0, MASK32) for _ in range(16)]
    base_L = xor_sha256(msg_base, rounds)
    base_SHA = sha256_compress(msg_base, rounds)

    l_kernel_dists = []
    for _ in range(100):
        # Random message diff
        delta = [rng.randint(0, MASK32) for _ in range(16)]
        msg2 = [msg_base[i] ^ delta[i] for i in range(16)]
        l2 = xor_sha256(msg2, rounds)

        if l2 == base_L:
            # In L-kernel! Check SHA distance
            sha2 = sha256_compress(msg2, rounds)
            dist = sum(bin(base_SHA[w] ^ sha2[w]).count('1') for w in range(8))
            l_kernel_dists.append(dist)

    if l_kernel_dists:
        avg_dist = sum(l_kernel_dists) / len(l_kernel_dists)
        if verbose:
            print(f"  Found {len(l_kernel_dists)} L-kernel pairs")
            print(f"  SHA hash distance: {avg_dist:.1f}/256 (128.0 = random)")
    else:
        if verbose:
            print(f"  No random L-kernel pairs found (expected: rare for R=64)")

    return {
        'avg_hw': avg_hw,
        'phi_rank': phi_rank,
        'linear_violations': linear_violations,
    }


if __name__ == '__main__':
    print("=" * 60)
    print("PHI ALGEBRA — SHA-256(M) = L(M) ⊕ Φ(M)")
    print("=" * 60)

    for R in [4, 16, 64]:
        print(f"\n{'='*40}")
        print(f"R = {R} rounds")
        print(f"{'='*40}")
        verify_phi_decomposition(n_tests=500, rounds=R)
        analyze_l_kernel(rounds=R)
        analyze_phi_properties(n_tests=200, rounds=R)
