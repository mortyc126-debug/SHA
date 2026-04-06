"""
OMEGA Solver: linearize the quadratic system using a known assignment,
then find the α-kernel (message-changing directions).

Strategy:
1. Build OMEGA system (linear + quadratic equations)
2. Given a known message M, compute full SHA-256 trace
3. Use trace to LINEARIZE quadratic equations around M
4. Solve linearized system → find kernel
5. Filter kernel to message-changing directions → α-kernel
6. For each α-kernel vector, XOR with M to get second preimage candidate
7. Verify candidate hashes to same value
"""

import random
from qt_solver.sha256_traced import (
    MASK32, get_bit, sha256_compress, sha256_compress_traced, IV,
)
from qt_solver.gf2 import gf2_solve, gf2_kernel, gf2_gaussian_eliminate, bitvec_weight


def omega_linearize_and_solve(R, msg, verbose=True):
    """
    Build OMEGA system for R rounds, linearize using msg's trace,
    find α-kernel.

    Instead of the full symbolic OmegaSystem (which is complex),
    we use a DIRECT approach:

    1. Trace SHA-256 on msg to get all intermediate values
    2. Build GF(2) Jacobian: how does each message bit affect each hash bit?
    3. The kernel of this Jacobian = directions that don't change the hash
    4. α-kernel = kernel intersected with message space

    This is equivalent to OMEGA linearization but much simpler to implement.

    Returns dict with α-kernel info.
    """
    trace = sha256_compress_traced(msg, R)
    target_hash = trace['hash']
    n_msg = 512  # 16 words × 32 bits

    if verbose:
        print(f"OMEGA linearize for R={R}")
        print(f"  Message: {[hex(w) for w in msg[:4]]}...")
        print(f"  Hash: {[hex(w) for w in target_hash[:4]]}...")

    # Build Jacobian by flipping each message bit and observing hash change
    # This gives us the EXACT linear approximation around msg
    hash_bits = 256  # 8 words × 32 bits
    jacobian_rows = []  # one row per hash bit, columns = message bits

    # Reference hash as bit vector
    ref_hash = 0
    for w in range(8):
        for b in range(32):
            if get_bit(target_hash[w], b):
                ref_hash |= (1 << (w * 32 + b))

    # For each message bit, compute hash difference
    msg_effects = [0] * n_msg  # msg_effects[i] = which hash bits flip when msg bit i flips

    for wi in range(16):
        for bi in range(32):
            msg_idx = wi * 32 + bi
            # Flip bit
            msg_copy = list(msg)
            msg_copy[wi] ^= (1 << bi)
            new_hash = sha256_compress(msg_copy, R)

            diff = 0
            for w in range(8):
                xd = target_hash[w] ^ new_hash[w]
                for b in range(32):
                    if get_bit(xd, b):
                        diff |= (1 << (w * 32 + b))
            msg_effects[msg_idx] = diff

    # Build hash-bit equations: for each hash bit, which message bits affect it?
    # Row per hash bit, column per message bit
    # We want: J * δM = δH = 0 (for second preimage, δH = 0)
    rows = []
    for hb in range(hash_bits):
        row = 0
        for mb in range(n_msg):
            if (msg_effects[mb] >> hb) & 1:
                row |= (1 << mb)
        rows.append(row)

    # Find kernel of J (directions where δH = 0)
    kernel = gf2_kernel(rows, n_msg)

    if verbose:
        echelon, pivots = gf2_gaussian_eliminate(list(rows), n_msg)
        rank = len(pivots)
        print(f"  Jacobian rank: {rank}/{hash_bits}")
        print(f"  Kernel dimension: {len(kernel)}")
        print(f"  Expected α-kernel: {max(0, 32*(16-R))}")

    return {
        'kernel': kernel,
        'alpha_dim': len(kernel),
        'expected_alpha': max(0, 32 * (16 - R)),
        'msg': msg,
        'hash': target_hash,
        'R': R,
        'jacobian_rows': rows,
    }


def find_second_preimage(R, msg=None, seed=42, verbose=True):
    """
    Find a second preimage for R-round SHA-256.

    For R <= 15: should find algebraic preimage via α-kernel.
    For R = 16: α-kernel = 0, no algebraic preimage.
    """
    if msg is None:
        rng = random.Random(seed)
        msg = [rng.randint(0, MASK32) for _ in range(16)]

    result = omega_linearize_and_solve(R, msg, verbose)
    kernel = result['kernel']

    if not kernel:
        if verbose:
            print(f"  α-kernel empty — no algebraic second preimage at R={R}")
        return result

    # Try each kernel vector as a message perturbation
    target_hash = result['hash']
    found = []

    for i, kvec in enumerate(kernel):
        # Build perturbed message
        msg2 = list(msg)
        for wi in range(16):
            word_bits = 0
            for bi in range(32):
                if (kvec >> (wi * 32 + bi)) & 1:
                    word_bits |= (1 << bi)
            msg2[wi] ^= word_bits

        if msg2 == msg:
            continue  # trivial (shouldn't happen if kernel is message-changing)

        # Verify hash
        hash2 = sha256_compress(msg2, R)
        hash_diff = sum(bin(target_hash[w] ^ hash2[w]).count('1') for w in range(8))

        if hash_diff == 0:
            found.append(msg2)
            if verbose and len(found) <= 3:
                print(f"  ✓ Second preimage #{len(found)}: msg differs in {bitvec_weight(kvec)} bits, hash MATCHES")
        else:
            if verbose and i < 3:
                print(f"  ✗ Kernel vec {i}: hash differs in {hash_diff} bits (linearization error)")

    result['preimages'] = found
    result['n_found'] = len(found)
    if verbose:
        print(f"  Total: {len(found)} second preimages from {len(kernel)} kernel vectors")
        # Note: linearization is approximate for R > ~5 due to quadratic terms
        # Exact preimages need full OMEGA system, not just Jacobian

    return result


def verify_omega_range(R_min=1, R_max=15, seeds=3, verbose=True):
    """
    Verify OMEGA finds second preimages for R=R_min..R_max.
    """
    results = {}
    for R in range(R_min, R_max + 1):
        if verbose:
            print(f"\n{'='*60}")
            print(f"R = {R}")
            print(f"{'='*60}")
        for seed in range(seeds):
            res = find_second_preimage(R, seed=seed + 100*R, verbose=verbose)
            key = (R, seed)
            results[key] = res
    return results


if __name__ == '__main__':
    print("OMEGA System Verification")
    print("=" * 60)
    verify_omega_range(R_min=1, R_max=16, seeds=2, verbose=True)
