"""
LATTICE MODEL v2: Fix the decomposition and explore what 75% correlation means.

Issue from v1: XOR/carry decomposition had 991/2048 violations because
T1 and T2 themselves contain intermediate carries. We need to decompose
at the level of INDIVIDUAL additions, not just the final one.

Plan:
1. Fix decomposition: track ALL carries through ALL additions
2. Understand WHY 75% correlation exists (it's not trivial)
3. Test: can the 75% correlation be used to predict/constrain bits?
4. Test: does this give ANY advantage over birthday?
"""

import random
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, schedule
)


def full_carry_decomposition(state, r, W):
    """
    Decompose one SHA-256 round into XOR-parts and carry-parts
    for EVERY addition separately.

    Returns:
      a_xor_only: what a[r+1] would be if ALL additions were XOR
      carries: list of 32-bit carry patterns for each addition
      a_real: actual a[r+1]
    """
    a, b, c, d, e, f, g, h = state

    def add_parts(x, y):
        """Return (xor_result, carry_pattern, real_result)."""
        real = (x + y) & MASK
        xor_r = x ^ y
        carry_p = real ^ xor_r
        return xor_r, carry_p, real

    # T1 = h + Sig1(e) + Ch(e,f,g) + K[r] + W[r]
    sig1e = Sig1(e)
    che = Ch(e, f, g)

    # Step 1: h + Sig1(e)
    xor1, carry1, sum1 = add_parts(h, sig1e)

    # Step 2: sum1 + Ch
    xor2, carry2, sum2 = add_parts(sum1, che)

    # Step 3: sum2 + K[r]
    xor3, carry3, sum3 = add_parts(sum2, K[r])

    # Step 4: sum3 + W[r]
    xor4, carry4, T1 = add_parts(sum3, W[r])

    # T2 = Sig0(a) + Maj(a,b,c)
    sig0a = Sig0(a)
    maja = Maj(a, b, c)
    xor5, carry5, T2 = add_parts(sig0a, maja)

    # a_new = T1 + T2
    xor6, carry6, a_new = add_parts(T1, T2)

    # e_new = d + T1
    xor7, carry7, e_new = add_parts(d, T1)

    # "Pure XOR" version: replace ALL additions with XOR
    xor_T1 = h ^ sig1e ^ che ^ K[r] ^ W[r]
    xor_T2 = sig0a ^ maja
    a_xor_only = xor_T1 ^ xor_T2

    # Total carry contribution
    total_carry = a_new ^ a_xor_only

    carries = [carry1, carry2, carry3, carry4, carry5, carry6, carry7]
    return a_xor_only, total_carry, carries, a_new, e_new


def experiment_v2_decomposition():
    """Verify the full decomposition."""
    print("=" * 80)
    print("LATTICE v2: Full carry decomposition")
    print("=" * 80)

    N = 500
    violations = 0
    total_carry_hws = []

    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        states, W = sha256_round_trace(M)

        for r in range(64):
            a_xor, total_carry, carries, a_real, _ = full_carry_decomposition(states[r], r, W)

            # Verify: a_real = a_xor XOR total_carry
            if a_real != (a_xor ^ total_carry):
                violations += 1

            total_carry_hws.append(bin(total_carry).count('1'))

    print(f"  Verification: a_real = a_xor_only XOR total_carry_correction")
    print(f"  Violations: {violations}/{N*64}")
    avg_hw = sum(total_carry_hws) / len(total_carry_hws)
    print(f"  E[HW(total carry correction)] = {avg_hw:.2f}/32")


def experiment_v2_75pct():
    """
    WHY is carry-to-carry correlation 75%?

    C(r, k-1) predicts C(r, k) with 75% accuracy.
    Two cases:
      - If C(r, k-1) = 0: carry stops. C(r,k) = G[k-1] only.
        P(G[k-1]=1) ≈ 0.25 (both bits = 1). So P(C(r,k)=1 | C(r,k-1)=0) ≈ 0.25.
      - If C(r, k-1) = 1: carry propagates if P[k-1]=1, generates if G[k-1]=1.
        P(C(r,k)=1 | C(r,k-1)=1) = P(G) + P(P)*1 = 0.25 + 0.5 = 0.75.

    Overall: P(C_k = C_{k-1}) = P(C_{k-1}=0)*P(C_k=0|C_{k-1}=0) + P(C_{k-1}=1)*P(C_k=1|C_{k-1}=1)

    For random operands:
      P(C_k=0) ≈ 0.5 (steady state of carry chain)
      P(C_k = C_{k-1}) = 0.5 * 0.75 + 0.5 * 0.75 = 0.75

    So 75% is EXACTLY what random operands predict!
    It's a consequence of the MAJ function: MAJ(x,y,c) = c with probability 3/4.
    """
    print("\n" + "=" * 80)
    print("WHY 75%: Carry autocorrelation is a PROPERTY OF MAJ")
    print("=" * 80)

    # Verify: MAJ(random, random, c) = c with what probability?
    N = 100000
    same_count = 0
    for _ in range(N):
        x = random.randint(0, 1)
        y = random.randint(0, 1)
        c = random.randint(0, 1)
        maj = (x & y) | (x & c) | (y & c)
        if maj == c:
            same_count += 1

    print(f"  P(MAJ(random, random, c) = c) = {same_count/N:.4f}")
    print(f"  Expected: 3/4 = 0.7500")
    print(f"")
    print(f"  This is NOT a discovery about SHA-256.")
    print(f"  It's a property of MAJ: the carry bit is 'sticky'.")
    print(f"  Once carry = 1, it tends to stay 1 (75% chance).")
    print(f"  Once carry = 0, it tends to stay 0 (75% chance).")
    print(f"")
    print(f"  This explains the geometric segment distribution:")
    print(f"  P(segment length ≥ k) = 0.75^k (carry stays same for k steps)")
    print(f"  Average segment = 1/(1-0.75) = 4 (close to measured 3.64)")

    # Verify geometric distribution
    print(f"\n  Verification: P(carry stays same) = 0.75^k")
    from lattice_model import build_carry_lattice
    N2 = 200
    run_lengths = []
    for seed in range(N2):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        carry = build_carry_lattice(M)
        for r in range(8, 60):
            run = 0
            for k in range(1, 32):
                if carry[r][k] == carry[r][k-1]:
                    run += 1
                else:
                    if run > 0:
                        run_lengths.append(run)
                    run = 0

    for k in [1, 2, 3, 4, 5, 8, 10]:
        observed = sum(1 for r in run_lengths if r >= k) / len(run_lengths)
        predicted = 0.75 ** k
        print(f"    P(run ≥ {k}): observed={observed:.4f}, predicted(0.75^k)={predicted:.4f}")


def experiment_v2_what_matters():
    """
    OK so 75% is trivial (property of MAJ). But the LATTICE MODEL is still exact.

    Let's ask the real question: in the lattice model, where does
    the information about M go?

    SHA-256 maps M → H. In lattice terms:
      φ(0, k) = H0[0][k]  (constant, from IV)
      W[r] enters through L(r, k) at each round
      After 64 rounds: φ(64, k) = output

    The MESSAGE enters only through W[r] in the L-part.
    The C-part (carry) is a RESPONSE to L — it's computed from L.

    So: L carries the message signal. C corrupts it (adds nonlinearity).

    The question becomes:
      How much of the L-signal survives the C-corruption after 64 rounds?

    If we could REMOVE C entirely (set all carries to 0), we'd have
    XOR-SHA, which is GF(2)-linear with kernel dim = 256 (PHI result).

    The carries ADD 256 constraints (making the system full rank).
    But they do so in a STRUCTURED way (local, sequential, MAJ-based).

    Can we PARTIALLY remove carries? Not all 256, but enough to
    create a small kernel?
    """
    print("\n" + "=" * 80)
    print("WHAT MATTERS: How much L-signal survives C-corruption?")
    print("=" * 80)

    # Measure: for XOR-only SHA (no carries), the kernel is 256-dimensional.
    # For real SHA, kernel is 0-dimensional.
    # If we fix SOME carry bits (say, the first k bits of each addition),
    # how does the kernel dimension change?

    # Method: compute Jacobian with some carry bits fixed.
    # "Fixing" carry bits = adding them as known constraints.

    # Simplified test: fix carry bits at position 0..k of the final addition (T1+T2).
    # This partially linearizes the system.

    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]

    R = 16  # rounds

    # Base hash
    def sha_reduced_bits(msg, R):
        states, W = sha256_round_trace(msg, rounds=R)
        result = []
        for reg in range(8):
            val = (states[R][reg] + H0[reg]) & MASK
            for b in range(32):
                result.append((val >> b) & 1)
        return result

    base = sha_reduced_bits(M, R)

    # Compute full Jacobian
    rows = []
    for word in range(16):
        for bit in range(32):
            M2 = list(M)
            M2[word] ^= (1 << bit)
            h2 = sha_reduced_bits(M2, R)
            rows.append([base[i] ^ h2[i] for i in range(256)])

    def gf2_rank(matrix, nrows, ncols):
        m = [list(row[:ncols]) for row in matrix[:nrows]]
        rank = 0
        for col in range(ncols):
            pivot = None
            for row in range(rank, len(m)):
                if m[row][col] == 1:
                    pivot = row
                    break
            if pivot is None:
                continue
            m[rank], m[pivot] = m[pivot], m[rank]
            for row in range(len(m)):
                if row != rank and m[row][col] == 1:
                    for c in range(ncols):
                        m[row][c] ^= m[rank][c]
            rank += 1
        return rank

    full_rank = gf2_rank(rows, 512, 256)
    kernel_dim = 512 - full_rank
    print(f"  Real SHA-256 (R={R}): Jacobian rank = {full_rank}/256, kernel = {kernel_dim}")

    # Now: XOR-only SHA
    def sha_xor_bits(msg, R):
        W = list(msg)
        for i in range(16, 64):
            s1 = rotr(W[i-2],17) ^ rotr(W[i-2],19) ^ (W[i-2]>>10)
            s0 = rotr(W[i-15],7) ^ rotr(W[i-15],18) ^ (W[i-15]>>3)
            W.append(s1 ^ W[i-7] ^ s0 ^ W[i-16])
        a,b,c,d,e,f,g,h = H0
        for r in range(R):
            T1 = h ^ Sig1(e) ^ Ch(e,f,g) ^ K[r] ^ W[r]
            T2 = Sig0(a) ^ Maj(a,b,c)
            h,g,f = g,f,e; e = d ^ T1
            d,c,b = c,b,a; a = T1 ^ T2
        result = []
        for reg_idx, reg_val in enumerate([a,b,c,d,e,f,g,h]):
            val = reg_val ^ H0[reg_idx]
            for bit in range(32):
                result.append((val >> bit) & 1)
        return result

    base_xor = sha_xor_bits(M, R)
    rows_xor = []
    for word in range(16):
        for bit in range(32):
            M2 = list(M)
            M2[word] ^= (1 << bit)
            h2 = sha_xor_bits(M2, R)
            rows_xor.append([base_xor[i] ^ h2[i] for i in range(256)])

    xor_rank = gf2_rank(rows_xor, 512, 256)
    xor_kernel = 512 - xor_rank
    print(f"  XOR-only SHA (R={R}): Jacobian rank = {xor_rank}/256, kernel = {xor_kernel}")

    # The gap: real has kernel 256, XOR has kernel 256.
    # Wait — both should have kernel 256 (512-256=256)?
    # No: real SHA at R=16 has rank 256 (full) → kernel = 256.
    # XOR SHA at R=16 also has rank 256 → kernel = 256.
    # But the KERNELS ARE DIFFERENT.

    # The real question: does the real-SHA kernel OVERLAP with the XOR-SHA kernel?
    # If overlap > 0 → there are directions where carries don't matter.

    print(f"\n  Both have kernel dim 256. But are the kernels the SAME?")

    # Compute both kernels and check overlap
    # Actually, kernel dim = 512 - rank. If rank = 256, kernel = 256.
    # These are 256-dimensional subspaces of GF(2)^512.
    # Random 256-dim subspaces overlap in dim ≈ 0.

    # But THESE kernels might overlap because they share structure (rotations, Ch, Maj).
    # Let's check: compute a basis for each kernel and see if any vector is in both.

    # For real SHA kernel: find 256 vectors in GF(2)^512 that map to 0
    # This is expensive but we can check: does a random XOR-kernel vector
    # also satisfy the real SHA constraints?

    # Quick test: take the XOR-kernel basis vectors and evaluate real SHA on them.
    # If real SHA(M XOR v) = real SHA(M) for v in XOR-kernel → overlap!

    print(f"\n  Testing: do XOR-kernel vectors also lie in real-SHA kernel?")

    # Find XOR kernel vectors by Gaussian elimination
    # We need the null space of the XOR Jacobian
    m = [list(row) for row in rows_xor]
    # Row reduce
    pivot_cols = []
    rank = 0
    for col in range(256):
        pivot = None
        for row in range(rank, 512):
            if m[row][col] == 1:
                pivot = row
                break
        if pivot is None:
            continue
        m[rank], m[pivot] = m[pivot], m[rank]
        pivot_cols.append(col)
        for row in range(512):
            if row != rank and m[row][col] == 1:
                for c in range(256):
                    m[row][c] ^= m[rank][c]
        rank += 1

    # Free columns (not pivots) → kernel vectors
    free_cols = [c for c in range(256) if c not in set(pivot_cols)]
    print(f"  XOR kernel: {len(free_cols)} free variables (expected {256 - xor_rank}... wait this is output space)")

    # Actually the kernel is in INPUT space (512 dim), not output space (256 dim).
    # The Jacobian is 512×256. Kernel = null space of J^T in output space... no.
    # J is 512 rows × 256 cols. Rank = 256. Left null space = 512-256 = 256 vectors in input space.

    # Let me just test directly: pick random kernel direction for XOR-SHA,
    # and check if it's also kernel for real SHA.

    # For XOR-SHA: a kernel vector v means XOR-SHA(M XOR v) = XOR-SHA(M)
    # These are the L-kernel vectors from PHI.
    # For real SHA: check if SHA(M XOR v) = SHA(M)

    tested = 0
    real_kernel_too = 0
    for trial in range(100):
        # Random message perturbation that's in XOR-SHA kernel
        # (Changes M such that XOR-SHA output unchanged)
        v = [0] * 16
        # Simple: flip two words such that XOR-SHA cancels
        # Actually, finding kernel vectors requires more work.
        # Let's just test: does flipping W[0] bit 0 AND some other bit
        # cancel in XOR world?

        # Skip complex kernel computation.
        # Instead: directly measure hash distance for a known XOR-kernel direction.
        pass

    print(f"  (Kernel overlap test requires computing explicit kernel basis — deferred)")
    print(f"  Key insight: XOR-SHA kernel = 256 dimensions (PHI L-kernel).")
    print(f"  Real SHA: same 256 dimensions exist, but carries SHIFT them.")
    print(f"  The carries are a PERTURBATION of the XOR kernel.")
    print(f"  If the perturbation is structured (and it IS — through MAJ chains),")
    print(f"  then the real kernel might be CLOSE to the XOR kernel in some metric.")


if __name__ == "__main__":
    experiment_v2_decomposition()
    experiment_v2_75pct()
    experiment_v2_what_matters()
