"""
PARALLEL: Three research threads at once.

THREAD A: Prove T3 (nilpotency) analytically.
THREAD B: BTE "spectrum" — eigenvalues of round operator over GF(2).
THREAD C: Schedule interaction with layers — quantify.
"""

import random, math
from collections import Counter


# ═══════════════════════════════════════════════════
# THREAD A: Analytical proof of T3 (Nilpotency)
# ═══════════════════════════════════════════════════

def prove_nilpotency():
    """
    THEOREM 3: C_y^n(x) = 0 for all x, y ∈ {0,1}^n.

    PROOF:
    C_y(x)[k] = carry into position k when computing x + y.
    C_y(x)[0] = 0 always.
    C_y(x)[k] = MAJ(x[k-1], y[k-1], C_y(x)[k-1])
              = x[k-1]·y[k-1] + x[k-1]·C[k-1] + y[k-1]·C[k-1]

    Key observation: C_y(x)[k] can be 1 ONLY IF at least 2 of
    {x[k-1], y[k-1], C[k-1]} are 1.

    Now apply C_y AGAIN: C_y(C_y(x)).
    Input = C_y(x), which has C_y(x)[0] = 0.
    So C_y(C_y(x))[1] = MAJ(C_y(x)[0], y[0], 0) = MAJ(0, y[0], 0) = 0.
    Regardless of y[0]!

    After 1st application: position 0 = 0.
    After 2nd application: positions 0,1 = 0.
    After kth application: positions 0..k-1 = 0.
    After nth application: ALL positions = 0.

    WHY? Because C_y(x)[k] depends on C[k-1] (the previous carry).
    If C[k-1] = 0, then C[k] = x[k-1] AND y[k-1] (just generate).
    The generate pattern of (C_y(x), y) is SPARSER than (x, y)
    because C_y(x) has more zeros (position 0 = 0 guaranteed).

    More precisely: after m applications, positions 0..m-1 are 0.
    At position m: C_y^m(x)[m] = MAJ(C_y^{m-1}(x)[m-1], y[m-1], C_y^m(x)[m-1])
                                = MAJ(0, y[m-1], 0)  [since positions 0..m-1 = 0]
                                = 0.

    Wait — that proves positions clear one at a time from the bottom!
    After m applications: pos 0..m-1 = 0, pos m = MAJ(0, y[m-1], 0) = 0.
    So after m+1 applications: pos 0..m = 0.
    By induction: after n applications, all n positions = 0.

    Actually: C_y(x)[0] = 0. So after 1 application, only pos 0 = 0.
    C_y(v)[1] where v[0] = 0: C_y(v)[1] = MAJ(v[0], y[0], C_y(v)[0]) = MAJ(0, y[0], 0) = 0.
    C_y(v)[2] = MAJ(v[1], y[1], C_y(v)[1]) = MAJ(v[1], y[1], 0) = v[1] AND y[1].
    This is NOT necessarily 0!

    Hmm. The argument breaks at position 2 because v[1] can be nonzero.
    Let me reconsider.

    After 1 application: v = C_y(x). v[0] = 0.
    After 2 applications: w = C_y(v). w[0] = 0, w[1] = MAJ(v[0], y[0], w[0]) = 0.
    After 3 applications: u = C_y(w). u[0] = 0, u[1] = 0,
        u[2] = MAJ(w[1], y[1], u[1]) = MAJ(0, y[1], 0) = 0.
    After k applications: positions 0..k-1 = 0.
    After n applications: ALL = 0.  QED!

    The key: each application zeros out ONE MORE position from the bottom,
    because that position's carry depends on the (already zero) position below.
    """
    print("=" * 80)
    print("THREAD A: Analytical proof of Nilpotency (T3)")
    print("=" * 80)

    print("""
    THEOREM 3 (Carry Nilpotency): C_y^n(x) = 0 for all x, y in {0,1}^n.

    PROOF BY INDUCTION on k (number of applications):

    Claim: After k applications of C_y, positions 0..k-1 are all 0.

    Base: k=1. C_y(x)[0] = 0 by definition (no carry into bit 0).

    Step: Assume after k applications, v = C_y^k(x) has v[0]=...=v[k-1]=0.
          Apply C_y once more: w = C_y(v).
          w[0] = 0 (always).
          For j = 1..k:
            w[j] = MAJ(v[j-1], y[j-1], w[j-1])
            v[j-1] = 0 (by induction, j-1 < k)
            w[j-1] = 0 (by this same argument, j-1 < k)
            → w[j] = MAJ(0, y[j-1], 0) = 0.

          So w[0] = ... = w[k] = 0. ✓

    After n applications: all n positions = 0.  ∎
    """)

    # Verify
    n = 8
    N = 2**n
    confirmed = True
    for y in range(N):
        for x in range(N):
            v = x
            for step in range(n):
                # Apply C_y
                c = 0
                new_v = 0
                for k in range(n):
                    new_v |= (c << k)
                    vk = (v >> k) & 1
                    yk = (y >> k) & 1
                    c = (vk & yk) | (vk & c) | (yk & c)
                v = new_v
            if v != 0:
                confirmed = False
                break
        if not confirmed:
            break

    print(f"    Exhaustive verification n=8: C_y^8(x) = 0 for all x,y: {confirmed}")


# ═══════════════════════════════════════════════════
# THREAD B: BTE Spectrum
# ═══════════════════════════════════════════════════

def bte_spectrum():
    """
    The round function F over GF(2) has a Jacobian J (256×256).
    Eigenvalues of J over GF(2): λ ∈ {0, 1}.
    Eigenvalue 1 multiplicity = dim of invariant subspace.
    Eigenvalue 0 multiplicity = dim of "killed" subspace.

    We found earlier: eigenvalue-1 dim varies 0-3 depending on state.
    But the CHARACTERISTIC POLYNOMIAL of J might be universal.

    For the BTE CLASS: what is the characteristic polynomial of
    the "skeleton" Jacobian (the 93.7% fixed part)?
    """
    print("\n" + "=" * 80)
    print("THREAD B: BTE Spectrum — characteristic polynomial of skeleton")
    print("=" * 80)

    # The skeleton = Jacobian with all carry-dependent entries set to their
    # MOST COMMON value (0 for entries that are usually 0, 1 for usually 1).
    # Actually: the "always-0" and "always-1" entries = the skeleton.
    # Variable entries → set to 0 (since they depend on state).

    # From bte_star.py: 61387 entries always same, 4149 variable.
    # The skeleton has rank 256 (always invertible).

    # Over GF(2): characteristic polynomial of 256×256 matrix.
    # det(J - λI) over GF(2). λ = 0 or 1.
    # det(J) = 1 if invertible (which it is — rank 256).
    # det(J - I) = 0 if eigenvalue 1 exists.

    # We already measured: eigenvalue-1 dim = 0 or 1 depending on state.
    # For skeleton (fixed part): let's compute.

    # But computing det of 256×256 is expensive symbolically.
    # Instead: compute rank(J - I) for the skeleton.

    # We need the skeleton matrix. Build it from the "always" values.
    # From bte_star.py experiment: we collected always_one and always_zero.
    # Let's rebuild quickly.

    from stage0_observe import sha256_round_trace, MASK, H0, K, Sig0, Sig1, Ch, Maj, schedule

    N = 50  # Fewer samples for speed (we already know 93.7% is fixed)

    always_one = [[True]*256 for _ in range(256)]
    always_zero = [[True]*256 for _ in range(256)]

    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        states, W = sha256_round_trace(M)

        r = 20
        a,b,c,d,e,f,g,h = states[r]
        T1 = (h + Sig1(e) + Ch(e,f,g) + K[r] + W[r]) & MASK
        T2 = (Sig0(a) + Maj(a,b,c)) & MASK
        base_out = [(T1+T2)&MASK, a, b, c, (d+T1)&MASK, e, f, g]

        base_bits = []
        for reg in base_out:
            for k in range(32):
                base_bits.append((reg >> k) & 1)

        for ri in range(8):
            for bi in range(32):
                ts = list(states[r])
                ts[ri] ^= (1 << bi)
                ta,tb,tc,td,te,tf,tg,th = ts
                tT1 = (th+Sig1(te)+Ch(te,tf,tg)+K[r]+W[r])&MASK
                tT2 = (Sig0(ta)+Maj(ta,tb,tc))&MASK
                tout = [(tT1+tT2)&MASK, ta, tb, tc, (td+tT1)&MASK, te, tf, tg]
                tbits = []
                for reg in tout:
                    for k in range(32):
                        tbits.append((reg >> k) & 1)

                idx = ri*32 + bi
                for j in range(256):
                    val = base_bits[j] ^ tbits[j]
                    if val == 0:
                        always_one[idx][j] = False
                    if val == 1:
                        always_zero[idx][j] = False

    # Build skeleton: always_one = 1, always_zero = 0, variable = 0
    skeleton = [[0]*256 for _ in range(256)]
    n_fixed = 0
    for i in range(256):
        for j in range(256):
            if always_one[i][j]:
                skeleton[i][j] = 1
                n_fixed += 1
            elif always_zero[i][j]:
                n_fixed += 1
            # else: variable → 0

    print(f"  Skeleton: {n_fixed}/{256*256} fixed entries")

    # Rank of skeleton
    m = [list(row) for row in skeleton]
    rank = 0
    for col in range(256):
        pivot = None
        for row in range(rank, 256):
            if m[row][col] == 1:
                pivot = row
                break
        if pivot is None:
            continue
        m[rank], m[pivot] = m[pivot], m[rank]
        for row in range(256):
            if row != rank and m[row][col] == 1:
                for c in range(256):
                    m[row][c] ^= m[rank][c]
        rank += 1

    print(f"  Skeleton rank: {rank}/256")

    # Rank of (skeleton - I)
    skel_I = [[skeleton[i][j] ^ (1 if i == j else 0) for j in range(256)] for i in range(256)]
    m2 = [list(row) for row in skel_I]
    rank2 = 0
    for col in range(256):
        pivot = None
        for row in range(rank2, 256):
            if m2[row][col] == 1:
                pivot = row
                break
        if pivot is None:
            continue
        m2[rank2], m2[pivot] = m2[pivot], m2[rank2]
        for row in range(256):
            if row != rank2 and m2[row][col] == 1:
                for c in range(256):
                    m2[row][c] ^= m2[rank2][c]
        rank2 += 1

    eig1_dim = 256 - rank2
    print(f"  Skeleton rank(J-I): {rank2}/256")
    print(f"  Skeleton eigenvalue-1 dim: {eig1_dim}")

    if eig1_dim > 0:
        print(f"  *** SKELETON HAS {eig1_dim} INVARIANT DIRECTIONS ***")
    else:
        print(f"  Skeleton has no eigenvalue-1 (fully mixing)")


# ═══════════════════════════════════════════════════
# THREAD C: Schedule × Layers interaction
# ═══════════════════════════════════════════════════

def schedule_layers():
    """
    The schedule W[16..63] = f(W[0..15]) is GF(2)-linear.
    How does it interact with the 4-layer structure?

    Schedule operates on WORDS (32 bits). Layer structure operates on BITS.
    The schedule links word i to words i-2, i-7, i-15, i-16.
    Each word has 32 bit-layers.

    Question: does the schedule PRESERVE layer structure?
    I.e., does schedule at bit k only involve bit k of input words?
    """
    print("\n" + "=" * 80)
    print("THREAD C: Schedule × Layer interaction")
    print("=" * 80)

    # Schedule: W[i] = sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16]
    # sig0(x) = ROTR(x,7) XOR ROTR(x,18) XOR (x >> 3)
    # sig1(x) = ROTR(x,17) XOR ROTR(x,19) XOR (x >> 10)

    # sig0(x)[k] = x[(k+7)%32] XOR x[(k+18)%32] XOR x[k+3] (if k+3 < 32, else 0)
    # This reads from bit positions (k+7)%32, (k+18)%32, k+3.
    # These are DIFFERENT layers!

    # So: schedule at bit k reads from layers (k+7)%32, (k+18)%32, k+3
    # of one word, plus bit k of three other words.

    # The schedule COUPLES layers through sig0 and sig1 rotations.
    # Specifically:
    #   sig0 couples layer k to layers (k+7)%32, (k+18)%32, k+3
    #   sig1 couples layer k to layers (k+17)%32, (k+19)%32, k+10

    print(f"  Schedule sig0(x)[k] reads from layers: (k+7)%32, (k+18)%32, k+3")
    print(f"  Schedule sig1(x)[k] reads from layers: (k+17)%32, (k+19)%32, k+10")
    print(f"  + bit k of W[i-7] and W[i-16] (same layer)")
    print(f"")

    # For layer 0: schedule reads from layers 7, 18, 3 (sig0) and 17, 19, 10 (sig1)
    # Plus layer 0 of W[i-7] and W[i-16].
    print(f"  Layer 0 schedule dependencies:")
    sig0_layers = [(0+7)%32, (0+18)%32, min(0+3, 31)]
    sig1_layers = [(0+17)%32, (0+19)%32, min(0+10, 31)]
    print(f"    From sig0: layers {sig0_layers}")
    print(f"    From sig1: layers {sig1_layers}")
    print(f"    From W[i-7], W[i-16]: layer 0 (same layer)")
    print(f"")

    # KEY: the schedule at layer 0 reads from layers 3, 7, 10, 17, 18, 19.
    # This means: layer 0 of W[16+] depends on OTHER layers of W[0..15].
    # The schedule BREAKS layer independence!

    # BUT: the schedule is about W values, not about state.
    # The state trajectory at layer 0 doesn't directly use schedule layer coupling.
    # State depends on W[r], and at layer 0: only W[r][0] matters (no carry).
    # W[r][0] = sig1(W[r-2])[0] + W[r-7][0] + sig0(W[r-15])[0] + W[r-16][0]
    # And sig0(W[r-15])[0] = W[r-15][7] XOR W[r-15][18] XOR W[r-15][3]

    # So: W[r][0] involves W[r-15] at bits 3, 7, 18 (different layers!)
    # This IS the layer coupling through schedule.

    # But on bit-0 level: the + in schedule = XOR at bit 0 (no carry).
    # So: W[r][0] = W[r-2][17] XOR W[r-2][19] XOR W[r-7][0] XOR W[r-15][7] XOR W[r-15][18] XOR ...

    # Wait: the schedule has ADDITION (not XOR). At bit 0:
    # (a + b)[0] = a[0] XOR b[0] (no carry). So yes, XOR at bit 0.

    # But sig0(W[r-15]) involves rotation which is XOR of bits at different positions.
    # sig0(x)[0] = x[7] XOR x[18] XOR (x >> 3)[0] = x[7] XOR x[18] XOR x[3]
    # These are bits 3, 7, 18 of x. NOT bit 0.

    # So: the layer-0 schedule equation for W[r][0] involves bits of OTHER layers.
    # This means: our layer-0 system is NOT closed in W-space.
    # It WOULD be closed if schedule only involved same-layer bits.
    # But sig0/sig1 in schedule couple different layers.

    # CONSEQUENCE: The 127-rank layer-0 system implicitly includes
    # constraints from other layers (through schedule sig0/sig1).
    # The "layer independence" we proved (MI=0) on OUTPUT is consistent
    # because schedule coupling is LINEAR and doesn't create bias.

    print(f"  CONSEQUENCE:")
    print(f"  Schedule couples layers through sig0/sig1 rotations.")
    print(f"  Layer 0 of W[16+] depends on layers 3,7,10,17,18,19 of earlier W.")
    print(f"  BUT: all coupling is LINEAR (XOR at bit 0, rotations = XOR).")
    print(f"  Linear coupling → no bias → MI = 0 on output. Consistent.")
    print(f"")
    print(f"  The schedule makes layers CONNECTED but not CORRELATED.")
    print(f"  Connection = structural (through linear dependencies).")
    print(f"  Correlation = statistical (through bias). None found.")
    print(f"")
    print(f"  INSIGHT: Schedule is the ONLY source of cross-layer coupling")
    print(f"  that spans all 64 rounds. Within round: Sig0/Sig1 couple layers.")
    print(f"  Across rounds: schedule couples layers additionally.")
    print(f"  Both are LINEAR → no exploitable correlation.")


if __name__ == "__main__":
    prove_nilpotency()
    bte_spectrum()
    schedule_layers()
