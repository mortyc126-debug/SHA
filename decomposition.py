"""
THE DECOMPOSITION: SHA-256(M) = L(M) ⊕ Φ(M)

L(M) = XOR-SHA(M) — all additions replaced by XOR. GF(2)-linear.
Φ(M) = SHA-256(M) ⊕ L(M) — the carry cocycle map.

Properties of L:
  - GF(2)-linear (verified: PHI algebra)
  - Kernel dim = 256 (from PHI)
  - Rank = 256 (full, verified)

Properties of Φ:
  - Nonlinear
  - Composed of per-round cocycles (T5)
  - Each round's carry: nilpotent (T3), binomial rank (T4)
  - Φ maps 512 → 256 bits

This decomposition passes through ALL 64 rounds.
Both L and Φ are defined for any number of rounds.

COLLISION through decomposition:
  L(M1) ⊕ Φ(M1) = L(M2) ⊕ Φ(M2)
  → L(M1 ⊕ M2) = Φ(M1) ⊕ Φ(M2)

  Left: linear in δM (known for fixed δM).
  Right: "carry difference" between two messages.

Let's verify this decomposition and study Φ.
"""

import random
from stage0_observe import (
    sha256_round_trace, MASK, H0, K,
    rotr, Sig0, Sig1, Ch, Maj, schedule
)


def sha256_real(M, R=64):
    states, _ = sha256_round_trace(M, rounds=R)
    return [(states[R][i] + H0[i]) & MASK for i in range(8)]


def sha256_xor(M, R=64):
    """XOR-only SHA-256: all + replaced by XOR."""
    W = list(M)
    for i in range(16, R):
        s1 = rotr(W[i-2],17) ^ rotr(W[i-2],19) ^ (W[i-2]>>10)
        s0 = rotr(W[i-15],7) ^ rotr(W[i-15],18) ^ (W[i-15]>>3)
        W.append(s1 ^ W[i-7] ^ s0 ^ W[i-16])
    a,b,c,d,e,f,g,h = H0
    for r in range(R):
        T1 = h ^ Sig1(e) ^ Ch(e,f,g) ^ K[r] ^ W[r]
        T2 = Sig0(a) ^ Maj(a,b,c)
        h,g,f = g,f,e; e = d ^ T1; d,c,b = c,b,a; a = T1 ^ T2
    return [a^H0[0],b^H0[1],c^H0[2],d^H0[3],e^H0[4],f^H0[5],g^H0[6],h^H0[7]]


def phi_map(M, R=64):
    """Φ(M) = SHA-256(M) ⊕ L(M). The carry cocycle map."""
    H_real = sha256_real(M, R)
    H_xor = sha256_xor(M, R)
    return [H_real[i] ^ H_xor[i] for i in range(8)]


def experiment_decomposition_verify():
    """Verify: SHA-256(M) = L(M) XOR Φ(M) for all rounds."""
    print("=" * 80)
    print("DECOMPOSITION: SHA-256(M) = L(M) ⊕ Φ(M)")
    print("=" * 80)

    N = 1000
    for R in [4, 16, 64]:
        violations = 0
        for seed in range(N):
            rng = random.Random(seed)
            M = [rng.randint(0, MASK) for _ in range(16)]
            H = sha256_real(M, R)
            L = sha256_xor(M, R)
            P = phi_map(M, R)

            for i in range(8):
                if H[i] != (L[i] ^ P[i]):
                    violations += 1

        print(f"  R={R:>2}: SHA = L ⊕ Φ verified: {violations}/{N*8} violations")


def experiment_phi_properties():
    """Study properties of Φ — the carry cocycle map."""
    print("\n" + "=" * 80)
    print("Φ PROPERTIES: The carry cocycle map")
    print("=" * 80)

    N = 500

    # 1. Is Φ balanced? (HW ≈ 128 for random M)
    hws = []
    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        P = phi_map(M)
        hw = sum(bin(p).count('1') for p in P)
        hws.append(hw)

    avg_hw = sum(hws) / N
    print(f"\n  Balance: E[HW(Φ)] = {avg_hw:.1f}/256 (random = 128)")

    # 2. Is Φ XOR-linear? (Φ(M1⊕M2) = Φ(M1) ⊕ Φ(M2)?)
    linear_violations = 0
    total = 0
    for seed in range(200):
        rng = random.Random(seed)
        M1 = [rng.randint(0, MASK) for _ in range(16)]
        M2 = [rng.randint(0, MASK) for _ in range(16)]
        M12 = [M1[i] ^ M2[i] for i in range(16)]

        P1 = phi_map(M1)
        P2 = phi_map(M2)
        P12 = phi_map(M12)

        for i in range(8):
            if P12[i] != (P1[i] ^ P2[i]):
                linear_violations += 1
        total += 8

    print(f"  XOR-linearity: {linear_violations}/{total} violations → {'NONLINEAR' if linear_violations > 0 else 'LINEAR'}")

    # 3. Φ rank (GF(2))
    rng = random.Random(42)
    M_base = [rng.randint(0, MASK) for _ in range(16)]
    P_base = phi_map(M_base)

    J = []
    for word in range(16):
        for bit in range(32):
            M_flip = list(M_base)
            M_flip[word] ^= (1 << bit)
            P_flip = phi_map(M_flip)
            row = []
            for i in range(8):
                diff = P_base[i] ^ P_flip[i]
                for k in range(32):
                    row.append((diff >> k) & 1)
            J.append(row)

    # GF(2) rank
    m = [list(row) for row in J]
    rank = 0
    for col in range(256):
        pivot = None
        for row in range(rank, 512):
            if m[row][col] == 1:
                pivot = row; break
        if pivot is None: continue
        m[rank], m[pivot] = m[pivot], m[rank]
        for row in range(512):
            if row != rank and m[row][col] == 1:
                for c in range(256): m[row][c] ^= m[rank][c]
        rank += 1

    print(f"  Φ Jacobian rank: {rank}/256 (same as SHA-256)")

    # 4. COLLISION EQUATION: L(δM) = Φ(M) ⊕ Φ(M⊕δM)
    # For fixed δM: left side = constant. Right side = function of M.
    # How does Φ(M) ⊕ Φ(M⊕δM) behave?
    print(f"\n  COLLISION EQUATION: L(δM) = Φ(M) ⊕ Φ(M⊕δM)")

    delta_bit = 0  # δM = flip bit 0 of W[0]
    rng_base = random.Random(42)
    M0 = [rng_base.randint(0, MASK) for _ in range(16)]
    L_delta = sha256_xor([M0[i] ^ (1 if i == 0 else 0) for i in range(16)])
    L_delta_val = [sha256_xor(M0)[i] ^ L_delta[i] for i in range(8)]

    # Check: for many M, does Φ(M) ⊕ Φ(M⊕δM) = L(δM)?
    matches = 0
    total_check = 0
    phi_diffs = []
    for seed in range(N):
        rng = random.Random(seed)
        M = [rng.randint(0, MASK) for _ in range(16)]
        M_delta = list(M); M_delta[0] ^= 1

        P1 = phi_map(M)
        P2 = phi_map(M_delta)
        phi_diff = [P1[i] ^ P2[i] for i in range(8)]

        L_M = sha256_xor(M)
        L_Md = sha256_xor(M_delta)
        L_diff = [L_M[i] ^ L_Md[i] for i in range(8)]

        # For collision: phi_diff should equal L_diff
        if phi_diff == L_diff:
            matches += 1
        total_check += 1

        phi_diffs.append(phi_diff)

    print(f"  Φ(M)⊕Φ(M⊕δ) = L(δM) for {matches}/{total_check} messages")

    if matches == 0:
        # How far is Φ-diff from L-diff?
        hw_diffs = []
        for seed in range(N):
            rng = random.Random(seed)
            M = [rng.randint(0, MASK) for _ in range(16)]
            M_delta = list(M); M_delta[0] ^= 1
            P1 = phi_map(M); P2 = phi_map(M_delta)
            L_M = sha256_xor(M); L_Md = sha256_xor(M_delta)
            phi_diff = [P1[i]^P2[i] for i in range(8)]
            L_diff = [L_M[i]^L_Md[i] for i in range(8)]
            hw = sum(bin(phi_diff[i]^L_diff[i]).count('1') for i in range(8))
            hw_diffs.append(hw)

        avg_dist = sum(hw_diffs)/len(hw_diffs)
        print(f"  Average distance Φ-diff to L-diff: {avg_dist:.1f}/256 (collision needs 0)")
        print(f"  (random = 128)")

    # 5. Is L(δM) constant (independent of M) for fixed δM?
    # YES — because L is linear! L(M⊕δM) = L(M) ⊕ L(δM).
    # So L-diff = L(δM) for ALL M. This is the POWER of the decomposition.
    print(f"\n  L(δM) is CONSTANT for fixed δM (because L is linear).")
    print(f"  → Collision = find M where Φ(M)⊕Φ(M⊕δM) = L(δM) = known constant.")
    print(f"  → This is a PREIMAGE problem for the function Ψ(M) = Φ(M)⊕Φ(M⊕δM).")


def experiment_phi_per_round():
    """How does Φ grow round by round? Is it structured at early rounds?"""
    print("\n" + "=" * 80)
    print("Φ PER ROUND: How does the carry map evolve?")
    print("=" * 80)

    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]

    print(f"\n  HW(Φ(M)) at each round count:")
    for R in [1, 2, 4, 8, 12, 16, 24, 32, 48, 64]:
        P = phi_map(M, R)
        hw = sum(bin(p).count('1') for p in P)
        # Also: is Φ "closer to 0" at early rounds?
        print(f"    R={R:>2}: HW(Φ) = {hw:>3}/256")


if __name__ == "__main__":
    experiment_decomposition_verify()
    experiment_phi_properties()
    experiment_phi_per_round()
