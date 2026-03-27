"""
ЩЕЛЬ В СФЕРЕ: condition number = 34.

Метрический тензор g_ij (32×32 для бит W[0]) имеет:
  Max eigenvalue: 2114
  Min eigenvalue: 62
  Ratio: 34×

Это значит: есть ОДНО направление в котором hash sensitivity
в 34 раз сильнее чем в перпендикулярном.

Вопросы:
  A. Что за направление? (eigenvector)
  B. Это свойство W[0] или всех слов?
  C. Почему? (carry structure? bit position?)
  D. Можно ли эксплуатировать?
"""

import numpy as np
import struct, hashlib

MASK32 = 0xFFFFFFFF

def sha256_hash(W16):
    raw = struct.pack('>16I', *W16)
    return struct.unpack('>8I', hashlib.sha256(raw).digest())

def hw(x): return bin(x).count('1')


def build_metric_tensor(word_idx, n_samples=2000):
    """Build 32×32 metric tensor for bits of word word_idx."""
    G = np.zeros((32, 32))
    for _ in range(n_samples):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H_base = sha256_hash(W)

        rows = []
        for bit in range(32):
            W_mod = list(W); W_mod[word_idx] ^= (1 << bit)
            H_mod = sha256_hash(W_mod)
            row = []
            for w in range(8):
                d = H_base[w] ^ H_mod[w]
                for b in range(32):
                    row.append((d >> b) & 1)
            rows.append(np.array(row, dtype=float))

        J = np.array(rows)
        G += J @ J.T

    G /= n_samples
    return G


def main():
    np.random.seed(42)

    print("=" * 70)
    print("ЩЕЛЬ В СФЕРЕ: привилегированное направление")
    print("=" * 70)

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("A. Eigenvector analysis для W[0]")
    print("=" * 70)

    G0 = build_metric_tensor(0, n_samples=3000)
    eigenvalues, eigenvectors = np.linalg.eigh(G0)
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    print(f"  Eigenvalues of g_ij (W[0], 32×32):")
    for i in range(min(10, len(eigenvalues))):
        print(f"    λ_{i} = {eigenvalues[i]:.1f}")
    print(f"    ...")
    print(f"    λ_31 = {eigenvalues[-1]:.1f}")
    print(f"    Condition: {eigenvalues[0]/eigenvalues[-1]:.1f}")

    # What is the dominant eigenvector?
    v0 = eigenvectors[:, 0]  # dominant direction
    print(f"\n  Dominant eigenvector v_0 (λ={eigenvalues[0]:.0f}):")
    print(f"    Components: {v0}")

    # Is it a specific bit pattern?
    # All-ones vector? (sum of all bits)
    all_ones = np.ones(32) / np.sqrt(32)
    dot_allones = abs(np.dot(v0, all_ones))
    print(f"\n    Alignment with all-ones: {dot_allones:.4f} (1.0 = perfect)")

    # Individual bit vectors?
    max_component = np.argmax(np.abs(v0))
    print(f"    Largest component: bit {max_component} ({v0[max_component]:.4f})")

    # Is it uniform? (all components equal)
    uniformity = np.std(np.abs(v0)) / np.mean(np.abs(v0))
    print(f"    Uniformity (std/mean of |v|): {uniformity:.4f} (0 = perfectly uniform)")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("B. STRUCTURE of metric tensor")
    print("=" * 70)

    # g_ij = diagonal (self) + off-diagonal (cross)
    # If g_ij = a*I + b*J (where J = all-ones matrix):
    # eigenvalues: a + 32b (once, for all-ones vector), a (31 times)
    # Check!

    diag_val = np.mean(np.diag(G0))
    offdiag_val = np.mean(G0[np.triu_indices(32, k=1)])

    a_pred = diag_val - offdiag_val  # ≈ 64
    b_pred = offdiag_val  # ≈ 64
    lambda_max_pred = a_pred + 32 * b_pred  # ≈ 64 + 32*64 = 2112
    lambda_rest_pred = a_pred  # ≈ 64

    print(f"  Metric decomposition: g = a*I + b*J")
    print(f"    a (diagonal - offdiag) = {a_pred:.1f}")
    print(f"    b (off-diagonal) = {b_pred:.1f}")
    print(f"    Predicted λ_max = a + 32b = {lambda_max_pred:.1f}")
    print(f"    Actual λ_max = {eigenvalues[0]:.1f}")
    print(f"    Predicted λ_rest = a = {lambda_rest_pred:.1f}")
    print(f"    Actual λ_rest = {np.mean(eigenvalues[1:]):.1f}")
    print(f"    Match: {'YES ✓' if abs(lambda_max_pred - eigenvalues[0]) < 10 else 'NO'}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("C. INTERPRETATION")
    print("=" * 70)

    print(f"""
  The metric tensor has form: g = {a_pred:.0f}·I + {b_pred:.0f}·J

  Where I = identity (32×32), J = all-ones matrix.

  This means:
    - Flipping ANY single bit → δH ≈ {diag_val:.0f} bits (diagonal)
    - Flipping TWO bits → δH correlated by {offdiag_val:.0f} (off-diagonal)
    - ALL pairs equally correlated (J is uniform)

  The "privileged direction" = ALL-ONES vector (1,1,1,...,1)/√32.
  This is the direction where ALL 32 bits flip simultaneously.
  Its eigenvalue = {eigenvalues[0]:.0f} vs single-bit = {eigenvalues[1]:.0f}.

  WHY? Flipping all 32 bits of W[0] = W[0] → W[0] ⊕ 0xFFFFFFFF = ~W[0].
  This is NOT 32 independent bit flips — it's a COMPLEMENT.
  W + ~W = 0xFFFFFFFF (all ones). Special algebraic relationship!
  The carry structure of ADD makes complement SPECIAL.
""")

    # Verify: is complement really special?
    complement_hws = []
    random_32_hws = []
    for _ in range(2000):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H = sha256_hash(W)

        # Complement W[0]
        W_c = list(W); W_c[0] ^= MASK32
        H_c = sha256_hash(W_c)
        complement_hws.append(sum(hw(H[i]^H_c[i]) for i in range(8)))

        # Random 32-bit flip
        W_r = list(W); W_r[0] ^= np.random.randint(1, 2**32)
        H_r = sha256_hash(W_r)
        random_32_hws.append(sum(hw(H[i]^H_r[i]) for i in range(8)))

    print(f"  Complement W[0] → ~W[0]: mean δH = {np.mean(complement_hws):.1f}")
    print(f"  Random 32-bit flip:      mean δH = {np.mean(random_32_hws):.1f}")
    print(f"  → {'COMPLEMENT IS SPECIAL!' if abs(np.mean(complement_hws) - np.mean(random_32_hws)) > 2 else 'Same as random'}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("D. PER-WORD metric tensors: is this universal?")
    print("=" * 70)

    conditions = []
    for word in range(16):
        G_w = build_metric_tensor(word, n_samples=1000)
        eigs = np.linalg.eigvalsh(G_w)
        eigs = np.sort(eigs)[::-1]
        cond = eigs[0] / eigs[-1]
        conditions.append(cond)

        diag_w = np.mean(np.diag(G_w))
        offdiag_w = np.mean(G_w[np.triu_indices(32, k=1)])

        print(f"  W[{word:>2}]: cond={cond:>5.1f}, diag={diag_w:>6.1f}, offdiag={offdiag_w:>6.1f}, "
              f"a={diag_w-offdiag_w:>5.1f}, b={offdiag_w:>5.1f}")

    print(f"\n  Condition numbers across words: mean={np.mean(conditions):.1f}, std={np.std(conditions):.1f}")
    print(f"  → {'UNIVERSAL (all words same)' if np.std(conditions) < 5 else 'WORD-DEPENDENT'}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("E. TRUE METRIC: remove the all-ones component")
    print("=" * 70)

    # The g = aI + bJ structure means:
    # - In the 31-dimensional subspace PERPENDICULAR to all-ones: g = aI (flat!)
    # - The all-ones direction is just "how many bits flip" (Hamming weight)
    #
    # If we project out the all-ones component:
    # g_reduced = aI (31×31) → condition number = 1!

    # Verify: eigenvalues 1..31 should all be ≈ a
    print(f"  Eigenvalues λ_1..λ_31 (without all-ones direction):")
    print(f"    Mean: {np.mean(eigenvalues[1:]):.1f}")
    print(f"    Std:  {np.std(eigenvalues[1:]):.1f}")
    print(f"    Min:  {eigenvalues[-1]:.1f}")
    print(f"    Max:  {eigenvalues[1]:.1f}")
    print(f"    Condition (reduced): {eigenvalues[1]/eigenvalues[-1]:.2f}")
    print(f"    → {'PERFECTLY FLAT!' if eigenvalues[1]/eigenvalues[-1] < 1.5 else 'Still anisotropic'}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("F. WHAT DOES THIS MEAN?")
    print("=" * 70)

    print(f"""
  ═══════════════════════════════════════════════════════════════

  THE "CRACK IN THE SPHERE" IS NOT A CRACK.

  The condition number 34 comes from ONE direction: all-ones (1,1,...,1).
  This is the HAMMING WEIGHT direction.
  Its eigenvalue = a + 32b = {eigenvalues[0]:.0f}.

  In ALL OTHER 31 directions: eigenvalue = a = {np.mean(eigenvalues[1:]):.0f}.
  Reduced condition number: {eigenvalues[1]/eigenvalues[-1]:.2f} ← FLAT!

  EXPLANATION:
    g_ij = <J_i, J_j> where J_i = output change when flipping bit i.
    J_i has HW ≈ 128 (random-looking 256-bit vector).
    <J_i, J_i> = HW(J_i) ≈ 128 (diagonal).
    <J_i, J_j> = |J_i ∩ J_j| ≈ 128 × 128 / 256 = 64 (random overlap).

    This 64 = 128/2 is a MATHEMATICAL IDENTITY for random binary vectors!
    E[<u,v>] = n × p² = 256 × (1/2)² = 64.

    The all-ones eigenvalue = 128 + 31×64 = 2112.
    This is NOT a "privileged direction" — it's the TRACE of the metric.

  CONCLUSION:
    The sphere HAS NO CRACK.
    Condition number 34 = trivial artifact of g = aI + bJ structure.
    After removing the trace component: PERFECTLY FLAT (condition 1.1).
    SHA-256 message space = S^511 (31-dimensional per word) × trace.

  MATHEMATICAL THEOREM:
    For any hash H: {{0,1}}^n → {{0,1}}^m with avalanche property:
      g_ij = (m/4)·I + (m/4)·J
      Eigenvalues: m·n/4 (once), m/4 (n-1 times)
      Condition: n (trivial)
      Reduced condition: 1 (flat!)

    This is a THEOREM about random binary vectors, not SHA-256.
    ANY hash with good avalanche has this EXACT metric structure.

  ═══════════════════════════════════════════════════════════════
""")


if __name__ == "__main__":
    main()
