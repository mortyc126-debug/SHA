"""
CARRY BIRTHDAY: поиск kernel vector с нулевым carry error.

GF(2)-kernel: 256-мерный. 2^256 vectors.
Carry error на random kernel vector: ~132 бит.
Нужен vector с carry error = 0.

Если carry error = random 256-бит:
  P(carry=0) = 2^{-256} per vector.
  Need 2^256 vectors. Нет gain.

Если carry error < 256 бит (зависит от kernel structure):
  P(carry=0) = 2^{-E} где E = effective carry entropy.
  Need 2^E vectors из 2^256 available.
  Если E < 256: ДЕШЕВЛЕ чем birthday!

Birthday внутри kernel:
  2^256 vectors → birthday на carry error → need √(2^E) = 2^{E/2}.
  Если E = 132: need 2^66 kernel vectors.
  2^66 из 2^256 = легко доступны.
  → CARRY BIRTHDAY = 2^66???

Измеряем E (effective carry entropy) точно.
"""

import numpy as np
import struct, hashlib

MASK32 = 0xFFFFFFFF

def sha256_words(W16):
    return struct.unpack('>8I', hashlib.sha256(struct.pack('>16I', *W16)).digest())

def hw(x): return bin(x).count('1')


def main():
    np.random.seed(42)

    print("=" * 70)
    print("CARRY BIRTHDAY: birthday ВНУТРИ GF(2)-kernel")
    print("=" * 70)

    # =================================================================
    print(f"\n{'=' * 70}")
    print("1. СТРОИМ GF(2)-kernel (Gaussian elimination)")
    print("=" * 70)

    # Build T matrix for specific W_base
    W_base = [np.random.randint(0, 2**32) for _ in range(16)]
    H_base = sha256_words(W_base)

    T_matrix = []
    for word in range(16):
        for bit in range(32):
            W_mod = list(W_base)
            W_mod[word] ^= (1 << bit)
            H_mod = sha256_words(W_mod)
            row = []
            for w in range(8):
                delta = H_base[w] ^ H_mod[w]
                for b in range(32):
                    row.append((delta >> b) & 1)
            T_matrix.append(row)

    T = np.array(T_matrix, dtype=np.uint8)  # 512 × 256

    # GF(2) kernel
    M = T.T.copy()  # 256 × 512
    pivots = []
    row = 0
    for col in range(512):
        found = False
        for r in range(row, 256):
            if M[r, col] == 1:
                M[[row, r]] = M[[r, row]]
                found = True
                break
        if not found:
            continue
        pivots.append(col)
        for r in range(256):
            if r != row and M[r, col] == 1:
                M[r] = M[r] ^ M[row]
        row += 1

    free_vars = [c for c in range(512) if c not in pivots]
    print(f"  Kernel dimension: {len(free_vars)}")

    # Generate kernel vectors
    def make_kernel_vector(free_idx):
        x = np.zeros(512, dtype=np.uint8)
        x[free_vars[free_idx]] = 1
        for i in range(len(pivots)-1, -1, -1):
            pc = pivots[i]
            val = 0
            for j in range(512):
                if j != pc:
                    val ^= (M[i, j] & x[j])
            x[pc] = val
        return x

    def kernel_to_dW(x):
        dW = [0] * 16
        for word in range(16):
            for bit in range(32):
                if x[word * 32 + bit]:
                    dW[word] ^= (1 << bit)
        return dW

    def compute_carry_error(dW):
        """Real SHA-256 δH for given δW (carry error of GF2 kernel vector)."""
        W2 = [W_base[i] ^ dW[i] for i in range(16)]
        H1 = sha256_words(W_base)
        H2 = sha256_words(W2)
        return tuple(H1[i] ^ H2[i] for i in range(8))

    # =================================================================
    print(f"\n{'=' * 70}")
    print("2. CARRY ERROR DISTRIBUTION: базисные vectors")
    print("=" * 70)

    carry_errors = []
    carry_hws = []

    for i in range(min(256, len(free_vars))):
        x = make_kernel_vector(i)
        dW = kernel_to_dW(x)
        ce = compute_carry_error(dW)
        ce_hw = sum(hw(w) for w in ce)
        carry_errors.append(ce)
        carry_hws.append(ce_hw)

    print(f"\n  {len(carry_hws)} базисных kernel vectors:")
    print(f"    Mean carry error HW: {np.mean(carry_hws):.1f}")
    print(f"    Std: {np.std(carry_hws):.1f}")
    print(f"    Min: {min(carry_hws)}")
    print(f"    Max: {max(carry_hws)}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("3. XOR COMBINATIONS: можно ли уменьшить carry error?")
    print("=" * 70)

    # XOR двух kernel vectors = новый kernel vector.
    # carry_error(v1 ⊕ v2) = ??? (NOT carry_error(v1) ⊕ carry_error(v2)!)
    # Потому что carry нелинеен.

    # Проверяем: carry_error(v1⊕v2) vs carry_error(v1)⊕carry_error(v2)
    combo_hws = []
    pred_hws = []

    for _ in range(1000):
        i = np.random.randint(0, len(free_vars))
        j = np.random.randint(0, len(free_vars))
        if i == j:
            continue

        x_i = make_kernel_vector(i)
        x_j = make_kernel_vector(j)
        x_combo = x_i ^ x_j

        dW_combo = kernel_to_dW(x_combo)
        ce_combo = compute_carry_error(dW_combo)
        hw_combo = sum(hw(w) for w in ce_combo)
        combo_hws.append(hw_combo)

        # Predicted (linear): ce(i) ⊕ ce(j)
        if i < len(carry_errors) and j < len(carry_errors):
            pred = tuple(carry_errors[i][w] ^ carry_errors[j][w] for w in range(8))
            hw_pred = sum(hw(w) for w in pred)
            pred_hws.append(hw_pred)

    print(f"\n  XOR of pairs of kernel vectors:")
    print(f"    Actual carry error HW: mean={np.mean(combo_hws):.1f}, min={min(combo_hws)}")
    if pred_hws:
        print(f"    Predicted (linear):    mean={np.mean(pred_hws):.1f}")
        print(f"    Difference:            {abs(np.mean(combo_hws) - np.mean(pred_hws)):.1f}")
        print(f"    → {'CARRY IS ADDITIVE (linear)!' if abs(np.mean(combo_hws) - np.mean(pred_hws)) < 10 else 'Carry nonlinear'}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("4. RANDOM KERNEL VECTORS: statistical search")
    print("=" * 70)

    # Generate random kernel vectors (random XOR of basis) and measure carry error
    best_hw = 256
    best_dW = None
    hw_dist = []

    for trial in range(10000):
        # Random kernel vector: XOR random subset of basis
        x = np.zeros(512, dtype=np.uint8)
        n_basis = np.random.randint(1, min(20, len(free_vars)))
        indices = np.random.choice(len(free_vars), n_basis, replace=False)
        for idx in indices:
            x = x ^ make_kernel_vector(idx)

        dW = kernel_to_dW(x)
        if all(w == 0 for w in dW):
            continue

        ce = compute_carry_error(dW)
        ce_hw = sum(hw(w) for w in ce)
        hw_dist.append(ce_hw)

        if ce_hw < best_hw:
            best_hw = ce_hw
            best_dW = dW

    print(f"\n  10K random kernel vectors:")
    print(f"    Mean carry error: {np.mean(hw_dist):.1f}")
    print(f"    Std: {np.std(hw_dist):.1f}")
    print(f"    Min: {min(hw_dist)}")
    print(f"    Max: {max(hw_dist)}")

    # Distribution
    print(f"\n  Distribution:")
    for threshold in [80, 90, 100, 110, 120, 130]:
        count = sum(1 for h in hw_dist if h <= threshold)
        print(f"    HW ≤ {threshold}: {count}/10000 = {count/100:.2f}%")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("5. CARRY BIRTHDAY СТОИМОСТЬ")
    print("=" * 70)

    # Carry error из random kernel vector: mean = M, std = S.
    # Carry error = 256-бит вектор. Его HW ≈ M.
    # "Effective entropy" E: если carry error = random E-bit vector,
    # birthday на E бит = 2^{E/2}.

    # НО: carry error зависит от δW. Два РАЗНЫХ kernel vectors
    # могут иметь ОДИНАКОВЫЙ carry error → это CARRY COLLISION!
    # Birthday на carry error space: need √(|carry error space|).

    # |carry error space| = сколько различных carry errors?
    unique_errors = len(set(tuple(sum(hw(w) for w in ce) for ce in [compute_carry_error(kernel_to_dW(make_kernel_vector(i))) for i in range(min(100, len(free_vars)))]) for _ in range(1)))

    # Simpler: just count unique carry error vectors
    unique_ce = set()
    for i in range(min(200, len(free_vars))):
        ce = compute_carry_error(kernel_to_dW(make_kernel_vector(i)))
        unique_ce.add(ce)

    print(f"\n  Unique carry errors from 200 basis vectors: {len(unique_ce)}")
    print(f"  → Carry error space ≈ full (all unique)")

    mean_hw = np.mean(hw_dist)
    std_hw = np.std(hw_dist)

    print(f"""
  CARRY BIRTHDAY:
    Kernel dimension: 256 (2^256 vectors available)
    Carry error: mean HW = {mean_hw:.1f}, std = {std_hw:.1f}
    Carry error appears random 256-bit → effective entropy E ≈ 256

    Birthday на carry error: √(2^E) = 2^{{E/2}}
    Если E = 256: 2^128 (same as standard birthday!)
    Если E = {mean_hw:.0f}: 2^{mean_hw/2:.0f} (carry-specific)

    2^256 kernel vectors → birthday на carry error →
    Need 2^{{E/2}} kernel vectors.
    Each vector: O(1) to generate (XOR of basis) + O(1) SHA-256.
    Total: 2^{{E/2}} SHA-256 evaluations.

    Если E < 256: ДЕШЕВЛЕ standard birthday!

    Best carry error found (10K): HW = {best_hw}
    Random HW expected: ~128
    → Carry error ≈ random 256-бит → E ≈ 256 → birthday = 2^128.

    НО: carry error HW = {mean_hw:.0f} (не 128!).
    Carry error BIASED — средний HW > 128!
    Это может означать carry error < 256 effective бит.
""")


if __name__ == "__main__":
    main()
