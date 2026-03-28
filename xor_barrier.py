#!/usr/bin/env python3
"""
XOR Barrier Attack — Axis 7
==============================

RAYON: "XOR has no controlling value → constant propagation fails"
Carry-web: "schedule linear over GF(2) → XOR structure opaque"

Both stop at the same wall. Attack it from 3 angles:

  X1: PARTIAL XOR CANCELLATION — which ΔW[0..15] give δW[r]=0 (XOR)
      for maximum number of schedule words r=16..63?
      GF(2) linear → solvable by Gaussian elimination.
      More zeros in schedule → more carry-free rounds → more RAYON cutoffs.

  X2: SCHEDULE WEIGHT MINIMIZATION — find ΔW[0..15] with minimum
      total HW(δW[16..63]) over XOR schedule.
      Low HW = fewer active bits = easier differential.
      This is a coding theory problem (minimum weight codeword).

  X3: COMBINED XOR+CARRY — XOR schedule gives δW[r] (exact over GF(2)).
      Carry correction = δW_real[r] - δW_xor[r] (nonlinear part).
      If carry correction is SMALL for specific δW[0..15] patterns →
      XOR analysis is approximately correct → RAYON can use it.
"""

import numpy as np
import time, sys

MASK = 0xFFFFFFFF

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&MASK
def sig0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def sig1(x): return rotr(x,17)^rotr(x,19)^(x>>10)

def hw(x): return bin(x&MASK).count('1')


# ============================================================
# X1: PARTIAL XOR CANCELLATION
# ============================================================

def experiment_X1(seed=21000):
    print("="*70)
    print("X1: PARTIAL XOR CANCELLATION — Maximize zeros in XOR schedule")
    print("="*70)

    # XOR schedule: δW[r] = sig1(δW[r-2]) ⊕ δW[r-7] ⊕ sig0(δW[r-15]) ⊕ δW[r-16]
    # This is LINEAR over GF(2)^32.
    # Build 1536×512 matrix (48 output words × 32 bits each, 16 input words × 32 bits)
    # Find input that maximizes number of zero output words.

    print(f"\n  Building GF(2) schedule matrix (1536×512)...")
    t0 = time.time()

    # sig0 and sig1 as 32×32 GF(2) matrices
    def func_to_matrix(f):
        """Convert 32-bit → 32-bit function to GF(2) matrix."""
        M = np.zeros((32, 32), dtype=np.int8)
        for b in range(32):
            out = f(1 << b)
            for ob in range(32):
                M[ob, b] = (out >> ob) & 1
        return M

    S0 = func_to_matrix(sig0)
    S1 = func_to_matrix(sig1)
    I32 = np.eye(32, dtype=np.int8)
    Z32 = np.zeros((32, 32), dtype=np.int8)

    # Schedule recurrence: W[i] = sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16]
    # In XOR: W[i] = sig1(W[i-2]) ⊕ W[i-7] ⊕ sig0(W[i-15]) ⊕ W[i-16]

    # Build expansion: for each output word W[16..63], express as function of W[0..15]
    # W[i] = Σ M_j × W[j] over GF(2), where M_j are 32×32 matrices

    # Dependency matrix per output word: coeff[i][j] = 32×32 GF(2) matrix
    # such that W[i] = ⊕_j coeff[i][j] × W[j]
    coeff = {}
    for j in range(16):
        coeff[(j, j)] = I32.copy()

    for i in range(16, 64):
        # W[i] = sig1(W[i-2]) ⊕ W[i-7] ⊕ sig0(W[i-15]) ⊕ W[i-16]
        for j in range(16):
            M = Z32.copy()
            if (i-2, j) in coeff:
                M = (M + S1 @ coeff[(i-2, j)]) % 2
            if (i-7, j) in coeff:
                M = (M + coeff[(i-7, j)]) % 2
            if (i-15, j) in coeff:
                M = (M + S0 @ coeff[(i-15, j)]) % 2
            if (i-16, j) in coeff:
                M = (M + coeff[(i-16, j)]) % 2
            if np.any(M):
                coeff[(i, j)] = M

    print(f"  Matrix built: {time.time()-t0:.1f}s")

    # For each schedule word i=16..63: which input words affect it?
    print(f"\n  --- Schedule word dependencies ---")
    for i in [16, 17, 18, 19, 20, 30, 47, 48, 63]:
        deps = [j for j in range(16) if (i, j) in coeff]
        hw_total = sum(np.sum(coeff[(i,j)]) for j in deps)
        print(f"    W[{i:2d}]: depends on {len(deps)} words, total GF2 weight = {hw_total}")

    # Find ΔW[0..15] that zeros maximum schedule words
    # Strategy: for each target word i, solve coeff[i] × ΔW = 0
    # The intersection of kernels gives ΔW that zeros ALL targets

    # Build full 1536×512 matrix
    L = np.zeros((48*32, 16*32), dtype=np.int8)
    for i in range(16, 64):
        row_base = (i-16)*32
        for j in range(16):
            col_base = j*32
            if (i, j) in coeff:
                L[row_base:row_base+32, col_base:col_base+32] = coeff[(i, j)]

    rank = np.linalg.matrix_rank(L.astype(float))
    print(f"\n  Full schedule matrix: {L.shape}, rank={rank}/512")
    print(f"  Kernel dimension: {512 - rank}")

    if rank < 512:
        print(f"  ★ Kernel exists! There are {2**(512-rank)} nonzero ΔW with some δW[16..63]=0")
    else:
        print(f"  Full rank → only ΔW=0 gives all zeros (expected)")

    # Find ΔW that zeros SPECIFIC words (e.g., W[17], W[19], W[21])
    # From methodology: δW[17]=δW[19]=δW[21]=0 when δW[0]=1, W[1..15]=0
    print(f"\n  --- Specific zero targets ---")
    for targets in [[17], [17,19], [17,19,21], [17,19,21,23]]:
        rows = []
        for t in targets:
            rows.extend(range((t-16)*32, (t-16)*32+32))
        L_sub = L[rows, :]
        rank_sub = np.linalg.matrix_rank(L_sub.astype(float))
        kernel_dim = 512 - rank_sub
        print(f"    Targets {targets}: rank={rank_sub}, kernel={kernel_dim} → 2^{kernel_dim} solutions")

    # HW distribution of schedule output for random nonzero ΔW
    print(f"\n  --- Schedule output HW for random ΔW ---")
    rng = np.random.RandomState(seed)
    hw_totals = []
    zero_word_counts = []
    for _ in range(1000):
        dw = rng.randint(0, 2, size=512).astype(np.int8)
        if not np.any(dw): continue
        out = L @ dw % 2
        hw_total = np.sum(out)
        n_zero_words = sum(1 for i in range(48) if not np.any(out[i*32:(i+1)*32]))
        hw_totals.append(hw_total)
        zero_word_counts.append(n_zero_words)

    print(f"    E[HW(δW[16..63])]: {np.mean(hw_totals):.1f} / 1536")
    print(f"    E[zero words]:     {np.mean(zero_word_counts):.2f} / 48")
    print(f"    Max zero words:    {max(zero_word_counts)}")
    print(f"    Expected random:   HW≈768, zero words≈0")

    return L, coeff


# ============================================================
# X2: MINIMUM WEIGHT CODEWORD
# ============================================================

def experiment_X2(L, seed=21001):
    print(f"\n{'='*70}")
    print("X2: MINIMUM WEIGHT CODEWORD — Lightest nonzero schedule output")
    print("="*70)

    # The set {L×x : x ∈ GF(2)^512} is a linear code.
    # Minimum distance = min HW of nonzero codeword.
    # If min_d is small → sparse differential through schedule.

    rng = np.random.RandomState(seed)

    # Strategy 1: random search
    print(f"\n  --- Random search (10K inputs) ---")
    min_hw = 1536
    min_zeros = 0
    best_input = None

    for trial in range(10000):
        # Sparse random input (few active bits)
        n_active = rng.randint(1, 33)
        dw = np.zeros(512, dtype=np.int8)
        active_bits = rng.choice(512, size=n_active, replace=False)
        dw[active_bits] = 1

        out = L @ dw % 2
        hw_out = int(np.sum(out))
        n_zeros = sum(1 for i in range(48) if not np.any(out[i*32:(i+1)*32]))

        if hw_out < min_hw:
            min_hw = hw_out
            min_zeros = n_zeros
            best_input = dw.copy()

    print(f"    Min HW(output): {min_hw} / 1536")
    print(f"    Zero words at min: {min_zeros} / 48")
    print(f"    Input HW: {int(np.sum(best_input))}")

    # Strategy 2: single-bit inputs (16×32 = 512 possibilities)
    print(f"\n  --- Single-bit ΔW inputs ---")
    single_hws = []
    single_zeros = []
    for b in range(512):
        dw = np.zeros(512, dtype=np.int8)
        dw[b] = 1
        out = L @ dw % 2
        hw_out = int(np.sum(out))
        n_zeros = sum(1 for i in range(48) if not np.any(out[i*32:(i+1)*32]))
        single_hws.append(hw_out)
        single_zeros.append(n_zeros)

    best_b = np.argmin(single_hws)
    word_idx = best_b // 32
    bit_idx = best_b % 32
    print(f"    Best single bit: W[{word_idx}][b{bit_idx}]")
    print(f"    HW(output): {single_hws[best_b]} / 1536")
    print(f"    Zero words: {single_zeros[best_b]} / 48")

    # Top-5 sparsest single bits
    top5 = sorted(range(512), key=lambda b: single_hws[b])[:5]
    print(f"\n    Top-5 sparsest single-bit ΔW:")
    for b in top5:
        w, bi = b//32, b%32
        print(f"      W[{w}][b{bi}]: HW={single_hws[b]}, zeros={single_zeros[b]}")

    # Statistics
    print(f"\n    Single-bit HW: mean={np.mean(single_hws):.1f}, min={min(single_hws)}, max={max(single_hws)}")
    print(f"    Single-bit zeros: mean={np.mean(single_zeros):.1f}, max={max(single_zeros)}")

    # Which schedule words are ZERO for the best single bit?
    dw_best = np.zeros(512, dtype=np.int8)
    dw_best[top5[0]] = 1
    out_best = L @ dw_best % 2
    zero_rounds = [i+16 for i in range(48) if not np.any(out_best[i*32:(i+1)*32])]
    print(f"\n    Zero schedule words for best bit: {zero_rounds}")

    return min_hw, single_hws


# ============================================================
# X3: XOR+CARRY COMBINED
# ============================================================

def experiment_X3(N=3000, seed=21002):
    print(f"\n{'='*70}")
    print("X3: XOR+CARRY COMBINED — How close is XOR to real schedule?")
    print(f"N={N}")
    print("="*70)

    rng = np.random.RandomState(seed)

    def schedule_real(W16):
        W=list(W16)+[0]*48
        for i in range(16,64): W[i]=(sig1(W[i-2])+W[i-7]+sig0(W[i-15])+W[i-16])&MASK
        return W

    def schedule_xor(W16):
        W=list(W16)+[0]*48
        for i in range(16,64): W[i]=sig1(W[i-2])^W[i-7]^sig0(W[i-15])^W[i-16]
        return W

    # For random messages: HW(W_real[r] XOR W_xor[r]) measures carry correction
    print(f"\n  --- Carry correction per schedule word ---")
    carry_hw = np.zeros(48)

    for _ in range(N):
        W16 = [int(rng.randint(0, 1<<32)) for _ in range(16)]
        Wr = schedule_real(W16)
        Wx = schedule_xor(W16)
        for r in range(16, 64):
            carry_hw[r-16] += hw(Wr[r] ^ Wx[r])

    carry_hw /= N
    print(f"  {'r':>3} | {'E[HW(carry_corr)]':>18} | {'fraction':>9}")
    print(f"  {'-'*3}-+-{'-'*18}-+-{'-'*9}")

    for r_off in [0, 1, 2, 3, 4, 10, 20, 30, 40, 47]:
        r = r_off + 16
        frac = carry_hw[r_off] / 32
        marker = " ★" if carry_hw[r_off] < 8 else ""
        print(f"  {r:3d} | {carry_hw[r_off]:18.2f} | {frac:9.3f}{marker}")

    print(f"\n  Mean carry correction: {np.mean(carry_hw):.2f} / 32")
    print(f"  This is {np.mean(carry_hw)/32*100:.1f}% of bits — {'small' if np.mean(carry_hw)<10 else 'large'}")

    # For DIFFERENTIAL: δW_real = δW_xor + carry_diff
    print(f"\n  --- Differential carry correction ---")
    print(f"  For M, M' with δW[0]=1 (1-bit diff): how close is δW_real to δW_xor?")

    diff_hw = np.zeros(48)
    for _ in range(N):
        W16 = [int(rng.randint(0, 1<<32)) for _ in range(16)]
        W16p = list(W16); W16p[0] ^= 1

        Wr1 = schedule_real(W16); Wr2 = schedule_real(W16p)
        Wx1 = schedule_xor(W16); Wx2 = schedule_xor(W16p)

        for r in range(16, 64):
            dw_real = Wr1[r] ^ Wr2[r]  # XOR difference in real schedule
            dw_xor = Wx1[r] ^ Wx2[r]   # XOR difference in XOR schedule
            diff_hw[r-16] += hw(dw_real ^ dw_xor)

    diff_hw /= N
    print(f"\n  {'r':>3} | {'E[HW(δW_real ⊕ δW_xor)]':>25} | {'match%':>7}")
    print(f"  {'-'*3}-+-{'-'*25}-+-{'-'*7}")

    near_zero = 0
    for r_off in range(48):
        r = r_off + 16
        match = 1 - diff_hw[r_off]/32
        if diff_hw[r_off] < 0.01:
            near_zero += 1
        if r_off < 10 or diff_hw[r_off] < 5 or r_off >= 45:
            marker = " ★★" if diff_hw[r_off] < 0.01 else (" ★" if diff_hw[r_off] < 5 else "")
            print(f"  {r:3d} | {diff_hw[r_off]:25.2f} | {match*100:6.1f}%{marker}")

    print(f"\n  Schedule words where XOR ≈ real diff (HW < 0.01): {near_zero}/48")
    print(f"  Mean differential correction: {np.mean(diff_hw):.2f} bits")

    if near_zero > 5:
        print(f"\n  ★★★ {near_zero} schedule words are EXACTLY XOR-predictable!")
        print(f"  For these words: RAYON can use XOR model with zero error.")
        print(f"  Carry correction = 0 → constant propagation applies perfectly.")
    elif np.mean(diff_hw) < 8:
        print(f"\n  ★ XOR model is approximate (mean error {np.mean(diff_hw):.1f} bits)")
    else:
        print(f"\n  XOR model inadequate (mean error {np.mean(diff_hw):.1f} bits)")

    return carry_hw, diff_hw


# ============================================================
# MAIN
# ============================================================

if __name__ == '__main__':
    print("XOR Barrier Attack — Axis 7")
    print(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")

    N = 3000
    if len(sys.argv) > 1 and sys.argv[1] == '--fast':
        N = 1500

    t_start = time.time()
    L, coeff = experiment_X1()
    min_hw, single_hws = experiment_X2(L)
    carry_hw, diff_hw = experiment_X3(N=N)

    total = time.time() - t_start
    print(f"\n{'='*70}")
    print(f"TOTAL TIME: {total:.1f}s")
    print(f"{'='*70}")

    print(f"""
{'='*70}
SYNTHESIS: XOR Barrier Attack (Axis 7)
{'='*70}

X1: SCHEDULE MATRIX
  Full rank = {np.linalg.matrix_rank(L.astype(float))}/512 → no universal zero patterns
  But specific target zeros possible (kernel of submatrices)

X2: MINIMUM WEIGHT
  Best single-bit: HW = {min(single_hws)}/1536
  Sparsest schedule differential → easiest carry correction

X3: XOR ≈ REAL?
  Near-zero correction rounds: {sum(1 for d in diff_hw if d < 0.01)}/48
  Mean carry correction: {np.mean(carry_hw):.1f} bits per word

KEY FOR RAYON:
  If carry correction is 0 for specific rounds → RAYON DFS works
  on those rounds. XOR-predictable rounds = RAYON-attackable rounds.
""")
