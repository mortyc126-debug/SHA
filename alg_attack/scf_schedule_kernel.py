#!/usr/bin/env python3
"""
SCF: Schedule Kernel Attack — соединяем forward Wang и backward chain
через GF(2)-ядро расписания.

Schedule: W[i] = σ1(W[i-2]) + W[i-7] + σ0(W[i-15]) + W[i-16]
Над GF(2) это ЛИНЕЙНОЕ отображение (512 → 1536 бит).

Идея:
1. Вычислить GF(2)-матрицу schedule: δW[0..15] → δW[16..63]
2. Найти δW[0..15] с минимальным HW(δW[48..63]) — "tail-sparse"
3. Forward Wang chain: δe[2..16]=0 (от выбора δW)
4. Backward chain: δe[r-3]=-δW[r], поэтому sparse tail → short backward gap
5. Оценить: сколько раундов перекрывается?
"""
import os, sys

MASK32 = 0xFFFFFFFF

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&MASK32
def sig0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def sig1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def add32(*args):
    s=0
    for x in args: s=(s+x)&MASK32
    return s
def hw(x): return bin(x).count('1')

K = [
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2,
]
IV = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
      0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]

def Sig0(x): return rotr(x,2)^rotr(x,13)^rotr(x,22)
def Sig1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def Ch(e,f,g): return (e&f)^(~e&g)&MASK32
def Maj(a,b,c): return (a&b)^(a&c)^(b&c)

def expand_xor(W16):
    """GF(2)-linear schedule (XOR only)."""
    W=list(W16)
    for i in range(16,64):
        W.append(sig1(W[i-2])^W[i-7]^sig0(W[i-15])^W[i-16])
    return W

def expand_real(W16):
    W=list(W16)
    for i in range(16,64):
        W.append(add32(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    return W


# ============================================================
# Build GF(2) schedule matrix (512 → 1536)
# ============================================================
def build_schedule_matrix():
    """Build the binary matrix M such that δW_xor[16..63] = M × δW[0..15].
    M is 1536×512 over GF(2)."""

    rows = []  # 1536 rows (48 words × 32 bits), each is 512-bit vector

    for out_word in range(16, 64):
        for out_bit in range(32):
            row = [0] * 512  # 16 words × 32 bits

            # Compute: flip each input bit, check this output bit
            for in_word in range(16):
                for in_bit in range(32):
                    W_zero = [0] * 16
                    W_zero[in_word] = 1 << in_bit
                    Wexp = expand_xor(W_zero)
                    if (Wexp[out_word] >> out_bit) & 1:
                        row[in_word * 32 + in_bit] = 1

            rows.append(row)

    return rows


def gf2_rank(matrix):
    m = [list(row) for row in matrix]
    nrows = len(m)
    ncols = len(m[0]) if m else 0
    rank = 0
    for col in range(ncols):
        pivot = -1
        for row in range(rank, nrows):
            if m[row][col] == 1:
                pivot = row
                break
        if pivot == -1:
            continue
        m[rank], m[pivot] = m[pivot], m[rank]
        for row in range(nrows):
            if row != rank and m[row][col] == 1:
                for c in range(ncols):
                    m[row][c] ^= m[rank][c]
        rank += 1
    return rank


# ============================================================
# EXP 1: Schedule matrix properties
# ============================================================
def exp1_schedule_matrix():
    print("="*70)
    print("EXP 1: GF(2) SCHEDULE MATRIX")
    print("="*70)

    M = build_schedule_matrix()
    print(f"  Matrix size: {len(M)} × {len(M[0])} (output bits × input bits)")

    # Rank of full matrix
    r = gf2_rank(M)
    print(f"  Rank: {r} / {min(len(M), len(M[0]))}")
    print(f"  Kernel dimension: {len(M[0]) - r}")

    # Rank of sub-matrices (tail of schedule)
    for start_word in [16, 32, 48, 56]:
        start_row = (start_word - 16) * 32
        sub = M[start_row:]
        sub_rank = gf2_rank(sub)
        n_rows = len(sub)
        print(f"  Rank of W[{start_word}..63] matrix: {sub_rank}/{min(n_rows,512)} ({n_rows} rows)")

    # Rank of JUST the tail (W[48..63])
    tail_start = (48 - 16) * 32  # row 1024
    tail = M[tail_start:]
    tail_rank = gf2_rank(tail)
    print(f"\n  ★ Rank of tail W[48..63]: {tail_rank}/512")
    print(f"    Kernel of tail: {512 - tail_rank} dimensions")
    if 512 - tail_rank > 0:
        print(f"    → {512 - tail_rank} free bits can be chosen without affecting W[48..63]!")

    return M


# ============================================================
# EXP 2: Find tail-sparse δW[0..15]
# ============================================================
def exp2_tail_sparse(N):
    print("\n" + "="*70)
    print("EXP 2: TAIL-SPARSE δW[0..15]")
    print("  Find δW[0..15] minimizing HW(δW[48..63]) in GF(2) schedule")
    print("="*70)

    import math

    # Random baseline
    baseline_tail_hw = []
    for _ in range(min(N, 1000)):
        dW = [0]*16
        dW[0] = int.from_bytes(os.urandom(4),'big')
        if dW[0] == 0: dW[0] = 1
        Wexp = expand_xor(dW)
        tail_hw = sum(hw(Wexp[r]) for r in range(48, 64))
        baseline_tail_hw.append(tail_hw)

    avg_base = sum(baseline_tail_hw)/len(baseline_tail_hw)
    min_base = min(baseline_tail_hw)
    print(f"\n  Baseline (random δW[0], δW[1..15]=0):")
    print(f"    Avg HW(δW[48..63]): {avg_base:.1f}")
    print(f"    Min HW(δW[48..63]): {min_base}")

    # SA search: optimize δW[0..15] for minimal tail
    best_dW = [0]*16
    best_dW[0] = 0x80000000  # Wang-style

    def tail_hw_score(dW):
        Wexp = expand_xor(dW)
        return sum(hw(Wexp[r]) for r in range(48, 64))

    best_score = tail_hw_score(best_dW)
    current_dW = list(best_dW)
    current_score = best_score

    for it in range(min(N*100, 50000)):
        trial = list(current_dW)
        word = int.from_bytes(os.urandom(1),'big') % 16
        bit = int.from_bytes(os.urandom(1),'big') % 32
        trial[word] ^= (1 << bit)

        # Keep at least one nonzero word
        if all(x == 0 for x in trial):
            continue

        sc = tail_hw_score(trial)

        T = max(0.01, 1.0 - it/(N*100))
        if sc < current_score or \
           math.exp(-(sc-current_score)/(T*5)) > int.from_bytes(os.urandom(4),'big')/(1<<32):
            current_dW = trial
            current_score = sc

        if current_score < best_score:
            best_score = current_score
            best_dW = list(current_dW)

    print(f"\n  SA-optimized δW[0..15]:")
    print(f"    Tail HW(δW[48..63]): {best_score}")
    print(f"    Reduction: {avg_base - best_score:.1f} bits ({(avg_base-best_score)/avg_base*100:.1f}%)")

    # Show schedule profile
    Wexp = expand_xor(best_dW)
    print(f"\n  Schedule HW profile (GF(2)):")
    print(f"  {'r':>3} | HW(δW[r]) | bar")
    for r in range(64):
        h = hw(Wexp[r])
        bar = "█" * min(h, 30)
        marker = " ←" if h == 0 else ""
        if r < 16 or r >= 48 or h == 0:
            print(f"  {r:3d} | {h:9d} | {bar}{marker}")

    zeros = sum(1 for r in range(16,64) if Wexp[r] == 0)
    tail_zeros = sum(1 for r in range(48,64) if Wexp[r] == 0)
    print(f"\n  Zero δW in schedule: {zeros}/48")
    print(f"  Zero δW in tail [48..63]: {tail_zeros}/16")

    # Also check: how sparse can we make W[56..63] (last 8)?
    def tail8_score(dW):
        Wexp = expand_xor(dW)
        return sum(hw(Wexp[r]) for r in range(56, 64))

    best8 = list(best_dW)
    best8_score = tail8_score(best8)
    current8 = list(best8)
    current8_score = best8_score

    for it in range(min(N*50, 30000)):
        trial = list(current8)
        word = int.from_bytes(os.urandom(1),'big') % 16
        bit = int.from_bytes(os.urandom(1),'big') % 32
        trial[word] ^= (1 << bit)
        if all(x == 0 for x in trial): continue

        sc = tail8_score(trial)
        T = max(0.01, 1.0 - it/(N*50))
        if sc < current8_score or \
           math.exp(-(sc-current8_score)/(T*3)) > int.from_bytes(os.urandom(4),'big')/(1<<32):
            current8 = trial
            current8_score = sc
        if current8_score < best8_score:
            best8_score = current8_score
            best8 = list(current8)

    print(f"\n  Last 8 rounds [56..63] optimization:")
    print(f"    Best HW(δW[56..63]): {best8_score}")

    return best_dW, best_score


# ============================================================
# EXP 3: Real schedule (with carry) — does carry help or hurt?
# ============================================================
def exp3_real_vs_xor_schedule(best_dW):
    print("\n" + "="*70)
    print("EXP 3: REAL vs XOR SCHEDULE — does carry change sparsity?")
    print("="*70)

    Wexp_xor = expand_xor(best_dW)
    # For real schedule, we need a base message
    W_base = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    W_modified = [(W_base[i] ^ best_dW[i]) for i in range(16)]

    Wexp_base = expand_real(W_base)
    Wexp_mod = expand_real(W_modified)

    print(f"\n  {'r':>3} | HW_xor  HW_real  diff")
    print("  " + "-"*40)
    total_xor = 0
    total_real = 0
    for r in range(48, 64):
        hw_xor = hw(Wexp_xor[r])
        hw_real = hw(Wexp_base[r] ^ Wexp_mod[r])
        total_xor += hw_xor
        total_real += hw_real
        diff = hw_real - hw_xor
        print(f"  {r:3d} | {hw_xor:7d} {hw_real:7d} {diff:+5d}")

    print(f"\n  Total tail [48..63]: XOR={total_xor}, Real={total_real}, diff={total_real-total_xor:+d}")


# ============================================================
# EXP 4: Combined attack estimate
# ============================================================
def exp4_combined_estimate(best_dW, best_tail_hw):
    print("\n" + "="*70)
    print("EXP 4: COMBINED ATTACK ESTIMATE")
    print("  Forward Wang (15 rounds) + Backward chain (sparse tail)")
    print("="*70)

    Wexp = expand_xor(best_dW)

    # Backward chain: δe[r-3] = -δW[r]
    # If HW(δW[r]) = h, then δe[r-3] has ~h nonzero bits
    # This δe propagates into Ch at next round, adding ~h/2 bits of noise

    # Count backward-controllable rounds
    backward_rounds = 0
    for r in range(63, 15, -1):
        if hw(Wexp[r]) == 0:
            backward_rounds += 1
        else:
            break

    # Forward Wang chain: 15 rounds typically (r=2..16)
    forward_rounds = 15

    gap = 64 - forward_rounds - backward_rounds
    print(f"\n  Forward Wang chain: {forward_rounds} rounds (r=2..16)")
    print(f"  Backward zero chain: {backward_rounds} rounds from end")
    print(f"  GAP: {gap} rounds")

    # Even with non-zero backward rounds, low-HW tail helps
    # Estimate: each nonzero δW[r] in backward creates ~HW(δW[r]) bits of noise
    # This noise compounds through Ch (degree 2)
    print(f"\n  Tail schedule profile (last 16 rounds):")
    backward_noise = 0
    for r in range(63, 47, -1):
        h = hw(Wexp[r])
        backward_noise += h
        status = "ZERO" if h == 0 else f"{h} bits noise"
        print(f"    r={r}: HW(δW)={h:3d} → {status}")

    print(f"\n  Total backward noise: {backward_noise} bits (vs random: ~256)")
    print(f"  Noise reduction: {(256-backward_noise)/256*100:.1f}%")

    # Attack cost estimate
    # Forward: Wang chain = O(1) (deterministic)
    # Barrier: birthday on δe[17] = O(2^32)
    # Backward: noise from schedule creates ~backward_noise bits of uncertainty
    # Each bit needs ~O(2^1) to resolve
    # Total additional cost: O(2^backward_noise) roughly

    print(f"\n  ATTACK COST ESTIMATE:")
    print(f"    Forward Wang chain: O(1)")
    print(f"    Barrier r=17: O(2^32)")
    print(f"    Backward noise resolution: ~O(2^{backward_noise//2}) (birthday)")
    print(f"    Classical birthday: O(2^128)")
    print(f"    Potential speedup: 2^{128 - 32 - backward_noise//2}")

    if 128 - 32 - backward_noise//2 > 0:
        print(f"\n  ★ Combined attack is FASTER than birthday by 2^{128-32-backward_noise//2}!")
    else:
        print(f"\n  Combined attack is NOT faster than birthday.")


# ============================================================
# EXP 5: Exhaustive search for schedule-zero tail
# ============================================================
def exp5_zero_tail_search():
    print("\n" + "="*70)
    print("EXP 5: SEARCH FOR δW WITH ZERO SCHEDULE TAIL")
    print("  Can any δW[0..15] give δW_xor[R..63]=0?")
    print("="*70)

    # For each starting round R, find if there's a δW making tail zero
    for R in [60, 56, 52, 48]:
        tail_size = 64 - R
        print(f"\n  Target: δW_xor[{R}..63] = 0 ({tail_size} words = {tail_size*32} bits)")

        # Build sub-matrix for just these output words
        sub_rows = []
        for out_word in range(R, 64):
            for out_bit in range(32):
                row = [0] * 512
                for in_word in range(16):
                    for in_bit in range(32):
                        W_zero = [0]*16
                        W_zero[in_word] = 1 << in_bit
                        Wexp = expand_xor(W_zero)
                        if (Wexp[out_word] >> out_bit) & 1:
                            row[in_word*32+in_bit] = 1
                sub_rows.append(row)

        r = gf2_rank(sub_rows)
        kernel_dim = 512 - r
        print(f"    Rank: {r}/512")
        print(f"    Kernel: {kernel_dim} dimensions")

        if kernel_dim > 0:
            print(f"    ★ {kernel_dim} free directions exist!")
            print(f"    These δW[0..15] produce ZERO δW_xor in rounds {R}..63")
        else:
            print(f"    No nontrivial solution (full rank)")


# ============================================================
if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 200

    print("="*70)
    print("SCF: SCHEDULE KERNEL ATTACK")
    print("  Forward Wang + GF(2) schedule kernel + backward chain")
    print("="*70)

    M = exp1_schedule_matrix()
    best_dW, best_tail_hw = exp2_tail_sparse(N)
    exp3_real_vs_xor_schedule(best_dW)
    exp4_combined_estimate(best_dW, best_tail_hw)
    exp5_zero_tail_search()

    print("\n" + "="*70)
    print("ФИНАЛЬНЫЙ ИТОГ")
    print("="*70)
    print("""
  Schedule kernel analysis:
  - GF(2) schedule is a LINEAR map 512 → 1536
  - Kernel of full schedule = 0 (rank 512, injective)
  - But TAIL sub-matrices may have kernels!
  - If kernel_dim(W[R..63]) > 0 for some R:
    → We can choose δW[0..15] that zeros the tail
    → Backward chain works for those rounds
    → Gap shrinks → attack improves

  Combined strategy:
    Forward Wang: 15 zero rounds (r=2..16) — FREE
    + Schedule kernel: zero or sparse tail — COMPUTABLE
    + Backward chain: zero rounds from end — CONDITIONAL on tail
    = Reduced gap → potentially below birthday bound
    """)
