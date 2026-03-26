#!/usr/bin/env python3
"""
SCF: Kernel + Carry Suppression в Schedule

Факт: ker(Schedule_GF2[52..63]) dim=128, даёт δW_xor=0 в хвосте.
Факт: carry добавляет ~13 бит/слово → δW_real ≠ 0.
Факт: schedule carry зависит ТОЛЬКО от message (нет state!).

Идея: для КОНКРЕТНОГО W_base, найти элемент ядра, где carry
в schedule[52..63] ТОЖЕ мал. Carry контролируем полностью
потому что знаем все операнды schedule.

Стратегия:
1. Извлечь ядро (dim=128 для R=52)
2. Для конкретного W_base: вычислить carry каждого базисного вектора
3. SA по 2^128 комбинациям → минимизировать РЕАЛЬНЫЙ tail HW
4. Перебирать W_base → найти пару (W_base, δW) с δW_real[52..63]≈0
"""
import os, sys, math

MASK32 = 0xFFFFFFFF

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&MASK32
def sig0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def sig1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def add32(*args):
    s=0
    for x in args: s=(s+x)&MASK32
    return s
def hw(x): return bin(x).count('1')

def expand_xor(W16):
    W=list(W16)
    for i in range(16,64):
        W.append(sig1(W[i-2])^W[i-7]^sig0(W[i-15])^W[i-16])
    return W

def expand_real(W16):
    W=list(W16)
    for i in range(16,64):
        W.append(add32(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    return W

def bits_to_words(bits):
    words = []
    for w in range(16):
        val = 0
        for b in range(32):
            if bits[w*32+b]:
                val |= (1 << b)
        words.append(val)
    return words


def extract_kernel_basis(R_start=52):
    n_out_bits = (64 - R_start) * 32
    n_in_bits = 512

    M = []
    for out_word in range(R_start, 64):
        for out_bit in range(32):
            row = [0] * n_in_bits
            for in_word in range(16):
                for in_bit in range(32):
                    W_test = [0]*16
                    W_test[in_word] = 1 << in_bit
                    Wexp = expand_xor(W_test)
                    if (Wexp[out_word] >> out_bit) & 1:
                        row[in_word*32 + in_bit] = 1
            M.append(row)

    m = [list(row) for row in M]
    n_rows, n_cols = len(m), n_in_bits
    pivot_cols = []
    row_idx = 0
    for col in range(n_cols):
        pivot = -1
        for r in range(row_idx, n_rows):
            if m[r][col] == 1:
                pivot = r; break
        if pivot == -1: continue
        m[row_idx], m[pivot] = m[pivot], m[row_idx]
        for r in range(n_rows):
            if r != row_idx and m[r][col] == 1:
                for c in range(n_cols): m[r][c] ^= m[row_idx][c]
        pivot_cols.append(col)
        row_idx += 1

    free_cols = [c for c in range(n_cols) if c not in pivot_cols]
    kernel_basis = []
    for fc in free_cols:
        vec = [0]*n_cols
        vec[fc] = 1
        for i, pc in enumerate(pivot_cols):
            if m[i][fc] == 1: vec[pc] = 1
        kernel_basis.append(vec)

    return kernel_basis


def real_tail_hw(dW_words, W_base, R_start=52):
    W_mod = [(W_base[i] ^ dW_words[i]) for i in range(16)]
    We1 = expand_real(W_base)
    We2 = expand_real(W_mod)
    return sum(hw(We1[r] ^ We2[r]) for r in range(R_start, 64))

def real_tail_per_word(dW_words, W_base, R_start=52):
    W_mod = [(W_base[i] ^ dW_words[i]) for i in range(16)]
    We1 = expand_real(W_base)
    We2 = expand_real(W_mod)
    return [hw(We1[r] ^ We2[r]) for r in range(R_start, 64)]


# ============================================================
# EXP 1: Joint optimization (W_base, kernel_combo)
# ============================================================
def exp1_joint_optimize(basis, R_start=52, N_outer=50, N_inner=10000):
    """For each W_base candidate, find best kernel combination."""
    print("="*70)
    print(f"EXP 1: JOINT (W_base × kernel) OPTIMIZATION")
    print(f"  Search {N_outer} base messages × {N_inner} kernel combos each")
    print("="*70)

    n_basis = len(basis)
    global_best_score = 9999
    global_best_W = None
    global_best_dW = None

    for outer in range(N_outer):
        W_base = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

        # SA over kernel combinations
        current_bits = list(basis[0])
        current_words = bits_to_words(current_bits)
        current_score = real_tail_hw(current_words, W_base, R_start)

        best_bits = list(current_bits)
        best_score = current_score

        for it in range(N_inner):
            trial_bits = list(current_bits)
            idx = int.from_bytes(os.urandom(2),'big') % n_basis
            for j in range(512):
                trial_bits[j] ^= basis[idx][j]
            if all(b==0 for b in trial_bits): continue

            trial_words = bits_to_words(trial_bits)
            trial_score = real_tail_hw(trial_words, W_base, R_start)

            T = max(0.01, 1.0 - it/N_inner)
            if trial_score < current_score or \
               math.exp(-(trial_score-current_score)/(T*2)) > \
               int.from_bytes(os.urandom(4),'big')/(1<<32):
                current_bits = trial_bits
                current_score = trial_score

            if current_score < best_score:
                best_score = current_score
                best_bits = list(current_bits)

        if best_score < global_best_score:
            global_best_score = best_score
            global_best_W = list(W_base)
            global_best_dW = bits_to_words(best_bits)
            print(f"  [{outer:3d}] New best: tail_HW = {global_best_score}")

    # Analyze best result
    print(f"\n  BEST RESULT: real_tail_HW = {global_best_score}")
    per_word = real_tail_per_word(global_best_dW, global_best_W, R_start)
    zeros = sum(1 for h in per_word if h == 0)
    print(f"  Per-word profile (r={R_start}..63):")
    for i, h in enumerate(per_word):
        r = R_start + i
        marker = " ★ ZERO!" if h == 0 else (" ☆ low" if h <= 5 else "")
        print(f"    r={r}: HW={h:3d}{marker}")
    print(f"  Zero words in real tail: {zeros}/{64-R_start}")

    return global_best_W, global_best_dW, global_best_score


# ============================================================
# EXP 2: Carry-aware basis scoring
# ============================================================
def exp2_carry_basis_analysis(basis, R_start=52, N_bases=20):
    """For multiple W_base, score each basis vector.
    Find: are some basis vectors carry-quiet across many bases?"""
    print("\n" + "="*70)
    print("EXP 2: CARRY-AWARE BASIS ANALYSIS")
    print("  Which basis vectors are carry-quiet across MANY W_base?")
    print("="*70)

    n = len(basis)
    # Score matrix: [basis_idx][w_base_idx] → real_tail_hw
    avg_scores = [0.0] * n

    for b_idx in range(N_bases):
        W_base = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        for v_idx in range(n):
            words = bits_to_words(basis[v_idx])
            sc = real_tail_hw(words, W_base, R_start)
            avg_scores[v_idx] += sc

    for v_idx in range(n):
        avg_scores[v_idx] /= N_bases

    # Sort by average
    ranked = sorted(range(n), key=lambda i: avg_scores[i])

    print(f"\n  Top 10 carry-quiet basis vectors (avg over {N_bases} bases):")
    for rank, idx in enumerate(ranked[:10]):
        print(f"    basis[{idx:3d}]: avg_tail_HW = {avg_scores[idx]:.1f}")

    print(f"\n  Bottom 3:")
    for idx in ranked[-3:]:
        print(f"    basis[{idx:3d}]: avg_tail_HW = {avg_scores[idx]:.1f}")

    overall_avg = sum(avg_scores)/n
    print(f"\n  Overall average: {overall_avg:.1f}")
    print(f"  Best single basis: {avg_scores[ranked[0]]:.1f}")
    print(f"  Spread: {avg_scores[ranked[-1]] - avg_scores[ranked[0]]:.1f}")

    return ranked, avg_scores


# ============================================================
# EXP 3: Focused search — use top basis vectors only
# ============================================================
def exp3_focused_search(basis, ranked, R_start=52, N_outer=100, N_inner=5000):
    """Use only top-K carry-quiet basis vectors for combination."""
    print("\n" + "="*70)
    print("EXP 3: FOCUSED SEARCH — top carry-quiet basis vectors")
    print("="*70)

    for top_k in [16, 32, 64]:
        top_basis = [basis[ranked[i]] for i in range(min(top_k, len(ranked)))]
        n_top = len(top_basis)

        best_score = 9999
        best_W = None

        for outer in range(N_outer):
            W_base = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

            current_bits = list(top_basis[0])
            current_score = real_tail_hw(bits_to_words(current_bits), W_base, R_start)
            local_best = current_score

            for it in range(N_inner):
                trial_bits = list(current_bits)
                idx = int.from_bytes(os.urandom(2),'big') % n_top
                for j in range(512):
                    trial_bits[j] ^= top_basis[idx][j]
                if all(b==0 for b in trial_bits): continue

                sc = real_tail_hw(bits_to_words(trial_bits), W_base, R_start)

                T = max(0.01, 1.0 - it/N_inner)
                if sc < current_score or \
                   math.exp(-(sc-current_score)/(T*2)) > \
                   int.from_bytes(os.urandom(4),'big')/(1<<32):
                    current_bits = trial_bits
                    current_score = sc
                if current_score < local_best:
                    local_best = current_score

            if local_best < best_score:
                best_score = local_best
                best_W = list(W_base)

        print(f"  Top-{top_k}: best real_tail_HW = {best_score}")


# ============================================================
# EXP 4: THE KEY — can we find δW_real = 0 for ANY word?
# ============================================================
def exp4_per_word_zero(basis, R_start=52, N_outer=200, N_inner=3000):
    """Can we make even ONE word of δW_real = 0 in the tail?"""
    print("\n" + "="*70)
    print("EXP 4: CAN WE ZERO EVEN ONE REAL SCHEDULE WORD?")
    print(f"  For r in [{R_start}..63]: find (W_base, kernel_combo) with δW_real[r]=0")
    print("="*70)

    n = len(basis)

    for target_r in range(63, R_start-1, -1):
        best_hw = 32

        for outer in range(N_outer):
            W_base = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

            current_bits = list(basis[outer % n])
            W_mod = [(W_base[i] ^ bits_to_words(current_bits)[i]) for i in range(16)]
            We1 = expand_real(W_base)
            We2 = expand_real(W_mod)
            current_hw_r = hw(We1[target_r] ^ We2[target_r])

            for it in range(N_inner):
                trial_bits = list(current_bits)
                idx = int.from_bytes(os.urandom(2),'big') % n
                for j in range(512): trial_bits[j] ^= basis[idx][j]
                if all(b==0 for b in trial_bits): continue

                tw = bits_to_words(trial_bits)
                Wm = [(W_base[i]^tw[i]) for i in range(16)]
                We2t = expand_real(Wm)
                hw_r = hw(We1[target_r] ^ We2t[target_r])

                if hw_r < current_hw_r:
                    current_bits = trial_bits
                    current_hw_r = hw_r

            if current_hw_r < best_hw:
                best_hw = current_hw_r

            if best_hw == 0:
                break

        marker = " ★★★ ZERO FOUND!" if best_hw == 0 else (" ★ low!" if best_hw <= 3 else "")
        print(f"  r={target_r}: best HW(δW_real) = {best_hw}{marker}")

        if best_hw == 0:
            print(f"    → δW_real[{target_r}] = 0 IS ACHIEVABLE through kernel!")


# ============================================================
if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 100

    print("="*70)
    print("SCF: KERNEL + CARRY SUPPRESSION IN SCHEDULE")
    print("="*70)

    R = 52
    print(f"\nExtracting kernel basis (R={R})...")
    basis = extract_kernel_basis(R)
    print(f"  Kernel dim: {len(basis)}")

    ranked, avg_scores = exp2_carry_basis_analysis(basis, R, N_bases=min(N//5+1, 20))
    exp3_focused_search(basis, ranked, R, N_outer=min(N, 50), N_inner=3000)
    exp1_joint_optimize(basis, R, N_outer=min(N, 30), N_inner=min(N*50, 10000))
    exp4_per_word_zero(basis, R, N_outer=min(N, 100), N_inner=2000)

    print("\n" + "="*70)
    print("ИТОГ")
    print("="*70)
    print("""
  Вопрос: можно ли обнулить РЕАЛЬНЫЙ schedule tail через ядро?

  GF(2)-ядро гарантирует δW_xor[52..63] = 0.
  Carry-correction δcc = δW_real делает tail ненулевым.
  Carry зависит от КОНКРЕТНОГО W_base.

  Если найдена пара (W_base, δW∈kernel) с δW_real[r]=0:
  → backward chain получает бесплатный ноль на раунде r-3
  → gap между forward и backward СУЖАЕТСЯ

  Если δW_real[r]=0 для нескольких последовательных r:
  → backward chain работает на этом участке
  → формируется мост от forward Wang к backward chain
    """)
