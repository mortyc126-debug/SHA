#!/usr/bin/env python3
"""
SCF Задача B: Carry-корреляции в message schedule.

Вопрос: если carry(W[52])=0, какова P(carry(W[53])=0)?
Если P >> 2^{-16} → carry коррелированы → обнуление хвоста дешевле.

Schedule: W[i] = σ1(W[i-2]) + W[i-7] + σ0(W[i-15]) + W[i-16]
Carry: cc[i] = (σ1(W[i-2]) + W[i-7] + σ0(W[i-15]) + W[i-16]) ⊕
               (σ1(W[i-2]) ⊕ W[i-7] ⊕ σ0(W[i-15]) ⊕ W[i-16])
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

def psi4(a,b,c,d):
    """Carry correction for 4-operand addition."""
    return add32(a,b,c,d) ^ (a^b^c^d)

def schedule_carry_profile(W16):
    """Compute carry correction HW for each schedule word."""
    W = list(W16)
    cc = [0]*64  # carry correction HW
    for i in range(16, 64):
        ops = [sig1(W[i-2]), W[i-7], sig0(W[i-15]), W[i-16]]
        carry = psi4(*ops)
        cc[i] = hw(carry)
        W.append(add32(*ops))
    return cc, W

# ============================================================
# EXP 1: Pairwise carry correlation
# ============================================================
def exp1_pairwise(N):
    print("="*70)
    print("EXP 1: PAIRWISE CARRY CORRELATION IN SCHEDULE")
    print("  P(HW(cc[j])≤t | HW(cc[i])≤t) vs P(HW(cc[j])≤t)")
    print("="*70)

    # Collect carry profiles
    profiles = []
    for _ in range(N):
        W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        cc, _ = schedule_carry_profile(W)
        profiles.append(cc)

    # Threshold: HW ≤ t means "low carry"
    for t in [0, 3, 5, 8]:
        print(f"\n  --- Threshold: HW(cc) ≤ {t} ---")

        # Marginal P(HW(cc[i]) ≤ t)
        marginals = {}
        for i in range(52, 64):
            count = sum(1 for p in profiles if p[i] <= t)
            marginals[i] = count / N

        print(f"  Marginal P(cc[i]≤{t}):")
        for i in range(52, 64):
            print(f"    r={i}: P={marginals[i]:.4f}")

        # Conditional P(cc[j]≤t | cc[i]≤t) for consecutive pairs
        print(f"\n  Conditional P(cc[j]≤{t} | cc[i]≤{t}):")
        print(f"  {'i→j':>6} | P(j|i)   P(j)    ratio   amplif")
        print("  " + "-"*50)

        for i in range(52, 63):
            j = i + 1
            both = sum(1 for p in profiles if p[i] <= t and p[j] <= t)
            cond_i = sum(1 for p in profiles if p[i] <= t)

            if cond_i > 0:
                p_j_given_i = both / cond_i
                p_j = marginals[j]
                ratio = p_j_given_i / p_j if p_j > 0 else 0
                amplif = "★" if ratio > 1.2 else ("▼" if ratio < 0.8 else "")
                print(f"  {i}→{j} | {p_j_given_i:.4f}  {p_j:.4f}  {ratio:.3f}x  {amplif}")

        # Also check non-consecutive: i, i+2
        print(f"\n  Gap=2: P(cc[j]≤{t} | cc[i]≤{t}):")
        for i in range(52, 62):
            j = i + 2
            both = sum(1 for p in profiles if p[i] <= t and p[j] <= t)
            cond_i = sum(1 for p in profiles if p[i] <= t)
            if cond_i > 0:
                p_j_given_i = both / cond_i
                p_j = marginals[j]
                ratio = p_j_given_i / p_j if p_j > 0 else 0
                amplif = "★" if ratio > 1.2 else ""
                print(f"  {i}→{j} | {p_j_given_i:.4f}  {p_j:.4f}  {ratio:.3f}x  {amplif}")


# ============================================================
# EXP 2: Joint probability of k consecutive zero-carry words
# ============================================================
def exp2_joint_zero(N):
    print("\n" + "="*70)
    print("EXP 2: JOINT PROBABILITY OF k CONSECUTIVE LOW-CARRY WORDS")
    print("  P(cc[52..52+k-1] all ≤ t)")
    print("="*70)

    profiles = []
    for _ in range(N):
        W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        cc, _ = schedule_carry_profile(W)
        profiles.append(cc)

    for t in [0, 3, 5, 8]:
        print(f"\n  --- Threshold ≤ {t} ---")
        print(f"  {'k':>3} | P(joint)    P(indep)    ratio    log2(P)")
        print("  " + "-"*55)

        for k in range(1, 13):
            joint = sum(1 for p in profiles
                       if all(p[52+j] <= t for j in range(k)))
            p_joint = joint / N if N > 0 else 0

            # Independent estimate
            p_indep = 1.0
            for j in range(k):
                marg = sum(1 for p in profiles if p[52+j] <= t) / N
                p_indep *= marg

            ratio = p_joint / p_indep if p_indep > 0 else 0
            import math as _m
            log2p = -999 if p_joint == 0 else (0 if p_joint >= 1 else
                     int(_m.log2(p_joint)))

            marker = " ★★" if ratio > 2 else (" ★" if ratio > 1.3 else "")
            print(f"  {k:3d} | {p_joint:.6f}  {p_indep:.6f}  {ratio:.3f}x  ~2^{log2p}{marker}")


# ============================================================
# EXP 3: Carry correlation STRUCTURE — what creates it?
# ============================================================
def exp3_structure(N):
    print("\n" + "="*70)
    print("EXP 3: WHAT CREATES CARRY CORRELATION?")
    print("  Schedule: W[i] depends on W[i-2], W[i-7], W[i-15], W[i-16]")
    print("  W[i+1] depends on W[i-1], W[i-6], W[i-14], W[i-15]")
    print("  Shared operand: W[i-15] appears in BOTH!")
    print("="*70)

    # How many operands are shared between W[i] and W[i+k]?
    print("\n  Operand overlap between W[i] and W[i+k]:")
    for k in range(1, 8):
        # W[i] uses: i-2, i-7, i-15, i-16
        # W[i+k] uses: i+k-2, i+k-7, i+k-15, i+k-16
        deps_i = {-2, -7, -15, -16}
        deps_ik = {k-2, k-7, k-15, k-16}
        shared = deps_i & deps_ik
        print(f"    k={k}: shared={len(shared)} operands: {shared if shared else '∅'}")

    # Measure: when operands have low HW, is carry lower?
    print("\n  Correlation between operand HW and carry HW:")

    low_op_low_cc = 0
    low_op_total = 0
    high_op_low_cc = 0
    high_op_total = 0

    for _ in range(N):
        W16 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W = list(W16)
        for i in range(16, 64):
            ops = [sig1(W[i-2]), W[i-7], sig0(W[i-15]), W[i-16]]
            total_op_hw = sum(hw(o) for o in ops)
            cc_hw = hw(psi4(*ops))
            W.append(add32(*ops))

            if i >= 52:
                if total_op_hw < 48:  # low operand HW (< 12 avg per op)
                    low_op_total += 1
                    if cc_hw <= 5: low_op_low_cc += 1
                else:
                    high_op_total += 1
                    if cc_hw <= 5: high_op_low_cc += 1

    p_low = low_op_low_cc / low_op_total if low_op_total > 0 else 0
    p_high = high_op_low_cc / high_op_total if high_op_total > 0 else 0
    print(f"    P(cc≤5 | low operand HW): {p_low:.4f} ({low_op_total} samples)")
    print(f"    P(cc≤5 | high operand HW): {p_high:.4f} ({high_op_total} samples)")
    if p_high > 0:
        print(f"    Ratio: {p_low/p_high:.2f}x")


# ============================================================
# EXP 4: W_base optimization for multi-word low carry
# ============================================================
def exp4_optimize_multi_low(N):
    print("\n" + "="*70)
    print("EXP 4: OPTIMIZE W[0..15] FOR MULTI-WORD LOW CARRY IN TAIL")
    print("="*70)

    import math

    def tail_score(W16, t=5):
        cc, _ = schedule_carry_profile(W16)
        return sum(1 for i in range(52, 64) if cc[i] <= t)

    best_count = 0
    best_W = None

    # Random search first
    for trial in range(min(N*10, 10000)):
        W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        count = tail_score(W)
        if count > best_count:
            best_count = count
            best_W = list(W)
            if count >= 4:
                print(f"  [{trial:5d}] {count}/12 words with cc≤5 in tail")

    print(f"\n  Random search: best = {best_count}/12 low-carry words")

    # SA optimization
    if best_W:
        current_W = list(best_W)
        current_count = best_count

        for it in range(min(N*50, 30000)):
            trial_W = list(current_W)
            w = int.from_bytes(os.urandom(1),'big') % 16
            b = int.from_bytes(os.urandom(1),'big') % 32
            trial_W[w] ^= (1 << b)

            tc = tail_score(trial_W)

            T = max(0.01, 1.0 - it/30000)
            if tc > current_count or \
               (tc == current_count and int.from_bytes(os.urandom(1),'big') < 128) or \
               math.exp((tc-current_count)/(T*0.5)) > int.from_bytes(os.urandom(4),'big')/(1<<32):
                current_W = trial_W
                current_count = tc

            if current_count > best_count:
                best_count = current_count
                best_W = list(current_W)
                print(f"  SA [{it:5d}] {best_count}/12 low-carry words")

        print(f"\n  SA optimized: {best_count}/12 low-carry words in tail")

        # Show profile
        cc, _ = schedule_carry_profile(best_W)
        print(f"  Tail carry profile:")
        for i in range(52, 64):
            marker = " ★" if cc[i] <= 5 else ""
            print(f"    r={i}: HW(cc)={cc[i]:2d}{marker}")

        # Joint probability estimate
        n_low = best_count
        if n_low >= 4:
            print(f"\n  ★ {n_low} consecutive-ish low-carry words found!")
            print(f"    If independent: P ≈ 2^{-16*n_low} per message")
            print(f"    Found in {min(N*10,10000)+min(N*50,30000)} tries")
            print(f"    → Effective cost: ~2^{int(math.log2(min(N*10,10000)+min(N*50,30000))+1)}")


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5000

    print("="*70)
    print("SCF ЗАДАЧА B: CARRY CORRELATION IN SCHEDULE")
    print("="*70)

    exp1_pairwise(N)
    exp2_joint_zero(N)
    exp3_structure(N)
    exp4_optimize_multi_low(min(N, 1000))

    print("\n" + "="*70)
    print("ИТОГ ЗАДАЧИ B")
    print("="*70)
