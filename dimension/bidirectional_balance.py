"""
БАЛАНС FORWARD/BACKWARD: ищем оптимальную точку.

Проблема: forward fix (16 раундов) портит schedule → backward слабеет.
Идея: что если фиксить МЕНЬШЕ forward (оставить W стабильнее) →
      backward получит БОЛЬШЕ контроля?

Варианты:
  A. Forward 16 + Backward 4 = 20 fixed (текущий)
  B. Forward 8 + Backward ? = ? fixed (меньше forward портежа)
  C. Forward 4 + Backward ? = ?
  D. Forward 0 + Backward ? = ?
  E. Backward-only (без forward fix)
  F. Partial forward: fix только a-chain или только e-chain

Также: что если forward fix ВЫБОРОЧНО — не все W, а только некоторые?
"""

import numpy as np
import struct, hashlib

MASK32 = 0xFFFFFFFF
K = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
    0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
    0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
    0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
    0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
    0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
]
IV = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]

def hw(x): return bin(x).count('1')
def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
def add32(x, y): return (x + y) & MASK32
def sub32(x, y): return (x - y) & MASK32
def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)

def R(state, W_r, r_idx):
    a, b, c, d, e, f, g, h = state
    T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e, f, g)), K[r_idx]), W_r)
    T2 = add32(Sigma0(a), Maj(a, b, c))
    return (add32(T1, T2), a, b, c, add32(d, T1), e, f, g)

def R_inv(sn, W_r, r_idx):
    ap, bp, cp, dp, ep, fp, gp, hp = sn
    a, b, c, e, f, g = bp, cp, dp, fp, gp, hp
    T2 = add32(Sigma0(a), Maj(a, b, c))
    T1 = sub32(ap, T2)
    d = sub32(ep, T1)
    h = sub32(sub32(sub32(sub32(T1, Sigma1(e)), Ch(e, f, g)), K[r_idx]), W_r)
    return (a, b, c, d, e, f, g, h)

def expand_schedule(W16):
    W = list(W16)
    for r in range(16, 64):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))
    return W

def compute_needed_W(state, target_next, r_idx):
    a, b, c, d, e, f, g, h = state
    T2 = add32(Sigma0(a), Maj(a, b, c))
    T1 = sub32(target_next[0], T2)
    return sub32(sub32(sub32(sub32(T1, h), Sigma1(e)), Ch(e, f, g)), K[r_idx])

def state_diff(s1, s2):
    return sum(hw(s1[i] ^ s2[i]) for i in range(8))


def measure_bidirectional(fwd_fix_depth, n_trials=300):
    """
    Измеряем: при forward fix глубиной fwd_fix_depth,
    какой backward контроль и общий gap?
    """
    results = []

    for trial in range(n_trials):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W1f = expand_schedule(W1)

        # Forward trace 1: полный
        s1 = tuple(IV)
        states1 = [s1]
        for r in range(64):
            s1 = R(s1, W1f[r], r)
            states1.append(s1)

        # Forward fix trace 2: меняем W[0], фиксим W[1..fwd_fix_depth-1]
        W2 = list(W1)
        W2[0] ^= (1 << (trial % 32))

        s2 = tuple(IV)
        s2 = R(s2, W2[0], 0)

        actual_fix = min(fwd_fix_depth, 15)
        for r in range(1, actual_fix + 1):
            W2[r] = compute_needed_W(s2, states1[r + 1], r)
            s2 = R(s2, W2[r], r)

        # Оставшиеся W[fix+1..15] = W1 (не трогаем!)
        W2f = expand_schedule(W2)

        # Continue forward
        for r in range(actual_fix + 1, 64):
            s2 = R(s2, W2f[r], r)

        # Forward δstate profile
        s2_fwd = tuple(IV)
        s2_fwd = R(s2_fwd, W2[0], 0)
        fwd_states = [tuple(IV), s2_fwd]
        for r in range(1, actual_fix + 1):
            s2_fwd = R(s2_fwd, W2[r], r)
            fwd_states.append(s2_fwd)
        for r in range(actual_fix + 1, 64):
            s2_fwd = R(s2_fwd, W2f[r], r)
            fwd_states.append(s2_fwd)

        fwd_profile = [state_diff(states1[r], fwd_states[r]) if r < len(fwd_states) else 256
                       for r in range(65)]

        # Backward from state1[64] with W2f schedule
        bwd_states = [None] * 65
        bwd_states[64] = states1[64]
        for r in range(63, -1, -1):
            bwd_states[r] = R_inv(bwd_states[r + 1], W2f[r], r)

        bwd_profile = [state_diff(states1[r], bwd_states[r]) for r in range(65)]

        # Measure: forward fixed rounds, backward controlled, gap
        fwd_fixed = sum(1 for r in range(actual_fix + 2) if fwd_profile[r] < 2)
        bwd_controlled = sum(1 for r in range(55, 65) if bwd_profile[r] < 50)

        # Schedule damage
        sched_damage = sum(hw(W1f[r] ^ W2f[r]) for r in range(16, 64))

        # HW(δH)
        H1 = tuple(add32(IV[i], states1[64][i]) for i in range(8))
        H2 = tuple(add32(IV[i], s2[i]) for i in range(8))
        dh = sum(hw(H1[i] ^ H2[i]) for i in range(8))

        results.append({
            'fwd_fixed': fwd_fixed,
            'bwd_controlled': bwd_controlled,
            'sched_damage': sched_damage,
            'dh': dh,
            'bwd_r63': bwd_profile[63],
            'bwd_r62': bwd_profile[62],
            'bwd_r61': bwd_profile[61],
            'bwd_r60': bwd_profile[60],
        })

    return results


def main():
    np.random.seed(42)

    print("=" * 70)
    print("БАЛАНС FORWARD/BACKWARD: поиск оптимума")
    print("=" * 70)

    print(f"\n  {'FwdDepth':>8} {'FwdFixed':>9} {'BwdCtrl':>8} {'SchedDmg':>9} "
          f"{'δr63':>5} {'δr62':>5} {'δr61':>5} {'δr60':>5} {'HW(δH)':>8} {'minδH':>6}")
    print(f"  {'-'*8} {'-'*9} {'-'*8} {'-'*9} {'-'*5} {'-'*5} {'-'*5} {'-'*5} {'-'*8} {'-'*6}")

    best_config = None
    best_bwd = 0

    for fwd_depth in [0, 1, 2, 3, 4, 6, 8, 10, 12, 15]:
        results = measure_bidirectional(fwd_depth, n_trials=300)

        mean_fwd = np.mean([r['fwd_fixed'] for r in results])
        mean_bwd = np.mean([r['bwd_controlled'] for r in results])
        mean_sd = np.mean([r['sched_damage'] for r in results])
        mean_dh = np.mean([r['dh'] for r in results])
        min_dh = min(r['dh'] for r in results)
        mean_r63 = np.mean([r['bwd_r63'] for r in results])
        mean_r62 = np.mean([r['bwd_r62'] for r in results])
        mean_r61 = np.mean([r['bwd_r61'] for r in results])
        mean_r60 = np.mean([r['bwd_r60'] for r in results])

        marker = ""
        if mean_bwd > best_bwd:
            best_bwd = mean_bwd
            best_config = fwd_depth
            marker = " ★"

        print(f"  {fwd_depth:>8} {mean_fwd:>8.1f} {mean_bwd:>7.1f} {mean_sd:>8.0f} "
              f"{mean_r63:>4.0f} {mean_r62:>4.0f} {mean_r61:>4.0f} {mean_r60:>4.0f} "
              f"{mean_dh:>7.1f} {min_dh:>5}{marker}")

    print(f"\n  ЛУЧШИЙ БАЛАНС: forward depth = {best_config}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print(f"ДЕТАЛЬНО: forward depth = {best_config}")
    print("=" * 70)

    results = measure_bidirectional(best_config, n_trials=500)

    # Distribution of backward control
    bwd_ctrls = [r['bwd_controlled'] for r in results]
    print(f"\n  Backward controlled rounds distribution:")
    for v in sorted(set(bwd_ctrls)):
        count = bwd_ctrls.count(v)
        bar = "█" * (count * 50 // len(results))
        print(f"    {v} rounds: {count:4d} ({count/len(results)*100:5.1f}%)  {bar}")

    # δH distribution
    dhs = [r['dh'] for r in results]
    print(f"\n  HW(δH) distribution:")
    print(f"    Mean:  {np.mean(dhs):.1f}")
    print(f"    Std:   {np.std(dhs):.1f}")
    print(f"    Min:   {min(dhs)}")
    print(f"    <100:  {sum(1 for d in dhs if d < 100)} / {len(dhs)}")
    print(f"    <110:  {sum(1 for d in dhs if d < 110)} / {len(dhs)}")

    # Schedule damage vs backward control
    sds = [r['sched_damage'] for r in results]
    bwds = [r['bwd_controlled'] for r in results]
    corr_sd_bwd = np.corrcoef(sds, bwds)[0, 1]
    print(f"\n  Corr(schedule_damage, bwd_controlled): {corr_sd_bwd:+.4f}")
    print(f"  → {'Schedule damage ОПРЕДЕЛЯЕТ backward!' if abs(corr_sd_bwd) > 0.2 else 'Слабая связь'}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("ИТОГ: ОПТИМАЛЬНЫЙ BIDIRECTIONAL")
    print("=" * 70)

    # Подсчёт: сколько мы выиграли
    total_fixed_best = np.mean([r['fwd_fixed'] + r['bwd_controlled'] for r in results])
    gap_best = 64 - total_fixed_best

    print(f"""
  Forward depth: {best_config}
  Forward fixed: {np.mean([r['fwd_fixed'] for r in results]):.1f} раундов
  Backward controlled: {np.mean([r['bwd_controlled'] for r in results]):.1f} раундов
  Total controlled: {total_fixed_best:.1f} раундов
  Gap: {gap_best:.1f} раундов

  Сравнение:
    Forward-only (depth=15): gap ≈ 48
    Bidirectional optimal:   gap ≈ {gap_best:.0f}
    Выигрыш: {48 - gap_best:.0f} раундов

  В единицах нашего измерения:
    Каждый контролируемый раунд backward =
    одна "тихая" зона в ткани (e-chain не диффундирует).
    {total_fixed_best:.0f}/64 раундов под контролем = {total_fixed_best/64*100:.0f}% ткани.
""")


if __name__ == "__main__":
    main()
