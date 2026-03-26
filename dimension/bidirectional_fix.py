"""
BIDIRECTIONAL FIX: forward 16 раундов + backward 37 раундов.

Forward fix: δstate[1..16] = 0 через подбор W[1..15] (проверено: 100%)
Backward fix: δstate[32..64] контролируем через backward (37 раундов)

Зазор: раунды 17-31 = 15 неконтролируемых раундов.
(было 48 при forward-only!)

Стратегия:
  1. Forward: W[0..15] → fix δstate[1..16]=0
  2. Backward: от target state[64], используя W[32..63],
     fix δstate[32..64]=0
  3. GAP: раунды 17-31 — schedule определяет δW[16..31]
     → δstate растёт → нужен birthday только на 15 раундов
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

def R_inv(state_next, W_r, r_idx):
    ap, bp, cp, dp, ep, fp, gp, hp = state_next
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
    ap = target_next[0]
    T2 = add32(Sigma0(a), Maj(a, b, c))
    T1 = sub32(ap, T2)
    return sub32(sub32(sub32(sub32(T1, h), Sigma1(e)), Ch(e, f, g)), K[r_idx])

def state_diff(s1, s2):
    return sum(hw(s1[i] ^ s2[i]) for i in range(8))


def bidirectional_experiment():
    np.random.seed(42)

    print("=" * 70)
    print("BIDIRECTIONAL FIX: forward + backward")
    print("=" * 70)

    # =================================================================
    print("\n" + "=" * 70)
    print("1. FORWARD FIX (r=0..16): повторяем проверенный метод")
    print("=" * 70)

    # Trace 1: W1 → forward → states1, H1
    # Trace 2: W2 (δW[0]=1 бит) → forward fix W2[1..15] → δstate[1..16]=0
    # Потом: forward W2 c schedule → state2[17..64], H2

    n_trials = 200
    meeting_points = []  # state на round 32 для обоих traces

    for trial in range(n_trials):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W1f = expand_schedule(W1)

        # Forward trace 1
        s1 = tuple(IV)
        states1 = [s1]
        for r in range(64):
            s1 = R(s1, W1f[r], r)
            states1.append(s1)
        H1 = tuple(add32(IV[i], states1[64][i]) for i in range(8))

        # Forward fix trace 2
        W2 = list(W1)
        W2[0] ^= (1 << (trial % 32))

        s2 = tuple(IV)
        s2 = R(s2, W2[0], 0)  # round 0 with δW
        for r in range(1, 16):
            W2[r] = compute_needed_W(s2, states1[r + 1], r)
            s2 = R(s2, W2[r], r)

        # Verify forward fix
        fwd_fix_ok = (s2 == states1[16])

        # Continue trace 2 with schedule from W2
        W2f = expand_schedule(W2)
        for r in range(16, 64):
            s2 = R(s2, W2f[r], r)

        H2 = tuple(add32(IV[i], s2[i]) for i in range(8))
        dh = sum(hw(H1[i] ^ H2[i]) for i in range(8))

        # Meeting point: δstate at round 32
        # Recompute state at r=32 for both
        s1_32 = states1[32]
        # For trace 2, recompute from states1[16] (=s2 at r=16)
        s2_r = states1[16]  # forward fix ensures this
        for r in range(16, 32):
            s2_r = R(s2_r, W2f[r], r)
        ds_32 = state_diff(s1_32, s2_r)

        meeting_points.append({
            'fwd_ok': fwd_fix_ok,
            'ds_32': ds_32,
            'dh': dh,
        })

    fwd_ok_pct = sum(1 for m in meeting_points if m['fwd_ok']) / n_trials * 100
    mean_ds32 = np.mean([m['ds_32'] for m in meeting_points])
    mean_dh = np.mean([m['dh'] for m in meeting_points])
    min_dh = min(m['dh'] for m in meeting_points)

    print(f"\n  Forward fix OK: {fwd_ok_pct:.0f}%")
    print(f"  δstate[32] (meeting point): {mean_ds32:.1f}")
    print(f"  δH: mean={mean_dh:.1f}, min={min_dh}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("2. BACKWARD TRACE: от state1[64], какой δstate на r=32?")
    print("=" * 70)

    # Backward trace: от state1[64], используя W2f[63..32]
    # Это backward R⁻¹: state1[64] → ... → state[32]
    # С ДРУГИМИ W (W2f вместо W1f): δ появляется

    backward_ds = []
    for trial in range(n_trials):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W1f = expand_schedule(W1)

        W2 = list(W1)
        W2[0] ^= (1 << (trial % 32))
        # Forward fix
        s_tmp = tuple(IV)
        s_tmp = R(s_tmp, W2[0], 0)
        states1_local = [tuple(IV)]
        s1_local = tuple(IV)
        for r in range(64):
            s1_local = R(s1_local, W1f[r], r)
            states1_local.append(s1_local)

        for r in range(1, 16):
            W2[r] = compute_needed_W(s_tmp, states1_local[r + 1], r)
            s_tmp = R(s_tmp, W2[r], r)

        W2f = expand_schedule(W2)

        # State[64] для trace 1
        s64_1 = states1_local[64]

        # Backward от s64_1 с W1f
        sb1 = s64_1
        for r in range(63, 31, -1):
            sb1 = R_inv(sb1, W1f[r], r)

        # Backward от s64_1 с W2f (ДРУГОЙ schedule!)
        sb2 = s64_1
        for r in range(63, 31, -1):
            sb2 = R_inv(sb2, W2f[r], r)

        ds = state_diff(sb1, sb2)
        backward_ds.append(ds)

    print(f"\n  Backward от ОДНОГО state[64], два schedule (W1 vs W2):")
    print(f"  δstate[32]: mean={np.mean(backward_ds):.1f}, min={min(backward_ds)}, max={max(backward_ds)}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("3. ЗАЗОР: δstate по раундам 16..32 (forward) и 32..64 (backward)")
    print("=" * 70)

    # Полный profile для одного trial
    W1 = [np.random.randint(0, 2**32) for _ in range(16)]
    W1f = expand_schedule(W1)

    W2 = list(W1)
    W2[0] ^= 0x80000000

    # Forward fix
    s1_all = [tuple(IV)]
    s = tuple(IV)
    for r in range(64):
        s = R(s, W1f[r], r)
        s1_all.append(s)

    s2_fwd = tuple(IV)
    s2_fwd = R(s2_fwd, W2[0], 0)
    for r in range(1, 16):
        W2[r] = compute_needed_W(s2_fwd, s1_all[r + 1], r)
        s2_fwd = R(s2_fwd, W2[r], r)

    W2f = expand_schedule(W2)

    # Continue forward trace 2
    s2_states = list(s1_all[:17])  # s2[0..16] = s1[0..16] by fix
    s2 = s2_fwd
    for r in range(16, 64):
        s2 = R(s2, W2f[r], r)
        s2_states.append(s2)

    # Backward trace от s1[64] with W2f
    s2_back = [None] * 65
    s2_back[64] = s1_all[64]  # start from SAME state[64]
    for r in range(63, -1, -1):
        s2_back[r] = R_inv(s2_back[r + 1], W2f[r], r)

    print(f"\n  {'Round':>6} {'Fwd δstate':>11} {'Bwd δstate':>11} {'Gap':>6}")
    print(f"  {'-'*6} {'-'*11} {'-'*11} {'-'*6}")

    for r in range(65):
        fwd_ds = state_diff(s1_all[r], s2_states[r]) if r < len(s2_states) else 0
        bwd_ds = state_diff(s1_all[r], s2_back[r]) if s2_back[r] else 0

        marker = ""
        if fwd_ds < 2 and bwd_ds < 2:
            marker = " ★ BOTH FIXED"
        elif fwd_ds < 2:
            marker = " ← fwd fixed"
        elif bwd_ds < 2:
            marker = " ← bwd fixed"

        gap = min(fwd_ds, bwd_ds)

        if r <= 20 or r >= 28 or r in [24, 25, 26]:
            print(f"  r={r:3d}  {fwd_ds:10d}  {bwd_ds:10d}  {gap:5d}{marker}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("4. MEETING POINT: где forward и backward сходятся?")
    print("=" * 70)

    # Ищем раунд где fwd_ds ≈ bwd_ds (минимальный суммарный gap)
    fwd_profile = []
    bwd_profile = []
    for r in range(65):
        fwd_ds = state_diff(s1_all[r], s2_states[r]) if r < len(s2_states) else 256
        bwd_ds = state_diff(s1_all[r], s2_back[r]) if s2_back[r] else 256
        fwd_profile.append(fwd_ds)
        bwd_profile.append(bwd_ds)

    # Meeting = where min(fwd, bwd) is minimized
    best_meeting = 65
    best_gap = 999
    for r in range(65):
        gap = max(fwd_profile[r], bwd_profile[r])  # both must be small
        combined = fwd_profile[r] + bwd_profile[r]
        if combined < best_gap:
            best_gap = combined
            best_meeting = r

    print(f"\n  Best meeting point: round {best_meeting}")
    print(f"  Forward δ at meeting: {fwd_profile[best_meeting]}")
    print(f"  Backward δ at meeting: {bwd_profile[best_meeting]}")
    print(f"  Combined: {best_gap}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("5. БУХГАЛТЕРИЯ BIDIRECTIONAL")
    print("=" * 70)

    # Forward: fix r=0..16 (δstate=0)
    # Backward: fix r=32..64 (δstate идёт только по e-chain, медленно)
    # Gap: r=17..31 (15 раундов)

    # В зазоре δstate растёт от 0 до ~128 за 15 раундов (forward side)
    # И от ~small до 128 (backward side)

    # Meeting at r=best_meeting: δstate = best_gap бит
    # Это 256 бит условий (state match at meeting)
    # Birthday на 256 бит = 2^128... но!

    # У нас 512 бит свободы (W[0..15])
    # 256 бит потрачено на forward fix (W[1..15])
    # Оставшиеся ~32 бит (W[0]) + backward вариации
    # Backward использует W[32..63] = schedule от W[0..15] — не free!

    fwd_fixed = sum(1 for r in range(17) if fwd_profile[r] < 2)
    bwd_fixed = sum(1 for r in range(33, 65) if bwd_profile[r] < 50)
    gap_rounds = 64 - fwd_fixed - bwd_fixed

    print(f"""
  FORWARD:  {fwd_fixed} раундов зафиксированы (δstate < 2)
  BACKWARD: {bwd_fixed} раундов контролируемы (δstate < 50)
  GAP:      {gap_rounds} неконтролируемых раундов

  При forward-only: gap = 48 раундов → 2^128
  При bidirectional: gap = {gap_rounds} раундов

  КЛЮЧЕВАЯ ПРОБЛЕМА: backward использует schedule(W[0..15])
  Forward fix МЕНЯЕТ W[1..15] → schedule ДРУГОЙ →
  backward trace ДРУГОЙ → backward fix не совместим с forward fix.

  Forward и backward СВЯЗАНЫ через одни и те же W[0..15].
  Нельзя оптимизировать их НЕЗАВИСИМО.

  НО: backward диффузия МЕДЛЕННЕЕ →
  даже при зависимости, backward сторона БОЛЕЕ КОНТРОЛИРУЕМА.
  Meeting point ближе к середине → gap МЕНЬШЕ → birthday ДЕШЕВЛЕ.
""")


if __name__ == "__main__":
    bidirectional_experiment()
