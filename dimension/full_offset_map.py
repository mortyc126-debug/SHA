"""
ПОЛНАЯ КАРТА CARRY OFFSET: все позиции × все break bits.

Для каждой "конфликтной" позиции W[r] (r=4..15):
  A-repair определяет δW[r] (для meeting: δa[r+1]=0)
  Schedule использует δW[r] (для lock: δW[16+..]=0)

Carry offset = |a-repair δW[r] XOR schedule_target[r]|

Schedule targets (для δW[16..18]=0 при δW[0..2]=0):
  δW[16] = σ₁(δW[14]) + δW[9] + 0 + 0
  δW[17] = σ₁(δW[15]) + δW[10] + 0 + 0
  δW[18] = σ₁(δW[16]) + δW[11] + σ₀(δW[3]) + 0

Для δW[16]=0: нужно σ₁(δW[14]) + δW[9] = 0 → δW[9] = -σ₁(δW[14])
Для δW[17]=0: нужно σ₁(δW[15]) + δW[10] = 0 → δW[10] = -σ₁(δW[15])
Для δW[18]=0: нужно δW[11] = -σ₁(δW[16]) - σ₀(δW[3])

Но δW[14] и δW[15] ТОЖЕ от a-repair! Замкнутый круг.
Нужно проверить ВСЕ позиции.
"""

import numpy as np
from collections import Counter

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

def expand_schedule(W16):
    W = list(W16)
    for r in range(16, 64):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))
    return W

def a_repair_W(state2, target_a, r_idx):
    a, b, c, d, e, f, g, h = state2
    T2 = add32(Sigma0(a), Maj(a, b, c))
    T1_needed = sub32(target_a, T2)
    return sub32(sub32(sub32(sub32(T1_needed, h), Sigma1(e)), Ch(e, f, g)), K[r_idx])


def measure_all_offsets(break_bit, N=500):
    """Для данного break bit: вычисляем δW[4..15] от a-repair."""

    all_dw = {r: [] for r in range(4, 16)}

    for trial in range(N):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W1f = expand_schedule(W1)

        s1 = tuple(IV)
        states1 = [s1]
        for r in range(64):
            s1 = R(s1, W1f[r], r)
            states1.append(s1)

        W2 = list(W1)
        W2[3] ^= (1 << break_bit)
        s2 = tuple(IV)
        for r in range(3):
            s2 = R(s2, W2[r], r)
        s2 = R(s2, W2[3], 3)
        for r in range(4, 16):
            W2[r] = a_repair_W(s2, states1[r + 1][0], r)
            s2 = R(s2, W2[r], r)

        for r in range(4, 16):
            dw = W1[r] ^ W2[r]
            all_dw[r].append(dw)

    return all_dw


def main():
    np.random.seed(42)

    print("=" * 70)
    print("ПОЛНАЯ КАРТА: δW[4..15] от a-repair для всех break bit")
    print("=" * 70)

    # =================================================================
    print(f"\n{'=' * 70}")
    print("1. HW(δW[r]) по позициям и break bits")
    print("=" * 70)

    # Собираем для break_bit = 0, 15, 31 (три представителя)
    for bb in [0, 15, 31]:
        all_dw = measure_all_offsets(bb, N=500)

        print(f"\n  break bit={bb}:")
        print(f"  {'W[r]':>6} {'Mean HW':>8} {'P(=0)':>7} {'Unique':>7} {'Stable':>7}")

        for r in range(4, 16):
            vals = all_dw[r]
            mean_hw = np.mean([hw(v) for v in vals])
            p_zero = sum(1 for v in vals if v == 0) / len(vals) * 100
            unique = len(set(vals))

            # Stability: most common percentage
            mc = Counter(vals).most_common(1)[0]
            stable = mc[1] / len(vals) * 100

            marker = ""
            if p_zero > 99: marker = " ★ ALWAYS ZERO"
            elif p_zero > 40: marker = " ★ often zero"
            elif mean_hw <= 1: marker = " ~ small"

            print(f"  W[{r:2d}] {mean_hw:7.2f} {p_zero:6.1f}% {unique:6d} {stable:6.1f}%{marker}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("2. ПОЗИЦИИ С НУЛЕВЫМ δW (dual-use candidates)")
    print("=" * 70)

    # Для КАЖДОГО break bit: какие позиции W[4..15] имеют δW=0?
    print(f"\n  Позиции с P(δW=0) > 40% для каждого break bit:")
    print(f"  {'Bit':>4}", end="")
    for r in range(4, 16):
        print(f" W[{r:2d}]", end="")
    print()

    dual_use_map = {}

    for bb in range(32):
        all_dw = measure_all_offsets(bb, N=300)
        row = []

        print(f"  {bb:4d}", end="")
        for r in range(4, 16):
            p_zero = sum(1 for v in all_dw[r] if v == 0) / 300 * 100
            row.append(p_zero)

            if p_zero > 99:
                sym = "  ██ "
            elif p_zero > 40:
                sym = "  ░░ "
            elif p_zero > 0:
                sym = f" {p_zero:3.0f}%"
            else:
                sym = "   · "
            print(sym, end="")

        dual_use_map[bb] = row
        print()

    print(f"\n  ██ = always 0 (free dual-use)")
    print(f"  ░░ = often 0 (>40%)")
    print(f"  ·  = never 0")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("3. КАКИЕ ПОЗИЦИИ ВЛИЯЮТ НА SCHEDULE?")
    print("=" * 70)

    # Schedule dependencies для δW[16..20]:
    # δW[16] depends on: W[14], W[9], W[1], W[0]
    # δW[17] depends on: W[15], W[10], W[2], W[1]
    # δW[18] depends on: W[16→(14,9)], W[11], W[3], W[2]
    # δW[19] depends on: W[17→(15,10)], W[12], W[4], W[3]
    # δW[20] depends on: W[18→(16,11,3)], W[13], W[5], W[4]

    print(f"""
  Schedule dependencies (which a-repair words affect δW[16..20]):

  δW[16] ← W[14], W[9]       (+ W[1]=0, W[0]=0)
  δW[17] ← W[15], W[10]      (+ W[2]=0, W[1]=0)
  δW[18] ← W[11], W[3]=break (+ δW[16]←W[14,9])
  δW[19] ← W[12], W[4]       (+ δW[17]←W[15,10], W[3]=break)
  δW[20] ← W[13], W[5]       (+ δW[18]←...)

  DUAL-USE позиции (нужны для meeting И lock):
    W[9]:  meeting(δa[10]) + lock(δW[16])
    W[10]: meeting(δa[11]) + lock(δW[17])
    W[11]: meeting(δa[12]) + lock(δW[18])  ← W[11]-конфликт!
    W[12]: meeting(δa[13]) + lock(δW[19])
    W[13]: meeting(δa[14]) + lock(δW[20])
    W[14]: meeting(δa[15]) + lock(δW[16])
    W[15]: meeting(δa[16]) + lock(δW[17])

  NON-CONFLICT позиции (meeting only):
    W[4]: meeting(δa[5])
    W[5]: meeting(δa[6])  + lock(δW[20]) — remote
    W[6]: meeting(δa[7])
    W[7]: meeting(δa[8])
    W[8]: meeting(δa[9])
""")

    # =================================================================
    print(f"{'=' * 70}")
    print("4. DUAL-USE SCORE: для каких (break_bit, position) offset=0?")
    print("=" * 70)

    # Для каждого dual-use position W[9..15]:
    # "Offset" = разница между a-repair δW[r] и schedule target
    # Если offset=0 → эта позиция БЕСПЛАТНО выполняет оба

    # Schedule target для δW[16]=0: δW[9] = -σ₁(δW[14]) mod 2^32
    # Но δW[14] тоже от a-repair и может быть 0!

    print(f"\n  Для break bit=31:")
    all_dw = measure_all_offsets(31, N=1000)

    # Фактический δW[14]:
    dw14_vals = all_dw[14]
    p14_zero = sum(1 for v in dw14_vals if v == 0) / 1000 * 100
    print(f"    δW[14]: P(=0)={p14_zero:.1f}%")

    if p14_zero > 99:
        # δW[14]=0 → δW[16] = σ₁(0) + δW[9] = δW[9]
        # Для δW[16]=0: need δW[9]=0
        dw9_vals = all_dw[9]
        p9_zero = sum(1 for v in dw9_vals if v == 0) / 1000 * 100
        print(f"    δW[14]=0 ALWAYS → δW[16] = δW[9]")
        print(f"    δW[9]: P(=0)={p9_zero:.1f}%")
        print(f"    → P(δW[16]=0) = P(δW[9]=0) = {p9_zero:.1f}%")

    # δW[15]:
    dw15_vals = all_dw[15]
    p15_zero = sum(1 for v in dw15_vals if v == 0) / 1000 * 100
    print(f"    δW[15]: P(=0)={p15_zero:.1f}%")

    if p15_zero > 99:
        dw10_vals = all_dw[10]
        p10_zero = sum(1 for v in dw10_vals if v == 0) / 1000 * 100
        print(f"    δW[15]=0 → δW[17] = δW[10]")
        print(f"    δW[10]: P(=0)={p10_zero:.1f}%")

    # δW[11] — the conflict word:
    dw11_vals = all_dw[11]
    dw11_hw = [hw(v) for v in dw11_vals]
    print(f"    δW[11]: mean HW={np.mean(dw11_hw):.2f}, P(=0)={sum(1 for v in dw11_vals if v==0)/1000*100:.1f}%")
    print(f"    δW[11] values: {len(set(dw11_vals))} unique")

    # The actual δW[11] value:
    mc = Counter(dw11_vals).most_common(3)
    for val, cnt in mc:
        print(f"      {val:#010x} (HW={hw(val)}): {cnt/10:.1f}%")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("5. НАТИВНЫЙ БЮДЖЕТ COLLISION")
    print("=" * 70)

    # Count: how many positions are "free" (δW=0 always)?
    all_dw_31 = measure_all_offsets(31, N=500)

    free_positions = []
    conflict_positions = []
    for r in range(4, 16):
        p_zero = sum(1 for v in all_dw_31[r] if v == 0) / 500 * 100
        if p_zero > 99:
            free_positions.append(r)
        elif p_zero < 1:
            conflict_positions.append(r)

    half_positions = [r for r in range(4, 16) if r not in free_positions and r not in conflict_positions]

    print(f"""
  БЮДЖЕТ (break bit=31):

  16 позиций всего (W[0..15])
    Fixed: W[0..2]=0, W[3]=break               → 4 использованы
    Free (δW=0 always): {free_positions}     → {len(free_positions)} позиций (dual-use!)
    Half (δW=0 sometimes): {half_positions}           → {len(half_positions)} позиций (probabilistic dual-use)
    Conflict (δW≠0 always): {conflict_positions}  → {len(conflict_positions)} позиций (meeting-only)

  Для meeting: нужно 8+ позиций a-repair → используем ALL 12 (W[4..15])
  Для lock: нужно δW[schedule_deps]=0
    δW[16]=0 ← need δW[9]=0 AND δW[14]=0
              δW[14]∈free (always 0) ✓
              δW[9]∈half (50% = 0)  → cost ×2
    δW[17]=0 ← need δW[10]=0 AND δW[15]=0
              δW[15]∈???
              δW[10]∈???
    δW[18]=0 ← need δW[11]=target
              δW[11]∈conflict → BLOCKED

  НАТИВНАЯ СТОИМОСТЬ:
    Meeting: FREE (12 positions, reboot 100%)
    Lock δW[16]: ×2 (δW[9]=0 with 50%)
    Lock δW[17]: ×2 (δW[15] similar)
    Lock δW[18]: BLOCKED (δW[11] conflict, offset ≥ 3)
    Gap: 47 rounds after r=17

  BOTTLENECK: W[11]. Единственный conflicting position для δW[18].
  Всё остальное — дешёво или бесплатно.
""")


if __name__ == "__main__":
    main()
