"""
КАРТА СТОИМОСТЕЙ АТАКИ в нашем измерении.

Для каждого break bit (0..31):
  1. Meeting cost (reboot round, probability)
  2. Schedule transparency (δW[16], δW[17], δW[18])
  3. Carry offset (HW, variants)
  4. Gap size (rounds after last BOTH=0)
  5. TOTAL COST estimate

Единицы: NOT bits. Нативные единицы нашего измерения:
  - Meeting cost = число попыток для reboot
  - Lock cost = число попыток для schedule transparency
  - Gap cost = birthday на δstate at meeting point
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

def state_diff(s1, s2):
    return sum(hw(s1[i] ^ s2[i]) for i in range(8))


def full_analysis(break_bit, N=2000):
    """Полный анализ для одного break bit."""
    results = {
        'reboot_round': [],
        'reboot_pct': 0,
        'dw16': [], 'dw17': [], 'dw18': [],
        'r17_pass': 0, 'r18_pass': 0,
        'last_zero': [],
        'dh': [],
        'offset_hw': 0,
        'offset_masks': [],
    }

    target_gf2 = sigma0(1 << break_bit)

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

        W2f = expand_schedule(W2)

        # Profile
        s2_full = tuple(IV)
        profile = []
        for r in range(64):
            w = W2[r] if r < 16 else W2f[r]
            s2_full = R(s2_full, w, r)
            profile.append(state_diff(states1[r + 1], s2_full))

        # Reboot
        for r in range(8, 20):
            if profile[r] == 0 and r > 6:
                results['reboot_round'].append(r + 1)
                break

        # Schedule
        results['dw16'].append(hw(W1f[16] ^ W2f[16]))
        results['dw17'].append(hw(W1f[17] ^ W2f[17]))
        results['dw18'].append(hw(W1f[18] ^ W2f[18]))

        # Pass rates
        if len(profile) > 17 and profile[16] == 0 and profile[17] == 0:
            results['r17_pass'] += 1
        if len(profile) > 18 and profile[16] == 0 and profile[17] == 0 and profile[18] == 0:
            results['r18_pass'] += 1

        # Last BOTH=0
        last_z = 0
        for r in range(64):
            if profile[r] == 0:
                last_z = r + 1
        results['last_zero'].append(last_z)

        # Hash diff
        H1 = tuple(add32(IV[i], states1[64][i]) for i in range(8))
        H2 = tuple(add32(IV[i], s2_full[i]) for i in range(8))
        dh = sum(hw(H1[i] ^ H2[i]) for i in range(8))
        results['dh'].append(dh)

        # Offset
        actual_dw11 = W1[11] ^ W2[11]
        offset = actual_dw11 ^ target_gf2
        results['offset_masks'].append(offset)

    # Aggregate
    results['reboot_pct'] = len(results['reboot_round']) / N * 100
    results['r17_pass'] = results['r17_pass'] / N * 100
    results['r18_pass'] = results['r18_pass'] / N * 100

    offset_counter = Counter(results['offset_masks'])
    main_offset = offset_counter.most_common(1)[0][0]
    results['offset_hw'] = hw(main_offset)
    results['offset_unique'] = len(offset_counter)

    return results


def main():
    np.random.seed(42)

    print("=" * 70)
    print("КАРТА СТОИМОСТЕЙ АТАКИ В НАШЕМ ИЗМЕРЕНИИ")
    print("=" * 70)

    # =================================================================
    print(f"\n{'=' * 70}")
    print("1. ПОЛНАЯ КАРТА: все 32 break бита")
    print("=" * 70)

    all_results = {}
    for bit in range(32):
        r = full_analysis(bit, N=1000)
        all_results[bit] = r

    print(f"\n  {'Bit':>4} {'Reboot':>7} {'δW16':>5} {'δW17':>5} {'δW18':>5} {'r17%':>6} {'r18%':>6} "
          f"{'LastZ':>6} {'OffHW':>6} {'minδH':>6}")
    print(f"  {'-'*4} {'-'*7} {'-'*5} {'-'*5} {'-'*5} {'-'*6} {'-'*6} {'-'*6} {'-'*6} {'-'*6}")

    for bit in range(32):
        r = all_results[bit]
        reboot = f"r={int(np.mean(r['reboot_round']))}" if r['reboot_round'] else "NONE"
        last_z = np.mean(r['last_zero'])
        print(f"  {bit:4d} {reboot:>7} {np.mean(r['dw16']):4.1f} {np.mean(r['dw17']):4.1f} "
              f"{np.mean(r['dw18']):4.1f} {r['r17_pass']:5.1f} {r['r18_pass']:5.1f} "
              f"{last_z:5.1f} {r['offset_hw']:5d} {min(r['dh']):5d}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("2. СТОИМОСТЬ В НАШЕМ ИЗМЕРЕНИИ")
    print("=" * 70)

    print(f"\n  Единицы нашего измерения:")
    print(f"    Meeting cost = 1/P(reboot)")
    print(f"    Lock cost = 1/P(δW[16..last]=0)")
    print(f"    Gap cost = birthday на δstate at meeting boundary")
    print(f"    TOTAL = Meeting × Lock × Gap")

    print(f"\n  {'Bit':>4} {'MeetCost':>9} {'LockCost':>9} {'GapRounds':>10} {'GapCost':>9} {'TOTAL':>12}")
    print(f"  {'-'*4} {'-'*9} {'-'*9} {'-'*10} {'-'*9} {'-'*12}")

    best_bit = -1
    best_total = 999

    for bit in range(32):
        r = all_results[bit]

        # Meeting cost
        if r['reboot_pct'] > 0:
            meet_cost_exp = np.log2(100 / r['reboot_pct'])
        else:
            meet_cost_exp = 20  # never found

        # Lock cost: P(schedule transparent to r=last_zero)
        # We use r17_pass as proxy for lock
        if r['r17_pass'] > 0:
            lock_cost_exp = np.log2(100 / r['r17_pass'])
        else:
            lock_cost_exp = 14  # not found in sample

        # Gap: from last_zero to r=64
        mean_last = np.mean(r['last_zero'])
        gap_rounds = max(64 - mean_last, 1)

        # Gap cost in our dimension:
        # Each gap round: δstate grows by ~32 bits
        # After 4 rounds: saturated at 128
        # Birthday on saturated state: 2^128
        # BUT: if gap < 4 rounds, not saturated!
        if gap_rounds <= 4:
            gap_state_bits = gap_rounds * 32
            gap_cost_exp = gap_state_bits / 2
        else:
            gap_cost_exp = 128  # saturated → standard birthday

        total_exp = meet_cost_exp + lock_cost_exp + gap_cost_exp

        if total_exp < best_total:
            best_total = total_exp
            best_bit = bit

        marker = ""
        if total_exp < 135: marker = " ★"
        if total_exp < 130: marker = " ★★"
        if total_exp < 120: marker = " ★★★"

        if bit < 5 or bit > 28 or marker or r['r17_pass'] > 10:
            print(f"  {bit:4d} 2^{meet_cost_exp:5.1f}  2^{lock_cost_exp:5.1f}  {gap_rounds:9.1f}  "
                  f"2^{gap_cost_exp:5.0f}   2^{total_exp:6.1f}{marker}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print(f"3. ЛУЧШАЯ АТАКА: bit={best_bit}")
    print("=" * 70)

    r = all_results[best_bit]
    print(f"""
  ОПТИМАЛЬНЫЙ BREAK BIT: {best_bit}

  Meeting:
    Reboot: {r['reboot_pct']:.1f}%
    Cost: 2^{np.log2(max(100/max(r['reboot_pct'],0.1), 1)):.1f}

  Schedule lock:
    δW[16] = {np.mean(r['dw16']):.1f} бит (P(=0) = {sum(1 for x in r['dw16'] if x==0)/len(r['dw16'])*100:.1f}%)
    δW[17] = {np.mean(r['dw17']):.1f} бит (P(=0) = {sum(1 for x in r['dw17'] if x==0)/len(r['dw17'])*100:.1f}%)
    δW[18] = {np.mean(r['dw18']):.1f} бит (P(=0) = {sum(1 for x in r['dw18'] if x==0)/len(r['dw18'])*100:.1f}%)
    r=17 pass: {r['r17_pass']:.1f}%
    r=18 pass: {r['r18_pass']:.1f}%

  Carry offset: HW={r['offset_hw']}, {r['offset_unique']} variants

  Gap: {64 - np.mean(r['last_zero']):.1f} rounds

  Hash: HW(δH) mean={np.mean(r['dh']):.1f}, min={min(r['dh'])}

  TOTAL COST: 2^{best_total:.1f}
  Standard birthday: 2^128
  GAIN: {128 - best_total:.1f} бит
""")

    # =================================================================
    print(f"{'=' * 70}")
    print("4. СВОДКА: нативная vs стандартная стоимость")
    print("=" * 70)

    print(f"""
  ┌─────────────────────────────────────────────────────┐
  │ КАРТА АТАК В НАШЕМ ИЗМЕРЕНИИ                        │
  ├─────────────────────────────────────────────────────┤
  │ Standard birthday:              2^128               │
  │                                                     │
  │ A-repair (bit=31):                                  │
  │   Meeting r=12 (free) + r17 pass 51%               │
  │   Gap 46 rounds → 2^128                             │
  │   Total: 2^128 (no gain on full collision)          │
  │                                                     │
  │ A-repair (bit={best_bit}, optimal):                     │
  │   Meeting: 2^{np.log2(max(100/max(all_results[best_bit]['reboot_pct'],0.1),1)):.1f}                                  │
  │   Lock: 2^{np.log2(max(100/max(all_results[best_bit]['r17_pass'],0.1),1)):.1f}                                    │
  │   Gap: 2^{128 if 64-np.mean(all_results[best_bit]['last_zero']) > 4 else int((64-np.mean(all_results[best_bit]['last_zero']))*16)}                                     │
  │   Total: 2^{best_total:.0f}                                  │
  │                                                     │
  │ GAIN: {128-best_total:.0f} бит                                     │
  │                                                     │
  │ Carry offset: ≥3 бит (архитектурный минимум)        │
  │ W[11]-конфликт: round vs schedule НЕСОВМЕСТИМЫ      │
  │ Schedule lock: единственный оставшийся барьер       │
  └─────────────────────────────────────────────────────┘
""")


if __name__ == "__main__":
    main()
