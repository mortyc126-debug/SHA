"""
CARRY OFFSET 4 БИТА: определяем и компенсируем.

Факт: a-repair δW[11] отличается от GF2 target РОВНО на 4 бита.
Всегда. 100K из 100K. Детерминистически.

Вопросы:
  1. КАКИЕ 4 бита? Всегда одни и те же позиции?
  2. Зависят ли позиции от W₁?
  3. Можно ли flip'нуть эти 4 бита в W₂[11] и сохранить a-repair?
  4. Если flip — что происходит с reboot?
"""

import numpy as np

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


def main():
    np.random.seed(42)

    target_gf2 = sigma0(1 << 31)

    print("=" * 70)
    print("CARRY OFFSET 4 БИТА: идентификация и компенсация")
    print("=" * 70)

    # =================================================================
    print(f"\n{'=' * 70}")
    print("1. КАКИЕ 4 бита? Позиции carry offset")
    print("=" * 70)

    N = 10000
    xor_masks = []  # actual_dw11 XOR target

    for trial in range(N):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W1f = expand_schedule(W1)
        s1 = tuple(IV)
        states1 = [s1]
        for r in range(64):
            s1 = R(s1, W1f[r], r)
            states1.append(s1)

        W2 = list(W1)
        W2[3] ^= (1 << 31)
        s2 = tuple(IV)
        for r in range(3):
            s2 = R(s2, W2[r], r)
        s2 = R(s2, W2[3], 3)
        for r in range(4, 12):
            W2[r] = a_repair_W(s2, states1[r + 1][0], r)
            s2 = R(s2, W2[r], r)

        actual_dw11 = W1[11] ^ W2[11]
        xor_mask = actual_dw11 ^ target_gf2
        xor_masks.append(xor_mask)

    # Анализируем XOR mask
    unique_masks = set(xor_masks)
    print(f"\n  Target GF2: {target_gf2:#010x}")
    print(f"  Unique XOR masks: {len(unique_masks)}")

    if len(unique_masks) <= 20:
        print(f"\n  ВСЕ маски:")
        from collections import Counter
        cnt = Counter(xor_masks)
        for mask, count in cnt.most_common(20):
            pct = count / N * 100
            bits = [i for i in range(32) if (mask >> i) & 1]
            print(f"    {mask:#010x} (HW={hw(mask)}, bits={bits}): {count:5d} ({pct:.1f}%)")
    else:
        # Частота каждой битовой позиции
        bit_freq = np.zeros(32)
        for mask in xor_masks:
            for b in range(32):
                if (mask >> b) & 1:
                    bit_freq[b] += 1
        bit_freq /= N

        print(f"\n  Частота бит в XOR mask:")
        for b in range(32):
            bar = "█" * int(bit_freq[b] * 50)
            marker = " ★" if bit_freq[b] > 0.9 else ""
            if bit_freq[b] > 0.01 or b < 5 or b > 27:
                print(f"    bit {b:2d}: {bit_freq[b]:.4f}  {bar}{marker}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("2. КОМПЕНСАЦИЯ: flip carry offset в W₂[11]")
    print("=" * 70)

    # Берём a-repair W₂[11], flip'аем carry offset → новый W₂[11]
    # Пересчитываем: reboot работает? schedule lock?

    compensated_reboot = 0
    compensated_dw18 = []
    compensated_profiles = []

    # Используем ПЕРВУЮ (или самую частую) маску
    if len(unique_masks) == 1:
        compensation_mask = list(unique_masks)[0]
    else:
        compensation_mask = Counter(xor_masks).most_common(1)[0][0]

    print(f"\n  Compensation mask: {compensation_mask:#010x} (HW={hw(compensation_mask)})")

    for trial in range(5000):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W1f = expand_schedule(W1)
        s1 = tuple(IV)
        states1 = [s1]
        for r in range(64):
            s1 = R(s1, W1f[r], r)
            states1.append(s1)

        # A-repair standard
        W2 = list(W1)
        W2[3] ^= (1 << 31)
        s2 = tuple(IV)
        for r in range(3):
            s2 = R(s2, W2[r], r)
        s2 = R(s2, W2[3], 3)
        for r in range(4, 11):
            W2[r] = a_repair_W(s2, states1[r + 1][0], r)
            s2 = R(s2, W2[r], r)

        # Round 11: a-repair PLUS compensation
        W2_11_arepair = a_repair_W(s2, states1[12][0], 11)
        W2_11_compensated = W2_11_arepair ^ compensation_mask
        W2[11] = W2_11_compensated
        s2 = R(s2, W2[11], 11)

        # Continue a-repair for r=12..15
        for r in range(12, 16):
            W2[r] = a_repair_W(s2, states1[r + 1][0], r)
            s2 = R(s2, W2[r], r)

        W2f = expand_schedule(W2)

        # Check reboot
        s2_check = tuple(IV)
        profile = []
        for r in range(64):
            w = W2[r] if r < 16 else W2f[r]
            s2_check = R(s2_check, w, r)
            profile.append(state_diff(states1[r + 1], s2_check))
        compensated_profiles.append(profile)

        if profile[11] == 0:  # reboot at r=12
            compensated_reboot += 1

        compensated_dw18.append(hw(W1f[18] ^ W2f[18]))

    mean_profile = np.mean(compensated_profiles, axis=0)
    both_pct = [np.mean([p[r] == 0 for p in compensated_profiles]) * 100 for r in range(64)]

    print(f"\n  С компенсацией (5000 trials):")
    print(f"  {'r':>4} {'δstate':>8} {'BOTH%':>7}")
    for r in range(22):
        marker = ""
        if both_pct[r] > 99: marker = " ★ZERO"
        elif both_pct[r] > 10: marker = " ~pass"
        print(f"  {r:4d} {mean_profile[r]:7.1f} {both_pct[r]:6.1f}%{marker}")

    print(f"\n  Reboot r=12: {compensated_reboot}/5000 = {compensated_reboot/50:.1f}%")
    print(f"  δW[18]: mean={np.mean(compensated_dw18):.1f}, P(=0)={sum(1 for x in compensated_dw18 if x==0)/50:.1f}%")

    # Без компенсации (baseline)
    print(f"\n  Сравнение:")
    print(f"    {'':>25} {'Reboot%':>8} {'δW[18]':>8} {'r=18 BOTH%':>11}")
    print(f"    {'A-repair (standard)':>25} {'100.0%':>8} {'7.7':>8} {'24.8%':>11}")
    print(f"    {'A-repair + compensation':>25} {compensated_reboot/50:>7.1f}% {np.mean(compensated_dw18):>7.1f} {both_pct[17]:>10.1f}%")

    # =================================================================
    if compensated_reboot == 0:
        print(f"\n{'=' * 70}")
        print("3. КОМПЕНСАЦИЯ УБИВАЕТ REBOOT. Анализ почему.")
        print("=" * 70)

        # Compensation flip'ает 4 бита W₂[11] → δa[12] ≠ 0
        # Потому что W₂[11] определяет a₂[12] через round function
        # Flip 4 бита W → flip ~4 бита T1 → flip ~4 бита a'
        print(f"""
  W₂[11] определяет a₂[12]:
    a₂[12] = T1(state₂[11], W₂[11]) + T2(state₂[11])

  A-repair: W₂[11] выбран чтобы a₂[12] = a₁[12] (δa=0)
  Compensation: flip 4 бита W₂[11] → T1 меняется → a₂[12] ≠ a₁[12]

  Flip в W входит в T1 = h + Σ₁(e) + Ch + K + W.
  δT1 = δW = compensation_mask (4 бита)
  δa' = δT1 + δT2 = compensation_mask + 0 = compensation_mask
  δa[12] = {hw(compensation_mask)} бит ← REBOOT РАЗРУШЕН

  Это ФУНДАМЕНТАЛЬНО: W[11] определяет И a[12] (round) И δW[18] (schedule).
  Нельзя flip W[11] для schedule без поломки round.
  Это W[11]-КОНФЛИКТ: round vs schedule привязаны к ОДНОМУ слову.
        """)

        # Но: δa[12] = 4 бита. Это МАЛО!
        # Может a-repair на r=12..15 может ПОЧИНИТЬ эти 4 бита?
        print(f"  НО: δa[12] = {hw(compensation_mask)} бита — МАЛО!")
        print(f"  A-repair на r=12 ПОЧИНИТ δa[13]=0?")
        print(f"  r=12 BOTH%={both_pct[11]:.1f}%, r=13 BOTH%={both_pct[12]:.1f}%")

        # Проверяем: что если a-repair продолжает ПОСЛЕ compensation?
        # r=12..15 a-repair уже включен → δa[13..16]=0
        # Reboot на r=12 сломан, но r=13+ a-repair ВОССТАНАВЛИВАЕТ!

        # Нужен ли reboot вообще? Или a-repair НАПРЯМУЮ фиксирует r=12..16?
        print(f"\n  Profile после compensation с a-repair r=12..15:")
        for r in range(10, 20):
            print(f"    r={r+1}: δstate={mean_profile[r]:.1f}, BOTH={both_pct[r]:.1f}%")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("4. НОВАЯ СТРАТЕГИЯ: compensation + extended a-repair")
    print("=" * 70)

    # Reboot сломан на r=12 (δa=4 бит).
    # Но a-repair на r=12 ВОССТАНАВЛИВАЕТ δa[13]=0.
    # Потом reboot ПОЗЖЕ: r=12+8=20? Или convergence быстрее?

    # Посмотрим convergence после compensation
    print(f"\n  Convergence после compensation (a-repair r=4..10, comp r=11, a-repair r=12..15):")

    # Reboot-like: когда δstate впервые = 0 ПОСЛЕ r=11?
    reboot_after = []
    for p in compensated_profiles:
        found = False
        for r in range(12, 25):
            if p[r] == 0:
                reboot_after.append(r + 1)
                found = True
                break
        if not found:
            reboot_after.append(-1)

    valid = [r for r in reboot_after if r > 0]
    print(f"    Re-convergence: {len(valid)}/5000 = {len(valid)/50:.1f}%")
    if valid:
        from collections import Counter
        cnt = Counter(valid)
        for r, c in sorted(cnt.items()):
            print(f"      r={r}: {c} ({c/50:.1f}%)")


if __name__ == "__main__":
    main()
