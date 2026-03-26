"""
FABRIC REPAIR: ремонт ткани — оптимизация под SCHEDULE LOCK.

A-repair: оптимизирует δa=0 → schedule = побочный эффект
Fabric repair: оптимизирует δW[16..]=0 → state = подстраиваем

Schedule над GF(2) (без carry):
  δW[16] = σ₁(δW[14]) ⊕ δW[9] ⊕ σ₀(δW[1]) ⊕ δW[0]
  δW[17] = σ₁(δW[15]) ⊕ δW[10] ⊕ σ₀(δW[2]) ⊕ δW[1]
  δW[18] = σ₁(δW[16]) ⊕ δW[11] ⊕ σ₀(δW[3]) ⊕ δW[2]

При δW[0..2]=0, δW[3]=break:
  δW[16] = σ₁(δW[14]) ⊕ δW[9]        → для =0: δW[9] = σ₁(δW[14])
  δW[17] = σ₁(δW[15]) ⊕ δW[10]       → для =0: δW[10] = σ₁(δW[15])
  δW[18] = σ₁(0) ⊕ δW[11] ⊕ σ₀(δW[3]) → для =0: δW[11] = σ₀(δW[3])

ТРИ уравнения — ТРИ определённых слова!
Остальные δW[4..8,12..15] — СВОБОДНЫ для state repair.

Fabric repair = schedule lock (GF2) + state repair (свободные слова).
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

def state_diff(s1, s2):
    return sum(hw(s1[i] ^ s2[i]) for i in range(8))

def a_repair_W(state2, target_a, r_idx):
    a, b, c, d, e, f, g, h = state2
    T2 = add32(Sigma0(a), Maj(a, b, c))
    T1_needed = sub32(target_a, T2)
    return sub32(sub32(sub32(sub32(T1_needed, h), Sigma1(e)), Ch(e, f, g)), K[r_idx])


def fabric_repair(W1, break_bit=31):
    """
    Fabric repair: schedule lock (GF2) + a-repair (свободные слова).

    Step 1: Вычисляем GF2-решение для δW[9,10,11] чтобы δW[16,17,18]=0
    Step 2: Используем свободные δW[4..8,12..15] для a-repair (δa=0)
    Step 3: Проверяем carry-шум
    """
    W1f = expand_schedule(W1)

    # Forward trace 1
    s1 = tuple(IV)
    states1 = [s1]
    for r in range(64):
        s1 = R(s1, W1f[r], r)
        states1.append(s1)

    # Break
    W2 = list(W1)
    W2[3] ^= (1 << break_bit)

    # STEP 1: GF(2) schedule solution
    # δW[0..2] = 0, δW[3] = (1 << break_bit)
    # Пока оставляем δW[4..8,12..15] = 0 (будем менять для a-repair)
    # Фиксируем δW[14,15] = 0 (свободный выбор)
    # Тогда:
    #   δW[9] = σ₁(δW[14]) = σ₁(0) = 0
    #   δW[10] = σ₁(δW[15]) = σ₁(0) = 0
    #   δW[11] = σ₀(δW[3]) = σ₀(1 << break_bit)

    delta_W3 = 1 << break_bit
    gf2_dw9 = 0   # σ₁(δW[14]) with δW[14]=0
    gf2_dw10 = 0  # σ₁(δW[15]) with δW[15]=0
    gf2_dw11 = sigma0(delta_W3)  # σ₀(δW[3])

    # W2 с GF2-оптимальными значениями
    # δW[9] = gf2_dw9, значит W2[9] = W1[9] ^ gf2_dw9
    W2[9] = W1[9] ^ gf2_dw9
    W2[10] = W1[10] ^ gf2_dw10
    W2[11] = W1[11] ^ gf2_dw11

    # δW[4..8,12..15] = используем для a-repair
    # Сначала: forward trace 2 до break
    s2 = tuple(IV)
    for r in range(3):
        s2 = R(s2, W2[r], r)
    s2 = R(s2, W2[3], 3)  # break

    # A-repair на свободных словах: r=4..8 (5 слов) + r=12..15 (4 слова)
    # r=9,10,11 — заняты GF2-решением
    for r in range(4, 9):  # свободные для a-repair
        W2[r] = a_repair_W(s2, states1[r + 1][0], r)
        s2 = R(s2, W2[r], r)

    # r=9,10,11 — GF2-фиксированные
    for r in [9, 10, 11]:
        s2 = R(s2, W2[r], r)

    # r=12..15 — свободные для a-repair
    for r in range(12, 16):
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

    dw = [hw(W1f[r] ^ W2f[r]) for r in range(64)]
    H1 = tuple(add32(IV[i], states1[64][i]) for i in range(8))
    H2 = tuple(add32(IV[i], s2_full[i]) for i in range(8))
    dh = sum(hw(H1[i] ^ H2[i]) for i in range(8))

    return profile, dw, dh, W2


def main():
    np.random.seed(42)

    print("=" * 70)
    print("FABRIC REPAIR: schedule lock + state repair")
    print("=" * 70)

    # =================================================================
    print(f"\n{'=' * 70}")
    print("1. GF(2) SCHEDULE SOLUTION")
    print("=" * 70)

    delta_W3 = 1 << 31
    gf2_dw11 = sigma0(delta_W3)
    print(f"\n  Break: δW[3] = 2^31 = {delta_W3:#010x}")
    print(f"  GF(2) solution:")
    print(f"    δW[9]  = σ₁(δW[14]) = σ₁(0) = 0")
    print(f"    δW[10] = σ₁(δW[15]) = σ₁(0) = 0")
    print(f"    δW[11] = σ₀(δW[3])  = σ₀(2^31) = {gf2_dw11:#010x} (HW={hw(gf2_dw11)})")
    print(f"\n  Это ГАРАНТИРУЕТ (в GF2): δW[16]=0, δW[17]=0, δW[18]=0")
    print(f"  Carry-шум может добавить ошибку.")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("2. FABRIC REPAIR vs A-REPAIR: сравнение")
    print("=" * 70)

    N = 5000
    fr_profiles = []
    fr_dws = []
    fr_dhs = []
    fr_meeting = []
    fr_locks = []

    for trial in range(N):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        profile, dw, dh, W2 = fabric_repair(W1)
        fr_profiles.append(profile)
        fr_dws.append(dw)
        fr_dhs.append(dh)

        # Meeting: longest consecutive BOTH=0
        both_chain = 0
        max_chain = 0
        for r in range(64):
            if profile[r] == 0:
                both_chain += 1
                max_chain = max(max_chain, both_chain)
            else:
                both_chain = 0
        fr_meeting.append(max_chain)

        # Schedule lock check
        lock_from = 64
        for r in range(16, 64):
            if dw[r] != 0:
                lock_from = r
                break
        fr_locks.append(lock_from)

    mean_profile = np.mean(fr_profiles, axis=0)
    mean_dw = np.mean(fr_dws, axis=0)
    both_pct = [np.mean([p[r] == 0 for p in fr_profiles]) * 100 for r in range(64)]

    print(f"\n  Fabric repair profile ({N} trials):")
    print(f"  {'r':>4} {'δstate':>8} {'BOTH%':>7} {'δW':>6}")
    for r in range(25):
        marker = ""
        if both_pct[r] > 99: marker = " ★ZERO"
        elif both_pct[r] > 10: marker = " ~pass"
        print(f"  {r:4d} {mean_profile[r]:7.1f} {both_pct[r]:6.1f}% {mean_dw[r]:5.1f}{marker}")

    print(f"\n  Schedule δW[16..20]:")
    for r in range(16, 25):
        p0 = sum(1 for d in fr_dws if d[r] == 0) / N * 100
        print(f"    δW[{r}]: {mean_dw[r]:.2f} бит, P(=0)={p0:.1f}%")

    print(f"\n  Meeting chain: mean={np.mean(fr_meeting):.1f}, max={max(fr_meeting)}")
    print(f"  HW(δH): mean={np.mean(fr_dhs):.1f}, min={min(fr_dhs)}")

    # Schedule lock
    lock_16 = sum(1 for l in fr_locks if l > 18) / N * 100
    lock_all = sum(1 for l in fr_locks if l >= 64) / N * 100
    print(f"\n  Schedule lock from r=16: P(δW[16..18]=0) = {lock_16:.2f}%")
    print(f"  Full schedule lock: P(δW[16..63]=0) = {lock_all:.2f}%")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("3. СРАВНЕНИЕ: fabric repair vs a-repair vs Wang")
    print("=" * 70)

    # A-repair baseline
    ar_dw16 = []
    ar_dw17 = []
    ar_dw18 = []
    ar_meeting = []

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
        for r in range(4, 16):
            W2[r] = a_repair_W(s2, states1[r + 1][0], r)
            s2 = R(s2, W2[r], r)

        W2f = expand_schedule(W2)
        ar_dw16.append(hw(W1f[16] ^ W2f[16]))
        ar_dw17.append(hw(W1f[17] ^ W2f[17]))
        ar_dw18.append(hw(W1f[18] ^ W2f[18]))

    print(f"\n  {'':>20} {'δW[16]':>8} {'δW[17]':>8} {'δW[18]':>8} {'P(18=0)':>9}")
    print(f"  {'A-repair':>20} {np.mean(ar_dw16):7.2f} {np.mean(ar_dw17):7.2f} {np.mean(ar_dw18):7.2f} "
          f"{sum(1 for x in ar_dw18 if x==0)/N*100:8.2f}%")
    print(f"  {'Fabric repair':>20} {mean_dw[16]:7.2f} {mean_dw[17]:7.2f} {mean_dw[18]:7.2f} "
          f"{sum(1 for d in fr_dws if d[18]==0)/N*100:8.2f}%")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("4. REBOOT: работает ли при fabric repair?")
    print("=" * 70)

    # Fabric repair: a-repair на r=4..8 и r=12..15, пропуск r=9,10,11
    # r=9,10,11 НЕ a-repaired → δa может быть ≠ 0 на этих раундах
    # → reboot может НЕ работать

    reboot_count = sum(1 for p in fr_profiles if p[11] == 0)
    meeting_12 = sum(1 for p in fr_profiles if p[11] == 0)

    print(f"\n  Reboot на r=12: {reboot_count}/{N} = {reboot_count/N*100:.1f}%")
    print(f"  (A-repair: 100%)")

    if reboot_count < N * 0.5:
        print(f"\n  ⚠ Reboot НАРУШЕН! Пропуск a-repair на r=9,10,11 убивает convergence.")
        print(f"  GF2-фиксация δW[9,10,11] конфликтует с δa=0 на этих раундах.")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("5. TRADE-OFF: schedule lock vs reboot")
    print("=" * 70)

    print(f"""
  A-REPAIR:
    ✓ Reboot 100% (δa=0 на r=4..15)
    ✓ Meeting r=12 бесплатен
    × δW[16]=0.5 бит, δW[17]=0.5 бит, δW[18]=7.7 бит
    × Schedule lock: 0%

  FABRIC REPAIR:
    {'✓' if reboot_count > N*0.9 else '×'} Reboot {reboot_count/N*100:.0f}%
    ✓ δW[16]={mean_dw[16]:.1f}, δW[17]={mean_dw[17]:.1f}, δW[18]={mean_dw[18]:.1f}
    {'✓' if lock_16 > 1 else '×'} Schedule lock r=16..18: {lock_16:.1f}%

  ВОПРОС: можно ли иметь ОБА?
    GF2-решение фиксирует δW[9,10,11].
    A-repair требует δW[9,10,11] для δa=0.
    Если GF2-значения СЛУЧАЙНО совпадают с a-repair →
    reboot + schedule lock.
    P(совпадение) = ???
""")

    # Проверяем: в скольких случаях GF2 и a-repair дают одинаковые δW[9,10,11]?
    match_count = 0
    for trial in range(N):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W1f = expand_schedule(W1)
        s1 = tuple(IV)
        states1 = [s1]
        for r in range(64):
            s1 = R(s1, W1f[r], r)
            states1.append(s1)

        # A-repair δW values
        W2_ar = list(W1)
        W2_ar[3] ^= (1 << 31)
        s2 = tuple(IV)
        for r in range(3):
            s2 = R(s2, W2_ar[r], r)
        s2 = R(s2, W2_ar[3], 3)
        for r in range(4, 16):
            W2_ar[r] = a_repair_W(s2, states1[r + 1][0], r)
            s2 = R(s2, W2_ar[r], r)

        ar_dw9 = W1[9] ^ W2_ar[9]
        ar_dw10 = W1[10] ^ W2_ar[10]
        ar_dw11 = W1[11] ^ W2_ar[11]

        gf2_dw9_needed = 0  # σ₁(δW[14]) with δW[14]=0 in a-repair → actually δW[14] from a-repair
        # Hmm, need to compute properly. δW[14] from a-repair = W1[14] ^ W2_ar[14]
        ar_dw14 = W1[14] ^ W2_ar[14]
        ar_dw15 = W1[15] ^ W2_ar[15]

        gf2_need_dw9 = sigma1(ar_dw14)   # need δW[9] = σ₁(δW[14]) for δW[16]=0 in GF2
        gf2_need_dw10 = sigma1(ar_dw15)  # need δW[10] = σ₁(δW[15])
        gf2_need_dw11 = sigma0(1 << 31)  # need δW[11] = σ₀(δW[3])

        if ar_dw9 == gf2_need_dw9 and ar_dw10 == gf2_need_dw10 and ar_dw11 == gf2_need_dw11:
            match_count += 1

    print(f"\n  P(a-repair δW[9,10,11] = GF2 solution): {match_count}/{N} = {match_count/N*100:.3f}%")
    if match_count > 0:
        print(f"  → СОВПАДЕНИЕ НАЙДЕНО! Reboot + schedule lock возможны одновременно!")
        print(f"  → Стоимость: 2^{-np.log2(match_count/N):.1f}")
    else:
        print(f"  → 0 совпадений из {N}. Нужно >2^{np.log2(N):.0f} попыток.")


if __name__ == "__main__":
    main()
