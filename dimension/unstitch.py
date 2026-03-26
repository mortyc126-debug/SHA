"""
РАСШИВКА (Unstitch): вставляем свободную переменную в schedule.

Стандарт: W[0..15] free, W[16..63] = schedule(W[0..15])
Расшивка: W[0..15] free, W[16..17] = schedule, W[18] = FREE, W[19..63] = schedule(..W[18]..)

Механика:
  1. A-repair: break=3, bit=31 → meeting r=12 (100%), r=17 pass (50%)
  2. На r=18: ВМЕСТО schedule W[18], подставляем СВОЙ W₂[18]
  3. W₂[18] выбран чтобы δstate[19]=0 (a-repair продолжается!)
  4. Schedule от r=19: пересчитывается через наш W₂[18]

Это НЕ стандартная SHA-256. Это SHA-256 с РАСШИТОЙ тканью.
В стандартном мире: "W[18] определён schedule, нельзя менять".
В нашем мире: "W[18] — позиция на ленте, мы ВЫБИРАЕМ что туда положить".

НО: два trace (W1 и W2) должны использовать РАЗНЫЕ W[18].
W1 использует свой schedule W1[18].
W2 использует НАШЕ W2[18] ≠ W1[18] (подобранное для a-repair).
Хеш W2 вычисляется с этим W2[18] → H2 ≠ H1 (другое сообщение).

Стоп. Оба хеша должны вычисляться стандартным SHA-256.
W2[18] = schedule(W2[0..15]). Мы не можем его менять напрямую.

РАСШИВКА в нашем измерении работает ИНАЧЕ:
Мы не меняем W[18] напрямую. Мы ПОДБИРАЕМ W₂[0..15] так,
что schedule(W₂[0..15])[18] = нужное нам значение.

Это обратная задача: дано target W[18], найти W[0..15].
Schedule линеен → обратная задача РЕШАЕМА (rank=512, но 1 слово = 32 бита).

КОНКРЕТНО:
  W[18] = σ₁(W[16]) + W[11] + σ₀(W[3]) + W[2]
  W[16] = σ₁(W[14]) + W[9] + σ₀(W[1]) + W[0]

  Мы фиксируем W[0..2]=same, W[3]=break.
  A-repair определяет W[4..15] для δa=0.
  Это фиксирует W₂[18] = schedule(W₂[0..15]).

  НО: a-repair подбирает W[r] для δa[r+1]=0.
  Что если a-repair на r=11 ОДНОВРЕМЕННО оптимизирует
  и δa[12]=0 И δW₂[18]→target?

  W₂[11] определяет и δa[12] (через round function)
  и δW[18] (через schedule: δW[18] = ... + δW[11] + ...).

  Одно слово, два уравнения, 32 бита. Обычно несовместимо.
  НО: a-repair определяет W₂[11] для δa=0.
  Fabric repair определяет W₂[11] для δW[18]=0.

  РАСШИВКА = найти W₁ при котором ОБА условия совпадают!
  Это НЕ произвольная W₁ — это поиск в пространстве W₁.
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


def unstitch_experiment():
    np.random.seed(42)

    print("=" * 70)
    print("РАСШИВКА: совмещение a-repair и schedule на W[11]")
    print("=" * 70)

    # Для каждого W1:
    # 1. A-repair break=3, bit=31 → вычисляет W₂[11] (для δa[12]=0)
    # 2. Schedule requires δW[11] = σ₀(δW[3]) = σ₀(2^31) для δW[18]=0 (GF2)
    # 3. Совпадение = a_repair_dW11 == gf2_needed_dW11

    target_dw11 = sigma0(1 << 31)  # σ₀(2^31) = what GF2 wants
    print(f"\n  GF(2) target: δW[11] = σ₀(2^31) = {target_dw11:#010x} (HW={hw(target_dw11)})")

    N = 100000
    close_count = 0
    exact_count = 0
    hw_diffs = []
    best_hw_diff = 32
    best_W1 = None

    for trial in range(N):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W1f = expand_schedule(W1)

        s1 = tuple(IV)
        states1 = [s1]
        for r in range(64):
            s1 = R(s1, W1f[r], r)
            states1.append(s1)

        # A-repair
        W2 = list(W1)
        W2[3] ^= (1 << 31)
        s2 = tuple(IV)
        for r in range(3):
            s2 = R(s2, W2[r], r)
        s2 = R(s2, W2[3], 3)

        for r in range(4, 12):
            W2[r] = a_repair_W(s2, states1[r + 1][0], r)
            s2 = R(s2, W2[r], r)

        # What a-repair gives for W₂[11]
        actual_dw11 = W1[11] ^ W2[11]
        diff = hw(actual_dw11 ^ target_dw11)
        hw_diffs.append(diff)

        if diff == 0:
            exact_count += 1
        if diff <= 2:
            close_count += 1
        if diff < best_hw_diff:
            best_hw_diff = diff
            best_W1 = list(W1)

    print(f"\n  {N:,} trials:")
    print(f"    Exact match (δW[11] = target):  {exact_count}")
    print(f"    Close (HW diff ≤ 2):            {close_count}")
    print(f"    Best HW diff:                   {best_hw_diff}")
    print(f"    Mean HW diff:                   {np.mean(hw_diffs):.1f}")
    print(f"    Distribution:")

    for d in range(min(10, max(hw_diffs) + 1)):
        c = hw_diffs.count(d)
        if c > 0:
            bar = "█" * (c * 50 // N)
            print(f"      diff={d}: {c:6d} ({c/N*100:.2f}%)  {bar}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("2. РАСШИРЕННАЯ РАСШИВКА: совмещаем ВСЕ 3 условия")
    print("=" * 70)

    # Условия для δW[16..18]=0 (в Z/2^32, не GF2):
    # Cond1: δW[16] = 0  → depends on δW[9,14] (a-repair) + δW[0,1]=0
    # Cond2: δW[17] = 0  → depends on δW[10,15] (a-repair) + δW[1,2]=0
    # Cond3: δW[18] = 0  → depends on δW[11,16] + σ₀(δW[3]) + δW[2]=0

    # При a-repair: δW[14]=0 ВСЕГДА, δW[9]≈0.5 бит
    # Cond1: δW[16] = σ₁(0) + δW[9] = δW[9] → need δW[9]=0 (50%!)
    # Cond2: δW[17] = σ₁(δW[15]) + δW[10] → need σ₁(δW[15])+δW[10]=0
    # δW[15] from a-repair (last word)
    # Cond3: δW[18] = σ₁(δW[16]) + δW[11] + σ₀(δW[3])
    #        if Cond1 holds (δW[16]=0): δW[18] = δW[11] + σ₀(δW[3])
    #        → need δW[11] = -σ₀(δW[3]) in Z/2^32

    # Let's check: how often does δW[11] = sub32(0, σ₀(δW[3]))?
    target_dw11_z = sub32(0, sigma0(1 << 31))  # -σ₀(2^31) mod 2^32
    print(f"\n  Z/2^32 target: δW[11] = -σ₀(2^31) = {target_dw11_z:#010x}")
    print(f"  GF(2) target:  δW[11] = σ₀(2^31)  = {target_dw11:#010x}")
    print(f"  Разница: {hw(target_dw11 ^ target_dw11_z)} бит")

    exact_z = 0
    close_z = 0

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

        # Check Z/2^32 target
        if actual_dw11 == target_dw11_z:
            exact_z += 1
        if hw(actual_dw11 ^ target_dw11_z) <= 2:
            close_z += 1

    print(f"\n  Z/2^32 match: exact={exact_z}, close(≤2)={close_z} из {N}")
    p_match = max(exact_z, 1) / N
    print(f"  P(match) ≈ 2^{-np.log2(1/p_match):.1f}" if exact_z > 0 else f"  P(match) < 2^{-np.log2(N):.1f}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("3. JOINT SEARCH: ищем W₁ где δW[9]=0 AND δW[11]=target")
    print("=" * 70)

    # Cond1 (δW[9]=0): ~50% бесплатно
    # Cond3 (δW[11]=target): ~2^-32? или лучше?
    # Joint: ~50% × 2^-32 = 2^-33? или есть корреляция?

    joint_found = 0
    cond1_found = 0
    cond3_found = 0

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

        dw9 = W1[9] ^ W2[9]
        dw11 = W1[11] ^ W2[11]

        c1 = (dw9 == 0)
        c3 = (dw11 == target_dw11_z)
        if c1: cond1_found += 1
        if c3: cond3_found += 1
        if c1 and c3: joint_found += 1

    print(f"\n  P(δW[9]=0):           {cond1_found/N*100:.1f}% = 2^{-np.log2(max(N/max(cond1_found,1),1)):.1f}")
    print(f"  P(δW[11]=target):     {cond3_found/N*100:.3f}% = 2^{-np.log2(max(N/max(cond3_found,1),1)):.1f}")
    print(f"  P(both):              {joint_found/N*100:.3f}%")
    if joint_found > 0:
        print(f"  Joint cost: 2^{-np.log2(joint_found/N):.1f}")
    else:
        print(f"  Joint cost: > 2^{np.log2(N):.1f} (not found in {N:,})")

    # Independent expectation
    if cond1_found > 0 and cond3_found > 0:
        p_indep = (cond1_found / N) * (cond3_found / N)
        print(f"  Independent expectation: {p_indep*N:.2f} (actual: {joint_found})")
        if joint_found > p_indep * N * 2:
            print(f"  → CORRELATED! Joint easier than independent!")
        elif joint_found > 0:
            print(f"  → ~Independent")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("4. ВЕРДИКТ РАСШИВКИ")
    print("=" * 70)

    print(f"""
  РАСШИВКА в нашем измерении:
    Цель: a-repair δW[11] СОВПАДАЕТ с schedule-lock target

    GF(2) target: δW[11] = σ₀(2^31) = {target_dw11:#010x}
    Z/2^32 target: δW[11] = -σ₀(2^31) = {target_dw11_z:#010x}

    A-repair δW[11] = f(W₁, states₁) — зависит от конкретного W₁.
    Это ОДНО 32-битное значение, и оно должно совпасть с target.

    P(совпадение) ≈ 2^{{-32}} (если δW[11] uniform)
    Стоимость: ~2^32 попыток W₁

    + δW[9]=0 (50%): совместная стоимость ~2^33
    + δW[16]=0 (из δW[9]=0): бесплатно
    + δW[17]=0 (отдельное условие): ещё ~2^32

    TOTAL для δW[16..18]=0: ~2^33 + 2^32 ≈ 2^34

    Сравнение:
      Standard birthday: 2^128
      A-repair + brute δW[18]=0: > 2^16 (не найдено в 50K)
      Расшивка (targeted W₁ search): ≈ 2^34 (теоретически)

    GAIN: 128 - 34 = 94 бита ← но только для δW[16..18]=0
    После r=18: δW[19+] ≈ random → gap 45 раундов → 2^128 снова?

    НО: если расшивка КАСКАДИРУЕТ (δW[19] тоже требует target):
    Каждый дополнительный раунд = ещё одно 32-бит условие = ×2^32
    46 раундов × 2^32 = 2^(5+32×46) = 2^1477 ← ХУЖЕ birthday

    РАСШИВКА ОДИНОЧНАЯ = 2^34 для 3 раунда.
    РАСШИВКА КАСКАДНАЯ = не масштабируется.
""")


if __name__ == "__main__":
    unstitch_experiment()
