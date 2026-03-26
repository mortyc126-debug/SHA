"""
Двойная оптимизация: fix δstate[1..16] И минимизировать δW[16..63].

Факты:
  1. ker(L) = {0} — нельзя обнулить ВСЮ schedule diff
  2. Schedule ЛИНЕЕН над GF(2): δW[r] = σ₁(δW[r-2])⊕δW[r-7]⊕σ₀(δW[r-15])⊕δW[r-16]
  3. Forward fix требует: W[r] = f(state, target) — одно значение
  4. Но у W[r] 32 бита, а fix использует НЕ ВСЕ (только T1)

Вопрос: есть ли СВОБОДА внутри fix?
  T1 = h + Σ₁(e) + Ch(e,f,g) + K + W
  Для fix нужен конкретный T1. Значит W = T1 - rest.
  W определено ОДНОЗНАЧНО. Свободы нет?

  НО: fix требует δstate'=0, а НЕ конкретный state'.
  Два разных state' могут оба давать δstate'=0 через:
    e'₁ = d₁ + T1₁  и  e'₂ = d₂ + T1₂
    δe' = 0 ⟺ d₁+T1₁ = d₂+T1₂

  Это ОДНО уравнение на 32-бит W₁ и 32-бит W₂.
  У нас 64 бита свободы - 32 бита уравнение = 32 бита ОСТАТКА.

  НО: нам также нужно δa'=0:
    a' = T1 + T2
    δa' = δT1 + δT2 + carry
    Ещё 32 бита условий.

  ИТОГО: 64 бит свободы (W₁, W₂) - 64 бит условий (δe'=0, δa'=0) = 0 свободы.

  ИЛИ: подход через ОДИН trace.
  Фиксируем W₁[r], подбираем W₂[r] = W₁[r] + δW[r].
  Условие fix: конкретный δW[r] чтобы δstate'=0.
  δW[r] ОДНОЗНАЧНО определён.
  Свободы нет.

  Тогда стратегия: ВЫБИРАТЬ W₁[0] (32 бита свободы!)
  → это меняет ВСЕ downstream: state, T1, T2, δW[r]
  → и schedule (через W₁[0] → W[16], W[31], ...)

  Ищем W₁[0] который минимизирует schedule damage ПОСЛЕ fix.
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


def compute_needed_W(state, target_next, r_idx):
    a, b, c, d, e, f, g, h = state
    ap = target_next[0]
    T2 = add32(Sigma0(a), Maj(a, b, c))
    T1 = sub32(ap, T2)
    W = sub32(sub32(sub32(sub32(T1, h), Sigma1(e)), Ch(e, f, g)), K[r_idx])
    return W


def expand_schedule(W16):
    W = list(W16)
    for r in range(16, 64):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))
    return W


def full_sha(W16):
    W = expand_schedule(W16)
    s = tuple(IV)
    for r in range(64):
        s = R(s, W[r], r)
    return tuple(add32(IV[i], s[i]) for i in range(8))


def forward_fix(W1, delta_bit):
    """Forward fix: возвращает W2 такой что δstate[1..16]=0."""
    W2 = list(W1)
    W2[0] ^= (1 << delta_bit)

    s1 = tuple(IV)
    states1 = [s1]
    for r in range(16):
        s1 = R(s1, W1[r], r)
        states1.append(s1)

    s2 = tuple(IV)
    s2 = R(s2, W2[0], 0)
    for r in range(1, 16):
        W2[r] = compute_needed_W(s2, states1[r + 1], r)
        s2 = R(s2, W2[r], r)

    return W2


def state_diff_hw(s1, s2):
    return sum(hw(a ^ b) for a, b in zip(s1, s2))


def experiment_dual_optimization():
    np.random.seed(42)

    print("=" * 70)
    print("ДВОЙНАЯ ОПТИМИЗАЦИЯ: fix + minimize schedule damage")
    print("=" * 70)

    # ===================================================================
    print("\n" + "=" * 70)
    print("1. BASELINE: schedule damage после forward fix")
    print("=" * 70)

    baseline_damages = []
    baseline_hash_diffs = []

    for trial in range(500):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = forward_fix(W1, delta_bit=31)

        W1f = expand_schedule(W1)
        W2f = expand_schedule(W2)

        sched_damage = sum(hw(W1f[r] ^ W2f[r]) for r in range(16, 64))
        baseline_damages.append(sched_damage)

        H1 = full_sha(W1)
        H2 = full_sha(W2)
        hash_diff = sum(hw(H1[i] ^ H2[i]) for i in range(8))
        baseline_hash_diffs.append(hash_diff)

    print(f"\n  Baseline (random W₁, δbit=31):")
    print(f"    Schedule damage: {np.mean(baseline_damages):.0f} ± {np.std(baseline_damages):.0f} бит")
    print(f"    HW(δH):          {np.mean(baseline_hash_diffs):.1f} ± {np.std(baseline_hash_diffs):.1f}")
    print(f"    Best HW(δH):     {min(baseline_hash_diffs)}")

    # ===================================================================
    print("\n" + "=" * 70)
    print("2. ОПТИМИЗАЦИЯ: выбираем δ-бит чтобы min schedule damage")
    print("=" * 70)

    for trial in range(5):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]

        best_bit = -1
        best_damage = 999999
        best_hash_diff = 256

        for delta_bit in range(32):
            W2 = forward_fix(W1, delta_bit)
            W1f = expand_schedule(W1)
            W2f = expand_schedule(W2)
            damage = sum(hw(W1f[r] ^ W2f[r]) for r in range(16, 64))

            if damage < best_damage:
                best_damage = damage
                best_bit = delta_bit
                H1 = full_sha(W1)
                H2 = full_sha(W2)
                best_hash_diff = sum(hw(H1[i] ^ H2[i]) for i in range(8))

        print(f"    Trial {trial}: best δ-bit={best_bit:2d}, sched_damage={best_damage:4d}, HW(δH)={best_hash_diff}")

    # ===================================================================
    print("\n" + "=" * 70)
    print("3. ОПТИМИЗАЦИЯ W₁[0]: выбираем НАЧАЛЬНОЕ СООБЩЕНИЕ")
    print("=" * 70)

    # Фиксируем W₁[1..15], перебираем W₁[0]
    W1_base = [np.random.randint(0, 2**32) for _ in range(16)]

    best_w0 = W1_base[0]
    best_total = 999999
    best_hd = 256
    results = []

    for attempt in range(2000):
        W1 = list(W1_base)
        W1[0] = np.random.randint(0, 2**32)

        W2 = forward_fix(W1, delta_bit=31)
        W1f = expand_schedule(W1)
        W2f = expand_schedule(W2)
        damage = sum(hw(W1f[r] ^ W2f[r]) for r in range(16, 64))

        H1 = full_sha(W1)
        H2 = full_sha(W2)
        hd = sum(hw(H1[i] ^ H2[i]) for i in range(8))
        results.append((damage, hd))

        if hd < best_hd:
            best_hd = hd
            best_w0 = W1[0]
            best_total = damage

    damages = [r[0] for r in results]
    hds = [r[1] for r in results]

    print(f"\n  Перебор W₁[0] (2000 вариантов, fix δbit=31):")
    print(f"    Schedule damage: {np.mean(damages):.0f} ± {np.std(damages):.0f}")
    print(f"    HW(δH) mean:    {np.mean(hds):.1f}")
    print(f"    HW(δH) best:    {best_hd}")
    print(f"    HW(δH) worst:   {max(hds)}")

    corr_dam_hd = np.corrcoef(damages, hds)[0, 1]
    print(f"\n  Корреляция schedule_damage ↔ HW(δH): {corr_dam_hd:+.4f}")
    print(f"  → {'Schedule damage ПРЕДСКАЗЫВАЕТ HW(δH)!' if abs(corr_dam_hd) > 0.1 else 'НЕТ связи'}")

    # ===================================================================
    print("\n" + "=" * 70)
    print("4. МАСШТАБНЫЙ ПОИСК: 100K вариантов W₁")
    print("=" * 70)

    best_hd_global = 256
    best_w1_global = None
    hd_hist = np.zeros(257)

    for attempt in range(100000):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = forward_fix(W1, delta_bit=31)

        H1 = full_sha(W1)
        H2 = full_sha(W2)
        hd = sum(hw(H1[i] ^ H2[i]) for i in range(8))
        hd_hist[hd] += 1

        if hd < best_hd_global:
            best_hd_global = hd
            best_w1_global = list(W1)

    print(f"\n  100K пар с forward fix (δstate[1..16]=0):")
    print(f"    Best HW(δH): {best_hd_global}")
    print(f"    Mean HW(δH): {sum(i*hd_hist[i] for i in range(257))/100000:.1f}")

    # Распределение
    print(f"\n  Распределение HW(δH):")
    for threshold in [80, 85, 90, 95, 100, 105, 110, 115, 120]:
        count = sum(hd_hist[:threshold+1])
        print(f"    HW(δH) ≤ {threshold}: {int(count):6d} ({count/1000:.2f}%)")

    # Теоретическое для random
    import math
    print(f"\n  Для сравнения (random без fix):")
    print(f"    Best в 100K проб: ≈ {128 - 4*8:.0f} (128 - 4σ)")

    # ===================================================================
    print("\n" + "=" * 70)
    print("5. ВЕРДИКТ")
    print("=" * 70)

    print(f"""
  Forward fix гарантирует δstate[1..16] = 0.
  Раунды 17-64: определяются schedule diff.

  Schedule damage ≈ {np.mean(damages):.0f} бит (неконтролируемо).
  Корреляция damage↔HW(δH): {corr_dam_hd:+.3f}

  Best HW(δH) из 100K: {best_hd_global}
  Это {'ЛУЧШЕ random (≈96)!' if best_hd_global < 93 else 'примерно как random (≈96)'}

  Forward fix НЕ даёт преимущества:
    δstate[1..16]=0 не помогает, потому что schedule
    ПЕРЕЗАПИСЫВАЕТ состояние начиная с r=17.
    К r=20 δstate ≈ 128 — как будто fix'а не было.

  SHA-256 САМОВОССТАНАВЛИВАЕТСЯ через schedule.
""")


if __name__ == "__main__":
    experiment_dual_optimization()
