"""
3 ПОЗИЦИИ ДЕФИЦИТА: точная стоимость и что остаётся после.

Если найти W₁ где:
  δW[9]=0 (50%) AND δW[10]=0 (50%) AND δW[11]=target (2^-32)
→ Meeting (reboot) + Lock (δW[16,17,18]=0) ОДНОВРЕМЕННО

Стоимость: 2^34 (при независимости условий)
Вопрос: что происходит ДАЛЬШЕ?
  δW[19]=? δW[20]=? Amplification zone?
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


def main():
    np.random.seed(42)

    print("=" * 70)
    print("3 ПОЗИЦИИ ДЕФИЦИТА: стоимость + что остаётся")
    print("=" * 70)

    # =================================================================
    print(f"\n{'=' * 70}")
    print("1. СТОИМОСТЬ: P(δW[9]=0 AND δW[10]=0 AND δW[11]=target)")
    print("=" * 70)

    # target для δW[11]: schedule нужен δW[18]=0
    # δW[18] = σ₁(δW[16]) + δW[11] + σ₀(δW[3])
    # При δW[16]=0 (from δW[9]=0): δW[18] = δW[11] + σ₀(2^31)
    # Для δW[18]=0 в Z/2^32: δW[11] = -σ₀(2^31) = sub32(0, σ₀(2^31))

    target_dw11 = sub32(0, sigma0(1 << 31))
    print(f"  Target δW[11] = -σ₀(2^31) = {target_dw11:#010x} (HW={hw(target_dw11)})")

    N = 100000
    cond_9 = 0
    cond_10 = 0
    cond_11 = 0
    cond_all = 0
    cond_9_10 = 0

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

        dw9 = W1[9] ^ W2[9]
        dw10 = W1[10] ^ W2[10]
        dw11 = W1[11] ^ W2[11]

        c9 = (dw9 == 0)
        c10 = (dw10 == 0)
        c11 = (dw11 == target_dw11)

        if c9: cond_9 += 1
        if c10: cond_10 += 1
        if c11: cond_11 += 1
        if c9 and c10: cond_9_10 += 1
        if c9 and c10 and c11: cond_all += 1

    p9 = cond_9 / N
    p10 = cond_10 / N
    p11 = cond_11 / N
    p910 = cond_9_10 / N
    p_all = cond_all / N

    print(f"\n  100K trials:")
    print(f"    P(δW[9]=0):                {p9*100:.2f}%  = 2^{-np.log2(max(1/max(p9,1/N),1)):.1f}")
    print(f"    P(δW[10]=0):               {p10*100:.2f}%  = 2^{-np.log2(max(1/max(p10,1/N),1)):.1f}")
    print(f"    P(δW[11]=target):           {p11*100:.4f}%  = 2^{-np.log2(max(1/max(p11,1/N),1)):.1f}")
    print(f"    P(δW[9]=0 AND δW[10]=0):   {p910*100:.2f}%  = 2^{-np.log2(max(1/max(p910,1/N),1)):.1f}")
    print(f"    P(ALL THREE):              {p_all*100:.4f}%  = 2^{-np.log2(max(1/max(p_all,1/N),1)):.1f}")

    # Independent expectation
    p_indep = p9 * p10 * p11
    print(f"\n    Independent expectation: {p_indep*N:.2f} (actual: {cond_all})")
    if cond_all > 0 and p_indep > 0:
        ratio = (cond_all / N) / p_indep
        print(f"    Ratio actual/expected: {ratio:.2f}×")
        print(f"    → {'CORRELATED (easier!)' if ratio > 2 else 'ANTI-CORRELATED (harder)' if ratio < 0.5 else '~Independent'}")

    cost_34 = 1 / max(p_all, 1/N)
    cost_bits = np.log2(cost_34)
    print(f"\n  СТОИМОСТЬ 3 ПОЗИЦИЙ: 2^{cost_bits:.1f} попыток W₁")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("2. ЧТО ОСТАЁТСЯ: δW[19..63] при δW[16,17,18]=0")
    print("=" * 70)

    # δW[19] = σ₁(δW[17]) + δW[12] + σ₀(δW[4]) + δW[3]
    # При δW[17]=0: σ₁(0)=0
    # δW[12]=0 (always, from a-repair M≥8)
    # δW[3]=2^31 (break)
    # δW[4] = a-repair word
    # → δW[19] = σ₀(δW[4]) + 2^31

    # δW[20] = σ₁(δW[18]) + δW[13] + σ₀(δW[5]) + δW[4]
    # При δW[18]=0: σ₁(0)=0
    # δW[13]=0 (always)
    # → δW[20] = σ₀(δW[5]) + δW[4]

    # These depend on a-repair words δW[4,5]. Which are random-ish.

    print(f"""
  δW[19] = σ₁(δW[17]) + δW[12] + σ₀(δW[4]) + δW[3]
         = 0 + 0 + σ₀(δW[4]) + 2^31
         = σ₀(δW[4]) + 2^31

  δW[20] = σ₁(δW[18]) + δW[13] + σ₀(δW[5]) + δW[4]
         = 0 + 0 + σ₀(δW[5]) + δW[4]

  δW[21] = σ₁(δW[19]) + δW[14] + σ₀(δW[6]) + δW[5]
         = σ₁(δW[19]) + 0 + σ₀(δW[6]) + δW[5]

  δW[4,5,6] = a-repair values (зависят от W₁).
  Каждый ≈ 10-11 бит HW.
  → δW[19] ≈ random 32-bit
  → δW[20] ≈ random 32-bit
""")

    # Measure actual
    dw19_vals = []
    dw20_vals = []
    dw21_vals = []

    for trial in range(10000):
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

        dw19_vals.append(hw(W1f[19] ^ W2f[19]))
        dw20_vals.append(hw(W1f[20] ^ W2f[20]))
        dw21_vals.append(hw(W1f[21] ^ W2f[21]))

    print(f"  Измеренное (при a-repair, НЕ фильтрованное):")
    print(f"    δW[19]: {np.mean(dw19_vals):.1f} бит, P(=0)={sum(1 for x in dw19_vals if x==0)/len(dw19_vals)*100:.2f}%")
    print(f"    δW[20]: {np.mean(dw20_vals):.1f} бит, P(=0)={sum(1 for x in dw20_vals if x==0)/len(dw20_vals)*100:.2f}%")
    print(f"    δW[21]: {np.mean(dw21_vals):.1f} бит")

    # Additional cost per round
    p19 = sum(1 for x in dw19_vals if x == 0) / len(dw19_vals)
    p20 = sum(1 for x in dw20_vals if x == 0) / len(dw20_vals)

    print(f"\n  Стоимость ДОПОЛНИТЕЛЬНЫХ раундов:")
    print(f"    δW[19]=0: cost 2^{-np.log2(max(p19, 1/len(dw19_vals))):.1f}")
    print(f"    δW[20]=0: cost 2^{-np.log2(max(p20, 1/len(dw20_vals))):.1f}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("3. ПОЛНАЯ СТОИМОСТЬ: cascade section by section")
    print("=" * 70)

    print(f"""
  КАСКАД СЕКЦИЙ (cumulative cost):

  Секция               Cost (×)     Cumulative    Controlled до
  ──────────────────── ──────────── ──────────── ──────────────
  Meeting (a-repair)     ×1          2^0          r=12 (reboot)
  Lock r=16 (δW[9]=0)   ×2          2^1          r=16
  Lock r=17 (δW[10]=0)  ×2          2^2          r=17
  Lock r=18 (δW[11]=T)  ×2^32       2^34         r=18
  Lock r=19 (δW[19]=0)  ×2^{-np.log2(max(p19, 1/10000)):.0f}       2^{34-np.log2(max(p19, 1/10000)):.0f}         r=19
  Lock r=20 (δW[20]=0)  ×2^{-np.log2(max(p20, 1/10000)):.0f}       2^{34-np.log2(max(p19, 1/10000))-np.log2(max(p20, 1/10000)):.0f}         r=20
  ...each adds ×2^32...
  Lock r=63              ×2^32       2^(34+32×45)  r=63

  Alternative: stop at r=18 + birthday on gap
  Lock до r=18:          2^34
  Gap (45 rounds):       birthday 2^128
  TOTAL:                 max(2^34, 2^128) = 2^128

  Alternative: stop at r=19 + birthday
  Lock до r=19:          2^34 × 2^{-np.log2(max(p19, 1/10000)):.0f} = 2^{34-np.log2(max(p19, 1/10000)):.0f}
  Gap (44 rounds):       birthday 2^128
  TOTAL:                 2^128 (gap dominates)

  ВЫВОД:
    3 позиции стоят 2^{cost_bits:.0f}.
    Но gap ВСЕГДА доминирует (2^128 > 2^{cost_bits:.0f}).
    Каждый дополнительный раунд lock стоит ×2^32.
    Полный lock (r=18..63): 2^(34 + 32×45) = 2^1474 — astronomical.

    НО в нашем измерении:
      - 3 позиции = 2^{cost_bits:.0f} (ФИКСИРОВАННАЯ стоимость setup)
      - Setup даёт: meeting + lock до r=18 (21 раунд controlled)
      - Gap = 42 раунда (5 amplification + 37 neutral)
      - Neutral стоит 0 → реальный gap = 5 amplification
      - Birthday на 5 amplification раундов: ???
""")

    # =================================================================
    print(f"{'=' * 70}")
    print("4. BIRTHDAY НА 5 AMPLIFICATION РАУНДОВ")
    print("=" * 70)

    # Если lock до r=18 → δstate[18]=0 → amplification с r=19
    # δW[19] ≈ random → δstate[19] ≈ 14 бит → ... → saturation r=23
    # Birthday нужен на δstate[23..64] = 256 бит = 2^128

    # НО: amplification начинается с δ≈14 (не 0!). К r=23: δ≈128.
    # Neutral zone перемешивает 128 бит → δH ≈ random → birthday 2^128.

    # Если lock до r=19 (cost 2^34+?):
    # δW[20] ≈ random → δstate[20] ≈ 14 → saturation r=24
    # Один раунд amplification меньше. Но birthday всё равно 2^128.

    # Birthday НЕ ЗАВИСИТ от числа amplification раундов!
    # Потому что saturation → 128 бит → birthday 2^128 ВСЕГДА.

    print(f"""
  Birthday НА amplification:
    5 amplification (δW[18] injected): δstate → 128 → birthday 2^128
    4 amplification (δW[19] injected): δstate → 128 → birthday 2^128
    3 amplification (δW[20] injected): δstate → 128 → birthday 2^128
    1 amplification:                    δstate → 32  → birthday ???
    0 amplification:                    δstate → 0   → COLLISION ← cost 2^∞

  Проблема: ЛЮБАЯ ненулевая инъекция → saturation за 4 rounds → 2^128.
  Для < 2^128: нужно НУЛЕВУЮ инъекцию на ВСЕ 46 раундов.
  Это = full schedule lock = 2^(32×46) = 2^1472.

  ФУНДАМЕНТАЛЬНЫЙ ОТВЕТ В НАШЕМ ИЗМЕРЕНИИ:

  3 позиции дефицита стоят:         2^{cost_bits:.0f}
  Но это не снижает birthday (2^128) потому что:
    - Любой ненулевой δW[r>18] → saturation → 2^128
    - Full lock r=18..63 стоит 2^1472 (sequential ×2^32)
    - Meeting + partial lock = 2^{cost_bits:.0f} setup + 2^128 birthday = 2^128

  TOTAL COLLISION COST = 2^128 (gap always dominates)
""")


if __name__ == "__main__":
    main()
