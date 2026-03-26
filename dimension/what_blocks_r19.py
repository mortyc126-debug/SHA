"""
ЧТО ИМЕННО БЛОКИРУЕТ R=19?

A-repair + schedule transparency:
  r=17: 19.2% проход (δW[16]=2 бит)
  r=18: 6.4% проход  (δW[17]=1.8 бит)
  r=19: 0% проход    (δW[18]=9.2 бит)  ← ЧТО ЗДЕСЬ?

Разбираем по компонентам в нашем измерении:
  1. Откуда приходит δ на r=19? (какие переходы его несут?)
  2. Через NODE или через PIPE?
  3. Из T1 или T2?
  4. Из δW или из δstate?
  5. Можно ли это ОБОЙТИ нативными средствами?
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

reg_names = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']


def main():
    np.random.seed(42)
    N = 1000

    print("=" * 70)
    print("ЧТО БЛОКИРУЕТ R=19? Анатомия в нашем измерении")
    print("=" * 70)

    # Собираем пары которые ПРОШЛИ r=17 (19.2%)
    passed_r17 = []
    all_pairs = []

    for trial in range(N * 50):  # нужно много для 19.2% pass rate
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W1f = expand_schedule(W1)

        s1 = tuple(IV)
        states1 = [s1]
        for r in range(64):
            s1 = R(s1, W1f[r], r)
            states1.append(s1)

        W2 = list(W1)
        W2[3] ^= (1 << (trial % 32))
        s2 = tuple(IV)
        for r in range(3):
            s2 = R(s2, W2[r], r)
        s2 = R(s2, W2[3], 3)
        for r in range(4, 16):
            W2[r] = a_repair_W(s2, states1[r + 1][0], r)
            s2 = R(s2, W2[r], r)

        W2f = expand_schedule(W2)

        # Check r=17: need to go through r=16 and r=17
        s2_r16 = R(s2, W2f[16], 16)
        s2_17 = R(s2_r16, W2f[17], 17)

        ds17 = sum(hw(states1[18][i] ^ s2_17[i]) for i in range(8))

        if ds17 == 0:
            # Continue to r=18, 19, 20
            s2_cont = s2_17
            decomp = []
            for r in range(18, 22):
                a2, b2, c2, d2, e2, f2, g2, h2 = s2_cont
                a1, b1, c1, d1, e1, f1, g1, h1 = states1[r]

                # δT1 components
                dh = hw(h2 ^ h1)
                dSig1 = hw(Sigma1(e2) ^ Sigma1(e1))
                dCh = hw(Ch(e2, f2, g2) ^ Ch(e1, f1, g1))
                dW = hw(W2f[r] ^ W1f[r])

                T1_1 = add32(add32(add32(add32(h1, Sigma1(e1)), Ch(e1, f1, g1)), K[r]), W1f[r])
                T1_2 = add32(add32(add32(add32(h2, Sigma1(e2)), Ch(e2, f2, g2)), K[r]), W2f[r])
                dT1 = hw(T1_1 ^ T1_2)

                T2_1 = add32(Sigma0(a1), Maj(a1, b1, c1))
                T2_2 = add32(Sigma0(a2), Maj(a2, b2, c2))
                dT2 = hw(T2_1 ^ T2_2)

                # δ по регистрам
                dreg = [hw(states1[r][i] ^ s2_cont[i]) for i in range(8)]

                s2_cont = R(s2_cont, W2f[r], r)
                s1_next = states1[r + 1]
                dreg_next = [hw(s1_next[i] ^ s2_cont[i]) for i in range(8)]
                ds_next = sum(dreg_next)

                decomp.append({
                    'r': r,
                    'dreg': dreg,
                    'dreg_next': dreg_next,
                    'ds_next': ds_next,
                    'dh': dh, 'dSig1': dSig1, 'dCh': dCh, 'dW': dW,
                    'dT1': dT1, 'dT2': dT2,
                })

            passed_r17.append(decomp)

        if len(passed_r17) >= 200:
            break

    print(f"\n  Собрано {len(passed_r17)} пар, прошедших r=17")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("1. ИСТОЧНИКИ δ НА r=18,19,20 (при r=17 BOTH=0)")
    print("=" * 70)

    for r_offset in range(4):  # r=18,19,20,21
        r = 18 + r_offset

        means = {}
        for key in ['dh', 'dSig1', 'dCh', 'dW', 'dT1', 'dT2', 'ds_next']:
            vals = [p[r_offset][key] for p in passed_r17]
            means[key] = np.mean(vals)

        dreg_means = np.mean([p[r_offset]['dreg'] for p in passed_r17], axis=0)
        dreg_next_means = np.mean([p[r_offset]['dreg_next'] for p in passed_r17], axis=0)

        print(f"\n  ─── Раунд {r} ───")
        print(f"  ВХОД δstate: {sum(dreg_means):.1f}")
        print(f"    δ по регистрам: ", end="")
        for i in range(8):
            marker = "█" if dreg_means[i] > 2 else "░" if dreg_means[i] > 0.5 else " "
            print(f"δ{reg_names[i]}={dreg_means[i]:.1f}{marker} ", end="")
        print()

        print(f"  КОМПОНЕНТЫ δT1:")
        print(f"    δh (pipe→NODE_e):     {means['dh']:.1f}")
        print(f"    δΣ₁(e):               {means['dSig1']:.1f}")
        print(f"    δCh(e,f,g):           {means['dCh']:.1f}")
        print(f"    δW[{r}] (schedule):    {means['dW']:.1f}  ← SCHEDULE")

        print(f"  РЕЗУЛЬТАТ:")
        print(f"    δT1: {means['dT1']:.1f}")
        print(f"    δT2: {means['dT2']:.1f}")
        print(f"    δstate[{r+1}]: {means['ds_next']:.1f}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("2. ДЕКОМПОЗИЦИЯ: что ВНОСИТ δ, что ПЕРЕДАЁТ?")
    print("=" * 70)

    # На r=18 (первый раунд после pass r=17):
    # state[18] = state[17] (оба = 0, pass!)
    # δ может прийти ТОЛЬКО от δW[18]

    print(f"""
  r=18 (первый раунд gap после pass r=17):
    δstate[17] = 0 (pass!)
    δstate[18] = f(δstate[17], δW[18])

    Поскольку δstate[17] = 0:
      δT2 = Σ₀(δa) + δMaj = Σ₀(0) + Maj(0,0,0) = 0
      δT1 = δh + δΣ₁(e) + δCh + K + δW[18]
           = 0 + 0 + 0 + 0 + δW[18]
           = δW[18]  ← ЕДИНСТВЕННЫЙ источник!

      δa[19] = δT1 + δT2 = δW[18] + 0 = δW[18]
      δe[19] = δd + δT1 = 0 + δW[18] = δW[18]

    δstate[19] определяется ПОЛНОСТЬЮ через δW[18].

  δW[18] = σ₁(δW[16]) + δW[11] + σ₀(δW[3]) + δW[2]
         = σ₁(δW[16]) + δW[11] + σ₀(δW[3]) + 0

    δW[2] = 0 (a-repair)
    δW[3] = 1 бит (break) → σ₀(δW[3]) = σ₀(small) ≈ 3 бита
    δW[11] = a-repair value (переменный)
    δW[16] = small (2 бит avg) → σ₁(δW[16]) ≈ 3 бита
""")

    # Измеряем реальный вклад каждого компонента δW[18]
    dw18_components = {'sig1_w16': [], 'w11': [], 'sig0_w3': [], 'w2': [], 'total': []}

    for trial in range(min(len(passed_r17) * 50, N * 10)):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W1)
        W2[3] ^= (1 << (trial % 32))

        s1 = tuple(IV)
        states1 = [s1]
        for r in range(20):
            s1 = R(s1, expand_schedule(W1)[r], r)
            states1.append(s1)

        s2 = tuple(IV)
        for r in range(3):
            s2 = R(s2, W2[r], r)
        s2 = R(s2, W2[3], 3)
        for r in range(4, 16):
            W2[r] = a_repair_W(s2, states1[r + 1][0], r)
            s2 = R(s2, W2[r], r)

        W1f = expand_schedule(W1)
        W2f = expand_schedule(W2)

        # Components of δW[18]
        c1 = hw(sigma1(W1f[16]) ^ sigma1(W2f[16]))  # σ₁(δW[16])
        c2 = hw(W1f[11] ^ W2f[11])                    # δW[11]
        c3 = hw(sigma0(W1f[3]) ^ sigma0(W2f[3]))      # σ₀(δW[3])
        c4 = hw(W1f[2] ^ W2f[2])                      # δW[2] (=0)
        total = hw(W1f[18] ^ W2f[18])

        dw18_components['sig1_w16'].append(c1)
        dw18_components['w11'].append(c2)
        dw18_components['sig0_w3'].append(c3)
        dw18_components['w2'].append(c4)
        dw18_components['total'].append(total)

    print(f"  КОМПОНЕНТЫ δW[18]:")
    print(f"    σ₁(δW[16]): {np.mean(dw18_components['sig1_w16']):.1f} бит")
    print(f"    δW[11]:     {np.mean(dw18_components['w11']):.1f} бит  ← a-repair word")
    print(f"    σ₀(δW[3]):  {np.mean(dw18_components['sig0_w3']):.1f} бит  ← break expansion")
    print(f"    δW[2]:      {np.mean(dw18_components['w2']):.1f} бит  (=0 by design)")
    print(f"    TOTAL δW[18]: {np.mean(dw18_components['total']):.1f} бит")
    print(f"    P(δW[18]=0): {sum(1 for x in dw18_components['total'] if x == 0)/len(dw18_components['total'])*100:.2f}%")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("3. ГЛАВНЫЙ ВИНОВНИК")
    print("=" * 70)

    # Какой компонент вносит больше всего?
    components = [
        ('σ₁(δW[16])', dw18_components['sig1_w16']),
        ('δW[11]', dw18_components['w11']),
        ('σ₀(δW[3])', dw18_components['sig0_w3']),
    ]

    components.sort(key=lambda x: np.mean(x[1]), reverse=True)

    print(f"\n  Ранжирование источников δW[18]:")
    for name, vals in components:
        m = np.mean(vals)
        p0 = sum(1 for v in vals if v == 0) / len(vals) * 100
        bar = "█" * int(m * 2)
        print(f"    {name:15s}: {m:5.1f} бит (P(=0)={p0:5.1f}%)  {bar}")

    worst = components[0][0]
    print(f"\n  ГЛАВНЫЙ ВИНОВНИК: {worst}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("4. МОЖНО ЛИ УСТРАНИТЬ ВИНОВНИКА?")
    print("=" * 70)

    # Для каждого компонента: можем ли мы его обнулить?
    print(f"""
  σ₁(δW[16]):
    δW[16] = σ₁(δW[14]) + δW[9] + σ₀(δW[1]) + δW[0]
    δW[0]=0, δW[1]=0 (a-repair design)
    δW[9] и δW[14] — из a-repair (подобраны для δa=0)
    → σ₁(δW[16]) = σ₁(σ₁(δW[14]) + δW[9])
    Можно обнулить если δW[14]=δW[9]=0
    НО: δW[9] и δW[14] — a-repair words, нужны для δa=0 на r=10,15

  δW[11]:
    a-repair word для r=11 (δa[12]=0)
    НУЖЕН для reboot — без него нет схождения
    Обнулить = потерять a-repair на r=12 = потерять reboot

  σ₀(δW[3]):
    δW[3] = break word (1 бит)
    σ₀ разворачивает: ROTR7 ⊕ ROTR18 ⊕ SHR3
    1 бит → 3 бита после σ₀
    Можно уменьшить: выбрать ПОЗИЦИЮ бита δW[3] чтобы σ₀ был минимальным!
""")

    # Тест: какая позиция бита break (δW[3]) минимизирует σ₀(δW[3])?
    print(f"  σ₀ expansion по позиции break-бита:")
    for bit in range(32):
        val = 1 << bit
        expanded = sigma0(val)
        h = hw(expanded)
        bar = "█" * h
        opt = " ★" if h <= 2 else ""
        if h <= 3 or bit < 5 or bit > 28:
            print(f"    bit {bit:2d}: σ₀(2^{bit}) = {expanded:#010x}, HW={h}  {bar}{opt}")


if __name__ == "__main__":
    main()
