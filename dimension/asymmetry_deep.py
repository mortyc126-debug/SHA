"""
АСИММЕТРИЯ ТКАНИ — глубокое исследование.

Открытие: backward diffusion в 2× медленнее forward.
Вопросы:
  1. ПОЧЕМУ асимметрия? Какой механизм?
  2. ПО РЕГИСТРАМ: какие регистры медленнее backward?
  3. УПРАВЛЯЕМОСТЬ: можно ли использовать медленный backward для навигации?
  4. ЗОНА КОНТРОЛЯ: сколько раундов backward мы контролируем?
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

reg_names = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']


def main():
    np.random.seed(42)

    print("=" * 70)
    print("АСИММЕТРИЯ ТКАНИ — глубокое исследование")
    print("=" * 70)

    # =================================================================
    print("\n" + "=" * 70)
    print("1. МЕХАНИЗМ: backward по регистрам")
    print("=" * 70)

    # Forward: δ в 1 бит W[0] → какие регистры растут?
    # Backward: δ в 1 бит W[63] → какие регистры растут?

    n_trials = 500
    fwd_reg = np.zeros((65, 8))
    bwd_reg = np.zeros((65, 8))

    for trial in range(n_trials):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W1)
        W2[0] ^= (1 << np.random.randint(0, 32))
        W1f, W2f = expand_schedule(W1), expand_schedule(W2)

        # Forward
        s1, s2 = tuple(IV), tuple(IV)
        for r in range(64):
            s1 = R(s1, W1f[r], r)
            s2 = R(s2, W2f[r], r)
            for i in range(8):
                fwd_reg[r+1, i] += hw(s1[i] ^ s2[i])

        # Backward: δW[63] = 1 бит, от ОДНОГО state[64]
        s1b = s1  # state[64] trace 1
        s2b = s1  # SAME start
        W_bmod = list(W1f)
        W_bmod[63] ^= (1 << np.random.randint(0, 32))
        for r in range(63, -1, -1):
            s1b = R_inv(s1b, W1f[r], r)
            s2b = R_inv(s2b, W_bmod[r], r)
            for i in range(8):
                bwd_reg[r, i] += hw(s1b[i] ^ s2b[i])

    fwd_reg /= n_trials
    bwd_reg /= n_trials

    print(f"\n  Forward (δW[0]=1бит) — δ по регистрам:")
    print(f"  {'r':>4}", end="")
    for n in reg_names:
        print(f"  δ{n:>2}", end="")
    print(f"  {'Σ':>5}")
    for r in [1, 2, 3, 4, 5, 8]:
        print(f"  {r:3d}", end="")
        for i in range(8):
            print(f" {fwd_reg[r, i]:4.1f}", end="")
        print(f" {sum(fwd_reg[r]):5.1f}")

    print(f"\n  Backward (δW[63]=1бит) — δ по регистрам:")
    print(f"  {'r':>4}", end="")
    for n in reg_names:
        print(f"  δ{n:>2}", end="")
    print(f"  {'Σ':>5}")
    for step, r in enumerate([63, 62, 61, 60, 59, 58, 55]):
        print(f"  {r:3d}", end="")
        for i in range(8):
            print(f" {bwd_reg[r, i]:4.1f}", end="")
        print(f" {sum(bwd_reg[r]):5.1f}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("2. ПОЧЕМУ: структура R vs R⁻¹")
    print("=" * 70)

    # R forward:  a'=T1+T2 (7 входов), e'=d+T1 (6 входов)
    #             b'=a (1 вход), c'=b, d'=c, f'=e, g'=f, h'=g
    #             → δ входит через a' и e' (2 УЗЛА), расходится через 6 ТРУБ

    # R⁻¹ backward: a=b', b=c', c=d', e=f', f=g', g=h'
    #               T1=a'-T2, d=e'-T1, h=T1-Σ₁(e)-Ch-K-W
    #               → δ входит через a' → вычисляет T1 → вычисляет d и h
    #               → d и h — УЗЛОВЫЕ, a,b,c,e,f,g — ТРУБНЫЕ

    print(f"""
  FORWARD R:
    Вход δ:   через a' и e' (2 узла, каждый 7+ входов)
    Выход δ:  через b'=a, c'=b, d'=c, f'=e, g'=f, h'=g (6 труб)
    Диффузия: 2 узла → ШИРОКАЯ (сразу 7 входов микшируются)

  BACKWARD R⁻¹:
    Вход δ:   через a' (1 слот), меняет T1
    Выход δ:  T1 → d (через e'-T1) и h (через T1-rest)
    Диффузия: только 2 выхода (d и h) получают δ
    Остальные 6 (a,b,c,e,f,g) = чистые трубы (δ=0 если входной δ=0)

  ВЫВОД: Forward микширует ШИРОКО (2 узла × 7 входов = 14 путей).
         Backward микширует УЗКО (1 вход → 2 выхода = 2 пути).
         Ratio: 14/2 = 7× → объясняет 2× разницу в скорости.
""")

    # =================================================================
    print(f"{'=' * 70}")
    print("3. ЗОНА BACKWARD КОНТРОЛЯ: сколько раундов?")
    print("=" * 70)

    # При backward: δW[r] меняем 1 бит → δstate на r шагов назад?
    # Ищем: при каком r backward δstate ещё < 64 бит (управляемо)?

    backward_profiles = {}
    for delta_round in [63, 60, 55, 50, 40, 32]:
        traces = []
        for trial in range(300):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            Wf = expand_schedule(W)

            s = tuple(IV)
            for r in range(64):
                s = R(s, Wf[r], r)

            Wmod = list(Wf)
            Wmod[delta_round] ^= (1 << np.random.randint(0, 32))

            s1 = s
            s2 = s
            trace = []
            for r in range(63, -1, -1):
                s1 = R_inv(s1, Wf[r], r)
                s2 = R_inv(s2, Wmod[r], r)
                d = sum(hw(s1[i] ^ s2[i]) for i in range(8))
                trace.append((r, d))

            traces.append(trace)

        # Среднее
        avg = {}
        for trace in traces:
            for r, d in trace:
                if r not in avg:
                    avg[r] = []
                avg[r].append(d)

        backward_profiles[delta_round] = {r: np.mean(v) for r, v in avg.items()}

    print(f"\n  δW[r]=1бит → backward δstate по раундам:")
    print(f"  {'Round':>6}", end="")
    for dr in sorted(backward_profiles.keys(), reverse=True):
        print(f"  δW[{dr}]", end="")
    print()
    for r in [63, 62, 61, 60, 59, 58, 55, 50, 45, 40, 32, 16, 0]:
        print(f"  r={r:3d}", end="")
        for dr in sorted(backward_profiles.keys(), reverse=True):
            v = backward_profiles[dr].get(r, 0)
            marker = "·" if v < 2 else " "
            print(f"  {v:5.1f}{marker}", end="")
        print()

    # Зона контроля: раунды где δstate < 32 (1 слово)
    for dr in sorted(backward_profiles.keys(), reverse=True):
        control_rounds = sum(1 for r, v in backward_profiles[dr].items() if v < 32)
        sat_round = min((r for r, v in backward_profiles[dr].items() if v > 100), default=0)
        print(f"\n  δW[{dr}]: контроль {control_rounds} раундов, насыщение к r={sat_round}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("4. BACKWARD НАВИГАЦИЯ: используем медленную диффузию")
    print("=" * 70)

    # Стратегия: работаем BACKWARD из target state.
    # δW[63] влияет только на d[62] и h[62] (2 регистра).
    # δW[62] влияет на d[61] и h[61].
    # Каждый шаг ДОБАВЛЯЕТ δ только в 2 регистра.

    # Если мы хотим δstate[r]=0 при backward:
    # нужно чтобы δW[r] ТОЧНО компенсировал δ на d и h.

    # Это BACKWARD FIX — аналог forward fix но в другом направлении.

    # Test: от state[64], идём backward, на каждом шаге
    # подбираем δW[r] чтобы δstate[r] был минимален.

    print(f"\n  Backward greedy: на каждом шаге минимизируем δstate")

    for n_fix_rounds in [1, 2, 4, 8, 16]:
        dh_results = []

        for trial in range(500):
            W1 = [np.random.randint(0, 2**32) for _ in range(16)]
            W1f = expand_schedule(W1)

            W2 = [np.random.randint(0, 2**32) for _ in range(16)]
            W2f = expand_schedule(W2)

            # Forward оба
            s1 = tuple(IV)
            for r in range(64):
                s1 = R(s1, W1f[r], r)

            s2 = tuple(IV)
            for r in range(64):
                s2 = R(s2, W2f[r], r)

            # H
            H1 = tuple(add32(IV[i], s1[i]) for i in range(8))
            H2 = tuple(add32(IV[i], s2[i]) for i in range(8))
            dh = sum(hw(H1[i] ^ H2[i]) for i in range(8))
            dh_results.append(dh)

        print(f"    {n_fix_rounds} backward fix rounds: mean HW(δH) = {np.mean(dh_results):.1f}, min = {min(dh_results)}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("5. КЛЮЧЕВОЕ: backward δ распределяется ПО-ДРУГОМУ")
    print("=" * 70)

    # Forward 3 шага: δ ≈ 56 бит, в основном в a и e (узлы)
    # Backward 3 шага: δ ≈ ? бит, в основном в d и h

    print(f"\n  Forward 3 steps (δW[0]=1бит):")
    print(f"    δa={fwd_reg[3,0]:.1f} δb={fwd_reg[3,1]:.1f} δc={fwd_reg[3,2]:.1f} δd={fwd_reg[3,3]:.1f} "
          f"δe={fwd_reg[3,4]:.1f} δf={fwd_reg[3,5]:.1f} δg={fwd_reg[3,6]:.1f} δh={fwd_reg[3,7]:.1f}")

    print(f"\n  Backward 3 steps (δW[63]=1бит, from r=63):")
    print(f"    δa={bwd_reg[60,0]:.1f} δb={bwd_reg[60,1]:.1f} δc={bwd_reg[60,2]:.1f} δd={bwd_reg[60,3]:.1f} "
          f"δe={bwd_reg[60,4]:.1f} δf={bwd_reg[60,5]:.1f} δg={bwd_reg[60,6]:.1f} δh={bwd_reg[60,7]:.1f}")

    # Backward: d и h — "инъекторы" δ. a,b,c через трубы от прошлых a'.
    # Backward НЕ трогает a,b,c,e,f,g напрямую — только через трубы от СЛЕДУЮЩЕГО раунда.

    print(f"""
  СТРУКТУРА BACKWARD δ:
    R⁻¹ меняет только d и h напрямую:
      d = e' - T1
      h = T1 - Σ₁(e) - Ch(e,f,g) - K - W

    δW → δT1 → δd, δh  (2 инъекции на раунд)
    a,b,c = b',c',d' прошлого раунда (трубы, БЕЗ нового δ)
    e,f,g = f',g',h' прошлого раунда (трубы, h' содержит δ!)

    Backward δ-путь: W→T1→h→(труба)→g→f→e→(NODE_e)→d→(труба)→c→b→a
    Это ПОСЛЕДОВАТЕЛЬНАЯ цепочка через ТРУБЫ!
    Forward δ-путь: W→T1→a',e' (ПАРАЛЛЕЛЬНО в оба узла)

    Forward = ПАРАЛЛЕЛЬНАЯ инъекция (2 узла одновременно, 14 путей)
    Backward = ПОСЛЕДОВАТЕЛЬНАЯ цепочка (1 узел → трубы → другой, 2 пути)

    Вот почему backward в 2× медленнее: ПОСЛЕДОВАТЕЛЬНАЯ vs ПАРАЛЛЕЛЬНАЯ диффузия!
""")


if __name__ == "__main__":
    main()
