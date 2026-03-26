"""
ИССЛЕДОВАНИЕ НАШЕГО МИРА — с нуля, без импорта старых законов.

Что мы НЕ знаем о нашем измерении:
  1. Как ведут себя БОЛЬШИЕ метки (δW с HW > 100)?
  2. Как выглядит state-пространство изнутри (через R⁻¹)?
  3. Есть ли "потоки" в ткани — направления наименьшего сопротивления?
  4. Симметрична ли ткань — или есть привилегированные области?

Подход: не тестируем гипотезы, а НАБЛЮДАЕМ. Как натуралисты в новом мире.
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

def sha256_words(W16):
    return struct.unpack('>8I', hashlib.sha256(struct.pack('>16I', *W16)).digest())

def state_hw(s1, s2):
    return sum(hw(s1[i] ^ s2[i]) for i in range(8))


def main():
    np.random.seed(42)

    print("=" * 70)
    print("ИССЛЕДОВАНИЕ НАШЕГО МИРА — наблюдение, не гипотезы")
    print("=" * 70)

    # =================================================================
    print("\n" + "=" * 70)
    print("НАБЛЮДЕНИЕ 1: Ландшафт при РАЗНЫХ размерах δW")
    print("  (от 1 бита до 256 бит)")
    print("=" * 70)

    for target_hw_dw in [1, 2, 4, 8, 16, 32, 64, 128, 192, 256]:
        dh_results = []
        for _ in range(3000):
            W1 = [np.random.randint(0, 2**32) for _ in range(16)]
            W2 = list(W1)

            # Генерируем δW с заданным HW
            bits_to_flip = np.random.choice(512, target_hw_dw, replace=False)
            for bit in bits_to_flip:
                word = bit // 32
                pos = bit % 32
                W2[word] ^= (1 << pos)

            H1 = sha256_words(W1)
            H2 = sha256_words(W2)
            dh = sum(hw(H1[i] ^ H2[i]) for i in range(8))
            dh_results.append(dh)

        m = np.mean(dh_results)
        s = np.std(dh_results)
        mn = min(dh_results)
        bar = "█" * int(m / 4)
        print(f"  HW(δW)={target_hw_dw:3d}: HW(δH) = {m:6.1f} ± {s:.1f}  min={mn:3d}  {bar}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("НАБЛЮДЕНИЕ 2: R⁻¹ navigation — движение из state-пространства")
    print("  Стартуем из target δH=0, идём НАЗАД через R⁻¹")
    print("=" * 70)

    # Идея: выбираем два РАЗНЫХ W, но целевой state[64] ОДИНАКОВЫЙ.
    # Идём назад R⁻¹ от одного state до state[0].
    # Смотрим: как расходятся два backward-trace?

    W1 = [np.random.randint(0, 2**32) for _ in range(16)]
    W1_full = list(W1)
    for r in range(16, 64):
        W1_full.append(add32(add32(add32(sigma1(W1_full[r-2]), W1_full[r-7]),
                       sigma0(W1_full[r-15])), W1_full[r-16]))

    # Forward: compute state[64]
    s = tuple(IV)
    for r in range(64):
        s = R(s, W1_full[r], r)
    state_64_target = s

    # Backward trace 1: с правильными W
    s_back = state_64_target
    for r in range(63, -1, -1):
        s_back = R_inv(s_back, W1_full[r], r)
    recovered_iv = s_back

    print(f"\n  Forward → state[64] → Backward → IV:")
    print(f"  Recovered IV == real IV: {recovered_iv == tuple(IV)}")

    # Backward trace 2: с ДРУГИМИ W (W2), от ТОГО ЖЕ state[64]
    W2 = [np.random.randint(0, 2**32) for _ in range(16)]
    W2_full = list(W2)
    for r in range(16, 64):
        W2_full.append(add32(add32(add32(sigma1(W2_full[r-2]), W2_full[r-7]),
                       sigma0(W2_full[r-15])), W2_full[r-16]))

    s_back2 = state_64_target
    backward_states = [s_back2]
    for r in range(63, -1, -1):
        s_back2 = R_inv(s_back2, W2_full[r], r)
        backward_states.append(s_back2)

    backward_states.reverse()  # теперь [state_0, ..., state_64]

    # Backward-IV — это state[0] при "неправильных" W
    fake_iv = backward_states[0]

    print(f"\n  Backward с W2 (чужие W) от того же state[64]:")
    print(f"  'Fake IV': ({', '.join(hex(x)[:10] for x in fake_iv)})")
    print(f"  Real IV:   ({', '.join(hex(x)[:10] for x in IV)})")
    print(f"  δIV HW:    {state_hw(fake_iv, tuple(IV))}")

    # Вопрос: если бы IV был СВОБОДЕН, это была бы collision!
    # state[64](W1, real_IV) = state[64](W2, fake_IV) = target
    # → H(W1) - IV = H(W2) - fake_IV
    # → Для collision: нужен fake_IV == real_IV

    # =================================================================
    print(f"\n{'=' * 70}")
    print("НАБЛЮДЕНИЕ 3: Backward divergence — как быстро расходятся?")
    print("=" * 70)

    # Два backward trace от одного state[64], с δW[63] = 1 бит
    results_by_round = []

    for trial in range(200):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        W_full = list(W)
        for r in range(16, 64):
            W_full.append(add32(add32(add32(sigma1(W_full[r-2]), W_full[r-7]),
                           sigma0(W_full[r-15])), W_full[r-16]))

        s = tuple(IV)
        for r in range(64):
            s = R(s, W_full[r], r)

        # Backward 1: with W_full
        s1 = s
        # Backward 2: with W_full but δW[63] = 1 bit
        W_mod = list(W_full)
        W_mod[63] ^= (1 << np.random.randint(0, 32))
        s2 = s

        trace = []
        for r in range(63, -1, -1):
            s1 = R_inv(s1, W_full[r], r)
            s2 = R_inv(s2, W_mod[r], r)
            d = state_hw(s1, s2)
            trace.append(d)

        results_by_round.append(trace)

    mean_trace = np.mean(results_by_round, axis=0)

    print(f"\n  δW[63]=1бит, backward от state[64]:")
    print(f"  {'Back step':>10} {'Round':>6} {'δstate':>8}")
    for step in [0, 1, 2, 3, 4, 5, 8, 16, 32, 48, 63]:
        r = 63 - step
        print(f"  step={step:3d}     r={r:3d}   {mean_trace[step]:7.1f}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("НАБЛЮДЕНИЕ 4: Backward SYMMETRY — forward vs backward diffusion")
    print("=" * 70)

    # Forward: δW[0]=1бит → δstate растёт 3→24→86→128
    # Backward: δW[63]=1бит → δstate растёт ???
    # Одинаковая ли скорость?

    forward_diff = np.zeros(65)
    backward_diff = np.zeros(65)

    for trial in range(500):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W1)
        W2[0] ^= (1 << np.random.randint(0, 32))

        W1f = list(W1)
        W2f = list(W2)
        for r in range(16, 64):
            W1f.append(add32(add32(add32(sigma1(W1f[r-2]), W1f[r-7]), sigma0(W1f[r-15])), W1f[r-16]))
            W2f.append(add32(add32(add32(sigma1(W2f[r-2]), W2f[r-7]), sigma0(W2f[r-15])), W2f[r-16]))

        # Forward
        s1, s2 = tuple(IV), tuple(IV)
        for r in range(64):
            s1 = R(s1, W1f[r], r)
            s2 = R(s2, W2f[r], r)
            forward_diff[r + 1] += state_hw(s1, s2)

        # Backward from state[64] of W1, using W1 vs W2
        s1b = s1  # state[64] of W1
        s2b = s1  # SAME start (like collision target)
        for r in range(63, -1, -1):
            s1b = R_inv(s1b, W1f[r], r)
            s2b = R_inv(s2b, W2f[r], r)
            backward_diff[r] += state_hw(s1b, s2b)

    forward_diff /= 500
    backward_diff /= 500

    print(f"\n  {'Round':>6} {'Forward δ':>10} {'Backward δ':>11} {'Ratio':>8}")
    for r in [0, 1, 2, 3, 4, 5, 8, 16, 32, 48, 60, 62, 63, 64]:
        if r <= 64:
            f = forward_diff[r]
            b = backward_diff[r] if r < 64 else 0
            ratio = f / max(b, 0.1) if b > 0.1 else 0
            print(f"  r={r:3d}   {f:9.1f}  {b:10.1f}  {ratio:7.2f}")

    # Асимметрия?
    fwd_speed = forward_diff[5]  # diffusion at round 5
    bwd_speed = backward_diff[59]  # backward diffusion "5 steps back"
    print(f"\n  Forward 5 steps:  δ = {fwd_speed:.1f}")
    print(f"  Backward 5 steps: δ = {bwd_speed:.1f}")
    print(f"  → {'СИММЕТРИЧНО' if abs(fwd_speed - bwd_speed) < 10 else 'АСИММЕТРИЯ!'}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("НАБЛЮДЕНИЕ 5: Ткань при НУЛЕВОМ δW — только разные IV")
    print("=" * 70)

    # Если IV свободен (free-start): collision значительно проще?
    # δW=0, но δIV≠0. Как быстро расходится?

    for delta_hw_iv in [1, 4, 16, 64, 128]:
        dh_results = []
        for _ in range(2000):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            W_full = list(W)
            for r in range(16, 64):
                W_full.append(add32(add32(add32(sigma1(W_full[r-2]), W_full[r-7]),
                               sigma0(W_full[r-15])), W_full[r-16]))

            IV1 = list(IV)
            IV2 = list(IV)
            bits = np.random.choice(256, delta_hw_iv, replace=False)
            for bit in bits:
                reg = bit // 32
                pos = bit % 32
                IV2[reg] ^= (1 << pos)

            s1 = tuple(IV1)
            s2 = tuple(IV2)
            for r in range(64):
                s1 = R(s1, W_full[r], r)
                s2 = R(s2, W_full[r], r)

            # H = state + IV (РАЗНЫЕ IV!)
            H1 = tuple(add32(IV1[i], s1[i]) for i in range(8))
            H2 = tuple(add32(IV2[i], s2[i]) for i in range(8))
            dh = sum(hw(H1[i] ^ H2[i]) for i in range(8))
            dh_results.append(dh)

        print(f"  HW(δIV)={delta_hw_iv:3d}: HW(δH) = {np.mean(dh_results):6.1f} ± {np.std(dh_results):.1f}  min={min(dh_results):3d}")


if __name__ == "__main__":
    main()
