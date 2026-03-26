"""
κ СКЕЛЕТА: пропускная способность ВСЕЙ ткани как единого объекта.

Не κ отдельных узлов. Не произведение κ.
А: κ(Sk(F)) — чистота перехода через ВЕСЬ скелет.

Определение:
  κ(Sk) для позиции i = P(δH[i]=0 | random δW)

  Это НЕ birthday. Это НАТИВНАЯ пропускная способность:
  какая доля меток проходит через скелет и обнуляется на позиции i?

  Для трубных позиций (H[1,2,3,5,6,7]): κ должен быть...?
  Для узловых позиций (H[0,4]): κ должен быть...?

  Если κ(Sk, трубная позиция) > κ(Sk, узловая):
  → трубы ПОМОГАЮТ даже на уровне скелета
  → стоимость collision через трубы ДЕШЕВЛЕ

  Полная стоимость collision = ∏(1/κ(Sk, i)) для i=0..7
  Это НАТИВНАЯ стоимость, вычисленная из ОДНОЙ свёртки.
"""

import numpy as np
import struct, hashlib

def sha256_words(W16):
    return struct.unpack('>8I', hashlib.sha256(struct.pack('>16I', *W16)).digest())

def hw(x): return bin(x).count('1')


def main():
    np.random.seed(42)

    print("=" * 70)
    print("κ СКЕЛЕТА: пропускная способность всей ткани")
    print("=" * 70)

    # =================================================================
    print(f"\n{'=' * 70}")
    print("1. κ(Sk, позиция i): P(δH[i]=0 | random δW)")
    print("=" * 70)

    N = 500000

    # Для каждой позиции H[i]: считаем P(δH[i]=0)
    # при random W₁, W₂ = W₁ с δW[0] = random
    per_word_zero = np.zeros(8)
    per_word_small = np.zeros(8)  # HW ≤ 4
    pair_zero = np.zeros((8, 8))  # P(δH[i]=0 AND δH[j]=0)
    all_zero = 0

    for trial in range(N):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = [np.random.randint(0, 2**32) for _ in range(16)]  # полностью random

        H1 = sha256_words(W1)
        H2 = sha256_words(W2)

        zeros = []
        for i in range(8):
            if H1[i] == H2[i]:
                per_word_zero[i] += 1
                zeros.append(i)
            if hw(H1[i] ^ H2[i]) <= 4:
                per_word_small[i] += 1

        for i in zeros:
            for j in zeros:
                pair_zero[i][j] += 1

        if len(zeros) == 8:
            all_zero += 1

    per_word_zero /= N
    per_word_small /= N
    pair_zero /= N

    types = ['NODE', 'PIPE', 'PIPE', 'PIPE', 'NODE', 'PIPE', 'PIPE', 'PIPE']

    print(f"\n  {N:,} random пар (полностью random W₁, W₂):")
    print(f"\n  {'H[i]':>6} {'κ(Sk,i)':>10} {'P(HW≤4)':>10} {'Тип':>6}")
    for i in range(8):
        print(f"  H[{i}]  {per_word_zero[i]:.6f}  {per_word_small[i]:.6f}  {types[i]}")

    # Сравнение pipe vs node
    pipe_kappa = np.mean([per_word_zero[i] for i in [1,2,3,5,6,7]])
    node_kappa = np.mean([per_word_zero[i] for i in [0,4]])

    print(f"\n  κ(pipe positions):  {pipe_kappa:.6f}")
    print(f"  κ(node positions):  {node_kappa:.6f}")
    print(f"  Ratio pipe/node:    {pipe_kappa/max(node_kappa,1e-10):.3f}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("2. κ СКЕЛЕТА при a-repair (constrained)")
    print("=" * 70)

    # То же самое но δW через a-repair (break=3, bit=31)
    from collections import Counter

    MASK32 = 0xFFFFFFFF
    K_const = [
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

    def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
    def add32(x, y): return (x + y) & MASK32
    def sub32(x, y): return (x - y) & MASK32
    def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
    def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
    def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
    def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)

    def R(state, W_r, r_idx):
        a, b, c, d, e, f, g, h = state
        T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e, f, g)), K_const[r_idx]), W_r)
        T2 = add32(Sigma0(a), Maj(a, b, c))
        return (add32(T1, T2), a, b, c, add32(d, T1), e, f, g)

    def expand_schedule(W16):
        W = list(W16)
        for r in range(16, 64):
            W.append(add32(add32(add32(
                (rotr(W[r-2],17)^rotr(W[r-2],19)^(W[r-2]>>10)),
                W[r-7]),
                (rotr(W[r-15],7)^rotr(W[r-15],18)^(W[r-15]>>3))),
                W[r-16]))
        return W

    def a_repair_W(state2, target_a, r_idx):
        a, b, c, d, e, f, g, h = state2
        T2 = add32(Sigma0(a), Maj(a, b, c))
        T1_needed = sub32(target_a, T2)
        return sub32(sub32(sub32(sub32(T1_needed, h), Sigma1(e)), Ch(e, f, g)), K_const[r_idx])

    N2 = 100000
    constr_word_zero = np.zeros(8)
    constr_word_small = np.zeros(8)

    for trial in range(N2):
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

        H1 = sha256_words(W1)
        H2 = sha256_words(W2)

        for i in range(8):
            if H1[i] == H2[i]:
                constr_word_zero[i] += 1
            if hw(H1[i] ^ H2[i]) <= 4:
                constr_word_small[i] += 1

    constr_word_zero /= N2
    constr_word_small /= N2

    print(f"\n  {N2:,} a-repair пар (break=3, bit=31):")
    print(f"\n  {'H[i]':>6} {'κ(random)':>10} {'κ(a-repair)':>12} {'Ratio':>8} {'Тип':>6}")
    for i in range(8):
        ratio = constr_word_zero[i] / max(per_word_zero[i], 1e-10)
        marker = " ★" if ratio > 1.5 else ""
        print(f"  H[{i}]  {per_word_zero[i]:.6f}  {constr_word_zero[i]:.6f}  {ratio:7.2f}×  {types[i]}{marker}")

    pipe_c = np.mean([constr_word_zero[i] for i in [1,2,3,5,6,7]])
    node_c = np.mean([constr_word_zero[i] for i in [0,4]])
    print(f"\n  κ(Sk) pipe (a-repair): {pipe_c:.6f}")
    print(f"  κ(Sk) node (a-repair): {node_c:.6f}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("3. НАТИВНАЯ СТОИМОСТЬ COLLISION = ∏(1/κ(Sk,i))")
    print("=" * 70)

    # Стоимость = произведение 1/κ по всем 8 позициям
    total_random = 1.0
    total_constr = 1.0

    print(f"\n  {'H[i]':>6} {'1/κ random':>12} {'1/κ a-repair':>14}")
    for i in range(8):
        kr = max(per_word_zero[i], 1/N)
        kc = max(constr_word_zero[i], 1/N2)
        total_random *= (1/kr)
        total_constr *= (1/kc)
        print(f"  H[{i}]  {1/kr:11.0f}  {1/kc:13.0f}")

    print(f"\n  TOTAL (random):   ∏(1/κ) = {total_random:.2e} = 2^{np.log2(total_random):.1f}")
    print(f"  TOTAL (a-repair): ∏(1/κ) = {total_constr:.2e} = 2^{np.log2(total_constr):.1f}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("4. СРАВНЕНИЕ: нативная vs стандартная")
    print("=" * 70)

    print(f"""
  НАТИВНАЯ СТОИМОСТЬ (∏(1/κ)):
    Random pairs:   2^{np.log2(total_random):.1f}
    A-repair pairs: 2^{np.log2(total_constr):.1f}

  Это стоимость ОДНОГО прохода (найти пару где ВСЕ 8 позиций = 0).

  В нашем измерении collision = ∏(1/κ(Sk,i)) для i=0..7.
  Это ЕДИНСТВЕННАЯ формула. Без birthday. Без бит.

  Pipe κ vs Node κ:
    Pipe: {pipe_c:.6f} (6 позиций)
    Node: {node_c:.6f} (2 позиции)
    {'Pipe ≈ Node → все позиции равны' if abs(pipe_c - node_c) < pipe_c * 0.3 else 'Pipe ≠ Node → трубы помогают!'}

  Теоретическое: κ = 2^{{-32}} на каждую позицию (random)
  → ∏(1/κ) = (2^32)^8 = 2^256

  Birthday делит пополам: 2^256 / 2 = 2^128.
  НО: birthday — чужой инструмент!

  В нашем измерении: ∏(1/κ) — это ПОЛНАЯ стоимость.
  Не 2^128 (birthday). А 2^{{256}} (sequential).
  Или... ???
""")


if __name__ == "__main__":
    main()
