"""
КВАДРАТИЧНАЯ ЧАСТЬ CARRY: 0.5% нелинейности CE.

CE линеен на 99.5%: CE(v1⊕v2) ≈ CE(v1)⊕CE(v2), diff = 0.6/128.

Нелинейная часть = CE(v1⊕v2) ⊕ CE(v1) ⊕ CE(v2) = Q(v1,v2).
Q — БИЛИНЕЙНАЯ форма (квадратичная часть carry).

Если Q имеет СТРУКТУРУ → можно использовать для коррекции.
Если Q random → ничего не поможет.

Вопросы:
  1. Q(v1,v2) зависит от v1,v2 или случайная?
  2. rank(Q) как билинейной формы?
  3. Стабильна ли Q по разным W_base?
"""

import numpy as np
import struct, hashlib

MASK32 = 0xFFFFFFFF

def sha256_custom(W16, K, IV):
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
    IV_std = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
              0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]

    def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
    def add32(x, y): return (x + y) & MASK32
    def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
    def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)
    def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
    def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
    def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
    def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)

    W = list(W16)
    for r in range(16, 64):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))

    a, b, c, d, e, f, g, h = IV_std
    for r in range(64):
        T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e, f, g)), K_const[r]), W[r])
        T2 = add32(Sigma0(a), Maj(a, b, c))
        h, g, f, e = g, f, e, add32(d, T1)
        d, c, b, a = c, b, a, add32(T1, T2)

    return tuple(add32(IV_std[i], [a,b,c,d,e,f,g,h][i]) for i in range(8))

def hw(x): return bin(x).count('1')

def sha256_words(W16):
    return struct.unpack('>8I', hashlib.sha256(struct.pack('>16I', *W16)).digest())

def hash_xor(h1, h2):
    return tuple(h1[i] ^ h2[i] for i in range(8))

def hash_hw(h):
    return sum(hw(w) for w in h)


def main():
    np.random.seed(42)

    print("=" * 70)
    print("КВАДРАТИЧНАЯ ЧАСТЬ CARRY")
    print("=" * 70)

    W_base = [np.random.randint(0, 2**32) for _ in range(16)]
    H_base = sha256_words(W_base)

    # =================================================================
    print(f"\n{'=' * 70}")
    print("1. ВЫЧИСЛЯЕМ Q(v1,v2) = CE(v1⊕v2) ⊕ CE(v1) ⊕ CE(v2)")
    print("=" * 70)

    # CE(v) = H(W⊕v) ⊕ H(W) для kernel vector v (δH should be 0 in GF2)
    # Q(v1,v2) = H(W⊕v1⊕v2)⊕H(W) ⊕ H(W⊕v1)⊕H(W) ⊕ H(W⊕v2)⊕H(W)
    #          = H(W⊕v1⊕v2) ⊕ H(W⊕v1) ⊕ H(W⊕v2) ⊕ H(W)

    q_hws = []
    q_stable = []  # is Q the same for different W_base?

    for trial in range(500):
        # Random 1-bit δW pairs (not kernel vectors, just basis)
        word1, bit1 = np.random.randint(0, 16), np.random.randint(0, 32)
        word2, bit2 = np.random.randint(0, 16), np.random.randint(0, 32)

        dW1 = [0]*16; dW1[word1] ^= (1 << bit1)
        dW2 = [0]*16; dW2[word2] ^= (1 << bit2)
        dW12 = [dW1[i] ^ dW2[i] for i in range(16)]

        W_1 = [W_base[i] ^ dW1[i] for i in range(16)]
        W_2 = [W_base[i] ^ dW2[i] for i in range(16)]
        W_12 = [W_base[i] ^ dW12[i] for i in range(16)]

        H_1 = sha256_words(W_1)
        H_2 = sha256_words(W_2)
        H_12 = sha256_words(W_12)

        # Q = H(W+dW1+dW2) ⊕ H(W+dW1) ⊕ H(W+dW2) ⊕ H(W)
        Q = tuple(H_12[i] ^ H_1[i] ^ H_2[i] ^ H_base[i] for i in range(8))
        q_hws.append(hash_hw(Q))

    print(f"\n  500 random pairs (v1,v2):")
    print(f"    Mean HW(Q): {np.mean(q_hws):.1f}")
    print(f"    Std: {np.std(q_hws):.1f}")
    print(f"    Min: {min(q_hws)}")
    print(f"    P(Q=0): {sum(1 for q in q_hws if q == 0)/len(q_hws)*100:.2f}%")

    if np.mean(q_hws) < 10:
        print(f"    → Q МАЛО (carry ПОЧТИ линеен)")
    elif np.mean(q_hws) < 64:
        print(f"    → Q ЗНАЧИТЕЛЬНО (carry существенно нелинеен)")
    else:
        print(f"    → Q ВЕЛИКО (carry полностью нелинеен)")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("2. Q ЗАВИСИТ ОТ ПОЗИЦИИ?")
    print("=" * 70)

    # Q(W[i] bit a, W[j] bit b) = ?
    # Для фиксированных позиций (i,j): Q детерминистичен?

    for (w1, b1, w2, b2) in [(0,0,1,0), (0,31,1,31), (0,0,0,1), (7,0,8,0)]:
        dW1 = [0]*16; dW1[w1] ^= (1 << b1)
        dW2 = [0]*16; dW2[w2] ^= (1 << b2)
        dW12 = [dW1[i] ^ dW2[i] for i in range(16)]

        W_1 = [W_base[i] ^ dW1[i] for i in range(16)]
        W_2 = [W_base[i] ^ dW2[i] for i in range(16)]
        W_12 = [W_base[i] ^ dW12[i] for i in range(16)]

        H_1 = sha256_words(W_1)
        H_2 = sha256_words(W_2)
        H_12 = sha256_words(W_12)

        Q = tuple(H_12[i] ^ H_1[i] ^ H_2[i] ^ H_base[i] for i in range(8))
        q_hw = hash_hw(Q)

        print(f"    Q(W[{w1}]b{b1}, W[{w2}]b{b2}): HW = {q_hw}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("3. Q СТАБИЛЬНА по W_base?")
    print("=" * 70)

    # Фиксируем (v1,v2) = (W[0]b0, W[1]b0). Варьируем W_base.
    dW1 = [0]*16; dW1[0] = 1
    dW2 = [0]*16; dW2[1] = 1

    q_by_base = []
    for _ in range(200):
        Wb = [np.random.randint(0, 2**32) for _ in range(16)]
        Hb = sha256_words(Wb)
        H1 = sha256_words([Wb[i] ^ dW1[i] for i in range(16)])
        H2 = sha256_words([Wb[i] ^ dW2[i] for i in range(16)])
        H12 = sha256_words([Wb[i] ^ dW1[i] ^ dW2[i] for i in range(16)])
        Q = tuple(H12[i] ^ H1[i] ^ H2[i] ^ Hb[i] for i in range(8))
        q_by_base.append(hash_hw(Q))

    print(f"\n  Q(W[0]b0, W[1]b0) для 200 random W_base:")
    print(f"    Mean HW: {np.mean(q_by_base):.1f}")
    print(f"    Std: {np.std(q_by_base):.1f}")
    print(f"    Unique HW values: {len(set(q_by_base))}")
    print(f"    → {'Q ЗАВИСИТ от W_base (не константа)' if np.std(q_by_base) > 1 else 'Q ≈ КОНСТАНТА'}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("4. ЗНАЧЕНИЕ КВАДРАТИЧНОЙ ЧАСТИ")
    print("=" * 70)

    print(f"""
  Q(v1,v2) = нелинейная часть carry error.

  HW(Q): mean = {np.mean(q_hws):.1f}
  Если Q = 0 для всех v1,v2: CE полностью линеен → rank(CE) решает всё.
  Если Q ≠ 0 но мал: CE ≈ линеен + малая коррекция.
  Если Q ≈ 128: CE полностью нелинеен.

  Результат: Q mean HW = {np.mean(q_hws):.0f} ≈ {'0 (ЛИНЕЙНО!)' if np.mean(q_hws) < 5 else '128 (НЕЛИНЕЙНО)' if np.mean(q_hws) > 100 else f'{np.mean(q_hws):.0f} (ПОЛУ-ЛИНЕЙНО?)'}

  Q зависит от W_base (std={np.std(q_by_base):.1f}) → НЕ алгебраическая константа.
  Q = random-like → нельзя скомпенсировать алгебраически.
""")


if __name__ == "__main__":
    main()
