"""
HIGHER-ORDER ZEROS: реальные или trivial?

При r=32, order-2: 6/2000 full 256-bit zeros.
Random function: ожидаем 0.

Возможные причины:
  A. BIT COLLISION: два случайных бита попали на один и тот же (word,bit)
     → δW = 0 → тривиальный ноль
  B. STRUCTURAL: реальная алгебраическая степень < 2
  C. BUG: ошибка в коде

Проверяем ВСЁ.
"""

import numpy as np
import struct, hashlib

MASK32 = 0xFFFFFFFF

def hw(x): return bin(x).count('1')
def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def add32(x, y): return (x + y) & MASK32
def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)
def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)

K = [
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
]
IV = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]


def sha256_r(W16, n_r):
    W = list(W16)
    for r in range(16, max(n_r, 16)):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))
    a,b,c,d,e,f,g,h = IV
    for r in range(n_r):
        T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K[r]), W[r])
        T2 = add32(Sigma0(a), Maj(a,b,c))
        h,g,f,e = g,f,e,add32(d,T1)
        d,c,b,a = c,b,a,add32(T1,T2)
    return tuple(add32(IV[i], [a,b,c,d,e,f,g,h][i]) for i in range(8))


def main():
    np.random.seed(42)

    print("=" * 70)
    print("HIGHER-ORDER ZEROS: Bug hunt")
    print("=" * 70)

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("CHECK A: Are the two random bits IDENTICAL?")
    print("=" * 70)

    # If bits = [(w1,b1), (w2,b2)] and w1==w2, b1==b2 → XOR cancels → δW=0

    # P(collision) = 1/512 per trial. In 2000 trials: ~4 expected.
    # THIS EXPLAINS THE ~3-6 zeros!

    N = 5000
    collision_count = 0
    zero_with_collision = 0
    zero_without_collision = 0

    for _ in range(N):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        w1, b1 = np.random.randint(0, 16), np.random.randint(0, 32)
        w2, b2 = np.random.randint(0, 16), np.random.randint(0, 32)

        is_collision = (w1 == w2 and b1 == b2)
        if is_collision:
            collision_count += 1

        # Compute order-2 differential
        xor_sum = [0] * 8
        for mask in range(4):
            Wm = list(W)
            if mask & 1: Wm[w1] ^= (1 << b1)
            if mask & 2: Wm[w2] ^= (1 << b2)
            H = sha256_r(Wm, 32)
            for i in range(8):
                xor_sum[i] ^= H[i]

        is_zero = all(x == 0 for x in xor_sum)
        if is_zero:
            if is_collision:
                zero_with_collision += 1
            else:
                zero_without_collision += 1

    total_zeros = zero_with_collision + zero_without_collision
    print(f"  {N} trials at r=32:")
    print(f"    Bit collisions (same (w,b)): {collision_count} (expected: {N/512:.0f})")
    print(f"    Total zeros: {total_zeros}")
    print(f"    Zeros WITH bit collision: {zero_with_collision}")
    print(f"    Zeros WITHOUT bit collision: {zero_without_collision}")
    print(f"    → {'ALL zeros from collisions!' if zero_without_collision == 0 else 'REAL zeros exist!'}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("CHECK B: Force DISTINCT bits, test again")
    print("=" * 70)

    for n_r in [16, 20, 24, 32, 48, 64]:
        zeros = 0
        N_test = 5000

        for _ in range(N_test):
            W = [np.random.randint(0, 2**32) for _ in range(16)]

            # FORCE distinct bits
            while True:
                w1, b1 = np.random.randint(0, 16), np.random.randint(0, 32)
                w2, b2 = np.random.randint(0, 16), np.random.randint(0, 32)
                if w1 != w2 or b1 != b2:
                    break

            xor_sum = [0] * 8
            for mask in range(4):
                Wm = list(W)
                if mask & 1: Wm[w1] ^= (1 << b1)
                if mask & 2: Wm[w2] ^= (1 << b2)
                H = sha256_r(Wm, n_r)
                for i in range(8):
                    xor_sum[i] ^= H[i]

            if all(x == 0 for x in xor_sum):
                zeros += 1

        print(f"  r={n_r:>2}, distinct bits: zeros = {zeros}/{N_test}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("CHECK C: Order-3 with GUARANTEED distinct bits")
    print("=" * 70)

    for n_r in [16, 20, 24, 32, 64]:
        zeros = 0
        N_test = 3000

        for _ in range(N_test):
            W = [np.random.randint(0, 2**32) for _ in range(16)]

            # 3 GUARANTEED distinct bits
            bits_set = set()
            bit_list = []
            while len(bit_list) < 3:
                w, b = np.random.randint(0, 16), np.random.randint(0, 32)
                if (w, b) not in bits_set:
                    bits_set.add((w, b))
                    bit_list.append((w, b))

            xor_sum = [0] * 8
            for mask in range(8):  # 2^3
                Wm = list(W)
                for j in range(3):
                    if (mask >> j) & 1:
                        Wm[bit_list[j][0]] ^= (1 << bit_list[j][1])
                H = sha256_r(Wm, n_r)
                for i in range(8):
                    xor_sum[i] ^= H[i]

            if all(x == 0 for x in xor_sum):
                zeros += 1

        print(f"  r={n_r:>2}, order-3 distinct: zeros = {zeros}/{N_test}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("CHECK D: per-word zeros (not full) with distinct bits")
    print("=" * 70)

    for n_r in [16, 24, 32, 64]:
        word_zeros = [0] * 8
        N_test = 5000

        for _ in range(N_test):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            while True:
                w1, b1 = np.random.randint(0, 16), np.random.randint(0, 32)
                w2, b2 = np.random.randint(0, 16), np.random.randint(0, 32)
                if w1 != w2 or b1 != b2: break

            xor_sum = [0] * 8
            for mask in range(4):
                Wm = list(W)
                if mask & 1: Wm[w1] ^= (1 << b1)
                if mask & 2: Wm[w2] ^= (1 << b2)
                H = sha256_r(Wm, n_r)
                for i in range(8):
                    xor_sum[i] ^= H[i]

            for i in range(8):
                if xor_sum[i] == 0:
                    word_zeros[i] += 1

        rates = [f"{z/N_test:.4f}" for z in word_zeros]
        expected = 1/2**32
        any_significant = any(z/N_test > 0.001 for z in word_zeros)
        print(f"  r={n_r:>2}: per-word zero rates = {rates}")
        print(f"       expected (random): {expected:.10f}")
        if any_significant:
            print(f"       ★ SIGNIFICANT!")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("ВЕРДИКТ")
    print("=" * 70)


if __name__ == "__main__":
    main()
