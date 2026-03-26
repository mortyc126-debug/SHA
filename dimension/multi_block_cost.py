"""
MULTI-BLOCK COLLISION: полная стоимость.

Block 1: birthday δH1[0,4,3,7]=0 (128 бит) → cost 2^64
  → δIV = (0,δ,δ,0,0,δ,δ,0) — только pipe-позиции ненулевые
  → Block 2 pipes чисты: δd[3]=0, δh[3]=0

Block 2: a-repair с δIV частично нулевым.
  Budget: 16(M2) + 8(IV partial) → meeting + lock.

Ключевой вопрос: при δIV[0,3,4,7]=0, чему равен δstate[3] в block 2?
Если мал → a-repair дешёвый → total < 2^128?
"""

import numpy as np
import struct, hashlib

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
IV_STD = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
          0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]

def hw(x): return bin(x).count('1')
def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def add32(x, y): return (x + y) & MASK32
def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)

def sha256_compress(W16, IV):
    W = list(W16)
    for r in range(16, 64):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))
    a,b,c,d,e,f,g,h = IV
    for r in range(64):
        T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K_const[r]), W[r])
        T2 = add32(Sigma0(a), Maj(a,b,c))
        h,g,f,e = g,f,e,add32(d,T1)
        d,c,b,a = c,b,a,add32(T1,T2)
    return tuple(add32(IV[i], [a,b,c,d,e,f,g,h][i]) for i in range(8))


def main():
    np.random.seed(42)

    print("=" * 70)
    print("MULTI-BLOCK COLLISION: полная стоимость")
    print("=" * 70)

    # =================================================================
    print(f"\n{'=' * 70}")
    print("1. СТОИМОСТЬ BLOCK 1 BIRTHDAY по разным conditions")
    print("=" * 70)

    # Генерируем пары M1 с birthday на разных подмножествах H1
    conditions = [
        ("δH1[7]=0", lambda h1a, h1b: h1a[7]==h1b[7], 32),
        ("δH1[3,7]=0", lambda h1a, h1b: h1a[3]==h1b[3] and h1a[7]==h1b[7], 64),
        ("δH1[0,7]=0", lambda h1a, h1b: h1a[0]==h1b[0] and h1a[7]==h1b[7], 64),
        ("δH1[0,4]=0", lambda h1a, h1b: h1a[0]==h1b[0] and h1a[4]==h1b[4], 64),
        ("δH1[0,3,4,7]=0", lambda h1a, h1b: h1a[0]==h1b[0] and h1a[3]==h1b[3] and h1a[4]==h1b[4] and h1a[7]==h1b[7], 128),
    ]

    for name, cond_fn, n_bits in conditions:
        cost = 2**(n_bits/2)
        print(f"  {name}: birthday {n_bits} бит → cost 2^{n_bits/2:.0f}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("2. δstate[r] В BLOCK 2 при разных δIV conditions")
    print("=" * 70)

    # Simulate: random δIV с разными zero-conditions
    N = 2000

    for name, zero_positions in [
        ("δIV[7]=0", [7]),
        ("δIV[3,7]=0", [3,7]),
        ("δIV[0,4]=0", [0,4]),
        ("δIV[0,3,4,7]=0", [0,3,4,7]),
        ("δIV[0,1,2,3,4,5,6,7]=0", list(range(8))),  # full collision on H1
    ]:
        dstate_by_round = np.zeros(65)
        count = 0

        for _ in range(N):
            # Random δIV with specified positions = 0
            IV_a = [np.random.randint(0, 2**32) for _ in range(8)]
            IV_b = list(IV_a)
            for i in range(8):
                if i not in zero_positions:
                    IV_b[i] = np.random.randint(0, 2**32)

            # Same M2 for both
            M2 = [np.random.randint(0, 2**32) for _ in range(16)]
            W = list(M2)
            for r in range(16, 64):
                W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))

            s_a = list(IV_a)
            s_b = list(IV_b)

            for r in range(64):
                a,b,c,d,e,f,g,h = s_a
                T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K_const[r]), W[r])
                T2 = add32(Sigma0(a), Maj(a,b,c))
                s_a = [add32(T1,T2),a,b,c,add32(d,T1),e,f,g]

                a,b,c,d,e,f,g,h = s_b
                T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K_const[r]), W[r])
                T2 = add32(Sigma0(a), Maj(a,b,c))
                s_b = [add32(T1,T2),a,b,c,add32(d,T1),e,f,g]

                ds = sum(hw(s_a[i] ^ s_b[i]) for i in range(8))
                dstate_by_round[r+1] += ds

            count += 1

        dstate_by_round /= count

        # Finalization
        # H = state[64] + IV. δH depends on both δstate[64] and δIV.
        print(f"\n  {name}:")
        print(f"    δstate[1]: {dstate_by_round[1]:.1f}")
        print(f"    δstate[3]: {dstate_by_round[3]:.1f}")
        print(f"    δstate[8]: {dstate_by_round[8]:.1f}")
        print(f"    δstate[16]: {dstate_by_round[16]:.1f}")
        print(f"    δstate[64]: {dstate_by_round[64]:.1f}")

        # Cost analysis
        n_zero = len(zero_positions)
        block1_cost = n_zero * 32 / 2  # birthday on n_zero × 32 bits
        remaining_bits = (8 - n_zero) * 32  # remaining δIV for finalization
        ds64 = dstate_by_round[64]

        print(f"    Block1 cost: 2^{block1_cost:.0f}")
        print(f"    Remaining δIV: {remaining_bits} бит")
        print(f"    δstate[64]: {ds64:.1f}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("3. ФИНАЛЬНАЯ БУХГАЛТЕРИЯ MULTI-BLOCK")
    print("=" * 70)

    print(f"""
  SINGLE-BLOCK:
    Budget: 12 free positions
    Deficit: 3 (meeting 8 + lock 7 > 12)
    Total: 2^128 (birthday C^4)

  MULTI-BLOCK (free-start на block 2):
    Block 1: birthday на δH1[pipe_positions] = 0
    Block 2: a-repair с custom δIV

    Вариант A: δH1[7]=0 (cost 2^16)
      Block 2: δstate[64] ≈ 128. No gain.
      Total: 2^16 + 2^128 = 2^128.

    Вариант B: δH1[0,4]=0 (cost 2^32)
      Block 2: δstate[3] has 4/8 regs = 0 (pipe cascade!)
      НО δstate[64] ≈ 128 (remaining δIV too large)
      Total: 2^32 + 2^128 = 2^128.

    Вариант C: δH1[0,1,2,3,4,5,6,7]=0 (= H1 collision, cost 2^128)
      Block 2: δIV=0 → standard single-block
      Total: 2^128 + single-block = 2^128. Circular!

    Вариант D: δH1[0,3,4,7]=0 (cost 2^64)
      Block 2: δd[3]=0, δh[3]=0 → pipe cascades clean
      δstate[3] = (δa, δb, δc, 0, δe, δf, δg, 0)
      Still 6 nonzero regs from δH1[1,2,5,6]
      δstate[64] ≈ 128. No gain.
      Total: 2^64 + 2^128 = 2^128.

  ПРОБЛЕМА: δIV на НЕПОГАШЕННЫХ позициях (H1[1,2,5,6]) = random.
  Они вносят full δ в block 2.
  Для collision в block 2: нужно погасить ВСЮ δIV = all H1 collision.
  → Multi-block = two serial H1-collisions. Cost = max(2^128, 2^128) = 2^128.

  Multi-block НЕ ПОМОГАЕТ потому что:
    Free-start на block 2 = collision на block 1.
    Мы заменяем одну collision другой.
    "Бюджет 24" — иллюзия: 8 extra positions = 8 more conditions.
""")


if __name__ == "__main__":
    main()
