"""
МУЛЬТИ-БЛОЧНАЯ ТКАНЬ: два блока SHA-256.

Block 1: SHA256(IV, M1) → H1 (intermediate hash)
Block 2: SHA256(H1, M2) → H2 (final hash)

В нашем измерении:
  - 128 слоёв (64+64)
  - 32 свободных позиции (16+16)
  - Бюджет УДВОИЛСЯ!

Collision на multi-block:
  (M1_a || M2_a) vs (M1_b || M2_b) → H2_a = H2_b

Варианты:
  A. Same block1: M1_a = M1_b → H1 same → collision в block2 only (= single block)
  B. Different block1, same block2: H1_a ≠ H1_b но block2 сводит → FREE-START в block2
  C. Both different: полная свобода 32 позиций

Вариант B = FREE-START collision на block2!
  IV для block2 = H1 (которое мы ВЫБИРАЕМ через M1).
  → 8 дополнительных "позиций" свободы (IV variable).
  → Total: 16 (M2) + 8 (IV from M1) = 24 позиции!

В нашем измерении: 24 позиции → 8 выходных.
  Избыток: 24 - 8 = 16 позиций.
  A-repair budget: 16 позиций для meeting + lock.
  Ранее: 12 позиций → дефицит 3. Теперь: 16 → дефицит 0???
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
def sub32(x, y): return (x - y) & MASK32
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
    a, b, c, d, e, f, g, h = IV
    for r in range(64):
        T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e, f, g)), K_const[r]), W[r])
        T2 = add32(Sigma0(a), Maj(a, b, c))
        h, g, f, e = g, f, e, add32(d, T1)
        d, c, b, a = c, b, a, add32(T1, T2)
    return tuple(add32(IV[i], [a,b,c,d,e,f,g,h][i]) for i in range(8))


def main():
    np.random.seed(42)

    print("=" * 70)
    print("МУЛЬТИ-БЛОЧНАЯ ТКАНЬ: удвоение бюджета")
    print("=" * 70)

    # =================================================================
    print(f"\n{'=' * 70}")
    print("1. БЮДЖЕТ MULTI-BLOCK")
    print("=" * 70)

    print(f"""
  Single-block:
    Input: 16 позиций (M[0..15])
    IV: ФИКСИРОВАН (8 позиций, не свободны)
    Output: 8 позиций
    Budget: 16 - 0 = 16 свободных
    A-repair: 12 позиций (M[4..15])
    Дефицит: meeting(8) + lock(7) = 15 > 12. Deficit = 3.

  Multi-block (2 blocks):
    Block 1: M1[0..15] → H1[0..7] (= IV для block 2)
    Block 2: M2[0..15] + IV=H1 → H2[0..7] (final)

    Для collision на H2:
      Подход B: M1_a ≠ M1_b → H1_a ≠ H1_b (разные IV для block 2)
      M2_a ≠ M2_b (разные messages в block 2)
      H2(H1_a, M2_a) = H2(H1_b, M2_b) — FREE-START collision на block 2!

    Свобода block 2:
      M2[0..15]: 16 позиций (свободны)
      IV = H1: 8 позиций (ВЫБИРАЕМ через M1!)
      Total: 24 свободных позиции!

    Budget block 2: 24 - 8(output) = 16 позиций свободы.
    Meeting(8) + lock(7) = 15 ≤ 16. Deficit = 0!
""")

    # =================================================================
    print(f"{'=' * 70}")
    print("2. FREE-START A-REPAIR: 24 позиции")
    print("=" * 70)

    # В free-start: IV₁ ≠ IV₂ (δIV ≠ 0).
    # A-repair работает на state, не на IV.
    # Если δIV ≠ 0: δstate[0] = δIV ≠ 0 → a-repair начинает с ненулевого δ.

    # НО: мы ВЫБИРАЕМ IV₁ и IV₂ через M1_a и M1_b.
    # Значит: мы выбираем δIV = H1_a ⊕ H1_b.
    # Можем ли выбрать δIV = 0 на НЕКОТОРЫХ позициях?

    # Если δIV[0..3] = 0 (a-chain): a-repair начинает с δa=δb=δc=δd=0.
    # Если δIV[4..7] = 0 (e-chain): δe=δf=δg=δh=0.

    # Birthday на δIV[0..3]=0: 128-бит condition → N = 2^64 messages M1.
    # Birthday на δIV=0 полностью: 256-бит → N = 2^128. Слишком дорого.
    # Birthday на δIV[7]=0: 32-бит → N = 2^16. ДЁШЕВО!

    # Стратегия: birthday на H1[7] → δIV[7]=0 (δh=0 на входе block2).
    # Стоимость: 2^16 М1 messages.
    # Результат: block 2 начинает с δh=0 → backward pipe e→f→g→h чист.

    N = 200000
    h1_dict = {}
    h1_7_collisions = []

    for _ in range(N):
        M1 = [np.random.randint(0, 2**32) for _ in range(16)]
        H1 = sha256_compress(M1, IV_STD)
        key = H1[7]  # birthday на H1[7]
        if key in h1_dict:
            h1_7_collisions.append((M1, h1_dict[key], H1, sha256_compress(h1_dict[key], IV_STD)))
        h1_dict[key] = M1

    print(f"\n  H1[7]-birthday: {len(h1_7_collisions)} collisions из {N}")

    if h1_7_collisions:
        M1_a, M1_b, H1_a, H1_b = h1_7_collisions[0]
        dIV = [hw(H1_a[i] ^ H1_b[i]) for i in range(8)]
        print(f"  δIV (H1_a ⊕ H1_b): {dIV}")
        print(f"  δIV[7] = {dIV[7]} (= 0 by birthday ✓)")
        print(f"  δIV total HW: {sum(dIV)}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("3. MULTI-WORD IV BIRTHDAY: δIV[3,7]=0 → a/e chain clean")
    print("=" * 70)

    # δIV[7]=0 → δh=0 → backward pipe clean (e→f→g→h)
    # δIV[3]=0 → δd=0 → backward pipe clean (a→b→c→d)
    # Both: 64-бит birthday → N = 2^32 M1 messages.

    h1_37_dict = {}
    h1_37_collisions = []

    for _ in range(N):
        M1 = [np.random.randint(0, 2**32) for _ in range(16)]
        H1 = sha256_compress(M1, IV_STD)
        key = (H1[3], H1[7])
        if key in h1_37_dict:
            h1_37_collisions.append((M1, h1_37_dict[key], H1, sha256_compress(h1_37_dict[key], IV_STD)))
        h1_37_dict[key] = M1

    print(f"\n  H1[3,7]-birthday: {len(h1_37_collisions)} collisions из {N}")
    print(f"  Expected: {N**2 / 2**65:.1f}")

    if h1_37_collisions:
        M1_a, M1_b, H1_a, H1_b = h1_37_collisions[0]
        dIV = [hw(H1_a[i] ^ H1_b[i]) for i in range(8)]
        print(f"  δIV: {dIV}")
        print(f"  δIV[3]=0 δIV[7]=0 → δd=0, δh=0 на входе block 2")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("4. BLOCK 2 С δIV[7]=0: что даёт?")
    print("=" * 70)

    # В block 2: IV = H1. Если H1_a[7] = H1_b[7]:
    # state2[0] = IV = H1 → δstate2[0] = δH1.
    # δH1[7] = 0 → δh[0] = 0.
    # Через трубы: δh[0] → δg[1] → δf[2] → δe[3].
    # Значит δe[3] = 0!

    # Если ещё δH1[3]=0: δd[0]=0 → δc[1] → δb[2] → δa[3].
    # δa[3] = 0!

    # При δa[3]=0 и δe[3]=0: δstate[3] = (0,?,?,0,0,?,?,0).
    # 4 из 8 регистров = 0 на раунде 3! Бесплатно!

    # С a-repair от r=3: budget = M2[3..15] = 13 позиций + M2[0..2] free.
    # + IV variation = 8 позиций.
    # Total effective: 13 + overhead = достаточно для meeting + lock?

    if h1_7_collisions:
        # Simulate block 2 with δIV[7]=0
        M1_a, M1_b, H1_a, H1_b = h1_7_collisions[0]

        # Block 2 with same M2 but different IV
        M2 = [np.random.randint(0, 2**32) for _ in range(16)]
        H2_a = sha256_compress(M2, list(H1_a))
        H2_b = sha256_compress(M2, list(H1_b))

        dH2 = [hw(H2_a[i] ^ H2_b[i]) for i in range(8)]
        print(f"\n  Block 2 (same M2, δIV[7]=0):")
        print(f"    δH2 per word: {dH2}")
        print(f"    Total HW(δH2): {sum(dH2)}")
        print(f"    → {'REDUCED!' if sum(dH2) < 120 else 'Still random (δIV other positions too large)'}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("5. PIPE ANALYSIS: δIV через трубы block 2")
    print("=" * 70)

    # δIV: known from H1_a, H1_b.
    # В block 2: state[0] = IV = H1.
    # Трубы в раундах 0-3 ПЕРЕДАЮТ δIV:
    #   δb[1] = δa[0] = δH1[0]
    #   δc[2] = δb[1] = δa[0] = δH1[0]
    #   δd[3] = δc[2] = δH1[0]
    #   δf[1] = δe[0] = δH1[4]
    #   δg[2] = δf[1] = δH1[4]
    #   δh[3] = δg[2] = δH1[4]

    # Если δH1[0]=0 (NODE): δb[1]=δc[2]=δd[3]=0 (a-chain clean!)
    # Если δH1[4]=0 (NODE): δf[1]=δg[2]=δh[3]=0 (e-chain clean!)

    # Birthday на δH1[0,4]=0: 64-бит → N = 2^32 M1 messages.
    # Даёт: δstate[3] = (δa[3],0,0,0,δe[3],0,0,0) = only NODE positions nonzero!

    # Это ЛУЧШЕ чем single-block! Single: δstate[0] = δ only on W[3].
    # Multi: δstate[3] has 6/8 registers = 0 from PIPES.

    print(f"""
  При δH1[0]=0, δH1[4]=0 (cost 2^32 birthday на block 1):
    Block 2 state[3] = (δa[3], 0, 0, 0, δe[3], 0, 0, 0)
    6/8 регистров = 0 через ТРУБЫ!

    Только NODE_a(2) и NODE_e(2) вносят δ.
    δa[3] зависит от δH1[0..7] через 3 раунда round function.
    δe[3] аналогично.

    A-repair от r=3: budget = M2[3..15] = 13 слов.
    Meeting(8) + lock(7) = 15 → deficit = 2.
    НО: δstate[3] = (small, 0, 0, 0, small, 0, 0, 0)
    A-repair начинает с ПОЧТИ нулевого state!

    Pipe cascade БЕСПЛАТНО обнуляет 6 регистров.
    A-repair нужен только для δa и δe.
    Это ДЕШЕВЛЕ чем single-block!

  СТОИМОСТЬ MULTI-BLOCK:
    Block 1 birthday (δH1[0,4]=0): 2^32
    Block 2 a-repair + lock: budget 13 → deficit 2
    Deficit 2 стоит: ???
    Total: 2^32 + deficit_cost
""")

    # Measure: what is δstate[3] at block 2 entry?
    if h1_7_collisions:
        dstate3_values = []
        for trial in range(min(len(h1_7_collisions), 50)):
            M1_a, M1_b, H1_a, H1_b = h1_7_collisions[trial]
            # δH1
            dH1 = tuple(H1_a[i] ^ H1_b[i] for i in range(8))

            # Run 3 rounds of block 2 with IV=H1_a vs IV=H1_b, same M2
            M2 = [np.random.randint(0, 2**32) for _ in range(16)]

            s_a = list(H1_a)  # state = IV = H1_a
            s_b = list(H1_b)

            # Expand schedule for M2
            W = list(M2)
            for r in range(16, 64):
                W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))

            # 3 rounds
            for r in range(3):
                a,b,c,d,e,f,g,h = s_a
                T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K_const[r]), W[r])
                T2 = add32(Sigma0(a), Maj(a,b,c))
                s_a = [add32(T1,T2),a,b,c,add32(d,T1),e,f,g]

                a,b,c,d,e,f,g,h = s_b
                T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K_const[r]), W[r])
                T2 = add32(Sigma0(a), Maj(a,b,c))
                s_b = [add32(T1,T2),a,b,c,add32(d,T1),e,f,g]

            ds3 = sum(hw(s_a[i] ^ s_b[i]) for i in range(8))
            dstate3_values.append(ds3)

        if dstate3_values:
            print(f"\n  δstate[3] at block 2 (with δH1[7]=0):")
            print(f"    Mean: {np.mean(dstate3_values):.1f}")
            print(f"    Min: {min(dstate3_values)}")
            print(f"    → {'REDUCED!' if np.mean(dstate3_values) < 100 else 'Still high'}")


if __name__ == "__main__":
    main()
