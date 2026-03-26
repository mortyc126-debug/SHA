"""
ПОИСК НУЛЕВОГО CARRY OFFSET: перебор всех 32 break бит.

Для break bit=31: carry offset = 0x91002000 (4 бита, КОНСТАНТА).
Вопрос: есть ли break bit где offset = 0?

Для каждого bit b:
  1. Break: δW[3] = 2^b
  2. A-repair: δW[11] = f(break) — детерминистический
  3. Schedule target: σ₀(2^b) — для δW[18]=0
  4. Carry offset: δW[11] XOR σ₀(2^b)

Ищем: HW(offset) = 0 или минимальный.
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

def state_diff(s1, s2):
    return sum(hw(s1[i] ^ s2[i]) for i in range(8))


def measure_offset(break_bit, n_samples=500):
    """Для данного break bit: вычисляем carry offset δW[11] vs schedule target."""
    from collections import Counter

    target_gf2 = sigma0(1 << break_bit)  # σ₀(δW[3]) для этого break bit
    offsets = []

    for trial in range(n_samples):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W1f = expand_schedule(W1)

        s1 = tuple(IV)
        states1 = [s1]
        for r in range(64):
            s1 = R(s1, W1f[r], r)
            states1.append(s1)

        W2 = list(W1)
        W2[3] ^= (1 << break_bit)
        s2 = tuple(IV)
        for r in range(3):
            s2 = R(s2, W2[r], r)
        s2 = R(s2, W2[3], 3)
        for r in range(4, 12):
            W2[r] = a_repair_W(s2, states1[r + 1][0], r)
            s2 = R(s2, W2[r], r)

        actual_dw11 = W1[11] ^ W2[11]
        offset = actual_dw11 ^ target_gf2
        offsets.append(offset)

    unique = len(set(offsets))
    most_common = Counter(offsets).most_common(1)[0]
    main_offset = most_common[0]
    main_pct = most_common[1] / n_samples * 100

    return hw(main_offset), main_offset, unique, main_pct


def main():
    np.random.seed(42)

    print("=" * 70)
    print("ПОИСК НУЛЕВОГО CARRY OFFSET: все 32 break бита")
    print("=" * 70)

    print(f"\n  {'Bit':>4} {'σ₀ HW':>6} {'Offset HW':>10} {'Offset':>12} {'Unique':>7} {'Stable%':>8} {'Status':>10}")
    print(f"  {'-'*4} {'-'*6} {'-'*10} {'-'*12} {'-'*7} {'-'*8} {'-'*10}")

    best_bit = -1
    best_hw = 32
    results = []

    for bit in range(32):
        offset_hw, offset_val, unique, stable_pct = measure_offset(bit, n_samples=500)

        sig0_hw = hw(sigma0(1 << bit))

        status = ""
        if offset_hw == 0:
            status = "★★★ ZERO!"
        elif offset_hw <= 2:
            status = "★★ CLOSE"
        elif offset_hw <= 4:
            status = "★"

        if offset_hw < best_hw:
            best_hw = offset_hw
            best_bit = bit

        results.append((bit, sig0_hw, offset_hw, offset_val, unique, stable_pct))
        print(f"  {bit:4d} {sig0_hw:5d} {offset_hw:9d} {offset_val:#012x} {unique:6d} {stable_pct:7.1f}% {status}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print(f"ЛУЧШИЙ BREAK BIT: {best_bit} (offset HW = {best_hw})")
    print("=" * 70)

    if best_hw == 0:
        print(f"\n  ★★★ НУЛЕВОЙ OFFSET НАЙДЕН на bit={best_bit}!")
        print(f"  A-repair δW[11] СОВПАДАЕТ с schedule target!")
        print(f"  Расшивка работает БЕСПЛАТНО!")
        print(f"  → Reboot + Schedule lock = ОДНОВРЕМЕННО")
        print(f"  → Meeting r=12 (free) + δW[16..18]=0 (free)")
        print(f"  → Collision cost = стоимость gap r=19..60")
    elif best_hw <= 2:
        print(f"\n  ★★ БЛИЗКИЙ OFFSET (HW={best_hw}) на bit={best_bit}")
        print(f"  Нужна компенсация {best_hw} бит")
        print(f"  Возможна через перебор 2^{best_hw} вариантов")
    else:
        print(f"\n  Минимальный offset = {best_hw} бит на bit={best_bit}")
        print(f"  Нулевой offset не найден ни для одного break bit")

    # Для лучшего бита: полный тест
    print(f"\n  Детальный тест bit={best_bit}:")
    bit = best_bit
    offset_hw, offset_val, unique, stable_pct = measure_offset(bit, n_samples=2000)
    print(f"    Offset: {offset_val:#010x} (HW={offset_hw})")
    print(f"    Unique offsets: {unique}")
    print(f"    Stability: {stable_pct:.1f}%")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("Z/2^32 OFFSET (вместо GF2)")
    print("=" * 70)

    # Тот же тест но target = -σ₀(2^b) mod 2^32 (Z/2^32)
    print(f"\n  Z/2^32 targets:")
    best_z_bit = -1
    best_z_hw = 32

    for bit in range(32):
        target_z = sub32(0, sigma0(1 << bit))
        offsets = []

        for trial in range(500):
            W1 = [np.random.randint(0, 2**32) for _ in range(16)]
            W1f = expand_schedule(W1)
            s1 = tuple(IV)
            states1 = [s1]
            for r in range(64):
                s1 = R(s1, W1f[r], r)
                states1.append(s1)

            W2 = list(W1)
            W2[3] ^= (1 << bit)
            s2 = tuple(IV)
            for r in range(3):
                s2 = R(s2, W2[r], r)
            s2 = R(s2, W2[3], 3)
            for r in range(4, 12):
                W2[r] = a_repair_W(s2, states1[r + 1][0], r)
                s2 = R(s2, W2[r], r)

            actual_dw11 = W1[11] ^ W2[11]
            offset = actual_dw11 ^ (target_z & MASK32)
            offsets.append(offset)

        from collections import Counter
        mc = Counter(offsets).most_common(1)[0]
        off_hw = hw(mc[0])

        if off_hw < best_z_hw:
            best_z_hw = off_hw
            best_z_bit = bit

        if off_hw <= 4 or bit == best_bit:
            print(f"    bit={bit:2d}: Z/2^32 offset HW={off_hw}  {mc[0]:#010x}"
                  f"{'  ★★★' if off_hw == 0 else '  ★★' if off_hw <= 2 else '  ★' if off_hw <= 4 else ''}")

    print(f"\n  Best Z/2^32: bit={best_z_bit}, offset HW={best_z_hw}")


if __name__ == "__main__":
    main()
