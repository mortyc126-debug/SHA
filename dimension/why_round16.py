"""
ПОЧЕМУ r=16? Связь границы безопасности с schedule.

Факты:
  - rank(CE) = 256 начиная с r=16 (SECURE)
  - rank(CE) = 224 при r=15 (32 collisions)
  - W[0..15] = свободные слова. W[16+] = schedule-expanded.
  - r=16 = ПЕРВЫЙ раунд с schedule-expanded W[16].

Совпадение? Или schedule expansion = причина security?

Тесты:
  1. SHA-256 БЕЗ schedule expansion (W[16..]=0): rank(CE) при r=16?
  2. SHA-256 с W[16..]=random (не schedule): rank(CE)?
  3. Где именно коллизионные δW живут? В W[0..7] или W[8..15]?
"""

import numpy as np
import struct

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

def hw(x): return bin(x).count('1')
def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def add32(x, y): return (x + y) & MASK32
def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)

def sha256_custom_schedule(W16, n_rounds, schedule_mode='standard'):
    """SHA-256 with configurable schedule."""
    W = list(W16)
    if schedule_mode == 'standard':
        for r in range(16, max(n_rounds, 16)):
            W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))
    elif schedule_mode == 'zero':
        W.extend([0] * max(0, n_rounds - 16))
    elif schedule_mode == 'copy':
        while len(W) < n_rounds:
            W.append(W[len(W) % 16])

    a,b,c,d,e,f,g,h = IV
    for r in range(n_rounds):
        T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K_const[r%64]), W[r])
        T2 = add32(Sigma0(a), Maj(a,b,c))
        h,g,f,e = g,f,e,add32(d,T1)
        d,c,b,a = c,b,a,add32(T1,T2)
    return tuple(add32(IV[i], [a,b,c,d,e,f,g,h][i]) for i in range(8))


def compute_ce_rank(n_rounds, W_base, schedule_mode='standard'):
    H_base = sha256_custom_schedule(W_base, n_rounds, schedule_mode)
    T_rows = []
    for word in range(16):
        for bit in range(32):
            W_mod = list(W_base)
            W_mod[word] ^= (1 << bit)
            H_mod = sha256_custom_schedule(W_mod, n_rounds, schedule_mode)
            row = []
            for w in range(8):
                d = H_base[w] ^ H_mod[w]
                for b in range(32): row.append((d >> b) & 1)
            T_rows.append(row)

    T = np.array(T_rows, dtype=np.uint8)
    rank_T = np.linalg.matrix_rank(T.astype(float))
    if rank_T < 256: return rank_T, None

    M = T.T.copy()
    pivots = []; row = 0
    for col in range(512):
        found = False
        for r in range(row, 256):
            if M[r, col] == 1:
                M[[row, r]] = M[[r, row]]; found = True; break
        if not found: continue
        pivots.append(col)
        for r in range(256):
            if r != row and M[r, col] == 1: M[r] ^= M[row]
        row += 1

    free_vars = [c for c in range(512) if c not in pivots]
    ces = []
    for fc in free_vars[:256]:
        x = np.zeros(512, dtype=np.uint8); x[fc] = 1
        for i in range(len(pivots)-1, -1, -1):
            pc = pivots[i]; val = np.uint8(0)
            for j in range(512):
                if j != pc: val ^= (M[i,j] & x[j])
            x[pc] = val
        dW = [0]*16
        for word in range(16):
            for bit in range(32):
                if x[word*32+bit]: dW[word] ^= (1<<bit)
        W2 = [W_base[i]^dW[i] for i in range(16)]
        H2 = sha256_custom_schedule(W2, n_rounds, schedule_mode)
        ce = []
        for w in range(8):
            d = H_base[w] ^ H2[w]
            for b in range(32): ce.append((d >> b) & 1)
        ces.append(ce)

    CE = np.array(ces[:256], dtype=np.uint8)
    return rank_T, np.linalg.matrix_rank(CE.astype(float))


def main():
    np.random.seed(42)

    print("=" * 70)
    print("ПОЧЕМУ r=16? Связь schedule и безопасности")
    print("=" * 70)

    W_base = [np.random.randint(0, 2**32) for _ in range(16)]

    # =================================================================
    print(f"\n{'=' * 70}")
    print("1. STANDARD vs NO SCHEDULE: rank(CE) по раундам")
    print("=" * 70)

    print(f"\n  {'Rounds':>6} {'Standard':>10} {'No sched':>10} {'Copy sched':>11}")
    print(f"  {'-'*6} {'-'*10} {'-'*10} {'-'*11}")

    for n_r in [8, 10, 12, 14, 15, 16, 17, 18, 20, 24, 32]:
        rT_s, rCE_s = compute_ce_rank(n_r, W_base, 'standard')
        rT_z, rCE_z = compute_ce_rank(n_r, W_base, 'zero')
        rT_c, rCE_c = compute_ce_rank(n_r, W_base, 'copy')

        def fmt(rT, rCE):
            if rCE is None: return f"T={rT}"
            return f"{rCE}"

        print(f"  {n_r:6d} {fmt(rT_s, rCE_s):>10} {fmt(rT_z, rCE_z):>10} {fmt(rT_c, rCE_c):>11}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("2. ГДЕ ЖИВУТ КОЛЛИЗИОННЫЕ δW?")
    print("=" * 70)

    # Для 8-round collision: δW = 1 бит. В каком слове?
    H_base = sha256_custom_schedule(W_base, 12, 'standard')

    collision_words = {w: 0 for w in range(16)}
    for word in range(16):
        for bit in range(32):
            dW = [0]*16; dW[word] = 1 << bit
            W2 = [W_base[i]^dW[i] for i in range(16)]
            H2 = sha256_custom_schedule(W2, 12, 'standard')
            if H_base == H2:
                collision_words[word] += 1

    print(f"\n  12-round collisions по входным словам (δW=1 бит):")
    for w in range(16):
        bar = "█" * collision_words[w]
        used = "← используется в round" if w < 12 else "← NOT used (W[12+])"
        print(f"    W[{w:2d}]: {collision_words[w]:3d} collisions  {bar}  {used}")

    # Повторяем для 15 rounds
    print(f"\n  15-round collisions по входным словам:")
    H_base_15 = sha256_custom_schedule(W_base, 15, 'standard')
    coll_15 = {w: 0 for w in range(16)}
    for word in range(16):
        for bit in range(32):
            dW = [0]*16; dW[word] = 1 << bit
            W2 = [W_base[i]^dW[i] for i in range(16)]
            H2 = sha256_custom_schedule(W2, 15, 'standard')
            if H_base_15 == H2:
                coll_15[word] += 1

    for w in range(16):
        bar = "█" * coll_15[w]
        print(f"    W[{w:2d}]: {coll_15[w]:3d} collisions  {bar}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("3. ОБЪЯСНЕНИЕ: почему r=16 = граница")
    print("=" * 70)

    total_12 = sum(collision_words.values())
    words_with_coll = [w for w in range(16) if collision_words[w] > 0]
    total_15 = sum(coll_15.values())
    words_15 = [w for w in range(16) if coll_15[w] > 0]

    print(f"""
  12-round: {total_12} collisions в словах {words_with_coll}
  15-round: {total_15} collisions в словах {words_15}

  PATTERN: collisions живут в W[r..15] (словах НЕ использованных до раунда r).

  SHA-256 r-round использует W[0..r-1] напрямую.
  W[r..15] НЕ входят в раунды 0..r-1 как W[r].
  НО: W[r..15] входят в SCHEDULE для W[16+].

  При r < 16: W[r..15] влияют на H ТОЛЬКО через schedule (W[16+]).
  Если r < 16: нет раундов использующих W[16+]
  → W[r..15] = "мёртвые" входы. Flip бит → carry error = 0.

  При r = 16: W[16] = schedule(W[0..15]) используется в раунде 16.
  Теперь δW[0..15] влияет через W[16]:
    δW[16] = σ₁(δW[14]) + δW[9] + σ₀(δW[1]) + δW[0]
  Carry в schedule ADD → rank(CE) прыгает!

  ВЫВОД: rank(CE) < 256 для r < 16 потому что
  W[r..15] = "мёртвые позиции" (не используются ни в раундах, ни в schedule).
  W[r..15] дают 32×(16-r) бит "бесплатных" δ = ker(CE).

  При r=16: все 16 слов используются → нет мёртвых → rank(CE)=256.
""")


if __name__ == "__main__":
    main()
