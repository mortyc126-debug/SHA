"""
ПОТОК ЭНТРОПИИ: информация через ткань.

В нашем измерении: ткань преобразует вход в выход.
Сколько "информации" о входе ОСТАЁТСЯ на каждом раунде?

Измеряем: mutual information proxy между W[0] и state[r].
Если MI высокая: state[r] "помнит" W[0].
Если MI ≈ 0: state[r] "забыл" W[0].

Также: сколько бит state[r] НЕЗАВИСИМЫ (effective entropy)?

Entropy flow = fundamental transport property нашей ткани.
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


def sha256_states(W16):
    """Return ALL intermediate states."""
    W = list(W16)
    for r in range(16, 64):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))
    a,b,c,d,e,f,g,h = IV
    states = [(a,b,c,d,e,f,g,h)]
    for r in range(64):
        T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K_const[r]), W[r])
        T2 = add32(Sigma0(a), Maj(a,b,c))
        h,g,f,e = g,f,e,add32(d,T1)
        d,c,b,a = c,b,a,add32(T1,T2)
        states.append((a,b,c,d,e,f,g,h))
    return states


def main():
    np.random.seed(42)

    print("=" * 70)
    print("ПОТОК ЭНТРОПИИ ЧЕРЕЗ ТКАНЬ")
    print("=" * 70)

    # =================================================================
    print(f"\n{'=' * 70}")
    print("1. ПАМЯТЬ: corr(W[0], state[r][0]) по раундам")
    print("=" * 70)

    # Для каждого раунда: корреляция между W[0] (input) и a[r] (state)
    N = 5000
    memory_by_round = []

    w0_lsb = []
    state_a_lsb = [[] for _ in range(65)]

    for _ in range(N):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        states = sha256_states(W)
        w0_lsb.append(W[0] & 0xFF)
        for r in range(65):
            state_a_lsb[r].append(states[r][0] & 0xFF)

    for r in range(65):
        if np.std(state_a_lsb[r]) > 0:
            corr = np.corrcoef(w0_lsb, state_a_lsb[r])[0, 1]
            memory_by_round.append(abs(corr) if not np.isnan(corr) else 0)
        else:
            memory_by_round.append(0)

    print(f"\n  |corr(W[0], a[r])| по раундам:")
    for r in [0, 1, 2, 3, 4, 5, 8, 12, 16, 20, 32, 48, 64]:
        bar = "█" * int(memory_by_round[r] * 50)
        print(f"    r={r:2d}: {memory_by_round[r]:.4f}  {bar}")

    # Memory half-life: round where corr drops below 0.01
    half_life = 64
    for r in range(65):
        if memory_by_round[r] < 0.01 and r > 0:
            half_life = r
            break

    print(f"\n  Memory half-life: r={half_life} (corr < 0.01)")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("2. EFFECTIVE ENTROPY state[r]: сколько бит независимы?")
    print("=" * 70)

    # Proxy: для 10K random W, сколько UNIQUE state[r][0] (a-register)?
    # Если все unique → 32 бит entropy. Если repeats → меньше.

    # Метод: вместо подсчёта unique (ограничен sample size),
    # измеряем collision rate. P(collision) = 1/2^H → H = -log2(P).
    # Для a[r]: P(a[r]_i = a[r]_j) для random i,j.

    for r in [0, 1, 2, 4, 8, 16, 32, 64]:
        N_ent = 5000
        a_values = []
        for _ in range(N_ent):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            states = sha256_states(W)
            a_values.append(states[r][0])

        # Collision count among pairs
        from collections import Counter
        cnt = Counter(a_values)
        collisions = sum(c * (c - 1) // 2 for c in cnt.values())
        total_pairs = N_ent * (N_ent - 1) // 2
        p_coll = collisions / total_pairs if total_pairs > 0 else 1

        if p_coll > 0 and p_coll < 1:
            est_entropy = -np.log2(p_coll)
        elif p_coll == 0:
            est_entropy = np.log2(total_pairs)  # lower bound
        else:
            est_entropy = 0

        unique = len(set(a_values))
        print(f"    r={r:2d}: unique={unique:>5}/{N_ent}, collisions={collisions}, "
              f"est. entropy ≈ {min(est_entropy, 32):.1f} bits (max 32)")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("3. INFORMATION INJECTION: W[r] → state[r] per round")
    print("=" * 70)

    # How much does W[r] add to state?
    # Proxy: corr(W[r], a[r+1]) — does W[r] influence state?

    print(f"\n  corr(W[r], a[r+1]) — injection per round:")
    for r_target in [0, 1, 2, 4, 8, 12, 15]:
        w_vals = []
        a_vals = []
        for _ in range(5000):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            states = sha256_states(W)
            w_vals.append(W[r_target] & 0xFF)
            a_vals.append(states[r_target + 1][0] & 0xFF)

        corr = abs(np.corrcoef(w_vals, a_vals)[0, 1])
        bar = "█" * int(corr * 50)
        print(f"    W[{r_target:2d}] → a[{r_target+1:2d}]: corr={corr:.4f}  {bar}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("4. ENTROPY FLOW: input → output")
    print("=" * 70)

    # Total information path: 512 bits input → 256 bits output.
    # Per-word: 32 bits W[r] → 256 bits state.
    # After 64 rounds: 512 bits compressed to 256.

    # Proxy: rank of Jacobian at each round
    W_base = [np.random.randint(0, 2**32) for _ in range(16)]

    print(f"\n  rank(Jacobian) of W → state[r]:")
    for r in [1, 2, 4, 8, 16, 32, 64]:
        states_base = sha256_states(W_base)
        T_rows = []
        for word in range(16):
            for bit in range(32):
                W_mod = list(W_base)
                W_mod[word] ^= (1 << bit)
                states_mod = sha256_states(W_mod)
                row = []
                for w in range(8):
                    d = states_base[r][w] ^ states_mod[r][w]
                    for b in range(32): row.append((d >> b) & 1)
                T_rows.append(row)

        T = np.array(T_rows, dtype=np.uint8)
        rank = np.linalg.matrix_rank(T.astype(float))
        bar = "█" * (rank // 8)
        print(f"    r={r:2d}: rank={rank:3d}/256  {bar}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("5. ENTROPY MODEL нашего измерения")
    print("=" * 70)

    print(f"""
  ПОТОК ЭНТРОПИИ В ТКАНИ SHA-256:

  INJECTION:
    W[r] вносит 32 бита на раунде r.
    Injection rate = 32 бит/раунд (r=0..15).
    Total injection = 512 бит (16 слов × 32).

  DIFFUSION:
    rank(W→state) растёт:
    r=1: ~32 (одно слово), r=8: ~256 (saturated).
    State entropy: saturates at 256 бит за ~8 раундов.

  MEMORY:
    corr(W[0], state[r]) → 0 за {half_life} раундов.
    Input information "забывается" — но НЕ теряется.
    Она ПЕРЕМЕШИВАЕТСЯ со state.

  COMPRESSION:
    512 бит вход → 256 бит выход.
    Compression ratio = 2:1.
    Finalization (state + IV): preserves 256 бит.

  В нашем измерении: entropy FLOW = pipe + node transport.
    Pipes: transmit entropy without loss (κ=1).
    Nodes: MIX entropy with new input (κ<1 for differential).

  Entropy is NEVER LOST — only MIXED.
  SHA-256 = entropy mixer, not entropy destroyer.
""")


if __name__ == "__main__":
    main()
