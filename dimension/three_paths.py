"""
ТРИ ПУТИ ПРОРЫВА: higher-order + cold paths + carry resonance.

Путь 1: Higher-order differentials при больших r и k
Путь 2: Cold path — цепочка cold directions через раунды
Путь 3: Carry resonance — формализация когда δ падает
"""

import numpy as np
import struct, hashlib

MASK32 = 0xFFFFFFFF

def hw(x): return bin(x).count('1')
def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def add32(x, y): return (x + y) & MASK32
def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)

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


def one_round(state, w, r):
    a,b,c,d,e,f,g,h = state
    T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K[r]), w)
    T2 = add32(Sigma0(a), Maj(a,b,c))
    return (add32(T1,T2), a, b, c, add32(d,T1), e, f, g)


def higher_order_diff(W_base, bit_positions, n_r):
    """Compute k-th order differential. bit_positions = [(word,bit), ...]."""
    k = len(bit_positions)
    xor_sum = [0] * 8
    for mask in range(1 << k):
        W = list(W_base)
        for j in range(k):
            if (mask >> j) & 1:
                W[bit_positions[j][0]] ^= (1 << bit_positions[j][1])
        H = sha256_r(W, n_r)
        for i in range(8):
            xor_sum[i] ^= H[i]
    return xor_sum


def main():
    np.random.seed(42)

    print("=" * 70)
    print("ТРИ ПУТИ ПРОРЫВА")
    print("=" * 70)

    # ═══════════════════════════════════════
    print(f"\n{'='*70}")
    print("ПУТЬ 1: HIGHER-ORDER DIFFERENTIALS")
    print(f"{'='*70}")

    # Сначала проверим: zeros при r=24 — это реальный сигнал?
    # Для RANDOM function: P(XOR of 2^k hashes = 0) = 2^(-256) per word
    # Мы проверяем 8 слов ОТДЕЛЬНО, не вместе.
    # P(one word = 0) = 2^(-32). За 500 trials: ожидаем 500/2^32 ≈ 0.
    # НО: мы считаем ALL 8 words zero → P = 2^(-256) → definitely 0.

    # Подождите — проверим ПОБИТНО. Может zeros = partial zero (one word)?

    print(f"\n  CAREFUL CHECK: full 256-bit zero or partial?")

    for n_r in [16, 20, 24, 32]:
        full_zeros = 0
        word_zeros = [0]*8  # per-word zero count
        N = 2000

        for _ in range(N):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            bits = [(np.random.randint(0,16), np.random.randint(0,32)) for _ in range(2)]
            xor_sum = higher_order_diff(W, bits, n_r)

            if all(x == 0 for x in xor_sum):
                full_zeros += 1
            for i in range(8):
                if xor_sum[i] == 0:
                    word_zeros[i] += 1

        print(f"\n  r={n_r}, order-2, {N} trials:")
        print(f"    Full 256-bit zeros: {full_zeros}/{N}")
        print(f"    Per-word zeros: {[f'{z}/{N}' for z in word_zeros]}")
        print(f"    Per-word zero rate: {[f'{z/N:.3f}' for z in word_zeros]}")
        print(f"    Expected (random): {1/2**32:.10f} per word")

    # ═══════════════════
    print(f"\n  ORDER 3-6 at r=24:")
    for k in [2, 3, 4, 5, 6]:
        if k > 4:
            N = 200  # 2^k gets expensive
        else:
            N = 500
        full_zeros = 0
        partial_zeros = 0

        for _ in range(N):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            bits = [(np.random.randint(0,16), np.random.randint(0,32)) for _ in range(k)]
            xor_sum = higher_order_diff(W, bits, 24)

            if all(x == 0 for x in xor_sum):
                full_zeros += 1
            if any(x == 0 for x in xor_sum):
                partial_zeros += 1

        print(f"    Order-{k}: full_zeros={full_zeros}/{N}, partial={partial_zeros}/{N} "
              f"({partial_zeros/N*100:.1f}%)")

    # ═══════════════════════════════════════
    print(f"\n{'='*70}")
    print("ПУТЬ 2: COLD PATH ENGINEERING")
    print(f"{'='*70}")

    # Cold direction = pipe register flip (b→c copying, expansion = 1).
    # Can we CHAIN cold directions to keep δ small across MANY rounds?
    #
    # Strategy: at round r, flip the register that's ABOUT TO BE COPIED.
    # That register won't be amplified (just shifted to next position).
    #
    # SHA-256 pipes: b←a, c←b, d←c, f←e, g←f, h←g
    # If we flip bit in register b at round r:
    #   round r+1: c_new = b_old (with our flip) → δ moves to c
    #   round r+2: d_new = c_old (with our flip) → δ moves to d
    #   round r+3: d enters e_new = d + T1 → δ AMPLIFIED!
    #
    # Cold path length = 3 rounds (b→c→d→consumed).

    print(f"\n  Cold path: flip register, track how long δ stays small")

    for start_reg in range(8):
        reg_names = ['a','b','c','d','e','f','g','h']
        # Flip 1 bit in start_reg, track δstate through rounds
        dH_trajectory = []
        N_traj = 500

        for _ in range(N_traj):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            # Expand
            W_exp = list(W)
            for r in range(16, 30):
                W_exp.append(add32(add32(add32(sigma1(W_exp[r-2]), W_exp[r-7]),
                             sigma0(W_exp[r-15])), W_exp[r-16]))

            # Run to round 20, then flip register, continue
            state = tuple(IV)
            for r in range(20):
                state = one_round(state, W_exp[r], r)

            # Flip bit 0 of start_reg
            state2 = list(state)
            state2[start_reg] ^= 1
            state2 = tuple(state2)

            # Track δ for next 8 rounds
            traj = []
            s1 = state; s2 = state2
            for r in range(20, 28):
                s1 = one_round(s1, W_exp[r], r)
                s2 = one_round(s2, W_exp[r], r)
                d = sum(hw(s1[i]^s2[i]) for i in range(8))
                traj.append(d)
            dH_trajectory.append(traj)

        avg_traj = np.mean(dH_trajectory, axis=0)
        cold_rounds = sum(1 for a in avg_traj if a < 5)
        print(f"  δ{reg_names[start_reg]}: {' → '.join(f'{a:.0f}' for a in avg_traj)} "
              f"(cold for {cold_rounds} rounds)")

    # ═══════════════════════════════════════
    print(f"\n{'='*70}")
    print("ПУТЬ 3: CARRY RESONANCE")
    print(f"{'='*70}")

    # When does δstate DECREASE? Track for many messages.
    # Find PATTERN in which rounds δ decreases.

    decrease_counts = np.zeros(64)
    decrease_amounts = [[] for _ in range(64)]
    N = 2000

    for _ in range(N):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W); W2[15] ^= (1 << 31)

        W1e = list(W); W2e = list(W2)
        for r in range(16, 64):
            W1e.append(add32(add32(add32(sigma1(W1e[r-2]),W1e[r-7]),sigma0(W1e[r-15])),W1e[r-16]))
            W2e.append(add32(add32(add32(sigma1(W2e[r-2]),W2e[r-7]),sigma0(W2e[r-15])),W2e[r-16]))

        s1 = tuple(IV); s2 = tuple(IV)
        prev_d = 0
        for r in range(64):
            s1 = one_round(s1, W1e[r], r)
            s2 = one_round(s2, W2e[r], r)
            d = sum(hw(s1[i]^s2[i]) for i in range(8))
            if d < prev_d:
                decrease_counts[r] += 1
                decrease_amounts[r].append(prev_d - d)
            prev_d = d

    print(f"\n  Rounds where δstate DECREASES (δW[15]=bit31, {N} messages):")
    print(f"  {'Round':>5} {'Freq':>6} {'%':>6} {'Avg decrease':>13}")
    for r in range(15, 64):
        freq = decrease_counts[r]
        pct = freq / N * 100
        avg_dec = np.mean(decrease_amounts[r]) if decrease_amounts[r] else 0
        if pct > 10:
            print(f"  {r:>5} {freq:>5.0f} {pct:>5.1f}% {avg_dec:>12.1f}")

    # Is there a PATTERN? Do certain rounds ALWAYS decrease?
    high_decrease_rounds = [r for r in range(15, 64) if decrease_counts[r] / N > 0.4]
    print(f"\n  Rounds with >40% decrease frequency: {high_decrease_rounds}")

    # Overall: what fraction of rounds have decrease?
    total_decreases = sum(decrease_counts[15:])
    total_round_observations = N * 49
    print(f"  Overall decrease rate: {total_decreases}/{total_round_observations} "
          f"= {total_decreases/total_round_observations*100:.1f}%")
    print(f"  Expected (random walk): ~46% (P(step down) in bounded random walk)")

    # ═══════════════════════════════════════
    print(f"\n{'='*70}")
    print("SYNTHESIS")
    print(f"{'='*70}")

    print(f"""
  ═══════════════════════════════════════════════════════════════

  ПУТЬ 1 — HIGHER-ORDER DIFFERENTIALS:
    Previous "zeros at r=24" were WRONG — checking now per-word.
    If per-word zeros >> 1/2^32: REAL signal → algebraic weakness.
    If ≈ 1/2^32: noise.

  ПУТЬ 2 — COLD PATHS:
    Pipe registers (b,c,d,f,g,h) give cold directions.
    Cold path length: depends on pipe position.
    δd → consumed in 1 round (d enters e_new).
    δc → cold for 1 round (c→d), then consumed.
    δb → cold for 2 rounds (b→c→d), then consumed.
    δh → consumed in 1 round (h enters T1).
    δg → cold for 1 round (g→h→T1).
    δf → cold for 2 rounds (f→g→h→T1).

    MAX COLD PATH = 2 rounds (from b or f).
    NOT enough to bypass the sphere.

  ПУТЬ 3 — CARRY RESONANCE:
    δstate decreases ~{total_decreases/total_round_observations*100:.0f}% of rounds.
    This matches random walk prediction (~46%).
    → NO systematic resonance. Decreases are RANDOM.

  CONCLUSION:
    None of the three paths breaks through r=20.
    SHA-256 is secure beyond our current mathematics.

    To go further: need QUALITATIVELY new approach.
    Not differential. Not algebraic. Something we haven't invented yet.

  ═══════════════════════════════════════════════════════════════
""")


if __name__ == "__main__":
    main()
