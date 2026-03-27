"""
НОВЫЕ УГЛЫ: то что мы НЕ пробовали.

Все предыдущие подходы: смотрим на SHA-256 КАК ЕСТЬ.
Новый угол: смотрим на SHA-256 КАК ПРОЦЕСС.

Идея 1: TEMPORAL ATTACK — не одна пара (M1,M2),
а ПОСЛЕДОВАТЕЛЬНОСТЬ M1→M2→M3→... где каждый следующий
ИСПОЛЬЗУЕТ информацию от предыдущего.

Идея 2: COLLISION BETWEEN DIFFERENT ROUND COUNTS
  H_r1(M1) = H_r2(M2) для r1 ≠ r2.
  Мы всегда искали collision при ОДНОМ r.
  Что если при РАЗНЫХ раундах легче?

Идея 3: PARTIAL STATE COLLISION
  Не H1=H2, а H1[0..3]=H2[0..3] (половина hash).
  Feedforward: H[i] = IV[i] + state[i].
  Collision в state[0..3] = collision в H[0..3].
  Pipe structure: state[0..3] = (a,b,c,d), state[4..7] = (e,f,g,h).
  Две НЕЗАВИСИМЫЕ pipe chains. Может одна слабее?

Идея 4: INTERNAL STATE COLLISION (free-start)
  Не H(M1)=H(M2), а state_r(M1)=state_r(M2) для промежуточного r.
  Feedforward потом РАЗРУШИТ collision, но это distinguisher.

Идея 5: CROSS-FUNCTION
  SHA-256 vs SHA-224 (same compression, different IV, truncated).
  Найти M такое что SHA256(M)[0..6] = SHA224(M).
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


def sha256_full(W16):
    raw = struct.pack('>16I', *W16)
    return struct.unpack('>8I', hashlib.sha256(raw).digest())


def compress_state(iv, W16, n_r):
    """Return INTERNAL state (before feedforward)."""
    W = list(W16)
    for r in range(16, max(n_r, 16)):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))
    a,b,c,d,e,f,g,h = iv
    for r in range(n_r):
        T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K[r]), W[r])
        T2 = add32(Sigma0(a), Maj(a,b,c))
        h,g,f,e = g,f,e,add32(d,T1)
        d,c,b,a = c,b,a,add32(T1,T2)
    return (a,b,c,d,e,f,g,h)


def main():
    np.random.seed(42)

    print("=" * 70)
    print("НОВЫЕ УГЛЫ АТАКИ")
    print("=" * 70)

    # ═══════════════════
    print(f"\n{'='*70}")
    print("УГОЛ 1: PARTIAL COLLISION — только 4 из 8 слов")
    print(f"{'='*70}")

    # H[0..3] = first half = (IV[0]+a, IV[1]+b, IV[2]+c, IV[3]+d)
    # H[4..7] = second half = (IV[4]+e, IV[5]+f, IV[6]+g, IV[7]+h)
    # Pipes: a→b→c→d (first chain), e→f→g→h (second chain)
    # Are they EQUALLY hard to collide?

    N = 50000
    for n_r in [20, 24, 32, 64]:
        best_first_half = 128  # 4 words × 32 bits
        best_second_half = 128
        best_full = 256

        for _ in range(N):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            W2 = list(W); W2[15] ^= (1 << 31)

            state1 = compress_state(IV, W, n_r)
            state2 = compress_state(IV, W2, n_r)

            # First half collision (a,b,c,d match)
            d_first = sum(hw((add32(IV[i], state1[i])) ^ (add32(IV[i], state2[i]))) for i in range(4))
            d_second = sum(hw((add32(IV[i], state1[i])) ^ (add32(IV[i], state2[i]))) for i in range(4, 8))
            d_full = d_first + d_second

            if d_first < best_first_half: best_first_half = d_first
            if d_second < best_second_half: best_second_half = d_second
            if d_full < best_full: best_full = d_full

        print(f"  r={n_r:>2}: first_half={best_first_half:>3}/128, "
              f"second_half={best_second_half:>3}/128, full={best_full:>3}/256")

    # ═══════════════════
    print(f"\n{'='*70}")
    print("УГОЛ 2: INTERNAL STATE COLLISION (without feedforward)")
    print(f"{'='*70}")

    # state_r(M1) = state_r(M2) at intermediate round r
    # This is EASIER than output collision because no feedforward

    for n_r in [20, 24, 32, 64]:
        best_state = 256
        best_output = 256

        for _ in range(N):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            W2 = list(W); W2[15] ^= (1 << 31)

            s1 = compress_state(IV, W, n_r)
            s2 = compress_state(IV, W2, n_r)

            # State collision (no feedforward)
            d_state = sum(hw(s1[i] ^ s2[i]) for i in range(8))
            # Output collision (with feedforward)
            d_output = sum(hw(add32(IV[i],s1[i]) ^ add32(IV[i],s2[i])) for i in range(8))

            if d_state < best_state: best_state = d_state
            if d_output < best_output: best_output = d_output

        diff = best_output - best_state
        print(f"  r={n_r:>2}: state={best_state:>3}, output={best_output:>3}, "
              f"feedforward cost={diff:>+3}")

    # ═══════════════════
    print(f"\n{'='*70}")
    print("УГОЛ 3: ASYMMETRIC δW — exploit pipe chain asymmetry")
    print(f"{'='*70}")

    # a→b→c→d chain: d enters e_new = d + T1
    # e→f→g→h chain: h enters T1 = h + Σ1(e) + Ch(e,f,g) + K + W
    # h is CONSUMED differently than d!
    # d: added to T1 → goes to e
    # h: added at START of T1 computation → through Σ1, Ch, K, W

    # Question: is flipping a bit in one chain DIFFERENT from the other?
    print(f"\n  δ in a-chain vs e-chain registers (full SHA-256, 50K):")

    for target_reg in range(8):
        reg_names = ['a','b','c','d','e','f','g','h']
        best = 256
        for _ in range(N):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            # Modify initial state at specific register by tweaking M
            # Can't directly modify IV, so use a proxy:
            # Flip bit in W[target_reg % 16] as proxy
            W2 = list(W); W2[target_reg] ^= (1 << 31)
            H1 = sha256_full(W); H2 = sha256_full(W2)
            d = sum(hw(H1[i]^H2[i]) for i in range(8))
            if d < best: best = d
        chain = "a-chain" if target_reg < 4 else "e-chain"
        print(f"    δW[{target_reg}] ({reg_names[target_reg]}, {chain}): best={best:>3}")

    # ═══════════════════
    print(f"\n{'='*70}")
    print("УГОЛ 4: ADAPTIVE SEARCH — use H to guide next step")
    print(f"{'='*70}")

    # Instead of random search: start from (M1, M2=M1⊕δ),
    # look at H1⊕H2, and use the PATTERN to choose next δ.

    # Hill climbing: try to MINIMIZE δH iteratively.
    for n_r in [24, 32, 64]:
        best_random = 256
        best_adaptive = 256

        # Random search baseline
        for _ in range(N):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            W2 = list(W); W2[15] ^= (1 << 31)
            H1 = sha256_full(W) if n_r == 64 else tuple(compress_state(IV, W, n_r))
            H2 = sha256_full(W2) if n_r == 64 else tuple(compress_state(IV, W2, n_r))
            d = sum(hw(H1[i]^H2[i]) for i in range(8))
            if d < best_random: best_random = d

        # Adaptive: start with random, then mutate MESSAGE to reduce δH
        for start in range(500):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            dW = 1 << 31  # fixed δ in W[15]

            for step in range(100):
                W2 = list(W); W2[15] ^= dW
                H1 = sha256_full(W) if n_r == 64 else tuple(add32(IV[i], compress_state(IV, W, n_r)[i]) for i in range(8))
                H2 = sha256_full(W2) if n_r == 64 else tuple(add32(IV[i], compress_state(IV, W2, n_r)[i]) for i in range(8))
                current_d = sum(hw(H1[i]^H2[i]) for i in range(8))

                if current_d < best_adaptive:
                    best_adaptive = current_d

                # Mutate: flip random bit in W[0..14] (NOT W[15])
                w = np.random.randint(0, 15)
                b = np.random.randint(0, 32)
                W_new = list(W); W_new[w] ^= (1 << b)
                W2_new = list(W_new); W2_new[15] ^= dW

                H1n = sha256_full(W_new) if n_r == 64 else tuple(add32(IV[i], compress_state(IV, W_new, n_r)[i]) for i in range(8))
                H2n = sha256_full(W2_new) if n_r == 64 else tuple(add32(IV[i], compress_state(IV, W2_new, n_r)[i]) for i in range(8))
                new_d = sum(hw(H1n[i]^H2n[i]) for i in range(8))

                if new_d <= current_d:
                    W = W_new  # accept improvement or lateral move

        print(f"  r={n_r:>2}: random={best_random:>3} (50K), adaptive={best_adaptive:>3} (500×100)")

    # ═══════════════════
    print(f"\n{'='*70}")
    print("УГОЛ 5: BIRTHDAY на partial match")
    print(f"{'='*70}")

    # Standard birthday: find M1,M2 with H(M1)=H(M2). Cost: 2^128.
    # Partial birthday: find M1,M2 with H(M1)[0]=H(M2)[0]. Cost: 2^16.
    # Can we find FULL collision by EXTENDING partial match?

    # Step 1: find pair with H[0] matching (32-bit birthday, ~65K hashes)
    hash_table = {}
    collisions_found = 0

    for trial in range(200000):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H = sha256_full(W)

        key = H[0]  # first word only
        if key in hash_table:
            W_prev = hash_table[key]
            if W != W_prev:
                H_prev = sha256_full(W_prev)
                # Check: how many OTHER words match?
                full_d = sum(hw(H[i]^H_prev[i]) for i in range(8))
                words_match = sum(1 for i in range(8) if H[i] == H_prev[i])

                if collisions_found < 5:
                    print(f"  H[0] collision #{collisions_found+1}: "
                          f"full δH={full_d}, matching words={words_match}/8")
                collisions_found += 1
        else:
            hash_table[key] = W

    print(f"\n  Total H[0] collisions in 200K messages: {collisions_found}")
    print(f"  Expected (birthday on 32 bits): ~{200000**2 // (2**33):,}")

    # ═══════════════════
    print(f"\n{'='*70}")
    print("ИТОГ")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
