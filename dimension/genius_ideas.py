"""
ГЕНИАЛЬНЫЕ, БЕЗУМНЫЕ, ЛОГИЧЕСКИЕ.

Все 27 подходов смотрели на SHA-256 СНАРУЖИ.
Новая идея: залезть ВНУТРЬ и СЛОМАТЬ изнутри.

БЕЗУМНАЯ 1: SELF-REFERENTIAL ATTACK
  Найти M такое что H(M) СОДЕРЖИТ инструкции для нахождения M'.
  H(M) интерпретируется как δW → M' = M ⊕ decode(H(M)).
  Если H(M) = H(M') → SELF-REFERENTIAL COLLISION.

БЕЗУМНАЯ 2: FIXED POINT
  Найти M такое что H(M) = M (первые 256 бит = hash).
  Или: H(M) = M[0..7] (hash = first 8 words of input).
  Fixed point = collision с САМИМ СОБОЙ в extended sense.

БЕЗУМНАЯ 3: BIRTHDAY МЕЖДУ РАЗНЫМИ ФУНКЦИЯМИ
  SHA256(M1) = SHA256(M2) трудно (2^128).
  Но: SHA256(M) = SHA256(rev(M))? Или SHA256(M) = SHA256(~M)?
  Structured pairs — не random birthday, а ALGEBRAIC relation.

ЛОГИЧЕСКАЯ 1: MEET-IN-THE-MIDDLE на schedule
  Schedule: W[0..15] → W[0..63].
  Разделить: W[0..7] → forward half, W[8..15] → backward half.
  Они ПЕРЕСЕКАЮТСЯ в W[16..23].
  MitM: forward вычисляет state[0..32] от W[0..7],
        backward вычисляет state[32..64] от W[8..15],
        match в state[32].

ЛОГИЧЕСКАЯ 2: EXPLOIT ADD OVERFLOW
  ADD mod 2^32: overflow = carry out of bit 31.
  Overflow = lost information. Each ADD loses ~0.5 bits.
  64 rounds × 7 ADDs = ~448 ADDs → ~224 bits lost.
  224 > 128 → enough information loss for collision?

ГЕНИАЛЬНАЯ: HASH AS EQUATION SYSTEM
  SHA-256(M) = H → 256 equations in 512 unknowns.
  Collision: SHA-256(M1) = SHA-256(M2) → 256 equations in 1024 unknowns.
  Solution space: 1024 - 256 = 768 dimensional.
  Can we SOLVE this system? Not brute force — ALGEBRA.
  Linearize: approximate as linear system, solve, use solution
  as starting point for nonlinear correction.
"""

import numpy as np
import struct, hashlib
import time

MASK32 = 0xFFFFFFFF

def hw(x): return bin(x).count('1')

def sha256_full(W16):
    raw = struct.pack('>16I', *W16)
    return struct.unpack('>8I', hashlib.sha256(raw).digest())

def sha256_r(W16, n_r):
    def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
    def add32(x, y): return (x + y) & MASK32
    def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
    def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)
    def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
    def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
    def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
    def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
    K = [0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
         0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
         0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
         0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
         0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
         0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
         0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
         0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2]
    IV = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]
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
    print("ГЕНИАЛЬНЫЕ, БЕЗУМНЫЕ, ЛОГИЧЕСКИЕ ИДЕИ")
    print("=" * 70)

    # ═══════════════════
    print(f"\n{'='*70}")
    print("БЕЗУМНАЯ 1: SELF-REFERENTIAL COLLISION")
    print("  M' = M ⊕ decode(H(M)). Ищем H(M) = H(M').")
    print(f"{'='*70}")

    # decode(H): use first 512 bits of hash as δW
    # (need 512 bits = 2 hashes? or truncate to modify only some words)
    # Simple: δW[i] = H[i % 8] for i=0..15

    best_self = 256
    for _ in range(100000):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H = sha256_full(W)

        # M' = M ⊕ H (H interpreted as δW)
        W2 = [W[i] ^ H[i % 8] for i in range(16)]

        # Is M' ≠ M? (H must be non-zero — always true)
        if W == W2:
            continue

        H2 = sha256_full(W2)
        dH = sum(hw(H[i] ^ H2[i]) for i in range(8))
        if dH < best_self:
            best_self = dH

    print(f"  100K self-referential pairs: best δH = {best_self}")
    print(f"  (128 = random, <128 = structure)")

    # ═══════════════════
    print(f"\n{'='*70}")
    print("БЕЗУМНАЯ 2: FIXED POINT H(M) = M[0..7]")
    print(f"{'='*70}")

    # Search: random M, check if H(M) ≈ M[0..7]
    best_fp = 256
    for _ in range(100000):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H = sha256_full(W)
        # Distance: HW(H[i] ⊕ W[i]) for i=0..7
        d = sum(hw(H[i] ^ W[i]) for i in range(8))
        if d < best_fp:
            best_fp = d

    print(f"  100K random: min HW(H(M) ⊕ M[0..7]) = {best_fp}")
    print(f"  Expected min (random): ~93")
    print(f"  → {'STRUCTURE' if best_fp < 85 else 'Random (no fixed point proximity)'}")

    # ═══════════════════
    print(f"\n{'='*70}")
    print("БЕЗУМНАЯ 3: STRUCTURED PAIRS")
    print("  SHA256(M) vs SHA256(transform(M))")
    print(f"{'='*70}")

    transforms = {
        'reverse words': lambda W: list(reversed(W)),
        'complement': lambda W: [w ^ MASK32 for w in W],
        'byte-swap': lambda W: [((w>>24)&0xFF)|((w>>8)&0xFF00)|((w<<8)&0xFF0000)|((w<<24)&0xFF000000) for w in W],
        'rotate all <<1': lambda W: [((w<<1)|(w>>31)) & MASK32 for w in W],
        'add 1 to all': lambda W: [(w+1) & MASK32 for w in W],
        'W[i] = W[15-i]': lambda W: [W[15-i] for i in range(16)],
        'W[i] ^= W[i+8%16]': lambda W: [W[i] ^ W[(i+8)%16] for i in range(16)],
        'negate all': lambda W: [(~w+1) & MASK32 for w in W],
    }

    N = 50000
    print(f"\n  {N} trials each:")
    for name, tf in transforms.items():
        best = 256
        for _ in range(N):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            W2 = tf(W)
            if W == W2: continue
            H1 = sha256_full(W); H2 = sha256_full(W2)
            d = sum(hw(H1[i]^H2[i]) for i in range(8))
            if d < best: best = d
        print(f"    {name:<25}: best δH = {best:>3} {'★' if best < 88 else ''}")

    # ═══════════════════
    print(f"\n{'='*70}")
    print("ЛОГИЧЕСКАЯ 1: ADD OVERFLOW = INFORMATION LOSS")
    print(f"{'='*70}")

    # Each ADD mod 2^32: carry out of bit 31 = LOST.
    # P(overflow) = P(sum ≥ 2^32) ≈ probability that carry=1.
    # For random operands: P(carry out) ≈ 50%.
    # Info lost per ADD: ~0.5 bits (entropy of carry out).
    # 64 rounds × 7 ADDs = 448 ADDs → ~224 bits lost.
    # Input: 512 bits. Output: 256 bits. Lost: 256 bits.
    # Overflow accounts for 224/256 = 87.5% of information loss!

    # Verify: count actual overflows
    overflow_counts = []
    N = 5000
    for _ in range(N):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        count = 0
        # Simulate SHA-256 counting overflows
        def add32_count(x, y):
            nonlocal count
            s = x + y
            if s >= (1 << 32):
                count += 1
            return s & MASK32

        def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
        def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
        def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
        def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
        def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
        def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
        def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)

        K = [0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
             0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174]
        IV = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]

        Ww = list(W)
        for r in range(16, 64):
            Ww.append(add32_count(add32_count(add32_count(sigma1(Ww[r-2]), Ww[r-7]), sigma0(Ww[r-15])), Ww[r-16]))
        a,b,c,d,e,f,g,h = IV
        for r in range(64):
            T1 = add32_count(add32_count(add32_count(add32_count(h, Sigma1(e)), Ch(e,f,g)), K[r%16]), Ww[r])
            T2 = add32_count(Sigma0(a), Maj(a,b,c))
            h,g,f,e = g,f,e,add32_count(d,T1)
            d,c,b,a = c,b,a,add32_count(T1,T2)
        overflow_counts.append(count)

    total_adds = 64 * 7 + 48 * 3  # round ADDs + schedule ADDs
    print(f"  Total ADDs per hash: ~{total_adds}")
    print(f"  Overflows per hash: mean={np.mean(overflow_counts):.0f}, std={np.std(overflow_counts):.0f}")
    print(f"  Overflow rate: {np.mean(overflow_counts)/total_adds*100:.1f}%")
    print(f"  Info lost: ~{np.mean(overflow_counts)*0.5:.0f} bits per hash")
    print(f"  (Input 512 → Output 256 = 256 bits lost)")
    print(f"  Overflow accounts for {np.mean(overflow_counts)*0.5/256*100:.0f}% of info loss")

    # ═══════════════════
    print(f"\n{'='*70}")
    print("ГЕНИАЛЬНАЯ: LINEARIZED SYSTEM + CORRECTION")
    print(f"{'='*70}")

    # SHA-256 ≈ Linear(M) + Noise(M)
    # Linear part: GF(2) Jacobian T (we computed this — rank 256)
    # Noise: CE matrix (carry errors)
    # Solve: T × δM = 0 (GF(2) kernel, dim=256)
    # Any δM in kernel: T × δM = 0 → δH ≈ CE × δM (carry errors only)
    # If CE × δM ≈ 0: near-collision!

    # We already know: CE is full rank (256) and CE tail = random.
    # BUT: what about the LINEAR system T?
    # T × δM = 0 has 2^256 solutions.
    # For each: actual δH = CE(δM) (nonlinear function of δM).
    # Can we find δM in kernel where CE(δM) is SMALL?

    # This is EXACTLY our GF(2)-kernel search from earlier.
    # Result: CE(kernel vectors) has same distribution as random δH.
    # Already verified — no advantage.

    # NEW TWIST: what about ITERATIVE linearization?
    # Round 1: linearize around current M, find δM.
    # Round 2: M' = M ⊕ δM, linearize around M', find δM'.
    # Repeat: converge to collision?

    print(f"\n  Iterative linearization (Newton-like method):")
    print(f"  Start: random M, target: H(M) = H(M ⊕ δM)")

    successes = 0
    for start in range(100):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H_target = sha256_full(W)

        # Iterative: try to find δM that makes H(M⊕δM) closer to H(M)
        current_dH = 256
        best_iter_dH = 256

        for iteration in range(20):
            # Random δM (can't compute true Jacobian inverse easily)
            # Use random search in GF(2) neighborhood
            best_local = current_dH
            for _ in range(500):
                dM = [0]*16
                # Flip 1-3 random bits
                for _ in range(np.random.randint(1, 4)):
                    w = np.random.randint(0, 16)
                    b = np.random.randint(0, 32)
                    dM[w] ^= (1 << b)

                W2 = [W[i] ^ dM[i] for i in range(16)]
                H2 = sha256_full(W2)
                dH = sum(hw(H_target[i] ^ H2[i]) for i in range(8))

                if dH < best_local:
                    best_local = dH
                    if dH < best_iter_dH:
                        best_iter_dH = dH
                        W = W2  # accept improvement

        if best_iter_dH < 100:
            successes += 1

    print(f"  100 starts × 20 iterations × 500 probes = 1M evals")
    print(f"  Best δH reached: varies by start")
    print(f"  Starts reaching δH < 100: {successes}/100")

    # Compare with pure random 1M evals
    best_pure = 256
    for _ in range(1000000):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = [W[i] ^ (np.random.randint(0, 2**32) if np.random.random() < 0.1 else 0) for i in range(16)]
        if W == W2: continue
        H1 = sha256_full(W); H2 = sha256_full(W2)
        dH = sum(hw(H1[i]^H2[i]) for i in range(8))
        if dH < best_pure: best_pure = dH

    print(f"  Pure random 1M: best δH = {best_pure}")

    # ═══════════════════
    print(f"\n{'='*70}")
    print("ИТОГ")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
