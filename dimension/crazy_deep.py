"""
ГЛУБОКИЕ БЕЗУМСТВА.

1. SELF-REFERENTIAL VERIFICATION: 92 = сигнал или шум? 1M проверка.

2. SHA-256 × SHA-256: две копии скрещиваются.
   Compute H1 = SHA256(M). Use H1 to MODIFY M → M'.
   Compute H2 = SHA256(M'). Use H2 to modify M → M''.
   Iterate. Ищем СХОДИМОСТЬ к fixed cycle.

3. PERIOD OF ITERATION: H → H(H) → H(H(H)) → ...
   При 32-bit truncation: period ≈ 2^16 (birthday).
   Но: SHA-256 iteration может иметь КОРОТКИЙ period для special M.

4. BACKWARD HASH: вычисляем SHA-256 В ОБРАТНОМ ПОРЯДКЕ раундов.
   SHA256_reverse(M) = rounds 63,62,...,0 вместо 0,1,...,63.
   SHA256(M) = SHA256_reverse(M) → SYMMETRIC COLLISION.
"""

import numpy as np
import struct, hashlib
import time

MASK32 = 0xFFFFFFFF

def hw(x): return bin(x).count('1')
def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def add32(x, y): return (x + y) & MASK32

def sha256_full(W16):
    raw = struct.pack('>16I', *W16)
    return struct.unpack('>8I', hashlib.sha256(raw).digest())


def main():
    np.random.seed(42)

    print("=" * 70)
    print("ГЛУБОКИЕ БЕЗУМСТВА")
    print("=" * 70)

    # ═══════════════════
    print(f"\n{'='*70}")
    print("1. SELF-REFERENTIAL: 1M verification")
    print(f"{'='*70}")

    # M' = M ⊕ H(M). Is δH consistently < 128?
    N = 500000
    self_dHs = []
    random_dHs = []

    t0 = time.time()
    for _ in range(N):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H = sha256_full(W)

        # Self-referential
        W_self = [W[i] ^ H[i % 8] for i in range(16)]
        H_self = sha256_full(W_self)
        self_dHs.append(sum(hw(H[i] ^ H_self[i]) for i in range(8)))

        # Random pair (same HW of δW for fair comparison)
        W_rand = [np.random.randint(0, 2**32) for _ in range(16)]
        H_rand = sha256_full(W_rand)
        random_dHs.append(sum(hw(H[i] ^ H_rand[i]) for i in range(8)))

    elapsed = time.time() - t0
    print(f"\n  {N} pairs in {elapsed:.1f}s:")
    print(f"    Self-ref:  mean={np.mean(self_dHs):.3f}, std={np.std(self_dHs):.3f}, "
          f"min={min(self_dHs)}")
    print(f"    Random:    mean={np.mean(random_dHs):.3f}, std={np.std(random_dHs):.3f}, "
          f"min={min(random_dHs)}")

    diff = np.mean(self_dHs) - np.mean(random_dHs)
    from scipy import stats
    t, p = stats.ttest_ind(self_dHs, random_dHs)
    print(f"    Difference: {diff:.3f} bits")
    print(f"    t-test: t={t:.2f}, p={p:.4f}")
    print(f"    → {'★ SIGNAL!' if p < 0.001 else 'Noise'}")

    # HW of δW for self-ref
    self_dW_hw = []
    for _ in range(10000):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H = sha256_full(W)
        W_self = [W[i] ^ H[i % 8] for i in range(16)]
        self_dW_hw.append(sum(hw(W[i] ^ W_self[i]) for i in range(16)))
    print(f"\n    HW(δW) for self-ref: mean={np.mean(self_dW_hw):.1f}")
    print(f"    (high HW → more mixing → δH closer to 128)")

    # ═══════════════════
    print(f"\n{'='*70}")
    print("2. SHA-256 × SHA-256: cross-iteration")
    print(f"{'='*70}")

    # Start with M. Iterate:
    # M₁ = M ⊕ H(M)
    # M₂ = M₁ ⊕ H(M₁)
    # ...
    # Look for: H(Mₙ) = H(Mₘ) for some n ≠ m (cycle in H-space)

    # Track: HW(H(Mₙ) ⊕ H(Mₙ₋₁)) per step. Does it decrease?

    print(f"\n  Cross-iteration trajectories (10 starts × 1000 steps):")

    for start_idx in range(10):
        M = [np.random.randint(0, 2**32) for _ in range(16)]
        H_prev = sha256_full(M)

        min_dH = 256
        convergence_step = -1

        for step in range(1000):
            # M_next = M ⊕ H(M) (folded: H applied to message positions)
            M_next = [M[i] ^ H_prev[i % 8] for i in range(16)]
            H_curr = sha256_full(M_next)

            dH = sum(hw(H_prev[i] ^ H_curr[i]) for i in range(8))
            if dH < min_dH:
                min_dH = dH
                convergence_step = step

            M = M_next
            H_prev = H_curr

        print(f"    Start {start_idx}: min δH(consecutive) = {min_dH} at step {convergence_step}")

    # ═══════════════════
    print(f"\n{'='*70}")
    print("3. ITERATION PERIOD: H^n(M) cycle detection")
    print(f"{'='*70}")

    # Truncate to 16 bits for tractability
    # H_trunc(M) = SHA256(M)[0] & 0xFFFF

    print(f"\n  16-bit truncated iteration: period detection")

    for trial in range(5):
        M = [np.random.randint(0, 2**32) for _ in range(16)]

        # Floyd cycle detection
        # Tortoise and hare on truncated hash
        def next_hash(W):
            H = sha256_full(W)
            W_next = list(W)
            W_next[0] = H[0] & 0xFFFF  # truncate and feed back
            return W_next

        tortoise = next_hash(list(M))
        hare = next_hash(next_hash(list(M)))

        steps = 1
        max_steps = 200000
        while tortoise[0] != hare[0] and steps < max_steps:
            tortoise = next_hash(tortoise)
            hare = next_hash(next_hash(hare))
            steps += 1

        if steps < max_steps:
            # Find cycle length
            hare = next_hash(tortoise)
            cycle_len = 1
            while hare[0] != tortoise[0]:
                hare = next_hash(hare)
                cycle_len += 1
            print(f"    Trial {trial}: cycle found at step {steps}, length = {cycle_len}")
            print(f"    Expected (16-bit): ~{int(2**8)} (sqrt(2^16))")
        else:
            print(f"    Trial {trial}: no cycle in {max_steps} steps")

    # ═══════════════════
    print(f"\n{'='*70}")
    print("4. BACKWARD SHA-256: reverse round order")
    print(f"{'='*70}")

    # Normal: rounds 0→63 with K[0]→K[63] and W[0]→W[63]
    # Reverse: rounds 63→0 with K[63]→K[0] and W[63]→W[0]
    # If H_forward(M) = H_backward(M) → symmetric collision!

    def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
    def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
    def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
    def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
    def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
    def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)

    K = [0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
         0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
         0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
         0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
         0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
         0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
         0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
         0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2]
    IV = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
          0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]

    def sha256_reverse(W16):
        """SHA-256 with rounds in reverse order."""
        W = list(W16)
        for r in range(16, 64):
            W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))
        a,b,c,d,e,f,g,h = IV
        for r in range(63, -1, -1):  # REVERSE
            T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K[r]), W[r])
            T2 = add32(Sigma0(a), Maj(a,b,c))
            h,g,f,e = g,f,e,add32(d,T1)
            d,c,b,a = c,b,a,add32(T1,T2)
        return tuple(add32(IV[i], [a,b,c,d,e,f,g,h][i]) for i in range(8))

    # Search: H_forward(M) close to H_reverse(M)?
    N = 50000
    forward_reverse_dHs = []
    best_fr = 256

    for _ in range(N):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H_fwd = sha256_full(W)
        H_rev = sha256_reverse(W)
        dH = sum(hw(H_fwd[i] ^ H_rev[i]) for i in range(8))
        forward_reverse_dHs.append(dH)
        if dH < best_fr: best_fr = dH

    print(f"\n  H_forward(M) vs H_reverse(M) ({N} messages):")
    print(f"    Mean δH: {np.mean(forward_reverse_dHs):.1f}")
    print(f"    Min δH:  {best_fr}")
    print(f"    → {'★ STRUCTURE' if np.mean(forward_reverse_dHs) < 126 or best_fr < 85 else 'Random'}")

    # ═══════════════════
    print(f"\n{'='*70}")
    print("5. HASH COMPOSITION: H(M) ⊕ H(reverse(M))")
    print(f"{'='*70}")

    # SHA256(M) ⊕ SHA256(reverse_words(M))
    # If this has LOW HW → words order doesn't matter much → structure
    compose_hws = []
    for _ in range(N):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        W_rev = list(reversed(W))
        H1 = sha256_full(W); H2 = sha256_full(W_rev)
        compose_hws.append(sum(hw(H1[i]^H2[i]) for i in range(8)))

    print(f"  H(M) ⊕ H(rev(M)): mean={np.mean(compose_hws):.1f}, min={min(compose_hws)}")

    # H(M) ⊕ H(M+1) where M+1 = increment W[0]
    inc_hws = []
    for _ in range(N):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W); W2[0] = (W2[0] + 1) & MASK32
        H1 = sha256_full(W); H2 = sha256_full(W2)
        inc_hws.append(sum(hw(H1[i]^H2[i]) for i in range(8)))

    print(f"  H(M) ⊕ H(M+1): mean={np.mean(inc_hws):.1f}, min={min(inc_hws)}")

    # ═══════════════════
    print(f"\n{'='*70}")
    print("ИТОГ")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
