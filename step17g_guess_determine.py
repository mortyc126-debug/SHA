#!/usr/bin/env python3
"""Step 17g: Guess-and-determine collision search for reduced SHA-256.

The NB approach works BACKWARDS from conditions:
1. Fix a differential path target
2. Compute required conditions on state bits
3. Choose message words to satisfy conditions

Key: for each round, W[t] directly determines a[t] and e[t].
Given the previous state, W[t] is DETERMINED by the target a[t] or e[t].

Strategy for 9-round collision (local collision at steps 0-8):
- DW[0]=+1, DW[1]=-1, DW[8]=correction
- Target: Da=+1, De=+1 after step 0; Da=0, De=0 after step 1; stay 0
- After step 0: diffs in (a,e) = (+1,+1)
- After step 1: diffs in (b,f) = (+1,+1), rest 0 if corrected
- After step 2: diffs in (c,g) = (+1,+1) from shift register
- ...
- After step 4: diffs exit d,h
- After step 8: all diffs gone

BUT: Da_new - De_new = T2_diff - d_diff is FIXED. If we want Da=0, De=0
we need T2_diff = d_diff. We showed P(this) = 0.

ALTERNATIVE TARGET: Don't require Da=0 at every step. Instead, let Da
follow its natural value, and only control De. The Da values shift through
b,c,d and exit. After 4 rounds of zero Da, all will be clear.

Better approach: use DW at steps 0-3 to control De (force to 0).
After step 3: De=0, but Db=Da0, Dc=Da1, Dd=Da2, Df=0, Dg=0, Dh=0.
Da3 is whatever it is.

After step 4 (DW=0): De is determined by Dd + DT1 = Da2 + DT1.
For De4=0: DT1 = -Da2 (need specific condition).

This is still probabilistic. Let me try a completely different approach:
COMPUTE W[0..15] that satisfies the collision condition directly.

For N-round collision: need hash(IV, W) = hash(IV, W')
where DW has a specific pattern.

For 9 rounds: output = IV + state_after_9_rounds.
Need: state1 = state2 (after 9 rounds with W vs W').

This gives 256-bit constraint with ~480 bits freedom.
Should have solutions. But finding them requires solving a nonlinear system.

Let's try: fix W[0] and DW pattern, then use Newton-like iteration
on W[1..15] to minimize the hash difference.
"""

import random
import math
MASK32 = 0xFFFFFFFF

def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def shr(x, n): return x >> n
def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ shr(x, 3)
def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ shr(x, 10)
def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
def hw(x): return bin(x).count('1')

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

def sha256_step(state, W_t, K_t):
    a, b, c, d, e, f, g, h = state
    T1 = (h + Sigma1(e) + Ch(e, f, g) + K_t + W_t) & MASK32
    T2 = (Sigma0(a) + Maj(a, b, c)) & MASK32
    return ((T1 + T2) & MASK32, a, b, c, (d + T1) & MASK32, e, f, g)

def expand_message(W16):
    W = list(W16)
    for t in range(16, 64):
        W.append((sigma1(W[t-2]) + W[t-7] + sigma0(W[t-15]) + W[t-16]) & MASK32)
    return W

def hash_n(W16, n):
    """Compute n-round reduced SHA-256 hash."""
    W = expand_message(W16)
    s = list(IV)
    for t in range(n):
        s = list(sha256_step(s, W[t], K[t]))
    return tuple((IV[i] + s[i]) & MASK32 for i in range(8))

def hash_diff_hw(W, W2, n):
    h1 = hash_n(W, n)
    h2 = hash_n(W2, n)
    return sum(hw(h1[i] ^ h2[i]) for i in range(8))

# ============================================================
# Approach: Gradient-like bit-flipping optimization
# ============================================================
def bitflip_optimization():
    """Start from random message, flip bits to minimize hash diff."""
    print("=" * 70)
    print("Gradient-like bit-flipping optimization")
    print("=" * 70)

    for n_rounds in [9, 13, 17, 21]:
        print(f"\n  === {n_rounds}-round collision search ===")

        best_overall = 999
        best_W = None

        for restart in range(100):
            W = [random.randint(0, MASK32) for _ in range(16)]
            W2 = list(W)
            W2[0] = (W[0] + 1) & MASK32

            current_hw = hash_diff_hw(W, W2, n_rounds)

            # Greedy bit-flipping: flip each bit of W[1..15] and keep if improves
            improved = True
            while improved:
                improved = False
                for word in range(1, 16):
                    for bit in range(32):
                        W_test = list(W)
                        W_test[word] ^= (1 << bit)
                        W2_test = list(W_test)
                        W2_test[0] = (W_test[0] + 1) & MASK32

                        test_hw = hash_diff_hw(W_test, W2_test, n_rounds)
                        if test_hw < current_hw:
                            W = W_test
                            W2 = W2_test
                            current_hw = test_hw
                            improved = True
                            if current_hw == 0:
                                break
                    if current_hw == 0:
                        break

            if current_hw < best_overall:
                best_overall = current_hw
                best_W = list(W)

            if current_hw == 0:
                print(f"    COLLISION at restart {restart}!")
                print(f"    W  = {[f'0x{w:08x}' for w in W]}")
                h = hash_n(W, n_rounds)
                print(f"    H  = {[f'0x{x:08x}' for x in h[:4]]}...")
                break

            if restart % 20 == 0:
                print(f"    restart {restart}: current best = {best_overall}")

        print(f"    Best after 100 restarts: HW = {best_overall}")

# ============================================================
# Approach 2: Stochastic hill climbing with temperature
# ============================================================
def stochastic_search():
    """Simulated annealing style search."""
    print("\n" + "=" * 70)
    print("Stochastic hill climbing")
    print("=" * 70)

    for n_rounds in [9, 13]:
        print(f"\n  === {n_rounds}-round collision search ===")

        best_overall = 999

        for restart in range(20):
            W = [random.randint(0, MASK32) for _ in range(16)]
            W2 = list(W)
            W2[0] = (W[0] + 1) & MASK32
            current_hw = hash_diff_hw(W, W2, n_rounds)

            temp = 10.0
            for step in range(50000):
                # Random perturbation: flip 1-3 bits
                W_test = list(W)
                n_flips = random.choice([1, 1, 1, 2, 2, 3])
                for _ in range(n_flips):
                    word = random.randint(1, 15)
                    bit = random.randint(0, 31)
                    W_test[word] ^= (1 << bit)

                W2_test = list(W_test)
                W2_test[0] = (W_test[0] + 1) & MASK32
                test_hw = hash_diff_hw(W_test, W2_test, n_rounds)

                # Accept if better, or with probability exp(-delta/temp)
                delta = test_hw - current_hw
                if delta < 0 or random.random() < math.exp(-delta / max(temp, 0.01)):
                    W = W_test
                    W2 = W2_test
                    current_hw = test_hw

                temp *= 0.9999

                if current_hw == 0:
                    print(f"    COLLISION at restart {restart}, step {step}!")
                    print(f"    W  = {[f'0x{w:08x}' for w in W[:8]]}...")
                    break

            if current_hw < best_overall:
                best_overall = current_hw

            if current_hw == 0:
                break

        print(f"    Best: HW = {best_overall}")

# ============================================================
# Approach 3: Meet-in-the-middle on internal state
# ============================================================
def mitm_search():
    """Split the computation at the midpoint and match."""
    print("\n" + "=" * 70)
    print("Meet-in-the-middle on internal state")
    print("=" * 70)

    # For 9-round collision: split at round 4
    # Forward from IV with W[0..3] → state at round 4
    # Backward from target (IV) with W[8..4] → state at round 4
    # If they match, collision found

    # But backward computation requires inverting the round function,
    # which is complex. Let's try a simpler approach:

    # Fix W[5..15] and vary W[1..4].
    # For each W[1..4], compute forward hash and store.
    # Look for collisions (same hash for different W[1..4] with DW[0]=+1).

    n_rounds = 9
    print(f"\n  === {n_rounds}-round partial MITM ===")

    # Fix most words, vary W[1..3]
    W_base = [random.randint(0, MASK32) for _ in range(16)]

    # Store: hash → W[1..3]
    table = {}
    N = 2**20  # 1M
    best_hw = 999

    for trial in range(N):
        W = list(W_base)
        # Randomize W[1], W[2], W[3]
        W[1] = random.randint(0, MASK32)
        W[2] = random.randint(0, MASK32)
        W[3] = random.randint(0, MASK32)

        W2 = list(W)
        W2[0] = (W[0] + 1) & MASK32

        h1 = hash_n(W, n_rounds)
        h2 = hash_n(W2, n_rounds)

        thw = sum(hw(h1[i] ^ h2[i]) for i in range(8))
        if thw < best_hw:
            best_hw = thw
            if thw == 0:
                print(f"    COLLISION at trial {trial}!")
                break

    print(f"    Best after {N} trials (varying W[1..3]): HW = {best_hw}")

    # Also try: fix everything, vary ONE word at a time
    for vary_word in [1, 5, 10, 14]:
        best_hw = 999
        for trial in range(N):
            W = list(W_base)
            W[vary_word] = random.randint(0, MASK32)
            W2 = list(W)
            W2[0] = (W[0] + 1) & MASK32

            h1 = hash_n(W, n_rounds)
            h2 = hash_n(W2, n_rounds)

            thw = sum(hw(h1[i] ^ h2[i]) for i in range(8))
            if thw < best_hw:
                best_hw = thw
                if thw == 0:
                    print(f"    COLLISION varying W[{vary_word}] at trial {trial}!")
                    break
        print(f"    Varying W[{vary_word}] only: best HW = {best_hw}")


if __name__ == "__main__":
    random.seed(42)

    bitflip_optimization()
    stochastic_search()
    mitm_search()

    print("\n" + "=" * 70)
    print("Step 17g Complete")
    print("=" * 70)
