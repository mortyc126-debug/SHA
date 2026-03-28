#!/usr/bin/env python3
"""
RAYON Per-Round — Axis 7.3
=============================

RAYON shows: constant propagation works when scheme has AND/OR chains.
SHA-256 round = h + Sig1(e) + Ch(e,f,g) + K + W.

Ch(e,f,g) = (e AND f) XOR (NOT e AND g) — has AND gates!
Modular addition = carry chains = AND chains!
Sig1 = XOR of rotations — LINEAR (no AND).

When δW[r]=0: round r gets same W for both messages.
The ONLY difference comes from state diff.
State diff propagates through Ch (AND-based) and addition (carry).

TEST: For a SINGLE round with known input state diff and W=0:
How many bits of output are DETERMINED by partial input assignment?
This is the RAYON "determination probability" per round.
"""

import numpy as np
import time, sys

MASK = 0xFFFFFFFF
H0 = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
      0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]
K = [0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
     0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
     0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
     0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
     0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
     0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
     0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
     0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2]

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&MASK
def Sig1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def Sig0(x): return rotr(x,2)^rotr(x,13)^rotr(x,22)
def Ch(e,f,g): return (e&f)^(~e&g)&MASK
def Maj(a,b,c): return (a&b)^(a&c)^(b&c)
def hw(x): return bin(x&MASK).count('1')

def one_round(a,b,c,d,e,f,g,h,W,r):
    """One SHA-256 round. Returns new (a,b,c,d,e,f,g,h)."""
    T1 = (h + Sig1(e) + Ch(e,f,g) + K[r] + W) & MASK
    T2 = (Sig0(a) + Maj(a,b,c)) & MASK
    return (T1+T2)&MASK, a, b, c, (d+T1)&MASK, e, f, g


def determination_test_round(n_bits, r, N_trials=500, seed=23000):
    """
    RAYON-style determination test for one SHA-256 round.

    Fix n_bits of 256-bit input state randomly.
    Check: is output bit determined?

    Measure: P(output determined after fixing k/256 state bits).
    """
    rng = np.random.RandomState(seed)

    # Single round: 256 input bits (a,b,c,d,e,f,g,h each 32)
    # + W[r] (32 bits) + K[r] (constant)
    # Output: 256 bits (new state)
    # We fix W[r] = known constant (no δW).
    # Vary state bits. Check output.

    total_input_bits = 256  # 8 words × 32 bits
    W_fixed = rng.randint(0, 1 << 32)  # arbitrary fixed W

    det_by_k = {}  # k → fraction determined

    for k in [32, 64, 96, 128, 160, 192, 224, 240, 248, 252, 254, 255]:
        n_det = 0
        for trial in range(N_trials):
            # Random base state
            state = [int(rng.randint(0, 1<<32)) for _ in range(8)]

            # Fix k bits (first k bits = fixed)
            # Free bits: bits k..255
            free_bits = total_input_bits - k

            if free_bits > 20:
                # Sample: check if ALL sampled completions give same output
                vals = set()
                n_samples = min(64, 1 << free_bits)
                for s in range(n_samples):
                    state_s = list(state)
                    # Flip random free bits
                    for fb in range(free_bits):
                        bit_pos = k + fb
                        word_idx = bit_pos // 32
                        bit_idx = bit_pos % 32
                        if rng.randint(0, 2):
                            state_s[word_idx] ^= (1 << bit_idx)
                    out = one_round(*state_s, W_fixed, r)
                    vals.add(out[4] & 1)  # check bit 0 of new_e
                    if len(vals) > 1:
                        break
                if len(vals) == 1:
                    n_det += 1
            else:
                # Exhaustive
                vals = set()
                for mask in range(1 << free_bits):
                    state_m = list(state)
                    for fb in range(free_bits):
                        bit_pos = k + fb
                        word_idx = bit_pos // 32
                        bit_idx = bit_pos % 32
                        if (mask >> fb) & 1:
                            state_m[word_idx] ^= (1 << bit_idx)
                    out = one_round(*state_m, W_fixed, r)
                    vals.add(out[4] & 1)
                    if len(vals) > 1:
                        break
                if len(vals) == 1:
                    n_det += 1

        det_by_k[k] = n_det / N_trials

    return det_by_k


def experiment_main(seed=23000):
    print("="*70)
    print("RAYON Per-Round DFS — Determination probability")
    print("="*70)

    # Test for specific rounds
    rounds_to_test = [0, 16, 17, 30, 47, 63]

    print(f"\n  --- P(output bit determined) after fixing k/256 state bits ---")
    print(f"  (Single round, fixed W, checking new_e[bit 0])")
    print()

    header = f"  {'k fixed':>8} |"
    for r in rounds_to_test:
        header += f" r={r:2d} |"
    print(header)
    print(f"  {'-'*8}-+" + "-+".join(['-'*6]*len(rounds_to_test)))

    all_results = {}
    for r in rounds_to_test:
        t0 = time.time()
        det = determination_test_round(256, r, N_trials=300, seed=seed+r)
        all_results[r] = det
        print(f"  ... r={r} done ({time.time()-t0:.1f}s)")

    # Print table
    for k in [128, 192, 224, 240, 248, 252, 254, 255]:
        line = f"  {k:4d}/256 |"
        for r in rounds_to_test:
            p = all_results[r].get(k, 0)
            marker = "*" if p > 0.6 else " "
            line += f" {p:.2f}{marker}|"
        print(line)

    # KEY: does the round constant K[r] affect determination?
    print(f"\n  --- Effect of K[r] on determination ---")
    print(f"  Hypothesis: large K[r] → more carries → more AND chains → higher determination")
    for r in rounds_to_test:
        k_ratio = K[r] / (1 << 32)
        det_at_240 = all_results[r].get(240, 0)
        print(f"    r={r:2d}: K/T={k_ratio:.3f}, P(det|240 fixed)={det_at_240:.3f}")

    # Ch analysis: how much does Ch help determination?
    print(f"\n  --- Ch(e,f,g) determination ---")
    print(f"  Ch = (e AND f) XOR (NOT e AND g)")
    print(f"  If e[bit]=0 → Ch[bit] = g[bit] (determined by g alone)")
    print(f"  If e[bit]=1 → Ch[bit] = f[bit] (determined by f alone)")
    print(f"  Ch is determined once e[bit] is known!")
    print(f"  → P(Ch determined | e known) = 1.0")
    print(f"  → 32 bits of Ch determined by 32 bits of e")
    print(f"  This is the AND-chain that RAYON exploits!")

    # Carry chain in addition
    print(f"\n  --- Carry chain in T1 = h + Sig1(e) + Ch(e,f,g) + K + W ---")
    print(f"  Addition of 5 values: carry from bit i → bit i+1")
    print(f"  If all inputs to bit i known → carry determined")
    print(f"  → LSB (bit 0) always determined (no carry in)")
    print(f"  → bit 1 determined if carry from bit 0 determined")
    print(f"  → cascade: determined from bit 0 upward")
    print(f"  RAYON uses this: fix lower bits first → propagate carries")

    return all_results


if __name__ == '__main__':
    print("RAYON Per-Round — Axis 7.3")
    print(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")

    results = experiment_main()

    print(f"""
{'='*70}
SYNTHESIS: RAYON Per-Round (Axis 7.3)
{'='*70}

Determination probability for SINGLE SHA-256 round:
  How many state bits must be fixed to determine output?

If P(det) high at k=240 (94% fixed) → RAYON DFS ε ≈ 0.06
If P(det) high at k=192 (75% fixed) → RAYON DFS ε ≈ 0.25
If P(det) high at k=128 (50% fixed) → RAYON DFS ε ≈ 0.50

RAYON SHA-1 result: ε=0.87 for one round (without schedule).
Our SHA-256: comparable or different?
""")
