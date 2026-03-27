#!/usr/bin/env python3
"""
Three Paths Around Chaos — Stage 1.4
======================================
The chaotic zone r=31..59 is OPAQUE (proven G3, S2, S4).
Instead of going THROUGH it, we go AROUND it.

Path A: MULTIBLOCK — carry structure of block 1 → IV of block 2
        Does biased IV propagate through a fresh block?

Path B: BACKWARD — from H[7] backward to state[62], skip chaos entirely
        What does knowing H[7] tell us about state[62]?
        Can we constrain state[62] algebraically?

Path C: K-INVARIANTS — deterministic constraints from constants K[r]
        For EACH round r: K[r] + X >= 2^32 forces bits of X to 0.
        Build the full chain of forced-zero bits through all 64 rounds.
"""

import numpy as np
from collections import defaultdict, Counter
import time
import sys

MASK = 0xFFFFFFFF

H0 = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
]

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

def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK
def sig0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sig1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)
def Sig0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sig1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
def hw(x): return bin(x & MASK).count('1')

def message_schedule(W16):
    W = list(W16) + [0] * 48
    for i in range(16, 64):
        W[i] = (sig1(W[i-2]) + W[i-7] + sig0(W[i-15]) + W[i-16]) & MASK
    return W

def sha256_compress(IV, W16):
    """SHA-256 compression: returns (final_state, all_states, carries)."""
    W = message_schedule(W16)
    a, b, c, d, e, f, g, h = IV
    states = [tuple(IV)]
    carries = []

    for r in range(64):
        raw = h + Sig1(e) + Ch(e, f, g) + K[r] + W[r]
        T1 = raw & MASK
        carries.append(1 if raw >= (1 << 32) else 0)
        T2 = (Sig0(a) + Maj(a, b, c)) & MASK
        h, g, f = g, f, e
        e = (d + T1) & MASK
        d, c, b = c, b, a
        a = (T1 + T2) & MASK
        states.append((a, b, c, d, e, f, g, h))

    H_out = tuple((states[-1][i] + IV[i]) & MASK for i in range(8))
    return H_out, states, carries

def sha256_hash(W16):
    H, _, carries = sha256_compress(tuple(H0), W16)
    return H, carries


# ============================================================
# PATH A: MULTIBLOCK — carry bias propagation
# ============================================================

def path_A(N=10000, seed=70):
    """
    Does carry bias from block 1 propagate into block 2?

    Block 1: W1[0..15] → H1 = compress(IV, W1)
    Block 2: W2[0..15] → H2 = compress(H1, W2)

    If H1 is biased (e.g. from carry[63]=0), does H2 inherit any bias?
    """
    print("=" * 70)
    print("PATH A: Multiblock Carry Propagation")
    print(f"N={N}, seed={seed}")
    print("=" * 70)

    rng = np.random.RandomState(seed)

    # Phase 1: Find W1 with carry[63]=0 via Hill Climbing
    print("\n--- Phase 1: Finding carry[63]=0 inputs ---")

    biased_IVs = []  # H1 from block1 where carry[63]=0
    normal_IVs = []  # H1 from block1 where carry[63]=1
    t0 = time.time()
    target = min(N // 20, 200)

    for attempt in range(N * 10):
        W0 = int(rng.randint(0, 1 << 32))

        # Hill Climb to minimize raw[63]
        best_W0 = W0
        _, _, best_c = sha256_compress(tuple(H0), [W0] + [0]*15)
        # Compute SS63 for scoring
        _, best_states, _ = sha256_compress(tuple(H0), [best_W0] + [0]*15)
        st62 = best_states[63]
        best_ss = st62[7] + Sig1(st62[4]) + Ch(st62[4], st62[5], st62[6])

        for step in range(150):
            b = rng.randint(0, 32)
            W0_try = best_W0 ^ (1 << b)
            _, states_try, _ = sha256_compress(tuple(H0), [W0_try] + [0]*15)
            st62_try = states_try[63]
            ss_try = st62_try[7] + Sig1(st62_try[4]) + Ch(st62_try[4], st62_try[5], st62_try[6])
            if ss_try < best_ss:
                best_ss = ss_try
                best_W0 = W0_try

        H1, c1 = sha256_hash([best_W0] + [0]*15)
        if c1[63] == 0:
            biased_IVs.append(H1)
        elif len(normal_IVs) < target:
            normal_IVs.append(H1)

        if len(biased_IVs) >= target and len(normal_IVs) >= target:
            break
        if (attempt + 1) % 2000 == 0:
            print(f"  {attempt+1}: biased={len(biased_IVs)}/{target}, normal={len(normal_IVs)}/{target}")

    n_biased = len(biased_IVs)
    n_normal = len(normal_IVs)

    print(f"  Found: {n_biased} biased IVs, {n_normal} normal IVs ({time.time()-t0:.1f}s)")

    if n_biased < 10:
        print("  Too few biased IVs, skipping multiblock test.")
        return

    # Phase 2: Run block 2 with random W2, compare carry profiles
    print(f"\n--- Phase 2: Block 2 with biased vs normal IV ---")

    def measure_block2(IVs, label, n_per_iv=5):
        carries_63 = []
        h7_bits = {b: [] for b in [28, 29, 30, 31]}
        carry_sums = []

        for iv in IVs[:min(500, len(IVs))]:
            for _ in range(n_per_iv):
                W2 = [int(rng.randint(0, 1 << 32)) for _ in range(16)]
                H2, _, c2 = sha256_compress(iv, W2)
                carries_63.append(c2[63])
                carry_sums.append(sum(c2))
                for b in [28, 29, 30, 31]:
                    h7_bits[b].append((H2[7] >> b) & 1)

        print(f"\n  {label} (N={len(carries_63)}):")
        print(f"    P(carry2[63]=0) = {1 - np.mean(carries_63):.5f}")
        print(f"    E[carry_sum2]   = {np.mean(carry_sums):.2f}")
        for b in [29, 30, 31]:
            p = np.mean(h7_bits[b])
            print(f"    P(H2[7][b{b}]=1) = {p:.4f} (delta={p-0.5:+.4f})")

        return carries_63, h7_bits

    c_biased, hb_biased = measure_block2(biased_IVs, "Biased IV (carry1[63]=0)")
    c_normal, hb_normal = measure_block2(normal_IVs, "Normal IV (carry1[63]=1)")

    # Compare
    print(f"\n  --- Comparison ---")
    p_b = 1 - np.mean(c_biased)
    p_n = 1 - np.mean(c_normal)
    print(f"  P(carry2[63]=0 | biased IV):  {p_b:.5f}")
    print(f"  P(carry2[63]=0 | normal IV):  {p_n:.5f}")
    print(f"  Ratio: {p_b/max(p_n, 1e-10):.3f}x")

    for b in [29, 30, 31]:
        pb = np.mean(hb_biased[b])
        pn = np.mean(hb_normal[b])
        print(f"  H2[7][b{b}]: biased={pb:.4f}, normal={pn:.4f}, delta={pb-pn:+.4f}")


# ============================================================
# PATH B: BACKWARD — from H to state[62]
# ============================================================

def path_B(N=10000, seed=71):
    """
    Backward analysis: what does H[7] tell us about internal state?

    H[i] = state[63][i] + IV[i]  →  state[63][i] = H[i] - IV[i]
    state[63] is fully determined by H. No information loss.

    But: can we go further back? state[62] from state[63]?
    Inverse round: given (a,b,c,d,e,f,g,h) after round 63,
    recover (a,b,c,d,e,f,g,h) before round 63.

    This requires W[63] — which is determined by W[0..15].
    """
    print("\n" + "=" * 70)
    print("PATH B: Backward Analysis — H → state[62] → state[60]")
    print(f"N={N}, seed={seed}")
    print("=" * 70)

    rng = np.random.RandomState(seed)

    # Backward round
    def inverse_round(state_after, r, W):
        a, b, c, d, e, f, g, h = state_after
        # b_after = a_before → a_before = b
        # c_after = b_before → b_before = c
        # d_after = c_before → c_before = d
        # f_after = e_before ... but e_after = d_before + T1
        # So: a_before = b, b_before = c, c_before = d
        #     e_before = f, f_before = g, g_before = h

        # T2 = Sig0(a_before) + Maj(a_before, b_before, c_before)
        #    = Sig0(b) + Maj(b, c, d)
        T2 = (Sig0(b) + Maj(b, c, d)) & MASK

        # T1 = a_after - T2
        T1 = (a - T2) & MASK

        # d_before = e_after - T1
        d_before = (e - T1) & MASK

        # h_before = T1 - Sig1(e_before) - Ch(e_before, f_before, g_before) - K[r] - W[r]
        #          = T1 - Sig1(f) - Ch(f, g, h) - K[r] - W[r]
        h_before = (T1 - Sig1(f) - Ch(f, g, h) - K[r] - W[r]) & MASK

        return (b, c, d, d_before, f, g, h, h_before)

    # Test inverse
    print("\n--- Verification: inverse round correctness ---")
    n_verify = 100
    errors = 0
    for _ in range(n_verify):
        W16 = [int(rng.randint(0, 1 << 32)) for _ in range(16)]
        H, states, carries = sha256_compress(tuple(H0), W16)
        W64 = message_schedule(W16)

        # Recover state[63] from H
        state63_recovered = tuple((H[i] - H0[i]) & MASK for i in range(8))
        if state63_recovered != states[64]:
            errors += 1

        # Inverse round 63 → state[62]
        state62_recovered = inverse_round(states[64], 63, W64)
        if state62_recovered != states[63]:
            errors += 1

        # Chain: inverse rounds 63 → 60
        state = states[64]
        for r in range(63, 59, -1):
            state = inverse_round(state, r, W64)
        if state != states[60]:
            errors += 1

    print(f"  Verified {n_verify} traces: {errors} errors")
    print(f"  Inverse round: {'CORRECT' if errors == 0 else 'ERRORS FOUND'}")

    # Key analysis: what does state[62] look like when carry[63]=0?
    print(f"\n--- State[62] structure when carry[63]=0 ---")

    state62_carry0 = []
    state62_normal = []

    for i in range(N):
        W16 = [int(rng.randint(0, 1 << 32))] + [0] * 15
        H, states, carries = sha256_compress(tuple(H0), W16)

        if carries[63] == 0:
            state62_carry0.append(states[63])  # state after round 62
        elif len(state62_normal) < 500:
            state62_normal.append(states[63])

    print(f"  carry[63]=0 samples: {len(state62_carry0)}")
    print(f"  carry[63]=1 samples: {len(state62_normal)}")

    if len(state62_carry0) >= 5:
        reg_names = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
        print(f"\n  {'reg':>5} | {'E[HW] c=0':>10} | {'E[HW] c=1':>10} | {'delta':>7} | {'E[val]/2^32 c=0':>16}")
        print("  " + "-" * 65)
        for i in range(8):
            vals_0 = [s[i] for s in state62_carry0]
            vals_1 = [s[i] for s in state62_normal]
            hw_0 = np.mean([hw(v) for v in vals_0])
            hw_1 = np.mean([hw(v) for v in vals_1])
            mean_0 = np.mean(vals_0) / (1 << 32)
            marker = " ***" if abs(hw_0 - hw_1) > 1.0 else ""
            print(f"  {reg_names[i]:>5} | {hw_0:10.2f} | {hw_1:10.2f} | {hw_0-hw_1:+7.2f} | {mean_0:16.4f}{marker}")

        # Which bits of state[62] are most constrained?
        print(f"\n  Per-bit bias of state[62] registers when carry[63]=0:")
        for reg_i, reg_name in enumerate(reg_names):
            biases = []
            for b in range(32):
                p = np.mean([(s[reg_i] >> b) & 1 for s in state62_carry0])
                biases.append(abs(p - 0.5))
            max_bias = max(biases)
            max_bit = biases.index(max_bias)
            if max_bias > 0.05:
                p_val = np.mean([(s[reg_i] >> max_bit) & 1 for s in state62_carry0])
                print(f"    {reg_name}[b{max_bit}]: P(=1|c=0)={p_val:.3f} (bias={max_bias:.3f}) ***")

    return state62_carry0, state62_normal


# ============================================================
# PATH C: K-INVARIANTS — forced bits from constants
# ============================================================

def path_C(N=10000, seed=72):
    """
    K-invariants: which bits are FORCED by the magnitude of K[r]?

    For carry[r]=0: raw[r] = h + Sig1(e) + Ch(e,f,g) + K[r] + W[r] < 2^32
    SS[r] = h + Sig1(e) + Ch(e,f,g)

    carry[r]=0  ⟺  SS[r] + K[r] + W[r] < 2^32
                ⟺  SS[r] < 2^32 - K[r] - W[r]

    When K[r] is large (close to 2^32), SS[r] must be VERY small.
    This forces high bits of SS components (h, Sig1, Ch) to 0.

    Build the full K-constraint map for all 64 rounds.
    """
    print("\n" + "=" * 70)
    print("PATH C: K-Invariants — Forced Bits from Constants")
    print(f"N={N}, seed={seed}")
    print("=" * 70)

    # Phase 1: Analytical — max SS allowed for carry=0 at each round
    print(f"\n--- Phase 1: K-constraint budget per round ---")
    print(f"  carry[r]=0 ⟺ SS[r] < 2^32 - K[r] - W[r]")
    print(f"  Worst case (W[r]=0): SS[r] < 2^32 - K[r]")
    print(f"  Best case (W[r]=2^32-1): SS[r] < 2^32 - K[r] - (2^32-1) = 1-K[r] (impossible if K[r]>0)")
    print()

    print(f"  {'r':>3} | {'K[r]':>12} | {'K[r]/2^32':>10} | {'budget':>12} | {'budget/2^32':>12} | {'forced_bits':>12}")
    print("  " + "-" * 75)

    k_budgets = []
    forced_bits_per_round = []
    for r in range(64):
        budget = ((1 << 32) - K[r]) & MASK  # max SS for carry=0 (W[r]=0 best case)
        k_ratio = K[r] / (1 << 32)
        b_ratio = budget / (1 << 32)

        # How many high bits of SS are forced to 0?
        forced = 0
        for b in range(31, -1, -1):
            if budget < (1 << (b + 1)):
                forced += 1
            else:
                break

        k_budgets.append(budget)
        forced_bits_per_round.append(forced)

        marker = " ***" if forced >= 2 else (" *" if forced >= 1 else "")
        print(f"  {r:3d} | {K[r]:12d} | {k_ratio:10.4f} | {budget:12d} | {b_ratio:10.4f} | {forced:12d}{marker}")

    # Summary
    total_forced = sum(forced_bits_per_round)
    rounds_with_forced = sum(1 for f in forced_bits_per_round if f > 0)
    print(f"\n  Rounds with forced bits (budget < 2^31): {rounds_with_forced}/64")
    print(f"  Total forced high bits: {total_forced}")
    print(f"  Max forced bits: {max(forced_bits_per_round)} at r={forced_bits_per_round.index(max(forced_bits_per_round))}")

    # Top-10 most constrained rounds
    top_rounds = sorted(range(64), key=lambda r: forced_bits_per_round[r], reverse=True)[:10]
    print(f"\n  Top-10 most constrained rounds:")
    for r in top_rounds:
        print(f"    r={r:2d}: {forced_bits_per_round[r]} forced bits, K=0x{K[r]:08x}, budget=0x{k_budgets[r]:08x}")

    # Phase 2: Experimental — verify forced bits
    print(f"\n--- Phase 2: Verify forced bits experimentally ---")
    rng = np.random.RandomState(seed)

    # For each highly constrained round, measure actual SS distribution
    for r in top_rounds[:5]:
        if forced_bits_per_round[r] < 1:
            continue

        ss_vals = []
        carry0_count = 0
        ss_carry0 = []

        for _ in range(N):
            W16 = [int(rng.randint(0, 1 << 32)) for _ in range(16)]
            _, states, carries = sha256_compress(tuple(H0), W16)

            # SS[r] = h + Sig1(e) + Ch(e,f,g) at state before round r
            st = states[r]
            a, b, c, d, e, f, g, h = st
            ss = h + Sig1(e) + Ch(e, f, g)  # without mod
            ss_vals.append(ss)

            if carries[r] == 0:
                carry0_count += 1
                ss_carry0.append(ss)

        print(f"\n  Round {r}: K=0x{K[r]:08x}, budget=0x{k_budgets[r]:08x}")
        print(f"    P(carry=0): {carry0_count/N:.5f}")
        print(f"    E[SS]: {np.mean(ss_vals)/2**32:.4f} × 2^32")
        if ss_carry0:
            print(f"    E[SS|carry=0]: {np.mean(ss_carry0)/2**32:.4f} × 2^32")
            print(f"    max(SS|carry=0): 0x{max(ss_carry0):x}")
            print(f"    max(SS|carry=0) < budget? {max(ss_carry0) < k_budgets[r] + (1<<32)}")

    # Phase 3: Chain of forced conditions
    print(f"\n--- Phase 3: K-invariant chain (carry=0 at round r forces what?) ---")

    # For the most constrained round (r=63): trace which registers are forced
    r_target = 63
    print(f"\n  Target round r={r_target}: K=0x{K[r_target]:08x}")
    print(f"  carry[{r_target}]=0 ⟹ SS[{r_target}] < {k_budgets[r_target]} (0x{k_budgets[r_target]:08x})")
    print(f"  SS = h[{r_target-1}] + Sig1(e[{r_target-1}]) + Ch(e[{r_target-1}],f[{r_target-1}],g[{r_target-1}])")

    # What about the CHAIN: carry[63]=0 → constraints on state[62]
    # → carry[62] more likely 0 → constraints on state[61] → ...
    print(f"\n  Chain propagation (conditional probabilities):")

    chain_data = defaultdict(lambda: {'total': 0, 'carry0': 0})

    for _ in range(N):
        W16 = [int(rng.randint(0, 1 << 32))] + [0] * 15
        _, states, carries = sha256_compress(tuple(H0), W16)

        # Build chain from r=63 backward
        if carries[63] == 0:
            for r2 in range(62, 55, -1):
                key = f"c{r2}|c63=0"
                chain_data[key]['total'] += 1
                if carries[r2] == 0:
                    chain_data[key]['carry0'] += 1

            # Double condition
            if carries[62] == 0:
                for r2 in range(61, 55, -1):
                    key = f"c{r2}|c63=c62=0"
                    chain_data[key]['total'] += 1
                    if carries[r2] == 0:
                        chain_data[key]['carry0'] += 1

    for key in sorted(chain_data.keys()):
        d = chain_data[key]
        if d['total'] > 0:
            p = d['carry0'] / d['total']
            print(f"    P({key.split('|')[0]}=0 | {key.split('|')[1]}) = {p:.4f} (N={d['total']})")

    return k_budgets, forced_bits_per_round


# ============================================================
# MAIN
# ============================================================

if __name__ == '__main__':
    print("Three Paths Around Chaos — Stage 1.4")
    print(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print()

    N_A = 8000
    N_B = 8000
    N_C = 8000

    if len(sys.argv) > 1 and sys.argv[1] == '--fast':
        N_A, N_B, N_C = 3000, 3000, 3000

    t_start = time.time()

    path_A(N=N_A)
    state62_c0, state62_c1 = path_B(N=N_B)
    k_budgets, forced_bits = path_C(N=N_C)

    total = time.time() - t_start

    print(f"\n{'='*70}")
    print(f"TOTAL TIME: {total:.1f}s")
    print(f"{'='*70}")

    print(f"""
{'='*70}
SYNTHESIS: Three Paths Around Chaos
{'='*70}

PATH A (Multiblock): Does biased IV from block 1 propagate?
  → If yes: chain attack accumulates bias over blocks
  → If no: each block is independent barrier

PATH B (Backward): From H to state[62] — what's constrained?
  → H[7] = g[62]+IV[7] → g[62] known exactly from H[7]
  → carry[63]=0 forces specific register bits
  → Inverse round is exact (verified)

PATH C (K-invariants): Constants create deterministic constraints
  → {sum(1 for f in forced_bits if f > 0)} rounds have forced high bits when carry=0
  → Total forced bits across all rounds: {sum(forced_bits)}
  → r=63: K=0xc67178f2, budget=0x{k_budgets[63]:08x}
     carry[63]=0 forces SS<{k_budgets[63]} — bits 31,30 of Ch MUST be 0
     (This is T_CH_INVARIANT from methodology!)

THE QUESTION: Which path (or combination) gives new mathematics?
""")
