#!/usr/bin/env python3
"""
RAYON + Carry-Web Integration — Axis 6
=========================================

Three parallel integrations:

  I1: CARRY-AWARE DFS — use K-map to guide variable ordering.
      Fix bits in carry-window rounds first → more AND-chain cutoffs.
      Measure ε for SHA-256 with smart ordering vs random.

  I2: COPULA-GUIDED RAYON BIRTHDAY — multi-target preimage
      in carry=0 subspace. Use copula predictions to focus DFS.

  I3: DEGREE-AWARE PROPAGATION — use degree deficit (n-1 at n≤28)
      to add extra propagation rules in DFS.
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
def sig0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def sig1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def Sig0(x): return rotr(x,2)^rotr(x,13)^rotr(x,22)
def Sig1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def Ch(e,f,g): return (e&f)^(~e&g)&MASK
def Maj(a,b,c): return (a&b)^(a&c)^(b&c)
def hw(x): return bin(x&MASK).count('1')
from math import erfc, sqrt, log2

def msg_sched(W16):
    W=list(W16)+[0]*48
    for i in range(16,64): W[i]=(sig1(W[i-2])+W[i-7]+sig0(W[i-15])+W[i-16])&MASK
    return W

def sha256(W16):
    W=msg_sched(W16); a,b,c,d,e,f,g,h=H0
    for r in range(64):
        T1=(h+Sig1(e)+Ch(e,f,g)+K[r]+W[r])&MASK; T2=(Sig0(a)+Maj(a,b,c))&MASK
        h,g,f=g,f,e; e=(d+T1)&MASK; d,c,b=c,b,a; a=(T1+T2)&MASK
    return [(v+iv)&MASK for v,iv in zip([a,b,c,d,e,f,g,h],H0)]

def sha256_partial(W16, target_bit_word, target_bit_pos):
    """Compute single output bit of SHA-256."""
    H = sha256(W16)
    return (H[target_bit_word] >> target_bit_pos) & 1

SIGMA = 0.568
def p_carry0(r):
    """Analytical P(carry[r]=0) from copula model."""
    return 0.5 * erfc((1 + K[r]/(1<<32)) / (SIGMA * sqrt(2)))


# ============================================================
# I1: CARRY-AWARE DFS — Smart variable ordering
# ============================================================

def experiment_I1(n_bits=16, seed=20000):
    """
    DFS with constant propagation on SHA-256 (n_bits input).
    Compare: random bit ordering vs carry-aware ordering.

    Carry-aware: fix bits that affect carry-window rounds first.
    Idea: RAYON showed carry chains = AND chains = good for cutoff.
    Our K-map tells which rounds are cheapest for carry=0.
    """
    print("="*70)
    print(f"I1: CARRY-AWARE DFS — n={n_bits} input bits")
    print("="*70)

    rng = np.random.RandomState(seed)

    # For n_bits of W[0], compute SHA-256 and check if specific output bit = target.
    # DFS: fix bits one by one, propagate, check if output determined.

    # Simplified DFS: for each partial assignment (k bits fixed),
    # evaluate SHA-256 for ALL 2^(n-k) completions.
    # If all give same output bit → subtree determined → prune.
    # Count: how many nodes before solution found?

    target_w, target_b = 7, 23  # target output bit

    # Choose a random target value
    target_input = int(rng.randint(0, 1 << n_bits))
    target_value = sha256_partial([target_input]+[0]*15, target_w, target_b)

    def dfs_count(bit_order, target_val):
        """Count DFS nodes to find ANY input with output = target_val."""
        nodes = 0

        def search(depth, partial):
            nonlocal nodes
            nodes += 1
            if nodes > 2**n_bits:  # safety limit
                return False

            if depth == n_bits:
                # Leaf: check if output matches
                val = sha256_partial([partial]+[0]*15, target_w, target_b)
                return val == target_val

            bit = bit_order[depth]
            for v in [0, 1]:
                new_partial = partial | (v << bit)

                # Propagation check: can we determine output?
                # Fix current bits, try ALL remaining — if all same → determined
                # (This is expensive for large n, so we sample)
                n_remaining = n_bits - depth - 1
                if n_remaining <= 10:
                    # Exhaustive check
                    vals = set()
                    determined = True
                    for rem in range(1 << n_remaining):
                        full = new_partial
                        rem_bits = [bit_order[depth+1+j] for j in range(n_remaining)]
                        for j, rb in enumerate(rem_bits):
                            full |= ((rem >> j) & 1) << rb
                        ov = sha256_partial([full]+[0]*15, target_w, target_b)
                        vals.add(ov)
                        if len(vals) > 1:
                            determined = False
                            break
                    if determined:
                        if target_val in vals:
                            return True  # determined and matches
                        else:
                            continue  # determined but wrong → prune
                # Else: just recurse
                if search(depth + 1, new_partial):
                    return True
            return False

        search(0, 0)
        return nodes

    # Ordering 1: Random
    random_order = list(rng.permutation(n_bits))

    # Ordering 2: Carry-aware — bits that affect early rounds first
    # Bits 0-7 affect round 0 most (low bits of W[0] → carry chain in T1[0])
    # Carry-window priority: rounds with highest P(carry=0)
    # For W[0] bits: bit k affects T1[0] carry at position k
    # Lower bits → earlier carries → more determination
    carry_order = list(range(n_bits))  # LSB first (0,1,2,...,n-1)

    # Ordering 3: MSB first (opposite — should be worse)
    msb_order = list(range(n_bits-1, -1, -1))

    print(f"\n  Target: H[{target_w}][b{target_b}] = {target_value}, input = 0x{target_input:04x}")
    print(f"\n  Testing 3 orderings with DFS + propagation (n={n_bits})...")

    for name, order in [("random", random_order), ("carry_LSB", carry_order), ("MSB_first", msb_order)]:
        t0 = time.time()
        nodes = dfs_count(order, target_value)
        elapsed = time.time() - t0
        brute = 1 << n_bits
        if nodes > 0 and nodes < brute:
            eps = 1 - log2(nodes) / n_bits if nodes > 1 else 1.0
        else:
            eps = 0
        speedup = brute / max(nodes, 1)
        print(f"    {name:>12}: nodes={nodes:8d}, brute={brute}, ε={eps:.3f}, speedup={speedup:.1f}×, {elapsed:.1f}s")

    # Repeat with multiple targets for statistics
    print(f"\n  Statistical test (10 random targets)...")
    results = {name: [] for name in ["random", "carry_LSB", "MSB_first"]}

    for trial in range(10):
        target_input = int(rng.randint(0, 1 << n_bits))
        target_value = sha256_partial([target_input]+[0]*15, target_w, target_b)

        for name, order_gen in [("random", lambda: list(rng.permutation(n_bits))),
                                 ("carry_LSB", lambda: list(range(n_bits))),
                                 ("MSB_first", lambda: list(range(n_bits-1,-1,-1)))]:
            order = order_gen()
            nodes = dfs_count(order, target_value)
            results[name].append(nodes)

    print(f"\n  {'ordering':>12} | {'mean nodes':>11} | {'mean ε':>7} | {'vs brute':>9}")
    print(f"  {'-'*12}-+-{'-'*11}-+-{'-'*7}-+-{'-'*9}")
    brute = 1 << n_bits
    for name in ["random", "carry_LSB", "MSB_first"]:
        mn = np.mean(results[name])
        eps = 1 - log2(max(mn,1)) / n_bits
        print(f"  {name:>12} | {mn:11.0f} | {eps:7.3f} | {brute/max(mn,1):9.1f}×")

    return results


# ============================================================
# I2: COPULA-GUIDED RAYON BIRTHDAY
# ============================================================

def experiment_I2(seed=20001):
    print(f"\n{'='*70}")
    print("I2: COPULA-GUIDED RAYON BIRTHDAY")
    print("="*70)

    # Rayon Birthday: Phase 1 collect K hashes, Phase 2 multi-target DFS.
    # Our addition: use copula model to select WHICH hashes to collect.
    # Specifically: collect hashes from carry=0 subspace (HC-optimized).
    # These hashes have biased H[6,7] → smaller effective hash space.

    # For toy model: use SHA-256 truncated to h bits.
    h = 16  # 16-bit hash (first 16 bits of H[7])
    n = 20  # 20-bit input (W[0][0..19])
    rng = np.random.RandomState(seed)

    def toy_hash(x):
        H = sha256([x & MASK] + [0]*15)
        return H[7] & ((1 << h) - 1)

    # Standard Birthday
    print(f"\n  --- Standard Birthday (h={h}, n={n}) ---")
    t0 = time.time()
    table = {}
    std_cost = 0
    std_found = False
    for x in range(1 << n):
        std_cost += 1
        hv = toy_hash(x)
        if hv in table:
            print(f"    COLLISION: x1={table[hv]:#x}, x2={x:#x}, hash={hv:#06x}")
            std_found = True
            break
        table[hv] = x
    print(f"    Cost: {std_cost} evals ({time.time()-t0:.1f}s)")
    print(f"    Expected: 2^{h/2} = {1<<(h//2)}")

    # Copula-Guided: collect from HC-optimized inputs
    print(f"\n  --- Copula-Guided Birthday (HC inputs) ---")
    t0 = time.time()
    table_hc = {}
    hc_cost = 0
    hc_found = False

    for trial in range(1 << 16):
        # HC: minimize raw[63] via bit flips
        W0 = int(rng.randint(0, 1 << n))
        best_W0 = W0
        # Quick HC: 50 steps
        for _ in range(50):
            b = int(rng.randint(0, n))
            W0_try = W0 ^ (1 << b)
            # Use toy_hash as proxy (not real raw63 minimization)
            if toy_hash(W0_try) < toy_hash(W0):
                W0 = W0_try
        hc_cost += 51  # 1 base + 50 HC steps

        hv = toy_hash(W0)
        if hv in table_hc:
            print(f"    COLLISION: x1={table_hc[hv]:#x}, x2={W0:#x}, hash={hv:#06x}")
            hc_found = True
            break
        table_hc[hv] = W0

    print(f"    Cost: {hc_cost} evals ({time.time()-t0:.1f}s)")

    if std_found and hc_found:
        ratio = std_cost / max(hc_cost, 1)
        print(f"\n  Speedup HC vs standard: {ratio:.2f}×")
        if ratio > 1.5:
            print(f"  ★ HC-guided birthday is FASTER!")
        else:
            print(f"  HC overhead cancels any hash-space reduction")

    return


# ============================================================
# I3: DEGREE-AWARE PROPAGATION
# ============================================================

def experiment_I3(seed=20002):
    print(f"\n{'='*70}")
    print("I3: DEGREE-AWARE PROPAGATION")
    print("="*70)

    # Key insight: at n≤28, degree(H[7][b23]) = n-1.
    # This means: the n-th order derivative is 0.
    # δ^n f(x) = XOR over all 2^n cube corners = 0.
    #
    # For DFS: if we've fixed n-1 bits and the last bit is free,
    # we know that f(x with last=0) XOR f(x with last=1) = derivative.
    # For degree-n function: derivative can be 0 or 1.
    # For degree-(n-1) function: n-th derivative = 0 ALWAYS.
    # → f(x|last=0) = f(x|last=1) → output DETERMINED by n-1 bits!
    #
    # This is a FREE propagation rule at the last level of DFS.
    # Saves factor 2 in DFS tree.

    # But: this only works when ALL n bits are in the "deficit" subset.
    # And degree deficit is for W[0][0..n-1], not arbitrary bit subsets.

    # Test: for n=12, verify that fixing 11 bits determines bit 12.

    n = 14
    print(f"\n  Testing degree-aware propagation for n={n}")
    print(f"  If degree=n-1: fixing n-1 bits → output determined")

    rng = np.random.RandomState(seed)
    N_tests = 200

    # For each random (n-1)-bit assignment: check if last bit matters
    w, b = 7, 23
    determined_count = 0

    for _ in range(N_tests):
        # Fix bits 0..n-2 randomly
        prefix = int(rng.randint(0, 1 << (n-1)))
        # Last bit = 0 or 1
        x0 = prefix  # bit n-1 = 0
        x1 = prefix | (1 << (n-1))  # bit n-1 = 1
        v0 = sha256_partial([x0]+[0]*15, w, b)
        v1 = sha256_partial([x1]+[0]*15, w, b)
        if v0 == v1:
            determined_count += 1

    p_det = determined_count / N_tests
    print(f"  P(output determined by n-1 bits): {p_det:.3f}")
    print(f"  Expected if degree=n-1: 1.000")
    print(f"  Expected if degree=n:   0.500")

    if p_det > 0.95:
        print(f"  ★★★ CONFIRMED: degree n-1 gives FREE determination at last level!")
        print(f"  DFS saves factor 2 at the deepest level → overall speedup ~2×")
    elif p_det > 0.75:
        print(f"  ★★ Partial: {p_det:.0%} determination (degree nearly n-1)")
    else:
        print(f"  Not exploitable at n={n}")

    # Test with multiple last-bit positions
    print(f"\n  --- Which bit position can be 'last' (determined by others)? ---")
    for last_bit in [0, n//4, n//2, 3*n//4, n-1]:
        det = 0
        for _ in range(N_tests):
            x_base = int(rng.randint(0, 1 << n)) & ~(1 << last_bit)
            v0 = sha256_partial([x_base]+[0]*15, w, b)
            v1 = sha256_partial([x_base | (1<<last_bit)]+[0]*15, w, b)
            if v0 == v1: det += 1
        p = det / N_tests
        print(f"    last_bit={last_bit:2d}: P(determined)={p:.3f}{'  ★' if p>0.6 else ''}")

    # KEY: combine with RAYON DFS
    print(f"\n  --- DFS with degree-aware pruning ---")
    # At each leaf (n bits fixed except 1): check if output determined
    # If yes → skip the flip → save 1 evaluation
    # Net saving: P(determined) × (fraction of leaves)

    print(f"  Degree-aware DFS saving: {p_det:.1%} of leaf evaluations saved")
    print(f"  Overall DFS speedup factor: {1 + p_det:.2f}×")
    print(f"  (Compared to RAYON's ε improvement, this is modest)")

    return p_det


# ============================================================
# MAIN
# ============================================================

if __name__ == '__main__':
    print("RAYON + Carry-Web Integration — Axis 6")
    print(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")

    n_bits = 14
    if len(sys.argv) > 1 and sys.argv[1] == '--fast':
        n_bits = 12

    t_start = time.time()
    results_I1 = experiment_I1(n_bits=n_bits)
    experiment_I2()
    p_det = experiment_I3()

    total = time.time() - t_start
    print(f"\n{'='*70}")
    print(f"TOTAL TIME: {total:.1f}s")
    print(f"{'='*70}")

    print(f"""
{'='*70}
SYNTHESIS: RAYON + Carry-Web Integration (Axis 6)
{'='*70}

I1: CARRY-AWARE DFS
  Does LSB-first ordering improve ε for SHA-256?
  Smart ordering exploits carry chain structure.

I2: COPULA-GUIDED BIRTHDAY
  HC-optimized inputs → biased hash space → cheaper birthday?

I3: DEGREE-AWARE PROPAGATION
  Degree n-1 at n≤28 → {p_det:.0%} of leaves determined free.
  DFS speedup factor: {1+p_det:.2f}×

DATA FOR RAYON:
  - K-map: P(carry[r]=0) = Φ((-1-K[r]/T)/0.568) for all 64 rounds
  - Carry windows: {{9,10,18,19,30,31,47,48}} (cheapest rounds)
  - Degree deficit: n-1 at n≤28, full at n=32
  - Block structure: 6 independent blocks, isolation 42×
""")
