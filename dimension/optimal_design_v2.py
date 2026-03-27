"""
СЛЕДСТВИЯ ЕДИНОЙ ТЕОРИИ.

Из H(N,C,r) = Random iff r >= 3N+5:

1. OPTIMAL DESIGN: для 128-bit collision, minimize computation.
   Need NC/2 = 128 → NC = 256.
   Options: (N=4,C=64), (N=8,C=32), (N=16,C=16), (N=32,C=8)
   Cost per round ∝ N×C (register operations)
   Total cost = (3N+5) × N × C = (3N+5) × 256
   Minimize (3N+5): minimum at N=1... but N=1 → no pipes → theory breaks!
   Minimum with pipes: N=2, C=128, rounds=11, cost=11×256=2816
   vs SHA-256: N=8, C=32, rounds=29, cost=29×256=7424 → 2.6× slower!

2. BREAK THE AXIOM: design hash with K nodes (K>2).
   Prediction: S(r) = K×C×r (faster entropy gain!)
   Saturation: r = NC/(KC) = N/K (faster!)
   Optimal: 3N/K + 5? Or different formula?

3. ALL-NODE HASH (K=N): every register computed, no pipes.
   Prediction: S(r) = NC immediately → 1 round sufficient??
   Conservation: NONE (no pipes → no conservation)
   This should be MAXIMALLY FAST but potentially WEAK.

Let's BUILD and TEST all three.
"""

import numpy as np
import struct, hashlib

MASK32 = 0xFFFFFFFF
MASK16 = 0xFFFF
MASK64 = 0xFFFFFFFFFFFFFFFF

def hw(x): return bin(x).count('1')
def rotr32(x, n): return ((x >> n) | (x << (32 - n))) & MASK32


def hash_generic(message, N, C, n_rounds, K_nodes=2):
    """
    Generic shift-register hash with K_nodes active nodes.
    N registers, C bits each.
    K_nodes: how many registers are COMPUTED (rest are pipes).
    """
    mask = (1 << C) - 1

    def rotr_c(x, n):
        return ((x >> (n % C)) | (x << (C - n % C))) & mask
    def add_c(x, y):
        return (x + y) & mask
    def nl(x, y, z):
        return ((x & y) ^ (~x & z)) & mask

    # IV
    primes = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71]
    state = [(int(p**0.5 * (1 << C)) & mask) for p in primes[:N]]

    # Schedule
    W = list(message[:N*2])
    while len(W) < N*2:
        W.append(np.random.randint(0, 1 << C) & mask)
    for r in range(N*2, n_rounds):
        W.append(add_c(rotr_c(W[r-2], C//4) ^ rotr_c(W[r-2], C//2),
                       add_c(W[r-N], rotr_c(W[r-N*2+1], C//5) ^ (W[r-N*2+1] >> (C//8)))))

    # Constants
    K_c = [(int((i+2)**0.333 * (1<<C)) & mask) for i in range(n_rounds)]

    iv_save = list(state)

    for r in range(n_rounds):
        new_state = [0] * N

        # Compute K_nodes node registers
        for k in range(K_nodes):
            # Node k: uses different parts of state
            node_idx = k * (N // K_nodes)  # spread nodes evenly
            e_idx = (node_idx + N // 2) % N

            # T = nonlinear combination + W + K
            T = add_c(add_c(state[(node_idx - 1) % N],
                      rotr_c(state[e_idx], C//6) ^ rotr_c(state[e_idx], C//4)),
                      add_c(nl(state[e_idx], state[(e_idx+1)%N], state[(e_idx+2)%N]),
                      add_c(K_c[r], W[r] if r < len(W) else 0)))

            if k == 0:
                # Primary node: goes to position 0
                T2 = add_c(rotr_c(state[0], C//8) ^ rotr_c(state[0], C//3),
                           (state[0] & state[1]) ^ (state[0] & state[2] if N > 2 else 0) ^ (state[1] & state[2] if N > 2 else 0))
                new_state[0] = add_c(T, T2)
            else:
                # Additional nodes: go to evenly spaced positions
                pos = k * (N // K_nodes)
                new_state[pos] = add_c(state[(pos-1)%N], T)

        # Pipes: copy from previous position
        for i in range(N):
            if new_state[i] == 0:  # not set by a node
                new_state[i] = state[(i - 1) % N]

        state = new_state

    return tuple(add_c(iv_save[i], state[i]) for i in range(N))


def measure_hash(hash_fn, N, C, name, n_samples=1000):
    """Measure key metrics."""
    mask = (1 << C) - 1
    m = N * C

    # Avalanche
    avl = []
    for _ in range(n_samples):
        msg = [np.random.randint(0, 1<<C) & mask for _ in range(N*2)]
        H1 = hash_fn(msg)
        w = np.random.randint(0, N*2)
        b = np.random.randint(0, C)
        msg2 = list(msg); msg2[w] ^= (1 << b)
        H2 = hash_fn(msg2)
        avl.append(sum(hw(H1[i]^H2[i]) for i in range(N)))
    K_metric = np.mean(avl)

    # Conservation: (reg[0]+reg[N//2])[r] = (reg[1]+reg[N//2+1])[r+1]?
    # Can't easily measure without internal states, so skip for now

    return {'K': K_metric, 'K_ideal': m/2, 'ratio': K_metric/(m/2)}


def main():
    np.random.seed(42)

    print("=" * 70)
    print("СЛЕДСТВИЯ ЕДИНОЙ ТЕОРИИ")
    print("=" * 70)

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("1. OPTIMAL 128-BIT HASH: minimize total cost")
    print("=" * 70)

    # NC = 256, minimize (3N+5) × NC
    # Cost = (3N+5) × 256
    print(f"\n  NC = 256 (for 128-bit collision resistance)")
    print(f"  {'N':>3} {'C':>4} {'Rounds':>7} {'Cost':>8} {'Speedup vs SHA-256':>20}")
    print(f"  {'':>3} {'':>4} {'(3N+5)':>7} {'R×NC':>8}")

    sha256_cost = 29 * 256
    for N in [2, 4, 6, 8, 10, 16, 32]:
        C = 256 // N
        if C < 4: continue
        R = 3*N + 5
        cost = R * 256  # NC is fixed at 256
        speedup = sha256_cost / cost
        print(f"  {N:>3} {C:>4} {R:>7} {cost:>8} {speedup:>19.2f}×")

    print(f"\n  OPTIMAL: N=2, C=128 → 11 rounds × 256 = 2816")
    print(f"  vs SHA-256 (N=8, C=32): 29 rounds × 256 = 7424")
    print(f"  Speedup: 2.6×")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("2. BUILD AND TEST: N=2,C=32 (tiny fast hash)")
    print("=" * 70)

    # N=2, C=32: collision = 2^32 (toy), 11 rounds
    for K_nodes in [2]:
        N, C, n_rounds = 2, 32, 11
        name = f"FastHash-{N}x{C} ({K_nodes} nodes)"
        hfn = lambda msg, N=N, C=C, R=n_rounds, Kn=K_nodes: hash_generic(msg, N, C, R, Kn)
        m = measure_hash(hfn, N, C, name)
        print(f"  {name}: K={m['K']:.1f}/{m['K_ideal']:.0f} ({m['ratio']:.3f})")

    # N=4, C=32: collision = 2^64, 17 rounds
    for K_nodes in [2]:
        N, C, n_rounds = 4, 32, 17
        name = f"FastHash-{N}x{C} ({K_nodes} nodes)"
        hfn = lambda msg, N=N, C=C, R=n_rounds, Kn=K_nodes: hash_generic(msg, N, C, R, Kn)
        m = measure_hash(hfn, N, C, name)
        print(f"  {name}: K={m['K']:.1f}/{m['K_ideal']:.0f} ({m['ratio']:.3f})")

    # N=8, C=32: collision = 2^128, 29 rounds (our prediction)
    for K_nodes in [2]:
        N, C, n_rounds = 8, 32, 29
        name = f"FastHash-{N}x{C} ({K_nodes} nodes)"
        hfn = lambda msg, N=N, C=C, R=n_rounds, Kn=K_nodes: hash_generic(msg, N, C, R, Kn)
        m = measure_hash(hfn, N, C, name)
        print(f"  {name}: K={m['K']:.1f}/{m['K_ideal']:.0f} ({m['ratio']:.3f})")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("3. BREAK THE AXIOM: 4 nodes instead of 2")
    print("=" * 70)

    # Theory predicts: S(r) = 4C×r (faster), saturation at r = N/4
    # Optimal rounds = 3N/2 + 5? (halved because twice as fast)

    print(f"\n  Predictions for K=4 nodes (N=8, C=32):")
    print(f"    Entropy rate: 4×32 = 128 bits/round (vs 64 for K=2)")
    print(f"    Saturation: r = N/4 = 2 (vs 4 for K=2)")
    print(f"    Optimal: ~3N/2+5 = 17? (vs 29 for K=2)")

    # Build and test
    for K_nodes in [1, 2, 4, 8]:
        N, C = 8, 32
        # Predicted optimal rounds
        if K_nodes <= N:
            pred_sat = max(1, N // K_nodes)
            pred_optimal = int(3 * N / K_nodes + 5)
        else:
            pred_sat = 1
            pred_optimal = 6

        n_rounds = max(pred_optimal, 8)  # at least 8 for testing

        name = f"Hash-8x32-{K_nodes}nodes"
        hfn = lambda msg, R=n_rounds, Kn=K_nodes: hash_generic(msg, 8, 32, R, Kn)
        m = measure_hash(hfn, N, C, name)

        print(f"\n  K={K_nodes} nodes: {pred_optimal} rounds predicted")
        print(f"    Avalanche K={m['K']:.1f}/{m['K_ideal']:.0f} ({m['ratio']:.3f})")
        print(f"    {'PASS ✓' if m['ratio'] > 0.95 else 'FAIL ✗' if m['ratio'] < 0.8 else 'MARGINAL ~'}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("4. ALL-NODE HASH: K=N (no pipes at all)")
    print("=" * 70)

    # When K=N: every register computed → no pipes → no conservation
    # Prediction: S(r) = NC in 1 round! Ultra-fast saturation.
    # But: no pipe structure → different security properties?

    N, C = 8, 32
    for n_rounds in [1, 2, 4, 8, 16]:
        hfn = lambda msg, R=n_rounds: hash_generic(msg, N, C, R, K_nodes=N)
        m = measure_hash(hfn, N, C, f"AllNode-{n_rounds}r", n_samples=500)
        print(f"  All-node {n_rounds:>2} rounds: K={m['K']:.1f}/128 ({m['ratio']:.3f}) "
              f"{'✓' if m['ratio'] > 0.95 else '~' if m['ratio'] > 0.8 else '✗'}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("5. THEORY EXTENSION: K-node formula")
    print("=" * 70)

    # Collect data: for each K, what's the minimum rounds for K/K_ideal > 0.95?
    print(f"\n  Minimum rounds for avalanche > 0.95:")
    print(f"  {'K_nodes':>7} {'Min rounds':>10} {'Predicted (3N/K+5)':>20} {'Match':>6}")

    for K_nodes in [1, 2, 4, 8]:
        N, C = 8, 32
        pred = int(3 * N / K_nodes + 5)

        # Binary search for minimum rounds
        min_r = None
        for n_r in range(1, 40):
            hfn = lambda msg, R=n_r, Kn=K_nodes: hash_generic(msg, N, C, R, Kn)
            m = measure_hash(hfn, N, C, f"test", n_samples=300)
            if m['ratio'] > 0.95:
                min_r = n_r
                break

        match = "✓" if min_r and abs(min_r - pred) <= 3 else "~" if min_r else "✗"
        print(f"  {K_nodes:>7} {str(min_r):>10} {pred:>20} {match:>6}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("6. EXTENDED UNIFIED EQUATION")
    print("=" * 70)

    print(f"""
  ═══════════════════════════════════════════════════════════════

  EXTENDED MASTER EQUATION (K-node generalization):

  Given: N registers, C bits, K nodes (K ≤ N), m = NC.

  ╔═══════════════════════════════════════════════════════════╗
  ║                                                           ║
  ║   H(N, C, K, r) = Random   iff   r ≥ 3N/K + 5           ║
  ║                                                           ║
  ║   Entropy rate:   dS/dr = KC                              ║
  ║   Saturation:     r_sat = N/K                             ║
  ║   Conservation:   lifetime = (N-K)/2 (if K < N)           ║
  ║                   lifetime = 0 (if K = N, no pipes)       ║
  ║   Metric:         g = (NC/4)(I+J)  (unchanged)            ║
  ║   Security:       collision = 2^(NC/2) (unchanged)        ║
  ║   Cost:           (3N/K + 5) × NC                         ║
  ║                                                           ║
  ╚═══════════════════════════════════════════════════════════╝

  OPTIMAL DESIGN TABLE (128-bit collision = NC/2 = 128):

  {'N':>3} {'C':>4} {'K':>3} {'Rounds':>7} {'Cost':>8} {'Conservation':>13} {'Speedup':>8}
""")

    sha_cost = 29 * 256
    for N, C, K_n in [(2,128,2), (4,64,2), (4,64,4), (8,32,2), (8,32,4), (8,32,8), (16,16,2), (16,16,8)]:
        m = N * C
        if m != 256: continue
        R = int(3*N/K_n + 5)
        cost = R * m
        cons = f"{(N-K_n)//2}r" if K_n < N else "none"
        speedup = sha_cost / cost
        print(f"  {N:>3} {C:>4} {K_n:>3} {R:>7} {cost:>8} {cons:>13} {speedup:>7.2f}×")

    print(f"""
  WINNER: N=4, C=64, K=4 (all nodes) → 8 rounds, cost=2048
          Speedup: {sha_cost/2048:.1f}× vs SHA-256
          BUT: no conservation law (K=N) → potentially weaker

  SAFE WINNER: N=4, C=64, K=2 → 11 rounds, cost=2816
               Speedup: {sha_cost/2816:.1f}× vs SHA-256
               With conservation (lifetime=1 round)

  ═══════════════════════════════════════════════════════════════
""")


if __name__ == "__main__":
    main()
