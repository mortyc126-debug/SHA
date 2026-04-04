"""
PROVABLE SECURITY DIRECTION: Can BTE theory prove lower bounds?

Sketch:
  T6: After R_H rounds, Hessian = random (0.5).
  → Function has no degree-2 structure.
  → Any algebraic attack needs degree ≥ 3.
  → XL at degree 3: cost ≥ C(N, 3) ≈ N^3 / 6.
  → For N = 512: cost ≥ 2^{27} (far below birthday).

This doesn't prove birthday hardness. But it proves:
  "No degree-2 algebraic attack works after R_H rounds."

Can we strengthen to: "No degree-d attack for any d"?
  From Hessian profile: D_k (k-th derivative) → random by round ~R_H + (k-2)×Δ.
  D_2 random at R_H ≈ 12.
  D_3 random at R ≈ 12-16 (measured: 0.47 at R=64, but what about R=12?).
  D_4 random at R ≈ ?

If ALL derivatives random by R_H + constant → function = random →
birthday bound follows.

Let's measure D_3 and D_4 at REDUCED rounds to find their transition points.
"""

import random
from stage0_observe import sha256_round_trace, MASK, H0, K, Sig0, Sig1, Ch, Maj


def sha256_hash_R(M, R):
    states, _ = sha256_round_trace(M, rounds=R)
    return [(states[R][i] + H0[i]) & MASK for i in range(8)]


def derivative_profile():
    """
    Measure D_2, D_3, D_4 at each round count.
    Find transition point for each.
    """
    print("=" * 80)
    print("DERIVATIVE TRANSITION PROFILE: D_2, D_3, D_4 vs rounds")
    print("=" * 80)

    rng = random.Random(42)
    M = [rng.randint(0, MASK) for _ in range(16)]

    print(f"\n  {'R':>3} | {'D2':>6} {'D3':>6} {'D4':>6} | {'D2 random?':>10} {'D3 random?':>10}")
    print("-" * 60)

    for R in [2, 4, 6, 8, 10, 12, 14, 16, 20, 24, 32, 64]:
        def f(msg):
            return sha256_hash_R(msg, R)[0] & 1

        N = 300
        d2_nz = 0; d3_nz = 0; d4_nz = 0; total = 0

        for trial in range(N):
            r = random.Random(trial * 1000 + R)
            flips = [(r.randint(0, 15), r.randint(0, 31)) for _ in range(4)]
            if len(set(flips)) < 4:
                continue

            def flip_M(M, active_flips):
                M2 = list(M)
                for w, b in active_flips:
                    M2[w] ^= (1 << b)
                return M2

            # All subsets
            vals = {}
            for mask_val in range(16):
                active = [flips[i] for i in range(4) if (mask_val >> i) & 1]
                vals[mask_val] = f(flip_M(M, active))

            # D2: subset size 2
            d2 = vals[0] ^ vals[1] ^ vals[2] ^ vals[3]  # {}, {0}, {1}, {0,1}
            # D3: subset size 3
            d3 = 0
            for s in range(8):
                d3 ^= vals[s]  # all subsets of {0,1,2}
            # D4: all 16 subsets
            d4 = 0
            for s in range(16):
                d4 ^= vals[s]

            if d2: d2_nz += 1
            if d3: d3_nz += 1
            if d4: d4_nz += 1
            total += 1

        p2 = d2_nz / total if total > 0 else 0
        p3 = d3_nz / total if total > 0 else 0
        p4 = d4_nz / total if total > 0 else 0

        r2 = "random" if abs(p2 - 0.5) < 0.1 else f"struct({p2:.2f})"
        r3 = "random" if abs(p3 - 0.5) < 0.1 else f"struct({p3:.2f})"

        print(f"  {R:>3} | {p2:>6.3f} {p3:>6.3f} {p4:>6.3f} | {r2:>10} {r3:>10}")

    print(f"""
    INTERPRETATION:
    D_k transition round = first R where D_k ≈ 0.5 (random).

    If D_2 transition ≈ D_3 transition ≈ D_4 transition ≈ R_H:
      → ALL derivatives random simultaneously at R_H
      → Function is RANDOM beyond R_H
      → This would PROVE birthday hardness for R > R_H

    If D_k transitions are DIFFERENT:
      → Higher-order structure persists longer
      → Need more rounds for higher-order randomness
      → Security depends on the HIGHEST-ORDER transition
    """)


if __name__ == "__main__":
    derivative_profile()
