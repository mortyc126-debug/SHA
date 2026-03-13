#!/usr/bin/env python3
"""
jacobian_gf2.py — GF(2) Jacobian analysis of the SHA-256 differential system.

Two experiments for each W0_seed:

  A) Full search: DW ∈ {0,1}^16 (2^16 = 65536 patterns)
     — measures P(non-trivial Sol_1 exists in full compact space)
     — computes 15×16 Jacobian J_A at DW=0

  B) Birthday slice: DW[0]=1 fixed, DW[1..15] ∈ {0,1}^15 (2^15 = 32768 patterns)
     — this is the context of the "9% phenomenon" from the methodology
     — computes 15×15 Jacobian J_B around reference DW=(1,0,...,0)
     — tests hypothesis: P(Sol_1 ≠ ∅) = P(F_ref ∈ Image(J_B)) = 2^{rank(J_B) - 15}

Theory for experiment B:
  F: {0,1}^15 → {0,1}^15,  F(x)_r = De_r(DW[0]=1, DW[1..15]=x) mod 2
  Random fn: E[|F^{-1}(0)|] = 1  →  P(sol) ≈ 1 - 1/e ≈ 63%
  SHA-256: P(sol) ≈ 9%  (from П-47..П-55 methodology)
  Hypothesis: rank(J_B) ≈ 11-12
    rank=12 → P = 2^{12-15} = 12.5%
    rank=11 → P = 2^{11-15} = 6.25%

Usage:
  python3 jacobian_gf2.py [n_seeds]    (default: 200)
"""

import numpy as np
import sys
import time

# ── SHA-256 constants (rounds 1..17 only) ─────────────────────────────────────
K17 = np.array([
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b,
    0x59f111f1, 0x923f82a4, 0xab1c5ed5, 0xd807aa98, 0x12835b01,
    0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7,
    0xc19bf174, 0xe49b69c1,
], dtype=np.uint32)

IV = np.array([
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
], dtype=np.uint32)


# ── SHA-256 primitives (numpy-vectorised) ─────────────────────────────────────
def R(x, n):
    n = np.uint32(n)
    return (x >> n) | (x << (np.uint32(32) - n))

SIG0 = lambda x: R(x, 2)  ^ R(x, 13) ^ R(x, 22)
SIG1 = lambda x: R(x, 6)  ^ R(x, 11) ^ R(x, 25)
sig0 = lambda x: R(x, 7)  ^ R(x, 18) ^ (x >> np.uint32(3))
sig1 = lambda x: R(x, 17) ^ R(x, 19) ^ (x >> np.uint32(10))
CH   = lambda e, f, g: (e & f) ^ (~e & g)
MAJ  = lambda a, b, c: (a & b) ^ (a & c) ^ (b & c)


def eval_de_batch(W0_seed: int, DW: np.ndarray) -> np.ndarray:
    """
    Evaluate De_r = e2_r - e1_r  for r = 3..17 over a batch of DW vectors.

    Parameters
    ----------
    W0_seed : int
        Base message W[0].  W[1..15] = 0.
        Second message: W'[0] = W0_seed + DW[:,0],  W'[i>0] = DW[:,i].
    DW : (N, 16) uint32

    Returns
    -------
    (N, 15) uint32  —  De_r values (mod 2^32), columns r=3..17
    """
    N  = DW.shape[0]
    W0 = np.uint32(W0_seed)

    def fill(v): return np.full(N, np.uint32(v), dtype=np.uint32)

    a1, b1, c1, d1, e1, f1, g1, h1 = [fill(v) for v in IV]
    a2, b2, c2, d2, e2, f2, g2, h2 = [fill(v) for v in IV]

    def step(a, b, c, d, e, f, g, h, w, K):
        T1 = h + SIG1(e) + CH(e, f, g) + K + w
        T2 = SIG0(a) + MAJ(a, b, c)
        return T1 + T2, a, b, c, d + T1, e, f, g   # new a, b, c, d, e, f, g, h

    # Round 1
    a1,b1,c1,d1,e1,f1,g1,h1 = step(a1,b1,c1,d1,e1,f1,g1,h1, fill(W0),         K17[0])
    a2,b2,c2,d2,e2,f2,g2,h2 = step(a2,b2,c2,d2,e2,f2,g2,h2, W0 + DW[:,0],     K17[0])

    # Round 2
    a1,b1,c1,d1,e1,f1,g1,h1 = step(a1,b1,c1,d1,e1,f1,g1,h1, fill(0),   K17[1])
    a2,b2,c2,d2,e2,f2,g2,h2 = step(a2,b2,c2,d2,e2,f2,g2,h2, DW[:,1],   K17[1])

    out = np.zeros((N, 15), dtype=np.uint32)
    for r in range(3, 18):
        if r <= 16:
            w1 = fill(0)
            w2 = DW[:, r - 1]
        else:   # r = 17: W[16] from message schedule
            # W[16]  = sig1(0)     + 0     + sig0(0)      + W0  = W0
            # W'[16] = sig1(DW14)  + DW9   + sig0(DW1)    + W0 + DW0
            w1 = fill(W0)
            w2 = W0 + DW[:,0] + sig0(DW[:,1]) + DW[:,9] + sig1(DW[:,14])

        a1,b1,c1,d1,e1,f1,g1,h1 = step(a1,b1,c1,d1,e1,f1,g1,h1, w1, K17[r-1])
        a2,b2,c2,d2,e2,f2,g2,h2 = step(a2,b2,c2,d2,e2,f2,g2,h2, w2, K17[r-1])

        out[:, r - 3] = e2 - e1   # mod 2^32 automatic in uint32

    return out


def parity_mat(De: np.ndarray) -> np.ndarray:
    """(N,15) uint32  →  (N,15) uint8, 1 iff De_r is odd."""
    return (De & np.uint32(1)).astype(np.uint8)


def rank_gf2(M: np.ndarray) -> int:
    """Rank of matrix M over GF(2) via Gaussian elimination."""
    A = np.asarray(M, dtype=np.uint8) & 1
    rows, cols = A.shape
    r = 0
    for col in range(cols):
        idx = np.where(A[r:, col])[0]
        if not len(idx):
            continue
        p = r + idx[0]
        if p != r:
            A[[r, p]] = A[[p, r]]
        mask = A[:, col].astype(bool)
        mask[r] = False
        A[mask] ^= A[r]
        r += 1
        if r == rows:
            break
    return r


def in_image_gf2(J: np.ndarray, b: np.ndarray) -> bool:
    """True iff b ∈ Image(J) over GF(2).  J: (m,n), b: (m,)."""
    J_aug = np.hstack([J, b.reshape(-1, 1)]).astype(np.uint8) & 1
    return rank_gf2(J_aug) == rank_gf2(J)


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    N_seeds = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    print(f"=== SHA-256 GF(2) Jacobian / Sol_1 analysis  (N_seeds={N_seeds}) ===\n")

    # ── Precompute DW pattern arrays ──────────────────────────────────────────

    # A: full {0,1}^16
    n_A   = 1 << 16
    idx_A = np.arange(n_A, dtype=np.uint32)
    DW_A  = np.zeros((n_A, 16), dtype=np.uint32)
    for i in range(16):
        DW_A[:, i] = (idx_A >> np.uint32(i)) & 1

    # B: birthday slice — DW[0]=1 fixed, DW[1..15] ∈ {0,1}^15
    n_B   = 1 << 15
    idx_B = np.arange(n_B, dtype=np.uint32)
    DW_B  = np.zeros((n_B, 16), dtype=np.uint32)
    DW_B[:, 0] = 1
    for i in range(1, 16):
        DW_B[:, i] = (idx_B >> np.uint32(i - 1)) & 1

    # Jacobian sample indices
    # A: e_i^(A) = DW_A[2^i]   (unit vector in direction i, DW[0]=0)
    e_idx_A = [1 << i for i in range(16)]   # indices into DW_A
    # B: at reference DW=(1,0,...,0) = DW_B[0]
    #    e_j^(B) = DW_B[2^{j-1}]  for j=1..15  (flip DW[j] relative to ref)
    e_idx_B = [1 << (j - 1) for j in range(1, 16)]   # indices into DW_B

    rng   = np.random.default_rng(0xC0DE_2026_03_13)
    seeds = rng.integers(0, 1 << 32, N_seeds).astype(np.uint32)

    # ── Accumulators ──────────────────────────────────────────────────────────
    A_nt        = 0               # non-trivial Sol_1 found (exp A)
    A_sol_sizes = []
    A_ranks     = {}

    B_sol       = 0               # any Sol_1 found (exp B)
    B_img       = 0               # F_ref ∈ Image(J_B)
    B_ranks     = {}
    # per-rank: [seeds, sol_count, img_count]
    B_rank_detail = {}

    t0 = time.time()
    for j, seed in enumerate(seeds):
        W0 = int(seed)

        # ══════════════ EXPERIMENT A ══════════════
        De_A  = eval_de_batch(W0, DW_A)               # (65536, 15) uint32
        par_A = parity_mat(De_A)                       # (65536, 15) uint8; 1=odd
        # Sol_1: all parities = 0 (all De_r even)
        sol_A = (par_A == 0).all(axis=1)              # (65536,) bool
        n_sol = int(sol_A.sum())
        has_nt = bool(sol_A[1:].any())                # exclude trivial DW=0

        # J_A[r-3, i] = par_A[e_i, r-3]  (parity of De_r under unit difference)
        J_A  = par_A[e_idx_A, :].T                    # (15, 16) uint8
        rk_A = rank_gf2(J_A)

        A_nt += has_nt
        A_sol_sizes.append(n_sol)
        A_ranks[rk_A] = A_ranks.get(rk_A, 0) + 1

        # ══════════════ EXPERIMENT B ══════════════
        De_B  = eval_de_batch(W0, DW_B)               # (32768, 15) uint32
        par_B = parity_mat(De_B)                       # (32768, 15) uint8
        sol_B = (par_B == 0).all(axis=1)              # True iff all De_r even
        has_sol_B = bool(sol_B.any())

        # Reference point: DW_B[0] = (1,0,...,0)
        F_ref = par_B[0]                              # (15,) parity at ref

        # J_B[r-3, j-1] = par_B[2^{j-1}, r-3]  XOR  F_ref[r-3]
        # = effect of flipping DW[j] from 0→1 on parity of De_r
        J_B_cols = [par_B[e_idx_B[j]] ^ F_ref for j in range(15)]
        J_B  = np.array(J_B_cols).T                   # (15, 15) uint8
        rk_B = rank_gf2(J_B)

        img_B = in_image_gf2(J_B, F_ref)             # F_ref ∈ Image(J_B)?

        B_sol += has_sol_B
        B_img += img_B
        B_ranks[rk_B] = B_ranks.get(rk_B, 0) + 1
        if rk_B not in B_rank_detail:
            B_rank_detail[rk_B] = [0, 0, 0]  # total, sol, img
        B_rank_detail[rk_B][0] += 1
        B_rank_detail[rk_B][1] += has_sol_B
        B_rank_detail[rk_B][2] += img_B

        # progress
        if (j + 1) % max(1, N_seeds // 10) == 0 or j == N_seeds - 1:
            pA  = 100 * A_nt    / (j + 1)
            pB  = 100 * B_sol   / (j + 1)
            pI  = 100 * B_img   / (j + 1)
            print(f"  [{j+1:5d}/{N_seeds}]  "
                  f"A: P(nt-sol)={pA:.1f}%  |"
                  f"  B: P(sol)={pB:.1f}%  P(F∈Im)={pI:.1f}%  "
                  f"ranks_B={dict(sorted(B_ranks.items()))}  "
                  f"t={time.time()-t0:.0f}s", flush=True)

    N = N_seeds
    print()
    print("=" * 72)
    print()

    # ── Report A ──────────────────────────────────────────────────────────────
    print("─── EXPERIMENT A: full DW ∈ {0,1}^16  (65536 patterns/seed) ───")
    print(f"  P(non-trivial Sol_1 ≠ ∅) = {A_nt}/{N} = {100*A_nt/N:.1f}%")
    print(f"  Expected random function  ≈ 86%   (2^16 inputs, 2^15 outputs, E[nt]≈2)")
    print(f"  Avg |Sol_1| (incl. trivial DW=0) = {sum(A_sol_sizes)/N:.3f}")
    print(f"  Expected random ≈ 2.000")
    print()
    print(f"  Rank of J_A (15×16) over GF(2):")
    for rk in sorted(A_ranks):
        ker = 1 << (16 - rk)
        print(f"    rank={rk:2d}:  {A_ranks[rk]:5d} seeds  ({100*A_ranks[rk]/N:4.1f}%)"
              f"  |ker(J_A)|=2^{16-rk}={ker}")

    print()
    print("─── EXPERIMENT B: DW[0]=1 fixed, DW[1..15] ∈ {0,1}^15  (32768/seed) ───")
    print(f"  P(Sol_1 ≠ ∅)             = {B_sol}/{N} = {100*B_sol/N:.1f}%")
    print(f"  Expected random function  ≈ 63%   (2^15 in, 2^15 out, E[sol]=1)")
    print(f"  Expected SHA-256          ≈  9%   (from methodology П-47..П-55)")
    print()
    print(f"  P(F_ref ∈ Image(J_B))    = {B_img}/{N} = {100*B_img/N:.1f}%")
    print(f"  [linear prediction; should ≈ P(Sol_1) if F is well-linearised]")
    print()
    print(f"  Rank of J_B (15×15) over GF(2):")
    print(f"  {'rank':>4}  {'seeds':>5}  {'%':>4}  "
          f"{'P(sol|rank)':>11}  {'P(F∈Im|rank)':>12}  "
          f"{'P_theory=2^(r-15)':>18}")
    for rk in sorted(B_rank_detail):
        tot, sol, img = B_rank_detail[rk]
        p_sol = 100 * sol / tot if tot else 0
        p_img = 100 * img / tot if tot else 0
        p_th  = 100 * (2 ** (rk - 15))
        print(f"  {rk:>4}  {tot:>5}  {100*tot/N:>3.0f}%  "
              f"{p_sol:>10.1f}%  {p_img:>11.1f}%  "
              f"{p_th:>10.2f}%")

    print()
    print("─── KEY CONCLUSION ───────────────────────────────────────────────────")
    print("Hypothesis: rank(J_B) explains the 9% phenomenon.")
    print("  rank=15 → P = 2^0   = 100% (fully solvable, always has solution)")
    print("  rank=14 → P = 2^-1  =  50%")
    print("  rank=13 → P = 2^-2  =  25%")
    print("  rank=12 → P = 2^-3  = 12.5%   ← matches ~9%?")
    print("  rank=11 → P = 2^-4  =  6.25%  ← matches ~9%?")
    print("  rank=10 → P = 2^-5  =  3.1%")
    print()
    print("If rank distribution is concentrated at 11-12 and")
    print("P(F∈Im|rank) ≈ 2^{rank-15}, the hypothesis is CONFIRMED.")

    elapsed = time.time() - t0
    print(f"\nTotal time: {elapsed:.1f}s  ({elapsed/N:.2f}s/seed)")


if __name__ == "__main__":
    main()
