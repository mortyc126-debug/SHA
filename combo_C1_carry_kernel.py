#!/usr/bin/env python3
"""
combo_C1_carry_kernel.py -- Carry-Algebra x Bilinear Rank-5 Kernel

Tests whether the carry identity D_add[L](x,y) = L(c) - 2(L(x) AND L(c))
simplifies or vanishes when perturbations are applied along kernel
directions of the bilinear form B_round = B_Ch + B_Maj.

Kernel directions (invisible to quadratic nonlinearity of Ch/Maj):
  - d alone        (index 3)
  - h alone        (index 7)
  - f XOR g        (indices 5,6 flipped simultaneously)
  - a XOR b XOR c  (indices 0,1,2 flipped simultaneously)

Non-kernel (control) directions:
  - a alone        (index 0)
  - e alone        (index 4)
  - f alone        (index 5)

References:
  - clifford_step3b.py: carry identity D_add[L](x,y) = L(c) - 2(L(x) AND L(c))
  - clifford_step1.py: bilinear form, rank-5 kernel {d, h, f^g, a^b^c}
"""

import random
import numpy as np
from collections import defaultdict

MASK = 0xFFFFFFFF
N = 10000

# ── SHA-256 primitives ──────────────────────────────────────────────

def rotr(x, n):
    return ((x >> n) | (x << (32 - n))) & MASK

def Sig0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sig1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g):  return ((e & f) ^ (~e & g)) & MASK
def Maj(a, b, c): return ((a & b) ^ (a & c) ^ (b & c)) & MASK
def add(x, y): return (x + y) & MASK
def sub(x, y): return (x - y) & MASK
def hw(x): return bin(x).count('1')

K0 = 0x428a2f98  # first round constant

def sha_one_round(state, W, K_r=K0):
    a, b, c, d, e, f, g, h = state
    T1 = add(add(add(add(h, Sig1(e)), Ch(e, f, g)), K_r), W)
    T2 = add(Sig0(a), Maj(a, b, c))
    return [add(T1, T2), a, b, c, add(d, T1), e, f, g]


# ── Perturbation definitions ────────────────────────────────────────

def perturb_state(state, name, delta):
    """Apply perturbation delta to state along direction 'name'.
    Returns perturbed state (XOR perturbation on the relevant registers)."""
    s = list(state)
    if name == "d":
        s[3] ^= delta
    elif name == "h":
        s[7] ^= delta
    elif name == "f^g":
        s[5] ^= delta
        s[6] ^= delta
    elif name == "a^b^c":
        s[0] ^= delta
        s[1] ^= delta
        s[2] ^= delta
    elif name == "a":
        s[0] ^= delta
    elif name == "e":
        s[4] ^= delta
    elif name == "f":
        s[5] ^= delta
    else:
        raise ValueError(f"Unknown direction: {name}")
    return s

KERNEL_DIRS    = ["d", "h", "f^g", "a^b^c"]
NONKERNEL_DIRS = ["a", "e", "f"]
ALL_DIRS = KERNEL_DIRS + NONKERNEL_DIRS

# ── Experiment ──────────────────────────────────────────────────────

print("=" * 72)
print("  C1: Carry-Algebra x Bilinear Rank-5 Kernel")
print("=" * 72)
print()
print(f"N = {N} random states, random 32-bit delta per trial")
print()

# Storage
hw_xor_diff   = {d: [] for d in ALL_DIRS}   # HW of XOR output diff
hw_add_diff   = {d: [] for d in ALL_DIRS}   # HW of additive output diff
carry_corr_hw = {d: [] for d in ALL_DIRS}   # HW of carry correction term

# Linearity test storage:
# D(x + d1 + d2) =? D(x + d1) XOR D(x + d2)
linearity_ok = {d: 0 for d in ALL_DIRS}
linearity_total = {d: 0 for d in ALL_DIRS}

random.seed(42)

for trial in range(N):
    state = [random.randint(0, MASK) for _ in range(8)]
    W = random.randint(0, MASK)
    out_orig = sha_one_round(state, W)

    delta = random.randint(1, MASK)  # nonzero random 32-bit

    for dname in ALL_DIRS:
        # ── (A) Single-perturbation output difference ────────────
        state_p = perturb_state(state, dname, delta)
        out_p = sha_one_round(state_p, W)

        # XOR difference
        xdiff = [out_p[i] ^ out_orig[i] for i in range(8)]
        hw_xor_diff[dname].append(sum(hw(x) for x in xdiff))

        # Additive difference (mod 2^32 per register)
        adiff = [sub(out_p[i], out_orig[i]) for i in range(8)]
        hw_add_diff[dname].append(sum(hw(x) for x in adiff))

        # ── (B) Carry correction magnitude ───────────────────────
        # For the registers that participate in modular addition in T1/T2,
        # measure the carry correction 2*(L(x) AND L(c)) for the Sigma
        # functions applied to the perturbed registers.
        #
        # We measure the full-round carry effect: for each output register,
        # compare additive diff to XOR diff. The "carry correction" is:
        #   additive_out - xor_out (when interpreted as integer diff vs XOR diff)
        # But more directly: for each output word, carry_correction =
        #   out_p[i] - out_orig[i] XOR (out_p[i] XOR out_orig[i])
        # i.e.  sub(out_p[i], out_orig[i]) XOR xdiff[i]
        # which captures how much the additive and XOR views disagree.
        total_carry_hw = 0
        for i in range(8):
            # The "carry correction" in the identity:
            # D_add = L(c) - 2(L(x) AND L(c))
            # The 2(L(x) AND L(c)) part is the difference between L(c) and D_add
            # We measure it as: L(c) - D_add = 2(L(x) AND L(c))
            # where L(c) would be the XOR-predicted diff if everything were linear.
            # For the full round this is an approximation; we measure the
            # discrepancy between xor diff and additive diff as carry magnitude.
            carry_term = adiff[i] ^ xdiff[i]
            total_carry_hw += hw(carry_term)
        carry_corr_hw[dname].append(total_carry_hw)

    # ── (C) Linearity test (every 5th trial to keep runtime manageable) ──
    if trial % 5 == 0:
        delta1 = random.randint(1, MASK)
        delta2 = random.randint(1, MASK)

        for dname in ALL_DIRS:
            s1 = perturb_state(state, dname, delta1)
            s2 = perturb_state(state, dname, delta2)
            s12 = perturb_state(perturb_state(state, dname, delta1), dname, delta2)
            # Note: for XOR perturbations, s12 = perturb(state, dname, delta1 ^ delta2)
            # but we do it this way to be explicit.

            o1 = sha_one_round(s1, W)
            o2 = sha_one_round(s2, W)
            o12 = sha_one_round(s12, W)

            # D(x+d1) XOR D(x+d2) vs D(x+d1+d2)
            # where D(x+d) = Round(x+d) XOR Round(x)
            diff1 = [o1[i] ^ out_orig[i] for i in range(8)]
            diff2 = [o2[i] ^ out_orig[i] for i in range(8)]
            diff12 = [o12[i] ^ out_orig[i] for i in range(8)]

            # Linearity check: diff12 == diff1 XOR diff2 ?
            predicted = [diff1[i] ^ diff2[i] for i in range(8)]
            if predicted == diff12:
                linearity_ok[dname] += 1
            linearity_total[dname] += 1


# ── Print results ───────────────────────────────────────────────────

print("-" * 72)
print("METRIC 1: Average Hamming Weight of XOR output difference")
print("-" * 72)
print(f"  {'Direction':<12} {'Mean HW':>10} {'Std':>10}   Type")
print()
for dname in ALL_DIRS:
    arr = hw_xor_diff[dname]
    tag = "KERNEL" if dname in KERNEL_DIRS else "non-kernel"
    print(f"  {dname:<12} {np.mean(arr):10.2f} {np.std(arr):10.2f}   {tag}")

kern_mean_xor = np.mean([np.mean(hw_xor_diff[d]) for d in KERNEL_DIRS])
nonk_mean_xor = np.mean([np.mean(hw_xor_diff[d]) for d in NONKERNEL_DIRS])
print()
print(f"  Kernel avg:     {kern_mean_xor:.2f}")
print(f"  Non-kernel avg: {nonk_mean_xor:.2f}")
print(f"  Ratio (non-kernel / kernel): {nonk_mean_xor / kern_mean_xor:.4f}")

print()
print("-" * 72)
print("METRIC 2: Average Hamming Weight of additive output difference")
print("-" * 72)
print(f"  {'Direction':<12} {'Mean HW':>10} {'Std':>10}   Type")
print()
for dname in ALL_DIRS:
    arr = hw_add_diff[dname]
    tag = "KERNEL" if dname in KERNEL_DIRS else "non-kernel"
    print(f"  {dname:<12} {np.mean(arr):10.2f} {np.std(arr):10.2f}   {tag}")

kern_mean_add = np.mean([np.mean(hw_add_diff[d]) for d in KERNEL_DIRS])
nonk_mean_add = np.mean([np.mean(hw_add_diff[d]) for d in NONKERNEL_DIRS])
print()
print(f"  Kernel avg:     {kern_mean_add:.2f}")
print(f"  Non-kernel avg: {nonk_mean_add:.2f}")
print(f"  Ratio (non-kernel / kernel): {nonk_mean_add / kern_mean_add:.4f}")

print()
print("-" * 72)
print("METRIC 3: Carry correction magnitude |2(L(x) AND L(c))|")
print("-" * 72)
print(f"  {'Direction':<12} {'Mean HW':>10} {'Std':>10}   Type")
print()
for dname in ALL_DIRS:
    arr = carry_corr_hw[dname]
    tag = "KERNEL" if dname in KERNEL_DIRS else "non-kernel"
    print(f"  {dname:<12} {np.mean(arr):10.2f} {np.std(arr):10.2f}   {tag}")

kern_mean_carry = np.mean([np.mean(carry_corr_hw[d]) for d in KERNEL_DIRS])
nonk_mean_carry = np.mean([np.mean(carry_corr_hw[d]) for d in NONKERNEL_DIRS])
print()
print(f"  Kernel avg:     {kern_mean_carry:.2f}")
print(f"  Non-kernel avg: {nonk_mean_carry:.2f}")
if kern_mean_carry > 0:
    print(f"  Ratio (non-kernel / kernel): {nonk_mean_carry / kern_mean_carry:.4f}")
else:
    print(f"  Kernel carry correction is ZERO!")

print()
print("-" * 72)
print("METRIC 4: XOR-linearity test  D(x+d1+d2) =? D(x+d1) XOR D(x+d2)")
print("-" * 72)
print(f"  {'Direction':<12} {'Pass':>8} / {'Total':>5}  {'Rate':>8}   Type")
print()
for dname in ALL_DIRS:
    ok = linearity_ok[dname]
    tot = linearity_total[dname]
    rate = ok / tot if tot > 0 else 0
    tag = "KERNEL" if dname in KERNEL_DIRS else "non-kernel"
    print(f"  {dname:<12} {ok:8d} / {tot:5d}  {rate:8.4f}   {tag}")

kern_lin = np.mean([linearity_ok[d] / linearity_total[d] for d in KERNEL_DIRS])
nonk_lin = np.mean([linearity_ok[d] / linearity_total[d] for d in NONKERNEL_DIRS])
print()
print(f"  Kernel linearity rate:     {kern_lin:.4f}")
print(f"  Non-kernel linearity rate: {nonk_lin:.4f}")

# ── Per-direction breakdown: which output registers are affected ────

print()
print("-" * 72)
print("METRIC 5: Per-register XOR diff HW (averaged over all trials)")
print("-" * 72)
print(f"  {'Direction':<12}", end="")
for reg in "abcdefgh":
    print(f" {reg:>7}", end="")
print("   Type")
print()

# Recompute per-register stats with a smaller sample for display
random.seed(42)
per_reg_hw = {d: [[] for _ in range(8)] for d in ALL_DIRS}
for trial in range(min(N, 5000)):
    state = [random.randint(0, MASK) for _ in range(8)]
    W = random.randint(0, MASK)
    out_orig = sha_one_round(state, W)
    delta = random.randint(1, MASK)
    for dname in ALL_DIRS:
        state_p = perturb_state(state, dname, delta)
        out_p = sha_one_round(state_p, W)
        for i in range(8):
            per_reg_hw[dname][i].append(hw(out_p[i] ^ out_orig[i]))

for dname in ALL_DIRS:
    tag = "KERNEL" if dname in KERNEL_DIRS else "non-kernel"
    print(f"  {dname:<12}", end="")
    for i in range(8):
        print(f" {np.mean(per_reg_hw[dname][i]):7.2f}", end="")
    print(f"   {tag}")

# ── Verdict ─────────────────────────────────────────────────────────

print()
print("=" * 72)
print("  VERDICT")
print("=" * 72)
print()

# Gather effect sizes
xor_effect = abs(kern_mean_xor - nonk_mean_xor) / max(kern_mean_xor, nonk_mean_xor)
add_effect = abs(kern_mean_add - nonk_mean_add) / max(kern_mean_add, nonk_mean_add)
carry_effect = abs(kern_mean_carry - nonk_mean_carry) / max(kern_mean_carry, nonk_mean_carry) if max(kern_mean_carry, nonk_mean_carry) > 0 else 0
lin_effect = abs(kern_lin - nonk_lin)

print(f"  XOR diff effect size (relative):     {xor_effect:.4f}")
print(f"  Add diff effect size (relative):     {add_effect:.4f}")
print(f"  Carry correction effect size (rel):  {carry_effect:.4f}")
print(f"  Linearity rate difference (abs):     {lin_effect:.4f}")
print()

# Thresholds for significance
# If kernel directions show >5% relative difference in any metric, ANOMALY
# If >15%, ALIVE. If <2% on all, DEAD.

max_effect = max(xor_effect, add_effect, carry_effect)
sig_linearity = lin_effect > 0.01  # linearity rates differ by >1%

if max_effect > 0.15 or (lin_effect > 0.05):
    verdict = "ALIVE"
    reasoning = (
        "Kernel directions show SIGNIFICANT differences from non-kernel.\n"
        "  The carry correction and/or output diffusion differs substantially\n"
        "  along kernel vs non-kernel directions. The quadratic-invisible\n"
        "  subspace interacts non-trivially with the carry algebra."
    )
elif max_effect > 0.05 or sig_linearity:
    verdict = "ANOMALY"
    reasoning = (
        "Kernel directions show MODERATE differences from non-kernel.\n"
        "  There is a measurable but modest effect of kernel structure\n"
        "  on carry propagation behavior."
    )
else:
    verdict = "DEAD"
    reasoning = (
        "Kernel directions show NO significant differences from non-kernel.\n"
        "  The carry correction term 2(L(x) AND L(c)) does NOT simplify\n"
        "  along kernel directions. The carry algebra overwhelms the\n"
        "  bilinear kernel structure."
    )

print(f"  *** {verdict} ***")
print()
print(f"  {reasoning}")
print()

# Additional context
print("-" * 72)
print("INTERPRETATION")
print("-" * 72)
print()
print("The kernel {d, h, f^g, a^b^c} is invisible to Ch/Maj quadratic forms.")
print("But the SHA-256 round also has MODULAR ADDITION which introduces")
print("its own nonlinearity via carry propagation.")
print()
print("If ALIVE:  Kernel directions reduce BOTH quadratic AND carry")
print("           nonlinearity -- the round function is closer to linear")
print("           in the kernel subspace. This is exploitable.")
print()
print("If DEAD:   Carry nonlinearity is independent of the bilinear kernel.")
print("           The kernel structure helps with Ch/Maj but not with")
print("           modular addition. Exploitation requires handling both")
print("           sources of nonlinearity separately.")
print()
print("If ANOMALY: Partial interaction. Some metrics show kernel advantage,")
print("            others do not. Further investigation needed.")
print()
