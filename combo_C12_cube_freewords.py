#!/usr/bin/env python3
"""
C12: Cube Attack × Free Words

Key insight: Free words W[12..15] (128 bits) don't affect barrier Da13
but DO affect De17 (the first barrier equation).

We apply Dinur-Shamir cube attacks to De17 viewed as a function of W[12..15].
If the superpoly has degree ≤ 3, we can solve barriers faster than brute-force 2^32.

Experiment:
- Fix W[0..11], DW[0] = 0x80000000
- De17 is a function of W[12], W[13], W[14], W[15]
- Apply cube attacks with small cubes (k=1..6 bits) from W[14] and cross-word cubes
- Test superpoly constancy and linearity
- Control: verify De13 does NOT depend on W[12..15]
"""

import random
import itertools

MASK = 0xFFFFFFFF

def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK
def Sig0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sig1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def sig0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sig1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)
def Ch(e, f, g):  return ((e & f) ^ (~e & g)) & MASK
def Maj(a, b, c): return ((a & b) ^ (a & c) ^ (b & c)) & MASK
def add(*args):
    s = 0
    for x in args:
        s = (s + x) & MASK
    return s

IV = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]
K = [0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
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
     0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2]

DW_BIT = 0x80000000

def schedule(W16):
    W = list(W16) + [0] * 48
    for i in range(16, 64):
        W[i] = add(sig1(W[i-2]), W[i-7], sig0(W[i-15]), W[i-16])
    return W

def sha_round_fn(state, W_r, K_r):
    a, b, c, d, e, f, g, h = state
    T1 = add(h, Sig1(e), Ch(e, f, g), K_r, W_r)
    T2 = add(Sig0(a), Maj(a, b, c))
    return [add(T1, T2), a, b, c, add(d, T1), e, f, g]

def compute_De(W16, num_rounds):
    """Compute De at given round: differential in e register."""
    W64 = schedule(W16)
    W16_f = list(W16)
    W16_f[0] ^= DW_BIT
    W64_f = schedule(W16_f)
    s_n = list(IV)
    s_f = list(IV)
    for r in range(num_rounds):
        s_n = sha_round_fn(s_n, W64[r], K[r])
        s_f = sha_round_fn(s_f, W64_f[r], K[r])
    return s_n[4] ^ s_f[4]  # e register differential

def compute_Da(W16, num_rounds):
    """Compute Da at given round: differential in a register."""
    W64 = schedule(W16)
    W16_f = list(W16)
    W16_f[0] ^= DW_BIT
    W64_f = schedule(W16_f)
    s_n = list(IV)
    s_f = list(IV)
    for r in range(num_rounds):
        s_n = sha_round_fn(s_n, W64[r], K[r])
        s_f = sha_round_fn(s_f, W64_f[r], K[r])
    return s_n[0] ^ s_f[0]  # a register differential

# ============================================================
# Fix W[0..11] randomly; W[12..15] are the free variables
# ============================================================
random.seed(0xC12C0BE)
W_fixed = [random.randint(0, MASK) for _ in range(16)]

def make_W16(w12, w13, w14, w15):
    W = list(W_fixed)
    W[12] = w12
    W[13] = w13
    W[14] = w14
    W[15] = w15
    return W

# ============================================================
# PART 0: Control test — verify De13 does NOT depend on W[12..15]
# ============================================================
print("=" * 72)
print("C12: CUBE ATTACK x FREE WORDS")
print("=" * 72)
print()
print("-" * 72)
print("PART 0: Control — De13 independence from W[12..15]")
print("-" * 72)

base_w12 = random.randint(0, MASK)
base_w13 = random.randint(0, MASK)
base_w14 = random.randint(0, MASK)
base_w15 = random.randint(0, MASK)

W_base = make_W16(base_w12, base_w13, base_w14, base_w15)
de13_base = compute_De(W_base, 13)
da13_base = compute_Da(W_base, 13)

n_indep_tests = 200
de13_varies = False
da13_varies = False
for _ in range(n_indep_tests):
    rw12 = random.randint(0, MASK)
    rw13 = random.randint(0, MASK)
    rw14 = random.randint(0, MASK)
    rw15 = random.randint(0, MASK)
    W_test = make_W16(rw12, rw13, rw14, rw15)
    de13_test = compute_De(W_test, 13)
    da13_test = compute_Da(W_test, 13)
    if de13_test != de13_base:
        de13_varies = True
    if da13_test != da13_base:
        da13_varies = True

print(f"  Da13 varies with W[12..15]? {da13_varies}")
print(f"  De13 varies with W[12..15]? {de13_varies}")
if not da13_varies:
    print("  -> CONFIRMED: Da13 independent of free words (expected)")
else:
    print("  -> ANOMALY: Da13 depends on free words!")
if de13_varies:
    print("  -> De13 depends on free words — checking which round it starts...")
    for rnd in range(13, 18):
        varies = False
        de_base = compute_De(W_base, rnd)
        for _ in range(50):
            W_t = make_W16(random.randint(0, MASK), random.randint(0, MASK),
                           random.randint(0, MASK), random.randint(0, MASK))
            if compute_De(W_t, rnd) != de_base:
                varies = True
                break
        print(f"    De{rnd} varies? {varies}")
print()

# ============================================================
# PART 1: Single-word cubes from W[14] on De17
# ============================================================
print("-" * 72)
print("PART 1: Cube attack on De17 — cubes from W[14]")
print("-" * 72)
print()

test_bits = [0, 8, 16, 24]  # representative output bits of De17

def eval_de17_bit(w12, w13, w14, w15, bit_j):
    """Evaluate single bit j of De17."""
    W16 = make_W16(w12, w13, w14, w15)
    de17 = compute_De(W16, 17)
    return (de17 >> bit_j) & 1

def cube_sum_w14(cube_positions, bit_j, w12, w13, w14_rest, w15):
    """
    Sum De17[bit_j] over all 2^k assignments to cube_positions in W[14].
    cube_positions: list of bit indices in W[14] to vary
    w14_rest: base value of W[14] (cube positions will be overwritten)
    Returns: XOR sum (0 or 1)
    """
    k = len(cube_positions)
    # Build mask for cube bits
    cube_mask = 0
    for p in cube_positions:
        cube_mask |= (1 << p)

    total = 0
    for assignment in range(1 << k):
        # Set cube bits according to assignment
        w14 = w14_rest & ~cube_mask  # clear cube bits
        for idx, p in enumerate(cube_positions):
            if assignment & (1 << idx):
                w14 |= (1 << p)
        total ^= eval_de17_bit(w12, w13, w14, w15, bit_j)
    return total

N_RANDOM = 100  # random points to evaluate superpoly

for bit_j in test_bits:
    print(f"  Output bit j={bit_j}:")
    for k in range(1, 7):
        # Choose random cube positions from W[14]
        cube_pos = sorted(random.sample(range(32), k))

        # Evaluate superpoly at N_RANDOM random points
        superpoly_vals = []
        for _ in range(N_RANDOM):
            rw12 = random.randint(0, MASK)
            rw13 = random.randint(0, MASK)
            rw14 = random.randint(0, MASK)
            rw15 = random.randint(0, MASK)
            sp_val = cube_sum_w14(cube_pos, bit_j, rw12, rw13, rw14, rw15)
            superpoly_vals.append(sp_val)

        unique_vals = set(superpoly_vals)
        n_ones = sum(superpoly_vals)
        is_constant = (len(unique_vals) == 1)

        if is_constant:
            verdict = f"CONSTANT={superpoly_vals[0]} -> degree(f) <= {k}"
        else:
            # Test linearity of superpoly
            # The superpoly is a function of remaining variables (w12, w13, w14_rest, w15)
            # Test: p(x ^ a) ^ p(x) ^ p(a) ^ p(0) = 0 for random x, a
            # Here "x" and "a" are the full remaining variable set
            lin_pass = 0
            lin_total = 100
            for _ in range(lin_total):
                # Random base points
                x_w12 = random.randint(0, MASK)
                x_w13 = random.randint(0, MASK)
                x_w14 = random.randint(0, MASK)
                x_w15 = random.randint(0, MASK)

                a_w12 = random.randint(0, MASK)
                a_w13 = random.randint(0, MASK)
                a_w14 = random.randint(0, MASK)
                a_w15 = random.randint(0, MASK)

                p_x = cube_sum_w14(cube_pos, bit_j, x_w12, x_w13, x_w14, x_w15)
                p_a = cube_sum_w14(cube_pos, bit_j, a_w12, a_w13, a_w14, a_w15)
                p_xa = cube_sum_w14(cube_pos, bit_j,
                                    x_w12 ^ a_w12, x_w13 ^ a_w13,
                                    x_w14 ^ a_w14, x_w15 ^ a_w15)
                p_0 = cube_sum_w14(cube_pos, bit_j, 0, 0, 0, 0)

                if p_xa ^ p_x ^ p_a ^ p_0 == 0:
                    lin_pass += 1

            if lin_pass == lin_total:
                verdict = f"LINEAR superpoly -> degree(f) = {k+1}"
            else:
                verdict = f"NONLINEAR superpoly (lin_pass={lin_pass}/{lin_total}) -> degree(f) > {k+1}"

        print(f"    k={k} cube_bits={cube_pos}: "
              f"superpoly vals={n_ones}/{N_RANDOM} ones | {verdict}")
    print()

# ============================================================
# PART 2: Cross-word cubes (W[12] + W[14])
# ============================================================
print("-" * 72)
print("PART 2: Cross-word cubes (bits from W[12] AND W[14]) on De17")
print("-" * 72)
print()

def cube_sum_cross(cube_w12_pos, cube_w14_pos, bit_j, w12_rest, w13, w14_rest, w15):
    """
    Sum De17[bit_j] over all assignments to cube bits in both W[12] and W[14].
    """
    k12 = len(cube_w12_pos)
    k14 = len(cube_w14_pos)
    k = k12 + k14

    mask12 = 0
    for p in cube_w12_pos:
        mask12 |= (1 << p)
    mask14 = 0
    for p in cube_w14_pos:
        mask14 |= (1 << p)

    total = 0
    for assignment in range(1 << k):
        w12 = w12_rest & ~mask12
        w14 = w14_rest & ~mask14
        for idx, p in enumerate(cube_w12_pos):
            if assignment & (1 << idx):
                w12 |= (1 << p)
        for idx, p in enumerate(cube_w14_pos):
            if assignment & (1 << (k12 + idx)):
                w14 |= (1 << p)
        total ^= eval_de17_bit(w12, w13, w14, w15, bit_j)
    return total

for bit_j in test_bits:
    print(f"  Output bit j={bit_j}:")
    for k in [2, 3, 4, 5, 6]:
        # Split cube: half from W[12], half from W[14]
        k12 = k // 2
        k14 = k - k12
        pos_w12 = sorted(random.sample(range(32), k12))
        pos_w14 = sorted(random.sample(range(32), k14))

        superpoly_vals = []
        for _ in range(N_RANDOM):
            rw12 = random.randint(0, MASK)
            rw13 = random.randint(0, MASK)
            rw14 = random.randint(0, MASK)
            rw15 = random.randint(0, MASK)
            sp_val = cube_sum_cross(pos_w12, pos_w14, bit_j,
                                    rw12, rw13, rw14, rw15)
            superpoly_vals.append(sp_val)

        unique_vals = set(superpoly_vals)
        n_ones = sum(superpoly_vals)
        is_constant = (len(unique_vals) == 1)

        if is_constant:
            verdict = f"CONSTANT={superpoly_vals[0]} -> degree(f) <= {k}"
        else:
            lin_pass = 0
            lin_total = 100
            for _ in range(lin_total):
                x_w12 = random.randint(0, MASK)
                x_w13 = random.randint(0, MASK)
                x_w14 = random.randint(0, MASK)
                x_w15 = random.randint(0, MASK)
                a_w12 = random.randint(0, MASK)
                a_w13 = random.randint(0, MASK)
                a_w14 = random.randint(0, MASK)
                a_w15 = random.randint(0, MASK)

                p_x = cube_sum_cross(pos_w12, pos_w14, bit_j,
                                     x_w12, x_w13, x_w14, x_w15)
                p_a = cube_sum_cross(pos_w12, pos_w14, bit_j,
                                     a_w12, a_w13, a_w14, a_w15)
                p_xa = cube_sum_cross(pos_w12, pos_w14, bit_j,
                                      x_w12 ^ a_w12, x_w13 ^ a_w13,
                                      x_w14 ^ a_w14, x_w15 ^ a_w15)
                p_0 = cube_sum_cross(pos_w12, pos_w14, bit_j,
                                     0, 0, 0, 0)

                if p_xa ^ p_x ^ p_a ^ p_0 == 0:
                    lin_pass += 1

            if lin_pass == lin_total:
                verdict = f"LINEAR superpoly -> degree(f) = {k+1}"
            else:
                verdict = f"NONLINEAR superpoly (lin_pass={lin_pass}/{lin_total}) -> degree(f) > {k+1}"

        print(f"    k={k} W12_bits={pos_w12} W14_bits={pos_w14}: "
              f"{n_ones}/{N_RANDOM} ones | {verdict}")
    print()

# ============================================================
# PART 3: De13 cube control test (should always be constant)
# ============================================================
print("-" * 72)
print("PART 3: Control — cube attack on De13 (should NOT depend on W[12..15])")
print("-" * 72)
print()

def eval_de13_bit(w12, w13, w14, w15, bit_j):
    W16 = make_W16(w12, w13, w14, w15)
    de13 = compute_De(W16, 13)
    return (de13 >> bit_j) & 1

de13_ctrl_pass = 0
de13_ctrl_total = 0
for bit_j in test_bits:
    for k in [1, 3, 5]:
        cube_pos = sorted(random.sample(range(32), k))
        superpoly_vals = []
        for _ in range(50):
            rw12 = random.randint(0, MASK)
            rw13 = random.randint(0, MASK)
            rw14 = random.randint(0, MASK)
            rw15 = random.randint(0, MASK)

            # Cube sum over W[14] bits
            cube_mask = 0
            for p in cube_pos:
                cube_mask |= (1 << p)
            total = 0
            for assignment in range(1 << k):
                w14 = rw14 & ~cube_mask
                for idx, p in enumerate(cube_pos):
                    if assignment & (1 << idx):
                        w14 |= (1 << p)
                total ^= eval_de13_bit(rw12, rw13, w14, rw15, bit_j)
            superpoly_vals.append(total)

        is_constant = (len(set(superpoly_vals)) == 1)
        de13_ctrl_total += 1
        if is_constant:
            de13_ctrl_pass += 1
        status = "CONSTANT (expected)" if is_constant else "VARIES (unexpected!)"
        print(f"  De13 bit={bit_j} k={k}: superpoly is {status}")

print(f"\n  De13 control: {de13_ctrl_pass}/{de13_ctrl_total} constant")
print()

# ============================================================
# PART 4: Summary and degree estimates
# ============================================================
print("=" * 72)
print("PART 4: Summary — Degree Estimates and Verdict")
print("=" * 72)
print()

# Re-run systematic degree estimation for each output bit
print("Systematic degree lower-bound for De17[j] as function of W[14]:")
print("(Using single-word cubes of increasing dimension)")
print()

for bit_j in test_bits:
    degree_lower = 0
    for k in range(1, 7):
        cube_pos = sorted(random.sample(range(32), k))

        # Check if superpoly is constant across many random points
        vals = set()
        for trial in range(N_RANDOM):
            rw12 = random.randint(0, MASK)
            rw13 = random.randint(0, MASK)
            rw14 = random.randint(0, MASK)
            rw15 = random.randint(0, MASK)
            sp = cube_sum_w14(cube_pos, bit_j, rw12, rw13, rw14, rw15)
            vals.add(sp)
            if len(vals) > 1:
                break  # Not constant, no need to continue

        if len(vals) > 1:
            degree_lower = k + 1
        else:
            # Constant superpoly: degree <= k
            print(f"  De17[{bit_j:2d}]: degree <= {k} "
                  f"(constant superpoly at k={k})")
            break
    else:
        print(f"  De17[{bit_j:2d}]: degree >= {degree_lower} "
              f"(non-constant superpoly even at k=6)")

print()

# Cross-word comparison
print("Cross-word cube degree comparison:")
for bit_j in [0, 16]:
    # Single-word k=4
    cube_w14 = sorted(random.sample(range(32), 4))
    sw_vals = set()
    for _ in range(N_RANDOM):
        rw12 = random.randint(0, MASK)
        rw13 = random.randint(0, MASK)
        rw14 = random.randint(0, MASK)
        rw15 = random.randint(0, MASK)
        sp = cube_sum_w14(cube_w14, bit_j, rw12, rw13, rw14, rw15)
        sw_vals.add(sp)
        if len(sw_vals) > 1:
            break
    sw_const = (len(sw_vals) == 1)

    # Cross-word k=4 (2 from W[12], 2 from W[14])
    cw12 = sorted(random.sample(range(32), 2))
    cw14 = sorted(random.sample(range(32), 2))
    cw_vals = set()
    for _ in range(N_RANDOM):
        rw12 = random.randint(0, MASK)
        rw13 = random.randint(0, MASK)
        rw14 = random.randint(0, MASK)
        rw15 = random.randint(0, MASK)
        sp = cube_sum_cross(cw12, cw14, bit_j, rw12, rw13, rw14, rw15)
        cw_vals.add(sp)
        if len(cw_vals) > 1:
            break
    cw_const = (len(cw_vals) == 1)

    print(f"  De17[{bit_j}] k=4: single-word={'CONST' if sw_const else 'VARIES'}, "
          f"cross-word={'CONST' if cw_const else 'VARIES'}")
    if cw_const and not sw_const:
        print(f"    -> ANOMALY: cross-word cubes find lower degree!")
    elif sw_const and not cw_const:
        print(f"    -> Single-word cubes more effective for this bit")
    elif sw_const and cw_const:
        print(f"    -> Both constant: degree <= 4")
    else:
        print(f"    -> Both vary: degree > 5")

print()

# Final verdict
print("-" * 72)
print("FINAL VERDICT")
print("-" * 72)
print()
print("De17 as a function of free words W[12..15]:")
print("  - If superpoly degree <= 3 at cube dim 4-6: barrier can be linearized")
print("    -> potential to solve De17=0 faster than 2^32 brute force")
print("  - If superpoly remains nonlinear even at dim 6: degree >= 7")
print("    -> cube attack at these dimensions does NOT break the barrier")
print("  - Cross-word cubes finding lower degree would indicate")
print("    inter-word algebraic coupling (exploitable structure)")
print()

# Count overall results
n_const = 0
n_linear = 0
n_nonlinear = 0
for bit_j in test_bits:
    for k in range(1, 7):
        cube_pos = sorted(random.sample(range(32), k))
        vals = set()
        for _ in range(60):
            rw12 = random.randint(0, MASK)
            rw13 = random.randint(0, MASK)
            rw14 = random.randint(0, MASK)
            rw15 = random.randint(0, MASK)
            sp = cube_sum_w14(cube_pos, bit_j, rw12, rw13, rw14, rw15)
            vals.add(sp)
            if len(vals) > 1:
                break
        if len(vals) == 1:
            n_const += 1
        else:
            # Quick linearity check
            lin_ok = True
            for _ in range(30):
                x12 = random.randint(0, MASK); x13 = random.randint(0, MASK)
                x14 = random.randint(0, MASK); x15 = random.randint(0, MASK)
                a12 = random.randint(0, MASK); a13 = random.randint(0, MASK)
                a14 = random.randint(0, MASK); a15 = random.randint(0, MASK)
                p_x = cube_sum_w14(cube_pos, bit_j, x12, x13, x14, x15)
                p_a = cube_sum_w14(cube_pos, bit_j, a12, a13, a14, a15)
                p_xa = cube_sum_w14(cube_pos, bit_j,
                                    x12^a12, x13^a13, x14^a14, x15^a15)
                p_0 = cube_sum_w14(cube_pos, bit_j, 0, 0, 0, 0)
                if p_xa ^ p_x ^ p_a ^ p_0 != 0:
                    lin_ok = False
                    break
            if lin_ok:
                n_linear += 1
            else:
                n_nonlinear += 1

total_tests = n_const + n_linear + n_nonlinear
print(f"Across {total_tests} cube tests (4 bits x 6 dimensions):")
print(f"  Constant superpolys:  {n_const:3d}  ({100*n_const/total_tests:.1f}%)")
print(f"  Linear superpolys:    {n_linear:3d}  ({100*n_linear/total_tests:.1f}%)")
print(f"  Nonlinear superpolys: {n_nonlinear:3d}  ({100*n_nonlinear/total_tests:.1f}%)")
print()

if n_const > total_tests * 0.3:
    print("VERDICT: ALIVE — significant fraction of constant superpolys")
    print("  -> cube attack can reduce De17 barrier complexity")
elif n_linear > total_tests * 0.3:
    print("VERDICT: ALIVE — significant fraction of linear superpolys")
    print("  -> linearization attack may reduce barrier below 2^32")
elif n_const + n_linear > 0:
    print("VERDICT: ANOMALY — some low-degree superpolys found but minority")
    print("  -> partial algebraic structure, but not enough for full attack")
else:
    print("VERDICT: DEAD — all superpolys are nonlinear at tested dimensions")
    print("  -> De17 behaves as high-degree function of free words")
    print("  -> cube attack at dim <= 6 does NOT break the 2^32 barrier")
