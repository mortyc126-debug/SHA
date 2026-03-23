"""
Verification of weapon_combined_v2 differential distinguisher.
Key question: is the 12-round result real or artifact?

Test plan:
1. Run differential KS-test against ACTUAL random oracle => should give p ~ uniform
2. Run against SHA-256 rounds 8, 10, 12 => check if p really underflows
3. Multiple seeds to rule out seed-dependent artifacts
4. Check SHA-256 implementation correctness via known test vector
"""

import numpy as np
from scipy import stats
import hashlib, struct, time

# === SHA-256 implementation (copied from weapon_combined_v2.py) ===
K = [0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
     0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174]
IV = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]

def rotr(x, n):
    return ((x >> n) | (x << (32 - n))) & 0xFFFFFFFF

def sha256_rounds(W, nr):
    w = list(W[:16])
    for i in range(16, max(nr, 16)):
        s0 = rotr(w[i-15], 7) ^ rotr(w[i-15], 18) ^ (w[i-15] >> 3)
        s1 = rotr(w[i-2], 17) ^ rotr(w[i-2], 19) ^ (w[i-2] >> 10)
        w.append((w[i-16] + s0 + w[i-7] + s1) & 0xFFFFFFFF)
    a, b, c, d, e, f, g, h = IV
    for i in range(nr):
        S1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25)
        ch = (e & f) ^ ((~e) & g) & 0xFFFFFFFF
        t1 = (h + S1 + ch + K[i] + w[i]) & 0xFFFFFFFF
        S0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22)
        maj = (a & b) ^ (a & c) ^ (b & c)
        t2 = (S0 + maj) & 0xFFFFFFFF
        h, g, f, e, d, c, b, a = g, f, e, (d + t1) & 0xFFFFFFFF, c, b, a, (t1 + t2) & 0xFFFFFFFF
    return [a, b, c, d, e, f, g, h]

def hw(x):
    return bin(x).count('1')

# === Test 0: Verify SHA-256 implementation against hashlib ===
def verify_implementation():
    """Check our reduced-round impl matches hashlib for structure."""
    print("=" * 60)
    print("TEST 0: SHA-256 implementation sanity check")
    print("=" * 60)

    # We can't directly compare reduced rounds with hashlib (which does 64),
    # but we can verify the round function is correct by checking that
    # 64 rounds + feed-forward matches hashlib.

    # Full 64-round implementation for verification
    K64 = [
        0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
        0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
        0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
        0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
        0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
        0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
        0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
        0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
    ]

    msg = b"abc"
    # Pad manually
    ml = len(msg) * 8
    padded = msg + b'\x80' + b'\x00' * (55 - len(msg)) + struct.pack('>Q', ml)
    W = list(struct.unpack('>16I', padded))

    # Extend schedule
    w = list(W)
    for i in range(16, 64):
        s0 = rotr(w[i-15], 7) ^ rotr(w[i-15], 18) ^ (w[i-15] >> 3)
        s1 = rotr(w[i-2], 17) ^ rotr(w[i-2], 19) ^ (w[i-2] >> 10)
        w.append((w[i-16] + s0 + w[i-7] + s1) & 0xFFFFFFFF)

    a, b, c, d, e, f, g, h = IV
    for i in range(64):
        S1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25)
        ch = (e & f) ^ ((~e) & g) & 0xFFFFFFFF
        t1 = (h + S1 + ch + K64[i] + w[i]) & 0xFFFFFFFF
        S0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22)
        maj = (a & b) ^ (a & c) ^ (b & c)
        t2 = (S0 + maj) & 0xFFFFFFFF
        h, g, f, e, d, c, b, a = g, f, e, (d + t1) & 0xFFFFFFFF, c, b, a, (t1 + t2) & 0xFFFFFFFF

    # Add feed-forward
    result = [(x + y) & 0xFFFFFFFF for x, y in zip([a,b,c,d,e,f,g,h], IV)]
    our_hash = ''.join(f'{x:08x}' for x in result)

    expected = hashlib.sha256(msg).hexdigest()

    match = our_hash == expected
    print(f"  Our hash:    {our_hash}")
    print(f"  hashlib:     {expected}")
    print(f"  Match: {'YES' if match else 'NO *** BUG ***'}")

    if not match:
        # Debug: check Ch precedence
        e_val, f_val, g_val = 0xAAAAAAAA, 0x55555555, 0xFFFFFFFF
        ch_weapon = (e_val & f_val) ^ ((~e_val) & g_val) & 0xFFFFFFFF
        ch_correct = ((e_val & f_val) ^ ((~e_val) & g_val)) & 0xFFFFFFFF
        print(f"  Ch precedence check: weapon={ch_weapon:#010x}, correct={ch_correct:#010x}")
        if ch_weapon != ch_correct:
            print("  *** Ch has OPERATOR PRECEDENCE BUG! ***")

    return match

# === Test 1: Null hypothesis - random oracle ===
def test_null_hypothesis(n_trials=20, n_pairs=10000):
    """Run differential KS-test against random oracle.
    Should give p-values uniformly distributed on [0,1]."""
    print("\n" + "=" * 60)
    print("TEST 1: Null hypothesis (random oracle)")
    print(f"  {n_trials} trials, {n_pairs} pairs each")
    print("=" * 60)

    pvals = []
    for trial in range(n_trials):
        hw_vals = np.empty(n_pairs, dtype=np.int32)
        for i in range(n_pairs):
            # Two independent random 32-bit values (simulating random oracle)
            r1 = np.random.randint(0, 0xFFFFFFFF + 1, dtype=np.uint64)
            r2 = np.random.randint(0, 0xFFFFFFFF + 1, dtype=np.uint64)
            diff = int(r1) ^ int(r2)
            hw_vals[i] = hw(diff)

        expected = stats.binom(32, 0.5)
        _, pval = stats.kstest(hw_vals, expected.cdf)
        pvals.append(pval)

    pvals = np.array(pvals)
    print(f"  p-values: min={pvals.min():.6e}, median={np.median(pvals):.4f}, max={pvals.max():.4f}")
    print(f"  p < 0.05: {np.sum(pvals < 0.05)}/{n_trials} (expected ~{n_trials*0.05:.0f})")
    print(f"  p < 0.01: {np.sum(pvals < 0.01)}/{n_trials} (expected ~{n_trials*0.01:.0f})")

    # KS test on p-values themselves (should be uniform)
    _, meta_p = stats.kstest(pvals, 'uniform')
    print(f"  Meta-KS (p-values uniform?): p = {meta_p:.4f}")

    return pvals

# === Test 2: SHA-256 differential test, multiple seeds ===
def test_sha256_differential(rounds_list=[8, 10, 12], n_pairs=10000, n_seeds=5):
    """Run differential KS-test on SHA-256 with multiple seeds."""
    print("\n" + "=" * 60)
    print("TEST 2: SHA-256 differential KS-test (multiple seeds)")
    print(f"  Rounds: {rounds_list}, {n_pairs} pairs, {n_seeds} seeds")
    print("=" * 60)

    delta = 0x80000000

    for nr in rounds_list:
        print(f"\n  --- Round {nr} ---")
        pvals = []
        for seed in range(n_seeds):
            np.random.seed(seed * 1000 + 42)
            hw_vals = np.empty(n_pairs, dtype=np.int32)
            for i in range(n_pairs):
                W = [int(x) for x in np.random.randint(0, 0xFFFFFFFF + 1, 16, dtype=np.uint64)]
                W2 = list(W)
                W2[0] = W[0] ^ delta
                out1 = sha256_rounds(W, nr)
                out2 = sha256_rounds(W2, nr)
                diff_e = out1[4] ^ out2[4]
                hw_vals[i] = hw(diff_e)

            expected = stats.binom(32, 0.5)
            ks_stat, pval = stats.kstest(hw_vals, expected.cdf)
            pvals.append(pval)

            mean_hw = np.mean(hw_vals)
            std_hw = np.std(hw_vals)
            print(f"    seed={seed}: p={pval:.6e}, mean_HW={mean_hw:.3f}, std_HW={std_hw:.3f}")

        pvals = np.array(pvals)
        all_sig = np.all(pvals < 0.01)
        print(f"    All seeds p < 0.01? {'YES => REAL SIGNAL' if all_sig else 'NO => may be artifact'}")

# === Test 3: Critical check - is HW distribution actually non-binomial? ===
def test_hw_distribution_detail(nr=12, n_pairs=50000):
    """Detailed look at HW distribution for round 12."""
    print("\n" + "=" * 60)
    print(f"TEST 3: HW distribution detail (round {nr}, {n_pairs} pairs)")
    print("=" * 60)

    delta = 0x80000000
    np.random.seed(42)

    hw_vals = np.empty(n_pairs, dtype=np.int32)
    for i in range(n_pairs):
        W = [int(x) for x in np.random.randint(0, 0xFFFFFFFF + 1, 16, dtype=np.uint64)]
        W2 = list(W)
        W2[0] = W[0] ^ delta
        out1 = sha256_rounds(W, nr)
        out2 = sha256_rounds(W2, nr)
        diff_e = out1[4] ^ out2[4]
        hw_vals[i] = hw(diff_e)

    expected = stats.binom(32, 0.5)

    # Show distribution vs expected
    print(f"\n  {'HW':>4} | {'Observed':>10} | {'Expected':>10} | {'Ratio':>8}")
    print("  " + "-" * 45)
    for k in range(10, 23):
        obs = np.sum(hw_vals == k)
        exp = expected.pmf(k) * n_pairs
        ratio = obs / exp if exp > 0 else 0
        flag = " ***" if abs(ratio - 1.0) > 0.05 else ""
        print(f"  {k:4d} | {obs:10d} | {exp:10.1f} | {ratio:8.4f}{flag}")

    mean_hw = np.mean(hw_vals)
    std_hw = np.std(hw_vals)
    print(f"\n  mean = {mean_hw:.4f} (expected 16.0000)")
    print(f"  std  = {std_hw:.4f} (expected {np.sqrt(32*0.25):.4f})")

    ks_stat, pval = stats.kstest(hw_vals, expected.cdf)
    print(f"  KS stat = {ks_stat:.6f}, p = {pval:.6e}")

# === Test 4: Check ALL output registers, not just e ===
def test_all_registers(nr=12, n_pairs=10000):
    """Check if the bias is specific to register e or affects all."""
    print("\n" + "=" * 60)
    print(f"TEST 4: All registers differential KS-test (round {nr})")
    print("=" * 60)

    delta = 0x80000000
    np.random.seed(42)

    reg_names = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
    hw_all = [np.empty(n_pairs, dtype=np.int32) for _ in range(8)]

    for i in range(n_pairs):
        W = [int(x) for x in np.random.randint(0, 0xFFFFFFFF + 1, 16, dtype=np.uint64)]
        W2 = list(W)
        W2[0] = W[0] ^ delta
        out1 = sha256_rounds(W, nr)
        out2 = sha256_rounds(W2, nr)
        for r in range(8):
            hw_all[r][i] = hw(out1[r] ^ out2[r])

    expected = stats.binom(32, 0.5)
    for r in range(8):
        ks_stat, pval = stats.kstest(hw_all[r], expected.cdf)
        mean_hw = np.mean(hw_all[r])
        flag = "***" if pval < 0.001 else ""
        print(f"  {reg_names[r]}: p={pval:.6e}, mean_HW={mean_hw:.3f} {flag}")

# === Test 5: No-difference control (W1 = W2, should give HW=0) ===
def test_no_diff_control(nr=12, n_pairs=100):
    """Sanity: with delta=0, diff should be exactly 0."""
    print("\n" + "=" * 60)
    print(f"TEST 5: No-difference control (round {nr})")
    print("=" * 60)
    np.random.seed(42)

    all_zero = True
    for i in range(n_pairs):
        W = [int(x) for x in np.random.randint(0, 0xFFFFFFFF + 1, 16, dtype=np.uint64)]
        out1 = sha256_rounds(W, nr)
        out2 = sha256_rounds(list(W), nr)
        for r in range(8):
            if out1[r] != out2[r]:
                all_zero = False
                print(f"  BUG: pair {i}, register {r}: {out1[r]:#010x} != {out2[r]:#010x}")
                break
        if not all_zero:
            break

    print(f"  All diffs zero: {'YES (correct)' if all_zero else 'NO *** BUG ***'}")

if __name__ == "__main__":
    t0 = time.time()

    impl_ok = verify_implementation()

    if not impl_ok:
        print("\n*** SHA-256 IMPLEMENTATION IS BUGGY - all results suspect ***\n")

    test_null_hypothesis(n_trials=20, n_pairs=5000)
    test_no_diff_control(nr=12)
    test_sha256_differential(rounds_list=[8, 10, 12], n_pairs=5000, n_seeds=3)
    test_hw_distribution_detail(nr=12, n_pairs=20000)
    test_all_registers(nr=12, n_pairs=5000)

    print(f"\n{'='*60}")
    print(f"Total verification time: {time.time()-t0:.1f}s")
