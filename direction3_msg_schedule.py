#!/usr/bin/env python3
"""
Direction 3: Deep analysis of the SHA-256 message schedule.

W[t] = sigma1(W[t-2]) + W[t-7] + sigma0(W[t-15]) + W[t-16]  for t >= 16

Studies:
  1. Carry density in message schedule additions
  2. Low-weight message schedule paths (hill climbing)
  3. Differential propagation: actual vs GF(2)-linear
  4. GF(2) linear recurrence rank analysis
  5. Carry predictability for 1-bit input differences
"""

import random
import time

MASK32 = 0xFFFFFFFF

def rotr(x, n):
    return ((x >> n) | (x << (32 - n))) & MASK32

def shr(x, n):
    return x >> n

def sigma0(x):
    return rotr(x, 7) ^ rotr(x, 18) ^ shr(x, 3)

def sigma1(x):
    return rotr(x, 17) ^ rotr(x, 19) ^ shr(x, 10)

def hw(x):
    """Hamming weight."""
    return bin(x).count('1')

def expand_schedule(W):
    """Expand W[0..15] to W[0..63] using SHA-256 message schedule."""
    W = list(W)
    for t in range(16, 64):
        W.append((sigma1(W[t-2]) + W[t-7] + sigma0(W[t-15]) + W[t-16]) & MASK32)
    return W

def expand_schedule_xor(W):
    """Expand using XOR instead of addition (no carries) - GF(2) version."""
    W = list(W)
    for t in range(16, 64):
        val = sigma1(W[t-2]) ^ W[t-7] ^ sigma0(W[t-15]) ^ W[t-16]
        W.append(val & MASK32)
    return W

def random_w():
    """Generate random 16-word message block."""
    return [random.getrandbits(32) for _ in range(16)]


def study1_carry_density(n_samples=5000):
    """Study 1: Carry density per addition in message schedule.
    Measures carry bits LOCAL to each addition: W[t] = A + B + C + D mod 2^32.
    Compare (A+B+C+D) mod 2^32 vs (A XOR B XOR C XOR D) at each step."""
    print("=" * 70)
    print("STUDY 1: Carry Density in Message Schedule Additions")
    print("=" * 70)
    print(f"Samples: {n_samples}")
    print()

    carry_bits = [0] * 48  # for t=16..63

    for _ in range(n_samples):
        W = list(random_w())
        for t in range(16, 64):
            a = sigma1(W[t-2])
            b = W[t-7]
            c = sigma0(W[t-15])
            d = W[t-16]
            actual = (a + b + c + d) & MASK32
            xor_val = a ^ b ^ c ^ d
            carry_diff = actual ^ xor_val
            carry_bits[t - 16] += hw(carry_diff)
            W.append(actual)

    print(f"{'Round t':>8}  {'Avg carry bits':>14}  {'Carry density':>13}  {'Bar'}")
    print("-" * 60)
    for idx in range(48):
        t = idx + 16
        avg = carry_bits[idx] / n_samples
        density = avg / 32.0
        bar = "#" * min(60, int(density * 40))
        if idx < 8 or idx >= 44 or idx % 4 == 0:
            print(f"{t:>8}  {avg:>14.3f}  {density:>13.4f}  {bar}")

    overall_avg = sum(carry_bits) / (48 * n_samples)
    print(f"\nOverall average carry bits per word: {overall_avg:.3f} ({overall_avg/32:.4f} density)")
    print()


def study2_low_weight_paths(n_iter=5000):
    """Study 2: Hill climbing for low-weight message schedule paths."""
    print("=" * 70)
    print("STUDY 2: Low-Weight Message Schedule Paths (Hill Climbing)")
    print("=" * 70)
    print(f"Iterations: {n_iter}")
    print()

    # Baseline: average total HW for random inputs
    n_baseline = 1000
    baseline_hws = []
    for _ in range(n_baseline):
        W = expand_schedule(random_w())
        total = sum(hw(W[t]) for t in range(16, 64))
        baseline_hws.append(total)
    baseline_avg = sum(baseline_hws) / len(baseline_hws)
    baseline_min = min(baseline_hws)

    # Hill climbing
    best_W = random_w()
    best_full = expand_schedule(best_W)
    best_hw = sum(hw(best_full[t]) for t in range(16, 64))

    # Try multiple restarts
    n_restarts = 10
    iters_per = n_iter // n_restarts
    global_best_hw = best_hw
    global_best_W = list(best_W)

    for restart in range(n_restarts):
        curr_W = random_w()
        curr_full = expand_schedule(curr_W)
        curr_hw = sum(hw(curr_full[t]) for t in range(16, 64))

        for _ in range(iters_per):
            # Pick random word and bit to flip
            word_idx = random.randint(0, 15)
            bit_idx = random.randint(0, 31)

            new_W = list(curr_W)
            new_W[word_idx] ^= (1 << bit_idx)
            new_full = expand_schedule(new_W)
            new_hw = sum(hw(new_full[t]) for t in range(16, 64))

            if new_hw <= curr_hw:
                curr_W = new_W
                curr_full = new_full
                curr_hw = new_hw

        if curr_hw < global_best_hw:
            global_best_hw = curr_hw
            global_best_W = list(curr_W)

    print(f"Random baseline avg total HW(W[16..63]): {baseline_avg:.1f}")
    print(f"Random baseline min total HW(W[16..63]): {baseline_min}")
    print(f"Best found via hill climbing:             {global_best_hw}")
    print(f"Reduction from avg baseline:              {baseline_avg - global_best_hw:.1f} "
          f"({(baseline_avg - global_best_hw)/baseline_avg*100:.1f}%)")
    print(f"Theoretical minimum (all zeros input):    ", end="")
    zero_W = expand_schedule([0]*16)
    zero_hw = sum(hw(zero_W[t]) for t in range(16, 64))
    print(f"{zero_hw}")
    print()

    # Show HW distribution of best solution
    best_full = expand_schedule(global_best_W)
    print("Best solution W[16..63] HW per word:")
    for t in range(16, 64):
        if t < 24 or t >= 60 or t % 8 == 0:
            print(f"  W[{t:2d}]: HW={hw(best_full[t]):2d}  {'*' * hw(best_full[t])}")
    print()


def study3_differential_propagation(n_samples=1000):
    """Study 3: Differential propagation - actual vs GF(2)-linear."""
    print("=" * 70)
    print("STUDY 3: Differential Propagation in Message Schedule")
    print("=" * 70)
    print(f"Samples per input word: {n_samples}")
    print()

    print(f"{'Input':>8}  {'Actual avg HW':>14}  {'GF2 avg HW':>11}  {'Ratio':>7}  {'Max dev':>8}")
    print("-" * 60)

    for i in range(16):
        actual_hw_total = 0
        gf2_hw_total = 0
        max_deviation = 0

        for _ in range(n_samples):
            W0 = random_w()
            bit = random.randint(0, 31)

            W0_prime = list(W0)
            W0_prime[i] ^= (1 << bit)

            W_actual = expand_schedule(W0)
            W_prime = expand_schedule(W0_prime)

            W_xor_base = expand_schedule_xor(W0)
            W_xor_prime = expand_schedule_xor(W0_prime)

            actual_diff_hw = 0
            gf2_diff_hw = 0
            for t in range(16, 64):
                ad = W_actual[t] ^ W_prime[t]
                gd = W_xor_base[t] ^ W_xor_prime[t]
                actual_diff_hw += hw(ad)
                gf2_diff_hw += hw(gd)

            actual_hw_total += actual_diff_hw
            gf2_hw_total += gf2_diff_hw
            dev = abs(actual_diff_hw - gf2_diff_hw)
            if dev > max_deviation:
                max_deviation = dev

        avg_actual = actual_hw_total / n_samples
        avg_gf2 = gf2_hw_total / n_samples
        ratio = avg_actual / avg_gf2 if avg_gf2 > 0 else float('inf')
        print(f"  W[{i:2d}]   {avg_actual:>14.2f}  {avg_gf2:>11.2f}  {ratio:>7.3f}  {max_deviation:>8d}")

    print()
    print("Ratio > 1 means carries amplify diffusion beyond linear prediction.")
    print("Ratio < 1 means carries sometimes cancel diffs (rare).")
    print()


def study4_gf2_rank():
    """Study 4: GF(2) linear recurrence - rank of dependency matrix."""
    print("=" * 70)
    print("STUDY 4: GF(2) Linear Recurrence Rank Analysis")
    print("=" * 70)
    print()

    # Build the 48*32 x 16*32 = 1536 x 512 binary dependency matrix
    # over GF(2). Each row is one output bit of W[16..63],
    # each column is one input bit of W[0..15].
    # We compute this by setting each input bit one at a time.

    n_out = 48 * 32  # 1536 output bits
    n_in = 16 * 32   # 512 input bits

    # Build matrix as list of n_out integers, each is a bitmask of n_in bits
    # Row r corresponds to bit (r%32) of W[16 + r//32]
    # We compute column by column: set one input bit, see which output bits flip.

    # More efficient: compute row by row as bitmask of which input bits affect it
    # Actually, column approach is simpler: for each input bit, run GF(2) schedule,
    # record which output bits are set.

    matrix = [0] * n_out  # each is an integer with up to 512 bits

    for col in range(n_in):
        word_idx = col // 32
        bit_idx = col % 32

        W = [0] * 16
        W[word_idx] = 1 << bit_idx

        W = expand_schedule_xor(W)

        for t in range(16, 64):
            val = W[t]
            for b in range(32):
                if val & (1 << b):
                    row = (t - 16) * 32 + b
                    matrix[row] |= (1 << col)

    # Gaussian elimination over GF(2) to find rank
    rank = 0
    mat = list(matrix)  # copy
    pivot_col = 0
    n_rows = len(mat)

    # Use standard GF(2) Gaussian elimination
    used_rows = [False] * n_rows
    for col in range(n_in):
        # Find a row with this column set
        pivot_row = -1
        for r in range(n_rows):
            if not used_rows[r] and (mat[r] & (1 << col)):
                pivot_row = r
                break
        if pivot_row == -1:
            continue
        rank += 1
        used_rows[pivot_row] = True
        pivot_val = mat[pivot_row]
        for r in range(n_rows):
            if r != pivot_row and (mat[r] & (1 << col)):
                mat[r] ^= pivot_val

    print(f"Dependency matrix dimensions: {n_out} rows x {n_in} cols")
    print(f"  (Each row = one output bit of W[16..63])")
    print(f"  (Each col = one input bit of W[0..15])")
    print(f"GF(2) rank: {rank} / {n_in}")
    print()

    if rank == n_in:
        print("RESULT: Full rank! All 512 input bits are linearly independent")
        print("in their GF(2) contribution to W[16..63].")
        print("This means the GF(2)-linear part of the schedule is a bijection")
        print("from input bits to a 512-dimensional subspace of the 1536-bit output.")
    else:
        print(f"RESULT: Rank deficiency detected. {n_in - rank} input bits are")
        print(f"linearly dependent over GF(2) in terms of schedule output.")

    # Per-round rank analysis
    print("\nPer-round cumulative GF(2) rank (W[16] through W[16+k]):")
    print(f"{'Through W[t]':>14}  {'Output bits':>12}  {'Rank':>6}  {'Max possible':>13}")
    for k_end in [0, 3, 7, 15, 23, 31, 39, 47]:
        sub_rows = (k_end + 1) * 32
        sub_mat = list(matrix[:sub_rows])
        r = 0
        used = [False] * sub_rows
        for c in range(n_in):
            pr = -1
            for rr in range(sub_rows):
                if not used[rr] and (sub_mat[rr] & (1 << c)):
                    pr = rr
                    break
            if pr == -1:
                continue
            r += 1
            used[pr] = True
            pv = sub_mat[pr]
            for rr in range(sub_rows):
                if rr != pr and (sub_mat[rr] & (1 << c)):
                    sub_mat[rr] ^= pv
        max_rank = min(n_in, sub_rows)
        t_val = 16 + k_end
        print(f"       W[{t_val:2d}]    {sub_rows:>12}  {r:>6}  {max_rank:>13}")

    print()


def study5_carry_predictability(n_samples=5000):
    """Study 5: Carry predictability for 1-bit diffs in W[0]."""
    print("=" * 70)
    print("STUDY 5: Carry Predictability for 1-bit Input Differences")
    print("=" * 70)
    print(f"Samples per bit position: {n_samples}")
    print()

    # Focus on flipping each bit of W[0], measuring P(carry at each bit of W[16])
    # W[16] = sigma1(W[14]) + W[9] + sigma0(W[1]) + W[0]
    # Only W[0] term is affected directly.

    # For each input bit position in W[0], measure the probability that
    # each bit of W[16] has a carry contribution
    # (i.e., actual_diff != xor_diff at that bit position)

    print("Flipping bits in W[0], measuring carry in W[16]:")
    print(f"{'Input bit':>10}  {'Avg carry bits':>14}  {'Most predictable':>17}  {'Least predictable':>18}")
    print("-" * 65)

    # Track per-bit carry probability for all 32 input bits
    bit_carry_probs = {}  # input_bit -> list of 32 probabilities

    for in_bit in range(32):
        carry_counts = [0] * 32  # per output bit of W[16]

        for _ in range(n_samples):
            W0 = random_w()
            W0_prime = list(W0)
            W0_prime[0] ^= (1 << in_bit)

            W = expand_schedule(W0)
            W_prime = expand_schedule(W0_prime)

            W_xor = expand_schedule_xor(W0)
            W_xor_prime = expand_schedule_xor(W0_prime)

            actual_diff = W[16] ^ W_prime[16]
            gf2_diff = W_xor[16] ^ W_xor_prime[16]
            carry_diff = actual_diff ^ gf2_diff

            for b in range(32):
                if carry_diff & (1 << b):
                    carry_counts[b] += 1

        probs = [c / n_samples for c in carry_counts]
        bit_carry_probs[in_bit] = probs
        avg_carry = sum(probs)

        # Most predictable = closest to 0 or 1
        predictabilities = [max(p, 1-p) for p in probs]
        most_pred_bit = max(range(32), key=lambda b: predictabilities[b])
        most_pred_val = predictabilities[most_pred_bit]
        least_pred_bit = min(range(32), key=lambda b: predictabilities[b])
        least_pred_val = predictabilities[least_pred_bit]

        if in_bit % 4 == 0 or in_bit < 4 or in_bit >= 28:
            print(f"{in_bit:>10}  {avg_carry:>14.3f}  "
                  f"bit {most_pred_bit:2d} (p={most_pred_val:.3f})  "
                  f"bit {least_pred_bit:2d} (p={least_pred_val:.3f})")

    # Overall summary
    print()
    print("Carry probability heatmap for W[16] (rows=input bit, cols=output bit):")
    print("  (showing probability that carry affects output bit, 0=never, 9=always)")
    print()
    header = "     " + "".join(f"{b%10}" for b in range(32))
    print(header)
    print("     " + "-" * 32)
    for in_bit in range(32):
        row = f"{in_bit:2d} | "
        for out_bit in range(32):
            p = bit_carry_probs[in_bit][out_bit]
            # Map to 0-9 scale
            c = int(p * 9.99)
            if c > 9:
                c = 9
            row += str(c)
        print(row)

    print()
    print("Legend: 0 = carry never occurs, 5 = carry occurs ~50%, 9 = carry always occurs")

    # Find the overall most and least predictable combinations
    print()
    all_pairs = []
    for ib in range(32):
        for ob in range(32):
            p = bit_carry_probs[ib][ob]
            entropy = -(p * (max(1e-10, p)) + (1-p) * max(1e-10, 1-p))  # not real entropy
            pred = max(p, 1-p)
            all_pairs.append((pred, ib, ob, p))
    all_pairs.sort()

    print("5 least predictable carry positions (closest to 50/50):")
    for pred, ib, ob, p in all_pairs[:5]:
        print(f"  Input bit {ib} -> W[16] bit {ob}: P(carry)={p:.4f} (predictability={pred:.4f})")

    print("\n5 most predictable carry positions:")
    for pred, ib, ob, p in all_pairs[-5:]:
        print(f"  Input bit {ib} -> W[16] bit {ob}: P(carry)={p:.4f} (predictability={pred:.4f})")

    print()


def main():
    t0 = time.time()
    print("SHA-256 Message Schedule Deep Analysis")
    print("=" * 70)
    print()

    random.seed(0x53484132353600)

    t1 = time.time()
    study1_carry_density(5000)
    print(f"[Study 1 completed in {time.time()-t1:.1f}s]\n")

    t2 = time.time()
    study2_low_weight_paths(5000)
    print(f"[Study 2 completed in {time.time()-t2:.1f}s]\n")

    t3 = time.time()
    study3_differential_propagation(1000)
    print(f"[Study 3 completed in {time.time()-t3:.1f}s]\n")

    t4 = time.time()
    study4_gf2_rank()
    print(f"[Study 4 completed in {time.time()-t4:.1f}s]\n")

    t5 = time.time()
    study5_carry_predictability(2000)
    print(f"[Study 5 completed in {time.time()-t5:.1f}s]\n")

    total = time.time() - t0
    print("=" * 70)
    print(f"TOTAL TIME: {total:.1f}s")
    if total < 45:
        print("PASS: Completed within 45-second limit.")
    else:
        print("FAIL: Exceeded 45-second limit!")

if __name__ == "__main__":
    main()
