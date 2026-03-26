#!/usr/bin/env python3
"""
SCF Carry Bitwise Decomposition — Per-bit carry analysis across SHA-256 rounds.

Goal: Decompose carry contributions BIT BY BIT to find hidden structure.

From ALG theory:
  - carry_correction cc = (a+b) XOR (a^b) for each addition
  - deg(cc[bit_i]) = i+1 (grows linearly with bit position)
  - Bit 0 is always carry-free (F3: carry-free LSB)
  - Carry survival P = 5/8 per bit position

Experiments:
  1. Per-bit carry profile (64 rounds x 32 bits x 4 carry sources)
  2. Bit-position heatmap (rounds x bits)
  3. Carry chain length per round
  4. Differential carry correlation
  5. Predictability window detection
"""

import os
import struct
import sys

MASK32 = 0xFFFFFFFF

# SHA-256 round constants
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
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
]

IV = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]

# --- Primitives ---
def rotr(x, n):
    return ((x >> n) | (x << (32 - n))) & MASK32

def shr(x, n):
    return x >> n

def Sig0(x):
    return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)

def Sig1(x):
    return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)

def sig0(x):
    return rotr(x, 7) ^ rotr(x, 18) ^ shr(x, 3)

def sig1(x):
    return rotr(x, 17) ^ rotr(x, 19) ^ shr(x, 10)

def Ch(e, f, g):
    return (e & f) ^ (~e & g) & MASK32

def Maj(a, b, c):
    return (a & b) ^ (a & c) ^ (b & c)

def add32(*args):
    s = 0
    for x in args:
        s = (s + x) & MASK32
    return s

def hw(x):
    return bin(x).count('1')

def bit(x, i):
    """Extract bit i from x."""
    return (x >> i) & 1

# --- Carry correction ---
def carry_correction_multi(*args):
    """cc = (a + b + c + ...) XOR (a ^ b ^ c ^ ...)"""
    arith = 0
    xor_sum = 0
    for x in args:
        arith = (arith + x) & MASK32
        xor_sum ^= x
    return arith ^ xor_sum

# --- Message schedule ---
def expand_schedule(W16):
    W = list(W16)
    for i in range(16, 64):
        W.append(add32(sig1(W[i-2]), W[i-7], sig0(W[i-15]), W[i-16]))
    return W

# --- Per-round carry decomposition ---
def round_carries(state, W_r, K_r):
    """
    Compute one SHA-256 round, returning per-carry-source corrections
    and the new state.
    """
    a, b, c, d, e, f, g, h = state

    # T1 = h + Sig1(e) + Ch(e,f,g) + K_r + W_r
    t1_operands = [h, Sig1(e), Ch(e, f, g), K_r, W_r]
    t1_real = add32(*t1_operands)
    t1_xor = 0
    for x in t1_operands:
        t1_xor ^= x
    cc_T1 = t1_real ^ t1_xor

    # T2 = Sig0(a) + Maj(a,b,c)
    t2_operands = [Sig0(a), Maj(a, b, c)]
    t2_real = add32(*t2_operands)
    t2_xor = Sig0(a) ^ Maj(a, b, c)
    cc_T2 = t2_real ^ t2_xor

    # e_new = d + T1
    e_real = add32(d, t1_real)
    cc_e = e_real ^ (d ^ t1_real)

    # a_new = T1 + T2
    a_real = add32(t1_real, t2_real)
    cc_a = a_real ^ (t1_real ^ t2_real)

    new_state = [a_real, a, b, c, e_real, e, f, g]

    return new_state, cc_T1, cc_T2, cc_e, cc_a

# --- Carry chain length measurement ---
def max_carry_chain(cc):
    """Find the longest consecutive run of 1-bits in cc (= longest carry chain)."""
    if cc == 0:
        return 0
    max_run = 0
    run = 0
    for i in range(32):
        if (cc >> i) & 1:
            run += 1
            if run > max_run:
                max_run = run
        else:
            run = 0
    return max_run

def avg_carry_chain(cc):
    """Average length of carry chain runs in cc."""
    if cc == 0:
        return 0.0
    runs = []
    run = 0
    for i in range(32):
        if (cc >> i) & 1:
            run += 1
        else:
            if run > 0:
                runs.append(run)
            run = 0
    if run > 0:
        runs.append(run)
    return sum(runs) / len(runs) if runs else 0.0


def random_w16():
    return [int.from_bytes(os.urandom(4), 'big') for _ in range(16)]


# ================================================================
# EXPERIMENT 1 & 2: Per-bit carry profile + heatmap
# ================================================================
def experiment_bitwise_profile(N=1000):
    """
    For each round (0..63) and bit position (0..31), measure P(cc[bit]=1)
    for all four carry sources: cc_T1, cc_T2, cc_e, cc_a.
    """
    # Accumulators: [round][bit] = count of times that bit was 1
    profile_T1 = [[0]*32 for _ in range(64)]
    profile_T2 = [[0]*32 for _ in range(64)]
    profile_e  = [[0]*32 for _ in range(64)]
    profile_a  = [[0]*32 for _ in range(64)]

    # Carry chain stats: [round] = (total_max, total_avg, count)
    chain_max_acc = [0]*64
    chain_avg_acc = [0.0]*64

    # Bit-0 violation counter
    bit0_violations = 0

    for trial in range(N):
        W = random_w16()
        Wexp = expand_schedule(W)
        state = list(IV)

        for r in range(64):
            state, cc_T1, cc_T2, cc_e, cc_a = round_carries(state, Wexp[r], K[r])

            # Per-bit accumulation
            for b in range(32):
                if bit(cc_T1, b): profile_T1[r][b] += 1
                if bit(cc_T2, b): profile_T2[r][b] += 1
                if bit(cc_e, b):  profile_e[r][b]  += 1
                if bit(cc_a, b):  profile_a[r][b]  += 1

            # Check F3: bit 0 must always be 0
            if bit(cc_T1, 0) or bit(cc_T2, 0) or bit(cc_e, 0) or bit(cc_a, 0):
                bit0_violations += 1

            # Carry chain stats (combined cc = cc_T1 for the 5-operand addition)
            chain_max_acc[r] += max_carry_chain(cc_T1)
            chain_avg_acc[r] += avg_carry_chain(cc_T1)

    return profile_T1, profile_T2, profile_e, profile_a, chain_max_acc, chain_avg_acc, bit0_violations


def print_bitwise_profile(profile_T1, profile_T2, profile_e, profile_a, N):
    """Print per-bit carry probabilities."""
    print("=" * 90)
    print("EXPERIMENT 1: Per-bit carry probability P(cc[bit]=1)")
    print("=" * 90)
    print()

    # Print summary for representative rounds and all 4 carry sources
    rep_rounds = [0, 1, 2, 3, 5, 10, 20, 40, 63]
    sources = [
        ("cc_T1", profile_T1),
        ("cc_T2", profile_T2),
        ("cc_e",  profile_e),
        ("cc_a",  profile_a),
    ]

    for src_name, profile in sources:
        print(f"\n--- {src_name} ---")
        print(f"{'r':>3} | bit0  bit1  bit2  bit3  bit4  bit7  bit15 bit23 bit31 | avg")
        print("-" * 80)

        sel_bits = [0, 1, 2, 3, 4, 7, 15, 23, 31]

        for r in rep_rounds:
            vals = [profile[r][b] / N for b in sel_bits]
            avg_all = sum(profile[r][b] for b in range(32)) / (32 * N)
            parts = " ".join(f"{v:.3f}" for v in vals)
            print(f"{r:3d} | {parts} | {avg_all:.3f}")

    # Verify F3: bit 0 always carry-free
    print()
    print("--- F3 Verification: bit 0 carry-free ---")
    all_ok = True
    for r in range(64):
        for src_name, profile in sources:
            if profile[r][0] != 0:
                print(f"  VIOLATION at round {r}, {src_name}: P(cc[0]=1) = {profile[r][0]/N:.4f}")
                all_ok = False
    if all_ok:
        print("  CONFIRMED: bit 0 = 0 across all rounds, all sources, all trials.")


def print_heatmap(profile_T1, N):
    """
    Print a 64x32 heatmap of P(cc_T1[bit]=1) using compact characters.
    Character encoding: ' '=0, .=<0.2, o=<0.4, O=<0.6, #=>=0.6
    """
    print()
    print("=" * 90)
    print("EXPERIMENT 2: Carry heatmap (rounds x bits) for cc_T1")
    print("  Key: ' '=0  .=<0.20  o=<0.40  O=<0.60  #=>=0.60")
    print("=" * 90)

    # Column header
    hdr = "r\\b|"
    for b in range(32):
        hdr += str(b % 10)
    hdr += "| (bit 0=LSB, 31=MSB)"
    print(hdr)
    print("---+" + "-"*32 + "+")

    for r in range(64):
        row = f"{r:2d} |"
        for b in range(32):
            p = profile_T1[r][b] / N
            if p == 0:
                ch = ' '
            elif p < 0.20:
                ch = '.'
            elif p < 0.40:
                ch = 'o'
            elif p < 0.60:
                ch = 'O'
            else:
                ch = '#'
            row += ch
        # Compute average probability for this round
        avg_p = sum(profile_T1[r][b] for b in range(32)) / (32 * N)
        row += f"| avg={avg_p:.3f}"
        print(row)

    print("---+" + "-"*32 + "+")

    # Per-bit average across all rounds
    print("\nPer-bit average P(cc_T1[bit]=1) across rounds 5-63 (steady state):")
    print("bit: ", end="")
    for b in range(32):
        print(f"{b:5d}", end="")
    print()
    print("P:   ", end="")
    for b in range(32):
        avg = sum(profile_T1[r][b] for r in range(5, 64)) / (59 * N)
        print(f"{avg:5.3f}", end="")
    print()

    # Carry wavefront analysis
    print("\nCarry wavefront analysis (first round where P(cc_T1[bit]>0.1)):")
    for b in [0, 1, 2, 3, 4, 5, 7, 10, 15, 20, 25, 31]:
        for r in range(64):
            if profile_T1[r][b] / N > 0.1:
                print(f"  bit {b:2d}: first active at round {r}")
                break
        else:
            print(f"  bit {b:2d}: never exceeds 0.1")


# ================================================================
# EXPERIMENT 3: Carry chain length per round
# ================================================================
def print_chain_lengths(chain_max_acc, chain_avg_acc, N):
    print()
    print("=" * 70)
    print("EXPERIMENT 3: Carry chain length per round (cc_T1)")
    print("=" * 70)
    print(f"{'r':>3} | {'avg_max':>8} {'avg_chain':>10} | bar (max chain)")
    print("-" * 60)

    for r in range(64):
        avg_max = chain_max_acc[r] / N
        avg_avg = chain_avg_acc[r] / N
        bar = '*' * int(avg_max + 0.5)

        if r < 10 or r % 4 == 0 or r == 63:
            print(f"{r:3d} | {avg_max:8.2f} {avg_avg:10.2f} | {bar}")

    # Steady-state averages
    ss_max = sum(chain_max_acc[10:64]) / (54 * N)
    ss_avg = sum(chain_avg_acc[10:64]) / (54 * N)
    print(f"\nSteady-state (r=10..63): avg max_chain = {ss_max:.2f}, avg avg_chain = {ss_avg:.2f}")
    print(f"Theoretical P(carry survival) = 5/8 = 0.625")
    print(f"Expected chain length (geometric) = 1/(1-5/8) = {1/(1-5/8):.2f}")


# ================================================================
# EXPERIMENT 4: Differential carry correlation
# ================================================================
def experiment_diff_carry(N=1000):
    """
    For pairs (W, W^delta), measure correlation of cc[bit_i] between
    the two messages. This tells us if carry at bit i 'remembers' the
    input difference.
    """
    # correlation[round][bit] = count of agreements (both 0 or both 1)
    agree = [[0]*32 for _ in range(64)]
    total = [[0]*32 for _ in range(64)]

    # Per-bit P(delta_cc[bit]=0) — probability that carry bit is SAME for both messages
    same_count = [[0]*32 for _ in range(64)]

    for trial in range(N):
        W1 = random_w16()
        W2 = list(W1)
        # Single-bit flip in W[0]
        flip_bit = trial % 32
        W2[0] ^= (1 << flip_bit)

        Wexp1 = expand_schedule(W1)
        Wexp2 = expand_schedule(W2)

        state1 = list(IV)
        state2 = list(IV)

        for r in range(64):
            state1, cc1_T1, _, _, _ = round_carries(state1, Wexp1[r], K[r])
            state2, cc2_T1, _, _, _ = round_carries(state2, Wexp2[r], K[r])

            diff_cc = cc1_T1 ^ cc2_T1

            for b in range(32):
                total[r][b] += 1
                if bit(diff_cc, b) == 0:
                    same_count[r][b] += 1

    return same_count, total


def print_diff_carry(same_count, total):
    print()
    print("=" * 90)
    print("EXPERIMENT 4: Differential carry correlation")
    print("  P(cc1[bit] == cc2[bit]) for (W, W^delta), delta = single-bit flip")
    print("  Value 0.50 = uncorrelated, >0.50 = correlated (carry remembers)")
    print("=" * 90)

    # Compact heatmap
    print(f"\n{'r':>3} | bit0  bit1  bit2  bit3  bit4  bit7  bit15 bit23 bit31 | avg")
    print("-" * 80)

    sel_bits = [0, 1, 2, 3, 4, 7, 15, 23, 31]
    rep_rounds = [0, 1, 2, 3, 4, 5, 7, 10, 15, 20, 30, 40, 50, 63]

    for r in rep_rounds:
        vals = []
        for b in sel_bits:
            p = same_count[r][b] / total[r][b] if total[r][b] > 0 else 0
            vals.append(p)
        avg_all = sum(same_count[r][b] for b in range(32)) / sum(total[r][b] for b in range(32))
        parts = " ".join(f"{v:.3f}" for v in vals)
        # Mark if significantly correlated
        marker = ""
        if avg_all > 0.55:
            marker = " << CORRELATED"
        elif avg_all > 0.52:
            marker = " < weak"
        print(f"{r:3d} | {parts} | {avg_all:.3f}{marker}")

    # Low-bit vs high-bit comparison
    print("\nLow bits (0-7) vs high bits (24-31) correlation across rounds:")
    print(f"{'r':>3} | low_avg  high_avg | gap")
    print("-" * 45)

    for r in rep_rounds:
        low_avg = sum(same_count[r][b] / total[r][b] for b in range(8)) / 8
        high_avg = sum(same_count[r][b] / total[r][b] for b in range(24, 32)) / 8
        gap = low_avg - high_avg
        marker = " << low bits more predictable" if gap > 0.02 else ""
        print(f"{r:3d} | {low_avg:.4f}  {high_avg:.4f} | {gap:+.4f}{marker}")


# ================================================================
# EXPERIMENT 5: Predictability window detection
# ================================================================
def experiment_predictability(N=1000):
    """
    For each bit position, find the round range where carry is 'predictable'
    (significantly different from 0.5 probability).
    """
    # Collect per-bit carry probability per round
    counts = [[0]*32 for _ in range(64)]

    for trial in range(N):
        W = random_w16()
        Wexp = expand_schedule(W)
        state = list(IV)

        for r in range(64):
            state, cc_T1, _, _, _ = round_carries(state, Wexp[r], K[r])
            for b in range(32):
                if bit(cc_T1, b):
                    counts[r][b] += 1

    return counts


def print_predictability(counts, N):
    print()
    print("=" * 90)
    print("EXPERIMENT 5: Predictability window — bit-position-dependent")
    print("  'Predictable' = |P(cc[bit]=1) - 0.5| > threshold")
    print("  Threshold = 0.05 (significant deviation from fair coin)")
    print("=" * 90)

    threshold = 0.05

    print(f"\n{'bit':>4} | {'window':>12} | {'min_P':>6} {'max_P':>6} {'avg_P(r5-63)':>12} | notes")
    print("-" * 80)

    for b in range(32):
        probs = [counts[r][b] / N for r in range(64)]

        # Find first round where |P - 0.5| < threshold (carry becomes unpredictable)
        predictable_end = 0
        for r in range(64):
            if abs(probs[r] - 0.5) > threshold:
                predictable_end = r + 1
            else:
                # Once unpredictable, check if it stays unpredictable
                stays_unpred = True
                for r2 in range(r, min(r + 5, 64)):
                    if abs(probs[r2] - 0.5) > threshold:
                        stays_unpred = False
                        break
                if stays_unpred:
                    break

        # Steady-state average
        avg_ss = sum(probs[5:64]) / 59 if len(probs) > 5 else 0
        min_p = min(probs)
        max_p = max(probs)

        notes = ""
        if b == 0:
            notes = "ALWAYS 0 (carry-free LSB)"
        elif avg_ss < 0.40:
            notes = "below 0.5 — asymmetric"
        elif avg_ss > 0.55:
            notes = "above 0.5 — carry-heavy"

        # Window display
        if b == 0:
            window = "all rounds"
        elif predictable_end <= 5:
            window = "r=0-4 only"
        else:
            window = f"r=0-{predictable_end-1}"

        print(f"{b:4d} | {window:>12} | {min_p:6.3f} {max_p:6.3f} {avg_ss:12.3f} | {notes}")

    # Summary: bit-position dependent structure
    print()
    print("--- Degree analysis ---")
    print("ALG theory predicts: deg(cc[bit_i]) = i+1")
    print("Bit 0: degree 1 (always 0, trivially predictable)")
    print("Bit 1: degree 2 (quadratic, low complexity)")
    print("Bit k: degree k+1 (complexity grows with position)")
    print()
    print("Implication: low-order bits of carry are algebraically simpler")
    print("and should remain predictable for more rounds.")


# ================================================================
# MAIN
# ================================================================
if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 1000

    print(f"SCF Carry Bitwise Decomposition (N={N} trials)")
    print(f"Analyzing 64 rounds x 32 bit positions x 4 carry sources")
    print()

    # --- Experiment 1 & 2: Per-bit profile + heatmap ---
    print("Running per-bit carry profile (Experiments 1, 2, 3)...")
    prof_T1, prof_T2, prof_e, prof_a, chain_max, chain_avg, bit0_viol = \
        experiment_bitwise_profile(N)

    print(f"  Bit-0 violations: {bit0_viol}")
    print()

    print_bitwise_profile(prof_T1, prof_T2, prof_e, prof_a, N)
    print_heatmap(prof_T1, N)

    # --- Experiment 3: Chain lengths ---
    print_chain_lengths(chain_max, chain_avg, N)

    # --- Experiment 4: Differential carry correlation ---
    print("\nRunning differential carry correlation (Experiment 4)...")
    same_count, total = experiment_diff_carry(N)
    print_diff_carry(same_count, total)

    # --- Experiment 5: Predictability window ---
    print("\nRunning predictability window analysis (Experiment 5)...")
    pred_counts = experiment_predictability(N)
    print_predictability(pred_counts, N)

    # ================================================================
    # FINAL SUMMARY
    # ================================================================
    print()
    print("=" * 90)
    print("SUMMARY: Key findings")
    print("=" * 90)

    # 1. F3 verification
    print(f"\n1. F3 (carry-free LSB): {'CONFIRMED' if bit0_viol == 0 else 'VIOLATED (' + str(bit0_viol) + ' violations)'}")

    # 2. Steady-state carry probability by bit position
    print("\n2. Steady-state P(cc_T1[bit]=1) by bit position (r=10..63):")
    print("   Low bits (0-3):  ", end="")
    for b in range(4):
        avg = sum(prof_T1[r][b] for r in range(10, 64)) / (54 * N)
        print(f"bit{b}={avg:.3f} ", end="")
    print()
    print("   Mid bits (14-17):", end="")
    for b in range(14, 18):
        avg = sum(prof_T1[r][b] for r in range(10, 64)) / (54 * N)
        print(f"bit{b}={avg:.3f} ", end="")
    print()
    print("   High bits (28-31):", end="")
    for b in range(28, 32):
        avg = sum(prof_T1[r][b] for r in range(10, 64)) / (54 * N)
        print(f"bit{b}={avg:.3f} ", end="")
    print()

    # 3. Carry chain
    ss_max = sum(chain_max[10:64]) / (54 * N)
    print(f"\n3. Average max carry chain (r=10..63): {ss_max:.2f} bits")
    print(f"   Theory: geometric(p=5/8) => E[max] grows with operand count")

    # 4. Differential correlation
    early_corr = sum(same_count[r][b] for r in range(5) for b in range(32)) / \
                 sum(total[r][b] for r in range(5) for b in range(32))
    late_corr = sum(same_count[r][b] for r in range(20, 64) for b in range(32)) / \
                sum(total[r][b] for r in range(20, 64) for b in range(32))
    print(f"\n4. Differential carry correlation:")
    print(f"   Early rounds (0-4): P(same) = {early_corr:.4f}")
    print(f"   Late rounds (20-63): P(same) = {late_corr:.4f}")
    print(f"   Decorrelation gap: {early_corr - late_corr:+.4f}")

    # 5. Key insight
    print(f"\n5. KEY QUESTION: Is there a bit-position-dependent predictability window?")
    print(f"   Answer: bit 0 is ALWAYS predictable (= 0).")
    low_bits_avg = sum(prof_T1[r][1] + prof_T1[r][2] for r in range(5)) / (2 * 5 * N)
    high_bits_avg = sum(prof_T1[r][30] + prof_T1[r][31] for r in range(5)) / (2 * 5 * N)
    print(f"   Early rounds (0-4): low bits P={low_bits_avg:.3f}, high bits P={high_bits_avg:.3f}")
    if abs(low_bits_avg - high_bits_avg) > 0.02:
        print(f"   => Low bits and high bits behave differently in early rounds!")
        print(f"   => Exploitable structure exists in bit-position space.")
    else:
        print(f"   => Similar behavior across bit positions in early rounds.")

    print()
