"""
Carry BIT PATTERN analysis for SHA-256.

Previous experiment (carry_profile_space.py) looked at carry COUNT per ADD.
This goes deeper: which specific bits carry? Are certain bit positions
deterministic or highly biased?

Tracks full 32-bit carry chains for key ADDs across rounds, measures
per-bit probabilities, differential patterns under delta-W, and
cross-round correlations.
"""

import numpy as np
import struct
import os

MASK32 = 0xFFFFFFFF
NUM_SAMPLES = 5000

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
IV = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]

def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def shr(x, n): return x >> n
def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ shr(x, 3)
def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ shr(x, 10)
def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ ((~e) & g) & MASK32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)


def add32_carry_bits(x, y):
    """Return (result, carry_bitvector) for x + y mod 2^32.
    carry_bitvector bit j = 1 means there was a carry OUT of position j."""
    result = (x + y) & MASK32
    carry_vec = 0
    c = 0
    for i in range(32):
        s = ((x >> i) & 1) + ((y >> i) & 1) + c
        c = s >> 1
        if c:
            carry_vec |= (1 << i)
    return result, carry_vec


def add32(x, y):
    return (x + y) & MASK32


def bitvec_to_array(v):
    """Convert 32-bit int to numpy array of bits."""
    return np.array([(v >> i) & 1 for i in range(32)], dtype=np.int8)


# ─── Key ADD labels ───
# In each round's compression step, the 7 ADDs are:
#   0: T1_step1 = h + Sigma1(e)
#   1: T1_step2 = T1_step1 + Ch(e,f,g)
#   2: T1_step3 = T1_step2 + K[r]
#   3: T1_final = T1_step3 + W[r]
#   4: T2       = Sigma0(a) + Maj(a,b,c)
#   5: e_new    = d + T1
#   6: a_new    = T1 + T2
ADD_NAMES = ["h+Sig1", "..+Ch", "..+K[r]", "..+W[r]", "Sig0+Maj", "d+T1", "T1+T2"]
KEY_ADDS = [3, 5, 6]  # T1_final, e_new, a_new


def sha256_round_carry_bits(W_input):
    """Run SHA-256 compression, return carry bit vectors for all 64x7 ADDs.
    Returns: carry_bits[round][add_idx] as 32-bit ints.
    """
    W = list(W_input[:16])
    for r in range(16, 64):
        s0 = sigma0(W[r-15])
        s1 = sigma1(W[r-2])
        t, _ = add32_carry_bits(s1, W[r-7])
        t, _ = add32_carry_bits(t, s0)
        t, _ = add32_carry_bits(t, W[r-16])
        W.append(t)

    a, b, c, d, e, f, g, h = IV

    carry_bits = []  # [round][7 adds]
    for r in range(64):
        round_carries = []

        # ADD 0: h + Sigma1(e)
        sig1_e = Sigma1(e)
        t1, cv0 = add32_carry_bits(h, sig1_e)
        round_carries.append(cv0)

        # ADD 1: + Ch(e,f,g)
        ch_val = Ch(e, f, g)
        t1, cv1 = add32_carry_bits(t1, ch_val)
        round_carries.append(cv1)

        # ADD 2: + K[r]
        t1, cv2 = add32_carry_bits(t1, K[r])
        round_carries.append(cv2)

        # ADD 3: + W[r]  (T1 final)
        t1, cv3 = add32_carry_bits(t1, W[r])
        round_carries.append(cv3)

        # ADD 4: Sigma0(a) + Maj(a,b,c)
        sig0_a = Sigma0(a)
        maj_val = Maj(a, b, c)
        t2, cv4 = add32_carry_bits(sig0_a, maj_val)
        round_carries.append(cv4)

        # ADD 5: e' = d + T1
        e_new, cv5 = add32_carry_bits(d, t1)
        round_carries.append(cv5)

        # ADD 6: a' = T1 + T2
        a_new, cv6 = add32_carry_bits(t1, t2)
        round_carries.append(cv6)

        carry_bits.append(round_carries)

        # Shift state
        h, g, f, e = g, f, e, e_new
        d, c, b, a = c, b, a, a_new

    return carry_bits


def random_W():
    data = os.urandom(64)
    return list(struct.unpack('>16I', data))


def main():
    np.set_printoptions(precision=4, linewidth=120, suppress=True)
    print("=" * 80)
    print("CARRY BIT PATTERN ANALYSIS — SHA-256")
    print(f"Samples: {NUM_SAMPLES}")
    print("=" * 80)

    # ─────────────────────────────────────────────────────
    # Section 1: Per-bit carry probability for key ADDs
    # ─────────────────────────────────────────────────────
    print("\n" + "=" * 80)
    print("SECTION 1: Per-bit carry probability P(carry[j]=1)")
    print("=" * 80)

    # Accumulate carry bit counts: [64 rounds][7 adds][32 bits]
    carry_accum = np.zeros((64, 7, 32), dtype=np.float64)

    # For differential analysis (section 3)
    delta_carry_accum = np.zeros((64, 7, 32), dtype=np.float64)

    # For cross-round correlation (section 5)
    # Store raw carry bit arrays for key adds
    raw_carries = {aid: np.zeros((NUM_SAMPLES, 64, 32), dtype=np.int8) for aid in KEY_ADDS}

    delta_W = 0x80000000  # flip MSB of W[0]

    for s in range(NUM_SAMPLES):
        W1 = random_W()
        cb1 = sha256_round_carry_bits(W1)

        W2 = list(W1)
        W2[0] ^= delta_W
        cb2 = sha256_round_carry_bits(W2)

        for r in range(64):
            for a_idx in range(7):
                bits1 = bitvec_to_array(cb1[r][a_idx])
                carry_accum[r, a_idx] += bits1

                bits_diff = bitvec_to_array(cb1[r][a_idx] ^ cb2[r][a_idx])
                delta_carry_accum[r, a_idx] += bits_diff

                if a_idx in KEY_ADDS:
                    raw_carries[a_idx][s, r] = bits1

    carry_prob = carry_accum / NUM_SAMPLES
    delta_prob = delta_carry_accum / NUM_SAMPLES

    # Print per-bit probabilities for key ADDs in selected rounds
    focus_rounds = [0, 1, 2, 3, 60, 61, 62, 63]
    for a_idx in KEY_ADDS:
        print(f"\n--- ADD {a_idx}: {ADD_NAMES[a_idx]} ---")
        for r in focus_rounds:
            probs = carry_prob[r, a_idx]
            # Find strongly biased bits (close to 0 or 1)
            near_zero = np.sum(probs < 0.05)
            near_one = np.sum(probs > 0.95)
            near_half = np.sum((probs > 0.4) & (probs < 0.6))
            biased = np.sum((probs < 0.1) | (probs > 0.9))

            print(f"  Round {r:2d}: "
                  f"near_0={near_zero:2d}  near_1={near_one:2d}  "
                  f"near_0.5={near_half:2d}  strongly_biased(>0.9 or <0.1)={biased:2d}")
            # Show the actual probabilities for bits 0-15 and 16-31
            print(f"    bits[ 0:15] = {probs[0:16]}")
            print(f"    bits[16:31] = {probs[16:32]}")

    # ─────────────────────────────────────────────────────
    # Section 2: Summary of bias across ALL rounds
    # ─────────────────────────────────────────────────────
    print("\n" + "=" * 80)
    print("SECTION 2: Bias summary across all rounds for key ADDs")
    print("=" * 80)

    for a_idx in KEY_ADDS:
        print(f"\n--- ADD {a_idx}: {ADD_NAMES[a_idx]} ---")
        print(f"  {'Round':>5s}  {'#det(p<.01|>.99)':>17s}  {'#biased(p<.1|>.9)':>18s}  "
              f"{'#unbiased(.4-.6)':>16s}  {'mean(p)':>7s}  {'std(p)':>7s}")
        for r in range(64):
            probs = carry_prob[r, a_idx]
            det = np.sum((probs < 0.01) | (probs > 0.99))
            biased = np.sum((probs < 0.1) | (probs > 0.9))
            unbiased = np.sum((probs > 0.4) & (probs < 0.6))
            print(f"  {r:5d}  {det:17d}  {biased:18d}  {unbiased:16d}  "
                  f"{np.mean(probs):7.4f}  {np.std(probs):7.4f}")

    # ─────────────────────────────────────────────────────
    # Section 3: Differential carry bit patterns
    # ─────────────────────────────────────────────────────
    print("\n" + "=" * 80)
    print("SECTION 3: Differential carry bit patterns (deltaW0 = 0x80000000)")
    print("Which carry BITS change between W1 and W2?")
    print("=" * 80)

    for a_idx in KEY_ADDS:
        print(f"\n--- ADD {a_idx}: {ADD_NAMES[a_idx]} ---")
        for r in focus_rounds:
            dp = delta_prob[r, a_idx]
            always_same = np.sum(dp < 0.01)
            always_diff = np.sum(dp > 0.99)
            sometimes = 32 - always_same - always_diff
            print(f"  Round {r:2d}: always_same={always_same:2d}  always_diff={always_diff:2d}  "
                  f"variable={sometimes:2d}")
            print(f"    P(delta_carry[j]=1): bits[ 0:15] = {dp[0:16]}")
            print(f"                         bits[16:31] = {dp[16:32]}")

    # Aggregate: how many bits never change across all rounds?
    print("\n--- Aggregate: bits where delta_carry is ALWAYS 0 across all 64 rounds ---")
    for a_idx in KEY_ADDS:
        all_zero = np.ones(32, dtype=bool)
        for r in range(64):
            all_zero &= (delta_prob[r, a_idx] < 0.01)
        print(f"  ADD {a_idx} ({ADD_NAMES[a_idx]}): {np.sum(all_zero)} bits always zero "
              f"(positions: {np.where(all_zero)[0].tolist()})")

    # ─────────────────────────────────────────────────────
    # Section 4: "46 bits deterministic" check for e_add in round 0
    # ─────────────────────────────────────────────────────
    print("\n" + "=" * 80)
    print("SECTION 4: Deterministic carry bits in round 0 (fixed IV)")
    print("How many carry bits are determined by W alone?")
    print("=" * 80)

    for a_idx in range(7):
        probs = carry_prob[0, a_idx]
        det_strict = np.sum((probs < 0.001) | (probs > 0.999))
        det_loose = np.sum((probs < 0.01) | (probs > 0.99))
        biased = np.sum((probs < 0.1) | (probs > 0.9))

        # Identify which bits are deterministic
        det_0_bits = np.where(probs < 0.001)[0].tolist()
        det_1_bits = np.where(probs > 0.999)[0].tolist()

        print(f"  ADD {a_idx} ({ADD_NAMES[a_idx]:>8s}): "
              f"deterministic(p<.001|>.999)={det_strict:2d}  "
              f"near-det(p<.01|>.99)={det_loose:2d}  "
              f"biased(p<.1|>.9)={biased:2d}")
        if det_0_bits:
            print(f"    Always 0: bits {det_0_bits}")
        if det_1_bits:
            print(f"    Always 1: bits {det_1_bits}")

    # Across all 7 ADDs in round 0: total deterministic bits
    total_det = 0
    total_biased = 0
    for a_idx in range(7):
        probs = carry_prob[0, a_idx]
        total_det += np.sum((probs < 0.001) | (probs > 0.999))
        total_biased += np.sum((probs < 0.01) | (probs > 0.99))

    print(f"\n  TOTAL across 7 ADDs in round 0:")
    print(f"    Deterministic (p<.001 or p>.999): {total_det} / {7*32} = {7*32} possible")
    print(f"    Near-deterministic (p<.01 or p>.99): {total_biased} / {7*32}")

    # Also check: how many across round 0 + round 1?
    total_det_r01 = 0
    for r in range(2):
        for a_idx in range(7):
            probs = carry_prob[r, a_idx]
            total_det_r01 += np.sum((probs < 0.001) | (probs > 0.999))
    print(f"    Deterministic in rounds 0-1: {total_det_r01} / {2*7*32}")

    # ─────────────────────────────────────────────────────
    # Section 5: Cross-round carry bit correlation
    # ─────────────────────────────────────────────────────
    print("\n" + "=" * 80)
    print("SECTION 5: Cross-round carry bit correlations")
    print("Correlation of carry[bit_j, round_r] with carry[bit_j, round_r+1]")
    print("=" * 80)

    for a_idx in KEY_ADDS:
        print(f"\n--- ADD {a_idx}: {ADD_NAMES[a_idx]} ---")
        data = raw_carries[a_idx]  # [samples, 64, 32]

        # Compute per-bit cross-round correlation for consecutive rounds
        corr_matrix = np.zeros((63, 32))
        for r in range(63):
            for j in range(32):
                v1 = data[:, r, j].astype(np.float64)
                v2 = data[:, r+1, j].astype(np.float64)
                std1 = np.std(v1)
                std2 = np.std(v2)
                if std1 > 0.01 and std2 > 0.01:
                    corr_matrix[r, j] = np.corrcoef(v1, v2)[0, 1]
                else:
                    corr_matrix[r, j] = np.nan  # one is constant

        # Summary stats
        valid = ~np.isnan(corr_matrix)
        if np.any(valid):
            abs_corr = np.abs(corr_matrix[valid])
            print(f"  Overall: mean|corr|={np.mean(abs_corr):.4f}  "
                  f"max|corr|={np.max(abs_corr):.4f}  "
                  f"fraction |corr|>0.1: {np.mean(abs_corr > 0.1):.4f}")

        # Show rounds with highest correlation
        print(f"  Top cross-round correlations (bit j, round r -> r+1):")
        flat_idx = np.argsort(np.abs(np.nan_to_num(corr_matrix)).ravel())[::-1]
        shown = 0
        for idx in flat_idx[:15]:
            r_idx = idx // 32
            j_idx = idx % 32
            val = corr_matrix[r_idx, j_idx]
            if not np.isnan(val) and abs(val) > 0.05:
                print(f"    bit {j_idx:2d}, round {r_idx:2d}->{r_idx+1:2d}: corr={val:+.4f}")
                shown += 1
        if shown == 0:
            print(f"    No significant correlations found (all |corr| < 0.05)")

        # Per-bit average correlation across all rounds
        mean_corr_per_bit = np.nanmean(np.abs(corr_matrix), axis=0)
        print(f"  Mean |corr| by bit position:")
        print(f"    bits[ 0:15] = {mean_corr_per_bit[0:16]}")
        print(f"    bits[16:31] = {mean_corr_per_bit[16:32]}")

    # ─────────────────────────────────────────────────────
    # Section 6: Carry bit entropy analysis
    # ─────────────────────────────────────────────────────
    print("\n" + "=" * 80)
    print("SECTION 6: Per-round carry bit entropy (bits of uncertainty)")
    print("Binary entropy H(p) for each carry bit, summed over 32 bits per ADD")
    print("Max possible: 32 bits (all carry bits perfectly random)")
    print("=" * 80)

    def binary_entropy(p):
        """H(p) = -p*log2(p) - (1-p)*log2(1-p), with 0*log(0)=0."""
        p = np.clip(p, 1e-12, 1 - 1e-12)
        return -p * np.log2(p) - (1 - p) * np.log2(1 - p)

    for a_idx in KEY_ADDS:
        print(f"\n--- ADD {a_idx}: {ADD_NAMES[a_idx]} ---")
        print(f"  {'Round':>5s}  {'entropy(bits)':>13s}  {'max_entropy':>11s}  {'det_bits':>8s}")
        for r in range(64):
            probs = carry_prob[r, a_idx]
            ent_per_bit = binary_entropy(probs)
            total_ent = np.sum(ent_per_bit)
            det = np.sum((probs < 0.01) | (probs > 0.99))
            print(f"  {r:5d}  {total_ent:13.4f}  {32.0:11.1f}  {det:8d}")

    # ─────────────────────────────────────────────────────
    # Final Summary
    # ─────────────────────────────────────────────────────
    print("\n" + "=" * 80)
    print("FINAL SUMMARY & INTERPRETATION")
    print("=" * 80)

    # Count total deterministic bits across all rounds and adds
    total_det_all = 0
    total_biased_all = 0
    total_bits = 64 * 7 * 32
    for r in range(64):
        for a_idx in range(7):
            probs = carry_prob[r, a_idx]
            total_det_all += np.sum((probs < 0.01) | (probs > 0.99))
            total_biased_all += np.sum((probs < 0.1) | (probs > 0.9))

    print(f"\n1. DETERMINISTIC CARRY BITS:")
    print(f"   Total near-deterministic (p<.01 or p>.99): {total_det_all} / {total_bits} "
          f"({100*total_det_all/total_bits:.1f}%)")
    print(f"   Total strongly biased (p<.1 or p>.9): {total_biased_all} / {total_bits} "
          f"({100*total_biased_all/total_bits:.1f}%)")

    # Early rounds vs late rounds
    for label, rng in [("Rounds 0-3", range(4)), ("Rounds 60-63", range(60, 64))]:
        det = 0
        for r in rng:
            for a_idx in range(7):
                probs = carry_prob[r, a_idx]
                det += np.sum((probs < 0.01) | (probs > 0.99))
        total = len(list(rng)) * 7 * 32
        print(f"   {label}: {det}/{total} deterministic ({100*det/total:.1f}%)")

    # Differential summary
    print(f"\n2. DIFFERENTIAL CARRY BITS (deltaW0 = MSB flip):")
    for a_idx in KEY_ADDS:
        total_always_same = 0
        total_always_diff = 0
        for r in range(64):
            dp = delta_prob[r, a_idx]
            total_always_same += np.sum(dp < 0.01)
            total_always_diff += np.sum(dp > 0.99)
        total = 64 * 32
        print(f"   ADD {a_idx} ({ADD_NAMES[a_idx]}): "
              f"always_same={total_always_same}/{total} ({100*total_always_same/total:.1f}%)  "
              f"always_diff={total_always_diff}/{total} ({100*total_always_diff/total:.1f}%)")

    # Cross-round correlation summary
    print(f"\n3. CROSS-ROUND CORRELATION:")
    for a_idx in KEY_ADDS:
        data = raw_carries[a_idx]
        all_corrs = []
        for r in range(63):
            for j in range(32):
                v1 = data[:, r, j].astype(np.float64)
                v2 = data[:, r+1, j].astype(np.float64)
                if np.std(v1) > 0.01 and np.std(v2) > 0.01:
                    all_corrs.append(abs(np.corrcoef(v1, v2)[0, 1]))
        if all_corrs:
            arr = np.array(all_corrs)
            print(f"   ADD {a_idx} ({ADD_NAMES[a_idx]}): "
                  f"mean|corr|={np.mean(arr):.4f}  "
                  f"max|corr|={np.max(arr):.4f}  "
                  f"|corr|>0.1: {np.mean(arr > 0.1)*100:.1f}%  "
                  f"|corr|>0.3: {np.mean(arr > 0.3)*100:.1f}%")

    # Entropy summary
    print(f"\n4. ENTROPY SUMMARY (carry uncertainty per ADD):")
    for a_idx in KEY_ADDS:
        ents = []
        for r in range(64):
            probs = carry_prob[r, a_idx]
            ent = np.sum(binary_entropy(probs))
            ents.append(ent)
        ents = np.array(ents)
        print(f"   ADD {a_idx} ({ADD_NAMES[a_idx]}): "
              f"mean={np.mean(ents):.2f}/32  "
              f"min={np.min(ents):.2f} (round {np.argmin(ents)})  "
              f"max={np.max(ents):.2f} (round {np.argmax(ents)})")

    print(f"\n5. KEY INSIGHT — Carry bit patterns vs carry counts:")
    print(f"   Previous experiment: 444/448 dimensions in carry count space.")
    print(f"   Bit-level view: each of those 448 counts decomposes into 32 bit decisions.")
    print(f"   Total bit-level space: {total_bits} dimensions.")
    print(f"   Near-deterministic: {total_det_all} ({100*total_det_all/total_bits:.1f}%)")
    print(f"   This means ~{total_bits - total_det_all} carry BITS are free,")
    print(f"   compared to ~444 carry COUNTS being free.")
    print(f"   Ratio: {(total_bits - total_det_all)/444:.1f}x more structure visible at bit level.")

    print("\nDone.")


if __name__ == "__main__":
    main()
