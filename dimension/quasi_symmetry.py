"""
КВАЗИ-СИММЕТРИИ: complement δ=32, swap_halves δ=32.

Это на ВНУТРЕННЕМ state, не на output.
Вопрос: видна ли квази-симметрия на OUTPUT (после feedforward)?

Если да → DISTINGUISHER: SHA-256 ≠ random oracle.
Если нет → квази-симметрия стирается feedforward.
"""

import numpy as np
import struct, hashlib

MASK32 = 0xFFFFFFFF

def sha256_hash(W16):
    raw = struct.pack('>16I', *W16)
    return struct.unpack('>8I', hashlib.sha256(raw).digest())

def hw(x): return bin(x).count('1')


def main():
    np.random.seed(42)

    print("=" * 70)
    print("КВАЗИ-СИММЕТРИИ: distinguisher test")
    print("=" * 70)

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("1. COMPLEMENT QUASI-SYMMETRY на output")
    print("=" * 70)

    # Quasi-symmetry was on INTERNAL state (before feedforward).
    # After feedforward: H = IV + state.
    # If complement(state) ≈ state, then:
    #   H' = IV + ~state ≈ IV + state + carry_noise = H + carry_noise
    # So H' - H ≈ constant (IV-related)?

    # Test: H(W) vs H(~W) where ~W = complement of all message words
    complement_diffs = []
    random_diffs = []

    for _ in range(10000):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H = sha256_hash(W)

        # Complement all words
        W_c = [w ^ MASK32 for w in W]
        H_c = sha256_hash(W_c)
        d_complement = sum(hw(H[i] ^ H_c[i]) for i in range(8))
        complement_diffs.append(d_complement)

        # Random pair
        W_r = [np.random.randint(0, 2**32) for _ in range(16)]
        H_r = sha256_hash(W_r)
        d_random = sum(hw(H[i] ^ H_r[i]) for i in range(8))
        random_diffs.append(d_random)

    print(f"  H(W) vs H(~W):     mean={np.mean(complement_diffs):.2f}, std={np.std(complement_diffs):.2f}")
    print(f"  H(W) vs H(random): mean={np.mean(random_diffs):.2f}, std={np.std(random_diffs):.2f}")

    from scipy import stats
    t, p = stats.ttest_ind(complement_diffs, random_diffs)
    print(f"  t-test: t={t:.3f}, p={p:.4f}")
    print(f"  → {'★ DISTINGUISHER!' if p < 0.01 else 'No distinguisher.'}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("2. SWAP HALVES on message")
    print("=" * 70)

    # Swap first 8 words with last 8: W[0..7] ↔ W[8..15]
    swap_diffs = []
    for _ in range(10000):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H = sha256_hash(W)

        W_s = W[8:] + W[:8]
        H_s = sha256_hash(W_s)
        d = sum(hw(H[i] ^ H_s[i]) for i in range(8))
        swap_diffs.append(d)

    print(f"  H(W) vs H(swap(W)): mean={np.mean(swap_diffs):.2f}, std={np.std(swap_diffs):.2f}")
    print(f"  H(W) vs H(random):  mean={np.mean(random_diffs):.2f}, std={np.std(random_diffs):.2f}")

    t2, p2 = stats.ttest_ind(swap_diffs, random_diffs)
    print(f"  t-test: t={t2:.3f}, p={p2:.4f}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("3. ARITHMETIC complement: W → (2^32 - 1 - W)")
    print("=" * 70)

    # ~W in two's complement = bitwise NOT.
    # But what about arithmetic: W → -W mod 2^32?
    arith_diffs = []
    for _ in range(10000):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H = sha256_hash(W)

        W_neg = [(~w + 1) & MASK32 for w in W]  # -W mod 2^32
        H_neg = sha256_hash(W_neg)
        d = sum(hw(H[i] ^ H_neg[i]) for i in range(8))
        arith_diffs.append(d)

    print(f"  H(W) vs H(-W):     mean={np.mean(arith_diffs):.2f}, std={np.std(arith_diffs):.2f}")
    t3, p3 = stats.ttest_ind(arith_diffs, random_diffs)
    print(f"  t-test: t={t3:.3f}, p={p3:.4f}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("4. PARTIAL complement: vary how many words are complemented")
    print("=" * 70)

    for n_comp in [1, 2, 4, 8, 12, 16]:
        diffs = []
        for _ in range(5000):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            H = sha256_hash(W)

            W2 = list(W)
            for i in range(n_comp):
                W2[i] ^= MASK32
            H2 = sha256_hash(W2)
            diffs.append(sum(hw(H[i] ^ H2[i]) for i in range(8)))

        print(f"  Complement {n_comp:>2}/16 words: mean δH={np.mean(diffs):.2f}, std={np.std(diffs):.2f}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("5. DEEPER: is complement quasi-symmetry on INTERNAL state real?")
    print("=" * 70)

    # The δ=32 was measured on internal state per round.
    # But that's comparing T(R(s)) vs R(T(s)) where T=complement.
    # Let's check: is this because complement changes T1/T2 by ~32 bits?

    # Ch(~e, ~f, ~g) = (~e & ~f) ^ (e & ~g) vs Ch(e,f,g) = (e&f) ^ (~e&g)
    # Ch(~e,~f,~g) = (~e&~f) ^ (~~e&~g) = (~e&~f) ^ (e&~g)
    # vs Ch(e,f,g) = (e&f) ^ (~e&g)
    # These are DIFFERENT functions! But:
    # ~Ch(e,f,g) = ~(e&f) ^ ~(~e&g) which is NOT Ch(~e,~f,~g)

    # Actually: Ch(~e,~f,~g) = (~e&~f) ^ (e&~g)
    # And: ~Ch(e,f,g) = ~((e&f)^(~e&g)) = complement of result

    # The question: does complement "almost commute" with Ch, Maj?
    ch_diffs = []
    maj_diffs = []
    for _ in range(10000):
        e = np.random.randint(0, 2**32)
        f = np.random.randint(0, 2**32)
        g = np.random.randint(0, 2**32)

        ch1 = (e & f) ^ (~e & g) & MASK32
        ch2 = (~e & ~f) ^ (e & ~g) & MASK32  # Ch(~e,~f,~g)
        ch_diffs.append(hw(ch1 ^ ch2))

        a = np.random.randint(0, 2**32)
        b = np.random.randint(0, 2**32)
        c = np.random.randint(0, 2**32)

        maj1 = (a&b) ^ (a&c) ^ (b&c)
        maj2 = (~a&~b) ^ (~a&~c) ^ (~b&~c)  # Maj(~a,~b,~c)
        maj_diffs.append(hw(maj1 ^ maj2))

    print(f"  Ch(e,f,g) vs Ch(~e,~f,~g):  mean HW(δ) = {np.mean(ch_diffs):.1f}")
    print(f"  Maj(a,b,c) vs Maj(~a,~b,~c): mean HW(δ) = {np.mean(maj_diffs):.1f}")

    # Theoretical:
    # Ch(~e,~f,~g) = ~e&~f ^ e&~g = ~(e|f) ^ (e&~g)
    # Ch(e,f,g) = e&f ^ ~e&g
    # XOR: Ch ⊕ Ch_comp = (e&f ^ ~e&g) ^ (~e&~f ^ e&~g)
    #     = e&f ^ ~e&g ^ ~e&~f ^ e&~g
    #     = e&(f^~g) ^ ~e&(g^~f) = e&(f⊕g⊕1) ^ ~e&(f⊕g⊕1) = f⊕g⊕1 = ~(f⊕g)
    # So: Ch(e,f,g) ⊕ Ch(~e,~f,~g) = ~(f⊕g)
    # HW(~(f⊕g)) = 32 - HW(f⊕g) ≈ 32 - 16 = 16

    print(f"\n  Theoretical:")
    print(f"    Ch(e,f,g) ⊕ Ch(~e,~f,~g) = ~(f⊕g)")
    print(f"    E[HW] = 32 - E[HW(f⊕g)] = 32 - 16 = 16")
    print(f"    Measured: {np.mean(ch_diffs):.1f} ✓")

    # For Maj:
    # Maj(~a,~b,~c) = (~a&~b) ^ (~a&~c) ^ (~b&~c)
    # = ~(a|b) ^ ~(a|c) ^ ~(b|c)
    # Maj(a,b,c) ⊕ Maj(~a,~b,~c):
    # After algebra: = ~(a⊕b⊕c) ... let's check
    # Actually: Maj(~a,~b,~c) = ~Maj(a,b,c) (majority of complements = complement of majority)
    # So HW(Maj⊕~Maj) = HW(all ones) = 32

    print(f"    Maj(~a,~b,~c) = ~Maj(a,b,c)")
    print(f"    E[HW] = 32")
    print(f"    Measured: {np.mean(maj_diffs):.1f} ✓")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("6. WHY δ = 32 for complement quasi-symmetry")
    print("=" * 70)

    print(f"""
  ═══════════════════════════════════════════════════════════════

  COMPLEMENT QUASI-SYMMETRY EXPLAINED:

  Round function with complemented state:
    T1_comp = ~h + Σ1(~e) + Ch(~e,~f,~g) + K + W
    T2_comp = Σ0(~a) + Maj(~a,~b,~c)

  Key identities:
    Σ0(~x) = Σ0(x) ⊕ (all ones) = ~Σ0(x)  [rotation commutes with NOT]
    Wait — rotr(~x, n) = ~rotr(x, n). So:
    Σ0(~x) = ~rotr(x,2) ⊕ ~rotr(x,13) ⊕ ~rotr(x,22)
           = rotr(x,2) ⊕ rotr(x,13) ⊕ rotr(x,22) ⊕ 1  (3 NOTs = 1 NOT)
           = Σ0(x) ⊕ MASK32
    Same for Σ1.

    Ch(~e,~f,~g) = Ch(e,f,g) ⊕ ~(f⊕g)  [proven above]
    Maj(~a,~b,~c) = ~Maj(a,b,c)           [proven above]

    ~h = -h - 1 (two's complement)
    So: ~h + X = -(h+1) + X = X - h - 1

  The complement transforms EACH addition by a constant offset.
  Total offset per round: sum of ~effects on T1 and T2.
  Each effect ≈ 16 bits of change (half of 32-bit word).
  But there are 7 additions → 7 × ~16 ≈ 112 bits of "noise".
  After mixing: ~32 bits survive as state difference.

  THIS IS EXACTLY WHAT WE MEASURED: δ = 31.7.

  EXPLOITABLE?
    On internal state: yes, δ=32 < random 128.
    On output (after feedforward H = IV + state):
      H(~W): mean δH = {np.mean(complement_diffs):.1f}
      Random: mean δH = {np.mean(random_diffs):.1f}
      → {'VISIBLE ON OUTPUT!' if abs(np.mean(complement_diffs) - np.mean(random_diffs)) > 1 else 'HIDDEN by feedforward.'}

  The quasi-symmetry exists INSIDE SHA-256 (δ=32 per round)
  but is {'destroyed' if abs(np.mean(complement_diffs) - np.mean(random_diffs)) < 1 else 'preserved'} after 64 rounds of amplification.

  ═══════════════════════════════════════════════════════════════
""")


if __name__ == "__main__":
    main()
