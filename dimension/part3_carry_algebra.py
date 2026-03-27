"""
ЧАСТЬ 3: Carry algebra — алгебра переносов.

SHA-256 round: 4 сложения mod 2^32 для T1, плюс T2.
Каждое сложение: c[i+1] = (a[i] & b[i]) | (c[i] & (a[i] ^ b[i]))
  — carry chain глубины 32.

За 1 раунд: ~128 бит carry.
За 64 раунда: ~8192 бит carry.
Стандартная криптография работает в XOR (GF2): ignores carry.
Мы работаем в Z/2^32: ВИДИМ carry.

НОВАЯ ОПЕРАЦИЯ: carry_trace
  Для двух сообщений M1, M2: отследить КАЖДЫЙ carry бит
  через все 64 раунда. Это 8192 бит дополнительной информации,
  которую XOR-анализ НЕ ВИДИТ.

Вопрос: содержат ли carry биты СТРУКТУРУ, которой нет в XOR?
"""

import numpy as np
import struct, hashlib

MASK32 = 0xFFFFFFFF

def hw(x): return bin(x).count('1')
def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)

K = [
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
]
IV = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]


def add32_with_carry(x, y):
    """ADD mod 2^32 + return carry bits."""
    result = (x + y) & MASK32
    # Carry at each bit position
    carry = 0
    carry_bits = 0
    for bit in range(32):
        a = (x >> bit) & 1
        b = (y >> bit) & 1
        s = a ^ b ^ carry
        carry = (a & b) | (carry & (a ^ b))
        carry_bits |= (carry << bit)
    return result, carry_bits


def sha256_with_carries(W16):
    """SHA-256 that also returns ALL carry bits from every addition."""
    W = list(W16)
    schedule_carries = []
    for r in range(16, 64):
        # W[r] = σ1(W[r-2]) + W[r-7] + σ0(W[r-15]) + W[r-16]
        s1 = sigma1(W[r-2])
        s0 = sigma0(W[r-15])
        t1, c1 = add32_with_carry(s1, W[r-7])
        t2, c2 = add32_with_carry(t1, s0)
        t3, c3 = add32_with_carry(t2, W[r-16])
        W.append(t3)
        schedule_carries.append((c1, c2, c3))

    a, b, c, d, e, f, g, h = IV
    round_carries = []

    for r in range(64):
        # T1 = h + Σ1(e) + Ch(e,f,g) + K[r] + W[r]
        sig1e = Sigma1(e)
        ch_efg = Ch(e, f, g)

        t1, c1 = add32_with_carry(h, sig1e)
        t2, c2 = add32_with_carry(t1, ch_efg)
        t3, c3 = add32_with_carry(t2, K[r])
        T1, c4 = add32_with_carry(t3, W[r])

        # T2 = Σ0(a) + Maj(a,b,c)
        sig0a = Sigma0(a)
        maj_abc = Maj(a, b, c)
        T2, c5 = add32_with_carry(sig0a, maj_abc)

        # a_new = T1 + T2
        a_new, c6 = add32_with_carry(T1, T2)
        # e_new = d + T1
        e_new, c7 = add32_with_carry(d, T1)

        round_carries.append((c1, c2, c3, c4, c5, c6, c7))

        h, g, f, e = g, f, e, e_new
        d, c, b, a = c, b, a, a_new

    # Final addition
    final_carries = []
    H = []
    for i, (iv, reg) in enumerate(zip(IV, [a,b,c,d,e,f,g,h])):
        result, carry = add32_with_carry(iv, reg)
        H.append(result)
        final_carries.append(carry)

    return tuple(H), round_carries, schedule_carries, final_carries


def main():
    np.random.seed(42)

    print("=" * 70)
    print("ЧАСТЬ 3: Carry algebra")
    print("=" * 70)

    # ═══════════════════════════════
    print(f"\n{'=' * 70}")
    print("ТЕСТ A: Carry bit statistics")
    print("=" * 70)

    # Сколько carry бит возникает per round?
    total_carries_per_round = []
    for _ in range(1000):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H, round_carries, _, _ = sha256_with_carries(W)

        for r in range(64):
            total_c = sum(hw(c) for c in round_carries[r])
            total_carries_per_round.append(total_c)

    total_carries = np.array(total_carries_per_round).reshape(1000, 64)
    print(f"  Carry bits per round:")
    print(f"    Mean: {np.mean(total_carries):.1f}")
    print(f"    Per addition: {np.mean(total_carries)/7:.1f} (7 additions per round)")
    print(f"    Expected (random ADD): ~15.5 per addition × 7 = ~108.5")

    # Per-round profile
    round_means = np.mean(total_carries, axis=0)
    print(f"\n  Round profile (mean carry bits):")
    for r in [0, 4, 8, 16, 32, 48, 63]:
        print(f"    r={r:2d}: {round_means[r]:.1f}")

    # Is carry count CONSTANT or varies?
    print(f"    Variance across rounds: {np.std(round_means):.2f}")

    # ═══════════════════════════════
    print(f"\n{'=' * 70}")
    print("ТЕСТ B: Carry CORRELATION between rounds")
    print("=" * 70)

    # Key question: are carries in round r CORRELATED with carries in round r+1?
    # For random: no correlation expected.
    # For SHA-256: carry propagation MAY create correlation.

    carry_per_round = np.mean(total_carries, axis=0)

    # Autocorrelation of carry count across rounds (per-message)
    autocorrs = []
    for msg_idx in range(1000):
        carry_seq = total_carries[msg_idx]
        carry_seq = carry_seq - np.mean(carry_seq)
        if np.std(carry_seq) > 0:
            ac = np.sum(carry_seq[:-1] * carry_seq[1:]) / np.sum(carry_seq**2)
            autocorrs.append(ac)

    print(f"  Carry autocorrelation (r → r+1):")
    print(f"    Mean: {np.mean(autocorrs):.4f}")
    print(f"    Std:  {np.std(autocorrs):.4f}")
    significant = abs(np.mean(autocorrs)) > 2 * np.std(autocorrs) / np.sqrt(len(autocorrs))
    print(f"    → {'★ SIGNIFICANT!' if significant else 'Not significant.'}")

    # ═══════════════════════════════
    print(f"\n{'=' * 70}")
    print("ТЕСТ C: Carry vs Output — предсказывают ли carries выход?")
    print("=" * 70)

    # If we know ALL carry bits, does that help predict the output?
    # Trivially yes (carries + XOR = full computation).
    # But: do PARTIAL carries help?

    # Test: does carry count in round 0 correlate with HW(H[0])?
    round0_carries = []
    h0_hws = []
    for _ in range(10000):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H, rc, _, _ = sha256_with_carries(W)
        round0_carries.append(sum(hw(c) for c in rc[0]))
        h0_hws.append(hw(H[0]))

    corr = np.corrcoef(round0_carries, h0_hws)[0, 1]
    print(f"  corr(carries_round0, HW(H[0])): {corr:.6f}")

    # Does TOTAL carries correlate with output?
    total_msg_carries = []
    total_output_hws = []
    for _ in range(10000):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H, rc, _, _ = sha256_with_carries(W)
        tc = sum(sum(hw(c) for c in carries) for carries in rc)
        total_msg_carries.append(tc)
        total_output_hws.append(sum(hw(h) for h in H))

    corr2 = np.corrcoef(total_msg_carries, total_output_hws)[0, 1]
    print(f"  corr(total_carries, HW(H)): {corr2:.6f}")
    print(f"  → {'★ CARRIES PREDICT OUTPUT!' if abs(corr2) > 0.02 else 'No predictive power.'}")

    # ═══════════════════════════════
    print(f"\n{'=' * 70}")
    print("ТЕСТ D: Carry difference — δcarry между двумя сообщениями")
    print("=" * 70)

    # For two messages M1, M2 with δW = 1 bit:
    # What's the carry DIFFERENCE at each round?
    # XOR-анализ видит δstate. Carry-анализ видит δcarry.
    # Есть ли в δcarry информация, которой нет в δstate?

    carry_diffs = []
    state_diffs = []

    for _ in range(2000):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W1); W2[0] ^= 1

        H1, rc1, _, _ = sha256_with_carries(W1)
        H2, rc2, _, _ = sha256_with_carries(W2)

        # δcarry per round
        dc_per_round = []
        for r in range(64):
            dc = sum(hw(rc1[r][i] ^ rc2[r][i]) for i in range(7))
            dc_per_round.append(dc)
        carry_diffs.append(dc_per_round)

    carry_diffs = np.array(carry_diffs)
    print(f"  δcarry per round (δW = 1 bit in W[0]):")
    for r in [0, 1, 2, 4, 8, 16, 32, 63]:
        print(f"    r={r:2d}: mean={np.mean(carry_diffs[:,r]):.1f}")

    # ═══════════════════════════════
    print(f"\n{'=' * 70}")
    print("ТЕСТ E: Carry как ДОПОЛНИТЕЛЬНАЯ ИНФОРМАЦИЯ")
    print("=" * 70)

    # KEY EXPERIMENT:
    # 1. Compute two near-collision messages (HW(δH) small)
    # 2. Look at their carry differences
    # 3. Is δcarry ALSO small? Or independent?
    #
    # If δcarry is correlated with δH → carry gives FREE information

    # Generate pairs with varying δH
    pairs_by_dH = {range_val: [] for range_val in [(90,100), (100,110), (110,120), (120,130), (130,140)]}

    for _ in range(20000):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W1)
        w = np.random.randint(0, 16)
        b = np.random.randint(0, 32)
        W2[w] ^= (1 << b)

        H1, rc1, _, fc1 = sha256_with_carries(W1)
        H2, rc2, _, fc2 = sha256_with_carries(W2)

        dH = sum(hw(H1[i]^H2[i]) for i in range(8))

        # δ in final carries
        dc_final = sum(hw(fc1[i]^fc2[i]) for i in range(8))

        for (lo, hi) in pairs_by_dH:
            if lo <= dH < hi:
                pairs_by_dH[(lo,hi)].append(dc_final)
                break

    print(f"  δH range → mean δ(final carries):")
    for (lo, hi) in sorted(pairs_by_dH.keys()):
        vals = pairs_by_dH[(lo,hi)]
        if vals:
            print(f"    δH ∈ [{lo},{hi}): mean δcarry = {np.mean(vals):.1f} (n={len(vals)})")

    print(f"\n  If δcarry correlates with δH → carries give EXTRA information")
    print(f"  If δcarry ≈ constant → carries are independent (no extra info)")

    # ═══════════════════════════════
    print(f"\n{'=' * 70}")
    print("ИТОГ ЧАСТИ 3")
    print("=" * 70)


if __name__ == "__main__":
    main()
