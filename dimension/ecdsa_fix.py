"""
ФИКС DISTINGUISHER: трассировка сигнала от H до k.

Сигнал: SHA-256 H имеет carry[63] bias (AUC=0.976 из методички).
Вопрос: ГДЕ сигнал теряется?

Цепочка: H → int(H) → int(H) mod N → k → k·G → r
На каждом шаге измеряем сигнал.

Дополнительно: используем правильный score из методички:
  - H[6][b31,b30,b29,b28] = (0,0,1,0) с P=0.263 (vs random 0.0625)
  - adj_xor_deviation (z=8.8 на 50K)
  - Cross-word asymmetry
"""

import hashlib
import os
import numpy as np

N_CURVE = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141


def sha256_bytes(data):
    return hashlib.sha256(data).digest()


def hw(x):
    return bin(x).count('1')


def h_to_words(h_bytes):
    """32 bytes → 8 words (big-endian, как SHA-256 output)."""
    words = []
    for i in range(8):
        w = int.from_bytes(h_bytes[i*4:(i+1)*4], 'big')
        words.append(w)
    return words


def score_h6_pattern(h_bytes):
    """T_CH_INVARIANT pattern: H[6][b31,b30,b29,b28] = (0,0,1,0).
    P(SHA-256) = 0.263, P(random) = 0.0625. Lift = 4.2×."""
    words = h_to_words(h_bytes)
    h6 = words[6]
    b31 = (h6 >> 31) & 1
    b30 = (h6 >> 30) & 1
    b29 = (h6 >> 29) & 1
    b28 = (h6 >> 28) & 1
    return (b31 == 0 and b30 == 0 and b29 == 1 and b28 == 0)


def score_adj_xor(h_bytes):
    """Adjacent bit XOR uniformity score."""
    x = int.from_bytes(h_bytes, 'big')
    adj_xor_count = 0
    for b in range(255):
        b0 = (x >> b) & 1
        b1 = (x >> (b + 1)) & 1
        adj_xor_count += b0 ^ b1
    return adj_xor_count / 255  # close to 0.5 for both, but variance differs


def score_word_asymmetry(h_bytes):
    """Cross-word HW asymmetry: |std(H[0..3]) - std(H[4..7])|."""
    words = h_to_words(h_bytes)
    hws = [hw(w) for w in words]
    a_path = np.std(hws[:4])
    e_path = np.std(hws[4:])
    return abs(a_path - e_path)


def score_h7_upper(h_bytes):
    """H[7] upper bits bias."""
    words = h_to_words(h_bytes)
    h7 = words[7]
    # T_H7_BIAS: E[H[7]/2^32] = 0.529 instead of 0.500
    return h7 / 2**32


def experiment():
    np.random.seed(42)
    privkey = os.urandom(32)

    print("=" * 70)
    print("ФИКС: ТРАССИРОВКА СИГНАЛА ОТ H ДО k")
    print("=" * 70)

    # =================================================================
    print("\n" + "=" * 70)
    print("1. СИГНАЛ В RAW HASH H (без mod N)")
    print("=" * 70)

    N_SAMPLES = 50000

    # SHA-256 hashes с фиксированным prefix (как уязвимый nonce)
    sha_hashes = []
    for i in range(N_SAMPLES):
        msg = f"tx_{i}".encode()
        h = sha256_bytes(privkey + msg)
        sha_hashes.append(h)

    # True random bytes
    rand_hashes = [os.urandom(32) for _ in range(N_SAMPLES)]

    # H[6] pattern test
    sha_h6 = sum(score_h6_pattern(h) for h in sha_hashes) / N_SAMPLES
    rand_h6 = sum(score_h6_pattern(h) for h in rand_hashes) / N_SAMPLES
    theo_h6 = 1/16  # 4 bits random

    print(f"\n  H[6] pattern (0,0,1,0):")
    print(f"    SHA-256:  {sha_h6:.4f}  (методичка: 0.263)")
    print(f"    Random:   {rand_h6:.4f}  (теория: {theo_h6:.4f})")
    print(f"    Lift:     {sha_h6/max(rand_h6,0.001):.2f}×")

    # Adj XOR
    sha_adj = np.mean([score_adj_xor(h) for h in sha_hashes[:5000]])
    rand_adj = np.mean([score_adj_xor(h) for h in rand_hashes[:5000]])
    print(f"\n  Adj XOR mean:")
    print(f"    SHA-256:  {sha_adj:.6f}")
    print(f"    Random:   {rand_adj:.6f}")

    # H[7] upper bias
    sha_h7 = np.mean([score_h7_upper(h) for h in sha_hashes])
    rand_h7 = np.mean([score_h7_upper(h) for h in rand_hashes])
    print(f"\n  E[H[7]/2^32]:")
    print(f"    SHA-256:  {sha_h7:.6f}  (методичка: 0.529)")
    print(f"    Random:   {rand_h7:.6f}  (теория: 0.500)")

    # Per-word HW
    print(f"\n  Mean HW по словам:")
    print(f"  {'Word':>6} {'SHA-256':>10} {'Random':>10} {'Δ':>8}")
    word_names = ['H[0]', 'H[1]', 'H[2]', 'H[3]', 'H[4]', 'H[5]', 'H[6]', 'H[7]']
    for w in range(8):
        sha_w = np.mean([hw(h_to_words(h)[w]) for h in sha_hashes[:10000]])
        rand_w = np.mean([hw(h_to_words(h)[w]) for h in rand_hashes[:10000]])
        delta = sha_w - rand_w
        print(f"  {word_names[w]:>6} {sha_w:10.4f} {rand_w:10.4f} {delta:+7.4f}")

    # =================================================================
    print("\n" + "=" * 70)
    print("2. ТРАССИРОВКА: где теряется сигнал?")
    print("=" * 70)

    # Step 0: Raw H
    # Step 1: int(H) — 256-бит integer
    # Step 2: int(H) mod N — mod secp256k1 order
    # Step 3: int(H) mod N → bytes — обратно в bytes

    sha_raw = sha_hashes[:10000]
    rand_raw = rand_hashes[:10000]

    # На каждом шаге: H[6] pattern survival
    # Step 0: raw hash
    s0_sha = sum(score_h6_pattern(h) for h in sha_raw) / len(sha_raw)
    s0_rand = sum(score_h6_pattern(h) for h in rand_raw) / len(rand_raw)

    # Step 1: int → mod N → back to bytes
    def mod_n_bytes(h):
        val = int.from_bytes(h, 'big')
        val_mod = val % N_CURVE
        return val_mod.to_bytes(32, 'big')

    sha_modn = [mod_n_bytes(h) for h in sha_raw]
    rand_modn = [mod_n_bytes(h) for h in rand_raw]

    s1_sha = sum(score_h6_pattern(h) for h in sha_modn) / len(sha_modn)
    s1_rand = sum(score_h6_pattern(h) for h in rand_modn) / len(rand_modn)

    # Step 2: truncate to different bit widths
    print(f"\n  H[6] pattern survival через трансформации:")
    print(f"  {'Шаг':<30} {'SHA':>8} {'Random':>8} {'Lift':>8} {'Signal':>8}")
    print(f"  {'-'*30} {'-'*8} {'-'*8} {'-'*8} {'-'*8}")

    lift0 = s0_sha / max(s0_rand, 0.001)
    lift1 = s1_sha / max(s1_rand, 0.001)

    print(f"  {'Raw H':30s} {s0_sha:8.4f} {s0_rand:8.4f} {lift0:7.2f}× {'✓' if lift0 > 1.5 else '✗'}")
    print(f"  {'int(H) mod N → bytes':30s} {s1_sha:8.4f} {s1_rand:8.4f} {lift1:7.2f}× {'✓' if lift1 > 1.5 else '✗'}")

    # H[7] bias survival
    h7_raw_sha = np.mean([score_h7_upper(h) for h in sha_raw])
    h7_raw_rand = np.mean([score_h7_upper(h) for h in rand_raw])
    h7_modn_sha = np.mean([score_h7_upper(h) for h in sha_modn])
    h7_modn_rand = np.mean([score_h7_upper(h) for h in rand_modn])

    print(f"\n  E[H[7]/2^32] survival:")
    print(f"  {'Raw H':30s} SHA={h7_raw_sha:.6f} rand={h7_raw_rand:.6f} Δ={h7_raw_sha-h7_raw_rand:+.6f}")
    print(f"  {'mod N → bytes':30s} SHA={h7_modn_sha:.6f} rand={h7_modn_rand:.6f} Δ={h7_modn_sha-h7_modn_rand:+.6f}")

    # Per-bit analysis: which bits survive mod N?
    print(f"\n  Per-bit bias survival (|P(bit=1) - 0.5|):")
    print(f"  {'Бит':>5} {'Raw SHA':>10} {'Raw Rand':>10} {'ModN SHA':>10} {'ModN Rand':>10}")

    for bit in [0, 7, 15, 16, 24, 28, 29, 30, 31,   # word-level bits
                224, 225, 226, 227, 228, 229, 230, 231,  # H[7] bits
                192, 193, 194, 195, 196, 197, 198, 199]:  # H[6] bits

        raw_sha_freq = np.mean([(int.from_bytes(h,'big') >> bit) & 1 for h in sha_raw])
        raw_rand_freq = np.mean([(int.from_bytes(h,'big') >> bit) & 1 for h in rand_raw])
        modn_sha_freq = np.mean([(int.from_bytes(h,'big') >> bit) & 1 for h in sha_modn])
        modn_rand_freq = np.mean([(int.from_bytes(h,'big') >> bit) & 1 for h in rand_modn])

        marker = ""
        if abs(raw_sha_freq - 0.5) > 0.01:
            marker = " ★raw"
        if abs(modn_sha_freq - 0.5) > 0.01:
            marker += " ★modn"

        if bit in [224,225,226,227,228,229,230,231, 192,193,194,195,196,197,198,199, 0, 31]:
            print(f"  {bit:5d} {abs(raw_sha_freq-0.5):10.4f} {abs(raw_rand_freq-0.5):10.4f} "
                  f"{abs(modn_sha_freq-0.5):10.4f} {abs(modn_rand_freq-0.5):10.4f}{marker}")

    # =================================================================
    print("\n" + "=" * 70)
    print("3. АЛЬТЕРНАТИВНЫЙ ПОДХОД: фиксированный prefix detection")
    print("=" * 70)

    # k_i = SHA256(FIXED_PREFIX || msg_i)
    # Все k_i имеют одинаковые первые 8 слов W[0..7].
    # Это создаёт КОРРЕЛЯЦИЮ между k_i, невидимую в random.

    # Test: pairwise correlation between hash outputs
    sha_pairs_hw = []
    rand_pairs_hw = []

    for i in range(0, 2000, 2):
        h1_sha = sha_hashes[i]
        h2_sha = sha_hashes[i+1]
        delta_sha = bytes(a ^ b for a, b in zip(h1_sha, h2_sha))
        sha_pairs_hw.append(sum(hw(b) for b in delta_sha))

        h1_rand = rand_hashes[i]
        h2_rand = rand_hashes[i+1]
        delta_rand = bytes(a ^ b for a, b in zip(h1_rand, h2_rand))
        rand_pairs_hw.append(sum(hw(b) for b in delta_rand))

    print(f"\n  Pairwise HW(H_i ⊕ H_j) — корреляция между выходами:")
    print(f"    SHA-256 (fixed prefix): mean={np.mean(sha_pairs_hw):.2f} ± {np.std(sha_pairs_hw):.2f}")
    print(f"    Random:                 mean={np.mean(rand_pairs_hw):.2f} ± {np.std(rand_pairs_hw):.2f}")
    print(f"    Δ: {np.mean(sha_pairs_hw) - np.mean(rand_pairs_hw):+.2f}")

    # Per-word pairwise
    print(f"\n  Per-word pairwise δ:")
    for w in range(8):
        sha_w_pairs = []
        rand_w_pairs = []
        for i in range(0, 2000, 2):
            w1s = h_to_words(sha_hashes[i])[w]
            w2s = h_to_words(sha_hashes[i+1])[w]
            sha_w_pairs.append(hw(w1s ^ w2s))

            w1r = h_to_words(rand_hashes[i])[w]
            w2r = h_to_words(rand_hashes[i+1])[w]
            rand_w_pairs.append(hw(w1r ^ w2r))

        delta = np.mean(sha_w_pairs) - np.mean(rand_w_pairs)
        marker = " ★" if abs(delta) > 0.1 else ""
        print(f"    {word_names[w]}: SHA={np.mean(sha_w_pairs):.3f} rand={np.mean(rand_w_pairs):.3f} Δ={delta:+.3f}{marker}")

    # =================================================================
    print("\n" + "=" * 70)
    print("4. CONDITIONAL TEST: фильтруем carry[63]-like хеши")
    print("=" * 70)

    # Из методички: carry[63]=0 → H[6] pattern (0,0,1,0) с P=0.263
    # Используем паттерн как ФИЛЬТР: берём только хеши где H[6] pattern = True
    sha_filtered = [h for h in sha_hashes if score_h6_pattern(h)]
    rand_filtered = [h for h in rand_hashes if score_h6_pattern(h)]

    print(f"\n  Фильтр H[6] pattern (0,0,1,0):")
    print(f"    SHA-256: {len(sha_filtered)}/{N_SAMPLES} = {len(sha_filtered)/N_SAMPLES*100:.2f}%")
    print(f"    Random:  {len(rand_filtered)}/{N_SAMPLES} = {len(rand_filtered)/N_SAMPLES*100:.2f}%")
    print(f"    Lift:    {(len(sha_filtered)/N_SAMPLES) / max(len(rand_filtered)/N_SAMPLES, 0.001):.3f}×")

    if len(sha_filtered) > 100 and len(rand_filtered) > 100:
        # Внутри отфильтрованных: есть ли дополнительный bias?
        sha_h7_filt = np.mean([score_h7_upper(h) for h in sha_filtered[:1000]])
        rand_h7_filt = np.mean([score_h7_upper(h) for h in rand_filtered[:1000]])
        print(f"\n  E[H[7]/2^32] внутри фильтра:")
        print(f"    SHA-256: {sha_h7_filt:.6f}")
        print(f"    Random:  {rand_h7_filt:.6f}")
        print(f"    Δ: {sha_h7_filt - rand_h7_filt:+.6f}")

        # H[7] bits 29,30 inside filter
        sha_b29 = np.mean([(h_to_words(h)[7] >> 29) & 1 for h in sha_filtered[:1000]])
        sha_b30 = np.mean([(h_to_words(h)[7] >> 30) & 1 for h in sha_filtered[:1000]])
        rand_b29 = np.mean([(h_to_words(h)[7] >> 29) & 1 for h in rand_filtered[:1000]])
        rand_b30 = np.mean([(h_to_words(h)[7] >> 30) & 1 for h in rand_filtered[:1000]])

        print(f"\n  H[7] bits внутри фильтра:")
        print(f"    P(b29=1): SHA={sha_b29:.4f} rand={rand_b29:.4f} Δ={sha_b29-rand_b29:+.4f}")
        print(f"    P(b30=1): SHA={sha_b30:.4f} rand={rand_b30:.4f} Δ={sha_b30-rand_b30:+.4f}")

    # =================================================================
    print("\n" + "=" * 70)
    print("5. ВЕРДИКТ: ГДЕ СИГНАЛ И КАК ЕГО ИСПОЛЬЗОВАТЬ")
    print("=" * 70)


if __name__ == "__main__":
    experiment()
