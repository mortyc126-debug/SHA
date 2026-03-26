"""
ECDSA Conditional Distinguisher — полный тест.

Что работает (из трассировки):
  1. H[6] pattern как фильтр для carry[63]≈0
  2. Внутри фильтра: H[7][b29] bias Δ=+0.043
  3. e-path asymmetry: H[4,5] pairwise Δ=+0.23

Тест:
  - 1000 уязвимых кошельков (k=SHA256(priv||msg))
  - 1000 безопасных кошельков (k=random)
  - Для каждого: N подписей
  - Conditional score → классификация
  - Accuracy по количеству подписей
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
    return [int.from_bytes(h_bytes[i*4:(i+1)*4], 'big') for i in range(8)]


def h6_filter(h_bytes):
    """True если H[6] upper nibble = 0b0010xxxx (pattern 0,0,1,0)."""
    h6 = h_to_words(h_bytes)[6]
    return ((h6 >> 28) & 0xF) == 0x2


def conditional_score(hash_list):
    """
    Conditional distinguisher:
      1. Фильтруем хеши через H[6] pattern
      2. Внутри фильтра: считаем H[7][b29,b30] bias
      3. Считаем e-path pairwise asymmetry
    Возвращаем composite score.
    """
    n = len(hash_list)
    if n == 0:
        return 0.5, {}

    # === Feature 1: H[6] filter rate ===
    filtered = [h for h in hash_list if h6_filter(h)]
    filter_rate = len(filtered) / n

    # === Feature 2: H[7] bias inside filter ===
    if len(filtered) >= 3:
        b29_freq = np.mean([(h_to_words(h)[7] >> 29) & 1 for h in filtered])
        b30_freq = np.mean([(h_to_words(h)[7] >> 30) & 1 for h in filtered])
        h7_bias = (b29_freq - 0.5) + (b30_freq - 0.5)
    else:
        h7_bias = 0.0
        b29_freq = 0.5
        b30_freq = 0.5

    # === Feature 3: E[H[7]/2^32] inside filter ===
    if len(filtered) >= 3:
        h7_mean = np.mean([h_to_words(h)[7] / 2**32 for h in filtered])
    else:
        h7_mean = 0.5

    # === Feature 4: e-path pairwise asymmetry ===
    if n >= 4:
        e_path_pairs = []
        a_path_pairs = []
        for i in range(0, min(n, 200) - 1, 2):
            w1 = h_to_words(hash_list[i])
            w2 = h_to_words(hash_list[i+1])
            e_path_pairs.append(np.mean([hw(w1[k] ^ w2[k]) for k in range(4, 8)]))
            a_path_pairs.append(np.mean([hw(w1[k] ^ w2[k]) for k in range(0, 4)]))
        asymmetry = np.mean(e_path_pairs) - np.mean(a_path_pairs)
    else:
        asymmetry = 0.0

    # === Feature 5: Unconditional H[7] mean ===
    h7_uncond = np.mean([h_to_words(h)[7] / 2**32 for h in hash_list])

    details = {
        'filter_rate': filter_rate,
        'n_filtered': len(filtered),
        'b29': b29_freq,
        'b30': b30_freq,
        'h7_bias': h7_bias,
        'h7_mean_cond': h7_mean,
        'h7_mean_uncond': h7_uncond,
        'asymmetry': asymmetry,
    }

    # Composite score: weighted sum
    # Каждая фича: (value - random_baseline) * weight
    score = 0.0
    score += (filter_rate - 0.0625) * 10     # H[6] filter excess
    score += h7_bias * 5                      # conditional H[7] bias
    score += (h7_mean - 0.5) * 3             # conditional H[7] mean
    score += asymmetry * 2                    # e-path asymmetry
    score += (h7_uncond - 0.5) * 3           # unconditional H[7]

    return score, details


def generate_wallet_signatures(n_sigs, mode='vulnerable'):
    """Генерируем n_sigs хешей (прокси для nonce)."""
    privkey = os.urandom(32)
    hashes = []
    for i in range(n_sigs):
        msg = f"tx_{i}_{os.urandom(4).hex()}".encode()
        if mode == 'vulnerable':
            h = sha256_bytes(privkey + msg)
        else:
            h = os.urandom(32)
        hashes.append(h)
    return hashes


def experiment():
    np.random.seed(42)

    print("=" * 70)
    print("CONDITIONAL DISTINGUISHER — ПОЛНЫЙ ТЕСТ")
    print("=" * 70)

    # =================================================================
    print("\n" + "=" * 70)
    print("1. CALIBRATION: score distribution для SHA vs Random")
    print("=" * 70)

    for n_sigs in [50, 200, 1000, 5000]:
        sha_scores = []
        rand_scores = []

        n_wallets = 200
        for _ in range(n_wallets):
            h_sha = generate_wallet_signatures(n_sigs, 'vulnerable')
            h_rand = generate_wallet_signatures(n_sigs, 'random')

            s_sha, _ = conditional_score(h_sha)
            s_rand, _ = conditional_score(h_rand)
            sha_scores.append(s_sha)
            rand_scores.append(s_rand)

        mean_sha = np.mean(sha_scores)
        mean_rand = np.mean(rand_scores)
        std_sha = np.std(sha_scores)
        std_rand = np.std(rand_scores)

        # Separability: (mean_sha - mean_rand) / sqrt(std²_sha + std²_rand)
        sep = abs(mean_sha - mean_rand) / max(np.sqrt(std_sha**2 + std_rand**2), 0.001)

        print(f"\n  N={n_sigs:5d}:")
        print(f"    SHA score:  {mean_sha:+.4f} ± {std_sha:.4f}")
        print(f"    Rand score: {mean_rand:+.4f} ± {std_rand:.4f}")
        print(f"    Separability: {sep:.3f}  {'← SEPARABLE!' if sep > 1.0 else ''}")

    # =================================================================
    print("\n" + "=" * 70)
    print("2. CLASSIFICATION: accuracy по числу подписей")
    print("=" * 70)

    print(f"\n  {'N_sigs':>7} {'Accuracy':>10} {'TPR':>8} {'TNR':>8} {'Sep':>8} {'Status':>10}")
    print(f"  {'-'*7} {'-'*10} {'-'*8} {'-'*8} {'-'*8} {'-'*10}")

    for n_sigs in [10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000]:
        n_wallets = 200  # 100 SHA + 100 random

        # Сначала определяем threshold на calibration set
        cal_sha = [conditional_score(generate_wallet_signatures(n_sigs, 'vulnerable'))[0]
                   for _ in range(50)]
        cal_rand = [conditional_score(generate_wallet_signatures(n_sigs, 'random'))[0]
                    for _ in range(50)]
        threshold = (np.mean(cal_sha) + np.mean(cal_rand)) / 2

        # Test set
        tp, tn, fp, fn = 0, 0, 0, 0
        for trial in range(n_wallets):
            is_sha = trial < n_wallets // 2

            if is_sha:
                hashes = generate_wallet_signatures(n_sigs, 'vulnerable')
            else:
                hashes = generate_wallet_signatures(n_sigs, 'random')

            score, _ = conditional_score(hashes)
            predicted_sha = score > threshold

            if is_sha and predicted_sha:
                tp += 1
            elif is_sha and not predicted_sha:
                fn += 1
            elif not is_sha and predicted_sha:
                fp += 1
            else:
                tn += 1

        accuracy = (tp + tn) / n_wallets * 100
        tpr = tp / max(tp + fn, 1) * 100
        tnr = tn / max(tn + fp, 1) * 100

        sep = abs(np.mean(cal_sha) - np.mean(cal_rand)) / max(
            np.sqrt(np.std(cal_sha)**2 + np.std(cal_rand)**2), 0.001)

        status = ""
        if accuracy > 70: status = "★"
        if accuracy > 80: status = "★★"
        if accuracy > 90: status = "★★★"

        print(f"  {n_sigs:7d} {accuracy:9.1f}% {tpr:7.1f}% {tnr:7.1f}% {sep:7.3f}  {status}")

    # =================================================================
    print("\n" + "=" * 70)
    print("3. ДЕТАЛЬНЫЙ АНАЛИЗ: какая фича работает?")
    print("=" * 70)

    n_sigs = 5000
    n_test = 100

    feature_scores = {k: {'sha': [], 'rand': []} for k in
                      ['filter_rate', 'h7_bias', 'h7_mean_cond', 'h7_mean_uncond', 'asymmetry']}

    for trial in range(n_test):
        h_sha = generate_wallet_signatures(n_sigs, 'vulnerable')
        h_rand = generate_wallet_signatures(n_sigs, 'random')

        _, d_sha = conditional_score(h_sha)
        _, d_rand = conditional_score(h_rand)

        for k in feature_scores:
            feature_scores[k]['sha'].append(d_sha.get(k, 0))
            feature_scores[k]['rand'].append(d_rand.get(k, 0))

    print(f"\n  Per-feature analysis (N={n_sigs}, {n_test} wallets each):")
    print(f"\n  {'Feature':<20} {'SHA':>10} {'Random':>10} {'Δ':>10} {'Sep':>8} {'Works?':>8}")
    print(f"  {'-'*20} {'-'*10} {'-'*10} {'-'*10} {'-'*8} {'-'*8}")

    for k in feature_scores:
        sha_vals = feature_scores[k]['sha']
        rand_vals = feature_scores[k]['rand']
        m_sha = np.mean(sha_vals)
        m_rand = np.mean(rand_vals)
        delta = m_sha - m_rand
        sep = abs(delta) / max(np.sqrt(np.std(sha_vals)**2 + np.std(rand_vals)**2), 0.001)
        works = "✓" if sep > 0.5 else "✗"

        print(f"  {k:<20} {m_sha:10.5f} {m_rand:10.5f} {delta:+9.5f} {sep:7.3f}    {works}")

    # =================================================================
    print("\n" + "=" * 70)
    print("4. ВЕРДИКТ")
    print("=" * 70)


if __name__ == "__main__":
    experiment()
