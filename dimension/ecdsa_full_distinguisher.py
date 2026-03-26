"""
Полный distinguisher на ECDSA nonce.

Проблема прошлого теста: использовали только 2 бита (carry[63] proxy).
Bias = 0.143 на 2 бита → теряется в шуме при mod N.

Полный подход:
  1. Intra-word bit correlation (carry chains → adjacent bit dependency)
  2. Cross-word asymmetry (a-path vs e-path разное число ADD)
  3. Bit-pair XOR statistics (SHA-256 ≠ random на уровне бит-пар)
  4. Hamming weight distribution (SHA-256 output HW не идеально нормален)
  5. Byte-level entropy

Собираем ВСЕ фичи → накапливаем статистику → классифицируем.
"""

import hashlib
import os
import numpy as np

# secp256k1 order
N_CURVE = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141


def sha256_bytes(data):
    return hashlib.sha256(data).digest()


def generate_nonce_vulnerable(privkey_bytes, msg):
    """k = SHA256(privkey || msg) mod N — УЯЗВИМЫЙ."""
    h = sha256_bytes(privkey_bytes + msg)
    k = int.from_bytes(h, 'big') % (N_CURVE - 1) + 1
    return k, h  # возвращаем и hash для анализа


def generate_nonce_random(msg):
    """k = random — БЕЗОПАСНЫЙ."""
    k = int.from_bytes(os.urandom(32), 'big') % (N_CURVE - 1) + 1
    h = k.to_bytes(32, 'big')
    return k, h


def extract_features(hash_bytes_list):
    """
    Извлекаем фичи из набора 32-байтных значений.
    Каждая фича — скаляр, накопленный по всем сэмплам.
    """
    n = len(hash_bytes_list)
    if n == 0:
        return {}

    features = {}

    # Преобразуем в массив бит
    all_ints = [int.from_bytes(h, 'big') for h in hash_bytes_list]

    # === 1. HW distribution ===
    hws = [bin(x).count('1') for x in all_ints]
    features['hw_mean'] = np.mean(hws)
    features['hw_std'] = np.std(hws)
    features['hw_skew'] = float(np.mean(((np.array(hws) - 128) / 8) ** 3))

    # === 2. Per-bit frequency (отклонение от 0.5) ===
    bit_freqs = np.zeros(256)
    for x in all_ints:
        for b in range(256):
            if (x >> b) & 1:
                bit_freqs[b] += 1
    bit_freqs /= n

    # Суммарное отклонение от 0.5
    features['bit_bias_sum'] = float(np.sum(np.abs(bit_freqs - 0.5)))
    features['bit_bias_max'] = float(np.max(np.abs(bit_freqs - 0.5)))

    # === 3. Adjacent bit XOR (carry chain signature) ===
    # SHA-256: adjacent bits коррелированы из-за carry
    adj_xor = np.zeros(255)
    for x in all_ints:
        for b in range(255):
            b0 = (x >> b) & 1
            b1 = (x >> (b + 1)) & 1
            adj_xor[b] += b0 ^ b1
    adj_xor /= n

    # Для random: P(adj XOR = 1) = 0.5 точно
    features['adj_xor_mean'] = float(np.mean(adj_xor))
    features['adj_xor_std'] = float(np.std(adj_xor))
    features['adj_xor_deviation'] = float(np.sum(np.abs(adj_xor - 0.5)))

    # === 4. Word-level analysis (8 слов по 32 бит) ===
    word_hws = np.zeros((n, 8))
    for idx, x in enumerate(all_ints):
        for w in range(8):
            word = (x >> (w * 32)) & 0xFFFFFFFF
            word_hws[idx, w] = bin(word).count('1')

    # Cross-word correlation
    word_means = np.mean(word_hws, axis=0)
    word_stds = np.std(word_hws, axis=0)

    # a-path words (H0-H3) vs e-path words (H4-H7) asymmetry
    a_path_std = np.mean(word_stds[:4])
    e_path_std = np.mean(word_stds[4:])
    features['word_asymmetry'] = float(abs(a_path_std - e_path_std))

    # Word-word correlation matrix
    if n > 10:
        corr_matrix = np.corrcoef(word_hws.T)
        # Среднее off-diagonal correlation
        mask = ~np.eye(8, dtype=bool)
        features['word_corr_mean'] = float(np.mean(np.abs(corr_matrix[mask])))
    else:
        features['word_corr_mean'] = 0.0

    # === 5. Byte entropy ===
    byte_counts = np.zeros(256)
    for h in hash_bytes_list:
        for b in h:
            byte_counts[b] += 1
    byte_probs = byte_counts / (n * 32)
    byte_probs = byte_probs[byte_probs > 0]
    features['byte_entropy'] = float(-np.sum(byte_probs * np.log2(byte_probs)))

    # === 6. Last word (H[7]) specific ===
    h7_values = [(x & 0xFFFFFFFF) for x in all_ints]
    h7_b29 = sum((v >> 29) & 1 for v in h7_values) / n
    h7_b30 = sum((v >> 30) & 1 for v in h7_values) / n
    features['h7_b29'] = h7_b29
    features['h7_b30'] = h7_b30
    features['h7_bias'] = abs(h7_b29 - 0.5) + abs(h7_b30 - 0.5)

    # === 7. Run length statistics ===
    run_lengths = []
    for x in all_ints[:min(n, 1000)]:
        current = (x >> 0) & 1
        length = 1
        for b in range(1, 256):
            bit = (x >> b) & 1
            if bit == current:
                length += 1
            else:
                run_lengths.append(length)
                current = bit
                length = 1
        run_lengths.append(length)

    features['run_mean'] = float(np.mean(run_lengths))
    features['run_max_mean'] = float(np.mean([max(1, r) for r in run_lengths]))

    return features


def compute_score(features_sha, features_rand):
    """Вычисляем score различия между двумя наборами фичей."""
    score = 0
    details = {}

    for key in features_sha:
        if key in features_rand:
            diff = abs(features_sha[key] - features_rand[key])
            # Нормируем по значению random (baseline)
            baseline = abs(features_rand[key]) + 1e-10
            normalized = diff / baseline
            score += normalized
            details[key] = (features_sha[key], features_rand[key], diff, normalized)

    return score, details


def experiment():
    np.random.seed(42)

    print("=" * 70)
    print("ПОЛНЫЙ DISTINGUISHER НА ECDSA NONCE")
    print("=" * 70)

    privkey = os.urandom(32)

    # =================================================================
    print("\n" + "=" * 70)
    print("1. СБОР ФИЧЕЙ: SHA-256 nonce vs Random nonce")
    print("=" * 70)

    for n_sigs in [100, 500, 1000, 5000, 10000]:
        # SHA-256 nonce (уязвимый)
        sha_hashes = []
        for i in range(n_sigs):
            msg = f"tx_{i}".encode()
            k, h = generate_nonce_vulnerable(privkey, msg)
            sha_hashes.append(h)

        # Random nonce
        rand_hashes = []
        for i in range(n_sigs):
            msg = f"tx_{i}".encode()
            k, h = generate_nonce_random(msg)
            rand_hashes.append(h)

        feat_sha = extract_features(sha_hashes)
        feat_rand = extract_features(rand_hashes)
        score, details = compute_score(feat_sha, feat_rand)

        print(f"\n  N = {n_sigs:,}:")
        print(f"    Total score: {score:.4f}")

        # Топ-5 различающих фичей
        sorted_details = sorted(details.items(), key=lambda x: x[1][3], reverse=True)
        print(f"    Топ-5 фичей:")
        for key, (sha_val, rand_val, diff, norm) in sorted_details[:5]:
            direction = "SHA>rand" if sha_val > rand_val else "SHA<rand"
            print(f"      {key:25s}: SHA={sha_val:.5f} rand={rand_val:.5f} "
                  f"Δ={diff:.5f} ({direction})")

    # =================================================================
    print("\n" + "=" * 70)
    print("2. BLIND TEST: отличить SHA-256 от random")
    print("=" * 70)

    # Генерируем 100 наборов: 50 SHA, 50 random
    # Для каждого — N подписей. Классифицируем.

    for n_sigs in [100, 500, 1000, 5000]:
        correct = 0
        total = 100

        # Сначала собираем reference из большой выборки
        ref_sha = []
        ref_rand = []
        for i in range(10000):
            msg = f"ref_{i}".encode()
            _, h_sha = generate_nonce_vulnerable(privkey, msg)
            _, h_rand = generate_nonce_random(msg)
            ref_sha.append(h_sha)
            ref_rand.append(h_rand)

        feat_ref_sha = extract_features(ref_sha)
        feat_ref_rand = extract_features(ref_rand)

        for trial in range(total):
            is_sha = trial < 50

            hashes = []
            for i in range(n_sigs):
                msg = f"test_{trial}_{i}".encode()
                if is_sha:
                    _, h = generate_nonce_vulnerable(privkey, msg)
                else:
                    _, h = generate_nonce_random(msg)
                hashes.append(h)

            feat = extract_features(hashes)

            # Score: ближе к SHA или к random?
            dist_sha = sum(abs(feat.get(k, 0) - feat_ref_sha.get(k, 0))
                          for k in feat_ref_sha)
            dist_rand = sum(abs(feat.get(k, 0) - feat_ref_rand.get(k, 0))
                           for k in feat_ref_rand)

            predicted_sha = dist_sha < dist_rand
            if predicted_sha == is_sha:
                correct += 1

        acc = correct / total * 100
        print(f"    N={n_sigs:5d}: accuracy = {acc:.1f}%  "
              f"{'← РАБОТАЕТ!' if acc > 65 else ''}")

    # =================================================================
    print("\n" + "=" * 70)
    print("3. ЧТО ОТЛИЧАЕТ SHA-256 HASH ОТ RANDOM?")
    print("=" * 70)

    # Большая выборка для точных измерений
    N_LARGE = 50000
    sha_large = []
    rand_large = []
    for i in range(N_LARGE):
        msg = f"large_{i}".encode()
        _, h_sha = generate_nonce_vulnerable(privkey, msg)
        sha_large.append(h_sha)
        rand_large.append(os.urandom(32))

    feat_sha = extract_features(sha_large)
    feat_rand = extract_features(rand_large)

    print(f"\n  50K сэмплов — полная карта отличий:")
    print(f"\n  {'Фича':<28} {'SHA-256':>10} {'Random':>10} {'Δ':>10} {'σ':>8}")
    print(f"  {'-'*28} {'-'*10} {'-'*10} {'-'*10} {'-'*8}")

    _, details = compute_score(feat_sha, feat_rand)
    for key, (sv, rv, diff, norm) in sorted(details.items(), key=lambda x: x[1][3], reverse=True):
        # Оценка σ (стандартное отклонение при random)
        sigma_est = abs(rv) * 0.01 + 0.001  # грубая оценка
        z_score = diff / sigma_est if sigma_est > 0 else 0
        significant = "★" if abs(z_score) > 3 else ""
        print(f"  {key:<28} {sv:10.5f} {rv:10.5f} {diff:10.5f} {z_score:7.1f} {significant}")

    # =================================================================
    print("\n" + "=" * 70)
    print("4. ВЕРДИКТ")
    print("=" * 70)


if __name__ == "__main__":
    experiment()
