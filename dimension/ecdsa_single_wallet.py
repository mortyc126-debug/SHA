"""
ECDSA distinguisher: ОДИН кошелёк, МНОГО подписей.

Сценарий реальной атаки:
  Атакующий видит кошелёк, делающий много транзакций.
  Вопрос: k = SHA256(priv||msg) или k = random?

  1. Собираем N хешей от ОДНОГО privkey
  2. Тестируем КАЖДЫЙ distinguisher feature
  3. Определяем порог N для уверенной детекции
  4. Без knowledge of privkey — только статистика выхода
"""

import hashlib
import os
import numpy as np
from scipy import stats as sp_stats

N_CURVE = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141

def sha256_bytes(data):
    return hashlib.sha256(data).digest()

def hw(x):
    return bin(x).count('1')

def h_to_words(h_bytes):
    return [int.from_bytes(h_bytes[i*4:(i+1)*4], 'big') for i in range(8)]


def collect_hashes(privkey, n, mode='sha'):
    """Собираем n хешей от одного privkey."""
    hashes = []
    for i in range(n):
        msg = f"payment_{i}_{os.urandom(8).hex()}".encode()
        if mode == 'sha':
            h = sha256_bytes(privkey + msg)
        else:
            h = os.urandom(32)
        hashes.append(h)
    return hashes


def compute_features_single(hash_list):
    """Полный набор фичей для одного набора хешей."""
    n = len(hash_list)
    ints = [int.from_bytes(h, 'big') for h in hash_list]
    words_list = [h_to_words(h) for h in hash_list]

    features = {}

    # 1. Adjacent bit XOR — самая сильная фича (z=8.8 на 50K)
    adj_counts = np.zeros(255)
    for x in ints:
        for b in range(255):
            adj_counts[b] += ((x >> b) & 1) ^ ((x >> (b+1)) & 1)
    adj_counts /= n
    features['adj_xor_deviation'] = float(np.sum(np.abs(adj_counts - 0.5)))
    features['adj_xor_variance'] = float(np.var(adj_counts))

    # 2. Per-word HW statistics
    word_hws = np.array([[hw(w[i]) for i in range(8)] for w in words_list])
    for i in range(8):
        features[f'word{i}_hw_mean'] = float(np.mean(word_hws[:, i]))
        features[f'word{i}_hw_std'] = float(np.std(word_hws[:, i]))

    # 3. Cross-word correlations
    if n > 20:
        corr = np.corrcoef(word_hws.T)
        features['cross_word_corr'] = float(np.mean(np.abs(corr[np.triu_indices(8, k=1)])))
    else:
        features['cross_word_corr'] = 0.0

    # 4. e-path vs a-path asymmetry
    features['ea_asymmetry'] = float(np.mean(word_hws[:, 4:]) - np.mean(word_hws[:, :4]))

    # 5. Bit frequency chi-squared
    bit_freq = np.zeros(256)
    for x in ints:
        for b in range(256):
            bit_freq[b] += (x >> b) & 1
    bit_freq /= n
    features['bit_chi2'] = float(np.sum((bit_freq - 0.5)**2 / 0.25 * n))

    # 6. Run length distribution
    all_runs = []
    for x in ints[:min(n, 2000)]:
        cur = x & 1
        length = 1
        for b in range(1, 256):
            if ((x >> b) & 1) == cur:
                length += 1
            else:
                all_runs.append(length)
                cur = 1 - cur
                length = 1
        all_runs.append(length)
    features['run_mean'] = float(np.mean(all_runs))
    features['run_var'] = float(np.var(all_runs))

    # 7. Byte distribution chi-squared
    byte_counts = np.zeros(256)
    for h in hash_list:
        for b in h:
            byte_counts[b] += 1
    expected = n * 32 / 256
    features['byte_chi2'] = float(np.sum((byte_counts - expected)**2 / expected))

    # 8. Serial correlation (lag-1 between consecutive hashes)
    if n > 2:
        h_ints = [int.from_bytes(h, 'big') % 2**32 for h in hash_list]
        features['serial_corr'] = float(np.corrcoef(h_ints[:-1], h_ints[1:])[0, 1])
    else:
        features['serial_corr'] = 0.0

    return features


def test_single_wallet():
    np.random.seed(42)

    print("=" * 70)
    print("ОДИН КОШЕЛЁК, МНОГО ПОДПИСЕЙ — single-wallet distinguisher")
    print("=" * 70)

    # Фиксируем ОДИН privkey
    privkey = os.urandom(32)

    # =================================================================
    print("\n" + "=" * 70)
    print("1. НАКОПЛЕНИЕ СИГНАЛА: один privkey, растущее N")
    print("=" * 70)

    for N in [100, 500, 1000, 5000, 10000, 50000]:
        h_sha = collect_hashes(privkey, N, 'sha')
        h_rand = collect_hashes(privkey, N, 'random')

        f_sha = compute_features_single(h_sha)
        f_rand = compute_features_single(h_rand)

        # Какие фичи отличаются?
        diffs = {}
        for k in f_sha:
            if k in f_rand:
                d = f_sha[k] - f_rand[k]
                diffs[k] = d

        # Топ-3 по абсолютной разнице
        top = sorted(diffs.items(), key=lambda x: abs(x[1]), reverse=True)[:3]

        print(f"\n  N={N:6d}:")
        for feat, diff in top:
            print(f"    {feat:25s}: SHA={f_sha[feat]:10.4f}  rand={f_rand[feat]:10.4f}  Δ={diff:+.4f}")

    # =================================================================
    print("\n" + "=" * 70)
    print("2. STATISTICAL TEST: p-value для каждой фичи")
    print("=" * 70)

    # Для каждой фичи: bootstrap p-value
    # Собираем 50K SHA хешей и 50K random, делим на подвыборки
    N_TOTAL = 50000
    h_sha_all = collect_hashes(privkey, N_TOTAL, 'sha')
    h_rand_all = collect_hashes(privkey, N_TOTAL, 'random')

    # Bootstrap: разбиваем на 100 подвыборок по 500
    n_bootstrap = 100
    sub_size = 500

    feature_distributions = {}

    for b in range(n_bootstrap):
        idx = np.random.choice(N_TOTAL, sub_size, replace=False)
        h_sub_sha = [h_sha_all[i] for i in idx]
        h_sub_rand = [h_rand_all[i] for i in idx]

        f_sha = compute_features_single(h_sub_sha)
        f_rand = compute_features_single(h_sub_rand)

        for k in f_sha:
            if k not in feature_distributions:
                feature_distributions[k] = {'sha': [], 'rand': []}
            feature_distributions[k]['sha'].append(f_sha[k])
            feature_distributions[k]['rand'].append(f_rand[k])

    print(f"\n  Bootstrap: {n_bootstrap} подвыборок × {sub_size} хешей")
    print(f"\n  {'Feature':<28} {'SHA mean':>10} {'Rand mean':>10} {'t-stat':>8} {'p-value':>10} {'Sig':>5}")
    print(f"  {'-'*28} {'-'*10} {'-'*10} {'-'*8} {'-'*10} {'-'*5}")

    significant_features = []
    for k in sorted(feature_distributions.keys()):
        sha_vals = feature_distributions[k]['sha']
        rand_vals = feature_distributions[k]['rand']

        t_stat, p_val = sp_stats.ttest_ind(sha_vals, rand_vals)
        m_sha = np.mean(sha_vals)
        m_rand = np.mean(rand_vals)

        sig = ""
        if p_val < 0.001: sig = "★★★"
        elif p_val < 0.01: sig = "★★"
        elif p_val < 0.05: sig = "★"

        if sig or k in ['adj_xor_deviation', 'byte_chi2', 'bit_chi2', 'ea_asymmetry']:
            print(f"  {k:<28} {m_sha:10.4f} {m_rand:10.4f} {t_stat:+7.2f} {p_val:10.6f} {sig}")
            if sig:
                significant_features.append((k, t_stat, p_val))

    # =================================================================
    print("\n" + "=" * 70)
    print("3. BLIND CLASSIFICATION: один кошелёк, разное N")
    print("=" * 70)

    # Для ОДНОГО privkey: N подписей → SHA или random?
    # Используем только significant features

    if not significant_features:
        # Используем все фичи с наибольшим t-stat
        all_feats = []
        for k in feature_distributions:
            sha_v = feature_distributions[k]['sha']
            rand_v = feature_distributions[k]['rand']
            t, p = sp_stats.ttest_ind(sha_v, rand_v)
            all_feats.append((k, abs(t), p))
        all_feats.sort(key=lambda x: x[1], reverse=True)
        significant_features = [(k, t, p) for k, t, p in all_feats[:5]]

    best_feat = significant_features[0][0] if significant_features else 'adj_xor_deviation'
    print(f"\n  Best feature: {best_feat} (t={significant_features[0][1]:.2f}, p={significant_features[0][2]:.6f})")

    # Threshold из bootstrap
    sha_vals = feature_distributions[best_feat]['sha']
    rand_vals = feature_distributions[best_feat]['rand']
    threshold = (np.mean(sha_vals) + np.mean(rand_vals)) / 2
    sha_higher = np.mean(sha_vals) > np.mean(rand_vals)

    print(f"  Threshold: {threshold:.4f} (SHA {'>' if sha_higher else '<'} rand)")

    # Тест на НОВЫХ данных от ТОГО ЖЕ кошелька
    print(f"\n  {'N':>7} {'Accuracy':>10} {'P(detect SHA)':>14} {'P(reject rand)':>16}")

    for N in [50, 100, 200, 500, 1000, 2000, 5000]:
        correct_sha = 0
        correct_rand = 0
        n_trials = 100

        for trial in range(n_trials):
            h_sha = collect_hashes(privkey, N, 'sha')
            h_rand = collect_hashes(privkey, N, 'random')

            f_sha = compute_features_single(h_sha)
            f_rand = compute_features_single(h_rand)

            if sha_higher:
                if f_sha[best_feat] > threshold: correct_sha += 1
                if f_rand[best_feat] <= threshold: correct_rand += 1
            else:
                if f_sha[best_feat] < threshold: correct_sha += 1
                if f_rand[best_feat] >= threshold: correct_rand += 1

        acc_sha = correct_sha / n_trials * 100
        acc_rand = correct_rand / n_trials * 100
        total_acc = (correct_sha + correct_rand) / (2 * n_trials) * 100

        status = ""
        if total_acc > 70: status = "★"
        if total_acc > 80: status = "★★"
        if total_acc > 90: status = "★★★"

        print(f"  {N:7d} {total_acc:9.1f}% {acc_sha:13.1f}% {acc_rand:15.1f}% {status}")

    # =================================================================
    print("\n" + "=" * 70)
    print("4. ВОСПРОИЗВОДИМОСТЬ: работает ли для ДРУГОГО privkey?")
    print("=" * 70)

    privkey2 = os.urandom(32)
    privkey3 = os.urandom(32)

    for pk_name, pk in [("Original", privkey), ("Privkey 2", privkey2), ("Privkey 3", privkey3)]:
        h_sha = collect_hashes(pk, 5000, 'sha')
        h_rand = collect_hashes(pk, 5000, 'random')

        f_sha = compute_features_single(h_sha)
        f_rand = compute_features_single(h_rand)

        pred_sha = (f_sha[best_feat] > threshold) if sha_higher else (f_sha[best_feat] < threshold)
        pred_rand = (f_rand[best_feat] <= threshold) if sha_higher else (f_rand[best_feat] >= threshold)

        print(f"  {pk_name:12s}: SHA={f_sha[best_feat]:.4f} → {'✓ detected' if pred_sha else '✗ missed'}  |  "
              f"Rand={f_rand[best_feat]:.4f} → {'✓ rejected' if pred_rand else '✗ false alarm'}")


if __name__ == "__main__":
    test_single_wallet()
