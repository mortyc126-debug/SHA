"""
Точный ALG distinguisher из методички — применяем к ECDSA.

Ключевое понимание:
  Distinguisher из методички сравнивает:
    SET A: SHA-256(random inputs) — N хешей
    SET B: random 32-byte strings — N строк
  И детектит что SET A имеет carry-fingerprint.

  Для ECDSA:
    SET A: SHA-256(privkey || msg_i) — N хешей от одного кошелька
    SET B: random 32-byte strings — N строк
  Вопрос: видим ли carry-fingerprint в SET A?

  Критически: SET A — это SHA-256 outputs! Значит fingerprint ДОЛЖЕН быть.
  Наша ошибка была в фичах. Используем ТОЧНЫЕ фичи из методички.
"""

import hashlib
import struct
import os
import math
import random
import numpy as np

def sha256_hash(msg_bytes):
    return hashlib.sha256(msg_bytes).digest()

def hash_to_words(h):
    return struct.unpack('>8I', h)

def bit(word, pos):
    return (word >> pos) & 1


def compute_alg_score(hash_list):
    """
    Точная копия ALG distinguisher score из методички.
    Вход: список 32-byte hash outputs.
    Выход: score (SHA-256 → высокий, random → низкий).
    """
    N = len(hash_list)
    words_list = [hash_to_words(h) for h in hash_list]

    # === Component 1: Bit correlation profile (carry fingerprint) ===
    # Adjacent bit XOR в первом слове H[0]
    xor_count = [0] * 31
    for h_words in words_list:
        w0 = h_words[0]
        for i in range(31):
            xor_count[i] += bit(w0, i) ^ bit(w0, i + 1)

    biases = [abs(xor_count[i] / N - 0.5) for i in range(31)]
    mean_bias = sum(biases) / len(biases)
    variance = sum((b - mean_bias)**2 for b in biases) / len(biases)
    std_bias = math.sqrt(variance)
    uniformity = std_bias / mean_bias if mean_bias > 0 else 0

    # === Component 2: Cross-word asymmetry ===
    word_hw = [[] for _ in range(8)]
    for h_words in words_list:
        for i in range(8):
            word_hw[i].append(bin(h_words[i]).count('1'))

    stds = []
    for i in range(8):
        mean = sum(word_hw[i]) / N
        std = math.sqrt(sum((x - mean)**2 for x in word_hw[i]) / N)
        stds.append(std)

    a_std = sum(stds[:4]) / 4
    e_std = sum(stds[4:]) / 4
    word_ratio = e_std / a_std if a_std > 0 else 1.0

    # === Combined score (точная формула из методички) ===
    score = uniformity * 1000 + abs(1.0 - word_ratio) * 10000

    return score, uniformity, word_ratio


def experiment():
    print("=" * 70)
    print("ТОЧНЫЙ ALG SCORE: SHA-256 vs Random vs ECDSA nonce")
    print("=" * 70)

    # =================================================================
    print("\n" + "=" * 70)
    print("1. CALIBRATION: score для SHA-256 vs Random (как в методичке)")
    print("=" * 70)

    for N in [1000, 5000, 10000, 50000, 200000]:
        # SHA-256 с random input (как в методичке)
        sha_hashes = [sha256_hash(random.getrandbits(440).to_bytes(55, 'big')) for _ in range(N)]
        # True random bytes
        rand_hashes = [os.urandom(32) for _ in range(N)]

        score_sha, uni_sha, wr_sha = compute_alg_score(sha_hashes)
        score_rand, uni_rand, wr_rand = compute_alg_score(rand_hashes)

        print(f"\n  N={N:>7,}:")
        print(f"    SHA-256: score={score_sha:8.2f}  uniformity={uni_sha:.6f}  word_ratio={wr_sha:.6f}")
        print(f"    Random:  score={score_rand:8.2f}  uniformity={uni_rand:.6f}  word_ratio={wr_rand:.6f}")
        print(f"    Δscore:  {score_sha - score_rand:+.2f}  {'← SHA detectable' if score_sha > score_rand * 1.5 else ''}")

    # =================================================================
    print("\n" + "=" * 70)
    print("2. ECDSA NONCE: score для k=SHA256(privkey||msg)")
    print("=" * 70)

    privkey = os.urandom(32)

    for N in [1000, 5000, 10000, 50000, 200000]:
        # ECDSA vulnerable nonce
        ecdsa_hashes = [sha256_hash(privkey + f"tx_{i}".encode()) for i in range(N)]
        # Random nonce
        rand_hashes = [os.urandom(32) for _ in range(N)]

        score_ecdsa, uni_ecdsa, wr_ecdsa = compute_alg_score(ecdsa_hashes)
        score_rand, uni_rand, wr_rand = compute_alg_score(rand_hashes)

        print(f"\n  N={N:>7,}:")
        print(f"    ECDSA:   score={score_ecdsa:8.2f}  uniformity={uni_ecdsa:.6f}  word_ratio={wr_ecdsa:.6f}")
        print(f"    Random:  score={score_rand:8.2f}  uniformity={uni_rand:.6f}  word_ratio={wr_rand:.6f}")
        print(f"    Δscore:  {score_ecdsa - score_rand:+.2f}  {'← ECDSA detectable' if score_ecdsa > score_rand * 1.5 else ''}")

    # =================================================================
    print("\n" + "=" * 70)
    print("3. BLIND TEST: 50 SHA + 50 Random → классификация")
    print("=" * 70)

    for N in [1000, 5000, 10000, 50000]:
        # Calibration
        cal_sha_scores = []
        cal_rand_scores = []
        for _ in range(20):
            sh = [sha256_hash(random.getrandbits(440).to_bytes(55, 'big')) for _ in range(N)]
            rh = [os.urandom(32) for _ in range(N)]
            cal_sha_scores.append(compute_alg_score(sh)[0])
            cal_rand_scores.append(compute_alg_score(rh)[0])

        threshold = (np.mean(cal_sha_scores) + np.mean(cal_rand_scores)) / 2

        # Test: 50 SHA-256 + 50 random
        tp, tn, fp, fn = 0, 0, 0, 0
        for trial in range(100):
            is_sha = trial < 50
            if is_sha:
                hashes = [sha256_hash(random.getrandbits(440).to_bytes(55, 'big')) for _ in range(N)]
            else:
                hashes = [os.urandom(32) for _ in range(N)]

            score = compute_alg_score(hashes)[0]
            predicted_sha = score > threshold

            if is_sha and predicted_sha: tp += 1
            elif is_sha: fn += 1
            elif predicted_sha: fp += 1
            else: tn += 1

        acc = (tp + tn) / 100 * 100
        print(f"    N={N:>6,}: accuracy={acc:.0f}% (TP={tp} TN={tn} FP={fp} FN={fn}) "
              f"thresh={threshold:.2f}  {'★★★' if acc > 90 else '★★' if acc > 70 else '★' if acc > 60 else ''}")

    # =================================================================
    print("\n" + "=" * 70)
    print("4. BLIND TEST: ECDSA nonce (фикс privkey)")
    print("=" * 70)

    for N in [1000, 5000, 10000, 50000]:
        # Same threshold from section 3 calibration
        cal_sha_scores = []
        cal_rand_scores = []
        for _ in range(20):
            sh = [sha256_hash(random.getrandbits(440).to_bytes(55, 'big')) for _ in range(N)]
            rh = [os.urandom(32) for _ in range(N)]
            cal_sha_scores.append(compute_alg_score(sh)[0])
            cal_rand_scores.append(compute_alg_score(rh)[0])
        threshold = (np.mean(cal_sha_scores) + np.mean(cal_rand_scores)) / 2

        tp, tn, fp, fn = 0, 0, 0, 0
        for trial in range(100):
            is_sha = trial < 50
            pk = os.urandom(32)  # разный privkey для каждого!
            if is_sha:
                hashes = [sha256_hash(pk + f"tx_{i}".encode()) for i in range(N)]
            else:
                hashes = [os.urandom(32) for _ in range(N)]

            score = compute_alg_score(hashes)[0]
            predicted_sha = score > threshold

            if is_sha and predicted_sha: tp += 1
            elif is_sha: fn += 1
            elif predicted_sha: fp += 1
            else: tn += 1

        acc = (tp + tn) / 100 * 100
        print(f"    N={N:>6,}: accuracy={acc:.0f}% (TP={tp} TN={tn} FP={fp} FN={fn}) "
              f"{'★★★' if acc > 90 else '★★' if acc > 70 else '★' if acc > 60 else ''}")

    # =================================================================
    print("\n" + "=" * 70)
    print("5. ВЕРДИКТ")
    print("=" * 70)


if __name__ == "__main__":
    experiment()
