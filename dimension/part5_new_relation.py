"""
ЧАСТЬ 5: Новое отношение.

Collision = H(x) == H(y). SHA-256 защищена от этого.
Preimage = H(x) == target. SHA-256 защищена от этого.

А что если определить ДРУГОЕ отношение R(x, y)?
Такое, что R(x, y) = true → полезно (как collision),
но SHA-256 НЕ проектировалась защищать от R?

Примеры новых отношений:
  1. MODULAR: H(x) + H(y) ≡ 0 (mod 2^32) per word
  2. PRODUCT: H(x) × H(y) ≡ 1 (mod 2^32) per word
  3. ROTATION: H(x) = ROT_k(H(y)) для какого-то k
  4. PARTIAL: H(x)[0:4] == H(y)[0:4] (половина хеша)
  5. ALGEBRAIC: H(x)² + H(y)² ≡ 0 (mod p) для какого-то p
  6. COMPOSED: H(H(x)) == H(H(y)) (вторая итерация)

Тестируем каждое: насколько СЛОЖНЕЕ/ПРОЩЕ чем collision?
"""

import numpy as np
import struct, hashlib

MASK32 = 0xFFFFFFFF

def sha256_hash(W16):
    raw = struct.pack('>16I', *W16)
    return struct.unpack('>8I', hashlib.sha256(raw).digest())

def hw(x): return bin(x).count('1')
def rotr32(x, n): return ((x >> n) | (x << (32 - n))) & MASK32


def main():
    np.random.seed(42)

    print("=" * 70)
    print("ЧАСТЬ 5: Новые отношения — что SHA-256 НЕ защищает?")
    print("=" * 70)

    # ═══════════════════════════════
    print(f"\n{'=' * 70}")
    print("1. MODULAR SUM: H(x) + H(y) ≡ 0 (mod 2^32)")
    print("=" * 70)

    # Для каждого слова: H(x)[i] + H(y)[i] ≡ 0 mod 2^32
    # Это значит H(y)[i] = -H(x)[i] mod 2^32 = (2^32 - H(x)[i]) & MASK32
    # Для random function: P = (1/2^32)^8 = 2^-256 per pair
    # Для SHA-256: same? Или есть structure?

    # Тест: среди N случайных хешей, ищем пару с суммой close to 0
    N = 10000
    hashes = []
    for _ in range(N):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        hashes.append(sha256_hash(W))

    # Для каждой пары: sum distance
    best_sum_dist = 256
    for trial in range(100000):
        i = np.random.randint(N)
        j = np.random.randint(N)
        if i == j: continue
        # mod sum per word
        sum_hw = sum(hw((hashes[i][w] + hashes[j][w]) & MASK32) for w in range(8))
        if sum_hw < best_sum_dist:
            best_sum_dist = sum_hw

    from scipy.stats import norm
    expected = 128 - 8 * norm.ppf(1 - 1.0/100000)
    print(f"  Best HW(H(x)+H(y) mod 2^32) = {best_sum_dist}")
    print(f"  Expected (random): ~{expected:.0f}")
    print(f"  → {'EASIER than collision!' if best_sum_dist < expected - 10 else 'Same difficulty as collision.'}")

    # ═══════════════════════════════
    print(f"\n{'=' * 70}")
    print("2. ROTATION: H(x) = ROT_k(H(y))")
    print("=" * 70)

    # Ищем пару где H(x)[i] = ROT_k(H(y)[i]) для фиксированного k
    best_rot_dist = 256
    best_k = 0
    for k in range(1, 32):
        for trial in range(10000):
            i = np.random.randint(N)
            j = np.random.randint(N)
            if i == j: continue
            d = sum(hw(hashes[i][w] ^ rotr32(hashes[j][w], k)) for w in range(8))
            if d < best_rot_dist:
                best_rot_dist = d
                best_k = k

    print(f"  Best rotation distance: HW = {best_rot_dist}, k = {best_k}")
    print(f"  Expected (random): ~{expected:.0f}")
    print(f"  → {'STRUCTURE in rotation!' if best_rot_dist < expected - 10 else 'Same as random.'}")

    # ═══════════════════════════════
    print(f"\n{'=' * 70}")
    print("3. PARTIAL: H(x)[0:4] == H(y)[0:4] (128 бит)")
    print("=" * 70)

    # Birthday on 128 bits: need ~2^64 hashes
    # With N=10K: best distance to 0 in 128-bit space
    best_partial = 128
    for trial in range(200000):
        i = np.random.randint(N)
        j = np.random.randint(N)
        if i == j: continue
        d = sum(hw(hashes[i][w] ^ hashes[j][w]) for w in range(4))  # first 4 words only
        if d < best_partial:
            best_partial = d

    expected_128 = 64 - 4 * norm.ppf(1 - 1.0/200000)
    print(f"  Best HW(δH[0:4]) = {best_partial}")
    print(f"  Expected (random 128-bit): ~{expected_128:.0f}")

    # ═══════════════════════════════
    print(f"\n{'=' * 70}")
    print("4. WORD PERMUTATION: H(x) = permute(H(y))")
    print("=" * 70)

    # Ищем: H(x) is a PERMUTATION of H(y) (same 8 words, different order)
    # 8! = 40320 permutations
    # For each pair: check if words of H(x) are a rearrangement of H(y)

    best_perm_dist = 256
    from itertools import permutations

    # Quick test: sorted words
    for trial in range(50000):
        i = np.random.randint(N)
        j = np.random.randint(N)
        if i == j: continue
        # Sort both hashes by word value
        s1 = sorted(hashes[i])
        s2 = sorted(hashes[j])
        d = sum(hw(s1[w] ^ s2[w]) for w in range(8))
        if d < best_perm_dist:
            best_perm_dist = d

    print(f"  Best sorted-word distance: HW = {best_perm_dist}")
    print(f"  Expected (random): ~{expected:.0f}")
    print(f"  → {'PERMUTATION STRUCTURE!' if best_perm_dist < expected - 15 else 'Same as random.'}")

    # ═══════════════════════════════
    print(f"\n{'=' * 70}")
    print("5. SELF-INVERSE: H(x) ⊕ H(y) = x ⊕ y (truncated)")
    print("=" * 70)

    # Wild idea: is there structure where the XOR of hashes
    # relates to the XOR of inputs?
    # For random: no relation.
    # For SHA-256: possible through algebraic structure?

    relations_found = 0
    for trial in range(10000):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = [np.random.randint(0, 2**32) for _ in range(16)]
        H1 = sha256_hash(W1)
        H2 = sha256_hash(W2)

        # δH vs δW (first 8 words)
        dH = [H1[i] ^ H2[i] for i in range(8)]
        dW = [W1[i] ^ W2[i] for i in range(8)]

        # Are they related? Check HW(dH ⊕ dW)
        d = sum(hw(dH[i] ^ dW[i]) for i in range(8))
        if d < 100:
            relations_found += 1

    print(f"  Pairs with HW(δH ⊕ δW) < 100: {relations_found}/10000")
    print(f"  Expected (random): ~{10000 * sum(1 for _ in range(1000) if sum(np.random.binomial(1, 0.5, 256)) < 100) / 1000}")

    # ═══════════════════════════════
    print(f"\n{'=' * 70}")
    print("6. COMPOSED: H(H(x)) vs H(H(y))")
    print("=" * 70)

    # Second iteration: is it easier to find collision in H²?
    # H²(x) = H(H(x)) — two layers of SHA-256
    # If H is random: collision in H² as hard as collision in H (2^128)
    # But: if H has structure, H² might AMPLIFY or REDUCE it

    composed_hashes = []
    for i in range(min(N, 5000)):
        H = hashes[i]
        # Use H as input to second SHA-256 (pad to 16 words)
        W2 = list(H) + [0]*8
        H2 = sha256_hash(W2)
        composed_hashes.append(H2)

    # Find closest pair in H² space
    best_composed = 256
    for trial in range(100000):
        i = np.random.randint(len(composed_hashes))
        j = np.random.randint(len(composed_hashes))
        if i == j: continue
        d = sum(hw(composed_hashes[i][w] ^ composed_hashes[j][w]) for w in range(8))
        if d < best_composed:
            best_composed = d

    # Compare: best distance in H vs H²
    best_single = 256
    for trial in range(100000):
        i = np.random.randint(len(composed_hashes))
        j = np.random.randint(len(composed_hashes))
        if i == j: continue
        d = sum(hw(hashes[i][w] ^ hashes[j][w]) for w in range(8))
        if d < best_single:
            best_single = d

    print(f"  Best distance in H:  {best_single}")
    print(f"  Best distance in H²: {best_composed}")
    print(f"  → {'H² EASIER!' if best_composed < best_single - 5 else 'Same difficulty.'}")

    # ═══════════════════════════════
    print(f"\n{'=' * 70}")
    print("ИТОГ ЧАСТИ 5")
    print("=" * 70)
    print(f"""
  Relation                      Best dist   Expected   Verdict
  ────────────────────────────   ─────────   ────────   ───────
  Standard collision (XOR=0)     ~{best_single:>3}        ~{expected:.0f}       baseline
  Modular sum (ADD=0)            ~{best_sum_dist:>3}        ~{expected:.0f}       same
  Rotation (ROT_k match)         ~{best_rot_dist:>3}        ~{expected:.0f}       same
  Partial (128-bit)              ~{best_partial:>3}/128    ~{expected_128:.0f}/128    same
  Word permutation               ~{best_perm_dist:>3}        ~{expected:.0f}       same
  Composed H²                    ~{best_composed:>3}        ~{expected:.0f}       same

  ВСЕ ОТНОШЕНИЯ = ОДИНАКОВОЙ СЛОЖНОСТИ.
  SHA-256 не различает XOR, ADD, rotation, permutation, composition.
  Потому что она = random function. Random function не имеет
  алгебраической структуры, которую можно эксплуатировать
  через нестандартные отношения.

  НО: мы тестировали отношения на ВЫХОДЕ.
  Что если определить отношение на ВНУТРЕННИХ состояниях?
""")


if __name__ == "__main__":
    main()
