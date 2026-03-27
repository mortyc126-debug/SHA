"""
ЧАСТЬ 1: Найти ОДНО новое свойство SHA-256.

Стандартные тесты видят:
  - XOR разницы (differential)
  - линейные корреляции (linear)
  - алгебраическую степень (higher-order)

Все они говорят: "SHA-256 = random".

Но SHA-256 ≠ random. У неё есть ADD mod 2^32.
ADD ≠ XOR. Carry = НАПРАВЛЕННАЯ операция (LSB → MSB).
Стандартная криптография работает в GF(2) — XOR мире.
SHA-256 живёт в Z/2^32 — модулярном мире.

Вопрос: есть ли свойство в Z/2^32, которое XOR-тесты не видят?
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
    print("ЧАСТЬ 1: ищем свойство, невидимое стандартным тестам")
    print("=" * 70)

    # ═══════════════════════════════════════════════
    # ТЕСТ A: Модулярное расстояние vs XOR расстояние
    # ═══════════════════════════════════════════════
    # XOR расстояние: HW(H1 ⊕ H2) — считает РАЗНЫЕ биты
    # MOD расстояние: |H1 - H2| mod 2^32 — считает АРИФМЕТИЧЕСКУЮ разницу
    #
    # Для random function: оба = случайные, некоррелированные
    # Для SHA-256: ADD создаёт СВЯЗЬ между XOR и MOD

    print(f"\n{'=' * 70}")
    print("ТЕСТ A: корреляция между XOR и MOD расстояниями")
    print("=" * 70)

    xor_dists = []
    mod_dists = []

    for _ in range(10000):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = [np.random.randint(0, 2**32) for _ in range(16)]
        H1 = sha256_hash(W1)
        H2 = sha256_hash(W2)

        xor_d = sum(hw(H1[i] ^ H2[i]) for i in range(8))
        mod_d = sum(min((H1[i] - H2[i]) % (2**32), (H2[i] - H1[i]) % (2**32)) for i in range(8))

        xor_dists.append(xor_d)
        mod_dists.append(mod_d)

    corr = np.corrcoef(xor_dists, mod_dists)[0, 1]
    print(f"  SHA-256: corr(XOR_dist, MOD_dist) = {corr:.6f}")

    # Сравнение: для ЧИСТО СЛУЧАЙНЫХ 256-бит чисел
    xor_r = []
    mod_r = []
    for _ in range(10000):
        A = [np.random.randint(0, 2**32) for _ in range(8)]
        B = [np.random.randint(0, 2**32) for _ in range(8)]
        xor_r.append(sum(hw(A[i]^B[i]) for i in range(8)))
        mod_r.append(sum(min((A[i]-B[i])%(2**32), (B[i]-A[i])%(2**32)) for i in range(8)))

    corr_r = np.corrcoef(xor_r, mod_r)[0, 1]
    print(f"  Random:  corr(XOR_dist, MOD_dist) = {corr_r:.6f}")
    print(f"  Разница: {abs(corr - corr_r):.6f}")

    # ═══════════════════════════════════════════════
    # ТЕСТ B: Carry propagation в выходе
    # ═══════════════════════════════════════════════
    # ADD mod 2^32: carry идёт от бита 0 к биту 31.
    # Это значит: бит 31 зависит от ВСЕХ нижних бит.
    # Бит 0 не зависит ни от чего (чистый XOR).
    #
    # Если SHA-256 сохраняет СЛЕДЫ carry:
    #   Бит 0 выхода = более "независимый"
    #   Бит 31 выхода = более "зависимый"

    print(f"\n{'=' * 70}")
    print("ТЕСТ B: бит 0 vs бит 31 — следы carry?")
    print("=" * 70)

    # Для каждого выходного бита: насколько он коррелирован
    # с СОСЕДНИМ выходным битом?
    bit_pair_corr = np.zeros(32)
    n_samples = 20000

    for _ in range(n_samples):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H = sha256_hash(W)

        # Для H[0]: корреляция между битом b и битом b+1
        for b in range(31):
            b0 = (H[0] >> b) & 1
            b1 = (H[0] >> (b+1)) & 1
            if b0 == b1:
                bit_pair_corr[b] += 1

    bit_pair_corr /= n_samples

    print(f"\n  Корреляция соседних бит в H[0]:")
    print(f"  (для random: все ≈ 0.500)")
    for b in range(31):
        deviation = abs(bit_pair_corr[b] - 0.5)
        star = " ★" if deviation > 0.01 else ""
        print(f"    bit({b},{b+1}): {bit_pair_corr[b]:.4f} (deviation: {deviation:.4f}){star}")

    # ═══════════════════════════════════════════════
    # ТЕСТ C: Модулярная сумма регистров
    # ═══════════════════════════════════════════════
    # В SHA-256: финальный хеш = IV + state[64] (mod 2^32 per word)
    # Это значит: H[i] = IV[i] + state[i] mod 2^32
    # Модулярная сумма всех H[i] = sum(IV) + sum(state) mod 2^32
    #
    # Есть ли что-то особенное в sum(H) mod 2^32?

    print(f"\n{'=' * 70}")
    print("ТЕСТ C: модулярная сумма H[0]+H[1]+...+H[7]")
    print("=" * 70)

    mod_sums = []
    for _ in range(50000):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H = sha256_hash(W)
        s = 0
        for h in H:
            s = (s + h) & MASK32
        mod_sums.append(s)

    # Для random: mod_sums uniform в [0, 2^32)
    # Тест: HW(mod_sum) должен быть Binomial(32, 0.5) = mean 16, std 2.83
    hw_sums = [hw(s) for s in mod_sums]
    print(f"  HW(sum(H)): mean={np.mean(hw_sums):.4f}, std={np.std(hw_sums):.4f}")
    print(f"  Expected:    mean=16.0000, std=2.8284")

    # Более тонкий тест: распределение mod_sum mod 256
    mod256 = [s % 256 for s in mod_sums]
    from collections import Counter
    counts = Counter(mod256)
    expected = 50000 / 256
    chi2 = sum((counts.get(i, 0) - expected)**2 / expected for i in range(256))
    print(f"  χ² test (sum mod 256): {chi2:.1f} (expected: 255 ± 22)")
    print(f"  → {'★ NON-UNIFORM!' if abs(chi2 - 255) > 50 else 'Uniform.'}")

    # ═══════════════════════════════════════════════
    # ТЕСТ D: XOR сумма vs MOD сумма
    # ═══════════════════════════════════════════════
    # XOR: H[0] ⊕ H[1] ⊕ ... ⊕ H[7]
    # MOD: (H[0] + H[1] + ... + H[7]) mod 2^32
    # Для random: некоррелированы
    # Для SHA-256: carry в ADD может создать связь

    print(f"\n{'=' * 70}")
    print("ТЕСТ D: XOR-fold vs MOD-fold")
    print("=" * 70)

    xor_folds = []
    mod_folds = []
    for _ in range(50000):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H = sha256_hash(W)

        xf = 0
        mf = 0
        for h in H:
            xf ^= h
            mf = (mf + h) & MASK32
        xor_folds.append(hw(xf))
        mod_folds.append(hw(mf))

    corr_fold = np.corrcoef(xor_folds, mod_folds)[0, 1]

    # Random baseline
    xor_fr = []
    mod_fr = []
    for _ in range(50000):
        vals = [np.random.randint(0, 2**32) for _ in range(8)]
        xf = 0; mf = 0
        for v in vals:
            xf ^= v
            mf = (mf + v) & MASK32
        xor_fr.append(hw(xf))
        mod_fr.append(hw(mf))

    corr_fold_r = np.corrcoef(xor_fr, mod_fr)[0, 1]

    print(f"  SHA-256: corr(HW(XOR-fold), HW(MOD-fold)) = {corr_fold:.6f}")
    print(f"  Random:  corr(HW(XOR-fold), HW(MOD-fold)) = {corr_fold_r:.6f}")
    diff = abs(corr_fold - corr_fold_r)
    print(f"  Разница: {diff:.6f}")
    print(f"  → {'★ РАЗЛИЧИЕ!' if diff > 0.01 else 'Одинаково.'}")

    # ═══════════════════════════════════════════════
    # ТЕСТ E: Bit-position bias при δW = +1 (арифметический)
    # ═══════════════════════════════════════════════
    # XOR flip: W[0] ^= 1 (меняет только бит 0)
    # ADD flip: W[0] += 1 (меняет бит 0 И может carry дальше)
    #
    # Вопрос: реагирует ли SHA-256 ОДИНАКОВО на XOR и ADD?

    print(f"\n{'=' * 70}")
    print("ТЕСТ E: XOR flip vs ADD flip")
    print("=" * 70)

    xor_flip_hws = []
    add_flip_hws = []

    for _ in range(10000):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H_base = sha256_hash(W)

        # XOR flip bit 0
        W_xor = list(W); W_xor[0] ^= 1
        H_xor = sha256_hash(W_xor)
        xor_flip_hws.append(sum(hw(H_base[i]^H_xor[i]) for i in range(8)))

        # ADD +1
        W_add = list(W); W_add[0] = (W_add[0] + 1) & MASK32
        H_add = sha256_hash(W_add)
        add_flip_hws.append(sum(hw(H_base[i]^H_add[i]) for i in range(8)))

    print(f"  XOR flip (W[0]^=1):  mean={np.mean(xor_flip_hws):.2f}, std={np.std(xor_flip_hws):.2f}")
    print(f"  ADD flip (W[0]+=1):  mean={np.mean(add_flip_hws):.2f}, std={np.std(add_flip_hws):.2f}")

    # Ключ: W[0]+=1 и W[0]^=1 совпадают когда бит 0 = 0 (нет carry)
    # Расходятся когда бит 0 = 1 (carry propagates)
    # Если SHA-256 полностью перемешивает carry → оба одинаковы
    # Если нет → разница

    from scipy import stats
    t, p = stats.ttest_ind(xor_flip_hws, add_flip_hws)
    print(f"  T-test: t={t:.3f}, p={p:.6f}")
    print(f"  → {'★ РАЗЛИЧАЮТСЯ!' if p < 0.01 else 'Неразличимы.'}")

    # Подробнее: разделяем по carry/no-carry
    carry_diffs = []
    no_carry_diffs = []
    for trial in range(10000):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H_base = sha256_hash(W)

        W_xor = list(W); W_xor[0] ^= 1
        W_add = list(W); W_add[0] = (W_add[0] + 1) & MASK32

        H_xor = sha256_hash(W_xor)
        H_add = sha256_hash(W_add)

        # δ между XOR-flip и ADD-flip
        d = sum(hw(H_xor[i] ^ H_add[i]) for i in range(8))

        if W[0] & 1:  # carry happens (bit 0 was 1)
            carry_diffs.append(d)
        else:  # no carry (bit 0 was 0, XOR=ADD)
            no_carry_diffs.append(d)

    print(f"\n  Когда carry НЕТ (bit0=0): XOR=ADD → δ should be 0")
    print(f"    Mean δ: {np.mean(no_carry_diffs):.2f} (expected: 0)")
    print(f"  Когда carry ЕСТЬ (bit0=1): XOR≠ADD → δ > 0")
    print(f"    Mean δ: {np.mean(carry_diffs):.2f}")
    print(f"    Min δ:  {min(carry_diffs)}")

    # ═══════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("ИТОГ ЧАСТИ 1")
    print("=" * 70)


if __name__ == "__main__":
    main()
