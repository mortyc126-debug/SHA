"""
Carry-Profile Space — физика в пространстве carry-профилей SHA-256.

Carry-профиль = вектор (carry_count для каждого ADD в каждом раунде).
64 раунда × 7 ADD = 448 компонентов.
Каждый компонент = число carry-бит (0..32).

Вопросы:
  1. Какова размерность пространства профилей?
  2. Какие компоненты детерминированы?
  3. Есть ли кластеры / аттракторы?
  4. Как δW влияет на профиль?
  5. Отличаются ли профили near-collision хешей?
"""

import numpy as np
import hashlib
import struct

MASK32 = 0xFFFFFFFF

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
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
def hw(x): return bin(x).count('1')

def add32_carry_count(x, y):
    result = (x + y) & MASK32
    c = 0
    count = 0
    for i in range(32):
        s = ((x >> i) & 1) + ((y >> i) & 1) + c
        c = s >> 1
        count += c
    return result, count

def add32(x, y): return (x + y) & MASK32


def get_carry_profile(W_input):
    """
    Полный carry-профиль SHA-256: для каждого ADD — число carry-бит.
    Возвращает вектор длины 448 (64 раунда × 7 ADD) + schedule.
    """
    W = list(W_input[:16])

    # Schedule carries (48 шагов × 3 ADD = 144 компонентов)
    sched_profile = []
    for r in range(16, 64):
        s0 = sigma0(W[r-15])
        s1 = sigma1(W[r-2])
        t1, c1 = add32_carry_count(s1, W[r-7])
        t2, c2 = add32_carry_count(t1, s0)
        t3, c3 = add32_carry_count(t2, W[r-16])
        sched_profile.extend([c1, c2, c3])
        W.append(t3)

    # Round carries (64 раунда × 7 ADD = 448 компонентов)
    a, b, c, d, e, f, g, h = IV
    round_profile = []

    for r in range(64):
        S1 = Sigma1(e)
        ch = Ch(e, f, g)
        S0 = Sigma0(a)
        maj = Maj(a, b, c)

        # T1 = h + S1 + ch + K[r] + W[r] (4 ADD)
        t1, c1 = add32_carry_count(h, S1)
        t2, c2 = add32_carry_count(t1, ch)
        t3, c3 = add32_carry_count(t2, K[r])
        T1, c4 = add32_carry_count(t3, W[r])

        # T2 = S0 + maj (1 ADD)
        T2, c5 = add32_carry_count(S0, maj)

        # e' = d + T1 (1 ADD)
        e_new, c6 = add32_carry_count(d, T1)

        # a' = T1 + T2 (1 ADD)
        a_new, c7 = add32_carry_count(T1, T2)

        round_profile.extend([c1, c2, c3, c4, c5, c6, c7])

        h, g, f, e = g, f, e, e_new
        d, c, b, a = c, b, a, a_new

    return np.array(round_profile), np.array(sched_profile)


def experiment_carry_profile_space(num_samples=2000):
    np.random.seed(42)

    print("=" * 70)
    print("CARRY-PROFILE SPACE — Физика carry-профилей SHA-256")
    print("=" * 70)

    # === 1. Собираем профили ===
    print("\n  Собираем carry-профили...")

    round_profiles = []
    sched_profiles = []

    for trial in range(num_samples):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        rp, sp = get_carry_profile(W)
        round_profiles.append(rp)
        sched_profiles.append(sp)

    R = np.array(round_profiles)  # (num_samples, 448)
    S = np.array(sched_profiles)  # (num_samples, 144)

    print(f"  Собрано {num_samples} профилей")
    print(f"  Round profile: {R.shape[1]} компонентов")
    print(f"  Schedule profile: {S.shape[1]} компонентов")

    # === 2. Статистика компонентов ===
    print("\n" + "=" * 70)
    print("1. СТАТИСТИКА CARRY-ПРОФИЛЯ")
    print("=" * 70)

    mean_R = np.mean(R, axis=0)
    std_R = np.std(R, axis=0)

    # Средний carry по типам ADD в раунде
    op_names = ['T1:h+S1', 'T1:+Ch', 'T1:+K', 'T1:+W', 'T2:S0+Maj', 'e:d+T1', 'a:T1+T2']
    print(f"\n  Средний carry по типам операций:")
    for i, name in enumerate(op_names):
        vals = R[:, i::7]  # каждый 7-й элемент
        print(f"    {name:12s}: mean={np.mean(vals):5.2f}, std={np.mean(np.std(vals, axis=0)):5.2f}")

    # Carry по раундам (суммарный)
    round_sums = np.array([np.sum(R[:, r*7:(r+1)*7], axis=1) for r in range(64)])  # (64, num_samples)
    mean_round = np.mean(round_sums, axis=1)
    std_round = np.std(round_sums, axis=1)

    print(f"\n  Суммарный carry по раундам:")
    for r in [0, 1, 2, 5, 10, 20, 32, 48, 63]:
        print(f"    r={r:2d}: {mean_round[r]:5.1f} ± {std_round[r]:4.1f}")

    # === 3. Дисперсия — какие компоненты детерминированы ===
    print("\n" + "=" * 70)
    print("2. ДЕТЕРМИНИРОВАННОСТЬ: ГДЕ ДИСПЕРСИЯ МИНИМАЛЬНА?")
    print("=" * 70)

    # Компоненты с наименьшей дисперсией = наиболее детерминированы
    low_var_threshold = 1.0
    determ_count = np.sum(std_R < low_var_threshold)
    print(f"\n  Компоненты с std < {low_var_threshold}: {determ_count} из {R.shape[1]} ({determ_count/R.shape[1]*100:.1f}%)")

    # По раундам
    print(f"\n  Детерминированность по раундам (доля компонентов с std < 3.0):")
    for r in range(0, 64, 4):
        round_stds = std_R[r*7:(r+1)*7]
        determ = np.sum(round_stds < 3.0)
        bar = "█" * (determ * 4)
        print(f"    r={r:2d}: {determ}/7 детерм.  {bar}")

    # === 4. PCA — истинная размерность ===
    print("\n" + "=" * 70)
    print("3. PCA: ИСТИННАЯ РАЗМЕРНОСТЬ ПРОСТРАНСТВА ПРОФИЛЕЙ")
    print("=" * 70)

    # Центрируем
    R_centered = R - np.mean(R, axis=0)

    # Ковариационная матрица (ограниченная для скорости)
    cov = np.cov(R_centered.T)

    # Собственные значения
    eigenvalues = np.linalg.eigvalsh(cov)
    eigenvalues = np.sort(eigenvalues)[::-1]

    # Сколько компонент нужно для 90%, 95%, 99% дисперсии
    total_var = np.sum(eigenvalues)
    cum_var = np.cumsum(eigenvalues) / total_var

    dim_90 = np.searchsorted(cum_var, 0.90) + 1
    dim_95 = np.searchsorted(cum_var, 0.95) + 1
    dim_99 = np.searchsorted(cum_var, 0.99) + 1

    print(f"\n  Всего компонентов: {R.shape[1]}")
    print(f"  Размерность для 90% дисперсии: {dim_90}")
    print(f"  Размерность для 95% дисперсии: {dim_95}")
    print(f"  Размерность для 99% дисперсии: {dim_99}")
    print(f"\n  Первые 20 собственных значений:")
    for i in range(min(20, len(eigenvalues))):
        pct = eigenvalues[i] / total_var * 100
        bar = "█" * int(pct * 2)
        print(f"    λ_{i:2d} = {eigenvalues[i]:8.2f} ({pct:5.2f}%)  {bar}")

    effective_dim = np.sum(eigenvalues > 0.01 * eigenvalues[0])
    print(f"\n  Эффективная размерность (λ > 1% от λ_max): {effective_dim}")

    # === 5. Дифференциальный профиль ===
    print("\n" + "=" * 70)
    print("4. ДИФФЕРЕНЦИАЛЬНЫЙ CARRY-ПРОФИЛЬ: δW → δProfile")
    print("=" * 70)

    diff_profiles = []
    for trial in range(500):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W1)
        W2[0] ^= 0x80000000

        rp1, _ = get_carry_profile(W1)
        rp2, _ = get_carry_profile(W2)
        diff_profiles.append(rp2 - rp1)

    D = np.array(diff_profiles)

    print(f"\n  δProfile при δW₀ = 1 бит MSB:")
    print(f"    Средний |δ| по компонентам: {np.mean(np.abs(D)):.3f}")
    print(f"    Макс. |δ| по компонентам:   {np.max(np.abs(np.mean(D, axis=0))):.3f}")

    # Какие раунды больше всего меняются?
    diff_by_round = np.zeros(64)
    for r in range(64):
        diff_by_round[r] = np.mean(np.abs(D[:, r*7:(r+1)*7]))

    print(f"\n  |δProfile| по раундам:")
    for r in [0, 1, 2, 3, 4, 5, 8, 16, 32, 48, 63]:
        bar = "█" * int(diff_by_round[r] * 5)
        print(f"    r={r:2d}: {diff_by_round[r]:5.2f}  {bar}")

    # PCA дифференциального профиля
    D_centered = D - np.mean(D, axis=0)
    cov_D = np.cov(D_centered.T)
    eig_D = np.sort(np.linalg.eigvalsh(cov_D))[::-1]
    total_var_D = np.sum(eig_D[eig_D > 0])
    cum_var_D = np.cumsum(eig_D[eig_D > 0]) / total_var_D if total_var_D > 0 else [0]

    dim_90_D = np.searchsorted(cum_var_D, 0.90) + 1
    dim_95_D = np.searchsorted(cum_var_D, 0.95) + 1

    print(f"\n  PCA дифференциального профиля:")
    print(f"    Размерность для 90%: {dim_90_D}")
    print(f"    Размерность для 95%: {dim_95_D}")

    # === 6. Near-collision: особый профиль? ===
    print("\n" + "=" * 70)
    print("5. NEAR-COLLISION: ОСОБЫЙ CARRY-ПРОФИЛЬ?")
    print("=" * 70)

    # Ищем пары с малым HW(δH)
    best_hw = 256
    best_profiles = None
    near_collision_profiles = []
    random_profiles = []

    for trial in range(1000):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W1)
        W2[0] ^= (1 << np.random.randint(0, 32))

        h1 = hashlib.sha256(struct.pack('>16I', *W1)).digest()
        h2 = hashlib.sha256(struct.pack('>16I', *W2)).digest()
        delta_h = bytes(a ^ b for a, b in zip(h1, h2))
        hw_h = sum(hw(b) for b in delta_h)

        rp1, _ = get_carry_profile(W1)
        rp2, _ = get_carry_profile(W2)
        dp = rp2 - rp1

        if hw_h < 100:
            near_collision_profiles.append(dp)
        else:
            random_profiles.append(dp)

        if hw_h < best_hw:
            best_hw = hw_h
            best_profiles = (rp1, rp2, dp, hw_h)

    NC = np.array(near_collision_profiles) if near_collision_profiles else np.zeros((1, 448))
    RN = np.array(random_profiles[:len(near_collision_profiles)]) if random_profiles else np.zeros((1, 448))

    print(f"\n  Near-collision (HW(δH) < 100): {len(near_collision_profiles)} пар")
    print(f"  Лучший HW(δH) найденный: {best_hw}")

    if len(near_collision_profiles) > 5 and len(random_profiles) > 5:
        mean_nc = np.mean(np.abs(NC), axis=0)
        mean_rn = np.mean(np.abs(RN), axis=0)

        # Сравнение по раундам
        print(f"\n  |δProfile| средний по раундам:")
        print(f"  {'Раунд':>6} {'Near-coll':>12} {'Random':>12} {'Ratio':>8}")
        for r in [0, 1, 2, 3, 5, 10, 20, 32, 48, 60, 63]:
            nc_round = np.mean(mean_nc[r*7:(r+1)*7])
            rn_round = np.mean(mean_rn[r*7:(r+1)*7])
            ratio = nc_round / max(rn_round, 0.001)
            marker = " ←" if ratio < 0.8 else ""
            print(f"  r={r:3d}  {nc_round:11.3f}  {rn_round:11.3f}  {ratio:7.3f}{marker}")

        # Суммарная энергия δProfile
        energy_nc = np.mean(np.sum(NC**2, axis=1))
        energy_rn = np.mean(np.sum(RN**2, axis=1))
        print(f"\n  'Энергия' δProfile (сумма квадратов):")
        print(f"    Near-collision: {energy_nc:.1f}")
        print(f"    Random:         {energy_rn:.1f}")
        print(f"    Ratio:          {energy_nc/max(energy_rn,1):.3f}")
    else:
        print("  Недостаточно near-collision пар для статистики")

    # === 7. Корреляция carry-профиля с выходом ===
    print("\n" + "=" * 70)
    print("6. CARRY-ПРОФИЛЬ → HASH: ЧТО ПРЕДСКАЗЫВАЕТ ЧТО?")
    print("=" * 70)

    # Корреляция компонентов профиля с битами хеша
    profiles_for_corr = []
    hash_hw_for_corr = []
    hash_lz_for_corr = []

    for trial in range(num_samples):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        rp, _ = get_carry_profile(W)
        h = hashlib.sha256(struct.pack('>16I', *W)).digest()
        hash_hw = sum(hw(b) for b in h)
        hash_lz = 0
        n = int.from_bytes(h, 'big')
        if n > 0:
            hash_lz = 256 - n.bit_length()

        profiles_for_corr.append(rp)
        hash_hw_for_corr.append(hash_hw)
        hash_lz_for_corr.append(hash_lz)

    P = np.array(profiles_for_corr)
    HW_arr = np.array(hash_hw_for_corr)

    # Корреляция каждой компоненты профиля с HW(hash)
    correlations = []
    for j in range(P.shape[1]):
        if np.std(P[:, j]) > 0:
            corr = np.corrcoef(P[:, j], HW_arr)[0, 1]
            correlations.append((j, corr))

    correlations.sort(key=lambda x: abs(x[1]), reverse=True)

    print(f"\n  Топ-10 компонентов профиля, коррелирующих с HW(hash):")
    for idx, (j, corr) in enumerate(correlations[:10]):
        r = j // 7
        op = j % 7
        print(f"    #{idx+1}: раунд {r:2d}, {op_names[op]:12s} — corr = {corr:+.4f}")

    max_corr = abs(correlations[0][1]) if correlations else 0
    print(f"\n  Макс. |корреляция|: {max_corr:.4f}")
    print(f"  → {'СТРУКТУРА!' if max_corr > 0.05 else 'Нет предсказательной силы'}")

    # Корреляция последних раундов с hash
    last_round_corrs = []
    for r in range(60, 64):
        for op in range(7):
            j = r * 7 + op
            if np.std(P[:, j]) > 0:
                corr = np.corrcoef(P[:, j], HW_arr)[0, 1]
                last_round_corrs.append((r, op_names[op], corr))

    print(f"\n  Корреляция carry (раунды 60-63) с HW(hash):")
    for r, name, corr in sorted(last_round_corrs, key=lambda x: abs(x[2]), reverse=True)[:10]:
        print(f"    r={r:2d} {name:12s}: corr = {corr:+.5f}")

    # === 8. Финальная сводка ===
    print("\n" + "=" * 70)
    print("СВОДКА: CARRY-PROFILE SPACE")
    print("=" * 70)

    print(f"""
  ПРОСТРАНСТВО:
    Номинальная размерность:  {R.shape[1]} компонентов (+ {S.shape[1]} schedule)
    PCA 90% дисперсии:        {dim_90} компонент
    PCA 95% дисперсии:        {dim_95} компонент
    PCA 99% дисперсии:        {dim_99} компонент
    Эффективная размерность:  {effective_dim}

  ДЕТЕРМИНИРОВАННОСТЬ:
    Компоненты с std < 1.0:   {determ_count} из {R.shape[1]}

  ДИФФЕРЕНЦИАЛЬНЫЙ:
    PCA δProfile 90%:         {dim_90_D} компонент
    PCA δProfile 95%:         {dim_95_D} компонент

  ПРЕДСКАЗАТЕЛЬНАЯ СИЛА:
    Макс. |corr(profile, HW(hash))|: {max_corr:.4f}
""")


if __name__ == "__main__":
    experiment_carry_profile_space(num_samples=2000)
