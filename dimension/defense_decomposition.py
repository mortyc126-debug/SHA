"""
Декомпозиция защитного механизма SHA-256.

Убираем компоненты по одному → измеряем 5 стен.
Цель: найти МИНИМАЛЬНОЕ ЯДРО, которое создаёт неуязвимость.

Варианты:
  A. Full SHA-256 (baseline)
  B. Без carry (ADD → XOR)
  C. Без Ch/Maj (заменяем на XOR)
  D. Без ротаций (Σ₀=Σ₁=σ₀=σ₁=identity)
  E. Без schedule expansion (W[16..63] = 0)
  F. Только 16 раундов
  G. Только 8 раундов

Для каждого варианта измеряем:
  Стена 1: Carry-нелинейность (HW(δH) distribution)
  Стена 2: Диффузия (раунды до 90% HW)
  Стена 3: Информация (corr(W[0], state[r]))
  Стена 4: ЦПТ-смешивание (std HW(δH))
  Стена 5: XOR-ADD несовместимость (corr XOR-δ vs ADD-δ)
"""

import numpy as np
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

def hw(x): return bin(x).count('1')
def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def shr(x, n): return x >> n


class SHA256Variant:
    """Конфигурируемый SHA-256 с возможностью отключения компонентов."""

    def __init__(self, name, use_carry=True, use_ch_maj=True,
                 use_rotations=True, use_schedule=True, num_rounds=64):
        self.name = name
        self.use_carry = use_carry
        self.use_ch_maj = use_ch_maj
        self.use_rotations = use_rotations
        self.use_schedule = use_schedule
        self.num_rounds = num_rounds

    def add(self, x, y):
        if self.use_carry:
            return (x + y) & MASK32
        else:
            return x ^ y  # XOR вместо ADD

    def sigma0(self, x):
        if self.use_rotations:
            return rotr(x, 7) ^ rotr(x, 18) ^ shr(x, 3)
        return x

    def sigma1(self, x):
        if self.use_rotations:
            return rotr(x, 17) ^ rotr(x, 19) ^ shr(x, 10)
        return x

    def Sigma0(self, x):
        if self.use_rotations:
            return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
        return x

    def Sigma1(self, x):
        if self.use_rotations:
            return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
        return x

    def Ch(self, e, f, g):
        if self.use_ch_maj:
            return (e & f) ^ (~e & g) & MASK32
        return e ^ f ^ g  # линейная замена

    def Maj(self, a, b, c):
        if self.use_ch_maj:
            return (a & b) ^ (a & c) ^ (b & c)
        return a ^ b ^ c  # линейная замена

    def compute(self, W_input):
        W = list(W_input[:16])
        if self.use_schedule:
            for r in range(16, 64):
                W.append(self.add(self.add(self.add(
                    self.sigma1(W[r-2]), W[r-7]),
                    self.sigma0(W[r-15])), W[r-16]))
        else:
            W.extend([0] * 48)

        a, b, c, d, e, f, g, h = IV
        states = [(a, b, c, d, e, f, g, h)]

        for r in range(self.num_rounds):
            S1 = self.Sigma1(e)
            ch = self.Ch(e, f, g)
            S0 = self.Sigma0(a)
            maj = self.Maj(a, b, c)
            T1 = self.add(self.add(self.add(self.add(h, S1), ch), K[r % 64]), W[r % len(W)])
            T2 = self.add(S0, maj)
            h, g, f, e = g, f, e, self.add(d, T1)
            d, c, b, a = c, b, a, self.add(T1, T2)
            states.append((a, b, c, d, e, f, g, h))

        H = [self.add(IV[i], states[-1][i]) for i in range(8)]
        return H, states


def measure_walls(variant, num_samples=1000):
    """Измеряем 5 стен для данного варианта."""
    np.random.seed(42)

    # === Стена 1: Диффузия (HW(δH) при 1 бит δW) ===
    dh_hws = []
    for _ in range(num_samples):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W1)
        W2[0] ^= 0x80000000
        H1, _ = variant.compute(W1)
        H2, _ = variant.compute(W2)
        dh = sum(hw(H1[i] ^ H2[i]) for i in range(8))
        dh_hws.append(dh)

    wall1_mean_hw = np.mean(dh_hws)
    wall1_std_hw = np.std(dh_hws)
    # Идеальный random oracle: mean=128, std≈8
    wall1_score = min(wall1_mean_hw / 128.0, 1.0)

    # === Стена 2: Скорость диффузии (раундов до 90%) ===
    target_hw = 0.9 * (min(variant.num_rounds, 64) / 64.0 * 128)
    if target_hw < 10:
        target_hw = 10
    round_90 = variant.num_rounds

    for test_rounds in range(1, variant.num_rounds + 1):
        test_var = SHA256Variant(
            f"test_{test_rounds}",
            use_carry=variant.use_carry,
            use_ch_maj=variant.use_ch_maj,
            use_rotations=variant.use_rotations,
            use_schedule=variant.use_schedule,
            num_rounds=test_rounds
        )
        hws = []
        for _ in range(200):
            W1 = [np.random.randint(0, 2**32) for _ in range(16)]
            W2 = list(W1)
            W2[0] ^= 0x80000000
            H1, _ = test_var.compute(W1)
            H2, _ = test_var.compute(W2)
            # Используем state, не H (чтобы видеть промежуточные)
            dh = sum(hw(H1[i] ^ H2[i]) for i in range(8))
            hws.append(dh)
        if np.mean(hws) >= target_hw:
            round_90 = test_rounds
            break

    wall2_rounds = round_90
    wall2_score = 1.0 - (wall2_rounds / variant.num_rounds)

    # === Стена 3: Информационная (corr W[0] → H[0]) ===
    w0_lsb = []
    h0_lsb = []
    for _ in range(num_samples):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H, _ = variant.compute(W)
        w0_lsb.append(W[0] & 0xFF)
        h0_lsb.append(H[0] & 0xFF)

    wall3_corr = abs(np.corrcoef(w0_lsb, h0_lsb)[0, 1])
    wall3_score = 1.0 - min(wall3_corr * 10, 1.0)  # 0 = плохо, 1 = хорошо

    # === Стена 4: ЦПТ (нормальность распределения HW(δH)) ===
    # Для random oracle: std(HW) ≈ 8. Если меньше — структура, если больше — тоже
    wall4_std = wall1_std_hw
    wall4_score = min(wall4_std / 8.0, 1.0)

    # === Стена 5: XOR-ADD несовместимость ===
    if variant.use_carry:
        xor_diffs = []
        add_diffs = []
        for _ in range(num_samples):
            W1 = [np.random.randint(0, 2**32) for _ in range(16)]
            W2 = list(W1)
            delta = np.random.randint(1, 2**32)
            W2[0] ^= delta
            H1, _ = variant.compute(W1)
            H2, _ = variant.compute(W2)
            xor_d = sum(hw(H1[i] ^ H2[i]) for i in range(8))
            add_d = hw((H1[0] - H2[0]) % (2**32))
            xor_diffs.append(xor_d)
            add_diffs.append(add_d)
        wall5_corr = abs(np.corrcoef(xor_diffs, add_diffs)[0, 1])
        wall5_score = 1.0 - min(wall5_corr * 5, 1.0)
    else:
        wall5_corr = 1.0  # без carry XOR и ADD — одно и то же
        wall5_score = 0.0

    # === Бонус: можно ли найти collision? ===
    # Пробуем birthday-подобный поиск на малом множестве
    collision_found = False
    hash_dict = {}
    for trial in range(min(num_samples * 5, 10000)):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H, _ = variant.compute(W)
        h_tuple = tuple(H)
        if h_tuple in hash_dict:
            if hash_dict[h_tuple] != tuple(W):
                collision_found = True
                break
        hash_dict[h_tuple] = tuple(W)

    return {
        'wall1_hw': wall1_mean_hw,
        'wall1_std': wall1_std_hw,
        'wall1_score': wall1_score,
        'wall2_rounds': wall2_rounds,
        'wall2_score': wall2_score,
        'wall3_corr': wall3_corr,
        'wall3_score': wall3_score,
        'wall4_std': wall4_std,
        'wall4_score': wall4_score,
        'wall5_corr': wall5_corr,
        'wall5_score': wall5_score,
        'collision': collision_found,
        'total_score': (wall1_score + wall2_score + wall3_score + wall4_score + wall5_score) / 5.0
    }


def main():
    print("=" * 78)
    print("ДЕКОМПОЗИЦИЯ ЗАЩИТНОГО МЕХАНИЗМА SHA-256")
    print("Убираем компоненты → смотрим какие стены рушатся")
    print("=" * 78)

    variants = [
        SHA256Variant("A: Full SHA-256",           True,  True,  True,  True,  64),
        SHA256Variant("B: Без carry (XOR)",         False, True,  True,  True,  64),
        SHA256Variant("C: Без Ch/Maj (линейн.)",    True,  False, True,  True,  64),
        SHA256Variant("D: Без ротаций",             True,  True,  False, True,  64),
        SHA256Variant("E: Без schedule",            True,  True,  True,  False, 64),
        SHA256Variant("F: 16 раундов",              True,  True,  True,  True,  16),
        SHA256Variant("G: 8 раундов",               True,  True,  True,  True,  8),
        SHA256Variant("H: Без carry+Ch/Maj",        False, False, True,  True,  64),
        SHA256Variant("I: Без carry+ротаций",       False, True,  False, True,  64),
        SHA256Variant("J: 8р без schedule",         True,  True,  True,  False, 8),
    ]

    results = {}
    for v in variants:
        print(f"\n  Анализ: {v.name}...", end="", flush=True)
        r = measure_walls(v, num_samples=500)
        results[v.name] = r
        print(f" done (score={r['total_score']:.3f})")

    # === Таблица результатов ===
    print("\n" + "=" * 78)
    print("РЕЗУЛЬТАТЫ: 5 СТЕН × 10 ВАРИАНТОВ")
    print("=" * 78)

    print(f"\n{'Вариант':<28} {'HW(δH)':>7} {'Diff.r':>6} {'Corr':>6} {'Std':>5} {'XA':>5} {'Score':>6} {'Coll':>5}")
    print(f"{'':<28} {'(128)':>7} {'(≤5)':>6} {'(~0)':>6} {'(~8)':>5} {'(~0)':>5} {'(1.0)':>6} {'':>5}")
    print("-" * 78)

    for name, r in results.items():
        coll_str = "YES!" if r['collision'] else "no"
        print(f"{name:<28} {r['wall1_hw']:>6.1f} {r['wall2_rounds']:>5d}  "
              f"{r['wall3_corr']:>5.3f} {r['wall4_std']:>5.1f} {r['wall5_corr']:>4.2f} "
              f"{r['total_score']:>6.3f} {coll_str:>5}")

    # === Анализ: какие компоненты КРИТИЧНЫ ===
    print("\n" + "=" * 78)
    print("АНАЛИЗ: КАКИЕ КОМПОНЕНТЫ КРИТИЧНЫ?")
    print("=" * 78)

    baseline = results["A: Full SHA-256"]

    print(f"\n  {'Компонент убран':<28} {'ΔScore':>8} {'ΔHW(δH)':>9} {'ΔCorr':>8} {'Collision':>10}")
    print("-" * 68)

    for name, r in results.items():
        if name == "A: Full SHA-256":
            continue
        ds = r['total_score'] - baseline['total_score']
        dhw = r['wall1_hw'] - baseline['wall1_hw']
        dcorr = r['wall3_corr'] - baseline['wall3_corr']
        coll = "COLLISION!" if r['collision'] else "safe"
        marker = " ← CRITICAL" if ds < -0.1 or r['collision'] else ""
        print(f"  {name:<28} {ds:>+7.3f} {dhw:>+8.1f} {dcorr:>+7.4f}  {coll:<10}{marker}")

    # === Какие стены падают при каждом удалении ===
    print("\n" + "=" * 78)
    print("МАТРИЦА: СТЕНА × КОМПОНЕНТ (что ломается)")
    print("=" * 78)

    wall_names = ['Диффузия', 'Скорость', 'Информ.', 'ЦПТ', 'XOR-ADD']
    wall_keys = ['wall1_score', 'wall2_score', 'wall3_score', 'wall4_score', 'wall5_score']

    header = f"{'Вариант':<28}"
    for wn in wall_names:
        header += f" {wn:>9}"
    print(header)
    print("-" * 78)

    for name, r in results.items():
        row = f"{name:<28}"
        for wk in wall_keys:
            val = r[wk]
            if val >= 0.8:
                symbol = "█████"  # стена стоит
            elif val >= 0.5:
                symbol = "███░░"  # повреждена
            elif val >= 0.2:
                symbol = "█░░░░"  # почти разрушена
            else:
                symbol = "░░░░░"  # разрушена
            row += f" {symbol:>9}"
        row += f"  {'COLL!' if r['collision'] else ''}"
        print(row)

    print(f"\n  █████ = стена стоит  ███░░ = повреждена  █░░░░ = критична  ░░░░░ = разрушена")

    # === Финальный вывод ===
    print("\n" + "=" * 78)
    print("МЕХАНИЗМ ЗАЩИТЫ SHA-256 — ФИНАЛЬНАЯ КАРТА")
    print("=" * 78)

    # Определяем критичные компоненты
    critical = []
    for name, r in results.items():
        if name == "A: Full SHA-256":
            continue
        if r['total_score'] < baseline['total_score'] - 0.1 or r['collision']:
            critical.append(name)

    non_critical = []
    for name, r in results.items():
        if name == "A: Full SHA-256":
            continue
        if r['total_score'] >= baseline['total_score'] - 0.1 and not r['collision']:
            non_critical.append(name)

    print(f"\n  КРИТИЧНЫЕ компоненты (удаление ломает защиту):")
    for c in critical:
        print(f"    × {c}")

    print(f"\n  УСТОЙЧИВЫЕ компоненты (удаление НЕ ломает):")
    for c in non_critical:
        print(f"    ✓ {c}")

    print(f"""
  МИНИМАЛЬНОЕ ЯДРО ЗАЩИТЫ:
  Компоненты, без которых SHA-256 ЛОМАЕТСЯ — это и есть механизм.
  Компоненты, которые можно убрать без потери — избыточная защита.

  Чтобы пробить SHA-256, нужно одновременно преодолеть
  ВСЕ критичные компоненты. Каждый из них создаёт свой
  независимый барьер.
""")


if __name__ == "__main__":
    main()
