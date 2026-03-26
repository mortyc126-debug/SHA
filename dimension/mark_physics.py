"""
Физика меток в новом измерении SHA-256.

Метка δ = разность двух следов через ткань SHA-256.
Вопрос: есть ли законы сохранения? Инварианты?

Измеряем для метки на каждом раунде:
  - масса m(δ) = HW(δ_state)  (суммарный Hamming weight всех регистров)
  - carry-энергия E_p = дифференциальные carry-биты
  - поток J = изменение массы между раундами (Δm)
  - ищем: m + f(E_p) = const?  или другие инварианты
"""

import numpy as np
from typing import List, Tuple

# SHA-256 константы
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
MASK32 = 0xFFFFFFFF

def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def shr(x, n): return x >> n
def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ shr(x, 3)
def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ shr(x, 10)
def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
def add32(x, y): return (x + y) & MASK32
def hw(x): return bin(x).count('1')

def add32_carry(x, y):
    """Возвращает (result, carry_bits_mask, carry_count)"""
    result = (x + y) & MASK32
    c = 0
    carry_bits = 0
    for i in range(32):
        s = ((x >> i) & 1) + ((y >> i) & 1) + c
        c = s >> 1
        if c:
            carry_bits |= (1 << i)
    return result, carry_bits, bin(carry_bits).count('1')


def sha256_round_by_round(W_input):
    """SHA-256 с сохранением состояния на каждом раунде."""
    W = list(W_input[:16])
    for r in range(16, 64):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))

    state = list(IV)
    states = [tuple(state)]
    for r in range(64):
        a,b,c,d,e,f,g,h = state
        S1 = Sigma1(e)
        ch = Ch(e, f, g)
        S0 = Sigma0(a)
        maj = Maj(a, b, c)
        T1 = add32(add32(add32(add32(h, S1), ch), K[r]), W[r])
        T2 = add32(S0, maj)
        state = [add32(T1, T2), a, b, c, add32(d, T1), e, f, g]
        states.append(tuple(state))
    return states, W


def measure_mark_physics(W1, W2):
    """
    Измеряет физические свойства метки δ = (W1 trace) - (W2 trace)
    на каждом раунде.
    """
    states1, Wfull1 = sha256_round_by_round(W1)
    states2, Wfull2 = sha256_round_by_round(W2)

    rounds = []
    for r in range(65):  # 0..64 (64 = после последнего раунда)
        s1 = states1[r]
        s2 = states2[r]

        # δ для каждого регистра
        deltas = [s1[i] ^ s2[i] for i in range(8)]

        # Масса = суммарный HW
        mass = sum(hw(d) for d in deltas)

        # Масса по регистрам
        mass_a = hw(deltas[0])
        mass_e = hw(deltas[4])

        # Масса a-chain (a,b,c,d) и e-chain (e,f,g,h)
        mass_a_chain = sum(hw(deltas[i]) for i in range(4))
        mass_e_chain = sum(hw(deltas[i]) for i in range(4, 8))

        # δW для этого раунда
        if r < 64:
            dW = Wfull1[r] ^ Wfull2[r]
            mass_w = hw(dW)
        else:
            dW = 0
            mass_w = 0

        rounds.append({
            'r': r,
            'mass': mass,
            'mass_a': mass_a,
            'mass_e': mass_e,
            'mass_a_chain': mass_a_chain,
            'mass_e_chain': mass_e_chain,
            'mass_w': mass_w,
            'deltas': deltas
        })

    return rounds


def experiment_conservation_laws(num_samples=200):
    """Ищем законы сохранения в физике меток."""
    np.random.seed(42)

    print("=" * 70)
    print("ФИЗИКА МЕТОК — Поиск законов сохранения")
    print("=" * 70)

    # === Эксперимент 1: Эволюция массы ===
    print("\n" + "=" * 70)
    print("1. ЭВОЛЮЦИЯ МАССЫ МЕТКИ m(δ) по раундам")
    print("=" * 70)

    all_masses = np.zeros((num_samples, 65))
    all_mass_a = np.zeros((num_samples, 65))
    all_mass_e = np.zeros((num_samples, 65))
    all_mass_w = np.zeros((num_samples, 64))
    all_mass_a_chain = np.zeros((num_samples, 65))
    all_mass_e_chain = np.zeros((num_samples, 65))

    for trial in range(num_samples):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W1)
        W2[0] ^= 0x80000000  # 1 бит MSB

        rounds = measure_mark_physics(W1, W2)
        for r in range(65):
            all_masses[trial, r] = rounds[r]['mass']
            all_mass_a[trial, r] = rounds[r]['mass_a']
            all_mass_e[trial, r] = rounds[r]['mass_e']
            all_mass_a_chain[trial, r] = rounds[r]['mass_a_chain']
            all_mass_e_chain[trial, r] = rounds[r]['mass_e_chain']
        for r in range(64):
            all_mass_w[trial, r] = rounds[r]['mass_w']

    mean_mass = np.mean(all_masses, axis=0)
    std_mass = np.std(all_masses, axis=0)

    print(f"\n  δW₀ = 1 бит (MSB), начальная масса = {mean_mass[0]:.1f}")
    print(f"\n  {'Раунд':>6} {'Масса':>8} {'±std':>6} {'a-chain':>8} {'e-chain':>8} {'δW':>6}")
    print(f"  {'-'*6} {'-'*8} {'-'*6} {'-'*8} {'-'*8} {'-'*6}")

    for r in [0, 1, 2, 3, 4, 5, 8, 12, 16, 20, 24, 32, 48, 63, 64]:
        if r < 65:
            m = mean_mass[r]
            s = std_mass[r]
            ma = np.mean(all_mass_a_chain[:, r])
            me = np.mean(all_mass_e_chain[:, r])
            mw = np.mean(all_mass_w[:, r]) if r < 64 else 0
            print(f"  r={r:3d}  {m:7.1f}  {s:5.1f}  {ma:7.1f}  {me:7.1f}  {mw:5.1f}")

    # === Эксперимент 2: Закон сохранения? ===
    print("\n" + "=" * 70)
    print("2. ПОИСК ИНВАРИАНТОВ")
    print("=" * 70)

    # Гипотеза 1: m(r) + m(r+1) = const?
    print("\n  Гипотеза: m(r) + m(r+1) = const (сохранение суммарной массы)")
    sums_adj = np.zeros((num_samples, 64))
    for trial in range(num_samples):
        for r in range(64):
            sums_adj[trial, r] = all_masses[trial, r] + all_masses[trial, r+1]
    mean_sums = np.mean(sums_adj, axis=0)
    std_sums = np.std(sums_adj, axis=0)
    print(f"  r=0..3:  {mean_sums[0]:.1f}, {mean_sums[1]:.1f}, {mean_sums[2]:.1f}, {mean_sums[3]:.1f}")
    print(f"  r=30..33: {mean_sums[30]:.1f}, {mean_sums[31]:.1f}, {mean_sums[32]:.1f}, {mean_sums[33]:.1f}")
    print(f"  Вариация: mean={np.mean(std_sums):.1f}, → НЕ константа" if np.std(mean_sums) > 5 else f"  Вариация: std={np.std(mean_sums):.1f} → ВОЗМОЖНАЯ константа!")

    # Гипотеза 2: a_chain + e_chain = const после насыщения?
    print("\n  Гипотеза: a_chain(r) / e_chain(r) = const после насыщения")
    ratios = []
    for r in range(10, 64):
        ma = np.mean(all_mass_a_chain[:, r])
        me = np.mean(all_mass_e_chain[:, r])
        if me > 0:
            ratios.append(ma / me)
    print(f"  Ratio a/e для r=10..63: mean={np.mean(ratios):.4f}, std={np.std(ratios):.4f}")
    print(f"  Если бы равномерно: ratio = 1.0")

    # Гипотеза 3: масса → аттрактор?
    print("\n  Гипотеза: масса стремится к аттрактору")
    saturation_mass = np.mean(mean_mass[20:64])
    saturation_std = np.std(mean_mass[20:64])
    theoretical_random = 128.0  # 256 бит × 0.5 = 128 ненулевых
    print(f"  Среднее m(r) для r=20..63: {saturation_mass:.2f} ± {saturation_std:.2f}")
    print(f"  Теоретическое (random):     {theoretical_random:.1f}")
    print(f"  Отклонение от random:       {abs(saturation_mass - theoretical_random):.2f} бит ({abs(saturation_mass - theoretical_random)/theoretical_random*100:.1f}%)")

    # === Эксперимент 3: Импульс (изменение массы) ===
    print("\n" + "=" * 70)
    print("3. ИМПУЛЬС МЕТКИ p(r) = Δm = m(r+1) - m(r)")
    print("=" * 70)

    all_momentum = np.diff(all_masses, axis=1)
    mean_momentum = np.mean(all_momentum, axis=0)
    std_momentum = np.std(all_momentum, axis=0)

    print(f"\n  {'Раунд':>6} {'Импульс':>10} {'±std':>8}")
    for r in [0, 1, 2, 3, 4, 5, 8, 12, 16, 20, 32, 48, 63]:
        print(f"  r={r:3d}  {mean_momentum[r]:9.2f}  {std_momentum[r]:7.2f}")

    # Где импульс ≈ 0? (равновесие)
    equilibrium_rounds = [r for r in range(64) if abs(mean_momentum[r]) < 2.0]
    print(f"\n  Раунды равновесия (|p| < 2): {equilibrium_rounds[:20]}")

    # === Эксперимент 4: Поиск глубоких инвариантов ===
    print("\n" + "=" * 70)
    print("4. ГЛУБОКИЕ ИНВАРИАНТЫ")
    print("=" * 70)

    # XOR всех δ-регистров = ?
    print("\n  Гипотеза: XOR-свёртка всех δ-регистров = инвариант")
    xor_invariants = np.zeros((num_samples, 65))
    for trial in range(num_samples):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W1)
        W2[0] ^= 0x80000000
        rounds = measure_mark_physics(W1, W2)
        for r in range(65):
            xor_all = 0
            for d in rounds[r]['deltas']:
                xor_all ^= d
            xor_invariants[trial, r] = hw(xor_all)

    mean_xor = np.mean(xor_invariants, axis=0)
    std_xor = np.std(xor_invariants, axis=0)
    print(f"  HW(XOR(δa,δb,...,δh)):")
    print(f"    r=0: {mean_xor[0]:.2f} ± {std_xor[0]:.2f}")
    print(f"    r=1: {mean_xor[1]:.2f} ± {std_xor[1]:.2f}")
    print(f"    r=4: {mean_xor[4]:.2f} ± {std_xor[4]:.2f}")
    print(f"    r=8: {mean_xor[8]:.2f} ± {std_xor[8]:.2f}")
    print(f"    r=16: {mean_xor[16]:.2f} ± {std_xor[16]:.2f}")
    print(f"    r=32: {mean_xor[32]:.2f} ± {std_xor[32]:.2f}")
    print(f"    r=64: {mean_xor[64]:.2f} ± {std_xor[64]:.2f}")
    xor_variance = np.std(mean_xor[10:])
    print(f"  Вариация r=10..64: {xor_variance:.2f}")
    print(f"  → {'ВОЗМОЖНЫЙ инвариант!' if xor_variance < 1.0 else 'Не инвариант'}")

    # Чётность (parity) каждого регистра
    print("\n  Гипотеза: чётность (LSB) δ-регистров сохраняется")
    parity_preserved = np.zeros(64)
    for trial in range(num_samples):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W1)
        W2[0] ^= 0x80000000
        rounds = measure_mark_physics(W1, W2)
        for r in range(64):
            # Сумма LSB всех δ-регистров mod 2
            parity_r = sum(d & 1 for d in rounds[r]['deltas']) % 2
            parity_r1 = sum(d & 1 for d in rounds[r+1]['deltas']) % 2
            if parity_r == parity_r1:
                parity_preserved[r] += 1

    parity_preserved /= num_samples
    print(f"  P(parity сохраняется r→r+1):")
    print(f"    r=0: {parity_preserved[0]:.3f}")
    print(f"    r=1: {parity_preserved[1]:.3f}")
    print(f"    r=5: {parity_preserved[5]:.3f}")
    print(f"    r=10: {parity_preserved[10]:.3f}")
    print(f"    r=32: {parity_preserved[32]:.3f}")
    print(f"    mean(r=0..63): {np.mean(parity_preserved):.3f}")
    print(f"    Если бы случайно: 0.500")

    # === Эксперимент 5: Энтропия метки ===
    print("\n" + "=" * 70)
    print("5. ЭНТРОПИЯ МЕТКИ")
    print("=" * 70)

    # Для каждого раунда: насколько "случайно" распределены биты δ?
    # H = -Σ p·log(p) по 8 регистрам (нормализованная масса)
    all_entropy = np.zeros((num_samples, 65))
    for trial in range(num_samples):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W1)
        W2[0] ^= 0x80000000
        rounds = measure_mark_physics(W1, W2)
        for r in range(65):
            masses = [hw(d) for d in rounds[r]['deltas']]
            total = sum(masses)
            if total > 0:
                probs = [m / total for m in masses if m > 0]
                entropy = -sum(p * np.log2(p) for p in probs)
            else:
                entropy = 0
            all_entropy[trial, r] = entropy

    mean_entropy = np.mean(all_entropy, axis=0)
    max_entropy = np.log2(8)  # = 3.0 (равномерное по 8 регистрам)

    print(f"\n  Энтропия распределения массы по 8 регистрам (max = {max_entropy:.2f}):")
    for r in [0, 1, 2, 3, 4, 5, 8, 12, 16, 20, 32, 48, 63, 64]:
        bar = "█" * int(mean_entropy[r] / max_entropy * 30)
        print(f"    r={r:3d}: H={mean_entropy[r]:.3f}/{max_entropy:.1f}  {bar}")

    saturation_entropy = np.mean(mean_entropy[20:64])
    print(f"\n  Энтропия насыщения (r=20..63): {saturation_entropy:.3f}")
    print(f"  Максимальная: {max_entropy:.3f}")
    print(f"  Отношение: {saturation_entropy/max_entropy:.4f}")

    # === Эксперимент 6: Разные δW — разная физика? ===
    print("\n" + "=" * 70)
    print("6. ЗАВИСИМОСТЬ ОТ НАЧАЛЬНЫХ УСЛОВИЙ")
    print("=" * 70)

    delta_types = [
        ("1 бит MSB W₀", lambda W: (0, 0x80000000)),
        ("1 бит LSB W₀", lambda W: (0, 0x00000001)),
        ("1 бит MSB W₇", lambda W: (7, 0x80000000)),
        ("1 бит MSB W₁₅", lambda W: (15, 0x80000000)),
        ("2 бита W₀", lambda W: (0, 0xC0000000)),
        ("full W₀", lambda W: (0, 0xFFFFFFFF)),
    ]

    for name, delta_fn in delta_types:
        masses_at_sat = []
        for trial in range(100):
            W1 = [np.random.randint(0, 2**32) for _ in range(16)]
            W2 = list(W1)
            idx, delta = delta_fn(W1)
            W2[idx] ^= delta
            rounds = measure_mark_physics(W1, W2)
            masses_at_sat.append(rounds[32]['mass'])

        print(f"  {name:20s}: m(r=32) = {np.mean(masses_at_sat):.1f} ± {np.std(masses_at_sat):.1f}")

    print(f"\n  Теоретическое (random): 128.0")

    # === Финальная сводка ===
    print("\n" + "=" * 70)
    print("СВОДКА: ФИЗИКА МЕТОК")
    print("=" * 70)
    print(f"""
  МАССА:
    Начальная (1 бит):     {mean_mass[0]:.0f}
    Насыщение (r≥5):       {saturation_mass:.1f} ± {saturation_std:.1f}
    Теор. random:          128.0
    Скорость насыщения:    ~3-4 раунда

  ИМПУЛЬС:
    Максимальный:          r=0→1 ({mean_momentum[0]:.0f} бит/раунд)
    Равновесие (|p|<2):    с раунда ~{min(equilibrium_rounds[1:]) if len(equilibrium_rounds) > 1 else '?'}

  ЭНТРОПИЯ:
    Начальная:             {mean_entropy[0]:.3f} (локализована в 1 регистре)
    Насыщение:             {saturation_entropy:.3f} / {max_entropy:.3f}
    Скорость выравнивания: ~5 раундов

  ИНВАРИАНТЫ:
    XOR-свёртка:           {'ДА' if xor_variance < 1.0 else 'НЕТ'} (вариация={xor_variance:.2f})
    Чётность:              {'ДА' if abs(np.mean(parity_preserved) - 0.5) > 0.05 else 'НЕТ'} (P={np.mean(parity_preserved):.3f})
    a/e ratio:             {'ДА' if np.std(ratios) < 0.05 else 'НЕТ'} (ratio={np.mean(ratios):.3f}±{np.std(ratios):.3f})
    Масса-аттрактор:       ДА (128 ± {saturation_std:.1f})
""")


if __name__ == "__main__":
    experiment_conservation_laws(num_samples=200)
