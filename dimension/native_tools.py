"""
НАТИВНЫЕ ИНСТРУМЕНТЫ НАШЕГО ИЗМЕРЕНИЯ.

Отказываемся от:
  × Bits (единица старого мира)
  × Correlation (Pearson — чужой инструмент)
  × Birthday (слепой перебор — варварство)
  × HW (Hamming weight — битовая мера)

Создаём:
  ✓ СВЯЗНОСТЬ (coupling) — нативная мера связи через ткань
  ✓ МОДЫ — нативные единицы измерения (не биты)
  ✓ НАВИГАЦИЯ — направленное движение по ткани (не перебор)
  ✓ ЭНЕРГИЯ МЕТКИ — мера в нашем пространстве
"""

import numpy as np
import struct, hashlib

MASK32 = 0xFFFFFFFF

def sha256_words(W16):
    return struct.unpack('>8I', hashlib.sha256(struct.pack('>16I', *W16)).digest())

def hw(x): return bin(x).count('1')


# =======================================================================
# ОПРЕДЕЛЕНИЕ 1: СВЯЗНОСТЬ (Coupling)
#
# В старом мире: corr(X, Y) = Pearson correlation (статистическая)
# В нашем мире: coupling(pos_i, pos_j) = чистота пути между ними
#
# Coupling = 1: позиции соединены ТРУБОЙ (δ передаётся без потерь)
# Coupling = 0: позиции не связаны (δ не доходит)
# Coupling ∈ (0,1): связаны через УЗЛЫ (δ частично искажается)
# =======================================================================

def measure_coupling(pos_i, pos_j, n_samples=5000):
    """
    Coupling между позициями хеша:
    coupling(i,j) = P(δH[j] мал | δH[i]=0)

    Если i и j связаны трубами → coupling высокий.
    Если через узлы → coupling низкий.
    """
    np.random.seed(42 + pos_i * 8 + pos_j)

    # Ищем пары с δH[pos_i]=0 и смотрим δH[pos_j]
    h_dict = {}
    conditional_dh = []
    unconditional_dh = []

    for trial in range(n_samples * 100):  # нужно много для birthday на 32 бит
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H = sha256_words(W)
        key = H[pos_i]

        if key in h_dict:
            H_prev = h_dict[key]
            dh_j = hw(H[pos_j] ^ H_prev[pos_j])
            conditional_dh.append(dh_j)
            if len(conditional_dh) >= 100:
                break
        else:
            h_dict[key] = H

        # Unconditional: случайная пара
        if trial % 1000 == 0 and trial > 0:
            W2 = [np.random.randint(0, 2**32) for _ in range(16)]
            H2 = sha256_words(W2)
            unconditional_dh.append(hw(H[pos_j] ^ H2[pos_j]))

    if not conditional_dh:
        return 0.0, 16.0, 16.0

    mean_cond = np.mean(conditional_dh)
    mean_uncond = np.mean(unconditional_dh) if unconditional_dh else 16.0

    # Coupling = (16 - mean_conditional) / 16
    # = 0 если conditional = 16 (random, не связаны)
    # = 1 если conditional = 0 (полная связь)
    coupling = max(0, (16 - mean_cond) / 16)

    return coupling, mean_cond, mean_uncond


# =======================================================================
# ОПРЕДЕЛЕНИЕ 2: ЭНЕРГИЯ МЕТКИ
#
# В старом мире: HW(δH) — число отличающихся бит
# В нашем мире: E(δ) — взвешенная мера через coupling
#
# E(δ) = Σ coupling_weight[i] × |δH[i]|
# Позиции с высоким coupling к уже обнулённым → дешевле
# =======================================================================

def mark_energy(delta_H, coupling_matrix):
    """
    Энергия метки в нашем измерении.

    delta_H: список из 8 значений δH[i]
    coupling_matrix: 8×8 матрица связности

    Энергия учитывает: если позиция i сильно связана с уже
    нулевой позицией j, то δH[i] "дешевле" обнулить.
    """
    n = len(delta_H)
    energy = 0.0

    for i in range(n):
        if delta_H[i] == 0:
            continue

        # Вес позиции = 1 - max coupling с нулевыми позициями
        max_coupling = 0.0
        for j in range(n):
            if j != i and delta_H[j] == 0:
                max_coupling = max(max_coupling, coupling_matrix[i][j])

        weight = 1.0 - max_coupling
        energy += weight * hw(delta_H[i])

    return energy


# =======================================================================
# ОПРЕДЕЛЕНИЕ 3: НАВИГАЦИЯ
#
# В старом мире: birthday = слепой перебор → ждём совпадения
# В нашем мире: навигация = направленное движение по ткани
#
# Алгоритм навигации:
#   1. Выбираем начальную метку δ
#   2. Вычисляем энергию E(δ) через ткань
#   3. Модифицируем δ в направлении МИНИМАЛЬНОЙ энергии
#   4. Используем coupling-матрицу для выбора направления
# =======================================================================

def navigate_fabric(coupling_matrix, n_attempts=50000):
    """
    Навигация по ткани: ищем метку с минимальной энергией.

    Вместо слепого перебора:
      - Используем coupling для выбора направления
      - Двигаемся вдоль СВЯЗНЫХ позиций
      - Обнуляем сначала САМЫЕ СВЯЗНЫЕ слова
    """
    np.random.seed(42)

    best_energy = float('inf')
    best_hw = 256
    best_delta_W = None

    # Определяем порядок обнуления: сначала самые связные
    # (те, чья стоимость снижается при обнулении соседей)
    word_connectivity = [sum(coupling_matrix[i]) for i in range(8)]
    priority = np.argsort(word_connectivity)[::-1]  # от самого связного

    for attempt in range(n_attempts):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]

        # НАПРАВЛЕННЫЙ поиск: варьируем δW в направлении
        # обнуления приоритетных позиций

        # Стратегия: генерируем δW, ориентируясь на coupling
        if attempt % 3 == 0:
            # Случайный 1-бит δ
            delta = [0] * 16
            delta[0] = 1 << np.random.randint(0, 32)
        elif attempt % 3 == 1:
            # δ в нескольких словах (multi-word)
            delta = [0] * 16
            for w in range(np.random.randint(1, 4)):
                delta[np.random.randint(0, 16)] = np.random.randint(1, 2**32)
        else:
            # δ подобранный под coupling structure
            delta = [0] * 16
            # Варьируем слова, соответствующие высокоприоритетным H
            target_h = priority[0]  # самое связное H
            # H[target] зависит от W[target] через round target
            delta[min(target_h, 15)] = 1 << np.random.randint(0, 32)

        W2 = [W1[i] ^ delta[i] for i in range(16)]
        if all(d == 0 for d in delta):
            continue

        H1 = sha256_words(W1)
        H2 = sha256_words(W2)

        delta_H = [H1[i] ^ H2[i] for i in range(8)]

        # Энергия в нашем измерении
        energy = mark_energy(delta_H, coupling_matrix)
        total_hw = sum(hw(d) for d in delta_H)

        if energy < best_energy:
            best_energy = energy
            best_hw = total_hw
            best_delta_W = delta

    return best_energy, best_hw, n_attempts


def main():
    np.random.seed(42)

    print("=" * 70)
    print("НАТИВНЫЕ ИНСТРУМЕНТЫ НАШЕГО ИЗМЕРЕНИЯ")
    print("=" * 70)

    # ===================================================================
    print("\n" + "=" * 70)
    print("1. МАТРИЦА СВЯЗНОСТИ (вместо correlation)")
    print("=" * 70)

    # Вычисляем coupling для всех пар H[i], H[j]
    print("\n  Вычисляем coupling matrix (нужно birthday на каждую пару)...")

    coupling_matrix = np.zeros((8, 8))
    mean_cond_matrix = np.zeros((8, 8))

    reg = ['H0(a)', 'H1(b)', 'H2(c)', 'H3(d)', 'H4(e)', 'H5(f)', 'H6(g)', 'H7(h)']
    types = ['NODE', 'PIPE', 'PIPE', 'PIPE', 'NODE', 'PIPE', 'PIPE', 'PIPE']

    for i in range(8):
        for j in range(8):
            if i == j:
                coupling_matrix[i][j] = 1.0
                mean_cond_matrix[i][j] = 0.0
            else:
                c, mc, mu = measure_coupling(i, j, n_samples=500)
                coupling_matrix[i][j] = c
                mean_cond_matrix[i][j] = mc

    print(f"\n  МАТРИЦА СВЯЗНОСТИ (coupling):")
    print(f"  {'':>8}", end="")
    for j in range(8):
        print(f" {reg[j]:>7}", end="")
    print()

    for i in range(8):
        print(f"  {reg[i]:>8}", end="")
        for j in range(8):
            if i == j:
                print(f"    --- ", end="")
            else:
                c = coupling_matrix[i][j]
                marker = "★" if c > 0.02 else " "
                print(f"  {c:.3f}{marker}", end="")
        print(f"  [{types[i]}]")

    print(f"\n  Coupling > 0.02 (★) = позиции СВЯЗАНЫ через ткань")

    # Pipe-chain verification
    print(f"\n  Pipe-chain coupling:")
    print(f"    H[7]→H[6] (e-chain): {coupling_matrix[7][6]:.4f}")
    print(f"    H[6]→H[5] (e-chain): {coupling_matrix[6][5]:.4f}")
    print(f"    H[3]→H[2] (a-chain): {coupling_matrix[3][2]:.4f}")
    print(f"    H[2]→H[1] (a-chain): {coupling_matrix[2][1]:.4f}")
    print(f"    H[7]→H[4] (cross):   {coupling_matrix[7][4]:.4f}")
    print(f"    H[3]→H[0] (cross):   {coupling_matrix[3][0]:.4f}")

    # ===================================================================
    print(f"\n{'=' * 70}")
    print("2. ЭНЕРГИЯ МЕТКИ (вместо HW)")
    print("=" * 70)

    # Сравниваем: HW vs Energy для случайных пар
    hw_list = []
    energy_list = []

    for _ in range(10000):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W1)
        W2[0] ^= (1 << np.random.randint(0, 32))
        H1 = sha256_words(W1)
        H2 = sha256_words(W2)
        delta_H = [H1[i] ^ H2[i] for i in range(8)]

        total_hw = sum(hw(d) for d in delta_H)
        energy = mark_energy(delta_H, coupling_matrix)

        hw_list.append(total_hw)
        energy_list.append(energy)

    print(f"\n  10K случайных пар:")
    print(f"    HW(δH):     mean={np.mean(hw_list):.1f} ± {np.std(hw_list):.1f}")
    print(f"    Energy(δ):   mean={np.mean(energy_list):.1f} ± {np.std(energy_list):.1f}")
    print(f"    Ratio E/HW:  {np.mean(energy_list)/np.mean(hw_list):.3f}")
    print(f"    Corr(E,HW):  {np.corrcoef(hw_list, energy_list)[0,1]:.4f}")

    if np.mean(energy_list) < np.mean(hw_list) * 0.95:
        saving = (1 - np.mean(energy_list)/np.mean(hw_list)) * 100
        print(f"\n    >>> Энергия МЕНЬШЕ HW на {saving:.1f}% — coupling снижает стоимость!")

    # ===================================================================
    print(f"\n{'=' * 70}")
    print("3. НАВИГАЦИЯ (вместо birthday)")
    print("=" * 70)

    # Навигация vs random search
    nav_energy, nav_hw, nav_attempts = navigate_fabric(coupling_matrix, n_attempts=50000)

    # Random baseline
    np.random.seed(99)
    rand_best_hw = 256
    for _ in range(50000):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W1)
        W2[0] ^= (1 << np.random.randint(0, 32))
        H1 = sha256_words(W1)
        H2 = sha256_words(W2)
        h = sum(hw(H1[i] ^ H2[i]) for i in range(8))
        if h < rand_best_hw:
            rand_best_hw = h

    print(f"\n  50K попыток:")
    print(f"    Навигация: best HW={nav_hw}, best Energy={nav_energy:.1f}")
    print(f"    Random:    best HW={rand_best_hw}")
    print(f"    → {'НАВИГАЦИЯ ЛУЧШЕ!' if nav_hw < rand_best_hw else 'Одинаково' if nav_hw == rand_best_hw else 'Random лучше'}")

    # ===================================================================
    print(f"\n{'=' * 70}")
    print("4. НАТИВНЫЕ ЕДИНИЦЫ ИЗМЕРЕНИЯ")
    print("=" * 70)

    # Сколько "мод ткани" нужно обнулить для collision?
    # Мода = независимое направление в coupling-пространстве

    eigenvalues = np.sort(np.linalg.eigvalsh(coupling_matrix))[::-1]
    significant = sum(1 for ev in eigenvalues if ev > 0.1 * eigenvalues[0])

    print(f"\n  Eigenvalues of coupling matrix:")
    for i, ev in enumerate(eigenvalues):
        bar = "█" * int(ev * 20)
        print(f"    Mode {i}: λ={ev:.4f}  {bar}")

    print(f"\n  Significant modes: {significant} / 8")
    print(f"\n  В нашем измерении:")
    print(f"    Старая единица: 256 бит = 2^128 birthday")
    print(f"    Наша единица:   {significant} мод × ? стоимость/мода")
    print(f"    Стоимость моды зависит от coupling — НЕ фиксирована!")


if __name__ == "__main__":
    main()
