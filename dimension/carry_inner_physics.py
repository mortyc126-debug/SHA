"""
Физика ВНУТРИ carry — 1D клеточный автомат.

Carry chain: c_i = MAJ(x_i, y_i, c_{i-1})
Это Rule 232 в классификации Wolfram.

Вопросы:
  1. Есть ли фазовые переходы в carry?
  2. Есть ли "скорость света" (макс. скорость распространения)?
  3. Есть ли интерференция carry-цепочек?
  4. Есть ли корреляции в ПРОСТРАНСТВЕ carry (не во времени)?
  5. Можно ли предсказать carry из частичной информации?
"""

import numpy as np
import hashlib
import struct

MASK32 = 0xFFFFFFFF

def hw(x): return bin(x).count('1')

def get_carry_chain(x, y):
    """Возвращает полную carry-цепочку сложения x + y."""
    c = 0
    carries = []
    for i in range(32):
        xi = (x >> i) & 1
        yi = (y >> i) & 1
        s = xi + yi + c
        c = s >> 1
        carries.append(c)
    return carries  # carries[i] = carry OUT of position i

def carry_chain_stats(x, y):
    """Анализ carry-цепочки: длины серий, паттерны."""
    chain = get_carry_chain(x, y)

    # Длины серий единиц (carry runs)
    runs = []
    current_run = 0
    for c in chain:
        if c == 1:
            current_run += 1
        else:
            if current_run > 0:
                runs.append(current_run)
            current_run = 0
    if current_run > 0:
        runs.append(current_run)

    # Число переключений (0→1, 1→0)
    transitions = sum(1 for i in range(1, 32) if chain[i] != chain[i-1])

    return {
        'chain': chain,
        'hw': sum(chain),
        'runs': runs,
        'max_run': max(runs) if runs else 0,
        'num_runs': len(runs),
        'transitions': transitions,
    }


def experiment_carry_inner(num_samples=10000):
    np.random.seed(42)

    print("=" * 70)
    print("ФИЗИКА ВНУТРИ CARRY — 1D клеточный автомат")
    print("=" * 70)

    # === 1. Базовая статистика carry-цепочки ===
    print("\n" + "=" * 70)
    print("1. БАЗОВАЯ СТАТИСТИКА CARRY-ЦЕПОЧКИ")
    print("=" * 70)

    all_hw = []
    all_max_run = []
    all_num_runs = []
    all_transitions = []
    all_chains = []

    for _ in range(num_samples):
        x = np.random.randint(0, 2**31) * 2 + np.random.randint(0, 2)
        y = np.random.randint(0, 2**31) * 2 + np.random.randint(0, 2)
        stats = carry_chain_stats(x, y)
        all_hw.append(stats['hw'])
        all_max_run.append(stats['max_run'])
        all_num_runs.append(stats['num_runs'])
        all_transitions.append(stats['transitions'])
        all_chains.append(stats['chain'])

    print(f"\n  Случайные x, y ∈ Z/2³²:")
    print(f"    HW(carry chain):      {np.mean(all_hw):.2f} ± {np.std(all_hw):.2f}  (из 32)")
    print(f"    Теор. P(carry=1):     0.5 → ожидаем HW = 16")
    print(f"    Макс. серия carry=1:  {np.mean(all_max_run):.2f} ± {np.std(all_max_run):.2f}")
    print(f"    Число серий:          {np.mean(all_num_runs):.2f} ± {np.std(all_num_runs):.2f}")
    print(f"    Переключения:         {np.mean(all_transitions):.2f} ± {np.std(all_transitions):.2f}")

    # === 2. Пространственная корреляция carry ===
    print("\n" + "=" * 70)
    print("2. ПРОСТРАНСТВЕННАЯ КОРРЕЛЯЦИЯ CARRY")
    print("=" * 70)

    # Корреляция carry[i] и carry[i+d] для разных d
    chains_arr = np.array(all_chains)  # (num_samples, 32)

    print(f"\n  Корреляция carry[i] vs carry[i+d]:")
    for d in [1, 2, 3, 4, 5, 8, 16]:
        corrs = []
        for i in range(32 - d):
            c1 = chains_arr[:, i]
            c2 = chains_arr[:, i + d]
            if np.std(c1) > 0 and np.std(c2) > 0:
                corrs.append(np.corrcoef(c1, c2)[0, 1])
        mean_corr = np.mean(corrs) if corrs else 0
        bar = "█" * int(abs(mean_corr) * 50)
        print(f"    d={d:2d}: corr = {mean_corr:+.4f}  {bar}")

    # === 3. Carry при МАЛОМ δ (дифференциальный) ===
    print("\n" + "=" * 70)
    print("3. ДИФФЕРЕНЦИАЛЬНЫЙ CARRY: δ = 1 бит")
    print("=" * 70)

    # x + y vs (x ⊕ δ) + y — как меняется carry chain при 1-бит δ?
    diff_positions = np.zeros(32)  # где carry меняется
    diff_counts = []

    for _ in range(num_samples):
        x = np.random.randint(0, 2**31) * 2 + np.random.randint(0, 2)
        y = np.random.randint(0, 2**31) * 2 + np.random.randint(0, 2)
        bit = np.random.randint(0, 32)
        x2 = x ^ (1 << bit)

        c1 = get_carry_chain(x, y)
        c2 = get_carry_chain(x2, y)

        diff = [c1[i] ^ c2[i] for i in range(32)]
        for i in range(32):
            diff_positions[i] += diff[i]
        diff_counts.append(sum(diff))

    diff_positions /= num_samples

    print(f"\n  P(carry[i] меняется) при δ = 1 случайный бит:")
    for i in range(32):
        bar = "█" * int(diff_positions[i] * 50)
        print(f"    бит {i:2d}: P = {diff_positions[i]:.4f}  {bar}")

    print(f"\n  Среднее число изменённых carry-бит: {np.mean(diff_counts):.2f} ± {np.std(diff_counts):.2f}")

    # === 3b. δ в конкретном бите ===
    print(f"\n  Зависимость от ПОЗИЦИИ δ-бита:")
    for delta_bit in [0, 4, 8, 15, 16, 24, 28, 31]:
        diffs_for_bit = []
        propagation_lengths = []
        for _ in range(num_samples):
            x = np.random.randint(0, 2**31) * 2 + np.random.randint(0, 2)
            y = np.random.randint(0, 2**31) * 2 + np.random.randint(0, 2)
            x2 = x ^ (1 << delta_bit)

            c1 = get_carry_chain(x, y)
            c2 = get_carry_chain(x2, y)

            diff = [c1[i] ^ c2[i] for i in range(32)]
            total_diff = sum(diff)
            diffs_for_bit.append(total_diff)

            # Длина распространения: от delta_bit до последнего изменённого carry
            if total_diff > 0:
                last_changed = max(i for i in range(32) if diff[i])
                propagation_lengths.append(last_changed - delta_bit)

        mean_prop = np.mean(propagation_lengths) if propagation_lengths else 0
        print(f"    δ бит {delta_bit:2d}: mean Δcarry = {np.mean(diffs_for_bit):5.2f}, propagation = {mean_prop:.1f} бит")

    # === 4. Фазовые переходы ===
    print("\n" + "=" * 70)
    print("4. ФАЗОВЫЕ ПЕРЕХОДЫ: CARRY VS ПЛОТНОСТЬ ЕДИНИЦ")
    print("=" * 70)

    # При какой плотности единиц в x,y carry ведёт себя по-разному?
    for density in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
        hws = []
        max_runs = []
        for _ in range(num_samples):
            # Генерим x с заданной плотностью единиц
            x = 0
            for b in range(32):
                if np.random.random() < density:
                    x |= (1 << b)
            y = 0
            for b in range(32):
                if np.random.random() < density:
                    y |= (1 << b)

            stats = carry_chain_stats(x, y)
            hws.append(stats['hw'])
            max_runs.append(stats['max_run'])

        expected_carry = 1 - (1 - density)**2 / (1 + density**2 - (1-density)**2 + 0.001)
        print(f"    density={density:.1f}: HW(carry)={np.mean(hws):5.1f}/32, max_run={np.mean(max_runs):4.1f}")

    # === 5. Каскадный carry: ADD(ADD(x,y),z) ===
    print("\n" + "=" * 70)
    print("5. КАСКАДНЫЙ CARRY: ИНТЕРФЕРЕНЦИЯ")
    print("=" * 70)

    # В SHA-256 NODE_a = T1 + T2 = (((h + Σ₁) + Ch) + K + W) + (Σ₀ + Maj)
    # Carry каскадируется. Интерферируют ли carry-цепочки?

    cascade_hw = []
    independent_hw = []

    for _ in range(num_samples):
        a = np.random.randint(0, 2**31) * 2 + np.random.randint(0, 2)
        b = np.random.randint(0, 2**31) * 2 + np.random.randint(0, 2)
        c = np.random.randint(0, 2**31) * 2 + np.random.randint(0, 2)

        # Каскад: (a + b) + c
        ab = (a + b) & MASK32
        carry_ab = get_carry_chain(a, b)
        carry_abc = get_carry_chain(ab, c)
        cascade_hw.append(sum(carry_abc))

        # Независимый: a + c (без b)
        carry_ac = get_carry_chain(a, c)
        independent_hw.append(sum(carry_ac))

    print(f"\n  Каскад (a+b)+c — carry второго ADD:")
    print(f"    HW(carry): {np.mean(cascade_hw):.2f} ± {np.std(cascade_hw):.2f}")
    print(f"  Независимый a+c:")
    print(f"    HW(carry): {np.mean(independent_hw):.2f} ± {np.std(independent_hw):.2f}")
    print(f"  Разница: {np.mean(cascade_hw) - np.mean(independent_hw):.4f}")
    print(f"  → {'ИНТЕРФЕРЕНЦИЯ!' if abs(np.mean(cascade_hw) - np.mean(independent_hw)) > 0.2 else 'Нет интерференции'}")

    # Carry overlap между двумя ADD в каскаде
    overlaps = []
    for _ in range(num_samples):
        a = np.random.randint(0, 2**31) * 2 + np.random.randint(0, 2)
        b = np.random.randint(0, 2**31) * 2 + np.random.randint(0, 2)
        c = np.random.randint(0, 2**31) * 2 + np.random.randint(0, 2)

        ab = (a + b) & MASK32
        c1 = get_carry_chain(a, b)
        c2 = get_carry_chain(ab, c)

        # Jaccard overlap
        both = sum(c1[i] & c2[i] for i in range(32))
        either = sum(c1[i] | c2[i] for i in range(32))
        overlaps.append(both / max(either, 1))

    print(f"\n  Carry overlap (Jaccard) между ADD₁ и ADD₂ в каскаде:")
    print(f"    Jaccard = {np.mean(overlaps):.4f}")
    print(f"    Random (независимые): ~0.333")
    print(f"    → {'СВЯЗАНЫ!' if abs(np.mean(overlaps) - 0.333) > 0.02 else 'Независимы'}")

    # === 6. Предсказуемость carry ===
    print("\n" + "=" * 70)
    print("6. ПРЕДСКАЗУЕМОСТЬ: МОЖНО ЛИ УГАДАТЬ CARRY ИЗ ЧАСТИЧНОЙ ИНФОРМАЦИИ?")
    print("=" * 70)

    # Знаем x[0..15] и y[0..15] — можем ли предсказать carry[16..31]?
    correct_predictions = 0
    total_predictions = 0

    for _ in range(num_samples):
        x = np.random.randint(0, 2**31) * 2 + np.random.randint(0, 2)
        y = np.random.randint(0, 2**31) * 2 + np.random.randint(0, 2)

        full_chain = get_carry_chain(x, y)

        # Предсказание: carry[15] (выход нижней половины) определяет
        # поведение carry в верхней половине?
        carry_15 = full_chain[15]

        # Наивное предсказание: если carry[15]=1, то carry[16..31] смещены к 1
        for i in range(16, 32):
            predicted = carry_15  # наивная экстраполяция
            actual = full_chain[i]
            if predicted == actual:
                correct_predictions += 1
            total_predictions += 1

    accuracy = correct_predictions / total_predictions
    print(f"\n  Предсказание carry[16..31] по carry[15]:")
    print(f"    Accuracy = {accuracy:.4f}")
    print(f"    Random baseline = 0.5000")
    print(f"    Gain = {(accuracy - 0.5) * 100:.2f}%")

    # Более умное предсказание: по x[0..15], y[0..15], carry[15]
    correct_smart = 0
    total_smart = 0

    for _ in range(num_samples):
        x = np.random.randint(0, 2**31) * 2 + np.random.randint(0, 2)
        y = np.random.randint(0, 2**31) * 2 + np.random.randint(0, 2)

        full_chain = get_carry_chain(x, y)

        # "Умное" предсказание: знаем carry[15] и бит x[16], y[16]
        # carry[16] = MAJ(x[16], y[16], carry[15]) — ТОЧНО!
        carry_15 = full_chain[15]
        x16 = (x >> 16) & 1
        y16 = (y >> 16) & 1
        predicted_16 = (x16 & y16) | (x16 & carry_15) | (y16 & carry_15)
        actual_16 = full_chain[16]

        if predicted_16 == actual_16:
            correct_smart += 1
        total_smart += 1

    accuracy_smart = correct_smart / total_smart
    print(f"\n  Предсказание carry[16] по x[16], y[16], carry[15] (формула MAJ):")
    print(f"    Accuracy = {accuracy_smart:.4f}")
    print(f"    → {'ТОЧНОЕ (carry вычислим!)' if accuracy_smart > 0.99 else 'Неточное'}")

    # А если мы НЕ знаем x[16], y[16]?
    correct_partial = 0
    total_partial = 0

    for _ in range(num_samples):
        x = np.random.randint(0, 2**31) * 2 + np.random.randint(0, 2)
        y = np.random.randint(0, 2**31) * 2 + np.random.randint(0, 2)

        full_chain = get_carry_chain(x, y)
        carry_15 = full_chain[15]

        # Без знания x[16], y[16]: лучший предиктор = carry[15]
        predicted = carry_15
        actual = full_chain[16]
        if predicted == actual:
            correct_partial += 1
        total_partial += 1

    accuracy_partial = correct_partial / total_partial
    print(f"\n  Предсказание carry[16] ТОЛЬКО по carry[15] (без x[16], y[16]):")
    print(f"    Accuracy = {accuracy_partial:.4f}")
    print(f"    Gain over random = {(accuracy_partial - 0.5) * 100:.2f}%")

    # === 7. ВЕРДИКТ ===
    print("\n" + "=" * 70)
    print("7. ВЕРДИКТ: ГДЕ СТРОИТЬ ФИЗИКУ?")
    print("=" * 70)

    print(f"""
  CARRY ИЗНУТРИ — что мы нашли:

  1. ПРОСТРАНСТВЕННАЯ КОРРЕЛЯЦИЯ: carry[i] и carry[i+1] сильно коррелированы
     (серии единиц / нулей). Carry — НЕ случайная последовательность.

  2. РАСПРОСТРАНЕНИЕ: δ в бите k влияет на carry[k..31], но НЕ на carry[0..k-1].
     Carry распространяется ТОЛЬКО ВВЕРХ (от LSB к MSB).
     "Скорость света" = 1 бит/позицию. Строгая причинность.

  3. КАСКАДНЫЙ CARRY: два ADD в каскаде — carries связаны или нет?

  4. ПРЕДСКАЗУЕМОСТЬ:
     - Зная x[i], y[i], carry[i-1] → carry[i] ТОЧНО (формула MAJ).
     - Зная ТОЛЬКО carry[i-1] → carry[i] с accuracy {accuracy_partial:.1%}.
     - Это +{(accuracy_partial - 0.5) * 100:.1f}% над random.

  КЛЮЧЕВОЙ ВЫВОД:
     Carry ВЫЧИСЛИМ, если знать (x, y, c_prev).
     В SHA-256 мы НЕ знаем промежуточные x, y — они зависят от ВСЕХ W.

     Барьер НЕ в carry как таковом.
     Барьер в том, что ВХОДЫ каждого ADD зависят от предыдущих ВЫХОДОВ.

     Это РЕКУРСИВНАЯ ЗАВИСИМОСТЬ:
       carry(round r) = f(state(r)) = f(g(state(r-1))) = f(g(h(state(r-2)))) ...

     Carry — симптом. Болезнь — рекурсивная глубина 64.
""")


if __name__ == "__main__":
    experiment_carry_inner(num_samples=10000)
