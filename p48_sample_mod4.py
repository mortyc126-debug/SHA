"""
П-48 (быстрый): случайная выборка 10000 δ ∈ {0,1}^15 для поиска mod-4 решений.
Вместо полного перебора 2^15 — вероятностный тест.
Если плотность решений > 0 → найдём за O(1/density) попыток.
Если решений нет → 10000 попыток дают P(miss) < 2^{-14} при density≥1.
"""
import random, sys, math
sys.path.insert(0, '/home/user/SHA')
from p48_mod4_barrier import sha256_state, MASK, find_seed_mod2

random.seed(42)


def f_mod4(W_base, seed_vals16, delta_bits):
    """Проверяет f(seed + 2*δ) ≡ 0 mod 4 с ранним выходом. delta_bits = int 0..2^15."""
    dW = list(seed_vals16)
    for j in range(15):
        dW[j + 1] = (seed_vals16[j + 1] + ((delta_bits >> j) & 1) * 2) & MASK

    # Ранний выход по первым 3 ограничениям
    for r in range(3, 6):
        s1 = sha256_state(W_base, r)
        W2 = list(W_base)
        for i in range(16):
            W2[i] = (W_base[i] + dW[i]) & MASK
        s2 = sha256_state(W2, r)
        if (s2[0] - s1[0]) % 4 != 0:
            return False

    W2 = list(W_base)
    for i in range(16):
        W2[i] = (W_base[i] + dW[i]) & MASK

    # Остальные Da
    for r in range(6, 17):
        s1 = sha256_state(W_base, r)
        s2 = sha256_state(W2, r)
        if (s2[0] - s1[0]) % 4 != 0:
            return False

    # De17
    s1 = sha256_state(W_base, 17)
    s2 = sha256_state(W2, 17)
    if (s2[4] - s1[4]) % 4 != 0:
        return False

    return True


def sample_mod4(W_base, seed_dict, n_samples=5000):
    """Случайная выборка n_samples δ ∈ {0,1}^15."""
    seed_vals = [seed_dict.get(i, 0) for i in range(16)]
    found = []
    for _ in range(n_samples):
        bits = random.randrange(1 << 15)
        if f_mod4(W_base, seed_vals, bits):
            found.append(bits)
            if len(found) >= 5:
                break
    return found


print("П-48: случайная выборка mod 4 (N=5000 на семя)")
print("=" * 60)
print()

mod4_found = 0
no_mod4 = 0
total_seeds = 0

for msg_idx in range(10):
    W_base = [random.randint(0, MASK) for _ in range(16)]
    msg_has = False

    for _ in range(6):
        dw0 = random.randint(1, MASK)
        seed = find_seed_mod2(W_base, dw0, n_trials=3000)
        if seed is None:
            continue
        total_seeds += 1

        found = sample_mod4(W_base, seed, n_samples=5000)
        if found:
            mod4_found += 1
            msg_has = True
            print(f"  msg={msg_idx+1} DW0=0x{dw0:08x}: MOD4 НАЙДЕНО! {len(found)} решений")
        else:
            no_mod4 += 1
            print(f"  msg={msg_idx+1} DW0=0x{dw0:08x}: 5000 попыток — нет mod 4")

    print()

print(f"ИТОГ: найдено={mod4_found}/{total_seeds}  нет_решений={no_mod4}")
density_upper = 3 / max(no_mod4 * 5000, 1)
log2_density = math.log2(density_upper) if density_upper > 0 else -float('inf')
if mod4_found == 0:
    print(f"Верхняя оценка плотности решений mod 4: < 3/{no_mod4*5000} ≈ 2^{log2_density:.1f}")
    print()
    print("ВЫВОД:")
    print("  Если решений mod 4 нет вообще → T_MOD4_BARRIER_ABSOLUTE.")
    print("  Если есть, но редкие → нужен более широкий поиск.")
else:
    frac = mod4_found / total_seeds
    print(f"  Плотность: {mod4_found}/{total_seeds} семян продолжаются mod 4.")
    print(f"  T_MOD4_SOLVABLE! Дальше → k-шаговый подъём до mod 2^32.")
