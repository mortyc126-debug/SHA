"""
П-48: быстрый исчерпывающий поиск mod 4.
Иерархическое просеивание с ранним выходом.
"""
import random
import sys
sys.path.insert(0, '/home/user/SHA')
from p48_mod4_barrier import sha256_state, MASK, find_seed_mod2

random.seed(42)


def Da_mod4(W_base, dW_vals16, r):
    W2 = list(W_base)
    for i in range(16):
        W2[i] = (W_base[i] + dW_vals16[i]) & MASK
    s1 = sha256_state(W_base, r)
    s2 = sha256_state(W2, r)
    return (s2[0] - s1[0]) % 4


def De17_mod4(W_base, dW_vals16):
    W2 = list(W_base)
    for i in range(16):
        W2[i] = (W_base[i] + dW_vals16[i]) & MASK
    s1 = sha256_state(W_base, 17)
    s2 = sha256_state(W2, 17)
    return (s2[4] - s1[4]) % 4


def exhaustive_mod4(W_base, seed_dict):
    """
    Полный перебор δ ∈ {0,1}^15 с ранним выходом по первым 3 ограничениям.
    seed_dict: {idx: delta_value} - базовое семя (f ≡ 0 mod 2).
    Ищем δ такое, что f(seed + 2*δ) ≡ 0 mod 4.
    """
    seed_vals = [seed_dict.get(i, 0) for i in range(16)]
    found = []

    for bits in range(1 << 15):
        dW = list(seed_vals)
        for j in range(15):
            dW[j + 1] = (seed_vals[j + 1] + ((bits >> j) & 1) * 2) & MASK

        # Ранний выход по первым трём
        if Da_mod4(W_base, dW, 3) != 0:
            continue
        if Da_mod4(W_base, dW, 4) != 0:
            continue
        if Da_mod4(W_base, dW, 5) != 0:
            continue

        # Остальные ограничения
        ok = True
        for r in range(6, 17):
            if Da_mod4(W_base, dW, r) != 0:
                ok = False
                break
        if not ok:
            continue
        if De17_mod4(W_base, dW) != 0:
            continue

        found.append(bits)
        if len(found) >= 3:
            break

    return found


print("П-48: исчерпывающий поиск mod 4")
print("=" * 60)
print("Ищем δ ∈ {0,1}^15 с f(seed + 2δ) ≡ 0 mod 4")
print()

mod4_found_total = 0
tested_seeds = 0

for msg_idx in range(8):
    W_base = [random.randint(0, MASK) for _ in range(16)]
    found_any = False

    for _ in range(6):
        dw0 = random.randint(1, MASK)
        seed = find_seed_mod2(W_base, dw0, n_trials=3000)
        if seed is None:
            continue
        tested_seeds += 1

        found = exhaustive_mod4(W_base, seed)

        if found:
            mod4_found_total += 1
            found_any = True
            print(f"  msg={msg_idx+1} DW0=0x{dw0:08x}: MOD4 НАЙДЕНО! {len(found)} решений")
        else:
            print(f"  msg={msg_idx+1} DW0=0x{dw0:08x}: нет решений mod 4")

    print()

print(f"ИТОГ: mod4_найдено={mod4_found_total}/{tested_seeds}")
if mod4_found_total == 0:
    print("T_MOD4_BARRIER_ABSOLUTE: система не разрешима mod 4 ни для одного семени.")
    print("T_GLOBAL_BARRIER_FINAL подтверждён исчерпывающим перебором.")
else:
    print("T_MOD4_SOLVABLE: решения mod 4 существуют!")
    print(f"Плотность: ~{mod4_found_total}/{tested_seeds} семян имеют продолжение.")
