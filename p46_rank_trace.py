"""Диагностика: меняется ли ранг J по ходу подъёма?"""
import random, sys
sys.path.insert(0, '/home/user/SHA')
from p46_dw0_hypothesis import (
    numerical_jacobian_15x15, rank_mod2, gauss_mod2,
    compute_all_constraints, MASK
)

random.seed(42)

for msg_idx in range(3):
    W_base = [random.randint(0, MASK) for _ in range(16)]
    print(f"\n=== Сообщение {msg_idx+1} ===")

    # Найдём DW0 с rank=15
    found_dw0 = None
    for dw0 in range(1, 1000):
        dW = {0: dw0}
        for j in range(1, 16):
            dW[j] = 0
        J = numerical_jacobian_15x15(W_base, dW)
        r = rank_mod2(J)
        if r == 15:
            found_dw0 = dw0
            break
    if found_dw0 is None:
        print("  DW0 с rank=15 не найден")
        continue

    print(f"  DW0={found_dw0} (rank=15 при DW1..DW15=0)")

    # Трассируем побитовый подъём, печатаем rank на каждом бите
    dW = {0: found_dw0}
    for j in range(1, 16):
        dW[j] = 0

    for bit in range(10):
        f_vec = compute_all_constraints(W_base, dW)
        norm = sum(1 for v in f_vec if v != 0)
        J = numerical_jacobian_15x15(W_base, dW)
        r = rank_mod2(J)
        rhs = [(-v) % 2 for v in f_vec]
        delta = gauss_mod2(J, rhs)

        # Показываем состояние
        print(f"  bit={bit}: norm={norm}  rank(J)={r}  delta={'OK' if delta else 'NONE'}")

        if delta is None:
            print(f"    --> SINGULAR at bit={bit}, остановка")
            break

        if norm == 0:
            print(f"    --> РЕШЕНИЕ НАЙДЕНО на bit={bit}!")
            break

        for j in range(15):
            var_idx = j + 1
            dW[var_idx] = (dW.get(var_idx, 0) + delta[j] * (1 << bit)) & MASK

    f_final = compute_all_constraints(W_base, dW)
    norm_final = sum(1 for v in f_final if v != 0)
    print(f"  Итог: norm={norm_final}")
