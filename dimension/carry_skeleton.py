"""
Carry Skeleton — Этап B: Построение carry-скелета SHA-256
в терминах нового измерения (Ткань-Переход-Слой)

Разделяем SHA-256 на:
  L(W) — линейная часть (XOR, ротации, σ₀, σ₁, линейные Ch/Maj)
  C(W) — carry-часть (все carry от ADD mod 2³²)

Строим карту: какие carries коррелированы, где плотность минимальна.
"""

import struct
import numpy as np
from typing import List, Tuple, Dict

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

IV = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
]

MASK32 = 0xFFFFFFFF


def rotr(x, n):
    return ((x >> n) | (x << (32 - n))) & MASK32

def shr(x, n):
    return x >> n

def sigma0(x):
    return rotr(x, 7) ^ rotr(x, 18) ^ shr(x, 3)

def sigma1(x):
    return rotr(x, 17) ^ rotr(x, 19) ^ shr(x, 10)

def Sigma0(x):
    return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)

def Sigma1(x):
    return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)

def Ch(e, f, g):
    return (e & f) ^ (~e & g) & MASK32

def Maj(a, b, c):
    return (a & b) ^ (a & c) ^ (b & c)


def add32_with_carry_info(x, y):
    """Сложение mod 2^32 с подсчётом carry-бит."""
    result = (x + y) & MASK32
    # carry на каждой битовой позиции
    carry_bits = 0
    c = 0
    for i in range(32):
        xi = (x >> i) & 1
        yi = (y >> i) & 1
        s = xi + yi + c
        c = s >> 1
        if c:
            carry_bits |= (1 << i)
    return result, carry_bits, bin(carry_bits).count('1')


def sha256_with_carry_tracking(W_input: List[int]):
    """
    SHA-256 с полным трекингом carry на каждом сложении.
    Возвращает:
      - hash
      - carry_map: словарь {(раунд, операция): (carry_bits, carry_count)}
      - round_states: состояния регистров на каждом раунде
    """
    # Schedule expansion с трекингом carry
    W = list(W_input[:16])
    schedule_carries = {}

    for r in range(16, 64):
        s0 = sigma0(W[r-15])
        s1 = sigma1(W[r-2])

        # W[r] = s1 + W[r-7] + s0 + W[r-16]  (3 сложения)
        tmp1, c1, n1 = add32_with_carry_info(s1, W[r-7])
        tmp2, c2, n2 = add32_with_carry_info(tmp1, s0)
        tmp3, c3, n3 = add32_with_carry_info(tmp2, W[r-16])

        schedule_carries[(r, 'add1')] = (c1, n1)
        schedule_carries[(r, 'add2')] = (c2, n2)
        schedule_carries[(r, 'add3')] = (c3, n3)
        W.append(tmp3)

    # Round function с трекингом carry
    a, b, c, d, e, f, g, h = IV
    round_carries = {}
    round_states = [(a, b, c, d, e, f, g, h)]

    for r in range(64):
        S1 = Sigma1(e)
        ch = Ch(e, f, g)
        S0 = Sigma0(a)
        maj = Maj(a, b, c)

        # T1 = h + S1 + ch + K[r] + W[r]  (4 сложения)
        t1_1, c1, n1 = add32_with_carry_info(h, S1)
        t1_2, c2, n2 = add32_with_carry_info(t1_1, ch)
        t1_3, c3, n3 = add32_with_carry_info(t1_2, K[r])
        T1,   c4, n4 = add32_with_carry_info(t1_3, W[r])

        round_carries[(r, 'T1_add1')] = (c1, n1)
        round_carries[(r, 'T1_add2')] = (c2, n2)
        round_carries[(r, 'T1_add3')] = (c3, n3)
        round_carries[(r, 'T1_add4')] = (c4, n4)

        # T2 = S0 + maj  (1 сложение)
        T2, c5, n5 = add32_with_carry_info(S0, maj)
        round_carries[(r, 'T2_add1')] = (c5, n5)

        # e' = d + T1  (1 сложение)
        e_new, c6, n6 = add32_with_carry_info(d, T1)
        round_carries[(r, 'e_add')] = (c6, n6)

        # a' = T1 + T2  (1 сложение)
        a_new, c7, n7 = add32_with_carry_info(T1, T2)
        round_carries[(r, 'a_add')] = (c7, n7)

        # Shift registers (ТРУБЫ — нулевая нелинейность)
        h = g
        g = f
        f = e
        e = e_new
        d = c
        c = b
        b = a
        a = a_new

        round_states.append((a, b, c, d, e, f, g, h))

    # Финализация
    final_carries = {}
    H = list(IV)
    for i, val in enumerate([a, b, c, d, e, f, g, h]):
        result, cf, nf = add32_with_carry_info(H[i], val)
        final_carries[('final', i)] = (cf, nf)
        H[i] = result

    carry_map = {**schedule_carries, **round_carries, **final_carries}
    return H, carry_map, round_states


def analyze_carry_skeleton(num_samples=100):
    """
    Анализ carry-скелета SHA-256.
    """
    np.random.seed(42)

    print("=" * 70)
    print("CARRY SKELETON — Анализ нелинейности SHA-256")
    print("Новое измерение: Ткань-Переход-Слой")
    print("=" * 70)

    # Статистика carry по операциям
    carry_stats = {}
    carry_correlation = {}

    all_carry_maps = []

    for trial in range(num_samples):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H, carry_map, states = sha256_with_carry_tracking(W)
        all_carry_maps.append(carry_map)

        for key, (cbits, count) in carry_map.items():
            if key not in carry_stats:
                carry_stats[key] = []
            carry_stats[key].append(count)

    # === 1. Среднее carry по типам операций ===
    print("\n" + "=" * 70)
    print("1. СРЕДНЯЯ CARRY-НЕЛИНЕЙНОСТЬ ПО ТИПАМ ОПЕРАЦИЙ")
    print("=" * 70)

    type_groups = {}
    for key, counts in carry_stats.items():
        if isinstance(key[0], int) and key[0] < 64:
            op_type = key[1]  # T1_add1, T1_add2, etc.
        elif isinstance(key[0], int) and key[0] >= 16:
            op_type = f"schedule_{key[1]}"
        else:
            op_type = "final"

        if op_type not in type_groups:
            type_groups[op_type] = []
        type_groups[op_type].extend(counts)

    for op_type in sorted(type_groups.keys()):
        counts = type_groups[op_type]
        mean_c = np.mean(counts)
        std_c = np.std(counts)
        print(f"  {op_type:20s}: mean={mean_c:.2f} bits, std={std_c:.2f}")

    # === 2. Carry по раундам ===
    print("\n" + "=" * 70)
    print("2. СУММАРНАЯ CARRY-НЕЛИНЕЙНОСТЬ ПО РАУНДАМ")
    print("=" * 70)

    round_totals = np.zeros(64)
    for r in range(64):
        round_ops = [k for k in carry_stats if isinstance(k[0], int) and k[0] == r and k[1].startswith('T')]
        round_ops += [k for k in carry_stats if isinstance(k[0], int) and k[0] == r and k[1] in ('e_add', 'a_add')]
        for key in round_ops:
            round_totals[r] += np.mean(carry_stats[key])

    # Зоны
    ctrl_zone = round_totals[:14]
    transition_zone = round_totals[14:32]
    wild_zone = round_totals[32:]

    print(f"\n  Контролируемая зона (r=0..13):  mean={np.mean(ctrl_zone):.1f} bits/round, total={np.sum(ctrl_zone):.0f} bits")
    print(f"  Переходная зона    (r=14..31): mean={np.mean(transition_zone):.1f} bits/round, total={np.sum(transition_zone):.0f} bits")
    print(f"  Дикая зона         (r=32..63): mean={np.mean(wild_zone):.1f} bits/round, total={np.sum(wild_zone):.0f} bits")
    print(f"\n  ВСЕГО раунды: {np.sum(round_totals):.0f} bits нелинейности")

    # Schedule carry
    sched_total = 0
    for r in range(16, 64):
        for add_i in ['add1', 'add2', 'add3']:
            key = (r, add_i)
            if key in carry_stats:
                sched_total += np.mean(carry_stats[key])
    print(f"  ВСЕГО schedule: {sched_total:.0f} bits нелинейности")
    print(f"  ВСЕГО (раунды + schedule + финал): {np.sum(round_totals) + sched_total + 8*np.mean([np.mean(carry_stats[('final',i)]) for i in range(8)]):.0f} bits")

    # === 3. Корреляция carry между раундами ===
    print("\n" + "=" * 70)
    print("3. КОРРЕЛЯЦИЯ CARRY МЕЖДУ РАУНДАМИ (инертный поток)")
    print("=" * 70)

    # Корреляция между NODE_e(r) и NODE_e(r+4) через трубу d
    correlations_e = []
    for r in range(60):
        key_r = (r, 'e_add')
        key_r4 = (r+4, 'e_add')
        if key_r in carry_stats and key_r4 in carry_stats:
            c_r = [all_carry_maps[t][key_r][0] for t in range(num_samples)]
            c_r4 = [all_carry_maps[t][key_r4][0] for t in range(num_samples)]
            # Побитовая корреляция
            overlap = [bin(c_r[t] & c_r4[t]).count('1') / max(1, bin(c_r[t] | c_r4[t]).count('1')) for t in range(num_samples)]
            correlations_e.append(np.mean(overlap))

    # Корреляция между NODE_a(r) и NODE_a(r+4)
    correlations_a = []
    for r in range(60):
        key_r = (r, 'a_add')
        key_r4 = (r+4, 'a_add')
        if key_r in carry_stats and key_r4 in carry_stats:
            c_r = [all_carry_maps[t][key_r][0] for t in range(num_samples)]
            c_r4 = [all_carry_maps[t][key_r4][0] for t in range(num_samples)]
            overlap = [bin(c_r[t] & c_r4[t]).count('1') / max(1, bin(c_r[t] | c_r4[t]).count('1')) for t in range(num_samples)]
            correlations_a.append(np.mean(overlap))

    print(f"\n  NODE_e(r) vs NODE_e(r+4) — корреляция carry через трубу d→e:")
    print(f"    mean={np.mean(correlations_e):.4f}, std={np.std(correlations_e):.4f}")
    print(f"    Если бы независимы: ~0.33 (Jaccard двух случайных 16-bit)")

    print(f"\n  NODE_a(r) vs NODE_a(r+4) — корреляция carry:")
    print(f"    mean={np.mean(correlations_a):.4f}, std={np.std(correlations_a):.4f}")

    # Корреляция между смежными раундами (r, r+1)
    correlations_adj = []
    for r in range(63):
        key_r = (r, 'e_add')
        key_r1 = (r+1, 'e_add')
        if key_r in carry_stats and key_r1 in carry_stats:
            c_r = [all_carry_maps[t][key_r][0] for t in range(num_samples)]
            c_r1 = [all_carry_maps[t][key_r1][0] for t in range(num_samples)]
            overlap = [bin(c_r[t] & c_r1[t]).count('1') / max(1, bin(c_r[t] | c_r1[t]).count('1')) for t in range(num_samples)]
            correlations_adj.append(np.mean(overlap))

    print(f"\n  NODE_e(r) vs NODE_e(r+1) — смежная корреляция:")
    print(f"    mean={np.mean(correlations_adj):.4f}")

    # === 4. Дифференциальный carry ===
    print("\n" + "=" * 70)
    print("4. ДИФФЕРЕНЦИАЛЬНЫЙ CARRY (метка δW)")
    print("=" * 70)

    # Сравниваем carry для W и W+δ при разных HW(δ)
    for hw_target in [1, 2, 4, 8, 16]:
        diff_carries = []
        for trial in range(num_samples):
            W1 = [np.random.randint(0, 2**32) for _ in range(16)]
            # δW с заданным HW в W[0]
            delta = 0
            bits = np.random.choice(32, hw_target, replace=False)
            for bit in bits:
                delta ^= (1 << bit)
            W2 = list(W1)
            W2[0] ^= delta

            _, cm1, _ = sha256_with_carry_tracking(W1)
            _, cm2, _ = sha256_with_carry_tracking(W2)

            # Считаем сколько carry-бит изменилось
            total_diff = 0
            for key in cm1:
                if key in cm2:
                    diff = cm1[key][0] ^ cm2[key][0]
                    total_diff += bin(diff).count('1')
            diff_carries.append(total_diff)

        print(f"  HW(δW₀)={hw_target:2d}: mean diff carry = {np.mean(diff_carries):.1f} bits, std={np.std(diff_carries):.1f}")

    # === 5. Маршрутизация — тихие раунды ===
    print("\n" + "=" * 70)
    print("5. МАРШРУТИЗАЦИЯ: ПОИСК ТИХИХ ПУТЕЙ")
    print("=" * 70)

    # Тестируем: при δW[0] = 0x80000000 (1 бит), δW[1..15] = 0
    # Какие раунды "тихие" (малый carry diff)?
    quiet_rounds = np.zeros(64)
    total_nonlin_per_round = np.zeros(64)

    for trial in range(num_samples):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W1)
        W2[0] ^= 0x80000000  # 1 бит в MSB

        _, cm1, states1 = sha256_with_carry_tracking(W1)
        _, cm2, states2 = sha256_with_carry_tracking(W2)

        for r in range(64):
            round_diff = 0
            for op in ['T1_add1', 'T1_add2', 'T1_add3', 'T1_add4', 'T2_add1', 'e_add', 'a_add']:
                key = (r, op)
                if key in cm1 and key in cm2:
                    diff = cm1[key][0] ^ cm2[key][0]
                    round_diff += bin(diff).count('1')
            total_nonlin_per_round[r] += round_diff
            if round_diff <= 2:
                quiet_rounds[r] += 1

    total_nonlin_per_round /= num_samples
    quiet_rounds /= num_samples

    print(f"\n  δW₀ = 0x80000000 (1 бит MSB), остальные δ = 0")
    print(f"\n  Зона контроля (r=0..13):")
    for r in range(14):
        bar = "█" * int(total_nonlin_per_round[r])
        quiet_pct = quiet_rounds[r] * 100
        print(f"    r={r:2d}: nonlin={total_nonlin_per_round[r]:5.1f} bits  quiet={quiet_pct:5.1f}%  {bar}")

    print(f"\n  Переходная зона (r=14..31):")
    for r in range(14, 32):
        bar = "█" * int(total_nonlin_per_round[r])
        quiet_pct = quiet_rounds[r] * 100
        print(f"    r={r:2d}: nonlin={total_nonlin_per_round[r]:5.1f} bits  quiet={quiet_pct:5.1f}%  {bar}")

    print(f"\n  Дикая зона (r=32..63):")
    for r in range(32, 64):
        bar = "█" * int(total_nonlin_per_round[r])
        quiet_pct = quiet_rounds[r] * 100
        print(f"    r={r:2d}: nonlin={total_nonlin_per_round[r]:5.1f} bits  quiet={quiet_pct:5.1f}%  {bar}")

    total_nonlin = np.sum(total_nonlin_per_round)
    total_quiet = np.sum(quiet_rounds > 0.5)
    print(f"\n  ИТОГО: {total_nonlin:.0f} бит дифф. carry-нелинейности")
    print(f"  Тихих раундов (>50% проб quiet): {total_quiet} из 64")

    # === 6. Бухгалтерия нового измерения ===
    print("\n" + "=" * 70)
    print("6. БУХГАЛТЕРИЯ НОВОГО ИЗМЕРЕНИЯ")
    print("=" * 70)

    struct_kernel = 256
    total_measured_nonlin = total_nonlin
    # schedule carry diff
    sched_diff = 0
    for trial_idx in range(num_samples):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W1)
        W2[0] ^= 0x80000000
        _, cm1, _ = sha256_with_carry_tracking(W1)
        _, cm2, _ = sha256_with_carry_tracking(W2)
        for r in range(16, 64):
            for add_i in ['add1', 'add2', 'add3']:
                key = (r, add_i)
                if key in cm1 and key in cm2:
                    diff = cm1[key][0] ^ cm2[key][0]
                    sched_diff += bin(diff).count('1')
    sched_diff /= num_samples

    total_effective = total_measured_nonlin + sched_diff

    print(f"\n  АКТИВ:")
    print(f"    Свобода входа (16 × 32):     512 бит")
    print(f"    Условия выхода (8 × 32):     256 бит")
    print(f"    Структурное ядро:             {struct_kernel} бит")

    print(f"\n  ПАССИВ (измеренная нелинейность при δW₀=1 бит):")
    print(f"    Carry раундов (дифф.):        {total_measured_nonlin:.0f} бит")
    print(f"    Carry schedule (дифф.):       {sched_diff:.0f} бит")
    print(f"    ИТОГО:                        {total_effective:.0f} бит")

    balance = struct_kernel - total_effective
    print(f"\n  БАЛАНС: {struct_kernel} - {total_effective:.0f} = {balance:.0f} бит")

    if balance > 0:
        print(f"\n  >>> ЯДРО НЕПУСТО: {balance:.0f} бит свободы")
        print(f"  >>> Оценка стоимости коллизии: ~2^{int((256 - balance) / 2)}")
    else:
        print(f"\n  >>> ЯДРО ПУСТО: дефицит {-balance:.0f} бит")
        print(f"  >>> SHA-256 защищена с запасом {-balance:.0f} бит сверх структурного ядра")
        print(f"  >>> Стоимость коллизии: 2^128 (birthday bound)")


if __name__ == "__main__":
    analyze_carry_skeleton(num_samples=100)
