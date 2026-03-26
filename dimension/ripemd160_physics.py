"""
RIPEMD-160 через призму нового измерения.
Сравнение с SHA-256: carry-скелет, физика меток, законы сохранения.

RIPEMD-160: 80 раундов, 5 регистров, ДВЕ параллельные цепочки (left/right),
финальное слияние. 160-бит выход.
"""

import numpy as np

MASK32 = 0xFFFFFFFF

# RIPEMD-160 константы
KL = [0x00000000, 0x5A827999, 0x6ED9EBA1, 0x8F1BBCDC, 0xA953FD4E]
KR = [0x50A28BE6, 0x5C4DD124, 0x6D703EF3, 0x7A6D76E9, 0x00000000]

RL = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
      7,4,13,1,10,6,15,3,12,0,9,5,2,14,11,8,
      3,10,14,4,9,15,8,1,2,7,0,6,13,11,5,12,
      1,9,11,10,0,8,12,4,13,3,7,15,14,5,6,2,
      4,0,5,9,7,12,2,10,14,1,3,8,11,6,15,13]

RR = [5,14,7,0,9,2,11,4,13,6,15,8,1,10,3,12,
      6,11,3,7,0,13,5,10,14,15,8,12,4,9,1,2,
      15,5,1,3,7,14,6,9,11,8,12,2,10,0,4,13,
      8,6,4,1,3,11,15,0,5,12,2,13,9,7,10,14,
      12,15,10,4,1,5,8,7,6,2,13,14,0,3,9,11]

SL = [11,14,15,12,5,8,7,9,11,13,14,15,6,7,9,8,
      7,6,8,13,11,9,7,15,7,12,15,9,11,7,13,12,
      11,13,6,7,14,9,13,15,14,8,13,6,5,12,7,5,
      11,12,14,15,14,15,9,8,9,14,5,6,8,6,5,12,
      9,15,5,11,6,8,13,12,5,12,13,14,11,8,5,6]

SR = [8,9,9,11,13,15,15,5,7,7,8,11,14,14,12,6,
      9,13,15,7,12,8,9,11,7,7,12,7,6,15,13,11,
      9,7,15,11,8,6,6,14,12,13,5,14,13,13,7,5,
      15,5,8,11,14,14,6,14,6,9,12,9,12,5,15,8,
      8,5,12,9,12,5,14,6,8,13,6,5,15,13,11,11]

IV_RIPEMD = [0x67452301, 0xEFCDAB89, 0x98BADCFE, 0x10325476, 0xC3D2E1F0]

def rotl(x, n):
    return ((x << n) | (x >> (32 - n))) & MASK32

def f_ripemd(j, x, y, z):
    if j < 16:   return x ^ y ^ z
    elif j < 32: return (x & y) | (~x & z) & MASK32
    elif j < 48: return (x | ~y & MASK32) ^ z
    elif j < 64: return (x & z) | (y & ~z & MASK32)
    else:        return x ^ (y | ~z & MASK32)

def add32(x, y): return (x + y) & MASK32
def hw(x): return bin(x).count('1')

def add32_carry_count(x, y):
    result = (x + y) & MASK32
    c = 0
    count = 0
    for i in range(32):
        s = ((x >> i) & 1) + ((y >> i) & 1) + c
        c = s >> 1
        if c: count += 1
    return result, count


def ripemd160_round_by_round(W_input):
    """RIPEMD-160 с сохранением состояния на каждом раунде."""
    W = list(W_input[:16])

    # Left chain
    al, bl, cl, dl, el = IV_RIPEMD
    left_states = [(al, bl, cl, dl, el)]
    left_carries = []

    for j in range(80):
        fval = f_ripemd(j, bl, cl, dl)
        # T = a + f + W[r] + K
        t1, c1 = add32_carry_count(al, fval)
        t2, c2 = add32_carry_count(t1, W[RL[j]])
        t3, c3 = add32_carry_count(t2, KL[j // 16])
        T = rotl(t3, SL[j])
        t4, c4 = add32_carry_count(T, el)

        left_carries.append(c1 + c2 + c3 + c4)

        al = el
        el = dl
        dl = rotl(cl, 10)
        cl = bl
        bl = t4
        left_states.append((al, bl, cl, dl, el))

    # Right chain
    ar, br, cr, dr, er = IV_RIPEMD
    right_states = [(ar, br, cr, dr, er)]
    right_carries = []

    for j in range(80):
        fval = f_ripemd(79 - j, br, cr, dr)
        t1, c1 = add32_carry_count(ar, fval)
        t2, c2 = add32_carry_count(t1, W[RR[j]])
        t3, c3 = add32_carry_count(t2, KR[j // 16])
        T = rotl(t3, SR[j])
        t4, c4 = add32_carry_count(T, er)

        right_carries.append(c1 + c2 + c3 + c4)

        ar = er
        er = dr
        dr = rotl(cr, 10)
        cr = br
        br = t4
        right_states.append((ar, br, cr, dr, er))

    return left_states, right_states, left_carries, right_carries


def analyze_ripemd160(num_samples=200):
    np.random.seed(42)

    print("=" * 70)
    print("RIPEMD-160 — Физика меток в новом измерении")
    print("=" * 70)

    # === Структура ткани ===
    print("\n  СТРУКТУРА ТКАНИ RIPEMD-160:")
    print(f"    Регистров:       5 (a,b,c,d,e)")
    print(f"    Раундов:         80")
    print(f"    Параллельных цепочек: 2 (left + right)")
    print(f"    Труб/раунд:      3 (a→e, d→e→d, c→rotl→d, b→c)")
    print(f"    Узлов/раунд:     1 (b' = rotl(a + f + W + K, s) + e)")
    print(f"    Трубы/Узлы:      3:1 = 75% инертно (как SHA-256!)")
    print(f"    Слов входа:      16")
    print(f"    Schedule:        НЕТ расширения! Только перестановка индексов")
    print(f"    Горлышко:        5 регистров × 2 цепочки = rank ≤ 10")
    print(f"    Выход:           160 бит (5 × 32)")

    # === Carry-скелет ===
    print("\n" + "=" * 70)
    print("1. CARRY-СКЕЛЕТ RIPEMD-160")
    print("=" * 70)

    all_left_carries = np.zeros((num_samples, 80))
    all_right_carries = np.zeros((num_samples, 80))

    for trial in range(num_samples):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        ls, rs, lc, rc = ripemd160_round_by_round(W)
        all_left_carries[trial] = lc
        all_right_carries[trial] = rc

    mean_lc = np.mean(all_left_carries, axis=0)
    mean_rc = np.mean(all_right_carries, axis=0)

    print(f"\n  Средний carry/раунд:")
    print(f"    Left chain:  {np.mean(mean_lc):.1f} бит/раунд")
    print(f"    Right chain: {np.mean(mean_rc):.1f} бит/раунд")
    print(f"    ВСЕГО: {(np.sum(mean_lc) + np.sum(mean_rc)):.0f} бит абсолютного carry")

    # === Физика меток ===
    print("\n" + "=" * 70)
    print("2. ФИЗИКА МЕТОК RIPEMD-160")
    print("=" * 70)

    all_mass_left = np.zeros((num_samples, 81))
    all_mass_right = np.zeros((num_samples, 81))

    for trial in range(num_samples):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W1)
        W2[0] ^= 0x80000000

        ls1, rs1, _, _ = ripemd160_round_by_round(W1)
        ls2, rs2, _, _ = ripemd160_round_by_round(W2)

        for r in range(81):
            delta_l = [ls1[r][i] ^ ls2[r][i] for i in range(5)]
            delta_r = [rs1[r][i] ^ rs2[r][i] for i in range(5)]
            all_mass_left[trial, r] = sum(hw(d) for d in delta_l)
            all_mass_right[trial, r] = sum(hw(d) for d in delta_r)

    mean_ml = np.mean(all_mass_left, axis=0)
    mean_mr = np.mean(all_mass_right, axis=0)

    print(f"\n  Эволюция массы (δW₀ = 1 бит MSB):")
    print(f"\n  {'Раунд':>6} {'Left':>8} {'Right':>8} {'Total':>8}")
    for r in [0, 1, 2, 3, 4, 5, 8, 12, 16, 24, 32, 48, 64, 79, 80]:
        if r <= 80:
            print(f"  r={r:3d}  {mean_ml[r]:7.1f}  {mean_mr[r]:7.1f}  {mean_ml[r]+mean_mr[r]:7.1f}")

    theoretical_random_l = 80.0  # 5 × 32 × 0.5
    theoretical_random_total = 160.0

    sat_l = np.mean(mean_ml[20:80])
    sat_r = np.mean(mean_mr[20:80])
    print(f"\n  Аттрактор left:  {sat_l:.1f} (теор. random = {theoretical_random_l:.0f})")
    print(f"  Аттрактор right: {sat_r:.1f}")
    print(f"  Аттрактор total: {sat_l + sat_r:.1f} (теор. = {theoretical_random_total:.0f})")

    # === Дифференциальный carry ===
    print("\n" + "=" * 70)
    print("3. ДИФФЕРЕНЦИАЛЬНЫЙ CARRY RIPEMD-160")
    print("=" * 70)

    diff_carry_per_round = np.zeros(80)
    for trial in range(num_samples):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W1)
        W2[0] ^= 0x80000000

        _, _, lc1, _ = ripemd160_round_by_round(W1)
        _, _, lc2, _ = ripemd160_round_by_round(W2)

        for r in range(80):
            diff_carry_per_round[r] += abs(lc1[r] - lc2[r])

    diff_carry_per_round /= num_samples

    print(f"\n  |ΔCarry| по раундам (left chain, δW₀=1 бит):")
    for r in [0, 1, 2, 3, 4, 5, 8, 16, 32, 48, 64, 79]:
        bar = "█" * int(diff_carry_per_round[r])
        print(f"    r={r:2d}: {diff_carry_per_round[r]:5.1f}  {bar}")

    # === Скорость диффузии ===
    print("\n" + "=" * 70)
    print("4. СКОРОСТЬ ДИФФУЗИИ: RIPEMD-160 vs SHA-256")
    print("=" * 70)

    # Сколько раундов до 90% насыщения?
    target_90 = theoretical_random_l * 0.9
    round_90_l = next((r for r in range(81) if mean_ml[r] >= target_90), 80)
    target_50 = theoretical_random_l * 0.5
    round_50_l = next((r for r in range(81) if mean_ml[r] >= target_50), 80)

    print(f"\n  RIPEMD-160 (left chain):")
    print(f"    50% насыщения: раунд {round_50_l}")
    print(f"    90% насыщения: раунд {round_90_l}")
    print(f"    Регистров: 5, раундов: 80")
    print(f"\n  SHA-256 (из предыдущего эксперимента):")
    print(f"    50% насыщения: раунд ~3")
    print(f"    90% насыщения: раунд ~5")
    print(f"    Регистров: 8, раундов: 64")

    # === Бухгалтерия ===
    print("\n" + "=" * 70)
    print("5. БУХГАЛТЕРИЯ НОВОГО ИЗМЕРЕНИЯ: RIPEMD-160")
    print("=" * 70)

    struct_kernel = 16 * 32 - 160  # 512 - 160 = 352
    total_abs_carry = np.sum(mean_lc) + np.sum(mean_rc)

    print(f"\n  АКТИВ:")
    print(f"    Свобода входа (16 × 32):     512 бит")
    print(f"    Условия выхода (5 × 32):     160 бит")
    print(f"    Структурное ядро:             {struct_kernel} бит")
    print(f"\n  ПАССИВ:")
    print(f"    Абсолютный carry (left+right): {total_abs_carry:.0f} бит")
    print(f"    Дифф. carry (оценка):         ~{total_abs_carry * 0.8:.0f} бит")
    print(f"\n  БАЛАНС: {struct_kernel} - {total_abs_carry:.0f} = {struct_kernel - total_abs_carry:.0f} бит")
    print(f"  Запас прочности: ×{total_abs_carry / struct_kernel:.1f}")
    print(f"\n  Для сравнения:")
    print(f"    SHA-256:    запас ×31.7")
    print(f"    RIPEMD-160: запас ×{total_abs_carry / struct_kernel:.1f}")
    print(f"    Birthday:   SHA-256 = 2^128, RIPEMD-160 = 2^80")


if __name__ == "__main__":
    analyze_ripemd160(num_samples=200)
