"""
ECDSA secp256k1 через призму нового измерения.

secp256k1: y² = x³ + 7 над F_p
p = 2²⁵⁶ - 2³² - 977

Безопасность = дискретный логарифм на эллиптической кривой.
"Collision" = найти k по k·G = Q (ECDLP).

Анализируем:
  - Структура ткани (point doubling + addition)
  - Carry-скелет модулярной арифметики F_p
  - Физика меток: как δ распространяется при scalar multiplication
"""

import numpy as np

# secp256k1 параметры
P = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
Gx = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
Gy = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8
A_COEFF = 0
B_COEFF = 7

def hw256(x):
    """Hamming weight of a 256-bit number."""
    return bin(x).count('1')

def mod_inv(a, p):
    """Modular inverse using extended Euclidean."""
    if a == 0:
        return 0
    g, x, _ = extended_gcd(a % p, p)
    if g != 1:
        return None
    return x % p

def extended_gcd(a, b):
    if a == 0:
        return b, 0, 1
    g, x, y = extended_gcd(b % a, a)
    return g, y - (b // a) * x, x

def point_add(x1, y1, x2, y2, p):
    """Elliptic curve point addition with carry tracking."""
    if x1 is None: return x2, y2, 0
    if x2 is None: return x1, y1, 0

    if x1 == x2 and y1 == y2:
        # Point doubling
        num = (3 * x1 * x1 + A_COEFF) % p
        den = (2 * y1) % p
    elif x1 == x2:
        return None, None, 0  # point at infinity
    else:
        num = (y2 - y1) % p
        den = (x2 - x1) % p

    inv = mod_inv(den, p)
    if inv is None:
        return None, None, 0

    lam = (num * inv) % p
    x3 = (lam * lam - x1 - x2) % p
    y3 = (lam * (x1 - x3) - y1) % p

    # "Carry" в модулярной арифметике = число бит,
    # отличающихся от "наивного" (без mod) результата
    # Приближение: HW промежуточных значений
    carry_estimate = hw256(lam) + hw256(x3) + hw256(y3)

    return x3, y3, carry_estimate

def scalar_mult(k, x, y, p):
    """Double-and-add scalar multiplication with per-step tracking."""
    rx, ry = None, None
    steps = []
    bits = bin(k)[2:]

    for i, bit in enumerate(bits):
        # Double
        rx, ry, carry_d = point_add(rx, ry, rx, ry, p)
        if bit == '1':
            # Add
            rx, ry, carry_a = point_add(rx, ry, x, y, p)
        else:
            carry_a = 0

        if rx is not None:
            steps.append({
                'step': i,
                'bit': int(bit),
                'x_hw': hw256(rx),
                'y_hw': hw256(ry),
                'carry_double': carry_d,
                'carry_add': carry_a,
            })

    return rx, ry, steps


def analyze_secp256k1(num_samples=50):
    np.random.seed(42)

    print("=" * 70)
    print("ECDSA secp256k1 — Физика в новом измерении")
    print("=" * 70)

    # === Структура ткани ===
    print("\n  СТРУКТУРА ТКАНИ secp256k1:")
    print(f"    Объект:          Точка (x, y) на кривой y²=x³+7")
    print(f"    Размер поля:     256 бит (p ≈ 2²⁵⁶)")
    print(f"    Операция:        Скалярное умножение k·G")
    print(f"    'Раундов':       256 (double-and-add, по числу бит k)")
    print(f"    Операции/раунд:  1 doubling + 0-1 addition")
    print(f"    На каждую операцию: 1 mod_inv + 3 mod_mul + mod_sub")
    print(f"    Schedule:        НЕТ (биты k — прямой вход)")
    print(f"    Горлышко:        2 координаты × 256 бит = 512 бит состояния")
    print(f"    Вход:            256 бит (скаляр k)")
    print(f"    Выход:           512 бит (точка Q = k·G)")

    # === 1. Carry-скелет ===
    print("\n" + "=" * 70)
    print("1. CARRY-СКЕЛЕТ secp256k1 (модулярная арифметика)")
    print("=" * 70)

    all_carries = []
    all_hw_x = []
    all_hw_y = []

    for trial in range(num_samples):
        k = int.from_bytes(np.random.bytes(32), 'big') % (N - 2) + 2
        _, _, steps = scalar_mult(k, Gx, Gy, P)

        carries = [s['carry_double'] + s['carry_add'] for s in steps]
        hw_xs = [s['x_hw'] for s in steps]
        hw_ys = [s['y_hw'] for s in steps]

        all_carries.append(carries)
        all_hw_x.append(hw_xs)
        all_hw_y.append(hw_ys)

    # Выровняем по длине (все ~256 шагов)
    min_len = min(len(c) for c in all_carries)
    carries_arr = np.array([c[:min_len] for c in all_carries])
    hw_x_arr = np.array([h[:min_len] for h in all_hw_x])
    hw_y_arr = np.array([h[:min_len] for h in all_hw_y])

    mean_carry = np.mean(carries_arr, axis=0)
    mean_hw_x = np.mean(hw_x_arr, axis=0)
    mean_hw_y = np.mean(hw_y_arr, axis=0)

    print(f"\n  Среднее HW промежуточных значений (≈carry) по шагам:")
    print(f"  {'Шаг':>6} {'Carry':>8} {'HW(x)':>8} {'HW(y)':>8}")
    for s in [0, 1, 2, 5, 10, 20, 50, 100, 150, 200, min_len-1]:
        if s < min_len:
            print(f"  s={s:3d}  {mean_carry[s]:7.1f}  {mean_hw_x[s]:7.1f}  {mean_hw_y[s]:7.1f}")

    total_carry = np.sum(mean_carry)
    print(f"\n  ВСЕГО carry за scalar_mult: {total_carry:.0f} бит")

    # === 2. Физика меток: δk → δQ ===
    print("\n" + "=" * 70)
    print("2. ФИЗИКА МЕТОК: как δk влияет на δQ")
    print("=" * 70)

    print(f"\n  Изменяем 1 бит в скаляре k → смотрим δ(Q)")

    bit_positions = [0, 1, 7, 15, 31, 63, 127, 255]
    for bit_pos in bit_positions:
        diffs_x = []
        diffs_y = []
        for trial in range(min(num_samples, 30)):
            k1 = int.from_bytes(np.random.bytes(32), 'big') % (N - 2) + 2
            k2 = k1 ^ (1 << bit_pos)
            if k2 <= 0 or k2 >= N:
                continue

            x1, y1, _ = scalar_mult(k1, Gx, Gy, P)
            x2, y2, _ = scalar_mult(k2, Gx, Gy, P)

            if x1 is not None and x2 is not None:
                dx = x1 ^ x2 if (x1 is not None and x2 is not None) else 0
                dy = y1 ^ y2 if (y1 is not None and y2 is not None) else 0
                diffs_x.append(hw256(dx))
                diffs_y.append(hw256(dy))

        if diffs_x:
            print(f"  δk бит {bit_pos:3d}: HW(δx)={np.mean(diffs_x):.1f}±{np.std(diffs_x):.1f}, HW(δy)={np.mean(diffs_y):.1f}±{np.std(diffs_y):.1f}")
        else:
            print(f"  δk бит {bit_pos:3d}: недостаточно сэмплов")

    # === 3. Сравнение трёх систем ===
    print("\n" + "=" * 70)
    print("3. СРАВНИТЕЛЬНАЯ ТАБЛИЦА: SHA-256 vs RIPEMD-160 vs secp256k1")
    print("=" * 70)

    print(f"""
  {'Параметр':<35} {'SHA-256':>12} {'RIPEMD-160':>12} {'secp256k1':>12}
  {'-'*35} {'-'*12} {'-'*12} {'-'*12}
  {'Тип':<35} {'хеш':>12} {'хеш':>12} {'ECDLP':>12}
  {'Раундов/шагов':<35} {'64':>12} {'80×2=160':>12} {'~256':>12}
  {'Регистров':<35} {'8':>12} {'5×2=10':>12} {'2 (x,y)':>12}
  {'Бит состояния':<35} {'256':>12} {'160×2=320':>12} {'512':>12}
  {'Бит входа':<35} {'512':>12} {'512':>12} {'256':>12}
  {'Бит выхода':<35} {'256':>12} {'160':>12} {'512':>12}
  {'Структурное ядро':<35} {'256':>12} {'352':>12} {'0*':>12}
  {'Труб (% инертных)':<35} {'75%':>12} {'60%':>12} {'0%':>12}
  {'Birthday bound':<35} {'2^128':>12} {'2^80':>12} {'2^128':>12}
  {'Тип нелинейности':<35} {'carry+Ch+Maj':>12} {'carry+f':>12} {'mod_inv':>12}

  * secp256k1: структурное ядро = 0, потому что вход (256) < выход (512).
    Безопасность — другая: не collision, а discrete log.
""")

    # === 4. Ткань secp256k1 в нашем измерении ===
    print("=" * 70)
    print("4. ТКАНЬ secp256k1 — ФУНДАМЕНТАЛЬНОЕ ОТЛИЧИЕ ОТ ХЕШЕЙ")
    print("=" * 70)

    print(f"""
  SHA-256 / RIPEMD-160 (хеши):
    Ткань = ОДНА структура, все W проходят через ОДНУ цепочку раундов.
    Коллизия = два входа, один выход.
    Барьер = нелинейность carry в фиксированной ткани.

  secp256k1 (ECDSA):
    Ткань = РАЗНАЯ для каждого k! Биты k определяют КАКИЕ операции выполняются.
    Бит=1: doubling + addition. Бит=0: только doubling.
    Ткань не фиксирована — она ЗАВИСИТ от входа.
    ECDLP = найти k по результату, когда СТРУКТУРА ТКАНИ неизвестна.

  В нашем измерении:
    Хеш:      Ткань фиксирована, ищем δ_in в ядре → carry-стена
    ECDSA:    Ткань переменная, ищем САМУ ТКАНЬ → комбинаторная стена

  Это РАЗНЫЕ типы защиты:
    Хеш:      термодинамический барьер (энтропия)
    ECDSA:    комбинаторный барьер (экспоненциально много тканей)
""")

    # === 5. Аттрактор масса ===
    print("=" * 70)
    print("5. АТТРАКТОРЫ: МАССА НА ВЫХОДЕ")
    print("=" * 70)

    hw_final_x = []
    hw_final_y = []
    for trial in range(num_samples):
        k = int.from_bytes(np.random.bytes(32), 'big') % (N - 2) + 2
        x, y, _ = scalar_mult(k, Gx, Gy, P)
        if x is not None:
            hw_final_x.append(hw256(x))
            hw_final_y.append(hw256(y))

    print(f"\n  HW финальной точки Q = k·G:")
    print(f"    HW(x): {np.mean(hw_final_x):.1f} ± {np.std(hw_final_x):.1f}  (теор. random = 128)")
    print(f"    HW(y): {np.mean(hw_final_y):.1f} ± {np.std(hw_final_y):.1f}  (теор. random = 128)")
    print(f"\n  Для сравнения:")
    print(f"    SHA-256 H:   128.3 ± 0.7  (аттрактор)")
    print(f"    secp256k1 Q: {np.mean(hw_final_x):.1f} ± {np.std(hw_final_x):.1f}")


if __name__ == "__main__":
    analyze_secp256k1(num_samples=50)
