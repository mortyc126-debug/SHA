"""
ТОПОЛОГИЯ ТКАНИ: циклы, неподвижные точки, орбиты.

В нашем измерении ткань — отображение F: вход → выход.
Если выход подать на вход (итерация): W → H(W) → H(H(W)) → ...

Вопросы:
  1. Есть ли НЕПОДВИЖНЫЕ ТОЧКИ: H(W) = W? (F(x) = x)
  2. Есть ли ЦИКЛЫ: H^k(W) = W для некоторого k?
  3. Длина цикла = ?
  4. Как быстро орбита сходится к циклу?
  5. Цикл = collision! (W и H^{k-1}(W) дают одинаковый H)

В нашем измерении: цикл — ЗАМКНУТЫЙ ПУТЬ через ткань.
Неподвижная точка — СТАЦИОНАРНЫЙ СЛЕД (ткань не меняет его).

Для итерации: нужно "согласовать" вход и выход.
Вход = 16 позиций (512 бит). Выход = 8 позиций (256 бит).
→ нельзя напрямую подать выход на вход (размеры разные).

Решение: padding. W = (H, 0x80000000, 0, ..., 0, 256) — стандартный padding.
Или: W[0..7] = H[0..7], W[8..15] = constant.
"""

import numpy as np
import struct, hashlib

MASK32 = 0xFFFFFFFF

def sha256_words(W16):
    return struct.unpack('>8I', hashlib.sha256(struct.pack('>16I', *W16)).digest())

def hw(x): return bin(x).count('1')


def h_to_w(H, padding_style='zero'):
    """Convert 8-word hash to 16-word message for iteration."""
    if padding_style == 'zero':
        return list(H) + [0] * 8
    elif padding_style == 'mirror':
        return list(H) + list(H)
    elif padding_style == 'const':
        return list(H) + [0x12345678] * 8
    elif padding_style == 'sha_pad':
        # Standard SHA-256 padding for 256-bit message
        return list(H) + [0x80000000, 0, 0, 0, 0, 0, 0, 0x00000100]


def iterate_sha(W_init, n_steps, padding='zero'):
    """Iterate: W → H(W) → H(H(W)) → ..."""
    W = list(W_init)
    orbit = [tuple(W[:8])]  # track first 8 words (hash part)

    for step in range(n_steps):
        H = sha256_words(W)
        W = h_to_w(H, padding)
        orbit.append(H)

    return orbit


def main():
    np.random.seed(42)

    print("=" * 70)
    print("ТОПОЛОГИЯ ТКАНИ: циклы и орбиты")
    print("=" * 70)

    # =================================================================
    print(f"\n{'=' * 70}")
    print("1. ОРБИТЫ: итерация SHA-256 (small scale)")
    print("=" * 70)

    # Итерируем с маленьким пространством: только 16 бит
    # W[0] = x (16 бит), W[1..15] = 0
    # H(W) → берём H[0] mod 2^16 → новый x

    def small_iterate(x, n_bits=16):
        mask = (1 << n_bits) - 1
        W = [x & mask] + [0] * 15
        H = sha256_words(W)
        return H[0] & mask

    # Floyd cycle detection
    def floyd(x0, n_bits=16):
        # Phase 1: find meeting point
        tortoise = small_iterate(x0, n_bits)
        hare = small_iterate(small_iterate(x0, n_bits), n_bits)
        steps = 1
        while tortoise != hare:
            tortoise = small_iterate(tortoise, n_bits)
            hare = small_iterate(small_iterate(hare, n_bits), n_bits)
            steps += 1
            if steps > 2**n_bits:
                return None, None, steps

        # Phase 2: find cycle start
        tortoise = x0
        mu = 0
        while tortoise != hare:
            tortoise = small_iterate(tortoise, n_bits)
            hare = small_iterate(hare, n_bits)
            mu += 1

        # Phase 3: find cycle length
        lam = 1
        hare = small_iterate(tortoise, n_bits)
        while tortoise != hare:
            hare = small_iterate(hare, n_bits)
            lam += 1

        return mu, lam, steps

    for n_bits in [8, 12, 16, 20]:
        print(f"\n  Пространство 2^{n_bits} = {2**n_bits}:")

        mus = []
        lams = []
        for trial in range(min(10, 2**(n_bits//2))):
            x0 = np.random.randint(0, 2**n_bits)
            mu, lam, steps = floyd(x0, n_bits)
            if mu is not None:
                mus.append(mu)
                lams.append(lam)

        if mus:
            print(f"    Tail (μ): mean={np.mean(mus):.0f}, range={min(mus)}-{max(mus)}")
            print(f"    Cycle (λ): mean={np.mean(lams):.0f}, range={min(lams)}-{max(lams)}")
            print(f"    Expected (random function): μ ≈ λ ≈ √(πN/8) ≈ {int(np.sqrt(np.pi * 2**n_bits / 8))}")

            # Отклонение от random?
            expected = np.sqrt(np.pi * 2**n_bits / 8)
            ratio_mu = np.mean(mus) / expected
            ratio_lam = np.mean(lams) / expected
            print(f"    Ratio μ/expected: {ratio_mu:.2f}")
            print(f"    Ratio λ/expected: {ratio_lam:.2f}")
            print(f"    → {'RANDOM-LIKE' if 0.5 < ratio_mu < 2.0 and 0.5 < ratio_lam < 2.0 else 'STRUCTURED!'}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("2. НЕПОДВИЖНЫЕ ТОЧКИ: H(W) = W?")
    print("=" * 70)

    # На малом пространстве: полный перебор
    for n_bits in [8, 12, 16]:
        fixed_points = []
        for x in range(2**n_bits):
            y = small_iterate(x, n_bits)
            if y == x:
                fixed_points.append(x)

        expected_fp = 1  # random function: ~1 fixed point on average
        print(f"\n  2^{n_bits}: {len(fixed_points)} fixed points (expected ≈{expected_fp})")
        if fixed_points and len(fixed_points) <= 10:
            print(f"    Values: {[hex(fp) for fp in fixed_points]}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("3. COLLISION ИЗ ЦИКЛА: ρ-метод в нашем измерении")
    print("=" * 70)

    # Цикл длины λ с хвостом μ:
    # x₀ → x₁ → ... → x_μ → ... → x_{μ+λ-1} → x_μ (цикл)
    # Collision: x_{μ-1} и x_{μ+λ-1} оба отображаются в x_μ.
    # Если x_{μ-1} ≠ x_{μ+λ-1}: COLLISION!

    # Стоимость: μ + λ шагов (нахождение цикла).
    # Для random function: μ + λ ≈ √(πN/4) ≈ √N.

    print(f"\n  ρ-метод (Floyd) для поиска collision:")

    for n_bits in [16, 20, 24]:
        x0 = np.random.randint(0, 2**n_bits)
        mu, lam, steps = floyd(x0, n_bits)

        if mu is not None:
            # Find the actual collision pair
            # x_{mu-1} and x_{mu+lam-1} map to x_mu
            x = x0
            for i in range(mu - 1):
                x = small_iterate(x, n_bits)
            x_before_cycle = x

            x2 = x0
            for i in range(mu + lam - 1):
                x2 = small_iterate(x2, n_bits)
            x_lam_before = x2

            y1 = small_iterate(x_before_cycle, n_bits)
            y2 = small_iterate(x_lam_before, n_bits)

            is_collision = (y1 == y2) and (x_before_cycle != x_lam_before)

            expected_cost = int(np.sqrt(np.pi * 2**n_bits / 4))

            print(f"\n  2^{n_bits}: μ={mu}, λ={lam}, total={mu+lam}")
            print(f"    Expected √N cost: {expected_cost}")
            print(f"    Ratio: {(mu+lam)/expected_cost:.2f}")
            print(f"    Collision: {is_collision} (x={hex(x_before_cycle)}, x'={hex(x_lam_before)}, H→{hex(y1)})")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("4. НАТИВНАЯ СТОИМОСТЬ ρ-COLLISION")
    print("=" * 70)

    print(f"""
  ρ-метод в нашем измерении:
    Итерация = СВЁРТКА (fold) ткани на себя.
    Цикл = ЗАМКНУТЫЙ ПУТЬ в ткани.
    Collision = две РАЗНЫЕ точки пути ведут в одну.

  Стоимость:
    N = размер пространства итерации.
    Для SHA-256 с padding: N = 2^256 (hash space).
    ρ-cost = √N = 2^128.

    ρ-метод = overlay (birthday) В НАШЕМ ИЗМЕРЕНИИ:
    оба дают √N = C^4 = 2^128.

  НО: ρ-метод использует O(1) ПАМЯТЬ (vs O(N) для birthday).
  В нашем измерении: O(1) = 1 позиция хранения vs N позиций.

  ρ-метод = НАИБОЛЕЕ ЭКОНОМНАЯ форма overlay.
  Стоимость та же (C^4), но РЕСУРСЫ меньше.
""")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("5. СТРУКТУРА ОРБИТ: random или нет?")
    print("=" * 70)

    # Если SHA-256 ≠ random function → μ, λ могут отличаться.
    # Проверяем: ratio μ/√N и λ/√N — если ≠ 1, есть структура.

    ratios_mu = []
    ratios_lam = []

    for trial in range(50):
        x0 = np.random.randint(0, 2**16)
        mu, lam, _ = floyd(x0, 16)
        if mu is not None:
            expected = np.sqrt(np.pi * 2**16 / 8)
            ratios_mu.append(mu / expected)
            ratios_lam.append(lam / expected)

    print(f"\n  50 орбит в 2^16 пространстве:")
    print(f"    μ/expected: mean={np.mean(ratios_mu):.2f}, std={np.std(ratios_mu):.2f}")
    print(f"    λ/expected: mean={np.mean(ratios_lam):.2f}, std={np.std(ratios_lam):.2f}")

    # Для random function: mean ratio ≈ 1.0
    # Если < 1: orbits SHORTER (easier collision)
    # Если > 1: orbits LONGER (harder collision)

    ratio_avg = (np.mean(ratios_mu) + np.mean(ratios_lam)) / 2

    print(f"\n    Average ratio: {ratio_avg:.2f}")
    if ratio_avg < 0.8:
        gain = 1 / ratio_avg
        print(f"    → ORBITS SHORTER! Collision {gain:.1f}× easier than random")
        print(f"    → Effective cost: 2^{128 - np.log2(gain):.1f} instead of 2^128")
    elif ratio_avg > 1.2:
        print(f"    → ORBITS LONGER! Collision harder than random")
    else:
        print(f"    → RANDOM-LIKE orbits. No structural advantage.")


if __name__ == "__main__":
    main()
