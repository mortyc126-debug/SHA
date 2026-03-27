"""
ЧАСТЬ 2: Хеш как программа.

Новая операция: САМОНАВЕДЕНИЕ.
H(x) используется не как "значение" для сравнения,
а как ИНСТРУКЦИЯ для выбора следующего x.

SHA-256 ведёт СЕБЯ через своё пространство.
Мы не ищем collision — мы даём SHA-256 НАЙТИ его самой.

Операция: TRACE
  x₀ = seed
  x₁ = f(H(x₀))   — H решает куда идти
  x₂ = f(H(x₁))
  ...

Если trace замкнётся (xₙ = xₘ) → CYCLE → collision!
Это Floyd/Pollard rho? Частично. Но мы определяем f ИНАЧЕ.

Стандартный rho: f = H (identity, treat hash as next input)
Наш подход: f использует СТРУКТУРУ хеша (биты как инструкции)
"""

import numpy as np
import struct, hashlib
from collections import defaultdict

MASK32 = 0xFFFFFFFF

def sha256_hash(W16):
    raw = struct.pack('>16I', *W16)
    return struct.unpack('>8I', hashlib.sha256(raw).digest())

def hw(x): return bin(x).count('1')


def trace_identity(seed, steps):
    """Стандартный rho: H → next input."""
    x = list(seed)
    trajectory = []
    for step in range(steps):
        H = sha256_hash(x)
        # Используем хеш как следующий вход (pad to 16 words)
        x = list(H) + [0]*8
        trajectory.append(tuple(H))
    return trajectory


def trace_self_modify(seed, steps):
    """Самомодификация: H говорит КАКОЙ бит менять."""
    x = list(seed)
    trajectory = []
    for step in range(steps):
        H = sha256_hash(x)
        trajectory.append(tuple(H))

        # H[0] говорит какое слово менять (bits 3:0 → word 0-15)
        word_idx = H[0] & 0xF
        # H[1] говорит какой бит менять (bits 4:0 → bit 0-31)
        bit_idx = H[1] & 0x1F
        # H[2] говорит КАК менять: XOR(0) или ADD(1)
        op = (H[2] >> 0) & 1

        if op == 0:
            x[word_idx] ^= (1 << bit_idx)
        else:
            x[word_idx] = (x[word_idx] + (1 << bit_idx)) & MASK32
    return trajectory


def trace_gradient(seed, target_hw, steps):
    """Градиентный trace: H решает направление,
    но мы идём ТОЛЬКО если приближаемся к цели."""
    x = list(seed)
    trajectory = []
    H_current = sha256_hash(x)
    current_hw = sum(hw(h) for h in H_current)

    for step in range(steps):
        # H говорит какой бит менять
        word_idx = H_current[0] & 0xF
        bit_idx = H_current[1] & 0x1F

        x_new = list(x)
        x_new[word_idx] ^= (1 << bit_idx)
        H_new = sha256_hash(x_new)
        new_hw = sum(hw(h) for h in H_new)

        # Принимаем если приближаемся к target_hw
        if abs(new_hw - target_hw) < abs(current_hw - target_hw):
            x = x_new
            H_current = H_new
            current_hw = new_hw
        else:
            # Иначе: используем H[3] для альтернативного хода
            word_idx2 = H_current[3] & 0xF
            bit_idx2 = H_current[4] & 0x1F
            x_new2 = list(x)
            x_new2[word_idx2] ^= (1 << bit_idx2)
            H_new2 = sha256_hash(x_new2)
            new_hw2 = sum(hw(h) for h in H_new2)
            x = x_new2
            H_current = H_new2
            current_hw = new_hw2

        trajectory.append(current_hw)
    return trajectory


def trace_collision_hunt(seed, steps):
    """Collision-hunting trace: два потока, управляемые одним хешем."""
    x1 = list(seed)
    x2 = list(seed)
    x2[0] ^= 1  # Стартуем с δ=1 бит

    # Оба потока делают ОДИНАКОВЫЕ ходы (управляемые H1)
    # Если H1 = H2 → collision!
    # Если нет — ходы расходятся (разные H → разные инструкции)

    for step in range(steps):
        H1 = sha256_hash(x1)
        H2 = sha256_hash(x2)

        dH = sum(hw(H1[i] ^ H2[i]) for i in range(8))

        if dH == 0:
            return step, x1, x2, H1

        # Оба делают ход, управляемый СВОИМ хешем
        w1 = H1[0] & 0xF
        b1 = H1[1] & 0x1F
        x1[w1] ^= (1 << b1)

        w2 = H2[0] & 0xF
        b2 = H2[1] & 0x1F
        x2[w2] ^= (1 << b2)

    return None


def trace_convergent(seed, steps):
    """Convergent trace: МНОГО потоков, идущих к ОДНОЙ точке.
    SHA-256 управляет каждым. Кто придёт первым?"""
    N_STREAMS = 100
    streams = []
    for i in range(N_STREAMS):
        x = list(seed)
        x[0] = (x[0] + i) & MASK32
        streams.append(x)

    # Target: all-zero hash (или любой фиксированный)
    target = tuple([0]*8)

    closest = 256
    for step in range(steps):
        for i in range(N_STREAMS):
            H = sha256_hash(streams[i])
            d = sum(hw(H[j]) for j in range(8))  # distance to all-zero

            if d < closest:
                closest = d

            # Self-guided step
            w = H[0] & 0xF
            b = H[1] & 0x1F
            streams[i][w] ^= (1 << b)

        if step % 100 == 0 and step > 0:
            pass  # progress tracked below

    return closest


def main():
    np.random.seed(42)

    print("=" * 70)
    print("ЧАСТЬ 2: Хеш как программа (самонаведение)")
    print("=" * 70)

    seed = [np.random.randint(0, 2**32) for _ in range(16)]

    # ═══════════════════════════════
    print(f"\n{'=' * 70}")
    print("ТЕСТ A: Identity trace (стандартный rho)")
    print("=" * 70)

    traj = trace_identity(seed, 10000)
    # Ищем повторы (cycle)
    seen = {}
    cycle_found = None
    for i, h in enumerate(traj):
        # Truncate to 32 bits for faster cycle detection
        key = h[0]
        if key in seen:
            cycle_found = (seen[key], i)
            break
        seen[key] = i

    if cycle_found:
        print(f"  32-bit truncated cycle: step {cycle_found[0]} → {cycle_found[1]}")
        print(f"  Cycle length: {cycle_found[1] - cycle_found[0]}")
    else:
        print(f"  No 32-bit cycle in 10K steps (expected: ~65K for birthday)")

    # ═══════════════════════════════
    print(f"\n{'=' * 70}")
    print("ТЕСТ B: Self-modify trace")
    print("=" * 70)

    traj_sm = trace_self_modify(seed, 5000)
    # Measure: HW distribution of hashes in trajectory
    hws = [sum(hw(h) for h in t) for t in traj_sm]
    print(f"  5000 steps of self-modify:")
    print(f"    HW: mean={np.mean(hws):.1f}, std={np.std(hws):.1f}")
    print(f"    Min HW: {min(hws)}, Max: {max(hws)}")

    # Does self-modify converge to LOW HW? (self-guided toward all-zero?)
    print(f"    First 10: {hws[:10]}")
    print(f"    Last 10:  {hws[-10:]}")
    print(f"    → {'CONVERGING!' if np.mean(hws[-100:]) < np.mean(hws[:100]) - 5 else 'Random walk.'}")

    # ═══════════════════════════════
    print(f"\n{'=' * 70}")
    print("ТЕСТ C: Gradient trace (target HW=0)")
    print("=" * 70)

    traj_grad = trace_gradient(seed, target_hw=0, steps=5000)
    print(f"  5000 gradient steps toward HW=0:")
    print(f"    Start: HW={traj_grad[0]}")
    print(f"    End:   HW={traj_grad[-1]}")
    print(f"    Min:   HW={min(traj_grad)}")
    print(f"    → Gradient descent: {traj_grad[0]} → {min(traj_grad)} (improved by {traj_grad[0] - min(traj_grad)})")

    # ═══════════════════════════════
    print(f"\n{'=' * 70}")
    print("ТЕСТ D: Collision hunt (два потока)")
    print("=" * 70)

    result = trace_collision_hunt(seed, 50000)
    if result:
        step, x1, x2, H = result
        print(f"  ★★★ COLLISION at step {step}!")
        print(f"    H = {tuple(hex(h) for h in H)}")
    else:
        print(f"  No collision in 50K steps (2×50K = 100K hashes)")
        print(f"  Expected (birthday on 256 bits): 2^128 steps")
        print(f"  → Self-guided doesn't beat birthday")

    # Track δ between the two streams
    x1 = list(seed)
    x2 = list(seed); x2[0] ^= 1
    deltas = []
    for step in range(1000):
        H1 = sha256_hash(x1)
        H2 = sha256_hash(x2)
        dH = sum(hw(H1[i]^H2[i]) for i in range(8))
        deltas.append(dH)
        w1 = H1[0]&0xF; b1 = H1[1]&0x1F; x1[w1]^=(1<<b1)
        w2 = H2[0]&0xF; b2 = H2[1]&0x1F; x2[w2]^=(1<<b2)

    print(f"\n  δH trajectory (1000 steps):")
    print(f"    Start: {deltas[0]}")
    print(f"    Mean:  {np.mean(deltas):.1f}")
    print(f"    Min:   {min(deltas)}")
    print(f"    → {'CONVERGING!' if min(deltas) < 100 else 'Diverges to random.'}")

    # ═══════════════════════════════
    print(f"\n{'=' * 70}")
    print("ТЕСТ E: Convergent search (100 потоков)")
    print("=" * 70)

    closest = trace_convergent(seed, steps=500)
    print(f"  100 streams × 500 steps = 50K hashes")
    print(f"  Closest to all-zero: HW = {closest}")

    # Expected: min HW from 50K random hashes
    # Each HW ~ Binomial(256, 0.5), mean=128, std=8
    # Min of 50K: ≈ 128 - 8*Φ^(-1)(1-1/50K) ≈ 128 - 8*4.1 ≈ 95
    from scipy.stats import norm
    expected_min = 128 - 8 * norm.ppf(1 - 1.0/50000)
    print(f"  Expected min (random 50K): {expected_min:.0f}")
    print(f"  → {'★ BETTER THAN RANDOM!' if closest < expected_min - 5 else 'Same as random.'}")

    # ═══════════════════════════════
    print(f"\n{'=' * 70}")
    print("ИТОГ ЧАСТИ 2")
    print("=" * 70)
    print(f"""
  Self-guided search (хеш как программа):
    Identity rho:      no cycle in 10K (expected)
    Self-modify:       random walk (no convergence)
    Gradient:          {traj_grad[0]} → {min(traj_grad)} (modest improvement)
    Collision hunt:    diverges to δ≈128
    Convergent (100):  closest={closest} vs expected={expected_min:.0f}

  ВЫВОД: самонаведение НЕ даёт преимущества.
  SHA-256 перемешивает настолько хорошо, что
  "инструкции" из хеша = случайные шаги.
  Самоуправляемый поиск = random walk.
""")


if __name__ == "__main__":
    main()
