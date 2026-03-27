"""
ОТКРЫТЫЕ ПРОБЛЕМЫ МАТЕМАТИКИ В НАШЕМ ИЗМЕРЕНИИ.

Переформулируем великие проблемы в терминах SHA-256 пространства.

1. P vs NP → "Обратимость": структура прообразов
2. Collatz → "Орбиты": итерация хеш-функции
3. Riemann → "Простые": особые точки и их распределение
4. Navier-Stokes → "Турбулентность": хаос в потоке информации
5. Hodge → "Волокна": алгебраическая структура прообразов
"""

import numpy as np
import struct, hashlib
from collections import Counter, defaultdict

MASK32 = 0xFFFFFFFF

def sha256_hash(W16):
    raw = struct.pack('>16I', *W16)
    return struct.unpack('>8I', hashlib.sha256(raw).digest())

def hw(x): return bin(x).count('1')

def hash_to_input(H):
    """Treat 256-bit hash as 16-word input (pad with zeros)."""
    return list(H) + [0]*8


def main():
    np.random.seed(42)

    # ═════════════════════════════════════════════════════════
    print("=" * 70)
    print("ПРОБЛЕМА 1: P vs NP в нашем измерении")
    print("=" * 70)
    print("""
  Формулировка: Существует ли СТРУКТУРА в прообразах SHA-256,
  которая позволяет найти прообраз быстрее чем 2^256?

  В нашем измерении это = "есть ли КОРОТКИЙ ПУТЬ назад?"
  Мы знаем: вперёд = 1 шаг (polynomial).
  Назад = ??? шагов.

  Тест: для truncated SHA-256 (8 бит), сколько РЕАЛЬНО нужно
  попыток чтобы найти прообраз? Меньше чем 256 (brute force)?
""")

    # Preimage search: find x such that H(x) = target
    # For 8-bit truncated: brute force = 256 tries on average
    # Can we do better using STRUCTURE?

    def sha256_8bit(W16):
        H = sha256_hash(W16)
        return H[0] & 0xFF

    # Method 1: Brute force
    target = 42
    attempts_bf = []
    for trial in range(1000):
        for attempt in range(1, 10000):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            if sha256_8bit(W) == target:
                attempts_bf.append(attempt)
                break

    print(f"  Brute force (8-bit target):")
    print(f"    Mean attempts: {np.mean(attempts_bf):.1f}")
    print(f"    Expected: 256.0")

    # Method 2: Structured search using Conservation Law
    # (a+e)[r] = (d+h)[r+3] — can we exploit this?
    # For preimage: we know H[0] (target), need to find W[0..15]
    # H[0] = IV[0] + a[64], so a[64] = H[0] - IV[0]
    # With conservation: (a+e)[60] = (d+h)[63]
    # But we don't know e[64], so this gives 1 equation with 2 unknowns

    # Method 3: Near-target search
    # Start with random W, if H is "close" to target, modify W locally
    attempts_near = []
    for trial in range(1000):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        best_diff = 8
        for attempt in range(1, 10000):
            H_val = sha256_8bit(W)
            if H_val == target:
                attempts_near.append(attempt)
                break
            # Flip one bit guided by current distance
            diff = hw(H_val ^ target)
            w = attempt % 16
            b = (attempt * 7) % 32
            W_new = list(W)
            W_new[w] ^= (1 << b)
            H_new = sha256_8bit(W_new)
            if hw(H_new ^ target) <= diff:
                W = W_new
        else:
            attempts_near.append(10000)

    print(f"  Near-target search:")
    print(f"    Mean attempts: {np.mean(attempts_near):.1f}")
    print(f"  → {'FASTER than brute force!' if np.mean(attempts_near) < np.mean(attempts_bf) * 0.8 else 'Same as brute force'}")

    # ═════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("ПРОБЛЕМА 2: COLLATZ в нашем измерении")
    print("=" * 70)
    print("""
  Collatz: n → n/2 (even) или 3n+1 (odd). Все сходятся к 1?

  В нашем измерении: x → H(x) → H(H(x)) → ...
  Вопросы:
    A. Все ли орбиты в конце концов зацикливаются?
    B. Какова длина цикла?
    C. Сколько различных циклов существует?
    D. Есть ли "1" — универсальный аттрактор?
""")

    # Iterate: x → H(x), treat H as next input
    def iterate_hash(seed, max_steps=10000):
        x = list(seed)
        trajectory = []
        seen = {}
        for step in range(max_steps):
            H = sha256_hash(x)
            # Truncate to make cycle detection feasible
            key = H[0]  # 32-bit fingerprint
            trajectory.append(key)
            if key in seen:
                cycle_start = seen[key]
                cycle_len = step - cycle_start
                return step, cycle_start, cycle_len, trajectory
            seen[key] = step
            x = list(H) + [0]*8  # H as next input
        return max_steps, -1, -1, trajectory

    print(f"\n  A. Orbit cycling (32-bit truncated, 100 seeds):")
    cycle_lens = []
    tail_lens = []
    for trial in range(100):
        seed = [np.random.randint(0, 2**32) for _ in range(16)]
        total, start, clen, _ = iterate_hash(seed, 100000)
        if clen > 0:
            cycle_lens.append(clen)
            tail_lens.append(start)

    found = len(cycle_lens)
    print(f"    Cycles found: {found}/100 (in 100K steps)")
    if cycle_lens:
        print(f"    Cycle length: mean={np.mean(cycle_lens):.0f}, min={min(cycle_lens)}, max={max(cycle_lens)}")
        print(f"    Tail length:  mean={np.mean(tail_lens):.0f}")
        # Birthday bound for 32-bit: sqrt(2^32) ≈ 65536
        print(f"    Expected (birthday on 32 bits): ~{int(2**16)}")
    else:
        print(f"    No cycles in 100K steps (expected: need ~2^16 = 65536)")

    # B. Do different seeds converge to SAME cycle?
    if len(cycle_lens) >= 2:
        # Check if cycles share the same "loop value"
        cycle_values = set()
        for trial in range(min(10, found)):
            seed = [np.random.randint(0, 2**32) for _ in range(16)]
            np.random.seed(trial)
            seed = [np.random.randint(0, 2**32) for _ in range(16)]
            total, start, clen, traj = iterate_hash(seed, 100000)
            if clen > 0:
                cycle_val = traj[start]  # first value in cycle
                cycle_values.add(cycle_val)

        print(f"\n  B. Distinct cycle attractors (from {min(10,found)} orbits): {len(cycle_values)}")
        print(f"    → {'ONE attractor (like Collatz!)' if len(cycle_values) == 1 else f'{len(cycle_values)} attractors'}")

    # ═════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("ПРОБЛЕМА 3: RIEMANN в нашем измерении")
    print("=" * 70)
    print("""
  Riemann: нули ζ(s) лежат на линии Re(s)=1/2.
  О распределении простых чисел.

  В нашем измерении: "простые" = точки с особыми свойствами.
  Определим "простое" хеш-значение:
    H "простое" если HW(H) = простое число.

  Вопрос: как распределены "простые хеши"?
  Есть ли аналог PNT (Prime Number Theorem)?
  Есть ли "нули" — значения с аномальной плотностью?
""")

    # Generate many hashes, classify by HW
    N = 100000
    hw_counts = Counter()
    for _ in range(N):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H = sha256_hash(W)
        total_hw = sum(hw(h) for h in H)
        hw_counts[total_hw] += 1

    # Is HW distribution binomial?
    from scipy.stats import binom, norm
    mean_hw = sum(k*v for k,v in hw_counts.items()) / N
    std_hw = np.sqrt(sum((k-mean_hw)**2 * v for k,v in hw_counts.items()) / N)

    print(f"\n  HW distribution of SHA-256 hashes (N={N}):")
    print(f"    Mean: {mean_hw:.2f} (expected: 128.0)")
    print(f"    Std:  {std_hw:.2f} (expected: 8.0)")

    # "Prime" hashes: HW = prime number
    def is_prime(n):
        if n < 2: return False
        for i in range(2, int(n**0.5)+1):
            if n % i == 0: return False
        return True

    prime_hws = {k: v for k, v in hw_counts.items() if is_prime(k)}
    non_prime_hws = {k: v for k, v in hw_counts.items() if not is_prime(k)}
    prime_count = sum(prime_hws.values())
    total = sum(hw_counts.values())

    print(f"\n  'Prime' hashes (HW = prime): {prime_count}/{total} = {prime_count/total*100:.1f}%")

    # What fraction of HW values in [100..156] are prime?
    prime_density = sum(1 for k in range(100, 157) if is_prime(k)) / 57
    print(f"  Prime density in [100,156]: {prime_density*100:.1f}%")
    print(f"  Expected 'prime' hash fraction: ~{prime_density*100:.1f}%")

    # More interesting: define "SHA-256 primes" differently
    # A hash H is "prime" if it cannot be written as H = H(H(x)) for any x
    # (i.e., it's not a "composite" = result of double-hashing)
    # This is hard to test... but we can check: among iterated hashes,
    # what fraction of 32-bit truncated values are "reachable"?

    print(f"\n  Alternative: 'reachable' values under iteration")
    reachable = set()
    for trial in range(1000):
        x = [np.random.randint(0, 2**32) for _ in range(16)]
        for step in range(100):
            H = sha256_hash(x)
            reachable.add(H[0] & 0xFFFF)  # 16-bit truncated
            x = list(H) + [0]*8

    print(f"    After 1000×100 iterations: {len(reachable)}/65536 16-bit values reached")
    print(f"    Coverage: {len(reachable)/65536*100:.1f}%")
    print(f"    → {'Full coverage (ergodic)' if len(reachable) > 60000 else f'Partial: {len(reachable)} values are attractors'}")

    # ═════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("ПРОБЛЕМА 4: NAVIER-STOKES в нашем измерении")
    print("=" * 70)
    print("""
  Navier-Stokes: существует ли гладкое решение для потока жидкости?
  Или поток может развить СИНГУЛЯРНОСТЬ (бесконечную скорость)?

  В нашем измерении: информация "течёт" через раунды.
  "Скорость" = rate of change (dState/dr).
  Вопрос: может ли "скорость" стать бесконечной (сингулярность)?
  Или поток ВСЕГДА гладкий?
""")

    # "Velocity" = HW(state[r] ⊕ state[r+1]) per round
    # "Acceleration" = change in velocity
    # Singularity = acceleration → infinity?

    from functools import reduce
    def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
    def add32(x, y): return (x + y) & MASK32
    def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
    def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
    def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
    def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
    def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
    def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)

    K_const = [
        0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
        0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
        0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
        0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
        0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
        0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
        0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
        0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
    ]
    IV_c = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
            0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]

    velocities_all = []
    accelerations_all = []
    max_accel = 0

    for _ in range(1000):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        W_exp = list(W)
        for r in range(16, 64):
            W_exp.append(add32(add32(add32(sigma1(W_exp[r-2]), W_exp[r-7]), sigma0(W_exp[r-15])), W_exp[r-16]))

        a,b,c,d,e,f,g,h = IV_c
        prev_state = (a,b,c,d,e,f,g,h)
        velocities = []

        for r in range(64):
            T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K_const[r]), W_exp[r])
            T2 = add32(Sigma0(a), Maj(a,b,c))
            h,g,f,e = g,f,e,add32(d,T1)
            d,c,b,a = c,b,a,add32(T1,T2)
            curr_state = (a,b,c,d,e,f,g,h)

            v = sum(hw(prev_state[i] ^ curr_state[i]) for i in range(8))
            velocities.append(v)
            prev_state = curr_state

        velocities_all.append(velocities)

        # Acceleration
        accels = [abs(velocities[i+1] - velocities[i]) for i in range(63)]
        accelerations_all.append(accels)
        if max(accels) > max_accel:
            max_accel = max(accels)

    # Average velocity profile
    avg_vel = np.mean(velocities_all, axis=0)
    avg_acc = np.mean(accelerations_all, axis=0)

    print(f"  Velocity profile (HW of state change per round):")
    print(f"    r=0: {avg_vel[0]:.1f}")
    print(f"    r=4: {avg_vel[4]:.1f}")
    print(f"    r=8: {avg_vel[8]:.1f}")
    print(f"    r=32: {avg_vel[32]:.1f}")
    print(f"    r=63: {avg_vel[63]:.1f}")
    print(f"    Overall mean: {np.mean(avg_vel):.1f}")
    print(f"    Max velocity: {max(avg_vel):.1f}")

    print(f"\n  Acceleration (change in velocity):")
    print(f"    Mean: {np.mean(avg_acc):.1f}")
    print(f"    Max ever observed: {max_accel}")
    print(f"    → {'SINGULARITY possible!' if max_accel > 200 else f'Bounded by {max_accel} (no singularity)'}")

    # Is velocity BOUNDED?
    all_vels = [v for vl in velocities_all for v in vl]
    print(f"\n  Velocity bounds (64K samples):")
    print(f"    Min: {min(all_vels)}, Max: {max(all_vels)}")
    print(f"    Theoretical max: 256 (all bits flip)")
    print(f"    → Flow is {'BOUNDED (no blow-up)' if max(all_vels) < 200 else 'approaches maximum'}")

    # ═════════════════════════════════════════════════════════
    print(f"\n{'=' * 70}")
    print("СВОДКА")
    print("=" * 70)


if __name__ == "__main__":
    main()
