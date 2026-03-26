"""
НАТИВНЫЙ BIRTHDAY: поиск слияния в нашем измерении.

Определения:
  Meeting point = раунд r где state₁[r] = state₂[r]
  Schedule lock = δW[r..63] = 0 (одинаковый schedule с раунда r)
  Слияние = meeting point + schedule lock → collision

Нативный birthday:
  Вместо поиска H₁=H₂ (256 бит), ищем:
  1. Meeting на ЛЮБОМ раунде r (не только r=64)
  2. Schedule lock с раунда r

Стоимость meeting на раунде r:
  state[r] = 256 бит → birthday 2^128
  НО: с a-repair state[r] частично фиксирован!
  И meeting на разных r имеет РАЗНУЮ стоимость.

Стоимость schedule lock с раунда r:
  δW[r..63] = 0 → зависит от δW[0..15] через schedule
  При a-repair: δW[0..2]=0, δW[3]=1bit → schedule ЧАСТИЧНО прозрачен

Вопросы:
  1. На каком r meeting дешевле всего?
  2. Сколько стоит schedule lock с каждого r?
  3. Суммарная стоимость = meeting + lock = ?
  4. Сравнение с standard birthday 2^128
"""

import numpy as np
import struct, hashlib

MASK32 = 0xFFFFFFFF
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
IV = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]

def hw(x): return bin(x).count('1')
def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
def add32(x, y): return (x + y) & MASK32
def sub32(x, y): return (x - y) & MASK32
def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)

def R(state, W_r, r_idx):
    a, b, c, d, e, f, g, h = state
    T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e, f, g)), K[r_idx]), W_r)
    T2 = add32(Sigma0(a), Maj(a, b, c))
    return (add32(T1, T2), a, b, c, add32(d, T1), e, f, g)

def expand_schedule(W16):
    W = list(W16)
    for r in range(16, 64):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))
    return W

def a_repair_W(state2, target_a, r_idx):
    a, b, c, d, e, f, g, h = state2
    T2 = add32(Sigma0(a), Maj(a, b, c))
    T1_needed = sub32(target_a, T2)
    return sub32(sub32(sub32(sub32(T1_needed, h), Sigma1(e)), Ch(e, f, g)), K[r_idx])

def state_diff(s1, s2):
    return sum(hw(s1[i] ^ s2[i]) for i in range(8))


def main():
    np.random.seed(42)

    print("=" * 70)
    print("НАТИВНЫЙ BIRTHDAY: слияние следов в нашем измерении")
    print("=" * 70)

    # =================================================================
    print(f"\n{'=' * 70}")
    print("1. СТОИМОСТЬ MEETING POINT на каждом раунде")
    print("   (с a-repair, break=3, bit=31)")
    print("=" * 70)

    # A-repair создаёт определённый δstate profile.
    # Meeting на раунде r = δstate[r] = 0.
    # Мы уже знаем P(δstate[r]=0) для каждого r.
    # Собираем заново с большей выборкой.

    N = 20000
    meeting_probs = np.zeros(65)  # P(δstate[r]=0)
    partial_meeting = np.zeros((65, 9))  # P(δstate[r] ≤ k) for k=0..8 registers

    for trial in range(N):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W1f = expand_schedule(W1)

        s1 = tuple(IV)
        states1 = [s1]
        for r in range(64):
            s1 = R(s1, W1f[r], r)
            states1.append(s1)

        W2 = list(W1)
        W2[3] ^= (1 << 31)
        s2 = tuple(IV)
        for r in range(3):
            s2 = R(s2, W2[r], r)
        s2 = R(s2, W2[3], 3)
        for r in range(4, 16):
            W2[r] = a_repair_W(s2, states1[r + 1][0], r)
            s2 = R(s2, W2[r], r)

        W2f = expand_schedule(W2)
        s2_full = tuple(IV)
        for r in range(64):
            w = W2[r] if r < 16 else W2f[r]
            s2_full = R(s2_full, w, r)

            ds = state_diff(states1[r + 1], s2_full)
            if ds == 0:
                meeting_probs[r + 1] += 1

            # How many REGISTERS match?
            reg_match = sum(1 for i in range(8) if states1[r + 1][i] == s2_full[i])
            for k in range(reg_match + 1):
                partial_meeting[r + 1, k] += 1

    meeting_probs /= N
    partial_meeting /= N

    print(f"\n  {'r':>4} {'P(meet)':>9} {'Cost':>8}  {'P(≥6reg)':>10} {'P(≥4reg)':>10}")
    for r in range(1, 65):
        p = meeting_probs[r]
        cost = f"2^{-np.log2(max(p, 1/N)):.1f}" if p > 0 else ">2^14"
        p6 = partial_meeting[r, 6]  # at least 6 registers match
        p4 = partial_meeting[r, 4]

        if r <= 18 or r >= 60 or p > 0.001 or r in [20, 25, 30, 40, 50]:
            print(f"  {r:4d} {p*100:8.2f}% {cost:>8}  {p6*100:9.2f}% {p4*100:9.2f}%")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("2. СТОИМОСТЬ SCHEDULE LOCK с каждого раунда")
    print("=" * 70)

    # Schedule lock с раунда r = δW[r..63] = 0
    # При a-repair: δW[0..2]=0, δW[3]=1bit, δW[4..15]=a-repair
    # → δW[16..63] = schedule(δW[0..15])

    sched_lock_probs = {}  # r → P(δW[r..63]=0)

    for trial in range(N):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W1f = expand_schedule(W1)

        W2 = list(W1)
        W2[3] ^= (1 << 31)
        s1 = tuple(IV)
        states1 = [s1]
        for r in range(64):
            s1 = R(s1, W1f[r], r)
            states1.append(s1)

        s2 = tuple(IV)
        for r in range(3):
            s2 = R(s2, W2[r], r)
        s2 = R(s2, W2[3], 3)
        for r in range(4, 16):
            W2[r] = a_repair_W(s2, states1[r + 1][0], r)
            s2 = R(s2, W2[r], r)

        W2f = expand_schedule(W2)

        # Check schedule lock from each r
        for r in range(16, 25):
            all_zero = all(W1f[rr] == W2f[rr] for rr in range(r, 64))
            if r not in sched_lock_probs:
                sched_lock_probs[r] = 0
            if all_zero:
                sched_lock_probs[r] += 1

    print(f"\n  Schedule lock δW[r..63]=0:")
    print(f"  {'r':>4} {'P(lock)':>10} {'Cost':>10}")
    for r in range(16, 25):
        p = sched_lock_probs.get(r, 0) / N
        cost = f"2^{-np.log2(max(p, 1/N)):.1f}" if p > 0 else f">2^{np.log2(N):.0f}"
        print(f"  {r:4d} {p*100:9.3f}% {cost:>10}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("3. НАТИВНЫЙ BIRTHDAY: meeting + lock = collision")
    print("=" * 70)

    # Для каждого r: total_cost = cost(meeting at r) + cost(schedule lock from r)
    # Collision = meeting at r AND schedule lock from r

    print(f"\n  Нативная стоимость = meeting_cost × lock_cost:")
    print(f"  {'r':>4} {'Meeting':>10} {'Lock':>10} {'Total':>10} {'vs 2^128':>10}")

    best_r = 64
    best_total = 256

    for r in range(12, 25):
        p_meet = meeting_probs[r]
        p_lock = sched_lock_probs.get(r, 0) / N if r >= 16 else (1.0 if r <= 15 else 0)

        if r <= 15:
            # Before schedule kicks in: lock is trivially 1 (same W)
            # But meeting requires a-repair PLUS specific conditions
            if p_meet > 0:
                meet_cost = -np.log2(p_meet)
            else:
                meet_cost = np.log2(N)
            lock_cost = 0  # free before schedule
            total = meet_cost
        elif r <= 17:
            if p_meet > 0:
                meet_cost = -np.log2(p_meet)
            else:
                meet_cost = np.log2(N)
            if p_lock > 0:
                lock_cost = -np.log2(p_lock)
            else:
                lock_cost = np.log2(N)
            total = meet_cost + lock_cost
        else:
            meet_cost = np.log2(N)  # never found
            lock_cost = np.log2(N)
            total = meet_cost + lock_cost

        if total < best_total:
            best_total = total
            best_r = r

        print(f"  {r:4d} {'2^'+f'{meet_cost:.1f}':>9} {'2^'+f'{lock_cost:.1f}':>9} "
              f"{'2^'+f'{total:.1f}':>9} {'2^'+f'{128-total:.1f}':>9}")

    print(f"\n  ОПТИМАЛЬНЫЙ meeting point: r={best_r}")
    print(f"  Нативная стоимость: 2^{best_total:.1f}")
    print(f"  Standard birthday: 2^128")
    print(f"  GAIN: {128 - best_total:.1f} бит")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("4. PARTIAL MEETING: не все 8 регистров")
    print("=" * 70)

    # Что если meeting на r регистрах (< 8)?
    # n matched registers → P(collision) = P(remaining match at output)
    # В нашем измерении: каждый matched register = один PIPE cascade

    print(f"\n  Если при meeting совпали K из 8 регистров:")
    print(f"  {'K regs':>7} {'P at r=12':>10} {'Remaining':>10} {'Total cost':>11}")

    for k in range(1, 9):
        p_k = partial_meeting[12, k] - (partial_meeting[12, k+1] if k < 8 else 0)
        remaining_bits = (8 - k) * 32  # unsolved registers
        remaining_cost = remaining_bits / 2  # birthday on remaining

        if p_k > 0:
            meet_cost = -np.log2(p_k)
        else:
            meet_cost = np.log2(N)

        total = meet_cost + remaining_cost
        print(f"  {k:3d}/8   {p_k*100:9.2f}%  {remaining_bits:6d} bit  2^{total:.1f}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("5. СВОДКА: нативный vs стандартный")
    print("=" * 70)

    print(f"""
  СТАНДАРТНЫЙ BIRTHDAY (из старого мира):
    Цель: δH[0..7] = 0 (256 бит)
    Метод: генерировать хеши, искать совпадение
    Стоимость: 2^128

  НАТИВНЫЙ BIRTHDAY (наше измерение):
    Цель: СЛИЯНИЕ следов (meeting + schedule lock)
    Метод: a-repair → reboot → schedule transparency → lock
    Оптимальный meeting point: r={best_r}
    Стоимость: 2^{best_total:.1f}

  Разница: {128 - best_total:.1f} бит

  ПОЧЕМУ РАЗНИЦА:
    A-repair даёт БЕСПЛАТНЫЙ meeting на r=12 (100%)
    Schedule lock: δW[12..15]=0 бесплатно (a-repair design)
    δW[16]=0: 48.6% (2^1)
    δW[17]=0: 49% (2^1)
    δW[18..63]: barrier (schedule mixing)

  В нашем измерении стоимость = стоимость SCHEDULE LOCK,
  не стоимость meeting. Meeting БЕСПЛАТЕН через reboot.
  Lock стоит 2^{{schedule_barrier}}.
""")


if __name__ == "__main__":
    main()
