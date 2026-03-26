"""
ПОЧЕМУ ЦИКЛЫ КОРОЧЕ? Структура орбит SHA-256.

Из topology.py: λ/expected = 0.62 (циклы 1.6× короче random).
Вопросы:
  1. Масштабируется ли эффект? (2^8 → 2^24: ratio stable?)
  2. Зависит ли от padding? (zero, mirror, const)
  3. Зависит ли от КАКОЙ позиции берём? (H[0] vs H[7])
  4. Связан ли с PIPE/NODE структурой?
"""

import numpy as np

MASK32 = 0xFFFFFFFF

def sha256_compress_one_word(x, n_bits, position=0):
    """SHA-256 iterate: x → H[position] mod 2^n_bits."""
    import struct, hashlib
    mask = (1 << n_bits) - 1
    W = [0] * 16
    W[0] = x & mask
    h = hashlib.sha256(struct.pack('>16I', *W)).digest()
    words = struct.unpack('>8I', h)
    return words[position] & mask


def floyd(x0, n_bits, position=0):
    f = lambda x: sha256_compress_one_word(x, n_bits, position)
    tortoise = f(x0)
    hare = f(f(x0))
    steps = 1
    while tortoise != hare:
        tortoise = f(tortoise)
        hare = f(f(hare))
        steps += 1
        if steps > 3 * (2**n_bits):
            return None, None
    # Cycle start
    tortoise = x0
    mu = 0
    while tortoise != hare:
        tortoise = f(tortoise)
        hare = f(hare)
        mu += 1
    # Cycle length
    lam = 1
    hare = f(tortoise)
    while tortoise != hare:
        hare = f(hare)
        lam += 1
    return mu, lam


def main():
    np.random.seed(42)

    print("=" * 70)
    print("СТРУКТУРА ЦИКЛОВ: масштабирование и зависимости")
    print("=" * 70)

    # =================================================================
    print(f"\n{'=' * 70}")
    print("1. МАСШТАБИРОВАНИЕ: ratio λ/√N по размеру пространства")
    print("=" * 70)

    print(f"\n  {'n_bits':>6} {'N':>10} {'√N':>8} {'μ mean':>8} {'λ mean':>8} {'μ/√N':>7} {'λ/√N':>7} {'Status':>10}")

    scaling_ratios = []

    for n_bits in [8, 10, 12, 14, 16, 18, 20]:
        N = 2**n_bits
        expected = np.sqrt(np.pi * N / 8)

        mus = []
        lams = []
        n_trials = min(30, max(3, 2**(n_bits // 3)))

        for trial in range(n_trials):
            x0 = np.random.randint(0, 2**n_bits)
            mu, lam = floyd(x0, n_bits, position=0)
            if mu is not None:
                mus.append(mu)
                lams.append(lam)

        if mus:
            r_mu = np.mean(mus) / expected
            r_lam = np.mean(lams) / expected
            scaling_ratios.append((n_bits, r_mu, r_lam))

            status = "SHORTER" if r_lam < 0.7 else "LONGER" if r_lam > 1.3 else "~random"
            print(f"  {n_bits:>6} {N:>10} {int(expected):>8} {np.mean(mus):>7.0f} {np.mean(lams):>7.0f} "
                  f"{r_mu:>6.2f} {r_lam:>6.2f} {status:>10}")

    # Trend
    if len(scaling_ratios) > 2:
        ns = [s[0] for s in scaling_ratios]
        rs = [s[2] for s in scaling_ratios]
        trend = np.polyfit(ns, rs, 1)
        print(f"\n  λ/√N trend: slope = {trend[0]:+.4f} per bit")
        print(f"  → {'CONVERGING to random (slope→0)' if abs(trend[0]) < 0.01 else 'PERSISTENT structure' if trend[0] < -0.01 else 'DIVERGING'}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("2. ЗАВИСИМОСТЬ ОТ ВЫХОДНОЙ ПОЗИЦИИ (H[0] vs H[7])")
    print("=" * 70)

    n_bits = 16
    for pos in range(8):
        mus = []
        lams = []
        for trial in range(20):
            x0 = np.random.randint(0, 2**n_bits)
            mu, lam = floyd(x0, n_bits, position=pos)
            if mu is not None:
                mus.append(mu)
                lams.append(lam)

        expected = np.sqrt(np.pi * 2**n_bits / 8)
        types = ['NODE', 'PIPE', 'PIPE', 'PIPE', 'NODE', 'PIPE', 'PIPE', 'PIPE']

        if mus:
            r_lam = np.mean(lams) / expected
            marker = " ★" if r_lam < 0.5 else ""
            print(f"    H[{pos}] ({types[pos]}): λ={np.mean(lams):.0f}, λ/√N={r_lam:.2f}{marker}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("3. ПОЛНЫЙ ПЕРЕБОР: все циклы в 2^16")
    print("=" * 70)

    # Найти ВСЕ циклы в пространстве 2^16
    n_bits = 16
    N = 2**n_bits

    # Compute the full function table
    print(f"\n  Computing full function table for 2^{n_bits}...")
    f_table = np.zeros(N, dtype=np.int32)
    for x in range(N):
        f_table[x] = sha256_compress_one_word(x, n_bits, 0)

    # Find all cycles using visited flags
    visited = np.zeros(N, dtype=np.uint8)
    cycles = []
    cycle_members = set()

    for start in range(N):
        if visited[start]:
            continue

        # Trace orbit
        path = []
        x = start
        while not visited[x]:
            visited[x] = 1
            path.append(x)
            x = f_table[x]

        # x is now either in a known cycle or we found a new one
        if x in cycle_members:
            continue  # connects to known cycle

        # Find cycle containing x (if x is in path)
        if x in path:
            cycle_start = path.index(x)
            cycle = path[cycle_start:]
            cycles.append(cycle)
            for c in cycle:
                cycle_members.add(c)

    print(f"  Total cycles found: {len(cycles)}")
    cycle_lengths = [len(c) for c in cycles]
    print(f"  Cycle lengths: {sorted(cycle_lengths)[:20]}{'...' if len(cycle_lengths) > 20 else ''}")
    print(f"  Min: {min(cycle_lengths)}, Max: {max(cycle_lengths)}, Mean: {np.mean(cycle_lengths):.1f}")
    print(f"  Total cycle members: {len(cycle_members)} / {N} ({len(cycle_members)/N*100:.1f}%)")

    # Expected for random function: √(πN/8) average cycle, ~1 main cycle
    print(f"  Expected (random): ~{int(np.sqrt(np.pi*N/8))} mean cycle length")

    # Fixed points
    fps = [c[0] for c in cycles if len(c) == 1]
    print(f"  Fixed points (cycle length 1): {len(fps)}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("4. CYCLE DISTRIBUTION: random vs SHA-256")
    print("=" * 70)

    # For random function on N elements:
    # Expected number of cycles: ~ln(N)/2 ≈ 5.5 for N=2^16
    # Expected max cycle: ~0.78√N ≈ 200

    expected_num_cycles = np.log(N) / 2
    expected_max_cycle = 0.78 * np.sqrt(N)

    print(f"\n  SHA-256 on 2^{n_bits}:")
    print(f"    Cycles: {len(cycles)} (expected random: {expected_num_cycles:.1f})")
    print(f"    Max cycle: {max(cycle_lengths)} (expected: {expected_max_cycle:.0f})")
    print(f"    Ratio cycles: {len(cycles) / expected_num_cycles:.2f}×")
    print(f"    Ratio max: {max(cycle_lengths) / expected_max_cycle:.2f}×")

    if len(cycles) > expected_num_cycles * 1.5:
        print(f"    → MORE CYCLES than random (fragmented orbits)")
    elif len(cycles) < expected_num_cycles * 0.5:
        print(f"    → FEWER CYCLES (more connected)")
    else:
        print(f"    → ~random cycle count")

    if max(cycle_lengths) < expected_max_cycle * 0.5:
        print(f"    → SHORTER max cycle (no dominant cycle)")
    elif max(cycle_lengths) > expected_max_cycle * 1.5:
        print(f"    → LONGER max cycle (dominant cycle exists)")


if __name__ == "__main__":
    main()
