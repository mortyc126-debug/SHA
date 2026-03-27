"""
ЦИКЛОВАЯ АНОМАЛИЯ: 482 значений на циклах vs ожидаемые 160 (3×).

Единственное отклонение SHA-256 от random function в полной карте.

Вопросы:
  A. Воспроизводится на других seed_base?
  B. Зависит ли от truncation bits?
  C. Сравнение с 100 random functions — насколько аномально?
  D. Если реально — ПОЧЕМУ?
"""

import numpy as np
import struct, hashlib
from collections import Counter

MASK32 = 0xFFFFFFFF

def sha256_hash(W16):
    raw = struct.pack('>16I', *W16)
    return struct.unpack('>8I', hashlib.sha256(raw).digest())

def hw(x): return bin(x).count('1')


def analyze_map(fwd_map, N):
    """Find cycles and compute statistics for a function map."""
    visited = set()
    cycles = []

    for start in range(N):
        if start in visited:
            continue
        path = []
        seen = {}
        val = start
        while val not in seen:
            if val in visited:
                break
            seen[val] = len(path)
            path.append(val)
            val = fwd_map.get(val, 0)
        if val in seen:
            cs = seen[val]
            cycle = path[cs:]
            cycles.append(cycle)
            visited.update(cycle)
        visited.update(path)

    image = set(fwd_map.values())
    total_on_cycles = sum(len(c) for c in cycles)
    cycle_lens = sorted([len(c) for c in cycles])

    return {
        'image_size': len(image),
        'n_cycles': len(cycles),
        'total_on_cycles': total_on_cycles,
        'cycle_lens': cycle_lens,
        'max_cycle': max(cycle_lens) if cycle_lens else 0,
    }


def build_sha_map(seed_base, trunc_bits):
    """Build complete SHA-256 map for truncated output."""
    mask = (1 << trunc_bits) - 1
    N = 1 << trunc_bits
    fwd = {}
    for val in range(N):
        W = list(seed_base)
        W[0] = val
        H = sha256_hash(W)
        fwd[val] = H[0] & mask
    return fwd, N


def main():
    np.random.seed(42)

    print("=" * 70)
    print("ЦИКЛОВАЯ АНОМАЛИЯ: 3× больше значений на циклах")
    print("=" * 70)

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("A. 20 разных seed_base, 16-bit truncation")
    print("=" * 70)

    sha_stats = []
    for trial in range(20):
        np.random.seed(trial * 137)
        seed = [np.random.randint(0, 2**32) for _ in range(16)]
        fwd, N = build_sha_map(seed, 16)
        stats = analyze_map(fwd, N)
        sha_stats.append(stats)

    print(f"  {'Trial':>5} {'#Cycles':>8} {'On cycles':>10} {'Max cycle':>10} {'Image%':>8}")
    for i, s in enumerate(sha_stats):
        print(f"  {i:>5} {s['n_cycles']:>8} {s['total_on_cycles']:>10} {s['max_cycle']:>10} {s['image_size']/N*100:>7.1f}%")

    sha_on_cycles = [s['total_on_cycles'] for s in sha_stats]
    sha_n_cycles = [s['n_cycles'] for s in sha_stats]

    print(f"\n  Summary (20 SHA-256 maps):")
    print(f"    On cycles: mean={np.mean(sha_on_cycles):.0f}, std={np.std(sha_on_cycles):.0f}")
    print(f"    #Cycles:   mean={np.mean(sha_n_cycles):.1f}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("B. 100 random functions for comparison")
    print("=" * 70)

    rand_stats = []
    N = 65536
    for trial in range(100):
        np.random.seed(10000 + trial)
        fwd = {v: np.random.randint(0, N) for v in range(N)}
        stats = analyze_map(fwd, N)
        rand_stats.append(stats)

    rand_on_cycles = [s['total_on_cycles'] for s in rand_stats]
    rand_n_cycles = [s['n_cycles'] for s in rand_stats]
    rand_max_cycle = [s['max_cycle'] for s in rand_stats]

    print(f"  Random functions (100 trials):")
    print(f"    On cycles: mean={np.mean(rand_on_cycles):.0f}, std={np.std(rand_on_cycles):.0f}")
    print(f"    #Cycles:   mean={np.mean(rand_n_cycles):.1f}")
    print(f"    Max cycle: mean={np.mean(rand_max_cycle):.0f}")

    print(f"\n  SHA-256 vs Random:")
    print(f"    On cycles: SHA={np.mean(sha_on_cycles):.0f} vs Random={np.mean(rand_on_cycles):.0f}")

    # Z-score
    z = (np.mean(sha_on_cycles) - np.mean(rand_on_cycles)) / np.std(rand_on_cycles)
    print(f"    Z-score: {z:.2f}")
    print(f"    → {'ANOMALOUS!' if abs(z) > 2 else 'Within normal range'}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("C. Distribution comparison")
    print("=" * 70)

    # Histogram of on-cycle counts
    print(f"\n  On-cycle count distribution:")
    for percentile in [5, 25, 50, 75, 95]:
        r_val = np.percentile(rand_on_cycles, percentile)
        print(f"    Random P{percentile}: {r_val:.0f}")

    sha_median = np.median(sha_on_cycles)
    # What percentile is SHA median in random distribution?
    sha_percentile = np.mean([1 for r in rand_on_cycles if r < sha_median]) * 100
    print(f"    SHA-256 median: {sha_median:.0f} (= P{sha_percentile:.0f} of random)")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("D. Vary truncation bits: 8, 10, 12, 14, 16")
    print("=" * 70)

    np.random.seed(42)
    seed_base = [np.random.randint(0, 2**32) for _ in range(16)]

    for bits in [8, 10, 12, 14, 16]:
        N_b = 1 << bits

        # SHA-256
        fwd_sha, _ = build_sha_map(seed_base, bits)
        sha_s = analyze_map(fwd_sha, N_b)

        # 10 random functions
        rand_oc = []
        for _ in range(10):
            fwd_r = {v: np.random.randint(0, N_b) for v in range(N_b)}
            rand_oc.append(analyze_map(fwd_r, N_b)['total_on_cycles'])

        expected = int(np.sqrt(np.pi * N_b / 8))

        print(f"  {bits}-bit (N={N_b:>6}): SHA on_cycles={sha_s['total_on_cycles']:>5}, "
              f"random={np.mean(rand_oc):>5.0f}±{np.std(rand_oc):>4.0f}, "
              f"theory={expected:>4}, "
              f"ratio={sha_s['total_on_cycles']/max(np.mean(rand_oc),1):.2f}x")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("E. WHY might SHA-256 have more on cycles?")
    print("=" * 70)

    # One possible explanation: our map is W[0] → H[0]&mask,
    # but W[1..15] are FIXED. This means the function is
    # f(x) = SHA256(x, seed[1..15])[0] & mask.
    #
    # This is NOT the same as iterating H→H→H.
    # Our original experiment ITERATED (H as next full input).
    # Here we only vary W[0]. Different setup!
    #
    # For true iteration: x is FULL hash, truncated for cycle detection.
    # For our map: x is a PARTIAL input.
    #
    # The difference matters: our map f: {0..N-1} → {0..N-1} is
    # a SPECIFIC random function (determined by seed[1..15]).
    # True iteration would be a DIFFERENT function.

    print(f"""
  NOTE: Our map is f(x) = SHA256(x || seed[1..15])[0] & mask.
  This is a FIXED random function, not iterated hashing.

  For a truly random function f: {{0..N-1}} → {{0..N-1}}:
    E[total on cycles] = sqrt(πN/8)
    Variance is HIGH: std ≈ 0.6 × mean

  The "3×" from earlier was ONE sample.
  With 20 samples: SHA-256 mean = {np.mean(sha_on_cycles):.0f}
  Random functions mean = {np.mean(rand_on_cycles):.0f}

  Z-score = {z:.2f}
""")

    # ═══════════════════
    print(f"{'=' * 70}")
    print("ВЕРДИКТ")
    print("=" * 70)


if __name__ == "__main__":
    main()
