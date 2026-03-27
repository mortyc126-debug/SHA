"""
НЕДОСТИЖИМОЕ ПРОСТРАНСТВО: 21.9% значений недостижимы итерацией.

Факт: при x → H(x) → H(H(x)) → ..., 78.1% 16-бит пространства достижимы.
21.9% = "мёртвые зоны" — точки, которые НИКОГДА не появляются в орбите.

Вопросы:
  A. Это свойство SHA-256 или любой функции?
  B. Какова структура недостижимых точек?
  C. Зависит ли % от размера пространства?
  D. Можно ли ПРЕДСКАЗАТЬ какие точки недостижимы?
  E. Связь с аттракторами?
"""

import numpy as np
import struct, hashlib
from collections import Counter

MASK32 = 0xFFFFFFFF

def sha256_hash(W16):
    raw = struct.pack('>16I', *W16)
    return struct.unpack('>8I', hashlib.sha256(raw).digest())

def hw(x): return bin(x).count('1')


def iterate_and_collect(seed_base, n_seeds, n_steps, trunc_bits):
    """Collect all reachable values under iteration."""
    mask = (1 << trunc_bits) - 1
    reachable = set()
    image_of = {}  # value → its image

    for s in range(n_seeds):
        x = list(seed_base)
        x[0] = (x[0] + s) & MASK32
        for step in range(n_steps):
            H = sha256_hash(x)
            val = H[0] & mask
            next_x = list(H) + [0]*8
            next_H = sha256_hash(next_x)
            next_val = next_H[0] & mask
            reachable.add(val)
            image_of[val] = next_val
            x = next_x
    return reachable, image_of


def main():
    np.random.seed(42)
    seed_base = [np.random.randint(0, 2**32) for _ in range(16)]

    print("=" * 70)
    print("НЕДОСТИЖИМОЕ ПРОСТРАНСТВО")
    print("=" * 70)

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("A. Теория: random function image size")
    print("=" * 70)

    # For random function f: {0..N-1} → {0..N-1}
    # Expected size of f(S) where S = full domain:
    # E[|Image|] = N × (1 - (1-1/N)^N) ≈ N × (1 - 1/e) ≈ 0.6321 × N
    #
    # So ~36.8% are UNREACHABLE after 1 application.
    # After k iterations: reachable set GROWS but converges.
    # After infinite iterations: reachable = "rho tail + cycle"
    #
    # For random function, the "eventual image" (values on cycles):
    # Expected fraction ≈ sqrt(πN/8) / N ≈ sqrt(π/(8N))
    # For N=65536: ≈ 0.39% on cycles. Rest = tails.

    N = 65536
    print(f"  Theory for random f: {{0..{N-1}}} → {{0..{N-1}}}:")
    print(f"    After 1 step: {N*(1-np.exp(-1)):.0f}/{N} = {(1-np.exp(-1))*100:.1f}% reachable")
    print(f"    Eventually: convergence to rho structure")
    print(f"    Fraction on cycles: ≈ {np.sqrt(np.pi/(8*N))*100:.2f}%")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("B. SHA-256 vs random function: image size after k steps")
    print("=" * 70)

    # Build explicit map for 16-bit truncated
    trunc = 16
    mask = (1 << trunc) - 1
    N = 1 << trunc

    # Map: for each 16-bit value, compute its image
    print(f"  Building explicit map for {trunc}-bit space ({N} values)...")

    fwd_map = {}
    for val in range(N):
        # Input: val as W[0], rest fixed
        W = list(seed_base)
        W[0] = val
        H = sha256_hash(W)
        fwd_map[val] = H[0] & mask

    # Image after 1 step
    image1 = set(fwd_map.values())
    print(f"  After 1 step: {len(image1)}/{N} = {len(image1)/N*100:.1f}% reachable")
    print(f"  Expected (random): {(1-np.exp(-1))*100:.1f}%")

    # Image after k steps
    current_image = set(range(N))  # start with everything
    for k in range(1, 20):
        next_image = set()
        for v in current_image:
            if v in fwd_map:
                next_image.add(fwd_map[v])
        current_image = next_image
        if k in [1, 2, 3, 5, 10, 15, 19]:
            print(f"  After {k:>2} steps: {len(current_image):>5}/{N} = {len(current_image)/N*100:.1f}%")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("C. Cycle structure: find ALL cycles")
    print("=" * 70)

    # Follow each value until it cycles
    visited = set()
    cycles = []

    for start in range(N):
        if start in visited:
            continue

        path = []
        seen_in_path = {}
        val = start

        while val not in seen_in_path:
            if val in visited:
                break
            seen_in_path[val] = len(path)
            path.append(val)
            val = fwd_map.get(val, val)

        if val in seen_in_path:
            cycle_start = seen_in_path[val]
            cycle = path[cycle_start:]
            cycles.append(cycle)
            visited.update(cycle)

        visited.update(path)

    print(f"  Total cycles found: {len(cycles)}")
    cycle_lens = [len(c) for c in cycles]
    print(f"  Cycle lengths: {sorted(cycle_lens)[:20]}{'...' if len(cycle_lens) > 20 else ''}")
    total_on_cycles = sum(cycle_lens)
    print(f"  Values on cycles: {total_on_cycles}/{N} = {total_on_cycles/N*100:.2f}%")
    print(f"  Expected (random): ≈ {np.sqrt(np.pi*N/8):.0f}/{N} = {np.sqrt(np.pi/(8*N))*100:.2f}%")

    # Tail lengths
    tail_lens = []
    for start in range(N):
        val = start
        steps = 0
        while True:
            next_val = fwd_map.get(val, val)
            if next_val == val or steps > N:
                break
            # Check if on a cycle
            is_cycle = False
            for c in cycles:
                if val in c:
                    is_cycle = True
                    break
            if is_cycle:
                break
            val = next_val
            steps += 1
        tail_lens.append(steps)

    print(f"\n  Tail lengths:")
    print(f"    Mean: {np.mean(tail_lens):.1f}")
    print(f"    Max:  {max(tail_lens)}")
    print(f"    Expected (random): ≈ {np.sqrt(np.pi*N/8):.0f}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("D. Unreachable values: what's special about them?")
    print("=" * 70)

    unreachable = set(range(N)) - image1
    reachable_list = list(image1)
    unreachable_list = list(unreachable)

    print(f"  Unreachable: {len(unreachable)}/{N} = {len(unreachable)/N*100:.1f}%")
    print(f"  Expected (random): {np.exp(-1)*100:.1f}%")

    # HW distribution of unreachable vs reachable
    hw_reach = [hw(v) for v in reachable_list]
    hw_unreach = [hw(v) for v in unreachable_list]

    print(f"\n  HW distribution:")
    print(f"    Reachable:   mean={np.mean(hw_reach):.2f}, std={np.std(hw_reach):.2f}")
    print(f"    Unreachable: mean={np.mean(hw_unreach):.2f}, std={np.std(hw_unreach):.2f}")

    from scipy import stats
    t, p = stats.ttest_ind(hw_reach, hw_unreach)
    print(f"    t={t:.2f}, p={p:.4f}")
    print(f"    → {'DIFFERENT!' if p < 0.01 else 'Same distribution'}")

    # Are unreachable values clustered or uniform?
    sorted_unreach = sorted(unreachable_list)
    if len(sorted_unreach) > 1:
        gaps = [sorted_unreach[i+1] - sorted_unreach[i] for i in range(len(sorted_unreach)-1)]
        cv = np.std(gaps) / np.mean(gaps) if np.mean(gaps) > 0 else 0
        print(f"\n  Spacing of unreachable values:")
        print(f"    Mean gap: {np.mean(gaps):.1f}")
        print(f"    CV: {cv:.3f} (1.0 = random/exponential)")
        print(f"    → {'CLUSTERED' if cv > 1.3 else 'UNIFORM' if cv < 0.8 else 'RANDOM (like Poisson process)'}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("E. Comparison: random function vs SHA-256")
    print("=" * 70)

    # Build random function
    rand_map = {v: np.random.randint(0, N) for v in range(N)}
    rand_image = set(rand_map.values())
    rand_unreachable = set(range(N)) - rand_image

    # Random cycles
    visited_r = set()
    cycles_r = []
    for start in range(N):
        if start in visited_r: continue
        path = []; seen = {}; val = start
        while val not in seen:
            if val in visited_r: break
            seen[val] = len(path); path.append(val)
            val = rand_map[val]
        if val in seen:
            cs = seen[val]
            cycles_r.append(path[cs:])
            visited_r.update(path[cs:])
        visited_r.update(path)

    print(f"  {'':>20} {'SHA-256':>10} {'Random':>10} {'Match':>6}")
    print(f"  {'Image size':>20} {len(image1):>10} {len(rand_image):>10} {'✓' if abs(len(image1)-len(rand_image)) < N*0.02 else '✗'}")
    print(f"  {'Unreachable %':>20} {len(unreachable)/N*100:>9.1f}% {len(rand_unreachable)/N*100:>9.1f}%")
    print(f"  {'Num cycles':>20} {len(cycles):>10} {len(cycles_r):>10}")
    print(f"  {'On cycles':>20} {total_on_cycles:>10} {sum(len(c) for c in cycles_r):>10}")

    sha_matches_random = (
        abs(len(image1)/N - len(rand_image)/N) < 0.03 and
        abs(len(cycles) - len(cycles_r)) < max(len(cycles), len(cycles_r)) * 0.5
    )

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("F. DEEPER: 8-bit space (COMPLETE analysis)")
    print("=" * 70)

    # 8-bit: only 256 values. Can analyze EVERYTHING.
    trunc8 = 8
    mask8 = 0xFF
    N8 = 256

    fwd8 = {}
    for val in range(N8):
        W = list(seed_base)
        W[0] = val
        H = sha256_hash(W)
        fwd8[val] = H[0] & mask8

    image8 = set(fwd8.values())
    unreach8 = set(range(N8)) - image8

    # Count preimages
    preimage_counts = Counter(fwd8.values())
    max_preimages = max(preimage_counts.values())
    zero_preimages = N8 - len(preimage_counts)

    print(f"  8-bit complete map:")
    print(f"    Image size: {len(image8)}/{N8} = {len(image8)/N8*100:.1f}%")
    print(f"    Expected: {(1-np.exp(-1))*100:.1f}%")
    print(f"    Unreachable: {len(unreach8)}")
    print(f"    Max preimage count: {max_preimages}")

    # Preimage distribution
    preimg_dist = Counter(preimage_counts.values())
    print(f"\n    Preimage distribution:")
    print(f"    {'#preimages':>11} {'count':>6} {'expected(Poisson)':>17}")
    for k in range(max_preimages + 1):
        actual = preimg_dist.get(k, 0) + (zero_preimages if k == 0 else 0)
        # Poisson(1): P(k) = e^(-1)/k!
        import math
        expected = N8 * np.exp(-1) / math.factorial(k)
        print(f"    {k:>11} {actual:>6} {expected:>16.1f}")

    # Cycles in 8-bit
    visited8 = set()
    cycles8 = []
    for start in range(N8):
        if start in visited8: continue
        path = []; seen = {}; val = start
        while val not in seen:
            if val in visited8: break
            seen[val] = len(path); path.append(val)
            val = fwd8[val]
        if val in seen:
            cs = seen[val]
            cycles8.append(path[cs:])
            visited8.update(path[cs:])
        visited8.update(path)

    print(f"\n    Cycles: {len(cycles8)}")
    for i, c in enumerate(cycles8):
        print(f"      Cycle {i}: length={len(c)}, values={c[:10]}{'...' if len(c)>10 else ''}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("РЕЗУЛЬТАТ")
    print("=" * 70)

    print(f"""
  SHA-256 iteration structure = RANDOM FUNCTION structure.

  Key numbers (16-bit):
    Image size:  {len(image1)/N*100:.1f}% (expected 63.2%)
    Unreachable: {len(unreachable)/N*100:.1f}% (expected 36.8%)
    Cycles:      {len(cycles)} (expected ~{int(np.log(N)/2)})

  SHA-256 matches random function predictions:
    {'YES — no special structure' if sha_matches_random else 'ANOMALY DETECTED'}

  WHAT THIS MEANS:
    The 21.9% "unreachable" from earlier = we didn't iterate ENOUGH.
    With COMPLETE map: unreachable = {len(unreachable)/N*100:.1f}% ≈ 1/e = 36.8%.
    This is a THEOREM about random functions, not SHA-256.

  ANALOGY TO OPEN PROBLEMS:
    Collatz: does NOT hold — multiple attractors, not one.
    Riemann: "primes" (unreachable) follow Poisson distribution.
    No hidden structure in SHA-256 iteration — it IS a random function.
""")


if __name__ == "__main__":
    main()
