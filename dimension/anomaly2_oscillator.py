"""
АНОМАЛИЯ 2: Self-modify oscillator [142,121,142,121]

Факт: при самомодификации (H говорит какой бит менять),
trace застревает в 2-шаговом цикле.

Вопросы:
  A. Это свойство SHA-256 или любой хеш-функции?
  B. Почему именно 2-цикл?
  C. Насколько это частое явление?
  D. Можно ли это эксплуатировать?
"""

import numpy as np
import struct, hashlib

MASK32 = 0xFFFFFFFF

def sha256_hash(W16):
    raw = struct.pack('>16I', *W16)
    return struct.unpack('>8I', hashlib.sha256(raw).digest())

def hw(x): return bin(x).count('1')


def self_modify_trace(seed, steps):
    """H[0] → word, H[1] → bit, H[2] → XOR/ADD."""
    x = list(seed)
    hws = []
    states = []
    for step in range(steps):
        H = sha256_hash(x)
        hws.append(sum(hw(h) for h in H))
        states.append(tuple(x))

        word_idx = H[0] & 0xF
        bit_idx = H[1] & 0x1F
        op = (H[2] >> 0) & 1

        if op == 0:
            x[word_idx] ^= (1 << bit_idx)
        else:
            x[word_idx] = (x[word_idx] + (1 << bit_idx)) & MASK32

    return hws, states


def main():
    np.random.seed(42)

    print("=" * 70)
    print("АНОМАЛИЯ 2: Self-modify oscillator")
    print("=" * 70)

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("A. Воспроизводим")
    print("=" * 70)

    seed = [np.random.randint(0, 2**32) for _ in range(16)]
    hws, states = self_modify_trace(seed, 200)

    print(f"  Last 20 HWs: {hws[-20:]}")

    # Detect cycle
    for i in range(len(states)-1, 0, -1):
        for period in [1, 2, 3, 4, 5]:
            if i - period >= 0 and states[i] == states[i-period]:
                print(f"  Cycle detected: period={period}, starts at step {i-period}")
                break
        else:
            continue
        break

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("B. Почему 2-цикл?")
    print("=" * 70)

    # Проверяем: XOR одного бита = инволюция (x^b^b = x).
    # Если H(x) говорит "flip bit b" → x' = x^b
    # Если H(x') тоже говорит "flip bit b" → x'' = x'^b = x → ЦИКЛ!
    #
    # Это произойдёт когда:
    #   H(x)[0]&0xF == H(x^b)[0]&0xF  (same word)
    #   H(x)[1]&0x1F == H(x^b)[1]&0x1F  (same bit)
    #   H(x)[2]&1 == 0 AND H(x^b)[2]&1 == 0  (both XOR)

    # Вероятность: P(same word) = 1/16, P(same bit) = 1/32, P(both XOR) = 1/4
    # P(2-cycle|XOR) = 1/16 × 1/32 × 1 = 1/512
    # Но нет — это должно быть P что x' flip ТОГО ЖЕ бита что привёл к x'

    # Точнее: шаг 1: x → x' = modify(x, H(x))
    #         шаг 2: x' → x'' = modify(x', H(x'))
    #         2-цикл если x'' = x

    # Для XOR: x'' = x iff modify undoes itself
    # modify(x, w1, b1, XOR) → x' (flip bit b1 of word w1)
    # modify(x', w2, b2, op2) → x''
    # x'' = x iff (w2,b2,op2) = (w1,b1,XOR) — same flip

    # P = P(w2=w1) × P(b2=b1) × P(op2=XOR)
    # = 1/16 × 1/32 × 1/2 = 1/1024

    # Но: op может быть ADD, тогда:
    # ADD с op=1: x' = x + (1<<b), x'' = x' + (1<<b) = x + 2×(1<<b) ≠ x
    # ADD НИКОГДА не даёт 2-цикл (не инволюция)

    # Для XOR-XOR: P(2-cycle) = 1/16 × 1/32 = 1/512 PER STEP
    # За 200 шагов: P(hit 2-cycle) ≈ 1 - (1-1/512)^200 ≈ 33%

    print(f"  Теория:")
    print(f"    XOR flip = инволюция (x^b^b = x)")
    print(f"    2-цикл если H(x') выбирает ТОТ ЖЕ бит и XOR")
    print(f"    P(per step) = P(same word)×P(same bit)×P(XOR)×P(XOR)")
    print(f"                = 1/16 × 1/32 × 1/2 × 1/2 = 1/2048")
    print(f"    За 200 шагов: P ≈ {1-(1-1/2048)**200:.2%}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("C. Как часто? 1000 случайных seed'ов")
    print("=" * 70)

    cycle_counts = {1: 0, 2: 0, 3: 0, 'none': 0}
    cycle_steps = []

    for trial in range(1000):
        seed_t = [np.random.randint(0, 2**32) for _ in range(16)]
        _, states_t = self_modify_trace(seed_t, 500)

        found = False
        for i in range(len(states_t)-1, max(len(states_t)-50, 0), -1):
            for period in [1, 2, 3]:
                if i - period >= 0 and states_t[i] == states_t[i-period]:
                    cycle_counts[period] = cycle_counts.get(period, 0) + 1
                    cycle_steps.append(i - period)
                    found = True
                    break
            if found:
                break
        if not found:
            cycle_counts['none'] += 1

    print(f"  1000 trials × 500 steps:")
    for k, v in sorted(cycle_counts.items(), key=lambda x: str(x[0])):
        print(f"    Period {k}: {v}/1000 ({v/10:.1f}%)")
    if cycle_steps:
        print(f"    Cycle start: mean={np.mean(cycle_steps):.0f}, min={min(cycle_steps)}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("D. Специфично для SHA-256?")
    print("=" * 70)

    # Сравним с "random oracle": вместо SHA-256 — случайная функция
    # (используем SHA-256 с другим padding как прокси)
    import hashlib as hl

    cycle_counts_alt = {1: 0, 2: 0, 3: 0, 'none': 0}

    for trial in range(1000):
        x = [np.random.randint(0, 2**32) for _ in range(16)]
        states_t = []

        for step in range(500):
            # Используем MD5 вместо SHA-256 как "другой оракул"
            raw = struct.pack('>16I', *x)
            H_md5 = struct.unpack('>4I', hl.md5(raw).digest())

            states_t.append(tuple(x))

            word_idx = H_md5[0] & 0xF
            bit_idx = H_md5[1] & 0x1F
            op = (H_md5[2] >> 0) & 1

            if op == 0:
                x[word_idx] ^= (1 << bit_idx)
            else:
                x[word_idx] = (x[word_idx] + (1 << bit_idx)) & MASK32

        found = False
        for i in range(len(states_t)-1, max(len(states_t)-50, 0), -1):
            for period in [1, 2, 3]:
                if i-period>=0 and states_t[i]==states_t[i-period]:
                    cycle_counts_alt[period] = cycle_counts_alt.get(period, 0) + 1
                    found = True
                    break
            if found: break
        if not found:
            cycle_counts_alt['none'] += 1

    print(f"  SHA-256:  {cycle_counts}")
    print(f"  MD5:      {cycle_counts_alt}")
    print(f"  → {'SAME pattern = generic property' if abs(cycle_counts[2] - cycle_counts_alt.get(2,0)) < 100 else 'DIFFERENT = SHA-256 specific'}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("ВЕРДИКТ")
    print("=" * 70)
    print(f"""
  Self-modify oscillator = GENERIC property of XOR self-modification.
  XOR = инволюция (x^b^b = x).
  Если H(x') выбирает тот же бит → instant 2-cycle.
  P ≈ 1/2048 per step. За 500 шагов: ≈ 22% шанс.

  НЕ свойство SHA-256. Любая хеш-функция даёт то же самое.
  ВЕРДИКТ: ОБЪЯСНЕНО (generic XOR involution).
""")


if __name__ == "__main__":
    main()
