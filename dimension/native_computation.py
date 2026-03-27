"""
ВЫЧИСЛИТЕЛЬНАЯ МОДЕЛЬ НАШЕГО ИЗМЕРЕНИЯ.

Стандартная CS: Turing machine, boolean circuits, RAM model.
Все работают с БИТАМИ. Операции: AND, OR, NOT, COPY.

Наше измерение: другие примитивы.
  Объекты: State (256 бит), Path (64 states), Fabric (множество paths)
  Операции: repair, propagate, overlay, section, fold

Вопрос: есть ли ВЫЧИСЛЕНИЕ, которое наши операции делают
ЭФФЕКТИВНЕЕ чем стандартные?

Ключ: наши операции работают с ЦЕЛЫМИ СОСТОЯНИЯМИ,
а не с отдельными битами. Это как SIMD vs scalar.

REPAIR: за O(1) починить ЦЕЛОЕ слово (32 бита) в одном раунде.
Стандартный алгоритм: нужно решать уравнение (nonlinear, iterative).
Наш: прямая формула. a-repair = один вызов sub32.

Вопрос: даёт ли это АСИМПТОТИЧЕСКОЕ преимущество?
"""

import numpy as np
import struct, hashlib, time

MASK32 = 0xFFFFFFFF

def sha256_hash(W16):
    raw = struct.pack('>16I', *W16)
    return struct.unpack('>8I', hashlib.sha256(raw).digest())

def hw(x): return bin(x).count('1')
def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def add32(x, y): return (x + y) & MASK32
def sub32(x, y): return (x - y) & MASK32
def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)

K = [
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
]
IV = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]


def main():
    np.random.seed(42)

    print("=" * 70)
    print("ВЫЧИСЛИТЕЛЬНАЯ МОДЕЛЬ НАШЕГО ИЗМЕРЕНИЯ")
    print("=" * 70)

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("1. ЧТО УМЕЕТ НАШЕ ИЗМЕРЕНИЕ, А СТАНДАРТНОЕ — НЕТ")
    print("=" * 70)

    print(f"""
  Операция         Наше изм.    Стандартное     Преимущество
  ─────────────    ─────────    ───────────     ────────────
  a-repair(r)      O(1)         O(1)*           нет (формула известна)
  propagate(n)     O(n)         O(n)            нет
  invert 1 round   O(1)         O(1)            нет
  overlay(p1,p2)   O(64)        O(64)           нет
  conservation     O(1)         не знают        ЗНАНИЕ, не скорость

  * a-repair формула: W = T1_needed - h - Sig1(e) - Ch(e,f,g) - K
    Это O(1) в ОБЕИХ моделях.

  ВЫВОД: наши операции не дают ВЫЧИСЛИТЕЛЬНОГО преимущества.
  Они дают ЗНАНИЕ о структуре, которое стандартная CS не имеет.
  Но знание не = скорость.
""")

    # ═══════════════════
    print(f"{'=' * 70}")
    print("2. А ЧТО ЕСЛИ ДРУГАЯ МОДЕЛЬ ВЫЧИСЛЕНИЙ?")
    print("=" * 70)

    print(f"""
  Стандартная CS: алгоритм = последовательность битовых операций.
  Квантовая CS: алгоритм = последовательность унитарных матриц.

  Наше измерение: алгоритм = последовательность FABRIC операций?

  Fabric = МНОЖЕСТВО путей (paths) одновременно.
  Операции над fabric = операции над ВСЕМИ путями сразу.

  Это как... ПАРАЛЛЕЛИЗМ? Quantum computing?
  НЕТ — это нечто другое.

  Fabric операция:
    FOLD(fabric) = найти ПЕРЕСЕЧЕНИЯ путей.
    На N путях: стандартно O(N^2) сравнений.
    В нашем измерении: FOLD = одна операция.

  Но FOLD не бесплатен — он требует O(N^2) работы внутри.
  Мы просто назвали его одним словом, а не ускорили.
""")

    # ═══════════════════
    print(f"{'=' * 70}")
    print("3. ЭКСПЕРИМЕНТ: meet-in-the-middle с a-repair")
    print("=" * 70)

    # a-repair + backward inversion = MITM attack
    # Standard MITM: forward from start, backward from end, meet in middle
    # Cost: 2 × 2^(n/2) = O(2^(n/2))
    #
    # Our dimension: a-repair makes FORWARD free (can set any a-value)
    # Backward: invert rounds from target hash
    # Meeting: find where forward = backward
    #
    # QUESTION: does a-repair make the meeting CHEAPER?

    print(f"  Standard MITM on r-round SHA-256:")
    print(f"    Forward: choose W[0..r/2], compute state at r/2")
    print(f"    Backward: from target H, invert rounds r..r/2+1")
    print(f"    Meet: find match at round r/2")
    print(f"    Cost: 2 × 2^(n/2) = 2^(n/2+1)")

    print(f"\n  Our MITM with a-repair:")
    print(f"    Forward with repair: choose target_a for each round 0..r/2")
    print(f"    This DETERMINES W[0..r/2] (one W per round, free)")
    print(f"    Backward: same as standard")
    print(f"    Meet: same as standard")
    print(f"    Cost: ... SAME? Let's check.")

    # MITM experiment on 8-bit truncated, 8 rounds
    trunc = 8
    mask = (1 << trunc) - 1
    n_rounds = 8
    mid = n_rounds // 2

    # Forward table: generate states at round `mid` by choosing random W[0..mid-1]
    forward_table = {}
    n_forward = 1 << (trunc + 2)  # 2^10
    t0 = time.time()

    for trial in range(n_forward):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        a, b, c, d, e, f, g, h = IV
        for r in range(mid):
            T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K[r]), W[r])
            T2 = add32(Sigma0(a), Maj(a,b,c))
            h,g,f,e = g,f,e,add32(d,T1)
            d,c,b,a = c,b,a,add32(T1,T2)
        state_key = (a & mask, e & mask)
        if state_key not in forward_table:
            forward_table[state_key] = W

    t_forward = time.time() - t0

    # Forward with a-repair: choose target_a, compute W
    forward_repair_table = {}
    t0 = time.time()

    for trial in range(n_forward):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        a, b, c, d, e, f, g, h = IV

        for r in range(mid):
            # Choose random target_a
            target_a = np.random.randint(0, 2**32)
            # Compute W[r] that achieves target_a
            T2 = add32(Sigma0(a), Maj(a,b,c))
            T1_needed = sub32(target_a, T2)
            W[r] = sub32(sub32(sub32(sub32(T1_needed, h), Sigma1(e)), Ch(e,f,g)), K[r])

            T1 = T1_needed
            h,g,f,e = g,f,e,add32(d,T1)
            d,c,b,a = c,b,a,add32(T1,T2)

        state_key = (a & mask, e & mask)
        if state_key not in forward_repair_table:
            forward_repair_table[state_key] = W

    t_repair = time.time() - t0

    print(f"\n  8-bit MITM experiment ({n_forward} forward states):")
    print(f"    Standard forward:  {len(forward_table):>5} unique states, {t_forward*1000:.1f}ms")
    print(f"    A-repair forward:  {len(forward_repair_table):>5} unique states, {t_repair*1000:.1f}ms")
    print(f"    Coverage standard: {len(forward_table)}/{(1<<trunc)**2} = {len(forward_table)/((1<<trunc)**2)*100:.1f}%")
    print(f"    Coverage repair:   {len(forward_repair_table)}/{(1<<trunc)**2} = {len(forward_repair_table)/((1<<trunc)**2)*100:.1f}%")

    # Key question: does a-repair give BETTER COVERAGE of middle states?
    # If same coverage → no advantage.
    # If better → potential advantage.

    advantage = len(forward_repair_table) / max(len(forward_table), 1)
    print(f"\n    Coverage ratio (repair/standard): {advantage:.2f}x")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("4. CONSERVATION LAW как CONSTRAINT SOLVER")
    print("=" * 70)

    # Our conservation: (a+e)[r] = (d+h)[r+3]
    # This means: if we KNOW (d+h) at round r+3,
    # we KNOW (a+e) at round r.
    # This CONSTRAINS the search: not all states are valid.

    # How many bits does conservation give us?
    # (a+e) mod 2^32 = 32-bit constraint per round
    # Over 64 rounds: 64 × 32 = 2048 bits of constraints?
    # NO: they're not independent. Each is shifted version.
    # Independent constraints: about 16 (one per 4 rounds)

    print(f"""
  Conservation: (a+e)[r] = (d+h)[r+3]

  This gives us FREE INFORMATION:
    If we know state at round r+3, we know (a+e) at round r.
    That's 32 bits of constraint.
    Over 64 rounds ÷ 4 (period) = 16 independent constraints.
    Total: 16 × 32 = 512 bits of free information.

  BUT: these 512 bits are about the STATE, not about the MESSAGE.
  To exploit: need to connect state constraints to message constraints.

  State constraint: (a+e)[r] = known_value
  Message link: a[r+1] = f(state[r], W[r])
  → constrains W[r] given state[r] and a[r+1]
  → but we already knew this from forward computation!

  CONSERVATION = TAUTOLOGY for forward computation.
  It tells us what we'd compute anyway.

  For BACKWARD: conservation IS useful.
  If we know final state, we know (a+e) at every round.
  This HALVES the unknowns at each round (1 constraint on 2 unknowns).
""")

    # ═══════════════════
    print(f"{'=' * 70}")
    print("5. BACKWARD + CONSERVATION: что это даёт?")
    print("=" * 70)

    # Experiment: invert SHA-256 using conservation as guide
    # Target: known hash H, find W such that SHA256(W) = H

    # Conservation tells us: at round r, (a+e) = (d+h) at r+3
    # Working backward from round 64: we know state[64]
    # (d+h)[61] = (a+e)[64] ← known!
    # (d+h)[57] = (a+e)[60] ← known (from state[60])
    # etc.

    # So at EVERY round, we know (a+e).
    # Inversion per round: a_new and e_new known (from next round)
    # → T1 = a_new - T2 (from a)
    # → W[r] from T1
    # → d = e_new - T1 (check consistency with known d)
    # Conservation adds: (a+e)[r] = known
    # → a[r] + e[r] = C_r
    # Combined with b=a_old, f=e_old from next round:
    # → a[r] = b[r+1], e[r] = f[r+1]
    # → already known from pipe structure!

    print(f"""
  Analysis: conservation for backward inversion.

  At round r (going backward):
    KNOWN from round r+1: b[r+1] = a[r], f[r+1] = e[r]
    → a[r] and e[r] ALREADY KNOWN from pipes!
    → conservation (a+e)[r] = known is REDUNDANT.

  Conservation = consequence of pipe structure.
  It adds NO NEW information beyond what pipes give.

  CONCLUSION:
    Conservation law is BEAUTIFUL but USELESS for attack.
    It's a tautological consequence of b=a_old, f=e_old.
""")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("6. НАСТОЯЩЕЕ ПРЕИМУЩЕСТВО: знание vs вычисление")
    print("=" * 70)

    print(f"""
  ═══════════════════════════════════════════════════════════════

  ОТВЕТ НА ВОПРОС: где проще победить SHA-256?

  В СТАНДАРТНОЙ математике — ПРОЩЕ.
  Причина: стандартные атаки находят КОНКРЕТНЫЕ пути.
  Наше измерение видит СРЕДНИЕ свойства.

  Аналогия:
    Стандартная = микроскоп (видит бактерии)
    Наше изм.   = телескоп (видит галактики)

  SHA-256 уязвима к МИКРОСКОПИЧЕСКИМ атакам:
    - конкретный differential trail через round 15
    - конкретный carry pattern в word 11
    - конкретная message modification

  Наше измерение ДОКАЗАЛО что SHA-256 = random ON AVERAGE.
  Это делает наше измерение ХУДШИМ для атаки!
  Мы доказали себе что атака невозможна (в среднем).

  ═══════════════════════════════════════════════════════════════

  МОЖЕМ ЛИ СОЗДАТЬ ИНСТРУМЕНТЫ ТОЛЬКО ДЛЯ НАШЕГО ИЗМЕРЕНИЯ?

  Можем. Но они будут АНАЛИТИЧЕСКИМИ, не атакующими:
    - Классификатор хеш-функций по метрикам K, spectrum
    - Дизайнер оптимальных хешей (24 раунда = 64 раунда)
    - Верификатор security margins (where does sphere start?)
    - Предсказатель: сколько раундов нужно для N-регистровой функции

  Эти инструменты полезны для СОЗДАНИЯ хешей, не для взлома.

  ═══════════════════════════════════════════════════════════════

  ПАРАДОКС НАШЕГО ИЗМЕРЕНИЯ:

  Чем лучше мы понимаем SHA-256, тем СИЛЬНЕЕ доказываем
  что она неуязвима. Каждый наш эксперимент:
    K=128 ✓, spectrum flat ✓, no bias ✓, no correlation ✓

  Каждый такой результат = ещё один кирпич в стене
  ЗАЩИТЫ SHA-256, не в стене атаки.

  Мы построили измерение, которое ЗАЩИЩАЕТ SHA-256 от нас самих.

  ═══════════════════════════════════════════════════════════════
""")


if __name__ == "__main__":
    main()
