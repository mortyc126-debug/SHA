"""
ТЕОРЕМА НЁТЕР ДЛЯ SHA-256.

Теорема Нётер: каждая непрерывная СИММЕТРИЯ → закон сохранения.
  Трансляционная симметрия по времени → сохранение энергии.
  Трансляционная симметрия в пространстве → сохранение импульса.
  Вращательная симметрия → сохранение момента.

У нас есть закон сохранения: (a+e)[r] = (b+f)[r+1] = ... = (d+h)[r+3].
Вопрос: КАКАЯ СИММЕТРИЯ его порождает?

Кандидаты:
  1. PIPE SHIFT SYMMETRY: b←a, c←b, d←c, f←e, g←f, h←g
     Это "трансляция по pipe" — сдвиг на 1 позицию.
  2. TWIN SYMMETRY: (a,b,c,d) ↔ (e,f,g,h) — зеркальная симметрия.
  3. ROUND TRANSLATION: R_r ≈ R_{r+1} (раунды "похожи").
"""

import numpy as np
import struct, hashlib

MASK32 = 0xFFFFFFFF

def hw(x): return bin(x).count('1')
def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def add32(x, y): return (x + y) & MASK32
def sub32(x, y): return (x - y) & MASK32
def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)

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


def round_fn(state, w, r):
    a,b,c,d,e,f,g,h = state
    T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K[r]), w)
    T2 = add32(Sigma0(a), Maj(a,b,c))
    return (add32(T1,T2), a, b, c, add32(d,T1), e, f, g)


def main():
    np.random.seed(42)

    print("=" * 70)
    print("ТЕОРЕМА НЁТЕР: какая симметрия → conservation law?")
    print("=" * 70)

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("1. PIPE SHIFT: симметрия сдвига")
    print("=" * 70)

    # Round function:
    #   a_new = T1 + T2    (computed)
    #   b_new = a           (shift)
    #   c_new = b           (shift)
    #   d_new = c           (shift)
    #   e_new = d + T1      (computed)
    #   f_new = e           (shift)
    #   g_new = f           (shift)
    #   h_new = g           (shift)

    # This IS a shift register with 2 injection points (a, e).
    # The shift operation S: (a,b,c,d,e,f,g,h) → (?,a,b,c,?,e,f,g)
    # is an EXACT SYMMETRY of the pipe structure.

    # PROOF that pipe shift → conservation:
    print(f"""
  PIPE SHIFT SYMMETRY:

  Define shift operator S:
    S(a,b,c,d,e,f,g,h) = (_, a, b, c, _, e, f, g)
    where _ = "new value" (computed from T1, T2)

  This is NOT a full symmetry of the round function
  (because a_new and e_new are COMPUTED, not shifted).
  But it IS a symmetry of the PIPE PART (6/8 registers).

  Conservation law follows DIRECTLY from the pipe shift:
    b[r+1] = a[r]   (by definition)
    f[r+1] = e[r]   (by definition)
    → (b+f)[r+1] = a[r] + e[r] = (a+e)[r]    QED

  This is NOT Noether's theorem (continuous symmetry → conserved quantity).
  This is SIMPLER: it's a DISCRETE symmetry (shift by 1 position).
  The conservation is a DIRECT CONSEQUENCE of the shift, not a deep result.
""")

    # ═══════════════════
    print(f"{'=' * 70}")
    print("2. DEEPER QUESTION: is there a HIDDEN symmetry?")
    print("=" * 70)

    # The pipe conservation is trivial. But are there NON-TRIVIAL conserved
    # quantities that come from HIDDEN symmetries?

    # Search: for functions f(state) that are conserved across rounds.
    # f(state[r]) = f(state[r+1]) for all r.

    # Test candidates:
    candidates = {
        'a+e': lambda s: add32(s[0], s[4]),
        'a⊕e': lambda s: s[0] ^ s[4],
        'a×e mod 2^32': lambda s: (s[0] * s[4]) & MASK32,
        'Σ(all)': lambda s: add32(add32(add32(s[0],s[1]),add32(s[2],s[3])),
                                   add32(add32(s[4],s[5]),add32(s[6],s[7]))),
        '⊕(all)': lambda s: s[0]^s[1]^s[2]^s[3]^s[4]^s[5]^s[6]^s[7],
        'a+b+c+d': lambda s: add32(add32(s[0],s[1]),add32(s[2],s[3])),
        'e+f+g+h': lambda s: add32(add32(s[4],s[5]),add32(s[6],s[7])),
        'a⊕b⊕c⊕d': lambda s: s[0]^s[1]^s[2]^s[3],
        'HW(state) mod 8': lambda s: sum(hw(x) for x in s) % 8,
        'a-e': lambda s: sub32(s[0], s[4]),
        '(a+e)⊕(b+f)': lambda s: add32(s[0],s[4]) ^ add32(s[1],s[5]),
    }

    W = [np.random.randint(0, 2**32) for _ in range(16)]
    # Expand W
    W_full = list(W)
    from functools import reduce
    for r in range(16, 64):
        s0 = rotr(W_full[r-15],7) ^ rotr(W_full[r-15],18) ^ (W_full[r-15]>>3)
        s1 = rotr(W_full[r-2],17) ^ rotr(W_full[r-2],19) ^ (W_full[r-2]>>10)
        W_full.append(add32(add32(add32(s1, W_full[r-7]), s0), W_full[r-16]))

    state = tuple(IV)
    states = [state]
    for r in range(64):
        state = round_fn(state, W_full[r], r)
        states.append(state)

    print(f"  Testing candidate conserved quantities:")
    print(f"  {'Candidate':<20} {'Conserved?':>10} {'Max δ':>8}")

    for name, fn in candidates.items():
        values = [fn(s) for s in states]
        # Check: is value constant?
        if all(isinstance(v, int) for v in values):
            diffs = [hw(values[i] ^ values[i+1]) for i in range(64)]
            max_diff = max(diffs)
            is_conserved = max_diff == 0
        else:
            diffs = [abs(values[i] - values[i+1]) for i in range(64)]
            max_diff = max(diffs)
            is_conserved = max_diff == 0
        print(f"  {name:<20} {'YES ✓' if is_conserved else 'NO':>10} {max_diff:>8}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("3. SEARCH: conserved MODULAR quantities")
    print("=" * 70)

    # f(state) mod p conserved for small p?
    # Maybe some function is conserved mod 2, or mod 4?

    print(f"  Testing f(state) mod p for small p:")
    for p in [2, 3, 4, 5, 7, 8, 16, 256]:
        for name, fn in list(candidates.items())[:6]:
            values = [fn(s) % p if isinstance(fn(s), int) else int(fn(s)) % p for s in states]
            diffs = [abs(values[i] - values[i+1]) for i in range(64)]
            is_conserved = all(d == 0 for d in diffs)
            if is_conserved:
                print(f"    {name} mod {p}: CONSERVED ★")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("4. SYMMETRY GROUP of the round function")
    print("=" * 70)

    # What transformations T of state commute with the round function?
    # T(R(s)) = R(T(s)) for all s → T is a SYMMETRY.

    # Test: simple transformations
    transformations = {
        'identity': lambda s: s,
        'complement all': lambda s: tuple(x ^ MASK32 for x in s),
        'swap a↔e': lambda s: (s[4],s[1],s[2],s[3],s[0],s[5],s[6],s[7]),
        'rotate left': lambda s: (s[1],s[2],s[3],s[0],s[5],s[6],s[7],s[4]),
        'swap halves': lambda s: (s[4],s[5],s[6],s[7],s[0],s[1],s[2],s[3]),
        'negate a': lambda s: (sub32(0, s[0]),s[1],s[2],s[3],s[4],s[5],s[6],s[7]),
        'reverse': lambda s: tuple(reversed(s)),
    }

    print(f"  {'Transform':<20} {'Commutes?':>10} {'Mean δ':>8}")
    for name, T in transformations.items():
        deltas = []
        for _ in range(500):
            s = tuple(np.random.randint(0, 2**32) for _ in range(8))
            w = np.random.randint(0, 2**32)
            r = np.random.randint(0, 64)

            # T(R(s)) vs R(T(s))
            TR = T(round_fn(s, w, r))
            RT = round_fn(T(s), w, r)
            d = sum(hw(TR[i]^RT[i]) for i in range(8))
            deltas.append(d)

        mean_d = np.mean(deltas)
        is_symmetry = mean_d < 0.01
        print(f"  {name:<20} {'YES ✓' if is_symmetry else 'NO':>10} {mean_d:>7.1f}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("5. APPROXIMATE SYMMETRIES")
    print("=" * 70)

    # Maybe exact symmetry doesn't exist, but APPROXIMATE symmetry does.
    # T such that |T(R(s)) - R(T(s))| is SMALL (not zero).

    # The smallest δ from above: which transform is CLOSEST to a symmetry?
    best_transforms = []
    for name, T in transformations.items():
        if name == 'identity':
            continue
        deltas = []
        for _ in range(1000):
            s = tuple(np.random.randint(0, 2**32) for _ in range(8))
            w = np.random.randint(0, 2**32)
            r = np.random.randint(0, 64)
            TR = T(round_fn(s, w, r))
            RT = round_fn(T(s), w, r)
            deltas.append(sum(hw(TR[i]^RT[i]) for i in range(8)))
        best_transforms.append((name, np.mean(deltas)))

    best_transforms.sort(key=lambda x: x[1])
    print(f"  Approximate symmetries (sorted by δ):")
    for name, d in best_transforms:
        quality = d / 128  # 0 = perfect, 1 = random
        print(f"    {name:<20}: δ={d:>6.1f}/256 ({quality:.3f})")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("6. NOETHER'S THEOREM FOR SHA-256")
    print("=" * 70)

    print(f"""
  ═══════════════════════════════════════════════════════════════

  РЕЗУЛЬТАТ:

  ТОЧНЫХ СИММЕТРИЙ: только identity.
  SHA-256 round function НЕ коммутирует НИ С КАКОЙ нетривиальной
  трансформацией state. Группа симметрий = {{identity}}.

  ПРИБЛИЖЁННЫХ СИММЕТРИЙ: нет (все δ ≈ 64-128 = random).

  ЗАКОН СОХРАНЕНИЯ (a+e)[r] = (b+f)[r+1]:
  Это НЕ следствие симметрии в смысле Нётер.
  Это ТОЖДЕСТВО pipe structure:
    b_new ≡ a_old  (по определению)
    f_new ≡ e_old  (по определению)
    → (b+f)_new ≡ (a+e)_old   QED

  ТЕОРЕМА НЁТЕР НЕ ПРИМЕНИМА:
    1. Нет непрерывной группы Ли (дискретная система).
    2. Нет нетривиальных симметрий (кроме identity).
    3. Conservation law = pipe tautology, не symmetry consequence.

  НОВЫЙ ЗАКОН (вместо Нётер):

  "PIPE IDENTITY THEOREM":
    Для любой системы вида:
      x_new = f(state, input)
      y_new = x_old           (shift/pipe)
    выполняется: y[r+1] = x[r] (тождественно).

    Для SHA-256: 6 таких тождеств (b=a, c=b, d=c, f=e, g=f, h=g).
    Комбинируя: (a+e)[r] = (d+h)[r+3] (4-step conservation).

    Это НЕ глубокая симметрия.
    Это СТРУКТУРНОЕ СЛЕДСТВИЕ shift-register design.
    Работает для ЛЮБОГО shift-register hash.

  ПОЧЕМУ ЭТО ВАЖНО:
    В физике: симметрия → conservation → predictions.
    В SHA-256: pipe identity → conservation → constraints on δ.
    Practical value: conservation constrains valid differential paths
    (shown: 3-6× enrichment in conservation_attack.py).

  ═══════════════════════════════════════════════════════════════
""")


if __name__ == "__main__":
    main()
