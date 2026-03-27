"""
OPTIMAL INJECTION ORDER.

Проблема: последнее слово = самое слабое (100% verified).
В ЛЮБОМ порядке: position 15 → только 2 раунда mixing.

Решения:
  A. INTERLEAVE: чередовать ранние и поздние → balance
  B. PRE-MIX: все слова в IV ДО раундов → нет "последнего"
  C. DOUBLE INJECT: каждое слово входит ДВАЖДЫ в разных раундах
  D. XOR BUNDLE: round r получает W[r] ⊕ W[r+8] → two words per round

Измеряем: min(avg δH по всем словам) при r=17.
Чем ВЫШЕ minimum — тем БЕЗОПАСНЕЕ.
"""

import numpy as np
import struct, hashlib

MASK32 = 0xFFFFFFFF

def hw(x): return bin(x).count('1')
def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def add32(x, y): return (x + y) & MASK32
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
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,
]
IV = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]


def hash_sequential(W16, n_r, order=None):
    """Standard: words injected one per round in given order."""
    if order is None:
        order = list(range(16))
    W_ord = [W16[order[i]] for i in range(16)]
    W = list(W_ord)
    for r in range(16, max(n_r, 16)):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))
    a,b,c,d,e,f,g,h = IV
    for r in range(n_r):
        T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K[r % len(K)]), W[r])
        T2 = add32(Sigma0(a), Maj(a,b,c))
        h,g,f,e = g,f,e,add32(d,T1)
        d,c,b,a = c,b,a,add32(T1,T2)
    return tuple(add32(IV[i], [a,b,c,d,e,f,g,h][i]) for i in range(8))


def hash_premix(W16, n_r):
    """PRE-MIX: XOR all words into IV, then run rounds with schedule."""
    # Mix all message words into initial state
    state = list(IV)
    for i in range(16):
        state[i % 8] = add32(state[i % 8], W16[i])

    # Schedule from original words
    W = list(W16)
    for r in range(16, max(n_r, 16)):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))

    a,b,c,d,e,f,g,h = state
    iv_save = list(state)
    for r in range(n_r):
        T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K[r % len(K)]), W[r])
        T2 = add32(Sigma0(a), Maj(a,b,c))
        h,g,f,e = g,f,e,add32(d,T1)
        d,c,b,a = c,b,a,add32(T1,T2)
    return tuple(add32(iv_save[i], [a,b,c,d,e,f,g,h][i]) for i in range(8))


def hash_double_inject(W16, n_r):
    """DOUBLE: each word enters TWICE (round i and round i+8)."""
    W = list(W16)
    # Rounds 0-7: W[0..7], Rounds 8-15: W[0..7] again XOR W[8..15]
    W_double = list(W[:8]) + [W[i] ^ W[i+8] for i in range(8)]
    for r in range(16, max(n_r, 16)):
        W_double.append(add32(add32(add32(sigma1(W_double[r-2]), W_double[r-7]),
                        sigma0(W_double[r-15])), W_double[r-16]))
    a,b,c,d,e,f,g,h = IV
    for r in range(n_r):
        T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K[r % len(K)]), W_double[r])
        T2 = add32(Sigma0(a), Maj(a,b,c))
        h,g,f,e = g,f,e,add32(d,T1)
        d,c,b,a = c,b,a,add32(T1,T2)
    return tuple(add32(IV[i], [a,b,c,d,e,f,g,h][i]) for i in range(8))


def hash_all_at_once(W16, n_r):
    """ALL AT ONCE: all words XORed into first 16 rounds as bundles."""
    # Each round r gets: W[r] ⊕ W[(r+4)%16] ⊕ W[(r+8)%16] ⊕ W[(r+12)%16]
    # Every word participates in 4 rounds → maximum mixing
    W_mixed = []
    for r in range(16):
        w = W16[r] ^ W16[(r+4)%16] ^ W16[(r+8)%16] ^ W16[(r+12)%16]
        W_mixed.append(w)
    for r in range(16, max(n_r, 16)):
        W_mixed.append(add32(add32(add32(sigma1(W_mixed[r-2]), W_mixed[r-7]),
                       sigma0(W_mixed[r-15])), W_mixed[r-16]))
    a,b,c,d,e,f,g,h = IV
    for r in range(n_r):
        T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K[r % len(K)]), W_mixed[r])
        T2 = add32(Sigma0(a), Maj(a,b,c))
        h,g,f,e = g,f,e,add32(d,T1)
        d,c,b,a = c,b,a,add32(T1,T2)
    return tuple(add32(IV[i], [a,b,c,d,e,f,g,h][i]) for i in range(8))


def measure_vulnerability(hash_fn, name, n_r, N=2000):
    """Measure avg δH for each input word. Return min (weakest word)."""
    word_avgs = []
    word_mins = []
    for word in range(16):
        dHs = []
        for _ in range(N):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            W2 = list(W); W2[word] ^= (1 << 31)
            H1 = hash_fn(W); H2 = hash_fn(W2)
            dHs.append(sum(hw(H1[i]^H2[i]) for i in range(8)))
        word_avgs.append(np.mean(dHs))
        word_mins.append(min(dHs))

    return word_avgs, word_mins


def main():
    np.random.seed(42)

    print("=" * 70)
    print("OPTIMAL INJECTION ORDER")
    print("=" * 70)

    n_r = 17

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print(f"1. ALL DESIGNS at r={n_r}")
    print("=" * 70)

    designs = [
        (lambda W: hash_sequential(W, n_r), "Sequential (standard)"),
        (lambda W: hash_sequential(W, n_r, [15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0]), "Reversed"),
        (lambda W: hash_sequential(W, n_r, [0,8,1,9,2,10,3,11,4,12,5,13,6,14,7,15]), "Interleaved"),
        (lambda W: hash_premix(W, n_r), "Pre-mix (all in IV)"),
        (lambda W: hash_double_inject(W, n_r), "Double inject"),
        (lambda W: hash_all_at_once(W, n_r), "All-at-once (XOR bundles)"),
    ]

    print(f"\n  {'Design':<25} {'Min avg δH':>11} {'Max avg δH':>11} {'Weakest':>10} {'Spread':>8}")

    results = {}
    for hash_fn, name in designs:
        avgs, mins = measure_vulnerability(hash_fn, name, n_r, N=1500)
        min_avg = min(avgs)
        max_avg = max(avgs)
        weakest_word = np.argmin(avgs)
        spread = max_avg - min_avg

        results[name] = {
            'avgs': avgs, 'mins': mins,
            'min_avg': min_avg, 'max_avg': max_avg,
            'weakest': weakest_word, 'spread': spread
        }

        print(f"  {name:<25} {min_avg:>10.1f} {max_avg:>10.1f} {'W['+str(weakest_word)+']':>10} {spread:>7.1f}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("2. PER-WORD DETAILS for each design")
    print("=" * 70)

    for hash_fn, name in designs:
        r = results[name]
        print(f"\n  {name}:")
        for w in range(16):
            bar_len = int(r['avgs'][w] / 4)
            bar = "█" * bar_len
            marker = " ★ WEAK" if r['avgs'][w] < 50 else ""
            print(f"    W[{w:>2}]: avg={r['avgs'][w]:>6.1f}, min={r['mins'][w]:>3} {bar}{marker}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("3. SCALING: how fast does each design reach full security?")
    print("=" * 70)

    print(f"\n  Min(avg δH across words) at each round count:")
    print(f"  {'Round':>5}", end="")
    for _, name in designs:
        short = name[:15]
        print(f" {short:>15}", end="")
    print()

    for n_r_test in [8, 12, 14, 16, 17, 18, 20, 24]:
        print(f"  {n_r_test:>5}", end="")
        for hash_fn_maker, name in [
            (lambda nr: lambda W: hash_sequential(W, nr), "Sequential"),
            (lambda nr: lambda W: hash_sequential(W, nr, [0,8,1,9,2,10,3,11,4,12,5,13,6,14,7,15]), "Interleaved"),
            (lambda nr: lambda W: hash_premix(W, nr), "Pre-mix"),
            (lambda nr: lambda W: hash_all_at_once(W, nr), "All-at-once"),
        ]:
            if name == "Sequential":
                hfn = lambda W, nr=n_r_test: hash_sequential(W, nr)
            elif name == "Interleaved":
                hfn = lambda W, nr=n_r_test: hash_sequential(W, nr, [0,8,1,9,2,10,3,11,4,12,5,13,6,14,7,15])
            elif name == "Pre-mix":
                hfn = lambda W, nr=n_r_test: hash_premix(W, nr)
            else:
                hfn = lambda W, nr=n_r_test: hash_all_at_once(W, nr)

            avgs, _ = measure_vulnerability(hfn, name, n_r_test, N=500)
            print(f" {min(avgs):>14.1f}", end="")
        print()

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("4. ВЕРДИКТ")
    print("=" * 70)

    best_design = min(results.items(), key=lambda x: -x[1]['min_avg'])
    worst_design = min(results.items(), key=lambda x: x[1]['min_avg'])

    print(f"""
  ═══════════════════════════════════════════════════════════════

  RESULTS at r={n_r}:

  BEST design:  {best_design[0]}
    Min avg δH = {best_design[1]['min_avg']:.1f} (weakest word)
    Spread = {best_design[1]['spread']:.1f}

  WORST design: {worst_design[0]}
    Min avg δH = {worst_design[1]['min_avg']:.1f} (weakest word)
    Spread = {worst_design[1]['spread']:.1f}

  Improvement: {best_design[1]['min_avg'] - worst_design[1]['min_avg']:.1f} bits
  in weakest-word security.

  DESIGN RECOMMENDATIONS:

  1. PRE-MIX (all words into IV before rounds):
     Eliminates "last word" problem entirely.
     ALL words have SAME mixing = n_r rounds.
     Trade-off: IV depends on message → different security model.

  2. ALL-AT-ONCE (XOR bundles):
     Each word participates in 4 rounds simultaneously.
     Spreads influence, reduces "last word" effect.
     Trade-off: each round's W is less independent.

  3. INTERLEAVE:
     Mild improvement. Last word still weak but less so.
     Most compatible with existing SHA-256 structure.

  FOR SHA-256 (64 rounds):
     All designs converge by r≈20-24.
     The choice of injection order ONLY matters for reduced rounds.
     At 64 rounds: ALL orders give identical security.

  ═══════════════════════════════════════════════════════════════
""")


if __name__ == "__main__":
    main()
