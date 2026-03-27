"""
НОВОЕ НАПРАВЛЕНИЕ 1: Закон сохранения.

В физике: каждая симметрия → закон сохранения (теорема Нётер).
Энергия сохраняется из-за однородности времени.
Импульс — из-за однородности пространства.

SHA-256: есть ли ВЕЛИЧИНА, которая НЕ МЕНЯЕТСЯ через раунды?

Кандидаты:
  A. HW(state) — "масса" (количество единиц)
  B. XOR-fold: a⊕b⊕c⊕d⊕e⊕f⊕g⊕h — "чётность"
  C. ADD-fold: (a+b+c+d+e+f+g+h) mod 2^32 — "энергия"
  D. Что-то с участием W[r] — "полная энергия" (state + input)
  E. Неизвестная функция f(state, W) — ищем

Если найдём — это НАСТОЯЩИЙ закон физики нашего измерения.
"""

import numpy as np

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
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
]
IV = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]


def sha256_all_states(W16):
    W = list(W16)
    for r in range(16, 64):
        W.append((sigma1(W[r-2]) + W[r-7] + sigma0(W[r-15]) + W[r-16]) & MASK32)

    a, b, c, d, e, f, g, h = IV
    states = [(a, b, c, d, e, f, g, h)]
    for r in range(64):
        T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K[r]), W[r])
        T2 = add32(Sigma0(a), Maj(a,b,c))
        h,g,f,e = g,f,e,add32(d,T1)
        d,c,b,a = c,b,a,add32(T1,T2)
        states.append((a,b,c,d,e,f,g,h))

    return states, W


def main():
    np.random.seed(42)

    print("=" * 70)
    print("ЗАКОН СОХРАНЕНИЯ: ищем инвариант SHA-256")
    print("=" * 70)

    W16 = [np.random.randint(0, 2**32) for _ in range(16)]
    states, W = sha256_all_states(W16)

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("A. HW(state) — сохраняется ли 'масса'?")
    print("=" * 70)

    hws = [sum(hw(r) for r in s) for s in states]
    print(f"  HW(state) по раундам:")
    print(f"    r=0:  {hws[0]}")
    print(f"    r=32: {hws[32]}")
    print(f"    r=64: {hws[64]}")
    print(f"    Std:  {np.std(hws):.1f}")
    print(f"  → {'CONSERVED' if np.std(hws) < 2 else 'NOT conserved'} (varies by ±{np.std(hws):.0f})")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("B. XOR-fold: a⊕b⊕c⊕d⊕e⊕f⊕g⊕h")
    print("=" * 70)

    xor_folds = []
    for s in states:
        xf = 0
        for r in s: xf ^= r
        xor_folds.append(xf)

    # Check if constant
    all_same = all(x == xor_folds[0] for x in xor_folds)
    print(f"  XOR-fold constant? {all_same}")
    if not all_same:
        # How much does it change?
        diffs = [hw(xor_folds[i] ^ xor_folds[i+1]) for i in range(len(xor_folds)-1)]
        print(f"  HW(δ XOR-fold) per round: mean={np.mean(diffs):.1f}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("C. ADD-fold: (a+b+c+d+e+f+g+h) mod 2^32")
    print("=" * 70)

    add_folds = []
    for s in states:
        af = 0
        for r in s: af = (af + r) & MASK32
        add_folds.append(af)

    all_same_add = all(x == add_folds[0] for x in add_folds)
    print(f"  ADD-fold constant? {all_same_add}")
    if not all_same_add:
        diffs_add = [hw(add_folds[i] ^ add_folds[i+1]) for i in range(len(add_folds)-1)]
        print(f"  HW(δ ADD-fold) per round: mean={np.mean(diffs_add):.1f}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("D. TOTAL ENERGY: state + cumulative W + cumulative K")
    print("=" * 70)

    # Round function: new_state = f(old_state, W[r], K[r])
    # Maybe: state_sum + W_sum + K_sum = constant?
    # state_sum = a+b+c+d+e+f+g+h mod 2^32

    # Let's track: sum(state) - sum(W[0..r]) - sum(K[0..r])
    state_sums = []
    for i, s in enumerate(states):
        ss = 0
        for r in s: ss = (ss + r) & MASK32
        state_sums.append(ss)

    W_cum = [0]
    K_cum = [0]
    for r in range(64):
        W_cum.append((W_cum[-1] + W[r]) & MASK32)
        K_cum.append((K_cum[-1] + K[r]) & MASK32)

    # E = state_sum - W_cumulative - K_cumulative
    energies = []
    for i in range(65):
        E = (state_sums[i] - W_cum[i] - K_cum[i]) & MASK32
        energies.append(E)

    all_same_E = all(e == energies[0] for e in energies)
    print(f"  E = sum(state) - sum(W) - sum(K) constant? {all_same_E}")
    if not all_same_E:
        diffs_E = [hw(energies[i] ^ energies[i+1]) for i in range(64)]
        print(f"  HW(δE) per round: mean={np.mean(diffs_E):.1f}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("E. ALGEBRAIC INVARIANT: a+e vs b+f, c+g, d+h")
    print("=" * 70)

    # In SHA-256: a_new = T1+T2, e_new = d+T1
    # b_new = a, f_new = e
    # So: a_new + e_new = (T1+T2) + (d+T1) = 2*T1 + T2 + d
    # And: b_new + f_new = a + e (previous!)
    # So: (b+f)[r+1] = (a+e)[r] — EXACT CONSERVATION!

    print(f"  Hypothesis: (b+f)[r+1] = (a+e)[r] mod 2^32")
    print(f"  Because: b_new = a_old, f_new = e_old")

    conserved = True
    for r in range(64):
        ae_r = (states[r][0] + states[r][4]) & MASK32   # a+e at round r
        bf_r1 = (states[r+1][1] + states[r+1][5]) & MASK32  # b+f at round r+1
        if ae_r != bf_r1:
            conserved = False
            break

    print(f"  Verified: {conserved}")
    if conserved:
        print(f"  ★ (a+e)[r] = (b+f)[r+1] — EXACT for all 64 rounds!")

    # Extend: (b+f)[r] = (c+g)[r+1]?
    print(f"\n  Extended pipe conservation:")
    for pair1, pair2, name in [
        ((1,5), (2,6), "(b+f)[r] = (c+g)[r+1]"),
        ((2,6), (3,7), "(c+g)[r] = (d+h)[r+1]"),
        ((0,4), (1,5), "(a+e)[r] = (b+f)[r+1]"),
    ]:
        ok = True
        for r in range(64):
            v1 = (states[r][pair1[0]] + states[r][pair1[1]]) & MASK32
            v2 = (states[r+1][pair2[0]] + states[r+1][pair2[1]]) & MASK32
            if v1 != v2: ok = False; break
        print(f"    {name}: {ok}")

    # XOR version
    print(f"\n  XOR version:")
    for pair1, pair2, name in [
        ((0,4), (1,5), "(a⊕e)[r] = (b⊕f)[r+1]"),
        ((1,5), (2,6), "(b⊕f)[r] = (c⊕g)[r+1]"),
        ((2,6), (3,7), "(c⊕g)[r] = (d⊕h)[r+1]"),
    ]:
        ok = True
        for r in range(64):
            v1 = states[r][pair1[0]] ^ states[r][pair1[1]]
            v2 = states[r+1][pair2[0]] ^ states[r+1][pair2[1]]
            if v1 != v2: ok = False; break
        print(f"    {name}: {ok}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("F. DEEP SEARCH: 4-round conserved quantity")
    print("=" * 70)

    # After 4 rounds: (a+e)[0] → (b+f)[1] → (c+g)[2] → (d+h)[3]
    # So: (a+e)[r] = (d+h)[r+3] — CONSERVED OVER 4 ROUNDS!

    print(f"  Hypothesis: (a+e)[r] = (d+h)[r+3]")
    ok_4 = True
    for r in range(61):
        ae_r = (states[r][0] + states[r][4]) & MASK32
        dh_r3 = (states[r+3][3] + states[r+3][7]) & MASK32
        if ae_r != dh_r3: ok_4 = False; break
    print(f"  Verified: {ok_4}")
    if ok_4:
        print(f"  ★ (a+e)[r] = (d+h)[r+3] — CONSERVED ACROSS 4 ROUNDS!")

    # What about (a+e)[r] = something at [r+4]?
    # (d+h)[r+3] → at r+4: d_new = c_old, h_new = g_old
    # So we need (c+g)[r+3] = (d+h)[r+4]... which IS our conservation
    # Chain: (a+e)[r] = (b+f)[r+1] = (c+g)[r+2] = (d+h)[r+3]
    # At r+4: d_new=c, h_new=g → but c,g come from r+3's b,f... which is r+2's a,e
    # Wait: at r+4: the value enters the NEW a and e computation
    # It gets MODIFIED by T1,T2 → conservation BREAKS

    print(f"\n  Full chain: (a+e)[r] → (b+f)[r+1] → (c+g)[r+2] → (d+h)[r+3]")
    print(f"  At r+4: d enters e_new = d+T1 → value CONSUMED by T1")
    print(f"  Conservation lives for exactly 4 rounds (pipe length)")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("G. VERIFY ON MULTIPLE MESSAGES")
    print("=" * 70)

    for trial in range(10):
        W_t = [np.random.randint(0, 2**32) for _ in range(16)]
        st, _ = sha256_all_states(W_t)

        # Check (a+e)[r] = (d+h)[r+3] for all r
        ok = True
        for r in range(61):
            if (st[r][0]+st[r][4])&MASK32 != (st[r+3][3]+st[r+3][7])&MASK32:
                ok = False; break

        # Check XOR version
        ok_xor = True
        for r in range(61):
            if st[r][0]^st[r][4] != st[r+3][3]^st[r+3][7]:
                ok_xor = False; break

        print(f"  Message {trial}: ADD conserved={ok}, XOR conserved={ok_xor}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("ЗАКОН СОХРАНЕНИЯ — РЕЗУЛЬТАТ")
    print("=" * 70)


if __name__ == "__main__":
    main()
