"""
PIPE/NODE RATIO: нужны ли трубы для безопасности?

SHA-256: 6 труб + 2 узла = 75%/25% per round.
Что если:
  A. 0 труб + 8 узлов (все NODE) — "full-mix"
  B. 4 трубы + 4 узла (50/50)
  C. 6 труб + 2 узла (standard SHA-256)
  D. 7 труб + 1 узел (minimal node)

В нашем измерении:
  Трубы = бесплатная передача δ (κ=1)
  Узлы = нелинейное преобразование (κ<1)

  Больше узлов = больше нелинейности per round?
  Или: трубы НУЖНЫ для propagation?
"""

import numpy as np
import struct

MASK32 = 0xFFFFFFFF
K_const = [0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
           0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
           0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
           0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174]
IV = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]

def hw(x): return bin(x).count('1')
def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def add32(x, y): return (x + y) & MASK32
def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)

def round_standard(state, W_r, r):
    """Standard: 6 pipes + 2 nodes."""
    a,b,c,d,e,f,g,h = state
    T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K_const[r%16]), W_r)
    T2 = add32(Sigma0(a), Maj(a,b,c))
    return (add32(T1,T2), a, b, c, add32(d,T1), e, f, g)

def round_all_node(state, W_r, r):
    """All NODE: each register recomputed independently."""
    a,b,c,d,e,f,g,h = state
    T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K_const[r%16]), W_r)
    T2 = add32(Sigma0(a), Maj(a,b,c))
    # Instead of pipes: ALL registers get NODE computation
    a_new = add32(T1, T2)
    b_new = add32(a, Sigma0(b))  # node instead of pipe
    c_new = add32(b, Maj(b,c,d))
    d_new = add32(c, Ch(c,d,e))
    e_new = add32(d, T1)
    f_new = add32(e, Sigma1(f))  # node instead of pipe
    g_new = add32(f, Ch(f,g,h))
    h_new = add32(g, Maj(g,h,a))
    return (a_new, b_new, c_new, d_new, e_new, f_new, g_new, h_new)

def round_half_node(state, W_r, r):
    """Half NODE (4 pipes + 4 nodes)."""
    a,b,c,d,e,f,g,h = state
    T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K_const[r%16]), W_r)
    T2 = add32(Sigma0(a), Maj(a,b,c))
    a_new = add32(T1, T2)     # NODE
    b_new = a                  # PIPE
    c_new = add32(b, Maj(b,c,d))  # NODE
    d_new = c                  # PIPE
    e_new = add32(d, T1)      # NODE
    f_new = e                  # PIPE
    g_new = add32(f, Ch(f,g,h))   # NODE
    h_new = g                  # PIPE
    return (a_new, b_new, c_new, d_new, e_new, f_new, g_new, h_new)

def round_min_node(state, W_r, r):
    """Minimal NODE (7 pipes + 1 node)."""
    a,b,c,d,e,f,g,h = state
    T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K_const[r%16]), W_r)
    T2 = add32(Sigma0(a), Maj(a,b,c))
    a_new = add32(T1, T2)  # only NODE
    return (a_new, a, b, c, d, e, f, g)  # 7 pipes


def hash_variant(W16, n_rounds, round_fn):
    a,b,c,d,e,f,g,h = IV
    state = (a,b,c,d,e,f,g,h)
    for r in range(n_rounds):
        state = round_fn(state, W16[r%16], r)
    return tuple(add32(IV[i], state[i]) for i in range(8))


def compute_metrics(round_fn, n_rounds):
    np.random.seed(42)
    W_base = [np.random.randint(0, 2**32) for _ in range(16)]
    H_base = hash_variant(W_base, n_rounds, round_fn)

    # Sensitivity
    sens = []
    for _ in range(500):
        w,b = np.random.randint(0,16), np.random.randint(0,32)
        W2 = list(W_base); W2[w]^=(1<<b)
        H2 = hash_variant(W2, n_rounds, round_fn)
        sens.append(sum(hw(H_base[i]^H2[i]) for i in range(8)))

    # Curvature
    Ks = []
    for _ in range(300):
        b1,b2 = np.random.choice(512, 2, replace=False)
        W1=list(W_base); W1[b1//32]^=(1<<(b1%32))
        W2=list(W_base); W2[b2//32]^=(1<<(b2%32))
        W12=list(W_base); W12[b1//32]^=(1<<(b1%32)); W12[b2//32]^=(1<<(b2%32))
        H1=hash_variant(W1,n_rounds,round_fn)
        H2=hash_variant(W2,n_rounds,round_fn)
        H12=hash_variant(W12,n_rounds,round_fn)
        nl = sum(hw((H_base[i]^H12[i])^((H_base[i]^H1[i])^(H_base[i]^H2[i]))) for i in range(8))
        Ks.append(nl)

    # rank(T)
    T_rows = []
    for word in range(16):
        for bit in range(32):
            W_mod = list(W_base); W_mod[word]^=(1<<bit)
            H_mod = hash_variant(W_mod, n_rounds, round_fn)
            row = []
            for w in range(8):
                d = H_base[w]^H_mod[w]
                for b in range(32): row.append((d>>b)&1)
            T_rows.append(row)
    T = np.array(T_rows, dtype=np.uint8)
    rT = np.linalg.matrix_rank(T.astype(float))

    return np.mean(sens), np.mean(Ks), rT


def main():
    np.random.seed(42)

    print("=" * 70)
    print("PIPE/NODE RATIO: роль труб в безопасности")
    print("=" * 70)

    variants = [
        ("Standard (6P+2N)", round_standard),
        ("All-node (0P+8N)", round_all_node),
        ("Half (4P+4N)", round_half_node),
        ("Min-node (7P+1N)", round_min_node),
    ]

    # =================================================================
    print(f"\n{'=' * 70}")
    print("1. МЕТРИКИ ПРИ r=16 (все слова активны)")
    print("=" * 70)

    print(f"\n  {'Variant':<22} {'Sens':>7} {'K(curv)':>8} {'rank(T)':>8}")
    for name, fn in variants:
        s, k, rT = compute_metrics(fn, 16)
        marker = " ★" if rT == 256 and k > 100 else ""
        print(f"  {name:<22} {s:6.1f} {k:7.1f} {rT:>7}{marker}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("2. СКОРОСТЬ ДИФФУЗИИ: rounds to saturation")
    print("=" * 70)

    for name, fn in variants:
        print(f"\n  {name}:")
        for n_r in [4, 8, 12, 16, 24]:
            s, k, rT = compute_metrics(fn, n_r)
            bar = "█" * int(s / 8)
            print(f"    r={n_r:2d}: Sens={s:5.1f} K={k:5.1f} rank={rT:3d}  {bar}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("3. ЗНАЧЕНИЕ PIPE/NODE RATIO")
    print("=" * 70)

    print(f"""
  В нашем измерении:

  ТРУБЫ (pipes):
    - Бесплатная передача δ (κ=1)
    - Создают SHIFT REGISTER cascade (a→b→c→d, e→f→g→h)
    - Обеспечивают REBOOT (a-repair → convergence за 8 раундов)
    - НО: не добавляют нелинейность

  УЗЛЫ (nodes):
    - Нелинейное преобразование (κ<1)
    - Carry + Ch/Maj = source of rank(CE)
    - НО: без труб δ не распространяется по регистрам

  ОПТИМУМ: SHA-256 (6P+2N = 75% труб) —
    2 узла per round ДОСТАТОЧНО для rank=256 за 16 раундов.
    6 труб обеспечивают ПОЛНУЮ propagation за 4 раунда.

  All-node: быстрее saturation? Или медленнее?
  Min-node: достаточно? Или нет?
""")


if __name__ == "__main__":
    main()
