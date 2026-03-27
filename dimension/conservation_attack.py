"""
CONSERVATION-GUIDED DIFFERENTIAL SEARCH.

Стандартный подход: ищем differential trail перебором.
Наш подход: conservation law ОГРАНИЧИВАЕТ пространство поиска.

Закон: (a+e)[r] = (b+f)[r+1] = (c+g)[r+2] = (d+h)[r+3]
Значит: если знаем δ(a+e) на раунде r, мы ТОЧНО знаем
δ(b+f) на r+1, δ(c+g) на r+2, δ(d+h) на r+3.

Это СОКРАЩАЕТ пространство поиска дифференциальных путей:
  Стандартно: 8 независимых δ-регистров × 32 бита = 256 бит/раунд
  С conservation: фактически 2 независимых (δa, δe) + 6 зависимых
  = 64 бита/раунд.

Reduction: 256 → 64 = в 4 раза меньше свободы = в 2^192 раз
меньше пространство поиска? Нет — не совсем, потому что зависимость
нелинейная. Но направление верное.

ПЛАН:
1. Используем conservation чтобы ПРЕДСКАЗАТЬ δ на раундах r+1..r+3
2. Ищем δW только для двух свободных регистров (a, e)
3. Проверяем: conservation-guided search быстрее random?
"""

import numpy as np
import struct, hashlib
import time

MASK32 = 0xFFFFFFFF

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


def sha256_states(W16, n_rounds):
    """Return all intermediate states."""
    W = list(W16)
    for r in range(16, max(n_rounds, 16)):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))
    states = []
    a,b,c,d,e,f,g,h = IV
    states.append((a,b,c,d,e,f,g,h))
    for r in range(n_rounds):
        T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K[r]), W[r])
        T2 = add32(Sigma0(a), Maj(a,b,c))
        h,g,f,e = g,f,e,add32(d,T1)
        d,c,b,a = c,b,a,add32(T1,T2)
        states.append((a,b,c,d,e,f,g,h))
    final = tuple(add32(IV[i], states[n_rounds][i]) for i in range(8))
    return states, W, final


def conservation_predict(states1, states2, r):
    """Use conservation law to predict δ at r+1, r+2, r+3 from r."""
    # (a+e)[r] = (b+f)[r+1] = (c+g)[r+2] = (d+h)[r+3]
    ae1 = add32(states1[r][0], states1[r][4])
    ae2 = add32(states2[r][0], states2[r][4])
    delta_ae = ae1 ^ ae2

    # At r+1: δ(b+f) should equal δ(a+e)[r]
    bf1 = add32(states1[r+1][1], states1[r+1][5])
    bf2 = add32(states2[r+1][1], states2[r+1][5])
    delta_bf = bf1 ^ bf2

    return delta_ae, delta_bf, hw(delta_ae), hw(delta_bf)


def main():
    np.random.seed(42)

    print("=" * 70)
    print("CONSERVATION-GUIDED DIFFERENTIAL SEARCH")
    print("=" * 70)

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("1. Conservation law верификация: δ(a+e) предсказывает δ(b+f)")
    print("=" * 70)

    # Show that conservation law gives EXACT prediction
    W1 = [np.random.randint(0, 2**32) for _ in range(16)]
    W2 = list(W1); W2[0] ^= 1

    states1, _, _ = sha256_states(W1, 24)
    states2, _, _ = sha256_states(W2, 24)

    print(f"  Verification: δ(a+e)[r] == δ(b+f)[r+1]?")
    for r in range(20):
        dae, dbf, hw_ae, hw_bf = conservation_predict(states1, states2, r)
        match = "✓" if dae == dbf else "✗"
        print(f"    r={r:>2}: HW(δ(a+e))={hw_ae:>3}, HW(δ(b+f)[r+1])={hw_bf:>3} {match}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("2. Conservation CONSTRAINS differential paths")
    print("=" * 70)

    # In a differential path (δa[r], δb[r], ..., δh[r]):
    # Conservation says: δ(a+e)[r] flows through pipes for 4 rounds.
    # This means: δb[r+1] + δf[r+1] is DETERMINED by δa[r] + δe[r].
    #
    # Standard search: choose all 8 δ-registers independently.
    # Conservation-aware: only choose δa[r] and δe[r], rest follows.
    #
    # But: b[r+1] = a[r] (exact copy), f[r+1] = e[r] (exact copy).
    # So δb[r+1] = δa[r] and δf[r+1] = δe[r].
    # This means δ(b+f)[r+1] = δ(a+e)[r] TRIVIALLY (it's the same sum).
    #
    # The NON-trivial part: at r+4, when d enters the T1 computation,
    # the conserved quantity CONSTRAINS what T1 can produce.

    print(f"""
  Conservation chain:
    (a+e)[r] → (b+f)[r+1] → (c+g)[r+2] → (d+h)[r+3]

  At round r+3: d and h carry the conserved quantity.
  At round r+4: d feeds into e_new = d + T1.
    → The conserved value ENTERS the nonlinear computation.

  Constraint: e_new[r+4] = d[r+4] + T1[r+4]
    where d[r+4] = c[r+3] = b[r+2] = a[r+1]
    and a[r+1] was computed with freedom (W[r+1]).

  So: e_new[r+4] depends on a[r+1] + T1[r+4].
    a[r+1] carries part of (a+e)[r].
    This LINKS round r+4 back to round r.

  The link is: e_new[r+4] partially determined by (a+e)[r].
""")

    # ═══════════════════
    print(f"{'=' * 70}")
    print("3. PRACTICAL USE: conservation-pruned search")
    print("=" * 70)

    # Strategy: for each candidate δW, compute (a+e) at each round.
    # Use conservation to predict: will this δW lead to SMALL δH?
    #
    # Key insight: if (a+e) is SMALL at some round r,
    # then (d+h) will be SMALL at r+3.
    # If (d+h) is small: d contributes less to e_new → less amplification.
    #
    # Search: find δW where (a+e) is small at CRITICAL rounds (r=16-19).

    n_rounds = 20
    n_trials = 10000

    # Standard search: random δW, measure δH
    t0 = time.time()
    best_random = 256
    random_hws = []
    for _ in range(n_trials):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W1)
        # Random 1-bit flip in random word
        w = np.random.randint(0, 16)
        b = np.random.randint(0, 32)
        W2[w] ^= (1 << b)

        st1, _, H1 = sha256_states(W1, n_rounds)
        st2, _, H2 = sha256_states(W2, n_rounds)
        dH = sum(hw(H1[i]^H2[i]) for i in range(8))
        random_hws.append(dH)
        if dH < best_random: best_random = dH

    t_random = time.time() - t0

    # Conservation-guided search: same trials, but FILTER based on
    # (a+e) magnitude at critical rounds
    t0 = time.time()
    best_guided = 256
    guided_hws = []
    filtered_count = 0

    for _ in range(n_trials):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W1)
        w = np.random.randint(0, 16)
        b = np.random.randint(0, 32)
        W2[w] ^= (1 << b)

        # EARLY CHECK: compute only first 17 rounds
        st1, _, _ = sha256_states(W1, 17)
        st2, _, _ = sha256_states(W2, 17)

        # Conservation filter: check (a+e) at round 16
        ae1_16 = add32(st1[16][0], st1[16][4])
        ae2_16 = add32(st2[16][0], st2[16][4])
        delta_ae_16 = hw(ae1_16 ^ ae2_16)

        # If δ(a+e) at r=16 is LOW → promising (conservation will keep it low for 3 more rounds)
        if delta_ae_16 < 12:  # threshold: low conserved difference
            filtered_count += 1
            # Compute full hash
            _, _, H1 = sha256_states(W1, n_rounds)
            _, _, H2 = sha256_states(W2, n_rounds)
            dH = sum(hw(H1[i]^H2[i]) for i in range(8))
            guided_hws.append(dH)
            if dH < best_guided: best_guided = dH

    t_guided = time.time() - t0

    print(f"  r={n_rounds}, {n_trials} trials:")
    print(f"  Random search:         best={best_random}, mean={np.mean(random_hws):.1f}, time={t_random:.1f}s")
    g_mean = f"{np.mean(guided_hws):.1f}" if guided_hws else "N/A"
    g_best = str(best_guided) if guided_hws else "N/A"
    print(f"  Conservation-guided:   best={g_best}, mean={g_mean}, "
          f"filtered={filtered_count}/{n_trials}, time={t_guided:.1f}s")

    if guided_hws:
        # Compare distributions
        print(f"\n  Distribution comparison:")
        for t in [100, 110, 115, 120]:
            r_pct = sum(1 for h in random_hws if h <= t) / len(random_hws) * 100
            g_pct = sum(1 for h in guided_hws if h <= t) / len(guided_hws) * 100 if guided_hws else 0
            print(f"    HW ≤ {t}: random={r_pct:.1f}%, guided={g_pct:.1f}%")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("4. MULTI-ROUND conservation filter")
    print("=" * 70)

    # Filter on MULTIPLE conservation quantities simultaneously
    best_multi = 256
    multi_hws = []
    multi_filtered = 0

    for _ in range(n_trials):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W1)
        w = np.random.randint(0, 16)
        b = np.random.randint(0, 32)
        W2[w] ^= (1 << b)

        st1, _, _ = sha256_states(W1, 18)
        st2, _, _ = sha256_states(W2, 18)

        # Check conservation quantities at rounds 14, 15, 16, 17
        total_conserved_hw = 0
        for r in [14, 15, 16, 17]:
            ae1 = add32(st1[r][0], st1[r][4])
            ae2 = add32(st2[r][0], st2[r][4])
            total_conserved_hw += hw(ae1 ^ ae2)

        # Filter: total conserved HW < threshold
        if total_conserved_hw < 40:
            multi_filtered += 1
            _, _, H1 = sha256_states(W1, n_rounds)
            _, _, H2 = sha256_states(W2, n_rounds)
            dH = sum(hw(H1[i]^H2[i]) for i in range(8))
            multi_hws.append(dH)
            if dH < best_multi: best_multi = dH

    print(f"  Multi-round filter (sum δ(a+e)[14..17] < 40):")
    print(f"    Filtered: {multi_filtered}/{n_trials}")
    print(f"    Best δH: {best_multi if multi_hws else 'N/A'}")
    if multi_hws:
        print(f"    Mean δH of filtered: {np.mean(multi_hws):.1f}")
        print(f"    vs random mean: {np.mean(random_hws):.1f}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("5. COMBINED: a-repair + conservation filter + schedule-aware")
    print("=" * 70)

    # The full pipeline:
    # 1. Choose δW[0] (1 bit)
    # 2. a-repair rounds 1-15
    # 3. Conservation filter on rounds 14-17
    # 4. If passes → compute full hash

    best_combined = 256
    combined_hws = []
    combined_filtered = 0
    total_combined = 20000

    t0 = time.time()
    for trial in range(total_combined):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W1)
        W2[0] ^= (1 << np.random.randint(0, 32))

        # a-repair
        a1,b1,c1,d1,e1,f1,g1,h1 = IV
        a2,b2,c2,d2,e2,f2,g2,h2 = IV

        T1=add32(add32(add32(add32(h1,Sigma1(e1)),Ch(e1,f1,g1)),K[0]),W1[0])
        T2=add32(Sigma0(a1),Maj(a1,b1,c1))
        h1,g1,f1,e1=g1,f1,e1,add32(d1,T1); d1,c1,b1,a1=c1,b1,a1,add32(T1,T2)

        T1=add32(add32(add32(add32(h2,Sigma1(e2)),Ch(e2,f2,g2)),K[0]),W2[0])
        T2=add32(Sigma0(a2),Maj(a2,b2,c2))
        h2,g2,f2,e2=g2,f2,e2,add32(d2,T1); d2,c2,b2,a2=c2,b2,a2,add32(T1,T2)

        for r in range(1, 16):
            T1_1=add32(add32(add32(add32(h1,Sigma1(e1)),Ch(e1,f1,g1)),K[r]),W1[r])
            T2_1=add32(Sigma0(a1),Maj(a1,b1,c1))
            a1_new=add32(T1_1,T2_1); e1_new=add32(d1,T1_1)
            T2_2=add32(Sigma0(a2),Maj(a2,b2,c2))
            T1_n=sub32(a1_new,T2_2)
            W2[r]=sub32(sub32(sub32(sub32(T1_n,h2),Sigma1(e2)),Ch(e2,f2,g2)),K[r])
            h1,g1,f1,e1=g1,f1,e1,e1_new; d1,c1,b1,a1=c1,b1,a1,a1_new
            h2,g2,f2,e2=g2,f2,e2,add32(d2,T1_n); d2,c2,b2,a2=c2,b2,a2,a1_new

        # Conservation filter: compute through round 18
        st1, _, _ = sha256_states(W1, 18)
        st2, _, _ = sha256_states(W2, 18)

        # Sum of conserved δ at rounds 15-17
        cons_hw = 0
        for r in [15, 16, 17]:
            ae1 = add32(st1[r][0], st1[r][4])
            ae2 = add32(st2[r][0], st2[r][4])
            cons_hw += hw(ae1 ^ ae2)

        if cons_hw < 30:  # pass filter
            combined_filtered += 1
            _, _, H1 = sha256_states(W1, n_rounds)
            _, _, H2 = sha256_states(W2, n_rounds)
            dH = sum(hw(H1[i]^H2[i]) for i in range(8))
            combined_hws.append(dH)
            if dH < best_combined: best_combined = dH

    t_combined = time.time() - t0

    print(f"  a-repair + conservation + search ({total_combined} trials, {t_combined:.1f}s):")
    print(f"    Pass filter: {combined_filtered}/{total_combined} ({combined_filtered/total_combined*100:.1f}%)")
    if combined_hws:
        print(f"    Best δH: {best_combined}")
        print(f"    Mean δH (filtered): {np.mean(combined_hws):.1f}")
    print(f"    Random best (same count): {best_random}")
    print(f"    Advantage: {best_random - best_combined if combined_hws else 0} bits")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("6. SCALING across rounds")
    print("=" * 70)

    for n_r in [18, 20, 22, 24]:
        best_c = 256
        best_r = 256
        for _ in range(5000):
            W1 = [np.random.randint(0, 2**32) for _ in range(16)]
            W2 = list(W1); W2[0] ^= (1 << np.random.randint(0, 32))

            # Random
            _, _, H1r = sha256_states(W1, n_r)
            _, _, H2r = sha256_states(W2, n_r)
            dH_r = sum(hw(H1r[i]^H2r[i]) for i in range(8))
            if dH_r < best_r: best_r = dH_r

            # a-repair + conservation
            a1,b1,c1,d1,e1,f1,g1,h1 = IV
            a2,b2,c2,d2,e2,f2,g2,h2 = IV
            T1=add32(add32(add32(add32(h1,Sigma1(e1)),Ch(e1,f1,g1)),K[0]),W1[0])
            T2=add32(Sigma0(a1),Maj(a1,b1,c1))
            h1,g1,f1,e1=g1,f1,e1,add32(d1,T1); d1,c1,b1,a1=c1,b1,a1,add32(T1,T2)
            T1=add32(add32(add32(add32(h2,Sigma1(e2)),Ch(e2,f2,g2)),K[0]),W2[0])
            T2=add32(Sigma0(a2),Maj(a2,b2,c2))
            h2,g2,f2,e2=g2,f2,e2,add32(d2,T1); d2,c2,b2,a2=c2,b2,a2,add32(T1,T2)
            W2c = list(W2)
            for r in range(1, 16):
                T1_1=add32(add32(add32(add32(h1,Sigma1(e1)),Ch(e1,f1,g1)),K[r]),W1[r])
                T2_1=add32(Sigma0(a1),Maj(a1,b1,c1))
                a1_new=add32(T1_1,T2_1); e1_new=add32(d1,T1_1)
                T2_2=add32(Sigma0(a2),Maj(a2,b2,c2))
                T1_n=sub32(a1_new,T2_2)
                W2c[r]=sub32(sub32(sub32(sub32(T1_n,h2),Sigma1(e2)),Ch(e2,f2,g2)),K[r])
                h1,g1,f1,e1=g1,f1,e1,e1_new; d1,c1,b1,a1=c1,b1,a1,a1_new
                h2,g2,f2,e2=g2,f2,e2,add32(d2,T1_n); d2,c2,b2,a2=c2,b2,a2,a1_new

            _, _, H1c = sha256_states(W1, n_r)
            _, _, H2c = sha256_states(W2c, n_r)
            dH_c = sum(hw(H1c[i]^H2c[i]) for i in range(8))
            if dH_c < best_c: best_c = dH_c

        adv = best_r - best_c
        print(f"  r={n_r:>2}: combined={best_c:>3}, random={best_r:>3}, advantage={adv:>+3}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("ИТОГ")
    print("=" * 70)


if __name__ == "__main__":
    main()
