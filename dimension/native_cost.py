"""
НАТИВНАЯ СТОИМОСТЬ: измеряем κ для каждого узла gap-зоны.

Теория нашего измерения:
  Gap cost = ∏ (1/κᵢ) для всех узлов инъекции
  κᵢ = P(метка δ проходит узел i БЕЗ увеличения)

Измеряем:
  Для каждого раунда r=18..25 (инъекция):
    κ(NODE_a, r) = P(δa[r+1] = δa[r])  — узел не усилил метку
    κ(NODE_e, r) = P(δe[r+1] = δe[r])  — узел не усилил метку

  Для saturated zone (r=25+):
    κ_sat = P(δstate[r+1] = δstate[r])  — плато, нет роста

  Общая стоимость = ∏ (1/κᵢ)
"""

import numpy as np

MASK32 = 0xFFFFFFFF
K = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
    0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
    0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
    0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
    0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
    0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
]
IV = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]

def hw(x): return bin(x).count('1')
def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
def add32(x, y): return (x + y) & MASK32
def sub32(x, y): return (x - y) & MASK32
def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)

def R(state, W_r, r_idx):
    a, b, c, d, e, f, g, h = state
    T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e, f, g)), K[r_idx]), W_r)
    T2 = add32(Sigma0(a), Maj(a, b, c))
    return (add32(T1, T2), a, b, c, add32(d, T1), e, f, g)

def expand_schedule(W16):
    W = list(W16)
    for r in range(16, 64):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))
    return W

def a_repair_W(state2, target_a, r_idx):
    a, b, c, d, e, f, g, h = state2
    T2 = add32(Sigma0(a), Maj(a, b, c))
    T1_needed = sub32(target_a, T2)
    return sub32(sub32(sub32(sub32(T1_needed, h), Sigma1(e)), Ch(e, f, g)), K[r_idx])


def main():
    np.random.seed(42)

    print("=" * 70)
    print("НАТИВНАЯ СТОИМОСТЬ: κ для каждого узла")
    print("=" * 70)

    N = 5000

    # Собираем per-round, per-register δ
    # При a-repair break=3, bit=31
    all_dreg = np.zeros((N, 65, 8))  # δ per register per round

    for trial in range(N):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W1f = expand_schedule(W1)
        s1 = tuple(IV)
        states1 = [s1]
        for r in range(64):
            s1 = R(s1, W1f[r], r)
            states1.append(s1)

        W2 = list(W1)
        W2[3] ^= (1 << 31)
        s2 = tuple(IV)
        for r in range(3):
            s2 = R(s2, W2[r], r)
        s2 = R(s2, W2[3], 3)
        for r in range(4, 16):
            W2[r] = a_repair_W(s2, states1[r + 1][0], r)
            s2 = R(s2, W2[r], r)

        W2f = expand_schedule(W2)
        s2_full = tuple(IV)
        for r in range(64):
            w = W2[r] if r < 16 else W2f[r]
            s2_full = R(s2_full, w, r)
            for i in range(8):
                all_dreg[trial, r + 1, i] = hw(states1[r + 1][i] ^ s2_full[i])

    # =================================================================
    print(f"\n{'=' * 70}")
    print("1. κ ПО УЗЛАМ: P(δ не растёт) для каждого раунда")
    print("=" * 70)

    reg_names = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']

    # κ для NODE_a: P(δa[r+1] ≤ δa[r])
    # κ для NODE_e: P(δe[r+1] ≤ δe[r])
    # Трубы: κ=1 (δb[r+1]=δa[r], всегда)

    print(f"\n  {'r':>4} {'κ(NODE_a)':>10} {'κ(NODE_e)':>10} {'κ(total)':>10} {'δstate':>8}")

    kappa_products_a = 1.0
    kappa_products_e = 1.0
    kappa_products_total = 1.0

    gap_kappas = []

    for r in range(17, 64):
        # κ(NODE_a) = P(δa[r+1] ≤ δa[r])
        da_same_or_less = sum(1 for t in range(N) if all_dreg[t, r+1, 0] <= all_dreg[t, r, 0])
        ka = da_same_or_less / N

        # κ(NODE_e) = P(δe[r+1] ≤ δe[r])
        de_same_or_less = sum(1 for t in range(N) if all_dreg[t, r+1, 4] <= all_dreg[t, r, 4])
        ke = de_same_or_less / N

        # κ(round) = P(δstate[r+1] ≤ δstate[r])
        ds_r = np.sum(all_dreg[:, r, :], axis=1)
        ds_r1 = np.sum(all_dreg[:, r+1, :], axis=1)
        kt = sum(1 for t in range(N) if ds_r1[t] <= ds_r[t]) / N

        mean_ds = np.mean(ds_r1)

        if r < 30:
            kappa_products_a *= max(ka, 1/N)
            kappa_products_e *= max(ke, 1/N)
            kappa_products_total *= max(kt, 1/N)

        gap_kappas.append((r, ka, ke, kt, mean_ds))

        if r <= 25 or r >= 60 or r in [30, 40, 50]:
            print(f"  {r:4d} {ka:9.4f} {ke:9.4f} {kt:9.4f} {mean_ds:7.1f}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("2. АЛЬТЕРНАТИВНАЯ κ: P(δ = 0) для каждого узла")
    print("=" * 70)

    # Более строгая мера: P(δ ОБНУЛЯЕТСЯ) на каждом раунде
    print(f"\n  {'r':>4} {'P(δa=0)':>10} {'P(δe=0)':>10} {'P(δstate=0)':>13} {'1/P':>10}")

    cumulative_cost = 0

    for r in range(17, 30):
        pa0 = sum(1 for t in range(N) if all_dreg[t, r+1, 0] == 0) / N
        pe0 = sum(1 for t in range(N) if all_dreg[t, r+1, 4] == 0) / N
        ps0 = sum(1 for t in range(N) if np.sum(all_dreg[t, r+1, :]) == 0) / N

        if ps0 > 0:
            cost = 1/ps0
            cost_log = np.log2(cost)
        else:
            cost = N  # lower bound
            cost_log = np.log2(N)

        cumulative_cost += cost_log

        print(f"  {r:4d} {pa0:9.4f} {pe0:9.4f} {ps0:12.4f} 2^{cost_log:.1f}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("3. НАТИВНЫЙ κ: P(δ УМЕНЬШАЕТСЯ) per round")
    print("=" * 70)

    # Вместо P(=0), считаем P(уменьшилось)
    # В нашем измерении: ткань "поглощает" δ если δ уменьшается

    print(f"\n  {'r':>4} {'P(δ↓)':>8} {'P(δ→)':>8} {'P(δ↑)':>8} {'Net':>8}  {'Ткань':>8}")

    for r in range(17, 35):
        ds_r = np.sum(all_dreg[:, r, :], axis=1)
        ds_r1 = np.sum(all_dreg[:, r+1, :], axis=1)

        p_down = sum(1 for t in range(N) if ds_r1[t] < ds_r[t]) / N
        p_same = sum(1 for t in range(N) if ds_r1[t] == ds_r[t]) / N
        p_up = sum(1 for t in range(N) if ds_r1[t] > ds_r[t]) / N

        net = p_down - p_up  # >0 = absorbing, <0 = amplifying
        tissue = "absorbs" if net > 0.05 else "neutral" if net > -0.05 else "amplifies"

        print(f"  {r:4d} {p_down:7.3f} {p_same:7.3f} {p_up:7.3f} {net:+7.3f}  {tissue:>8}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("4. НАТИВНАЯ СТОИМОСТЬ COLLISION")
    print("=" * 70)

    # Meeting: free (reboot 100%)
    meeting_cost = 1  # 1 try

    # Lock: ×4 (δW[9]=0 AND δW[10]=0)
    lock_cost = 4

    # Gap: measured by per-round κ
    # Phase 1 (injection, r=18..21): δ grows from ~8 to 128
    # Phase 2 (saturation, r=22..63): δ ≈ 128, circulating

    # Cost of injection phase: need δ NOT to grow → P ≈ product of κ
    # If δstate[17]=0 and δW[18]≠0:
    #   P(δstate[18]=0) = measured above
    #   P(δstate[19]=0 | δstate[18]=0) = ???

    # Conditional: given δstate[17]=0 (from lock), P(δstate[r]=0)?
    # Filter trials where δstate[17]=0
    locked_trials = [t for t in range(N) if np.sum(all_dreg[t, 17, :]) == 0]
    print(f"\n  Trials с δstate[17]=0: {len(locked_trials)}/{N}")

    if locked_trials:
        print(f"\n  Conditional P(δstate[r]=0 | lock до r=17):")
        for r in range(18, 30):
            p0 = sum(1 for t in locked_trials if np.sum(all_dreg[t, r, :]) == 0) / len(locked_trials)
            cost_str = f"2^{-np.log2(max(p0, 1/len(locked_trials))):.1f}" if p0 > 0 else f">2^{np.log2(len(locked_trials)):.0f}"
            print(f"    r={r}: P(δstate=0)={p0:.4f} → cost={cost_str}")

    # TOTAL
    print(f"\n  НАТИВНАЯ СТОИМОСТЬ:")
    print(f"    Meeting:  ×{meeting_cost}")
    print(f"    Lock:     ×{lock_cost}")

    # Gap cost from conditional
    if locked_trials:
        # P(δstate[64]=0 | lock до r=17)
        # = P(full collision | meeting + lock)
        p_full = sum(1 for t in locked_trials if np.sum(all_dreg[t, 64, :]) == 0) / len(locked_trials)
        p_full_str = f"2^{-np.log2(max(p_full, 1/len(locked_trials))):.1f}" if p_full > 0 else f">2^{np.log2(len(locked_trials)):.0f}"
        print(f"    Gap:      {p_full_str} (conditional on lock)")

        total = meeting_cost * lock_cost / max(p_full, 1/len(locked_trials))
        print(f"    TOTAL:    {meeting_cost} × {lock_cost} × {p_full_str}")
        print(f"              = 2^{np.log2(total):.1f} cascade sections")
    else:
        print(f"    Gap:      unmeasured (no locked trials)")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("5. НАТИВНЫЕ vs СТАНДАРТНЫЕ ЕДИНИЦЫ")
    print("=" * 70)

    print(f"""
  Наши единицы:                Стандартные:
  ─────────────────────────── ──────────────
  Секция                       Bit
  Позиция (1 из 16)            32 бита
  κ (пропускная способность)   Probability
  Cascade cost = ∏(1/κᵢ)       Birthday = 2^(n/2)

  Meeting = 0 секций            ≡ 0 бит (бесплатно в обоих)
  Lock = 4 секции               ≡ 2 бита
  Gap = ∏(1/κᵢ) секций          ≡ ??? бит ← ЗДЕСЬ РАЗНИЦА

  Если κ gap-узлов > 2^{{-4}}: gap cost < 2^{{128}}
  Если κ gap-узлов = 2^{{-4}}: gap cost = 2^{{4×94}} ≫ 2^128

  Ключ: κ НЕ КОНСТАНТА. Она зависит от СОСТОЯНИЯ δ.
  В инъекции (δ мал): κ высокая (δ может не вырасти)
  В сатурации (δ=128): κ ≈ 0.5 (монетка — может ↑ или ↓)
""")


if __name__ == "__main__":
    main()
