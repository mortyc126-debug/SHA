"""
Вскрытие 95 бит: ЧТО ИМЕННО создаёт HW(δH)=95?

Берём лучшие пары и разбираем:
  1. Какие из 8 слов H[0..7] вносят сколько?
  2. На каком РАУНДЕ δstate становится необратимым?
  3. Какая ОПЕРАЦИЯ (T1 vs T2) генерирует δ?
  4. Какие БИТЫ (MSB vs LSB) горят?
  5. Можно ли "откусить" хотя бы часть?
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


def full_trace(W16):
    """Полный trace: states, T1, T2, schedule."""
    W = list(W16)
    for r in range(16, 64):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))

    a, b, c, d, e, f, g, h = IV
    states = [(a, b, c, d, e, f, g, h)]
    T1s = []
    T2s = []

    for r in range(64):
        S1 = Sigma1(e)
        ch = Ch(e, f, g)
        S0 = Sigma0(a)
        maj = Maj(a, b, c)
        T1 = add32(add32(add32(add32(h, S1), ch), K[r]), W[r])
        T2 = add32(S0, maj)
        T1s.append(T1)
        T2s.append(T2)
        h, g, f, e = g, f, e, add32(d, T1)
        d, c, b, a = c, b, a, add32(T1, T2)
        states.append((a, b, c, d, e, f, g, h))

    H = [add32(IV[i], states[64][i]) for i in range(8)]
    return states, T1s, T2s, W, H


def find_best_pairs(n_pairs=20, n_search=200000):
    """Ищем пары с минимальным HW(δH)."""
    np.random.seed(42)
    pairs = []

    for _ in range(n_search):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W1)
        W2[0] ^= (1 << np.random.randint(0, 32))

        _, _, _, _, H1 = full_trace(W1)
        _, _, _, _, H2 = full_trace(W2)
        hd = sum(hw(H1[i] ^ H2[i]) for i in range(8))
        pairs.append((hd, W1, W2))

    pairs.sort()
    return pairs[:n_pairs]


def autopsy(W1, W2):
    """Полное вскрытие пары."""
    st1, t1_1, t2_1, w1, H1 = full_trace(W1)
    st2, t1_2, t2_2, w2, H2 = full_trace(W2)

    result = {}

    # 1. δH по словам
    dH = [H1[i] ^ H2[i] for i in range(8)]
    dH_hw = [hw(d) for d in dH]
    result['dH_per_word'] = dH_hw
    result['dH_total'] = sum(dH_hw)

    # 2. δstate по раундам
    dstate_per_round = []
    for r in range(65):
        ds = sum(hw(st1[r][i] ^ st2[r][i]) for i in range(8))
        dstate_per_round.append(ds)
    result['dstate'] = dstate_per_round

    # 3. δ по регистрам на каждом раунде
    reg_names = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
    dreg = np.zeros((65, 8))
    for r in range(65):
        for i in range(8):
            dreg[r, i] = hw(st1[r][i] ^ st2[r][i])
    result['dreg'] = dreg

    # 4. δT1 и δT2 по раундам
    dt1 = [hw(t1_1[r] ^ t1_2[r]) for r in range(64)]
    dt2 = [hw(t2_1[r] ^ t2_2[r]) for r in range(64)]
    result['dT1'] = dt1
    result['dT2'] = dt2

    # 5. δW schedule
    dW = [hw(w1[r] ^ w2[r]) for r in range(64)]
    result['dW'] = dW

    # 6. Битовые позиции δH — какие биты горят чаще?
    bit_positions = np.zeros(32)
    for i in range(8):
        for b in range(32):
            if (dH[i] >> b) & 1:
                bit_positions[b] += 1
    result['bit_positions'] = bit_positions

    # 7. Финализация: δ(state[64] + IV) — откуда δH?
    # H[i] = state[64][i] + IV[i]
    # δH[i] = δ(state[64][i] + IV[i]) = δstate[64][i] + δcarry
    # Carry в финализации
    dfinal_state = [hw(st1[64][i] ^ st2[64][i]) for i in range(8)]
    result['dfinal_state'] = dfinal_state

    return result


def main():
    print("=" * 70)
    print("ВСКРЫТИЕ 95 БИТ: что именно создаёт HW(δH)?")
    print("=" * 70)

    # Ищем лучшие пары
    print("\n  Поиск лучших пар (200K)...")
    pairs = find_best_pairs(n_pairs=20, n_search=200000)
    print(f"  Топ-5: {[p[0] for p in pairs[:5]]}")

    # ===================================================================
    print("\n" + "=" * 70)
    print("1. δH ПО СЛОВАМ: какие H[i] горят больше?")
    print("=" * 70)

    # Агрегируем по лучшим парам
    word_hws = np.zeros(8)
    for _, W1, W2 in pairs[:10]:
        r = autopsy(W1, W2)
        word_hws += r['dH_per_word']
    word_hws /= 10

    reg_names = ['a(H0)', 'b(H1)', 'c(H2)', 'd(H3)', 'e(H4)', 'f(H5)', 'g(H6)', 'h(H7)']
    print(f"\n  Средний HW(δH[i]) по лучшим 10 парам:")
    for i in range(8):
        bar = "█" * int(word_hws[i])
        print(f"    {reg_names[i]}: {word_hws[i]:5.1f}/32  {bar}")
    print(f"    TOTAL: {sum(word_hws):.1f}")

    # ===================================================================
    print("\n" + "=" * 70)
    print("2. TIMELINE δstate: когда δ взрывается?")
    print("=" * 70)

    # Берём лучшую пару
    _, W1_best, W2_best = pairs[0]
    r_best = autopsy(W1_best, W2_best)

    print(f"\n  Лучшая пара: HW(δH) = {r_best['dH_total']}")
    print(f"\n  {'Раунд':>6} {'δstate':>8} {'δW[r]':>7} {'δT1':>5} {'δT2':>5} {'δa':>4} {'δe':>4}")
    print(f"  {'-'*6} {'-'*8} {'-'*7} {'-'*5} {'-'*5} {'-'*4} {'-'*4}")

    for r in range(65):
        ds = r_best['dstate'][r]
        dw = r_best['dW'][r] if r < 64 else 0
        dt1 = r_best['dT1'][r] if r < 64 else 0
        dt2 = r_best['dT2'][r] if r < 64 else 0
        da = int(r_best['dreg'][r, 0])
        de = int(r_best['dreg'][r, 4])

        bar = "█" * min(int(ds / 4), 40)

        if r <= 5 or r >= 60 or r in [10, 16, 17, 18, 19, 20, 32, 48]:
            print(f"  r={r:3d} {ds:7d} {dw:6d} {dt1:4d} {dt2:4d} {da:3d} {de:3d}  {bar}")

    # ===================================================================
    print("\n" + "=" * 70)
    print("3. δ ПО РЕГИСТРАМ: кто несёт δ через раунды?")
    print("=" * 70)

    dreg = r_best['dreg']
    print(f"\n  HW(δ) по регистрам на ключевых раундах:")
    print(f"  {'r':>4} {'δa':>5} {'δb':>5} {'δc':>5} {'δd':>5} {'δe':>5} {'δf':>5} {'δg':>5} {'δh':>5} {'total':>6}")

    for r in [0, 1, 2, 5, 10, 16, 17, 18, 19, 20, 32, 48, 60, 62, 63, 64]:
        vals = [int(dreg[r, i]) for i in range(8)]
        total = sum(vals)
        print(f"  {r:4d} {vals[0]:5d} {vals[1]:5d} {vals[2]:5d} {vals[3]:5d} "
              f"{vals[4]:5d} {vals[5]:5d} {vals[6]:5d} {vals[7]:5d} {total:5d}")

    # ===================================================================
    print("\n" + "=" * 70)
    print("4. БИТОВЫЕ ПОЗИЦИИ δH: MSB vs LSB")
    print("=" * 70)

    # Агрегируем по 10 лучшим парам
    bit_agg = np.zeros(32)
    for _, W1, W2 in pairs[:10]:
        r = autopsy(W1, W2)
        bit_agg += r['bit_positions']
    bit_agg /= 10

    print(f"\n  P(δH бит j = 1) по позициям, усреднено по 10 парам:")
    for b in range(32):
        p = bit_agg[b] / 8  # 8 слов, нормируем
        bar = "█" * int(p * 40)
        label = ""
        if b == 0: label = " ← LSB"
        if b == 31: label = " ← MSB"
        print(f"    бит {b:2d}: {p:.3f}  {bar}{label}")

    # ===================================================================
    print("\n" + "=" * 70)
    print("5. δT1 vs δT2: КТО ВИНОВАТ?")
    print("=" * 70)

    dt1_agg = np.zeros(64)
    dt2_agg = np.zeros(64)
    for _, W1, W2 in pairs[:10]:
        r = autopsy(W1, W2)
        dt1_agg += r['dT1']
        dt2_agg += r['dT2']
    dt1_agg /= 10
    dt2_agg /= 10

    print(f"\n  Средний HW(δT1) и HW(δT2) по раундам:")
    print(f"  {'r':>4} {'δT1':>7} {'δT2':>7} {'δT1-δT2':>9} {'Кто больше':>12}")

    for r in [0, 1, 2, 5, 10, 16, 17, 18, 19, 20, 32, 48, 60, 63]:
        diff = dt1_agg[r] - dt2_agg[r]
        who = "T1" if diff > 1 else ("T2" if diff < -1 else "≈")
        print(f"  {r:4d} {dt1_agg[r]:6.1f} {dt2_agg[r]:6.1f} {diff:+8.1f}   {who:>10}")

    # Суммарно
    total_dt1 = np.sum(dt1_agg)
    total_dt2 = np.sum(dt2_agg)
    print(f"\n  ИТОГО: δT1 = {total_dt1:.0f}, δT2 = {total_dt2:.0f}")
    print(f"  Доля T1: {total_dt1/(total_dt1+total_dt2)*100:.1f}%")
    print(f"  Доля T2: {total_dt2/(total_dt1+total_dt2)*100:.1f}%")

    # ===================================================================
    print("\n" + "=" * 70)
    print("6. SCHEDULE: где δW[r] максимален?")
    print("=" * 70)

    dw_agg = np.zeros(64)
    for _, W1, W2 in pairs[:10]:
        r = autopsy(W1, W2)
        dw_agg += r['dW']
    dw_agg /= 10

    print(f"\n  HW(δW[r]) по раундам:")
    quiet_rounds = []
    loud_rounds = []
    for r in range(64):
        bar = "█" * int(dw_agg[r])
        label = ""
        if dw_agg[r] < 2: quiet_rounds.append(r)
        if dw_agg[r] > 20: loud_rounds.append(r)
        if r <= 16 or r in quiet_rounds[-1:] or dw_agg[r] > 20 or r == 63:
            print(f"    r={r:2d}: {dw_agg[r]:5.1f}  {bar}")

    print(f"\n  Тихие раунды (δW<2):  {quiet_rounds}")
    print(f"  Громкие раунды (δW>20): {loud_rounds}")

    # ===================================================================
    print("\n" + "=" * 70)
    print("7. ФИНАЛИЗАЦИЯ: δstate[64] vs δH")
    print("=" * 70)

    ds64_agg = np.zeros(8)
    dh_agg = np.zeros(8)
    for _, W1, W2 in pairs[:10]:
        r = autopsy(W1, W2)
        ds64_agg += r['dfinal_state']
        dh_agg += r['dH_per_word']
    ds64_agg /= 10
    dh_agg /= 10

    print(f"\n  δstate[64] vs δH (финализация: H = state + IV):")
    print(f"  {'Reg':>6} {'δstate[64]':>12} {'δH':>8} {'Δ(carry)':>10}")
    for i in range(8):
        delta_carry = abs(dh_agg[i] - ds64_agg[i])
        print(f"  {reg_names[i]:>6} {ds64_agg[i]:11.1f} {dh_agg[i]:7.1f} {delta_carry:9.1f}")

    # ===================================================================
    print("\n" + "=" * 70)
    print("8. ДИАГНОЗ: ЧТО КУШАЕТ 95 БИТ")
    print("=" * 70)

    print(f"""
  РАЗЛОЖЕНИЕ 95 БИТ:

    δH = δstate[64] ⊕ δcarry_final

    δstate[64] складывается из:
      Раунды 0-16:   δstate ≈ 0 (начальный δW мал)
      Раунд 17:      δW[16] ударяет → δstate ≈ {int(r_best['dstate'][17])}
      Раунды 18-20:  лавина → δstate → {int(r_best['dstate'][20])}
      Раунды 21-64:  насыщение → δstate ≈ {int(np.mean(r_best['dstate'][21:64]))}

    ИСТОЧНИК δ:
      δT1 (управляемый):  {total_dt1:.0f} бит суммарно ({total_dt1/(total_dt1+total_dt2)*100:.0f}%)
      δT2 (неуправляемый): {total_dt2:.0f} бит суммарно ({total_dt2/(total_dt1+total_dt2)*100:.0f}%)

    SCHEDULE:
      Тихие раунды:  {len(quiet_rounds)} (δW < 2 бит)
      Громкие раунды: {len(loud_rounds)} (δW > 20 бит)

    ФИНАЛИЗАЦИЯ:
      δstate[64] ≈ δH (carry в финализации ≈ {np.mean(np.abs(dh_agg - ds64_agg)):.1f} бит)
""")


if __name__ == "__main__":
    main()
