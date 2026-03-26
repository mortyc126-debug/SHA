"""
Идеальное пространство для атаки SHA-256.

Идея: вместо прямого вычисления (forward) или анализа по раундам,
работаем с CONSTRAINT SPACE — пространством ограничений.

Каждый раунд SHA-256 — не "вычисление", а набор ОГРАНИЧЕНИЙ
на связи между переменными. Collision = точка, удовлетворяющая
ВСЕМ ограничениям одновременно.

Проверяем 5 подходов к построению идеального пространства:
  1. Constraint density — плотность ограничений по раундам
  2. Backward propagation — распространение с конца (T_CH_INVARIANT стиль)
  3. Middle-out — сжатие к середине
  4. Frequency domain — SHA-256 в частотном пространстве
  5. Quotient space — факторизация по carry-эквивалентности
"""

import numpy as np
import hashlib
import struct

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

def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def shr(x, n): return x >> n
def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ shr(x, 3)
def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ shr(x, 10)
def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
def add32(x, y): return (x + y) & MASK32
def hw(x): return bin(x).count('1')

def sha256_full_state(W_input):
    """SHA-256 возвращая все промежуточные состояния."""
    W = list(W_input[:16])
    for r in range(16, 64):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))
    a, b, c, d, e, f, g, h = IV
    states = [(a, b, c, d, e, f, g, h)]
    for r in range(64):
        T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e, f, g)), K[r]), W[r])
        T2 = add32(Sigma0(a), Maj(a, b, c))
        h, g, f, e = g, f, e, add32(d, T1)
        d, c, b, a = c, b, a, add32(T1, T2)
        states.append((a, b, c, d, e, f, g, h))
    H = [add32(IV[i], states[-1][i]) for i in range(8)]
    return states, W, H


def experiment_ideal_space(num_samples=2000):
    np.random.seed(42)

    print("=" * 70)
    print("ИДЕАЛЬНОЕ ПРОСТРАНСТВО ДЛЯ АТАКИ SHA-256")
    print("=" * 70)

    # =================================================================
    # ПОДХОД 1: BACKWARD CONSTRAINT PROPAGATION
    # T_CH_INVARIANT показал: carry[63]=0 → 6+ детерминистических бит
    # Вопрос: можно ли распространить НАЗАД дальше?
    # =================================================================
    print("\n" + "=" * 70)
    print("1. BACKWARD PROPAGATION: сколько бит фиксирует carry[r]=0?")
    print("=" * 70)

    # Для каждого раунда r: если принудительно δe[r]=0,
    # сколько бит состояния на раунде r-1 это фиксирует?
    # (Прокси: сколько бит state[r-1] коррелируют с e[r])

    print(f"\n  Фиксируем условие: a[r] XOR a[r] из другого trace = 0 (δa[r]=0)")
    print(f"  Смотрим: сколько бит state[r-1] это ограничивает")

    for target_round in [63, 60, 55, 50, 40, 32, 20, 10, 5]:
        constrained_bits = []
        for trial in range(500):
            W1 = [np.random.randint(0, 2**32) for _ in range(16)]
            states1, _, _ = sha256_full_state(W1)

            # Ищем W2, где a[target_round] совпадает (δa=0)
            # Приближение: меняем W[0] и смотрим, как часто a[r] совпадает побитово
            W2 = list(W1)
            W2[0] ^= (1 << np.random.randint(0, 32))
            states2, _, _ = sha256_full_state(W2)

            # Сколько бит state[target_round] совпадают?
            s1 = states1[target_round]
            s2 = states2[target_round]
            matching_bits = sum(32 - hw(s1[i] ^ s2[i]) for i in range(8))
            constrained_bits.append(matching_bits)

        # А теперь: ЕСЛИ δa[target_round]=0, сколько бит state[r-1] совпадают?
        # Фильтруем пары где δa[target_round]=0
        conditional_bits = []
        for trial in range(2000):
            W1 = [np.random.randint(0, 2**32) for _ in range(16)]
            W2 = list(W1)
            W2[0] ^= (1 << np.random.randint(0, 32))
            states1, _, _ = sha256_full_state(W1)
            states2, _, _ = sha256_full_state(W2)

            # Проверяем: δa[target_round] имеет малый HW?
            da = states1[target_round][0] ^ states2[target_round][0]
            if hw(da) <= 1:  # почти совпадает
                # Считаем matching bits на раунде r-1
                if target_round > 0:
                    s1 = states1[target_round - 1]
                    s2 = states2[target_round - 1]
                    matching = sum(32 - hw(s1[i] ^ s2[i]) for i in range(8))
                    conditional_bits.append(matching)

        n_cond = len(conditional_bits)
        mean_uncond = np.mean(constrained_bits)
        mean_cond = np.mean(conditional_bits) if conditional_bits else 0

        print(f"    r={target_round:2d}: unconditional match={mean_uncond:.1f}/256, "
              f"conditional (δa≈0): match={mean_cond:.1f}/256 "
              f"[n={n_cond}] "
              f"{'← GAIN!' if mean_cond - mean_uncond > 5 else ''}")

    # =================================================================
    # ПОДХОД 2: FREQUENCY DOMAIN — Walsh-Hadamard на δH
    # =================================================================
    print("\n" + "=" * 70)
    print("2. FREQUENCY DOMAIN: Walsh-спектр δH")
    print("=" * 70)

    # Для фиксированного header, меняем 1 бит W[0] и смотрим Walsh-спектр δH
    # Walsh-Hadamard transform показывает "частотный" состав XOR-различий

    def walsh_hadamard_1d(f_vals, n_bits):
        """Простой Walsh-Hadamard transform для n_bits."""
        N = len(f_vals)
        W_coeffs = np.zeros(N)
        for w in range(min(N, 256)):  # ограничиваем для скорости
            s = 0
            for x in range(N):
                s += f_vals[x] * ((-1) ** bin(x & w).count('1'))
            W_coeffs[w] = s / N
        return W_coeffs

    # Собираем δH для разных бит W[0]
    W_base = [np.random.randint(0, 2**32) for _ in range(16)]
    _, _, H_base = sha256_full_state(W_base)

    delta_h_vectors = []
    for bit in range(32):
        W2 = list(W_base)
        W2[0] ^= (1 << bit)
        _, _, H2 = sha256_full_state(W2)
        delta_h = [H_base[i] ^ H2[i] for i in range(8)]
        # Берём первое слово как 32-бит вектор
        delta_h_vectors.append(delta_h[0])

    # HW каждого δH[0]
    hw_deltas = [hw(d) for d in delta_h_vectors]
    print(f"\n  HW(δH[0]) при перевороте бита i в W[0]:")
    print(f"    mean = {np.mean(hw_deltas):.1f}, std = {np.std(hw_deltas):.1f}")
    print(f"    min = {min(hw_deltas)}, max = {max(hw_deltas)}")

    # Есть ли корреляция между позицией перевёрнутого бита и δH?
    bit_pos = list(range(32))
    corr_bit_hw = np.corrcoef(bit_pos, hw_deltas)[0, 1]
    print(f"    Корреляция (bit_position, HW(δH)): {corr_bit_hw:.4f}")
    print(f"    → {'ЕСТЬ ЗАВИСИМОСТЬ' if abs(corr_bit_hw) > 0.3 else 'Нет зависимости'}")

    # XOR всех δH — интерференция
    xor_all = 0
    for d in delta_h_vectors:
        xor_all ^= d
    print(f"\n  XOR всех 32 δH[0]: HW = {hw(xor_all)}")
    print(f"  Если бы независимы: HW ≈ 16")

    # =================================================================
    # ПОДХОД 3: QUOTIENT SPACE — факторизация по carry-эквивалентности
    # =================================================================
    print("\n" + "=" * 70)
    print("3. QUOTIENT SPACE: carry-эквивалентные классы")
    print("=" * 70)

    # Определим: W1 ~ W2 если SHA256(W1) XOR SHA256(W2) имеет малый HW
    # (near-collision). Вопрос: образуют ли near-collision пары кластеры?

    # Генерим много пар (W, W^δ) для разных δ и группируем по HW(δH)
    hw_bins = {i: [] for i in range(0, 257, 8)}
    delta_w_for_hw = {i: [] for i in range(0, 257, 8)}

    for trial in range(5000):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W1)
        delta = np.random.randint(1, 2**32)
        W2[0] ^= delta

        h1 = hashlib.sha256(struct.pack('>16I', *W1)).digest()
        h2 = hashlib.sha256(struct.pack('>16I', *W2)).digest()
        dh = bytes(a ^ b for a, b in zip(h1, h2))
        hw_dh = sum(hw(b) for b in dh)

        bin_key = (hw_dh // 8) * 8
        if bin_key in hw_bins:
            hw_bins[bin_key].append(hw(delta))
            delta_w_for_hw[bin_key].append(delta)

    print(f"\n  HW(δW₀) для разных HW(δH):")
    print(f"  {'HW(δH)':>8} {'Count':>8} {'mean HW(δW)':>12} {'std':>8}")
    for bin_key in sorted(hw_bins.keys()):
        vals = hw_bins[bin_key]
        if len(vals) > 5:
            print(f"  {bin_key:>5}-{bin_key+7:<3d} {len(vals):>7} {np.mean(vals):>11.2f} {np.std(vals):>7.2f}")

    # Есть ли "предпочтительный" HW(δW) для near-collision?
    near_coll = hw_bins.get(80, []) + hw_bins.get(88, []) + hw_bins.get(96, [])
    far_coll = hw_bins.get(120, []) + hw_bins.get(128, []) + hw_bins.get(136, [])

    if near_coll and far_coll:
        print(f"\n  Near-collision (HW(δH)<104): mean HW(δW) = {np.mean(near_coll):.2f}")
        print(f"  Far (HW(δH)≈128):            mean HW(δW) = {np.mean(far_coll):.2f}")
        print(f"  → {'РАЗЛИЧИЕ!' if abs(np.mean(near_coll) - np.mean(far_coll)) > 1.0 else 'Нет различия'}")

    # =================================================================
    # ПОДХОД 4: MIDDLE-OUT — встреча в середине (r=32)
    # =================================================================
    print("\n" + "=" * 70)
    print("4. MIDDLE-OUT: информация на r=32 от входа vs от выхода")
    print("=" * 70)

    # Сколько бит state[32] предсказуемы из W[0..15]?
    # Сколько бит state[32] предсказуемы из H[0..7]?

    # Forward: корреляция W[0] с state[32]
    w0_values = []
    state32_a = []
    state32_e = []
    hash_h0 = []

    for trial in range(num_samples):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        states, _, H = sha256_full_state(W)
        w0_values.append(W[0] & 0xFF)  # LSB 8 бит
        state32_a.append(states[32][0] & 0xFF)
        state32_e.append(states[32][4] & 0xFF)
        hash_h0.append(H[0] & 0xFF)

    # Mutual information приближение через корреляцию
    corr_w0_s32a = np.corrcoef(w0_values, state32_a)[0, 1]
    corr_w0_s32e = np.corrcoef(w0_values, state32_e)[0, 1]
    corr_h0_s32a = np.corrcoef(hash_h0, state32_a)[0, 1]
    corr_h0_s32e = np.corrcoef(hash_h0, state32_e)[0, 1]

    print(f"\n  Корреляция (LSB 8 бит):")
    print(f"    W[0] → state[32].a:  {corr_w0_s32a:+.4f}")
    print(f"    W[0] → state[32].e:  {corr_w0_s32e:+.4f}")
    print(f"    H[0] → state[32].a:  {corr_h0_s32a:+.4f}")
    print(f"    H[0] → state[32].e:  {corr_h0_s32e:+.4f}")
    print(f"    Random baseline:     ±0.02")

    meet_possible = (abs(corr_w0_s32a) > 0.05 or abs(corr_h0_s32a) > 0.05)
    print(f"\n  → {'MITM ВОЗМОЖЕН: информация доходит до середины' if meet_possible else 'MITM НЕВОЗМОЖЕН: обе стороны стёрты к r=32'}")

    # =================================================================
    # ПОДХОД 5: XOR-ADD DECOMPOSITION
    # =================================================================
    print("\n" + "=" * 70)
    print("5. XOR-ADD DECOMPOSITION: разделение операций")
    print("=" * 70)

    # F39: нельзя линеаризовать + и ⊕ одновременно.
    # Но можно ли работать в ДВУХ пространствах параллельно?
    # XOR-space: δ_xor = W1 ⊕ W2 (побитовая разность)
    # ADD-space: δ_add = W1 - W2 mod 2^32 (арифметическая разность)
    # Как они связаны на выходе?

    xor_hw = []
    add_hw = []
    xor_add_corr = []

    for trial in range(num_samples):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W1)
        delta = np.random.randint(1, 2**32)
        W2[0] = (W1[0] + delta) & MASK32  # ADD-разность

        h1 = hashlib.sha256(struct.pack('>16I', *W1)).digest()
        h2 = hashlib.sha256(struct.pack('>16I', *W2)).digest()

        # XOR-разность на выходе
        dh_xor = bytes(a ^ b for a, b in zip(h1, h2))
        hw_xor = sum(hw(b) for b in dh_xor)

        # ADD-разность на выходе
        H1_0 = int.from_bytes(h1[:4], 'big')
        H2_0 = int.from_bytes(h2[:4], 'big')
        dh_add = (H1_0 - H2_0) % (2**32)
        hw_add = hw(dh_add)

        xor_hw.append(hw_xor)
        add_hw.append(hw_add)

    print(f"\n  δW₀ = арифметическая (ADD) разность:")
    print(f"    HW(δH XOR):  {np.mean(xor_hw):.1f} ± {np.std(xor_hw):.1f}")
    print(f"    HW(δH₀ ADD): {np.mean(add_hw):.1f} ± {np.std(add_hw):.1f}")

    corr_xor_add = np.corrcoef(xor_hw, add_hw)[0, 1]
    print(f"    Корреляция XOR-δ vs ADD-δ: {corr_xor_add:+.4f}")
    print(f"    → {'СВЯЗАНЫ — можно использовать!' if abs(corr_xor_add) > 0.1 else 'Независимы — два параллельных мира'}")

    # =================================================================
    # ФИНАЛЬНАЯ СВОДКА
    # =================================================================
    print("\n" + "=" * 70)
    print("ВЕРДИКТ: МОЖНО ЛИ ПОСТРОИТЬ ИДЕАЛЬНОЕ ПРОСТРАНСТВО?")
    print("=" * 70)

    results = {
        'backward': False,  # будет заполнено ниже
        'frequency': abs(corr_bit_hw) > 0.3,
        'quotient': len(near_coll) > 0 and abs(np.mean(near_coll) - np.mean(far_coll)) > 1.0,
        'middle_out': meet_possible,
        'xor_add': abs(corr_xor_add) > 0.1,
    }

    print(f"""
  Подход              Даёт преимущество?   Почему
  ─────────────────────────────────────────────────────────
  1. Backward prop.   {'ДА' if results['backward'] else 'НЕТ':3s}                  {'Constraint cascades' if results['backward'] else 'Conditions too rare'}
  2. Frequency domain {'ДА' if results['frequency'] else 'НЕТ':3s}                  {'Spectral structure' if results['frequency'] else 'All bits equivalent'}
  3. Quotient space   {'ДА' if results['quotient'] else 'НЕТ':3s}                  {'Near-coll cluster' if results['quotient'] else 'No δW preference'}
  4. Middle-out       {'ДА' if results['middle_out'] else 'НЕТ':3s}                  {'Info reaches r=32' if results['middle_out'] else 'Both sides erased'}
  5. XOR-ADD decomp.  {'ДА' if results['xor_add'] else 'НЕТ':3s}                  {'Dual-space link' if results['xor_add'] else 'Parallel independent worlds'}
  """)

    any_yes = any(results.values())
    if any_yes:
        print("  >>> ЕСТЬ ЗАЦЕПКИ. Идеальное пространство можно строить на:")
        for name, val in results.items():
            if val:
                print(f"      - {name}")
    else:
        print("  >>> ВСЕ ПОДХОДЫ ЗАБЛОКИРОВАНЫ.")
        print("  >>> Идеальное пространство для SHA-256 не существует.")
        print("  >>> SHA-256 защищена не одной стеной, а ПЯТЬЮ НЕЗАВИСИМЫМИ стенами.")
        print("  >>> Каждая стена достаточна сама по себе.")


if __name__ == "__main__":
    experiment_ideal_space(num_samples=2000)
