"""
ТКАНЕВАЯ КОЛЛИЗИЯ — нативная операция нашего измерения.

Старая математика: collision = δH[0..7] = 0 (256 бит, 2^128)
Наше измерение:    fabric collision = δH[0]=δH[4]=0 (64 бит, 2^32)

H[0] и H[4] — NODE-слова (проходят через узлы, κ≈0.05)
H[1,2,3,5,6,7] — PIPE-слова (трубы, κ=1, структурированы)

Эксперимент:
  1. Найти тканевую коллизию (δH[0]=δH[4]=0) — birthday 2^32
  2. Анализировать pipe-residual (δH[1,2,3,5,6,7])
  3. Pipe-residual структурирован? Можно сжать?
"""

import numpy as np
import struct
import hashlib
import time

def sha256_words(W16):
    h = hashlib.sha256(struct.pack('>16I', *W16)).digest()
    return struct.unpack('>8I', h)

def hw(x): return bin(x).count('1')

MASK32 = 0xFFFFFFFF


def find_fabric_collisions(n_target=50):
    """Birthday на H[0]||H[4] (64 бит). Ожидаем ~2^32 попыток."""
    print(f"  Ищем тканевые коллизии (δH[0]=δH[4]=0)...")

    # Для ускорения: используем multi-target birthday
    # Генерим много хешей, группируем по (H[0], H[4])
    np.random.seed(42)
    h04_dict = {}
    collisions = []
    t0 = time.time()

    batch = 0
    while len(collisions) < n_target:
        batch += 1
        N = 500000
        for i in range(N):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            H = sha256_words(W)
            key = (H[0], H[4])

            if key in h04_dict:
                W_prev, H_prev = h04_dict[key]
                if W != W_prev:
                    collisions.append((W, H, W_prev, H_prev))
                    if len(collisions) >= n_target:
                        break
            else:
                h04_dict[key] = (W, H)

        elapsed = time.time() - t0
        total_hashes = batch * N
        print(f"    Batch {batch}: {total_hashes:,} хешей, {len(collisions)} коллизий ({elapsed:.1f}с)")

        if batch > 20:
            break

    return collisions


def analyze_fabric_collisions(collisions):
    """Анализируем pipe-residual в тканевых коллизиях."""

    print(f"\n  Анализ {len(collisions)} тканевых коллизий:")
    print(f"\n  {'Слово':<10} {'Mean HW(δ)':>12} {'Std':>8} {'Min':>6} {'Max':>6} {'Тип':>8}")
    print(f"  {'-'*10} {'-'*12} {'-'*8} {'-'*6} {'-'*6} {'-'*8}")

    reg = ['H[0]', 'H[1]', 'H[2]', 'H[3]', 'H[4]', 'H[5]', 'H[6]', 'H[7]']
    types = ['NODE', 'PIPE', 'PIPE', 'PIPE', 'NODE', 'PIPE', 'PIPE', 'PIPE']

    per_word = [[] for _ in range(8)]
    for W1, H1, W2, H2 in collisions:
        for k in range(8):
            per_word[k].append(hw(H1[k] ^ H2[k]))

    for k in range(8):
        m = np.mean(per_word[k])
        s = np.std(per_word[k])
        mn = min(per_word[k])
        mx = max(per_word[k])
        print(f"  {reg[k]:<10} {m:11.2f} {s:7.2f} {mn:5d} {mx:5d} {types[k]:>8}")

    # Pipe-residual total
    pipe_totals = []
    for W1, H1, W2, H2 in collisions:
        pipe_hw = sum(hw(H1[k] ^ H2[k]) for k in [1,2,3,5,6,7])
        pipe_totals.append(pipe_hw)

    print(f"\n  Pipe-residual (δH[1,2,3,5,6,7]):")
    print(f"    Mean: {np.mean(pipe_totals):.1f} / 192")
    print(f"    Std:  {np.std(pipe_totals):.1f}")
    print(f"    Min:  {min(pipe_totals)}")

    return per_word, pipe_totals


def analyze_pipe_structure(collisions):
    """Pipe-residual СТРУКТУРИРОВАН? (трубная теория)."""
    print(f"\n  Pipe-residual: δH[1]=δa[63], δH[2]=δa[62], δH[3]=δa[61]")
    print(f"  Pipe-residual: δH[5]=δe[63], δH[6]=δe[62], δH[7]=δe[61]")
    print(f"\n  Если структурирован: δH[3] должен ПРЕДСКАЗЫВАТЬ δH[2] и δH[1]")
    print(f"  (потому что a[61]→a[62]→a[63] через NODE_a)")

    # Корреляция между pipe-словами
    data = {k: [] for k in range(8)}
    for W1, H1, W2, H2 in collisions:
        for k in range(8):
            data[k].append(hw(H1[k] ^ H2[k]))

    print(f"\n  Корреляция HW(δH[i]) ↔ HW(δH[j]):")
    pipe_indices = [1, 2, 3, 5, 6, 7]
    pipe_names = ['H1(a63)', 'H2(a62)', 'H3(a61)', 'H5(e63)', 'H6(e62)', 'H7(e61)']

    corr_data = np.array([[data[k][i] for k in pipe_indices] for i in range(len(collisions))])
    if len(collisions) > 5:
        corr = np.corrcoef(corr_data.T)
        print(f"  {'':>10}", end="")
        for name in pipe_names:
            print(f" {name:>9}", end="")
        print()
        for i, name in enumerate(pipe_names):
            print(f"  {name:>10}", end="")
            for j in range(len(pipe_names)):
                if i == j:
                    print(f"     ---  ", end="")
                else:
                    c = corr[i, j]
                    marker = "*" if abs(c) > 0.2 else " "
                    print(f"  {c:+.3f}{marker}", end="")
            print()

    # Pipe-chain: δH[3]=δa[61], и a[62]=NODE_a(61), a[63]=NODE_a(62)
    # Если NODE_a почти линеен: δa[62] ≈ f(δa[61]), δa[63] ≈ f(δa[62])
    # → δH[2] ≈ f(δH[3]), δH[1] ≈ f(δH[2])

    # Проверяем: HW(δH[3]) предсказывает HW(δH[2])?
    if len(collisions) > 10:
        corr_32 = np.corrcoef(data[3], data[2])[0, 1]
        corr_21 = np.corrcoef(data[2], data[1])[0, 1]
        corr_76 = np.corrcoef(data[7], data[6])[0, 1]
        corr_65 = np.corrcoef(data[6], data[5])[0, 1]

        print(f"\n  Chain correlations:")
        print(f"    a-chain: H[3]→H[2]: {corr_32:+.3f}  H[2]→H[1]: {corr_21:+.3f}")
        print(f"    e-chain: H[7]→H[6]: {corr_76:+.3f}  H[6]→H[5]: {corr_65:+.3f}")

        has_structure = abs(corr_32) > 0.15 or abs(corr_76) > 0.15
        print(f"\n  >>> {'PIPE-RESIDUAL СТРУКТУРИРОВАН!' if has_structure else 'Pipe-residual случаен'}")

        if has_structure:
            # Если структурирован: эффективный размер < 192 бит
            # Оценка через PCA
            from numpy.linalg import eigvalsh
            cov = np.cov(corr_data.T)
            eigs = sorted(eigvalsh(cov), reverse=True)
            total_var = sum(eigs)
            cum_var = np.cumsum(eigs) / total_var
            dim_90 = np.searchsorted(cum_var, 0.90) + 1
            dim_95 = np.searchsorted(cum_var, 0.95) + 1

            print(f"\n  PCA pipe-residual:")
            print(f"    Dimension for 90% variance: {dim_90} (из 6)")
            print(f"    Dimension for 95% variance: {dim_95} (из 6)")

            effective_bits = dim_90 * 32
            print(f"    Effective dimensionality: {effective_bits} бит (из 192)")
            print(f"    Compression: {192 - effective_bits} бит!")
            return effective_bits

    return 192  # no compression


def experiment():
    np.random.seed(42)

    print("=" * 70)
    print("ТКАНЕВАЯ КОЛЛИЗИЯ — нативная операция нашего измерения")
    print("Fabric collision: δH[0]=δH[4]=0 (NODE-слова)")
    print("=" * 70)

    # Этап 1: найти тканевые коллизии
    collisions = find_fabric_collisions(n_target=30)

    if len(collisions) < 5:
        print("  Недостаточно коллизий для анализа")
        return

    # Этап 2: анализ pipe-residual
    per_word, pipe_totals = analyze_fabric_collisions(collisions)

    # Этап 3: структура pipe-residual
    effective_bits = analyze_pipe_structure(collisions)

    # Этап 4: стоимость в нашем измерении
    print(f"\n{'=' * 70}")
    print("СТОИМОСТЬ В НАШЕМ ИЗМЕРЕНИИ")
    print("=" * 70)

    print(f"""
  ТКАНЕВАЯ КОЛЛИЗИЯ (δH[0]=δH[4]=0):
    Стоимость: birthday на 64 бит = 2^32
    Pipe-residual: {np.mean(pipe_totals):.0f} бит из 192
    Effective dimensions: {effective_bits} бит

  ПОЛНАЯ COLLISION (δH[0..7]=0):
    Из тканевой: нужно обнулить pipe-residual
    Pipe-residual: {effective_bits} эффективных бит
    Birthday на {effective_bits} бит: 2^{effective_bits // 2}
    + стоимость тканевой: 2^32

    ИТОГО: 2^{max(effective_bits // 2, 32)}

  СТАНДАРТ: 2^128

  ВЫИГРЫШ: {128 - max(effective_bits // 2, 32)} бит
  (если pipe-residual структурирован)
""")


if __name__ == "__main__":
    experiment()
