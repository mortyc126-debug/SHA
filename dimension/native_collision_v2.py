"""
Тканевая коллизия v2: начинаем с H[0]-only (32 бит, 2^16 birthday).
Потом расширяем.
"""

import numpy as np
import struct, hashlib, time

def sha256_words(W16):
    h = hashlib.sha256(struct.pack('>16I', *W16)).digest()
    return struct.unpack('>8I', h)

def hw(x): return bin(x).count('1')

def main():
    np.random.seed(42)

    print("=" * 70)
    print("ТКАНЕВАЯ КОЛЛИЗИЯ v2: поэтапное наращивание")
    print("=" * 70)

    for target_words, n_max in [
        ([0], 500000),           # 32 бит, 2^16
        ([4], 500000),           # 32 бит, 2^16
        ([0, 4], 10000000),      # 64 бит, 2^32 — долго
        ([7], 500000),           # 32 бит, pipe word
    ]:
        label = "+".join(f"H[{w}]" for w in target_words)
        bits = len(target_words) * 32
        expected_N = int(2 ** (bits / 2) * 1.2)

        print(f"\n  --- {label} ({bits} бит, need ≈{expected_N:,}) ---")

        h_dict = {}
        collisions = []
        t0 = time.time()

        for i in range(min(n_max, expected_N * 10)):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            H = sha256_words(W)
            key = tuple(H[w] for w in target_words)

            if key in h_dict:
                W_prev, H_prev = h_dict[key]
                if W != W_prev:
                    collisions.append((W, H, W_prev, H_prev))
                    if len(collisions) >= 50:
                        break
            else:
                h_dict[key] = (W, H)

        elapsed = time.time() - t0
        print(f"    {i+1:,} хешей, {len(collisions)} коллизий, {elapsed:.1f}с")

        if not collisions:
            print(f"    Нет коллизий (нужно больше)")
            continue

        # Анализ pipe-residual
        reg = ['H[0]', 'H[1]', 'H[2]', 'H[3]', 'H[4]', 'H[5]', 'H[6]', 'H[7]']
        types = ['NODE', 'PIPE', 'PIPE', 'PIPE', 'NODE', 'PIPE', 'PIPE', 'PIPE']

        per_word = {k: [] for k in range(8)}
        for W1, H1, W2, H2 in collisions:
            for k in range(8):
                per_word[k].append(hw(H1[k] ^ H2[k]))

        print(f"\n    {'Слово':<8} {'HW(δ)':>8} {'Тип':>6}")
        for k in range(8):
            m = np.mean(per_word[k])
            marker = " ← ZERO" if m < 0.5 else ""
            print(f"    {reg[k]:<8} {m:7.2f}  {types[k]:>5}{marker}")

        # Корреляции между pipe-словами
        if len(collisions) >= 10:
            # a-chain: H[3]→H[2]→H[1]
            c32 = np.corrcoef(per_word[3], per_word[2])[0, 1] if len(per_word[3]) > 3 else 0
            c21 = np.corrcoef(per_word[2], per_word[1])[0, 1] if len(per_word[2]) > 3 else 0
            c76 = np.corrcoef(per_word[7], per_word[6])[0, 1] if len(per_word[7]) > 3 else 0
            c65 = np.corrcoef(per_word[6], per_word[5])[0, 1] if len(per_word[6]) > 3 else 0

            print(f"\n    Chain correlations:")
            print(f"      a: H[3]↔H[2]={c32:+.3f}  H[2]↔H[1]={c21:+.3f}")
            print(f"      e: H[7]↔H[6]={c76:+.3f}  H[6]↔H[5]={c65:+.3f}")

            any_sig = max(abs(c32), abs(c21), abs(c76), abs(c65)) > 0.15
            print(f"      → {'СТРУКТУРА!' if any_sig else 'Нет структуры'}")

            # Cross-chain: H[3]↔H[7] (a[61]↔e[61])
            c37 = np.corrcoef(per_word[3], per_word[7])[0, 1] if len(per_word[3]) > 3 else 0
            print(f"      Cross: H[3]↔H[7]={c37:+.3f} (a[61]↔e[61])")

    print(f"\n{'=' * 70}")
    print("ВЫВОД")
    print("=" * 70)


if __name__ == "__main__":
    main()
