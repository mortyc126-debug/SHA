"""
АНОМАЛИЯ 4: Bias W[7][15] → H[0][31] = 0.5345

Факт: бит 15 слова W[7] коррелирует с битом 31 слова H[0]
с вероятностью 0.5345 вместо 0.500. Bias = 3.45%.

Вопросы:
  A. Воспроизводится ли?
  B. Это единственная такая пара или есть ещё?
  C. Если да — это свойство SHA-256 или артефакт малой выборки?
  D. Масштаб: сколько нужно samples для уверенности?
"""

import numpy as np
import struct, hashlib

MASK32 = 0xFFFFFFFF

def sha256_hash(W16):
    raw = struct.pack('>16I', *W16)
    return struct.unpack('>8I', hashlib.sha256(raw).digest())


def measure_bias(in_word, in_bit, out_word, out_bit, n_samples):
    """P(H[out_word][out_bit] = W[in_word][in_bit])."""
    match = 0
    for _ in range(n_samples):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H = sha256_hash(W)
        ib = (W[in_word] >> in_bit) & 1
        ob = (H[out_word] >> out_bit) & 1
        if ib == ob:
            match += 1
    return match / n_samples


def main():
    np.random.seed(42)

    print("=" * 70)
    print("АНОМАЛИЯ 4: Bias W[7][15] -> H[0][31] = 0.5345")
    print("=" * 70)

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("A. Воспроизводим с разными N")
    print("=" * 70)

    for N in [1000, 5000, 10000, 50000, 200000]:
        bias = measure_bias(7, 15, 0, 31, N)
        # 95% CI for binomial: ±1.96*sqrt(p(1-p)/N)
        ci = 1.96 * np.sqrt(0.25 / N)
        sig = abs(bias - 0.5) / (np.sqrt(0.25 / N))
        print(f"  N={N:>7}: P = {bias:.4f} (±{ci:.4f}), σ = {sig:.1f}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("B. Проверяем ВСЕ пары (W[w][b] → H[0][31])")
    print("=" * 70)

    # Scan all input bits → H[0][31]
    N = 20000
    biases = []
    anomalies = []

    for w in range(16):
        for b in range(32):
            match = 0
            for trial in range(N):
                W = [np.random.randint(0, 2**32) for _ in range(16)]
                H = sha256_hash(W)
                ib = (W[w] >> b) & 1
                ob = (H[0] >> 31) & 1
                if ib == ob: match += 1
            p = match / N
            biases.append((w, b, p))
            sigma = abs(p - 0.5) / np.sqrt(0.25 / N)
            if sigma > 3.0:
                anomalies.append((w, b, p, sigma))

    # Stats
    all_p = [x[2] for x in biases]
    print(f"  512 input bits → H[0][31] (N={N} each):")
    print(f"    Mean P: {np.mean(all_p):.5f}")
    print(f"    Std P:  {np.std(all_p):.5f}")
    print(f"    Expected std: {np.sqrt(0.25/N):.5f}")
    print(f"    Max |P-0.5|: {max(abs(p-0.5) for _,_,p in biases):.5f}")
    print(f"    Anomalies (>3σ): {len(anomalies)}/512")
    print(f"    Expected false positives (>3σ): {512*0.0027:.1f}")

    if anomalies:
        print(f"\n    Anomalous pairs:")
        for w, b, p, sig in sorted(anomalies, key=lambda x: -x[3]):
            print(f"      W[{w}][{b}] → H[0][31]: P={p:.4f}, σ={sig:.1f}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("C. Проверяем W[7][15] → ВСЕ выходные биты")
    print("=" * 70)

    biases_out = []
    for ow in range(8):
        for ob in range(32):
            match = 0
            for trial in range(N):
                W = [np.random.randint(0, 2**32) for _ in range(16)]
                H = sha256_hash(W)
                ib = (W[7] >> 15) & 1
                obit = (H[ow] >> ob) & 1
                if ib == obit: match += 1
            p = match / N
            sigma = abs(p - 0.5) / np.sqrt(0.25 / N)
            biases_out.append((ow, ob, p, sigma))

    all_p_out = [x[2] for x in biases_out]
    anom_out = [x for x in biases_out if x[3] > 3.0]
    print(f"  W[7][15] → all 256 output bits (N={N}):")
    print(f"    Mean P: {np.mean(all_p_out):.5f}")
    print(f"    Anomalies (>3σ): {len(anom_out)}/256")
    print(f"    Expected: {256*0.0027:.1f}")
    if anom_out:
        for ow, ob, p, sig in sorted(anom_out, key=lambda x: -x[3])[:5]:
            print(f"      W[7][15] → H[{ow}][{ob}]: P={p:.4f}, σ={sig:.1f}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("D. GLOBAL scan: max bias across random pairs")
    print("=" * 70)

    # Sample 200 random (in,out) pairs, find max bias
    max_sigma = 0
    max_pair = None
    for _ in range(200):
        iw = np.random.randint(0, 16)
        ib = np.random.randint(0, 32)
        ow = np.random.randint(0, 8)
        ob = np.random.randint(0, 32)

        match = 0
        for trial in range(N):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            H = sha256_hash(W)
            i_bit = (W[iw] >> ib) & 1
            o_bit = (H[ow] >> ob) & 1
            if i_bit == o_bit: match += 1
        p = match / N
        sigma = abs(p - 0.5) / np.sqrt(0.25 / N)
        if sigma > max_sigma:
            max_sigma = sigma
            max_pair = (iw, ib, ow, ob, p)

    print(f"  200 random (in,out) pairs, N={N}:")
    print(f"    Max σ: {max_sigma:.2f}")
    print(f"    Pair: W[{max_pair[0]}][{max_pair[1]}] → H[{max_pair[2]}][{max_pair[3]}]: P={max_pair[4]:.4f}")
    print(f"    Expected max σ from 200 normal: ~3.0")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("ВЕРДИКТ")
    print("=" * 70)


if __name__ == "__main__":
    main()
