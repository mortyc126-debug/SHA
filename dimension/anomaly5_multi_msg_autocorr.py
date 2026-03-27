"""
АНОМАЛИЯ 5: Multi-message autocorrelation = 0.0995

Факт: HW(H(M₁)), HW(H(M₂)), HW(H(M₃))... для ПОСЛЕДОВАТЕЛЬНЫХ
сообщений (M₁, M₂=M₁+1, ...) показывает автокорреляцию ~0.10.

Вопросы:
  A. Воспроизводится ли?
  B. Что значит "последовательные"? (increment W[0]? random?)
  C. Это свойство SHA-256 или артефакт increment?
"""

import numpy as np
import struct, hashlib

MASK32 = 0xFFFFFFFF

def sha256_hash(W16):
    raw = struct.pack('>16I', *W16)
    return struct.unpack('>8I', hashlib.sha256(raw).digest())

def hw(x): return bin(x).count('1')


def main():
    np.random.seed(42)

    print("=" * 70)
    print("АНОМАЛИЯ 5: Multi-message autocorrelation = 0.0995")
    print("=" * 70)

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("A. Sequential messages (W[0] += 1)")
    print("=" * 70)

    W_base = [np.random.randint(0, 2**32) for _ in range(16)]
    hws_seq = []
    for i in range(10000):
        W = list(W_base)
        W[0] = (W[0] + i) & MASK32
        H = sha256_hash(W)
        hws_seq.append(sum(hw(h) for h in H))

    seq = np.array(hws_seq, dtype=float)
    seq -= np.mean(seq)
    ac1 = np.sum(seq[:-1] * seq[1:]) / np.sum(seq**2)

    print(f"  Sequential W[0]+=1 (10K messages):")
    print(f"    Mean HW(H): {np.mean(hws_seq):.2f}")
    print(f"    Autocorr lag=1: {ac1:.4f}")

    # Lags
    for lag in [1, 2, 3, 5, 10, 100]:
        ac = np.sum(seq[:-lag] * seq[lag:]) / np.sum(seq**2)
        print(f"    Lag {lag:>3}: {ac:.4f}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("B. Random messages (no sequential relation)")
    print("=" * 70)

    hws_rand = []
    for _ in range(10000):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H = sha256_hash(W)
        hws_rand.append(sum(hw(h) for h in H))

    seq_r = np.array(hws_rand, dtype=float)
    seq_r -= np.mean(seq_r)
    ac1_r = np.sum(seq_r[:-1] * seq_r[1:]) / np.sum(seq_r**2)

    print(f"  Random messages (10K):")
    print(f"    Autocorr lag=1: {ac1_r:.4f}")
    print(f"  → {'DIFFERENT' if abs(ac1 - ac1_r) > 0.02 else 'SAME'}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("C. Bit-level: sequential W[0]+=1 changes only LOW bits")
    print("=" * 70)

    # W[0]+=1: if low bits are 1111, carry flips many bits
    # Most of the time: only bit 0 changes → HW(δW) = 1 or 2
    # This means: H(W) and H(W+1) differ by ~128 bits (full avalanche)
    # But the INPUT is "almost the same" → is there carry bias?

    # Count: how many bits of W actually change for +1?
    dW_hws = []
    for i in range(10000):
        w0 = (W_base[0] + i) & MASK32
        w1 = (w0 + 1) & MASK32
        dW_hws.append(hw(w0 ^ w1))

    print(f"  HW(W[0] ⊕ (W[0]+1)) distribution:")
    print(f"    Mean: {np.mean(dW_hws):.2f}")
    from collections import Counter
    cnt = Counter(dW_hws)
    for k in sorted(cnt.keys()):
        print(f"      HW={k}: {cnt[k]/100:.1f}%")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("D. Control: MD5 sequential")
    print("=" * 70)

    hws_md5 = []
    for i in range(10000):
        W = list(W_base)
        W[0] = (W[0] + i) & MASK32
        raw = struct.pack('>16I', *W)
        H_md5 = hashlib.md5(raw).digest()
        hws_md5.append(sum(hw(b) for b in H_md5))

    seq_md5 = np.array(hws_md5, dtype=float)
    seq_md5 -= np.mean(seq_md5)
    ac1_md5 = np.sum(seq_md5[:-1] * seq_md5[1:]) / np.sum(seq_md5**2)

    print(f"  MD5 sequential W[0]+=1:")
    print(f"    Autocorr lag=1: {ac1_md5:.4f}")

    # SHA-512
    hws_512 = []
    for i in range(10000):
        W = list(W_base)
        W[0] = (W[0] + i) & MASK32
        raw = struct.pack('>16I', *W)
        h512 = hashlib.sha512(raw).digest()
        hws_512.append(sum(hw(b) for b in h512))

    seq_512 = np.array(hws_512, dtype=float)
    seq_512 -= np.mean(seq_512)
    ac1_512 = np.sum(seq_512[:-1] * seq_512[1:]) / np.sum(seq_512**2)

    print(f"  SHA-512 sequential:")
    print(f"    Autocorr lag=1: {ac1_512:.4f}")

    print(f"\n  Comparison:")
    print(f"    SHA-256 sequential: {ac1:.4f}")
    print(f"    SHA-256 random:     {ac1_r:.4f}")
    print(f"    MD5 sequential:     {ac1_md5:.4f}")
    print(f"    SHA-512 sequential: {ac1_512:.4f}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("E. Bigger sample for significance")
    print("=" * 70)

    # 100K sequential
    hws_big = []
    for i in range(100000):
        W = list(W_base)
        W[0] = (W[0] + i) & MASK32
        H = sha256_hash(W)
        hws_big.append(sum(hw(h) for h in H))

    seq_big = np.array(hws_big, dtype=float)
    seq_big -= np.mean(seq_big)
    ac1_big = np.sum(seq_big[:-1] * seq_big[1:]) / np.sum(seq_big**2)
    # SE of autocorrelation ≈ 1/sqrt(N)
    se = 1.0 / np.sqrt(100000)

    print(f"  100K sequential:")
    print(f"    Autocorr lag=1: {ac1_big:.5f}")
    print(f"    SE: {se:.5f}")
    print(f"    σ: {abs(ac1_big)/se:.1f}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("ВЕРДИКТ")
    print("=" * 70)


if __name__ == "__main__":
    main()
