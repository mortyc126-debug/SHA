"""
АНОМАЛИЯ 3: Schedule autocorrelation lag=1 = 0.228

Факт: HW(W[r]) коррелирует с HW(W[r+1]) в expanded schedule.
W[r] = σ1(W[r-2]) + W[r-7] + σ0(W[r-15]) + W[r-16]

Вопросы:
  A. Это свойство schedule формулы? (W[r] и W[r+1] делят W[r-1], W[r-6], W[r-14])
  B. Или свойство ADD (carry)?
  C. Специфично для SHA-256 или для любой рекуррентной формулы?
"""

import numpy as np
import struct, hashlib

MASK32 = 0xFFFFFFFF

def hw(x): return bin(x).count('1')
def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)


def sha256_schedule(W16):
    W = list(W16)
    for r in range(16, 64):
        W.append((sigma1(W[r-2]) + W[r-7] + sigma0(W[r-15]) + W[r-16]) & MASK32)
    return W


def main():
    np.random.seed(42)

    print("=" * 70)
    print("АНОМАЛИЯ 3: Schedule autocorrelation = 0.228")
    print("=" * 70)

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("A. Воспроизводим")
    print("=" * 70)

    autocorrs = []
    for _ in range(2000):
        W16 = [np.random.randint(0, 2**32) for _ in range(16)]
        W = sha256_schedule(W16)

        hws = [hw(w) for w in W[16:]]  # only expanded part
        seq = np.array(hws, dtype=float)
        seq -= np.mean(seq)
        if np.std(seq) > 0:
            ac = np.sum(seq[:-1] * seq[1:]) / np.sum(seq**2)
            autocorrs.append(ac)

    print(f"  SHA-256 schedule autocorrelation (expanded W[16..63]):")
    print(f"    Mean: {np.mean(autocorrs):.4f}")
    print(f"    Std:  {np.std(autocorrs):.4f}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("B. Shared inputs — W[r] и W[r+1] делят операнды")
    print("=" * 70)

    # W[r]   = σ1(W[r-2])  + W[r-7]  + σ0(W[r-15]) + W[r-16]
    # W[r+1] = σ1(W[r-1])  + W[r-6]  + σ0(W[r-14]) + W[r-15]
    #
    # Shared: NONE directly! All indices shift by 1.
    # But: W[r-1] was input to W[r] computation?
    # No — W[r] uses W[r-2], not W[r-1].
    #
    # Indirect: W[r-7] and W[r-6] are ADJACENT → correlated if THEY are correlated.
    # This is recursive!

    print(f"  W[r]   uses: W[r-2], W[r-7],  W[r-15], W[r-16]")
    print(f"  W[r+1] uses: W[r-1], W[r-6],  W[r-14], W[r-15]")
    print(f"  Direct overlap: W[r-15] appears in W[r+1] as σ0 input")
    print(f"                  W[r-15] appears in W[r] as σ0 input")
    print(f"  Wait — W[r] uses W[r-15], W[r+1] uses W[r-14]. NO direct overlap!")
    print(f"  But W[r-7] and W[r-6] are adjacent → if THEY correlate, it propagates.")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("C. Control: simple LFSR (XOR recurrence, no ADD)")
    print("=" * 70)

    # Replace ADD with XOR: W[r] = σ1(W[r-2]) ⊕ W[r-7] ⊕ σ0(W[r-15]) ⊕ W[r-16]
    autocorrs_xor = []
    for _ in range(2000):
        W16 = [np.random.randint(0, 2**32) for _ in range(16)]
        W = list(W16)
        for r in range(16, 64):
            W.append(sigma1(W[r-2]) ^ W[r-7] ^ sigma0(W[r-15]) ^ W[r-16])

        hws = [hw(w) for w in W[16:]]
        seq = np.array(hws, dtype=float)
        seq -= np.mean(seq)
        if np.std(seq) > 0:
            ac = np.sum(seq[:-1] * seq[1:]) / np.sum(seq**2)
            autocorrs_xor.append(ac)

    print(f"  XOR schedule (no carry): autocorr = {np.mean(autocorrs_xor):.4f}")
    print(f"  ADD schedule (SHA-256):  autocorr = {np.mean(autocorrs):.4f}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("D. Control: generic ADD recurrence")
    print("=" * 70)

    # W[r] = W[r-1] + W[r-2] (Fibonacci-like)
    autocorrs_fib = []
    for _ in range(2000):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        for r in range(16, 64):
            W.append((W[r-1] + W[r-2]) & MASK32)

        hws = [hw(w) for w in W[16:]]
        seq = np.array(hws, dtype=float)
        seq -= np.mean(seq)
        if np.std(seq) > 0:
            ac = np.sum(seq[:-1] * seq[1:]) / np.sum(seq**2)
            autocorrs_fib.append(ac)

    print(f"  Fibonacci ADD: autocorr = {np.mean(autocorrs_fib):.4f}")

    # W[r] = W[r-2] + W[r-7] + W[r-15] + W[r-16] (same taps, no sigma)
    autocorrs_plain = []
    for _ in range(2000):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        for r in range(16, 64):
            W.append((W[r-2] + W[r-7] + W[r-15] + W[r-16]) & MASK32)

        hws = [hw(w) for w in W[16:]]
        seq = np.array(hws, dtype=float)
        seq -= np.mean(seq)
        if np.std(seq) > 0:
            ac = np.sum(seq[:-1] * seq[1:]) / np.sum(seq**2)
            autocorrs_plain.append(ac)

    print(f"  Plain ADD (same taps, no σ): autocorr = {np.mean(autocorrs_plain):.4f}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("E. Deeper lags")
    print("=" * 70)

    # Autocorrelation at various lags
    for lag in [1, 2, 3, 5, 7, 10, 15, 16]:
        acs = []
        for _ in range(2000):
            W16 = [np.random.randint(0, 2**32) for _ in range(16)]
            W = sha256_schedule(W16)
            hws = np.array([hw(w) for w in W[16:]], dtype=float)
            hws -= np.mean(hws)
            if np.std(hws) > 0 and len(hws) > lag:
                ac = np.sum(hws[:-lag] * hws[lag:]) / np.sum(hws**2)
                acs.append(ac)
        print(f"    Lag {lag:>2}: autocorr = {np.mean(acs):.4f}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("ВЕРДИКТ")
    print("=" * 70)


if __name__ == "__main__":
    main()
