"""
ТРИ ДИКИЕ ИДЕИ.

1. STATE INTERFERENCE: 3-4 слова с δ, оптимизированных чтобы
   их вклады в state ГАСИЛИ друг друга (destructive interference).

2. INCREMENTAL BIRTHDAY: сначала k=3 match (325 пар),
   потом ИЗ НИХ ищем k=4, потом k=5... bootstrap к collision.

3. SPECTRAL RESONANCE: подаём белый шум, слушаем отклик.
   FFT(H(M_1), H(M_2), ...) может показать резонансы.
"""

import numpy as np
import struct, hashlib
import time

MASK32 = 0xFFFFFFFF

def hw(x): return bin(x).count('1')

def sha256_full(W16):
    raw = struct.pack('>16I', *W16)
    return struct.unpack('>8I', hashlib.sha256(raw).digest())


def main():
    np.random.seed(42)

    print("=" * 70)
    print("ТРИ ДИКИЕ ИДЕИ")
    print("=" * 70)

    # ═══════════════════
    print(f"\n{'='*70}")
    print("1. STATE INTERFERENCE: multi-word δ optimized for cancellation")
    print(f"{'='*70}")

    # Strategy: try MANY random multi-word δ patterns.
    # For each: measure δH. Keep the BEST pattern.
    # Then: use that pattern for massive search.

    # Phase 1: find best δ PATTERN (which words, which bits)
    print(f"\n  Phase 1: search for best δ pattern (50K patterns × 50 evals)")

    best_pattern = None
    best_avg = 256
    t0 = time.time()

    for _ in range(50000):
        # Random pattern: 1-4 words, 1 bit each
        n_words = np.random.randint(1, 5)
        pattern = {}
        for _ in range(n_words):
            w = np.random.randint(0, 16)
            b = np.random.randint(0, 32)
            pattern[w] = pattern.get(w, 0) ^ (1 << b)
        # Remove zero entries
        pattern = {w: v for w, v in pattern.items() if v != 0}
        if not pattern:
            continue

        # Quick eval: 50 random messages
        total_dH = 0
        for _ in range(50):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            W2 = list(W)
            for w, v in pattern.items():
                W2[w] ^= v
            H1 = sha256_full(W); H2 = sha256_full(W2)
            total_dH += sum(hw(H1[i]^H2[i]) for i in range(8))

        avg = total_dH / 50
        if avg < best_avg:
            best_avg = avg
            best_pattern = dict(pattern)

    elapsed = time.time() - t0
    print(f"  Best pattern found in {elapsed:.1f}s:")
    print(f"    δ words: {best_pattern}")
    total_hw = sum(hw(v) for v in best_pattern.values())
    print(f"    Total HW(δ): {total_hw}")
    print(f"    Avg δH: {best_avg:.1f}")

    # Phase 2: massive search with best pattern
    print(f"\n  Phase 2: 100K search with best pattern")
    best_dH = 256
    for _ in range(100000):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W)
        for w, v in best_pattern.items():
            W2[w] ^= v
        H1 = sha256_full(W); H2 = sha256_full(W2)
        dH = sum(hw(H1[i]^H2[i]) for i in range(8))
        if dH < best_dH:
            best_dH = dH

    print(f"  Best δH with optimized pattern: {best_dH}")

    # Compare: random 1-bit, 100K
    best_random = 256
    for _ in range(100000):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W)
        W2[np.random.randint(0, 16)] ^= (1 << np.random.randint(0, 32))
        H1 = sha256_full(W); H2 = sha256_full(W2)
        dH = sum(hw(H1[i]^H2[i]) for i in range(8))
        if dH < best_random: best_random = dH

    print(f"  Best δH random 1-bit:           {best_random}")
    print(f"  Advantage: {best_random - best_dH:+d}")

    # ═══════════════════
    print(f"\n{'='*70}")
    print("2. INCREMENTAL BIRTHDAY: bootstrap from partial match")
    print(f"{'='*70}")

    # Step 1: collect MANY messages
    # Step 2: find pairs matching on k LSBs
    # Step 3: among those, find ones matching on MORE bits

    N_MSGS = 300000
    print(f"\n  Collecting {N_MSGS} hashes...")
    t0 = time.time()

    messages = []
    hashes = []
    for _ in range(N_MSGS):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H = sha256_full(W)
        messages.append(W)
        hashes.append(H)

    elapsed = time.time() - t0
    print(f"  Collected in {elapsed:.1f}s")

    # Step 2: bucket by k LSBs, find pairs
    for k in [2, 3, 4]:
        mask_k = (1 << k) - 1
        buckets = {}

        for idx in range(N_MSGS):
            key = tuple(hashes[idx][i] & mask_k for i in range(8))
            if key not in buckets:
                buckets[key] = []
            buckets[key].append(idx)

        # Find pairs within buckets
        pairs = []
        for key, indices in buckets.items():
            if len(indices) >= 2:
                for i in range(min(len(indices), 10)):
                    for j in range(i+1, min(len(indices), 10)):
                        idx1, idx2 = indices[i], indices[j]
                        full_dH = sum(hw(hashes[idx1][w] ^ hashes[idx2][w]) for w in range(8))
                        pairs.append(full_dH)

        if pairs:
            print(f"\n  k={k}-bit match: {len(pairs)} pairs")
            print(f"    Full δH: mean={np.mean(pairs):.1f}, min={min(pairs)}, max={max(pairs)}")
            # Expected: each pair has 8*(32-k) random bits = 8*(32-k)/2 expected δH
            expected_dH = 8 * (32 - k) / 2
            print(f"    Expected δH (random remaining bits): {expected_dH:.0f}")
            print(f"    Actual vs expected: {np.mean(pairs) - expected_dH:+.1f}")

            # Among k-bit matches: find pairs with ADDITIONAL bits matching
            extra_match = 0
            for d in pairs:
                if d < expected_dH - 10:  # significantly better than expected
                    extra_match += 1
            print(f"    Pairs with δH < expected-10: {extra_match}/{len(pairs)}")
        else:
            print(f"\n  k={k}: no pairs found in {N_MSGS} messages")

    # Step 3: can we CHAIN? k=2 pair → refine to k=3?
    print(f"\n  CHAINING: k=2 pairs → check for k=3 match among them")
    mask2 = 0x3; mask3 = 0x7

    bucket2 = {}
    for idx in range(N_MSGS):
        key2 = tuple(hashes[idx][i] & mask2 for i in range(8))
        if key2 not in bucket2:
            bucket2[key2] = []
        bucket2[key2].append(idx)

    # Find k=2 pairs
    pairs_2 = []
    for key, indices in bucket2.items():
        if len(indices) >= 2:
            for i in range(min(len(indices), 5)):
                for j in range(i+1, min(len(indices), 5)):
                    pairs_2.append((indices[i], indices[j]))
        if len(pairs_2) > 50000:
            break

    # Among k=2 pairs: how many are also k=3 matches?
    k3_matches = 0
    for idx1, idx2 in pairs_2[:10000]:
        if all((hashes[idx1][i] & mask3) == (hashes[idx2][i] & mask3) for i in range(8)):
            k3_matches += 1

    print(f"  k=2 pairs: {len(pairs_2)}")
    print(f"  Of those, k=3 matches: {k3_matches} ({k3_matches/min(len(pairs_2),10000)*100:.2f}%)")
    print(f"  Expected (independent): {1/256*100:.2f}%")
    print(f"  → {'CORRELATION!' if k3_matches/min(len(pairs_2),10000) > 1.5/256 else 'Independent (no bootstrap advantage)'}")

    # ═══════════════════
    print(f"\n{'='*70}")
    print("3. SPECTRAL RESONANCE: FFT of hash sequence")
    print(f"{'='*70}")

    # Treat SHA-256 as a SYSTEM: input = M, output = H(M).
    # Feed STRUCTURED input (counter: M=0,1,2,...), observe output.
    # FFT of output sequence → frequency response of the system.
    # Peaks in FFT → resonant frequencies → exploitable structure.

    N_SEQ = 10000
    seed = [0]*16

    # Generate sequence: H(0), H(1), H(2), ...
    seq_bit0 = []  # bit 0 of H[0]
    seq_hw = []    # HW of full hash

    for n in range(N_SEQ):
        W = list(seed); W[0] = n
        H = sha256_full(W)
        seq_bit0.append((H[0] >> 0) & 1)
        seq_hw.append(sum(hw(h) for h in H))

    # FFT of bit 0 sequence
    bit0_arr = np.array(seq_bit0, dtype=float) - 0.5  # center
    fft_bit0 = np.abs(np.fft.fft(bit0_arr))**2 / N_SEQ
    # FFT of HW sequence
    hw_arr = np.array(seq_hw, dtype=float) - 128
    fft_hw = np.abs(np.fft.fft(hw_arr))**2 / N_SEQ

    # Find peaks
    top_freqs_bit = np.argsort(fft_bit0[1:N_SEQ//2])[::-1][:5] + 1
    top_freqs_hw = np.argsort(fft_hw[1:N_SEQ//2])[::-1][:5] + 1

    print(f"\n  Sequential input M=0,1,2,...{N_SEQ-1}:")
    print(f"\n  Bit 0 of H[0] — top 5 spectral peaks:")
    for f in top_freqs_bit:
        power = fft_bit0[f]
        # Expected for white noise: ~1/N
        ratio = power * N_SEQ
        print(f"    freq={f:>5} (period={N_SEQ/f:.1f}): power={power:.4f}, "
              f"ratio vs noise={ratio:.1f}× {'★ PEAK' if ratio > 5 else ''}")

    print(f"\n  HW(hash) — top 5 spectral peaks:")
    for f in top_freqs_hw:
        power = fft_hw[f]
        ratio = power * N_SEQ / np.mean(fft_hw[1:N_SEQ//2])
        print(f"    freq={f:>5} (period={N_SEQ/f:.1f}): power={power:.2f}, "
              f"ratio={ratio:.1f}× {'★ PEAK' if ratio > 5 else ''}")

    # Flatness of spectrum (1.0 = perfect white noise)
    flatness_bit = np.std(fft_bit0[1:N_SEQ//2]) / np.mean(fft_bit0[1:N_SEQ//2])
    flatness_hw = np.std(fft_hw[1:N_SEQ//2]) / np.mean(fft_hw[1:N_SEQ//2])

    print(f"\n  Spectral flatness (1.0 = white noise):")
    print(f"    Bit 0: {flatness_bit:.3f}")
    print(f"    HW:    {flatness_hw:.3f}")
    print(f"    → {'FLAT (no resonance)' if flatness_bit < 1.2 and flatness_hw < 1.2 else 'STRUCTURE DETECTED ★'}")

    # Compare with TRUE white noise
    noise = np.random.randn(N_SEQ)
    fft_noise = np.abs(np.fft.fft(noise))**2 / N_SEQ
    flatness_noise = np.std(fft_noise[1:N_SEQ//2]) / np.mean(fft_noise[1:N_SEQ//2])
    print(f"    White noise reference: {flatness_noise:.3f}")

    # ═══════════════════
    print(f"\n{'='*70}")
    print("ИТОГ")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
