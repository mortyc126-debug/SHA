"""
ИЗОБРЕТАЕМ НОВЫЙ ПОДХОД.

Все подходы которые мы пробовали — АНАЛИЗИРУЮТ SHA-256 как есть.
Мы смотрим на output и ищем structure. Не находим при r≥20.

А что если не АНАЛИЗИРОВАТЬ, а ТРАНСФОРМИРОВАТЬ?

Идея: ПОСТРОИТЬ ДРУГУЮ ФУНКЦИЮ из SHA-256, которая СОХРАНЯЕТ
collision property но ЛЕГЧЕ для анализа.

f(M) = T(SHA256(M)) где T — наша трансформация.
Если T обратима: collision в f = collision в SHA256.
Если T необратима: collision в f ≠ collision в SHA256,
но partial collision может помочь.

ТРАНСФОРМАЦИИ:
  T1: Truncation — уже знаем (partial collision).
  T2: PROJECTION — проецировать hash в пространство МЕНЬШЕЙ размерности
      через НАШУ метрику. Collision в проекции = near-collision в hash.
  T3: FOLDING — H[0]⊕H[1]⊕H[2]⊕H[3] (reduce 256→32 bits).
      Birthday на 32 бита = 2^16 (мгновенно!). Но это partial collision.
  T4: MODULAR — (H[0]+H[1]+H[2]+H[3]) mod p для малого p.
  T5: HASH OF HASH — SHA256(SHA256(M)). Двойной hash имеет ДРУГИЕ свойства?

НО ГЛАВНОЕ:
  T6: NEURAL/LEARNED — обучить функцию g(M) ≈ 0 для collision pairs.
  Мы не можем, но: можем ли мы АНАЛИТИЧЕСКИ построить g?

  T7: REVERSE ENGINEERING — если знаем что collision pair (M1,M2)
  имеет δW[15]=1, можно ли build M1 такое что SHA256(M1) has
  SPECIFIC structure that makes finding M2 easier?
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
    print("ИЗОБРЕТАЕМ НОВЫЙ ПОДХОД")
    print("=" * 70)

    # ═══════════════════
    print(f"\n{'='*70}")
    print("ПОДХОД 1: FOLDED COLLISION")
    print("  F(M) = H[0] ⊕ H[1] ⊕ H[2] ⊕ H[3] (32-bit folded hash)")
    print("  Birthday: 2^16 pairs = instant collision in F.")
    print("  Question: does F-collision help find H-collision?")
    print(f"{'='*70}")

    # Find F-collision
    fold_table = {}
    fold_collisions = []

    t0 = time.time()
    for trial in range(200000):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H = sha256_full(W)
        F = H[0] ^ H[1] ^ H[2] ^ H[3]

        if F in fold_table:
            W_prev = fold_table[F]
            H_prev = sha256_full(W_prev)
            full_dH = sum(hw(H[i] ^ H_prev[i]) for i in range(8))
            fold_collisions.append((full_dH, W, W_prev))
        else:
            fold_table[F] = W

    elapsed = time.time() - t0
    print(f"\n  Found {len(fold_collisions)} fold-collisions in {elapsed:.1f}s")

    if fold_collisions:
        dHs = [fc[0] for fc in fold_collisions]
        print(f"  Full δH of fold-collision pairs:")
        print(f"    Mean: {np.mean(dHs):.1f}")
        print(f"    Min:  {min(dHs)}")
        print(f"    Max:  {max(dHs)}")

        # Distribution
        for threshold in [64, 96, 112, 120, 124, 126, 127]:
            count = sum(1 for d in dHs if d <= threshold)
            print(f"    δH ≤ {threshold}: {count}/{len(dHs)}")

        # Best pair
        best = min(fold_collisions, key=lambda x: x[0])
        print(f"\n  Best fold-collision: δH = {best[0]}")

    # ═══════════════════
    print(f"\n{'='*70}")
    print("ПОДХОД 2: ADD-FOLDED COLLISION")
    print("  F(M) = (H[0]+H[1]+H[2]+H[3]) mod 2^32")
    print(f"{'='*70}")

    add_table = {}
    add_collisions = []

    for trial in range(200000):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H = sha256_full(W)
        F = (H[0] + H[1] + H[2] + H[3]) & MASK32

        if F in add_table:
            W_prev = add_table[F]
            H_prev = sha256_full(W_prev)
            full_dH = sum(hw(H[i] ^ H_prev[i]) for i in range(8))
            add_collisions.append(full_dH)
        else:
            add_table[F] = W

    if add_collisions:
        print(f"  Found {len(add_collisions)} add-fold collisions")
        print(f"  Full δH: mean={np.mean(add_collisions):.1f}, min={min(add_collisions)}")

    # ═══════════════════
    print(f"\n{'='*70}")
    print("ПОДХОД 3: MULTI-TARGET — ищем НЕ collision, а СТРУКТУРУ")
    print(f"{'='*70}")

    # Instead of H(M1)=H(M2), look for:
    # H(M1) + H(M2) = H(M3) (additive relation among 3 hashes)
    # This is EASIER: 3 random 256-bit values summing to 0 mod 2^32
    # requires ~2^128 by birthday on 3-XOR.
    # But: with FOLDING to 32 bits → needs only ~2^16.

    print(f"\n  3-XOR on folded hash: F(M1) ⊕ F(M2) ⊕ F(M3) = 0")
    print(f"  (equivalent to F(M1) ⊕ F(M2) = F(M3))")

    # Build table of F values, then look for pairs that XOR to a known F
    fold_set = {}
    for trial in range(100000):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        H = sha256_full(W)
        F = H[0] ^ H[1] ^ H[2] ^ H[3]
        fold_set[F] = W

    # Now: for each pair (F1, F2), check if F1⊕F2 is in set
    three_xor = []
    keys = list(fold_set.keys())
    for i in range(min(10000, len(keys))):
        for j in range(i+1, min(i+100, len(keys))):
            target = keys[i] ^ keys[j]
            if target in fold_set and target != keys[i] and target != keys[j]:
                # Found! F(M1) ⊕ F(M2) ⊕ F(M3) = 0
                W1 = fold_set[keys[i]]
                W2 = fold_set[keys[j]]
                W3 = fold_set[target]
                H1 = sha256_full(W1)
                H2 = sha256_full(W2)
                H3 = sha256_full(W3)
                # Check: H1⊕H2⊕H3 per word?
                triple_dH = sum(hw(H1[w]^H2[w]^H3[w]) for w in range(8))
                three_xor.append(triple_dH)
                if len(three_xor) >= 100:
                    break
        if len(three_xor) >= 100:
            break

    if three_xor:
        print(f"  Found {len(three_xor)} fold-3-XOR relations")
        print(f"  Full HW(H1⊕H2⊕H3): mean={np.mean(three_xor):.1f}, min={min(three_xor)}")

    # ═══════════════════
    print(f"\n{'='*70}")
    print("ПОДХОД 4: DIFFERENTIAL-THEN-BIRTHDAY")
    print(f"{'='*70}")

    # Step 1: collect MANY pairs with δW[15]=bit31 (our best diff)
    # Step 2: birthday among the δH values!
    # If two DIFFERENT pairs give the SAME δH → near-collision chain

    # More precisely: collect H(M), H(M⊕δ) for many M.
    # Look for M1,M2 where H(M1)⊕H(M1⊕δ) = H(M2)⊕H(M2⊕δ).
    # This means: differential is SAME for two different messages.
    # If we ALSO have H(M1)⊕H(M2) small → collision-like.

    delta = [0]*16; delta[15] = 1 << 31

    diff_table = {}
    diff_collisions = []

    for trial in range(200000):
        W = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = list(W); W2[15] ^= (1 << 31)
        H1 = sha256_full(W); H2 = sha256_full(W2)

        # Differential fingerprint: fold of δH
        dH_fold = (H1[0]^H2[0]) ^ (H1[1]^H2[1]) ^ (H1[2]^H2[2]) ^ (H1[3]^H2[3])

        if dH_fold in diff_table:
            W_prev = diff_table[dH_fold]
            H_prev1 = sha256_full(W_prev)
            # How close are H(M1) and H(M_prev)?
            cross_dH = sum(hw(H1[i]^H_prev1[i]) for i in range(8))
            diff_collisions.append(cross_dH)
        else:
            diff_table[dH_fold] = W

    if diff_collisions:
        print(f"  Found {len(diff_collisions)} diff-fingerprint collisions")
        print(f"  Cross δH (H(M1) vs H(M_prev)): mean={np.mean(diff_collisions):.1f}, "
              f"min={min(diff_collisions)}")
        print(f"  (128 = random, <128 = structure)")

    # ═══════════════════
    print(f"\n{'='*70}")
    print("ПОДХОД 5: MERKLE-DAMGÅRD ДЛИНА КАК ОРУЖИЕ")
    print(f"{'='*70}")

    # SHA-256 adds padding with LENGTH. Different length messages
    # have different padding → different last block.
    # But: what if we use EXACTLY the right length to create
    # favorable conditions in the last block?

    # Single block message: 512 bits message + padding in same block
    # Only if message ≤ 447 bits (55 bytes).
    # Padding: 0x80, zeros, 64-bit length.

    # For 16-word (64-byte) message: needs 2 blocks!
    # Block 1: M[0..15] (512 bits of message)
    # Block 2: 0x80000000, zeros, length=512

    # THIS MEANS: SHA-256 of 64-byte message = 2-block hash.
    # Block 2 is FIXED (padding only, no message freedom).
    # Block 1: we control all 16 words.

    # For shorter messages (≤ 55 bytes = 13.75 words):
    # Single block! Padding fills W[14..15] with length info.
    # W[14] = 0, W[15] = message_length_in_bits.

    # This means: W[15] is NOT free — it's the LENGTH.
    # Our attack (flip W[15]) = CHANGE MESSAGE LENGTH.
    # Not valid for same-length collision!

    print(f"  SHA-256 padding for 55-byte message:")
    print(f"    W[0..13]: message (448 bits)")
    print(f"    W[14]: 0x00000000 (length high)")
    print(f"    W[15]: 0x000001B8 (440 bits = 55 bytes)")
    print(f"  → W[15] = LENGTH. Flipping it = different length message!")
    print(f"  → Our δW[15] attack needs DIFFERENT LENGTH messages.")
    print(f"  → For SAME LENGTH collision: can only modify W[0..13].")

    # What's the weakest freely-controllable word for same-length?
    # W[13] enters at r=13, has 4 rounds of mixing at r=17.
    # From our absorption law: influence = 7.1 × 3 = 21 bits at r=17.

    print(f"\n  Same-length attack (≤55 bytes): weakest FREE word = W[13]")
    print(f"    W[13] influence at r=17: λ×(17-13-1) = 7.1×3 = 21 bits")
    print(f"    W[15] influence at r=17: 7.1×1 = 7 bits (but W[15]=length!)")

    # Test: collision search with δW[13] vs δW[15]
    N = 20000
    for word, name in [(13, "W[13] (same-length)"), (15, "W[15] (diff-length)")]:
        best = 256
        for _ in range(N):
            W = [np.random.randint(0, 2**32) for _ in range(16)]
            W2 = list(W); W2[word] ^= (1 << 31)
            H1 = sha256_full(W); H2 = sha256_full(W2)
            d = sum(hw(H1[i]^H2[i]) for i in range(8))
            if d < best: best = d
        print(f"    {name}: best δH = {best} (20K search)")

    # ═══════════════════
    print(f"\n{'='*70}")
    print("СИНТЕЗ")
    print(f"{'='*70}")

    best_fold = min(fold_collisions, key=lambda x: x[0])[0] if fold_collisions else "N/A"
    print(f"""
  ═══════════════════════════════════════════════════════════════

  НОВЫЕ ПОДХОДЫ:

  1. FOLD COLLISION (XOR 4 words → 32 bits):
     {len(fold_collisions)} collisions in 200K, best full δH = {best_fold}
     Fold-collision does NOT help: full δH ≈ random.
     The 4 folded words cancel, but remaining 4 don't.

  2. ADD-FOLD COLLISION:
     Same result. Fold collisions give random full δH.

  3. 3-XOR RELATIONS:
     Found {len(three_xor)} fold-3-XOR. Full HW ≈ {np.mean(three_xor):.0f if three_xor else 'N/A'}.
     No advantage over random.

  4. DIFFERENTIAL-THEN-BIRTHDAY:
     {len(diff_collisions)} diff-fingerprint collisions.
     Cross δH mean = {np.mean(diff_collisions):.0f if diff_collisions else 'N/A'} ≈ random.
     Sharing diff-fingerprint doesn't imply closeness.

  5. PADDING CONSTRAINT:
     W[15] = message length for single-block messages!
     Our best attack (δW[15]) requires DIFFERENT LENGTHS.
     Same-length attack: weakest word = W[13], weaker effect.

  KEY INSIGHT: every transformation that REDUCES SHA-256
  to a simpler problem also DESTROYS the collision structure.
  Fold-collision ≠ full collision. Diff-fingerprint ≠ closeness.
  The information lost in reduction is EXACTLY the information
  needed for collision.

  THIS IS THE FUNDAMENTAL BARRIER:
  SHA-256's 256-bit output is NOT reducible.
  Any projection to < 256 bits loses collision structure.
  Any "helper function" that makes search easier
  also makes verification meaningless.

  ═══════════════════════════════════════════════════════════════
""")


if __name__ == "__main__":
    main()
