"""CE tail: FAST version. Precompute basis, then XOR."""

import numpy as np
import struct, hashlib

def sha256_full(W16):
    return struct.unpack('>8I', hashlib.sha256(struct.pack('>16I', *W16)).digest())

def hw(x): return bin(x).count('1')
MASK32 = 0xFFFFFFFF

def main():
    np.random.seed(42)
    print("=" * 70)
    print("CE TAIL (fast): precomputed basis")
    print("=" * 70)

    W_base = [np.random.randint(0, 2**32) for _ in range(16)]
    H_base = sha256_full(W_base)

    # Build T, kernel, basis CE
    T_rows = []
    for word in range(16):
        for bit in range(32):
            W_mod = list(W_base); W_mod[word] ^= (1 << bit)
            H_mod = sha256_full(W_mod)
            row = []
            for w in range(8):
                d = H_base[w] ^ H_mod[w]
                for b in range(32): row.append((d >> b) & 1)
            T_rows.append(row)

    T = np.array(T_rows, dtype=np.uint8)
    M = T.T.copy()
    pivots = []; row = 0
    for col in range(512):
        found = False
        for r in range(row, 256):
            if M[r,col]==1: M[[row,r]]=M[[r,row]]; found=True; break
        if not found: continue
        pivots.append(col)
        for r in range(256):
            if r!=row and M[r,col]==1: M[r]^=M[row]
        row += 1
    free_vars = [c for c in range(512) if c not in pivots]

    # PRECOMPUTE all 256 basis dW and their H
    print(f"  Precomputing {len(free_vars)} basis vectors + SHA-256...")
    basis_dW = []
    basis_H = []
    for fc in free_vars[:256]:
        x = np.zeros(512, dtype=np.uint8); x[fc] = 1
        for i in range(len(pivots)-1, -1, -1):
            pc = pivots[i]; val = np.uint8(0)
            for j in range(512):
                if j!=pc: val ^= (M[i,j]&x[j])
            x[pc] = val
        dW = [0]*16
        for word in range(16):
            for bit in range(32):
                if x[word*32+bit]: dW[word]^=(1<<bit)
        basis_dW.append(dW)
        W2 = [W_base[i]^dW[i] for i in range(16)]
        basis_H.append(sha256_full(W2))

    # Basis CE HW
    basis_ce_hw = [sum(hw(H_base[i]^basis_H[k][i]) for i in range(8)) for k in range(len(basis_dW))]
    print(f"  Basis CE HW: mean={np.mean(basis_ce_hw):.1f}, min={min(basis_ce_hw)}, max={max(basis_ce_hw)}")

    # Random XOR combinations (FAST: just XOR dW, then 1 SHA per combo)
    print(f"\n  Sampling 10000 random kernel vectors (XOR of basis)...")
    combo_hws = []
    for trial in range(10000):
        n_xor = np.random.randint(2, min(20, len(basis_dW)))
        indices = np.random.choice(len(basis_dW), n_xor, replace=False)
        dW_combo = [0]*16
        for idx in indices:
            for w in range(16): dW_combo[w] ^= basis_dW[idx][w]
        if all(w==0 for w in dW_combo): continue
        W2 = [W_base[i]^dW_combo[i] for i in range(16)]
        H2 = sha256_full(W2)
        ce_hw = sum(hw(H_base[i]^H2[i]) for i in range(8))
        combo_hws.append(ce_hw)

    print(f"\n  CE HW distribution (10K random kernel combos):")
    print(f"    Mean: {np.mean(combo_hws):.1f}")
    print(f"    Std:  {np.std(combo_hws):.1f}")
    print(f"    Min:  {min(combo_hws)}")
    print(f"    Max:  {max(combo_hws)}")

    # Tail
    for t in [80, 90, 95, 100, 105, 110]:
        c = sum(1 for h in combo_hws if h <= t)
        print(f"    HW ≤ {t}: {c}/10000 = {c/100:.2f}%")

    # Random baseline
    rand_hws = []
    for _ in range(10000):
        W2 = [np.random.randint(0, 2**32) for _ in range(16)]
        H2 = sha256_full(W2)
        rand_hws.append(sum(hw(H_base[i]^H2[i]) for i in range(8)))

    print(f"\n  Random δW baseline (10K):")
    print(f"    Mean: {np.mean(rand_hws):.1f}, Min: {min(rand_hws)}")

    print(f"\n  COMPARISON:")
    print(f"    {'Source':<25} {'Mean':>6} {'Min':>5} {'P(≤100)':>8}")
    for name, hws in [("Random δW", rand_hws), ("GF2-kernel (basis)", basis_ce_hw), ("GF2-kernel (combos)", combo_hws)]:
        p = sum(1 for h in hws if h <= 100)/len(hws)*100
        print(f"    {name:<25} {np.mean(hws):>5.1f} {min(hws):>4} {p:>7.2f}%")

    # Is there a difference?
    gain = np.mean(rand_hws) - np.mean(combo_hws)
    tail_gain = min(rand_hws) - min(combo_hws)
    print(f"\n  Mean gain: {gain:+.1f} bits")
    print(f"  Tail gain: {tail_gain:+d} bits")
    print(f"  → {'GF2-kernel has BETTER TAIL!' if tail_gain > 5 else 'GF2-kernel ≈ random tail.'}")


if __name__ == "__main__":
    main()
