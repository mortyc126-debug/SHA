"""
ХВОСТ CE: распределение carry error для GF2-kernel vectors.

rank(CE) = 256 = СРЕДНИЙ carry error full rank.
Но: РАСПРЕДЕЛЕНИЕ HW carry error по kernel vectors может иметь ХВОСТ.
Если существуют vectors с HW(CE) << 128: реальные атаки ИСПОЛЬЗУЮТ ИХ.

Это мост между:
  Нашей теорией (average-case: rank=256)
  и реальными атаками (worst-case: specific low-CE trails)
"""

import numpy as np
import struct, hashlib

MASK32 = 0xFFFFFFFF
K_const = [
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
def add32(x, y): return (x + y) & MASK32
def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)

def sha256_full(W16):
    return struct.unpack('>8I', hashlib.sha256(struct.pack('>16I', *W16)).digest())

def sha256_r(W16, n_r):
    W = list(W16)
    for r in range(16, max(n_r, 16)):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))
    a,b,c,d,e,f,g,h = IV
    for r in range(n_r):
        T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e,f,g)), K_const[r]), W[r])
        T2 = add32(Sigma0(a), Maj(a,b,c))
        h,g,f,e = g,f,e,add32(d,T1)
        d,c,b,a = c,b,a,add32(T1,T2)
    return tuple(add32(IV[i], [a,b,c,d,e,f,g,h][i]) for i in range(8))


def build_kernel_and_ce(hash_fn, W_base):
    H_base = hash_fn(W_base)
    T_rows = []
    for word in range(16):
        for bit in range(32):
            W_mod = list(W_base)
            W_mod[word] ^= (1 << bit)
            H_mod = hash_fn(W_mod)
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
    return M, pivots, free_vars, H_base, W_base, hash_fn


def sample_ce_distribution(M, pivots, free_vars, H_base, W_base, hash_fn, N_samples=5000):
    """Sample carry error HW from random GF2-kernel vectors."""
    ce_hws = []

    for trial in range(N_samples):
        # Random kernel vector: XOR random subset of basis
        n_basis = np.random.randint(1, min(30, len(free_vars)))
        indices = np.random.choice(len(free_vars), n_basis, replace=False)

        x = np.zeros(512, dtype=np.uint8)
        for idx in indices:
            fc = free_vars[idx]
            v = np.zeros(512, dtype=np.uint8); v[fc] = 1
            for i in range(len(pivots)-1, -1, -1):
                pc = pivots[i]; val = np.uint8(0)
                for j in range(512):
                    if j != pc: val ^= (M[i,j] & v[j])
                v[pc] = val
            x = x ^ v

        dW = [0]*16
        for word in range(16):
            for bit in range(32):
                if x[word*32+bit]: dW[word] ^= (1<<bit)

        if all(w == 0 for w in dW): continue

        W2 = [W_base[i]^dW[i] for i in range(16)]
        H2 = hash_fn(W2)
        ce_hw = sum(hw(H_base[i]^H2[i]) for i in range(8))
        ce_hws.append(ce_hw)

    return ce_hws


def main():
    np.random.seed(42)

    print("=" * 70)
    print("ХВОСТ CE: распределение carry error")
    print("=" * 70)

    # =================================================================
    print(f"\n{'=' * 70}")
    print("1. CE DISTRIBUTION для full 64-round SHA-256")
    print("=" * 70)

    W_base = [np.random.randint(0, 2**32) for _ in range(16)]
    M, piv, fv, Hb, Wb, fn = build_kernel_and_ce(sha256_full, W_base)

    print(f"  Kernel dim: {len(fv)}")
    print(f"  Sampling 5000 random kernel vectors...")

    hws_64 = sample_ce_distribution(M, piv, fv, Hb, Wb, fn, 5000)

    print(f"\n  CE HW distribution (64-round):")
    print(f"    Mean: {np.mean(hws_64):.1f}")
    print(f"    Std:  {np.std(hws_64):.1f}")
    print(f"    Min:  {min(hws_64)}")
    print(f"    Max:  {max(hws_64)}")

    # Tail analysis
    for threshold in [90, 100, 110, 120]:
        count = sum(1 for h in hws_64 if h <= threshold)
        print(f"    HW ≤ {threshold}: {count}/5000 = {count/50:.2f}%")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("2. CE DISTRIBUTION для 24-round (наш OptimalHash)")
    print("=" * 70)

    fn24 = lambda W: sha256_r(W, 24)
    M24, piv24, fv24, Hb24, Wb24, _ = build_kernel_and_ce(fn24, W_base)
    print(f"  Kernel dim: {len(fv24)}")

    hws_24 = sample_ce_distribution(M24, piv24, fv24, Hb24, Wb24, fn24, 5000)

    print(f"\n  CE HW distribution (24-round):")
    print(f"    Mean: {np.mean(hws_24):.1f}")
    print(f"    Std:  {np.std(hws_24):.1f}")
    print(f"    Min:  {min(hws_24)}")

    for threshold in [90, 100, 110, 120]:
        count = sum(1 for h in hws_24 if h <= threshold)
        print(f"    HW ≤ {threshold}: {count}/5000 = {count/50:.2f}%")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("3. СРАВНЕНИЕ ХВОСТОВ: 24r vs 64r vs random")
    print("=" * 70)

    # Random comparison: 5000 random δW
    hws_rand = []
    for _ in range(5000):
        dW = [np.random.randint(0, 2**32) for _ in range(16)]
        W2 = [W_base[i]^dW[i] for i in range(16)]
        H2 = sha256_full(W2)
        Hb_full = sha256_full(W_base)
        hws_rand.append(sum(hw(Hb_full[i]^H2[i]) for i in range(8)))

    print(f"\n  {'Source':<20} {'Mean':>6} {'Std':>6} {'Min':>5} {'P(≤100)':>8}")
    for name, hws in [("Random δW", hws_rand), ("GF2-kernel 64r", hws_64), ("GF2-kernel 24r", hws_24)]:
        p100 = sum(1 for h in hws if h <= 100) / len(hws) * 100
        print(f"  {name:<20} {np.mean(hws):>5.1f} {np.std(hws):>5.1f} {min(hws):>4} {p100:>7.2f}%")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("4. ЗНАЧЕНИЕ ХВОСТА")
    print("=" * 70)

    gain_64 = np.mean(hws_rand) - np.mean(hws_64)
    gain_24 = np.mean(hws_rand) - np.mean(hws_24)

    print(f"""
  CE HW DISTRIBUTION:
    Random δW:       mean={np.mean(hws_rand):.1f}, min={min(hws_rand)}
    GF2-kernel 64r:  mean={np.mean(hws_64):.1f}, min={min(hws_64)}
    GF2-kernel 24r:  mean={np.mean(hws_24):.1f}, min={min(hws_24)}

  GAIN from GF2-kernel (vs random):
    64r: {gain_64:+.1f} bits mean ({"BETTER" if gain_64 > 2 else "SAME"})
    24r: {gain_24:+.1f} bits mean

  TAIL (min HW):
    Random: {min(hws_rand)} (from 5000 samples)
    GF2-64: {min(hws_64)}
    GF2-24: {min(hws_24)}

  {"★ GF2-kernel gives LOWER min than random!" if min(hws_64) < min(hws_rand) - 5 else "GF2-kernel ≈ random (no tail advantage)."}

  BRIDGE TO REAL ATTACKS:
    Real attacks find δW with HW(δH) << 128 through SPECIFIC trails.
    Our GF2-kernel: HW(δH) = {min(hws_64)} minimum (from 5000).
    Differential cryptanalysis: HW(δH) = 0 (collision) at r≈28.

    Gap: our best {min(hws_64)} vs real 0.
    Real attacks use STRUCTURE (not GF2-kernel) to find low-HW δW.
    Our dimension captures AVERAGE (rank=256) not WORST-CASE (specific trail).
""")


if __name__ == "__main__":
    main()
