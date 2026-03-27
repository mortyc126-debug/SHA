"""
АНОМАЛИЯ 7: a-repair даёт HW(δH) = 119-121, consistently < 128.

Факт: когда мы forcing δa=0 для r=1..15, state сходится к δ=0 к r=8,
но после schedule rounds (r=16+) финальный HW(δH) ≈ 119-121.
Random δW даёт HW(δH) ≈ 128. Разница 7-9 бит.

Вопросы:
  A. Воспроизводится на разных seed?
  B. Statistically significant?
  C. ПОЧЕМУ ниже 128? Что a-repair оставляет в state?
  D. Можно ли усилить эффект?
"""

import numpy as np
import struct, hashlib

MASK32 = 0xFFFFFFFF

def hw(x): return bin(x).count('1')
def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def add32(x, y): return (x + y) & MASK32
def sub32(x, y): return (x - y) & MASK32
def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)

K = [
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
]
IV = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]


def sha256_full(W16):
    raw = struct.pack('>16I', *W16)
    return struct.unpack('>8I', hashlib.sha256(raw).digest())


def a_repair(W_base, dW0):
    """Force δa=0 for rounds 1..15. Return (W_prime, HW(δH))."""
    W_prime = list(W_base)
    W_prime[0] ^= dW0

    a1,b1,c1,d1,e1,f1,g1,h1 = IV
    a2,b2,c2,d2,e2,f2,g2,h2 = IV

    # Round 0: both with original W
    T1_1 = add32(add32(add32(add32(h1, Sigma1(e1)), Ch(e1,f1,g1)), K[0]), W_base[0])
    T2_1 = add32(Sigma0(a1), Maj(a1,b1,c1))
    h1,g1,f1,e1 = g1,f1,e1,add32(d1,T1_1)
    d1,c1,b1,a1 = c1,b1,a1,add32(T1_1,T2_1)

    T1_2 = add32(add32(add32(add32(h2, Sigma1(e2)), Ch(e2,f2,g2)), K[0]), W_prime[0])
    T2_2 = add32(Sigma0(a2), Maj(a2,b2,c2))
    h2,g2,f2,e2 = g2,f2,e2,add32(d2,T1_2)
    d2,c2,b2,a2 = c2,b2,a2,add32(T1_2,T2_2)

    # Rounds 1..15: force δa=0
    for r in range(1, 16):
        T1_1 = add32(add32(add32(add32(h1, Sigma1(e1)), Ch(e1,f1,g1)), K[r]), W_base[r])
        T2_1 = add32(Sigma0(a1), Maj(a1,b1,c1))
        a1_new = add32(T1_1, T2_1)

        # Force a2_new = a1_new
        T2_2 = add32(Sigma0(a2), Maj(a2,b2,c2))
        T1_2_needed = sub32(a1_new, T2_2)
        W_r_prime = sub32(sub32(sub32(sub32(T1_2_needed, h2), Sigma1(e2)), Ch(e2,f2,g2)), K[r])
        W_prime[r] = W_r_prime

        e1_new = add32(d1, T1_1)
        h1,g1,f1,e1 = g1,f1,e1,e1_new
        d1,c1,b1,a1 = c1,b1,a1,a1_new

        e2_new = add32(d2, T1_2_needed)
        h2,g2,f2,e2 = g2,f2,e2,e2_new
        d2,c2,b2,a2 = c2,b2,a2,a1_new

    H1 = sha256_full(W_base)
    H2 = sha256_full(W_prime)
    dH = sum(hw(H1[i] ^ H2[i]) for i in range(8))
    return W_prime, dH


def main():
    np.random.seed(42)

    print("=" * 70)
    print("АНОМАЛИЯ 7: a-repair gives HW(δH) ≈ 119-121")
    print("=" * 70)

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("A. 1000 random seeds × 32 bit positions")
    print("=" * 70)

    repair_hws = []
    random_hws = []

    for trial in range(1000):
        W_base = [np.random.randint(0, 2**32) for _ in range(16)]
        H_base = sha256_full(W_base)

        # a-repair with δW[0] = random single bit
        bit = np.random.randint(0, 32)
        dW0 = 1 << bit
        _, dH_repair = a_repair(W_base, dW0)
        repair_hws.append(dH_repair)

        # Random: just flip bit in W[0], no repair
        W_rand = list(W_base)
        W_rand[0] ^= dW0
        H_rand = sha256_full(W_rand)
        dH_rand = sum(hw(H_base[i] ^ H_rand[i]) for i in range(8))
        random_hws.append(dH_rand)

    from scipy import stats
    t, p = stats.ttest_ind(repair_hws, random_hws)

    print(f"  a-repair: mean={np.mean(repair_hws):.2f}, std={np.std(repair_hws):.2f}")
    print(f"  Random:   mean={np.mean(random_hws):.2f}, std={np.std(random_hws):.2f}")
    print(f"  Difference: {np.mean(random_hws) - np.mean(repair_hws):.2f} bits")
    print(f"  t={t:.2f}, p={p:.2e}")
    print(f"  → {'SIGNIFICANT!' if p < 0.01 else 'Not significant.'}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("B. Distribution comparison")
    print("=" * 70)

    for threshold in [100, 105, 110, 115, 120]:
        r_pct = sum(1 for h in repair_hws if h <= threshold) / len(repair_hws) * 100
        n_pct = sum(1 for h in random_hws if h <= threshold) / len(random_hws) * 100
        print(f"  HW ≤ {threshold}: repair={r_pct:.1f}%, random={n_pct:.1f}%")

    print(f"\n  Min: repair={min(repair_hws)}, random={min(random_hws)}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("C. WHY? What δ remains after repair converges?")
    print("=" * 70)

    # After a-repair, state converges to δ=0 by round 8.
    # But W_prime ≠ W_base for rounds 1..15.
    # Schedule: W[16+] = f(W[0..15]).
    # δW_prime ≠ 0 for W[1..15] → δW[16+] ≠ 0 → δ re-enters.

    # Measure: HW(δW) for schedule words
    W_base = [np.random.randint(0, 2**32) for _ in range(16)]
    W_prime, dH = a_repair(W_base, 1)

    # Expand both
    def expand(W16):
        W = list(W16)
        for r in range(16, 64):
            W.append((sigma1(W[r-2]) + W[r-7] + sigma0(W[r-15]) + W[r-16]) & MASK32)
        return W

    W_exp = expand(W_base)
    W_exp_p = expand(W_prime)

    print(f"  δW schedule after a-repair:")
    total_dW = 0
    for r in range(64):
        d = hw(W_exp[r] ^ W_exp_p[r])
        total_dW += d
        if r < 20 or d == 0:
            print(f"    δW[{r:>2}] HW = {d:>2}")

    print(f"  Total δW bits: {total_dW}")

    # Count non-zero δW in schedule
    nonzero_sched = sum(1 for r in range(16, 64) if W_exp[r] != W_exp_p[r])
    print(f"  Non-zero δW[16..63]: {nonzero_sched}/48")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("D. Control: random δW with SAME HW(δW)")
    print("=" * 70)

    # a-repair modifies W[0..15]. Total δW has specific HW.
    # What if we use RANDOM δW with same total HW?
    # If same → the bias is from HW(δW), not from structure.

    # Measure total HW(δW) from repair
    repair_dW_hws = []
    repair_dH_list = []
    random_matched_dH = []

    for trial in range(500):
        W_base = [np.random.randint(0, 2**32) for _ in range(16)]

        bit = np.random.randint(0, 32)
        W_prime, dH_r = a_repair(W_base, 1 << bit)

        dW_total_hw = sum(hw(W_base[i] ^ W_prime[i]) for i in range(16))
        repair_dW_hws.append(dW_total_hw)
        repair_dH_list.append(dH_r)

        # Random δW with same total HW
        # Create random δW by flipping dW_total_hw random bits
        W_rand = list(W_base)
        bits_to_flip = np.random.choice(512, dW_total_hw, replace=False)
        for bf in bits_to_flip:
            W_rand[bf // 32] ^= (1 << (bf % 32))
        H_base = sha256_full(W_base)
        H_rand = sha256_full(W_rand)
        dH_matched = sum(hw(H_base[i] ^ H_rand[i]) for i in range(8))
        random_matched_dH.append(dH_matched)

    print(f"  a-repair total HW(δW): mean={np.mean(repair_dW_hws):.1f}")
    print(f"  a-repair HW(δH):       mean={np.mean(repair_dH_list):.2f}")
    print(f"  Random matched HW(δH): mean={np.mean(random_matched_dH):.2f}")

    t2, p2 = stats.ttest_ind(repair_dH_list, random_matched_dH)
    print(f"  t={t2:.2f}, p={p2:.2e}")
    print(f"  → {'REPAIR IS SPECIAL' if p2 < 0.01 else 'Same as random with matched HW'}")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("E. Control: e-repair (force δe=0)")
    print("=" * 70)

    e_repair_hws = []
    for trial in range(500):
        W_base = [np.random.randint(0, 2**32) for _ in range(16)]
        W_prime = list(W_base)
        W_prime[0] ^= (1 << np.random.randint(0, 32))

        a1,b1,c1,d1,e1,f1,g1,h1 = IV
        a2,b2,c2,d2,e2,f2,g2,h2 = IV

        T1_1 = add32(add32(add32(add32(h1,Sigma1(e1)),Ch(e1,f1,g1)),K[0]),W_base[0])
        T2_1 = add32(Sigma0(a1),Maj(a1,b1,c1))
        h1,g1,f1,e1 = g1,f1,e1,add32(d1,T1_1)
        d1,c1,b1,a1 = c1,b1,a1,add32(T1_1,T2_1)

        T1_2 = add32(add32(add32(add32(h2,Sigma1(e2)),Ch(e2,f2,g2)),K[0]),W_prime[0])
        T2_2 = add32(Sigma0(a2),Maj(a2,b2,c2))
        h2,g2,f2,e2 = g2,f2,e2,add32(d2,T1_2)
        d2,c2,b2,a2 = c2,b2,a2,add32(T1_2,T2_2)

        for r in range(1, 16):
            T1_1 = add32(add32(add32(add32(h1,Sigma1(e1)),Ch(e1,f1,g1)),K[r]),W_base[r])
            e1_new = add32(d1, T1_1)

            # Force e2_new = e1_new
            T1_2_needed = sub32(e1_new, d2)
            W_r_prime = sub32(sub32(sub32(sub32(T1_2_needed,h2),Sigma1(e2)),Ch(e2,f2,g2)),K[r])
            W_prime[r] = W_r_prime

            T2_1 = add32(Sigma0(a1),Maj(a1,b1,c1))
            h1,g1,f1,e1 = g1,f1,e1,e1_new
            d1,c1,b1,a1 = c1,b1,a1,add32(T1_1,T2_1)

            T2_2 = add32(Sigma0(a2),Maj(a2,b2,c2))
            h2,g2,f2,e2 = g2,f2,e2,add32(d2,T1_2_needed)
            d2,c2,b2,a2 = c2,b2,a2,add32(T1_2_needed,T2_2)

        H1 = sha256_full(W_base)
        H2 = sha256_full(W_prime)
        e_repair_hws.append(sum(hw(H1[i]^H2[i]) for i in range(8)))

    print(f"  e-repair HW(δH): mean={np.mean(e_repair_hws):.2f}, std={np.std(e_repair_hws):.2f}")
    print(f"  a-repair HW(δH): mean={np.mean(repair_dH_list):.2f}")
    print(f"  → a-repair {'BETTER' if np.mean(repair_dH_list) < np.mean(e_repair_hws) - 1 else 'SAME as'} e-repair")

    # ═══════════════════
    print(f"\n{'=' * 70}")
    print("ВЕРДИКТ")
    print("=" * 70)


if __name__ == "__main__":
    main()
