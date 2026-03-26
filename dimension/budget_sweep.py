"""
ПОЛНЫЙ ПЕРЕБОР БЮДЖЕТА: M слов meeting + (12-M) слов lock.

Budget: W[4..15] = 12 позиций.
  M позиций → a-repair (W[4..3+M])
  12-M позиций → schedule lock (W[4+M..15] = W1)

Для каждого M=0..12:
  Meeting: a-repair на r=4..3+M
  Lock: δW[4+M..15] = 0 (SET to W1)
  Schedule: δW[16..] depends on all δW[0..15]

Ключевые метрики:
  Reboot? (meeting at some round)
  δW[16..20] (schedule transparency)
  Last BOTH=0
  HW(δH)
"""

import numpy as np

MASK32 = 0xFFFFFFFF
K = [
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
def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
def add32(x, y): return (x + y) & MASK32
def sub32(x, y): return (x - y) & MASK32
def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ (x >> 3)
def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ (x >> 10)

def R(state, W_r, r_idx):
    a, b, c, d, e, f, g, h = state
    T1 = add32(add32(add32(add32(h, Sigma1(e)), Ch(e, f, g)), K[r_idx]), W_r)
    T2 = add32(Sigma0(a), Maj(a, b, c))
    return (add32(T1, T2), a, b, c, add32(d, T1), e, f, g)

def expand_schedule(W16):
    W = list(W16)
    for r in range(16, 64):
        W.append(add32(add32(add32(sigma1(W[r-2]), W[r-7]), sigma0(W[r-15])), W[r-16]))
    return W

def a_repair_W(state2, target_a, r_idx):
    a, b, c, d, e, f, g, h = state2
    T2 = add32(Sigma0(a), Maj(a, b, c))
    T1_needed = sub32(target_a, T2)
    return sub32(sub32(sub32(sub32(T1_needed, h), Sigma1(e)), Ch(e, f, g)), K[r_idx])


def test_budget(M, N=1000):
    """Test M words a-repair + (12-M) words lock."""
    results = {
        'reboot': 0, 'reboot_round': [],
        'dw16': [], 'dw17': [], 'dw18': [], 'dw19': [], 'dw20': [],
        'max_both': [], 'dh': [], 'last_zero': [],
    }

    for trial in range(N):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W1f = expand_schedule(W1)
        s1 = tuple(IV)
        states1 = [s1]
        for r in range(64):
            s1 = R(s1, W1f[r], r)
            states1.append(s1)

        W2 = list(W1)
        W2[3] ^= (1 << 31)

        s2 = tuple(IV)
        for r in range(3):
            s2 = R(s2, W2[r], r)
        s2 = R(s2, W2[3], 3)

        # A-repair on W[4..3+M]
        repair_end = min(4 + M, 16)
        for r in range(4, repair_end):
            W2[r] = a_repair_W(s2, states1[r + 1][0], r)
            s2 = R(s2, W2[r], r)

        # Lock on W[repair_end..15] = W1
        for r in range(repair_end, 16):
            W2[r] = W1[r]
            s2 = R(s2, W2[r], r)

        W2f = expand_schedule(W2)

        # Profile
        s2_check = tuple(IV)
        profile = []
        for r in range(64):
            w = W2[r] if r < 16 else W2f[r]
            s2_check = R(s2_check, w, r)
            profile.append(sum(hw(states1[r+1][i] ^ s2_check[i]) for i in range(8)))

        # Metrics
        for r in range(5, 20):
            if profile[r] == 0 and r > 6:
                results['reboot'] += 1
                results['reboot_round'].append(r + 1)
                break

        results['dw16'].append(hw(W1f[16] ^ W2f[16]))
        results['dw17'].append(hw(W1f[17] ^ W2f[17]))
        results['dw18'].append(hw(W1f[18] ^ W2f[18]))
        results['dw19'].append(hw(W1f[19] ^ W2f[19]))
        results['dw20'].append(hw(W1f[20] ^ W2f[20]))

        max_b = 0
        chain = 0
        for r in range(64):
            if profile[r] == 0:
                chain += 1
                max_b = max(max_b, chain)
            else:
                chain = 0
        results['max_both'].append(max_b)

        lz = 0
        for r in range(64):
            if profile[r] == 0:
                lz = r + 1
        results['last_zero'].append(lz)

        H1 = tuple(add32(IV[i], states1[64][i]) for i in range(8))
        H2 = tuple(add32(IV[i], s2_check[i]) for i in range(8))
        results['dh'].append(sum(hw(H1[i] ^ H2[i]) for i in range(8)))

    return results


def main():
    np.random.seed(42)

    print("=" * 70)
    print("ПОЛНЫЙ ПЕРЕБОР БЮДЖЕТА: M meeting + (12-M) lock")
    print("=" * 70)

    N = 2000

    print(f"\n  {'M':>3} {'Repair':>8} {'Lock':>6} {'Reboot%':>8} {'δW16':>5} {'δW17':>5} "
          f"{'δW18':>5} {'δW19':>5} {'δW20':>5} {'Chain':>6} {'LastZ':>6} {'minδH':>6}")
    print(f"  {'-'*3} {'-'*8} {'-'*6} {'-'*8} {'-'*5} {'-'*5} "
          f"{'-'*5} {'-'*5} {'-'*5} {'-'*6} {'-'*6} {'-'*6}")

    all_results = {}

    for M in range(13):  # 0..12
        r = test_budget(M, N)
        all_results[M] = r

        repair_range = f"W[4..{3+M}]" if M > 0 else "none"
        lock_range = f"W[{4+M}..15]" if M < 12 else "none"
        reboot_pct = r['reboot'] / N * 100

        marker = ""
        # Best = high reboot + low δW[18]
        if reboot_pct > 50 and np.mean(r['dw18']) < 8:
            marker = " ★★★"
        elif reboot_pct > 50 or np.mean(r['dw18']) < 10:
            marker = " ★★"
        elif reboot_pct > 10 or np.mean(r['dw18']) < 12:
            marker = " ★"

        print(f"  {M:3d} {repair_range:>8} {lock_range:>6} {reboot_pct:7.1f}% "
              f"{np.mean(r['dw16']):4.1f} {np.mean(r['dw17']):4.1f} "
              f"{np.mean(r['dw18']):4.1f} {np.mean(r['dw19']):4.1f} {np.mean(r['dw20']):4.1f} "
              f"{np.mean(r['max_both']):5.1f} {np.mean(r['last_zero']):5.1f} "
              f"{min(r['dh']):5d}{marker}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("ОПТИМАЛЬНОЕ РАЗБИЕНИЕ")
    print("=" * 70)

    # Score = reboot_pct + (16 - dw18) * 5 + last_zero
    best_M = -1
    best_score = -999

    for M in range(13):
        r = all_results[M]
        reboot = r['reboot'] / N * 100
        dw18 = np.mean(r['dw18'])
        last_z = np.mean(r['last_zero'])
        chain = np.mean(r['max_both'])

        # Score: we want BOTH reboot AND low schedule damage
        score = (reboot / 10) + (16 - dw18) * 3 + last_z + chain * 2

        if score > best_score:
            best_score = score
            best_M = M

    r = all_results[best_M]
    print(f"\n  Лучший M = {best_M} (repair W[4..{3+best_M}], lock W[{4+best_M}..15])")
    print(f"    Reboot: {r['reboot']/N*100:.1f}%")
    print(f"    δW[16]: {np.mean(r['dw16']):.1f}")
    print(f"    δW[17]: {np.mean(r['dw17']):.1f}")
    print(f"    δW[18]: {np.mean(r['dw18']):.1f}")
    print(f"    Chain: {np.mean(r['max_both']):.1f}")
    print(f"    Last BOTH=0: {np.mean(r['last_zero']):.1f}")
    print(f"    HW(δH): mean={np.mean(r['dh']):.1f}, min={min(r['dh'])}")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("TRADE-OFF ВИЗУАЛИЗАЦИЯ")
    print("=" * 70)

    print(f"\n  M  Meeting────────────  Lock────────────")
    for M in range(13):
        r = all_results[M]
        reboot = r['reboot'] / N * 100
        dw18 = np.mean(r['dw18'])

        meet_bar = "█" * int(reboot / 10)
        lock_bar = "█" * int(max(0, 16 - dw18))

        print(f"  {M:2d} [{meet_bar:<10}] {reboot:5.1f}%  [{lock_bar:<16}] δW18={dw18:.1f}")


if __name__ == "__main__":
    main()
