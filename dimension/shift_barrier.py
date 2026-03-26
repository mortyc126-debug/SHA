"""
СДВИГ БАРЬЕРА: выбор break point → σ₀(break) входит в ДРУГОЙ δW[r].

Schedule: δW[r] = σ₁(δW[r-2]) + δW[r-7] + σ₀(δW[r-15]) + δW[r-16]

σ₀(break) входит в δW[break+15]. Значит:
  Break=3: σ₀(break) → δW[18]     ← текущий barrier
  Break=4: σ₀(break) → δW[19]
  Break=5: σ₀(break) → δW[20]     ← pushed 2 rounds!
  Break=6: σ₀(break) → δW[21]
  Break=7: σ₀(break) → δW[22]     ← pushed 4 rounds!

Одновременно: break через σ₁ входит в δW[break+2]:
  Break=3: σ₁(break) → δW[5] (внутри a-repair, поглощается)
  Break=5: σ₁(break) → δW[7] (внутри a-repair)
  Break=7: σ₁(break) → δW[9] (внутри a-repair)

И break напрямую: δW[break+16] и δW[break+7]:
  Break=3: δW[19]=...+δW[3], δW[10]=...+δW[3]... нет, break=W[3], not free word
  Actually: break word is at position break_r. It enters δW[break_r+16] directly.
  Break=3: δW[3+16]=δW[19] gets +δW[3]
  Break=5: δW[5+16]=δW[21] gets +δW[5]
  Break=7: δW[7+16]=δW[23] gets +δW[7]

Каждый сдвиг break на +2 → все schedule contributions сдвигаются на +2!
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

def state_diff(s1, s2):
    return sum(hw(s1[i] ^ s2[i]) for i in range(8))


def test_break_point(break_r, break_bit=31, N=3000):
    """Test a-repair with break at different rounds."""
    meeting_chain_max = []
    dw_profiles = np.zeros((N, 64))
    profile_data = np.zeros((N, 65))
    dhs = []

    for trial in range(N):
        W1 = [np.random.randint(0, 2**32) for _ in range(16)]
        W1f = expand_schedule(W1)

        s1 = tuple(IV)
        states1 = [s1]
        for r in range(64):
            s1 = R(s1, W1f[r], r)
            states1.append(s1)

        W2 = list(W1)
        W2[break_r] ^= (1 << break_bit)

        s2 = tuple(IV)
        for r in range(break_r):
            s2 = R(s2, W2[r], r)
        s2 = R(s2, W2[break_r], break_r)

        for r in range(break_r + 1, 16):
            W2[r] = a_repair_W(s2, states1[r + 1][0], r)
            s2 = R(s2, W2[r], r)

        W2f = expand_schedule(W2)

        s2_full = tuple(IV)
        for r in range(64):
            w = W2[r] if r < 16 else W2f[r]
            s2_full = R(s2_full, w, r)
            profile_data[trial, r + 1] = state_diff(states1[r + 1], s2_full)

        for r in range(64):
            dw_profiles[trial, r] = hw(W1f[r] ^ W2f[r])

        H1 = tuple(add32(IV[i], states1[64][i]) for i in range(8))
        H2 = tuple(add32(IV[i], s2_full[i]) for i in range(8))
        dhs.append(sum(hw(H1[i] ^ H2[i]) for i in range(8)))

        # Max consecutive BOTH=0
        chain = 0
        mc = 0
        for r in range(65):
            if profile_data[trial, r] == 0:
                chain += 1
                mc = max(mc, chain)
            else:
                chain = 0
        meeting_chain_max.append(mc)

    return profile_data, dw_profiles, dhs, meeting_chain_max


def main():
    np.random.seed(42)

    print("=" * 70)
    print("СДВИГ БАРЬЕРА: break point → schedule barrier position")
    print("=" * 70)

    # Теория: σ₀(δW[break]) входит в δW[break+15]
    print(f"\n  Теоретический сдвиг:")
    for br in range(2, 11):
        sig0_entry = br + 15
        sig1_entry = br + 2
        direct_entry = br + 16
        direct7_entry = br + 7
        print(f"    break={br}: σ₀→δW[{sig0_entry}]  σ₁→δW[{sig1_entry}]  "
              f"direct→δW[{direct_entry}]  +7→δW[{direct7_entry}]")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("ТЕСТ: break r=2..10, bit=31, schedule transparency")
    print("=" * 70)

    results = {}

    for break_r in [2, 3, 4, 5, 6, 7, 8, 9, 10]:
        N = 3000
        profiles, dw_prof, dhs, chains = test_break_point(break_r, 31, N)

        mean_dw = np.mean(dw_prof, axis=0)
        mean_profile = np.mean(profiles, axis=0)
        both_pct = [np.mean(profiles[:, r] == 0) * 100 for r in range(65)]

        # Reboot round: first r > break+3 where BOTH=0 = 100%
        reboot = None
        for r in range(break_r + 4, 20):
            if both_pct[r] > 99:
                reboot = r
                break

        # Schedule transparency: consecutive δW=0 starting from r=16
        first_nonzero = 16
        for r in range(16, 25):
            if mean_dw[r] > 1:
                first_nonzero = r
                break

        # Transparent rounds
        transparent = first_nonzero - 16

        # Max meeting chain
        mean_chain = np.mean(chains)
        max_chain = max(chains)

        # r=17 pass rate
        r17_pass = sum(1 for i in range(N) if profiles[i, 17] == 0) / N * 100

        results[break_r] = {
            'reboot': reboot,
            'transparent': transparent,
            'first_nonzero': first_nonzero,
            'mean_chain': mean_chain,
            'max_chain': max_chain,
            'r17_pass': r17_pass,
            'min_dh': min(dhs),
            'mean_dw': mean_dw,
            'both_pct': both_pct,
        }

        marker = ""
        if transparent >= 2 and reboot:
            marker = " ★"
        if transparent >= 3:
            marker = " ★★"
        if transparent >= 4:
            marker = " ★★★"

        print(f"\n  break={break_r}: reboot={reboot}, transparent={transparent} (δW=0 for r=16..{first_nonzero-1}), "
              f"r17pass={r17_pass:.1f}%, chain={mean_chain:.1f}/{max_chain}{marker}")

        # Schedule detail
        print(f"    δW: ", end="")
        for r in range(16, 25):
            d = mean_dw[r]
            sym = "·" if d < 1 else "░" if d < 4 else "█"
            print(f"[{r}]={d:.1f}{sym} ", end="")
        print()

    # =================================================================
    print(f"\n{'=' * 70}")
    print("СВОДКА: break point → barrier shift")
    print("=" * 70)

    print(f"\n  {'Break':>5} {'Reboot':>7} {'Transp':>7} {'Barrier':>8} {'r17%':>6} {'Chain':>6} {'minδH':>6}")
    for br in sorted(results):
        r = results[br]
        barrier = r['first_nonzero']
        print(f"  {br:5d} {str(r['reboot']):>7} {r['transparent']:>7} r={barrier:>5} "
              f"{r['r17_pass']:>5.1f} {r['mean_chain']:>5.1f} {r['min_dh']:>5}")

    # Best
    best_br = max(results, key=lambda x: results[x]['transparent'] + (1 if results[x]['reboot'] else 0) * 5)
    print(f"\n  ЛУЧШИЙ: break={best_br}")
    r = results[best_br]
    print(f"    Reboot: r={r['reboot']}")
    print(f"    Schedule transparent: {r['transparent']} rounds (r=16..{r['first_nonzero']-1})")
    print(f"    Barrier at: δW[{r['first_nonzero']}]")
    print(f"    r=17 pass: {r['r17_pass']:.1f}%")

    # =================================================================
    print(f"\n{'=' * 70}")
    print("АНТИАРХИТЕКТУРНОЕ РЕШЕНИЕ")
    print("=" * 70)

    print(f"""
  SHA-256 архитектура: schedule deps at offsets -2, -7, -15, -16.
  σ₀(δW[break]) входит в δW[break+15].
  Break shift → barrier shift. Каждый +1 в break = +1 в barrier.

  Break=3: barrier at δW[18] (3 culprits: σ₁(δW[16]) + δW[11] + σ₀(δW[3]))
  Break=5: barrier at δW[20] (σ₀(δW[3])=0, σ₀(δW[5]) входит позже)
  Break=7: barrier at δW[22]

  НО: a-repair range сужается (break позже → меньше слов для repair).
  Break=3: a-repair r=4..15 (12 слов)
  Break=5: a-repair r=6..15 (10 слов)
  Break=7: a-repair r=8..15 (8 слов)

  Trade-off: barrier pushed +2 за каждые -2 слова a-repair.
  Reboot нужен 8 раундов → минимум 8 слов a-repair → break ≤ 7.

  ОПТИМУМ: break = max, при котором reboot ещё работает.
""")


if __name__ == "__main__":
    main()
