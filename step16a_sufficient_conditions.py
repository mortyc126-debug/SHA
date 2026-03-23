#!/usr/bin/env python3
"""Step 16a: Wang-style sufficient conditions for SHA-256.

For a chosen differential path DW[0]=0x80000000,
derive the per-round conditions on state registers
that guarantee the differential propagates correctly.

Wang's approach:
1. Choose a differential characteristic (which bits differ each round)
2. Derive sufficient conditions on intermediate values
3. Use message modification to satisfy conditions
"""

import random
MASK32 = 0xFFFFFFFF

def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK32
def shr(x, n): return x >> n
def sigma0(x): return rotr(x, 7) ^ rotr(x, 18) ^ shr(x, 3)
def sigma1(x): return rotr(x, 17) ^ rotr(x, 19) ^ shr(x, 10)
def Sigma0(x): return rotr(x, 2) ^ rotr(x, 13) ^ rotr(x, 22)
def Sigma1(x): return rotr(x, 6) ^ rotr(x, 11) ^ rotr(x, 25)
def Ch(e, f, g): return (e & f) ^ (~e & g) & MASK32
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)
def hw(x): return bin(x).count('1')

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

def expand_message(W16):
    W = list(W16)
    for t in range(16, 64):
        W.append((sigma1(W[t-2]) + W[t-7] + sigma0(W[t-15]) + W[t-16]) & MASK32)
    return W

def sha256_round(state, W_t, K_t):
    a, b, c, d, e, f, g, h = state
    T1 = (h + Sigma1(e) + Ch(e, f, g) + K_t + W_t) & MASK32
    T2 = (Sigma0(a) + Maj(a, b, c)) & MASK32
    return ((T1 + T2) & MASK32, a, b, c, (d + T1) & MASK32, e, f, g)

# ============================================================
# Part 1: Empirical differential path extraction
# ============================================================
def extract_differential_path(N=100000):
    """Find the most common per-bit differential at each round."""
    print("=" * 70)
    print("Part 1: Extracting empirical differential path")
    print("=" * 70)

    DW0 = 0x80000000

    # For each round, track per-bit XOR differential frequencies
    # bit_freq[round][reg][bit] = count of times this bit differs
    bit_freq = {}
    for t in range(25):
        bit_freq[t] = [[0]*32 for _ in range(8)]

    for _ in range(N):
        W = [random.randint(0, MASK32) for _ in range(16)]
        W2 = list(W)
        W2[0] ^= DW0

        Wexp = expand_message(W)
        W2exp = expand_message(W2)

        s = list(IV)
        s2 = list(IV)
        for t in range(25):
            s = list(sha256_round(s, Wexp[t], K[t]))
            s2 = list(sha256_round(s2, W2exp[t], K[t]))
            for reg in range(8):
                xd = s[reg] ^ s2[reg]
                for bit in range(32):
                    if (xd >> bit) & 1:
                        bit_freq[t][reg][bit] += 1

    print(f"\nPer-round state XOR differential (P(bit differs) for registers a,e):")
    print(f"{'Round':>5} | {'Da HW (avg)':>11} | {'De HW (avg)':>11} | {'Da top bits':>30}")
    print("-" * 80)

    diff_path = {}  # round -> (Da_mask, De_mask) for high-probability bits

    for t in range(25):
        da_hw = sum(1 for b in range(32) if bit_freq[t][0][b] > N*0.4)
        de_hw = sum(1 for b in range(32) if bit_freq[t][4][b] > N*0.4)

        # Build masks of bits that differ with P > 0.9 (deterministic)
        da_det = 0
        de_det = 0
        da_high = []
        for b in range(32):
            pa = bit_freq[t][0][b] / N
            pe = bit_freq[t][4][b] / N
            if pa > 0.9:
                da_det |= (1 << b)
            if pe > 0.9:
                de_det |= (1 << b)
            if pa > 0.4:
                da_high.append(f"b{b}({pa:.2f})")

        diff_path[t] = (da_det, de_det)

        # Show avg HW including probabilistic bits
        da_avg = sum(bit_freq[t][0][b]/N for b in range(32))
        de_avg = sum(bit_freq[t][4][b]/N for b in range(32))

        det_info = f"det_Da=0x{da_det:08x}(HW{hw(da_det)}) det_De=0x{de_det:08x}(HW{hw(de_det)})"
        print(f"{t:5d} | {da_avg:11.2f} | {de_avg:11.2f} | {det_info}")

    return diff_path

# ============================================================
# Part 2: Sufficient conditions from Ch and Maj
# ============================================================
def derive_ch_conditions():
    """Derive conditions on e,f,g for Ch differential to behave."""
    print("\n" + "=" * 70)
    print("Part 2: Ch/Maj sufficient conditions (per-bit)")
    print("=" * 70)

    # Ch(e,f,g) = (e AND f) XOR (NOT e AND g)
    # DCh = Ch(e+De, f+Df, g+Dg) - Ch(e, f, g)
    #
    # For XOR differential De (ignoring carries for now):
    # If De[i]=0: DCh[i] = 0 (no change at bit i)
    # If De[i]=1, Df[i]=0, Dg[i]=0:
    #   DCh[i] = f[i] XOR g[i]  (flips between selecting f and g)
    #   For DCh[i]=0 we need f[i]=g[i]
    #
    # If De[i]=1, Df[i]=1, Dg[i]=0:
    #   Old: e[i]*f[i] + (1-e[i])*g[i]
    #   New: (1-e[i])*(1-f[i]) + e[i]*g[i]
    #   DCh[i] = depends...

    print("\nCh truth table with single-bit De=1, Df=Dg=0:")
    print(f"  e f g | Ch | Ch' (e'=1-e) | DCh")
    print("-" * 40)
    for e in range(2):
        for f in range(2):
            for g in range(2):
                ch = (e & f) ^ ((1-e) & g)
                ch2 = ((1-e) & f) ^ (e & g)
                dch = ch ^ ch2
                cond = "f=g needed for DCh=0" if dch else "OK"
                print(f"  {e} {f} {g} |  {ch} |      {ch2}      |  {dch}  {cond}")

    print("\n  => Sufficient condition for DCh[i]=0 when De[i]=1:")
    print("     f[i] = g[i]")
    print("  => This is a condition on the STATE, not the message!")

    print("\nMaj truth table with single-bit Da=1, Db=Dc=0:")
    print(f"  a b c | Maj | Maj' (a'=1-a) | DMaj")
    print("-" * 42)
    for a in range(2):
        for b in range(2):
            for c in range(2):
                maj = (a & b) ^ (a & c) ^ (b & c)
                maj2 = ((1-a) & b) ^ ((1-a) & c) ^ (b & c)
                dmaj = maj ^ maj2
                cond = "b=c needed" if dmaj else "OK"
                print(f"  {a} {b} {c} |  {maj}  |      {maj2}       |  {dmaj}   {cond}")

    print("\n  => Sufficient condition for DMaj[i]=0 when Da[i]=1:")
    print("     b[i] = c[i]")

# ============================================================
# Part 3: Count conditions needed per round
# ============================================================
def count_conditions(diff_path, N=50000):
    """Count how many sufficient conditions are needed and satisfiable."""
    print("\n" + "=" * 70)
    print("Part 3: Condition count per round")
    print("=" * 70)

    DW0 = 0x80000000

    # For each round, count:
    # - Number of bits where De differs (need f[i]=g[i])
    # - Number of bits where Da differs (need b[i]=c[i])
    # - How often these conditions are naturally satisfied

    cond_stats = {}

    for _ in range(N):
        W = [random.randint(0, MASK32) for _ in range(16)]
        W2 = list(W)
        W2[0] ^= DW0

        Wexp = expand_message(W)
        W2exp = expand_message(W2)

        s = list(IV)
        s2 = list(IV)

        for t in range(20):
            s = list(sha256_round(s, Wexp[t], K[t]))
            s2 = list(sha256_round(s2, W2exp[t], K[t]))

            a, b, c, d, e, f, g, h = s
            da_xor = s[0] ^ s2[0]
            de_xor = s[4] ^ s2[4]

            # Ch conditions: for bits where De differs, need f=g
            ch_cond_needed = hw(de_xor)
            ch_cond_sat = 0
            for bit in range(32):
                if (de_xor >> bit) & 1:
                    if ((f >> bit) & 1) == ((g >> bit) & 1):
                        ch_cond_sat += 1

            # Maj conditions: for bits where Da differs, need b=c
            maj_cond_needed = hw(da_xor)
            maj_cond_sat = 0
            for bit in range(32):
                if (da_xor >> bit) & 1:
                    if ((b >> bit) & 1) == ((c >> bit) & 1):
                        maj_cond_sat += 1

            if t not in cond_stats:
                cond_stats[t] = {'ch_need': [], 'ch_sat': [], 'maj_need': [], 'maj_sat': []}
            cond_stats[t]['ch_need'].append(ch_cond_needed)
            cond_stats[t]['ch_sat'].append(ch_cond_sat)
            cond_stats[t]['maj_need'].append(maj_cond_needed)
            cond_stats[t]['maj_sat'].append(maj_cond_sat)

    print(f"\n{'Round':>5} | {'Ch conds':>8} | {'Ch sat%':>7} | {'Maj conds':>9} | {'Maj sat%':>8} | {'Total unsatisfied':>17}")
    print("-" * 65)
    for t in range(20):
        cs = cond_stats[t]
        ch_n = sum(cs['ch_need']) / N
        ch_s = sum(cs['ch_sat']) / N
        maj_n = sum(cs['maj_need']) / N
        maj_s = sum(cs['maj_sat']) / N
        ch_pct = (ch_s / ch_n * 100) if ch_n > 0 else 100
        maj_pct = (maj_s / maj_n * 100) if maj_n > 0 else 100
        unsat = (ch_n - ch_s) + (maj_n - maj_s)
        print(f"{t:5d} | {ch_n:8.1f} | {ch_pct:6.1f}% | {maj_n:9.1f} | {maj_pct:7.1f}% | {unsat:17.1f}")

# ============================================================
# Part 4: Which message words can fix which conditions?
# ============================================================
def message_to_condition_map():
    """Map which W[i] bits affect which state conditions."""
    print("\n" + "=" * 70)
    print("Part 4: Message word -> state condition dependency")
    print("=" * 70)

    # In round t: W[t] directly affects a_t and e_t
    # a_t = T1 + T2, where T1 = h + Sigma1(e) + Ch(e,f,g) + K[t] + W[t]
    # e_t = d + T1
    #
    # So W[t] has direct linear influence on a_t and e_t
    # Flipping W[t] bit j flips a_t and e_t (plus carry effects)
    #
    # But conditions at round t involve state from round t-1:
    # f_t = e_{t-1}, g_t = e_{t-2}, b_t = a_{t-1}, c_t = a_{t-2}
    #
    # So conditions at round t+1 (f=g for Ch) require e_t = e_{t-1}
    # at specific bit positions. We can influence e_t via W[t].

    print("\n  Round t state mapping:")
    print("  a[t] = T1 + T2  (depends on W[t])")
    print("  e[t] = d[t-1] + T1  (depends on W[t])")
    print("  b[t] = a[t-1]  (depends on W[t-1])")
    print("  c[t] = a[t-2]  (depends on W[t-2])")
    print("  f[t] = e[t-1]  (depends on W[t-1])")
    print("  g[t] = e[t-2]  (depends on W[t-2])")
    print()
    print("  => W[t] can directly modify a[t] and e[t]")
    print("  => To fix Ch condition at round t+1 (f[t+1]=g[t+1]):")
    print("     Need e[t] = e[t-1] at specific bits")
    print("     W[t] controls e[t], W[t-1] controls e[t-1]")
    print("  => To fix Maj condition at round t+1 (b[t+1]=c[t+1]):")
    print("     Need a[t] = a[t-1] at specific bits")
    print("     W[t] controls a[t], W[t-1] controls a[t-1]")
    print()
    print("  KEY INSIGHT: For rounds 0-15, W[t] is a FREE variable")
    print("  For rounds 16+, W[t] is determined by message schedule")
    print("  => Direct message modification works for rounds 0-16")
    print("  => After round 16, need indirect or probabilistic methods")

    # Verify: flipping W[t] bit j, how does e[t] change?
    print("\n  Verifying W[t] -> e[t] sensitivity:")
    for t_test in [0, 5, 10, 15]:
        sensitivities = []
        for bit in range(32):
            same_count = 0
            N = 1000
            for _ in range(N):
                W = [random.randint(0, MASK32) for _ in range(16)]
                Wf = list(W)
                Wf[t_test] ^= (1 << bit)

                Wexp = expand_message(W)
                Wfexp = expand_message(Wf)

                s = list(IV)
                sf = list(IV)
                for t in range(t_test + 1):
                    s = list(sha256_round(s, Wexp[t], K[t]))
                    sf = list(sha256_round(sf, Wfexp[t], K[t]))

                # Check if e[t_test] changed at bit position
                de = s[4] ^ sf[4]
                if (de >> bit) & 1:
                    same_count += 1
            sensitivities.append(same_count / N)

        avg_sens = sum(sensitivities) / 32
        high_sens = sum(1 for s in sensitivities if s > 0.9)
        print(f"    W[{t_test:2d}]: avg bit sensitivity = {avg_sens:.3f}, "
              f"bits with >90% sensitivity: {high_sens}/32")


if __name__ == "__main__":
    random.seed(42)

    diff_path = extract_differential_path(N=50000)
    derive_ch_conditions()
    count_conditions(diff_path, N=20000)
    message_to_condition_map()

    print("\n" + "=" * 70)
    print("Step 16a Complete")
    print("=" * 70)
