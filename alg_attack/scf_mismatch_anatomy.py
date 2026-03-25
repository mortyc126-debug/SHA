#!/usr/bin/env python3
"""
SCF: АНАТОМИЯ MISMATCH — какие биты, какие раунды, какие слова.

Цель: открытие структуры mismatch → создание точного инструмента.

1. Bit-position profile: какие биты mismatch'а чаще ненулевые?
2. Influence map: бит W[word][bit] → Δ(mismatch) на каждом раунде
3. Positive interference: биты помогающие НЕСКОЛЬКИМ раундам
4. Coupling: корреляция mismatch между соседними раундами
"""
import os, sys

MASK32 = 0xFFFFFFFF
K = [
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2,
]
IV = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
      0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&MASK32
def Sig0(x): return rotr(x,2)^rotr(x,13)^rotr(x,22)
def Sig1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def sig0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def sig1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def Ch(e,f,g): return (e&f)^(~e&g)&MASK32
def Maj(a,b,c): return (a&b)^(a&c)^(b&c)
def add32(*args):
    s=0
    for x in args: s=(s+x)&MASK32
    return s
def sub32(a,b): return (a-b)&MASK32
def hw(x): return bin(x).count('1')

def expand_real(W16):
    W=list(W16)
    for i in range(16,64):
        W.append(add32(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    return W

def sha_all_states(W16):
    W=expand_real(W16); s=list(IV); states=[list(s)]
    for r in range(64):
        a,b,c,d,e,f,g,h=s
        T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add32(Sig0(a),Maj(a,b,c))
        s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
        states.append(list(s))
    return states, W

def compute_mismatch(W1, W2):
    """Mismatch[r] = required_W2[r] XOR actual_W2[r] for gap rounds."""
    states1, Wexp1 = sha_all_states(W1)
    states2, Wexp2 = sha_all_states(W2)
    mismatches = {}
    for r in range(17, 52):
        s1r = states1[r]; s2r = states2[r]
        T1_1 = sub32(states1[r+1][4], s1r[3])
        T1_2_req = add32(T1_1, sub32(s1r[3], s2r[3]))
        W2_req = sub32(sub32(sub32(sub32(T1_2_req,
                  s2r[7]), Sig1(s2r[4])), Ch(s2r[4],s2r[5],s2r[6])), K[r])
        mismatches[r] = W2_req ^ Wexp2[r]
    return mismatches


# ============================================================
# Discovery 1: Bit-position profile of mismatch
# ============================================================
def disc1_bit_profile(N):
    print("="*70)
    print("DISCOVERY 1: BIT-POSITION PROFILE OF MISMATCH")
    print("  P(mismatch[r][bit]=1) for each bit position 0..31")
    print("="*70)

    bit_counts = [[0]*32 for _ in range(35)]  # 35 gap rounds × 32 bits

    for _ in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = list(W1); W2[0] ^= int.from_bytes(os.urandom(4),'big')
        if W2[0]==W1[0]: W2[0]^=1

        mm = compute_mismatch(W1, W2)
        for ri, r in enumerate(range(17, 52)):
            for b in range(32):
                if (mm[r] >> b) & 1:
                    bit_counts[ri][b] += 1

    # Average across rounds
    avg_per_bit = [sum(bit_counts[ri][b] for ri in range(35)) / (35*N)
                   for b in range(32)]

    print(f"\n  P(mismatch[bit]=1) averaged over rounds (N={N}):")
    print(f"  {'bit':>4} | P(=1)  | deviation from 0.5")
    print("  " + "-"*40)
    for b in range(32):
        p = avg_per_bit[b]
        dev = p - 0.5
        bar = "+" * int(abs(dev)*200) if dev > 0 else "-" * int(abs(dev)*200)
        marker = " ★" if abs(dev) > 0.03 else ""
        print(f"  {b:4d} | {p:.4f} | {dev:+.4f} {bar}{marker}")

    # Key question: is bit 0 special? (carry-free)
    print(f"\n  Bit 0 (carry-free LSB): P = {avg_per_bit[0]:.4f}")
    print(f"  Bit 1 (carry-biased):   P = {avg_per_bit[1]:.4f}")
    print(f"  Average bit 2-31:       P = {sum(avg_per_bit[2:])/30:.4f}")


# ============================================================
# Discovery 2: Influence map — which W[word][bit] helps most?
# ============================================================
def disc2_influence_map(N):
    print("\n" + "="*70)
    print("DISCOVERY 2: INFLUENCE MAP — Δ(total_mismatch) per input bit")
    print("  Flip W2[word][bit], measure change in total mismatch")
    print("="*70)

    W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    W2 = list(W1); W2[0] ^= 0x80000000

    base_mm = compute_mismatch(W1, W2)
    base_total = sum(hw(base_mm[r]) for r in range(17, 52))

    # For each word and sampled bits
    influence = {}  # (word, bit) → avg Δmismatch

    print(f"\n  Base mismatch: {base_total}")
    print(f"  {'word':>5} | best_bit Δmismatch | worst_bit | avg")
    print("  " + "-"*55)

    for word in range(16):
        word_influences = []
        for bit in range(32):
            W2_flip = list(W2)
            W2_flip[word] ^= (1 << bit)
            mm_flip = compute_mismatch(W1, W2_flip)
            flip_total = sum(hw(mm_flip[r]) for r in range(17, 52))
            delta = flip_total - base_total
            word_influences.append((delta, bit))
            influence[(word, bit)] = delta

        word_influences.sort()
        best = word_influences[0]
        worst = word_influences[-1]
        avg = sum(d for d, _ in word_influences) / 32

        marker = " ★" if best[0] < -5 else ""
        print(f"  W[{word:2d}] | bit {best[1]:2d}: {best[0]:+4d}      | "
              f"bit {worst[1]:2d}: {worst[0]:+4d}  | {avg:+.1f}{marker}")

    # Top 10 most helpful bits globally
    print(f"\n  Top 10 most HELPFUL bit flips:")
    all_influences = sorted(influence.items(), key=lambda x: x[1])
    for (w, b), delta in all_influences[:10]:
        print(f"    W[{w}][{b:2d}]: Δmismatch = {delta:+d}")

    # Top 10 most HARMFUL
    print(f"\n  Top 10 most HARMFUL bit flips:")
    for (w, b), delta in all_influences[-10:]:
        print(f"    W[{w}][{b:2d}]: Δmismatch = {delta:+d}")

    return all_influences


# ============================================================
# Discovery 3: Positive interference — helps multiple rounds
# ============================================================
def disc3_interference(N):
    print("\n" + "="*70)
    print("DISCOVERY 3: POSITIVE INTERFERENCE")
    print("  Bits that reduce mismatch on MULTIPLE rounds simultaneously")
    print("="*70)

    W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    W2 = list(W1); W2[0] ^= 0x80000000

    base_mm = compute_mismatch(W1, W2)
    base_per_round = {r: hw(base_mm[r]) for r in range(17, 52)}

    # For each input bit: count rounds where it HELPS (reduces mismatch)
    print(f"\n  {'word':>5} {'bit':>4} | rounds_helped | rounds_hurt | net | Δtotal")
    print("  " + "-"*65)

    best_helpers = []

    for word in range(16):
        for bit in [0, 7, 15, 23, 31]:  # Sample 5 bits per word
            W2_flip = list(W2)
            W2_flip[word] ^= (1 << bit)
            mm_flip = compute_mismatch(W1, W2_flip)

            helped = 0; hurt = 0; total_delta = 0
            for r in range(17, 52):
                d = hw(mm_flip[r]) - hw(base_mm[r])
                total_delta += d
                if d < 0: helped += 1
                if d > 0: hurt += 1

            net = helped - hurt
            best_helpers.append((net, total_delta, word, bit, helped, hurt))

            if net > 5 or total_delta < -10:
                print(f"  W[{word:2d}] {bit:4d} | {helped:13d} | {hurt:11d} | {net:+3d} | {total_delta:+4d} ★")

    best_helpers.sort(reverse=True)
    print(f"\n  Top 5 by net rounds helped:")
    for net, td, w, b, h, hu in best_helpers[:5]:
        print(f"    W[{w}][{b:2d}]: net={net:+d}, Δtotal={td:+d}, helped={h}, hurt={hu}")


# ============================================================
# Discovery 4: Round coupling — mismatch correlation
# ============================================================
def disc4_coupling(N):
    print("\n" + "="*70)
    print("DISCOVERY 4: ROUND COUPLING")
    print("  Correlation of mismatch between adjacent rounds")
    print("  High correlation → they change TOGETHER → can't fix independently")
    print("="*70)

    # Collect per-round mismatch HW vectors
    vectors = []
    for _ in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = list(W1); W2[0] ^= int.from_bytes(os.urandom(4),'big')
        if W2[0]==W1[0]: W2[0]^=1
        mm = compute_mismatch(W1, W2)
        vectors.append([hw(mm[r]) for r in range(17, 52)])

    # Correlation between round r and round r+k
    print(f"\n  Correlation(HW(mm[r]), HW(mm[r+k])):")
    for k in [1, 2, 3, 4, 7]:
        corrs = []
        for ri in range(35-k):
            xs = [v[ri] for v in vectors]
            ys = [v[ri+k] for v in vectors]
            mx = sum(xs)/len(xs); my = sum(ys)/len(ys)
            cov = sum((x-mx)*(y-my) for x,y in zip(xs,ys))/len(xs)
            sx = (sum((x-mx)**2 for x in xs)/len(xs))**0.5
            sy = (sum((y-my)**2 for y in ys)/len(ys))**0.5
            if sx > 0 and sy > 0:
                corrs.append(cov/(sx*sy))
        avg_corr = sum(corrs)/len(corrs) if corrs else 0
        marker = " ★" if abs(avg_corr) > 0.1 else ""
        print(f"    k={k}: avg_corr = {avg_corr:+.4f}{marker}")

    # XOR coupling: are specific BITS of mismatch correlated across rounds?
    print(f"\n  Bit-level coupling (bit 0 of round r vs bit 0 of round r+1):")
    for test_bit in [0, 1, 15, 31]:
        match_count = 0; total = 0
        for _ in range(N):
            W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
            W2 = list(W1); W2[0] ^= int.from_bytes(os.urandom(4),'big')
            if W2[0]==W1[0]: W2[0]^=1
            mm = compute_mismatch(W1, W2)
            for r in range(17, 51):
                b_r = (mm[r] >> test_bit) & 1
                b_r1 = (mm[r+1] >> test_bit) & 1
                if b_r == b_r1: match_count += 1
                total += 1
        p = match_count/total if total else 0
        marker = " ★" if abs(p-0.5) > 0.02 else ""
        print(f"    bit {test_bit:2d}: P(same) = {p:.4f}{marker}")


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 200

    disc1_bit_profile(N)
    influences = disc2_influence_map(N)
    disc3_interference(N)
    disc4_coupling(N)

    print("\n" + "="*70)
    print("ИТОГ: ЧТО ОТКРЫЛИ → ЧТО СТРОИМ")
    print("="*70)
    print("""
  Открытие → Инструмент:

  1. Bit profile → Bit-weighted solver (если есть лёгкие биты → целимся в них)
  2. Influence map → Greedy bit-flip solver (начинаем с самых полезных бит)
  3. Interference → Multi-round optimizer (используем биты с positive net)
  4. Coupling → Decoupled solver (если coupling < 0.1 → раунды независимы)
    """)
