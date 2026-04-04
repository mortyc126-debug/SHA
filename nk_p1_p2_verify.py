#!/usr/bin/env python3
"""
Nova Cryptarithmetica — P1 + P2 verification

P1: Shield potential S(r) for post-barrier rounds (r=17..31)
    How much of the △ from δa gets ABSORBED vs PASSED at each round
    AFTER the Wang cascade ends?

P2: Minimum cost through r=17..20 in GPK terms
    How many P-positions (undetermined carry bits) exist at barrier rounds?
    What's the minimum over different δW[0] choices?
"""

import random
from collections import defaultdict

MASK32 = 0xFFFFFFFF

K = [0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2]
IV = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&MASK32
def sig0(x): return rotr(x,2)^rotr(x,13)^rotr(x,22)
def sig1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def ssig0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def ssig1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def ch(e,f,g): return (e&f)^(~e&g)&MASK32
def maj(a,b,c): return (a&b)^(a&c)^(b&c)
def hw(x): return bin(x&MASK32).count('1')

def run_sha(M, R):
    W = list(M)+[0]*(64-len(M))
    for i in range(16,64): W[i]=(ssig1(W[i-2])+W[i-7]+ssig0(W[i-15])+W[i-16])&MASK32
    a,b,c,d,e,f,g,h = IV
    states = [{'a':a,'b':b,'c':c,'d':d,'e':e,'f':f,'g':g,'h':h}]
    T1s, T2s = [], []
    for r in range(R):
        T1=(h+sig1(e)+ch(e,f,g)+K[r]+W[r])&MASK32
        T2=(sig0(a)+maj(a,b,c))&MASK32
        T1s.append(T1); T2s.append(T2)
        h,g,f=g,f,e; e=(d+T1)&MASK32; d,c,b=c,b,a; a=(T1+T2)&MASK32
        states.append({'a':a,'b':b,'c':c,'d':d,'e':e,'f':f,'g':g,'h':h})
    return states, T1s, T2s, W

def wang_cascade(W_base, dW0, R):
    W1 = list(W_base)
    DW = [0]*16; DW[0] = dW0
    for step in range(min(R-1,15)):
        wi = step+1
        W2 = [(W1[i]+DW[i])&MASK32 for i in range(16)]
        s1,_,_,_ = run_sha(W1, step+2)
        s2,_,_,_ = run_sha(W2, step+2)
        de = (s2[-1]['e']-s1[-1]['e'])&MASK32
        DW[wi] = (-de)&MASK32
    W2 = [(W1[i]+DW[i])&MASK32 for i in range(16)]
    return W1, W2, DW

def gpk_at_bit(a_bit, b_bit):
    if a_bit==1 and b_bit==1: return 'G'
    if a_bit==0 and b_bit==0: return 'K'
    return 'P'

def count_gpk(a_val, b_val):
    """Count G, P, K symbols and return P-count (= undetermined carry bits)."""
    g=p=k=0
    for bit in range(32):
        ab=(a_val>>bit)&1; bb=(b_val>>bit)&1
        if ab==1 and bb==1: g+=1
        elif ab==0 and bb==0: k+=1
        else: p+=1
    return g, p, k


# ================================================================
# P1: Shield potential for POST-BARRIER rounds
#
# Extend Wang cascade to 16 rounds, then let SHA run freely
# for rounds 17-31 (no more adaptive DW).
# Measure: at each post-barrier round, how many bit-differences
# in δstate get ABSORBED vs PROPAGATED?
# ================================================================

print("=" * 70)
print("P1: Shield potential for post-barrier rounds (r=17..31)")
print("=" * 70)

R_TOTAL = 32
N_PAIRS = 300
random.seed(42)

# Per-round statistics
p1_stats = defaultdict(lambda: {
    'hw_da': [], 'hw_de': [],
    'hw_da_xor': [], 'hw_de_xor': [],
    'new_delta_bits': [],    # bits that were 0 in δstate[r] but 1 in δstate[r+1]
    'killed_delta_bits': [], # bits that were 1 in δstate[r] but 0 in δstate[r+1]
    'survive_rate_e': [],    # P(δe[r] bit survived to δe[r+1])
    'birth_rate_e': [],      # P(new δe bit born at r+1)
})

for trial in range(N_PAIRS):
    W_base = [random.randint(0,MASK32) for _ in range(16)]
    W1, W2, DW = wang_cascade(W_base, 1, min(R_TOTAL, 16))

    s1,t1_1,t2_1,Wexp1 = run_sha(W1, R_TOTAL)
    s2,t1_2,t2_2,Wexp2 = run_sha(W2, R_TOTAL)

    for r in range(R_TOTAL):
        st1 = s1[r+1]; st2 = s2[r+1]  # state AFTER round r
        da = (st2['a']-st1['a'])&MASK32
        de = (st2['e']-st1['e'])&MASK32
        da_xor = st1['a']^st2['a']
        de_xor = st1['e']^st2['e']

        p1_stats[r]['hw_da'].append(hw(da))
        p1_stats[r]['hw_de'].append(hw(de))
        p1_stats[r]['hw_da_xor'].append(hw(da_xor))
        p1_stats[r]['hw_de_xor'].append(hw(de_xor))

        if r > 0:
            prev_st1 = s1[r]; prev_st2 = s2[r]
            # Full state XOR diff
            prev_diff = 0
            curr_diff = 0
            for reg in ['a','b','c','d','e','f','g','h']:
                prev_diff |= prev_st1[reg] ^ prev_st2[reg]  # any bit differs
                curr_diff |= st1[reg] ^ st2[reg]

            # For e-register specifically
            prev_de = prev_st1['e'] ^ prev_st2['e']
            curr_de = st1['e'] ^ st2['e']

            survived = hw(prev_de & curr_de)  # bits that were diff AND still diff
            born = hw(curr_de & ~prev_de)      # bits that are new diffs
            killed = hw(prev_de & ~curr_de)    # bits that were diff but now same

            total_prev = hw(prev_de)
            p1_stats[r]['survive_rate_e'].append(survived / max(total_prev, 1))
            p1_stats[r]['birth_rate_e'].append(born)
            p1_stats[r]['killed_delta_bits'].append(killed)
            p1_stats[r]['new_delta_bits'].append(born)

print(f"\nN={N_PAIRS} Wang pairs, R={R_TOTAL}")
print(f"\n{'r':>3} | {'HW(δa)':>6} {'HW(δe)':>6} | {'δe_xor':>6} | {'survive':>7} {'born':>5} {'killed':>6} | shield_e")

for r in range(R_TOTAL):
    da = sum(p1_stats[r]['hw_da'])/N_PAIRS
    de = sum(p1_stats[r]['hw_de'])/N_PAIRS
    de_xor = sum(p1_stats[r]['hw_de_xor'])/N_PAIRS

    if r > 0 and p1_stats[r]['survive_rate_e']:
        surv = sum(p1_stats[r]['survive_rate_e'])/N_PAIRS
        born = sum(p1_stats[r]['birth_rate_e'])/N_PAIRS
        killed = sum(p1_stats[r]['killed_delta_bits'])/N_PAIRS

        # Shield = fraction of prev-round δe bits that got killed
        prev_de = sum(p1_stats[r-1]['hw_de_xor'])/N_PAIRS
        shield = killed / max(prev_de, 0.01)

        zone = "CASCADE" if r <= 16 and de < 0.5 else ("BARRIER" if r == 17 else "POST")
        if r <= 1: zone = "INIT"
        print(f"{r:>3} | {da:>6.1f} {de:>6.1f} | {de_xor:>6.1f} | {surv:>7.3f} {born:>5.1f} {killed:>6.1f} | {shield:>7.3f}  [{zone}]")
    else:
        print(f"{r:>3} | {da:>6.1f} {de:>6.1f} | {de_xor:>6.1f} |       -     -      - |       -")


# ================================================================
# P2: Minimum P-positions at barrier rounds (r=17..20)
#
# For each δW[0], the Wang cascade produces specific δW[1..15].
# Schedule then determines δW[16..19].
# At round 17: e_new = d + T1, and T1 includes W[16].
# GPK of (d, T1) tells us how many bits are P (undetermined carry).
#
# Lower P-count = cheaper path through the barrier.
# ================================================================

print()
print("=" * 70)
print("P2: P-positions (carry cost) at barrier rounds r=17..20")
print("=" * 70)

N_P2 = 500
random.seed(123)

# For each trial: compute GPK of the e-computation at rounds 17-20
# and count P-positions (undetermined carry = the "cost" in bits)

p2_stats = defaultdict(lambda: {'p_count': [], 'g_count': [], 'k_count': []})

for trial in range(N_P2):
    W_base = [random.randint(0,MASK32) for _ in range(16)]
    W1, W2, DW = wang_cascade(W_base, 1, 20)

    s1,t1_1,t2_1,Wexp1 = run_sha(W1, 21)
    s2,t1_2,t2_2,Wexp2 = run_sha(W2, 21)

    for r in range(17, 21):
        if r >= len(t1_1): break

        # e_new = d + T1. GPK of this addition.
        d1 = s1[r]['d']; d2 = s2[r]['d']
        T1_v1 = t1_1[r]; T1_v2 = t1_2[r]

        # GPK for message 1
        g1, p1, k1 = count_gpk(d1, T1_v1)
        # GPK for message 2
        g2, p2, k2 = count_gpk(d2, T1_v2)

        # The DIFFERENTIAL GPK: how many positions have DIFFERENT GPK?
        diff_gpk = 0
        p_in_diff = 0  # P-positions among differing GPK
        for bit in range(32):
            gpk1 = gpk_at_bit((d1>>bit)&1, (T1_v1>>bit)&1)
            gpk2 = gpk_at_bit((d2>>bit)&1, (T1_v2>>bit)&1)
            if gpk1 != gpk2:
                diff_gpk += 1
                if gpk1 == 'P' or gpk2 == 'P':
                    p_in_diff += 1

        p2_stats[r]['p_count'].append(p1)  # P-positions in msg1
        p2_stats[r]['g_count'].append(g1)
        p2_stats[r]['k_count'].append(k1)

print(f"\nN={N_P2} pairs")
print(f"\n{'r':>3} | {'E[G]':>5} {'E[P]':>5} {'E[K]':>5} | {'min_P':>5} {'max_P':>5} | {'P<8':>4} {'P<12':>5}")

for r in range(17, 21):
    if not p2_stats[r]['p_count']: continue
    ps = p2_stats[r]['p_count']
    gs = p2_stats[r]['g_count']
    ks = p2_stats[r]['k_count']

    ep = sum(ps)/len(ps)
    eg = sum(gs)/len(gs)
    ek = sum(ks)/len(ks)
    minp = min(ps)
    maxp = max(ps)
    p_lt8 = sum(1 for x in ps if x < 8)
    p_lt12 = sum(1 for x in ps if x < 12)

    print(f"{r:>3} | {eg:>5.1f} {ep:>5.1f} {ek:>5.1f} | {minp:>5} {maxp:>5} | {p_lt8:>4} {p_lt12:>5}")


# ================================================================
# P2 extension: Does δW[0] choice affect P-count at barrier?
# ================================================================

print()
print("=" * 70)
print("P2 ext: P-count at r=17 vs δW[0] bit position")
print("=" * 70)

N_P2B = 200
random.seed(456)

bit_results = []

for bit in [0, 6, 15, 16, 25, 31]:  # sample of positions
    dW0 = 1 << bit
    p_counts_r17 = []

    for trial in range(N_P2B):
        W_base = [random.randint(0,MASK32) for _ in range(16)]
        W1, W2, DW = wang_cascade(W_base, dW0, 18)
        s1,t1_1,_,_ = run_sha(W1, 18)
        s2,t1_2,_,_ = run_sha(W2, 18)

        d1 = s1[17]['d']; T1_v1 = t1_1[17]
        _, p1, _ = count_gpk(d1, T1_v1)
        p_counts_r17.append(p1)

    avg_p = sum(p_counts_r17)/N_P2B
    min_p = min(p_counts_r17)
    max_p = max(p_counts_r17)
    bit_results.append((bit, avg_p, min_p, max_p))

print(f"\nN={N_P2B} pairs per bit")
print(f"\n{'bit':>3} | {'E[P]':>5} | {'min':>4} {'max':>4} | note")
for bit, avg, mn, mx in sorted(bit_results, key=lambda x: x[1]):
    note = "← Wang standard" if bit == 0 else ""
    if avg == min(r[1] for r in bit_results): note += " ★ BEST"
    print(f"{bit:>3} | {avg:>5.1f} | {mn:>4} {mx:>4} | {note}")


# ================================================================
# P2 deeper: What determines P-count?
# Hypothesis: P-count at r=17 depends on HW(δW[16])
# because δW[16] = σ₁(δW[14]) + δW[9] + σ₀(δW[1]) + δW[0]
# ================================================================

print()
print("=" * 70)
print("P2 deep: P-count at r=17 vs HW(δW[16])")
print("=" * 70)

N_P2C = 500
random.seed(789)

hw_dw16_vs_p = defaultdict(list)

for trial in range(N_P2C):
    W_base = [random.randint(0,MASK32) for _ in range(16)]
    W1, W2, DW = wang_cascade(W_base, 1, 18)

    # Compute δW[16] from schedule
    Wexp1 = list(W1)+[0]*48
    Wexp2 = list(W2)+[0]*48
    for i in range(16,64):
        Wexp1[i]=(ssig1(Wexp1[i-2])+Wexp1[i-7]+ssig0(Wexp1[i-15])+Wexp1[i-16])&MASK32
        Wexp2[i]=(ssig1(Wexp2[i-2])+Wexp2[i-7]+ssig0(Wexp2[i-15])+Wexp2[i-16])&MASK32

    dW16 = Wexp1[16] ^ Wexp2[16]
    hw_dw16 = hw(dW16)

    s1,t1_1,_,_ = run_sha(W1, 18)
    d1 = s1[17]['d']; T1_v1 = t1_1[17]
    _, p_count, _ = count_gpk(d1, T1_v1)

    bucket = (hw_dw16 // 4) * 4
    hw_dw16_vs_p[bucket].append(p_count)

print(f"\nN={N_P2C}")
print(f"\n{'HW(δW16)':>9} | {'E[P]':>5} | {'min':>4} {'max':>4} | {'n':>4}")
for bucket in sorted(hw_dw16_vs_p.keys()):
    ps = hw_dw16_vs_p[bucket]
    if ps:
        print(f"{bucket:>5}-{bucket+3:<3} | {sum(ps)/len(ps):>5.1f} | {min(ps):>4} {max(ps):>4} | {len(ps):>4}")


# ================================================================
# SUMMARY
# ================================================================

print()
print("=" * 70)
print("SUMMARY")
print("=" * 70)

print("""
P1 (Shield potential post-barrier):
  Rounds 1-16: δe=0 (Wang cascade), shield = N/A (no δe to shield)
  Round 17+: δe ≠ 0, shield analysis measures bit survival/birth/death

P2 (P-positions at barrier):
  GPK analysis of d+T1 at rounds 17-20 gives P-count =
  number of carry bits that must be "guessed" for the path.

  Lower P-count = cheaper differential path.

Next steps:
  - Correlate P-count with actual collision probability
  - Find W configurations that minimize P at barrier
  - Connect to Λ (total path cost) formula from NK §3
""")
