#!/usr/bin/env python3
"""
SCF: SHADOW ALGEBRA — алгебра дифференциальных паттернов SHA-256.
Построена С НУЛЯ из структуры SHA.

Элемент: 8-бит вектор P = (δa≠0, δb≠0, ..., δh≠0)
  P[i] = 1 если регистр i имеет ненулевой дифференциал
  P[i] = 0 если δreg[i] = 0

Shift register constraints:
  P_b[r+1] = P_a[r]  (b←a)
  P_c[r+1] = P_b[r] = P_a[r-1]
  P_d[r+1] = P_c[r] = P_a[r-2]
  P_f[r+1] = P_e[r]
  P_g[r+1] = P_f[r] = P_e[r-1]
  P_h[r+1] = P_g[r] = P_e[r-2]

Nonlinear nodes:
  P_a[r+1] = 1 if δ(T1+T2) ≠ 0
  P_e[r+1] = 1 if δ(d+T1) ≠ 0

T2 depends on a,b,c → P_a, P_b, P_c
T1 depends on h,e,f,g,W → P_h, P_e, P_f, P_g, P_W
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

def expand_real(W16):
    W=list(W16)
    for i in range(16,64):
        W.append(add32(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    return W

def sha_states(W16):
    W=expand_real(W16); s=list(IV); states=[list(s)]
    for r in range(64):
        a,b,c,d,e,f,g,h=s
        T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add32(Sig0(a),Maj(a,b,c))
        s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
        states.append(list(s))
    return states

def pattern(st1, st2):
    """8-bit shadow: which registers differ?"""
    return sum((1 if st1[i]!=st2[i] else 0) << i for i in range(8))

def pattern_str(p):
    regs = "abcdefgh"
    return "".join(regs[i] if (p>>i)&1 else "." for i in range(8))


# ============================================================
# EXP 1: Measure empirical shadow trajectories
# ============================================================
def exp1_shadow_trajectories(N):
    print("="*70)
    print("EXP 1: EMPIRICAL SHADOW TRAJECTORIES")
    print("  Pattern = 8 bits: which registers have δ≠0")
    print("="*70)

    from collections import Counter
    pattern_at_round = [Counter() for _ in range(65)]

    for trial in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = list(W1); W2[0] ^= 0x80000000

        st1 = sha_states(W1)
        st2 = sha_states(W2)

        for r in range(65):
            p = pattern(st1[r], st2[r])
            pattern_at_round[r][p] += 1

    print(f"\n  {'r':>3} | top patterns (showing most common)")
    print("  " + "-"*60)
    for r in [0,1,2,3,4,5,6,10,16,17,20,30,40,50,60,64]:
        top = pattern_at_round[r].most_common(3)
        parts = []
        for p, count in top:
            parts.append(f"{pattern_str(p)}({count})")
        n_unique = len(pattern_at_round[r])
        print(f"  {r:3d} | {' '.join(parts):40s} | {n_unique} unique")


# ============================================================
# EXP 2: Shadow transition rules (empirical)
# ============================================================
def exp2_transitions(N):
    print("\n" + "="*70)
    print("EXP 2: SHADOW TRANSITION RULES")
    print("  P(pattern[r+1] | pattern[r]) — what follows what?")
    print("="*70)

    from collections import Counter, defaultdict
    transitions = defaultdict(Counter)  # {pattern_r: Counter(pattern_r+1)}

    for trial in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = list(W1); W2[0] ^= int.from_bytes(os.urandom(4),'big')
        if W2[0]==W1[0]: W2[0]^=1

        st1 = sha_states(W1)
        st2 = sha_states(W2)

        for r in range(20, 60):  # Focus on gap region
            p_r = pattern(st1[r], st2[r])
            p_r1 = pattern(st1[r+1], st2[r+1])
            transitions[p_r][p_r1] += 1

    # Analyze: from pattern 0x00 (all zero), what happens?
    print(f"\n  From all-zero pattern (........):")
    if 0 in transitions:
        total = sum(transitions[0].values())
        for p, count in transitions[0].most_common(5):
            print(f"    → {pattern_str(p)}: {count}/{total} = {count/total:.3f}")
    else:
        print(f"    Never observed")

    # From pattern 0xFF (all nonzero)
    print(f"\n  From all-nonzero (abcdefgh):")
    if 0xFF in transitions:
        total = sum(transitions[0xFF].values())
        for p, count in transitions[0xFF].most_common(5):
            print(f"    → {pattern_str(p)}: {count}/{total} = {count/total:.3f}")

    # KEY: which patterns can transition TO all-zero?
    print(f"\n  Which patterns transition TO all-zero?")
    to_zero = []
    for p_from, counts in transitions.items():
        if 0 in counts:
            total = sum(counts.values())
            prob = counts[0] / total
            to_zero.append((p_from, prob, counts[0], total))

    if to_zero:
        to_zero.sort(key=lambda x: -x[1])
        for p, prob, cnt, tot in to_zero[:10]:
            print(f"    {pattern_str(p)} → ........: P={prob:.4f} ({cnt}/{tot})")
    else:
        print(f"    NONE — all-zero is unreachable from gap patterns!")

    # How many unique patterns exist in the gap?
    all_patterns = set()
    for p_from in transitions:
        all_patterns.add(p_from)
        for p_to in transitions[p_from]:
            all_patterns.add(p_to)
    print(f"\n  Unique patterns observed in gap (r=20..60): {len(all_patterns)} / 256")

    # Pattern with most zeros
    print(f"\n  Patterns with most zero registers:")
    for p in sorted(all_patterns, key=lambda x: bin(x).count('1')):
        n_zero = 8 - bin(p).count('1')
        if n_zero >= 3:
            total_obs = sum(transitions[p].values()) if p in transitions else 0
            print(f"    {pattern_str(p)} ({n_zero} zeros): observed {total_obs} times")


# ============================================================
# EXP 3: Theoretical shadow transition (deterministic part)
# ============================================================
def exp3_deterministic():
    print("\n" + "="*70)
    print("EXP 3: DETERMINISTIC SHADOW RULES")
    print("  From shift register alone (no probabilistic part)")
    print("="*70)

    # Shift register gives us:
    # b[r+1] = a[r] → P_b' = P_a
    # c[r+1] = b[r] → P_c' = P_b = P_a[-1]
    # d[r+1] = c[r] → P_d' = P_c = P_a[-2]
    # f[r+1] = e[r] → P_f' = P_e
    # g[r+1] = f[r] → P_g' = P_f = P_e[-1]
    # h[r+1] = g[r] → P_h' = P_g = P_e[-2]
    #
    # Nonlinear nodes:
    # P_a' = ? (depends on T1+T2, always 1 if ANY input ≠ 0)
    # P_e' = ? (depends on d+T1, always 1 if P_d=1 OR P_T1≠0)
    #
    # Conservative rule (worst case):
    # P_a' = 1 if P_a OR P_b OR P_c OR P_e OR P_f OR P_g OR P_h OR P_W
    # P_e' = 1 if P_d OR P_e OR P_f OR P_g OR P_h OR P_W
    #
    # Simplification: once any bit is 1, it spreads everywhere.
    # QUESTION: can we have P_a'=0? Only if δT1=0 AND δT2=0.
    # δT2=0 requires δa=δb=δc=0 (from Maj)
    # δT1=0 requires δh=δe=δf=δg=0 AND δW=0 AND δ(carry)=0
    #   → basically ALL registers zero AND δW=0

    print("""
  Deterministic rules:
    P_b' = P_a          (shift)
    P_c' = P_b          (shift)
    P_d' = P_c          (shift)
    P_f' = P_e          (shift)
    P_g' = P_f          (shift)
    P_h' = P_g          (shift)

  Nonlinear nodes (probabilistic):
    P_a' = 0 ONLY IF: P_a=P_b=P_c=0 (δT2=0)
                  AND: P_h=P_e=P_f=P_g=0 AND P_W=0 (δT1=0)
           → P_a' = 0 requires FULL ZERO + δW=0

    P_e' = 0 ONLY IF: P_d=0 (δd=0)
                  AND: P_h=P_e=P_f=P_g=0 AND P_W=0 (δT1=0)
           → P_e' = 0 requires δd=0 AND most registers zero + δW=0

  CONCLUSION:
    In the shadow algebra, once pattern becomes all-1s (0xFF),
    it STAYS all-1s unless δW[r]=0 for that round.

    The ONLY way to create zeros is through δW[r]=0.
    This links shadow algebra directly to schedule kernel!
    """)

    # Compute: for each 8-bit pattern, what's the deterministic next?
    # (assuming δW≠0, which is true in the gap)
    print("  Deterministic transition (assuming δW≠0):")
    print("  Pattern → deterministic part of next pattern")
    print("  " + "-"*40)

    for p in range(256):
        Pa = (p>>0)&1; Pb = (p>>1)&1; Pc = (p>>2)&1; Pd = (p>>3)&1
        Pe = (p>>4)&1; Pf = (p>>5)&1; Pg = (p>>6)&1; Ph = (p>>7)&1

        # Shift register (deterministic)
        Pb_new = Pa
        Pc_new = Pb
        Pd_new = Pc
        Pf_new = Pe
        Pg_new = Pf
        Ph_new = Pg

        # Nonlinear: with δW≠0, T1 has δW contribution → δT1≠0 almost surely
        # → Pa_new = 1 (from T1+T2, T1≠0)
        # → Pe_new = 1 (from d+T1, T1≠0)
        Pa_new = 1  # Almost always 1 when δW≠0
        Pe_new = 1  # Almost always 1 when δW≠0

        p_new = Pa_new | (Pb_new<<1) | (Pc_new<<2) | (Pd_new<<3) | \
                (Pe_new<<4) | (Pf_new<<5) | (Pg_new<<6) | (Ph_new<<7)

        # Only show interesting patterns
        n_zero_in = 8 - bin(p).count('1')
        n_zero_out = 8 - bin(p_new).count('1')
        if n_zero_in >= 2 or n_zero_out >= 2:
            print(f"    {pattern_str(p)} → {pattern_str(p_new)}"
                  f"  (zeros: {n_zero_in} → {n_zero_out})")

    # KEY FINDING: count rounds to reach all-1s from any start
    print(f"\n  Rounds to all-1s from each pattern (δW≠0):")
    for p in range(256):
        if p == 0xFF: continue
        current = p
        for step in range(10):
            Pa = (current>>0)&1; Pb = (current>>1)&1; Pc = (current>>2)&1
            Pe = (current>>4)&1; Pf = (current>>5)&1; Pg = (current>>6)&1
            Pb_n=Pa; Pc_n=Pb; Pd_n=Pc; Pf_n=Pe; Pg_n=Pf; Ph_n=Pg
            current = 1|(Pb_n<<1)|(Pc_n<<2)|(Pd_n<<3)|0x10|(Pf_n<<5)|(Pg_n<<6)|(Ph_n<<7)
            if current == 0xFF:
                n_z = 8-bin(p).count('1')
                if n_z >= 4:
                    print(f"    {pattern_str(p)} → all-1s in {step+1} rounds")
                break


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 500

    print("="*70)
    print("SCF: SHADOW ALGEBRA — новая алгебра из первых принципов")
    print("="*70)

    exp3_deterministic()
    exp1_shadow_trajectories(min(N, 500))
    exp2_transitions(min(N, 1000))
