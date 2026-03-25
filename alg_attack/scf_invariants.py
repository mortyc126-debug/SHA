#!/usr/bin/env python3
"""
SCF: Каталог инвариантов SHA-256 через 64 раунда.
Ищем что НЕ разрушается в хаосе.
"""
import os, sys, math

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
def hw(x): return bin(x).count('1')
def bit(x,i): return (x>>i)&1

def sha_all_rounds(W16):
    """Return state at every round (0..64)."""
    W = list(W16)
    for i in range(16,64):
        W.append(add32(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    states = []
    s = list(IV)
    states.append(list(s))
    for r in range(64):
        a,b,c,d,e,f,g,h = s
        T1 = add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2 = add32(Sig0(a),Maj(a,b,c))
        s = [add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
        states.append(list(s))
    return states

# ============================================================
# Experiment 1: Parity invariants
# ============================================================
def test_parity_invariants(N):
    """Test if parity of any register or XOR combination is conserved."""
    print("="*70)
    print("EXP 1: PARITY INVARIANTS")
    print("  Test: P(parity(reg[r]) == parity(reg[0])) for each register")
    print("="*70)

    # Track parity correlation: register i, round r
    # parity(x) = hw(x) % 2
    match_counts = [[0]*65 for _ in range(8)]  # 8 regs × 65 rounds

    for _ in range(N):
        W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        states = sha_all_rounds(W)
        par0 = [hw(states[0][i]) % 2 for i in range(8)]
        for r in range(65):
            for i in range(8):
                if hw(states[r][i]) % 2 == par0[i]:
                    match_counts[i][r] += 1

    print(f"\n  P(parity match with r=0) — 0.5 = random")
    regs = ['a','b','c','d','e','f','g','h']
    print(f"{'r':>3} | " + " ".join(f"{regs[i]:>6}" for i in range(8)))
    print("-"*65)
    for r in [0,1,2,3,4,5,10,20,32,64]:
        vals = [match_counts[i][r]/N for i in range(8)]
        markers = ["" if abs(v-0.5)<0.03 else "*" for v in vals]
        print(f"{r:3d} | " + " ".join(f"{vals[i]:5.3f}{markers[i]}" for i in range(8)))

    # XOR combinations: a^e, a^b, etc
    print(f"\n  Cross-register parity: P(parity(X[r]^Y[r]) stable)")
    pairs = [(0,4,'a^e'), (0,1,'a^b'), (4,5,'e^f'), (0,2,'a^c'), (4,6,'e^g')]
    xor_match = {name: [0]*65 for _,_,name in pairs}

    for _ in range(N):
        W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        states = sha_all_rounds(W)
        for i,j,name in pairs:
            par0 = hw(states[0][i] ^ states[0][j]) % 2
            for r in range(65):
                if hw(states[r][i] ^ states[r][j]) % 2 == par0:
                    xor_match[name][r] += 1

    print(f"{'r':>3} | " + " ".join(f"{name:>7}" for _,_,name in pairs))
    print("-"*50)
    for r in [0,1,5,10,32,64]:
        vals = [xor_match[name][r]/N for _,_,name in pairs]
        print(f"{r:3d} | " + " ".join(f"{v:7.3f}" for v in vals))


# ============================================================
# Experiment 2: Hamming weight mod k invariants
# ============================================================
def test_hw_mod_invariants(N):
    print("\n" + "="*70)
    print("EXP 2: HAMMING WEIGHT MOD K")
    print("  Test: bias in HW(a[r]) mod k and HW(a[r])+HW(e[r]) mod k")
    print("="*70)

    for k in [2, 4]:
        # Track distribution of HW(a[r]) mod k
        dist_a = [[0]*k for _ in range(65)]
        dist_e = [[0]*k for _ in range(65)]
        dist_ae = [[0]*k for _ in range(65)]

        for _ in range(N):
            W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
            states = sha_all_rounds(W)
            for r in range(65):
                dist_a[r][hw(states[r][0]) % k] += 1
                dist_e[r][hw(states[r][4]) % k] += 1
                dist_ae[r][(hw(states[r][0]) + hw(states[r][4])) % k] += 1

        print(f"\n  --- mod {k} ---")
        print(f"  HW(a[r]) mod {k} distribution (expected uniform = {1/k:.3f}):")
        print(f"  {'r':>3} | " + " ".join(f"  ={i}" for i in range(k)) + " | chi2")
        for r in [0,1,5,10,32,64]:
            probs = [dist_a[r][i]/N for i in range(k)]
            expected = N/k
            chi2 = sum((dist_a[r][i]-expected)**2/expected for i in range(k))
            sig = "***" if chi2 > 10.83 else ("**" if chi2 > 6.63 else ("*" if chi2 > 3.84 else ""))
            print(f"  {r:3d} | " + " ".join(f"{p:.3f}" for p in probs) + f" | {chi2:6.2f} {sig}")

        print(f"\n  HW(a[r])+HW(e[r]) mod {k}:")
        print(f"  {'r':>3} | " + " ".join(f"  ={i}" for i in range(k)) + " | chi2")
        for r in [0,1,5,10,32,64]:
            probs = [dist_ae[r][i]/N for i in range(k)]
            expected = N/k
            chi2 = sum((dist_ae[r][i]-expected)**2/expected for i in range(k))
            sig = "***" if chi2 > 10.83 else ("**" if chi2 > 6.63 else ("*" if chi2 > 3.84 else ""))
            print(f"  {r:3d} | " + " ".join(f"{p:.3f}" for p in probs) + f" | {chi2:6.2f} {sig}")


# ============================================================
# Experiment 3: Modular invariants (a[r] mod 2^k)
# ============================================================
def test_mod_invariants(N):
    print("\n" + "="*70)
    print("EXP 3: MODULAR INVARIANTS — a[r] mod 2^k")
    print("  Test: does a[r] mod 2^k retain info about a[0] mod 2^k?")
    print("="*70)

    for k in [1, 2, 3, 4]:
        mod = 1 << k
        match_a = [0]*65
        match_e = [0]*65

        for _ in range(N):
            W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
            states = sha_all_rounds(W)
            a0_mod = states[0][0] % mod
            e0_mod = states[0][4] % mod
            for r in range(65):
                if states[r][0] % mod == a0_mod:
                    match_a[r] += 1
                if states[r][4] % mod == e0_mod:
                    match_e[r] += 1

        expected = 1.0/mod
        print(f"\n  mod 2^{k} (expected match = {expected:.4f}):")
        print(f"  {'r':>3} | P(a[r]≡a[0]) P(e[r]≡e[0]) | a_dev   e_dev")
        for r in [0,1,2,3,5,10,20,32,64]:
            pa = match_a[r]/N
            pe = match_e[r]/N
            da = (pa - expected) / max(0.001, (expected*(1-expected)/N)**0.5)
            de = (pe - expected) / max(0.001, (expected*(1-expected)/N)**0.5)
            sig_a = "***" if abs(da)>3 else ("**" if abs(da)>2 else "")
            sig_e = "***" if abs(de)>3 else ("**" if abs(de)>2 else "")
            print(f"  {r:3d} |     {pa:.4f}       {pe:.4f}   | {da:+5.1f}{sig_a:3s} {de:+5.1f}{sig_e:3s}")


# ============================================================
# Experiment 4: Differential invariants
# ============================================================
def test_differential_invariants(N):
    print("\n" + "="*70)
    print("EXP 4: DIFFERENTIAL INVARIANTS")
    print("  Given (W, W⊕δ), what's predictable in δstate after saturation?")
    print("="*70)

    # Track per-bit predictability of differential
    bit_predict = [[0]*32 for _ in range(65)]  # round × bit → count where δa[bit]=0

    # Track HW(δa) - HW(δe) bias
    hw_diff = [[] for _ in range(65)]

    # Track δa ⊕ δe
    da_xor_de_hw = [[] for _ in range(65)]

    for trial in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = list(W1)
        W2[0] ^= (1 << (trial % 32))

        st1 = sha_all_rounds(W1)
        st2 = sha_all_rounds(W2)

        for r in range(65):
            da = st1[r][0] ^ st2[r][0]
            de = st1[r][4] ^ st2[r][4]

            for b in range(32):
                if bit(da, b) == 0:
                    bit_predict[r][b] += 1

            hw_diff[r].append(hw(da) - hw(de))
            da_xor_de_hw[r].append(hw(da ^ de))

    # Report bit predictability
    print(f"\n  Per-bit P(δa[bit]=0) — 0.5 = random, >0.5 = predictable")
    print(f"  {'r':>3} | bit0   bit1   bit8  bit15  bit24  bit31 | avg")
    print("-"*65)
    for r in [0,1,2,3,5,10,20,32,64]:
        bits_to_show = [0,1,8,15,24,31]
        vals = [bit_predict[r][b]/N for b in bits_to_show]
        avg = sum(bit_predict[r][b]/N for b in range(32))/32
        print(f"  {r:3d} | " + " ".join(f"{v:.3f} " for v in vals) + f"| {avg:.3f}")

    # HW(δa) - HW(δe) bias
    print(f"\n  HW(δa) - HW(δe) bias (0 = symmetric)")
    print(f"  {'r':>3} | mean    std    | deviation")
    for r in [0,1,2,5,10,20,32,64]:
        if hw_diff[r]:
            m = sum(hw_diff[r])/len(hw_diff[r])
            v = (sum((x-m)**2 for x in hw_diff[r])/len(hw_diff[r]))**0.5
            sig = abs(m)/(v/len(hw_diff[r])**0.5) if v > 0 else 0
            marker = "***" if sig > 3 else ("**" if sig > 2 else "")
            print(f"  {r:3d} | {m:+6.2f}  {v:5.2f}  | {sig:5.1f}σ {marker}")

    # δa ⊕ δe HW
    print(f"\n  HW(δa ⊕ δe) — expected 16 if independent")
    print(f"  {'r':>3} | mean   std")
    for r in [0,1,2,5,10,20,32,64]:
        if da_xor_de_hw[r]:
            m = sum(da_xor_de_hw[r])/len(da_xor_de_hw[r])
            v = (sum((x-m)**2 for x in da_xor_de_hw[r])/len(da_xor_de_hw[r]))**0.5
            print(f"  {r:3d} | {m:5.2f}  {v:4.2f}")


# ============================================================
# Experiment 5: Shift register verification
# ============================================================
def test_shift_register(N):
    print("\n" + "="*70)
    print("EXP 5: SHIFT REGISTER INVARIANTS (absolute, P=1)")
    print("  Verify: b[r]=a[r-1], c[r]=a[r-2], d[r]=a[r-3]")
    print("         f[r]=e[r-1], g[r]=e[r-2], h[r]=e[r-3]")
    print("="*70)

    violations = 0
    checks = 0
    for _ in range(min(N, 1000)):
        W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        states = sha_all_rounds(W)
        for r in range(3, 65):
            checks += 6
            if states[r][1] != states[r-1][0]: violations += 1  # b=a[-1]
            if states[r][2] != states[r-2][0]: violations += 1  # c=a[-2]
            if states[r][3] != states[r-3][0]: violations += 1  # d=a[-3]
            if states[r][5] != states[r-1][4]: violations += 1  # f=e[-1]
            if states[r][6] != states[r-2][4]: violations += 1  # g=e[-2]
            if states[r][7] != states[r-3][4]: violations += 1  # h=e[-3]

    print(f"\n  {violations} violations / {checks} checks → {'✓ ABSOLUTE INVARIANT' if violations==0 else '✗ BROKEN'}")


# ============================================================
# Experiment 6: Low-bit structure survival
# ============================================================
def test_lowbit_structure(N):
    print("\n" + "="*70)
    print("EXP 6: LOW-BIT STRUCTURE SURVIVAL")
    print("  Do low bits of state retain more structure than high bits?")
    print("="*70)

    # Measure: correlation of bit i of a[r] with W[0] bit i
    corr_by_bit = [[0]*32 for _ in range(65)]

    for trial in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = list(W1)
        b_flip = trial % 32
        W2[0] ^= (1 << b_flip)

        st1 = sha_all_rounds(W1)
        st2 = sha_all_rounds(W2)

        for r in range(65):
            da = st1[r][0] ^ st2[r][0]
            for b in range(32):
                if bit(da, b) == 0:
                    corr_by_bit[r][b] += 1

    print(f"\n  P(δa[bit]=0) by bit position group (0.5=random, 1.0=unaffected)")
    print(f"  {'r':>3} | low(0-7) mid(8-15) mid(16-23) high(24-31) | gap(low-high)")
    print("-"*70)
    for r in [0,1,2,3,4,5,6,10,20,32,64]:
        groups = [
            sum(corr_by_bit[r][b]/N for b in range(0,8))/8,
            sum(corr_by_bit[r][b]/N for b in range(8,16))/8,
            sum(corr_by_bit[r][b]/N for b in range(16,24))/8,
            sum(corr_by_bit[r][b]/N for b in range(24,32))/8,
        ]
        gap = groups[0] - groups[3]
        marker = " ★" if gap > 0.02 else ""
        print(f"  {r:3d} |   {groups[0]:.4f}    {groups[1]:.4f}     {groups[2]:.4f}      {groups[3]:.4f}  | {gap:+.4f}{marker}")


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5000

    print(f"SCF: Каталог инвариантов SHA-256 (N={N})")
    print()

    test_shift_register(N)
    test_parity_invariants(N)
    test_hw_mod_invariants(N)
    test_mod_invariants(N)
    test_differential_invariants(N)
    test_lowbit_structure(N)

    print("\n" + "="*70)
    print("КАТАЛОГ НАЙДЕННЫХ ИНВАРИАНТОВ")
    print("="*70)
    print("""
    АБСОЛЮТНЫЕ (P=1):
      I1. Shift register: b[r]=a[r-1], ..., h[r]=e[r-3]
      I2. Carry-free LSB: cc[bit0] = 0 always (from carry bitwise exp)

    СТАТИСТИЧЕСКИЕ (bias > 0):
      I3. Parity — to be determined from results above
      I4. HW mod k — to be determined
      I5. Modular — to be determined
      I6. Low-bit advantage — to be determined

    NULL RESULTS (no invariant found):
      These are equally important — they confirm chaos is real.
    """)
