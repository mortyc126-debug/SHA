#!/usr/bin/env python3
"""
SCF Задача A: Carleman Linearization SHA-256.

Идея: SHA round = degree 2 (Ch, Maj). Квадратичная система может быть
"поднята" в линейную через введение новых переменных для мономов.

x_new = Ax + B(x⊗x)  →  [x; y=x⊗x] = [A B; C D][x; y]

Если SHA-round квадратична, Carleman-лифт делает её ЛИНЕЙНОЙ.
64 раунда = 64 шага линейной системы = A^64.

Вопрос: какова ЭФФЕКТИВНАЯ размерность Carleman-пространства?
Если мала → линейные методы (eigenvalues, Jordan form) применимы.

Ключевое наблюдение: SHA имеет только 2 нелинейных операции на раунд:
  Ch(e,f,g) = ef ⊕ ēg  (3 квадратичных монома)
  Maj(a,b,c) = ab ⊕ ac ⊕ bc  (3 квадратичных монома)
Всего: 6 новых квадратичных мономов на раунд (не C(256,2)=32640!)
"""
import os, sys

MASK32 = 0xFFFFFFFF

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

# ============================================================
# PART 1: Count effective nonlinear monomials per round
# ============================================================
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
K_CONST = K

def count_monomials():
    """
    In GF(2), one SHA round introduces these quadratic monomials:
    Ch(e,f,g)[i] = e[i]·f[i] ⊕ e[i]·g[i] ⊕ g[i]
      → 2 new monomials per bit: e[i]f[i], e[i]g[i]
      → 64 new monomials for 32 bits of Ch
    Maj(a,b,c)[i] = a[i]·b[i] ⊕ a[i]·c[i] ⊕ b[i]·c[i]
      → 3 new monomials per bit: a[i]b[i], a[i]c[i], b[i]c[i]
      → 96 new monomials for 32 bits

    But: after shift register, b=a[-1], c=a[-2], f=e[-1], g=e[-2]
    So Ch(e,f,g) = Ch(e[r], e[r-1], e[r-2]):
      → monomials: e[r][i]·e[r-1][i], e[r][i]·e[r-2][i]
    Maj(a,b,c) = Maj(a[r], a[r-1], a[r-2]):
      → monomials: a[r][i]·a[r-1][i], a[r][i]·a[r-2][i], a[r-1][i]·a[r-2][i]

    Total per round: 5 types × 32 bits = 160 monomials.
    Over 64 rounds: 160 × 64 = 10240 monomials.
    BUT: many are REUSED (same a[r],e[r] appear in multiple rounds).
    """
    print("="*70)
    print("PART 1: MONOMIAL COUNT PER ROUND (XOR-SHA, GF(2))")
    print("="*70)

    # Track unique monomials across rounds
    # Each monomial = (type, round1, round2, bit)
    # where type = 'ee', 'aa', or 'ae' (cross-term doesn't exist in SHA)

    monomials = set()

    for r in range(64):
        # Ch uses e[r], e[r-1], e[r-2] (through f=e[-1], g=e[-2])
        for i in range(32):
            monomials.add(('e*e', r, r-1, i))  # e[r]·e[r-1]
            monomials.add(('e*e', r, r-2, i))  # e[r]·e[r-2]

        # Maj uses a[r], a[r-1], a[r-2]
        for i in range(32):
            monomials.add(('a*a', r, r-1, i))   # a[r]·a[r-1]
            monomials.add(('a*a', r, r-2, i))   # a[r]·a[r-2]
            monomials.add(('a*a', r-1, r-2, i)) # a[r-1]·a[r-2]

    print(f"  Total unique monomials: {len(monomials)}")
    print(f"  Linear variables: 256 (state) + 512 (message) = 768")
    print(f"  Carleman-lifted dimension: {768 + len(monomials)}")

    # Per-round accumulation
    acc = set()
    print(f"\n  Cumulative monomials by round:")
    print(f"  {'r':>3} | new monomials | cumulative | lifted dim")
    print("  " + "-"*55)
    for r in range(64):
        new = set()
        for i in range(32):
            new.add(('e*e', r, r-1, i))
            new.add(('e*e', r, r-2, i))
            new.add(('a*a', r, r-1, i))
            new.add(('a*a', r, r-2, i))
            new.add(('a*a', r-1, r-2, i))
        before = len(acc)
        acc |= new
        added = len(acc) - before
        if r < 10 or r % 8 == 0 or r == 63:
            print(f"  {r:3d} |         {added:5d} |      {len(acc):5d} | {768+len(acc):6d}")

    return len(monomials)


# ============================================================
# PART 2: Simulate Carleman dimension numerically
# ============================================================
def simulate_carleman_dim(N):
    """
    Run XOR-SHA (no carry) and track which degree-2 monomials
    actually appear in the output as functions of input bits.

    Method: perturb pairs of input bits, check if output changes.
    If flipping bit i AND bit j changes output differently from
    flipping each alone → there's a degree-2 monomial involving i,j.
    """
    print("\n" + "="*70)
    print("PART 2: EFFECTIVE CARLEMAN DIMENSION (numerical)")
    print("  How many input-bit-pairs affect the output?")
    print("="*70)

    def xor_sha(W16):
        W = list(W16)
        for i in range(16,64):
            W.append(sig1(W[i-2])^W[i-7]^sig0(W[i-15])^W[i-16])
        s = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
             0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]
        for r in range(64):
            a,b,c,d,e,f,g,h = s
            T1 = h^Sig1(e)^Ch(e,f,g)^(K_CONST[r] if 'K_CONST' in dir() else 0)^W[r]
            T2 = Sig0(a)^Maj(a,b,c)
            s = [T1^T2,a,b,c,d^T1,e,f,g]
        return s

    # Sample a random base
    W_base = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    H_base = xor_sha(W_base)

    # Count nonlinear interactions between input bits
    # For efficiency, sample random pairs
    n_pairs_tested = 0
    n_nonlinear = 0

    # First-order: single bit flips
    first_order = {}
    for word in range(16):
        for bit_pos in range(32):
            W_flip = list(W_base)
            W_flip[word] ^= (1 << bit_pos)
            H_flip = xor_sha(W_flip)
            dH = tuple(H_base[i] ^ H_flip[i] for i in range(8))
            first_order[(word, bit_pos)] = dH

    # Second-order: pairs of bit flips
    # If δ²H = δH(i,j) ⊕ δH(i) ⊕ δH(j) ≠ 0 → nonlinear interaction
    n_sample = min(N * 10, 5000)
    for trial in range(n_sample):
        w1 = int.from_bytes(os.urandom(1),'big') % 16
        b1 = int.from_bytes(os.urandom(1),'big') % 32
        w2 = int.from_bytes(os.urandom(1),'big') % 16
        b2 = int.from_bytes(os.urandom(1),'big') % 32
        if (w1, b1) >= (w2, b2):
            continue

        W_flip12 = list(W_base)
        W_flip12[w1] ^= (1 << b1)
        W_flip12[w2] ^= (1 << b2)
        H_flip12 = xor_sha(W_flip12)

        dH_12 = tuple(H_base[i] ^ H_flip12[i] for i in range(8))
        dH_1 = first_order[(w1, b1)]
        dH_2 = first_order[(w2, b2)]

        # Second-order differential
        d2H = tuple(dH_12[i] ^ dH_1[i] ^ dH_2[i] for i in range(8))
        n_pairs_tested += 1

        if any(x != 0 for x in d2H):
            n_nonlinear += 1

    frac = n_nonlinear / n_pairs_tested if n_pairs_tested > 0 else 0
    total_pairs = 512 * 511 // 2
    estimated_nonlinear = int(frac * total_pairs)

    print(f"  Pairs tested: {n_pairs_tested}")
    print(f"  Nonlinear interactions: {n_nonlinear} ({frac:.4f})")
    print(f"  Estimated total nonlinear pairs: {estimated_nonlinear} / {total_pairs}")
    print(f"  Estimated Carleman dimension: {512 + estimated_nonlinear}")

    return estimated_nonlinear


# ============================================================
# PART 3: Reduced-round Carleman — does dimension grow or saturate?
# ============================================================
def carleman_by_rounds(N):
    """Measure effective Carleman dimension for R=1,2,...,20 rounds."""
    print("\n" + "="*70)
    print("PART 3: CARLEMAN DIMENSION VS ROUNDS")
    print("="*70)

    def xor_sha_R(W16, R):
        W = list(W16)
        for i in range(16,64):
            W.append(sig1(W[i-2])^W[i-7]^sig0(W[i-15])^W[i-16])
        s = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
             0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]
        for r in range(R):
            a,b,c,d,e,f,g,h = s
            T1 = h^Sig1(e)^Ch(e,f,g)^K_CONST[r]^W[r]
            T2 = Sig0(a)^Maj(a,b,c)
            s = [T1^T2,a,b,c,d^T1,e,f,g]
        return s

    W_base = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    n_sample = min(N, 2000)

    print(f"  {'R':>3} | {'nonlinear%':>10} | {'est_dim':>8} | growth")
    print("  " + "-"*50)

    prev_dim = 512

    for R in [1, 2, 3, 4, 5, 6, 8, 10, 16, 32, 64]:
        H_base = xor_sha_R(W_base, R)

        # First-order
        first_order = {}
        for word in range(16):
            for bp in range(32):
                W_f = list(W_base)
                W_f[word] ^= (1 << bp)
                H_f = xor_sha_R(W_f, R)
                first_order[(word,bp)] = tuple(H_base[i]^H_f[i] for i in range(8))

        # Second-order sample
        n_nl = 0
        n_tested = 0
        for trial in range(n_sample):
            w1 = int.from_bytes(os.urandom(1),'big') % 16
            b1 = int.from_bytes(os.urandom(1),'big') % 32
            w2 = int.from_bytes(os.urandom(1),'big') % 16
            b2 = int.from_bytes(os.urandom(1),'big') % 32
            if (w1,b1) >= (w2,b2): continue

            W_f12 = list(W_base)
            W_f12[w1] ^= (1<<b1)
            W_f12[w2] ^= (1<<b2)
            H_f12 = xor_sha_R(W_f12, R)

            d2 = tuple(
                (H_base[i]^H_f12[i]) ^ first_order[(w1,b1)][i] ^ first_order[(w2,b2)][i]
                for i in range(8))
            n_tested += 1
            if any(x!=0 for x in d2): n_nl += 1

        frac = n_nl/n_tested if n_tested > 0 else 0
        total_pairs = 512*511//2
        est_nl = int(frac * total_pairs)
        est_dim = 512 + est_nl
        growth = est_dim / prev_dim if prev_dim > 0 else 0

        marker = ""
        if R <= 3 and frac < 0.01: marker = " ← LINEAR!"
        elif growth < 1.1: marker = " ← SATURATED"

        print(f"  {R:3d} | {frac:10.4f} | {est_dim:8d} | {growth:.2f}x{marker}")
        prev_dim = est_dim


# ============================================================
# PART 4: Real SHA vs XOR-SHA Carleman comparison
# ============================================================
def carleman_real_vs_xor(N):
    print("\n" + "="*70)
    print("PART 4: REAL SHA vs XOR-SHA — Carleman dimension")
    print("="*70)

    def sha_R(W16, R, use_xor=False):
        if use_xor:
            W = list(W16)
            for i in range(16,64):
                W.append(sig1(W[i-2])^W[i-7]^sig0(W[i-15])^W[i-16])
        else:
            W = list(W16)
            for i in range(16,64):
                W.append(add32(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))

        s = list([0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
                  0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19])
        for r in range(R):
            a,b,c,d,e,f,g,h = s
            if use_xor:
                T1 = h^Sig1(e)^Ch(e,f,g)^K_CONST[r]^W[r]
                T2 = Sig0(a)^Maj(a,b,c)
                s = [T1^T2,a,b,c,d^T1,e,f,g]
            else:
                T1 = add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
                T2 = add32(Sig0(a),Maj(a,b,c))
                s = [add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
        return s

    W_base = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    n_sample = min(N, 1500)

    print(f"  {'R':>3} | {'XOR nonlin%':>11} {'Real nonlin%':>12} | {'XOR dim':>8} {'Real dim':>9}")
    print("  " + "-"*60)

    for R in [1, 2, 3, 4, 6, 10, 20, 64]:
        results = {}
        for mode_name, use_xor in [("xor", True), ("real", False)]:
            H_base = sha_R(W_base, R, use_xor)
            fo = {}
            for w in range(16):
                for bp in range(32):
                    Wf = list(W_base); Wf[w]^=(1<<bp)
                    Hf = sha_R(Wf, R, use_xor)
                    fo[(w,bp)] = tuple(H_base[i]^Hf[i] for i in range(8))

            n_nl, n_t = 0, 0
            for trial in range(n_sample):
                w1=int.from_bytes(os.urandom(1),'big')%16
                b1=int.from_bytes(os.urandom(1),'big')%32
                w2=int.from_bytes(os.urandom(1),'big')%16
                b2=int.from_bytes(os.urandom(1),'big')%32
                if (w1,b1)>=(w2,b2): continue
                Wf=list(W_base); Wf[w1]^=(1<<b1); Wf[w2]^=(1<<b2)
                Hf=sha_R(Wf, R, use_xor)
                d2=tuple((H_base[i]^Hf[i])^fo[(w1,b1)][i]^fo[(w2,b2)][i] for i in range(8))
                n_t += 1
                if any(x!=0 for x in d2): n_nl += 1

            frac = n_nl/n_t if n_t>0 else 0
            est = 512 + int(frac * 512*511//2)
            results[mode_name] = (frac, est)

        xf, xd = results["xor"]
        rf, rd = results["real"]
        gap = " ★" if abs(xf-rf) > 0.05 else ""
        print(f"  {R:3d} | {xf:11.4f} {rf:12.4f} | {xd:8d} {rd:9d}{gap}")


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 500

    print("="*70)
    print("SCF ЗАДАЧА A: CARLEMAN LINEARIZATION OF SHA-256")
    print("="*70)

    n_mono = count_monomials()
    simulate_carleman_dim(N)
    carleman_by_rounds(N)
    carleman_real_vs_xor(N)

    print("\n" + "="*70)
    print("ИТОГ ЗАДАЧИ A")
    print("="*70)
    print(f"""
  Теоретически:
    SHA round: degree 2 (Ch + Maj)
    Уникальные мономы: {n_mono} (за 64 раунда)
    Carleman dimension: 768 + {n_mono} = {768+n_mono}

  Если Carleman-dim < 2^20 → линейная система решаема!
  Если Carleman-dim насыщается → SHA "компактна" в lifted space.
  Если Carleman-dim растёт экспоненциально → подъём не помогает.
    """)
