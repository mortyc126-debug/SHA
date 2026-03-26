#!/usr/bin/env python3
"""
SCF: Matching Condition Dimension Analysis
Вопрос: δπ_⊕ = δπ_Ψ — это 256 бит ограничений или меньше?

Если matching space имеет rank < 256 → есть свободные направления.
Свободные направления = пространство, где collision проще.
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

def expand_schedule(W16):
    W=list(W16)
    for i in range(16,64):
        W.append(add32(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    return W

def expand_schedule_xor(W16):
    W=list(W16)
    for i in range(16,64):
        W.append(sig1(W[i-2])^W[i-7]^sig0(W[i-15])^W[i-16])
    return W

def sha_real(W16):
    W=expand_schedule(W16); s=list(IV)
    for r in range(64):
        a,b,c,d,e,f,g,h=s
        T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add32(Sig0(a),Maj(a,b,c))
        s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
    return s

def sha_xor(W16):
    W=expand_schedule_xor(W16); s=list(IV)
    for r in range(64):
        a,b,c,d,e,f,g,h=s
        T1=h^Sig1(e)^Ch(e,f,g)^K[r]^W[r]
        T2=Sig0(a)^Maj(a,b,c)
        s=[T1^T2,a,b,c,d^T1,e,f,g]
    return s

def state_to_bits(state):
    """Convert 8×32-bit state to 256-bit vector (list of 0/1)."""
    bits = []
    for word in state:
        for i in range(32):
            bits.append((word >> i) & 1)
    return bits

def bits_to_state(bits):
    state = []
    for w in range(8):
        val = 0
        for i in range(32):
            val |= bits[w*32+i] << i
        state.append(val)
    return state

# ============================================================
# GF(2) rank computation (Gaussian elimination)
# ============================================================
def gf2_rank(matrix):
    """Compute rank of binary matrix over GF(2).
    matrix = list of rows, each row is a list of 0/1."""
    if not matrix:
        return 0
    m = [list(row) for row in matrix]
    nrows = len(m)
    ncols = len(m[0])
    rank = 0
    for col in range(ncols):
        # Find pivot
        pivot = -1
        for row in range(rank, nrows):
            if m[row][col] == 1:
                pivot = row
                break
        if pivot == -1:
            continue
        # Swap
        m[rank], m[pivot] = m[pivot], m[rank]
        # Eliminate
        for row in range(nrows):
            if row != rank and m[row][col] == 1:
                for c in range(ncols):
                    m[row][c] ^= m[rank][c]
        rank += 1
    return rank


# ============================================================
# EXP 1: Jacobian rank of the matching map M(W) = π_⊕(W) ⊕ π_Ψ(W)
# ============================================================
def exp1_matching_jacobian(N_samples=50):
    """
    M(W) = F_xor(W) ⊕ (F_real(W) ⊕ F_xor(W)) = F_real(W)... wait.

    Actually: matching condition is about DIFFERENTIALS.
    For collision F(x) = F(y), with x = W, y = W⊕δ:
      δF_real = 0
      ⟺ δF_xor ⊕ δΨ = 0
      ⟺ δF_xor = δΨ

    So we need the Jacobian of the map δ → (δF_xor(W,δ), δΨ(W,δ))
    evaluated at a specific W.

    The matching constraint is δF_xor = δΨ, i.e., δF_xor ⊕ δΨ = 0.
    This is a 256-bit constraint on δ (512 bits).

    Effective dimension of matching = rank of J_{matching}.
    If rank < 256 → there are null directions → easier collision.
    """
    print("="*70)
    print("EXP 1: JACOBIAN RANK OF MATCHING CONDITION")
    print("  M(δ) = δF_xor(W,δ) ⊕ δΨ(W,δ)")
    print("  Collision ⟺ M(δ) = 0")
    print("  rank(J_M) = effective number of constraints")
    print("="*70)

    ranks = []

    for sample in range(N_samples):
        W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

        # Build Jacobian: flip each of 512 input bits, measure 256-bit output
        # But 512×256 is big. Sample: use only W[0] (32 bits) for speed.
        # Jacobian of M restricted to δW[0].

        jacobian_rows = []

        for word_idx in range(min(16, 4)):  # First 4 words
            for bit_idx in range(32):
                delta_W = list(W)
                delta_W[word_idx] ^= (1 << bit_idx)

                # δF_xor
                xor1 = sha_xor(W)
                xor2 = sha_xor(delta_W)
                d_xor = [xor1[i] ^ xor2[i] for i in range(8)]

                # δΨ = δF_real ⊕ δF_xor
                real1 = sha_real(W)
                real2 = sha_real(delta_W)
                d_real = [real1[i] ^ real2[i] for i in range(8)]
                d_psi = [d_real[i] ^ d_xor[i] for i in range(8)]

                # M(δ) = δF_xor ⊕ δΨ = δF_real (trivially!)
                # Wait — this is just δF_real. The rank of J_M = rank of J_F.
                # That's always 256 (full rank) — we know this.

                # The REAL question: what is the rank of J_{δΨ} alone?
                # And: rank of the map δ → (δF_xor, δΨ) as a 512-dim output?
                # If δF_xor and δΨ are correlated, the joint rank < 512.

                row_xor = state_to_bits(d_xor)
                row_psi = state_to_bits(d_psi)
                # Joint: 512-bit row
                joint_row = row_xor + row_psi

                jacobian_rows.append(joint_row)

        # Rank of joint Jacobian (input: 128 bits from W[0..3], output: 512 bits)
        n_input = len(jacobian_rows)
        rank_joint = gf2_rank(jacobian_rows)

        # Rank of δF_xor alone
        xor_rows = [row[:256] for row in jacobian_rows]
        rank_xor = gf2_rank(xor_rows)

        # Rank of δΨ alone
        psi_rows = [row[256:] for row in jacobian_rows]
        rank_psi = gf2_rank(psi_rows)

        ranks.append((rank_joint, rank_xor, rank_psi, n_input))

        if sample < 5:
            print(f"\n  Sample {sample}: {n_input} input bits")
            print(f"    rank(J_joint) = {rank_joint} / {n_input}")
            print(f"    rank(J_xor)   = {rank_xor} / {n_input}")
            print(f"    rank(J_psi)   = {rank_psi} / {n_input}")
            print(f"    Redundancy    = {rank_xor + rank_psi - rank_joint}")
            print(f"      (>0 means δF_xor and δΨ share information!)")

    avg_joint = sum(r[0] for r in ranks)/len(ranks)
    avg_xor = sum(r[1] for r in ranks)/len(ranks)
    avg_psi = sum(r[2] for r in ranks)/len(ranks)
    avg_redundancy = sum(r[1]+r[2]-r[0] for r in ranks)/len(ranks)

    print(f"\n  Summary ({N_samples} samples, {ranks[0][3]} input bits each):")
    print(f"    E[rank(J_joint)] = {avg_joint:.1f}")
    print(f"    E[rank(J_xor)]   = {avg_xor:.1f}")
    print(f"    E[rank(J_psi)]   = {avg_psi:.1f}")
    print(f"    E[redundancy]    = {avg_redundancy:.1f}")
    print(f"      Redundancy = bits shared between XOR and Ψ Jacobians")

    return avg_redundancy


# ============================================================
# EXP 2: Per-bit correlation between δF_xor and δΨ
# ============================================================
def exp2_per_bit_correlation(N):
    print("\n" + "="*70)
    print("EXP 2: PER-BIT CORRELATION between δF_xor and δΨ")
    print("  For each output bit: corr(δF_xor[bit], δΨ[bit])")
    print("  High correlation → redundancy → fewer real constraints")
    print("="*70)

    # For each of 256 output bits, count agreement
    agree = [0]*256
    total = [0]*256

    for trial in range(N):
        W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        delta_W = list(W)
        delta_W[0] ^= (1 << (trial % 32))

        xor1 = sha_xor(W); xor2 = sha_xor(delta_W)
        real1 = sha_real(W); real2 = sha_real(delta_W)

        d_xor = state_to_bits([xor1[i]^xor2[i] for i in range(8)])
        d_real = state_to_bits([real1[i]^real2[i] for i in range(8)])
        d_psi = [d_real[i]^d_xor[i] for i in range(256)]

        for b in range(256):
            total[b] += 1
            if d_xor[b] == d_psi[b]:
                agree[b] += 1

    corr = [agree[b]/total[b] for b in range(256)]

    # Summary by register
    regs = ['a','b','c','d','e','f','g','h']
    print(f"\n  Per-register avg correlation P(δF_xor[bit] = δΨ[bit]):")
    print(f"  (0.5 = independent, >0.5 = correlated, <0.5 = anti-correlated)")
    for reg_idx in range(8):
        reg_corr = corr[reg_idx*32:(reg_idx+1)*32]
        avg = sum(reg_corr)/32
        mn = min(reg_corr)
        mx = max(reg_corr)
        marker = " ★" if abs(avg - 0.5) > 0.02 else ""
        print(f"    {regs[reg_idx]}: avg={avg:.4f} min={mn:.3f} max={mx:.3f}{marker}")

    # By bit position across all registers
    print(f"\n  By bit position (averaged across registers):")
    print(f"  {'bit':>4} | corr")
    for b in range(0, 32, 4):
        avg_b = sum(corr[reg*32+b] for reg in range(8))/8
        bar = "█" * int(abs(avg_b-0.5)*200)
        print(f"  {b:4d} | {avg_b:.4f} {bar}")

    overall_avg = sum(corr)/256
    print(f"\n  Overall average: {overall_avg:.4f}")
    print(f"  Deviation from 0.5: {overall_avg-0.5:+.4f}")

    return corr


# ============================================================
# EXP 3: Effective matching dimension via sampling
# ============================================================
def exp3_effective_dimension(N):
    print("\n" + "="*70)
    print("EXP 3: EFFECTIVE MATCHING DIMENSION")
    print("  Sample random δ, compute M(δ) = δF_xor ⊕ δΨ = δF_real")
    print("  Measure: how many independent M(δ) vectors in 256-dim space?")
    print("  Full rank = 256 → no shortcut. Rank < 256 → shortcut exists.")
    print("="*70)

    for n_words in [1, 2, 4, 8, 16]:
        # Fix a random W, sample random deltas affecting n_words
        W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

        rows = []
        for trial in range(min(N*2, 512)):
            delta_W = list(W)
            for w in range(n_words):
                delta_W[w] ^= int.from_bytes(os.urandom(4),'big')

            real1 = sha_real(W); real2 = sha_real(delta_W)
            d_real = state_to_bits([real1[i]^real2[i] for i in range(8)])
            rows.append(d_real)

            if len(rows) >= 300:
                break

        rank = gf2_rank(rows)
        n_samples = len(rows)
        print(f"  δW affects {n_words} words: rank = {rank}/{min(256, n_samples)} "
              f"({n_samples} samples)")


# ============================================================
# EXP 4: Ψ-kernel — directions where δΨ = 0
# ============================================================
def exp4_psi_kernel(N):
    print("\n" + "="*70)
    print("EXP 4: Ψ-KERNEL — directions where δΨ ≈ 0")
    print("  If δΨ = 0 for some δ, then collision reduces to δF_xor = 0")
    print("  δF_xor = 0 is a QUADRATIC system (much simpler!)")
    print("="*70)

    # For a fixed W, search for δ where HW(δΨ) is minimal
    W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

    best_hw_psi = 256
    best_delta_desc = ""

    # Strategy A: Single-bit deltas
    for word in range(16):
        for bit in range(32):
            delta_W = list(W)
            delta_W[word] ^= (1 << bit)

            real1 = sha_real(W); real2 = sha_real(delta_W)
            xor1 = sha_xor(W); xor2 = sha_xor(delta_W)

            d_real = [real1[i]^real2[i] for i in range(8)]
            d_xor = [xor1[i]^xor2[i] for i in range(8)]
            d_psi = [d_real[i]^d_xor[i] for i in range(8)]

            hw_psi = sum(hw(x) for x in d_psi)

            if hw_psi < best_hw_psi:
                best_hw_psi = hw_psi
                best_delta_desc = f"W[{word}] bit {bit}"

    print(f"\n  Single-bit δ: best HW(δΨ) = {best_hw_psi} ({best_delta_desc})")

    # Strategy B: SA search for low δΨ
    current_W2 = list(W)
    current_W2[0] ^= 1  # Start with small delta

    def compute_dpsi(W1, W2):
        r1 = sha_real(W1); r2 = sha_real(W2)
        x1 = sha_xor(W1); x2 = sha_xor(W2)
        return sum(hw((r1[i]^r2[i])^(x1[i]^x2[i])) for i in range(8))

    best_score = compute_dpsi(W, current_W2)
    best_W2 = list(current_W2)

    for it in range(min(N*20, 20000)):
        trial = list(current_W2)
        word = int.from_bytes(os.urandom(1),'big') % 16
        trial[word] ^= (1 << (int.from_bytes(os.urandom(1),'big') % 32))

        # Don't allow δ=0
        if trial == W:
            continue

        score = compute_dpsi(W, trial)

        T = max(0.01, 1.0 - it/(N*20))
        if score < best_score or \
           (score == best_score and int.from_bytes(os.urandom(1),'big') < 128) or \
           math.exp(-(score-best_score)/(T*5)) > int.from_bytes(os.urandom(4),'big')/(1<<32):
            current_W2 = trial
            if score < best_score:
                best_score = score
                best_W2 = list(trial)

    print(f"  SA search: best HW(δΨ) = {best_score}")

    # How many words differ?
    n_diff_words = sum(1 for i in range(16) if W[i] != best_W2[i])
    total_hw_delta = sum(hw(W[i]^best_W2[i]) for i in range(16))
    print(f"    δ affects {n_diff_words} words, total HW(δ) = {total_hw_delta}")

    # What is δF_xor for this near-kernel δ?
    x1 = sha_xor(W); x2 = sha_xor(best_W2)
    d_xor_hw = sum(hw(x1[i]^x2[i]) for i in range(8))
    r1 = sha_real(W); r2 = sha_real(best_W2)
    d_real_hw = sum(hw(r1[i]^r2[i]) for i in range(8))

    print(f"    HW(δF_real) = {d_real_hw}")
    print(f"    HW(δF_xor)  = {d_xor_hw}")
    print(f"    HW(δΨ)      = {best_score}")
    print(f"    Note: δF_real = δF_xor ⊕ δΨ")

    if best_score < 100:
        print(f"\n  ★ Near-kernel found! {256-best_score} of 256 Ψ-bits are zero.")
        print(f"    In those {256-best_score} positions, collision = XOR-collision!")

    # Distribution of δΨ for random deltas
    rand_psi = []
    for _ in range(min(N, 1000)):
        dW = list(W)
        dW[int.from_bytes(os.urandom(1),'big')%16] ^= int.from_bytes(os.urandom(4),'big')
        if dW == W: continue
        rand_psi.append(compute_dpsi(W, dW))

    if rand_psi:
        print(f"\n  Distribution of HW(δΨ) for random δ:")
        print(f"    avg = {sum(rand_psi)/len(rand_psi):.1f}")
        print(f"    min = {min(rand_psi)}, max = {max(rand_psi)}")
        print(f"    SA improvement over random min: {min(rand_psi) - best_score} bits")


# ============================================================
# EXP 5: Matching condition structure
# ============================================================
def exp5_matching_structure(N):
    print("\n" + "="*70)
    print("EXP 5: MATCHING STRUCTURE — where does δπ_⊕ = δπ_Ψ fail?")
    print("  For random δ, count per-bit: P(δF_xor[b] = δΨ[b])")
    print("  Bits where this is high → 'easy' constraints")
    print("  Bits where this is ~0.5 → 'hard' constraints")
    print("="*70)

    # Count per-bit match probability
    match_count = [0]*256
    total = 0

    for trial in range(min(N, 1000)):
        W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        dW = list(W)
        dW[trial%16] ^= (1 << (trial%32))

        r1=sha_real(W); r2=sha_real(dW)
        x1=sha_xor(W); x2=sha_xor(dW)

        d_xor = state_to_bits([x1[i]^x2[i] for i in range(8)])
        d_psi = state_to_bits([(r1[i]^r2[i])^(x1[i]^x2[i]) for i in range(8)])

        total += 1
        for b in range(256):
            if d_xor[b] == d_psi[b]:
                match_count[b] += 1

    probs = [match_count[b]/total for b in range(256)]

    # Count "easy" and "hard" bits
    easy = sum(1 for p in probs if p > 0.6)
    hard = sum(1 for p in probs if 0.45 < p < 0.55)
    very_easy = sum(1 for p in probs if p > 0.7)

    print(f"\n  256 output bits classified:")
    print(f"    Very easy (P > 0.7):  {very_easy}")
    print(f"    Easy (P > 0.6):       {easy}")
    print(f"    Hard (0.45 < P < 0.55): {hard}")
    print(f"    Total bits:           256")

    print(f"\n  Effective constraint dimension ≈ {hard}")
    print(f"  'Free' dimensions ≈ {256 - hard}")

    if very_easy > 10:
        print(f"\n  ★ {very_easy} bits are 'very easy' — matching nearly automatic!")
        print(f"    These bits have δF_xor ≈ δΨ with P > 0.7")
        print(f"    Only need to solve for the remaining {256-very_easy} bits")

    # Per-register breakdown
    regs = ['a','b','c','d','e','f','g','h']
    print(f"\n  Per-register breakdown:")
    for reg in range(8):
        reg_probs = probs[reg*32:(reg+1)*32]
        easy_r = sum(1 for p in reg_probs if p > 0.6)
        hard_r = sum(1 for p in reg_probs if 0.45 < p < 0.55)
        avg_r = sum(reg_probs)/32
        print(f"    {regs[reg]}: avg_match={avg_r:.3f}, easy={easy_r}/32, hard={hard_r}/32")


# ============================================================
if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 200

    print("="*70)
    print("SCF: MATCHING CONDITION DIMENSION ANALYSIS")
    print("  Is δπ_⊕ = δπ_Ψ really 256 bits, or less?")
    print("="*70)

    redundancy = exp1_matching_jacobian(N_samples=min(N//4, 30))
    corr = exp2_per_bit_correlation(min(N, 500))
    exp3_effective_dimension(N)
    exp4_psi_kernel(N)
    exp5_matching_structure(N)

    print("\n" + "="*70)
    print("ФИНАЛЬНЫЙ СИНТЕЗ")
    print("="*70)
    print(f"""
  1. Jacobian redundancy = {redundancy:.1f} bits
     (bits shared between δF_xor and δΨ Jacobians)

  2. Per-bit correlation δF_xor vs δΨ:
     overall = {sum(corr)/256:.4f}
     (0.5 = independent → 256 real constraints)
     (>0.5 = correlated → fewer real constraints)

  3. ВЫВОД:
     Если redundancy > 0 → matching < 256 constraints → shortcut
     Если redundancy ≈ 0 → matching = 256 → no shortcut → birthday bound

  Формула:
     Effective collision cost = 2^(matching_dimension / 2)
     If matching_dim = 256 → cost = 2^128 (birthday)
     If matching_dim = 240 → cost = 2^120 (8 bits saved!)
     If matching_dim = 200 → cost = 2^100 (28 bits saved!)
    """)
