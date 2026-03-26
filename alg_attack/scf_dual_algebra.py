#!/usr/bin/env python3
"""
SCF: Dual Algebra Framework (DAF)
Третья алгебра для SHA-256 — видит одновременно GF(2) и Z_{2^32}.

ИДЕЯ: Наблюдатель O_⊕ видит L+Q, но не C (carry = шум).
      Наблюдатель O_+ видит + линейно, но не Σ (rotations = шум).
      Наблюдатель O_Ψ работает в пространстве ПАРЫ (x_⊕, x_+)
      и использует Ψ как СВЯЗУЮЩИЙ оператор.

Формально:
  Состояние в DAF: s_dual = (s_gf2, s_int, Ψ_coupling)
  где Ψ_coupling[r] = Ψ(T1_operands) на каждом раунде

Вопрос: в dual representation, какова ЭФФЕКТИВНАЯ степень SHA?
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
def psi(a,b): return ((a+b)&MASK32)^(a^b)
def psi_multi(*args):
    ar,xr = 0,0
    for x in args: ar=(ar+x)&MASK32; xr^=x
    return ar^xr

def expand_schedule(W16):
    W=list(W16)
    for i in range(16,64):
        W.append(add32(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    return W

def expand_schedule_xor(W16):
    """Schedule with XOR instead of +"""
    W=list(W16)
    for i in range(16,64):
        W.append(sig1(W[i-2])^W[i-7]^sig0(W[i-15])^W[i-16])
    return W


# ============================================================
# DUAL STATE: track (state_real, state_xor, Ψ_coupling) simultaneously
# ============================================================
def sha_dual_trace(W16):
    """
    Run SHA-256 in three parallel universes:
    1. Real (arithmetic): standard SHA
    2. XOR (GF(2)): all + replaced by ⊕
    3. Ψ-coupling: the difference between 1 and 2

    Returns per-round: (state_real, state_xor, psi_state, psi_T1, psi_T2)
    """
    W_real = expand_schedule(W16)
    W_xor = expand_schedule_xor(W16)

    s_real = list(IV)
    s_xor = list(IV)

    trace = []

    for r in range(64):
        # --- Real universe ---
        a,b,c,d,e,f,g,h = s_real
        T1_real = add32(h, Sig1(e), Ch(e,f,g), K[r], W_real[r])
        T2_real = add32(Sig0(a), Maj(a,b,c))
        psi_T1 = psi_multi(h, Sig1(e), Ch(e,f,g), K[r], W_real[r])
        psi_T2 = psi_multi(Sig0(a), Maj(a,b,c))
        psi_e = psi(d, T1_real)
        psi_a = psi(T1_real, T2_real)
        s_real_new = [add32(T1_real,T2_real),a,b,c,add32(d,T1_real),e,f,g]

        # --- XOR universe ---
        ax,bx,cx,dx,ex,fx,gx,hx = s_xor
        T1_xor = hx ^ Sig1(ex) ^ Ch(ex,fx,gx) ^ K[r] ^ W_xor[r]
        T2_xor = Sig0(ax) ^ Maj(ax,bx,cx)
        s_xor_new = [T1_xor^T2_xor,ax,bx,cx,dx^T1_xor,ex,fx,gx]

        # --- Ψ coupling: real ⊕ xor ---
        psi_state = [s_real_new[i] ^ s_xor_new[i] for i in range(8)]

        trace.append({
            'r': r,
            's_real': list(s_real_new),
            's_xor': list(s_xor_new),
            'psi_state': psi_state,
            'psi_T1': psi_T1,
            'psi_T2': psi_T2,
            'psi_e': psi_e,
            'psi_a': psi_a,
            'hw_psi_state': sum(hw(x) for x in psi_state),
            'hw_psi_T1': hw(psi_T1),
            'hw_psi_T2': hw(psi_T2),
        })

        s_real = s_real_new
        s_xor = s_xor_new

    return trace


# ============================================================
# EXP 1: Ψ-coupling profile — how does the coupling evolve?
# ============================================================
def exp1_psi_profile(N):
    print("="*70)
    print("EXP 1: Ψ-COUPLING PROFILE")
    print("  Track HW(Ψ_state) = HW(state_real ⊕ state_xor) per round")
    print("  This measures how far the two universes diverge")
    print("="*70)

    avg_hw = [0.0]*64
    avg_T1 = [0.0]*64
    avg_T2 = [0.0]*64

    for _ in range(N):
        W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        trace = sha_dual_trace(W)
        for r in range(64):
            avg_hw[r] += trace[r]['hw_psi_state']
            avg_T1[r] += trace[r]['hw_psi_T1']
            avg_T2[r] += trace[r]['hw_psi_T2']

    for r in range(64):
        avg_hw[r] /= N
        avg_T1[r] /= N
        avg_T2[r] /= N

    print(f"\n  {'r':>3} | {'HW(Ψ_state)':>12} {'HW(Ψ_T1)':>10} {'HW(Ψ_T2)':>10} | {'phase':>10}")
    print("  " + "-"*60)
    for r in range(64):
        if r < 10 or r % 5 == 0 or r == 63:
            phase = "GROWTH" if avg_hw[r] < 120 else "SATURATED"
            if avg_hw[r] < 20: phase = "SMALL"
            print(f"  {r:3d} | {avg_hw[r]:12.1f} {avg_T1[r]:10.1f} {avg_T2[r]:10.1f} | {phase}")

    # When does Ψ saturate?
    for r in range(64):
        if avg_hw[r] > 120:
            print(f"\n  Ψ-coupling saturates at r={r} (HW > 120)")
            break

    return avg_hw, avg_T1, avg_T2


# ============================================================
# EXP 2: Ψ-coupling differential — does Ψ carry info about δW?
# ============================================================
def exp2_psi_differential(N):
    print("\n" + "="*70)
    print("EXP 2: Ψ-COUPLING DIFFERENTIAL")
    print("  For (W, W⊕δ): does δΨ_state carry info about δW?")
    print("  If yes → Ψ-space contains exploitable structure")
    print("="*70)

    corr_by_round = [0.0]*64
    hw_delta_psi = [[] for _ in range(64)]

    for trial in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = list(W1)
        W2[0] ^= (1 << (trial % 32))

        trace1 = sha_dual_trace(W1)
        trace2 = sha_dual_trace(W2)

        for r in range(64):
            # δΨ = Ψ(W1) ⊕ Ψ(W2)
            delta_psi = [trace1[r]['psi_state'][i] ^ trace2[r]['psi_state'][i] for i in range(8)]
            hw_dp = sum(hw(x) for x in delta_psi)
            hw_delta_psi[r].append(hw_dp)

            # Also compare with δstate_real
            delta_real = [trace1[r]['s_real'][i] ^ trace2[r]['s_real'][i] for i in range(8)]
            hw_dr = sum(hw(x) for x in delta_real)

            # Correlation: does δΨ predict δstate?
            if hw_dr > 0:
                # Normalized overlap
                overlap = sum(hw(delta_psi[i] & delta_real[i]) for i in range(8))
                corr_by_round[r] += overlap / hw_dr

    for r in range(64):
        corr_by_round[r] /= N

    print(f"\n  {'r':>3} | {'E[HW(δΨ)]':>11} {'σ':>6} | {'corr(δΨ,δreal)':>15} | {'status':>10}")
    print("  " + "-"*55)
    for r in range(64):
        if r < 10 or r % 5 == 0 or r == 63:
            m = sum(hw_delta_psi[r])/len(hw_delta_psi[r])
            s = (sum((x-m)**2 for x in hw_delta_psi[r])/len(hw_delta_psi[r]))**0.5
            c = corr_by_round[r]
            status = "STRUCTURED" if c > 0.1 else "RANDOM"
            if c > 0.3: status = "STRONG ★"
            print(f"  {r:3d} | {m:11.1f} {s:6.1f} | {c:15.3f} | {status}")


# ============================================================
# EXP 3: Ψ decomposition — is Ψ itself decomposable?
# ============================================================
def exp3_psi_decompose(N):
    print("\n" + "="*70)
    print("EXP 3: Ψ DECOMPOSITION")
    print("  Ψ_total = Ψ_T1 ⊕ Ψ_T2 ⊕ Ψ_e ⊕ Ψ_a (4 sources)")
    print("  Which source dominates? Are they independent?")
    print("="*70)

    # Track per-source contributions
    src_names = ['Ψ_T1', 'Ψ_T2', 'Ψ_e', 'Ψ_a']
    src_hw = {name: [0.0]*64 for name in src_names}

    # Track pairwise correlations
    pair_corr = {}
    for i in range(4):
        for j in range(i+1,4):
            pair_corr[(src_names[i],src_names[j])] = [0.0]*64

    for _ in range(N):
        W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        trace = sha_dual_trace(W)
        for r in range(64):
            vals = [trace[r]['psi_T1'], trace[r]['psi_T2'],
                    trace[r]['psi_e'], trace[r]['psi_a']]
            hws = [hw(v) for v in vals]
            for i,name in enumerate(src_names):
                src_hw[name][r] += hws[i]
            # Pairwise: overlap
            for i in range(4):
                for j in range(i+1,4):
                    overlap = hw(vals[i] & vals[j])
                    expected = hws[i]*hws[j]/32.0
                    pair_corr[(src_names[i],src_names[j])][r] += (overlap - expected)

    for name in src_names:
        for r in range(64):
            src_hw[name][r] /= N
    for key in pair_corr:
        for r in range(64):
            pair_corr[key][r] /= N

    print(f"\n  {'r':>3} | {'Ψ_T1':>6} {'Ψ_T2':>6} {'Ψ_e':>6} {'Ψ_a':>6} | {'T1%':>5} {'T2%':>5}")
    print("  " + "-"*50)
    for r in range(64):
        if r < 10 or r % 5 == 0 or r == 63:
            vals = [src_hw[name][r] for name in src_names]
            total = sum(vals) if sum(vals) > 0 else 1
            print(f"  {r:3d} | {vals[0]:6.1f} {vals[1]:6.1f} {vals[2]:6.1f} {vals[3]:6.1f} | "
                  f"{100*vals[0]/total:5.1f} {100*vals[1]/total:5.1f}")

    # Independence test
    print(f"\n  Pairwise excess overlap (0 = independent, >0 = correlated):")
    print(f"  {'r':>3} | " + " ".join(f"{k[0][-2:]}/{k[1][-2:]}" for k in sorted(pair_corr.keys())))
    for r in [0, 5, 10, 20, 63]:
        vals = [pair_corr[k][r] for k in sorted(pair_corr.keys())]
        print(f"  {r:3d} | " + " ".join(f"{v:+5.2f}" for v in vals))


# ============================================================
# EXP 4: Dual state rank — information content of Ψ-space
# ============================================================
def exp4_psi_information(N):
    print("\n" + "="*70)
    print("EXP 4: Ψ INFORMATION CONTENT")
    print("  How many bits of info does Ψ_state carry about W?")
    print("  Method: flip each W bit, measure HW(δΨ_state)")
    print("  If HW(δΨ)≈128 → Ψ is random (no info)")
    print("  If HW(δΨ)<128 → Ψ retains structure")
    print("="*70)

    # For each W word position and bit, measure δΨ at r=64
    results = {}
    for word in range(16):
        for bit in range(0, 32, 8):  # Sample every 8th bit
            delta_psi_list = []
            for _ in range(min(N, 200)):
                W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
                W2 = list(W1)
                W2[word] ^= (1 << bit)

                t1 = sha_dual_trace(W1)
                t2 = sha_dual_trace(W2)

                # δΨ at final round
                dpsi = [t1[63]['psi_state'][i] ^ t2[63]['psi_state'][i] for i in range(8)]
                delta_psi_list.append(sum(hw(x) for x in dpsi))

            avg = sum(delta_psi_list)/len(delta_psi_list)
            results[(word, bit)] = avg

    print(f"\n  E[HW(δΨ_state[64])] by (word, bit) — expected 128 if random:")
    print(f"  {'word':>5} | {'bit0':>6} {'bit8':>6} {'bit16':>6} {'bit24':>6} | {'avg':>6}")
    print("  " + "-"*50)
    for word in range(16):
        vals = [results.get((word, b), 0) for b in [0, 8, 16, 24]]
        avg = sum(vals)/4
        marker = " ★" if abs(avg - 128) > 5 else ""
        print(f"  W[{word:2d}] | {vals[0]:6.1f} {vals[1]:6.1f} {vals[2]:6.1f} {vals[3]:6.1f} | {avg:6.1f}{marker}")


# ============================================================
# EXP 5: The key question — Ψ-guided collision search
# ============================================================
def exp5_psi_guided_search(N):
    print("\n" + "="*70)
    print("EXP 5: Ψ-GUIDED COLLISION SEARCH")
    print("  Strategy: minimize HW(Ψ_state) instead of HW(δstate)")
    print("  If Ψ=0, then real=xor, and XOR-SHA is simpler!")
    print("="*70)

    # Can we find W where Ψ_state is small?
    def psi_state_hw(W16):
        trace = sha_dual_trace(W16)
        return trace[63]['hw_psi_state']

    # SA to minimize Ψ_state at round 64
    best_W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    best_score = psi_state_hw(best_W)
    current_W = list(best_W)
    current_score = best_score
    init_score = best_score

    for it in range(min(N*10, 10000)):
        trial_W = list(current_W)
        word_idx = int.from_bytes(os.urandom(1),'big') % 16
        trial_W[word_idx] ^= (1 << (int.from_bytes(os.urandom(1),'big') % 32))

        trial_score = psi_state_hw(trial_W)

        T = max(0.01, 1.0 - it/(N*10))
        if trial_score < current_score or \
           math.exp(-(trial_score-current_score)/(T*5)) > \
           int.from_bytes(os.urandom(4),'big')/(1<<32):
            current_W = trial_W
            current_score = trial_score

        if current_score < best_score:
            best_score = current_score
            best_W = list(current_W)

    # Baseline
    rand_scores = []
    for _ in range(min(N, 500)):
        W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        rand_scores.append(psi_state_hw(W))
    rand_avg = sum(rand_scores)/len(rand_scores)
    rand_min = min(rand_scores)

    print(f"\n  HW(Ψ_state[64]) — how different are real and XOR universes?")
    print(f"    Random:    avg={rand_avg:.1f}, min={rand_min}")
    print(f"    Optimized: {best_score}")
    print(f"    Reduction: {(rand_avg-best_score)/rand_avg*100:.1f}%")

    if best_score < rand_min:
        print(f"    ★ BELOW RANDOM MINIMUM — Ψ is partially controllable!")
    else:
        print(f"    Within random range — Ψ not controllable at r=64")

    # Show the Ψ profile for optimized W
    trace_opt = sha_dual_trace(best_W)
    print(f"\n  Ψ profile for optimized W:")
    print(f"  {'r':>3} | {'HW(Ψ)':>6}")
    for r in range(0, 64, 4):
        print(f"  {r:3d} | {trace_opt[r]['hw_psi_state']:6d}")


# ============================================================
# EXP 6: Dual differential — δ in Ψ-space
# ============================================================
def exp6_dual_differential(N):
    print("\n" + "="*70)
    print("EXP 6: DUAL DIFFERENTIAL — collision in Ψ-space")
    print("  Find (W1, W2) where Ψ(W1) = Ψ(W2) at round 64")
    print("  This means: carry patterns are IDENTICAL")
    print("  Then δstate = δstate_xor (purely GF(2) problem!)")
    print("="*70)

    best_hw = 256
    best_pair = None

    for trial in range(min(N, 2000)):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = list(W1)
        # Small perturbation
        W2[trial % 16] ^= (1 << (trial % 32))

        t1 = sha_dual_trace(W1)
        t2 = sha_dual_trace(W2)

        # δΨ at round 64
        delta_psi = [t1[63]['psi_state'][i] ^ t2[63]['psi_state'][i] for i in range(8)]
        hw_dpsi = sum(hw(x) for x in delta_psi)

        if hw_dpsi < best_hw:
            best_hw = hw_dpsi
            best_pair = (W1, W2, trial)

    print(f"\n  Best HW(δΨ_state[64]): {best_hw} (out of 256 bits)")
    print(f"  Expected random: ~128")

    if best_hw < 100:
        print(f"  ★ Significant Ψ-collision: {best_hw} bits differ in carry space")
        print(f"    This means {256-best_hw} bits have IDENTICAL carry patterns")
    else:
        print(f"  No significant Ψ-collision found")


# ============================================================
# MAIN
# ============================================================
if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 300

    print("="*70)
    print("SCF: DUAL ALGEBRA FRAMEWORK (DAF)")
    print("Третья алгебра: пространство Ψ-связи между GF(2) и Z_{2^32}")
    print("="*70)

    hw_profile, _, _ = exp1_psi_profile(min(N, 500))
    exp2_psi_differential(min(N, 300))
    exp3_psi_decompose(min(N, 500))
    exp4_psi_information(N)
    exp5_psi_guided_search(N)
    exp6_dual_differential(N)

    print("\n" + "="*70)
    print("СИНТЕЗ: DUAL ALGEBRA FRAMEWORK")
    print("="*70)
    print("""
  ФОРМАЛИЗАЦИЯ:

  Пусть F: {0,1}^512 → {0,1}^256 — SHA-256 compress.

  Определим три проекции:
    π_⊕(x) = F_xor(x)     — XOR-вселенная (линейная при r≤3)
    π_+(x) = F_real(x)     — Real-вселенная (стандартная SHA)
    π_Ψ(x) = π_+(x) ⊕ π_⊕(x)  — Ψ-связь (carry-вселенная)

  Тогда: F(x) = π_⊕(x) ⊕ π_Ψ(x)

  Collision в F:   F(x) = F(y)
  ⟺ π_⊕(x) ⊕ π_Ψ(x) = π_⊕(y) ⊕ π_Ψ(y)
  ⟺ δπ_⊕ = δπ_Ψ

  Т.е. collision ⟺ XOR-differential = Ψ-differential!

  Это РАЗДЕЛЯЕТ задачу на две:
    1. Найти δπ_⊕ (GF(2) — квадратичная система)
    2. Найти δπ_Ψ (carry — марковская система)
    3. Потребовать δπ_⊕ = δπ_Ψ (matching condition)

  ВОПРОС: Легче ли решать 1+2+3 чем исходную задачу?
    """)
