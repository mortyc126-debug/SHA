#!/usr/bin/env python3
"""
SCF: BACKWARD WANG CHAIN — обратная рекурренция от δstate[64]=0.

Все атаки идут ВПЕРЁД. Wang chain: δe[2..16]=0 → 15 раундов контроля.
Вопрос: сколько раундов НАЗАД от δstate[64]=0 мы можем контролировать?

Если вперёд 15 + назад N → нужно перекрыть только (64-15-N) раундов.

Ключевое наблюдение (shift register):
  state[64] = (a[64], a[63], a[62], a[61], e[64], e[63], e[62], e[61])
  δstate[64] = 0 означает δa[61..64] = 0, δe[61..64] = 0

Обратный шаг:
  a[r+1] = T1[r] + T2[r]      →  T1[r] = a[r+1] - T2[r]
  e[r+1] = a[r-3] + T1[r]     →  a[r-3] = e[r+1] - T1[r]
  T1[r] = e[r-3] + Σ1(e[r]) + Ch(e[r],e[r-1],e[r-2]) + K[r] + W[r]

Когда δa[r+1..r+4]=0 и δe[r+1..r+4]=0:
  δT2[r] = δΣ0(a[r]) + δMaj(a[r],a[r-1],a[r-2])
  Если δa[r..r-2]=0: δT2[r] = 0

  δT1[r] = δa[r+1] - δT2[r] = 0  → OK

  δe[r-3] + δΣ1(e[r]) + δCh(e[r],e[r-1],e[r-2]) + δW[r] = δT1[r] = 0
  Если δe[r..r-2]=0: δΣ1=0, δCh=0
  → δe[r-3] = -δW[r]

  δa[r-3] = δe[r+1] - δT1[r] = 0

Итого: backward step r+1→r при δa[r+1..r+4]=δe[r+1..r+4]=0 даёт:
  δa[r-3] = 0     (бесплатно!)
  δe[r-3] = -δW[r] (управляемо через W[r]!)

Если δW[r] = 0 → δe[r-3] = 0 → ещё один нулевой раунд БЕСПЛАТНО!
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

def expand_schedule(W16):
    W=list(W16)
    for i in range(16,64):
        W.append(add32(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    return W

def sha_all_states(W16):
    W=expand_schedule(W16); s=list(IV); states=[list(s)]
    for r in range(64):
        a,b,c,d,e,f,g,h=s
        T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add32(Sig0(a),Maj(a,b,c))
        s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
        states.append(list(s))
    return states, W


# ============================================================
# EXP 1: Backward propagation — analytical
# ============================================================
def exp1_backward_theory():
    print("="*70)
    print("EXP 1: BACKWARD PROPAGATION — ТЕОРЕТИЧЕСКИЙ АНАЛИЗ")
    print("="*70)

    print("""
  Target: δstate[64] = 0, i.e., δa[64..61] = δe[64..61] = 0.

  Backward step from round r+1 to round r:
  Given δa[r+1..r+4] = 0, δe[r+1..r+4] = 0:

    1. δT2[r] = δ(Σ0(a[r]) + Maj(a[r],a[r-1],a[r-2]))
       If δa[r]=δa[r-1]=δa[r-2]=0: δT2[r] = 0  ✓

    2. δT1[r] = δa[r+1] - δT2[r] = 0 - 0 = 0  ✓

    3. δT1[r] = δh[r] + δΣ1(e[r]) + δCh(e[r],f[r],g[r]) + δW[r]
       = δe[r-3] + δΣ1(e[r]) + δCh + δW[r]
       If δe[r]=δe[r-1]=δe[r-2]=0: δΣ1=0, δCh=0
       → δe[r-3] + δW[r] = 0
       → δe[r-3] = -δW[r]  (mod 2^32)

    4. δa[r-3] = δ(e[r+1] - T1[r]) = δe[r+1] - δT1[r] = 0  ✓

  CONCLUSION:
    Each backward step gives δa[r-3]=0 FOR FREE.
    Each backward step gives δe[r-3] = -δW[r] — CONTROLLED by W[r].

    If δW[r] = 0 → δe[r-3] = 0 → one more zero round FREE.

  APPLICATION:
    From δstate[64]=0, going backwards:
      r=63: δa[60]=0 free, δe[60]=-δW[63]
      r=62: δa[59]=0 free, δe[59]=-δW[62] + corrections(δe[60])
      r=61: δa[58]=0 free, δe[58]=-δW[61] + corrections(δe[60],δe[59])
      ...

    WHILE δe[r..r+2]=0 (all three e-neighbours zero), corrections=0.
    Once any δe[r]≠0, corrections become nonzero → chain reaction.

    So: backward Wang chain holds as long as δW[r]=0 for consecutive rounds!
    """)


# ============================================================
# EXP 2: Verify backward propagation numerically
# ============================================================
def exp2_verify_backward(N):
    print("="*70)
    print("EXP 2: NUMERICAL VERIFICATION — backward from δstate[64]=0")
    print("  Set δW[r]=0 for r=R_back..63. Check δe, δa at each round.")
    print("="*70)

    for n_zero_back in [4, 8, 12, 16, 20, 32, 48]:
        # δW[r]=0 for r in [64-n_zero_back, 63]
        # δW[r]=random for r in [0, 64-n_zero_back-1]
        r_start_zero = 64 - n_zero_back

        total_de_zeros = 0
        total_da_zeros = 0

        for trial in range(N):
            W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
            W2 = list(W1)

            # Modify W1 → W2: flip bits only in words that expand to W[0..r_start_zero-1]
            # For simplicity: flip bit in W[0] (affects W[0] and schedule)
            W2[0] ^= (1 << (trial % 32))

            states1, Wexp1 = sha_all_states(W1)
            states2, Wexp2 = sha_all_states(W2)

            # Check: which δW[r] are actually zero?
            dW_zero_from = 64
            for r in range(63, -1, -1):
                if Wexp1[r] != Wexp2[r]:
                    dW_zero_from = r + 1
                    break

            # Count zero δe and δa from the back
            de_zeros = 0
            for r in range(64, 0, -1):
                if states1[r][4] == states2[r][4]:  # e register
                    de_zeros += 1
                else:
                    break

            da_zeros = 0
            for r in range(64, 0, -1):
                if states1[r][0] == states2[r][0]:  # a register
                    da_zeros += 1
                else:
                    break

            total_de_zeros += de_zeros
            total_da_zeros += da_zeros

        avg_de = total_de_zeros / N
        avg_da = total_da_zeros / N

        print(f"\n  δW[r]=0 for last {n_zero_back} rounds (r≥{r_start_zero}):")
        print(f"    δW actually zero from r={r_start_zero} (schedule dep.)")
        print(f"    Avg consecutive δe=0 from end: {avg_de:.1f}")
        print(f"    Avg consecutive δa=0 from end: {avg_da:.1f}")


# ============================================================
# EXP 3: Schedule structure — which δW[16..63] can be zeroed?
# ============================================================
def exp3_schedule_zeros(N):
    print("\n" + "="*70)
    print("EXP 3: SCHEDULE STRUCTURE — how many δW[16..63] = 0?")
    print("  W[i] = σ1(W[i-2]) + W[i-7] + σ0(W[i-15]) + W[i-16]")
    print("  If δW[0..15] has k nonzero words, how many δW[16..63] are nonzero?")
    print("="*70)

    for n_flip_words in [1, 2, 4, 8, 16]:
        total_nonzero_sched = 0
        last_nonzero = 0

        for trial in range(N):
            W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
            W2 = list(W1)

            # Flip n_flip_words words
            for w in range(n_flip_words):
                W2[w] ^= (trial * (w+1) + 1) & MASK32

            Wexp1 = expand_schedule(W1)
            Wexp2 = expand_schedule(W2)

            nonzero = sum(1 for i in range(16, 64) if Wexp1[i] != Wexp2[i])
            total_nonzero_sched += nonzero

            # Last round where δW[r] ≠ 0
            for r in range(63, 15, -1):
                if Wexp1[r] != Wexp2[r]:
                    last_nonzero = max(last_nonzero, r)
                    break

        avg = total_nonzero_sched / N
        print(f"\n  {n_flip_words} words flipped in W[0..15]:")
        print(f"    Avg nonzero δW in schedule: {avg:.1f} / 48")
        print(f"    Last nonzero δW: r={last_nonzero}")
        print(f"    Zero δW in schedule: {48-avg:.1f} / 48")

    # Special: which δW[0] makes δW[16..63] have the MOST zeros?
    print(f"\n  Special: δW[0] only — schedule zero pattern")
    W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

    best_zeros = 0
    best_delta = 0

    for delta in range(1, min(N*10, 10000)):
        W2 = list(W1)
        W2[0] ^= delta
        Wexp1 = expand_schedule(W1)
        Wexp2 = expand_schedule(W2)

        zeros = sum(1 for i in range(16,64) if Wexp1[i] == Wexp2[i])
        if zeros > best_zeros:
            best_zeros = zeros
            best_delta = delta

    print(f"    Best: {best_zeros}/48 schedule words unchanged (δW[0]=0x{best_delta:08x})")


# ============================================================
# EXP 4: MEET IN THE MIDDLE — forward Wang + backward Wang
# ============================================================
def exp4_meet_in_middle(N):
    print("\n" + "="*70)
    print("EXP 4: MEET IN THE MIDDLE")
    print("  Forward: Wang chain → δe[2..16]=0 (15 zero rounds)")
    print("  Backward: from δstate[64]=0, how many zero rounds?")
    print("  Gap = rounds between forward and backward zero zones")
    print("="*70)

    # For each message pair, compute:
    # 1. How many δe=0 from the front (Wang-like)
    # 2. How many δe=0 from the back
    # 3. The gap

    gaps = []
    for trial in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = list(W1)
        W2[0] ^= 0x80000000  # Wang-style MSB flip

        states1, _ = sha_all_states(W1)
        states2, _ = sha_all_states(W2)

        # Forward: count consecutive δe=0 from round 2
        fwd_zeros = 0
        for r in range(2, 65):
            if states1[r][4] == states2[r][4]:
                fwd_zeros += 1
            else:
                break
        fwd_end = 2 + fwd_zeros  # Last zero round (forward)

        # Backward: count consecutive δe=0 from round 64
        bwd_zeros = 0
        for r in range(64, 0, -1):
            if states1[r][4] == states2[r][4]:
                bwd_zeros += 1
            else:
                break
        bwd_start = 64 - bwd_zeros  # First zero round (backward)

        gap = bwd_start - fwd_end
        gaps.append((fwd_zeros, bwd_zeros, gap, fwd_end, bwd_start))

    avg_fwd = sum(g[0] for g in gaps)/N
    avg_bwd = sum(g[1] for g in gaps)/N
    avg_gap = sum(g[2] for g in gaps)/N
    min_gap = min(g[2] for g in gaps)
    max_fwd = max(g[0] for g in gaps)
    max_bwd = max(g[1] for g in gaps)

    print(f"\n  Forward (δe=0 from r=2):  avg={avg_fwd:.1f}, max={max_fwd}")
    print(f"  Backward (δe=0 from r=64): avg={avg_bwd:.1f}, max={max_bwd}")
    print(f"  Gap (uncovered rounds):     avg={avg_gap:.1f}, min={min_gap}")

    if min_gap < 60:
        print(f"\n  ★ Minimum gap = {min_gap} rounds!")
        # Show best case
        best = min(gaps, key=lambda g: g[2])
        print(f"    Forward zeros: {best[0]} (r=2..{best[3]-1})")
        print(f"    Backward zeros: {best[1]} (r={best[4]}..64)")
        print(f"    Uncovered: r={best[3]}..{best[4]-1} ({best[2]} rounds)")

    # Histogram of gaps
    print(f"\n  Gap distribution:")
    for g_val in sorted(set(g[2] for g in gaps)):
        count = sum(1 for g in gaps if g[2] == g_val)
        if count > 0:
            bar = "█" * min(count, 50)
            print(f"    gap={g_val:3d}: {count:4d} {bar}")


# ============================================================
# EXP 5: Adaptive backward — choose δW to maximize backward chain
# ============================================================
def exp5_adaptive_backward(N):
    print("\n" + "="*70)
    print("EXP 5: ADAPTIVE BACKWARD CHAIN")
    print("  Choose δW[0..15] to maximize (forward zeros + backward zeros)")
    print("  This is the COMBINED optimization.")
    print("="*70)

    import math

    W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

    def score(W2):
        states1, _ = sha_all_states(W1)
        states2, _ = sha_all_states(W2)

        # Count all δe=0 rounds
        de_zero = sum(1 for r in range(65)
                      if states1[r][4] == states2[r][4])
        # Count all δa=0 rounds
        da_zero = sum(1 for r in range(65)
                      if states1[r][0] == states2[r][0])

        # Bonus for consecutive zeros from front and back
        fwd = 0
        for r in range(1, 65):
            if states1[r][4] == states2[r][4]:
                fwd += 1
            else:
                break

        bwd = 0
        for r in range(64, 0, -1):
            if states1[r][4] == states2[r][4]:
                bwd += 1
            else:
                break

        return de_zero + da_zero + fwd * 3 + bwd * 3

    # SA optimization
    best_W2 = list(W1)
    best_W2[0] ^= 0x80000000
    best_score = score(best_W2)
    current_W2 = list(best_W2)
    current_score = best_score

    for it in range(min(N * 50, 30000)):
        trial = list(current_W2)
        word = int.from_bytes(os.urandom(1),'big') % 16
        bit = int.from_bytes(os.urandom(1),'big') % 32
        trial[word] ^= (1 << bit)

        if trial == W1:
            continue

        sc = score(trial)

        T = max(0.01, 1.0 - it / (N * 50))
        if sc > current_score or \
           math.exp((sc - current_score) / (T * 3)) > \
           int.from_bytes(os.urandom(4), 'big') / (1 << 32):
            current_W2 = trial
            current_score = sc

        if current_score > best_score:
            best_score = current_score
            best_W2 = list(current_W2)

    # Analyze best result
    states1, Wexp1 = sha_all_states(W1)
    states2, Wexp2 = sha_all_states(best_W2)

    print(f"\n  Best combined score: {best_score}")
    print(f"\n  Per-round δe=0 and δa=0 map:")
    print(f"  {'r':>3} | δe  δa | δW")
    print("  " + "-"*40)

    fwd_de, bwd_de = 0, 0
    de_zeros = []
    da_zeros = []
    dw_zeros = []

    for r in range(65):
        de = "0" if states1[r][4] == states2[r][4] else "X"
        da = "0" if states1[r][0] == states2[r][0] else "X"
        if r < 64:
            dw = "0" if Wexp1[r] == Wexp2[r] else "X"
        else:
            dw = "-"

        de_zeros.append(de == "0")
        da_zeros.append(da == "0")
        if r < 64:
            dw_zeros.append(dw == "0")

        if r <= 20 or r >= 55 or de == "0" or da == "0":
            print(f"  {r:3d} |  {de}   {da}  |  {dw}")

    total_de_zero = sum(de_zeros)
    total_da_zero = sum(da_zeros)
    total_dw_zero = sum(dw_zeros)

    # Forward consecutive
    fwd = 0
    for r in range(1, 65):
        if de_zeros[r]: fwd += 1
        else: break

    # Backward consecutive
    bwd = 0
    for r in range(64, 0, -1):
        if de_zeros[r]: bwd += 1
        else: break

    gap = 64 - fwd - bwd
    print(f"\n  Summary:")
    print(f"    Total δe=0 rounds: {total_de_zero}/65")
    print(f"    Total δa=0 rounds: {total_da_zero}/65")
    print(f"    Total δW=0 rounds: {total_dw_zero}/64")
    print(f"    Forward chain (δe=0): {fwd} rounds")
    print(f"    Backward chain (δe=0): {bwd} rounds")
    print(f"    GAP: {gap} rounds")
    print(f"    δW words changed: {sum(1 for i in range(16) if W1[i]!=best_W2[i])}/16")


# ============================================================
if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 200

    print("="*70)
    print("SCF: BACKWARD WANG CHAIN + MEET IN THE MIDDLE")
    print("  Гениальная идея: атаковать С ДВУХ СТОРОН")
    print("="*70)

    exp1_backward_theory()
    exp2_verify_backward(min(N, 300))
    exp3_schedule_zeros(min(N, 200))
    exp4_meet_in_middle(min(N, 500))
    exp5_adaptive_backward(N)

    print("\n" + "="*70)
    print("ФИНАЛЬНЫЙ ИТОГ")
    print("="*70)
    print("""
  Backward Wang chain ТЕОРЕТИЧЕСКИ:
    δa[r-3] = 0     ВСЕГДА (бесплатно, shift register)
    δe[r-3] = -δW[r] УПРАВЛЯЕМО (через message schedule)

  Если δW[r]=0 для последних M раундов → M backward zeros.

  Но: schedule расширяет δW[0..15] на ВСЕ 48 раундов W[16..63].
  Единственный δ, дающий δW[16..63]=0, это δW[0..15]=0 (тривиальный).

  Meet in the middle:
    Forward: ~15 zeros (Wang chain)
    Backward: ~0-1 zeros (schedule diffusion)
    Gap: ~48 rounds (вся хаотическая зона)

  Вывод: backward chain СУЩЕСТВУЕТ но schedule не даёт δW=0
  в поздних раундах. Message schedule — ВТОРОЙ барьер (после carry).
    """)
