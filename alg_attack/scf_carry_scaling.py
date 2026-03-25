#!/usr/bin/env python3
"""
SCF: Carry Suppression Scaling — можно ли подавить carry через ВСЕ 64 раунда?

Вопрос: если Strategy 4 даёт cc=0 для r=0..15 (свободные W),
можно ли выбрать W[0..15] так, чтобы schedule W[16..63]
тоже давал подавление carry?

Три подхода:
1. Измерить: как carry в r=16..63 зависит от W[0..15]?
2. SA-оптимизация: минимизировать carry по ВСЕМ 64 раундам
3. Гибрид: обнулить carry в r=0..15 (adaptive), затем оценить r=16..63
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

def psi_multi(*args):
    arith, xor_s = 0, 0
    for x in args:
        arith = (arith+x)&MASK32
        xor_s ^= x
    return arith ^ xor_s

def psi(a,b): return ((a+b)&MASK32)^(a^b)

def expand_schedule(W16):
    W=list(W16)
    for i in range(16,64):
        W.append(add32(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    return W

def total_carry_profile(W16):
    """Compute carry HW at each of 64 rounds."""
    W = expand_schedule(W16)
    s = list(IV)
    profile = []
    for r in range(64):
        a,b,c,d,e,f,g,h = s
        cc_T1 = hw(psi_multi(h, Sig1(e), Ch(e,f,g), K[r], W[r]))
        cc_T2 = hw(psi_multi(Sig0(a), Maj(a,b,c)))
        T1 = add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2 = add32(Sig0(a),Maj(a,b,c))
        cc_e = hw(psi(d, T1))
        cc_a = hw(psi(T1, T2))
        profile.append(cc_T1 + cc_T2 + cc_e + cc_a)
        s = [add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
    return profile

def total_carry_score(W16, weight_late=1.0):
    """Weighted total carry across all rounds.
    weight_late > 1 emphasizes later rounds."""
    profile = total_carry_profile(W16)
    score = 0
    for r in range(64):
        w = 1.0 if r < 16 else weight_late
        score += profile[r] * w
    return score

# ============================================================
# EXP 1: Baseline — how much carry varies across random messages?
# ============================================================
def exp1_baseline(N):
    print("="*70)
    print("EXP 1: BASELINE — carry variance across random messages")
    print("="*70)

    profiles = []
    total_scores = []
    for _ in range(N):
        W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        p = total_carry_profile(W)
        profiles.append(p)
        total_scores.append(sum(p))

    avg_total = sum(total_scores)/N
    min_total = min(total_scores)
    max_total = max(total_scores)
    std_total = (sum((x-avg_total)**2 for x in total_scores)/N)**0.5

    print(f"\n  Total carry across 64 rounds (N={N}):")
    print(f"    Mean:  {avg_total:.1f}")
    print(f"    Std:   {std_total:.1f}")
    print(f"    Min:   {min_total}")
    print(f"    Max:   {max_total}")
    print(f"    Range: {max_total - min_total}")

    # Per-phase averages
    phases = [(0,16,"r=0..15 (free W)"), (16,32,"r=16..31"), (32,48,"r=32..47"), (48,64,"r=48..63")]
    print(f"\n  Per-phase carry:")
    for r_s, r_e, name in phases:
        avg = sum(sum(p[r_s:r_e]) for p in profiles) / N
        mn = min(sum(p[r_s:r_e]) for p in profiles)
        print(f"    {name}: avg={avg:.1f}, min={mn}")

    # Key: what is the natural carry per round?
    per_round_avg = [sum(p[r] for p in profiles)/N for r in range(64)]
    print(f"\n  Per-round carry (avg):")
    print(f"  {'r':>3} | carry")
    for r in range(0, 64, 4):
        print(f"  {r:3d} | {per_round_avg[r]:.1f}")

    return avg_total, min_total, per_round_avg

# ============================================================
# EXP 2: SA optimization — minimize carry across ALL 64 rounds
# ============================================================
def exp2_sa_global(N_iter, weight_late=2.0):
    print("\n" + "="*70)
    print(f"EXP 2: SA — minimize carry across ALL 64 rounds (weight_late={weight_late})")
    print("="*70)

    best_W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    best_score = total_carry_score(best_W, weight_late)
    current_W = list(best_W)
    current_score = best_score
    init_score = best_score

    history = []

    for it in range(N_iter):
        trial_W = list(current_W)

        # Mutation
        r_type = int.from_bytes(os.urandom(1),'big') % 5
        word_idx = int.from_bytes(os.urandom(1),'big') % 16

        if r_type == 0:
            # Flip single bit
            bit = int.from_bytes(os.urandom(1),'big') % 32
            trial_W[word_idx] ^= (1 << bit)
        elif r_type == 1:
            # Small XOR
            trial_W[word_idx] ^= int.from_bytes(os.urandom(1),'big')
        elif r_type == 2:
            # Set to 0 (minimal carry contribution)
            trial_W[word_idx] = 0
        elif r_type == 3:
            # Set to power of 2
            trial_W[word_idx] = 1 << (int.from_bytes(os.urandom(1),'big') % 32)
        else:
            # Random
            trial_W[word_idx] = int.from_bytes(os.urandom(4),'big')

        trial_score = total_carry_score(trial_W, weight_late)

        T = max(0.01, 1.0 - it/N_iter)
        if trial_score < current_score or \
           math.exp(-(trial_score - current_score)/(T*20)) > \
           int.from_bytes(os.urandom(4),'big')/(1<<32):
            current_W = trial_W
            current_score = trial_score

        if current_score < best_score:
            best_score = current_score
            best_W = list(current_W)

        if it % (N_iter//10) == 0:
            history.append((it, current_score, best_score))

    print(f"\n  Initial score: {init_score:.1f}")
    print(f"  Final score:   {best_score:.1f}")
    print(f"  Reduction:     {(init_score-best_score)/init_score*100:.1f}%")

    # Show optimized carry profile
    opt_profile = total_carry_profile(best_W)
    rand_profile_avg = [0]*64

    # Compare with average random
    for _ in range(200):
        W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        p = total_carry_profile(W)
        for r in range(64):
            rand_profile_avg[r] += p[r]
    for r in range(64):
        rand_profile_avg[r] /= 200

    print(f"\n  {'r':>3} | {'optimized':>9} {'random':>7} {'reduction':>10}")
    print("  " + "-"*40)
    for r in range(0, 64, 4):
        opt_phase = sum(opt_profile[r:r+4])
        rand_phase = sum(rand_profile_avg[r:r+4])
        red = (rand_phase-opt_phase)/rand_phase*100 if rand_phase>0 else 0
        bar = "█" * max(0, int(red/5))
        print(f"  {r:3d} | {opt_phase:9.1f} {rand_phase:7.1f} {red:+9.1f}% {bar}")

    return best_W, best_score, opt_profile

# ============================================================
# EXP 3: Schedule carry — how much carry is in W[16..63] itself?
# ============================================================
def exp3_schedule_carry(N):
    print("\n" + "="*70)
    print("EXP 3: Schedule carry — carry WITHIN message expansion")
    print("  W[i] = σ1(W[i-2]) + W[i-7] + σ0(W[i-15]) + W[i-16]")
    print("  How much carry is in this addition?")
    print("="*70)

    schedule_carry = [0.0]*64
    for _ in range(N):
        W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        Wexp = expand_schedule(W)
        for i in range(16, 64):
            ops = [sig1(Wexp[i-2]), Wexp[i-7], sig0(Wexp[i-15]), Wexp[i-16]]
            cc = psi_multi(*ops)
            schedule_carry[i] += hw(cc)

    for i in range(16, 64):
        schedule_carry[i] /= N

    print(f"\n  Schedule carry per expanded word (expected ~16 for random):")
    print(f"  {'i':>3} | carry_hw | bar")
    print("  " + "-"*40)
    for i in range(16, 64, 2):
        bar = "█" * int(schedule_carry[i])
        print(f"  {i:3d} | {schedule_carry[i]:8.1f} | {bar}")

    avg = sum(schedule_carry[16:])/48
    print(f"\n  Average schedule carry: {avg:.1f} (expected 16 for random add)")
    print(f"  Schedule carry is {'controllable' if avg < 14 else 'essentially random'}")

    return schedule_carry

# ============================================================
# EXP 4: Can we reduce schedule carry by choosing W[0..15]?
# ============================================================
def exp4_schedule_optimize(N_iter=10000):
    print("\n" + "="*70)
    print("EXP 4: Optimize W[0..15] to minimize SCHEDULE carry (W[16..63])")
    print("="*70)

    def schedule_carry_total(W16):
        Wexp = expand_schedule(W16)
        total = 0
        for i in range(16, 64):
            ops = [sig1(Wexp[i-2]), Wexp[i-7], sig0(Wexp[i-15]), Wexp[i-16]]
            total += hw(psi_multi(*ops))
        return total

    best_W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    best_sc = schedule_carry_total(best_W)
    current_W = list(best_W)
    current_sc = best_sc
    init_sc = best_sc

    for it in range(N_iter):
        trial_W = list(current_W)
        word_idx = int.from_bytes(os.urandom(1),'big') % 16
        r_type = int.from_bytes(os.urandom(1),'big') % 4

        if r_type == 0:
            trial_W[word_idx] ^= (1 << (int.from_bytes(os.urandom(1),'big') % 32))
        elif r_type == 1:
            trial_W[word_idx] = 0
        elif r_type == 2:
            trial_W[word_idx] ^= int.from_bytes(os.urandom(2),'big')
        else:
            trial_W[word_idx] = int.from_bytes(os.urandom(4),'big')

        trial_sc = schedule_carry_total(trial_W)

        T = max(0.01, 1.0 - it/N_iter)
        if trial_sc < current_sc or \
           math.exp(-(trial_sc-current_sc)/(T*10)) > \
           int.from_bytes(os.urandom(4),'big')/(1<<32):
            current_W = trial_W
            current_sc = trial_sc

        if current_sc < best_sc:
            best_sc = current_sc
            best_W = list(current_W)

    # Random baseline
    rand_scores = []
    for _ in range(500):
        W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        rand_scores.append(schedule_carry_total(W))
    rand_avg = sum(rand_scores)/len(rand_scores)
    rand_min = min(rand_scores)

    print(f"\n  Random schedule carry: avg={rand_avg:.1f}, min={rand_min}")
    print(f"  Optimized schedule carry: {best_sc}")
    print(f"  Reduction from avg: {(rand_avg-best_sc)/rand_avg*100:.1f}%")
    print(f"  Reduction from min: {(rand_min-best_sc)/rand_min*100:.1f}%")

    # Show the optimized W
    print(f"\n  Optimized W[0..15]:")
    for i in range(16):
        print(f"    W[{i:2d}] = 0x{best_W[i]:08x} (HW={hw(best_W[i]):2d})")

    return best_W, best_sc

# ============================================================
# EXP 5: The key experiment — combined compression + schedule optimization
# ============================================================
def exp5_combined(N_iter=15000):
    print("\n" + "="*70)
    print("EXP 5: COMBINED — minimize total carry (compression + schedule)")
    print("  This is the FULL carry suppression scaling experiment.")
    print("="*70)

    def combined_score(W16):
        """Total carry across all 64 rounds of SHA-256."""
        return sum(total_carry_profile(W16))

    best_W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    best_sc = combined_score(best_W)
    current_W = list(best_W)
    current_sc = best_sc
    init_sc = best_sc

    for it in range(N_iter):
        trial_W = list(current_W)
        word_idx = int.from_bytes(os.urandom(1),'big') % 16
        r_type = int.from_bytes(os.urandom(1),'big') % 5

        if r_type == 0:
            trial_W[word_idx] ^= (1 << (int.from_bytes(os.urandom(1),'big') % 32))
        elif r_type == 1:
            trial_W[word_idx] = 0
        elif r_type == 2:
            trial_W[word_idx] ^= int.from_bytes(os.urandom(2),'big')
        elif r_type == 3:
            trial_W[word_idx] = int.from_bytes(os.urandom(4),'big')
        else:
            # Swap two words
            other = int.from_bytes(os.urandom(1),'big') % 16
            trial_W[word_idx], trial_W[other] = trial_W[other], trial_W[word_idx]

        trial_sc = combined_score(trial_W)

        T = max(0.01, 1.0 - it/N_iter)
        if trial_sc < current_sc or \
           math.exp(-(trial_sc-current_sc)/(T*15)) > \
           int.from_bytes(os.urandom(4),'big')/(1<<32):
            current_W = trial_W
            current_sc = trial_sc

        if current_sc < best_sc:
            best_sc = current_sc
            best_W = list(current_W)

    # Baseline
    rand_scores = []
    for _ in range(500):
        W = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        rand_scores.append(combined_score(W))
    rand_avg = sum(rand_scores)/len(rand_scores)
    rand_min = min(rand_scores)

    opt_profile = total_carry_profile(best_W)

    print(f"\n  Random total carry: avg={rand_avg:.1f}, min={rand_min}")
    print(f"  Optimized total carry: {best_sc}")
    print(f"  Reduction from avg: {(rand_avg-best_sc)/rand_avg*100:.1f}%")
    print(f"  Reduction from random min: {(rand_min-best_sc)/rand_min*100:.1f}%")

    # Per-phase
    print(f"\n  Per-phase profile (optimized vs random avg):")
    phases = [(0,16), (16,32), (32,48), (48,64)]
    for rs, re in phases:
        opt_phase = sum(opt_profile[rs:re])
        rand_phase = sum(sum(total_carry_profile(
            [int.from_bytes(os.urandom(4),'big') for _ in range(16)])[rs:re]
        ) for _ in range(100)) / 100
        red = (rand_phase-opt_phase)/rand_phase*100 if rand_phase > 0 else 0
        print(f"    r={rs:2d}..{re-1:2d}: opt={opt_phase:5.0f} rand={rand_phase:5.0f} reduction={red:+.1f}%")

    # Per-round detail for optimized
    print(f"\n  Optimized carry per round:")
    print(f"  {'r':>3} | carry | bar")
    for r in range(0, 64, 2):
        c = opt_profile[r]
        bar = "█" * min(c, 30)
        print(f"  {r:3d} | {c:5d} | {bar}")

    # How many zero-carry rounds?
    zero_rounds = sum(1 for c in opt_profile if c == 0)
    low_rounds = sum(1 for c in opt_profile if c <= 5)
    print(f"\n  Zero-carry rounds: {zero_rounds}/64")
    print(f"  Low-carry (≤5) rounds: {low_rounds}/64")

    return best_W, best_sc, opt_profile


# ============================================================
# MAIN
# ============================================================
if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 500

    avg_total, min_total, per_round = exp1_baseline(min(N, 500))
    best_W_sa, score_sa, profile_sa = exp2_sa_global(N_iter=min(N*20, 15000))
    schedule_carry = exp3_schedule_carry(min(N, 500))
    best_W_sch, sc_sch = exp4_schedule_optimize(N_iter=min(N*20, 15000))
    best_W_comb, sc_comb, profile_comb = exp5_combined(N_iter=min(N*30, 20000))

    print("\n" + "="*70)
    print("ФИНАЛЬНЫЙ ИТОГ: CARRY SCALING")
    print("="*70)
    print(f"""
  Random carry (64 rounds):           avg={avg_total:.0f}, min={min_total}
  SA-optimized (compression focus):   {score_sa:.0f}
  Schedule-optimized:                 {sc_sch}
  COMBINED (full optimization):       {sc_comb}

  Reduction from random avg:          {(avg_total-sc_comb)/avg_total*100:.1f}%
  Reduction from random min:          {(min_total-sc_comb)/min_total*100:.1f}%

  ВЫВОД: Если reduction > 30% — carry suppression масштабируется.
         Если reduction ~ 10-15% — это шум, не масштабируется.
         Если zero-carry rounds > 5 — есть структурный контроль.
    """)
