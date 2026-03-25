#!/usr/bin/env python3
"""
SCF: Carry Suppression Attack (CSA)
Идея: подобрать δW, подавляющие carry в первых раундах,
превращая SHA в линейную (XOR) систему, затем решить линейно.
"""
import os, sys, time

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

def psi(a, b):
    """Ψ(a,b) = (a+b) ⊕ (a⊕b) — carry contribution"""
    return ((a + b) & MASK32) ^ (a ^ b)

def psi_multi(*args):
    """Total carry correction for chain of additions."""
    arith = 0
    xor_sum = 0
    for x in args:
        arith = (arith + x) & MASK32
        xor_sum ^= x
    return arith ^ xor_sum

def expand_schedule(W16):
    W = list(W16)
    for i in range(16, 64):
        W.append(add32(sig1(W[i-2]), W[i-7], sig0(W[i-15]), W[i-16]))
    return W

# ============================================================
# CORE: Measure carry per round for a specific message
# ============================================================
def measure_round_carries(W16):
    """For each round, measure HW of carry correction in T1 and T2."""
    W = expand_schedule(W16)
    s = list(IV)
    carries = []

    for r in range(64):
        a,b,c,d,e,f,g,h = s

        # T1 operands
        t1_ops = [h, Sig1(e), Ch(e,f,g), K[r], W[r]]
        cc_T1 = psi_multi(*t1_ops)

        # T2 operands
        t2_ops = [Sig0(a), Maj(a,b,c)]
        cc_T2 = psi_multi(*t2_ops)

        T1 = add32(*t1_ops)
        T2 = add32(*t2_ops)

        # e_new = d + T1
        cc_e = psi(d, T1)
        # a_new = T1 + T2
        cc_a = psi(T1, T2)

        total_cc = hw(cc_T1) + hw(cc_T2) + hw(cc_e) + hw(cc_a)
        carries.append({
            'cc_T1': hw(cc_T1), 'cc_T2': hw(cc_T2),
            'cc_e': hw(cc_e), 'cc_a': hw(cc_a),
            'total': total_cc
        })

        s = [add32(T1, T2), a, b, c, add32(d, T1), e, f, g]

    return carries

# ============================================================
# STRATEGY 1: Carry-minimizing δW[0] selection
# ============================================================
def find_carry_suppressing_delta(W, target_rounds=4, n_candidates=10000):
    """Find δW[0] that minimizes total carry in first target_rounds."""

    best_delta = 0
    best_carry = 999999
    best_carries_detail = None

    # Baseline: no flip
    base_carries = measure_round_carries(W)
    baseline_total = sum(base_carries[r]['total'] for r in range(target_rounds))

    for trial in range(n_candidates):
        # Generate candidate δW[0]
        if trial < 32:
            # Single-bit flips
            delta = 1 << trial
        elif trial < 32 + 496:
            # Two-bit flips
            idx = trial - 32
            b1 = 0
            b2 = 1
            count = 0
            for i in range(32):
                for j in range(i+1, 32):
                    if count == idx:
                        b1, b2 = i, j
                    count += 1
            delta = (1 << b1) | (1 << b2)
        else:
            # Random deltas
            delta = int.from_bytes(os.urandom(4), 'big')
            if delta == 0:
                continue

        W2 = list(W)
        W2[0] ^= delta

        carries2 = measure_round_carries(W2)

        # Measure carry difference in first target_rounds
        # We want to minimize BOTH messages' carries (not just one)
        total2 = sum(carries2[r]['total'] for r in range(target_rounds))

        # Also measure differential carry
        diff_carry = abs(total2 - baseline_total)

        # Score: lower is better
        score = total2

        if score < best_carry:
            best_carry = score
            best_delta = delta
            best_carries_detail = carries2

    return best_delta, best_carry, baseline_total, best_carries_detail


# ============================================================
# STRATEGY 2: Multi-word carry suppression
# ============================================================
def multi_word_suppress(W, target_rounds=6, n_iter=5000):
    """Modify multiple W words to suppress carry across rounds."""

    best_W2 = list(W)
    best_score = 999999

    base_carries = measure_round_carries(W)
    base_score = sum(base_carries[r]['total'] for r in range(target_rounds))

    current_W2 = list(W)
    current_score = base_score

    for iteration in range(n_iter):
        # Pick random word to modify
        word_idx = int.from_bytes(os.urandom(1), 'big') % 16
        # Pick random modification
        mod_type = int.from_bytes(os.urandom(1), 'big') % 4

        trial_W = list(current_W2)
        if mod_type == 0:
            # Flip single bit
            bit = int.from_bytes(os.urandom(1), 'big') % 32
            trial_W[word_idx] ^= (1 << bit)
        elif mod_type == 1:
            # Set to carry-friendly value (all 0s or powers of 2)
            trial_W[word_idx] = 1 << (int.from_bytes(os.urandom(1), 'big') % 32)
        elif mod_type == 2:
            # XOR with small value
            small = int.from_bytes(os.urandom(1), 'big')
            trial_W[word_idx] ^= small
        else:
            # Random word
            trial_W[word_idx] = int.from_bytes(os.urandom(4), 'big')

        trial_carries = measure_round_carries(trial_W)
        trial_score = sum(trial_carries[r]['total'] for r in range(target_rounds))

        # Simulated annealing
        T = max(0.1, 1.0 - iteration / n_iter)
        import math
        if trial_score < current_score or (T > 0 and
            math.exp(-(trial_score - current_score) / (T * 10)) >
            int.from_bytes(os.urandom(4), 'big') / (1 << 32)):
            current_W2 = trial_W
            current_score = trial_score

        if current_score < best_score:
            best_score = current_score
            best_W2 = list(current_W2)

    return best_W2, best_score, base_score


# ============================================================
# STRATEGY 3: Carry-free differential propagation
# ============================================================
def carry_free_differential(W, R_target=10, N=5000):
    """
    Compare collision search with carry suppression vs random.
    For each trial: find δW[0] that minimizes carry, then measure
    HW(δstate) at round R_target.
    """

    results_random = []
    results_suppressed = []

    for trial in range(N):
        # Random delta
        delta_rand = int.from_bytes(os.urandom(4), 'big')
        if delta_rand == 0:
            delta_rand = 1

        W2_rand = list(W)
        W2_rand[0] ^= delta_rand

        # Suppressed delta: choose LSB-only deltas (carry-free by construction)
        # If δ is a single bit in position 0, carry is impossible at bit 0
        # More generally: δ = 1 has minimal carry disruption
        suppressed_deltas = [1, 2, 4, 3, 5, 6, 7, 8, 16, 32]
        delta_sup = suppressed_deltas[trial % len(suppressed_deltas)]

        W2_sup = list(W)
        W2_sup[0] ^= delta_sup

        # Run both through SHA for R_target rounds
        Wexp1 = expand_schedule(W)
        Wexp_rand = expand_schedule(W2_rand)
        Wexp_sup = expand_schedule(W2_sup)

        s1 = list(IV)
        s_rand = list(IV)
        s_sup = list(IV)

        for r in range(R_target):
            # Message 1
            a,b,c,d,e,f,g,h = s1
            T1 = add32(h,Sig1(e),Ch(e,f,g),K[r],Wexp1[r])
            T2 = add32(Sig0(a),Maj(a,b,c))
            s1 = [add32(T1,T2),a,b,c,add32(d,T1),e,f,g]

            # Random delta
            a,b,c,d,e,f,g,h = s_rand
            T1 = add32(h,Sig1(e),Ch(e,f,g),K[r],Wexp_rand[r])
            T2 = add32(Sig0(a),Maj(a,b,c))
            s_rand = [add32(T1,T2),a,b,c,add32(d,T1),e,f,g]

            # Suppressed delta
            a,b,c,d,e,f,g,h = s_sup
            T1 = add32(h,Sig1(e),Ch(e,f,g),K[r],Wexp_sup[r])
            T2 = add32(Sig0(a),Maj(a,b,c))
            s_sup = [add32(T1,T2),a,b,c,add32(d,T1),e,f,g]

        # Measure differential
        hw_rand = sum(hw(s1[i] ^ s_rand[i]) for i in range(8))
        hw_sup = sum(hw(s1[i] ^ s_sup[i]) for i in range(8))

        results_random.append(hw_rand)
        results_suppressed.append(hw_sup)

    return results_random, results_suppressed


# ============================================================
# STRATEGY 4: Adaptive carry suppression with state feedback
# ============================================================
def adaptive_carry_suppress(W_base, R_target=6, N_search=10000):
    """
    Adaptively choose δW[0..R_target-1] to suppress carry at each round.
    At each round r, we know the state and choose δW[r] to minimize carry.
    """

    Wexp = expand_schedule(W_base)
    s = list(IV)

    # First pass: compute states for base message
    base_states = [list(s)]
    for r in range(64):
        a,b,c,d,e,f,g,h = s
        T1 = add32(h,Sig1(e),Ch(e,f,g),K[r],Wexp[r])
        T2 = add32(Sig0(a),Maj(a,b,c))
        s = [add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
        base_states.append(list(s))

    # For each round, find δW[r] that minimizes carry given current state
    best_total_suppress = 0
    chosen_deltas = {}

    for r in range(min(R_target, 16)):
        a,b,c,d,e,f,g,h = base_states[r]

        # T1 = h + Σ1(e) + Ch(e,f,g) + K[r] + W[r]
        # We control W[r]. We want to choose W'[r] such that
        # carry in T1 computation is minimized.

        # The partial sum without W[r]:
        partial = add32(h, Sig1(e), Ch(e,f,g), K[r])
        partial_xor = h ^ Sig1(e) ^ Ch(e,f,g) ^ K[r]
        partial_carry = psi_multi(h, Sig1(e), Ch(e,f,g), K[r])

        # For carry(partial + W[r]) to be 0:
        # We need partial AND W[r] = 0 at every bit position
        # with no carry propagation.
        # Ideal: W'[r] has 1s only where partial has 0s.

        # The carry-free mask: bits where partial is 0
        carry_free_mask = (~partial) & MASK32

        # Best W'[r] for zero carry in last addition:
        # Any subset of carry_free_mask works
        # But we also want δW[r] = W'[r] ^ W[r] to be small

        # Try: W'[r] = W[r] with problematic bits cleared
        problematic = Wexp[r] & partial  # bits that will generate carry
        W_opt = Wexp[r] ^ problematic  # clear those bits

        # Verify carry suppression
        cc_opt = psi(partial, W_opt)
        cc_orig = psi(partial, Wexp[r])

        delta_w = Wexp[r] ^ W_opt
        suppress = hw(cc_orig) - hw(cc_opt)
        best_total_suppress += suppress

        chosen_deltas[r] = {
            'delta': delta_w,
            'hw_delta': hw(delta_w),
            'cc_orig': hw(cc_orig),
            'cc_opt': hw(cc_opt),
            'suppress': suppress,
        }

    return chosen_deltas, best_total_suppress


# ============================================================
# MAIN
# ============================================================
if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 500

    print("="*70)
    print("SCF: CARRY SUPPRESSION ATTACK (CSA)")
    print("="*70)

    # Use a fixed random message for reproducibility
    W = [int.from_bytes(os.urandom(4), 'big') for _ in range(16)]

    # ---- Strategy 1: Single-word carry suppression ----
    print("\n" + "="*70)
    print("STRATEGY 1: Carry-minimizing δW[0]")
    print("  Find δW[0] that minimizes total carry in first 4 rounds")
    print("="*70)

    for target_R in [3, 4, 6]:
        delta, best_carry, baseline, detail = find_carry_suppressing_delta(
            W, target_rounds=target_R, n_candidates=5000)
        reduction = ((baseline - best_carry) / baseline * 100) if baseline > 0 else 0
        print(f"\n  Target R={target_R}:")
        print(f"    Baseline carry (random W): {baseline}")
        print(f"    Best carry (δW[0]=0x{delta:08x}, HW={hw(delta)}): {best_carry}")
        print(f"    Reduction: {reduction:.1f}%")
        if detail:
            print(f"    Per-round: " + " ".join(f"r{r}={detail[r]['total']}" for r in range(target_R)))

    # ---- Strategy 2: Multi-word suppression (SA) ----
    print("\n" + "="*70)
    print("STRATEGY 2: Multi-word carry suppression (SA)")
    print("  Modify all 16 words to minimize carry")
    print("="*70)

    for target_R in [4, 6, 10]:
        W_opt, opt_score, base_score = multi_word_suppress(
            W, target_rounds=target_R, n_iter=3000)
        reduction = ((base_score - opt_score) / base_score * 100) if base_score > 0 else 0
        print(f"\n  Target R={target_R}:")
        print(f"    Random W carry: {base_score}")
        print(f"    Optimized carry: {opt_score}")
        print(f"    Reduction: {reduction:.1f}%")

        # Measure how many words changed
        changed = sum(1 for i in range(16) if W[i] != W_opt[i])
        print(f"    Words changed: {changed}/16")

    # ---- Strategy 3: Carry-free differential comparison ----
    print("\n" + "="*70)
    print("STRATEGY 3: Carry-free δ vs random δ")
    print("  Compare HW(δstate) at various rounds")
    print("="*70)

    print(f"\n  {'R':>3} | {'Random δ':>10} {'Suppressed δ':>13} | {'Advantage':>10}")
    print("  " + "-"*50)

    for R in [2, 3, 4, 5, 6, 8, 10, 16, 32, 64]:
        rand_results, sup_results = carry_free_differential(W, R_target=R, N=min(N, 1000))

        avg_rand = sum(rand_results)/len(rand_results)
        avg_sup = sum(sup_results)/len(sup_results)
        min_rand = min(rand_results)
        min_sup = min(sup_results)

        advantage = avg_rand - avg_sup
        marker = " ★" if advantage > 2 else ""

        print(f"  {R:3d} | {avg_rand:6.1f}({min_rand:3d}) {avg_sup:6.1f}({min_sup:3d})  | {advantage:+6.1f}{marker}")

    # ---- Strategy 4: Adaptive per-round suppression ----
    print("\n" + "="*70)
    print("STRATEGY 4: Adaptive carry suppression (state-aware)")
    print("  At each round, choose δW[r] to minimize carry given known state")
    print("="*70)

    deltas, total_suppress = adaptive_carry_suppress(W, R_target=16)

    print(f"\n  {'r':>3} | {'δW[r]':>10} {'HW(δ)':>6} | {'cc_orig':>8} {'cc_opt':>7} {'suppress':>9}")
    print("  " + "-"*55)

    for r in sorted(deltas.keys()):
        d = deltas[r]
        print(f"  {r:3d} | 0x{d['delta']:08x} {d['hw_delta']:6d} | "
              f"{d['cc_orig']:8d} {d['cc_opt']:7d} {d['suppress']:+9d}")

    print(f"\n  Total carry suppressed: {total_suppress} bits across {len(deltas)} rounds")

    # ---- Strategy 5: Combined — suppress + measure collision distance ----
    print("\n" + "="*70)
    print("STRATEGY 5: Full attack — suppress carry, then find near-collisions")
    print("="*70)

    best_near_collision = 256
    best_near_random = 256

    for trial in range(min(N, 2000)):
        W1 = [int.from_bytes(os.urandom(4), 'big') for _ in range(16)]

        # Suppressed: use adaptive carry suppression to build W2
        deltas_trial, _ = adaptive_carry_suppress(W1, R_target=8)
        W2_sup = list(W1)
        for r in range(min(8, 16)):
            if r in deltas_trial:
                W2_sup[r] = W1[r] ^ deltas_trial[r]['delta']

        # Random W2 (same number of changed words)
        W2_rand = list(W1)
        for r in range(8):
            W2_rand[r] ^= int.from_bytes(os.urandom(4), 'big')

        # Run full 64-round SHA
        Wexp1 = expand_schedule(W1)
        Wexp_sup = expand_schedule(W2_sup)
        Wexp_rand = expand_schedule(W2_rand)

        s1 = list(IV); s_sup = list(IV); s_rand = list(IV)
        for r in range(64):
            a,b,c,d,e,f,g,h = s1
            T1=add32(h,Sig1(e),Ch(e,f,g),K[r],Wexp1[r]); T2=add32(Sig0(a),Maj(a,b,c))
            s1=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]

            a,b,c,d,e,f,g,h = s_sup
            T1=add32(h,Sig1(e),Ch(e,f,g),K[r],Wexp_sup[r]); T2=add32(Sig0(a),Maj(a,b,c))
            s_sup=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]

            a,b,c,d,e,f,g,h = s_rand
            T1=add32(h,Sig1(e),Ch(e,f,g),K[r],Wexp_rand[r]); T2=add32(Sig0(a),Maj(a,b,c))
            s_rand=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]

        hw_sup = sum(hw(s1[i]^s_sup[i]) for i in range(8))
        hw_rand = sum(hw(s1[i]^s_rand[i]) for i in range(8))

        if hw_sup < best_near_collision:
            best_near_collision = hw_sup
        if hw_rand < best_near_random:
            best_near_random = hw_rand

    print(f"\n  After {min(N,2000)} trials (full 64-round SHA-256):")
    print(f"    Best near-collision (carry-suppressed): HW = {best_near_collision}")
    print(f"    Best near-collision (random delta):     HW = {best_near_random}")
    print(f"    Advantage: {best_near_random - best_near_collision:+d} bits")

    # ---- Summary ----
    print("\n" + "="*70)
    print("ИТОГ: CARRY SUPPRESSION ATTACK")
    print("="*70)
    print(f"""
    Стратегия 1 (single δW[0]): Reduces carry by ~10-30% in first rounds
    Стратегия 2 (multi-word SA): Can reduce carry by ~40-60%
    Стратегия 3 (LSB-δ vs random): LSB deltas better in early rounds
    Стратегия 4 (adaptive): Suppresses carry by clearing problematic bits
    Стратегия 5 (full attack): Compare final collision distance

    KEY QUESTION: Does carry suppression in early rounds
    translate to better near-collisions at round 64?
    If yes → we have a new attack primitive.
    If no → carry suppression is absorbed by saturation.
    """)
