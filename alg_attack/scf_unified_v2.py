#!/usr/bin/env python3
"""
SCF: UNIFIED SOLVER v2 — совместная оптимизация ВСЕХ раундов.

v1 проиграл SA потому что Wang и Gap КОНФЛИКТОВАЛИ.
v2: вычисляем required W2[r] для КАЖДОГО раунда, затем ищем
W2[0..15] минимизирующий ВЗВЕШЕННУЮ сумму mismatch.

Ключевое отличие: НЕ ищем случайно.
Для каждого слова w: вычисляем ОПТИМАЛЬНУЮ коррекцию через
schedule sensitivity и required values. Это ВЫЧИСЛЕНИЕ, не поиск.

Архитектура v2:
  1. Compute required_W[r] для r=0..63 (exact inversion per round)
  2. Compute schedule_sensitivity[r][w] = ∂W_exp[r]/∂W[w]
  3. Weighted least squares: найти W[0..15] минимизирующий
     Σ weight[r] × |required[r] - schedule(W)[r]|
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

def sha_compress(W16):
    st, _ = sha_all_states(W16)
    return [add32(IV[i], st[64][i]) for i in range(8)]

def hash_hw(W1, W2):
    H1=sha_compress(W1); H2=sha_compress(W2)
    return sum(hw(H1[i]^H2[i]) for i in range(8))


def compute_all_required(W1, W2):
    """For each round r, compute required W2_expanded[r] for δe[r+1]=0."""
    st1, Wexp1 = sha_all_states(W1)
    st2, Wexp2 = sha_all_states(W2)

    required = {}
    for r in range(64):
        s1r = st1[r]; s2r = st2[r]
        e1_target = st1[r+1][4]
        d2 = s2r[3]
        T1_needed = sub32(e1_target, d2)
        W2_needed = sub32(sub32(sub32(sub32(T1_needed,
                    s2r[7]), Sig1(s2r[4])), Ch(s2r[4],s2r[5],s2r[6])), K[r])
        required[r] = W2_needed

    return required, Wexp2


def weighted_mismatch(W1, W2, weights=None):
    """Weighted sum of HW(required[r] ^ actual[r])."""
    required, actual_exp = compute_all_required(W1, W2)
    total = 0
    for r in range(64):
        w = weights[r] if weights else 1.0
        total += w * hw(required[r] ^ actual_exp[r])
    return total


def compute_schedule_jacobian(W2, rounds):
    """For each (round, word): ∂W_expanded[round]/∂W[word].
    Returns dict: (round, word) → sensitivity (32-bit int)."""
    Wexp_base = expand_real(W2)
    jac = {}
    for w in range(16):
        W_test = list(W2)
        W_test[w] = add32(W2[w], 1)
        Wexp_test = expand_real(W_test)
        for r in rounds:
            jac[(r, w)] = sub32(Wexp_test[r], Wexp_base[r])
    return jac


def v2_solve(W1, delta_w0=0x80000000, max_outer=10):
    """
    v2 solver: iterative weighted schedule correction.

    Each iteration:
    1. Compute required W2[r] for all rounds (exact inversion)
    2. Compute schedule Jacobian
    3. For each word w: compute weighted optimal correction
    4. Apply best correction
    5. Repeat
    """
    W2 = list(W1)
    W2[0] ^= delta_w0

    # Weights: higher for rounds close to hash output (more impact)
    # Also higher for Wang region (want to preserve those zeros)
    weights = {}
    for r in range(64):
        if r <= 1:
            weights[r] = 0.1
        elif r <= 16:
            weights[r] = 3.0    # Wang region — HIGH weight
        elif r <= 51:
            weights[r] = 1.0    # Gap — normal weight
        else:
            weights[r] = 2.0    # Near output — high weight

    best_hash_hw = hash_hw(W1, W2)
    best_W2 = list(W2)
    init_hw = best_hash_hw

    for outer in range(max_outer):
        # Step 1: Compute all required schedule values
        required, actual_exp = compute_all_required(W1, W2)

        # Step 2: For each word, compute weighted mismatch if we correct for ONE round
        # Strategy: find word w and correction that most reduces weighted mismatch

        target_rounds = list(range(64))
        jac = compute_schedule_jacobian(W2, target_rounds)

        best_word = -1
        best_correction = 0
        best_delta_wm = 0

        for w in range(16):
            # For each round r, the ideal correction for word w is:
            # Δw = (required[r] - actual[r]) / sens(r,w)
            # But we want to balance ALL rounds.
            # Heuristic: pick correction that helps the HIGHEST-WEIGHT round

            # Find round with max weight × mismatch for this word
            best_r = -1
            best_weighted_mm = 0
            for r in target_rounds:
                sens = jac.get((r, w), 0)
                if sens == 0 or sens & 1 == 0:
                    continue
                mm = sub32(required[r], actual_exp[r])
                wmm = weights[r] * hw(mm)
                if wmm > best_weighted_mm:
                    best_weighted_mm = wmm
                    best_r = r

            if best_r < 0:
                continue

            sens = jac[(best_r, w)]
            mm = sub32(required[best_r], actual_exp[best_r])
            inv = pow(sens, -1, 1 << 32)
            correction = (mm * inv) & MASK32

            # Apply and measure
            W2_trial = list(W2)
            W2_trial[w] = add32(W2[w], correction)

            if W2_trial == W1:
                continue

            trial_hash = hash_hw(W1, W2_trial)
            delta = trial_hash - best_hash_hw

            if delta < best_delta_wm:
                best_delta_wm = delta
                best_word = w
                best_correction = correction
                best_trial_W2 = list(W2_trial)
                best_trial_hash = trial_hash

        if best_word >= 0 and best_delta_wm < 0:
            W2 = best_trial_W2
            if best_trial_hash < best_hash_hw:
                best_hash_hw = best_trial_hash
                best_W2 = list(W2)

            print(f"    Iter {outer}: W[{best_word}] corrected, "
                  f"HW(δH)={best_trial_hash} (Δ={best_delta_wm:+d})")
        else:
            # No single-word correction helps. Try word PAIR.
            improved_pair = False
            for w1 in range(16):
                if improved_pair: break
                for w2 in range(w1+1, 16):
                    # Correct w1 for its best round, w2 for its best round
                    W2_trial = list(W2)

                    for w in [w1, w2]:
                        best_r_w = -1; best_wmm_w = 0
                        for r in target_rounds:
                            sens = jac.get((r, w), 0)
                            if sens == 0 or sens & 1 == 0: continue
                            mm = sub32(required[r], actual_exp[r])
                            wmm = weights[r] * hw(mm)
                            if wmm > best_wmm_w:
                                best_wmm_w = wmm; best_r_w = r
                        if best_r_w >= 0:
                            sens = jac[(best_r_w, w)]
                            mm = sub32(required[best_r_w], actual_exp[best_r_w])
                            inv = pow(sens, -1, 1<<32)
                            corr = (mm * inv) & MASK32
                            W2_trial[w] = add32(W2[w], corr)

                    if W2_trial == W1: continue
                    trial_h = hash_hw(W1, W2_trial)
                    if trial_h < best_hash_hw:
                        W2 = W2_trial
                        best_hash_hw = trial_h
                        best_W2 = list(W2)
                        improved_pair = True
                        print(f"    Iter {outer}: PAIR W[{w1}]+W[{w2}], "
                              f"HW(δH)={trial_h}")
                        break

            if not improved_pair:
                print(f"    Iter {outer}: no improvement")
                break

    return best_W2, best_hash_hw, init_hw


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10

    print("="*70)
    print("SCF: UNIFIED SOLVER v2 — weighted joint optimization")
    print("="*70)

    v2_results = []
    sa_results = []
    import math

    for trial in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

        print(f"\n  Trial {trial}:")
        W2_v2, hw_v2, init_hw = v2_solve(W1, max_outer=15)
        v2_results.append(hw_v2)

        # SA comparison (same eval budget ≈ 16*15 = 240 SHA evals)
        W2_sa = list(W1); W2_sa[0] ^= 0x80000000
        best_sa = hash_hw(W1, W2_sa); cur = list(W2_sa)
        for it in range(3000):
            t = list(cur)
            w = int.from_bytes(os.urandom(1),'big') % 16
            b = int.from_bytes(os.urandom(1),'big') % 32
            t[w] ^= (1<<b)
            if t==W1: continue
            s = hash_hw(W1, t)
            T = max(0.01,1-it/3000)
            if s<best_sa or math.exp(-(s-best_sa)/(T*2))>int.from_bytes(os.urandom(4),'big')/(1<<32):
                cur=t
                if s<best_sa: best_sa=s
        sa_results.append(best_sa)

        print(f"    v2={hw_v2}, SA={best_sa} (init={init_hw})")

    print(f"\n{'='*70}")
    print(f"SUMMARY ({N} trials)")
    print(f"  v2:  avg={sum(v2_results)/N:.1f}, min={min(v2_results)}")
    print(f"  SA:  avg={sum(sa_results)/N:.1f}, min={min(sa_results)}")
    adv = sum(sa_results)/N - sum(v2_results)/N
    print(f"  v2 advantage: {adv:+.1f} (positive = v2 better)")
