#!/usr/bin/env python3
"""
SCF: GLOBAL MISMATCH SOLVER — все раунды одновременно.

Carry Solver ломался потому что корректировал ПО ОДНОМУ раунду.
Новый подход: для КАЖДОГО раунда r вычислить required W2[r].
Затем найти W2[0..15], чей schedule БЛИЖЕ ВСЕГО к required.

Это LEAST SQUARES над Z_{2^32}:
  minimize Σ_r || schedule(W2)[r] - required[r] ||^2

16 переменных (W2[0..15]), 48 условий (W[16..63]).
Schedule ПОЧТИ линеен → линейное приближение + итерация.
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


def compute_required_schedule(W1, W2):
    """For pair (W1, W2), compute what W2's schedule SHOULD be
    to make δe[r+1]=0 at every round r."""
    states1, Wexp1 = sha_all_states(W1)
    states2, Wexp2 = sha_all_states(W2)

    required = list(Wexp2)  # Start with actual, override gap rounds

    for r in range(17, 52):
        s1r = states1[r]; s2r = states2[r]
        # T1_1 = e1[r+1] - d1[r]
        T1_1 = sub32(states1[r+1][4], s1r[3])
        # Required T1_2 = T1_1 + d1 - d2
        T1_2_req = add32(T1_1, sub32(s1r[3], s2r[3]))
        # Required W2[r]
        W2_req = sub32(sub32(sub32(sub32(T1_2_req,
                  s2r[7]), Sig1(s2r[4])), Ch(s2r[4],s2r[5],s2r[6])), K[r])
        required[r] = W2_req

    return required, Wexp2


def global_mismatch_score(W1, W2):
    """Total HW of mismatch between required and actual schedule."""
    required, actual = compute_required_schedule(W1, W2)
    return sum(hw(required[r] ^ actual[r]) for r in range(17, 52))


def global_mismatch_per_round(W1, W2):
    required, actual = compute_required_schedule(W1, W2)
    return [(r, hw(required[r] ^ actual[r])) for r in range(17, 52)]


# ============================================================
# Tool 1: Gradient descent on mismatch (word-by-word)
# ============================================================
def gradient_mismatch(W1, W2_init, max_iter=200):
    """For each word, try +1/-1 perturbation, pick best."""
    W2 = list(W2_init)
    score = global_mismatch_score(W1, W2)
    best_score = score
    best_W2 = list(W2)

    for iteration in range(max_iter):
        improved = False
        for word in range(16):
            # Try +1
            W2_up = list(W2); W2_up[word] = add32(W2[word], 1)
            s_up = global_mismatch_score(W1, W2_up)

            # Try -1
            W2_dn = list(W2); W2_dn[word] = sub32(W2[word], 1)
            s_dn = global_mismatch_score(W1, W2_dn)

            if s_up < score:
                W2 = W2_up; score = s_up; improved = True
            if s_dn < score:
                W2 = W2_dn; score = s_dn; improved = True

        if score < best_score:
            best_score = score
            best_W2 = list(W2)

        if not improved:
            # Random kick
            w = int.from_bytes(os.urandom(1),'big') % 16
            W2[w] ^= (1 << (int.from_bytes(os.urandom(1),'big') % 32))
            score = global_mismatch_score(W1, W2)

        if iteration < 5 or iteration % 50 == 0:
            print(f"    Iter {iteration}: mismatch={best_score}")

    return best_W2, best_score


# ============================================================
# Tool 2: SA on global mismatch
# ============================================================
def sa_global_mismatch(W1, W2_init, max_iter=5000):
    W2 = list(W2_init)
    score = global_mismatch_score(W1, W2)
    best_score = score; best_W2 = list(W2)

    for it in range(max_iter):
        trial = list(W2)
        w = int.from_bytes(os.urandom(1),'big') % 16
        b = int.from_bytes(os.urandom(1),'big') % 32
        trial[w] ^= (1 << b)

        s = global_mismatch_score(W1, trial)
        T = max(0.01, 1.0 - it/max_iter)
        if s < score or math.exp(-(s-score)/(T*5)) > int.from_bytes(os.urandom(4),'big')/(1<<32):
            W2 = trial; score = s
        if score < best_score:
            best_score = score; best_W2 = list(W2)

    return best_W2, best_score


# ============================================================
# Tool 3: Schedule-aware Newton — use schedule Jacobian
# ============================================================
def schedule_jacobian_correct(W1, W2_init, target_r, max_iter=10):
    """Use schedule Jacobian to correct mismatch at specific round.
    Key: compute ∂W_expanded[target_r]/∂W[word] for each word."""
    W2 = list(W2_init)

    for iteration in range(max_iter):
        required, actual = compute_required_schedule(W1, W2)
        mismatch = sub32(required[target_r], actual[target_r])
        if mismatch == 0:
            return W2, True

        # Schedule Jacobian at target_r
        best_word = -1; best_residual = hw(mismatch)

        for word in range(16):
            W2_test = list(W2)
            W2_test[word] = add32(W2[word], 1)
            Wexp_test = expand_real(W2_test)
            Wexp_base = expand_real(W2)

            sens = sub32(Wexp_test[target_r], Wexp_base[target_r])
            if sens == 0 or sens & 1 == 0:
                continue

            inv = pow(sens, -1, 1 << 32)
            correction = (mismatch * inv) & MASK32

            W2_corr = list(W2)
            W2_corr[word] = add32(W2[word], correction)
            Wexp_corr = expand_real(W2_corr)
            new_mismatch = sub32(required[target_r], Wexp_corr[target_r])

            if hw(new_mismatch) < best_residual:
                best_residual = hw(new_mismatch)
                best_word = word
                best_W2 = list(W2_corr)

        if best_word >= 0:
            W2 = best_W2
        else:
            break

    return W2, False


# ============================================================
# MAIN: Combined pipeline
# ============================================================
if __name__ == '__main__':
    print("="*70)
    print("SCF: GLOBAL MISMATCH SOLVER")
    print("="*70)

    W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    W2 = list(W1); W2[0] ^= 0x80000000

    # Baseline
    baseline = global_mismatch_score(W1, W2)
    print(f"\n  Baseline mismatch: {baseline} bits ({baseline/35:.1f}/round)")

    # Tool 1: Gradient
    print(f"\n  Tool 1: Gradient descent on global mismatch")
    W2_grad, score_grad = gradient_mismatch(W1, list(W2), max_iter=100)
    H1 = sha_compress(W1); H2g = sha_compress(W2_grad)
    hw_grad = sum(hw(H1[i]^H2g[i]) for i in range(8))
    print(f"  Result: mismatch={score_grad} ({score_grad/35:.1f}/r), hash HW={hw_grad}")

    # Tool 2: SA
    print(f"\n  Tool 2: SA on global mismatch")
    W2_sa, score_sa = sa_global_mismatch(W1, list(W2), max_iter=10000)
    H2s = sha_compress(W2_sa)
    hw_sa = sum(hw(H1[i]^H2s[i]) for i in range(8))
    print(f"  Result: mismatch={score_sa} ({score_sa/35:.1f}/r), hash HW={hw_sa}")

    # Tool 3: Sequential schedule Newton
    print(f"\n  Tool 3: Sequential schedule Newton")
    W2_sn = list(W2)
    n_solved = 0
    for tr in range(17, 52, 4):
        W2_sn, success = schedule_jacobian_correct(W1, W2_sn, tr, max_iter=5)
        if success: n_solved += 1

    score_sn = global_mismatch_score(W1, W2_sn)
    H2sn = sha_compress(W2_sn)
    hw_sn = sum(hw(H1[i]^H2sn[i]) for i in range(8))
    print(f"  Solved {n_solved}/9 target rounds")
    print(f"  Result: mismatch={score_sn} ({score_sn/35:.1f}/r), hash HW={hw_sn}")

    # Tool 4: COMBINED — Newton per round then SA globally
    print(f"\n  Tool 4: COMBINED (Newton + SA)")
    W2_comb = list(W2_sn)  # Start from Newton result
    W2_comb, score_comb = sa_global_mismatch(W1, W2_comb, max_iter=10000)
    H2c = sha_compress(W2_comb)
    hw_comb = sum(hw(H1[i]^H2c[i]) for i in range(8))
    print(f"  Result: mismatch={score_comb} ({score_comb/35:.1f}/r), hash HW={hw_comb}")

    # Per-round profile of best result
    best_W2 = min([(score_grad, W2_grad), (score_sa, W2_sa),
                    (score_sn, W2_sn), (score_comb, W2_comb)])[1]
    best_total = global_mismatch_score(W1, best_W2)

    print(f"\n  Best overall mismatch: {best_total} bits")
    profile = global_mismatch_per_round(W1, best_W2)
    zeros = sum(1 for r, h in profile if h == 0)
    low = sum(1 for r, h in profile if h <= 5)
    print(f"  Zero-mismatch rounds: {zeros}/35")
    print(f"  Low-mismatch (≤5): {low}/35")

    if zeros > 0:
        print(f"  ★ Rounds with exact schedule match:")
        for r, h in profile:
            if h == 0:
                print(f"    r={r}: mismatch=0 ★★★")

    # Random baseline comparison
    rand_mismatches = []
    for _ in range(100):
        Wr = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        rand_mismatches.append(global_mismatch_score(W1, Wr))
    avg_rand = sum(rand_mismatches)/len(rand_mismatches)
    min_rand = min(rand_mismatches)

    print(f"\n  Comparison:")
    print(f"    Random: avg={avg_rand:.0f}, min={min_rand}")
    print(f"    Our best: {best_total}")
    print(f"    Advantage: {min_rand - best_total:+d} bits vs random min")
