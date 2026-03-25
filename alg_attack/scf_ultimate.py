#!/usr/bin/env python3
"""
SCF: ULTIMATE SOLVER — все открытия в одном инструменте.

Впитывает:
1. Adaptive δ screening (512 кандидатов → top-K)
2. v4 structure solver (exact inversion + δe/δa targeting)
3. Wang chain (sens=1 для ранних раундов)
4. Schedule kernel awareness (128 dims)
5. Multi-start (top-3 δ параллельно)
6. SA hybrid polish
7. Iterative recomputation (fresh states)

Архитектура:
  OUTER: screen δ candidates → top-K
  MID:   for each δ → Wang chain (if applicable) + structure solve
  INNER: SA polish from structured result
  RETURN: best across all starts
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

def sha_full(W16):
    W=expand_real(W16); s=list(IV); states=[list(s)]
    for r in range(64):
        a,b,c,d,e,f,g,h=s
        T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add32(Sig0(a),Maj(a,b,c))
        s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
        states.append(list(s))
    return [add32(IV[i],s[i]) for i in range(8)], states, W

def sha_hash(W16):
    return sha_full(W16)[0]

def hash_diff(H1, W2):
    H2 = sha_hash(W2)
    return sum(hw(H1[i]^H2[i]) for i in range(8))


# ============================================================
# Component 1: DELTA SCREENER — find top-K starting deltas
# ============================================================
def screen_deltas(W1, H1, n_candidates=512, top_k=5):
    """Screen single-bit deltas, return top-K by raw hash diff."""
    candidates = []
    for word in range(16):
        for bit in range(32):
            W2 = list(W1); W2[word] ^= (1 << bit)
            d = hash_diff(H1, W2)
            candidates.append((d, word, bit))
    candidates.sort()
    return candidates[:top_k]


# ============================================================
# Component 2: STRUCTURE SOLVER — exact inversion per round
# ============================================================
def structure_step(W1, W2, H1):
    """One structure-guided correction step.
    Returns improved W2 or None if no improvement found."""
    H1_states = sha_full(W1)[1]
    H2_full, st2, Wexp2 = sha_full(W2)
    current_hw = sum(hw(H1[i] ^ H2_full[i]) for i in range(8))

    best_hw = current_hw
    best_W2 = None

    # Try corrections for rounds 50-63 (closest to hash)
    for w in range(16):
        W2_sens = list(W2); W2_sens[w] = add32(W2[w], 1)
        Wexp_sens = expand_real(W2_sens)

        for r in range(63, 49, -1):
            sens = sub32(Wexp_sens[r], Wexp2[r])
            if sens == 0 or sens & 1 == 0:
                continue

            inv = pow(sens, -1, 1 << 32)
            s1r = H1_states[r]; s2r = st2[r]; s1r1 = H1_states[r+1]

            # Try δe → 0
            T1_1 = sub32(s1r1[4], s1r[3])
            T1_2_need = add32(T1_1, sub32(s1r[3], s2r[3]))
            W_need_e = sub32(sub32(sub32(sub32(T1_2_need,
                       s2r[7]),Sig1(s2r[4])),Ch(s2r[4],s2r[5],s2r[6])),K[r])
            mm_e = sub32(W_need_e, Wexp2[r])
            if mm_e != 0:
                corr = (mm_e * inv) & MASK32
                W2_c = list(W2); W2_c[w] = add32(W2[w], corr)
                if W2_c != W1:
                    trial = hash_diff(H1, W2_c)
                    if trial < best_hw:
                        best_hw = trial; best_W2 = list(W2_c)

            # Try δa → 0
            T2_2 = add32(Sig0(s2r[0]), Maj(s2r[0],s2r[1],s2r[2]))
            T1_2_need_a = sub32(s1r1[0], T2_2)
            W_need_a = sub32(sub32(sub32(sub32(T1_2_need_a,
                       s2r[7]),Sig1(s2r[4])),Ch(s2r[4],s2r[5],s2r[6])),K[r])
            mm_a = sub32(W_need_a, Wexp2[r])
            if mm_a != 0:
                corr = (mm_a * inv) & MASK32
                W2_c = list(W2); W2_c[w] = add32(W2[w], corr)
                if W2_c != W1:
                    trial = hash_diff(H1, W2_c)
                    if trial < best_hw:
                        best_hw = trial; best_W2 = list(W2_c)

            break  # One round per word

    return best_W2, best_hw


# ============================================================
# Component 3: WANG CHAIN — adaptive δW for early rounds
# ============================================================
def wang_chain_adaptive(W1, W2_init, protected_word=-1):
    """Try to set δW[r] to cancel δe at early rounds.
    Protected_word: don't modify this word (keeps δ nontrivial)."""
    W2 = list(W2_init)
    for r in range(1, 16):
        if r == protected_word:
            continue
        st1 = sha_full(W1)[1]
        st2 = sha_full(W2)[1]
        de = sub32(st2[r+1][4], st1[r+1][4])
        if de != 0 and hw(de) < 20:
            W2_try = list(W2)
            W2_try[r] = sub32(W2[r], de)
            if W2_try != W1:  # Don't allow trivial
                W2 = W2_try
    return W2


# ============================================================
# Component 4: SA POLISH — local refinement
# ============================================================
def sa_polish(W1, W2_init, H1, budget, fixed_word=-1, fixed_bit=-1):
    W2 = list(W2_init)
    best = hash_diff(H1, W2)
    best_W2 = list(W2)
    cur = list(W2); cur_s = best

    for it in range(budget):
        t = list(cur)
        w = int.from_bytes(os.urandom(1),'big') % 16
        b = int.from_bytes(os.urandom(1),'big') % 32
        if w == fixed_word and b == fixed_bit: continue
        t[w] ^= (1 << b)
        if t == W1: continue
        s = hash_diff(H1, t)
        T = max(0.01, 1 - it/budget)
        if s < cur_s or math.exp(-(s-cur_s)/(T*2)) > int.from_bytes(os.urandom(4),'big')/(1<<32):
            cur = t; cur_s = s
            if s < best:
                best = s; best_W2 = list(t)
    return best_W2, best


# ============================================================
# ULTIMATE SOLVER — all components combined
# ============================================================
def ultimate_solve(W1, total_budget=3000):
    H1 = sha_hash(W1)
    evals = 0

    # Phase 1: Screen deltas (cost: 512 evals)
    top_deltas = screen_deltas(W1, H1, top_k=5)
    evals += 512

    budget_per_delta = (total_budget - 512) // len(top_deltas)

    best_overall_hw = 256
    best_overall_W2 = None

    for raw_hw, d_word, d_bit in top_deltas:
        W2 = list(W1)
        W2[d_word] ^= (1 << d_bit)

        # Phase 2: Wang chain attempt (protect the delta word)
        W2 = wang_chain_adaptive(W1, W2, protected_word=d_word)
        evals += 30

        # Phase 3: Structure solver (iterative)
        struct_budget = budget_per_delta // 3
        for _ in range(min(struct_budget // 32, 15)):
            improved_W2, improved_hw = structure_step(W1, W2, H1)
            evals += 32  # ~16 words × 2 targets
            if improved_W2 and improved_W2 != W1 and improved_hw < hash_diff(H1, W2):
                W2 = improved_W2
            else:
                break

        # Phase 4: SA polish
        sa_budget = budget_per_delta - (evals - 512) // len(top_deltas)
        sa_budget = max(sa_budget, 100)
        W2, final_hw = sa_polish(W1, W2, H1, sa_budget, d_word, d_bit)
        evals += sa_budget

        if final_hw < best_overall_hw:
            best_overall_hw = final_hw
            best_overall_W2 = list(W2)

    return best_overall_W2, best_overall_hw, evals


# ============================================================
# BENCHMARK
# ============================================================
if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    budget = int(sys.argv[2]) if len(sys.argv) > 2 else 3000

    print("="*70)
    print(f"SCF: ULTIMATE SOLVER (budget={budget})")
    print("="*70)

    ult_results = []
    sa_results = []

    for trial in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        H1 = sha_hash(W1)

        # Ultimate solver
        W2_ult, hw_ult, ev_ult = ultimate_solve(W1, total_budget=budget)
        ult_results.append(hw_ult)

        # Pure SA with same budget
        W2_sa = list(W1); W2_sa[0] ^= 0x80000000
        best_sa = hash_diff(H1, W2_sa); cur = list(W2_sa)
        for it in range(budget):
            t = list(cur)
            w = int.from_bytes(os.urandom(1),'big') % 16
            b = int.from_bytes(os.urandom(1),'big') % 32
            t[w] ^= (1<<b)
            if t == W1: continue
            s = hash_diff(H1, t)
            T = max(0.01, 1-it/budget)
            if s < best_sa or math.exp(-(s-best_sa)/(T*2)) > int.from_bytes(os.urandom(4),'big')/(1<<32):
                cur = t
                if s < best_sa: best_sa = s
        sa_results.append(best_sa)

        winner = "ULT ★" if hw_ult < best_sa else ("SA" if best_sa < hw_ult else "tie")
        print(f"  Trial {trial:2d}: Ultimate={hw_ult:3d}  SA={best_sa:3d}  → {winner}")

    print(f"\n{'='*70}")
    avg_u = sum(ult_results)/N; avg_s = sum(sa_results)/N
    ult_wins = sum(1 for i in range(N) if ult_results[i] < sa_results[i])
    sa_wins = sum(1 for i in range(N) if sa_results[i] < ult_results[i])

    print(f"SUMMARY ({N} trials, budget={budget})")
    print(f"  Ultimate: avg={avg_u:.1f}  min={min(ult_results)}")
    print(f"  SA:       avg={avg_s:.1f}  min={min(sa_results)}")
    print(f"  Ultimate wins: {ult_wins}  SA wins: {sa_wins}  ties: {N-ult_wins-sa_wins}")
    print(f"  Advantage: {avg_s-avg_u:+.1f} (positive = Ultimate better)")
