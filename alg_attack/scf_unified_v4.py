#!/usr/bin/env python3
"""
SCF: UNIFIED SOLVER v4 — v3 improvements + δa targeting + hybrid SA.

v3 = SA при 540 evals. v4 цель: ПРЕВЗОЙТИ SA.

Улучшения:
1. Target BOTH δe AND δa (два нелинейных узла → 2× info)
2. Score = hash_hw напрямую (не mismatch proxy)
3. After structure phase → short SA polish (hybrid)
4. Multiple correction candidates per iteration → pick BEST by hash
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
    """Returns (hash, all_states, expanded_W)."""
    W=expand_real(W16); s=list(IV); states=[list(s)]
    for r in range(64):
        a,b,c,d,e,f,g,h=s
        T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add32(Sig0(a),Maj(a,b,c))
        s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
        states.append(list(s))
    H = [add32(IV[i],s[i]) for i in range(8)]
    return H, states, W

def hash_hw_fast(H1, W2):
    H2,_,_ = sha_full(W2)
    return sum(hw(H1[i]^H2[i]) for i in range(8))


def invert_for_de_zero(s1_at_r, s2_at_r, s1_after_r, r):
    """Required W2[r] for δe[r+1]=0."""
    T1_1 = sub32(s1_after_r[4], s1_at_r[3])
    T1_2_needed = add32(T1_1, sub32(s1_at_r[3], s2_at_r[3]))
    return sub32(sub32(sub32(sub32(T1_2_needed,
           s2_at_r[7]),Sig1(s2_at_r[4])),Ch(s2_at_r[4],s2_at_r[5],s2_at_r[6])),K[r])

def invert_for_da_zero(s1_at_r, s2_at_r, s1_after_r, r):
    """Required W2[r] for δa[r+1]=0. a_new = T1+T2."""
    a1_target = s1_after_r[0]
    T2_2 = add32(Sig0(s2_at_r[0]), Maj(s2_at_r[0],s2_at_r[1],s2_at_r[2]))
    T1_2_needed = sub32(a1_target, T2_2)
    return sub32(sub32(sub32(sub32(T1_2_needed,
           s2_at_r[7]),Sig1(s2_at_r[4])),Ch(s2_at_r[4],s2_at_r[5],s2_at_r[6])),K[r])


def v4_solve(W1, delta_w0=0x80000000, struct_iters=25, sa_iters=300):
    W2 = list(W1); W2[0] ^= delta_w0

    H1, st1_all, Wexp1 = sha_full(W1)
    best_hw = hash_hw_fast(H1, W2)
    best_W2 = list(W2)
    evals = 1

    # === PHASE 1: Structure-guided corrections ===
    for iteration in range(struct_iters):
        H2, st2_all, Wexp2 = sha_full(W2)
        current_hw = sum(hw(H1[i]^H2[i]) for i in range(8))
        evals += 1

        if current_hw == 0:
            print(f"    [{iteration}] ★★★ COLLISION!")
            return W2, 0, evals

        # Generate candidates: for each (word, round, target=e|a)
        candidates = []

        for w in range(16):
            W2_sens = list(W2); W2_sens[w] = add32(W2[w], 1)
            Wexp_sens = expand_real(W2_sens)

            for r in range(55, max(15, 55-iteration*2), -1):
                sens = sub32(Wexp_sens[r], Wexp2[r])
                if sens == 0 or sens & 1 == 0:
                    continue

                inv = pow(sens, -1, 1<<32)

                # Target δe[r+1]=0
                W2_needed_e = invert_for_de_zero(
                    st1_all[r], st2_all[r], st1_all[r+1], r)
                mm_e = sub32(W2_needed_e, Wexp2[r])
                if mm_e != 0:
                    corr = (mm_e * inv) & MASK32
                    W2_c = list(W2); W2_c[w] = add32(W2[w], corr)
                    if W2_c != W1:
                        candidates.append((w, r, 'e', W2_c))

                # Target δa[r+1]=0
                W2_needed_a = invert_for_da_zero(
                    st1_all[r], st2_all[r], st1_all[r+1], r)
                mm_a = sub32(W2_needed_a, Wexp2[r])
                if mm_a != 0:
                    corr = (mm_a * inv) & MASK32
                    W2_c = list(W2); W2_c[w] = add32(W2[w], corr)
                    if W2_c != W1:
                        candidates.append((w, r, 'a', W2_c))

                break  # One round per word

        # Evaluate top candidates (limit evals)
        best_cand_hw = current_hw
        best_cand = None

        for w, r, target, W2_c in candidates[:32]:
            ch = hash_hw_fast(H1, W2_c); evals += 1
            if ch < best_cand_hw:
                best_cand_hw = ch
                best_cand = (w, r, target, W2_c)

        if best_cand and best_cand_hw < current_hw:
            w, r, target, W2_c = best_cand
            W2 = W2_c
            if best_cand_hw < best_hw:
                best_hw = best_cand_hw
                best_W2 = list(W2)
            if iteration < 5 or best_cand_hw < 95:
                print(f"    [{iteration}] W[{w}] r={r} δ{target}→0, HW={best_cand_hw} (Δ={best_cand_hw-current_hw:+d})")
        else:
            # Escape: random word perturbation
            w = iteration % 16
            W2[w] = add32(W2[w], iteration*7 + 3)
            if W2 == W1: W2[w] = add32(W2[w], 1)

    struct_best = best_hw
    struct_evals = evals

    # === PHASE 2: SA polish from structured result ===
    W2_sa = list(best_W2)
    sa_best = best_hw

    for it in range(sa_iters):
        t = list(W2_sa)
        w = int.from_bytes(os.urandom(1),'big') % 16
        b = int.from_bytes(os.urandom(1),'big') % 32
        t[w] ^= (1<<b)
        if t == W1: continue
        s = hash_hw_fast(H1, t); evals += 1
        T = max(0.01, 1-it/sa_iters)
        if s < sa_best or math.exp(-(s-sa_best)/(T*2)) > int.from_bytes(os.urandom(4),'big')/(1<<32):
            W2_sa = t
            if s < sa_best:
                sa_best = s
                if s < best_hw:
                    best_hw = s
                    best_W2 = list(W2_sa)

    if best_hw < struct_best:
        print(f"    SA polish: {struct_best} → {best_hw}")

    return best_W2, best_hw, evals


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10

    print("="*70)
    print("SCF: UNIFIED SOLVER v4 — δe+δa targeting + SA hybrid")
    print("="*70)

    v4_results = []; sa_results = []; v4_evals_list = []

    for trial in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        init = hash_hw_fast(sha_full(W1)[0], [W1[0]^0x80000000]+W1[1:])

        print(f"\n  Trial {trial} (init={init}):")
        W2, hw_v4, ev4 = v4_solve(W1, struct_iters=20, sa_iters=500)
        v4_results.append(hw_v4); v4_evals_list.append(ev4)

        # Pure SA with same total evals
        W2_sa = list(W1); W2_sa[0] ^= 0x80000000
        H1 = sha_full(W1)[0]
        best_sa = hash_hw_fast(H1, W2_sa); cur = list(W2_sa)
        for it in range(ev4):
            t=list(cur); w=int.from_bytes(os.urandom(1),'big')%16
            b=int.from_bytes(os.urandom(1),'big')%32
            t[w]^=(1<<b)
            if t==W1: continue
            s=hash_hw_fast(H1,t)
            T=max(0.01,1-it/ev4)
            if s<best_sa or math.exp(-(s-best_sa)/(T*2))>int.from_bytes(os.urandom(4),'big')/(1<<32):
                cur=t
                if s<best_sa: best_sa=s
        sa_results.append(best_sa)

        winner = "v4 ★" if hw_v4 < best_sa else ("SA" if best_sa < hw_v4 else "tie")
        print(f"    v4={hw_v4}, SA={best_sa} ({ev4} evals) → {winner}")

    print(f"\n{'='*70}")
    print(f"FINAL SUMMARY ({N} trials)")
    avg_v4 = sum(v4_results)/N; avg_sa = sum(sa_results)/N
    avg_ev = sum(v4_evals_list)/N
    v4_wins = sum(1 for i in range(N) if v4_results[i] < sa_results[i])
    sa_wins = sum(1 for i in range(N) if sa_results[i] < v4_results[i])
    ties = N - v4_wins - sa_wins

    print(f"  v4:  avg={avg_v4:.1f}, min={min(v4_results)}")
    print(f"  SA:  avg={avg_sa:.1f}, min={min(sa_results)}")
    print(f"  Avg evals: {avg_ev:.0f}")
    print(f"  v4 wins: {v4_wins}, SA wins: {sa_wins}, ties: {ties}")
    print(f"  v4 advantage: {avg_sa-avg_v4:+.1f}")
