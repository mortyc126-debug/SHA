#!/usr/bin/env python3
"""
SCF: UNIFIED SOLVER v3 — structure-guided с iterative recomputation.

v2 проблемы: коррекция W[word] меняет states → required меняется.
v3: после КАЖДОЙ коррекции ПЕРЕСЧИТЫВАЕМ required для всех раундов.

Плюс: больший eval budget, multi-word corrections, структурные веса
основанные на РАССТОЯНИИ до hash output (ближе → важнее).
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

def sha_compress(W16):
    W=expand_real(W16); s=list(IV)
    for r in range(64):
        a,b,c,d,e,f,g,h=s
        T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add32(Sig0(a),Maj(a,b,c))
        s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
    return [add32(IV[i],s[i]) for i in range(8)]

def hash_hw(W1, W2):
    H1=sha_compress(W1); H2=sha_compress(W2)
    return sum(hw(H1[i]^H2[i]) for i in range(8))


def v3_solve(W1, delta_w0=0x80000000, max_iter=50):
    """
    v3: iterative structure-guided correction.

    Each iteration:
    1. Run FULL SHA for both messages → get states + expanded schedule
    2. For EACH round r, compute required W2[r] (exact inversion)
    3. For EACH word w, compute: if I correct w for round r_best,
       what happens to the HASH? (not mismatch — HASH directly)
    4. Apply best (word, correction) pair
    5. Repeat with FRESH states
    """
    W2 = list(W1)
    W2[0] ^= delta_w0

    best_hw = hash_hw(W1, W2)
    best_W2 = list(W2)
    evals = 0

    for iteration in range(max_iter):
        # Fresh computation every iteration
        H1 = sha_compress(W1)
        current_hw = sum(hw(H1[i] ^ sha_compress(W2)[i]) for i in range(8))
        evals += 2

        if current_hw == 0:
            print(f"    Iter {iteration}: ★★★ COLLISION!")
            break

        # For each word: try correcting for the round that helps hash MOST
        best_word = -1
        best_trial_hw = current_hw
        best_trial_W2 = None

        for w in range(16):
            # Compute schedule sensitivity for this word
            W2_test = list(W2)
            W2_test[w] = add32(W2[w], 1)
            Wexp_base = expand_real(W2)
            Wexp_test = expand_real(W2_test)

            # Find rounds where this word has invertible sensitivity
            for r in range(63, -1, -1):  # Start from end (closest to hash)
                sens = sub32(Wexp_test[r], Wexp_base[r])
                if sens == 0 or sens & 1 == 0:
                    continue

                # Compute required W2[r] for δe[r+1]=0
                # Need both messages' states at round r
                W_exp1 = expand_real(W1); s1 = list(IV)
                for rr in range(r+1):
                    a,b,c,d,e,f,g,h=s1
                    T1=add32(h,Sig1(e),Ch(e,f,g),K[rr],W_exp1[rr])
                    T2=add32(Sig0(a),Maj(a,b,c))
                    s1=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]

                s2 = list(IV)
                for rr in range(r+1):
                    a,b,c,d,e,f,g,h=s2
                    T1=add32(h,Sig1(e),Ch(e,f,g),K[rr],Wexp_base[rr])
                    T2=add32(Sig0(a),Maj(a,b,c))
                    s2=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]

                # Mismatch for this round
                mm = sub32(Wexp_base[r], Wexp_base[r])  # Will be computed properly
                # Actually use exact inversion
                e1_next = s1[4]  # e of state after round r
                d2_at_r_prev = s2[3] if r > 0 else IV[3]
                # Recompute states at r (before the round)
                s1_at_r = list(IV); s2_at_r = list(IV)
                for rr in range(r):
                    a,b,c,d,e,f,g,h=s1_at_r
                    T1=add32(h,Sig1(e),Ch(e,f,g),K[rr],W_exp1[rr])
                    T2=add32(Sig0(a),Maj(a,b,c))
                    s1_at_r=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
                    a,b,c,d,e,f,g,h=s2_at_r
                    T1=add32(h,Sig1(e),Ch(e,f,g),K[rr],Wexp_base[rr])
                    T2=add32(Sig0(a),Maj(a,b,c))
                    s2_at_r=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]

                T1_1_needed = sub32(s1[4], s1_at_r[3])
                T1_2_needed = add32(T1_1_needed, sub32(s1_at_r[3], s2_at_r[3]))
                W2_needed = sub32(sub32(sub32(sub32(T1_2_needed,
                            s2_at_r[7]),Sig1(s2_at_r[4])),
                            Ch(s2_at_r[4],s2_at_r[5],s2_at_r[6])),K[r])

                schedule_mismatch = sub32(W2_needed, Wexp_base[r])
                if schedule_mismatch == 0:
                    continue

                # Correction
                inv = pow(sens, -1, 1 << 32)
                correction = (schedule_mismatch * inv) & MASK32

                W2_corr = list(W2)
                W2_corr[w] = add32(W2[w], correction)

                if W2_corr == W1:
                    continue

                trial_hw = hash_hw(W1, W2_corr)
                evals += 1

                if trial_hw < best_trial_hw:
                    best_trial_hw = trial_hw
                    best_word = w
                    best_trial_W2 = list(W2_corr)
                    best_round = r

                break  # Only try best round per word

        if best_word >= 0 and best_trial_hw < current_hw:
            W2 = best_trial_W2
            if best_trial_hw < best_hw:
                best_hw = best_trial_hw
                best_W2 = list(W2)

            if iteration < 10 or iteration % 10 == 0 or best_trial_hw < 95:
                print(f"    Iter {iteration}: W[{best_word}] for r={best_round}, "
                      f"HW={best_trial_hw} (Δ={best_trial_hw-current_hw:+d})")
        else:
            # No single-word helps → try random word perturbation
            w = iteration % 16
            W2[w] = add32(W2[w], iteration + 1)
            if W2 == W1:
                W2[w] = add32(W2[w], 1)

    return best_W2, best_hw, evals


if __name__ == '__main__':
    import math
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10

    print("="*70)
    print("SCF: UNIFIED SOLVER v3")
    print("  Structure-guided + iterative recomputation")
    print("="*70)

    v3_results = []
    sa_results = []

    for trial in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        init = hash_hw(W1, [W1[0]^0x80000000]+W1[1:])

        print(f"\n  Trial {trial} (init={init}):")
        W2, hw_v3, evals_v3 = v3_solve(W1, max_iter=30)
        v3_results.append(hw_v3)

        # SA with same evals
        W2_sa = list(W1); W2_sa[0] ^= 0x80000000
        best_sa = hash_hw(W1, W2_sa); cur = list(W2_sa)
        for it in range(evals_v3):
            t=list(cur); w=int.from_bytes(os.urandom(1),'big')%16
            b=int.from_bytes(os.urandom(1),'big')%32
            t[w]^=(1<<b)
            if t==W1: continue
            s=hash_hw(W1,t)
            T=max(0.01,1-it/max(evals_v3,1))
            if s<best_sa or math.exp(-(s-best_sa)/(T*2))>int.from_bytes(os.urandom(4),'big')/(1<<32):
                cur=t
                if s<best_sa: best_sa=s
        sa_results.append(best_sa)

        print(f"    v3={hw_v3} ({evals_v3} evals), SA={best_sa} ({evals_v3} evals)")

    print(f"\n{'='*70}")
    print(f"SUMMARY ({N} trials, SAME eval budget)")
    avg_v3 = sum(v3_results)/N
    avg_sa = sum(sa_results)/N
    print(f"  v3: avg={avg_v3:.1f}, min={min(v3_results)}")
    print(f"  SA: avg={avg_sa:.1f}, min={min(sa_results)}")
    print(f"  v3 advantage: {avg_sa-avg_v3:+.1f}")
