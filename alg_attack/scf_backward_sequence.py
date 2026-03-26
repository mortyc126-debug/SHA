#!/usr/bin/env python3
"""
SCF: BACKWARD SEQUENCE CONSTRUCTOR — строим a-sequence назад от цели.

Ключ: a[r+1] = T1[r] + T2[r], где T1 и T2 вычислимы из a[r..r-7].
Обращение: из a[r+1] и a[r..r-6] вычислить W[r].

W[r] = a[r+1] - T2[r] - h[r] - Σ1(e[r]) - Ch(e[r],f[r],g[r]) - K[r]

Всё выражается через a-sequence (D-coupling theorem).

Алгоритм:
1. Зафиксировать цель: a2[53..64] = a1[53..64]
2. Для каждого r от 63 к 53: вычислить НУЖНЫЙ W2[r]
3. Нужный W2[r] должен совпасть с schedule(W2[0..15])
4. Это 12 уравнений на 16 неизвестных → rank=12 → 4 свободных слова
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

def sha_full_a(W16):
    """Return a-sequence + all states + expanded W."""
    W=expand_real(W16); s=list(IV); states=[list(s)]
    a_seq=[s[0]]
    for r in range(64):
        a,b,c,d,e,f,g,h=s
        T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add32(Sig0(a),Maj(a,b,c))
        s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
        states.append(list(s))
        a_seq.append(s[0])
    H=[add32(IV[i],s[i]) for i in range(8)]
    return a_seq, states, W, H


def required_W_from_states(state_before, target_a_next, r):
    """Given state at round r, compute W[r] needed to get target a[r+1].
    a[r+1] = T1 + T2.  T2 = Σ0(a)+Maj(a,b,c). T1 = target - T2.
    W = T1 - h - Σ1(e) - Ch(e,f,g) - K[r]."""
    a,b,c,d,e,f,g,h = state_before
    T2 = add32(Sig0(a), Maj(a,b,c))
    T1 = sub32(target_a_next, T2)
    return sub32(sub32(sub32(sub32(T1, h), Sig1(e)), Ch(e,f,g)), K[r])


# ============================================================
# EXP 1: Backward inversion — required W[53..63] for collision
# ============================================================
def exp1_backward_required(N):
    print("="*70)
    print("EXP 1: BACKWARD SEQUENCE INVERSION")
    print("  For collision: a2[53..64] = a1[53..64]")
    print("  Compute required W2[52..63] from backward inversion")
    print("  Compare with actual schedule W2[52..63]")
    print("="*70)

    for trial in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = list(W1); W2[0] ^= 0x80000000

        a1, st1, Wexp1, H1 = sha_full_a(W1)
        a2, st2, Wexp2, H2 = sha_full_a(W2)

        # For each round 52..63: compute required W2[r] so that
        # a2[r+1] = a1[r+1] GIVEN the actual state2 at round r.
        print(f"\n  Trial {trial}:")
        print(f"  {'r':>3} | {'W2_required':>12} {'W2_actual':>12} {'mismatch_HW':>12}")
        print("  " + "-"*45)

        total_mm = 0
        for r in range(52, 64):
            W_req = required_W_from_states(st2[r], a1[r+1], r)
            W_act = Wexp2[r]
            mm = hw(W_req ^ W_act)
            total_mm += mm
            print(f"  {r:3d} | 0x{W_req:08x} 0x{W_act:08x} {mm:12d}")

        print(f"  Total mismatch: {total_mm} bits ({total_mm/12:.1f}/round)")

        # Now the KEY question: this uses state2 at round r.
        # But state2 ALSO depends on W2. Changing W2[52..63]
        # would change state2[53..64] → different required W.
        # This is the COUPLED problem.
        # BUT: if we fix state2[0..52], then state2[52] is fixed,
        # and we only need to solve for W2[52..63] given state2[52].


# ============================================================
# EXP 2: How many free words remain after fixing a[53..64]?
# ============================================================
def exp2_free_words():
    print("\n" + "="*70)
    print("EXP 2: FREE WORDS — schedule vs backward constraints")
    print("  Schedule: W[r] = σ1(W[r-2])+W[r-7]+σ0(W[r-15])+W[r-16]")
    print("  Backward: W[r] = required(state, target_a)")
    print("  Intersection: how many W[0..15] satisfy both?")
    print("="*70)

    # Schedule for r=52..63 depends on specific W[0..15] words:
    # W[52] = f(W[50],W[45],W[37],W[36])
    #        = f(f(W[48],...), W[45], f(W[35],...), f(W[34],...))
    # Everything traces back to W[0..15].

    # Count: how many of W[0..15] affect W[52..63]?
    W_base = [0]*16
    affected = set()
    for word in range(16):
        for r in range(52, 64):
            W_test = [0]*16; W_test[word] = 1
            Wexp_test = expand_real(W_test)
            if Wexp_test[r] != 0:
                affected.add(word)
                break

    print(f"  W[0..15] words affecting W[52..63]: {sorted(affected)}")
    print(f"  Count: {len(affected)}/16")
    print(f"  Free words (not affecting tail): {sorted(set(range(16))-affected)}")

    # More precisely: schedule Jacobian for rounds 52..63
    W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    Wexp1 = expand_real(W1)

    # J[r][word] = sensitivity of W_expanded[r] to W[word]
    J = {}
    for r in range(52, 64):
        for word in range(16):
            W_test = list(W1)
            W_test[word] = add32(W1[word], 1)
            Wexp_test = expand_real(W_test)
            sens = sub32(Wexp_test[r], Wexp1[r])
            J[(r, word)] = sens

    # How many words have INVERTIBLE sensitivity for each round?
    print(f"\n  Schedule sensitivity (odd = invertible):")
    for r in range(52, 64):
        invertible = [w for w in range(16) if J[(r,w)] & 1]
        print(f"    W[{r}]: {len(invertible)} invertible words: {invertible[:8]}...")


# ============================================================
# EXP 3: Solve backward → check if solution exists
# ============================================================
def exp3_backward_solve(N):
    print("\n" + "="*70)
    print("EXP 3: BACKWARD SOLVE — fix a[53..64], solve for W[0..15]")
    print("  12 target values → 12 equations on schedule → solve")
    print("="*70)

    for trial in range(N):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = list(W1); W2[0] ^= 0x80000000

        a1, st1, Wexp1, H1 = sha_full_a(W1)
        a2, st2, Wexp2, H2 = sha_full_a(W2)

        # Iterative: for each round 52..63, correct ONE word of W2
        # to make schedule(W2)[r] = required(state2[r], a1[r+1])

        W2_work = list(W2)
        corrections = 0

        for r in range(52, 64):
            # Recompute everything (states change after each correction)
            a2w, st2w, Wexp2w, _ = sha_full_a(W2_work)

            target = a1[r+1]
            W_req = required_W_from_states(st2w[r], target, r)
            W_act = Wexp2w[r]
            mm = sub32(W_req, W_act)

            if mm == 0:
                continue  # Already correct

            # Find word to correct
            best_word = -1; best_residual = 32
            for word in range(16):
                W_test = list(W2_work)
                W_test[word] = add32(W2_work[word], 1)
                Wexp_test = expand_real(W_test)
                sens = sub32(Wexp_test[r], Wexp2w[r])

                if sens & 1:
                    inv = pow(sens, -1, 1<<32)
                    corr = (mm * inv) & MASK32
                    W_corr = list(W2_work)
                    W_corr[word] = add32(W2_work[word], corr)
                    if W_corr == W1: continue

                    Wexp_c = expand_real(W_corr)
                    residual = hw(sub32(Wexp_c[r], W_req))
                    if residual < best_residual:
                        best_residual = residual
                        best_word = word
                        best_W2 = list(W_corr)

            if best_word >= 0:
                W2_work = best_W2
                corrections += 1

        # Check result
        a2f, _, _, H2f = sha_full_a(W2_work)

        # How many a-values match at the end?
        a_match = sum(1 for r in range(53, 65) if a1[r] == a2f[r])
        hash_hw_val = sum(hw(H1[i]^H2f[i]) for i in range(8))
        n_diff = sum(1 for i in range(16) if W1[i] != W2_work[i])

        print(f"  Trial {trial}: corrections={corrections}, "
              f"a_match={a_match}/12, HW(δH)={hash_hw_val}, "
              f"words_diff={n_diff}")

        if a_match >= 8:
            print(f"    ★ {a_match} a-values matched!")
        if hash_hw_val < 100:
            print(f"    ★ HW < 100!")


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    exp1_backward_required(2)
    exp2_free_words()
    exp3_backward_solve(N)
