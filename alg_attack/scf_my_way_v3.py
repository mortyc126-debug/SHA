#!/usr/bin/env python3
"""
15/16 сходятся. W[15] застревает. Но КУДА сходятся 15?
К W1 (тривиально) или к ДРУГИМ значениям (нетривиально)?

Я проверю: после итерации, W2[0..14] = W1[0..14]?
Если да — тривиально (как раньше).
Если нет — у нас 15 НЕТРИВИАЛЬНЫХ фиксированных точек.
"""
import os, sys

M = 0xFFFFFFFF

def R(x,n): return ((x>>n)|(x<<(32-n)))&M
def S0(x): return R(x,2)^R(x,13)^R(x,22)
def S1(x): return R(x,6)^R(x,11)^R(x,25)
def s0(x): return R(x,7)^R(x,18)^(x>>3)
def s1(x): return R(x,17)^R(x,19)^(x>>10)
def ch(e,f,g): return (e&f)^(~e&g)&M
def mj(a,b,c): return (a&b)^(a&c)^(b&c)
def A(*a):
    s=0
    for x in a: s=(s+x)&M
    return s
def D(a,b): return (a-b)&M
def H(x): return bin(x).count('1')

C = [
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2,
]
I = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
     0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]

def expand(w):
    e=list(w)
    for i in range(16,64): e.append(A(s1(e[i-2]),e[i-7],s0(e[i-15]),e[i-16]))
    return e

def run(w):
    e=expand(w); s=list(I)
    for r in range(64):
        a,b,c,d,ee,f,g,h=s
        t1=A(h,S1(ee),ch(ee,f,g),C[r],e[r])
        t2=A(S0(a),mj(a,b,c))
        s=[A(t1,t2),a,b,c,A(d,t1),ee,f,g]
    return tuple(A(I[i],s[i]) for i in range(8))

def F(w1, w2):
    e1=expand(w1); e2=expand(w2)
    s1=list(I); s2=list(I)
    req=[0]*16
    for r in range(64):
        a1,b1,c1,d1,e1r,f1,g1,h1=s1
        a2,b2,c2,d2,e2r,f2,g2,h2=s2
        t1_1=A(h1,S1(e1r),ch(e1r,f1,g1),C[r],e1[r])
        w2_need=D(D(D(D(t1_1,h2),S1(e2r)),ch(e2r,f2,g2)),C[r])
        if r<16: req[r]=w2_need
        t1_2=A(h2,S1(e2r),ch(e2r,f2,g2),C[r],e2[r])
        t2_1=A(S0(a1),mj(a1,b1,c1))
        t2_2=A(S0(a2),mj(a2,b2,c2))
        s1=[A(t1_1,t2_1),a1,b1,c1,A(d1,t1_1),e1r,f1,g1]
        s2=[A(t1_2,t2_2),a2,b2,c2,A(d2,t1_2),e2r,f2,g2]
    return req


# =====================================================
# ПРОВЕРКА: тривиальная или нетривиальная сходимость?
# =====================================================
def check_trivial(N):
    print("="*60)
    print("ПРОВЕРКА: W2[0..14] = W1[0..14]?")
    print("  Если да → тривиально (W2 = W1)")
    print("  Если нет → НЕТРИВИАЛЬНАЯ фиксированная точка!")
    print("="*60)

    for trial in range(N):
        w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        w2 = list(w1); w2[0]=(w2[0]+1)&M

        # Итерируем до сходимости
        for it in range(20):
            req = F(w1, w2)
            w2_new = list(req)
            if w2_new == w1:
                w2_new[it%16]=(w2_new[it%16]+it+1)&M
            w2 = w2_new

        # Проверяем: W2[0..14] = W1[0..14]?
        same = [w2[i] == w1[i] for i in range(16)]
        n_same = sum(same)
        n_converged_to_req = sum(1 for i in range(16) if w2[i] == req[i])

        # Ключевой вопрос: converged W2[i] = W1[i]?
        converged_trivial = sum(1 for i in range(15) if w2[i] == w1[i])
        converged_nontrivial = 15 - converged_trivial

        diff_words = [(i, w1[i], w2[i]) for i in range(16) if w1[i] != w2[i]]

        if trial < 10:
            print(f"\n  Trial {trial}:")
            print(f"    W2 = W1 on {n_same}/16 words")
            print(f"    W2 = required on {n_converged_to_req}/16 words")
            print(f"    Converged AND trivial (=W1): {converged_trivial}/15")
            print(f"    Converged AND nontrivial:    {converged_nontrivial}/15")
            if diff_words:
                print(f"    Different words:")
                for i, v1, v2 in diff_words[:5]:
                    print(f"      W[{i:2d}]: W1=0x{v1:08x} W2=0x{v2:08x} δ=0x{v1^v2:08x} HW={H(v1^v2)}")

            if converged_nontrivial > 0:
                print(f"    ★★★ {converged_nontrivial} NONTRIVIAL converged words!")
                # Verify hash
                h1=run(w1); h2=run(w2)
                dh=sum(H(h1[i]^h2[i]) for i in range(8))
                print(f"    Hash diff: HW={dh}")


# =====================================================
# ЧИСТЫЙ ТЕСТ: запустить F(w1, w2) когда w2 отличается ТОЛЬКО W[0]
# =====================================================
def clean_test(N):
    print(f"\n{'='*60}")
    print("ЧИСТЫЙ ТЕСТ: w2 = w1 кроме W[0]")
    print("  Если F(w1,w2)[r] = w1[r] для r=1..14 → тривиально")
    print("  Если F(w1,w2)[r] ≠ w1[r] для какого-то r → нетривиально!")
    print("="*60)

    for trial in range(N):
        w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        w2 = list(w1); w2[0] = (w1[0] + 1) & M  # Только W[0] отличается

        req = F(w1, w2)  # Одна итерация, не цикл

        # Проверяем: req[r] = w1[r] для r=0..15?
        matches_w1 = [req[i] == w1[i] for i in range(16)]
        matches_w2 = [req[i] == w2[i] for i in range(16)]

        if trial < 10:
            print(f"\n  Trial {trial}:")
            for i in range(16):
                m1 = "=W1" if matches_w1[i] else "≠W1"
                m2 = "=W2" if matches_w2[i] else "≠W2"
                delta = H(req[i] ^ w1[i])
                status = ""
                if matches_w1[i]: status = "→ TRIVIAL"
                elif matches_w2[i]: status = "→ STAYS"
                else: status = f"→ NEW (δ from W1: {delta} bits)"
                print(f"      req[{i:2d}] {m1} {m2} {status}")


# =====================================================
# ГЛУБОКИЙ ТЕСТ: что если начать с W2 ДАЛЕКО от W1?
# =====================================================
def deep_test(N):
    print(f"\n{'='*60}")
    print("ГЛУБОКИЙ ТЕСТ: W2 ДАЛЕКО от W1")
    print("  Если F всё равно сходится к W1 → тривиально навсегда")
    print("  Если F идёт к ДРУГОЙ точке → нетривиальный аттрактор!")
    print("="*60)

    for trial in range(N):
        w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        # W2 ПОЛНОСТЬЮ случайный (далеко от W1)
        w2 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

        for it in range(30):
            req = F(w1, w2)
            match_w1 = sum(1 for i in range(16) if req[i] == w1[i])
            match_req = sum(1 for i in range(16) if req[i] == w2[i])

            if it < 5 or it % 10 == 0:
                print(f"  [{it:2d}] req=W1: {match_w1}/16, req=W2: {match_req}/16, "
                      f"Σδ(req,W1)={sum(H(req[i]^w1[i]) for i in range(16))}")

            w2 = list(req)

        # Итог: к чему сошлось?
        final_dist_w1 = sum(H(w2[i]^w1[i]) for i in range(16))
        h1=run(w1); h2=run(w2)
        dh=sum(H(h1[i]^h2[i]) for i in range(8))
        is_trivial = w2 == w1

        print(f"  → Final: dist_from_W1={final_dist_w1}, dH={dh}, trivial={is_trivial}")
        if not is_trivial and dh < 128:
            print(f"  ★ NONTRIVIAL ATTRACTOR! dist={final_dist_w1}, dH={dh}")


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    clean_test(N)
    check_trivial(N)
    deep_test(min(N, 3))
