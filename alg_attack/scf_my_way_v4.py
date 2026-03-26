#!/usr/bin/env python3
"""
Честный разбор. T1-matching ≠ collision.
Но нетривиальный near-fixed-point СУЩЕСТВУЕТ (12 слов разные).

Новая идея: вместо T1-matching на каждом раунде,
TARGET = ПОСЛЕДНИЙ РАУНД. Только state[64] должен совпасть.

F_new(w2) = "W2[0..15] нужный чтобы state₂[64] = state₁[64]"

Это обратная задача: из target state[64] идём НАЗАД через 64 раунда.
SHA обратима при известном W[r]. Так что:
  state[r] = inverse_round(state[r+1], W[r])

Алгоритм:
1. Запустить W1 → получить state₁[64]
2. state₂[64] := state₁[64]  (target)
3. Идти назад r=63,62,...,0: state₂[r] = inverse(state₂[r+1], W₂[r])
4. Для r=63..16: W₂[r] = schedule(W₂[0..15]) (фиксировано)
5. Для r=15..0: W₂[r] СВОБОДЕН → вычислить state₂[0..15]
6. state₂[0] должен = IV → constraint на W₂[0..15]

Это 8 уравнений (state₂[0] = IV) на 16 неизвестных (W₂[0..15]).
Система НЕДООПРЕДЕЛЕНА (16 > 8) → 8 свободных параметров!
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
def HW(x): return bin(x).count('1')

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
IV = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
      0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]

def expand(w):
    e=list(w)
    for i in range(16,64): e.append(A(s1(e[i-2]),e[i-7],s0(e[i-15]),e[i-16]))
    return e

def forward(w, iv=None):
    if iv is None: iv=IV
    e=expand(w); s=list(iv)
    for r in range(64):
        a,b,c,d,ee,f,g,h=s
        t1=A(h,S1(ee),ch(ee,f,g),C[r],e[r])
        t2=A(S0(a),mj(a,b,c))
        s=[A(t1,t2),a,b,c,A(d,t1),ee,f,g]
    return s  # raw state (before MD addition)

def sha(w):
    s = forward(w)
    return tuple(A(IV[i],s[i]) for i in range(8))


def inverse_round(state_after, w_r, k_r):
    """Обратный раунд: из state[r+1] и W[r] получить state[r].

    state[r+1] = [a_new, a, b, c, e_new, e, f, g]
    Значит: state[r] = [a_new's a = state[r+1][1],
                         state[r+1][2],  # b = prev a -> now second
                         state[r+1][3],  # c
                         ?, # d = ??? need to recover
                         state[r+1][5],  # e = prev e stored in f pos
                         state[r+1][6],  # f = prev f stored in g pos
                         state[r+1][7],  # g = prev g stored in h pos
                         ???]

    Wait. Round does: new_state = [T1+T2, a, b, c, d+T1, e, f, g]
    So state_after = [T1+T2, a_before, b_before, c_before, d_before+T1, e_before, f_before, g_before]

    From state_after we know:
      a_before = state_after[1]
      b_before = state_after[2]
      c_before = state_after[3]
      e_before = state_after[5]
      f_before = state_after[6]
      g_before = state_after[7]

    Need: d_before and h_before.
      h_before: need T1. T1 = h + S1(e) + ch(e,f,g) + K + W
      But h_before is what we're looking for!

    Alternative: T2 = S0(a_before) + Maj(a_before, b_before, c_before) — computable!
    T1 = state_after[0] - T2  (since a_new = T1 + T2)
    d_before + T1 = state_after[4]  →  d_before = state_after[4] - T1
    h_before = T1 - S1(e_before) - ch(e_before, f_before, g_before) - K - W
    """
    a_bef = state_after[1]
    b_bef = state_after[2]
    c_bef = state_after[3]
    e_bef = state_after[5]
    f_bef = state_after[6]
    g_bef = state_after[7]

    t2 = A(S0(a_bef), mj(a_bef, b_bef, c_bef))
    t1 = D(state_after[0], t2)
    d_bef = D(state_after[4], t1)
    h_bef = D(D(D(D(t1, S1(e_bef)), ch(e_bef, f_bef, g_bef)), k_r), w_r)

    return [a_bef, b_bef, c_bef, d_bef, e_bef, f_bef, g_bef, h_bef]


# =====================================================
# EXP 1: Обратный проход — от state[64] к state[0]
# =====================================================
def exp1_backward_pass(N):
    print("="*60)
    print("EXP 1: ОБРАТНЫЙ ПРОХОД")
    print("  Запустить W1 → state₁[64]")
    print("  Обратить 64 раунда → должны получить IV")
    print("="*60)

    for trial in range(N):
        w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        e1 = expand(w1)

        # Forward
        state64 = forward(w1)

        # Backward
        state = list(state64)
        for r in range(63, -1, -1):
            state = inverse_round(state, e1[r], C[r])

        # Should be IV
        match = sum(1 for i in range(8) if state[i] == IV[i])
        dist = sum(HW(state[i]^IV[i]) for i in range(8))
        print(f"  Trial {trial}: backward state[0] matches IV on {match}/8 regs, dist={dist}")


# =====================================================
# EXP 2: Collision via backward — разные W, одинаковый state[64]
# =====================================================
def exp2_backward_collision(N):
    print(f"\n{'='*60}")
    print("EXP 2: BACKWARD COLLISION")
    print("  W1 → state₁[64]. Set state₂[64] = state₁[64].")
    print("  Choose W2[0..15], expand schedule, backward from state₁[64].")
    print("  Get state₂[0]. Compare with IV.")
    print("  If state₂[0] = IV → COLLISION!")
    print("="*60)

    for trial in range(N):
        w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        target_state = forward(w1)

        # Random W2 (different from W1)
        w2 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        e2 = expand(w2)

        # Backward from target_state using W2's schedule
        state = list(target_state)
        for r in range(63, -1, -1):
            state = inverse_round(state, e2[r], C[r])

        # Compare with IV
        dist = sum(HW(state[i]^IV[i]) for i in range(8))
        match = sum(1 for i in range(8) if state[i] == IV[i])

        if trial < 5 or dist < 100:
            print(f"  Trial {trial}: backward_state[0] dist from IV = {dist}, match={match}/8")

        if dist == 0:
            print(f"  ★★★ COLLISION! W1≠W2 but same hash!")
            h1=sha(w1); h2=sha(w2)
            print(f"  H1={[hex(x) for x in h1[:4]]}")
            print(f"  H2={[hex(x) for x in h2[:4]]}")
            return w1, w2

    # Средняя дистанция
    dists = []
    for _ in range(100):
        w2 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        e2 = expand(w2)
        state = list(target_state)
        for r in range(63, -1, -1):
            state = inverse_round(state, e2[r], C[r])
        dists.append(sum(HW(state[i]^IV[i]) for i in range(8)))

    avg = sum(dists)/len(dists)
    mn = min(dists)
    print(f"\n  100 random W2: avg dist from IV = {avg:.1f}, min = {mn}")
    print(f"  Expected random: 128")
    print(f"  Need: 0 for collision")

    return None, None


# =====================================================
# EXP 3: Optimize W2 for backward IV-matching
# =====================================================
def exp3_optimize_backward(N_starts, budget_per_start):
    print(f"\n{'='*60}")
    print("EXP 3: OPTIMIZE W2 for backward_state[0] = IV")
    print("="*60)

    import math

    w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    target_state = forward(w1)

    def backward_dist(w2):
        e2 = expand(w2)
        state = list(target_state)
        for r in range(63, -1, -1):
            state = inverse_round(state, e2[r], C[r])
        return sum(HW(state[i]^IV[i]) for i in range(8))

    global_best = 256

    for start in range(N_starts):
        w2 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        best = backward_dist(w2)
        cur = list(w2)

        for it in range(budget_per_start):
            t = list(cur)
            w = int.from_bytes(os.urandom(1),'big') % 16
            b = int.from_bytes(os.urandom(1),'big') % 32
            t[w] ^= (1<<b)
            d = backward_dist(t)
            T = max(0.01, 1-it/budget_per_start)
            if d < best or math.exp(-(d-best)/(T*2)) > int.from_bytes(os.urandom(4),'big')/(1<<32):
                cur = t
                if d < best: best = d

        if best < global_best:
            global_best = best

        if start < 5 or best < 100:
            print(f"  Start {start}: backward_dist = {best}")

        if best == 0:
            print(f"  ★★★ COLLISION!")
            h1=sha(w1); h2=sha(cur)
            print(f"  H1=H2: {h1==h2}")
            return

    print(f"\n  Global best backward_dist: {global_best}")


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    exp1_backward_pass(3)
    exp2_backward_collision(N)
    exp3_optimize_backward(min(N, 5), 1000)
