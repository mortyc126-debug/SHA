#!/usr/bin/env python3
"""
ТРИ ПОТОКА НАВСТРЕЧУ: вперёд + назад + середина.

Поток 1 (вперёд): break r=3, a-repair → BOTH=0 на r=1-3,12-17.
Поток 2 (назад): обратный SHA от target state[64].
  SHA обратима! state[r] = inv_round(state[r+1], W[r]).
  Backward a-repair: фиксируем a-sequence с конца.

Поток 3 (середина): свободный state[32]. Мост.

Я построю: forward a-repair (r=0→16) + backward a-repair (r=63→48).
Измерю: сколько BOTH=0 раундов с каждой стороны.
Потом: gap = оставшиеся раунды в середине.
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

def step(state, w, k):
    a,b,c,d,e,f,g,h = state
    t1 = A(h,S1(e),ch(e,f,g),k,w)
    t2 = A(S0(a),mj(a,b,c))
    return [A(t1,t2),a,b,c,A(d,t1),e,f,g]

def inv_step(state_after, w, k):
    """Обратный раунд."""
    a_bef = state_after[1]
    b_bef = state_after[2]
    c_bef = state_after[3]
    e_bef = state_after[5]
    f_bef = state_after[6]
    g_bef = state_after[7]
    t2 = A(S0(a_bef), mj(a_bef,b_bef,c_bef))
    t1 = D(state_after[0], t2)
    d_bef = D(state_after[4], t1)
    h_bef = D(D(D(D(t1, S1(e_bef)), ch(e_bef,f_bef,g_bef)), k), w)
    return [a_bef,b_bef,c_bef,d_bef,e_bef,f_bef,g_bef,h_bef]

def sha(w):
    e=expand(w); s=list(IV)
    for r in range(64):
        s = step(s, e[r], C[r])
    return tuple(A(IV[i],s[i]) for i in range(8))


def three_way(w1):
    """Три потока навстречу."""
    e1 = expand(w1)

    # Поток 1 вперёд: все state₁
    states1 = [list(IV)]
    s = list(IV)
    for r in range(64):
        s = step(s, e1[r], C[r])
        states1.append(list(s))

    # TARGET: state₁[64] (для collision hash)
    target = list(states1[64])

    # === FORWARD a-repair (break r=3) ===
    w2_fwd = list(w1)
    w2_fwd[3] = A(w1[3], 1)  # break

    s2_fwd = list(IV)
    for r in range(3):
        s2_fwd = step(s2_fwd, e1[r], C[r])
    # Break round
    s2_fwd = step(s2_fwd, w2_fwd[3], C[3])
    # a-repair r=4..15
    for r in range(4, 16):
        a2,b2,c2,d2,e2,f2,g2,h2 = s2_fwd
        a_target = states1[r+1][0]
        t2_2 = A(S0(a2),mj(a2,b2,c2))
        t1_need = D(a_target, t2_2)
        w2_fwd[r] = D(D(D(D(t1_need,h2),S1(e2)),ch(e2,f2,g2)),C[r])
        s2_fwd = step(s2_fwd, w2_fwd[r], C[r])

    # Forward BOTH=0 check
    e2_fwd = expand(w2_fwd)
    s1c=list(IV); s2c=list(IV)
    fwd_both0 = []
    fwd_states2 = [list(IV)]
    for r in range(64):
        s1c = step(s1c, e1[r], C[r])
        s2c = step(s2c, e2_fwd[r], C[r])
        fwd_states2.append(list(s2c))
        if s1c == s2c:
            fwd_both0.append(r+1)

    # Последний forward BOTH=0 round
    fwd_reach = max(fwd_both0) if fwd_both0 else 0

    # === BACKWARD a-repair ===
    # Идём от target state[64] назад.
    # Backward: state[r] = inv_step(state[r+1], W[r]).
    # W₁[r] для r≥16 = schedule (фиксирован).
    # Для backward a-repair: используем W₁[r] и проверяем
    # сколько раундов state₂ = state₁ при обратном проходе.

    # Backward: если state₂[64] = state₁[64] (target),
    # и W₂ = W₁ для r=63..R_back, то state₂[R_back] = state₁[R_back].
    # Это тривиально! Потому что одинаковый W + одинаковый state → одинаковый.

    # ДЛЯ НЕТРИВИАЛЬНОСТИ: backward break на каком-то раунде.
    # Break backward: W₂[60] ≠ W₁[60]. Потом backward a-repair 59→48.

    # Но W₂[60] для r≥16 = schedule(W₂[0..15]). W₂[0..15] уже
    # ОПРЕДЕЛЕНЫ forward a-repair! Значит W₂[60] ФИКСИРОВАН.
    # Нельзя делать backward break на schedule rounds!

    # ОДНАКО: W₂[0..15] определяют schedule₂.
    # schedule₂[r] ≠ schedule₁[r] для r≥16.
    # Значит: backward от target с schedule₂ АВТОМАТИЧЕСКИ
    # даёт другую backward trajectory.

    # Проверю: сколько раундов BACKWARD state₂ = state₁?
    # (идя от r=64 назад с schedule₂)

    back_both0 = []
    s1_back = list(states1[64])
    s2_back = list(states1[64])  # Start from SAME target!

    for r in range(63, -1, -1):
        s1_back = inv_step(s1_back, e1[r], C[r])
        s2_back = inv_step(s2_back, e2_fwd[r], C[r])
        if s1_back == s2_back:
            back_both0.append(r)

    back_reach = min(back_both0) if back_both0 else 64

    # Gap
    gap = back_reach - fwd_reach if back_reach > fwd_reach else 0

    h1=sha(w1); h2=sha(w2_fwd)
    dh=sum(HW(h1[i]^h2[i]) for i in range(8))

    return fwd_both0, back_both0, fwd_reach, back_reach, gap, dh


def run(N):
    print("="*60)
    print("THREE-WAY: forward + backward a-repair")
    print("="*60)

    for trial in range(N):
        w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        fwd, bwd, fr, br, gap, dh = three_way(w1)

        print(f"\n  Trial {trial}:")
        print(f"    Forward BOTH=0: {fwd[:15]} (reach r={fr})")
        print(f"    Backward BOTH=0: {sorted(bwd)[:15]} (reach r={br})")
        print(f"    GAP: r={fr}..{br} = {gap} rounds")
        print(f"    dH = {dh}")

        if gap < 30:
            print(f"    ★ Gap < 30!")
        if len(fwd) + len(bwd) > 20:
            print(f"    ★ {len(fwd)+len(bwd)} total BOTH=0 rounds!")


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    run(N)
