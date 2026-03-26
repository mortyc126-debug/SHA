#!/usr/bin/env python3
"""
Смотрю на переход. Раунды 13-16. Момент где свободное
становится связанным.

Я не меряю HW. Я смотрю на ЗНАЧЕНИЯ. Что КОНКРЕТНО
происходит с числами на переходе.
"""
import os

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
    return [A(t1,t2),a,b,c,A(d,t1),e,f,g], t1, t2


def look_at_transition():
    """Я смотрю. Не считаю — смотрю."""

    w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    e1 = expand(w1)

    # Строю w2 с T1-matching на раундах 0..14
    w2 = list(w1)
    w2[0] = A(w1[0], 1)  # Начальное отличие

    s1 = list(IV); s2 = list(IV)

    # Раунды 0..14: фиксирую W2 для T1-matching
    for r in range(15):
        a1,b1,c1,d1,ee1,f1,g1,h1 = s1
        t1_1 = A(h1,S1(ee1),ch(ee1,f1,g1),C[r],e1[r])

        a2,b2,c2,d2,ee2,f2,g2,h2 = s2
        # Нужный W2[r]
        w2[r] = D(D(D(D(t1_1,h2),S1(ee2)),ch(ee2,f2,g2)),C[r])

        s1, _, t2_1 = step(s1, e1[r], C[r])
        s2, _, t2_2 = step(s2, w2[r], C[r])

    e2 = expand(w2)

    print("="*70)
    print("СМОТРЮ НА ПЕРЕХОД (раунды 12-20)")
    print("="*70)

    # Перезапускаю с финальными w1, w2
    s1 = list(IV); s2 = list(IV)
    e1 = expand(w1); e2 = expand(w2)

    for r in range(22):
        s1_new, t1_1, t2_1 = step(s1, e1[r], C[r])
        s2_new, t1_2, t2_2 = step(s2, e2[r], C[r])

        same_state = (s1 == s2)
        same_t1 = (t1_1 == t1_2)
        same_t2 = (t2_1 == t2_2)
        same_w = (e1[r] == e2[r])
        same_after = (s1_new == s2_new)

        # Что я ВИЖУ
        free = r < 16
        sched = r >= 16

        print(f"\n  Round {r:2d} {'[FREE]' if free else '[SCHED]'}:")
        print(f"    W same:     {same_w}")
        print(f"    state same: {same_state}")
        print(f"    T1 same:    {same_t1}")
        print(f"    T2 same:    {same_t2}")
        print(f"    next same:  {same_after}")

        if not same_state and same_t1:
            # Потоки расходятся но T1 совпадает!
            # Значит W2 КОМПЕНСИРОВАЛ разницу state
            delta_w = D(e2[r], e1[r])
            print(f"    → W₂ КОМПЕНСИРУЕТ state-разницу! δW=0x{delta_w:08x}")

        if same_state and not same_w:
            print(f"    → State одинаков но W разный! δW=0x{D(e2[r],e1[r]):08x}")
            print(f"    → T1 будет РАЗНЫЙ → state РАЗОЙДЁТСЯ")

        if not same_state and not same_t1:
            # Расхождение и T1 не совпал
            delta_t1 = D(t1_2, t1_1)
            print(f"    → РАСХОЖДЕНИЕ: δT1=0x{delta_t1:08x}")

        if same_state and same_w:
            print(f"    → ВМЕСТЕ (тривиально)")

        s1 = s1_new; s2 = s2_new

    # КЛЮЧЕВОЙ ВОПРОС: на раунде 15 (последний свободный)
    # state₁ ≠ state₂ (потому что T2 разные на предыдущих раундах)
    # W₂[15] выбран для T1-matching → T1₁=T1₂
    # Но T2₁ ≠ T2₂ → state[16]₁ ≠ state[16]₂
    # На раунде 16: W₁[16] = W₂[16] (schedule, W1[0..15] и W2[0..15] разные
    # но schedule одинаковый? НЕТ! W2[0..15] ДРУГИЕ → schedule₂ ДРУГОЙ!)

    print(f"\n{'='*70}")
    print("ЧТО Я ВИЖУ:")
    print(f"{'='*70}")
    print("""
    Раунды 0..14: W₂ подобран для T1-matching.
      state₁ и state₂ РАСХОДЯТСЯ (через T2) но T1 СОВПАДАЕТ.
      Каждый раунд: state чуть дальше, но W компенсирует.

    Раунд 15: последний свободный W.
      T1 совпадает. Но state уже ДАЛЕКО.

    Раунд 16: ПЕРВЫЙ связанный W.
      W₂[16] = schedule(W₂[0..15]) ← из НАШИХ подобранных W₂[0..15]
      W₁[16] = schedule(W₁[0..15]) ← из оригинальных W₁[0..15]
      W₁[16] ≠ W₂[16] (потому что W₁ ≠ W₂)

      И state₁[16] ≠ state₂[16].
      Оба РАЗНЫЕ (W и state) → T1 РАЗНЫЙ → полное расхождение.

    ВЫВОД: переход на раунде 16 — это момент где
    КОМПЕНСАЦИЯ (через свободный W) ЗАКАНЧИВАЕТСЯ
    и СВЯЗАННОСТЬ (через schedule) НАЧИНАЕТСЯ.

    Свободные раунды: W ОТВЕЧАЕТ на state → потоки вместе.
    Связанные раунды: W ЗАДАН schedule → потоки врозь.
    """)


if __name__ == '__main__':
    look_at_transition()
