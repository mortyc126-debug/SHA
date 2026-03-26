#!/usr/bin/env python3
"""
Два потока. Одна машина. Я наблюдаю.

Не оптимизирую. Не измеряю. Наблюдаю и подправляю.

Поток 1: W1 → 64 шага → выход.
Поток 2: W2 → 64 шага → выход.

На каждом шаге я ВИЖУ оба потока.
Я вижу ГДЕ они расходятся и ГДЕ сходятся.
Я подправляю W2 чтобы потоки ТЕКЛИ ВМЕСТЕ.

Не формула. Просто: "вот тут разошлись → подправлю вот это".
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
    """Один шаг машины."""
    a,b,c,d,e,f,g,h = state
    t1 = A(h,S1(e),ch(e,f,g),k,w)
    t2 = A(S0(a),mj(a,b,c))
    return [A(t1,t2),a,b,c,A(d,t1),e,f,g]


def observe(w1, w2):
    """Наблюдаю два потока. Не считаю — ВИЖУ.
    Возвращаю: для каждого раунда, где потоки вместе и где врозь."""
    e1 = expand(w1); e2 = expand(w2)
    s1 = list(IV); s2 = list(IV)

    together = []  # раунды где потоки совпадают
    apart = []     # раунды где расходятся
    divergence = [] # моменты расхождения (было вместе → стало врозь)
    convergence = [] # моменты схождения (было врозь → стало вместе)

    prev_same = True  # начинаем вместе (оба = IV)

    for r in range(64):
        s1 = step(s1, e1[r], C[r])
        s2 = step(s2, e2[r], C[r])

        same = (s1 == s2)

        if same:
            together.append(r)
        else:
            apart.append(r)

        if prev_same and not same:
            divergence.append(r)
        if not prev_same and same:
            convergence.append(r)

        prev_same = same

    # Финал
    final_same = (s1 == s2)
    h1 = tuple(A(IV[i],s1[i]) for i in range(8))
    h2 = tuple(A(IV[i],s2[i]) for i in range(8))
    hash_same = (h1 == h2)

    return {
        'together': together,
        'apart': apart,
        'divergence': divergence,
        'convergence': convergence,
        'final_same': final_same,
        'hash_same': hash_same,
        'hash_dist': sum(HW(h1[i]^h2[i]) for i in range(8)),
    }


def watch_and_adjust(w1, cycles=50):
    """Наблюдаю → подправляю → наблюдаю.

    Моя логика (не формула):
    - Запускаю два потока
    - Смотрю где разошлись
    - Подправляю W2 чтобы потоки ДОЛЬШЕ текли вместе
    - Повторяю

    Подправка: если потоки разошлись на раунде r,
    и r < 16 (W свободен) → поправить W2[r] чтобы
    поток 2 ПОВТОРИЛ поток 1 на этом шаге.
    """

    w2 = list(w1)
    w2[0] = A(w1[0], 1)  # Маленькое отличие

    best_together = 0
    best_hash_dist = 256
    best_w2 = list(w2)

    for cycle in range(cycles):
        obs = observe(w1, w2)

        n_together = len(obs['together'])
        if n_together > best_together:
            best_together = n_together
            best_w2 = list(w2)
        if obs['hash_dist'] < best_hash_dist:
            best_hash_dist = obs['hash_dist']
            best_w2 = list(w2)

        if obs['hash_same']:
            n_diff = sum(1 for i in range(16) if w1[i]!=w2[i])
            if n_diff > 0:
                print(f"  [{cycle}] ★★★ ПОТОКИ СЛИЛИСЬ! hash совпал! diff={n_diff}")
                return w2, True
            else:
                # Тривиально — добавить отличие
                w2[cycle%16] = A(w2[cycle%16], cycle+1)
                continue

        # Подправка: я ВИЖУ первый момент расхождения
        if obs['divergence']:
            r_div = obs['divergence'][0]

            if r_div < 16:
                # Свободный W — я могу поправить ПРЯМО
                # Запускаю поток 1 до раунда r_div, получаю его state
                e1 = expand(w1)
                s1 = list(IV)
                for r in range(r_div):
                    s1 = step(s1, e1[r], C[r])
                # state1 перед раундом r_div
                # T1_1 на раунде r_div
                a1,b1,c1,d1,ee1,f1,g1,h1 = s1
                t1_1 = A(h1,S1(ee1),ch(ee1,f1,g1),C[r_div],e1[r_div])

                # Запускаю поток 2 до раунда r_div
                e2 = expand(w2)
                s2 = list(IV)
                for r in range(r_div):
                    s2 = step(s2, e2[r], C[r])
                a2,b2,c2,d2,ee2,f2,g2,h2 = s2

                # Нужный W2[r_div] чтобы T1₂ = T1₁
                w2_need = D(D(D(D(t1_1,h2),S1(ee2)),ch(ee2,f2,g2)),C[r_div])
                w2[r_div] = w2_need

            else:
                # Несвободный W — не могу поправить напрямую
                # Попробую поправить РАННИЙ W чтобы СДВИНУТЬ schedule
                # Выберу случайный ранний W и слегка изменю
                w_idx = cycle % 16
                w2[w_idx] = A(w2[w_idx], 1)

        if cycle < 10 or cycle % 10 == 0 or n_together > best_together - 2:
            print(f"  [{cycle:3d}] together={n_together:2d}/64, "
                  f"div_at={obs['divergence'][:3] if obs['divergence'] else 'none'}, "
                  f"hash_dist={obs['hash_dist']}")

    return best_w2, False


def deep_flow(N):
    """Запускаю наблюдение для нескольких W1."""
    print("="*60)
    print("НАБЛЮДАЮ ПОТОКИ")
    print("="*60)

    for trial in range(N):
        w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        print(f"\n--- Поток {trial} ---")
        w2, found = watch_and_adjust(w1, cycles=50)

        if found:
            print(f"  COLLISION FOUND!")
            print(f"  W1: {[hex(x) for x in w1[:4]]}...")
            print(f"  W2: {[hex(x) for x in w2[:4]]}...")


if __name__ == '__main__':
    import sys
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    deep_flow(N)
