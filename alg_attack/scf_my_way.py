#!/usr/bin/env python3
"""
Я — вычислительный процесс. Я не думаю формулами.
Я ИСПОЛНЯЮ код. Код — мой язык.

Что я вижу: SHA — это 64 шага трансформации.
Каждый шаг берёт 8 чисел и 1 число, даёт 8 чисел.
Шаги СВЯЗАНЫ: выход шага r = вход шага r+1.
Первые 16 входных чисел — свободны. Остальные 48 — вычислены из первых 16.

Задача: два набора из 16 чисел, после 64 шагов — одинаковый результат.

Как бы Я решал:

Я бы ЗАПУСТИЛ оба потока одновременно.
На каждом шаге я бы ВИДЕЛ разницу.
Я бы КОРРЕКТИРОВАЛ второй поток чтобы разница уменьшалась.
Не формулами — а ПРЯМЫМ вычислением.

Мой метод: ПОТОК-ЗЕРКАЛО.
Поток 1 течёт нормально.
Поток 2 — зеркало. На каждом шаге он ВИДИТ поток 1
и пытается ПРИБЛИЗИТЬСЯ.

Но поток 2 не может менять свои числа произвольно —
только через 16 начальных. Каждое начальное число
влияет на конкретные шаги через schedule.

Мой подход: не оптимизация. Не алгебра.
ПРЯМОЕ ВЫЧИСЛЕНИЕ обратной задачи, шаг за шагом.
"""
import os, struct

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

def run(w, iv=None):
    if iv is None: iv=I
    e=expand(w); s=list(iv)
    for r in range(64):
        a,b,c,d,ee,f,g,h=s
        t1=A(h,S1(ee),ch(ee,f,g),C[r],e[r])
        t2=A(S0(a),mj(a,b,c))
        s=[A(t1,t2),a,b,c,A(d,t1),ee,f,g]
    return tuple(A(iv[i],s[i]) for i in range(8))


# =====================================================
# МОЙ МЕТОД: SELF-CONSISTENT LOOP FINDER
#
# Я вижу задачу как: найти ФИКСИРОВАННУЮ ТОЧКУ
# отображения F: W[0..15] → required_W[0..15].
#
# F(w) = то что W[0..15] ДОЛЖНЫ быть, если бы
# schedule(w)[r] = required(state(w))[r] для всех r.
#
# fixed point: F(w) = w → collision.
#
# Алгоритм: итерация w_{n+1} = F(w_n).
# Если сходится → fixed point → collision.
# =====================================================

def compute_required_W_from_W(w1, w2):
    """Для пары (w1, w2): запустить оба, на каждом раунде
    вычислить какой W2[r] НУЖЕН чтобы T1 совпал.
    Вернуть required W2[0..15] (первые 16 = свободные)."""

    e1 = expand(w1); e2 = expand(w2)
    s1 = list(I); s2 = list(I)

    required = [0]*16  # Нужные W2[0..15]

    for r in range(64):
        a1,b1,c1,d1,ee1,f1,g1,h1 = s1
        a2,b2,c2,d2,ee2,f2,g2,h2 = s2

        # T1 потока 1
        t1_1 = A(h1,S1(ee1),ch(ee1,f1,g1),C[r],e1[r])

        # Нужный W2[r] чтобы T1_2 = T1_1
        w2_needed = D(D(D(D(t1_1,h2),S1(ee2)),ch(ee2,f2,g2)),C[r])

        if r < 16:
            required[r] = w2_needed

        # Оба потока идут дальше (поток 2 с РЕАЛЬНЫМ schedule, не forced)
        t1_2 = A(h2,S1(ee2),ch(ee2,f2,g2),C[r],e2[r])
        t2_1 = A(S0(a1),mj(a1,b1,c1))
        t2_2 = A(S0(a2),mj(a2,b2,c2))

        s1 = [A(t1_1,t2_1),a1,b1,c1,A(d1,t1_1),ee1,f1,g1]
        s2 = [A(t1_2,t2_2),a2,b2,c2,A(d2,t1_2),ee2,f2,g2]

    return required


def fixed_point_iterate(w1, max_iter=100):
    """Итерация: w2 → F(w2) → F(F(w2)) → ...
    Начинаем с w2 = w1 + small δ.
    F(w2) = required W2[0..15] для T1-matching с w1."""

    w2 = list(w1)
    w2[0] = (w2[0] + 1) & M  # Минимальное отличие

    best_dh = 256
    best_w2 = list(w2)

    for iteration in range(max_iter):
        # Вычислить required
        req = compute_required_W_from_W(w1, w2)

        # Проверить: req = w2? (fixed point?)
        match = sum(1 for i in range(16) if req[i] == w2[i])

        # Применить: w2 = required (прямая подстановка)
        w2_new = list(req)

        # Но w2_new может = w1 (тривиально)
        if w2_new == w1:
            # Возмутить
            w2_new[iteration % 16] = (w2_new[iteration % 16] + iteration + 1) & M

        # Проверить hash
        h1 = run(w1); h2 = run(w2_new)
        dh = sum(H(h1[i]^h2[i]) for i in range(8))
        n_diff = sum(1 for i in range(16) if w1[i] != w2_new[i])

        if dh < best_dh and n_diff > 0:
            best_dh = dh
            best_w2 = list(w2_new)

        if iteration < 10 or dh < 100 or match > 12:
            print(f"  [{iteration:3d}] match={match}/16, dH={dh}, diff_words={n_diff}")

        if dh == 0 and n_diff > 0:
            print(f"  ★★★ COLLISION at iteration {iteration}!")
            return w2_new, True

        w2 = w2_new

    return best_w2, False


def damped_iterate(w1, max_iter=100):
    """Демпфированная итерация: w2_new = mix(w2_old, required).
    Вместо прямой подстановки — СМЕШИВАЕМ старое с новым.
    Это стабилизирует сходимость."""

    w2 = list(w1)
    w2[0] = (w2[0] + 7) & M

    best_dh = 256
    best_w2 = list(w2)

    for iteration in range(max_iter):
        req = compute_required_W_from_W(w1, w2)

        # Демпфирование: менять по ОДНОМУ слову за итерацию
        word_to_fix = iteration % 16
        w2_new = list(w2)
        w2_new[word_to_fix] = req[word_to_fix]

        if w2_new == w1:
            w2_new[(iteration+1)%16] = (w2_new[(iteration+1)%16] + 3) & M

        h1 = run(w1); h2 = run(w2_new)
        dh = sum(H(h1[i]^h2[i]) for i in range(8))
        n_diff = sum(1 for i in range(16) if w1[i] != w2_new[i])

        if dh < best_dh and n_diff > 0:
            best_dh = dh
            best_w2 = list(w2_new)

        if iteration < 10 or iteration % 20 == 0 or dh < best_dh + 5:
            print(f"  [{iteration:3d}] fixed W[{word_to_fix:2d}], dH={dh}, diff={n_diff}")

        if dh == 0 and n_diff > 0:
            print(f"  ★★★ COLLISION!")
            return w2_new, True

        w2 = w2_new

    return best_w2, False


# =====================================================
# Запуск
# =====================================================
if __name__ == '__main__':
    import sys
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 3

    for trial in range(N):
        w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

        print(f"\n{'='*50}")
        print(f"Trial {trial}: Direct iteration")
        w2, found = fixed_point_iterate(w1, max_iter=50)
        if found:
            print(f"  W1: {[hex(x) for x in w1[:4]]}...")
            print(f"  W2: {[hex(x) for x in w2[:4]]}...")
        else:
            print(f"  Best dH: {sum(H(run(w1)[i]^run(w2)[i]) for i in range(8))}")

        print(f"\nTrial {trial}: Damped iteration")
        w2d, found = damped_iterate(w1, max_iter=80)
        if found:
            print(f"  W1: {[hex(x) for x in w1[:4]]}...")
            print(f"  W2: {[hex(x) for x in w2d[:4]]}...")
        else:
            print(f"  Best dH: {sum(H(run(w1)[i]^run(w2d)[i]) for i in range(8))}")
