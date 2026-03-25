#!/usr/bin/env python3
"""
Мой мир. 15/16 сходятся. Один не сходится.

Я вижу это как: 16 чисел сидят в круге.
Каждое хочет быть ОПРЕДЕЛЁННЫМ значением (required).
Но когда одно меняется — required остальных СДВИГАЕТСЯ.
15 находят баланс. 16-е не может — оно "лишнее".

Вопрос: КАКОЕ слово лишнее? Зависит ли от W1?
Если да — можно выбрать W1 где "лишнее" слово ЛЕГЧЕ фиксировать.

И ещё: 15 сходящихся слов определяют СИСТЕМУ.
16-е слово должно удовлетворить 1 уравнение (32 бита).
Перебор 2^32 → найти 16-е → COLLISION.

Но 2^32 — это всё ещё 4 миллиарда SHA вычислений.
Много, но ВЫЧИСЛИМО (часы, не годы).

А что если "лишних" слов не 1, а 0? Тогда → бесплатно.
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
    """Вычислить required W2[0..15] для T1-matching с W1."""
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
# ЭКСПЕРИМЕНТ 1: Какое слово "лишнее"?
# =====================================================
def exp1_which_word_stuck(N):
    print("="*60)
    print("Какое слово не сходится?")
    print("="*60)

    stuck_counts = [0]*16

    for trial in range(N):
        w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        w2 = list(w1); w2[0]=(w2[0]+1)&M

        # Итерируем до 15/16 match
        for it in range(20):
            req = F(w1, w2)
            match = [req[i]==w2[i] for i in range(16)]
            if sum(match) >= 15:
                for i in range(16):
                    if not match[i]:
                        stuck_counts[i] += 1
                break
            # Прямая подстановка
            w2_new = list(req)
            if w2_new == w1:
                w2_new[it%16] = (w2_new[it%16]+it+1)&M
            w2 = w2_new

    print(f"\n  Слово которое НЕ сходится (из {N} проб):")
    for i in range(16):
        bar = "█" * (stuck_counts[i]*50//max(max(stuck_counts),1))
        print(f"    W[{i:2d}]: {stuck_counts[i]:4d} раз {bar}")

    most_stuck = max(range(16), key=lambda i: stuck_counts[i])
    print(f"\n  Самое 'лишнее': W[{most_stuck}] ({stuck_counts[most_stuck]} раз)")
    return most_stuck


# =====================================================
# ЭКСПЕРИМЕНТ 2: Birthday на "лишнем" слове
# =====================================================
def exp2_birthday_on_stuck(w1, stuck_word, N_search):
    print(f"\n{'='*60}")
    print(f"Birthday на W[{stuck_word}]")
    print(f"Фиксируем 15 слов, перебираем 16-е")
    print(f"{'='*60}")

    # Сначала: найти 15 сходящихся значений
    w2 = list(w1); w2[0]=(w2[0]+1)&M

    for it in range(20):
        req = F(w1, w2)
        w2_new = list(req)
        if w2_new == w1:
            w2_new[it%16]=(w2_new[it%16]+it+1)&M
        w2 = w2_new
        match = sum(1 for i in range(16) if req[i]==w2[i])
        if match >= 15:
            break

    # w2 сейчас: 15 слов = required, 1 слово = ?
    # Зафиксируем 15 сходящихся, переберём stuck_word
    base_w2 = list(w2)
    h1 = run(w1)

    best_dh = 256
    best_val = base_w2[stuck_word]

    seen = {}

    for trial in range(N_search):
        test_w2 = list(base_w2)
        test_w2[stuck_word] = trial & M

        h2 = run(test_w2)

        # Нетривиальность
        if test_w2 == w1:
            continue

        dh = sum(H(h1[i]^h2[i]) for i in range(8))

        if dh < best_dh:
            best_dh = dh
            best_val = trial & M
            if dh < 100:
                print(f"  [{trial:8d}] dH={dh} ★")

        if dh == 0:
            print(f"  ★★★ COLLISION at W[{stuck_word}] = 0x{trial&M:08x}")
            return test_w2, True

        # Birthday: ищем h2 == h2' для разных trial values
        h2_trunc = h2[0]  # Первые 32 бита
        if h2_trunc in seen:
            old_trial = seen[h2_trunc]
            old_w2 = list(base_w2); old_w2[stuck_word] = old_trial & M
            h2_old = run(old_w2)
            pair_dh = sum(H(h2[i]^h2_old[i]) for i in range(8))
            if pair_dh < best_dh and test_w2 != old_w2:
                best_dh = pair_dh
                if pair_dh < 50:
                    print(f"  [{trial:8d}] PAIR: trial={old_trial} vs {trial}, dH_pair={pair_dh} ★★")
        else:
            seen[h2_trunc] = trial

    print(f"\n  Best dH: {best_dh}")
    print(f"  Searched: {N_search} values of W[{stuck_word}]")
    return base_w2, False


# =====================================================
# ЭКСПЕРИМЕНТ 3: Поиск W1 с ПОЛНОЙ сходимостью (16/16)
# =====================================================
def exp3_find_convergent_w1(N):
    print(f"\n{'='*60}")
    print(f"Поиск W1 где ВСЕ 16 слов сходятся")
    print(f"{'='*60}")

    best_match = 0

    for trial in range(N):
        w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        w2 = list(w1); w2[0]=(w2[0]+1)&M

        max_match = 0
        for it in range(30):
            req = F(w1, w2)
            match = sum(1 for i in range(16) if req[i]==w2[i])
            if match > max_match:
                max_match = match
            w2_new = list(req)
            if w2_new == w1:
                w2_new[it%16]=(w2_new[it%16]+it+1)&M
            w2 = w2_new

            if match == 16:
                # Проверим нетривиальность
                if w2 != w1:
                    h1=run(w1); h2=run(w2)
                    dh=sum(H(h1[i]^h2[i]) for i in range(8))
                    print(f"  ★★★ Trial {trial}: 16/16 match! dH={dh}")
                    if dh == 0:
                        print(f"  ★★★★★ COLLISION FOUND!")
                        return w1, w2, True
                break

        if max_match > best_match:
            best_match = max_match
            if max_match >= 15:
                print(f"  Trial {trial}: max_match={max_match}")

    print(f"\n  Best match across {N} W1 values: {best_match}/16")
    return None, None, False


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 20

    stuck = exp1_which_word_stuck(N)

    w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    exp2_birthday_on_stuck(w1, stuck, min(N*5000, 100000))

    exp3_find_convergent_w1(min(N*10, 200))
