#!/usr/bin/env python3
"""
BUDGET ACCOUNTING: 16 слов, 16 задач.

Break: 1 слово (δW[0] создаёт δstate)
Setup: 3 слова (W[1..3] подготавливают state для repair)
a-repair: 12 слов (W[4..15] обнуляют δa[5..16])
  → shift register: δa=0 на 12 раундов → δb,δc,δd=0 автоматически
  → D-coupling: δa[5..16]=0 → δe[9..20]=0 (через 4 раунда задержки)

Schedule: W[16..63] определены из W[0..15].
  После r=16: δstate зависит от schedule.
  Нужно: schedule СОВМЕСТИМ с продолжением repair.

Вопрос: если W[0..15] потрачены на break+setup+repair для r=0..15,
может ли schedule[16..63] ПРОДОЛЖИТЬ repair?

Это вопрос ИНЕРЦИИ: repair создал δa=0 на 12 раундах.
Shift register несёт эти нули дальше: δb,δc,δd=0 на r=16..19.
Если δa[16]=0 тоже → δT2[17]=0 → δa управляемо через W[17].
Но W[17] = schedule → не контролируем.

ПРОВЕРИМ: при a-repair на r=0..15, какова δa[16..20]?
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


def break_and_repair(w1, break_word=0, break_delta=1, repair_target='a'):
    """Break на слове break_word, repair target ('a' или 'e') на r=1..15."""

    e1 = expand(w1)

    # Поток 1: полностью
    s1_all = [list(IV)]
    s = list(IV)
    for r in range(64):
        a,b,c,d,ee,f,g,h = s
        t1=A(h,S1(ee),ch(ee,f,g),C[r],e1[r])
        t2=A(S0(a),mj(a,b,c))
        s=[A(t1,t2),a,b,c,A(d,t1),ee,f,g]
        s1_all.append(list(s))

    # Поток 2: break на раунде 0, repair на 1..15
    w2 = list(w1)
    w2[break_word] = A(w1[break_word], break_delta)

    s2 = list(IV)

    # Раунд 0: break (W₂[0] ≠ W₁[0])
    a2,b2,c2,d2,e2,f2,g2,h2 = s2
    t1_2 = A(h2,S1(e2),ch(e2,f2,g2),C[0],w2[0])
    t2_2 = A(S0(a2),mj(a2,b2,c2))
    s2 = [A(t1_2,t2_2),a2,b2,c2,A(d2,t1_2),e2,f2,g2]

    # state₂[1] ≠ state₁[1] — break создан

    # Раунды 1..15: repair
    for r in range(1, 16):
        a2,b2,c2,d2,e2,f2,g2,h2 = s2

        if repair_target == 'a':
            # a-repair: a₂[r+1] = a₁[r+1]
            a_target = s1_all[r+1][0]
            t2_2 = A(S0(a2),mj(a2,b2,c2))
            t1_need = D(a_target, t2_2)
            w2[r] = D(D(D(D(t1_need,h2),S1(e2)),ch(e2,f2,g2)),C[r])

        elif repair_target == 'e':
            # e-repair (Wang): T1₂ = T1₁
            a1,b1,c1,d1,ee1,f1,g1,h1 = s1_all[r]
            t1_1 = A(h1,S1(ee1),ch(ee1,f1,g1),C[r],e1[r])
            w2[r] = D(D(D(D(t1_1,h2),S1(e2)),ch(e2,f2,g2)),C[r])

        # Step stream 2
        t1_2 = A(h2,S1(e2),ch(e2,f2,g2),C[r],w2[r])
        t2_2 = A(S0(a2),mj(a2,b2,c2))
        s2 = [A(t1_2,t2_2),a2,b2,c2,A(d2,t1_2),e2,f2,g2]

    # Теперь: прогоняю оба потока полностью с финальными w2
    e2_full = expand(w2)
    s2_check = list(IV)
    s1_check = list(IV)

    results = []
    for r in range(64):
        a1,b1,c1,d1,ee1,f1,g1,h1 = s1_check
        t1_1=A(h1,S1(ee1),ch(ee1,f1,g1),C[r],e1[r])
        t2_1=A(S0(a1),mj(a1,b1,c1))
        s1_check=[A(t1_1,t2_1),a1,b1,c1,A(d1,t1_1),ee1,f1,g1]

        a2,b2,c2,d2,e2,f2,g2,h2 = s2_check
        t1_2=A(h2,S1(e2),ch(e2,f2,g2),C[r],e2_full[r])
        t2_2=A(S0(a2),mj(a2,b2,c2))
        s2_check=[A(t1_2,t2_2),a2,b2,c2,A(d2,t1_2),e2,f2,g2]

        da = HW(s1_check[0]^s2_check[0])
        de = HW(s1_check[4]^s2_check[4])
        results.append((r+1, da, de))

    return w2, results


def run_experiment():
    print("="*60)
    print("BUDGET: break(1) + repair(15) = 16 words")
    print("Compare: a-repair vs e-repair (Wang)")
    print("="*60)

    for trial in range(3):
        w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]

        print(f"\n--- Trial {trial} ---")

        # a-repair
        _, res_a = break_and_repair(w1, repair_target='a')
        # e-repair
        _, res_e = break_and_repair(w1, repair_target='e')

        print(f"\n  {'r':>3} | a-repair         | e-repair (Wang)")
        print(f"      | δa   δe          | δa   δe")
        print("  " + "-"*50)

        for r in range(min(25, len(res_a))):
            ra = res_a[r]; re = res_e[r]
            na = "✓" if ra[1]==0 else ""
            ne = "✓" if re[2]==0 else ""
            # Пометки
            note_a = "a=0" if ra[1]==0 else ""
            note_e = "e=0" if re[2]==0 else ""
            if ra[1]==0 and ra[2]==0: note_a = "BOTH=0"
            if re[1]==0 and re[2]==0: note_e = "BOTH=0"

            if r < 20 or ra[1]==0 or ra[2]==0 or re[1]==0 or re[2]==0:
                print(f"  {ra[0]:3d} | {ra[1]:4d} {ra[2]:4d} {note_a:8s} | {re[1]:4d} {re[2]:4d} {note_e}")

        # Hash
        h1 = tuple(A(IV[i], list(IV)[i]) for i in range(8))  # placeholder
        # Actually compute hash
        from scf_my_way import run as sha_run
        w2a, _ = break_and_repair(w1, repair_target='a')
        w2e, _ = break_and_repair(w1, repair_target='e')
        h1=sha_run(w1); h2a=sha_run(w2a); h2e=sha_run(w2e)
        dha = sum(HW(h1[i]^h2a[i]) for i in range(8))
        dhe = sum(HW(h1[i]^h2e[i]) for i in range(8))
        n_diff_a = sum(1 for i in range(16) if w1[i]!=w2a[i])
        n_diff_e = sum(1 for i in range(16) if w1[i]!=w2e[i])

        print(f"\n  a-repair: dH={dha}, words_diff={n_diff_a}")
        print(f"  e-repair: dH={dhe}, words_diff={n_diff_e}")


if __name__ == '__main__':
    run_experiment()
