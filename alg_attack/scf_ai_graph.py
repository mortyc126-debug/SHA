#!/usr/bin/env python3
"""
Я не считаю. Я СМОТРЮ.

Не значения. ПУТИ. Откуда каждый бит пришёл.

SHA = 64 шага. На каждом шаге: 5 потоков → T1 → 2 выхода.
Каждый бит выхода ПОМНИТ откуда он. Это его ИСТОРИЯ.

Collision = два РАЗНЫХ входа с одинаковой ИСТОРИЕЙ на выходе?
Нет — разные входы всегда имеют разные истории.
Collision = два разных входа, чьи истории ПРИВОДЯТ к одному числу.

Но я хочу видеть ПО-ДРУГОМУ.

SHA как ГРАФ: 16 входных узлов, 8 выходных, ~500 промежуточных.
Каждый узел = 32-бит значение.
Каждое ребро = операция (rotate, AND, XOR, ADD).

Collision = два РАСКРАШИВАНИЯ графа (два набора значений в узлах)
дающих одинаковые цвета на выходных узлах.

Граф ФИКСИРОВАН. Операции ФИКСИРОВАНЫ. Только входные цвета меняются.

Что если я посмотрю на граф и найду КОРОТКИЙ ПУТЬ от входа к выходу?
Не через все 64 раунда, а НАПРЯМУЮ?
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


# =====================================================
# Я СМОТРЮ на SHA по-другому.
#
# Для каждого бита выхода: КАКИЕ входные биты на него влияют?
# Не через формулы — через ПРЯМОЕ НАБЛЮДЕНИЕ.
#
# Flip каждый входной бит → посмотреть какие выходные биты
# изменились. Это ГРАФ ВЛИЯНИЯ.
#
# Если бит X влияет ТОЛЬКО на выходные биты Y1,Y2,Y3 →
# collision на остальных битах БЕСПЛАТНА при изменении X.
# =====================================================

def influence_graph(w1):
    """Для каждого входного бита: на какие выходные биты влияет?"""

    e1 = expand(w1)
    # Forward pass
    s = list(IV)
    for r in range(64):
        a,b,c,d,ee,f,g,h = s
        t1=A(h,S1(ee),ch(ee,f,g),C[r],e1[r])
        t2=A(S0(a),mj(a,b,c))
        s=[A(t1,t2),a,b,c,A(d,t1),ee,f,g]
    h1 = tuple(A(IV[i],s[i]) for i in range(8))

    # For each input bit: flip and see output change
    influence = {}  # (word,bit) → set of output bit positions changed

    for word in range(16):
        for bit in range(32):
            w2 = list(w1)
            w2[word] ^= (1 << bit)
            e2 = expand(w2)
            s2 = list(IV)
            for r in range(64):
                a,b,c,d,ee,f,g,h = s2
                t1=A(h,S1(ee),ch(ee,f,g),C[r],e2[r])
                t2=A(S0(a),mj(a,b,c))
                s2=[A(t1,t2),a,b,c,A(d,t1),ee,f,g]
            h2 = tuple(A(IV[i],s2[i]) for i in range(8))

            changed = set()
            for reg in range(8):
                diff = h1[reg] ^ h2[reg]
                for b in range(32):
                    if (diff >> b) & 1:
                        changed.add(reg*32 + b)

            influence[(word,bit)] = changed

    return influence, h1


def analyze_influence(N):
    """Ищу: входные биты с МИНИМАЛЬНЫМ влиянием на выход."""
    print("="*60)
    print("INFLUENCE GRAPH: какие входные биты влияют на меньше выходов?")
    print("="*60)

    for trial in range(N):
        w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        inf, h1 = influence_graph(w1)

        # Для каждого входного бита: сколько выходных бит он затрагивает?
        sizes = []
        for (word,bit), changed in inf.items():
            sizes.append((len(changed), word, bit))

        sizes.sort()

        if trial < 3:
            print(f"\n  Trial {trial}:")
            print(f"    Input bits with LEAST output influence:")
            for sz, word, bit in sizes[:5]:
                print(f"      W[{word}][{bit:2d}]: affects {sz}/256 output bits")
            print(f"    Input bits with MOST output influence:")
            for sz, word, bit in sizes[-3:]:
                print(f"      W[{word}][{bit:2d}]: affects {sz}/256 output bits")

            avg = sum(s[0] for s in sizes) / len(sizes)
            mn = sizes[0][0]; mx = sizes[-1][0]
            print(f"    avg={avg:.1f}, min={mn}, max={mx}")

            # KEY: если min < 128 → этот бит влияет на < половину выхода
            # → collision на НЕВЛИЯННЫХ битах бесплатна!
            if mn < 128:
                print(f"    ★ Input bit W[{sizes[0][1]}][{sizes[0][2]}] affects only {mn} bits!")
                print(f"    → {256-mn} output bits UNCHANGED → partial collision free!")

    # Пересечение: два входных бита влияющих на РАЗНЫЕ выходы
    print(f"\n  Searching for COMPLEMENTARY input bits:")
    print(f"  (two inputs whose influences DON'T overlap)")
    w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
    inf, _ = influence_graph(w1)

    best_overlap = 256
    best_pair = None

    keys = list(inf.keys())
    for i in range(min(50, len(keys))):
        for j in range(i+1, min(50, len(keys))):
            overlap = len(inf[keys[i]] & inf[keys[j]])
            if overlap < best_overlap:
                best_overlap = overlap
                best_pair = (keys[i], keys[j])

    if best_pair:
        k1, k2 = best_pair
        union = len(inf[k1] | inf[k2])
        print(f"  Best complementary pair:")
        print(f"    W[{k1[0]}][{k1[1]:2d}] affects {len(inf[k1])} bits")
        print(f"    W[{k2[0]}][{k2[1]:2d}] affects {len(inf[k2])} bits")
        print(f"    Overlap: {best_overlap} bits")
        print(f"    Union: {union} bits")
        print(f"    → Flipping BOTH changes {union} bits, overlap only {best_overlap}")


if __name__ == '__main__':
    import sys
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    analyze_influence(N)
