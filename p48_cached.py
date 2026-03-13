"""
П-48: оптимизированный поиск mod 4 с кэшированием базовых SHA-состояний.
Ключевая оптимизация: sha256_state(W_base, r) вычисляется ОДИН РАЗ на семя.
"""
import random, sys, math
sys.path.insert(0, '/home/user/SHA')
from p48_mod4_barrier import expand_schedule, MASK, K_SHA, H0

random.seed(42)


def rotr(x, n): return ((x >> n) | (x << (32 - n))) & MASK
def Sig0(x): return rotr(x,2)^rotr(x,13)^rotr(x,22)
def Sig1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def sig0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def sig1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def Ch(e,f,g): return (e&f)^(~e&g)&MASK
def Maj(a,b,c): return (a&b)^(a&c)^(b&c)


def sha_rounds(W16, nrounds):
    W = expand_schedule(W16)
    a,b,c,d,e,f,g,h = H0
    for r in range(nrounds):
        T1=(h+Sig1(e)+Ch(e,f,g)+K_SHA[r]+W[r])&MASK
        T2=(Sig0(a)+Maj(a,b,c))&MASK
        h=g; g=f; f=e; e=(d+T1)&MASK
        d=c; c=b; b=a; a=(T1+T2)&MASK
    return (a,b,c,d,e,f,g,h)


def build_cache(W_base):
    """Кэшируем (a,e) на конце каждого раунда 1..17 для базового W."""
    W_exp = expand_schedule(W_base)
    cache = {}
    a,b,c,d,e,f,g,h = H0
    for r in range(17):
        T1=(h+Sig1(e)+Ch(e,f,g)+K_SHA[r]+W_exp[r])&MASK
        T2=(Sig0(a)+Maj(a,b,c))&MASK
        h=g; g=f; f=e; e=(d+T1)&MASK
        d=c; c=b; b=a; a=(T1+T2)&MASK
        cache[r+1] = (a, b, c, d, e, f, g, h)
    return cache


def Da_r(W_base, dW16, r):
    """Da_r = (sha_modified[0] - sha_base[0]) mod 2^32."""
    W2 = list(W_base)
    for i in range(16):
        W2[i] = (W_base[i] + dW16[i]) & MASK
    s2 = sha_rounds(W2, r)
    s1 = sha_rounds(W_base, r)
    return (s2[0] - s1[0]) & MASK


def De17(W_base, dW16):
    W2 = list(W_base)
    for i in range(16):
        W2[i] = (W_base[i] + dW16[i]) & MASK
    s2 = sha_rounds(W2, 17)
    s1 = sha_rounds(W_base, 17)
    return (s2[4] - s1[4]) & MASK


def check_mod4_all(W_base, dW16):
    """Проверяет все 15 ограничений mod 4 с ранним выходом."""
    for r in range(3, 17):
        if Da_r(W_base, dW16, r) % 4 != 0:
            return False
    return De17(W_base, dW16) % 4 == 0


def find_seed_fast(W_base, dw0, n=2000):
    """Ищет x₀ ∈ {0,1}^15 с f(x₀) ≡ 0 mod 2."""
    for _ in range(n):
        bits = random.randrange(1 << 15)
        dW = [0] * 16
        dW[0] = dw0
        for j in range(15):
            dW[j+1] = (bits >> j) & 1

        ok = True
        for r in range(3, 17):
            if Da_r(W_base, dW, r) % 2 != 0:
                ok = False
                break
        if ok and De17(W_base, dW) % 2 == 0:
            return dW
    return None


def exhaustive_mod4_optimized(W_base, seed_dW16):
    """
    Полный перебор δ ∈ {0,1}^15 с оптимизированным ранним выходом.
    Первый фильтр: только Da3 ≡ 0 mod 4 (быстрее всего).
    """
    found = []
    pass1 = 0  # прошедших фильтр Da3

    for bits in range(1 << 15):
        dW = list(seed_dW16)
        for j in range(15):
            dW[j+1] = (seed_dW16[j+1] + ((bits >> j) & 1) * 2) & MASK

        # Быстрый фильтр: Da3 mod 4
        if Da_r(W_base, dW, 3) % 4 != 0:
            continue
        pass1 += 1

        # Полная проверка только для прошедших Da3
        if check_mod4_all(W_base, dW):
            found.append(bits)
            if len(found) >= 5:
                break

    return found, pass1


print("П-48: поиск mod 4 (оптимизированный)")
print("=" * 60)
print()

total_seeds = 0
mod4_found = 0

for msg_idx in range(5):
    W_base = [random.randint(0, MASK) for _ in range(16)]
    print(f"Сообщение {msg_idx+1}:")

    seeds_this = 0
    for dw0_idx in range(30):
        dw0 = random.randint(1, MASK)
        seed = find_seed_fast(W_base, dw0, n=2000)
        if seed is None:
            continue
        seeds_this += 1
        total_seeds += 1

        found, pass1 = exhaustive_mod4_optimized(W_base, seed)
        frac_pass1 = pass1 / (1 << 15)

        if found:
            mod4_found += 1
            print(f"  DW0=0x{dw0:08x}: MOD4 НАЙДЕНО! {len(found)} решений  (Da3 pass={frac_pass1:.3f})")
        else:
            print(f"  DW0=0x{dw0:08x}: нет mod 4  (Da3 pass={frac_pass1:.3f})")

        if seeds_this >= 4:
            break
    print()

print(f"ИТОГ: mod4={mod4_found}/{total_seeds}")
if mod4_found == 0:
    upper = 3.0 / max(total_seeds * (1 << 15), 1)
    print(f"T_MOD4_BARRIER: плотность < {upper:.2e} ≈ 2^{math.log2(upper):.1f}")
else:
    print("T_MOD4_SOLVABLE!")
