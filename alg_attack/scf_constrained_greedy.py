#!/usr/bin/env python3
"""
SCF: CONSTRAINED GREEDY — запрещаем тривиальный путь.

Открытие: единственный минимум mismatch = δW=0 (тривиальный).
Решение: зафиксировать ЧАСТЬ δW и оптимизировать остальное.

Стратегия: фиксируем δW[0] = 0x80000000 (Wang-style, неотменяемый).
Оптимизируем W2[1..15] для минимизации mismatch.
Greedy solver НЕ МОЖЕТ отменить δW[0] → ищет НЕТРИВИАЛЬНЫЙ путь.
"""
import os, sys, math

MASK32 = 0xFFFFFFFF
K = [
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

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&MASK32
def Sig0(x): return rotr(x,2)^rotr(x,13)^rotr(x,22)
def Sig1(x): return rotr(x,6)^rotr(x,11)^rotr(x,25)
def sig0(x): return rotr(x,7)^rotr(x,18)^(x>>3)
def sig1(x): return rotr(x,17)^rotr(x,19)^(x>>10)
def Ch(e,f,g): return (e&f)^(~e&g)&MASK32
def Maj(a,b,c): return (a&b)^(a&c)^(b&c)
def add32(*args):
    s=0
    for x in args: s=(s+x)&MASK32
    return s
def sub32(a,b): return (a-b)&MASK32
def hw(x): return bin(x).count('1')

def expand_real(W16):
    W=list(W16)
    for i in range(16,64):
        W.append(add32(sig1(W[i-2]),W[i-7],sig0(W[i-15]),W[i-16]))
    return W

def sha_all_states(W16):
    W=expand_real(W16); s=list(IV); states=[list(s)]
    for r in range(64):
        a,b,c,d,e,f,g,h=s
        T1=add32(h,Sig1(e),Ch(e,f,g),K[r],W[r])
        T2=add32(Sig0(a),Maj(a,b,c))
        s=[add32(T1,T2),a,b,c,add32(d,T1),e,f,g]
        states.append(list(s))
    return states, W

def sha_compress(W16):
    st, _ = sha_all_states(W16)
    return [add32(IV[i], st[64][i]) for i in range(8)]

def hash_hw(W1, W2):
    H1 = sha_compress(W1); H2 = sha_compress(W2)
    return sum(hw(H1[i]^H2[i]) for i in range(8))


def constrained_greedy(W1, fixed_delta_word=0, fixed_delta_val=0x80000000,
                       max_rounds=20):
    """Greedy solver with FIXED δW[word] — can't undo the perturbation."""
    W2 = list(W1)
    W2[fixed_delta_word] ^= fixed_delta_val

    best_hw = hash_hw(W1, W2)
    best_W2 = list(W2)
    current_hw = best_hw

    for round_num in range(max_rounds):
        # Compute influence for each ALLOWED bit (not fixed_delta bits)
        influences = []
        for word in range(16):
            for bit in range(32):
                # Skip: can't flip the fixed delta bit
                if word == fixed_delta_word and bit == 31 and fixed_delta_val == 0x80000000:
                    continue

                W2_flip = list(W2)
                W2_flip[word] ^= (1 << bit)

                # Ensure δW stays nontrivial
                if W2_flip == W1:
                    continue

                flip_hw = hash_hw(W1, W2_flip)
                influences.append((flip_hw - current_hw, word, bit, flip_hw))

        influences.sort()

        # Flip best helpful bit
        if influences and influences[0][0] < 0:
            delta, word, bit, new_hw = influences[0]
            W2[word] ^= (1 << bit)
            current_hw = new_hw
            if current_hw < best_hw:
                best_hw = current_hw
                best_W2 = list(W2)

            if round_num < 10 or round_num % 5 == 0:
                n_diff = sum(1 for i in range(16) if W1[i]!=W2[i])
                print(f"  Round {round_num:2d}: HW(δH)={current_hw:3d}, "
                      f"flipped W[{word}][{bit:2d}] (Δ={delta:+d}), "
                      f"words_diff={n_diff}")

            if current_hw == 0:
                print(f"  ★★★ COLLISION AT ROUND {round_num}!")
                break
        else:
            print(f"  Round {round_num:2d}: no helpful bits, HW={current_hw}")
            break

    return best_W2, best_hw


def multi_start_constrained(N_starts):
    """Run constrained greedy from multiple starting points."""
    print("="*70)
    print("MULTI-START CONSTRAINED GREEDY")
    print(f"  Fixed: δW[0] = 0x80000000 (MSB flip, can't undo)")
    print(f"  Free: W2[1..15] optimized via greedy bit-flip")
    print("="*70)

    global_best = 256
    global_best_W1 = None
    global_best_W2 = None

    for start in range(N_starts):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        print(f"\n  --- Start {start} ---")
        W2_opt, best_hw = constrained_greedy(W1, max_rounds=20)

        if best_hw < global_best:
            global_best = best_hw
            global_best_W1 = list(W1)
            global_best_W2 = list(W2_opt)
            print(f"  ★ NEW GLOBAL BEST: HW = {global_best}")

    # SA comparison
    print(f"\n  SA comparison (same constraint):")
    for _ in range(3):
        W1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        W2 = list(W1); W2[0] ^= 0x80000000
        best_sa = hash_hw(W1, W2)
        cur = list(W2)
        for it in range(5000):
            t = list(cur)
            w = int.from_bytes(os.urandom(1),'big') % 16
            b = int.from_bytes(os.urandom(1),'big') % 32
            if w==0 and b==31: continue  # Don't undo
            t[w] ^= (1<<b)
            if t==W1: continue
            s = hash_hw(W1, t)
            T = max(0.01, 1-it/5000)
            if s < best_sa or math.exp(-(s-best_sa)/(T*2)) > int.from_bytes(os.urandom(4),'big')/(1<<32):
                cur=t
                if s < best_sa: best_sa = s
        print(f"    SA: HW = {best_sa}")

    print(f"\n  GLOBAL BEST (greedy): HW = {global_best}")
    return global_best_W1, global_best_W2, global_best


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    W1, W2, best = multi_start_constrained(N)

    if W1 and W2:
        n_diff = sum(1 for i in range(16) if W1[i]!=W2[i])
        print(f"\n  Best pair: {n_diff} words differ, HW(δH)={best}")
        print(f"  W1[0] = 0x{W1[0]:08x}")
        print(f"  W2[0] = 0x{W2[0]:08x}")
        print(f"  δW[0] = 0x{W1[0]^W2[0]:08x}")
