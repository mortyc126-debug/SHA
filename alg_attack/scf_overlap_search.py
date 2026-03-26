#!/usr/bin/env python3
"""
OVERLAP SEARCH: два входных бита с максимальным overlap влияния.

Если flip(bit_A) меняет множество S_A выходных бит,
и flip(bit_B) меняет множество S_B,
то flip(A и B вместе) меняет S_A △ S_B = (S_A ∪ S_B) \ (S_A ∩ S_B).
Размер: |S_A| + |S_B| - 2×|S_A ∩ S_B|.

Для COLLISION: |S_A △ S_B| = 0 → |S_A ∩ S_B| = (|S_A|+|S_B|)/2.
Если |S_A|=|S_B|=128 → нужно overlap = 128 → S_A = S_B.

Но S_A = S_B означает: два входных бита влияют на ОДИНАКОВЫЕ выходы.
XOR двух flips = flip(A)⊕flip(B) на КАЖДОМ affected бите.
Если flip(A) и flip(B) дают ОДИНАКОВОЕ изменение → XOR = 0 → collision!

НО: flip(A) и flip(B) дают РАЗНЫЕ значения δ (не просто 1).
Overlap = "оба бита ВЛИЯЮТ на этот выход" не значит "одинаково".

Нужно: не overlap, а СОВПАДЕНИЕ δ-ЗНАЧЕНИЙ.
δH(flip A ⊕ flip B) = δH(flip A) ⊕ δH(flip B)  (над GF(2)!)

Collision: δH(flip A) = δH(flip B) → δH(A⊕B) = 0!

Это BIRTHDAY на δH-ЗНАЧЕНИЯХ 512 single-bit flips!
512 flip'ов, каждый даёт 256-bit δH.
Birthday на 512 элементах: P(match) ≈ 512²/(2×2²⁵⁶) ≈ 0.
Нет.

НО: δH не random 256 бит. HW(δH) ≈ 128. И они КЛАСТЕРИЗОВАНЫ.
Closest pair = 91 (раньше). При N=512: P(HW(δ(δ))=0) ≈ 0.

ДРУГОЙ подход: не два single-bit flip'а, а ПОДМНОЖЕСТВО из K бит.
2^K подмножеств. Каждое → δH. Birthday на 2^K элементах.
Для 2^K ≈ 2^{128}: K = 128 → birthday O(2^{64}).
Но 2^{128} подмножеств → 2^{128} SHA evaluations. Не лучше.

НО: я вижу ДРУГОЕ. Influence graph показал: overlap = 44.
Net change = |A|+|B| - 2×overlap.
Если overlap БОЛЬШЕ → net change МЕНЬШЕ.
Ищу pair с МАКСИМАЛЬНЫМ overlap → МИНИМАЛЬНЫМ net change.
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

def sha(w):
    e=expand(w); s=list(IV)
    for r in range(64):
        a,b,c,d,ee,f,g,h=s
        t1=A(h,S1(ee),ch(ee,f,g),C[r],e[r])
        t2=A(S0(a),mj(a,b,c))
        s=[A(t1,t2),a,b,c,A(d,t1),ee,f,g]
    return tuple(A(IV[i],s[i]) for i in range(8))


def compute_all_dH(w1):
    """Для каждого single-bit flip: вычислить δH как 256-bit int."""
    h1 = sha(w1)
    h1_int = 0
    for i in range(8): h1_int |= h1[i] << (i*32)

    dH_map = {}
    for word in range(16):
        for bit in range(32):
            w2 = list(w1); w2[word] ^= (1<<bit)
            h2 = sha(w2)
            h2_int = 0
            for i in range(8): h2_int |= h2[i] << (i*32)
            dH_map[(word,bit)] = h1_int ^ h2_int

    return dH_map


def search_overlap(N):
    print("="*60)
    print("OVERLAP + δH MATCHING")
    print("="*60)

    for trial in range(N):
        w1 = [int.from_bytes(os.urandom(4),'big') for _ in range(16)]
        dH_map = compute_all_dH(w1)

        keys = list(dH_map.keys())
        dH_vals = [dH_map[k] for k in keys]

        # 1. Closest δH pair (XOR distance)
        best_dist = 256; best_pair = None

        for i in range(len(keys)):
            for j in range(i+1, len(keys)):
                xor_dH = dH_vals[i] ^ dH_vals[j]
                dist = HW(xor_dH)
                if dist < best_dist:
                    best_dist = dist
                    best_pair = (keys[i], keys[j])

        print(f"\n  Trial {trial}:")
        if best_pair:
            k1, k2 = best_pair
            hw1 = HW(dH_map[k1]); hw2 = HW(dH_map[k2])
            print(f"    Closest δH pair: W[{k1[0]}][{k1[1]:2d}] (HW={hw1}) "
                  f"vs W[{k2[0]}][{k2[1]:2d}] (HW={hw2})")
            print(f"    XOR distance: {best_dist}")
            print(f"    → flip BOTH: δH has HW={best_dist} (this IS the near-collision!)")

            if best_dist < 90:
                print(f"    ★ Below 90!")

            # Verify: actually flip both and check
            w2_both = list(w1)
            w2_both[k1[0]] ^= (1 << k1[1])
            w2_both[k2[0]] ^= (1 << k2[1])
            h1 = sha(w1); h2 = sha(w2_both)
            actual_dh = sum(HW(h1[i]^h2[i]) for i in range(8))
            print(f"    Actual δH(flip both): HW={actual_dh}")
            print(f"    Theory (XOR of single flips): HW={best_dist}")
            print(f"    Match: {'YES' if actual_dh == best_dist else 'NO — nonlinear!'}")
            print(f"    Difference: {actual_dh - best_dist:+d}")

        # 2. Triplets: find 3 flips where XOR of all 3 δH is smallest
        best_triple = 256; best_triple_keys = None
        # Sample random triples
        for _ in range(5000):
            i = int.from_bytes(os.urandom(2),'big') % len(keys)
            j = int.from_bytes(os.urandom(2),'big') % len(keys)
            k = int.from_bytes(os.urandom(2),'big') % len(keys)
            if i==j or j==k or i==k: continue
            xor3 = dH_vals[i] ^ dH_vals[j] ^ dH_vals[k]
            d = HW(xor3)
            if d < best_triple:
                best_triple = d
                best_triple_keys = (keys[i], keys[j], keys[k])

        if best_triple_keys:
            k1,k2,k3 = best_triple_keys
            print(f"\n    Best triple: W[{k1[0]}][{k1[1]}], W[{k2[0]}][{k2[1]}], W[{k3[0]}][{k3[1]}]")
            print(f"    XOR₃ distance: {best_triple}")

            w3 = list(w1)
            w3[k1[0]]^=(1<<k1[1]); w3[k2[0]]^=(1<<k2[1]); w3[k3[0]]^=(1<<k3[1])
            h1=sha(w1); h3=sha(w3)
            actual3 = sum(HW(h1[i]^h3[i]) for i in range(8))
            print(f"    Actual δH(flip triple): HW={actual3}")
            print(f"    Linear prediction: HW={best_triple}")
            print(f"    Nonlinear error: {actual3 - best_triple:+d}")


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    search_overlap(N)
